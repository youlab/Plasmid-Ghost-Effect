# Keio_NGS_data_processing

How to go from **raw amplicon reads** (NCBI SRA) to the **strain barcode count tables** that the
analysis scripts in [Keio_Community/](../Keio_Community/) use.

Each *E. coli* Keio strain in the synthetic communities (Comm57 and Comm87) carries a unique
chromosomal barcode. One round of PCR adds a well-specific index pair (i5 + i7) to the barcode
amplicon, all wells are pooled, and the pool is sequenced. Processing therefore means:
**demultiplex by index pair → count strain barcodes per well → merge into one table.**

---

## 1. Get the data

Raw reads are deposited under NCBI BioProject
[PRJNA1360409](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA1360409) (study SRP644022).
The two runs `NGS of E.coli Keio Comm57 / Comm87` are the Keio barcode-amplicon libraries; the remaining runs in the BioProject
belong to sink isolate community experiments in the paper.

The deposited reads are already separated by community, so each run is processed independently:

```bash
# SRA Toolkit
prefetch SRR39821239 SRR39821238
fasterq-dump --split-files SRR39821239   # -> SRR39821239_1.fastq (R1), SRR39821239_2.fastq (R2)
fasterq-dump --split-files SRR39821238
```

Inside a run, reads from all wells are still mixed — that is what step 3 takes care of.

## 2. What is in this folder

| File | Purpose |
|---|---|
| `Galaxy-Workflow-Plasmid_Barcode_1-cycle_PCR_with_PhiX_Filter_Dec_2024.ga` | Galaxy workflow: reads → per-well barcode counts |
| `i5_16.txt` | 16 i5 (P5) index sequences, 8 bp, `<id>\t<sequence>` |
| `i7.txt` | 24 i7 (P7) index sequences, 8 bp, `<id>\t<sequence>` |
| `Comm57_barcodes.txt` | 60 strain barcodes expected in Comm57; IDs are Keio strain numbers |
| `Comm87_barcodes.txt` | 94 strain barcodes expected in Comm87; IDs are Keio strain numbers |
| `Metadata.xlsx` | Plate map: which sample sits in each i5 × i7 well |
| `Merge_Count.py` | Merges the per-well count files into one `.csv` (+ a QC heatmap) |

Each strain barcode is 52 bp: an 18 bp variable barcode, the 16 bp constant linker
`CCTCAGGGTCACTAGG`, and a second 18 bp variable barcode.

### Where each tag sits in the reads (1-based, inclusive)

```
R1:  [ 1-52  strain barcode ] ............ [ 87-94  i7 index ] ...
R2:  ........................ ........... [ 87-94  i5 index ] ...
```

Because the insert is short, both reads run through the adapter into the sample index, which is why
the indices are read out of the middle of R1/R2 rather than from separate index reads.

### Which wells belong to which community

* **Comm87** = i5 9–16 × i7 13–24 (i5 14–16 unused)
* **Comm57** = all remaining wells

Every populated cell of `Metadata.xlsx` names the sample in that well (community, time point,
plasmid, and dilution rate, 1:100 or 1:500). Cells shaded green are the communities that received
the antibiotic pulse between day 2 and day 3.

## 3. Count barcodes with the Galaxy workflow

Upload the workflow to a [Galaxy](https://usegalaxy.org) instance and run it **once per community**.

Inputs:

| Workflow input | Use |
|---|---|
| Read 1 / Read 2 | the run's R1 / R2 FASTQ |
| phiX Genome | phiX174 reference FASTA (spike-in) |
| i5.txt | `i5_16.txt` |
| i7.txt | `i7.txt` |
| Barcodes (Barcodes.txt) | `Comm57_barcodes.txt` **or** `Comm87_barcodes.txt`, matching the run |

What it does:

1. **BBDuk** removes phiX spike-in reads from the pair.
2. Trims R2 to positions 87–94 → i5 tag, pastes it onto the front of R1.
3. **Barcode Splitter** on that tag (5′ end, ≤1 mismatch) → one dataset per i5.
4. Trims positions 95–102 of the joined read (= R1 87–94) → i7 tag, pastes it on too.
5. **Barcode Splitter** on the i7 tag (≤1 mismatch) → one dataset per i5 × i7 well.
6. Trims to the 52 bp strain barcode and splits it against the community barcode file
   (≤2 mismatches). The **"Summary: split barcodes"** output is the count table for that well.

Each summary is a two-column tabular file listing every expected barcode, plus `unmatched` and
`total` rows:

```
# Barcode   Count
1           3201
11          74
...
unmatched   50
total       6102
```

## 4. Merge the per-well counts

Download the summaries and lay them out as `<i5>/<i7>.tabular` under a `raw_data/` folder, i.e. one
subfolder per i5 index containing one file per i7 index:

```
<working folder>/
├── raw_data/
│   ├── 1/  1.tabular … 24.tabular
│   ├── 2/  1.tabular … 24.tabular
│   └── …
└── processed_data/
```

Then run:

```bash
python Merge_Count.py
```

Outputs (into `processed_data/`):

* `merged_counts.csv` — rows = strain barcodes, columns = `<i5>-<i7>` wells, values = read counts.
  Missing wells are filled with zeros.
* `merged_counts.png` — heatmap of log₁₀ total reads per well (P7 × P5), a quick check that every
  expected well was sequenced deeply enough.

Two things to know before running:

* `Merge_Count.py` reads the header from `./raw_data/1/1.tabular` but the loop reads
  `<i5>/<i7>.tabular` — set both to the same prefix (or run from inside `raw_data/`) or every well
  silently comes back as zeros. [Keio_Community/Comm57/Merge_Count.py](../Keio_Community/Comm57/Merge_Count.py)
  is the same script and has the same mismatch; only
  [Keio_Community/Comm87/Merge_Count.py](../Keio_Community/Comm87/Merge_Count.py) is fixed, using a
  single `path` variable and a well range already restricted to Comm87 (i5 9–16 × i7 13–24) — the
  easier starting point when reproducing that community.
* Both copies under `Keio_Community/` have the `fig.savefig` line commented out, so they show the
  QC heatmap without writing `merged_counts.png`. The copy in this folder writes it.

The script imports [natsort](https://pypi.org/project/natsort/), which is included in
[environment.yml](../environment.yml).

## 5. Downstream analysis

`merged_counts.csv` is the entry point for the community-composition analyses:
`organize_data.py` reshapes it into the per-experiment Excel files, and `Community_Dynamics_*.py`
turns those into the strain-level dynamics figures. See the **Keio_Community** section of the
[main README](../README.md) for the full chain.
