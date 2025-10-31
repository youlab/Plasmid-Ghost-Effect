import numpy as np
import matplotlib.pyplot as plt
from scipy import stats

import numpy as np
from scipy import stats

def report_welch_t(a, b, group_names=("Group A", "Group B"), alpha=0.05):
    """
    Format a Welch's independent-samples t-test report:
    - Descriptives (mean±SD, n)
    - Welch t, df, exact p
    - Mean difference with 95% CI
    - Cohen's d and Hedges' g (small-sample correction)
    NaNs are ignored.
    """
    a = np.asarray(a, dtype=float)
    b = np.asarray(b, dtype=float)

    # Remove NaNs for all calculations
    a = a[np.isfinite(a)]
    b = b[np.isfinite(b)]

    n1, n2 = len(a), len(b)
    m1, m2 = np.mean(a), np.mean(b)
    s1, s2 = np.std(a, ddof=1), np.std(b, ddof=1)

    # Welch's t-test (two-sided by default)
    t_stat, p_value = stats.ttest_ind(a, b, equal_var=False, nan_policy="omit")

    # Welch-Satterthwaite df
    v1, v2 = s1**2, s2**2
    se = np.sqrt(v1/n1 + v2/n2)
    df_num = (v1/n1 + v2/n2)**2
    df_den = (v1**2 / (n1**2 * (n1 - 1))) + (v2**2 / (n2**2 * (n2 - 1)))
    df = df_num / df_den

    # Mean difference CI
    md = m1 - m2
    t_crit = stats.t.ppf(1 - alpha/2, df)
    ci_lo, ci_hi = md - t_crit*se, md + t_crit*se

    # Effect size: Cohen's d (pooled SD) and Hedges' g correction
    sp2 = ((n1 - 1) * v1 + (n2 - 1) * v2) / (n1 + n2 - 2)
    sp = np.sqrt(sp2)
    d = md / sp if sp > 0 else np.nan
    J = 1 - 3 / (4*(n1 + n2) - 9)  # Hedges' small-sample correction
    g = J * d if np.isfinite(d) else np.nan

    # p-value formatting (APA-ish)
    if p_value < 0.001:
        p_str = "p<.001"
    else:
        p_str = f"p={p_value:.3f}".replace("0.", ".")

    # Assemble the sentence
    direction = "higher" if m1 > m2 else "lower"
    text = (
        f"{group_names[0]} (n={n1}) had a {direction} mean than {group_names[1]}: "
        f"{m1:.1f}±{s1:.1f} vs {m2:.1f}±{s2:.1f}. "
        f"Welch’s t-test: t({df:.2f})={t_stat:.2f}, {p_str}, "
        f"mean difference={md:.1f} [95% CI: {ci_lo:.1f}, {ci_hi:.1f}], "
        f"Cohen’s d={d:.2f} (Hedges’ g={g:.2f})."
    )
    return p_value, text


colors = ['#8da0cb', '#fc8d62']

def p_to_star(p):
    if p < 0.001:
        return "***"
    elif p < 0.01:
        return "**"
    elif p < 0.05:
        return "*"
    else:
        return "ns"

def loglinear_crossing_time(
        t: np.ndarray,
        mu: np.ndarray,
        se: np.ndarray,
        mu_star: float
    ):
    """
    Log-linear interpolation on ln(mu) with delta-method SE,
    robust to NaNs in mu or sd (rows with NaNs are skipped).
    """
    # 0. Remove rows with NaN in mu or sd --------------------------------------
    nan_idx = np.isnan(mu)
    t   = t[~nan_idx]
    mu  = mu[~nan_idx]
    se  = se[~nan_idx]

    if len(t) < 2:
        #print("not enough points left")
        return np.nan, np.nan                      # not enough points left

    # 1. Move to log space ---------------------------------
    g     = np.log(mu)
    se_g  = se / mu                               # delta: σ/μ
    g_star = np.log(mu_star)

    # 2. Locate the bracketing segment ----------------------------------------
    idx = np.where(g <= g_star)[0]                # first point below threshold
    if len(idx) == 0 or idx[0] == 0:
        #print("threshold never reached")
        return t[-1], np.nan                     # threshold never reached

    i  = idx[0] - 1                               # segment [i, i+1]
    dt = t[i + 1] - t[i]
    N  = g[i] - g_star
    D  = g[i] - g[i + 1]

    t_cross = t[i] + dt * N / D                   # log-linear interpolation

    # 3. Delta-method SE -------------------------------------------------------
    se_i, se_ip1 = se_g[i], se_g[i + 1]
    var_t  = (dt**2 / D**4) * ((g_star - g[i + 1])**2 * se_i**2 +
                               N**2 * se_ip1**2)
    se_cross = np.sqrt(var_t)
    return t_cross, se_cross

# Comm87 + R388
fig,ax=plt.subplots(1,1,figsize=(2.05, 1.9))

abundance = np.load(f"./processed_data/R388_mean.npy")
se_abundance = np.load(f"./processed_data/R388_se.npy")

# define a hard floor for plasmid abundance (0.1%), does not affect the estimation of half-lives, but as a quick fix
# for performing log-linear interpolation
zero_idx = (abundance < 0.1)
abundance[zero_idx] = 0.1
time = np.array([0, 2, 3, 7, 10,15, 20,25, 30,36])

n_ab = abundance.shape[0]  # with / without antibiotic pulse
bio_rep = abundance.shape[1]  # number of biological replicates

for j in range(n_ab):
    zorder = -j
    for k in range(bio_rep):
        ax.plot(time, abundance[j,k,:], color=colors[j], linewidth=1.5, zorder=zorder)
        ax.scatter(time, abundance[j,k,:], s=60, color=colors[j], linewidth=1, edgecolors='black', zorder=zorder)
        ax.errorbar(time, abundance[j,k,:], se_abundance[j,k,:], marker='None', linewidth=0, elinewidth=1, capsize=2, color="k", zorder=zorder)
ax.fill_between(x=[2, 3], y1=[1, 1], y2=[150, 150], color="#808080", zorder=-20, alpha=0.3)#, lw=0.2)
ax.text(x=0.9,y=0.2,s="Comm87\nR388",ha="right",transform=ax.transAxes)
ax.set_xlim([-1, 37])
ax.set_xticks([0, 18,36])
ax.set_ylim([1, 150])
ax.set_yticks([1, 10, 100])
ax.set_yscale("log")
ax.set_xlabel("time (days)")
ax.set_ylabel("P%",rotation=0,va="center",ha="right")
fig.subplots_adjust(left=0.25, right=0.95, bottom=0.23, top=0.95)
fig.savefig("./figures/Comm87_plasmid_dynamics.png",dpi=300)
fig.savefig("./figures/Comm87_plasmid_dynamics.svg")

fig2, ax2 = plt.subplots(1, 1, figsize=(2.05, 1.9))
HL = np.zeros((n_ab, bio_rep))
SE_HL = np.zeros((n_ab, bio_rep))
# For LB group: calculate half-life from day 0
P0 = 25
t_start = 0
for k in range(bio_rep):
    t, t_se = loglinear_crossing_time(time[t_start:], abundance[0, k, t_start:], se_abundance[0, k, t_start:], P0/2)
    HL[0, k] = t
    SE_HL[0, k] = t_se
# For LB+Trim group: calculate half-life from day 3
P0 = 100
t_start = 2
for k in range(bio_rep):
    t, t_se = loglinear_crossing_time(time[t_start:], abundance[1, k, t_start:], se_abundance[1, k, t_start:], P0/2)
    HL[1, k] = t - 3  # start counting after the antibiotic pulse
    SE_HL[1, k] = t_se
p_value, report = report_welch_t(HL[1,:], HL[0,:], group_names=("+Trim","LB"))
print(f"Comm87+R388 Welch's t-test")
print(report)

# HL: shape (n_ab, bio_rep)
# SE_HL: shape (n_ab, bio_rep)
n_ab, bio_rep = HL.shape

# --- summary stats ---
means = HL.mean(axis=1)
sems  = HL.std(axis=1, ddof=1) / np.sqrt(bio_rep)  # SEM across biological replicates
# --- plot ---
x = np.arange(n_ab)

# bars + error bars
error_kw = {'elinewidth': 1.2, 'capsize': 3, 'capthick': 1.2, 'color': 'blue'}
bars = ax2.bar(x, means, yerr=sems, width=0.6, error_kw=error_kw, lw=1.2, facecolor="None", edgecolor="k", zorder=2)

# jittered points (strip-like)
rng = np.random.default_rng(42)  # set seed for reproducible jitter
for j in range(n_ab):
    jitter_scale = 0.1 if j==0 else 0.04
    x_j = rng.normal(loc=x[j], scale=jitter_scale, size=bio_rep)
    ax2.scatter(x_j, HL[j], s=60, linewidths=1, facecolor=colors[j], edgecolors='k', zorder=3)
    ax2.errorbar(x_j, HL[j], SE_HL[j], marker='o', markersize=0, elinewidth=1.2, linewidth=0,
                 capsize=2, capthick=1.2, color="k", zorder=4)
ax2.text(0.5, np.max(HL)*1.08, p_to_star(p_value),ha='center', va='bottom', fontsize=12)
ax2.plot(x,np.ones(2)*np.max(HL)*1.13,c="k",lw=1)

ax2.set_ylim([0,np.max(HL)*1.3])
ax2.set_xticks(x)
ax2.set_xticklabels(["no pulse","+Trim"])
ax2.set_ylabel(r"$\tau_{1/2}$ (days)")
fig2.subplots_adjust(left=0.25, right=0.95, bottom=0.15, top=0.87)
fig2.savefig("./figures/Comm87_half_life.png",dpi=300)
fig2.savefig("./figures/Comm87_half_life.svg")
plt.show()
