import subprocess
import sys

plasmids = ["pSC101", "colE1", "pUC"]
days = [3, 20]

for day in days:
    for plasmid in plasmids:
        print(f"\n{'='*60}")
        print(f"Running: day={day}, plasmid={plasmid}")
        print(f"{'='*60}\n")
        
        # Read the original file
        with open("GFP_evolved_fitness.py", "r") as f:
            content = f.read()
        
        # Modify the day and plasmid variables
        import re
        content_modified = re.sub(r'day = \d+', f'day = {day}', content)
        content_modified = re.sub(r'plasmid = "[^"]*"', f'plasmid = "{plasmid}"', content_modified)
        
        # Write to a temporary file
        with open("_temp_fitness.py", "w") as f:
            f.write(content_modified)
        
        # Run the temporary file
        result = subprocess.run([sys.executable, "_temp_fitness.py"], cwd=".")
        
        if result.returncode != 0:
            print(f"Error running day={day}, plasmid={plasmid}")
        else:
            print(f"Successfully completed day={day}, plasmid={plasmid}")

print("\n" + "="*60)
print("All fitness calculations completed!")
print("="*60)
