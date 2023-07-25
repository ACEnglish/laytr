import os
import subprocess

# Path to the TSV file
tsv_file = "v3.1-GRCh38-subset.tsv"

# Output BED file
output_bed_file = "v3.1-GRCh38-subset-multiintersect.bed"

# Create a temporary directory to store intermediate BED files
temp_dir = "temp"
os.makedirs(temp_dir, exist_ok=True)

# Read the TSV file and extract stratification IDs and relative paths
stratifications = []
with open(tsv_file, "r") as f:
    for line in f:
        strat_id, path = line.strip().split("\t")
        stratifications.append((strat_id, path))

# Generate individual BED files for each stratification
bed_files = []
for strat_id, path in stratifications:
    temp_bed_file = os.path.join(temp_dir, f"{strat_id}.bed")
    subprocess.run(["gunzip", "-c", path], stdout=open(temp_bed_file, "w"), check=True)
    bed_files.append(temp_bed_file)

# Run bedtools multiintersect
cmd = ["bedtools", "multiinter", "-i"] + bed_files + ["-header", "-wa"]
subprocess.run(cmd, stdout=open(output_bed_file, "w"), check=True)

# Clean up temporary BED files
for file in bed_files:
    os.remove(file)
os.rmdir(temp_dir)

print("All stratifications have been merged into a single BED file:", output_bed_file)
