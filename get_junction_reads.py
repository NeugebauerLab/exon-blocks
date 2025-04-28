
import pysam
import pandas as pd
import os
import argparse

def parse_args():
    parser = argparse.ArgumentParser(description="Calculate junction coverage (RPM) and percentage from BAM files.")
    parser.add_argument("--bam_dir", type=str, required=True, help="Directory containing BAM files.")
    parser.add_argument("--junctions_csv", type=str, required=True, help="CSV file with junction coordinates.")
    parser.add_argument("--output_csv", type=str, required=True, help="Output CSV file for results.")
    parser.add_argument("--tolerance", type=int, default=5, help="Tolerance for junction matching (default: 5bp).")
    return parser.parse_args()

def adjust_chrom_name(chrom, bam_contigs):
    """
    Adjust chromosome name to match the BAM file format.
    - If 'chrX', 'chrY', 'chrM' are in CSV but BAM uses 'X', 'Y', 'M', remove 'chr'.
    - Keep 'chr1', 'chr2', etc., unchanged.
    """
    if chrom.startswith("chr"):
        chrom_base = chrom[3:]  # Remove "chr" prefix
        if chrom_base in {"X", "Y", "M"} and chrom_base in bam_contigs:
            return chrom_base  # Convert "chrX" → "X", "chrY" → "Y", "chrM" → "M"
    return chrom  # Keep "chr1", "chr2", etc.

def count_junction_reads(bam_file, junctions, tolerance):
    """
    Count reads spanning each junction in a BAM file with a tolerance.
    Tracks:
    - Total reads spanning the junction region
    - Reads containing the junction (`N` operation in CIGAR)
    """
    bam = pysam.AlignmentFile(bam_file, "rb")
    
    # Get the list of reference contigs in the BAM file
    bam_contigs = set(bam.references)

    junction_counts = []  # Reads containing the junction (N in CIGAR)
    total_counts = []     # Reads spanning the junction region

    for _, junction in junctions.iterrows():
        chrom = str(junction["chr"])  # Ensure it's a string

        # Adjust naming convention for chromosomes that do not have "chr" prefix in BAM
        if chrom.startswith("chr"):
            chrom_base = chrom[3:]  # Remove "chr" prefix
            if chrom_base in {"M", "X", "Y"} and chrom_base in bam_contigs:
                chrom = chrom_base  # Convert "chrM" -> "M", "chrX" -> "X", "chrY" -> "Y"
            elif chrom not in bam_contigs and chrom_base in bam_contigs:
                chrom = chrom_base  # Convert "chr1" -> "1" if needed

        # Check if chromosome exists in BAM file
        if chrom not in bam_contigs:
            print(f"Warning: Chromosome {chrom} not found in BAM. Skipping.")
            junction_counts.append(0)
            total_counts.append(0)
            continue  # Skip this junction if chromosome isn't found

        start = junction["end_ups_exon"]
        end = junction["start_dns_exon"]

        junction_count = 0
        total_count = 0

        for read in bam.fetch(chrom, start - tolerance, end + tolerance):
            if read.cigartuples:
                ref_start = read.reference_start
                spans_junction = False
                contains_junction = False
                
                for op, length in read.cigartuples:
                    if op == 3:  # 'N' operation (splicing)
                        junction_start = ref_start
                        junction_end = ref_start + length
                        if abs(junction_start - start) <= tolerance and abs(junction_end - end) <= tolerance:
                            contains_junction = True
                    if op in {0, 2, 3}:
                        ref_start += length

                total_count += 1  # Any read spanning the junction region is counted
                if contains_junction:
                    junction_count += 1  # Only reads with an 'N' in the CIGAR string are counted here

        junction_counts.append(junction_count)
        total_counts.append(total_count)

    bam.close()
    return junction_counts, total_counts


def calculate_rpm(counts, total_reads):
    """
    Convert raw counts to Reads Per Million (RPM).
    """
    return [(count / total_reads) * 1e6 if total_reads > 0 else 0 for count in counts]

def calculate_junction_percent(junction_counts, total_counts):
    """
    Calculate the percentage of reads spanning the junction region that contain the junction.
    """
    return [(junc / total * 100 if total > 0 else 0) for junc, total in zip(junction_counts, total_counts)]

def main():
    args = parse_args()
    
    # Load junction coordinates
    junctions = pd.read_csv(args.junctions_csv)
    
    # BAM files
    bam_files = sorted([os.path.join(args.bam_dir, f) for f in os.listdir(args.bam_dir) if f.endswith(".bam")])
    bam_files_kd = [f for f in bam_files if "dano" in f]
    bam_files_ctr = [f for f in bam_files if "dmso" in f]

    # Prepare results
    results = junctions.copy()
    total_reads = {}

    for bam_file in bam_files:
        # Get total reads in the BAM file
        bam = pysam.AlignmentFile(bam_file, "rb")
        total_reads[bam_file] = sum(1 for _ in bam.fetch(until_eof=True))
        bam.close()

    # Process KD samples
    kd_junction_percentages = []
    for i, bam_file in enumerate(bam_files_kd):
        print(f"Processing KD sample: {bam_file}")
        junction_counts, total_counts = count_junction_reads(bam_file, junctions, args.tolerance)
        rpms = calculate_rpm(junction_counts, total_reads[bam_file])
        junction_percentages = calculate_junction_percent(junction_counts, total_counts)
        
        results[f"kd_sample_{i + 1}"] = rpms
        results[f"kd_junction_percent_{i + 1}"] = junction_percentages
        kd_junction_percentages.append(junction_percentages)

    # Process CTR samples
    ctr_junction_percentages = []
    for i, bam_file in enumerate(bam_files_ctr):
        print(f"Processing CTR sample: {bam_file}")
        junction_counts, total_counts = count_junction_reads(bam_file, junctions, args.tolerance)
        rpms = calculate_rpm(junction_counts, total_reads[bam_file])
        junction_percentages = calculate_junction_percent(junction_counts, total_counts)

        results[f"ctr_sample_{i + 1}"] = rpms
        results[f"ctr_junction_percent_{i + 1}"] = junction_percentages
        ctr_junction_percentages.append(junction_percentages)

    # Calculate mean RPMs
    kd_columns = [f"kd_sample_{i + 1}" for i in range(len(bam_files_kd))]
    ctr_columns = [f"ctr_sample_{i + 1}" for i in range(len(bam_files_ctr))]
    results["mean_kd_rpm"] = results[kd_columns].mean(axis=1)
    results["mean_ctr_rpm"] = results[ctr_columns].mean(axis=1)

    # Calculate mean junction percentages
    kd_percent_columns = [f"kd_junction_percent_{i + 1}" for i in range(len(bam_files_kd))]
    ctr_percent_columns = [f"ctr_junction_percent_{i + 1}" for i in range(len(bam_files_ctr))]
    results["mean_percent_jr_kd"] = results[kd_percent_columns].mean(axis=1)
    results["mean_percent_jr_ctr"] = results[ctr_percent_columns].mean(axis=1)

    # Save results
    results.to_csv(args.output_csv, index=False)
    print(f"Results saved to {args.output_csv}")

if __name__ == "__main__":
    main()


