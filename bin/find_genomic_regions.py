#!/usr/bin/env python3

import argparse
import gffutils
from collections import defaultdict
from intervaltree import Interval, IntervalTree
import os
import sys
import gzip
import matplotlib.pyplot as plt
import numpy as np

def parse_arguments():
    parser = argparse.ArgumentParser(
        description='Find genomic regions from VCF variants'
    )
    parser.add_argument(
        '--vcf',
        required=True,
        help='Input VCF file (can be gzipped)'
    )
    parser.add_argument(
        '--gtf',
        required=True,
        help='GTF annotation file'
    )
    parser.add_argument(
        '--fai',
        required=True,
        help='FASTA index (.fai) file'
    )
    parser.add_argument(
        '--cancer_type',
        required=True,
        help='Cancer type identifier'
    )
    parser.add_argument(
        '--output',
        required=True,
        help='Output TSV file for variant regions'
    )
    parser.add_argument(
        '--summary',
        required=True,
        help='Output TSV file for summary statistics'
    )
    parser.add_argument(
        '--plot',
        required=True,
        help='Output PNG file for density plot'
    )
    parser.add_argument(
        '--db_name',
        default='gencode_grch38.db',
        help='Name for gffutils database file (default: gencode_grch38.db)'
    )
    
    return parser.parse_args()

def get_genome_length_from_fai(fai_path):
    """
    Calculates total genome length by summing chromosome lengths from a .fai file and
    normalizes chromosome names to "chrX" format.
    """
    total_length = 0
    print(f"Reading genome lengths from .fai file: {fai_path}")
    try:
        with open(fai_path, 'r') as f:
            for line in f:
                parts = line.strip().split('\t')
                if len(parts) >= 2:
                    chrom_name = parts[0]
                    # Normalize chrom name from FAI to match "chrX" if needed, consistent with GTF/VCF processing
                    if not chrom_name.startswith("chr") and chrom_name.lower() != "mt": # Handle 'MT' vs 'chrM'
                        chrom_name = "chr" + chrom_name
                    elif chrom_name.lower() == "mt":
                        chrom_name = "chrM"
                        
                    # Only sum standard chromosomes (1-22, X, Y, M) to avoid unplaced contigs if desired
                    if chrom_name.startswith("chr") and (chrom_name[3:].isdigit() or chrom_name[3:] in ['X', 'Y', 'M']):
                         total_length += int(parts[1])
    except FileNotFoundError:
        print(f"Error: .fai file not found at {fai_path}. Please check the path.", file=sys.stderr)
        sys.exit(1)
    except Exception as e:
        print(f"Error reading .fai file: {e}", file=sys.stderr)
        sys.exit(1)
    print(f"Total genome length from .fai: {total_length} bp")
    return total_length

def classify_variant_region(chrom, pos_1based, cds_trees, utr5_trees, utr3_trees, exon_trees, gene_trees):
    """
    Classify a variant into genomic regions based on its 1-based position.
    Priority order: CDS > 5'UTR > 3'UTR > Exon > Intron > Intergenic
    """
    # Convert 1-based VCF position to 0-based half-open interval for IntervalTree query
    # A single base at 1-based 'pos_1based' is represented as [pos_1based - 1, pos_1based)
    query_start_0based = pos_1based - 1
    query_end_0based = pos_1based # End is exclusive in intervaltree

    # Check for overlaps in priority order
    in_cds = chrom in cds_trees and cds_trees[chrom].overlaps(query_start_0based, query_end_0based)
    in_5utr = chrom in utr5_trees and utr5_trees[chrom].overlaps(query_start_0based, query_end_0based)
    in_3utr = chrom in utr3_trees and utr3_trees[chrom].overlaps(query_start_0based, query_end_0based)
    in_exon = chrom in exon_trees and exon_trees[chrom].overlaps(query_start_0based, query_end_0based)
    in_gene = chrom in gene_trees and gene_trees[chrom].overlaps(query_start_0based, query_end_0based)

    if in_cds:
        return "CDS"
    elif in_5utr:
        return "5'UTR"
    elif in_3utr:
        return "3'UTR"
    elif in_exon:
        return "Exon"
    elif in_gene: # If it's in a gene but not in CDS/UTRs/Exon, it's an intron
        return "Intron"
    else: # If it's not in any gene, it's intergenic
        return "Intergenic"

def accumulate_lengths(trees, label, region_lengths):
    """Accumulate non-redundant lengths for a given feature type."""
    for chrom, tree in trees.items():
        # Using tree.merge_overlaps() provides a set of non-overlapping intervals
        # then sum their lengths.
        # This is crucial for accurate region lengths, especially for exons/UTRs/CDS.
        merged_intervals = tree.copy() # Make a copy to avoid modifying original tree
        merged_intervals.merge_overlaps() # Merge intervals in place within the copy
        for interval in merged_intervals:
            region_lengths[label] += interval.end - interval.begin

def main():
    args = parse_arguments()
    
    # Assign arguments to variables
    vcf_file = args.vcf
    gtf_file = args.gtf
    fasta_fai_file = args.fai
    cancer_type = args.cancer_type
    output_file = args.output
    summary_output_file = args.summary
    output_plot_path = args.plot
    db_name = args.db_name
    
    print(f"Processing cancer type: {cancer_type}")
    print(f"VCF file: {vcf_file}")
    print(f"GTF file: {gtf_file}")
    print(f"FAI file: {fasta_fai_file}")
    
    # Calculate total genome length once at the beginning
    TOTAL_GENOME_LENGTH = get_genome_length_from_fai(fasta_fai_file)
    if TOTAL_GENOME_LENGTH == 0:
        print("Warning: Total genome length from FAI is 0. This might lead to division by zero or incorrect normalization. Exiting.", file=sys.stderr)
        sys.exit(1)
    
    # Build GTF DB 
    print(f"Building gffutils database from {gtf_file}...")
    try:
        if not os.path.exists(db_name):
            db = gffutils.create_db(
                gtf_file,
                dbfn=db_name,
                force=False,
                disable_infer_transcripts=True,
                disable_infer_genes=True,
                keep_order=True
            )
            print("Database built successfully.")
        else:
            db = gffutils.FeatureDB(db_name)
            print("Using existing gffutils database.")
    except Exception as e:
        print(f"Error building/loading gffutils DB: {e}", file=sys.stderr)
        sys.exit(1)
    
    # Build Interval Trees 
    utr5_trees = defaultdict(IntervalTree)
    utr3_trees = defaultdict(IntervalTree)
    exon_trees = defaultdict(IntervalTree)
    cds_trees = defaultdict(IntervalTree)
    gene_trees = defaultdict(IntervalTree) # To identify introns (within gene, not exon/UTR/CDS)
    
    print("Loading features into IntervalTrees...")
    try:
        for feature in db.all_features():
            chrom = feature.seqid
            start = feature.start
            end = feature.end
    
            # Normalize chromosome names to "chrX" format if not already
            if chrom.lower() == "mt":
                chrom = "chrM"
            elif not chrom.startswith("chr"):
                chrom = "chr" + chrom    
    
            # Convert 1-based inclusive GTF coordinates to 0-based half-open for IntervalTree
            # [start, end] (1-based inclusive) becomes [start-1, end) (0-based half-open)
            interval = Interval(start - 1, end)
    
            if feature.featuretype == 'five_prime_utr':
                utr5_trees[chrom].add(interval)
            elif feature.featuretype == 'three_prime_utr':
                utr3_trees[chrom].add(interval)
            elif feature.featuretype == 'exon':
                exon_trees[chrom].add(interval)
            elif feature.featuretype == 'CDS':
                cds_trees[chrom].add(interval)
            elif feature.featuretype == 'gene':
                gene_trees[chrom].add(interval)
    
        print("Feature loading complete.")
    except Exception as e:
        print(f"Error loading features: {e}", file=sys.stderr)
        sys.exit(1)
    
    # Process VCF 
    print(f"Processing VCF file {vcf_file} and writing to {output_file}...")
    
    region_counts = defaultdict(int)
    
    try:
        # Use gzip.open in 'rt' mode to read gzipped file as text
        with gzip.open(vcf_file, "rt") as vcf_in, open(output_file, "w") as out_fh:
            out_fh.write("CHROM\tPOS\tREF\tALT\tREGION\n")
    
            for line in vcf_in:
                if line.startswith('#'): # Skip header and meta-info lines
                    continue
                
                parts = line.strip().split('\t')
                # Ensure line has enough columns before accessing indices
                if len(parts) < 5: 
                    continue
    
                chrom = parts[0]
                pos_1based = int(parts[1]) # VCF POS is 1-based
                ref = parts[3]
                alt = parts[4]
    
                # Normalize chromosome names for VCF variants too
                if chrom.lower() == "mt":
                    chrom = "chrM"
                elif not chrom.startswith("chr"):
                    chrom = "chr" + chrom
    
                region = classify_variant_region(chrom, pos_1based, cds_trees, utr5_trees, utr3_trees, exon_trees, gene_trees)
                region_counts[region] += 1
                out_fh.write(f"{chrom}\t{pos_1based}\t{ref}\t{alt}\t{region}\n")
    
    except Exception as e:
        print(f"VCF processing error: {e}", file=sys.stderr)
        sys.exit(1)
    
    # Compute Region Lengths
    print("Computing region lengths for normalization...")
    
    region_lengths = defaultdict(int)
    
    accumulate_lengths(cds_trees, "CDS", region_lengths)
    accumulate_lengths(utr5_trees, "5'UTR", region_lengths)
    accumulate_lengths(utr3_trees, "3'UTR", region_lengths)
    accumulate_lengths(exon_trees, "Exon", region_lengths)
    accumulate_lengths(gene_trees, "Gene", region_lengths) # Accumulate total gene length
    
    # Calculate total annotated gene length for intergenic calculation
    total_gene_annotated_length = region_lengths.get("Gene", 0)
    
    # Intergenic length = Total Genome Length - Total Annotated Gene Length
    region_lengths["Intergenic"] = TOTAL_GENOME_LENGTH - total_gene_annotated_length
    if region_lengths["Intergenic"] < 0:
        print(f"Warning: Calculated intergenic length was negative ({region_lengths['Intergenic']}). Setting to 0. "
              "This might indicate an issue with genome length or annotation coverage.", file=sys.stderr)
        region_lengths["Intergenic"] = 0
    
    # Calculate the TRUE Intron length
    total_exonic_length = 0
    for chrom, tree in exon_trees.items():
        merged_exons = tree.copy()
        merged_exons.merge_overlaps()
        for interval in merged_exons:
            total_exonic_length += interval.end - interval.begin
    
    # The total intron length is the total gene length minus the total exonic length
    region_lengths["Intron"] = total_gene_annotated_length - total_exonic_length
    
    # Handle cases where this might be negative due to annotation inconsistencies or overlapping genes
    if region_lengths["Intron"] < 0:
        print(f"Warning: Calculated intron length was negative ({region_lengths['Intron']}). Setting to 0. "
              "This might indicate annotation overlaps or inconsistencies.", file=sys.stderr)
        region_lengths["Intron"] = 0
    
    # Generate Summary TSV 
    print(f"Generating summary file: {summary_output_file}...")
    try:
        with open(summary_output_file, "w") as sum_fh:
            sum_fh.write("REGION\tNUMBER OF VARIANTS\tTOTAL REGION LENGTH (bp)\n")
            
            # Ensure consistent order for summary as well
            regions_order = ["CDS", "5'UTR", "3'UTR", "Exon", "Intron", "Intergenic"]
            
            for region in regions_order:
                count = region_counts.get(region, 0)
                length = region_lengths.get(region, 0)
                sum_fh.write(f"{region}\t{count}\t{length}\n")
        print("Summary file generated successfully.")
    except Exception as e:
        print(f"Error generating summary file: {e}", file=sys.stderr)
        sys.exit(1)
    
    # Normalization and Plotting 
    regions_order = ["CDS", "5'UTR", "3'UTR", "Exon", "Intron", "Intergenic"]
    norm_vals = []
    
    # Total variant count for plot title
    total_variant_count = sum(region_counts.get(region, 0) for region in regions_order)
    
    print("\nNormalized Variant Densities:")
    for region in regions_order:
        count = region_counts.get(region, 0)
        length = region_lengths.get(region, 0) 
        
        if length <= 0:
            norm = 0.0
            print(f"  {region}: Count = {count}, Length = {length} (cannot normalize, length is zero or negative). Density = {norm:.2e} per kb")
        else:
            norm = count / (length/1000)
            print(f"  {region}: Count = {count}, Length = {length}, Normalized = {norm:.2e} per kb")
        
        norm_vals.append(norm)
    
    # Assign different colors on region types
    region_colors = {
        "CDS": "#e41a1c",         # red
        "5'UTR": "#ff7f00",       # orange
        "3'UTR": "#ff7f00",       # orange
        "Exon": "#377eb8",        # blue
        "Intron": "#999999",      # gray
        "Intergenic": "#4daf4a"   # green
    }
    bar_colors = [region_colors[r] for r in regions_order]
    
    plt.figure(figsize=(10, 6))
    bars = plt.bar(regions_order, norm_vals, color=bar_colors)
    
    plt.ylabel("Normalized Variant Density (per bp, log scale)")
    plt.xlabel("Genomic Region")
    plt.title(f"Normalized Variant Density in {cancer_type} (n={total_variant_count} variants)")
    plt.yscale("log")
    plt.grid(True, axis="y", linestyle="--", alpha=0.7)
    plt.tight_layout()
    
    plt.savefig(output_plot_path)
    print("Completed.")
    print(f"\nPlot saved to: {output_plot_path}")
    print(f"Output written to: {output_file}")
    print(f"Summary written to: {summary_output_file}")

if __name__ == "__main__":
    main()

