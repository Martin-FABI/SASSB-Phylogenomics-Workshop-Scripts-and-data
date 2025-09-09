#!/usr/bin/env python3
# Martin Coetzee, University of Pretoria, developed this approach.
# This is a work in progress, and correct results are not guaranteed.
# The user should not distribute this script or publish it.
# The user should contact martin.coetzee@up.ac.za should he/she want to use it for their research.

import os
import numpy as np
from Bio import AlignIO
from Bio.Phylo.TreeConstruction import DistanceCalculator
import argparse
from scipy import stats
import matplotlib.pyplot as plt
from pathlib import Path

def parse_arguments():
    """Parse command line arguments with sensible defaults for testing."""
    parser = argparse.ArgumentParser(description='Detect outlier sequences in multiple sequence alignments.')
    parser.add_argument('--input_dir', default='.', help='Directory containing FASTA alignment files (default: current directory)')
    parser.add_argument('--output_dir', default='outlier_results', help='Directory to save results (default: outlier_results)')
    parser.add_argument('--z_threshold', type=float, default=3.0, 
                        help='Z-score threshold for outlier detection (default: 3.0)')
    parser.add_argument('--plot', action='store_true', help='Generate distance distribution plots')
    parser.add_argument('--show_plots', action='store_true', help='Display plots interactively')
    return parser.parse_args()

def load_alignment(file_path):
    """Load a multiple sequence alignment from a FASTA file."""
    try:
        # Auto-detect format - FASTA is most common but this allows flexibility
        alignment = AlignIO.read(file_path, "fasta")
        return alignment
    except Exception as e:
        print(f"Error loading alignment from {file_path}: {e}")
        return None

def calculate_distance_matrix(alignment, model="identity"):
    """Calculate the distance matrix for the alignment using specified model."""
    calculator = DistanceCalculator(model)
    try:
        # Calculate the distance matrix
        distance_matrix = calculator.get_distance(alignment)
        return distance_matrix
    except Exception as e:
        print(f"Error calculating distance matrix: {e}")
        return None

def analyze_distances(distance_matrix, sequence_ids, z_threshold=3.0):
    """
    Analyze the distance matrix to identify potential outlier sequences.
    Returns a list of (sequence_id, avg_distance, z_score) tuples for outliers.
    """
    n_sequences = len(sequence_ids)
    
    # Calculate average distance for each sequence to all others
    avg_distances = []
    for i in range(n_sequences):
        # Extract distances from one sequence to all others (excluding self-comparison)
        seq_distances = [distance_matrix[i, j] for j in range(n_sequences) if i != j]
        avg_distance = np.mean(seq_distances)
        avg_distances.append(avg_distance)
    
    # Calculate Z-scores for average distances
    z_scores = stats.zscore(avg_distances)
    
    # Identify outliers based on Z-score threshold
    outliers = []
    for i, (avg_dist, z_score) in enumerate(zip(avg_distances, z_scores)):
        if z_score > z_threshold:
            outliers.append((sequence_ids[i], avg_dist, z_score))
    
    return outliers, avg_distances, z_scores

def plot_distance_distribution(avg_distances, z_scores, sequence_ids, outliers, alignment_name, output_dir, show_plot=False):
    """Generate a plot showing the distribution of average distances and outliers."""
    plt.figure(figsize=(12, 8))
    
    # Scatter plot of average distances
    plt.scatter(range(len(avg_distances)), avg_distances, alpha=0.7)
    
    # Highlight outliers
    outlier_indices = [sequence_ids.index(outlier[0]) for outlier in outliers]
    if outlier_indices:
        plt.scatter(outlier_indices, [avg_distances[i] for i in outlier_indices], 
                   color='red', s=100, label='Outliers')
    
    # Add labels and other visual elements
    plt.axhline(y=np.mean(avg_distances), color='green', linestyle='-', 
               label=f'Mean distance: {np.mean(avg_distances):.4f}')
    plt.axhline(y=np.mean(avg_distances) + 3 * np.std(avg_distances), color='orange', 
               linestyle='--', label='3σ threshold')
    
    plt.xlabel('Sequence Index')
    plt.ylabel('Average Distance to Other Sequences')
    plt.title(f'Distance Distribution - {alignment_name}')
    plt.legend()
    plt.grid(True, alpha=0.3)
    
    # Add sequence IDs for outliers as annotations
    for i in outlier_indices:
        plt.annotate(sequence_ids[i], (i, avg_distances[i]), 
                    textcoords="offset points", xytext=(0,10), ha='center')
    
    # Add a small plot showing z-score distribution
    ax2 = plt.axes([0.15, 0.15, 0.3, 0.3])  # [left, bottom, width, height]
    ax2.hist(z_scores, bins=10, alpha=0.7)
    ax2.axvline(x=3.0, color='red', linestyle='--', label='Z=3')
    ax2.set_title('Z-score Distribution')
    
    # Save the plot
    plot_path = os.path.join(output_dir, f"{alignment_name}_distance_plot.png")
    plt.savefig(plot_path)
    
    # Show plot if requested
    if show_plot:
        plt.show()
    else:
        plt.close()
    
    return plot_path

def process_alignment_file(file_path, output_dir, z_threshold=3.0, generate_plot=False, show_plot=False):
    """Process a single alignment file and detect outliers."""
    alignment = load_alignment(file_path)
    if alignment is None:
        return None
    
    # Extract alignment information
    alignment_name = os.path.basename(file_path).split('.')[0]
    sequence_ids = [record.id for record in alignment]
    
    # Calculate distance matrix
    distance_matrix = calculate_distance_matrix(alignment)
    if distance_matrix is None:
        return None
    
    # Find potential outliers
    outliers, avg_distances, z_scores = analyze_distances(
        distance_matrix, sequence_ids, z_threshold)
    
    # Generate report
    results = {
        "alignment_name": alignment_name,
        "file_path": file_path,
        "num_sequences": len(alignment),
        "alignment_length": alignment.get_alignment_length(),
        "outliers": outliers,
        "plot_path": None
    }
    
    # Generate plot if requested
    if generate_plot:
        plot_path = plot_distance_distribution(
            avg_distances, z_scores, sequence_ids, outliers, 
            alignment_name, output_dir, show_plot)
        results["plot_path"] = plot_path
    
    return results

def write_report(results, output_dir):
    """Write a summary report of all processed alignments."""
    report_path = os.path.join(output_dir, "outlier_report.txt")
    with open(report_path, 'w') as f:
        f.write("OUTLIER SEQUENCE ANALYSIS REPORT\n")
        f.write("===============================\n\n")
        
        found_outliers = False
        for result in results:
            if result is None:
                continue
                
            f.write(f"Alignment: {result['alignment_name']}\n")
            f.write(f"File: {result['file_path']}\n")
            f.write(f"Sequences: {result['num_sequences']}\n")
            f.write(f"Alignment Length: {result['alignment_length']}\n")
            
            if result['outliers']:
                found_outliers = True
                f.write("\nPotential outlier sequences:\n")
                for seq_id, avg_dist, z_score in result['outliers']:
                    f.write(f"  - {seq_id}: Avg Distance = {avg_dist:.4f}, Z-score = {z_score:.4f}\n")
                
                if result['plot_path']:
                    f.write(f"\nDistance plot saved to: {result['plot_path']}\n")
            else:
                f.write("\nNo outlier sequences detected.\n")
            
            f.write("\n" + "-"*50 + "\n\n")
        
        if not found_outliers:
            f.write("No outlier sequences were detected in any alignments.\n")
    
    print(f"Report saved to {report_path}")
    return report_path

def setup_matplotlib():
    """Configure matplotlib to work in various environments."""
    # Try to set a backend that works in most environments
    try:
        # For systems with a GUI
        import matplotlib
        if os.environ.get('DISPLAY'):
            matplotlib.use('TkAgg')  # Good for Linux with display
        elif os.name == 'nt':  # Windows
            matplotlib.use('TkAgg')
        elif os.name == 'posix':  # macOS
            matplotlib.use('MacOSX')
        else:
            matplotlib.use('Agg')  # Fallback for headless systems
        
        print(f"Using matplotlib backend: {matplotlib.get_backend()}")
    except Exception as e:
        print(f"Warning: Could not set matplotlib backend: {e}")
        print("Plots will be saved but may not display interactively.")

def main():
    """Main function to process all alignment files in the current directory."""
    args = parse_arguments()
    
    # Setup matplotlib if we want to show plots
    if args.show_plots:
        setup_matplotlib()
    
    # Create output directory if it doesn't exist
    os.makedirs(args.output_dir, exist_ok=True)
    
    # Get list of FASTA files in current directory
    input_files = []
    for ext in ['.fasta', '.fa', '.faa', '.fna', '.aln', '.fas']:
        input_files.extend(list(Path(args.input_dir).glob(f'*{ext}')))
    
    if not input_files:
        print(f"No FASTA files found in {args.input_dir}")
        return
    
    print(f"Found {len(input_files)} alignment files to process.")
    
    # Process each alignment file
    results = []
    for file_path in input_files:
        print(f"Processing {file_path}...")
        result = process_alignment_file(
            file_path, args.output_dir, 
            z_threshold=args.z_threshold, 
            generate_plot=args.plot,
            show_plot=args.show_plots
        )
        results.append(result)
    
    # Write summary report
    write_report(results, args.output_dir)
    print("Analysis complete!")
    
    # Print instructions for viewing saved plots
    if args.plot and not args.show_plots:
        print("\nPlots have been saved to the output directory but were not displayed.")
        print(f"You can find them in: {os.path.abspath(args.output_dir)}")
        print("To view plots interactively during analysis, run with --show_plots flag")

if __name__ == "__main__":
    main()
