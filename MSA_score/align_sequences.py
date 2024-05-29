import matplotlib.pyplot as plt
import numpy as np
from Bio import SeqIO, AlignIO

def parse_blast_results(blast_file):
    """Parses a BLAST results file and returns a dictionary with (query_id, subject_id) as keys
       and (alignment_array, sequence, query_start, query_end) as values.
    """
    blast_dict = {}
    with open(blast_file, "r") as file:
        for line in file:
            if line.startswith("#"):  # Skip comment lines
                continue
            parts = line.strip().split("\t")
            query_id = parts[0]
            subject_id = parts[1]
            query_start = int(parts[3]) - 1 # Convert to 0-based index
            query_end = int(parts[4]) - 1  # Convert to 0-based index
            sequence = parts[5]  # Assuming the aligned sequence is in this column
            reference = parts[6]

            # Create an empty alignment array of the length of the query alignment
            alignment_length = query_end - query_start + 1
            alignment_array = [""] * alignment_length

            key = (query_id, subject_id)
            value = (alignment_array, sequence, query_start, query_end, reference)
            blast_dict[key] = value
    return blast_dict

def compare_sequences(blast_results):
    """Compares sequences from BLAST results with the reference sequences and fills the alignment array."""
    for (query_id, subject_id), (alignment_array, query_sequence, query_start, query_end, reference_sequence) in blast_results.items():
        ref_index = 0
        query_index = 0

        while query_index < len(query_sequence) and ref_index < len(reference_sequence):
            if query_sequence[query_index] == reference_sequence[ref_index]:
                alignment_array[ref_index] = query_sequence[query_index]
                query_index += 1
            else:
                alignment_array[ref_index] = '-'
            ref_index += 1


def extract_alignment_arrays(blast_results):
    """Creates a new dictionary with (query_id, subject_id) as keys and alignment_array as values."""
    alignment_dict = {}
    for key, (alignment_array, _, _, _, _) in blast_results.items():
        alignment_dict[key] = alignment_array
    return alignment_dict

def parse_msa(msa_file):
    """Parses a MSA file and returns a dictionary with subject IDs as keys and arrays of bases as values."""
    alignment = AlignIO.read(msa_file, "fasta")
    msa_dict = {}
    
    for record in alignment:
        subject_id = record.id
        sequence = str(record.seq).upper()
        sequence_array = list(sequence)
        match_count_array = [0] * len(sequence)  # Initialize the match count array with zeros
        msa_dict[subject_id] = (sequence_array, match_count_array)
    
    return msa_dict

def update_match_count(blast_results, msa_dict):
    """Updates the match count array in the msa_dict based on the comparison with BLAST results."""
    for (query_id, subject_id), (alignment_array, query_sequence, query_start, query_end, reference) in blast_results.items():
        if subject_id not in msa_dict:
            print(f"Subject ID {subject_id} not found in MSA.")
            continue

        msa_sequence_array, match_count_array = msa_dict[subject_id]

        # Find the start position in the MSA sequence corresponding to the query start position
        msa_index = 0
        non_gap_count = 0
        while non_gap_count < query_start and msa_index < len(msa_sequence_array):
            if msa_sequence_array[msa_index] != '-':
                non_gap_count += 1
            msa_index += 1

        # Now msa_index is at the start position of the query sequence in the MSA sequence
        #ref_index = query_start
        query_index = 0
        while msa_index < len(msa_sequence_array) and query_index < len(query_sequence):
            # Increment match count if characters match or if both are gaps
            if query_sequence[query_index] == msa_sequence_array[msa_index]:
                match_count_array[msa_index] += 1
                query_index += 1
            msa_index += 1
            #ref_index += 1
    
def extract_match_counts(msa_dict):
    """Extracts just the keys and the match count arrays from the msa_dict."""
    match_counts = {}
    for subject_id, (_, match_count_array) in msa_dict.items():
        match_counts[subject_id] = match_count_array
    return match_counts

def create_merged_counts(msa_dict):
    """Creates a merged counts array that accumulates the counts from the MSA dictionary."""
    msa_length = len(list(msa_dict.values())[0][1]) # gets the number of bases in the MSA from the first ref sequence in the MSA
    merged_counts = [0] * msa_length

    for _, (_, match_count_array) in msa_dict.items():
        for i in range(msa_length):
            merged_counts[i] += match_count_array[i]
    
    return merged_counts

def plot_counts(merged_counts, output_file):
    """Plots the merged counts and saves the plot to a file."""
    x = np.arange(len(merged_counts))
    
    plt.figure(figsize=(10, 6))
    plt.bar(x, merged_counts, color='skyblue', edgecolor='black')
    plt.title('Bar Plot of Merged Counts')
    plt.xlabel('Index')
    plt.ylabel('Count')
    
    plt.savefig(output_file)
    plt.close()


def main(blast_file, msa_file, output_file):
    blast_results = parse_blast_results(blast_file)
    compare_sequences(blast_results)
    alignment_dict = extract_alignment_arrays(blast_results)
    msa_sequences = parse_msa(msa_file)
    update_match_count(blast_results, msa_sequences)
    msa_matches = extract_match_counts(msa_sequences)
    merged_counts = create_merged_counts(msa_sequences)

    # Testing the results 
    print("\nAlignment Arrays:")
    for key, alignment_array in alignment_dict.items():
        print(f"{key}: {alignment_array}")

    print("\nMSA Results:")
    for key, value in msa_matches.items():
        print(f"{key}: {value}")

    print("\nMerged Counts:")
    print(merged_counts)

    # Generate and save the plot
    plot_counts(merged_counts, output_file)
    print(f"Plot saved to {output_file}")


if __name__ == "__main__":
    import argparse
    
    parser = argparse.ArgumentParser(description="Parse BLAST and MSA result files.")
    parser.add_argument("blast_file", help="Input BLAST file")
    parser.add_argument("msa_file", help="Input MSA file")
    parser.add_argument("output_file", help="Output file for the plot")

    
    args = parser.parse_args()
    
    main(args.blast_file, args.msa_file, args.output_file)