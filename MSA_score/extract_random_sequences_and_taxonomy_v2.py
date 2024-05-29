import random
from Bio import SeqIO

def pick_random_sequences(input_fasta, num_sequences):
    sequences = list(SeqIO.parse(input_fasta, "fasta"))
    
    if len(sequences) < num_sequences:
        raise ValueError(f"Input file only contains {len(sequences)} sequences, but {num_sequences} sequences are requested.")
    
    random_sequences = random.sample(sequences, num_sequences)
    return random_sequences

def write_sequences_to_fasta(sequences, output_fasta):
    with open(output_fasta, "w") as output_handle:
        SeqIO.write(sequences, output_handle, "fasta")

def extract_taxonomy_data(taxonomy_file, sequence_fasta, output_taxonomy_file):
    sequence_ids = {record.id for record in SeqIO.parse(sequence_fasta, "fasta")}
    taxonomy_data = {}
    
    with open(taxonomy_file, "r") as tax_file:
        for line in tax_file:
            parts = line.strip().split("\t")
            seq_id = parts[0]
            if seq_id in sequence_ids:
                taxonomy_data[seq_id] = line.strip()
    
    with open(output_taxonomy_file, "w") as out_file:
        for seq_id in sequence_ids:
            if seq_id in taxonomy_data:
                out_file.write(taxonomy_data[seq_id] + "\n")

def main(input_fasta, taxonomy_file, output_fasta, output_taxonomy_file, num_sequences):
    random_sequences = pick_random_sequences(input_fasta, num_sequences)
    
    write_sequences_to_fasta(random_sequences, output_fasta)
    extract_taxonomy_data(taxonomy_file, output_fasta, output_taxonomy_file)
    
    print(f"Successfully wrote {num_sequences} sequences to {output_fasta} and corresponding taxonomy data to {output_taxonomy_file}")

if __name__ == "__main__":
    import argparse
    
    parser = argparse.ArgumentParser(description="Extract random sequences and corresponding taxonomy data.")
    parser.add_argument("input_fasta", help="Input FASTA file")
    parser.add_argument("taxonomy_file", help="Taxonomy file")
    parser.add_argument("output_fasta", help="Output FASTA file for random sequences")
    parser.add_argument("output_taxonomy_file", help="Output file for corresponding taxonomy data")
    parser.add_argument("num_sequences", type=int, help="Number of sequences to extract")
    
    args = parser.parse_args()
    
    main(args.input_fasta, args.taxonomy_file, args.output_fasta, args.output_taxonomy_file, args.num_sequences)
