import sys
import csv
import argparse
from collections import defaultdict

import numpy as np
import pyrodigal
import Bio.SeqIO

import frames
import code

def nucleotid_count(sequence):
    result = {nt: 0 for nt in "ATCGN"}
    for nt in sequence:
        result[nt] += 1
    return result

def calculate_gc_content(sequence):
    nts = nucleotid_count(sequence)
    gc_count = nts['C'] + nts['G']
    return (gc_count / len(sequence)) if sequence else 0

def main():
    parser = argparse.ArgumentParser(
        prog='fasta_stats.py',
        description='print stats for a fasta file'
    )
    parser.add_argument('files', metavar='FILE', type=str, nargs='+', help="input file, default stdin")
    args = parser.parse_args()

    tsv_writer = csv.writer(sys.stdout, delimiter="\t")

    for file in args.files:
        row = []

        print(file, file=sys.stderr)
        print("reading fasta from " + file, file=sys.stderr)
        global sequences
        sequences = list(Bio.SeqIO.parse(file, "fasta"))

        # TODO extend to multi chromosomes
        if len(sequences) != 1:
            print("expecting exactly one sequence in fasta, found len(sequences)")
            sys.exit(1)

        row.append(sequences[0].id)

        global sequence
        sequence = str(sequences[0].seq)

        print("finding genes", file=sys.stderr)
        gene_finder = pyrodigal.GeneFinder()
        gene_finder.train(sequence)

        global genes
        genes = gene_finder.find_genes(sequence)

        print("calculating statistics", file=sys.stderr)

        sequence_lengths = [g.end - g.begin + 1 for g in genes]

        # Whole genome length
        genome_length = len(sequence)
        print(f"Length of Genome: {genome_length}", file=sys.stderr)
        row.append(genome_length)

        # Total Length of all Genes
        genes_total_length = sum(sequence_lengths)
        print(f"Total length of all genes: {genes_total_length} bases", file=sys.stderr)
        row.append(genes_total_length)

        # Number of Genes
        genes_length = len(genes)
        print(f"Number of genes: {genes_length}", file=sys.stderr)

        # Genes Length Statistics
        min_length = min(sequence_lengths)
        max_length = max(sequence_lengths)
        quartiles_length = np.percentile(sequence_lengths, [25, 50, 75])
        print(f"Length Min/Quartiles/Max: {min_length}, {quartiles_length[0]}, {quartiles_length[1]}, {quartiles_length[2]}, {max_length}", file=sys.stderr)
        row += [min_length, *quartiles_length, max_length]

        # GC Content
        genome_gc = calculate_gc_content(sequence)
        print(f"Whole Genome GC Conent: {genome_gc*100:0.4f}%", file=sys.stderr)
        row.append(genome_gc)

        gc_contens = [ calculate_gc_content(g.sequence()) for g in genes ]
        min_gc = min(gc_contens)
        max_gc = max(gc_contens)
        quartiles_gc = np.percentile(gc_contens, [25, 50, 75])
        print(f"Length Min/Quartiles/Max: {min_gc*100:0.4f}%, {quartiles_gc[0]*100:0.4f}%, {quartiles_gc[1]*100:0.4f}%, {quartiles_gc[2]*100:0.4f}%, {max_gc*100:0.4f}%", file=sys.stderr)
        row += [min_gc, *quartiles_gc, max_gc]

        # Coding Density
        coding_density = genes_total_length / (genome_length * 2) # Count Forward and Backward Strand
        print(f"Coding Density: {coding_density*100: 0.4}%", file=sys.stderr)
        row.append(coding_density)

        # Start Codonsstart_codons
        start_codons = defaultdict(int)
        for g in genes:
            start_codons[g.start_type] += 1
        
        row.append(start_codons.get("ATG", 0))
        row.append(start_codons.get("GTG", 0))
        row.append(start_codons.get("TTG", 0))
        row.append(start_codons.get("Edge", 0))
        start_codon_count = sum(start_codons.values())
        start_codons = sorted(list(start_codons.items()), key=lambda x: -x[1])
        start_codons = ", ".join(f"{codon}: {float(count)/start_codon_count*100:0.4f}%" for codon, count in start_codons)
        print(f"Start Codon Distribution: {start_codons}", file=sys.stderr)

        # codon usage in genome
        codon_counts = defaultdict(int)
        for seq in frames.frames(sequence):
            for (_, codon) in seq:
                codon_counts[codon] += 1
        codon_count = sum(codon_counts.values())
        codon_counts_row = sorted(codon_counts.items(), key=lambda x: code.codon_sort_key(x[0]))
        row += [ float(c) / codon_count for _, c in codon_counts_row]
        codon_counts_print = sorted(codon_counts.items(), key=lambda x: -x[1])
        codon_counts_print = ", ".join(f"{codon}: {float(count) / codon_count * 100:0.4f}%" for codon, count in codon_counts_print)
        print(f"Codon Distribution Whole Genome: {codon_counts_print}", file=sys.stderr)

        # codon usage in genes
        codon_counts = defaultdict(int)
        for gene in genes:
            for (_, codon) in frames.codons(gene.sequence()):
                codon_counts[codon] += 1
        codon_count = sum(codon_counts.values())
        codon_counts_row = sorted(codon_counts.items(), key=lambda x: code.codon_sort_key(x[0]))
        row += [ float(c) / codon_count for _, c in codon_counts_row]
        codon_counts_print = sorted(codon_counts.items(), key=lambda x: -x[1])
        codon_counts_print = ", ".join(f"{codon}: {float(count) / codon_count * 100:0.4f}%" for codon, count in codon_counts_print)
        print(f"Codon Distribution Genes: {codon_counts_print}", file=sys.stderr)

        tsv_writer.writerow(row)


if __name__ == "__main__":
    main()