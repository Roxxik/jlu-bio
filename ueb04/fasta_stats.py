from code import code, codon_sort_key

import argparse
import numpy as np
import matplotlib.pyplot as plt

def read_fasta(file_path):
  with open(file_path, 'r') as file:
    sequences = {None: ""}
    header = None
    for line in file:
      line = line.strip()
      if line[0] == ">":
        header = line[1:].split(" ")[0] # remove > and only take the id (until first space)
        sequences[header] = ""
      else:
        sequences[header] += line
  if len(sequences[None]) != 0:
    print("Warning: sequence without header")
  del sequences[None]
  return sequences

def calculate_gc_content(sequence):
  gc_count = sum(1 for nt in sequence if nt in ['G', 'C'])
  return (gc_count / len(sequence)) * 100 if sequence else 0

def print_fasta_statistics(fasta_dict):
  sequence_lengths = [len(seq) for seq in fasta_dict.values()]

  # Number of sequences
  num_sequences = len(fasta_dict)
  print(f"Number of sequences: {num_sequences}")

  # Total size
  total_size = sum(sequence_lengths)
  print(f"Total size of all sequences: {total_size} bases")

  # Size statistics
  min_size = min(sequence_lengths)
  max_size = max(sequence_lengths)
  quartiles = np.percentile(sequence_lengths, [25, 50, 75])
  print(f"Min sequence size: {min_size}")
  print(f"Quartiles of sequence size: {quartiles}")
  print(f"Max sequence size: {max_size}")

  # GC content
  combined_gc_content = calculate_gc_content(''.join(fasta_dict.values()))
  print(f"GC content of all sequences combined: {combined_gc_content:.2f}%")

  individual_gc_contents = [calculate_gc_content(seq) for seq in fasta_dict.values()]
  min_gc = min(individual_gc_contents)
  max_gc = max(individual_gc_contents)
  quartiles_gc = np.percentile(individual_gc_contents, [25, 50, 75])
  print(f"Min GC content: {min_gc:.2f}%")
  print(f"Quartiles of GC content: {quartiles_gc}")
  print(f"Max GC content: {max_gc:.2f}%")

  return sequence_lengths, individual_gc_contents

def plot_box_plots(sequence_lengths, gc_contents):
  fig, axs = plt.subplots(1, 2, figsize=(12, 6))

  # Box plot for sequence lengths
  axs[0].boxplot(sequence_lengths)
  axs[0].set_title("Distribution of Sequence Lengths")
  axs[0].set_ylabel("Length")

  # Box plot for GC content
  axs[1].boxplot(gc_contents)
  axs[1].set_title("Distribution of GC Content")
  axs[1].set_ylabel("GC Content (%)")

  plt.tight_layout()
  plt.show()

def plot_box_plots_with_stats(sequence_lengths, gc_contents):
  fig, axs = plt.subplots(1, 2, figsize=(12, 6))

  # Function to add statistics text
  def add_stats(ax, data):
    stats = np.percentile(data, [25, 50, 75])
    whiskers = np.percentile(data, [5, 95])
    ax.text(1.1, max(data), f"Max: {max(data):.2f}", verticalalignment='top')
    ax.text(1.1, stats[2], f"Q3: {stats[2]:.2f}", verticalalignment='center')
    ax.text(1.1, stats[1], f"Median: {stats[1]:.2f}", verticalalignment='center')
    ax.text(1.1, stats[0], f"Q1: {stats[0]:.2f}", verticalalignment='center')
    ax.text(1.1, min(data), f"Min: {min(data):.2f}", verticalalignment='bottom')
    ax.text(1.1, whiskers[1], f"95%: {whiskers[1]:.2f}", verticalalignment='center')
    ax.text(1.1, whiskers[0], f"5%: {whiskers[0]:.2f}", verticalalignment='center')

  # Box plot for sequence lengths
  axs[0].boxplot(sequence_lengths, showmeans=True)
  axs[0].set_yscale('log')
  axs[0].set_title("Distribution of Sequence Lengths")
  axs[0].set_ylabel("Length")
  add_stats(axs[0], sequence_lengths)

  # Box plot for GC content
  axs[1].boxplot(gc_contents, showmeans=True)
  axs[1].set_title("Distribution of GC Content")
  axs[1].set_ylabel("GC Content (%)")
  add_stats(axs[1], gc_contents)

  plt.tight_layout()
  plt.show()

def plot_correlation(sequence_lengths, gc_contents):
  # Calculate Pearson correlation coefficient
  correlation_coefficient = np.corrcoef(sequence_lengths, gc_contents)[0, 1]

  plt.figure(figsize=(8, 6))
  plt.scatter(sequence_lengths, gc_contents, s=1)
  plt.xscale('log')  # Log scale for sequence length if needed
  plt.xlabel('Sequence Length (Log Scale)')
  plt.ylabel('GC Content (%)')
  plt.title(f'Correlation between Sequence Length and GC Content\nPearson Correlation: {correlation_coefficient:.2f}')
  plt.grid(True)
  plt.show()

def calculate_codon_frequency(sequences):
    codon_frequency = {}
    for seq in sequences:
        for i in range(0, len(seq) - 2, 3):
            codon = seq[i:i+3]
            if codon in codon_frequency:
                codon_frequency[codon] += 1
            else:
                codon_frequency[codon] = 1
    return codon_frequency

def plot_codon_frequency(codon_freq):
   # Sorting codons for better visualization
  sorted_codons = sorted(codon_freq.keys())
  sorted_freqs = [codon_freq[codon] for codon in sorted_codons]

  plt.figure(figsize=(15, 5))
  plt.bar(sorted_codons, sorted_freqs, color='skyblue')
  plt.xlabel('Codon')
  plt.ylabel('Frequency')
  plt.title('Codon Usage Frequency')
  plt.xticks(rotation=90)  # Rotate labels for better readability
  plt.show()

# Example data
sequence_lengths = [100, 200, 150, 120, 180, 220, 140, 160]
gc_contents = [40, 50, 45, 60, 55, 65, 42, 58]

def group_codon_frequency_by_amino_acid(sequences):
    aa_codon_frequency = {}
    for seq in sequences:
        for i in range(0, len(seq) - 2, 3):
            codon = seq[i:i+3]
            aa = code.get(codon, None)  # Get the amino acid encoded by the codon
            if aa:
                if aa not in aa_codon_frequency:
                    aa_codon_frequency[aa] = {}
                aa_codon_frequency[aa][codon] = aa_codon_frequency[aa].get(codon, 0) + 1
    return aa_codon_frequency

def plot_codons_grouped(aa_codon_freq):
  # Prepare data for plotting
  amino_acids = sorted(aa_codon_freq.keys())
  grouped_data = {codon: [aa_codon_freq[aa].get(codon, 0) for aa in amino_acids] for codon in code.values()}

  # Number of bars for each amino acid
  n_bars = len(amino_acids)

  # Positions of the bars on the x-axis
  ind = np.arange(n_bars)

  # Width of each bar
  bar_width = 0.35

  plt.figure(figsize=(15, 8))

  # Create bars for each codon
  for i, (codon, freqs) in enumerate(grouped_data.items()):
      plt.bar(ind + i * bar_width, freqs, width=bar_width, label=codon)

  # Add labels, title, and legend
  plt.xlabel('Amino Acid')
  plt.ylabel('Frequency')
  plt.title('Codon Usage Grouped by Amino Acid')
  plt.xticks(ind + bar_width / 2, amino_acids)
  plt.legend(title="Codon")

  plt.show()

def plot_codons_grouped(aa_codon_freq):
  amino_acids = sorted(aa_codon_freq.keys())
  
  # Number of bars for each amino acid
  n_bars = len(amino_acids)
  
  # Flatten the data for plotting
  # Flatten the data for plotting
  flat_data = []
  for aa, codons in aa_codon_freq.items():
      for codon, freq in codons.items():
          flat_data.append((f"{codon} ({aa})", freq))

  # Sort data for consistent ordering
  flat_data.sort(key=lambda x: codon_sort_key(x[0]))

  # Separate labels and frequencies
  labels, freqs = zip(*flat_data)

  # Create a color map for amino acids
  amino_acids = sorted(set(aa for aa, _ in aa_codon_freq.items()))
  colors = plt.cm.tab20(np.linspace(0, 1, len(amino_acids)))
  aa_color_map = {aa: colors[i] for i, aa in enumerate(amino_acids)}

  # Assign colors to each bar based on amino acid
  bar_colors = [aa_color_map[label.split('(')[-1].strip(')')] for label in labels]

  plt.figure(figsize=(15, 6))
  plt.bar(labels, freqs, color=bar_colors)
  plt.xlabel('Codon (Amino Acid)')
  plt.ylabel('Frequency')
  plt.title('Codon Usage Frequency by Amino Acid')
  plt.xticks(rotation=90)  # Rotate labels for readability
  plt.show()


def main():
  parser = argparse.ArgumentParser(
    prog='fasta_stats.py',
    description='print stats for a fasta file'
  )
  parser.add_argument('-i', '--input', default='/dev/fd/0', help="input file, default stdin")
  args = parser.parse_args()

  print("reading sequences")
  sequences = read_fasta(args.input)
  print("evaluating stats")
  sequence_lengths, gc_contents = print_fasta_statistics(sequences)
  #plot_box_plots_with_stats(sequence_lengths, gc_contents)
  plot_correlation(sequence_lengths, gc_contents)
  #codon_freq = calculate_codon_frequency(sequences.values())
  #print(codon_freq)
  #plot_codon_frequency(codon_freq)
  #groups = group_codon_frequency_by_amino_acid(sequences.values())
  #print(groups)
  #plot_codons_grouped(groups)
if __name__ == "__main__":
  main()