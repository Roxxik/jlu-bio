#!/usr/bin/env python3

import sys
import argparse
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

def get_similarity(seq1, seq2):
  length = 0
  matches = 0
  for a,b in zip(seq1, seq2):
    length += 1
    if a == b != '-':
      matches += 1
  return matches / length if length != 0 else 0

# only call this function with exactly two sequences
def calculate_mutation_rates(sequences, window_length):
  seq1, seq2 = sequences.values()
  mutation_rates = []
  for i in range(0, len(seq1) - window_length + 1):
    window1 = seq1[i:i+window_length]
    window2 = seq2[i:i+window_length]
    similarity = get_similarity(window1, window2)
    mutation_rates.append(1 - similarity)
  return mutation_rates

def plot_mutation_rates(mutation_rates, window_length, file_name):
  plt.figure(figsize=(10, 4))
  plt.plot(mutation_rates, label=f'Window Length: {window_length}')
  plt.xlabel('Position')
  plt.ylabel('Mutation Rate')
  plt.title('Mutation Rates Across the Sequences')
  plt.legend()
  plt.savefig(file_name.replace('.fasta', f'_{window_length}.png'))
  plt.close()

def main():
  # parse args
  parser = argparse.ArgumentParser(
    prog='mutation',
    description='compare mutation of two aligned sequences, plot the result'
  )
  parser.add_argument('fasta_file', help="a fasta file with exactly two aligned sequences")
  parser.add_argument('-l', '--length', default=500, help="the window length to compute compute the mutation in")
  args = parser.parse_args()

  fasta_file = args.fasta_file
  window_length = int(args.length)

  # read file
  if not fasta_file.endswith('.fasta'):
    print("Please provide a fasta file with .fasta ending", file=sys.stderr)
    sys.exit(1)
  sequences = read_fasta(fasta_file)

  if len(sequences) != 2:
    print("There are not exactly two sequences in the fasta file", file=sys.stderr)
    sys.exit(1)

  # print whole genome stats
  print(f"whole sequence mutation: {(1 - get_similarity(*sequences.values()))*100:.3f}%")

  # calculate windows
  mutation_rates = calculate_mutation_rates(sequences, window_length)

  # plot windows
  plot_mutation_rates(mutation_rates, window_length, fasta_file)
  
  # gather stats on windows
  max_rate = max(mutation_rates)
  min_rate = min(mutation_rates)
  max_pos = mutation_rates.index(max_rate)
  min_pos = mutation_rates.index(min_rate)

  # print stats on windows
  print(f'max mutation rate: {max_rate*100:>7.3f}% at position: {max_pos:>7}')
  print(f'min mutation rate: {min_rate*100:>7.3f}% at position: {min_pos:>7}')

if __name__ == "__main__":
  main()