import os
import os.path
import sys

def codons(seq): # str -> [(int, str)]
    for i in range(0, len(seq), 3):
        if i+3 <= len(seq):
            yield i, seq[i:i+3]

# generate reading frames, return index of start codon and index of end codon
def reading_frames(seq): # [(int, str)] -> iter (int, int)
    for i in range(len(seq)):
        if seq[i][1] == "ATG":
            for j in range(i, len(seq)):
                if seq[j][1] in ["TAG", "TGA", "TAA"]:
                    yield i, j
                    break

def print_frames(sequences):
    for c, seq in enumerate(sequences):
        for i,j in reading_frames(seq):
            if abs(i - j) < 30:
                continue
            print(f"seq {c}: {seq[i][0]: >5} {seq[j][0]}: ", end="")
            for k in range(i,j+1):
                print(f"{seq[k][1]} ", end="")
            print()


translation_table = str.maketrans('ATGC', 'TACG')
def complement_dna(seq): # str -> str
    return seq.translate(translation_table)

def complementary_sequence(seq): # [(int, str)] -> [(int, str)]
    return [ (i + 2, complement_dna(c[::-1])) for i,c in seq[::-1] ]

def main():
    # declare vars to be available interactively
    global sequence
    global sequences

    with open("Xanthomonas.fna") as f:
        lines = f.readlines()
        lines = lines[1:] # cut out header

        sequence = "".join(line[:-1] for line in lines) # remove newlines
        sequences = [
            sequence, 
            sequence[1:], 
            sequence[2:]
        ]
        sequences = [ list(codons(s)) for s in sequences ]
        sequences = sequences + [ complementary_sequence(s) for s in sequences ]

        print_frames(sequences)

if __name__ == "__main__":
    main()