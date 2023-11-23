def codons(seq): # str -> [(int, str)]
    for i in range(0, len(seq), 3):
        if i+3 <= len(seq):
            yield i, seq[i:i+3]

translation_table = str.maketrans('ATGC', 'TACG')
def complement_dna(seq): # str -> str
    return seq.translate(translation_table)

def complementary_sequence(seq): # [(int, str)] -> [(int, str)]
    return [ (i + 2, complement_dna(c[::-1])) for i,c in seq[::-1] ]


def frames(sequence):
    sequences = [
        sequence, 
        sequence[1:], 
        sequence[2:]
    ]
    sequences = [ list(codons(s)) for s in sequences ]
    sequences = sequences + [ complementary_sequence(s) for s in sequences ]
    return sequences