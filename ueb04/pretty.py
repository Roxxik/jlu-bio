import sys

def filtered(seq):
    for c in seq:
        if c not in " \n":
            yield c

seq = sys.stdin.read()
seq = filtered(seq)

# Table 11
code_rev = {
    "F": ["TTY"], # Phenylanalin
    "L": ["TTR", "CTN"], # Leucin
    "I": ["ATH"], # Isoleucin
    "M": ["ATG"], # Methionin
    "V": ["GTN"], # Valin
    "S": ["TCN", "AGY"], # Serin
    "P": ["CCN"], # Prolin
    "T": ["ACN"], # Threonin
    "A": ["GCN"], # Alanin
    "Y": ["TAY"], # Tyrosin
    "*": ["TAR", "TGA"], # STOP
    "H": ["CAY"], # Histidin
    "Q": ["CAR"], # Glutamin
    "N": ["AAY"], # Asparagin
    "K": ["AAR"], # Lysin
    "D": ["GAY"], # Asparaginsäure
    "E": ["GAR"], # Glutaminsäure
    "C": ["TGY"], # Cystein
    "W": ["TGG"], # Tryptophan
    "R": ["CGN", "AGR"], # Arginin
    "G": ["GGN"] # Glycin
}

translation = {
    "N": "ATGC",
    "R": "AG",
    "Y": "TC",
    "H": "ATC"
}

def expand_codons(code_rev, translation):
    # Initialize an empty dictionary to hold the expanded codon to amino acid mappings
    codon_to_aa = {}

    # Iterate over the amino acid codes and their corresponding codon patterns
    for aa, patterns in code_rev.items():
        # For each pattern, replace the special letters with the corresponding nucleotides
        for pattern in patterns:
            # Generate all possible combinations for the codon pattern
            expanded_patterns = [pattern]
            for char, replacements in translation.items():
                new_expanded_patterns = []
                # Replace each special character with each of its possible nucleotides
                for expanded_pattern in expanded_patterns:
                    if char in expanded_pattern:
                        for replacement in replacements:
                            new_expanded_patterns.append(expanded_pattern.replace(char, replacement, 1))
                    else:
                        new_expanded_patterns.append(expanded_pattern)
                expanded_patterns = new_expanded_patterns

            # Add the expanded patterns to the codon_to_aa dictionary
            for codon in expanded_patterns:
                codon_to_aa[codon] = aa

    return codon_to_aa

code = expand_codons(code_rev, translation)

n = 48

try:
    i = 0
    print(f"{i:>5}: ", end="")
    codons = list()
    while True:
        codon = ""
        
        # first base
        a = next(seq)
        codon += a
        i += 1
        print (a, end="")

        # second base
        b = next(seq)
        codon += b
        i += 1
        print (b, end="")

        # third base
        c = next(seq)
        codon += c
        i += 1
        print (c, end="")

        codons.append(codon)

        if i % n == 0:
            print("  " + "".join(code.get(codon, 'X') for codon in codons), end="")
            codons = list()
            print()
            print(f"{i:>5}: ", end="")
        else:
            print(" ", end="")

except StopIteration:
    pass
print((" " * (int ((n - i % n) * (4/3.0) - 0.5) + 2)) + "".join(code.get(codon, 'X') for codon in codons), end="")
print()