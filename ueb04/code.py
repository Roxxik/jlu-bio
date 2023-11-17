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
    "N": "TCAG",
    "R": "AG",
    "Y": "TC",
    "H": "TCA"
}

def shuffle(s):
    return ''.join([s[1], s[0], s[2]])

def codon_sort_key(codon):
    table = str.maketrans('TCAG', 'ABCD')
    return shuffle(codon.translate(table))

def expand_codons(code_rev, translation):
    code = dict()
    for aa, patterns in code_rev.items():
        for pattern in patterns:
            codons = [pattern]
            for char, replacements in translation.items():
                new_codons = []
                for codon in codons:
                    if char in codon:
                        for replacement in replacements:
                            new_codons.append(codon.replace(char, replacement, 1))
                    else:
                        new_codons.append(codon)
                codons = new_codons

            for codon in codons:
                code[codon] = aa
    result = dict(sorted(code.items(), key=lambda x: codon_sort_key(x[0])))
    return result

code = expand_codons(code_rev, translation)