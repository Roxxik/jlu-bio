import sys
from collections import defaultdict
import argparse

def reverse_dict(d):
    res = defaultdict(list)
    for key, val in d.items():
        res[val].append(key)
    return res


def filtered(seq):
    for c in seq:
        if c not in " \n":
            yield c


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

# generator that consumes input lines of nucleotides
# and prints them to a file in a table. with indices and amino acids
# no validation is done in here. every charater is used as a nucleotide
def tabulator(output_file):
    try:
        line = yield
    except StopIteration:
        return
    seq = iter(line)
    try:
        i = 0
        codons = list()
        while True:
            codon = ""

            # first base
            a = None
            while a is None:
                try:
                    a = next(seq)
                except StopIteration:
                    try:
                        line = yield
                    except StopIteration:
                        return
                    seq = iter(line)
            # if we actually got a first base, print in between stuff
            # either the rows header (last rows amino acids and this rows offset)
            # or a space between tripletts
            if i % n == 0:
                if i != 0:
                    print("  " + "".join(code.get(codon, 'X') for codon in codons), end="", file=output_file)
                    print(file=output_file)
                    codons = list()
                print(f"{i:>5}: ", end="", file=output_file)
            else:
                print(" ", end="", file=output_file)

            codon += a
            i += 1
            print(a, end="", file=output_file)

            # second base
            b = None
            while b is None:
                try:
                    b = next(seq)
                except StopIteration:
                    try:
                        line = yield
                    except StopIteration:
                        return
                    seq = iter(line)
            codon += b
            i += 1
            print(b, end="", file=output_file)

            # third base
            c = None
            while c is None:
                try:
                    c = next(seq)
                except StopIteration:
                    try:
                        line = yield
                    except StopIteration:
                        return
                    seq = iter(line)
            codon += c
            i += 1
            print(c, end="", file=output_file)

            codons.append(codon)
    finally:
        # only do anything at all if we are in the midst of a line
        if i % n != 0:
            j = i
            # finish the last triplett
            k = (2 - ((j-1) % 3))
            print(" " * k, end="", file=output_file)
            j += k
            while j % n != 0:
                print("    ", end="", file=output_file)
                j += 3
        if i != 0 and i % n not in (1,2):
            print("  " + "".join(code.get(codon, 'X') for codon in codons), end="", file=output_file)
        if i != 0:
            print(file=output_file)

def mktabulator(output_file):
    generator = tabulator(output_file)
    next(generator) # start it up
    return generator

def stoptabulator(gen):
    try:
        gen.throw(StopIteration)
    except StopIteration:
        pass

def main():
    parser = argparse.ArgumentParser(
        prog='pretty.py',
        description='prettify nucleotid sequences'
    )
    parser.add_argument('-i', '--input', default='/dev/fd/0', help="input file, default stdin")
    parser.add_argument('-o', '--output', default='/dev/fd/1', help="output file, default stdout")
    parser.add_argument('-n', '--blocks', default=10, help="the number of tripletts in a row")
    args = parser.parse_args()
    global n
    n = int(args.blocks) * 3


    with open(args.input) as input_file, open(args.output, "w") as output_file:
        gen = mktabulator(output_file)
        for line in input_file:
            line = line.strip()
            if line.startswith(">"):
                stoptabulator(gen)
                print(line, file=output_file)
                gen = mktabulator(output_file)
            else:
                gen.send(line)
        stoptabulator(gen)

if __name__ == "__main__":
    main()