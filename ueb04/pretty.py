import argparse
from code import code

# generator that consumes input lines of nucleotides
# and prints them to a file in a table. with indices and amino acids
# no validation is done in here. every charater is used as a nucleotide
def tabulator(output_file, begin=None, dir=None):
    if begin is None:
        begin = 0
    if dir is None:
        dir = 1
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
                print(f"{begin + dir*i:0>10}: ", end="", file=output_file)
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

def mktabulator(*args, **kwargs):
    generator = tabulator(*args, **kwargs)
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
                # try to read prodigal header
                parts = line.split("#")
                if len(parts) == 5:
                    id, begin, end, dir, notes = parts
                    begin = int(begin.strip())
                    end = int(end.strip())
                    dir = int(dir.strip())
                    if dir == 1:
                        gen = mktabulator(output_file, begin, dir)
                    else:
                        gen = mktabulator(output_file, end, dir)
                else:
                    gen = mktabulator(output_file)
            else:
                gen.send(line)
        stoptabulator(gen)

if __name__ == "__main__":
    main()