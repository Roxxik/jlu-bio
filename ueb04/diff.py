# compare getorf output with prodigal output
# write file with prodigal id and getorf id, write - if not found
# a gene is identified by start, stop, and read direction
# getorf writes [start - stop] regardless of read direction, if start > stop the gene is on the complementary strand

# call getorf as $ getorf -sequence file.fna -outseq file.getorf -table 11 -find 3
# call prodigal as $ prodigal -i file.fna -f gff -o file.gff -d nuc_file
import argparse
from collections import defaultdict

def main():
    global result
    parser = argparse.ArgumentParser(
        prog='diff',
        description='compare reading frames'
    )
    parser.add_argument('getorf_file', help="output of getorf")
    parser.add_argument('gff_file', help="output of prodigal")
    args = parser.parse_args()

    # getorf header format: 
    # >ID [start - end] other info
    getorf = dict()
    with open(args.getorf_file) as f:
        for i, line in enumerate(f.readlines()):
            if line[0] != ">": # only look at header lines
                continue
            id, begin, _, end, *_ = line.split(" ")
            id = id[1:]
            begin = int(begin[1:])
            end = int(end[:-1])
            dir = 1
            if begin > end:
                end, begin = begin, end
                dir = -1
                begin -= 3 # getorf does not output the stop codon, prodigal does
            else:
                end += 3
            start, stop = begin, end
            getorf[(stop, dir)] = (id, start, stop)
            # print(f"{id}\t{begin}\t{end}\t{' +-'[dir]}")
            """
            print((id, start, stop, dir))
            if i > 10:
                return
            """

    prodigal = dict()
    with open(args.gff_file) as f:
        for i, line in enumerate(f.readlines()):
            if line[0] == '#': # skip comment lines
                continue
            genome_id, _, _, begin, end, _, dir, _, notes = line.split("\t")
            gene_id, *_ = notes.split(";")


            id = genome_id + "/" + gene_id[3:]
            begin = int(begin)
            end = int(end)
            if dir == "+":
                start, stop = begin, end
            else:
                start, stop = end, begin
            dir += "1"
            dir = int(dir)

            prodigal[(stop, dir)] = (id, start, stop)
            print(f"{id}\t{begin}\t{end}\t{' +-'[dir]}")
        return

    merged = defaultdict(dict)
    for k,v in getorf.items():
        merged[k]["getorf"] = v
    for k,v in prodigal.items():
        merged[k]["prodigal"] = v

    result = list()
    for i, (k, v) in enumerate(merged.items()):
        _, dir = k
        getorf_id, start_g, stop_g = v.get("getorf", ("@", -1, -1))
        prodigal_id, start_p, stop_p = v.get("prodigal", ("@", -1, -1))

        print(f"{getorf_id}\t{prodigal_id}\t{' +-'[dir]}\t{start_g}\t{stop_g}\t{start_p}\t{stop_p}")
        #print(f"{getorf_id}\t{prodigal_id}\t{begin}\t{end}\t{' +-'[dir]}\t{start_g}\t{stop_g}\t{start_p}\t{stop_p}")
        #print(f"{getorf_id}\t{prodigal_id}\t{start_g if start_g != -1 else start_p}\t{stop}\t{' +-'[dir]}")
        #result.append((getorf_id, prodigal_id, start_g, start_p, stop, dir))

if __name__ == "__main__":
    main()