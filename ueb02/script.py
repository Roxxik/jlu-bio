import argparse
from collections import defaultdict
import csv
import os
import os.path
import subprocess
import sys

e_epsilon = 1e-20

# expects a standard blastp output from diamond
def read_matches(match_file):
    matches = dict()
    with open(match_file) as f:
        reader = csv.reader(f, delimiter="\t")
        for row in reader:
            protein = row[0]
            match = row[1]
            e_value = float(row[10])
            score = row[11]
            if e_value > e_epsilon:
                continue
            if protein not in matches or score < matches[protein][2]:
                matches[protein] = (match, e_value, score)
    return matches

def main():
    parser = argparse.ArgumentParser(
        prog='script',
        description='compare proteomes'
    )
    parser.add_argument('target', help="One of the filenames in the proteomes folder")
    args = parser.parse_args()

    target = args.target
    target_file = os.path.basename(target)

    target_name, target_extension = os.path.splitext(target_file)

    if target_extension != '.faa':
        print("only works with .faa file endings", file=sys.stderr)
        return

    print("target: " + target_name)

    try:
        os.mkdir("out")
    except FileExistsError:
        pass

    # compare target to each reference
    matches_reciprocal = defaultdict(dict)
    references = []
    for file in os.listdir('proteomes'):
        reference = os.path.join('proteomes', file)
        if target == reference:
            continue
        references.append(reference)

        # query target in reference
        print("##########")
        print("query target in reference")
        print("##########")
        matches1 = os.path.join("out", "matches1.tsv")
        subprocess.call([
            "/vol/software/bin/diamond-2.0.0", 
            "blastp", 
            "-d", 
            reference,
            "-q",
            target,
            "-o",
            matches1
        ])

        # read forward matches
        matches = read_matches(matches1)

        # query reference in target
        matches2 = os.path.join("out", "matches2.tsv")
        subprocess.call([
            "/vol/software/bin/diamond-2.0.0", 
            "blastp", 
            "-d", 
            target,
            "-q",
            reference,
            "-o",
            matches2
        ])

        # read reverse matches
        matches_reverse = read_matches(matches2)

        # for each protein in target find the reciprocal best match in reference
        for protein in matches:
            match = matches[protein][0]
            try:
                reciprocal = matches_reverse[match][0]
            except KeyError:
                continue
            if protein == reciprocal:
                matches_reciprocal[protein][reference] = match
    
    # write output
    with open(os.path.join("out", "reciprocal.tsv"), "w") as f:
        writer = csv.writer(f, delimiter="\t")
        for protein in matches:
            row = [protein]
            for reference in references:
                if reference in matches_reciprocal[protein]:
                    row.append(matches_reciprocal[protein][reference])
                else:
                    row.append("no reciprocal")
            writer.writerow(row)


if __name__ == "__main__":
    main()