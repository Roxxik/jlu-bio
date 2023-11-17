import argparse
from collections import defaultdict
import csv
import os
import os.path
import sys

e_epsilon = 1e-20

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
    parser.add_argument('-o', '--output', help="Output filename", required=True)
    args = parser.parse_args()

    output_file = args.output

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

        reference_name, reference_extension = os.path.splitext(file)
        if reference_extension != '.faa':
            continue

        references.append(reference)

        matches_forward_file = os.path.join("out", f"{target_name}_v_{reference_name}.tsv")
        matches = read_matches(matches_forward_file)

        matches_reverse_file = os.path.join("out", f"{reference_name}_v_{target_name}.tsv")
        matches_reverse = read_matches(matches_reverse_file)

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
    with open(output_file, "w") as f:
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