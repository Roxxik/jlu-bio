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