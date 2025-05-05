import subprocess

import io_helpers

data_folder = "../data"

paths = []
for bin in range(1,10):
    paths.append(f"{data_folder}/bin.{bin}/bin.{bin}.fna")

print("Paths: ")
for path in paths:
    print(path)
print("")

sequences = []
for index, path in enumerate(paths):
    sequence = io_helpers.get_sequence(path)
    if sequence:
        element = (f"sequence_{index}", sequence)
        sequences.append(element)

if False:
    for index, element in enumerate(sequences):
        print(f"Element {index}: ({element[0]}, {element[1][1:20]}...)")
        print(f"Sequences: {type(sequences)}")
        print(f"Element: {type(element)}")
        print(f"Header: {type(element[0])}")
        print(f"Header: {type(element[1])}")

io_helpers.write_fasta(sequences, "test.fasta")
