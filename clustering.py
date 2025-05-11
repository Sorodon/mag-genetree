# TODO: add docstrings

from typing import List, Set
import subprocess
import re

import fasta as fs
import io_helpers as io

def merge_fastas(fastas: List[fs.Fasta], names: List[str] = [], sep: str = "_", names_sep: str = "-") -> str:
    if len(names) != 0 and len(fastas) != len(names):
        raise ValueError("names and fastas do not match in length!")
    strings = []
    for index, fasta in enumerate(fastas):
        for seq in fasta:
            if len(names) != 0:
                strings.append(f"{names[index]}{names_sep}{seq.sequence}")
            else:
                strings.append(seq.sequence)
    return sep.join(strings)

def concat_fastas(fastas: List[fs.Fasta], names: List[str] = [], names_sep: str = "-") -> fs.Fasta :
    if len(names) != 0 and len(fastas) != len(names):
        raise ValueError("names and fastas do not match in length!")
    result = fs.Fasta()
    for index, fasta in enumerate(fastas):
        for seq in fasta:
            if len(names) != 0:
                new_header = f">{names[index]}{names_sep}{seq.header[1:]}"
                new_Sequence = fs.Sequence(header=new_header, sequence=seq.sequence)
            result.add(new_Sequence)
    return result

def build_clusters(text: List[str]) -> List[Set[str]]:
    result = []
    pattern = r"^(bin\.[0-9]{1,2}-\w{6}_\d{5})\t(bin\.[0-9]{1,2}-\w{6}_\d{5})"
    
    for line in text:
        match = re.search(pattern, line)
        if match:
            item1, item2 = match.group(1), match.group(2)
            if item1 != item2:
                connection = False
                for cluster in result:
                    if item1 in cluster or item2 in cluster:
                        connection = True
                        cluster.update({item1, item2})
                        break  # Avoid redundant checks
                if not connection:
                    result.append({item1, item2})
    return result

if __name__ == "__main__":
    THRESHOLD = 30
    DIAMOND_EXECUTABLE = "./diamond/diamond"
    DATA_FILES = [
        "../data/bin.1/bin.1.faa",
        "../data/bin.2/bin.2.faa",
        "../data/bin.3/bin.3.faa",
        "../data/bin.4/bin.4.faa",
        "../data/bin.5/bin.5.faa",
        "../data/bin.6/bin.6.faa",
        "../data/bin.7/bin.7.faa",
        "../data/bin.8/bin.8.faa",
        "../data/bin.9/bin.9.faa",
        "../data/bin.10/bin.10.faa",
        "../data/bin.11/bin.11.faa",
        "../data/bin.12/bin.12.faa",
        "../data/bin.13/bin.13.faa",
        "../data/bin.14/bin.14.faa",
        "../data/bin.15/bin.15.faa",
        "../data/bin.16/bin.16.faa",
        "../data/bin.17/bin.17.faa",
    ]
    IDENTIFIERS = [
        "bin.1",
        "bin.2",
        "bin.3",
        "bin.4",
        "bin.5",
        "bin.6",
        "bin.7",
        "bin.8",
        "bin.9",
        "bin.10",
        "bin.11",
        "bin.12",
        "bin.13",
        "bin.14",
        "bin.15",
        "bin.16",
        "bin.17",
    ]
    
    fastas = []
    for file in DATA_FILES:
        fasta = fs.Fasta()
        fasta.read(file)
        fastas.append(fasta)
    fasta = concat_fastas(fastas, names=IDENTIFIERS)
    fasta.write("temp")
    command = [DIAMOND_EXECUTABLE, "cluster", "-d", "temp.fasta", "-o", "clusters.tsv", "--approx-id", str(THRESHOLD), "-M", "64G"]
    subprocess.run(command)

    connections = io.read_lines("clusters.tsv")
    clusters = build_clusters(connections)

    file_string = "\n".join([str(cluster) for cluster in clusters])
    io.write_file("clusters.txt", file_string)
