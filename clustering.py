# vim: set foldmethod=marker :
# TODO: add docstrings

from typing import List, Set, Optional
import subprocess
import re
import argparse
import os

import fasta as fs
import io_helpers as io

def merge_fastas( # {{{
        fastas: List[fs.Fasta],
        names: List[str] = [],
        sep: str = "_",
        names_sep: str = "-"
) -> str:
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
# }}}

def concat_fastas( # {{{
        fastas: List[fs.Fasta],
        names: List[str] = [], 
        names_sep: str = "-"
) -> fs.Fasta :
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
# }}}

def build_clusters( # {{{
        text: List[str]
) -> List[Set[str]]:
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
# }}}

def main( # {{{
        data_file: str,
        diamond: str,
        threshold: int,
        output: Optional[str] = None,
        verbose: bool = False
):
    data = io.parse_csv(
        data_file,
        columns = 2,
        sep = ",",
        pattern = r"([a-zA-Z0-9/_\.]*)"
    )
    fastas, names = [], []
    for file, identifier in data:
        fasta = fs.Fasta()
        fasta.read(file)
        fastas.append(fasta)
        names.append(identifier)
    fasta = concat_fastas(fastas, names=names)
    fasta.write("concat")
    command = [
        diamond,
        "cluster",
        "-d",
        "concat.fasta",
        "-o",
        "diamond_out.tsv",
        "--approx-id",
        str(threshold),
        "-M",
        "64G"
    ]
    subprocess.run(
        command,
        stdout=None if verbose else subprocess.DEVNULL,
        stderr=None if verbose else subprocess.DEVNULL
    )

    connections = io.read_file("diamond_out.tsv", lines=True)
    clusters = build_clusters(connections)

    try:
        os.remove("diamond_out.tsv")
        os.remove("concat.fasta")
    except Exception as e:
        print(f"There was an error removing temporary files: {e}")
    file_string = "\n".join([str(cluster) for cluster in clusters])
    if output:
        io.write_file(output, file_string)
    else:
        return clusters
# }}}

if __name__ == "__main__":
    parser = argparse.ArgumentParser(prog="clustering.py")

    parser.add_argument(
        "DATA",
        help="Path to a csv-file containing the data files and identifiers [file, name]",
        type=str
    )
    parser.add_argument(
        "-d",
        "--diamond",
        metavar = "PATH",
        help = "Path to the diamond executable, defaults to './diamond/diamond'",
        type = str
    )
    parser.add_argument(
        "-o",
        "--output",
        metavar = "FILE",
        help = "The output file to save the new clusters to. Omitting will return to stdout.",
        type = str
    )
    parser.add_argument(
        "-t",
        "--threshold",
        metavar = "THRESHOLD",
        help = "The minimal similarity between protein sequences to be clustered, default 90.",
        type = int
    )
    parser.add_argument(
        "-v",
        "--verbose",
        action = "store_true",
        help = "Set to show the output of diamond."
    )

    args = parser.parse_args()

    # Defaults
    THRESHOLD = 30 if not args.threshold else args.threshold
    DIAMOND = "./diamond/diamond" if not args.diamond else args.diamond-executable

    output = main(
        data_file = args.DATA,
        diamond = DIAMOND,
        threshold = THRESHOLD,
        output = args.output,
        verbose = args.verbose
    )
    if output:
        for cluster in output:
            print(cluster)
