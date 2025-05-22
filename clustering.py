# vim: set foldmethod=marker:
# vim: set foldclose=all foldlevel=0:
# vim: set foldenable: 

# TODO: add docstrings

from typing import List, Set, Optional
import subprocess
import re
import argparse
import os

import fasta as fs
import io_helpers as io

def flatten_fastas( # {{{
        fastas: List[fs.Fasta],
        names: List[str] = [],
        sep: str = "_",
        names_sep: str = "-"
) -> str:
    """
    Creates a single string of a list of fastas with optional identifiers
    """
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
    """
    Turns a list of fastas into a single fasta
    """
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

def read_clusters( # {{{
    filepath: str,
) -> List[Set[str]]:
    file_content = io.read_file(filepath, lines=True)
    clusters = []
    for line in file_content:
        cluster = line.strip()[1:-1]
        proteins = [protein[1:-1] for protein in cluster.split(", ")]
        clusters.append(set(proteins))
    return clusters
# }}}

def write_clusters( # {{{
    filepath: str,
    clusters: List[Set[str]]
) -> None:
    file_string = "\n".join([str(cluster) for cluster in clusters])
    io.write_file(filepath, file_string)
# }}}

def translate_proteins( # {{{
    clusters: List[Set[str]],
    database: fs.Fasta,
    sep: str = "-"
) -> List[Set[str]]:
    clusters_new = []
    for cluster in clusters:
        cluster_new = set()
        for protein in cluster:
            bin, _, id = protein.partition(sep)
            for sequence in database.search(id):
                cluster_new.add(f"{bin}-{sequence.header}")
        clusters_new.append(cluster_new)
    return clusters_new
# }}}

def main( # {{{
        data_file: str,
        diamond: str,
        threshold: int,
        output: Optional[str] = None,
        verbose: bool = False
):
    # Get data paths and identifiers from csv file
    data = io.parse_csv(
        data_file,
        columns = 2,
        sep = ",",
        pattern = r"([a-zA-Z0-9/_\.]*)"
    )

    # Create big List of fasta including bin names
    fastas, names = [], []
    for file, identifier in data:
        fasta = fs.Fasta()
        fasta.read(file)
        fastas.append(fasta)
        names.append(identifier)
    fasta = concat_fastas(fastas, names=names)
    
    # save to temporary file and run diamond on it
    fasta.write("diamond_in")
    command = [
        diamond,
        "cluster",
        "-d",
        "diamond_in.fasta",
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

    # read diamond output and build cluster list from it
    connections = io.read_file("diamond_out.tsv", lines=True)
    clusters = build_clusters(connections)
    clusters = translate_proteins(clusters, fasta)

    # clean temporary files
    try:
        os.remove("diamond_out.tsv")
        os.remove("diamond_in.fasta")
    except Exception as e:
        print(f"There was an error removing temporary files: {e}")

    # write to file or return
    if output:
        write_clusters(output, clusters)
    else:
        return clusters
# }}}

if __name__ == "__main__": # {{{
    parser = argparse.ArgumentParser(prog="clustering.py") # {{{

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
    # }}}

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
# }}}
