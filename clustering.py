# vim: set foldmethod=marker:
# vim: set foldclose=all foldlevel=0:
# vim: set foldenable: 

# TODO: add docstrings

from typing import List, Set, Optional
import subprocess
import re
import argparse
import os
import time

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
    if names != []:
        for index, fasta in enumerate(fastas):
            fasta.redit(edit=(r"^>(.+)$", rf">{names[index]}{names_sep}\1"))
    [result.add(fasta) for fasta in fastas]
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

def main( # {{{
        data_file: str,
        diamond_path: str,
        threshold: int,
        verbose: bool = False
):
    data = io.parse_csv(
        data_file,
        sep = ",",
    )
    # Create big List of fasta including bin names
    fastas = [(fasta := fs.Fasta()).read(file) or fasta for file, _ in data]
    names = [name for _, name in data]
    fasta = concat_fastas(fastas, names=names)

    start_time = time.time()
    clusters = diamond(fasta, threshold, verbose=verbose, executable=diamond_path)
    end_time = time.time()
    if verbose: print(f"Diamond took: {(end_time-start_time):.4f}s")

    start_time = time.time()
    clusters = grow_clusters(clusters)
    end_time = time.time()
    if verbose: print(f"Growing took: {(end_time-start_time):.4f}s")

    start_time = time.time()
    clusters = purge_clusters(clusters)
    end_time = time.time()
    # [print(cluster) for cluster in clusters]
    if verbose: print(f"Purging took: {(end_time-start_time):.4f}s")

    start_time = time.time()
    clusters = parse_clusters(clusters, fasta)
    end_time = time.time()
    if verbose: print(f"Parsing took: {(end_time-start_time):.4f}s")

    return clusters
# }}}

def grow_clusters( #{{{
    clusters: List[Set[str]]
) -> List[Set[str]]:
    """
    Merges all clusters with overlapping elements until only disjoint Sets (Clusters) are left
    """
    while True:
        new_clusters = []
        merged = False
        while clusters:
            seed = clusters.pop(0)
            for cluster in clusters:
                if not seed.isdisjoint(cluster):
                    seed = seed.union(cluster)
                    clusters.remove(cluster)
                    merged = True
            new_clusters.append(seed)
        clusters = new_clusters
        if not merged:
            break
    return clusters
#}}}

def parse_clusters( #{{{
    clusters: List[Set[str]],
    reference: fs.Fasta
) -> List[fs.Fasta]:
    """
    Turn a list of clusters into a list of their respective fasta objects.
    Raises a ValueError if an object is ambiguous
    """
    result = []
    for cluster in clusters:
        new_cluster = fs.Fasta()
        for sequence in cluster:
            search_result = reference.search(sequence)
            if len(search_result) != 1:
                raise ValueError(f"The search in the reference Fasta was faulty (found: {len(search_result)})")
            else:
                new_cluster.add(search_result[0])
        result.append(new_cluster)
    return result
#}}}

def diamond( #{{{
    fasta: fs.Fasta,
    threshold: int,
    executable: str = "./diamond/diamond",
    verbose:bool = False
) -> List[fs.Fasta]:
    fasta.write("diamond_in")
    command = [
        executable,
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
    
    diamond_out = io.read_file("diamond_out.tsv", lines=True)

    # build the clusters
    clusters = []
    for line in diamond_out:
        protein1, _, protein2 = line.partition("\t")
        clusters.append(set([protein1, protein2]))
    # clusters = grow_clusters(clusters)

    # clean temporary files
    try:
        os.remove("diamond_out.tsv")
        os.remove("diamond_in.fasta")
    except Exception as e:
        print(f"There was an error removing temporary files: {e}")

    return clusters
#}}}

def purge_clusters( #{{{
    clusters:List[Set[str]],
    min:int = 2
) -> List[Set[str]]: 
    """
    Remove all clusters with less than a specified amount of sequences
    """
    return [cluster for cluster in clusters if len(cluster) >= min]
#}}}

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
        help = "Set to show the output of diamond and additional timing information."
    )

    args = parser.parse_args()
    # }}}

    # Defaults
    THRESHOLD = 30 if not args.threshold else args.threshold
    DIAMOND = "./diamond/diamond" if not args.diamond else args.diamond-executable

    clusters = main(
        data_file = args.DATA,
        diamond_path = DIAMOND,
        threshold = THRESHOLD,
        verbose = args.verbose
    )
    [print(cluster) for cluster in clusters]
# }}}
