# vim: set foldmethod=marker:
# vim: set foldclose=all foldlevel=0:
# vim: set foldenable: 

from typing import List, Set, Optional
import subprocess
import re
import argparse
import os
import time

import fasta as fs
import io_helpers as io

def concat_fastas( # {{{
        fastas: List[fs.Fasta],
        names: List[str] = [], 
        names_sep: str = "-"
) -> fs.Fasta :
    """
    Turns a list of fastas into a single fasta
    Args:
        fastas (List[Fasta]): A list of fasta objects to be merged.
        names: (List[str]): A list of names to prepend to fasta headers to identify them later. Omiting will disable this feature.
        names_sep (str): A string to separate name from header when marking fastas with names. Defaults to '-'.
    Returns:
        Fasta: The merged Fasta object.
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

def main( # {{{
    data_file: str,
    executable: str,
    threshold: int,
    method:str = "cluster",
    verbose: bool = False,
    nopurge:bool = False
):
    """
    Main entrypoint into the clustering module
    Args:
        data_file (str): The data file to be processed (should contain comma separated fasta paths and names (e.g. 'some/fasta.faa,id1')).
        executable (str): The path to the diamond executable.
        threshold (int): The identity threshold to apply when clustering in percent. Should range from 100 to 0.
        method (str): Can be either 'cluster' or 'linclust' depending on the preferred clustering method. Defaults to 'cluster'.
        verbose (bool): Whether to print additional info like runtimes of different steps. Defaults to False.
        nopurge (bool): Whether to keep singluar clusters before parsing. Defaults to False.
    Returns:
        List[Fasta]: A list of Fasta objects, each one being one cluster.
    """
    data = io.parse_csv(
        data_file,
        sep = ",",
    )
    # Create big List of fasta including bin names
    fastas = [(fasta := fs.Fasta()).read(file) or fasta for file, _ in data]
    names = [name for _, name in data]
    fasta = concat_fastas(fastas, names=names)

    # Run diamond on it
    start_time = time.time()
    clusters = diamond(fasta, threshold, executable=executable)
    if verbose: print(f"Diamond took: {(time.time()-start_time):.4f}s")

    # Turn links into an actual list
    start_time = time.time()
    clusters = grow_clusters(clusters)
    end_time = time.time()
    if verbose: print(f"Growing took: {(end_time-start_time):.4f}s")

    # purge clusters with less than two members
    if not nopurge:
        start_time = time.time()
        clusters = purge_clusters(clusters)
        end_time = time.time()
        if verbose: print(f"Purging took: {(end_time-start_time):.4f}s")

    # Turn locus tags into actual sequences
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
    Args:
         clusters (List[Set[str]]): A list containing sets with identifier strings.
    Returns:
        List[Set[str]]: The resulting list of sets. All sets will be disjoint.
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
    Args:
        clusters (List[Set[str]]): A list of sets with sequence identifiers.
        reference (Fasta): A single Fasta object containing all sequences (with identifiers).
    Returns:
        List[Fasta]: A list of proper Fasta objects, each being one cluster.
    Raises:
        ValueError: If an identifier is ambiguous when looking it up in <reference>.
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
    method:str = "cluster",
    verbose:bool = False
) -> List[set]:
    """
    Run diamond on a fasta returning clusters
    Args:
        fasta (Fasta): The input Fasta object
        threshold (int): The identity treshold in percent.
        executable (str): Path of the diamond executable. Defaults to './diamond/diamond'.
        method (str): Either 'cluster' or 'linclust' depending on the preferred clustering method. Defaults to 'cluster'.
        verbose (bool): Whether to allow the output of diamond on stdout. Defaults to False.
    Returns:
        List[set]: A list of sets of cardinality 2, each containing two sequence identifiers.
    """
    fasta.write("diamond_in")
    command = [
        executable,
        method,
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
    Args:
        clusters (List[Set[str]]): A list of sets of sequence identifiers
        min (int): The minimal set size to keep. Defaults to 2.
    Returns:
        List[Set[str]]: The <clusters> input with all sets below the <min> threshold size removed.
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
    [print(cluster, "\n") for cluster in clusters]
# }}}
