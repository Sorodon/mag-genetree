# vim: set foldmethod=marker:
# vim: set foldclose=all foldlevel=0:
# vim: set foldenable: 

# TODO: add docstrings

from typing import List, Set, Tuple, Union
import subprocess, os

import bakta_table as bt
import clustering as cl
import fasta as fs
import io_helpers as io


def verify_clusters( #{{{
    clusters: List[Set[str]],
    lookup_table: bt.Bakta_table,
) -> Union[List[Tuple[Set[str], Set[str]]], dict]:
    results = []
    stats = {"total": len(clusters)}
    for cluster in clusters:
        seen = set()
        for protein in cluster:
            bin, _, locus_tag = protein.partition("-")
            lookup = lookup_table.find(locus_tag, "Locus Tag")
            for match in lookup:
                id = bt.get_uniprot(match, level=50)
                if id is not None:
                    seen.add(id)
        results.append((cluster, seen))
#        stats[count := str(len(seen))] = stats.get(count, 0) + 1
        stats[str(len(seen))] = stats.get(str(len(seen)), 0) + 1
    return results, stats
#}}}

if __name__ == "__main__":
    data = io.parse_csv("data.csv", sep=",")
    fastas, names = [], []
    for file, identifier in data:
        fasta = fs.Fasta()
        fasta.read(file)
        fastas.append(fasta)
        names.append(identifier)
    fasta = cl.concat_fastas(fastas, names=names)
    fasta.write("diamond_in")
    command = [
        "./diamond/diamond",
        "cluster",
        "-d",
        "diamond_in.fasta",
        "-o",
        "diamond_out.tsv",
        "--approx-id",
        "95",
        "-M",
        "64G"
    ]
    verbose = False
    subprocess.run(
        command,
        stdout=None if verbose else subprocess.DEVNULL,
        stderr=None if verbose else subprocess.DEVNULL
    )

    # read diamond output and build cluster list from it
    connections = io.read_file("diamond_out.tsv", lines=True)
    clusters = cl.build_clusters(connections)

    # clean temporary files
    try:
        os.remove("diamond_out.tsv")
        os.remove("diamond_in.fasta")
    except Exception as e:
        print(f"There was an error removing temporary files: {e}")


    print("Diamond is done!\n")

    print("Creating lookup table\n")
    bakta_table = bt.Bakta_table()
    bakta_table.read("../data/bin.1/bin.1.tsv", skip=5)
    bakta_table.read("../data/bin.2/bin.2.tsv", skip=5)
    bakta_table.read("../data/bin.3/bin.3.tsv", skip=5)
    bakta_table.read("../data/bin.4/bin.4.tsv", skip=5)
    bakta_table.read("../data/bin.5/bin.5.tsv", skip=5)
    bakta_table.read("../data/bin.6/bin.6.tsv", skip=5)
    bakta_table.read("../data/bin.7/bin.7.tsv", skip=5)
    bakta_table.read("../data/bin.8/bin.8.tsv", skip=5)
    bakta_table.read("../data/bin.9/bin.9.tsv", skip=5)
    bakta_table.read("../data/bin.10/bin.10.tsv", skip=5)
    bakta_table.read("../data/bin.11/bin.11.tsv", skip=5)
    bakta_table.read("../data/bin.12/bin.12.tsv", skip=5)
    bakta_table.read("../data/bin.13/bin.13.tsv", skip=5)
    bakta_table.read("../data/bin.14/bin.14.tsv", skip=5)
    bakta_table.read("../data/bin.15/bin.15.tsv", skip=5)
    bakta_table.read("../data/bin.16/bin.16.tsv", skip=5)
    bakta_table.read("../data/bin.17/bin.17.tsv", skip=5)

    print(f"Number of clusters: {len(clusters)}")
    filtered, stats = verify_clusters(clusters, bakta_table)
    [print(f"{key}: {value} ({value/len(clusters)*100:.2f}%)") for key,value in stats.items() if key != "total"]
