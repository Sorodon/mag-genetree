# vim: set foldmethod=marker:
# vim: set foldclose=all foldlevel=0:
# vim: set foldenable: 

# TODO: add docstrings

import argparse, os

import clustering as cl
import filtering as fl
import bakta_table as bt
import io_helpers as io

def main(
    data_file,
    diamond_path,
    threshold,
    verbose:bool = False,
    alignment_path:str = None,
    method:str = "cluster",
    size_threshold = 2,
    gaps_threshold = 1,
    uniref_lookup = None,
    uniref50_threshold = float('inf'),
    uniref90_threshold = float('inf'),
    uniref100_threshold = float('inf'),
    out_file = None
):
    
    # Step 1: Creating clusters with diamond
    clusters = cl.main(
        data_file = data_file,
        diamond_path = diamond_path,
        threshold = threshold,
        method = method,
        verbose = verbose
    )
    if verbose: print(f">>> Diamond found {len(clusters)} clusters")

    if verbose: print(">>> Start Aligning Clusters")
    # Step 2: Aligning
    clusters = [cluster.align() for cluster in clusters]
    
    # Step 3: Filtering
    # 3.1 By size
    if verbose: print(">>> Start filtering by size")
    clusters = fl.filter_size(clusters, size_threshold)
    # 3.2 By gaps
    if verbose: print(">>> Start filtering by gaps")
    clusters = fl.filter_gaps(clusters, gaps_threshold)
    # 3.3 By UniRef ID
    if uniref_lookup:
        paths = [row[0] for row in io.parse_csv(uniref_lookup)]
        lookup = bt.Bakta_table()
        lookup.read(paths, skip=5)
        # 3.3.1 UniRef100
        if verbose: print(">>> Start filtering by UniRef100 IDs")
        clusters = fl.filter_uniref(clusters, lookup, uniref100_threshold, level=100)
        # 3.3.2 UniRef90
        if verbose: print(">>> Start filtering by UniRef90 IDs")
        clusters = fl.filter_uniref(clusters, lookup, uniref90_threshold, level=90)
        # 3.3.3 UniRef50
        if verbose: print(">>> Start filtering by UniRef50 IDs")
        clusters = fl.filter_uniref(clusters, lookup, uniref50_threshold, level=50)

    if verbose: print(f">>> Filtering done ({len(clusters)} Clusters left)")
    
    # Saving Alginments
    if alignment_path:
        for index, cluster in enumerate(clusters):
            cluster.write(os.path.join(alignment_path, f"{index:0{len(str(len(clusters)))}}"))
    if verbose: print(f">>> Alignments have been saved to {alignment_path}")


    # Step 4: Distance matrices
    if verbose: print(">>> Start calculating distance matrices")
    for cluster in clusters:
        cluster.cd()

    # Step 5: Trees
    if verbose: print(">>> Start calculating distance matrices")
    trees = [cluster.distmat.upgma() for cluster in clusters]

    # Step 6: Output
    if out_file:
        if verbose: print(">>> Now writing to file")
        io.write_file(out_file, "\n".join(trees))
    else:
        [print(tree) for tree in trees]

if __name__ == "__main__":
    parser = argparse.ArgumentParser(prog="main.py") # {{{

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
        "--out",
        metavar = "FILE",
        help = "The output file to save the trees to. Omitting will return to stdout.",
        type = str
    )
    parser.add_argument(
        "-a",
        "--alignment_out",
        metavar = "FOLDER",
        help = "The path where to save the alignment for each cluster as separate fasta. Not saved if omitted.",
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
        "--linclust",
        help = "Use linclust instead of cluster mode for DIAMOND",
        action = "store_true"
    )
    parser.add_argument(
        "-v",
        "--verbose",
        action = "store_true",
        help = "Set to show more detailed output"
    )
    parser.add_argument(
        "--size",
        metavar = "CLUSTER_SIZE",
        help = "The minimal cluster size to keep, defaults to 2",
        type = int
    )
    parser.add_argument(
        "--gaps",
        metavar = "GAP_RATE",
        help = "The max gap rate to keep",
        type = float
    )
    parser.add_argument(
        "--lookup",
        metavar = "FILE",
        help = "A csv file containing paths to all necessary bakta files in column one",
        type = str
    )
    parser.add_argument(
        "--ur50",
        metavar = "THRESHOLD",
        help = "The max amount of distinct UniRef50 IDs per cluster",
        type = int
    )
    parser.add_argument(
        "--ur90",
        metavar = "THRESHOLD",
        help = "The max amount of distinct UniRef90 IDs per cluster",
        type = int
    )
    parser.add_argument(
        "--ur100",
        metavar = "THRESHOLD",
        help = "The max amount of distinct UniRef100 IDs per cluster",
        type = int
    )

    args = parser.parse_args()
    # }}}

    # Setting defaults
    THRESHOLD = 90 if not args.threshold else args.threshold
    DIAMOND = "./diamond/diamond" if not args.diamond else args.diamond

    params = {"data_file": args.DATA}
    params["diamond_path"] = DIAMOND
    params["threshold"] = THRESHOLD
    if args.alignment_out: params["alignment_path"] = args.alignment_out
    if args.linclust: params["method"] = "linclust"
    if args.verbose: params["verbose"] = args.verbose
    if args.size: params["size_threshold"] = args.size
    if args.gaps: params["gaps_threshold"] = args.gaps
    if args.lookup: params["uniref_lookup"] = args.lookup
    if args.ur50: params["uniref50_threshold"] = args.ur50
    if args.ur90: params["uniref90_threshold"] = args.ur90
    if args.ur100: params["uniref100_threshold"] = args.ur100
    if args.out: params["out_file"] = args.out

    main(**params)
