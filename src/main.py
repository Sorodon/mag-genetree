# vim: set foldmethod=marker:
# vim: set foldclose=all foldlevel=0:
# vim: set foldenable: 

# TODO: add docstrings

import argparse, os, time

import clustering as cl
import filtering as fl
import bakta_table as bt
import io_helpers as io

def main(
    data_file,
    executable,
    threshold,
    verbose:bool = False,
    timing:bool = False,
    alignment_path:str = None,
    method:str = "cluster",
    size_threshold = 3,
    gaps_threshold = 1,
    length_threshold = 1,
    uniref_lookup = None,
    uniref50_threshold = float('inf'),
    uniref90_threshold = float('inf'),
    uniref100_threshold = float('inf'),
    out_file = None
):
    start_time_total = time.time()   
    # Step 1: Creating clusters with diamond
    if verbose: print(">>> Start Clustering")
    time_clustering = time.time()
    clusters = cl.main(
        data_file = data_file,
        executable = executable,
        threshold = threshold,
        method = method,
        verbose = timing
    )
    if verbose: print(f"Clustering found {len(clusters)} clusters")
    if timing: print(f"Clustering took: {time.time()-time_clustering:.4f}s")

    # Step 2: Filter by size
    if verbose: print(f">>> Start filtering by size ({len(clusters)} clusters left)")
    time_filter_size = time.time()
    clusters = fl.filter_size(clusters, size_threshold)
    if timing: print(f"Size filtering took: {time.time()-time_filter_size:.4f}s")

    # Step 3: Aligning matrices
    # if verbose: print(f">>> Start Aligning Clusters")
    # time_aligning = time.time()
    # clusters = [cluster.align() for cluster in clusters]
    # if timing: print(f"Aligning took: {time.time()-time_aligning:.4f}s")

    # Step 4: Distance matrices
    # if verbose: print(">>> Start calculating distance matrices")
    # time_distmatrices = time.time()
    # for cluster in clusters:
    #     cluster.cd()
    # if timing: print(f"Calculating distance matrices took: {time.time()-time_distmatrices:.4f}s")

    # Step 3: Clustalo
    if verbose: print(">>> Start Alignment and Distance matrix calculation")
    time_clustalo = time.time()
    clusters = [cluster.clustalo() for cluster in clusters]
    if timing: print(f"Aligning took: {time.time()-time_clustalo:.4f}s")
    
    # Step 4: Filter by gaps
    if verbose: print(f">>> Start filtering by gaps ({len(clusters)} clusters left)")
    time_filter_gaps =  time.time()
    clusters = fl.filter_gaps(clusters, threshold=gaps_threshold, absolute=False, average=True)
    if timing: print(f"Gap filtering took: {time.time()-time_filter_gaps:.4f}s")

    # Step 5: Filter by UniRef ID
    if uniref_lookup:
        if verbose: print(f">>> Start filtering by UniRefIDs ({len(clusters)} clusters left)")
        time_filter_uniref = time.time()
        paths = [row[0] for row in io.parse_csv(uniref_lookup)]
        lookup = bt.Bakta_table()
        lookup.read(paths, skip=5)
        # 3.3.1 UniRef100
        if verbose: print(f">>> Start filtering by UniRef100 IDs ({len(clusters)} clusters left)")
        time_filter_uniref_100 = time.time()
        clusters = fl.filter_uniref(clusters, lookup, uniref100_threshold, level=100)
        if timing: print(f"UniRef100 filtering took {time.time()-timed-filter_uniref_100:.4f}s")
        # 3.3.2 UniRef90
        if verbose: print(f">>> Start filtering by UniRef90 IDs ({len(clusters)} clusters left)")
        time_filter_uniref_90 = time.time()
        clusters = fl.filter_uniref(clusters, lookup, uniref90_threshold, level=90)
        if timing: print(f"UniRef90 filtering took {time.time()-timed-filter_uniref_90:.4f}s")
        # 3.3.3 UniRef50
        if verbose: print(f">>> Start filtering by UniRef50 IDs ({len(clusters)} clusters left)")
        time_filter_uniref_50 = time.time()
        clusters = fl.filter_uniref(clusters, lookup, uniref50_threshold, level=50)
        if timing: print(f"UniRef50 filtering took {time.time()-timed-filter_uniref_50:.4f}s")
        if timing: print(f"UniRef filtering took {time.time()-timed-filter_uniref:.4f}s")

    # Step 6: Filter by Length
    if verbose: print(f">>> Start filtering by length ({len(clusters)} clusters left)")
    time_filter_length = time.time()
    clusters = fl.filter_length(clusters, threshold=length_threshold)
    if timing: print(f"Length filtering took: {time.time()-time_filter_length:.4f}s")

    if verbose: print(f">>> Filtering done ({len(clusters)} Clusters left)")
    
    # Step 7: Saving Alginments
    if alignment_path:
        for index, cluster in enumerate(clusters):
            cluster.write(os.path.join(alignment_path, f"{index:0{len(str(len(clusters)))}}"))
    if verbose: print(f">>> Alignments have been saved to {alignment_path}")

    # Step 8: Trees
    if verbose: print(">>> Start calculating trees")
    time_trees = time.time()
    trees = [cluster.distmat.upgma() for cluster in clusters]
    if timing: print(f"Calculating Trees took: {time.time()-time_trees:.4f}s")

    # Step 9: Output
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
        "--diamond_path",
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
        "--timing",
        action = "store_true",
        help = "Whether to show detailed timing information"
    )
    parser.add_argument(
        "--size",
        metavar = "CLUSTER_SIZE",
        help = "The minimal cluster size to keep, defaults to 3",
        type = int
    )
    parser.add_argument(
        "--gaps",
        metavar = "GAP_RATE",
        help = "The max gap rate to keep",
        type = float
    )
    parser.add_argument(
        "--length",
        metavar = "LENGTH_DIFFERENCE",
        help = "The max difference in sequence length (relative) to keep",
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
    DIAMOND = "./diamond/diamond" if not args.diamond_path else args.diamond_path
    executable = DIAMOND

    params = {"data_file": args.DATA}
    params["executable"] = executable
    params["threshold"] = THRESHOLD
    if args.alignment_out: params["alignment_path"] = args.alignment_out
    if args.linclust: params["method"] = "linclust"
    if args.verbose: params["verbose"] = args.verbose
    if args.timing: params["timing"] = args.timing
    if args.size: params["size_threshold"] = args.size
    if args.gaps: params["gaps_threshold"] = args.gaps
    if args.length: params["length_threshold"] = args.length
    if args.lookup: params["uniref_lookup"] = args.lookup
    if args.ur50: params["uniref50_threshold"] = args.ur50
    if args.ur90: params["uniref90_threshold"] = args.ur90
    if args.ur100: params["uniref100_threshold"] = args.ur100
    if args.out: params["out_file"] = args.out

    main(**params)
