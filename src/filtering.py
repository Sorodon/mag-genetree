# vim: set foldmethod=marker:
# vim: set foldclose=all foldlevel=0:
# vim: set foldenable: 

from typing import List, Set, Union, Tuple

import fasta as fs
import bakta_table as bt

## Number of Members
def filter_size( #{{{
    clusters: List[fs.Fasta],
    threshold:int
    ) -> List[fs.Fasta]:
    """
    Filter a list of clusters (Fasta objects) by their size.
    Args:
        clusters (List[Fasta]): The list to filter.
        threshold (int): The minimal size to keep.
    Returns:
        List[Fasta]: A list of all Fasta objects that passed the filter.
    """
    return [fasta for fasta in clusters if len(fasta)>=threshold]
#}}}

## Gap count (absolute and fraction)
def filter_gaps( #{{{
    clusters:List[fs.Fasta],
    threshold:int,
    absolute:bool = False,
    average:bool = True
) -> List[fs.Fasta]:
    """
    Filter a list of clusters (Fasta objects) by their gaps.
    Args:
        clusters (List[Fasta]): The list to filter.
        threshold (int): The maximal number of gaps to keep.
        absolute (bool): Whether <threshold> is to be understood as an absolute value. Defaults to False.
        averarage (bool): Whether <threshold> is to be understood as the average across sequences instead of a sum (only for <absolute>=True). Defaults to True.
    Returns:
        List[Fasta]: A list of all Fasta objects that passed the filter.
    """
    return [fasta for fasta in clusters if fasta.count(symbol="-", absolute=absolute, average=average)<=threshold]
#}}}

## UniRefID
def filter_uniref( #{{{
    clusters:List[fs.Fasta],
    lookup: bt.Bakta_table,
    threshold:int = 1,
    level:int = 100,
    stats:bool = False,
    accept_missing:bool = True,
    sep:str = "-"
) -> List[fs.Fasta] | Tuple[List[fs.Fasta], dict]:
    """
    Filter a list of clusters (Fasta objects) by their number of distict UniRef IDs.
    Args:
        clusters (List[Fasta]): The list to filter.
        lookup (Bakta_table): A Bakta_table object containing all concerned sequences.
        threshold (int): The maximal allowed number of distict IDs per cluster. Defaults to 1.
        level (int): The UniRef ID level to work with (e.g. 90). Defaults to 100.
        stats (bool): Whether to output a dictionary with stats in addition to the filtered list. Defaults to False.
        accept_missing (bool): Whether to treat missing IDs as non-distict. Defaults to True.
        sep (str): The separator between identifier and locus tag in the sequence header if their is one (set to '' if their is none). Defaults to '-'.
    Returns:
        List[Fasta] or Tuple[List[Fasta], dict]: A list of all Fasta objects that passed the filter or a Tuple containing the former and a stat dictionary.
    """
    result = []
    stat_dict = {}
    for cluster in clusters:
        ids = set()
        missing_id = False
        for seq in cluster:
            if sep:
                _, _, tag = seq.header.partition(sep)
            else:
                tag = seq.header
            id = lookup.get_uniref(tag.split(' ', 1)[0], level)
            if id:
                ids.add(id)
            else:
                missing_id = True
        if len(ids) <= threshold and (accept_missing or not missing_id):
            result.append(cluster)
        if missing_id:
            stat_dict[f"{len(ids)}m"] = stat_dict.get(f"{len(ids)}m", 0) + 1
        else:
            stat_dict[f"{len(ids)}"] = stat_dict.get(f"{len(ids)}", 0) + 1
    if stats:
        return (result, stat_dict)
    else:
        return result
#}}}

## Filter Length difference
def filter_length( #{{{
    clusters: List[fs.Fasta],
    threshold:float
    ) -> List[fs.Fasta]:
    """
    Filter a list of clusters (Fasta objects) by their length.
    Args:
        clusters (List[Fasta]): The list to filter.
        threshold (int): The maximal percentage difference in length to keep.
    Returns:
        List[Fasta]: A list of all Fasta objects that passed the filter.
    """
    result = []
    for fasta in clusters:
        ratio = len(min(fasta.sequences)) / len(max(fasta.sequences))
        if ratio >= 1-threshold:
            result.append(fasta)
    return result
#}}}
