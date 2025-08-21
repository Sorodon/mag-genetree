# vim: set foldmethod=marker:
# vim: set foldclose=all foldlevel=0:
# vim: set foldenable: 

# TODO: add pydoc string
from typing import List, Set, Union, Tuple

import fasta as fs
import bakta_table as bt

## Number of Members
def filter_size( #{{{
    clusters: List[fs.Fasta],
    threshold:int
    ) -> List[fs.Fasta]:
    return [fasta for fasta in clusters if len(fasta)>=threshold]
#}}}

## Gap count (absolute and fraction)
def filter_gaps( #{{{
    clusters:List[fs.Fasta],
    threshold:int,
    absolute:bool = False,
    average:bool = True
) -> List[fs.Fasta]:
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
