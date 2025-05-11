# TODO: add docstrings

from typing import List

import fasta as fs

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

