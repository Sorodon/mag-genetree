# vim: set foldmethod=marker:
# vim: set foldclose=all foldlevel=0:
# vim: set foldenable: 

# TODO: add docstrings

from typing import List, Set, Optional, Union

import fasta as fs
import io_helpers as io

class Bakta_table: #{{{
    table: List[dict]
    index: int

    def __init__( #{{{
        self, 
        table: List[dict] = []
    ):
        self.table = table
        self.index = 0
    #}}}

    def __getitem__( #{{{
        self,
        args: Union[tuple, int, slice]
    ) -> list[dict]:
        if isinstance(args, tuple):
            key, value = args
            result = []
            for entry in self.table:
                if key is not None and entry.get(key) == value:
                    result.append(entry)
                elif value in entry.values():
                    result.append(entry)
            return result
        elif isinstance(args, int):
            return self.table[args]
        elif isinstance(args, slice):
            return self.table[args]
    #}}}

    def __str__( #{{{
        self
    ):
        return str(self.table)
    #}}}

    def __iter__(self): #{{{
        return self
    #}}}

    def __next__(self): #{{{
        if self.index < len(self.table):
            value = self.table[self.index]
            self.index += 1
            return value
        else:
            raise StopIteration
    #}}}

    def __len__(self): #{{{
        return len(self.table)
    #}}}

    def read( #{{{
        self, 
        filepath,
        sep = "\t",
        skip = 0
    ):
        file = io.parse_csv(
            filepath = filepath,
            sep = sep,
            header_row = True,
            skip = skip
        )
        self.table = self.table + file
    #}}}

    def find( #{{{
        self,
        value: str,
        key: str = None
    ) -> List[dict]:
        if key is None:
            return [entry for entry in self.table if any(value in v for v in entry.values())]
        elif isinstance(key, list):
            result = []
            for entry in self.table:
                for keyx, valuex in entry.items():
                    if value in valuex and keyx in key:
                        result.append(entry)
            return result
        else:
            result = []
            for entry in self.table:
                for keyx, valuex in entry.items():
                    if value in valuex and keyx in key:
                        result.append(entry)
            return result
    #}}}
#}}}

if __name__ == "__main__":
    bakta = {
        "sequence_id": "sequence_id",
        "typ": "typ",
        "start": "start",
        "stop": "stop",
        "strand": "+",
        "locus": "locus",
        "tag": "tag",
        "gene": "gene",
        "product": "product",
        "dbxrefs": "dbxrefs"
    }
    bakta_table = Bakta_table()

    bakta_table.read("../data/bin.1/bin.1.tsv", skip=5)
    [print(bakta_table[i]) for i in range(0,20)]
    print(len(bakta_table))
    print(bakta_table["DbXrefs", "SO:0001217, UniRef:UniRef50_UPI00260DA7F0"])
    print("\n")
    print(bakta_table.find("UniRef:UniRef50_UPI00260CD7CB"))
    print(bakta_table.find("UniRef:UniRef50_UPI00260CD7CB", key=["DbXrefs"]))
    print(bakta_table.find("UniRef:UniRef50_UPI00260CD7CB", "DbXrefs"))
