# vim: set foldmethod=marker:
# vim: set foldclose=all foldlevel=0:
# vim: set foldenable: 

from typing import List, Set, Optional, Union
import re

import fasta as fs
import io_helpers as io

class Bakta_table: #{{{
    table: List[dict]
    index: int

    def __init__( #{{{
        self, 
        table: List[dict] = []
    ):
        """
        Create a Bakta_table object
        Args:
            table (List[dict]): A list of dictionaries parsed from a bakta tsv. Defaults to an empty list.
        Returns:
            None
        """
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
        self.index = 0
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
        filepath:Union[str, List[str]],
        sep:str = "\t",
        skip:int = 0
    ):
        """
        Read a bakta tsv into the table, expanding it
        Args:
            filepath (str | List[str]): Filepath(s) to be read
            sep (str): The separator for the file(s). Defaults to '\\t'.
            skip (int): The number of lines to skip before the header row. Defaults to 0.
        Returns:
            None
        """
        if isinstance(filepath, str):
            filepath = [filepath]
        for entry in filepath:
            file = io.parse_csv(
                filepath = entry,
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
        """
        Lookup an entry in the Bakta_table object
        Args:
            value (str): The value to search for.
            key (str, List[str]): The column(s) to look in. Omitting causes the method to look through all columns.
        Returns:
            List[dict]: A list of dictionaries, each being one entry from the Bakta_table object where <value> was found.
        """
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

    def get_uniref( #{{{
        self,
        locus_tag: str,
        level: int
    ) -> Optional[str]:
        """
        A more specialized version of the find method to look specifically for unirefids
        Args:
            locus_tag (str): The locus_tag (identifier) to look for
            level (int): The UniRef ID level to look for (e.g 90).
        Returns:
            str or None: Returns the ID as string or None if no ID was found.
        """
        # Find the entry based on the locus tag
        entry = self.find(locus_tag, "locus tag")

        # Ensure there's only one entry
        if len(entry) > 1:
            raise ValueError(f"Multiple entries for <{locus_tag}>")
        elif len(entry) == 0:
            return None  # No entry found
        
        # Extract the single entry
        entry = entry[0]
        
        # Define the regex pattern for the UniRef ID
        pattern = rf"UniRef:UniRef{level}_([a-zA-Z0-9]+)"
        
        # Check if the 'locus tag' field exists and extract the UniRef ID
        locus_tag_field = entry.get('dbxrefs')
        if locus_tag_field:
            match = re.search(pattern, locus_tag_field)
            if match:
                return match.group(1)
        
        # Return None if no match is found
        return None
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

    bakta_table.read("../data/bin.2/bin.2.tsv", skip=5)
    [print(bakta_table[i]) for i in range(0,20)]
    print(len(bakta_table))
    print(bakta_table["DbXrefs", "SO:0001217, UniRef:UniRef50_UPI00260DA7F0"])
    print("\n")
    annotation = bakta_table.find("OJFFMF_00010", "Locus Tag")[0]
    print(annotation)
    match = get_uniprot(annotation, level=50)
    print(match)
