# vim: set foldmethod=marker:
# vim: set foldclose=all foldlevel=0:
# vim: set foldenable: 

from __future__ import annotations
from typing import List, Tuple, Union
from io import StringIO
import re
import subprocess
import os
import copy

import io_helpers as io

class Sequence: #{{{
    header: str
    sequence: str

    def __init__(self, header: str, sequence: str):
        """
        Create a Sequence object
        Args:
            header (str): The sequence header used to identify the sequence
            sequence (str): The sequence itself.
        """
        self.header = header
        self.sequence = sequence
        self.index = 0

    def __str__(self):
        return self.sequence

    def __repr__(self):
        return "\n".join([str(self.header), str(self.sequence)])

    def __len__(self):
        return len(self.sequence)

    def __lt__(self, other):
        if not isinstance(other, Sequence):
            return NotImplemented
        return len(self) < len(other)

    def __eq__(self, other):
        return self.sequence == other.sequence

    def __iter__(self):
        return self

    def __next__(self):
        if self.index >= len(self.sequence):
            raise StopIteration
        char = self.sequence[self.index]
        self.index += 1
        return char

    def redit(self, edit: Tuple[str, str], field: str = "header"):
        """
        Edit the sequence using regular expressions
        Args:
            edit (Tuple[str, str]): A tuple containing a find and a replace regular expression.
            field (str): The field to edit. Can be either 'header' or 'sequence'. Defaults to 'header'.
        Returns:
            None
        Raises:
            ValueError: If a value other that 'header' or 'sequence' is provided for <field>.
        """
        if field == "header":
            self.header = re.sub(edit[0], edit[1], self.header)
        elif field == "sequence":
            self.sequence = re.sub(edit[0], edit[1], self.sequence)
        else:
            raise ValueError(f"Valid fields: header, sequence (provided: {field})")

    def count(self, symbol:str="-", relative:bool=False) -> Union[int,float]:
        """
        Count the occurences of a specified symbol (or substring) in the sequence
        Args:
            symbol (str): The symbol in question. Defaults to '-'.
            relative (bool): Whether to return the fraction of sequence occupied by <symbol>. Defaults to False.
        Returns:
            int or float: Either the absolute count of <symbol> or the fraction of <symbol> in the sequence.
        """
        total = len([occurence for occurence in self.sequence if occurence == symbol])
        return total if not relative else total/len(self.sequence)

#}}}

class Fasta: #{{{
    sequences: List[Sequence]
    distmat: Distmat

    def __init__( #{{{
        self,
        sequences: List[Sequence] = None
    ) -> None:
        """
        Create a Fasta object
        Args:
            sequences (List[Sequence]): A list of Sequence objects. Defaults to None.
        Returns:
            None
        """
        self.distmat = None
        if sequences is None:
            self.sequences = []
        else:
            self.sequences = sequences
        self.index = 0
    #}}}

    def __str__( #{{{
        self
    ):
        return "\n".join(repr(sequence) for sequence in self.sequences)
    #}}}

    def __getitem__( #{{{
        self,
        key
    ):
        return self.sequences[key]
    #}}}

    def __len__( #{{{
        self
    ) -> int:
        return len(self.sequences)
    #}}}

    def __iter__( #{{{
        self
    ):
        self.index = 0
        return self
    #}}}

    def __next__( #{{{
        self
    ):
        if self.index >= len(self.sequences):
            raise StopIteration
        sequence = self.sequences[self.index]
        self.index += 1
        return sequence
    #}}}

    def add( #{{{
        self,
        new: Union['Sequence', 'Fasta', List['Sequence']]
    ) -> None:
        """
        Add (a) new Sequence(s) to the Fasta object.
        Args:
            new (Sequence | Fast | List[Sequence]): The content to be added to the Fasta object.
        Returns:
            None
        """
        if isinstance(new, Sequence):
            self.sequences.append(new)
        elif isinstance(new, list):
            for sequence in new:
                self.add(sequence)
        elif isinstance(new, Fasta):
            self.add(new.sequences)
    #}}}

    def delete( #{{{
        self,
        index: int
    ) -> None:
        """
        Delete a sequence from the Fasta object by index.
        Args:
            index (int): The index of the sequence to delete.
        Returns:
            None
        """
        del self.sequences[index]
    #}}}

    def get( #{{{
        self,
        index:int
    ) -> Sequence:
        """
        Get a sequence from the Fasta object by index:
        Args:
            index (int): The index of the Sequence object to return.
        Returns:
            Sequence: The Sequence object at the specified index.
        """
        return self.sequences[index]
    #}}}

    def search( #{{{
        self,
        search:str,
        in_seq:bool=False,
        regex:bool=False
    ) -> List[Sequence]:
        """
        Search for sequences inside the Fasta object.
        Args:
            search (str): The searchstring.
            in_seq (bool): Whether to search inside the sequence (instead of its header). Defaults to False.
            regex (bool): Whether to understand <search> as regular expression. Defaults to False.
        Return:
            List[Sequence]: A list of all found Sequences matching the search criteria.
        """
        result = []
        for sequence in self.sequences:
            if not regex:
                if (search in sequence.sequence if in_seq else search in sequence.header):
                    result.append(sequence)
            else:
                if re.search(search, sequence.sequence if in_seq else sequence.header):
                    result.append(sequence)
        return result
    #}}}

    def redit( #{{{
        self,
        edit: Tuple[str, str],
        field: str = "header"
    ) -> None:
        """
        Edit the Sequences in the Fasta object using regular expressions.
        This just calls the redit function of every single Sequence in the Fasta object
        Args:
            edit (Tuple[str, str]): A tuple containting a search and a replace expression.
            field (str): The field to edit. Either 'header' or 'sequence'. Defaults to 'header'.
        Returns:
            None
        """
        for sequence in self.sequences:
            sequence.redit(edit=edit, field=field)
    #}}}

    def write( #{{{
        self,
        filename: str,
        line_length=0
    ) -> None:
        """
        Write the Fasta object to a fasta file.
        Args:
            filename (str): The filename of the file to be written. The extension will always be '.fasta'.
            line_length (int): The maxmimum line length before breaking when writing the fasta file. Omitting will disable linebreaks.
        Returns:
            None
        """
        write_string = ""
        for seq in self.sequences:
            header = seq.header
            sequence_string = seq.sequence
            if line_length > 0:
                split_sequence_lines = []
                for i in range(0, len(sequence_string), line_length):
                    chunk = sequence_string[i:i + line_length]
                    split_sequence_lines.append(chunk)
                split_sequence = "\n".join(split_sequence_lines)
            else:
                split_sequence = sequence_string
            write_string += f"{header}\n{split_sequence}\n"
        io.write_file(f"{filename}.fasta", content=write_string)
    #}}}

    def read( #{{{
        self,
        input_file: str,
        identifier:str = "",
        from_file:bool = True
    ) -> None:
        """
        Read a fasta file into the Fasta object. If the the object is non-empty, the file content gets appended to the current object.
        Args:
            input_file (str): The path to the input file.
            identifier (str): An identifier to mark read headers with (gets prepended). Empty by default.
            from_file (bool): Wheter <input_file> refers to a file or contains the fasta string directly. Defaults to True.
        """
        if from_file:
            file = io.read_file(input_file, lines=True)
        else:
            file = [line.strip() for line in input_file.splitlines()]
        current_header = None
        current_sequence = ""
        for line in file:
            if line.startswith(">"):
                if current_header is not None:
                    current_header = f">{identifier}{current_header[1:]}"
                    sequence = Sequence(header=current_header, sequence=current_sequence)
                    self.sequences.append(sequence)
                current_header = line.strip()
                current_sequence = ""
            else:
                current_sequence += line.strip()
        if current_header is not None:
            current_header = f">{identifier}{current_header[1:]}"
            sequence = Sequence(header=current_header, sequence=current_sequence)
            self.sequences.append(sequence)
    #}}}

    def align( #{{{
        self,
    ) -> Fasta:
        """
        Align the Fasta object using clustalo subprocess.
        Does not overwrite the original object.
        Args:
            None
        Returns:
            Fasta: A new Fasta object with aligned sequences.
        """
        command = ["clustalo", "-i", "-"]
        process = subprocess.Popen(
            command,
            stdin = subprocess.PIPE,
            stdout = subprocess.PIPE,
            stderr = subprocess.PIPE,
        )
        stdout, _ = process.communicate(input=str(self).encode())

        result = Fasta()
        result.read(
            input_file=stdout.decode("utf-8").replace("\\n", "\n"),
            from_file=False
        )
        return result
    #}}}

    def distmat( #{{{
        self
    ) -> Distmat:
        """
        Get the Fasta objects distance matrix.
        Should only be called if the matrix has been created before.
        Args:
            None
        Returns:
            Distmat: The distance matrix for this Fasta object.
        """
        return Distmat(self)
    #}}}

    def count( #{{{
        self,
        symbol: str = '-',
        absolute:bool = False,
        average:bool = True,
    ) -> int | float:
        """
        Count occuences of a specified symbol over all sequences.
        Args:
            symbol (str): The symbol (or substring) in question. Defaults to '-'.
            absolute (bool): Whether to count absolute occuences or relative occurences. Defaults to False.
            average (bool): Whether to return summed or averaged occuence counts. Obsolete if <absolute> is set to False. Defaults to True.
        Returns:
            int or float: Either the absolute sum or average of sequences or the average ratio over all sequences.
        """
        if absolute and not average:
            return sum([sequence.count(symbol=symbol, relative=False) for sequence in self.sequences])
        elif absolute and average:
            return sum([sequence.count(symbol=symbol, relative=False) for sequence in self.sequences])/len(self)
        elif not absolute:
            return sum([sequence.count(symbol=symbol, relative=True) for sequence in self.sequences])/len(self)
    #}}}

    def cd( #{{{
        self
    ) -> None:
        """
        [C]alculate a [d]istance matrix for this Fasta object.
        Args:
            None
        Returns:
            None
        """
        self.distmat = Distmat(self)
    #}}}

    def clustalo( #{{{
        self
    ) -> Fasta:
        """
        Combination of the align and cd methods.
        Creates a new aligned Fasta object with calculated distance matrix.
        Does in no way overwrite the current Fasta object.
        Args:
            None
        Returns:
            Fasta: An aligned Fasta object with distance matrix attribute.
        """
        result = Fasta()
        # Alignment
        command = ["clustalo", "--full", "--force", "--distmat-out=matrix.temp", "-i", "-"]
        process = subprocess.Popen(
            command,
            stdin = subprocess.PIPE,
            stdout = subprocess.PIPE,
            stderr = subprocess.PIPE,
        )
        stdout, _ = process.communicate(input=str(self).encode())
        result.read(
            input_file=stdout.decode("utf-8").replace("\\n", "\n"),
            from_file=False
        )

        # Distance Matrix
        labels = []
        matrix = []

        file = io.read_file("matrix.temp")
        for line in file.splitlines()[1:len(result)+1]:
            parts = line.split()
            labels.append(parts[0])
            matrix.append([float(x) for x in parts[1:]])

        result.distmat = Distmat(matrix=matrix, labels=labels)

        try:
            os.remove("matrix.temp")
        except Exception as e:
            print(f"There was an error removing temporary files: {e}")

        return result
    #}}}
#}}}

class Distmat: #{{{
    matrix: List[List[float]]
    labels: List[str]

    def __init__( #{{{
        self,
        matrix,
        labels=None
    ) -> None:
        """
        Create a distance matrix object
        Args:
            matrix (int or Fasta, List[List]): The source of the matrix. Can either be an integer (specifying the size for an empty matrix), a Fasta to create it from or a predefined one.
            labels (List[str]): A list of labels for the elements in the matrix. Calculated automatically if <matrix> is a Fasta object. Omitting uses Letters (A, B, ..., AA, AB, ...).
        Returns:
            None
         """
        # Check if the matrix is a tuple (x, y) and initialize a zero matrix
        if isinstance(matrix, int):
            if matrix <= 0:
                raise ValueError("Matrix dimension must be a positive integer.")
            matrix = [[0.0 for _ in range(matrix)] for _ in range(matrix)]
        elif isinstance(matrix, Fasta):
            self.matrix, self.labels = self._from_fasta(matrix)
            return
        else:
            # Check if the matrix is quadratic (number of rows == number of columns)
            if not all(len(row) == len(matrix) for row in matrix):
                raise ValueError("The given matrix is not quadratic (must have the same number of rows and columns).")
        
        self.matrix = matrix

        # Generate labels if none are provided
        if labels is None:
            size = len(matrix)
            self.labels = self._generate_labels(size)
        else:
            if len(labels) != len(matrix):
                raise ValueError("The number of labels must match the size of the matrix.")
            self.labels = labels
    #}}}

    def _from_fasta( #{{{
        self,
        fasta
    ):
        command = ["clustalo", "--full", "--force", "--distmat-out=/dev/stdout", "-o", "/dev/null", "-i", "-"]
        process = subprocess.Popen(
            command,
            stdin = subprocess.PIPE,
            stdout = subprocess.PIPE,
            stderr = subprocess.PIPE,
        )
        stdout, _ = process.communicate(input=str(fasta).encode())

        labels = []
        matrix = []

        for line in stdout.decode("utf-8").splitlines()[1:len(fasta)+1]:
            parts = line.split()
            labels.append(parts[0])
            matrix.append([float(x) for x in parts[1:]])

        return matrix, labels
    #}}}

    def _generate_labels( #{{{
        self,
        size
    ):
        labels = []
        alphabet = "ABCDEFGHIJKLMNOPQRSTUVWXYZ"
        length = 1

        while len(labels) < size:
            for i in range(len(alphabet) ** length):
                label = ""
                num = i
                for _ in range(length):
                    label = alphabet[num % len(alphabet)] + label
                    num //= len(alphabet)
                labels.append(label)
                if len(labels) == size:
                    break
            length += 1

        return labels
    #}}}

    def __str__( #{{{
        self
    ):
        max_val = max(max(row) for row in self.matrix)
        value_padding = len(str(max_val)) + 1

        label_padding = max(len(label) for label in self.labels)

        padding = max(value_padding, label_padding)

        header = " " * (label_padding + 1) + " ".join(f"{label:>{padding}}" for label in self.labels)

        rows = []
        for label, row in zip(self.labels, self.matrix):
            row_str = " ".join(f"{val:>{padding}}" for val in row)
            rows.append(f"{label:<{label_padding}} {row_str}")

        return f"{header}\n" + "\n".join(rows)
    #}}}

    def __iter__( #{{{
        self
    ):
        self._iterator = (item for row in self.matrix for item in row)
        return self._iterator
    #}}}

    def __next__( #{{{
        self
    ):
        if self._iterator is None:
            self.__iter__()
        return next(self._iterator)
    #}}}

    def __len__( #{{{
        self
    ):
        return sum(len(row) for row in self.matrix)
    #}}}

    def __eq__( #{{{
        self,
        other
    ):
        if not isinstance(other, Distmat):
            return False
        return self.matrix == other.matrix
    #}}}

    def __lt__( #{{{
        self,
        other
    ):
        if not isinstance(other, Distmat):
            return NotImplemented
        return sum(sum(row) for row in self.matrix) < sum(sum(row) for row in other.matrix)
    #}}}

    def smallest( #{{{
        self,
        matrix:List[List[float]]
    ) -> Tuple[int, int]:
        """
        Return the smallest value from a matrix. Excludes self-distances (the diagonal).
        Does not work on the object it's attached to but rather a provided matrix.
        Args:
            matrix (List[List[float]]): The matrix to work on.
        Returns:
            Tuple[int, int]: The coordinates of the smallest value inside the matrix.
        """
        result = (-1, -1)
        current_min = float("inf")
        for i, row in enumerate(matrix):
            for j, col in enumerate(row):
                if col < current_min and i != j:
                    current_min = col
                    result = (i, j)
        return result
    #}}}

    def _join_cells( #{{{
            _,
            matrix: List[List[float]],
            labels: List[str],
            a: int,
            b: int,
            distances: bool = False
    ):
        if b < a: 
            a, b = b, a

        # Calculate branch lengths for a and b
        branch_length_a = matrix[a][b] / 2
        branch_length_b = matrix[a][b] / 2

        # Create a new label for the merged cluster
        if distances:
            new_label = f"({labels[a]}:{branch_length_a},{labels[b]}:{branch_length_b})"
        else:
            new_label = f"({labels[a]},{labels[b]})"

        # Create a new labels list to avoid modifying the original
        new_labels = labels[:]
        new_labels.append(new_label)
        del new_labels[b]
        del new_labels[a]

        a_dists = matrix[a]
        b_dists = matrix[b]
        new_dists = []

        # Row
        for i in range(0, len(a_dists)):
            new_dists.append((a_dists[i] + b_dists[i]) / 2)
        new_dists.append(0.0)
        matrix.append(new_dists)

        # Column
        for i in range(0, len(matrix) - 1):
            matrix[i].append(new_dists[i])

        # Deleting
        for row in range(0, len(matrix)):
            del matrix[row][b]
            del matrix[row][a]
        del matrix[b]
        del matrix[a]
        return matrix, new_labels
    #}}}

    def upgma( #{{{
            self,
            distances: bool = False,
    ) -> str:
        """
        Perform UPGMA on the matrix, creating a tree in Newick format.
        Args:
            distances (bool): Whether to include distances. Currently not implemented. Defaults to False.
        Returns:
            str: The string of the created tree in Newick format.
        """
        # Create deep copies of matrix and labels
        matrix = copy.deepcopy(self.matrix)
        labels = copy.deepcopy(self.labels)
        
        while len(labels) > 1:
            a, b = self.smallest(matrix)
            matrix, labels = self._join_cells(matrix, labels, a, b, distances)
        return labels[0]
    #}}}

#}}}
