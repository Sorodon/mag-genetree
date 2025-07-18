# vim: set foldmethod=marker:
# vim: set foldclose=all foldlevel=0:
# vim: set foldenable: 

# TODO: add pydoc string
from __future__ import annotations
from typing import List, Tuple, Union
from Bio import Phylo
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
        if field == "header":
            self.header = re.sub(edit[0], edit[1], self.header)
        elif field == "sequence":
            self.sequence = re.sub(edit[0], edit[1], self.sequence)
        else:
            raise ValueError(f"Valid fields: header, sequence (provided: {field})")
#}}}

class Fasta: #{{{
    sequences: List[Sequence]
    distmat: Distmat

    def __init__(self, sequences: List[Sequence] = None):
        if sequences is None:
            self.sequences = []
            self.distmat = None
        else:
            self.sequences = sequences
            self.distmat = Distmat(self)
        self.index = 0

    def __str__(self):
        return "\n".join(repr(sequence) for sequence in self.sequences)

    def __getitem__(
        self,
        key
    ):
        return self.sequences[key]

    def __len__(self):
        return len(self.sequences)

    def __iter__(self):
        self.index = 0
        return self

    def __next__(self):
        if self.index >= len(self.sequences):
            raise StopIteration
        sequence = self.sequences[self.index]
        self.index += 1
        return sequence

    def add(self, new: Union['Sequence', 'Fasta', List['Sequence']]) -> None:
        if isinstance(new, Sequence):
            self.sequences.append(new)
        elif isinstance(new, list):
            for sequence in new:
                self.add(sequence)
        elif isinstance(new, Fasta):
            self.add(new.sequences)
        self.distmat = Distmat(self)


    def delete(self, index: int) -> None:
        del self.sequences[index]
        self.distmat = Distmat(self)

    def get(self, index:int) -> Sequence:
        return self.sequences[index]

    def search(self, search:str, in_seq:bool=False, regex:bool=False) -> List[Sequence]:
        result = []
        for sequence in self.sequences:
            if not regex:
                if (search in sequence.sequence if in_seq else search in sequence.header):
                    result.append(sequence)
            else:
                if re.search(search, sequence.sequence if in_seq else sequence.header):
                    result.append(sequence)
        return result

    def redit(self, edit: Tuple[str, str], field: str = "header"):
        for sequence in self.sequences:
            sequence.redit(edit=edit, field=field)
        self.distmat = Distmat(self)

    def write(self, filename: str, line_length=0) -> None:
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


    def read(self, input_file: str, identifier:str = "", from_file:bool = True) -> None:
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
        self.distmat = Distmat(self)

    def align( #{{{
        self,
    ) -> Fasta:
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
    ):
        return Distmat(self)
    #}}}

    def count( #{{{
        self,
        symbol: str = '-',
        stats: bool = False
    ) -> int | list[float]:
        if not stats:
            # Count occurrences of `symbol` across all sequences
            amount = 0
            for sequence in self.sequences:
                amount += sequence.sequence.count(symbol)
            return amount
        else:
            # Determine the length of the shortest sequence
            if not self.sequences:
                return []

            min_length = min(len(sequence.sequence) for sequence in self.sequences)
            num_sequences = len(self.sequences)
            fractions = []

            # Count occurrences of `symbol` for each column up to the shortest sequence length
            for col in range(min_length):
                count_symbol = sum(
                    1 for sequence in self.sequences if sequence.sequence[col] == symbol
                )
                fractions.append(count_symbol / num_sequences)

            return fractions
    #}}}

#}}}

class Distmat: #{{{
    matrix: List[List[float]]
    labels: List[str]

    def __init__(self, matrix, labels=None):
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

    def _from_fasta(self, fasta):
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

        for line in stdout.decode("utf-8").splitlines()[1:len(fasta)]:
            parts = line.split()
            labels.append(parts[0])
            matrix.append([float(x) for x in parts[1:]])

        return (matrix, labels)

    def _generate_labels(self, size):
        # Generate labels like A, B, ..., Z, AA, AB, ..., ZZ, AAA, ...
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

    def __str__(self):
        # Determine the padding dynamically based on the largest number in the matrix
        max_val = max(max(row) for row in self.matrix)
        value_padding = len(str(max_val)) + 1

        # Determine the padding based on the longest label
        label_padding = max(len(label) for label in self.labels)

        # Final padding is the maximum of value_padding and label_padding
        padding = max(value_padding, label_padding)

        # Create the header with column labels
        header = " " * (label_padding + 1) + " ".join(f"{label:>{padding}}" for label in self.labels)

        # Create rows with row labels and values
        rows = []
        for label, row in zip(self.labels, self.matrix):
            row_str = " ".join(f"{val:>{padding}}" for val in row)
            rows.append(f"{label:<{label_padding}} {row_str}")

        # Combine header and rows into the final string
        return f"{header}\n" + "\n".join(rows)

    def __iter__(self):
    
        self._iterator = (item for row in self.matrix for item in row)
        return self._iterator

    def __next__(self):
        if self._iterator is None:
            self.__iter__()
        return next(self._iterator)

    def __len__(self):
        return sum(len(row) for row in self.matrix)

    def __eq__(self, other):
        if not isinstance(other, Distmat):
            return False
        return self.matrix == other.matrix

    def __lt__(self, other):
        if not isinstance(other, Distmat):
            return NotImplemented
        return sum(sum(row) for row in self.matrix) < sum(sum(row) for row in other.matrix)

    def smallest(self, matrix:List[List[float]]) -> Tuple[int, int]:
        result = (-1, -1)
        current_min = float("inf")
        for i, row in enumerate(matrix):
            for j, col in enumerate(row):
                if col < current_min and i != j:
                    current_min = col
                    result = (i, j)
        return result

    def _join_cells(
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

    def upgma(
            self,
            distances: bool = False,
            draw:bool = False
    ) -> str:
        # Create deep copies of matrix and labels
        matrix = copy.deepcopy(self.matrix)
        labels = copy.deepcopy(self.labels)
        
        while len(labels) > 1:
            a, b = self.smallest(matrix)
            matrix, labels = self._join_cells(matrix, labels, a, b, distances)
        if draw:
            Phylo.draw(Phylo.read(StringIO(labels[0]), "newick"))
        else:
            return labels[0]

#}}}
