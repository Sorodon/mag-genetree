# TODO: add pydoc string

import io_helpers
from typing import List

class Sequence:
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

    def sequence(self) -> str:
        return self.sequence

    def header(self) -> str:
        return self.header

class Fasta:
    sequences: List[Sequence]

    def __init__(self, sequences: List[Sequence]):
        self.sequences = sequences
        self.index = 0

    def __str__(self):
        return "\n".join(repr(sequence) for sequence in self.sequences)

    def __iter__(self):
        return self

    def __next__(self):
        if self.index >= len(self.sequences):
            raise StopIteration
        sequence = self.sequences[self.index]
        self.index += 1
        return sequence

    def add(self, sequence: Sequence) -> None:
        self.sequences.append(sequence)

    def delete(self, index: int) -> None:
        del self.sequences[index]

    def get(self, index:int) -> Sequence:
        return self.sequences[index]

    def write(self, filename: str, line_length=0) -> None:
        write_string = ""
        for sequence in self.sequences:
            header = sequence.header
            sequence_string = sequence.sequence
            if line_length > 0:
                split_sequence = "\n".join(sequence[i:i + line_length] for i in range(0, len(sequence_string), line_length))
            else:
                split_sequence = sequence_string
            write_string += f"{header}\n{split_sequence}\n"
        io_helpers.write_file(f"{filename}.fasta", write_string)
