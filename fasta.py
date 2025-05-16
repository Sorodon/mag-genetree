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

class Fasta:
    sequences: List[Sequence]

    def __init__(self, sequences: List[Sequence] = None):
        if sequences is None:
            self.sequences = []
        else:
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

    def __len__(self):
        return len(self.sequences)

    def add(self, sequence: Sequence) -> None:
        self.sequences.append(sequence)

    def delete(self, index: int) -> None:
        del self.sequences[index]

    def get(self, index:int) -> Sequence:
        return self.sequences[index]

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
        io_helpers.write_file(f"{filename}.fasta", content=write_string)


    def read(self, filename: str) -> None:
        file = io_helpers.read_file(filename, lines=True)
        current_header = None
        current_sequence = ""
        for line in file:
            if line.startswith(">"):
                if current_header is not None:
                    sequence = Sequence(header=current_header, sequence=current_sequence)
                    self.sequences.append(sequence)
                current_header = line.strip()
                current_sequence = ""
            else:
                current_sequence += line.strip()
        if current_header is not None:
            sequence = Sequence(header=current_header, sequence=current_sequence)
            self.sequences.append(sequence)
