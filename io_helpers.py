from typing import List, Tuple, Optional

def read_file(
    filepath: str
) -> Optional[List[str]]:
    """
    Read contents from a file

    Args:
        filepath (str): The path to the file from which should be read

    Returns:
        List[str] or None: A list of lines in the file or None if there was an error reading from the file
    """
    try:
        with open(filepath, 'r') as file:
            content = file.readlines()
            return content
    except PermissionError:
        print(f"You do not have permission to write the file at {filepath}.")
    except IOError:
        print("An I/O error occured.")
    except Exception as e:
        print(f"An unspecified error has occured while opening the file at {filepath}: \n {e}")

def write_file(
    filepath: str,
    contents: str
) -> None:
    """
    Write contents to a file
    
    Args:
        filepath (str): The path of the file to be written
        contents (str): The contents to be written.
    """
    try:
        with open(filepath, 'w') as file:
            file.write(contents)
    except PermissionError:
        print(f"You do not have permission to write the file at {filepath}.")
    except FileExistsError:
        print(f"There is already a file at {filepath}. Aborting.")
    except IOError:
        print("An I/O error occured.")
    except Exception as e:
        print(f"An unspecified error has occured while opening the file at {filepath}: \n {e}")

def parse_fasta(
        filepath: str
) -> List[Tuple[str,str]]:
    """
    Get a list of seqences from a fasta file
    
    Args:
        filepath (str): A filepath pointing to the fasta file

    Returns:
        List[Tuple[str,str]]: A list of Tuples each containing a header and sequence
    """
    content = read_file(filepath)
    sequences = []
    if not content:
        return sequences

    current_header = None
    current_sequence = []
    for line in content:
        if line.startswith(">"):
            if current_header is not None:
                sequences.append((current_header, "".join(current_sequence)))
            current_header = line
            current_group = []
        else:
            current_group.append(line)
    return sequences

def write_fasta(
        filepath: str,
        content: List[Tuple[str, str]]
) -> None:
    write_string = ""
    for element in content:
        header = element[0]
        sequence = element[1]
        split_sequence = "\n".join(sequence[i:i + line_length] for i in range(0, len(sequence), line_length))
        write_string += f"{header}\n{split_sequence}\n"
    write_file(filepath, write_string)

if __name__ == "__main__":
    None
