from typing import List, Tuple, Optional

def read_lines(
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

def read_file(
    filepath: str
) -> Optional[List[str]]:
    """
    Read contents from a file

    Args:
        filepath (str): The path to the file from which should be read

    Returns:
        str or None: The file contents as string or None if there was an error reading from the file.
    """
    try:
        with open(filepath, 'r') as file:
            content = file.read()
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

def write_lines(
    filepath: str,
    contents: List[str]
) -> None:
    """
    Write contents to a file
    
    Args:
        filepath (str): The path of the file to be written
        contents (List[str]): The lines to be written as list.
    """
    try:
        with open(filepath, 'w') as file:
            file.writelines(contents)
    except PermissionError:
        print(f"You do not have permission to write the file at {filepath}.")
    except FileExistsError:
        print(f"There is already a file at {filepath}. Aborting.")
    except IOError:
        print("An I/O error occured.")
    except Exception as e:
        print(f"An unspecified error has occured while opening the file at {filepath}: \n {e}")

if __name__ == "__main__":
    None
