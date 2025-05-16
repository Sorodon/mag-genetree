from typing import List, Tuple, Optional, Union
import re

def read_file(
    filepath: str,
    lines: bool = False
) -> Optional[Union[str, List[str]]]:
    """
    Read contents from a file

    Args:
        filepath (str): The path to the file from which should be read
        lines (bool, optional): Whether to return a list of strings instead of one big string.

    Returns:
        str or List[str] or None: A string, a list of lines in the file or None if there was an error reading from the file
    """
    try:
        with open(filepath, 'r') as file:
            if lines:
                content = file.readlines()
            else:
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
    content: Union[str, List[str]]
) -> None:
    """
    Write content to a file
    
    Args:
        filepath (str): The path of the file to be written
        content (str, List[str]): The content to be written. Either a sing€ string or a list of strings.
    """
    try:
        with open(filepath, 'w') as file:
            if isinstance(content, list):
                file.writelines(content)
            else:
                file.write(content)
    except PermissionError:
        print(f"You do not have permission to write the file at {filepath}.")
    except FileExistsError:
        print(f"There is already a file at {filepath}. Aborting.")
    except IOError:
        print("An I/O error occured.")
    except Exception as e:
        print(f"An unspecified error has occured while opening the file at {filepath}: \n {e}")

def parse_csv(
        filepath: str,
        columns: int,
        sep: str = ",",
        pattern = r"([a-zA-z0-9_]*)"
) -> List[Tuple[str, ...]]:
    """
    Parse a csv (or similar file) and return a list of tuples (containing columns)

    Args:
        filepath (str): The path to the csv file
        columns (str): The number of columns to expect
        sep (str, optional): The separator for columns (defaults to ',')
        pattern (str, optional): The pattern matching single fields. Should contain one capture group. Defaults to r'([a-zA-z0-9_])'

    Returns:
        List[Tuple[str, ...]]: A list of tuples, one tuple per line with as many strings as defined in columns. Lines not matching the pattern get ignored.
    """
    contents = read_file(filepath=filepath, lines=True)
    full_pattern = f"^{sep.join([pattern for _ in range(columns)])}$"
    result = []
    for line in contents:
        match_obj = re.fullmatch(full_pattern, line.strip())
        if match_obj:
            result.append(tuple(match_obj.group(i) for i in range(1,columns+1)))
    return result

if __name__ == "__main__":
    None
