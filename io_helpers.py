# vim: set foldmethod=marker:
# vim: set foldclose=all foldlevel=0:
# vim: set foldenable: 

from typing import List, Tuple, Optional, Union
import re

def read_file( #{{{
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
                content = [line.strip() for line in file.readlines()]
            else:
                content = file.read().strip()
            return content
    except PermissionError:
        print(f"You do not have permission to write the file at {filepath}.")
    except IOError:
        print("An I/O error occured.")
    except Exception as e:
        print(f"An unspecified error has occured while opening the file at {filepath}: \n {e}")
#}}}

def write_file( #{{{
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
#}}}

def parse_csv( #{{{
        filepath: str,
        sep: str = ",",
        header_row: bool = False,
        skip: int = 0
) -> List[List[str]]:
    """
    Parse a csv (or similar file) and return a list of lists (containing columns)

    Args:
        filepath (str): The path to the csv file
        sep (str, optional): The separator for columns. (defaults to ',')
        header_row (bool): If True, will output a list of dictionaries with named column fields using the first row as names. (defaults to False)
        skip (int): Skips the specified amount of rows from the top of the file. (defaults to 0)

    Returns:
        List[List[str, ...]]: A list of lists, one list per line with as many strings as defined in columns.
    """
    contents = read_file(filepath=filepath, lines=True)
    contents = contents[skip:]
    result = []
    for line in contents:
        columns = []
        while sep in line:
            column, _, line = line.partition(sep)
            columns.append(column)
        columns.append(line)
        result.append(columns)
    if header_row:
        dicts = []
        header = result.pop(0)
        print(header)
        for line in result:
            resulting_dict = {}
            for index, field in enumerate(header):
                try:
                    resulting_dict[field] = line[index]
                except IndexError:
                    None
            dicts.append(resulting_dict)
        return dicts
    return result
#}}}

if __name__ == "__main__":
    file = parse_csv(
        "../data/bin.2/bin.2.tsv",
        sep = "\t",
        skip = 5,
        header_row = True
    )
    [print(line) for line in file]
