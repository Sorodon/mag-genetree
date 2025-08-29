from io import StringIO
from Bio import Phylo
import matplotlib.pyplot as plt
from typing import Optional

def draw(
    newick:str,
    mode:str = 'show',
    path:Optional[str] = None
) -> None:
    """
    Draw a tree provided in Newick format.
    Args:
        newick (str): The tree in Newick notation
        mode (str): Either 'show', 'save' or 'ascii'. 'show' opens the tree directly, while 'save' and 'ascii' save the tree either as image or ascii representation to a file defined in <path>. Defaults to 'show'.
        path (str | None): The path to save the tree when <mode> is set to either 'save' or 'ascii'. Defaults to None.
    Returns:
        None
    """
    tree = Phylo.read(StringIO(newick), 'newick')
    match mode:
        case 'show':
            Phylo.draw(tree)
        case 'save':
            try:
                Phylo.draw(tree, do_show=False)
                plt.savefig(path)
            except Exception as e:
                print(f"There was an error saving the file: {e}")
        case 'ascii':
            try:
                with open(path, 'w') as file:
                    Phylo.draw_ascii(tree, file=file)
            except Exception as e:
                print(f"There was an error saving the file: {e}")
