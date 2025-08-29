# Pipeline for gene tree calculation from annotated long read MAGs and their visualization
A pipeline to turn annotated Bakta fasta-files into trees in newick format with options to output alignments and trees as images or ascii-render.

## Usage
`python src/main.py <Options> DATA`

Positional arguments:
- `DATA`: Path to a csv-file containing the data files and identifiers [file, name]

| Option            | Shorthand | Parameter    | Default             | Description                                                             |
| ----------------- | --------- | ------------ | ------------------- | ----------------------------------------------------------------------- |
| `--help`          | `-h`      |              |                     | show this help message and exit                                         |
| `--diamond_path`  |           | `FILE`       | `./diamond/diamond` | Path to the diamond executable                                          |
| `--out`           | `-o`      | `FILE`       | stdout              | The output file to save the trees to                                    |
| `--alignment_out` | `-a`      | `FOLDER`     | Not saved           | The path where to save the alignment for each cluster as separate fasta |
| `--images`        |           | `FOLDER`     | Not saved           | The path where to save the image for each tree                          |
| `--ascii`         |           | `FOLDER`     | Not saved           | The path where to save the ascii render for each tree                   |
| `--threshold`     | `-t`      | `THRESHOLD`  | `90`                | The minimal similarity between protein sequences to be clustered        |
| `--linclust`      |           |              |                     | Use linclust instead of cluster mode for DIAMOND                        |
| `--verbose`       | `-v`      |              |                     | Set to show more detailed output                                        |
| `--timing`        |           |              |                     | Whether to show detailed timing information                             |
| `--size`          |           | `SIZE`       | `3`                 | The minimal cluster size to keep                                        |
| `--gaps`          |           | `GAP_RATE`   | `1`                 | The max gap rate to keep                                                |
| `--length`        |           | `DIFFERENCE` | `1`                 | The max difference in sequence length (relative) to keep                |
| `--lookup`        |           | `FILE`       |                     | A csv file containing paths to all necessary bakta files in column one  |
| `--ur50`          |           | `THRESHOLD`  | $\infty$            | The max amount of distinct UniRef50 IDs per cluster                     |
| `--ur90`          |           | `THRESHOLD`  | $\infty$            | The max amount of distinct UniRef90 IDs per cluster                     |
| `--ur100`         |           | `THRESHOLD`  | $\infty$            | The max amount of distinct UniRef100 IDs per cluster                    |
| `--nopurge`       |           |              |                     | Do not purge singleton clusters before parsing them                     |

