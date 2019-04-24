# Protfarm

Protfarm is a Python package designed for managing large amounts of parallel
sequencing protein engineering datasets. It relies heavily on the functionality
of the [pepars](https://github.com/GradinaruLab/pepars) package, and mostly just
provides a database and API for storing, querying, and analyzing multiple
protein engineering experiments, with multiple samples per experiment. 

## Prerequisites

Protfarm relies on the pepars package, so this must be installed before
installing protfarm. Refer to the
[pepars repo](https://github.com/GradinaruLab/pepars) for installation
instructions.

## Installation

You can then install protfarm via pip:
```
pip install git+https://github.com/GradinaruLab/protfarm.git
```
Note: use pip or pip3, whichever is associated with your Python3

## GUI

If you prefer to work with a GUI instead of in Jupyter notebooks, Protfarm also
has a graphical user interface, which is available as a separate package, here:
[https://github.com/GradinaruLab/protfarm-gui](
https://github.com/GradinaruLab/protfarm-gui)

## Terminology

Before diving into the functionality, it is important to establish some
terminology that will be used throughout the instructions, and the package.

- **Data Path**: Protfarm expects all the data for all experiments to be
contained within a single folder. This folder will largely be managed by
Protfarm, so requires little user interaction (other than dropping raw FASTQ
files in, or getting exported data out). It is recommended to create a folder
just for this purpose, and not manually put any data/analysis here. **There are
no guarantees that Protfarm won't delete/manipulate data in this folder!**
- **Experiment**: Protfarm organizes data into "Experiments". An experiment can
consist of any number of samples and rounds of data, but every sample in an
experiment is expected to have the same mutation strategy. For example, if you
mutate two separate regions of a protein as part of two different libraries,
all the data associated with these two regions should be bundled together into
two separate experiments.
- **Sample/Library**: The terms sample and library are used interchangeably
(sorry), and refer to a set of FASTQ files that should be bundled together for
determining variant counts. Typically, a sample/library has a one-to-one
correlation with a particular index in a sequencing run.

## Functionality

The easiest way to see the functionality of Protfarm is to follow the examples
in the following order:
1. ```workspace/workspace_initialization.ipynb```: Set up a new workspace from
scratch
2. ```alignment/alignment.ipynb```: Do an alignment of some FASTQ files against
a template
3. ```alignment/export_alignment.ipnyb```: Export the variant counts to a CSV
file
4. ```analysis/sequence_counts.ipynb```: Investigate the sequence counts to see
library diversity
5. ```analysis/export_enrichment.ipynb```: Calculate the enrichment of one
sample over another, and export it
6. ```analysis/amino_acid_heatmap.ipynb```: Look at the distribution of amino
acids, normalized by their intrinsic bias
7. ```analysis/collapsing_sequences.ipynb```: Collapse sequences that are
likely to be PCR errors, and compare the library diversity

## License

License information can be found in the LICENSE file