# AA table visualization

## Installation:

1. You will need a copy of conda in order to manage dependencies. We strongly
recommend using 'mamba' derivatives of conda, as these provide built-in strict
version control when downloading packages, have more robust package solving, and
provide quicker package download times. The easiest way to obtain mamba is with
micromamba, available here:
https://mamba.readthedocs.io/en/latest/installation/micromamba-installation.html
The relevant command is here:
```bash
"${SHELL}" <(curl -L micro.mamba.pm/install.sh)
```
2. You can obtain a copy of this git repository with this command:
```bash
git clone https://github.com/simkin-bioinformatics/AA_table_visualization.git
```
3. After changing directory to the cloned AA_table_visualization folder,
you can cd into the folder (if you're unsure whether you're in the correct
folder, make sure the folder you cd into contains an environment.yaml file) and
run this command to build a conda environment that contains all package
dependencies (you can substitute micromamba or conda instead of mamba depending
on which type of conda you have):
```bash
mamba env create -f environment.yaml
```

4. (optional - useful if your system times out or throws errors during the
static image graphing step). Install chrome for plotly by activating your
environment and running the included install script:
```bash
mamba activate aa_table_visualization
```
```bash
plotly_get_chrome
```

## Usage:

You'll need to have coverage, alternate, and reference AA tables copied into the
github folder that you cloned, along with a metadata table that includes columns
with sample names, latitudes, longitudes, and location names that you would like
to use to group your samples together (e.g. geographic region, district, or
health facility names)

1. activate the environment with:
```bash
mamba activate aa_table_visualization
```

2. Launch the jupyter notebook with:
```
jupyter lab variant_graphing.ipynb
```

3. make sure to open the jupyter notebook in a web browser (e.g. using the link
included in the output messages), and open the variant_graphing.ipynb file (e.g.
by double clicking it). Follow the instructions in the notebook carefully. An
example usage visualizing a dataset from Tanzania is pre-filled out in the
notebook.
