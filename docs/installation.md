# Installation

## PyPI installation

Install the core package from PyPI:

```bash
pip install repseq
```

Optional functionality is distributed as extras:

```bash
pip install "repseq[clustering]"
pip install "repseq[pgen]"
pip install "repseq[tcrdist]"
pip install "repseq[r]"
```

For development:

```bash
git clone https://github.com/mmjmike/repseq
cd repseq
pip install -e ".[dev]"
```

## Library update

To update an editable development checkout, use `git pull -p` inside its folder.

<br>For instance, if `~/soft` directory is used:
    <br>`mkdir ~/soft`
    <br>`cd ~/soft`
    <br>`git clone https://github.com/mmjmike/repseq`
    <br>`cd ~/soft/repseq`
    <br>`git pull -p`

## Setting up the environment

For conda-based installation, use the conda package once it is available from
the selected channel:

```bash
conda install -c conda-forge repseq
```

For local development with all historical analysis dependencies, you can still
set up the full Conda environment.

* The .yml environemnt file (`main_repseq.yml`) should be inside 'repseq' folder after the library installation or an update
* Enter this folder. For instance, if `~/soft` directory is used: `cd ~/soft/repseq`
* Make sure Conda is working: You should see (base) or another environment name on the left side of the command prompt. If it’s not active, run: `conda activate bash`
* Create the environment and install dependencies from the .yml file: `conda env create -n main_repseq -f main_repseq.yml` (-n is used to specify the environment name)
* Wait for all dependencies to finish installing. Shortly after that, the new environment will appear in the environment list in Jupyter Hub.

In Jupyter Hub, make sure to select the installed environment from the menu in the upper-right corner.
For clustering and calculating statistics, it's recommended to run the Jupyter Hub server in a short session with 32 processors to fully utilize parallel computations during clustering.
