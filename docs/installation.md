# Installation

!!! note "PyPI"
    Currently, a more straightforward installation using PyPI (The Python Package Index) isn't available, however, this feature is planned for future releases. 

## Library installation

Clone the library into a designated directory: `git clone https://github.com/mmjmike/repseq`

## Library update

To update the library, use `git pull -p` inside its folder

<br>For instance, if `~/soft` directory is used:
    <br>`mkdir ~/soft`
    <br>`cd ~/soft`
    <br>`git clone https://github.com/mmjmike/repseq`
    <br>`cd ~/soft/repseq`
    <br>`git pull -p`

## Setting up the environment

To resolve all dependencies, you need to set up a Conda environment.

* The .yml environemnt file (`main_repseq.yml`) should be inside 'repseq' folder after the library installation or an update
* Enter this folder. For instance, if `~/soft` directory is used: `cd ~/soft/repseq`
* Make sure Conda is working: You should see (base) or another environment name on the left side of the command prompt. If it’s not active, run: `conda activate bash`
* Create the environment and install dependencies from the .yml file: `conda env create -n main_repseq -f main_repseq.yml` (-n is used to specify the environment name)
* Wait for all dependencies to finish installing. Shortly after that, the new environment will appear in the environment list in Jupyter Hub.

In Jupyter Hub, make sure to select the installed environment from the menu in the upper-right corner.
For clustering and calculating statistics, it's recommended to run the Jupyter Hub server in a short session with 32 processors to fully utilize parallel computations during clustering.