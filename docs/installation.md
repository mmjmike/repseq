# Installation

## Library installation

* Clone the library into a designated directory: `git clone https://github.com/mmjmike/repseq`
* To update the library use `git pull -p` inside its folder

<br>For instance, if `~/soft` directory is used:
    <br>`mkdir ~/soft`
    <br>`cd ~/soft`
    <br>`git clone https://github.com/mmjmike/repseq`
    <br>`cd ~/soft/repseq`
    <br>`git pull -p`

## Setting up the environment

To resolve all dependencies, you need to set up a Conda environment.

* Download the file main_repseq.yml, for example, into the folder `~/resources/conda_envs`
* Enter that folder: `cd ~/resources/conda_envs`
* Make sure Conda is working: You should see (base) or another environment name on the left side of the command prompt. If it’s not active, run: `conda activate bash`
* Create the environment and install dependencies from the .yml file: `conda env create -n main_repseq -f main_repseq.yml` (-n is used to specify the environment name)
* Wait for all dependencies to finish installing. Shortly after that, the new environment will appear in the environment list in Jupyter Hub.

In Jupyter Hub, make sure to select the installed environment from the menu in the upper-right corner.
For clustering and calculating statistics, it's recommended to run the Jupyter Hub server in a short session with 32 processors to fully utilize parallel computations during clustering.