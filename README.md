# pyrk
Analysis framework for R(J/Psi) analysis (forked from R(K) analysis)

## Installation

To ease dependency tracking this package is meant to run in a conda environment (python 3)
First install miniconda following the instructions [here](https://docs.conda.io/en/latest/miniconda.html).

Then clone this package and install the environment:
```
git clone --recursive git@github.com:friti/pyrk.git
cd pyrk
conda create --name pyrk python=3.7
conda activate pyrk
conda install six
pip install -r env.pip
conda install -c conda-forge root
```

You are ready to go! Remember to activate the pyrk environment at each login!
