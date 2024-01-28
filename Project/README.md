
# How to setup PyPulseq with KomaMRI with Pulse Sequence example

## Installation
For setting up the PyPulseq installation, we use a Conda environment with python 3.10.12 for compatibility. As of January 2024, the latest version of PyPulseq is compatible with this version.

Run the following commands where Conda is installed
```
conda create --name pypulseq -python=3.10.12
conda activate pypulseq
pip install pypulseq
```
Now that we have installed PyPulseq, let's continue with KomaMRI. We are following the guide at this link: [https://cncastillo.github.io/KomaMRI.jl/stable/getting-started/](https://cncastillo.github.io/KomaMRI.jl/stable/getting-started/). First, we download the appropriate Julia version, as of now it is 1.8.0 and can be found here: [https://julialang.org/downloads/](https://julialang.org/downloads/).

If you are on a unix-based system extract and move the installation to appropriate destination and add to path:
```
tar -xvzf julia-1.8.0-linux-x86_64.tar.gz
mv julia-1.8.0 ~
export PATH="~/julia-1.8.0/bin:$PATH"

