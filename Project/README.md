
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
```
```
julia> ]
(@v1.8) pkg> add KomaMRI
julia> using KomaMRI
```

To start the KomaUI with defalult settings run:
`
julia> KomaUI()
`

If you want to customize the UI settings, `launch_komaui.jl` shows an example. 

The latest Julia version was not compatible with KomaMRI but 1.8 worked well. 

For Pulseq, just download the github repository and add it to your path in Matlab

`git clone git@github.com:pulseq/pulseq.git`


## To test out the pipeline 
Run `spiral.m` to generate an example spiral sequence from Pulseq. Simulate the sequence with `simulation.jl` and do a re-gridding and reconstruction with `gridding.m`


## Project files
The implemented four-shot sequence is in `four_shot.m` and the corresponding gridding which is currently not working is found in `gridding_four_shot.m`.
