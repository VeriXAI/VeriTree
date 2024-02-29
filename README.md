# Safety Verification of Decision-Tree Policies in Continuous Time

This is the repeatability package for the NeurIPS 2023 paper of the same title.
To cite the work, please use this citation:

```
@inproceedings{SchillingLDL23,
  author       = {Christian Schilling and
                  Anna Lukina and
                  Emir Demirovi{\'{c}} and
                  Kim Guldstrand Larsen},
  title        = {Safety verification of decision-tree policies in continuous time},
  booktitle    = {{NeurIPS}},
  volume       = {36},
  pages        = {14750--14769},
  publisher    = {Curran Associates, Inc.},
  year         = {2023},
  url          = {https://proceedings.neurips.cc//paper_files/paper/2023/hash/2f89a23a19d1617e7fb16d4f7a049ce2-Abstract-Conference.html}
}
```


# Installation

The following has been tested on an Ubuntu Linux computer, but it should work
across several major platforms.

Install a Julia compiler for your system.
You can download one [here](https://julialang.org/).
The instructions below assume you have added the `julia` command to your `PATH`
variable.


# Running

Navigate to this folder in the terminal.
Run Julia with the following command.

```bash
julia --project=. experiments/neurips2023_experiments.jl
```

The script tells you that it will write all plots to the folder `plots`.

The script also tells you that verification runs are executed twice.
This is because Julia is just-in-time compiled, so the first run compiles
everything; this is a common standard for just-in-time-compiled languages.
Note that the scripts only run the benchmarks once.

Plotting of the quadrotor results takes a long time (~15 minutes).
