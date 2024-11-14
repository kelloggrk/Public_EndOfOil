# Code and replication package for "The End Of Oil" by Ryan Kellogg


## Overview
---------------------------------------
There are two components to this code package. The first is a script aimed at users who want to set their own input parameters and run their own specifications of the paper's model. The second component is a set of scripts that replicates all results, across all specifications, in the paper and appendix.

For both components, all code uses the Julia programming language. If you are unfamiliar with Julia, an excellent overview and set of tutorials aimed at economists is provided in [Quantitative Economics with Julia](https://julia.quantecon.org/intro.html), by Jesse Perla, Thomas J. Sargent, and John Stachurski.

To either run a custom version of the model or replicate the entire paper, you should first clone this repository to your machine, e.g. to a directory such as `/Users/kelloggr/EndOfOil` or `C:/Work/EndOfOil`. 

The rest of this README file is composed of three sections:
1. Software and package requirements
2. How to use the provided script `CustomEndOfOilModel.jl` to run your own version of the model
3. How to replicate all results and the paper.


## Software and package requirements
---------------------------------------
To install Julia, see [this link](https://docs.julialang.org/en/v1/manual/installation/) for installation instructions for all common operating systems. I strongly recommend using Julia in conjunction with [Visual Studio Code](https://code.visualstudio.com/) as an integrated development environment.

Once Julia is installed, you must install the necessary packages via the commands below (these commands are also included as commented lines in `code/analysis/EOOmodel.jl`):

import Pkg; Pkg.add("Parameters")

import Pkg; Pkg.add("Roots")

import Pkg; Pkg.add("QuadGK")

import Pkg; Pkg.add("NLsolve")

import Pkg; Pkg.add("BenchmarkTools")

import Pkg; Pkg.add("JLD2")

import Pkg; Pkg.add("FileIO")

import Pkg; Pkg.add("Printf")

import Pkg; Pkg.add("Plots")


## Running a custom specification of the model
---------------------------------------

### File
To run a single, custom specification of the model, open `code/analysis/CustomEndOfOilModel.jl`. This file will let you set any parameter input --- and potentially change the number of producing regions if desired --- and then simulate a full set of results, along with plots corresponding to figures 2 and 3 in the paper.

### User directory
The first block of code in `CustomEndOfOilModel.jl` sets the working directory to `EndOfOil/code/analysis` and then loads the model by reading in the scripts `EOOmodel.jl`, `EOOmodel_AKS.jl`, and `PlotFunctions.jl` from this directory. 

If the Julia environment is loaded by opening `CustomEndOfOilModel.jl` to launch Visual Studio Code (e.g., in Windows, by double clicking on `CustomEndOfOilModel.jl` in a Windows folder explorer), then running this first block of code will automatically set the working directory to the path to `EndOfOil/code/analysis` on your machine. If not (for instance, if you load a Julia environment and then open `CustomEndOfOilModel.jl`), you will need to manually set the working directory, `dircurrent`, to your path.

### Input parameters
The second block of code is where users can enter input parameters. The defaults match the reference case model from the paper. If you do not change any of this code, then `CustomEndOfOilModel.jl` will reproduce the reference case results.

#### Flags
The input parameter code block begins with three flags that affect how the inputs are used:
1. Does the model include investment? Setting `p_AKS = true` includes investment dynamics in the model. Set this flag to `false` if you desire to run an extraction model without any investment dynamics. Note that in this case, the lambda parameters will be ignored.
2. Should the cost slope parameters gamma be re-calibrated before running the model? Setting `p_estslopes = true` will run code that calibrates the gamma parameters so that, given the other parameter inputs, simulated baseline 2023 investment and production match actual 2023 investment and production.
    * Calibrating the slope parameters can require 30-60 minutes of run time, depending on the model and the quality of the initial guess. The default is `p_estslopes = false` so that users can quickly run the reference case model if desired without modifying any code.
3. Should the model run the unanticipated demand decline scenario? Setting `p_unant = true` (the default) will run the unanticipated decline scenario. Because this scenario requires re-solving the model in every time step, it can require ~30 minutes to run
    * Running the model with both `p_estslopes = false` and `p_unant = false` is fast (less than one minute typically). But then you will not be able to compare outcomes under the anticipated vs unanticipated decline.

#### Parameter inputs
After the three flags, you may modify any of the parameter inputs. Some notes:
* It is possible to change the number of regions by changing the lengths of the vectors of region-specific inputs. All vectors must have the same length.
* The drilling cost slope parameters `p_gamma` are used directly if the flag `p_estslopes = false`. If `p_estslopes = true`, then these values are used as the initial guess in the algorithm that searches for the slope parameters that equate simulated investment to actual investment.
  - If you change any of the other input parameters and set `p_estslopes = false`, then the model will generically not reproduce 2023 actual investment at baseline
  - To help users set good initial guesses of the gamma parameters, I have provided further below the values of gamma that correspond to every calibrated model specification in the paper.
* `mu_0_bbl_guess` are guesses of the initial values of baseline mu0 (in $/bbl). These are used in the calibration algorithm.
  - To help set good initial guesses of the mu0, I have provided further below the values of baseline mu0 that correspond to every calibrated model specification in the paper.
  - If `p_estslopes = true` and the initial guesses of gamma and mu0 are poor, the calibration routine may fail to converge.


### Values for the drilling marginal cost slopes gamma and the initial baseline shadow values mu0
The list below provides the values of gamma and baseline mu0 for every specification used in the paper. For any of these specifications, if the gamma values are used as inputs for `p_gamma` in the `CustomEndOfOilModel.jl` script, and if all other inputs correspond to that specification, then the model that is created will replicate the paper's results. Using the values of mu0 for `mu0_bbl_guess` will speed computation if that particular specification is used.

For custom specifications, this list of values can serve as a guide to values of `p_gamma` and `mu0_bbl_guess` that are good initial guesses for setting a baseline version of the model that will quickly converge during calibration.

The units for mu0 are always in dollars/bbl. The units for gamma are in dollars/bbl per mmbbl/d of investment for the models that include investment dynamics. The units are dollars/bbl per mmbbl/d of production for models without investment dynamics.
* Note that as lambda increases, the values of gamma decrease because the amount of first period investment, in mmbbl/d of capacity created, required to maintain baseline production growth increases mechanically with lambda.

Reference case:
* p_gamma = [30.3947, 53.8906, 50.0906, 36.2644]
* mu0_bbl_guess = [30.8797, 3.59828, 2.7732, 2.37266]

All resources have an 8% production decline:
* p_gamma = [30.3435, 53.8401, 50.0472, 136.798]
* mu0_bbl_guess = [30.8941, 3.64174, 2.80363, 4.54921]

All resources have an 6% production decline:
* p_gamma = [36.9182, 69.9794, 65.3401, 175.279]
* mu0_bbl_guess = [34.5792, 5.01951, 3.78779, 6.34039]

All resources have an 30% production decline:
* p_gamma = [9.37405, 14.1581, 13.0799, 36.5041]
* mu0_bbl_guess = [23.5702, 1.79032, 1.4393, 2.20034]

No investment
* p_gamma = [1.65876, 2.33872, 2.1565, 6.03627]
* mu0_bbl_guess = [21.0618, 1.35088, 1.09661, 1.66459]

High reserves
* p_gamma = [56.9581, 55.5353, 51.2496, 37.6623]
* mu0_bbl_guess = [9.4789, 0.296036, 0.203891, 0.225045]

Low reserves
* p_gamma = [16.2907, 51.5204, 46.0072, 32.7213]
* mu0_bbl_guess = [44.0237, 8.62873, 9.70362, 7.64653]

All resources discount at 9%
* p_gamma = [68.7178, 53.7133, 49.8421, 36.3766]
* mu0_bbl_guess = [1.56015, 3.28459, 2.53182, 2.09518]

All resources discount at 3%
* p_gamma = [33.6981, 31.5315, 28.5105, 18.4386]
* mu0_bbl_guess = [30.2892, 37.3502, 33.6956, 27.8425]

No investment, low reserves
* p_gamma = [1.081, 2.27388, 2.06612, 5.59849]
* mu0_bbl_guess = [29.9883, 3.32314, 3.66968, 5.3507]

No investment, all resources discount at 3%
* p_gamma = [1.65874, 1.55372, 1.3852, 3.33266]
* mu0_bbl_guess = [21.0621, 25.2289, 23.0549, 24.429]

Russia in core OPEC, 9% discounting
* p_gamma = [21.3528, 83.1755, 49.917, 36.4503]
* mu0_bbl_guess = [2.263, 2.50552, 2.45056, 2.00456]

No market power
* p_gamma = [73.7888, 53.7988, 49.979, 36.2637]
* mu0_bbl_guess = [31.642, 3.55639, 2.74792, 2.32999]

High demand elasticity
* p_gamma = [42.458, 53.9129, 50.0939, 36.3957]
* mu0_bbl_guess = [26.8851, 3.22843, 2.42841, 2.10341]

Low demand elasticity
* p_gamma = [13.4203, 53.9187, 50.1669, 36.072]
* mu0_bbl_guess = [36.7308, 4.18742, 3.31684, 2.81209]

No demand growth
* p_gamma = [34.4185, 56.8003, 52.7156, 36.7353]
* mu0_bbl_guess = [29.2022, 3.13555, 2.40437, 2.17779]

Higher demand growth
* p_gamma = [22.6495, 47.9863, 44.7794, 35.1561]
* mu0_bbl_guess = [34.8461, 4.8601, 3.78892, 2.9178]

High cost function intercepts
* p_gamma = [11.3146, 42.7438, 38.0886, 25.6751]
* mu0_bbl_guess = [29.4419, 3.33609, 2.62382, 2.17272]


### Running the model
Once all inputs are specified, the code in the block below the heading `RUN MODEL. NO NEED FOR USERS TO EDIT THIS CODE` will run the calibration (if `p_estslopes = true`) and then run the model. The comments below this block of code provide a dictionary for the output variables and figures.




## Replicating the full paper
---------------------------

### Initial downloads:
Users interested in replicating the paper should download to their machines [this publicly-accessible (zipped) folder](https://www.dropbox.com/scl/fo/gjg0pgc4077i4278iag3p/AD7QCuqWESxDTOYNpFYY1Gk?rlkey=cpzu5a4ghrl001qzug0qmuz5z&st=nf54jh1c&dl=0) (in addition to cloning this repository to their machines as discussed at the top of this README). I recommend downloading this folder to a directory such as `/Users/kelloggr/EndOfOil_data` or `C:/Users/kelloggr/Dropbox/EndOfOil_data`. This folder contains two sub-folders:
* `intermediatefiles` holds intermediate data files (.jld2 files that contain calibrated parameter sets and simulation results) that facilitate the paper's analysis.
* `CartoonFigures` includes two figures that are included in the paper's theory section that are not produced by the Julia scripts.


### Directory specification:
* So that the Julia scripts locate your local copy of `intermediatefiles`, you must edit `code/analysis/intfilepath.jl`. This file should contain a single line that contains your full path to this directory. For instance, my own path is:
```
dirint = "C:/Users/kelloggr/Dropbox/EndOfOil/intermediatefiles/"
```
  - If you are working collaboratively via GitHub, I recommend adding `intfilepath.jl` to `.gitignore` so that you and your collaborators can have different paths.
* So that the shell scripts locate your local files, you must specify the `REPODIR` and `DBDIR` variables in `EndOfOil_bash_file.sh`. These point to your local root repo directory and the directory holding `intermediatefiles` (on my own machine, this folder lives on dropbox)
  - `REPODIR` should look something like `REPODIR=C:/Work/EndOfOil`
  - `DBDIR` should look something like `DBDIR="C:/Users/kelloggr/Dropbox/EndOfOil_data`
    * Note: do NOT include a white space on either side of the equal sign in any of these expressions
* `EndOfOil_bash_file.sh` also specifies the number of threads to use when running the simulation scripts in parallel. The default is `export JULIA_NUM_THREADS=4`. You can set a number that is appropriate for your machine.


### Running the scripts
* The `EndOfOil_bash_file.sh` shell script is located in the root repo directory. Executing this script will: (1) delete all intermediate data and results; (2) copy the files in `CartoonFigures` into `REPODIR/output/figures`; (3) conduct all of the calibration and simulations; (4) output all figures, tables, and single number tex files into subfolders of `REPODIR/output`; and (5) compile the paper and appendix in `REPODIR/paper`.
  - The deletion of the intermediate data and results files ensures that the paper's results are fully replicated from scratch and that there are no hidden, improper file dependencies.
  - Users who only wish to run part of the code (e.g. the simulations but not the calibration) should NOT run the script that deletes the intermediate data and results files.
* To execute `EndOfOil_bash_file.sh`, I recommend opening your bash shell, changing the directory to your local repository, and then using the following command
```
bash -x EndOfOil_bash_file.sh |& tee EndOfOil_bash_file_out.txt
```
  - This command will log output and any error messages to `EndOfOil_bash_file_out.txt`
* The entire script takes roughly one day to execute on my desktop workstation.
