#!/bin/sh

#-------------------------------------------------------------------------------
# Name:        EndOfOil_bash_file.sh
# Purpose:     Calls every file from calibration through final simulations,
#              and compilation of the paper
#
# Created:     8 November, 2024
#-------------------------------------------------------------------------------


# ---------------------------------------------------------------------------
# Execution
# Users can run this script by opening a bash/unix shell, changing the directory
# to their local repository ($REPODIR below), and then using this command:
#
# bash -x EndOfOil_bash_file.sh |& tee EndOfOil_bash_file_out.txt
#
# ---------------------------------------------------------------------------


# USER INPUT NEEDED HERE FOR ENVIRONMENT VARIABLES:
# USER-SPECIFIC PATHS TO THIS SCRIPT'S DIRECTORY AND THE DATA DIRECTORY
# USER-SPECIFIC OPERATING SYSTEMS
REPODIR=C:/Work/EndOfOil
DBDIR="C:/Users/kelloggr/Dropbox/EndOfOil"


# EXPORT VARIABLES TO SUBSCRIPTS
export REPODIR
export DBDIR


# EXPORT NUMBER OF THREADS TO JULIA
export JULIA_NUM_THREADS=4


# CLEAR INTERMEDIATE RESULTS AND OUTPUT FOLDERS
bash -x $REPODIR/EndOfOil_clear_folders.sh |& tee EndOfOil_clear_folders_out.txt


# RUN THE JULIA CALIBRATION SCRIPTS
julia $REPODIR/code/analysis/calibration.jl
julia $REPODIR/code/analysis/calibrationout.jl

# RUN THE JULIA SIMULATION SCRIPTS
julia $REPODIR/code/analysis/runmodels_altmp.jl
julia $REPODIR/code/analysis/runmodels_altdec.jl
julia $REPODIR/code/analysis/runmodels_n.jl
julia $REPODIR/code/analysis/runmodels_altres.jl
julia $REPODIR/code/analysis/runmodels_altdall1.jl
julia $REPODIR/code/analysis/runmodels_altdall2.jl
julia $REPODIR/code/analysis/runmodels_altdemdec.jl
julia $REPODIR/code/analysis/runmodels_altelast.jl
julia $REPODIR/code/analysis/resultsout.jl


# COMPILE THE PAPER
bash -x $REPODIR/EndOfOil_paper_script.sh |& tee EndOfOil_paper_script_out.txt

# clean up log files
rm *.log

exit

