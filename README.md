
## Cloud-Based Tide Prediction Software

### What is this software?

This software is used to predict the sea surface-height anomaly associated with the
baroclinic (internal) tide. The user inputs a filename containing a list of
longtiude, latitude, time coordinates, and the software will output the predicted
tidal elevation from the HRET 8.1 baroclinic tide model. For detailed usage
instructions, see below.

### Who is this software for?

This software is useful for researchers who would like to predict baroclinic sea
level anomaly or baroclinic surface pressure associated with the main ocean
tides, M2, S2, K1, and O1, as well as the annual modulates of M2, denoted
MA2 and MB2. The predictions are based on a model for the phase-locked tides
inferred from the 30-year record of satellite altimetry. Thus, the predictions
will not include modulations related to refraction or other interactions
with mesoscale processes, the so-called non-phase-locked, nonstationary, or
incoherent tide. The amplitude of the non-phase-locked tide is thought to
exceed the phase-locked tide throughout the tropics, so the predictions
offered by this program may differ by several centimeters from the true baroclinic
sea level anomaly at a particular instant.

### What is the status of this software?

This software is under development. The prediction of tidal sea surface-height anomaly
is routine, but I expect to add new capabilities for predicting subsurface
properties (i.e., tidal currents and density anomalies). The main limitation at the moment
is the provisioning of enough space in the cloud for the input data files. Roughly 10GB are
required to store the vertical modes and other information needed to convert HRET 8.1 into
subsurface predictions.

### How to get more information?

Email the author: [Edward D. Zaron](mailto:edward.d.zaron@oregonstate.edu)

Read about HRET 8.1: [Baroclinic Tidal Sea Level from Exact-Repeat Altimetry](http://dx.doi.org/10.1175/JPO-D-18-0127.1)

Learn about tides and tidal analysis: [Coastal Tides](http://refmar.shom.fr/sea_level_news_2013/2013_t4/ouvage-reference-sur-theorie-et-pratique-maree-francais-et-anglais) by Bernard Simon.

## Quickstart

### Installation

- Login to the SMCE cloud environment (Jupyterlab).

- Open a Terminal and checkout this software:
```
git clone https://github.com/ezaron/SMCE.git SMCE
```

- Open a Terminal and install the Julia language:
```
cd
mkdir -p opt/src
cd opt/src
wget https://julialang-s3.julialang.org/bin/linux/x64/1.8/julia-1.8.0-linux-x86_64.tar.gz
tar -xvzf julia-1.8.0-linux-x86_64.tar.gz
mkdir -p opt/bin
mv julia-1.8.0/bin/julia ../bin
export PATH=/home/jovyan/opt/bin:${PATH}
```

- Get Julia ready for use:
```
cd
cd SMCE
export PATH=/home/jovyan/opt/bin:${PATH}
julia ./setup.jl
```

### Options for using the software:

- Open the .ipynb in this directory, and launch an interactive session.

- Open a Terminal and run the software, `julia ./test.jl <PATHSPEC>`.

- Open a Julia notebook: click on File > New Launcher > Notebook > Julia to launch a Julia session and experiment with the `test.jl` script.

### Command-line usage:

```
julia test.jl <PATHSPEC>
```
1. If _PATHSPEC_ ends in ".nc" (i.e., if it is the full path to a NetCDF file)
   then we compute baroclinic tide SLA predictions and write them in a file 
   called `PATHSPEC_hret.nc`.
   
2. If _PATHSPEC_ ends in ".txt" then we treat the contents of _PATHSPEC_
   as a list of files and generate predictions, as above.
   
3. If _PATHSPEC_ ends in "/" then we compute tide predictions for all
   the files with names ending in ".nc" in this directory, as above.

This program will generate the output files, mentioned above, in addition to simple plots
of the input coordinates (latitude,longitude), useful for basic sanity checking.

### Caveats:

This program attempts to read and parse the input files to extract
the latitude, longitude, and time at which predictions are requested.
It has enough logic to handle cases, such as the simulated SWOT
data files, where the two-dimensional (latitude,longitude) arrays
must be matched with a one-dimensional time array.

The most important thing to verify, as a user, is that the
reference time for the time variable has been parsed correctly.
The start and end times of each input file, in YYYY-MM-DD HH:MM:SS
format, are output to the console and also included on the figures
created by this program.
