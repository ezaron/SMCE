
## HRET14 Tide Prediction Software

### What is this software?

This software is used to predict the sea surface-height anomaly or ocean surface
current associated with the
baroclinic (internal) tide. The user inputs a filename containing a list of
longtiude, latitude, time coordinates, and the software will output the predicted
tidal elevation and current from the HRET14 baroclinic tide model. Input files in different formats
are supported (RADS, SWOT, and GDP), but the user can modify these for other input
format specifications.
For detailed usage instructions, see below.

### Who is this software for?

This software is for researchers who would like to predict baroclinic sea
level anomaly (baroclinic surface pressure) or baroclinic surface currents
associated with the ocean tides, M2, S2, N2, K1, and O1.
The predictions are based on a model for the phase-locked tides
obtained by mapping historical satellite altimeter data and
hourly velocity estimates from surface drifters. Thus, the predictions
will not include modulations related to time-variable
refraction or other interactions with the environment.

### How to get more information?

Email the author: [Edward D. Zaron](mailto:edward.d.zaron@oregonstate.edu)

Read about the development of HRET14: [manuscript in review at JTech](https://oregonstate.box.com/s/gssl5uqysucqw2f3q7tfkzylpijba2ny)

Read about the previous version of HRET (version 8.1): [Baroclinic Tidal Sea Level from Exact-Repeat Altimetry](http://dx.doi.org/10.1175/JPO-D-18-0127.1)

Learn about tides and tidal analysis: [Coastal Tides](http://refmar.shom.fr/sea_level_news_2013/2013_t4/ouvage-reference-sur-theorie-et-pratique-maree-francais-et-anglais) by Bernard Simon.

## Quickstart for a cloud environment or Linux commandline installation

### Installation

- Login to the SMCE cloud environment (Jupyterlab).

- Checkout the most recent version of this software and data files:
```
fossil clone https://ingria.ceoas.oregonstate.edu/fossil/SMCE
```

- Alternately, if you already have the HRET14 data files, you may grab the source code
from github. Open a Terminal and checkout this software:
```
git clone https://github.com/ezaron/SMCE.git SMCE
```
The HRET14 data files are located at
```
https://ingria.ceoas.oregonstate.edu/fossil/SMCE/dir?ci=tip&name=HRET14
```

- Open a Terminal and install the Julia language:
```
cd
mkdir -p opt/src
cd opt/src
wget https://julialang-s3.julialang.org/bin/linux/x64/1.8/julia-1.8.0-linux-x86_64.tar.gz
tar -xvzf julia-1.8.0-linux-x86_64.tar.gz
mkdir -p opt/bin
ln -s /home/jovyan/opt/src/julia-1.8.0/bin/julia /home/jovyan/opt/bin
export PATH=/home/jovyan/opt/bin:${PATH}
```

- Get Julia ready for use:
```
cd
cd SMCE
export PATH=/home/jovyan/opt/bin:${PATH}
julia ./setup.jl
```

- Move the sample data files into place. The sample data files are located in the './SWOTdata', './RADSdata', and './GDPdata'
directories of this [fossil repository](https://ingria.ceoas.oregonstate.edu/fossil/SMCE/dir?ci=tip)

Within the SMCE environment the sample files can be copied into your user directory as follows:
```
mkdir /home/jovyan/DEMO_FILES
cp /efs/SWOT_shared/data/SWOT_SIMULATED_L2_KARIN_SSH_ECCO_LLC4320_SCIENCE_V1/SWOT_L2_LR_SSH_Expert_018_290_20121112T003212_20121112T012339_DG10_01.nc /home/jovyan/DEMO_FILES/
```

### Options for using the software:

- Open the .ipynb in this directory, and launch an interactive session. This will generate predictions for the files under /home/jovyan/DEMO_FILES/.

- Open a Terminal and run the software, `julia ./test.jl <PATHSPEC>`.

- Open a Julia notebook: click on File > New Launcher > Notebook > Julia to launch a Julia session and experiment with the `test.jl` script.

### Command-line usage:

```
julia test.jl <PATHSPEC>
```
1. If _PATHSPEC_ ends in ".nc" (i.e., if it is the full path to a NetCDF file)
   then we compute baroclinic tide SLA predictions and write them in a file 
   called `PATHSPEC_hret.nc`. The SWOT sample data, RADS passfile, and GDP data
   file formats are currently supported.
   
2. If _PATHSPEC_ ends in ".txt" then we treat the contents of _PATHSPEC_
   as a list of files and generate predictions, as above.
   
3. If _PATHSPEC_ ends in "/" then we compute tide predictions for all
   the files with names ending in ".nc" in this directory, as above.

4. If _PATHSPEC_ is omitted, then the default value is used, _PATHSPEC=/home/jovyan/DEMO_FILES/_.

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
format, are output to the console for you to check. The times are also
written prominently on the sample figures created by this program.
