# Introduction
## Overview
This project contains the PSLab code used by the IceCube collaboration for some point-source (PS) analyses, and a script example of how to run a PS analysis. The data used in this project was collected by IceCube across 10 years (2008-2018), and can be downloaded [here](https://icecube.wisc.edu/data-releases/2021/01/all-sky-point-source-icecube-data-years-2008-2018/). The data release also contains information about the detector response (effective area, reconstructed angular error and energy estimation) and the kinematic neutrino/muon angle, that will be needed for the signal simulations. A detailed documentation about the aforementioned data release can be found ([here](https://arxiv.org/abs/2101.09836)). There is no need to download the files of the data release to run this analysis, as they can be found in the folder *data_release* of this project. The description of the method for point-source searches in IceCube can be found [here](https://www.sciencedirect.com/science/article/pii/S0927650508000303?via%3Dihub) (time-independent searches) and [here](https://www.sciencedirect.com/science/article/pii/S0927650510000241?via%3Dihub) (time-dependent searches).

## Description of Files in This Project
This project contains the following files:
  * **data_release** (folder): it contains the 10-yr IceCube data and the information on the detector response, as mentioned in the previous section.
  * **psLab** (folder): it contains the PSLab code, with the relevant classes to run a point-source search.
  * **README.md**: this readme file, with an introduction to the codes and instructions to run the search.
  * **catalogue.txt**: the file containing a toy catalogue of sources to analyse.
  * **createListOfTimesForScrambling.py**: a python script needed to create a list of uptimes of the detector, that will be used in the scramble method to produce background pseudo-experiments (see below).
  * **make_datafiles.py**: a python script to read the information of the data stored in the *data_release* folder, and produce data files in a suitable format for the analysis codes (see below).
  * **detector_response.py**: a python file containing the classes to read and interpret the information on the detector response store in the *data_release* folder.
  * **start.sh**: a bash script to load the environment. **IMPORTANT**: to properly run this project, open the file and change the *_LAB_MAIN_DIR* variable (line 7) to point to the *full path* of your psLab folder.
  * **PSTimeDepAna.py**: the main python script to run a time-dependent point-source analysis.
  * **stackingAnalysis.py**: the main python script to run a stacking analysis.
  * **SimpleAnalysisStack.C** and **SimpleAnalysisPS.C**: auxiliary C macros called by *PSTimeDepAna.py* and *stackingAnalysis.py* respectively, to iterate the background scramble, signal injection and likelihood maximization.
  * **load_ark.C**: an auxiliary C macro called by *PSTimeDepAna.py* and *stackingAnalysis.py* to set some basic important information for the analysis.
  * **tree_loader.C**: an auxiliary C macro called by *PSTimeDepAna.py* and *stackingAnalysis.py* to load the data.
  
## Archival Searches Based on PSLab
The PSLab code was used by the IceCube collaboration for the following analyses:
  1. *First Neutrino Point-Source Results From the 22-String IceCube Detector*, the IceCube collaboration, 2009, [link](https://iopscience.iop.org/article/10.1088/0004-637X/701/1/L47)
  1. *Time-Integrated Searches for Point-like Sources of Neutrinos with the 40-String IceCube Detector*, the IceCube collaboration, 2011, [link](https://iopscience.iop.org/article/10.1088/0004-637X/732/1/18)
  2. *Time-Dependent Searches for Point Sources of Neutrinos with the 40-String and 22-String Configurations of IceCube*, the IceCube collaboration, 2012, [link](https://iopscience.iop.org/article/10.1088/0004-637X/744/1/1)
  3. *Search for Time-independent Neutrino Emission from Astrophysical Sources with 3 yr of IceCube Data*, the IceCube collaboration, 2013, [link](https://iopscience.iop.org/article/10.1088/0004-637X/779/2/132)
  4. *Search for time-dependent neutrino sources with IceCube data from 2008 to 2012*, the IceCube collaboration, 2015, [link](https://iopscience.iop.org/article/10.1088/0004-637X/807/1/46)
  5. *Search for correlations between the arrival directions of IceCube neutrino events and ultrahigh-energy cosmic rays detected by the Pierre Auger Observatory and the Telescope Array*, the IceCube collaboration, 2016, [link](https://iopscience.iop.org/article/10.1088/1475-7516/2016/01/037)
  6. *Search for annihilating dark matter in the Sun with 3 years of IceCube data*, the IceCube collaboration, 2017, [link](https://link.springer.com/article/10.1140/epjc/s10052-017-4689-9)
  7. *Neutrino emission from the direction of the blazar TXS 0506+056 prior to the IceCube-170922A alert*, the IceCube collaboration, 2018, [link](https://www.science.org/doi/10.1126/science.aat2890)
  8. *Time-integrated Neutrino Source Searches with 10 Years of IceCube Data*, the IceCube collaboration, 2020, [link](https://journals.aps.org/prl/abstract/10.1103/PhysRevLett.124.051103)
  9. *A Search for IceCube Events in the Direction of ANITA Neutrino Candidates*, the IceCube collaboration, 2020, [link](https://iopscience.iop.org/article/10.3847/1538-4357/ab791d)
  10. *A Search for Time-dependent Astrophysical Neutrino Emission with IceCube Data from 2012 to 2017*, the IceCube collaboration, 2021, [link](https://iopscience.iop.org/article/10.3847/1538-4357/abe7e6)
  11. *Search for Multi-flare Neutrino Emissions in 10 yr of IceCube Data from a Catalog of Sources*, the IceCube collaboration, 2021, [link](https://iopscience.iop.org/article/10.3847/2041-8213/ac2c7b)
  
For internal analyses, IceCube adopts a signal injection method based on Monte Carlo events, thus different from the method implemented for the release of this code.

# Software Requirements
## ROOT v6.18.04
The ROOT software can be downloaded from [here](https://root.cern/releases/release-61804/). Make sure to download the version 6.18.04 to run this analysis: the compatility of subsequent ROOT versions has not been tested

## Python 3.8.10
Python 3.8.10 must be used, and the following modules are required:
  1. numpy
  2. pandas
  3. os
  4. random
  5. argparse
  6. ROOT
  7. array
  8. sys
  9. copy
  10. math
  11. time
  12. astropy

## ERFA and PAL
The local/equatorial coordinate conversion is based on [PAL](https://github.com/Starlink/pal) (Positional Astronomy Library). Notice that PAL requires that either the SOFA C library or the [ERFA](https://github.com/liberfa/erfa) library are installed.<br />
To install PAL and ERFA, follow the instructions provided [here](https://docs.icecube.aq/icetray/main/projects/cmake/tools/pal.html) to download and install ERFA and PAL.<br />
You can specify the path where ERFA and PAL must be installed by running the configure file with the `--prefix=<user_defined_path>`. In this case, remember to run the PAL configure file with the option `--with-erfa=<path_to_erfa>`, specifying the path where ERFA was installed. In case you do not install the PAL library in the default library path, remember to export the path of `libpal.so` to your LIBPATH variable.<br />
The description paper for PAL is: ["PAL: A Positional Astronomy Library"](https://ui.adsabs.harvard.edu/abs/2013ASPC..475..307J/abstract), Jenness, T. & Berry, D. S., in Astronomical Data Anaysis Software and Systems XXII, Friedel, D. N. (ed), ASP Conf. Ser. 475, p307.

# Prepare Files for Analysis
Some preliminary files must be run _only once_ to generate the files needed for the analysis (i.e. the data file in an appropriate format, the lists of times used for the scramble method to generate background pseudo-experiments).

## Data Files
To generate the data file in the appropriate format, a ROOT tree, and a file containing the start/stop times of the IceCube uptimes, run the following command:

```
python make_datafiles.py @inDir data_release/
```

The script accepts the following parameters:
  * @inDir (required): path to the directory containing the files of the data release.
  * @season (not required): the name of the IceCube season to generate. If not provided, the script generates all the seasons by default.
Open the script and read the help descriptor for more details.

The script output is:
  1. A set of files (ROOT trees) in the folder *data_release*. Each file corresponds to nearly one year of IceCube data, and it is named after the number of strings of the detector (40, 59, 79, 86), and where needed after the year of the data.
  2. A set of files (ROOT trees) in the folder *data_release/GRL*. Each file contains the start and stop moment of the uptime intervals of IceCube in each year, and also the duration of this interval (called *livetime* in the files and equal to the difference between the stop and start times). These files will be used to produce a list of possible uptimes to be used in the scramble method when producing background pseudo-experiments (see next section).

## List of Good Times for Scramble
This part is needed to generate a file containing a list of possible uptimes of the detector. These times are used to produce background pseudo-realizations of the data through the so-called scrmable method: events are assigned a random time from this list, and their equatorial coordinates (right ascension and declination) are corrected accordingly, assuming fixed local coordinates (azimuth, zenith). To produce such files, run the following command:

```
python createListOfTimesForScrambling.py -seed 123
```

The script accepts the following parameters:
  * -seed (required): the seed of the random number generator. It can be any positive integer number.
  * -season (not required): the IceCube season of which the script must generate the possible uptimes. If not provided, the script will run across all the seasons.
  * -scale (not required): a scale factor for the number of random times to be generated in each uptime interval. Defaul is 1000 (do not use less than 1000)

The script output is a list of .txt files in the folder *ListOfTimes*. Each file contains the possible uptimes of the detector for that season, that will be used by the scramble method.

# Running the PS Analysis
This project contains two example codes to run a simple search for point-like neutrino sources:
  1. **PSTimeDepAna.py**: it looks for astrophysical neutrino emission from a single direction in the sky.
  2. **stackingAnalysis.py**: it looks for astrophysical neutrino emission by stacking all the sources in the file *catalogue.txt*.
   
## Preliminary: Load the Environment
To load the environment, run the following command:
```
./start.sh
```
**Important**: before running the above command, make sure to replace the value of the variable *_LAB_MAIN_DIR* (line 7) with the *full path* to your PSLab folder.

Make sure that your environment is using the correct python version (3.8.10) and root version (6.18.04).

## Point-Source Time-Dependent Analysis
This search uses the script **PSTimeDepAna.py** to look for Gaussian-shaped temporal clusters of astrophysical neutrinos from a single direction in the sky.

The code takes as input the following parameters (see alse the argparse in the file for other details):
 * @srcra (required): the Right Ascension (in degrees) of the source to test
 * @srcdec(required): the declination (in degrees) of the source to test
 * @ns (not required, default 0): the (integer) number of signal events to inject when generating pseudo-experiments (ns=0 for background, ns>0 for signal)
 * @gamma (not required, default 2): the spectral index of the astrophysical neutrino flux to use for the signal injection (if ns>0)
 * @sigmaT (not required, default 10): the standard deviation (in days) of the Gaussian flare to be injected (if ns>0)
 * @upLimSigma (not required, default 200): the maximum sigmaT (in days) allowed for the fit
 * @lowLimSigma (not required): the minimum sigmaT (in days) allowed for the fit
 * @nInj (not required, default 10): the number of pseudo-experiments to generate
 * @tInj (not required, default 57000): the mean value (in Modified Julian Date) of the Gaussian flare to be injected (if ns>0)
 * @seed (required): the seed of the random number generator (should be a positive integer)
 * @nSigTrunc (not required, default 4): the number of standard deviation at which the Gaussian flare is truncated
 * @TSthr (not required, default 2): the test statistics threshold that a cluster of events must pass to be considered a flare. This is only relevant for multi-flare search (default is single-flare): to enable multi-flare searches, set the flag `isSingleFlare_` to 'False' in `psLab/llhTimeDep/public/MultiGaussAnalysisFn.h`
 * @unbl (not required, default 0): use the code to test the unblinded data (unbl=1) or to generate pseudo-experiments (unbl=0, default choice)
 
The code searches for the number of signal neutrino candidates ns, the index of the energy spectrum gamma, the central time of the flare, and the flare duration sigmaT that maximize the likelihood function of the analysis.

### Background and Signal Pseudo-Experiments

To generate pseudo-experiments for point-source searches, run the following command:
```
python PSTimeDepAna.py @srcra 180 @srcdec 0.0 @ns 0 @gamma 2 @sigmaT 100 @upLimSigma 200 @nInj 100 @tInj 57000 @seed 1 @unbl 0
```
Replace @srcra and @srcdec with the desired Right Ascension and declination of the source, @ns with the desired number of signal neutrinos to inject (ns=0 for **background** pseudo-experiments, ns>0 for **signal** pseudo-experiments), @gamma with the desired index of the energy spectrum to inject (relevant only if ns>0), @tInj with the desired flare time to inject (relevant only if ns>0) and @sigmaT with its standard deviation (relevant only if ns>0).

This command is used to produce background scrambles of the data (i.e. background pseudo-experiments), to inject astrophysical signal events (if ns>0), and to maximize the likelihood function. It returns a ROOT file containing the best-fit parameterss and the test statistics of each pseudo-experiments.

### Unblinded Data

To run the point-source analysis on unblinded data, use the following command:
```
python PSTimeDepAna.py @srcra 180 @srcdec 0.0 @upLimSigma 200 @nInj 1 @seed 1 @unbl 1
```
Replace @srcra and @srcdec with the desired Right Ascension and declination of the source.

This command is used to perform a likelihood scan on unblinded data. It returns a ROOT file containing the best-fit number of signal candidates and the test statistics.

## Stacking Analysis
This search uses the **stackingAnalysis.py** to search for stacking neutrino emission from the sources listed in *catalogue.txt*. The source fluence is used as a weight in the stacking process.

The code takes as input the following parameters (see alse the argparse in the file for other details):
  * @ns: the (integer) number of signal events to inject when generating pseudo-experiments (ns=0 for background, ns>0 for signal)
  * @gamma: the (fixed) spectral index of the astrophysical neutrino flux to use for the signal injection (if any) and for the likelihood scan
  * @nInj: the number of pseudo-experiments to generate
  * @s: the seed of the random number generator (should be a positive integer)
  * @hemi: the hemisphere to use for the analysis ('north' or 'south'). Given the very different IceCube sensitivity in the two hemispheres, it is worth running separate analyses for northern and southern sources.
  * @unbl: use the code to test the unblinded data (unbl=1) or to generate pseudo-experiments (unbl=0, default choice)

The code searches for the number of signal neutrino candidates that maximizes the stacking likelihood function of the analysis.

### Background and Signal Pseudo-Experiments

To generate pseudo-experiments for single-source searches, run the following command:
```
python stackingAnalysis.py @ns 0 @gamma 2 @nInj 100 @s 1 @hemi north @unbl 0
```
Replace @ns with the desired number of signal neutrinos to inject (ns=0 for **background** pseudo-experiments, ns>0 for **signal** pseudo-experiments), @gamma with the desired spectral index to use in the likelihood scan, @nInj with the desired number of pseudo-experiments to generate, @s with the desired seed for the random number generator, and @hemi with the desired hemisphere to scan ('north' or 'south'). The *ns* parameter must be intended as a stacked sum of signal neutrinos from all the sources.

This command is used to produce background scrambles of the data (i.e. background pseudo-experiments), to inject astrophysical signal events (if ns>0), and to maximize the likelihood function. It returns a ROOT file containing the best-fit stacked number of signal candidates and the test statistics of each pseudo-experiments.

### Unblinded Data

To run the stacking analysis on unblinded data, use the following command:
```
python stackingAnalysis.py @ns 0 @gamma 2 @nInj 1 @s 1 @hemi north @unbl 1
```
Replace @gamma with the desired spectral index to use in the likelihood scan, and @hemi with the desired hemisphere.

This command is used to perform a likelihood scan on unblinded data. It returns a ROOT file containing the best-fit stacked  number of signal candidates and the test statistics.

# Acknowledgements
The PSLab code included in this repository results from the contribution and effort of several people. The main credits for this project go to Teresa Montaruli, who introduced the likelihood idea and supervised the development and maintenance of the code, and Chad Finley, who is the creator of the code and wrote the main classes used for the point-source analyses. We also acknowledge the contributions of Juan Antonio Aguilar, Mike Baker, Anastasia Barbano, Asen Christov, Jon Dumm, Francesco Lucarelli, and Mohamed Rameez.

We also thank the many members of the IceCube collaboration, for their feedback to the developement of this code during the analysis reviews.

The developement of the PSLab code has been supported by the following agencies and institutions: U.S. National Science Foundation-Physics Division, Wisconsin Research Foundation, Swiss National Science Foundation (SNSF), Université de Genève.
