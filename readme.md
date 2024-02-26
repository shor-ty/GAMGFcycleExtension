# GAMG Extension and modularity extension
This repository uses extends the standard GAMG solver within OpenFOAM by the F-Cycle. Hence, the user can specify if one wants to use the common V-cycle or the new implemented F-cycle. Tests showed that the F-cycle does not have any benefit in terms of speedup to the simulation cases but it is also stated that the F-cycle benefit is mainly related to stiff problems. 

The library also extends the GAMG library in a more modular way offering a better extension for other cycles and tests.

## OpenFOAM Versions
* OpenFOAM-v2212
* OpenFOAM-v2306
* OpenFOAM-v2312

## How to use the library?
1. Clone the library into a directory you want
> git clone https://github.com/shor-ty/GAMGExtension.git
2. Enter the directory
> cd GAMGExtension
3. Source your OpenFOAM version (as an example):
> source <pathToYourOpenFOAMInstallation>/etc/bashrc
4. Compile the new library
> wmake libso
5. Add the following line to the __controlDict__ of your simulation case
> libs (GAMGFcycleExtension);
6. Go to your __fvSolution__ file and change the solver name from GAMG to GAMGExtension
7. Add the new keyword `cycleMode  Vcycle;` or `cycleMode Fcycle;` (VCycle is default)
8. Done

## Credit
This repository was re-build and modified by Tobias Holzmann. The base idea however, came from Minhao Xu, see gitlab issue #2054.
Based on his intension, Tobias Holzmann rebuild the code on basis of the latest ESI-OpenFOAM version.
