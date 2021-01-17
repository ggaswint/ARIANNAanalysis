<img src="https://pbs.twimg.com/media/EkYqKtbWsAYwkGC?format=jpg&name=large"/>

# ARIANNAanalysis

Plotting and Data Analysis scripts for the ARIANNA collaboration. Anyone looking to study large data sets from the ARIANNA detectors such as triggering, livetime, or essentially any measurable, and or looking to extend the ARIANNA neutrino reconstruction techniques would benefit from the scripts within this package. See *Installation Instructions and Additional Requirements* at the end of this text.

# Example Figures
These plots show some of the figures produced from the scripts within this package. They have also been published at: https://iopscience.iop.org/article/10.1088/1748-0221/15/09/P09039

<img src="https://pbs.twimg.com/media/Er90lZ-VgAIQKtF?format=jpg&name=large"/>

<img src="https://pbs.twimg.com/media/Er90oEhXAAIlPat?format=jpg&name=large"/>

# Software

**angularPolarPlot.py:**
Takes the saved .npy RF angular data from the script *getAngularReconstructionData.py* or *getAngularReconstructionDataNurInput.py* and makes a polar plot of the reconstructed RF angular directions.

**getAnglesFromSPICE.py**
Used to find match the arrival and launch angular direction given the depth of the SPICEcore pulser for ARIANNA station 51 reciever

**getAngularReconstructionData**
Takes the .root file that contains the 2018 SPICEcore data from station 51 and stores the reconstructed the angular direction in a .npy file format. Requires pyROOT to be installed. The requirement for pyROOT is being phased out of the ARIANNA data analysis, and so it recommended that instead on uses the .nur files which are simply conversions from .root to .nur. These .nur files do not require additional software to be installed. A second script (*getAngularReconstructionDataNurInput.py*) is provided for the .nur input file type.


# Installation Instructions and Additional Requirements

Required packages include NuRadioReco and NuRadioMC along with all their dependencies. These two packages and a list of their required dependencies are found on GitHub at:
https://github.com/nu-radio/NuRadioReco
and
https://github.com/nu-radio/NuRadioMC .
To install, follow the instruction on the wiki page found within the ReadMe file.

The ARIANNA data is stored as ROOT files. NuRadioReco has its own file type, ".nur", along with a module to convert .root files to .nur files. If you can find and work exclusively with the .nur files then ROOT is unnecessary. The scripts below come in two versions when applicable (one version that reads .nur input files and the other for .root input files). However, I highly recommend the installation of ROOT and the other required software for ARIANNA data processing (instruction below). One reason is that the python astropy package is often being upgraded which can cause AttributeErrors when using .nur files that were produced from earlier versions of NuRadioReco. The best way to fix this is to remake the .nur files from the original .root files.

# Installation of ROOT on Unix type computers with python3

The ARIANNA data processing software (snowShovel) was constructed to work with ROOT. Both pieces of software will need to be installed to manage ARIANNA data at the lowest level. This process can be tedious and often bug out multiple times before success. I provide a standard command sequence that should work for most Unix type systems. Note that snowShovel6 along with root 6.18.00 are two verified versions that work together and with NuRadioReco/python3. If you do not have python3 set as your default, DO IT! Don't live in the past.

This particular sequence of commands proved successful for my Ubuntu 18.04.3 machine.


$ git clone --branch v6-22-00-patches https://github.com/root-project/root.git root_src
$ mkdir root_build root_install && cd root_build
$ cmake -DCMAKE_INSTALL_PREFIX=../root_install ../root_src # && check cmake configuration output for warnings or errors
$ cmake --build . -- install -j4 # if you have 4 cores available for compilation
$ source ../root_install/bin/thisroot.sh # or thisroot.{fish,csh}

Clone the repo
    $ cd # got to home directoy
    $ mkdir root-6.18.00
    $ git clone --branch v6-18-00-patches https://github.com/root-project/root.git root_src

Make a directory for building

    $ mkdir build
    $ cd build

Run cmake and make

    $ cmake ../root
    $ make -j8

Setup and run ROOT

    $ source bin/thisroot.sh
    $ root
