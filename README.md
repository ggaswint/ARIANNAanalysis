<img src="https://pbs.twimg.com/media/EkYqKtbWsAYwkGC?format=jpg&name=large"/>

# ARIANNAanalysis

Plotting and Data Analysis scripts for the ARIANNA collaboration. Anyone looking to study large data sets from the ARIANNA detectors such as triggering, livetime, or essentially any measurable, and or looking to extend the ARIANNA neutrino reconstruction techniques would benefit from the scripts within this package. See *Installation Instructions and Additional Requirements* at the end of this text and do not forget to add a link in .bashrc or .bash_profile (example at very end of this text).

# Example Figures
These plots show some of the figures produced from the scripts within this package. They have also been published at: https://iopscience.iop.org/article/10.1088/1748-0221/15/09/P09039

<img src="https://pbs.twimg.com/media/Er90lZ-VgAIQKtF?format=jpg&name=large"/>

<img src="https://pbs.twimg.com/media/Er90oEhXAAIlPat?format=jpg&name=large"/>

# Software

**angularPolarPlot.py:**
Takes the saved .npy RF angular data from the script *getAngularReconstructionData.py* or *getAngularReconstructionDataNurInput.py* and makes a polar plot of the reconstructed RF angular directions.

**getAnglesFromSPICE.py**
Used to match the arrival and launch angular direction given the depth of the SPICEcore pulser for an ARIANNA station 51 receiver. Imported by other scripts in ARIANNAanalysis.

**getAngularReconstructionDataRootInput** and **getAngularReconstructionDataNurInput**
Reconstructs events arrival direction. In particular it is setup with and example that takes the .root file or equivalent .nur file that contains the 2018 SPICEcore data from station 51 and stores the reconstructed angular direction in a .npy file format. Using root input file requires pyROOT to be installed. The requirement for pyROOT is being phased out of the ARIANNA data analysis, and so it recommended that instead use the .nur files which are simply conversions from .root to .nur. These .nur files do not require additional software to be installed. However it is still recommended to have ROOT installed in order to make these file conversions and have some backwards compatibility. See **writeNurFileFromRootFile.py** for a conversion script. A second script (*getAngularReconstructionDataNurInput.py*) is provided for the .nur input file type.

**getElectricFieldDataAtReciever.py**
Reconstructs the electric field at the ARIANNA station. In particule using the SPICE data from 2018 which involves matching date stamps to pulser depths and getting the expected or reconstructed arrival direction from this. This script sort of builds on top of getAngularReconstructionData in the sense that it first requires a reconstruction of the RF direction. This is necessary to deconvolve out the antenna responses. Note that their is no ROOT input file type example. To see how to convert this to a ROOT input instead of nur input, check out differences between **getAngularReconstructionData** and **getAngularReconstructionDataNurInput**.


# Installation Instructions and Additional Requirements

Required packages include NuRadioReco and NuRadioMC along with all their dependencies. These two packages and a list of their required dependencies are found on GitHub at:
https://github.com/nu-radio/NuRadioReco
and
https://github.com/nu-radio/NuRadioMC .
To install, follow the instruction on the wiki page found within the ReadMe file.

The ARIANNA data is stored as ROOT files. NuRadioReco has its own file type, ".nur", along with a module to convert .root files to .nur files. If you can find and work exclusively with the .nur files then ROOT is unnecessary. The scripts below come in two versions when applicable (one version that reads .nur input files and the other for .root input files). However, I highly recommend the installation of ROOT and the other required software for ARIANNA data processing (instruction below). One reason is that the python astropy package is often being upgraded which can cause AttributeErrors when using .nur files that were produced from earlier versions of NuRadioReco. The best way to fix this is to remake the .nur files from the original .root files.

# Installation of ROOT on Unix type computers with python3

The ARIANNA data processing software (snowShovel) was constructed to work with ROOT. Both pieces of software will need to be installed to manage ARIANNA data at the lowest level. This process can be tedious and often bug out multiple times before success. I provide a standard command sequence that should work for most Unix type systems. Note that snowShovel6 along with root 6.18.00 are two verified versions that work together and with NuRadioReco/python3. If you do not have python3 set as your default, DO IT! Don't live in the past.

This particular sequence of commands proved successful for my Ubuntu 18.04.3 machine. Hopefully copy paste into the terminal proves successful for you as well.


$ git clone --branch v6-22-00-patches https://github.com/root-project/root.git root_src
$ mkdir root_build root_install && cd root_build
$ cmake -DCMAKE_INSTALL_PREFIX=../root_install ../root_src # && check cmake configuration output for warnings or errors
$ cmake --build . -- install -j4 # if you have 4 cores available for compilation
$ source ../root_install/bin/thisroot.sh # or thisroot.{fish,csh}

Make root folder and clone GitHub repo

    $ cd
    $ mkdir root-6.18.00 && cd root-6.18.00
    $ git clone --branch v6-18-00-patches https://github.com/root-project/root.git root_src

Make the build file and insure it is using python3! Note that *-DPYTHON_EXECUTABLE=/path/to/desired/python* needs to be modified with your correct python path to python3. If your system uses python3 by default (i.e. python foo.py runs in python3) then you should be able to get the path with typing *which python* in the command line.

    $ cmake -DPYTHON_EXECUTABLE=/path/to/desired/python -DCMAKE_INSTALL_PREFIX=../root_install ../root_src

Make build file. Note -j4 means use 4 cores for build. Should be fine.

    $ cmake --build . -- install -j4

You need to tell the terminal about root by "sourceing" it

    $ source ../root_install/bin/thisroot.sh

A better way would be to edit your .bashrc or .bash_profile located in the home directory. At the end of this text I give my .bashrc extra lines to be used. To test root, first simply type root in the terminal. If a root GUI pops up, then it is installed successfully. Type *.q* to quit. Who wants to use this GUI? ew. Now lets see if pyROOT is working. Make a python file and type import ROOT. Execute this with *python filename.py*. If this works then you are all set with the ROOT installation. If either of these are broken, time to google or ask someone for help or simply retry the process (I know sounds crazy but i've seen crazy things in my life).

In the extras folder there are two files called .rootrc and .rootlogon.C. Save these files as is to your home directory. They are crucial for getting MACROS to work between ROOT and snowShovel (the final step).

# Installing snowShovel6

snowShovel is the first framework built for ARIANNA. It contains scripts to communicate with the stations in Antarctica, and for debugging local stations. scripts/online will be a very familiar destination for any scientist working directly with the ARIANNA stations. snowShovel6 is a private code base, as with it comes the power to control the ARIANNA experiment. To download snowShovel type:

    $ svn co https://arianna.ps.uci.edu/svn/repos/snowShovel/trunk snowShovel

Note that a username and password is required, and insure that this link is to the latest version (version 6). Checkout https://arianna.ps.uci.edu/mediawiki/index.php/Local_DAQ_Instructions for detailed instructions on setting up snowShovel. Note that the root and snowShovel packages on that link are outdated and will not work with NuRadioReco, however the setup instructions for snowShovel are identical.

# .bashrc extra links

Here is an example of the extra lines to save in a .bashrc file. Note that */home/geoffrey* will need to be changed to whatever your home directoy is, and insure that all the paths are correct for your setup. Note that their is an additional link to ARIANNAanalysis.

    export HOME="/home/geoffrey"
    export SNS="/home/geoffrey/snowShovel6/snowShovel"
    export SNSscripts="/home/geoffrey/snowShovel6/snowShovel/scripts"
    export ROOTSYS="/home/geoffrey/root-6.18.00/root_install"
    export Nu="/home/geoffrey/NuRadioReco"
    export NuM="/home/geoffrey/NuRadioMC"
    export Radio="/home/geoffrey/radiotools"
    export ARIANNAanalysis="/home/geoffrey/ARIANNAanalysis/ARIANNAanalysis"

    source ${ROOTSYS}/bin/thisroot.sh

    export LD_LIBRARY_PATH=${ROOTSYS}/lib:${LD_LIBRARY_PATH}
    export PYTHONPATH=${ROOTSYS}/lib:${PYTHONPATH}
    export LD_LIBRARY_PATH=${SNS}/lib:${SNS}:${LD_LIBRARY_PATH}
    export PYTHONPATH=${HOME}:${PYTHONPATH}
    export PYTHONPATH=${SNS}/lib:${PYTHONPATH}
    export PYTHONPATH=${SNS}/scripts:${PYTHONPATH}
    export PYTHONPATH=${SNS}:${PYTHONPATH}
    export PYTHONPATH=${Nu}/:${PYTHONPATH}
    export PYTHONPATH=${NuM}/:${PYTHONPATH}
    export PYTHONPATH=${Radio}/:${PYTHONPATH}
    export PYTHONPATH=${ARIANNAanalysis}/:${PYTHONPATH}
    export ROOT_INCLUDE_PATH=${SNS}/include:${ROOTSYS}/include/

Do not forget to either restart the command line or source this file with *source .bashrc*
