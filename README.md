<img src="https://pbs.twimg.com/media/EkYqKtbWsAYwkGC?format=jpg&name=large" align="right"  height="50"/>

# ARIANNAanalysis

Plotting and Data Analysis scripts for the ARIANNA collaboration.


Required packages include NuRadioReco and NuRadioMC along with all their dependencies. These two packages and a list of required dependencies are found on GitHub at:
https://github.com/nu-radio/NuRadioReco
and
https://github.com/nu-radio/NuRadioMC
To install, follow the instruction on the wiki page found within the ReadMe file.

The ARIANNA data is stored as ROOT files. NuRadioReco has its own file type, ".nur", along with a module to convert .root files to .nur files. If you can find and work exclusively with the .nur files then ROOT is unnecessary. The scripts below come in two versions when applicable (one version that reads .nur input files and the other for .root input files). I highly recommend the installing of ROOT and the other required software for ARIANNA data processing.

Note that the python astropy package is often being upgraded which can cause AttributeErrors when using .nur files that were produced from earlier versions of NuRadioReco. To fix this, the .nur files will need to be remade using the latest version of NuRadioReco.

**angularPolarPlot.py:**
Takes the saved .npy RF angular data from the script *getAngularReconstructionData.py* or *getAngularReconstructionDataNurInput.py* and makes a polar plot of the reconstructed RF angular directions.

**getAnglesFromSPICE.py**
Used to find match the arrival and launch angular direction given the depth of the SPICEcore pulser for ARIANNA station 51 reciever

**getAngularReconstructionData**
Takes the .root file that contains the 2018 SPICEcore data from station 51 and stores the reconstructed the angular direction in a .npy file format. Requires pyROOT to be installed. The requirement for pyROOT is being phased out of the ARIANNA data analysis, and so it recommended that instead on uses the .nur files which are simply conversions from .root to .nur. These .nur files do not require additional software to be installed. A second script (*getAngularReconstructionDataNurInput.py*) is provided for the .nur input file type.
