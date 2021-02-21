from NuRadioReco.utilities import units
from NuRadioMC.EvtGen.generator import generate_eventlist_cylinder
import h5py
import numpy as np
import os

PathToARIANNAanalysis = os.environ['ARIANNAanalysis']

if __name__== "__main__":
    # define simulation volume
    zmin = -2.7 * units.km  # the ice sheet at South Pole is 2.7km deep
    zmax = 0 * units.km
    rmin = 0 * units.km
    rmax = 2 * units.km

    volume = {"fiducial_rmin": rmin, "fiducial_rmax" : rmax, "fiducial_zmin" : zmin, "fiducial_zmax" : zmax}
    # generate one event list at 1e19 eV with 10 neutrinos

    saveFile = PathToARIANNAanalysis + '/data/tenNeutrinosAt1e19.hdf5'
    generate_eventlist_cylinder(saveFile, 10, 1e19 * units.eV, 1e19 * units.eV, volume)
