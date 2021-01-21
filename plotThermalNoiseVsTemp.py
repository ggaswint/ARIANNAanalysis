import numpy as np
import matplotlib.pyplot as plt
from NuRadioReco.utilities import units
from scipy import constants
import os

PathToARIANNAanalysis = os.environ['ARIANNAanalysis']



K_B = constants.k
def thermalNoisePower(temp, bandwidth):
    return bandwidth * K_B * temp

def thermalNoiseVrms50ohm(temp, bandwidth):
    R = 50.0 # Ohms, common impedence
    # Note some definitions remove the factor of 4, though I think this is technically wrong
    return (4* bandwidth * K_B * temp * R)**0.5 * units.V

BW = (500.0 - 50.0) * units.MHz / units.Hz
T = 250 # Kelvin
print("Vrms for bandwidth " + str(BW*units.Hz/units.MHz) + " MHz and temp: " + str(T) + " K is: " + str(thermalNoiseVrms50ohm(T,BW)/units.mV) + " mV")
BW = (800.0 - 50.0) * units.MHz / units.Hz
T = 300 # Kelvin
print("Vrms for bandwidth " + str(BW*units.Hz/units.MHz) + " MHz and temp: " + str(T) + " K is: " + str(thermalNoiseVrms50ohm(T,BW)/units.mV) + " mV")


temperatures = np.linspace(200,400,100) # Kelvin
bandwidth = (500.0 - 50.0) * units.MHz / units.Hz
func = np.vectorize(thermalNoiseVrms50ohm)

fig1, ax1 = plt.subplots(1, 1, sharex=False, figsize=(5, 5))
ax1.plot(temperatures,func(temperatures,bandwidth))


ax1.set_ylabel(r'Vrms [V]')
ax1.set_xlabel(r'T [K]')

fig1.tight_layout()
fig1.savefig(PathToARIANNAanalysis + '/plots/thermalNoiseVrmsPerTemp.png')
fig1.savefig(PathToARIANNAanalysis + '/plots/thermalNoiseVrmsPerTemp.pdf')

plt.show()
