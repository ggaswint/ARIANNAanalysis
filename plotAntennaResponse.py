import numpy as np
import matplotlib.pyplot  as plt
from NuRadioReco.utilities import units
from NuRadioReco.detector import antennapattern
import os

PathToARIANNAanalysis = os.environ['ARIANNAanalysis']

antenna_provider = antennapattern.AntennaPatternProvider()

def plotAntennas(ax,zen_arrival,azi_arrival,antennaType):
	ff = np.linspace(0,1,1000)

	LPDA = antenna_provider.load_antenna_pattern(antennaType)
	VEL = LPDA.get_antenna_response_vectorized(ff,zen_arrival*units.deg, azi_arrival*units.deg, 180*units.deg, 0, 90.0*units.deg, 8.21*units.deg)

	ax[0].semilogy(ff,np.abs(VEL['theta']),label = 'eTheta')
	ax[1].semilogy(ff,np.abs(VEL['phi']),label = 'ePhi')
	ax[0].legend()
	ax[1].legend()
	ax[0].set_ylabel('|Antenna Effective Height| [m]')
	ax[0].set_xlabel('frequency [GHz]')
	ax[1].set_xlabel('frequency [GHz]')

fig1, ax1 = plt.subplots(1,2,figsize=(12, 8),sharex='col',sharey='row')
plotAntennas(ax1,135,312,"createLPDA_100MHz_InfFirn")
fig1.tight_layout()
fig1.savefig(PathToARIANNAanalysis + '/plots/lpdaResponse.png')
fig1.savefig(PathToARIANNAanalysis + '/plots/lpdaResponse.pdf')

fig2, ax2 = plt.subplots(1,2,figsize=(12, 8))
plotAntennas(ax2,135,312,"bicone_v8_InfFirn")
fig2.tight_layout()
fig2.savefig(PathToARIANNAanalysis + '/plots/dipoleResponse.png')
fig2.savefig(PathToARIANNAanalysis + '/plots/dipoleResponse.pdf')

plt.show()
