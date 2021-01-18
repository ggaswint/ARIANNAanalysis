import matplotlib.pyplot  as plt
import numpy as np
import os
from scipy import signal
PathToARIANNAanalysis = os.environ['ARIANNAanalysis']


# Uses 2019 data instead of 2018 data because the transmitters signal amplitude changed over time, and the expected efield was given in 2019. Another SPICEcore run was performed in 2019 and therefor this data was used
file_names = ['expectedEfieldSpice2019.npy','predictedEfieldSpice2019.npy'] #singles 2019. Second file taken from when pulser was at 1200m depth, corresponds to launch angle used for first file (anechoic data)

# 1/r scaling factors
factor_anechoic = 3.0 # Seperation distance between transmitter and LPDA receiver for anechoic test
factor_spice = 1256.6524228044634 # Seperation distance between transmitter and LPDA receiver for this SPICE event

alphas = [1.0,1.0]
label = ['anechoic','SPICE']

fig, ax = plt.subplots(2, 1,figsize=(9, 6)) # test
file = PathToARIANNAanalysis + '/data/'+file_names[0]
data = np.load(file,allow_pickle=True)
etheta_anechoic = data[1]
ephi_anechoic = data[2]
x_anechoic = data[0]

file = PathToARIANNAanalysis + '/data/'+file_names[1]
data2 = np.load(file,allow_pickle=True)

etheta_anechoic = etheta_anechoic*factor_anechoic
ephi_anechoic = ephi_anechoic*factor_anechoic
data2[1] = data2[1]*factor_spice
x_spice = data2[0]
etheta_spice = data2[1]
ephi_spice = data2[2]*factor_spice


time0 = np.argmax(signal.correlate(etheta_anechoic, etheta_anechoic))
time_c = np.argmax(signal.correlate(etheta_anechoic, etheta_spice))
delta_t = x_spice[1] - x_spice[0]

time_offset = (time0-time_c)*delta_t

print(np.max(etheta_anechoic)/np.max(etheta_spice))
print(np.max(etheta_spice)/np.max(etheta_anechoic))

time_offset = x_anechoic[np.argmax(etheta_anechoic)] - x_spice[np.argmax(etheta_spice)]

ax[0].plot(x_anechoic,etheta_anechoic,label=label[0],alpha=alphas[0],color='#1f77b4')
ax[0].plot(x_spice+time_offset-4.951-.3,etheta_spice,label=label[1],alpha=alphas[0],color='#ff7f0e')

ax[0].set_ylabel('Electric Field @ 1m [V/m]')
ax[0].set_xlabel('time [ns]')
delta_t1 = x_anechoic[1] - x_anechoic[0]
freqs = np.fft.rfftfreq(len(etheta_anechoic), delta_t1)
ax[1].plot(freqs*1e3, 50.0*np.abs(np.fft.rfft(etheta_anechoic, norm='ortho')))  #50 to go from "50GHz bins" to GHz
delta_t2 = x_spice[1] - x_spice[0]
freqs = np.fft.rfftfreq(len(etheta_spice), delta_t2)
ax[1].plot(freqs*1e3, 50.0*np.abs(np.fft.rfft(etheta_spice, norm='ortho')))
ax[1].tick_params('y', colors='k')
ax[1].grid(False)
ax[1].set_xlabel('[MHz]')
ax[1].set_ylabel('[V*GHz/m]')
ax[1].set_xlim(80.0,300.0)
ax[0].set_xlim(280.0,420.0)
ax[0].legend()



fig.tight_layout()
save = PathToARIANNAanalysis + '/plots/efieldExpectedVPredictedAtPulserDepth1200m.png'
fig.savefig(save)
save = PathToARIANNAanalysis + '/plots/efieldExpectedVPredictedAtPulserDepth1200m.pdf'
fig.savefig(save)



plt.show()
