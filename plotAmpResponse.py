import numpy as np
import matplotlib.pyplot  as plt
import NuRadioReco.detector.ARIANNA.analog_components as ac
import os

PathToARIANNAanalysis = os.environ['ARIANNAanalysis']

fig, ax = plt.subplots(1,1,figsize=(12, 8),sharex='col',sharey='row')
ff = np.linspace(0.01,1.2,1000)

amp_response = '100'
ax.plot(ff,20*np.log10(np.abs(ac.amplifier_response[amp_response]['gain'](ff))),linewidth=3,linestyle='solid',label='100-series',color='darkblue')

amp_response = '200'
ax.plot(ff,20*np.log10(np.abs(ac.amplifier_response[amp_response]['gain'](ff))),linewidth=3,linestyle='dotted',label='200-series',color='darkorange')

amp_response = '300'
ax.plot(ff,20*np.log10(np.abs(ac.amplifier_response[amp_response]['gain'](ff))),linewidth=3,linestyle='dashed',label='300-series',color='darkgreen')

ax.set_ylim(0,75)
ax.set_xlim(0.0,1.2)
ax.legend()
ax.set_ylabel('gain [dB]')
ax.set_xlabel('frequency [GHz]')

fig.tight_layout()
fig.savefig(PathToARIANNAanalysis + '/plots/ampResponse.png')
fig.savefig(PathToARIANNAanalysis + '/plots/ampResponse.pdf')

plt.show()
