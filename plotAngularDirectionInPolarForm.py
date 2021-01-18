import os
import matplotlib.pyplot  as plt
import numpy as np
PathToARIANNAanalysis = os.environ['ARIANNAanalysis']

def angle_data(data):
    Angles = np.asarray(data)
    Zen = []
    for i in range(len(Angles)):
        Zen.append(float(Angles[i][0]))
    Azi = []
    for i in range(len(Angles)):
        Azi.append(float(Angles[i][1]))
    return Zen, Azi


file = PathToARIANNAanalysis + '/data/cc_lpdas_stn51.npy'
data = np.load(file,allow_pickle=True,encoding='bytes')
times = data[0]
zen, azi = angle_data(data[1])

fig, ax = plt.subplots(1, 1,figsize=(11, 7),sharex=True)
ax.scatter(azi, zen, s=50,alpha=0.25)
ax.set_ylabel(r'RF zenith [$\degree$]')
ax.set_xlabel(r'RF azimuth [$\degree$]')
ax.grid(True)
ax.tick_params(axis='both', which='both', direction='in')
#ax.set_xlim(300,320)   # feel free to constrict the axis with this line
fig.tight_layout()

save = PathToARIANNAanalysis + '/plots/angularDirectionScatter.png'
fig.savefig(save)
# pdf versions of the figures are nice for LaTeX documents as they preserve resolution better than png files when zooming
save = PathToARIANNAanalysis + '/plots/angularDirectionScatter.pdf'
fig.savefig(save)

### I wouldn't use this polar plot, but you can implement it like this
# Radius is zenith, angle is azimuth
fig = plt.figure()
ax = fig.add_subplot(111, projection='polar')
ax.plot(np.deg2rad(np.asarray(azi)), zen, 'o', MarkerSize=8,alpha=0.25)
save = PathToARIANNAanalysis + '/plots/angularDirectionPolar.png'
fig.savefig(save)
save = PathToARIANNAanalysis + '/plots/angularDirectionPolar.pdf'
fig.savefig(save)


plt.show()
