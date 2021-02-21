import matplotlib.pyplot  as plt
import numpy as np
from radiotools import helper as hp
from radiotools import stats as stat
import helpfulUtils as hu
from mpl_toolkits.axes_grid1 import make_axes_locatable
from scipy.stats import norm
import os
PathToARIANNAanalysis = os.environ['ARIANNAanalysis']

def angle_data(data):
    Angle = np.asarray(data)# L1 data data
    Zen = []
    for i in range(len(Angle)):
        Zen.append(float(Angle[i][0]))
    Azi = []
    for i in range(len(Angle)):
        Azi.append(float(Angle[i][1]))
    return Zen, Azi

def getData(file, dataExpected, reverse=False):
    set1 = np.load(file,allow_pickle=True,encoding='bytes')
    depth = set1[0]
    data = set1[1]
    Zen, Azi = angle_data(data)
    ZenE = []
    for i in range(len(depth)):
        idx = hu.find_nearest(np.asarray(dataExpected[0]), depth[i])
        ZenE.append(dataExpected[1][idx])
    return np.asarray(Zen), np.asarray(Azi), np.asarray(ZenE), np.ones(len(Azi))*312.448284, np.asarray(depth)

def gaussianHist(ax,data,color,label,pos,line,label2):
    (mu, sigma) = norm.fit(data)
    fig2, ax2 = plt.subplots(1, 1)
    n, bins, patches = ax.hist(data,linestyle=line,bins=np.arange(-1.5, 3.6, 0.1),edgecolor=color,fill=False,label=label2,histtype='step',orientation='horizontal')#,weights=weights) histype = bar
    if label2 == 'dipoles':
        ax.text(0.99,0.98,label,horizontalalignment='right',verticalalignment='top',transform=ax.transAxes,color='darkred')
    if label2 == 'lpdas':
        ax.text(0.99,0.90,label,horizontalalignment='right',verticalalignment='top',transform=ax.transAxes,color='midnightblue')
    plt.close(fig2)


# Note this is now a rayleigh distribution
def getMeanSTDStr(data):
    mean = np.mean(data)
    textstr = r"$\mu$ = %.2g$^{\circ}$, $\sigma_{68\%%}$ = %.2g$^{\circ}$" % (mean,stat.quantile_1d(data, np.ones_like(data), 0.68))
    return textstr


def aveError(depth,data):
    depth = np.asarray(depth)
    depth = np.round(depth.astype(float), -1)
    #depth = np.unique(np.round(depth.astype(float), -1))
    means = []
    stds = []
    depths = []
    for d in np.unique(np.round(depth.astype(float), -1)):
        mask = depth==d
        means.append(data[mask].mean())
        stds.append(data[mask].std())
        depths.append(d)
    return means, stds, depths

def get_single_angle(zenith_reco,azimuth_reco,zenith_exp,azimuth_exp):
    vspherical_to_cartesian = np.vectorize(hp.spherical_to_cartesian,otypes=[np.ndarray])
    vget_angle = np.vectorize(hp.get_angle)
    v1 = vspherical_to_cartesian(zenith_reco, azimuth_reco)
    v2 = vspherical_to_cartesian(zenith_exp, azimuth_exp)
    return vget_angle(v1, v2)

# A few depths of interest that marks when the reflection coefficient becomes 0.5 or 0.1 or is TIR
R_50 = 937.63763764
R_10 = 1180.98098098
R_TIR = 918.91891892


max_diff = 10.0
datafile = PathToARIANNAanalysis + '/data/expectedArrivalDirectionSpice2018smooth.npy'
dataExpected = np.load(datafile,allow_pickle=True)


fig, ax = plt.subplots(1, 1,figsize=(11, 7),sharex=True)

file = PathToARIANNAanalysis + '/data/reconstructedAngularDirectionFromDipoles1180mDepthAndBelowWithBandpassFilter80to300MHz.npy'
file2 = PathToARIANNAanalysis + '/data/reconstructedAngularDirectionFromLPDAs938mDepthAndBelowWithBandpassFilter80to300MHz.npy'

# R for reconstructed, E for expected
zenR_d, AziR_d, ZenE_d, AziE_d, depth_d = getData(file,dataExpected)
zenR_l, AziR_l, ZenE_l, AziE_l, depth_l = getData(file2,dataExpected)

angles_d = get_single_angle(zenR_d,AziR_d,ZenE_d,AziE_d)
angles_l = get_single_angle(zenR_l,AziR_l,ZenE_l,AziE_l)

means_d, errors_d, depths_d = aveError(depth_d,np.asarray(angles_d))
means_l, errors_l, depths_l = aveError(depth_l,np.asarray(angles_l))


size_marker = 40
ax.scatter(depth_d,angles_d,s=size_marker,marker='s',color='red',alpha=0.25)
ax.scatter(depth_l,angles_l,s=size_marker,marker='>',color='deepskyblue',alpha=0.25)

ax.scatter(depths_d,means_d,s=size_marker,marker='s',color='darkred',label='Dipoles ave.')
ax.scatter(depths_l,means_l,s=size_marker,marker='>',color='midnightblue',label='LPDAs ave.')

ax.fill_between([900,946],[17,17],facecolor='none',hatch='x', edgecolor='black',linewidth=0.0)
ax.fill_between([900,946],[-17,-17],facecolor='none',hatch='x', edgecolor='black',linewidth=0.0)
ax.fill_between([1233,1270.1],[2.5,2.5],facecolor='none',hatch='x', edgecolor='black',linewidth=0.0)
ax.fill_between([1233,1270.1],[-17,-17],facecolor='none',hatch='x', edgecolor='black',linewidth=0.0)
ax.fill_between([1520,1547],[17,17],facecolor='none',hatch='x', edgecolor='black',linewidth=0.0)
ax.fill_between([1520,1547],[-17,-17],facecolor='none',hatch='x', edgecolor='black',linewidth=0.0)

ax.axhline(y=0.0,color='black')
ax.set_ylabel('$\Delta$ anglular direction [$\degree$]')

ax.set_ylim(-0.49,3.0)
ax.set_xlim(800,1700)

ax.plot([1180,1180],[-5,5],color='darkred')
ax.plot([938,938],[-5,5],color='midnightblue')

ax.annotate(s='', xy=(1250.0,2.7), xytext=(1180,2.7),ha='center',va='top',color='darkred',arrowprops=dict(arrowstyle='->',color='darkred'))
ax.text(0.53,0.98,s=r'R$\leq$0.1',horizontalalignment='right',verticalalignment='top',transform=ax.transAxes,color='darkred')

ax.annotate(s='', xy=(1008.0,2.7), xytext=(938,2.7),ha='center',va='top',color='midnightblue',arrowprops=dict(arrowstyle='->',color='midnightblue'))
ax.text(0.26,0.98,s=r'R$\leq$0.5',horizontalalignment='right',verticalalignment='top',transform=ax.transAxes,color='midnightblue')



fig.tight_layout()
fig.subplots_adjust(wspace = 0, hspace = 0.02)


divider = make_axes_locatable(ax)
axHisty0 = divider.append_axes("right", 2.2, pad=0.1,sharey=ax)

mask_dipole = (depth_d >= 1180.0)
mask_lpda = (depth_l >= 938.0)

print('stats without depth cuts:')
print('Lpdas')
textstr = getMeanSTDStr(angles_l)
print(textstr)
print('Dipoles')
textstr = getMeanSTDStr(angles_d)
print(textstr)

gaussianHist(axHisty0,angles_l[mask_lpda],'midnightblue',getMeanSTDStr(angles_l[mask_lpda]),[-4.35,160],'--','lpdas')
gaussianHist(axHisty0,angles_d[mask_dipole],'darkred',getMeanSTDStr(angles_d[mask_dipole]),[0,165],'-','dipoles')


ax.text(0.99,0.98,r'$\blacksquare$ dipoles',horizontalalignment='right',verticalalignment='top',transform=ax.transAxes,color='darkred')
ax.text(0.99,0.90,r'$\blacktriangleright$ lpdas',horizontalalignment='right',verticalalignment='top',transform=ax.transAxes,color='midnightblue')

ax.xaxis.set_ticks(np.arange(800,1700,100))


axHisty0.tick_params(labelleft=False)
axHisty0.xaxis.set_ticks(np.arange(250,1250,250))

axHisty0.set_xlabel('Number of entries')


save = PathToARIANNAanalysis + '/plots/spaceAngularSpiceData.pdf'
fig.savefig(save)
save = PathToARIANNAanalysis + '/plots/spaceAngularSpiceData.png'
fig.savefig(save)

plt.show()
