import matplotlib.pyplot  as plt
import numpy as np

import helpfulUtils as hu
from mpl_toolkits.axes_grid1 import make_axes_locatable
from scipy.stats import norm
from radiotools import stats as stat
import os
PathToARIANNAanalysis = os.environ['ARIANNAanalysis']

expectedAzi = 312.448284

def angle_data(data):
    Angle = np.asarray(data)# L1 data data
    Zen = []
    for i in range(len(Angle)):
        Zen.append(float(Angle[i][0]))
    Azi = []
    for i in range(len(Angle)):
        Azi.append(float(Angle[i][1]))
    return Zen, Azi


def getDataDiff(file, dataExpected, reverse=False):
    # Finds the difference between reconstructed and expected angular directions
    # Reverse should be true if the event set is describing a pulser that was being raised versus being lowered into a borehole
    # when full is set to true, the start and stop values are ignored and the entire depth range is taken into consideration
    dataSet = np.load(file,allow_pickle=True,encoding='bytes')
    depth = dataSet[0]
    data = dataSet[1]
    Zen, Azi = angle_data(data)
    ZenExpected = []
    AziExpected = []
    for i in range(len(depth)):
        idx = hu.find_nearest(np.asarray(dataExpected[0]), depth[i])
        ZenExpected.append(dataExpected[1][idx])
        AziExpected.append(dataExpected[2][idx])
    deltaZen = (np.asarray(Zen) - np.asarray(ZenExpected)).astype(float)
    deltaAzi = (np.asarray(Azi) - np.asarray(AziExpected)).astype(float)
    return deltaZen,deltaAzi,depth

def gaussianHist(ax,data,color,label,pos,line,label2):
    (mu, sigma) = norm.fit(data)
    fig2, ax2 = plt.subplots(1, 1)
    n, bins, patches = ax.hist(data,linestyle=line,bins=np.arange(-1.5, 1.6, 0.2),edgecolor=color,fill=False,label=label2,histtype='step',orientation='horizontal')
    if label2 == 'dipoles':
        ax.text(0.99,0.98,label,horizontalalignment='right',verticalalignment='top',transform=ax.transAxes,color='darkred',fontsize=11)
    if label2 == 'lpdas':
        ax.text(0.99,0.90,label,horizontalalignment='right',verticalalignment='top',transform=ax.transAxes,color='midnightblue',fontsize=11)
    plt.close(fig2)

def getMeanSTDStr(data):
    mean = np.mean(data)
    std = np.std(data)
    textstr = r"mean = %.2g$^{\circ}$, STD = %.2g$^{\circ}$" % (round(mean,2),round(std,2)) # The few outlier events make the STD much larger, changed to 68% for better measurement of error
    tweights = np.ones_like(data)
    textstr = r"mean = %.2g$^{\circ}$, $\sigma_{68\%%}$=%.2g$^{\circ}$" % (round(mean,2), stat.quantile_1d(data,tweights,0.68))

    tweights = np.ones_like(data)
    q1 = stat.quantile_1d(data, tweights, 0.16)
    q2 = stat.quantile_1d(data, tweights, 0.84)
    median = stat.median(data, tweights)
    textstr = "$\mathrm{median}:%.2g^{+%.2g}_{-%.2g}$" % (round(median,2), np.abs(median - q2),np.abs(median - q1))

    return textstr


def aveError(depth,data):
    depth = np.round(depth.astype(float), -1)
    means = []
    stds = []
    depths = []
    for d in np.unique(np.round(depth.astype(float), -1)):
        mask = depth==d
        means.append(data[mask].mean())
        stds.append(data[mask].std())
        depths.append(d)
    return np.asarray(means), np.asarray(stds), np.asarray(depths)


# A few depths of interest that marks when the reflection coefficient becomes 0.5 or 0.1 or is TIR
R_50 = 937.63763764
R_10 = 1180.98098098
R_TIR = 918.91891892


color = ['C1','C2','C3','C4','C5','C6','C7','C8']
datafile = PathToARIANNAanalysis + '/data/expectedArrivalDirectionSpice2018smooth.npy'
dataExpected = np.load(datafile,allow_pickle=True)
tmpAzi = []
for i in range(len(dataExpected[0])):
    tmpAzi.append(expectedAzi)
dataExpected = np.asarray([dataExpected[0],dataExpected[1],np.asarray(tmpAzi)])

correctForTilt = True
if correctForTilt:
    import getLaunchAndArrivalAnglesFromSPICEwithNewTiltMeasurements2020 as tilt2020

    zen, azi, depths2 = tilt2020.getOffsetForAngularData()
    fig4, ax4 = plt.subplots(1, 1)
    ax4.plot(depths2, zen,linewidth=3,linestyle='solid',label='zenith')
    ax4.plot(depths2, azi,linewidth=2,linestyle='dotted',label='azimuth')
    ax4.legend()
    ax4.set_xlabel(r'Z [m]')
    ax4.set_ylabel(r'$\Delta$ receive azimuth from no tilt [$^{\circ}$]')
    ax4.set_xlim(980,1700.0)
    fig4.tight_layout()
    fig4.savefig(PathToARIANNAanalysis + '/plots/changeInRecieveAngleFor2020TiltProfile.png')
    fig4.savefig(PathToARIANNAanalysis + '/plots/changeInRecieveAngleFor2020TiltProfile.pdf')
    for i in range(len(dataExpected[0])):
        idx = hu.find_nearest(depths2, dataExpected[0][i])
        dataExpected[1][i] += zen[idx]
        dataExpected[2][i] += azi[idx]

fig, ax = plt.subplots(2, 1,figsize=(11, 7),sharex=True)


fileDipoles = PathToARIANNAanalysis + '/data/reconstructedAngularDirectionFromDipoles1180mDepthAndBelowWithBandpassFilter80to300MHz.npy'
fileLPDAs = PathToARIANNAanalysis + '/data/reconstructedAngularDirectionFromLPDAs938mDepthAndBelowWithBandpassFilter80to300MHz.npy'

deltaZenDipoles, deltaAziDipoles, depth = getDataDiff(fileDipoles,dataExpected)
deltaZenLPDAs, deltaAziLPDAs, depth = getDataDiff(fileLPDAs,dataExpected)

means, errors, depths = aveError(depth,deltaZenDipoles)
means2, errors2, depths2 = aveError(depth,deltaAziDipoles)
means3, errors3, depths3 = aveError(depth,np.asarray(deltaZenLPDAs).astype(float))
means4, errors4, depths4 = aveError(depth,np.asarray(deltaAziLPDAs).astype(float))

size_marker = 40
ax[0].scatter(depth,deltaZenDipoles,s=size_marker,marker='s',color='red',alpha=0.25)
ax[1].scatter(depth,deltaAziDipoles,s=size_marker,marker='s',color='red',alpha=0.25)
ax[0].scatter(depth,deltaZenLPDAs,s=size_marker,marker='>',color='deepskyblue',alpha=0.25)
ax[1].scatter(depth,deltaAziLPDAs,s=size_marker,marker='>',color='deepskyblue',alpha=0.25)

ax[0].scatter(depths,means,s=size_marker,marker='s',color='darkred',label='Dipoles ave.')
ax[1].scatter(depths2,means2,s=size_marker,marker='s',color='darkred',label='Dipoles ave.')
ax[0].scatter(depths3,means3,s=size_marker,marker='>',color='midnightblue',label='LPDAs ave.')
ax[1].scatter(depths4,means4,s=size_marker,marker='>',color='midnightblue',label='LPDAs ave.')

mask1700 = (depths > 1680.0) & (depths < 1710)
mask1000 = (depths > 990.0) & (depths < 1010)
print(errors[mask1700])
print(errors[mask1000])

mask1700 = (depths2 > 1680.0) & (depths2 < 1710)
mask1000 = (depths2 > 990.0) & (depths2 < 1010)
print(errors2[mask1700])
print(errors2[mask1000])

mask1700 = (depths3 > 1680.0) & (depths3 < 1710)
mask1000 = (depths3 > 990.0) & (depths3 < 1010)
print(errors3[mask1700])
print(errors3[mask1000])

mask1700 = (depths4 > 1680.0) & (depths4 < 1710)
mask1000 = (depths4 > 990.0) & (depths4 < 1010)
print(errors4[mask1700])
print(errors4[mask1000])
print(1/0)


mask_space1 = (np.asarray(depth) <= 920.0)
mask_space2 = (np.asarray(depth) >= 920.0) & (np.asarray(depth) <= 1250.0)
mask_space3 = (np.asarray(depths) >= 1250.0) & (np.asarray(depths) <= 1530.0)
mask_space4 = (np.asarray(depth) >= 1530.0)
depth = np.asarray(depth)
depths = np.asarray(depths)
vals = []
vals.append(depth[mask_space1][-1])
vals.append(depth[mask_space2][0])
vals.append(depth[mask_space2][-1])
vals.append(depths[mask_space3][0])
vals.append(depths[mask_space3][-1])
vals.append(depth[mask_space4][0])

ax[0].fill_between([vals[0],vals[1]],[17,17],color='black',alpha=0.15,linewidth=0.0)
ax[0].fill_between([vals[0],vals[1]],[-17,-17],color='black',alpha=0.15,linewidth=0.0)
ax[0].fill_between([vals[2],vals[3]],[17,17],color='black',alpha=0.15,linewidth=0.0)
ax[0].fill_between([vals[2],vals[3]],[-17,-17],color='black',alpha=0.15,linewidth=0.0)
ax[0].fill_between([vals[4],vals[5]],[17,17],color='black',alpha=0.15,linewidth=0.0)
ax[0].fill_between([vals[4],vals[5]],[-17,-17],color='black',alpha=0.15,linewidth=0.0)
ax[1].fill_between([vals[0],vals[1]],[17,17],color='black',alpha=0.15,linewidth=0.0)
ax[1].fill_between([vals[0],vals[1]],[-17,-17],color='black',alpha=0.15,linewidth=0.0)
ax[1].fill_between([vals[2],vals[3]],[17,17],color='black',alpha=0.15,linewidth=0.0)
ax[1].fill_between([vals[2],vals[3]],[-17,-17],color='black',alpha=0.15,linewidth=0.0)
ax[1].fill_between([vals[4],vals[5]],[17,17],color='black',alpha=0.15,linewidth=0.0)
ax[1].fill_between([vals[4],vals[5]],[-17,-17],color='black',alpha=0.15,linewidth=0.0)

ax[0].axhline(y=0.0,color='black')
ax[0].set_ylabel('$\Delta$ zenith angle [$\degree$]')

ax[1].axhline(y=0.0,color='black')
ax[1].set_ylabel(r'$\Delta$ azimuth angle [$\degree$]')
ax[1].set_xlabel('pulser depth [m]')

ax[0].set_ylim(-1.49,1.49)
ax[0].set_xlim(800,1700)
ax[1].set_ylim(-1.49,1.49)
ax[1].set_xlim(800,1700)


ax[0].plot([1180,1180],[-5,5],color='darkred')
ax[1].plot([1180,1180],[-5,5],color='darkred')
ax[0].plot([938,938],[-5,5],color='midnightblue')
ax[1].plot([938,938],[-5,5],color='midnightblue')

ax[0].annotate(s='', xy=(1250.0,1.23), xytext=(1180,1.23),ha='center',va='top',color='darkred',arrowprops=dict(arrowstyle='->',color='darkred'))
ax[0].text(0.53,0.98,s=r'R$\leq$0.1',horizontalalignment='right',verticalalignment='top',transform=ax[0].transAxes,color='darkred')
ax[1].annotate(s='', xy=(1250.0,1.23), xytext=(1180,1.23),ha='center',va='top',color='darkred',arrowprops=dict(arrowstyle='->',color='darkred'))
ax[1].text(0.53,0.98,s=r'R$\leq$0.1',horizontalalignment='right',verticalalignment='top',transform=ax[1].transAxes,color='darkred')

ax[0].annotate(s='', xy=(1008.0,1.23), xytext=(938,1.23),ha='center',va='top',color='midnightblue',arrowprops=dict(arrowstyle='->',color='midnightblue'))
ax[0].text(0.26,0.98,s=r'R$\leq$0.5',horizontalalignment='right',verticalalignment='top',transform=ax[0].transAxes,color='midnightblue')
ax[1].annotate(s='', xy=(1008.0,1.23), xytext=(938,1.23),ha='center',va='top',color='midnightblue',arrowprops=dict(arrowstyle='->',color='midnightblue'))
ax[1].text(0.26,0.98,s=r'R$\leq$0.5',horizontalalignment='right',verticalalignment='top',transform=ax[1].transAxes,color='midnightblue')



fig.tight_layout()
fig.subplots_adjust(wspace = 0, hspace = 0.02)


divider = make_axes_locatable(ax[0])
axHisty0 = divider.append_axes("right", 2.2, pad=0.1,sharey=ax[0])

divider = make_axes_locatable(ax[1])
axHisty1 = divider.append_axes("right", 2.2, pad=0.1,sharey=ax[1],sharex=axHisty0)

mask_dipole = (depth >= 1180.0)
mask_lpda = (depth >= 938.0)

print('stats without depth cuts:')
print('zen: lpdas')
textstr3 = r"$\mu$ = %.2g$^{\circ}$ +- %.2g$^{\circ}$" % (round(np.mean(deltaZenLPDAs),4),round(np.std(deltaZenLPDAs),4))
print(textstr3)
print('zen: dipoles')
textstr3 = r"$\mu$ = %.2g$^{\circ}$ +- %.2g$^{\circ}$" % (round(np.mean(deltaZenDipoles),4),round(np.std(deltaZenDipoles),4))
print(textstr3)
print('azi: lpdas')
textstr3 = r"$\mu$ = %.2g$^{\circ}$ +- %.2g$^{\circ}$" % (round(np.mean(deltaAziLPDAs),4),round(np.std(deltaAziLPDAs),4))
print(textstr3)
print('azi: dipoles')
textstr3 = r"$\mu$ = %.2g$^{\circ}$ +- %.2g$^{\circ}$" % (round(np.mean(deltaAziDipoles),4),round(np.std(deltaAziDipoles),4))
print(textstr3)

gaussianHist(axHisty0,deltaZenLPDAs[mask_lpda],'midnightblue',getMeanSTDStr(deltaZenLPDAs[mask_lpda]),[-4.35,160],'--','lpdas')
gaussianHist(axHisty0,deltaZenDipoles[mask_dipole],'darkred',getMeanSTDStr(deltaZenDipoles[mask_dipole]),[0,165],'-','dipoles')

gaussianHist(axHisty1,deltaAziLPDAs[mask_lpda],'midnightblue',getMeanSTDStr(deltaAziLPDAs[mask_lpda]),[0.2,150],'--','lpdas')
gaussianHist(axHisty1,deltaAziDipoles[mask_dipole],'darkred',getMeanSTDStr(deltaAziDipoles[mask_dipole]),[-4.25,140],'-','dipoles')

ax[0].text(0.99,0.98,r'$\blacksquare$ dipoles',horizontalalignment='right',verticalalignment='top',transform=ax[0].transAxes,color='darkred')
ax[0].text(0.99,0.90,r'$\blacktriangleright$ lpdas',horizontalalignment='right',verticalalignment='top',transform=ax[0].transAxes,color='midnightblue')
ax[0].text(0.99,0.82,'Comm. Period',horizontalalignment='right',verticalalignment='top',transform=ax[0].transAxes,color='grey')

ax[0].xaxis.set_ticks(np.arange(800,1700,100))
ax[1].xaxis.set_ticks(np.arange(800,1700,100))



axHisty0.tick_params(labelleft=False)
axHisty1.tick_params(labelleft=False)
axHisty1.xaxis.set_ticks(np.arange(250,750,250))
axHisty1.set_xlabel('Number of entries')

save = PathToARIANNAanalysis + '/plots/angularSpiceData.pdf'
fig.savefig(save)
save = PathToARIANNAanalysis + '/plots/angularSpiceData.png'
fig.savefig(save)



plt.show()
