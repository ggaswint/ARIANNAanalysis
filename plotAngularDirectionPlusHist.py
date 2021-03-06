import matplotlib.pyplot  as plt
import numpy as np

import helpfulUtils as hu
from mpl_toolkits.axes_grid1 import make_axes_locatable
from scipy.stats import norm
import os
PathToARIANNAanalysis = os.environ['ARIANNAanalysis']



cut_val = 1.49 # events larger than this uncertainty are cut for plotting purposes (outside plot limits)


def getStartStop(depth, start_depth, end_depth, reverse=False):
    first = True
    second = True
    start=0
    end=0
    if reverse:
        depth = depth[::-1]
    for i in range(len(depth)):
        if depth[i] < end_depth and first:
            start = i
            first = False
        if depth[i] < start_depth and second:
            end = i
            second = False
    return start,end

def angle_data(data):
    Angle = np.asarray(data)# L1 data data
    Zen = []
    for i in range(len(Angle)):
        Zen.append(float(Angle[i][0]))
    Azi = []
    for i in range(len(Angle)):
        Azi.append(float(Angle[i][1]))
    return Zen, Azi


def getDataDiff(file, start_depth, stop_depth, full, dataE, ys, reverse=False):
    set1 = np.load(file,allow_pickle=True,encoding='bytes')
    depth = set1[0]
    data = set1[1]
    start, end = getStartStop(depth,start_depth,stop_depth,reverse)
    if full:
        start = 0
        end = len(depth)
    Zen, Azi = angle_data(data)
    ZenE = []
    for i in range(len(depth)):
        idx = hu.find_nearest(np.asarray(dataE[0]), depth[i])
        ZenE.append(ys[idx])
    gg_exp_diff = (np.asarray(Zen) - np.asarray(ZenE)).astype(float)
    gg_exp_diff2 = (np.asarray(Azi) - 312.448284).astype(float)
    return gg_exp_diff,gg_exp_diff2,start,end,depth

def gaussianHist(ax,gg_exp_diff,color,label,pos,line,label2):
    gg_exp_diff2 = []
    for j in range(len(gg_exp_diff)):
        if np.abs(gg_exp_diff[j]) < cut_val:
            gg_exp_diff2.append(gg_exp_diff[j])

    (mu, sigma) = norm.fit(gg_exp_diff)

    fig2, ax2 = plt.subplots(1, 1)
    n, bins, patches = ax.hist(gg_exp_diff,linestyle=line,bins=np.arange(-1.5, 1.6, 0.2),edgecolor=color,fill=False,label=label2,histtype='step',orientation='horizontal')#,weights=weights) histype = bar

    if label2 == 'dipoles':
        ax.text(0.99,0.98,label,horizontalalignment='right',verticalalignment='top',transform=ax.transAxes,color='darkred',fontsize=11)

    if label2 == 'lpdas':
        ax.text(0.99,0.90,label,horizontalalignment='right',verticalalignment='top',transform=ax.transAxes,color='midnightblue',fontsize=11)

    plt.close(fig2)

def getMeanSTDStr(data):
    gg_exp_diff2 = []
    for j in range(len(data)):
        if np.abs(data[j]) < cut_val:
            gg_exp_diff2.append(data[j])
    data = gg_exp_diff2
    mean = np.mean(data)
    std = np.std(data)
    textstr1 = r"$\mu$ = %.2g" % (round(mean,2))
    textstr2 = r"$\sigma$ = %.2g" % (round(std,2))
    #return [textstr1, textstr2]
    textstr3 = r"mean = %.2g$^{\circ}$, STD = %.2g$^{\circ}$" % (round(mean,2),round(std,2))
    return textstr3


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


R_50 = 937.63763764
R_10 = 1180.98098098
R_TIR = 918.91891892


color = ['C1','C2','C3','C4','C5','C6','C7','C8']
max_diff = 10.0
datafile = PathToARIANNAanalysis + '/data/expectedArrivalDirectionSpice2018smooth.npy'
dataE = np.load(datafile,allow_pickle=True)
ys = np.asarray(dataE[1])


fig, ax = plt.subplots(2, 1,figsize=(11, 7),sharex=True)


file = PathToARIANNAanalysis + '/data/reconstructedAngularDirectionFromDipoles1180mDepthAndBelowWithBandpassFilter80to300MHz.npy'
file2 = PathToARIANNAanalysis + '/data/reconstructedAngularDirectionFromLPDAs938mDepthAndBelowWithBandpassFilter80to300MHz.npy'

gg_exp_diff, gg_exp_diff2, start, end, depth = getDataDiff(file,500,1400,False,dataE,ys)
gg_exp_diff3, gg_exp_diff4, start, end, depth = getDataDiff(file2,1000,1800,False,dataE,ys)

means, errors, depths = aveError(depth,gg_exp_diff)
means2, errors2, depths2 = aveError(depth,gg_exp_diff2)
means3, errors3, depths3 = aveError(depth,np.asarray(gg_exp_diff3).astype(float))
means4, errors4, depths4 = aveError(depth,np.asarray(gg_exp_diff4).astype(float))

size_marker = 40
ax[0].scatter(depth,gg_exp_diff,s=size_marker,marker='s',color='red',alpha=0.25)
ax[1].scatter(depth,gg_exp_diff2,s=size_marker,marker='s',color='red',alpha=0.25)
ax[0].scatter(depth,gg_exp_diff3,s=size_marker,marker='>',color='deepskyblue',alpha=0.25)
ax[1].scatter(depth,gg_exp_diff4,s=size_marker,marker='>',color='deepskyblue',alpha=0.25)

ax[0].scatter(depths,means,s=size_marker,marker='s',color='darkred',label='Dipoles ave.')
ax[1].scatter(depths2,means2,s=size_marker,marker='s',color='darkred',label='Dipoles ave.')
ax[0].scatter(depths3,means3,s=size_marker,marker='>',color='midnightblue',label='LPDAs ave.')
ax[1].scatter(depths4,means4,s=size_marker,marker='>',color='midnightblue',label='LPDAs ave.')

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
textstr3 = r"$\mu$ = %.2g$^{\circ}$ +- %.2g$^{\circ}$" % (round(np.mean(gg_exp_diff3),4),round(np.std(gg_exp_diff3),4))
print(textstr3)
print('zen: dipoles')
textstr3 = r"$\mu$ = %.2g$^{\circ}$ +- %.2g$^{\circ}$" % (round(np.mean(gg_exp_diff),4),round(np.std(gg_exp_diff),4))
print(textstr3)
print('azi: lpdas')
textstr3 = r"$\mu$ = %.2g$^{\circ}$ +- %.2g$^{\circ}$" % (round(np.mean(gg_exp_diff4),4),round(np.std(gg_exp_diff4),4))
print(textstr3)
print('azi: dipoles')
textstr3 = r"$\mu$ = %.2g$^{\circ}$ +- %.2g$^{\circ}$" % (round(np.mean(gg_exp_diff2),4),round(np.std(gg_exp_diff2),4))
print(textstr3)

gaussianHist(axHisty0,gg_exp_diff3[mask_lpda],'midnightblue',getMeanSTDStr(gg_exp_diff3[mask_lpda]),[-4.35,160],'--','lpdas')
gaussianHist(axHisty0,gg_exp_diff[mask_dipole],'darkred',getMeanSTDStr(gg_exp_diff[mask_dipole]),[0,165],'-','dipoles')

gaussianHist(axHisty1,gg_exp_diff4[mask_lpda],'midnightblue',getMeanSTDStr(gg_exp_diff4[mask_lpda]),[0.2,150],'--','lpdas')
gaussianHist(axHisty1,gg_exp_diff2[mask_dipole],'darkred',getMeanSTDStr(gg_exp_diff2[mask_dipole]),[-4.25,140],'-','dipoles')

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
