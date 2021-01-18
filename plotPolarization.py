import os
import numpy as np
import matplotlib.pyplot  as plt
from scipy.interpolate import interp1d

PathToARIANNAanalysis = os.environ['ARIANNAanalysis']

def gaussianHist(ax,data,color,label,line):
    n, bins, patches = ax.hist(data,linestyle=line,bins=np.arange(-15.0, 15.1, 1.0),edgecolor=color,fill=False,histtype='step')#,weights=weights) histype = bar
    ax.text(1.0,0.99,label,horizontalalignment='right',verticalalignment='top',transform=ax.transAxes,color='midnightblue')

def getMeanSTDStr(data):
    mean = np.mean(data)
    std = np.std(data)
    textstr = r"mean = %.2g$^{\circ}$, STD = %.2g$^{\circ}$" % (round(mean,2),round(std,2))
    return textstr

def aveError(depth,data):
    depth = np.asarray(depth)
    depth = np.round(depth.astype(float), -1)
    means = []
    stds = []
    depths = []
    for d in np.unique(np.round(depth.astype(float), -1)):
        mask = depth==d
        means.append(data[mask].mean())
        stds.append(data[mask].std())
        depths.append(d)
    return means, stds, depths

def plotAnechoic(ax):
    file = PathToARIANNAanalysis + '/data/anechoicPolarizationPrediction.npy' # 1,2,4,6
    data = np.load(file,allow_pickle=True)
    # Data ===> Depths, PolReco + sigma, PolReco - sigma, PolReco
    ax.fill_between(data[0], data[1], data[2],color='darkorange',label='Anechoic no Atten', alpha=.25)
    return data[0], data[3]


def plotPolarization():
    fig, ax = plt.subplots(1, 1,figsize=(11, 7))
    fig2, ax2 = plt.subplots(1, 1,figsize=(5,5))

    file = PathToARIANNAanalysis + '/data/SPICE_Dec30_2018_EFieldDataLPDA_AtTransmitter.npy'
    data = np.load(file,allow_pickle=True)
    # data ===> depths, polarization, max amplitudes theta, max amplitudes phi,time at max theta, time at max phi
    data0 = np.delete(data[0],720) # A bad event 721,722,725,2002,2143,2367
    data1 = np.delete(data[1],720)


    means, errors, depths = aveError(data0,np.rad2deg(np.arctan(data1)))
    depths = np.asarray(depths)
    means = np.asarray(means)
    errors = np.asarray(errors)

    # Find gray shading regions, where no data lies.
    shading_inputs = []
    mask_space1 = (depths <= 920.0)
    mask_space2 = (depths >= 920.0) & (depths <= 1250.0)
    mask_space3 = (depths >= 1250.0) & (depths <= 1530.0)
    mask_space4 = (depths >= 1530.0)
    shading_inputs.append(depths[mask_space1][-1])
    shading_inputs.append(depths[mask_space2][0])
    shading_inputs.append(depths[mask_space2][-1])
    shading_inputs.append(depths[mask_space3][0])
    shading_inputs.append(depths[mask_space3][-1])
    shading_inputs.append(depths[mask_space4][0])

    # Plot data
    ax.errorbar(depths,means,errors,fmt='d',MarkerSize=8.0,color='midnightblue')

    file = PathToARIANNAanalysis + '/data/SPICE_Dec30_2018_EFieldDataLPDA_AtTransmitter_SystematicAntennaStudy.npy'
    SPICE_data_sysStudy = np.load(file,allow_pickle=True)
    # SPICE_data_sysStudy ===> depths, polData + 68%, polData - 68%
    mask_space1 = (np.asarray(SPICE_data_sysStudy[0]) <= 920.0)
    mask_space2 = (np.asarray(SPICE_data_sysStudy[0]) >= 920.0) & (np.asarray(SPICE_data_sysStudy[0]) <= 1250.0)
    mask_space3 = (np.asarray(SPICE_data_sysStudy[0]) >= 1250.0) & (np.asarray(SPICE_data_sysStudy[0]) <= 1530.0)
    mask_space4 = (np.asarray(SPICE_data_sysStudy[0]) >= 1530.0)
    ax.fill_between(SPICE_data_sysStudy[0][mask_space1], SPICE_data_sysStudy[1][mask_space1], SPICE_data_sysStudy[2][mask_space1],color='deepskyblue', alpha=.25)
    ax.fill_between(SPICE_data_sysStudy[0][mask_space2], SPICE_data_sysStudy[1][mask_space2], SPICE_data_sysStudy[2][mask_space2],color='deepskyblue', alpha=.25)
    ax.fill_between(SPICE_data_sysStudy[0][mask_space3], SPICE_data_sysStudy[1][mask_space3], SPICE_data_sysStudy[2][mask_space3],color='deepskyblue', alpha=.25)
    ax.fill_between(SPICE_data_sysStudy[0][mask_space4], SPICE_data_sysStudy[1][mask_space4], SPICE_data_sysStudy[2][mask_space4],color='deepskyblue', alpha=.25)

    ax.set_xlabel('pulser depth [m]')
    ax.set_ylabel(r'polarization $[^\circ]$')

    #ax.fill_between([800,900], 0, 30, color='grey',alpha=0.5)
    ax.set_ylim(0,30)
    ax.set_xlim(800,1700)

    ax.text(0.99,0.98,'SPICE',horizontalalignment='right',verticalalignment='top',transform=ax.transAxes,color='midnightblue')
    ax.text(0.99,0.94,'Anechoic',horizontalalignment='right',verticalalignment='top',transform=ax.transAxes,color='darkorange')
    ax.text(0.99,0.90,'Comm. Period',horizontalalignment='right',verticalalignment='top',transform=ax.transAxes,color='grey')
    ax.annotate(s='', xy=(1008.0,28.2), xytext=(938,28.2),ha='center',va='top',color='midnightblue',arrowprops=dict(arrowstyle='->',color='midnightblue'))
    ax.text(0.27,0.98,s=r'R$\leq$0.5',horizontalalignment='right',verticalalignment='top',transform=ax.transAxes,color='midnightblue')
    ax.plot([938,938],[-5,35],color='midnightblue')

    ax.fill_between([shading_inputs[0],shading_inputs[1]],[60,60],color='black',alpha=0.15,linewidth=0.0)
    ax.fill_between([shading_inputs[0],shading_inputs[1]],[-17,-17],color='black',alpha=0.15,linewidth=0.0)
    ax.fill_between([shading_inputs[2],shading_inputs[3]],[60,60],color='black',alpha=0.15,linewidth=0.0)
    ax.fill_between([shading_inputs[2],shading_inputs[3]],[-17,-17],color='black',alpha=0.15,linewidth=0.0)
    ax.fill_between([shading_inputs[4],shading_inputs[5]],[60,60],color='black',alpha=0.15,linewidth=0.0)
    ax.fill_between([shading_inputs[4],shading_inputs[5]],[-17,-17],color='black',alpha=0.15,linewidth=0.0)


    mask_spice = (data0 >= 938.0)
    spice_depths = data0[mask_spice]
    spice_data = np.rad2deg(np.arctan(np.asarray(data1)))[mask_spice]

    AnechoicDepths, AnechoicPolarizations = plotAnechoic(ax)

    func = interp1d(AnechoicDepths,AnechoicPolarizations)
    anechoic_data = func(spice_depths)

    gaussianHist(ax2,spice_data-anechoic_data,'midnightblue',getMeanSTDStr(spice_data-anechoic_data),'-')
    ax2.set_xlabel(r'polarization [$\degree$]')
    ax2.set_ylabel(r'Number of events')

    fig.tight_layout()
    fig.savefig(PathToARIANNAanalysis + '/plots/polReco.png')
    fig.savefig(PathToARIANNAanalysis + '/plots/polReco.pdf')

    fig2.tight_layout()
    fig2.savefig(PathToARIANNAanalysis + '/plots/polRecoHist.png')
    fig2.savefig(PathToARIANNAanalysis + '/plots/polRecoHist.pdf')


def main():
    plotPolarization()
    plt.show()

if __name__== "__main__":
    main()
