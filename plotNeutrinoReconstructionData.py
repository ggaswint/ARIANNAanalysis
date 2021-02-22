import numpy as np
from NuRadioReco.utilities import units, fft, trace_utilities
import matplotlib.pyplot as plt
import h5py
from radiotools import stats as stat
from scipy.stats import norm
import matplotlib.mlab as mlab
from scipy.optimize import curve_fit
import pylab as py
from radiotools import plthelpers as ph
from radiotools import helper as hp
import numpy.polynomial.polynomial as poly
from colour import Color
from radiotools import helper as hp
import os
PathToARIANNAanalysis = os.environ['ARIANNAanalysis']

cherenkov_angle = np.rad2deg(np.arccos(1./1.78))

def gaussianHist(ax,data,bins,color,line,density,weights,useWeights):
    if useWeights:
        n, bins, patches = ax.hist(data,linestyle=line,bins=bins,edgecolor=color,fill=False,histtype='step',density=density,weights=weights)
    else:
        n, bins, patches = ax.hist(data,linestyle=line,bins=bins,edgecolor=color,fill=False,histtype='step',density=density)
    return n, bins, patches

def get_single_angle(zenith_reco,azimuth_reco,zenith_exp,azimuth_exp):
    vspherical_to_cartesian = np.vectorize(hp.spherical_to_cartesian,otypes=[np.ndarray])
    vget_angle = np.vectorize(hp.get_angle)
    v1 = vspherical_to_cartesian(zenith_reco, azimuth_reco)
    v2 = vspherical_to_cartesian(zenith_exp, azimuth_exp)
    return vget_angle(v1, v2)

def parseData(data):
    SNRs = []
    weights = []
    triggers = []
    nu_zenith_sim = []
    nu_zenith_reco = []
    nu_azimuth_sim = []
    nu_azimuth_reco = []
    shower_energy_sim = []
    shower_energy_reco = []

    # Turn this one when using multiple reconstructed events
    # Assuming that multiple events have been reconstructed and the resulting files concatenated together to make arrays for each of these keys
    if 0:
       SNRs.append(data.item().get("SNRs")[0])
       weights.append(data.item().get("weight")[0])
       triggers.append(data.item().get("triggers")[0])
       nu_zenith_sim.append(data.item().get("nu_zenith_sim")[0])
       nu_zenith_reco.append(data.item().get("nu_zenith_reco")[0])
       nu_azimuth_sim.append(data.item().get("nu_azimuth_sim")[0])
       nu_azimuth_reco.append(data.item().get("nu_azimuth_reco")[0])
       shower_energy_sim.append(data.item().get("shower_energy_sim")[0])
       shower_energy_reco.append(data.item().get("shower_energy_reco")[0])

    # This is for a single event, done for ARIANNAanalysis repo example, but the plots will be weak because only 1 data point in the example file
    if 1:
       SNRs.append(data.item().get("SNRs"))
       weights.append(data.item().get("weight"))
       triggers.append(data.item().get("triggers"))
       nu_zenith_sim.append(data.item().get("nu_zenith_sim"))
       nu_zenith_reco.append(data.item().get("nu_zenith_reco"))
       nu_azimuth_sim.append(data.item().get("nu_azimuth_sim"))
       nu_azimuth_reco.append(data.item().get("nu_azimuth_reco"))
       shower_energy_sim.append(data.item().get("shower_energy_sim"))
       shower_energy_reco.append(data.item().get("shower_energy_reco"))

    SNRs = np.asarray(SNRs)
    weights = np.asarray(weights)
    triggers = np.asarray(triggers)
    nu_zenith_sim = np.asarray(nu_zenith_sim)
    nu_zenith_reco = np.asarray(nu_zenith_reco)
    nu_azimuth_sim = np.asarray(nu_azimuth_sim)
    nu_azimuth_reco = np.asarray(nu_azimuth_reco)
    shower_energy_sim = np.asarray(shower_energy_sim)
    shower_energy_reco = np.asarray(shower_energy_reco)

    return SNRs, weights, triggers, nu_zenith_sim ,nu_zenith_reco ,nu_azimuth_sim, nu_azimuth_reco, shower_energy_sim, shower_energy_reco

def parseSNRs(SNRs):
    SNRsMax = []
    SNRsMaxLPDA = []
    SNRsMinLPDA = []
    SNRsCh0 = []
    SNRsCh1 = []
    SNRsCh2 = []
    SNRsCh3 = []
    SNRsCh4 = []

    for i in range(len(SNRs)):
        SNRsMax.append(max(SNRs[i]))
        SNRsMaxLPDA.append(max(SNRs[i][0:4]))
        SNRsMinLPDA.append(min(SNRs[i][0:4]))
        SNRsCh0.append(SNRs[i][0])
        SNRsCh1.append(SNRs[i][1])
        SNRsCh2.append(SNRs[i][2])
        SNRsCh3.append(SNRs[i][3])
        SNRsCh4.append(SNRs[i][4])

    SNRsMax = np.asarray(SNRsMax)
    SNRsMaxLPDA = np.asarray(SNRsMaxLPDA)
    SNRsMinLPDA = np.asarray(SNRsMinLPDA)
    SNRsCh0 = np.asarray(SNRsCh0)
    SNRsCh1 = np.asarray(SNRsCh1)
    SNRsCh2 = np.asarray(SNRsCh2)
    SNRsCh3 = np.asarray(SNRsCh3)
    SNRsCh4 = np.asarray(SNRsCh4)

    return SNRsMax, SNRsMaxLPDA, SNRsMinLPDA, SNRsCh0, SNRsCh1, SNRsCh2, SNRsCh3, SNRsCh4

def getMeanSTDStrRaylieghWeightedHist(data,weights):
    textstr = r"$\sigma_{68\%%}$=%.2g$^{\circ}$" % (stat.quantile_1d(data,weights,0.68))
    return textstr

def getMeanSTDStrWeighted(data, weights):
    data = data.astype('float')
    tweights = np.ones_like(data)
    median = stat.median(data, tweights)
    textstr = "$\mathrm{median} = %.2g^{+%.2g}_{-%.2g}$" % (median, (stat.quantile_1d(data,weights,0.84)-median),(median-stat.quantile_1d(data,weights,0.16)))
    return textstr

def populatePlots(inputFile,color):
    data = np.load(inputFile,allow_pickle=True)
    SNRs, weights, triggers, nu_zenith_sim ,nu_zenith_reco ,nu_azimuth_sim, nu_azimuth_reco, shower_energy_sim, shower_energy_reco  = parseData(data)

    # example to check if trigger type was true
    #if triggers["LPDA_2of4_100Hz"].has_triggered():
    #   print("do something")

    SNRsMax, SNRsMaxLPDA, SNRsMinLPDA, SNRsCh0, SNRsCh1, SNRsCh2, SNRsCh3, SNRsCh4 = parseSNRs(SNRs)

    mask = SNRsMinLPDA >= 0.0 # example mask can adjust, currently has no effect becasue SNR is always >= 0 unless it does not exist at all

    print("weighted percentage: " + str(round(100.0*np.sum(weights[mask])/np.sum(weights),1)))

    label = getMeanSTDStrWeighted(np.log10(shower_energy_reco[mask]/shower_energy_sim[mask]), weights[mask])
    n, bins, patches = gaussianHist(ax1,np.log10(shower_energy_reco[mask]/shower_energy_sim[mask]),np.arange(-1.0, 1.0, 0.01),color,'-',True,weights[mask],True)
    ax1.text(1.0,0.99,label,horizontalalignment='right',verticalalignment='top',fontsize='small',transform=ax1.transAxes,color=color)

    n, bins, patches = gaussianHist(ax2,nu_zenith_sim[mask]-nu_zenith_reco[mask],np.arange(-20., 20.0, 0.5),color,'-',False,weights[mask],True)
    label = getMeanSTDStrWeighted(nu_zenith_sim[mask]-nu_zenith_reco[mask],weights[mask])
    ax2.text(1.0,0.99,label,horizontalalignment='right',verticalalignment='top',transform=ax2.transAxes,color=color)

    n, bins, patches = gaussianHist(ax3,nu_azimuth_sim[mask]-nu_azimuth_reco[mask],np.arange(-20., 20.0, 0.5),color,'-',False,weights[mask],True)
    label = getMeanSTDStrWeighted(nu_azimuth_sim[mask]-nu_azimuth_reco[mask],weights[mask])
    ax3.text(1.0,0.99,label,horizontalalignment='right',verticalalignment='top',transform=ax3.transAxes,color=color)

    #ax30.text(1.0,0.99 - (0)*0.07,label,horizontalalignment='right',verticalalignment='top',transform=ax30.transAxes,color='black')
    angles = get_single_angle(nu_zenith_reco[mask]*units.deg,nu_azimuth_reco[mask]*units.deg,nu_zenith_sim[mask]*units.deg,nu_azimuth_sim[mask]*units.deg)
    n, bins, patches = gaussianHist(ax4,angles/units.deg,np.arange(0.0, 15.1, 0.5),color,'-',False,weights[mask],True)
    label = getMeanSTDStrRaylieghWeightedHist(angles/units.deg, weights[mask])
    ax4.text(1.0,0.99,label,horizontalalignment='right',verticalalignment='top',transform=ax4.transAxes,color=color)



if __name__ == "__main__":
    fig1, ax1 = plt.subplots(1, 1, sharex=False, figsize=(5, 5))
    fig2, ax2 = plt.subplots(1, 1, sharex=False, figsize=(5, 5))
    fig3, ax3 = plt.subplots(1, 1, sharex=False, figsize=(5, 5))
    fig4, ax4 = plt.subplots(1, 1, sharex=False, figsize=(5, 5))

    extraSave = 'evt0'
    dataFile = PathToARIANNAanalysis + "/data/reconstructedNeutrinoProperties_0.npy"
    populatePlots(dataFile,"blue")

    ax1.set_xlabel(r'$log_{10}(E_{rec}/E_{true})$')
    ax1.set_ylabel(r'normalized triggered events')
    fig1.tight_layout()
    fig1.savefig(PathToARIANNAanalysis + '/plots/Energy_'+extraSave+'.png')

    ax2.set_xlabel(r'$\Delta$ nu zenith (sim - reco) [$\degree$]')
    ax2.set_ylabel('num events')
    fig2.tight_layout()
    fig2.savefig(PathToARIANNAanalysis + '/plots/ZenHist_'+extraSave+'.png')

    ax3.set_xlabel(r'$\Delta$ nu azimuth (sim - reco) [$\degree$]')
    ax3.set_ylabel('num events')
    fig3.tight_layout()
    fig3.savefig(PathToARIANNAanalysis + '/plots/AziHist_'+extraSave+'.png')

    ax4.set_xlabel(r'$\Delta$$\psi_{\nu}$ [$\degree$]')
    ax4.set_ylabel(r'triggered events')
    fig4.tight_layout()
    fig4.savefig(PathToARIANNAanalysis + '/plots/SpaceAngle_'+extraSave+'.png')
