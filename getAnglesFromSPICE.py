import numpy as np
import matplotlib.pyplot as plt
import helpfulUtils as hu
import os
PathToARIANNAanalysis = os.environ['ARIANNAanalysis']

def findZenith(depth): # Uses raytracer solutions from NuRadioMC
    # the data structure for this input file is: [depth (m), launch angle, arival rangle, predicted arrival angle, reflected angle]
    datafile = PathToARIANNAanalysis + '/data/AngleVsDepthSP_650m_C.csv'
    data = np.loadtxt(datafile, skiprows=1, delimiter=',').T
    depths = data[0]
    index = hu.find_nearest(depths, depth)
    return [depths[index],data[1][index],data[2][index],data[3][index]]

def findZenAziFromSpiceDataLPDAs(depth): # Found from cross correlator on SPICE data
    # the data structure for this input file is: [depth (m), launch angle, arival rangle, predicted arrival angle]
    data = np.load(PathToARIANNAanalysis + '/data/reconstructedAnglesLPDAs.npy',allow_pickle=True,encoding='bytes')
    depths = data[0]
    index = hu.find_nearest(depths, depth)
    return [depths[index],[],[],data[1][index]]

def plotZenith():
    datafile = PathToARIANNAanalysis + '/data/AngleVsDepthSP_650m_C.csv'
    data = np.loadtxt(datafile, skiprows=1, delimiter=',').T

    fig, ax = plt.subplots(1, 1,figsize=(12, 5))
    ax.plot(data[0],data[2],'o',label='arrival zen.')
    ax.plot(data[0],data[1],'o',label='Launch zen.')
    ax.plot(data[0],data[3],'o',label='arrival zen. calced')
    ax.set_ylabel(r'RF zenith [$\degree$]')
    plt.axvline(x=400.0,label='shadow zone',color='black')
    plt.legend()
    ax.set_xlabel('Pulser depth (m)')

    save = PathToARIANNAanalysis + '/plots/arrivalAndLaunchAnglesForSpiceData.png'
    fig.savefig(save)
    # pdf versions of the figures are nice for LaTeX documents as they preserve resolution better than png files when zooming
    save = PathToARIANNAanalysis + '/plots/arrivalAndLaunchAnglesForSpiceData.pdf'
    fig.savefig(save)

    plt.show()


def main():
    plotZenith()

if __name__== "__main__":
    main()
