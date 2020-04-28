import numpy as np
import matplotlib.pyplot as plt
from matplotlib import mlab
from pylab import *
import helpful_defs as hd
import os

PathToARIANNAanalysis = os.getcwd()

def findZenith(depth): # Uses raytracer solutions from NuRadioMC
# depth (m), launch angle, arival rangle, predicted arrival angle, reflected angle
	datafile = PathToARIANNAanalysis + '/data/AngleVsDepthSP_650m_C.csv'
	data = np.loadtxt(datafile, skiprows=1, delimiter=',').T
	depths = data[0]
	index = hd.find_nearest(depths, depth)
	return [depths[index],data[1][index],data[2][index],data[3][index]]
	# returns depths , Launch angle, arrival zenith angle (bad at shallow depths), predicted arrival zenith angle

def findZenAziFromSpiceDataLPDAs(depth): # Found from cross correlator on SPICE data
# depth (m), launch angle, arival rangle, predicted arrival angle
	data = np.load(PathToARIANNAanalysis + '/data/reconstructedAnglesLPDAs.npy',allow_pickle=True,encoding='bytes')
	depths = data[0]
	index = hd.find_nearest(depths, depth)
	return [depths[index],[],[],data[1][index]]
	# returns depths , [], [], arrival angle [zenith, azimuth]

def plotZenith():
	datafile = PathToARIANNAanalysis + '/data/AngleVsDepthSP_650m_C.csv'
	data = np.loadtxt(datafile, skiprows=1, delimiter=',').T

	fig, ax = plt.subplots(1, 1,figsize=(12, 5))
	ax.plot(data[0],data[2],'o',label='arrival zen.')
	ax.plot(data[0],data[1],'o',label='Launch zen.')
	ax.plot(data[0],data[3],'o',label='arrival zen. calced')
	ax.set_ylabel('Zen')
	plt.axvline(x=400.0,label='shadow_zone',color='black')
	plt.legend()
	ax.set_xlabel('Pulser depth (m)')
	plt.show()


def main():
	plotZenith()

if __name__== "__main__":
	main()
