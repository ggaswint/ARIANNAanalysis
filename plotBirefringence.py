import numpy as np
import matplotlib.pyplot as plt
import os

PathToARIANNAanalysis = os.environ['ARIANNAanalysis']

def aveError(depth,data):
	depth = np.asarray(depth)
	depth = np.round(depth.astype(float), -1)
	means = []
	stds = []
	depths = []
	for d in np.unique(depth):
		mask = depth==d
		n = len(data[mask])
		means.append(data[mask].mean())
		stds.append(data[mask].std()/np.sqrt(n-1))
		depths.append(d)
	return means, np.asarray(stds), depths

def plotBirefringence(file_LPDA, file_Dipole):
	fig, ax = plt.subplots(1, 1,figsize=(8, 5)) # test
	data_LPDA = np.load(file_LPDA,allow_pickle=True)
	data_Dipole = np.load(file_Dipole,allow_pickle=True)
	means_LPDA, stds_LPDA, depths = aveError(data_LPDA[0],data_LPDA[4])
	means_Dipole, stds_Dipole, depths = aveError(data_Dipole[0],data_Dipole[4])
	ax.scatter(depths,np.asarray(means_Dipole)-np.asarray(means_LPDA),s=50,alpha=0.75,color='midnightblue',label='SPICE')
	ax.set_ylabel('time of (vpol - hpol) [ns]')
	ax.set_xlabel('depths [m]')

	# NOTE: theory data is defined as hpol (LPDA) - vpol (Dipole)
	data_theory = np.loadtxt(PathToARIANNAanalysis + '/data/BirefringenceTheory.txt', delimiter=' ')
	data_theory = np.asarray(data_theory).T
	ax.scatter(data_theory[0],-data_theory[1],s=50,alpha=0.75,color='orange',label='Theory')
	ax.set_xlim(800,1700)
	ax.legend()

	fig.tight_layout()
	fig.savefig(PathToARIANNAanalysis + '/plots/Birefringence.png')
	fig.savefig(PathToARIANNAanalysis + '/plots/Birefringence.pdf')
	plt.show()

def main():
	file_LPDA = PathToARIANNAanalysis + '/data/SPICE_Dec30_2018_EFieldDataLPDA.npy'
	file_Dipole = PathToARIANNAanalysis + '/data/SPICE_Dec30_2018_EFieldDataDipole.npy'

	plotBirefringence(file_LPDA, file_Dipole)

if __name__== "__main__":
	main()
