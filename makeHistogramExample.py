import numpy as np
import matplotlib.pyplot as plt
from radiotools import stats as stat

def gaussianHist(ax,data,bins,color='blue',line='-',label='',weights=[],useWeights=False,density=False):
	#(mu, sigma) = norm.fit(gg_exp_diff)
	if useWeights:
		n, bins, patches = ax.hist(data,linestyle=line,bins=bins,edgecolor=color,fill=False,label=label,histtype='step',density=density,weights=weights)#density=True,weights=weights) histype = bar
	else:
		n, bins, patches = ax.hist(data,linestyle=line,bins=bins,edgecolor=color,fill=False,label=label,histtype='step',density=density)
	return n, bins, patches

def getStatsStrRaylieghWeighted(data,weights): # Use for one sided distirbutions
	data = np.asarray(data)
	textstr = r"$\sigma_{68\%%}$=%.2g$^{\circ}$" % (stat.quantile_1d(data,weights,0.68))
	return textstr

def getStatsStrWeighted(data, weights):
	data = np.asarray(data)
	median = stat.median(data, tweights)
	textstr = "$\mathrm{median} = %.2g^{+%.2g}_{-%.2g}$" % (median, (stat.quantile_1d(data,weights,0.84)-median),(median-stat.quantile_1d(data,weights,0.16)))
	return textstr

if __name__ == "__main__":
    bins = np.arange(-10.0, 10.1, 1.0)
    data = np.random.normal(0, 1.0, 1000)
    tweights = np.ones_like(data)
    fig, ax = plt.subplots(1, 1, sharex=False, figsize=(5, 5))
    n, bins, patches = gaussianHist(ax,data,bins,color='black',line='-',label=r'Gaussian Norm',weights=tweights,useWeights=True,density=False)
    label = getStatsStrWeighted(data,tweights)
    ax.text(1.0,0.99,label,horizontalalignment='right',verticalalignment='top',transform=ax.transAxes,color='black')

    plt.show()
