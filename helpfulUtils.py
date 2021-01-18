import numpy as np
from NuRadioReco.utilities import units

def CalcL1(a):
    ct = np.array(a)**2
    l1 = np.max(ct)/(np.sum(ct)-np.max(ct))
    return l1

def FindL1(fx, len_trace):
	b = []
	for i in range(len(fx)):
		b.append(fx[i]/(np.sqrt(len_trace)))
	return CalcL1(b)

# used in getAnglesFromSPICE.py
def find_nearest(array,value):
	nparr = np.asarray(array)
	return (np.abs(nparr-value)).argmin()

def fit(xs,ys,n=1):
	p, V = np.polyfit(xs, ys, 1, cov=True)
	slope = p[0]
	intercept = p[1]
	error = np.sqrt(V[0][0])
	abline_values = [slope * k + intercept for k in xs]
	label = format(slope, '.1e') + " + " + format(error, '.1e')
	fit_data = [slope,error]
	#n = 1
	t = np.linspace(xs[0], xs[len(xs)-1], len(xs))
	# Matrix with rows 1, t, t**2, ...:
	TT = np.vstack([t**(n-i) for i in range(n+1)]).T
	yi = np.dot(TT, p)  # matrix multiplication calculates the polynomial values
	C_yi = np.dot(TT, np.dot(V, TT.T)) # C_y = TT*C_z*TT.T
	sig_yi = np.sqrt(np.diag(C_yi))  # Standard deviations are sqrt of diagonal
	# Do the plotting:
	return abline_values,label,fit_data,t,yi+sig_yi,yi-sig_yi
