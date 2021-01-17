import numpy as np
from NuRadioReco.utilities import units

def mergeArrays(arr1, arr2, arr4, arr5):
	n1 = len(arr1)
	n2 = len(arr2)
	arr3 = [None] * (n1 + n2)
	arr6 = [None] * (n1 + n2)
	i = 0
	j = 0
	k = 0

	# Traverse both array
	while i < n1 and j < n2:

		# Check if current element
		# of first array is smaller
		# than current element of
		# second array. If yes,
		# store first array element
		# and increment first array
		# index. Otherwise do same
		# with second array
		if arr1[i] < arr2[j]:
			arr3[k] = arr1[i]
			arr6[k] = arr4[i]
			k = k + 1
			i = i + 1
		else:
			arr3[k] = arr2[j]
			arr6[k] = arr5[j]
			k = k + 1
			j = j + 1


	# Store remaining elements
	# of first array
	while i < n1:
		arr3[k] = arr1[i]
		arr6[k] = arr4[i]
		k = k + 1
		i = i + 1

	# Store remaining elements
	# of second array
	while j < n2:
		arr3[k] = arr2[j]
		arr6[k] = arr5[j]
		k = k + 1
		j = j + 1
	#print("Array after merging")
	#for i in range(n1 + n2):
	#    print(str(arr3[i]), end = " ")
	return arr3, arr6

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

def bandPassFilter_EField(evt, station, det, passband=[55 * units.MHz, 1000 * units.MHz],
		filter_type='rectangular',):

	channels = station.get_channels()
	for channel in channels:
		frequencies = channel.get_electric_field().get_frequencies()
		trace_fft = channel.get_electric_field().get_frequency_spectrum()

		if(filter_type == 'rectangular'):
			trace_fft[np.where(frequencies < passband[0])] = 0.
			trace_fft[np.where(frequencies > passband[1])] = 0.
		elif(filter_type == 'butter10'):
			b, a = scipy.signal.butter(10, passband, 'bandpass', analog=True)
			w, h = scipy.signal.freqs(b, a, frequencies)
			trace_fft *= h
		elif(filter_type == 'butter10abs'):
			b, a = scipy.signal.butter(10, passband, 'bandpass', analog=True)
			w, h = scipy.signal.freqs(b, a, frequencies)
			trace_fft *= np.abs(h)
		channel.get_electric_field().set_frequency_spectrum(trace_fft, channel.get_sampling_rate())

def chi2fit(xpoints, ypoints):
	xpoints = np.array(xpoints)
	ypoints = np.array(ypoints)

	# define some intermediate terms (analytic solution to
	# neyman's chisquare definition, assuming sigma=sqrt(N))
	A = np.sum(xpoints/ypoints)
	B = np.sum(1/ypoints)
	C = len(ypoints)
	D = np.sum(xpoints**2/ypoints)
	E = np.sum(xpoints)

	# point estimates
	slope = (E*B - C*A) / (D*B - A**2)
	offset = (D*C - E*A) / (D*B - A**2)
	Vaa = B/(D*B - A**2)
	Vbb = D/(D*B - A**2)
	Vab = -A/(D*B - A**2)

	yfit = np.poly1d([slope, offset])
	yerrs = lambda x: np.sqrt(Vbb + (np.asarray(x)**2)*Vaa + 2*np.asarray(x)*Vab)

	return [slope,offset], yfit, yerrs


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
