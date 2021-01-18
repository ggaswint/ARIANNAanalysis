import numpy as np
from NuRadioReco.utilities import units
import matplotlib.pyplot as plt
import os
PathToARIANNAanalysis = os.environ['ARIANNAanalysis']

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


def findZenith(depth):
# depth (m), launch angle, arival rangle, predicted arrival angle, reflected angle
    datafile = PathToARIANNAanalysis + '/data/AngleVsDepthSP_650m_C.csv'
    data = np.loadtxt(datafile, skiprows=1, delimiter=',').T
    depths = data[0]
    index = find_nearest(depths, depth)
    return [depths[index],data[1][index],data[2][index],data[3][index]]

def findZenAziFromSpiceDataLPDAs(depth):
# depth (m), launch angle, arival rangle, predicted arrival angle
    datafile = PathToARIANNAanalysis + '/data/angleDataLPDASpulserGoingDownDepth750mDownFilterRectangular80to300MHz.npy'
    data = np.load(datafile,allow_pickle=True,encoding='bytes')
    depths = data[0]
    index = find_nearest(depths, depth)
    return [depths[index],[],[],data[1][index]]

def findZenAziFromSpiceDataLPDAsNoCableAdjustment(depth):
# depth (m), launch angle, arival rangle, predicted arrival angle
    datafile = PathToARIANNAanalysis + '/data/angleDataLPDASpulserGoingupDepth750mDownFilterRectangular80to300MHz.npy'
    data = np.load(datafile,allow_pickle=True,encoding='bytes')
    depths = data[0]
    index = find_nearest(depths, depth)
    return [depths[index],[],[],data[1][index]]

def findZenithSmooth(depth):
# depth (m), launch angle, arival rangle, predicted arrival angle
    datafile = PathToARIANNAanalysis + '/data/expectedArrivalDirectionSpice2018smooth.npy'
    data = np.load(datafile)
    depths = data[0]
    index = find_nearest(depths, depth)
    return [depths[index],[],[],data[1][index]]

def findZenithSmooth2(depth):
# depth (m), launch angle, arival rangle, predicted arrival angle
    datafile = PathToARIANNAanalysis + '/data/expectedArrivalDirectionSpice2018smoothNsurfaceAt1_353.npy'
    data = np.load(datafile)
    depths = data[0]
    index = find_nearest(depths, depth)
    return [depths[index],data[2][index],[],data[1][index],data[3][index]]

def findZenithSmoothLargeDepthsRange(depth):
# depth (m), launch angle, arival rangle, predicted arrival angle
    datafile = PathToARIANNAanalysis + '/data/expectedArrivalDirectionSpice2018LargerDepthRangeSmoothNsurfaceAt1_353.npy'
    data = np.load(datafile)
    depths = data[0]
    index = find_nearest(depths, depth)
    return [depths[index],data[2][index],[],data[1][index]]

def findZenithN1_378AtSurface(depth):
# depth (m), launch angle, arival rangle, predicted arrival angle
    datafile = PathToARIANNAanalysis + '/data/expectedArrivalDirectionSpice2018smoothNsurfaceAt1_378.npy'
    data = np.load(datafile)
    depths = data[0]
    index = find_nearest(depths, depth)
    return [depths[index],[],[],data[1][index]]

def findZenithTilt1degreeTowards(depth):
# depth (m), launch angle, arival rangle, predicted arrival angle
    datafile = PathToARIANNAanalysis + '/data/expectedArrivalDirectionSpice2018smoothWithTilt1degreeTowardsStation.npy'
    data = np.load(datafile)
    depths = data[0]
    index = find_nearest(depths, depth)
    return [depths[index],[],[],data[1][index]]

def findZenithTilt05degreeTowards(depth):
# depth (m), launch angle, arival rangle, predicted arrival angle
    datafile = PathToARIANNAanalysis + '/data/expectedArrivalDirectionSpice2018smoothWithTiltHalfDegreeTowardsStation.npy'
    data = np.load(datafile)
    depths = data[0]
    index = find_nearest(depths, depth)
    return [depths[index],[],[],data[1][index]]


def plotZenith():
    dataFile = PathToARIANNAanalysis + '/data/expectedArrivalDirectionSpice2018LargerDepthRangeSmoothNsurfaceAt1_353.npy'
    data = np.load(dataFile)

    fig, ax = plt.subplots(1, 1,figsize=(12, 5))
    ax.plot(data[0],data[1],'o',label='arrival zen.')
    ax.plot(data[0],np.asarray(data[2]),'o',label='Launch zen.')
    #ax.plot(data[0],data[3],'o',label='arrival zen. calced')
    mask = np.abs(data[2] - 15.0) < 0.01

    ax.set_ylabel('Zen')
    plt.axvline(x=400.0,label='shadow_zone',color='black')
    plt.legend()
    ax.set_xlabel('Pulser depth (m)')
    save = PathToARIANNAanalysis + "/plots/ArrivalVLaunchAnglesWiderRangeOfDepths.png"
    fig.savefig(save)
    save = PathToARIANNAanalysis + "/plots/ArrivalVLaunchAnglesWiderRangeOfDepths.pdf"
    fig.savefig(save)
    plt.show()


def main():
    depth = 100#m
    print('finding arrival angle at depth: ' + str(depth))
    print('launch angle: ' + str(findZenith(depth)[1]))
    print('arrival angle: ' + str(findZenith(depth)[2]))
    plotZenith()

if __name__== "__main__":
    main()
