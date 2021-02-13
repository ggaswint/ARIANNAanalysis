import matplotlib.pyplot as plt
import numpy as np
from NuRadioMC.SignalProp import analyticraytracing as ray
from NuRadioMC.utilities import medium
from NuRadioReco.utilities import units
import os

PathToARIANNAanalysis = os.environ['ARIANNAanalysis']

def getX(tilt, deltaY):
    return np.abs(np.tan(tilt) * deltaY)

def getXvaluesFromSpiceTild(z_space):
    z_space = np.flip(z_space)
    x = []
    xOffset1 = 0
    xOffset2 = 0
    xOffset3 = 0
    zPrev = -1400.0
    mask = (z_space >= -1600.0) & (z_space < -1400.0)
    numZValues = len(z_space[mask])
    angleChangesBetween1400And1600 = np.linspace(1.0,2.25,numZValues)
    angleIdx = 0
    for z in z_space:
        if z >= -1000.0:
            x.append(0)
        elif z >= -1200.0:
            xShift = getX(np.deg2rad(0.75),np.abs(z+1000.0))
            x.append(xShift)
            xOffset1 = x[len(x)-1]
        elif z >= -1400.0:
            xShift = getX(np.deg2rad(1.0),np.abs(z+1200.0))
            x.append(xShift + xOffset1)
            xOffset2 = x[len(x)-1]
        elif z >= -1600.0:
            xShift = getX(np.deg2rad(angleChangesBetween1400And1600[angleIdx]),np.abs(z-zPrev))
            x.append(xShift + xOffset2)
            xOffset2 = x[len(x)-1]
            zPrev = z
            angleIdx += 1
        else:
            xShift = getX(np.deg2rad(2.25),np.abs(z+1600.0))
            x.append(xShift + xOffset2)
    return np.flip(np.asarray(x))

def getPerpiniduclarXchange(spiceTiltLocations,r):
    azimuthAngleChanges = []
    new_x = []
    for x in spiceTiltLocations:
        new_x.append(np.sqrt(x**2 + r**2))
        azimuthAngleChanges.append(np.arctan(x/r))
    return -1*np.asarray(new_x), np.rad2deg(np.asarray(azimuthAngleChanges))

def getAnglesFromSpice(x_space,z_space):
    x2 = [0., -1.0]  # ARIANNA antanna
    ice = medium.southpole_simple()
    z_0_best = ice.z_0
    z_0 = 80.0
    delta_n = 1.78-1.353

    #x_space = xy_space - 653.8
    thetas = np.zeros(len(z_space))
    ice.z_0 = z_0
    ice.delta_n = delta_n
    ice.n_ice = 1.78
    r = ray.ray_tracing_2D(ice)
    recieveAngles = []
    for i, z in enumerate(z_space):
        x1 = [x_space[i] * units.m , z * units.m]
        #x1[0] = x_space[i] * units.m
        solution = r.find_solutions(x1, x2)
        recieveAngle = r.get_receive_angle(x1,x2,solution[0]['C0'])
        recieveAngles.append(np.rad2deg(recieveAngle))
    return np.asarray(recieveAngles)


if __name__ == "__main__":
    bottom = -1750.0
    top = -425.0
    r = 653.804524 # meters
    phi = 312.448284

    fig1, ax1 = plt.subplots(1, 1)
    z_space = np.linspace(bottom,top,int(np.abs(bottom - top))+1)
    x_space = getXvaluesFromSpiceTild(z_space)
    azimuth_space = np.linspace(phi,phi,int(np.abs(bottom - top))+1)

    ax1.fill_between([0,20],[-1750,-1750],color='skyblue')
    ax1.axhline(y=-200.0,color='black',label='firn boundary')
    ax1.plot(x_space,z_space,color='red',label=r'$\Delta X$ from center')
    ax1.set_xlabel(r'X [m]')
    ax1.set_ylabel(r'Z [m]')
    ax1.set_ylim(bottom,top)
    ax1.set_xlim(0,max(x_space))
    fig1.tight_layout()
    fig1.savefig(PathToARIANNAanalysis + '/plots/spiceTiltFromSurface.png')
    fig1.savefig(PathToARIANNAanalysis + '/plots/spiceTiltFromSurface.pdf')


    fig2, ax2 = plt.subplots(1, 1)

    z_space = np.linspace(bottom,top,int(np.abs(bottom - top))+1)
    x_space = -1*np.linspace(r,r,int(np.abs(bottom - top))+1)
    changesInSpiceHoleRelativeToSurface = getXvaluesFromSpiceTild(z_space)
    recieveAnglesNoTilt = getAnglesFromSpice(x_space,z_space)

    x_space = -r + changesInSpiceHoleRelativeToSurface
    recieveAngles = getAnglesFromSpice(x_space,z_space)
    ax2.plot(-z_space / units.m, recieveAngles-recieveAnglesNoTilt,linewidth=3,linestyle='solid',label='towards')

    x_space = -r - changesInSpiceHoleRelativeToSurface
    recieveAngles = getAnglesFromSpice(x_space,z_space)
    ax2.plot(-z_space / units.m, recieveAngles-recieveAnglesNoTilt,linewidth=2,linestyle='dotted',label='away')

    x_space, azimuthAngleChanges = getPerpiniduclarXchange(changesInSpiceHoleRelativeToSurface,r)
    recieveAngles = getAnglesFromSpice(x_space,z_space)
    ax2.plot(-z_space / units.m, recieveAngles-recieveAnglesNoTilt,linewidth=1,linestyle='dashed',label='perpindicular')

    ax2.legend()
    ax2.set_xlabel(r'Z [m]')
    ax2.set_ylabel(r'$\Delta$ receive zenith from no tilt [$^{\circ}$]')
    ax2.set_xlim(980,1700.0)
    fig2.tight_layout()
    fig2.savefig(PathToARIANNAanalysis + '/plots/changeInRecieveAngleZenithWithTilt.png')
    fig2.savefig(PathToARIANNAanalysis + '/plots/changeInRecieveAngleZenithWithTilt.pdf')


    fig3, ax3 = plt.subplots(1, 1)
    ax3.plot(-z_space / units.m, azimuth_space-azimuth_space,linewidth=3,linestyle='solid',label='towards')
    ax3.plot(-z_space / units.m, azimuth_space-azimuth_space,linewidth=2,linestyle='dotted',label='away')
    ax3.plot(-z_space / units.m, azimuthAngleChanges,linewidth=1,linestyle='dashed',label='perpindicular')
    ax3.legend()
    ax3.set_xlabel(r'Z [m]')
    ax3.set_ylabel(r'$\Delta$ receive azimuth from no tilt [$^{\circ}$]')
    ax3.set_xlim(980,1700.0)
    fig3.tight_layout()
    fig3.savefig(PathToARIANNAanalysis + '/plots/changeInRecieveAngleAzimuthWithTilt.png')
    fig3.savefig(PathToARIANNAanalysis + '/plots/changeInRecieveAngleAzimuthWithTilt.pdf')

    plt.show()
