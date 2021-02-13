import os
import numpy as np
import matplotlib.pyplot  as plt
import scipy
from NuRadioMC.SignalProp import analyticraytracing as ray
from NuRadioMC.utilities import medium
from NuRadioReco.utilities import units

PathToARIANNAanalysis = os.environ['ARIANNAanalysis']


def getExpectedReflected(z,xless=False):
    delta_n = 1.78-1.353
    x1 = [-653.8 * units.m, -800. * units.m]  # SPICE at 800m
    x2 = [0., -1.0]  # ARA antanna
    if xless:
        x1[0] = -680.0 * units.m
    ice = medium.southpole_simple()
    z_0 = 80.0
    ice.z_0 = z_0
    ice.delta_n = delta_n
    ice.n_ice = 1.78
    r = ray.ray_tracing_2D(ice)
    x1[1] = z * units.m
    solution = r.find_solutions(x1, x2)
    data = r.get_path(x1,x2,solution[1]['C0'])
    return data

def getExpected(z,xless=False):
    delta_n = 1.78-1.353
    x1 = [-653.8 * units.m, -800. * units.m]  # SPICE at 800m
    x2 = [0., -1.0]  # ARA antanna
    if xless:
        x1[0] = -680.0 * units.m
    ice = medium.southpole_simple()
    z_0 = 80.0
    ice.z_0 = z_0
    ice.delta_n = delta_n
    ice.n_ice = 1.78
    r = ray.ray_tracing_2D(ice)
    x1[1] = z * units.m
    solution = r.find_solutions(x1, x2)
    data = r.get_path(x1,x2,solution[0]['C0'])
    return data

def main():
    fig, ax = plt.subplots(1, 1,figsize=(5,7))
    data = getExpected(-1000.0)
    plt.plot(data[0],data[1])
    data = getExpectedReflected(-1000.0)
    plt.plot(data[0],data[1])
    data = getExpected(-1700.0)
    plt.plot(data[0],data[1],'--')
    data = getExpected(-418.0)
    plt.plot(data[0],data[1],'--')
    plt.fill_between([-800,0],[-1750,-1750],color='skyblue')
    plt.fill_between([-800,0],[-200,-200],color='aliceblue')
    data = getExpected(-440.0,xless=True)
    ys = scipy.zeros(len(data[1]))
    ys[:] = 10.0
    ax.fill_between(data[0], data[1], ys, color='grey',alpha=0.5)
    plt.xlim(-680.0,0)
    plt.ylim(-1750,10)
    plt.ylabel('depth [m]')
    plt.xlabel('distance [m]')

    plt.annotate('ARIANNA\n  STN', xy=(0, 0), xytext=(10, -240),size=16)
    plt.annotate('', xy=(0, 0), xytext=(100, -100),arrowprops=dict(facecolor='black',width=0.5,headwidth=5))

    ax.set_aspect('equal')

    ax.axvline(x=-653.8,color='black')
    fig.tight_layout()
    fig.savefig(PathToARIANNAanalysis + '/plots/rayTracingSouthPole.png')
    fig.savefig(PathToARIANNAanalysis + '/plots/rayTracingSouthPole.pdf')

    plt.show()

if __name__== "__main__":
    main()
