import numpy as np
import matplotlib.pyplot  as plt
from scipy.signal import hilbert
from radiotools import helper as hp
from scipy.interpolate import interp1d
import csv
from NuRadioMC.SignalProp import analyticraytracing as ray
from NuRadioMC.utilities import medium
from NuRadioReco.utilities import fft as FFTs
from NuRadioReco.utilities import units
from NuRadioReco.utilities import ice
from scipy.interpolate import UnivariateSpline
from scipy import signal
from NuRadioReco.detector import antennapattern # 180, irrelevant, 90, 0 ,90, [0,90 or -90]
import sys
import helpfulUtils as hu
from radiotools import plthelpers as php
from mpl_toolkits.axes_grid1 import make_axes_locatable
from scipy.stats import norm
from matplotlib.ticker import NullFormatter
from radiotools import stats as stat
import os
PathToARIANNAanalysis = os.environ['ARIANNAanalysis']

antenna_provider = antennapattern.AntennaPatternProvider()

lowers = [180,80,135,80]
uppers = [260,300,300,135]

lower_bound = 80
upper_bound = 300

color = ['C1','C2','C3','C4','C5','C6','C7','C8']


#lower_bound = lowers[0]
#upper_bound = uppers[0]

title = str(lower_bound) + 'MHz - ' + str(upper_bound) + 'MHz'


def gaussianHist(ax,gg_exp_diff,color,label,pos,line,label2):
    (mu, sigma) = norm.fit(gg_exp_diff)
    n, bins, patches = ax.hist(gg_exp_diff,linestyle=line,bins=np.arange(-15.0, 15.1, 1.0),edgecolor=color,fill=False,label=label2,histtype='step')#,weights=weights) histype = bar
    ax.text(1.0,0.99,label,horizontalalignment='right',verticalalignment='top',transform=ax.transAxes,color='midnightblue')


def getMedianStr(data):
    tweights = np.ones_like(data)
    from radiotools import stat as stat
    q1 = stat.quantile_1d(data, tweights, 0.16)
    q2 = stat.quantile_1d(data, tweights, 0.84)
    median = stat.median(data, tweights)
    textstr = "$\mathrm{median} = %.2g^{+%.2g}_{-%.2g}$" % (median, np.abs(median - q2),np.abs(median - q1))
    return textstr

def getMeanSTDStr(data):
    #gg_exp_diff2 = []
    #for j in range(len(data)):
    #    if np.abs(data[j]) < cut_val:
    #        gg_exp_diff2.append(data[j])
    #data = gg_exp_diff2
    mean = np.mean(data)
    std = np.std(data)
    textstr1 = r"$\mu$ = %.2g" % (mean)
    textstr2 = r"$\sigma$ = %.2g" % (std)
    #return [textstr1, textstr2]
    textstr3 = r"mean = %.2g$^{\circ}$, STD = %.2g$^{\circ}$" % (round(mean,2),round(std,2))
    return textstr3


def readData(my_data,debug=False):
    my_data.T[0][17:]
    x_arr = my_data.T[0][17:]
    y_arr = my_data.T[1][17:]
    delta_t = x_arr[1] - x_arr[0]

    if debug == True:
        f, ax = plt.subplots(2,1,figsize=(12, 15))
        ax[0].plot(x_arr, y_arr*1e3)
        ax[0].grid(False)
        ax[0].set_xlabel('[ns]')
        ax[0].set_ylabel('[mV]')

        freqs = np.fft.rfftfreq(len(y_arr), delta_t)

        ax[1].plot(freqs*1e-9, np.abs(np.fft.rfft(y_arr, norm='ortho')))
        ax[1].tick_params('y', colors='k')
        ax[1].grid(False)
        ax[1].set_xlabel('[GHz]')
        #ax[2].set_xlim(0.085,0.5)
        suptitle = 'vpol - raw ' + title
        f.suptitle(suptitle)
    return x_arr, y_arr, delta_t

def preprocess(etheta,ephi,x_theta,x_phi,delta_t,depth):
    ######################################## bandpass filter

    delta_t = x_theta[1]-x_theta[0]
    f_temp = np.fft.rfftfreq(len(etheta), delta_t)
    ftheta = f_temp
    etheta_tilda = np.fft.rfft(etheta, norm='ortho')

    #mask = f_temp*1e-6 < lower_bound
    #etheta_tilda[mask] = 0
    #mask = f_temp*1e-6 > upper_bound
    #etheta_tilda[mask] = 0

    delta_f = f_temp[1]-f_temp[0]
    Nt = (len(f_temp)-1)*2
    delta_t_test = 1.0/(delta_f*Nt)
    t = np.arange(Nt)*delta_t_test
    etheta = np.fft.irfft(etheta_tilda, norm='ortho')
    x_theta = t

    f_temp = np.fft.rfftfreq(len(ephi), delta_t)
    fphi = f_temp
    ephi_tilda = np.fft.rfft(ephi, norm='ortho')

    #mask = f_temp*1e-6 < lower_bound
    #ephi_tilda[mask] = 0
    #mask = f_temp*1e-6 > upper_bound
    #ephi_tilda[mask] = 0


    delta_f = f_temp[1]-f_temp[0]
    Nt = (len(f_temp)-1)*2
    delta_t_test = 1.0/(delta_f*Nt)
    t = np.arange(Nt)*delta_t_test
    ephi = np.fft.irfft(ephi_tilda, norm='ortho')
    x_phi = t

    ###### pad
    amount = 1000
    etheta = np.pad(etheta, (amount,amount), 'constant', constant_values=(0, 0))

    delta_t = x_theta[1] - x_theta[0]
    times_pad = []
    for i in range(amount):
       times_pad.append(x_theta[0] - delta_t*(amount-1-i))
    for i in range(len(x_theta)):
       times_pad.append(x_theta[i] + delta_t)
    for i in range(amount):
       times_pad.append(x_theta[len(x_theta)-1] + delta_t*(i+2))
    x_theta = times_pad


    amount = 1000
    ephi = np.pad(ephi, (amount,amount), 'constant', constant_values=(0, 0))

    delta_t = x_phi[1] - x_phi[0]
    times_pad = []
    for i in range(amount):
       times_pad.append(x_phi[0] - delta_t*(amount-1-i))
    for i in range(len(x_phi)):
       times_pad.append(x_phi[i] + delta_t)
    for i in range(amount):
       times_pad.append(x_phi[len(x_phi)-1] + delta_t*(i+2))
    x_phi = times_pad


    ############################ interpolate
    file = PathToARIANNAanalysis + '/data/AntennasResponseFrequenciesLPDA_VEL_180_0_ds_air.npy'
    ff = np.load(file,allow_pickle=True)
    ff = np.linspace(0,1.0,500)
    LPDA = antenna_provider.load_antenna_pattern("createLPDA_100MHz_InfAir")
    VEL = LPDA.get_antenna_response_vectorized(ff,90*units.deg, 0, 90*units.deg, 0, 180*units.deg, 0)
    VEL2 = LPDA.get_antenna_response_vectorized(ff, 90*units.deg, 0, 90*units.deg, 0, 90*units.deg, 90*units.deg)
    efield_antenna_factor = [[VEL['theta'],VEL['phi']],[VEL2['theta'],VEL2['phi']]]

    freq = ff*1e9



    delta_f = freq[1]-freq[0]
    Nt = (len(freq)-1)*2
    delta_t_test = 1.0/(delta_f*Nt)
    VEL_times = np.arange(Nt)*delta_t_test

    func = interp1d(x_theta,etheta)
    etheta = func(VEL_times)
    etheta_tilda = np.fft.rfft(etheta, norm='ortho')

    delta_t = VEL_times[1] - VEL_times[0]
    freqs = np.fft.rfftfreq(len(etheta), delta_t)

    func = interp1d(x_phi,ephi)
    ephi = func(VEL_times)
    ephi_tilda = np.fft.rfft(ephi, norm='ortho')

    delta_t = VEL_times[1] - VEL_times[0]
    freqs = np.fft.rfftfreq(len(ephi), delta_t)

    V = [ephi_tilda,etheta_tilda]

    delta_f = freq[1]-freq[0]
    Nt = (len(freq)-1)*2
    delta_t_test = 1.0/(delta_f*Nt)
    t = np.arange(Nt)*delta_t_test
    etheta = np.fft.irfft(etheta_tilda, norm='ortho')
    x_theta = t



    delta_f = freq[1]-freq[0]
    Nt = (len(freq)-1)*2
    delta_t_test = 1.0/(delta_f*Nt)
    t = np.arange(Nt)*delta_t_test
    ephi = np.fft.irfft(ephi_tilda, norm='ortho')
    x_phi = t

    ####### deconvolve manually # debug

    if 0:
        E1 = np.zeros_like(V[0])
        E2 = np.zeros_like(V[1])
        E1 = etheta_tilda/efield_antenna_factor[0][1]
        E2 = ephi_tilda/efield_antenna_factor[1][0]

    if 0:
        etheta_direct = etheta_tilda/efield_antenna_factor[0][1]
        ephi_direct = ephi_tilda/efield_antenna_factor[1][0]

        mask = freqs*1e-6 < lower_bound
        etheta_direct[mask] = 0
        mask = freqs*1e-6 > upper_bound
        etheta_direct[mask] = 0

        mask = freqs*1e-6 < lower_bound
        ephi_direct[mask] = 0
        mask = freqs*1e-6 > upper_bound
        ephi_direct[mask] = 0

        delta_f = freq[1]-freq[0]
        Nt = (len(freq)-1)*2
        delta_t_test = 1.0/(delta_f*Nt)
        t = np.arange(Nt)*delta_t_test
        e_times = t
        etheta2 = np.fft.irfft(etheta_direct, norm='ortho')
        ephi2 = np.fft.irfft(ephi_direct, norm='ortho')



    ################################## deconvolve antenna
    # efield_antenna_factor[0][0] = h_theta_first
    # efield_antenna_factor[0][1] = h_phi_first
    # efield_antenna_factor[1][0] = h_theta_2nd
    # efield_antenna_factor[1][1] = h_phi_second
    # V[0] = Vco
    # V[1] = Vcross 2nd antenna
    # E0 = etheta
    # E1 = ephi
    if 1:
        n_frequencies = len(V[0])
        denom = (efield_antenna_factor[0][0] * efield_antenna_factor[1][1] - efield_antenna_factor[0][1] * efield_antenna_factor[1][0])
        mask = np.abs(denom) != 0
        # solving for electric field using just two orthorgonal antennas
        E1 = np.zeros_like(V[0])
        E2 = np.zeros_like(V[1])
        E1[mask] = (V[0] * efield_antenna_factor[1][1] - V[1] * efield_antenna_factor[0][1])[mask] / denom[mask]
        E2[mask] = (V[1] - efield_antenna_factor[1][0] * E1)[mask] / efield_antenna_factor[1][1][mask]

    if 0:
        n_frequencies = len(V[0])
        denom = (efield_antenna_factor[1][0] * efield_antenna_factor[0][1] - efield_antenna_factor[1][1] * efield_antenna_factor[0][0])
        mask = np.abs(denom) != 0
        # solving for electric field using just two orthorgonal antennas
        E1 = np.zeros_like(V[0])
        E2 = np.zeros_like(V[1])
        E1[mask] = (V[0] * efield_antenna_factor[0][1] - V[1] * efield_antenna_factor[1][1])[mask] / denom[mask]
        E2[mask] = (V[1] - efield_antenna_factor[0][0] * E1)[mask] / efield_antenna_factor[0][1][mask]


    #### frequency dependent propagation attenuation
    ice = medium.southpole_simple()
    z_0 = 80.0
    delta_n = 1.78-1.353
    ice.z_0 = z_0
    ice.delta_n = delta_n
    ice.n_ice = 1.78
    r = ray.ray_tracing_2D(ice)
    x1 = [-653.8 * units.m, depth * units.m]  # SPICE at 800m
    x2 = [0.* units.m, -1.0* units.m]  # ARA antanna
    solution = r.find_solutions(x1, x2)
    C_0 = 1.0 # solution['C0']
    max_detector_freq = 1*units.GHz

    freq_shift = 1.78
    freqs = freqs / freq_shift

    attenuation = r.get_attenuation_along_path(x1, x2, C_0, freq*units.Hz, max_detector_freq)
    #print(attenuation)
    #E1 = attenuation*E1
    #E2 = attenuation*E2

    etheta = np.fft.irfft(E1, norm='ortho')
    ephi = np.fft.irfft(E2, norm='ortho')

    #freq_shift = np.sqrt(1.78)

    mask = freqs*1e-6 < lower_bound
    E1[mask] = 0
    mask = freqs*1e-6 > upper_bound
    E1[mask] = 0

    mask = freqs*1e-6 < lower_bound
    E2[mask] = 0
    mask = freqs*1e-6 > upper_bound
    E2[mask] = 0


    delta_f = freqs[1]-freqs[0]
    Nt = (len(freqs)-1)*2
    delta_t_test = 1.0/(delta_f*Nt)
    t = np.arange(Nt)*delta_t_test
    e_times = t
    etheta = np.fft.irfft(E1, norm='ortho')
    ephi = np.fft.irfft(E2, norm='ortho')

    delta_t = e_times[1]-e_times[0]
    f_temp = np.fft.rfftfreq(len(etheta), delta_t)
    etheta_tilda = np.fft.rfft(etheta, norm='ortho')

    #np.argmax(etheta)
    #plt.plot(x_theta*1e9,etheta)
    #plt.plot(e_times*1e9,etheta)
    #plt.show()

    return etheta, ephi,e_times# x_theta



def plot_electric_fields(etheta,ephi,x_theta,angle,set_mltp,i):
    exact = True

    e_times = x_theta
    h_etheta = hilbert(etheta)
    h_ephi = hilbert(ephi)
    h3 = np.sqrt(np.abs(h_etheta)**2 + np.abs(h_ephi)**2)
    fwhm = hp.get_FWHM_hilbert(h3)

    if exact:
        mltp = int(set_mltp*1e-9/(e_times[1]-e_times[0]))
        width = fwhm[1]-fwhm[0]
        fwhm = hp.get_FWHM_hilbert(np.abs(h_etheta))
        mid_fwhm = fwhm[0] + int((fwhm[1] - fwhm[0])/2)

        IW = int(mid_fwhm+mltp/2) - int(mid_fwhm-mltp/2)
        noise_start = int(1.1*int(mid_fwhm+mltp/2))
        signal_start = int(mid_fwhm+mltp/2)
        #max_etheta = np.sqrt(np.abs(np.sum((np.abs(etheta[int(mid_fwhm-mltp/2):int(mid_fwhm+mltp/2)]))**2) - np.sum((np.abs(etheta[noise_start:noise_start+IW]))**2)))
        #max_ephi = np.sqrt(np.abs(np.sum((np.abs(ephi[int(mid_fwhm-mltp/2):int(mid_fwhm+mltp/2)]))**2) - np.sum((np.abs(ephi[noise_start:noise_start+IW]))**2)))

        max_etheta = np.sqrt(np.abs(np.sum((np.abs(etheta[signal_start-IW:signal_start]))**2) - np.sum((np.abs(etheta[noise_start:noise_start+IW]))**2)))
        max_ephi = np.sqrt(np.abs(np.sum((np.abs(ephi[signal_start-IW:signal_start]))**2) - np.sum((np.abs(ephi[noise_start:noise_start+IW]))**2)))

    else:
        theta_idx = [[181, 205, 115, 161],[181, 205, 115, 161],[175, 199, 109, 155],[175, 196, 104, 150],[175, 196, 107, 151],[165, 186, 99, 145],[190, 209, 86, 154]]
        phi_idx = [[181, 204, 115, 161],[181, 204, 115, 161],[175, 201, 109, 155],[170, 190, 104, 150],[170, 185, 107, 151],[165, 185, 99, 145],[185, 204, 86, 154]]
        width = fwhm[1]-fwhm[0]
        width1 = int(0.6*width)
        idx = np.argmax(etheta)
        theta_idx = np.asarray(theta_idx[i])*50
        phi_idx = np.asarray(phi_idx[i])*50
        max_etheta = np.sqrt(np.sum((etheta[theta_idx[0]:theta_idx[1]])**2) - ((theta_idx[1]-theta_idx[0])/(theta_idx[3]-theta_idx[2])) * np.sum((etheta[theta_idx[2]:theta_idx[3]])**2))
        max_ephi = np.sqrt(np.sum((ephi[phi_idx[0]:phi_idx[1]])**2) - ((phi_idx[1]-phi_idx[0])/(phi_idx[3]-phi_idx[2])) * np.sum((ephi[phi_idx[2]:phi_idx[3]])**2))
    if 0:
        f, ax = plt.subplots(2,1,figsize=(12, 15))
        ax[0].plot(e_times*1e9,etheta)
        #ax[0].plot(e_times,np.abs(h_etheta))
        ax[1].plot(e_times*1e9,ephi)
        #ax[1].plot(e_times,np.abs(h_ephi))
        ax[0].plot(e_times[signal_start-IW:signal_start]*1e9,etheta[signal_start-IW:signal_start],color='red')
        ax[1].plot(e_times[signal_start-IW:signal_start]*1e9,ephi[signal_start-IW:signal_start],color='red')
        ax[0].plot(e_times*1e9,np.abs(h_etheta),color='orange',alpha=0.25)
        ax[1].plot(e_times*1e9,np.abs(h_ephi),color='orange',alpha=0.25)
        ax[0].plot(e_times[noise_start:noise_start+IW]*1e9,etheta[noise_start:noise_start+IW],color='green')
        ax[1].plot(e_times[noise_start:noise_start+IW]*1e9,ephi[noise_start:noise_start+IW],color='green')
        ax[0].set_xlabel('[ns]')
        ax[0].set_ylabel('eTheta [v/m]')
        ax[1].set_xlabel('[ns]')
        ax[1].set_ylabel('ePhi [v/m]')
        ax[1].set_ylim(-1.1*np.max(ephi),1.1*np.max(ephi))
        ax[0].set_xlim(200,700)
        ax[1].set_xlim(200,700)
        #save = PathToARIANNAanalysis + '/data/expectedEfieldSpice2019noBandPassFilter'
        #np.save(savefile,[e_times*1e9,etheta,ephi])
        title = 'angle: ' + str(int(angle)) + ' pol: ' + str(np.rad2deg(np.arctan(np.asarray(max_ephi/max_etheta))))
        f.suptitle(title)

        save = PathToARIANNAanalysis + '/data/electricFieldExpectedSpiceData_'+str(lower_bound)+'_'+str(upper_bound)+'_'+str(int(i))+ ".png"

        f.savefig(save)
        plt.close(f)
        #print(savefile)
        #plt.show()


    polarization = float(np.rad2deg(np.arctan(max_ephi/max_etheta)))

    return [90.0-angle, polarization]


def run(lower_bound,upper_bound,sets):
    angles = [0,15,30,45,60,75,90]
    depths = [-200.0,-500.0,-1000.0,-1200.0,-1700.0,-1700.0,-1700.0]
    depths = [-171.8321655,-216.84784785,-317.1782894,-581.13113113,-1089.93993994,-1692.69402736,-1692.69402736]

    ################################################################################
    #   vpol
    ################################################################################
    main_x_array = []
    first = True

    x_ave_all_v = []
    y_ave_all_v = []
    delta_t_all_v = []
    x_all_v = []
    y_all_v = []
    delta_ts_all_v = []
    angles_nums = [1,61,11,21,31,41,51]
    for angle_num in angles_nums:
        x_ave = []
        y_ave = []
        delta_ts = []
        for file_num in range(10):
            num = str(angle_num+file_num)
            while len(num) < 4:
                num = '0' + num
            file_name = PathToARIANNAanalysis + '/data/AnechoicChamberData/tek'+num+'CH1.csv'
            my_data = np.genfromtxt(file_name, delimiter=',')
            x, y, delta_t = readData(my_data)
            mean = np.mean(y[0:800])
            y = np.asarray(y) - mean
            x_ave.append(x)
            y_ave.append(y)
            if first:
                first = False
                main_x_array = x
            delta_ts.append(delta_t)
            x_all_v.append(x)
            y_all_v.append(y)
            delta_ts_all_v.append(delta_t)
        x_ave_all_v.append(sum(x_ave)/10.0)
        y_ave_all_v.append(sum(y_ave)/10.0)
        delta_t_all_v.append(delta_ts)

    if 0:
        f, ax = plt.subplots(len(x_ave_all_v),2,figsize=(12, 15),sharey='col',sharex='col')
        for i in range(len(x_ave_all_v)):
            label=str(angles[i]) + r' $^\circ$'
            ax[i][0].plot(x_ave_all_v[i]*1e9, y_ave_all_v[i]*1e3,label=label)
            ax[i][0].grid(False)
            ax[i][0].set_xlim(-50.0,200.0)

            delta_t = delta_t_all_v[i][0]
            freqs = np.fft.rfftfreq(len(y_ave_all_v[i]), delta_t)

            ax[i][1].plot(freqs*1e-9, np.abs(np.fft.rfft(y_ave_all_v[i], norm='ortho')))
            ax[i][1].tick_params('y', colors='k')
            ax[i][1].grid(False)
            ax[i][1].set_xlim(0.0,0.5)
            ax[i][0].legend(loc='right')

        ax[3][1].set_xlabel('[GHz]')
        ax[3][0].set_xlabel('[ns]')
        ax[3][0].set_ylabel('[mV]')
        suptitle = 'vpol - raw ' + title
        f.suptitle(suptitle)
        save = PathToARIANNAanalysis + '/plots/anechoic_vpol_raw'+title+".png"
        f.savefig(save)




    ################################################################################
    #   hpol
    ################################################################################

    x_ave_all_h = []
    y_ave_all_h = []
    delta_t_all_h = []
    x_all_h = []
    y_all_h = []
    delta_ts_all_h = []
    angles_nums = [71,81,91,101,111,121,131]
    for angle_num in angles_nums:
        x_ave = []
        y_ave = []
        delta_ts = []
        for file_num in range(10):
            num = str(angle_num+file_num)
            while len(num) < 4:
                num = '0' + num
            file_name = PathToARIANNAanalysis + '/data/AnechoicChamberData/tek'+num+'CH1.csv'
            my_data = np.genfromtxt(file_name, delimiter=',')
            x, y, delta_t = readData(my_data)
            mean = np.mean(y[0:800])
            y = np.asarray(y) - mean
            func_volt = interp1d(x,y,fill_value=0)
            y = func_volt(main_x_array)
            x = main_x_array
            x_ave.append(x)
            y_ave.append(y)
            delta_ts.append(delta_t)
            x_all_h.append(x)
            y_all_h.append(y)
            delta_ts_all_h.append(delta_t)
        x_ave_all_h.append(sum(x_ave)/10.0)
        y_ave_all_h.append(sum(y_ave)/10.0)
        delta_t_all_h.append(delta_ts)

    if 0:
        f, ax = plt.subplots(len(x_ave_all_h),2,figsize=(12, 15),sharey='col',sharex='col')
        for i in range(len(x_ave_all_h)):
            label=str(angles[i]) + r' $^\circ$'
            ax[i][0].plot(x_ave_all_h[i]*1e9, y_ave_all_h[i]*1e3,label=label)
            ax[i][0].grid(False)
            ax[i][0].set_xlim(-50.0,200.0)

            delta_t = delta_t_all_h[i][0]
            freqs = np.fft.rfftfreq(len(y_ave_all_h[i]), delta_t)

            ax[i][1].plot(freqs*1e-9, np.abs(np.fft.rfft(y_ave_all_h[i], norm='ortho')))
            ax[i][1].tick_params('y', colors='k')
            ax[i][1].grid(False)
            ax[i][1].set_xlim(0.0,0.5)
            ax[i][0].legend(loc='right')

        ax[3][1].set_xlabel('[GHz]')
        ax[3][0].set_xlabel('[ns]')
        ax[3][0].set_ylabel('[mV]')
        suptitle = 'hpol - raw ' + title
        f.suptitle(suptitle)
        save = PathToARIANNAanalysis + '/plots/anechoic_hpol_raw'+title+".png"
        f.savefig(save)


    if 0: # test efield reconstruction
        file_name = PathToARIANNAanalysis + '/data/AnechoicChamberData/tek0021CH1.csv'
        my_data = np.genfromtxt(file_name, delimiter=',')
        x2, y2, delta_t2 = readData(my_data)

        file_name = PathToARIANNAanalysis + '/data/AnechoicChamberData/tek0101CH1.csv'
        my_data = np.genfromtxt(file_name, delimiter=',')
        x, y, delta_t = readData(my_data)
        func_volt = interp1d(x,y,fill_value=0)
        y = func_volt(x2)
        x = x2

        etheta = y2
        ephi = y
        delta_t = delta_t2
        x_theta = x2
        x_phi = x

        etheta, ephi, x_theta = preprocess(etheta,ephi,x_theta,x_phi,delta_t)

        fig1, ax1 = plt.subplots(1,2,figsize=(12, 15))
        fig2, ax2 = plt.subplots(1,1,figsize=(12, 15))
        plot_electric_fields(etheta, ephi, x_theta,angles,ax1,ax2)
        plt.show()
        print(1/0)


    ethetas = []
    ephis = []
    x_thetas = []
    for j in range(len(depths)):
        for i in range(10):
            etheta = y_all_v[j*10 + i]
            ephi = y_all_h[j*10 + i]
            delta_t = delta_ts_all_v[j*10 + i]
            x_theta2 = x_all_v[j*10 + i]
            x_phi = x_all_h[j*10 + i]
            etheta, ephi, x_theta2 = preprocess(etheta,ephi,x_theta2,x_phi,delta_t,depths[j])
            sample_rate = int(89.5/((x_theta2[1]-x_theta2[0])*1e9))  # sqrt(1.78),  3/1 * 256 , channelResampler look at
            sample_rate = int(40.0/((x_theta2[1]-x_theta2[0])*1e9))  # 1.78,  3/1 * 256 , channelResampler look at
            etheta, x_theta =signal.resample(etheta,sample_rate*len(x_theta2),x_theta2)
            ephi, x_theta =signal.resample(ephi,sample_rate*len(x_theta2),x_theta2)
            #print('delta t: ' + str(x_theta[1]-x_theta[0]))
            ethetas.append(etheta)
            ephis.append(ephi)
            x_thetas.append(x_theta)

    plot_data_pol_all = []
    for k in range(len(sets)):
        plot_data_pol = []
        for j in range(len(depths)):
            for i in range(10):
                   plot_data = plot_electric_fields(np.asarray(ethetas[j*10 + i]), np.asarray(ephis[j*10 + i]), np.asarray(x_thetas[j*10 + i]),angles[j],sets[k],j*10 + i)
                   plot_data_pol.append(plot_data)
        plot_data_pol_all.append(plot_data_pol)

    return np.asarray(plot_data_pol_all[0])

def process_polData(pol_data):
    pol_data = np.asarray(pol_data).T
    ave_x = []
    ave_y = []
    sigma = []
    q1 = []
    q2 = []
    from radiotools import stats as stat
    for j in range(7):
        x = []
        y = []
        for i in range(10):
            x.append(pol_data[0][j*10 + i])
            y.append(pol_data[1][j*10 + i])
        ave_x.append(sum(x)/10.0)
        ave_y.append(sum(y)/10.0)
        sigma.append(np.std(y))
        tweights = np.ones_like(y)
        q1.append(stat.quantile_1d(y, tweights, 0.16))
        q2.append(stat.quantile_1d(y, tweights, 0.84))
    q1 = np.asarray(ave_y) - np.asarray(q1)
    q2 = np.asarray(q2) - np.asarray(ave_y)
    pol_data = [ave_x, ave_y]
    pol_data = np.asarray(pol_data).T
    return pol_data, sigma, q1, q2


def aveError(depth,data):
    depth = np.asarray(depth)
    depth = np.round(depth.astype(float), -1)
    #depth = np.unique(np.round(depth.astype(float), -1))
    means = []
    stds = []
    depths = []
    for d in np.unique(np.round(depth.astype(float), -1)):
        mask = depth==d
        means.append(data[mask].mean())
        stds.append(data[mask].std())
        depths.append(d)
    return means, stds, depths


def aveError2(depth,data):
    depth = np.asarray(depth)
    depth = np.round(depth.astype(float), -1)
    #depth = np.unique(np.round(depth.astype(float), -1))
    means = []
    stds = []
    depths = []
    q1s = []
    q2s = []
    for d in np.unique(np.round(depth.astype(float), -1)):
        mask = depth==d
        means.append(data[mask].mean())
        stds.append(data[mask].std())
        depths.append(d)
        ave = sum(data[mask])/len(data[mask])
        tweights = np.ones_like(data[mask])
        q1 = stat.quantile_1d(data[mask], tweights, 0.16)
        q2 = stat.quantile_1d(data[mask], tweights, 0.84)
        median = data[mask].mean()
        q1s.append(np.abs(median-q1))
        q2s.append(np.abs(median-q2))
    return means, q1s,q2s, depths



depths = np.arange(200,2501,1)
print(depths)
Zen_a = []
Zen_l = []
for depth in depths:
    Zen = hu.findZenithSmoothLargeDepthsRange(depth)[3]*units.deg
    Zen_a.append(Zen)
    Zen = hu.findZenithSmoothLargeDepthsRange(depth)[1]*units.deg
    Zen_l.append(Zen)
Azi = 312.448284*units.deg

delta_zen = Zen_a[1] - Zen_a[0]
pairs = []
angles = [90.0*units.deg,105.0*units.deg,120.0*units.deg,135.0*units.deg,150.0*units.deg,165.0*units.deg,180.0*units.deg]
for i in range(len(Zen_a)):
    for angle in angles:
        if np.abs((Zen_a[i]/units.deg-angle/units.deg)) < 0.5:
            pairs.append([str(Zen_l[i]/units.deg),str(Zen_a[i]/units.deg),str(depths[i])])


fig3, ax3 = plt.subplots(1,1,figsize=(4, 4))
fig, ax = plt.subplots(1, 1,figsize=(11, 7)) # average rectilinear 80-300
sets = [[17.5*4]]#[[875]]
labels = ['1.3','1.0','0.5']
file_ex = ['1','2','3']
print(title)
for bound in range(len(sets)):
    #lower_bound = lowers[1]
    #upper_bound = uppers[1]
    label = labels[bound] + 'FWHMa'
    pol_data = run(lower_bound,upper_bound,sets[bound])
    pol_data_all = np.copy(pol_data)
    pol_data_all = pol_data_all.T
    pol_data_all[0] = pol_data_all[0][::-1]
    pol_data_all[1] = pol_data_all[1][::-1]
    ax3.plot(pol_data_all[0][10:],pol_data_all[1][10:],'o',color=color[bound],alpha=0.25)

    pol_data, pol_errors, pol_err_m, pol_err_p = process_polData(pol_data)

    plot_data_pol = pol_data

    pol_data = pol_data.T
    pol_data[0] = pol_data[0][::-1]
    pol_data[1] = pol_data[1][::-1]
    pol_errors = pol_errors[::-1]
    pol_err_m = pol_err_m[::-1]
    pol_err_p = pol_err_p[::-1]



    #spl = UnivariateSpline(pol_data[0],pol_data[1])
    #xs = np.linspace(0, 90, 1000)
    ##plt.plot(xs, spl(xs), 'g', lw=3)
    #spl.set_smoothing_factor(0.5)
#
    #func = interp1d(xs,spl(xs))
    func_p1sigma = interp1d(pol_data[0],np.asarray(pol_data[1])+np.asarray(pol_err_p))
    interp_pol_data_y_p1sigma = func_p1sigma((np.asarray(Zen_l)/units.deg)[::-1])
    interp_pol_data_x_p1sigma = (np.asarray(Zen_l)/units.deg)[::-1]

    func_m1sigma = interp1d(pol_data[0],np.asarray(pol_data[1])-np.asarray(pol_err_m))
    interp_pol_data_y_m1sigma = func_m1sigma((np.asarray(Zen_l)/units.deg)[::-1])
    interp_pol_data_x_m1sigma = (np.asarray(Zen_l)/units.deg)[::-1]


    func = interp1d(pol_data[0],pol_data[1])
    interp_pol_data_y = func((np.asarray(Zen_l)/units.deg)[::-1])
    interp_pol_data_x = (np.asarray(Zen_l)/units.deg)[::-1]



    depths_cut2 = [1021.0875000000003, 1034.7041666666667, 1048.0833333333323, 1064.8666666666666, 1079.8291666666669, 1095.6624999999995, 1111.9472086588537, 1130.722208658854, 1150.7805684407554, 1167.6166697184242, 1183.7277852376303, 1201.7250172932945, 1219.1666666666665, 1254.4166666666672, 1290.4999999999995, 1305.2291666666656, 1317.4999999999995, 1330.8749999999995, 1345.25, 1360.7124999999996, 1376.7999999999997, 1392.2500000000007, 1406.950000000001, 1419.025000000001, 1431.1000000000008, 1445.5321350097659, 1459.1666829427081, 1471.6071594238279, 1484.6428568522138, 1497.6190450032555, 1510.0595357259115, 1537.678588867187, 1566.2499999999995, 1578.5714111328125, 1588.33332824707, 1598.035720825195, 1609.4444478352866, 1623.777803548177, 1637.4444478352862, 1647.7222290039062, 1657.5000081380206, 1667.9444569905597, 1680.000000000001, 1691.7777669270838, 1698.5000091552733]
    mask = (depths >= 1021.0875000000003) & (depths <= 1698.5000091552733)
    mask = (depths >= 938.0) & (depths <= 1700.0)

    depths_cut = [837.8611094156904, 850.8305506388351, 864.1375015258782, 879.0958343505853, 916.7972178141276, 956.3624959309902, 973.1152760823572, 988.3916666666669, 1003.9874999999997, 1021.0083333333334, 1034.7041666666667, 1048.0833333333323, 1064.8666666666666, 1079.8291666666669, 1095.6624999999995, 1111.9472086588537, 1130.722208658854, 1150.7805684407554, 1167.6166697184242, 1183.7277852376303, 1201.7250172932945, 1219.1666666666665, 1254.4166666666672, 1290.4999999999995, 1305.2291666666656, 1317.4999999999995, 1330.8749999999995, 1345.25, 1360.7124999999996, 1376.7999999999997, 1392.2500000000007, 1406.950000000001, 1419.025000000001, 1431.1000000000008, 1445.5321350097659, 1459.1666829427081, 1471.6071594238279, 1484.6428568522138, 1497.6190450032555, 1510.0595357259115, 1537.678588867187, 1566.2499999999995, 1578.5714111328125, 1588.33332824707, 1598.035720825195, 1609.4444478352866, 1623.777803548177, 1637.4444478352862, 1647.7222290039062, 1657.5000081380206, 1667.9444569905597, 1680.000000000001, 1691.7777669270838, 1698.5000091552733]
    #mask = (depths > 837.85) & (depths < 1698.5000091552733)

    depths2 = np.linspace(800.0,2500.0,1000)#depths_cut2

    x_interp_p = np.asarray(interp_pol_data_x_p1sigma[::-1][mask]).tolist()
    y_interp_p = np.asarray(interp_pol_data_y_p1sigma[::-1][mask]).tolist()
    x_interp_m = np.asarray(interp_pol_data_x_m1sigma[::-1][mask]).tolist()
    y_interp_m = np.asarray(interp_pol_data_y_m1sigma[::-1][mask]).tolist()

    x_interp_p1sigma = np.asarray(interp_pol_data_x_p1sigma[::-1][mask]).tolist()
    y_interp_p1sigma = np.asarray(interp_pol_data_y_p1sigma[::-1][mask]).tolist()

    func_p1sigma = interp1d(depths,interp_pol_data_y_p1sigma[::-1])
    interp_pol_data_y_p1sigma = func_p1sigma(depths2)
    interp_pol_data_x_p1sigma = depths2

    x_interp_m1sigma = np.asarray(interp_pol_data_x_m1sigma[::-1][mask]).tolist()
    y_interp_m1sigma = np.asarray(interp_pol_data_y_m1sigma[::-1][mask]).tolist()

    func_m1sigma = interp1d(depths,interp_pol_data_y_m1sigma[::-1])
    interp_pol_data_y_m1sigma = func_m1sigma(depths2)
    interp_pol_data_x_m1sigma = depths2



    x_interp = np.asarray(interp_pol_data_x[::-1][mask]).tolist()
    y_interp = np.asarray(interp_pol_data_y[::-1][mask]).tolist()

    func = interp1d(depths,interp_pol_data_y[::-1])
    interp_pol_data_y = func(depths2)
    interp_pol_data_x = depths2

    #spl = UnivariateSpline(depths_cut2,interp_pol_data_y)
    #xs = depths_cut2
    #plt.plot(xs, spl(xs), 'g', lw=3)
    #spl.set_smoothing_factor(0.5)
    #print(pol_data[1][2])
    #mask = np.abs((np.asarray(interp_pol_data_y) - pol_data[1][1])) <= 0.001
    #print(interp_pol_data_x[mask])
    #print(interp_pol_data_y)
    #print(1/0)

    #1772.29458918
    ax.errorbar([1772.29458918],[pol_data[1][1]],yerr=[[pol_err_m[1]], [pol_err_p[1]]],fmt='d',MarkerSize=8.0,color='darkorange',label='Anechoic')
    ax.errorbar([1102.0],[pol_data[1][2]],yerr=[[pol_err_m[2]], [pol_err_p[2]]],fmt='d',MarkerSize=8.0,color='darkorange',label='Anechoic')
    ax.fill_between(interp_pol_data_x, interp_pol_data_y_p1sigma, interp_pol_data_y_m1sigma,color='darkorange',label='Anechoic', alpha=.25)
    #np.save('/home/geoffrey/ARIANNA/SpiceData/anechoic_data_9FWHM',[interp_pol_data_x, interp_pol_data_y])
    #np.save('/home/geoffrey/ARIANNA/SpiceData/anechoic_data_line_4FWHM_atten',[interp_pol_data_x, interp_pol_data_y_p1sigma, interp_pol_data_y_m1sigma, interp_pol_data_y])
    #print(max(interp_pol_data_y))
    #print(1/0)


    chan_sets = ['0123']
    up_zen = []
    up_azi = []
    down_zen = []
    down_azi = []

    for num,chans in enumerate(chan_sets):
        ##file = '/home/geoffrey/ARIANNA/SpiceData/polData_dec30_800m_constTime_newAveref_sameChan_noE_single_test2_'+str(lower_bound)+'_'+str(upper_bound)+'.npy' # 1,2,4,6
        #file = '/home/geoffrey/ARIANNA/ARIANNAanalysis/ARIANNAanalysis/data/SPICE_Dec30_2018_EFieldDataLPDA_AtTransmitter.npy' # 1,2,4,6
        ##label = 'chans: '+chans
        #data = np.load(file,allow_pickle=True)
        #label = 'Spice'
        #ax.plot(data[0],np.rad2deg(np.arctan(np.asarray(data[1]))),'o',color=color[5],alpha=0.25)#,label=label)
        #means, errors, depths = aveError(data[0],np.rad2deg(np.arctan(np.asarray(data[1]))))
        #ax.errorbar(depths,means,errors,fmt='d',MarkerSize=10.0,color='purple',label='chan -5 degree zenith')
        #ax.set_ylabel(r'$\arctan{\frac{\phi}{\theta}}$ $[^\circ]$')
        #ax.text(0.99,0.90,'Spice; chan0 extra -0.5m depth',horizontalalignment='right',verticalalignment='top',transform=ax.transAxes,color='purple')

        shading_inputs = []

        file = PathToARIANNAanalysis + '/data/polarizationReconstructedDataSpice2018Experiment.npy' # 1,2,4,6
        #file = '/home/geoffrey/ARIANNA/polData_dec30_800m_studyRotationAllH_80_300_azi0_zen0_amount_2_5_v2.npy' # 1,2,4,6
        #label = 'chans: '+chans
        data = np.load(file,allow_pickle=True)
        data0 = np.delete(data[0],721)
        data1 = np.delete(data[1],721)

        label = 'Spice'
        mask_space1 = (np.asarray(data0) <= 920.0)
        mask_space2 = (np.asarray(data0) >= 920.0) & (np.asarray(data0) <= 1250.0)
        mask_space3 = (np.asarray(data0) >= 1250.0) & (np.asarray(data0) <= 1530.0)
        mask_space4 = (np.asarray(data0) >= 1530.0)

        #ax.plot(data0,np.rad2deg(np.arctan(np.asarray(data1))),'o',color='deepskyblue',alpha=0.25)#,label=label)
        #ax.plot(data[0],np.rad2deg(np.arctan(np.asarray(data[1]))),'--',color=color[1])#,label=label)
        mask_spice = (np.asarray(data0) >= 938.0)
        spice_depths = np.asarray(data0)[mask_spice]
        spice_data = np.rad2deg(np.arctan(np.asarray(data1)))[mask_spice]
        means, errors, depths = aveError(data0,np.rad2deg(np.arctan(np.asarray(data1))))
        mask_space1 = (np.asarray(depths) <= 920.0)
        mask_space2 = (np.asarray(depths) >= 920.0) & (np.asarray(depths) <= 1250.0)
        mask_space3 = (np.asarray(depths) >= 1250.0) & (np.asarray(depths) <= 1530.0)
        mask_space4 = (np.asarray(depths) >= 1530.0)
        depths = np.asarray(depths)
        means = np.asarray(means)
        errors = np.asarray(errors)
        shading_inputs.append(depths[mask_space1][-1])
        shading_inputs.append(depths[mask_space2][0])
        shading_inputs.append(depths[mask_space2][-1])
        shading_inputs.append(depths[mask_space3][0])
        shading_inputs.append(depths[mask_space3][-1])
        shading_inputs.append(depths[mask_space4][0])
        print(shading_inputs)

        mask_test = (np.asarray(depths) >= 1000.0) & (np.asarray(depths) <= 1010.0)
        print('errors: ' + str(errors[mask_test]))
        print('means: ' + str(means[mask_test]))
        print('depths: ' + str(depths[mask_test]))
        #means, q1,q2, depths = aveError2(data0,np.rad2deg(np.arctan(np.asarray(data1))))
        #print('68: ' + str(np.asarray(q1)[mask_test]) + ' ' + str(np.asarray(q1)[mask_test]))
        #print(1/0)

        ax.errorbar(depths[mask_space1],means[mask_space1],errors[mask_space1],fmt='d',MarkerSize=8.0,color='midnightblue',label=label)
        ax.errorbar(depths[mask_space2],means[mask_space2],errors[mask_space2],fmt='d',MarkerSize=8.0,color='midnightblue',label=label)
        ax.errorbar(depths[mask_space3],means[mask_space3],errors[mask_space3],fmt='d',MarkerSize=8.0,color='midnightblue',label=label)
        ax.errorbar(depths[mask_space4],means[mask_space4],errors[mask_space4],fmt='d',MarkerSize=8.0,color='midnightblue',label=label)


        #ax.set_ylabel(r'$\arctan{\frac{f_{\phi}}{f_{\theta}}}$ $[^\circ]$')
        ax.set_ylabel(r'polarization $[^\circ]$')
        file = PathToARIANNAanalysis + '/data/expectedPolarizationUpperQuantile68.npy'
        data_q1 = np.load(file,allow_pickle=True)
        file = PathToARIANNAanalysis + '/data/expectedPolarizationLowerQuantile68.npy'
        data_q2 = np.load(file,allow_pickle=True)
        mask_space1 = (np.asarray(data_q1[2]) <= 920.0)
        mask_space2 = (np.asarray(data_q1[2]) >= 920.0) & (np.asarray(data_q1[2]) <= 1250.0)
        mask_space3 = (np.asarray(data_q1[2]) >= 1250.0) & (np.asarray(data_q1[2]) <= 1530.0)
        mask_space4 = (np.asarray(data_q1[2]) >= 1530.0)
        ax.fill_between(data_q1[2][mask_space1], np.asarray(data_q2).T[1].T[0][mask_space1], np.asarray(data_q2).T[2].T[0][mask_space1],color='deepskyblue', alpha=.25)
        ax.fill_between(data_q1[2][mask_space2], np.asarray(data_q2).T[1].T[0][mask_space2], np.asarray(data_q2).T[2].T[0][mask_space2],color='deepskyblue', alpha=.25)
        ax.fill_between(data_q1[2][mask_space3], np.asarray(data_q2).T[1].T[0][mask_space3], np.asarray(data_q2).T[2].T[0][mask_space3],color='deepskyblue', alpha=.25)
        ax.fill_between(data_q1[2][mask_space4], np.asarray(data_q2).T[1].T[0][mask_space4], np.asarray(data_q2).T[2].T[0][mask_space4],color='deepskyblue', alpha=.25)

    ax.set_xlabel('pulser depth [m]')

    #ax.fill_between([800,900], 0, 30, color='grey',alpha=0.5)
    ax.set_ylim(0,30)
    ax.set_xlim(800,1700)

    y_data = spice_data
    y_mean = np.mean(y_data)
    y_sigma = np.std(y_data)


    xs2 = np.linspace(800.0,1700.0,1000)#depths_cut2
    y2 = func(spice_depths)
    fig2, ax2 = plt.subplots(1, 1,figsize=(5,5))
    gaussianHist(ax2,y_data-y2,'midnightblue',getMeanSTDStr(y_data-y2),[-4.35,160],'-','')
    ax2.set_xlabel(r'polarization [$\degree$]')
    ax2.set_ylabel(r'Number of events')
    fig2.tight_layout()



    ax3.errorbar(pol_data[0][1:],pol_data[1][1:],pol_errors[1:],fmt='-d',MarkerSize=10.0,color='black')
    ax3.set_xlabel(r'Launch angle [$^\circ$]')
    ax3.set_ylabel(r'Polarization [$^\circ$]')



    ax3.axis('tight')

    ax.text(0.99,0.98,'SPICE',horizontalalignment='right',verticalalignment='top',transform=ax.transAxes,color='midnightblue')
    ax.text(0.99,0.94,'Anechoic',horizontalalignment='right',verticalalignment='top',transform=ax.transAxes,color='darkorange')
    ax.text(0.99,0.90,'Comm. Period',horizontalalignment='right',verticalalignment='top',transform=ax.transAxes,color='grey')
    ax.annotate(s='', xy=(1008.0,28.2), xytext=(938,28.2),ha='center',va='top',color='midnightblue',arrowprops=dict(arrowstyle='->',color='midnightblue'))
    ax.text(0.27,0.98,s=r'R$\leq$0.5',horizontalalignment='right',verticalalignment='top',transform=ax.transAxes,color='midnightblue')
    ax.plot([938,938],[-5,35],color='midnightblue')

    ax.fill_between([shading_inputs[0],shading_inputs[1]],[60,60],color='black',alpha=0.15,linewidth=0.0)
    ax.fill_between([shading_inputs[0],shading_inputs[1]],[-17,-17],color='black',alpha=0.15,linewidth=0.0)
    ax.fill_between([shading_inputs[2],shading_inputs[3]],[60,60],color='black',alpha=0.15,linewidth=0.0)
    ax.fill_between([shading_inputs[2],shading_inputs[3]],[-17,-17],color='black',alpha=0.15,linewidth=0.0)
    ax.fill_between([shading_inputs[4],shading_inputs[5]],[60,60],color='black',alpha=0.15,linewidth=0.0)
    ax.fill_between([shading_inputs[4],shading_inputs[5]],[-17,-17],color='black',alpha=0.15,linewidth=0.0)


    plot_data_pol = np.asarray(plot_data_pol).T
    ax3.plot(plot_data_pol[0][1:],plot_data_pol[1][1:], 'b', lw=3,alpha=0.5,color='black',label='Anechoic')

    ax3.legend()
    ax3.fill_between(x_interp, y_interp_p, y_interp_m, color='green',alpha=0.5)


    fig.savefig(PathToARIANNAanalysis + '/plots/polReco.png')
    fig.savefig(PathToARIANNAanalysis + '/plots/polReco.pdf')
    fig2.savefig(PathToARIANNAanalysis + '/plots/polRecoHist.png')
    fig2.savefig(PathToARIANNAanalysis + '/plots/polRecoHist.pdf')
    fig3.savefig(PathToARIANNAanalysis + '/plots/polRecoExp.png')
    fig3.savefig(PathToARIANNAanalysis + '/plots/polRecoExp.pdf')



plt.show()
