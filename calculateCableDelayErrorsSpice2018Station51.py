import os
import numpy as np
import ROOT
import subprocess
from optparse import OptionParser
from matplotlib import cm
import math
import matplotlib.pyplot  as plt
from scipy import signal
import logging
import copy
from NuRadioReco.detector import detector
from NuRadioReco.utilities import units
from NuRadioReco.modules.io.snowshovel import readARIANNADataCalib as CreadARIANNAData
from NuRadioReco.modules import channelResampler as CchannelResampler
from NuRadioReco.modules.ARIANNA import hardwareResponseIncorporator as ChardwareResponseIncorporator
from NuRadioReco.modules import correlationDirectionFitter as CcorrelationDirectionFitter
from NuRadioReco.modules import channelTimeWindow as cTWindow
import NuRadioReco.modules.voltageToEfieldConverter
import NuRadioReco.modules.channelSignalReconstructor
import NuRadioReco.modules.channelBandPassFilter
import NuRadioReco.modules.channelStopFilter
from NuRadioReco.framework.parameters import stationParameters as stnp
import logging
from scipy.interpolate import interp1d
from NuRadioReco.utilities import geometryUtilities as geo_utl
from NuRadioReco.utilities import ice
import helpfulUtils as hu
logger = logging.getLogger('plotDeconvolvedEvents')
logging.basicConfig(level=logging.DEBUG)

PathToARIANNAanalysis = os.environ['ARIANNAanalysis']

det = detector.Detector(source='sql',assume_inf=False)  # establish mysql connection

cTW = cTWindow.channelTimeWindow()
readARIANNAData = CreadARIANNAData.readARIANNAData()
channelResampler = CchannelResampler.channelResampler()
channelResampler.begin(debug=False)
channelBandPassFilter = NuRadioReco.modules.channelBandPassFilter.channelBandPassFilter()
hardwareResponseIncorporator = ChardwareResponseIncorporator.hardwareResponseIncorporator()
hardwareResponseIncorporator.begin(debug=False)
correlationDirectionFitter = CcorrelationDirectionFitter.correlationDirectionFitter()
channelSignalReconstructor = NuRadioReco.modules.channelSignalReconstructor.channelSignalReconstructor()
voltageToEfieldConverter = NuRadioReco.modules.voltageToEfieldConverter.voltageToEfieldConverter()
channelStopFilter = NuRadioReco.modules.channelStopFilter.channelStopFilter()
cTW.begin()

expData = np.load(PathToARIANNAanalysis + '/data/dec30_timeVSdepth.npy')
timesRun = expData[0]
depthsRun = expData[1]
f_data = interp1d(timesRun, depthsRun)

def getDepths30(time):
    depth = 0
    time = str(time)[11:]
    time_fixed = (int(time[:2]) + int(time[3:5])/60.0 + int(time[6:8])/3600.0)
    if time_fixed < 10:
        time_fixed += 24.0
    return -f_data(time_fixed)


def get_array_of_channels_voltages(station, use_channels, det, zenith, azimuth,
                          antenna_pattern_provider, time_domain=False):
    time_shifts = np.zeros(len(use_channels))
    station_id = station.get_id()
    for iCh, channel in enumerate(station.iter_channels()):
        channel_id = iCh

        antenna_position = det.get_relative_position(station_id, channel_id)
        # determine refractive index of signal propagation speed between antennas
        refractive_index = 1.353#ice.get_refractive_index(1, site)  # if signal comes from above, in-air propagation speed
        if(zenith > 0.5 * np.pi):
            refractive_index = 1.353#ice.get_refractive_index(antenna_position[2], site)  # if signal comes from below, use refractivity at antenna position
        time_shift = -geo_utl.get_time_delay_from_direction(zenith, azimuth, antenna_position, n=refractive_index)
        time_shift += channel.get_trace_start_time()
        time_shifts[iCh] = time_shift

    delta_t = time_shifts.max() - time_shifts.min()
    tmin = time_shifts.min()
    tmax = time_shifts.max()
    trace_length = station.get_channel(0).get_times()[-1] - station.get_channel(0).get_times()[0]

    traces = []
    n_samples = None
    for iCh, channel in enumerate(station.iter_channels()):
        tstart = delta_t - (time_shifts[iCh] - tmin)
        tstop = tmax - time_shifts[iCh] - delta_t + trace_length
        iStart = int(round(tstart * channel.get_sampling_rate()))
        iStop = int(round(tstop * channel.get_sampling_rate()))
        if(n_samples is None):
            n_samples = iStop - iStart
            if(n_samples % 2):
                n_samples -= 1

        trace = copy.copy(channel.get_trace())  # copy to not modify data structure
        trace = trace[iStart:(iStart + n_samples)]
        base_trace = NuRadioReco.framework.base_trace.BaseTrace()  # create base trace class to do the fft with correct normalization etc.
        base_trace.set_trace(trace, channel.get_sampling_rate())
        traces.append(base_trace)
    times = traces[0].get_times()  # assumes that all channels have the same sampling rate
    if(time_domain):  # save time domain traces first to avoid extra fft
        V_timedomain = np.zeros((len(use_channels), len(times)))
        for iCh, trace in enumerate(traces):
            V_timedomain[iCh] = trace.get_trace()

    return V_timedomain


def findTimes(file):

    ConsecutiveEntryNum = 0
    n_events = readARIANNAData.begin([file])

    event_count = 0
    times = [[],[],[],[],[],[],[],[]]
    depths = []
    for evt in readARIANNAData.run():
        for station in evt.get_stations():
            depth = float(getDepths30(station.get_station_time()))
            if station.has_triggered() and depth > 1180.98:
                det.update(station.get_station_time())
                channelStopFilter.run(evt,station,det)
                channelBandPassFilter.run(evt, station, det, passband=[80 * units.MHz, 300 * units.MHz], filter_type='rectangular')
                hardwareResponseIncorporator.run(evt, station, det, sim_to_data=False)
                channelSignalReconstructor.run(evt,station,det)
                channelResampler.run(evt, station, det, sampling_rate=50*units.GHz)
                cTW.run(evt, station, det, window_function='hanning')
                polarization = 0
                depths.append(float(depth))
                Zen = hu.findZenith(depth)[2]*units.deg
                Azi = 312.448284*units.deg
                station.set_parameter(stnp.zenith,Zen)
                station.set_parameter(stnp.azimuth,Azi)

                station_id = station.get_id()

                n_index = 1.353
                sampling_rate = station.get_channel(0).get_sampling_rate()
                chan = 6
                traces = get_array_of_channels_voltages(station, [0,1,2,3,4,5,6,7], det, Zen, Azi, voltageToEfieldConverter.antenna_provider, time_domain=True)

                mltp = [1,-1,-1,1,1,1,1,1]
                time0 = np.argmax(signal.correlate(traces[chan], traces[chan]))
                for iCh, trace in enumerate(traces):
                    time_c = signal.correlate(traces[chan], mltp[iCh]*trace)
                    times[iCh].append((time0-np.argmax(time_c))/sampling_rate)


                event_count+=1

                channelResampler.run(evt, station, det, sampling_rate=1*units.GHz)

    return depths, times

def main():


    file = PathToARIANNAanalysis + '/data/Spice_750mDown_Dec30_2018_idl_10dB.root'


    depths, times_per_channel = findTimes(file)
    np.save(PathToARIANNAanalysis + '/data/reCalculatedCableDelaysSpice2018.npy',[depths, times_per_channel])


if __name__== "__main__":
    main()
