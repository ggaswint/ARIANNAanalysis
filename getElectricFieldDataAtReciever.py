import os
import numpy as np
import math
import matplotlib.pyplot  as plt
from scipy import signal
from scipy.interpolate import interp1d
from scipy.signal import hilbert
from radiotools import helper as hp
import getAnglesFromSPICE as getAngles # part of ARIANNAnalysis
import getInsituCableDelays as getTimeDelays # part of ARIANNAnalysis

from NuRadioReco.utilities import units
from NuRadioReco.modules import channelResampler as CchannelResampler
from NuRadioReco.modules.ARIANNA import hardwareResponseIncorporator as ChardwareResponseIncorporator
from NuRadioReco.modules import channelTimeWindow as cTWindow
import NuRadioReco.modules.voltageToEfieldConverter
import NuRadioReco.modules.channelSignalReconstructor
import NuRadioReco.modules.channelBandPassFilter
import NuRadioReco.modules.channelStopFilter
from NuRadioReco.framework.parameters import stationParameters as stnp
import NuRadioReco.modules.electricFieldBandPassFilter
from NuRadioReco.modules.io import NuRadioRecoio
from NuRadioReco.detector import detector_sys_uncertainties
from NuRadioReco.detector import detector

import logging
logger = logging.getLogger('GetElectricFieldDataAtReciever')
logging.basicConfig(level=logging.WARNING)
logger.disabled = True

PathToARIANNAanalysis = os.environ['ARIANNAanalysis']

channelResampler = CchannelResampler.channelResampler()
channelResampler.begin(debug=False)
channelBandPassFilter = NuRadioReco.modules.channelBandPassFilter.channelBandPassFilter()
hardwareResponseIncorporator = ChardwareResponseIncorporator.hardwareResponseIncorporator()
hardwareResponseIncorporator.begin(debug=False)
channelSignalReconstructor = NuRadioReco.modules.channelSignalReconstructor.channelSignalReconstructor()
voltageToEfieldConverter = NuRadioReco.modules.voltageToEfieldConverter.voltageToEfieldConverter()
electricFieldBandPassFilter = NuRadioReco.modules.electricFieldBandPassFilter.electricFieldBandPassFilter()
channelStopFilter = NuRadioReco.modules.channelStopFilter.channelStopFilter()
cTW = cTWindow.channelTimeWindow()
cTW.begin(debug=False)

#det = detector_sys_uncertainties.DetectorSysUncertainties(source='sql',assume_inf=False)  # establish mysql connection
# An example of using a json detector as input. Can convert sql to jason using NuRadioReco in the detector folder there is a sql_to_json python script
det = detector.Detector(source='json',json_filename=PathToARIANNAanalysis + '/data/arianna_detector_db.json',assume_inf=False)

expData = np.load(PathToARIANNAanalysis + '/data/SPICE_dec30_timeVSdepth.npy')
timesRun = expData[0]
depthsRun = expData[1]
f_data = interp1d(timesRun, depthsRun)

#converts time stamp to fractional hour so can be used to acces what the depth of the pulser was per second.
def getDepthsFromTimes(time):
    depth = 0
    time = str(time)[11:]
    time_fixed = (int(time[:2]) + int(time[3:5])/60.0 + int(time[6:8])/3600.0)
    if time_fixed < 10:
        time_fixed += 24.0
    return -float(f_data(time_fixed))

def findElectricFieldProperties(nurFile, chans, force_Polarization=False):
    template = NuRadioRecoio.NuRadioRecoio(nurFile)
    event_count = 0
    polData = []
    depths = []
    ampsT = []
    ampsP = []
    EtimeT = []
    EtimeP = []
    lower = 80 # Filter MHz
    upper = 300 # Filter MHz
    for evt in template.get_events():
        for station in evt.get_stations():
            depth = getDepthsFromTimes(station.get_station_time())
            if station.has_triggered() and depth > 800.0:# and bad_event not in [721,722,725,2002,2143,2367]:
                det.update(station.get_station_time())
                #sim_station = NuRadioReco.framework.sim_station.SimStation(station.get_id())
                #sim_station._particle_type = 'test'
                #station.set_sim_station(sim_station)
                print(station.get_sim_station())
                station.set_is_neutrino()
                channelStopFilter.run(evt,station,det)
                channelBandPassFilter.run(evt, station, det, passband=[lower * units.MHz, upper * units.MHz], filter_type='rectangular')
                hardwareResponseIncorporator.run(evt, station, det, sim_to_data=False)
                channelSignalReconstructor.run(evt,station,det)
                channelResampler.run(evt, station, det, sampling_rate=50*units.GHz)
                Zen, Azi = np.asarray(getAngles.findZenAziFromSpiceDataLPDAs(depth)[3])*units.deg
                station.set_parameter(stnp.zenith,Zen)
                station.set_parameter(stnp.azimuth,Azi)
                sampling_rate = station.get_channel(0).get_sampling_rate()
                #cTW.run(evt, station, det, window_function='hanning')

                # time delays between channels
                time_offsets = getTimeDelays.getAdditionalInsituCableDelaysFromSpice()
                for channel in station.iter_channels():
                    channel.set_trace_start_time(channel.get_trace_start_time()-time_offsets[channel.get_id()])


                voltageToEfieldConverter.run(evt, station, det, use_channels=chans, force_Polarization=force_Polarization)
                electricFieldBandPassFilter.run(evt, station, det, passband=[lower * units.MHz, upper * units.MHz], filter_type='rectangular')

                etheta = station.get_electric_fields()[0].get_trace()[1]
                ephi = station.get_electric_fields()[0].get_trace()[2]
                e_times = station.get_electric_fields()[0].get_times()
                h_etheta = hilbert(station.get_electric_fields()[0].get_trace()[1])
                h_ephi = hilbert(station.get_electric_fields()[0].get_trace()[2])
                h3 = np.sqrt(np.abs(h_etheta)**2 + np.abs(h_ephi)**2)
                fwhm = hp.get_FWHM_hilbert(h3)

                IW = int(sampling_rate*70.0) # length of window
                mid_fwhm = fwhm[0] + int((fwhm[1] - fwhm[0])/2) # Center of FWHM
                noise_idx = int(1.1*int(mid_fwhm+IW/2)) # Noise start
                signal_idx = int(mid_fwhm+IW/2) # signal end

                max_etheta = np.sqrt(np.abs(np.sum((np.abs(etheta[signal_idx-IW:signal_idx]))**2) - np.sum((np.abs(etheta[noise_idx:noise_idx+IW]))**2)))
                max_ephi = np.sqrt(np.abs(np.sum((np.abs(ephi[signal_idx-IW:signal_idx]))**2) - np.sum((np.abs(ephi[noise_idx:noise_idx+IW]))**2)))

                time_idx_theta = np.argmax(np.abs(h_etheta))
                time_idx_phi = np.argmax(np.abs(h_ephi))

                depths.append(depth)
                polData.append(max_ephi/max_etheta)
                ampsT.append(np.max(etheta))
                ampsP.append(np.max(ephi))
                EtimeT.append(e_times[time_idx_theta])
                EtimeP.append(e_times[time_idx_phi])


                channelResampler.run(evt, station, det, sampling_rate=1*units.GHz)
            event_count += 1

    return depths, polData, ampsT, ampsP, EtimeT, EtimeP

def main():
    print('Please replace file from below to the data you would like to analyze')
    file = PathToARIANNAanalysis + '/data/Spice_750mDown_Dec30_2018_idl_10dB.nur'
    print(file)

    depths, polData, ampsT, ampsP, EtimeT, EtimeP = findElectricFieldProperties(file, [0,1,2,3])
    np.save(PathToARIANNAanalysis + '/data/SPICE_Dec30_2018_EFieldDataLPDA',[depths, polData, ampsT, ampsP, EtimeT, EtimeP])

    depths, polData, ampsT, ampsP, EtimeT, EtimeP = findElectricFieldProperties(file, [4,5,6,7],'eTheta')
    np.save(PathToARIANNAanalysis + '/data/SPICE_Dec30_2018_EFieldDataDipole',[depths, polData, ampsT, ampsP, EtimeT, EtimeP])


if __name__== "__main__":
    main()
