from NuRadioReco.detector import detector
from NuRadioReco.detector import detector_sys_uncertainties
from NuRadioReco.utilities import units
from NuRadioReco.modules import channelResampler as CchannelResampler
from NuRadioReco.modules.ARIANNA import hardwareResponseIncorporator as ChardwareResponseIncorporator
from NuRadioReco.modules import channelTimeWindow as cTWindow
import NuRadioReco.modules.channelSignalReconstructor
import NuRadioReco.modules.channelBandPassFilter
import NuRadioReco.modules.channelStopFilter
import NuRadioReco.modules.correlationDirectionFitter
from NuRadioReco.framework.parameters import stationParameters as stnp
from NuRadioReco.modules.io import NuRadioRecoio
import numpy as np
from scipy.interpolate import interp1d
import os
PathToARIANNAanalysis = os.environ['ARIANNAanalysis']

import logging
logger = logging.getLogger('plotDeconvolvedEvents')
logging.basicConfig(level=logging.WARNING)

channelResampler = CchannelResampler.channelResampler()
channelResampler.begin(debug=False)
channelBandPassFilter = NuRadioReco.modules.channelBandPassFilter.channelBandPassFilter()
hardwareResponseIncorporator = ChardwareResponseIncorporator.hardwareResponseIncorporator()
hardwareResponseIncorporator.begin(debug=False)
channelSignalReconstructor = NuRadioReco.modules.channelSignalReconstructor.channelSignalReconstructor()
correclationDirectionFitter = NuRadioReco.modules.correlationDirectionFitter.correlationDirectionFitter()
correclationDirectionFitter.begin(debug=False)
channelStopFilter = NuRadioReco.modules.channelStopFilter.channelStopFilter()
cTW = cTWindow.channelTimeWindow()
cTW.begin(debug=False)

# The SysUncertainties version simply allows one to change the orientation and position of the antennas on the fly. Only useful for studying systematic effects due to uncertainty in antennas.
# May have issues with certain modules so if you ever run into a problem, try switching to the regular detector
#det = detector_sys_uncertainties.DetectorSysUncertainties(source='sql',assume_inf=False)
det = detector.Detector(source='sql',assume_inf=False)

expData = np.load(PathToARIANNAanalysis + '/data/dec30_timeVSdepth.npy')
timesRun = expData[0]
depthsRun = expData[1]
f_data = interp1d(timesRun, depthsRun)

def getDepthsDec30SPICErun(time):
    # convert time to a hour unit for comparing with SPICE run time stamps
	time = str(time)[11:]
	time_fixed = (int(time[:2]) + int(time[3:5])/60.0 + int(time[6:8])/3600.0)
	if time_fixed < 10: # Times on the 31st of December are written as hours greater than 24
		time_fixed += 24.0
	return -float(f_data(time_fixed))

def getAngularData(nurFile,channel_pairs):
    template = NuRadioRecoio.NuRadioRecoio(nurFile)
    direction_plot = []
    times = []
    depths = []
    for evt in template.get_events():
        for station in evt.get_stations():
            time = station.get_station_time()
            if station.has_triggered():
                det.update(station.get_station_time())
                depths.append(float(getDepthsDec30SPICErun(station.get_station_time())))
                station.set_is_neutrino()
                l1s = []
                channelStopFilter.run(evt,station,det)
                #channelBandPassFilter.run(evt, station, det, passband=[80 * units.MHz, 300 * units.MHz], filter_type='butterabs',order=10)
                channelBandPassFilter.run(evt, station, det, passband=[80 * units.MHz, 300 * units.MHz], filter_type='rectangular')
                hardwareResponseIncorporator.run(evt, station, det, sim_to_data=False)
                channelSignalReconstructor.run(evt,station,det)
                channelResampler.run(evt, station, det, sampling_rate=50*units.GHz)
                cTW.run(evt, station, det, window_function='hanning')
                # Uncomment below to add more precision to the cable delays from the SPICE data
				#time_offsets = [-1.3367996742671011,-0.7025223759153785,-0.1587785016286645,0.0,0.22120196238757153,-0.06695943120033458,0.0,-0.935954487989886] # calculated additional cable delays found from SPICE data using getInsituCableDelays.py
				#for channel in station.iter_channels():
				#   channel.set_trace_start_time(channel.get_trace_start_time()-time_offsets[channel.get_id()])
                correclationDirectionFitter.run(evt,station,det,n_index=1.353,AziLim=[309 * units.deg, 315 * units.deg],ZenLim=[120 * units.deg, 160 * units.deg],channel_pairs=channel_pairs)  #AziLim=[309 * units.deg, 315 * units.deg]
                times.append(time)
                direction_plot.append([station.get_parameter(stnp.zenith)/units.deg,station.get_parameter(stnp.azimuth)/units.deg])
                channelResampler.run(evt, station, det, sampling_rate=1*units.GHz)
                print("Reconstructed Angular Directions: " + str([station.get_parameter(stnp.zenith)/units.deg,station.get_parameter(stnp.azimuth)/units.deg]))
    return depths, direction_plot, times

def main():

    file = PathToARIANNAanalysis + '/data/Spice_750mDown_Dec30_2018_idl_10dB.nur'

    depths, direction_plot, times = getAngularData(file,((0, 2), (1, 3)))
    np.save(PathToARIANNAanalysis + '/data/cc_lpdas_stn51',[depths, direction_plot, times])

    depths, direction_plot, times = getAngularData(file,((4, 6), (5, 7)))
    np.save(PathToARIANNAanalysis + '/data/cc_dipoles_stn51',[depths, direction_plot, times])


if __name__== "__main__":
    main()
