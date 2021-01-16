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
import os

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
det = detector_sys_uncertainties.DetectorSysUncertainties(source='sql',assume_inf=False)  # establish mysql connection

PathToARIANNAanalysis = os.getcwd()

def printHeaderDetailsPerEvent(nurFile,channel_pairs):
	template = NuRadioRecoio.NuRadioRecoio(nurFile)
	direction_plot = []
	times = []
	for evt in template.get_events():
		for station_object in evt.get_stations():
			time = station_object.get_station_time()
			if station_object.has_triggered():
				det.update(station_object.get_station_time())
				l1s = []
				channelStopFilter.run(evt,station_object,det)
				#channelBandPassFilter.run(evt, station_object, det, passband=[80 * units.MHz, 300 * units.MHz], filter_type='butterabs',order=10)
				channelBandPassFilter.run(evt, station_object, det, passband=[80 * units.MHz, 300 * units.MHz], filter_type='rectangular')
				hardwareResponseIncorporator.run(evt, station_object, det, sim_to_data=False)
				channelSignalReconstructor.run(evt,station_object,det)
				channelResampler.run(evt, station_object, det, sampling_rate=50*units.GHz)
				cTW.run(evt, station_object, det, window_function='hanning')
				correclationDirectionFitter.run(evt,station_object,det,n_index=1.353,AziLim=[309 * units.deg, 315 * units.deg],ZenLim=[120 * units.deg, 160 * units.deg],channel_pairs=channel_pairs)  #AziLim=[309 * units.deg, 315 * units.deg]
				times.append(time)
				direction_plot.append([station_object.get_parameter(stnp.zenith)/units.deg,station_object.get_parameter(stnp.azimuth)/units.deg])
				channelResampler.run(evt, station_object, det, sampling_rate=1*units.GHz)
				print("Reconstructed Angular Directions: " + str([station_object.get_parameter(stnp.zenith)/units.deg,station_object.get_parameter(stnp.azimuth)/units.deg]))
	return times, direction_plot

def main():

	file = PathToARIANNAanalysis + '/data/Spice_750mDown_Dec30_2018_idl_10dB.nur'

	times, angles = printHeaderDetailsPerEvent(file,((0, 2), (1, 3)))
	np.save(PathToARIANNAanalysis + '/data/cc_lpdas_stn51',[times,angles])

	times, angles = printHeaderDetailsPerEvent(file,((4, 6), (5, 7)))
	np.save(PathToARIANNAanalysis + '/data/cc_dipoles_stn51',[times,angles])


if __name__== "__main__":
	main()
