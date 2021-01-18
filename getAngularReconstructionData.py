from NuRadioReco.detector import detector
from NuRadioReco.detector import detector_sys_uncertainties
from NuRadioReco.utilities import units
from NuRadioReco.modules.io.snowshovel import readARIANNADataCalib as CreadARIANNAData
from NuRadioReco.modules import channelResampler as CchannelResampler
from NuRadioReco.modules.ARIANNA import hardwareResponseIncorporator as ChardwareResponseIncorporator
from NuRadioReco.modules import channelTimeWindow as cTWindow
import NuRadioReco.modules.channelSignalReconstructor
import NuRadioReco.modules.channelBandPassFilter
import NuRadioReco.modules.channelStopFilter
import NuRadioReco.modules.correlationDirectionFitter
from NuRadioReco.framework.parameters import stationParameters as stnp
import numpy as np
import os
PathToARIANNAanalysis = os.environ['ARIANNAanalysis']

import logging
logger = logging.getLogger('plotDeconvolvedEvents')
logging.basicConfig(level=logging.WARNING)

readARIANNAData = CreadARIANNAData.readARIANNAData()
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

def getAngularData(file,channel_pairs):
    n_events = readARIANNAData.begin([file])
    direction_plot = []
    times = []
    for evt in readARIANNAData.run():
        for station in evt.get_stations():
            time = station.get_station_time()
            if station.has_triggered():
                det.update(station.get_station_time())
                station.set_is_neutrino()
                l1s = []
                channelStopFilter.run(evt,station,det)
                #channelBandPassFilter.run(evt, station, det, passband=[80 * units.MHz, 300 * units.MHz], filter_type='butterabs',order=10)
                channelBandPassFilter.run(evt, station, det, passband=[80 * units.MHz, 300 * units.MHz], filter_type='rectangular')
                hardwareResponseIncorporator.run(evt, station, det, sim_to_data=False)
                channelSignalReconstructor.run(evt,station,det)
                channelResampler.run(evt, station, det, sampling_rate=50*units.GHz)
                cTW.run(evt, station, det, window_function='hanning')
                correclationDirectionFitter.run(evt,station,det,n_index=1.353,AziLim=[309 * units.deg, 315 * units.deg],ZenLim=[120 * units.deg, 160 * units.deg],channel_pairs=channel_pairs)  #AziLim=[309 * units.deg, 315 * units.deg]
                times.append(time)
                direction_plot.append([station.get_parameter(stnp.zenith)/units.deg,station.get_parameter(stnp.azimuth)/units.deg])
                channelResampler.run(evt, station, det, sampling_rate=1*units.GHz)
                print("Reconstructed Angular Directions: " + str([station.get_parameter(stnp.zenith)/units.deg,station.get_parameter(stnp.azimuth)/units.deg]))
    return times, direction_plot

def main():

    file = PathToARIANNAanalysis + '/data/Spice_750mDown_Dec30_2018_idl_10dB.root'

    times, angles = getAngularData(file,((0, 2), (1, 3)))
    np.save(PathToARIANNAanalysis + '/data/cc_lpdas_stn51_rootInput',[times,angles])

    times, angles = getAngularData(file,((4, 6), (5, 7)))
    np.save(PathToARIANNAanalysis + '/data/cc_dipoles_stn51_rootInput',[times,angles])


if __name__== "__main__":
    main()
