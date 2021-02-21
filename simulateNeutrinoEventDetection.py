import argparse
# import detector simulation modules
import NuRadioReco.modules.trigger.highLowThreshold
import NuRadioReco.modules.channelBandPassFilter
import NuRadioReco.modules.triggerTimeAdjuster
from NuRadioReco.utilities import units
from NuRadioMC.simulation import simulation

# initialize detector sim modules
highLowThreshold = NuRadioReco.modules.trigger.highLowThreshold.triggerSimulator()
channelBandPassFilter = NuRadioReco.modules.channelBandPassFilter.channelBandPassFilter()

triggerTimeAdjuster = NuRadioReco.modules.triggerTimeAdjuster.triggerTimeAdjuster()
triggerTimeAdjuster.begin(pre_trigger_time=155.*units.ns)

thresholds = {
  '2/4_100Hz': 3.9498194908011524,
  '2/4_10mHz': 4.919151494949084,
  '2/6_100Hz': 4.04625348733533,
  '2/6_10mHz': 5.015585491483261,
  'fhigh': 0.15,
  'flow': 0.08
  }

passband_low = {}
passband_high = {}
for channel_id in range(0, 5):
	passband_low[channel_id] = [1 * units.MHz, 500 * units.MHz]
	passband_high[channel_id] = [80 * units.MHz, 1000 * units.GHz]
for channel_id in range(5, 9):
	passband_low[channel_id] = [1 * units.MHz, thresholds['fhigh']]
	passband_high[channel_id] = [thresholds['flow'], 800 * units.GHz]

passband_low[9] = [1 * units.MHz, 730 * units.MHz]
passband_high[9] = [100 * units.MHz, 800 * units.GHz]


class mySimulation(simulation.simulation):

	def _detector_simulation_filter_amp(self, evt, station, det):
		channelBandPassFilter.run(evt, station, det,
								  passband=passband_low, filter_type="butter", order=10)
		channelBandPassFilter.run(evt, station, det,
								  passband=passband_high, filter_type="butter", order=5)

	def _detector_simulation_trigger(self, evt, station, det):
		# run a high/low trigger on the 4 downward pointing LPDAs
		threshold_high = {}
		threshold_low = {}
		for channel_id in det.get_channel_ids(station.get_id()):
			threshold_high[channel_id] = 3.9 * self._Vrms_per_channel[station.get_id()][channel_id]
			threshold_low[channel_id] = -3.9 * self._Vrms_per_channel[station.get_id()][channel_id]

		highLowThreshold.run(evt, station, det,
									threshold_high=threshold_high,
									threshold_low=threshold_low,
									coinc_window=40 * units.ns,
									triggered_channels=[5, 6, 7, 8],  # select the LPDA channels
									number_concidences=2,  # 2/4 majority logic
									trigger_name='LPDA_2of4_3.9sigma')
		triggerTimeAdjuster.run(evt, station, det)


parser = argparse.ArgumentParser(description='Run NuRadioMC simulation')
parser.add_argument('inputfilename', type=str,
					help='path to NuRadioMC input event list')
parser.add_argument('detectordescription', type=str,
					help='path to file containing the detector description')
parser.add_argument('config', type=str,
					help='NuRadioMC yaml config file')
parser.add_argument('outputfilename', type=str,
					help='hdf5 output filename')
parser.add_argument('outputfilenameNuRadioReco', type=str, nargs='?', default=None,
					help='outputfilename of NuRadioReco detector sim file')
args = parser.parse_args()

sim = mySimulation(inputfilename=args.inputfilename,
							outputfilename=args.outputfilename,
							detectorfile=args.detectordescription,
							outputfilenameNuRadioReco=args.outputfilenameNuRadioReco,
							config_file=args.config,
							default_detector_station=1)
sim.run()

# Example for running this script:
# python simulateNeutrinoEventDetection.py data/tenNeutrinosAt1e19.hdf5 configurations/ARIANNA_4LPDA_1dipole.json configurations/simulateNeutrinosConfig.yaml data/triggeredNeutrinoEvents.hdf5 data/triggeredNeutrinoEvents.nur
