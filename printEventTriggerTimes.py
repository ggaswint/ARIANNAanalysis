from NuRadioReco.detector import detector
from NuRadioReco.modules.io.snowshovel import readARIANNADataCalib as CreadARIANNAData
import os

PathToARIANNAanalysis = os.environ['ARIANNAanalysis']
readARIANNAData = CreadARIANNAData.readARIANNAData()
det = detector.Detector(source='sql',assume_inf=False)  # establish mysql connection

def printTimes(file):
	n_events = readARIANNAData.begin([file])
	event_count = 0
	for evt in readARIANNAData.run():
		for station_object in evt.get_stations():
			print_str = 'event: ' + str(event_count) + ' time: ' + str(station_object.get_station_time())
			print(print_str)
			event_count += 1

def main():
    file = PathToARIANNAanalysis + '/data/Spice_750mDown_Dec30_2018_idl_10dB.root'
    printTimes(file)

if __name__== "__main__":
	main()
