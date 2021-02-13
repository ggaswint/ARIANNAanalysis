from NuRadioReco.modules.io import NuRadioRecoio
import numpy as np
import os
import datetime
PathToARIANNAanalysis = os.environ['ARIANNAanalysis']

def getAngularData(nurFile,channel_pairs):
    template = NuRadioRecoio.NuRadioRecoio(nurFile)
    direction_plot = []
    times = []
    for evt in template.get_events():
        for station in evt.get_stations():
            time = str(station.get_station_time())
            t1 = time.split('T')
            t11 = t1[0].split('-')
            t12 = t1[1].split(':')
            t22 = t12[2].split('.')

            tdt = datetime.datetime(int(t11[0]),int(t11[1]),int(t11[2]),int(t12[0]),int(t12[1]),int(t22[0]),int(t22[1])*1000)
            unix = (tdt - datetime.datetime(1970,1,1)).total_seconds()

            if station.has_triggered():
                times.append(unix)
    return times

def main():

    file = PathToARIANNAanalysis + '/data/Spice_750mDown_Dec30_2018_idl_10dB.nur'

    times = getAngularData(file,((0, 2), (1, 3)))
    np.save(PathToARIANNAanalysis + '/data/unixTimeStampsSpiceRunDec302018ARIANNAstn51TriggeredEvents',[times])

if __name__== "__main__":
    main()
