from NuRadioReco.detector import detector
from NuRadioReco.modules.io import NuRadioRecoio
import os
import matplotlib.pyplot as plt
import numpy as np
PathToARIANNAanalysis = os.environ['ARIANNAanalysis']

det = detector.Detector(source='sql',assume_inf=False)

def plotRawTrace(nurFile):
    template = NuRadioRecoio.NuRadioRecoio(nurFile)
    for evt in template.get_events():
        for station in evt.get_stations():
            if station.has_triggered():
                det.update(station.get_station_time())

                fig1, ax1 = plt.subplots(8, 2, sharex=False, figsize=(12, 20))
                for iCh, channel in enumerate(station.iter_channels()):
                    ax1[iCh][0].plot(channel.get_times(), channel.get_trace(),linewidth=2)
                    delta_t = channel.get_times()[1] - channel.get_times()[0]
                    freqs = np.fft.rfftfreq(len(channel.get_trace()), delta_t)
                    ax1[iCh][1].plot(freqs, np.abs(np.fft.rfft(channel.get_trace(), norm='ortho')),linewidth=2)
                fig1.tight_layout()
                plt.show()


def main():
    file = PathToARIANNAanalysis + '/data/Spice_750mDown_Dec30_2018_idl_10dB.nur'
    plotRawTrace(file)


if __name__== "__main__":
    main()
