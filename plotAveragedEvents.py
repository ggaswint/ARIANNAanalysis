import os
import numpy as np
import ROOT
import matplotlib.pyplot  as plt
import copy
from radiotools import helper as hp

from NuRadioReco.modules.io.snowshovel import readARIANNADataCalib as CreadARIANNAData

PathToARIANNAanalysis = os.environ['ARIANNAanalysis']

readARIANNAData = CreadARIANNAData.readARIANNAData()


color = ['C1','C2','C3','C4','C5','C6','C7','C8']

def plotting_tools(ax1,ax2,time,trace,chan_id,color):
    ax1.plot(time, trace*1e3, label=chan_id, color=color, linewidth = 1)
    delta_t = time[1] - time[0]
    freqs = np.fft.rfftfreq(len(trace), delta_t)
    ax2.plot(freqs, np.abs(np.fft.rfft(trace, norm='ortho')), color=color, linewidth = 1)
    ax2.set_xlim(0.0,0.5)
    if chan_id == 7:
        ax1.set_xlabel('[ns]')
        ax1.set_ylabel('[mV]')
        ax2.set_xlabel('[GHz]')
    else:
        ax1.xaxis.set_ticks([])
        ax2.xaxis.set_ticks([])

def correlate(t1,t2):
    xcorr = hp.get_normalized_xcorr(t1, t2)
    pos = np.argmax(xcorr)
    return xcorr, pos

def averageEvents(file):
    n_events = readARIANNAData.begin([file])

    event_count = 0
    traces = []
    traces_first_event = []
    for evt in readARIANNAData.run():
        for station in evt.get_stations():
            #################################################################################
            if event_count == 0:
                times = copy.copy(station.get_channel(0).get_times())
                for iCh, chan in enumerate(station.iter_channels()):
                    traces.append(chan.get_trace())
                    traces_first_event.append(chan.get_trace())

            else:
                for iCh, chan in enumerate(station.iter_channels()):
                    xcorr, pos = correlate(traces_first_event[iCh], chan.get_trace())
                    traces[iCh] += np.roll(chan.get_trace(), pos + len(traces_first_event[iCh]) + 1)

            event_count += 1

    traces = np.asarray(traces)/event_count
    return traces,times

def main():

    file = PathToARIANNAanalysis + '/data/Spice_750mDown_Dec30_2018_idl_10dB.root'

    traces, times = averageEvents(file)

    f, ax = plt.subplots(8,2,figsize=(16, 8))
    for iCh in range(len(traces)):
        plotting_tools(ax[iCh][0],ax[iCh][1],times,traces[iCh],iCh,color=color[iCh])

    f.tight_layout()
    f.savefig(PathToARIANNAanalysis + '/plots/averagedTraces.png')
    f.savefig(PathToARIANNAanalysis + '/plots/averagedTraces.pdf')
    plt.show()

if __name__== "__main__":
    main()
