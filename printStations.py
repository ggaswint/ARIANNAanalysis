from NuRadioReco.detector import detector
from astropy.time import Time
from NuRadioReco.utilities import units

det = detector.Detector(source='sql',assume_inf=False)  # establish mysql connection

def printStations():
    det.update(Time("2018-12-30 23:26:58"))
    station = 51
    channel = 0
    det.get_number_of_channels(station)
    det.get_antenna_type(station,channel)
    #    get anetnna orientation prints the following four values as a list
    #    * orientation theta: boresight direction (zenith angle, 0deg is the zenith, 180deg is straight down)
    #    * orientation phi: boresight direction (azimuth angle counting from East counterclockwise)
    #    * rotation theta: rotation of the antenna, is perpendicular to 'orientation', for LPDAs: vector in plane of tines pointing away from connector
    #    * rotation phi: rotation of the antenna, is perpendicular to 'orientation', for LPDAs: vector in plane of tines pointing away from connector
    # Note if first argument is 180 or pi, then antenna is pointing down, last argument is the horizontal rotation also known as azimuth or phi
    det.get_antenna_orientation(station,channel)

    # stations 51 and 61 are at the SouthPole, the rest are in Moore's Bay
    stations = [14,15,17,18,19,30,32,51,52,61] # excludes 50 because it does not have station details for the time given above
    for station in stations:
        print('station: ' + str(station))
        for channel in range(det.get_number_of_channels(station)):
            if abs(det.get_antenna_orientation(station,channel)[0]/units.deg - 180) < 10.0:
                print('channel: ' + str(channel) + ' type: ' + str(det.get_antenna_type(station,channel)) +' antenna is pointing down')
            else:
                print('channel: ' + str(channel) + ' type: ' + str(det.get_antenna_type(station,channel)))
            print('orientation: ' + str(det.get_antenna_orientation(station,channel)/units.deg) + 'degrees')
        print('\n')



def main():
    printStations()

if __name__== "__main__":
    main()
