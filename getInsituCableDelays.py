def getAdditionalInsituCableDelaysFromSpice():
	# First 4 time offsets are lpdas: Spice750mDownData but only uses 938mDownData(reflection interferences are less dominate) with smooth ice profile, rectangular filter 80-300, with NuRadio default time window
	# Last 4 time offsets are dipoles: Spice750mDownData but only uses 1180mDownData(reflection interferences are less dominate) with smooth ice profile, rectangular filter 80-300, with NuRadio default time window
	time_offsets = [-1.357185244587009,-0.7436246992782678,-0.1495028067361668,0.0,0.216011258544431,-0.1058818737270876,0.0,-0.8744174557431041]
	return time_offsets

def main():
	print('Time delays from SPICE 2018 are: ' + str(getAdditionalInsituCableDelaysFromSpice()))

if __name__== "__main__":
	main()
