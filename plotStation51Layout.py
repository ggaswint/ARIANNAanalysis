import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
from matplotlib.ticker import NullFormatter
import numpy as np
import os

PathToARIANNAanalysis = os.environ['ARIANNAanalysis']

# uses 0,0 bottom left corner and 1,1 top right corner

#https://www.researchgate.net/publication/256538503_South_Pole_Glacial_Climate_Reconstruction_from_Multi-Borehole_Laser_Particulate_Stratigraphy
ice_flow = 320 # Degrees
Direction_from_Spice = 312.448284 # degrees
# 0.05 = 1m
scale = 0.05
side_len = 6 #m sides of square for station
dipole_channel_text_spacing = 0.03
lpda_channel_text_spacing = 0.05
meter_scale_location_below_square = 0.15
bottom_left_square = [0.1,0.2]
axis_text_spacing = 0.05
lpda_length = 2 #m
dipole_color = 'blue'
lpda_color = 'red'
square_side = side_len *scale
axis_arrow_len = square_side/1.5


ice_flow = 40.0
Direction_to_Spice = Direction_from_Spice + 180 + 90 - 360
print('angle between ice flow and spice is: ' + str(ice_flow-Direction_to_Spice+180))

bottom_right_square = [bottom_left_square[0] + square_side, bottom_left_square[1]]
top_left_square = [bottom_left_square[0], bottom_left_square[1] + square_side]
top_right_square = [bottom_left_square[0] + square_side, bottom_left_square[1] + square_side]
bottom_middle = (bottom_right_square[0] + bottom_left_square[0])/2.0
top_middle = (top_right_square[0] + top_left_square[0])/2.0
right_middle = (top_right_square[1] + bottom_right_square[1])/2.0
left_middle = (top_left_square[1] + bottom_left_square[1])/2.0
start_meter = [bottom_left_square[0] , bottom_left_square[1] - meter_scale_location_below_square]
axis_origin = [bottom_right_square[0] + square_side,right_middle]


fig, ax = plt.subplots(1, 1,figsize=(9, 9))

### Dipoles
plt.plot(top_left_square[0],top_left_square[1],'o',color=dipole_color)
plt.text(top_left_square[0] - dipole_channel_text_spacing,top_left_square[1] + dipole_channel_text_spacing,'Ch4',rotation=45,ha='center',va='center')
plt.plot(top_right_square[0],top_right_square[1],'o',color=dipole_color)
plt.text(top_right_square[0] + dipole_channel_text_spacing,top_right_square[1] + dipole_channel_text_spacing,'Ch5',rotation=-45,ha='center',va='center')
plt.plot(bottom_left_square[0],bottom_left_square[1],'o',color=dipole_color)
plt.text(bottom_left_square[0] - dipole_channel_text_spacing,bottom_left_square[1] - dipole_channel_text_spacing,'Ch7',rotation=-45,ha='center',va='center')
plt.plot(bottom_right_square[0],bottom_right_square[1],'o',color=dipole_color)
plt.text(bottom_right_square[0] + dipole_channel_text_spacing,bottom_right_square[1] - dipole_channel_text_spacing,'Ch6',rotation=45,ha='center',va='center')

# Lpdas
plt.annotate(s='', xy=(top_middle + lpda_length*scale/2.0,top_left_square[1]), xytext=(top_middle - lpda_length*scale/2.0,top_left_square[1]), arrowprops=dict(arrowstyle='-',color=lpda_color))
plt.text(top_middle,top_left_square[1] + lpda_channel_text_spacing,'Ch0',rotation=0,ha='center',va='center')
plt.annotate(s='', xy=(bottom_right_square[0],right_middle + lpda_length*scale/2.0), xytext=(bottom_right_square[0],right_middle - lpda_length*scale/2.0), arrowprops=dict(arrowstyle='-',color=lpda_color))
plt.text(bottom_right_square[0] + lpda_channel_text_spacing/1.5,right_middle,'Ch1',rotation=-90,ha='center',va='center')
plt.annotate(s='', xy=(bottom_middle + lpda_length*scale/2.0,bottom_left_square[1]), xytext=(bottom_middle - lpda_length*scale/2.0,bottom_left_square[1]), arrowprops=dict(arrowstyle='-',color=lpda_color))
plt.text(bottom_middle,bottom_left_square[1] - lpda_channel_text_spacing,'Ch2',rotation=0,ha='center',va='center')
plt.annotate(s='', xy=(bottom_left_square[0],left_middle + lpda_length*scale/2.0), xytext=(bottom_left_square[0],left_middle - lpda_length*scale/2.0), arrowprops=dict(arrowstyle='-',color=lpda_color))
plt.text(bottom_left_square[0] - lpda_channel_text_spacing/1.5,left_middle,'Ch3',rotation=90,ha='center',va='center')

#label LPDA downfacing
plt.annotate(s='Downward facing LPDA', xy=(top_middle,top_left_square[1] + lpda_channel_text_spacing*2), xytext=(top_middle,top_left_square[1] + lpda_channel_text_spacing*4),ha='center',va='center', arrowprops=dict(arrowstyle='->'))

#label dipole
plt.annotate(s='Dipole', xy=(top_right_square[0] + dipole_channel_text_spacing*2,top_right_square[1] + dipole_channel_text_spacing*2), xytext=(top_right_square[0] + dipole_channel_text_spacing*6,top_right_square[1] + dipole_channel_text_spacing*6),ha='center',va='center', arrowprops=dict(arrowstyle='->'))


# meter scale
plt.annotate(s='', xy=(start_meter[0] + side_len*scale,start_meter[1]), xytext=(start_meter[0],start_meter[1]), arrowprops=dict(arrowstyle='<->'))
plt.annotate(s='6m', xy=(start_meter[0] + (side_len/2.0)*scale,start_meter[1] - 0.04),ha='center')


# Axis labels
plt.annotate(s='', xy=(axis_origin[0] + axis_arrow_len,axis_origin[1]), xytext=(axis_origin[0],axis_origin[1]),ha='center',va='center', arrowprops=dict(arrowstyle='->'))
plt.text(axis_origin[0] + axis_arrow_len,axis_origin[1]-axis_text_spacing,'Grid North',rotation=0,ha='center',va='center')
plt.annotate(s='', xy=(axis_origin[0],axis_origin[1] - axis_arrow_len), xytext=(axis_origin[0],axis_origin[1]),ha='center',va='center', arrowprops=dict(arrowstyle='->'))
plt.text(axis_origin[0] + axis_text_spacing,axis_origin[1] - axis_arrow_len -axis_text_spacing,'South Pole Station',rotation=0,ha='center',va='center')
plt.annotate(s='', xy=(axis_origin[0] + axis_arrow_len*np.cos(np.deg2rad(ice_flow)),axis_origin[1] + axis_arrow_len*np.sin(np.deg2rad(ice_flow))), xytext=(axis_origin[0],axis_origin[1]),ha='center',va='center', arrowprops=dict(arrowstyle='->',color='brown'))
plt.text(axis_origin[0] + axis_arrow_len*np.cos(np.deg2rad(ice_flow)) - axis_text_spacing,axis_origin[1] + axis_arrow_len*np.sin(np.deg2rad(ice_flow)),'Ice Flow',rotation=ice_flow,ha='center',va='center',color='brown')
plt.annotate(s='', xy=(axis_origin[0] + axis_arrow_len*np.cos(np.deg2rad(Direction_to_Spice)),axis_origin[1] + axis_arrow_len*np.sin(np.deg2rad(Direction_to_Spice))), xytext=(axis_origin[0],axis_origin[1]),ha='center',va='center', arrowprops=dict(arrowstyle='->',color='purple'))
plt.text(axis_origin[0] + axis_arrow_len*np.cos(np.deg2rad(Direction_to_Spice)),axis_origin[1] + axis_arrow_len*np.sin(np.deg2rad(Direction_to_Spice)) + axis_text_spacing,'SPICE',rotation=(Direction_to_Spice-180-7),ha='center',va='center',color='purple')

extra_spacing_ice_flow_origion = [axis_origin[0]+axis_arrow_len*np.cos(np.deg2rad(ice_flow))/2.0,axis_origin[1]+axis_arrow_len*np.sin(np.deg2rad(ice_flow))/2.0]
extra_spacing_spice_origion = [axis_origin[0]+axis_arrow_len*np.cos(np.deg2rad(Direction_to_Spice))/2.0,axis_origin[1]+axis_arrow_len*np.sin(np.deg2rad(Direction_to_Spice))/2.0]


plt.xlim(0,1)
plt.ylim(0,1)


plt.grid(False)
ax.xaxis.set_major_formatter(NullFormatter())
ax.yaxis.set_major_formatter(NullFormatter())
ax.axis('off')

fig.savefig(PathToARIANNAanalysis + '/plots/station51Layout.png')
fig.savefig(PathToARIANNAanalysis + '/plots/station51Layout.pdf')

plt.show()
