import numpy as np
import matplotlib.pyplot as plt
from colour import Color
import os

PathToARIANNAanalysis = os.environ['ARIANNAanalysis']

# This nice chunk of code that returns an array of hex color values
# Matplotlib bugs out when the hex values are not length 6, very rarily occurs
# This is why the additional for loop is done to verify the length of the colors
# NOTE: therefor the length of colors may be smaller than numColors
# For numColors = 1000, there are about 5 hex codes that do not have length 6 and therefor the returned array is length 995 or so
def getColorGradient(numColors):
    numColors = 1000
    red = Color("red")
    colorsHex = list(red.range_to(Color("purple"),numColors))
    colors = []
    for i in range(len(colorsHex)):
       colors.append(str(colorsHex[i]))
    index = 0
    for i in range(len(colors)):
       if len(colors[index]) < 6:
          colors.pop(index)
          index-=1
       index+=1
    return colors


def linear(x, slope, yIntercept):
    return slope * x + yIntercept

def quad(x, width):
    return width * x * x

colors = getColorGradient(1000)

yIntercept = 0.0
func = np.vectorize(linear)
xs = np.linspace(0,400,100) # Kelvin

fig1, ax1 = plt.subplots(1, 1, sharex=False, figsize=(5, 5))
for i in range(len(colors)):
    ax1.plot(xs,func(xs,i,yIntercept),color=colors[i],linewidth=0.5)

ax1.set_ylabel(r'X')
ax1.set_xlabel(r'Y')
fig1.tight_layout()
fig1.savefig(PathToARIANNAanalysis + '/plots/colorGradientBeauty.png')
fig1.savefig(PathToARIANNAanalysis + '/plots/colorGradientBeauty.pdf')


func = np.vectorize(quad)
xs = np.linspace(-10,10,100) # Kelvin

fig1, ax1 = plt.subplots(1, 1, sharex=False, figsize=(5, 5))
for i in range(len(colors)):
    ax1.plot(xs,func(xs,i+1),color=colors[i],linewidth=0.5)

ax1.set_ylabel(r'X')
ax1.set_xlabel(r'Y')
fig1.tight_layout()
fig1.savefig(PathToARIANNAanalysis + '/plots/colorGradientBeauty2.png')
fig1.savefig(PathToARIANNAanalysis + '/plots/colorGradientBeauty2.pdf')

plt.show()
