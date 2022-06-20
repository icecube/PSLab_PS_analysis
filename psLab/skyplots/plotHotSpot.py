#!/usr/bin/env python
import numpy as np
import pylab 
from numpy.random import uniform, seed
from pylab import *

from matplotlib.projections.geo import HammerAxes

import  mpl_toolkits.axes_grid.angle_helper as angle_helper
from mpl_toolkits.axes_grid.grid_finder import  GridFinder
from mpl_toolkits.axes_grid.axislines import  AxisArtistHelper, GridHelperBase, AxisArtist
from mpl_toolkits.axes_grid.parasite_axes import SubplotHost, ParasiteAxesAuxTrans
from matplotlib.transforms import Affine2D

from  mpl_toolkits.axes_grid.grid_helper_curvelinear import GridHelperCurveLinear

import ROOT

cdict = {'blue': ((0.0, 0.0, 1.0),
                  (0.05, 1.0, 1.0),
                  (0.4, 1.0, 1.0),
                  (0.6, 1.0, 1.0),
                  (0.7, 0.2, 0.2),
                  (1.0, 0.0, 0.0)),
         'green': ((0.0, 0.0, 1.0),
                   (0.05, 1.0, 1.0),
                   (0.5, 0.0416, 0.0416),
                   (0.6, 0.0, 0.0),
                   (0.8, 0.5, 0.5),
                   (1.0, 1.0, 1.0)),
         'red': ((0.0, 0.0, 1.0),
                 (0.05, 1.0, 1.0),
                 (0.5, 0.0416, 0.0416),
                 (0.6, 0.0416, 0.0416),
                 (0.7, 1.0, 1.0),
                 (1.0, 1.0, 1.0))}

my_cmap = matplotlib.colors.LinearSegmentedColormap('root_colormap',cdict,256)


f = ROOT.TFile("../macro_llh/ic79/ExtendedSky/IC79_IC59_IC40_AllSky_E2_Sirin_Sigma_5.root")

if (f.IsOpen()):
    print "Openning root file "
    roothistfine = ROOT.TH2D(f.Get("hExtendedSkyCoarse"))
else:
    print "Unable to open root file ", self.file
    
try:
    if not roothistfine.InheritsFrom("TH2"):
        raise TypeError("%s does not inherit from TH2" % hist)
except:
    raise TypeError("%s is not a ROOT object" % hist)


fevent = open("resources/IC40_events_noPoles_ASCII.txt","r")

evxcoord = []
evycoord = []

for line in fevent:
    evxcoord.append(float(line.split()[0]))
    evycoord.append(float(line.split()[1]))


fevent59 = open("resources/IC59_events_noPoles_ASCII.txt","r")

evxcoord59 = []
evycoord59 = []

for line in fevent59:
    evxcoord59.append(float(line.split()[0]))
    evycoord59.append(float(line.split()[1]))
    

fevent79 = open("resources/IC79_events_noPoles_ASCII.txt","r")

evxcoord79 = []
evycoord79 = []

for line in fevent79:
    evxcoord79.append(float(line.split()[0]))
    evycoord79.append(float(line.split()[1]))
    


array = np.array

nbinsx = nx = roothistfine.GetNbinsX()
nbinsy = ny = roothistfine.GetNbinsY()

maxbin = roothistfine.GetMaximumBin()
maxval = roothistfine.GetBinContent(maxbin)

i, j = ROOT.Long(0), ROOT.Long(0)

roothistfine.GetBinWithContent2(maxval, i, j, 1, 3600, 1, 1700, 0.01)

#ra_max = roothistfine.GetXaxis().GetBinCenter(i)
#dec_max = roothistfine.GetYaxis().GetBinCenter(j)

#ra_max = 219.25
#dec_max = -38.75

ra_max = 254
dec_max = -58.24


#Check the y limits, we need from -90 to 90

ymin = roothistfine.GetYaxis().GetBinLowEdge(1)
ymax = roothistfine.GetYaxis().GetBinLowEdge(ny + 1)
binwidth = roothistfine.GetXaxis().GetBinWidth(1)

underbins = int((ymin + 90) / binwidth)
overbins = int((90 - ymax) / binwidth)

print "Underbins in y-axis: ", underbins
print "Overbins in y-axis: ", overbins


content = array([array([roothistfine.GetBinContent(i, j)
                        for i in range(1, nx + 1)])
                 for j in range(1 - underbins, ny + 1 + overbins)])

xedges = array([roothistfine.GetXaxis().GetBinLowEdge(i)
                for i in range(1, nx + 2)])


yedges = array([roothistfine.GetYaxis().GetBinLowEdge(i)
                for i in range(1 - underbins, ny + 2 + overbins)])

xedgesrad = xedges * np.pi/180. - np.pi
yedgesrad = yedges * np.pi/180.


x  = array([(xedges[i+1] + xedges[i])/2
                for i in range(nx)])


y   = array([(yedges[i+1] + yedges[i])/2
             for i in range(ny + overbins + underbins)])


xrad = x * np.pi/180. - np.pi
yrad = y * np.pi/180.

fig = figure()

fig.canvas.set_window_title('Sky Map')

'''
tr = HammerAxes.HammerTransform().transform
inv_tr = HammerAxes.InvertedHammerTransform().transform

extreme_finder = angle_helper.ExtremeFinderCycle(20, 20, lon_cycle = None, lat_cycle = None, lon_minmax = None, lat_minmax = None)

grid_locator1 = angle_helper.LocatorDMS(12)

tick_formatter1 = angle_helper.FormatterDMS()

grid_helper = GridHelperCurveLinear((tr, inv_tr))
'''




#ax1 = SubplotHost(fig, 1, 1, 1, grid_helper=grid_helper)
#plt = fig.add_subplot(ax1)

plt = fig.add_subplot(111)

skymap = plt.pcolormesh(x, y, content, antialiased=True, cmap=get_cmap(my_cmap), shading='gouraud', vmin=0, vmax = 1.5)

ic40 = plt.scatter(evxcoord, evycoord, s = 15, marker='x')
ic59 = plt.scatter(evxcoord59, evycoord59, s = 15, marker='o',c='black')
ic79 = plt.scatter(evxcoord79, evycoord79, s = 15, marker='^',c='black')

plt.set_xlim(ra_max - 10, ra_max + 10)
plt.set_ylim(dec_max - 10, dec_max + 10)

plt.legend((ic40,ic59, ic79),('IC40 events', 'IC59 events', 'IC79 events'), loc='upper right', scatterpoints=1)

zbar = colorbar(skymap, orientation="vertical")
zbar.set_label("-log$_{10}$ p", fontsize=15)
plt.set_xlabel("RA",fontsize=15)
plt.set_ylabel("$\delta$",fontsize=15)

plt.annotate('PRELIMINARY', xy=(.9, .1),
              xycoords='axes fraction',
              horizontalalignment='right', verticalalignment='bottom',
              fontsize=25, color='red')

#ax1.grid(True)
grid(True)

savefig("hotspot_north.png")

show()
