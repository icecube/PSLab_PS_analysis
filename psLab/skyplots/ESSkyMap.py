import numpy as np
import pylab 
from numpy.random import uniform, seed
from pylab import *

import ROOT


""" Color map similar to what is being used in ROOT """
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

class ESSkyMap:
    """ A python class to draw skymaps """
    
    def __init__(self):

        self.cmap = my_cmap
        self.divlong = 6
        self.divlat = 6
        self.orientation = "horizontal"
        self.ztitle = "-log$_{10}$ p" 
        self.gridcolor = "blue"
        self.gridwidth = 0.7
        self.gridstyle = ":"
        self.galfile = False
        self.galcoord = []
        self.lcolor = "black"
        self.lstyle = "-"
        self.lwidth = 1
        self.evfile = False
        self.evxcoord = []
        self.evycoord = []
        self.evmarker = "o"
        self.ymin = -85
        self.ymax = 85
        self.yticks = []
        self.sourcelist =  ""
        self.sourcenames = []
        self.sourcera = []
        self.sourcedec = []
        self.annotate = False
        self.vmax = 6


        """ Values for the maximum significance """
        
        self.max = 0
        self.ra_max = 0
        self.dec_max = 0

        
    def SetZorientation(self, orientation):
        """ Set the orientation of the z bar [\"horizontal\"/\"vertical\"]. """
        
        if orientation in ["horizontal", "vertical"]:
            self.orientation = orientation
        else:
            print "Not defined orientation %s " % orientation
            print "Setting default to horizontal"
            self.orientation = "horizontal"

    
    def SetCMap(self,cmap):
        """ Method to choose the heat color map. Type SetCMap("") to retrieve the list of color maps available. """
        
        if cmap in [m for m in cm.datad]:
            self.cmap = cmap
        else:
            print "No colormap %s found" %cmap
            print "List of color maps are: "
            for m in cm.datad:
                if not m.endswith("_r"):
                    print " = ", m
                    
            print " = root_colormap (default) "


    def SetZtitle(self, ztitle):
        """ Set title of the z bar. """
        self.ztitle = ztitle
        

    def SetGridColor(self, color):
        self.gridcolor = color
        
    def AddGalPlane(self, galfile, lcolor="black", lstyle="-", lwidth=0.7):
        """ Add file to Galactic Plane coordinates. Also define the style of the line using matplotlib standards."""
        
        self.galcoord = []
        self.lcolor = lcolor
        self.lstyle = lstyle
        self.lwidth = lwidth
        
        
        f = open(galfile, "r")
        
        for line in f.readlines():
            self.galcoord.append( [float(line.split()[0]), float(line.split()[1])] )

        self.galcoord.sort()

    def AddSourceList(self, file, Annotate=False):
        self.sourcelist = file
        self.annotate = Annotate
        
        f = open(self.sourcelist)

        for line in f:
            name = line.split()[0]
            ra = (float(line.split()[1]))
            dec = float(line.split()[2])

            self.sourcenames.append(name.replace("_", " "))
            self.sourcera.append(ra)
            self.sourcedec.append(dec)



    def AddEvents(self, evfile):
        """ Add file with raDeg and decDeg of the events. """
        self.evxcoord = []
        self.evycoord = []

        fevent = open(evfile, "r")

        for line in fevent:
            
            self.evxcoord.append(float(line.split()[0]))
            self.evycoord.append(float(line.split()[1]))

        

            
    def AddHistogramFile(self, file, ClockWise=False):
        """ Method to add a ROOT histrogram file. The file must contain a TH2D histogram with name hAllSkyCoarse """
        """ You can decide from 0 to 24h (ClockWise=True) or 24h to 0 (ClokWise=False (default))"""
        
        self.file = file
        
        f = ROOT.TFile(self.file)

        if (f.IsOpen()):
            roothist = ROOT.TH2D(f.Get("hExtendedSkyCoarse"))
        else:
            print "Unable to open root file ", self.file

        try:
            if not roothist.InheritsFrom("TH2"):
                raise TypeError("%s does not inherit from TH2" % hist)
        except:
            raise TypeError("%s is not a ROOT object" % hist)


        array = np.array

        self.nbinsx = nx = roothist.GetNbinsX()
        self.nbinsy = ny = roothist.GetNbinsY()


        maxbin = roothist.GetMaximumBin()
        self.maxval = roothist.GetBinContent(maxbin)
        i, j = ROOT.Long(0), ROOT.Long(0)
        
        roothist.GetBinWithContent2(self.maxval, i, j, 1, 3600, 1, 1700, 0.01)

        self.ra_max = roothist.GetXaxis().GetBinCenter(i)
        self.dec_max = roothist.GetYaxis().GetBinCenter(j)
        
        


        #Check the y limits, we need from -90 to 90

        ymin = roothist.GetYaxis().GetBinLowEdge(1)
        ymax = roothist.GetYaxis().GetBinLowEdge(ny + 1)
        binwidth = roothist.GetXaxis().GetBinWidth(1)

        underbins = int((ymin + 90) / binwidth)
        overbins = int((90 - ymax) / binwidth)

        self.ymin = ymin
        self.ymax = ymax
        
        print "Underbins in y-axis: ", underbins
        print "Overbins in y-axis: ", overbins
        

        self.content = array([array([roothist.GetBinContent(i, j)
                                     for i in range(1, nx + 1)])
                              for j in range(1 - underbins, ny + 1 + overbins)])

        '''
        self.content = array([array([0
                                     for i in range(1, nx + 1)])
                              for j in range(1 - underbins, ny + 1 + overbins)])
        '''
        
    
        self.xedges = array([roothist.GetXaxis().GetBinLowEdge(i)
                             for i in range(1, nx + 2)])
        
        
        self.yedges = array([roothist.GetYaxis().GetBinLowEdge(i)
                             for i in range(1 - underbins, ny + 2 + overbins)])
        
        self.xedgesrad = self.xedges * np.pi/180. - np.pi
        self.yedgesrad = self.yedges * np.pi/180.
        
        
        self.x      = array([(self.xedges[i+1] + self.xedges[i])/2
                             for i in range(nx)])


        self.y      = array([(self.yedges[i+1] + self.yedges[i])/2
                             for i in range(ny + overbins + underbins)])
        
        
        self.xrad = self.x * np.pi/180. - np.pi
        self.yrad = self.y * np.pi/180.


    def set_longitude_grid(self, div_long):
        """ Set number of division in longitud """

        self.divlong = div_long
        
        


    def set_latitude_grid(self, div_lat):
        """ Set number of division in latitude """

        self.divlat = div_lat
        
        
        
    def Plot(self, proj, clockwise=False):

        rc('grid', color=self.gridcolor, linewidth=self.gridwidth, linestyle=self.gridstyle)


        
        #if self.orientation == "horizontal":
        #    w, h = figaspect(.8)
        #else:
        #    w, h = figaspect(.5)
        
        
        #w, h = (10, 7)
        w, h = (15, 10)

        
        fig = figure(figsize=(w,h))
            
        fig.canvas.set_window_title('Sky Map')

        

        if proj not in ["aitoff", "hammer", "lambert", "mollweide", "rectangular"]:

            print "No projection %s defined" % proj
            print "The list of accepted projections is:" % proj
            for p in ["aitoff", "hammer", "lambert", "mollweide", "rectangular"]:
                print p

            print "Setting default as Hammer-Aitoff"
            proj = "hammer"
            

        if(proj == "rectangular"):
            plt = fig.add_subplot(111)

            skymap = plt.pcolormesh(self.x, self.y, self.content, antialiased=True, cmap=get_cmap(self.cmap), shading='gouraud', vmin=0, vmax = self.vmax)
                            
            #skymap = plt.imshow(self.content, interpolation='nearest',
            #                    extent=[self.xedges[0], self.xedges[-1],
            #                            self.yedges[0], self.yedges[-1]],
            #                    aspect='auto', cmap=get_cmap(self.cmap))
            

            yticks = [self.ymin]
            yticks = np.append(yticks, np.arange(-90 + 180/self.divlat, 90, 180/self.divlat))
            yticks = np.append(yticks, [self.ymax])
            
            plt.set_yticks(yticks)

            ylabels = ["%+i$^{\circ}$" % int(l) for l in yticks]

            for n, label in enumerate(ylabels):
                if label == "+0$^{\circ}$":
                    ylabels[n] = "0$^{\circ}$"
                    
                    
            plt.set_yticklabels(ylabels, fontsize=15)

            xticks = [self.xedges[0]]
            xticks = np.append(xticks, np.arange(0, 360, 360/self.divlong))
            xticks = np.append(xticks, self.xedges[-1])
            
            if clockwise == True:
                plt.set_xlim(self.xedges[0], self.xedges[-1])
            else:
                plt.set_xlim(self.xedges[-1], self.xedges[0])
                xticks = xticks[::-1]
            
                
            plt.set_xticks(xticks)

            for label in plt.get_xticklabels():
                label.set_fontsize(15)

                
        else:
            
            plt = fig.add_subplot(111, projection=proj)

            plt.set_longitude_grid_ends(self.ymax)

            plt.set_latitude_grid(180/self.divlat)            
            plt.set_longitude_grid(360/self.divlong)
            
            yticks = [self.ymin * np.pi/180.]

            yticks = np.append(yticks, plt.get_yticks())
            yticks = np.append(yticks, [self.ymax * np.pi/180.])
            
            plt.set_yticks(yticks)



            ylabels = ["%+i$^{\circ}$" % math.ceil(l * 180/np.pi) for l in yticks]

            #Remove the 0 label
            
            for n, label in enumerate(ylabels):
                if label == "+0$^{\circ}$":
                    ylabels[n] = ""
            
            plt.set_yticklabels(ylabels, fontsize=15, family='cursive', weight='ultralight')

            plt.set_xticklabels([],fontsize =0)

            alignment = {'horizontalalignment':'center', 'verticalalignment':'center'}

            xmargin = np.pi/20

            if proj == "lambert":
                xmargin = xmargin + 1.9 * np.pi
                
            if clockwise == True:
                skymap = plt.pcolormesh(self.xrad, self.yrad, self.content, antialiased=True, cmap=get_cmap(self.cmap), shading='gouraud', vmin=0, vmax=self.vmax)
                
                plt.text(-np.pi - xmargin, 0, "0h", size="large", weight='semibold',**alignment)
                plt.text(np.pi + xmargin, 0, "24h", size="large", weight='semibold', **alignment)

                
            else:
                skymap = plt.pcolormesh(self.xrad[::-1], self.yrad, self.content, antialiased=True, cmap=get_cmap(self.cmap), shading='gouraud', vmin=0, vmax=self.vmax)

                
                plt.text(-np.pi - xmargin, 0, "24h",size="large", weight='semibold', **alignment)
                plt.text(np.pi + xmargin, 0, "0h",size="large", weight='semibold', **alignment)

                
                

        if (len(self.evxcoord) is not 0):

            if proj == "rectangular":
                xcoord = self.evxcoord
                ycoord = self.evycoord
            else:
                xcoord = [(x - 180) * np.pi/180. for x in self.evxcoord]
                ycoord =  [y*np.pi/180. for y in self.evycoord]

                if clockwise == False:
                    xcoord = [(360 - x - 180)*np.pi/180. for x in self.evxcoord]


                    

            
            xlims = plt.get_xlim()
            plt.scatter(xcoord, ycoord, s=0.05,marker='x')
            plt.set_xlim(xlims)

        
        if (len(self.galcoord) is not 0):
            
            if proj == "rectangular":
                xcoord = [coord[0] for coord in self.galcoord if coord[0]]
                ycoord = [coord[1] for coord in self.galcoord if coord[0]]

            else:
                xcoord = [(coord[0] - 180)* np.pi/180. for coord in self.galcoord]
                ycoord = [coord[1]*np.pi/180. for coord in self.galcoord]

                if clockwise == False:
                    xcoord = [(360 - coord[0] - 180)* np.pi/180. for coord in self.galcoord]
                    


            # Plotting the galactic plane changes the limit
            
            xlims = plt.get_xlim()
                    
            plt.plot(xcoord, ycoord, linestyle=self.lstyle, c=self.lcolor, linewidth=self.lwidth)
            
            plt.set_xlim(xlims)
                    
            
        if (len(self.sourcenames) is not 0):

            if proj is not "rectangular":

                if clockwise == False:
                    self.sourcera = [(360 - ra - 180)*np.pi/180 for ra in self.sourcera]
                else:
                    self.sourcera = [(ra - 180)*np.pi/180 for ra in self.sourcera]

                self.sourcedec = [d * np.pi/180 for d in self.sourcedec]

            xlims = plt.get_xlim()
                
            plt.scatter(self.sourcera, self.sourcedec, lw=2, facecolor='white', marker='o', s=60)



            if self.annotate == True:
                for i, name in enumerate(self.sourcenames):
                    xyt = (-20, -20)

                    if name.startswith("Crab"):
                        xyt = (2, 20)

                    if name.startswith("4C"):
                        xyt = (-20, 20)

                    if name.startswith("Cyg"):
                        xyt = (-30, 20)

                    if name.startswith("Srg"):
                        xyt = (-20, -20)

                    if name.startswith("IC4"):
                        xyt = (-20, 20)

                    if name.startswith("1ES"):
                        xyt = (20, 20)

                    if name.startswith("Mrk"):
                        xyt = (20, 20)

                    if name.startswith("Cas"):
                        xyt = (20, -10)

                    if name.startswith("MGRO"):
                        xyt = (-20, -30)


                    aa = annotate(name, xy=(self.sourcera[i],self.sourcedec[i]),  xycoords='data',
                                  xytext=xyt, textcoords='offset points',
                                  size=12,
                                  bbox=dict(boxstyle="round", fc="0.7", ec="none", alpha=0.5),
                                  arrowprops=dict(arrowstyle="wedge,tail_width=0.7",
                                                  fc="0.7", ec="none",
                                                  #   relpos=(0.2, 0.8),
                                                  connectionstyle="arc3,rad=-0.9"),
                                  )
    
                plt.set_xlim(xlims)
                

        if len(self.yticks) is not 0:
            plt.set_yticks(self.yticks)
        
        zbar = colorbar(skymap, orientation=self.orientation)

        zbar.set_label(self.ztitle, fontsize=25)
        
        for label in zbar.ax.get_xticklabels():
            label.set_fontsize(15)
            
        for label in zbar.ax.get_yticklabels():
            label.set_fontsize(15)
            
            
        grid(True)

                
        print "Fine Grid Hottest Spot:"
        print "-log10(p) = ",self.maxval
        print "Ra: ",self.ra_max
        print "Dec:",self.dec_max        
        
        
