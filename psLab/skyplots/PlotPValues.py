from SkyMap import *
import sys, glob, dashi, numpy

dashi.visual()

basedir='/net/user/jfeintzeig/2013/psLab/macro_llh/ic86/ScrambledSkyMaps/pvalues/'
OutDir='/net/user/jfeintzeig/2013/psLab/macro_llh/ic86/FinalPlots/'

UNBLINDING=True
if UNBLINDING:
	map_name=sys.argv[1]

def LoadScrambledPValues(files):
	'''
	files: list of text files containing scrambled pvalues to load
	'''
	PValues=[]
	for file in files:
		f=open(file)
		f_pvalues=numpy.array(f.readlines(),dtype=float)
		PValues=numpy.concatenate((PValues,f_pvalues))
		f.close()
	print "loaded %i p-values from %i files" % (len(PValues),len(files))
	return PValues

def MakePrettyScrambledPlot():
        pylab.xlabel(r'-log$_{10}$(p-value)',fontsize=18)
        pylab.ylabel('Fraction of trials',fontsize=18)
	pylab.grid(True)
	pylab.figtext(0.65,0.37,'ICECUBE\nPRELIMINARY',color='red',fontsize=18,
		fontweight='bold')
	pylab.legend(prop={'size':'16'})

def PlotScrambledPValues(arr,name):
        bins=numpy.linspace(3,9,50)
        hist=dashi.factory.hist1d(arr,bins)
	hist.bincontent/=float(len(arr))
        hist.line(linewidth=2)
	MakePrettyScrambledPlot()
        pylab.savefig(OutDir+name)
        #pylab.show()
	pylab.clf()
        return hist

def PlotUnblindedPValue(scrambled_arr,scrambled_hist,hemisphere,map_name=None):
	# plot scrambles
	scrambled_hist.line(linewidth=2,color='black',label='Scrambled trials')
	y=1.1*max(scrambled_hist.bincontent)
	if hemisphere=='North':
		sign='\geq'
	if hemisphere=='South':
		sign='<'

	# get hottest spot from map
	if map_name:
		s=SkyMap()
		s.AddHistogramFile(map_name)
		if hemisphere=='North':
			HottestSpot=s.maxvalN
		elif hemisphere=='South':
			HottestSpot=s.maxvalS
		else:
			sys.exit("What kinda hemisphere is %s?" % (hemisphere))
	else:
		HottestSpot=4.7
		pylab.text(HottestSpot+0.06,y,'FAKE PVALUE!',fontsize=16)

	# plot
	pylab.axvline(HottestSpot,linestyle='solid',linewidth=3,color='red',
		label='Observed p-value')
	PostTrialsPValue=sum(scrambled_arr>=HottestSpot)/float(len(scrambled_arr))
	pylab.figtext(0.525,0.65,'Post-trials p-value: %5.3f'
		% (PostTrialsPValue),fontsize=17)
	pylab.title(r'%sern Sky ($\delta %s -5^{\circ}$) Hottest Spot'
		% (hemisphere,sign),fontsize=20)

	x=scrambled_hist.bincenters
	MakePrettyScrambledPlot()
	pylab.savefig(OutDir+'Unblinded_%s.png' % (hemisphere))
	pylab.show()
	pylab.clf()

PValues={}
Hists={}

for item in ['North','South']:
	# get list of files
	files=glob.glob('%s*%s.txt' % (basedir,item))
	# open files, get array of p-values
	PValues[item]=LoadScrambledPValues(files)
	# plot scrambled p-values
	Hists[item]=PlotScrambledPValues(PValues[item],'Scrambled_%s' % (item))

	# overplot unblinded p-value
	if UNBLINDING:
		PlotUnblindedPValue(PValues[item],Hists[item],item,map_name)
	else:
		PlotUnblindedPValue(PValues[item],Hists[item],item)

