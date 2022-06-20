from SkyMap import *
import sys, glob, numpy

indir=sys.argv[1]
outname=sys.argv[2]

OutDir='/net/user/jfeintzeig/2013/psLab/macro_llh/ic86/ScrambledSkyMaps/pvalues/'

files=glob.glob(indir+'*.root')

print "Getting hottest spots from %i files" % (len(files))

def GetHottestSpots(file):
	s=SkyMap()
	s.AddHistogramFile(file)
	HS_N=s.maxvalN
	HS_S=s.maxvalS
	return HS_N, HS_S

HotSpots={'North':[],'South':[]}
Hist={}

i=0
for file in files:
	N, S = GetHottestSpots(file)
	HotSpots['North']+=[N]
	HotSpots['South']+=[S]
	i+=1
	if i%10==0:
		print "already read %i/%i maps" % (i, len(files))

for item in HotSpots:
	HotSpots[item]=numpy.array(HotSpots[item])
	out=open(OutDir+'PValues_%s_%s.txt' % (outname,item),'a')
	for val in HotSpots[item]:
		out.write('%5.10f\n' % (val))

	out.close()
