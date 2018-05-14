#!/USR/bin/python2.7

from optparse import OptionParser
import collections
import matplotlib
import matplotlib.gridspec as gridspec
import matplotlib.pyplot as plt
from scipy import *
from math import *
import numpy as numpy

parser = OptionParser()
parser.add_option('-f', "--filename", help="Filename of file for processing", action="store")
parser.add_option('-m', "--minFilter", help="Minimum Value for Contour Plot Filter", action="store")
parser.add_option('-t', "--threshold", help="Threshold Intensity for Picking MZ in Breakdown Graph", action="store")
parser.add_option('-o', "--outputFile", help="Name of Output File without extension", action="store")
parser.add_option('-b', "--breakDown", help="Enable/Disable Breakdown Plots", action="store_true", default=False)
parser.add_option('-z', "--zoom", help="Enable/Disable Zoom", action="store_true", default=False)

parser.add_option('-c', "--collision_start", help="Start index for collision", action="store")
parser.add_option('-k', "--collision_end", help="End index for collision", action="store")

options, args = parser.parse_args()

filename=options.filename

startCollInd = options.collision_start
endCollInd = options.collision_end

#http://scipy-cookbook.readthedocs.io/items/SignalSmooth.html
def smooth(x, window_len=5, window='hanning'):
   x = numpy.array(x)
   if x.ndim != 1:
      raise ValueError, "smooth only accepts 1 dimension arrays."

   if x.size < window_len:
      raise ValueError, "Input vector needs to be bigger than window size."


   if window_len<3:
      return x


   if not window in ['flat', 'hanning', 'hamming', 'bartlett', 'blackman']:
      raise ValueError, "Window is on of 'flat', 'hanning', 'hamming', 'bartlett', 'blackman'"


   s=numpy.r_[x[window_len-1:0:-1],x,x[-2:-window_len-1:-1]]
   #print(len(s))
   if window == 'flat': #moving average
      w=numpy.ones(window_len,'d')
   else:
      w=eval('numpy.'+window+'(window_len)')

   y=numpy.convolve(w/w.sum(),s,mode='valid')
   return list(y)
                                                          






def anon(x):
   #Used for testing lambda functions
   print x

def GraphCleanUp(ax):
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.get_xaxis().tick_bottom()
    ax.get_yaxis().tick_left()

def RoundToTen(val):
    return int(ceil(val/10.0))*10

def RoundToTenFloor(val):
    return int(floor(val/10.0))*10

def RoundToHundred(val):
    return int(ceil(val/100.0))*100

def round1sig(yTick):
    for i in range(0, len(yTick)):
        #stackoverflow: 3410976
        if (yTick[i] == 0):
            continue
        yTick[i] = round(yTick[i], -int(floor(log10(abs(yTick[i])))))
    return yTick

def Valueround1sig(yTick):
    return round(yTick, -int(floor(log10(abs(yTick)))))

file = open(filename, 'r')
#lines = file.readlines()

def normalizeMSSpectrum(mzList, maxValue):
    relInt = []
    for i in mzList:
       relInt.append((float(i)/float(maxValue))*100)
       #print ((float(i)/float(maxValue))*100)
    return relInt

#Loop Var initialization
arr = {} #Time to ColToSpec
ColToSpec = {} #Collision Energy to Spectrum
starttime = 0.0
collenergy = 0
mz = []
intens = []

mzActive = 0 #Switches to see if the mz have been collected
intensActive = 0
contourArr = {}
minVal = 100000000000
maxVal = 0

#To make threshold % based
minInt = 100000000000
maxInt = 0

writeCollEnergy = 0
CECounter =0
startCETrack = startCollInd
endCETrack = endCollInd 
switch_CETrack=0


mzToggle = 0
intensToggle = 0
collectToggle = 0

#Parse through text file to find all the relevant parameters for the experiment
with open(filename, 'r') as MSFileObject:
   for line in MSFileObject:
      SpecDict = {}
      if "scan start time" in line:
         starttime = float(line.split(', ')[1])
      if "collision energy" in line:
         collenergy = int(line.split(', ')[1])
      if mzToggle == 1:
         mz = line.split("] ")[1].split()
         mzToggle = 0
      if "m/z array" in line:
         #mz = lines[i+1].split("] ")[1].split()
         mzActive = 1
         mzToggle = 1
      if intensToggle == 1:
         intens = line.split("] ")[1].split()
         intensToggle = 0
      if "intensity array" in line:
         #intens = lines[i+1].split("] ")[1].split()
         intensActive = 1
         intensToggle = 1
      if ((mzActive & intensActive) != 0):
         SpecDict = collections.OrderedDict(zip(mz, intens))
         for j in SpecDict:
            SpecDict[j] = int(float(SpecDict[j])) #int cast intensity
            if(SpecDict[j] > maxInt):
               maxInt = SpecDict[j]
            if(SpecDict[j] < minInt):
               minInt = SpecDict[j]
            if(int(float(j)) > maxVal):
               maxVal = int(float(j))
            elif(int(float(j)) < minVal):
               minVal = int(float(j))
               print minVal

         contourArr[collenergy] = SpecDict
         arr[(starttime,collenergy)] = SpecDict
         mzActive = 0 #Switches to see if the mz have been collected
         intensActive = 0
         mz = []
         intens = []



print "SORT"
#Dictionary with key of (time, collision energy): ((mz: intensity), (mz: intensity)....)
tArr = collections.OrderedDict(sorted(arr.items()))


#################################
#SUMMING Spectra information

print "SUM"
c = collections.Counter()
for i in arr:
    c.update(tArr[i]) #Summing values of the same m/z

summed = collections.OrderedDict(sorted(c.items(), key = lambda x: float(x[0])))

#################################
#Constructing Z axis matrix for Contour Plot
print "Z ORDER"
contourPlotArr = collections.OrderedDict(sorted(contourArr.items()))

#################################3
#MAKE Breakdown Graphs
#ie. Reconstructed SIM traces of the Reactants, Intermediates and Products

print "BREAKDOWN"
threshold = options.threshold
threshold = float(maxInt-minInt)*(float(threshold)/100.0)
threshold += minInt
print threshold
breakDownMZ = []
localMaxima = 0
for i in summed:
    if(summed[i] > int(threshold)) and (int(float(i)) not in breakDownMZ):
        if(summed[i] > localMaxima):
            localMaxima = summed[i]
        else:
            #Remove local maximas that are located on the side of large local maximas ie 123.3333: 12, 123.4545: 13, 12,5444: 12, 12.6333: 126, 12.7333: 12
            #Essentially, try to remove peaks in peaks
            breakDownMZ.append(int(float(i)))
    else:
        localMaxima = 0

listofAllTrace = {}
for mz in breakDownMZ:
    listofAllTrace[mz] = []

#################################
#BINNING DATA TO UNIT MASS RESOLUTION
print "BINNING"
binArr = {}
for i in range(minVal, maxVal+1):
   binArr[i] = 0

TraceMaxInt = 0
ZMatrix = []
contourz = contourPlotArr.values()
for i in contourz:
   for j in i:
      binArr[int(float(j))] += i[j]
   ZMatrix.append(binArr.values())
   for mz in breakDownMZ:
       listofAllTrace[mz].append(binArr[mz])
       #print str(max(binArr)) + " " + str(TraceMaxInt)
       if(max(binArr) > TraceMaxInt):
           TraceMaxInt = max(binArr) 
   binArr = dict.fromkeys(binArr, 0)

#################################
#AT THIS POINT BASICALLY, c has the summed m/z and arr has the time, collision energy and associated spectrum

print "PLOT"
newArr = collections.OrderedDict(sorted(tArr.items()))

showBreakdown = options.breakDown
line_colours = ('BlueViolet', 'Crimson', 'ForestGreen', 'Indigo', 'Tomato', 'Maroon')

minFilter = int(options.minFilter)
print array(ZMatrix)
levels = arange(minFilter, max(array(ZMatrix).max(axis=1)), 6)

subplotVal = 1
if(showBreakdown == True):
    subplotVal = 2

plt.figure(figsize=(8,9))
matplotlib.rcParams.update({'font.size':16})
matplotlib.rcParams['xtick.direction'] = 'out'
matplotlib.rcParams['ytick.direction'] = 'out'

gs1 = gridspec.GridSpec(subplotVal,2,width_ratios=[3,1], height_ratios=[1,3])
gs1.update(wspace=0.025, hspace=0.025)

print "CONTOUR"
ax2 = plt.subplot(gs1[subplotVal])
ContourPlot = ax2.contour(binArr.keys(), contourPlotArr.keys(), ZMatrix, levels, colors = 'k')
ContouryTick = arange(0, Valueround1sig(max(contourPlotArr.keys()))+10, 10)

print "SETUP TICKS"

step = 10 ##SET TOGGLABLE

#Labels the smallest value of the axis and then continue labelling with the arange increment ie 50, 750, 1500, 2250 rather than 50, 800, 1550, 2300 (tl;dr 50 does not affect the arange values)
ax2.set_yticks(arange(0, max(contourPlotArr.keys()), step))
ContourxTick = arange(0, Valueround1sig(max(binArr.keys()))+750, 750)
ax2.set_ylabel('Collision Energy (V)')
ax2.set_xlabel('m/z', style='italic')
GraphCleanUp(ax2)

print "SUMMED SPECTRA"
ax1 = plt.subplot(gs1[0], sharex=ax2)
SummedFullScan = ax1.plot(summed.keys(), normalizeMSSpectrum(summed.values(), max(summed.values())))
#plt.title('Energy Dependent-ESI MS/MS Plot', y=1.08)
plt.setp(ax1.get_xticklabels(), visible=False)
ax1.set_ylabel('%')
#stackoverflow: 11244514
ax1.set_yticks([0, 50, 100]) #Relative Intensity`
GraphCleanUp(ax1)
ax1.spines['bottom'].set_visible(False)
ax1.xaxis.set_ticks_position('none')


xMinVal = 0
xMaxVal = int(maxVal/500) * 500
ax2.set_xlim(xmin=xMinVal, xmax=xMaxVal)
ax1.set_xlim(xmin=xMinVal, xmax=xMaxVal)

print "ZOOM"
zoom = options.zoom 
if (zoom):
    StartRange = min(listofAllTrace.keys())
    EndRange = max(listofAllTrace.keys())
    print EndRange
    Range = (EndRange-StartRange)/10
    xMinVal = StartRange - Range
    xMaxVal = EndRange + Range
    print RoundToHundred(EndRange)
    xtickRange = arange(RoundToHundred(StartRange)-100, RoundToHundred(EndRange)+100, 100)
    print xtickRange
    ax2.set_xticks(xtickRange)

    ax2.set_xlim(xmin=xMinVal, xmax=xMaxVal)
    ax1.set_xlim(xmin=xMinVal, xmax=xMaxVal)

print "BREAKDOWN GRAPH"
showBreakdown = options.breakDown
windowSize = 11 #Toggle
if(showBreakdown == True):
    ax3 = plt.subplot(gs1[3], sharey=ax2)
    maxOfAllTraces = 0
    for i in listofAllTrace:
        testMaxVal = max(listofAllTrace[i])
        if (testMaxVal > maxOfAllTraces):
            maxOfAllTraces = testMaxVal 
    for mz in breakDownMZ:
        
        ax3.plot(smooth(normalizeMSSpectrum(listofAllTrace[mz], maxOfAllTraces), windowSize), smooth(contourPlotArr.keys(), windowSize)) 
    plt.setp(ax3.get_yticklabels(), visible=False)
    ax3.set_xlabel('%')
    ax3.set_xlim(xmin=0, xmax=100)
    ax3.set_xticks([100])
    GraphCleanUp(ax3)
print "SAVEFIG"
plt.savefig('{OUTFILE}.png'.format(OUTFILE = options.outputFile), bbox_inches='tight')
