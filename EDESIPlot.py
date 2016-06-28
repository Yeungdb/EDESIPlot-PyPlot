#!/usr/bin/python2.7

from optparse import OptionParser
import collections
import plotly
import plotly.graph_objs as go
import plotly.tools as tls
import plotly.plotly as py

parser = OptionParser()
parser.add_option('-f', "--filename", help="Filename of file for processing", action="store")
parser.add_option('-w', "--width", help="Width of line", action="store")
parser.add_option('-m', "--minFilter", help="Minimum Value for Contour Plot Filter", action="store")
parser.add_option('-s', "--fontSize", help="Font Size for Plot", action="store") 
parser.add_option('-t', "--threshold", help="Threshold Intensity for Picking MZ in Breakdown Graph", action="store")
parser.add_option('-o', "--outputFile", help="Name of Output File (enter without .html)", action="store")
parser.add_option('-k', '--tickDist', help="Distance between each tick", action="store")
parser.add_option('-r', '--maxRange', help="Max value in the X Axis", action="store")
parser.add_option('-b', "--breakDown", help="Enable/Disable Breakdown Plots", action="store_true", default=False)
parser.add_option('-p', "--publication", help="Styles the plot to publication standards", action="store_true", default=False)
options, args = parser.parse_args()

filename=options.filename

def anon(x):
   print x


def GraphCleanUp(ax):
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.get_xaxis().tick_bottom()
    ax.get_yaxis().tick_left()

file = open(filename, 'r')
lines = file.readlines()

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
minVal = 1000
maxVal = 0
maxIntensity = 0

for i in range(67, len(lines)):
    SpecDict = {}
    if "scan start time" in lines[i]:
        starttime = float(lines[i].split(', ')[1])
    if "collision energy" in lines[i]:
        collenergy = int(lines[i].split(', ')[1])
    if "m/z array" in lines[i]:
        mz = lines[i+1].split("] ")[1].split()
        mzActive = 1
    if "intensity array" in lines[i]:
        intens = lines[i+1].split("] ")[1].split()
        intensActive = 1
    if ((mzActive & intensActive) != 0):
        SpecDict = collections.OrderedDict(zip(mz, intens))
        for j in SpecDict:
            SpecDict[j] = int(SpecDict[j]) #int cast intensity
            if(SpecDict[j] > maxIntensity):
               maxIntensity = SpecDict[j]
            j = float(j) #float cast m/z
            if(j > maxVal):
               maxVal = int(j)
            elif(j < minVal):
               minVal = int(j)
        contourArr[collenergy] = SpecDict
        arr[(starttime,collenergy)] = SpecDict
        mzActive = 0 #Switches to see if the mz have been collected
        intensActive = 0
        mz = []
        intens = []



#Dictionary with key of (time, collision energy): ((mz: intensity), (mz: intensity)....)
tArr = collections.OrderedDict(sorted(arr.items()))



#################################
#SUMMING Spectra information

c = collections.Counter()
for i in arr:
    c.update(tArr[i]) #Summing values of the same m/z

summed = collections.OrderedDict(sorted(c.items(), key = lambda x: float(x[0])))

#################################
#Constructing Z axis matrix for Contour Plot
contourPlotArr = collections.OrderedDict(sorted(contourArr.items()))



#################################3
#MAKE Breakdown Graphs


threshold = options.threshold
breakDownMZ = []
localMaxima = 0
for i in summed:
    if(summed[i] > int(threshold)) and (int(float(i)) not in breakDownMZ):
        if(summed[i] > localMaxima):
            localMaxima = summed[i]
        else:
            #Remove local maximas that are located on the side of large local maximas ie 123.3333: 12, 123.4545: 13, 12,5444: 12, 12.6333: 126, 12.7333: 12
            #print str(i) + " " + str(summed[i]) 
            breakDownMZ.append(int(float(i)))
    else:
        localMaxima = 0
#print breakDownMZ

listofAllTrace = {}
for mz in breakDownMZ:
    listofAllTrace[mz] = []




#################################
#BINNING DATA TO UNIT MASS RESOLUTION

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

newArr = collections.OrderedDict(sorted(tArr.items()))


import matplotlib
import matplotlib.gridspec as gridspec
import matplotlib.pyplot as plt
from scipy import *
from math import log10, floor

def round1sig(yTick):
    for i in range(0, len(yTick)):
        #stackoverflow: 3410976
        if (yTick[i] == 0):
            continue
        yTick[i] = round(yTick[i], -int(floor(log10(abs(yTick[i])))))
    return yTick

def Valueround1sig(yTick):
    return round(yTick, -int(floor(log10(abs(yTick)))))

showBreakdown = options.breakDown
line_colours = ('BlueViolet', 'Crimson', 'ForestGreen', 'Indigo', 'Tomato', 'Maroon')

minFilter = int(options.minFilter)
#print minFilter
levels = arange(minFilter, max(array(ZMatrix).max(axis=1)), 6)
#line_widths = (1, 1.5, 2, 2.5, 3, 3.5)

subplotVal = 1
if(showBreakdown == True):
    subplotVal = 2

plt.figure(figsize=(8,9))
matplotlib.rcParams.update({'font.size':16})
matplotlib.rcParams['xtick.direction'] = 'out'
matplotlib.rcParams['ytick.direction'] = 'out'

gs1 = gridspec.GridSpec(subplotVal,2,width_ratios=[3,1], height_ratios=[1,3])
gs1.update(wspace=0.025, hspace=0.025)


ax2 = plt.subplot(gs1[subplotVal])
ContourPlot = ax2.contour(binArr.keys(), contourPlotArr.keys(), ZMatrix, levels, colors = 'k')
ContouryTick = arange(0, Valueround1sig(max(contourPlotArr.keys()))+10, 10)

#Labels the smallest value of the axis and then continue labelling with the arange increment ie 50, 750, 1500, 2250 rather than 50, 800, 1550, 2300 (tl;dr 50 does not affect the arange values)
ax2.set_yticks(concatenate(([Valueround1sig(min(contourPlotArr.keys()))],ContouryTick[:-1]), axis=0))
#print ContouryTick[:-1]
ContourxTick = arange(0, Valueround1sig(max(binArr.keys()))+750, 750)
#print ContourxTick[:-1]
#ax2.set_xticks(concatenate(([Valueround1sig(min(binArr.keys()))],ContourxTick[:-1]), axis=0))
ax2.set_ylabel('Collision Energy (V)')
ax2.set_xlabel('m/z', style='italic')
GraphCleanUp(ax2)


def normalizeMSSpectrum(mzList, maxValue):
    relInt = []
    for i in mzList:
       relInt.append((float(i)/float(maxValue))*100)
       #print ((float(i)/float(maxValue))*100)
    return relInt

ax1 = plt.subplot(gs1[0], sharex=ax2)
#SummedFullScan = ax1.plot(summed.keys(), summed.values()) #Absolute Intensity
SummedFullScan = ax1.plot(summed.keys(), normalizeMSSpectrum(summed.values(), max(summed.values())))
plt.title('Energy Dependent-ESI MS/MS Plot', y=1.08)
plt.setp(ax1.get_xticklabels(), visible=False)
#ax1.set_ylabel('Intensity') #Absolute Intensity
ax1.set_ylabel('Intensity (%)')
#stackoverflow: 11244514
#ax1yTick = arange(0, Valueround1sig(max(summed.values()))+10000, 10000) #Absolute Intensity
ax1.set_yticks([0, 50, 100]) #Relative Intensity`
GraphCleanUp(ax1)
ax1.spines['bottom'].set_visible(False)
ax1.xaxis.set_ticks_position('none')

def RoundToTen(val):
    import math
    return int(math.ceil(val/10.0))*10

def RoundToTenFloor(val):
    import math
    return int(math.floor(val/10.0))*10

def RoundToHundred(val):
    import math
    return int(math.ceil(val/100.0))*100

zoom = True
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
    #for i in range(0, len(xtickRange)):
    #   xtickRange[i] = RoundToTenFloor(xtickRange[i])
    #print xtickRange
    ax2.set_xticks(xtickRange)

    #xMinVal = RoundToTenFloor(xMinVal)
    #xMaxVal = RoundToTenFloor(xMaxVal)
    ax2.set_xlim(xmin=xMinVal, xmax=xMaxVal)
    ax1.set_xlim(xmin=xMinVal, xmax=xMaxVal)

showBreakdown = options.breakDown
if(showBreakdown == True):
    ax3 = plt.subplot(gs1[3], sharey=ax2)
    maxOfAllTraces = 0
    for i in listofAllTrace:
        testMaxVal = max(listofAllTrace[i])
        if (testMaxVal > maxOfAllTraces):
            maxOfAllTraces = testMaxVal 
    for mz in breakDownMZ:
        #ax3.plot(listofAllTrace[mz], contourPlotArr.keys()) #Absolute Intensity
        ax3.plot(normalizeMSSpectrum(listofAllTrace[mz], maxOfAllTraces), contourPlotArr.keys()) 
    plt.setp(ax3.get_yticklabels(), visible=False)
    #BreakDownInt = arange(0, TraceMaxInt, 1500)
    ax3.set_xlabel('Intensity (%)')
    #print TraceMaxInt
    #ax3.set_xticks(BreakDownInt)
    ax3.set_xticks([' ', ' ',100])
    GraphCleanUp(ax3)

plt.savefig('../{OUTFILE}.png'.format(OUTFILE = options.outputFile), bbox_inches='tight')
