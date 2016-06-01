#!/usr/bin/python2.7

from optparse import OptionParser
import collections
import plotly
import plotly.graph_objs as go
import plotly.tools as tls

parser = OptionParser()
parser.add_option('-f', "--filename", help="Filename of file for processing", action="store")
parser.add_option('-w', "--width", help="Width of line", action="store")
parser.add_option('-m', "--minFilter", help="Minimum Value for Contour Plot Filter", action="store")
parser.add_option('-t', "--threshold", help="Threshold Intensity for Picking MZ in Breakdown Graph", action="store")
parser.add_option('-o', "--outputFile", help="Name of Output File (enter without .html)", action="store")
options, args = parser.parse_args()

filename=options.filename

def anon(x):
   print x


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
            print str(i) + " " + str(summed[i]) 
            breakDownMZ.append(int(float(i)))
    else:
        localMaxima = 0
print breakDownMZ

listofAllTrace = {}
for mz in breakDownMZ:
    listofAllTrace[mz] = []




#################################
#BINNING DATA TO UNIT MASS RESOLUTION

binArr = {}
for i in range(minVal, maxVal+1):
   binArr[i] = 0


ZMatrix = []
contourz = contourPlotArr.values()
for i in contourz:
   for j in i:
      binArr[int(float(j))] += i[j]
   ZMatrix.append(binArr.values())
   for mz in breakDownMZ:
       listofAllTrace[mz].append(binArr[mz])
   binArr = dict.fromkeys(binArr, 0)

#################################
#AT THIS POINT BASICALLY, c has the summed m/z and arr has the time, collision energy and associated spectrum

newArr = collections.OrderedDict(sorted(tArr.items()))



##################################
#PLOTLY PLOTS

#Summed Scan Full Spectrum
trace0 = go.Scatter(
   x = summed.keys(),
   y = summed.values(),
   xaxis = 'x1',
   yaxis = 'y2',
   mode = 'lines+text',
   name = "Summed Mass Spectrum " + filename,
   line = dict(
      color = ('rgb(205, 12, 24)'),
      width = int(options.width)
   )
)

minFilter = options.minFilter

#Contour Plot for EDESI (Collision Energy vs m/z)
trace1 = go.Contour(
   z = ZMatrix, 
   x = binArr.keys(),
   y = contourPlotArr.keys(),
   yaxis='y1',
   name = "Energy dependent Electrospray Ionization",
   autocontour=False,
   showscale=True,
   contours=dict(
      coloring='lines',
      start = minFilter,
      end=maxVal,
      size = 200
   )
)

data = [trace0, trace1]
for mz in breakDownMZ:
    trace = go.Scatter(
       x = listofAllTrace[mz],
       name = "" + str(mz),
       xaxis='x2',
       text="" + str(mz) + "m/z",
       textposition = 'top',
       line = dict(
          width = int(options.width)
       )
    )
    data.append(trace)

layout = go.Layout(
        showlegend=False,
        xaxis=dict(
            title='m/z',
            domain=[0, 0.78],
            showgrid=False,
            zeroline=False
            ),
        yaxis=dict(
            title='Collision Energy(V)',
            domain=[0, 0.78],
            showgrid=False,
            zeroline=False
            ),
        margin=dict(
            t=50,
            r=100
            ),
        xaxis2=dict(
            title='Intensity',
            domain=[0.78,1],
            showgrid=False,
            zeroline=False
        ),
        yaxis2=dict(
            title='Intensity',
            domain=[0.78,1],
            showgrid=False,
            zeroline=False
        )
    )


fig = go.Figure(data=data, layout=layout)

plotly.offline.plot(fig, validate=False, filename='{OUTFILE}.html'.format(OUTFILE=options.outputFile))
