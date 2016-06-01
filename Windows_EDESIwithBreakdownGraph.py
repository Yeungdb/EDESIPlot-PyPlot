#!/usr/bin/python
from tkFileDialog import askopenfilename, askdirectory
import tkSimpleDialog
from Tkinter import *

MAINDIRECTORY = ""

class MyDialog(tkSimpleDialog.Dialog):
    def body(self, master):
        #Label(master, text="Line Width for Mass Spectrum (pixel: 1-10):").grid(row=0)
        Label(master, text="Name of Output File").grid(row=0)
        Label(master, text="Minimum Intensity Counts for Contour Filter: ").grid(row=1)
        Label(master, text="Threshold Intensity for Summed Mass Spectrum: ").grid(row=2)
        self.e1 = Entry(master)
        self.e2 = Entry(master)
        self.e3 = Entry(master)
        
        self.e1.grid(row=0, column=1)
        self.e2.grid(row=1, column=1)
        self.e3.grid(row=2, column=1)
        return self.e1
    
    def apply(self):
        outfileName = str(self.e1.get())
        minFilter = int(self.e2.get())
        threshold = int(self.e3.get())
        self.result = outfileName, minFilter, threshold
        
def getFile():
    root = Tk()
    root.update()
    filename = askopenfilename()
    root.destroy()
    return filename
    
def getDirectory():
    root = Tk()
    root.update()
    if(MAINDIRECTORY == ""):
        directory = askdirectory()
    else:
        directory = MAINDIRECTORY
    root.destroy()
    return directory

def getParams():
    root = Tk()
    result = MyDialog(root).result
    root.destroy()
    return result
    
from optparse import OptionParser

parser = OptionParser()
parser.add_option('-f', "--filename", help="Filename of file for processing", action="store")
parser.add_option('-w', "--width", help="Width of line", action="store")
parser.add_option('-m', "--minFilter", help="Minimum Value for Contour Plot Filter", action="store")
parser.add_option('-t', "--threshold", help="Threshold Intensity for Picking MZ in Breakdown Graph", action="store")
parser.add_option('-o', "--outputFile", help="Name of Output File (enter without .html)", action="store")
options, args = parser.parse_args()

filename = getFile()
paramGetter = getParams()

def anon(x):
   print x

import collections
#file = open("C:\Users\Darien\Desktop\UVic\Karlee\MSDataTxt\KB-0173.txt", 'r')
#filename = "/home/darien/KB-0173.txt"

file = open(filename, 'r')
lines = file.readlines()
#print lines[67]

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
            #contourArr[(collenergy, float(j))] = SpecDict[j]
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
#print "Specta Arr"
tArr = collections.OrderedDict(sorted(arr.items()))
#SUMMING Spectra information

c = collections.Counter()
for i in arr:
    c.update(tArr[i]) #Summing values of the same m/z

contourPlotArr = collections.OrderedDict(sorted(contourArr.items()))
#print contourPlotArr
#x, y = zip(*contourPlotArr.keys())
#contourx = contourPlotArr.keys()
#print "X"
#print contourx


#################################
#Sum the data for Summed Full Scan

summed = collections.OrderedDict(sorted(c.items(), key = lambda x: float(x[0])))
#print "summ"


#################################3
#MAKE Breakdown Graphs


threshold = paramGetter[2]
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
#print maxVal
#print minVal
for i in range(minVal, maxVal+1):
   binArr[i] = 0


ZMatrix = []
contourz = contourPlotArr.values()
for i in contourz:
   for j in i:
      binArr[int(float(j))] += i[j]
      #print j
   ZMatrix.append(binArr.values())
   for mz in breakDownMZ:
       listofAllTrace[mz].append(binArr[mz])
   binArr = dict.fromkeys(binArr, 0)
   #print len(i)
#print ZMatrix

#print "breakdown"
#print breakDownMZ
#print binArr

#print "Y"
#print contoury

print "LIST OF ALL TRACE"
print listofAllTrace

#print "Z"
#print contourz
#print summed
#################################
#AT THIS POINT BASICALLY, c has the summed m/z and arr has the time, collision energy and associated spectrum

newArr = collections.OrderedDict(sorted(tArr.items()))



##################################
#PLOTLY PLOTS
import plotly
import plotly.graph_objs as go
import plotly.tools as tls

#SUMMED SPECTRA
#print "KEYS"
#print summed.keys() 
#print "VALUES"
#print summed.values() 

lineWidth = 1

trace0 = go.Scatter(
   x = summed.keys(),
   y = summed.values(),
   xaxis = 'x1',
   yaxis = 'y2',
   mode = 'lines+text',
   name = "Summed Mass Spectrum " + filename,
   line = dict(
      color = ('rgb(205, 12, 24)'),
      width = lineWidth
   )
)

#print "X"
#print binArr.keys()
#print "Y"
#print contourPlotArr.keys()


minFilter = paramGetter[1]
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
          width = lineWidth
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

#------------------------------------------------------------------------------------------

#data = [trace0, trace1]
#layout = go.Layout(
#      yaxis=dict(
#         domain=[0,0.5]
#      ),
#      yaxis2=dict(
#         domain=[05,1]
#      )
#   )

#fig = go.Figure(data=data, layout=layout)

plotly.offline.plot(fig, validate=False, filename='{directory}/{filename}.html'.format(directory=getDirectory(), filename=paramGetter[0]))
