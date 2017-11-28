# EDESI Plot
## Darien Yeung ([@Yeungdb](https://twitter.com/Yeungdb))


## Installation Guide

There a few things that needs to be installed before it can run properly.

1) [Canopy](https://store.enthought.com/downloads/#default) (This is Free)

2) [Proteowizard](http://proteowizard.sourceforge.net/downloads.shtml) (This is Free)


## Operation (Running) Steps

1) Open Canopy and then open EDESIwithBreakdownGraph.py.

2) Then open MSConvert inside Proteowizard:

   (Windows 7) Go to the start menu > All Programs   >   Proteowizard    >    MSConvert

   (Windows 8 and above) Go to the start menu  > Type in MSConvert

3) Find your .raw file and using MSConvert, convert the .raw into .txt and remember where you saved the .txt file

4) Go back to Canopy and press the Green Play Button.

5) First it will ask you for the .txt file, find the txt file and open it with the popup window.

6) It will ask you for:

   Name of the Output File: (Type in the name you want your file to be)

   Minimum Filter Value of the Contour Plot:  (This is the smallest threshold value that the contour plot will show)

      (For example, if you typed in 20, peak in a spectra with a peak intensity higher than 20 counts will be coloured in the contour plot)

   Threshold Intensity for Summed Mass Spectra: 

      (This is one is for the break down graph. The value you enter say 1150, any peak in the summed spectra above 1150 counts will be selected for plotting in the breakdown graph)

7) Once you have filled out the values, press enter and let the program run (~10 seconds)

8) Once its finished running, it will ask you for a folder to save in. Select which folder you want your EDESI plot to be saved in.

9) It will automatically open the EDESI plot that it generated from the txt file.

