# Marcus Andersson - Percolator Project
# Takes data from IntegrationTest_Percolator_SpeedGraphData.py and output graph images.
# Input: Any number of data files, but at least 1.

pathToOutputData = "@pathToOutputData@"
pathToTestScripts = "@pathToTestScripts@"

import sys

import pathlib
import datetime
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd
import ast 
from matplotlib.cm import get_cmap
import matplotlib.lines as mlines
from matplotlib.ticker import FormatStrFormatter
import matplotlib.ticker as mtick
import math

markers = ["o","v","^","<",">","1","s","p","P","*","H","x","d","|","+"]
outputFilename = "outGraph"

#This function is used to get the number of zeroes for the smallest x-value. 
#This enables the x-ticks length to be reduced by writing them as multiplications of 10^(exponent).
def getMinExponent(xStepValue):
    exponent = 1
    while xStepValue % 10 == 0 and xStepValue > 0:
        exponent *= 10
        xStepValue /= 10
    return exponent

def addPlotData(data, name, y_values, x_values):
    d = {}
    for i in range (0, len(x_values)):
        d[x_values[i]] = y_values[i]
    data[name] = d

def getPlotData(files, fileNames, xSteps):
    dataWall = {}
    dataCPU = {}
    for i in range(0,len(files)):
        f = open(files[i], "r")
        lines = f.read().splitlines()
        f.close()
        runTimes = ast.literal_eval(lines[0]) 
        wallSeconds = [item[1] for item in runTimes]
        cpuSeconds = [item[0] for item in runTimes]
        addPlotData(dataWall, fileNames[i], wallSeconds, xSteps)
        addPlotData(dataCPU, fileNames[i], cpuSeconds, xSteps)
    return dataWall, dataCPU

def getPlotColors():
    colorScheme = "Dark2"
    cmap = get_cmap(colorScheme)  # type: matplotlib.colors.ListedColormap
    return cmap.colors 

def setBarPlotLegends(data, colors, counterStart = 0):
    counter = counterStart
    legendColors = {}
    for keyName in data.keys():
        legendColors[keyName] = colors[counter]
        counter += 1      
    labels = list(legendColors.keys())
    handles = [plt.Rectangle((0,0),1,1, color=legendColors[label]) for label in labels]
    labels = [label.replace("_", " ") for label in labels]
    plt.legend(handles, labels, loc="upper left")

def setLinePlotLegends(data, colors):
    counter = 0
    legendColors = {}
    legendMarkers = {}
    for keyName in data.keys():
        legendColors[keyName] = colors[counter]
        legendMarkers[keyName] = markers[counter]
        counter += 1      
    labels = list(legendColors.keys())
    handles = [plt.Line2D([0,0],[0,1], color=legendColors[label], marker=legendMarkers[label], linestyle='-') for label in labels]
    labels = [label.replace("_", " ") for label in labels]
    plt.legend(handles, labels)

def setDescriptiveText(minExponent, yTitleName):
    plt.xlabel('Number of PSMs x$10^{}$'.format(len(str(minExponent))-1), fontsize='x-large')
    plt.ylabel(yTitleName, fontsize='x-large')

def prepareBarChart(data):
    colors = getPlotColors()
    df = pd.DataFrame(data)

    ax = df.plot(kind='bar', width=0.4, color=colors, zorder=3)
    ax.grid(zorder=0, axis='y')

    setBarPlotLegends(data, colors)
    plt.setp(plt.xticks()[1], rotation=0)

def renderLines(data, colors):
    dataFrames = []
    for keyName in data.keys():
        tmpData = {keyName: data[keyName]}
        dataFrames.append(pd.DataFrame(tmpData))
    prevAx = None
    counter = 0
    for df in dataFrames:
        if prevAx is not None:
            prevAx = df.plot(color=colors[counter], marker=markers[counter], ax=prevAx)
        else:
            prevAx = df.plot(color=colors[counter], marker=markers[counter])
        counter = (counter + 1) % min(len(markers), len(colors))

def formatHorizontalLine(data):
    ax = plt.gca()
    ax.xaxis.set_major_formatter(FormatStrFormatter('%g'))
    numXBins = len(data[next(iter(data))])
    plt.locator_params(axis = 'x', nbins = numXBins)

def prepareLineChart(data):
    colors = getPlotColors()
    renderLines(data, colors)

    setLinePlotLegends(data, colors)
    plt.setp(plt.xticks()[1], rotation=0)
    plt.grid(True)

    formatHorizontalLine(data)
    return plt

def getFiles():
    files = sys.argv[1:]
    numFiles = len(files)
    fileNames = [Path(item).stem for item in files]
    print("Files: " + str(fileNames))
    return files, fileNames

def getXValues(file):
    f = open(file, "r")
    lines = f.read().splitlines()
    f.close()
    xSteps = ast.literal_eval(lines[1])
    minExponent = getMinExponent(xSteps[0])
    xSteps = [int(item/minExponent) for item in xSteps]
    return xSteps, minExponent

def savePlot(name):
    plt.tight_layout()
    plt.savefig("{}/{}.png".format(pathToOutputData, name), dpi=200)
    plt.close()

def makeBarChart(data, minExponent, yTitleName, outputFilename):
    prepareBarChart(data)
    setDescriptiveText(minExponent, yTitleName)
    savePlot(outputFilename + "Bar")

def makeLineChart(data, minExponent, yTitleName, outputFilename):
    prepareLineChart(data)
    setDescriptiveText(minExponent, yTitleName)
    savePlot(outputFilename + "Line")

def getMaxAndMinPercentage(data):
    minVal = sys.maxsize
    maxVal = -sys.maxsize
    for value in data.values():
        percentages = list(value.values())
        maxVal = max(maxVal, math.ceil(max(percentages)))
        minVal = min(minVal, min(percentages))
    return minVal, maxVal

def nextClosePowerOf2(x):  
    if x >= 1:
        x = int(x)
        return 1 if x == 0 else 1 << (x-1).bit_length()

    minDivisionOf2 = 1
    while minDivisionOf2 >= x:
        minDivisionOf2 /= 2
    return minDivisionOf2

def setRelativePercentages(data, referenceValues):
    for key in data.keys():
        relativeValues = data[key].values()
        percentages = [y/x for x, y in zip(referenceValues, relativeValues)] #Assume it's impossible for any run-time to be 0.
        percentages = [ round(elem, 3) for elem in percentages ]
        data[key].update(zip(data[key], percentages))

def addRelativeLines(data, xTicks, colors):
    counter = 1
    for key in data.keys():
        plt.semilogy(xTicks, data[key].values(), color=colors[counter], marker=markers[counter], base=2)
        counter += 1

def extractReferenceValuesAndXTicks(tmpData):
    dataIter = iter(tmpData) 
    firstKey = next(dataIter)
    referenceValues = tmpData[firstKey].values()
    xTicks = list(tmpData[firstKey].keys())
    del tmpData[firstKey]
    return referenceValues, xTicks

def isPowerOf2(n): #https://stackoverflow.com/a/57025941
    n = int(n)
    return (n & (n-1) == 0) and n != 0

def plotIntermediateGridlinesSemiLogYBase2(min_value, max_value):
    plt.axhline(y=1, color='silver', linestyle='--', linewidth=0.5, zorder=0)
    yLineLocation = 1
    multiplier = 1
    while yLineLocation <= max_value:
        yLineLocation += 0.1*multiplier
        yLineLocation = round(yLineLocation, 6)
        plt.axhline(y=yLineLocation, color='silver', linestyle='--', linewidth=0.5, zorder=0)
        if yLineLocation.is_integer() and isPowerOf2(yLineLocation):
            multiplier = 1 << (multiplier).bit_length()
    yLineLocationUpwards = 1
    yLineLocationDownwards = 1
    multiplier = 1
    divisor = 2
    #The scaling changes for numbers below 1. For values higher than 1 the scale goes in terms of 100% but for values below 1, like the range 1 to 0.5, "50%" becomes 100%.
    #Thus, we need more spacing between lines to show margins with 10% spacing. This is why we double the distance with "0.2" instead of "0.1" like it's done above.
    while yLineLocationDownwards >= min_value:
        yLineLocationUpwards += 0.2*multiplier
        yLineLocationDownwards -= (0.2*multiplier)/divisor
        yLineLocationUpwards = round(yLineLocationUpwards, 6)
        yLineLocationDownwards = round(yLineLocationDownwards, 6)
        plt.axhline(y=yLineLocationDownwards, color='silver', linestyle='--', linewidth=0.5, zorder=0)
        if yLineLocationUpwards.is_integer() and isPowerOf2(yLineLocationUpwards):
            multiplier = 1 << (multiplier).bit_length()
            divisor = 1 << (divisor*2).bit_length()


def makeRelativeGraph(data, minExponent, yTitleName, outputFilename):
    colors = getPlotColors()
    tmpData = data.copy() #Assume first file is used as a "reference" for other files. I.e. how speed changes relative to the first file.
    referenceValues, xTicks = extractReferenceValuesAndXTicks(tmpData)
    setRelativePercentages(tmpData, referenceValues)

    min_value, max_value = getMaxAndMinPercentage(tmpData)
    min_value = nextClosePowerOf2(min_value)
    max_value = nextClosePowerOf2(max_value)
    max_value = max(max_value, 2)
    min_value = min(min_value, 0.5)

    plotIntermediateGridlinesSemiLogYBase2(min_value, max_value)
    ax = plt.gca()
    ax.grid(True, which="both", zorder=0)
    setBarPlotLegends(tmpData, colors, 1)
    addRelativeLines(tmpData, xTicks, colors)
    ax.set_ylim(bottom=min_value, top=max_value)

    setDescriptiveText(minExponent, yTitleName)
    ax.yaxis.set_major_formatter(mtick.ScalarFormatter())
    ax.yaxis.get_major_formatter().set_scientific(False)
    ax.yaxis.get_major_formatter().set_useOffset(False)
    ax.yaxis.set_major_formatter(FormatStrFormatter('%g'))
    ax.set_xticks(xTicks)
    plt.setp(plt.xticks()[1], rotation=0)

    savePlot(outputFilename + "Relative")

    
files, fileNames = getFiles()
xSteps, minExponent = getXValues(files[0]) #Assumes all plots share the same step size for the horizontal line, x.
dataWall, dataCPU = getPlotData(files, fileNames, xSteps)

makeBarChart(dataWall, minExponent,  'Wall Time in Seconds', outputFilename + "Wall")
makeBarChart(dataCPU, minExponent,   'CPU Time in Seconds', outputFilename + "CPU")
makeLineChart(dataWall, minExponent, 'Wall Time in Seconds', outputFilename + "Wall")
makeLineChart(dataCPU, minExponent,  'CPU Time in Seconds', outputFilename + "CPU")

makeRelativeGraph(dataWall, minExponent, "Fraction of Original {} Time".format("Wall"), outputFilename + "Wall")
makeRelativeGraph(dataCPU, minExponent, "Fraction of Original {} Time".format("CPU"), outputFilename + "CPU")



