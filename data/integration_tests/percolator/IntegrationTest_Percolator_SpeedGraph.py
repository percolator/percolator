# Marcus Andersson - Percolator Project
# Takes data from IntegrationTest_Percolator_SpeedGraphData.py and outputs graph images.
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


markers = ["o","v","^","<",">","1","s","p","P","*","H","x","d","|","+"]
outputFilename = "outGraph"

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

def setBarPlotLegends(data, colors):
    colorCounter = 0
    legendColors = {}
    for keyName in data.keys():
        legendColors[keyName] = colors[colorCounter]
        colorCounter += 1      
    labels = list(legendColors.keys())
    handles = [plt.Rectangle((0,0),1,1, color=legendColors[label]) for label in labels]
    plt.legend(handles, labels)

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
    plt.legend(handles, labels)

def setDescriptiveText(minExponent, titleName):
    plt.suptitle(titleName, fontsize='xx-large')
    plt.xlabel('Number of PSM Lines x$10^{}$'.format(str(len(str(minExponent))-1)), fontsize='x-large')
    plt.ylabel('Time in Seconds', fontsize='x-large')

def getBarChart(data):
    colors = getPlotColors()
    df = pd.DataFrame(data)

    ax = df.plot(kind='bar', width=0.4, color=colors, zorder=3)
    ax.grid(zorder=0)

    setBarPlotLegends(data, colors)
    plt.setp(plt.xticks()[1], rotation=0)
    return plt

def renderLines(data, colors):
    dataFrames = []
    for keyName in data.keys():
        tmpData = {keyName: data[keyName]}
        print(tmpData)
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

def getLineChart(data):
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

def makeBarChart(data, minExponent, titleName, outputFilename):
    chart = getBarChart(data)
    setDescriptiveText(minExponent, titleName)
    plt.tight_layout()
    chart.savefig("{}/{}.png".format(pathToOutputData, outputFilename + "Bar"))
    chart.close()

def makeLineChart(data, minExponent, titleName, outputFilename):
    chart = getLineChart(data)
    setDescriptiveText(minExponent, titleName)
    plt.tight_layout()
    chart.savefig("{}/{}.png".format(pathToOutputData, outputFilename + "Line"))
    chart.close()


files, fileNames = getFiles()
xSteps, minExponent = getXValues(files[0]) #Assumes all plots share the same step size for the horizontal line, x.
dataWall, dataCPU = getPlotData(files, fileNames, xSteps)

makeBarChart(dataWall, minExponent,  'Wall Clock Speed', outputFilename + "Wall")
makeBarChart(dataCPU, minExponent,   'CPU Clock Speed', outputFilename + "CPU")
makeLineChart(dataWall, minExponent, 'Wall Clock Speed', outputFilename + "Wall")
makeLineChart(dataCPU, minExponent,  'CPU Clock Speed', outputFilename + "CPU")


