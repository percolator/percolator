# Marcus Andersson - Percolator Project
# Script that tests the performances of percolator for different data input sizes. Produces a graph image when finished.

import os
import sys
import re

from argparse import ArgumentParser
import time
import contextlib
import subprocess
from subprocess import Popen, PIPE
from decimal import Decimal
import pathlib
import datetime
from pathlib import Path
import urllib
from urllib.request import urlopen, Request
from urllib.error import URLError, HTTPError
import tarfile

import matplotlib.pyplot as plt
import pandas as pd


success = True


scriptDirectory = pathlib.Path(__file__).parent.absolute()
testScriptPath = str(scriptDirectory) + "/" + "IntegrationTest_Percolator_Speed.py"
testDataFolderPath = str(scriptDirectory) + "/" + "testData"
testDataFilePath = ""

path = Path(testDataFolderPath)
path.mkdir(parents=False, exist_ok=True)

def getArguments():
    parser = ArgumentParser(description="Measure speed of perculator application with varying input and creates a graph.")
    parser._action_groups.pop()
    required = parser.add_argument_group('required arguments')
    optional = parser.add_argument_group('optional arguments')

    required.add_argument('-d','--data', type=str, metavar='', required=True, help="Path to input-data used by percolator.")
    optional.add_argument('-py','--python_version', type=str, metavar='', default="python3" , required=False, help="Python version to call other python modules from.")
    optional.add_argument('-r','--runs', type=int, metavar='', default=3, required=False, help="Number of tests to evaluate.")
    optional.add_argument('-f','--flags', type=str, default = "" , metavar='', required=False, help="Flags to use with the executable.")
    optional.add_argument('-a','--await_tests', default=False, action="store_true", required=False, help="Wait for others processes running this script to finish before starting.")
    optional.add_argument('-c','--comments', type=str, default="", metavar='', required=False, help="Comments regarding test details.")
    optional.add_argument('-fname','--filename', type=str, default="barFigure", metavar='', required=False, help="Name of produced charts.")
    return parser.parse_args()

def makeReducedFile(fileNameOriginal, fileNameReduced, numLines):
    fileReduced = open(fileNameReduced, "w")
    counter = 0
    with open(fileNameOriginal) as f:
        for line in f:
            counter += 1
            if (counter > numLines):
                break
            fileReduced.write(line)
    fileReduced.close()

from itertools import (takewhile,repeat)
def rawincount(filename): ####https://stackoverflow.com/a/27518377/7202012
    f = open(filename, 'rb')
    bufgen = takewhile(lambda x: x, (f.raw.read(1024*1024) for _ in repeat(None)))
    return sum( buf.count(b'\n') for buf in bufgen )

def getTemporaryInDataFileName(number):
    return "psmInData" + str(number) + ".tab"

def getTerminalCommand(args, data):
    result = args.python_version
    result += " " + testScriptPath
    result += " -d " + data
    if len(args.flags) > 0:
        result += " -f \"" + args.flags + "\""
    result += " -r " + str(args.runs)
    if args.await_tests is True:
        result += " -a "
    print("RES: " + result)
    return result
    
def addElement(dict, key, value):
    if key not in dict:
        dict[key] = []
    dict[key].append(value)

def getChartLinesShortened(tmpNumLines):
    tmpNumLines /= 1000000
    return int(tmpNumLines)

def getBarChart(y_valuesWall, y_valuesCPU, x_numLines):
    d1 = {}
    d2 = {}
    data = {}
    for i in range (0,len(x_numLines)):
        d1[x_numLines[i]] = y_valuesWall[i]
        d2[x_numLines[i]] = y_valuesCPU[i]
    data['Wall'] = d1
    data['CPU'] = d2
    df = pd.DataFrame(data)
    df.plot(kind='bar', width=0.4, color=['darkgray','gray'])
    colors = {'Wall':'darkgray', 'CPU':'gray'}         
    labels = list(colors.keys())
    handles = [plt.Rectangle((0,0),1,1, color=colors[label]) for label in labels]
    plt.legend(handles, labels)
    plt.tight_layout()
    plt.setp(plt.xticks()[1], rotation=0)
    return plt


args = getArguments()
testDataFilePath = args.data

#################
numLines = 5000000
numDataFiles = 0
numLinesMax = rawincount(testDataFilePath)
for i in range(0, int(numLinesMax/numLines)):
    numDataFiles += 1
    makeReducedFile(testDataFilePath, testDataFolderPath + "/" + getTemporaryInDataFileName(i+1), numLines*numDataFiles)
#################
runTimes = []
for i in range(1, numDataFiles+1):
    terminalCommand = getTerminalCommand(args, testDataFolderPath + "/" + getTemporaryInDataFileName(i))
    p = subprocess.Popen(terminalCommand, stderr=PIPE , stdout=PIPE , text=True, shell=True)
    print("Running subprocess: " + str(p.pid))
    p.wait()
    result = p.stdout.read()
    print(p.stderr.read())
    result = result.splitlines()
    for item in result:
        if "Min" in item:
            runTime = re.findall(r"[-+]?\d*\.\d+|\d+", item)
            runTime = [float(i) for i in runTime]
            runTimes.append(runTime)
#################
y_valuesWall = []
y_valuesCPU = []
x_numLines = []
for i in range(1, numDataFiles+1):
    numHandledLines = 0
    if(i  == numDataFiles):
        numHandledLines = numLinesMax
    else:
        numHandledLines = i*numLines
    x_numLines.append(numHandledLines)
    wallSec = runTimes[i-1][1]
    cpuSec = runTimes[i-1][0]
    y_valuesWall.append(wallSec)
    y_valuesCPU.append(cpuSec)
x_numLines = [getChartLinesShortened(item) for item in x_numLines]
##############################
chart = getBarChart(y_valuesWall, y_valuesCPU, x_numLines)
chart.savefig(str(scriptDirectory) + "/" + args.filename + ".png")
chart.close()


if success==True:
 print("...TEST SUCCEEDED")
 exit(0)
else:
 print("...TEST FAILED")
 exit(1)