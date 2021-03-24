# Marcus Andersson - Percolator Project
# Script that tests the performances of percolator for different data input sizes. Produces a file containing data for graphs.
# psutil is imported in the case that the script should wait for other integration-test scripts to finish.

pathToOutputData = "/home/ekvall/percolator/tests/integration_tests/percolator/testData"
pathToTestScripts = "/home/ekvall/percolator/tests/integration_tests/percolator"

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

import matplotlib.pyplot as plt
import pandas as pd
from itertools import (takewhile,repeat)

testScriptPath = pathToTestScripts + "/IntegrationTest_Percolator_Speed.py"
#The amount of input lines that also determine the step size between tests. Tests are made with increasing input size determined by number of PSMs.
minNumLines = -1
processText = 'IntegrationTest_Percolator_SpeedGraphData.py'


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
    optional.add_argument('-xstep','--x_step_size', type=int, default=50000, metavar='', required=False, help="Amount of input lines that also determine the step size between tests.")
    return parser.parse_args()

def makeReducedFile(fileNameOriginal, fileNameReduced, numLines):
    if Path(fileNameReduced).is_file():
        return
    fileReduced = open(fileNameReduced, "w")
    counter = 0
    with open(fileNameOriginal) as f:
        for line in f:
            counter += 1
            if counter > numLines:
                break
            fileReduced.write(line)
    fileReduced.close()

def countLinesInFile(filename): ####https://stackoverflow.com/a/27518377/7202012
    f = open(filename, 'rb')
    bufgen = takewhile(lambda x: x, (f.raw.read(1024*1024) for _ in repeat(None)))
    return sum( buf.count(b'\n') for buf in bufgen )

def getTemporaryInDataFileName(number, dataFilename):
    return dataFilename + str(number) + ".tab"

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

def generateInputFiles(args):
    maxNumLines = countLinesInFile(args.data)
    dataFilename = Path(args.data).stem
    numDataFiles = 1

    intialNumLines = getInitialNumLines()
    makeReducedFile(args.data, "{}/{}".format(pathToOutputData, getTemporaryInDataFileName(intialNumLines, dataFilename)), intialNumLines)

    for i in range(0, int(maxNumLines/minNumLines)):
        fileNumber = minNumLines*numDataFiles
        numDataFiles += 1
        makeReducedFile(args.data, "{}/{}".format(pathToOutputData, getTemporaryInDataFileName(fileNumber, dataFilename)), fileNumber)
    if numDataFiles == 0:
        print("Error, data input file must contain at least {} lines of PSMs.".format(minNumLines))
        exit(1)
    return numDataFiles

def runCommandAndWait(terminalCommand, runTimes):
    p = subprocess.Popen(terminalCommand, stderr=PIPE , stdout=PIPE , text=True, shell=True)
    print("Running subprocess: " + str(p.pid))
    p.wait()
    result = p.stdout.read()
    print(p.stderr.read())
    result = result.splitlines()
    for item in result:
        if "Min" in item:
            #This regex script depends on the output from IntegrationTest_Percolator_Speed.py to be in a very specific format.
            runTime = re.findall(r"[+]?\d*\.\d+|\d+", item)
            runTime = [float(i) for i in runTime]
            runTimes.append(runTime)

def getInitialNumLines():
    return 10**(len(str(minNumLines)) -1 )

def getRunTimes(args):
    runTimes = []
    dataFilename = Path(args.data).stem
    intialNumLines = getInitialNumLines()
    terminalCommand = getTerminalCommand(args, "{}/{}".format(pathToOutputData, getTemporaryInDataFileName(intialNumLines, dataFilename)) )
    runCommandAndWait(terminalCommand, runTimes)
    for i in range(1, numDataFiles):
        terminalCommand = getTerminalCommand(args, "{}/{}".format(pathToOutputData, getTemporaryInDataFileName(minNumLines*i, dataFilename)) )
        runCommandAndWait(terminalCommand, runTimes)
    return runTimes

def getXNumLines(numDataFiles):
    x_numLines = []
    intialNumLines = getInitialNumLines()
    x_numLines.append(intialNumLines)
    for i in range(1, numDataFiles):
        x_numLines.append(i*minNumLines)
    return x_numLines

def saveData(args, runTimes, x_numLines):
    f = open(pathToOutputData + "/" + args.filename, "w")
    f.write(str(runTimes) + '\n' )
    f.write(str(x_numLines) + '\n' )
    f.write(args.comments + '\n' )
    f.close()

#waitForProcess(pid) could cause the computer to wait indefinitely if stuck in a race condition. It is assumed that this will not occur normally however.
#Example of eternal wait scenario: A user spawns 2 processes, let's call them A and B. B waits for A to finish, but when A is finished another process spawns with the same pid as A had. B may fail to notice this.
def waitForProcess(pid):
    numIterations = 0
    delayMultiplier = 1
    while(True):
        if psutil.pid_exists(pid):
            if numIterations == 0 or numIterations > 60*delayMultiplier:
                delayMultiplier = delayMultiplier * 2
                numIterations = 1
                print("Waiting for PID {} to finish. ({})".format(pid, datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")))
            time.sleep(45)
        else:
            return
        numIterations += 1

def waitForProcesses(pids):
    if pids is None:
        return
    for p in pids:
        waitForProcess(p)

def waitForOtherTests():
    thisPid = os.getpid()
    parentPid = psutil.Process(thisPid).ppid()
    allProcesses = []
    for p in psutil.process_iter():
        if thisPid == p.pid or parentPid == p.pid:
            continue
        if processText in p.name() or processText in ' '.join(p.cmdline()):
            allProcesses.append(p.pid)
    waitForProcesses(allProcesses)

args = getArguments()
if args.await_tests is True:
    import psutil
    waitForOtherTests()

minNumLines = args.x_step_size
numDataFiles = generateInputFiles(args)

runTimes = getRunTimes(args)
x_numLines = getXNumLines(numDataFiles)
saveData(args, runTimes, x_numLines)

print("...TEST FINISHED")


