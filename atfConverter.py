#!/usr/bin/env python3
# Names: Joseph Olivier 
# Program Description: This program will look in a directory with an arbitrary number of ATF files. It will group and plot the different 
# the different files according to user preference. Right now, plots are only for square-wave voltammetry because that is what we primarily
# use in the lab. It will plot lines with different colors. It will create a folder according to current datetime and save graphs there.
# After looking at derivative graphs, I decided not to include 2nd derivative graphs. This is because I felt like they did not produce anymore
# useful information than did the 1st derivative graphs.
#
# Improvements over Current Practice: Currently, all graphs are generated through excel. In addition, in order to find current information,
# the current parsing program looks at specific timepoints in the square-wave graph in order to derive the difference between the top of the
# square-wave and the bottom of the square-wave. However, there are oftentimes "spikes" in data which can lead to extremely inaccurate measurements
# when using this method. My new data parsing method takes the median value over a range of 100 in order to ensure that a spike is not selected.
#
# Future Directions: This progam will be a work in progress for some time. I will be adding a GUI so a user is able to input the exact 
# program parameters in place on the command line. This will allow users who are not experts at the command line to use this program.
# I will also add support for different types of waves. Right now, default support is for square wave voltammetry as that is what we use.
# Estimated Completion with all functionality: several months
import matplotlib
matplotlib.use('Agg')
class CommandLine():
    def __init__(self, inOpts=None):
        import argparse
        self.parser = argparse.ArgumentParser(description="This program is intended to generate graphs based on the Voltage : Current of a nanopipette in solution.")
        self.parser.add_argument('-hz', '--Hertz', default=100,help="Amount of measurements per second with a default of 100")
        self.parser.add_argument('-s','--Steps',default=9,help="Amount of steps in program. Default = 8")
        self.parser.add_argument('-p', '--PulseWidth',default=100,help="Sets pulse width of steps. Default = 100 ms")
        self.parser.add_argument('-fd','--FirstDuration',default=1000,help="Sets the duration of each pulse at each Voltage. Default = 1000 ms")
        self.parser.add_argument('r','--Repeats',default=7,help="Number of repeats for each well. Default = 7")
        if inOpts is None:
            self.args = self.parser.parse_args()
        else:
            self.args = self.parser.parse_args(inOpts)

class GraphMaker:
	
	def __init__(self,listofCurrents,graphOutput):
		self.masterList = listofCurrents
		self.graphOutput = graphOutput
		self.AverageList = []
		self.allYPoints = []
		self.outputList = []
	
	def grapher(self):
		import matplotlib.pyplot as plt
		import sys
		import numpy as np
		
		
		plt.figure(figsize=(4,2))

		#Blank and sample, to be changed in the future
		numberOfGroups = 2
		panel1=plt.axes([0.2,0.2,1/4,1/2])
		panel1.tick_params(axis='both',which='both',\
                   bottom='on', labelbottom='on',\
                   left='on', labelleft='on',\
                   right='off', labelright='off',\
                   top='off', labeltop='off')
		panel2=plt.axes([0.6,0.2,1/4,1/2])
		panel2.tick_params(axis='both',which='both',\
                   bottom='on', labelbottom='on',\
                   left='on', labelleft='on',\
                   right='off', labelright='off',\
                   top='off', labeltop='off')
		overallMax = -float('inf')
		overallMin = float('inf')
		derivativeMax = -float('inf')
		derivativeMin = float('inf')
		firstDerivativeX = [-1000,-800,-600,-400,-200,0,200,400]
		numberOfFiles = len(self.masterList)
		lengthOfGroups = numberOfFiles/numberOfGroups
		#Iterate through each file in the listOfFiles
		averageY = [[] for i in range(8)]
		averageYDeriv = [[] for i in range(8)]
		
		for h in range(0,len(self.masterList)):
			print(h)
			if h != 0:
				prevXVals = xVals
				prevYVals = yVals
			xVals = []
			yVals = []
			file = self.masterList[h]
			firstDerivativeY = [0]
			for i in range(0,len(file)-1):
				tple = file[i]
				nextTple = file[i+1]
				#Normalize value due to step size
				currentVal = file[i][2]*(50/abs(tple[1][1]-tple[0][1]))
				if i != 0:
					print(h,i,file[i-1][2],file[i-1])
					firstDerivativeY.append(currentVal/((file[i-1][2]+.01)*(200/abs(tple[1][1]-tple[0][1])))-1)
				xVals.append(tple[1][1])
				yVals.append(currentVal)
				#Normalization to make sure equal units between points
				if currentVal<overallMin:
					overallMin = currentVal
				elif currentVal>overallMax:
					overallMax = currentVal
				if firstDerivativeY[-1]<derivativeMin:
					derivativeMin = firstDerivativeY[-1]
				elif firstDerivativeY[-1]>derivativeMax:
					derivativeMax = firstDerivativeY[-1]
				averageY[i].append(currentVal)
				if i == 0:
					averageYDeriv[i].append(0)
				if i != 0:
					averageYDeriv[i].append(currentVal/((file[i-1][2]+0.01)*(200/abs(tple[1][1]-tple[0][1])))-1)
				panel1.plot(tple[1][1],currentVal,'-o',ms=2,mfc='black', mew=0,linewidth=1,zorder=500)
				#print(tple[1][1],currentVal)
				#print(file[i][2])
				panel2.plot(firstDerivativeX[i],firstDerivativeY[-1],'-o',ms=2,mfc='black', mew=0,linewidth=1,zorder=500)

				if i == len(file)-2:
					saveColor1 = h
					nextVal = nextTple[2]*(50/abs(nextTple[1][1]-nextTple[0][1]))
					firstDerivativeY.append(nextVal/(file[i-1][2]*(200/abs(tple[1][1]-tple[0][1])))-1)
					panel1.plot(nextTple[1][1],nextVal,'-o',ms=2,mfc='black', mew=0,linewidth=1,zorder=500)
					panel2.plot(firstDerivativeX[i+1],firstDerivativeY[-1],'-o',ms=2,mfc='black', mew=0,linewidth=1,zorder=500)
					xVals.append(nextTple[1][1])
					yVals.append(nextVal)
					averageY[i+1].append(nextVal)
					averageYDeriv[i+1].append(nextVal/(file[i-1][2]*(200/abs(tple[1][1]-tple[0][1])))-1)
					panel1.plot(xVals,yVals,color=(.1*(h+1)%1,1/(h+1)%1,.1*(h+1)%1,.5))

					panel2.plot(firstDerivativeX,firstDerivativeY,color=(.1*(h+1)%1,1/(h+1)%1,.1*(h+1)%1,.5))
					if currentVal<overallMin:
						overallMin = currentVal
					elif currentVal>overallMax:
						overallMax = currentVal
					if firstDerivativeY[-1]<derivativeMin:
						derivativeMin = firstDerivativeY[-1]
					elif firstDerivativeY[-1]>derivativeMax:
						derivativeMax = firstDerivativeY[-1]
			
			#if h == lengthOfGroups-1:
			"""TO DO: Make a helper function to combine all of these lines"""
			if h == 1:
				newAverageY = []
				newAverageDerivY = []
				for i in range(len(averageY)):
					newAverageY.append(np.mean(averageY[i]))
					newAverageDerivY.append(np.mean(averageYDeriv[i]))
				self.AverageList.append((newAverageY,newAverageDerivY))
				self.allYPoints.append((averageY,averageYDeriv))

				panel1.set_xlim(-620,-180)
				panel1.set_ylim(overallMin-50,overallMax+50)
				panel1.set_xlabel("Voltage (mV)")
				panel1.set_ylabel("Current (pA)")
				#panel1.set_xticks([-1000,-600,-100,400])
				panel1.set_title("Blank Current to Voltage Graph",size=7)
				panel1.grid(b=True, which='major', color='k', linestyle='-')
				panel1.grid(b=True, which='minor', color='r', linestyle='-', alpha=0.2)
				panel1.minorticks_on()
				
				panel2.set_xlim(-620,-180)
				panel2.set_ylim(overallMin-50,overallMax+50)
				panel2.set_xlabel("Voltage (mV)")
				panel2.set_ylabel("Derivative")
				panel2.set_title("First Derivatives",size=7)
				#panel2.set_xticks([-1000,-600,-100,400])
				panel2.grid(b=True, which='major', color='k', linestyle='-')
				panel2.grid(b=True, which='minor', color='r', linestyle='-', alpha=0.2)
				panel2.minorticks_on()
				plt.savefig(self.graphOutput+"/data_blank.pdf")
				plt.clf()
				averageY = [[] for i in range(8)]
				averageYDeriv = [[] for i in range(8)]

				panel1=plt.axes([0.2,0.2,1/4,1/2])
				panel1.tick_params(axis='both',which='both',\
                   bottom='on', labelbottom='on',\
                   left='on', labelleft='on',\
                   right='off', labelright='off',\
                   top='off', labeltop='off')
				panel2=plt.axes([0.6,0.2,1/4,1/2])
				panel2.tick_params(axis='both',which='both',\
                   bottom='on', labelbottom='on',\
                   left='on', labelleft='on',\
                   right='off', labelright='off',\
                   top='off', labeltop='off')

			#elif (h-1) % lengthOfGroups == 0:
			elif h==3:
				newAverageY = []
				newAverageDerivY = []
				for i in range(len(averageY)):
					newAverageY.append(np.mean(averageY[i]))
					newAverageDerivY.append(np.mean(averageYDeriv[i]))
				self.AverageList.append((newAverageY,newAverageDerivY))
				self.allYPoints.append((averageY,averageYDeriv))
				panel1.set_xlim(-620,-180)
				panel1.set_ylim(overallMin-50,overallMax+50)
				panel1.set_xlabel("Voltage (mV)")
				panel1.set_ylabel("Current (pA)")
				#panel1.set_xticks([-1000,-600,-100,400])
				panel1.set_title("Sample Current to Voltage Graph",size=7)
				panel1.grid(b=True, which='major', color='k', linestyle='-')
				panel1.grid(b=True, which='minor', color='r', linestyle='-', alpha=0.2)
				panel1.minorticks_on()

				panel2.set_xlim(-620,-180)
				panel2.set_ylim(overallMin-50,overallMax+50)
				panel2.set_xlabel("Voltage (mV)")
				panel2.set_ylabel("Derivative")
				panel2.set_title("First Derivatives",size=7)
				#panel2.set_xticks([-1000,-600,-100,400])
				panel2.grid(b=True, which='major', color='k', linestyle='-')
				panel2.grid(b=True, which='minor', color='r', linestyle='-', alpha=0.2)
				panel2.minorticks_on()
				plt.savefig(self.graphOutput+"/data_sample.pdf")
				plt.clf()

				averageY = [[] for i in range(8)]
				averageYDeriv = [[] for i in range(8)]
				panel1=plt.axes([0.2,0.2,1/4,1/2])
				panel1.tick_params(axis='both',which='both',\
                   bottom='on', labelbottom='on',\
                   left='on', labelleft='on',\
                   right='off', labelright='off',\
                   top='off', labeltop='off')
				panel2=plt.axes([0.6,0.2,1/4,1/2])
				panel2.tick_params(axis='both',which='both',\
                   bottom='on', labelbottom='on',\
                   left='on', labelleft='on',\
                   right='off', labelright='off',\
                   top='off', labeltop='off')
				panel1.grid(True, which='both')
				panel1.minorticks_on
				panel2.grid(True, which='both')
				panel2.minorticks_on
		
		panel1.set_xlim(-620,-180)
		panel1.set_ylim(overallMin-50,overallMax+50)
		panel1.set_xlabel("Voltage (mV)")
		panel1.set_ylabel("Current (pA)")
		panel1.set_title("Combined Current to Voltage",size=7)
		#panel1.set_xticks([-1000,-600,-100,400])
		panel1.grid(b=True, which='major', color='k', linestyle='-')
		panel1.grid(b=True, which='minor', color='r', linestyle='-', alpha=0.2)
		panel1.minorticks_on()

		panel2.set_xlim(-620,-180)
		panel2.set_ylim(derivativeMin-.1,derivativeMax+.1)
		panel2.set_xlabel("Voltage (mV)")
		panel2.set_ylabel("Derivative")
		panel2.set_title("Combined First Derivatives",size=7)
		#panel2.set_xticks([-1000,-600,-100,400])
		panel2.grid(b=True, which='major', color='k', linestyle='-')
		panel2.grid(b=True, which='minor', color='r', linestyle='-', alpha=0.2)
		panel2.minorticks_on()
		xValues = [-600,-550,-500,-450,-400,-350,-300,-250]

		for i in range(len(self.AverageList)):
			currentList = self.AverageList[i]
			listofY = []
			listofYDeriv = []
			for j in range(len(currentList[0])):
				listofY.append(currentList[0][j])
				listofYDeriv.append(currentList[1][j])
			if i == 1:
				panel1.plot(xValues,listofY,color=(.1*(1+saveColor1)%1,1/(1+saveColor1)%1,.1*(1+saveColor1)%1,.5))
				panel2.plot(xValues,listofYDeriv,color=(.1*(saveColor1+1)%1,1/(saveColor1+1)%1,.1*(saveColor1+1)%1,.5))
			else:
				panel1.plot(xValues,listofY,color=(.1*(1+i%1),1/(1+i)%1,.1*(1+i)%1,.5))
				panel2.plot(xValues,listofYDeriv,color=(.1*(i+1)%1,1/(i+1)%1,.1*(i+1)%1,.5))
		lowestMin = float('inf')
		highestMax = -float('inf')
		
		for i in range(len(self.allYPoints)):
			currentList = self.allYPoints[i]
			for j in range(len(currentList)):
				currentGraph = currentList[j]
				pan1 = currentList[0]
				pan2 = currentList[1]
				for k in range(len(pan1)):
					amountOfDots = len(pan1[0])
					val1 = ""
					val2 = ""
					val11 = ""
					val21 = ""
					for l in range(amountOfDots):
						currentDots = pan1[l]
						if pan1[k][l] > highestMax:
							highestMax = pan1[k][l]
						elif pan1[k][l] < lowestMin:
							lowestMin = pan1[k][l]
						panel1.plot(xValues[k],pan1[k][l],'-o',ms=2,mfc='black', mew=0,linewidth=1,zorder=500)
						panel2.plot(xValues[k],pan2[k][l],'-o',ms=2,mfc='black', mew=0,linewidth=1,zorder=500)
						if l == 0:
							val1 = pan1[k][l]
							val2 = pan2[k][l]
						else:
							val11 = pan1[k][l]
							val21 = pan2[k][l]
					panel1.plot([xValues[k],xValues[k]],[val1,val11],color=(0,0,0,.5))
					panel2.plot([xValues[k],xValues[k]],[val2,val21],color=(0,0,0,.5))
		panel1.set_ylim(lowestMin-50,highestMax+50)
		plt.savefig(self.graphOutput+"/data_combined.pdf")

class DataParser:

	def __init__(self,listofFiles,graphOutput):
		self.graphOutput = graphOutput
		self.masterList = listofFiles

	def parser(self):
		import numpy as np
		listofCurrents = []

		#Iterates through all of the files that is in the currentWorkingDirectory/ATFFiles
		counter = 0
		for wholeFile in self.masterList:
			internalFileList = []
			#Assuming 9 steps, will change with GUI. 
			for i in range(0,8):
				firstStart = 1411
				PulseWidth=100
				#Value between FULL steps
				offset = 1000
				#Finds position of median in sorted list (depends on pulse width)
				median = int((firstStart+PulseWidth-firstStart)/2)
				#Sorts the currents at the top of the square wave
				sortedbyCurrentTop = sorted(wholeFile[firstStart+offset*i:(firstStart+100)+offset*i],key=lambda current:current[1])
				#Sorts the currents at the bottom of the square wave
				sortedbyCurrentBottom = sorted(wholeFile[firstStart+offset*i+100:(firstStart+100)+offset*i+150],key=lambda current:current[1])
				differenceInStep = sortedbyCurrentTop[median][1] - sortedbyCurrentBottom[median][1]
				internalFileList.append((sortedbyCurrentTop[int(PulseWidth/2)][1:],sortedbyCurrentBottom[int(PulseWidth/2)][1:],differenceInStep))
			listofCurrents.append(internalFileList)
		GraphMaker(listofCurrents,self.graphOutput).grapher()

class FileReader:
	def __init__(self,fileToRead):
		self.file = open(fileToRead)
	
	def reader(self):
		listofLines = []
		#Skip the header lines of the file
		for i in range(0,11):
			self.file.readline()
		#Convert all items to floats in file instead of strings.
		for line in self.file:
			line = line.strip().split()
			line[0] = float(line[0])
			line[1] = float(line[1])
			line[2] = float(line[2])
			listofLines.append(line)
		return listofLines

def main():
	#TO-DO:
	#Split FileReader and main function into separate file.
	import os
	import datetime
	#Location where the ATF files should be located (.../CWD/ATFFiles)
	refATFin = os.getcwd() + '/ATFFiles'
	#Location where the output directory will be located (.../CWD/graphs_current_date_to_second)
	graphOutput = os.getcwd() + '/graphs' + datetime.datetime.now().strftime('_%Y_%b_%d_%H_%M_%S')
	os.makedirs(graphOutput,0o777)
	FileList = []
	#Iterates through the ATF File folder and find all of the names of the ATF files.
	#Adds ATF file names to list, parses data using graphOutput Function.
	for dirname, dirnames, filenames in os.walk(refATFin):
		for filename in filenames:
			startfile = os.path.join(refATFin,filename)
			currentFile = FileReader(startfile).reader()
			FileList.append(currentFile)
	#Send list of files to dataParser
	DataParser(FileList,graphOutput).parser()

main()