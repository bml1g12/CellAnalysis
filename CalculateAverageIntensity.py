"""
MIT License

Copyright (c) 2018 Benjamin Mark Lowe

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
"""
#======================
#Version 1.0 
#Jython Script to calculate the average intensity across many cells which have been 
#labelled using the ImageJ MorphoLibJ Segmentation plugin. 
#See seperate documentation for usage instructions.
#======================
####User Parameters####
#Note: imageJ label number can be visualised using Image>Color>Color Picker's value
#By default background is set to largest label, but can be manually set here:
manually_assign_backgroundlayer_to_label=None #(None/Int) Set to None to auto-select the background as the label with largest area.							  
warning_threshold=0.01 			# If the normalised intensity of a cell (R_a=(C1-B1)/(C2-B2)) is below this threshold, throw a warning at the end of the script. 
show_table=False 				# Show a visual summary of results (True/False)
for_each_cell_calculate="Mean" 	# ("Mean"/"Median"). Sets what method is used for calculating C1,B1,C2 and B2 WITHIN each cell
######################

from ij import IJ, WindowManager, ImagePlus
from java.util import ArrayList
from ij.gui import GenericDialog
from ij.measure import ResultsTable
from ij.plugin import PlugIn, Duplicator
from inra.ijpb.measure import IntensityMeasures
from inra.ijpb.label.LabelImages import keepLargestLabel
from inra.ijpb.label import LabelImages
import operator
from math import sqrt
import sys
import os
import csv

assert for_each_cell_calculate == "Mean" or for_each_cell_calculate == "Median", "Error: for_each_cell_calculate must be set to either 'Mean' or 'Median'"

def median(lst):
    quotient, remainder = divmod(len(lst), 2)
    if remainder:
        return sorted(lst)[quotient]
    return float(sum(sorted(lst)[quotient - 1:quotient + 1]) / 2)

def standard_deviation(lst, population=True):
    """Calculates the standard deviation for a list of numbers."""
    num_items = len(lst)
    mean = sum(lst) / num_items
    differences = [x - mean for x in lst]
    sq_differences = [d ** 2 for d in differences]
    ssd = sum(sq_differences)
 
    if population is True:
        #print('This is POPULATION standard deviation.')
        variance = ssd / num_items
    else:
        #print('This is SAMPLE standard deviation.')
        variance = ssd / (num_items - 1)
    sd = sqrt(variance)

    return sd

def load_environment(hardcoded_data_name=None, hardcoded_label_name=None):
	"""Load required files. By default uses the GUI, otherwise filepaths can be manually specified"""

	if hardcoded_data_name and hardcoded_label_name:
		print("Loading {} and {} from hard coded selection".format(hardcoded_data_name,hardcoded_label_name))
		rawimg=IJ.openImage(hardcoded_data_name)
		labelimg=IJ.openImage(hardcoded_label_name)
		rawimg_name=rawimg.getTitle()
		labelimg_name=labelimg.getTitle()
		
		print("Data filename is: "+str(rawimg_name))
		print("Labels filename is: "+str(labelimg_name))
		print("Which has these labels: {}".format( list(LabelImages.findAllLabels( labelimg ))))

	else:
		print("Loading images from GUI selection")
		nbima = WindowManager.getImageCount()
		
		print("There are this {} windows open".format(nbima))
	
		if nbima < 2 :
			IJ.error( "Intensity Measures input error", \
					"ERROR: At least two images need to be open to run "\
							+ "Intensity Measures" )
			return 0
	
		names=[]
		for i in xrange(nbima):
			names.append(WindowManager.getImage(i + 1).getTitle())
	    	
		
		gd = GenericDialog( "Measure 3D" )
		gd.addChoice( "Input Data", names, names[0] )
		gd.addChoice( "Labels", names, names[1] )
		gd.showDialog()
		if gd.wasOKed():
			inputIndex = gd.getNextChoiceIndex()
			labelsIndex = gd.getNextChoiceIndex()
	
		rawimg=WindowManager.getImage(inputIndex+1)
		labelimg=WindowManager.getImage(labelsIndex+1)
		rawimg_name=rawimg.getTitle()
		labelimg_name=labelimg.getTitle()
		print("Data filename is: "+str(rawimg_name))
		print("Which has {} time frames.".format(rawimg.getNFrames()))
		print("Labels filename is: "+str(labelimg_name))
		print("Which has these labels: {}".format( list(LabelImages.findAllLabels( labelimg ))))
		

	#largestlabelimg=keepLargestLabel(labelimg) # can be used to make a label of just the background
	if rawimg.getWidth() != labelimg.getWidth() or rawimg.getHeight() != labelimg.getHeight():
		IJ.error( "Intensity Measures 2D/3D input error", "Error: input"
				+ " and label images must have the same size" )

	return rawimg, labelimg, rawimg_name, labelimg_name

def summarise_input(rawimg, labelimg):
	"""This function takes as input: 'rawimg' (the data) and 'labelimg' (the cell boundary cartoon) 
	Then using the z=1, channel=1 frame, produces a summary table for inspection. 
	It also calculates which label is the background, assuming it is the largest cell.
	It returns the following:
	(e.g. if there are three labels labeled 17, 20, 41. Where "20" is the background. This function will return:
	1. A list of all the labels  (e.g. [17,20,41])
	2. The number of labels (e.g. 3)
	3. The index of the background (e.g. 1) 
	4. The label name of the background (e.g. 20)
	"""
	#Take a snapshot of the image at z=1 and c=1 for this analysis
	inputimg=Duplicator().run(rawimg, 1, 1, 1, 1, 1, 1); #ImagePlus imp,   int firstC, int lastC,    int firstZ, int lastZ,   int firstT, int lastT)
	results = ArrayList()
	im = IntensityMeasures( inputimg, labelimg )
	results.add( im.getMean() )
	results.add( im.getStdDev() )
	results.add( im.getNumberOfVoxels())
	results.add( im.getMin() )
	results.add( im.getMax() )
	results.add( im.getMedian()  )
	results.add( im.getMode()  )
	
	mergedTable = ResultsTable()
	numLabels = results.get(0).getCounter()

	###Create a dictionary to store data###
	d={}
	d["label"]=[]
	for i in xrange(results.size()): #for each heading (mean, std. dev. etc.)
		measure = results.get( i ).getColumnHeading( 0 )
		d[measure]=[]
	######################################
	
	for i in xrange(numLabels):
		mergedTable.incrementCounter()
		label = results.get( 0 ).getLabel( i ) #obtains the 0-indexed ith label, regardless of its string-name.
		d["label"].append(label)
		mergedTable.addLabel(label)
		for j in xrange(results.size()):
			measure = results.get( j ).getColumnHeading( 0 )
			value = results.get( j ).getValue( measure, i )
			mergedTable.addValue( measure, value )
			d[measure].append(value)

	if show_table:
		mergedTable.show( inputimg.getShortTitle() +"-intensity-measurements" )
	


	###Ensure labels file is in the correct format: ###

    #Labels sometimes have gaps (e.g. labels=[4,40,82]  is possible). 
    #The Python script stores them in a python list, and accesses them by “python indexes” (i.e. their order, starting with 0)  
    #In this example, label 4 would have a python index of 0 and label 40 would have a python index of 1 etc.
  
	tmp=map(int, d["label"]) #convert label numbers (strings) to integers
	assert sorted(tmp) == tmp, "FATAL ERROR: The labels provided are not in numerical order, \
								whereas this script was written assuming they are. \
								If this error occurs, it means the script needs editing"
    ###################################################
    
	if manually_assign_backgroundlayer_to_label:
		background_label_index=tmp.index(manually_assign_backgroundlayer_to_label)
		print("The background has been manually selected as label {} (i.e. python index {})".format(manually_assign_backgroundlayer_to_label, background_label_index)) 
	else:
		background_label_index, background_number_of_voxels = max(enumerate(d["NumberOfVoxels"]), key=operator.itemgetter(1))
		print("The auto-selected background is at label {} (i.e. python index {})".format(d["label"][background_label_index], background_label_index)) 

	return d["label"], numLabels, background_label_index, d["label"][background_label_index]

def build_database(numLabels, rawimg, labelimg):
	"""Return a database 'df' which contains all the mean intensity data for each cell 
	across both channels and all z-positions. 

	df is a Python dictionary which has two keys - df['channel1'] and df['channel2']
	Each channel's data is stored as another dictionary with each key being the label of the cell.
	e.g. df['channel1'][0] will provide the mean intensity of the first cell in channel 1. 
	Note that the label name is unknown but can be found with function: LabelImages.findAllLabels(labelimg)[0] where 0 is the index] 
	e.g. df['channel2'][background_label_index] will provide the mean intensity of the background in channel 2
    """
	
	results = ArrayList()
	nSlices=rawimg.getNSlices()
	print("Creating a database of over both channels and over {} z-slices...".format(nSlices))

	nFrames=rawimg.getNFrames()
	#Create a dictionary with key channel1 whose value is a dictionary of key time_index (0,1,2...nFrames) whose value if an empty dictionary
	df={}
	df["channel1"]=dict([i for i in zip(xrange(nFrames),   [{} for i in xrange(nFrames)]    )]) 
	df["channel2"]=dict([i for i in zip(xrange(nFrames),   [{} for i in xrange(nFrames)]    )])
	for t_index in xrange(nFrames):
		for label_index in xrange(numLabels):
			df["channel1"][t_index][label_index]=[] #The list will contain a list average intensity for each cell in that channel at time=t_index
			df["channel2"][t_index][label_index]=[] #The list will contain a list average intensity for each cell in that channel at time=t_index

	for t_index in xrange(nFrames):
		for c in [1,2]:#iterate over channels
			for z in xrange(nSlices):#iterate over z-dimension
				inputimg=Duplicator().run(rawimg, c, c, z, z, t_index+1, t_index+1)
				im = IntensityMeasures( inputimg, labelimg )
				for label_index in xrange(numLabels): #Add mean intensity for each cell
					if c == 1:
						if for_each_cell_calculate == "Mean":
							df["channel1"][t_index][label_index].append(im.getMean().getValue("Mean",label_index)) #Store mean intensity for the label at label_index at a z position of z 
						elif for_each_cell_calculate == "Median":
							df["channel1"][t_index][label_index].append(im.getMedian().getValue("Median",label_index)) #Store median intensity for the label at label_index at a z position of z 
					elif c == 2:
						if for_each_cell_calculate == "Mean":
							df["channel2"][t_index][label_index].append(im.getMean().getValue("Mean",label_index)) #Store mean intensity for the label at label_index at a z position of z 
						elif for_each_cell_calculate == "Median":
							df["channel2"][t_index][label_index].append(im.getMedian().getValue("Median",label_index)) #Store median intensity for the label at label_index at a z position of z 

	return df

def process_jjmfile(csvwriter, manual_data_name=None, manual_label_name=None, hardcoded_data_name=None):
	
	rawimg, labelimg, rawimg_name, labelimg_name=load_environment(manual_data_name, manual_label_name)
	
	print("Loading first frame to show a tabular summary (and to find the background)...") 
	labels, numLabels, background_label_index, background_label = summarise_input(rawimg, labelimg) 
	
	df=build_database(numLabels, rawimg, labelimg)
	
	print("Database Built. Now calculating the intensity...")
	for t_index in xrange(rawimg.getNFrames()):
		if rawimg.getNFrames() > 1:
			print("====Analysing Time Frame {} ====".format(t_index+1))
			
		r_tot=0
		per_cell_intensity_list=[]
		cell_label_index_list=[i for i in xrange(numLabels)]
		cell_label_index_list.remove(background_label_index)
	
		for cell_label_index in cell_label_index_list: 
			zpos_index, maxintensity=max(enumerate(df["channel1"][t_index][cell_label_index]), key=operator.itemgetter(1))
			r_cell=(df["channel1"][t_index][cell_label_index][zpos_index] - df["channel1"][t_index][background_label_index][zpos_index])/(df["channel2"][t_index][cell_label_index][zpos_index] - df["channel2"][t_index][background_label_index][zpos_index])
			per_cell_intensity_list.append(r_cell)
	
		N=len(per_cell_intensity_list)
		r_mean=sum(per_cell_intensity_list)/float(N)
		r_median=median(per_cell_intensity_list)
		r_std=standard_deviation(per_cell_intensity_list, population=False)
		
		print("Analysis Complete!")
		print("Mean:{}\nStd:{}\nN:{}\nMedian:{}".format(r_mean, r_std, N, r_median))
		
		print("Each cell's average intensity is shown below:")
		print(per_cell_intensity_list)
		nwarnings=0
		for i, cell in enumerate(per_cell_intensity_list):
			if cell < warning_threshold:
				print("WARNING: Cell Label {} (i.e. python index {}) has an low intensity of: {}".format(labels[i], i, cell)) 
				nwarnings+=1

		if hardcoded_data_name:
			filenamestring=hardcoded_data_name
		else:
			filenamestring=t_index
			
		csvwriter.writerow([filenamestring,r_mean,r_std,r_median,N,nwarnings,per_cell_intensity_list])
		print(" ")
		
	return r_mean,r_std,r_median,N,nwarnings, per_cell_intensity_list




#==================
######
#This block of code can be used instead of process_jjmfile() to process a series of .lsm images
#And output the results to a .csv file in the directory of this script
######

path='C:\\CellAnalysis\\Analysis1\\' #path at which more than one .lsm files can be found
f=open(os.path.join(os.path.dirname(os.path.realpath(sys.argv[0])), 'summary.csv'), 'wb') 
csvwriter = csv.writer(f)
csvwriter.writerow(["filename","mean","std","median","N","N_warnings","cell_intensities"])
number_of_files=15
for i in xrange(number_of_files):
	hardcoded_data_name="HeLa30uMtimelapse"+str(i+1)+".lsm"
	process_jjmfile(csvwriter, os.path.join(path, hardcoded_data_name), os.path.join(path, "HeLa30uMtimelapse1_label_wc.tif"), hardcoded_data_name)
f.close()
print("Script Complete.")
#==================


"""

#==================
######
#This block of code can be used to process a single .lsm file (with 1 or many time steps)
#And output the results to a .csv file in the directory of this script
######
#For single .lsm file 
path=os.path.dirname(os.path.realpath(sys.argv[0])) #open a file in path of this script
f=open(os.path.join(path, 'summary.csv'), 'wb') 
csvwriter = csv.writer(f)
csvwriter.writerow(["time-step","mean","std","median","N","N_warnings","cell_intensities"])
process_jjmfile(csvwriter)
f.close()
print("Script Complete.")
#==================
"""
