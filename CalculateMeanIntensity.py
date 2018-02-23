from ij import IJ, WindowManager, ImagePlus
from java.util import ArrayList
from ij.gui import GenericDialog
from ij.measure import ResultsTable
from ij.plugin import PlugIn, Duplicator
from inra.ijpb.measure import IntensityMeasures
from inra.ijpb.label.LabelImages import keepLargestLabel
from inra.ijpb.label import LabelImages
import operator

######User Parameters########
window_number_of_the_data=1   #1 = first
window_number_of_the_labels=2 #2 = second
#By default background is set to largest label, but can be manually set here:
manually_assign_backgroundlayer_to_label=None #Based on imageJ label number (e.g. Image>Color>Color Picker's value)
#############################


nbima = WindowManager.getImageCount()
print("There are this many windows open: "+str(nbima))
rawimg=WindowManager.getImage(window_number_of_the_data)
labelimg=WindowManager.getImage(window_number_of_the_labels)
rawimg_name=rawimg.getShortTitle()
labelimg_name=labelimg.getShortTitle()
print("Data filename is: "+str(rawimg_name))
print("Labels filename is: "+str(labelimg_name))
labels = LabelImages.findAllLabels( labelimg )
print("Which has these labels...")
print(labels)

#largestlabelimg=keepLargestLabel(labelimg) # can be used to make a label of just the background
if rawimg.getWidth() != labelimg.getWidth() or rawimg.getHeight() != labelimg.getHeight():
	IJ.error( "Intensity Measures 2D/3D input error", "Error: input"
			+ " and label images must have the same size" );


def summarise_input(rawimg, labelimg):
	"""This function takes as input: 'rawimg' (the data) and 'labelimg' (the cell boundary cartoon) 
	Then using the z=1, channel=1 frame, produces a summary table for inspection. 
	It also calculates which label is the background, assuming it is the largest cell.
	It returns the following:
	1. The number of labels
	2. The index of the background (i.e. the label number of the background minus 1)
	"""
	#Take a snapshot of the frame at z=1 and c=1:
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
		
	####Find Background Label and Summarise Data for First C=1, Z=1 slide in a table.
	
	mergedTable = ResultsTable()
	numLabels = results.get(0).getCounter()
	d={}
	d["label"]=[]
	for i in xrange(results.size()):
		measure = results.get( i ).getColumnHeading( 0 )
		d[measure]=[]
		
	for i in xrange(numLabels):
		mergedTable.incrementCounter()
		label = results.get( 0 ).getLabel( i )
		d["label"].append(label)
		mergedTable.addLabel(label)
		for j in xrange(results.size()):
			measure = results.get( j ).getColumnHeading( 0 )
			value = results.get( j ).getValue( measure, i )
			mergedTable.addValue( measure, value )
			d[measure].append(value)
	mergedTable.show( inputimg.getShortTitle() +"-intensity-measurements" )
	
	background_label_index, background_number_of_voxels = max(enumerate(d["NumberOfVoxels"]), key=operator.itemgetter(1))
	print("The auto-selected background is at label {} (i.e. python index {})".format(d["label"][background_label_index], background_label_index)) 


	###############Ensure labels file is in the correct format: ########################
	tmp=map(int, d["label"]) #convert label numbers (strings) to integers
	assert sorted(tmp) == tmp, "FATAL ERROR: The labels provided are not in numerical order, \
								whereas this script was written assuming they are. \
								If this error occurs, it means the script needs editing"
    ####################################################################################
	if manually_assign_backgroundlayer_to_label:
		background_label_index=manually_assign_backgroundlayer_to_label-1

		
	
	return numLabels, background_label_index

def build_database(numLabels, rawimg, labelimg):
	"""Return a database 'df' which contains all the mean intensity data for each cell 
	across both channels and all z-positions. 

	df is a Python dictionary which has two keys - df['channel1'] and df['channel2']
	Each channel's data is stored as another dictionary with each key being the label of the cell.
	e.g. df['channel1'][0] will provide the mean intensity of the first cell (label=1) in channel 1 
	e.g. df['channel2'][background_label_index] will provide the mean intensity of the background in channel 2
    """
	
	results = ArrayList()
	nSlices=rawimg.getNSlices()
	print("Creating a database of over both channels and over {} z-slices...".format(nSlices))
	
	df={"channel1":{}, "channel2":{}}
	
	for label in xrange(numLabels):
		df["channel1"][label]=[] #assumes labels are 0,1,2... etc. with no gaps. The list will be a list mean intensity for each cell in that channel
		df["channel2"][label]=[] #assumes labels are 0,1,2... etc. with no gaps. The list will be a list mean intensity for each cell in that channel
		
	for c in [1,2]:#iterate over channels
		for z in xrange(nSlices):#iterate over z-dimension
			inputimg=Duplicator().run(rawimg, c, c, z, z, 1, 1)
			im = IntensityMeasures( inputimg, labelimg )
			for label in xrange(numLabels): #Add mean intensity for each cell
				if c == 1:
					df["channel1"][label].append(im.getMean().getValue("Mean",label)) #mean for the label called label at z slice called z
				elif c == 2:
					df["channel2"][label].append(im.getMean().getValue("Mean",label)) #mean for the label called label at z slice called z

	return df

print("Loading first frame to show a tabular summary (and to find the background)...") 
numLabels, background_label_index = summarise_input(rawimg, labelimg) 
df=build_database(numLabels, rawimg, labelimg)
print("Database Built. Now calculating the intensity...")

r_tot=0
cell_label_list=[i for i in xrange(numLabels)]
cell_label_list.remove(background_label_index)

for cell_label in cell_label_list: 
	zpos_index, maxintensity=max(enumerate(df["channel1"][cell_label]), key=operator.itemgetter(1))
	r_cell=(df["channel1"][cell_label][zpos_index] - df["channel1"][background_label_index][zpos_index])/(df["channel2"][cell_label][zpos_index] - df["channel2"][background_label_index][zpos_index])
	r_tot+=r_cell

r_mean=r_tot/float(len(cell_label_list))

print("Analysis Complete! Result={}.".format(r_mean))

