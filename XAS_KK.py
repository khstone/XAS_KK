from numpy import *
from scipy import *
import pylab

#Create an array to link atom name with Z
periodic_table = ['h', 'he', 'li', 'be', 'b', 'c', 'n', 'o', 'f', 'ne', 'na', 
'mg', 'al', 'si', 'p', 's', 'cl', 'ar', 'k', 'ca', 'sc', 'ti', 'v', 'cr', 'mn', 
'fe', 'co', 'ni', 'cu', 'zn', 'ga', 'ge', 'as', 'se', 'br', 'kr', 'rb', 'sr', 
'y', 'zr', 'nb', 'mo', 'tc', 'ru', 'rh', 'pd', 'ag', 'cd', 'in', 'sn', 'sb', 
'te', 'i', 'xe', 'cs', 'ba', 'la', 'ce', 'pr', 'nd', 'pm', 'sm', 'eu', 'gd', 
'tb', 'dy', 'ho', 'er', 'tm', 'yb', 'lu', 'hf', 'ta', 'w', 're', 'os', 'ir', 
'pt', 'au', 'hg', 'tl', 'pb', 'bi', 'po', 'at', 'rn', 'fr', 'ra', 'ac', 'th', 
'pa', 'u']

# Change this to reflect where the scattering factor files are on your computer
atompath = "C:/Users/Kevin/Dropbox/Projects/Python Stuff/ScatteringFactors/"

#Display introductory information
print "Welcome to XAS Normalizer"
print "This program will normalize XAS data"
print "Please specify a scan file to work with (include path if not in this directory):"
raw_name = raw_input()

#open specified data file
raw_file=open(raw_name, 'r')
raw_data=[]

#read in the first 12 lines (header information)
#for i in range (0,0):
#	line = raw_file.readline()

#read lines until we find one that is not empty
line=raw_file.readline()
while line.isspace():
	line=raw_file.readline()
	
#assign counter names to a variable
counter_names = line.split('\t')

#read in the first two lines of data, I need two lines to make a 2-dimensional array
line=raw_file.readline()
temp=line.split()
l1=[]
for a in temp:
	b=float(a)
	l1.append(b)
line=raw_file.readline()
temp=line.split()
l2=[]
for a in temp:
	b=float(a)
	l2.append(b)
	
#write the first two lines to the array data
data=array([l1,l2])

#loop through the rest of the data lines adding them to the array
while line:
	line=raw_file.readline()
	temp=line.split()
	l1=[]
	for a in temp:
		b= float(a)
		l1.append(b)
	if shape(l1)==shape(data[1]):   #to make sure we dont try to append the empty last line
		data=append(data,[l1],axis=0)
raw_file.close()
#data is now an array of the data entries, i.e. data[point number][counter]

#Open file and extract f2 and f1 data from Henke
print "Please enter the element abbreviation you want to work with (in lowercase only):"
atomtype = raw_input()
Z = periodic_table.index(atomtype) + 1  #since the indexing starts at 0
atomname = atompath + atomtype + ".nff"
atom_file = open(atomname, 'r')
line = atom_file.readline()
#read in the first two lines of data, I need two lines to make a 2-dimensional array
line=atom_file.readline()
temp=line.split()
l1=[]
for a in temp:
	b=float(a)
	l1.append(b)
line=atom_file.readline()
temp=line.split()
l2=[]
for a in temp:
	b=float(a)
	l2.append(b)
	
#write the first two lines to the array data
atom_data=array([l1,l2])

#loop through the rest of the data lines adding them to the array
while line:
	line=atom_file.readline()
	temp=line.split()
	l1=[]
	for a in temp:
		b= float(a)
		l1.append(b)
	if shape(l1)==shape(atom_data[1]):   #to make sure we dont try to append the empty last line
		atom_data=append(atom_data,[l1],axis=0)
atom_file.close()
#atom_data is now an array of the data entries, i.e. data[point number][counter]

#Make a plot of the specified counters
for i in range (0, size(counter_names)):
	print '[', i, '] - ', counter_names[i]
print "Enter the counter number to use for x:"
x_num = int(raw_input())
print "Enter the counter number to use for y:"
y_num = int(raw_input())
print "Enter the counter number to use for Io:"
I0 = int(raw_input())

#Pick the low energy linear fitting range by clicking on the plot itself
print "Click on the points to use for linear fitting"
pylab.plot(data[:,x_num],1*data[:,y_num]/data[:,I0])
linpoints = pylab.ginput(2)
pylab.close()
low_E = linpoints[0][0]
high_E = linpoints[1][0]
if low_E > high_E:
	low_E,high_E=high_E,low_E
num_points = size(data)/size(counter_names)
lin_x=[]
lin_y=[]
for i in range (0, num_points):
	test_x = data[i][x_num]
	if low_E <= test_x <= high_E:
		lin_x = append(lin_x, test_x)
		lin_y = append(lin_y, 1*data[i][y_num]/data[i][I0])
line = polyfit(lin_x, lin_y, 1)
#subtract the linear extrapolation from the data set
new_y = array(1*data[:,y_num]/data[:,I0] - line[0]*data[:,x_num] - line[1])

#Prompt for point to use to normalize to 1
print "Pick the energy to normalize to 1:"
pylab.plot(data[:,x_num], new_y[:])
set1 = pylab.ginput(1)
pylab.close()
for i in range (0, num_points):
	if data[i,x_num] >= set1[0][0]:
		break
multfactor = 1/new_y[i]
new_y = new_y*multfactor

#Find the edge using the TEY data range
min_E = data[:,x_num].min()
max_E = data[:,x_num].max()
end_pre_edge = high_E
start_pre_edge = min_E
start_post_edge = set1[0][0]
end_post_edge = max_E
atom_points = size(atom_data)/3

#Find the fit to the pre and post edge regions
pre_x=[]
pre_y=[]
post_x=[]
post_y=[]
f2_x=[]
f2_y=[]
for i in range (0, atom_points):
	test_x = atom_data[i,0]
	if start_pre_edge <= test_x <= end_pre_edge:
		pre_x = append(pre_x, test_x)
		pre_y = append(pre_y, atom_data[i,2])
	if start_post_edge - 25 <= test_x <= end_post_edge + 40:
		post_x = append(post_x, test_x)
		post_y = append(post_y, atom_data[i,2])
	if start_pre_edge <= test_x <= end_post_edge:
		f2_x = append(f2_x, test_x)
		f2_y = append(f2_y, atom_data[i,2])
log_pre_x = log(pre_x)
log_pre_y = log(pre_y)
log_pre_edge = polyfit(log_pre_x, log_pre_y, 1)
log_post_x = log(post_x)
log_post_y = log(post_y)
log_post_edge = polyfit(log_post_x, log_post_y, 3)

#Plot the edge region of the Henke data and the linear fits to the pre and post edges
pre_edge = exp(log_pre_edge[0]*log(f2_x[:]) + log_pre_edge[1])
post_edge =  exp(log(f2_x[:])*log_post_edge[0] + log_post_edge[1])
pylab.plot(f2_x[:], f2_y[:])
pylab.plot(f2_x[:], pre_edge[:])
pylab.plot(f2_x[:], post_edge[:])

#Subtract the pre_edge extrapolation from the post_edge extrapolation
pre_edge = exp(log_pre_edge[0]*log(data[:,x_num]) + log_pre_edge[1])
post_edge =  exp(log(data[:,x_num])*log_post_edge[0] + log_post_edge[1])
post_edge = post_edge - pre_edge

#multiply data by the subtracted post-edge extrapolation
new_y = new_y * post_edge

#add back in the pre-edge extrapolation
new_y = new_y + pre_edge
post_edge = post_edge + pre_edge
pylab.figure()
pylab.plot(f2_x[:], f2_y[:], 'k-')
pylab.plot(data[:,x_num], pre_edge[:], 'r.')
pylab.plot(data[:,x_num], post_edge[:], 'r.')
pylab.plot(data[:,x_num], new_y[:])
pylab.xlabel('Energy (eV)')
pylab.ylabel('f2')

#write an xye formatted file with the f2 data
xye_name = raw_name[:-4]
xye_name += '_f2.xye'
xye_file=open(xye_name, 'w')
x = array(data[:,x_num])
y = array(new_y)
xye = (x,y)
for i in range (0, num_points):
	x_string = str(xye[0][i])
	y_string = str(xye[1][i])
	xye_file.write(x_string)
	xye_file.write('\t')
	xye_file.write(y_string)
	xye_file.write('\n')
xye_file.close()

f1_y = []
f1_x = []

#create full f2_y which will be the measured data stitched onto the full henke data set
full_f2_x = []
full_f2_y = []
for i in range (0, atom_points):
	if atom_data[i,0] >= min_E:
		start_i = i
		break
	full_f2_x = append(full_f2_x, atom_data[i,0])
	full_f2_y = append(full_f2_y, atom_data[i,2])
for i in range (0, num_points):
	full_f2_x = append(full_f2_x, x[i])
	full_f2_y = append(full_f2_y, y[i])
	if x[i] >= set1[0][0]:
		extra = i
		break
for i in range (start_i, atom_points):
	if atom_data[i,0] <= max(x):
		stop_i = i + extra
	if atom_data[i,0] >= set1[0][0]:
		full_f2_x = append(full_f2_x, atom_data[i,0])
		full_f2_y = append(full_f2_y, atom_data[i,2])

Zstar = Z - (Z/82.5)**2.37

#Kramers Kronig for the full_f2 data
full_f1_y = []
full_points = size(full_f2_x)
for i in range(0, full_points):
	f1_temp = 0
	for j in range (0, full_points):
		if j == 0:
			dE = 0.5*(full_f2_x[j+1] - full_f2_x[j])
		if j == full_points:
				dE = 0.5*(full_f2_x[j] - full_f2_x[j-1])
		if 0 < j < (full_points-1):  #not sure why full_points - 1, but just atom_points doesnt seem to work
				dE = 0.5*(full_f2_x[j+1] - full_f2_x[j-1])
		if full_f2_x[j] != full_f2_x[i]:
			f1_temp += dE*(-2/ pi )*((full_f2_x[j]*full_f2_y[j] - full_f2_x[i]*full_f2_y[i])/(full_f2_x[j]**2 - full_f2_x[i]**2))
		if full_f2_x[j] == full_f2_x[i]:
			f1_temp += dE*(-2/pi)*(full_f2_y[i]/(2*full_f2_x[i]))
	full_f1_y = append(full_f1_y, f1_temp + Zstar)

#Hilbert transform to do the KK transform
#hilbert_f2_y = -1*(hilbert(full_f2_y))

#plot f2, and both version of f1
pylab.figure()
pylab.plot(full_f2_x[:], full_f1_y[:], 'r-')
pylab.plot(full_f2_x, full_f2_y[:], 'b-')
#pylab.plot(full_f2_x, hilbert_f2_y, 'g-')
#pylab.axis([680, 760, -40, 60])

#spit out a file with Energy, F1, F2 for the edge
xff_name = raw_name[:-4]
xff_name += '.xff'
xff_file=open(xff_name, 'w')
#x = array(data[:,x_num])
for i in range (start_i, stop_i):
	x_string = str(full_f2_x[i])
	f1_string = str(full_f1_y[i])
	f2_string = str(full_f2_y[i])
	xff_file.write(x_string)
	xff_file.write('\t')
	xff_file.write(f1_string)
	xff_file.write('\t')
	xff_file.write(f2_string)
	xff_file.write('\n')
xff_file.close()

pylab.show()

