#
# reads Bruker MR data (Paravision v5.1)
# vizualizes raw acquisition data (FID files)
# this version is for the FLOWMAP method only
#
# ----- VERSION HISTORY -----
#
# Version 0.1 - 08, August 2017
#       - 1st public github Release
#
# ----- REQUIREMENTS ----- 
#
#    This program was developed under Python Version 2.7
#    with the following additional libraries: 
#    - numpy
#    - matplotlib.pylab
#

from __future__ import print_function
import sys
import os
import numpy as np
import pylab as pl


TK_installed=True
try: from tkFileDialog import askopenfilename # Python 2
except: 
	try: from tkinter.filedialog import askopenfilename; # Python3
	except: TK_installed=False
try: import Tkinter as tk; # Python2
except: 
	try: import tkinter as tk; # Python3
	except: TK_installed=False
if not TK_installed:
	print ('ERROR: tkinter not installed')
	print ('       on Linux try "yum install tkinter"')
	print ('       on MacOS install ActiveTcl from:')
	print ('       http://www.activestate.com/activetcl/downloads')
	sys.exit(2)

def ReadParamFile(filepath):
    global OrigFilename;
    #Read a Bruker MRI experiment's method or acqp file to a dictionary.
    param_dict = {}
    with open(filepath, "r") as f:
        while True:
            line = f.readline()
            if not line:
                break
            # when line contains parameter
            if line.startswith('$$ /'):
			    OrigFilename = line[line.find('/nmr/')+5:]
			    OrigFilename = OrigFilename[0:len(OrigFilename)-8]
			    OrigFilename = OrigFilename.replace(".", "_")
			    OrigFilename = OrigFilename.replace("/", "_")
            if line.startswith('##$'):
                (param_name, current_line) = line[3:].split('=') # split at "="
                # if current entry (current_line) is arraysize
                if current_line[0:2] == "( " and current_line[-3:-1] == " )":
                    value = ParseArray(f, current_line)
                # if current entry (current_line) is struct/list
                elif current_line[0] == "(" and current_line[-3:-1] != " )":
                    # if neccessary read in multiple lines
                    while current_line[-2] != ")":
                        current_line = current_line[0:-1] + f.readline()
                    # parse the values to a list
                    value = [ParseSingleValue(x)
                             for x in current_line[1:-2].split(', ')]
                # otherwise current entry must be single string or number
                else:
                    value = ParseSingleValue(current_line)
                # save parsed value to dict
                param_dict[param_name] = value
    return param_dict
	
def ParseArray(current_file, line):
    # extract the arraysize and convert it to numpy
    line = line[1:-2].replace(" ", "").split(",")
    arraysize = np.array([int(x) for x in line])
    # then extract the next line
    vallist = current_file.readline().split()
    # if the line was a string, then return it directly
    try:
        float(vallist[0])
    except ValueError:
        return " ".join(vallist)
    # include potentially multiple lines
    while len(vallist) != np.prod(arraysize):
        vallist = vallist + current_file.readline().split()
    # try converting to int, if error, then to float
    try:
        vallist = [int(x) for x in vallist]
    except ValueError:
        vallist = [float(x) for x in vallist]
    # convert to numpy array
    if len(vallist) > 1:
        return np.reshape(np.array(vallist), arraysize)
    # or to plain number
    else:
        return vallist[0]

def ParseSingleValue(val):
    try: # check if int
        result = int(val)
    except ValueError:
        try: # then check if float
            result = float(val)
        except ValueError:
            # if not, should  be string. Remove  newline character.
            result = val.rstrip('\n')
    return result	
	
	
#general initialization stuff  
space=' '; slash='/'; 
if sys.platform=="win32": slash='\\' # not really needed, but looks nicer ;)
Program_name = os.path.basename(sys.argv[0]); 
if Program_name.find('.')>0: Program_name = Program_name[:Program_name.find('.')]
python_version = str(sys.version_info[0])+'.'+str(sys.version_info[1])+'.'+str(sys.version_info[2])
# sys.platform = [linux2, win32, cygwin, darwin, os2, os2emx, riscos, atheos, freebsd7, freebsd8]
if sys.platform=="win32": os.system("title "+Program_name)
	
#TK initialization       
TKwindows = tk.Tk(); TKwindows.withdraw() #hiding tkinter window
TKwindows.update()
# the following tries to disable showing hidden files/folders under linux
try: TKwindows.tk.call('tk_getOpenFile', '-foobarz')
except: pass
try: TKwindows.tk.call('namespace', 'import', '::tk::dialog::file::')
except: pass
try: TKwindows.tk.call('set', '::tk::dialog::file::showHiddenBtn', '1')
except: pass
try: TKwindows.tk.call('set', '::tk::dialog::file::showHiddenVar', '0')
except: pass
TKwindows.update()
	
#intercatively choose input FID file
FIDfile = askopenfilename(title="Choose Bruker FID file", filetypes=[("FID files","fid")])
if FIDfile == "": print ('ERROR: No FID input file specified'); sys.exit(2)
FIDfile = os.path.abspath(FIDfile) 
TKwindows.update()

#read FID 
with open(FIDfile, "r") as f: FIDrawdata= np.fromfile(f, dtype=np.int32) 
FIDrawdata_CPX = FIDrawdata[0::2] + 1j * FIDrawdata[1::2]
FIDrawdata = 0 #free memory

#read method file
METHODfile=os.path.dirname(FIDfile)+slash+'method'
METHODdata=ReadParamFile(METHODfile)

#check for not implemented stuff
if METHODdata["Method"] != "FLOWMAP" or METHODdata["FlowMode"] != "VelocityMapping":
	print ('ERROR: Recon only implemented for FLOWMAP VelocityMapping'); sys.exit(1)
if METHODdata["PVM_SpatDimEnum"] != "3D" or METHODdata["FlowEncodingDirection"] != "AllDirections":
	print ('ERROR: Recon only implemented for 3D acquisition with flow incoding in all directions'); sys.exit(1)
if METHODdata["FlowEncLoop"] !=4:
	print ('ERROR: ops, expected flow encoding loop = 4 '); sys.exit(1)
if METHODdata["PVM_NSPacks"] != 1:
	print ('ERROR: Recon only implemented 1 package'); sys.exit(1)
if METHODdata["PVM_NRepetitions"] != 1:
	print ('ERROR: Recon only implemented 1 repetition'); sys.exit(1)
if METHODdata["PVM_EncPpiAccel1"] != 1 or METHODdata["PVM_EncPftAccel1"] != 1 or \
   METHODdata["PVM_EncZfAccel1"] != 1 or METHODdata["PVM_EncZfAccel2"] != 1 or \
   METHODdata["PVM_EncTotalAccel"] != 1 or METHODdata["PVM_EncNReceivers"] != 1:
	print ('ERROR: Recon for parallel acquisition not implemented'); sys.exit(1)

#reshape FID data according to dimensions from method file
#"order="F" means Fortran style order as by BRUKER conventions
dim=METHODdata["PVM_EncMatrix"]
try: FIDrawdata_CPX = FIDrawdata_CPX.reshape(dim[0],4,dim[1],dim[2], order="F")
except: print ('ERROR: k-space data reshape failed (dimension problem)'); sys.exit(1)
print (FIDrawdata_CPX.shape)
	
#view k-space
EchoPosition_raw=METHODdata["PVM_EchoPosition"]
EchoPosition_raw=50-(50-EchoPosition_raw)
EchoPosition=int(EchoPosition_raw/100.*dim[0])
# show the resulting images
print(dim)
pl.figure()
pl.subplot(4,3,1)
pl.imshow(abs(FIDrawdata_CPX[:,0,:,dim[2]/2]))
pl.subplot(4,3,2)
pl.imshow(abs(FIDrawdata_CPX[:,0,dim[1]/2,:]))
pl.subplot(4,3,3)
pl.imshow(abs(FIDrawdata_CPX[EchoPosition,0,:,:]))
pl.subplot(4,3,4)
pl.imshow(abs(FIDrawdata_CPX[:,1,:,dim[2]/2]))
pl.subplot(4,3,5)
pl.imshow(abs(FIDrawdata_CPX[:,1,dim[1]/2,:]))
pl.subplot(4,3,6)
pl.imshow(abs(FIDrawdata_CPX[EchoPosition,1,:,:]))
pl.subplot(4,3,7)
pl.imshow(abs(FIDrawdata_CPX[:,2,:,dim[2]/2]))
pl.subplot(4,3,8)
pl.imshow(abs(FIDrawdata_CPX[:,2,dim[1]/2,:]))
pl.subplot(4,3,9)
pl.imshow(abs(FIDrawdata_CPX[EchoPosition,2,:,:]))
pl.subplot(4,3,10)
pl.imshow(abs(FIDrawdata_CPX[:,3,:,dim[2]/2]))
pl.subplot(4,3,11)
pl.imshow(abs(FIDrawdata_CPX[:,3,dim[1]/2,:]))
pl.subplot(4,3,12)
pl.imshow(abs(FIDrawdata_CPX[EchoPosition,3,:,:]))
pl.show()
	
	