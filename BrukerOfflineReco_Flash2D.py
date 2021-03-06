#
# reads Bruker MR data (Paravision v5.1)
# reconstructs images from raw acquisition data (FID files)
# this version is for the FLASH 3D method only
#
# ----- VERSION HISTORY -----
#
# Version 0.1 - 08, August 2017
#       - 1st public github Release
#
# ----- LICENSE -----                 
#
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#    For more detail see the GNU General Public License.
#    <http://www.gnu.org/licenses/>.
#
#    THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
#    IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
#    FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
#    AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
#    LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
#    OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
#    THE SOFTWARE.
#
# ----- REQUIREMENTS ----- 
#
#    This program was developed under Python Version 2.7
#    with the following additional libraries: 
#    - numpy
#    - nibabel
#    - matplotlib.pylab (optional)
#

from __future__ import print_function
try: import win32gui, win32console
except: pass #silent
from math import ceil, floor
import sys
import os
import numpy as np
import nibabel as nib


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
                (param_name, current_line) = line[3:].split('=') #split at "="
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

def smooth(x,window_len):
    w=np.hanning(window_len)
    s=np.r_[2*x[0]-x[window_len-1::-1],x,2*x[-1]-x[-1:-window_len:-1]]
    w=np.hanning(window_len)
    y=np.convolve(w/w.sum(),s,mode='same')
    return y[window_len:-window_len+1]  
    
def FFT2D (array):
    for k in range(0,array.shape[1]): array[:,k,:,:] = np.fft.fft(array[:,k,:,:], axis=(0))
    for k in range(0,array.shape[1]): array[:,k,:,:] = np.fft.fft(array[:,k,:,:], axis=(1))        
    return array 

def iFFT2D (array):
    for k in range(0,array.shape[1]): array[:,k,:,:] = np.fft.ifft(array[:,k,:,:], axis=(0))
    for k in range(0,array.shape[1]): array[:,k,:,:] = np.fft.ifft(array[:,k,:,:], axis=(1))   
    return array 
    
def RG_to_voltage(RG):
    return np.power(10,11.995/20.) * np.power(RG,19.936/20.)
    # this comes from the Bruker provided conversion list below
    #
    #
    # Receiver,   Gain Equivalent Voltage Gain [dB],        V_out/V_in
    #             20*log (V_out/V_in)
    #
    # 2050,       78,	                                    7943.282347	
    # 1820,       77,	                                    7079.457844	
    # 1620,       76,	                                    6309.573445	
    # 1440,       75,	                                    5623.413252	
    # 1290,       74,	                                    5011.872336	
    # 1150,       73,	                                    4466.835922	
    # 1030,       72,	                                    3981.071706	
    # 912,        71,	                                    3548.133892	
    # 812,        70,	                                    3162.27766	
    # 724,        69,	                                    2818.382931	
    # 645,        68,	                                    2511.886432	
    # 575,        67,	                                    2238.721139	
    # 512,        66,	                                    1995.262315	
    # 456,        65,	                                    1778.27941	
    # 406,        64,	                                    1584.893192	
    # 362,        63,	                                    1412.537545	
    # 322,        62,	                                    1258.925412	
    # 287,        61,	                                    1122.018454	
    # 256,        60,	                                    1000	    
    # 228,        59,	                                    891.2509381	
    # 203,        58,	                                    794.3282347	
    # 181,        57,	                                    707.9457844	
    # 161,        56,	                                    630.9573445	
    # 144,        55,	                                    562.3413252	
    # 128,        54,	                                    501.1872336	
    # 114,        53,	                                    446.6835922	
    # 101,        52,	                                    398.1071706	
    # 90.5,       51,	                                    354.8133892	
    # 80.6,       50,	                                    316.227766	
    # 71.8,       49,	                                    281.8382931	
    # 64,         48,	                                    251.1886432	
    # 57,         47,	                                    223.8721139	
    # 50.8,       46,	                                    199.5262315	
    # 45.2,       45,	                                    177.827941	
    # 40.3,       44,	                                    158.4893192	
    # 36,         43,	                                    141.2537545	
    # 32,         42,	                                    125.8925412	
    # 28.5,       41,	                                    112.2018454	
    # 25.4,       40,	                                    100	        
    # 22.6,       39,	                                    89.12509381	
    # 20.2,       38,	                                    79.43282347	
    # 18,         37,	                                    70.79457844	
    # 16,         36,	                                    63.09573445	
    # 14.2,       35,	                                    56.23413252	
    # 12.7,       34,	                                    50.11872336	
    # 11.3,       33,	                                    44.66835922	
    # 10,         32,	                                    39.81071706	
    # 9,    	  31,	                                    35.48133892	
    # 8,          30,	                                    31.6227766	
    # 7.12,       29,	                                    28.18382931	
    # 6.35,       28,	                                    25.11886432	
    # 5.6,        27,	                                    22.38721139	
    # 5,          26,	                                    19.95262315	
    # 4.5,        25,	                                    17.7827941	
    # 4,          24,	                                    15.84893192	
    # 3.56,       23,	                                    14.12537545	
    # 3.2,        22,	                                    12.58925412	
    # 2.8,        21,	                                    11.22018454	
    # 2.56,       20,	                                    10	        
    # 2.25,       19,	                                    8.912509381	
    # 2,          18,	                                    7.943282347	
    # 1.78,       17,	                                    7.079457844	
    # 1.6,        16,	                                    6.309573445	
    # 1.4,        15,	                                    5.623413252	
    # 1.28,       14,	                                    5.011872336	
    # 1.12,       13,	                                    4.466835922	
    # 1,          12,	                                    3.981071706	
    # 0.89,       11,	                                    3.548133892	
    # 0.8,        10,	                                    3.16227766	
    # 0.7,        9,                                        2.818382931	
    # 0.64,       8,	                                    2.511886432	
    # 0.56,       7,                                        2.238721139	
    # 0.5,        6,	                                    1.995262315	
    # 0.44,       5,	                                    1.77827941	
    # 0.4,        4,	                                    1.584893192	
    # 0.35,       3,	                                    1.412537545	
    # 0.32,       2,	                                    1.258925412	
    # 0.28,       1,	                                    1.122018454	
    # 0.25,       0,	                                    1	  

    
#general initialization stuff  
space=' '; slash='/'; 
if sys.platform=="win32": slash='\\' # not really needed, but looks nicer ;)
Program_name = os.path.basename(sys.argv[0]); 
if Program_name.find('.')>0: Program_name = Program_name[:Program_name.find('.')]
python_version=str(sys.version_info[0])+'.'+str(sys.version_info[1])+'.'+str(sys.version_info[2])
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
try: win32gui.SetForegroundWindow(win32console.GetConsoleWindow())
except: pass #silent

#read FID 
with open(FIDfile, "rb") as f: FIDrawdata= np.fromfile(f, dtype=np.int32) 
FIDrawdata_CPX = FIDrawdata[0::2] + 1j * FIDrawdata[1::2]
FIDrawdata = 0 #free memory

#read acqp file
ACQPfile=os.path.dirname(FIDfile)+slash+'acqp'
ACQPdata=ReadParamFile(ACQPfile)

#read method file
METHODfile=os.path.dirname(FIDfile)+slash+'method'
METHODdata=ReadParamFile(METHODfile)

#compatibility with Paravision 6
METHODdata["Method"] = METHODdata["Method"].replace ("<Bruker:","")
METHODdata["Method"] = METHODdata["Method"].replace (">","")


#check for not implemented stuff
if  not(METHODdata["Method"] == "FLASH" or METHODdata["Method"] == "FISP") or METHODdata["PVM_SpatDimEnum"] != "2D":
    print ('ERROR: Recon only implemented for FLASH/FISP 2D method'); 
    sys.exit(1)
if METHODdata["PVM_NSPacks"] != 1:
    print ('ERROR: Recon only implemented 1 package'); 
    sys.exit(1) 
if METHODdata["PVM_EncPpiAccel1"] != 1 or METHODdata["PVM_EncNReceivers"] != 1 or\
   METHODdata["PVM_EncZfAccel1"] != 1:
    print ('ERROR: Recon for parallel acquisition not implemented'); 
    sys.exit(1)

#start
print ('Starting recon')    

#reshape FID data according to dimensions from method file
#"order="F" means Fortran style order as by BRUKER conventions
repetitions=METHODdata["PVM_NRepetitions"]
dim=METHODdata["PVM_EncMatrix"]
dim0 = dim[0]; dim0_mod_128 = dim0%128
if dim0_mod_128!=0: dim0=(int(dim0/128)+1)*128 # Bruker sets readout point to a multiple of 128
if  METHODdata["Method"] == "FLASH":
    dim=[dim[0],METHODdata["PVM_SPackArrNSlices"],dim[1],repetitions]# insert slice dimension  and repetitions
    try: FIDrawdata_CPX = FIDrawdata_CPX.reshape(dim0,dim[1],dim[2],dim[3], order="F")
    except: print ('ERROR: k-space data reshape failed (dimension problem)'); sys.exit(1)
elif METHODdata["Method"] == "FISP":
    dim=[dim[0],dim[1],METHODdata["PVM_SPackArrNSlices"],repetitions]# insert slice dimension and repetitions
    try: FIDrawdata_CPX = FIDrawdata_CPX.reshape(dim0,dim[1],dim[2],dim[3], order="F")
    except: print ('ERROR: k-space data reshape failed (dimension problem)'); sys.exit(1)
    FIDrawdata_CPX = np.transpose (FIDrawdata_CPX, axes=(0,2,1,3))
    dim = FIDrawdata_CPX.shape
else:
    print ('Error: unknown method',METHODdata["Method"])
    sys.exit(1)
if dim0 != dim[0]: FIDrawdata_CPX = FIDrawdata_CPX[0:dim[0],:,:,:]  
  
#partial phase acquisition - add zeros
if METHODdata["PVM_EncPftAccel1"] != 1:
   zeros_ = np.zeros (shape=(dim[0],int(dim[1]*(float(METHODdata["PVM_EncPftAccel1"])-1.)),dim[2],dim[3]))
   FIDrawdata_CPX = np.append (FIDrawdata_CPX, zeros_,axis=1)
   dim=FIDrawdata_CPX.shape
   
#reorder data
FIDdata_tmp=np.empty(shape=(dim[0],dim[1],dim[2],dim[3]),dtype=np.complex64)
FIDdata=np.empty(shape=(dim[0],dim[1],dim[2],dim[3]),dtype=np.complex64)
order1=METHODdata["PVM_EncSteps1"]+dim[2]/2                             
for i in range(0,order1.shape[0]): FIDdata_tmp[:,:,order1[i],:]=FIDrawdata_CPX[:,:,i,:]
FIDrawdata_CPX = 0 #free memory  
order2=METHODdata["PVM_ObjOrderList"]
if dim[1]>1: # more than one slice
   for i in range(0,order2.shape[0]): FIDdata[:,order2[i],:,:]=FIDdata_tmp[:,i,:,:]
else: # only one slice
   FIDdata=FIDdata_tmp
FIDdata_tmp = 0 #free memory  
print('.', end='') #progress indicator

# apply FOV offsets = (linear phase in k-space)
PackArrPhase1Offset=METHODdata["PVM_SPackArrPhase1Offset"]
realFOV = METHODdata["PVM_Fov"]*METHODdata["PVM_AntiAlias"]
phase_step1 = +2.*np.pi*float(PackArrPhase1Offset)/float(realFOV[1])
mag = np.abs(FIDdata[:,:,:,:]); ph = np.angle(FIDdata[:,:,:,:])
for i in range(0,FIDdata.shape[2]): ph[:,:,i,:] -= float(i-int(FIDdata.shape[2]/2))*phase_step1
FIDdata [:,:,:,:] = mag * np.exp(1j*ph)
print('.', end='') #progress indicator

#zero fill
zero_fill=2.
SpatResol=METHODdata["PVM_SpatResol"]/zero_fill
res1=METHODdata["PVM_SPackArrSliceDistance"]
if dim[1]<=1: res1=METHODdata["PVM_SliceThick"] # only one slice
SpatResol=[SpatResol[0],res1,SpatResol[1]]# insert slice dimension
FIDdata_ZF = np.zeros(shape=(int(dim[0]*zero_fill),dim[1],
                             int(dim[2]*zero_fill),dim[3]),dtype=np.complex64)
dim0start=int(dim[0]*(zero_fill-1)/2)
dim2start=int(dim[2]*(zero_fill-1)/2)
FIDdata_ZF[dim0start:dim0start+dim[0],:,dim2start:dim2start+dim[2],:] = \
    FIDdata[0:dim[0],:,0:dim[2],:]
FIDdata=FIDdata_ZF;
FIDdata_ZF = 0 #free memory 
dim=FIDdata.shape
print('.', end='') #progress indicator

#roll partial echo (at Bruker aka echo position)
EchoPosition_raw=METHODdata["PVM_EchoPosition"]
EchoPosition_raw=50-(50-EchoPosition_raw)/zero_fill
EchoPosition=int(EchoPosition_raw/100.*dim[0])
if METHODdata["Method"] == "FISP" and METHODdata["ssfp"] == "ECHO":
      EchoPosition=dim[0]-int(EchoPosition_raw/100.*dim[0])      
FIDdata=np.roll(FIDdata, dim[0]/2-EchoPosition, axis=(0))

#find borders in case of partial echo and/or phase encoding
nz = np.asarray(np.nonzero (FIDdata))
first_x=np.amin(nz[0,:]); last_x=np.amax(nz[0,:])
first_z=np.amin(nz[2,:]); last_z=np.amax(nz[2,:])
#calculate % increase of resolution if partial fourrier recon used
percentual_inc_x=float(last_x+first_x+1-dim[0])/float(last_x-first_x)*100.
percentual_inc_z=float(last_z+first_z+1-dim[2])/float(last_z-first_z)*100.
print('.', end='') #progress indicator
print('.', end='') #progress indicator

min_percentual=10. # if the potential increase in resolution is less than this % then don't even try
if abs(percentual_inc_x)>min_percentual  or abs(percentual_inc_z)>min_percentual:
    #low pass filter for phase correction (function: 1-hanning^2)
    percentage = 10 # center only (lowpass)
    FIDlowpass = np.empty(shape=FIDdata.shape,dtype=np.complex64)
    FIDlowpass [:,:,:,:] = FIDdata [:,:,:,:]
    npoints_x = int(float(dim[0]/zero_fill)*percentage/100.)
    hanning_x = np.zeros(shape=(dim[0]),dtype=np.float32)
    x_ = np.linspace (- np.pi/2.,np.pi/2.,num=2*npoints_x+1)
    hanning_x [int(dim[0]/2)-npoints_x:int(dim[0]/2)+npoints_x+1] = 1-np.power(np.sin(x_),4)
    FIDlowpass[:,:,:,:] *= hanning_x [:,None,None,None]
    npoints_z = int(float(dim[2]/zero_fill)*percentage/100.)
    hanning_z = np.zeros(shape=(dim[2]),dtype=np.float32)
    z_ = np.linspace (-np.pi/2.,np.pi/2.,num=2*npoints_z+1)
    hanning_z [int(dim[2]/2)-npoints_z:int(dim[2]/2)+npoints_z+1] = 1-np.power(np.sin(z_),4)
    FIDlowpass[:,:,:,:] *= hanning_z [None,None,:,None]
    print('.', end='') #progress indicator
    #FFT lowpass data
    FIDlowpass = np.fft.fftshift(FIDlowpass, axes=(0,2))
    FIDlowpass = FFT2D(FIDlowpass)
    print('.', end='') #progress indicator
    #FFT actual data
    FIDdata = np.fft.fftshift(FIDdata, axes=(0,2))
    FIDdata = FFT2D(FIDdata)
    print('.', end='') #progress indicator
    # subtract phase difference from actual
    FIDlowpass = FIDdata/FIDlowpass # use this phase
    FIDdata = np.abs(FIDdata) * np.exp(1j*np.angle(FIDlowpass)) #here
    FIDlowpass = 0 # free memory
    #inverse FFT
    FIDdata = iFFT2D(FIDdata)
    FIDdata = np.fft.fftshift(FIDdata, axes=(0,2))    
    print('.', end='') #progress indicator

    # copy complex conjugates
    percentage = 5 # mix conjugate with original
    if percentual_inc_x>min_percentual: # dimension 0 points missing at the beginning        
        npoints_x = int(float(dim[0]/zero_fill)*percentage/100.)
        compl_conjugate_x = np.conj(FIDdata[dim[0]-first_x-npoints_x:dim[0],:,:,:])
        compl_conjugate_x = compl_conjugate_x[::-1,:,::-1,:] # reverse array
        compl_conjugate_x=np.roll(compl_conjugate_x, 1, axis=(2)) #symetry point in dim/2
        hanning_x = np.zeros(shape=(compl_conjugate_x.shape[0]),dtype=np.float32)
        x_ = np.linspace (1./(npoints_x-1.)*np.pi/2.,(1.-1./(npoints_x-1))*np.pi/2.,num=npoints_x)
        hanning_x [compl_conjugate_x.shape[0]-npoints_x:compl_conjugate_x.shape[0]] = np.power(np.sin(x_),2)
        FIDdata[1:first_x+1+npoints_x,:,:,:] *= hanning_x [:,None,None,None]
        hanning_x = 1.- hanning_x
        compl_conjugate_x *= hanning_x [:,None,None,None]
        #print (FIDdata[first_x+npoints_x,dim[1]/2,dim[2]/2])
        #print (compl_conjugate_x[compl_conjugate_x.shape[0]-1,dim[1]/2,dim[2]/2])
        FIDdata[1:first_x+1+npoints_x,:,:,:] += compl_conjugate_x[:,:,:,:]
        first_x=dim[0]-last_x
        compl_conjugate_x = 0 # free memory
    elif -1.*percentual_inc_x>min_percentual: # dimension 0 points missing at the end 
        npoints_x = int(float(dim[0]/zero_fill)*percentage/100.)       
        compl_conjugate_x = np.conj(FIDdata[1:dim[0]-last_x+npoints_x,:,:,:])
        compl_conjugate_x = compl_conjugate_x[::-1,:,::-1,:] # reverse array
        compl_conjugate_x=np.roll(compl_conjugate_x, 1, axis=(2)) #symetry point in dim/2
        hanning_x = np.zeros(shape=(compl_conjugate_x.shape[0]),dtype=np.float32)
        x_ = np.linspace (1./(npoints_x-1.)*np.pi/2.,(1.-1./(npoints_x-1))*np.pi/2.,num=npoints_x)
        hanning_x [compl_conjugate_x.shape[0]-npoints_x:compl_conjugate_x.shape[0]] = np.power(np.sin(x_),2)
        hanning_x = hanning_x [::-1] #reverse array
        FIDdata[last_x+1-npoints_x:dim[0],:,:,:] *= hanning_x [:,None,None,None]
        hanning_x = 1.- hanning_x
        compl_conjugate_x *= hanning_x [:,None,None,None]
        #print (FIDdata[last_x+1-npoints_x,dim[1]/2,dim[2]/2])
        #print (compl_conjugate_x[0,dim[1]/2,dim[2]/2])
        FIDdata[last_x+1-npoints_x:dim[0],:,:,:] += compl_conjugate_x[:,:,:,:]       
        last_x=dim[0]-first_x
        compl_conjugate_x = 0 # free memory     
    if percentual_inc_z>min_percentual: # dimension 2 points missing at the beginning
        npoints_z = int(float(dim[2]/zero_fill)*percentage/100.)  
        compl_conjugate_z = np.conj(FIDdata[:,:,dim[2]-first_z-npoints_z:dim[2],:])
        compl_conjugate_z = compl_conjugate_z[::-1,:,::-1,:] # reverse array
        compl_conjugate_z=np.roll(compl_conjugate_z, 1, axis=(0)) #symetry point in dim/2
        hanning_z = np.zeros(shape=(compl_conjugate_z.shape[2]),dtype=np.float32)
        y_ = np.linspace (1./(npoints_z-1.)*np.pi/2.,(1.-1./(npoints_z-1))*np.pi/2.,num=npoints_z)
        hanning_z [compl_conjugate_z.shape[2]-npoints_z:compl_conjugate_z.shape[2]] = np.power(np.sin(y_),2)
        FIDdata[:,:,1:first_z+1+npoints_z,:] *= hanning_z [None,None,:,None]
        hanning_z = 1.- hanning_z
        compl_conjugate_z *= hanning_z [None,None,:,None]
        FIDdata[:,:,1:first_z+1+npoints_z,:] += compl_conjugate_z[:,:,:,:]
        first_z=dim[2]-last_z
        compl_conjugate_z = 0 # free memory         
    elif -1.*percentual_inc_z>min_percentual: # dimension 2 points missing at the end
        npoints_z = int(float(dim[2]/zero_fill)*percentage/100.)        
        compl_conjugate_z = np.conj(FIDdata[:,:,1:dim[2]-last_z+npoints_z,:])
        compl_conjugate_z = compl_conjugate_z[::-1,:,::-1,:] # reverse array
        compl_conjugate_z=np.roll(compl_conjugate_z, 1, axis=(0)) #symetry point in dim/2
        hanning_z = np.zeros(shape=(compl_conjugate_z.shape[2]),dtype=np.float32)
        x_ = np.linspace (1./(npoints_z-1.)*np.pi/2.,(1.-1./(npoints_z-1))*np.pi/2.,num=npoints_z)
        hanning_z [compl_conjugate_z.shape[2]-npoints_z:compl_conjugate_z.shape[2]] = np.power(np.sin(x_),2)
        hanning_z = hanning_z [::-1] #reverse array
        FIDdata[:,:,last_z+1-npoints_z:dim[2],:] *= hanning_z [None,None,:,None]
        hanning_z = 1.- hanning_z
        compl_conjugate_z *= hanning_z [None,None,:,None]        
        FIDdata[:,:,last_z+1-npoints_z:dim[2],:] += compl_conjugate_z[:,:,:,:]       
        last_z=dim[2]-first_z 
        compl_conjugate_z = 0 # free memory
    print('.', end='') #progress indicator 
  
#Hanning filter
percentage = 10.
npoints_x = int(float(dim[0]/zero_fill)*percentage/100.)
hanning_x = np.zeros(shape=(dim[0]),dtype=np.float32)
x_ = np.linspace (1./(npoints_x-1.)*np.pi/2.,(1.-1./(npoints_x-1))*np.pi/2.,num=npoints_x)
hanning_x [first_x:first_x+npoints_x] = np.power(np.sin(x_),2)
hanning_x [first_x+npoints_x:last_x-npoints_x+1] = 1
x_ = x_[::-1] # reverse x_
hanning_x [last_x-npoints_x+1:last_x+1] = np.power(np.sin(x_),2)
#print (hanning_x.shape, hanning_x)
FIDdata[:,:,:,:] *= hanning_x [:,None,None,None]
npoints_z = int(float(dim[2]/zero_fill)*percentage/100.)
hanning_z = np.zeros(shape=(dim[2]),dtype=np.float32)
z_ = np.linspace (1./(npoints_z-1.)*np.pi/2.,(1.-1./(npoints_z-1))*np.pi/2.,num=npoints_z)
hanning_z [first_z:first_z+npoints_z] = np.power(np.sin(z_),2)
hanning_z [first_z+npoints_z:last_z-npoints_z+1] = 1
z_ = z_[::-1] # reverse z_
hanning_z [last_z-npoints_z+1:last_z+1] = np.power(np.sin(z_),2)
#print (hanning_z.shape, hanning_z)
FIDdata[:,:,:,:] *= hanning_z [None,None,:,None]
print('.', end='') #progress indicator      

#FFT
IMGdata=FIDdata
FIDdata = 0 #free memory 
IMGdata = np.fft.fftshift(IMGdata, axes=(0,2))
IMGdata = FFT2D(IMGdata)
IMGdata = np.fft.fftshift(IMGdata, axes=(0,2))          
print('.', end='') #progress indicator

#throw out antialiasing
crop=METHODdata["PVM_AntiAlias"]
dim0start=int((dim[0]-dim[0]/crop[0])/2)
dim2start=int((dim[2]-dim[2]/crop[1])/2)
dim0end = int(dim0start+dim[0]/crop[0])
dim2end = int(dim2start+dim[2]/crop[1])
IMGdata = IMGdata[dim0start:dim0end,:,dim2start:dim2end,:]
dim=IMGdata.shape
print('.', end='') #progress indicator

#permute dimensions
#worx for PVM_SPackArrSliceOrient=axial
#this way results are comparabled to ImageJ's BrukerOpener plugin
if METHODdata["PVM_SPackArrSliceOrient"] == "axial":
    SpatResol_perm = np.empty(shape=(3))
    SpatResol_perm[0]=SpatResol[0]
    SpatResol_perm[1]=SpatResol[2]
    SpatResol_perm[2]=SpatResol[1]
    IMGdata = np.transpose (IMGdata, axes=(0,2,1,3))
    temp=SpatResol_perm[0]
    SpatResol_perm[0]=SpatResol_perm[1]
    SpatResol_perm[1]=temp
    if METHODdata["PVM_SPackArrReadOrient"] == "A_P":
        temp=SpatResol_perm[0]
        SpatResol_perm[0]=SpatResol_perm[1]
        SpatResol_perm[1]=temp
        IMGdata = np.rot90(IMGdata, k=1, axes=(0,1)) # rotate (axial A_P)
    elif METHODdata["PVM_SPackArrReadOrient"] == "L_R":
        IMGdata = IMGdata[::-1,:,:,:] # flip axis (axial L_R)
    else:
        SpatResol_perm=SpatResol    
        print ('Warning: unknown Orientation',METHODdata["PVM_SPackArrSliceOrient"],
                METHODdata["PVM_SPackArrReadOrient"]);
        print ('         resulting images may be rotated incorrectly');
elif METHODdata["PVM_SPackArrSliceOrient"] == "sagittal":
    if METHODdata["PVM_SPackArrReadOrient"] == "A_P":
        SpatResol_perm = np.empty(shape=(3))
        SpatResol_perm[0]=SpatResol[0]
        SpatResol_perm[1]=SpatResol[2]
        SpatResol_perm[2]=SpatResol[1]    
        IMGdata = np.transpose (IMGdata, axes=(0,2,1,3))
        IMGdata = IMGdata[::-1,:,:,:] # flip axis
    elif METHODdata["PVM_SPackArrReadOrient"] == "H_F":
        SpatResol_perm = np.empty(shape=(3))
        SpatResol_perm[0]=SpatResol[0]
        SpatResol_perm[1]=SpatResol[2]
        SpatResol_perm[2]=SpatResol[1]    
        IMGdata = np.transpose (IMGdata, axes=(0,2,1,3))
        IMGdata = IMGdata[::-1,:,:,:] # flip axis         
    else:
        SpatResol_perm=SpatResol
        print ('Warning: unknown Orientation',METHODdata["PVM_SPackArrSliceOrient"],
                METHODdata["PVM_SPackArrReadOrient"]);
        print ('         resulting images may be rotated incorrectly');
elif METHODdata["PVM_SPackArrSliceOrient"] == "coronal":
    if METHODdata["PVM_SPackArrReadOrient"] == "L_R":
        SpatResol_perm = np.empty(shape=(3))
        SpatResol_perm[0]=SpatResol[0]
        SpatResol_perm[1]=SpatResol[2]
        SpatResol_perm[2]=SpatResol[1]
        IMGdata = np.transpose (IMGdata, axes=(0,2,1,3))        
        IMGdata = IMGdata[::-1,:,::-1,::-1] # flip axis 
    elif METHODdata["PVM_SPackArrReadOrient"] == "H_F":
        SpatResol_perm = np.empty(shape=(3))
        SpatResol_perm[0]=SpatResol[0]
        SpatResol_perm[1]=SpatResol[2]
        SpatResol_perm[2]=SpatResol[1]
        IMGdata = np.transpose (IMGdata, axes=(0,2,1,3))
        IMGdata = np.rot90(IMGdata, k=1, axes=(0,1)) # rotate (90 deg)          
    else:
        SpatResol_perm=SpatResol
        print ('Warning: unknown Orientation',METHODdata["PVM_SPackArrSliceOrient"],
                METHODdata["PVM_SPackArrReadOrient"]);
        print ('         resulting images may be rotated incorrectly');                    
else:
    SpatResol_perm=SpatResol
    print ('Warning: unknown Orientation',METHODdata["PVM_SPackArrSliceOrient"],
            METHODdata["PVM_SPackArrReadOrient"]);
    print ('         resulting images may be rotated incorrectly'); 
print('.', end='') #progress indicator

#find noise mask threshold from histogram
#n_points=IMGdata.shape[0]*IMGdata.shape[1]*IMGdata.shape[2]
#steps=int(n_points/1000); start=0; fin=np.max(np.abs(IMGdata[:,:,:]))
#xbins =  np.linspace(start,fin,steps)
#ybins, binedges = np.histogram(np.abs(IMGdata[:,:,:]), bins=xbins)
#ybins = np.resize (ybins,len(xbins)); ybins[len(ybins)-1]=0
#ybins = smooth(ybins,steps/20)
#--- old code find minimum ---
#start=ybins.argmax()
#i=start;minx=0;miny=ybins[start]
#while i<len(ybins):
#    i+=1
#    if ybins[i]<=miny: miny=ybins[i]; minx=i; 
#    else: i=len(ybins);
#threshold=xbins[minx]
#--- new code find FWHM ---
#start=ybins.argmax(); i=start
#while i<len(ybins):
#    i+=1
#    if ybins[i]<np.max(ybins)/2: 
#        threshold=xbins[start]+4.0*(xbins[i]-xbins[start]);
#        i=len(ybins);        
#mask =  abs(IMGdata [:,:,:]) > threshold  
#enable the following view histogram plot
#print ('\nThreshold = %.2e' % threshold)
#import pylab; pylab.plot(xbins,ybins, linewidth=1.5); pylab.draw();
#pylab.show(block=False); os.system("pause"); pylab.close(); 

# use noise in all 8 corners to establish threshold
N=10 # use 10% at the corners of the FOV
tresh=np.empty(shape=8,dtype=np.float)
avg=np.empty(shape=8,dtype=np.float)
std=np.empty(shape=8,dtype=np.float)
xstart=0; xend=int(ceil(float(IMGdata.shape[0])/float(N)))
ystart=0; yend=int(ceil(float(IMGdata.shape[1])/float(N)))
zstart=0; zend=int(ceil(float(IMGdata.shape[2])/float(N)))
arr=np.abs(IMGdata[xstart:xend,ystart:yend,zstart:zend,:])
avg[0]=np.mean(arr)
std[0]=np.std(arr)
tresh[0]=avg[0] + 4*std[0]
xstart=int(floor(float(IMGdata.shape[0])-float(IMGdata.shape[0])/float(N))); xend=IMGdata.shape[0]
ystart=0; yend=int(ceil(float(IMGdata.shape[1])/float(N)))
zstart=0; zend=int(ceil(float(IMGdata.shape[2])/float(N)))
arr=np.abs(IMGdata[xstart:xend,ystart:yend,zstart:zend,:])
avg[1]=np.mean(arr)
std[1]=np.std(arr)
tresh[1]=avg[1] + 4*std[1]
xstart=0; xend=int(ceil(float(IMGdata.shape[0])/float(N)))
ystart=int(floor(float(IMGdata.shape[1])-float(IMGdata.shape[1])/float(N))); yend=IMGdata.shape[1]
zstart=0; zend=int(ceil(float(IMGdata.shape[2])/float(N)))
arr=np.abs(IMGdata[xstart:xend,ystart:yend,zstart:zend,:])
avg[2]=np.mean(arr)
std[2]=np.std(arr)
tresh[2]=avg[2] + 4*std[2]
xstart=int(floor(float(IMGdata.shape[0])-float(IMGdata.shape[0])/float(N))); xend=IMGdata.shape[0]
ystart=int(floor(float(IMGdata.shape[1])-float(IMGdata.shape[1])/float(N))); yend=IMGdata.shape[1]
zstart=0; zend=int(ceil(float(IMGdata.shape[2])/float(N)))
arr=np.abs(IMGdata[xstart:xend,ystart:yend,zstart:zend,:])
avg[3]=np.mean(arr)
std[3]=np.std(arr)
tresh[3]=avg[3] + 4*std[3]
xstart=0; xend=int(ceil(float(IMGdata.shape[0])/float(N)))
ystart=0; yend=int(ceil(float(IMGdata.shape[1])/float(N)))
zstart=int(floor(float(IMGdata.shape[2])-float(IMGdata.shape[2])/float(N))); zend=IMGdata.shape[2]
arr=np.abs(IMGdata[xstart:xend,ystart:yend,zstart:zend,:])
avg[4]=np.mean(arr)
std[4]=np.std(arr)
tresh[4]=avg[4] + 4*std[4]
xstart=int(floor(float(IMGdata.shape[0])-float(IMGdata.shape[0])/float(N))); xend=IMGdata.shape[0]
ystart=0; yend=int(ceil(float(IMGdata.shape[1])/float(N)))
zstart=int(floor(float(IMGdata.shape[2])-float(IMGdata.shape[2])/float(N))); zend=IMGdata.shape[2]
arr=np.abs(IMGdata[xstart:xend,ystart:yend,zstart:zend])
avg[5]=np.mean(arr)
std[5]=np.std(arr)
tresh[5]=avg[5] + 4*std[5]
xstart=0; xend=int(ceil(float(IMGdata.shape[0])/float(N)))
ystart=int(floor(float(IMGdata.shape[1])-float(IMGdata.shape[1])/float(N))); yend=IMGdata.shape[1]
zstart=int(floor(float(IMGdata.shape[2])-float(IMGdata.shape[2])/float(N))); zend=IMGdata.shape[2]
arr=np.abs(IMGdata[xstart:xend,ystart:yend,zstart:zend,:])
avg[6]=np.mean(arr)
std[6]=np.std(arr)
tresh[6]=avg[6] + 4*std[6]
xstart=int(floor(float(IMGdata.shape[0])-float(IMGdata.shape[0])/float(N))); xend=IMGdata.shape[0]
ystart=int(floor(float(IMGdata.shape[1])-float(IMGdata.shape[1])/float(N))); yend=IMGdata.shape[1]
zstart=int(floor(float(IMGdata.shape[2])-float(IMGdata.shape[2])/float(N))); zend=IMGdata.shape[2]
arr=np.abs(IMGdata[xstart:xend,ystart:yend,zstart:zend,:])
avg[7]=np.mean(arr)
std[7]=np.std(arr)
tresh[7]=avg[7] + 4*std[7]
threshold=np.min(tresh)
mask =  abs(IMGdata [:,:,:,:]) > threshold

# calculate phase difference
if repetitions > 1:
   IMGdata_rawPH = np.zeros(shape=IMGdata.shape,dtype=np.complex64)
   for i in range(0,dim[3]): 
      IMGdata_rawPH[:,:,:,i] = IMGdata[:,:,:,i]/IMGdata[:,:,:,0]
else:
   IMGdata_rawPH = IMGdata
      
#transform to int
ReceiverGain = ACQPdata["RG"] # RG is a simple attenuation FACTOR, NOT in dezibel (dB) unit !!!
n_Averages = METHODdata["PVM_NAverages"]
IMGdata_ABS = np.abs(IMGdata)/RG_to_voltage(ReceiverGain)/n_Averages; 
max_ABS = np.amax(IMGdata_ABS);
IMGdata_ABS *= 32767./max_ABS
IMGdata_ABS = IMGdata_ABS.astype(np.int16)
print('.', end='') #progress indicator
IMGdata_PH  = np.angle(IMGdata_rawPH)*mask; # use this to mask out background noise
max_PH = np.pi; 
IMGdata_PH *= 32767./max_PH
IMGdata_PH = IMGdata_PH.astype(np.int16)
# set max/min in 0,0,0/1,1,1 corners (this save the venc parameterin the image)
IMGdata_PH [0,0,0,:] = 32767
if IMGdata_PH.shape[0]>1 and IMGdata_PH.shape[1]>1 and IMGdata_PH.shape[2]>1: IMGdata_PH [1,1,1,:] = -32767 
print('.', end='') #progress indicator

#save NIFTI
aff = np.eye(4)
aff[0,0] = SpatResol_perm[0]*1000; aff[0,3] = -(IMGdata.shape[0]/2)*aff[0,0]
aff[1,1] = SpatResol_perm[1]*1000; aff[1,3] = -(IMGdata.shape[1]/2)*aff[1,1]
aff[2,2] = SpatResol_perm[2]*1000; aff[2,3] = -(IMGdata.shape[2]/2)*aff[2,2]
NIFTIimg_ABS = nib.Nifti1Image(IMGdata_ABS, aff)
NIFTIimg_ABS.header.set_slope_inter(max_ABS/32767.,0)
NIFTIimg_ABS.header.set_xyzt_units(3, 8)
NIFTIimg_ABS.set_sform(aff, code=0)
NIFTIimg_ABS.set_qform(aff, code=1)
NIFTIimg_ABS_masked = nib.Nifti1Image(IMGdata_ABS*mask, aff)
NIFTIimg_ABS_masked.header.set_slope_inter(max_ABS/32767.,0)
NIFTIimg_ABS_masked.header.set_xyzt_units(3, 8)
NIFTIimg_ABS_masked.set_sform(aff, code=0)
NIFTIimg_ABS_masked.set_qform(aff, code=1)
NIFTIimg_PH  = nib.Nifti1Image(IMGdata_PH, aff)
NIFTIimg_PH.header.set_slope_inter(max_PH/32767.,0)
NIFTIimg_PH.header.set_xyzt_units(3, 8)
NIFTIimg_PH.set_sform(aff, code=0)
NIFTIimg_PH.set_qform(aff, code=1)
#write
try:
    print('.', end='') #progress indicator
    nib.save(NIFTIimg_ABS, os.path.join(os.path.dirname(FIDfile),OrigFilename+'_MAGNT.nii.gz'))
    print('.', end='') #progress indicator
    nib.save(NIFTIimg_ABS_masked, os.path.join(os.path.dirname(FIDfile),OrigFilename+'_MAG_m.nii.gz'))    
    print('.', end='') #progress indicator
    nib.save(NIFTIimg_PH , os.path.join(os.path.dirname(FIDfile),OrigFilename+'_PHASE.nii.gz'))
except:
    print ('\nERROR:  problem while writing results'); sys.exit(1)
print ('\nSuccessfully written output files '+OrigFilename+'_MAGNT/PHASE.nii.gz')   

#end
if sys.platform=="win32": os.system("pause") # windows
else: 
    #os.system('read -s -n 1 -p "Press any key to continue...\n"')
    import termios
    print("Press any key to continue...")
    fd = sys.stdin.fileno()
    oldterm = termios.tcgetattr(fd)
    newattr = termios.tcgetattr(fd)
    newattr[3] = newattr[3] & ~termios.ICANON & ~termios.ECHO
    termios.tcsetattr(fd, termios.TCSANOW, newattr)
    try: result = sys.stdin.read(1)
    except IOError: pass
    finally: termios.tcsetattr(fd, termios.TCSAFLUSH, oldterm)