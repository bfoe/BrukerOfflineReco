#
# reads Bruker MR data (Paravision v5.1)
# reconstructs images from raw acquisition data (FID files)
# this version is for the MSME 2D method only
#
# outputs 4D NIFTI containing all individual echo images
# and a 3D NIFTI containing the average of all echoes
# does not calculate T2 maps, use "NIFTI_T2map" for that
#
# ----- VERSION HISTORY -----
#
# Version 0.1 - 19, February 2019
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
#    - scipy
#    - nibabel
#

from __future__ import print_function
try: import win32gui, win32console
except: pass #silent
from math import ceil, floor
import sys
import os
import numpy as np
import nibabel as nib
if getattr( sys, 'frozen', False ): # running as pyinstaller bundle
   from scipy_extract import median_filter
   from scipy_extract import gaussian_filter 
else: # running native python
   from scipy.ndimage import median_filter 
   from scipy.ndimage import gaussian_filter 
   

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
    
def FFT2D (array):
    for k in range(0,array.shape[1]): array[:,k,:,:] = np.fft.fft(array[:,k,:,:], axis=(0))
    for k in range(0,array.shape[1]): array[:,k,:,:] = np.fft.fft(array[:,k,:,:], axis=(1))        
    return array 

def iFFT2D (array):
    for k in range(0,array.shape[1]): array[:,k,:,:] = np.fft.ifft(array[:,k,:,:], axis=(0))
    for k in range(0,array.shape[1]): array[:,k,:,:] = np.fft.ifft(array[:,k,:,:], axis=(1))   
    return array     
    
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
with open(FIDfile, "r") as f: FIDrawdata= np.fromfile(f, dtype=np.int32) 
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
if  METHODdata["Method"] != "MSME"  or METHODdata["PVM_SpatDimEnum"] != "2D":
    print ('ERROR: Recon only implemented for MSME 2D method'); 
    sys.exit(1)
if METHODdata["PVM_NSPacks"] != 1:
    print ('ERROR: Recon only implemented 1 package'); 
    sys.exit(1) 
if METHODdata["PVM_NRepetitions"] != 1:
    print ('ERROR: Recon only implemented 1 repetition'); 
    sys.exit(1)    
if METHODdata["PVM_EncPpiAccel1"] != 1 or METHODdata["PVM_EncNReceivers"] != 1 or\
   METHODdata["PVM_EncZfAccel1"] != 1:
    print ('ERROR: Recon for parallel acquisition not implemented'); 
    sys.exit(1)
if METHODdata["constNEchoes"] != "Yes"]:
    print ('ERROR: non constNEchoes not implemented'); 
    sys.exit(1)
if METHODdata["PVM_NEchoImages"] != METHODdata["NEchoes"]:
    print ('ERROR: unequal #echoes'); 
    sys.exit(1)         
NEchoes = METHODdata["NEchoes"]
TEs = METHODdata["EffectiveTE"].astype(np.float32)

#start
print ('Starting recon')    

#reshape FID data according to dimensions from method file
#"order="F" means Fortran style order as by BRUKER conventions
dim=METHODdata["PVM_EncMatrix"]
dim0 = dim[0]; dim0_mod_128 = dim0%128
if dim0_mod_128!=0: dim0=(int(dim0/128)+1)*128 # Bruker sets readout point to a multiple of 128
dim=[dim[0],METHODdata["PVM_SPackArrNSlices"],dim[1]]# insert slice dimension
try: FIDrawdata_CPX = FIDrawdata_CPX.reshape(dim0,NEchoes,dim[1],dim[2], order="F")
except: print ('ERROR: k-space data reshape failed (dimension problem)'); sys.exit(1)
if dim0 != dim[0]: FIDrawdata_CPX = FIDrawdata_CPX[0:dim[0],:,:]    

#transpose such that echo dimension is last
FIDrawdata_CPX = np.transpose (FIDrawdata_CPX, axes=(0,2,3,1))
dim=FIDrawdata_CPX.shape

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


#calc phase correction
#
# the reason for this is that eddy currents produce a
# shift in the images along (mainly) in readout direction
# which depends on the echo number, as a result uncorrected
# images do not allign varying the echo number
# this is uncritical for the average magnitude image
# but does produce bad artifacts when trying to fit T2 maps
#
# this filter corrcts for slow varying phase error in k-space
# which in image space relates basically to correcting the 
# object's position error
# (linear phase in time domain = shift in frequency domain)
#
#
for j in range (3):
    AVG = np.average(FIDdata[:,:,:,:],axis=3)   #complex average all echoes
    PH = FIDdata[:,:,:,:]/AVG[:,:,:,None]       #complex division = phase difference
    PH = np.angle(PH)                           #calc phase angle   
    for i in range(dim[3]):    
       PH[:,:,:,i] = median_filter  (PH[:,:,:,i], size = (5,5,5))    #median filter 1   
       PH[:,:,:,i] = median_filter  (PH[:,:,:,i], size = (5,5,5))    #median filter 2   
       PH[:,:,:,i] = gaussian_filter(PH[:,:,:,i], sigma=10, truncate=3)
       print('.', end='') #progress indicator    
    FIDdata [:,:,:,:] = np.abs(FIDdata) * np.exp(1j*(np.angle(FIDdata)-PH)) #apply correction
    #save images to check visually
    #aff = np.eye(4)
    #IMG = nib.Nifti1Image(PH, aff)
    #nib.save(IMG, os.path.join(os.path.dirname(FIDfile),OrigFilename+'_testPH'+str(j+1)+'.nii.gz'))
AVG = 0; PH = 0 #free memory


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
IMGdata[:,:,:,:] = FFT2D(IMGdata[:,:,:,:])
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
else:
    SpatResol_perm=SpatResol
    print ('Warning: unknown Orientation',METHODdata["PVM_SPackArrSliceOrient"],
            METHODdata["PVM_SPackArrReadOrient"]);
    print ('         resulting images may be rotated incorrectly'); 
print('.', end='') #progress indicator


# use noise in all 8 corners to establish threshold
N=10 # use 10% at the corners of the FOV
tresh=np.empty(shape=8,dtype=np.float)
avg=np.empty(shape=8,dtype=np.float)
std=np.empty(shape=8,dtype=np.float)
xstart=0; xend=int(IMGdata.shape[0]/N)
ystart=0; yend=int(IMGdata.shape[1]/N)
zstart=0; zend=int(ceil(float(IMGdata.shape[2])/float(N)))
arr=np.abs(IMGdata[xstart:xend,ystart:yend,zstart:zend,:])
avg[0]=np.mean(arr)
std[0]=np.std(arr)
tresh[0]=avg[0] + 4*std[0]
xstart=int(IMGdata.shape[0]-IMGdata.shape[0]/N); xend=IMGdata.shape[0]
ystart=0; yend=int(IMGdata.shape[1]/N)
zstart=0; zend=int(ceil(float(IMGdata.shape[2])/float(N)))
arr=np.abs(IMGdata[xstart:xend,ystart:yend,zstart:zend,:])
avg[1]=np.mean(arr)
std[1]=np.std(arr)
tresh[1]=avg[1] + 4*std[1]
xstart=0; xend=int(IMGdata.shape[0]/N)
ystart=int(IMGdata.shape[1]-IMGdata.shape[1]/N); yend=IMGdata.shape[1]
zstart=0; zend=int(ceil(float(IMGdata.shape[2])/float(N)))
arr=np.abs(IMGdata[xstart:xend,ystart:yend,zstart:zend,:])
avg[2]=np.mean(arr)
std[2]=np.std(arr)
tresh[2]=avg[2] + 4*std[2]
xstart=int(IMGdata.shape[0]-IMGdata.shape[0]/N); xend=IMGdata.shape[0]
ystart=int(IMGdata.shape[1]-IMGdata.shape[1]/N); yend=IMGdata.shape[1]
zstart=0; zend=int(ceil(float(IMGdata.shape[2])/float(N)))
arr=np.abs(IMGdata[xstart:xend,ystart:yend,zstart:zend,:])
avg[3]=np.mean(arr)
std[3]=np.std(arr)
tresh[3]=avg[3] + 4*std[3]
xstart=0; xend=int(IMGdata.shape[0]/N)
ystart=0; yend=int(IMGdata.shape[1]/N)
zstart=int(floor(float(IMGdata.shape[2])-float(IMGdata.shape[2])/float(N))); zend=IMGdata.shape[2]
arr=np.abs(IMGdata[xstart:xend,ystart:yend,zstart:zend,:])
avg[4]=np.mean(arr)
std[4]=np.std(arr)
tresh[4]=avg[4] + 4*std[4]
xstart=int(IMGdata.shape[0]-IMGdata.shape[0]/N); xend=IMGdata.shape[0]
ystart=0; yend=int(IMGdata.shape[1]/N)
zstart=int(floor(float(IMGdata.shape[2])-float(IMGdata.shape[2])/float(N))); zend=IMGdata.shape[2]
arr=np.abs(IMGdata[xstart:xend,ystart:yend,zstart:zend])
avg[5]=np.mean(arr)
std[5]=np.std(arr)
tresh[5]=avg[5] + 4*std[5]
xstart=0; xend=int(IMGdata.shape[0]/N)
ystart=int(IMGdata.shape[1]-IMGdata.shape[1]/N); yend=IMGdata.shape[1]
zstart=int(floor(float(IMGdata.shape[2])-float(IMGdata.shape[2])/float(N))); zend=IMGdata.shape[2]
arr=np.abs(IMGdata[xstart:xend,ystart:yend,zstart:zend,:])
avg[6]=np.mean(arr)
std[6]=np.std(arr)
tresh[6]=avg[6] + 4*std[6]
xstart=int(IMGdata.shape[0]-IMGdata.shape[0]/N); xend=IMGdata.shape[0]
ystart=int(IMGdata.shape[1]-IMGdata.shape[1]/N); yend=IMGdata.shape[1]
zstart=int(floor(float(IMGdata.shape[2])-float(IMGdata.shape[2])/float(N))); zend=IMGdata.shape[2]
arr=np.abs(IMGdata[xstart:xend,ystart:yend,zstart:zend,:])
avg[7]=np.mean(arr)
std[7]=np.std(arr)
tresh[7]=avg[7] + 4*std[7]
threshold=np.min(tresh)
mask =  abs(IMGdata [:,:,:,:]) > threshold
      
#transform to int
ReceiverGain = ACQPdata["RG"] # RG is a simple attenuation FACTOR, NOT in dezibel (dB) unit !!!
n_Averages = METHODdata["PVM_NAverages"]
IMGdata_ABS = np.abs(IMGdata)/ReceiverGain/n_Averages; 
max_ABS = np.amax(IMGdata_ABS);
IMGdata_ABS *= 32767./max_ABS
IMGdata_ABS = IMGdata_ABS.astype(np.int16)
print('.', end='') #progress indicator
IMGdata_PH  = np.angle(IMGdata)*mask; # use this to mask out background noise
max_PH = np.pi; 
IMGdata_PH *= 32767./max_PH
IMGdata_PH = IMGdata_PH.astype(np.int16)
# set max/min in 0,0,0/1,1,1 corners (this save the venc parameterin the image)
IMGdata_PH [0,0,0,:] = 32767
if IMGdata_PH.shape[0]>1 and IMGdata_PH.shape[1]>1 and IMGdata_PH.shape[2]>1: IMGdata_PH [1,1,1,:] = -32767 
print('.', end='') #progress indicator

# calc sum of all echoes
#IMGdata_AVG = np.average(np.abs(IMGdata[:,:,:,1::2]),axis=3) #odd
#IMGdata_AVG = np.average(np.abs(IMGdata[:,:,:,::2]),axis=3)  #even
#IMGdata_AVG = np.average(np.abs(IMGdata),axis=3)  #all
IMGdata_AVG = np.abs(np.average(IMGdata,axis=3))  #complex average all echoes
IMGdata_AVG = IMGdata_AVG/ReceiverGain/n_Averages;
max_AVG = np.amax(IMGdata_AVG);
IMGdata_AVG *= 32767./max_AVG
IMGdata_AVG = IMGdata_AVG.astype(np.int16)
print('.', end='') #progress indicator

#save NIFTI
aff = np.eye(4)
aff[0,0] = SpatResol_perm[0]*1000; aff[0,3] = -(IMGdata.shape[0]/2)*aff[0,0]
aff[1,1] = SpatResol_perm[1]*1000; aff[1,3] = -(IMGdata.shape[1]/2)*aff[1,1]
aff[2,2] = SpatResol_perm[2]*1000; aff[2,3] = -(IMGdata.shape[2]/2)*aff[2,2]
NIFTIimg_AVG = nib.Nifti1Image(IMGdata_AVG, aff)
NIFTIimg_AVG.header.set_slope_inter(max_AVG/32767.,0)
NIFTIimg_AVG.header.set_xyzt_units(3, 8)
NIFTIimg_AVG.set_sform(aff, code=0)
NIFTIimg_AVG.set_qform(aff, code=1)
NIFTIimg_ABS = nib.Nifti1Image(IMGdata_ABS, aff)
NIFTIimg_ABS.header['descrip'] = "TE1 = "+str(TEs[0])
NIFTIimg_ABS.header.set_slope_inter(max_ABS/32767.,0)
NIFTIimg_ABS.header.set_xyzt_units(3, 8)
NIFTIimg_ABS.set_sform(aff, code=0)
NIFTIimg_ABS.set_qform(aff, code=1)
NIFTIimg_ABS_masked = nib.Nifti1Image(IMGdata_ABS*mask, aff)
NIFTIimg_ABS_masked.header['descrip'] = "TE1 = "+str(TEs[0])
NIFTIimg_ABS_masked.header.set_slope_inter(max_ABS/32767.,0)
NIFTIimg_ABS_masked.header.set_xyzt_units(3, 8)
NIFTIimg_ABS_masked.set_sform(aff, code=0)
NIFTIimg_ABS_masked.set_qform(aff, code=1)
NIFTIimg_PH  = nib.Nifti1Image(IMGdata_PH, aff)
NIFTIimg_PH.header['descrip'] = "TE1 = "+str(TEs[0])
NIFTIimg_PH.header.set_slope_inter(max_PH/32767.,0)
NIFTIimg_PH.header.set_xyzt_units(3, 8)
NIFTIimg_PH.set_sform(aff, code=0)
NIFTIimg_PH.set_qform(aff, code=1)
#write
try:
    print('.', end='') #progress indicator
    nib.save(NIFTIimg_AVG, os.path.join(os.path.dirname(FIDfile),OrigFilename+'_AVG.nii.gz')) 
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