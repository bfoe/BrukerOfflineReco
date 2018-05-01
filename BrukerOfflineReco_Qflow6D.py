#
# reads Bruker MR data (Paravision v5.1)
# reconstructs images from raw acquisition data (FID files)
# this version is for the FLOWMAP method only
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
#    - scipy
#

from __future__ import print_function
try: import win32gui, win32console
except: pass #silent
import sys
import os
import numpy as np
import nibabel as nib
if getattr( sys, 'frozen', False ): # running as pyinstaller bundle
   from scipy_extract import zoom
   from scipy_extract import median_filter
   from scipy_extract import gaussian_filter
else: # running native python
   from scipy.ndimage import zoom 
   from scipy.ndimage import median_filter 
   from scipy.ndimage import gaussian_filter    
pyfftw_installed = True
try: 
    import pyfftw
    test = pyfftw.interfaces.numpy_fft.fft(np.asarray([0,1,2,3]))
except: pyfftw_installed = False

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

def smooth(x,window_len):
    w=np.hanning(window_len)
    s=np.r_[2*x[0]-x[window_len-1::-1],x,2*x[-1]-x[-1:-window_len:-1]]
    w=np.hanning(window_len)
    y=np.convolve(w/w.sum(),s,mode='same')
    return y[window_len:-window_len+1]      

def is_powerof_2(n):
    return bool(n and not (n&(n-1)))    
def next_powerof_2(n):
    list = [1 << i for i in range(17)]  # 1..65536
    i=0
    while n>list[i]: i += 1
    return (list[i])

def FFT3D (array):
    if pyfftw_installed: 
       array = pyfftw.interfaces.numpy_fft.fftn(array, axes=(0,1,2))
    else: 
        for k in range(0,array.shape[1]): array[:,k,:] = np.fft.fft(array[:,k,:], axis=(0))
        for i in range(0,array.shape[0]): array[i,:,:] = np.fft.fft(array[i,:,:], axis=(0))
        for i in range(0,array.shape[0]): array[i,:,:] = np.fft.fft(array[i,:,:], axis=(1))          
    return array 

def iFFT3D (array):
    if pyfftw_installed: 
       array = pyfftw.interfaces.numpy_fft.ifftn(array, axes=(0,1,2))
    else: 
        for k in range(0,array.shape[1]): array[:,k,:] = np.fft.ifft(array[:,k,:], axis=(0))
        for i in range(0,array.shape[0]): array[i,:,:] = np.fft.ifft(array[i,:,:], axis=(0))
        for i in range(0,array.shape[0]): array[i,:,:] = np.fft.ifft(array[i,:,:], axis=(1))          
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
    
#intercatively choose input FID files (flow and zero flow reference)
FIDfile = askopenfilename(title="Choose Bruker FID file (main acq with flow)", filetypes=[("FID files","fid")])
if FIDfile == "": print ('ERROR: No FID input file specified'); sys.exit(2)
FIDfile = os.path.abspath(FIDfile)
haveRef = True 
FIDfileRef = askopenfilename(title="Choose Bruker FID file (reference zero flow)", filetypes=[("FID files","fid")])
if FIDfileRef == "": haveRef = False
else: FIDfileRef = os.path.abspath(FIDfileRef) 
TKwindows.update()
try: win32gui.SetForegroundWindow(win32console.GetConsoleWindow())
except: pass #silent

#read FID (main & ref)
with open(FIDfile, "r") as f: FIDrawdata= np.fromfile(f, dtype=np.int32) 
FIDrawdata_CPX = FIDrawdata[0::2] + 1j * FIDrawdata[1::2]
if haveRef:
    with open(FIDfileRef, "r") as f: FIDrawdata= np.fromfile(f, dtype=np.int32) 
    FIDrawdata_CPX_REF = FIDrawdata[0::2] + 1j * FIDrawdata[1::2]
FIDrawdata = 0 #free memory

#read acqp file (main & ref)
if haveRef:
    ACQPfileRef=os.path.dirname(FIDfileRef)+slash+'acqp'
    ACQPdataRef=ReadParamFile(ACQPfileRef)
ACQPfile=os.path.dirname(FIDfile)+slash+'acqp'
ACQPdata=ReadParamFile(ACQPfile)

#read method file
if haveRef:
    METHODfileRef=os.path.dirname(FIDfileRef)+slash+'method'
    METHODdataRef=ReadParamFile(METHODfileRef)
METHODfile=os.path.dirname(FIDfile)+slash+'method'
METHODdata=ReadParamFile(METHODfile)

#check for not implemented stuff
if METHODdata["Method"] != "FLOWMAP" or METHODdata["FlowMode"] != "VelocityMapping":
    print ('ERROR: Recon only implemented for FLOWMAP VelocityMapping'); 
    sys.exit(1)
if METHODdata["PVM_SpatDimEnum"] !="3D" or METHODdata["FlowEncodingDirection"] != "AllDirections":
    print ('ERROR: Recon only implemented for 3D acquisition'); 
    print ('       with flow incoding in all directions');
    sys.exit(1)
if METHODdata["FlowEncLoop"] !=4:
    print ('ERROR: ops, expected flow encoding loop = 4 '); 
    sys.exit(1)
if METHODdata["PVM_NSPacks"] != 1:
    print ('ERROR: Recon only implemented 1 package'); 
    sys.exit(1)
if METHODdata["PVM_NRepetitions"] != 1:
    print ('ERROR: Recon only implemented 1 repetition'); 
    sys.exit(1)
if METHODdata["PVM_EncPpiAccel1"] != 1 or METHODdata["PVM_EncPftAccel1"] != 1 or \
   METHODdata["PVM_EncZfAccel1"] != 1 or METHODdata["PVM_EncZfAccel2"] != 1 or \
   METHODdata["PVM_EncTotalAccel"] != 1 or METHODdata["PVM_EncNReceivers"] != 1:
    print ('ERROR: Recon for parallel acquisition not implemented'); sys.exit(1)


#check for identical parameters in main and reference
if haveRef:
    #sequence parameters
    if METHODdata["Method"] != METHODdataRef["Method"] or METHODdata["FlowMode"] != METHODdataRef["FlowMode"]:
        print ('ERROR: Parameter mismatch beteen Main and Reference files (Acquisition Method)'); 
        sys.exit(1)
    if METHODdata["PVM_EchoTime"] != METHODdataRef["PVM_EchoTime"]:
        print ('ERROR: Parameter mismatch beteen Main and Reference files (TE)'); 
        sys.exit(1)
    if METHODdata["PVM_RepetitionTime"] != METHODdataRef["PVM_RepetitionTime"]:
        print ('ERROR: Parameter mismatch beteen Main and Reference files (TR)'); 
        sys.exit(1)
    if METHODdata["FlowRange"] != METHODdataRef["FlowRange"]:
        print ('ERROR: Parameter mismatch beteen Main and Reference files (venc)'); 
        sys.exit(1)  
    if METHODdata["PVM_EffSWh"] != METHODdataRef["PVM_EffSWh"]:
        print ('ERROR: Parameter mismatch beteen Main and Reference files (Bandwidth)'); 
        sys.exit(1)  
    #geometry parameters
    if METHODdata["PVM_SPackArrSliceOrient"] != METHODdataRef["PVM_SPackArrSliceOrient"]:
        print ('ERROR: Parameter mismatch beteen Main and Reference files (Slice orientation)'); 
        sys.exit(1)
    if METHODdata["PVM_SPackArrReadOrient"] != METHODdataRef["PVM_SPackArrReadOrient"]:
        print ('ERROR: Parameter mismatch beteen Main and Reference files (Readout direction)'); 
        sys.exit(1)
    if METHODdata["PVM_SPackArrPhase1Offset"] != METHODdataRef["PVM_SPackArrPhase1Offset"]:
        print ('ERROR: Parameter mismatch beteen Main and Reference files (offset Phase)'); 
        sys.exit(1)
    if METHODdata["PVM_SPackArrSliceOffset"] != METHODdataRef["PVM_SPackArrSliceOffset"]:
        print ('ERROR: Parameter mismatch beteen Main and Reference files (offset Slice)'); 
        sys.exit(1)
    if METHODdata["PVM_SPackArrReadOffset"] != METHODdataRef["PVM_SPackArrReadOffset"]:
        print ('ERROR: Parameter mismatch beteen Main and Reference files (offset Read)'); 
        sys.exit(1)    
    if str(METHODdata["PVM_Fov"]) != str(METHODdataRef["PVM_Fov"]):
        print ('ERROR: Parameter mismatch beteen Main and Reference files (FOV)'); 
        sys.exit(1)     
    if str(METHODdata["PVM_AntiAlias"]) != str(METHODdataRef["PVM_AntiAlias"]):
        print ('ERROR: Parameter mismatch beteen Main and Reference files (AntiAlias)'); 
        sys.exit(1)
    # angulations     
    if str(ACQPdata["ACQ_grad_matrix"]) != str(ACQPdataRef["ACQ_grad_matrix"]):
        print ('ERROR: Parameter mismatch beteen Main and Reference files (Angulations)'); 
        sys.exit(1)    

# read input from keyboard
if not haveRef:
    xcor=0; ycor=0; zcor=0
    OK=False
    while not OK:
        dummy = raw_input("Enter Flow correction X [cm/s]: ")
        if dummy == '': dummy=0
        try: xcor = float(dummy); OK=True
        except: print ("Input Error")
    OK=False
    while not OK:
        dummy = raw_input("Enter Flow correction Y [cm/s]: ")
        if dummy == '': dummy=0    
        try: ycor = float(dummy); OK=True
        except: print ("Input Error")
    OK=False
    while not OK:
        dummy = raw_input("Enter Flow correction Z [cm/s]: ")
        if dummy == '': dummy=0    
        try: zcor = float(dummy); OK=True
        except: print ("Input Error")

        
#start
print ('Starting recon')    

#reshape FID data according to dimensions from method file
#"order="F" means Fortran style order as by BRUKER conventions
dim=METHODdata["PVM_EncMatrix"]
dim0 = dim[0]; dim0_mod_128 = dim0%128
if dim0_mod_128!=0: dim0=(int(dim0/128)+1)*128 # Bruker sets readout point to a multiple of 128
try: FIDrawdata_CPX = FIDrawdata_CPX.reshape(dim0,4,dim[1],dim[2], order="F")
except: print ('ERROR: k-space data reshape failed (dimension problem)'); sys.exit(1)
if dim0 != dim[0]: FIDrawdata_CPX = FIDrawdata_CPX[0:dim[0],:,:,:]

# reshape reference
if haveRef:
    dimRef=METHODdataRef["PVM_EncMatrix"]
    dim0 = dimRef[0]; dim0_mod_128 = dim0%128
    if dim0_mod_128!=0: dim0=(int(dim0/128)+1)*128 # Bruker sets readout point to a multiple of 128
    try: FIDrawdata_CPX_REF = FIDrawdata_CPX_REF.reshape(dim0,4,dimRef[1],dimRef[2], order="F")
    except: print ('ERROR: reference k-space data reshape failed (dimension problem)'); sys.exit(1)
    if dim0 != dimRef[0]: FIDrawdata_CPX_REF = FIDrawdata_CPX_REF[0:dimRef[0],:,:,:]

#reorder data (main)
FIDdata_tmp=np.empty(shape=(dim[0],4,dim[1],dim[2]),dtype=np.complex64)
FIDdata=np.empty(shape=(dim[0],4,dim[1],dim[2]),dtype=np.complex64)
order1=METHODdata["PVM_EncSteps1"]+dim[1]/2                             
for i in range(0,dim[1]): FIDdata_tmp[:,:,order1[i],:]=FIDrawdata_CPX[:,:,i,:]
FIDrawdata_CPX = 0 #free memory  
order2=METHODdata["PVM_EncSteps2"]+dim[2]/2                             
for i in range(0,dim[2]): FIDdata[:,:,:,order2[i]]=FIDdata_tmp[:,:,:,i]
FIDdata_tmp = 0 #free memory  
print('.', end='') #progress indicator

#reorder data (reference)
if haveRef:
    FIDdata_tmp=np.empty(shape=(dimRef[0],4,dimRef[1],dimRef[2]),dtype=np.complex64)
    FIDdataRef=np.empty(shape=(dimRef[0],4,dimRef[1],dimRef[2]),dtype=np.complex64)
    order1=METHODdataRef["PVM_EncSteps1"]+dimRef[1]/2                             
    for i in range(0,dimRef[1]): FIDdata_tmp[:,:,order1[i],:]=FIDrawdata_CPX_REF[:,:,i,:]
    FIDrawdata_CPX_REF = 0 #free memory  
    order2=METHODdataRef["PVM_EncSteps2"]+dimRef[2]/2                             
    for i in range(0,dimRef[2]): FIDdataRef[:,:,:,order2[i]]=FIDdata_tmp[:,:,:,i]
    FIDdata_tmp = 0 #free memory  
    print('.', end='') #progress indicator

#Hanning filter (main)
percentage = 5.
npoints_x = int(float(dim[0])*percentage/100.)
hanning_x = np.empty(shape=(dim[0]),dtype=np.float32)
x_ = np.linspace (1./(npoints_x-1.)*np.pi/2.,(1.-1./(npoints_x-1))*np.pi/2.,num=npoints_x)
hanning_x [0:npoints_x] = np.power(np.sin(x_),2)
hanning_x [npoints_x:hanning_x.shape[0]-npoints_x] = 1
x_ = x_[::-1] # reverse x_
hanning_x [hanning_x.shape[0]-npoints_x:hanning_x.shape[0]] = np.power(np.sin(x_),2)
#print (hanning_x)
FIDdata[:,:,:,:] *= hanning_x [:,None,None,None]
npoints_y = int(float(dim[1])*percentage/100.)
hanning_y = np.empty(shape=(dim[1]),dtype=np.float32)
y_ = np.linspace (1./(npoints_y-1.)*np.pi/2.,(1.-1./(npoints_y-1))*np.pi/2.,num=npoints_y)
hanning_y [0:npoints_y] = np.power(np.sin(y_),2)
hanning_y [npoints_y:hanning_y.shape[0]-npoints_y] = 1
y_ = y_[::-1] # reverse y_
hanning_y [hanning_y.shape[0]-npoints_y:hanning_y.shape[0]] = np.power(np.sin(y_),2)
#print (hanning_y)
FIDdata[:,:,:,:] *= hanning_y [None,None,:,None]
npoints_z = int(float(dim[2])*percentage/100.)
hanning_z = np.empty(shape=(dim[2]),dtype=np.float32)
z_ = np.linspace (1./(npoints_z-1.)*np.pi/2.,(1.-1./(npoints_z-1))*np.pi/2.,num=npoints_z)
hanning_z [0:npoints_z] = np.power(np.sin(z_),2)
hanning_z [npoints_z:hanning_z.shape[0]-npoints_z] = 1
z_ = z_[::-1] # reverse z_
hanning_z [hanning_z.shape[0]-npoints_z:hanning_z.shape[0]] = np.power(np.sin(z_),2)
#print (hanning_z)
FIDdata[:,:,:,:] *= hanning_z [None,None,None,:]
print('.', end='') #progress indicator

#for reference dimension < main imensions fill in zeros
if haveRef:
    if not np.array_equal (dimRef,dim):
        if dimRef[0]>dim[0]: print ('ERROR: in reference dimension 1'); sys.exit(1)
        if dimRef[1]>dim[1]: print ('ERROR: in reference dimension 2'); sys.exit(1) 
        if dimRef[2]>dim[2]: print ('ERROR: in reference dimension 3'); sys.exit(1)  
        FIDdataRef_=np.zeros(shape=FIDdata.shape,dtype=np.complex64)
        dim0start=int((dim[0]-dimRef[0])/2)
        dim1start=int((dim[1]-dimRef[1])/2)
        dim2start=int((dim[2]-dimRef[2])/2)
        FIDdataRef_[dim0start:dim0start+dimRef[0],:,dim1start:dim1start+dimRef[1],dim2start:dim2start+dimRef[2]] = \
            FIDdataRef[0:dimRef[0],:,0:dimRef[1],0:dimRef[2]]
        FIDdataRef=FIDdataRef_
        FIDdataRef_ = 0 #free memory
        dimRef=dim

#modified Hanning filter (reference)
if haveRef:
    percentage_full = 20. 
    #find real echo positions
    max_echo_pos_x = np.argmax(np.sum(np.abs(FIDdataRef[:,:,:,:]),axis=(1,2,3)))
    max_echo_pos_y = np.argmax(np.sum(np.abs(FIDdataRef[:,:,:,:]),axis=(0,1,3)))
    max_echo_pos_z = np.argmax(np.sum(np.abs(FIDdataRef[:,:,:,:]),axis=(0,1,2)))
    #print ('echo maximum position found at', max_echo_pos_x, max_echo_pos_y, max_echo_pos_z)
    npoints_x = int(float(dimRef[0])*percentage_full/100.)
    npoints_x = int(npoints_x/2.)*2+1 # make it odd
    hanning_x = np.zeros(shape=(dimRef[0]),dtype=np.float32)
    x_ = np.linspace (0,np.pi,num=npoints_x)
    hanning_x [int(dimRef[0]/2)-npoints_x/2:int(dimRef[0]/2)+npoints_x/2+1] = 1-np.power(np.cos(x_),4)
    #print (hanning_x)
    hanning_x = np.roll (hanning_x, max_echo_pos_x-np.argmax(hanning_x)) # set to real echo position
    FIDdataRef[:,:,:,:] *= hanning_x [:,None,None,None]
    npoints_y = int(float(dimRef[1])*percentage_full/100.)
    npoints_y = int(npoints_y/2.)*2+1 # make it odd
    hanning_y = np.zeros(shape=(dimRef[1]),dtype=np.float32)
    y_ = np.linspace (0,np.pi,num=npoints_y)
    hanning_y [int(dimRef[1]/2)-npoints_y/2:int(dimRef[1]/2)+npoints_y/2+1] = 1-np.power(np.cos(y_),4)
    #print (hanning_y)
    hanning_y = np.roll (hanning_y, max_echo_pos_y-np.argmax(hanning_y)) # set to real echo position
    FIDdataRef[:,:,:,:] *= hanning_y [None,None,:,None]
    npoints_z = int(float(dimRef[2])*percentage_full/100.)
    npoints_z = int(npoints_z/2.)*2+1 # make it odd
    hanning_z = np.zeros(shape=(dimRef[2]),dtype=np.float32)
    z_ = np.linspace (0,np.pi,num=npoints_z)
    hanning_z [int(dimRef[2]/2)-npoints_z/2:int(dimRef[2]/2)+npoints_z/2+1] = 1-np.power(np.cos(z_),4)
    #print (hanning_z)
    hanning_z = np.roll (hanning_z, max_echo_pos_z-np.argmax(hanning_z)) # set to real echo position
    FIDdataRef[:,:,:,:] *= hanning_z [None,None,None,:]
    print('.', end='') #progress indicator


# apply FOV offsets = (linear phase in k-space) main
PackArrPhase1Offset=METHODdata["PVM_SPackArrPhase1Offset"]
SPackArrSliceOffset=METHODdata["PVM_SPackArrSliceOffset"]
realFOV = METHODdata["PVM_Fov"]*METHODdata["PVM_AntiAlias"]
phase_step1 = +2.*np.pi*float(PackArrPhase1Offset)/float(realFOV[1])
phase_step2 = -2.*np.pi*float(SPackArrSliceOffset)/float(realFOV[2])
mag = np.abs(FIDdata[:,:,:,:]); ph = np.angle(FIDdata[:,:,:,:])
for i in range(0,FIDdata.shape[2]): ph[:,:,i,:] -= float(i-int(FIDdata.shape[2]/2))*phase_step1
for j in range(0,FIDdata.shape[3]): ph[:,:,:,j] -= float(j-int(FIDdata.shape[3]/2))*phase_step2
FIDdata [:,:,:,:] = mag * np.exp(1j*ph)
mag= 0; ph=0 # free memory 
print('.', end='') #progress indicator

# apply FOV offsets = (linear phase in k-space) reference
if haveRef:
    PackArrPhase1Offset=METHODdataRef["PVM_SPackArrPhase1Offset"]
    SPackArrSliceOffset=METHODdataRef["PVM_SPackArrSliceOffset"]
    realFOV = METHODdataRef["PVM_Fov"]*METHODdataRef["PVM_AntiAlias"]
    phase_step1 = +2.*np.pi*float(PackArrPhase1Offset)/float(realFOV[1])
    phase_step2 = -2.*np.pi*float(SPackArrSliceOffset)/float(realFOV[2])
    mag = np.abs(FIDdataRef[:,:,:,:]); ph = np.angle(FIDdataRef[:,:,:,:])
    for i in range(0,FIDdataRef.shape[2]): ph[:,:,i,:] -= float(i-int(FIDdataRef.shape[2]/2))*phase_step1
    for j in range(0,FIDdataRef.shape[3]): ph[:,:,:,j] -= float(j-int(FIDdataRef.shape[3]/2))*phase_step2
    FIDdataRef [:,:,:,:] = mag * np.exp(1j*ph)
    mag= 0; ph=0 # free memory
    print('.', end='') #progress indicator

#zero fill (main)
zero_fill=2.
SpatResol=METHODdata["PVM_SpatResol"]/zero_fill
FIDdata_ZF=np.zeros(shape=(int(dim[0]*zero_fill),4,int(dim[1]*zero_fill),
                           int(dim[2]*zero_fill)),dtype=np.complex64)
dim0start=int(dim[0]*(zero_fill-1)/2)
dim1start=int(dim[1]*(zero_fill-1)/2)
dim2start=int(dim[2]*(zero_fill-1)/2)
FIDdata_ZF[dim0start:dim0start+dim[0],:,dim1start:dim1start+dim[1],dim2start:dim2start+dim[2]] = \
    FIDdata[0:dim[0],:,0:dim[1],0:dim[2]]
FIDdata=FIDdata_ZF;
FIDdata_ZF = 0 #free memory
dim=dim*zero_fill; dim=dim.astype(int) 

#zero fill (reference)
if haveRef:
    FIDdata_ZF=np.zeros(shape=(int(dimRef[0]*zero_fill),4,int(dimRef[1]*zero_fill),
                               int(dimRef[2]*zero_fill)),dtype=np.complex64)
    dim0start=int(dimRef[0]*(zero_fill-1)/2)
    dim1start=int(dimRef[1]*(zero_fill-1)/2)
    dim2start=int(dimRef[2]*(zero_fill-1)/2)
    FIDdata_ZF[dim0start:dim0start+dimRef[0],:,dim1start:dim1start+dimRef[1],dim2start:dim2start+dimRef[2]] = \
        FIDdataRef[0:dimRef[0],:,0:dimRef[1],0:dimRef[2]]
    FIDdataRef=FIDdata_ZF;
    FIDdata_ZF = 0 #free memory
    dimRef=dimRef*zero_fill; dimRef=dimRef.astype(int) 

#FFT (main)
EchoPosition_raw=METHODdata["PVM_EchoPosition"]
EchoPosition_raw=50-(50-EchoPosition_raw)/zero_fill
EchoPosition=int(EchoPosition_raw/100.*dim[0])
IMGdata=FIDdata
FIDdata = 0 #free memory
IMGdata=np.roll(IMGdata, -EchoPosition, axis=(0))
IMGdata = np.fft.fftshift(IMGdata, axes=(2))
IMGdata = np.fft.fftshift(IMGdata, axes=(3)); print('.', end='') #progress indicator
for i in range(0,IMGdata.shape[1]):
    IMGdata [:,i,:,:] = FFT3D(IMGdata[:,i,:,:])
    print('.', end='') #progress indicator
IMGdata = np.fft.fftshift(IMGdata, axes=(0))
IMGdata = np.fft.fftshift(IMGdata, axes=(2))
IMGdata = np.fft.fftshift(IMGdata, axes=(3))
print('.', end='') #progress indicator

#FFT (reference)
if haveRef:
    EchoPosition_raw=METHODdataRef["PVM_EchoPosition"]
    EchoPosition_raw=50-(50-EchoPosition_raw)/zero_fill
    EchoPosition=int(EchoPosition_raw/100.*dimRef[0])
    IMGdataRef=FIDdataRef
    FIDdataRef = 0 #free memory
    IMGdataRef=np.roll(IMGdataRef, -EchoPosition, axis=(0))
    IMGdataRef = np.fft.fftshift(IMGdataRef, axes=(2))
    IMGdataRef = np.fft.fftshift(IMGdataRef, axes=(3)); print('.', end='') #progress indicator
    for i in range(0,IMGdataRef.shape[1]):
        IMGdataRef [:,i,:,:] = FFT3D(IMGdataRef[:,i,:,:])
        print('.', end='') #progress indicator
    IMGdataRef = np.fft.fftshift(IMGdataRef, axes=(0))
    IMGdataRef = np.fft.fftshift(IMGdataRef, axes=(2))
    IMGdataRef = np.fft.fftshift(IMGdataRef, axes=(3))
    print('.', end='') #progress indicator

        
#throw out antialiasing (main)
crop=METHODdata["PVM_AntiAlias"]
dim0start=int((dim[0]-dim[0]/crop[0])/2)
dim1start=int((dim[1]-dim[1]/crop[1])/2)
dim2start=int((dim[2]-dim[2]/crop[2])/2)
dim0end = int(dim0start+dim[0]/crop[0])
dim1end = int(dim1start+dim[1]/crop[1])
dim2end = int(dim2start+dim[2]/crop[2])
IMGdata = IMGdata[dim0start:dim0end,:,dim1start:dim1end,dim2start:dim2end]
dim[0] = dim0end-dim0start
dim[1] = dim1end-dim1start
dim[2] = dim2end-dim2start
print('.', end='') #progress indicator

#throw out antialiasing (reference)
if haveRef:
    crop=METHODdataRef["PVM_AntiAlias"]
    dim0start=int((dimRef[0]-dimRef[0]/crop[0])/2)
    dim1start=int((dimRef[1]-dimRef[1]/crop[1])/2)
    dim2start=int((dimRef[2]-dimRef[2]/crop[2])/2)
    dim0end = int(dim0start+dimRef[0]/crop[0])
    dim1end = int(dim1start+dimRef[1]/crop[1])
    dim2end = int(dim2start+dimRef[2]/crop[2])
    IMGdataRef = IMGdataRef[dim0start:dim0end,:,dim1start:dim1end,dim2start:dim2end]
    dimRef[0] = dim0end-dim0start
    dimRef[1] = dim1end-dim1start
    dimRef[2] = dim2end-dim2start
    print('.', end='') #progress indicator



    
    
#Phase correction 
if haveRef: # using reference, this is what actually does the trick
    # 1) Magnitude correction (with veeeerrry heavy filtering to avoid SNR degradation)
    temp = np.abs(IMGdataRef)
    base = np.average(temp, axis=(1))
    mag_cor = np.zeros (shape=IMGdataRef.shape)
    mag_cor[:,:,:,:] = temp[:,:,:,:]/base[:,None,:,:]
    temp = 0; base = 0  # free memory    
    mag_cor = zoom(mag_cor[:,:,:,:],[0.25,1,0.25,0.25],order=1) # downsample
    print('.', end='') #progress indicator
    mag_cor = median_filter(mag_cor, size = (7,1,7,7)) # median filter
    print('.', end='') #progress indicator
    s = 3; w = 9; t = (((w - 1)/2)-0.5)/s
    mag_cor[:,0,:,:] = gaussian_filter(mag_cor[:,0,:,:], sigma=s, truncate=t)
    mag_cor[:,1,:,:] = gaussian_filter(mag_cor[:,1,:,:], sigma=s, truncate=t)
    mag_cor[:,2,:,:] = gaussian_filter(mag_cor[:,2,:,:], sigma=s, truncate=t)
    mag_cor[:,3,:,:] = gaussian_filter(mag_cor[:,3,:,:], sigma=s, truncate=t)  
    print('.', end='') #progress indicator
    mag_cor = zoom(mag_cor[:,:,:,:],[4,1,4,4],order=1) # upsample    
    print('.', end='') #progress indicator
    # 2) Phase Correction - complex division is equivalent to phase subtraction
    temp = IMGdata/IMGdataRef; 
    ph = np.angle(temp [:,:,:,:])
    temp = 0 # free memory
    #Put things together again
    mag = np.abs(IMGdata [:,:,:,:])/mag_cor; 
    IMGdata [:,:,:,:] = mag * np.exp(1j*ph)
    print('.', end='') #progress indicator
    mag_cor = 0 # free memory
    IMGdataRef = 0 # free memory
else: # without reference enable manual 0 order correction
    venc=METHODdata["FlowRange"]
    mag = np.abs(IMGdata [:,0,:,:]); ph = np.angle(IMGdata [:,0,:,:])
    ph += -1.0*xcor/venc*np.pi/4. -1.0*ycor/venc*np.pi/4. -1.0*zcor/venc*np.pi/4.
    IMGdata [:,0,:,:] = mag * np.exp(1j*ph)
    print('.', end='') #progress indicator
    mag = np.abs(IMGdata [:,1,:,:]); ph = np.angle(IMGdata [:,1,:,:])
    ph += +1.0*xcor/venc*np.pi/4. +1.0*ycor/venc*np.pi/4. -1.0*zcor/venc*np.pi/4.
    IMGdata [:,1,:,:] = mag * np.exp(1j*ph)
    print('.', end='') #progress indicator
    mag = np.abs(IMGdata [:,2,:,:]); ph = np.angle(IMGdata [:,2,:,:])
    ph += +1.0*xcor/venc*np.pi/4. -1.0*ycor/venc*np.pi/4. +1.0*zcor/venc*np.pi/4.
    IMGdata [:,2,:,:] = mag * np.exp(1j*ph)
    print('.', end='') #progress indicator
    mag = np.abs(IMGdata [:,3,:,:]); ph = np.angle(IMGdata [:,3,:,:])
    ph += -1.0*xcor/venc*np.pi/4. +1.0*ycor/venc*np.pi/4. +1.0*zcor/venc*np.pi/4.
    IMGdata [:,3,:,:] = mag * np.exp(1j*ph)
    print('.', end='') #progress indicator
    
    
# Hadamard decoding
#
# from https://www.ncbi.nlm.nih.gov/pubmed/1790361
# Simultaneous acquisition of phase-contrast angiograms and stationary-tissue
# images with Hadamard encoding of flow-induced phase shifts.
# Dumoulin CL, Souza SP, Darrow RD, Pelc NJ, Adams WJ, Ash SA.
# J Magn Reson Imaging. 1991 Jul-Aug;1(4):399-404.
# DOI: 10.1002/jmri.1880010403
#
# Exp0 = -1 -1 -1
# Exp1 = +1 +1 -1
# Exp2 = +1 -1 +1
# Exp3 = -1 +1 +1
#
# Stationary =  Exp0 + Exp1 + Exp2 + Exp3
# Flow_X     = -Exp0 + Exp1 + Exp2 - Exp3
# Flow_Y     = -Exp0 + Exp1 - Exp2 + Exp3
# Flow_Z     = -Exp0 - Exp1 + Exp2 + Exp3
#
#
# to recon magnitude images take the complex difference (+-)
IMGdata_decoded_ABS = np.empty(shape=(dim[0],5,dim[1],dim[2]),dtype=np.float32)
IMGdata_decoded_ABS [:,0,:,:] = np.abs(    IMGdata [:,0,:,:] + IMGdata [:,1,:,:] + IMGdata [:,2,:,:] + IMGdata [:,3,:,:])
IMGdata_decoded_ABS [:,1,:,:] = np.abs( -1*IMGdata [:,0,:,:] + IMGdata [:,1,:,:] + IMGdata [:,2,:,:] - IMGdata [:,3,:,:])
IMGdata_decoded_ABS [:,2,:,:] = np.abs( -1*IMGdata [:,0,:,:] + IMGdata [:,1,:,:] - IMGdata [:,2,:,:] + IMGdata [:,3,:,:])
IMGdata_decoded_ABS [:,3,:,:] = np.abs( -1*IMGdata [:,0,:,:] - IMGdata [:,1,:,:] + IMGdata [:,2,:,:] + IMGdata [:,3,:,:])
IMGdata_decoded_ABS [:,4,:,:] = np.sqrt(np.power(IMGdata_decoded_ABS[:,1,:,:],2) +
                                        np.power(IMGdata_decoded_ABS[:,2,:,:],2) +
                                        np.power(IMGdata_decoded_ABS[:,3,:,:],2))
# to recon phase images take the complex division (*/)
IMGdata_decoded_PH = np.empty(shape=(dim[0],5,dim[1],dim[2]),dtype=np.float32)
IMGdata_decoded_PH [:,0,:,:] = -1.0*np.angle(   IMGdata [:,0,:,:] * IMGdata [:,1,:,:] * IMGdata [:,2,:,:] * IMGdata [:,3,:,:])
IMGdata_decoded_PH [:,1,:,:] = -1.0*np.angle( 1/IMGdata [:,0,:,:] * IMGdata [:,1,:,:] * IMGdata [:,2,:,:] / IMGdata [:,3,:,:])
IMGdata_decoded_PH [:,2,:,:] = -1.0*np.angle( 1/IMGdata [:,0,:,:] * IMGdata [:,1,:,:] / IMGdata [:,2,:,:] * IMGdata [:,3,:,:])
IMGdata_decoded_PH [:,3,:,:] = -1.0*np.angle( 1/IMGdata [:,0,:,:] / IMGdata [:,1,:,:] * IMGdata [:,2,:,:] * IMGdata [:,3,:,:])
IMGdata_decoded_PH [:,4,:,:] = np.sqrt(np.power(IMGdata_decoded_PH[:,1,:,:],2) +
                                        np.power(IMGdata_decoded_PH[:,2,:,:],2) +
                                        np.power(IMGdata_decoded_PH[:,3,:,:],2))

#permute dimensions
#worx for PVM_SPackArrSliceOrient=sagittal, PVM_SPackArrReadOrient="H_F"
#this way results are comparabled to ImageJ's BrukerOpener plugin
if METHODdata["PVM_SPackArrSliceOrient"] == "sagittal" and METHODdata["PVM_SPackArrReadOrient"] == "H_F":
    SpatResol_perm = np.empty(shape=(3))
    SpatResol_perm[0]=SpatResol[1]
    SpatResol_perm[1]=SpatResol[0]
    SpatResol_perm[2]=SpatResol[2]
    IMGdata_decoded_ABS = np.transpose (IMGdata_decoded_ABS, axes=(2,1,0,3))
    IMGdata_decoded_ABS = np.rot90(IMGdata_decoded_ABS, k=2, axes=(0,3)) # k=2 is a 180 degree rotation
    IMGdata_decoded_PH = np.transpose (IMGdata_decoded_PH, axes=(2,1,0,3))
    IMGdata_decoded_PH = np.rot90(IMGdata_decoded_PH, k=2, axes=(0,3)) # k=2 is a 180 degree rotation
    #The following is to get the vector vizualization right in Paraview
    dummy = np.zeros(shape=(IMGdata_decoded_PH.shape[0],IMGdata_decoded_PH.shape[2],IMGdata_decoded_PH.shape[3]),dtype=np.float32)
    dummy [:] = IMGdata_decoded_PH[:,1,:,:] # save flow X component
    IMGdata_decoded_PH[:,1,:,:] = IMGdata_decoded_PH[:,2,:,:] # X = Y
    IMGdata_decoded_PH[:,2,:,:] = dummy # Y = X
    IMGdata_decoded_PH[:,1,:,:] *= -1 # X invert sign    
elif METHODdata["PVM_SPackArrSliceOrient"] == "sagittal" and METHODdata["PVM_SPackArrReadOrient"] == "A_P":
    SpatResol_perm = np.empty(shape=(3))
    SpatResol_perm[0]=SpatResol[0]
    SpatResol_perm[1]=SpatResol[1]
    SpatResol_perm[2]=SpatResol[2]
    IMGdata_decoded_ABS = IMGdata_decoded_ABS[::-1,:,:,:] # reverse x  
    IMGdata_decoded_ABS = IMGdata_decoded_ABS[:,:,:,::-1] # reverse y      
    IMGdata_decoded_PH = IMGdata_decoded_PH[::-1,:,:,:] # reverse x  
    IMGdata_decoded_PH = IMGdata_decoded_PH[:,:,:,::-1] # reverse y    
    #The following is to get the vector vizualization right in Paraview
    IMGdata_decoded_PH[:,1,:,:] *= -1 # X invert sign
elif METHODdata["PVM_SPackArrSliceOrient"] == "coronal" and METHODdata["PVM_SPackArrReadOrient"] == "H_F":
    SpatResol_perm = np.empty(shape=(3))
    SpatResol_perm[0]=SpatResol[1]
    SpatResol_perm[1]=SpatResol[0]
    SpatResol_perm[2]=SpatResol[2]
    IMGdata_decoded_ABS = np.transpose (IMGdata_decoded_ABS, axes=(2,1,0,3))
    IMGdata_decoded_ABS = IMGdata_decoded_ABS[::-1,:,:,:] # flip axis
    IMGdata_decoded_ABS = IMGdata_decoded_ABS[:,:,:,::-1] # flip axis    
    IMGdata_decoded_PH = np.transpose (IMGdata_decoded_PH, axes=(2,1,0,3))
    IMGdata_decoded_PH = IMGdata_decoded_PH[::-1,:,:,:] # flip axis
    IMGdata_decoded_PH = IMGdata_decoded_PH[:,:,:,::-1] # flip axis
    #The following is to get the vector vizualization right in Paraview
    dummy = np.zeros(shape=(IMGdata_decoded_PH.shape[0],IMGdata_decoded_PH.shape[2],IMGdata_decoded_PH.shape[3]),dtype=np.float32)
    dummy [:] = IMGdata_decoded_PH[:,1,:,:] # save flow X component
    IMGdata_decoded_PH[:,1,:,:] = IMGdata_decoded_PH[:,2,:,:] # X = Y
    IMGdata_decoded_PH[:,2,:,:] = dummy # Y = X
    IMGdata_decoded_PH[:,1,:,:] *= -1 # X invert sign
elif METHODdata["PVM_SPackArrSliceOrient"] == "axial" and METHODdata["PVM_SPackArrReadOrient"] == "A_P":
    IMGdata_decoded_ABS = np.rot90(IMGdata_decoded_ABS, k=1, axes=(0,2)) # rotate (axial A_P)
    IMGdata_decoded_PH = np.rot90(IMGdata_decoded_PH, k=1, axes=(0,2)) # rotate (axial A_P)
    SpatResol_perm = np.empty(shape=(3))    
    SpatResol_perm[0] = SpatResol[1]
    SpatResol_perm[1] = SpatResol[0]
    SpatResol_perm[2] = SpatResol[2]
    #The following is to get the vector vizualization right in Paraview
    dummy = np.zeros(shape=(IMGdata_decoded_PH.shape[0],IMGdata_decoded_PH.shape[2],IMGdata_decoded_PH.shape[3]),dtype=np.float32)
    dummy [:] = IMGdata_decoded_PH[:,1,:,:] # save flow X component
    IMGdata_decoded_PH[:,1,:,:] = IMGdata_decoded_PH[:,2,:,:] # X = Y
    IMGdata_decoded_PH[:,2,:,:] = dummy # Y = X
    IMGdata_decoded_PH[:,1,:,:] *= -1 # X invert sign 
    IMGdata_decoded_PH[:,3,:,:] *= -1 # Z invert sign    
else:
    SpatResol_perm = SpatResol
    print ('Warning: unknown Orientation',METHODdata["PVM_SPackArrSliceOrient"], 
            METHODdata["PVM_SPackArrReadOrient"]);
    print ('         resulting images may be rotated incorrectly'); 
print('.', end='') #progress indicator

#find noise mask threshold from histogram
#image_number = 4 # 0 is static, 4 is flow
#n_points=IMGdata_decoded_ABS.shape[0]*IMGdata_decoded_ABS.shape[2]*IMGdata_decoded_ABS.shape[3]
#steps=int(n_points/1000); start=0; fin=np.max(IMGdata_decoded_ABS [:,image_number,:,:])
#xbins =  np.linspace(start,fin,steps)
#ybins, binedges = np.histogram(IMGdata_decoded_ABS [:,image_number,:,:], bins=xbins)
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
#mask =  IMGdata_decoded_ABS [:,image_number,:,:] > threshold  
#enable the following view histogram plot
#print ('\nThreshold = %.2e' % threshold)
#import pylab; pylab.plot(xbins,ybins, linewidth=1.2); pylab.draw();
#pylab.show(block=False); os.system("pause"); pylab.close(); 

# use noise in all 8 corners to establish threshold
image_number = 4 # 0 is static, 4 is flow
N=10 # use 10% at the corners of the FOV
std_fac = 6 # how many standard deviations to add
tresh=np.empty(shape=8,dtype=np.float)
avg=np.empty(shape=8,dtype=np.float)
std=np.empty(shape=8,dtype=np.float)
xstart=0; xend=int(IMGdata_decoded_ABS.shape[0]/N)
ystart=0; yend=int(IMGdata_decoded_ABS.shape[2]/N)
zstart=0; zend=int(IMGdata_decoded_ABS.shape[3]/N)
arr=IMGdata_decoded_ABS[xstart:xend,image_number,ystart:yend,zstart:zend]
avg[0]=np.mean(arr)
std[0]=np.std(arr)
tresh[0]=avg[0] + std_fac*std[0]
xstart=int(IMGdata_decoded_ABS.shape[0]-IMGdata_decoded_ABS.shape[0]/N); xend=IMGdata_decoded_ABS.shape[0]
ystart=0; yend=int(IMGdata_decoded_ABS.shape[2]/N)
zstart=0; zend=int(IMGdata_decoded_ABS.shape[3]/N)
arr=IMGdata_decoded_ABS[xstart:xend,image_number,ystart:yend,zstart:zend]
avg[1]=np.mean(arr)
std[1]=np.std(arr)
tresh[1]=avg[1] + std_fac*std[1]
xstart=0; xend=int(IMGdata_decoded_ABS.shape[0]/N)
ystart=int(IMGdata_decoded_ABS.shape[2]-IMGdata_decoded_ABS.shape[2]/N); yend=IMGdata_decoded_ABS.shape[2]
zstart=0; zend=int(IMGdata_decoded_ABS.shape[3]/N)
arr=IMGdata_decoded_ABS[xstart:xend,image_number,ystart:yend,zstart:zend]
avg[2]=np.mean(arr)
std[2]=np.std(arr)
tresh[2]=avg[2] + std_fac*std[2]
xstart=int(IMGdata_decoded_ABS.shape[0]-IMGdata_decoded_ABS.shape[0]/N); xend=IMGdata_decoded_ABS.shape[0]
ystart=int(IMGdata_decoded_ABS.shape[2]-IMGdata_decoded_ABS.shape[2]/N); yend=IMGdata_decoded_ABS.shape[2]
zstart=0; zend=int(IMGdata_decoded_ABS.shape[3]/N)
arr=IMGdata_decoded_ABS[xstart:xend,image_number,ystart:yend,zstart:zend]
avg[3]=np.mean(arr)
std[3]=np.std(arr)
tresh[3]=avg[3] + std_fac*std[3]
xstart=0; xend=int(IMGdata_decoded_ABS.shape[0]/N)
ystart=0; yend=int(IMGdata_decoded_ABS.shape[2]/N)
zstart=int(IMGdata_decoded_ABS.shape[3]-IMGdata_decoded_ABS.shape[2]/N); zend=IMGdata_decoded_ABS.shape[3]
arr=IMGdata_decoded_ABS[xstart:xend,image_number,ystart:yend,zstart:zend]
avg[4]=np.mean(arr)
std[4]=np.std(arr)
tresh[4]=avg[4] + std_fac*std[4]
xstart=int(IMGdata_decoded_ABS.shape[0]-IMGdata_decoded_ABS.shape[0]/N); xend=IMGdata_decoded_ABS.shape[0]
ystart=0; yend=int(IMGdata_decoded_ABS.shape[2]/N)
zstart=int(IMGdata_decoded_ABS.shape[3]-IMGdata_decoded_ABS.shape[3]/N); zend=IMGdata_decoded_ABS.shape[3]
arr=IMGdata_decoded_ABS[xstart:xend,image_number,ystart:yend,zstart:zend]
avg[5]=np.mean(arr)
std[5]=np.std(arr)
tresh[5]=avg[5] + std_fac*std[5]
xstart=0; xend=int(IMGdata_decoded_ABS.shape[0]/N)
ystart=int(IMGdata_decoded_ABS.shape[2]-IMGdata_decoded_ABS.shape[2]/N); yend=IMGdata_decoded_ABS.shape[2]
zstart=int(IMGdata_decoded_ABS.shape[3]-IMGdata_decoded_ABS.shape[3]/N); zend=IMGdata_decoded_ABS.shape[3]
arr=IMGdata_decoded_ABS[xstart:xend,image_number,ystart:yend,zstart:zend]
avg[6]=np.mean(arr)
std[6]=np.std(arr)
tresh[6]=avg[6] + std_fac*std[6]
xstart=int(IMGdata_decoded_ABS.shape[0]-IMGdata_decoded_ABS.shape[0]/N); xend=IMGdata_decoded_ABS.shape[0]
ystart=int(IMGdata_decoded_ABS.shape[2]-IMGdata_decoded_ABS.shape[2]/N); yend=IMGdata_decoded_ABS.shape[2]
zstart=int(IMGdata_decoded_ABS.shape[3]-IMGdata_decoded_ABS.shape[3]/N); zend=IMGdata_decoded_ABS.shape[3]
arr=IMGdata_decoded_ABS[xstart:xend,image_number,ystart:yend,zstart:zend]
avg[7]=np.mean(arr)
std[7]=np.std(arr)
tresh[7]=avg[7] + std_fac*std[7]
threshold=np.min(tresh)
mask =  IMGdata_decoded_ABS [:,image_number,:,:] > threshold


#transform to int
ReceiverGain = ACQPdata["RG"] # RG is a simple attenuation FACTOR, NOT in dezibel (dB) unit !!!
n_Averages = METHODdata["PVM_NAverages"]
IMGdata_decoded_ABS /= ReceiverGain*n_Averages; 
max_ABS = np.amax(IMGdata_decoded_ABS);
IMGdata_decoded_ABS *= 32767./max_ABS
IMGdata_decoded_ABS = IMGdata_decoded_ABS.astype(np.int16)
print('.', end='') #progress indicator
IMGdata_decoded_PH [:,0,:,:] *= mask*32767./np.pi
IMGdata_decoded_PH [:,1,:,:] *= mask*32767./np.pi
IMGdata_decoded_PH [:,2,:,:] *= mask*32767./np.pi
IMGdata_decoded_PH [:,3,:,:] *= mask*32767./np.pi
IMGdata_decoded_PH [:,4,:,:] *= mask*32767./(np.sqrt(np.square(np.pi)*3.)) # 5.44
IMGdata_decoded_PH = IMGdata_decoded_PH.astype(np.int16)
# set max/min in 0,0,0/1,1,1 corners (this save the venc parameterin the image)
IMGdata_decoded_PH [0,:,0,0] = 32767
IMGdata_decoded_PH [1,0:3,1,1] = -32767  
print('.', end='') #progress indicator

#createNIFTI's
venc=METHODdata["FlowRange"]
aff = np.eye(4)
aff[0,0] = SpatResol_perm[0]*1000; aff[0,3] = -(IMGdata_decoded_ABS.shape[0]/2)*aff[0,0]
aff[1,1] = SpatResol_perm[1]*1000; aff[1,3] = -(IMGdata_decoded_ABS.shape[2]/2)*aff[1,1]
aff[2,2] = SpatResol_perm[2]*1000; aff[2,3] = -(IMGdata_decoded_ABS.shape[3]/2)*aff[2,2]
#write Magnitude static
NIFTIimg = nib.Nifti1Image(IMGdata_decoded_ABS[:,0,:,:], aff)
NIFTIimg.header['sform_code']=1
NIFTIimg.header['qform_code']=1
NIFTIimg.header.set_slope_inter(max_ABS/32767.,0)
try: nib.save(NIFTIimg, os.path.join(os.path.dirname(FIDfile),OrigFilename+'_MAGNT_Static.nii.gz'))
except: print ('\nERROR:  problem while writing results'); sys.exit(1)
print('.', end='') #progress indicator
#write Magnitude flow
NIFTIimg = nib.Nifti1Image(IMGdata_decoded_ABS[:,4,:,:], aff)
NIFTIimg.header['sform_code']=1
NIFTIimg.header['qform_code']=1
NIFTIimg.header.set_slope_inter(max_ABS/32767.,0)
try: nib.save(NIFTIimg, os.path.join(os.path.dirname(FIDfile),OrigFilename+'_MAGNT_Flow.nii.gz'))
except: print ('\nERROR:  problem while writing results'); sys.exit(1)
print('.', end='') #progress indicator
#write Phase flow X
NIFTIimg = nib.Nifti1Image(IMGdata_decoded_PH[:,1,:,:], aff)
NIFTIimg.header['sform_code']=1
NIFTIimg.header['qform_code']=1
NIFTIimg.header.set_slope_inter(venc/32767.,0)
try: nib.save(NIFTIimg, os.path.join(os.path.dirname(FIDfile),OrigFilename+'_PHASE_Flow_X.nii.gz'))
except: print ('\nERROR:  problem while writing results'); sys.exit(1)
print('.', end='') #progress indicator
#write Phase flow Y
NIFTIimg = nib.Nifti1Image(IMGdata_decoded_PH[:,2,:,:], aff)
NIFTIimg.header['sform_code']=1
NIFTIimg.header['qform_code']=1
NIFTIimg.header.set_slope_inter(venc/32767.,0)
try: nib.save(NIFTIimg, os.path.join(os.path.dirname(FIDfile),OrigFilename+'_PHASE_Flow_Y.nii.gz'))
except: print ('\nERROR:  problem while writing results'); sys.exit(1)
print('.', end='') #progress indicator
#write Phase flow Z
NIFTIimg = nib.Nifti1Image(IMGdata_decoded_PH[:,3,:,:], aff)
NIFTIimg.header['sform_code']=1
NIFTIimg.header['qform_code']=1
NIFTIimg.header.set_slope_inter(venc/32767.,0)
try: nib.save(NIFTIimg, os.path.join(os.path.dirname(FIDfile),OrigFilename+'_PHASE_Flow_Z.nii.gz'))
except: print ('\nERROR:  problem while writing results'); sys.exit(1)
print('.', end='') #progress indicator
print('.', end='') #progress indicator
#write Phase flow ALL
NIFTIimg = nib.Nifti1Image(IMGdata_decoded_PH[:,4,:,:], aff)
NIFTIimg.header['sform_code']=1
NIFTIimg.header['qform_code']=1
NIFTIimg.header.set_slope_inter(np.sqrt(3)*venc/32767.,0)
try: nib.save(NIFTIimg, os.path.join(os.path.dirname(FIDfile),OrigFilename+'_PHASE_Flow.nii.gz'))
except: print ('\nERROR:  problem while writing results'); sys.exit(1)
print('.', end='') #progress indicator
#writesucess    
print ('\nSuccessfully written output files '+OrigFilename+'_*.nii.gz') 

#write flowvolume results
if METHODdata["PVM_SPackArrSliceOrient"] == "sagittal" and METHODdata["PVM_SPackArrReadOrient"] == "H_F":
    with open(os.path.join(os.path.dirname(FIDfile),OrigFilename+'_FlowVolumes.txt'), "w") as text_file:
        text_file.write("Flow Volumes per slice (X):\n")
        for i in range(0,IMGdata_decoded_PH.shape[2]): # in our data shape[2] is the main flow direction
            flowvol = np.sum(IMGdata_decoded_PH[:,1,i,:])
            flowvol *= 10.*venc/32767. # venc is in cm/s, multiply by 10. to get this in mm/s
            flowvol *= SpatResol_perm[0]*SpatResol_perm[2] # multiply with inplane spatial resolution, result is in mm^3/s
            flowvol /= 1000. # convert mm^3/s ---> ml/s
            text_file.write("Slice %d:\t%0.2f\tml/s\n" % (i, flowvol))
        text_file.write("\n")
elif METHODdata["PVM_SPackArrSliceOrient"] == "sagittal" and METHODdata["PVM_SPackArrReadOrient"] == "A_P":
    with open(os.path.join(os.path.dirname(FIDfile),OrigFilename+'_FlowVolumes.txt'), "w") as text_file:
        text_file.write("Flow Volumes per slice (X):\n")
        for i in range(0,IMGdata_decoded_PH.shape[2]): # in our data shape[2] is the main flow direction
            flowvol = np.sum(IMGdata_decoded_PH[:,2,i,:])
            flowvol *= 10.*venc/32767. # venc is in cm/s, multiply by 10. to get this in mm/s
            flowvol *= SpatResol_perm[0]*SpatResol_perm[2] # multiply with inplane spatial resolution, result is in mm^3/s
            flowvol /= 1000. # convert mm^3/s ---> ml/s
            text_file.write("Slice %d:\t%0.2f\tml/s\n" % (i, flowvol))
        text_file.write("\n")               
elif METHODdata["PVM_SPackArrSliceOrient"] == "coronal" and METHODdata["PVM_SPackArrReadOrient"] == "H_F":
    with open(os.path.join(os.path.dirname(FIDfile),OrigFilename+'_FlowVolumes.txt'), "w") as text_file:    
        text_file.write("Flow Volumes per slice (X):\n")
        for i in range(0,IMGdata_decoded_PH.shape[2]): # in our data shape[2] is the main flow direction
            flowvol = np.sum(IMGdata_decoded_PH[:,1,i,:])
            flowvol *= 10.*venc/32767. # venc is in cm/s, multiply by 10. to get this in mm/s
            flowvol *= SpatResol_perm[0]*SpatResol_perm[2] # multiply with inplane spatial resolution, result is in mm^3/s
            flowvol /= 1000. # convert mm^3/s ---> ml/s
            text_file.write("Slice %d:\t%0.2f\tml/s\n" % (i, flowvol))
        text_file.write("\n")        
elif METHODdata["PVM_SPackArrSliceOrient"] == "axial" and METHODdata["PVM_SPackArrReadOrient"] == "A_P":
    with open(os.path.join(os.path.dirname(FIDfile),OrigFilename+'_FlowVolumes.txt'), "w") as text_file:
        text_file.write("Flow Volumes per slice (X):\n")
        for i in range(0,IMGdata_decoded_PH.shape[3]): # in our data shape[2] is the main flow direction
            flowvol = np.sum(IMGdata_decoded_PH[:,3,:,i])
            flowvol *= 10.*venc/32767. # venc is in cm/s, multiply by 10. to get this in mm/s
            flowvol *= SpatResol_perm[0]*SpatResol_perm[1] # multiply with inplane spatial resolution, result is in mm^3/s
            flowvol /= 1000. # convert mm^3/s ---> ml/s
            text_file.write("Slice %d:\t%0.2f\tml/s\n" % (i, flowvol))
        text_file.write("\n")
else:
    try: os.remove (os.path.join(os.path.dirname(FIDfile),OrigFilename+'_FlowVolumes.txt'))
    except: pass #silent
    print ('Warning: textfile with flow values not written (unknown Orientation)');   
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