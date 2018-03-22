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
#    - matplotlib.pylab (optional)
#

from __future__ import print_function
try: import win32gui, win32console
except: pass #silent
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

# read input from keyboard
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

#reorder data
FIDdata_tmp=np.empty(shape=(dim[0],4,dim[1],dim[2]),dtype=np.complex64)
FIDdata=np.empty(shape=(dim[0],4,dim[1],dim[2]),dtype=np.complex64)
order1=METHODdata["PVM_EncSteps1"]+dim[1]/2                             
for i in range(0,dim[1]): FIDdata_tmp[:,:,order1[i],:]=FIDrawdata_CPX[:,:,i,:]
FIDrawdata_CPX = 0 #free memory  
order2=METHODdata["PVM_EncSteps2"]+dim[2]/2                             
for i in range(0,dim[2]): FIDdata[:,:,:,order2[i]]=FIDdata_tmp[:,:,:,i]
FIDdata_tmp = 0 #free memory  
print('.', end='') #progress indicator

#Hanning filter
percentage = 5.
npoints_x = int(float(dim[0])*percentage/100.)
hanning_x = np.empty(shape=(dim[0]),dtype=np.float32)
x_ = np.linspace (1./(npoints_x-1.)*np.pi/2.,(1.-1./(npoints_x-1))*np.pi/2.,num=npoints_x)
hanning_x [0:npoints_x] = np.power(np.sin(x_),2)
hanning_x [npoints_x:hanning_x.shape[0]-npoints_x] = 1
x_ = x_[::-1] # reverse x_
hanning_x [hanning_x.shape[0]-npoints_x:hanning_x.shape[0]] = np.power(np.sin(x_),2)
#print (hanning_x)
FIDdata[:,:,:] *= hanning_x [:,None,None,None]
npoints_y = int(float(dim[1])*percentage/100.)
hanning_y = np.empty(shape=(dim[1]),dtype=np.float32)
y_ = np.linspace (1./(npoints_y-1.)*np.pi/2.,(1.-1./(npoints_y-1))*np.pi/2.,num=npoints_y)
hanning_y [0:npoints_y] = np.power(np.sin(y_),2)
hanning_y [npoints_y:hanning_y.shape[0]-npoints_y] = 1
y_ = y_[::-1] # reverse y_
hanning_y [hanning_y.shape[0]-npoints_y:hanning_y.shape[0]] = np.power(np.sin(y_),2)
#print (hanning_x)
FIDdata[:,:,:] *= hanning_y [None,None,:,None]
npoints_z = int(float(dim[2])*percentage/100.)
hanning_z = np.empty(shape=(dim[2]),dtype=np.float32)
z_ = np.linspace (1./(npoints_z-1.)*np.pi/2.,(1.-1./(npoints_z-1))*np.pi/2.,num=npoints_z)
hanning_z [0:npoints_z] = np.power(np.sin(z_),2)
hanning_z [npoints_z:hanning_z.shape[0]-npoints_z] = 1
z_ = z_[::-1] # reverse z_
hanning_z [hanning_z.shape[0]-npoints_z:hanning_z.shape[0]] = np.power(np.sin(z_),2)
#print (hanning_x)
FIDdata[:,:,:] *= hanning_z [None,None,None,:]
print('.', end='') #progress indicator

# apply FOV offsets = (linear phase in k-space)
PackArrPhase1Offset=METHODdata["PVM_SPackArrPhase1Offset"]
SPackArrSliceOffset=METHODdata["PVM_SPackArrSliceOffset"]
realFOV = METHODdata["PVM_Fov"]*METHODdata["PVM_AntiAlias"]
phase_step1 = +2.*np.pi*float(PackArrPhase1Offset)/float(realFOV[1])
phase_step2 = -2.*np.pi*float(SPackArrSliceOffset)/float(realFOV[2])
mag = np.abs(FIDdata[:,:,:,:]); ph = np.angle(FIDdata[:,:,:,:])
for i in range(0,FIDdata.shape[2]): ph[:,:,i,:] -= float(i-int(FIDdata.shape[2]/2))*phase_step1
for j in range(0,FIDdata.shape[3]): ph[:,:,:,j] -= float(j-int(FIDdata.shape[3]/2))*phase_step2
FIDdata [:,:,:,:] = mag * np.exp(1j*ph)
print('.', end='') #progress indicator

#zero fill
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

#FFT (individually by axis, to save memory)
EchoPosition_raw=METHODdata["PVM_EchoPosition"]
EchoPosition_raw=50-(50-EchoPosition_raw)/zero_fill
EchoPosition=int(EchoPosition_raw/100.*dim[0])
IMGdata=FIDdata
FIDdata = 0 #free memory
Memory_OK=True
try: dummy=np.empty(IMGdata.shape, dtype=np.complex128); dummy=0 # try to allocate more memory
except: Memory_OK=False
if Memory_OK:
    IMGdata=np.roll(IMGdata, -EchoPosition, axis=(0))
    IMGdata = np.fft.fftshift(IMGdata, axes=(2))
    IMGdata = np.fft.fftshift(IMGdata, axes=(3)); print('.', end='') #progress indicator
    IMGdata = np.fft.fft(IMGdata, axis=0); print('.', end='') #progress indicator
    IMGdata = np.fft.fft(IMGdata, axis=2); print('.', end='') #progress indicator
    IMGdata = np.fft.fft(IMGdata, axis=3); print('.', end='') #progress indicator
    IMGdata = np.fft.fftshift(IMGdata, axes=(0))
    IMGdata = np.fft.fftshift(IMGdata, axes=(2))
    IMGdata = np.fft.fftshift(IMGdata, axes=(3))
else: # for large datasets use loops
    for k in range(0,IMGdata.shape[2]):
        if k%30==0: print(',', end='') #progress indicator
        IMGdata[:,:,k,:] = np.roll(IMGdata[:,:,k,:], -EchoPosition, axis=(0))
        IMGdata[:,:,k,:] = np.fft.fft(IMGdata[:,:,k,:], axis=(0))
        IMGdata[:,:,k,:] = np.fft.fftshift(IMGdata[:,:,k,:], axes=(0))
    for i in range(0,IMGdata.shape[0]):
        if i%30==0: print(',', end='') #progress indicator
        IMGdata[i,:,:,:] = np.fft.fftshift(IMGdata[i,:,:,:], axes=(1))
        IMGdata[i,:,:,:] = np.fft.fft(IMGdata[i,:,:,:], axis=(1))
        IMGdata[i,:,:,:] = np.fft.fftshift(IMGdata[i,:,:,:], axes=(1))
    for i in range(0,IMGdata.shape[0]):
        if i%30==0: print(',', end='') #progress indicator
        IMGdata[i,:,:,:] = np.fft.fftshift(IMGdata[i,:,:,:], axes=(2))
        IMGdata[i,:,:,:] = np.fft.fft(IMGdata[i,:,:,:], axis=(2))
        IMGdata[i,:,:,:] = np.fft.fftshift(IMGdata[i,:,:,:], axes=(2))   
print('.', end='') #progress indicator
        
#throw out antialiasing
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

#Phase offset correction
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
elif METHODdata["PVM_SPackArrSliceOrient"] == "sagittal" and METHODdata["PVM_SPackArrReadOrient"] == "A_P":
    SpatResol_perm = np.empty(shape=(3))
    SpatResol_perm[0]=SpatResol[0]
    SpatResol_perm[1]=SpatResol[1]
    SpatResol_perm[2]=SpatResol[2]
    IMGdata_decoded_ABS = IMGdata_decoded_ABS[::-1,:,:,:] # reverse x  
    IMGdata_decoded_ABS = IMGdata_decoded_ABS[:,:,:,::-1] # reverse y      
    IMGdata_decoded_PH = IMGdata_decoded_PH[::-1,:,:,:] # reverse x  
    IMGdata_decoded_PH = IMGdata_decoded_PH[:,:,:,::-1] # reverse y          
elif METHODdata["PVM_SPackArrSliceOrient"] == "coronal" and METHODdata["PVM_SPackArrReadOrient"] == "H_F":
    SpatResol_perm = np.empty(shape=(3))
    SpatResol_perm[0]=SpatResol[1]
    SpatResol_perm[1]=SpatResol[0]
    SpatResol_perm[2]=SpatResol[2]
    IMGdata_decoded_ABS = np.transpose (IMGdata_decoded_ABS, axes=(2,1,0,3))
    IMGdata_decoded_ABS = IMGdata_decoded_ABS[::-1,:,:,:] # flip axis
    IMGdata_decoded_PH = np.transpose (IMGdata_decoded_PH, axes=(2,1,0,3))
    IMGdata_decoded_PH = IMGdata_decoded_PH[::-1,:,:,:] # flip axis
elif METHODdata["PVM_SPackArrSliceOrient"] == "axial" and METHODdata["PVM_SPackArrReadOrient"] == "A_P":
    IMGdata_decoded_ABS = np.rot90(IMGdata_decoded_ABS, k=1, axes=(0,2)) # rotate (axial A_P)
    IMGdata_decoded_PH = np.rot90(IMGdata_decoded_PH, k=1, axes=(0,2)) # rotate (axial A_P)
    SpatResol_perm = np.empty(shape=(3))    
    SpatResol_perm[0] = SpatResol[1]
    SpatResol_perm[1] = SpatResol[0]
    SpatResol_perm[2] = SpatResol[2]
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
IMGdata_decoded_ABS /= ReceiverGain; 
max_ABS = np.amax(IMGdata_decoded_ABS);
IMGdata_decoded_ABS *= 32767./max_ABS
IMGdata_decoded_ABS = IMGdata_decoded_ABS.astype(np.int16)
print('.', end='') #progress indicator
IMGdata_decoded_PH [:,0,:,:] *= mask*32767./np.pi
IMGdata_decoded_PH [:,1,:,:] *= mask*32767./np.pi
IMGdata_decoded_PH [:,2,:,:] *= mask*32767./np.pi
IMGdata_decoded_PH [:,3,:,:] *= mask*32767./np.pi
IMGdata_decoded_PH [:,4,:,:] *= mask*32767./(np.sqrt(np.power(np.pi,2)*3)) # 5.44
IMGdata_decoded_PH = IMGdata_decoded_PH.astype(np.int16)
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
