#
# reads Bruker MR data (Paravision v5.1)
# reconstructs images from raw acquisition data (FID files)
# this version is for the FLASH 3D method only
#
# ----- VERSION HISTORY -----
#
# Version 0.1a - 21, November 2017
#       - this a customized version of the generic BrukerOfflineReco_Flash3D.py
#       - diff1: does not recon phase imges (due to memory limitations)
#       - diff2: does output a noise masked magnitude image
#       - diff3: does some quantitative measurements(undisclosed for IP reasons) 
#       - if you don't know whta this is about this version is probably NOT for you
#       - in this case probably better get the generic version BrukerOfflineReco_Flash3D.py
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

#read method file
METHODfile=os.path.dirname(FIDfile)+slash+'method'
METHODdata=ReadParamFile(METHODfile)

#check for not implemented stuff
if  not(METHODdata["Method"] == "FLASH" or METHODdata["Method"] == "FISP" or METHODdata["Method"] =="GEFC") or METHODdata["PVM_SpatDimEnum"] != "3D":
    print ('ERROR: Recon only implemented for FLASH/FISP 3D method'); 
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
    print ('ERROR: Recon for parallel acquisition not implemented'); 
    sys.exit(1)

#start
print ('Starting recon')    

#reshape FID data according to dimensions from method file
#"order="F" means Fortran style order as by BRUKER conventions
dim=METHODdata["PVM_EncMatrix"]
dim0 = dim[0]; dim0_mod_128 = dim0%128
if dim0_mod_128!=0: dim0=(int(dim0/128)+1)*128 # Bruker sets readout point to a multiple of 128
try: FIDrawdata_CPX = FIDrawdata_CPX.reshape(dim0,dim[1],dim[2], order="F")
except: print ('ERROR: k-space data reshape failed (dimension problem)'); sys.exit(1)
if dim0 != dim[0]: FIDrawdata_CPX = FIDrawdata_CPX[0:dim[0],:,:]

#reorder data
FIDdata_tmp=np.empty(shape=(dim[0],dim[1],dim[2]),dtype=np.complex64)
FIDdata=np.empty(shape=(dim[0],dim[1],dim[2]),dtype=np.complex64)
order1=METHODdata["PVM_EncSteps1"]+dim[1]/2                             
for i in range(0,dim[1]): FIDdata_tmp[:,order1[i],:]=FIDrawdata_CPX[:,i,:]
FIDrawdata_CPX = 0 #free memory  
order2=METHODdata["PVM_EncSteps2"]+dim[2]/2                             
for i in range(0,dim[2]): FIDdata[:,:,order2[i]]=FIDdata_tmp[:,:,i]
FIDdata_tmp = 0 #free memory  
print('.', end='') #progress indicator

#Hanning filter
percentage = 0.
npoints_x = int(float(dim[0])*percentage/100.)
hanning_x = np.empty(shape=(dim[0]),dtype=np.float32)
x_ = np.linspace (1./(npoints_x-1.)*np.pi/2.,(1.-1./(npoints_x-1))*np.pi/2.,num=npoints_x)
hanning_x [0:npoints_x] = np.power(np.sin(x_),2)
hanning_x [npoints_x:hanning_x.shape[0]-npoints_x] = 1
x_ = x_[::-1] # reverse x_
hanning_x [hanning_x.shape[0]-npoints_x:hanning_x.shape[0]] = np.power(np.sin(x_),2)
#print (hanning_x)
FIDdata[:,:,:] *= hanning_x [:,None,None]
npoints_y = int(float(dim[1])*percentage/100.)
hanning_y = np.empty(shape=(dim[1]),dtype=np.float32)
y_ = np.linspace (1./(npoints_y-1.)*np.pi/2.,(1.-1./(npoints_y-1))*np.pi/2.,num=npoints_y)
hanning_y [0:npoints_y] = np.power(np.sin(y_),2)
hanning_y [npoints_y:hanning_y.shape[0]-npoints_y] = 1
y_ = y_[::-1] # reverse y_
hanning_y [hanning_y.shape[0]-npoints_y:hanning_y.shape[0]] = np.power(np.sin(y_),2)
#print (hanning_x)
FIDdata[:,:,:] *= hanning_y [None,:,None]
npoints_z = int(float(dim[2])*percentage/100.)
hanning_z = np.empty(shape=(dim[2]),dtype=np.float32)
z_ = np.linspace (1./(npoints_z-1.)*np.pi/2.,(1.-1./(npoints_z-1))*np.pi/2.,num=npoints_z)
hanning_z [0:npoints_z] = np.power(np.sin(z_),2)
hanning_z [npoints_z:hanning_z.shape[0]-npoints_z] = 1
z_ = z_[::-1] # reverse z_
hanning_z [hanning_z.shape[0]-npoints_z:hanning_z.shape[0]] = np.power(np.sin(z_),2)
#print (hanning_x)
FIDdata[:,:,:] *= hanning_z [None,None,:]
print('.', end='') #progress indicator

# apply FOV offsets = (linear phase in k-space)
PackArrPhase1Offset=METHODdata["PVM_SPackArrPhase1Offset"]
SPackArrSliceOffset=METHODdata["PVM_SPackArrSliceOffset"]
realFOV = METHODdata["PVM_Fov"]*METHODdata["PVM_AntiAlias"]
phase_step1 = +2.*np.pi*float(PackArrPhase1Offset)/float(realFOV[1])
phase_step2 = -2.*np.pi*float(SPackArrSliceOffset)/float(realFOV[2])
mag = np.abs(FIDdata[:,:,:]); ph = np.angle(FIDdata[:,:,:])
for i in range(0,FIDdata.shape[1]): ph[:,i,:] -= float(i-int(FIDdata.shape[1]/2))*phase_step1
for j in range(0,FIDdata.shape[2]): ph[:,:,j] -= float(j-int(FIDdata.shape[2]/2))*phase_step2
FIDdata [:,:,:] = mag * np.exp(1j*ph)
print('.', end='') #progress indicator

#zero fill
zero_fill=2.
SpatResol=METHODdata["PVM_SpatResol"]/zero_fill
FIDdata_ZF = np.empty(shape=(int(dim[0]*zero_fill),int(dim[1]*zero_fill),
                             int(dim[2]*zero_fill)),dtype=np.complex64)
dim0start=int(dim[0]*(zero_fill-1)/2)
dim1start=int(dim[1]*(zero_fill-1)/2)
dim2start=int(dim[2]*(zero_fill-1)/2)
FIDdata_ZF[dim0start:dim0start+dim[0],dim1start:dim1start+dim[1],dim2start:dim2start+dim[2]] = \
    FIDdata[0:dim[0],0:dim[1],0:dim[2]]
FIDdata=FIDdata_ZF;
FIDdata_ZF = 0 #free memory 
dim=FIDdata.shape
print('.', end='') #progress indicator

#FFT (individually by axis, to save memory)
EchoPosition_raw=METHODdata["PVM_EchoPosition"]
EchoPosition_raw=50-(50-EchoPosition_raw)/zero_fill
EchoPosition=int(EchoPosition_raw/100.*dim[0])
IMGdata=FIDdata
FIDdata = 0 #free memory 
#do FFT with for loop  to save memory
for k in range(0,IMGdata.shape[1]):
    if k%30==0: print(',', end='') #progress indicator
    IMGdata[:,k,:] = np.roll(IMGdata[:,k,:], -EchoPosition, axis=(0))
    IMGdata[:,k,:] = np.fft.fft(IMGdata[:,k,:], axis=(0))
    IMGdata[:,k,:] = np.fft.fftshift(IMGdata[:,k,:], axes=(0))
for i in range(0,IMGdata.shape[0]):
    if i%30==0: print(',', end='') #progress indicator
    IMGdata[i,:,:] = np.fft.fftshift(IMGdata[i,:,:], axes=(0))
    IMGdata[i,:,:] = np.fft.fft(IMGdata[i,:,:], axis=(0))
    IMGdata[i,:,:] = np.fft.fftshift(IMGdata[i,:,:], axes=(0))
for i in range(0,IMGdata.shape[0]):
    if i%30==0: print(',', end='') #progress indicator
    IMGdata[i,:,:] = np.fft.fftshift(IMGdata[i,:,:], axes=(1))
    IMGdata[i,:,:] = np.fft.fft(IMGdata[i,:,:], axis=(1))
    IMGdata[i,:,:] = np.fft.fftshift(IMGdata[i,:,:], axes=(1))       
print('.', end='') #progress indicator

#throw out antialiasing
crop=METHODdata["PVM_AntiAlias"]
dim0start=int((dim[0]-dim[0]/crop[0])/2)
dim1start=int((dim[1]-dim[1]/crop[1])/2)
dim2start=int((dim[2]-dim[2]/crop[2])/2)
dim0end = int(dim0start+dim[0]/crop[0])
dim1end = int(dim1start+dim[1]/crop[1])
dim2end = int(dim2start+dim[2]/crop[2])
IMGdata = IMGdata[dim0start:dim0end,dim1start:dim1end,dim2start:dim2end]
print('.', end='') #progress indicator

#permute dimensions
#worx for PVM_SPackArrSliceOrient=sagittal, PVM_SPackArrReadOrient="H_F"
#this way results are comparabled to ImageJ's BrukerOpener plugin
if METHODdata["PVM_SPackArrSliceOrient"] == "sagittal":
    if METHODdata["PVM_SPackArrReadOrient"] == "H_F":
        SpatResol_perm = np.empty(shape=(3))
        SpatResol_perm[0]=SpatResol[1]
        SpatResol_perm[1]=SpatResol[0]
        SpatResol_perm[2]=SpatResol[2]
        IMGdata = np.transpose (IMGdata, axes=(1,0,2))
        IMGdata = np.rot90(IMGdata, k=2, axes=(0, 2)) # k=2 is a 180 degree rotation
    else:
        SpatResol_perm=SpatResol    
        print ('Warning: unknown Orientation',METHODdata["PVM_SPackArrSliceOrient"],
                METHODdata["PVM_SPackArrReadOrient"]);
        print ('         resulting images may be rotated incorrectly');
elif METHODdata["PVM_SPackArrSliceOrient"] == "axial":
    if METHODdata["PVM_SPackArrReadOrient"] == "L_R":
        SpatResol_perm = SpatResol
        IMGdata = IMGdata[::-1,:,:] # flip axis (axial L_R)
        IMGdata = IMGdata[:,:,::-1] # flip axis (axial L_R)
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
zstart=0; zend=int(IMGdata.shape[2]/N)
arr=np.abs(IMGdata[xstart:xend,ystart:yend,zstart:zend])
avg[0]=np.mean(arr)
std[0]=np.std(arr)
tresh[0]=avg[0] + 4*std[0]
xstart=int(IMGdata.shape[0]-IMGdata.shape[0]/N); xend=IMGdata.shape[0]
ystart=0; yend=int(IMGdata.shape[1]/N)
zstart=0; zend=int(IMGdata.shape[2]/N)
arr=np.abs(IMGdata[xstart:xend,ystart:yend,zstart:zend])
avg[1]=np.mean(arr)
std[1]=np.std(arr)
tresh[1]=avg[1] + 4*std[1]
xstart=0; xend=int(IMGdata.shape[0]/N)
ystart=int(IMGdata.shape[1]-IMGdata.shape[1]/N); yend=IMGdata.shape[1]
zstart=0; zend=int(IMGdata.shape[2]/N)
arr=np.abs(IMGdata[xstart:xend,ystart:yend,zstart:zend])
avg[2]=np.mean(arr)
std[2]=np.std(arr)
tresh[2]=avg[2] + 4*std[2]
xstart=int(IMGdata.shape[0]-IMGdata.shape[0]/N); xend=IMGdata.shape[0]
ystart=int(IMGdata.shape[1]-IMGdata.shape[1]/N); yend=IMGdata.shape[1]
zstart=0; zend=int(IMGdata.shape[2]/N)
arr=np.abs(IMGdata[xstart:xend,ystart:yend,zstart:zend])
avg[3]=np.mean(arr)
std[3]=np.std(arr)
tresh[3]=avg[3] + 4*std[3]
xstart=0; xend=int(IMGdata.shape[0]/N)
ystart=0; yend=int(IMGdata.shape[1]/N)
zstart=int(IMGdata.shape[2]-IMGdata.shape[2]/N); zend=IMGdata.shape[2]
arr=np.abs(IMGdata[xstart:xend,ystart:yend,zstart:zend])
avg[4]=np.mean(arr)
std[4]=np.std(arr)
tresh[4]=avg[4] + 4*std[4]
xstart=int(IMGdata.shape[0]-IMGdata.shape[0]/N); xend=IMGdata.shape[0]
ystart=0; yend=int(IMGdata.shape[1]/N)
zstart=int(IMGdata.shape[2]-IMGdata.shape[2]/N); zend=IMGdata.shape[2]
arr=np.abs(IMGdata[xstart:xend,ystart:yend,zstart:zend])
avg[5]=np.mean(arr)
std[5]=np.std(arr)
tresh[5]=avg[5] + 4*std[5]
xstart=0; xend=int(IMGdata.shape[0]/N)
ystart=int(IMGdata.shape[1]-IMGdata.shape[1]/N); yend=IMGdata.shape[1]
zstart=int(IMGdata.shape[2]-IMGdata.shape[2]/N); zend=IMGdata.shape[2]
arr=np.abs(IMGdata[xstart:xend,ystart:yend,zstart:zend])
avg[6]=np.mean(arr)
std[6]=np.std(arr)
tresh[6]=avg[6] + 4*std[6]
xstart=int(IMGdata.shape[0]-IMGdata.shape[0]/N); xend=IMGdata.shape[0]
ystart=int(IMGdata.shape[1]-IMGdata.shape[1]/N); yend=IMGdata.shape[1]
zstart=int(IMGdata.shape[2]-IMGdata.shape[2]/N); zend=IMGdata.shape[2]
arr=np.abs(IMGdata[xstart:xend,ystart:yend,zstart:zend])
avg[7]=np.mean(arr)
std[7]=np.std(arr)
tresh[7]=avg[7] + 4*std[7]
# old code trying to remove outliers and average
#avg1=np.mean(tresh); std1=np.std(tresh)
#remove1 = np.where(tresh>avg1+2*std1)
#avg=np.delete(avg,remove1)
#std=np.delete(std,remove1)
#tresh=np.delete(tresh,remove1)
#remove2 = np.where(tresh<avg1-2*std1)
#avg=np.delete(avg,remove2)
#std=np.delete(std,remove2)
#tresh=np.delete(tresh,remove2)
#threshold=np.mean(tresh)
# new code, just take the minimum
threshold=np.min(tresh)
mask =  abs(IMGdata [:,:,:]) > threshold

# SNR estimation
xstart=int((1/2.-1/8.)*IMGdata.shape[0]); xend=int((1/2.+1/8.)*IMGdata.shape[0])
xstart=int((1/2.-1/8.)*IMGdata.shape[1]); xend=int((1/2.+1/8.)*IMGdata.shape[1])
xstart=int((1/2.-1/8.)*IMGdata.shape[2]); xend=int((1/2.+1/8.)*IMGdata.shape[2])
arr=np.abs(IMGdata[xstart:xend,ystart:yend,zstart:zend])
SNR=np.amax(arr)/np.mean(avg)

#transform to int
IMGdata_ABS = np.abs(IMGdata);
IMGdata =0 # free memory 
maximum_ABS = np.amax(IMGdata_ABS)
IMGdata_ABS *= 32767./maximum_ABS 
IMGdata_ABS = IMGdata_ABS.astype(np.int16)
print('.', end='') #progress indicator

#save NIFTI
aff = np.eye(4)
aff[0,0] = SpatResol_perm[0]*1000; aff[0,3] = -(IMGdata_ABS.shape[0]/2)*aff[0,0]
aff[1,1] = SpatResol_perm[1]*1000; aff[1,3] = -(IMGdata_ABS.shape[1]/2)*aff[1,1]
aff[2,2] = SpatResol_perm[2]*1000; aff[2,3] = -(IMGdata_ABS.shape[2]/2)*aff[2,2]
NIFTIimg_ABS = nib.Nifti1Image(IMGdata_ABS, aff)
NIFTIimg_ABS.header['sform_code']=1
NIFTIimg_ABS.header['qform_code']=1
NIFTIimg_ABS.header.set_slope_inter(1,0)
#write unmasked
try:
    print('.', end='') #progress indicator
    nib.save(NIFTIimg_ABS, os.path.join(os.path.dirname(FIDfile),OrigFilename+'_MAGNT.nii.gz'))
except:
    print ('\nERROR:  problem while writing results'); sys.exit(1)
NIFTIimg_ABS = nib.Nifti1Image(IMGdata_ABS*mask, aff)
NIFTIimg_ABS.header['sform_code']=1
NIFTIimg_ABS.header['qform_code']=1
NIFTIimg_ABS.header.set_slope_inter(1,0)
#write masked
try:
    print('.', end='') #progress indicator
    nib.save(NIFTIimg_ABS, os.path.join(os.path.dirname(FIDfile),OrigFilename+'_MAGNT_masked.nii.gz'))
except:
    print ('\nERROR:  problem while writing results'); sys.exit(1)
print ('\nSuccessfully written output files '+OrigFilename+'_MAGNT*.nii.gz')

#calculate volumes
vol=np.empty(shape=IMGdata_ABS.shape[0],dtype=np.float)
vol_percentage=np.empty(shape=IMGdata_ABS.shape[0],dtype=np.float)
for i in range(0,IMGdata_ABS.shape[0]):
   arr=IMGdata_ABS[i,:,:]*mask[i,:,:]
   arr=arr[arr.nonzero()]
   vol[i] = np.count_nonzero(IMGdata_ABS[i,:,:]*mask[i,:,:])
   # local normalization over N slices
   N=int(2*zero_fill)
   start=max(i-N,0)
   end=min(i+N,IMGdata_ABS.shape[0])      
   maximum=np.amax(IMGdata_ABS[start:end,:,:]*mask[start:end,:,:])
   if arr.size>0: partial_volume=np.average(arr)/maximum
   else:partial_volume=0
   # calculate the volume
   vol[i] *= SpatResol_perm[0]*SpatResol_perm[1]*SpatResol_perm[2] # result is in mm^3/s
   vol[i] /= 1000. # convert mm^3/s ---> ml/s
   vol[i] *= partial_volume
   vol_percentage[i] = vol[i]/(11.*(SpatResol_perm[0]/10.))*100. # cores have 11cm^2 axial area, resolution is in mm

#write volume results
if METHODdata["PVM_SPackArrSliceOrient"] == "axial" and METHODdata["PVM_SPackArrReadOrient"] == "L_R":
    with open(os.path.join(os.path.dirname(FIDfile),OrigFilename+'_Volumes.txt'), "w") as text_file:
        start = int((1/2.-1/8.)*IMGdata_ABS.shape[0])
        end   = int((1/2.+1/8.)*IMGdata_ABS.shape[0])
        vol_avg = np.mean(vol_percentage[start:end])
        vol_std = np.std(vol_percentage[start:end])
        print("Volume:    %0.2f%% %s %0.2f" % (vol_avg, chr(241), vol_std))
        print("SNR:       %0.2f" % (SNR))
        print("Threshold: %0.2f%%" % (threshold/maximum_ABS*100))
        text_file.write("Volume average:\t%0.2f\t%%\t%0.2f\t%%\n" % (vol_avg,vol_std))
        text_file.write("estimated SNR: \t%0.2f\n" % (SNR))
        text_file.write("Threshold:     \t%0.2f%%\n\n\n" % (threshold/maximum_ABS*100))
        text_file.write("Volumes per slice (X):\n")
        for i in range(0,IMGdata_ABS.shape[0]):
            text_file.write("Slice %d:\t%0.6f\tml\t%0.2f\t%%\n" % (i, vol[i],vol_percentage[i]))
        text_file.write("\n")        
else:
    try: os.remove (os.path.join(os.path.dirname(FIDfile),OrigFilename+'_Volumes.txt'))
    except: pass #silent
    print ('Warning: textfile with volume measures not written (unknown Orientation)');   
#end


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