#
# calculates IVIM parameters
#  slow diffusion time constant (aka ADC)
#  fast diffusion time constant 
#  perfusion fraction
#      
#
# ----- VERSION HISTORY -----
#
# Version 0.1 - 20, September 2019
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
# to do
# - calculate individual threshold for each image type B(B0, Bs)
# - calc kernel for result filtering based on resolution, to work for 2D and 3D
# - test interpolation 2x
# - only enter curvefit if not all points are zero
# - parallelize
# - use current parameter maps as initial condition for true bi-exponential fit
# - bayasian fit (?)



from __future__ import print_function
try: import win32gui, win32console
except: pass #silent
from math import ceil, floor
import sys
import os
import warnings
import numpy as np
import nibabel as nib
if getattr( sys, 'frozen', False ): # running as pyinstaller bundle
   from scipy_extract import zoom
   from scipy_extract import median_filter
   from scipy_extract import gaussian_filter 
   from scipy_extract import curve_fit 
else: # running native python
   from scipy.ndimage import zoom 
   from scipy.ndimage import median_filter 
   from scipy.ndimage import gaussian_filter 
   from scipy.optimize import curve_fit   

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

#define fit functions
def expdecay(x,a,D):
    return a* np.exp(-x*D)
def FIT (x,y): 
    bad_fit=False
    initial_conditions=[max(y), 1/((max(x)+min(x))/2)]
    try: 
        pars,covar = curve_fit(expdecay,x,y,p0=initial_conditions,maxfev=100*(len(x)+1))
        if (np.size(pars)==1):  bad_fit=True
        if (np.size(covar)==1): bad_fit=True
    except RuntimeError: bad_fit=True
    if bad_fit: A=0;Aerr=0;D=0;Derr=0                           #fit failed
    elif (pars[0]<0) or (pars[1]<0): A=0;Aerr=0;D=0;Derr=0      #does not make sense
    elif  covar[1,1]>0.5*pars[1]: A=0;Aerr=0;D=0;Derr=0         #error too large    
    #elif  covar[0,0]>50*pars[0]: A=0;Aerr=0;D=0;Derr=0          #error too large        
    else: A=pars[0]; Aerr=covar[0,0];D=pars[1]; Derr=covar[1,1] #fit OK
    return A, Aerr, D, Derr



    
#general initialization stuff  
Program_name = os.path.basename(sys.argv[0]); 
if Program_name.find('.')>0: Program_name = Program_name[:Program_name.find('.')]
python_version=str(sys.version_info[0])+'.'+str(sys.version_info[1])+'.'+str(sys.version_info[2])
# sys.platform = [linux2, win32, cygwin, darwin, os2, os2emx, riscos, atheos, freebsd7, freebsd8]
if sys.platform=="win32": os.system("title "+Program_name)
warnings.filterwarnings("ignore")
    
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
    
#intercatively choose input NIFTI file
InputFile = askopenfilename(title="Choose NIFTI file", filetypes=[("NIFTI files",('MAGNT.nii.gz'))])
if InputFile=="":print ('ERROR: No input file specified'); sys.exit(2)
InputFile = os.path.abspath(InputFile)
TKwindows.update()
try: win32gui.SetForegroundWindow(win32console.GetConsoleWindow())
except: pass #silent
dirname  = os.path.dirname(InputFile)
basename = os.path.basename(InputFile)
basename = os.path.splitext(basename)[0]
basename = basename[0:basename.rfind('_MAG')]

print ('Reading NIFTI file'); sys.stdout.flush()
img = nib.load(InputFile)
IMGdata = img.get_data().astype(np.float32)
SpatResol = np.asarray(img.header.get_zooms())
Shape = np.asarray(img.header.get_data_shape())
affine = img.affine
xyzt_units1  = img.header.get_xyzt_units()[0]
xyzt_units2  = img.header.get_xyzt_units()[1]
sform = int(img.header['sform_code'])
qform = int(img.header['qform_code'])
bvec_str = str(img.header['descrip']).replace("b'","").replace("'","")
if bvec_str[0:4] != "B = ":
   print ("ERROR: unable to identify B values from NIFTI header"); sys.exit(1)
bvec_str = bvec_str[4:].split(',')
bvec = np.zeros (len(bvec_str), dtype=np.int)
for i in range (len(bvec_str)):
   try: bvec[i] = float(bvec_str[i])
   except: print ("ERROR: unable to parse B-values in NIFTY header"); sys.exit(1)
bvec_str = np.array2string(bvec,max_line_width=1000)
bvec_str = bvec_str.replace('.]','').replace(']','').replace('[','')
bvec_str = bvec_str.replace('. ',' ').replace('   ',' ').replace('  ',' ')
print ("B's ="+bvec_str); sys.stdout.flush()


#check some stuff
if len(IMGdata.shape) != 4:
    print ('ERROR: 4D NIFTI expected, 3D NIFTI found, this is probably not an IVIM file'); 
    sys.exit(1)
if bvec.shape[0]!=IMGdata.shape[3]:
    print ("ERROR: number of B's in header unequal data dimension "); 
    sys.exit(1)
if np.unique(bvec).shape[0]<3:
    print ("ERROR: need at least 3 unique B's"); 
    sys.exit(1)
if np.amin(IMGdata)<0 or abs(np.amax(IMGdata)-np.pi)<0.2:
    print ("ERROR: this looks like a Phase Image"); 
    sys.exit(1)
    
# check if already masked
N=10 # use 10% at the corners of the FOV
tresh=np.empty(shape=8,dtype=np.float)
avg=np.empty(shape=8,dtype=np.float)
std=np.empty(shape=8,dtype=np.float)
xstart=0; xend=int(IMGdata.shape[0]/N)
ystart=0; yend=int(IMGdata.shape[1]/N)
zstart=0; zend=int(ceil(float(IMGdata.shape[2])/float(N)))
arr=np.abs(IMGdata[xstart:xend,ystart:yend,zstart:zend,:])
if np.count_nonzero(arr==0)>0.2*arr.size: #masked input image not OK
    print ("ERROR: this looks like an image with noise mask"); 
    sys.exit(1)    

    
# --------------------------- start doing something --------------------------

'''
#Increase resolution
print ('\nIncrease resolution 2x'); sys.stdout.flush()
dim = IMGdata.shape
zoomed = np.zeros((IMGdata.shape[0]*2,IMGdata.shape[1]*2,IMGdata.shape[2],IMGdata.shape[3]),dtype=np.float32) 
for i in range(IMGdata.shape[3]):
   print ('.',end=''); sys.stdout.flush() #progress indicator
   zoomed[:,:,:,i] = zoom(IMGdata[:,:,:,i],[2,2,1],order=2)
zoomed[zoomed<0]=0 #remove negative values
IMGdata = zoomed
'''

# calculate mask
# use noise in all 8 corners to establish threshold
N=10 # use 10% at the corners of the FOV
std_factor = 4 # thresh = avg + std_factor*std
thresh=np.empty(shape=8,dtype=np.float)
avg=np.empty(shape=8,dtype=np.float)
std=np.empty(shape=8,dtype=np.float)
xstart=0; xend=int(IMGdata.shape[0]/N)
ystart=0; yend=int(IMGdata.shape[1]/N)
zstart=0; zend=int(ceil(float(IMGdata.shape[2])/float(N)))
arr=np.abs(IMGdata[xstart:xend,ystart:yend,zstart:zend,:])
avg[0]=np.mean(arr)
std[0]=np.std(arr)
thresh[0]=avg[0] + std_factor*std[0]
xstart=int(IMGdata.shape[0]-IMGdata.shape[0]/N); xend=IMGdata.shape[0]
ystart=0; yend=int(IMGdata.shape[1]/N)
zstart=0; zend=int(ceil(float(IMGdata.shape[2])/float(N)))
arr=np.abs(IMGdata[xstart:xend,ystart:yend,zstart:zend,:])
avg[1]=np.mean(arr)
std[1]=np.std(arr)
thresh[1]=avg[1] + std_factor*std[1]
xstart=0; xend=int(IMGdata.shape[0]/N)
ystart=int(IMGdata.shape[1]-IMGdata.shape[1]/N); yend=IMGdata.shape[1]
zstart=0; zend=int(ceil(float(IMGdata.shape[2])/float(N)))
arr=np.abs(IMGdata[xstart:xend,ystart:yend,zstart:zend,:])
avg[2]=np.mean(arr)
std[2]=np.std(arr)
thresh[2]=avg[2] + std_factor*std[2]
xstart=int(IMGdata.shape[0]-IMGdata.shape[0]/N); xend=IMGdata.shape[0]
ystart=int(IMGdata.shape[1]-IMGdata.shape[1]/N); yend=IMGdata.shape[1]
zstart=0; zend=int(ceil(float(IMGdata.shape[2])/float(N)))
arr=np.abs(IMGdata[xstart:xend,ystart:yend,zstart:zend,:])
avg[3]=np.mean(arr)
std[3]=np.std(arr)
thresh[3]=avg[3] + std_factor*std[3]
xstart=0; xend=int(IMGdata.shape[0]/N)
ystart=0; yend=int(IMGdata.shape[1]/N)
zstart=int(floor(float(IMGdata.shape[2])-float(IMGdata.shape[2])/float(N))); zend=IMGdata.shape[2]
arr=np.abs(IMGdata[xstart:xend,ystart:yend,zstart:zend,:])
avg[4]=np.mean(arr)
std[4]=np.std(arr)
thresh[4]=avg[4] + std_factor*std[4]
xstart=int(IMGdata.shape[0]-IMGdata.shape[0]/N); xend=IMGdata.shape[0]
ystart=0; yend=int(IMGdata.shape[1]/N)
zstart=int(floor(float(IMGdata.shape[2])-float(IMGdata.shape[2])/float(N))); zend=IMGdata.shape[2]
arr=np.abs(IMGdata[xstart:xend,ystart:yend,zstart:zend])
avg[5]=np.mean(arr)
std[5]=np.std(arr)
thresh[5]=avg[5] + std_factor*std[5]
xstart=0; xend=int(IMGdata.shape[0]/N)
ystart=int(IMGdata.shape[1]-IMGdata.shape[1]/N); yend=IMGdata.shape[1]
zstart=int(floor(float(IMGdata.shape[2])-float(IMGdata.shape[2])/float(N))); zend=IMGdata.shape[2]
arr=np.abs(IMGdata[xstart:xend,ystart:yend,zstart:zend,:])
avg[6]=np.mean(arr)
std[6]=np.std(arr)
thresh[6]=avg[6] + std_factor*std[6]
xstart=int(IMGdata.shape[0]-IMGdata.shape[0]/N); xend=IMGdata.shape[0]
ystart=int(IMGdata.shape[1]-IMGdata.shape[1]/N); yend=IMGdata.shape[1]
zstart=int(floor(float(IMGdata.shape[2])-float(IMGdata.shape[2])/float(N))); zend=IMGdata.shape[2]
arr=np.abs(IMGdata[xstart:xend,ystart:yend,zstart:zend,:])
avg[7]=np.mean(arr)
std[7]=np.std(arr)
thresh[7]=avg[7] + std_factor*std[7]
mask_threshold=np.min(thresh)
mask =  IMGdata [:,:,:,:] > mask_threshold
print ('.',end=''); sys.stdout.flush()
#filter mask
mask = median_filter (mask, size = (3,3,1,1)) #filter in 2D spatial dimension
mask = mask > 0.8  # takes out the 0.5 case of the median filter
print ('.',end=''); sys.stdout.flush()
#appy mask
IMGdata = IMGdata*mask
print ('.',end=''); sys.stdout.flush()

# reshape flatten
dim = IMGdata.shape
IMGdata = np.reshape(IMGdata, (dim[0]*dim[1]*dim[2],dim[3]))
print ('.',end=''); sys.stdout.flush()

#fit slow diffusion component (aka ADC)
print ('\nFitting slow diffusion component'); sys.stdout.flush()
min_b_index = int(ceil(2./3.*bvec.shape[0]))                # use only high b-values
if min_b_index>bvec.shape[0]-4: min_b_index=bvec.shape[0]-4 # at least 4
Dslow_clip = 1./(bvec[min_b_index]*5.)
if Dslow_clip<3e-3: Dslow_clip=3e-3   # water is ~2e-3 (50% safety margin) 
print ("Using B's",str(bvec[min_b_index:]).replace('[','').replace(']',''))
print ("Dslow clip is",Dslow_clip*1e3)
data_Dslow = np.zeros(dtype=np.float32, shape=(dim[0]*dim[1]*dim[2]))
data_Aslow = np.zeros(dtype=np.float32, shape=(dim[0]*dim[1]*dim[2]))
progress_tag = int(dim[0]*dim[1]*dim[2]/70)
for i in range(dim[0]*dim[1]*dim[2]):
   if i%progress_tag==0: print ('.',end=''); sys.stdout.flush()
   nz = np.nonzero (IMGdata [i,min_b_index:])     
   if nz[0].shape[0]>=3: #need at least 3 unmasked points
      B_temp = bvec[min_b_index:][nz]
      IMGdata_temp = IMGdata [i,min_b_index:][nz]    
      data_Aslow[i], Aerr, data_Dslow [i], ADCerr = FIT (B_temp, IMGdata_temp)
print (''); sys.stdout.flush()

#subtract long component
for i in range(dim[0]*dim[1]*dim[2]):
    IMGdata [i,:] -= expdecay(bvec,data_Aslow[i],data_Dslow [i])
#reapply threshold 
#CAUTION: this also remove negative points
# if the fitting of th slow compontent overestimates the amplitude
# and/or underestimates the diffusion coefficient this will remove points
# that could potientially be still meaningfull when fitted
# with "exp(-x*D)-offset" (we currently don't use an offset)
IMGdata[IMGdata<mask_threshold]=0 

#fit fast diffusion component
print ('Fitting fast diffusion component'); sys.stdout.flush()
max_b_index = int(floor(2./3.*bvec.shape[0]))           # use only low b-values
if max_b_index<3: max_b_index=3                         # at least 4 (first index is 0)
if max_b_index>=bvec.shape[0]: max_b_index=bvec.shape[0]-1 # all
Dfast_thresh = 1./(bvec[max_b_index-1]*2.)
if Dfast_thresh<3e-3: Dfast_thresh=3e-3   # water is ~2e-3 (50% safety margin)
print ("Using B's",str(bvec[0:max_b_index]).replace('[','').replace(']',''))
print ("Dfast threshold is",Dfast_thresh*1e3)
data_Dfast = np.zeros(dtype=np.float32, shape=(dim[0]*dim[1]*dim[2]))
data_Afast = np.zeros(dtype=np.float32, shape=(dim[0]*dim[1]*dim[2]))
progress_tag = int(dim[0]*dim[1]*dim[2]/70)
for i in range(dim[0]*dim[1]*dim[2]):
   if i%progress_tag==0: print ('.',end=''); sys.stdout.flush()
   nz = np.nonzero (IMGdata [i,0:max_b_index])     
   if nz[0].shape[0]>=3: #need at least 3 unmasked points
      B_temp = bvec[0:max_b_index][nz]
      IMGdata_temp = IMGdata [i,0:max_b_index][nz]    
      data_Afast[i], Arr, data_Dfast [i], Derr = FIT (B_temp, IMGdata_temp)
print (''); sys.stdout.flush()

#calculate Perfusion Fraction
data_Pfrac = np.zeros(dtype=np.float32, shape=(dim[0]*dim[1]*dim[2]))
nz = np.nonzero (data_Afast+data_Aslow)
data_Pfrac[nz] = data_Afast[nz]/(data_Afast[nz]+data_Aslow[nz])*100.

#reshape to original
data_Dslow = np.reshape(data_Dslow,(dim[0],dim[1],dim[2]))
data_Dfast = np.reshape(data_Dfast,(dim[0],dim[1],dim[2]))
data_Pfrac = np.reshape(data_Pfrac,(dim[0],dim[1],dim[2]))
IMGdata    = np.reshape(IMGdata,   (dim[0],dim[1],dim[2],dim[3]))

#clip slow (aka ADC)
data_Dslow[data_Dslow>Dslow_clip]=Dslow_clip

#clip fast 
data_Dfast[data_Dfast<Dfast_thresh]=0
data_Pfrac[data_Dfast<Dfast_thresh]=0

#filter slow
data_Dslow = median_filter (data_Dslow, size = (3,3,1))
for i in range(data_Dslow.shape[2]): 
  data_Dslow[:,:,i] = gaussian_filter(data_Dslow[:,:,i], sigma=0.7, truncate=3)
#filter fast
data_Dfast = median_filter (data_Dfast, size = (3,3,1))
for i in range(data_Dfast.shape[2]): 
  data_Dfast[:,:,i] = gaussian_filter(data_Dfast[:,:,i], sigma=0.7, truncate=3)
#filter Perfusion Fraction
data_Pfrac = median_filter (data_Pfrac, size = (3,3,1))
for i in range(data_Pfrac.shape[2]): 
  data_Pfrac[:,:,i] = gaussian_filter(data_Pfrac[:,:,i], sigma=0.7, truncate=3)  

'''
#filter slow
data_Dslow = median_filter (data_Dslow, size = (5,5,1))
for i in range(data_Dslow.shape[2]): 
  data_Dslow[:,:,i] = gaussian_filter(data_Dslow[:,:,i], sigma=1.2, truncate=3)
#filter fast
data_Dfast = median_filter (data_Dfast, size = (5,5,1))
for i in range(data_Dfast.shape[2]): 
  data_Dfast[:,:,i] = gaussian_filter(data_Dfast[:,:,i], sigma=1.2, truncate=3)
#filter Perfusion Fraction
data_Pfrac = median_filter (data_Pfrac, size = (5,5,1))
for i in range(data_Pfrac.shape[2]): 
  data_Pfrac[:,:,i] = gaussian_filter(data_Pfrac[:,:,i], sigma=1.2, truncate=3)  
#decrease resolution
print ('Decrease resolution 2x'); sys.stdout.flush()
data_Dslow  = zoom(data_Dslow,[0.5,0.5,1],order=1)
data_Dfast  = zoom(data_Dfast,[0.5,0.5,1],order=1)
data_Pfrac  = zoom(data_Pfrac ,[0.5,0.5,1],order=1)   
'''

#convert ADC unit to 10^-3 mm^2/s
#(we suppose that Bs are in mm^2/s)
data_Dslow *= 1e3 
#usual in-vivo values:
#    white matter: 0.6-0.8
#    grey matter: 0.8-1.0
#    CSF: 3.0-3.5
#    restriction (e.g. tumor, stroke) <0.7 
data_Dfast *= 1e3 # use same unit
         
#transform to int
max_Dslow = np.amax(data_Dslow);
data_Dslow *= 32767./max_Dslow
data_Dslow = data_Dslow.astype(np.int16)
max_Dfast = np.amax(data_Dfast);
data_Dfast *= 32767./max_Dfast
data_Dfast = data_Dfast.astype(np.int16)
max_Pfrac = np.amax(data_Pfrac);
data_Pfrac *= 32767./max_Pfrac
data_Pfrac = data_Pfrac.astype(np.int16)
max_Residual = np.amax(IMGdata);
IMGdata *= 32767./max_Residual
IMGdata = IMGdata.astype(np.int16)

#Write slow diffusion component (aka ADC)
ADC_img = nib.Nifti1Image(data_Dslow, affine)
ADC_img.header.set_slope_inter(max_Dslow/32767.,0)
ADC_img.header.set_xyzt_units(xyzt_units1,xyzt_units2)
ADC_img.set_sform(affine, sform)
ADC_img.set_qform(affine, qform)
try: 
   nib.save(ADC_img, os.path.join(dirname,basename+'_Dslow.nii.gz'))
except: print ('ERROR:  problem while writing results'); exit(1)

#Write fast diffusion component
Dslow_img = nib.Nifti1Image(data_Dfast, affine)
Dslow_img.header.set_slope_inter(max_Dfast/32767.,0)
Dslow_img.header.set_xyzt_units(xyzt_units1,xyzt_units2)
Dslow_img.set_sform(affine, sform)
Dslow_img.set_qform(affine, qform)
try: 
   nib.save(Dslow_img, os.path.join(dirname,basename+'_Dfast.nii.gz'))
except: print ('ERROR:  problem while writing results'); exit(1)

#Write Perfusion Fraction
Pfrac_img = nib.Nifti1Image(data_Pfrac, affine)
Pfrac_img.header.set_slope_inter(max_Pfrac/32767.,0)
Pfrac_img.header.set_xyzt_units(xyzt_units1,xyzt_units2)
Pfrac_img.set_sform(affine, sform)
Pfrac_img.set_qform(affine, qform)
try: 
   nib.save(Pfrac_img, os.path.join(dirname,basename+'_Pfrac.nii.gz'))
except: print ('ERROR:  problem while writing results'); exit(1)

'''
#Write Residual after subtraction of long diffusion component
Residual_img = nib.Nifti1Image(IMGdata, affine)
Residual_img.header.set_slope_inter(max_Residual/32767.,0)
Residual_img.header.set_xyzt_units(xyzt_units1,xyzt_units2)
Residual_img.set_sform(affine, sform)
Residual_img.set_qform(affine, qform)
try: 
   nib.save(Residual_img, os.path.join(dirname,basename+'_Residual.nii.gz'))
except: print ('ERROR:  problem while writing results'); exit(1)
'''

print ('Successfully written output files'); sys.stdout.flush()

os.system("pause") # windows
