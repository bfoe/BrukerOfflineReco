#
# calculates IVIM parameters
#  slow diffusion coefficient (aka ADC)
#  fast diffusion coefficient 
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


from __future__ import print_function
try: import win32gui, win32console
except: pass #silent
from math import ceil, floor
import sys
import os
import fnmatch
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

#define mono-exponential fit function
def expdecay(x,a,D):
    if a<=0 or D<=0: result = 0
    else: result = a*np.exp(-x*D)
    return result
def FIT (x,y,initial_conditions): 
    bad_fit=False
    try: pars,covar = curve_fit(expdecay,x,y,p0=initial_conditions)
    except RuntimeError: bad_fit=True    
    if bad_fit: A=0;D=0       #fit failed
    else: A=pars[0];D=pars[1] #fit OK
    return A, D
    
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
bval_str = str(img.header['descrip']).replace("b'","").replace("'","")
if bval_str[0:4] != "B = ":
   print ("Unable to identify B values from NIFTI header")
   # alternatively try to find bval file
   bval_filename = ""
   for root, dirs, files in os.walk(dirname):
      for name in files:
         if fnmatch.fnmatch(name, basename+"*.bval") and root==dirname: 
            bval_filename=name
   if bval_filename!="": #found
      print ("Reading B's from file",bval_filename)
      with open(os.path.join(dirname,bval_filename), "r") as f:
         line = f.readline() #read just one line
         bval_str = "B = "+line.replace(" ",",")         
   else: # not found, nothing else to do, sorry
      print ("ERROR: unable to identify B values, no suitable .bval file found"); sys.exit(1)
bval_str = bval_str[4:].split(',')
bval = np.zeros (len(bval_str), dtype=np.float32)
for i in range (len(bval_str)):
   try: bval[i] = float(bval_str[i])
   except: print ("ERROR: unable to parse B-values in NIFTY header"); sys.exit(1)
bval_str = np.array2string(bval.astype(np.int16),max_line_width=1000)
bval_str = bval_str.replace('.]','').replace(']','').replace('[','')
bval_str = bval_str.replace('. ',' ').replace('   ',' ').replace('  ',' ')

#check some stuff
if len(IMGdata.shape) != 4:
    print ('ERROR: 4D NIFTI expected, 3D NIFTI found, this is probably not an IVIM file'); 
    sys.exit(1)
if bval.shape[0]!=IMGdata.shape[3]:
    print ("ERROR: number of B's in header unequal data dimension "); 
    sys.exit(1)
if np.unique(bval).shape[0]<7:
    print ("ERROR: need at least 7 unique B's"); 
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

#sort in ascending b-value order (just in case)
order = np.argsort (bval)
bval = bval[order]
IMGdata = IMGdata [:,:,:,order]
print ("B's ="+bval_str); sys.stdout.flush()



# --------------------------- start doing something --------------------------

# calculate mask
# use noise in all 8 corners to establish threshold
N=10 # use 10% at the corners of the FOV
std_factor = 2 # thresh = avg + std_factor*std
mask = np.zeros(shape=IMGdata.shape,dtype=np.float)
mask_all = np.zeros(shape=IMGdata.shape[0:3],dtype=np.int16)
mask_threshold=np.zeros(shape=IMGdata.shape[3],dtype=np.float)
for i in range (IMGdata.shape[3]):
    thresh=np.empty(shape=8,dtype=np.float)
    avg=np.empty(shape=8,dtype=np.float)
    std=np.empty(shape=8,dtype=np.float)
    xstart=0; xend=int(IMGdata.shape[0]/N)
    ystart=0; yend=int(IMGdata.shape[1]/N)
    zstart=0; zend=int(ceil(float(IMGdata.shape[2])/float(N)))
    arr=np.abs(IMGdata[xstart:xend,ystart:yend,zstart:zend,i])
    avg[0]=np.mean(arr)
    std[0]=np.std(arr)
    thresh[0]=avg[0] + std_factor*std[0]
    xstart=int(IMGdata.shape[0]-IMGdata.shape[0]/N); xend=IMGdata.shape[0]
    ystart=0; yend=int(IMGdata.shape[1]/N)
    zstart=0; zend=int(ceil(float(IMGdata.shape[2])/float(N)))
    arr=np.abs(IMGdata[xstart:xend,ystart:yend,zstart:zend,i])
    avg[1]=np.mean(arr)
    std[1]=np.std(arr)
    thresh[1]=avg[1] + std_factor*std[1]
    xstart=0; xend=int(IMGdata.shape[0]/N)
    ystart=int(IMGdata.shape[1]-IMGdata.shape[1]/N); yend=IMGdata.shape[1]
    zstart=0; zend=int(ceil(float(IMGdata.shape[2])/float(N)))
    arr=np.abs(IMGdata[xstart:xend,ystart:yend,zstart:zend,i])
    avg[2]=np.mean(arr)
    std[2]=np.std(arr)
    thresh[2]=avg[2] + std_factor*std[2]
    xstart=int(IMGdata.shape[0]-IMGdata.shape[0]/N); xend=IMGdata.shape[0]
    ystart=int(IMGdata.shape[1]-IMGdata.shape[1]/N); yend=IMGdata.shape[1]
    zstart=0; zend=int(ceil(float(IMGdata.shape[2])/float(N)))
    arr=np.abs(IMGdata[xstart:xend,ystart:yend,zstart:zend,i])
    avg[3]=np.mean(arr)
    std[3]=np.std(arr)
    thresh[3]=avg[3] + std_factor*std[3]
    xstart=0; xend=int(IMGdata.shape[0]/N)
    ystart=0; yend=int(IMGdata.shape[1]/N)
    zstart=int(floor(float(IMGdata.shape[2])-float(IMGdata.shape[2])/float(N))); zend=IMGdata.shape[2]
    arr=np.abs(IMGdata[xstart:xend,ystart:yend,zstart:zend,i])
    avg[4]=np.mean(arr)
    std[4]=np.std(arr)
    thresh[4]=avg[4] + std_factor*std[4]
    xstart=int(IMGdata.shape[0]-IMGdata.shape[0]/N); xend=IMGdata.shape[0]
    ystart=0; yend=int(IMGdata.shape[1]/N)
    zstart=int(floor(float(IMGdata.shape[2])-float(IMGdata.shape[2])/float(N))); zend=IMGdata.shape[2]
    arr=np.abs(IMGdata[xstart:xend,ystart:yend,zstart:zend,i])
    avg[5]=np.mean(arr)
    std[5]=np.std(arr)
    thresh[5]=avg[5] + std_factor*std[5]
    xstart=0; xend=int(IMGdata.shape[0]/N)
    ystart=int(IMGdata.shape[1]-IMGdata.shape[1]/N); yend=IMGdata.shape[1]
    zstart=int(floor(float(IMGdata.shape[2])-float(IMGdata.shape[2])/float(N))); zend=IMGdata.shape[2]
    arr=np.abs(IMGdata[xstart:xend,ystart:yend,zstart:zend,i])
    avg[6]=np.mean(arr)
    std[6]=np.std(arr)
    thresh[6]=avg[6] + std_factor*std[6]
    xstart=int(IMGdata.shape[0]-IMGdata.shape[0]/N); xend=IMGdata.shape[0]
    ystart=int(IMGdata.shape[1]-IMGdata.shape[1]/N); yend=IMGdata.shape[1]
    zstart=int(floor(float(IMGdata.shape[2])-float(IMGdata.shape[2])/float(N))); zend=IMGdata.shape[2]
    arr=np.abs(IMGdata[xstart:xend,ystart:yend,zstart:zend,i])
    avg[7]=np.mean(arr)
    std[7]=np.std(arr)
    thresh[7]=avg[7] + std_factor*std[7]
    mask_threshold[i]=np.min(thresh)
    mask[:,:,:,i] = IMGdata [:,:,:,i] > mask_threshold[i]
    IMGdata [:,:,:,i] *= mask[:,:,:,i] # apply mask    
print ('.',end=''); sys.stdout.flush()
mask_all = np.sum(mask, axis=3)
mask_all = mask_all>=2 # need at leas 2 good points
mask_all = mask_all.astype(np.int16)

# reshape flatten
dim = IMGdata.shape
IMGdata  = np.reshape(IMGdata, (dim[0]*dim[1]*dim[2],dim[3]))
mask_all = np.reshape(mask_all,(dim[0]*dim[1]*dim[2]))
print ('.',end=''); sys.stdout.flush()

# vector reduction
temp_IMGdata = IMGdata [mask_all>0,:]

#fit slow diffusion component (aka ADC)
print ('\nFitting slow diffusion component'); sys.stdout.flush()
min_b_index = int(ceil(2./3.*bval.shape[0]))                # use only high b-values
if min_b_index>bval.shape[0]-4: min_b_index=bval.shape[0]-4 # at least 4
Dslow_clip = 1./(bval[min_b_index]*5.)
if Dslow_clip<3e-3: Dslow_clip=3e-3   # water is ~2e-3 (50% safety margin) 
print ("Using B's",str(bval[min_b_index:].astype(np.int16))\
       .replace('[','').replace(']','').replace('   ',' ').replace('  ',' '))
print ("Dslow clip is",Dslow_clip*1e3)
temp_Aslow = np.zeros(shape=temp_IMGdata.shape[0], dtype=np.float32)
temp_Dslow = np.zeros(shape=temp_IMGdata.shape[0], dtype=np.float32)
progress_tag = int(temp_IMGdata.shape[0]/70)
for i in range(temp_IMGdata.shape[0]):
   if i%progress_tag==0: print ('.',end=''); sys.stdout.flush()
   nz = np.nonzero (temp_IMGdata [i,min_b_index:])     
   if nz[0].shape[0]>=3: #permit down to minimum of 3 unmasked points 
      B_nz = bval[min_b_index:][nz]
      IMGdata_nz = temp_IMGdata [i,min_b_index:][nz] 
      initial_conditions=[max(IMGdata_nz), 1./bval[-1]]
      temp_Aslow[i], temp_Dslow [i] = FIT (B_nz, IMGdata_nz, initial_conditions)
print (''); sys.stdout.flush() 

#subtract long component
IMGdata_residual = np.zeros(shape=temp_IMGdata.shape, dtype=np.float32)
for j in range(temp_IMGdata.shape[0]):
    IMGdata_residual[j,:] = temp_IMGdata[j,:] - expdecay(bval,temp_Aslow[j],temp_Dslow [j])
for i in range(temp_IMGdata.shape[1]):   
    #reapply threshold 
    #CAUTION: this also remove negative points
    # if the fitting of th slow compontent overestimates the amplitude
    # and/or underestimates the diffusion coefficient this will remove points
    # that could potientially be still meaningfull when fitted
    # with "exp(-x*D)-offset" (we currently don't use an offset)
    IMGdata_residual[:,i][IMGdata_residual[:,i]<mask_threshold[i]] = 0 

#fit fast diffusion component
print ('Fitting fast diffusion component'); sys.stdout.flush()
max_b_index = int(floor(2./3.*bval.shape[0]))           # use only low b-values
if max_b_index<3: max_b_index=3                         # at least 4 (first index is 0)
if max_b_index>=bval.shape[0]: max_b_index=bval.shape[0]-1 # all
Dfast_thresh = 1./(bval[max_b_index-1]*1.)
if Dfast_thresh<3e-3: Dfast_thresh=3e-3   # water is ~2e-3 (50% safety margin)
print ("Using B's",str(bval[0:max_b_index].astype(np.int16))\
       .replace('[','').replace(']','').replace('   ',' ').replace('  ',' '))
print ("Dfast threshold is",Dfast_thresh*1e3)
temp_Afast = np.zeros(shape=temp_IMGdata.shape[0], dtype=np.float32)
temp_Dfast = np.zeros(shape=temp_IMGdata.shape[0], dtype=np.float32)
progress_tag = int(IMGdata_residual.shape[0]/70)
for i in range(IMGdata_residual.shape[0]):
   if i%progress_tag==0: print ('.',end=''); sys.stdout.flush()
   nz = np.nonzero (IMGdata_residual [i,0:max_b_index])     
   if nz[0].shape[0]>=4: #need at least 4 unmasked points
      B_nz = bval[0:max_b_index][nz]
      IMGdata_nz = IMGdata_residual [i,0:max_b_index][nz] 
      initial_conditions=[max(IMGdata_nz), 1./bval[1]]      
      temp_Afast[i], temp_Dfast [i] = FIT (B_nz, IMGdata_nz, initial_conditions)
print (''); sys.stdout.flush()

#calculate Perfusion Fraction
temp_Pfrac = np.zeros(dtype=np.float32, shape=temp_IMGdata.shape[0])
nz = np.nonzero (temp_Afast+temp_Aslow)
temp_Pfrac[nz] = temp_Afast[nz]/(temp_Afast[nz]+temp_Aslow[nz])*100.

#undo vector reduction
data_Dslow = np.zeros((dim[0]*dim[1]*dim[2]),dtype=np.float32)
data_Dfast = np.zeros((dim[0]*dim[1]*dim[2]),dtype=np.float32)
data_Pfrac = np.zeros((dim[0]*dim[1]*dim[2]),dtype=np.float32)
data_Dslow[mask_all>0] = temp_Dslow[:] # undo vector reduction 
data_Dfast[mask_all>0] = temp_Dfast[:] # undo vector reduction
data_Pfrac[mask_all>0] = temp_Pfrac[:] # undo vector reduction     

#reshape to original
data_Dslow  = np.reshape(data_Dslow, (dim[0],dim[1],dim[2]))
data_Dfast  = np.reshape(data_Dfast, (dim[0],dim[1],dim[2]))
data_Pfrac  = np.reshape(data_Pfrac, (dim[0],dim[1],dim[2]))

#clip slow (aka ADC)
data_Dslow [data_Dslow >Dslow_clip]=Dslow_clip

#clip fast 
data_Dfast [data_Dfast <Dfast_thresh]=0
data_Pfrac [data_Dfast <Dfast_thresh]=0

#calculate kernel that serves both 2D and 3D acquisitions
second=np.sort(SpatResol[0:3])[1]; #2nd smalest spatial res 
kernel = np.round(float(second)/SpatResol[0:3],0).astype(np.int) # 2nd smalest element is index 1
kernel_sizes=np.asarray([1,3,5,7])
sigmas = np.round(kernel_sizes/3.*0.7,1); sigmas[kernel_sizes==1]=0
kernel[kernel>len(kernel_sizes)-1]=len(kernel_sizes)-1 # limit to index size
sigma_=np.zeros(shape=kernel.shape, dtype=np.float32)
sigma_[0] = sigmas [kernel[0]] # translate index
sigma_[1] = sigmas [kernel[1]] # translate index
sigma_[2] = sigmas [kernel[2]] # translate index
kernel[0] = kernel_sizes[kernel[0]] # translate index
kernel[1] = kernel_sizes[kernel[1]] # translate index
kernel[2] = kernel_sizes[kernel[2]] # translate index

#median filter 
print ("Median filter results with kernel",str(kernel).replace("[","").replace("]",""))
data_Dslow = median_filter (data_Dslow, size=kernel)
data_Dfast = median_filter (data_Dfast, size=kernel)
data_Pfrac = median_filter (data_Pfrac, size=kernel)

#gaussian filter
print ("Gauss  filter results with sigma ",str(sigma_).replace("[","").replace("]","").replace("0. ","0"))
data_Dslow = gaussian_filter(data_Dslow, sigma=sigma_, truncate=3)
data_Dfast = gaussian_filter(data_Dfast, sigma=sigma_, truncate=3)
data_Pfrac = gaussian_filter(data_Pfrac, sigma=sigma_, truncate=3) 

#clip fast again
data_Dfast [data_Dfast <Dfast_thresh]=0
data_Pfrac [data_Dfast <Dfast_thresh]=0

#convert ADC unit to 10^-3 mm^2/s
#(we suppose that Bs are in mm^2/s)
data_Dslow  *= 1e3 
#usual in-vivo values:
#    white matter: 0.6-0.8
#    grey matter: 0.8-1.0
#    CSF: 3.0-3.5
#    restriction (e.g. tumor, stroke) <0.7 
data_Dfast  *= 1e3 # use same unit
         
#transform to int
max_Dslow  = np.amax(data_Dslow);
data_Dslow *= 32767./max_Dslow
data_Dslow = data_Dslow.astype(np.int16)
max_Dfast  = np.amax(data_Dfast);
data_Dfast *= 32767./max_Dfast
data_Dfast = data_Dfast.astype(np.int16)
max_Pfrac  = np.amax(data_Pfrac);
data_Pfrac *= 32767./max_Pfrac
data_Pfrac = data_Pfrac.astype(np.int16)

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

print ('Successfully written output files'); sys.stdout.flush()

os.system("pause") # windows
