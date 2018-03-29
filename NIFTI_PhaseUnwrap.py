#
# reads 2 or more NIFTI files that contain 
# phase data from a PCA acqusition with different venc values 
# 
# information from the higer venc values is used to
# correct phase wraps in the lower venc images
# result is saved back to NIFTI format
#
# ----- VERSION HISTORY -----
#
# Version 0.1 - 30, March 2018
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
#


from __future__ import print_function
try: import win32gui, win32console
except: pass #silent
import math
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

def combine_arrays(arrays):
    #
    # Generate a cartesian product of input arrays.
    #   Input:   list of 1-D arrays.
    #   Returns: 2-D array of shape (M, len(arrays)) containing cartesian product.
    #
    # Example
    # >>> combine_arrays(([1, 2, 3], [4, 5], [6, 7]))
    # array([[1, 4, 6],[1, 4, 7],[1, 5, 6],[1, 5, 7],[2, 4, 6],[2, 4, 7],
    #        [2, 5, 6],[2, 5, 7],[3, 4, 6],[3, 4, 7],[3, 5, 6],[3, 5, 7]])
    #
    arrays = [np.asarray(a) for a in arrays]
    shape = (len(x) for x in arrays)
    ix = np.indices(shape, dtype=int)
    ix = ix.reshape(len(arrays), -1).T
    for n, arr in enumerate(arrays):
        ix[:, n] = arrays[n][ix[:, n]]
    return ix    
    
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
    
#intercatively choose input FID files
nfiles=0
answer="dummy"
FIDfile=np.array([])
while answer!="":
   answer = askopenfilename(title="Choose NIFTI file "+str(nfiles+1)+" (press cancel to end)", filetypes=[("NIFTI files",('*.nii','*.NII','*.nii.gz','*.NII.GZ'))])
   if answer!="":
        answer = os.path.abspath(answer)
        FIDfile = np.append(FIDfile, answer)
        nfiles+=1
if nfiles==0: print ('ERROR: No input file specified'); sys.exit(2)
if nfiles==1: print ('ERROR: Need at least 2 files'); sys.exit(2)
TKwindows.update()
try: win32gui.SetForegroundWindow(win32console.GetConsoleWindow())
except: pass #silent

print ('Reading file 1')
img = nib.load(FIDfile[0])
shape0 = img.shape
IMGdata_decoded_PH = img.get_data().astype(np.float32)
IMGdata_decoded_PH = np.expand_dims(IMGdata_decoded_PH, axis=(3))
for i in range (1,nfiles):
    print ('Reading file',i+1)
    img = nib.load(FIDfile[i])
    if img.shape != shape0: print ('ERROR: incompatible file dimensions'); sys.exit(2)
    IMGdata_decoded_PH = np.append(IMGdata_decoded_PH, np.expand_dims(img.get_data().astype(np.float32), axis=(3)),axis=(3))

print ('Starting Phase unwrap')
#calc parameters for phase unwrap
mask =  IMGdata_decoded_PH [:,:,:,:] > 0
venc_arr = np.amax (IMGdata_decoded_PH[:,:,:,:], axis=(0,1,2))
max_venc = np.amax(venc_arr)
niter_arr = np.zeros (shape=(nfiles),dtype=np.int32)
list_of_arrays = []
for file_nr in range (0,nfiles): 
    niter_arr[file_nr] = math.ceil(max_venc/venc_arr[file_nr])
    list_of_arrays.append (np.asarray(range (-1*niter_arr[file_nr]+1, niter_arr[file_nr])))
PhUnwrap_combinations = combine_arrays (list_of_arrays)
i=0
while i<PhUnwrap_combinations.shape[0]:
    signs = np.unique (np.sign(PhUnwrap_combinations[i,:]))
    if (1 in signs and -1 in signs):
        PhUnwrap_combinations = np.delete (PhUnwrap_combinations, i, axis=(0))
    else: i +=1
PhUnwrap_combinations = PhUnwrap_combinations.astype (float)   
n_combinations = PhUnwrap_combinations.shape[0]
for file_nr in range (0,nfiles):
    PhUnwrap_combinations[:,file_nr] *= venc_arr[file_nr]
# force combinations increasing with decreasing venc
i=0
while i<PhUnwrap_combinations.shape[0]:
    remove = False
    for file_nr in range (0,nfiles-1):
        if venc_arr[file_nr] > venc_arr[file_nr+1]:
            if np.abs(PhUnwrap_combinations[i,file_nr])>np.abs(PhUnwrap_combinations[i,file_nr+1]): remove = True
        else:
            if np.abs(PhUnwrap_combinations[i,file_nr])<np.abs(PhUnwrap_combinations[i,file_nr+1]): remove = True
    if remove: PhUnwrap_combinations = np.delete (PhUnwrap_combinations, i, axis=(0))
    else: i +=1
    
#PhaseUnwrap
Img_PH_flow = np.zeros (shape=(IMGdata_decoded_PH.shape[0],IMGdata_decoded_PH.shape[1],IMGdata_decoded_PH.shape[2]),dtype=np.float32)  
for x in range (0, Img_PH_flow.shape[0]):
   if x%4 == 0: print('.', end='') #progress indicator
   for y in range (0, Img_PH_flow.shape[1]):
     for z in range (0, Img_PH_flow.shape[2]):
        nz = np.nonzero (mask[x,y,z,:].astype(int))
        local_nfiles = np.asarray(nz).shape[1]
        if   local_nfiles == 0: Img_PH_flow [x,y,z] = 0
        elif local_nfiles == 1: IMGdata_decoded_PH[x,y,z,nz]
        else:      
           local_combinations=np.unique(PhUnwrap_combinations[:,nz], axis=(0))
           local_combinations=local_combinations.reshape((local_combinations.shape[0],local_nfiles)) #eliminate an empty dimension
           local_max_venc = np.amax(venc_arr[nz])
           index_local_max_venc = np.argmax(venc_arr[nz])
           i=0
           while i<local_combinations.shape[0]:
               if local_combinations[i,index_local_max_venc] != 0:
                   local_combinations = np.delete (local_combinations, i, axis=(0))
               else: i +=1              
           local_n_combinations=local_combinations.shape[0]
           PhUnwrap_try = np.zeros (shape=(local_n_combinations,local_nfiles),dtype=np.float32)           
           for i in range (0,local_n_combinations):
              PhUnwrap_try[i,:] = IMGdata_decoded_PH[x,y,z,nz]+local_combinations[i,:]            
           avg = np.average (PhUnwrap_try[:,:],axis=(1))
           chi = np.sqrt(np.sum(np.square(PhUnwrap_try[:,:]-avg[:,None]),axis=(1)))/local_nfiles   
           min_chi_index = np.argmin (chi)
           Img_PH_flow [x,y,z] = avg[min_chi_index]
print (' ')

#transform to int
Max_PH_flow   = np.amax (np.abs(Img_PH_flow));
Img_PH_flow *= 32767./Max_PH_flow
Img_PH_flow   = Img_PH_flow.astype(np.int16)


#make results folder
dirname = os.path.abspath(os.path.dirname(FIDfile[0])+slash+'..'+slash+'PhaseUnwrap')
new_dirname = dirname
i=0
while os.path.exists(new_dirname):
   i+=1
   new_dirname = dirname+'('+str(i)+')'
try: os.makedirs(new_dirname)
except: print ('ERROR: unable to make folder', new_dirname); sys.exit(2)
print ("Saving results")
#write averaged fid file
aff = img.get_affine()
unit_xyz, unit_t = img.header.get_xyzt_units()
if unit_xyz == 'unknown': unit_xyz=0
if unit_t   == 'unknown': unit_t=0
img_SoS = nib.Nifti1Image(Img_PH_flow, aff)
img_SoS.header.set_xyzt_units(unit_xyz, unit_t)
img_SoS.set_sform(aff, code=0)
img_SoS.set_qform(aff, code=1)
img_SoS.header.set_slope_inter(Max_PH_flow/32767.,0)
new_filename=os.path.basename(FIDfile[0])
new_filename = new_filename[0:new_filename.rfind('_')]+'_PhaseUnwrap.nii.gz'
nib.save(img_SoS, new_dirname+slash+new_filename)   
#write logfile      
with open(new_dirname+slash+'Logfile.txt', "w") as logfile:
    logfile.write('Output NIFTI file is the phase unwrapping result from:\n')
    for i in range (0,nfiles):
       logfile.write(FIDfile[i]+'\n')
print ("done\n")  
     
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