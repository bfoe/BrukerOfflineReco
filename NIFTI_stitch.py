#
# This tool reads two NIFTI files and stiches them together
#      
#
# ----- VERSION HISTORY -----
#
# Version 0.1 - 08, August 2018
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
    
#intercatively choose input NIFTI files
InputFile1 = askopenfilename(title="Choose 1st NIFTI file", filetypes=[("NIFTI files",('*.nii','*.NII','*.nii.gz','*.NII.GZ'))])
if InputFile1=="":print ('ERROR: No input file specified'); sys.exit(2)
InputFile1 = os.path.abspath(InputFile1)
InputFile2 = askopenfilename(title="Choose 2nd NIFTI file", filetypes=[("NIFTI files",('*.nii','*.NII','*.nii.gz','*.NII.GZ'))])
if InputFile2=="":print ('ERROR: No input file specified'); sys.exit(2)
InputFile2 = os.path.abspath(InputFile2)
TKwindows.update()
try: win32gui.SetForegroundWindow(win32console.GetConsoleWindow())
except: pass #silent
dirname  = os.path.dirname(InputFile1)
basename = os.path.basename(InputFile1)
filename = os.path.splitext(InputFile1)[0]
extension = os.path.splitext(InputFile1)[1].lower()
filename = filename[0:filename.rfind('*.nii.gz')]
filename_connected = filename+'_stitched.nii.gz'


print ('Reading NIFTI files')
img1 = nib.load(InputFile1)
data1 = img1.get_data().astype(np.float32)
SpatResol1 = np.asarray(img1.header.get_zooms())
img2 = nib.load(InputFile2)
data2 = img2.get_data().astype(np.float32)
SpatResol2 = np.asarray(img2.header.get_zooms())

#some checks
if not np.array_equal(data1.shape,data2.shape): print ('ERROR: image dimension mismatch'); sys.exit(2)
if not np.array_equal(SpatResol1,SpatResol2): print ('ERROR: image resolution mismatch'); sys.exit(2)
directions = np.argsort(data1.shape)
directions = directions[::-1] # decreasing order
if directions[0]!=1: print ('ERROR: largest dimension is not index 1, not implemented'); sys.exit(2)

#initialize match search
overlap=int(0.05*data1.shape[1]) #5%
overlap=int(overlap/2)*2 # make it even
print ('Overlap is ',overlap)
stitch_search_range = int(0.20*data1.shape[1]) #25%
print ('Shift search range is 0..'+str(2*stitch_search_range))
roll1_search_range  = int(0.05*data1.shape[0]) #5%
print ('Roll1 search range is -'+str(roll1_search_range)+'..'+str(roll1_search_range))
roll2_search_range  = int(0.05*data1.shape[0]) #5%
print ('Roll2 search range is -'+str(roll2_search_range)+'..'+str(roll2_search_range))
goodness = np.zeros ((2*stitch_search_range, 2*roll1_search_range+1, 2*roll2_search_range+1),dtype=np.float64)
#find match
print ('optimizing ',end='')
for i in range (0,stitch_search_range):
   print ('.',end='')
   for j in range (0,2*roll1_search_range+1):
       for k in range (0,2*roll2_search_range+1):
           #0
           test1 = data1[:,data1.shape[1]-i-overlap:data1.shape[1]-i,:]          
           test2 = data2[:,i:i+overlap,:]
           test2 = np.roll(test2,j-roll1_search_range,axis=0)
           test2 = np.roll(test2,k-roll2_search_range,axis=2)
           test1=test1.flatten().astype(np.float64); test2=test2.flatten().astype(np.float64)
           goodness [2*i,j,k] = np.correlate(test1,test2)
           #1
           test1 = data1[:,data1.shape[1]-i-overlap:data1.shape[1]-i,:]          
           test2 = data2[:,i+1:i+overlap+1,:] # add 1
           test2 = np.roll(test2,j-roll1_search_range,axis=0)
           test2 = np.roll(test2,k-roll2_search_range,axis=2)          
           test1=test1.flatten().astype(np.float64); test2=test2.flatten().astype(np.float64)
           goodness [2*i+1,k] = np.correlate(test1,test2)           
print ('')
max_=0; max_i=0; max_j=0; max_k=0
for i in range (0,2*stitch_search_range):
   for j in range (0,2*roll1_search_range+1):
       for k in range (0,2*roll2_search_range+1):
           if goodness [i,j,k]>max_: 
              max_ = goodness [i,j,k]
              max_i=i
              max_j=j-roll1_search_range 
              max_k=k-roll2_search_range
print ('Optimal shift found at',max_i)
print ('Optimal roll1 found at',max_j)
print ('Optimal roll2 found at',max_k)

              
#delete some points
#data1 = data1[:,0:data1.shape[1]+1,:] #full
#data2 = data2[:,0:data1.shape[1]+1,:] #full
#data1 = data1[:,0:data1.shape[1]+1-78,:] #OK
#data2 = data2[:,50:data2.shape[1]+1,:]   #OK 
#data2 = np.roll (data2,2, axis=0)        #OK 
#data2 = np.roll (data2,5, axis=2)        #OK
crop1 = int(max_i/2)
crop2 = max_i-crop1
crop1 += int(overlap/2)
crop2 += int(overlap/2)
#data1 = data1[:,0:data1.shape[1]+1-crop1,:] #looks worse
data1 = data1[:,0:data1.shape[1]-crop1,:] #looks better
data2 = data2[:,crop2:data2.shape[1]+1,:]   
data2 = np.roll (data2,max_j, axis=0)
data2 = np.roll (data2,max_k, axis=2)

#concat
data = np.concatenate ((data1, data2), axis=1)
    
#transform to int 
max_data = np.amax(data);
data *= 32767./max_data
data = data.astype(np.int16)
#save NIFTI
print ('Writing output files')
aff = np.eye(4)
aff[0,0] = SpatResol1[0]*1000; aff[0,3] = -(data.shape[0]/2)*aff[0,0]
aff[1,1] = SpatResol1[1]*1000; aff[1,3] = -(data.shape[1]/2)*aff[1,1]
aff[2,2] = SpatResol1[2]*1000; aff[2,3] = -(data.shape[2]/2)*aff[2,2]
NIFTIimg = nib.Nifti1Image(data, aff)
NIFTIimg.header.set_slope_inter(max_data/32767.,0)
NIFTIimg.header.set_xyzt_units(3, 8)
NIFTIimg.set_sform(aff, code=0)
NIFTIimg.set_qform(aff, code=1)
try:
    nib.save(NIFTIimg, os.path.join(dirname,filename_connected))
except:
    print ('\nERROR:  problem while writing results'); sys.exit(1)
print ('If the result looks wrong, try choosing the input files in inverse order')   
print ('\nSuccessfully written output file')   
        


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