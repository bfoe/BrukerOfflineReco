#
# calculates flow volume in [ml/s] along the main flow direction
# (for coronal H_F acquisitions the X component)
#      
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
#

from __future__ import print_function
try: import win32gui, win32console
except: pass #silent
import math
import sys
import os
import numpy as np
import nibabel as nib
'''
if getattr( sys, 'frozen', False ): # running as pyinstaller bundle
   from scipy_extract import label   
else: # running native python
   from scipy.ndimage import label 
'''

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
    
#intercatively choose input NIFTI file
InputFile = askopenfilename(title="Choose NIFTI file", filetypes=[("NIFTI files",('*.nii','*.NII','*.nii.gz','*.NII.GZ'))])
if InputFile=="":print ('ERROR: No input file specified'); sys.exit(2)
InputFile = os.path.abspath(InputFile)
TKwindows.update()
try: win32gui.SetForegroundWindow(win32console.GetConsoleWindow())
except: pass #silent

print ('Reading NIFTI file')
img = nib.load(InputFile)
data = img.get_data().astype(np.float32)
SpatResol = img.header.get_zooms()
'''
#mask out low flow values
threshold=np.max(data)*0.01 # 1%
mask =  data [:,:,:] > threshold

# clear up mask: leave only the largest cluster of connected points
print ('Calculating ...')
s = [[[1,1,1],[1,1,1],[1,1,1]], [[1,1,1],[1,1,1],[1,1,1]], [[1,1,1],[1,1,1],[1,1,1]]]
labeled_mask, num_clusters = label(mask, structure=s)
unique, counts = np.unique(labeled_mask, return_counts=True)
max_count=0
for i in range(0,unique.shape[0]): # find the largest nonzero count
    if counts[i]>max_count and unique[i]!=0: max_count=counts[i]
remove_labels = unique[np.where(counts<max_count)] # leave only the largest cluster of connected points
remove_indices = np.where(np.isin(labeled_mask,remove_labels))
mask[remove_indices] = 0

# apply mask
data[:,:,:] *= mask [:,:,:]
'''

#write flowvolume results
print ('Writing flow rate file')
dirname  = os.path.dirname(InputFile)
basename = os.path.basename(InputFile)
filename = os.path.splitext(InputFile)[0]
with open(os.path.join(dirname,filename+'_FlowVolumes.txt'), "w") as text_file:    
    text_file.write("Flow Volumes per slice (X):\n")
    for i in range(0,data.shape[0]): # in our data shape[2] is the main flow direction
        flowvol = np.sum(data[i,:,:])
        flowvol = abs(flowvol)
        flowvol *= 10. # venc is in cm/s, multiply by 10. to get this in mm/s
        flowvol *= SpatResol[1]/1000.*SpatResol[2]/1000. # multiply with inplane spatial resolution, result is in mm^3/s
        flowvol /= 1000. # convert mm^3/s ---> ml/s
        text_file.write("Slice %d:\t%0.2f\tml/s\n" % (i, flowvol))
    text_file.write("\n")        


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