#
# This tool reads a NIFTI image, analyzes the histogram to etablish a threshold that includes
# 2% of the total amount of voxels with the highest signal intensity and filters out any
# isolated clusters leaving only the single largest cluster of interconnected voxels.
# The resulting image is again written in NIFTI format adding "conected" to the filename.
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
if getattr( sys, 'frozen', False ): # running as pyinstaller bundle
   from scipy_extract import label   
else: # running native python
   from scipy.ndimage import label 


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
    
#intercatively choose input NIFTI file
InputFile = askopenfilename(title="Choose (masked) NIFTI file", filetypes=[("NIFTI files",('MAG_m.nii.gz','MAGNT_masked.nii.gz'))])
if InputFile=="":print ('ERROR: No input file specified'); sys.exit(2)
InputFile = os.path.abspath(InputFile)
TKwindows.update()
try: win32gui.SetForegroundWindow(win32console.GetConsoleWindow())
except: pass #silent
dirname  = os.path.dirname(InputFile)
basename = os.path.basename(InputFile)
filename = os.path.splitext(InputFile)[0]
extension = os.path.splitext(InputFile)[1].lower()
filename = filename[0:filename.rfind('MAG_m')]
filename = filename[0:filename.rfind('_MAGNT_masked')]
filename_connected = filename+'_Connected.nii.gz'


print ('Reading NIFTI file')
img = nib.load(InputFile)
data = img.get_data().astype(np.float32)
SpatResol = np.asarray(img.header.get_zooms())


print ('Analyzing histogram')
#find largest extension of the volume)
directions = np.argsort(data.shape)
directions = directions[::-1] # decreasing order
if   directions[0]==0: d=data[data.shape[0]/4:data.shape[0]*3/4,:,:]
elif directions[0]==1: d=data[:,data.shape[1]/4:data.shape[1]*3/4,:]
else:                  d=data[:,:,data.shape[2]/4:data.shape[2]*3/4]
#calculate histogram
n_points=d.shape[0]*d.shape[1]*d.shape[2]
steps=int(n_points/10000); start=0; fin=np.max(data [:,:,:])
xbins =  np.linspace(start,fin,steps)
ybins, binedges = np.histogram(d[:,:,:], bins=xbins)
ybins = np.resize (ybins,len(xbins)); ybins[len(ybins)-1]=0
xbins = xbins[1:-1] # through out first point (lots of counts)
ybins = ybins[1:-1] # through out first point (lots of counts)
nz = np.nonzero (ybins)
xbins = xbins[nz] # through out histogram points with zero count (basically the initial points)
ybins = ybins[nz] # through out histogram points with zero count (basically the initial points)
ybins = smooth(ybins,ybins.shape[0]/20)
##find minimum
#start=ybins.argmax()
#i=start;x_min=0;y_min=ybins[start]
#while i<len(ybins):
#    i+=1
#    if ybins[i]<=y_min: y_min=ybins[i]; x_min=i; 
#    else: i=len(ybins);
#minimum=xbins[x_min]
#print ('Minimum   =',minimum)
##find maximum
#start=x_min
#i=start;x_max=0;y_max=ybins[start]
#while i<len(ybins):
#    i+=1
#    if ybins[i]>y_max: y_max=ybins[i]; x_max=i; 
#    else: i=len(ybins);
#maximum=xbins[x_max]
#print ('Maximum   =',maximum)
#threshold=maximum/3.
i=len(ybins)-1
points=0; threshold = 0
while i>0:
   points += ybins[i]
   threshold = xbins[i]   
   if points >  n_points*0.02: #2%
      i=1
   i -= 1      
print ('Threshold = %0.2f' % threshold) 
mask =  data [:,:,:] > threshold
mask = mask.astype(np.int16)
##write histogram results
#with open(os.path.join(dirname,filename+'_Histogram.txt'), "w") as text_file:    
#    for i in range(0,xbins.shape[0]):     
#        text_file.write("%e\t %e\n" % (xbins[i], ybins[i]))
#    text_file.write("\n")      

# clear up mask: leave only the largest cluster of connected points
print ('Extracting larges connected cluster')
s = [[[1,1,1],[1,1,1],[1,1,1]], [[1,1,1],[1,1,1],[1,1,1]], [[1,1,1],[1,1,1],[1,1,1]]]
labeled_mask, num_clusters = label(mask, structure=s)
unique, counts = np.unique(labeled_mask, return_counts=True)
max_count=0
for i in range(0,unique.shape[0]): # find the largest nonzero count
    if counts[i]>max_count and unique[i]!=0: max_count=counts[i]
remove_labels = unique[np.where(counts<max_count)] # leave only the largest cluster of connected points
remove_indices = np.where(np.isin(labeled_mask,remove_labels))
mask[remove_indices] = 0    
data2=data*mask
    
#transform to int 
max_data = np.amax(data);
data *= 32767./max_data
data = data.astype(np.int16)
max_data2 = np.amax(data2);
data2 *= 32767./max_data2
data2 = data2.astype(np.int16)
#save NIFTI's
print ('Writing output files')
aff = np.eye(4)
aff[0,0] = SpatResol[0]; aff[0,3] = -(data.shape[0]/2)*aff[0,0]
aff[1,1] = SpatResol[1]; aff[1,3] = -(data.shape[1]/2)*aff[1,1]
aff[2,2] = SpatResol[2]; aff[2,3] = -(data.shape[2]/2)*aff[2,2]
NIFTIimg = nib.Nifti1Image(data2, aff)
NIFTIimg.header.set_slope_inter(max_data2/32767.,0)
NIFTIimg.header.set_xyzt_units(3, 8)
NIFTIimg.set_sform(aff, code=0)
NIFTIimg.set_qform(aff, code=1)
try:
    nib.save(NIFTIimg, os.path.join(dirname,filename_connected))
except:
    print ('\nERROR:  problem while writing results'); sys.exit(1)
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