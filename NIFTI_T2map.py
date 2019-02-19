#
# calculates maps of T2 decay
#      
#
# ----- VERSION HISTORY -----
#
# Version 0.1 - 20, February 2019
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
import binascii
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
    
#intercatively choose input NIFTI file
InputFile = askopenfilename(title="Choose NIFTI file", filetypes=[("NIFTI files",('*.nii','*.NII','*.nii.gz','*.NII.GZ'))])
if InputFile=="":print ('ERROR: No input file specified'); sys.exit(2)
InputFile = os.path.abspath(InputFile)
TKwindows.update()
try: win32gui.SetForegroundWindow(win32console.GetConsoleWindow())
except: pass #silent
dirname  = os.path.dirname(InputFile)
basename = os.path.basename(InputFile)
filename = os.path.splitext(InputFile)[0]
extension = os.path.splitext(InputFile)[1].lower()


print ('Reading NIFTI file')
img = nib.load(InputFile)
IMGdata = img.get_data().astype(np.float32)
SpatResol = np.asarray(img.header.get_zooms())
TE_str = img.header['descrip'] 
try: TEs = np.fromstring(binascii.a2b_uu(str(TE_str)[6:]),dtype=np.int16)/10.
except: print ("ERROR: parsing TE information from header"); sys.exit(1)

#check some stuff
if len(IMGdata.shape) != 4:
    print ('ERROR: 4D NIFTI expected, 3D NIFTI found, this is probably not a multiecho file'); 
    sys.exit(1)
if TEs.shape[0]!=IMGdata.shape[3]:
    print ("ERROR: number of TE's in header unequal data dimension "); 
    sys.exit(1)
if np.unique(TEs).shape[0]<7:
    print ("ERROR: need at least 7 unique TE's"); 
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
if np.count_nonzero(arr==0)<0.2*arr.size: 
    # use noise in all 8 corners to establish threshold
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
    IMGdata = IMGdata*mask
    
print ('number of  zero points', np.count_nonzero(IMGdata==0))
print ('total number of points', IMGdata.size)
print ('data shape', IMGdata.shape)    
print ("TE's = ",TEs)


    
#do stuff here


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