#
# calculates the histogram of voxels along 50% of the largest image extension
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
import zlib
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

#read NIFTI
print ('Reading NIFTI file')
img = nib.load(InputFile)
data = img.get_data().astype(np.float32)
SpatResol = np.asarray(img.header.get_zooms())

#find main flow direction (suposed to be the largest extension of the volume)
directions = np.argsort(data.shape)
directions = directions[::-1] # decreasing order
#change order according to direction with larges extension
data = np.transpose(data, directions)
SpatResol = SpatResol [directions]

#calculate histogram
FOV_start = int(1/4. * data.shape[0])
FOV_end   = int(3/4. * data.shape[0])
n_points=np.prod(data.shape)/20
steps =np.sqrt(n_points) 
start = np.min(data[np.nonzero(data[:,:,:])])
fin   = np.max(data[FOV_start:FOV_end,:,:])
xbins =  np.linspace(start,fin,steps)
ybins, binedges = np.histogram(data[FOV_start:FOV_end,:,:], bins=xbins)
ybins = np.resize (ybins,len(xbins)); ybins[len(ybins)-1]=0
ybins_ln = np.zeros(ybins.shape[0],dtype=np.float32)
ybins_ln[ybins>0] = np.log(ybins[ybins>0])
ybins = smooth(ybins,21)
ybins_ln = smooth(ybins_ln,21)

#write histogram results
print ('Writing histogram file')
with open(os.path.join(dirname,filename+'_Histogram.txt'), "w") as text_file:    
    for i in range(0,xbins.shape[0]):     
        text_file.write("%e\t %e\t %e\n" % (xbins[i], ybins[i], ybins_ln[i]))
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