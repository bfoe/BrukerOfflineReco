#
# reads 2 or more NIFTI files
# and calculates the Sum of Squares
# result is saved back to NIFTI format
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
from math import floor
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
img0 = nib.load(FIDfile[0])
data0 = img0.get_data().astype(np.float32)
data0 = np.square(data0)
for i in range (1,nfiles):
    print ('Reading file',i+1)
    img = nib.load(FIDfile[i])
    if img.shape != img0.shape: print ('ERROR: incompatible file dimensions'); sys.exit(2)
    data = img.get_data().astype(np.float32)
    data = np.square(data)
    data0 += data
data=0 # free memory    
data0 = np.sqrt(data0)/np.sqrt(nfiles)
max_ = np.amax(data0)
slope=max_/32767.
data0 /= slope
data0 = data0.astype(np.int16)

#make results folder
dirname = os.path.abspath(os.path.dirname(FIDfile[0])+slash+'..'+slash+'SumOfSquaresResult')
new_dirname = dirname
i=0
while os.path.exists(new_dirname):
   i+=1
   new_dirname = dirname+'('+str(i)+')'
try: os.makedirs(new_dirname)
except: print ('ERROR: unable to make folder', new_dirname); sys.exit(2)
print ("Saving results")
#write averaged fid file
img_SoS = nib.Nifti1Image(data0, img0.get_affine())
img_SoS.header['sform_code']=1
img_SoS.header['qform_code']=1
img_SoS.header.set_slope_inter(slope,0)
new_filename=os.path.basename(FIDfile[0])
new_filename = new_filename[0:new_filename.rfind('_')]+'_SumOfSquares.nii.gz'
nib.save(img_SoS, new_dirname+slash+new_filename)   
#write logfile      
with open(new_dirname+slash+'Logfile.txt', "w") as logfile:
    logfile.write('Result fid file is the complex average of:\n')
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