#
# reduces the resolution of a NIFTI image by a factor of 2
#
# ----- VERSION HISTORY -----
#
# Version 0.1 - 27, September 2018
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
import sys
import os
import numpy as np
import nibabel as nib

pywin32_installed=True
try: import win32console, win32gui, win32con
except: pywin32_installed=True

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
    exit(2)

       
#general initialization stuff
space=' '; slash='/'; 
if sys.platform=="win32": slash='\\' # not really needed, but looks nicer ;)
Program_name = os.path.basename(sys.argv[0]); 
if Program_name.find('.')>0: Program_name = Program_name[:Program_name.find('.')]
python_version=str(sys.version_info[0])+'.'+str(sys.version_info[1])+'.'+str(sys.version_info[2])
# sys.platform = [linux2, win32, cygwin, darwin, os2, os2emx, riscos, atheos, freebsd7, freebsd8]
if sys.platform=="win32": 
    os.system("title "+Program_name)    
   
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

#interactively choose NIFTI file
FIDfile1 = askopenfilename(title="Choose moving NIFTI file", filetypes=[("NIFTI files",('*.nii','*.NII','*.nii.gz','*.NII.GZ'))])
if FIDfile1 == "": print ('ERROR: NIFTI file not specified'); exit(2)
FIDfile1 = os.path.abspath(FIDfile1) 


# Set Output filenames
dirname  = os.path.dirname(FIDfile1)
basename = os.path.basename(FIDfile1)
basename = basename[0:basename.rfind('.nii.gz')]
outfile = basename+'_downsampled.nii.gz'


#read NIFTI file
print ('Read  NIFTI file')
img = nib.load(FIDfile1)
data = img.get_data().astype(np.float32)
SpatResol = np.asarray(img.header.get_zooms())
Shape = np.asarray(img.header.get_data_shape())
affine = img.affine
xyzt_units1  = img.header.get_xyzt_units()[0]
xyzt_units2  = img.header.get_xyzt_units()[1]
sform = int(img.header['sform_code'])
qform = int(img.header['qform_code'])

#downsample
data = data[::2,::2,::2]
affine [:,0:3] *= 2

#rewrite NIFTI
#transform to int
max_ = np.amax(data)
data *= 32767./max_
data = data.astype(np.int16)
#Write NIFTI
img = nib.Nifti1Image(data, affine)
img.header.set_slope_inter(max_/32767.,0)
img.header.set_xyzt_units(xyzt_units1,xyzt_units2)
img.set_sform(affine, sform)
img.set_sform(affine, qform)
try: nib.save(img, os.path.join(dirname,outfile))
except: print ('ERROR:  problem while writing results'); exit(1)
print ('Successfully written output file '+outfile)
      
#end
exit(0)

