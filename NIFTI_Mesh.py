#
# This tool reads a previousely binarized image in NIFTI format 
# and extracts a mesh, result written in STL format
#      
#
# ----- VERSION HISTORY -----
#
# Version 0.1 - 1, November 2018
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
#    - vtk
#

from __future__ import print_function
try: import win32gui, win32console
except: pass #silent
import math
import sys
import os
import warnings
import numpy as np
import nibabel as nib
import vtk
from vtk.util import numpy_support


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
InputFile = askopenfilename(title="Choose (masked) NIFTI file", filetypes=[("NIFTI files",('.nii.gz'))])
if InputFile=="":print ('ERROR: No input file specified'); sys.exit(2)
InputFile = os.path.abspath(InputFile)
TKwindows.update()
try: win32gui.SetForegroundWindow(win32console.GetConsoleWindow())
except: pass #silent
dirname  = os.path.dirname(InputFile)
basename = os.path.basename(InputFile)
basename = basename[0:basename.rfind('.nii.gz')]



print ('Reading NIFTI file')

img = nib.load(InputFile)
data = img.get_data().astype(np.float32)
shape = data.shape
affine = img.affine
Origin = np.asarray((affine[0,3],affine[1,3],affine[2,3]))
SpatResol = np.asarray(img.header.get_zooms())

# convert to VTK order
data = np.transpose (data, [2,1,0]).copy() 

#check if binarized
data = data/np.amax(data)
if np.unique(data).shape[0]!=2: print ('ERROR: Inputfile is not binary (containes values != 0,1)'); sys.exit(2)   

# --------- VTK code starts here --------

# convert numpy array to VTK image algorithm output
warnings.filterwarnings("ignore")
VTK_data = numpy_support.numpy_to_vtk(num_array=data.flatten(), deep=True, array_type=vtk.VTK_INT)
img_vtk = vtk.vtkImageData()
img_vtk.SetDimensions(shape)
img_vtk.SetSpacing(SpatResol)
img_vtk.SetOrigin(Origin)
img_vtk.GetPointData().SetScalars(VTK_data)
cast = vtk.vtkImageCast ()
cast.SetInputData(img_vtk)
cast.Update()

#VTK read NIFTI
#reader = vtk.vtkNIFTIImageReader()
#reader.SetFileName(InputFile)
#reader.Update()

print ('Creating Mesh')
close_filter = vtk.vtkImageOpenClose3D()
close_filter.SetInputConnection(cast.GetOutputPort())
close_filter.SetOpenValue(0)
close_filter.SetCloseValue(1)
close_filter.SetKernelSize(3,3,3)
close_filter.Update()

mesh = vtk.vtkDiscreteMarchingCubes()
mesh.SetInputConnection(close_filter.GetOutputPort())
mesh.GenerateValues(1, 1, 1)
mesh.Update()

#only largest connected
conn = vtk.vtkPolyDataConnectivityFilter()
conn.SetInputConnection(mesh.GetOutputPort())
conn.SetExtractionModeToLargestRegion()
conn.Update()

#smooth
smooth = vtk.vtkSmoothPolyDataFilter()
smooth.SetInputConnection(conn.GetOutputPort())
smooth.SetNumberOfIterations(2)
smooth.SetRelaxationFactor(0.5)
smooth.FeatureEdgeSmoothingOff()
smooth.BoundarySmoothingOff()
normals = vtk.vtkPolyDataNormals()
normals.SetInputConnection(smooth.GetOutputPort())
normals.FlipNormalsOn()
mapper = vtk.vtkPolyDataMapper()
mapper.SetInputConnection(normals.GetOutputPort())
actor = vtk.vtkActor()
actor.SetMapper(mapper)
actor.GetProperty().SetInterpolationToFlat()
smooth.Update()

#write STL
writer = vtk.vtkSTLWriter()
writer.SetInputConnection(smooth.GetOutputPort())
writer.SetFileTypeToBinary()
writer.SetFileName(basename+".stl")
writer.Write()       


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