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
import vtk
from vtk.util import numpy_support
import multiprocessing as mp

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

def redirect_vtk_messages ():
    ow = vtk.vtkOutputWindow()
    ow.SendToStdErrOn()
    #1
    #errOut = vtk.vtkFileOutputWindow()
    #errOut.SetFileName(os.path.join(dirname,'VTK_errors.log')
    #vtkStdErrOut = vtk.vtkOutputWindow()
    #vtkStdErrOut.SetInstance(errOut)
    #2
    #log = vtk.vtkFileOutputWindow()
    #log.SetFlush(1)
    #log.SetFileName(os.path.join(dirname,'VTK_errors.log'))
    #log.SetInstance(log)         
    
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
   

# --------- VTK code starts here --------
redirect_vtk_messages ()

#set number of cores   
cores=mp.cpu_count()-1; cores = max (1,cores)
print ('Vtk multithreading set to %d cores ' % cores)   
vtk.vtkMultiThreader.SetGlobalMaximumNumberOfThreads(cores)

#VTK read NIFTI
print ('Reading NIFTI image')
reader = vtk.vtkNIFTIImageReader()
reader.SetFileName(InputFile)
reader.Update()

#convert to float
tofloat = vtk.vtkImageCast()
tofloat.SetInputConnection(reader.GetOutputPort())
tofloat.SetOutputScalarTypeToFloat()
tofloat.Update()

#check binary
vtk_data = tofloat.GetOutput().GetPointData().GetScalars()
numpy_data = numpy_support.vtk_to_numpy(vtk_data)
u = np.unique(numpy_data).shape[0]
if u!=2: 
   print ('ERROR: Inputfile is not binary')
   sys.exit(2)
del u; del numpy_data; del vtk_data # free memory

#get current resolution and set interpolation factor
spacing  = tofloat.GetOutput().GetSpacing()
interpolation = 2.0
spacingX = spacing[0]/interpolation
spacingY = spacing[1]/interpolation 
spacingZ = spacing[2]/interpolation

print ('Interpolating image')
interp = vtk.vtkImageSincInterpolator()
interp.SetWindowFunctionToLanczos()
resize = vtk.vtkImageReslice()
resize.SetOutputSpacing((spacingX, spacingY, spacingZ))
resize.SetInterpolator(interp)
resize.SetInputConnection(tofloat.GetOutputPort())
resize.Update()
del tofloat # free memory

print ('Smoothing image')
imagesmooth = vtk.vtkImageGaussianSmooth()
imagesmooth.SetInputConnection(resize.GetOutputPort())
imagesmooth.SetRadiusFactors(1.2,1.2,1.2) # very light smoothing
imagesmooth.Update()
del resize # free memory

#normalize to 1
normalize = vtk.vtkImageMathematics()
normalize.SetInputConnection(imagesmooth.GetOutputPort())
normalize.SetOperationToMultiplyByK()
normalize.SetConstantK(1./imagesmooth.GetOutput().GetScalarRange()[1])
normalize.Update()
del imagesmooth # free memory

print ('Thresholding image')
thresh = vtk.vtkImageThreshold()
thresh.SetInputConnection(normalize.GetOutputPort())
thresh.ThresholdByUpper(0.45) # 0.5 theoreticaly
thresh.SetInValue(1)
thresh.SetOutValue(0)
thresh.SetOutputScalarTypeToShort()
thresh.Update()
del normalize # free memory

#convert to short
toshort = vtk.vtkImageCast()
toshort.SetInputConnection(thresh.GetOutputPort())
toshort.SetOutputScalarTypeToShort()
toshort.Update()
del thresh # free memory

print ('Apllying closing filter')
close_filter = vtk.vtkImageOpenClose3D()
close_filter.SetInputConnection(toshort.GetOutputPort())
close_filter.SetOpenValue(0)
close_filter.SetCloseValue(1)
close_filter.SetKernelSize(7,7,7)
close_filter.Update()
del toshort # free memory

print ('Creating Mesh    ',end='')
mesh = vtk.vtkDiscreteMarchingCubes()
mesh.SetInputConnection(close_filter.GetOutputPort())
mesh.GenerateValues(1, 1, 1)
mesh.Update()
ncells = mesh.GetOutput().GetNumberOfCells()
if ncells>0: print(" --> Found %d cells" % ncells)
else: print("\nERROR: zero cells in mesh"); sys.exit(1)
del close_filter # free memory

print ('Removing isolated',end='')
conn = vtk.vtkPolyDataConnectivityFilter()
conn.SetInputConnection(mesh.GetOutputPort())
conn.SetExtractionModeToLargestRegion()
conn.Update()
ncells = conn.GetOutput().GetNumberOfCells()
if ncells>0: print(" --> Remaining %d cells" % ncells)
else: print("\nERROR: zero cells in mesh"); sys.exit(1)
del mesh # free memory

smooth1 = vtk.vtkSmoothPolyDataFilter()
smooth1.SetInputConnection(conn.GetOutputPort())
smooth1.SetNumberOfIterations(7)
smooth1.SetRelaxationFactor(0.15)
smooth1.FeatureEdgeSmoothingOff()
smooth1.BoundarySmoothingOff()
smooth1.Update()
del conn # free memory

print ('Decimating mesh  ',end='')
reduction = 1.-1./interpolation
decimate = vtk.vtkQuadricDecimation()
decimate.SetInputConnection(smooth1.GetOutputPort())
decimate.SetTargetReduction(reduction)
decimate.Update()
ncells = decimate.GetOutput().GetNumberOfCells()
if ncells>0: print(" --> Remaining %d cells" % ncells)
else: print("\nERROR: zero cells in mesh"); sys.exit(1)
del smooth1 # free memory

smooth2 = vtk.vtkSmoothPolyDataFilter()
smooth2.SetInputConnection(decimate.GetOutputPort())
smooth2.SetNumberOfIterations(7)
smooth2.SetRelaxationFactor(0.15)
smooth2.FeatureEdgeSmoothingOff()
smooth2.BoundarySmoothingOff()
smooth2.Update()
del decimate # free memory

writer = vtk.vtkSTLWriter()
writer.SetInputConnection(smooth2.GetOutputPort())
writer.SetFileTypeToBinary()
writer.SetFileName(os.path.join(dirname,basename+'.stl'))
writer.Write()       
del smooth2 # free memory
del writer  # free memory

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