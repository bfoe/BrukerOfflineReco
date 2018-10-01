#
# registers two NIFTI images using elastix (http://elastix.isi.uu.nl/)
# with the following strategy:
#  - Metric "AdvancedMattesMutualInformation"
#  - Optimizer "AdaptiveStochasticGradientDescent"
#  - Transforms: 
#    1) permutation of dimensions
#    2) rough estimation of rotation along largest dimension (rotx)
#    3) Transform "EulerTransform"
#    4) Transform "BSplineTransform"
#
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
#    - elastix executables:
#       1) elastix.exe 
#       2) transformix.exe
#       3) ANNlib-4.9.dll
#      from:
#      https://github.com/SuperElastix/elastix/releases/tag/4.9.0
#


from __future__ import print_function
import sys
import os
import signal
import datetime
import shutil
import subprocess
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

def signal_handler(signal, frame): 
    lprint ('\nUser abort detected\n')
    exit(0)    

def wait_for_keypress():
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
    
def reenable_close_window():
    if pywin32_installed:
        try: # reenable console windows close button (useful if called command line or batch file)
            hwnd = win32console.GetConsoleWindow()
            hMenu = win32gui.GetSystemMenu(hwnd, False)
            win32gui.EnableMenuItem(hMenu, win32con.SC_CLOSE, win32con.MF_ENABLED)        
        except: pass #silent

def exit (code):
    try: logfile.close()
    except: pass
    wait_for_keypress()  
    try: shutil.rmtree(tempdir)
    except: pass # silent
    reenable_close_window()     
    sys.exit(code)        

def lprint (text):
    text = str(text)
    print (text)
    try:
       logfile.write(text+'\n')
       logfile.flush()
    except: pass
            
def permutations():
    dirs  = np.asarray ([[1,2,3],[1,3,2]])
    signs = np.asarray ([[1,1,1],[1,1,-1],[1,-1,1],[-1,1,1],[-1,-1,1],[-1,1,-1],[1,-1,-1],[-1,-1,-1]])
    result = np.zeros ((dirs.shape[0]*signs.shape[0],3), dtype=int)
    for i in range(dirs.shape[0]):
        for j in range(signs.shape[0]):
              result [i*signs.shape[0]+j,:]=dirs[i,:]*signs[j,:]
    return result
    
def invert_transpose(transp):
    #do some test to guarantee correct results
    try:test = np.asarray(transp)
    except: lprint('Error inverting transpose array (1)'); exit(0)
    if not np.array_equal(np.sort(test),range(test.shape[0])): 
        lprint('Error inverting transpose array (2)'); exit(0)
    #actually invert the array
    result=[]
    for i in range(len(transp)): result.append(transp.index(i))
    return result

def deletefile(filename):
    try: os.remove(filename)
    except: pass
    
def copyfile(filename1,filename2):
    try: shutil.copy2 (filename1, filename2)
    except: pass  
    
def run (command):
    process = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)   
    (stdout, stderr) = process.communicate()
    if process.returncode!=0:
       lprint ('WARNING: Subprocess call returned with error (see logfile for details)')
       logfile.write (command)   
       logfile.write (stderr)
       logfile.write (stdout)
       
def make_rotx_parameters():
    f = open(os.path.join(tempdir,'rotx_parameters.txt'), "w")
    f.write ('(FixedInternalImagePixelType "float")\n')
    f.write ('(FixedImageDimension 3)\n')
    f.write ('(MovingInternalImagePixelType "float")\n')
    f.write ('(MovingImageDimension 3)\n')
    f.write ('(ResultImagePixelType "float")\n')
    f.write ('(Registration "MultiResolutionRegistration")\n')
    f.write ('(FixedImagePyramid "FixedSmoothingImagePyramid")\n')
    f.write ('(MovingImagePyramid "MovingSmoothingImagePyramid")\n')
    f.write ('(Interpolator "BSplineInterpolator")\n')
    f.write ('(Metric "AdvancedMattesMutualInformation")\n')
    f.write ('(Optimizer "FullSearch")\n')
    f.write ('(ResampleInterpolator "FinalBSplineInterpolator")\n')
    f.write ('(Resampler "DefaultResampler")\n')
    f.write ('(Transform "EulerTransform")\n')
    f.write ('(NumberOfResolutions 1)\n')
    f.write ('(ImagePyramidSchedule 0)\n')
    f.write ('(AutomaticScalesEstimation "true")\n')
    # if necessary uncomment following line and change preceeding to "false"
    # f.write ('(Scales 100000 100000 100000 100 100 100)\n')
    f.write ('(AutomaticTransformInitialization "true")\n')
    f.write ('(AutomaticTransformInitializationMethod "GeometricalCenter")\n')
    f.write ('(HowToCombineTransforms "Compose")\n')
    f.write ('(UseDirectionCosines "false")\n')
    f.write ('(FullSearchSpace0  "rotation_x" 0 -0.8 0.8 0.2)\n')
    # if necessary uncomment following line(s) and comment preceeding    
    # f.write ('(FullSearchSpace0 "rotation_x" 0 -0.8 0.8 0.2 \
    #           "translation_y" 4 -2 2 0.5 "translation_z" 5 -2 2 0.5)\n')       
    f.write ('(NumberOfHistogramBins 64)\n')    
    f.write ('(FixedLimitRangeRatio 0.0)\n')
    f.write ('(MovingLimitRangeRatio 0.0)\n')
    f.write ('(FixedKernelBSplineOrder 1)\n')
    f.write ('(MovingKernelBSplineOrder 3)\n')
    f.write ('(WriteTransformParametersEachIteration "false")\n')
    f.write ('(WriteTransformParametersEachResolution "false")\n')
    f.write ('(WriteResultImage "false")\n')
    f.write ('(ResultImageFormat "nii.gz")\n')
    f.write ('(CompressResultImage "true")\n')
    f.write ('(ShowExactMetricValue "false")\n')
    f.write ('(ErodeMask "true")\n')
    f.write ('(UseDifferentiableOverlap "false")\n')
    f.write ('(ImageSampler "RandomCoordinate")\n')
    f.write ('(NewSamplesEveryIteration "false")\n')
    f.write ('(BSplineInterpolationOrder 1)\n')
    f.write ('(FinalBSplineInterpolationOrder 3)\n')
    f.write ('(DefaultPixelValue 0)\n')
    f.write ('(CheckNumberOfSamples "true")\n')
    f.write ('(UseMultiThreadingForMetrics "true")\n')
    f.write ('(NumberOfSpatialSamples 5000)\n')
    f.write ('(FixedImageBSplineInterpolationOrder 1)\n')
    f.write ('(UseRandomSampleRegion "true")\n')
    f.write ('(NumberOfFixedHistogramBins 64)\n')
    f.write ('(NumberOfMovingHistogramBins 64)\n')    
    f.write ('(UseFastAndLowMemoryVersion "true")\n')
    f.write ('(UseJacobianPreconditioning "false")\n')
    f.write ('(FiniteDifferenceDerivative "false")\n')
    f.write ('\n')    
    f.close() 
  
def make_rigid_parameters():
    f = open(os.path.join(tempdir,'rigid_parameters.txt'), "w")    
    f.write ('(NumberOfSpatialSamples 1000 2000 4000)\n')
    f.write ('(MaximumNumberOfIterations 100 200 400)\n')      
    f.write ('(UseRandomSampleRegion "true")\n') 
    f.write ('(AutomaticScalesEstimation "true")\n') 
    f.write ('(AutomaticParameterEstimation "true")\n')
    f.write ('(MaximumStepLength 4 2 1)\n')
    f.write ('(MovingImageDimension 3)\n')
    f.write ('(HowToCombineTransforms "Compose")\n')
    f.write ('(Transform "EulerTransform")\n')
    f.write ('(FixedImageDimension 3)\n')
    f.write ('(NewSamplesEveryIteration "true")\n')
    f.write ('(Interpolator "BSplineInterpolator")\n')
    f.write ('(FinalBSplineInterpolationOrder 3)\n')
    f.write ('(WriteResultImage "false")\n')
    f.write ('(ResultImageFormat "nii.gz")\n')
    f.write ('(CompressResultImage "true")\n')
    f.write ('(Resampler "DefaultResampler")\n')
    f.write ('(NumberOfHistogramBins 32 64 64)\n')
    f.write ('(Optimizer "AdaptiveStochasticGradientDescent")\n')
    f.write ('(MovingImagePyramid "MovingRecursiveImagePyramid")\n')
    f.write ('(FixedImagePyramid "FixedRecursiveImagePyramid")\n') 
    f.write ('(Metric "AdvancedMattesMutualInformation")\n')
    f.write ('(AutomaticTransformInitialization "true")\n')
    f.write ('(ImageSampler "RandomCoordinate")\n')
    f.write ('(DefaultPixelValue "0")\n')
    f.write ('(UseDirectionCosines "false")\n')
    f.write ('(NumberOfResolutions "3")\n')
    f.write ('(ImagePyramidSchedule 3 2 1)\n') 
    f.write ('(Registration "MultiResolutionRegistration")\n')
    f.write ('(ResultImagePixelType "float")\n')
    f.write ('(MovingInternalImagePixelType "float")\n')
    f.write ('(FixedInternalImagePixelType "float")\n')
    f.write ('(BSplineInterpolationOrder "1")\n')
    f.write ('(ResampleInterpolator "FinalBSplineInterpolator")\n')
    f.write ('\n')    
    f.close()
    
def make_bspline_parameters():
    f = open(os.path.join(tempdir,'bspline_parameters.txt'), "w")    
    f.write ('(NumberOfSpatialSamples 1000 2000 4000)\n')
    f.write ('(MaximumNumberOfIterations 100 200 400)\n')
    f.write ('(UseRandomSampleRegion "true")\n') 
    f.write ('(AutomaticScalesEstimation "true")\n') 
    f.write ('(AutomaticParameterEstimation "true")\n')    
    f.write ('(MovingImageDimension 3)\n')
    f.write ('(Transform "BSplineTransform")\n')
    f.write ('(FixedImageDimension 3)\n')
    f.write ('(NewSamplesEveryIteration "true")\n')
    f.write ('(FixedImagePyramid "FixedRecursiveImagePyramid")\n')
    f.write ('(FinalBSplineInterpolationOrder 3)\n')
    f.write ('(WriteResultImage "false")\n')
    f.write ('(ResultImageFormat "nii.gz")\n')
    f.write ('(CompressResultImage "true")\n')
    f.write ('(Resampler "DefaultResampler")\n')
    f.write ('(FinalGridSpacingInVoxels 50 50 50)\n')
    f.write ('(NumberOfHistogramBins 32 64 64)\n')
    f.write ('(Optimizer "AdaptiveStochasticGradientDescent")\n')
    f.write ('(MovingImagePyramid "MovingRecursiveImagePyramid")\n')
    f.write ('(FixedInternalImagePixelType "float")\n')
    f.write ('(BSplineInterpolationOrder "1")\n')
    f.write ('(HowToCombineTransforms "Compose")\n')
    f.write ('(ImageSampler "RandomCoordinate")\n')
    f.write ('(DefaultPixelValue "0")\n')
    f.write ('(UseDirectionCosines "false")\n')
    f.write ('(NumberOfResolutions 3)\n')
    f.write ('(ImagePyramidSchedule 3 2 1)\n')    
    f.write ('(Metric "AdvancedMattesMutualInformation")\n')
    f.write ('(Registration "MultiResolutionRegistration")\n')
    f.write ('(ResultImagePixelType "float")\n')
    f.write ('(MovingInternalImagePixelType "float")\n')
    f.write ('(Interpolator "BSplineInterpolator")\n')
    f.write ('(ResampleInterpolator "FinalBSplineInterpolator")\n')
    f.write ('\n')    
    f.close()
      
def get_metric():
    result=0   
    with open(os.path.join(tempdir,'elastix.log')) as f:
        lines = f.readlines()
    for l in range (len(lines)):
          if lines[l].startswith('Final metric value  = '):
             try: result = float(lines[l].split('=')[1])
             except: pass
    return result

          
#general initialization stuff
space=' '; slash='/'; 
if sys.platform=="win32": slash='\\' # not really needed, but looks nicer ;)
Program_name = os.path.basename(sys.argv[0]); 
if Program_name.find('.')>0: Program_name = Program_name[:Program_name.find('.')]
try: resourcedir = sys._MEIPASS+slash # when on PyInstaller 
except: resourcedir = os.path.abspath(os.path.dirname(sys.argv[0]))+slash; 
python_version=str(sys.version_info[0])+'.'+str(sys.version_info[1])+'.'+str(sys.version_info[2])
# sys.platform = [linux2, win32, cygwin, darwin, os2, os2emx, riscos, atheos, freebsd7, freebsd8]
if sys.platform=="win32": 
    os.system("title "+Program_name)
    if pywin32_installed:
        try: # disable console windows close button (substitutes catch shell exit under linux)
            hwnd = win32console.GetConsoleWindow()
            hMenu = win32gui.GetSystemMenu(hwnd, False)
            win32gui.EnableMenuItem(hMenu, win32con.SC_CLOSE, win32con.MF_GRAYED)          
        except: pass #silent 
signal.signal(signal.SIGINT, signal_handler)  # keyboard interrupt
signal.signal(signal.SIGTERM, signal_handler) # kill/shutdown
if  'SIGHUP' in dir(signal): signal.signal(signal.SIGHUP, signal_handler)  # shell exit (linux)
#check for external executables
if not os.path.isfile(os.path.join(resourcedir,'elastix.exe')):
    print ('ERROR:  Elastix executable not found '); exit(1)
if not os.path.isfile(os.path.join(resourcedir,'transformix.exe')):
    print ('ERROR:  Transformix executable not found '); exit(1)    
if not os.path.isfile(os.path.join(resourcedir,'ANNlib-4.9.dll')):
    print ('ERROR:  Elastix DLL not found '); exit(1)    
   
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

#interactively choose input NIFTI files
FIDfile1 = askopenfilename(title="Choose moving NIFTI file", filetypes=[("NIFTI files",('*.nii','*.NII','*.nii.gz','*.NII.GZ'))])
if FIDfile1 == "": print ('ERROR: 1st input file not specified'); exit(2)
FIDfile1 = os.path.abspath(FIDfile1) 
FIDfile2 = askopenfilename(title="Choose fixed reference NIFTI file", filetypes=[("NIFTI files",('*.nii','*.NII','*.nii.gz','*.NII.GZ'))])
if FIDfile2 == "": print ('ERROR: 2nd input file not specified'); exit(2)
FIDfile2 = os.path.abspath(FIDfile2) 
TKwindows.update()
try: win32gui.SetForegroundWindow(win32console.GetConsoleWindow())
except: pass #silent
FIDfile1=str(FIDfile1)
FIDfile2=str(FIDfile2)

# Set Output filenames
dirname  = os.path.dirname(FIDfile1)
basename = os.path.basename(FIDfile1)
basename = basename[0:basename.rfind('.nii.gz')]
outfile = basename+'_Registered.nii.gz'
logname = basename+'_Register.log'
try: logfile = open(os.path.join(dirname,logname), "w")
except: print ('ERROR opening logfile'); exit(2)
logfile.write ('Registration of\n')
logfile.write(FIDfile1+'\nto\n'+FIDfile2+'\n\n');logfile.flush()

# make tempdir
timestamp=datetime.datetime.now().strftime("%Y%m%d%H%M%S")
tempname='.temp'+timestamp
tempdir=os.path.join(dirname,tempname)
if os.path.isdir(tempdir): # this should never happen
    lprint ('ERROR:  Problem creating temp dir (already exists)'); exit(1) 
try: os.mkdir (tempdir)
except: lprint ('ERROR:  Problem creating temp dir: '+tempdir); exit(1)

lprint ('Read  NIFTI files')
#read moving NIFTI file
img_moving = nib.load(FIDfile1)
data_moving = img_moving.get_data().astype(np.float32)
SpatResol_moving = np.asarray(img_moving.header.get_zooms())
Shape_moving = np.asarray(img_moving.header.get_data_shape())
del img_moving # free memory
#find and fix main directions
FOV_moving = data_moving.shape*SpatResol_moving
directions_moving = np.argsort(FOV_moving)
directions_moving = directions_moving[::-1] # decreasing order
if directions_moving[0]==1: transpose_moving = [1,2,0]
elif directions_moving[0]==2: transpose_moving = [2,0,1]
else: transpose_moving = [0,1,2]
lprint ('Moving image transposition: '+np.array2string(np.asarray(transpose_moving)+1))
data_moving = np.transpose (data_moving, axes=transpose_moving)
SpatResol_moving = SpatResol_moving[transpose_moving]
Shape_moving = Shape_moving[transpose_moving]

permutation_arr = permutations()
for i in range (permutation_arr.shape[0]):
    data_moving_perm = np.transpose (data_moving, axes=np.abs(permutation_arr[i])-1).copy()
    if permutation_arr[i,0]<0: data_moving_perm[:,:,:]=data_moving_perm[::-1,:,:]
    if permutation_arr[i,1]<0: data_moving_perm[:,:,:]=data_moving_perm[:,::-1,:]
    if permutation_arr[i,2]<0: data_moving_perm[:,:,:]=data_moving_perm[:,:,::-1] 
    SpatResol_moving_perm = SpatResol_moving[np.abs(permutation_arr[i])-1]
    Shape_moving_perm = Shape_moving[np.abs(permutation_arr[i])-1]
    #transform to int
    max_moving_perm = np.amax(data_moving_perm)
    data_moving_perm *= 32767./max_moving_perm
    data_moving_perm = data_moving_perm.astype(np.int16)    
    #write NIFTI files
    aff = np.eye(4)
    aff[0,0] = -SpatResol_moving_perm[0]; aff[0,3] = -(Shape_moving_perm[0]/2)*aff[0,0]
    aff[1,1] = -SpatResol_moving_perm[1]; aff[1,3] = -(Shape_moving_perm[1]/2)*aff[1,1]
    aff[2,2] = SpatResol_moving_perm[2]; aff[2,3] = -(Shape_moving_perm[2]/2)*aff[2,2]
    img_moving_perm = nib.Nifti1Image(data_moving_perm, aff)
    img_moving_perm.header.set_slope_inter(max_moving_perm/32767.,0)
    img_moving_perm.header.set_xyzt_units(3, 8)   
    img_moving_perm.set_sform(aff, code=1)
    img_moving_perm.set_qform(aff, code=1)
    infile_moving_perm=os.path.join(tempdir,'moving'+str(i)+'.nii')
    nib.save(img_moving_perm, infile_moving_perm)
#free memory
del img_moving_perm; del data_moving_perm; del data_moving

#read fixed NIFTI file
img_fixed = nib.load(FIDfile2)
hdr_fixed = img_fixed.header; affine_fixed = img_fixed.affine #save header info
data_fixed = img_fixed.get_data().astype(np.float32)
SpatResol_fixed = np.asarray(img_fixed.header.get_zooms())
Shape_fixed = np.asarray(img_fixed.header.get_data_shape())
del img_fixed # free memory
#find and fix main directions
FOV_fixed = data_fixed.shape*SpatResol_fixed
directions_fixed = np.argsort(FOV_fixed)
directions_fixed = directions_fixed[::-1] # decreasing order
if directions_fixed[0]==1: transpose_fixed = [1,2,0]
elif directions_fixed[0]==2: transpose_fixed = [2,0,1]
else: transpose_fixed = [0,1,2]
lprint ('Fixed  image transposition: '+np.array2string(np.asarray(transpose_fixed)+1))
data_fixed  = np.transpose (data_fixed, axes=transpose_fixed)
SpatResol_fixed = SpatResol_fixed[transpose_fixed]
Shape_fixed = Shape_fixed[transpose_fixed]
#transform to int
max_fixed = np.amax(data_fixed)
data_fixed *= 32767./max_fixed
data_fixed = data_fixed.astype(np.int16)
#write NIFTI file
aff = np.eye(4)
aff[0,0] = -SpatResol_fixed[0]; aff[0,3] = -(Shape_fixed[0]/2)*aff[0,0]
aff[1,1] = -SpatResol_fixed[1]; aff[1,3] = -(Shape_fixed[1]/2)*aff[1,1]
aff[2,2] = SpatResol_fixed[2]; aff[2,3] = -(Shape_fixed[2]/2)*aff[2,2]
img_fixed = nib.Nifti1Image(data_fixed, aff)
img_fixed.header.set_slope_inter(max_fixed/32767.,0)
img_fixed.header.set_xyzt_units(3, 8)
img_fixed.set_sform(aff, code=1)
img_fixed.set_qform(aff, code=1)
infile_fixed=os.path.join(tempdir,'fixed.nii')
nib.save(img_fixed, infile_fixed)
#free memory
del img_fixed; del data_fixed

#get started    
make_rigid_parameters()                        
permutation_arr = permutations()
metric = np.zeros(shape=permutation_arr.shape[0], dtype=np.float32)
make_rotx_parameters()
lprint ('Start registration')
for i in range (permutation_arr.shape[0]):
    infile_moving_perm=os.path.join(tempdir,'moving'+str(i)+'.nii')
    copyfile (os.path.join(tempdir,'affine_'+str(i+1)+'_transform.txt'),
              os.path.join(tempdir,'affine_transform.txt'))     
    command = '"'+os.path.join(resourcedir,'elastix.exe')+'" '
    command +='-m "'+infile_moving_perm+'" '
    command +='-f "'+infile_fixed+'" '
    command +='-out "'+tempdir+'" '
    command +='-p  "'+os.path.join(tempdir,'rotx_parameters.txt')+'" '      
    run (command)
    copyfile (os.path.join(tempdir,'TransformParameters.0.txt'),
              os.path.join(tempdir,'rotx_transform.txt'))    
    command = '"'+os.path.join(resourcedir,'elastix.exe')+'" '
    command +='-m "'+infile_moving_perm+'" '
    command +='-f "'+infile_fixed+'" '
    command +='-out "'+tempdir+'" '
    command +='-p  "'+os.path.join(tempdir,'rigid_parameters.txt')+'" '
    command +='-t0 "'+os.path.join(tempdir,'rotx_transform.txt')+'" '  
    run (command)
    metric[i]=get_metric()
    lprint ('Permutation  %10s %2d: %f' % (np.array2string(permutation_arr[i]), i, metric[i]))
    copyfile (os.path.join(tempdir,'rotx_transform.txt'),
              os.path.join(tempdir,'rotx_'+str(i+1)+'_transform.txt'))
    copyfile (os.path.join(tempdir,'TransformParameters.0.txt'),
              os.path.join(tempdir,'rigid_'+str(i+1)+'_transform.txt'))   
#save best for later use                
optimum = np.argmin(metric)
copyfile (os.path.join(tempdir,'rigid_'+str(optimum+1)+'_transform.txt'),
          os.path.join(tempdir,'rigid_transform.txt'))
copyfile (os.path.join(tempdir,'rotx_'+str(optimum+1)+'_transform.txt'),
          os.path.join(tempdir,'rotx_transform.txt'))                         
#cleanup    
for i in range (permutation_arr.shape[0]):   
    deletefile(os.path.join(tempdir,'rigid_' +str(i+1)+'_transform.txt'))
    deletefile(os.path.join(tempdir,'rotx_' +str(i+1)+'_transform.txt'))        
lprint ('Best  permutation is at %2d: %f' % (optimum, metric[optimum]))    
    
#bspline registration
lprint ('Start Bspline registration')
infile_moving_perm=os.path.join(tempdir,'moving'+str(optimum)+'.nii') 
make_bspline_parameters()
command = '"'+os.path.join(resourcedir,'elastix.exe')+'" '
command +='-m "'+infile_moving_perm+'" '
command +='-f "'+infile_fixed+'" '
command +='-out "'+tempdir+'" '
command +='-p  "'+os.path.join(tempdir,'bspline_parameters.txt')+'" '
command +='-t0 "'+os.path.join(tempdir,'rigid_transform.txt')+'" '
run (command)
copyfile (os.path.join(tempdir,'TransformParameters.0.txt'),
          os.path.join(tempdir,'bspline_transform.txt'))  
#apply transforms
lprint ('Apply transforms') 
command = '"'+os.path.join(resourcedir,'transformix.exe')+'" '
command +='-in "'+infile_moving_perm+'" '
command +='-out "'+tempdir+'" '
command +='-tp "'+os.path.join(tempdir,'bspline_transform.txt')+'" '
run (command)
#rewrite NIFTI
img = nib.load(os.path.join(tempdir,'result.nii.gz'))
data = img.get_data().astype(np.float32)
SpatResol = np.asarray(img.header.get_zooms())
Shape = np.asarray(img.header.get_data_shape())
#fix directions
transpose_fixed_inv = invert_transpose(transpose_fixed)
lprint ('Inverse fixed image transposition: '+np.array2string(np.asarray(transpose_fixed_inv)+1))
data = np.transpose (data, axes=transpose_fixed_inv)
SpatResol = SpatResol[transpose_fixed_inv]
Shape = Shape[transpose_fixed_inv]
#transform to int
max_ = np.amax(data)
data *= 32767./max_
data = data.astype(np.int16)
#Write NIFTI
img = nib.Nifti1Image(data, affine_fixed)
img.header.set_slope_inter(max_/32767.,0)
img.header.set_xyzt_units(hdr_fixed.get_xyzt_units()[0],hdr_fixed.get_xyzt_units()[1])
img.set_sform(affine_fixed, int(hdr_fixed['sform_code']))
img.set_sform(affine_fixed, int(hdr_fixed['qform_code']))
try: nib.save(img, os.path.join(dirname,outfile))
except: lprint ('ERROR:  problem while writing results'); exit(1)
lprint ('Successfully written output file '+outfile)

#get transforms
with open(os.path.join(tempdir,'rotx_transform.txt')) as f: content1 = f.readlines()
with open(os.path.join(tempdir,'rigid_transform.txt')) as f: content2 = f.readlines()
with open(os.path.join(tempdir,'bspline_transform.txt')) as f: content3 = f.readlines()
#write transforms
f = open(os.path.join(dirname,basename+'_Register.transform'), "w")
f.write ('//------Transform_file_start------(rotx_transform.txt)\n')        
for item in content1: f.write (item.replace (tempdir, '.'))
f.write ('//------Transform_file_start------(rigid_transform.txt)\n')        
for item in content2: f.write (item.replace (tempdir, '.')) 
f.write ('//------Transform_file_start------(bspline_transform.txt)\n')    
for item in content3: f.write (item.replace (tempdir, '.'))
#write additional stuff   
f.write ('//------Transform_file_start------(additional_parameters.txt)\n')
f.write ('Transpose_Moving=%s,%s,%s\n' % (transpose_moving[0]+1,transpose_moving[1]+1,transpose_moving[2]+1))
f.write ('Transpose_Fixed=%s,%s,%s\n' %  (transpose_fixed[0]+1, transpose_fixed[1]+1, transpose_fixed[2]+1))
f.write ('Transpose_Best=%s,%s,%s\n' %  (permutation_arr[optimum,0], permutation_arr[optimum,1], permutation_arr[optimum,2]))
f.write ('Affine_Matrix=%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s\n' %\
         (affine_fixed[0,0],affine_fixed[0,1],affine_fixed[0,2],affine_fixed[0,3],
          affine_fixed[1,0],affine_fixed[1,1],affine_fixed[1,2],affine_fixed[1,3],
          affine_fixed[2,0],affine_fixed[2,1],affine_fixed[2,2],affine_fixed[2,3],
          affine_fixed[3,0],affine_fixed[3,1],affine_fixed[3,2],affine_fixed[3,3]))
f.write ('Sform_Code=%s\n' % hdr_fixed['sform_code'])
f.write ('Qform_Code=%s\n' % hdr_fixed['qform_code'])
f.write ('xyzt_units1=%s\n' % hdr_fixed.get_xyzt_units()[0])
f.write ('xyzt_units2=%s\n' % hdr_fixed.get_xyzt_units()[1])
f.write ('//------Transform_file_start------(full_header.txt)\n')
f.write (str(hdr_fixed))
f.write ('\n')
f.close()
      
#end
exit(0)

