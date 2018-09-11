#
# registers two NIFTI images: 3D Rigid Body, Mutual Information, RegularStepGradientDescentOptimizer
#
# ----- VERSION HISTORY -----
#
# Version 0.1 - 11, April 2018
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
#    - itk 
#


from __future__ import print_function
try: import win32gui, win32console
except: pass #silent
from math import floor
import sys
import os
import numpy as np
import nibabel as nib
import new # required for ITK work with pyinstaller
import itk
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

def smooth(x,window_len):
    w=np.hanning(window_len)
    s=np.r_[2*x[0]-x[window_len-1::-1],x,2*x[-1]-x[-1:-window_len:-1]]
    w=np.hanning(window_len)
    y=np.convolve(w/w.sum(),s,mode='same')
    return y[window_len:-window_len+1]  

def lprint (text):
    text = str(text)
    print (text);
    logfile.write(text+'\n')
    logfile.flush()

def itk_GetDirectionArray(itkImage):
    arr = np.zeros ((3,3),dtype=np.float32)
    vnl_matrix = itkImage.GetDirection().GetVnlMatrix()
    for i in range(3):
        for j in range(3):
            arr[i,j]=vnl_matrix.get(i,j)
    return arr

def itk_GetiDirectionArray(itkImage):
    arr = np.zeros ((3,3),dtype=np.float32)
    vnl_matrix = itkImage.GetInverseDirection().GetVnlMatrix()
    for i in range(3):
        for j in range(3):
            arr[i,j]=vnl_matrix.get(i,j)
    return arr    

def itk_SetDirectionArray (itkImage,arr):
    for i in range(3):
        for j in range(3):
            itkImage.GetDirection().GetVnlMatrix().set(i,j,arr[i,j])

def itk_SetiDirectionArray (itkImage,arr):
    for i in range(3):
        for j in range(3):
            itkImage.GetInverseDirection().GetVnlMatrix().set(i,j,arr[i,j])            

    
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
nfiles=0
answer="dummy"
FIDfile1 = askopenfilename(title="Choose moving NIFTI file (A)", filetypes=[("NIFTI files",('*.png','*.nii','*.NII','*.nii.gz','*.NII.GZ'))])
if FIDfile1 == "": print ('ERROR: X input file not specified'); sys.exit(2)
FIDfile1 = os.path.abspath(FIDfile1) 
FIDfile2 = askopenfilename(title="Choose reference NIFTI file (B)", filetypes=[("NIFTI files",('*.png','*.nii','*.NII','*.nii.gz','*.NII.GZ'))])
if FIDfile2 == "": print ('ERROR: Y input file not specified'); sys.exit(2)
FIDfile2 = os.path.abspath(FIDfile2) 
TKwindows.update()
try: win32gui.SetForegroundWindow(win32console.GetConsoleWindow())
except: pass #silent
FIDfile1=str(FIDfile1)
FIDfile2=str(FIDfile2)


# Set Output filenames
dirname  = os.path.dirname(FIDfile1)
Outfile = os.path.basename(FIDfile1);
Outfile = Outfile[0:Outfile.rfind('.nii.gz')]+'_RigidBodyRegistered.nii.gz'
logname = os.path.basename(FIDfile1);
logname = logname[0:logname.rfind('.nii.gz')]+'_RigidBodyRegistration.log'
try: logfile = open(os.path.join(dirname,logname), "w")
except: print ('ERROR opening logfile'); sys.exit(2)
logfile.write ('3D Rigid Body registration\nof\n')
logfile.write(FIDfile1+'\nto\n'+FIDfile2+'\n\n');logfile.flush()

lprint ('Reading NIFTI files')
img_moving = nib.load(FIDfile1)
data_moving = img_moving.get_data().astype(np.float32)
SpatResol_moving = np.asarray(img_moving.header.get_zooms())
img_fixed = nib.load(FIDfile2)
data_fixed = img_fixed.get_data().astype(np.float32)
SpatResol_fixed = np.asarray(img_fixed.header.get_zooms())
    
# ------------------       ITK code starts here --------------------

# convert Numpy Array to ITK 
movingImage = itk.GetImageFromArray(data_moving[:,:,:])
itk_SetDirectionArray  (movingImage, np.asarray([[-1,0,0],[0,-1,0],[0,0,1]]))
itk_SetiDirectionArray (movingImage, np.asarray([[-1,0,0],[0,-1,0],[0,0,1]]))
movingImage.SetSpacing(SpatResol_moving.tolist())
Origin_moving = SpatResol_moving*data_moving.shape/2
Origin_moving [2] *= -1 # duno why this is needed
movingImage.SetOrigin(Origin_moving.tolist())
fixedImage  = itk.GetImageFromArray(data_fixed[:,:,:])   
itk_SetDirectionArray  (fixedImage, np.asarray([[-1,0,0],[0,-1,0],[0,0,1]]))
itk_SetiDirectionArray (fixedImage, np.asarray([[-1,0,0],[0,-1,0],[0,0,1]]))
fixedImage.SetSpacing(SpatResol_moving.tolist())
Origin_fixed = SpatResol_fixed*data_fixed.shape/2
Origin_fixed [2] *= -1 # duno why this is needed
fixedImage.SetOrigin(Origin_fixed.tolist())
#lprint (movingImage.GetLargestPossibleRegion().GetSize())
#lprint (movingImage.GetSpacing())
#lprint (movingImage.GetOrigin())
#lprint (itk_GetDirectionArray(movingImage))
#lprint (itk_GetiDirectionArray(movingImage))
#lprint (fixedImage)
#lprint ('')
#lprint (fixedImage.GetLargestPossibleRegion().GetSize())
#lprint (fixedImage.GetSpacing())
#lprint (fixedImage.GetOrigin())
#lprint (itk_GetDirectionArray(fixedImage))
#lprint (itk_GetiDirectionArray(fixedImage))
#lprint (fixedImage)
#lprint ('')

#
#  Define data types
#
FixedImageType   = itk.Image[itk.F, 3]
MovingImageType  = itk.Image[itk.F, 3]
#TransformType    = itk.AffineTransform[itk.D, 3]
TransformType    = itk.VersorRigid3DTransform[itk.D]
OptimizerType    = itk.RegularStepGradientDescentOptimizerv4[itk.D]
RegistrationType = itk.ImageRegistrationMethodv4[FixedImageType,MovingImageType]
#MetricType       = itk.MeanSquaresImageToImageMetricv4[FixedImageType,MovingImageType]
#MetricType       = itk.CorrelationImageToImageMetricv4[FixedImageType,MovingImageType]
MetricType       = itk.MattesMutualInformationImageToImageMetricv4[FixedImageType,MovingImageType]

'''
#  Read the fixed and moving images
lprint ('Reading NIFTI files')
movingImageReader = itk.ImageFileReader[MovingImageType].New()
fixedImageReader  = itk.ImageFileReader[FixedImageType].New()
movingImageReader.SetFileName(FIDfile1)
fixedImageReader.SetFileName(FIDfile2)
movingImageReader.Update()
fixedImageReader.Update()
movingImage = movingImageReader.GetOutput()
fixedImage  = fixedImageReader.GetOutput()
if np.min(fixedImage.GetSpacing())<10.:
   fixedImage.SetSpacing(fixedImage.GetSpacing()*1000) # undo ITK scale to xyz_units=um
if np.min(movingImage.GetSpacing())<10.:
   movingImage.SetSpacing(movingImage.GetSpacing()*1000) # undo ITK scale to xyz_units=um
#lprint (movingImage.GetLargestPossibleRegion().GetSize())
#lprint (movingImage.GetSpacing())
#lprint (movingImage.GetOrigin())
#lprint (itk_GetDirectionArray(movingImage))
#lprint (itk_GetiDirectionArray(movingImage))
#lprint (fixedImage)
#lprint ('')
#lprint (fixedImage.GetLargestPossibleRegion().GetSize())
#lprint (fixedImage.GetSpacing())
#lprint (fixedImage.GetOrigin())
#lprint (itk_GetDirectionArray(fixedImage))
#lprint (itk_GetiDirectionArray(fixedImage))
#lprint (fixedImage)
#lprint ('')
'''

#  Instantiate the classes for the registration framework
registration = RegistrationType.New()
imageMetric  = MetricType.New()
transform    = TransformType.New()
optimizer    = OptimizerType.New()
registration.SetMetric(imageMetric)
registration.SetInitialTransform(transform)
registration.SetOptimizer(optimizer)
registration.SetFixedImage(fixedImage)
registration.SetMovingImage(movingImage)

#  Define optimizer parameters
optimizer.SetLearningRate(1) # 1
optimizer.SetMinimumStepLength(0.001) #0.001
optimizer.SetRelaxationFactor(0.5)
optimizer.SetNumberOfIterations(500)

# One level registration process without shrinking and smoothing.
registration.SetNumberOfLevels(1)
registration.SetSmoothingSigmasPerLevel([0])
registration.SetShrinkFactorsPerLevel([1])

# Iteration Observer
def iterationUpdate():
    #print('.', end='') # just print a progress indicator
    currentParameter = registration.GetOutput().Get().GetParameters()
    lprint ("%f %0.4f %0.4f %0.4f %0.4f %0.4f %0.4f %0.4f %0.4f" % ( optimizer.GetValue(),
                                  currentParameter.GetElement(0),
                                  currentParameter.GetElement(1), 
                                  currentParameter.GetElement(2),
                                  currentParameter.GetElement(3),
                                  currentParameter.GetElement(4),
                                  currentParameter.GetElement(5),
                                  currentParameter.GetElement(6),    
                                  currentParameter.GetElement(7)))

iterationCommand = itk.PyCommand.New()
iterationCommand.SetCommandCallable(iterationUpdate)
optimizer.AddObserver(itk.IterationEvent(),iterationCommand)

#  Start the registration process
lprint ("Starting registration")
registration.Update()


# Get the final parameters of the transformation
finalParameters = registration.GetOutput().Get().GetParameters()
lprint (" ")
lprint ("Final Registration Parameters ")
lprint ("%f %f %f %f %f %f %f %f %f" % (optimizer.GetValue(),
                                  finalParameters.GetElement(0),
                                  finalParameters.GetElement(1), 
                                  finalParameters.GetElement(2),
                                  finalParameters.GetElement(3),
                                  finalParameters.GetElement(4),
                                  finalParameters.GetElement(5),
                                  finalParameters.GetElement(6),
                                  finalParameters.GetElement(7)))

# Now, we use the final transform for resampling the moving image.
resampler = itk.ResampleImageFilter[MovingImageType,FixedImageType].New()
resampler.SetTransform(registration.GetTransform())
resampler.SetInput(movingImage)
region = fixedImage.GetLargestPossibleRegion()
resampler.SetSize(region.GetSize())
resampler.SetOutputOrigin(fixedImage.GetOrigin())
resampler.SetOutputSpacing(fixedImage.GetSpacing())
resampler.SetOutputDirection(fixedImage.GetDirection())
resampler.SetDefaultPixelValue(0)
resampler.Update()
OutputImageType  = itk.Image[itk.F, 3]
outputCast = itk.CastImageFilter[FixedImageType, OutputImageType].New()
outputCast.SetInput(resampler.GetOutput())
outputCast.Update()

# transform ITK images to numpy arrays
img = outputCast.GetOutput()
data = itk.GetArrayViewFromImage(img)
data = np.transpose (data, axes=(2,1,0))
Resolution_out = img.GetSpacing()

#---------------- ITK ends --------------------


#write NIFTI
lprint ('Writing registration result ')
max_data = np.amax(data)
data *= 32767./max_data;
data = data.astype(np.int16)
#createNIFTI's
aff = np.eye(4)
aff[0,0] = Resolution_out[0]; aff[0,3] = -(data.shape[0]/2)*aff[0,0]
aff[1,1] = Resolution_out[1]; aff[1,3] = -(data.shape[1]/2)*aff[1,1]
aff[2,2] = Resolution_out[2]; aff[2,3] = -(data.shape[2]/2)*aff[2,2]
NIFTIimg = nib.Nifti1Image(data[:,:,:], aff)
NIFTIimg.header.set_slope_inter(max_data/32767.,0)
NIFTIimg.header.set_xyzt_units(3, 8)
NIFTIimg.set_sform(aff, code=0)
NIFTIimg.set_qform(aff, code=1)
try: nib.save(NIFTIimg, os.path.join(dirname,Outfile))
except: lprint ('\nERROR:  problem while writing results'); sys.exit(1)
lprint ("\ndone\n")
logfile.close()
 
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