#
# reads 3 NIFTI files that contain X,Y,Z components
# combines the input to VTK vector format and regularizes
# the vector field with the Gradient Vector Flow (GVF) algorithm
# as implemented in the Insight Segmentation and Registration Toolkit (ITK)
#    https://itk.org/  
# which is based on: 
#    Xu et al., "Snakes, Shapes, and Gradient Vector Flow", 
#    IEEE Transactions on Image Processing, vol. 7, No. 3, Mar. 1998, pp. 359-369
#    DOI: 10.1109/83.661186 
#    https://ieeexplore.ieee.org/stamp/stamp.jsp?tp=&arnumber=661186
# the original and regularized vector images are saved in MHA format
# the individual X,Y,Z components resulting from the regularization 
# are saved in NIFTI format. ForConvenience, the magnitude of the regularized 
# components is saved too
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
#


from __future__ import print_function
try: import win32gui, win32console
except: pass #silent
from math import floor
import sys
import os
import numpy as np
import nibabel as nib
import itk
# to test pyinstaller
'''
ImportProblem = False
try: from itk import ImageFileReader
except: print ('Error importing ITK module ImageFileReader'); ImportProblem = True
try: from itk import ImageFileWriter
except: print ('Error importing ITK module ImageFileWriter'); ImportProblem = True
try: from itk import ComposeImageFilter
except: print ('Error importing ITK module ComposeImageFilter'); ImportProblem = True
try: from itk import GradientVectorFlowImageFilter
except: print ('Error importing ITK module GradientVectorFlowImageFilter'); ImportProblem = True
try: from itk import VectorIndexSelectionCastImageFilter
except: print ('Error importing ITK module VectorIndexSelectionCastImageFilter'); ImportProblem = True
try: from itk import GetArrayViewFromImage
except: print ('Error importing ITK module GetArrayViewFromImage'); ImportProblem = True
if ImportProblem: sys.exit(2)
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

#intercatively choose input NIFTI files
nfiles=0
answer="dummy"
FIDfile1 = askopenfilename(title="Choose NIFTI file X component", filetypes=[("NIFTI files",('*.nii','*.NII','*.nii.gz','*.NII.GZ'))])
if FIDfile1 == "": print ('ERROR: X input file not specified'); sys.exit(2)
FIDfile1 = os.path.abspath(FIDfile1) 
FIDfile2 = askopenfilename(title="Choose NIFTI file Y component", filetypes=[("NIFTI files",('*.nii','*.NII','*.nii.gz','*.NII.GZ'))])
if FIDfile2 == "": print ('ERROR: Y input file not specified'); sys.exit(2)
FIDfile2 = os.path.abspath(FIDfile2) 
FIDfile3 = askopenfilename(title="Choose NIFTI file Z component", filetypes=[("NIFTI files",('*.nii','*.NII','*.nii.gz','*.NII.GZ'))])
if FIDfile3 == "": print ('ERROR: Z input file not specified'); sys.exit(2)
FIDfile3 = os.path.abspath(FIDfile3) 
TKwindows.update()
try: win32gui.SetForegroundWindow(win32console.GetConsoleWindow())
except: pass #silent
FIDfile1=str(FIDfile1)
FIDfile2=str(FIDfile2)
FIDfile3=str(FIDfile3)

# Set Output filenames
logname = os.path.basename(FIDfile1);
logname = logname[0:logname.rfind('_')];logname = logname[0:logname.rfind('_')]+'-Regularization.log'
Vectorfile = os.path.basename(FIDfile1);
Vectorfile = Vectorfile[0:Vectorfile.rfind('_')]+'Vector.mha'
rVectorfile = os.path.basename(FIDfile1);
rVectorfile = rVectorfile[0:rVectorfile.rfind('_')]+'Vector-regularized.mha'
OutXfile = os.path.basename(FIDfile1);
OutXfile = OutXfile[0:OutXfile.rfind('.nii.gz')]+'-regularized.nii.gz'
OutYfile = os.path.basename(FIDfile2);
OutYfile = OutYfile[0:OutYfile.rfind('.nii.gz')]+'-regularized.nii.gz'
OutZfile = os.path.basename(FIDfile3);
OutZfile = OutZfile[0:OutZfile.rfind('.nii.gz')]+'-regularized.nii.gz'
Out_file = os.path.basename(FIDfile1); # SumOfSquares
Out_file = Out_file[0:Out_file.rfind('_')]+'-regularized.nii.gz'      

#make results folder
dirname = os.path.abspath(os.path.dirname(FIDfile1)+slash+'..'+slash+'VectorRegularizationResult')
new_dirname = dirname
i=0
while os.path.exists(new_dirname):
   i+=1
   new_dirname = dirname+'('+str(i)+')'
try: os.makedirs(new_dirname)
except: print ('ERROR: unable to make folder', new_dirname); sys.exit(2)



#write logfile     
with open(os.path.join(new_dirname,logname), "w") as logfile:
    logfile.write('Result VTK file is combination of:\n')
    logfile.write(FIDfile1+'\n')
    logfile.write(FIDfile2+'\n')
    logfile.write(FIDfile3+'\n')


    
# ------------------       ITK code starts here --------------------

print ('Reading NIFTI Images')
image_X = itk.imread(FIDfile1); image_X.Update()
image_Y = itk.imread(FIDfile2); image_Y.Update()
image_Z = itk.imread(FIDfile3); image_Z.Update()

print ('Writing Vector Image')
composer = itk.ComposeImageFilter[itk.Image.F3, itk.Image.VF33].New()
composer.SetInput(0, image_X)
composer.SetInput(1, image_Y)
composer.SetInput(2, image_Z)
writer = itk.ImageFileWriter[itk.Image.VF33].New(composer.GetOutput())
writer.SetFileName(os.path.join(new_dirname,Vectorfile))
writer.UseCompressionOn ()
writer.Update()

print ('Regularizing Vector Image')
#image_X = itk.MedianImageFilter (image_X, Radius = 1) #only good forimages without much detail, e.g. phantoms
#image_Y = itk.MedianImageFilter (image_Y, Radius = 1) #only good forimages without much detail, e.g. phantoms
#image_Z = itk.MedianImageFilter (image_Z, Radius = 1) #only good forimages without much detail, e.g. phantoms
composer = itk.ComposeImageFilter[itk.Image.F3, itk.Image.VF33].New()
composer.SetInput(0, image_X)
composer.SetInput(1, image_Y)
composer.SetInput(2, image_Z)
# IterationNum MUST be = 1
# txx file says:
#   m_TimeStep = 0.001; m_NoiseLevel = 200; (lines 31-32)
#   m_TimeStep = 0.2/m_NoiseLevel; (line 56)
# from: http://itk-users.7.n7.nabble.com/units-of-SetNoiseLevel-in-itkGradientVectorFlowImageFilter-td20115.html
# higher NoiseLevel smooth more
# lower Timestep idem
# for smoother images increase IterationNum to 2-5
# for phantoms with less degree of details IterationNum = 10 is OK
filtered_image = itk.GradientVectorFlowImageFilter.New(composer.GetOutput(), IterationNum=1, NoiseLevel=2000.0, NumberOfThreads = 2)
filtered_image.Update()

print ('Writing regularized Vector Image')
writer = itk.ImageFileWriter[itk.Image.VF33].New(filtered_image)
writer.SetFileName(os.path.join(new_dirname,rVectorfile))
writer.UseCompressionOn ()
writer.Update()

# transform ITK images to numpy arrays
Xcomponent_itk = itk.VectorIndexSelectionCastImageFilter(filtered_image, Index = 0)
Xcomponent_np = itk.GetArrayViewFromImage(Xcomponent_itk)
Ycomponent_itk = itk.VectorIndexSelectionCastImageFilter(filtered_image, Index = 1)
Ycomponent_np = itk.GetArrayViewFromImage(Ycomponent_itk)
Zcomponent_itk = itk.VectorIndexSelectionCastImageFilter(filtered_image, Index = 2)
Zcomponent_np = itk.GetArrayViewFromImage(Zcomponent_itk)
SpatResol = np.zeros(shape=(3),dtype=np.float32)
SpatResol[0] = Xcomponent_itk.GetSpacing()[0]
SpatResol[1] = Xcomponent_itk.GetSpacing()[1]
SpatResol[2] = Xcomponent_itk.GetSpacing()[2]



# ------------------       ITK code ends here --------------------


# still some flips needed :(
Xcomponent_np = np.transpose (Xcomponent_np, axes=(2,1,0))
Ycomponent_np = np.transpose (Ycomponent_np, axes=(2,1,0))
Zcomponent_np = np.transpose (Zcomponent_np, axes=(2,1,0))
SpatResol_perm = np.zeros(shape=(3),dtype=np.float32)
SpatResol_perm[0]=SpatResol[2]
SpatResol_perm[1]=SpatResol[1]
SpatResol_perm[2]=SpatResol[0]

print ('Writing regularized X,Y,Z Components ', end='')
image_SOS = np.sqrt(np.square(Xcomponent_np) + np.square(Ycomponent_np) + np.square(Zcomponent_np))
max_SOS = np.amax(image_SOS)
image_SOS *= 32767./max_SOS; image_SOS = image_SOS.astype(np.int16)
max_X = np.amax(Xcomponent_np);
max_Y = np.amax(Ycomponent_np);
max_Z = np.amax(Zcomponent_np);
max_ALL = max (max_X,max_Y,max_Z)
Xcomponent_np *= 32767./max_ALL; Xcomponent_np = Xcomponent_np.astype(np.int16)
Ycomponent_np *= 32767./max_ALL; Ycomponent_np = Ycomponent_np.astype(np.int16)
Zcomponent_np *= 32767./max_ALL; Zcomponent_np = Zcomponent_np.astype(np.int16)

#createNIFTI's
aff = np.eye(4)
aff[0,0] = SpatResol_perm[0]; aff[0,3] = -(Xcomponent_np.shape[0]/2)*aff[0,0]
aff[1,1] = SpatResol_perm[1]; aff[1,3] = -(Xcomponent_np.shape[1]/2)*aff[1,1]
aff[2,2] = SpatResol_perm[2]; aff[2,3] = -(Xcomponent_np.shape[2]/2)*aff[2,2]
#write Phase flow X
NIFTIimg = nib.Nifti1Image(Xcomponent_np[:,:,:], aff)
NIFTIimg.header['sform_code']=1
NIFTIimg.header['qform_code']=1
NIFTIimg.header.set_slope_inter(max_ALL/32767.,0)
try: nib.save(NIFTIimg, os.path.join(new_dirname,OutXfile))
except: print ('\nERROR:  problem while writing results'); sys.exit(1)
print('.', end='') #progress indicator
#write Phase flow Y
NIFTIimg = nib.Nifti1Image(Ycomponent_np[:,:,:], aff)
NIFTIimg.header['sform_code']=1
NIFTIimg.header['qform_code']=1
NIFTIimg.header.set_slope_inter(max_ALL/32767.,0)
try: nib.save(NIFTIimg, os.path.join(new_dirname,OutYfile))
except: print ('\nERROR:  problem while writing results'); sys.exit(1)
print('.', end='') #progress indicator
#write Phase flow Z
NIFTIimg = nib.Nifti1Image(Zcomponent_np[:,:,:], aff)
NIFTIimg.header['sform_code']=1
NIFTIimg.header['qform_code']=1
NIFTIimg.header.set_slope_inter(max_ALL/32767.,0)
try: nib.save(NIFTIimg, os.path.join(new_dirname,OutZfile))
except: print ('\nERROR:  problem while writing results'); sys.exit(1)
print('.', end='') #progress indicator
#write Phase flow Z
NIFTIimg = nib.Nifti1Image(image_SOS[:,:,:], aff)
NIFTIimg.header['sform_code']=1
NIFTIimg.header['qform_code']=1
NIFTIimg.header.set_slope_inter(max_SOS/32767.,0)
try: nib.save(NIFTIimg, os.path.join(new_dirname,Out_file))
except: print ('\nERROR:  problem while writing results'); sys.exit(1)
print('.', end='') #progress indicator


print ("\ndone\n") 
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