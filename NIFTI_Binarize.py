#
# This is basically a wrapper arround some ANTs executables
#
# the objective of this tool is to binarize a grayscale NIFTI image
# output is a written to a NIFTI file that containes only 0 and 1 values.
# the method used is segmentation with the kmeans algorithm
# using the "Atropos" tool from the ANTs library (http://stnava.github.io/ANTs/)
#
#
# ----- VERSION HISTORY -----
#
# Version 0.1 - 01, November 2018
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
#    - ANTs executables:
#       1) N4BiasFieldCorrection.exe
#       2) Atropos.exe
#      from:
#      https://github.com/ANTsX/ANTs/releases/tag/v2.1.0
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
if not os.path.isfile(os.path.join(resourcedir,'Atropos.exe')):
    print ('ERROR:  Atropos executable not found '); exit(1)
if not os.path.isfile(os.path.join(resourcedir,'N4BiasFieldCorrection.exe')):
    print ('ERROR:  N4BiasFieldCorrection executable not found '); exit(1)       
   
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

#interactively choose input NIFTI file
FIDfile = askopenfilename(title="Choose NIFTI file", filetypes=[("NIFTI files",('*.nii','*.NII','*.nii.gz','*.NII.GZ'))])
if FIDfile == "": print ('ERROR: No input file specified'); exit(2)
FIDfile = os.path.abspath(FIDfile) 
FIDfile=str(FIDfile)
TKwindows.update()
try: win32gui.SetForegroundWindow(win32console.GetConsoleWindow())
except: pass #silent


# Set Output filenames
dirname  = os.path.dirname(FIDfile)
basename = os.path.basename(FIDfile)
basename = basename[0:basename.rfind('.nii.gz')]
outfile = basename+'_Binarized.nii.gz'
logname = basename+'_Binarize.log'
try: logfile = open(os.path.join(dirname,logname), "w")
except: print ('ERROR opening logfile'); exit(2)
logfile.write ('Binarization of\n')
logfile.write(FIDfile+'\n\n');logfile.flush()

# make tempdir
timestamp=datetime.datetime.now().strftime("%Y%m%d%H%M%S")
tempname='.temp'+timestamp
tempdir=os.path.join(dirname,tempname)
if os.path.isdir(tempdir): # this should never happen
    lprint ('ERROR:  Problem creating temp dir (already exists)'); exit(1) 
try: os.mkdir (tempdir)
except: lprint ('ERROR:  Problem creating temp dir: '+tempdir); exit(1)

#temp filenames
infile0 = os.path.join(tempdir,'input.nii.gz')
infile1 = os.path.join(tempdir,'input_biascorrected.nii.gz')
maskfile = os.path.join(tempdir,'mask.nii.gz')

#copy NIFTI file to tempdir
copyfile (FIDfile, infile0)

#create mask file (just a file with ones)
img = nib.load(infile0)
data = img.get_data().astype(np.float32)
SpatResol = np.asarray(img.header.get_zooms())
Shape = np.asarray(img.header.get_data_shape())
affine = img.affine
xyzt_units1  = img.header.get_xyzt_units()[0]
xyzt_units2  = img.header.get_xyzt_units()[1]
sform = int(img.header['sform_code'])
qform = int(img.header['qform_code'])
data.fill(1); data = data.astype(np.int16)
img = nib.Nifti1Image(data, affine)
img.header.set_slope_inter(1,0)
img.header.set_xyzt_units(xyzt_units1,xyzt_units2)
img.set_sform(affine, sform)
img.set_qform(affine, qform)
try: nib.save(img, maskfile)
except: print ('ERROR:  Problem writing temporary mask file'); exit(1)
del img; del data # free memory

#bias field correction  
lprint ('Start bias field correction')     
command = '"'+os.path.join(resourcedir,'N4BiasFieldCorrection.exe')+'" '
command +='-d 3 '
command +='-i "'+infile0+'" '
command +='-o "'+infile1+'" '
run (command)

#segmentation  
lprint ('Start binarization')     
command = '"'+os.path.join(resourcedir,'Atropos.exe')+'" '
command +='-d 3 '
command +='-a "'+infile1+'" '
command +='-i kmeans[4] -p Socrates[1] '
command +='-x "'+maskfile+'" '
command +='-m [0.1,1x1x1] '
command +='-o ["'
command +=os.path.join(tempdir,'Segmentation.nii.gz')+'","'
command +=os.path.join(tempdir,'SegmentationPosteriors%d.nii.gz')+'"] '
run (command)

#read segmentation result
img = nib.load(os.path.join(tempdir,'SegmentationPosteriors4.nii.gz'))
data = img.get_data().astype(np.float32) 

#threshold binarize
data /= np.amax(data)
threshold = 0.5
data[data>threshold]=1
data[data<=threshold]=0
#data = data.astype(np.uint8)
data = data.astype(np.int16)

#rewrite NIFTI
img = nib.Nifti1Image(data, affine)
img.header.set_slope_inter(1,0)
img.header.set_xyzt_units(xyzt_units1,xyzt_units2)
img.set_sform(affine, sform)
img.set_qform(affine, qform)
try: nib.save(img, outfile)
except: print ('ERROR:  Problem writing output file'); exit(1)
lprint ('Successfully written output file '+outfile)

#end
exit(0)

