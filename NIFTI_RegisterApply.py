#
# applies a transform previously calculated with the "NIFTI_Register" tool
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
if not os.path.isfile(os.path.join(resourcedir,'transformix.exe')):
    print ('ERROR:  Transformix executable not found '); exit(1)    
if not os.path.isfile(os.path.join(resourcedir,'ANNlib-4.9.dll')):
    print ('ERROR:  Transformix DLL not found '); exit(1)

   
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

#intercatively choose NIFTI and transform files
nfiles=0
answer="dummy"
FIDfile1 = askopenfilename(title="Choose moving NIFTI file", filetypes=[("NIFTI files",('*.nii','*.NII','*.nii.gz','*.NII.GZ'))])
if FIDfile1 == "": print ('ERROR: NIFTI file not specified'); exit(2)
FIDfile1 = os.path.abspath(FIDfile1) 
TransformFile = askopenfilename(title="Choose transformfile", filetypes=[("Transform files",('*Register.transform'))])
if TransformFile == "": print ('ERROR: Transform file not specified'); exit(2)
TransformFile = os.path.abspath(TransformFile) 
TKwindows.update()
try: win32gui.SetForegroundWindow(win32console.GetConsoleWindow())
except: pass #silent
FIDfile1=str(FIDfile1)
TransformFile=str(TransformFile)

# Set Output filenames
dirname  = os.path.dirname(FIDfile1)
basename = os.path.basename(FIDfile1)
basename = basename[0:basename.rfind('.nii.gz')]
outfile = basename+'_Transformed.nii.gz'
logname = basename+'_Transform.log'
try: logfile = open(os.path.join(dirname,logname), "w")
except: print ('ERROR opening logfile'); exit(2)
logfile.write ('Transform of\n')
logfile.write(FIDfile1+'\nwith\n'+TransformFile+'\n\n');logfile.flush()

# make tempdir
timestamp=datetime.datetime.now().strftime("%Y%m%d%H%M%S")
tempname='.temp'+timestamp
tempdir=os.path.join(dirname,tempname)
if os.path.isdir(tempdir): # this should never happen
    lprint ('ERROR:  Problem creating temp dir (already exists)'); exit(1) 
try: os.mkdir (tempdir)
except: lprint ('ERROR:  Problem creating temp dir: '+tempdir); exit(1)

#read Transpose file
with open(TransformFile) as f: content = f.readlines()
if not content[0].startswith('//------Transform_file_start------('):
    lprint ("Error: Choosen transform file doesn't look correct"); exit(2)
for i in range(len(content)):
    if content[i].startswith('//------Transform_file_start------('):
       s = content[i]
       filename = s[s.rfind('(')+1:-1].strip(')\n')
       try: f.close()
       except: pass
       f = open(os.path.join(tempdir,filename), "w")
    else: # actually write the respective files
       s = content[i]
       if s.startswith('(InitialTransformParametersFileName'):
          s = s.replace('.\\',tempdir+'\\')
       f.write (s)       
f.close()
with open(os.path.join(tempdir,'additional_parameters.txt')) as f: content = f.readlines()
moreparameters_dict = {}
for item in content:
    (param_name, current_line) = item.split('=') #split at "="
    param_name = param_name.strip()
    current_line = current_line.strip()
    value = ParseSingleValue(current_line)
    moreparameters_dict[param_name] = value
t = moreparameters_dict['Transpose_Moving'].split(',')
transpose_moving = [int(i)-1 for i in t]
t = moreparameters_dict['Transpose_Fixed'].split(',')
transpose_fixed = [int(i)-1 for i in t]
t = moreparameters_dict['Transpose_Best'].split(',')
transpose_best = [int(i) for i in t]
t = moreparameters_dict['Affine_Matrix'].split(',')
affine_fixed = [float(i) for i in t]
affine_fixed = np.asarray(affine_fixed).reshape((4,4))
sform = moreparameters_dict['Sform_Code']
qform = moreparameters_dict['Qform_Code']
xyzt_units1 = moreparameters_dict['xyzt_units1']
xyzt_units2 = moreparameters_dict['xyzt_units2']

#read moving NIFTI file
lprint ('Read  NIFTI file')
img_moving = nib.load(FIDfile1)
data_moving = img_moving.get_data().astype(np.float32)
SpatResol_moving = np.asarray(img_moving.header.get_zooms())
Shape_moving = np.asarray(img_moving.header.get_data_shape())
img_moving = 0 # free memory
#fix main directions
lprint ('Moving image transposition: '+np.array2string(np.asarray(transpose_moving)+1))
data_moving = np.transpose (data_moving, axes=transpose_moving)
SpatResol_moving = SpatResol_moving[transpose_moving]
Shape_moving = Shape_moving[transpose_moving]
lprint ('Moving image with best permutation: '+str(transpose_best))
data_moving = np.transpose (data_moving, axes=np.abs(transpose_best)-1)
if transpose_best[0]<0: data_moving[:,:,:]=data_moving[::-1,:,:]
if transpose_best[1]<0: data_moving[:,:,:]=data_moving[:,::-1,:]
if transpose_best[2]<0: data_moving[:,:,:]=data_moving[:,:,::-1] 
SpatResol_moving = SpatResol_moving[np.abs(transpose_best)-1]
Shape_moving = Shape_moving[np.abs(transpose_best)-1]
#transform to int
max_moving = np.amax(data_moving)
data_moving *= 32767./max_moving
data_moving = data_moving.astype(np.int16)  
#write NIFTI file
aff = np.eye(4)
aff[0,0] = -SpatResol_moving[0]; aff[0,3] = -(Shape_moving[0]/2)*aff[0,0]
aff[1,1] = -SpatResol_moving[1]; aff[1,3] = -(Shape_moving[1]/2)*aff[1,1]
aff[2,2] =  SpatResol_moving[2]; aff[2,3] = -(Shape_moving[2]/2)*aff[2,2]
img_moving = nib.Nifti1Image(data_moving, aff)
img_moving.header.set_slope_inter(max_moving/32767.,0)
img_moving.header.set_xyzt_units(3, 8)
img_moving.set_sform(aff, code=1)
img_moving.set_qform(aff, code=1)
nib.save(img_moving, os.path.join(tempdir,'moving.nii'))
#free memory
img_moving = 0; data_moving = 0
 
#apply transforms
lprint ('Apply transforms') 
command = '"'+os.path.join(resourcedir,'transformix.exe')+'" '
command +='-in "'+os.path.join(tempdir,'moving.nii')+'" '
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
img.header.set_xyzt_units(xyzt_units1,xyzt_units2)
img.set_sform(affine_fixed, sform)
img.set_sform(affine_fixed, qform)
try: nib.save(img, os.path.join(dirname,outfile))
except: lprint ('ERROR:  problem while writing results'); exit(1)
lprint ('Successfully written output file '+outfile)
      
#end
exit(0)

