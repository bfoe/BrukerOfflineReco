#
# reads Bruker MR data (Paravision v5.1)
#
# average 2 or more fid files
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
#    - matplotlib.pylab (optional)
#


from __future__ import print_function
try: import win32gui, win32console
except: pass #silent
from math import floor
import sys
import os
import shutil
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

def ReadParamFile(filepath):
    global OrigFilename;
    #Read a Bruker MRI experiment's method or acqp file to a dictionary.
    param_dict = {}
    with open(filepath, "r") as f:
        while True:
            line = f.readline()
            if not line:
                break
            # when line contains parameter
            if line.startswith('$$ /'):
                OrigFilename = line[line.find('/nmr/')+5:]
                OrigFilename = OrigFilename[0:len(OrigFilename)-8]
                OrigFilename = OrigFilename.replace(".", "_")
                OrigFilename = OrigFilename.replace("/", "_")
            if line.startswith('##$'):
                (param_name, current_line) = line[3:].split('=') #split at "="
                # if current entry (current_line) is arraysize
                if current_line[0:2] == "( " and current_line[-3:-1] == " )":
                    value = ParseArray(f, current_line)
                # if current entry (current_line) is struct/list
                elif current_line[0] == "(" and current_line[-3:-1] != " )":
                    # if neccessary read in multiple lines
                    while current_line[-2] != ")":
                        current_line = current_line[0:-1] + f.readline()
                    # parse the values to a list
                    value = [ParseSingleValue(x)
                             for x in current_line[1:-2].split(', ')]
                # otherwise current entry must be single string or number
                else:
                    value = ParseSingleValue(current_line)
                # save parsed value to dict
                param_dict[param_name] = value
    return param_dict
    
def ParseArray(current_file, line):
    # extract the arraysize and convert it to numpy
    line = line[1:-2].replace(" ", "").split(",")
    arraysize = np.array([int(x) for x in line])
    # then extract the next line
    vallist = current_file.readline().split()
    # if the line was a string, then return it directly
    try:
        float(vallist[0])
    except ValueError:
        return " ".join(vallist)
    # include potentially multiple lines
    while len(vallist) != np.prod(arraysize):
        vallist = vallist + current_file.readline().split()
    # try converting to int, if error, then to float
    try:
        vallist = [int(x) for x in vallist]
    except ValueError:
        vallist = [float(x) for x in vallist]
    # convert to numpy array
    if len(vallist) > 1:
        return np.reshape(np.array(vallist), arraysize)
    # or to plain number
    else:
        return vallist[0]

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
FIDfile1 = askopenfilename(title="Choose first Bruker FID file", filetypes=[("FID files","fid")])
if FIDfile1 == "": print ('ERROR: First FID input file not specified'); sys.exit(2)
FIDfile1 = os.path.abspath(FIDfile1) 
FIDfile2 = askopenfilename(title="Choose second Bruker FID file", filetypes=[("FID files","fid")])
if FIDfile2 == "": print ('ERROR: Second FID input file not specified'); sys.exit(2)
FIDfile2 = os.path.abspath(FIDfile2) 
TKwindows.update()
try: win32gui.SetForegroundWindow(win32console.GetConsoleWindow())
except: pass #silent

#read n_Averages from method files
METHODfile1 = os.path.dirname(FIDfile1)+slash+'method'
METHODdata1 = ReadParamFile(METHODfile1)
n_Averages1=METHODdata1["PVM_NAverages"]
METHODfile2 = os.path.dirname(FIDfile2)+slash+'method'
METHODdata2 = ReadParamFile(METHODfile2)
n_Averages2=METHODdata2["PVM_NAverages"]

#read ReceiverGain from acqp files
ACQPfile1 = os.path.dirname(FIDfile1)+slash+'acqp'
ACQPdata1 = ReadParamFile(ACQPfile1)
ReceiverGain1=ACQPdata1["RG"]
ACQPfile2 = os.path.dirname(FIDfile2)+slash+'acqp'
ACQPdata2 = ReadParamFile(ACQPfile2)
ReceiverGain2=ACQPdata2["RG"]
   
#read FIDs and subtract
FIDdata=0
print ("Reading file 1")
with open(FIDfile1, "r") as f: input_data1= np.fromfile(f, dtype=np.int32) 
input_data1 = input_data1.astype(float)/float(ReceiverGain1)/float(n_Averages1)
print ("Reading file 2")
with open(FIDfile2, "r") as f: input_data2= np.fromfile(f, dtype=np.int32) 
input_data2 = input_data2.astype(float)/float(ReceiverGain2)/float(n_Averages2)
FIDdata = input_data1 - input_data2
input_data1 = 0 # free memory
input_data2 = 0 # free memory
max_result = np.amax(np.abs(FIDdata))
if max_result == 0: print ('ERROR: Result is zero (subtracting identical images ?)'); sys.exit(2)
new_RG = int(floor(2147483647./max_result))
if new_RG<1: new_RG=1
FIDdata *= new_RG # to maintain maximum dynamic range
FIDdata = FIDdata.astype(np.int32)

#make results folder
dirname = os.path.abspath(os.path.dirname(FIDfile1)+slash+'..'+slash+'SubtractionResult')
new_dirname = dirname
i=0
while os.path.exists(new_dirname):
   i+=1
   new_dirname = dirname+'('+str(i)+')'
try: os.makedirs(new_dirname)
except: print ('ERROR: unable to make folder', new_dirname); sys.exit(2)
print ("Saving results")
#write averaged fid file   
FIDdata.tofile(new_dirname+slash+"fid")
#open/modify/write method file
with open(os.path.dirname(FIDfile1)+slash+'method') as method_file: 
    method_data = method_file.readlines()
with open(new_dirname+slash+'method', "w") as method_file:
    for i in range(0,len(method_data)):
       if method_data[i].startswith('$$ /') and method_data[i].endswith('/method\n'):
          method_data[i] = method_data[i][0:len(method_data[i])-8]+'_sbtr/method\n'
       if method_data[i].startswith('##$PVM_NAverages'):
          method_data[i] = '##$PVM_NAverages=1\n'          
       method_file.write(method_data[i])
#open/modify/write acqp file
with open(os.path.dirname(FIDfile1)+slash+'acqp')   as acqp_file:   
    acqp_data = acqp_file.readlines()
with open(new_dirname+slash+'acqp', "w") as acqp_file:
   for i in range(0,len(acqp_data)):
      if acqp_data[i][0:6] == '##$RG=': acqp_file.write('##$RG='+str(new_RG)+'\n')      
      else: acqp_file.write(acqp_data[i])
#write logfile      
with open(new_dirname+slash+'Logfile.txt', "w") as logfile:
    logfile.write('Result fid file is the complex average of:\n')
    for i in range (0,nfiles):
       logfile.write(FIDfile1+'\n')
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