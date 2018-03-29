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
FIDfile=np.array([])
while answer!="":
   answer = askopenfilename(title="Choose FID file "+str(nfiles+1)+" (press cancel to end)", filetypes=[("FID files","fid")])
   if answer!="":
        answer = os.path.abspath(answer)
        if nfiles==0: 
            FIDfile = np.append(FIDfile, answer)
            nfiles+=1
        elif os.path.getsize(answer) == os.path.getsize(FIDfile[0]):
            FIDfile = np.append(FIDfile, answer)
            nfiles+=1
        else:
            print ('WARNING: rejecting FID input file (size mismatch), try again')
if nfiles==0: print ('ERROR: No FID input file specified'); sys.exit(2)
if nfiles==1: print ('ERROR: Need at least 2 files'); sys.exit(2)
TKwindows.update()
try: win32gui.SetForegroundWindow(win32console.GetConsoleWindow())
except: pass #silent

#read n_Averages from method files
n_Averages=np.empty(shape=(nfiles),dtype=int)
for i in range (0,nfiles):
   METHODfile = os.path.dirname(FIDfile[i])+slash+'method'
   METHODdata = ReadParamFile(METHODfile)
   n_Averages[i]=METHODdata["PVM_NAverages"]

#read ReceiverGain from acqp files
ReceiverGain=np.empty(shape=(nfiles),dtype=int)
for i in range (0,nfiles):
   ACQPfile = os.path.dirname(FIDfile[i])+slash+'acqp'
   ACQPdata = ReadParamFile(ACQPfile)
   ReceiverGain[i]=ACQPdata["RG"]
   
#read FIDs and sum
FIDdata=0
for i in range (0,nfiles):
    print ("Reading file", i)
    with open(FIDfile[i], "r") as f: input_data= np.fromfile(f, dtype=np.int32) 
    input_data = input_data.astype(float)/float(ReceiverGain[i])/float(n_Averages[i])
    FIDdata += input_data
input_data = 0 # free memory
FIDdata /= nfiles
new_RG = int(floor(2147483647./np.amax(np.abs(FIDdata))))
if new_RG<1: new_RG=1
FIDdata *= new_RG # to maintain maximum dynamic range
FIDdata = FIDdata.astype(np.int32)

#make results folder
dirname = os.path.abspath(os.path.dirname(FIDfile[0])+slash+'..'+slash+'AveragingResult')
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
with open(os.path.dirname(FIDfile[0])+slash+'method') as method_file: 
    method_data = method_file.readlines()
with open(new_dirname+slash+'method', "w") as method_file:
    for i in range(0,len(method_data)):
       if method_data[i].startswith('$$ /') and method_data[i].endswith('/method\n'):
          method_data[i] = method_data[i][0:len(method_data[i])-8]+'_sum/method\n'
       method_file.write(method_data[i])
#open/modify/write acqp file
with open(os.path.dirname(FIDfile[0])+slash+'acqp')   as acqp_file:   
    acqp_data = acqp_file.readlines()
with open(new_dirname+slash+'acqp', "w") as acqp_file:
   for i in range(0,len(acqp_data)):
      if acqp_data[i][0:6] == '##$RG=': acqp_file.write('##$RG='+str(new_RG)+'\n')      
      else: acqp_file.write(acqp_data[i])
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