#
# calculates flow volume in [ml/s] along the main flow direction
# (for coronal H_F acquisitions the X component)
#      
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
#

from __future__ import print_function
try: import win32gui, win32console
except: pass #silent
import math
import sys
import os
import zlib
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
    
#intercatively choose input NIFTI file
InputFile = askopenfilename(title="Choose NIFTI file", filetypes=[("NIFTI files",('*.mha','*.nii','*.NII','*.nii.gz','*.NII.GZ'))])
if InputFile=="":print ('ERROR: No input file specified'); sys.exit(2)
InputFile = os.path.abspath(InputFile)
TKwindows.update()
try: win32gui.SetForegroundWindow(win32console.GetConsoleWindow())
except: pass #silent
dirname  = os.path.dirname(InputFile)
basename = os.path.basename(InputFile)
filename = os.path.splitext(InputFile)[0]
extension = os.path.splitext(InputFile)[1].lower()

if extension != ".mha":
    print ('Reading NIFTI file')
    img = nib.load(InputFile)
    data = img.get_data().astype(np.float32)
    SpatResol = img.header.get_zooms()
else:
    print ('Reading MHA file')
    #read MHA header 
    end_header=False
    header_dict = {}
    with open(InputFile, "rb") as f:
        while not end_header:
            line = f.readline()    
            (param_name, current_line) = line.split('=') #split at "=" and strip of spaces
            param_name = param_name.strip()
            current_line = current_line.strip()
            value = ParseSingleValue(current_line)
            header_dict[param_name] = value
            if param_name == 'ElementDataFile': end_header=True
        rawdata = f.read()
        
    #extract relevant parameters from header and check for not implemented stuff
    try: objecttype = header_dict["ObjectType"]
    except: print ('ERROR: Parameter "ObjectType" not found in MHA header'); sys.exit(2);
    if objecttype !='Image': print ('ERROR: ObjectType must be "Image"'); sys.exit(2);
    try: ndim = header_dict["NDims"]
    except: print ('ERROR: Parameter "NDims" not found in MHA header'); sys.exit(2);
    if ndim !=3: print ('ERROR: Parameter "NDims"<>3 not implemented'); sys.exit(2);
    try: binarydata = header_dict["BinaryData"]
    except: print ('ERROR: Parameter "BinaryData" not found in MHA header'); sys.exit(2);
    if binarydata !='True': print ('ERROR: only format with BinaryData implemented'); sys.exit(2);
    try: order = header_dict["BinaryDataByteOrderMSB"]
    except: print ('Warning: Parameter "BinaryDataByteOrderMSB" not found assuming "False"'); order='False'
    if order !='False': print ('ERROR: only format with BinaryDataByteOrderMSB=Flase implemented'); sys.exit(2);
    try: compressed = header_dict["CompressedData"]
    except: print ('Warning: Parameter "CompressedData" not found assuming "False"'); compressed='False'
    if compressed =='True': compressed=True
    else: compressed=False
    try: Resolution = header_dict["ElementSpacing"]
    except: print ('ERROR: Parameter "ElementSpacing" not found in MHA header'); sys.exit(2);
    try:
        SpatResol = np.zeros ((3), dtype=np.float32)
        Resolution  = Resolution.split()
        SpatResol[0] = float (Resolution[0])
        SpatResol[1] = float (Resolution[1])
        SpatResol[2] = float (Resolution[2])
    except: print ('ERROR: Problem parsing parameter "ElementSpacing"'); sys.exit(2);
    try: dims = header_dict["DimSize"]
    except: print ('ERROR: Parameter "DimSize" not found in MHA header'); sys.exit(2);
    try:
        dims  = dims.split()
        dim1 = int (dims[0])
        dim2 = int (dims[1])
        dim3 = int (dims[2])
    except: print ('ERROR: Problem parsing parameter "DimSize"'); sys.exit(2);
    try: veclen = header_dict["ElementNumberOfChannels"]
    except: print ('ERROR: Parameter "ElementNumberOfChannels" not found in MHA header'); sys.exit(2);
    if veclen !=3: print ('ERROR: Parameter "ElementNumberOfChannels"<>3 not implemented'); sys.exit(2);
    try: datatype = header_dict["ElementType"]
    except: print ('ERROR: Parameter "ElementType" not found in MHA header'); sys.exit(2);
    if datatype !='MET_FLOAT': print ('ERROR: ElementType must be "MET_FLOAT"'); sys.exit(2);
    try: datalocation = header_dict["ElementDataFile"]
    except: print ('ERROR: Parameter "ElementDataFile" not found in MHA header'); sys.exit(2);
    if datalocation !='LOCAL': print ('ERROR: Parameter "ElementDataFile" must be "LOCAL"'); sys.exit(2);
    print('.', end='') #progress indicator
    # paramters that are ignored: TransformMatrix, Offset, CenterOfRotation, AnatomicalOrientation, CompressedDataSize

    #decode binary string to floats
    if compressed: rawdata = zlib.decompress(rawdata); print('.', end='') #progress indicator
    if (len(rawdata) % 4) > 0:
        print ("Warning: Data length not a multiple of 4, truncating ....")
        length = int(len(rawdata)/4.0)*4
        rawdata = rawdata[0:length]
    if (len(rawdata)) > dim1*dim2*dim3*veclen*4:
        print ("Warning: Data length larger than expected, truncating ....")
        rawdata = rawdata[0:int(dim1*dim2*dim3*veclen*4)]    
    if (len(rawdata)) < dim1*dim2*dim3*veclen*4:
        print ("ERROR: Data length less than expected")
        sys.exit(2)
    data = np.fromstring (rawdata, dtype=np.float32)
    print('.', end='') #progress indicator    
    data = data.reshape(dim3,dim2,dim1,veclen)
    print('.', end='') #progress indicator
    #find main flow component
    flow_components=np.sum(data[:,:,:,:],axis=(0,1,2))
    main_component = np.argmax(flow_components)
    data = data [:,:,:,main_component]
   
   
#find main flow direction (suposed to be the largest extension of the volume)
flow_directions = np.argsort(data.shape)
flow_directions = flow_directions[::-1] # decreasing order
data = np.transpose(data, flow_directions)

#write flowvolume results
print ('Writing flow rate file')

with open(os.path.join(dirname,filename+'_FlowVolumes.txt'), "w") as text_file:    
    text_file.write("Flow Volumes per slice (X):\n")
    for i in range(0,data.shape[0]): # in our data shape[2] is the main flow direction
        flowvol = np.sum(data[i,:,:])
        flowvol = abs(flowvol)
        flowvol *= 10. # venc is in cm/s, multiply by 10. to get this in mm/s
        flowvol *= SpatResol[1]/1000.*SpatResol[2]/1000. # multiply with inplane spatial resolution, result is in mm^3/s
        flowvol /= 1000. # convert mm^3/s ---> ml/s
        text_file.write("Slice %d:\t%0.2f\tml/s\n" % (i, flowvol))
    text_file.write("\n")        


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