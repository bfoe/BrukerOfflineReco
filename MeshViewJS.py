#
# This tool reads a STL mesh file and generates a simple 
# JavaScript HTML file for visualization 
#      
#
# ----- VERSION HISTORY -----
#
# Version 0.1 - 5, November 2018
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
#    - scipy   (called by trimesh)
#    - trimesh ( https://github.com/mikedh/trimesh )
#    - MeshViewJS_template.html (from this repository)
#

from __future__ import print_function
try: import win32gui, win32console
except: pass #silent
import numpy as np
import sys
import os
import trimesh

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
try: resourcedir = sys._MEIPASS+slash # when on PyInstaller 
except: resourcedir = os.path.abspath(os.path.dirname(sys.argv[0]))+slash; 
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
    
#intercatively choose input STL file
InputFile = askopenfilename(title="Choose STL Mesh file", filetypes=[("Mesh files",('.stl'))])
if InputFile=="":print ('ERROR: No input file specified'); sys.exit(2)
InputFile = os.path.abspath(InputFile)
TKwindows.update()
try: win32gui.SetForegroundWindow(win32console.GetConsoleWindow())
except: pass #silent
dirname  = os.path.dirname(InputFile)
basename = os.path.basename(InputFile)
basename = basename[0:basename.rfind('.stl')] 
template = 'MeshViewJS_template.html'   
    
# Load mesh from file
mesh = trimesh.load(InputFile)
vertices = np.asarray(mesh.vertices)
faces = np.asarray(mesh.faces)
print ('Found %d vertices ' % vertices.shape[0])
print ('Found %d faces    ' % faces.shape[0])
p1=int(vertices.shape[0]/20) # for progress indicator
p2=int(faces.shape[0]/20)    # for progress indicator
vertices /= 1000. # rescale um to mm
xmin=np.amin(vertices[:,0])*1.1; xmax=np.amax(vertices[:,0])*1.1
ymin=np.amin(vertices[:,1])*1.1; ymax=np.amax(vertices[:,1])*1.1
zmin=np.amin(vertices[:,2])*1.1; zmax=np.amax(vertices[:,2])*1.1
if max(abs(xmin),abs(xmax))>400: print ('ERROR: extension in X direction exceed maximum');sys.exit(2)
if max(abs(ymin),abs(ymax))>400: print ('ERROR: extension in Y direction exceed maximum');sys.exit(2)
if max(abs(zmin),abs(zmax))>400: print ('ERROR: extension in Z direction exceed maximum');sys.exit(2)

# read template
try: 
   with open(os.path.join(resourcedir,template)) as f: content = f.readlines()
except: print ('ERROR: reading template file '+template);sys.exit(2)

# write html
print ('Writing html file ',end='')
try: f = open(os.path.join(dirname,basename+'.html'), "w")
except: print ('ERROR: opening output file for write');sys.exit(2)
for i in range(len(content)):
    if content[i].startswith('"viewpoint":'):
       f.write ('"viewpoint":[0,1,0],\n')
    elif content[i].startswith('"xMin":'):
       f.write('"xMin":%.2f,"yMin":%.2f,"zMin":%.2f,\n' % (xmin, ymin, zmin))
    elif content[i].startswith('"xMax":'):
       f.write('"xMax":%.2f,"yMax":%.2f,"zMax":%.2f};\n' % (xmax, ymax, zmax))          
    elif content[i].startswith('"vertices":['):
       f.write('"vertices":[')
       for j in range (vertices.shape[0]):
          f.write (np.array2string(vertices[j],precision=2, separator=',').replace(' ',''))
          if j<vertices.shape[0]-1: f.write(',')
          if (j+1)%10 == 0: f.write('\n')
          if j%p1 == 0: print ('.',end='') # progress bar
       f.write('],\n')   
    elif content[i].startswith('"faces":['):
       f.write('"faces":[')
       for j in range (faces.shape[0]):
          f.write (np.array2string(faces[j],separator=',').replace(' ',''))
          if j<faces.shape[0]-1: f.write(',')
          if (j+1)%10 == 0: f.write('\n')
          if j%p2 == 0: print ('.',end='') # progress bar          
       f.write('],\n')   
    else: f.write(content[i])    
f.close()  
print ('\ndone\n')

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