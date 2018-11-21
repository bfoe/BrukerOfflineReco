#
# This tool reconstructs a surface mesh from a list o spheres
# in the following format: 
#        1.0, 1.0, 1.0, 2.0
#        1.0,-1.0,-1.0, 1.5
#       -1.0, 1.0,-1.0, 1.5
#       -1.0,-1.0, 1.0, 1.5
# the first three columns are the X,Y,Z components of the sphere's origin
# the last column contains the radius of the sphere
# the file should be named <something>.sphrs
# output mesh is written STL format
#      
#
# ----- VERSION HISTORY -----
#
# Version 0.1 - 15, November 2018
#       - 1st public github Release
#
# ----- LICENSE -----                 
#
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the MIT License <https://opensource.org/licenses/MIT>
#
# ----- REQUIREMENTS ----- 
#
#    This program was developed under Python Version 2.7
#    with the following additional libraries: 
#    - numpy
#    - vtk 
#    - numexpr (optional, speeds up some things)
#    - SSDRecon.exe executable from
#      http://www.cs.jhu.edu/~misha/Code/PoissonRecon/Version10.04/
#

from __future__ import print_function
try: import win32gui, win32console
except: pass #silent
import math
import sys
import os
import numpy as np
import subprocess
import vtk
from vtk.util import numpy_support
import multiprocessing as mp

numexpr_installed = True
try: import numexpr as ne
except: numexpr_installed = False


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
        print ('WARNING: Subprocess call returned with error')
        print (command)   
        print (stderr)
        print (stdout)

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
#check for external executables
if not os.path.isfile(os.path.join(resourcedir,'SSDRecon.exe')):
    print ('ERROR:  SSDRecon executable not found '); sys.exit(1)
    
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
InputFile = askopenfilename(title="Choose (masked) NIFTI file", filetypes=[("List of spheres files",('.sphrs'))])
if InputFile=="":print ('ERROR: No input file specified'); sys.exit(2)
InputFile = os.path.abspath(InputFile)
TKwindows.update()
try: win32gui.SetForegroundWindow(win32console.GetConsoleWindow())
except: pass #silent
dirname  = os.path.dirname(InputFile)
basename = os.path.basename(InputFile)
basename = basename[0:basename.rfind('.sphrs')]

#read sphere file
print ('Reading input file') 
data = np.genfromtxt(InputFile, comments='#', delimiter=',')
if data.shape[0]==0: print ("ERROR no data found in input file"); sys.exit(2)
if len(data.shape)==1: # just one sphere
    data = np.reshape(data, (-1, data.shape[0]))    
if data.shape[1]!=4: print ("ERROR inconsistent data in input file"); sys.exit(2)
if len(data.shape)!=2: print ("ERROR inconsistent dimension of input file"); sys.exit(2)
min_radius = np.min (data[:,3])       

# ------------------------- VTK code starts here -----------------------------

#redirect VTK messages
ow = vtk.vtkOutputWindow()
ow.SendToStdErrOn()

#set number of cores   
cores=mp.cpu_count()-1; cores = max (1,cores)
print ('Vtk multithreading set to %d cores ' % cores)   
vtk.vtkMultiThreader.SetGlobalMaximumNumberOfThreads(cores)

def sphere(params): #input [X,Y,Z,Radius] #output mesh
    #
    # vtkSphereSource gives a geodesic aproximation
    # which does not look good, instead ...
    #
    # we aproximate the mesh for the sphere by an icosahedron
    # increasing the resolution with appropriate subdivisions
    #
    resolution = 1+int(math.ceil(params[3]/min_radius)**0.5)
    radius = -1.414*params[3]    
    origin = params[0:3]
    ico = vtk.vtkPlatonicSolidSource()
    ico.SetSolidTypeToIcosahedron()
    subdivide = vtk.vtkLoopSubdivisionFilter()        # set radius = -1.414*params[3] (baseline)
    #subdivide = vtk.vtkButterflySubdivisionFilter()  # set radius = -1.000*params[3] (similar, leaves more faces - bad)
    #subdivide = vtk.vtkLinearSubdivisionFilter()     # set radius = -1.414*params[3] (slightly faster, but bad result)
    subdivide.SetNumberOfSubdivisions(resolution)
    subdivide.SetInputConnection(ico.GetOutputPort())
    scale_transform = vtk.vtkTransform()
    scale_transform.Scale(radius,radius,radius)
    scale_transform.Translate(origin[0]/radius,origin[1]/radius,origin[2]/radius)    
    scaled = vtk.vtkTransformPolyDataFilter()
    scaled.SetInputConnection(subdivide.GetOutputPort())
    scaled.SetTransform(scale_transform)
    scaled.Update()
    #ncells = scaled.GetOutput().GetNumberOfCells();
    #print(ncells)
    return scaled.GetOutput()
    '''
    # simple solution with vtkSphereSource
    resolution = 15+int(math.ceil(params[3]/min_radius))#**1.0)
    source = vtk.vtkSphereSource()
    source.SetCenter(params[0],params[1],params[2])
    source.SetRadius(params[3])
    source.SetThetaResolution (resolution)
    source.SetPhiResolution (resolution)
    source.Update()
    #ncells = source.GetOutput().GetNumberOfCells();
    #print(ncells)
    return source.GetOutput()
    '''

#create mesh from spheres
print ('Creating mesh ',end='')
prog_dec = int(data.shape[0]/20)+1
mesh=sphere(data[0,:])
for k in range (1,data.shape[0]):
    if k%prog_dec==0: print ('.',end='') # progress indicator
    new_sphere = sphere(data[k,:])    
    append = vtk.vtkAppendPolyData()
    append.AddInputData(mesh)
    append.AddInputData(new_sphere)
    append.Update()
    mesh = append.GetOutput()
print ('')    
ncells = mesh.GetNumberOfCells(); npoints = mesh.GetNumberOfPoints()     
if ncells>0: print("Mesh has %d cells and %d points" % (ncells, npoints))
else: print("\nERROR: zero cells in mesh"); sys.exit(1)


#find interior points
print ('Identifying interior points ',end='')
prog_dec = int(data.shape[0]/20)+1
point_remove_array = np.zeros ((0),dtype=np.int64)
allpoints = np.asarray(mesh.GetPoints().GetData())
if numexpr_installed: #40% faster
    ne.set_num_threads(cores)
    for k in range (0,data.shape[0]):
        if k%prog_dec==0: print ('.',end='') # progress indicator
        sphere_origin = data[k,0:3]
        sphere_radius_sq = (0.95*data[k,3])**2.
        allpoints_diff = ne.evaluate("allpoints-sphere_origin")
        x = allpoints_diff[:,0]; y = allpoints_diff[:,1]; z = allpoints_diff[:,2]
        dist_sq = ne.evaluate("x**2+y**2+z**2") 
        temp_array = np.asarray(np.where(ne.evaluate("dist_sq<sphere_radius_sq"))).flatten()
        point_remove_array = np.append(point_remove_array,temp_array)        
else:
    for k in range (0,data.shape[0]):
        if k%prog_dec==0: print ('.',end='') # progress indicator
        sphere_origin = data[k,0:3]
        sphere_radius_sq = (0.95*data[k,3])**2.
        allpoints_diff = allpoints - sphere_origin
        dist_sq = allpoints_diff[:,0]**2.+allpoints_diff[:,1]**2.+allpoints_diff[:,2]**2. 
        temp_array = np.asarray(np.where(dist_sq<sphere_radius_sq)).flatten()
        point_remove_array = np.append(point_remove_array,temp_array)
point_remove_array = np.unique (point_remove_array)
print ('')    

#find all associated cells
print ('Identifying interior cells ',end='')
prog_dec = int(point_remove_array.shape[0]/20)+1
cell_remove_list =[]
for k in range (point_remove_array.shape[0]):
    if k%prog_dec==0: print ('.',end='') # progress indicator
    idlist=vtk.vtkIdList()   
    mesh.GetPointCells(point_remove_array[k], idlist)        
    for l in range (idlist.GetNumberOfIds()):   
        cell_remove_list.append(idlist.GetId(l))
cell_remove_array=np.unique(np.asarray(cell_remove_list))
print ('')

#delete and remove the identified cells
print ('Removing interior cells ',end='')
prog_dec = int(cell_remove_array.shape[0]/20)+1
for k in range (cell_remove_array.shape[0]):
    mesh.DeleteCell(cell_remove_array[k])
    if k%prog_dec==0: print ('.',end='') # progress indicator
print ('')
mesh.RemoveDeletedCells()

#delete unused points
clean1 = vtk.vtkCleanPolyData()
clean1.SetInputData(mesh)
clean1.Update()
mesh = clean1.GetOutput()
ncells = mesh.GetNumberOfCells(); npoints = mesh.GetNumberOfPoints()     
if ncells>0: print("Mesh has %d cells and %d points" % (ncells, npoints))
else: print("\nERROR: zero cells in mesh"); sys.exit(1)

#get points into numpy array
allpoints = np.asarray(mesh.GetPoints().GetData())

#calculate point normals and put into numpy array
Normals= vtk.vtkPolyDataNormals()
Normals.SetInputData(mesh)
Normals.ComputeCellNormalsOff()
Normals.ComputePointNormalsOn()
Normals.SplittingOff() 
Normals.ConsistencyOff()
Normals.AutoOrientNormalsOff()
Normals.Update()
allnormals = np.asarray(Normals.GetOutput().GetPointData().GetNormals())

#just checking
if allpoints.shape!=allnormals.shape:
    print ('Error: Soething went wrong badly while calculation point normals')
    sys.exit(1)

#save point cloud
pointcloud_file=os.path.join(dirname,basename+'_pointcloud.ply')
f = open(pointcloud_file,"w") 
f.write ("ply\n")
f.write ("format ascii 1.0\n")
f.write ("comment homebrew generated\n")
f.write ("element vertex %d\n" % allpoints.shape[0])
f.write ("property float x\n")
f.write ("property float y\n")
f.write ("property float z\n")
f.write ("property float nx\n")
f.write ("property float ny\n")
f.write ("property float nz\n")
f.write ("element face 0\n")
f.write ("property list uchar int vertex_indices\n")
f.write ("end_header\n")
for k in range (allpoints.shape[0]):
    f.write ("%f %f %f %f %f %f\n" % (allpoints[k,0],allpoints[k,1],allpoints[k,2],
                                allnormals[k,0],allnormals[k,1],allnormals[k,2]))                             
f.close()

#run SSDRecon.exe
print ('Reconstructing surface mesh with SSDRecon') 
PLYfile = os.path.join(dirname,basename+'.ply')
command = '"'+os.path.join(resourcedir,'SSDRecon.exe')+'" '
command +='--in "'+pointcloud_file+'" '
command +='--out "'+PLYfile+'" '
command +='--depth 10' 
run (command)
deletefile (pointcloud_file) 

#read PLY
reader = vtk.vtkPLYReader()
reader.SetFileName(PLYfile)
reader.Update()
deletefile (PLYfile)
ncells  = reader.GetOutput().GetNumberOfCells(); 
npoints = reader.GetOutput().GetNumberOfPoints()     
if ncells>0: print("Mesh has %d cells and %d points" % (ncells, npoints))
else: print("\nERROR: zero cells in mesh"); sys.exit(1)

smooth1 = vtk.vtkSmoothPolyDataFilter()
smooth1.SetInputConnection(reader.GetOutputPort())
smooth1.SetNumberOfIterations(7)
smooth1.SetRelaxationFactor(0.1)
smooth1.FeatureEdgeSmoothingOff()
smooth1.BoundarySmoothingOff()
smooth1.Update()

print ('Decimating mesh')
decimate = vtk.vtkQuadricDecimation()
decimate.SetInputConnection(smooth1.GetOutputPort())
decimate.SetTargetReduction(0.5)
decimate.Update()
ncells  = decimate.GetOutput().GetNumberOfCells(); 
npoints = decimate.GetOutput().GetNumberOfPoints()     
if ncells>0: print("Mesh has %d cells and %d points" % (ncells, npoints))
else: print("\nERROR: zero cells in mesh"); sys.exit(1)

smooth2 = vtk.vtkSmoothPolyDataFilter()
smooth2.SetInputConnection(decimate.GetOutputPort())
smooth2.SetNumberOfIterations(7)
smooth2.SetRelaxationFactor(0.1)
smooth2.FeatureEdgeSmoothingOff()
smooth2.BoundarySmoothingOff()
smooth2.Update()

writer = vtk.vtkSTLWriter()
writer.SetInputConnection(smooth2.GetOutputPort())
writer.SetFileTypeToBinary()
writer.SetFileName(os.path.join(dirname,basename+'.stl'))
writer.Write()       


print("done")
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