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
#    - numexpr
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
import numexpr as ne


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

def sphere2points(params, min_radius): #input [X,Y,Z,Radius],min #output array[[0..#points],[x,y,z,nx,ny,nz]]
    vtk.vtkMultiThreader.SetGlobalMaximumNumberOfThreads(1) # parallelized otherwise
    resolution = 1+int(math.ceil(params[3]/min_radius)**0.5)  
    ico = vtk.vtkPlatonicSolidSource()
    ico.SetSolidTypeToIcosahedron() 
    subdivide = vtk.vtkLoopSubdivisionFilter()
    subdivide.SetNumberOfSubdivisions(resolution)
    subdivide.SetInputConnection(ico.GetOutputPort())
    subdivide.Update()   
    points  = np.asarray(subdivide.GetOutput().GetPoints().GetData()) 
    x = points[:,0]; y = points[:,1]; z = points[:,2]       
    r = ne.evaluate("sqrt(x**2 + y**2 + z**2)")   
    points = np.divide (points [:,:], r[:,None])
    normals = points    
    points = points*params[3] + params[0:3] # scale and translation   
    result = np.concatenate ((points, normals), axis=1)
    return result        

def worker_spheres(start,end,data,prog_dec,min_radius):
    points = np.zeros ((0,6),dtype=np.float32)
    for k in range (start, end):  
        points = np.concatenate ((points, sphere2points(data[k,:],min_radius)), axis=0)
        if k%prog_dec==0: print ('.',end='') # progress indicator          
    return points  

def worker_remove_interior_points (points,data,prog_dec):  
    ne.set_num_threads(1) # parallelized here
    for k in range (data.shape[0]):   
        sphere_r_sq = (0.99*data[k,3])**2.
        sphere_x = data[k,0]; sphere_y = data[k,1]; sphere_z = data[k,2]
        points_x = points[:,0]; points_y = points[:,1]; points_z = points[:,2]
        dist_sq = ne.evaluate("(points_x-sphere_x)**2 + (points_y-sphere_y)**2 + (points_z-sphere_z)**2")
        points = points[dist_sq>sphere_r_sq,:]
        if k%prog_dec==0: print(".", end='')          
    return points
    
    
    
if __name__ == '__main__':
    mp.freeze_support() #required for pyinstaller
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

    #redirect VTK messages
    ow = vtk.vtkOutputWindow()
    ow.SendToStdErrOn()
        
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


    #set number of cores   
    cores=mp.cpu_count()-1; cores = max (1,cores)
    cores_act=min(cores,data.shape[0]) # don't allocate unnecessary cores
    print ('Multithreading set to %d cores ' % cores_act)
    p = mp.Pool(cores_act)
    return_vals=[]        
    #create mesh from spheres
    print ('Creating mesh')
    prog_dec = int(data.shape[0]/(160-cores_act))+1
    for i in range (cores_act):    
        workpiece=int(math.ceil(float(data.shape[0])/float(cores_act)))
        start = i*workpiece
        end   = start+workpiece
        if end > data.shape[0]: end = data.shape[0]
        return_vals.append(p.apply_async(worker_spheres, args = (start, end, data, prog_dec, min_radius)))
    p.close()
    p.join()
    allpoints = np.zeros ((0,6),dtype=np.float32)
    for i in range(cores_act):
        allpoints = np.concatenate ((allpoints, return_vals[i].get()), axis=0) # return values 
    print ('')     
    print("Mesh has %d points" % allpoints.shape[0])

    
    cores_act=min(cores,allpoints.shape[0]) # don't allocate unnecessary cores
    print ('Multithreading set to %d cores ' % cores_act)
    p = mp.Pool(cores_act)
    return_vals=[]        
    #sucessively removing interior points
    print ('Removing interior points')
    prog_dec = int(data.shape[0]/(160-cores_act)*cores_act)+1
    for i in range (cores_act):
        workpiece=int(math.ceil(float(allpoints.shape[0])/float(cores_act)))
        start = i*workpiece
        end   = start+workpiece
        if end > allpoints.shape[0]: end = allpoints.shape[0]
        return_vals.append(p.apply_async(worker_remove_interior_points, args = (allpoints[start:end,:], data, prog_dec)))
    p.close()
    p.join()
    allpoints = np.zeros ((0,6),dtype=np.float32)
    for i in range(cores_act):
        allpoints = np.concatenate ((allpoints, return_vals[i].get()), axis=0) # return values 
    print ('')     
    print("Mesh has %d points" % allpoints.shape[0])        

    
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
                                    allpoints[k,3],allpoints[k,4],allpoints[k,5]))                             
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

    #set VTK multiprocessing
    vtk.vtkMultiThreader.SetGlobalMaximumNumberOfThreads(cores)
    
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