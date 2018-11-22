#
# This tool extracts node data from and saves as textfile fin format 
#        1.0, 1.0, 1.0, 2.0
#        1.0,-1.0,-1.0, 1.5
#       -1.0, 1.0,-1.0, 1.5
#       -1.0,-1.0, 1.0, 1.5
# the first three columns are the X,Y,Z components of the sphere's origin
# the last column contains the radius of the sphere
#      
#
# ----- VERSION HISTORY -----
#
# Version 0.1 - 21, November 2018
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
#    - networkx 
#

from __future__ import print_function
import sys
import os
import math
import gzip
import json
import numpy as np
import networkx as nx
import multiprocessing as mp

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
    
                                         
def loadGraph(filepath):
    with gzip.GzipFile(filepath, 'rb') as file:
        json_str = file.read().decode('utf-8').replace('\0', '')
    data = json.loads(json_str[json_str.find('{')::], encoding='utf-8')
    G = nx.readwrite.json_graph.node_link_graph(data)
    return G

def get_data(nodes, i, prog_dec):
    if i%prog_dec==0: print ('.',end='') # progress indicator
    x = nodes.values()[i]['metadata']['node_coordinates']['x']
    y = nodes.values()[i]['metadata']['node_coordinates']['y']
    z = nodes.values()[i]['metadata']['node_coordinates']['z']
    r = nodes.values()[i]['metadata']['node_squared_radius']
    # in some json files instead of ['node_squared_radius'] use ['node_radius']
    # then comment the below line: data[:,3] = np.sqrt(data[:,3])
    return x,y,z,r;        
    
def worker_get_data(start,end,nodes,prog_dec):    
    return np.asarray(map(lambda i: get_data(nodes,i,prog_dec), range(start,end)))   


    
if __name__ == '__main__':
    mp.freeze_support() #required for pyinstaller 
    cores=mp.cpu_count()-1; cores = max (1,cores) #set number of cores
    
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
    InputFile = askopenfilename(title="Choose (masked) NIFTI file", filetypes=[("List of spheres files",('.json'))])
    if InputFile=="":print ('ERROR: No input file specified'); sys.exit(2)
    InputFile = os.path.abspath(InputFile)
    TKwindows.update()
    try: win32gui.SetForegroundWindow(win32console.GetConsoleWindow())
    except: pass #silent
    dirname  = os.path.dirname(InputFile)
    basename = os.path.basename(InputFile)
    basename = basename[0:basename.rfind('.json')]
    
      
    print ('Loading json file')
    G = loadGraph(InputFile);
    print ('Creating node dictionary')
    nodes = G.node_dict_factory(G.nodes(data=True));
    print ('Cleaning up node dictionary')
    empty_keys = [k for k,v in nodes.iteritems() if not v]
    for k in empty_keys: del nodes[k] # resolve issue of networkx v2.2

    cores_act=min(cores,len(nodes)) # just in case (don't allocate unnecessary cores)
    print ('Multithreading set to %d cores ' % cores_act)   
    print ('Extracting nodes ',end='')
    p = mp.Pool(cores_act)
    return_vals=[]    
    prog_dec = int(len(nodes)/60)+1    
    for i in range(cores_act):
        workpiece=int(math.ceil(float(len(nodes))/float(cores_act)))
        start = i*workpiece
        end   = start+workpiece
        if end > len(nodes): end = len(nodes)       
        return_vals.append(p.apply_async(worker_get_data, args = (start, end, nodes, prog_dec)))
    p.close()
    p.join()
    data = np.zeros((0,4),dtype=int)
    for i in range(0,cores_act):
        data = np.concatenate ((data, return_vals[i].get()), axis=0) # return values

    print ('\nComputing sqrt(r)')
    data[:,3] = np.sqrt(data[:,3])
    print ('Sorting')
    idx = np.argsort(data[:,3])
    idx = idx [::-1] # decreasing order
    data = data[idx,:]
    print ('Writing results')
    np.savetxt (os.path.join(dirname,basename+'.sphrs'),data, delimiter=',', newline='\n')
    print ('done')    
        
        
