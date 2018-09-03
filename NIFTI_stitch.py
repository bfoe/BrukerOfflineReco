#
# This tool reads two NIFTI files and stiches them together along the largest dimension
# using simple correlation as similarity criterium and (bruteforce) searching over a 
# subset of possible translations (3 degrees  of freedom)
#      
#
# ----- VERSION HISTORY -----
#
# Version 0.1 - 08, August 2018
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
#    - scipy
#    - nibabel
#

from __future__ import print_function
try: import win32gui, win32console
except: pass #silent
import math
import sys
import os
import numpy as np
import nibabel as nib
import multiprocessing as mp
from multiprocessing import Pool



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
  
def smooth(x,window_len):
    w=np.hanning(window_len)
    s=np.r_[2*x[0]-x[window_len-1::-1],x,2*x[-1]-x[-1:-window_len:-1]]
    w=np.hanning(window_len)
    y=np.convolve(w/w.sum(),s,mode='same')
    return y[window_len:-window_len+1]   
    
def anlyze_histogram (d):
    d=d.flatten()     
    n_points=d.shape[0]
    steps=int(np.sqrt(n_points)); start=0; fin=np.max(d)
    xbins =  np.linspace(start,fin,steps)
    ybins, binedges = np.histogram(d, bins=xbins)
    ybins = np.resize (ybins,len(xbins)); ybins[len(ybins)-1]=0
    xbins = xbins[1:-1] # through out first point (lots of counts)
    ybins = ybins[1:-1] # through out first point (lots of counts)
    nz = np.nonzero (ybins)
    xbins = xbins[nz] # through out histogram points with zero count (basically the initial points)
    ybins = ybins[nz] # through out histogram points with zero count (basically the initial points)
    ybins = smooth(ybins,ybins.shape[0]/20)
    #find minimum
    start=ybins.argmax()
    i=start;x_min=0;y_min=ybins[start]
    while i<len(ybins):
        i+=1
        if ybins[i]<=y_min: y_min=ybins[i]; x_min=i; 
        else: i=len(ybins);
    minimum=xbins[x_min]
    #find maximum
    start=x_min
    i=start;x_max=0;y_max=ybins[start]
    while i<len(ybins):
        i+=1
        if ybins[i]>y_max: y_max=ybins[i]; x_max=i; 
        else: i=len(ybins);
    maximum=xbins[x_max]
    i=len(ybins)-1
    points=0; threshold = 0
    while i>0:
       points += ybins[i]
       threshold = xbins[i]   
       if points >  n_points*0.02: #2%
          i=1
       i -= 1     
    return minimum, maximum, threshold

def worker(data1,data2,stitch_start,stitch_end,overlap,roll1_search_range,roll2_search_range):
    goodness = np.zeros ((2*(stitch_end-stitch_start), 2*roll1_search_range+1, 2*roll2_search_range+1),dtype=np.float64)
    for i in range (stitch_start,stitch_end):
       print ('.',end='')
       for j in range (0,2*roll1_search_range+1):
           for k in range (0,2*roll2_search_range+1):
               #0
               test1 = data1[:,data1.shape[1]-i-overlap:data1.shape[1]-i,:]          
               test2 = data2[:,i:i+overlap,:]
               test2 = np.roll(test2,j-roll1_search_range,axis=0)
               test2 = np.roll(test2,k-roll2_search_range,axis=2)
               test1=test1.flatten().astype(np.float64); test2=test2.flatten().astype(np.float64)
               goodness [2*(i-stitch_start),j,k] = np.correlate(test1,test2)
               #1
               test1 = data1[:,data1.shape[1]-i-overlap:data1.shape[1]-i,:]          
               test2 = data2[:,i+1:i+overlap+1,:] # add 1
               test2 = np.roll(test2,j-roll1_search_range,axis=0)
               test2 = np.roll(test2,k-roll2_search_range,axis=2)          
               test1=test1.flatten().astype(np.float64); test2=test2.flatten().astype(np.float64)
               goodness [2*(i-stitch_start)+1,k] = np.correlate(test1,test2)           
    return goodness

    
if __name__ == '__main__':
    mp.freeze_support()    
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
    InputFile1 = askopenfilename(title="Choose 1st NIFTI file", filetypes=[("NIFTI files",('*.nii','*.NII','*.nii.gz','*.NII.GZ'))])
    if InputFile1=="":print ('ERROR: No input file specified'); sys.exit(2)
    InputFile1 = os.path.abspath(InputFile1)
    InputFile2 = askopenfilename(title="Choose 2nd NIFTI file", filetypes=[("NIFTI files",('*.nii','*.NII','*.nii.gz','*.NII.GZ'))])
    if InputFile2=="":print ('ERROR: No input file specified'); sys.exit(2)
    InputFile2 = os.path.abspath(InputFile2)
    TKwindows.update()
    try: win32gui.SetForegroundWindow(win32console.GetConsoleWindow())
    except: pass #silent
    dirname  = os.path.dirname(InputFile1)
    basename = os.path.basename(InputFile1)
    filename = os.path.splitext(InputFile1)[0]
    filename = os.path.splitext(filename)[0]
    filename_connected = filename+'_Stitched.nii.gz'


    print ('Reading NIFTI files')
    img1 = nib.load(InputFile1)
    data1 = img1.get_data().astype(np.float32)
    SpatResol1 = np.asarray(img1.header.get_zooms())
    img2 = nib.load(InputFile2)
    data2 = img2.get_data().astype(np.float32)
    SpatResol2 = np.asarray(img2.header.get_zooms())


    # histogram correction
    min1, max1, high1 = anlyze_histogram(data1)
    min2, max2, high2 = anlyze_histogram(data2)
    factor = high2/high1
    data2 /= factor
    print ('Histogram correction: %0.2f' % factor)


    #some checks
    if not np.array_equal(data1.shape,data2.shape): print ('ERROR: image dimension mismatch'); sys.exit(2)
    if not np.array_equal(SpatResol1,SpatResol2): print ('ERROR: image resolution mismatch'); sys.exit(2)
    directions = np.argsort(data1.shape)
    directions = directions[::-1] # decreasing order
    if directions[0]!=1: print ('ERROR: largest dimension is not index 1, not implemented'); sys.exit(2)

    #initialize match search
    overlap=int(0.05*data1.shape[1]) #5%
    overlap=int(overlap/2)*2 # make it even
    print ('Overlap is ',overlap)
    stitch_search_range = int(0.25*data1.shape[1]) #25%
    print ('Shift search range is 0..'+str(2*stitch_search_range))
    roll1_search_range  = int(0.05*data1.shape[0]) #5%
    print ('Roll1 search range is -'+str(roll1_search_range)+'..'+str(roll1_search_range))
    roll2_search_range  = int(0.05*data1.shape[0]) #5%
    print ('Roll2 search range is -'+str(roll2_search_range)+'..'+str(roll2_search_range))
    #find match

    cores=mp.cpu_count()
    print ('optimizing using',cores,'cores ',end='')
    p = Pool(cores)
    return_vals=[]
    for i in range(0,cores):
        workpiece=int(math.ceil(float(stitch_search_range)/float(cores)))
        stitch_start = i*workpiece
        stitch_end   = stitch_start+workpiece
        if stitch_end > stitch_search_range: stitch_end = stitch_search_range       
        return_vals.append(p.apply_async(worker, args = (data1,data2,stitch_start,stitch_end,overlap,roll1_search_range,roll2_search_range)))
    p.close()
    p.join() 
    goodness = np.zeros ((0, 2*roll1_search_range+1, 2*roll2_search_range+1),dtype=np.float64)
    for i in range(0,cores):   
        goodness = np.concatenate ((goodness, return_vals[i].get()), axis=0)
    print ('')
    max_=0; max_i=0; max_j=0; max_k=0
    for i in range (0,2*stitch_search_range):
       for j in range (0,2*roll1_search_range+1):
           for k in range (0,2*roll2_search_range+1):
               if goodness [i,j,k]>max_: 
                  max_ = goodness [i,j,k]
                  max_i=i
                  max_j=j-roll1_search_range 
                  max_k=k-roll2_search_range
    print ('Optimal shift found at',max_i)
    print ('Optimal roll1 found at',max_j)
    print ('Optimal roll2 found at',max_k)
                  
    #initialize join
    crop1 = int(max_i/2)
    crop2 = max_i-crop1
    crop1 += int(overlap/2)
    crop2 += int(overlap/2)
    data2 = np.roll (data2,max_j, axis=0)
    data2 = np.roll (data2,max_k, axis=2)
    data = np.zeros ((data1.shape[0],data1.shape[1]+data2.shape[1]-crop1-crop2,data1.shape[2]),dtype=np.float32)
    hanning1 = np.zeros(overlap,dtype=np.float32)
    hanning2 = np.zeros(overlap,dtype=np.float32)
    x_ = np.linspace (0,np.pi/2,num=overlap+2); x_ = x_[1:overlap+1]
    # do join
    data [:,0:data1.shape[1]-crop1,:] = data1[:,0:data1.shape[1]-crop1,:]
    data [:,data1.shape[1]-crop1:-1,:] = data2[:,crop2:-1,:]
    test1 = data1[:,data1.shape[1]-crop1-overlap/2:data1.shape[1]-crop1+overlap/2,:]
    test2 = data2[:,crop2-overlap/2:crop2+overlap/2,:]
    hanning1 [0:overlap] = np.power(np.cos(x_),2)
    hanning2 = 1-hanning1
    data [:,data1.shape[1]-crop1-overlap/2:data1.shape[1]-crop1+overlap/2,:] = test1[:,:,:]*hanning1[None,:,None] + test2[:,:,:]*hanning2[None,:,None]

        
    #transform to int 
    max_data = np.amax(data);
    data *= 32767./max_data
    data = data.astype(np.int16)
    #save NIFTI
    print ('Writing output files')
    aff = np.eye(4)
    aff[0,0] = SpatResol1[0]*1000; aff[0,3] = -(data.shape[0]/2)*aff[0,0]
    aff[1,1] = SpatResol1[1]*1000; aff[1,3] = -(data.shape[1]/2)*aff[1,1]
    aff[2,2] = SpatResol1[2]*1000; aff[2,3] = -(data.shape[2]/2)*aff[2,2]
    NIFTIimg = nib.Nifti1Image(data, aff)
    NIFTIimg.header.set_slope_inter(max_data/32767.,0)
    NIFTIimg.header.set_xyzt_units(3, 8)
    NIFTIimg.set_sform(aff, code=0)
    NIFTIimg.set_qform(aff, code=1)
    try:
        nib.save(NIFTIimg, os.path.join(dirname,filename_connected))
    except:
        print ('\nERROR:  problem while writing results'); sys.exit(1)
    print ('If the result looks wrong, try choosing the input files in inverse order')   
    print ('\nSuccessfully written output file')   
            


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