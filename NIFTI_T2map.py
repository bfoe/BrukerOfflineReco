#
# calculates maps of T2 decay from 4D NIFTI
#      
#
# ----- VERSION HISTORY -----
#
# Version 0.1 - 20, February 2019
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
from math import ceil, floor
import sys
import os
import math
import binascii
import warnings
import numpy as np
import nibabel as nib
import multiprocessing as mp
if getattr( sys, 'frozen', False ): # running as pyinstaller bundle
   from scipy_extract import zoom
   from scipy_extract import median_filter
   from scipy_extract import gaussian_filter 
   from scipy_extract import curve_fit 
else: # running native python
   from scipy.ndimage import zoom 
   from scipy.ndimage import median_filter 
   from scipy.ndimage import gaussian_filter 
   from scipy.optimize import curve_fit   


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

#define fit functions
def expdecay(x,a,T2):
    return a* np.exp(-x/T2)
def FIT (x,y,T2_clip,R2_clip): 
    bad_fit1=False
    initial_conditions=[np.max(y), (np.max(x)+np.min(x))/2]
    try: 
        pars1,covar1 = curve_fit(expdecay,x,y,p0=initial_conditions,maxfev=100*(len(x)+1))
        if  (np.size(pars1)==1): bad_fit1=True
        if (np.size(covar1)==1): bad_fit1=True
    except RuntimeError: bad_fit1=True        
    if bad_fit1: T2=0;T2err=0; A=0; Aerr=0                                            #fit failed
    elif pars1[1]>T2_clip: T2=T2_clip; T2err=covar1[1,1]; A=pars1[0]; Aerr=covar1[0,0] #above clip value
    elif pars1[1]<R2_clip: T2=R2_clip; T2err=covar1[1,1]; A=pars1[0]; Aerr=covar1[0,0] #above clip value
    else: T2=pars1[1]; T2err=covar1[1,1]; A=pars1[0]; Aerr=covar1[0,0]                 #fit OK
    return T2, T2err, A, Aerr

def worker_curvefit(TE,IMGdata,p,T2_clip,R2_clip):
    warnings.filterwarnings("ignore")
    data_T2map = np.zeros((IMGdata.shape[0]),dtype=np.float32)
    data_R2map = np.zeros((IMGdata.shape[0]),dtype=np.float32)    
    for i in range(IMGdata.shape[0]):
       if i%p==0: print ('.',end='')
       #filter     
       zeros = np.nonzero (IMGdata[i,:]==0)[0]
       zeros = np.append  (zeros,IMGdata.shape[1])      
       IMGdata[i,zeros[0]:IMGdata.shape[1]]  = 0
       #extract
       nz = np.nonzero (IMGdata [i,:])     
       TE_temp = TE[nz]
       IMGdata_temp = IMGdata [i,:][nz]
       #fit
       if TE_temp.shape[0]>5:
          data_T2map [i], T2err, A, Aerr = FIT (TE_temp, IMGdata_temp, T2_clip,R2_clip)
          if data_T2map[i]>0: data_R2map [i] = 1000./data_T2map [i]
    return data_T2map, data_R2map

def worker_zoom(IMGdata,p):
    z=2
    zoomed = np.zeros((IMGdata.shape[0]*z,IMGdata.shape[1]*z,IMGdata.shape[2],IMGdata.shape[3]),dtype=np.float32)
    for i in range(IMGdata.shape[3]):
       if i%p==0: print ('.',end='')
       zoomed[:,:,:,i] = zoom(IMGdata[:,:,:,i],[z,z,1],order=2)
    return zoomed
    
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
    InputFile = askopenfilename(title="Choose NIFTI file", filetypes=[("NIFTI files",('MAGNT.nii.gz'))])
    if InputFile=="":print ('ERROR: No input file specified'); sys.exit(2)
    InputFile = os.path.abspath(InputFile)
    TKwindows.update()
    try: win32gui.SetForegroundWindow(win32console.GetConsoleWindow())
    except: pass #silent
    dirname  = os.path.dirname(InputFile)
    basename = os.path.basename(InputFile)
    basename = os.path.splitext(basename)[0]
    basename = basename[0:basename.rfind('_MAG')]
    outfile1 = basename+'_T2map.nii.gz'
    outfile2 = basename+'_R2map.nii.gz'

    print ('Reading NIFTI file')
    img = nib.load(InputFile)
    IMGdata = img.get_data().astype(np.float32)
    SpatResol = np.asarray(img.header.get_zooms())
    Shape = np.asarray(img.header.get_data_shape())
    affine = img.affine
    xyzt_units1  = img.header.get_xyzt_units()[0]
    xyzt_units2  = img.header.get_xyzt_units()[1]
    sform = int(img.header['sform_code'])
    qform = int(img.header['qform_code'])
    TE_str = img.header['descrip'] 
    try: TE = np.fromstring(binascii.a2b_uu(str(TE_str)[6:]),dtype=np.int16)/100.
    except: print ("ERROR: parsing TE information from header"); sys.exit(1)


    #check some stuff
    if len(IMGdata.shape) != 4:
        print ('ERROR: 4D NIFTI expected, 3D NIFTI found, this is probably not a multiecho file'); 
        sys.exit(1)
    if TE.shape[0]!=IMGdata.shape[3]:
        print ("ERROR: number of TE's in header unequal data dimension "); 
        sys.exit(1)
    if np.unique(TE).shape[0]<5:
        print ("ERROR: need at least 5 unique TE's"); 
        sys.exit(1)
    if np.amin(IMGdata)<0 or abs(np.amax(IMGdata)-np.pi)<0.2:
        print ("ERROR: this looks like a Phase Image"); 
        sys.exit(1)
        
    # check if already masked
    N=10 # use 10% at the corners of the FOV
    tresh=np.empty(shape=8,dtype=np.float)
    avg=np.empty(shape=8,dtype=np.float)
    std=np.empty(shape=8,dtype=np.float)
    xstart=0; xend=int(IMGdata.shape[0]/N)
    ystart=0; yend=int(IMGdata.shape[1]/N)
    zstart=0; zend=int(ceil(float(IMGdata.shape[2])/float(N)))
    arr=np.abs(IMGdata[xstart:xend,ystart:yend,zstart:zend,:])
    if np.count_nonzero(arr==0)>0.2*arr.size: #masked input image not OK
        print ("ERROR: this looks like an image with noise mask"); 
        sys.exit(1)    

        
    # --------------------------- start doing something --------------------------

    #correct first echo
    #
    #problem description:
    # In 2D sequences used for T2 quantification the signal amplitude
    # of the first echo often is considerably below what would be expected for
    # an exponential decay. Actually it's not the first echo that is too low, 
    # but the following echoes that are higher, this due to an additional
    # contribution arising from stimulated echos at the borders of the slice 
    # selection profile where the refocusing angle is below 180. RF pulses
    # wit well behaved slice profiles like "hermitian" minimize the problem, 
    # not so well behaved "gauss" pulses aggravate the effect. 3D sequences
    # intrinsically avoid this kind of problem altogether. T2* sequences do not
    # refocus echos by RF, so no problem there either even not if T2* is 2D.
    #
    #elaborate solutions to this problem:
    # https://www.ncbi.nlm.nih.gov/pubmed/22012743
    # https://www.ncbi.nlm.nih.gov/pubmed/24648387
    # 
    #How the below algorithm worx:
    # the solutions proposed in the above publications require information of
    # the applied pulse sequence and it's parameters. We implemented a generic
    # solution based on a simple correction of the amplitude of the first echo,
    # which consists of the following steps:
    # 1) analysis of the histogram of the last echo image to find the 1% 
    #    highest intensity points (these are are the regions with the longest T2).
    # 2) calculate the T2 fit from the average of these points.
    # 3) back calculate the signal of the first echo, compare to the actual 
    #    value and calculate a correction factor.
    # 4) apply the correction to the first echo image only.
    # in step 2) we also estimate the normalized RMS error and only apply 4)
    # under certain quality constraints.    
    #
    #
    # analyze histogram to find 1% higest intensity points of last echo image
    n_points=IMGdata.shape[0]*IMGdata.shape[1]*IMGdata.shape[2]
    steps=np.sqrt(n_points); start=0; fin=np.max(IMGdata [:,:,:,-1])
    xbins =  np.linspace(start,fin,steps)
    ybins, binedges = np.histogram(IMGdata [:,:,:,-1], bins=xbins)
    i=len(ybins)-1
    points=0; thresh = 0
    while i>0:
       points += ybins[i]       
       thresh = xbins[i]   
       if points >  n_points*0.01: #1%
          i=1
       i -= 1
    #get longT2 points       
    longT2_idx = IMGdata [:,:,:,-1]>thresh    
    longT2_data = np.average(IMGdata[longT2_idx,:],axis=0)
    #fit longT2 points
    longT2_value, T2err, A, Aerr = FIT (TE[1:-1], longT2_data[1:-1], np.max(TE)*100, np.min(TE)*0.01)
    #estimate normalized root mean square error
    rmse=0
    for i in range(1,TE.shape[0]):
       rmse += ((expdecay(TE[i],A,longT2_value)-longT2_data[i])/longT2_data[i])**2
    rmse = np.sqrt(rmse/(TE.shape[0]-1))
    #do correction
    first_echo_expected = expdecay(TE[0],A,longT2_value)
    corr = first_echo_expected/longT2_data[0]    
    if rmse>0.05: #5% RMS error limit
       print ('1st echo correction skipped (RMS error too high)')
    elif corr<=1:
       print ('1st echo correction skipped (correction factor < 1.0)')    
    else: #do correction
       print ('1st echo correction factor is %0.2f' % corr)
       IMGdata [:,:,:,0] *= corr

    #increase resolution 2x
    #IMGdata = zoom(IMGdata,[2,2,1,1],order=2) #simple non-mp code
    cores_min = min (cores,IMGdata.shape[3])    
    print ('Increase resolution 2x using',cores_min,'cores')
    p = mp.Pool(cores_min)
    return_vals=[]
    progress_tag = int(math.ceil(IMGdata.shape[3]/70.))
    for i in range(cores_min):
        workpiece=int(math.ceil(float(IMGdata.shape[3])/float(cores_min)))
        start = i*workpiece
        end   = start+workpiece
        if end > IMGdata.shape[3]: end = IMGdata.shape[3]  
        return_vals.append(p.apply_async(worker_zoom, args = (IMGdata[:,:,:,start:end], progress_tag)))
    p.close()
    p.join()    
    #get results  
    IMGdata = return_vals[0].get()     
    for i in range(1,cores_min):
        IMGdata = np.concatenate ((IMGdata, return_vals[i].get()),axis=3)
    print ('') 

    # calculate mask
    # use noise in all 8 corners to establish threshold
    N=10 # use 10% at the corners of the FOV
    std_factor = 2 # thresh = avg + std_factor*std
    thresh=np.empty(shape=8,dtype=np.float)
    avg=np.empty(shape=8,dtype=np.float)
    std=np.empty(shape=8,dtype=np.float)
    xstart=0; xend=int(IMGdata.shape[0]/N)
    ystart=0; yend=int(IMGdata.shape[1]/N)
    zstart=0; zend=int(ceil(float(IMGdata.shape[2])/float(N)))
    arr=np.abs(IMGdata[xstart:xend,ystart:yend,zstart:zend,:])
    avg[0]=np.mean(arr)
    std[0]=np.std(arr)
    thresh[0]=avg[0] + std_factor*std[0]
    xstart=int(IMGdata.shape[0]-IMGdata.shape[0]/N); xend=IMGdata.shape[0]
    ystart=0; yend=int(IMGdata.shape[1]/N)
    zstart=0; zend=int(ceil(float(IMGdata.shape[2])/float(N)))
    arr=np.abs(IMGdata[xstart:xend,ystart:yend,zstart:zend,:])
    avg[1]=np.mean(arr)
    std[1]=np.std(arr)
    thresh[1]=avg[1] + std_factor*std[1]
    xstart=0; xend=int(IMGdata.shape[0]/N)
    ystart=int(IMGdata.shape[1]-IMGdata.shape[1]/N); yend=IMGdata.shape[1]
    zstart=0; zend=int(ceil(float(IMGdata.shape[2])/float(N)))
    arr=np.abs(IMGdata[xstart:xend,ystart:yend,zstart:zend,:])
    avg[2]=np.mean(arr)
    std[2]=np.std(arr)
    thresh[2]=avg[2] + std_factor*std[2]
    xstart=int(IMGdata.shape[0]-IMGdata.shape[0]/N); xend=IMGdata.shape[0]
    ystart=int(IMGdata.shape[1]-IMGdata.shape[1]/N); yend=IMGdata.shape[1]
    zstart=0; zend=int(ceil(float(IMGdata.shape[2])/float(N)))
    arr=np.abs(IMGdata[xstart:xend,ystart:yend,zstart:zend,:])
    avg[3]=np.mean(arr)
    std[3]=np.std(arr)
    thresh[3]=avg[3] + std_factor*std[3]
    xstart=0; xend=int(IMGdata.shape[0]/N)
    ystart=0; yend=int(IMGdata.shape[1]/N)
    zstart=int(floor(float(IMGdata.shape[2])-float(IMGdata.shape[2])/float(N))); zend=IMGdata.shape[2]
    arr=np.abs(IMGdata[xstart:xend,ystart:yend,zstart:zend,:])
    avg[4]=np.mean(arr)
    std[4]=np.std(arr)
    thresh[4]=avg[4] + std_factor*std[4]
    xstart=int(IMGdata.shape[0]-IMGdata.shape[0]/N); xend=IMGdata.shape[0]
    ystart=0; yend=int(IMGdata.shape[1]/N)
    zstart=int(floor(float(IMGdata.shape[2])-float(IMGdata.shape[2])/float(N))); zend=IMGdata.shape[2]
    arr=np.abs(IMGdata[xstart:xend,ystart:yend,zstart:zend])
    avg[5]=np.mean(arr)
    std[5]=np.std(arr)
    thresh[5]=avg[5] + std_factor*std[5]
    xstart=0; xend=int(IMGdata.shape[0]/N)
    ystart=int(IMGdata.shape[1]-IMGdata.shape[1]/N); yend=IMGdata.shape[1]
    zstart=int(floor(float(IMGdata.shape[2])-float(IMGdata.shape[2])/float(N))); zend=IMGdata.shape[2]
    arr=np.abs(IMGdata[xstart:xend,ystart:yend,zstart:zend,:])
    avg[6]=np.mean(arr)
    std[6]=np.std(arr)
    thresh[6]=avg[6] + std_factor*std[6]
    xstart=int(IMGdata.shape[0]-IMGdata.shape[0]/N); xend=IMGdata.shape[0]
    ystart=int(IMGdata.shape[1]-IMGdata.shape[1]/N); yend=IMGdata.shape[1]
    zstart=int(floor(float(IMGdata.shape[2])-float(IMGdata.shape[2])/float(N))); zend=IMGdata.shape[2]
    arr=np.abs(IMGdata[xstart:xend,ystart:yend,zstart:zend,:])
    avg[7]=np.mean(arr)
    std[7]=np.std(arr)
    thresh[7]=avg[7] + std_factor*std[7]
    mask_treshold=np.min(thresh)  
    
    # reshape flatten
    dim = IMGdata.shape
    IMGdata = np.reshape(IMGdata, (dim[0]*dim[1]*dim[2],dim[3]))  

    #apply mask
    mask =  IMGdata > mask_treshold 
    IMGdata = IMGdata*mask

    # vector reduction
    mask = np.sum(mask, axis=1)
    mask = mask >= 5 # at least 5 good echoes
    IMGdata = IMGdata [mask,:]
    
    #T2map calculation  
    print ('Start fitting using', cores, 'cores')
    #set up multiprocessing
    p = mp.Pool(cores)
    return_vals=[]
    T2_clip_factor=5    # factor of max(TE) where to start clipping T2
    R2_clip_factor=1/3. # factor of min(TE) where to start clipping R2
    T2_clip = np.max(TE)*T2_clip_factor
    R2_clip = np.min(TE)*R2_clip_factor
    progress_tag = IMGdata.shape[0]/70
    for i in range(cores):
        workpiece=int(math.ceil(float(IMGdata.shape[0])/float(cores)))
        start = i*workpiece
        end   = start+workpiece
        if end > IMGdata.shape[0]: end = IMGdata.shape[0]  
        return_vals.append(p.apply_async(worker_curvefit, args = (TE, IMGdata[start:end,:], progress_tag, T2_clip, R2_clip)))
    p.close()
    p.join()    
    #get results  
    temp_T2map = np.zeros((0),dtype=np.float32)
    temp_R2map = np.zeros((0),dtype=np.float32)        
    for i in range(cores):
        temp_T2map = np.concatenate ((temp_T2map, return_vals[i].get()[0]))
        temp_R2map = np.concatenate ((temp_R2map, return_vals[i].get()[1]))
    print ('')        

    #undo vector reduction
    data_T2map = np.zeros((dim[0]*dim[1]*dim[2]),dtype=np.float32)
    data_R2map = np.zeros((dim[0]*dim[1]*dim[2]),dtype=np.float32)     
    data_T2map[mask] = temp_T2map[:] # undo vector reduction 
    data_R2map[mask] = temp_R2map[:] # undo vector reduction    
    
    #reshape to original
    data_T2map = np.reshape(data_T2map, (dim[0],dim[1],dim[2]))
    data_R2map = np.reshape(data_R2map, (dim[0],dim[1],dim[2]))

    #filter 2D
    print ('Filter T2 map')
    for i in range(dim[2]):
       data_T2map[:,:,i] = median_filter  (data_T2map[:,:,i], size = (5,5))
       data_T2map[:,:,i] = gaussian_filter(data_T2map[:,:,i], sigma=0.7)
       data_R2map[:,:,i] = median_filter  (data_R2map[:,:,i], size = (5,5))
       data_R2map[:,:,i] = gaussian_filter(data_R2map[:,:,i], sigma=0.7)

    #decrease resolution
    print ('Decrease resolution 2x')
    data_T2map = zoom(data_T2map,[0.5,0.5,1],order=1)
    data_R2map = zoom(data_R2map,[0.5,0.5,1],order=1)
      
    #transform to int
    max_T2 = np.amax(data_T2map);
    data_T2map *= 32767./max_T2
    data_T2map = data_T2map.astype(np.int16)
    max_R2 = np.amax(data_R2map);
    data_R2map *= 32767./max_R2
    data_R2map = data_R2map.astype(np.int16)

    #Write NIFTIs
    T2_img = nib.Nifti1Image(data_T2map, affine)
    T2_img.header.set_slope_inter(max_T2/32767.,0)
    T2_img.header.set_xyzt_units(xyzt_units1,xyzt_units2)
    T2_img.set_sform(affine, sform)
    T2_img.set_qform(affine, qform)
    R2_img = nib.Nifti1Image(data_R2map, affine)
    R2_img.header.set_slope_inter(max_R2/32767.,0)
    R2_img.header.set_xyzt_units(xyzt_units1,xyzt_units2)
    R2_img.set_sform(affine, sform)
    R2_img.set_qform(affine, qform)
    try: 
       nib.save(T2_img, os.path.join(dirname,outfile1))
       nib.save(R2_img, os.path.join(dirname,outfile2))   
    except: print ('ERROR:  problem while writing results'); exit(1)
    print ('Successfully written output files')

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