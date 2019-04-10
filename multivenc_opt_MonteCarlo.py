#
# This tool calculates an optimized set of venc values 
# for a multivenc PCA experiment
# Otimization is done with Monte Carlo (simulated anealing)
# The optimization criterium is to maintain a flat signal curve
# over a specified range of possible velocity values
# for the sum of squares of the joint magnitude flow image
#
#
# ----- VERSION HISTORY -----
#
# Version 0.1 - 10, April 2019
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


from __future__ import print_function
import sys
import os
import math
import random
import signal
import numpy as np

#setup constants
venc_min = 6  #in cm/s
venc_max = 100 #in cm/s
n_vencs = 7
precision = 1e-0

# this defines the abort criterium of the inner loop     
max_n_nochange = 500*int(np.sqrt(1/precision))

# this defines the abort criterium of the outer loop
# with the outer loop we try to avoid local minima  
max_m_nochange = 100 
                       
# this defines the #points in the velocity range venc_min..venc_max
# that will be optimized
#
# start with 10000 and decrease until geting unstable
# unstability is seen by 
# a) different optimized venc resutls  
# b) higher signal variation                      
n_velocities = 1000


'''
venc_min venc_max n_vencs precision max_n max_m n_vel Signal_Variation  Optimum_vencs
5        100      4       1e-0      500   200   1000  0.13465219630223  28   36   54   100
5        100      5       1e-0      500   200   1000  0.11045434389920  22   27   36   53   100
5        100      6       1e-0      500   200   1000  0.09838418700704  18   21   27   35   53  100

1        100      7       1e-0      500   200   1000  0.20199153239086   7    8   18   21   26   34   100 (bad)
2        100      7       1e-0      500   200   1000  0.16776329037842   8   18   21   26   35   52   100 (reasonable)
3        100      7       1e-0      500   200   1000  0.13409696245506  16   18   22   27   36   54   100 (OK)
4        100      7       1e-0      500   200   1000  0.10617855081583  16   18   22   27   36   54   100 (OK)
5        100      7       1e-0      500   200   1000  0.09058257241690  16   18   22   27   36   54   100 (OK)
6        100      7       1e-0      500   200   1000  0.08486399313936  16   18   22   27   36   54   100 (OK)

1        100      7       1e-1      500   200   1000  0.19101596431492   0.1  6.9  7.3 10.3 36.8 53.8 100 (bad)
2        100      7       1e-1      500   200   1000  0.15976038921775  10.5 11.5 12.7 27.1 35.3 52.4 100 (reasonable)
3        100      7       1e-1      500   200   1000  0.12439550809659  15.4 17.7 21.1 26.3 35.0 52.3 100 (OK)
4        100      7       1e-1      500   200   1000  0.09583690004737  15.4 17.7 21.1 26.3 35.0 52.5 100 (OK)
5        100      7       1e-1      500   200   1000  0.07990270758201  15.5 17.8 21.2 26.5 35.2 52.7 100 (OK) 
6        100      7       1e-1      500   200   1000  0.07431882303518  15.5 17.8 21.2 26.5 35.2 52.8 100 (OK)
'''

                      
#
# this is the function that is to be minimized
# to find the best venc combination
#
def opt_min (vencs,velocities):
    signal = 0
    # sum of squares
    for i in range (vencs.shape[0]):  
        signal += np.square(np.sin(velocities/vencs[i]*np.pi))
    signal = np.sqrt(signal)
    return np.std(signal)

def randrange_float(start, stop, step):
    return np.round(random.randint(0, int((stop - start) / step)) * step + start, decimals=digits)    

def signal_handler(signal, frame): 
    print ('\nUser abort detected\n')        
    try: print ("Total iterations : ",total_iter)
    except: pass # silent
    try: print ("Signal variation : ",final_opt_value)
    except: print ("No result yet")
    try: print ("Optimum vencs    : ",final_opt_venc)
    except: print ("No result yet")
    print ("")
    wait_for_keypress()
    sys.exit(0)    

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

#setup signal handler
signal.signal(signal.SIGINT, signal_handler)  # keyboard interrupt
signal.signal(signal.SIGTERM, signal_handler) # kill/shutdown
if  'SIGHUP' in dir(signal): signal.signal(signal.SIGHUP, signal_handler)  # shell exit (linux)


digits = int(math.ceil(-math.log10(precision)))
velocities = np.linspace (venc_min,venc_max,n_velocities, endpoint=True)
m=0
m_last=0
total_iter=0
final_opt_value = sys.float_info.max
while m-m_last<max_m_nochange:
    m += 1
    vencs = np.linspace (venc_min,venc_max,n_vencs, endpoint=True)
    vencs = np.round(vencs, decimals=digits)
    vencs = np.unique(vencs)
    #print ('velocities =', velocities)
    #print (0, opt_min (vencs,velocities), vencs)     
    n=0
    n_last=0
    opt_value = sys.float_info.max
    progress_indicator = int(max_n_nochange/5.)
    while n-n_last<max_n_nochange:
        total_iter += 1
        n += 1    
        #MonteCarlo modify a random venc by a random value
        venc_index = random.randint(0, n_vencs-2)  # allow all but last to be modified  
        rnd_value = sys.float_info.max    # set an invald value
        p=0; p_max = 100
        while (vencs[venc_index]+rnd_value in vencs or\
              vencs[venc_index]+rnd_value<=0 or\
              vencs[venc_index]+rnd_value>venc_max) and\
              p<p_max:
           p += 1              
           rnd_value = randrange_float(-vencs[venc_index],vencs[venc_index],precision)        
        if p<p_max: # valid random value found
           temp_vencs = np.copy(vencs)          
           temp_vencs[venc_index] += rnd_value       
           #check if improoved        
           opt_min_value = opt_min (temp_vencs,velocities)
           if opt_min_value<opt_value:
              n_last = n
              temp_vencs = np.sort(temp_vencs)
              opt_venc  = temp_vencs
              opt_value = opt_min_value              
              vencs = temp_vencs # next attempt start from here   
              #print (n, opt_value, opt_venc)
        else: n = max_n_nochange+n_last+1 # no valid random value found, abort the whole attempt    
        if n%progress_indicator==0: print ('.',end='') # progress indicator
    if opt_value<final_opt_value:
        m_last = m  
        final_opt_venc  = opt_venc
        final_opt_value = opt_value
        #print (n, final_opt_value, final_opt_venc)
        print ('|',end='') # progress indicator
    else: print (';',end='') # progress indicator
    
print ("")         
print ("Total iterations : ",total_iter)
print ("Signal variation : ",final_opt_value)
print ("Optimum vencs    : ",final_opt_venc)
print ("")

wait_for_keypress()