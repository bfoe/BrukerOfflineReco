#
# This tool calculates an optimized set of venc values 
# for a multivenc PCA experiment
# The optimization criterion is to maintain a flat signal curve
# over a specified range of possible velocity values
# for the sum of squares of the joint magnitude flow image
#
# This version sweeps the entire parameter space, which
# for a high number of venc's and/or high precision
# can be very time consuming
# Basically this theoretical approach used only to check if 
# the Monte Carlo version of this program finds the global optimum
# For any practical application use the alternative Monte Carlo 
# method implemented in "multivenc_opt_MonteCarlo.py" which is 
# much faster
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
import numpy as np



#setup constants
venc_min = 6  #in cm/s
venc_max = 100 #in cm/s
n_vencs = 7
venc_step = 1
n_velocities = 1000

'''
venc_min   venc_max     n_vencs     venc_step   n_velocities    variation   vencs
4          100          7           2           1000            0.1149937   16 18 22 28 36 54 100
6          100          7           2           1000            0.0956260   16 18 22 28 36 54 100

CAUTION: the below line with venc_step=1 take several days to calculate            
6          100          7           1           1000            0.08        16 18 22 27 36 54 100 (hope this is the result)
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


n=0
end = False
opt_value = sys.float_info.max
venc_temp = np.zeros(n_vencs,dtype=np.float32)
venc_index = np.zeros(n_vencs,dtype=int)
vencs_shape = int((venc_max-venc_min)/venc_step)+1
vencs = np.linspace (venc_min,venc_max,vencs_shape, endpoint=True)
velocities = np.linspace (venc_min,venc_max,n_velocities, endpoint=True)
#print (vencs)
#print (velocities)
for i in range (n_vencs): venc_index[i]=vencs.shape[0]-i-1
while not end:
      n += 1
      #check if improoved
      opt_min_value = opt_min (vencs[venc_index],velocities)
      if opt_min_value<opt_value:
          opt_value = opt_min_value
          opt_venc = vencs[venc_index]
      #decrement venc      
      venc_index[n_vencs-1] -= 1      
      for i in range (1,n_vencs-1):
         if venc_index[n_vencs-i]<i-1: 
            venc_index[n_vencs-i-1] -= 1
            if i==n_vencs-3: print ('-',end='') #progress
            elif i==n_vencs-2: print ('|')      #progress
      for i in range (n_vencs-2,0,-1):
         if venc_index[n_vencs-i]<i-1:        
            venc_index[n_vencs-i] = venc_index[n_vencs-i-1]-1
      #detect end      
      if venc_index[n_vencs-1]<0: end = True 

print ("")         
print ("Total iterations : ",n)
print ("Signal variation : ",opt_min (opt_venc,velocities))
print ("Optimum vencs    : ",np.sort (opt_venc) )

