import numpy as numpy
from BayesicFitting import Tools
from BayesicFitting import NonLinearModel

__author__ = "BF"
__year__ = 2019
__license__ = "GPL3"
__version__ = "0.9"
__maintainer__ = "BF"
__status__ = "Development"

#  *
#  * This file is part of the BayesicFitting package.
#  *
#  * BayesicFitting is free software: you can redistribute it and/or modify
#  * it under the terms of the GNU Lesser General Public License as
#  * published by the Free Software Foundation, either version 3 of
#  * the License, or ( at your option ) any later version.
#  *
#  * BayesicFitting is distributed in the hope that it will be useful,
#  * but WITHOUT ANY WARRANTY; without even the implied warranty of
#  * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
#  * GNU Lesser General Public License for more details.
#  *
#  * The GPL3 license can be found at <http://www.gnu.org/licenses/>.
#  *
#  * A JAVA version of this code was part of the Herschel Common
#  * Science System (HCSS), also under GPL3.
#  *
#  *    2007 - 2014 Do Kester, SRON (Java code)
#  *    2016 - 2017 Do Kester

class IVIM_Model( NonLinearModel ):
    """
    Exponential Model.
    .. math::
        f( x:p ) = p_0*((1-p_1)*np.exp(-x*p_2) + p_1*np.exp(-x*(p_2+1/p_3)))

    where: 
    p_0 = amplitude, p_1 = perfusion fraction, 
    p_2 = slow diffusion coefficient
    p_3 = inverse fast diffusion coefficient
    As always x = input.

    The parameters are initialized at {1.0, 0.01, 1e-3, 1/1e-2}. It is a non-linear model.

    included ABS to avoid negative results

    Attributes from Model
    ---------------------
        npchain, parameters, stdevs, xUnit, yUnit

    Attributes from FixedModel
    --------------------------
        npmax, fixed, parlist, mlist

    Attributes from BaseModel
    --------------------------
        npbase, ndim, priors, posIndex, nonZero, tiny, deltaP, parNames

    Examples
    --------
    >>> em = IVIM_Model( )
    >>> print( em.getNumberOfParameters( ) )
    4

    Category:    mathematics/Fitting

    """
    def __init__( self, copy=None, **kwargs ):
        """
        Exponential model.
        <br>
        Number of parameters is 4.

        Parameters
        ----------
        copy : IVIM_Model
            to be copied

        """
        param = [1.0, 0.01, 1e-3, 1./1e-2]
        #param = [1.0, 0.01, 1e-3, 1e-2]
        names = ["amplitude","perfusion fraction",
                 "slow diffusion coefficient","inverse fast diffusion coefficient"]
        super( IVIM_Model, self ).__init__( 4, copy=copy, params=param,
                names=names, **kwargs )

    def copy( self ):
        """ Copy method.  """
        return IVIM_Model( copy=self )
        
    def baseResult( self, xdata, params ):
        """
        Returns the result of the model function.

        Parameters
        ----------
        xdata : array_like
            value at which to calculate the result
        params : array_like
            values for the parameters

        """

        params[0]=abs(params[0])
        params[1]=abs(params[1])
        params[2]=abs(params[2])
        params[3]=abs(params[3])   
        e1 = numpy.exp(-xdata* params[2]) 
        e2 = numpy.exp(-xdata*(params[2]+1./params[3]))  
        #e2 = numpy.exp(-xdata*(params[2]+params[3]))        
        return params[0]*((1-params[1])*e1 + params[1]*e2)        

    def basePartial( self, xdata, params, parlist=None ):
        """
        Returns the partials at the input value.

        Parameters
        ----------
        xdata : array_like
            value at which to calculate the result
        params : array_like
            values for the parameters
        parlist : array_like
            list of indices active parameters (or None for all)

        """
        np = self.npbase if parlist is None else len( parlist )
        partial = numpy.ndarray( ( Tools.length( xdata ), np ) )

        #ABS embedded in formula
        #e1 = numpy.exp(-xdata* abs(params[2])) 
        #e2 = numpy.exp(-xdata*(abs(params[2])+1./abs(params[3])))
        #parts = { 0 : ( lambda: params[0]/abs(params[0]) * (1-abs(params[1]))*e1 + abs(params[1])*e2 ),
        #          1 : ( lambda: params[1]/abs(params[1]) * abs(params[0])*(e2-e1) ),
        #          2 : ( lambda: params[2]/abs(params[2]) * -abs(params[0])*xdata*(abs(params[1])*e2 + (1-abs(params[1]))*e1) ),
        #          3 : ( lambda: params[3]  *  abs(params[0])*abs(params[1])*xdata*e2/abs(params[3])**3 )}

        params[0]=abs(params[0])
        params[1]=abs(params[1])
        params[2]=abs(params[2])
        params[3]=abs(params[3])        
        e1 = numpy.exp(-xdata*params[2]) 
        e2 = numpy.exp(-xdata*(params[2]+1./params[3]))
        #e2 = numpy.exp(-xdata*(params[2]+params[3]))
        
        parts = { 0 : ( lambda: (1-params[1])*e1 + params[1]*e2 ),
                  1 : ( lambda: params[0]*(e2-e1) ),
                  2 : ( lambda: -params[0]*xdata*(params[1]*e2 + (1-params[1])*e1) ),
                  3 : ( lambda: params[0]*params[1]*xdata*e2/params[3]**2 )}
                  #3 : ( lambda: params[0]*params[1]*xdata*e2 )}                  
        if parlist is None :
            parlist = range( self.npmax )

        for k,kp in enumerate( parlist ) :
            partial[:,k] = parts[kp]()

        return partial

    def baseDerivative( self, xdata, params ):
        """
        Returns the derivative df/dx at the input value.

        Parameters
        ----------
        xdata : array_like
            value at which to calculate the result
        params : array_like
            values for the parameters

        """
        
        params[0]=abs(params[0])
        params[1]=abs(params[1])
        params[2]=abs(params[2])
        params[3]=abs(params[3])        
        e1 = numpy.exp(-xdata* params[2]) 
        e2 = numpy.exp(-xdata*(params[2]+1./params[3]))
        #e2 = numpy.exp(-xdata*(params[2]+params[3]))        
        return -params[0]*(params[1]*(params[2]+1./params[3])*e2 + (1-params[1])*params[2]*e1)
        #return -params[0]*(params[1]*(params[2]+params[3])*e2 + (1-params[1])*params[2]*e1)        

    def baseName( self ):
        """
        Returns a string representation of the model.

        """
        return str( "IVIM: f( x:p ) = p_0*((1-p_1)*exp(-x*p_2) + p_1*exp(-x*(p_2+1/p_3)))" )

    def baseParameterUnit( self, k ):
        """
        Return the unit of the indicated parameter.
        Parameters: k    parameter number.

        """
        if k == 0:
            return self.yUnit
        return self.xUnit


