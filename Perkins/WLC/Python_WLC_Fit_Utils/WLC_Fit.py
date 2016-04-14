# force floating point division. Can still use integer with //
from __future__ import division
# This file is used for importing the common utilities classes.
import numpy as np
import matplotlib.pyplot as plt
import sys


from scipy.optimize import curve_fit
import GeneralUtil.python.GenUtilities as pGenUtil

class WLC_DEF:
    """
    Class defining defaults for inputs. 

    See Wang, 1997.
    """
    L0 = 650e-9 # meters
    Lp = 50e-9 # meters
    K0 = 1200e-12 # Newtons
    kbT = 4.1e-21 # 4.1 pN * nm = 4.1e-21 N*m

class WLC_MODELS:
    """
    Class definining valid models.
    """
    EXTENSIBLE_WANG_1997 = 0
    INEXTENSIBLE_BOUICHAT_1999 = 1


class WlcParamsToVary:
    """
    Class to keep track of what to vary
    """
    def __init__(self,VaryL0=True,VaryLp=False,VaryK0=False):
        """
        Args: 
           VaryL0: If true, contour length is allowed to freely vary
           VaryLp: If true, persistence length is allowed to freely vary
           VaryK0: If true, bulk modulus K0 is allowed to vary
        """
        self.VaryL0 = VaryL0
        self.VaryLp = VaryLp
        self.VaryK0 = VaryK0
    def GetFittingFunctionToCall(self,ToCall,**kwargs):
        """
        Gets the functions to use for fitting
        
        Args: 
           function to call. must follow signature of WlcNonExtensible
        """
        # XXX could make this more efficient ? ...
        if (self.VaryL0 and self.VaryLp and self.VaryK0):
            return lambda x,L0,Lp,K0 : ToCall(ext=x,L0=L0,Lp=Lp,K0=K0,**kwargs)
        elif (self.VaryL0 and self.VaryLp):
            return lambda x,L0,Lp : ToCall(ext=x,L0=L0,Lp=Lp,**kwargs)
        elif (self.VaryL0 and self.VaryK0):
            return lambda x,L0,K0 : ToCall(ext=x,L0=L0,K0=K0,**kwargs)
        elif (self.VaryLp and self.VaryK0):
            return lambda x,Lp,K0 : ToCall(ext=x,Lp=L0,K0=K0,**kwargs)
        elif (self.VaryLp):
            return lambda x,Lp : ToCall(ext=x,Lp=Lp,**kwargs)
        elif (self.VaryL0):
            return lambda x,L0 : ToCall(ext=x,L0=L0,**kwargs)
        elif (self.VaryK0):
            return lambda x,K0 : ToCall(ext=x,K0=K0,**kwargs)
    
class WlcParamValues:
    """
    Class to record parameter values given to a fit or gotten from the same
    """
    def __init__(self,kbT=WLC_DEF.kbT,
                 L0=WLC_DEF.L0,Lp=WLC_DEF.Lp,K0=WLC_DEF.K0):
                 
        """
        Args:
            kbT,Lp,L0 : see WlcPolyCorrect. Initial guesses
            K0: see WlcExtensible. Note this is ignored for non-extensible 
            models
        """
        self.L0 = L0
        self.Lp = Lp
        self.K0 = K0
        self.kbT = kbT
    def GetParamsInOrder(self):
        """
        Conveniene function, gets the parameters in the conventional order
        """
        return self.kbT,self.L0,self.Lp,self.K0

class WlcFitInfo:
    def __init__(self,Model=WLC_MODELS.EXTENSIBLE_WANG_1997,
                 ParamVals=WlcParamValues(),nIters=10,
                 rtol=1e-2,VaryObj=WlcParamsToVary()):
        """
        Args:
        Model: which model, from WLC_MODELS, to use
        ParamValues: values of the parameters for the fit (e.g. initial guesses,
        or final resuts )

        rTol: for extensible models, the relative tolerance between sucessive
        fits before quitting

        nIters: for extensible models, the maximum number of iterations.
        VaryObj: which parameters should be varied for the 
        """
        self.Model = Model
        self.ParamVals = ParamVals
        self.ParamsVaried = VaryObj
        self.nIters = nIters
        self.rtol = rtol
    
def WlcPolyCorrect(kbT,Lp,l):
    """
    From "Estimating the Persistence Length of a Worm-Like Chain Molecule ..."
    C. Bouchiat, M.D. Wang, et al.
    Biophysical Journal Volume 76, Issue 1, January 1999, Pages 409-413

web.mit.edu/cortiz/www/3.052/3.052CourseReader/38_BouchiatBiophysicalJ1999.pdf

    Args:
        kbT : the thermal energy in units of [ForceOutput]/Lp
        Lp  : the persisence length, sensible units of length
        l   : is either extension/Contour=z/L0 Length (inextensible) or   
        z/L0 - F/K0, where f is the force and K0 is the bulk modulus. See 
        Bouchiat, 1999 equation 13
    Returns:
        Model-predicted value for the force
    """
    # parameters taken from paper cited above
    a0=0 
    a1=0
    a2=-.5164228
    a3=-2.737418
    a4=16.07497
    a5=-38.87607
    a6=39.49949
    a7=-14.17718
    #http://docs.scipy.org/doc/numpy/reference/generated/numpy.polyval.html
    #If p is of length N, this function returns the value:a
    # p[0]*x**(N-1) + p[1]*x**(N-2) + ... + p[N-2]*x + p[N-1]
    # note: a0 and a1 are zero, including them for easy of use of polyval.
    polyValCoeffs = [a7,a6,a5,a4,a3,a2,a1,a0]
    inner = 1/(4*(1-l)**2) -1/4 + l + np.polyval(polyValCoeffs,l)
    return kbT/Lp * inner

def WlcNonExtensible(ext,kbT,Lp,L0,*args,**kwargs):
    """
    Gets the non-extensible model for WLC polymers, given and extension

    Args:
        kbT : see WlcPolyCorrect
        Lp: see WlcPolyCorrect
        L0 : contour length, units of ext
        ext: the extension, same units as L0
        *args: extra arguments, ignored (so we can all use the same call sig)
    Returns:
        see WlcPolyCorrect
    """
    return WlcPolyCorrect(kbT,Lp,ext/L0)

def WlcExtensible(ext,kbT,Lp,L0,K0,ForceGuess):
    """
    Fits to the (recursively defined) extensible model. Note this will need to
    be called several times to converge properly

    Args: 
        kbT,Lp,L0,ext : see WlcPolyCorrect
        K0: bulk modulus, units of Force
        ForceGuess: the 'guess' for the force
    Returns:
        see WlcPolyCorrect
    """
    # get the non-extensible model
    return WlcPolyCorrect(kbT,Lp,ext/L0-ForceGuess/K0)

def WlcFit(ext,force,WlcOptions=WlcFitInfo()):
    """
    General fiting function.

    Fits to a WLC model using the given parameters. If Extensible is set to 
    true, then fits an extensible model, subject to running nIters times, 
    or until the relative error between fits is less than rtol, whichever
    is first

    Args:
       ext: the (experimental) extension, 1D array of size N 
       force: the (experimental) force, 1D array of size N 
       Options: WlcFitInfo Object, giving the options for the fit
    Returns: 
       XXX The fitted values
    """
    model = WlcOptions.Model
    kbT,L0,Lp,K0 = WlcOptions.ParamVals.GetParamsInOrder()
    # p0 record our initial guesses; what does the user want use to fit?
    p0 = []
    fixed = dict(kbT=kbT)
    toVary = WlcOptions.ParamsVaried
    # determing what to vary
    if (toVary.VaryL0):
        p0.append(L0)
    else:
        fixed.update(dict(L0=L0))
    if (toVary.VaryLp):
        p0.append(Lp)
    else:
        fixed.update(dict(Lp=Lp))
    if (toVary.VaryK0):
        p0.append(K0)
    else:
        fixed.update(dict(K0=K0))
    # figure out what the model is
    if (model == WLC_MODELS.EXTENSIBLE_WANG_1997):
        # initially, use non-extensible for extensible model, as a first guess
        func = WlcNonExtensible
    elif (model == WLC_MODELS.INEXTENSIBLE_BOUICHAT_1999):
        func = WlcNonExtensible
    else:
        raise TypeError("Didnt recognize model {:s}".format(model))
    # POST: have the functions we want
    mFittingFunc = toVary.GetFittingFunctionToCall(func,**fixed)
    # note: we use p0 as the initial guess for the parameter values
    params,paramsStd,predicted = pGenUtil.GenFit(ext,force,mFittingFunc,p0=p0)
    if (model == WLC_MODELS.EXTENSIBLE_WANG_1997):
        secondFunc = WlcExtensible
        rtol = WlcOptions.rtol
        nIters = WlcOptions.nIters
        for i in range(nIters):
            # get the previous array
            prev = predicted.copy()
            # get rid of bad elements
            prev[np.where(prev<0)] =0
            fixed.update(dict(ForceGuess=prev))
            mFittingFunc = toVary.GetFittingFunctionToCall(secondFunc,**fixed)
            params,paramsStd,predicted = pGenUtil.GenFit(ext,force,
                                                         mFittingFunc,p0=p0)
            close = np.allclose(predicted,prev, rtol=rtol, atol=0,
                                equal_nan=False)
            if (close):
                # then we are close enough to our final result!
                break
    return predicted


def NonExtensibleWlcFit(ext,force,VaryL0=True,VaryLp=False,**kwargs):
    """
    Non-extensible version of the WLC fit. By default, varies the contour length
    to get the fit. Uses Bouichat, 1999 (see aboce) , by default

    Args:
        ext,force : see WlcFit
        VaryL0,VaryLp : see WlcParamsToVary, boolean if we should vary
        **kwargs: passed directly to WlcParamValues (ie: initial guesses)
    Returns:
        see WlcFit
    """
    model = WLC_MODELS.INEXTENSIBLE_BOUICHAT_1999
    mVals = WlcParamValues(**kwargs)
    toVary = WlcParamsToVary(VaryL0=VaryL0,VaryLp=VaryLp)
    # create the informaiton to pass on to the fitter
    mInfo = WlcFitInfo(Model=model,ParamVals=mVals,VaryObj=toVary)
    # call the fitter
    return WlcFit(ext,force,mInfo)

def ExtensibleWlcFit(ext,force,VaryL0=True,VaryLp=False,VaryK0=False,
                     nIters=10,rtol=1e-2,**kwargs):
    """
    extensible version of the WLC fit. By default, varies the contour length
    to get the fit. Uses Bouichat, 1999 (see aboce) , by default

    Args:
        ext,force : see WlcFit
        VaryL0,VaryLp : see WlcParamsToVary, boolean if we should vary
        nIters: number of iterations for the recursively defined function,
        before breaking

        rtol: the relative tolerance
        **kwargs: passed directly to WlcParamValues (ie: initial guesses)
    Returns:
        see WlcFit
    """
    model = WLC_MODELS.EXTENSIBLE_WANG_1997
    mVals = WlcParamValues(**kwargs)
    toVary = WlcParamsToVary(VaryL0=VaryL0,VaryLp=VaryLp,VaryK0=VaryK0)
    mInfo = WlcFitInfo(Model=model,ParamVals=mVals,VaryObj=toVary,
                       nIters=nIters,rtol=rtol)
    return WlcFit(ext,force,mInfo)
        
