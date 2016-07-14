# force floating point division. Can still use integer with //
from __future__ import division
# This file is used for importing the common utilities classes.
import numpy as np
import matplotlib.pyplot as plt
import sys

sys.path.append("../../../../../")
from FitUtil.WormLikeChain.Python.Code import WLC_Fit
from GeneralUtil.python import PlotUtilities as pPlotUtil

def SafeCast(x):
    if (type(x) is not int) and (type(x) is not float):
        return np.array(list(x)).astype(np.complex128)
    else:
        return complex(x)

def Power(x,y):
    return np.power(SafeCast(x),SafeCast(y))

def Complex(r,i):
    return r + (1j * i)

def Sqrt(x):
    return np.sqrt(SafeCast(x))

def ExtensionPerForce(kbT,Lp,L0,K0,F):
    ToRet = L0*(F/K0 - (-9*kbT - 4*F*Lp)/(12.*kbT) + 
     ((1 + Complex(0,1)*Sqrt(3))*(-9*Power(kbT,2) + 24*F*kbT*Lp - 
          16*Power(F,2)*Power(Lp,2)))/
      (24.*kbT*Power(-243*Power(kbT,3) + 108*F*Power(kbT,2)*Lp - 
          144*Power(F,2)*kbT*Power(Lp,2) + 64*Power(F,3)*Power(Lp,3) + 
          12*Sqrt(3)*Sqrt(135*Power(kbT,6) - 108*F*Power(kbT,5)*Lp + 
             144*Power(F,2)*Power(kbT,4)*Power(Lp,2) - 
             64*Power(F,3)*Power(kbT,3)*Power(Lp,3)),0.3333333333333333)) - 
     ((1 - Complex(0,1)*Sqrt(3))*Power(-243*Power(kbT,3) + 108*F*Power(kbT,2)*Lp - 
          144*Power(F,2)*kbT*Power(Lp,2) + 64*Power(F,3)*Power(Lp,3) + 
          12*Sqrt(3)*Sqrt(135*Power(kbT,6) - 108*F*Power(kbT,5)*Lp + 
             144*Power(F,2)*Power(kbT,4)*Power(Lp,2) - 
             64*Power(F,3)*Power(kbT,3)*Power(Lp,3)),0.3333333333333333))/(24.*kbT))
    # can have complex extension; we just want the real part. 
    return np.real(ToRet)

def run():
    Force = np.linspace(0,250,1000) * 1e-12
    Force += + np.random.rand(Force.size) * 10e-12
    kbT = 4.1e-21
    Lp = 10e-9
    L0 = 650e-9
    K0 = 1000e-12
    ExtSixFifty = ExtensionPerForce(kbT,Lp,L0,K0,Force)
    ExtFiveHundred = ExtensionPerForce(kbT,Lp,500e-9,K0,Force)
    plt.plot(ExtSixFifty,Force)
    plt.plot(ExtFiveHundred,Force,color='g',linestyle='--')
    plt.show()

if __name__ == "__main__":
    run()
