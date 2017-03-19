# force floating point division. Can still use integer with //
from __future__ import division
# This file is used for importing the common utilities classes.
import numpy as np
import matplotlib.pyplot as plt


from FitUtil.EnergyLandscapes.InverseWeierstrass.Python.Code import \
    InverseWeierstrass
from GeneralUtil.python import CheckpointUtilities as pCheckUtil
from GeneralUtil.python import PlotUtilities
from FitUtil.FitUtils.Python import FitUtil as pFitUtil
from Research.Perkins.AnalysisUtil.ForceExtensionAnalysis import \
    FEC_Util,FEC_Plot

class BoundsObj:
    def __init__(self,bounds_folded_nm,bounds_unfolded_nm,
                 bounds_transition_nm,force_one_half_N):
        self.bounds_folded_nm =bounds_folded_nm
        self.bounds_unfolded_nm = bounds_unfolded_nm
        self.bounds_transition_nm =bounds_transition_nm
        self.force_one_half_N = force_one_half_N

class TiltedLandscape:
    def __init__(self,landscape,bounds,kT=4.1e-21):
        """
        Creates a new tilted landscape object

        Args:
            landscape: the IWT landscape object (from InverseWeierstrass)
            bounds: the IWT_Util.BoundsObj
            kT: 1/beta, assumed constant
        """
        bounds_folded_nm = bounds.bounds_folded_nm
        bounds_transition_nm = bounds.bounds_transition_nm
        bounds_unfolded_nm = bounds.bounds_unfolded_nm
        f_one_half = bounds.force_one_half_N 
        self.kT = kT
        self.Landscape_kT =  landscape.EnergyLandscape/kT
        self.Tilted_kT = self.Landscape_kT - \
                         (landscape.Extensions*f_one_half)/kT
        landscape_ext_nm = landscape.Extensions * 1e9
        self.pred_fold_x,self.pred_fold_y,self.coeffs_fold = \
            FitToRegion(landscape_ext_nm,self.Tilted_kT,bounds_folded_nm)
        self.pred_tx_x,self.pred_tx_y,self.coeffs_tx = \
                FitToRegion(landscape_ext_nm,self.Tilted_kT,
                            bounds_transition_nm)
        self.pred_unfold_x,self.pred_unfold_y,self.coeffs_unfold = \
                FitToRegion(landscape_ext_nm,self.Tilted_kT,bounds_unfolded_nm)
        # get the energy landscapes in kT
        # fit a second order to the tilted one (easy to find transition states)
        # get the offsets for the SHO
        self.x0_tx = ExtensionOffsetFromCoeffs(self.coeffs_tx)
        self.x0_fold = ExtensionOffsetFromCoeffs(self.coeffs_fold)
        self.x0_unfold = ExtensionOffsetFromCoeffs(self.coeffs_unfold)
        # DeltaX dagger ais the distance between the transition and
        # folded state. See after equation 5:
        """
        Dudko, O. K., Hummer, G. & Szabo, A. 
        Theory, analysis, and interpretation of single-molecule force 
        spectroscopy experiments. PNAS 105, 1575515760 (2008).
        """
        self.DeltaXDagger = self.x0_tx-self.x0_fold
        """
        DeltaG_dagger is (ibid)
        "the apparent free-energy of activation in the absence of an external
        force." (ie: the height of the energy barrier without force)
        """
        self.ext_idx = np.argmin(np.abs(landscape_ext_nm - self.x0_tx))
        self.MinG = min(self.Landscape_kT)
        self.DeltaGDagger = self.Landscape_kT[self.ext_idx]-self.MinG
        self.landscape_ext_nm = landscape_ext_nm
        self.Offset = np.percentile(self.Tilted_kT,10)
        self.OffsetTilted_kT = self.Tilted_kT-self.Offset

        
def FitToRegion(x,y,bounds_x):
    """
    Fits to a bounded region

    Args:
        x,y: values to fit
        bounds_x: only fit within these bounds
    Returns:
        tuple of <predicted x, predicted y, coefficients>
    """
    GoodIdx = np.where( ( x >= min(bounds_x)) &
                        ( x <= max(bounds_x)) )
    pred_x = x[GoodIdx]
    pred_y = y[GoodIdx]
    coeffs = np.polyfit(x=pred_x,y=pred_y,deg=2)
    return pred_x,pred_y,coeffs
        
def ToIWTObject(o):
    obj = InverseWeierstrass.FEC_Pulling_Object(Time=o.Time,
                                                Extension=o.Separation,
                                                Force=o.Force,
                                                SpringConstant=o.SpringConstant,
                                                Velocity=o.Velocity)
    obj.SetWork(obj.CalculateForceCummulativeWork())
    return obj

def ToIWTObjects(TimeSepForceObjects):
    """
    Converts TimeSepForceObjects to InverseWeierstrass objects

    Args:
        TimeSepForceObjects: list of TimeSepForceObjects to transform
    """
    Objs = [ToIWTObject(o) for o in TimeSepForceObjects]
    return Objs

def ReadInAllFiles(FileNames,Limit,**kwargs):
    """
    Given a list of pxp files, reads them all into a list as 
    TimeSepForce Objcts

    Args:
        FileNames: List of .pxp full paths to data
        Limit: maximum number of curves to return
        kwargs: passed directly to ReadInData
    """
    toRet = []
    for f in FileNames:
        toRet.extend(FEC_Util.ReadInData(f,Limit=Limit,**kwargs))
        # see if we are done
    # only return the limited number we want
    return toRet[:Limit]


def GetIWTObj(Base,FullNames,Force,Limit=150,
              PastZeroExt=60e-9,FilterToMeters=0.25e-9):
    """
    Given files, returns a list of IWT Objects for landscape reconstruction

    Args:
       Base: where to save the cache
       FullNames: ReadInAllFiles 
       Force: if true, force re-creation of the cache
       Limit: maximum number of curves to use 
       PastZeroExt: How much past the surface touchoff to analyze
       FilterToMeters: for determining the surface location, position resolution
       for a savitsky golay filter.
    Returns:
       list of IWT objects, list of TimeSepForce for the retraction, list of 
       TimeSepForce for just the desired amount past touchoff.
    """ 
    mObjs = pCheckUtil.getCheckpoint(Base + "cache.pkl",ReadInAllFiles,
                                     Force,FullNames,Limit)
    ApproachList,RetractList = FEC_Util.BreakUpIntoApproachAndRetract(mObjs)
    # filter all the retractions to the resolution specified,
    # based on the average velocity and the data sampling frequency.
    GetFilter = lambda x: max(3,
                              int((FilterToMeters/x.Velocity)*x.Frequency))
    Touchoff = [FEC_Util.GetFECPullingRegion(r,FilterPoints=GetFilter(r),
                                             MetersAfterTouchoff=PastZeroExt,
                                             Correct=True)
                for r in RetractList]
    # get the IWT transform objects
    IwtObjects = ToIWTObjects(Touchoff)
    return IwtObjects,RetractList,Touchoff

def GetObjectsAndIWT(Base,FullName,Force,NumBins=125,Limit=150):
    """
    Get the IWT and landscape associated with all the force extension curves

    Args:
        See GetIWTObj
    Returns:
        See GetIWTObj, plus it returns a LandscapeObj
    """
    IwtObjects,RetractList,Touchoff = GetIWTObj(Base,FullName,Force,
                                                Limit=Limit)
    # get the IWT
    LandscapeObj = InverseWeierstrass.FreeEnergyAtZeroForce(IwtObjects,
                                                            NumBins=NumBins)
    return IwtObjects,RetractList,Touchoff,LandscapeObj

def toNano(x):
    return x * 1e9


def toPn(x):
    return x * 1e12


def GetAllExtensionsAndForceAndPlot(RetractList,Touchoff,IwtObjects,Base):
    """
    Returns all the normalized forces and extnsions, after plotting them
    in the location given by base

    Args:
        RetractList: list of TimeSepForce objects with just the retract
        Touchoff: list of TimeSepForceObjects *zero-offset* in x and y
        IwtObjects: Touchoff, converted to IWT
        Base: where to save
    Returns:
        tuple of <stage Z,force>
    """
    NPlots = 3
    ext = []
    force = []
    MaxForce = max([np.max(t.Force) for t in Touchoff])
    MinForce = min([np.min(t.Force) for t in Touchoff])
    MaxX = max([t.Zsnsr[-1] for t in Touchoff])
    MaxWork = max([np.max(t.Work) for t in IwtObjects])
    # convert all the maxes to plot-friendly units
    MaxX_nm = MaxX * 1e9
    MaxWork_kbT = MaxWork/(4.1e-21)
    MaxForce_pN = MaxForce * 1e12
    MinForce_pN = MinForce * 1e12
    # get the limits
    ForceLim_pN = [MinForce_pN,MaxForce_pN]
    XLim_nm = [0,MaxX_nm]
    WorkLim_kbT = [0, MaxWork_kbT]
    for i,(Retract,Touch) in enumerate(zip(RetractList,Touchoff)):
        fig = PlotUtilities.figure(figsize=(8,12))
        ForZ = Retract
        RetractZ = Retract.Zsnsr
        RetractZ -= np.min(RetractZ)
        plt.subplot(NPlots,1,1)
        # normalize and flip the force, XXX move to utility...
        ForceRetractPlot = toPn(Retract.Force)
        ForceRetractPlot *= -1
        N = ForceRetractPlot.size
        fraction = 0.2
        ForceRetractPlot -= np.median(ForceRetractPlot[-int(fraction*N):])
        plt.plot(toNano(RetractZ),ForceRetractPlot,alpha=0.3)
        PlotUtilities.lazyLabel("Z stage Position (nm), Absolute","Force (pN)",
                            "Determining Work and Force for FEC")
        plt.subplot(NPlots,1,2)
        Z = Touch.Zsnsr
        plt.plot(toNano(Z),toPn(Touch.Force),alpha=0.3)
        plt.xlim(XLim_nm)
        # force bounds in pN
        plt.ylim(ForceLim_pN)
        PlotUtilities.lazyLabel("","Force (pN)","")
        plt.subplot(NPlots,1,3)
        plt.plot(toNano(Z),IwtObjects[i].Work/(4.1e-21),
                 alpha=0.3)
        plt.xlim(XLim_nm)
        plt.ylim(WorkLim_kbT)
        PlotUtilities.lazyLabel("Z stage Position (nm), relative to touchoff",
                            "Work (kbT)","")
        PlotUtilities.savefig(fig,Base + "{:d}.png".format(i))
        ext.extend(toNano(Z))
        force.extend(toPn(Touch.Force))
    return ext,force

def ForceExtensionHistograms(ext,force,nBins=100,AddAverage=True):
    """
    Makes a 2-d force histogram (ext,force)
    
    Args:
        ext: list of extensions
        force: list of forces
        nBins: how many bins to use
        AddAverage: if true, add average at each bin
    """
    # make a heat map, essentially
    counts, xedges, yedges, Image = plt.hist2d(ext, force,
                                               bins=nBins,cmap='afmhot')
    if (AddAverage):
        x_bins = xedges[:-1]
        y_bins = yedges[:-1]
        bindiff = np.median(np.diff(x_bins))
        for i in range(x_bins.size):
            N = sum(counts[i,:])
            average = (sum(counts[i,:] * y_bins))/N
            label = "Avg binned force" if i == 0 else ""
            plt.plot(x_bins[i]+bindiff/2,average,'go',label=label)
    PlotUtilities.lazyLabel("Separation [nm]",
                            "Force [pN]",
                            "Two-Dimensional Force-Separation Histogram",
                            frameon=True)
    cbar = plt.colorbar()
    cbar.set_label('# in (Force,Separation) Bin', labelpad=10,rotation=270)


def EnergyLandscapePlot(LandscapeObj,FOneHalf=8e-12,
                        ZoomBoundsMeters=[22e-9,30e-9],
                        NumPointsAround=4,
                        stiffness_pN_per_nm=4):
    """
    Plots the enegry landscape  and tilted landscape

    Args:
        LandscapeObj: return from InverseWeierstrass.FreeEnergyAtZeroForce(
        FOneHalf: what to tilt by
    """
    plt.subplot(3,1,1)
    NanoExt =toNano(LandscapeObj.Extensions)
    FreeEnergyEq = LandscapeObj.EnergyLandscape
    plt.plot(NanoExt,FreeEnergyEq * LandscapeObj.Beta)
    PlotUtilities.lazyLabel("","G0",
                            "Reconstructed Energy Landscape")
    plt.subplot(3,1,2)
    TiltedEnergy = (FreeEnergyEq - LandscapeObj.Extensions * FOneHalf)
    TiltedEnergy -= (TiltedEnergy[0])
    EnergyByBeta = TiltedEnergy * LandscapeObj.Beta
    plt.plot(NanoExt,EnergyByBeta)
    PlotUtilities.lazyLabel("Molecular Extension (nm)","G at F-1/2 (kT)",
                            "")
    plt.subplot(3,1,3)
    # zoom on in a specific reigon, just eyeballing
    # the bounds based on the 'full' landscape plot
    ZoomNm = np.array(ZoomBoundsMeters)* 1e9
    WhereZoom = np.where( (NanoExt > min(ZoomNm)) &
                          (NanoExt < max(ZoomNm)))
    if (WhereZoom[0].size > 0):
        EnergyZoom = EnergyByBeta[WhereZoom]
        ExtZoom = NanoExt[WhereZoom]
        # Zero out everything
        MinIdx = np.argmin(EnergyZoom)
        EnergyZoom -= np.min(EnergyZoom)
        ExtZoom -= ExtZoom[MinIdx]
        # fit a parabola to the bottom of the well
        IdxStart = max(0,MinIdx-NumPointsAround)
        IdxEnd = min(ExtZoom.size,MinIdx+NumPointsAround)
        IdxSlice = slice(IdxStart,IdxEnd)
        FitX = ExtZoom[IdxSlice]
        FitY = EnergyZoom[IdxSlice]
        coeffs = np.polyfit(FitX,FitY,deg=2)
        # remove the linear term, which is the second
        coeffs[1] = 0
        xinterp = np.linspace(FitX[0],FitX[-1])
        vals = np.polyval(coeffs,x=xinterp)
        curvature_kbT_per_nm = coeffs[0]
        # note: the well is in kT/nm, so we convert to pN/nm
        # by multuplying by 4.1
        curvature_pN_per_nm = curvature_kbT_per_nm * 4.1
        plt.plot(ExtZoom,EnergyZoom)
        plt.plot(xinterp,vals,color='g',linewidth=4.0,linestyle='--',
                 label="Landscape SHO ({:.1f} pN/nm)".\
                 format(curvature_pN_per_nm))
        plt.plot(xinterp,stiffness_pN_per_nm*xinterp**2+max(vals),
                 label="Cantilever SHO ({:.1f} pN/nm)".\
                 format(stiffness_pN_per_nm))
        plt.ylim([0,max(EnergyZoom)])
        PlotUtilities.lazyLabel("Molecular Extension (nm)","G at F-1/2 (kT)",
                                "",frameon=True)

def set_separation_velocity_by_first_frac(iwt_data,fraction_for_vel):
    """
    Sets the velocity and offset of the given iwt_object by the first
    fraction [0,1] points in iwt_data.Time

    Args:
        see set_separation_velocity_by_first_num, except:
        fraction_for_vel: the fraction [0,1] to use for the fitting
    Returns:
        see set_separation_velocity_by_first_num
    """
    Num = int(np.ceil(iwt_data.Time.size * fraction_for_vel))
    return set_separation_velocity_by_first_num(iwt_data,Num)
                                
                                
def set_separation_velocity_by_first_num(iwt_data,num):
    """
    Sets the velocity and offset of the given iwt_object by the first
    num points in the separation vs time curve

    Args:
        iwt_data: the data to use
        num: the number of points to use
    Returns:
        nothing, but sets the iwt_data offset and velocity
    """
    time_slice = iwt_data.Time[:num]
    sep_slice = iwt_data.Extension[:num]
    coeffs = np.polyfit(x=time_slice,y=sep_slice,deg=1)
    # XXX could just get slope from all, then get offset from np.percentile
    velocity = coeffs[0]
    offset = coeffs[1]
    # adjust the Z function for the fitted velocity and time
    iwt_data.SetVelocityAndOffset(offset,velocity)


def split_into_iwt_objects(d,idx_end_of_unfolding=None,idx_end_of_folding=None,
                           fraction_for_vel=0.5,flip_forces=False):
    """
    given a 'raw' TimeSepForce object, gets the approach and retract 
    as IWT objects, accounting for the velocity and offset of the separation

    Args:
        d: Single TimeSepForce object to split. A single retract/approach
        idx_end_of_unfolding: where the unfolding stops. If not given, we
        assume it happens directly in the middle (ie: default is no 'padding').

        idx_end_of_folding: where unfolding stops. If not given, we assume
        it happens at exactly twice where the folding stops
    
        fraction_for_vel: fit this much of the retract/approach
        separation versus time to determine the true velocity
    returns:
        tuple of <unfolding,refolding> IWT Object
    """
    if (idx_end_of_unfolding is None):
        idx_end_of_unfolding = int(np.ceil(d.Force.size/2))
    if (idx_end_of_folding is None):
        idx_end_of_folding = 2 * idx_end_of_unfolding
    if (flip_forces):
        d.Force *= -1
    # get the unfolding and unfolds
    slice_unfolding = slice(0,idx_end_of_unfolding)
    unfold_tmp = FEC_Util.MakeTimeSepForceFromSlice(d,slice_unfolding)
    slice_folding = slice(idx_end_of_unfolding,idx_end_of_folding)
    fold_tmp = FEC_Util.MakeTimeSepForceFromSlice(d,slice_folding)
    # convert all the unfolding objects to IWT data
    try:
        IwtData = ToIWTObject(unfold_tmp)
        IwtData_fold = ToIWTObject(fold_tmp)
    except AttributeError:
        # Rob messes with the notes
        IwtData = RobTimeSepForceToIWT(unfold_tmp,ZFunc=None)
        IwtData_fold = RobTimeSepForceToIWT(fold_tmp,ZFunc=None)
    # switch the velocities of all ToIWTObject folding objects..
    # set the velocity and Z functions
    set_separation_velocity_by_first_frac(IwtData,fraction_for_vel)
    set_separation_velocity_by_first_frac(IwtData_fold,fraction_for_vel)
    return IwtData,IwtData_fold


def RobTimeSepForceToIWT(o,ZFunc):
    """
    converts a Rob-Walder style pull into a FEC_Pulling_Object

    Args:
         o: TimeSepForce object with Robs meta information
         ZFunc: the z function (schedule) passed along
    Returns:
         properly initialized FEC_Pulling_Object for use in IWT
    """
    # spring constant should be in N/m
    k = o.Meta.__dict__["K"]
    velocity = o.Meta.__dict__["RetractVelocity"]
    Obj = InverseWeierstrass.FEC_Pulling_Object(Time=o.Time,
                                                Extension=o.Separation,
                                                Force=o.Force,
                                                SpringConstant=k,
                                                Velocity=velocity,
                                                ZFunc=ZFunc)
    Obj.SetWork(Obj.CalculateForceCummulativeWork())
    return Obj

        
def ExtensionOffsetFromCoeffs(coeffs):
    """
    Gets the location of the center of the 'well' from the polynomial coeffs

    Args:
        coeffs: the polynomial coefficients, higherst first, from np.polyfit
    Returns:
        x0, from k/2 * (x-x0)**2, fit 
    """
    return -coeffs[1]/(2*coeffs[0])



        
