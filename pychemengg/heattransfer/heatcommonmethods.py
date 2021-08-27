# -*- coding: utf-8 -*-

# This file is part of pyChemEngg python package.
 
# PyChemEngg: A python-based framework to promote problem solving and critical
# thinking in chemical engineering.

# Copyright (c) 2021 Harvinder Singh Gill <profhsgill@gmail.com>

# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:

# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.

# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.

"""Module contains functions commonly used in heat transfer.
"""

import numpy as np
from scipy import interpolate
from scipy.interpolate import interp1d
import bisect
import warnings
warnings.simplefilter('always')

__all__ = ["calc_LMTD", "calc_internalenergychange",
           "calc_overallheattransfercoefficient", "calc_Re"]

def calc_internalenergychange(mass=None, specificheat=None, deltaT=None):
    r"""To compute change in internal energy due to temperature change.
    
    
    Parameters
    ----------
    mass : `int or float`
        Mass of object (solid, liquid or gas).
    specificheat : `int or float`
        Specific heat of object.
    deltaT : `int or float`
        Change in temperature of object.

    Returns
    -------
    internal energy change: `int or float`
        Internal energy change due to change in temperature of object.
    
    
    Notes
    -----
    The following formula is used:         
    
    .. math::
        
       \Delta U = mC_p\Delta T
        
           
    *where:*
        
        m = change in internal energy
        
        :math:`C_p` = specific heat
        
        :math:`\Delta T` = change in temperature
        
        :math:`\Delta U` = change in internal energy

    
    Examples
    --------
    First import the module **heatcommonmethods**.
       
    >>> from pychemengg.heattransfer import heatcommonmethods as hcm 
    >>> hcm.calc_internalenergychange(mass=4.686, specificheat=0.395,
                                      deltaT=150-100)
    92.5485

    
    References
    ----------
    [1] G. F. Nellis and S. A. Klein, "Introduction to Engineering 
    Heat Transfer", 1st Edition. Cambridge University Press, 2021.
    
    [2] Y. A. Cengel and A. J. Ghajar, "Heat And Mass Transfer
    Fundamentals and Applications", 6th Edition. New York, McGraw Hill
    Education, 2020.
    
    [3] T.L. Bergman, A. S. Lavine, F. P. Incropera, D. P. Dewitt, 
    "Fundamentals of Heat and Mass Transfer", 7th Edition, 
    John Wiley, 2011.
    
    """    
    return mass * specificheat * deltaT


def calc_Re(characteristic_length=None, velocity=None, density=None, viscosity=None):
    r"""To compute Reynolds number for flowing fluid.
    
    
    Parameters
    ----------
    characteristic_length : `int or float`
        Relevant length for computing Reynolds number (see Notes).
    velocity : `int or float`
        Velocity of fluid.
    density : `int or float`
        Density of fluid.
    viscosity : `int or float`
        Viscosity of fluid.

    Returns
    -------
    Reynolds number : `int or float`
        Reynolds number for a flowing fluid.
    
    
    Notes
    -----
    The following formula is used:         
    
    .. math::
        
       Re = \frac {L_c V \rho} {\mu}
       
           
    *where:*
        
        :math:`L_c`  *= characteristic length*
            for example:
                :math:`L_c` *= internal diameter of pipe/tube for fluid flowing*
                in a circular tube/pipe
                :math:`L_c` *= external diameter for fluid flowing
                outside a circular tube/pipe*
                :math:`L_c` *= hydraulic diameter for fluid flowing
                in a non-circular tube/pipe*

        *V = fluid velocity*
        
        :math:`\rho` *= density of fluid*
        
        :math:`\mu` *= viscosity of fluid*

    
    Examples
    --------
    First import the module **heatcommonmethods**.
       
    >>> from pychemengg.heattransfer import heatcommonmethods as hcm 
    >>> hcm.calc_Re(characteristic_length=0.015, velocity=6.43,
                    density=1.059, viscosity=2.008e-5)
    5086.680776892429

    
    References
    ----------
    [1] G. F. Nellis and S. A. Klein, "Introduction to Engineering 
    Heat Transfer", 1st Edition. Cambridge University Press, 2021.
    
    [2] Y. A. Cengel and A. J. Ghajar, "Heat And Mass Transfer
    Fundamentals and Applications", 6th Edition. New York, McGraw Hill
    Education, 2020.
    
    [3] T.L. Bergman, A. S. Lavine, F. P. Incropera, D. P. Dewitt, 
    "Fundamentals of Heat and Mass Transfer", 7th Edition, 
    John Wiley, 2011.
    """ 
    return characteristic_length * velocity * density / viscosity


def calc_Pr(viscosity=None, specificheat=None, thermalconductivity=None):
    r"""To compute Prandtl number for a fluid.
    
    
    Parameters
    ----------
    viscosity : `int or float`
        Viscosity of fluid.
    specificheat : `int or float`
        Specific heat of fluid.
    thermalconductivity : `int or float`
        Thermal conductivity of fluid.


    Returns
    -------
    Prandtl number : `int or float`
        Prandtl number for a fluid.
    
    
    Notes
    -----
    The following formula is used:         
    
    .. math::
        
       Pr = \frac {\mu C_p} {k}
       
           
    *where:*
        
        :math:`\mu` *= viscosity of fluid*        
        
        :math:`C_p`  *= specific heat*
        
        *k = thermal conductivity*

    
    Examples
    --------
    First import the module **heatcommonmethods**.
       
    >>> from pychemengg.heattransfer import heatcommonmethods as hcm 
    >>> hcm.calc_Pr(viscosity=0.8374, specificheat=1880, thermalconductivity=0.145)
    10857.324137931037

    
    References
    ----------
    [1] G. F. Nellis and S. A. Klein, "Introduction to Engineering 
    Heat Transfer", 1st Edition. Cambridge University Press, 2021.
    
    [2] Y. A. Cengel and A. J. Ghajar, "Heat And Mass Transfer
    Fundamentals and Applications", 6th Edition. New York, McGraw Hill
    Education, 2020.
    
    [3] T.L. Bergman, A. S. Lavine, F. P. Incropera, D. P. Dewitt, 
    "Fundamentals of Heat and Mass Transfer", 7th Edition, 
    John Wiley, 2011.
    """ 
    return  viscosity * specificheat / thermalconductivity


def calc_LMTD(deltaT1=None, deltaT2=None):
    r"""To compute log mean temperature difference (LMTD) 
    
    Parameters
    ----------
    deltaT1 : `int or float`
        Temperature difference.
    deltaT2 : `int or float`
        Temperature difference.
        

    Returns
    -------
    LMTD : `int or float`
        Reynolds number for a flowing fluid.
    
    
    Notes
    -----
    The following formula is used.
    
    .. math::
        
        LMTD = \cfrac {\Delta T_1 - \Delta T_2} {ln[ \Delta T_1 / \Delta T_2 }
        
           
    *where:*
        
        :math:`\Delta T_1` *= temperature difference #1*
        
        :math:`\Delta T_2` *= temperature difference #2*
        
        :math:`\Delta T_1` and :math:`\Delta T_2` vary on the type 
        of problem under consideration

    Internal to this function, the temperature differences are converted into 
    complex numbers to allow for the use of this function in situations
    where temperature must be iteratively solved using solvers
    such as 'fsolve' from scipy. The value that is returned is the real part.
    
    
    Examples
    --------
    First import the module **heatcommonmethods**.
       
    >>> from pychemengg.heattransfer import heatcommonmethods as hcm 
    >>> hcm.calc_LMTD(deltaT1=90-40, deltaT2=65-20)
    47.456107905149494

    
    References
    ----------
    [1] G. F. Nellis and S. A. Klein, "Introduction to Engineering 
    Heat Transfer", 1st Edition. Cambridge University Press, 2021.
    
    [2] Y. A. Cengel and A. J. Ghajar, "Heat And Mass Transfer
    Fundamentals and Applications", 6th Edition. New York, McGraw Hill
    Education, 2020.
    
    [3] T.L. Bergman, A. S. Lavine, F. P. Incropera, D. P. Dewitt, 
    "Fundamentals of Heat and Mass Transfer", 7th Edition, 
    John Wiley, 2011.
    """ 
    if deltaT1 == deltaT2:
        LMTD = deltaT1
        return LMTD
    else:    
        deltaT1 = deltaT1 + 0j
        deltaT2 = deltaT2 + 0j
        LMTD = (deltaT1-deltaT2) / np.log(deltaT1/deltaT2)
        LMTD = LMTD.real
        return LMTD


def calc_overallheattransfercoefficient(inner_heattransfercoefficient=None,
                          outer_heattransfercoefficient=None,
                          inner_foulingfactor=None,
                          outer_foulingfactor=None,
                          resistanceof_wall=None,
                          inner_area=None,
                          outer_area=None,
                          is_basedon_outerarea=False,
                          is_basedon_innerarea=False):
    r"""To compute overall heat transfer coefficient (U) from resistance
    values
    
    Parameters
    ----------
    inner_heattransfercoefficient : `int or float`
        Heat transfer coefficient on inside/oneside of object.
    outer_heattransfercoefficient : `int or float`
        Heat transfer coefficient on outside/otherside of object.
    inner_foulingfactor: `int or float`
        Fouling factor on inside/oneside of object.
    outer_foulingfactor: `int or float`
        Fouling factor on outside/otherside of object.
    resistanceof_wall: `int or float`,
        Resistance of conduction of object material (based on shape).
    inner_area: `int or float`
        Surface area of inside/oneside of object.
    outer_area: `int or float`
        Surface area of outside/otherside of object.
    is_basedon_outerarea : `bool`, default = False
        Specification whether overall heat transfer coefficient (U) 
        is based on outside/oneside surface area.
    is_basedon_innerarea : `bool`, default = False
        Specification whether overall heat transfer coefficient (U) 
        is based on inside/oneside surface area.
        

    Returns
    -------
    Overall heat transfer coefficient (U) : `int or float`
        If 'is_basedon_innerarea == True' 'U' based on inner/oneside area
        is returned
        
        If 'is_basedon_outerarea == True' 'U' based on inner/oneside area
        is returned
    
    
    Notes
    -----
    The following formula is used.
    
    .. math::
        
        \frac{1}{U_iA_i} = \frac{1}{U_oA_o} = 
        
        R_{convection, in} + R_{fouling, in} + R_{conduction}
        + R_{fouling, out} + R_{convection, out}
        
           
    *where:*
        
        :math:`U_i` *= overall heat transfer coefficient based on
        inner area*
        
        :math:`U_o` *= overall heat transfer coefficient based on
        outer area*
        
        :math:`A_i` *= inner surface area of heat transfer*
        
        :math:`A_o` *= outer area of heat transfer*
        
        :math:`R_{*}` *= resistance of conduction, convection or fouling*
    

    **What if some resistances are zero**
    
    *Simply leave their respective keywords at default value of 'None'.*
    
        
    See Also
    --------
    pychemengg.heattransfer.steadystate.Slab.resistanceof_cond
    pychemengg.heattransfer.steadystate.Slab.resistanceof_conv
    pychemengg.heattransfer.steadystate.Slab.resistanceof_fouling
    pychemengg.heattransfer.steadystate.Cylinder.resistanceof_cond
    pychemengg.heattransfer.steadystate.Cylinder.resistanceof_conv
    pychemengg.heattransfer.steadystate.Cylinder.resistanceof_fouling
    pychemengg.heattransfer.steadystate.Sphere.resistanceof_cond
    pychemengg.heattransfer.steadystate.Sphere.resistanceof_conv
    pychemengg.heattransfer.steadystate.Sphere.resistanceof_fouling
   
    
    Examples
    --------
    First import the module **heatcommonmethods**.
       
    >>> from pychemengg.heattransfer import heatcommonmethods as hcm 
    >>> hcm.calc_overallheattransfercoefficient(inner_heattransfercoefficient=800,
                              outer_heattransfercoefficient=1200,
                              inner_foulingfactor=0.0004,
                              outer_foulingfactor=0.0001,
                              resistanceof_wall=0.0025,
                              inner_area=0.0471,
                              outer_area=0.0597,
                              is_basedon_outerarea=True)
    315.06138931035497

    
    References
    ----------
    [1] G. F. Nellis and S. A. Klein, "Introduction to Engineering 
    Heat Transfer", 1st Edition. Cambridge University Press, 2021.
    
    [2] Y. A. Cengel and A. J. Ghajar, "Heat And Mass Transfer
    Fundamentals and Applications", 6th Edition. New York, McGraw Hill
    Education, 2020.
    
    [3] T.L. Bergman, A. S. Lavine, F. P. Incropera, D. P. Dewitt, 
    "Fundamentals of Heat and Mass Transfer", 7th Edition, 
    John Wiley, 2011.
    """
    resistanceof_conv_inside = 1/(inner_heattransfercoefficient*inner_area) if inner_heattransfercoefficient is not None else 0
    resistanceof_fouling_inside = inner_foulingfactor/inner_area if inner_foulingfactor is not None else 0
    resistanceof_wall = resistanceof_wall if resistanceof_wall is not None else 0
    resistanceof_conv_outside = 1/(outer_heattransfercoefficient*outer_area) if outer_heattransfercoefficient is not None else 0
    resistanceof_fouling_outside = outer_foulingfactor/outer_area if outer_foulingfactor is not None else 0
    
    U_times_area_inverse = (resistanceof_conv_inside
                  + resistanceof_fouling_inside
                  + resistanceof_wall
                  + resistanceof_conv_outside
                  + resistanceof_fouling_outside)
    if (is_basedon_outerarea==True) and (is_basedon_innerarea==False):
        return ((U_times_area_inverse)**(-1)) / outer_area
    if (is_basedon_outerarea==False) and (is_basedon_innerarea==True):
        return ((U_times_area_inverse)**(-1)) / inner_area
