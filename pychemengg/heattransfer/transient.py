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
"""
Module for transient/unsteady state heat transfer.

"""
import numpy as np
from scipy.optimize import brentq, fsolve
from scipy.special import j0, j1, erfc

__all__ = ["LumpedSystem", "NonLumpedSlab", "NonLumpedCylinder",
           "NonLumpedSphere", "SemiInfinite"]


sin=np.sin
cos=np.cos


class LumpedSystem():
    r""" Model for lumped system analysis.


    Parameters
    ----------
    surfacearea : `int or float`
        Surface area of solid object.
    volume : `int or float`
        Volume of solid object.
    density : `int or float`
        Density of solid object.
    specificheat : `int or float`
        Specific heat of solid object.
    thermalconductivity : `int or float`
        Thermal conductivity of solid object.
    heattransfercoefficient : `int or float`
        Heat transfer coefficient between solid object and surrounding.
    T_infinity : `int or float`
        Temperature of surroundings.
    T_initial : `int or float`
        Temperature of solid object at time = 0.
    
    
    Attributes
    ----------
    See "Parameters". All parameters are attributes. Additional attributes are listed below.
    mass : `int or float`
        Mass of solid object computed as (volume * density) of solid object.
    characteristiclength : `int or float`
        Characteristic length of object computed as (volume/surface area) of object.
     
    

    Examples
    --------
    First import the module **transient**
    
    Units used in this example: SI system
    
    However, any consistent units can be used
    
    >>> from pychemengg.heattransfer import transient
    >>> plate = transient.LumpedSystem(thermalconductivity=180,
                                       density=2800,
                                       specificheat=880,
                                       T_initial=700,
                                       T_infinity=15,
                                       heattransfercoefficient=53,
                                       surfacearea=2*1,
                                       volume=1*2e-2)
    # This will create an instance of 'LumpedSystem' with a name 'plate' 

    """
    
    def __init__(self, surfacearea=None,
                 volume=None,
                 density=None,
                 specificheat=None,
                 thermalconductivity=None,
                 heattransfercoefficient=None,
                 T_infinity=None,
                 T_initial=None):
        # assign
        self.surfacearea = surfacearea
        self.volume = volume
        self.density=density
        self.specificheat=specificheat
        self.thermalconductivity = thermalconductivity
        self.heattransfercoefficient=heattransfercoefficient
        self.T_infinity=T_infinity
        self.T_initial=T_initial
        # calculate
        self.mass = self.volume * self.density
        self.characteristiclength = self.volume/self.surfacearea
    
    def calc_Bi(self):
        r"""Computes Biot number.
        
        
        Parameters
        ----------
        `None_required` : 'None'
            Attributes that are already defined are used in calculation.
                                
        
        Returns
        -------
        Bi : `int or float`
            Biot number
        
        
        Notes
        -----
        Biot number is calculated using the following formula.
            
        .. math::
            Bi = \frac {h L_{c}} {k}
        
        *where:*
        
            *h = heat transfer coefficient*
        
            *k = thermal conductivity of solid object*
            
            :math:`L_c` *= characteristic length of solid object*
            
            *Bi = Biot number*
                     
        
        Examples
        --------
        First import the module **transient**
        
        Units used in this example: SI system
        
        However, any consistent units can be used
        
        >>> from pychemengg.heattransfer import transient
        >>> plate = transient.LumpedSystem(thermalconductivity=180, density=2800, specificheat=880, T_initial=700, T_infinity=15, heattransfercoefficient=53, surfacearea=2*1, volume=1*2e-2)
        # This will create an instance of 'LumpedSystem' with a name 'plate'
        # Next call calc_Bi
        >>> plate.calc_Bi()
        0.0029444444444444444
        
        
        References
        ----------
        [1] G. F. Nellis and S. A. Klein, "Introduction to Engineering 
        Heat Transfer", 1st Edition. Cambridge University Press, 2021.
        
        [2] Y. A. Cengel and A. J. Ghajar, "Heat And Mass Transfer
        Fundamentals and Applications", 6th Edition. New York, McGraw Hill
        Education, 2020.      
        """        
        self.Bi = (self.heattransfercoefficient
                           * self.characteristiclength
                           / self.thermalconductivity)
        return self.Bi
    
    def temperature_of_solid_at_time_t(self, time=None):
        r"""Temperature of solid object at a given time = t.
        
        
        Parameters
        ----------
        time : 'int or float'
            Time instant from begining of process, at which temperature
            of solid object is to be found.
                                
        
        Returns
        -------
        temperature : `int or float`
            Temperature of solid object at time = t 
        
        
        Notes
        -----
        Temperature of solid object at time = t is calculated using the following formula:
            
        .. math::
            T(t) = T_{infinity} + (T_{initial} - T_{infinity}) e^{-bt}
        
        *where:*
        
            :math:`T_{infinity} = temperature of surrounding fluid`
        
            :math:`T_{initial} = intitial temperature of solid object`
            
            :math:`b = \frac{hA_s}{\rho V C_p}`
            
            t = time at which temperature is to be computed
            
                where:
                
                h =  heat transfer coefficient
                
                :math:`A_s = surface area of solid object`
                
                :math:`\rho` = density of solid object
                
                V = volume of solid object
               
                :math:`C_p` = specific heat of solid object
                     
        
        Examples
        --------
        First import the module **transient**
        
        Units used in this example: SI system
        
        However, any consistent units can be used
        
        >>> from pychemengg.heattransfer import transient
        >>> plate = transient.LumpedSystem(thermalconductivity=180, density=2800, specificheat=880, T_initial=700, T_infinity=15, heattransfercoefficient=53, surfacearea=2*1, volume=1*2e-2)
        # This will create an instance of 'LumpedSystem' with a name 'plate'
        # Let temperature at time = 60 s after start of the process be needed.
        # Next call the following
        >>> plate.temperature_of_solid_at_time_t(time=60s)
        617.0619799301729
        
        
        References
        ----------
        [1] G. F. Nellis and S. A. Klein, "Introduction to Engineering 
        Heat Transfer", 1st Edition. Cambridge University Press, 2021.
        
        [2] Y. A. Cengel and A. J. Ghajar, "Heat And Mass Transfer
        Fundamentals and Applications", 6th Edition. New York, McGraw Hill
        Education, 2020.   
        """
        b = ((self.heattransfercoefficient * self.surfacearea)
             / (self.density*self.volume*self.specificheat)) # unit of b is 1/time
        solidtemp_at_time_t = (self.T_infinity
                                + (self.T_initial-self.T_infinity)*np.exp(-b*time+0j))
        return solidtemp_at_time_t.real
    
    def heatrateof_conv_at_time_t(self, time=None):
        r"""Heat rate of convection between object and surroundings at a given time = t.
        
        
        Parameters
        ----------
        time : 'int or float'
            Time instant from begining of process, at which heat rate is to be found.
                    
        
        Returns
        -------
        heat rate of convection : `int or float; Positive: Heat is gained by object, Negative: Heat is lost by object`
            Heat rate of convection between solid object and surroundings at time = t.
        
        
        Notes
        -----
        Heat rate is calculated using the following formula:
            
        .. math::
            q_{t} = h A_s (T_{infinity} - T_{t})
            
        *where:*
        
            t = time at which temperature is to be computed    
        
            h = heat transfer coefficient    
        
            :math:`T_{infinity} = temperature of surrounding fluid`
        
            :math:`T_{t} = temperature of solid object at time = t`
            
            :math:`A_s = surface area of solid object`
            
            :math:`q_{t} = heat rate at time = t`
                     
        
        Examples
        --------
        First import the module **transient**
        
        Units used in this example: SI system
        
        However, any consistent units can be used
        
        >>> from pychemengg.heattransfer import transient
        >>> plate = transient.LumpedSystem(thermalconductivity=180, density=2800, specificheat=880, T_initial=700, T_infinity=15, heattransfercoefficient=53, surfacearea=2*1, volume=1*2e-2)
        # This will create an instance of 'LumpedSystem' with a name 'plate'
        # Let temperature at time = 60 s after start of the process be needed.
        # Next call the following
        >>> plate.heatrateof_conv_at_time_t(time=60)
        -63818.56987259833
        # negative value indicates heat is being lost by the solid object
        
        
        References
        ----------
        [1] G. F. Nellis and S. A. Klein, "Introduction to Engineering 
        Heat Transfer", 1st Edition. Cambridge University Press, 2021.
        
        [2] Y. A. Cengel and A. J. Ghajar, "Heat And Mass Transfer
        Fundamentals and Applications", 6th Edition. New York, McGraw Hill
        Education, 2020.  
        """
        q_rate = (self.heattransfercoefficient * self.surfacearea
                *(self.temperature_of_solid_at_time_t(time=time) - self.T_infinity))
        return q_rate

    def totalheat_transferred_during_interval_t(self, time=None):
        r"""Heat transferred between solid object and surroundings during 
        time interval = 0 to t.
        
        
        Parameters
        ----------
        time : 'int or float'
            Time-limit after start of process for which
            heat transferred is to be computed.
                    
        
        Returns
        -------
        total heat transferred : `int or float; Positive: Heat is gained by object, Negative: Heat is lost by object`
            Total heat transferred between object and
            surroundings during interval 0 to t
        
        
        Notes
        -----
        Total heat  transferred in interval 0 to t is calculated using the
        following formula:
            
        .. math::
            q_{0 to t} = m C_p (T_{t} - T_{inintial})
            
        *where:*
        
            t = time marking the interval [0, t] for which heat
            transferred is to be computed    
        
            m = mass of object
            
            :math:`C_{p} = specific heat of object`
        
            :math:`T_{t} = temperature of object at time = t`
        
            :math:`T_{initial} = temperature of object at time = 0`
            
            :math:`q_{0 to t} = heat transferred in interval [0, t]`
                     
        
        Examples
        --------
        First import the module **transient**
        
        Units used in this example: SI system
        
        However, any consistent units can be used
        
        >>> from pychemengg.heattransfer import transient
        >>> plate = transient.LumpedSystem(thermalconductivity=180, density=2800, specificheat=880, T_initial=700, T_infinity=15, heattransfercoefficient=53, surfacearea=2*1, volume=1*2e-2)
        # This will create an instance of 'LumpedSystem' with a name 'plate'
        # Let heat transferred in time = 60 s after start of the process be needed.
        # Next call the following
        >>> plate.totalheat_transferred_during_interval_t(time=60)
        -4087185.6290410785
        # negative value indicates heat is being lost by the solid object
        
        
        References
        ----------
        [1] G. F. Nellis and S. A. Klein, "Introduction to Engineering 
        Heat Transfer", 1st Edition. Cambridge University Press, 2021.
        
        [2] Y. A. Cengel and A. J. Ghajar, "Heat And Mass Transfer
        Fundamentals and Applications", 6th Edition. New York, McGraw Hill
        Education, 2020.   
        """
        qtotal = (self.mass * self.specificheat
                 * (self.temperature_of_solid_at_time_t(time=time) - self.T_initial))
        return qtotal
    
    def maxheattransferpossible(self):
        r"""Maximum possible heat transfer between solid object and surroundings.
        
        
        Parameters
        ----------
        `None_required` : 'None'
            This class takes no parameters for instance creation.
     
        
         Returns
         -------
         maximum heat transfer possible: `int or float; Positive: Heat is gained by object, Negative: Heat is lost by object`
             Maximum heat transfer posssible between object and
             surroundings
        
        
        Notes
        -----
        Maximum heat transfer possible between solid object and surroundings 
        is calculated using the following formula. This is based on the assumption
        that final object temperature will eventually reach surrounding temperature
        of :math:`T_{infinity}`
            
        .. math::
            q_{max} = m C_p (T_{infinity} - T_{initial})
            
        *where:*  
        
            m = mass of solid object
            
            :math:`C_{p} = specific heat of solid object`
        
            :math:`T_{infinity} = temperature of surrounding, which the solid object will eventually attain`
        
            :math:`T_{initial} = temperature of solid object at time = initial`
            
            :math:`q_{max} = max heat transfer possible`
                     
        
        Examples
        --------
        First import the module **transient**
        
        Units used in this example: SI system
        
        However, any consistent units can be used
        
        >>> from pychemengg.heattransfer import transient
        >>> plate = transient.LumpedSystem(thermalconductivity=180, density=2800, specificheat=880, T_initial=700, T_infinity=15, heattransfercoefficient=53, surfacearea=2*1, volume=1*2e-2)
        # This will create an instance of 'LumpedSystem' with a name 'plate'
        # Next call the following
        >>> plate.maxheattransferpossible()
        -33756800.0
        # negative value indicates heat is being lost by the solid object
        
        
        References
        ----------
        [1] G. F. Nellis and S. A. Klein, "Introduction to Engineering 
        Heat Transfer", 1st Edition. Cambridge University Press, 2021.
        
        [2] Y. A. Cengel and A. J. Ghajar, "Heat And Mass Transfer
        Fundamentals and Applications", 6th Edition. New York, McGraw Hill
        Education, 2020.       
        """
        qtotal_max = self.mass * self.specificheat * (self.T_infinity - self.T_initial)
        return qtotal_max



class NonLumpedSlab():    
    r""" Model for nonlumped analysis of rectangular solid object.


    Parameters
    ----------
    thickness : `int or float`
        Thickness of solid object
    surfacearea : `int or float`
        Surface area of solid object.
    volume : `int or float`
        Volume of solid object.
    density : `int or float`
        Density of solid object.
    specificheat : `int or float`
        Specific heat of solid object.
    thermalconductivity : `int or float`
        Thermal conductivity of solid object.
    thermaldiffusivity : `int or float`
        Thermal diffusivity of solid object.
    heattransfercoefficient : `int or float`
        Heat transfer coefficient between solid object and surrounding.
    T_infinity : `int or float`
        Temperature of surroundings.
    T_initial : `int or float`
        Temperature of solid object at time = 0.
    
    
    
    Attributes
    ----------
    See "Parameters". All parameters are attributes. Additional attributes are listed below.
    
    mass : `int or float`
        Mass of solid object computed as (volume * density) of solid object.
 
    

    Examples
    --------
    First import the module **transient**
    
    Units used in this example: SI system
    
    However, any consistent units can be used
    
    >>> from pychemengg.heattransfer import transient
    >>> plate = transient.NonLumpedSlab(thickness=4e-2, surfacearea=1, volume=1*4e-2, density=8530, specificheat=380, thermalconductivity=110, thermaldiffusivity=None, heattransfercoefficient=120, T_infinity=500, T_initial=20)
    # This will create an instance of 'NonLumpedSlab' with a name 'plate' 

    """
    def __init__(self, thickness=None,
                 surfacearea=None,
                 volume=None,
                 density=None,
                 specificheat=None,
                 thermalconductivity=None,
                 thermaldiffusivity=None,
                 heattransfercoefficient=None,
                 T_infinity=None,
                 T_initial=None):
        # assign
        self.thickness = thickness
        self.surfacearea=surfacearea
        self.volume=volume
        self.density=density
        self.specificheat=specificheat
        self.thermalconductivity = thermalconductivity
        self.heattransfercoefficient=heattransfercoefficient
        self.T_infinity=T_infinity
        self.T_initial=T_initial
        # calculate mass
        if self.density is not None:
            self.mass = self.volume * self.density
        # calculate thermal diffusivity
        if (self.density is not None) and (self.specificheat is not None):
            self.thermaldiffusivity = self.thermalconductivity/self.density/self.specificheat
        else:
            if thermaldiffusivity is not None:
                self.thermaldiffusivity = thermaldiffusivity

        
    def calc_Bi(self):
        r"""Computes Biot number.
        
        
        Parameters
        ----------
        `None_required` : 'None'
            Attributes that are already defined are used in calculation.
                                
        
        Returns
        -------
        Bi : `int or float`
            Biot number
        
        
        Notes
        -----
        Biot number is calculated using the following formula.
            
        .. math::
            Bi = \frac {h L_{c}} {k}
        
        *where:*
        
            *h = heat transfer coefficient*
        
            *k = thermal conductivity of solid object*
            
            :math:`L_c` *= characteristic length of solid object = thickness/2*
            
            *Bi = Biot number*
                     
        
        Examples
        --------
        First import the module **transient**
        
        Units used in this example: SI system
        
        However, any consistent units can be used
        
        >>> from pychemengg.heattransfer import transient
        >>> plate = transient.NonLumpedSlab(thickness=4e-2, surfacearea=1, volume=1*4e-2, density=8530, specificheat=380, thermalconductivity=110, thermaldiffusivity=None, heattransfercoefficient=120, T_infinity=500, T_initial=20)
        # This will create an instance of 'NonLumpedSlab' with a name 'plate' 
        # Next call calc_Bi
        >>> plate.calc_Bi()
        0.021818181818181816
        
        
        References
        ----------
        [1] G. F. Nellis and S. A. Klein, "Introduction to Engineering 
        Heat Transfer", 1st Edition. Cambridge University Press, 2021.
        
        [2] Y. A. Cengel and A. J. Ghajar, "Heat And Mass Transfer
        Fundamentals and Applications", 6th Edition. New York, McGraw Hill
        Education, 2020.       
        """

        
        self.Bi = (self.heattransfercoefficient
                           * self.thickness/2
                           / self.thermalconductivity)
        return self.Bi

    
    def calc_Fo(self, time=None):
        r"""Computes Fourier number.
        
        
        Parameters
        ----------
        time : 'int or float'
            Time at which temperature or heat transfer is to be evaluated.
                                
        
        Returns
        -------
        Fo : `int or float`
            Fourier number
        
        
        Notes
        -----
        Fourier number is calculated using the following formula.
            
        .. math::
            Fo = \frac {\alpha t} {L^2}
        
        *where:*
        
            :math:`\alpha` *= thermal diffusivity*
        
            *t = time at which temperature or heat transfer is to be evaluated*
            
            *L = (slab thickness)/2*
            
            *Fo = Fourier number*
                     
        
        Examples
        --------
        First import the module **transient**
        
        Units used in this example: SI system
        
        However, any consistent units can be used
        
        >>> from pychemengg.heattransfer import transient
        >>> plate = transient.NonLumpedSlab(thickness=4e-2, surfacearea=1, volume=1*4e-2, density=8530, specificheat=380, thermalconductivity=110, thermaldiffusivity=None, heattransfercoefficient=120, T_infinity=500, T_initial=20)
        # This will create an instance of 'NonLumpedSlab' with a name 'plate' 
        # Next call calc_Fo assuming temperature is required at 7 min
        >>> plate.calc_Fo(time=7*60)
        35.63275128031097
        
        
        References
        ----------
        [1] G. F. Nellis and S. A. Klein, "Introduction to Engineering 
        Heat Transfer", 1st Edition. Cambridge University Press, 2021.
        
        [2] Y. A. Cengel and A. J. Ghajar, "Heat And Mass Transfer
        Fundamentals and Applications", 6th Edition. New York, McGraw Hill
        Education, 2020.       
        """
        self.Fo = self.thermaldiffusivity*time/(self.thickness/2)**2
        return self.Fo
    
    
    def calc_eigenvalues(self, numberof_eigenvalues_desired=10):
        r"""Computes eigen values of characteristic equation for Slab geometry.
        
        
        Parameters
        ----------
        numberof_eigenvalues_desired : 'int or float' (default = 10)
            Number of eigen values desired for the characteristic equation.
                                
        
        Returns
        -------
        eigenvalues : `np.array of int or float`
            Eigen values
        
        
        Notes
        -----
        Eigen values are calculated as roots of the following equation.
            
        .. math::
            x_n tan(x_n) - Bi = 0 , n = 1 \hspace{2pt} to \hspace{2pt} \infty
            
        *where:*
        
            :math:`x_n` *= nth eigen value*        
            
            *Bi = Biot number*
                     
        
        Examples
        --------
        First import the module **transient**
        
        Units used in this example: SI system
        
        However, any consistent units can be used
        
        >>> from pychemengg.heattransfer import transient
        >>> plate = transient.NonLumpedSlab(thickness=4e-2, surfacearea=1, volume=1*4e-2, density=8530, specificheat=380, thermalconductivity=110, thermaldiffusivity=None, heattransfercoefficient=120, T_infinity=500, T_initial=20)
        # This will create an instance of 'NonLumpedSlab' with a name 'plate' 
        # Next call calc_Bi
        >>> plate.calc_Bi()
        0.021818181818181816
        # Let first 5 eigen values be required
        >>> plate.calc_eigenvalues(numberof_eigenvalues_desired=5)
        array([ 0.14717481,  3.1485222 ,  6.28665585,  9.42709237, 12.56810661])
        
        
        References
        ----------
        [1] G. F. Nellis and S. A. Klein, "Introduction to Engineering 
        Heat Transfer", 1st Edition. Cambridge University Press, 2021.
        
        [2] Y. A. Cengel and A. J. Ghajar, "Heat And Mass Transfer
        Fundamentals and Applications", 6th Edition. New York, McGraw Hill
        Education, 2020.       
        """
        slab_eigenfunction = lambda x, Bi: x*np.tan(x)-Bi
        slab_eigenvalues = _get_eigenvalues(slab_eigenfunction, Bi=self.Bi,
                                          numberof_eigenvalues_desired=numberof_eigenvalues_desired)
        self.eigenvalues = np.array(slab_eigenvalues)
        return self.eigenvalues

        
    def temperature_of_solid_at_time_t(self, time=None, positioninsolid=None):
        r"""Calculates temperature of solid object at a given time = t and position = x.
        
        
        Parameters
        ----------
        time : `int or float`
            Time instant from begining of process, at which temperature
            of solid object is to be found.
        positioninsolid : `int or float`
            Distance measured from center of rectangular object where temperature is to be found.
                                
        
        Returns
        -------
        temperature : `int or float`
            Temperature of solid object at time = t and position = x.
        
        
        Notes
        -----
        Temperature of solid object at time = t and position = x is calculated using the following formula:
            
        .. math::
            T(t) = T_{infinity} + (T_{initial} - T_{infinity}) \displaystyle\sum_{i=1}^\infty \cfrac{4sin(\lambda_i)}{2 \lambda_i + sin(2 \lambda_i)} e^{- \lambda_i^2 \tau} cos(\lambda_i x/L)
        
        *where:*
        
            :math:`T_{infinity} = temperature of surrounding fluid`
        
            :math:`T_{initial} = intitial temperature of solid object`
            
            :math:`\lambda_n = nth eigen value of x_n tan(x_n) - Bi = 0 , n = 1 \hspace{2pt} to \hspace{2pt} \infty`
            
            :math:`\tau = Fourier number`
            
            x = distance from center of solid slab where temperature  is required
            
            L = thickness/2
                     
        
        Examples
        --------
        First import the module **transient**
        
        Units used in this example: SI system
        
        However, any consistent units can be used
        
        >>> from pychemengg.heattransfer import transient
        >>> plate = transient.NonLumpedSlab(thickness=4e-2, surfacearea=1, volume=1*4e-2, density=8530, specificheat=380, thermalconductivity=110, thermaldiffusivity=None, heattransfercoefficient=120, T_infinity=500, T_initial=20)
        # This will create an instance of 'NonLumpedSlab' with a name 'plate'
        # Let temperature at time = 7 min after start of the process be required.
        # Next call the following
        >>> plate.calc_Bi()
        0.021818181818181816
        >>> plate.calc_Fo(time=7*60)
        35.63275128031097
        plate.calc_eigenvalues()
        array([ 0.14717481,  3.1485222 ,  6.28665585,  9.42709237, 12.56810661,
               15.70935213, 18.85071334, 21.99214066, 25.13360932, 28.27510552])
        >>> plate.temperature_of_solid_at_time_t(time=7*60, positioninsolid=plate.thickness/2)
        279.76430920417204
        
        
        References
        ----------
        [1] G. F. Nellis and S. A. Klein, "Introduction to Engineering 
        Heat Transfer", 1st Edition. Cambridge University Press, 2021.
        
        [2] Y. A. Cengel and A. J. Ghajar, "Heat And Mass Transfer
        Fundamentals and Applications", 6th Edition. New York, McGraw Hill
        Education, 2020.       
        """
        term1 = 4*np.sin(self.eigenvalues)
        term2 = 2*self.eigenvalues + np.sin(2*self.eigenvalues)
        term3 = np.exp(-np.power(self.eigenvalues,2) * self.Fo)
        term4 = np.cos(self.eigenvalues*positioninsolid/(self.thickness/2))
        theta = np.sum(term1/term2*term3*term4)
        self.solidtemp_at_time_t = (self.T_infinity
                                + (self.T_initial-self.T_infinity)*theta)
        return self.solidtemp_at_time_t
    

    def heatrateof_conv_at_time_t(self, time=None):
        r"""Heat rate of convection between object and surroundings at a given time = t.
        
        
        Parameters
        ----------
        time : `int or float`
            Time instant from begining of process, at which heat rate is to be found.
                    
        
        Returns
        -------
        heat rate of convection : `int or float ; Positive: Heat is gained by object, Negative: Heat is lost by object`
            Heat rate of convection between solid object and surroundings at time = t.
        
        
        Notes
        -----
        Heat rate of convection is calculated using the following formula:
            
        .. math::
            q_{t} = h A_s (T_{infinity} - T_{t})
            
        *where:*
        
            t = time at which temperature is to be computed    
        
            h = heat transfer coefficient    
        
            :math:`T_{infinity} = temperature of surrounding fluid`
        
            :math:`T_{t} = temperature of surface of solid object at time = t`
            
            :math:`A_s` *= surface area of solid object*
            
            :math:`q_{t}` *= heat rate of convection at time = t*
                     
        
        Examples
        --------
        First import the module **transient**
        
        Units used in this example: SI system
        
        However, any consistent units can be used
        
        >>> from pychemengg.heattransfer import transient
        >>> plate = transient.NonLumpedSlab(thickness=4e-2, surfacearea=1, volume=1*4e-2, density=8530, specificheat=380, thermalconductivity=110, thermaldiffusivity=None, heattransfercoefficient=120, T_infinity=500, T_initial=20)
        # This will create an instance of 'NonLumpedSlab' with a name 'plate'
        # Let temperature at time = 7 min after start of the process be required.
        # Next call the following
        >>> plate.calc_Bi()
        0.021818181818181816
        >>> plate.calc_Fo(time=7*60)
        35.63275128031097
        plate.calc_eigenvalues()
        array([ 0.14717481,  3.1485222 ,  6.28665585,  9.42709237, 12.56810661,
               15.70935213, 18.85071334, 21.99214066, 25.13360932, 28.27510552])
        >>> plate.temperature_of_solid_at_time_t(time=7*60, positioninsolid=plate.thickness/2)
        279.76430920417204
        # Next call the following
        >>> plate.heatrateof_conv_at_time_t(time=7*60)
        26428.282895499357
        # Positive sign indicates gain of heat by the solid object
        
        
        References
        ----------
        [1] G. F. Nellis and S. A. Klein, "Introduction to Engineering 
        Heat Transfer", 1st Edition. Cambridge University Press, 2021.
        
        [2] Y. A. Cengel and A. J. Ghajar, "Heat And Mass Transfer
        Fundamentals and Applications", 6th Edition. New York, McGraw Hill
        Education, 2020.       
        """
        qrate = (self.heattransfercoefficient * self.surfacearea
                *(self.T_infinity - self.temperature_of_solid_at_time_t(time=None, positioninsolid=self.thickness/2)))
        # For convection, surface temperature is required, therefore let
        # positioninsolid = surface position = thickness/2 because origin is in middle
        return qrate


    def totalheat_transferred_during_interval_t(self, time=None):
        r"""Heat transferred between solid object and surroundings during 
        time interval = 0 to t.
        
        
        Parameters
        ----------
        time : 'int or float'
            Time-limit after start of process for which
            heat transferred is to be computed.
                    
        
        Returns
        -------
        total heat transferred : `int or float; Positive: Heat is gained by object, Negative: Heat is lost by object`
            Total heat transferred between object and
            surroundings during interval 0 to t
        
        
        Notes
        -----
        Total heat  transferred in interval 0 to t is calculated using the
        following formula:
            
        .. math::
            q_{0 to t} = q_{max} (1 - \displaystyle\sum_{i=1}^\infty \cfrac{4Sin( \lambda_i)}{2 \lambda_i + Sin(2 \lambda_i)} \frac{Sin( \lambda_i)}{\lambda_i} e^{- \lambda_n^2 \tau}
            
        *where:*
        
            :math:`\lambda_i = i^{th} eigen value of x_i tan(x_i) - Bi = 0 , i = 1 \hspace{2pt} to \hspace{2pt} \infty`
            
            :math:`\tau = Fourier number`
                        
            :math:`q_{max} = maximum heat transfer possible between object and surroundings`
       

        
        See Also
        ----------
        pychemengg.heattransfer.transient.NonLumpedSlab.maxheattransferpossible
                
                     
        
        Examples
        --------
        First import the module **transient**
        
        Units used in this example: SI system
        
        However, any consistent units can be used
        
        >>> from pychemengg.heattransfer import transient
        >>> plate = transient.NonLumpedSlab(thickness=4e-2, surfacearea=1, volume=1*4e-2, density=8530, specificheat=380, thermalconductivity=110, thermaldiffusivity=None, heattransfercoefficient=120, T_infinity=500, T_initial=20)
        # This will create an instance of 'NonLumpedSlab' with a name 'plate'
        # Let temperature at time = 7 min after start of the process be required.
        # Next call the following
        >>> plate.calc_Bi()
        0.021818181818181816
        >>> plate.calc_Fo(time=7*60)
        35.63275128031097
        plate.calc_eigenvalues()
        array([ 0.14717481,  3.1485222 ,  6.28665585,  9.42709237, 12.56810661,
               15.70935213, 18.85071334, 21.99214066, 25.13360932, 28.27510552])
        >>> plate.totalheat_transferred_during_interval_t(time=7*60)
        33472028.92491645
        # Positive value means heat is gained by the object.
        
        
        References
        ----------
        [1] G. F. Nellis and S. A. Klein, "Introduction to Engineering 
        Heat Transfer", 1st Edition. Cambridge University Press, 2021.
        
        [2] Y. A. Cengel and A. J. Ghajar, "Heat And Mass Transfer
        Fundamentals and Applications", 6th Edition. New York, McGraw Hill
        Education, 2020.       
        """
        term1 = 4*np.sin(self.eigenvalues)
        term2 = 2*self.eigenvalues + np.sin(2*self.eigenvalues)
        term3 = np.exp(-np.power(self.eigenvalues,2) * self.Fo)        
        term4 = np.sin(self.eigenvalues)/self.eigenvalues
        normalized_heatamount = 1 - np.sum(term1/term2*term3*term4)
        heattransferred = self.maxheattransferpossible() * normalized_heatamount
        return heattransferred  
    
    def maxheattransferpossible(self):
        r"""Maximum possible heat transfer between solid object and surroundings.
        
        
        Parameters
        ----------
        `None_required` : 'None'
            Attributes that are already defined are used in calculation.
     
        
         Returns
         -------
         maximum heat transfer possible: `int or float; Positive: Heat is gained by object, Negative: Heat is lost by object`
             Maximum heat transfer posssible between object and surroundings.
        
        
        Notes
        -----
        Maximum heat transfer possible between solid object and surroundings 
        is calculated using the following formula. This is based on the assumption
        that final object temperature will eventually reach surrounding temperature
        of :math:`T_{infinity}`
            
        .. math::
            q_{max} = m C_p (T_{infinity} - T_{initial})
            
        *where:*  
        
            m = mass of solid object
            
            :math:`C_{p} = specific heat of solid object`
        
            :math:`T_{infinity} = temperature of surrounding, which the solid object will eventually attain`
        
            :math:`T_{initial} = temperature of solid object at time = initial`
            
            :math:`q_{max} = max heat transfer possible`
                     
        
        Examples
        --------
        First import the module **transient**
        
        Units used in this example: SI system
        
        However, any consistent units can be used
        
        >>> from pychemengg.heattransfer import transient
        >>> plate = transient.NonLumpedSlab(thickness=4e-2, surfacearea=1, volume=1*4e-2, density=8530, specificheat=380, thermalconductivity=110, thermaldiffusivity=None, heattransfercoefficient=120, T_infinity=500, T_initial=20)
        # This will create an instance of 'NonLumpedSlab' with a name 'plate'
        >>> plate.maxheattransferpossible()
        62234880.0
        # negative value indicates heat is being lost by the solid object
        
        
        References
        ----------
        [1] G. F. Nellis and S. A. Klein, "Introduction to Engineering 
        Heat Transfer", 1st Edition. Cambridge University Press, 2021.
        
        [2] Y. A. Cengel and A. J. Ghajar, "Heat And Mass Transfer
        Fundamentals and Applications", 6th Edition. New York, McGraw Hill
        Education, 2020.       
        """
        qtotal_max = self.mass * self.specificheat * (self.T_infinity - self.T_initial)
        return qtotal_max
   
    
class NonLumpedCylinder():
    def __init__(self, radius=None,
                 surfacearea=None,
                 volume=None,
                 density=None,
                 specificheat=None,
                 thermalconductivity=None,
                 thermaldiffusivity=None,
                 heattransfercoefficient=None,
                 T_infinity=None,
                 T_initial=None):
        # assign
        self.radius = radius
        self.surfacearea=surfacearea
        self.volume=volume
        self.density=density
        self.specificheat=specificheat
        self.thermalconductivity = thermalconductivity
        self.heattransfercoefficient=heattransfercoefficient
        self.T_infinity=T_infinity
        self.T_initial=T_initial
        # calculate mass
        if self.density is not None:
            self.mass = self.volume * self.density
        # calculate thermal diffusivity
        if (self.density is not None) and (self.specificheat is not None):
            self.thermaldiffusivity = self.thermalconductivity/self.density/self.specificheat
        else:
            if thermaldiffusivity is not None:
                self.thermaldiffusivity = thermaldiffusivity

        
    def calc_Bi(self):        
        self.Bi = (self.heattransfercoefficient
                           * self.radius
                           / self.thermalconductivity)
        return self.Bi

    
    def calc_Fo(self, time=None):
        self.Fo = self.thermaldiffusivity * time / (self.radius)**2
        return self.Fo
    
    
    def calc_eigenvalues(self, numberof_eigenvalues_desired=None):
        cylinder_eigenfunction = lambda x, Bi: x*j1(x)/j0(x)-Bi
        cylinder_eigenvalues = _get_eigenvalues(cylinder_eigenfunction, Bi=self.Bi,
                                          numberof_eigenvalues_desired=numberof_eigenvalues_desired)
        self.eigenvalues = np.array(cylinder_eigenvalues)
        return self.eigenvalues

        
    def temperature_of_solid_at_time_t(self, time=None, rposition_tofindtemp=None):
        term1 = 2/self.eigenvalues*j1(self.eigenvalues)
        term2 = np.power(j0(self.eigenvalues),2) + np.power(j1(self.eigenvalues),2)
        term3 = np.exp(-np.power(self.eigenvalues,2) * self.Fo)        
        term4 = j0(self.eigenvalues*rposition_tofindtemp/self.radius)
        theta = np.sum(term1/term2*term3*term4)
        solidtemp_at_time_t = (self.T_infinity
                                + (self.T_initial-self.T_infinity)*theta)
        self.solidtemp_at_time_t  = solidtemp_at_time_t 
        return solidtemp_at_time_t

    
    def totalheat_transferred_during_interval_t(self):
        term1 = 2/self.eigenvalues*j1(self.eigenvalues)
        term2 = np.power(j0(self.eigenvalues),2) + np.power(j1(self.eigenvalues),2)
        term3 = np.exp(-np.power(self.eigenvalues,2) * self.Fo)        
        term4 = 2*j1(self.eigenvalues)/self.eigenvalues
        normalized_heatamount = 1 - np.sum(term1/term2*term3*term4)
        heattransferred = self.maxheattransferpossible() * normalized_heatamount
        return heattransferred  
    

    def heatrateof_conv_at_time_t(self):
        qrate = (self.heattransfercoefficient * self.surfacearea
                *(self.solidtemp_at_time_t - self.T_infinity))
        return qrate

    
    def maxheattransferpossible(self):
        qtotal_max = self.mass * self.specificheat * (self.T_infinity - self.T_initial)
        return qtotal_max


class NonLumpedSphere():
    def __init__(self, radius=None,
                 surfacearea=None,
                 volume=None,
                 density=None,
                 specificheat=None,
                 thermalconductivity=None,
                 thermaldiffusivity=None,
                 heattransfercoefficient=None,
                 T_infinity=None,
                 T_initial=None):
        # assign
        self.radius = radius
        self.surfacearea=surfacearea
        self.volume=volume
        self.density=density
        self.specificheat=specificheat
        self.thermalconductivity = thermalconductivity
        self.heattransfercoefficient=heattransfercoefficient
        self.T_infinity=T_infinity
        self.T_initial=T_initial
        # calculate
        if self.density is not None:
            self.mass = self.volume * self.density
        if (self.density is not None) and (self.specificheat is not None):
            self.thermaldiffusivity = self.thermalconductivity/self.density/self.specificheat
        else:
            self.thermaldiffusivity = thermaldiffusivity

        
    def calc_Bi(self):        
        self.Bi = (self.heattransfercoefficient
                           * self.radius
                           / self.thermalconductivity)
        return self.Bi

    
    def calc_Fo(self, time=None):
        self.Fo = self.thermaldiffusivity * time / (self.radius)**2
        return self.Fo
    
    
    def calc_eigenvalues(self, numberof_eigenvalues_desired=None):
        sphere_eigenfunction = lambda x,Bi: 1-x/np.tan(x)-Bi
        sphere_eigenvalues = _get_eigenvalues(sphere_eigenfunction, Bi=self.Bi,
                                          numberof_eigenvalues_desired=numberof_eigenvalues_desired)
        self.eigenvalues = np.array(sphere_eigenvalues)
        return self.eigenvalues

        
    def temperature_of_solid_at_time_t(self, time=None, rposition_tofindtemp=None):
        term1 = 4*(np.sin(self.eigenvalues)-self.eigenvalues*np.cos(self.eigenvalues))
        term2 = 2*self.eigenvalues - np.sin(2*self.eigenvalues)
        term3 = np.exp(-np.power(self.eigenvalues,2) * self.Fo)
        if rposition_tofindtemp == 0:        
            term4 = 1
        else:
            term4 = (np.sin(self.eigenvalues*rposition_tofindtemp/self.radius)
                     /(self.eigenvalues*rposition_tofindtemp/self.radius))
            
        theta = np.sum(term1/term2*term3*term4)
        solidtemp_at_time_t = (self.T_infinity
                                + (self.T_initial-self.T_infinity)*theta)
        self.solidtemp_at_time_t  = solidtemp_at_time_t 
        return solidtemp_at_time_t

    
    def totalheat_transferred_during_interval_t(self):
        term1 = 4*(np.sin(self.eigenvalues)-self.eigenvalues*np.cos(self.eigenvalues))
        term2 = 2*self.eigenvalues - np.sin(2*self.eigenvalues)
        term3 = np.exp(-np.power(self.eigenvalues,2) * self.Fo)        
        term4 = 3*((np.sin(self.eigenvalues-self.eigenvalues*np.cos(self.eigenvalues)))
                   / np.power(self.eigenvalues,3))
        normalized_heatamount = 1 - np.sum(term1/term2*term3*term4)
        heattransferred = self.maxheattransferpossible() * normalized_heatamount
        return heattransferred  
    

    def heatrateof_conv_at_time_t(self):
        qrate = (self.heattransfercoefficient * self.surfacearea
                *(self.solidtemp_at_time_t - self.T_infinity))
        return qrate

    
    def maxheattransferpossible(self):
        qtotal_max = self.mass * self.specificheat * (self.T_infinity - self.T_initial)
        return qtotal_max
         
                
def _get_eigenvalues(func, Bi=5, numberof_eigenvalues_desired=1):
    sol = []
    increment_leftrange = 1e-5
    incrementvalue = 0.1
    left=0.0  + increment_leftrange
    right=left
    for root_count in range(numberof_eigenvalues_desired):
        roots_found = False
        while roots_found == False:
            try:
                a = brentq(func, left, right, Bi)
                if abs(func(a, Bi)) <= 1e-7:
                    sol.append(a)
                    roots_found = True
                else:
                    left=right
                    right=left
            except Exception as e:
                if str(e) == "f(a) and f(b) must have different signs":
                    right = right + incrementvalue
        left = right
        right = left
    return sol


class SemiInfinite():
    
    def __init__(self, option=None,
                 distance_fromsurface=None,
                 time=None,
                 thermaldiffusivity=None,
                 thermalconductivity=None,
                 specificheat=None,
                 density=None,
                 T_initial=None,
                 T_infinity=None,
                 T_surface=None,
                 surfaceheatflux=None,
                 surfaceheattransfercoefficient=None,
                 surfaceenergypulse=None):
        self.option = option
        self.distance_fromsurface = distance_fromsurface
        self.time = time
        self.thermalconductivity = thermalconductivity
        self.specificheat=specificheat
        self.density=density
        self.T_initial = T_initial
        self.T_infinity = T_infinity
        self.T_surface = T_surface
        self.surfaceheatflux = surfaceheatflux
        self.surfaceheattransfercoefficient = surfaceheattransfercoefficient
        self.surfaceenergypulse = surfaceenergypulse
        # calculate
        if (self.density is not None) and (self.specificheat is not None):
            self.thermaldiffusivity = self.thermalconductivity/self.density/self.specificheat
        else:
            self.thermaldiffusivity = thermaldiffusivity
  
    def temp_at_given_distance_and_time(self):
        if self.option == "specified_surfacetemperature":
            theta = erfc(self.distance_fromsurface/2/np.power(self.thermaldiffusivity*self.time, 0.5))
            temp_at_given_x_and_time = self.T_initial + theta * (self.T_surface-self.T_initial)
        if self.option == "specified_surfaceheatflux":
            term1 = np.power(4*self.thermaldiffusivity*self.time/np.pi,0.5)
            term2 = np.exp(np.power(-self.distance_fromsurface,2)/4/self.thermaldiffusivity/self.time)
            term3 = self.distance_fromsurface*erfc(self.distance_fromsurface/2/np.power(self.thermaldiffusivity*self.time,0.5))
            temp_at_given_x_and_time = self.T_initial + self.surfaceheatflux/self.thermalconductivity*(term1*term2-term3)
        if self.option == "specified_surfaceconvection":
            term0 = self.distance_fromsurface/2/np.power(self.thermaldiffusivity*self.time, 0.5)
            term1 = erfc(term0)
            term2 = self.surfaceheattransfercoefficient*self.distance_fromsurface/self.thermalconductivity
            term3 = np.power(self.surfaceheattransfercoefficient, 2)*self.thermaldiffusivity*self.time/np.power(self.thermalconductivity,2)
            term4 = self.surfaceheattransfercoefficient*np.power(self.thermaldiffusivity*self.time, 0.5)/self.thermalconductivity
            theta = term1 - np.exp(term2+term3)*erfc(term0+term4)
            temp_at_given_x_and_time = self.T_initial + theta * (self.T_infinity-self.T_initial)
        if self.option == "specified_surfaceenergypulse":
            term1 = self.surfaceenergypulse/self.thermalconductivity
            term2 = np.power(np.pi*self.time/self.thermaldiffusivity, 0.5)
            term3 = np.exp(- np.power(self.distance_fromsurface,2)/4/self.thermaldiffusivity/self.time)
            temp_at_given_x_and_time = self.T_initial + term1/term2*term3
        self.temp_at_given_x_and_time = temp_at_given_x_and_time
        return self.temp_at_given_x_and_time
    
    def heatflux_at_surface (self):
        if self.option == "specified_surfacetemperature":
            term1 = np.power(np.pi*self.thermaldiffusivity*self.time, 0.5)
            heatflux = self.thermalconductivity * (self.T_surface - self.T_initial)/term1
            return heatflux
            
    def contacting_temp(self, other):
        self_param = np.power(self.thermalconductivity*self.density*self.specificheat, 0.5)
        other_param = np.power(other.thermalconductivity*other.density*other.specificheat, 0.5)
        term1 = self_param*self.T_initial + other_param*other.T_initial
        term2 = self_param + other_param
        self.contact_temp = other.contact_temp = term1/term2
        return self.contact_temp
        


if __name__ == '__main__':
    
    
# LUMPED SYSTEM EXAMPLES
#------------------------
# Ghajjar page 255
# lumped system

    heattransfercoefficient_fxn = lambda vel: 33*vel**0.8
    def findvel(vel):
        time = 10/.04
        plate = LumpedSystem(thermalconductivity=180,
                             density=2800,
                             specificheat=880,
                             T_initial=700,
                             T_infinity=15,
                             heattransfercoefficient=heattransfercoefficient_fxn(vel[0]),
                             surfacearea=2*1,
                             volume=1*2e-2)
        return (plate.temperature_of_solid_at_time_t(time))-50
    guess_value_vel = [1]
    velocity_needed = fsolve(findvel, guess_value_vel)
    heattransfercoefficient_needed = heattransfercoefficient_fxn(velocity_needed[0])
    print(f"velocity required = {velocity_needed[0]: 0.1f} m/s")
    print(f"heat_coeff needed = {heattransfercoefficient_needed: 0.0f} W/m2K")
    # check
    plate = LumpedSystem(thermalconductivity=180,
                         density=2800,
                         specificheat=880,
                         T_initial=700,
                         T_infinity=15,
                         heattransfercoefficient=heattransfercoefficient_needed,
                         surfacearea=2*1,
                         volume=1*2e-2)
    print(f"calculated characteristic length = {plate.characteristiclength}" \
          " and book value is = 0.01 m")
    print(f"calculated biot number = {plate.Bi(): 0.4f}" \
          " and valuye of biot number from book = 0.0163")
    


# NON-LUMPED SYSTEM EXAMPLES
#-----------------------------    
    
    #Ghajjar 4-3 6th edition
    plate=NonLumpedSlab(thickness=4e-2, surfacearea=1,
                volume=1*4e-2, density=8530,
                specificheat=380, thermalconductivity=110,
                thermaldiffusivity=None,
                heattransfercoefficient=120, T_infinity=500,
                T_initial=20)
    plate.Bi()
    plate.Fo(time=7*60)
    plate.calc_eigenvalues(numberof_eigenvalues_desired=10)
    plate.temperature_of_solid_at_time_t(time=7*60, positioninsolid=plate.thickness/2)

    
    #Ghajjar 4-4 6th edition
    cylinder=NonLumpedCylinder(radius=10e-2, surfacearea=1,
                      T_initial=600, volume=np.pi*10e-2**2*1,
                      T_infinity=200, density=7900,
                      thermaldiffusivity=None,
                      specificheat=477, heattransfercoefficient=80,
                      thermalconductivity=14.9)
    cylinder.Bi()
    cylinder.Fo(time=45*60)
    cylinder.calc_eigenvalues(numberof_eigenvalues_desired=10)
    cylinder.temperature_of_solid_at_time_t(time=45*60, rposition_tofindtemp=0)
    cylinder.maxheattransferpossible()
    print("heat transferred_",cylinder.totalheat_transferred_during_interval_t())
    

# Ghajjar 4-67 Problem
# Cooling tomatos
# This Problem is not properly defined
# The values of T0=10 and ts=7.1 cannot simultaneoulsy be true
    # def problem4_67(heattransfercoefficient_guess):
    #     tomato=NonLumpedSphere(radius=.04,
    #                   surfacearea=4*np.pi*.04**2,
    #                   volume=4/3*np.pi*0.04**3,
    #                   density=999,
    #                   specificheat=3.99e3,
    #                   thermaldiffusivity=0.141e-6,
    #                   thermaldiffusivity=None,
    #                   T_initial=30,
    #                   T_infinity=7,
    #                   thermalconductivity=0.59,
    #                   heattransfercoefficient=heattransfercoefficient_guess[0])
    #     # print("mass", tomato.mass)
    #     # print("volume", tomato.volume)
    #     # print("density", tomato.density)
    #     # print("thermaldiffusivity =", tomato.thermaldiffusivity)
    #     print("biot =",tomato.Bi())
    #     print("fourier =", tomato.Fo(time=2*3600))
    #     print(tomato.calc_eigenvalues(numberof_eigenvalues_desired=1))
    #     print("eigen =", tomato.eigenvalues)
    #     temp_center = tomato.temperature_of_solid_at_time_t(time=2*3600, rposition_tofindtemp=0)
    #     temp_surface = tomato.temperature_of_solid_at_time_t(time=2*3600, rposition_tofindtemp=0.04)
    #     return abs(temp_center-10)+abs(temp_surface-7.1)
    # tomato = fsolve(problem4_67, 100)

    
    # Ghajjar 4-68 Problem
    # Cooling tomatos
    def problem4_68(time):
        potato=NonLumpedSphere(radius=.0275,
                      surfacearea=4*np.pi*.0275**2,
                      volume=4/3*np.pi*0.0275**3,
                      density=None,
                      specificheat=None,
                      thermaldiffusivity=0.14e-6,
                      T_initial=8,
                      T_infinity=97,
                      thermalconductivity=0.6,
                      heattransfercoefficient=1400)
        potato.Bi()
        potato.Fo(time=time[0])
        potato.calc_eigenvalues(numberof_eigenvalues_desired=1)
        potato.eigenvalues
        temp_center = potato.temperature_of_solid_at_time_t(time=time[0], rposition_tofindtemp=0)
        return temp_center-70
     
    potatocooktime = fsolve(problem4_68, 1000)
    print(f"time for center of potato to reach 70 C = {potatocooktime[0]: 0.0f} s and ans from book = 1068s")
       
                      

    #Bergman Page 305
    wall = NonLumpedSlab(thickness = 80e-3,
                surfacearea=1,
                volume=1*80e-3,
                density=7832,
                specificheat=434,
                thermalconductivity=63.9,
#               thermaldiffusivity=None, # this is already by default none,
# so it need not be declared
                heattransfercoefficient=500,
                T_infinity=60,
                T_initial=-20)
    wall.Bi()
    wall.Fo(time=8*60)
    wall.calc_eigenvalues(numberof_eigenvalues_desired=1)
    wall.temperature_of_solid_at_time_t(time=8*60, positioninsolid=0)
    wall.temperature_of_solid_at_time_t(time=8*60, positioninsolid=40e-3)
    # print(wall.heatrateof_conv_at_time_t())
    print(wall.totalheat_transferred_during_interval_t())
    
    
    
# # Bergman
# # Problem 5.38
#     def problem5_38(time):
#         wall = NonLumpedSlab(thickness = 0.1,
#                     surfacearea=1,
#                     volume=1*0.1,
#                     density=7830,
#                     specificheat=550,
#                     thermalconductivity=48,
#                     heattransfercoefficient=250,
#                     T_infinity=800,
#                     T_initial=200)
#         wall.Bi()
#         wall.Fo(time=time[0])
#         eigens = wall.calc_eigenvalues(numberof_eigenvalues_desired=1)
#         # print("eign", eigens)
#         T0 = wall.temperature_of_solid_at_time_t(time=time[0], positioninsolid=0)
#         return T0-550
    
#     time_solution = fsolve(problem5_38, 1)
#     print(time_solution[0], " s and from book = 861 s ")
    
# # Bergman Problem 5_40
    
    
#     def problem5_40(time):
#         wall = NonLumpedSlab(thickness=0.3,
#                 surfacearea=1,
#                 volume=1*0.3,
#                 density=2600,
#                 specificheat=1000,
#                 thermalconductivity=1.5,
#                 heattransfercoefficient=100,
#                 T_infinity=950,
#                 T_initial=20)
#         wall.Bi()
#         wall.Fo(time=time[0])
#         eigens = wall.calc_eigenvalues(numberof_eigenvalues_desired=5)
#         # print("eign", eigens)
#         T0 = wall.temperature_of_solid_at_time_t(time=time[0], positioninsolid=0)
#         return T0-750
    
#     time_solution = fsolve(problem5_40, 1)
#     print(time_solution[0], " s while from book = 33,800 s")


# # Problem 5_54 Bergman
#     def problem5_54(time):
#         bearings = NonLumpedCylinder(radius=0.1/2,
#                         surfacearea=np.pi*0.1*1,
#                         volume=np.pi/4*0.1*0.1*1,
#                         density=7900,
#                         thermalconductivity=19,
#                         specificheat=546,
#                         heattransfercoefficient=500,
#                         T_infinity=30,
#                         T_initial=500)
#         print("biot=",bearings.Bi())
#         fourier = bearings.Fo(time=time[0])
#         print("fourier = ", fourier)
#         print(bearings.thermaldiffusivity)
#         eigens = bearings.calc_eigenvalues(numberof_eigenvalues_desired=1)
#         print("eign", eigens)
#         T0 = bearings.temperature_of_solid_at_time_t(time=time[0],rposition_tofindtemp=0)
#         qrate = bearings.totalheat_transferred_during_interval_t()
#         print("qrate", qrate)
#         return T0-50
    
#     time_solution = fsolve(problem5_54, 1000)
#     print(time_solution[0], " s while from book = 33,800 s")
    
#Rathore
#Page 405 Example 6.22
    def example6_22(time):
        egg = NonLumpedSphere(radius=5e-2/2,
                     surfacearea=4*np.pi*(5e-2/2)**2,
                     volume=4/3*np.pi*(5e-2)**3,
                     density=None,
                     specificheat=None,
                     thermalconductivity=0.6,
                     thermaldiffusivity=0.14e-6,
                     heattransfercoefficient=1200,
                     T_initial=2,
                     T_infinity=100)
        egg.Bi()
        egg.Fo(time=time[0])
        egg.calc_eigenvalues(numberof_eigenvalues_desired=1)
        temp_center = egg.temperature_of_solid_at_time_t(time=time[0], rposition_tofindtemp=0)
        return temp_center-75
         
    eggcooktime = fsolve(example6_22, 1000)
    print(f"time for center of egg to reach 75 C = {eggcooktime[0]: 0.0f} s and ans from book = 962s ")
        
    


# Ghajjar 5t edition and 6th Edition
# Example 4-6
    def example4_6(burydistance):
        pipe = SemiInfinite(option="specified_surfacetemperature",
                            T_surface=-10,
                            T_initial=15,
                            thermalconductivity=0.4,
                            thermaldiffusivity=0.15e-6,
                            time = 90*24*3600,
                            distance_fromsurface=burydistance[0])
        return pipe.temp_at_given_distance_and_time()-0
    burydist = fsolve(example4_6, 1)
    print(f"burydistance = {burydist[0]: 0.2f} m and ans from book = 0.80 m")

# Ghajjar 6th edition
# Example 4-7

    wood = SemiInfinite(option="specified_surfaceheatflux",
                        time = 20*60,
                        T_initial=20,
                        surfaceheatflux=1250,
                        thermalconductivity=0.159,
                        thermaldiffusivity=1.75e-7,
                        distance_fromsurface=0)
    temp_wood = wood.temp_at_given_distance_and_time()
    print(f"temp of wood surface from code = {temp_wood: 0.0f} C and ans from book is 149 C")
    
    aluminum = SemiInfinite(option="specified_surfaceheatflux",
                        time = 20*60,
                        T_initial=20,
                        surfaceheatflux=1250,
                        thermalconductivity=237,
                        thermaldiffusivity=9.71e-5,
                        distance_fromsurface=0)
    temp_aluminum = aluminum.temp_at_given_distance_and_time()
    print(f"temp of aluminum surface from code = {temp_aluminum: 0.0f} C and ans from book is 22 C")
    
    # Ghajjar 6th edition
    # Example 4-7
    
    plate = SemiInfinite(option="specified_surfaceconvection",
                         T_infinity=-70,
                         T_initial=10,
                         thermalconductivity=16.3,
                         specificheat=500,
                         density=8000,
                         surfaceheattransfercoefficient=300,
                         distance_fromsurface=0.01,
                         time=30*60)
    tempatbolttip = plate.temp_at_given_distance_and_time()
    print(f"temp at bolt tip from code = {tempatbolttip: 0.1f} C and from book = -40.7 C")
    
# Ghajjar 5th Edn
# Problem 4_99
    def problem4_99(energypulse):
        slab = SemiInfinite(option="specified_surfaceenergypulse",
                            thermalconductivity=63.9,
                            thermaldiffusivity=18.8e-6,
                            time=30,
                            distance_fromsurface=25e-3,
                            T_initial=20,
                            surfaceenergypulse=energypulse[0])
        temp = slab.temp_at_given_distance_and_time()

        return temp-130
    energy = fsolve(problem4_99, 100)
    print(f"energy from code = {energy[0]: 0.3e} and from book = 2.076e7 ")
    
# Ghajjar 5th Edn
# Problem 4_102

    human = SemiInfinite(option="specified_surfacetemperature",
                        thermalconductivity=1,
                        density=1,
                        specificheat=1.1e3**2,
                        T_initial=32)
    
    aluminum = SemiInfinite(option="specified_surfacetemperature",
                        thermalconductivity=1,
                        density=1,
                        specificheat=24e3**2,
                        T_initial=20)
    temp_human_alumium = human.contacting_temp(aluminum)
    print(f"contacting temp of human and aluminum = {temp_human_alumium: 0.1f} C and from book = 20.5 C ")
    
    wood = SemiInfinite(option="specified_surfacetemperature",
                        thermalconductivity=1,
                        density=1,
                        specificheat=0.38e3**2,
                        T_initial=20)
    temp_human_alumium = human.contacting_temp(wood)
    print(f"contacting temp of human and aluminum = {temp_human_alumium: 0.1f} C and from book = 28.9 C ")
                        


                     



# # slab_eigenfunction = lambda x,Bi: x*np.tan(x)-Bi
# # slablambdas = _get_eigenvalues(slab_eigenfunction, Bi=0.9, numberof_eigenvalues_desired=10)

# sphere_eigenfunction = lambda x,Bi: 1-x/np.tan(x)-Bi
# spherelambdas = _get_eigenvalues(sphere_eigenfunction, Bi=64.2, numberof_eigenvalues_desired=10)

# # cylinder_eigenfunction = lambda x, Bi: x*j1(x)/j0(x)-Bi
# # cylinderlambdas = _get_eigenvalues(cylinder_eigenfunction, Bi=0.9, numberof_eigenvalues_desired=10)        
# # wall = NonLumpedSlab(thickness = 4e-2, area=1, thermalconductivity=110)
# # wall.Bi = 100
# # print(wall.lambdaval())
