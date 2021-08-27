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
Module for steady state heat transfer via conduction,
convection, and radiation.

"""
from scipy.optimize import fsolve
from scipy.constants import Stefan_Boltzmann
import math
from matplotlib import pyplot as plt
import numpy as np


__all__ = ["Slab", "Cylinder", "CylinderBase", "Sphere"]

class Slab(object):
    r""" Models a rectangular object


    Parameters
    ----------
    thickness : `int or float`
        Thickness of rectangular slab
    area : `int or float`
        Area of rectangular slab (if area is not known put area = 1.0)
    thermalconductivity : `int or float`
        Thermal conductivity of rectangular slab   
 
    
    Attributes
    ----------
    thickness : `int or float`
        Thickness of rectangular slab
    area : `int or float`
        Area of rectangular slab
    thermalconductivity : `int or float`
        Thermal conductivity of rectangular slab  
    

    Examples
    --------
    First import the module **steadystate**
    
    Units used in this example: SI system
    
    However, any consistent units can be used
    
    >>> from pychemengg.heattransfer import steadystate as ss 
    >>> innerWall = ss.Slab(thickness=0.13, thermalconductivity=0.7, area=1.0)
    # This will create an instance of 'Slab' with a name 'innerwall' 

    """
    def __init__(self, thickness=None, area=None, thermalconductivity=None):
        self.thickness = thickness
        self.area = area
        self.thermalconductivity = thermalconductivity
    
    def heatrateof_cond(self, dT=None):
        r"""Computes heat rate of conduction for a rectangular object
        
        
        Parameters
        ----------
        dT : `int or float`
            Temperature difference between two faces of a slab
                
        
        Returns
        -------
        heatrate : `int or float`
            rate of heat transfer by conduction
        
        
        Notes
        -----
        Heat rate of conduction is calculated using the Fourier's Law
            
        .. math::
            Q (heatrate) = k A \frac{\Delta T}{\Delta x}
        
        *where:*
        
            *k = thermal conductivity*
        
            *A = area of heat transfer*
             
            :math:`\Delta T` *= temperature difference*
             
            :math:`\Delta x` *= slab thickness*
        
        
        Examples
        --------
        First import the module **steadystate**
        
        Units used in this example: SI system
        
        However, any consistent units can be used
        
        >>> from pychemengg.heattransfer import steadystate as ss
        >>> wall = ss.Slab(thickness=0.13, thermalconductivity=0.7, area=1.0)
        >>> wall.heatrateof_cond(dT=100)
        538.4615384615385
        
        
        References
        ----------
        [1] Yunus A. Cengel and Afshin J. Ghajar, "Heat And Mass Transfer
        Fundamentals and Applications", 6th Edition. New York, McGraw Hill
        Education, 2020.       
        """
        return self.thermalconductivity * self.area * dT / self.thickness
    
    def heatrateof_conv(self, heattransfercoefficient=None, dT=None):
        r"""Computes heat rate of convection for a rectangular object

        
        Parameters
        ----------
        heattransfercoefficient : `int or float`
            Heat transfer coefficient *`h`* for the slab surface
        dT : `int or float`
            Temperature difference between slab surface and surrounding fluid
                
        
        Returns
        -------
        heatrate : `int or float`
            Rate of heat transfer by convection
        
        
        Notes
        -----
        Heat rate of convection is calculated using the Newton's Law
            
        .. math::
            Q (heatrate) = h A \Delta T
        
        *where:*
        
            *h = heat transfer coefficient*
        
            *A = area of heat transfer*
             
            :math:`\Delta T` *= temperature difference*
             
        
        Examples
        --------
        First import the module **steadystate**
        
        Units used in this example: SI system
        
        However, any consistent units can be used
        
        >>> from pychemengg.heattransfer import steadystate as ss
        >>> wall = ss.Slab(thickness=0.13, thermalconductivity=0.7, area=1.0)
        >>> wall.heatrateof_conv(heattransfercoefficient=134.0, dT=60)
        8040.0
        
        
        References
        ----------
        [1] Yunus A. Cengel and Afshin J. Ghajar, "Heat And Mass Transfer
        Fundamentals and Applications", 6th Edition. New York, McGraw Hill
        Education, 2020        
        """

        return heattransfercoefficient * self.area * dT
    
    def heatrateof_rad(self, T_infinity=None, T_surface=None, emissivity=None):
        r"""Computes heat rate of radiation for a rectangular object

        
        Parameters
        ----------
        T_infinity : `int or float`
            Temperature of surroundings in **absolute temperature units**
        T_surface : `int or float`
            Temperature of slab surface in **absolute temperature units**
        emissivity : `int or float`
            Emissivity of the slab
        
        
        Returns
        -------
        heatrate : `int or float` (returns a positive value)
            Rate of heat transfer by radiation
        
        
        Notes
        -----
        Heat rate of radiation is calculated using the Stefan-Boltzmann law
            
        .. math::
            Q (heatrate) = \sigma \epsilon A (T_{infinity}^4 - T_{surface}^4)
        
        *where:*
        
            :math:`\sigma` *= Stefan-Boltzmann constant*
            
            :math:`\epsilon` *= emissivity of object*
        
            *A = area of heat transfer*
             
            :math:`T_{infinity}` *= absolute temperature of surroundings*
            
            :math:`T_{surface}` *= absolute temperature of surface*
             
        
        Examples
        --------
        First import the module **steadystate**
        
        Units used in this example: SI system
        
        However, any consistent units can be used
        
        >>> from pychemengg.heattransfer import steadystate as ss
        >>> wall = ss.Slab(area=1.4)
        >>> wall.heatrateof_rad(T_infinity=10+273, T_surface=30+273, emissivity=0.95)
        151.93639338614008
        
        
        References
        ----------
        [1] Yunus A. Cengel and Afshin J. Ghajar, "Heat And Mass Transfer
        Fundamentals and Applications", 6th Edition. New York, McGraw Hill
        Education, 2020        
        """        
        return abs(Stefan_Boltzmann*emissivity*self.area*(T_infinity**4 - T_surface**4))
        
    def resistanceof_cond(self):
        r"""Computes resistance of conduction for a rectangular object

        
        Parameters
        ----------
        `None_required` : 'None'
            Uses attributes defined during instance creation

                
        Returns
        -------
        resistance : `int or float`
            Conduction resistance

               
        Notes
        -----
        The following formula is used:
            
        .. math::
            R_{conduction} = \frac{L}{kA}
        
        *where:*
        
            *L = thickness*
            
            *k = thermal conductivity*
            
            *A = surface area of heat transfer*
            
            :math:`R_{conduction}` *= conduction resistance*
            
        
        Examples
        --------
        First import the module **steadystate**
        
        Units used in this example: SI system
        
        However, any consistent units can be used
        
        >>> from pychemengg.heattransfer import steadystate as ss
        >>> wall = ss.Slab(thickness=0.3, thermalconductivity=0.9, area=15)
        >>> wall.resistanceof_cond()
        0.022222222222222223
        
        
        References
        ----------
        [1] Yunus A. Cengel and Afshin J. Ghajar, "Heat And Mass Transfer
        Fundamentals and Applications", 6th Edition. New York, McGraw Hill
        Education, 2020        
        """        
        return self.thickness / (self.thermalconductivity * self.area)

    def resistanceof_conv(self, heattransfercoefficient=None):
        r"""Computes resistance of convection for a rectangular object

        
        Parameters
        ----------
        heattransfercoefficient : `int or float`
            Heat transfer coefficient *`h`* for the slab surface

        
        Returns
        -------
        resistance : `int or float`
            Convection resistance
        
        
        Notes
        -----
        The following formula is used:         
        
        .. math::
            R_{convection} = \frac{1}{hA}
        
        *where:*
        
            *h = heat transfer coefficient*
                    
            *A = area of heat transfer*
            
            :math:`R_{convection}` *= convection resistance*
            
        
        Examples
        --------
        First import the module **steadystate**
        
        Units used in this example: SI system
        
        However, any consistent units can be used
        
        >>> from pychemengg.heattransfer import steadystate as ss
        >>> wall = ss.Slab(thickness=0.3, thermalconductivity=0.9, area=0.25*1)
        >>> wall.resistanceof_conv(heattransfercoefficient=10)
        0.4

        
        References
        ----------
        [1] Yunus A. Cengel and Afshin J. Ghajar, "Heat And Mass Transfer
        Fundamentals and Applications", 6th Edition. New York, McGraw Hill
        Education, 2020
        
        """
        return 1/(heattransfercoefficient * self.area)
    
    def resistanceof_fouling(self, foulingfactor=None):
        r"""Computes resistance of fouling for a rectangular object
        
        
        Parameters
        ----------
        foulingfactor : `int or float`
            Fouling factor :math:`R_f` for the slab surface
            
            typical units are :math:`m^2` K/W        
        
        Returns
        -------
        resistance : `int or float`
            Fouling resistance
        
        
        Notes
        -----
        The following formula is used:         
        
        .. math::
            R_{fouling} = \frac{R_f}{A}
        
        *where:*
        
            :math:`R_f` *= fouling factor*
                    
            *A = fouled area of heat transfer*
            
            :math:`R_{fouling}` *= fouling resistance*
            
        
        Examples
        --------
        First import the module **steadystate**
        
        Units used in this example: SI system
        
        However, any consistent units can be used
        
        >>> from pychemengg.heattransfer import steadystate as ss
        >>> wall = ss.Slab(thickness=0.3, thermalconductivity=0.9, area=0.25*1)
        >>> wall.resistanceof_fouling(foulingfactor=0.0007)
        0.0028
        
        
        References
        ----------
        [1] Yunus A. Cengel and Afshin J. Ghajar, "Heat And Mass Transfer
        Fundamentals and Applications", 6th Edition. New York, McGraw Hill
        Education, 2020
        
        """
        return foulingfactor / self.area


class Cylinder(object):
    r""" Models a cylindrical object
    
    
    Parameters
    ----------
    length : `int or float`
        Length of cylindrical object (if length is not known put length = 1.0)
    inner_radius: `int or float`
        Inner radius of cylindrical object
    outer_radius: `int or float`
        Outer radius of cylindrical object    
    thermalconductivity : `int or float`
        Thermal conductivity of cylindrical object   
    
    
    Attributes
    ----------
    length : `int or float`
        Length of cylindrical object (if length is not known put length = 1.0)
    inner_radius: `int or float`
        Inner radius of cylindrical object
    outer_radius: `int or float`
        Outer radius of cylindrical object    
    thermalconductivity : `int or float`
        Thermal conductivity of cylindrical object
    
    
    Examples
    --------
    First import the module **steadystate**
    
    Units used in this example: SI system
    
    However, any consistent units can be used
    
    >>> from pychemengg.heattransfer import steadystate as ss
    >>> insulation = Cylinder(length=5, inner_radius=1.5e-3, outer_radius=3.5e-3, thermalconductivity=0.15)
    # This will create an instance of 'Cylinder' with a name 'insulation'  
    """
    
    def __init__(self, length=None, inner_radius=None, outer_radius=None, thermalconductivity=None):
        self.length = length
        self.inner_radius = inner_radius
        self.outer_radius = outer_radius
        self.thermalconductivity = thermalconductivity
    
    def area(self, radius=None):
        r"""Computes surface area of a cylinder
        
        
        Parameters
        ----------
        radius : `int or float`
            Radius at which surface area is to be computed
                
        
        Returns
        -------
        area : `int or float`
            Surface area of cylinder
        
        
        Notes
        -----
        Surface area of cylinder is computed using:
            
        .. math::
            A = 2\pi rL
                    
        *where:*
        
            *r = radius*
        
            *L = length*
             
            *A = surface area*
        
        
        Examples
        --------
        First import the module **steadystate**
        
        Units used in this example: SI system
        
        However, any consistent units can be used
        
        >>> from pychemengg.heattransfer import steadystate as ss

        >>> insulation = Cylinder(length=5, inner_radius=1.5e-3, outer_radius=3.5e-3, thermalconductivity=0.15)
        # This will create an instance of 'Cylinder' with a name 'insulation'
        
        >>> insulation.area(radius=insulation.outer_radius)
        0.10995574287564275
        # This computes surface area @ outer radius
        
        >>> insulation.area(radius=insulation.inner_radius)
        0.047123889803846894
        # This computes surface area @ inner radius
        
        >>> insulation.area(radius=1.5e-3)
        0.047123889803846894
        # This will also compute surface area @ inner radius, which is = 1.5e-3.
        # Here the radius is entered as a number directly, rather than as an attribute.
        """
        return 2 * math.pi * radius * self.length
    
    def volume(self, radius=None):
        r"""Computes volume of a cylinder
        
        
        Parameters
        ----------
        radius : `int or float`
            Radius at which volume is to be computed
                
        
        Returns
        -------
        volume : `int or float`
            Volume of cylinder
        
        
        Notes
        -----
        Volume of cylinder is computed using:
            
        .. math::
            V = \pi r^2 L
                    
        *where:*
        
            *r = radius*
            
            *L = length*
             
            *V = volume*
        
        
        Examples
        --------
        First import the module **steadystate**
        
        Units used in this example: SI system
        
        However, any consistent units can be used
        
        >>> insulation = Cylinder(length=5, inner_radius=1.5e-3, outer_radius=3.5e-3, thermalconductivity=0.15)
        This will create an instance of 'Cylinder' with a name 'insulation'
        
        >>> insulation.volume(radius=insulation.outer_radius)
        0.00019242255003237485
        # This computes volume @ outer radius
        
        >>> insulation.volume(radius=insulation.inner_radius)
        3.534291735288517e-05
        # This computes volume @ inner radius
        
        >>>  insulation.volume(radius=1.5e-3)
        3.534291735288517e-05
        # This will also compute volume @ inner radius, which is = 1.5e-3.
        # Here the radius is entered as a number directly, rather than as an attribute.
        """
        return math.pi * radius**2 * self.length
    
    def heatrateof_cond(self, dT=None):
        r"""Computes heat rate of conduction for a cylindrical object
        
        
        Parameters
        ----------
        dT : `int or float`
            Temperature difference between two surfaces of the cylinder
                
        
        Returns
        -------
        heatrate : `int or float`
            Rate of heat transfer by conduction
        
        
        Notes
        -----
        The following formula is used:
            
        .. math::
            Q (heatrate) = \frac{\Delta T}{R_{conduction}}
        
        *where:*
             
            :math:`\Delta T` *= temperature difference*
            
            :math:`R_{conduction}` *= conduction resistance* given by
                         
                :math:`R_{conduction} = \cfrac{ln(r_o/r_i)}{2\pi kL}`
        
                :math:`r_o` *= outer radius of cylinder*
                
                :math:`r_i` *= inner radius of cylinder*
            
                *L = length*
                
                *k = thermal conductivity*
        
        
        Examples
        --------
        First import the module **steadystate**
        
        Units used in this example: SI system
        
        However, any consistent units can be used
        
        >>> from pychemengg.heattransfer import steadystate as ss
        >>> insulation = Cylinder(length=5, inner_radius=1.5e-3, outer_radius=3.5e-3, thermalconductivity=0.15)
        >>> insulation.heatrateof_cond(dT=75)
        417.12506315941755
        
        References
        ----------
        [1] Yunus A. Cengel and Afshin J. Ghajar, "Heat And Mass Transfer
        Fundamentals and Applications", 6th Edition. New York, McGraw Hill
        Education, 2020        
        """
        return dT/self.resistanceof_cond()
    
    def heatrateof_conv(self, heattransfercoefficient=None, radius=None, dT=None):
        r"""Computes heat rate of convection for a cylindrical object
        
        
        Parameters
        ----------
        heattransfercoefficient : `int or float`
            Heat transfer coefficient *`h`* for the cylinder surface
        radius : `int or float`
            Radius of cylinder where convective heat transfer rate is to be computed
        dT : `int or float`
            Temperature difference between cylinder surface and surrounding fluid
                
        
        Returns
        -------
        heatrate : `int or float`
            Rate of heat transfer by convection
        
        
        Notes
        -----
        Heat rate of convection is calculated using the Newton's Law
            
        .. math::
            Q (heatrate) = h A \Delta T
        
        *where:*
        
            *h = heat transfer coefficient*
        
            *A = area of heat transfer*
             
            :math:`\Delta T` *= temperature difference*
             
        
        Examples
        --------
        
        First import the module **steadystate**
        
        Units used in this example: SI system
        
        However, any consistent units can be used
        
        >>> from pychemengg.heattransfer import steadystate as ss
        >>> insulation = Cylinder(length=5, inner_radius=1.5e-3, outer_radius=3.5e-3, thermalconductivity=0.15)
        >>> insulation.heatrateof_conv(heattransfercoefficient=12, radius=insulation.outer_radius, dT=90.6-30)
        79.9598162191674
        # NOTE: Mathematical expressions (such as: 90.6 - 30) can be used as an argument
        
        
        References
        ----------
        [1] Yunus A. Cengel and Afshin J. Ghajar, "Heat And Mass Transfer
        Fundamentals and Applications", 6th Edition. New York, McGraw Hill
        Education, 2020        
        """
        return heattransfercoefficient * self.area(radius=radius) * dT
    
    def heatrateof_rad(self, radius=None, T_infinity=None, T_surface=None, emissivity=None):
        r"""Computes heat rate of radiation for a cylindrical object
        
        
        Parameters
        ----------
        radius : `int or float`
            Radius of cylinder where radiation heat transfer rate is to be computed
        T_infinity : `int or float`
            Temperature of surroundings in **absolute temperature units**
        T_surface : `int or float`
            Temperature of cylinder surface in **absolute temperature units**
        emissivity : `int or float`
            Emissivity of the cylinder
        
        
        Returns
        -------
        heatrate : `int or float` (returns a positive value)
            Rate of heat transfer by radiation
        
        
        Notes
        -----
        Heat rate of radiation is calculated using the Stefan-Boltzmann law
            
        .. math::
            Q (heatrate) = \sigma \epsilon A (T_{infinity}^4 - T_{surface}^4)
        
        *where:*
        
            :math:`\sigma` *= Stefan-Boltzmann constant*
            
            :math:`\epsilon` *= emissivity of object*
        
            *A = area of heat transfer*
             
            :math:`T_{infinity}^4` *= absolute temperature of surroundings*
            
            :math:`T_{surface}^4` *= absolute temperature of surface*
             
        
        Examples
        --------
        
        First import the module **steadystate**
        
        Units used in this example: SI system
        
        However, any consistent units can be used
        
        >>> from pychemengg.heattransfer import steadystate as ss
        >>> insulation = Cylinder(length=5, inner_radius=1.5e-3, outer_radius=3.5e-3, thermalconductivity=0.15)
        >>> insulation.heatrateof_rad(radius=insulation.outer_radius, T_infinity=30+273, T_surface=90.6+273, emissivity=0.95)
        53.60018341250309
        
        
        References
        ----------
        [1] Yunus A. Cengel and Afshin J. Ghajar, "Heat And Mass Transfer
        Fundamentals and Applications", 6th Edition. New York, McGraw Hill
        Education, 2020        
        """
        return abs(Stefan_Boltzmann * emissivity * self.area(radius=radius) * (T_infinity**4 - T_surface**4))
    
    def resistanceof_cond(self):
        r"""Computes resistance of conduction for a cylindrical object
        
        
        Parameters
        ----------
        `None_required` : 'None'
            Uses attributes defined during instance creation
        
        Returns
        -------
        resistance : `int or float`
            Conduction resistance
                
        Notes
        -----
        The following formula is used:
            
        .. math::
            
            R_{conduction} = \cfrac{ln(r_o/r_i)}{2\pi kL}
        
        *where:*
        
            :math:`r_o` *= outer radius of cylinder*
            
            :math:`r_i` *= inner radius of cylinder*
        
            *L = thickness*
            
            *k = thermal conductivity*
            
            :math:`R_{conduction}` *= conduction resistance*
            
        
        Examples
        --------
        
        First import the module **steadystate**
        
        Units used in this example: SI system
        
        However, any consistent units can be used
        
        >>> from pychemengg.heattransfer import steadystate as ss
        >>> insulation = Cylinder(length=5, inner_radius=1.5e-3, outer_radius=3.5e-3, thermalconductivity=0.15)
        >>> insulation.resistanceof_cond()
        0.1798021903357468
        
        
        References
        ----------
        [1] Yunus A. Cengel and Afshin J. Ghajar, "Heat And Mass Transfer
        Fundamentals and Applications", 6th Edition. New York, McGraw Hill
        Education, 2020        
        """
        # Here complex numbers are used.
        # This is because when 'fsolve' is called, say, to find 'thickness' or
        # 'temperature', then this (i.e. resistanceof_cond) method can get called in the process.
        # In that case, 'fsolve', during internal iterations
        # may put 'inner_radius' and 'outer_radius' as negative, 
        # but math.log cannot handle negative numbers. So instead, here np.log is used.
        # But for np.log to handle negative numbers, the argument must be complex, therefore,
        # internal to this method the radii are transfromed to complex numbers. 
        inner_radius = self.inner_radius + 0j
        outer_radius = self.outer_radius + 0j
        value =  (np.log(outer_radius / inner_radius)) / (2 * math.pi * self.thermalconductivity * self.length)
        return value.real # return only real part
    
    def resistanceof_conv(self, heattransfercoefficient=None, radius=None):
        r"""Computes resistance of convection for a cylindrical object
        
        
        Parameters
        ----------
        heattransfercoefficient : `int or float`
            Heat transfer coefficient *`h`* for the cylindrical surface
        radius : `int or float`
            Radius of cylinder where convective heat transfer is to be computed
        
        
        Returns
        -------
        resistance : `int or float`
            Convection resistance
        
        
        Notes
        -----
        The following formula is used:         
        
        .. math::
            R_{convection} = \frac{1}{hA}
        
        *where:*
        
            *h = heat transfer coefficient*
                    
            *A = area of heat transfer*
            
            :math:`R_{convection}` *= convection resistance*
            
        
        Examples
        --------
        First import the module **steadystate**
        
        Units used in this example: SI system
        
        However, any consistent units can be used
        
        >>> from pychemengg.heattransfer import steadystate as ss
        >>> insulation = Cylinder(length=5, inner_radius=1.5e-3, outer_radius=3.5e-3, thermalconductivity=0.15)
        >>> insulation.resistanceof_conv(heattransfercoefficient=12, radius=insulation.outer_radius)
        0.7578806813899779
        
        
        References
        ----------
        [1] Yunus A. Cengel and Afshin J. Ghajar, "Heat And Mass Transfer
        Fundamentals and Applications", 6th Edition. New York, McGraw Hill
        Education, 2020
        
        """
        return 1/(heattransfercoefficient * self.area(radius=radius))
    
    def resistanceof_fouling(self, foulingfactor=None, radius=None):
        r"""Computes resistance of fouling for a cylindrical object
        
        
        Parameters
        ----------
        foulingfactor : `int or float`
            Fouling factor :math:`R_f` for the cylindrical surface
            
            typical units are :math:`m^2` K/W        
        
        Returns
        -------
        resistance : `int or float`
            Fouling resistance
        
        
        Notes
        -----
        The following formula is used:         
        
        .. math::
            R_{fouling} = \frac{R_f}{A}
        
        *where:*
        
            :math:`R_f` *= fouling factor*
                    
            *A = fouled area of heat transfer*
            
            :math:`R_{fouling}` *= fouling resistance*
            
        
        Examples
        --------
        First import the module **steadystate**
        
        Units used in this example: SI system
        
        However, any consistent units can be used
        
        >>> from pychemengg.heattransfer import steadystate as ss
        >>> insulation = Cylinder(length=5, inner_radius=1.5e-3, outer_radius=3.5e-3, thermalconductivity=0.15)
        >>> insulation.resistanceof_fouling(foulingfactor=0.0007, radius=insulation.outer_radius)
        0.006366197723675814
        
        
        References
        ----------
        [1] Yunus A. Cengel and Afshin J. Ghajar, "Heat And Mass Transfer
        Fundamentals and Applications", 6th Edition. New York, McGraw Hill
        Education, 2020
        
        """        
        return foulingfactor / self.area(radius=radius)
        
    
    
class CylinderBase(object):
    r""" Models the base of a solid circular cylinder
    
    
    Parameters
    ----------
    radius: `int or float`
        Radius of circular base of cylindrical object    
    thermalconductivity : `int or float`
        Thermal conductivity of cylindrical object   
    
    
    Attributes
    ----------
    radius : `int or float`
        Radius of circular base of cylindrical object    
    thermalconductivity : `int or float`
        Thermal conductivity of cylindrical object
    area : `int or float`
        Area of circular base
    
    
    Examples
    --------
    First import the module **steadystate**
    
    Units used in this example: SI system
    
    However, any consistent units can be used
    
    >>> discface = ss.CylinderBase(radius=1.5e-3, thermalconductivity=0.15)
    # This will create an instance of 'CylinderBase' with a name 'discface'  
    """
    
    def __init__(self, radius=None, thermalconductivity=None):
        self.radius = radius
        self.thermalconductivity = thermalconductivity
        self.area =math.pi * self.radius**2
    
    def heatrateof_conv(self, heattransfercoefficient=None, dT=None):
        r"""Computes heat rate of convection for the base of a cylindrical object
        
        
        Parameters
        ----------
        heattransfercoefficient : `int or float`
            Heat transfer coefficient *`h`* for the cylinder base
        radius : `int or float`
            Radius of cylinder base where convective heat transfer rate is to be computed
        dT : `int or float`
            Temperature difference between surface of cylinder base and surrounding fluid
                
        
        Returns
        -------
        heatrate : `int or float`
            Rate of heat transfer by convection
        
        
        Notes
        -----
        Heat rate of convection is calculated using the Newton's Law
            
        .. math::
            Q (heatrate) = h A \Delta T
        
        *where:*
        
            *h = heat transfer coefficient*
        
            *A = area of heat transfer (which is the base circular area)*
             
            :math:`\Delta T` *= temperature difference*
             
        
        Examples
        --------
        
        First import the module **steadystate**
        
        Units used in this example: SI system
        
        However, any consistent units can be used
        
        >>> from pychemengg.heattransfer import steadystate as ss
        >>> disc = ss.CylinderBase(radius=0.1, thermalconductivity=401)
        >>> disc.heatrateof_conv(heattransfercoefficient=132, dT=34)
        140.99467829310993
        
        
        References
        ----------
        [1] Yunus A. Cengel and Afshin J. Ghajar, "Heat And Mass Transfer
        Fundamentals and Applications", 6th Edition. New York, McGraw Hill
        Education, 2020        
        """
        return heattransfercoefficient * self.area * dT
    
    def heatrateof_rad(self, T_infinity=None, T_surface=None, emissivity=None):
        r"""Computes heat rate of radiation for the base of a cylindrical object
        
        
        Parameters
        ----------
        radius : `int or float`
            Radius of cylinder base where radiation heat transfer rate is to be computed
        T_infinity : `int or float`
            Temperature of surroundings in **absolute temperature units**
        T_surface : `int or float`
            Temperature of cylinder base surface in **absolute temperature units**
        emissivity : `int or float`
            Emissivity of the cylinder base
        
        
        Returns
        -------
        heatrate : `int or float` (returns a positive value)
            Rate of heat transfer by radiation
        
        
        Notes
        -----
        Heat rate of radiation is calculated using the Stefan-Boltzmann law
            
        .. math::
            Q (heatrate) = \sigma \epsilon A (T_{infinity}^4 - T_{surface}^4)
        
        *where:*
        
            :math:`\sigma` *= Stefan-Boltzmann constant*
            
            :math:`\epsilon` *= emissivity of object*
        
            *A = area of heat transfer*
             
            :math:`T_{infinity}^4` *= absolute temperature of surroundings*
            
            :math:`T_{surface}^4` *= absolute temperature of surface*
             
        
        Examples
        --------
        
        First import the module **steadystate**
        
        Units used in this example: SI system
        
        However, any consistent units can be used
        
        >>> from pychemengg.heattransfer import steadystate as ss
        >>> disc = ss.CylinderBase(radius=0.1, thermalconductivity=401)
        >>> disc.heatrateof_rad(T_infinity=300, T_surface=550, emissivity=0.9)
        133.72195405218372
        
        
        References
        ----------
        [1] Yunus A. Cengel and Afshin J. Ghajar, "Heat And Mass Transfer
        Fundamentals and Applications", 6th Edition. New York, McGraw Hill
        Education, 2020        
        """  
        return abs(Stefan_Boltzmann * emissivity * self.area * (T_infinity**4 - T_surface**4))
    
    def resistanceof_conv(self, heattransfercoefficient=None):
        r"""Computes resistance of convection for the base of a cylindrical object
        
        
        Parameters
        ----------
        heattransfercoefficient : `int or float`
            Heat transfer coefficient *`h`* for the base surface of cylindrical object
        
        
        Returns
        -------
        resistance : `int or float`
            Convection resistance
        
        
        Notes
        -----
        The following formula is used:         
        
        .. math::
            R_{convection} = \frac{1}{hA}
        
        *where:*
        
            *h = heat transfer coefficient*
                    
            *A = area of heat transfer*
            
            :math:`R_{convection}` *= convection resistance*
            
        
        Examples
        --------
        First import the module **steadystate**
        
        Units used in this example: SI system
        
        However, any consistent units can be used
        
        >>> from pychemengg.heattransfer import steadystate as ss
        >>> disc = ss.CylinderBase(radius=0.1, thermalconductivity=401)
        >>> disc.resistanceof_conv(heattransfercoefficient=132)
        0.24114385316953837
        
        
        References
        ----------
        [1] Yunus A. Cengel and Afshin J. Ghajar, "Heat And Mass Transfer
        Fundamentals and Applications", 6th Edition. New York, McGraw Hill
        Education, 2020
        
        """
        return 1/(heattransfercoefficient * self.area)
    
    def resistanceof_fouling(self, foulingfactor=None):
        r"""Computes resistance of fouling for the base of cylindrical object
        
        
        Parameters
        ----------
        foulingfactor : `int or float`
            Fouling factor :math:`R_f` for the cylindrical base surface
            
            typical units are :math:`m^2` K/W        
        
        Returns
        -------
        resistance : `int or float`
            Fouling resistance
        
        
        Notes
        -----
        The following formula is used:         
        
        .. math::
            R_{fouling} = \frac{R_f}{A}
        
        *where:*
        
            :math:`R_f` *= fouling factor*
                    
            *A = fouled area of heat transfer*
            
            :math:`R_{fouling}` *= fouling resistance*
            
        
        Examples
        --------
        First import the module **steadystate**
        
        Units used in this example: SI system
        
        However, any consistent units can be used
        
        >>> from pychemengg.heattransfer import steadystate as ss
        >>> disc = ss.CylinderBase(radius=0.1, thermalconductivity=401)
        >>> disc.resistanceof_fouling(foulingfactor=0.0007)
        0.022281692032865345
        
        
        References
        ----------
        [1] Yunus A. Cengel and Afshin J. Ghajar, "Heat And Mass Transfer
        Fundamentals and Applications", 6th Edition. New York, McGraw Hill
        Education, 2020
        
        """  
        return foulingfactor / self.area
        
    
class Sphere(object):
    r""" Models a spherical object
    
    
    Parameters
    ----------
    inner_radius: `int or float`
        Inner radius of spherical object
    outer_radius: `int or float`
        Outer radius of spherical object    
    thermalconductivity : `int or float`
        Thermal conductivity of spherical object   
    
    
    Attributes
    ----------
    inner_radius: `int or float`
        Inner radius of spherical object
    outer_radius: `int or float`
        Outer radius of spherical object    
    thermalconductivity : `int or float`
        Thermal conductivity of spherical object 
    
    
    Examples
    --------
    First import the module **steadystate**
    
    Units used in this example: SI system
    
    However, any consistent units can be used
    
    >>> from pychemengg.heattransfer import steadystate as ss
    >>> shellLead = ss.Sphere(inner_radius=0.25, outer_radius=0.30, thermalconductivity=35.3)
    # This will create an instance of 'Sphere' with a name 'shellLead'  
    """
    
    def __init__(self, inner_radius=None, outer_radius=None, thermalconductivity=None):
        self.inner_radius = inner_radius
        self.outer_radius = outer_radius
        self.thermalconductivity = thermalconductivity
    
    def area(self, radius=None):
        r"""Computes surface area of a sphere
        
        
        Parameters
        ----------
        radius : `int or float`
            Radius at which surface area is to be computed
                
        
        Returns
        -------
        area : `int or float`
            Surface area of sphere
        
        
        Notes
        -----
        Surface area of sphere is computed using:
            
        .. math::
            A = 4 \pi r^2
                    
        *where:*
        
            *r = radius*
             
            *A = surface area*
        
        
        Examples
        --------
        First import the module **steadystate**
        
        Units used in this example: SI system
        
        However, any consistent units can be used
        
        >>> from pychemengg.heattransfer import steadystate as ss
        
        >>> shellLead = ss.Sphere(inner_radius=0.25, outer_radius=0.30, thermalconductivity=35.3)
        # This will create an instance of 'Sphere' with a name 'shellLead'
        
        >>> shellLead.area(radius=shellLead.outer_radius)
        1.1309733552923256
        # This computes surface area @ outer radius
        
        >>> shellLead.area(radius=shellLead.inner_radius)
        0.7853981633974483
        # This computes surface area @ inner radius
        
        >>> shellLead.area(radius=0.25)
        0.7853981633974483
        # This will also compute surface area @ inner radius, which is = 0.25
        # Here the radius is entered as a number directly, rather than as an attribute.
        """
        return 4 * math.pi * radius * radius   
    
    def volume(self, radius=None):
        r"""Calculates volume of a sphere
        
        Parameters
        ----------
        radius : `int or float`
            Radius at which volume is to be computed
                
        
        Returns
        -------
        volume : `int or float`
            Volume of sphere
        
        
        Notes
        -----
        Volume of sphere is computed using:
            
        .. math::
            
            V = \frac{4}{3} \pi r^3
                    
        *where:*
        
            *r = radius*
             
            *V = volume*
        
        
        Examples
        --------
        
        First import the module **steadystate**
        
        Units used in this example: SI system
        
        However, any consistent units can be used
        
        >>> from pychemengg.heattransfer import steadystate as ss
        
        >>> shellLead = ss.Sphere(inner_radius=0.25, outer_radius=0.30, thermalconductivity=35.3)
        # This will create an instance of 'Sphere' with a name 'shellLead'
        
        >>> shellLead.volume(radius=shellLead.outer_radius)
        0.11309733552923253
        # This computes volume @ outer radius
        
        >>> shellLead.volume(radius=shellLead.inner_radius)
        0.06544984694978735
        # This computes volume @ inner radius
        
        >>>  shellLead.volume(radius=0.25)
        0.06544984694978735
        # This will also compute volume @ inner radius, which is = 0.25.
        # Here the radius is entered as a number directly, rather than as an attribute.
        """
        return 4/3 * math.pi * radius**3
    
    def heatrateof_cond(self, dT=None):
        r"""Computes heat rate of conduction for a spherical object
        
        
        Parameters
        ----------
        dT : `int or float`
            Temperature difference between two surfaces of the sphere
                
        
        Returns
        -------
        heatrate : `int or float`
            Rate of heat transfer by conduction
        
        
        Notes
        -----
        The following formula is used:
            
        .. math::
            Q (heatrate) = \frac{\Delta T}{R_{conduction}}
        
        *where:*
             
            :math:`\Delta T` *= temperature difference*
            
            :math:`R_{conduction}` *= conduction resistance* given by
                         
                :math:`R_{conduction} = \cfrac{(r_o - r_i)}{4\pi k r_o r_i}`
        
                :math:`r_o` *= outer radius of sphere*
                
                :math:`r_i` *= inner radius of sphere*
                
                *k = thermal conductivity*
        
        
        Examples
        --------
        First import the module **steadystate**
        
        Units used in this example: SI system
        
        However, any consistent units can be used
        
        >>> from pychemengg.heattransfer import steadystate as ss
        >>> shellSS = ss.Sphere(inner_radius=.30, outer_radius=.31, thermalconductivity=15.1)
        >>> shellSS.heatrateof_cond(dT=75)
        132352.15690308428
        
        References
        ----------
        [1] Yunus A. Cengel and Afshin J. Ghajar, "Heat And Mass Transfer
        Fundamentals and Applications", 6th Edition. New York, McGraw Hill
        Education, 2020        
        """
        return dT/self.resistanceof_cond()
    
    def heatrateof_conv(self, heattransfercoefficient=None, radius=None, dT=None):
        r"""Computes heat rate of convection for a spherical object
        
        
        Parameters
        ----------
        heattransfercoefficient : `int or float`
            Heat transfer coefficient *`h`* for the spherical surface
        radius : `int or float`
            Radius of sphere where convective heat transfer rate is to be computed
        dT : `int or float`
            Temperature difference between spherical surface and surrounding fluid
                
        
        Returns
        -------
        heatrate : `int or float`
            Rate of heat transfer by convection
        
        
        Notes
        -----
        Heat rate of convection is calculated using the Newton's Law
            
        .. math::
            Q (heatrate) = h A \Delta T
        
        *where:*
        
            *h = heat transfer coefficient*
        
            *A = area of heat transfer*
             
            :math:`\Delta T` *= temperature difference*
             
        
        Examples
        --------
        
        First import the module **steadystate**
        
        Units used in this example: SI system
        
        However, any consistent units can be used
        
        >>> from pychemengg.heattransfer import steadystate as ss
        >>> shellSS = ss.Sphere(inner_radius=.30, outer_radius=.31, thermalconductivity=15.1)
        >>> shellSS.heatrateof_conv(heattransfercoefficient=500, radius=shellSS.outer_radius, dT=25)
        15095.352700498957
        
        References
        ----------
        [1] Yunus A. Cengel and Afshin J. Ghajar, "Heat And Mass Transfer
        Fundamentals and Applications", 6th Edition. New York, McGraw Hill
        Education, 2020        
        """
        return heattransfercoefficient * self.area(radius=radius) * dT
    
    def heatrateof_rad(self, radius=None, T_infinity=None, T_surface=None, emissivity=None):
        r"""Computes heat rate of radiation for a spherical object
        
        
        Parameters
        ----------
        radius : `int or float`
            Radius of sphere where radiation heat transfer rate is to be computed
        T_infinity : `int or float`
            Temperature of surroundings in **absolute temperature units**
        T_surface : `int or float`
            Temperature of sphere's surface in **absolute temperature units**
        emissivity : `int or float`
            Emissivity of the sphere
        
        
        Returns
        -------
        heatrate : `int or float` (returns a positive value)
            Rate of heat transfer by radiation
        
        
        Notes
        -----
        Heat rate of radiation is calculated using the Stefan-Boltzmann law
            
        .. math::
            Q (heatrate) = \sigma \epsilon A (T_{infinity}^4 - T_{surface}^4)
        
        *where:*
        
            :math:`\sigma` *= Stefan-Boltzmann constant*
            
            :math:`\epsilon` *= emissivity of object*
        
            *A = area of heat transfer*
             
            :math:`T_{infinity}^4` *= absolute temperature of surroundings*
            
            :math:`T_{surface}^4` *= absolute temperature of surface*
             
        
        Examples
        --------
        
        First import the module **steadystate**
        
        Units used in this example: SI system
        
        However, any consistent units can be used
        
        >>> from pychemengg.heattransfer import steadystate as ss
        >>> shellSS = ss.Sphere(inner_radius=.30, outer_radius=.31, thermalconductivity=15.1)
        >>> shellSS.heatrateof_rad(radius=shellSS.outer_radius, T_infinity=30+273, T_surface=90.6+273, emissivity=0.95)
        588.6831572504626
        
        
        References
        ----------
        [1] Yunus A. Cengel and Afshin J. Ghajar, "Heat And Mass Transfer
        Fundamentals and Applications", 6th Edition. New York, McGraw Hill
        Education, 2020        
        """
        return abs(Stefan_Boltzmann*emissivity*self.area(radius=radius)*(T_infinity**4 - T_surface**4))

    def resistanceof_cond(self):
        r"""Computes resistance of conduction for a spherical object
        
        
        Parameters
        ----------
        `None_required` : 'None'
            Uses attributes defined during instance creation
        
        Returns
        -------
        resistance : `int or float`
            Conduction resistance
                
        Notes
        -----
        The following formula is used:
            
        .. math::
            
            R_{conduction} = \cfrac{(r_o - r_i)}{4\pi k r_o r_i}
        
        *where:*
        
            :math:`r_o` *= outer radius of sphere*
            
            :math:`r_i` *= inner radius of sphere*
            
            *k = thermal conductivity*
            
            :math:`R_{conduction}` *= conduction resistance*
            
        
        Examples
        --------
        
        First import the module **steadystate**
        
        Units used in this example: SI system
        
        However, any consistent units can be used
        
        >>> from pychemengg.heattransfer import steadystate as ss
        >>> shellSS = ss.Sphere(inner_radius=.30, outer_radius=.31, thermalconductivity=15.1)
        >>> shellSS.resistanceof_cond()
        0.0005666700245385441
        
        
        References
        ----------
        [1] Yunus A. Cengel and Afshin J. Ghajar, "Heat And Mass Transfer
        Fundamentals and Applications", 6th Edition. New York, McGraw Hill
        Education, 2020        
        """
        return (self.outer_radius - self.inner_radius) / (4 * math.pi * self.thermalconductivity * self.outer_radius * self.inner_radius)
    
    def resistanceof_conv(self, heattransfercoefficient=None, radius=None):
        r"""Computes resistance of convection for a spherical object
        
        
        Parameters
        ----------
        heattransfercoefficient : `int or float`
            Heat transfer coefficient *`h`* for the spherical surface
        radius : `int or float`
            Radius of sphere where convective heat transfer is to be computed
        
        
        Returns
        -------
        resistance : `int or float`
            Convection resistance
        
        
        Notes
        -----
        The following formula is used:         
        
        .. math::
            R_{convection} = \frac{1}{hA}
        
        *where:*
        
            *h = heat transfer coefficient*
                    
            *A = area of heat transfer*
            
            :math:`R_{convection}` *= convection resistance*
            
        
        Examples
        --------
        First import the module **steadystate**
        
        Units used in this example: SI system
        
        However, any consistent units can be used
        
        >>> from pychemengg.heattransfer import steadystate as ss
        >>> shellSS = ss.Sphere(inner_radius=.30, outer_radius=.31, thermalconductivity=15.1)
        >>> shellSS.resistanceof_conv(heattransfercoefficient=500, radius=shellSS.outer_radius)
        0.0016561388459094206
        
        
        References
        ----------
        [1] Yunus A. Cengel and Afshin J. Ghajar, "Heat And Mass Transfer
        Fundamentals and Applications", 6th Edition. New York, McGraw Hill
        Education, 2020
        
        """
        return 1/(heattransfercoefficient * self.area(radius=radius))

    def resistanceof_fouling(self, foulingfactor=None, radius=None):
        r"""Computes resistance of fouling for a spherical object
        
        
        Parameters
        ----------
        foulingfactor : `int or float`
            Fouling factor :math:`R_f` for the spherical surface
            
            typical units are :math:`m^2` K/W        
        
        Returns
        -------
        resistance : `int or float`
            Fouling resistance
        
        
        Notes
        -----
        The following formula is used:         
        
        .. math::
            R_{fouling} = \frac{R_f}{A}
        
        *where:*
        
            :math:`R_f` *= fouling factor*
                    
            *A = fouled area of heat transfer*
            
            :math:`R_{fouling}` *= fouling resistance*
            
        
        Examples
        --------
        First import the module **steadystate**
        
        Units used in this example: SI system
        
        However, any consistent units can be used
        
        >>> from pychemengg.heattransfer import steadystate as ss
        >>> shellSS = ss.Sphere(inner_radius=.30, outer_radius=.31, thermalconductivity=15.1)
        >>> shellSS.resistanceof_fouling(foulingfactor=0.0007, radius=shellSS.outer_radius)
        0.0005796485960682973
        
        
        References
        ----------
        [1] Yunus A. Cengel and Afshin J. Ghajar, "Heat And Mass Transfer
        Fundamentals and Applications", 6th Edition. New York, McGraw Hill
        Education, 2020
        
        """
        return foulingfactor / self.area(radius=radius)
