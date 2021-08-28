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
        
            *k = thermal conductivity*
            
            :math:`L_c` *= characteristic length*
            
            *Bi = Biot number*
                     
        
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
        # Next call calc_Bi
        >>> plate.calc_Bi()
        0.0029444444444444444
        
        
        References
        ----------
        [1] Yunus A. Cengel and Afshin J. Ghajar, "Heat And Mass Transfer
        Fundamentals and Applications", 6th Edition. New York, McGraw Hill
        Education, 2020.       
        """        
        self.Bi = (self.heattransfercoefficient
                           * self.characteristiclength
                           / self.thermalconductivity)
        return self.Bi
    
    def temp_of_solid_at_time_t(self, time=None):
        
        r"""Temperature of solid at a given time = t.
        
        
        Parameters
        ----------
        time : 'int or float'
            Time from start of process at which temperature is required to be found.
                                
        
        Returns
        -------
        temperature : `int or float`
            Temperature of object at time = t 
        
        
        Notes
        -----
        Temperature of object at time = t is calculated using the following formula:
            
        .. math::
            T(t) = T_{infinity} + (T_i - T_{infinity}) e^{-bt}
        
        *where:*
        
            *h = heat transfer coefficient*
        
            *k = thermal conductivity*
            
            :math:`L_c` *= characteristic length*
            
            *Bi = Biot number*
                     
        
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
        # Next call calc_Bi
        >>> plate.calc_Bi()
        0.0029444444444444444
        
        
        References
        ----------
        [1] Yunus A. Cengel and Afshin J. Ghajar, "Heat And Mass Transfer
        Fundamentals and Applications", 6th Edition. New York, McGraw Hill
        Education, 2020.       
        """
        b = ((self.heattransfercoefficient * self.surfacearea)
             / (self.density*self.volume*self.specificheat)) # unit is 1/time
        solidtemp_at_time_t = (self.T_infinity
                                + (self.T_initial-self.T_infinity)*np.exp(-b*time+0j))
        return solidtemp_at_time_t.real
    
    def heatrateof_conv_at_time_t(self):
        q_rate = (self.heattransfercoefficient * self.surfacearea
                *(self.solidtemp_at_time_t - self.T_infinity))
        return q_rate

    def totalheat_transferred_during_interval_t(self):
        qtotal = (self.mass * self.specificheat
                 * (self.solidtemp_at_time_t - self.T_initial))
        return qtotal
    
    def maxheattransferpossible(self):
        qtotal_max = self.mass * self.specificheat * (self.T_infinity - self.T_initial)
        return qtotal_max





class NonLumpedSlab():
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
        self.Bi = (self.heattransfercoefficient
                           * self.thickness/2
                           / self.thermalconductivity)
        return self.Bi

    
    def Fo(self, time=None):
        self.Fo = self.thermaldiffusivity*time/(self.thickness/2)**2
        return self.Fo
    
    
    def calc_eigenvalues(self, numberof_eigenvalues_desired=None):
        slab_eigenfunction = lambda x, Bi: x*np.tan(x)-Bi
        slab_eigenvalues = _get_eigenvalues(slab_eigenfunction, Bi=self.Bi,
                                          numberof_eigenvalues_desired=numberof_eigenvalues_desired)
        self.eigenvalues = np.array(slab_eigenvalues)
        return self.eigenvalues

        
    def temp_of_solid_at_time_t(self, time=None, xposition_tofindtemp=None):
        term1 = 4*np.sin(self.eigenvalues)
        term2 = 2*self.eigenvalues + np.sin(2*self.eigenvalues)
        term3 = np.exp(-np.power(self.eigenvalues,2) * self.Fo)
        term4 = np.cos(self.eigenvalues*xposition_tofindtemp/(self.thickness/2))
        theta = np.sum(term1/term2*term3*term4)
        self.solidtemp_at_time_t = (self.T_infinity
                                + (self.T_initial-self.T_infinity)*theta)
        return self.solidtemp_at_time_t
    

    def heatrateof_conv_at_time_t(self):
        qrate = (self.heattransfercoefficient * self.surfacearea
                *(self.solidtemp_at_time_t - self.T_infinity))
        return qrate


    def totalheat_transferred_during_interval_t(self, time=None, xposition_tofindtemp=None):
        term1 = 4*np.sin(self.eigenvalues)
        term2 = 2*self.eigenvalues + np.sin(2*self.eigenvalues)
        term3 = np.exp(-np.power(self.eigenvalues,2) * self.Fo)        
        term4 = np.sin(self.eigenvalues)/self.eigenvalues
        normalized_heatamount = 1 - np.sum(term1/term2*term3*term4)
        heattransferred = self.maxheattransferpossible() * normalized_heatamount
        return heattransferred  
    
    def maxheattransferpossible(self):
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

    
    def Fo(self, time=None):
        self.Fo = self.thermaldiffusivity * time / (self.radius)**2
        return self.Fo
    
    
    def calc_eigenvalues(self, numberof_eigenvalues_desired=None):
        cylinder_eigenfunction = lambda x, Bi: x*j1(x)/j0(x)-Bi
        cylinder_eigenvalues = _get_eigenvalues(cylinder_eigenfunction, Bi=self.Bi,
                                          numberof_eigenvalues_desired=numberof_eigenvalues_desired)
        self.eigenvalues = np.array(cylinder_eigenvalues)
        return self.eigenvalues

        
    def temp_of_solid_at_time_t(self, time=None, rposition_tofindtemp=None):
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

    
    def Fo(self, time=None):
        self.Fo = self.thermaldiffusivity * time / (self.radius)**2
        return self.Fo
    
    
    def calc_eigenvalues(self, numberof_eigenvalues_desired=None):
        sphere_eigenfunction = lambda x,Bi: 1-x/np.tan(x)-Bi
        sphere_eigenvalues = _get_eigenvalues(sphere_eigenfunction, Bi=self.Bi,
                                          numberof_eigenvalues_desired=numberof_eigenvalues_desired)
        self.eigenvalues = np.array(sphere_eigenvalues)
        return self.eigenvalues

        
    def temp_of_solid_at_time_t(self, time=None, rposition_tofindtemp=None):
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
        return (plate.temp_of_solid_at_time_t(time))-50
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
    plate.temp_of_solid_at_time_t(time=7*60, xposition_tofindtemp=plate.thickness/2)

    
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
    cylinder.temp_of_solid_at_time_t(time=45*60, rposition_tofindtemp=0)
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
    #     temp_center = tomato.temp_of_solid_at_time_t(time=2*3600, rposition_tofindtemp=0)
    #     temp_surface = tomato.temp_of_solid_at_time_t(time=2*3600, rposition_tofindtemp=0.04)
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
        temp_center = potato.temp_of_solid_at_time_t(time=time[0], rposition_tofindtemp=0)
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
    wall.temp_of_solid_at_time_t(time=8*60, xposition_tofindtemp=0)
    wall.temp_of_solid_at_time_t(time=8*60, xposition_tofindtemp=40e-3)
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
#         T0 = wall.temp_of_solid_at_time_t(time=time[0], xposition_tofindtemp=0)
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
#         T0 = wall.temp_of_solid_at_time_t(time=time[0], xposition_tofindtemp=0)
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
#         T0 = bearings.temp_of_solid_at_time_t(time=time[0],rposition_tofindtemp=0)
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
        temp_center = egg.temp_of_solid_at_time_t(time=time[0], rposition_tofindtemp=0)
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