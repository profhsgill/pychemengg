# -*- coding: utf-8 -*-

# This file is part of PyChemEngg python package.
 
# PyChemEngg: A python framework to promote problem solving and critical
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

"""Module to compute physical properties of air.

"""
import numpy as np
import os

def _loadair_data():
    __location__ = os.path.realpath(
                        os.path.join(os.getcwd(), os.path.dirname(__file__)))
    air = []
    file_to_open = "data_airproperties.txt"
    with open(os.path.join(__location__, file_to_open),"r") as air_file:
        for line in air_file:
            air.append([float(eval(x)) if x != "nan" else np.nan for x in line.strip().split()])
    air_file.close()
    air_data = np.array(air)
    return air_data

def density(T=None):
    r""" Provides density of air at a temperature T.
    
    
    Parameters
    ----------
    T : `int or float`
        Temperature in 'Celsius' at which density is required.


    Returns
    -------
    density : `int or float`
       Density (kg/m3) at temperature T .
   
    
    Notes
    -----
    Look up table adapted from ref [1].
    
    Linear interpolation is performed when the temperature lies
    between tabulated entries.

    
    Examples
    --------
    First import the module **airproperties**.
       
    >>> from pychemengg.physicalproperties import airproperties as ap 
    >>> ap.density(T=42.5)
    1.1179999999999999

    
    References
    ----------
    [1] Yunus A. Cengel and Afshin J. Ghajar, "Heat And Mass Transfer
    Fundamentals and Applications", 6th Edition. New York, McGraw Hill
    Education, 2020.
    
    """
    air_data = _loadair_data()
    return np.interp(T, air_data[:,0], air_data[:,1])

def specificheat(T=None):
    r""" Provides specific heat of air at a temperature T.
    
    
    Parameters
    ----------
    T : `int or float`
        Temperature in 'Celsius' at which specific heat is required.


    Returns
    -------
    specific heat : `int or float`
        Specific heat (J/kg K) at temperature T.
   
    
    Notes
    -----
    Look up table adapted from ref [1].
    
    Linear interpolation is performed when the temperature lies
    between tabulated entries.

    
    Examples
    --------
    First import the module **airproperties**.
       
    >>> from pychemengg.physicalproperties import airproperties as ap 
    >>> ap.specificheat(T=42.5)
    1007.0

    
    References
    ----------
    [1] Yunus A. Cengel and Afshin J. Ghajar, "Heat And Mass Transfer
    Fundamentals and Applications", 6th Edition. New York, McGraw Hill
    Education, 2020.
    
    """
    air_data = _loadair_data()
    return np.interp(T, air_data[:,0], air_data[:,2])

def thermalconductivity(T=None):
    r""" Provides thermal conductivity of air at a temperature T.
    
    
    Parameters
    ----------
    T : `int or float`
        Temperature in 'Celsius' at which thermal conductivity is required.


    Returns
    -------
    Thermal conductivity : `int or float`
        Thermal conductivity (W/mK) at temperature T.
   
    
    Notes
    -----
    Look up table adapted from ref [1].
    
    Linear interpolation is performed when the temperature lies
    between tabulated entries.

    
    Examples
    --------
    First import the module **airproperties**.
       
    >>> from pychemengg.physicalproperties import airproperties as ap 
    >>> ap.thermalconductivity(T=42.5)
    0.026805000000000002

    
    References
    ----------
    [1] Yunus A. Cengel and Afshin J. Ghajar, "Heat And Mass Transfer
    Fundamentals and Applications", 6th Edition. New York, McGraw Hill
    Education, 2020.
    
    """
    
    air_data = _loadair_data()
    return np.interp(T, air_data[:,0], air_data[:,3])

def thermaldiffusivity(T=None):
    r""" Provides thermal diffusivity of air at a temperature T.
    
    
    Parameters
    ----------
    T : `int or float`
        Temperature in 'Celsius' at which thermal diffusivity is required.


    Returns
    -------
    Thermal diffusivity: `int or float`
        Thermal diffusivity (m2/s2) at temperature T.
   
    
    Notes
    -----
    Look up table adapted from ref [1].
    
    Linear interpolation is performed when the temperature lies
    between tabulated entries.

    
    Examples
    --------
    First import the module **airproperties**.
       
    >>> from pychemengg.physicalproperties import airproperties as ap 
    >>> ap.thermaldiffusivity(T=42.5)
    2.3810000000000004e-05

    
    References
    ----------
    [1] Yunus A. Cengel and Afshin J. Ghajar, "Heat And Mass Transfer
    Fundamentals and Applications", 6th Edition. New York, McGraw Hill
    Education, 2020.
    
    """
    air_data = _loadair_data()
    return np.interp(T, air_data[:,0], air_data[:,4])

def viscosity(T=None):
    r""" Provides viscosity of air at a temperature T
    
    
    Parameters
    ----------
    T : `int or float`
        Temperature in 'Celsius' at which viscosity (dynamic) is required.


    Returns
    -------
    Viscosity: `int or float`
        Dynamic viscosity (kg/ms) at temperature T.
   
    
    Notes
    -----
    Look up table adapted from ref [1].
    
    Linear interpolation is performed when the temperature lies
    between tabulated entries.

    
    Examples
    --------
    First import the module **airproperties**.
       
    >>> from pychemengg.physicalproperties import airproperties as ap 
    >>> ap.viscosity(T=42.5)
    1.9295e-05

    
    References
    ----------
    [1] Yunus A. Cengel and Afshin J. Ghajar, "Heat And Mass Transfer
    Fundamentals and Applications", 6th Edition. New York, McGraw Hill
    Education, 2020.
    """
    air_data = _loadair_data()
    return np.interp(T, air_data[:,0], air_data[:,5])


def kinematicviscosity(T=None):
    r""" Provides kinematicviscosity of air at a temperature T
    
    
    Parameters
    ----------
    T : `int or float`
        Temperature in 'Celsius' at which kinematic viscosity is required.


    Returns
    -------
    Kinematicviscosity: `int or float`
        Kinematic viscosity (m2/s) at temperature T.
   
    
    Notes
    -----
    Look up table adapted from ref [1].
    
    Linear interpolation is performed when the temperature lies
    between tabulated entries.

    
    Examples
    --------
    First import the module **airproperties**.
       
    >>> from pychemengg.physicalproperties import airproperties as ap 
    >>> ap.kinematicviscosity(T=42.5)
    1.726e-05

    
    References
    ----------
    [1] Yunus A. Cengel and Afshin J. Ghajar, "Heat And Mass Transfer
    Fundamentals and Applications", 6th Edition. New York, McGraw Hill
    Education, 2020.
    """
    air_data = _loadair_data()
    return np.interp(T, air_data[:,0], air_data[:,6])


def prandtlnumber(T=None):
    r""" Provides Prandtl number of air at a temperature T
    
    
    Parameters
    ----------
    T : `int or float`
        Temperature in 'Celsius' at which Prandtl number is required.


    Returns
    -------
    Prandtl number: `int or float`
        Prandtl number at temperature T.
   
    
    Notes
    -----
    Look up table adapted from ref [1].
    
    Linear interpolation is performed when the temperature lies
    between tabulated entries.

    
    Examples
    --------
    First import the module **airproperties**.
       
    >>> from pychemengg.physicalproperties import airproperties as ap 
    >>> ap.prandtlnumber(T=42.5)
    0.7248

    
    References
    ----------
    [1] Yunus A. Cengel and Afshin J. Ghajar, "Heat And Mass Transfer
    Fundamentals and Applications", 6th Edition. New York, McGraw Hill
    Education, 2020.
    """
    air_data = _loadair_data()
    return np.interp(T, air_data[:,0], air_data[:,7])

# Details of table loaded from .txt file
# col_0 = 'Temp (C)'
# col_1 = 'Density kg/m3'
# col_2 = 'Specific heat J/kg K'
# col_3 = 'Thermal conductivity  W/mK'
# col_4 = 'Thermal diffusivity m2/s2'
# col_5 = 'Dynamic viscosity kg/m s'
# col_6 = 'Kinematic viscosity m2/s'
# col_7 = 'Prandtl Number'