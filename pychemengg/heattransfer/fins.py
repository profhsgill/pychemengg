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

"""Calculates fin efficiency 
"""

import math
import scipy.special

__all__ = ["Fin"]


class Fin(object):
    r""" To compute fin efficiency of different fin types.

    Parameters
    ----------
    length : `int or float`
        Length of fin.
    width : `int or float`
        Width of fin.
    thickness : `int or float`
        Thickness of fin.
    diameter : `int or float`
        Diameter of fin.
    inner_radius : `int or float`
        Inner radius of annular fin.
    outer_radius : `int or float`
        Outer radius of annular fin.
    heattransfercoefficient : `int or float`
        Heat transfer coefficient between fin and surroundings.
    thermalconductivity : `int or float`
        Thermal conductivity of fin material.   
 
    
    Attributes
    ----------
    See "Parameters". All parameters are attributes. Additional attributes are listed below.
    efficiency : `int or float`
        Fin efficiency.
    surfacearea : `int or float`
        Surface area of given fin geometry.
    
        
    Examples
    --------
    First import the module **fins**.
       
    >>> from heattransfer import fins as fins 
    >>> fin1 = fins.Fin(keyword1=..., keyword2=..., ...)
    # This will create an instance of class 'Fin' and assign it to the
    variable 'fin1'.
    # Methods of the class 'Fin' can then be called like so :-
    # fin1.cylindrical()
    # or
    # fin1.rectangular()

    """
    
    def __init__ (self, length=None, width=None, thickness=None,
                  diameter=None, inner_radius=None, outer_radius=None,
                  heattransfercoefficient=None, thermalconductivity=None):
        self.length = length # fin length
        self.width = width # fin width
        self.thickness = thickness # fin thickness
        self.diameter = diameter # fin diameter
        self.inner_radius = inner_radius # for annular fin: fin inner radius
        self.outer_radius = outer_radius # for annular fin: fin outer radius
        self.heattransfercoefficient = heattransfercoefficient # heat transfer coefficient
        self.thermalconductivity = thermalconductivity # thermal conductivity of fin

    def cylindrical(self):
        r""" Computes fin efficiency for a cylindrical fin.

        Parameters
        ----------
        `None_required` : 'None'
            Attributes that are already defined are used in calculation.
                
        
        Returns
        -------
        tuple containing efficiency and surfacearea :
            (efficiency : `int or float`, surfacearea : `int or float`)
        
        
        Notes
        -----
        The following formula is used:         
        
        .. math::
            
            \eta = \frac {\tanh (mL)} {mL} \\[15pt]
            
            \\[15pt]
            
            A_{surface area} = Perimeter \hspace{2pt} L

        *where:*
            
            L = fin length
            
            m = :math:`\sqrt {\frac {h \hspace{2pt} Perimeter} {k \hspace{2pt} A_{cross-section}}}`
            
            h = heat transfer coefficient
            
            Perimeter = fin perimeter
            
            :math:`A_{cross-section}` = cross section area of fin
            
            :math:`A_{surface area}` = surface area of fin
        
        
        Examples
        --------
        First import the module **fins**.
           
        >>> from pychemengg.heattransfer import fins as fins
        >>> fin1 = fins.Fin(length=0.06,
                             diameter=5e-3,
                             heattransfercoefficient=35,
                             thermalconductivity=8.3)
        >>> fin1.cylindrical()
        (0.2864128219415127, 0.000942477796076938)
        

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
        # def cylindrical(self, length=1, diameter=1, heattransfercoefficient=1, thermalconductivity=1):
        try:
            self.perimeter = math.pi*self.diameter
            self.crosssectionarea = math.pi/4*self.diameter**2 # Area of cross section
            self.surfacearea = self.perimeter * self.length # surface area of fin
            self.mL = math.sqrt(self.heattransfercoefficient * self.perimeter / self.thermalconductivity / self.crosssectionarea) * self.length
            self.efficiency = math.tanh(self.mL)/self.mL
            return self.efficiency, self.surfacearea
        except TypeError:
            print("Check your 'Fin()' object. You maynot have defined all required parameters.")
            print("\ncylindrical fin requires: length=someValue, diameter=someValue,, heattransfercoefficient=someValue, thermalconductivity=someValue")

    def rectangular(self):
        r""" Computes fin efficiency for a rectangular fin.
        

        Parameters
        ----------
        `None_required` : 'None'
            Attributes that are already defined are used in calculation.
                
        
        Returns
        -------
        tuple containing efficiency and surfacearea :
            (efficiency : `int or float`, surfacearea : `int or float`)
        
        
        
        Examples
        --------
        First import the module **fins**.
           
        >>> from pychemengg.heattransfer import fins as fins
        >>> fin1 = fins.Fin (length=0.064, width=1, thickness=0.008, heattransfercoefficient=25, thermalconductivity=210)
        >>> fin1.rectangular()
        (0.9609578814207821, 0.129024)
        

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
        # def rectangular(self, length=1, width=1, thickness=1, heattransfercoefficient=1, thermalconductivity=1):
        try:
            self.perimeter = 2 * (self.width + self.thickness); # perimeter
            self.crosssectionarea = self.width * self.thickness; # Area of cross section
            self.surfacearea = self.perimeter * self.length; # surface area of fin
            self.mL = math.sqrt(self.heattransfercoefficient * self.perimeter / self.thermalconductivity / self.crosssectionarea) * self.length;
            self.efficiency = math.tanh(self.mL)/self.mL;      
            return self.efficiency, self.surfacearea
        except TypeError:
            print("Check your 'Fin()' object. You maynot have defined all required parameters.")
            print("\nrectangular fin requires: length=someValue, width=SomeValue, thickness=someValue, heattransfercoefficient=someValue, thermalconductivity=someValue")

    def straighttriangular(self):
        r""" Computes fin efficiency for a straighttriangular fin.
        

        Parameters
        ----------
        `None_required` : 'None'
            Attributes that are already defined are used in calculation.
                
        
        Returns
        -------
        tuple containing efficiency and surfacearea :
            (efficiency : `int or float`, surfacearea : `int or float`)
        
        
        
        Examples
        --------
        First import the module **fins**.
           
        >>> from pychemengg.heattransfer import fins as fins
        >>> fin1 = fins.Fin(length=0.055, width=0.11, thickness=0.004, heattransfercoefficient=25, thermalconductivity=236)
        >>> fin1.straighttriangular()
        (0.9275971984358753, 0.012107997357118972)
        

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
        # def straighttriangular(self, length=1, width=1, thickness=1, heattransfercoefficient=1, thermalconductivity=1):
        try:
            self.surfacearea = 2 * self.width * math.sqrt(self.length**2 + (self.thickness/2)**2) #surface area of fin
            self.mL = math.sqrt(2 * self.heattransfercoefficient / self.thermalconductivity / self.thickness) * self.length
            self.efficiency = scipy.special.iv(1,2*self.mL) / (self.mL*scipy.special.iv(0,2*self.mL))
            return self.efficiency, self.surfacearea
        except TypeError:
            print("Check your 'Fin()' object. You maynot have defined all required parameters.")
            print("\nstraighttriangular fin requires: length=someValue, width=SomeValue, thickness=someValue, heattransfercoefficient=someValue, thermalconductivity=someValue")

    def straightparabolic(self):
        r""" Computes fin efficiency for a straightparabolic fin.
        

        Parameters
        ----------
        `None_required` : 'None'
            Attributes that are already defined are used in calculation.
                
        
        Returns
        -------
        tuple containing efficiency and surfacearea :
            (efficiency : `int or float`, surfacearea : `int or float`)
        
        
        Examples
        --------
        First import the module **fins**.
           
        >>> from pychemengg.heattransfer import fins as fins
        >>> fin1 = fins.Fin(length=0.03, width=0.2, thickness=0.001, heattransfercoefficient=10, thermalconductivity=125)
        >>> fin1.straightparabolic()
        (0.8867652295764488, 0.012002221851998751)
        

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
    # def straightparabolic(self, length=1, width=1, thickness=1, heattransfercoefficient=1, thermalconductivity=1):
        try:
            self.C1 = math.sqrt(1 + (self.thickness/self.length)**2)
            self.surfacearea = self.width * (self.C1 *self.length + 
                                    self.length**2 / self.thickness * math.log(self.thickness/self.length + self.C1)) #surface area of fin
            self.mL = math.sqrt(2 * self.heattransfercoefficient / self.thermalconductivity / self.thickness) * self.length
            self.efficiency = 2 / (math.sqrt(4*self.mL**2 + 1) + 1)
            return self.efficiency, self.surfacearea
        except:
            print("Check your 'Fin()' object. You maynot have defined all required parameters.")
            print("\nstraightparabolic fin requires: length=someValue, width=SomeValue, thickness=someValue, heattransfercoefficient=someValue, thermalconductivity=someValue")
            

    def pintriangular(self):
        r""" Computes fin efficiency for a pintriangular fin.
        

        Parameters
        ----------
        `None_required` : 'None'
            Attributes that are already defined are used in calculation.
                
        
        Returns
        -------
        tuple containing efficiency and surfacearea :
            (efficiency : `int or float`, surfacearea : `int or float`)
        
                
        Examples
        --------
        First import the module **fins**.
           
        >>> from pychemengg.heattransfer import fins as fins
        >>> fin1 = fins.Fin(length=0.015, diameter=0.0015, heattransfercoefficient=50, thermalconductivity=70)
        fin1.pintriangular()
        >>> fin1.pintriangular()
        (0.9354407776522699, 3.538706842238283e-05)
        

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
    # def pintriangular(self, length=1, diameter=1, heattransfercoefficient=1, thermalconductivity=1):
        try:
            besseli = scipy.special.iv # renaming scipy function for clarity 
            self.surfacearea = math.pi * self.diameter / 2 * math.sqrt(self.length**2 + (self.diameter/2)**2) #surface area of fin
            self.mL = math.sqrt(4*self.heattransfercoefficient/self.thermalconductivity/self.diameter) * self.length
            self.efficiency = 2 * besseli(2, 2*self.mL) / (self.mL*besseli(1,2*self.mL));
            return self.efficiency, self.surfacearea    
        except:
            print("Check your 'Fin()' object. You maynot have defined all required parameters.")
            print("\npintriangular fin requires: length=someValue, diameter=SomeValue, heattransfercoefficient=someValue, thermalconductivity=someValue")
 
    
    def rectangularannular(self):
        r""" Computes fin efficiency for a rectangularannular fin.

        Parameters
        ----------
        `None_required` : 'None'
            Attributes that are already defined are used in calculation.
                
        
        Returns
        -------
        tuple containing efficiency and surfacearea :
            (efficiency : `int or float`, surfacearea : `int or float`)
        
             
        Examples
        --------
        First import the module **fins**.
           
        >>> from pychemengg.heattransfer import fins as fins
        >>> fin1 = fins.Fin(inner_radius=0.025, outer_radius=0.045, thickness=0.001, heattransfercoefficient=50, thermalconductivity=210)
        >>> fin1.rectangularannular()
        (0.922097299521443, 0.00879645943005142)
        

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
    # def rectangularannular(self, inner_radius=1, outer_radius=1, thickness=1, heattransfercoefficient=1, thermalconductivity=1):
        try:
            self.surfacearea = 2 * math.pi * (self.outer_radius**2 - self.inner_radius**2) #surface area of fin
            mro = math.sqrt(2 * self.heattransfercoefficient / self.thermalconductivity / self.thickness) * self.outer_radius
            mri = math.sqrt(2 * self.heattransfercoefficient / self.thermalconductivity / self.thickness) * self.inner_radius
            part1 = 2 * mri/(mro**2 - mri**2)
            besselk = scipy.special.kv # renaming scipy function for clarity
            besseli = scipy.special.iv # renaming scipy function for clarity
            part2 = besselk(1,mri)*besseli(1,mro) - besseli(1,mri)*besselk(1,mro)
            part3 = besseli(0,mri)*besselk(1,mro) + besselk(0,mri)*besseli(1,mro)
            self.efficiency = part1 * part2/part3
            return self.efficiency, self.surfacearea
        except:
            print("Check your 'Fin()' object. You maynot have defined all required parameters.")
            print("\nrectangularannular fin requires: inner_radius=someValue, outer_radius=someValue, thickness=someValue, heattransfercoefficient=someValue, thermalconductivity=someValue")