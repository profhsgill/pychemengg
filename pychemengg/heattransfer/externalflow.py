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

"""Module for external flow heat transfer

This module contains Nusselt number correlations for heat transfer between
a solid and a fluid subjected to convective flow over the external surfaces
of the object.


under convective flow on the surfaces of the objects.

The module defines classes for the following solid geometries:
    
* Plate
* Cylinder
* Sphere
* Bank of tubes
"""

import numpy as np
from scipy import interpolate
from scipy.interpolate import interp1d
import bisect
import warnings
import pychemengg.heattransfer.heatcommonmethods as hcm
warnings.simplefilter('always')

__all__ = ["Plate", "Cylinder", "Sphere", "TubeBank"]



class Plate:
    r""" Models external flow over flat plat.


    Parameters
    ----------
    `None_required` : 'None'
        This class takes no parameters for instance creation.
 
    
    Attributes
    ----------
    `None_required` : 'None'
        This class does not expose any instance attributes.
    

    Examples
    --------
    First import the module **externalflow**.
       
    >>> from pychemengg.heattransfer import externalflow as extflow 
    >>> plateobject = extflow.Plate()
    # This will create an instance of the class 'Plate'.
    # Methods of the class 'Plate' can then be called like so :-
    # plateobject.method(kwarg1=x, ...)

    """

    # ================
    # Isothermal cases
    # ================
    def Nu_laminar_local(self, Re=None, Pr=None):
        r""" Local Nusselt number for isothermal laminar flow on flat plate.
        
        
        Parameters
        ----------
        Re : `int or float`
            Reynolds number at location/length 'x' on plate.
        Pr : `int or float`
            Prandtl number for the fluid.
        

        Returns
        -------
        Nu : `int or float`
            Local Nusselt number at location/length 'x' on plate.
       
        
        Notes
        -----
        The following formula is used:         
        
        .. math::
            
            Nu = 0.332 Re^{0.5} Pr^{1/3}
            
               
        *where:*
            
            :math:`Pr > 0.6`
            
            :math:`Re < 5 * 10^5`
            
            Fluid properties are at film temp (:math:`T_{film}`):
            
                :math:`T_{film} = (T_{infinity} + T_{surface})/2`
                
                :math:`T_{infinity}` *= temperature of fluid away from surface*
                
                :math:`T_{surface}` *= temperature of surface*
        
            
        Warnings
        --------
        A Nusselt number is returned based on the equation
        even if parameters (such as *Re*, *Pr*) do not fall in their
        respective allowable range limits (see above under 'Notes').
        However, if this happens, a warning is issued.

        
        Examples
        --------
        First import the module **externalflow**.
           
        >>> from pychemengg.heattransfer import externalflow as extflow 
        >>> plateobject = extflow.Plate()
        # This will create an instance of the class 'Plate'.
        # Then call the method like so :-
        >>> plateobject.Nu_laminar_local(Re=4.024e4, Pr=2962)
        956.4495954202185

        
        References
        ----------
        [1] Yunus A. Cengel and Afshin J. Ghajar, "Heat And Mass Transfer
        Fundamentals and Applications", 6th Edition. New York, McGraw Hill
        Education, 2020.
        
        """ 
        Nu = 0.332 * np.power(Re, 0.5) * np.power(Pr, 1/3)
        
        if not (Re<5e5): warnings.warn("Re is not within range.")
        if not (Pr>0.6): warnings.warn("Pr is not within range.")
        return Nu


        
    
    def Nu_laminar_average(self, Re=None, Pr=None):
        r""" Average Nusselt number for isothermal laminar flow on flat plate.
        
        
        Parameters
        ----------
        Re : `int or float`
            Reynolds number for full length 'L' of plate.
        Pr : `int or float`
            Prandtl number for the fluid.
        

        Returns
        -------
        Nu : `int or float`
            Average Nusselt number for laminar flow on flat plate.
        
        
        Notes
        -----
        The following formula is used:         
        
        .. math::
            
            Nu = 0.664 Re^{0.5} Pr^{1/3}
            
               
        *where:*
            
            :math:`Pr > 0.6`
            
            :math:`Re < 5 * 10^5`
            
            Fluid properties are at film temp (:math:`T_{film}`):
            
                :math:`T_{film} = (T_{infinity} + T_{surface})/2`
                
                :math:`T_{infinity}` *= temperature of fluid away from surface*
                
                :math:`T_{surface}` *= temperature of surface*
            
        
        Warnings
        --------
        A Nusselt number is returned based on the equation
        even if parameters (such as *Re*, *Pr*) do not fall in their
        respective allowable range limits (see above under 'Notes').
        However, if this happens, a warning is issued.
        
        
        Examples
        --------
        First import the module **externalflow**.
           
        >>> from pychemengg.heattransfer import externalflow as extflow 
        >>> plateobject = extflow.Plate()
        # This will create an instance of the class 'Plate'.
        # Then call the method like so :-
        >>> plateobject.Nu_laminar_average(Re=4.024e4, Pr=2962)
        1912.899190840437

        
        References
        ----------
        [1] Yunus A. Cengel and Afshin J. Ghajar, "Heat And Mass Transfer
        Fundamentals and Applications", 6th Edition. New York, McGraw Hill
        Education, 2020.
        
        """

        Nu = 0.664 * np.power(Re, 0.5) * np.power(Pr, 1/3)

        if not (Re<5e5): warnings.warn("Re is not within range.")
        if not (Pr>0.6): warnings.warn("Pr is not within range.")
        return Nu

        
    def Nu_turbulent_local(self, Re=None, Pr=None):
        r"""Local Nusselt number for isothermal turbulent flow on flat plate.
        
        
        Parameters
        ----------
        Re : `int or float`
            Reynolds number at location/length 'x' on plate.
        Pr : `int or float`
            Prandtl number for the fluid.
        

        Returns
        -------
        Nu : `int or float`
            Local Nusselt number at location/length 'x' on plate.
        
        
        Notes
        -----
        The following formula is used:         
        
        .. math::
            
            Nu = 0.0296 Re^{0.8} Pr^{1/3}
            
               
        *where:*
            
            :math:`0.6 \eqslantless Pr \eqslantless 60`
            
            :math:`5*10^5 \eqslantless Re_x \eqslantless 10^7`
            
            Fluid properties are at film temp (:math:`T_{film}`):
            
                :math:`T_{film} = (T_{infinity} + T_{surface})/2`
                
                :math:`T_{infinity}` *= temperature of fluid away from surface*
                
                :math:`T_{surface}` *= temperature of surface*
            
        
        Warnings
        --------
        A Nusselt number is returned based on the equation
        even if parameters (such as *Re*, *Pr*) do not fall in their
        respective allowable range limits (see above under 'Notes').
        However, if this happens, a warning is issued.
        
        
        Examples
        --------
        First import the module **externalflow**.
           
        >>> from pychemengg.heattransfer import externalflow as extflow 
        >>> plateobject = extflow.Plate()
        # This will create an instance of the class 'Plate'.
        # Then call the method like so :-
        >>> plateobject.Nu_laminar_average(Re=4.024e4, Pr=2962)
        1912.899190840437

        
        References
        ----------
        [1] Yunus A. Cengel and Afshin J. Ghajar, "Heat And Mass Transfer
        Fundamentals and Applications", 6th Edition. New York, McGraw Hill
        Education, 2020.
        
        """
        Nu = 0.0296 * np.power(Re, 0.8) * np.power(Pr, 1/3)

        if not (5e5<=Re<=1e7): warnings.warn("Re is not within range.")
        if not (0.6<=Pr<=60): warnings.warn("Pr is not within range.")
        return Nu
    
    def Nu_turbulent_average(self, Re=None, Pr=None):
        r""" Average Nusselt number for isothermal turbulent flow on flat plate.
        
        
        Parameters
        ----------
        Re : `int or float`
            Reynolds number for full length 'L' of plate.
        Pr : `int or float`
            Prandtl number for the fluid.
        

        Returns
        -------
        Nu : `int or float`
            Average Nusselt number for turbulent flow on flat plate.
        
        
        Notes
        -----
        The following formula is used:         
        
        .. math::
            
            Nu = 0.037 Re^{0.8} Pr^{1/3}
            
               
        *where:*
            
            :math:`0.6 \eqslantless Pr \eqslantless 60`
            
            :math:`5*10^5 \eqslantless Re_x \eqslantless 10^7`
            
            Fluid properties are at film temp (:math:`T_{film}`):
            
                :math:`T_{film} = (T_{infinity} + T_{surface})/2`
                
                :math:`T_{infinity}` *= temperature of fluid away from surface*
                
                :math:`T_{surface}` *= temperature of surface*
            
        
        Warnings
        --------
        A Nusselt number is returned based on the equation
        even if parameters (such as *Re*, *Pr*) do not fall in their
        respective allowable range limits (see above under 'Notes').
        However, if this happens, a warning is issued.
        
        
        Examples
        --------
        First import the module **externalflow**.
           
        >>> from pychemengg.heattransfer import externalflow as extflow 
        >>> plateobject = extflow.Plate()
        # This will create an instance of the class 'Plate'.
        # Then call the method like so :-
        >>> plateobject.Nu_turbulent_average(Re=8.875e5, Pr=0.7387)
        1918.1936389098291

        
        References
        ----------
        [1] Yunus A. Cengel and Afshin J. Ghajar, "Heat And Mass Transfer
        Fundamentals and Applications", 6th Edition. New York, McGraw Hill
        Education, 2020.
        
        """

        Nu = 0.037 * np.power(Re, 0.8) * np.power(Pr, 1/3)
        
        if not (5e5<=Re<=1e7): warnings.warn("Re is not within range.")
        if not (0.6<=Pr<=60): warnings.warn("Pr is not within range.")
        return Nu
    
    def Nu_mixed_average(self, Re=None, Pr=None, Re_critical=5e5):
        r""" Average Nusselt number for isothermal 'laminar + turbulent' flow on flat plate.
        
        
        Parameters
        ----------
        Re : `int or float`
            Reynolds number for full length 'L' of plate
        Pr : `int or float` for the fluid.
            Prandtl number
        Re_critical : `int or float`
            Critical Reynolds number indicating transition between *laminar <--> turbulent*.
            Default value = :math:`5x10^5`.
        

        Returns
        -------
        Nu : `int or float`
            Average Nusselt number for 'laminar + turbulent' flow on flat plate.
        
        
        Notes
        -----
        The following formula is used:         
        
        .. math::
            
            Nu = (0.037 Re^{0.8} – (0.037 Re_{critical}^{4/5} – 0.664Re_{critical}^{1/2}))Pr^{1/3}
            
               
        *where:*
            
            :math:`0.6 \eqslantless Pr \eqslantless 60`
            
            :math:`5*10^5 \eqslantless Re_x \eqslantless 10^7`
            
            :math:`Re_{critical}` : is the number where fluid transitions from
            *laminar* to *turbulent*. This number can change from case to case, 
            although it is often = :math:`5*10^5`. For example, by adding
            roughness to plate surface, the transition to turbulent flow
            can occur at lower Reynolds numbers.
            
            Fluid properties are at film temp (:math:`T_{film}`):
            
                :math:`T_{film} = (T_{infinity} + T_{surface})/2`
                
                :math:`T_{infinity}` *= temperature of fluid away from surface*
                
                :math:`T_{surface}` *= temperature of surface*
            
        
        Warnings
        --------
        A Nusselt number is returned based on the equation
        even if parameters (such as *Re*, *Pr*) do not fall in their
        respective allowable range limits (see above under 'Notes').
        However, if this happens, a warning is issued.
        
        
        Examples
        --------
        First import the module **externalflow**.
           
        >>> from pychemengg.heattransfer import externalflow as extflow 
        >>> plateobject = extflow.Plate()
        # This will create an instance of the class 'Plate'.
        # Then call the method like so :-
        >>> plateobject.Nu_mixed_average(Re=610626, Pr=0.7073, Re_critical=5e5)
        625.480561991628

        
        References
        ----------
        [1] Theodore L. Bergman, Adrienne S. Lavine,
        Frank P. Incropera. David P. DeWitt,
        "Fundamentals of Heat and Mass Transfer", 8th Edition, Wiley.
        
        """

        A = 0.037*np.power(Re_critical, 4/5) - 0.664*np.power(Re_critical, 1/2)
        Nu = (0.037 * np.power(Re, 0.8) - A) * np.power(Pr, 1/3)

        if not (5e5<=Re<=1e7): warnings.warn("Re is not within range.")
        if not (0.6<=Pr<=60): warnings.warn("Pr is not within range.")
        return Nu
    
    # # =================================================
    # # Unheated starting section : Isothermal conditions
    # # =================================================   
    
    def Nu_unheated_laminar_local(self, Re=None, Pr=None, unheatedlength=None, xlocation_tofindNu=None):
        r""" Local Nusselt number for isothermal laminar flow on flat plate
        with unheated starting length.
        
        
        Parameters
        ----------
        Re : `int or float`
            Reynolds number at location/length 'x' on plate.
        Pr : `int or float`
            Prandtl number for the fluid.
        unheatedlength : `int or float`
            Initial length of plate that is not participating in heat transfer.
        xlocation_tofindNu : `int or float`
            x position to find local Nusselt number
        

        Returns
        -------
        Nu : `int or float`
            Local Nusselt number at location/length 'x' on plate.
 
        
        Notes
        -----
        The following formula is used:         
        
        .. math::
            
            Nu = \cfrac {Nu_{for\zeta=0}} {[1-(\zeta/x)^{3/4}]^{1/3}}
            
               
        *where:*
            
            :math:`Pr > 0.6`
            
            :math:`Re < 5 * 10^5`
            
            :math:`\zeta` = length of plate that is not heated
            
            x = position along length of plate where local Nusselt number is
            to be computed
            
            Fluid properties are at film temp (:math:`T_{film}`):
            
                :math:`T_{film} = (T_{infinity} + T_{surface})/2`
                
                :math:`T_{infinity}` *= temperature of fluid away from surface*
                
                :math:`T_{surface}` *= temperature of surface*
            
        
        Warnings
        --------
        A Nusselt number is returned based on the equation
        even if parameters (such as *Re*, *Pr*) do not fall in their
        respective allowable range limits (see above under 'Notes').
        However, if this happens, a warning is issued.
        
        
        Examples
        --------
        First import the module **externalflow**.
           
        >>> from pychemengg.heattransfer import externalflow as extflow 
        >>> plateobject = extflow.Plate()
        # This will create an instance of the class 'Plate'.
        # Then call the method like so :-
        >>> plateobject.Nu_unheated_laminar_local(Re=61062, Pr=0.7073, unheatedlength=0.1, xlocation_tofindNu=0.5)
        82.28737639413838

        
        References
        ----------
        [1] Yunus A. Cengel and Afshin J. Ghajar, "Heat And Mass Transfer
        Fundamentals and Applications", 6th Edition. New York, McGraw Hill
        Education, 2020.
        """
        num = Plate.Nu_laminar_local(Re=Re, Pr=Pr)
        denom = np.power((1-np.power(unheatedlength/xlocation_tofindNu, 3/4)),1/3)
        Nu = num/denom

        if not (Re<5e5): warnings.warn("Re is not within range.")
        if not (Pr>0.6): warnings.warn("Pr is not within range.")
        return Nu
    
    def Nu_unheated_laminar_average(self, Re=None, Pr=None, unheatedlength=None, length=None):
        r""" Average Nusselt number for isothermal laminar flow on flat plate
        with unheated starting length.
        
        
        Parameters
        ----------
        Re : `int or float`
            Reynolds number on plate (plate length to be used to compute Re).
        Pr : `int or float`
            Prandtl number for the fluid.
        unheatedlength : `int or float`
            Initial length of plate that is not participating in heat transfer.
        length : `int or float`
            Overall length of plate.
        

        Returns
        -------
        Nu : `int or float`
            Average Nusselt number over the plate.
 
        
        Notes
        -----
        The following formula is used:         
        
        .. math::
            
            Nu = \left( Nu_{for \zeta=0} \right) \left( \frac{L}{L-\zeta} \right)
            \left[1- \left(\frac{\zeta}{L} \right)^{3/4} \right]^{2/3}
                       
               
        *where:*
            
            :math:`Pr > 0.6`
            
            :math:`Re < 5 * 10^5`
            
            :math:`\zeta` = Length of plate that is not heated
            
            L = Overall length of plate over which average Nusselt number is
            to be computed
            
            Fluid properties are at film temp (:math:`T_{film}`):
            
                :math:`T_{film} = (T_{infinity} + T_{surface})/2`
                
                :math:`T_{infinity}` *= temperature of fluid away from surface*
                
                :math:`T_{surface}` *= temperature of surface*
            
        
        Warnings
        --------
        A Nusselt number is returned based on the equation
        even if parameters (such as *Re*, *Pr*) do not fall in their
        respective allowable range limits (see above under 'Notes').
        However, if this happens, a warning is issued.
        
        
        Examples
        --------
        First import the module **externalflow**.
           
        >>> from pychemengg.heattransfer import externalflow as extflow 
        >>> plateobject = extflow.Plate()
        # This will create an instance of the class 'Plate'.
        # Then call the method like so :-
        >>> plateobject.Nu_unheated_laminar_average(Re=61062, Pr=0.7073, unheatedlength=0.1, length=0.5)
        144.1942769849126

        
        References
        ----------
        [1] Yunus A. Cengel and Afshin J. Ghajar, "Heat And Mass Transfer
        Fundamentals and Applications", 6th Edition. New York, McGraw Hill
        Education, 2020.
        
        """
        num1 = Plate.Nu_laminar_average(Re=Re, Pr=Pr)
        num2 = np.power((1-np.power(unheatedlength/length, 3/4)),2/3)
        Nu = num1 * (length/(length-unheatedlength)) * num2

        if not (Re<5e5): warnings.warn("Re is not within range.")
        if not (Pr>0.6): warnings.warn("Pr is not within range.")
        return Nu
       
    def Nu_unheated_turbulent_local(self, Re=None, Pr=None, unheatedlength=None, xlocation_tofindNu=None):
        r""" Local Nusselt number for isothermal turbulent flow on flat plate
        with unheated starting length.
        
        
        Parameters
        ----------
        Re : `int or float`
            Reynolds number at location/length 'x' on plate.
        Pr : `int or float`
            Prandtl number for the fluid.
        unheatedlength : `int or float`
            Initial length of plate that is not participating in heat transfer.
        xlocation_tofindNu : `int or float`
            x position to find local Nusselt number
        

        Returns
        -------
        Nu : `int or float`
            Local Nusselt number at location/length 'x' on plate.
 
        
        Notes
        -----
        The following formula is used:         
        
        .. math::
            
            Nu = \cfrac {Nu_{for\zeta=0}} {[1-(\zeta/x)^{9/10}]^{1/9}}
            
               
        *where:*
            
            :math:`0.6 \eqslantless Pr \eqslantless 60`
            
            :math:`5*10^5 \eqslantless Re_x \eqslantless 10^7`
            
            :math:`\zeta` = length of plate that is not heated
            
            x = position along length of plate where local Nusselt number is
            to be computed
            
            Fluid properties are at film temp (:math:`T_{film}`):
            
                :math:`T_{film} = (T_{infinity} + T_{surface})/2`
                
                :math:`T_{infinity}` *= temperature of fluid away from surface*
                
                :math:`T_{surface}` *= temperature of surface*
            
        
        Warnings
        --------
        A Nusselt number is returned based on the equation
        even if parameters (such as *Re*, *Pr*) do not fall in their
        respective allowable range limits (see above under 'Notes').
        However, if this happens, a warning is issued.
        
        
        Examples
        --------
        First import the module **externalflow**.
           
        >>> from pychemengg.heattransfer import externalflow as extflow 
        >>> plateobject = extflow.Plate()
        # This will create an instance of the class 'Plate'.
        # Then call the method like so :-
        >>> plateobject.Nu_unheated_turbulent_local(Re=610620, Pr=0.7073, unheatedlength=0.1, xlocation_tofindNu=0.5)
        1155.3088951877935
        
        
        References
        ----------
        [1] Yunus A. Cengel and Afshin J. Ghajar, "Heat And Mass Transfer
        Fundamentals and Applications", 6th Edition. New York, McGraw Hill
        Education, 2020.
        
        """
        num = Plate.Nu_turbulent_local(Re=Re, Pr=Pr)
        denom = np.power((1-np.power(unheatedlength/xlocation_tofindNu, 9/10)),1/9)
        Nu = num/denom

        if not (5e5<=Re<=1e7): warnings.warn("Re is not within range.")
        if not (0.6<=Pr<=60): warnings.warn("Pr is not within range.")
        return Nu

    def Nu_unheated_turbulent_average(self, Re=None, Pr=None, unheatedlength=None, length=None):
        r""" Average Nusselt number for isothermal turbulent flow on flat plate
        with unheated starting length.
        
        
        Parameters
        ----------
        Re : `int or float`
            Reynolds number on plate (plate length to be used to compute Re).
        Pr : `int or float`
            Prandtl number for the fluid.
        unheatedlength : `int or float`
            Initial length of plate that is not participating in heat transfer.
        length : `int or float`
            Overall length of plate.
        

        Returns
        -------
        Nu : `int or float`
            Average Nusselt number over the plate.
 
        
        Notes
        -----
        The following formula is used:         
        
        .. math::
            
            Nu = \left( Nu_{for \zeta=0} \right) \left( \frac{L}{L-\zeta} \right)
            \left[1- \left(\frac{\zeta}{L} \right)^{9/10} \right]^{8/9}
            
               
        *where:*
            
            :math:`0.6 \eqslantless Pr \eqslantless 60`
            
            :math:`5*10^5 \eqslantless Re_x \eqslantless 10^7`
                
            :math:`\zeta` = Length of plate that is not heated
            
            L = Overall length of plate over which average Nusselt number is
            to be computed
            
            Fluid properties are at film temp (:math:`T_{film}`):
            
                :math:`T_{film} = (T_{infinity} + T_{surface})/2`
                
                :math:`T_{infinity}` *= temperature of fluid away from surface*
                
                :math:`T_{surface}` *= temperature of surface*
            
        
        Warnings
        --------
        A Nusselt number is returned based on the equation
        even if parameters (such as *Re*, *Pr*) do not fall in their
        respective allowable range limits (see above under 'Notes').
        However, if this happens, a warning is issued.
        
        
        Examples
        --------
        First import the module **externalflow**.
           
        >>> from pychemengg.heattransfer import externalflow as extflow 
        >>> plateobject = extflow.Plate()
        # This will create an instance of the class 'Plate'.
        # Then call the method like so :-
        >>> plateobject.Nu_unheated_turbulent_average(Re=610620, Pr=0.7073, unheatedlength=0.1, length=0.5)
        1381.0927382916545

        
        References
        ----------
        [1] Yunus A. Cengel and Afshin J. Ghajar, "Heat And Mass Transfer
        Fundamentals and Applications", 6th Edition. New York, McGraw Hill
        Education, 2020.
        """

        num1 = Plate.Nu_turbulent_average(Re=Re, Pr=Pr)
        num2 = np.power((1-np.power(unheatedlength/length, 9/10)),8/9)
        Nu = num1 * (length/(length-unheatedlength)) * num2

        if not (5e5<=Re<=1e7): warnings.warn("Re is not within range.")
        if not (0.6<=Pr<=60): warnings.warn("Pr is not within range.")
        return Nu

     
        
    # =============================
    # Uniform Heat Flux Conditions
    # =============================    
    
    def Nu_uniformflux_laminar_local(self, Re=None, Pr=None):
        r""" Local Nusselt number for laminar flow on flat plate with uniform heat flux.
        
        
        Parameters
        ----------
        Re : `int or float`
            Reynolds number at location/length 'x' on plate.
        Pr : `int or float` for the fluid.
            Prandtl number
        

        Returns
        -------
        Nu : `int or float`
            Local Nusselt number at location/length 'x' on plate.
 
        
        Notes
        -----
        The following formula is used:         
        
        .. math::
            
            Nu = 0.453 Re^{0.5} Pr^{1/3}
            
               
        *where:*
            
            :math:`Pr > 0.6`
            
            :math:`Re < 5 * 10^5`
            
            Fluid properties are at film temp (:math:`T_{film}`):
            
                :math:`T_{film} = (T_{infinity} + T_{surface})/2`
                
                :math:`T_{infinity}` *= temperature of fluid away from surface*
                
                :math:`T_{surface}` *= temperature of surface*
            
        
        Warnings
        --------
        A Nusselt number is returned based on the equation
        even if parameters (such as *Re*, *Pr*) do not fall in their
        respective allowable range limits (see above under 'Notes').
        However, if this happens, a warning is issued.
        
        
        Examples
        --------
        First import the module **externalflow**.
           
        >>> from pychemengg.heattransfer import externalflow as extflow 
        >>> plateobject = extflow.Plate()
        # This will create an instance of the class 'Plate'.
        # Then call the method like so :-
        >>> plateobject.Nu_uniformflux_laminar_local(Re=61062, Pr=0.7073)
        99.73592054582787

        
        References
        ----------
        [1] Yunus A. Cengel and Afshin J. Ghajar, "Heat And Mass Transfer
        Fundamentals and Applications", 6th Edition. New York, McGraw Hill
        Education, 2020
        
        """

        Nu = 0.453 * np.power(Re, 0.5) * np.power(Pr, 1/3)

        if not (Re<5e5): warnings.warn("Re is not within range.")
        if not (Pr>0.6): warnings.warn("Pr is not within range.")
        return Nu
    
    
    def Nu_uniformflux_turbulent_local(self, Re=None, Pr=None):
        r""" Local Nusselt number for turbulent flow on flat plate with uniform heat flux.
        
        
        Parameters
        ----------
        Re : `int or float`
            Reynolds number at location/length 'x' on plate.
        Pr : `int or float`
            Prandtl number for the fluid.
        

        Returns
        -------
        Nu : `int or float`
            Local Nusselt number at location/length 'x' on plate.
 
        
        Notes
        -----
        The following formula is used:         
        
        .. math::
            
            Nu = 0.0308 Re^{0.8} Pr^{1/3}
            
               
        *where:*
            
            :math:`0.6 \eqslantless Pr \eqslantless 60`
            
            :math:`5*10^5 \eqslantless Re_x \eqslantless 10^7`
            
            Fluid properties are at film temp (:math:`T_{film}`):
            
                :math:`T_{film} = (T_{infinity} + T_{surface})/2`
                
                :math:`T_{infinity}` *= temperature of fluid away from surface*
                
                :math:`T_{surface}` *= temperature of surface*
            
        
        Warnings
        --------
        A Nusselt number is returned based on the equation
        even if parameters (such as *Re*, *Pr*) do not fall in their
        respective allowable range limits (see above under 'Notes').
        However, if this happens, a warning is issued.
        
        
        Examples
        --------
        First import the module **externalflow**.
        >>> from pychemengg.heattransfer import externalflow as extflow 
        >>> plateobject = extflow.Plate()
        # This will create an instance of the class 'Plate'.
        # Then call the method like so :-
        >>> plateobject.Nu_uniformflux_turbulent_local(Re=610620, Pr=0.7073)
        1166.9047896712561
        
        
        References
        ----------
        [1] Yunus A. Cengel and Afshin J. Ghajar, "Heat And Mass Transfer
        Fundamentals and Applications", 6th Edition. New York, McGraw Hill
        Education, 2020.
        
        """
        Nu = 0.0308 * np.power(Re, 0.8) * np.power(Pr, 1/3)

        if not (5e5<=Re<=1e7): warnings.warn("Re is not within range.")
        if not (0.6<=Pr<=60): warnings.warn("Pr is not within range.")
        return Nu

    
        
    def Nu_uniformflux_laminar_average(self, Re=None, Pr=None):
        r""" Average Nusselt number for laminar flow on flat plate with uniform heat flux.
        
        
        Parameters
        ----------
        Re : `int or float`
            Reynolds number for full length 'L' of plate.
        Pr : `int or float`
            Prandtl number for the fluid.
        

        Returns
        -------
        Nu : `int or float`
            Average Nusselt number for laminar flow on flat plate.
        
        
        Notes
        -----
        The following formula is used:         
        
        .. math::
            
            Nu = 0.664 Re^{0.5} Pr^{1/3}
            
               
        *where:*
            
            :math:`Pr > 0.6`
            
            :math:`Re < 5 * 10^5`
            
            Fluid properties are at film temp (:math:`T_{film}`):
            
                :math:`T_{film} = (T_{infinity} + T_{surface})/2`
                
                :math:`T_{infinity}` *= temperature of fluid away from surface*
                
                :math:`T_{surface}` *= temperature of surface*
            
        
        Warnings
        --------
        A Nusselt number is returned based on the equation
        even if parameters (such as *Re*, *Pr*) do not fall in their
        respective allowable range limits (see above under 'Notes').
        However, if this happens, a warning is issued.
        
        
        Examples
        --------
        First import the module **externalflow**.
           
        >>> from pychemengg.heattransfer import externalflow as extflow 
        >>> plateobject = extflow.Plate()
        # This will create an instance of the class 'Plate'.
        # Then call the method like so :-
        >>> plateobject.Nu_uniformflux_laminar_average(Re=4.024e4, Pr=2962)
        1912.899190840437

        
        References
        ----------
        [1] Yunus A. Cengel and Afshin J. Ghajar, "Heat And Mass Transfer
        Fundamentals and Applications", 6th Edition. New York, McGraw Hill
        Education, 2020.
        
        """
        return Plate.Nu_laminar_average(Re=Re, Pr=Pr)
    
    
    def Nu_uniformflux_turbulent_average(self, Re=None, Pr=None):
        r""" Average Nusselt number for turbulent flow on flat plate with
        uniform heat flux.

        Parameters
        ----------
        Re : `int or float`
            Reynolds number for full length 'L' of plate.
        Pr : `int or float`
            Prandtl number for the fluid.
        

        Returns
        -------
        Nu : `int or float`
            Average Nusselt number for turbulent flow on flat plate.
        
        
        Notes
        -----
        The following formula is used:         
        
        .. math::
            
            Nu = 0.037 Re^{0.8} Pr^{1/3}
            
               
        *where:*
            
            :math:`0.6 \eqslantless Pr \eqslantless 60`
            
            :math:`5*10^5 \eqslantless Re_x \eqslantless 10^7`
            
            Fluid properties are at film temp (:math:`T_{film}`):
            
                :math:`T_{film} = (T_{infinity} + T_{surface})/2`
                
                :math:`T_{infinity}` *= temperature of fluid away from surface*
                
                :math:`T_{surface}` *= temperature of surface*
            
        
        Warnings
        --------
        A Nusselt number is returned based on the equation
        even if parameters (such as *Re*, *Pr*) do not fall in their
        respective allowable range limits (see above under 'Notes').
        However, if this happens, a warning is issued.
        
        
        Examples
        --------
        First import the module **externalflow**.
           
        >>> from pychemengg.heattransfer import externalflow as extflow 
        >>> plateobject = extflow.Plate()
        # This will assign the class 'Plate' to the variable 'plateobject'.
        # Then call the method like so :-
        >>> plateobject.Nu_uniformflux_turbulent_average(Re=8.875e5, Pr=0.7387)
        1918.1936389098291

        
        References
        ----------
        [1] Yunus A. Cengel and Afshin J. Ghajar, "Heat And Mass Transfer
        Fundamentals and Applications", 6th Edition. New York, McGraw Hill
        Education, 2020.
        
        """
        return Plate.Nu_turbulent_average(Re=Re, Pr=Pr)
    

class Cylinder:

    r""" Models a circular cylinder.
    
    
    Parameters
    ----------
    `None_required` : 'None'
        This class takes no parameters for instance creation.
    
    
    Attributes
    ----------
    `None_required` : 'None'
        This class does not expose any instance attributes.
    
    
    Examples
    --------
    First import the module **externalflow**.
       
    >>> from pychemengg.heattransfer import externalflow as extflow 
    >>> cyl = extflow.Cylinder
    # This will create an instance of the class 'Cylinder'.
    # Methods of the class 'Cylinder' can then be called like so :-
    # cyl.method(kwarg1=x, ...)
    
    """
    
    def Nu_churchill_bernstein(self, Re=None, Pr=None):
        r""" Average Nusselt number for flow over circular cylinder.

        Parameters
        ----------
        Re : `int or float`
            Reynolds number for fluid flow over the cylinder.
        Pr : `int or float`
            Prandtl number for the fluid.
        

        Returns
        -------
        Nu : `int or float`
            Average Nusselt number for turbulent flow on circular cylinder.
        
        
        Notes
        -----
        The following formula is used:         
        
        .. math::
            
            Nu = 0.3 + \cfrac {0.62 Re^{1/2} Pr^{1/3}} {\left[ {1 + (\frac{0.4}{Pr})^{2/3}} \right]^{1/4}}
            \left[ 1 + \left( \frac {Re} {282000} \right)^{5/8} \right]^{4/5}
            
               
        *where:*
            
            :math:`Re Pr > 0.2`
            
            Fluid properties are at film temp (:math:`T_{film}`):
            
                :math:`T_{film} = (T_{infinity} + T_{surface})/2`
                
                :math:`T_{infinity}` *= temperature of fluid away from surface*
                
                :math:`T_{surface}` *= temperature of surface*
            
        
        Warnings
        --------
        A Nusselt number is returned based on the equation
        even if parameters (such as *Re*, *Pr*) do not fall in their
        respective allowable range limits (see above under 'Notes').
        However, if this happens, a warning is issued.
        
        
        Examples
        --------
        First import the module **externalflow**.
           
        >>> from pychemengg.heattransfer import externalflow as extflow 
        >>> cyl = extflow.Cylinder()
        # This will create an instance of the class 'Cylinder'.
        # Then call the method like so :-
        >>> cyl.Nu_churchill_bernstein(Re=4.219e4, Pr=0.7202)
        124.44556447378223

        
        References
        ----------
        [1] Yunus A. Cengel and Afshin J. Ghajar, "Heat And Mass Transfer
        Fundamentals and Applications", 6th Edition. New York, McGraw Hill
        Education, 2020.
        
        """
        term1 = 0.62 * np.power(Re, 1/2) * np.power(Pr, 1/3)
        term2 = np.power(1+np.power(0.4/Pr,2/3),1/4)
        # if 20000< Re <400000:
        #     term3 = 1+np.power(Re/282000,1/2)
        # elif 100<Re<=20000 or 400000<= Re <1e7:
        #     term3 = np.power(1+np.power(Re/282000,5/8),4/5)
        term3 = np.power(1+np.power(Re/282000,5/8),4/5)
        Nu = 0.3 + term1/term2*term3

        if not (Re*Pr >= 0.2): warnings.warn("Product of Re and Pr is not within range.")
        return Nu
    
    def Nu_hilpert_knudsen_katz(self, Re=None, Pr=None):
        r""" Average Nusselt number for flow over circular cylinder.

        Parameters
        ----------
        Re : `int or float`
            Reynolds number for fluid flow over the cylinder.
        Pr : `int or float`
            Prandtl number for the fluid.
        

        Returns
        -------
        Nu : `int or float`
            Average Nusselt number for laminar/turbulent flow on circular cylinder.
        
        
        Notes
        -----
        The following formula is used:         
        
        .. math::
            
            Nu = C Re^{m} Pr^{1/3} 
            
               
        *where:*
            
            C, m : depend on Re
            
            :math:`0.4 \eqslantless Re \eqslantless 400000`
            
            :math:`Pr \eqslantgtr 0.7`
            
            Fluid properties are at film temp (:math:`T_{film}`):
            
                :math:`T_{film} = (T_{infinity} + T_{surface})/2`
                
                :math:`T_{infinity}` *= temperature of fluid away from surface*
                
                :math:`T_{surface}` *= temperature of surface*
            
        
        Warnings
        --------
        A Nusselt number is returned based on the equation
        even if parameters (such as *Re*, *Pr*) do not fall in their
        respective allowable range limits (see above under 'Notes').
        However, if this happens, a warning is issued.
        
        
        Examples
        --------
        First import the module **externalflow**.
           
        >>> from pychemengg.heattransfer import externalflow as extflow 
        >>> cyl = extflow.Cylinder()
        # This will create an instance of the class 'Cylinder'.
        # Then call the method like so :-
        >>> cyl.Nu_hilpert_knudsen_katz(Re=4.219e4, Pr=0.7202)
        127.97991118302727
        
        References
        ----------
        [1] Yunus A. Cengel and Afshin J. Ghajar, "Heat And Mass Transfer
        Fundamentals and Applications", 6th Edition. New York, McGraw Hill
        Education, 2020.
        
        """
        if Re <4: c, m = 0.989, 0.330
        if 4<= Re <40: c, m = 0.911, 0.385
        if 40<= Re <4000: c, m = 0.683, 0.466
        if 4000<= Re <40000: c, m = 0.193, 0.618
        if 40000<= Re : c, m = 0.027, 0.805
        Nu = c * np.power(Re,m) * np.power(Pr, 1/3)

        if not (Pr >=0.7): warnings.warn("Pr is not within range.")
        if not (0.4<=Re<=400000): warnings.warn("Re is not within range.")
        return Nu


class Sphere:
    r""" Models a sphere.
    
    
    Parameters
    ----------
    `None_required` : 'None'
        This class takes no parameters for instance creation.
    
    
    Attributes
    ----------
    `None_required` : 'None'
        This class does not expose any instance attributes .
    
    
    Examples
    --------
    First import the module **externalflow**.
       
    >>> from pychemengg.heattransfer import externalflow as extflow 
    >>> sphere = extflow.Sphere
    # This will create an instance of the class 'Sphere'.
    # Methods of the class 'Sphere' can then be called like so :-
    # sphere.method(kwarg1=x, ...)
    
    """
    
    def Nu_whitaker(self, Re=None, Pr=None, viscosity_surface=None, viscosity_infinity=None):
        r""" Average Nusselt number for flow over sphere.
        
        Parameters
        ----------
        Re : `int or float`
            Reynolds number for fluid flow over the sphere.
        Pr : `int or float`
            Prandtl number for the fluid.
        viscosity_surface : `int or float`
            Viscosity of fluid at surface temperature.
        viscosity_infinity : `int or float`
            Viscosity of fluid at fluid temperature far from surface.
        
        
        Returns
        -------
        Nu : `int or float`
            Average Nusselt number for turbulent flow on sphere.
        
        
        Notes
        -----
        The following formula is used:         
        
        .. math::
            
            Nu = 2 + [0.4 Re^{1/2} + 0.06 Re^{2/3}] Pr^{0.4} 
            \left( \cfrac {\mu_{infinity}} {\mu_{surface}} \right)^{1/4} 
            
               
        *where:*
            
            :math:`3.5 \eqslantless Re \eqslantless 8x10^4`
            
            :math:`0.7 \eqslantless Pr \eqslantless 380`
            
            :math:`1.0 \eqslantless \cfrac {\mu_{infinity}} {\mu_{surface}} \eqslantless 3.2`
            
            Fluid properties are at film temp (:math:`T_{film}`):
            
                :math:`T_{film} = (T_{infinity} + T_{surface})/2`
                
                :math:`T_{infinity}` *= temperature of fluid away from surface*
                
                :math:`T_{surface}` *= temperature of surface*
            
        
        Warnings
        --------
        A Nusselt number is returned based on the equation
        even if parameters (such as *Re*, *Pr*) do not fall in their
        respective allowable range limits (see above under 'Notes').
        However, if this happens, a warning is issued.
        
        
        Examples
        --------
        First import the module **externalflow**.
           
        >>> from pychemengg.heattransfer import externalflow as extflow 
        >>> sphere = extflow.Sphere()
        # This will create an instance of the class 'Sphere'.
        # Then call the method like so :-
        >>> sphere.Nu_whitaker(Re=3059, Pr=0.708, viscosity_surface=2.075e-5, viscosity_infinity=1.8462e-5)
        124.44556447378223

        
        References
        ----------
        [1] Yunus A. Cengel and Afshin J. Ghajar, "Heat And Mass Transfer
        Fundamentals and Applications", 6th Edition. New York, McGraw Hill
        Education, 2020.
        
        """
        a = 0.4*np.power(Re, 1/2)
        b = 0.06*np.power(Re,2/3)
        c = np.power(Pr, 0.4)
        d = np.power(viscosity_infinity/viscosity_surface, 1/4)
        Nu = 2 + (a+b)*c*d

        if not (5e5<=Re<=1e7): warnings.warn("Re is not within range.")
        if not (0.6<=Pr<=60): warnings.warn("Pr is not within range.")
        return Nu


class TubeBank:
    r""" Models a bank of tubes.


    Parameters
    ----------
    config : `str`
        Tube pitch arrangement ("inline" or "staggered"; Default:"inline").
    totalrows : `int or float`
        Number of rows. Row: A line of tubes in direction perpendicular to
        the flow of fluid over the tubes.
    tubes_per_row : `int or float`
        Number of tubes in each row.
    transverse_pitch : `int or float`
        Center-to-center distance between consecutive tubes within
        a row measured in a direction perpendicular to external flow of fluid.
    longitudanal_pitch : `int or float`
        Center-to-center distance between tubes in consecutive rows measured
        in a direction parallel to external flow of fluid .
    length : `int or float`
        Length of tubes (If not provided, assume = 1).
    outer_tubediameter : `int or float`
        Diameter of tubes.
    T_infinity : `int or float`
        Temperature of external fluid before it enters tube bank.
    velocity_infinity : `int or float`
        Velocity of external fluid  before it enters tube bank.
        
  
    
    Attributes
    ----------
    See "Parameters". Additional attributes are listed below.
    diagonal_pitch : `int or float`
        Diagonal distance between tubes in adjacent rows.
    viscosity : `int or float`
        Viscosity of external fluid at mean temperature of fluid.
    specificheat : `int or float`
        Specific heat of external fluid at mean temperature of fluid.
    thermalconductivity : `int or float`
        Thermal conductivity of external fluid at mean temperature of fluid.
    density_surface : `int or float`, optional
        Density of fluid at surface temperature of tubes. Optional if 
        Prandtl number at surface temperature will not be computed.
    viscosity_surface : `int or float`, optional
        Viscosity of fluid at surface temperature of tubes. Optional if 
        Prandtl number at surface temperature will not be computed.
    specificheat_surface : `int or float`, optional
        Specific heat of fluid at surface temperature of tubes. Optional if 
        Prandtl number at surface temperature will not be computed.
    thermalconductivity_surface : `int or float`, optional
        Thermal conductivity of fluid  at surface temperature of tubes. Optional if 
        Prandtl number at surface temperature will not be computed.
    T_out : `int or float`
        Temperature of external fluid at outlet of tube bank.
    T_surface : `int or float`
        Temperature of tube surface.
    Pr : `int or float`
        Prandtl number for the fluid.
    Re : `int or float`
        Reynolds number for the fluid.
    Nu : `int or float`
        Nusselt number for the fluid-tube bank heat exchange.
    maxvelocity : `int or float`
        Maximum velocity of external fluid in tube bank.
    Pr_surface : `int or float`
        Prandtl number at surface temperature of tubes.
    LMTD : `int or float`
        Log mean temperature difference between fluid "in"/"out" 
        and tube surface temperatures.
    pressuredrop : `int or float`
        Pressure drop when external fluid flows in tube bank.
  

    Examples
    --------
    First import the module **externalflow**.
    
    Units used in this example: SI system.
    
    However, any consistent units can be used.
    
    >>> from pychemengg.heattransfer import externalflow as extflow 
    >>> tubebundle = extflow.TubeBank(config="staggered", totalrows=7, tubes_per_row=8,
                             transverse_pitch=31.3e-3, longitudanal_pitch=34.3e-3,
                             length=1, outer_tubediameter=16.4e-3,
                             T_infinity=15, velocity_infinity=6)
    # This will create an instance of 'TubeBank'.
    """
    
    def __init__(self, config="inline",
                 totalrows=None,
                 tubes_per_row=None,
                 transverse_pitch=None,
                 longitudanal_pitch=None,
                 length=None,
                 outer_tubediameter=None,
                 T_infinity=None,
                 velocity_infinity=None):
        # assign user input attributes
        self.config=config
        self.totalrows = totalrows
        self.tubes_per_row = tubes_per_row
        self.transverse_pitch = transverse_pitch
        self.longitudanal_pitch = longitudanal_pitch
        self.length = length
        self.outer_tubediameter = outer_tubediameter
        self.T_infinity = T_infinity
        self.velocity_infinity = velocity_infinity
        # compute additional attributes based on user input
        self.diagonal_pitch = np.sqrt(self.longitudanal_pitch**2 + (self.transverse_pitch/2)**2)

    def set_fluid_properties(self, density=None,
                             viscosity=None,
                             specificheat=None,
                             thermalconductivity=None,
                             density_surface=None,
                             viscosity_surface=None,
                             specificheat_surface=None,
                             thermalconductivity_surface=None):
        r"""Use this to input fluid properties.
        
        
        Parameters
        ----------
        viscosity : `int or float`
            Viscosity of external fluid at mean temperature of fluid.
        specificheat : `int or float`
            Specific heat of external fluid at mean temperature of fluid.
        thermalconductivity : `int or float`
            Thermal conductivity of external fluid at mean temperature of fluid.
        density_surface : `int or float`, optional
            Density of fluid at surface temperature of tubes. Optional if 
            Prandtl number at surface temperature will not be computed.
        viscosity_surface : `int or float`, optional
            Viscosity of fluid at surface temperature of tubes. Optional if 
            Prandtl number at surface temperature will not be computed.
        specificheat_surface : `int or float`, optional
            Specific heat of fluid at surface temperature of tubes. Optional if 
            Prandtl number at surface temperature will not be computed.
        thermalconductivity_surface : `int or float`, optional
            Thermal conductivity of fluid  at surface temperature of tubes. Optional if 
            Prandtl number at surface temperature will not be computed.
                
        
        Returns
        -------
        `None` : 'None'
            This method simply assigns different fluid properties to keyword 
            attributes. It does not return anything.
        
        
        Notes
        -----
        All properties are at :math:`T_{mean} = (T_{infinity} + T_{out})/2`.
        Exception: some properties with *_surface* are at tube surface temperature
        :math:`T_{surface}`. 
        
        
        Examples
        --------
        First import the module **externalflow**.
        
        Units used in this example: SI system.
        
        However, any consistent units can be used.
        
        >>> from pychemengg.heattransfer import externalflow as extflow 
        >>> tubebundle = extflow.TubeBank(config="staggered", totalrows=7, tubes_per_row=8,
                                 transverse_pitch=31.3e-3, longitudanal_pitch=34.3e-3,
                                 length=1, outer_tubediameter=16.4e-3,
                                 T_infinity=15, velocity_infinity=6)
        # This will create an instance of 'TubeBank'.
        # Next set fluid properties like so:-
        >>> tubebundle.set_fluid_properties(thermalconductivity=0.027484,
                                         density=1.108152,
                                         viscosity=1.92152e-05,
                                         specificheat=1007.6)
        >>>
      
        """
        self.density = density
        self.viscosity = viscosity
        self.specificheat = specificheat
        self.thermalconductivity = thermalconductivity
        self.density_surface = density_surface
        self.viscosity_surface = viscosity_surface
        self.specificheat_surface = specificheat_surface
        self.thermalconductivity_surface = thermalconductivity_surface
    
    def set_T_out(self, T_out=None):
        r"""Use this to input external fluid outlet temperature from tube bank.
        
        
        Parameters
        ----------
        T_out : `int or float`
            Temperature of fluid at outlet of tube bank.
                
        
        Returns
        -------
        `None` : 'None'
            This method simply assigns fluid outlet temperature
            to keyword attribute. It does not return anything.

                
        Examples
        --------
        First import the module **externalflow**.
        
        Units used in this example: SI system.
        
        However, any consistent units can be used.
        
        >>> from pychemengg.heattransfer import externalflow as extflow 
        >>> tubebundle = extflow.TubeBank(config="staggered", totalrows=7, tubes_per_row=8,
                                 transverse_pitch=31.3e-3, longitudanal_pitch=34.3e-3,
                                 length=1, outer_tubediameter=16.4e-3,
                                 T_infinity=15, velocity_infinity=6)
        # This will create an instance of 'TubeBank'.
        # Next set fluid outlet temperatue like so:-
        >>> tubebundle.set_T_out(114.5)
        >>>
    
        """
        self.T_out = T_out
    
    def set_T_surface(self, T_surface=None):
        r"""Use this to input tube surface temperature.
        
        
        Parameters
        ----------
        T_surface : `int or float`
            Temperature of tube surface.
                
        
        Returns
        -------
        `None` : 'None'
            This method simply assign tube surface temperature
            to keyword attribute. It does not return anything.
        
                
        Examples
        --------
        First import the module **externalflow**.
        
        Units used in this example: SI system.
        
        However, any consistent units can be used.
        
        >>> from pychemengg.heattransfer import externalflow as extflow 
        >>> tubebundle = extflow.TubeBank(config="staggered", totalrows=7, tubes_per_row=8,
                                 transverse_pitch=31.3e-3, longitudanal_pitch=34.3e-3,
                                 length=1, outer_tubediameter=16.4e-3,
                                 T_infinity=15, velocity_infinity=6)
        # This will create an instance of 'TubeBank'.
        # Next set tube surface temperatue like so:-
        >>> tubebundle.set_T_surface(70)
        >>>
       
        """
        self.T_surface = T_surface
    
    def set_Pr(self, Pr):
        r"""Use this to input fluid Prandtl number at :math:`T_{mean}`.
        
        
        Parameters
        ----------
        Pr : `int or float`
            Prandtl number of external fluid.
                
        
        Returns
        -------
        `None` : 'None'
            This method simply assign Prandtl number
            to keyword attribute. It does not return anything.
        
        
        See Also
        --------
        calc_Pr :
            Prandtl number can alternatively 
            be computed using the method 'calc_Pr'.
        
                
        Examples
        --------
        First import the module **externalflow**.
        
        Units used in this example: SI system.
        
        However, any consistent units can be used.
        
        >>> from pychemengg.heattransfer import externalflow as extflow 
        >>> tubebundle = extflow.TubeBank(config="staggered", totalrows=7, tubes_per_row=8,
                                 transverse_pitch=31.3e-3, longitudanal_pitch=34.3e-3,
                                 length=1, outer_tubediameter=16.4e-3,
                                 T_infinity=15, velocity_infinity=6)
        # This will create an instance of 'TubeBank'.
        # Next set Prandtl number like so:-
        >>> tubebundle.set_Pr(0.705)
        >>>
      
        """
        self.Pr = Pr
    
    def set_Re(self, Re):
        r"""Use this to input Reynolds number.
        
        
        Parameters
        ----------
        Re : `int or float`
            Reynolds number of external fluid flow.
                
        
        Returns
        -------
        `None` : 'None'
            This method simply assign Reynolds number
            to keyword attribute. It does not return anything.
        
        
        See Also
        --------
        calc_Re : 
            Reynolds number can alternatively 
            be computed using the method 'calc_Re'.
        
                
        Examples
        --------
        First import the module **externalflow**.
        
        Units used in this example: SI system.
        
        However, any consistent units can be used.
        
        >>> from pychemengg.heattransfer import externalflow as extflow 
        >>> tubebundle = extflow.TubeBank(config="staggered", totalrows=7, tubes_per_row=8,
                                 transverse_pitch=31.3e-3, longitudanal_pitch=34.3e-3,
                                 length=1, outer_tubediameter=16.4e-3,
                                 T_infinity=15, velocity_infinity=6)
        # This will create an instance of 'TubeBank'.
        # Next set Reynolds number like so:-
        >>> tubebundle.set_Re(13943)
        >>>
       
        """
        self.Re = Re

    def set_maxvelocity(self, maxvelocity):
        r"""Use this to input maximum velocity of fluid in tube bank.
        
        
        Parameters
        ----------
        maxvelocity : `int or float`
            Maximum velocity of external fluid flow.
                
        
        Returns
        -------
        `None` : 'None'
            This method simply assign maximum velocity
            to keyword attribute. It does not return anything.
        
    
        See Also
        --------
        calc_maxvelocity : 
            Maximum velocity can alternatively 
            be computed using the method 'calc_maxvelocity'.
        
                
        Examples
        --------
        First import the module **externalflow**.
        
        Units used in this example: SI system.
        
        However, any consistent units can be used.
        
        >>> from pychemengg.heattransfer import externalflow as extflow 
        >>> tubebundle = extflow.TubeBank(config="staggered", totalrows=7, tubes_per_row=8,
                                 transverse_pitch=31.3e-3, longitudanal_pitch=34.3e-3,
                                 length=1, outer_tubediameter=16.4e-3,
                                 T_infinity=15, velocity_infinity=6)
        # This will create an instance of 'TubeBank'.
        # Next set Reynolds number like so:-
        >>> tubebundle.set_maxvelocity(12.6)
        >>>
        
        """
        self.maxvelocity = maxvelocity
        
    def set_Pr_surface(self, Pr_surface):
        r"""Use this to input fluid Prandtl number at :math:`T_{surface}`.
        
        
        Parameters
        ----------
        Pr : `int or float`
            Prandtl number of external fluid.
                
        
        Returns
        -------
        `None` : 'None'
            This method simply assign Prandtl number
            to keyword attribute. It does not return anything.
        
        
        See Also
        --------
        calc_Pr_surface : 
            Prandtl number at surface temperature can alternatively 
            be computed using the method 'calc_Pr_surface'.
        
                
        Examples
        --------
        First import the module **externalflow**.
        
        Units used in this example: SI system.
        
        However, any consistent units can be used.
        
        >>> from pychemengg.heattransfer import externalflow as extflow 
        >>> tubebundle = extflow.TubeBank(config="staggered", totalrows=7, tubes_per_row=8,
                                 transverse_pitch=31.3e-3, longitudanal_pitch=34.3e-3,
                                 length=1, outer_tubediameter=16.4e-3,
                                 T_infinity=15, velocity_infinity=6)
        # This will create an instance of 'TubeBank'.
        # Next set Prandtl number like so:-
        >>> tubebundle.set_Pr_surface(0.701)
        >>>
      
        """
        self.Pr_surface = Pr_surface
    
    def calc_Pr(self):
        r"""Use this to calculate fluid Prandtl number at :math:`T_{mean}`.
        
        
        Parameters
        ----------
        `None_required` : 'None'
            Attributes that are already defined are used in calculation.
                
        
        Returns
        -------
        Pr : `int or float`
            Prandtl number for the fluid.
        
        Notes
        -----
        The following formula is used:         
        
        .. math:
            Pr = \frac {\mu C_p} {k}
               
        *where:*
        
            :math:`\mu` = viscosity of fluid
            
            :math:`C_p` = specific heat of fluid
            
            k =  thermal conductivity
            
            Pr = Prandtl number
            
            Fluid properties are at mean fluid temp (:math:`T_{mean}`):
            
                :math:`T_{mean} = (T_{infinity} + T_{out})/2`
                
                :math:`T_{infinity}` *= inlet temperature of fluid before it enters tube bank*
                
                :math:`T_{surface}` *= outlet temperature of fluid as it exits tube bank*
        
        
        See Also
        --------
        set_Pr : 
            If Prandtl number is known it can alternatively be set using the
            method 'set_Pr'.
        
                
        Examples
        --------
        First import the module **externalflow**.
        
        Units used in this example: SI system.
        
        However, any consistent units can be used.
        
        >>> from pychemengg.heattransfer import externalflow as extflow
        >>> tubebundle = extflow.TubeBank(config="staggered", totalrows=7, tubes_per_row=8,
                                 transverse_pitch=31.3e-3, longitudanal_pitch=34.3e-3,
                                 length=1, outer_tubediameter=16.4e-3,
                                 T_infinity=15, velocity_infinity=6)
        >>> tubebundle.set_fluid_properties(thermalconductivity=0.027484,
                                            density=1.108152,
                                            viscosity=1.92152e-05,
                                            specificheat=1007.6)
        >>> tubebundle.calc_Pr()
        0.7044547926066074
        >>>
        
        References
        ----------
        [1] Yunus A. Cengel and Afshin J. Ghajar, "Heat And Mass Transfer
        Fundamentals and Applications", 6th Edition. New York, McGraw Hill
        Education, 2020.
        """
        Pr = self.viscosity * self.specificheat / self.thermalconductivity
        self.Pr = Pr
        return self.Pr
    
    def calc_Pr_surface(self):
        r"""Use this to calculate fluid Prandtl number at :math:`T_{surface}`.
        
        
        Parameters
        ----------
        `None_required` : 'None'
            Attributes that are already defined are used in calculation.
                
        
        Returns
        -------
        Pr : `int or float`
            Prandtl number for the fluid.
        
        Notes
        -----
        The following formula is used:         
        
        .. math:
            Pr = \frac {\mu C_p} {k}
               
        *where:*
        
            :math:`\mu` = viscosity of fluid
            
            :math:`C_p` = specific heat of fluid
            
            k =  thermal conductivity
            
            Pr = Prandtl number
            
            Fluid properties are at surface temperature of tubes (:math:`T_{surface}`):

        
        See Also
        --------
        set_Pr_surface : 
            If Prandtl number at surface temperature is known it can
            alternatively be set using the method 'set_Pr_surface'.
        
                
        Examples
        --------
        First import the module **externalflow**.
        
        Units used in this example: SI system.
        
        However, any consistent units can be used.
        
        >>> from pychemengg.heattransfer import externalflow as extflow
        >>> tubebundle = extflow.TubeBank(config="staggered", totalrows=7, tubes_per_row=8,
                                 transverse_pitch=31.3e-3, longitudanal_pitch=34.3e-3,
                                 length=1, outer_tubediameter=16.4e-3,
                                 T_infinity=15, velocity_infinity=6)
        >>> tubebundle.set_fluid_properties(thermalconductivity_surface=0.029482,
                                            density_surface=1.018296,
                                            viscosity_surface=2.04896e-05,
                                            specificheat_surface=1008.7)
        >>> tubebundle.calc_Pr_surface()
        0.7010331565022726
        >>>

        References
        ----------
        [1] Yunus A. Cengel and Afshin J. Ghajar, "Heat And Mass Transfer
        Fundamentals and Applications", 6th Edition. New York, McGraw Hill
        Education, 2020.
        """
        Pr_surface = self.viscosity_surface * self.specificheat_surface / self.thermalconductivity_surface
        self.Pr_surface = Pr_surface
        return self.Pr_surface
    
    def calc_maxvelocity(self):
        r"""Use this to calculate maximum velocity of fluid in tube bank.
        
        Parameters
        ----------
        `None_required` : 'None'
            Attributes that are already defined are used in calculation.
                
        
        Returns
        -------
        maxvelocity : `int or float`
            Maximum velocity of the fluid in the tube bank.
        
        Notes
        -----
        The following formula are used:         
        
        if config == "inline":
        
        .. math::
            
            V_{max} = \left( \frac {S_T} {S_T - D} \right) V_\infty
        
        if config == "staggered":
            
            the greater of the following two:
            
         .. math:: 
             
            V_{max} = \left( \frac {S_T} {S_T - D} \right) V_\infty\\[15pt]
            
            
            |or|\\[15pt]
            

            V_{max} = \left( \frac {S_T} {2(S_D - D)} \right) V_\infty
               
        *where:*
        
            :math:`S_T` = transverse pitch
            
            :math:`S_D` = diagonal pitch
            
            :math:`D` = outside diameter of tubes
            
            :math:`V_\infty` = velocity of fluid at inlet of tube bank
            
            :math:`V_{max}` = maximum velocity of fluid in tube bank

        
        See Also
        --------
        set_maxvelocity : 
            If maximum velocity is known it can
            alternatively be set using the method 'set_maxvelocity'.
        
                
        Examples
        --------
        First import the module **externalflow**.
        
        Units used in this example: SI system.
        
        However, any consistent units can be used.
        
        >>> from pychemengg.heattransfer import externalflow as extflow
        >>> tubebundle = extflow.TubeBank(config="staggered", totalrows=7, tubes_per_row=8,
                                 transverse_pitch=31.3e-3, longitudanal_pitch=34.3e-3,
                                 length=1, outer_tubediameter=16.4e-3,
                                 T_infinity=15, velocity_infinity=6)
        >>> tubebundle.calc_maxvelocity()
        >>> 12.604026845637584
        >>>

        References
        ----------
        [1] Yunus A. Cengel and Afshin J. Ghajar, "Heat And Mass Transfer
        Fundamentals and Applications", 6th Edition. New York, McGraw Hill
        Education, 2020.
       """
        if self.config == "inline":
            self.maxvelocity = ((self.transverse_pitch)
                                / (self.transverse_pitch - self.outer_tubediameter)
                                * (self.velocity_infinity))
            return self.maxvelocity
        
        if self.config == "staggered":
            term1 = 2*(self.diagonal_pitch - self.outer_tubediameter)
            if term1 < (self.transverse_pitch - self.outer_tubediameter):
                self.maxvelocity = self.transverse_pitch / term1 * self.velocity_infinity
            if term1 > (self.transverse_pitch - self.outer_tubediameter):
                self.maxvelocity = ((self.transverse_pitch)
                                    / (self.transverse_pitch - self.outer_tubediameter)
                                    * (self.velocity_infinity))
            return self.maxvelocity
    
    
    def calc_Re(self):
        r"""Use this to calculate Reynolds number for external fluid flow in the tube bank.
        
        Parameters
        ----------
        `None_required` : 'None'
            Attributes that are already defined are used in calculation.
                
        
        Returns
        -------
        Re : `int or float`
            Reynold number for external fluid flow in the tube bank.
        
        Notes
        -----
        The following formula is used:         
        
        .. math::
            
            Re = \frac {D V_{max} \rho} {\mu}

        *where:*
            
            :math:`D` = outside diameter of tubes
            
            :math:`V_{max}` = maximum velocity of fluid in tube bank
            
            :math:`\rho` = fluid density
            
            :math:`\mu` = fluid viscosity
            
            Re = Reynolds number
            
            Fluid properties are at mean fluid temp (:math:`T_{mean}`):
            
                :math:`T_{mean} = (T_{infinity} + T_{out})/2`
                
                :math:`T_{infinity}` *= inlet temperature of fluid before it enters tube bank*
                
                :math:`T_{surface}` *= outlet temperature of fluid as it exits tube bank*

        
        See Also
        --------
        set_Re : 
            If Reynolds number is known it can
            alternatively be set using the method 'set_Re'.
        
                
        Examples
        --------
        First import the module **externalflow**.
        
        Units used in this example: SI system.
        
        However, any consistent units can be used.
        
        >>> from pychemengg.heattransfer import externalflow as extflow
        >>> tubebundle = extflow.TubeBank(config="staggered", totalrows=7, tubes_per_row=8,
                                 transverse_pitch=31.3e-3, longitudanal_pitch=34.3e-3,
                                 length=1, outer_tubediameter=16.4e-3,
                                 T_infinity=15, velocity_infinity=6)
        >>> tubebundle.set_fluid_properties(thermalconductivity=0.027484,
                                            density=1.108152,
                                            viscosity=1.92152e-05,
                                            specificheat=1007.6)
        >>> tubebundle.calc_Re()
        11920.860149026319


        References
        ----------
        [1] Yunus A. Cengel and Afshin J. Ghajar, "Heat And Mass Transfer
        Fundamentals and Applications", 6th Edition. New York, McGraw Hill
        Education, 2020.
       """
        self.Re = self.outer_tubediameter * self.maxvelocity * self.density / self.viscosity
        return self.Re
    
    
    def calc_Nu(self):
        r"""Use this to calculate Nusselt number for external fluid flow in the tube bank.
        
        Parameters
        ----------
        `None_required` : 'None'
            Attributes that are already defined are used in calculation.
                
        
        Returns
        -------
        Re : `int or float`
            Nusselt number for external fluid flow in the tube bank.
        
        
        Notes
        -----
        The following correlation by Zukauskas ([1]), and reproduced in 
        [2] was used.
        
        .. math::
            
            Nu = F C Re^{m} Pr^{n}  {\left(\frac{Pr}{Pr_{surface}}\right)}^{0.25}
            
               
        *where:*
            
            C, m, n : depend on Re and tube arrangement ("inline" vs "staggered")
            
            F : correction factor :math:`\eqslantless 1`
            *"F"* depends on the number of tube rows.
            If tube rows < 16, F < 1.0. 
            A look up table can be used to compute and interpolate *F* as
            a function of tube rows,  however, the values are only valid
            for Re > 1000.
            
            :math:`0.7 \eqslantless Pr \eqslantless 500`
            
            :math:`0 \eqslantless Re_x \eqslantless 2*10^6`
            
            Nu = Nusselt number
            
            Fluid properties are at mean fluid temp (:math:`T_{mean}`):
            
                :math:`T_{mean} = (T_{infinity} + T_{out})/2`
                
                :math:`T_{infinity}` *= inlet temperature of fluid before it enters tube bank*
                
                :math:`T_{out}` *= outlet temperature of fluid as it exits tube bank*
        
        
        Warnings
        --------
        A Nusselt number is returned based on the equation
        even if parameters (such as *Re*, *Pr*) do not fall in their
        respective allowable range limits (see above under 'Notes').
        However, if this happens, a warning is issued.
        
        
        Examples
        --------
        First import the module **externalflow**.
        
        Units used in this example: SI system.
        
        However, any consistent units can be used.
        
        >>> from pychemengg.heattransfer import externalflow as extflow
        >>> tubebundle = extflow.TubeBank(config="staggered", totalrows=7, tubes_per_row=8,
                                         transverse_pitch=31.3e-3, longitudanal_pitch=34.3e-3,
                                         length=1, outer_tubediameter=16.4e-3,
                                         T_infinity=15, velocity_infinity=6)
        >>> tubebundle.set_fluid_properties(thermalconductivity=0.027484, 
                                        density=1.108152, 
                                        viscosity=1.92152e-05, 
                                        specificheat=1007.6, 
                                        thermalconductivity_surface=0.029482,
                                        density_surface=1.018296,
                                        viscosity_surface=2.04896e-05,
                                        specificheat_surface=1008.7)
        >>> tubebundle.calc_maxvelocity()
        12.604026845637584
        >>> tubebundle.calc_Re()
        11920.860149026319
        >>> tubebundle.calc_Pr()
        0.7044547926066074
        >>> tubebundle.calc_Pr_surface()
        0.7010331565022726
        >>> tubebundle.calc_Nu()
        81.26988497431279
        
        References
        ----------
        [1] A. Zukauskas, "Heat Transfer From Tubes in Cross Flow." 
        In "Handbook of Single Phase Convective Heat Transfer."
        S. Kakac, R.K. Shah, and W. Aung (Eds.), 
        New York: Wiley Interscienc, 1987.
        
        [2] Yunus A. Cengel and Afshin J. Ghajar, "Heat And Mass Transfer
        Fundamentals and Applications", 6th Edition. New York, McGraw Hill 
        Education, 2020.
        
        """
        Re = self.Re
        Pr = self.Pr
        Pr_surface = self.Pr_surface
        
        if self.config == "inline":
            correction_factor = [0.70, 0.80, 0.86, 0.90, 0.93, 0.96, 0.98, 0.99, 1.0]
            tube_rownumber = [1, 2, 3, 4, 5, 7, 10, 13, 16]
            compute_C2 = interpolate.interp1d(tube_rownumber, correction_factor)
            if self.totalrows <= 16:
                C2 = compute_C2(self.totalrows)         
            if 0<= Re <1e2:
                Nu = 0.9 * Re**0.4 * Pr**0.36 * (Pr/Pr_surface)**0.25 * C2
            elif 1e2<= Re <1e3: 
                Nu = 0.52 * Re**0.5 * Pr**0.36 * (Pr/Pr_surface)**0.25 * C2
            elif 1e3<= Re <2e5:
                Nu = 0.27 * Re**0.63 * Pr**0.36 * (Pr/Pr_surface)**0.25 * C2
            elif 2e5<= Re <2e6:
                Nu = 0.333 * Re**0.8 * Pr**0.4 * (Pr/Pr_surface)**0.25 * C2
            
        if self.config == "staggered":
            correction_factor = [0.64, 0.76, 0.84, 0.89, 0.93, 0.96, 0.98, 0.99, 1.0]
            tube_rownumber = [1, 2, 3, 4, 5, 7, 10, 13, 16]
            compute_C2 = interpolate.interp1d(tube_rownumber, correction_factor)
            if self.totalrows<=16:
                C2 = compute_C2(self.totalrows)               
            if 0<= Re <500:
                Nu = 1.04 * Re**0.4 * Pr**0.36 * (Pr/Pr_surface)**0.25 * C2
            elif 500<= Re <1e3: 
                Nu = 0.71 * Re**0.5 * Pr**0.36 * (Pr/Pr_surface)**0.25 * C2
            elif 1e3<= Re <2e5:
                Nu = (0.35 * (self.transverse_pitch/self.longitudanal_pitch)**0.2
                        * Re**0.6 * Pr**0.36 * (Pr/Pr_surface)**0.25 * C2)
            elif 2e5<= Re <2e6:
                Nu = (0.031 * (self.transverse_pitch/self.longitudanal_pitch)**0.2
                      * Re**0.8 * Pr**0.36 * (Pr/Pr_surface)**0.25 * C2)
        self.Nu = Nu
        if not self.Re  > 1e3: warnings.warn("Re is not within range for use of correction factor data.")
        if not self.Re < 2e6: warnings.warn("Re is not within range for Nusselt number correlations.")
        if not 0.7<=Pr<=500: warnings.warn("Pr is not within range for Nusselt number correlations.")
        return self.Nu
    
    def calc_LMTD(self):
        r"""Use this to calculate log mean temperature difference (LMTD).
        
        Parameters
        ----------
        `None_required` : 'None'
            Attributes that are already defined are used in calculation.
                
        
        Returns
        -------
        LMTD : `int or float`
            Log mean temperature difference between :math:`T_{infinity}`,
            :math:`T_{out}`, and :math:`T_{surface}`.
        
        
        Notes
        -----
        The following formula is used.
        
        .. math::
            
            LMTD = \cfrac {(T_{surface} - T_{infinity}) - (T_{surface} - T_{out})}
            {ln[(T_{surface} - T_{infinity})/(T_{surface} - T_{out})]}
            
               
        *where:*
            
            :math:`T_{infinity}` *= inlet temperature of fluid before it enters tube bank*
            
            :math:`T_{out}` *= outlet temperature of fluid as it exits tube bank*
            
            :math:`T_{surface}` *= surface temperature of tubes in tube bank*
        
        
        See Also
        --------
        heatcommonmethods.LMTD : 
            This function is used to compute LMTD for the tube bank.
        
        
        Examples
        --------
        First import the module **externalflow**.
        
        Units used in this example: SI system.
        
        However, any consistent units can be used.
        
        >>> from pychemengg.heattransfer import externalflow as extflow
        >>> tubebundle = extflow.TubeBank(config="staggered", totalrows=7, tubes_per_row=8,
                                         transverse_pitch=31.3e-3, longitudanal_pitch=34.3e-3,
                                         length=1, outer_tubediameter=16.4e-3,
                                         T_infinity=15, velocity_infinity=6)
        >>> tubebundle.set_T_surface(70)
        >>> tubebundle.set_T_out(70-44.5)
        >>> tubebundle.calc_LMTD()
        49.56477500080964
        # Temperatures are in Celsius
        
        
        References
        ----------        
        [1] Yunus A. Cengel and Afshin J. Ghajar, "Heat And Mass Transfer
        Fundamentals and Applications", 6th Edition. New York, McGraw Hill 
        Education, 2020.
        
        """
        deltaT1 = self.T_surface-self.T_infinity
        deltaT2 = self.T_surface-self.T_out
        self.LMTD = hcm.calc_LMTD(deltaT1=deltaT1, deltaT2=deltaT2)
        return self.LMTD
    
    def calc_pressuredrop(self):
        r"""Use this to calculate pressure drop for fluid flow in the tube bank.
        
        Parameters
        ----------
        `None_required` : 'None'
            Attributes that are already defined are used in calculation.
                
        
        Returns
        -------
        pressure drop : `int or float`
            Pressure drop for external fluid flow in the tube bank.
        
        
        Notes
        -----
        The following correlation by Zukauskas and Ulinskas([1]) was used.
        This correlation is used in textbooks as a graphical plot.
        For example in [2].
        
        .. math::
            
           \Delta P = Eu \cfrac{\rho V_{max}^2}{2} N_L
            
               
        *where:*
            
            :math:`\Delta P` = Pressure drop
            
            Eu = Euler number
            
            :math:`\rho` = density of fluid
            
            :math:`V_{max}` = maximum fluid velocity in tube bank
            
            :math:`N_L` = number of tube rows
            
            Fluid properties are at mean fluid temp (:math:`T_{mean}`):
            
                :math:`T_{mean} = (T_{infinity} + T_{out})/2`
                
                :math:`T_{infinity}` *= inlet temperature of fluid before it enters tube bank*
                
                :math:`T_{out}` *= outlet temperature of fluid as it exits tube bank*
        
        
        Warnings
        --------
        A Nusselt number is returned based on the equation
        even if parameters (such as *Re*, *Pr*) do not fall in their
        respective allowable range limits (see above under 'Notes').
        However, if this happens, a warning is issued.
        
        
        Examples
        --------
        First import the module **externalflow**.
        
        Units used in this example: SI system.
        
        However, any consistent units can be used.
        
        >>> from pychemengg.heattransfer import externalflow as extflow
        >>> tubebundle = extflow.TubeBank(config="staggered", totalrows=7, tubes_per_row=8,
                                         transverse_pitch=31.3e-3, longitudanal_pitch=34.3e-3,
                                         length=1, outer_tubediameter=16.4e-3,
                                         T_infinity=15, velocity_infinity=6)
        >>> tubebundle.set_fluid_properties(thermalconductivity=0.027484, 
                                        density=1.108152, 
                                        viscosity=1.92152e-05, 
                                        specificheat=1007.6, 
                                        thermalconductivity_surface=0.029482,
                                        density_surface=1.018296,
                                        viscosity_surface=2.04896e-05,
                                        specificheat_surface=1008.7)
        >>> tubebundle.calc_maxvelocity()
        12.604026845637584
        >>> tubebundle.calc_Re()
        11920.860149026319
        >>> tubebundle.calc_Pr()
        0.7044547926066074
        >>> tubebundle.calc_Pr_surface()
        0.7010331565022726
        >>> tubebundle.calc_Nu()
        81.26988497431279
        >>> tubebundle.calc_pressuredrop()
        211.3578929319506
        
        References
        ----------
        [1] A. Zukauskas, "Heat Transfer From Tubes in Cross Flow." 
        In "Handbook of Single Phase Convective Heat Transfer."
        S. Kakac, R.K. Shah, and W. Aung (Eds.), 
        New York: Wiley Interscienc, 1987.
        
        [2] Yunus A. Cengel and Afshin J. Ghajar, "Heat And Mass Transfer
        Fundamentals and Applications", 6th Edition. New York, McGraw Hill 
        Education, 2020.
        
        """
        Eu_k1_ratio, k1 = self._calc_Eu_k1_values()
        frictionfactor = Eu_k1_ratio
        correctionfactor = k1
        self.pressuredrop = (self.totalrows
                             * frictionfactor
                             * correctionfactor
                             * self.density
                             * self.maxvelocity**2
                             / 2.0)
        return self.pressuredrop
        
    def _calc_Eu_k1_values(self): #Eu: Euler number
        if self.config == "inline":
            pitch_dia_ratio = self.longitudanal_pitch/self.outer_tubediameter
            pitch_dia_ratio_limits = [1.25, 1.5, 2]
        if self.config=="staggered":
            pitch_dia_ratio = self.transverse_pitch/self.outer_tubediameter
            pitch_dia_ratio_limits = [1.25, 1.5, 2, 2.5]
        
        # Step 1) calculate Eu_to_k1 ratio at interval limits
        
        # first find interval limits for pitch to dia ratio of tube bank
        # for example for "inline" if pitch_dia_ratio = 1.34
        # then looking at 'inline'  [1.25, 1.5, 2] limits
        # 1.34 will lie between 1.25 and 1.5
        # so we first need to find the interval limits that bound the given 
        # pitch_dia_ratio of the tube bank
        # print("pitch dia ratio", pitch_dia_ratio)
        leftlimit, rightlimit, outofrange = TubeBank._get_interval_limits(pitch_dia_ratio_limits, pitch_dia_ratio) 
        # print("left limit=", leftlimit, "right limit=", rightlimit)
        
        # next find Eu_k1_ratio at the left and right interval limits
        left_Eu_k1_ratio = self._get_Eu_k1_ratio( Re=self.Re, pitch_dia_ratio=pitch_dia_ratio_limits[leftlimit] )
        right_Eu_k1_ratio = self._get_Eu_k1_ratio( Re=self.Re, pitch_dia_ratio=pitch_dia_ratio_limits[rightlimit] )
        # print("left value=", pitch_dia_ratio_limits[leftlimit], "right value=", pitch_dia_ratio_limits[rightlimit] )
        
        # next create a function to interpolate 
        # x: interval limits
        # y: Eu_k1_ratio
        x = [pitch_dia_ratio_limits[leftlimit], pitch_dia_ratio_limits[rightlimit]]
        y = [left_Eu_k1_ratio, right_Eu_k1_ratio]
        func = interp1d(x,y)

        # next use the function and interpolate to get Eu_k1_ratio at given pitch_dia_ratio
        if not outofrange:
            Eu_k1_ratio = func(pitch_dia_ratio)
        if outofrange:
            Eu_k1_ratio = left_Eu_k1_ratio
            warnings.warn("Pitch-to-diameter ratio is out of interpolation range, so pressure drop computed at a boundary")
        # print("Eu_k1_ratio =", Eu_k1_ratio)
        
        # Step 2) get k1 geometry parameter
        # k1 is available for specific Reynolds numbers and must be interpolated
        # for the Re relevant to a given problem
        #
        if self.config == "inline":
            Re_with_k1_formulas = [1e3, 1e4, 1e5, 1e6, 1e7]
        if self.config=="staggered":
            Re_with_k1_formulas = [1e2, 1e3, 1e4, 1e5, 1e6]

        # first find interval limits
        # for example for "inline" if Re = 12234
        # then looking at 'inline'  [1e3, 1e4, 1e5, 1e6, 1e7] case
        # 12234 will lie between 1e4 and 1e5
        # so we first need to find the interval limits that bound the given 
        # Re of the tube bank
        leftlimit, rightlimit, outofrange = TubeBank._get_interval_limits(Re_with_k1_formulas, self.Re) 
        # print("left limit=", leftlimit, "right limit=", rightlimit)
        
        # next find k1 at these interval limits
        left_k1_value = self._get_k1_parameter( Re=Re_with_k1_formulas[leftlimit] )
        right_k1_value = self._get_k1_parameter( Re=Re_with_k1_formulas[rightlimit] )
        # print("left limit = ", Re_with_k1_formulas[leftlimit], "right limit=", Re_with_k1_formulas[rightlimit]  )
        # print("left k1 =", left_k1_value, "right k1", right_k1_value)

        # # next create a function to interpolate
        # # create a log-log function 
        # x: Re
        # y: k1 value for Re
        x = np.log10( [Re_with_k1_formulas[leftlimit], Re_with_k1_formulas[rightlimit]] )
        y = np.log10( [left_k1_value, right_k1_value] ) 
        func = interp1d( x,y )

        # use function and interpolate to get Eu_k1_ratio at given pitch_dia_ratio
        if not outofrange:
            k1 = np.power(10, func(np.log10(self.Re)))
        if outofrange:
            k1 = left_k1_value
            warnings.warn("Pitch-to-diameter ratio is out of interpolation range, so pressure drop computed at a boundary")
        # print("k1_ratio =", k1)
        return Eu_k1_ratio, k1
    

    def _get_Eu_k1_ratio(self, Re=None, pitch_dia_ratio=None):
        r""" Eu_k1 ratio is what textbooks refer to as 'friction factor'
        
        Add docstrings later on how this is computed.
        
        References
        ----------
        [1] A. Zukauskas, "Heat Transfer From Tubes in Cross Flow." 
        In "Handbook of Single Phase Convective Heat Transfer."
        S. Kakac, R.K. Shah, and W. Aung (Eds.), 
        New York: Wiley Interscienc, 1987.
        
        """
        
        if self.config=="inline":
            if pitch_dia_ratio==1.25 and 3<=Re<2e3:
                Eu_k1_ratio = (0.272 + 0.207e3/Re
                               + 0.102e3/Re**2
                               - 0.286e3/Re**3)
            elif pitch_dia_ratio==1.25 and 2e3<=Re<2e6:
                Eu_k1_ratio = (0.267 + 0.249e4/Re
                               - 0.927e7/Re**2
                               - 0.1e11/Re**3)
            elif pitch_dia_ratio==1.5 and 3<=Re<2e3:
                Eu_k1_ratio = (0.263 + 0.867e2/Re
                               - 0.202e1/Re**2)
            elif pitch_dia_ratio==1.5 and 2e3<=Re<2e6:
                Eu_k1_ratio = (0.235 + 0.197e4/Re
                               - 0.124e8/Re**2
                               + 0.312e11/Re**3
                               - 0.274e14/Re**4)
            elif pitch_dia_ratio==2.0 and 7<=Re<800:
                Eu_k1_ratio = (0.188 + 0.566e2/Re
                               - 0.646e3/Re**2
                               + 0.601e4/Re**3
                               - 0.183e5/Re**4)
            elif pitch_dia_ratio==2.0 and 800<=Re<2e6:
                Eu_k1_ratio = (0.247 - 0.595e-6*Re
                               + 0.15e-11*Re**2
                               - 0.137e-17*Re**3
                               + 0.396e-24/Re**4)
            elif pitch_dia_ratio==2.5 and 600<=Re<2e5:
                Eu_k1_ratio = (0.177 - 0.311e-6*Re
                               + 0.117e-11*Re**2)
            else:
                Eu_k1_ratio = None
                print("Eu_k1_ratio in 'inline' configuration could not be computed")
            return Eu_k1_ratio
    
        if self.config=="staggered":
            if pitch_dia_ratio==1.25 and 3<=Re<1e3:
                Eu_k1_ratio = (0.795 + 0.247e3/Re
                               + 0.335e3/Re**2
                               - 0.155e4/Re**3
                               + 0.241e4/Re**4)
            elif pitch_dia_ratio==1.25 and 1e3<=Re<2e6:
                Eu_k1_ratio = (0.245 + 0.339e4/Re
                               - 0.984e7/Re**2
                               + 0.132e11/Re**3
                               - 0.599e13/Re**4)            
            elif pitch_dia_ratio==1.5 and 3<=Re<1e3:
                Eu_k1_ratio = (0.683 + 0.111e3/Re
                               - 0.973e2/Re**2
                               + 0.426e3/Re**3
                               - 0.574e3/Re**4)
            elif pitch_dia_ratio==1.5 and 1e3<=Re<2e6:
                Eu_k1_ratio = (0.203 + 0.248e4/Re
                               - 0.758e7/Re**2
                               + 0.104e11/Re**3
                               - 0.482e13/Re**4)                       
            elif pitch_dia_ratio==2.0 and 7<=Re<1e2:
                Eu_k1_ratio = (0.713 + 0.448e2/Re
                               - 0.126e3/Re**2
                               - 0.582e3/Re**3)
            elif pitch_dia_ratio==2.0 and 1e2<=Re<1e4:
                Eu_k1_ratio = (0.343 + 0.303e3/Re
                               - 0.717e5/Re**2
                               + 0.88e7/Re**3
                               - 0.38e9/Re**4)            
            elif pitch_dia_ratio==2.0 and 1e4<=Re<2e6:
                Eu_k1_ratio = (0.162 + 0.181e4/Re
                               + 0.792e8/Re**2
                               - 0.165e13/Re**3
                               + 0.872e16/Re**4)
            elif pitch_dia_ratio==2.5 and 1e2<=Re<5e3:
                Eu_k1_ratio = (0.33 + 0.989e2/Re
                               - 0.148e5/Re**2
                               + 0.192e7/Re**3
                               - 0.862e8/Re**4)             
            elif pitch_dia_ratio==2.5 and 5e3<=Re<2e6:
                Eu_k1_ratio = (0.119 + 0.498e4/Re
                               - 0.507e8/Re**2
                               + 0.251e12/Re**3
                               - 0.463e15/Re**4)
            else:
                Eu_k1_ratio = None
                print("Eu_k1_ratio in 'staggered' configuration could not be computed")
            return Eu_k1_ratio

    
    def _get_k1_parameter(self, Re=None): # k1 is the geometry parameter
        r""" k1 is what textbooks refer to as 'correction factor'
        
        Add docstrings later on how this is computed.
        
        References
        ----------
        [1] A. Zukauskas, "Heat Transfer From Tubes in Cross Flow." 
        In "Handbook of Single Phase Convective Heat Transfer."
        S. Kakac, R.K. Shah, and W. Aung (Eds.), 
        New York: Wiley Interscienc, 1987.
        
        """
        a = self.transverse_pitch/self.outer_tubediameter
        b = self.longitudanal_pitch/self.outer_tubediameter
        tol = 1e-4 # to compare whether a ~= b
        
        if self.config=="inline":
            if np.isclose(a,b,tol): # if a~b, k1=1
                k1 = 1
                return k1
            
            if not np.isclose(a,b,tol): # if a != b
                param = (a-1)/(b-1)
                if 0.06< param <6:
                    if Re==1e3: k1 = 1.009*(param)**(-0.744)
                    if Re==1e4: k1 = 1.007*(param)**(-0.655)
                    if Re==1e5: k1 = 1.004*(param)**(-0.539)
                    if Re==1e6: k1 = 1.218 - 0.297*(param) + 0.0265*(param)**20
                    if Re==1e7: k1 = 1
                else:
                    k1 = None
                    print("pitch ratios are out of range to compute k1: geometry parameter in 'inline' case")
                return k1

        if self.config=="staggered":
            if np.isclose(a/b, 1.155, tol): # if a/b ~= 1.155: for equilateral triangle, k1=1
                k1 = 1
                return k1
            
            if not np.isclose(a/b, 1.155, tol): # if a/b != 1.155
                param = a/b
                print("param = ", param)
                if 0.5<param<1.2 and Re==1e3:
                    print(" in 1")
                    k1 = (param)**(-0.048)
                elif 0.45<param<3.5 and Re==1e4:
                    print(" in 2")
                    k1 = (1.28 - 0.708/param + 0.55/param**2 - 0.113/param**3)
                    print("k1 in 2 = ", k1)
                elif 0.45<param<3.5 and Re==1e5:
                    print(" in 3")
                    k1 = (2.016 - 1.675*param + 0.948*param**2
                          - 0.234*param**3 + 0.021*param**4)
                    print("k1 in 3 = ", k1)
                elif 0.45<param<3.5 and Re==1e6:
                    print(" in 4")
                    k1 = (2.016 - 1.675*param + 0.948*param**2
                          - 0.234*param**3 + 0.021*param**4)
                elif 1.25<param<3.5 and Re==1e2:
                    print(" in51")
                    k1 = 0.93*(param)**(0.48)
                elif 1.25<param<3.5 and Re==1e3:
                    print(" in 6")
                    k1 = 0.951*(param)**(0.284)
                else:
                    k1 = None
                    print(" in 7")
                    print("pitch ratios are out of range to compute geometry parameter")
  
                return k1 
                
    
    def _get_interval_limits(interval_limits_list, num): 
        r""" This is used by _get_k1_parameter and _get_Eu_k1_ratio
        
        Add docstrings later on how this is computed.
        
        References
        ----------
        [1] A. Zukauskas, "Heat Transfer From Tubes in Cross Flow." 
        In "Handbook of Single Phase Convective Heat Transfer."
        S. Kakac, R.K. Shah, and W. Aung (Eds.), 
        New York: Wiley Interscienc, 1987.
        
        """
        right = bisect.bisect_left(interval_limits_list, num)
        if right==0: # this means num is less than lowest interval limit
            outofrange = True
            left = right # Use lowest available left-most interval limit as the answer
        elif right==len(interval_limits_list): # this means num is greater than largest interval limit
            outofrange = True
            left = right = len(interval_limits_list)-1 # Use right-most largest available interval limit as the answer
        else:
            left = right-1
            outofrange = False
        return left, right, outofrange             
