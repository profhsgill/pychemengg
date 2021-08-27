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
Contains correlations to find Nusselt number
for INTERNAL flow

"""
import numpy as np
from scipy import interpolate
from scipy.optimize import fsolve
import warnings
warnings.simplefilter('always')

__all__ = ["CircularTube"]

class CircularTube:
    r""" For heat transfer from fluid flow through a circular tube.


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
    First import the module **internalflow**.
       
    >>> from pychemengg.heattransfer import internalflow as intflow 
    >>> tube = intflow.CircularTube
    # This will assign the class 'CircularTube' to the variable 'tube'.
    # Methods of the class 'CircularTube' can then be called like so :-
    # tube.method(kwarg1=x, ...)

    """
    
    def Nu_uniformheatflux_laminar(self):
        r""" Nusselt number for laminar flow in tube with constant heat flux.
        
        
        Parameters
        ----------
        `None_required` : 'None'
        

        Returns
        -------
        Nu : `float`
            Nusselt number = 4.36
       
        
        Notes
        -----
        The following formula is used:         
        
        .. math::
            
            Nu = 4.36
    

        Examples
        --------
        First import the module **internalflow**.
           
        >>> from pychemengg.heattransfer import internalflow as intflow 
        >>> tube = intflow.CircularTube()
        >>> tube.Nu_uniformheatflux_laminar()
        4.36
        
        References
        ----------
        [1] Yunus A. Cengel and Afshin J. Ghajar, "Heat And Mass Transfer
        Fundamentals and Applications", 6th Edition. New York, McGraw Hill
        Education, 2020.
        
        """
        return 4.36
    
    def Nu_constantsurfacetemp_laminar(self):
        r""" Nusselt number for laminar flow in tube (constant tube surface
        temperature).
        
        
        Parameters
        ----------
        `None_required` : 'None'
        

        Returns
        -------
        Nu : `float`
            Nusselt number = 3.36
       
        
        Notes
        -----
        The following formula is used:         
        
        .. math::
            
            Nu = 3.36
        

        Examples
        --------
        First import the module **internalflow**.
           
        >>> from pychemengg.heattransfer import internalflow as intflow 
        >>> tube = intflow.CircularTube()
        >>> tube.Nu_constantsurfacetemp_laminar()
        3.36
        
        References
        ----------
        [1] Yunus A. Cengel and Afshin J. Ghajar, "Heat And Mass Transfer
        Fundamentals and Applications", 6th Edition. New York, McGraw Hill
        Education, 2020.
        """
        return 3.66
    
    def Nu_thermallydeveloping_laminar_edwards(self, Re=None, Pr=None,
                                               length=None, diameter=None):
        r"""Nusselt number for thermally developing laminar flow and
        constant tube surface temperature.
        
        Parameters
        ----------
        Re : `int or float`
            Reynolds number for tube/pipe
        Pr : `int or float`
            Prandtl number for the fluid.
        length : `int or float`
            Length of tube
        diameter : `int or float`
            Diameter of tube
        

        Returns
        -------
        Nu : `int or float`
            Average Nusselt number
       
        
        Notes
        -----
        The following formula is used:         
        
        .. math::
            
            Nu = 3.66 + \cfrac {0.065 (D/L) Re Pr}
                               {1 + 0.04[(D/L) Re Pr]^{2/3}}
            
               
        *where:*
        
            Re = Reynolds number
            
            Pr = Prandtl number
            
            D = Pipe inner diameter
            
            L = Pipe length
            
            :math:`Re < 2300` (for laminar flow)
            
            All fluid properties are at bulk mean temp (:math:`T_{b}`):
            
                :math:`T_{b} = (T_{fluid,in} + T_{fluid,out})/2`
                
                except :math:`\mu_{surface}`
                

        
        Warnings
        --------
        A Nusselt number is returned based on the equation
        even if parameters (such as *Re*, *Pr*) do not fall in their
        respective allowable range limits (see above under 'Notes').
        However, if this happens, a warning is issued.
                

        
        Examples
        --------
        First import the module **internalflow**.
           
        >>> from pychemengg.heattransfer import internalflow as intflow 
        >>> tube = intflow.CircularTube()
        >>> tube.Nu_thermallydeveloping_laminar_edwards(Re=77.16, Pr=28750,
                                              length=1500, diameter=0.4)
        13.729061335098443

        
        References
        ----------
        [1] Yunus A. Cengel and Afshin J. Ghajar, "Heat And Mass Transfer
        Fundamentals and Applications", 6th Edition. New York, McGraw Hill
        Education, 2020.
        
        """
        term1 = diameter/length*Re*Pr
        Nu = 3.66 + (0.065*term1)/(1+0.04*np.power(term1, 2/3))
        if not (Re < 2300): warnings.warn("Flow is not laminar. Re is not within range. ")
        return Nu
    
    def Nu_thermallyandhydrodynamicallydeveloping_laminar_siedertate(self, Re=None, Pr=None, length=None,
                                        diameter=None, viscosity_bulk=None, viscosity_surface=None):
        r"""Nusselt number for thermally developing and hydrodynamically 
        developing laminar flow and constant tube surface temperature.
        
        
        Parameters
        ----------
        Re : `int or float`
            Reynolds number for tube/pipe
        Pr : `int or float`
            Prandtl number for the fluid.
        length : `int or float`
            Length of tube
        diameter : `int or float`
            Diameter of tube
        viscosity_bulk : `int or float`
            fluid viscosity at bulk temperature
        viscosity_surface : `int or float`
            fluid viscosity at surface temperature

        Returns
        -------
        Nu : `int or float`
            Average Nusselt number
       
        
        Notes
        -----
        The following formula is used:         
        
        .. math::
            
            Nu = 1.86 \left(\frac{Re Pr D} {L}\right)^{1/3} 
                      \left(\frac{\mu_{bulk}} {\mu_{surface}}\right)^{0.14}
            
               
        *where:*
            Re = Reynolds number
            
            Pr = Prandtl number
            
            D = Pipe inner diameter
            
            L = Pipe length
            
            :math:`\mu_{bulk}` = viscosity at bulk temperature
            
            :math:`\mu_{surface}` = viscosity at surface temperature
            
            :math:`Re < 2300` (for laminar flow)
            
            :math:`0.6 \eqslantless Pr \eqslantless 5`
            
            :math:`0.0044 \eqslantless \left(\frac{\mu_{bulk}} {\mu_{surface}}\right) \eqslantless 9.75`
            
            All fluid properties are at bulk mean temp (:math:`T_{b}`):
            
                :math:`T_{b} = (T_{fluid,in} + T_{fluid,out})/2`
                
                except :math:`\mu_{surface}` is at surface temperature
                

        
        Warnings
        --------
        A Nusselt number is returned based on the equation
        even if parameters (such as *Re*, *Pr*) do not fall in their
        respective allowable range limits (see above under 'Notes').
        However, if this happens, a warning is issued.
        
        
        Examples
        --------
        First import the module **internalflow**.
           
        >>> from pychemengg.heattransfer import internalflow as intflow 
        >>> tube = intflow.CircularTube()
        >>> tube.Nu_thermallyandhydrodynamicallydeveloping_laminar_siedertate(Re=1390, Pr=0.7228,
                                        length=0.1, diameter=0.005, viscosity_bulk=1.963,
                                        viscosity_surface=2.420)
        6.664822962495227

        
        References
        ----------
        [1] Yunus A. Cengel and Afshin J. Ghajar, "Heat And Mass Transfer
        Fundamentals and Applications", 6th Edition. New York, McGraw Hill
        Education, 2020.        
        """
        Nu = 1.86 * np.power(Re*Pr*diameter/length, 1/3) * np.power(viscosity_bulk/viscosity_surface, 0.14)
        if not (Re < 2300): warnings.warn("Flow is not laminar. Re is not within range. ")
        if not (0.6<=Pr<=5): warnings.warn("Pr is not within range.")
        if not (0.0044<=(viscosity_bulk/viscosity_surface)<=9.75):
                warnings.warn("Ratio of viscosity_bulk to viscosity_surface"
                              " is not within range.")
        return Nu
    
    def Nu_dittusboelter(self, Re=None, Pr=None, thermal_ID="heated"):
        r"""Nusselt number for turbulent flow constant tube surface temperature.
        
        
        Parameters
        ----------
        Re : `int or float`
            Reynolds number for tube/pipe
        Pr : `int or float`
            Prandtl number for the fluid.
        thermal_ID : `str`
            Identifier for whether fluid is being "heated" or "cooled".
            Default = "heated"
        

        Returns
        -------
        Nu : `int or float`
            Average Nusselt number
       
        
        Notes
        -----
        The following formula is used:         
        
        .. math::
            
            Nu = 0.023 Re^{0.8} Pr^{n}
        
        *where:*
        
            n = 0.3 when fluid is "cooled"
            
              = 0.4 when fluid is "heated"
            
            Re = Reynolds number
            
            Pr = Prandtl number
            
            :math:`Re > 10000` (for turbulent flow)
            
            All fluid properties are at bulk mean temp (:math:`T_{b}`):
            
                :math:`T_{b} = (T_{fluid,in} + T_{fluid,out})/2`
                

        
        Warnings
        --------
        A Nusselt number is returned based on the equation
        even if parameters (such as *Re*, *Pr*) do not fall in their
        respective allowable range limits (see above under 'Notes').
        However, if this happens, a warning is issued.
        
        
        Examples
        --------
        First import the module **internalflow**.
           
        >>> from pychemengg.heattransfer import internalflow as intflow 
        >>> tube = intflow.CircularTube()
        >>> tube.Nu_dittusboelter(Re=10760, Pr=4.32, thermal_ID="heated")
        69.40186669836865

        
        References
        ----------
        [1] Yunus A. Cengel and Afshin J. Ghajar, "Heat And Mass Transfer
        Fundamentals and Applications", 6th Edition. New York, McGraw Hill
        Education, 2020.        
        """
        if thermal_ID == "heated":
            Nu = 0.023 * np.power(Re, 0.8) * np.power(Pr, 0.4)
        if thermal_ID == "cooled":
            Nu = 0.023 * np.power(Re, 0.8) * np.power(Pr, 0.3)
        if not (Re>=10000): warnings.warn("Flow is not turbulent. Re is not within range.")
        return Nu
        
    def Nu_siedertate(self, Re=None, Pr=None, viscosity_bulk=None, viscosity_surfaceace=None):
        r"""Nusselt number for turbulent flow constant tube surface temperature.
        Use when temperature differences are large.
        
        
        Parameters
        ----------
        Re : `int or float`
            Reynolds number for tube/pipe
        Pr : `int or float`
            Prandtl number for the fluid.
        viscosity_bulk : `int or float`
            fluid viscosity at bulk temperature
        viscosity_surface : `int or float`
            fluid viscosity at surface temperature
        

        Returns
        -------
        Nu : `int or float`
            Average Nusselt number
       
        
        Notes
        -----
        The following formula is used:         
        
        .. math:: 
            Nu = 0.027 Re^{0.8} Pr^{1/3} \left( \frac {\mu_{bulk}} {\mu_{surface}} \right)^{0.14}
        
        *where:*
            Re = Reynolds number
            
            Pr = Prandtl number
            
            :math:`\mu_{bulk}` = viscosity at bulk temperature
            
            :math:`\mu_{surface}` = viscosity at surface temperature
            
            :math:`Re > 10000` (for turbulent flow)
            
            :math:`0.7 \eqslantless Pr \eqslantless 16700`
                        
            All fluid properties are at bulk mean temp (:math:`T_{b}`):
            
                :math:`T_{b} = (T_{fluid,in} + T_{fluid,out})/2`
                
                except :math:`\mu_{surface}` is at surface temperature
                
        
        Warnings
        --------
        A Nusselt number is returned based on the equation
        even if parameters (such as *Re*, *Pr*) do not fall in their
        respective allowable range limits (see above under 'Notes').
        However, if this happens, a warning is issued.
        
        
        Examples
        --------
        First import the module **internalflow**.
           
        >>> from pychemengg.heattransfer import internalflow as intflow 
        >>> tube = intflow.CircularTube()
        >>> tube.Nu_siedertate(self, Re=2.04e5, Pr=3.02,
                               viscosity_bulk=4.71,
                               viscosity_surfaceace=2.82)
        741.7512113817678

        
        References
        ----------
        [1] Yunus A. Cengel and Afshin J. Ghajar, "Heat And Mass Transfer
        Fundamentals and Applications", 6th Edition. New York, McGraw Hill
        Education, 2020.
        [2] M. N. Ozisik, "Heat Transfer A Basic Approach", McGraw Hill, 1985.
        """
        Nu = 0.027 * np.power(Re, 0.8) * np.power(Pr, 1/3) * np.power(viscosity_bulk/viscosity_surfaceace, 0.14)
        if not (Re>=10000): warnings.warn("Flow is not turbulent. Re is not within range.")
        if not (0.7<=Pr<=16700): warnings.warn("Pr is not within range.")
        return Nu
     
    def Nu_gnielinski(self, Re=None, Pr=None):
        r"""Nusselt number for turbulent flow constant tube surface temperature.
        A more accurate correlation.
        
        
        Parameters
        ----------
        Re : `int or float`
            Reynolds number for tube/pipe
        Pr : `int or float`
            Prandtl number for the fluid.
        

        Returns
        -------
        Nu : `int or float`
            Average Nusselt number
       
        
        Notes
        -----
        The following formula is used:         
        
        .. math::
            
            Nu = \cfrac{(f/8) (Re-1000) Pr}{ 1+12.7(f/8)^{0.5} (Pr^{2/3}-1)}
            
        
        *where:*
            Re = Reynolds number
            
            Pr = Prandtl number
            
            f = friction factor (can be determined from different correlations.
                Here following 'Petukhov' correlation is used.)
            
            :math:`f=(0.790 \ln Re - 1.64)^{-2}` (Petukhov correlation)
            
            :math:`3*10^3\eqslantless Re \eqslantless 5*10^6` (for turbulent flow)
            
            :math:`0.5 \eqslantless Pr \eqslantless 2000`
                        
            All fluid properties are at bulk mean temp (:math:`T_{b}`):
            
                :math:`T_{b} = (T_{fluid,in} + T_{fluid,out})/2`
                
        
        Warnings
        --------
        A Nusselt number is returned based on the equation
        even if parameters (such as *Re*, *Pr*) do not fall in their
        respective allowable range limits (see above under 'Notes').
        However, if this happens, a warning is issued.
        
        
        Examples
        --------
        First import the module **internalflow**.
           
        >>> from pychemengg.heattransfer import internalflow as intflow 
        >>> tube = intflow.CircularTube()
        >>> tube.Nu_siedertate(Re=2.04e5, Pr=3.02,
                               viscosity_bulk=4.71,
                               viscosity_surfaceace=2.82)
        741.7512113817678

        
        References
        ----------
        [1] Yunus A. Cengel and Afshin J. Ghajar, "Heat And Mass Transfer
        Fundamentals and Applications", 6th Edition. New York, McGraw Hill
        Education, 2020.
        [2] M. N. Ozisik, "Heat Transfer A Basic Approach", McGraw Hill, 1985.
        """
        f = np.power(0.79*np.log(Re) - 1.64, -2)
        num = (f/8)*(Re-1000)*Pr
        denom = 1 + 12.7*np.power(f/8, 0.5)*(np.power(Pr,2/3)-1)
        Nu = num/denom
        if not (3e3<=Re<=5e6): warnings.warn("Re is not within range.")
        if not (0.5<=Pr<=2000): warnings.warn("Pr is not within range.")
        return Nu

    def Nu_annulus_laminar(self, outerannulus_diameter=None,
                           innerannulus_diameter=None,
                           findNUat_ID="innerannulus_diameter"):
        
        r"""Nusselt number for laminar flow in annulus (one surface is
        isothermal and other adiabatic)
        
        
        Parameters
        ----------
        outerannulus_diameter : `int or float`
            Outer annulus diameter for fluid flow.
        innerannulus_diameter : `int or float`
            Inner annulus diameter for fluid flow.
        findNUat_ID : `str`
            Annulus diameter where  Nu is to be found. 
            Default = "innerannulus_diameter", 
            second option = "outerannulus_diameter"
        

        Returns
        -------
        Nu : `int or float`
            Average Nusselt number
       
        
        Notes
        -----
        A look table is provided in textbooks with ratio of
        outerannulus_diameter/innerannulus_diameter in one column and
        Nu_at_innerannulus_diameter and Nu_at_outerannulus_diameter in other
        two columns. This table is implemented and values are interpolated
        as required.
            
        
        Examples
        --------
        First import the module **internalflow**.
           
        >>> from pychemengg.heattransfer import internalflow as intflow 
        >>> tube = intflow.CircularTube()
        >>> tube.Nu_annulus_laminar(outerannulus_diameter=0.1,
                                    innerannulus_diameter=0.025,
                                    findNUat_ID="innerannulus_diameter")
        741.7512113817678

        
        References
        ----------
        [1] Yunus A. Cengel and Afshin J. Ghajar, "Heat And Mass Transfer
        Fundamentals and Applications", 6th Edition. New York, McGraw Hill
        Education, 2020.
        """
        di_do_vals = np.array([0, 0.05, 0.1, 0.25, 0.5, 1])
        Nu_i_vals = np.array([np.nan, 17.46, 11.56, 7.37, 5.74, 4.86])
        Nu_o_vals = np.array([3.66, 4.06, 4.11, 4.23, 4.43, 4.86])
        fxn_Nu_i = interpolate.interp1d(di_do_vals, Nu_i_vals)
        fxn_Nu_o = interpolate.interp1d(di_do_vals, Nu_o_vals)
        if findNUat_ID=="innerannulus_diameter":
            Nu = fxn_Nu_i(innerannulus_diameter/outerannulus_diameter) 
        if findNUat_ID=="outerannulus_diameter":
            Nu = fxn_Nu_o(innerannulus_diameter/outerannulus_diameter)   
        if np.isnan(Nu):
            print("For inner to outer annulus diameter ratio of "
                  "< 0.05, data for 'Nu' is not available")
        else:
            return float(Nu) # Nu is an np.array of float.
                             # Convert to Python float type before returning.
    
    def Nu_liquidmetal_isothermal(self, Re=None, Pr_surface=None):
        r"""Nusselt number for liquid metals (isothermal case).
        
        Parameters
        ----------
        Re : `int or float`
            Reynolds number.
        Pr_surface : `int or float`
            Prandtl number at surface temperature.
        

        Returns
        -------
        Nu : `int or float`
            Average Nusselt number
       
        
        Notes
        -----
        The following formula is used:         
        
        .. math::
            
            Nu = 4.8 + 0.0156 Re^{0.85} Pr_{surface}^{0.93}
        
        *where:*
            Re = Reynolds number
            
            :math:`Pr_{surface}` = Prandtl number at surface temperature
            
            :math:`10^4\eqslantless Re \eqslantless 10^6` (for turbulent flow)
            
            :math:`0.004 \eqslantless Pr \eqslantless 0.01`
                        
            All fluid properties are at bulk mean temp (:math:`T_{b}`):
            
                :math:`T_{b} = (T_{fluid,in} + T_{fluid,out})/2`
                
                except :math:`\Pr_{surface}` is at surface temperature
            
        
        Examples
        --------
        First import the module **internalflow**.
           
        >>> from pychemengg.heattransfer import internalflow as intflow 
        >>> tube = intflow.CircularTube()
        >>> tube. Nu_liquidmetal_isothermal(Re=13570, Pr_surface=0.0119)
        5.624282620227201

        
        References
        ----------
        [1] Yunus A. Cengel and Afshin J. Ghajar, "Heat And Mass Transfer
        Fundamentals and Applications", 6th Edition. New York, McGraw Hill
        Education, 2020.
        """
        Nu = 4.8 + 0.0156*np.power(Re,0.85)*np.power(Pr_surface, 0.93)
        if not (1e4<=Re<=1e6): warnings.warn("Re is not within range.")
        return Nu
    
    def Nu_liquidmetal_constflux(self, Re=None, Pr_surface=None):
        r"""Nusselt number for liquid metals (constant flux case).
        
        Parameters
        ----------
        Re : `int or float`
            Reynolds number.
        Pr_surface : `int or float`
            Prandtl number at surface temperature.
        

        Returns
        -------
        Nu : `int or float`
            Average Nusselt number.
       
        
        Notes
        -----
        The following formula is used:         
        
        .. math::
            
            Nu = 6.3 + 0.0167 Re^{0.85} Pr_{surface}^{0.93}
        
        *where:*
            Re = Reynolds number
            
            :math:`Pr_{surface}` = Prandtl number at surface temperature
            
            :math:`10^4\eqslantless Re \eqslantless 10^6` (for turbulent flow)
            
            :math:`0.004 \eqslantless Pr \eqslantless 0.01`
                        
            All fluid properties are at bulk mean temp (:math:`T_{b}`):
            
                :math:`T_{b} = (T_{fluid,in} + T_{fluid,out})/2`
                
                except :math:`\Pr_{surface}` is at surface temperature
            
        
        Examples
        --------
        First import the module **internalflow**.
           
        >>> from pychemengg.heattransfer import internalflow as intflow 
        >>> tube = intflow.CircularTube()
        >>> tube.Nu_liquidmetal_constflux(Re=13570, Pr_surface=0.0119)
        7.182405112679119

        
        References
        ----------
        [1] Yunus A. Cengel and Afshin J. Ghajar, "Heat And Mass Transfer
        Fundamentals and Applications", 6th Edition. New York, McGraw Hill
        Education, 2020.
        """
        Nu = 6.3 + 0.0167*np.power(Re,0.85)*np.power(Pr_surface, 0.93)
        if not (1e4<=Re<=1e6): warnings.warn("Re is not within range.")
        return Nu
    
    def pressuredrop(self, length=None, diameter=None, frictionfactor=None,
                     density=None, velocity=None):
        r"""Computes pressure drop for flow in a pipe.
        
        Parameters
        ----------
        length : `int or float`
            Length of pipe.
        diameter : `int or float`
            Diameter of pipe.
        frictionfactor: `int or float`
            Friction factor of flowing fluid in pipe.
        density : `int or float`
            Density of fluid.
        velocity : `int or float`
            Velocity of fluid in pipe.
            
        
        Returns
        -------
        pressure drop : `int or float`
            Pressure drop from flowing fluid in pipe
       
        
        Notes
        -----
        The following formula is used:         
        
        .. math::
            
            \Delta P = f \frac{L}{D} \frac{\rho V_{avg}^2}{2}
        
        *where:*
            
            :math:`\rho` = density
            
            :math:`V_{avg}` = average velocity
            
            L = pipe length
            
            D = pipe inner diameter
            
            f = friction factor
            
            :math:`\Delta P` = pressure drop in pipe

                        
            All fluid properties are at bulk mean temp (:math:`T_{b}`):
            
                :math:`T_{b} = (T_{fluid,in} + T_{fluid,out})/2`

            
        
        Examples
        --------
        First import the module **internalflow**.
           
        >>> from pychemengg.heattransfer import internalflow as intflow 
        >>> tube = intflow.CircularTube()
        >>> tube.pressuredrop(length=200, diameter=0.3, frictionfactor=0.1006,
                         density=888.1, velocity=2)
        119123.81333333332

        
        References
        ----------
        [1] Yunus A. Cengel and Afshin J. Ghajar, "Heat And Mass Transfer
        Fundamentals and Applications", 6th Edition. New York, McGraw Hill
        Education, 2020.
        """
        return frictionfactor *length/diameter/2*density*velocity**2
    
    def frictionfactor_laminar(self, Re=None):
        r"""Computes friction factor for laminar flow in a pipe.
        
        Parameters
        ----------
        Re : `int or float`
            Reynolds number
        
            
        Returns
        -------
        friction factor : `int or float`
            Friction factor for flowing fluid in pipe
       
        
        Notes
        -----
        The following formula is used:         
        
        .. math::
            
            f = \frac{64}{Re}
        
        *where:*
            
            Re = Reynolds number
            
            f = friction factor
            
            :math:`Re < 2300` (for laminar flow)
                        
            All fluid properties are at bulk mean temp (:math:`T_{b}`):
            
                :math:`T_{b} = (T_{fluid,in} + T_{fluid,out})/2`

        
        Warnings
        --------
        A friction factor is returned based on the equation
        even if parameters (such as *Re*) do not fall in their
        respective allowable range limits (see above under 'Notes').
        However, if this happens, a warning is issued.
        
        
        Examples
        --------
        First import the module **internalflow**.
           
        >>> from pychemengg.heattransfer import internalflow as intflow 
        >>> tube = intflow.CircularTube()
        >>> tube.frictionfactor_laminar(636)
        0.10062893081761007

        
        References
        ----------
        [1] Yunus A. Cengel and Afshin J. Ghajar, "Heat And Mass Transfer
        Fundamentals and Applications", 6th Edition. New York, McGraw Hill
        Education, 2020.
        """
        if not (Re<2300): warnings.warn("Flow is not laminar."
                                        " Re is not within range.")
        return 64/Re
    
    
    def frictionfactor_turbulent(self, Re=None, diameter=None, 
                                 roughness=None, surfacetype=None):
        r"""Computes friction factor for turbulent flow in a pipe.
        
        Parameters
        ----------
        Re : `int or float`
            Reynolds number.
        diameter : `int or float`
            Inner diameter of pipe.
        roughness : `int or float`
            Roughness of pipe.
        surfacetype : `str`
            Indicates if surface is to be considered
            "rough" or "smooth". At least one of the two must be entered.
            
        Returns
        -------
        friction factor : `int or float`
            Friction factor for flowing fluid in pipe
       
        
        Notes
        -----
        The following formula is used:         
        
        for "smooth" surface
            
        .. math::
            
            f=(0.790 \ln Re - 1.64)^{-2} \hspace{5pt}(Petukhov\hspace{5pt}correlation)
            
        for "rough surface"
        
        .. math::
            
            \frac{1}{\sqrt{f}} = ( -2.0 \log_{10} \left( \frac{\epsilon /D}{3.7}
                        + \frac{2.51}{Re \sqrt{f}}\right) 
                                  
                                  \hspace{5pt}(Colebrook\hspace{5pt}correlation)
        
        *where:*
            
            Re = Reynolds number
            
            f = friction factor
            
            :math:`\epsilon` = surface roughness of inside of pipes
                        
            All fluid properties are at bulk mean temp (:math:`T_{b}`):
            
                :math:`T_{b} = (T_{fluid,in} + T_{fluid,out})/2`

        
        Warnings
        --------
        A friction factor is returned based on the equation
        even if parameters (such as *Re*) do not fall in their
        respective allowable range limits (see above under 'Notes').
        However, if this happens, a warning is issued.
        
        
        Examples
        --------
        First import the module **internalflow**.
           
        >>> from pychemengg.heattransfer import internalflow as intflow 
        >>> tube = intflow.CircularTube()
        >>> tube.frictionfactor_turbulent(Re=126400, diameter=2/12, 
                                     roughness=7e-6, surfacetype="rough")
        0.017397627070796194

        
        References
        ----------
        [1] Yunus A. Cengel and Afshin J. Ghajar, "Heat And Mass Transfer
        Fundamentals and Applications", 6th Edition. New York, McGraw Hill
        Education, 2020.
        """
        correctchoice = ["rough", "smooth"]
        if surfacetype not in correctchoice:
            print("You must enter 'surfacetype as 'smooth' OR 'rough'")
            frictionfactor = None
        elif surfacetype == "smooth": # Petukhov equation
            frictionfactor = np.power((0.79*np.log(Re) - 1.64),-2)
        elif surfacetype == "rough":
            def colebrook(f):
                f=f+0.j
                LHS = 1/np.power(f, 0.5)
                RHS = -2*np.log10(roughness/diameter/3.7 + 2.51/Re/np.power(f,0.5))
                return (LHS-RHS).real
            f = fsolve(colebrook, 0.05)
            frictionfactor = f[0]
        if not (Re>2300): warnings.warn("Flow is not turbulent."
                                        " Re is not within range.")
        return frictionfactor
        
    def hydrodynamic_entrylength_laminar(self, Re=None, diameter=None):
        r"""Computes laminar hydrodynamic entry length.
        
        Parameters
        ----------
        Re : `int or float`
            Reynolds number.
        diameter : `int or float`
            Inner diameter of pipe.
        
            
        Returns
        -------
        entry length : `int or float`
            Hydrodynamic entry length.
       
        
        Notes
        -----
        The following formula is used:         
        
        .. math::
            
            L_{h} = 0.05 Re D

        
        *where:*
            
            Re = Reynolds number
            
            D = internal diameter of pipe
            
            :math:`L_h` = hydrodynamic entry length
                        
            All fluid properties are at bulk mean temp (:math:`T_{b}`):
            
                :math:`T_{b} = (T_{fluid,in} + T_{fluid,out})/2`

        
        Warnings
        --------
        A friction factor is returned based on the equation
        even if parameters (such as *Re*) do not fall in their
        respective allowable range limits (see above under 'Notes').
        However, if this happens, a warning is issued.
        
        
        Examples
        --------
        First import the module **internalflow**.
           
        >>> from pychemengg.heattransfer import internalflow as intflow 
        >>> tube = intflow.CircularTube()
        >>> tube.hydrodynamic_entrylength_laminar(Re=636, diameter=0.3)
        9.54

        
        References
        ----------
        [1] Yunus A. Cengel and Afshin J. Ghajar, "Heat And Mass Transfer
        Fundamentals and Applications", 6th Edition. New York, McGraw Hill
        Education, 2020.
        """
        if not (Re<2300): warnings.warn("Flow is not laminar."
                                        " Re is not within range.")
        return 0.05*Re*diameter
    
    def hydrodynamic_entrylength_turbulent(self, diameter=None):
        r"""Computes turbulent hydrodynamic entry length.
        
        Parameters
        ----------
        diameter : `int or float`
            Inner diameter of pipe.
        
            
        Returns
        -------
        entry length : `int or float`
            Hydrodynamic entry length.
       
        
        Notes
        -----
        The following formula is used:         
        
        .. math::
            
            L_{h} = 10 D

        
        *where:*
            
            D = internal diameter of pipe
            
            :math:`L_h` = hydrodynamic entry length
        
        
        Examples
        --------
        First import the module **internalflow**.
           
        >>> from pychemengg.heattransfer import internalflow as intflow 
        >>> tube = intflow.CircularTube()
        >>> tube.hydrodynamic_entrylength_turbulent(diameter=0.03)
        0.3

        
        References
        ----------
        [1] Yunus A. Cengel and Afshin J. Ghajar, "Heat And Mass Transfer
        Fundamentals and Applications", 6th Edition. New York, McGraw Hill
        Education, 2020.
        """
        return 10*diameter
    
    def thermal_entrylength_laminar(self, Re=None, Pr=None, diameter=None):
        r"""Computes laminar thermal entry length.
        
        Parameters
        ----------
        Re : `int or float`
            Reynolds number.
        Pr : `int or float`
            Prandtl number.    
        diameter : `int or float`
            Inner diameter of pipe.
        
            
        Returns
        -------
        entry length : `int or float`
            Thermal entry length.
       
        
        Notes
        -----
        The following formula is used:         
        
        .. math::
            
            L_{t} = 0.05 Re Pr D

        
        *where:*
            
            Re = Reynolds number
            
            Pr = Prandtl number
            
            D = internal diameter of pipe
            
            :math:`L_t` = thermal entry length
                        
            All fluid properties are at bulk mean temp (:math:`T_{b}`):
            
                :math:`T_{b} = (T_{fluid,in} + T_{fluid,out})/2`

        
        Warnings
        --------
        A friction factor is returned based on the equation
        even if parameters (such as *Re*) do not fall in their
        respective allowable range limits (see above under 'Notes').
        However, if this happens, a warning is issued.
        
        
        Examples
        --------
        First import the module **internalflow**.
           
        >>> from pychemengg.heattransfer import internalflow as intflow 
        >>> tube = intflow.CircularTube()
        >>> tube.thermal_entrylength_laminar(Re=636, Pr=10863, diameter=0.3)
        103633.02

        
        References
        ----------
        [1] Yunus A. Cengel and Afshin J. Ghajar, "Heat And Mass Transfer
        Fundamentals and Applications", 6th Edition. New York, McGraw Hill
        Education, 2020.
        """
        if not (Re<2300): warnings.warn("Flow is not laminar."
                                        " Re is not within range.")
        return 0.05*Re*Pr*diameter
    
    def thermal_entrylength_turbulent(self, diameter=None):
        r"""Computes turbulent thermal entry length.
        
        Parameters
        ----------
        diameter : `int or float`
            Inner diameter of pipe.
        
            
        Returns
        -------
        entry length : `int or float`
            Hydrodynamic entry length.
       
        
        Notes
        -----
        The following formula is used:         
        
        .. math::
            
            L_{t} = 10 D

        
        *where:*
            
            D = internal diameter of pipe
            
            :math:`L_t` = thermal entry length
        
        
        Examples
        --------
        First import the module **internalflow**.
           
        >>> from pychemengg.heattransfer import internalflow as intflow 
        >>> tube = intflow.CircularTube()
        >>> tube.thermal_entrylength_turbulent(diameter=0.03)
        0.3

        
        References
        ----------
        [1] Yunus A. Cengel and Afshin J. Ghajar, "Heat And Mass Transfer
        Fundamentals and Applications", 6th Edition. New York, McGraw Hill
        Education, 2020.
        """
        return 10*diameter
