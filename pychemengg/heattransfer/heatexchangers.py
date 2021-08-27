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
Module for heat exchanger analysis.

"""

import math
from scipy.optimize import fsolve


__all__ = ["FCorrectionFactor", "EffNTU", "calc_overallheattransfercoefficient_fromNTU"]

class FCorrectionFactor:
    r""" Contains functions to compute 'F' correction factor for LMTD
    for heat exchangers.


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
    First import the module **heatexchangers**.
       
    >>> from pychemengg.heattransfer import heatexchangers as hx 
    >>> hx1 = hx.FCorrectionFactor
    # This will assign the class 'FCorrectionFactor' to the
    # variable 'hx1'.
    # Methods of the class 'FCorrectionFactor' can then be called like so :-
    # hx1.methodname(kwarg1=x, ... etc)

    """
    
    def oneshell2ntubepasses(self, T_tubein=None, T_tubeout=None, T_shellin=None, T_shellout=None):
        r""" To find 'F' correction factor for LMTD.
        Use when shell passes = 1 and tube passes = 2, 4, 6, 8, etc

        
        Parameters
        ----------
        T_tubein : `int or float`
            Temperature of tube side (or cold) fluid at inlet.
        T_tubeout : `int or float`
            Temperature of tube side (or cold) fluid at outlet.
        T_shellin : `int or float`
            Temperature of shell side (or hot) fluid at inlet.
        T_shellout : `int or float`
            Temperature of shell side (or hot) fluid at outlet.        
        

        Returns
        -------
        F : `int or float`
            'F' correction factor for LMTD.
       
        
        Raises
        ------
        ValueError with message 'math domain error'
            If negative number is being passed to the 'log' function
        
        
        Notes
        -----
        Equation number (6) of reference [1] is used.

        This is the same equation that is used to generate plots presented
        in most textbooks.        
        
               
        Examples
        --------
        First import the module **heatexchangers**.
           
        >>> from pychemengg.heattransfer import heatexchangers as hx 
        >>> hx1 = hx.FCorrectionFactor()
        >>> hx1.oneshell2ntubepasses(T_shellin=300, T_shellout=200,
                                     T_tubein=100, T_tubeout=200)
        0.8022781617244771

        

        References
        ----------
        [1] R. A. Bowman, A. C. Mueller, and W. M. Nagle, "Mean Temperature 
        Difference in Design", Transactions of the ASME 62, 1940, pp:283-294.        .

        """
        try:
            R, P = FCorrectionFactor.calc_R_P(self, T_tubein=T_tubein,
                                          T_tubeout=T_tubeout,
                                          T_shellin=T_shellin, 
                                          T_shellout=T_shellout)
            ta = math.sqrt(R**2 + 1)
            tb = 2/P - 1 - R
            if R == 1:
                tt = P / (math.log(10)*(1-P))
            elif R != 1:
                arg1 = (1-P) / (1-P*R)
                tt = 1/(R-1) * math.log10(arg1)
            tc = ta * tt
            td = (tb + ta) / (tb - ta)
            f1_2 = tc/math.log10(td)
            return f1_2
        except ValueError as err:
            if str(err) == "math domain error":
                print(err)
                print("Most likely a negative number is being passed to the"
                      " 'log' function.")
                print("A heat exchanger of this type cannot be designed"
                      " for this particular temperature combination")

        
        
    def twoshell4ntubepasses(self, T_tubein=None, T_tubeout=None, T_shellin=None, T_shellout=None):
        r"""To find 'F' correction factor for LMTD.
        Use when shell passes = 2 and tube passes = 4, 8, 12, etc

        
        Parameters
        ----------
        T_tubein : `int or float`
            Temperature of tube side (or cold) fluid at inlet.
        T_tubeout : `int or float`
            Temperature of tube side (or cold) fluid at outlet.
        T_shellin : `int or float`
            Temperature of shell side (or hot) fluid at inlet.
        T_shellout : `int or float`
            Temperature of shell side (or hot) fluid at outlet.        
        

        Returns
        -------
        F : `int or float`
            'F' correction factor for LMTD.
       
        
        Raises
        ------
        ValueError with message 'math domain error'
            If negative number is being passed to the 'log' function
        
        
        Notes
        -----
        Equation number (8) of reference [1] is used.
        
        This is the same equation used to generate plots presented
        in most textbooks.         
        
               
        Examples
        --------
        First import the module **heatexchangers**.
           
        >>> from pychemengg.heattransfer import heatexchangers as hx 
        >>> hx1 = hx.FCorrectionFactor()
        >>> hx1.twoshell4ntubepasses(T_shellin=300, T_shellout=200,
                                     T_tubein=100, T_tubeout=200)
        0.9568453972970873

        

        References
        ----------
        [1] R. A. Bowman, A. C. Mueller, and W. M. Nagle, "Mean Temperature 
        Difference in Design", Transactions of the ASME 62, 1940, pp:283-294.        .

        """
        try:
            R, P = FCorrectionFactor.calc_R_P(self, T_tubein=T_tubein,
                                              T_tubeout=T_tubeout,
                                              T_shellin=T_shellin,
                                              T_shellout=T_shellout)
            ta = math.sqrt(R**2 + 1)
            tb = 2/P - 1 - R
            tc = 2/P * math.sqrt((1-P) * (1-P*R))
            if R == 1:
                tt = P / (math.log(10)*(1-P))
                tt = tt/2
            elif R != 1:
                arg1 = (1-P) / (1-P*R)
                tt = 1/(R-1) * math.log10(arg1)
                tt = tt/2
            tf = ta * tt
            td = (tb + tc + ta) / (tb + tc - ta)
            f2_2 = tf/math.log10(td)
            return f2_2
        except ValueError as err:
            if str(err) == "math domain error":
                print(err)
                print("Most likely a negative number is being passed to the"
                      " 'log' function.")
                print("A heat exchanger of this type cannot be designed"
                      " for this particular temperature combination")

    
    
    def calc_R_P(self, T_tubein=None, T_tubeout=None, T_shellin=None, T_shellout=None):
        r"""To find 'R' and 'P' parameters needed to compute 'F' correction 
        factors. Could be useful to manually look up tables
        of F correction factors in textbooks.

        
        Parameters
        ----------
        T_tubein : `int or float`
            Temperature of tube side (or cold) fluid at inlet.
        T_tubeout : `int or float`
            Temperature of tube side (or cold) fluid at outlet.
        T_shellin : `int or float`
            Temperature of shell side (or hot) fluid at inlet.
        T_shellout : `int or float`
            Temperature of shell side (or hot) fluid at outlet.        
        

        Returns
        -------
        tuple containg R-parameter and P-parameter :
            (R-parameter : `int or float`, P-parameter : `int or float`)
       
        
        Notes
        -----
        The following formulas defined in reference [1] are used:
        
        .. math::
            
            R = \frac {T_{shellin} - T_{shellout}}  {T_{tubeout} - T_{tubein}}
            
            P = \frac {T_{tubeout} - T_{tubein}} {T_{shellin} - T_{tubein}}
        
               
        Examples
        --------
        First import the module **heatexchangers**.
           
        >>> from pychemengg.heattransfer import heatexchangers as hx 
        >>> hx1 = hx.FCorrectionFactor()
        >>> hx1.calc_R_P(T_shellin=300, T_shellout=200,
                                     T_tubein=100, T_tubeout=200)
        (1.0, 0.5)

        

        References
        ----------
        [1] R. A. Bowman, A. C. Mueller, and W. M. Nagle, "Mean Temperature 
        Difference in Design", Transactions of the ASME 62, 1940, pp:283-294.        .

        """
        R = (T_shellin - T_shellout) / (T_tubeout - T_tubein)
        P = (T_tubeout - T_tubein) / (T_shellin - T_tubein)
        return R, P


    def singlepass_crossflow_bothunmixed(self, T_tubein=None, T_tubeout=None,
                                         T_shellin=None, T_shellout=None):
        r"""To find 'F' correction factor for LMTD.
        Use for single pass cross flow with both fluids unmixed.

        
        Parameters
        ----------
        T_tubein : `int or float`
            Temperature of tube side (or cold) fluid at inlet.
        T_tubeout : `int or float`
            Temperature of tube side (or cold) fluid at outlet.
        T_shellin : `int or float`
            Temperature of shell side (or hot) fluid at inlet.
        T_shellout : `int or float`
            Temperature of shell side (or hot) fluid at outlet.        
        

        Returns
        -------
        F : `int or float`
            'F' correction factor for LMTD.
       
        
        Notes
        -----
        Equation number (10) of reference [1] is used.
        
        This is the same equation that is used to generate plots presented
        in most textbooks.         
        
               
        Examples
        --------
        First import the module **heatexchangers**.
           
        >>> from pychemengg.heattransfer import heatexchangers as hx 
        >>> hx1 = hx.FCorrectionFactor()
        >>> hx1.singlepass_crossflow_bothunmixed(T_shellin=300, T_shellout=200,
                                                 T_tubein=100, T_tubeout=200)
        0.8945911509910063

        

        References
        ----------
        [1] R. A. Bowman, A. C. Mueller, and W. M. Nagle, "Mean Temperature 
        Difference in Design", Transactions of the ASME 62, 1940, pp:283-294.        .

        """
        try:
            T1 = T_shellin
            T2 = T_shellout
            t1 = T_tubein
            t2 = T_tubeout
            p = (T1-T2)/(T1-t1)
            q = (t2-t1)/(T1-t1)
            if p!=q:
                r0 = (p-q)/math.log((1-q)/(1-p)) # r0 stand for countercurrent configuration
            if p==q: # i.e. R = 1
                r0 = (T1-t2)/(T1-t1) # r0 stand for countercurrent configuration
                                     # when R=1 , deltaTm = LMTD = T1-t2 
                                     # so this term is in numerator
            def find_r(r): # this is equation 10 from Bowman, 1940
                n = 40
                summation = 0
                fact = math.factorial
                for u in range(n):
                    for v in range(n):
                        term = (-1)**(u+v) * fact(u+v)/fact(u)/fact(v)/fact(u+1)/fact(v+1)
                        term = term * math.pow(p/r, u) * math.pow(q/r, v)
                        summation = summation + term
                return r-summation
            ans = fsolve(find_r, 0.5)
            F = ans[0]/r0
            # print("p= ",p, " q= ", q, " F= ", F)
            return F
        except ValueError as err:
            if str(err) == "math domain error":
                print(err)
                print("Most likely a negative number is being passed to the"
                      " 'log' function.")
                print("A heat exchanger of this type cannot be designed"
                      " for this particular temperature combination")
    
    def singlepass_crossflow_oneunmixed(self, T_tubein=None, T_tubeout=None,
                                        T_shellin=None, T_shellout=None):
        r"""To find 'F' correction factor for LMTD.
        Use for single pass cross flow with one fluid unmixed.

        
        Parameters
        ----------
        T_tubein : `int or float`
            Temperature of tube side (or cold) fluid at inlet.
        T_tubeout : `int or float`
            Temperature of tube side (or cold) fluid at outlet.
        T_shellin : `int or float`
            Temperature of shell side (or hot) fluid at inlet.
        T_shellout : `int or float`
            Temperature of shell side (or hot) fluid at outlet.        
        

        Returns
        -------
        F : `int or float`
            'F' correction factor for LMTD.
       
        
        Notes
        -----
        Equation number (11) of reference [1] is used.
        
        This is the same equation that is used to generate plots presented
        in most textbooks.         
        
               
        Examples
        --------
        First import the module **heatexchangers**.
           
        >>> from pychemengg.heattransfer import heatexchangers as hx 
        >>> hx1 = hx.FCorrectionFactor()
        >>> hx1.singlepass_crossflow_oneunmixed(T_shellin=300, T_shellout=200,
                                                 T_tubein=100, T_tubeout=200)
        0.8464626304853572

        

        References
        ----------
        [1] R. A. Bowman, A. C. Mueller, and W. M. Nagle, "Mean Temperature 
        Difference in Design", Transactions of the ASME 62, 1940, pp:283-294.        .

        """
        # this is equation 11 from Bowman, 1940
        try:
            T1 = T_shellin
            T2 = T_shellout
            t1 = T_tubein
            t2 = T_tubeout
            P = (T1-T2)/(T1-t1)
            Q = (t2-t1)/(T1-t1)
            if P!=Q:
                r0 = (P-Q)/math.log((1-Q)/(1-P))
            if P==Q:
                r0 = (T1-t2)/(T1-t1)
            b = 1-(Q/P*math.log(1/(1-P)))
            R = Q/math.log(1/(b)) 
            F = R/r0
            return F
        except ValueError as err:
            if str(err) == "math domain error":
                print(err)
                print("Most likely a negative number is being passed to the"
                      " 'log' function.")
                print("A heat exchanger of this type cannot be designed"
                      " for this particular temperature combination")
    
    def singlepass_crossflow_bothmixed(self, T_tubein=None, T_tubeout=None,
                                       T_shellin=None, T_shellout=None):
        r"""To find 'F' correction factor for LMTD.
        Use for single pass cross flow with both fluids mixed.

        
        Parameters
        ----------
        T_tubein : `int or float`
            Temperature of tube side (or cold) fluid at inlet.
        T_tubeout : `int or float`
            Temperature of tube side (or cold) fluid at outlet.
        T_shellin : `int or float`
            Temperature of shell side (or hot) fluid at inlet.
        T_shellout : `int or float`
            Temperature of shell side (or hot) fluid at outlet.        
        

        Returns
        -------
        F : `int or float`
            'F' correction factor for LMTD.
       
        
        Notes
        -----
        Equation number (12) of reference [1] is used.
        
        This is the same equation that is used to generate plots presented
        in most textbooks.         
        
               
        Examples
        --------
        First import the module **heatexchangers**.
           
        >>> from pychemengg.heattransfer import heatexchangers as hx 
        >>> hx1 = hx.FCorrectionFactor()
        >>> hx1.singlepass_crossflow_bothmixed(T_shellin=300, T_shellout=200,
                                                 T_tubein=100, T_tubeout=200)
        0.7959050946318332

        

        References
        ----------
        [1] R. A. Bowman, A. C. Mueller, and W. M. Nagle, "Mean Temperature 
        Difference in Design", Transactions of the ASME 62, 1940, pp:283-294.        .

        """
        # this is equation 12 from Bowman, 1940
        try:
            T1 = T_shellin
            T2 = T_shellout
            t1 = T_tubein
            t2 = T_tubeout
            P = (T1-T2)/(T1-t1)
            Q = (t2-t1)/(T1-t1)
            if P!=Q:
                r0 = (P-Q)/math.log((1-Q)/(1-P))
            if P==Q:
                r0 = (T1-t2)/(T1-t1)
            def findR(R): # this is equation 10 from Bowman, 1940
                ta = (P/R)/(1-math.exp(-P/R))
                tb = (Q/R)/(1-math.exp(-Q/R))
                return R*(ta + tb - 1) - 1
            ans = fsolve(findR, 0.5)
            F = ans[0]/r0
            return F
        except ValueError as err:
            if str(err) == "math domain error":
                print(err)
                print("Most likely a negative number is being passed to the"
                      " 'log' function.")
                print("A heat exchanger of this type cannot be designed"
                      " for this particular temperature combination")


############################
# Following functions need more work before adding to API
############################

    def _find_P_equivalent (self, x, R=None, P=None, shellpass_count=None):
        """ 
            VERIFY AGAIN BEFORE ADDING TO PUBLIC API
        """
        n = shellpass_count
        if R != 1:
            term1 = ((1-x*R) / (1-x))**n
            Pnew = (1-term1) / (R-term1)
            return Pnew-P
        if R == 1:
            Pnew = (x * n) / (x*n - x + 1)
            return Pnew-P 

            
    
        
    def _nshell2ntube(self, T_tubein=None, T_tubeout=None, T_shellin=None,
                     T_shellout=None, shellpass_count=None):
        """ Use when shell passes >= 3
            and tube passes = even multiples of shell pass
            example : 3 shell - 6, 12, 18 ... tube passes
                      4 shell - 8, 12, 16 ... tube passes
                      ... and so on
                      
            VERIFY AGAIN BEFORE ADDING TO PUBLIC API
            
        """
        try:
            n = shellpass_count
            R, P = FCorrectionFactor.calc_R_P(self, T_tubein=T_tubein,
                                              T_tubeout=T_tubeout,
                                              T_shellin=T_shellin,
                                              T_shellout=T_shellout)
            # First we convert P for the 'n' shell '2N' tube exchanger
            # into 'P' that is equivalent to '1' shell '2' tube exchanger                    
            psolved = fsolve(lambda x: FCorrectionFactor._find_P_equivalent(x, R=R, P=P, shellpass_count=n), 0.5)
            pequivalent = psolved[0]
            print("equivalent P =", pequivalent)
            
            # Now we take this P-equivalent and use 1 shell - 2 tube formula to find 'F'
            # To do this we find equivalent T_tubein and T_tubeout values that will correspond 
            # to the exchanger with equivalent P
            # For this we use R, pequivalent, T_shellin and T_shellout and find equivalent T_tubein and equivalent T_tubeout
            equivalenttemps = FCorrectionFactor._get_t1_t2_for_given_P_R (T1=T_shellin, T2=T_shellout, P=pequivalent, R=R)
            t1equivalent = equivalenttemps[2]
            t2equivalent = equivalenttemps[3]
            # Now cal 1-2 exachanger with these equivalent temps
            fn_2n = FCorrectionFactor.oneshell2ntubepasses(T_tubein=t1equivalent,
                                                    T_tubeout=t2equivalent,
                                                    T_shellin=T_shellin,
                                                    T_shellout=T_shellout)
            return fn_2n
        except:
            return "Impossible"

    def _oneshell3tube(self, T_tubein=None, T_tubeout=None, T_shellin=None, T_shellout=None):
        """Based on equation by Fischer
        "Mean Temperature Difference Correction in Multipass Exchangers
        Industrial and Engineering Chemistry, 1938
        vol 30, No 4, page 377-383
        
        VERIFY AGAIN BEFORE ADDING TO PUBLIC API
        
        """
        try:
            R, P = FCorrectionFactor.calc_R_P(self, T_tubein=T_tubein,
                                              T_tubeout=T_tubeout,
                                              T_shellin=T_shellin,
                                              T_shellout=T_shellout)
            R1 = 1/R
            T1 = T_shellin
            T2 = T_shellout
            t1 = T_tubein
            t2 = T_tubeout
            
            if R1 == 1:
                LMTD_countercurrent = abs(T1-t2)
            if R1 !=1:
                LMTD_countercurrent = ((T1-t2)-(T2-t1)) / math.log((T1-t2)/(T2-t1))
           
            def solvefor_thetam(thetaM):
                # thetaM is the true mean temperature from which 'F' can be computed
                # F = thetaM / LMTD_countercurrent
                if R1 != 1:
                    lambda_1 = math.sqrt(9 - 4*R1*(1-R1))
                    m1L =1/2 * (3+lambda_1) * (T1-T2) / (3*thetaM)
                    m2L = 1/2 * (3-lambda_1) * (T1-T2) / (3*thetaM)
                    denom_1 = math.exp(m1L) - math.exp(m2L)
                    A = ((T1-t2)- math.exp(m2L)*(T2-t1)) / denom_1/(1-R1)
                    B = (math.exp(m1L)*(T2-t1) - (T1-t2)) / denom_1/(1-R1)
                    LHS = A*m1L + B*m2L
                    RHS_1 = (T1-T2)/thetaM * (T2 - t1 + math.exp((t2-t1)/3/thetaM) * (T1-t2))
                    RHS_2 = math.exp((t2-t1)/3/thetaM) * (A*m1L*math.exp(m1L) + B*m2L*math.exp(m2L))
                    RHS = RHS_1 - RHS_2
                    return LHS-RHS
                if R1 == 1:
                    # This case was derived by Prof. Harvinder Singh Gill
                    # Equation 12 in the paper by Fischer was adapted by putting R1=1
                    # and the entire solution was reworked for this case
                    # Fischer provided data for R1=1 in the form of a Table (Table - I)
                    # However, Fischer did not show the work nor the equation that
                    # he had used to generate Table I for the case R = 1
                    bLa = 1/9 * (t2-t1) * (t2-T1) / thetaM
                    eal = math.exp((T1-T2)/thetaM)
                    B = (T1-T2+bLa) / (eal-1)
                    A = T2-B
                    termA = math.exp((t2-t1)/3/thetaM) * (3*t2 - 3*T1 - (t2-T1)/3 + 3*B*eal)
                    LHS = (T1-T2) * (3*T2 - 3*t1 - termA)
                    RHS = 3*(T1-T2)*B - 1/3*(t2-t1)*(t2-T1)
                    return LHS-RHS
            ans = fsolve(solvefor_thetam, LMTD_countercurrent)
            F = ans[0]/LMTD_countercurrent
            return F
        except:
            return "Impossible"


    def _nshell3ntube(self, T_tubein=None, T_tubeout=None, T_shellin=None,
                     T_shellout=None, shellpass_count=None):
        """ Use when shell passes >= 2
            and tube passes = multiples of '3' of shell pass
            example : 2 shell-6 tube,
                      3 shell-9 tube,
                      4 shell-12 tube passes
                      ... and so on
                      
            VERIFY AGAIN BEFORE ADDING TO PUBLIC API
            
        """       
        try:
            n = shellpass_count
            R, P = FCorrectionFactor.calc_R_P(self, T_tubein=T_tubein, T_tubeout=T_tubeout, T_shellin=T_shellin, T_shellout=T_shellout)
            R1 = 1/R
            T1 = T_shellin
            T2 = T_shellout
            t1 = T_tubein
            t2 = T_tubeout
            # First we convert P for the 'n' shell '3N' tube exchanger
            # into 'P' that is equivalent to '1' shell '3' tube exchanger           
            psolved = fsolve(FCorrectionFactor._find_P_equivalent, 0.5, R, P, n)
            pequivalent = psolved[0]
            print("equivalent P =", pequivalent)
            
            # Now we take this P-equivalent and use 1 shell - 3 tube formula to find 'F'
            # To do this we find equivalent T_tubein and T_tubeout values that will correspond 
            # to the exchanger with equivalent P
            # For this we use R, pequivalent, T_shellin and T_shellout and find equivalent T_tubein and equivalent T_tubeout
            equivalenttemps = FCorrectionFactor._get_t1_t2_for_given_P_R (T1=T_shellin, T2=T_shellout, P=pequivalent, R=R)
            t1equivalent  = equivalenttemps[2]
            t2equivalent  = equivalenttemps[3]
            # Now cal 1-3 exachanger with these equivalent temps
            fn_2n = FCorrectionFactor._oneshell3tube(T_tubein=t1equivalent, T_tubeout=t2equivalent , T_shellin=T_shellin, T_shellout=T_shellout)
            return fn_2n
        except:
            return "Impossible"

    def _get_t1_t2_for_given_P_R (self, T1=1, T2=1, P=1, R=1):
        t1 = T1-(T1-T2)/P/R
        t2 = P*(T1-t1)+t1
        Rcalc = (T1-T2)/(t2-t1)
        Pcalc = (t2-t1)/(T1-t1)
        return T1,T2,t1,t2,Rcalc,Pcalc


    
    
    def _twopass_crossflow_shellmixed_tubeunmixed(self, T_tubein=None, T_tubeout=None,
                                                 T_shellin=None, T_shellout=None):
        # this is equation 13 from Bowman, 1940
        T1 = T_shellin
        T2 = T_shellout
        t1 = T_tubein
        t2 = T_tubeout
        P = (T1-T2)/(T1-t1)
        Q = (t2-t1)/(T1-t1)
        if P!=Q:
            r0 = (P-Q)/math.log((1-Q)/(1-P))
        if P==Q:
            r0 = (T1-t2)/(T1-t1)
        ta = math.sqrt((1-Q)/(1-P))
        qOverP = Q/P
        tb = 1 - qOverP * math.log((ta-qOverP)/(1-qOverP))
        R = Q/(2* math.log((1)/(tb)))
        F = R/r0
        return F


class EffNTU(object):
    r""" Contains functions to compute 'effectiveness'
    and NTU for heat exchangers.


    Parameters
    ----------
    Cmin : `int or float`
        Minimum heat capacity rate of fluids in heat exchanger.
        :math:`C_{min} = \dot{m} c_p` =  smaller of the two heat 
        capacity rates for the two fluids in the heat exchanger.
            
            :math:`\dot{m}` = mass flow rate of the fluid
            
            :math:`c_p` = specific heat of the fluid
    Cmax : `int or float`
        Maximum heat capacity rate of fluids in heat exchanger.
        :math:`C_{max} = \dot{m} c_p` =  Larger of the two heat l
        capacity rates for the two fluids in the heat exchanger
            
            :math:`\dot{m}` = mass flow rate of the fluid
            
            :math:`c_p` = specific heat of the fluid
    Cratio : `int or float`
        Cmin/Cmax
    effectiveness : `int or float or str`
        Efficiency of the heat exchanger.    
        Input 'int or float' value if effectiveness values is known, or 
        input "?" if it has to be computed. Default = "?".  
    NTU : `int or float or str`
        Number of transfer units for the heat exchanger. 
        Input 'int or float' value if NTU values is known, or 
        input "?" if it has to be computed. Default = "?".
    
    
    Notes
    -----
    1. At least one of effectiveness or NTU must be entered 
    as 'int or float' and the other will get computed.
    
    2. Textbooks contain explicit formulas to compute
    i) effectiveness and ii) NTU for different exchanger types.
    
    3. In this 'class' definition, the following exchanger types
    can be computed with respect to 'effectiveness' or 'NTU''
        
    - Double pipe - parallel flow
    
    - Double pipe - counter flow
    
    - Shell and tube
    
        * one shell pass and 2, 4, 6 etc tube passes
        * n shell passes and 2n, 4n, etc tube passes
        
    - Single pass counter flow:
        
        * both fluids unmixed
        * Cmax fluid mixed, and Cmin unmixed
        * Cmin fluid mixed, and Cmax unmixed
        
    - All exhangers when Cratio = Cmin/Cmax = 0
        
    4. Logic:
        
        Explicit formula for computing 'effectiveness' for each of the above
        exchanger types is implemented in a separate method. 
        When that method is called, it internally calls another
        method ("_get_eff_or_ntu"). This method checks which of
        "effectiveness" or "NTU" is missing. If "effectiveness" is
        missing, and "NTU" is provided, the explict equation is used
        to compute effectiveness from "NTU". If "NTU" is missing, the same
        explict equation for effectiveness becomes an implicit equation for "NTU",
        and it is solved using 'fsolve' of scipy to find "NTU". A
        seperate explict equation for calculating "NTU" from "effectiveness"
        is not implemented. The advantage of this approach is that for the
        case of "Single pass counter flow with both fluid unmixed", there
        is no explict relationship for "NTU" reported in textbooks, but
        there is one for computing "effectiveness". So in this manner by using
        'fsolve', the "NTU"" for the case of "Single pass counter flow
        with both fluid unmixed" can also be computed.
    
    Attributes
    ----------
    See "Parameters". All parameters are attributes.
    

    Examples
    --------
    First import the module **heatexchangers**.
       
    >>> from pychemengg.heattransfer import heatexchangers as hx 
    >>> hx1 = hx.EffNTU(Cmin=0.639, Cmax=0.836, NTU=0.854)
    # This will create an instance of EffNTU and assign it to "hx1".
    # Methods of the class 'EffNTU' can then be called.
    # In this particular example, NTU is assigned = 0.854.
    # By default effectiveness = "?".
    # When a heat exchanger type is next called, it will automatically
    # compute the effectiveness for that exchanger type
    # For example
    >>> hx1.shelltube(shellpass_count=1, tubepass_count=8)
    0.4621565341983249
    # This is the effectiveness.

    """

    def __init__(self, Cmin=None, Cmax=None, NTU="?", effectiveness="?"):
        # C: heat capacity rate = specificheat * massrate
        self.Cmin = Cmin
        self.Cmax = Cmax
        self.Cratio = Cmin/Cmax
        self.effectiveness = effectiveness
        self.NTU = NTU

    def _get_eff_or_ntu(self, func):
        r""" Determines which of "effectiveness" or "NTU" is provided and
        computes the unknown quantity.
        
        
        Parameters
        ----------
        func : `function object`
            function of the form: effectiveness = func(NTU)
        
        
        Returns
        -------
        effectiveness or NTU of heat exchanger : both as `int or float`
            
            
        Notes
        -----
        1. Function - 'func' is passed to this method.
        2. A check is made whether "effectiveness" or "NTU" is unknown.
        3. If "effectiveness" is unknown, the function object, "func"
        is used to compute "effectiveness" from known "NTU" because
        func is explicit for "effectiveness" in the form :
            effectiveness = func(NTU).
        4. If "NTU" is unknown, the same function object, "func" is
        used to solve for "NTU" with known effectiveness by setting
        the equation: func(x) - effectiveness = 0.
        """
        if self.effectiveness == "?" and self.NTU != "?":
            # Case: effectiveness is unknown but NTU is known
            self.effectiveness = func(self.NTU) # computes effectiveness from func
            return self.effectiveness
        if self.NTU == "?" and self.effectiveness != "?":
            # Case: NTU is unknown but effectiveness is known
            guess = 1
            NTUValue = fsolve(lambda x: func(x) - self.effectiveness, guess)
            # solve for NTU using fsolve solver of scipy
            self.NTU = NTUValue[0]
            return self.NTU
        

    def doublepipe_parallelflow(self):
        r""" Computes "effectiveness" or "NTU" for double pipe parallel flow
        heat exchanger.
        """
        def effectiveness(NTU):
            num = (1 - math.exp( - NTU*(1 + self.Cratio)))
            denom = 1 + self.Cratio
            effvalue = num/denom
            return effvalue
        return self._get_eff_or_ntu(effectiveness)
    
    
    def doublepipe_counterflow(self):
        r""" Computes "effectiveness" or "NTU" for double pipe counter flow
        heat exchanger.
        """
        def effectiveness(NTU):
            if self.Cratio < 1:
                num = (1 - math.exp( - NTU*(1 - self.Cratio)))
                denom = (1 - self.Cratio*math.exp( - NTU*(1 - self.Cratio)))
                effvalue = num/denom
                return effvalue
            elif self.Cratio ==1:
                effvalue = NTU/(1+self.NTU)
                return effvalue
        return self._get_eff_or_ntu(effectiveness)
    
    
    def shelltube(self, shellpass_count = None, tubepass_count = None):
        r""" Computes "effectiveness" or "NTU" for shell and tube
        heat exchanger.
        """
        def effectiveness(NTU):
            if tubepass_count%2 == 0:
                self.shellpass_count = shellpass_count
                self.tubepass_count = tubepass_count
                term1 = math.sqrt(1 + self.Cratio**2)
                term2 = math.exp(- NTU * term1)
                b = 2 * (1 + self.Cratio + term1 * (1+term2) / (1-term2))**(-1)
                if self.shellpass_count == 1:
                    effvalue = b
                    return effvalue
                if self.shellpass_count > 1:
                    n = self.shellpass_count
                    if (self.tubepass_count%self.shellpass_count)%2 == 0:
                        c = (1 - b*self.Cratio)/(1 - b)
                        d = (c**n - 1)/(c**n - self.Cratio)
                        effvalue = d
                        return effvalue
                    else:
                        print("Tube pass must be a multiple of shell pass")
                        print("Tube passes must be in multiples of '2'")
                        return None
            else:
                print("Tube passes must be in multiples of '2'")
                return None
        return self._get_eff_or_ntu(effectiveness)


    def crossflow_bothfluids_unmixed(self):
        r""" Computes "effectiveness" or "NTU" for cross flow heat exchanger
        with both fluid unmixed.
        """
        def effectiveness(NTU):
            effvalue = 1-math.exp((1/self.Cratio) 
                                    * NTU**0.22
                                    * (math.exp(- self.Cratio * NTU**0.78)-1))
            return effvalue
        return self._get_eff_or_ntu(effectiveness)            


    def crossflow_Cmin_unmixed(self):
        r""" Computes "effectiveness" or "NTU" for cross flow heat exchanger
        with both fluid with Cmin unmixed.
        """
        def effectiveness(NTU):
            effvalue = (1/self.Cratio) * ( 1-math.exp(- self.Cratio*(1-math.exp(-NTU))))
            return effvalue
        return self._get_eff_or_ntu(effectiveness) 
    
    
    def crossflow_Cmax_unmixed(self):
        r""" Computes "effectiveness" or "NTU" for cross flow heat exchanger
        with both fluid with Cmax unmixed.
        """        
        def effectiveness(NTU):
            effvalue = 1-math.exp(- 1/self.Cratio*(1-math.exp(-self.Cratio*NTU)))
            return effvalue
        return self._get_eff_or_ntu(effectiveness)
    
    
    def Cratio_is_zero(self):
        r""" Computes "effectiveness" or "NTU" for any exchanger with
        Cmin/Cmax = 0
        """
        def effectiveness(NTU):
            effvalue = 1-math.exp(-NTU)
            return effvalue
        return self._get_eff_or_ntu(effectiveness)
    
def calc_overallheattransfercoefficient_fromNTU(NTU=None, area=None, Cmin=None):
    r"""To compute overall heat transfer coefficient (U) from NTU
    
    Parameters
    ----------
    NTU : `int or float`
        Number of transfer units of heat exchanger.
    area : `int or float`
        Surface area of heat transfer for exchanger.
    Cmin: `int or float`
        Minimum heat capacity rate.
        

    Returns
    -------
    Overall heat transfer coefficient (U) : `int or float`
        
    
    
    Notes
    -----
    The following formula is used.
    
    .. math::
        U = \frac {NTU \hspace{2pt} C_{min}} {A_{surfacearea}}
        
           
    *where:*
        
        NTU = number of transfer units for the heat exchanger
        
        :math:`C_{min} = \dot{m} c_p` =  smaller of the two
        capacity rates for the two fluids in the heat exchanger
        
        
        :math:`\dot{m}` = mass flow rate of the fluid
        
        :math:`c_p` = specific heat of the fluid
        
        :math:`A_{surfacearea}` : Heat transfer surface area of the
        exchanger
        
        U = overall heat transfer coefficient
    
    
        
    See Also
    --------
    pychemengg.heattransfer.heatcommonmethods.calc_overallheattransfercoefficient
    
    
    Examples
    --------
    First import the module **heatexchangers**.
       
    >>> from pychemengg.heattransfer import heatexchangers as hx 
    >>> hx.calc_overallheattransfercoefficient_fromNTU(NTU=0.651, Cmin=5.02e3, area=5.11)
    639.5342465753424

    
    References
    ----------
    [1] Yunus A. Cengel and Afshin J. Ghajar, "Heat And Mass Transfer
    Fundamentals and Applications", 6th Edition. New York, McGraw Hill
    Education, 2020.
    """
    return NTU * Cmin / area
