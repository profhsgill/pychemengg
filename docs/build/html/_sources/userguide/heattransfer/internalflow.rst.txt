.. _internalflow:

Internal flow
=============

Introduction
-------------

The module ``internalflow`` contains functions related to fluid flow inside circular
pipes/tubes and annular spaces. It contains methods to compute:

1. Nusselt number 
2. Entry length, friction factor, and pressure drop and other quantities

.. container:: custom

    **How to use**
    
    It is recommended that the module be imported
    as ``from pychemengg import internalflow as intflow``

The following examples demonstrate how the module `internalflow``
can be used to solve heat transfer problems.

Examples
--------

Example 1: Internal flow ``circular tube, entry lengths``
.........................................................

``Example 1.`` Engine oil is to be cooled from 100 to 80 C by passing through
a tube maintained at a uniform surface temperature of 40 C. Oil flows
at a velocity of 0.03 ms/. Tube inner diameter is 2 cm. What is the
desired tube length. ``Ans: 2.67 m``

    .. code-block:: python

        # EXAMPLE 1
        from pychemengg.heattransfer import internalflow as intflow
        import math
        # Identify properties at mean buld fluid temp = (100+80)/2.
        # All units are in SI system
        density = 846
        specificheat = 2176
        kinematicviscosity=2.81e-5
        thermalconductivity = 0.138
        viscosity = kinematicviscosity * density
        velocity = 0.03
        diameter = 2e-2
        
        # Model the tube as CircularTube from module internalflow
        tube=intflow.CircularTube()
        # Note, do not use () after CircularTube
        # To use correlations in internalflow, Reynolds number and Prandtl number is needed
        # This is provided in module 'commonmethods'
        from pychemengg.heattransfer import commonmethods as hcm 
        Re = hcm.calc_Re(characteristic_length=2e-2, velocity=velocity, density=density, viscosity=viscosity)
        Pr = hcm.calc_Pr(viscosity=kinematicviscosity*density, specificheat=specificheat, thermalconductivity=thermalconductivity)
        # Print Re to see whether flow is laminar or turbulent
        print(f"Reynolds number = {Re: 0.0f}")
        # Re = 21, which is < 2300 and flow is laminar

        # NOTE:
        # The user needs to examine intermediate results and accordingly
        # select the path forward. 

        # Now check if flow is developed or not.
        # Use correlation for laminar flow.
        thermalentrylength = tube.thermal_entrylength_laminar(Re=Re, Pr=Pr, diameter=diameter)
        print(f"Thermal entry length = {thermalentrylength}")
        hydrolength = tube.hydrodynamic_entrylength_laminar(Re=Re, diameter=diameter)
        print(f"Hydrodynamic entry length = {hydrolength}")

        # hydrolength = 0.021 m, which means flow should be developed hydrodynamically.
        # thermalentrylength = 8 m, and since tube length is not known,
        # it is not feasible to conclude that the flow is thermally developed.
        # Therefore, thermally developing flow relation should be used to find Nu.
        # The following method can be used.
        # Nu_thermallydeveloping_laminar_edwards(Pr=None, length=None, diameter=None).
        # This requires 'length' as input, and this is precisely what has to be
        # computed. This indicates a function should be created that takes 'length'
        # as input and returns some 'equation' that can be solved for 'length'

        # Create function
        # Note that other variables needed such as Pr, diameter are not part of
        # function definition below. Based on Python variable 'scope' the function
        # is able to reach out to access these other variable.
        def func(length_guess):
            Nu = tube.Nu_thermallydeveloping_laminar_edwards(Re=Re, Pr=Pr, length=length_guess, diameter=diameter)
            
            # What equation can be formulated to return via this function
            # so that 'length' can be computed.
            # Enery balance can be used
            # change in internal energy = Energy gained via convection
            # m Cp deltaT_internalenergy = h A delta T_conv
            # These deltaTs are not the same
            # deltaT_internalenergy = In - out
            # deltaT_conv = LMTD
            # Internal energy change function is in module commonmethods
            # hcm.calc_internalenergychange(mass=mass, specificheat=specificheat, deltaT=deltaT)
            areacrosssection = math.pi/4 * diameter**2
            massflowrate = density * velocity * areacrosssection
            rate_deltainternalenergy = hcm.calc_internalenergychange(mass=massflowrate, specificheat=specificheat, deltaT=100-80)
            h = Nu*thermalconductivity/diameter
            deltaT1 = 100-40
            deltaT2 = 80-40
            deltaT_conv = hcm.calc_LMTD(deltaT1=deltaT1, deltaT2=deltaT2)
            heatrate_conv = h * (math.pi * diameter * length_guess) * (deltaT_conv)
            equation = rate_deltainternalenergy - heatrate_conv
            return equation
        from scipy.optimize import fsolve
        guess_length = 1
        solution = fsolve(func, guess_length)
        print(f"Length of tube needed = {solution[0]: 0.2f} m")
        # Notice that L = 2.67 < 8, thus the original assumption
        # of using equation of Nu for thermally developing flow was correct

        # PRINTED OUTPUT
        Length of tube needed =  2.67 m

