.. _externalflow:

External flow
=============

Introduction
-------------

The module ``externalflow`` contains methods related to heat transfer
between a fluid and objects, wherein the fluid is subjected to convective
flow over the external surfaces of the objects.

The module has methods for the following objects:

1. Cylinder
2. Plate
3. Sphere
4. Tube Banks
   
.. container:: custom

    **How to use**
    
    It is recommended that the module be imported
    as ``from pychemengg import externalflow as extflow``

The following examples demonstrate how the module `externalflow``
can be used to solve heat transfer problems.

Examples
--------

Example 1: Tube bank ``heat rate``
....................................................

``Example 1.`` A gas is to be heated from 20 C using a tube bank where tubes are at a
surface temperature of 120 C. The gas enters at 4.5 m/s velocity and 1 atm.
The tubes are arranged in "inline" configuration. Tubes have an outer
diameter of 1.5 cm, and longitudanal and transverse pitches are 5 cm.
There are 6 rows in the direction of gas flow, and there are 10 tubes
per row. Calculate the heat rate per unit length of tubes and pressure
drop in the tube bank. ``Ans: Heat rate=2.56e4 W``

.. note::
    To solve this problem, an iterative procedure is required, where,
    the outlet temperature is assumed, next mean bulk temperature is
    calculated, then gas properties are determined and finally outlet
    temperature is computed, which is then verified against the guess value.
    However, for demonstration purposes, only the final iteration
    is shown with bulk temperature = 25 C

Use the following gas properties:

At 25 C: density = 1.184 kg/m3, specific heat = 1007J/kgK, 
thermal conductivity = 0.02551 W/m K, viscosity = 1.849e-5 kg/ms,

Pr at 120 C = 0.7073

At 20 C (inlet): density = 1.204 kg/m3

    .. code-block:: python

        # EXAMPLE 1
        from pychemengg.heattransfer import externalflow as extflow
        import math

        # Model the system with TubeBank.
        bank = extflow.TubeBank(config="inline", totalrows=6,
                                tubes_per_row=10, transverse_pitch=5e-2,
                                longitudanal_pitch=5e-2, length=1,
                                outer_tubediameter=1.5e-2, T_infinity=20,
                                velocity_infinity=4.5)

        bank.set_fluid_properties(density=1.184, viscosity=1.849e-5, specificheat=1007,
                            thermalconductivity=0.02551, density_surface=None,
                            viscosity_surface=None, specificheat_surface=None,
                            thermalconductivity_surface=None)

        bank.set_Pr_surface(0.7073)

        # Calculate max velocity.
        maxvelocity = bank.calc_maxvelocity()

        #Calculate Reynolds number.
        Re = bank.calc_Re()

        #Calculate Prandtl number
        Pr = bank.calc_Pr()

        #Calculate Nusselt number.
        Nu = bank.calc_Nu()

        # Calculate heat transfer coefficient.
        h = Nu * bank.thermalconductivity / bank.outer_tubediameter

        # Heat transfer using convection formula.
        # q_conv =  h * area * LMTD

        area = (math.pi * bank.outer_tubediameter * bank.length
                * bank.tubes_per_row * bank.totalrows)
        from pychemengg.heattransfer import commonmethods as hcm
        LMTD = hcm.calc_LMTD((120-25), (120-25))
        print(LMTD)

        q_conv = h * area * LMTD
        print(f"Heat transferred = {q_conv: .3e} W/m")

        # Find T_out using : q_internalenergy = q_conv
        # calc_internalenergychange() is available in commonmethods.
        # To find q_internalenergy mass is required.
        # Use mass flow rate = volumetric flowrate * density at T_in
        density_T_in = 1.204
        area_crosssection_at_inlet = bank.tubes_per_row * bank.transverse_pitch
        massin =  density_T_in * bank.velocity_infinity * area_crosssection_at_inlet
        T_out = 20 + q_conv/massin/bank.specificheat
        # T_out turns to be = 29.4 C
        # T_bulk  = (20 + 29.4)/2 = 24.7, which is close to 25 C

        # PRINTED OUTPUT
        Heat transferred =  2.563e+04 W/m

