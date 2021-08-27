.. _heatexchangers:

Heat Exchangers
================

Introduction
-------------

The module ``heatexchangers`` contains methods to perform heat transfer
calculations related to heat exchangers.

The module has methods for the following:

1. Calculation of F factors
2. Calculation of effectiveness and NTU
3. Calculation of overall heat coefficient from NTU
   
.. container:: custom

    **How to use**
    
    It is recommended that the module be imported
    as ``from pychemengg import heatexchangers as hx``

The following examples demonstrate how the module `heatexchangers``
can be used to solve heat transfer problems.

Examples
--------

Example 1: Counter flow ``LMTD``
.................................

``Example 1.`` Hot fluid (Cp=2.09 kJ/kg-K) flows through a counter flow
heat exchanger at a rate of 0.63 kg/s. It enter at 193 C and leaves
at 65 C. Cold fluid (Cp=1.67 kJ/kgK) exits at 149 C at a rate of
1 kg/s. What area is required if the overall heat transfer
coefficient based on the inside area is 0.7 kW/m2K. ``Ans: Area = 8.5 m2``

    .. code-block:: python

        # EXAMPLE 1
        from pychemengg.heattransfer import heatexchangers as hx
        from pychemengg.heattransfer import heatcommonmethods as hcm

        # Use change in internal energy to compute heat transfered
        heattransfer_hot = hcm.calc_internalenergychange(mass=0.63, specificheat=2.09e3, deltaT=193-65)

        # heattransfer_hot = heattransfer_cold, use this to find inlet
        # temperature of cold fluid

        T_cold_in = 149 - heattransfer_hot/1.67e3

        # compute LMTD
        deltaT1 = 65-T_cold_in
        deltaT2 = 193-149
        LMTD = hcm.calc_LMTD(deltaT1=deltaT1, deltaT2=deltaT2)
        correctionfactor = 1 # Because counterflow
        # compute area using convection equation
        area = heattransfer_hot/0.7e3/LMTD
        print(f"Area required = {area: 0.1f} m2")
        
        # PRINTED OUTPUT
        Area required =  8.5 m2


Example 2: Double-pipe ``Effectiveness-NTU``
............................................

``Example 2.`` A fluid (Cp = 4.18 kJ/kgK) enters a parallel flow, double-pipe heat
exchanger at 40 C at 0.75 kg/s. It is heated by a second fluid (Cp = 1.581 kJ/kgK)
flowing at a rate of 1.5 kg/s with inlet temperature 115 C. If the area is 13 m2
and overall heat transfer coefficient is 0.205 kW/m2K, find
the total heat rate. ``Ans: heat rate = 87.17 kW``

    .. code-block:: python

        # EXAMPLE 2
        from pychemengg.heattransfer import heatexchangers as hx

        # Use effectiveness-NTU method.
        Ccold = 0.75*4.18e3
        Chot = 1.5*1.581e3
        # Print both to see which is Cmin
        print(f" Chot = {Chot} and Ccold = {Ccold}")
        # This produces:  Chot = 2371.5 and Ccold = 3135.0
        # Based on this, assign Cmin and Cmax.
        Cmin = Chot
        Cmax = Ccold
        # Calculate NTU = UA/Cmin
        NTU = .205e3 * 13/Cmin
        # Next given this NTU find effectiveness.
        # First create instance of EffNTU class
        doublepipe = hx.EffNTU(Cmin=Cmin, Cmax=Cmax, NTU=NTU, effectiveness="?")
        # The keyword - 'effectiveness' is assigned the string "?".
        # This alerts the function that effectiveness is to be computed.
        
        # Then call the appropriate exchanger type to compute effectiveness.
        effectiveness = doublepipe.doublepipe_parallelflow()
        
        heatratemax = Cmin*(115-40)
        heatrate_actual = effectiveness * heatratemax
        print(f"Heat rate = {heatrate_actual: 0.3e} W")


        # PRINTED OUTPUT
        Heat rate =  8.719e+04 W


Example 3: Comparison of cross flow and shell and tube ``F correction factors``
................................................................................

``Example 3.`` An exchanger is desired for operation with T\ :sub:`hot,in`\  = 400 C,
T\ :sub:`hot,out`\  = 130 C, and T\ :sub:`cold,in`\  = 25 C. The following
is known about the system:

**Hot fluid**: mass flow rate = 2 kg/s, specific heat = 2000 J/kg k

**Cold fluid**: mass flow rate = 6.857 kg/s, specific heat = 1050 J/kgK

**Overall heat transfer coefficient** = 150 W/m2K

Compare i) shell and tube heat exchanger design with ii) cross flow exchanger
with both fluids unmixed. Find the suitable configurations.

``Ans: 2 Shell - 4 tube: F = 0.93, area = 48.64 m2; Cross flow: F = 0.85, area = 53.8 m2``

    .. code-block:: python

        # EXAMPLE 3
        from pychemengg.heattransfer import heatexchangers as hx
        from pychemengg.heattransfer import heatcommonmethods as hcm

        # Use energy balance to find T_cold_out
        mass_hot = 2
        specificheat_hot = 2000
        T_hot_in = 400
        T_hot_out = 130
        deltaT_hot = T_hot_in - T_hot_out

        hot_internalenergychange = hcm.calc_internalenergychange(mass=mass_hot, specificheat=specificheat_hot, deltaT=deltaT_hot)

        mass_cold = 6.857
        specificheat_cold = 1050
        T_cold_in =  25
        T_cold_out = T_cold_in + hot_internalenergychange/mass_cold/specificheat_cold

        # case: shell-tube
        # Calculate LMTD.
        deltaT_1 = T_hot_in - T_cold_out
        deltaT_2 = T_hot_out - T_cold_in
        LMTD = hcm.calc_LMTD(deltaT1=deltaT_1, deltaT2=deltaT_2)
        print(f"LMTD = {LMTD: 0.1f} C")
        # Calculate F correction factor.
        exchanger = hx.FCorrectionFactor()
        oneshelltube_F_factor = exchanger.oneshell2ntubepasses(T_tubein=T_cold_in, T_tubeout=T_cold_out, T_shellin=T_hot_in, T_shellout=T_hot_out)
        print(f"1-Shell 2-Tube F Factor = {oneshelltube_F_factor: 0.2f}")
        # The F factor is 0.58, which is quite low.
        # Therefore the 2 shell-4 pass configuration can be examined.
        twoshelltube_F_factor = exchanger.twoshell4ntubepasses(T_tubein=T_cold_in, T_tubeout=T_cold_out, T_shellin=T_hot_in, T_shellout=T_hot_out)
        print(f"2-Shell 4-Tube F Factor = {twoshelltube_F_factor: 0.2f}")
        # The F factor for 2 shell and 4 tube is 0.93.
        # Calculate area for this using convection heat rate equation and internal energy change
        # hot_internalenergychange = U * A * LMTD * F
        overallheattransfercoefficient = 150 
        area_shelltube = hot_internalenergychange/overallheattransfercoefficient/LMTD/twoshelltube_F_factor
        print(f"Area needed for 2 shell-4 tube heat exchanger = {area_shelltube:0.2f} m2")

        # case: cross flow both fluids unmixed
        # Calculate F factor.
        crossflow_F_factor = exchanger.singlepass_crossflow_bothunmixed(T_tubein=T_cold_in, T_tubeout=T_cold_out, T_shellin=T_hot_in, T_shellout=T_hot_out)
        print(f"Cross flow F factor= {crossflow_F_factor: 0.2f}")
        # Calculate area as was done for shell-tube case
        area_crossflow = hot_internalenergychange/overallheattransfercoefficient/LMTD/crossflow_F_factor
        print(f"Area needed for cross flow heat exchanger = {area_crossflow:0.2f} m2")



        # PRINTED OUTPUT
        LMTD =  157.5 C
        1-Shell 2-Tube F Factor =  0.58
        2-Shell 4-Tube F Factor =  0.93
        Area needed for 2 shell-4 tube heat exchanger = 49.32 m2
        Cross flow F factor=  0.85
        Area needed for cross flow heat exchanger = 53.70 m2

