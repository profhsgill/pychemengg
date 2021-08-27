.. _fins:

Fins
=====

Introduction
-------------

The module ``fins`` contains functions that can be used to find the

1. fin-efficiency, and
2. surface area 

of different fin geometries.

The geometries available are:

1. cylindrical (fin with constant circular cross section)
2. rectangular (fin with constant rectangular cross section)
3. rectangularannular (circular fin around a cylindrical shape)
4. pintriangular (concial fin)
5. straightparabolic (rectangular base with parabolic cross-section profile)

.. note::
    The formulas for fin efficiency assume an **adiabatic** fin tip case. If the
    fin tip is not adiabatic, the appropriate corrected fin length (Lc) must be used
    for the "length" keyword.
    
    The corrected lengths (Lc) and corrected outer radius (roc) for the different geometries are:

    .. csv-table:: Corrected fin lengths
        :header: "Geometry", "Corrected length"
        :widths: 20, 20

        "rectangular", "Lc = L + t/2 (L = physical length, t = thickness, Lc = corrected length)"
        "cylindrical", "Lc = L + D/4 (L = physical length, D = diameter, , Lc = corrected length)"
        "rectangularannular", "roc = ro + t/2 (ro = outer radius, t = thickness, , roc = corrected outer radius)"
        "pintriangular", "Lc = L (no change, L = physical length, Lc = corrected length)"
        "straightparabolic", "Lc = L (no change, L = physical length, Lc = corrected length)"


.. container:: custom

    **How to use**
    
    It is recommended that the module be imported
    as ``from pychemengg import fins as fins``

The following examples demonstrate how the module ``fins``
can be used to solve heat transfer problems.

Examples
--------

Example 1: Fins ``cylindrical, single fin``
............................................

``Example 1.`` A rod of diameter D = 2 cm, length= 25 cm, and thermal
conductivity k = 50 W/mC is exposed to ambient temperature air at
T\ :sub:`infinity`\  = 20 C with a heat transfer coefficient h = 64 W/m2C.
If one end of the rod is maintained at 120 C, find the heat loss from the
rod. ``Ans: heat loss = 25.1 kW``

    .. code-block:: python

        # EXAMPLE 1
        from pychemengg.heattransfer import fins as fins
        # The rod is similar to a single fin.
        # Model it as a cylindrical fin.
        # The rod tip is assumed to be adiabatic, so length = actual length.
        rod = fins.Fin(length=25e-2, diameter=2e-2, heattransfercoefficient=64, thermalconductivity=50)
        # There are other parameters for the 'Fin' class, however, only the relevant ones
        # need to be provided. The rest can be left as 'None' by default.
        # Next call on the cylindrical method to get efficiency and surface area.
        rod.cylindrical()
        # The above method will run and create attributes 'efficiency' and 'surfacearea'.
        max_heatrate = rod.heattransfercoefficient * rod.surfacearea * (120-20)
        actual_heatrate = max_heatrate * rod.efficiency
        print(f"Heat loss = {actual_heatrate: 0.1f} W")

        # PRINTED OUTPUT
        Heat loss =  25.1 W

====

Example 2: Fins ``rectangular annular, many fins``
..................................................

``Example 2.`` 
Annular fins with inner radius of 12.5 mm, outer radius of 22.5 mm and thickness
of 1 mm are attached to a cylinder at 200 fins per meter. The cylinder surface is
maintained at 250°C. The surrounding air temperature is 25°C. Consider the convective
heat transfer coefficient (h) as 25 W/m2K and the thermal conductivity (k) of the fin
as 240 W/mK. 

1. Find heat loss per fin. ``Ans: heat loss = 13 W``

2. Heat loss per unit length of cylinder. ``Ans: heat loss = 2961 W/m``

    .. code-block:: python

        # EXAMPLE 2
        from pychemengg.heattransfer import fins as fins
        # Model fin as a rectangularannular fin.
        # The fin tip is not given to be adiabatic.
        # So outer radius should be corrected as correctedradius = outerfinradius + thickness/2
        correctedouter_radius = 22.5e-3 + 1e-3/2
        annularfin = fins.Fin(inner_radius=12.5e-3, outer_radius=correctedouter_radius, thickness=1e-3, heattransfercoefficient=25, thermalconductivity=240)
        # There are other parameters for the 'Fin' class, however, only the relevant ones
        # need to be provided. The rest can be left as 'None' by default.
        # Even if irrelevant ones are provided, the method only uses the relavant parameters.
        # Next call on the rectangularannular method to get efficiency and surface area.
        annularfin.rectangularannular()
        # The above method will run and create attributes 'efficiency' and 'surfacearea'.
        max_heatrateperfin = annularfin.heattransfercoefficient * annularfin.surfacearea * (250-25)
        actualheatrateperfin = max_heatrateperfin * annularfin.efficiency
        print(f"Heat loss from single fin = {actualheatrateperfin: 0.1f} W")
        # To compute heat loss from 1 m long finned cylinder
        # Total heat loss from 1m cylinder = heat loss from fins on 1m + heat loss from surface not covered by fins
        fincount=200
        total_finheatloss = fincount * actualheatrateperfin
        import math
        barecylindersurfacearea = math.pi * (12.5e-3 * 2) * 1 # area = pi * D * L
        baseareaof_singlefin = math.pi * (12.5e-3 *2) * 1e-3 # here L = 1e-3
        areacoveredby_allfins = 200 * baseareaof_singlefin
        heatlossfrom_area_without_fins = 25 * (barecylindersurfacearea-areacoveredby_allfins) * (250-25)
        totalheatloss = total_finheatloss + heatlossfrom_area_without_fins
        print(f"Heat loss for 1 m cylinder = {totalheatloss: 0.0f} W/m")

        # PRINTED OUTPUT
        Heat loss from single fin =  13.0 W
        Heat loss for 1 m cylinder =  2961 W/m

