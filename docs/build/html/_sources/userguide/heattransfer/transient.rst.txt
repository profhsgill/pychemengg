.. _transient:

Transient Heat Transfer
===========================

Introduction
-------------

The module ``transient`` contains functions that can be used to solve
unsteady steady state heat transfer problems.

The following three cases are covered in this module:

1. Lumped system
2. Nonlumped system

    * Large wall

    * Long Cylinder

    * Sphere

3. Semi infinite solids
    
    The following cases of boundary conditions are covered:

        * Surface temperature is specified

        * Heatflux is specified

        * Surface convection is specified

        * Energy pulse is specified


.. container:: custom

    **How to use**
    
    It is recommended that the module be imported
    as ``from pychemengg import transient as transient``.

The following examples demonstrate how the module ``transient``
can be used to solve heat transfer problems.

Examples
--------

Example 1: Lumped system ``rectangular``
................................................

``Example 1.`` Plates (1 cm thick) made of metal are heated by passing them
through an enclosed space maintained at 800 C. The plates enter
the oven at 20 C and stay in the oven for 2 min. Find the temperature
of the plates as they exit the enclosure. Use the heat transfer
coefficient as 200 W/m2K and the following properties for the metal plates.

Metal plates:

    k = 180 W/mK
    
    :math:`\rho` = 2800 kg/m3
    
    :math:`c_p` = 880 J/kg K
    
    ``Ans: Temperature of plates at exit = 689 C``

    .. code-block:: python

        # EXAMPLE 1
        from pychemengg.heattransfer import transient as transient
        # Start by assuming lumped system analysis can be applied
        # Create an instance of LumpedSystem
        # Since plate area is not given consider one side area = '1'
        # Since heat will transfer from both sides, surface area = 2*1
        # Volume then equals = Area of one side*thickness
        plate = transient.LumpedSystem(surfacearea=2, volume=1*1e-2, density=2800, specificheat=880, thermalconductivity=180, heattransfercoefficient=200, T_infinity=800, T_initial=20)
        # Compute Biot number to check lumped system assumption
        biotnumber = plate.calc_Bi()
        print(f"Biot number = {biotnumber}")
        # The result is Biot number = 0.005555555555555556 which is < 0.1
        # Therefore lumped system assumption is valid
        # Temperature of plates as they exit can be computed as temperature at time = 2 min
        exittemperature = plate.calc_temperature_of_solid_at_time_t(time=2*60)
        print(f"Temperature of plates at exit = {exittemperature: 0.0f} C")

        # PRINTED OUTPUT
        Temperature of plates at exit =  689 C

====

Example 2: Lumped system ``cylindrical``
.............................................

``Example 2.`` A long copper rod of diameter 2 cm is exposed to air at 20 C
with a heat transfer coefficient of 200 W/m2K. If the initial temperature
of the rod is 100 C, how long will it take for the rod to reach an average
temperature of 20 C. Take the following properties for copper:

Copper:

    k = 401 W/mK
    
    :math:`\rho` = 8933 kg/m3
    
    :math:`c_p` = 385 J/kg K
    
``Ans: Time = 238 s``

.. code-block:: python

    # EXAMPLE 2
    from pychemengg.heattransfer import transient as transient
    # Start by assuming lumped system analysis can be applied
    # Create an instance of LumpedSystem
    # All units are in SI system
    import math
    diameter = 2e-2
    # Assume length = 1
    length = 1
    area = math.pi * diameter * length
    volume = math.pi/4 *diameter**2 * length
    rod = transient.LumpedSystem(surfacearea=area, volume=volume, density=8933, specificheat=385, thermalconductivity=401, heattransfercoefficient=200, T_infinity=20, T_initial=100)
    # Compute Biot number to check lumped system assumption
    biotnumber = rod.calc_Bi()
    print(f"Biot number = {biotnumber}")
    # The result is Biot number = 0.0024937655860349127 which is < 0.1
    # Therefore lumped system assumption is valid
    # Time for rod to reach a temperature of 20 C can be found as follows
    # Set a function that accepts 'time' and returns 'temperature' equation
    
    def func (time):
        temp = rod.calc_temperature_of_solid_at_time_t(time=time[0])
        # This temp is desired to be = 20 C
        equation = temp-25
        return equation

    # Now solve this function using a solver
    from scipy.optimize import fsolve 
    # Typically these import statements are placed at top of code.
    guess_time = 1
    # Guess is required to solve, and this is a random value.
    # User can change it and the result should be the same
    solution = fsolve(func, guess_time)
    timetaken = solution[0]
    # Because output of fsolve is an array, use [0] to get the value.
    print(f"Time = {timetaken: 0.0f} s")

    # PRINTED OUTPUT
    Time =  238 s

====

Example 3: Non lumped system ``cylindrical``
...................................................

``Example 3.`` A cylindrical wooden log measuring 10 cm in diameter is at a uniform
temperature of 15 C. It is exposed to hot gases at 550 C in a fireplace with a
heat transfer coefficient of 13.6 W/m2K. If the ingnition temperature is 420 C find the
time needed for the log to ignite i) using one-term approximation, ii)
using 10 terms of series solution. Use the following properties for the wooden log:

    k = 0.17 W/mK
 
    :math:`\alpha` = 1.28e-7 m2/s
    
``Ans: i) Time = 2771 s; Time = 2966 s``

.. code-block:: python

    # EXAMPLE 3
    from pychemengg.heattransfer import transient as transient
    # Problem statement asks that non lumped model be used.
    # Create an instance of non lumped system.
    # All units are in SI system.
    import math
    diameter = 10e-2
    # Assume length = 1
    length = 1
    area = math.pi * diameter * length
    volume = math.pi/4 *diameter**2 * length
    wood = transient.NonLumpedCylinder(radius=diameter/2, surfacearea=area, volume=volume, density=None, specificheat=None, thermalconductivity=0.17, thermaldiffusivity=1.28e-7, heattransfercoefficient=13.6, T_infinity=550, T_initial=15)
    biotnumber = wood.calc_Bi()
    print(f"Biot number = {biotnumber}")
    # The following gets printed to screen.
    # Biot number = 4.0
    # Case 1: use of one term approximation.
    # Set up a function that accepts 'time' and returns
    # an equation that can be solved using a solver.

    def func(time, *number_of_terms):
        # Here 'time' is guess values
        # and 'number_of_terms' are the terms to consider for infinite series solution
        fouriernumber =  wood.calc_Fo(time=time[0])
        # calculate first eigen value
        eigenvalue = wood.calc_eigenvalues(numberof_eigenvalues_desired=number_of_terms[0])
        temperature = wood.calc_temperature_of_solid_at_time_t(rposition_tofindtemp=wood.radius)
        # at some time 't', this temperature should be = 420
        # Therefore equation becomes, equation = temperature-420 = 0
        equation = temperature-420
        return equation

    # Now solve this function using a solver
    from scipy.optimize import fsolve 
    # Typically these import statements are placed at top of code.
    guess_time = 1
    # Guess is required to solve, and this is a random value.
    # User can change it and the result should be the same
    solution = fsolve(func, guess_time, (1,))
    timetaken_1 = solution[0]
    # Because output of fsolve is an array, use [0] to get the value.
    print(f"Time = {timetaken_1: 0.0f} s")
    # The following gets printed to screen
    # Time =  2771 s
    # Check the Fourier number
    print(f"Fourier number = {wood.calc_Fo(time=timetaken_1)}")
    # This prints the following to screen
    # Fourier number = 0.1418840952127394
    # Fourier number is not greater than 0.2, thus one term approximation
    # is not valid.

    # Repeat the above calculation but this time use 10 terms of the series solution


    # Now solve with 10 terms
    solution = fsolve(func, guess_time, (10,))
    timetaken_10 = solution[0]
    # Because output of fsolve is an array, use [0] to get the value.
    print(f"Time = {timetaken_10: 0.0f} s")
    # The following gets printed to screen
    # Time =  2966 s

    # Now print the solutions for both cases:
    print(f"Case 1: 1 term approximation gives time = {timetaken_1: 0.0f} s as answer")
    print(f"Case 2: 10 term approximation gives time = {timetaken_10: 0.0f} s as answer")

    # PRINTED OUTPUT
    Case 1: 1 term approximation gives time =  2771 s as answer
    Case 2: 10 term approximation gives time =  2966 s as answer


====

Example 4: Semi infinite ``boundary condition: surface temperature is specified``
...................................................................................

``Example 4.`` Tennis court made of clay is initially at a uniform tempertaure
of 55 C. Suddenly due to rain the surface temperature becomes 25 C.
Assume temperature of clay surface is maintained at 25 C. What is the
temperature of the tennis court 3 cm deep and the heat flux that has been
transferred after 60 min. Use the following properties for clay:

    k = 0.062 W/mK
 
    :math:`\rho` = 2115 kg/m3

    :math:`c_p` = 920 J/kg K
    
``Ans: Temperature = 53.6 C, Flux = 98 W/m2``

.. code-block:: python

    # EXAMPLE 4
    from pychemengg.heattransfer import transient as transient
    # All units are in SI system.
    # Create an instance of semi infinite system.
    # Only relevant keywords are to be input.
    clay = transient.SemiInfinite(boundarycondition="surfacetemperature_specified", xposition_tofindtemp=3e-2, time=60*60, density=2115, specificheat=920, thermalconductivity=0.062, thermaldiffusivity=None, constantsurfacetemperature=25, heattransfercoefficient=None, heatflux=None, energypulse=None, T_infinity=None, T_initial=55)
    # Call the method calc_temperature.
    # Depending on the boundary condition, the method selects the
    # appropriate equation to compute the temperature at given 'x' and 'time'
    temperature = clay.calc_temperature()
    # Call the method calc_heatflux_forconstantsurfacetemperature()
    flux = clay.calc_heatflux_forconstantsurfacetemperature()
    # Print the results.
    print(f"Temperature at 3cm depth in clay = {temperature: 0.1f} C")
    print(f"Flux in 60 min = {flux: 0.0f} W/m2. Negative sign indicates loss of heat.")
    
    # PRINTED OUTPUT
    Temperature at 3 cm depth in clay =  53.6 C
    Flux in 60 min = -98 W/m2. Negative sign indicates loss of heat.


====

Example 5: Semi infinite ``boundary condition: heat flux is specified``
........................................................................

``Example 5.`` In a room designed to test experimental fuels, a particular fuel
is being burnt and as a result the brick walls of the room are subjected
to a constant flux of 20,000 W/m2. If the initial wall temperature was 15 C,
find the temperature in the wall at 1 cm depth after 1 h of exposure to
constant flux. Take the following properties for the wall:

    k = 1 W/mK
 
    :math:`\alpha` = 5.08e-7 m2/s
    
``Ans: Temperature = 793 C``

.. code-block:: python

    # EXAMPLE 5
    from pychemengg.heattransfer import transient as transient
    # All units are in SI system.
    # Create an instance of semi infinite system.
    # Only relevant keywords are to be input.
    wall = transient.SemiInfinite(boundarycondition="heatflux_specified", xposition_tofindtemp=1e-2, time=60*60, density=None, specificheat=None, thermalconductivity=1.0, thermaldiffusivity=5.08e-7, constantsurfacetemperature=None, heattransfercoefficient=None, heatflux=20000, energypulse=None, T_infinity=None, T_initial=15)
    # Call the method calc_temperature.
    # Depending on the boundary condition, the method selects the
    # appropriate equation to compute the temperature at given 'x' and 'time'
    temperature = wall.calc_temperature()
    # Print the results.
    print(f"Temperature at 1 cm depth in wall = {temperature: 0.0f} C")
   
    # PRINTED OUTPUT
    Temperature at 1 cm depth in wall =  793 C


====

Example 6: Semi infinite ``boundary condition: surface convection is specified``
..................................................................................

``Example 6.`` Insulation bricks are being tested. They must not exceed a temperature
of 450 C when exposed to hot gases at a temperature of 550 C for 5 minutes.
If the initial temperature of insulation bricks is 25 C, comment whether the bricks
will reach their safe operational limit of 450 C. Take the following properties
for the bricks:

    k = 0.17 W/mK
 
    :math:`\alpha` = 1.28e-7 m2/s

    heat transfer coefficient = 35 W/m2 K
    
``Ans: Surface temperature after 5 min = 360 C``

.. code-block:: python

    # EXAMPLE 6
    from pychemengg.heattransfer import transient as transient
    # All units are in SI system.
    # Create an instance of semi infinite system.
    # Only relevant keywords are to be input.
    brick = transient.SemiInfinite(boundarycondition="surfaceconvection_specified", xposition_tofindtemp=0, time=5*60, density=None, specificheat=None, thermalconductivity=0.17, thermaldiffusivity=1.28e-7, constantsurfacetemperature=None, heattransfercoefficient=35, heatflux=None, energypulse=None, T_infinity=550, T_initial=25)
    # Call the method calc_temperature.
    # Depending on the boundary condition, the method selects the
    # appropriate equation to compute the temperature at given 'x' and 'time'
    temperature = brick.calc_temperature()
    # Print the results.
    print(f"Temperature at brick surface = {temperature: 0.0f} C")
    # Comment: Since the temperature = 360 C < operational limit of 450 C
    # the conditions are safe.
   
    # PRINTED OUTPUT
    Temperature at brick surface = 360 C

    ====

Example 7: Semi infinite ``boundary condition: energy pulse is specified``
.............................................................................

``Example 7.`` Laser with a certain energy pulse output is hitting a large object.
The object is at a uniform initial temperature of 20 C. After 30 s of expsoure,
the temperature of the object at a depth of 25 mm from surface of laser expsoure
is found to be 130 C. Find i) the amount of energy output from the laser, ii)
and temperature of the object at 25 mm depth after 60 s of exposure.
Take the following properties for the object :

    k = 63.9 W/mK
 
    :math:`\alpha` = 18.8e-6 m2/s

    heat transfer coefficient = 35 W/m2 K
    
``Ans: i) Energy pulse = 2.076e7 J/m2, ii) Temperature at 25 cm @ 60s = 109 C``

.. code-block:: python

    # EXAMPLE 7
    from pychemengg.heattransfer import transient as transient
    # All units are in SI system.
    # Create an instance of semi infinite system.
    # Only relevant keywords are to be input.
    # Since energy pulse is required to create the instance and this is also
    # the quantity that has to be found, a function must be created to solve the problem.
    # The function should return an equation that can be solved with a solver

    def func(epulse):
        # Create instance
        object_1 = transient.SemiInfinite(boundarycondition="energypulse_specified", xposition_tofindtemp=25e-3, time=30, density=None, specificheat=None, thermalconductivity=63.9, thermaldiffusivity=18.8e-6, constantsurfacetemperature=None, heattransfercoefficient=None, heatflux=None, energypulse=epulse[0], T_infinity=None, T_initial=20)
        # call on method to find temperature
        temperature = object_1.calc_temperature()
        # This temperature must be 130 C at 30 s.
        # Thus, equation becomes : temperature - 130 = 0.
        equation = temperature - 130
        # Return this equation
        return equation

    # Now solve this function using a solver
    from scipy.optimize import fsolve 
    # Typically these import statements are placed at top of code.
    guess_epulse = 1e3
    # Guess is required to solve, and this is a random value.
    # User can change it and the result should be the same
    solution = fsolve(func, guess_epulse)
    epulse = solution[0]
    # Because output of fsolve is an array, use [0] to get the value.
    # Temperature at 25 mm deoth and 60 s can be computed as follows:
    object_1 = transient.SemiInfinite(boundarycondition="energypulse_specified", xposition_tofindtemp=25e-3, time=60, density=None, specificheat=None, thermalconductivity=63.9, thermaldiffusivity=18.8e-6, constantsurfacetemperature=None, heattransfercoefficient=None, heatflux=None, energypulse=epulse, T_infinity=None, T_initial=20)
    temperature_60s = object_1.calc_temperature()
    # print results
    print(f"Energy pulse generated by laser = {epulse: 0.3e} W/m2")
    print(f"Temperature at 25 mm depth and 60 s of expsoure = {temperature_60s: 0.0f} C")
 
    # PRINTED OUTPUT
    Energy pulse generated by laser =  2.076e+07 W/m2
    Temperature at 25 mm depth and 60 s of expsoure =  109 C

