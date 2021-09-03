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

Example 1: Lumped system ``rectangular plate``
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

Example 2: Lumped system ``cylindrical rod``
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

Example 3: Non lumped system ``rectangular plate``
...................................................

``Example 3.`` A cylindrical wooden log measuring 10 cm in diameter is at a uniform
temperature of 15 C. It is exposed to hot gases at 550 C in a fireplace with a
heat transfer coefficient of 13.6 W/m2K. If the ingnition temperature is 420 C find the
time needed for the log to ignite i) using one-term approximation, ii)
using 10 terms of series solution. Use the following properties for the wooden log:

    k = 0.17 W/mK
 
    :math:`\alpha` = 1.28e-7 m2/s
    
``Ans: Time = 2771 s``

.. code-block:: python

    # EXAMPLE 3
    from pychemengg.heattransfer import transient as transient
    # Problem statement asks that non lumped model be used
    # Create an instance of non lumped system
    # All units are in SI system
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
    # an equation that can be solved using a solver
    def func(time):
        fouriernumber =  wood.calc_Fo(time=time[0])
        # calculate first eigen value
        eigenvalue = wood.calc_eigenvalues(numberof_eigenvalues_desired=1)
        temperature = wood.calc_temperature_of_solid_at_time_t(rposition_tofindtemp=wood.radius)
        # at some time 't', this temperature should be = 420
        equation = temperature-420
        return equation

    # Now solve this function using a solver
    from scipy.optimize import fsolve 
    # Typically these import statements are placed at top of code.
    guess_time = 1
    # Guess is required to solve, and this is a random value.
    # User can change it and the result should be the same
    solution = fsolve(func, guess_time)
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

    # Repeat the above but this time use 10 terms of the series solution
    # Start with redefining the function 'func' so that 10 terms are computed 
    def func(time):
        fouriernumber =  wood.calc_Fo(time=time[0])
        # calculate first 10 eigen values
        eigenvalue = wood.calc_eigenvalues(numberof_eigenvalues_desired=10)
        temperature = wood.calc_temperature_of_solid_at_time_t(rposition_tofindtemp=wood.radius)
        # at some time 't', this temperature should be = 420
        equation = temperature-420
        return equation

    # Now solve this function using a solver
    from scipy.optimize import fsolve 
    # Typically these import statements are placed at top of code.
    guess_time = 1
    # Guess is required to solve, and this is a random value.
    # User can change it and the result should be the same
    solution = fsolve(func, guess_time)
    timetaken_10 = solution[0]
    # Because output of fsolve is an array, use [0] to get the value.
    print(f"Time = {timetaken_10: 0.0f} s")
    # The following gets printed to screen
    # Time =  2966 s

    # Now print the solutions for both cases:
    print(f"Case 1: 1 term approximation gives time = {timetaken_1} s as answer")
    print(f"Case 2: 10 term approximation gives time = {timetaken_10} s as answer")

    # PRINTED OUTPUT
    Case 1: 1 term approximation gives time =  2771 s as answer
    Case 2: 10 term approximation gives time =  2966 s as answer


    