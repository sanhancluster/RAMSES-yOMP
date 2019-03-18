Ricarda Beckmann, April 2018
====================

The files in this patch provide the initial conditions for an isothermal
sphere of a given radius, with a black hole sink particle at the centre.
The sphere is set up to collapse onto the black hole and be entirely
accreted by it over the course of ~ 20 kyr for the default initial
condtions.

The initial conditions in this file are based on the isothermal
sphere in Federrath2010, Section 3.4.

-------
INPUT PARAMETERS:

added via gravity_params, in the order
* T_sphere: temperature of the sphere in K
* R0: Truncation radius of the sphere in units of 5E16 cm
* rho_R0: Density at the truncation radius, in units of 3.83E-18 g/cm**3
* rho_out: Density outside the sphere in units of rho_R0

To reproduce Federrath2010 exactly, set:  gravity_params= 10,1,1,1E-2