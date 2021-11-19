# rebound_tools_js
Tools from the REBOUND N-body code (written in C) translated to JavaScript

For Hanno Rein's planetary dynamics course (November 2021), here are REBOUND tools translated from C to JavaScript:

reb_tools_orbit_to_particle_err()
This function takes in a orbital parameters (a, e, inc, Omega, omega, f) and translates them into cartesian coordinates (x, y, z, vx, vy, vz).

reb_tools_particle_to_orbit_err()
This function is the reverse coordinate transformation of the previous function.

Source code: 
https://github.com/hannorein/rebound/blob/main/src/tools.c#L838
https://github.com/hannorein/rebound/blob/main/src/tools.c#L950
