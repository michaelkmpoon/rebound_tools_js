const TINY = 1e-308         // Close to smallest representable floating point number, used for orbit calculation
const MIN_INC = 1e-8		///< Below this inclination, the broken angles pomega and theta equal the corresponding 
							///< unbroken angles to within machine precision, so a practical boundary for planar orbits
							//
const MIN_ECC = 1e-8        ///< Below this eccentricity, corrections at order e^2 are below machine precision, so we use
                            ///< stable expressions accurate to O(e) for the mean longitude below for near-circular orbits.

function reb_tools_orbit_to_particle_err(G, primary, m, a, e, inc, Omega, omega, f, errObj) {
    if (e == 1) {
        // JS doesn't have pointers so we pass in errObj as an object.
        errObj.code = 1 // Can't initialize a radial orbit with orbital elements
        return {}       // return an empty object
    }
    if (e < 0.) {
        errObj.code = 2 // Eccetricity must be greater or equal to zero.
        return {}
    }
    if (e > 1.) {
        if (a > 0.) {
            errObj.code = 3 // Bound orbit (a > 0) must have e < 1.
            return {}
        }
    }
    else {
        if (a < 0.){
            errObj.code = 4 // Unbound orbit (a < 0) must have e > 1.
            return {}
        } 
    }
    if (e*Math.cos(f) < -1.) {
        /* Unbound orbit can't have f set beyond range
           allowed by the asymptotes set by the parabola. */
        errObj.code = 5 
        return {}
        }
    if (primary.m < TINY) {
        errObj.code = 6 // Primary has no mass.
        return {}   
    }
    const p = {
        m: 0,
        x: 0,
        y: 0,
        z: 0,
        vx: 0,
        vy: 0,
        vz: 0,
        ax: 0,
        ay: 0,
        az: 0,
    }

    p.m = m
    const r = a*(1-e*e)/(1 + e*Math.cos(f));
    const v0 = Math.sqrt(G*(m+primary.m)/a/(1.-e*e)); // in this form it works for elliptical and hyperbolic orbits

    const cO = Math.cos(Omega)
    const sO = Math.sin(Omega)
    const co = Math.cos(omega)
    const so = Math.sin(omega)
    const cf = Math.cos(f)
    const sf = Math.sin(f)
    const ci = Math.cos(inc)
    const si = Math.sin(inc)

    // Murray & Dermott Eq 2.122
    p.x = primary.x + r*(cO*(co*cf-so*sf) - sO*(so*cf+co*sf)*ci)
    p.y = primary.y + r*(sO*(co*cf-so*sf) + cO*(so*cf+co*sf)*ci)
    p.z = primary.z + r*(so*cf+co*sf)*si

    // Murray & Dermott Eq. 2.36 after applying the 3 rotation matrices from Sec. 2.8 to the velocities in the orbital plane
    p.vx = primary.vx + v0*((e+cf)*(-ci*co*sO - cO*so) - sf*(co*cO - ci*so*sO))
    p.vy = primary.vy + v0*((e+cf)*(ci*co*cO - sO*so)  - sf*(co*sO + ci*so*cO))
    p.vz = primary.vz + v0*((e+cf)*co*si - sf*si*so)

    p.ax = 0	
    p.ay = 0 	
    p.az = 0

    return p
}

// returns acos(num/denom), using disambiguator to tell which quadrant to return.  
// will return 0 or pi appropriately if num is larger than denom by machine precision
// and will return 0 if denom is exactly 0.
function acos2(num, denom, disambiguator){
    const cosine = num/denom
    if (cosine > -1. && cosine < 1.) {
        val = Math.acos(cosine)
        if (disambiguator < 0.) {
            val = -val
        }
    }
    else {
        if (cosine <= -1) {
            val = Math.PI
        }
        else {
            val = 0.
        }
    }
    return val
}

function reb_tools_mod2pi(f){
    const pi2 = 2.*Math.PI;
    return (pi2 + (f % pi2)) % pi2
}
  
function reb_tools_particle_to_orbit_err(G, p, primary, errObj) {
    const o = {
        a: 0,
        e: 0,
        inc: 0,
        Omega: 0,
        omega: 0,
        pomega: 0,
        f: 0,
        M: 0,
        E: 0,
        theta: 0,
        T: 0,
        d: 0,
        h: 0,
        l: 0,
        v: 0,
        n: 0,
        P: 0,
        ax: 0,
        ay: 0,
        az: 0,
    }
    if (primary.m < TINY) {
        errObj.code = 1 // Primary has no mass.
        return {}   
    }
    const mu = G*(p.m+primary.m)
    const dx = p.x - primary.x
    const dy = p.y - primary.y
    const dz = p.z - primary.z
    const dvx = p.vx - primary.vx
    const dvy = p.vy - primary.vy
    const dvz = p.vz - primary.vz
    o.d = Math.sqrt ( dx*dx + dy*dy + dz*dz )

    const vsquared = dvx*dvx + dvy*dvy + dvz*dvz
    o.v = Math.sqrt(vsquared)
    const vcircsquared = mu/o.d
    o.a = -mu/( vsquared - 2.*vcircsquared )	// semi major axis

    const hx = (dy*dvz - dz*dvy) 					//angular momentum vector
    const hy = (dz*dvx - dx*dvz)
    const hz = (dx*dvy - dy*dvx)
    o.h = Math.sqrt ( hx*hx + hy*hy + hz*hz )		// abs value of angular momentum

    const vdiffsquared = vsquared - vcircsquared	
    if(o.d <= TINY){							
        errObj = 2									// particle is on top of primary
        return {}
    }
    const vr = (dx*dvx + dy*dvy + dz*dvz)/o.d
    const rvr = o.d*vr
    const muinv = 1./mu

    const ex = muinv*( vdiffsquared*dx - rvr*dvx )
    const ey = muinv*( vdiffsquared*dy - rvr*dvy )
    const ez = muinv*( vdiffsquared*dz - rvr*dvz )
    o.e = Math.sqrt( ex*ex + ey*ey + ez*ez )		// eccentricity
    o.n = o.a/Math.abs(o.a)*Math.sqrt(Math.abs(mu/(o.a*o.a*o.a)))	// mean motion (negative if hyperbolic)
    o.P = 2*Math.PI/o.n									// period (negative if hyperbolic)


    o.inc = acos2(hz, o.h, 1.)			// cosi = dot product of h and z unit vectors.  Always in [0,pi], so pass dummy disambiguator
    // will = 0 if h is 0.

    const nx = -hy							// vector pointing along the ascending node = zhat cross h
    const ny =  hx		
    const n = Math.sqrt( nx*nx + ny*ny )

    // Omega, pomega and theta are measured from x axis, so we can always use y component to disambiguate if in the range [0,pi] or [pi,2pi]
    o.Omega = acos2(nx, n, ny)			// cos Omega is dot product of x and n unit vectors. Will = 0 if i=0.

    if(o.e < 1.){
        ea = acos2(1.-o.d/o.a, o.e, vr)// from definition of eccentric anomaly.  If vr < 0, must be going from apo to peri, so ea = [pi, 2pi] so ea = -acos2(cosea)
        o.M = ea - o.e*Math.sin(ea)			// mean anomaly (Kepler's equation)
    }
    else{
        ea = Math.acosh((1.-o.d/o.a)/o.e)
        if(vr < 0.){                    // Approaching pericenter, so eccentric anomaly < 0.
            ea = -ea
        }
        o.M = o.e*Math.sinh(ea) - ea          // Hyperbolic Kepler's equation
    }

    // in the near-planar case, the true longitude is always well defined for the position, and pomega for the pericenter if e!= 0
    // we therefore calculate those and calculate the remaining angles from them
    if(o.inc < MIN_INC || o.inc > Math.PI - MIN_INC){	// nearly planar.  Use longitudes rather than angles referenced to node for numerical stability.
        o.theta = acos2(dx, o.d, dy)		// cos theta is dot product of x and r vectors (true longitude). 
        o.pomega = acos2(ex, o.e, ey)		// cos pomega is dot product of x and e unit vectors.  Will = 0 if e=0.

        if(o.inc < Math.PI/2.){
            o.omega = o.pomega - o.Omega
            o.f = o.theta - o.pomega
            if(o.e > MIN_ECC){              // pomega well defined
                o.l = o.pomega + o.M
            }
            else{                           // when e << 1 and pomega ill defined, use l = theta+(M-f). M-f is O(e) so well behaved
                o.l = o.theta - 2.*o.e*Math.sin(o.f) // M-f from Murray & Dermott Eq 2.93. This way l->theta smoothly as e->0
            }
        }
        else{
            o.omega = o.Omega - o.pomega
            o.f = o.pomega - o.theta
            if(o.e > MIN_ECC){              // pomega well defined
                o.l = o.pomega - o.M
            }
            else{                           // when e << 1 and pomega ill defined, use l = theta+(M-f). M-f is O(e) so well behaved
                o.l = o.theta + 2.*o.e*Math.sin(o.f) // M-f from Murray & Dermott Eq 2.93 (retrograde changes sign). This way l->theta smoothly as e->0
            }
        }

    }
    // in the non-planar case, we can't calculate the broken angles from vectors like above.  omega+f is always well defined, and omega if e!=0
    else{
        const wpf = acos2(nx*dx + ny*dy, n*o.d, dz)	// omega plus f.  Both angles measured in orbital plane, and always well defined for i!=0.
        o.omega = acos2(nx*ex + ny*ey, n*o.e, ez)
        if(o.inc < Math.PI/2.){
            o.pomega = o.Omega + o.omega
            o.f = wpf - o.omega
            o.theta = o.Omega + wpf
            if(o.e > MIN_ECC){              // pomega well defined
                o.l = o.pomega + o.M
            }
            else{                           // when e << 1 and pomega ill defined, use l = theta+(M-f). M-f is O(e) so well behaved
                o.l = o.theta - 2.*o.e*Math.sin(o.f) // M-f from Murray & Dermott Eq 2.93. This way l->theta smoothly as e->0
            }
        }
        else{
            o.pomega = o.Omega - o.omega
            o.f = wpf - o.omega
            o.theta = o.Omega - wpf
            if(o.e > MIN_ECC){              // pomega well defined
                o.l = o.pomega - o.M
            }
            else{                           // when e << 1 and pomega ill defined, use l = theta+(M-f). M-f is O(e) so well behaved
                o.l = o.theta + 2.*o.e*Math.sin(o.f) // M-f from Murray & Dermott Eq 2.93 (retrograde changes sign). This way l->theta smoothly as e->0
            }
        }
    }

    // move some of the angles into [0,2pi) range
    o.f = reb_tools_mod2pi(o.f)
    o.l = reb_tools_mod2pi(o.l)
    o.M = reb_tools_mod2pi(o.M)
    o.theta = reb_tools_mod2pi(o.theta)
    o.omega = reb_tools_mod2pi(o.omega)
    return o;
}

// Test cases

// Test primary
const primary = {
    m: 1, x: 1, y: 1, z: 1, vx: 0, vy: 0, vz: 0,
}

// Test orbit
const test_o = {
    m: 0.01, a: 1, e: 0.1, inc: 1., Omega: 0.1, omega: 0.1, f: 0.1,
}

// Test case 1: simple orbit
console.log('orbit properties: a=', test_o.a, ' e=', test_o.e, ' inc=', test_o.inc, 'Omega=', test_o.Omega, ' omega=', test_o.omega, ' f=', test_o.f)
p = reb_tools_orbit_to_particle_err(G=1., primary, m=test_o.m, a=test_o.a, e=test_o.e, inc=test_o.inc, Omega=test_o.Omega, omega=test_o.omega, f=test_o.f, errObj={code : 0})
o = reb_tools_particle_to_orbit_err(G=1., p, primary, errObjerrObj={code : 0})
console.log('orbit properties: a=', o.a, ' e=', o.e, ' inc=', o.inc, 'Omega=', o.Omega, ' omega=', o.omega, ' f=', o.f)




// Test case 2: e > 1
test_o.e = 10.
p = reb_tools_orbit_to_particle_err(G=1., primary, m=test_o.m, a=test_o.a, e=test_o.e, inc=test_o.inc, Omega=test_o.Omega, omega=test_o.omega, f=test_o.f, errObj={code : 0})
console.log(errObj.code)