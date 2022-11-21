/*
 * Copyright 2013 Daniel Warner <contact@danrw.com>
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 * http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */


#include <Observer.h>
#include <SGP4.h>
#include <Util.h>
#include <CoordTopocentric.h>
#include <CoordGeodetic.h>
#include <SolarPosition.h>

#include <cmath>
#include <iostream>
#include <list>

/*
 * http://dkami.umcs.pl/how-to-calculate-flares-from-satellites/
 */

static void
Magnitude(libsgp4::Vector& in_v)
{
	using std::pow;
	using std::sqrt;
	in_v.w = sqrt(
		pow(in_v.x, 2) +
		pow(in_v.y, 2) +
		pow(in_v.z, 2)
	);
}

static double
Dot(const libsgp4::Vector& in_v1, const libsgp4::Vector& in_v2)
{
	return (in_v1.x * in_v2.x) +
		(in_v1.y * in_v2.y) +
		(in_v1.z * in_v2.z);
}

static void
Cross(const libsgp4::Vector& in_v1, 
	const libsgp4::Vector& in_v2,
	libsgp4::Vector& out_v)
{
	out_v.x = in_v1.y * in_v2.z - in_v1.z * in_v2.y;
	out_v.y = in_v1.x * in_v2.z - in_v1.z * in_v2.x;
	out_v.z = in_v1.x * in_v2.y - in_v1.y * in_v2.x;
	Magnitude(out_v);
}

static void
Scalar_Multiply(double in_d, 
		const libsgp4::Vector& in_v, 
		libsgp4::Vector& out_v)
{
	using std::abs;
	out_v.x = in_d * in_v.x;
	out_v.y = in_d * in_v.y;
	out_v.z = in_d * in_v.z;
	out_v.w = abs(in_d) * in_v.w;
}

static void
Vec_Sub(const libsgp4::Vector& in_v1, 
	const libsgp4::Vector& in_v2, 
	libsgp4::Vector& out_v3)
{
	out_v3.x = in_v1.x - in_v2.x;
	out_v3.y = in_v1.y - in_v2.y;
	out_v3.z = in_v1.z - in_v2.z;
	Magnitude(out_v3);
}

static void 
Normalize(libsgp4::Vector& in_v)
{
	/* Normalize a vector */
	if(in_v.w == 0.0) return; // no divide by zero
	in_v.x /= in_v.w;
	in_v.y /= in_v.w;
	in_v.z /= in_v.w;
}

static double 
Angle(const libsgp4::Vector& in_v1, const libsgp4::Vector& in_v2)
{
	using std::acos;
	/* Calculates the angle between vectors v1 and v2 */
	libsgp4::Vector v1 = in_v1;
	libsgp4::Vector v2 = in_v2;
	Magnitude(v1);
	Magnitude(v2);
	return (acos(Dot(v1, v2) / (v1.w * v2.w)));
}

static double 
reflection_angle(
	const libsgp4::Vector& in_sun_pos, // In ECI reference frame
	const libsgp4::Vector& in_sat_pos, // In ECI reference frame
	const libsgp4::Vector& in_sat_vel, // In ECI reference frame
	const libsgp4::Vector& in_obs,     // In ECI reference frame
	double in_alpha, // reflection surface angle in xz plane
	double in_beta,  // refelction surface angle in xy plane 
	double *sun_sat_obs)
{
	double angle = 0.0;
	double PP[3][3];
	libsgp4::Vector xx, yy, zz;
	libsgp4::Vector n1, n2, n3, reflected, reflected_back;

	in_alpha = libsgp4::Util::DegreesToRadians(in_alpha);
	in_beta  = libsgp4::Util::DegreesToRadians(in_beta);

	Scalar_Multiply(1.0, in_sat_vel, xx);
	Normalize(xx);
	Scalar_Multiply(1.0, in_sat_pos, zz);
	Normalize(zz);

	Cross(xx, zz, yy);

	double temp1, temp2, nx, ny, nz;
	double sin_alpha = std::sin(in_alpha);
	double cos_alpha = std::cos(in_alpha);
	double sin_beta  = std::sin(in_beta);
	double cos_beta  = std::cos(in_beta);

	// x of mirror
	nx = 1.0; ny = 0.0; nz = 0.0;
	n1.x = nx * cos_alpha - nz * sin_alpha;
	n1.y = ny;
	n1.z = nx * sin_alpha + nz * cos_alpha;
	temp1 = n1.x; temp2 = n1.y;
	n1.x = temp1 * cos_beta + temp2 * sin_beta;
	n1.y = temp1 * -sin_beta + temp2 * cos_beta;
	PP[0][0] = n1.x * xx.x + n1.y * yy.x + n1.z * zz.x;
	PP[1][0] = n1.x * xx.y + n1.y * yy.y + n1.z * zz.y;
	PP[2][0] = n1.x * xx.z + n1.y * yy.z + n1.z * zz.z;

	// y of mirror
	nx = 0.0; ny = 1.0; nz = 0.0;
	n2.x = nx * cos_alpha - nz * sin_alpha;
	n2.y = ny;
	n2.z = nx * sin_alpha + nz * cos_alpha;
	temp1 = n2.x; temp2 = n2.y;
	n2.x = temp1 * cos_beta + temp2 * sin_beta;
	n2.y = temp1 * -sin_beta + temp2 * cos_beta;
	PP[0][1] = n2.x * xx.x + n2.y * yy.x + n2.z * zz.x;
	PP[1][1] = n2.x * xx.y + n2.y * yy.y + n2.z * zz.y;
	PP[2][1] = n2.x * xx.z + n2.y * yy.z + n2.z * zz.z;

	// z of mirror
	nx = 0.0; ny = 0.0; nz = 1.0;
	n3.x = nx * cos_alpha - nz * sin_alpha;
	n3.y = ny;
	n3.z = nx * sin_alpha + nz * cos_alpha;
	temp1 = n3.x; temp2 = n3.y;
	n3.x = temp1 * cos_beta + temp2 * sin_beta;
	n3.y = temp1 * -sin_beta + temp2 * cos_beta;
	PP[0][2] = n3.x * xx.x + n3.y * yy.x + n3.z * zz.x;
	PP[1][2] = n3.x * xx.y + n3.y * yy.y + n3.z * zz.y;
	PP[2][2] = n3.x * xx.z + n3.y * yy.z + n3.z * zz.z;

	// Calculate difference between sat vpos vec and obs vec
	libsgp4::Vector sat_obs;
	Vec_Sub(in_sat_pos, in_obs, sat_obs);
	Normalize(sat_obs);

	// Convert the obs-sat vector to the mirror coord system
	reflected.x = 
		sat_obs.x * PP[0][0] + 
		sat_obs.y * PP[1][0] + 
		sat_obs.z * PP[2][0];
	reflected.y = 
		sat_obs.x * PP[0][1] + 
		sat_obs.y * PP[1][1] + 
		sat_obs.z * PP[2][1];
	reflected.z = 
		sat_obs.x * PP[0][2] + 
		sat_obs.y * PP[1][2] + 
		sat_obs.z * PP[2][2];

	reflected.x *= -1.0; // mirror observer location

	if(reflected.x >= 0) {
		// reflection in x must be positive direction
		reflected_back.x = 
			reflected.x * PP[0][0] + 
			reflected.y * PP[0][1] + 
			reflected.z * PP[0][2];
		reflected_back.y = 
			reflected.x * PP[1][0] + 
			reflected.y * PP[1][1] + 
			reflected.z * PP[1][2];
		reflected_back.z = 
			reflected.x * PP[2][0] + 
			reflected.y * PP[2][1] + 
			reflected.z * PP[2][2];
		reflected_back.w = 0.0;
		libsgp4::Vector sun_sat;
		Vec_Sub(in_sun_pos, in_sat_pos, sun_sat);
		Normalize(sun_sat);
		angle = Angle(sun_sat, reflected_back);
		if(sun_sat_obs) {
			*sun_sat_obs = M_PI - Angle(sat_obs, sun_sat);
		}
		return angle;
	}

	if(sun_sat_obs) {
		*sun_sat_obs = 0.0;
	}
	return 100.0; // nothing calculated;
}

/*
 * Satellite visibility is based on:-
 * see : Visually Observing Earth Satellites By Dr. T.S. Kelso
 * at  :  https://celestrak.org/columns/v03n01/
 */
static
bool is_illuminated(
	const libsgp4::Eci& sat, 
	const libsgp4::Eci& sun)
{
    using std::abs;
    using std::pow;
    using std::asin;
    using std::acos;
    using std::sqrt;

    bool rval = true;
    libsgp4::Vector satvec = sat.Position();
    libsgp4::Vector sunvec = sun.Position();
    
    double psun = sqrt(
		      pow(satvec.x - sunvec.x, 2) + 
		      pow(satvec.y - sunvec.y, 2) + 
		      pow(satvec.z - sunvec.z, 2) 
		  );
    double pearth = sqrt(
		        pow(satvec.x, 2) + 
   		        pow(satvec.y, 2) + 
		        pow(satvec.z, 2)
		    );
    double theta_e = asin((12742. / 2.) / pearth);
    double theta_s = asin((1391000. / 2.) / psun);
    double dotproduct_peps = 
	    	abs(
			(satvec.x * sunvec.x) + 
			(satvec.y * sunvec.y) + 
			(satvec.z * sunvec.z)
		);
    double theta = acos(dotproduct_peps / (psun * pearth));
    if(theta_e > theta_s && theta < (theta_e - theta_s)) {
        rval = false; // Satellite is not lit
    }
    else if(abs(theta_e - theta_s) < theta && theta < (theta_e + theta_s)) {
        rval = true; // Satellite is penumbral  
    }
    else if (theta_s > theta_e && theta < (theta_s - theta_e)) {
        rval = true; // Satellite is eclipsed
    }

    return rval;
}

struct SatPoint
{
	bool is_lit;
	double azm;
	double alt;
	double sunalt;
	double refl_angle;
	double sun_sat_obs;
	libsgp4::DateTime at_time;
};

static std::list<struct SatPoint>
build_list(const libsgp4::CoordGeodetic& geo,
		libsgp4::SGP4& sgp4,
		libsgp4::DateTime& start_date,
		bool in_only_lit = true)
{
	std::list<struct SatPoint> slist;
	libsgp4::SolarPosition sp;
	libsgp4::Observer obs(geo);
	libsgp4::DateTime current_time = start_date;
	libsgp4::DateTime end_time = current_time + libsgp4::TimeSpan(8, 0, 0);
	for( ; current_time < end_time ; current_time = current_time + libsgp4::TimeSpan(0, 0, 10)) {
		libsgp4::Eci obs_eci(current_time, geo);
		libsgp4::Eci sat_eci = sgp4.FindPosition(current_time);
		libsgp4::CoordTopocentric topo = obs.GetLookAngle(sat_eci);
		if(libsgp4::Util::RadiansToDegrees(topo.elevation) <= 1.0) continue; 
		libsgp4::Eci sun_eci = sp.FindPosition(current_time);
		bool is_lit = is_illuminated(sat_eci, sun_eci);
		if(in_only_lit && !is_lit) continue;
		libsgp4::CoordTopocentric suntopo = obs.GetLookAngle(sun_eci);
		SatPoint r;
		r.azm = topo.azimuth;
		r.alt = topo.elevation;
		r.sunalt = suntopo.elevation;
		r.at_time = current_time;
		r.is_lit = is_lit;
		r.refl_angle = reflection_angle(
			sun_eci.Position(),
			sat_eci.Position(),
			sat_eci.Velocity(),
			obs_eci.Position(), 
			-90.0, 0.0,
			&r.sun_sat_obs
		);

		slist.push_back(r);
	}
	return slist;
}

int main()
{
    std::stringstream oss;
    libsgp4::CoordGeodetic geo(56.21, -3.0026, 0.05);
    libsgp4::Tle tle("BLUEWALKER3",
        "1 53807U 22111AL  22324.21613640  .00001955  00000-0  10857-3 0  9990",
        "2 53807  53.2014 307.5499 0014123 131.7034 228.5173 15.18600551 10671"
	);
    libsgp4::SGP4 sgp4(tle);

    //std::cout << tle << std::endl;

    //libsgp4::DateTime start_date = libsgp4::DateTime::Now(true);
    //start_date = start_date - libsgp4::TimeSpan(4, 0, 0);

    libsgp4::DateTime start_date = libsgp4::DateTime(2022, 11, 29, 17, 0, 0);

    std::cout << "Start time: " << start_date << std::endl;

    std::list<struct SatPoint> slist = build_list(geo, sgp4, start_date);

    if(slist.begin() == slist.end()) {
	oss << "No passes found" << std::endl;	
    }
    else {
	std::list<struct SatPoint>::const_iterator itor = slist.begin();
	for( ; itor != slist.end(); itor++) {
		oss << "Time: " << itor->at_time << " "
			<< "Azm: " << std::setw(6) << std::setprecision(1) << std::fixed << libsgp4::Util::RadiansToDegrees(itor->azm) << " "
			<< "Alt: " << std::setw(6) << std::setprecision(1) << std::fixed << libsgp4::Util::RadiansToDegrees(itor->alt) << " "
			<< "Sun: " << std::setw(6) << libsgp4::Util::RadiansToDegrees(itor->sunalt) << " "
			<< "R:" << std::setw(6) << libsgp4::Util::RadiansToDegrees(itor->refl_angle) << " "
			<< "S:" << std::setw(6) << libsgp4::Util::RadiansToDegrees(itor->sun_sat_obs) << " "
			<< "Lit: " << (itor->is_lit == true ? "1":"0")
			<< std::endl;
	}
	std::cout << oss.str();
	std::cout << "Found in total: " << slist.size() << std::endl;
    }
    
    return 0;
}
