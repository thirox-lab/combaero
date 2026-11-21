#pragma once

// The zero-drive convention of the network's mass-flow rows (#1193, #1209,
// #1256).
//
// A flow law driven by a pressure drop -- an isentropic jet (orifice, hole,
// film row), a Fanno duct -- rises like sqrt(drive) from zero: its slope is
// infinite at zero drive, and evaluated in floating point just above zero
// the nozzle law is a round-off staircase (its 1 - p/p0 loses the drive's
// digits). Newton then reads a row that is blind to its own drive, or sees a
// slope that flips between 0 and -22 on a 1e-10 Pa change of a start
// pressure.
//
// Below a drive floor such a law is replaced by an odd cubic in
// t = drive / floor through zero that matches the law's value G0 at the
// floor and, through sigma, its slope:
//   h(t) = (3 - sigma)/2 t + (sigma - 1)/2 t^3,   h(1) = 1, h'(1) = sigma,
// sigma the law's local exponent (1/2 for a jet; larger for a laminar duct,
// where friction rises as the flow falls). The flux is G0 h(t): finite and
// positive in slope at zero (h'(0) = (3 - sigma)/2), monotone for sigma < 3,
// and the exact law from the floor up. At exactly zero drive the flux is 0
// and its slope is the cubic's h'(0) G0 / floor.
//
// The Fanno duct's floor is kFannoFlowDriveFloor (fanno_mach.h, 1 Pa since
// #1017), with sigma matched to the duct's slope there; the jet kernels' is
// kJetDriveFloor below, with sigma the constant kJetFloorExponent.

namespace combaero {

// Drive [Pa] below which the compressible jet kernels (the nozzle orifice and
// Jobson's) take the cubic. What it guarantees: below it there is no
// staircase at any pressure and the slope is finite (1.25 G(floor) / floor at
// zero). Above it the jet law is main's, and so is its round-off: the
// nozzle's staircase (flat treads in G against the drive, measured in the
// #1259 review) ends near 8.7e-7 Pa at 1 bar, 3 decades below the floor;
// near 9.7e-4 Pa at 100 bar, at the floor; and at 1000 bar it continues past
// the floor to ~1e-2 Pa, with ~3e-4 relative noise in G there (pre-existing,
// #1261). Jobson's kernel is well conditioned: below 1e-6 Pa even at 1000 bar.
// Inside the floor the flow is not the jet's (h(t) < sqrt(t)), so the floor is
// kept small. It is not tuned to a test: the CR-3630 D.2 B(5,4,1) mc/mj 1.01
// cold solve, whose iterate crosses a reversing jet's zero drop, passes at
// floors from 1e-4 to 1e-1 Pa and stalls only at 1 Pa, and the Brayton start
// converges at every floor from 1e-4 to 1 Pa.
constexpr double kJetDriveFloor = 1.0e-3;

// The jet cubic's exponent: a constant, so the band's partials are exact (a
// sigma matched to the law's slope moves with the state, and its own partials
// would need the law's second derivatives). The flow is then continuous at the
// floor but its slope steps there: the nozzle's by 4e-7 at 1 bar and 6.1e-5 at
// 100 bar; Jobson's, with an approach, by up to 27% (see jobson_eval).
constexpr double kJetFloorExponent = 0.5;

// Largest local exponent the cubic accepts: h is monotone for sigma < 3.
constexpr double kDriveFloorSigmaMax = 2.9;

struct DriveFloorCubic {
  double h = 0.0;   // G / G0
  double dh = 0.0;  // dh/dt
};

// h(t) and h'(t) of the sub-floor cubic, t = drive / floor in [0, 1].
inline DriveFloorCubic drive_floor_cubic(double t, double sigma) {
  DriveFloorCubic c;
  c.h = 0.5 * (3.0 - sigma) * t + 0.5 * (sigma - 1.0) * t * t * t;
  c.dh = 0.5 * (3.0 - sigma) + 1.5 * (sigma - 1.0) * t * t;
  return c;
}

}  // namespace combaero
