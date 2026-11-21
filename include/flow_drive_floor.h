#pragma once

// The zero-drive convention of the network's mass-flow rows (#1193, #1209,
// #1256).
//
// A flow law driven by a pressure drop -- an isentropic jet (orifice, hole,
// film row), a Fanno duct -- rises like sqrt(drive) from zero: its slope is
// infinite at zero drive, and evaluated in floating point just above zero
// the law is a round-off staircase (the nozzle's 1 - p/p0 carries a few
// digits at a drive of 1e-14 of the pressure). Newton then reads a row that
// is blind to its own drive, or sees a slope that flips between 0 and -22 on
// a 1e-10 Pa change of a start pressure.
//
// Below a drive floor such a law is replaced by an odd cubic in
// t = drive / floor through zero that matches the law's value G0 AND slope
// at the floor:
//   h(t) = (3 - sigma)/2 t + (sigma - 1)/2 t^3,   h(1) = 1, h'(1) = sigma,
// with sigma = (dG/d drive) floor / G0 the law's local exponent there (1/2
// for a jet; larger for a laminar duct, where friction rises as the flow
// falls). The flux is G0 h(t): C1 at the floor, finite and positive in slope
// at zero (h'(0) = (3 - sigma)/2), monotone for sigma < 3, and the exact law
// from the floor up. At exactly zero drive the flux is 0 and its slope is the
// cubic's h'(0) G0 / floor.
//
// The Fanno duct's floor is kFannoFlowDriveFloor (fanno_mach.h, 1 Pa since
// #1017); the jet kernels' is kJetDriveFloor below.

namespace combaero {

// Drive [Pa] below which the compressible jet kernels (the nozzle orifice and
// Jobson's) take the cubic. 1e-3 Pa: the nozzle's staircase ends near 2e-9
// Pa at 1 bar (~2e-14 of P0), so the floor is 1e5x above it there and 1e3x
// at 100 bar, and the floor's flow and slope carry ~1e-8 relative round-off
// at 1 bar. It is kept that small because inside it the flow is not the
// jet's (h(t) < sqrt(t)), and a cold solve whose iterate crosses a reversing
// jet's zero drop sees that change: at a 1 Pa floor the CR-3630 D.2
// B(5,4,1) mc/mj 1.01 network took another path there and stalled (#1256).
// The slope at zero is 1.25 G(floor) / floor, ~32x a 1 Pa floor's: finite.
constexpr double kJetDriveFloor = 1.0e-3;

// Largest local exponent the cubic accepts: h is monotone for sigma < 3.
constexpr double kDriveFloorSigmaMax = 2.9;

struct DriveFloorCubic {
  double h = 0.0;   // G / G0
  double dh = 0.0;  // dh/dt
};

// h(t) and h'(t) of the sub-floor cubic, t = drive / floor in [0, 1], sigma
// the law's local exponent at the floor. The callers hold sigma fixed in
// their partials; each states what that costs.
inline DriveFloorCubic drive_floor_cubic(double t, double sigma) {
  DriveFloorCubic c;
  c.h = 0.5 * (3.0 - sigma) * t + 0.5 * (sigma - 1.0) * t * t * t;
  c.dh = 0.5 * (3.0 - sigma) + 1.5 * (sigma - 1.0) * t * t;
  return c;
}

}  // namespace combaero
