/**
 * Standalone Node.js physics verification suite for script.js.
 * Validates Hermite-Gaussian intensity equations, nodal zeros, peak locations,
 * power conservation across all 4 spatial modes, and beam waist scaling.
 */

"use strict";

const assert = require("assert");
const {
  MODE_DEFINITIONS,
  hermitePolynomial,
  calculateHermiteGaussianIntensity,
  getTheoreticalPeakRatio,
  getSliceCoordinates,
  computeBeamParameters
} = require("./script.js");

function approxEqual(actual, expected, tol = 1e-6, label = "") {
  const diff = Math.abs(actual - expected);
  assert(
    diff <= tol,
    `${label}: expected ${expected}, got ${actual} (diff=${diff} > tol=${tol})`
  );
}

console.log("Running Laser Beam Profiling Physics Verification...\n");

// 1. Verify Hermite polynomials
approxEqual(hermitePolynomial(0, 2.5), 1, 1e-12, "H_0(2.5)");
approxEqual(hermitePolynomial(1, 2.5), 5.0, 1e-12, "H_1(2.5)");
approxEqual(hermitePolynomial(2, 1.5), 4 * 1.5 * 1.5 - 2, 1e-12, "H_2(1.5)");

// 2. Verify TEM00 intensity at origin and 1/e^2 waist radius
const w0 = 1.2;
const I0 = 150;
approxEqual(
  calculateHermiteGaussianIntensity(0, 0, w0, I0, 0, 0),
  I0,
  1e-9,
  "TEM00 center peak"
);
approxEqual(
  calculateHermiteGaussianIntensity(w0, 0, w0, I0, 0, 0),
  I0 * Math.exp(-2),
  1e-9,
  "TEM00 1/e^2 intensity at r = w0"
);

// 3. Verify nodal zeros and exact peak locations for TEM01, TEM10, TEM11
const lobeOffset = w0 / Math.SQRT2;

// TEM01: nodal line y = 0, peak at (0, ±w0/sqrt(2))
approxEqual(calculateHermiteGaussianIntensity(0.5, 0, w0, I0, 0, 1), 0, 1e-12, "TEM01 nodal line y=0");
approxEqual(
  calculateHermiteGaussianIntensity(0, lobeOffset, w0, I0, 0, 1),
  I0 * 2 * Math.exp(-1),
  1e-9,
  "TEM01 peak at (0, +w0/sqrt(2))"
);

// TEM10: nodal line x = 0, peak at (±w0/sqrt(2), 0)
approxEqual(calculateHermiteGaussianIntensity(0, -0.7, w0, I0, 1, 0), 0, 1e-12, "TEM10 nodal line x=0");
approxEqual(
  calculateHermiteGaussianIntensity(-lobeOffset, 0, w0, I0, 1, 0),
  I0 * 2 * Math.exp(-1),
  1e-9,
  "TEM10 peak at (-w0/sqrt(2), 0)"
);

// TEM11: nodal lines x = 0 and y = 0, peak at (±w0/sqrt(2), ±w0/sqrt(2))
approxEqual(calculateHermiteGaussianIntensity(0, lobeOffset, w0, I0, 1, 1), 0, 1e-12, "TEM11 nodal x=0");
approxEqual(calculateHermiteGaussianIntensity(lobeOffset, 0, w0, I0, 1, 1), 0, 1e-12, "TEM11 nodal y=0");
approxEqual(
  calculateHermiteGaussianIntensity(lobeOffset, -lobeOffset, w0, I0, 1, 1),
  I0 * 4 * Math.exp(-2),
  1e-9,
  "TEM11 peak at (+w0/sqrt(2), -w0/sqrt(2))"
);

// 4. Verify 2D numerical integration (power conservation across all 4 modes)
// Integrate over [-5*w0, +5*w0]^2 using midpoint rule
const testW0 = 1.0; // mm = 0.1 cm
const testI0 = 100.0; // W/cm^2
const expectedPowerW = 0.5 * Math.PI * (testW0 * 0.1) ** 2 * testI0;
const limit = 5.0 * testW0;
const steps = 200;
const dMm = (2 * limit) / steps;
const dCm2 = (dMm * 0.1) ** 2;

for (const key of Object.keys(MODE_DEFINITIONS)) {
  const { m, n } = MODE_DEFINITIONS[key];
  let integralW = 0;
  for (let i = 0; i < steps; i++) {
    const x = -limit + (i + 0.5) * dMm;
    for (let j = 0; j < steps; j++) {
      const y = -limit + (j + 0.5) * dMm;
      integralW += calculateHermiteGaussianIntensity(x, y, testW0, testI0, m, n) * dCm2;
    }
  }
  approxEqual(integralW, expectedPowerW, 1e-4, `${key} total integrated power (W)`);
  console.log(`  [PASS] ${key}: Integrated Power = ${integralW.toFixed(6)} W (Expected ${expectedPowerW.toFixed(6)} W), Peak Ratio = ${getTheoreticalPeakRatio(m, n).toFixed(4)}`);
}

// 5. Verify derived parameters and beam waist scaling
const paramsSmall = computeBeamParameters({ mode: "TEM00", w0: 0.5, I0: 100, wavelengthNm: 632.8 });
const paramsLarge = computeBeamParameters({ mode: "TEM00", w0: 1.0, I0: 100, wavelengthNm: 632.8 });
approxEqual(paramsLarge.rayleighRangeM / paramsSmall.rayleighRangeM, 4.0, 1e-9, "Rayleigh range scales with w0^2");
approxEqual(paramsSmall.divergenceMrad / paramsLarge.divergenceMrad, 2.0, 1e-9, "Divergence scales inversely with w0");
approxEqual(paramsLarge.totalPowerW, 0.5 * Math.PI * 0.01 * 100, 1e-9, "Power formula check");

const tem11Params = computeBeamParameters({ mode: "TEM11", w0: 1.0, I0: 100, wavelengthNm: 632.8 });
assert.strictEqual(tem11Params.mx2, 3, "TEM11 Mx^2 should be 3");
assert.strictEqual(tem11Params.my2, 3, "TEM11 My^2 should be 3");

// 6. Verify slice coordinates
assert.deepStrictEqual(getSliceCoordinates(1, 1, 1.0, "center"), { xSlice: 0, ySlice: 0 });
const peakSlice = getSliceCoordinates(1, 1, 1.0, "peak");
approxEqual(peakSlice.xSlice, Math.SQRT1_2, 1e-12, "TEM11 peak slice x");
approxEqual(peakSlice.ySlice, Math.SQRT1_2, 1e-12, "TEM11 peak slice y");

console.log("\nAll mathematical verification tests PASSED!");
