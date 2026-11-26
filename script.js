/**
 * ============================================================================
 * Laser Beam Profiling and Spatial Mode Selection
 * ============================================================================
 * Educational simulation implementing Hermite-Gaussian (TEMmn) laser beam
 * transverse intensity profiles, 2D heatmaps, interactive 3D surface plots,
 * horizontal/vertical cross-sections, and multi-mode comparison analysis.
 *
 * Units used consistently across the simulation:
 *  - Transverse coordinates (x, y) : millimeters (mm)
 *  - Beam waist (w0)               : millimeters (mm)
 *  - Reference intensity (I0)      : Watts per square centimeter (W/cm^2)
 *  - Wavelength (lambda)           : nanometers (nm)
 * ============================================================================
 */

"use strict";

/* --------------------------------------------------------------------------
   1. MODE DEFINITIONS & DEFAULT STATE
   -------------------------------------------------------------------------- */
const MODE_DEFINITIONS = {
  TEM00: {
    label: "TEM00",
    m: 0,
    n: 0,
    lobes: 1,
    nodalLines: "None",
    color: "#2563eb",
    description:
      "Fundamental Gaussian mode (m=0, n=0) with a single central peak, circular symmetry, and diffraction-limited beam quality (M² = 1)."
  },
  TEM01: {
    label: "TEM01",
    m: 0,
    n: 1,
    lobes: 2,
    nodalLines: "Horizontal (y = 0)",
    color: "#dc2626",
    description:
      "First-order vertical mode (m=0, n=1) featuring two symmetric lobes along the y-axis separated by a dark horizontal nodal line at y = 0."
  },
  TEM10: {
    label: "TEM10",
    m: 1,
    n: 0,
    lobes: 2,
    nodalLines: "Vertical (x = 0)",
    color: "#059669",
    description:
      "First-order horizontal mode (m=1, n=0) featuring two symmetric lobes along the x-axis separated by a dark vertical nodal line at x = 0."
  },
  TEM11: {
    label: "TEM11",
    m: 1,
    n: 1,
    lobes: 4,
    nodalLines: "x = 0 and y = 0",
    color: "#9333ea",
    description:
      "Higher-order quadrant mode (m=1, n=1) with four equal-intensity lobes divided by perpendicular nodal lines along x = 0 and y = 0."
  }
};

const DEFAULT_STATE = {
  mode: "TEM00",
  w0: 1.0, // Beam waist radius in mm
  I0: 100.0, // Reference TEM00 peak intensity in W/cm^2
  wavelengthNm: 632.8, // He-Ne default wavelength in nm
  sliceMode: "peak", // "peak" (through peak lobe) or "center" (x=0, y=0)
  scaleMode: "absolute", // "absolute" (W/cm^2) or "normalized" (I / I_peak)
  colormap: "inferno", // "inferno" | "laser" | "viridis"
  showGaussianRef: true,
  compareAxis: "peakSlice", // "peakSlice" | "horizontal" | "diagonal"
  azimuthDeg: 38,
  elevationDeg: 30
};

// Fixed physical half-width (in mm) for 2D/3D plots and cross-sections so
// changing beam waist w0 visibly expands or focuses the laser spot.
const PLOT_HALF_RANGE_MM = 5.0;

// Application state (cloned from defaults)
let appState = { ...DEFAULT_STATE };

// Chart.js instances
let horizontalChart = null;
let verticalChart = null;
let comparisonChart = null;

// Cache for 2D heatmap hover inspection
let lastHeatmapGeometry = null;

/* --------------------------------------------------------------------------
   2. MATHEMATICAL & PHYSICS ENGINE (PURE FUNCTIONS)
   -------------------------------------------------------------------------- */

/**
 * Computes factorial k! for small non-negative integers.
 * @param {number} k
 * @returns {number}
 */
function factorial(k) {
  let result = 1;
  for (let i = 2; i <= k; i++) {
    result *= i;
  }
  return result;
}

/**
 * Evaluates the physicist's Hermite polynomial H_order(u)
 * using the standard three-term recurrence relation:
 *   H_0(u) = 1
 *   H_1(u) = 2u
 *   H_k(u) = 2u * H_{k-1}(u) - 2(k-1) * H_{k-2}(u)
 * @param {number} order - Non-negative integer mode index (m or n)
 * @param {number} u - Dimensionless scaled coordinate
 * @returns {number}
 */
function hermitePolynomial(order, u) {
  if (order === 0) return 1;
  if (order === 1) return 2 * u;

  let hPrev2 = 1;
  let hPrev1 = 2 * u;
  for (let k = 2; k <= order; k++) {
    const hCurr = 2 * u * hPrev1 - 2 * (k - 1) * hPrev2;
    hPrev2 = hPrev1;
    hPrev1 = hCurr;
  }
  return hPrev1;
}

/**
 * Calculates the 2D Hermite-Gaussian intensity I_mn(x, y) at the beam waist (z = 0)
 * under equal total optical power normalization across modes:
 *
 *   I_mn(x, y) = [I0 / (2^(m+n) * m! * n!)]
 *                * [H_m(sqrt(2)*x/w0) * H_n(sqrt(2)*y/w0)]^2
 *                * exp(-2*(x^2 + y^2) / w0^2)
 *
 * @param {number} x - Horizontal coordinate in mm
 * @param {number} y - Vertical coordinate in mm
 * @param {number} w0 - Beam waist radius in mm
 * @param {number} I0 - Reference TEM00 peak intensity in W/cm^2
 * @param {number} m - Horizontal mode index
 * @param {number} n - Vertical mode index
 * @returns {number} Intensity in W/cm^2
 */
function calculateHermiteGaussianIntensity(x, y, w0, I0, m, n) {
  const sqrt2 = Math.SQRT2;
  const ux = (sqrt2 * x) / w0;
  const uy = (sqrt2 * y) / w0;

  const hx = hermitePolynomial(m, ux);
  const hy = hermitePolynomial(n, uy);

  const normFactor = 1 / (Math.pow(2, m + n) * factorial(m) * factorial(n));
  const gaussianEnvelope = Math.exp((-2 * (x * x + y * y)) / (w0 * w0));

  return I0 * normFactor * hx * hx * hy * hy * gaussianEnvelope;
}

/**
 * Returns the exact analytical peak intensity ratio (I_peak / I0) for modes
 * with m, n in {0, 1}.
 *  - index 0 -> peak factor 1
 *  - index 1 -> peak factor 2 * e^(-1) ≈ 0.73575888 (occurs at ±w0/sqrt(2))
 * @param {number} m
 * @param {number} n
 * @returns {number}
 */
function getTheoreticalPeakRatio(m, n) {
  const factor1D = (idx) => {
    if (idx === 0) return 1;
    if (idx === 1) return 2 * Math.exp(-1);
    // Numerical fallback for any higher order if ever extended
    let maxVal = 0;
    const norm = 1 / (Math.pow(2, idx) * factorial(idx));
    for (let u = 0; u <= 4; u += 0.001) {
      const h = hermitePolynomial(idx, u);
      const val = norm * h * h * Math.exp(-u * u);
      if (val > maxVal) maxVal = val;
    }
    return maxVal;
  };
  return factor1D(m) * factor1D(n);
}

/**
 * Determines the coordinate (xSlice, ySlice) where cross-sections are taken.
 * @param {number} m
 * @param {number} n
 * @param {number} w0
 * @param {string} sliceMode - "peak" or "center"
 * @returns {{xSlice: number, ySlice: number}}
 */
function getSliceCoordinates(m, n, w0, sliceMode) {
  if (sliceMode === "center") {
    return { xSlice: 0, ySlice: 0 };
  }
  // For index 1, peak lobe lies at +w0 / sqrt(2)
  const xPeak = m > 0 ? w0 / Math.SQRT2 : 0;
  const yPeak = n > 0 ? w0 / Math.SQRT2 : 0;
  return { xSlice: xPeak, ySlice: yPeak };
}

/**
 * Computes all derived physical parameters for the current laser configuration.
 * @param {object} state
 * @returns {object}
 */
function computeBeamParameters(state) {
  const modeInfo = MODE_DEFINITIONS[state.mode];
  const { m, n } = modeInfo;
  const { w0, I0, wavelengthNm } = state;

  const peakRatio = getTheoreticalPeakRatio(m, n);
  const peakIntensity = I0 * peakRatio;

  // Total optical power P = (1/2) * pi * w0^2 * I0
  // Convert w0 from mm to cm (1 mm = 0.1 cm) so P is in Watts.
  const w0Cm = w0 * 0.1;
  const totalPowerW = 0.5 * Math.PI * w0Cm * w0Cm * I0;

  // 1/e^2 diameter of fundamental Gaussian envelope
  const diameterMm = 2 * w0;

  // Full Width at Half Maximum (FWHM) of TEM00 Gaussian profile
  const fwhmMm = w0 * Math.sqrt(2 * Math.LN2);

  // Wavelength in mm (1 nm = 1e-6 mm)
  const wavelengthMm = wavelengthNm * 1e-6;

  // Rayleigh range z_R = pi * w0^2 / lambda (converted from mm to meters)
  const rayleighRangeM = (Math.PI * w0 * w0) / wavelengthMm / 1000;

  // Far-field half-angle divergence theta = lambda / (pi * w0) in milliradians
  const divergenceMrad = (wavelengthMm / (Math.PI * w0)) * 1000;

  // Beam quality factors for Hermite-Gaussian modes
  const mx2 = 2 * m + 1;
  const my2 = 2 * n + 1;

  return {
    modeInfo,
    peakRatio,
    peakIntensity,
    totalPowerW,
    diameterMm,
    fwhmMm,
    rayleighRangeM,
    divergenceMrad,
    mx2,
    my2
  };
}

/* --------------------------------------------------------------------------
   3. SCIENTIFIC COLORMAPS (0.0 -> 1.0 -> [R, G, B])
   -------------------------------------------------------------------------- */
const COLORMAP_STOPS = {
  inferno: [
    [0.0, [10, 8, 22]],
    [0.2, [66, 10, 104]],
    [0.45, [165, 44, 96]],
    [0.7, [237, 105, 37]],
    [0.9, [251, 193, 39]],
    [1.0, [252, 255, 164]]
  ],
  laser: [
    [0.0, [8, 8, 12]],
    [0.25, [90, 0, 18]],
    [0.55, [200, 16, 36]],
    [0.8, [255, 92, 48]],
    [0.95, [255, 200, 140]],
    [1.0, [255, 255, 255]]
  ],
  viridis: [
    [0.0, [68, 1, 84]],
    [0.25, [59, 82, 139]],
    [0.5, [33, 145, 140]],
    [0.75, [94, 201, 98]],
    [1.0, [253, 231, 37]]
  ]
};

/**
 * Linearly interpolates a normalized value t in [0, 1] into an RGB triplet.
 * @param {number} t
 * @param {string} paletteName
 * @returns {[number, number, number]}
 */
function sampleColormap(t, paletteName) {
  const clamped = Math.max(0, Math.min(1, t));
  const stops = COLORMAP_STOPS[paletteName] || COLORMAP_STOPS.inferno;

  for (let i = 0; i < stops.length - 1; i++) {
    const [t0, c0] = stops[i];
    const [t1, c1] = stops[i + 1];
    if (clamped >= t0 && clamped <= t1) {
      const f = (clamped - t0) / (t1 - t0 || 1);
      return [
        Math.round(c0[0] + f * (c1[0] - c0[0])),
        Math.round(c0[1] + f * (c1[1] - c0[1])),
        Math.round(c0[2] + f * (c1[2] - c0[2]))
      ];
    }
  }
  return stops[stops.length - 1][1];
}

/* --------------------------------------------------------------------------
   4. HTML5 CANVAS RENDERERS (2D HEATMAP, COLORBAR, MINI PREVIEWS, 3D SURFACE)
   -------------------------------------------------------------------------- */

/**
 * Renders a 2D beam intensity heatmap onto any square canvas element.
 * @param {HTMLCanvasElement} canvas
 * @param {object} options
 */
function drawHeatmapToCanvas(canvas, options) {
  const {
    modeKey,
    w0,
    I0,
    colormap,
    showAxes = true,
    showOverlayGuides = true,
    sliceMode = "peak"
  } = options;

  const ctx = canvas.getContext("2d");
  const width = canvas.width;
  const height = canvas.height;
  const { m, n } = MODE_DEFINITIONS[modeKey];
  const peakIntensity = I0 * getTheoreticalPeakRatio(m, n);

  // Clear background
  ctx.fillStyle = "#0f172a";
  ctx.fillRect(0, 0, width, height);

  // Define plot area margins for axes if enabled
  const padLeft = showAxes ? 44 : 6;
  const padRight = showAxes ? 14 : 6;
  const padTop = showAxes ? 14 : 6;
  const padBottom = showAxes ? 36 : 6;

  const plotW = width - padLeft - padRight;
  const plotH = height - padTop - padBottom;

  // Create pixel buffer for the intensity map
  const imgData = ctx.createImageData(plotW, plotH);
  const data = imgData.data;

  for (let py = 0; py < plotH; py++) {
    // Top of canvas corresponds to +PLOT_HALF_RANGE_MM, bottom to -PLOT_HALF_RANGE_MM
    const yMm = PLOT_HALF_RANGE_MM - (py / (plotH - 1)) * (2 * PLOT_HALF_RANGE_MM);

    for (let px = 0; px < plotW; px++) {
      const xMm = -PLOT_HALF_RANGE_MM + (px / (plotW - 1)) * (2 * PLOT_HALF_RANGE_MM);
      const intensity = calculateHermiteGaussianIntensity(xMm, yMm, w0, I0, m, n);
      const norm = peakIntensity > 0 ? intensity / peakIntensity : 0;
      const [r, g, b] = sampleColormap(norm, colormap);

      const idx = (py * plotW + px) * 4;
      data[idx] = r;
      data[idx + 1] = g;
      data[idx + 2] = b;
      data[idx + 3] = 255;
    }
  }

  ctx.putImageData(imgData, padLeft, padTop);

  // Helper converters between mm and canvas pixels
  const mmToPxX = (xMm) =>
    padLeft + ((xMm + PLOT_HALF_RANGE_MM) / (2 * PLOT_HALF_RANGE_MM)) * plotW;
  const mmToPxY = (yMm) =>
    padTop + ((PLOT_HALF_RANGE_MM - yMm) / (2 * PLOT_HALF_RANGE_MM)) * plotH;

  if (showOverlayGuides) {
    ctx.save();
    ctx.beginPath();
    ctx.rect(padLeft, padTop, plotW, plotH);
    ctx.clip();

    // Draw dashed 1/e^2 beam waist circle (radius = w0)
    const cx = mmToPxX(0);
    const cy = mmToPxY(0);
    const radiusPx = (w0 / (2 * PLOT_HALF_RANGE_MM)) * plotW;

    ctx.strokeStyle = "rgba(255, 255, 255, 0.72)";
    ctx.lineWidth = 1.4;
    ctx.setLineDash([5, 4]);
    ctx.beginPath();
    ctx.arc(cx, cy, radiusPx, 0, Math.PI * 2);
    ctx.stroke();

    // Draw active horizontal & vertical slice guide lines
    const { xSlice, ySlice } = getSliceCoordinates(m, n, w0, sliceMode);
    ctx.strokeStyle = "rgba(56, 189, 248, 0.55)";
    ctx.lineWidth = 1;
    ctx.setLineDash([3, 3]);

    // Horizontal slice line at y = ySlice
    const sliceYPx = mmToPxY(ySlice);
    ctx.beginPath();
    ctx.moveTo(padLeft, sliceYPx);
    ctx.lineTo(padLeft + plotW, sliceYPx);
    ctx.stroke();

    // Vertical slice line at x = xSlice
    const sliceXPx = mmToPxX(xSlice);
    ctx.beginPath();
    ctx.moveTo(sliceXPx, padTop);
    ctx.lineTo(sliceXPx, padTop + plotH);
    ctx.stroke();

    ctx.restore();
  }

  // Border around plot
  ctx.strokeStyle = "#475569";
  ctx.lineWidth = 1;
  ctx.strokeRect(padLeft, padTop, plotW, plotH);

  // Draw axis ticks and labels
  if (showAxes) {
    ctx.fillStyle = "#cbd5e1";
    ctx.font = "11px Segoe UI, sans-serif";

    const ticks = [-4, -2, 0, 2, 4];
    ctx.textAlign = "center";
    ctx.textBaseline = "top";
    for (const t of ticks) {
      const x = mmToPxX(t);
      ctx.beginPath();
      ctx.moveTo(x, padTop + plotH);
      ctx.lineTo(x, padTop + plotH + 4);
      ctx.strokeStyle = "#94a3b8";
      ctx.stroke();
      ctx.fillText(t.toString(), x, padTop + plotH + 6);
    }

    ctx.textAlign = "right";
    ctx.textBaseline = "middle";
    for (const t of ticks) {
      const y = mmToPxY(t);
      ctx.beginPath();
      ctx.moveTo(padLeft - 4, y);
      ctx.lineTo(padLeft, y);
      ctx.strokeStyle = "#94a3b8";
      ctx.stroke();
      ctx.fillText(t.toString(), padLeft - 6, y);
    }

    // Axis titles
    ctx.fillStyle = "#e2e8f0";
    ctx.font = "600 11px Segoe UI, sans-serif";
    ctx.textAlign = "center";
    ctx.textBaseline = "bottom";
    ctx.fillText("x (mm)", padLeft + plotW / 2, height - 3);

    ctx.save();
    ctx.translate(12, padTop + plotH / 2);
    ctx.rotate(-Math.PI / 2);
    ctx.textAlign = "center";
    ctx.textBaseline = "top";
    ctx.fillText("y (mm)", 0, 0);
    ctx.restore();

    // Save geometry for mouse hover inspection
    lastHeatmapGeometry = {
      padLeft,
      padTop,
      plotW,
      plotH,
      width,
      height
    };
  }
}

/**
 * Draws the vertical calibrated colorbar beside the main 2D heatmap.
 */
function renderColorbar(peakIntensity, scaleMode, colormap) {
  const cbCanvas = document.getElementById("colorbarCanvas");
  const ctx = cbCanvas.getContext("2d");
  const w = cbCanvas.width;
  const h = cbCanvas.height;

  for (let py = 0; py < h; py++) {
    const norm = 1 - py / (h - 1);
    const [r, g, b] = sampleColormap(norm, colormap);
    ctx.fillStyle = `rgb(${r}, ${g}, ${b})`;
    ctx.fillRect(0, py, w, 1);
  }

  const maxLabel = document.getElementById("colorbarMax");
  maxLabel.textContent =
    scaleMode === "normalized" ? "1.00" : peakIntensity.toFixed(1);
}

/**
 * Renders an interactive 3D-style intensity surface plot using pure HTML5 Canvas.
 * Uses a painter's algorithm over a 36x36 spatial mesh with perspective projection.
 */
function renderSurface3D() {
  const canvas = document.getElementById("surface3dCanvas");
  const ctx = canvas.getContext("2d");
  const w = canvas.width;
  const h = canvas.height;

  const { m, n } = MODE_DEFINITIONS[appState.mode];
  const { w0, I0, colormap, azimuthDeg, elevationDeg } = appState;
  const peakIntensity = I0 * getTheoreticalPeakRatio(m, n);

  // Clear background with subtle laboratory gradient
  const bgGrad = ctx.createLinearGradient(0, 0, 0, h);
  bgGrad.addColorStop(0, "#0b1120");
  bgGrad.addColorStop(1, "#1e293b");
  ctx.fillStyle = bgGrad;
  ctx.fillRect(0, 0, w, h);

  const gridN = 36;
  const az = (azimuthDeg * Math.PI) / 180;
  const el = (elevationDeg * Math.PI) / 180;

  const cosAz = Math.cos(az);
  const sinAz = Math.sin(az);
  const cosEl = Math.cos(el);
  const sinEl = Math.sin(el);

  // Projects normalized 3D coordinates (nx, ny in [-1, 1], nz in [0, 1]) to 2D canvas
  function projectPoint(nx, ny, nz) {
    // Center height around 0.4 for balanced framing
    const zCentered = (nz - 0.35) * 1.15;

    // Rotate around vertical Z-axis by azimuth
    const rx = nx * cosAz - ny * sinAz;
    const ry = nx * sinAz + ny * cosAz;

    // Tilt around X-axis by elevation
    const yTilt = ry * cosEl - zCentered * sinEl;
    const zTilt = ry * sinEl + zCentered * cosEl;

    // Weak perspective projection
    const cameraDist = 4.2;
    const scale = (cameraDist / (cameraDist + yTilt)) * 128;

    return {
      sx: w * 0.5 + rx * scale,
      sy: h * 0.56 - zTilt * scale,
      depth: yTilt
    };
  }

  // Evaluate grid vertices
  const vertices = [];
  for (let i = 0; i <= gridN; i++) {
    const row = [];
    const nx = -1 + (2 * i) / gridN;
    const xMm = nx * PLOT_HALF_RANGE_MM;
    for (let j = 0; j <= gridN; j++) {
      const ny = -1 + (2 * j) / gridN;
      const yMm = ny * PLOT_HALF_RANGE_MM;
      const intensity = calculateHermiteGaussianIntensity(xMm, yMm, w0, I0, m, n);
      const norm = peakIntensity > 0 ? intensity / peakIntensity : 0;
      const proj = projectPoint(nx, ny, norm);
      row.push({ ...proj, norm });
    }
    vertices.push(row);
  }

  // Build quads and sort back-to-front (painter's algorithm)
  const quads = [];
  for (let i = 0; i < gridN; i++) {
    for (let j = 0; j < gridN; j++) {
      const p00 = vertices[i][j];
      const p10 = vertices[i + 1][j];
      const p11 = vertices[i + 1][j + 1];
      const p01 = vertices[i][j + 1];

      const avgDepth = (p00.depth + p10.depth + p11.depth + p01.depth) * 0.25;
      const avgNorm = (p00.norm + p10.norm + p11.norm + p01.norm) * 0.25;

      quads.push({ p00, p10, p11, p01, avgDepth, avgNorm });
    }
  }

  quads.sort((a, b) => b.avgDepth - a.avgDepth);

  // Draw base reference bounding box lines at z = 0
  const b00 = projectPoint(-1, -1, 0);
  const b10 = projectPoint(1, -1, 0);
  const b11 = projectPoint(1, 1, 0);
  const b01 = projectPoint(-1, 1, 0);
  const topZ = projectPoint(-1, -1, 1);

  ctx.strokeStyle = "rgba(148, 163, 184, 0.35)";
  ctx.lineWidth = 1;
  ctx.beginPath();
  ctx.moveTo(b00.sx, b00.sy);
  ctx.lineTo(b10.sx, b10.sy);
  ctx.lineTo(b11.sx, b11.sy);
  ctx.lineTo(b01.sx, b01.sy);
  ctx.closePath();
  ctx.moveTo(b00.sx, b00.sy);
  ctx.lineTo(topZ.sx, topZ.sy);
  ctx.stroke();

  // Draw sorted surface quads
  for (const q of quads) {
    const [r, g, b] = sampleColormap(q.avgNorm, colormap);
    ctx.fillStyle = `rgb(${r}, ${g}, ${b})`;
    ctx.strokeStyle = "rgba(15, 23, 42, 0.35)";
    ctx.lineWidth = 0.6;

    ctx.beginPath();
    ctx.moveTo(q.p00.sx, q.p00.sy);
    ctx.lineTo(q.p10.sx, q.p10.sy);
    ctx.lineTo(q.p11.sx, q.p11.sy);
    ctx.lineTo(q.p01.sx, q.p01.sy);
    ctx.closePath();
    ctx.fill();
    ctx.stroke();
  }

  // Draw axis labels on 3D canvas
  const xLabelPos = projectPoint(0, -1.16, 0);
  const yLabelPos = projectPoint(1.16, 0, 0);
  const zLabelPos = projectPoint(-1.08, -1.08, 1.05);

  ctx.fillStyle = "#cbd5e1";
  ctx.font = "600 11px Segoe UI, sans-serif";
  ctx.textAlign = "center";
  ctx.fillText("x (±5 mm)", xLabelPos.sx, xLabelPos.sy);
  ctx.fillText("y (±5 mm)", yLabelPos.sx, yLabelPos.sy);
  ctx.fillText("I(x,y)", zLabelPos.sx, zLabelPos.sy);

  // Update angle text readout
  const angleReadout = document.getElementById("angleReadout");
  angleReadout.textContent = `Azimuth ${Math.round(azimuthDeg)}°, Elevation ${Math.round(elevationDeg)}°`;
}

/* --------------------------------------------------------------------------
   5. CHART.JS CROSS-SECTION & COMPARISON GRAPHS
   -------------------------------------------------------------------------- */

/**
 * Generates sample coordinates from -PLOT_HALF_RANGE_MM to +PLOT_HALF_RANGE_MM.
 * @param {number} numPoints
 * @returns {number[]}
 */
function generateAxisPoints(numPoints = 121) {
  const pts = [];
  for (let i = 0; i < numPoints; i++) {
    const val = -PLOT_HALF_RANGE_MM + (i / (numPoints - 1)) * (2 * PLOT_HALF_RANGE_MM);
    pts.push(Number(val.toFixed(3)));
  }
  return pts;
}

/**
 * Draws a simple fallback message if Chart.js is unavailable (e.g. offline without cache).
 */
function drawChartFallback(canvasId, message) {
  const canvas = document.getElementById(canvasId);
  if (!canvas) return;
  const ctx = canvas.getContext("2d");
  ctx.clearRect(0, 0, canvas.width, canvas.height);
  ctx.fillStyle = "#64748b";
  ctx.font = "13px Segoe UI, sans-serif";
  ctx.textAlign = "center";
  ctx.fillText(message, canvas.width / 2, canvas.height / 2);
}

/**
 * Updates both Horizontal I(x) and Vertical I(y) cross-section Chart.js graphs.
 */
function updateCrossSectionCharts() {
  const coords = generateAxisPoints(121);
  const labels = coords.map((v) => v.toFixed(2));

  const { mode, w0, I0, sliceMode, scaleMode, showGaussianRef } = appState;
  const { m, n, color, label } = MODE_DEFINITIONS[mode];
  const { xSlice, ySlice } = getSliceCoordinates(m, n, w0, sliceMode);
  const peakIntensity = I0 * getTheoreticalPeakRatio(m, n);

  const isNorm = scaleMode === "normalized";
  const yAxisTitle = isNorm ? "Normalized Intensity (I / I_peak)" : "Intensity (W/cm²)";
  const yMax = isNorm ? 1.05 : Number((I0 * 1.08).toFixed(2));

  // Horizontal slice values: vary x, fix y = ySlice
  const horizValues = coords.map((x) => {
    const val = calculateHermiteGaussianIntensity(x, ySlice, w0, I0, m, n);
    return isNorm ? (peakIntensity > 0 ? val / peakIntensity : 0) : val;
  });

  // Vertical slice values: fix x = xSlice, vary y
  const vertValues = coords.map((y) => {
    const val = calculateHermiteGaussianIntensity(xSlice, y, w0, I0, m, n);
    return isNorm ? (peakIntensity > 0 ? val / peakIntensity : 0) : val;
  });

  // Fundamental TEM00 reference curve along center axis
  const refValues = coords.map((coord) => {
    const val = calculateHermiteGaussianIntensity(coord, 0, w0, I0, 0, 0);
    return isNorm ? val / I0 : val;
  });

  // Update descriptive subtitles
  document.getElementById("horizSubtitle").innerHTML = `Slice along <em>y</em> = ${ySlice.toFixed(2)} mm`;
  document.getElementById("vertSubtitle").innerHTML = `Slice along <em>x</em> = ${xSlice.toFixed(2)} mm`;

  if (typeof Chart === "undefined") {
    drawChartFallback("horizontalChart", "Chart.js library not loaded (check network connection).");
    drawChartFallback("verticalChart", "Chart.js library not loaded (check network connection).");
    return;
  }

  const buildDatasets = (primaryData) => {
    const sets = [
      {
        label: `${label} Profile`,
        data: primaryData,
        borderColor: color,
        backgroundColor: `${color}22`,
        borderWidth: 2.5,
        fill: true,
        pointRadius: 0,
        tension: 0.25
      }
    ];
    if (showGaussianRef && mode !== "TEM00") {
      sets.push({
        label: "TEM00 Gaussian Ref.",
        data: refValues,
        borderColor: "#64748b",
        borderDash: [6, 4],
        borderWidth: 1.75,
        fill: false,
        pointRadius: 0,
        tension: 0.25
      });
    }
    return sets;
  };

  const buildOptions = (xTitle) => ({
    responsive: true,
    maintainAspectRatio: false,
    animation: { duration: 150 },
    interaction: { mode: "index", intersect: false },
    plugins: {
      legend: { position: "top", labels: { boxWidth: 14, font: { size: 11 } } },
      tooltip: {
        callbacks: {
          title: (items) => `${xTitle.split(" ")[0]} = ${items[0].label} mm`,
          label: (ctx) => {
            const unit = isNorm ? "" : " W/cm²";
            return `${ctx.dataset.label}: ${ctx.parsed.y.toFixed(3)}${unit}`;
          }
        }
      }
    },
    scales: {
      x: {
        title: { display: true, text: xTitle, font: { size: 11, weight: "600" } },
        ticks: { maxTicksLimit: 9, font: { size: 10 } }
      },
      y: {
        min: 0,
        max: yMax,
        title: { display: true, text: yAxisTitle, font: { size: 11, weight: "600" } },
        ticks: { font: { size: 10 } }
      }
    }
  });

  // Update or create Horizontal Chart
  const horizCtx = document.getElementById("horizontalChart").getContext("2d");
  if (horizontalChart) {
    horizontalChart.data.labels = labels;
    horizontalChart.data.datasets = buildDatasets(horizValues);
    horizontalChart.options = buildOptions("Horizontal Position x (mm)");
    horizontalChart.update();
  } else {
    horizontalChart = new Chart(horizCtx, {
      type: "line",
      data: { labels, datasets: buildDatasets(horizValues) },
      options: buildOptions("Horizontal Position x (mm)")
    });
  }

  // Update or create Vertical Chart
  const vertCtx = document.getElementById("verticalChart").getContext("2d");
  if (verticalChart) {
    verticalChart.data.labels = labels;
    verticalChart.data.datasets = buildDatasets(vertValues);
    verticalChart.options = buildOptions("Vertical Position y (mm)");
    verticalChart.update();
  } else {
    verticalChart = new Chart(vertCtx, {
      type: "line",
      data: { labels, datasets: buildDatasets(vertValues) },
      options: buildOptions("Vertical Position y (mm)")
    });
  }
}

/**
 * Updates the 4 mini mode preview heatmaps, overlay comparison chart,
 * and quantitative mode comparison table.
 */
function updateComparisonPanel() {
  const { w0, I0, colormap, scaleMode, compareAxis, mode } = appState;
  const isNorm = scaleMode === "normalized";
  const modeKeys = ["TEM00", "TEM01", "TEM10", "TEM11"];

  // 1. Render mini heatmaps and highlight active mode card
  for (const key of modeKeys) {
    const miniCanvas = document.getElementById(`miniCanvas${key}`);
    const { m, n } = MODE_DEFINITIONS[key];
    const peak = I0 * getTheoreticalPeakRatio(m, n);

    drawHeatmapToCanvas(miniCanvas, {
      modeKey: key,
      w0,
      I0,
      colormap,
      showAxes: false,
      showOverlayGuides: false
    });

    const peakSpan = document.getElementById(`miniPeak${key}`);
    peakSpan.textContent = isNorm
      ? `Ratio ${(peak / I0).toFixed(3)}`
      : `${peak.toFixed(1)} W/cm²`;

    const card = miniCanvas.closest(".mini-mode-card");
    card.classList.toggle("active-compare", key === mode);
  }

  // 2. Populate comparison table
  const tbody = document.getElementById("comparisonTableBody");
  tbody.innerHTML = "";
  for (const key of modeKeys) {
    const info = MODE_DEFINITIONS[key];
    const ratio = getTheoreticalPeakRatio(info.m, info.n);
    const peak = I0 * ratio;

    const tr = document.createElement("tr");
    if (key === mode) tr.classList.add("active-row");
    tr.innerHTML = `
      <td><strong>${info.label}</strong></td>
      <td>(${info.m}, ${info.n})</td>
      <td>${info.lobes}</td>
      <td>${info.nodalLines}</td>
      <td>${ratio.toFixed(4)}</td>
      <td>${peak.toFixed(2)}</td>
      <td>${2 * info.m + 1}, ${2 * info.n + 1}</td>
    `;
    tbody.appendChild(tr);
  }

  // 3. Render overlay comparison Chart.js graph
  if (typeof Chart === "undefined") {
    drawChartFallback("comparisonChart", "Chart.js library not loaded (check network connection).");
    return;
  }

  const coords = generateAxisPoints(121);
  const labels = coords.map((v) => v.toFixed(2));
  const yAxisTitle = isNorm ? "Normalized Intensity (I / I_peak)" : "Intensity (W/cm²)";
  const yMax = isNorm ? 1.05 : Number((I0 * 1.08).toFixed(2));

  let xAxisTitle = "Coordinate r (mm)";
  const datasets = modeKeys.map((key) => {
    const { m, n, color, label } = MODE_DEFINITIONS[key];
    const peak = I0 * getTheoreticalPeakRatio(m, n);

    const data = coords.map((s) => {
      let x = s;
      let y = 0;

      if (compareAxis === "peakSlice") {
        // Slice along the principal axis of each mode through its peak lobe:
        // TEM01 varies along y (at x=0); TEM00, TEM10, TEM11 vary along x (at y=y_peak)
        xAxisTitle = "Position Along Principal Lobe Axis (mm)";
        if (key === "TEM01") {
          x = 0;
          y = s;
        } else {
          x = s;
          y = n > 0 ? w0 / Math.SQRT2 : 0;
        }
      } else if (compareAxis === "horizontal") {
        xAxisTitle = "Horizontal Position x at y = 0 (mm)";
        x = s;
        y = 0;
      } else {
        // Diagonal y = x (converted to radial distance s along diagonal)
        xAxisTitle = "Diagonal Coordinate d (x = y = d / √2) (mm)";
        x = s / Math.SQRT2;
        y = s / Math.SQRT2;
      }

      const intensity = calculateHermiteGaussianIntensity(x, y, w0, I0, m, n);
      return isNorm ? (peak > 0 ? intensity / peak : 0) : intensity;
    });

    return {
      label,
      data,
      borderColor: color,
      borderWidth: key === mode ? 3 : 2,
      borderDash: key === "TEM01" && compareAxis === "peakSlice" ? [5, 3] : [],
      fill: false,
      pointRadius: 0,
      tension: 0.25
    };
  });

  const compOptions = {
    responsive: true,
    maintainAspectRatio: false,
    animation: { duration: 150 },
    interaction: { mode: "index", intersect: false },
    plugins: {
      legend: { position: "top", labels: { boxWidth: 14, font: { size: 11 } } }
    },
    scales: {
      x: {
        title: { display: true, text: xAxisTitle, font: { size: 11, weight: "600" } },
        ticks: { maxTicksLimit: 9, font: { size: 10 } }
      },
      y: {
        min: 0,
        max: yMax,
        title: { display: true, text: yAxisTitle, font: { size: 11, weight: "600" } },
        ticks: { font: { size: 10 } }
      }
    }
  };

  const compCtx = document.getElementById("comparisonChart").getContext("2d");
  if (comparisonChart) {
    comparisonChart.data.labels = labels;
    comparisonChart.data.datasets = datasets;
    comparisonChart.options = compOptions;
    comparisonChart.update();
  } else {
    comparisonChart = new Chart(compCtx, {
      type: "line",
      data: { labels, datasets },
      options: compOptions
    });
  }
}

/* --------------------------------------------------------------------------
   6. DOM READOUT & FULL SIMULATION REFRESH
   -------------------------------------------------------------------------- */

/**
 * Updates all parameter readouts, mode descriptions, and live analysis text.
 */
function updateParameterReadouts() {
  const params = computeBeamParameters(appState);
  const { modeInfo } = params;

  document.getElementById("modeQuickDesc").textContent = modeInfo.description;
  document.getElementById("heatmapModeTag").innerHTML = `TEM<sub>${modeInfo.m}${modeInfo.n}</sub>`;

  document.getElementById("paramPeakIntensity").textContent = `${params.peakIntensity.toFixed(2)} W/cm²`;
  document.getElementById("paramNormPeak").textContent = params.peakRatio.toFixed(3);
  document.getElementById("paramPower").textContent = `${params.totalPowerW.toFixed(3)} W`;
  document.getElementById("paramDiameter").textContent = `${params.diameterMm.toFixed(2)} mm`;
  document.getElementById("paramFwhm").textContent = `${params.fwhmMm.toFixed(3)} mm`;
  document.getElementById("paramRayleigh").textContent = `${params.rayleighRangeM.toFixed(3)} m`;
  document.getElementById("paramDivergence").textContent = `${params.divergenceMrad.toFixed(3)} mrad`;
  document.getElementById("paramM2").textContent = `${params.mx2}, ${params.my2}`;

  // Update dynamic Results & Discussion snapshot
  const liveBox = document.getElementById("liveAnalysisBox");
  liveBox.innerHTML = `
    <strong>Current Simulation Snapshot (${modeInfo.label}):</strong>
    With <em>w</em><sub>0</sub> = ${appState.w0.toFixed(2)} mm, <em>I</em><sub>0</sub> = ${appState.I0.toFixed(0)} W/cm&sup2;,
    and &lambda; = ${appState.wavelengthNm.toFixed(1)} nm, the beam carries <strong>${params.totalPowerW.toFixed(3)} W</strong>
    total power and reaches a peak irradiance of <strong>${params.peakIntensity.toFixed(2)} W/cm&sup2;</strong>
    (${(params.peakRatio * 100).toFixed(1)}% of fundamental TEM<sub>00</sub> peak).
    Its Rayleigh range is <strong>${params.rayleighRangeM.toFixed(3)} m</strong> with a far-field divergence of
    <strong>${params.divergenceMrad.toFixed(3)} mrad</strong>.
  `;
}

/**
 * Re-computes and re-renders every visual component on the page.
 * @param {string} statusMessage
 */
function runSimulation(statusMessage = "Updated") {
  updateParameterReadouts();

  const { m, n } = MODE_DEFINITIONS[appState.mode];
  const peakIntensity = appState.I0 * getTheoreticalPeakRatio(m, n);

  // Main 2D Heatmap & Colorbar
  const heatmapCanvas = document.getElementById("heatmapCanvas");
  drawHeatmapToCanvas(heatmapCanvas, {
    modeKey: appState.mode,
    w0: appState.w0,
    I0: appState.I0,
    colormap: appState.colormap,
    showAxes: true,
    showOverlayGuides: true,
    sliceMode: appState.sliceMode
  });
  renderColorbar(peakIntensity, appState.scaleMode, appState.colormap);

  // 3D Surface
  renderSurface3D();

  // Chart.js graphs & comparison panel
  updateCrossSectionCharts();
  updateComparisonPanel();

  // Status badge feedback
  const statusEl = document.getElementById("statusIndicator");
  statusEl.textContent = statusMessage;
}

/* --------------------------------------------------------------------------
   7. EVENT HANDLERS & INTERACTIVE BINDINGS
   -------------------------------------------------------------------------- */

/**
 * Synchronizes UI input controls with appState (used during init & reset).
 */
function syncControlsFromState() {
  document.querySelectorAll(".mode-btn").forEach((btn) => {
    const isActive = btn.dataset.mode === appState.mode;
    btn.classList.toggle("active", isActive);
    btn.setAttribute("aria-pressed", String(isActive));
  });

  document.getElementById("beamWaistRange").value = appState.w0.toFixed(2);
  document.getElementById("beamWaistInput").value = appState.w0.toFixed(2);

  document.getElementById("intensityRange").value = appState.I0.toFixed(0);
  document.getElementById("intensityInput").value = appState.I0.toFixed(0);

  document.getElementById("wavelengthSelect").value = appState.wavelengthNm.toFixed(1);
  document.getElementById("sliceModeSelect").value = appState.sliceMode;
  document.getElementById("scaleModeSelect").value = appState.scaleMode;
  document.getElementById("colormapSelect").value = appState.colormap;
  document.getElementById("showGaussianRef").checked = appState.showGaussianRef;
  document.getElementById("compareAxisSelect").value = appState.compareAxis;
}

/**
 * Clamps a numeric input within [min, max] safely.
 */
function clampValue(val, min, max, fallback) {
  if (Number.isNaN(val)) return fallback;
  return Math.max(min, Math.min(max, val));
}

function setupEventListeners() {
  // Spatial Mode buttons
  document.querySelectorAll(".mode-btn").forEach((btn) => {
    btn.addEventListener("click", () => {
      appState.mode = btn.dataset.mode;
      syncControlsFromState();
      runSimulation(`Mode: ${appState.mode}`);
    });
  });

  // Clicking a mini-mode card in Comparison section selects that mode
  document.querySelectorAll(".mini-mode-card").forEach((card) => {
    card.addEventListener("click", () => {
      appState.mode = card.dataset.compareMode;
      syncControlsFromState();
      runSimulation(`Selected ${appState.mode}`);
    });
  });

  // Beam Waist slider & numeric input
  const waistRange = document.getElementById("beamWaistRange");
  const waistInput = document.getElementById("beamWaistInput");

  waistRange.addEventListener("input", () => {
    appState.w0 = parseFloat(waistRange.value);
    waistInput.value = appState.w0.toFixed(2);
    runSimulation("Waist Updated");
  });

  waistInput.addEventListener("change", () => {
    appState.w0 = clampValue(parseFloat(waistInput.value), 0.2, 2.5, DEFAULT_STATE.w0);
    waistRange.value = appState.w0.toFixed(2);
    waistInput.value = appState.w0.toFixed(2);
    runSimulation("Waist Updated");
  });

  // Laser Intensity slider & numeric input
  const intensityRange = document.getElementById("intensityRange");
  const intensityInput = document.getElementById("intensityInput");

  intensityRange.addEventListener("input", () => {
    appState.I0 = parseFloat(intensityRange.value);
    intensityInput.value = appState.I0.toFixed(0);
    runSimulation("Intensity Updated");
  });

  intensityInput.addEventListener("change", () => {
    appState.I0 = clampValue(parseFloat(intensityInput.value), 10, 500, DEFAULT_STATE.I0);
    intensityRange.value = appState.I0.toFixed(0);
    intensityInput.value = appState.I0.toFixed(0);
    runSimulation("Intensity Updated");
  });

  // Wavelength, slice mode, scale mode, colormap, Gaussian ref
  document.getElementById("wavelengthSelect").addEventListener("change", (e) => {
    appState.wavelengthNm = parseFloat(e.target.value);
    runSimulation("Wavelength Updated");
  });

  document.getElementById("sliceModeSelect").addEventListener("change", (e) => {
    appState.sliceMode = e.target.value;
    runSimulation("Slice Updated");
  });

  document.getElementById("scaleModeSelect").addEventListener("change", (e) => {
    appState.scaleMode = e.target.value;
    runSimulation("Scale Updated");
  });

  document.getElementById("colormapSelect").addEventListener("change", (e) => {
    appState.colormap = e.target.value;
    runSimulation("Palette Updated");
  });

  document.getElementById("showGaussianRef").addEventListener("change", (e) => {
    appState.showGaussianRef = e.target.checked;
    runSimulation("Reference Toggled");
  });

  document.getElementById("compareAxisSelect").addEventListener("change", (e) => {
    appState.compareAxis = e.target.value;
    updateComparisonPanel();
  });

  // Action Buttons: Calculate, Reset, Compare Modes
  document.getElementById("calculateBtn").addEventListener("click", () => {
    // Re-read numeric fields in case user typed without blurring
    appState.w0 = clampValue(parseFloat(waistInput.value), 0.2, 2.5, appState.w0);
    appState.I0 = clampValue(parseFloat(intensityInput.value), 10, 500, appState.I0);
    syncControlsFromState();
    runSimulation("Calculated");
  });

  document.getElementById("resetBtn").addEventListener("click", () => {
    appState = { ...DEFAULT_STATE };
    syncControlsFromState();
    runSimulation("Reset to Defaults");
  });

  document.getElementById("compareBtn").addEventListener("click", () => {
    updateComparisonPanel();
    const compSection = document.getElementById("comparison-section");
    compSection.scrollIntoView({ behavior: "smooth", block: "start" });
    compSection.classList.add("highlight-pulse");
    setTimeout(() => compSection.classList.remove("highlight-pulse"), 1400);
    document.getElementById("statusIndicator").textContent = "Comparing All Modes";
  });

  // 2D Heatmap hover coordinate & intensity inspector
  const heatmapCanvas = document.getElementById("heatmapCanvas");
  const hoverReadout = document.getElementById("hoverReadout");

  heatmapCanvas.addEventListener("mousemove", (e) => {
    if (!lastHeatmapGeometry) return;
    const rect = heatmapCanvas.getBoundingClientRect();
    const scaleX = heatmapCanvas.width / rect.width;
    const scaleY = heatmapCanvas.height / rect.height;
    const px = (e.clientX - rect.left) * scaleX;
    const py = (e.clientY - rect.top) * scaleY;

    const { padLeft, padTop, plotW, plotH } = lastHeatmapGeometry;
    if (px < padLeft || px > padLeft + plotW || py < padTop || py > padTop + plotH) {
      hoverReadout.innerHTML = "Move cursor inside the plot area to inspect local intensity.";
      return;
    }

    const xMm = -PLOT_HALF_RANGE_MM + ((px - padLeft) / plotW) * (2 * PLOT_HALF_RANGE_MM);
    const yMm = PLOT_HALF_RANGE_MM - ((py - padTop) / plotH) * (2 * PLOT_HALF_RANGE_MM);

    const { m, n } = MODE_DEFINITIONS[appState.mode];
    const localI = calculateHermiteGaussianIntensity(xMm, yMm, appState.w0, appState.I0, m, n);
    const peakI = appState.I0 * getTheoreticalPeakRatio(m, n);
    const normI = peakI > 0 ? localI / peakI : 0;

    hoverReadout.innerHTML =
      `<em>x</em> = <strong>${xMm.toFixed(2)} mm</strong>, ` +
      `<em>y</em> = <strong>${yMm.toFixed(2)} mm</strong> &nbsp;|&nbsp; ` +
      `<em>I</em>(<em>x</em>,<em>y</em>) = <strong>${localI.toFixed(2)} W/cm&sup2;</strong> ` +
      `(Norm: <strong>${normI.toFixed(3)}</strong>)`;
  });

  heatmapCanvas.addEventListener("mouseleave", () => {
    hoverReadout.innerHTML =
      "Move cursor over heatmap to inspect (<em>x</em>, <em>y</em>) and local intensity.";
  });

  // 3D Surface drag-to-rotate controls
  const surfaceCanvas = document.getElementById("surface3dCanvas");
  let isDragging = false;
  let lastMouseX = 0;
  let lastMouseY = 0;

  surfaceCanvas.addEventListener("mousedown", (e) => {
    isDragging = true;
    lastMouseX = e.clientX;
    lastMouseY = e.clientY;
  });

  window.addEventListener("mousemove", (e) => {
    if (!isDragging) return;
    const dx = e.clientX - lastMouseX;
    const dy = e.clientY - lastMouseY;
    lastMouseX = e.clientX;
    lastMouseY = e.clientY;

    appState.azimuthDeg = (appState.azimuthDeg + dx * 0.6) % 360;
    appState.elevationDeg = clampValue(appState.elevationDeg + dy * 0.45, 10, 82, 30);
    renderSurface3D();
  });

  window.addEventListener("mouseup", () => {
    isDragging = false;
  });

  // Touch support for 3D rotation on mobile/tablets
  surfaceCanvas.addEventListener(
    "touchstart",
    (e) => {
      if (e.touches.length === 1) {
        isDragging = true;
        lastMouseX = e.touches[0].clientX;
        lastMouseY = e.touches[0].clientY;
      }
    },
    { passive: true }
  );

  surfaceCanvas.addEventListener(
    "touchmove",
    (e) => {
      if (!isDragging || e.touches.length !== 1) return;
      const dx = e.touches[0].clientX - lastMouseX;
      const dy = e.touches[0].clientY - lastMouseY;
      lastMouseX = e.touches[0].clientX;
      lastMouseY = e.touches[0].clientY;

      appState.azimuthDeg = (appState.azimuthDeg + dx * 0.6) % 360;
      appState.elevationDeg = clampValue(appState.elevationDeg + dy * 0.45, 10, 82, 30);
      renderSurface3D();
    },
    { passive: true }
  );

  surfaceCanvas.addEventListener("touchend", () => {
    isDragging = false;
  });

  document.getElementById("resetViewBtn").addEventListener("click", () => {
    appState.azimuthDeg = DEFAULT_STATE.azimuthDeg;
    appState.elevationDeg = DEFAULT_STATE.elevationDeg;
    renderSurface3D();
  });
}

/* --------------------------------------------------------------------------
   8. INITIALIZATION & TEST EXPORTS
   -------------------------------------------------------------------------- */
if (typeof document !== "undefined") {
  document.addEventListener("DOMContentLoaded", () => {
    syncControlsFromState();
    setupEventListeners();
    runSimulation("Ready");
  });
}

// Export pure physics functions when imported in Node.js for verification tests
if (typeof module !== "undefined" && module.exports) {
  module.exports = {
    MODE_DEFINITIONS,
    factorial,
    hermitePolynomial,
    calculateHermiteGaussianIntensity,
    getTheoreticalPeakRatio,
    getSliceCoordinates,
    computeBeamParameters
  };
}
