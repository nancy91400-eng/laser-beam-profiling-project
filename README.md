# Laser Beam Profiling and Spatial Mode Selection

An interactive B.Tech mini project built with **HTML5, CSS3, JavaScript, HTML5 Canvas, and Chart.js** to simulate and visualize Hermite-Gaussian ($\text{TEM}_{mn}$) laser beam intensity distributions and spatial modes.

> **Educational Note:** This is a mathematical physics simulation at the beam waist ($z = 0$) based on paraxial wave optics, not a real hardware laser experiment.

---

## Project Features

1. **Laser Beam Profiling**
   - **2D Intensity Heatmap (HTML5 Canvas):** High-resolution pixel-level rendering with selectable scientific colormaps (Inferno, Laser Crimson, Viridis), $1/e^2$ beam waist circle overlay, slice guides, and live mouse-hover $(x, y)$ intensity inspection.
   - **3D-Style Intensity Surface (HTML5 Canvas):** Interactive perspective mesh plot with click/touch drag-to-rotate viewing angles and reset view.
   - **Horizontal & Vertical Cross-Sections (Chart.js):** Live $I(x)$ and $I(y)$ intensity graphs with **Through Peak Lobe** vs. **Center Axis** slicing, **Absolute / Normalized** scaling, and optional $\text{TEM}_{00}$ Gaussian reference overlay.

2. **Spatial Mode Selection & Comparison**
   - Supports **$\text{TEM}_{00}$**, **$\text{TEM}_{01}$**, **$\text{TEM}_{10}$**, and **$\text{TEM}_{11}$** modes.
   - **Compare Modes** panel with 4 synchronized mini-heatmaps, an overlay Chart.js comparison graph, and a quantitative mode characteristics table.

3. **Interactive Controls & Calculated Parameters**
   - Controls for **Beam Waist ($w_0$)** ($0.20\text{–}2.50\text{ mm}$), **Laser Intensity ($I_0$)** ($10\text{–}500\text{ W/cm}^2$), **Wavelength ($\lambda$)**, **Calculate**, **Reset**, and **Compare Modes** buttons.
   - Real-time calculation of Peak Intensity ($I_{\text{peak}}$), Normalized Peak Ratio, Total Power ($P$), $1/e^2$ Diameter ($2w_0$), $\text{TEM}_{00}$ FWHM, Rayleigh Range ($z_R$), Far-Field Divergence ($\theta$), and Beam Quality Factors ($M_x^2, M_y^2$).

4. **Academic Documentation Included**
   - Built-in educational sections covering Introduction, Objectives, Working Principle, Mathematical Model, Types of Spatial Modes, Applications, Results & Discussion (with live snapshot), Conclusion, and References.

---

## Mathematical Model Summary

At the beam waist ($z = 0$), the transverse intensity distribution of mode $\text{TEM}_{mn}$ at constant total optical power $P = \frac{1}{2}\pi w_0^2 I_0$ is:

$$I_{mn}(x, y) = \frac{I_0}{2^{m+n}\,m!\,n!} \left[ H_m\!\left(\frac{\sqrt{2}\,x}{w_0}\right) H_n\!\left(\frac{\sqrt{2}\,y}{w_0}\right) \right]^2 \exp\!\left(-\frac{2(x^2 + y^2)}{w_0^2}\right)$$

- **$\text{TEM}_{00}$ ($m=0, n=0$):** $I_{\text{peak}} = I_0$ at $(0, 0)$
- **$\text{TEM}_{10}$ ($m=1, n=0$):** $I_{\text{peak}} = 2e^{-1}I_0 \approx 0.7358\,I_0$ at $(\pm w_0/\sqrt{2}, 0)$
- **$\text{TEM}_{01}$ ($m=0, n=1$):** $I_{\text{peak}} = 2e^{-1}I_0 \approx 0.7358\,I_0$ at $(0, \pm w_0/\sqrt{2})$
- **$\text{TEM}_{11}$ ($m=1, n=1$):** $I_{\text{peak}} = 4e^{-2}I_0 \approx 0.5413\,I_0$ at $(\pm w_0/\sqrt{2}, \pm w_0/\sqrt{2})$

---

## Project Files

```text
laser-beam-profiling-project/
├── index.html          # Dashboard UI, canvases, charts, and 9 educational sections
├── style.css           # Responsive scientific dashboard styling
├── script.js           # Physics engine, 2D/3D canvas renderers, Chart.js logic
├── verify-physics.js   # Node.js verification script for physics & beam equations
└── README.md           # Project documentation
```

---

## How to Run Locally

1. **Direct Browser Launch:**
   Double-click `index.html` (or open it in Chrome, Edge, or Firefox). Note: internet access is used only for the Chart.js CDN script; all physics calculations and 2D/3D Canvas visualizations run locally.
2. **Local Static Server (Recommended):**
   ```bash
   npx serve .
   ```
   or with Python:
   ```bash
   python -m http.server 8000
   ```
   Then open `http://localhost:8000` in your browser.
3. **Run Physics Verification Tests (Optional, requires Node.js):**
   ```bash
   node verify-physics.js
   ```
