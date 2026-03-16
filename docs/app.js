// ======= Responsiveness ==========
const toggleBtn = document.getElementById("controls-toggle");
const controlsPanel = document.getElementById("controls-panel");

toggleBtn.addEventListener("click", () => {
  const isOpen = controlsPanel.classList.toggle("expanded");
  toggleBtn.classList.toggle("open", isOpen);
});

window.addEventListener("resize", () => {
  if (window.innerWidth > 768) {
    controlsPanel.classList.remove("expanded");
    controlsPanel.style.maxHeight = "";
    controlsPanel.style.overflow = "";
    controlsPanel.style.padding = "";
  }
});

// ======= Constants & Setup ==========
const SPECIES = [
  "H2",
  "O2",
  "N2",
  "H2O",
  "CO2",
  "CO",
  "OH",
  "H",
  "O",
  "N",
  "NO",
  "NO2",
  "HO2",
  "H2O2",
  "CH4",
  "Ar",
  "C3H8",
  "NH3",
];
const N = SPECIES.length;
const C_ELEM = 5;
const MOL_WEIGHTS = [
  2.01588, 31.9988, 28.0134, 18.01528, 44.0095, 28.0101, 17.00734, 1.00794,
  15.9994, 14.0067, 30.0061, 46.0055, 33.00674, 34.01468, 16.04246, 39.948,
  44.09562, 17.03052,
];

const COLORS = [
  "#3b82f6",
  "#ef4444",
  "#64748b",
  "#06b6d4",
  "#f97316",
  "#8b5cf6",
  "#22c55e",
  "#ec4899",
  "#eab308",
  "#14b8a6",
  "#f43f5e",
  "#a855f7",
  "#84cc16",
  "#6366f1",
  "#f59e0b",
  "#78716c",
  "#0ea5e9",
  "#d946ef",
];

const DISPLAY = {
  H2: "H\u2082",
  O2: "O\u2082",
  N2: "N\u2082",
  H2O: "H\u2082O",
  CO2: "CO\u2082",
  CO: "CO",
  OH: "OH",
  H: "H",
  O: "O",
  N: "N",
  NO: "NO",
  NO2: "NO\u2082",
  HO2: "HO\u2082",
  H2O2: "H\u2082O\u2082",
  CH4: "CH\u2084",
  Ar: "Ar",
  C3H8: "C\u2083H\u2088",
  NH3: "NH\u2083",
};

const IDX = {};
SPECIES.forEach((s, i) => (IDX[s] = i));

const ATOM_MATRIX = [
  [0, 2, 0, 0, 0],
  [0, 0, 2, 0, 0],
  [0, 0, 0, 2, 0],
  [0, 2, 1, 0, 0],
  [1, 0, 2, 0, 0],
  [1, 0, 1, 0, 0],
  [0, 1, 1, 0, 0],
  [0, 1, 0, 0, 0],
  [0, 0, 1, 0, 0],
  [0, 0, 0, 1, 0],
  [0, 0, 1, 1, 0],
  [0, 0, 2, 1, 0],
  [0, 1, 2, 0, 0],
  [0, 2, 2, 0, 0],
  [1, 4, 0, 0, 0],
  [0, 0, 0, 0, 1],
  [3, 8, 0, 0, 0],
  [0, 3, 0, 1, 0],
];

const STOICH_O2 = { CH4: 2, H2: 0.5, C3H8: 5 };

const PRESETS = {
  ch4_air: { fuel: "CH4", moles: 1, air: true },
  h2_air: { fuel: "H2", moles: 1, air: true },
  c3h8_air: { fuel: "C3H8", moles: 1, air: true },
  h2_o2: { fuel: "H2", moles: 2, air: false },
  ch4_o2: { fuel: "CH4", moles: 1, air: false },
};

let wasmReady = false;
let state = null;
let solverMode = "fixed-t";

// ======= Event Listeners ==========
document.querySelectorAll(".mode-toggle button").forEach((btn) => {
  btn.addEventListener("click", () => {
    document
      .querySelectorAll(".mode-toggle button")
      .forEach((b) => b.classList.remove("active"));
    btn.classList.add("active");
    solverMode = btn.dataset.mode;

    const isAdiabatic = solverMode === "adiabatic";
    document
      .getElementById("temp-group")
      .classList.toggle("hidden", isAdiabatic);
    document
      .getElementById("tfeed-group")
      .classList.toggle("hidden", !isAdiabatic);
    document
      .getElementById("summary-cards")
      .classList.toggle("hidden", !isAdiabatic);
    runSolver(true);
  });
});

document.getElementById("temp").addEventListener("input", (e) => {
  document.getElementById("temp-val").textContent = e.target.value + " K";
  runSolver(false);
});

document.getElementById("tfeed").addEventListener("input", (e) => {
  document.getElementById("tfeed-val").textContent = e.target.value + " K";
  runSolver(true);
});

document.getElementById("pressure").addEventListener("input", (e) => {
  document.getElementById("pressure-val").textContent =
    parseFloat(e.target.value).toFixed(1) + " bar";
  runSolver(true);
});

document.getElementById("phi").addEventListener("input", (e) => {
  document.getElementById("phi-val").textContent =
    "\u03C6 = " + parseFloat(e.target.value).toFixed(2);
  runSolver(true);
});

document.getElementById("preset").addEventListener("change", (e) => {
  document
    .getElementById("custom-feed")
    .classList.toggle("visible", e.target.value === "custom");
  runSolver(true);
});

document.querySelectorAll("#custom-feed input").forEach((inp) => {
  inp.addEventListener("input", () => runSolver(true));
});

// ======= State & Helpers ==========
function resetState(feed, totalAtoms) {
  const n = new Float64Array(N);
  let nT = 0;

  for (let i = 0; i < N; i++) {
    let possible = true;
    for (let j = 0; j < C_ELEM; j++) {
      if (ATOM_MATRIX[i][j] > 0 && totalAtoms[j] <= 1e-11) {
        possible = false;
        break;
      }
    }
    if (possible) {
      n[i] = feed[i] + 0.1;
    } else {
      n[i] = 1e-12;
    }
    nT += n[i];
  }

  state = {
    n: n,
    nT: nT,
    lambda: new Float64Array(C_ELEM).fill(1.0),
    nu: 1.0,
  };
}

function buildFeed(presetId, phi) {
  const feed = new Float64Array(N);

  if (presetId === "custom") {
    document.querySelectorAll("#custom-feed input").forEach((inp) => {
      const sp = inp.dataset.sp;
      if (sp && IDX[sp] !== undefined) {
        feed[IDX[sp]] = parseFloat(inp.value) || 0;
      }
    });
    return feed;
  }

  const preset = PRESETS[presetId];
  const stoichO2 = STOICH_O2[preset.fuel] * preset.moles;
  const actualO2 = stoichO2 / phi;

  feed[IDX[preset.fuel]] = preset.moles;
  feed[IDX.O2] = actualO2;

  if (preset.air) {
    feed[IDX.N2] = actualO2 * 3.76;
    feed[IDX.Ar] = actualO2 * 0.045;
  }

  return feed;
}

function computeTotalAtoms(feed) {
  const atoms = new Float64Array(C_ELEM);
  for (let i = 0; i < N; i++) {
    for (let j = 0; j < C_ELEM; j++) {
      atoms[j] += feed[i] * ATOM_MATRIX[i][j];
    }
  }
  return atoms;
}

// ======= Core Solver Dispatch ==========
function runSolver(doReset = false) {
  if (!wasmReady) return;
  if (solverMode === "adiabatic") {
    runAdiabatic();
  } else {
    runFixedT(doReset);
  }
}

function runFixedT(doReset) {
  const T = parseFloat(document.getElementById("temp").value);
  const P = parseFloat(document.getElementById("pressure").value);
  const P_ref = 1.0;
  const phi = parseFloat(document.getElementById("phi").value);
  const feed = buildFeed(document.getElementById("preset").value, phi);

  const totalFeed = feed.reduce((a, b) => a + b, 0);
  if (totalFeed < 1e-10) return;

  const totalAtoms = computeTotalAtoms(feed);
  for (let j = 0; j < C_ELEM; j++) {
    if (totalAtoms[j] < 1e-12) totalAtoms[j] = 1e-12;
  }

  if (doReset || !state) {
    resetState(feed, totalAtoms);
  }

  for (let j = 0; j < C_ELEM; j++) totalAtoms[j] /= totalFeed;

  const dblSize = 8;
  const nPtr = Module._malloc(N * dblSize);
  const nTPtr = Module._malloc(dblSize);
  const lamPtr = Module._malloc(C_ELEM * dblSize);
  const nuPtr = Module._malloc(dblSize);
  const atomPtr = Module._malloc(C_ELEM * dblSize);

  const HEAPF64 = new Float64Array(
    Module.wasmMemory ? Module.wasmMemory.buffer : HEAP8.buffer,
  );

  HEAPF64.set(state.n, nPtr / dblSize);
  HEAPF64[nTPtr / dblSize] = state.nT;
  HEAPF64.set(state.lambda, lamPtr / dblSize);
  HEAPF64[nuPtr / dblSize] = state.nu;
  HEAPF64.set(totalAtoms, atomPtr / dblSize);

  const result = Module._solve(
    T,
    P,
    P_ref,
    nPtr,
    nTPtr,
    lamPtr,
    nuPtr,
    atomPtr,
    N,
    C_ELEM,
  );

  state.n = new Float64Array(HEAPF64.buffer, nPtr, N).slice();
  state.nT = HEAPF64[nTPtr / dblSize];
  state.lambda = new Float64Array(HEAPF64.buffer, lamPtr, C_ELEM).slice();
  state.nu = HEAPF64[nuPtr / dblSize];

  Module._free(nPtr);
  Module._free(nTPtr);
  Module._free(lamPtr);
  Module._free(nuPtr);
  Module._free(atomPtr);
  updateUI(result, T, null, null, null);
}

function runAdiabatic() {
  const T_feed = parseFloat(document.getElementById("tfeed").value);
  const P = parseFloat(document.getElementById("pressure").value);
  const P_ref = 1.0;
  const phi = parseFloat(document.getElementById("phi").value);
  const feed = buildFeed(document.getElementById("preset").value, phi);

  const totalFeed = feed.reduce((a, b) => a + b, 0);
  if (totalFeed < 1e-10) return;

  const totalAtoms = computeTotalAtoms(feed);
  for (let j = 0; j < C_ELEM; j++) {
    if (totalAtoms[j] < 1e-12) totalAtoms[j] = 1e-12;
  }

  const feedNorm = new Float64Array(N);
  for (let i = 0; i < N; i++) feedNorm[i] = feed[i] / totalFeed;
  for (let j = 0; j < C_ELEM; j++) totalAtoms[j] /= totalFeed;

  const dblSize = 8;
  const intSize = 4;
  const feedPtr = Module._malloc(N * dblSize);
  const nPtr = Module._malloc(N * dblSize);
  const nTPtr = Module._malloc(dblSize);
  const gammaPtr = Module._malloc(dblSize);
  const atomPtr = Module._malloc(C_ELEM * dblSize);
  const numIterPtr = Module._malloc(intSize);

  const wasmBuf = Module.wasmMemory ? Module.wasmMemory.buffer : HEAP8.buffer;
  const HEAPF64 = new Float64Array(wasmBuf);
  const HEAP32 = new Int32Array(wasmBuf);

  HEAPF64.set(feedNorm, feedPtr / dblSize);
  HEAPF64.set(totalAtoms, atomPtr / dblSize);

  const T_guess = 2200.0;

  const T_ad = Module._solve_adiabatic_temperature(
    T_feed,
    T_guess,
    P,
    P_ref,
    feedPtr,
    nPtr,
    nTPtr,
    gammaPtr,
    atomPtr,
    N,
    C_ELEM,
    numIterPtr,
  );

  const n_out = new Float64Array(HEAPF64.buffer, nPtr, N).slice();
  const nT_out = HEAPF64[nTPtr / dblSize];
  const gamma_out = HEAPF64[gammaPtr / dblSize];
  const iter_out = HEAP32[numIterPtr / intSize];

  Module._free(feedPtr);
  Module._free(nPtr);
  Module._free(nTPtr);
  Module._free(gammaPtr);
  Module._free(atomPtr);
  Module._free(numIterPtr);

  state = {
    n: n_out,
    nT: nT_out,
    lambda: new Float64Array(C_ELEM).fill(0),
    nu: 0,
  };

  updateUI(T_ad > 0 ? 0 : -1, T_ad, T_ad, gamma_out, iter_out);
}

// ======= UI ==========
function updateUI(result, T, T_ad, gamma, iter) {
  const statusEl = document.getElementById("status");

  if (result >= 0) {
    statusEl.className = "status ok";
    if (solverMode === "adiabatic" && T_ad > 0) {
      statusEl.textContent = `Converged in ${iter} outer iters \u2014 T_ad = ${T_ad.toFixed(1)} K`;
    } else {
      statusEl.textContent = `Converged in ${result} iters @ T = ${T} K`;
    }
  } else if (result === -1) {
    statusEl.className = "status err";
    statusEl.textContent = "Singular Jacobian";
    return;
  } else {
    statusEl.className = "status err";
    statusEl.textContent = "Did not converge";
    return;
  }

  if (solverMode === "adiabatic" && T_ad > 0) {
    document.getElementById("card-tad").textContent = T_ad.toFixed(1);
    document.getElementById("card-gamma").textContent =
      gamma > 0 ? gamma.toFixed(4) : "\u2014";
  }

  const nT = state.nT;
  let totalMass = 0;
  for (let i = 0; i < N; i++) totalMass += state.n[i] * MOL_WEIGHTS[i];

  const rows = [];
  for (let i = 0; i < N; i++) {
    const xi = state.n[i] / nT;
    const wi = (state.n[i] * MOL_WEIGHTS[i]) / totalMass;
    rows.push({ idx: i, name: SPECIES[i], ni: state.n[i], xi, wi });
  }
  rows.sort((a, b) => b.xi - a.xi);

  const chart = document.getElementById("chart");
  const maxX = Math.max(...rows.map((r) => r.xi), 1e-10);

  chart.innerHTML = rows
    .filter((r) => r.xi > 1e-8)
    .map((r) => {
      const pct = ((r.xi / maxX) * 100).toFixed(2);
      const color = COLORS[r.idx];
      const xStr =
        r.xi > 0.001 ? (r.xi * 100).toFixed(2) + "%" : r.xi.toExponential(2);
      return `<div class="bar-row">
  <span class="bar-label">${DISPLAY[r.name]}</span>
  <div class="bar-track">
    <div class="bar-fill" style="width:${pct}%;background:${color}"></div>
  </div>
  <span class="bar-value">${xStr}</span>
</div>`;
    })
    .join("");

  const tbody = document.querySelector("#results tbody");
  tbody.innerHTML = rows
    .filter((r) => r.ni > 1e-10)
    .map((r) => {
      const niStr = r.ni > 0.001 ? r.ni.toFixed(6) : r.ni.toExponential(3);
      const xiStr =
        r.xi > 0.001 ? (r.xi * 100).toFixed(4) + "%" : r.xi.toExponential(3);
      const wiStr =
        r.wi > 0.001 ? (r.wi * 100).toFixed(4) + "%" : r.wi.toExponential(3);
      return `<tr>
  <td>${DISPLAY[r.name]}</td>
  <td>${niStr}</td>
  <td>${xiStr}</td>
  <td>${wiStr}</td>
</tr>`;
    })
    .join("");
}

// ======= Global WASM Hook ==========
var Module = {
  onRuntimeInitialized: function() {
    Module._init();
    wasmReady = true;
    document.getElementById("status").className = "status ok";
    document.getElementById("status").textContent = "WASM ready";
    runSolver(true);
  },
};
