// Phase 13 Part 3: shared frontend helpers. No build step, no framework (see
// docs/ARCHITECTURE.md) -- plain fetch() against the FastAPI JSON API, plain DOM updates.

const LWMOM = (() => {
  async function postJson(url, body) {
    const resp = await fetch(url, {
      method: "POST",
      headers: { "Content-Type": "application/json" },
      body: JSON.stringify(body),
    });
    const data = await resp.json().catch(() => ({}));
    if (!resp.ok) {
      const err = new Error(data.detail?.message || data.detail || resp.statusText);
      err.body = data;
      throw err;
    }
    return data;
  }

  async function postForm(url, formData) {
    const resp = await fetch(url, { method: "POST", body: formData });
    const data = await resp.json().catch(() => ({}));
    if (!resp.ok) {
      const err = new Error(data.detail?.message || data.detail || resp.statusText);
      err.body = data;
      throw err;
    }
    return data;
  }

  // Polls GET /api/jobs/{id} every `intervalMs` until status leaves queued/running.
  // Calls onTick(job) on every poll, then resolves with the final job.
  function pollJob(jobId, { intervalMs = 500, onTick = () => {} } = {}) {
    return new Promise((resolve, reject) => {
      const tick = async () => {
        try {
          const resp = await fetch(`/api/jobs/${jobId}`);
          const job = await resp.json();
          onTick(job);
          if (job.status === "done" || job.status === "failed" || job.status === "cancelled") {
            resolve(job);
          } else {
            setTimeout(tick, intervalMs);
          }
        } catch (err) {
          reject(err);
        }
      };
      tick();
    });
  }

  // PHASE13C section 3 ("Words"): validation-label CODES stay in the raw job result (and in the
  // `title` tooltip below, for anyone who wants them) but the screen shows a plain-English note --
  // never a snake_case/internal identifier as the only visible text.
  const LABEL_DESCRIPTIONS = {
    LOCAL_PERIODICITY_APPROXIMATION: () => "Estimated from individual cells — confirm in HFSS before fabrication.",
    FEED_REGION_MASKED: () => "Feed point excluded.",
    MESH_NOT_CONVERGED: () => "Fast accuracy — switch to Accurate in Advanced for final designs.",
    MESH_CONVERGED: () => "Accurate (mesh-converged).",
    ZERO_THICKNESS_METAL_MODEL: () => "Modeled with zero-thickness metal (standard for this tool).",
    REFERENCE_AMPLITUDE_ASSUMED_UNIFORM: () => "No magnitude file — uniform amplitude assumed.",
    EXPERIMENTAL: () => "Experimental.",
    CELLS_PER_FRINGE: (v) => `Cells are large for this wave (${v} per fringe). Try a smaller cell size for a cleaner pattern.`,
    PHASE_COVERAGE_SHORTFALL: (v) => `Some cells can't reach the needed value (worst error ${v}).`,
    REFERENCE_WAVE_RESCALED_TO_BETA_BAR: (v) => `Reference wave adjusted to the board's own wavenumber (${v}).`,
    REFERENCE_WAVE_NOT_RESCALED: (v) => `Using the raw reference wave — possible error: ${v}.`,
    REFERENCE_FIELD_FREQUENCY_MISMATCH: () => "The uploaded file's frequency doesn't match this design — check Step 1.",
    STEEP_GRADIENT: () => "Phase changes sharply from cell to cell — pattern accuracy may suffer.",
    HFSS_RESIDUAL_AT_RESONANCE: (v) => `Small unexplained gap vs HFSS near resonance (${v}).`,
    DISPERSION_FROM_HFSS: () => "Dispersion data imported from HFSS.",
    GEOMETRY_INVERSION_LOCKED: () => "Geometry inversion isn't available in this build.",
    STEP_EXPORT_LOCKED: () => "STEP export isn't available in this build.",
    HFSS_VALIDATED_OFF_RESONANCE: () => "Validated against HFSS off resonance.",
  };

  // PHASE13D item 2 ("Words"): label params are formatted in Python (labels.py) with a bare
  // "deg" suffix (e.g. "162deg") -- the raw code stays exactly that (title tooltip, matching
  // labels.py 1:1), but the visible text gets the degree sign like every other angle in the UI.
  function _formatLabelParam(param) {
    return param == null ? param : param.replace(/deg/g, "°");
  }

  function describeLabel(raw) {
    const m = /^([A-Z0-9_]+)(?:\((.*)\))?$/.exec(raw || "");
    if (!m) return raw;
    const [, code, param] = m;
    const fn = LABEL_DESCRIPTIONS[code];
    return fn ? fn(_formatLabelParam(param)) : raw;
  }

  function renderLabels(container, labels, warnings) {
    container.innerHTML = "";
    (labels || []).forEach((l) => {
      const span = document.createElement("span");
      span.className = "label";
      span.textContent = describeLabel(l);
      span.title = l;
      container.appendChild(span);
    });
    (warnings || []).forEach((w) => {
      const span = document.createElement("span");
      span.className = "label warn";
      span.textContent = w;
      container.appendChild(span);
    });
  }

  function markActiveNav() {
    const path = window.location.pathname;
    document.querySelectorAll("nav a").forEach((a) => {
      if (a.getAttribute("href") === path) a.classList.add("active");
    });
  }

  // Phase 17 hygiene: a visible version/commit footer so screenshots and bug reports are traceable.
  async function renderVersionFooter() {
    let info = {};
    try {
      const resp = await fetch("/api/version");
      info = await resp.json();
    } catch {
      // offline/startup race -- footer just stays blank, never blocks the page.
    }
    const footer = document.createElement("footer");
    footer.className = "version-footer";
    const commitPart = info.commit ? ` · ${info.commit}` : "";
    footer.textContent = `HMSA Designer · v${info.version || "unknown"}${commitPart}`;
    document.body.appendChild(footer);
  }

  return { postJson, postForm, pollJob, renderLabels, describeLabel, markActiveNav, renderVersionFooter };
})();

document.addEventListener("DOMContentLoaded", LWMOM.markActiveNav);
document.addEventListener("DOMContentLoaded", LWMOM.renderVersionFooter);
