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

  function renderLabels(container, labels, warnings) {
    container.innerHTML = "";
    (labels || []).forEach((l) => {
      const span = document.createElement("span");
      span.className = "label";
      span.textContent = l;
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

  return { postJson, postForm, pollJob, renderLabels, markActiveNav };
})();

document.addEventListener("DOMContentLoaded", LWMOM.markActiveNav);
