# SFPPy — TODO / Roadmap

*Maintained by Olivier Vitrac, PhD, HDR — kept identical between the development
and public trees. Items are scientific/engineering tasks, ordered by priority.*

---

## High priority

- [ ] **Expose the real-units physical engine through a public survey entry point.**
      `survey/chain_engine.py` (N-layer × K-step real-units solver) and
      `survey/physical_query.py` (physical-primary query) are shipped but currently
      reachable only via the two-step / multilayer survey path. Add a documented,
      flag-gated public API and a worked example.
- [ ] **Publish the bilayer / two-step survey classes.** The public `Survey` covers
      the monolayer single-step path; the multilayer and chained-contact survey
      drivers remain to be packaged for general use.
- [ ] **Tests for the new engine modules.** Add unit tests for `chain_engine`
      (numerical identity vs. the bilayer two-step reference), `physical_query`
      (key parity with the Fo-surface query), and `survey.utils.cf_at_user_grid`.

## Medium priority

- [ ] **Regenerate the API HTML documentation on every release** (`docs/`), so the
      published reference always matches the shipped code.
- [ ] **Consolidate the survey cache schema** and document its on-disk layout
      (`.survey_cache/`), including versioned keys.
- [ ] **Survey module versioning policy.** `survey/__init__.py` and
      `survey/CHANGELOG.md` are now aligned at `0.4.0`; keep them in lock-step on
      every survey-facing change.

## Low priority / housekeeping

- [ ] Reconcile the documentation sources (knowledge-base indexing) between the
      development and public trees.
- [ ] Audit the bundled PubChem / ToxTree caches for staleness and provide a
      refresh utility.
- [ ] Type-checking (mypy) and style (black/isort) sweep over the `survey/` package.

---

## Recently completed

- [x] **v1.6 — physical migration engine.** Constant sparse Jacobian passed to the
      BDF integrator (cold-start speed-up); `C0eq` guardrail and unit-`C0eq` resume
      fix across the physical chain; canonical equilibrium helpers consolidated in
      `survey/workers.py`.
