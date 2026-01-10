# SFPPy Studio Changelog

All notable changes to SFPPy Studio are documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [0.3.0] - 2026-01-09

### Added
- **Favicon endpoint** - Serves the official SFPPy logo (`/favicon.ico`)
- **Service worker** - Offline capability with intelligent caching strategies
  - CacheFirst for static assets and substance data
  - NetworkFirst for jobs and sessions
  - NetworkOnly for simulations and live searches
- **Substance search deduplication** - Prevents duplicate results from local and PubChem sources
- **Thumbnail caching endpoints** - Local caching for molecule structure images
  - Tiny thumbnails (`/api/substances/thumbnail/{cid}.png`)
  - Large thumbnails with fallback chain
- **SML units display** - Shows "mg/kg" units consistently in all UI locations
- **Search query preservation** - Substance names now respect user's original search term

### Fixed
- **Service worker scope error** - Added `/service-worker.js` endpoint with `Service-Worker-Allowed` header
- **Substance ID format mismatch** - Backend now uses `pubchem_{cid}` format matching frontend
- **404 on assignments endpoint** - Fixed ID format causing substance assignment failures
- **Service worker precache URLs** - Corrected API paths:
  - `/api/assembly/polymers` (was `/api/polymers/list`)
  - `/api/food/simulants` (was `/api/simulants/list`)
  - `/api/substances/common` (was `/api/substances/categories`)
- **CAS array handling** - Frontend now properly handles CAS numbers as arrays
- **Substance name selection** - Search results now prioritize names matching user's query

### Changed
- **Substance detail endpoint** - Returns best matching name based on user query
- **PubChem search** - Improved name scoring to prefer recognizable chemical names

## [0.2.0] - 2025-12-XX

### Added
- **Session management** - In-memory session storage with 24h expiration
- **Substance-to-layer assignments** - C0 matrix with D/k override support
- **Contact steps** - Multi-step temperature/duration profiles
- **Layer configuration** - Polymer selection with thickness units
- **D/k computation endpoints** - Temperature-dependent property calculation
  - Single polymer: `/api/sessions/compute/diffusivity`
  - Batch polymers: `/api/sessions/compute/diffusivity/batch`
  - Temperature range: `/api/sessions/compute/diffusivity/temperature`
- **ToxTree integration** - Cramer classification, TTC values, structural alerts
- **Regulatory data** - EU 10/2011, US FCN, CN GB9685 authorization status
- **Session file format** - `.sfppy.json` for scenario persistence
- **Example sessions** - Pre-built scenarios for testing
- **Concentration unit conversion** - mg/kg, ppm, ppb, g/kg, kg/m³ support
- **Contextual help API** - `/api/help/context/{tab}/{element}`

### Changed
- **Substances module** - Enhanced with migrantToxtree data loading
- **Assembly tab** - Substance assignment matrix with auto-compute toggles

## [0.1.0] - 2025-11-XX

### Added
- **Initial release** - FastAPI-based web application
- **Core simulation interface** - Layer assembly, food/contact conditions
- **PubChem integration** - Substance search and property lookup
- **Polymer database** - Common food contact materials (LDPE, HDPE, PP, PET, etc.)
- **Food simulant database** - EU 10/2011 compliant simulants
- **Geometry calculations** - Cylinder, box, sphere shapes with V/S ratio
- **Migration simulation** - Patankar solver integration
- **Job management** - Background simulation execution
- **Results visualization** - Chart.js migration curves
- **Dark mode** - Theme toggle support
- **Responsive design** - Tailwind CSS styling

---

## Roadmap

### Phase 1 - Examples & Validation (Current)
- [ ] Replay Example 1: Monolayer LDPE with BHT antioxidant
- [ ] Replay Example 2: Bilayer PET/LDPE functional barrier
- [ ] Replay Example 3: Trilayer ABA set-off study
- [ ] Compare Studio results with CLI reference outputs
- [ ] Validate 2×tcontact rule and √t sampling

### Phase 2 - Kinetic Fitting
- [ ] Inverse problem: fit D from experimental CF(t) data
- [ ] Import experimental kinetics (CSV, Excel)
- [ ] Least-squares fitting with uncertainty quantification
- [ ] Confidence intervals and goodness-of-fit metrics

### Phase 3 - RAG + Chat Interface
- [ ] Connect to regulatory knowledge base
- [ ] Natural language parameter input
- [ ] Automated compliance guidance

### Phase 4 - Scenario Management
- [ ] Load/save complete simulation state
- [ ] Scenario library with presets
- [ ] Compare multiple scenarios side-by-side

### Phase 5 - Advanced Features
- [ ] Real-time collaboration (WebSocket)
- [ ] Batch simulation campaigns
- [ ] Monte Carlo uncertainty propagation
- [ ] PDF report generation

---

**Author:** Olivier Vitrac, PhD, HDR
**Organization:** INRAE + Generative Simulation
**Contact:** olivier.vitrac@gmail.com
