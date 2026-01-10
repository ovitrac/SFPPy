"""
SFPPy Studio - Main FastAPI Application

Comprehensive web application for food contact migration analysis.

Author: Olivier Vitrac, PhD, HDR
Organization: INRAE + Generative Simulation

ROADMAP / TODO
==============

Phase 1 - Examples & Validation (Current)
------------------------------------------
- [ ] Replay Example 1: Monolayer LDPE with BHT antioxidant
- [ ] Replay Example 2: Bilayer PET/LDPE functional barrier
- [ ] Replay Example 3: Trilayer ABA set-off study
- [ ] Compare Studio results with CLI reference outputs
- [ ] Validate 2*tcontact rule and sqrt time sampling

Phase 2 - Kinetic Fitting
-------------------------
- [ ] Inverse problem: fit D from experimental CF(t) data
- [ ] Import experimental kinetics (CSV, Excel)
- [ ] Least-squares fitting with uncertainty quantification
- [ ] Confidence intervals on fitted parameters
- [ ] Residual analysis and goodness-of-fit metrics
- [ ] Export fitted parameters for compliance reports

Phase 3 - RAG + Chat Interface
-------------------------------
- [ ] Connect to regulatory knowledge base (EU 10/2011, FDA, GB 9685)
- [ ] Natural language queries for SML lookup
- [ ] Substance identification from partial names/synonyms
- [ ] Automated compliance guidance based on simulation results
- [ ] Citation and traceability for regulatory decisions
- [ ] Integration with RAGGAE or similar RAG pipeline
- [ ] **Chat interface to auto-feed simulation parameters**
- [ ] Parse natural language into layer/substance/conditions config
- [ ] "Simulate BHT migration from 100µm LDPE at 40°C for 10 days"
- [ ] Conversational refinement of parameters

Phase 4 - Scenario Replay & State Management
--------------------------------------------
- [ ] **Refeed interface from saved scenarios**
- [ ] Load/save complete simulation state (JSON export/import)
- [ ] Replay examples 1-5 with one-click restore
- [ ] Scenario library with user-defined presets
- [ ] Compare multiple scenarios side-by-side
- [ ] Undo/redo for parameter changes
- [ ] Session persistence across browser refresh

Phase 5 - Design Review (Required)
----------------------------------
- [ ] Review current state management architecture
- [ ] Evaluate state.* structure for extensibility
- [ ] API design for state serialization/deserialization
- [ ] UI/UX review for parameter input flow
- [ ] Ensure bidirectional sync: UI ↔ State ↔ Backend
- [ ] Document data flow for chat → parameters → simulation

Phase 6 - Advanced Features
---------------------------
- [ ] Real-time collaboration (WebSocket)
- [ ] Batch simulation campaigns
- [ ] Monte Carlo uncertainty propagation
- [ ] PDF report generation with compliance statements
- [ ] API for external integrations
"""

import os
import sys
from pathlib import Path
from datetime import datetime

from fastapi import FastAPI, Request
from fastapi.responses import HTMLResponse, JSONResponse
from fastapi.staticfiles import StaticFiles
from fastapi.templating import Jinja2Templates

# Add parent paths for imports
sys.path.insert(0, str(Path(__file__).parent.parent.parent))

from studio.version import __version__, get_version_info

# Import routes
from studio.app.routes import assembly, food, substances, simulation, jobs, config, exports, sessions, fitting

# Create FastAPI app
app = FastAPI(
    title="SFPPy Studio",
    description="Comprehensive Food Contact Migration Analysis",
    version=__version__,
    docs_url="/docs",
    redoc_url="/redoc",
)

# Static files and templates
STUDIO_DIR = Path(__file__).parent.parent
APP_DIR = Path(__file__).parent


# Favicon - serve the SFPPy logo
@app.get("/favicon.ico")
async def favicon():
    """Serve the SFPPy SVG logo as favicon."""
    from fastapi.responses import FileResponse
    logo_path = STUDIO_DIR.parent / "docs" / "assets" / "SFPPy.svg"
    if logo_path.exists():
        return FileResponse(path=str(logo_path), media_type="image/svg+xml")
    # Fallback to simple green "S" favicon if logo not found
    from fastapi.responses import Response
    svg = '<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 32 32"><circle cx="16" cy="16" r="14" fill="#10B981"/><text x="16" y="21" text-anchor="middle" font-size="14" font-weight="bold" fill="white">S</text></svg>'
    return Response(content=svg, media_type="image/svg+xml")


# Service Worker with proper scope header (must be before static files mount)
@app.get("/service-worker.js")
async def service_worker():
    """Serve service worker from root with proper scope header."""
    from fastapi.responses import FileResponse
    sw_path = APP_DIR / "static" / "js" / "service-worker.js"
    return FileResponse(
        path=str(sw_path),
        media_type="application/javascript",
        headers={"Service-Worker-Allowed": "/"}
    )


app.mount("/static", StaticFiles(directory=APP_DIR / "static"), name="static")
templates = Jinja2Templates(directory=APP_DIR / "templates")

# Include API routes
app.include_router(sessions.router, prefix="/api/sessions", tags=["Sessions"])
app.include_router(assembly.router, prefix="/api/assembly", tags=["Assembly"])
app.include_router(food.router, prefix="/api/food", tags=["Food & Conditions"])
app.include_router(substances.router, prefix="/api/substances", tags=["Substances"])
app.include_router(simulation.router, prefix="/api/simulation", tags=["Simulation"])
app.include_router(jobs.router, prefix="/api/jobs", tags=["Jobs"])
app.include_router(config.router, prefix="/api/config", tags=["Configuration"])
app.include_router(exports.router, prefix="/api/export", tags=["Export"])
app.include_router(fitting.router, prefix="/api/fitting", tags=["Fitting"])


# ========== WEB PAGES ==========

@app.get("/", response_class=HTMLResponse)
async def index(request: Request):
    """Main application page."""
    version_info = get_version_info()
    return templates.TemplateResponse("index.html", {
        "request": request,
        "version": version_info,
        "current_year": datetime.now().year,
    })


# ========== API INFO ==========

@app.get("/api/info")
async def api_info():
    """Get application information."""
    return JSONResponse({
        "success": True,
        "app": "SFPPy Studio",
        "description": "Food Contact Migration Analysis",
        **get_version_info(),
        "timestamp": datetime.utcnow().isoformat() + "Z",
    })


@app.get("/api/health")
async def health_check():
    """Health check endpoint."""
    return JSONResponse({
        "status": "healthy",
        "timestamp": datetime.utcnow().isoformat() + "Z",
    })


# ========== HELP API ==========

@app.get("/api/help/context/{tab}/{element}")
async def get_help(tab: str, element: str):
    """Get contextual help for UI elements."""
    # Help content database (simplified - could be loaded from YAML)
    help_db = {
        "assembly": {
            "layer_count": {
                "title": "Number of Layers",
                "short": "Set the number of layers in your packaging structure.",
                "long": "Layer 1 is always in contact with food. Add more layers to model multilayer packaging. Up to 10 layers are supported.",
                "tips": ["Start with 1 layer for simple scenarios", "Use 2-3 layers for functional barrier studies"],
            },
            "polymer": {
                "title": "Polymer Type",
                "short": "Select the polymer material for this layer.",
                "long": "The polymer type determines default diffusion and partition properties. Common polymers include LDPE, HDPE, PP, PET, and PS.",
                "tips": ["PET acts as a good functional barrier", "LDPE has high diffusivity"],
            },
            "thickness": {
                "title": "Layer Thickness",
                "short": "Set the thickness of this layer.",
                "long": "Thickness affects migration kinetics. Thicker layers take longer to reach equilibrium. Common units: µm for films, mm for sheets.",
                "tips": ["Films: 10-100 µm", "Sheets: 0.5-3 mm", "Containers: 1-5 mm"],
            },
            "D_override": {
                "title": "Diffusivity Override",
                "short": "Override the computed diffusion coefficient.",
                "long": "The diffusion coefficient D (m²/s) describes how fast molecules move through the polymer. Use this when you have experimental values or the model doesn't fit your specific case.",
                "tips": ["Typical range: 1e-16 to 1e-12 m²/s", "Use scientific notation: 1.2e-14"],
                "references": ["Piringer et al., Food Additives & Contaminants, 2008"],
            },
            "k_override": {
                "title": "Partition Coefficient Override",
                "short": "Override the computed partition coefficient.",
                "long": "The partition coefficient k describes the equilibrium distribution between polymer and food. k = C_polymer / C_food. Use experimental values when available.",
                "tips": ["k > 1: substance prefers polymer", "k < 1: substance prefers food"],
            },
            "substances": {
                "title": "Substance Assignments",
                "short": "Assign substances to layers with initial concentrations.",
                "long": "For each substance and layer combination, specify: (1) whether the substance is present (checkbox), (2) initial concentration C0 (default 1000 mg/kg), (3) units. D (diffusivity) and k (partition) are auto-computed based on polymer and temperature; override values if needed.",
                "tips": [
                    "Check box to assign substance to layer",
                    "Default C0 = 1000 mg/kg when assigned",
                    "Orange border = override value active",
                    "k0 (food partition) is set in Food tab",
                ],
            },
        },
        "food": {
            "category": {
                "title": "Food Category",
                "short": "Select the type of food in contact.",
                "long": "Food category determines which simulant to use for migration testing (EU 10/2011). Categories include aqueous, acidic, alcoholic, fatty, and dry foods.",
                "tips": ["Fatty foods: use olive oil or ethanol 95%", "Aqueous foods: use water or acetic acid 3%"],
            },
            "geometry": {
                "title": "Packaging Geometry",
                "short": "Define the 3D shape of the packaging.",
                "long": "The geometry determines the contact surface area and food volume. Volume/Surface ratio affects migration kinetics.",
                "tips": ["Enter dimensions in mm", "Volume and surface are calculated automatically"],
            },
            "contact_time": {
                "title": "Contact Duration",
                "short": "Set how long food is in contact with packaging.",
                "long": "Contact time is crucial for migration assessment. EU standard test is 10 days at specified temperature. Use actual shelf life for worst-case scenarios.",
                "tips": ["EU test: 10 days", "Shelf life: weeks to months"],
            },
            "k0": {
                "title": "Food Partition Coefficient (k0)",
                "short": "Partition coefficient between polymer and food.",
                "long": "k0 describes the equilibrium distribution of a substance between polymer and food simulant. K_F/P = k/k0 is the effective food/polymer partition coefficient. Values depend on food type (fatty vs. aqueous) and substance polarity (logP).",
                "tips": [
                    "k0 > 1: substance prefers polymer (less migration)",
                    "k0 < 1: substance prefers food (more migration)",
                    "Fatty foods: lipophilic substances migrate more",
                    "Aqueous foods: hydrophilic substances migrate more",
                ],
            },
            "simulant": {
                "title": "Food Simulant",
                "short": "EU 10/2011 standard simulant for migration testing.",
                "long": "Food simulants are standardized test media defined by EU 10/2011. They mimic different food types for migration testing. Selection affects k0 calculation and migration kinetics.",
                "tips": [
                    "Fatty foods: use ethanol 95% or olive oil",
                    "Aqueous foods: use water or acetic acid 3%",
                    "Alcoholic beverages: use ethanol 10-50%",
                    "Dry foods: use Tenax",
                ],
            },
        },
        "substances": {
            "search": {
                "title": "Substance Search",
                "short": "Search for chemical substances in databases.",
                "long": "Search by common name, CAS number, or PubChem CID. Results are cached locally for fast access. Common food contact additives are pre-cached.",
                "tips": ["Try partial names: 'Irganox'", "Use CAS for exact match: '2082-79-3'"],
            },
            "C0": {
                "title": "Initial Concentration",
                "short": "Set the starting concentration in the layer.",
                "long": "C0 (mg/kg polymer) is the amount of substance present in the layer before any migration occurs. This is typically from manufacturing.",
                "tips": ["Typical range: 100-10000 mg/kg", "Check supplier data sheets"],
            },
            "SML": {
                "title": "Specific Migration Limit",
                "short": "EU regulatory limit for this substance.",
                "long": "SML is the maximum allowed concentration in food (mg/kg food). Defined in EU 10/2011 for authorized substances. Exceeding SML means non-compliance.",
                "tips": ["Check: CF(final) < SML for compliance"],
            },
        },
    }

    if tab in help_db and element in help_db[tab]:
        return JSONResponse({
            "success": True,
            "element": element,
            **help_db[tab][element],
        })
    else:
        return JSONResponse({
            "success": False,
            "error": f"No help available for {tab}/{element}",
        }, status_code=404)
