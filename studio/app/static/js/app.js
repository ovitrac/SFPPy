/**
 * SFPPy Studio - Main Application JavaScript
 *
 * Handles UI interactions, API calls, and state management.
 */

// ========== STATE ==========
const state = {
    sessionId: null,  // Session ID for backend state management
    layers: [{ index: 1, polymer: 'LDPE', thickness: 100, thickness_unit: 'um' }],
    steps: [{ index: 1, temperature_C: 40, duration: 10, duration_unit: 'days', with_food: true, type: 'storage', setoff_type: 'stacked' }],
    food: { category: 'fatty', simulant: 'ethanol', matrixType: 'liquid' },
    geometry: { shape: 'cylinder', dimensions: { radius: 50, height: 100 } },
    substances: [],
    // C0 matrix: { substanceId: { layerIndex: { value, unit }, ... }, ... }
    c0Matrix: {},
    // C0 unit preference per substance: { substanceId: 'mg/kg' }
    c0Units: {},
    // CF0: Initial concentration in food per substance: { substanceId: { value, unit } }
    CF0: {},
    currentJobId: null,
    migrationChart: null,
    profileChart: null,  // Concentration profile chart
    selectedLayerIndex: 1,  // Currently selected layer for editing
    materialsCache: null,   // Cache for discovered materials
    previewedSubstance: null, // Currently previewed substance for adding
    concentrationUnits: ['mg/kg', 'ppm', 'ppb', 'g/kg', 'kg/kg', 'ng/kg'],  // Available units
    defaultC0Unit: 'mg/kg',
    // Restart functionality
    restartEnabled: false,
    restartSource: null,  // { filename, session } for restart from previous simulation
    // Session tracking
    loadedSessionName: null,  // Name of loaded session for display
};

// Material Tg database (from layer.py)
const MATERIAL_TG = {
    // Rubbery/Low Tg (-130 to 0°C)
    'LDPE': -130, 'LLDPE': -120, 'HDPE': -100,
    'PP': -10, 'PPrubber': -20, 'oPP': -10,
    // Mid-range Tg (0 to 70°C)
    'PA6': 50, 'PA66': 70, 'PVAc': 30, 'wPET': 46,
    // Glassy/High Tg (>70°C) - Barrier
    'gPET': 76, 'PET': 76, 'rPET': 76, 'PBT': 66, 'PEN': 120,
    'PS': 100, 'HIPS': 95, 'PMMA': 105, 'SBS': 90,
    'rigidPVC': 87, 'plasticizedPVC': 20,
    // Adhesives
    'AdhesiveEVA': -30, 'AdhesiveVAE': 10, 'AdhesivePVAC': 35,
    'AdhesiveAcrylate': -20, 'AdhesivePU': -50,
    'AdhesiveNaturalRubber': -70, 'AdhesiveSyntheticRubber': -50,
    // Paper & Board
    'Paper': 200, 'Cardboard': 200,
    // Gas
    'air': null,
};

// ========== API HELPERS ==========
async function apiGet(endpoint) {
    const response = await fetch(`/api${endpoint}`);
    return response.json();
}

async function apiPost(endpoint, data) {
    const response = await fetch(`/api${endpoint}`, {
        method: 'POST',
        headers: { 'Content-Type': 'application/json' },
        body: JSON.stringify(data),
    });
    return response.json();
}

async function apiPut(endpoint, data) {
    const response = await fetch(`/api${endpoint}`, {
        method: 'PUT',
        headers: { 'Content-Type': 'application/json' },
        body: JSON.stringify(data),
    });
    return response.json();
}

async function apiDelete(endpoint) {
    const response = await fetch(`/api${endpoint}`, { method: 'DELETE' });
    return response.json();
}

// ========== SESSION MANAGEMENT ==========
async function initSession() {
    try {
        // Create a new session for this page load
        const result = await apiPost('/sessions/create', { name: 'SFPPy Studio Session' });
        if (result.success) {
            state.sessionId = result.session.id;
            console.log('Session created:', state.sessionId);
            // Add default layer to session
            await syncLayersToBackend();
            // Add default step to session
            await syncStepsToBackend();
        }
    } catch (e) {
        console.warn('Failed to create session, running in offline mode:', e);
    }
}

async function syncLayersToBackend() {
    if (!state.sessionId) return;
    try {
        // Clear existing layers and re-add
        const existing = await apiGet(`/sessions/${state.sessionId}/layers`);
        for (const layer of existing.layers || []) {
            await apiDelete(`/sessions/${state.sessionId}/layers/${layer.index}`);
        }
        // Add current layers
        for (const layer of state.layers) {
            await apiPost(`/sessions/${state.sessionId}/layers`, {
                polymer: layer.polymer,
                thickness: layer.thickness,
                thickness_unit: layer.thickness_unit,
            });
        }
    } catch (e) {
        console.warn('Failed to sync layers:', e);
    }
}

async function syncStepsToBackend() {
    if (!state.sessionId) return;
    try {
        // Clear existing steps and re-add
        const existing = await apiGet(`/sessions/${state.sessionId}/steps`);
        for (const step of existing.steps || []) {
            await apiDelete(`/sessions/${state.sessionId}/steps/${step.index}`);
        }
        // Add current steps
        for (const step of state.steps) {
            await apiPost(`/sessions/${state.sessionId}/steps`, {
                temperature: step.temperature_C,
                temperature_unit: 'degC',
                duration: step.duration,
                duration_unit: step.duration_unit,
                with_food: step.with_food,
            });
        }
    } catch (e) {
        console.warn('Failed to sync steps:', e);
    }
}

// ========== SESSION FILE MANAGEMENT ==========
async function loadExampleSessions() {
    try {
        const data = await apiGet('/sessions/files/examples');
        const container = document.getElementById('example-sessions-list');
        if (!container) return;

        if (data.examples && data.examples.length > 0) {
            container.innerHTML = data.examples.map(ex => {
                // Extract base filename without .sfppy.json extension
                const exampleId = ex.filename.replace('.sfppy.json', '').replace('.json', '');
                return `
                <div class="example-session-card p-3 border border-gray-200 dark:border-gray-600 rounded-lg
                            hover:border-sfppy-green hover:bg-green-50 dark:hover:bg-green-900/20
                            cursor-pointer transition-all" data-example="${exampleId}">
                    <div class="flex justify-between items-start">
                        <div class="flex-1 min-w-0">
                            <h4 class="font-medium text-gray-800 dark:text-white text-sm truncate">${ex.name}</h4>
                            <p class="text-xs text-gray-500 dark:text-gray-400 mt-1 line-clamp-2">${ex.description || 'No description'}</p>
                        </div>
                        <div class="ml-2 flex-shrink-0">
                            <span class="text-xs text-gray-400">${ex.filename}</span>
                        </div>
                    </div>
                    <div class="mt-2 flex gap-2 text-xs">
                        <span class="px-2 py-0.5 bg-blue-100 dark:bg-blue-900/30 text-blue-700 dark:text-blue-400 rounded">
                            ${ex.layers_count || '?'} layers
                        </span>
                        <span class="px-2 py-0.5 bg-orange-100 dark:bg-orange-900/30 text-orange-700 dark:text-orange-400 rounded">
                            ${ex.steps_count || '?'} steps
                        </span>
                        <span class="px-2 py-0.5 bg-purple-100 dark:bg-purple-900/30 text-purple-700 dark:text-purple-400 rounded">
                            ${ex.substances_count || '?'} migrants
                        </span>
                    </div>
                </div>`;
            }).join('');

            // Add click handlers
            container.querySelectorAll('.example-session-card').forEach(card => {
                card.addEventListener('click', () => loadExampleSession(card.dataset.example));
            });
        } else {
            container.innerHTML = '<div class="text-gray-500 dark:text-gray-400 text-sm p-3 text-center">No example sessions found</div>';
        }
    } catch (e) {
        console.error('Failed to load example sessions:', e);
        const container = document.getElementById('example-sessions-list');
        if (container) {
            container.innerHTML = '<div class="text-red-500 dark:text-red-400 text-sm p-3 text-center">Failed to load examples</div>';
        }
    }
}

async function loadExampleSession(exampleName) {
    try {
        console.log('Loading example session:', exampleName);
        const data = await apiGet(`/sessions/files/load/${encodeURIComponent(exampleName)}`);
        console.log('API response:', data);
        if (data.success) {
            console.log('Applying session to state...');
            applySessionToState(data.session);
            console.log('Session applied, showing banner...');
            showSessionLoadedBanner(data.session.metadata?.name || exampleName);
            showToast(`Loaded: ${data.session.metadata?.name || exampleName}`, 'success');
        } else {
            showToast(`Failed to load: ${data.error || 'Unknown error'}`, 'error');
        }
    } catch (e) {
        console.error('Failed to load example session:', e);
        console.error('Error stack:', e.stack);
        showToast(`Failed to load: ${e.message || e}`, 'error');
    }
}

function applySessionToState(session) {
    // Apply layers
    if (session.layers && session.layers.length > 0) {
        state.layers = session.layers.map((l, i) => ({
            index: l.index || i + 1,
            polymer: l.polymer || 'LDPE',
            thickness: l.thickness?.value || 100,
            thickness_unit: l.thickness?.unit || 'um',
        }));
        setLayerCount(state.layers.length);
    }

    // Apply contact steps
    if (session.contact_steps && session.contact_steps.length > 0) {
        state.steps = session.contact_steps.map((s, i) => ({
            index: s.index || i + 1,
            temperature_C: s.temperature?.value || 25,
            duration: s.duration?.value || 10,
            duration_unit: s.duration?.unit || 'days',
            with_food: s.with_food_contact !== false,
            type: s.type || 'storage',
        }));
        updateContactStepsUI();
    }

    // Apply food settings
    if (session.food) {
        state.food.category = session.food.category || 'fatty';
        state.food.simulant = session.food.simulant || 'ethanol';
        state.food.matrixType = session.food.texture || 'liquid';

        // Load CF0 (initial concentration in food)
        if (session.food.CF0) {
            state.CF0 = {};
            for (const [subId, data] of Object.entries(session.food.CF0)) {
                state.CF0[subId] = {
                    value: typeof data === 'object' ? data.value : data,
                    unit: typeof data === 'object' ? (data.unit || 'mg/kg') : 'mg/kg',
                };
            }
        } else {
            state.CF0 = {};
        }
    }

    // Apply geometry - convert all dimensions to mm
    if (session.geometry) {
        state.geometry.shape = session.geometry.shape || 'cylinder';
        if (session.geometry.dimensions) {
            // Helper to convert dimension to mm
            const toMm = (dim) => {
                if (!dim) return undefined;
                const val = dim.value || dim;
                const unit = (dim.unit || 'mm').toLowerCase();
                if (unit === 'cm') return val * 10;
                if (unit === 'm') return val * 1000;
                return val;  // assume mm
            };
            state.geometry.dimensions = {
                radius: toMm(session.geometry.dimensions.radius) || 50,
                height: toMm(session.geometry.dimensions.height) || 100,
                length: toMm(session.geometry.dimensions.length),
                width: toMm(session.geometry.dimensions.width),
            };
        }
    }

    // Apply substances
    if (session.substances && session.substances.length > 0) {
        state.substances = session.substances.map(s => ({
            id: s.id,
            name: s.properties?.name || s.lookup_name || s.id,
            cid: s.cid,
            cas: s.properties?.cas,
            mw: s.properties?.mw,
            logP: s.properties?.logP,
            formula: s.properties?.formula,
            SML: s.SML,
            color: s.color,
        }));

        // Clear and rebuild C0 matrix
        state.c0Matrix = {};
        session.layers?.forEach(layer => {
            layer.substances?.forEach(sub => {
                if (!state.c0Matrix[sub.substance_id]) {
                    state.c0Matrix[sub.substance_id] = {};
                }
                state.c0Matrix[sub.substance_id][layer.index] = {
                    value: sub.C0?.value || 0,
                    unit: sub.C0?.unit || 'mg/kg',
                };
            });
        });
    }

    // Update simulation name
    const simName = document.getElementById('sim-name');
    if (simName && session.metadata?.name) {
        simName.value = session.metadata.name;
    }

    // Update geometry UI
    if (session.geometry) {
        const shapeSelect = document.getElementById('geometry-shape');
        if (shapeSelect) {
            shapeSelect.value = session.geometry.shape || 'cylinder';
            // Trigger change event to update dimension inputs visibility
            shapeSelect.dispatchEvent(new Event('change'));
        }

        // Update dimension inputs after shape change (use converted mm values from state)
        setTimeout(() => {
            if (state.geometry.dimensions) {
                const dims = state.geometry.dimensions;
                // Set dimension values based on shape (already in mm)
                if (dims.radius) {
                    const radiusInput = document.getElementById('dim-radius');
                    if (radiusInput) radiusInput.value = dims.radius;
                }
                if (dims.height) {
                    const heightInput = document.getElementById('dim-height');
                    if (heightInput) heightInput.value = dims.height;
                }
                if (dims.length) {
                    const lengthInput = document.getElementById('dim-length');
                    if (lengthInput) lengthInput.value = dims.length;
                }
                if (dims.width) {
                    const widthInput = document.getElementById('dim-width');
                    if (widthInput) widthInput.value = dims.width;
                }
            }
        }, 100);
    }

    // Update food UI
    if (session.food) {
        const categorySelect = document.getElementById('food-category');
        if (categorySelect && session.food.category) {
            categorySelect.value = session.food.category;
        }
        const simulantSelect = document.getElementById('food-simulant');
        if (simulantSelect && session.food.simulant) {
            simulantSelect.value = session.food.simulant;
        }
        const textureSelect = document.getElementById('food-texture');
        if (textureSelect && session.food.texture) {
            textureSelect.value = session.food.texture;
        }
    }

    // Refresh UI
    updateLayerVisualization();
    updateLayerTable();
    updateSubstancesList();
    updateAssemblySubstances();
    updateSimulateSummary();

    // Navigate to Substances tab after loading session
    const substancesTab = document.querySelector('.tab-btn[data-tab="substances"]');
    if (substancesTab) {
        substancesTab.click();
    }
}

function updateContactStepsUI() {
    const container = document.getElementById('contact-steps');
    if (!container) return;

    // Remove existing step cards except the first template
    const existingCards = container.querySelectorAll('.step-card');
    existingCards.forEach((card, i) => {
        if (i > 0) card.remove();
    });

    // Update or add steps
    state.steps.forEach((step, i) => {
        if (i === 0) {
            // Update first step
            const card = container.querySelector('.step-card[data-step="1"]');
            if (card) {
                const tempInput = card.querySelector('.step-temp');
                const durationInput = card.querySelector('.step-duration');
                const durationUnit = card.querySelector('.step-duration-unit');
                const typeSelect = card.querySelector('.step-type');
                const foodCheck = card.querySelector('.step-with-food');

                if (tempInput) tempInput.value = step.temperature_C;
                if (durationInput) durationInput.value = step.duration;
                if (durationUnit) durationUnit.value = step.duration_unit;
                if (typeSelect) typeSelect.value = step.type;
                if (foodCheck) foodCheck.checked = step.with_food;
            }
        } else {
            // Add additional steps
            addContactStepCard(step);
        }
    });
}

function addContactStepCard(stepData) {
    const container = document.getElementById('contact-steps');
    const stepIndex = stepData.index || state.steps.length;

    const card = document.createElement('div');
    card.className = 'border dark:border-gray-700 rounded-lg p-3 mb-2 step-card';
    card.dataset.step = stepIndex;

    card.innerHTML = `
        <div class="flex items-center justify-between mb-2">
            <span class="font-medium text-sm text-gray-800 dark:text-white">Step ${stepIndex}</span>
            <div class="flex items-center gap-2">
                <label class="flex items-center space-x-1">
                    <input type="checkbox" class="rounded text-sfppy-green step-with-food" data-step="${stepIndex}" ${stepData.with_food ? 'checked' : ''}>
                    <span class="text-xs text-gray-600 dark:text-gray-400">Food contact</span>
                </label>
                <button class="text-red-500 hover:text-red-700 text-xs remove-step" data-step="${stepIndex}">✕</button>
            </div>
        </div>
        <div class="grid grid-cols-3 gap-2">
            <div>
                <label class="text-xs text-gray-600 dark:text-gray-400">Type</label>
                <select class="input-field text-sm step-type" data-step="${stepIndex}">
                    <option value="storage" ${stepData.type === 'storage' ? 'selected' : ''}>Storage</option>
                    <option value="hotfill" ${stepData.type === 'hotfill' ? 'selected' : ''}>Hot Fill</option>
                    <option value="transport" ${stepData.type === 'transport' ? 'selected' : ''}>Transport</option>
                    <option value="setoff" ${stepData.type === 'setoff' ? 'selected' : ''}>Set-off</option>
                </select>
            </div>
            <div>
                <label class="text-xs text-gray-600 dark:text-gray-400">Temperature</label>
                <div class="flex">
                    <input type="number" class="flex-1 input-field rounded-r-none text-sm step-temp"
                           data-step="${stepIndex}" value="${stepData.temperature_C || 25}" step="any">
                    <select class="border-l-0 input-field rounded-l-none w-14 text-sm step-temp-unit" data-step="${stepIndex}">
                        <option value="C">°C</option>
                        <option value="K">K</option>
                    </select>
                </div>
            </div>
            <div>
                <label class="text-xs text-gray-600 dark:text-gray-400">Duration</label>
                <div class="flex">
                    <input type="number" class="flex-1 input-field rounded-r-none text-sm step-duration"
                           data-step="${stepIndex}" value="${stepData.duration || 10}" min="0" step="any">
                    <select class="border-l-0 input-field rounded-l-none w-16 text-sm step-duration-unit" data-step="${stepIndex}">
                        <option value="days" ${stepData.duration_unit === 'days' ? 'selected' : ''}>days</option>
                        <option value="hours" ${stepData.duration_unit === 'hours' ? 'selected' : ''}>h</option>
                        <option value="minutes" ${stepData.duration_unit === 'minutes' ? 'selected' : ''}>min</option>
                        <option value="months" ${stepData.duration_unit === 'months' ? 'selected' : ''}>mo</option>
                    </select>
                </div>
            </div>
        </div>
    `;

    container.appendChild(card);

    // Add event listeners for all inputs
    card.querySelector('.step-type')?.addEventListener('change', (e) => {
        const step = state.steps.find(s => s.index === stepIndex);
        if (step) step.type = e.target.value;
    });

    card.querySelector('.step-temp')?.addEventListener('input', (e) => {
        const step = state.steps.find(s => s.index === stepIndex);
        if (step) step.temperature_C = parseFloat(e.target.value) || 25;
    });

    card.querySelector('.step-duration')?.addEventListener('input', (e) => {
        const step = state.steps.find(s => s.index === stepIndex);
        if (step) step.duration = parseFloat(e.target.value) || 1;
    });

    card.querySelector('.step-duration-unit')?.addEventListener('change', (e) => {
        const step = state.steps.find(s => s.index === stepIndex);
        if (step) step.duration_unit = e.target.value;
    });

    card.querySelector('.step-with-food')?.addEventListener('change', (e) => {
        const step = state.steps.find(s => s.index === stepIndex);
        if (step) step.with_food = e.target.checked;
    });

    // Add remove handler
    card.querySelector('.remove-step')?.addEventListener('click', () => {
        state.steps = state.steps.filter(s => s.index !== stepIndex);
        card.remove();
        // Re-index remaining steps
        state.steps.forEach((s, i) => s.index = i + 1);
    });
}

function extractInitialCF(sourceSession) {
    // Extract CF values from previous simulation results for restart
    const cf = {};
    if (sourceSession?.results?.substances) {
        sourceSession.results.substances.forEach(sub => {
            if (sub.CF_at_tcontact !== undefined) {
                cf[sub.substance_id] = sub.CF_at_tcontact;
            }
        });
    }
    return cf;
}

function buildCurrentSession() {
    // Build session object from current state
    return {
        version: "1.0",
        metadata: {
            name: document.getElementById('sim-name')?.value || "SFPPy Studio Session",
            description: "Session exported from SFPPy Studio",
            author: "SFPPy Studio",
            created: new Date().toISOString(),
            modified: new Date().toISOString(),
        },
        geometry: {
            shape: state.geometry.shape,
            dimensions: {
                // All dimensions stored in mm for consistency
                radius: state.geometry.dimensions.radius ? { value: state.geometry.dimensions.radius, unit: "mm" } : undefined,
                height: state.geometry.dimensions.height ? { value: state.geometry.dimensions.height, unit: "mm" } : undefined,
                length: state.geometry.dimensions.length ? { value: state.geometry.dimensions.length, unit: "mm" } : undefined,
                width: state.geometry.dimensions.width ? { value: state.geometry.dimensions.width, unit: "mm" } : undefined,
            },
        },
        substances: state.substances.map(s => ({
            id: s.id || s.name?.toLowerCase().replace(/\s+/g, '_') || `sub_${Date.now()}`,
            source: s.cid ? "pubchem" : "custom",
            lookup_name: s.name,
            cid: s.cid,
            properties: {
                name: s.name,
                cas: s.cas,
                mw: s.mw,
                logP: s.logP,
                formula: s.formula,
            },
            SML: s.SML,
            color: s.color,
        })),
        layers: state.layers.map(l => ({
            index: l.index,
            polymer: l.polymer,
            thickness: { value: l.thickness, unit: l.thickness_unit },
            name: `Layer ${l.index} - ${l.polymer}`,
            substances: Object.entries(state.c0Matrix)
                .filter(([subId, layers]) => layers[l.index]?.value > 0)
                .map(([subId, layers]) => ({
                    substance_id: subId,
                    C0: { value: layers[l.index].value, unit: layers[l.index].unit || 'mg/kg' },
                    D_auto: true,
                    k_auto: true,
                })),
        })),
        food: {
            category: state.food.category || "realfood",
            texture: state.food.matrixType || "liquid",
            affinity: state.food.category?.includes('fat') ? "fat" : "aqueous",
            simulant: state.food.simulant,
            // CF0: initial concentration in food per substance
            CF0: Object.keys(state.CF0).length > 0 ? Object.fromEntries(
                Object.entries(state.CF0).map(([subId, data]) => [subId, {
                    value: data.value,
                    unit: data.unit || 'mg/kg',
                }])
            ) : undefined,
        },
        contact_steps: state.steps.map(s => ({
            index: s.index,
            name: `Step ${s.index}`,
            type: s.type,
            temperature: { value: s.temperature_C, unit: "degC" },
            duration: { value: s.duration, unit: s.duration_unit },
            with_food_contact: s.with_food,
            setoff_type: s.setoff_type || 'stacked',
        })),
        simulation_config: {
            solver: "senspatankar",
            time_factor: 2.0,
            n_time_points: 1000,
            n_space_points: 200,
        },
        restart_state: state.restartEnabled && state.restartSource ? {
            enabled: true,
            from_step: parseInt(document.getElementById('restart-from-step')?.value || '1'),
            initial_CF: extractInitialCF(state.restartSource.session),
            source_session: state.restartSource.filename,
        } : undefined,
    };
}

function exportSessionAsJSON() {
    const session = buildCurrentSession();
    const blob = new Blob([JSON.stringify(session, null, 2)], { type: 'application/json' });
    const url = URL.createObjectURL(blob);
    const a = document.createElement('a');
    a.href = url;
    a.download = `${session.metadata.name.replace(/\s+/g, '_')}.sfppy.json`;
    document.body.appendChild(a);
    a.click();
    document.body.removeChild(a);
    URL.revokeObjectURL(url);
}

async function loadSessionFromFile(file) {
    try {
        const text = await file.text();
        const session = JSON.parse(text);
        applySessionToState(session);
        showSessionLoadedBanner(session.metadata?.name || file.name);
        showToast(`Loaded: ${session.metadata?.name || file.name}`, 'success');
    } catch (e) {
        console.error('Failed to parse session file:', e);
        showToast('Failed to load session file: Invalid JSON', 'error');
    }
}

// Show session loaded banner indicator
function showSessionLoadedBanner(sessionName) {
    // Update or create the session loaded indicator
    let banner = document.getElementById('session-loaded-banner');

    if (!banner) {
        // Create banner after the Session Files section header
        const sessionSection = document.querySelector('#tab-simulate > div:first-child');
        if (sessionSection) {
            const h2 = sessionSection.querySelector('h2');
            if (h2) {
                banner = document.createElement('div');
                banner.id = 'session-loaded-banner';
                banner.className = 'mb-3 p-2 bg-green-100 dark:bg-green-900/30 border border-green-300 dark:border-green-700 rounded-lg flex items-center justify-between';
                h2.insertAdjacentElement('afterend', banner);
            }
        }
    }

    if (banner) {
        banner.innerHTML = `
            <div class="flex items-center gap-2 text-green-800 dark:text-green-200">
                <span class="text-lg">✓</span>
                <span class="font-medium text-sm">Session loaded:</span>
                <span class="text-sm">${sessionName}</span>
            </div>
            <button onclick="clearSessionBanner()" class="text-green-600 dark:text-green-400 hover:text-green-800 text-sm px-2">×</button>
        `;
        banner.classList.remove('hidden');
    }

    // Store loaded session name in state
    state.loadedSessionName = sessionName;
}

function clearSessionBanner() {
    const banner = document.getElementById('session-loaded-banner');
    if (banner) {
        banner.classList.add('hidden');
    }
    state.loadedSessionName = null;
}

function initSessionFileHandlers() {
    // Load session button
    const loadBtn = document.getElementById('load-session-btn');
    const fileInput = document.getElementById('session-file-input');
    if (loadBtn && fileInput) {
        loadBtn.addEventListener('click', () => fileInput.click());
        fileInput.addEventListener('change', (e) => {
            if (e.target.files && e.target.files[0]) {
                loadSessionFromFile(e.target.files[0]);
            }
        });
    }

    // Export JSON button
    const exportBtn = document.getElementById('export-session-json-btn');
    if (exportBtn) {
        exportBtn.addEventListener('click', exportSessionAsJSON);
    }

    // Save session button (placeholder - could save to backend)
    const saveBtn = document.getElementById('save-session-btn');
    if (saveBtn) {
        saveBtn.addEventListener('click', () => {
            // For now, just export as JSON
            exportSessionAsJSON();
        });
    }

    // Restart options toggle
    const restartCheckbox = document.getElementById('restart-enabled');
    const restartOptions = document.getElementById('restart-options');
    if (restartCheckbox && restartOptions) {
        restartCheckbox.addEventListener('change', () => {
            restartOptions.classList.toggle('hidden', !restartCheckbox.checked);
            state.restartEnabled = restartCheckbox.checked;
        });
    }

    // Restart browse button
    const restartBrowseBtn = document.getElementById('restart-browse-btn');
    if (restartBrowseBtn) {
        restartBrowseBtn.addEventListener('click', () => {
            // Create a temporary file input for restart session
            const input = document.createElement('input');
            input.type = 'file';
            input.accept = '.json,.sfppy.json';
            input.onchange = (e) => {
                if (e.target.files && e.target.files[0]) {
                    loadRestartSourceFile(e.target.files[0]);
                }
            };
            input.click();
        });
    }
}

async function loadRestartSourceFile(file) {
    try {
        const text = await file.text();
        const session = JSON.parse(text);

        // Check if session has results to restart from
        if (!session.results || !session.results.substances) {
            alert('Selected session has no simulation results to restart from.');
            return;
        }

        // Store restart data in state
        state.restartSource = {
            filename: file.name,
            session: session,
        };

        // Update UI - show filename of loaded session
        const sourceInput = document.getElementById('restart-source');
        if (sourceInput) {
            sourceInput.value = file.name;
        }

        // Note: The step dropdown shows CURRENT simulation's steps (where to restart TO),
        // not the loaded session's steps. The loaded session provides initial CF values.
        // updateRestartStepOptions() handles this via updateSimulateSummary().

        // Show info about loaded session
        const sourceSteps = session.contact_steps?.length || 0;
        const sourceSubstances = session.results?.substances?.length || 0;
        showToast(`Loaded restart source: ${session.metadata?.name || file.name} (${sourceSteps} steps, ${sourceSubstances} substances)`, 'success');
    } catch (e) {
        console.error('Failed to parse restart source file:', e);
        alert('Failed to load restart source file: Invalid JSON');
    }
}

// ========== D/K COMPUTATION ==========
async function computeDiffusivity(substanceName, polymer, temperature_C = 40) {
    try {
        const params = new URLSearchParams({
            substance: substanceName,
            polymer: polymer,
            temperature_C: temperature_C.toString(),
        });
        const result = await apiGet(`/sessions/compute/diffusivity?${params}`);
        if (result.success) {
            return result.properties;
        }
    } catch (e) {
        console.warn('Failed to compute diffusivity:', e);
    }
    return null;
}

async function computeDiffusivityBatch(substanceName, polymers = ['LDPE', 'PP', 'gPET'], temperature_C = 40) {
    try {
        const params = new URLSearchParams({
            substance: substanceName,
            polymers: polymers.join(','),
            temperature_C: temperature_C.toString(),
        });
        const result = await apiGet(`/sessions/compute/diffusivity/batch?${params}`);
        if (result.success) {
            return result.results;
        }
    } catch (e) {
        console.warn('Failed to compute batch diffusivity:', e);
    }
    return null;
}

// ========== DENSITY FUNCTIONS ==========

/**
 * Fetch layer density from API with temperature dependence
 * Uses patankar layer.density(T) method
 */
async function fetchLayerDensity(polymerCode, temperature_C = 25) {
    try {
        const params = new URLSearchParams({ temperature_C: temperature_C.toString() });
        const result = await apiGet(`/assembly/materials/${polymerCode}/density?${params}`);
        if (result.success) {
            return {
                rho: result.rho,
                unit: result.rho_unit,
                source: result.source,
                is_temperature_dependent: result.is_temperature_dependent
            };
        }
    } catch (e) {
        console.warn(`Failed to fetch density for ${polymerCode}:`, e);
    }
    // Fallback defaults
    const defaults = {
        'LDPE': 920, 'LLDPE': 920, 'HDPE': 960, 'PP': 910,
        'gPET': 1380, 'wPET': 1380, 'rPET': 1380, 'PET': 1380,
        'PS': 1050, 'PMMA': 1190, 'PA6': 1130, 'PA66': 1140
    };
    return { rho: defaults[polymerCode] || 1000, unit: 'kg/m³', source: 'default', is_temperature_dependent: false };
}

/**
 * Fetch food/simulant density from API
 */
async function fetchFoodDensity(simulantCode, temperature_C = 25) {
    try {
        const params = new URLSearchParams({ temperature_C: temperature_C.toString() });
        const result = await apiGet(`/food/simulants/density/${simulantCode}?${params}`);
        if (result.success) {
            // API returns density_kg_m3 (T-dependent) and density_25C_kg_m3 (reference)
            const hasTempDependence = result.temp_coeff && result.temp_coeff !== 0;
            return {
                rho: result.density_kg_m3,
                unit: result.unit || 'kg/m³',
                name: result.code,
                is_temperature_dependent: hasTempDependence
            };
        }
    } catch (e) {
        console.warn(`Failed to fetch food density for ${simulantCode}:`, e);
    }
    // Fallback defaults
    const defaults = {
        'water': 997, 'ethanol_10': 982, 'ethanol_20': 969,
        'ethanol_50': 914, 'ethanol_95': 789, 'ethanol': 789,
        'olive_oil': 920, 'vegetable_oil': 920, 'iso_octane': 690
    };
    return { rho: defaults[simulantCode] || 1000, unit: 'kg/m³', name: simulantCode, is_temperature_dependent: false };
}

/**
 * Fetch densities for all layers and update display
 */
async function fetchAndUpdateLayerDensities() {
    const temperature = state.steps[0]?.temperature_C || 25;

    for (const layer of state.layers) {
        const densityData = await fetchLayerDensity(layer.polymer, temperature);
        layer.rho_computed = densityData.rho;
        layer.rho_source = densityData.source;
        layer.rho_T_dependent = densityData.is_temperature_dependent;

        // Update display in table
        const rhoCell = document.querySelector(`.layer-rho-display[data-layer="${layer.index}"]`);
        if (rhoCell) {
            rhoCell.textContent = Math.round(densityData.rho);
            rhoCell.title = `ρ at ${temperature}°C (${densityData.source})`;
        }
    }

    // Update detail panel if showing
    if (state.selectedLayerIndex) {
        updateLayerDensityDisplay(state.selectedLayerIndex);
    }
}

/**
 * Update density display in layer detail panel
 */
function updateLayerDensityDisplay(layerIndex) {
    const layer = state.layers.find(l => l.index === layerIndex);
    if (!layer) return;

    const temperature = state.steps[0]?.temperature_C || 25;
    const densityValueEl = document.getElementById('layer-density-value');
    const densityNoteEl = document.getElementById('layer-density-note');
    const rhoComputedEl = document.getElementById('layer-rho-computed');

    if (densityValueEl) {
        const rho = layer.rho_computed || '-';
        densityValueEl.textContent = `${typeof rho === 'number' ? Math.round(rho) : rho} kg/m³`;
    }

    if (densityNoteEl) {
        if (layer.rho_T_dependent) {
            densityNoteEl.textContent = `Temperature-dependent: ρ(${temperature}°C) from patankar`;
        } else {
            densityNoteEl.textContent = `Source: ${layer.rho_source || 'default'}`;
        }
    }

    if (rhoComputedEl) {
        rhoComputedEl.textContent = layer.rho_computed ? `${Math.round(layer.rho_computed)} kg/m³` : '-';
    }
}

/**
 * Update food density display in Food tab
 */
async function updateFoodDensityDisplay() {
    const temperature = state.steps[0]?.temperature_C || 25;
    // Use state.food.simulant (the actual selection), not state.foodConfig.simulant
    const simulant = state.food?.simulant || 'ethanol';

    const densityData = await fetchFoodDensity(simulant, temperature);

    // Update display
    const densityValueEl = document.getElementById('food-density-value');
    const simulantNameEl = document.getElementById('food-simulant-name');
    const densityNoteEl = document.getElementById('food-density-note');

    if (densityValueEl) {
        densityValueEl.textContent = Math.round(densityData.rho);
    }

    if (simulantNameEl) {
        simulantNameEl.textContent = densityData.name || simulant;
    }

    if (densityNoteEl) {
        if (densityData.is_temperature_dependent) {
            densityNoteEl.textContent = `Temperature-dependent: ρ(${temperature}°C) = ρ₂₅ × (1 - α × (T - 25))`;
        } else {
            densityNoteEl.textContent = `Reference density at 25°C`;
        }
    }

    // Store for use in simulation
    state.foodConfig = state.foodConfig || {};
    state.foodConfig.density = densityData.rho;
}

// ========== UTILITY FUNCTIONS ==========
function formatScientific(value, precision = 2) {
    if (value === null || value === undefined) return '-';
    if (value === 0) return '0';
    const exp = Math.floor(Math.log10(Math.abs(value)));
    const mantissa = value / Math.pow(10, exp);
    return `${mantissa.toFixed(precision)}e${exp}`;
}

function formatNumber(value, precision = 2) {
    if (value === null || value === undefined) return '-';
    if (Math.abs(value) < 1e-6 || Math.abs(value) > 1e6) {
        return formatScientific(value, precision);
    }
    return value.toFixed(precision);
}

/**
 * Get structure image URL - uses PubChem thumbnails as primary (fast, small)
 * @param {number} cid - PubChem CID
 * @param {string} size - 's' for small thumbnail, 'l' for large
 * @returns {string} Image URL
 */
function getStructureImageUrl(cid, size = 's') {
    if (!cid) return '';
    // Use PubChem thumbnails as primary (fast, browser-cached)
    return `https://pubchem.ncbi.nlm.nih.gov/image/imgsrv.fcgi?cid=${cid}&t=${size}`;
}

/**
 * Get local SFPPy structure image URL
 * @param {number} cid - PubChem CID
 * @returns {string} Local image URL
 */
function getLocalStructureImageUrl(cid) {
    if (!cid) return '';
    return `/api/substances/structure/${cid}.png`;
}

/**
 * Set image source - PubChem thumbnail primary, local SFPPy as fallback
 * @param {HTMLImageElement} img - Image element
 * @param {number} cid - PubChem CID
 * @param {string} size - 's' for small, 'l' for large
 */
function setStructureImage(img, cid, size = 's') {
    if (!cid || !img) return;

    const pubchemUrl = `https://pubchem.ncbi.nlm.nih.gov/image/imgsrv.fcgi?cid=${cid}&t=${size}`;
    const localUrl = `/api/substances/structure/${cid}.png`;

    // Store CID for click handler
    img.dataset.cid = cid;

    // Primary: PubChem thumbnail (fast, cached by browser/CDN)
    img.src = pubchemUrl;
    img.onerror = function() {
        // Fallback to local SFPPy cache if PubChem unavailable (offline)
        if (this.src !== localUrl) {
            console.log('[Image] PubChem unavailable, using local cache for CID:', cid);
            this.src = localUrl;
            this.onerror = null; // Prevent infinite loop
        }
    };

    // Add click handler to show enlarged local image in modal
    if (!img.dataset.clickHandlerAttached) {
        img.style.cursor = 'pointer';
        img.title = 'Click to enlarge';
        img.addEventListener('click', function(e) {
            e.stopPropagation();
            showStructureModal(cid);
        });
        img.dataset.clickHandlerAttached = 'true';
    }
}

/**
 * Show structure image in a floating modal
 * @param {number} cid - PubChem CID
 */
function showStructureModal(cid) {
    if (!cid) return;

    // Create or get modal
    let modal = document.getElementById('structure-modal');
    if (!modal) {
        modal = document.createElement('div');
        modal.id = 'structure-modal';
        modal.className = 'fixed inset-0 z-50 hidden items-center justify-center bg-black bg-opacity-50';
        modal.innerHTML = `
            <div class="bg-white dark:bg-gray-800 rounded-lg shadow-xl p-4 max-w-lg mx-4 relative">
                <button id="structure-modal-close" class="absolute top-2 right-2 text-gray-500 hover:text-gray-700 dark:text-gray-400 dark:hover:text-gray-200 text-2xl font-bold">&times;</button>
                <div class="text-center">
                    <img id="structure-modal-img" class="max-w-full max-h-[70vh] mx-auto" alt="Structure">
                    <p id="structure-modal-cid" class="mt-2 text-sm text-gray-600 dark:text-gray-400"></p>
                    <p id="structure-modal-source" class="text-xs text-gray-500 dark:text-gray-500"></p>
                </div>
            </div>
        `;
        document.body.appendChild(modal);

        // Close handlers
        modal.addEventListener('click', function(e) {
            if (e.target === modal) {
                closeStructureModal();
            }
        });
        document.getElementById('structure-modal-close').addEventListener('click', closeStructureModal);
        document.addEventListener('keydown', function(e) {
            if (e.key === 'Escape' && !modal.classList.contains('hidden')) {
                closeStructureModal();
            }
        });
    }

    const img = document.getElementById('structure-modal-img');
    const cidLabel = document.getElementById('structure-modal-cid');
    const sourceLabel = document.getElementById('structure-modal-source');

    // Try local SFPPy image first (higher quality for modal)
    const localUrl = `/api/substances/structure/${cid}.png`;
    const pubchemLargeUrl = `https://pubchem.ncbi.nlm.nih.gov/image/imgsrv.fcgi?cid=${cid}&t=l`;

    img.src = localUrl;
    sourceLabel.textContent = 'Source: SFPPy Local Cache';

    img.onerror = function() {
        if (this.src !== pubchemLargeUrl) {
            console.log('[Modal] Local not found, using PubChem large for CID:', cid);
            this.src = pubchemLargeUrl;
            sourceLabel.textContent = 'Source: PubChem';
            this.onerror = null;
        }
    };

    cidLabel.textContent = `CID: ${cid}`;

    // Show modal
    modal.classList.remove('hidden');
    modal.classList.add('flex');
}

/**
 * Close the structure modal
 */
function closeStructureModal() {
    const modal = document.getElementById('structure-modal');
    if (modal) {
        modal.classList.add('hidden');
        modal.classList.remove('flex');
    }
}

// ========== TAB NAVIGATION ==========
function initTabs() {
    const tabs = document.querySelectorAll('.tab-btn');
    const panels = document.querySelectorAll('.tab-panel');

    tabs.forEach(tab => {
        tab.addEventListener('click', () => {
            // Update buttons
            tabs.forEach(t => t.classList.remove('active'));
            tab.classList.add('active');

            // Update panels - handle both tab-* and panel-* naming conventions
            const tabName = tab.dataset.tab;
            const targetId = tabName === 'help' ? 'panel-help' : `tab-${tabName}`;
            panels.forEach(p => p.classList.add('hidden'));
            const targetPanel = document.getElementById(targetId);
            if (targetPanel) targetPanel.classList.remove('hidden');

            // Load data for specific tabs
            if (tabName === 'food') loadFoodData();
            if (tabName === 'substances') loadSubstancesData();
            if (tabName === 'simulate') updateSimulateSummary();
            if (tabName === 'jobs') loadJobsList();
            if (tabName === 'config') loadConfig();

            // Update contextual help bar for this tab
            updateHelpBarForTab(tabName);
        });
    });
}

// ========== HELP PANEL ==========
function initHelp() {
    const helpPanel = document.getElementById('help-panel');
    const helpToggle = document.getElementById('help-toggle');
    const closeHelp = document.getElementById('close-help');
    const helpContent = document.getElementById('help-content');

    // Toggle both help panel and help bar
    helpToggle?.addEventListener('click', () => {
        helpPanel?.classList.toggle('collapsed');

        // Also toggle help bar visibility
        const helpBar = document.getElementById('help-bar');
        if (helpBar) {
            if (helpBarVisible) {
                hideHelpBar();
            } else {
                // Re-enable help bar (clear dismissal)
                localStorage.removeItem('helpBarDismissed');
                showHelpBar();
            }
        }
    });

    closeHelp?.addEventListener('click', () => {
        helpPanel?.classList.add('collapsed');
    });

    // Help icons
    document.querySelectorAll('.help-icon').forEach(icon => {
        icon.addEventListener('click', async (e) => {
            const helpKey = e.target.dataset.help;
            if (!helpKey) return;

            const [tab, element] = helpKey.split('/');
            try {
                const data = await apiGet(`/help/context/${tab}/${element}`);
                if (data.success) {
                    helpPanel.classList.remove('collapsed');
                    helpContent.innerHTML = `
                        <h4 class="font-semibold text-gray-800 mb-2">${data.title}</h4>
                        <p class="text-sm text-gray-600 mb-3">${data.long || data.short}</p>
                        ${data.tips ? `
                            <h5 class="text-xs font-medium text-gray-700 mb-1">Tips:</h5>
                            <ul class="text-xs text-gray-600 list-disc list-inside">
                                ${data.tips.map(t => `<li>${t}</li>`).join('')}
                            </ul>
                        ` : ''}
                    `;
                }
            } catch (e) {
                console.error('Failed to load help:', e);
            }
        });
    });

    // Initialize help bar
    initHelpBar();
}

// ========== HELP BAR SYSTEM ==========
const HELP_TIPS = {
    substances: [
        { icon: '🧪', text: 'Start by searching for substances in PubChem using CAS numbers, IUPAC names, or common names like "BHT" or "Irganox 1076".' },
        { icon: '📋', text: 'Click on common additives to quickly add them to your simulation (antioxidants, plasticizers, UV stabilizers).' },
        { icon: '⚖️', text: 'Check the SML column for regulatory Specific Migration Limits from EU 10/2011, US FDA, or Chinese GB 9685.' },
        { icon: '🔬', text: 'Cramer class (I/II/III) is from ToxTree. Class I = low concern (TTC: 30 µg/kg/day), Class III = high concern (TTC: 1.5 µg/kg/day).' },
        { icon: '⚠️', text: 'For NIAS (no SML), use TTC-based limits. Max Exposure = TTC × 60 kg ÷ 1000 gives mg/kg food. Check for structural alerts!' },
        { icon: '☢️', text: 'GENOTOXIC alerts (purple badge) require TTC = 0.0025 µg/kg bw/day → 0.00015 mg/kg food. Individual risk assessment mandatory.' },
        { icon: '🟠', text: 'Orange alerts are informational (skin sensitization, etc.) - relevant for cosmetics/homecare. Cramer TTC still applies for food contact.' },
    ],
    assembly: [
        { icon: '🧱', text: 'Build your packaging from food-contact layer (Layer 1) outward. Layer 1 is always in contact with food.' },
        { icon: '📏', text: 'Set layer thicknesses and initial concentrations (C₀) for each substance in the matrix below.' },
        { icon: '📦', text: 'Choose the packaging geometry to calculate the correct surface-to-volume ratio for migration.' },
        { icon: '⚡', text: 'Diffusivity (D) and partition coefficient (k) are auto-computed using the Piringer model. Override if you have measured values.' },
    ],
    food: [
        { icon: '🍽️', text: 'Select a food category and simulant. EU simulants: water (aqueous), 10% ethanol (alcohol), 95% ethanol (fatty/dry), olive oil (fatty).' },
        { icon: '⏱️', text: 'Define contact steps with temperature and duration. Use multiple steps for scenarios like hot-fill → storage.' },
        { icon: '📦', text: 'Set-off mode simulates stacking/rolling before food contact (no food in this step). Useful for realistic worst-case scenarios.' },
        { icon: '🌡️', text: 'Higher temperatures increase diffusivity exponentially. A hot-fill at 80°C can dominate migration even if brief.' },
    ],
    simulate: [
        { icon: '▶️', text: 'Review the summary counters (layers, steps, substances) before running. All must be > 0.' },
        { icon: '✅', text: 'Click "Validate" to check your configuration for errors before running the simulation.' },
        { icon: '📂', text: 'Load example sessions to see complete configurations for common scenarios.' },
        { icon: '🔄', text: 'Use "Restart" to continue from a previous simulation\'s final state for multi-stage contact scenarios.' },
    ],
    results: [
        { icon: '📊', text: 'The kinetics chart shows CF (food concentration) vs time. Final CF is compared against SML.' },
        { icon: '📈', text: 'Concentration profiles show how C evolves within the packaging layers over time.' },
        { icon: '✅', text: 'Green badge = compliant (CF < SML). Red badge = exceeds SML. Yellow = no SML defined.' },
        { icon: '📄', text: 'Export results as PDF reports, CSV data files, or publication-quality figures.' },
    ],
    jobs: [
        { icon: '📁', text: 'All simulations are saved here with unique job IDs for reproducibility.' },
        { icon: '🔄', text: 'Click on a job to reload its results. Use "Rerun" to repeat with modified parameters.' },
        { icon: '💾', text: 'Download session files (.json) to share configurations or archive for later.' },
    ],
    config: [
        { icon: '⚙️', text: 'Adjust solver settings like mesh resolution and time steps for accuracy vs. speed tradeoff.' },
        { icon: '🔢', text: 'Higher mesh resolution (Nx) improves accuracy for multilayer systems but increases computation time.' },
        { icon: '🖥️', text: 'System info shows SFPPy version and available modules for debugging.' },
    ],
    fitting: [
        { icon: '📈', text: 'Fit D (diffusivity) and k (partition coefficient) from experimental migration data.' },
        { icon: '🧪', text: 'Generate synthetic data with known parameters to test the fitting algorithm.' },
        { icon: '📁', text: 'Upload CSV files with time (days) and CF (concentration in food) columns.' },
    ],
    help: [
        { icon: '❓', text: 'This tab provides an overview of the SFPPy Studio workflow and features.' },
        { icon: '📧', text: 'Questions or issues? Contact olivier.vitrac@gmail.com or visit the GitHub repository.' },
        { icon: '📜', text: 'SFPPy is open-source under the MIT License. Results are provided AS IS.' },
    ],
};

let currentHelpTipIndex = 0;
let currentHelpTab = 'substances';
let helpBarVisible = false;

function initHelpBar() {
    const helpBar = document.getElementById('help-bar');
    const toggleBtn = document.getElementById('help-bar-toggle');
    const closeBtn = document.getElementById('help-bar-close');
    const prevBtn = document.getElementById('help-bar-prev');
    const nextBtn = document.getElementById('help-bar-next');
    const firstBtn = document.getElementById('help-bar-first');
    const lastBtn = document.getElementById('help-bar-last');

    if (!helpBar) return;

    // Check if user has dismissed help bar before
    const helpBarDismissed = localStorage.getItem('helpBarDismissed');
    if (!helpBarDismissed) {
        // Show help bar after a short delay
        setTimeout(() => {
            showHelpBar();
        }, 1000);
    }

    // Toggle button (floating 💡 button) - always toggles
    toggleBtn?.addEventListener('click', () => {
        if (helpBarVisible) {
            hideHelpBar();
        } else {
            showHelpBar();
        }
    });

    // Close/minimize button (▼) - just hides, can be reopened
    closeBtn?.addEventListener('click', () => {
        hideHelpBar();
        // Don't set localStorage - user can easily reopen with toggle button
    });

    // Navigation buttons
    firstBtn?.addEventListener('click', () => {
        currentHelpTipIndex = 0;
        updateHelpBarContent();
    });

    prevBtn?.addEventListener('click', () => {
        const tips = HELP_TIPS[currentHelpTab] || HELP_TIPS.substances;
        currentHelpTipIndex = (currentHelpTipIndex - 1 + tips.length) % tips.length;
        updateHelpBarContent();
    });

    nextBtn?.addEventListener('click', () => {
        const tips = HELP_TIPS[currentHelpTab] || HELP_TIPS.substances;
        currentHelpTipIndex = (currentHelpTipIndex + 1) % tips.length;
        updateHelpBarContent();
    });

    lastBtn?.addEventListener('click', () => {
        const tips = HELP_TIPS[currentHelpTab] || HELP_TIPS.substances;
        currentHelpTipIndex = tips.length - 1;
        updateHelpBarContent();
    });
}

function showHelpBar() {
    const helpBar = document.getElementById('help-bar');
    const toggleBtn = document.getElementById('help-bar-toggle');
    if (!helpBar) return;

    helpBar.style.display = 'block';
    // Trigger animation
    setTimeout(() => {
        helpBar.classList.remove('translate-y-full');
    }, 10);
    helpBarVisible = true;

    // Update content and counter
    updateHelpBarContent();

    // Update toggle button state
    if (toggleBtn) {
        toggleBtn.classList.add('ring-2', 'ring-white', 'ring-offset-2');
        toggleBtn.title = 'Hide Help Tips';
    }
}

function hideHelpBar() {
    const helpBar = document.getElementById('help-bar');
    const toggleBtn = document.getElementById('help-bar-toggle');
    if (!helpBar) return;

    helpBar.classList.add('translate-y-full');
    setTimeout(() => {
        helpBar.style.display = 'none';
    }, 300);
    helpBarVisible = false;

    // Update toggle button state
    if (toggleBtn) {
        toggleBtn.classList.remove('ring-2', 'ring-white', 'ring-offset-2');
        toggleBtn.title = 'Show Help Tips';
    }
}

function updateHelpBarForTab(tabName) {
    currentHelpTab = tabName;
    currentHelpTipIndex = 0;
    updateHelpBarContent();
}

function updateHelpBarContent() {
    const tips = HELP_TIPS[currentHelpTab] || HELP_TIPS.substances;
    const tip = tips[currentHelpTipIndex];

    const iconEl = document.getElementById('help-bar-icon');
    const textEl = document.getElementById('help-bar-text');
    const counterEl = document.getElementById('help-bar-counter');

    if (iconEl && tip) iconEl.textContent = tip.icon;
    if (textEl && tip) textEl.textContent = tip.text;
    if (counterEl) counterEl.textContent = `${currentHelpTipIndex + 1}/${tips.length}`;
}

// ========== ASSEMBLY TAB ==========
function initAssembly() {
    // Layer count buttons
    document.querySelectorAll('.layer-count-btn').forEach(btn => {
        btn.addEventListener('click', () => {
            const count = parseInt(btn.dataset.count);
            setLayerCount(count);

            // Update button styles
            document.querySelectorAll('.layer-count-btn').forEach(b => {
                b.classList.remove('bg-blue-600', 'text-white');
                b.classList.add('bg-white');
            });
            btn.classList.add('bg-blue-600', 'text-white');
            btn.classList.remove('bg-white');
        });
    });

    // Initialize first layer listeners
    initLayerListeners(1);

    // Validate button
    document.getElementById('assembly-validate').addEventListener('click', validateAssembly);

    // Reset button
    document.getElementById('assembly-reset').addEventListener('click', resetAssembly);

    // Initialize C0 profile histogram toggle
    initC0ProfileToggle();
}

async function setLayerCount(count) {
    // Adjust state
    while (state.layers.length < count) {
        state.layers.push({
            index: state.layers.length + 1,
            polymer: 'LDPE',
            thickness: 100,
            thickness_unit: 'um',
        });
    }
    state.layers = state.layers.slice(0, count);

    // Update layer table UI
    updateLayerTable();
    updateLayerVisualization();
    updateC0Matrix();

    // Sync to backend
    await syncLayersToBackend();
}

function updateLayerTable() {
    const tbody = document.getElementById('layer-table-body');
    if (!tbody) return;

    // Get current temperature for density calculation
    const temperature = state.steps[0]?.temperature_C || 25;

    tbody.innerHTML = state.layers.map(layer => {
        const tg = MATERIAL_TG[layer.polymer];
        const tgDisplay = tg !== null && tg !== undefined ? `${tg}°C` : '-';
        const isSelected = layer.index === state.selectedLayerIndex;
        // Get cached density or use default
        const rhoDisplay = layer.rho_computed ? Math.round(layer.rho_computed) : '-';

        return `
            <tr class="${isSelected ? 'selected' : ''} cursor-pointer layer-row" data-layer="${layer.index}">
                <td class="font-medium ${layer.index === 1 ? 'text-sfppy-green' : ''} layer-index-cell">${layer.index}</td>
                <td onclick="event.stopPropagation()">
                    <select class="input-field text-sm layer-polymer" data-layer="${layer.index}">
                        ${getMaterialOptionsHTML(layer.polymer)}
                    </select>
                </td>
                <td onclick="event.stopPropagation()">
                    <div class="flex">
                        <input type="number" class="flex-1 input-field rounded-r-none text-sm layer-thickness"
                               data-layer="${layer.index}" value="${layer.thickness}" min="0.001" step="any">
                        <select class="border-l-0 input-field rounded-l-none w-14 text-xs layer-thickness-unit" data-layer="${layer.index}">
                            <option value="um" ${layer.thickness_unit === 'um' ? 'selected' : ''}>µm</option>
                            <option value="mm" ${layer.thickness_unit === 'mm' ? 'selected' : ''}>mm</option>
                            <option value="nm" ${layer.thickness_unit === 'nm' ? 'selected' : ''}>nm</option>
                        </select>
                    </div>
                </td>
                <td class="text-xs layer-tg-display" data-layer="${layer.index}">${tgDisplay}</td>
                <td class="text-xs layer-rho-display" data-layer="${layer.index}" title="ρ at ${temperature}°C">${rhoDisplay}</td>
                <td><span class="override-badge computed" id="layer-${layer.index}-override-badge">Computed</span></td>
            </tr>
        `;
    }).join('');

    // Fetch densities asynchronously and update display
    fetchAndUpdateLayerDensities();

    // Re-init listeners for all layers
    state.layers.forEach(layer => initLayerTableListeners(layer.index));

    // Add row click listeners for selection (on the row itself)
    tbody.querySelectorAll('.layer-row').forEach(row => {
        row.addEventListener('click', (e) => {
            // Only select if clicking on the row itself or index cell
            if (e.target.tagName === 'TR' || e.target.classList.contains('layer-index-cell') ||
                e.target.classList.contains('layer-tg-display') || e.target.classList.contains('override-badge')) {
                const layerIndex = parseInt(row.dataset.layer);
                selectLayer(layerIndex);
            }
        });
    });

    // Update assembly substances matrix
    updateAssemblySubstances();
}

// Update the assembly tab substance assignments matrix
// Shows substance×layer matrix with presence toggle, C0, unit, D, k per layer
// D and k have auto-compute checkboxes - when checked, values are recalculated on polymer/temp change
async function updateAssemblySubstances() {
    const noSubstancesDiv = document.getElementById('assembly-no-substances');
    const matrixDiv = document.getElementById('assembly-substances-matrix');
    const tbody = document.getElementById('assembly-substances-body');

    if (!noSubstancesDiv || !matrixDiv || !tbody) return;

    if (state.substances.length === 0) {
        noSubstancesDiv.classList.remove('hidden');
        matrixDiv.classList.add('hidden');
        renderC0ProfileHistogram(); // Show "no substances" message in histogram
        return;
    }

    noSubstancesDiv.classList.add('hidden');
    matrixDiv.classList.remove('hidden');

    // Build rows: one row per substance × layer combination
    const rows = [];
    const temperature = state.steps[0]?.temperature_C || 40;

    for (const substance of state.substances) {
        const subId = substance.id || substance.cid || substance.name;
        const subName = substance.name || subId;
        const smlValue = substance.SML || substance.regulatory?.EU?.SML || null;

        // For each layer, show a row
        for (const layer of state.layers) {
            const layerIndex = layer.index;

            // Get current assignment (if any)
            const c0Data = state.c0Matrix[subId]?.[layerIndex];
            const isAssigned = c0Data !== undefined && (typeof c0Data === 'object' ? c0Data.value > 0 : c0Data > 0);
            const c0Value = typeof c0Data === 'object' ? c0Data.value : (c0Data || 0);
            const c0Unit = state.c0Units[subId] || state.defaultC0Unit;

            // Check auto-compute flags (default: true = auto-compute)
            const D_auto = substance[`D_auto_L${layerIndex}`] !== false;
            const k_auto = substance[`k_auto_L${layerIndex}`] !== false;

            // Get override values (only used if auto is disabled)
            const D_override = substance[`D_override_L${layerIndex}`] || null;
            const k_override = substance[`k_override_L${layerIndex}`] || null;

            // Compute D/k from backend (for display and when auto is enabled)
            let D_computed = null;
            let k_computed = null;

            if (isAssigned) {
                try {
                    const props = await computeDiffusivity(subName, layer.polymer, temperature);
                    if (props) {
                        D_computed = props.D;
                        k_computed = props.k;
                        // Store computed values on substance for reference
                        substance[`D_computed_L${layerIndex}`] = D_computed;
                        substance[`k_computed_L${layerIndex}`] = k_computed;
                    }
                } catch (e) {
                    // Use previously cached values if available
                    D_computed = substance[`D_computed_L${layerIndex}`] || null;
                    k_computed = substance[`k_computed_L${layerIndex}`] || null;
                }
            }

            // Determine display values
            let D_display = D_computed ? formatScientific(D_computed, 2) : '-';
            let k_display = k_computed ? k_computed.toFixed(3) : '-';

            // If auto is disabled and override is set, show override
            const D_value_display = (!D_auto && D_override) ? formatScientific(parseFloat(D_override), 2) : D_display;
            const k_value_display = (!k_auto && k_override) ? parseFloat(k_override).toFixed(3) : k_display;

            // Calculate partition coefficient KF/P = k/k0
            const k0_value = substance.k0_computed || substance.k0_override || 1.0;
            const k_current = (!k_auto && k_override) ? parseFloat(k_override) : k_computed;
            const K_FP = k_current && k0_value ? (k_current / k0_value).toFixed(3) : '-';
            const k_tooltip = k_auto
                ? `k = ${k_display} (auto)\\nk0 = ${k0_value.toFixed(3)}\\nKF/P = k/k0 = ${K_FP}`
                : `k = ${k_value_display} (manual)\\nk0 = ${k0_value.toFixed(3)}\\nKF/P = k/k0 = ${K_FP}`;

            const rowId = `${subId}_L${layerIndex}`;

            rows.push(`
                <tr class="${isAssigned ? 'bg-green-50 dark:bg-green-900/20' : ''}" data-row-id="${rowId}">
                    <td class="border border-gray-300 dark:border-gray-600 px-2 py-1.5 text-left">
                        <div class="flex items-center gap-2">
                            <input type="checkbox" class="rounded text-sfppy-green substance-layer-toggle"
                                   data-substance="${subId}" data-layer="${layerIndex}"
                                   ${isAssigned ? 'checked' : ''}
                                   onchange="toggleSubstanceLayer('${subId}', ${layerIndex}, this.checked)">
                            <span class="font-medium text-sm truncate max-w-32 text-gray-800 dark:text-gray-200" title="${subName}">${subName}</span>
                        </div>
                    </td>
                    <td class="border border-gray-300 dark:border-gray-600 px-2 py-1 text-center">
                        <span class="text-xs text-gray-600 dark:text-gray-400">L${layerIndex}</span>
                        <span class="text-xs font-medium text-gray-800 dark:text-gray-200">${layer.polymer}</span>
                    </td>
                    <td class="border border-gray-300 dark:border-gray-600 px-2 py-1 text-center">
                        <input type="number" class="w-20 text-center text-xs p-1 border rounded assembly-c0-input bg-white dark:bg-gray-800 text-gray-800 dark:text-gray-200 border-gray-300 dark:border-gray-600
                                      ${isAssigned ? '' : 'bg-gray-100 dark:bg-gray-700 text-gray-400 dark:text-gray-500'}"
                               data-substance="${subId}" data-layer="${layerIndex}"
                               value="${isAssigned ? c0Value : ''}"
                               placeholder="${isAssigned ? '1000' : '-'}"
                               ${isAssigned ? '' : 'disabled'}
                               min="0" step="any"
                               onchange="updateC0Value('${subId}', ${layerIndex}, this.value)">
                    </td>
                    <td class="border border-gray-300 dark:border-gray-600 px-2 py-1 text-center">
                        <select class="text-xs p-1 border rounded assembly-unit-select bg-white dark:bg-gray-700 text-gray-800 dark:text-gray-200 border-gray-300 dark:border-gray-600
                                      ${isAssigned ? '' : 'bg-gray-100 dark:bg-gray-700 text-gray-400 dark:text-gray-500'}"
                                data-substance="${subId}" data-layer="${layerIndex}"
                                ${isAssigned ? '' : 'disabled'}
                                onchange="updateC0Unit('${subId}', this.value)">
                            ${state.concentrationUnits.map(unit =>
                                `<option value="${unit}" ${unit === c0Unit ? 'selected' : ''}>${unit}</option>`
                            ).join('')}
                        </select>
                    </td>
                    <td class="border border-gray-300 dark:border-gray-600 px-2 py-1 text-center">
                        <div class="flex flex-col items-center gap-0.5">
                            <div class="flex items-center gap-1">
                                <input type="checkbox" class="w-3 h-3 rounded text-sfppy-green assembly-D-auto"
                                       data-substance="${subId}" data-layer="${layerIndex}"
                                       ${D_auto ? 'checked' : ''}
                                       ${isAssigned ? '' : 'disabled'}
                                       title="Auto-compute D"
                                       onchange="toggleDAutoCompute('${subId}', ${layerIndex}, this.checked)">
                                <input type="text" class="w-16 text-center text-[10px] p-0.5 border rounded font-mono assembly-D-input bg-white dark:bg-gray-800 text-gray-800 dark:text-gray-200 border-gray-300 dark:border-gray-600
                                              ${isAssigned ? '' : 'bg-gray-100 dark:bg-gray-700 text-gray-400 dark:text-gray-500'}
                                              ${!D_auto ? 'border-orange-400 bg-orange-50 dark:bg-orange-900/20' : ''}"
                                       data-substance="${subId}" data-layer="${layerIndex}"
                                       value="${D_auto ? '' : (D_override || '')}"
                                       placeholder="${D_display}"
                                       ${isAssigned && !D_auto ? '' : 'disabled'}
                                       title="${D_auto ? 'Auto: ' + D_display : 'Manual override'}"
                                       onchange="updateDOverride('${subId}', ${layerIndex}, this.value)">
                            </div>
                            <span class="text-[9px] text-gray-400">${D_auto ? 'auto' : 'manual'}</span>
                        </div>
                    </td>
                    <td class="border border-gray-300 dark:border-gray-600 px-2 py-1 text-center">
                        <div class="flex flex-col items-center gap-0.5">
                            <div class="flex items-center gap-1">
                                <input type="checkbox" class="w-3 h-3 rounded text-sfppy-green assembly-k-auto"
                                       data-substance="${subId}" data-layer="${layerIndex}"
                                       ${k_auto ? 'checked' : ''}
                                       ${isAssigned ? '' : 'disabled'}
                                       title="Auto-compute k"
                                       onchange="toggleKAutoCompute('${subId}', ${layerIndex}, this.checked)">
                                <input type="number" class="w-14 text-center text-xs p-0.5 border rounded assembly-k-input bg-white dark:bg-gray-800 text-gray-800 dark:text-gray-200 border-gray-300 dark:border-gray-600
                                              ${isAssigned ? '' : 'bg-gray-100 dark:bg-gray-700 text-gray-400 dark:text-gray-500'}
                                              ${!k_auto ? 'border-orange-400 bg-orange-50 dark:bg-orange-900/20' : ''}"
                                       data-substance="${subId}" data-layer="${layerIndex}"
                                       value="${k_auto ? '' : (k_override || '')}"
                                       placeholder="${k_display}"
                                       ${isAssigned && !k_auto ? '' : 'disabled'}
                                       min="0.001" step="0.01"
                                       title="${k_tooltip}"
                                       onchange="updateKOverride('${subId}', ${layerIndex}, this.value)">
                            </div>
                            <span class="text-[9px] text-gray-400">${k_auto ? 'auto' : 'manual'}</span>
                        </div>
                    </td>
                    <td class="border border-gray-300 dark:border-gray-600 px-2 py-1 text-center text-xs
                               ${smlValue ? 'text-green-600 dark:text-green-400 font-medium' : 'text-gray-400 dark:text-gray-500'}">
                        ${smlValue ? `${smlValue}` : '-'}
                    </td>
                </tr>
            `);
        }
    }

    tbody.innerHTML = rows.join('');

    // Update no-assignment warnings
    updateNoSubstanceAssignedWarnings();

    // Update C0 profile histogram
    renderC0ProfileHistogram();
}

// Toggle substance presence in a layer
async function toggleSubstanceLayer(substanceId, layerIndex, isChecked) {
    if (!state.c0Matrix[substanceId]) {
        state.c0Matrix[substanceId] = {};
    }

    if (isChecked) {
        // Assign with default C0 = 1000 mg/kg
        state.c0Matrix[substanceId][layerIndex] = {
            value: 1000,
            unit: state.c0Units[substanceId] || state.defaultC0Unit,
        };

        // Sync to backend (ignore errors - local state is authoritative)
        if (state.sessionId) {
            apiPost(`/sessions/${state.sessionId}/assignments`, {
                substance_id: substanceId,
                layer_index: layerIndex,
                C0: 1000,
                C0_unit: state.c0Units[substanceId] || state.defaultC0Unit,
            }).catch(() => {});  // Silently ignore - session may not exist
        }
    } else {
        // Remove assignment
        delete state.c0Matrix[substanceId][layerIndex];

        // Sync to backend (silently ignore - local state is authoritative)
        if (state.sessionId) {
            apiDelete(`/sessions/${state.sessionId}/assignments/${substanceId}/${layerIndex}`).catch(() => {});
        }
    }

    updateLayerVisualization();
    updateC0Matrix();
    await updateAssemblySubstances();
    updateNoSubstanceAssignedWarnings();
}

// Check if any substance is assigned to any layer with C0 > 0
function hasAnySubstanceAssigned() {
    for (const substanceId of Object.keys(state.c0Matrix)) {
        for (const layerIndex of Object.keys(state.c0Matrix[substanceId])) {
            const c0Data = state.c0Matrix[substanceId][layerIndex];
            const c0Value = typeof c0Data === 'object' ? c0Data.value : c0Data;
            if (c0Value > 0) {
                return true;
            }
        }
    }
    return false;
}

// Update warning displays in Assembly and Simulate tabs
function updateNoSubstanceAssignedWarnings() {
    const assemblyWarning = document.getElementById('assembly-no-assignment-warning');
    const simulateWarning = document.getElementById('simulate-no-assignment-warning');

    const hasSubstances = state.substances.length > 0;
    const hasAssignment = hasAnySubstanceAssigned();

    // Show warning if: substances exist but none are assigned
    const showWarning = hasSubstances && !hasAssignment;

    if (assemblyWarning) {
        assemblyWarning.classList.toggle('hidden', !showWarning);
    }
    if (simulateWarning) {
        simulateWarning.classList.toggle('hidden', !showWarning);
    }
}

// ========== C0 PROFILE HISTOGRAM ==========
// Track the current view mode
state.c0ProfileMode = 'layer'; // 'layer' or 'substance'

// Color palette for substances/layers (semi-transparent)
const HISTOGRAM_COLORS = [
    'rgba(59, 130, 246, 0.7)',   // blue
    'rgba(239, 68, 68, 0.7)',    // red
    'rgba(34, 197, 94, 0.7)',    // green
    'rgba(249, 115, 22, 0.7)',   // orange
    'rgba(168, 85, 247, 0.7)',   // purple
    'rgba(236, 72, 153, 0.7)',   // pink
    'rgba(20, 184, 166, 0.7)',   // teal
    'rgba(234, 179, 8, 0.7)',    // yellow
];

function renderC0ProfileHistogram() {
    const canvas = document.getElementById('c0-profile-canvas');
    if (!canvas) return;

    const ctx = canvas.getContext('2d');
    const width = canvas.width;
    const height = canvas.height;
    const padding = { top: 20, right: 20, bottom: 35, left: 50 };

    // Clear canvas
    ctx.clearRect(0, 0, width, height);

    // Check if we have data
    if (state.substances.length === 0 || state.layers.length === 0) {
        ctx.fillStyle = '#999';
        ctx.font = '12px sans-serif';
        ctx.textAlign = 'center';
        ctx.fillText('No substances assigned', width / 2, height / 2);
        return;
    }

    // Collect data based on mode
    const data = [];
    const labels = [];
    let maxC0 = 0;

    if (state.c0ProfileMode === 'layer') {
        // Group by layer, stacked bars for each substance
        state.layers.forEach(layer => {
            const layerData = [];
            state.substances.forEach((sub, subIdx) => {
                const subId = sub.id || sub.cid || sub.name;
                const c0 = getC0Value(subId, layer.index);
                layerData.push({ value: c0, color: HISTOGRAM_COLORS[subIdx % HISTOGRAM_COLORS.length], name: sub.name });
                if (c0 > maxC0) maxC0 = c0;
            });
            data.push(layerData);
            labels.push(`L${layer.index} (${layer.polymer})`);
        });
    } else {
        // Group by substance, stacked bars for each layer
        state.substances.forEach((sub, subIdx) => {
            const subId = sub.id || sub.cid || sub.name;
            const subData = [];
            state.layers.forEach((layer, layIdx) => {
                const c0 = getC0Value(subId, layer.index);
                subData.push({ value: c0, color: HISTOGRAM_COLORS[layIdx % HISTOGRAM_COLORS.length], name: `L${layer.index}` });
                if (c0 > maxC0) maxC0 = c0;
            });
            data.push(subData);
            // Truncate long substance names
            const name = sub.name.length > 12 ? sub.name.substring(0, 10) + '...' : sub.name;
            labels.push(name);
        });
    }

    if (maxC0 === 0) {
        ctx.fillStyle = '#999';
        ctx.font = '12px sans-serif';
        ctx.textAlign = 'center';
        ctx.fillText('All concentrations are zero', width / 2, height / 2);
        return;
    }

    // Calculate bar dimensions
    const chartWidth = width - padding.left - padding.right;
    const chartHeight = height - padding.top - padding.bottom;
    const barGroupWidth = chartWidth / data.length;
    const barWidth = barGroupWidth * 0.6;
    const barGap = (barGroupWidth - barWidth) / 2;

    // Draw Y axis
    ctx.strokeStyle = '#666';
    ctx.lineWidth = 1;
    ctx.beginPath();
    ctx.moveTo(padding.left, padding.top);
    ctx.lineTo(padding.left, height - padding.bottom);
    ctx.lineTo(width - padding.right, height - padding.bottom);
    ctx.stroke();

    // Draw Y axis labels
    ctx.fillStyle = '#666';
    ctx.font = '10px sans-serif';
    ctx.textAlign = 'right';
    ctx.textBaseline = 'middle';
    const ySteps = 4;
    for (let i = 0; i <= ySteps; i++) {
        const value = (maxC0 * i / ySteps);
        const y = height - padding.bottom - (chartHeight * i / ySteps);
        const label = value >= 1000 ? `${(value/1000).toFixed(1)}k` : value.toFixed(0);
        ctx.fillText(label, padding.left - 5, y);

        // Draw grid line
        ctx.strokeStyle = '#eee';
        ctx.beginPath();
        ctx.moveTo(padding.left, y);
        ctx.lineTo(width - padding.right, y);
        ctx.stroke();
    }

    // Draw Y axis title
    ctx.save();
    ctx.translate(12, height / 2);
    ctx.rotate(-Math.PI / 2);
    ctx.fillStyle = '#666';
    ctx.font = '10px sans-serif';
    ctx.textAlign = 'center';
    ctx.fillText('C₀ (mg/kg)', 0, 0);
    ctx.restore();

    // Draw bars (stacked)
    data.forEach((group, groupIdx) => {
        const x = padding.left + groupIdx * barGroupWidth + barGap;
        let yOffset = 0;

        group.forEach(bar => {
            if (bar.value > 0) {
                const barHeight = (bar.value / maxC0) * chartHeight;
                const y = height - padding.bottom - barHeight - yOffset;

                // Draw bar
                ctx.fillStyle = bar.color;
                ctx.fillRect(x, y, barWidth, barHeight);

                // Draw bar border
                ctx.strokeStyle = bar.color.replace('0.7', '1');
                ctx.lineWidth = 1;
                ctx.strokeRect(x, y, barWidth, barHeight);

                // Only stack in layer mode (show total per layer)
                if (state.c0ProfileMode === 'layer') {
                    yOffset += barHeight;
                }
            }
        });

        // Draw X label
        ctx.fillStyle = '#666';
        ctx.font = '9px sans-serif';
        ctx.textAlign = 'center';
        ctx.textBaseline = 'top';
        ctx.fillText(labels[groupIdx], x + barWidth / 2, height - padding.bottom + 5);
    });

    // Draw legend
    const legendItems = state.c0ProfileMode === 'layer' ? state.substances : state.layers;
    const legendX = padding.left + 5;
    let legendY = padding.top;

    ctx.font = '9px sans-serif';
    ctx.textAlign = 'left';
    legendItems.forEach((item, idx) => {
        const name = state.c0ProfileMode === 'layer'
            ? (item.name.length > 15 ? item.name.substring(0, 13) + '...' : item.name)
            : `L${item.index}`;
        ctx.fillStyle = HISTOGRAM_COLORS[idx % HISTOGRAM_COLORS.length];
        ctx.fillRect(legendX, legendY, 10, 10);
        ctx.fillStyle = '#666';
        ctx.fillText(name, legendX + 14, legendY + 8);
        legendY += 12;
    });
}

function initC0ProfileToggle() {
    const byLayerBtn = document.getElementById('c0-profile-by-layer');
    const bySubstanceBtn = document.getElementById('c0-profile-by-substance');

    if (byLayerBtn) {
        byLayerBtn.addEventListener('click', () => {
            state.c0ProfileMode = 'layer';
            byLayerBtn.classList.add('bg-sfppy-green', 'text-white');
            byLayerBtn.classList.remove('bg-white', 'dark:bg-gray-700', 'text-gray-700', 'dark:text-gray-300');
            bySubstanceBtn.classList.remove('bg-sfppy-green', 'text-white');
            bySubstanceBtn.classList.add('bg-white', 'dark:bg-gray-700', 'text-gray-700', 'dark:text-gray-300');
            renderC0ProfileHistogram();
        });
    }

    if (bySubstanceBtn) {
        bySubstanceBtn.addEventListener('click', () => {
            state.c0ProfileMode = 'substance';
            bySubstanceBtn.classList.add('bg-sfppy-green', 'text-white');
            bySubstanceBtn.classList.remove('bg-white', 'dark:bg-gray-700', 'text-gray-700', 'dark:text-gray-300');
            byLayerBtn.classList.remove('bg-sfppy-green', 'text-white');
            byLayerBtn.classList.add('bg-white', 'dark:bg-gray-700', 'text-gray-700', 'dark:text-gray-300');
            renderC0ProfileHistogram();
        });
    }
}

// Toggle D auto-compute for substance×layer
async function toggleDAutoCompute(substanceId, layerIndex, isAuto) {
    const sub = state.substances.find(s => (s.id || s.cid || s.name) === substanceId);
    if (!sub) return;

    if (isAuto) {
        // Enable auto-compute - clear override
        sub[`D_auto_L${layerIndex}`] = true;
        delete sub[`D_override_L${layerIndex}`];
    } else {
        // Disable auto-compute - keep current computed value as override
        sub[`D_auto_L${layerIndex}`] = false;
        const computed = sub[`D_computed_L${layerIndex}`];
        if (computed) {
            sub[`D_override_L${layerIndex}`] = computed.toExponential(2);
        }
    }

    await updateAssemblySubstances();
}

// Toggle k auto-compute for substance×layer
async function toggleKAutoCompute(substanceId, layerIndex, isAuto) {
    const sub = state.substances.find(s => (s.id || s.cid || s.name) === substanceId);
    if (!sub) return;

    if (isAuto) {
        // Enable auto-compute - clear override
        sub[`k_auto_L${layerIndex}`] = true;
        delete sub[`k_override_L${layerIndex}`];
    } else {
        // Disable auto-compute - keep current computed value as override
        sub[`k_auto_L${layerIndex}`] = false;
        const computed = sub[`k_computed_L${layerIndex}`];
        if (computed) {
            sub[`k_override_L${layerIndex}`] = computed;
        }
    }

    await updateAssemblySubstances();
}

// Update D override for substance×layer
function updateDOverride(substanceId, layerIndex, value) {
    const sub = state.substances.find(s => (s.id || s.cid || s.name) === substanceId);
    if (!sub) return;

    if (value && value.trim() !== '') {
        sub[`D_override_L${layerIndex}`] = value.trim();
    } else {
        delete sub[`D_override_L${layerIndex}`];
    }
}

// Update k override for substance×layer
function updateKOverride(substanceId, layerIndex, value) {
    const sub = state.substances.find(s => (s.id || s.cid || s.name) === substanceId);
    if (!sub) return;

    if (value && parseFloat(value) > 0) {
        sub[`k_override_L${layerIndex}`] = parseFloat(value);
    } else {
        delete sub[`k_override_L${layerIndex}`];
    }
}

// Called when polymer changes - recalculates D/k for substances with auto-compute enabled
async function onPolymerChange(layerIndex) {
    // Refresh the assembly substances matrix to recalculate D/k
    await updateAssemblySubstances();
}

function getMaterialOptionsHTML(selectedValue) {
    const groups = {
        'Polyethylenes (Rubbery)': [
            ['LDPE', 'LDPE - Low Density PE'],
            ['LLDPE', 'LLDPE - Linear Low Density PE'],
            ['HDPE', 'HDPE - High Density PE'],
        ],
        'Polypropylenes': [
            ['PP', 'PP - Polypropylene'],
            ['PPrubber', 'PP Rubber'],
            ['oPP', 'oPP - Oriented PP'],
        ],
        'Polyesters (Barrier)': [
            ['gPET', 'gPET - Glassy PET (barrier)'],
            ['rPET', 'rPET - Rubbery PET'],
            ['wPET', 'wPET - Plasticized PET'],
            ['PBT', 'PBT - Polybutylene Terephthalate'],
            ['PEN', 'PEN - Polyethylene Naphthalate'],
        ],
        'Polyamides (Barrier)': [
            ['PA6', 'PA6 - Nylon 6'],
            ['PA66', 'PA66 - Nylon 6,6'],
        ],
        'Styrenic Polymers': [
            ['PS', 'PS - Polystyrene'],
            ['HIPS', 'HIPS - High Impact PS'],
            ['SBS', 'SBS - Styrene-Butadiene-Styrene'],
        ],
        'Other Polymers': [
            ['PMMA', 'PMMA - Acrylic'],
            ['PVAc', 'PVAc - Polyvinyl Acetate'],
            ['rigidPVC', 'Rigid PVC'],
            ['plasticizedPVC', 'Plasticized PVC'],
        ],
        'Adhesives': [
            ['AdhesiveEVA', 'Adhesive EVA'],
            ['AdhesiveVAE', 'Adhesive VAE'],
            ['AdhesivePVAC', 'Adhesive PVAC'],
            ['AdhesiveAcrylate', 'Adhesive Acrylate'],
            ['AdhesivePU', 'Adhesive PU'],
            ['AdhesiveNaturalRubber', 'Adhesive Natural Rubber'],
            ['AdhesiveSyntheticRubber', 'Adhesive Synthetic Rubber'],
        ],
        'Paper & Board': [
            ['Paper', 'Paper'],
            ['Cardboard', 'Cardboard'],
        ],
        'Gases': [
            ['air', 'Air (gap)'],
        ],
    };

    return Object.entries(groups).map(([groupName, options]) => `
        <optgroup label="${groupName}">
            ${options.map(([value, label]) =>
                `<option value="${value}" ${value === selectedValue ? 'selected' : ''}>${label}</option>`
            ).join('')}
        </optgroup>
    `).join('');
}

function initLayerTableListeners(index) {
    const tbody = document.getElementById('layer-table-body');
    if (!tbody) return;

    const row = tbody.querySelector(`tr[data-layer="${index}"]`);
    if (!row) return;

    // Polymer change
    const polymerSelect = row.querySelector('.layer-polymer');
    if (polymerSelect) {
        polymerSelect.addEventListener('change', async (e) => {
            state.layers[index - 1].polymer = e.target.value;
            state.layers[index - 1].material = e.target.value;  // Also set material for consistency
            updateLayerVisualization();
            updateTgDisplay(index);
            // Update C0 matrix to reflect new polymer name
            updateC0Matrix();
            // Recalculate D/k for substances with auto-compute enabled
            await onPolymerChange(index);
        });
    }

    // Thickness change
    const thicknessInput = row.querySelector('.layer-thickness');
    if (thicknessInput) {
        thicknessInput.addEventListener('input', (e) => {
            state.layers[index - 1].thickness = parseFloat(e.target.value) || 0;
            updateLayerVisualization();
        });
    }

    const thicknessUnitSelect = row.querySelector('.layer-thickness-unit');
    if (thicknessUnitSelect) {
        thicknessUnitSelect.addEventListener('change', (e) => {
            state.layers[index - 1].thickness_unit = e.target.value;
        });
    }
}

function updateTgDisplay(index) {
    const polymer = state.layers[index - 1].polymer;
    const tg = MATERIAL_TG[polymer];
    const tgDisplay = tg !== null && tg !== undefined ? `${tg}°C` : '-';

    const tgCell = document.querySelector(`.layer-tg-display[data-layer="${index}"]`);
    if (tgCell) {
        tgCell.textContent = tgDisplay;
    }
}

function createLayerCard(index) {
    return `
        <div class="layer-card" data-layer="${index}">
            <div class="flex items-center justify-between mb-3">
                <h3 class="font-medium text-gray-800">
                    <span class="text-blue-600">Layer ${index}</span>
                    ${index === 1 ? '<span class="text-xs text-gray-500 ml-2">(Food Contact)</span>' : ''}
                </h3>
                <div class="flex items-center space-x-2">
                    <span class="override-badge computed" id="layer-${index}-override-badge">Computed</span>
                </div>
            </div>
            <div class="grid grid-cols-1 md:grid-cols-3 gap-4">
                <div>
                    <label class="block text-sm text-gray-600 mb-1">Polymer</label>
                    <select class="w-full border rounded px-3 py-2 layer-polymer" data-layer="${index}">
                        <option value="LDPE">LDPE - Low Density PE</option>
                        <option value="HDPE">HDPE - High Density PE</option>
                        <option value="PP">PP - Polypropylene</option>
                        <option value="PET">PET - Polyethylene Terephthalate</option>
                        <option value="gPET">gPET - Glassy PET</option>
                        <option value="PS">PS - Polystyrene</option>
                        <option value="PA6">PA6 - Nylon 6</option>
                        <option value="EVOH">EVOH - Ethylene Vinyl Alcohol</option>
                    </select>
                </div>
                <div>
                    <label class="block text-sm text-gray-600 mb-1">Thickness</label>
                    <div class="flex">
                        <input type="number" class="flex-1 border rounded-l px-3 py-2 layer-thickness"
                               data-layer="${index}" value="100" min="0.1" step="1">
                        <select class="border-t border-r border-b rounded-r px-2 py-2 bg-gray-50 layer-thickness-unit" data-layer="${index}">
                            <option value="um">µm</option>
                            <option value="mm">mm</option>
                            <option value="nm">nm</option>
                        </select>
                    </div>
                </div>
                <div>
                    <label class="block text-sm text-gray-600 mb-1">Temperature</label>
                    <div class="flex">
                        <input type="number" class="flex-1 border rounded-l px-3 py-2 layer-temp"
                               data-layer="${index}" value="25" step="1">
                        <span class="border-t border-r border-b rounded-r px-3 py-2 bg-gray-50">°C</span>
                    </div>
                </div>
            </div>
            <details class="mt-4 border-t pt-4">
                <summary class="cursor-pointer text-sm font-medium text-gray-700 hover:text-blue-600">
                    Parameter Overrides (D, k, k0)
                </summary>
                <div class="mt-3 grid grid-cols-1 md:grid-cols-3 gap-4">
                    <div>
                        <label class="flex items-center space-x-2 mb-1">
                            <input type="checkbox" class="rounded layer-D-use-computed" data-layer="${index}" checked>
                            <span class="text-sm text-gray-600">Use computed D</span>
                        </label>
                        <input type="text" class="w-full border rounded px-3 py-2 text-sm layer-D-override"
                               data-layer="${index}" placeholder="e.g., 1.2e-14" disabled>
                        <span class="text-xs text-gray-500">m²/s</span>
                    </div>
                    <div>
                        <label class="flex items-center space-x-2 mb-1">
                            <input type="checkbox" class="rounded layer-k-use-computed" data-layer="${index}" checked>
                            <span class="text-sm text-gray-600">Use computed k</span>
                        </label>
                        <input type="text" class="w-full border rounded px-3 py-2 text-sm layer-k-override"
                               data-layer="${index}" placeholder="e.g., 1.0" disabled>
                        <span class="text-xs text-gray-500">dimensionless</span>
                    </div>
                    <div>
                        <label class="flex items-center space-x-2 mb-1">
                            <input type="checkbox" class="rounded layer-k0-use-computed" data-layer="${index}" checked>
                            <span class="text-sm text-gray-600">Use computed k0</span>
                        </label>
                        <input type="text" class="w-full border rounded px-3 py-2 text-sm layer-k0-override"
                               data-layer="${index}" placeholder="e.g., 1.0" disabled>
                        <span class="text-xs text-gray-500">dimensionless</span>
                    </div>
                </div>
            </details>
        </div>
    `;
}

function initLayerListeners(index) {
    const card = document.querySelector(`.layer-card[data-layer="${index}"]`);
    if (!card) return;

    // Polymer change
    card.querySelector('.layer-polymer').addEventListener('change', async (e) => {
        state.layers[index - 1].polymer = e.target.value;
        updateLayerVisualization();
        // Recalculate D/k for substances with auto-compute enabled
        await onPolymerChange(index);
    });

    // Thickness change
    card.querySelector('.layer-thickness').addEventListener('input', (e) => {
        state.layers[index - 1].thickness = parseFloat(e.target.value) || 0;
        updateLayerVisualization();
    });

    card.querySelector('.layer-thickness-unit').addEventListener('change', (e) => {
        state.layers[index - 1].thickness_unit = e.target.value;
    });

    // Temperature change
    card.querySelector('.layer-temp').addEventListener('input', (e) => {
        state.layers[index - 1].temperature_C = parseFloat(e.target.value) || 25;
    });

    // C0 (initial concentration) change
    const c0Input = card.querySelector('.layer-C0');
    if (c0Input) {
        c0Input.addEventListener('input', (e) => {
            state.layers[index - 1].C0 = parseFloat(e.target.value) || 0;
            updateLayerVisualization();  // Update concentration histogram
        });
    }

    const c0UnitSelect = card.querySelector('.layer-C0-unit');
    if (c0UnitSelect) {
        c0UnitSelect.addEventListener('change', (e) => {
            state.layers[index - 1].C0_unit = e.target.value;
        });
    }

    // Override checkboxes
    ['D', 'k', 'k0'].forEach(param => {
        const checkbox = card.querySelector(`.layer-${param}-use-computed`);
        const input = card.querySelector(`.layer-${param}-override`);

        if (checkbox && input) {
            checkbox.addEventListener('change', () => {
                input.disabled = checkbox.checked;
                updateOverrideBadge(index);
            });
        }
    });
}

function updateOverrideBadge(index) {
    const card = document.querySelector(`.layer-card[data-layer="${index}"]`);
    const badge = document.getElementById(`layer-${index}-override-badge`);
    if (!card || !badge) return;

    const hasOverride = ['D', 'k', 'k0'].some(param => {
        const checkbox = card.querySelector(`.layer-${param}-use-computed`);
        return checkbox && !checkbox.checked;
    });

    if (hasOverride) {
        badge.textContent = 'Override';
        badge.classList.remove('computed');
        badge.classList.add('override');
    } else {
        badge.textContent = 'Computed';
        badge.classList.remove('override');
        badge.classList.add('computed');
    }
}

function updateLayerVisualization() {
    const stack = document.getElementById('layer-stack');
    const concBar = document.getElementById('concentration-bar');

    // Calculate total thickness for proportions
    const totalThickness = state.layers.reduce((sum, l) => {
        return sum + convertThickness(l.thickness, l.thickness_unit);
    }, 0);

    // Find max C0 across all layers (from C0 matrix)
    let maxC0 = 0.01;
    state.layers.forEach(layer => {
        const layerC0 = getTotalC0ForLayer(layer.index);
        if (layerC0 > maxC0) maxC0 = layerC0;
    });

    // Minimum visual width for thin layers (as % of total)
    const MIN_WIDTH_PERCENT = 8;

    // Build layer stack HTML
    let stackHtml = `
        <div class="food-block">
            <span class="text-lg">🍏</span>
            <span class="text-xs">Food</span>
        </div>
    `;

    // Build concentration bar HTML
    let concHtml = `<div class="conc-block" style="flex-grow: 1; min-width: 50px;"></div>`;

    state.layers.forEach(layer => {
        const thickness = convertThickness(layer.thickness, layer.thickness_unit);
        let widthPercent = totalThickness > 0 ? (thickness / totalThickness * 100) : 100;

        // Ensure minimum width for very thin layers
        widthPercent = Math.max(widthPercent, MIN_WIDTH_PERCENT);

        const polymerClass = `polymer-${layer.polymer}` || 'polymer-default';
        const isSelected = layer.index === state.selectedLayerIndex;
        const tg = MATERIAL_TG[layer.polymer];
        const tgLabel = tg !== null && tg !== undefined ? `${tg}°C` : '';

        // Calculate concentration bar height from C0 matrix (sum of all substances in layer)
        const totalC0 = getTotalC0ForLayer(layer.index);
        const concHeight = maxC0 > 0 ? (totalC0 / maxC0 * 100) : 0;
        const c0Display = totalC0 > 0 ? (totalC0 >= 1000 ? `${(totalC0/1000).toFixed(1)}k` : totalC0.toFixed(0)) : '0';

        stackHtml += `
            <div class="layer-block ${polymerClass} ${isSelected ? 'selected' : ''}"
                 style="flex-grow: ${widthPercent};"
                 data-layer="${layer.index}"
                 onclick="selectLayer(${layer.index})">
                <span>L${layer.index}</span>
                <span style="font-size: 0.6rem;">${layer.polymer}</span>
                <span class="thickness-label">${formatThickness(layer.thickness, layer.thickness_unit)}</span>
                ${tgLabel ? `<span class="tg-indicator">${tgLabel}</span>` : ''}
            </div>
        `;

        // Concentration bar for this layer
        concHtml += `
            <div class="conc-block" style="flex-grow: ${widthPercent};" data-layer="${layer.index}">
                <div class="conc-fill" style="height: ${concHeight}%;"></div>
                ${totalC0 > 0 ? `<span class="conc-label">${c0Display}</span>` : ''}
            </div>
        `;
    });

    stack.innerHTML = stackHtml;
    if (concBar) {
        concBar.innerHTML = concHtml;
    }

    // Update layer card visibility - show only selected layer card
    highlightSelectedLayerCard();
}

function getTotalC0ForLayer(layerIndex) {
    let total = 0;
    Object.values(state.c0Matrix).forEach(layerC0s => {
        if (layerC0s[layerIndex]) {
            // C0 can be stored as object {value, unit} or direct number
            const c0Data = layerC0s[layerIndex];
            const c0Value = typeof c0Data === 'object' ? c0Data.value : c0Data;
            if (typeof c0Value === 'number' && !isNaN(c0Value)) {
                total += c0Value;
            }
        }
    });
    return total;
}

function formatThickness(value, unit) {
    if (value >= 1000) return `${(value / 1000).toFixed(1)}k ${unit}`;
    if (value < 1) return `${value.toFixed(3)} ${unit}`;
    return `${value.toFixed(0)} ${unit}`;
}

function selectLayer(index) {
    state.selectedLayerIndex = index;
    updateLayerVisualization();
    updateLayerTable();
    updateLayerDetailPanel(index);
}

function updateLayerDetailPanel(index) {
    const detailNum = document.getElementById('detail-layer-num');
    const detailContact = document.getElementById('detail-layer-contact');

    if (detailNum) detailNum.textContent = index;
    if (detailContact) {
        detailContact.textContent = index === 1 ? '(Food Contact)' : '';
    }

    // Update override checkboxes data-layer attributes
    const panel = document.getElementById('layer-detail-panel');
    if (panel) {
        panel.querySelectorAll('[data-layer]').forEach(el => {
            el.setAttribute('data-layer', index);
        });
    }

    // Update density display for this layer
    updateLayerDensityDisplay(index);
}

function highlightSelectedLayerCard() {
    // Update table row selection
    document.querySelectorAll('#layer-table-body tr').forEach(row => {
        const layerIdx = parseInt(row.dataset.layer);
        if (layerIdx === state.selectedLayerIndex) {
            row.classList.add('selected');
        } else {
            row.classList.remove('selected');
        }
    });
}

// ========== C0 MATRIX FUNCTIONS ==========
function updateC0Matrix() {
    const container = document.getElementById('c0-matrix-container');
    const table = document.getElementById('c0-matrix');
    if (!container || !table) return;

    if (state.substances.length === 0) {
        container.classList.add('hidden');
        return;
    }

    container.classList.remove('hidden');

    // Build header with layer columns
    const thead = table.querySelector('thead tr');
    thead.innerHTML = '<th>Substance</th><th>Unit</th>' +
        state.layers.map(l => `<th>L${l.index} (${l.polymer})</th>`).join('') +
        '<th>SML</th>';

    // Build body with substance rows
    const tbody = table.querySelector('tbody');
    tbody.innerHTML = state.substances.map(sub => {
        const subId = sub.id || sub.cid || sub.name;
        const currentUnit = state.c0Units[subId] || state.defaultC0Unit;
        return `
            <tr>
                <td class="font-medium text-left text-gray-800 dark:text-gray-200">${sub.name}</td>
                <td>
                    <select class="c0-unit-select text-xs p-1 border rounded bg-white dark:bg-gray-700 text-gray-800 dark:text-gray-200 border-gray-300 dark:border-gray-600"
                            data-substance="${subId}"
                            onchange="updateC0Unit('${subId}', this.value)">
                        ${state.concentrationUnits.map(unit =>
                            `<option value="${unit}" ${unit === currentUnit ? 'selected' : ''}>${unit}</option>`
                        ).join('')}
                    </select>
                </td>
                ${state.layers.map(layer => {
                    const c0Value = getC0Value(subId, layer.index);
                    return `
                        <td>
                            <input type="number" class="c0-matrix-input"
                                   data-substance="${subId}" data-layer="${layer.index}"
                                   value="${c0Value}" min="0" step="any"
                                   onchange="updateC0Value('${subId}', ${layer.index}, this.value)">
                        </td>
                    `;
                }).join('')}
                <td class="text-xs ${sub.SML ? 'text-green-600 dark:text-green-400' : 'text-gray-400 dark:text-gray-500'}">${sub.SML || '-'}</td>
            </tr>
        `;
    }).join('');

    // Update no-assignment warnings
    updateNoSubstanceAssignedWarnings();
}

function getC0Value(substanceId, layerIndex) {
    if (!state.c0Matrix[substanceId]) return 0;
    const c0Data = state.c0Matrix[substanceId][layerIndex];
    if (!c0Data) return 0;
    // Handle both old format (number) and new format ({ value, unit })
    return typeof c0Data === 'object' ? c0Data.value : c0Data;
}

function getC0Unit(substanceId, layerIndex) {
    if (!state.c0Matrix[substanceId]) return state.defaultC0Unit;
    const c0Data = state.c0Matrix[substanceId][layerIndex];
    if (!c0Data) return state.c0Units[substanceId] || state.defaultC0Unit;
    return typeof c0Data === 'object' ? c0Data.unit : state.defaultC0Unit;
}

async function updateC0Value(substanceId, layerIndex, value) {
    if (!state.c0Matrix[substanceId]) {
        state.c0Matrix[substanceId] = {};
    }

    const numValue = parseFloat(value) || 0;
    const unit = state.c0Units[substanceId] || state.defaultC0Unit;

    state.c0Matrix[substanceId][layerIndex] = {
        value: numValue,
        unit: unit,
    };

    // Sync to backend (silently ignore - local state is authoritative)
    if (state.sessionId && numValue > 0) {
        apiPost(`/sessions/${state.sessionId}/assignments`, {
            substance_id: substanceId,
            layer_index: layerIndex,
            C0: numValue,
            C0_unit: unit,
        }).catch(() => {});
    }

    updateLayerVisualization();
    renderC0ProfileHistogram();
    // Don't call updateAssemblySubstances here to avoid recursion - it's called from where C0 changes originate
}

async function updateC0Unit(substanceId, unit) {
    state.c0Units[substanceId] = unit;

    // Update all assignments for this substance
    if (state.c0Matrix[substanceId]) {
        for (const layerIndex of Object.keys(state.c0Matrix[substanceId])) {
            const c0Data = state.c0Matrix[substanceId][layerIndex];
            if (typeof c0Data === 'object') {
                c0Data.unit = unit;
            } else {
                state.c0Matrix[substanceId][layerIndex] = { value: c0Data, unit };
            }

            // Sync to backend (silently ignore - local state is authoritative)
            const numValue = getC0Value(substanceId, parseInt(layerIndex));
            if (state.sessionId && numValue > 0) {
                apiPost(`/sessions/${state.sessionId}/assignments`, {
                    substance_id: substanceId,
                    layer_index: parseInt(layerIndex),
                    C0: numValue,
                    C0_unit: unit,
                }).catch(() => {});
            }
        }
    }
}

function convertThickness(value, unit) {
    const factors = { nm: 1e-9, um: 1e-6, mm: 1e-3, cm: 1e-2, m: 1 };
    return value * (factors[unit] || 1e-6);
}

async function validateAssembly() {
    const assemblyData = {
        name: 'Current Assembly',
        layers: state.layers.map(l => ({
            index: l.index,
            polymer: l.polymer,
            thickness: l.thickness,
            thickness_unit: l.thickness_unit,
            temperature_C: l.temperature_C,
        })),
    };

    try {
        const result = await apiPost('/assembly/validate', assemblyData);
        if (result.valid) {
            alert(`Assembly valid! Total thickness: ${result.total_thickness_um.toFixed(1)} µm`);
        } else {
            alert(`Validation errors:\n${result.errors.join('\n')}`);
        }
    } catch (e) {
        console.error('Validation failed:', e);
        alert('Validation failed. Check console for details.');
    }
}

function resetAssembly() {
    setLayerCount(1);
    document.querySelector('.layer-count-btn[data-count="1"]').click();
}

// ========== FOOD TAB ==========

// Food type definitions with icon, simulant mapping, and description
const FOOD_TYPES = [
    { code: 'fatty_oil', icon: '🫒', name: 'Fatty (oil)', simulant: 'olive_oil', category: 'fatty', desc: 'Olive oil simulant' },
    { code: 'fatty_eth', icon: '🧈', name: 'Fatty (ethanol)', simulant: 'ethanol', category: 'fatty', desc: 'Ethanol 95% simulant' },
    { code: 'aqueous', icon: '💧', name: 'Aqueous', simulant: 'water', category: 'aqueous', desc: 'Water simulant' },
    { code: 'acidic', icon: '🍋', name: 'Acidic', simulant: 'acetic_acid', category: 'acidic', desc: 'Acetic acid 3%' },
    { code: 'alcohol_10', icon: '🍺', name: 'Alcohol 10%', simulant: 'ethanol_10', category: 'alcoholic', desc: 'Beer, wine' },
    { code: 'alcohol_20', icon: '🍷', name: 'Alcohol 20%', simulant: 'ethanol_20', category: 'alcoholic', desc: 'Fortified wine' },
    { code: 'alcohol_50', icon: '🥃', name: 'Alcohol 50%', simulant: 'ethanol_50', category: 'alcoholic', desc: 'Spirits' },
    { code: 'dry', icon: '🍪', name: 'Dry', simulant: 'tenax', category: 'dry', desc: 'Tenax simulant' },
];

// Food texture definitions
const FOOD_TEXTURES = [
    { code: 'liquid', icon: '💧', name: 'Liquid', desc: 'Beverages, oils' },
    { code: 'semisolid', icon: '🍮', name: 'Semi-solid', desc: 'Yogurt, sauces' },
    { code: 'solid', icon: '🧱', name: 'Solid', desc: 'Cheese, meat' },
];

async function loadFoodData() {
    // Create food type + simulant icons
    const foodTypeContainer = document.getElementById('food-type-icons');
    if (foodTypeContainer) {
        foodTypeContainer.innerHTML = FOOD_TYPES.map(ft => `
            <button class="food-type-btn flex flex-col items-center p-2 border-2 rounded-lg
                           bg-white dark:bg-gray-700 hover:border-sfppy-green transition-colors min-w-[70px]
                           ${state.food.simulant === ft.simulant ? 'border-sfppy-green bg-green-50 dark:bg-green-900/30' : 'border-gray-300 dark:border-gray-600'}"
                    data-simulant="${ft.simulant}" data-category="${ft.category}" title="${ft.desc}">
                <span class="text-2xl">${ft.icon}</span>
                <span class="text-[10px] text-gray-600 dark:text-gray-400 text-center leading-tight mt-1">${ft.name}</span>
            </button>
        `).join('');

        foodTypeContainer.querySelectorAll('.food-type-btn').forEach(btn => {
            btn.addEventListener('click', () => {
                state.food.simulant = btn.dataset.simulant;
                state.food.category = btn.dataset.category;

                // Update visual selection
                foodTypeContainer.querySelectorAll('.food-type-btn').forEach(b => {
                    b.classList.remove('border-sfppy-green', 'bg-green-50', 'dark:bg-green-900/30');
                    b.classList.add('border-gray-300', 'dark:border-gray-600');
                });
                btn.classList.remove('border-gray-300', 'dark:border-gray-600');
                btn.classList.add('border-sfppy-green', 'bg-green-50', 'dark:bg-green-900/30');

                // Recalculate k0
                updateFoodK0Matrix();

                // Update food density display
                updateFoodDensityDisplay();
            });
        });
    }

    // Create texture icons
    const textureContainer = document.getElementById('food-texture-icons');
    if (textureContainer) {
        textureContainer.innerHTML = FOOD_TEXTURES.map(tx => `
            <button class="food-texture-btn flex flex-col items-center p-2 border-2 rounded-lg
                           bg-white dark:bg-gray-700 hover:border-sfppy-green transition-colors min-w-[70px]
                           ${state.food.matrixType === tx.code ? 'border-sfppy-green bg-green-50 dark:bg-green-900/30' : 'border-gray-300 dark:border-gray-600'}"
                    data-texture="${tx.code}" title="${tx.desc}">
                <span class="text-2xl">${tx.icon}</span>
                <span class="text-[10px] text-gray-600 dark:text-gray-400 text-center leading-tight mt-1">${tx.name}</span>
            </button>
        `).join('');

        textureContainer.querySelectorAll('.food-texture-btn').forEach(btn => {
            btn.addEventListener('click', () => {
                state.food.matrixType = btn.dataset.texture;

                // Update visual selection
                textureContainer.querySelectorAll('.food-texture-btn').forEach(b => {
                    b.classList.remove('border-sfppy-green', 'bg-green-50', 'dark:bg-green-900/30');
                    b.classList.add('border-gray-300', 'dark:border-gray-600');
                });
                btn.classList.remove('border-gray-300', 'dark:border-gray-600');
                btn.classList.add('border-sfppy-green', 'bg-green-50', 'dark:bg-green-900/30');
            });
        });
    }

    // Load condition presets
    try {
        const presetData = await apiGet('/food/conditions/presets');
        const container = document.getElementById('condition-presets');
        container.innerHTML = presetData.presets.map(p => `
            <button class="px-3 py-2 border border-gray-300 dark:border-gray-600 rounded text-sm
                           text-gray-700 dark:text-gray-200 bg-white dark:bg-gray-700
                           hover:bg-gray-50 dark:hover:bg-gray-600 hover:border-sfppy-green
                           transition-colors" data-preset="${p.code}">
                ${p.name}
            </button>
        `).join('');

        container.querySelectorAll('button').forEach(btn => {
            btn.addEventListener('click', () => {
                applyConditionPreset(btn.dataset.preset, presetData.presets);
            });
        });
    } catch (e) {
        console.error('Failed to load presets:', e);
    }

    // Init geometry
    initGeometry();

    // Update k0 matrix
    updateFoodK0Matrix();

    // Update food density display
    updateFoodDensityDisplay();
}

// Simulant polarity index for k0 computation (lower = more polar)
// Based on patankar.food simulant classes
const SIMULANT_POLARITY = {
    'water': 0.0,           // Most polar
    'acetic_acid': 0.1,     // Slightly less polar
    'ethanol_10': 0.3,
    'ethanol_20': 0.5,
    'ethanol_50': 0.7,
    'ethanol': 0.9,         // Ethanol 95%
    'olive_oil': 1.0,       // Non-polar (fatty)
    'tenax': 0.5,           // Dry (intermediate)
};

// Simulant display names
const SIMULANT_NAMES = {
    'water': 'Water',
    'acetic_acid': 'Acetic acid 3%',
    'ethanol_10': 'Ethanol 10%',
    'ethanol_20': 'Ethanol 20%',
    'ethanol_50': 'Ethanol 50%',
    'ethanol': 'Ethanol 95%',
    'olive_oil': 'Olive oil',
    'tenax': 'Tenax (dry)',
};

// Compute k0 from logP and simulant
// k0 is the partition coefficient between polymer and food
// For lipophilic substances (high logP), k0 is lower in fatty simulants (substance migrates more)
// For hydrophilic substances (low logP), k0 is lower in aqueous simulants
function computeK0(logP, simulant) {
    if (logP === null || logP === undefined) {
        return 1.0; // Default if logP unknown
    }

    const polarity = SIMULANT_POLARITY[simulant] || 0.5;

    // Simple model: k0 relates substance lipophilicity to simulant polarity
    // High logP + fatty simulant (polarity~1) → low k0 (migrates easily)
    // High logP + aqueous simulant (polarity~0) → high k0 (stays in polymer)
    // Low logP + aqueous simulant → low k0 (migrates easily)
    // Low logP + fatty simulant → high k0 (stays in polymer)

    // Empirical formula based on Piringer approach
    const logP_float = parseFloat(logP);

    // For fatty simulants (polarity > 0.7): lipophilic substances migrate more
    // For aqueous simulants (polarity < 0.3): hydrophilic substances migrate more
    let k0;
    if (polarity > 0.7) {
        // Fatty simulant: k0 decreases with logP (more migration for lipophilic)
        k0 = Math.pow(10, -0.5 * logP_float + 2);
    } else if (polarity < 0.3) {
        // Aqueous simulant: k0 decreases with -logP (more migration for hydrophilic)
        k0 = Math.pow(10, 0.5 * logP_float + 1);
    } else {
        // Intermediate: balanced
        k0 = Math.pow(10, Math.abs(logP_float - polarity * 5) * 0.2 + 0.5);
    }

    // Clamp to reasonable range
    return Math.max(0.01, Math.min(1000, k0));
}

// Update the Food tab k0 matrix
function updateFoodK0Matrix() {
    const noSubstancesDiv = document.getElementById('food-no-substances');
    const matrixDiv = document.getElementById('food-k0-matrix');
    const tbody = document.getElementById('food-k0-body');

    if (!noSubstancesDiv || !matrixDiv || !tbody) return;

    if (state.substances.length === 0) {
        noSubstancesDiv.classList.remove('hidden');
        matrixDiv.classList.add('hidden');
        return;
    }

    noSubstancesDiv.classList.add('hidden');
    matrixDiv.classList.remove('hidden');

    // Build rows for each substance
    const rows = [];
    const simulant = state.food.simulant || 'ethanol';
    const simulantName = SIMULANT_NAMES[simulant] || simulant;

    for (const substance of state.substances) {
        const subId = substance.id || substance.cid || substance.name;
        const subName = substance.name || subId;
        const logP = substance.logP || null;

        // Compute k0 from simulant and logP
        const k0_computed = computeK0(logP, simulant);
        substance.k0_computed = k0_computed; // Store for reference

        // Check auto-compute flag (default: true)
        const k0_auto = substance.k0_auto !== false;
        const k0_override = substance.k0_override || null;
        const k0_value = (k0_auto || !k0_override) ? k0_computed : k0_override;

        // Get k from first assigned layer (or default)
        let k_value = substance.k_computed || 1.0;

        // Compute K_FP = k / k0
        const K_FP = k_value / k0_value;

        // Get CF0 (initial concentration in food)
        const cf0Data = state.CF0[subId] || { value: 0, unit: 'mg/kg' };
        const cf0Value = cf0Data.value || 0;

        rows.push(`
            <tr data-substance="${subId}">
                <td class="border border-gray-300 dark:border-gray-600 px-2 py-1.5 text-left">
                    <span class="font-medium text-sm">${subName}</span>
                </td>
                <td class="border border-gray-300 dark:border-gray-600 px-2 py-1 text-center text-xs">
                    ${logP !== null ? parseFloat(logP).toFixed(2) : '-'}
                </td>
                <td class="border border-gray-300 dark:border-gray-600 px-2 py-1 text-center">
                    <div class="flex flex-col items-center gap-0.5">
                        <div class="flex items-center gap-1">
                            <input type="checkbox" class="w-3 h-3 rounded text-sfppy-green food-k0-auto"
                                   data-substance="${subId}"
                                   ${k0_auto ? 'checked' : ''}
                                   title="Auto-compute k0"
                                   onchange="toggleK0AutoCompute('${subId}', this.checked)">
                            <input type="number" class="w-14 text-center text-xs p-0.5 border rounded food-k0-input
                                          ${!k0_auto ? 'border-orange-400 bg-orange-50 dark:bg-orange-900/20' : ''}"
                                   data-substance="${subId}"
                                   value="${k0_auto ? '' : (k0_override || '')}"
                                   placeholder="${k0_computed.toFixed(2)}"
                                   ${k0_auto ? 'disabled' : ''}
                                   min="0.001" step="0.01"
                                   title="${k0_auto ? 'Auto: ' + k0_computed.toFixed(2) : 'Manual override'}"
                                   onchange="updateK0Override('${subId}', this.value)">
                        </div>
                        <span class="text-[9px] text-gray-400">${k0_auto ? 'auto' : 'manual'}</span>
                    </div>
                </td>
                <td class="border border-gray-300 dark:border-gray-600 px-2 py-1 text-center text-xs font-medium text-sfppy-green">
                    ${K_FP.toFixed(3)}
                </td>
                <td class="border border-gray-300 dark:border-gray-600 px-2 py-1 text-center">
                    <input type="number" class="w-16 text-center text-xs p-0.5 border rounded food-cf0-input
                                  ${cf0Value > 0 ? 'border-blue-400 bg-blue-50 dark:bg-blue-900/20' : 'border-gray-300 dark:border-gray-600'}"
                           data-substance="${subId}"
                           value="${cf0Value || ''}"
                           placeholder="0"
                           min="0" step="any"
                           title="Initial concentration in food (mg/kg)"
                           onchange="updateCF0('${subId}', this.value)">
                </td>
                <td class="border border-gray-300 dark:border-gray-600 px-2 py-1 text-center text-xs text-gray-600 dark:text-gray-400">
                    ${simulantName}
                </td>
            </tr>
        `);
    }

    tbody.innerHTML = rows.join('');
}

// Toggle k0 auto-compute for substance
function toggleK0AutoCompute(substanceId, isAuto) {
    const sub = state.substances.find(s => (s.id || s.cid || s.name) === substanceId);
    if (!sub) return;

    if (isAuto) {
        sub.k0_auto = true;
        delete sub.k0_override;
    } else {
        sub.k0_auto = false;
        // Keep computed value as override
        sub.k0_override = sub.k0_computed || 1.0;
    }

    updateFoodK0Matrix();
}

// Update k0 override for substance
function updateK0Override(substanceId, value) {
    const sub = state.substances.find(s => (s.id || s.cid || s.name) === substanceId);
    if (!sub) return;

    if (value && parseFloat(value) > 0) {
        sub.k0_override = parseFloat(value);
    } else {
        delete sub.k0_override;
    }

    // Refresh the k0 matrix to update K_FP
    updateFoodK0Matrix();
}

// Called when simulant changes - recalculates k0 for all substances
function onSimulantChange() {
    updateFoodK0Matrix();
}

// Update CF0 (initial concentration in food) for substance
function updateCF0(substanceId, value) {
    const numValue = parseFloat(value) || 0;

    if (numValue > 0) {
        state.CF0[substanceId] = {
            value: numValue,
            unit: 'mg/kg',
        };
    } else {
        delete state.CF0[substanceId];
    }

    // Refresh the matrix to update styling
    updateFoodK0Matrix();
}

// Geometry shape definitions with icons
const GEOMETRY_SHAPES = [
    { code: 'cylinder', icon: '🥫', name: 'Can/Jar', desc: 'Cylindrical container' },
    { code: 'bottle', icon: '🍼', name: 'Bottle', desc: 'Bottle with neck' },
    { code: 'rectangle', icon: '📦', name: 'Box/Tray', desc: 'Rectangular container' },
    { code: 'pouch', icon: '🧴', name: 'Pouch', desc: 'Flexible pouch/sachet' },
    { code: 'sphere', icon: '⚽', name: 'Sphere', desc: 'Spherical container' },
];

function initGeometry() {
    const shapeSelect = document.getElementById('geometry-shape');
    const dimContainer = document.getElementById('geometry-dimensions');
    const iconContainer = document.getElementById('geometry-shape-icons');

    const shapeDimensions = {
        cylinder: [
            { name: 'radius', label: 'Radius (mm)', default: 50 },
            { name: 'height', label: 'Height (mm)', default: 100 },
        ],
        bottle: [
            { name: 'body_radius', label: 'Body Radius (mm)', default: 40 },
            { name: 'body_height', label: 'Body Height (mm)', default: 200 },
            { name: 'neck_radius', label: 'Neck Radius (mm)', default: 15 },
            { name: 'neck_height', label: 'Neck Height (mm)', default: 50 },
        ],
        rectangle: [
            { name: 'length', label: 'Length (mm)', default: 200 },
            { name: 'width', label: 'Width (mm)', default: 150 },
            { name: 'height', label: 'Height (mm)', default: 50 },
        ],
        pouch: [
            { name: 'width', label: 'Width (mm)', default: 150 },
            { name: 'height', label: 'Height (mm)', default: 200 },
            { name: 'thickness', label: 'Thickness (mm)', default: 10 },
        ],
        sphere: [
            { name: 'radius', label: 'Radius (mm)', default: 50 },
        ],
    };

    function updateDimensions() {
        const shape = shapeSelect.value;
        const dims = shapeDimensions[shape] || [];

        dimContainer.innerHTML = dims.map(d => `
            <div class="flex items-center space-x-2">
                <label class="w-32 text-sm text-gray-600 dark:text-gray-400">${d.label}</label>
                <input type="number" class="flex-1 border dark:border-gray-600 rounded px-2 py-1 geom-dim
                                           bg-white dark:bg-gray-700 text-gray-800 dark:text-white"
                       data-dim="${d.name}" value="${d.default}" min="0.1" step="1">
            </div>
        `).join('');

        // Update state
        state.geometry.shape = shape;
        state.geometry.dimensions = {};
        dims.forEach(d => {
            state.geometry.dimensions[d.name] = d.default;
        });

        // Add listeners
        dimContainer.querySelectorAll('.geom-dim').forEach(input => {
            input.addEventListener('input', () => {
                state.geometry.dimensions[input.dataset.dim] = parseFloat(input.value) || 0;
                calculateGeometry();
            });
        });

        calculateGeometry();
    }

    // Create geometry shape icons
    if (iconContainer) {
        iconContainer.innerHTML = GEOMETRY_SHAPES.map(gs => `
            <button class="geom-shape-btn flex flex-col items-center p-2 border-2 rounded-lg
                           bg-white dark:bg-gray-700 hover:border-sfppy-green transition-colors min-w-[65px]
                           ${state.geometry.shape === gs.code ? 'border-sfppy-green bg-green-50 dark:bg-green-900/30' : 'border-gray-300 dark:border-gray-600'}"
                    data-shape="${gs.code}" title="${gs.desc}">
                <span class="text-2xl">${gs.icon}</span>
                <span class="text-[10px] text-gray-600 dark:text-gray-400 text-center leading-tight mt-1">${gs.name}</span>
            </button>
        `).join('');

        iconContainer.querySelectorAll('.geom-shape-btn').forEach(btn => {
            btn.addEventListener('click', () => {
                const shape = btn.dataset.shape;
                shapeSelect.value = shape;

                // Update visual selection
                iconContainer.querySelectorAll('.geom-shape-btn').forEach(b => {
                    b.classList.remove('border-sfppy-green', 'bg-green-50', 'dark:bg-green-900/30');
                    b.classList.add('border-gray-300', 'dark:border-gray-600');
                });
                btn.classList.remove('border-gray-300', 'dark:border-gray-600');
                btn.classList.add('border-sfppy-green', 'bg-green-50', 'dark:bg-green-900/30');

                updateDimensions();
            });
        });
    }

    shapeSelect.addEventListener('change', updateDimensions);
    updateDimensions();
}

async function calculateGeometry() {
    try {
        // Validate geometry data before sending
        if (!state.geometry.shape) {
            state.geometry.shape = 'cylinder';
        }
        if (!state.geometry.dimensions || typeof state.geometry.dimensions !== 'object') {
            state.geometry.dimensions = { radius: 50, height: 100 };
        }

        // Clean dimensions - ensure numeric values
        const cleanDims = {};
        for (const [key, val] of Object.entries(state.geometry.dimensions)) {
            if (val !== undefined && val !== null) {
                cleanDims[key] = typeof val === 'object' ? (val.value || 0) : parseFloat(val) || 0;
            }
        }

        const payload = {
            shape: state.geometry.shape,
            dimensions: cleanDims,
        };

        const result = await apiPost('/food/geometry/calculate', payload);
        if (result && result.success !== false) {
            const vol = result.volume_cm3 ?? result.volume ?? 0;
            const surf = result.surface_cm2 ?? result.surface ?? 0;
            const ratio = result.vs_ratio_cm ?? (surf > 0 ? vol / surf : 0);

            // Store in state for report access
            state.geometry.volume_cm3 = vol;
            state.geometry.surface_cm2 = surf;
            state.geometry.vs_ratio_cm = ratio;

            document.getElementById('geom-volume').textContent = `${vol.toFixed(2)} cm³`;
            document.getElementById('geom-surface').textContent = `${surf.toFixed(2)} cm²`;
            document.getElementById('geom-vs-ratio').textContent = `${ratio.toFixed(3)} cm`;
        }
    } catch (e) {
        console.warn('Geometry calculation failed:', e);
        // Set defaults on error
        state.geometry.volume_cm3 = null;
        state.geometry.surface_cm2 = null;
        state.geometry.vs_ratio_cm = null;
        document.getElementById('geom-volume').textContent = '-- cm³';
        document.getElementById('geom-surface').textContent = '-- cm²';
        document.getElementById('geom-vs-ratio').textContent = '-- cm';
    }
}

function applyConditionPreset(code, presets) {
    const preset = presets.find(p => p.code === code);
    if (!preset) return;

    // Clear existing steps
    const container = document.getElementById('contact-steps');
    container.innerHTML = '';
    state.steps = [];

    // Add preset steps
    preset.steps.forEach((s, i) => {
        const step = {
            index: i + 1,
            temperature_C: s.temperature_C || 25,
            duration: s.duration,
            duration_unit: s.unit || 'days',
            with_food: s.with_food !== false,
            type: s.type || 'storage',
        };
        state.steps.push(step);
        addStepCard(step);
    });
}

function addStepCard(step) {
    const container = document.getElementById('contact-steps');
    // Ensure setoff_type has a default value
    if (!step.setoff_type) step.setoff_type = 'stacked';

    // Determine card border color based on mode
    const borderClass = step.with_food
        ? 'border-gray-200 dark:border-gray-600'
        : 'border-amber-300 dark:border-amber-600';

    const html = `
        <div class="border-2 rounded-lg p-4 mb-3 step-card ${borderClass}" data-step="${step.index}">
            <div class="flex items-center justify-between mb-3">
                <div class="flex items-center space-x-3">
                    <span class="font-medium">Step ${step.index}</span>
                    <span class="step-mode-badge px-2 py-0.5 text-xs rounded ${step.with_food ? 'bg-green-100 text-green-800' : 'bg-amber-100 text-amber-800'}">
                        ${step.with_food ? '🍽️ Food Contact' : '📦 Set-off'}
                    </span>
                </div>
                <div class="flex items-center space-x-4">
                    ${step.index > 1 ? `<button class="text-red-500 hover:text-red-700 remove-step" data-step="${step.index}">Remove</button>` : ''}
                </div>
            </div>

            <!-- Contact Mode Selector -->
            <div class="mb-3 p-2 bg-gray-50 dark:bg-gray-700 rounded flex items-center space-x-4">
                <span class="text-sm font-medium text-gray-700 dark:text-gray-300">Contact Mode:</span>
                <label class="flex items-center space-x-2 cursor-pointer">
                    <input type="radio" name="contact-mode-${step.index}" class="step-contact-mode" data-step="${step.index}" value="food" ${step.with_food ? 'checked' : ''}>
                    <span class="text-sm text-gray-600 dark:text-gray-400">🍽️ Food Contact</span>
                </label>
                <label class="flex items-center space-x-2 cursor-pointer">
                    <input type="radio" name="contact-mode-${step.index}" class="step-contact-mode" data-step="${step.index}" value="setoff" ${!step.with_food ? 'checked' : ''}>
                    <span class="text-sm text-gray-600 dark:text-gray-400">📦 Set-off (no food)</span>
                </label>
                <!-- Set-off type selector (shown only in set-off mode) -->
                <div class="setoff-type-inline" style="display: ${step.with_food ? 'none' : 'flex'}; align-items: center; margin-left: auto;">
                    <span class="text-sm text-gray-600 dark:text-gray-400 mr-2">Config:</span>
                    <select class="border rounded px-2 py-1 text-sm step-setoff-type bg-white dark:bg-gray-600" data-step="${step.index}">
                        <option value="stacked" ${step.setoff_type === 'stacked' ? 'selected' : ''}>Stacked</option>
                        <option value="rolled" ${step.setoff_type === 'rolled' ? 'selected' : ''}>Rolled</option>
                    </select>
                </div>
            </div>

            <div class="grid grid-cols-3 gap-4">
                <div>
                    <label class="text-sm text-gray-600 dark:text-gray-400">Type</label>
                    <select class="w-full border rounded px-3 py-2 step-type dark:bg-gray-700" data-step="${step.index}">
                        <option value="storage" ${step.type === 'storage' ? 'selected' : ''}>Storage</option>
                        <option value="hotfill" ${step.type === 'hotfill' ? 'selected' : ''}>Hot Fill</option>
                        <option value="transport" ${step.type === 'transport' ? 'selected' : ''}>Transport</option>
                    </select>
                </div>
                <div>
                    <label class="text-sm text-gray-600 dark:text-gray-400">Temperature</label>
                    <div class="flex">
                        <input type="number" class="flex-1 border rounded-l px-3 py-2 step-temp dark:bg-gray-700"
                               data-step="${step.index}" value="${step.temperature_C}">
                        <span class="border-t border-r border-b rounded-r px-3 py-2 bg-gray-50 dark:bg-gray-600">°C</span>
                    </div>
                </div>
                <div>
                    <label class="text-sm text-gray-600 dark:text-gray-400">Duration</label>
                    <div class="flex">
                        <input type="number" class="flex-1 border rounded-l px-3 py-2 step-duration dark:bg-gray-700"
                               data-step="${step.index}" value="${step.duration}" min="0.1">
                        <select class="border-t border-r border-b rounded-r px-2 py-2 bg-gray-50 dark:bg-gray-600 step-duration-unit" data-step="${step.index}">
                            <option value="minutes" ${step.duration_unit === 'minutes' ? 'selected' : ''}>min</option>
                            <option value="hours" ${step.duration_unit === 'hours' ? 'selected' : ''}>hours</option>
                            <option value="days" ${step.duration_unit === 'days' ? 'selected' : ''}>days</option>
                            <option value="months" ${step.duration_unit === 'months' ? 'selected' : ''}>months</option>
                        </select>
                    </div>
                </div>
            </div>
            <!-- Set-off info banner (shown when in set-off mode) -->
            <div class="setoff-info-banner mt-3 p-2 bg-amber-50 dark:bg-amber-900/30 border border-amber-200 dark:border-amber-700 rounded text-sm text-amber-800 dark:text-amber-300" style="display: ${step.with_food ? 'none' : 'block'}">
                <strong>📦 Set-off mode:</strong> Simulates migration between stacked/rolled packaging layers <em>before</em> food contact (periodic boundary conditions). No food/simulant involved in this step.
            </div>
        </div>
    `;
    container.insertAdjacentHTML('beforeend', html);

    // Add listeners
    initStepListeners(step.index);
}

function initStepListeners(index) {
    const card = document.querySelector(`.step-card[data-step="${index}"]`);
    if (!card) return;

    card.querySelector('.step-type')?.addEventListener('change', (e) => {
        state.steps[index - 1].type = e.target.value;
    });

    card.querySelector('.step-temp')?.addEventListener('input', (e) => {
        state.steps[index - 1].temperature_C = parseFloat(e.target.value) || 25;
    });

    card.querySelector('.step-duration')?.addEventListener('input', (e) => {
        state.steps[index - 1].duration = parseFloat(e.target.value) || 1;
    });

    card.querySelector('.step-duration-unit')?.addEventListener('change', (e) => {
        state.steps[index - 1].duration_unit = e.target.value;
    });

    // Contact mode radio buttons (Food Contact vs Set-off)
    card.querySelectorAll('.step-contact-mode').forEach(radio => {
        radio.addEventListener('change', (e) => {
            const withFood = e.target.value === 'food';
            state.steps[index - 1].with_food = withFood;
            toggleSetoffUI(card, withFood);
        });
    });

    // Legacy checkbox support (if still present in HTML template)
    card.querySelector('.step-with-food')?.addEventListener('change', (e) => {
        state.steps[index - 1].with_food = e.target.checked;
        toggleSetoffUI(card, e.target.checked);
    });

    card.querySelector('.step-setoff-type')?.addEventListener('change', (e) => {
        state.steps[index - 1].setoff_type = e.target.value;
    });

    card.querySelector('.remove-step')?.addEventListener('click', () => {
        removeStep(index);
    });
}

function toggleSetoffUI(card, withFood) {
    // Show/hide set-off configuration based on contact mode

    // Update card border color
    if (withFood) {
        card.classList.remove('border-amber-300', 'dark:border-amber-600');
        card.classList.add('border-gray-200', 'dark:border-gray-600');
    } else {
        card.classList.remove('border-gray-200', 'dark:border-gray-600');
        card.classList.add('border-amber-300', 'dark:border-amber-600');
    }

    // Update mode badge
    const badge = card.querySelector('.step-mode-badge');
    if (badge) {
        badge.className = `step-mode-badge px-2 py-0.5 text-xs rounded ${withFood ? 'bg-green-100 text-green-800' : 'bg-amber-100 text-amber-800'}`;
        badge.textContent = withFood ? '🍽️ Food Contact' : '📦 Set-off';
    }

    // Show/hide set-off type selector (inline version)
    const setoffInline = card.querySelector('.setoff-type-inline');
    if (setoffInline) {
        setoffInline.style.display = withFood ? 'none' : 'flex';
    }

    // Legacy: show/hide block-level set-off container
    const setoffContainer = card.querySelector('.setoff-type-container');
    if (setoffContainer) {
        setoffContainer.style.display = withFood ? 'none' : 'block';
    }

    // Show/hide set-off info banner
    const setoffBanner = card.querySelector('.setoff-info-banner');
    if (setoffBanner) {
        setoffBanner.style.display = withFood ? 'none' : 'block';
    }
}

function removeStep(index) {
    state.steps = state.steps.filter(s => s.index !== index);
    // Re-index
    state.steps.forEach((s, i) => s.index = i + 1);

    // Rebuild UI
    const container = document.getElementById('contact-steps');
    container.innerHTML = '';
    state.steps.forEach(s => addStepCard(s));
}

function initFoodTab() {
    document.getElementById('add-step').addEventListener('click', () => {
        const newStep = {
            index: state.steps.length + 1,
            temperature_C: 25,
            duration: 10,
            duration_unit: 'days',
            with_food: true,
            type: 'storage',
            setoff_type: 'stacked',
        };
        state.steps.push(newStep);
        addStepCard(newStep);
    });

    // Rebuild all step cards from state to use the new UI structure
    // This replaces the static HTML template with the dynamic JavaScript version
    const container = document.getElementById('contact-steps');
    container.innerHTML = '';
    state.steps.forEach(s => addStepCard(s));

    // Food type and texture icons are initialized in loadFoodData()
}

// ========== SUBSTANCES TAB ==========
async function loadSubstancesData() {
    // Update C0 matrix to reflect latest layer polymers from Assembly tab
    updateC0Matrix();

    // Load common substances
    try {
        const data = await apiGet('/substances/common');
        const container = document.getElementById('common-substances');
        container.innerHTML = data.substances.map(s => `
            <button class="common-sub-btn p-2 border rounded text-left text-sm hover:border-blue-500"
                    data-substance="${s.id}">
                <div class="font-medium">${s.name}</div>
                <div class="text-xs text-gray-500">${s.category}</div>
            </button>
        `).join('');

        container.querySelectorAll('.common-sub-btn').forEach(btn => {
            btn.addEventListener('click', () => addSubstance(btn.dataset.substance, data.substances));
        });
    } catch (e) {
        console.error('Failed to load substances:', e);
    }
}

async function addSubstance(id, allSubstances) {
    const sub = allSubstances.find(s => s.id === id);
    if (!sub || state.substances.find(s => s.id === id)) return;

    const newSubstance = {
        ...sub,
        layer_index: 1,
        C0: sub.typical_C0?.LDPE || 1000,
        SML: sub.SML || null,
    };

    state.substances.push(newSubstance);

    // Initialize C0 matrix for this substance
    const subId = newSubstance.id || newSubstance.cid || newSubstance.name;
    const firstLayerPolymer = state.layers[0]?.polymer || 'LDPE';
    state.c0Matrix[subId] = {};
    state.c0Matrix[subId][1] = {
        value: sub.typical_C0?.[firstLayerPolymer] || 1000,
        unit: state.defaultC0Unit,
    };
    state.c0Units[subId] = state.defaultC0Unit;

    // Update all UI components
    updateSubstancesList();
    updateC0Matrix();
    updateLayerVisualization();
    await updateAssemblySubstances();
    updateFoodK0Matrix();
}

function updateSubstancesList() {
    const container = document.getElementById('selected-substances');

    if (state.substances.length === 0) {
        container.innerHTML = `
            <div class="text-gray-500 dark:text-gray-400 text-sm p-3 border border-dashed dark:border-gray-600 rounded text-center">
                No substances added. Search above or select from common substances.
            </div>
        `;
        return;
    }

    // Simplified view: just substance selection with basic info
    // D/k assignments are done in Assembly tab, k0 in Food tab
    container.innerHTML = state.substances.map((s, i) => {
        const subId = s.id || s.cid || s.name;

        // Build regulatory flags - check both 'authorized' flag and presence of SML/FCM
        let regFlags = '';
        const reg = s.regulatory || {};
        // EU: authorized if explicitly set OR has SML/FCM number
        if (reg.EU?.authorized || reg.EU?.SML || reg.EU?.FCM || s.SML) {
            regFlags += `<span class="reg-flag reg-eu" title="EU 10/2011 Authorized">🇪🇺</span>`;
        }
        // US: authorized if explicitly set OR has FCM_No/CFR reference
        if (reg.US?.authorized || reg.US?.FCM_No || reg.US?.CFR) {
            regFlags += `<span class="reg-flag reg-us" title="FDA Listed">🇺🇸</span>`;
        }
        // CN: authorized if explicitly set OR has SML/FCA_No
        if (reg.CN?.authorized || reg.CN?.SML || reg.CN?.FCA_No) {
            regFlags += `<span class="reg-flag reg-cn" title="CN GB 9685 Listed">🇨🇳</span>`;
        }
        if (!regFlags) {
            regFlags = `<span class="text-gray-400 text-xs">Not listed</span>`;
        }

        // Cramer class and ToxTree display - Enhanced for NIAS
        let cramerDisplay = '';
        let alertsDisplay = '';
        let ttcDisplay = '';

        const toxtree = s.toxtree || {};
        const cramerClass = s.cramer_class || toxtree.cramer_class;
        const cramerValue = s.cramer_value || toxtree.cramer_value;
        const ttc = s.ttc || toxtree.ttc;
        const cfTtc = s.cf_ttc || toxtree.cf_ttc;
        const hasAlerts = toxtree.has_alerts;
        const nAlerts = toxtree.nalerts || 0;

        if (cramerClass) {
            let cramerColor = 'text-gray-500';
            let cramerBg = 'bg-gray-100 dark:bg-gray-700';
            if (cramerValue === 1) { cramerColor = 'text-green-700'; cramerBg = 'bg-green-100 dark:bg-green-900'; }
            else if (cramerValue === 2) { cramerColor = 'text-yellow-700'; cramerBg = 'bg-yellow-100 dark:bg-yellow-900'; }
            else if (cramerValue === 3) { cramerColor = 'text-red-700'; cramerBg = 'bg-red-100 dark:bg-red-900'; }

            // Build tooltip with TTC info
            let tooltip = `Cramer ${cramerClass}`;
            if (ttc) tooltip += ` | TTC: ${ttc} µg/kg bw/day`;
            if (cfTtc) tooltip += ` | Max: ${cfTtc.toFixed(3)} mg/kg food`;
            else if (ttc) tooltip += ` | Max: ${(ttc * 60 / 1000).toFixed(3)} mg/kg food`;

            cramerDisplay = `<span class="px-1.5 py-0.5 rounded text-xs font-medium ${cramerBg} ${cramerColor}" title="${tooltip}">${cramerClass}</span>`;

            // Show TTC as separate badge if no SML (NIAS case)
            if (!s.SML && ttc) {
                const maxExp = cfTtc || (ttc * 60 / 1000);
                ttcDisplay = `<span class="px-1.5 py-0.5 rounded text-xs bg-blue-100 dark:bg-blue-900 text-blue-700 dark:text-blue-300" title="TTC-based limit for 60 kg adult">TTC: ${maxExp.toFixed(3)} mg/kg</span>`;
            }
        }

        // Structural alerts badge - with genotoxic detection
        const hasGenotoxic = toxtree.has_genotoxic;
        const genotoxCfTtc = toxtree.genotox_cf_ttc;

        if (hasGenotoxic) {
            // GENOTOXIC - Critical warning with stricter threshold
            alertsDisplay = `<span class="px-1.5 py-0.5 rounded text-xs font-bold bg-purple-600 text-white animate-pulse" title="GENOTOXIC ALERT: TTC = 0.0025 µg/kg bw/day → ${genotoxCfTtc?.toFixed(6) || '0.00015'} mg/kg food">☢️ GENOTOXIC</span>`;
            // Override TTC display for genotoxic
            if (!s.SML) {
                ttcDisplay = `<span class="px-1.5 py-0.5 rounded text-xs bg-purple-100 dark:bg-purple-900 text-purple-700 dark:text-purple-300 font-medium" title="Genotoxic TTC threshold for 60 kg adult">TTC: ${genotoxCfTtc?.toFixed(6) || '0.00015'} mg/kg</span>`;
            }
        } else if (hasAlerts || nAlerts > 0) {
            // Non-genotoxic alerts - informational (cosmetics/homecare)
            alertsDisplay = `<span class="px-1.5 py-0.5 rounded text-xs font-medium bg-orange-100 dark:bg-orange-900 text-orange-700 dark:text-orange-300" title="${nAlerts} informational alert(s) - cosmetics/homecare relevance">⚠️ ${nAlerts}</span>`;
        } else if (cramerClass) {
            // Only show "No alerts" if we have ToxTree data
            alertsDisplay = `<span class="px-1.5 py-0.5 rounded text-xs bg-green-50 dark:bg-green-900/30 text-green-600 dark:text-green-400" title="No structural alerts">✓</span>`;
        }

        // CAS display - handle array
        let casDisplay = '-';
        if (Array.isArray(s.cas)) {
            casDisplay = s.cas.slice(0, 2).join(', ');
            if (s.cas.length > 2) casDisplay += '...';
        } else if (s.cas) {
            casDisplay = s.cas;
        }

        return `
            <div class="substance-detail-card sub-card" data-index="${i}" data-id="${subId}">
                <div class="flex items-start gap-3">
                    ${s.cid ? `<img src="${getStructureImageUrl(s.cid)}"
                                   class="w-14 h-14 molecule-img flex-shrink-0 cursor-pointer"
                                   title="Click to enlarge"
                                   data-cid="${s.cid}"
                                   onclick="showStructureModal(${s.cid})"
                                   onerror="this.onerror=null; this.src='/api/substances/structure/${s.cid}.png';">` : ''}
                    <div class="flex-1 min-w-0">
                        <div class="flex items-start justify-between">
                            <div>
                                <h4 class="font-medium text-gray-800 dark:text-white truncate">${s.name}</h4>
                                <p class="text-xs text-gray-500 dark:text-gray-400">
                                    CAS: ${casDisplay} ${s.mw ? `| MW: ${parseFloat(s.mw).toFixed(2)}` : ''}
                                    ${s.logP ? `| logP: ${parseFloat(s.logP).toFixed(1)}` : ''}
                                </p>
                            </div>
                            <button class="text-red-500 hover:text-red-700 text-sm remove-substance" data-index="${i}">✕</button>
                        </div>
                        <div class="flex items-center gap-2 mt-2 flex-wrap">
                            <div class="flex items-center gap-1">${regFlags}</div>
                            ${s.SML ? `<span class="regulatory-badge authorized">SML: ${s.SML} mg/kg</span>` : ''}
                            ${cramerDisplay}
                            ${alertsDisplay}
                            ${ttcDisplay}
                        </div>
                    </div>
                </div>
                <p class="text-xs text-gray-400 dark:text-gray-500 mt-2 italic">
                    Set concentrations in Assembly tab, k0 in Food tab.
                </p>
            </div>
        `;
    }).join('');

    // Add listeners
    container.querySelectorAll('.remove-substance').forEach(btn => {
        btn.addEventListener('click', async () => {
            const idx = parseInt(btn.dataset.index);
            const subId = state.substances[idx].id || state.substances[idx].cid || state.substances[idx].name;

            // Remove from C0 matrix
            delete state.c0Matrix[subId];
            delete state.c0Units[subId];

            // Remove from substances array
            state.substances.splice(idx, 1);

            updateSubstancesList();
            updateC0Matrix();
            updateLayerVisualization();
            await updateAssemblySubstances();
            updateFoodK0Matrix();
        });
    });
}

function initSubstancesTab() {
    const searchInput = document.getElementById('substance-search');
    const searchBtn = document.getElementById('search-btn');
    const resultsContainer = document.getElementById('search-results');

    searchBtn.addEventListener('click', async () => {
        const query = searchInput.value.trim();
        if (!query) return;

        searchBtn.textContent = '...';
        searchBtn.disabled = true;

        try {
            const data = await apiGet(`/substances/search?q=${encodeURIComponent(query)}`);
            if (data.results.length > 0) {
                resultsContainer.classList.remove('hidden');
                resultsContainer.innerHTML = data.results.map(r => `
                    <div class="p-3 hover:bg-gray-50 dark:hover:bg-gray-700 cursor-pointer border-b dark:border-gray-700 search-result flex items-center gap-3"
                         data-id="${r.id}" data-source="${r.source}" data-cid="${r.cid || ''}">
                        ${r.cid ? `<img src="${getStructureImageUrl(r.cid)}"
                                       class="w-12 h-12 molecule-img cursor-pointer"
                                       title="Click to enlarge"
                                       data-cid="${r.cid}"
                                       onclick="event.stopPropagation(); showStructureModal(${r.cid})"
                                       onerror="this.onerror=null; this.src='/api/substances/structure/${r.cid}.png';">` : ''}
                        <div class="flex-1 min-w-0">
                            <div class="font-medium text-gray-800 dark:text-gray-200 truncate">${r.name}</div>
                            <div class="text-xs text-gray-500 dark:text-gray-400">
                                ${r.cas ? `CAS: ${r.cas}` : ''} ${r.mw ? `| MW: ${parseFloat(r.mw).toFixed(2)}` : ''} | ${r.source}
                            </div>
                        </div>
                        <button class="btn-secondary text-xs preview-substance-btn" data-id="${r.id}">Details</button>
                    </div>
                `).join('');

                // Click on row to add quickly
                resultsContainer.querySelectorAll('.search-result').forEach(div => {
                    div.addEventListener('click', async (e) => {
                        // Ignore if clicked on Details button
                        if (e.target.classList.contains('preview-substance-btn')) return;

                        const result = data.results.find(r => r.id === div.dataset.id);
                        if (result) {
                            await addSubstanceToState(result);
                            resultsContainer.classList.add('hidden');
                            searchInput.value = '';
                        }
                    });
                });

                // Details button to show preview
                resultsContainer.querySelectorAll('.preview-substance-btn').forEach(btn => {
                    btn.addEventListener('click', async (e) => {
                        e.stopPropagation();
                        const resultDiv = btn.closest('.search-result');
                        const result = data.results.find(r => r.id === resultDiv.dataset.id);
                        if (result) {
                            await showSubstancePreview(result);
                        }
                    });
                });
            } else {
                resultsContainer.classList.remove('hidden');
                resultsContainer.innerHTML = '<div class="p-3 text-gray-500 dark:text-gray-400">No results found</div>';
            }
        } catch (e) {
            console.error('Search failed:', e);
            resultsContainer.classList.remove('hidden');
            resultsContainer.innerHTML = `<div class="p-3 text-red-500">Search failed: ${e.message}</div>`;
        } finally {
            searchBtn.textContent = '🔍 Search';
            searchBtn.disabled = false;
        }
    });

    searchInput.addEventListener('keypress', (e) => {
        if (e.key === 'Enter') searchBtn.click();
    });

    // Add previewed substance button
    document.getElementById('add-previewed-substance')?.addEventListener('click', async () => {
        if (state.previewedSubstance) {
            await addSubstanceToState(state.previewedSubstance);
            document.getElementById('substance-preview').classList.add('hidden');
            document.getElementById('search-results').classList.add('hidden');
            document.getElementById('substance-search').value = '';
        }
    });
}

async function showSubstancePreview(basicResult) {
    const previewDiv = document.getElementById('substance-preview');
    if (!previewDiv) return;

    // Try to get detailed info from /substances/detail endpoint
    // Prefer CAS for lookup (most reliable), then name, then try CID
    // Note: migrantToxtree doesn't support direct CID lookups - needs name or CAS
    let detailData = null;
    const queryParam = basicResult.cas || basicResult.name || (basicResult.cid ? String(basicResult.cid) : null);
    console.log('[ToxTree Debug] showSubstancePreview called with:', basicResult);
    console.log('[ToxTree Debug] Using query param:', queryParam);

    try {
        const response = await apiGet(`/substances/detail/${encodeURIComponent(queryParam)}`);
        console.log('[ToxTree Debug] Detail response:', response);
        if (response.success) {
            detailData = response.substance;
            console.log('[ToxTree Debug] Detail data toxtree:', detailData?.toxtree);
        }
    } catch (e) {
        console.warn('Could not fetch detailed data:', e);
    }

    // Merge basic result with detailed data
    const substance = { ...basicResult, ...detailData };
    state.previewedSubstance = substance;

    // Update preview UI
    const cid = substance.cid;
    const structImg = document.getElementById('preview-structure');
    if (cid) {
        structImg.style.display = '';  // Show image
        setStructureImage(structImg, cid, 'l');  // Use large size for preview
    } else {
        structImg.src = '';
        structImg.style.display = 'none';
    }
    document.getElementById('preview-name').textContent = substance.name || 'Unknown';

    // Alternative names / synonyms with ellipsis
    let altNames = '';
    const synonyms = substance.synonyms || [];
    if (synonyms.length > 0) {
        // Show synonyms with ellipsis to indicate more exist
        altNames = synonyms.join(', ') + '...';
    } else if (substance.iupac_name && substance.iupac_name !== substance.name) {
        altNames = substance.iupac_name;
    }
    document.getElementById('preview-alt-names').textContent = altNames;

    // CAS - handle array of CAS numbers
    const casValue = substance.cas;
    let casDisplay = '-';
    if (Array.isArray(casValue)) {
        casDisplay = casValue.join(', ');
    } else if (casValue) {
        casDisplay = casValue;
    }
    document.getElementById('preview-cas').textContent = casDisplay;

    // Properties
    document.getElementById('preview-mw').textContent = substance.mw
        ? `${parseFloat(substance.mw).toFixed(2)} g/mol`
        : '-';
    document.getElementById('preview-logp').textContent = substance.logP
        ? parseFloat(substance.logP).toFixed(2)
        : '-';
    document.getElementById('preview-formula').textContent = substance.formula || '-';

    // Regulatory badges
    const badgesDiv = document.getElementById('preview-badges');
    let badgesHtml = '';

    const reg = substance.regulatory || {};
    const hasSML = substance.SML || reg.EU?.SML;

    // EU Regulatory - check authorized flag, SML, or FCM
    if (reg.EU?.authorized || reg.EU?.SML || reg.EU?.FCM || hasSML) {
        badgesHtml += `<span class="regulatory-badge authorized">🇪🇺 EU Authorized</span>`;
        const smlVal = reg.EU?.SML || substance.SML;
        if (smlVal) {
            badgesHtml += `<span class="regulatory-badge authorized">SML: ${smlVal} mg/kg</span>`;
        }
        if (reg.EU?.FCM) {
            badgesHtml += `<span class="regulatory-badge" title="FCM No: ${reg.EU.FCM}">FCM ${reg.EU.FCM}</span>`;
        }
    } else if (substance.regulatory) {
        badgesHtml += `<span class="regulatory-badge not-listed">🇪🇺 Not Listed</span>`;
    }

    // US Regulatory
    if (reg.US?.authorized || reg.US?.FCM_No || reg.US?.CFR) {
        badgesHtml += `<span class="regulatory-badge authorized">🇺🇸 FDA Listed</span>`;
    }

    // CN Regulatory
    if (reg.CN?.authorized || reg.CN?.SML || reg.CN?.FCA_No) {
        badgesHtml += `<span class="regulatory-badge authorized">🇨🇳 CN Listed</span>`;
        if (reg.CN?.SML) {
            badgesHtml += `<span class="regulatory-badge" title="China SML">🇨🇳SML: ${reg.CN.SML}</span>`;
        }
    }

    if (!badgesHtml) {
        badgesHtml = '<span class="regulatory-badge not-listed">Status unknown</span>';
    }
    badgesDiv.innerHTML = badgesHtml;

    // Cramer class / ToxTree - Enhanced for NIAS assessment
    const cramerSpan = document.getElementById('preview-cramer');
    const ttcSpan = document.getElementById('preview-ttc');
    const cfTtcSpan = document.getElementById('preview-cfttc');
    const alertsSection = document.getElementById('preview-alerts-section');
    const alertsBadge = document.getElementById('preview-alerts-badge');
    const alertsText = document.getElementById('preview-alerts-text');

    console.log('[ToxTree Debug] substance.toxtree:', substance.toxtree);

    if (substance.toxtree?.cramer_class) {
        console.log('[ToxTree Debug] Displaying ToxTree data');
        const cramerClass = substance.toxtree.cramer_class;
        const cramerValue = substance.toxtree.cramer_value;
        const ttc = substance.toxtree.ttc;
        const cfTtc = substance.toxtree.cf_ttc;
        const hasAlerts = substance.toxtree.has_alerts;
        const nAlerts = substance.toxtree.nalerts || 0;

        // Show Cramer class
        cramerSpan.textContent = cramerClass;

        // Color by class
        if (cramerValue === 1) {
            cramerSpan.className = 'tox-alert tox-class-I ml-1 font-medium text-green-600';
        } else if (cramerValue === 2) {
            cramerSpan.className = 'tox-alert tox-class-II ml-1 font-medium text-yellow-600';
        } else if (cramerValue === 3) {
            cramerSpan.className = 'tox-alert tox-class-III ml-1 font-medium text-red-600';
        } else {
            cramerSpan.className = 'tox-alert ml-1 font-medium';
        }

        // TTC value
        if (ttc !== null && ttc !== undefined) {
            ttcSpan.textContent = `${ttc} µg/kg bw/day`;
        } else {
            ttcSpan.textContent = '-';
        }

        // CF_TTC - Max exposure for 60 kg adult (mg/kg food)
        if (cfTtc !== null && cfTtc !== undefined) {
            cfTtcSpan.textContent = `${cfTtc.toFixed(4)} mg/kg food`;
            cfTtcSpan.title = 'Maximum concentration in food assuming 60 kg adult, 1 kg food/day';
        } else if (ttc !== null && ttc !== undefined) {
            // Calculate: TTC (µg/kg bw/day) × 60 kg / 1000 = mg/kg food
            const computed = ttc * 60 / 1000;
            cfTtcSpan.textContent = `${computed.toFixed(4)} mg/kg food`;
            cfTtcSpan.title = 'Computed: TTC × 60 kg ÷ 1000';
        } else {
            cfTtcSpan.textContent = '-';
        }

        // Structural alerts display - with genotoxic detection
        const hasGenotoxic = substance.toxtree.has_genotoxic;
        const genotoxTtc = substance.toxtree.genotox_ttc;
        const genotoxCfTtc = substance.toxtree.genotox_cf_ttc;
        const alertsList = substance.toxtree.alerts || [];

        if (hasGenotoxic) {
            // GENOTOXIC ALERT - Critical warning with stricter threshold
            alertsSection.classList.remove('hidden');
            alertsBadge.textContent = `☢️ GENOTOXIC`;
            alertsBadge.className = 'px-2 py-0.5 rounded text-xs font-bold bg-purple-600 text-white animate-pulse';
            alertsText.innerHTML = `<span class="text-purple-700 dark:text-purple-300 font-medium">TTC: ${genotoxTtc} µg/kg bw/day → ${genotoxCfTtc.toFixed(6)} mg/kg food</span>`;

            // Update the CF_TTC display to show genotoxic threshold
            if (cfTtcSpan) {
                cfTtcSpan.innerHTML = `<span class="line-through text-gray-400">${(cfTtc || ttc * 60 / 1000).toFixed(4)}</span> → <span class="text-purple-600 dark:text-purple-400 font-bold">${genotoxCfTtc.toFixed(6)} mg/kg</span>`;
                cfTtcSpan.title = 'Genotoxic threshold applies (0.0025 µg/kg bw/day)';
            }
        } else if (hasAlerts || nAlerts > 0) {
            // Non-genotoxic alerts - informational (cosmetics/homecare relevance)
            alertsSection.classList.remove('hidden');
            const genotoxCount = alertsList.filter(a => a.is_genotoxic).length;
            const infoCount = nAlerts - genotoxCount;
            alertsBadge.textContent = `⚠️ ${nAlerts} Alert${nAlerts !== 1 ? 's' : ''}`;
            alertsBadge.className = 'px-2 py-0.5 rounded text-xs font-medium bg-orange-100 text-orange-700 dark:bg-orange-900 dark:text-orange-300';
            alertsText.textContent = 'Informational alerts (cosmetics/homecare relevance)';
        } else {
            alertsSection.classList.remove('hidden');
            alertsBadge.textContent = '✓ No Alerts';
            alertsBadge.className = 'px-2 py-0.5 rounded text-xs font-medium bg-green-100 text-green-700 dark:bg-green-900 dark:text-green-300';
            alertsText.textContent = 'No structural alerts identified';
        }
    } else {
        cramerSpan.textContent = '-';
        cramerSpan.className = 'tox-alert ml-1 font-medium';
        if (ttcSpan) ttcSpan.textContent = '-';
        if (cfTtcSpan) cfTtcSpan.textContent = '-';
        if (alertsSection) alertsSection.classList.add('hidden');
    }

    // Show preview panel
    previewDiv.classList.remove('hidden');
}

async function addSubstanceToState(result) {
    const existingId = result.id || result.cid || result.name;
    if (state.substances.find(s => (s.id || s.cid || s.name) === existingId)) {
        showToast('Substance already added', 'warning');
        return;
    }

    // Always fetch full details including regulatory data from migrantToxtree
    const query = result.cid || result.cas || result.name;
    try {
        const detailResult = await apiGet(`/substances/detail/${encodeURIComponent(query)}`);
        if (detailResult.success && detailResult.substance) {
            const sub = detailResult.substance;
            result = {
                ...result,
                id: sub.cid ? `pubchem_${sub.cid}` : (result.id || sub.name),
                name: sub.name || result.name,
                cid: sub.cid || result.cid,
                cas: sub.cas || result.cas,
                mw: sub.mw || result.mw,
                logP: sub.logP || result.logP,
                formula: sub.formula,
                smiles: sub.smiles,
                cramer_class: sub.cramer_class,
                cramer_value: sub.cramer_value,
                ttc: sub.ttc,
                cf_ttc: sub.cf_ttc,
                regulatory: sub.regulatory,
                toxtree: sub.toxtree,
                SML: sub.SML || sub.regulatory?.EU?.SML,
            };
        }
    } catch (e) {
        console.warn('Failed to fetch full substance details:', e);
    }

    // Add to backend session if available
    if (state.sessionId) {
        try {
            const backendResult = await apiPost(`/sessions/${state.sessionId}/substances`, { query });
            if (backendResult.success) {
                console.log('Substance added to backend session');
            }
        } catch (e) {
            console.warn('Failed to add substance to backend:', e);
        }
    }

    const newSubstance = {
        ...result,
        id: existingId,
        SML: result.SML || result.regulatory?.EU?.SML || null,
    };

    // Compute D/k for the first layer's polymer
    const firstLayerPolymer = state.layers[0]?.polymer || 'LDPE';
    const temperature = state.steps[0]?.temperature_C || 40;

    try {
        const props = await computeDiffusivity(newSubstance.name, firstLayerPolymer, temperature);
        if (props) {
            newSubstance.D_computed = props.D;
            newSubstance.D0 = props.D0;
            newSubstance.k_computed = props.k;
            newSubstance.k0_computed = props.k0;
        }
    } catch (e) {
        console.warn('Failed to compute D/k for substance:', e);
    }

    state.substances.push(newSubstance);

    // Initialize C0 matrix for this substance (default: first layer gets typical C0)
    const subId = newSubstance.id;
    state.c0Matrix[subId] = {};
    state.c0Matrix[subId][1] = {
        value: result.typical_C0?.[firstLayerPolymer] || 1000,
        unit: state.defaultC0Unit,
    };
    state.c0Units[subId] = state.defaultC0Unit;

    // Sync assignment to backend (silently ignore - local state is authoritative)
    if (state.sessionId) {
        const c0Data = state.c0Matrix[subId][1];
        apiPost(`/sessions/${state.sessionId}/assignments`, {
            substance_id: subId,
            layer_index: 1,
            C0: c0Data.value,
            C0_unit: c0Data.unit,
        }).catch(() => {});
    }

    updateSubstancesList();
    updateC0Matrix();
    updateLayerVisualization();
    await updateAssemblySubstances();
    updateFoodK0Matrix();

    // Force DOM update for matrix visibility
    const matrixDiv = document.getElementById('assembly-substances-matrix');
    if (matrixDiv) {
        matrixDiv.style.display = 'none';
        matrixDiv.offsetHeight; // Force reflow
        matrixDiv.style.display = '';
    }
}

// ========== SIMULATE TAB ==========
function updateSimulateSummary() {
    document.getElementById('sim-layer-count').textContent = state.layers.length;
    document.getElementById('sim-step-count').textContent = state.steps.length;
    document.getElementById('sim-substance-count').textContent = state.substances.length;

    // Update restart step options based on current steps
    updateRestartStepOptions();
}

/**
 * Update the restart step dropdown based on current simulation steps.
 * The dropdown shows which step to restart FROM (using loaded session's final state).
 * Options are: "Last step" or specific step numbers from 1 to state.steps.length.
 */
function updateRestartStepOptions() {
    const stepSelect = document.getElementById('restart-from-step');
    if (!stepSelect) return;

    // Build options: show actual steps from current simulation
    const options = state.steps.map(s => {
        const stepLabel = s.with_food ? `${s.temperature_C}°C, ${s.duration} ${s.duration_unit}`
                                       : `${s.temperature_C}°C, ${s.duration} ${s.duration_unit} (set-off)`;
        return `<option value="${s.index}">Step ${s.index}: ${stepLabel}</option>`;
    });

    // Add "Last step" option at the end (most common use case)
    if (state.steps.length > 0) {
        const lastStep = state.steps[state.steps.length - 1];
        options.push(`<option value="last" selected>Last step (Step ${lastStep.index})</option>`);
    }

    stepSelect.innerHTML = options.join('');
}

async function initSimulateTab() {
    // Initialize restart step options based on current steps
    updateRestartStepOptions();

    // Load example sessions
    await loadExampleSessions();

    // Initialize session file handlers
    initSessionFileHandlers();

    // Load presets (only if container exists - may have been removed)
    try {
        const container = document.getElementById('simulation-presets');
        if (container) {
            const data = await apiGet('/simulation/presets');
            container.innerHTML = data.presets.map(p => `
                <button class="px-3 py-2 border border-gray-300 dark:border-gray-600 rounded text-sm
                               text-gray-700 dark:text-gray-200 bg-white dark:bg-gray-700
                               hover:bg-gray-50 dark:hover:bg-gray-600 hover:border-sfppy-green
                               transition-colors" data-preset="${p.code}">
                    ${p.name}
                </button>
            `).join('');

            container.querySelectorAll('button').forEach(btn => {
                btn.addEventListener('click', () => {
                    loadSimulationPreset(btn.dataset.preset, data.presets);
                });
            });
        }
    } catch (e) {
        console.error('Failed to load simulation presets:', e);
    }

    // Validate button
    document.getElementById('validate-all-btn').addEventListener('click', validateAllConfig);

    // Run button
    document.getElementById('run-simulation-btn').addEventListener('click', runSimulation);
}

function loadSimulationPreset(code, presets) {
    const preset = presets.find(p => p.code === code);
    if (!preset || !preset.config) return;

    // Apply layers
    if (preset.config.layers) {
        state.layers = preset.config.layers.map((l, i) => ({
            index: i + 1,
            polymer: l.polymer || 'LDPE',
            thickness: l.thickness || 100,
            thickness_unit: l.thickness_unit || 'um',
            temperature_C: 25,
            C0: l.C0 || 0,
        }));
        setLayerCount(state.layers.length);
    }

    // Apply steps
    if (preset.config.steps) {
        state.steps = preset.config.steps.map((s, i) => ({
            index: i + 1,
            temperature_C: s.temperature_C || 25,
            duration: s.duration || 10,
            duration_unit: s.duration_unit || 'days',
            with_food: s.with_food !== false,
            type: s.type || 'storage',
        }));
    }

    // Apply food
    if (preset.config.food) {
        state.food.category = preset.config.food.category || 'fatty';
    }

    document.getElementById('sim-name').value = preset.name;
    updateSimulateSummary();
    alert(`Loaded preset: ${preset.name}`);
}

async function validateAllConfig() {
    const statusDiv = document.getElementById('validation-status');

    // First try session validation if session exists
    if (state.sessionId) {
        try {
            const result = await apiGet(`/sessions/${state.sessionId}/validate`);
            statusDiv.classList.remove('hidden');

            if (result.valid) {
                statusDiv.className = 'mb-6 p-4 rounded border bg-green-50 border-green-200 dark:bg-green-900/30 dark:border-green-700';
                statusDiv.innerHTML = `
                    <div class="text-green-800 dark:text-green-400 font-medium">✓ Configuration Valid</div>
                    <div class="text-sm text-green-600 dark:text-green-300">
                        ${result.summary.layers_count} layers, ${result.summary.steps_count} steps,
                        ${result.summary.substances_count} substances, ${result.summary.assignments_count} assignments
                    </div>
                    ${result.warnings.length > 0 ? `<div class="mt-2 text-yellow-700 dark:text-yellow-400 text-sm">⚠️ ${result.warnings.join(', ')}</div>` : ''}
                `;
            } else {
                statusDiv.className = 'mb-6 p-4 rounded border bg-red-50 border-red-200 dark:bg-red-900/30 dark:border-red-700';
                statusDiv.innerHTML = `
                    <div class="text-red-800 dark:text-red-400 font-medium">✗ Validation Failed</div>
                    <ul class="text-sm text-red-600 dark:text-red-300 list-disc list-inside">
                        ${result.errors.map(e => `<li>${e}</li>`).join('')}
                    </ul>
                `;
            }
            return;
        } catch (e) {
            console.warn('Session validation failed, falling back to config validation:', e);
        }
    }

    // Fallback to config-based validation
    const config = buildSimulationConfig();
    try {
        const result = await apiPost('/simulation/validate', config);
        statusDiv.classList.remove('hidden');

        if (result.valid) {
            statusDiv.className = 'mb-6 p-4 rounded border bg-green-50 border-green-200 dark:bg-green-900/30 dark:border-green-700';
            statusDiv.innerHTML = `
                <div class="text-green-800 dark:text-green-400 font-medium">✓ Configuration Valid</div>
                <div class="text-sm text-green-600 dark:text-green-300">
                    ${result.layer_count} layers, ${result.step_count} steps, ${result.total_duration_days?.toFixed(1) || 'N/A'} days total
                </div>
                ${result.warnings?.length > 0 ? `<div class="mt-2 text-yellow-700 dark:text-yellow-400 text-sm">⚠️ ${result.warnings.join(', ')}</div>` : ''}
            `;
        } else {
            statusDiv.className = 'mb-6 p-4 rounded border bg-red-50 border-red-200 dark:bg-red-900/30 dark:border-red-700';
            statusDiv.innerHTML = `
                <div class="text-red-800 dark:text-red-400 font-medium">✗ Validation Failed</div>
                <ul class="text-sm text-red-600 dark:text-red-300 list-disc list-inside">
                    ${result.errors.map(e => `<li>${e}</li>`).join('')}
                </ul>
            `;
        }
    } catch (e) {
        console.error('Validation failed:', e);
        statusDiv.classList.remove('hidden');
        statusDiv.className = 'mb-6 p-4 rounded border bg-red-50 border-red-200 dark:bg-red-900/30 dark:border-red-700';
        statusDiv.innerHTML = '<div class="text-red-800 dark:text-red-400">Validation request failed</div>';
    }
}

function buildSimulationConfig() {
    // Build substances with their layer assignments
    const substances = state.substances.map(sub => {
        const subId = sub.id || sub.cid || sub.name;
        const layer_assignments = {};

        // Get C0 and D/k for each layer from c0Matrix and substance data
        state.layers.forEach(layer => {
            const c0Data = state.c0Matrix[subId]?.[layer.index];
            if (c0Data && c0Data.value > 0) {
                layer_assignments[layer.index] = {
                    C0: c0Data.value,
                    C0_unit: c0Data.unit || 'mg/kg',
                    D: sub[`D_L${layer.index}`] || sub.D_computed,
                    k: sub[`k_L${layer.index}`] || sub.k_computed,
                };
            }
        });

        // Get CF0 for this substance
        const cf0Data = state.CF0[subId];

        // Get CAS - handle array (common from PubChem)
        let cas = sub.cas;
        if (Array.isArray(cas)) {
            cas = cas[0];  // Use first CAS number
        }

        return {
            id: subId,
            name: sub.name,
            cas: cas || null,                    // CAS for fallback lookup
            cid: sub.cid || null,                // PubChem CID for fallback lookup
            mw: sub.mw || null,                  // Molecular weight
            logP: sub.logP || null,              // Log partition coefficient
            SML: sub.SML || sub.regulatory?.EU?.SML || null,
            CF0: cf0Data?.value || 0,            // Initial concentration in food
            CF0_unit: cf0Data?.unit || 'mg/kg',
            layer_assignments: layer_assignments,
        };
    }).filter(sub => Object.keys(sub.layer_assignments).length > 0);

    return {
        name: document.getElementById('sim-name').value || 'Simulation',
        layers: state.layers.map(l => ({
            index: l.index,
            polymer: l.polymer,
            thickness: l.thickness,
            thickness_unit: l.thickness_unit,
            C0: l.C0 || 0,
        })),
        steps: state.steps.map(s => ({
            index: s.index,
            temperature_C: s.temperature_C,
            duration: s.duration,
            duration_unit: s.duration_unit,
            with_food: s.with_food,
            setoff_type: s.setoff_type || 'stacked',
        })),
        geometry: state.geometry,
        food: {
            ...state.food,
            // Include CF0 per substance for simulation
            CF0: Object.keys(state.CF0).length > 0 ? state.CF0 : undefined,
        },
        substances: substances.length > 0 ? substances : null,
    };
}

async function runSimulation() {
    const progressDiv = document.getElementById('simulation-progress');
    const progressBar = document.getElementById('progress-bar');
    const progressText = document.getElementById('progress-text');

    progressDiv.classList.remove('hidden');
    progressBar.style.width = '10%';
    progressText.textContent = 'Starting simulation...';

    const config = buildSimulationConfig();

    try {
        const result = await apiPost('/simulation/run?async_mode=false', config);
        state.currentJobId = result.job_id;

        progressBar.style.width = '100%';

        if (result.status === 'failed') {
            // Simulation failed - show error
            progressText.textContent = 'Simulation failed';
            progressBar.classList.add('bg-red-600');
            progressBar.classList.remove('bg-sfppy-green');
            showToast(result.error || 'Simulation failed', 'error');
            console.error('Simulation failed:', result.error);
            return;
        }

        progressText.textContent = result.status === 'completed' ? 'Completed!' : `Status: ${result.status}`;

        if (result.results) {
            displayResults(result.results);
            // Switch to results tab
            document.querySelector('.tab-btn[data-tab="results"]').click();
        }
    } catch (e) {
        console.error('Simulation failed:', e);
        progressText.textContent = 'Simulation failed';
        progressBar.classList.add('bg-red-600');
        progressBar.classList.remove('bg-sfppy-green');
        showToast(e.message || 'Network error', 'error');
    }
}

// ========== RESULTS TAB ==========

// Chart configuration state
const chartConfig = {
    timeScale: 'linear',  // 'linear' or 'sqrt'
    concScale: 'linear',  // 'linear' or 'log'
    showSML: true,
    showTcontact: true,
    showGrid: true,
};

// Store current results for chart updates
let currentResults = null;

function displayResults(results) {
    currentResults = results;
    document.getElementById('no-results').classList.add('hidden');
    document.getElementById('results-content').classList.remove('hidden');

    // Contact time info
    const tcontactDays = results.tcontact_days || results.time_days?.[results.time_days.length - 1] / 2;

    // Compliance summary - handle multiple substances
    const complianceDiv = document.getElementById('compliance-summary');
    const substances = results.substances || [];

    if (substances.length > 1) {
        // Multiple substances - show detailed table
        const allCompliant = results.compliant;
        const complianceClass = allCompliant
            ? 'mb-6 p-4 rounded-lg bg-green-100 dark:bg-green-900/30 border border-green-300 dark:border-green-700'
            : 'mb-6 p-4 rounded-lg bg-red-100 dark:bg-red-900/30 border border-red-300 dark:border-red-700';

        complianceDiv.className = complianceClass;
        complianceDiv.innerHTML = `
            <div class="${allCompliant ? 'text-green-800 dark:text-green-400' : 'text-red-800 dark:text-red-400'} font-bold text-lg mb-3">
                ${allCompliant ? '✓ ALL SUBSTANCES COMPLIANT' : '✗ NON-COMPLIANT'}
            </div>
            <div class="text-sm text-gray-600 dark:text-gray-300 mb-2">
                Contact time: ${tcontactDays?.toFixed(1) || '--'} days | Simulation: ${(tcontactDays * 2)?.toFixed(1) || '--'} days (2×tcontact)
            </div>
            <table class="w-full text-sm mt-2">
                <thead>
                    <tr class="border-b dark:border-gray-600">
                        <th class="text-left py-1">Substance</th>
                        <th class="text-right py-1">CF @ t<sub>contact</sub></th>
                        <th class="text-right py-1">SML</th>
                        <th class="text-right py-1">Margin</th>
                        <th class="text-center py-1">Status</th>
                    </tr>
                </thead>
                <tbody>
                    ${substances.map(sub => `
                        <tr class="border-b dark:border-gray-700">
                            <td class="py-1 flex items-center gap-2">
                                <span style="background:${sub.color}; width:12px; height:12px; border-radius:2px; display:inline-block;"></span>
                                ${sub.name}
                            </td>
                            <td class="text-right py-1">${sub.CF_at_tcontact?.toFixed(3) || '--'} mg/kg</td>
                            <td class="text-right py-1">${sub.SML_mg_kg?.toFixed(1) || 'N/A'} mg/kg</td>
                            <td class="text-right py-1">${sub.margin_percent?.toFixed(1) || '--'}%</td>
                            <td class="text-center py-1">
                                ${sub.compliant
                                    ? '<span class="text-green-600 dark:text-green-400">✓</span>'
                                    : '<span class="text-red-600 dark:text-red-400">✗</span>'}
                            </td>
                        </tr>
                    `).join('')}
                </tbody>
            </table>
        `;
    } else {
        // Single substance (backward compatibility)
        if (results.compliant !== undefined) {
            if (results.compliant) {
                complianceDiv.className = 'mb-6 p-4 rounded-lg bg-green-100 dark:bg-green-900/30 border border-green-300 dark:border-green-700';
                complianceDiv.innerHTML = `
                    <div class="text-green-800 dark:text-green-400 font-bold text-lg">✓ COMPLIANT</div>
                    <div class="text-green-700 dark:text-green-300 text-sm">Migration below SML limit</div>
                    <div class="text-sm text-gray-600 dark:text-gray-400 mt-1">
                        Contact time: ${tcontactDays?.toFixed(1) || '--'} days
                    </div>
                `;
            } else {
                complianceDiv.className = 'mb-6 p-4 rounded-lg bg-red-100 dark:bg-red-900/30 border border-red-300 dark:border-red-700';
                complianceDiv.innerHTML = `
                    <div class="text-red-800 dark:text-red-400 font-bold text-lg">✗ NON-COMPLIANT</div>
                    <div class="text-red-700 dark:text-red-300 text-sm">Migration exceeds SML limit</div>
                    <div class="text-sm text-gray-600 dark:text-gray-400 mt-1">
                        Contact time: ${tcontactDays?.toFixed(1) || '--'} days
                    </div>
                `;
            }
        }
    }

    // Key values - generate cards per substance
    displayKeyValueCards(results);

    // Initialize chart options if not already done
    initChartOptions();

    // Draw charts
    drawMigrationChart(results);
    drawConcentrationProfileChart(results);
}

/**
 * Display key value cards for each substance with color-coded borders
 */
function displayKeyValueCards(results) {
    const container = document.getElementById('results-key-values-container');
    if (!container) return;

    const substances = results.substances || [];

    // If no substances, use backward-compatible single-substance display
    if (substances.length === 0) {
        const cfAtTcontact = results.final_CF_mg_kg ?? results.CF_mg_kg?.[results.CF_mg_kg?.length - 1] ?? 0;
        const eqCF = results.equilibrium_CF_mg_kg ?? cfAtTcontact;
        const eqPercent = eqCF && cfAtTcontact ? ((cfAtTcontact / eqCF) * 100).toFixed(1) : '--';

        container.innerHTML = `
            <div class="grid grid-cols-4 gap-2">
                <div class="bg-gray-50 dark:bg-gray-700 rounded-lg p-2 text-center">
                    <p class="text-xs text-gray-600 dark:text-gray-400">CF @ t<sub>contact</sub></p>
                    <p class="text-lg font-bold text-gray-800 dark:text-white">${cfAtTcontact?.toFixed(3) || '--'}</p>
                    <p class="text-xs text-gray-500">mg/kg</p>
                </div>
                <div class="bg-gray-50 dark:bg-gray-700 rounded-lg p-2 text-center">
                    <p class="text-xs text-gray-600 dark:text-gray-400">SML</p>
                    <p class="text-lg font-bold text-gray-800 dark:text-white">${results.SML_mg_kg?.toFixed(1) || '--'}</p>
                    <p class="text-xs text-gray-500">mg/kg</p>
                </div>
                <div class="bg-gray-50 dark:bg-gray-700 rounded-lg p-2 text-center">
                    <p class="text-xs text-gray-600 dark:text-gray-400">Margin</p>
                    <p class="text-lg font-bold text-gray-800 dark:text-white">${results.margin_percent?.toFixed(1) || '--'}</p>
                    <p class="text-xs text-gray-500">%</p>
                </div>
                <div class="bg-gray-50 dark:bg-gray-700 rounded-lg p-2 text-center">
                    <p class="text-xs text-gray-600 dark:text-gray-400">Equilibrium</p>
                    <p class="text-lg font-bold text-gray-800 dark:text-white">${eqPercent}</p>
                    <p class="text-xs text-gray-500">%</p>
                </div>
            </div>
        `;
        return;
    }

    // Multiple substances - generate one row of cards per substance
    let html = '<div class="space-y-3">';

    substances.forEach((sub, i) => {
        const cfAtTcontact = sub.CF_at_tcontact ?? sub.CF_mg_kg?.[sub.CF_mg_kg?.length - 1] ?? 0;
        const eqCF = sub.equilibrium_CF_mg_kg ?? cfAtTcontact;
        const eqPercent = eqCF && cfAtTcontact ? ((cfAtTcontact / eqCF) * 100).toFixed(1) : '--';
        const color = sub.color || `hsl(${i * 137}, 70%, 50%)`;
        const isCompliant = sub.compliant !== false;
        const bgClass = isCompliant ? 'bg-gray-50 dark:bg-gray-700' : 'bg-red-50 dark:bg-red-900/30';

        html += `
            <div class="rounded-lg overflow-hidden" style="border-left: 4px solid ${color};">
                <div class="flex items-center gap-2 px-3 py-1.5 bg-gray-100 dark:bg-gray-700/50">
                    <span style="background:${color}; width:10px; height:10px; border-radius:2px; display:inline-block;"></span>
                    <span class="text-sm font-medium text-gray-800 dark:text-white">${sub.name}</span>
                    <span class="ml-auto text-xs ${isCompliant ? 'text-green-600 dark:text-green-400' : 'text-red-600 dark:text-red-400'}">
                        ${isCompliant ? '✓ Compliant' : '✗ Exceeds SML'}
                    </span>
                </div>
                <div class="grid grid-cols-4 gap-2 p-2 ${bgClass}">
                    <div class="text-center">
                        <p class="text-xs text-gray-600 dark:text-gray-400">CF @ t<sub>c</sub></p>
                        <p class="text-base font-bold" style="color: ${color}">${cfAtTcontact?.toFixed(3) || '--'}</p>
                        <p class="text-xs text-gray-500">mg/kg</p>
                    </div>
                    <div class="text-center">
                        <p class="text-xs text-gray-600 dark:text-gray-400">SML</p>
                        <p class="text-base font-bold text-gray-800 dark:text-white">${sub.SML_mg_kg?.toFixed(1) || 'N/A'}</p>
                        <p class="text-xs text-gray-500">mg/kg</p>
                    </div>
                    <div class="text-center">
                        <p class="text-xs text-gray-600 dark:text-gray-400">Margin</p>
                        <p class="text-base font-bold ${sub.margin_percent < 20 ? 'text-orange-600' : 'text-gray-800 dark:text-white'}">${sub.margin_percent?.toFixed(1) || '--'}</p>
                        <p class="text-xs text-gray-500">%</p>
                    </div>
                    <div class="text-center">
                        <p class="text-xs text-gray-600 dark:text-gray-400">Equilibrium</p>
                        <p class="text-base font-bold text-gray-800 dark:text-white">${eqPercent}</p>
                        <p class="text-xs text-gray-500">%</p>
                    </div>
                </div>
            </div>
        `;
    });

    html += '</div>';
    container.innerHTML = html;
}

// Store original dimensions for fullscreen restoration
const chartOriginalDimensions = {};

function initChartOptions() {
    // Check if already initialized
    if (document.getElementById('chart-time-scale')) return;

    // Add chart controls to the container in the HTML
    const controlsContainer = document.getElementById('chart-controls-container');
    if (!controlsContainer) return;

    const controlsHtml = `
        <div id="chart-controls" class="flex flex-wrap items-center gap-3 text-sm">
            <div class="flex items-center gap-1">
                <label class="text-gray-600 dark:text-gray-300">Time:</label>
                <select id="chart-time-scale" class="border dark:border-gray-600 rounded px-2 py-1 text-xs bg-white dark:bg-gray-700">
                    <option value="linear">Linear</option>
                    <option value="sqrt">√t</option>
                </select>
            </div>
            <div class="flex items-center gap-1">
                <label class="text-gray-600 dark:text-gray-300">Conc:</label>
                <select id="chart-conc-scale" class="border dark:border-gray-600 rounded px-2 py-1 text-xs bg-white dark:bg-gray-700">
                    <option value="linear">Linear</option>
                    <option value="log">Log</option>
                </select>
            </div>
            <label class="flex items-center gap-1 text-gray-600 dark:text-gray-300">
                <input type="checkbox" id="chart-show-sml" checked class="rounded w-3 h-3">
                SML
            </label>
            <label class="flex items-center gap-1 text-gray-600 dark:text-gray-300">
                <input type="checkbox" id="chart-show-tcontact" checked class="rounded w-3 h-3">
                t<sub>c</sub>
            </label>
            <label class="flex items-center gap-1 text-gray-600 dark:text-gray-300">
                <input type="checkbox" id="chart-show-grid" checked class="rounded w-3 h-3">
                Grid
            </label>
        </div>
    `;

    controlsContainer.innerHTML = controlsHtml;

    // Add event listeners
    document.getElementById('chart-time-scale').addEventListener('change', (e) => {
        chartConfig.timeScale = e.target.value;
        if (currentResults) drawMigrationChart(currentResults);
    });

    document.getElementById('chart-conc-scale').addEventListener('change', (e) => {
        chartConfig.concScale = e.target.value;
        if (currentResults) drawMigrationChart(currentResults);
    });

    document.getElementById('chart-show-sml').addEventListener('change', (e) => {
        chartConfig.showSML = e.target.checked;
        if (currentResults) drawMigrationChart(currentResults);
    });

    document.getElementById('chart-show-tcontact').addEventListener('change', (e) => {
        chartConfig.showTcontact = e.target.checked;
        if (currentResults) drawMigrationChart(currentResults);
    });

    document.getElementById('chart-show-grid')?.addEventListener('change', (e) => {
        chartConfig.showGrid = e.target.checked;
        if (currentResults) drawMigrationChart(currentResults);
    });

    document.getElementById('chart-fullscreen')?.addEventListener('click', () => {
        toggleChartFullscreen('migration-chart-wrapper', 'migration');
    });
}

// Zoom control for charts
function zoomChart(chartType, action) {
    let chart;
    if (chartType === 'migration') {
        chart = state.migrationChart;
    } else if (chartType === 'profile') {
        chart = state.profileChart;
    } else if (chartType === 'fitting-data') {
        chart = fittingState.dataChart;
    } else if (chartType === 'fitting-result') {
        chart = fittingState.resultChart;
    }

    if (!chart) return;

    switch (action) {
        case 'in':
            chart.zoom(1.2);
            break;
        case 'out':
            chart.zoom(0.8);
            break;
        case 'reset':
            chart.resetZoom();
            break;
    }
}

function toggleChartFullscreen(wrapperId, chartType) {
    const wrapper = document.getElementById(wrapperId);
    if (!wrapper) return;

    const chartContainer = wrapper.querySelector('.chart-container');
    const isFullscreen = wrapper.classList.contains('chart-fullscreen-wrapper');

    if (isFullscreen) {
        // Exit fullscreen - restore original dimensions
        wrapper.classList.remove('chart-fullscreen-wrapper');
        document.body.style.overflow = '';

        // Restore original height
        if (chartOriginalDimensions[wrapperId] && chartContainer) {
            chartContainer.style.height = chartOriginalDimensions[wrapperId].height;
        }

        // Force chart resize after DOM update
        setTimeout(() => {
            if (chartType === 'migration' && state.migrationChart) {
                state.migrationChart.resize();
            } else if (chartType === 'profile' && state.profileChart) {
                state.profileChart.resize();
            }
        }, 50);

    } else {
        // Store original dimensions before going fullscreen
        if (chartContainer) {
            chartOriginalDimensions[wrapperId] = {
                height: chartContainer.style.height || '280px'
            };
        }

        // Enter fullscreen
        wrapper.classList.add('chart-fullscreen-wrapper');
        document.body.style.overflow = 'hidden';

        // Force chart resize after DOM update
        setTimeout(() => {
            if (chartType === 'migration' && state.migrationChart) {
                state.migrationChart.resize();
            } else if (chartType === 'profile' && state.profileChart) {
                state.profileChart.resize();
            }
        }, 50);

        // Add escape key listener
        const escHandler = (e) => {
            if (e.key === 'Escape') {
                toggleChartFullscreen(wrapperId, chartType);
                document.removeEventListener('keydown', escHandler);
            }
        };
        document.addEventListener('keydown', escHandler);
    }

    // Update fullscreen button text
    const fsBtn = wrapper.querySelector('#chart-fullscreen, #profile-fullscreen');
    if (fsBtn) {
        fsBtn.textContent = isFullscreen ? '⛶' : '✕';
        fsBtn.title = isFullscreen ? 'Toggle Fullscreen' : 'Exit Fullscreen (Esc)';
    }
}

function drawMigrationChart(results) {
    const ctx = document.getElementById('migration-chart').getContext('2d');

    if (state.migrationChart) {
        state.migrationChart.destroy();
    }

    const substances = results.substances || [];
    const timeDays = results.time_days || [];
    const tcontactDays = results.tcontact_days || timeDays[Math.floor(timeDays.length / 2)];

    // Transform time if sqrt scale selected
    let xData = [...timeDays];
    let xLabel = 'Time (days)';

    if (chartConfig.timeScale === 'sqrt') {
        xData = timeDays.map(t => Math.sqrt(t));
        xLabel = '√Time (√days)';
    }

    // Build datasets and collect CF values for legend
    const datasets = [];
    const legendData = [];  // For custom legend with CF values

    if (substances.length > 0) {
        // Multiple substances - each with its own color and SML line
        substances.forEach((sub, i) => {
            const cfAtTcontact = sub.CF_at_tcontact ?? sub.CF_mg_kg?.[sub.CF_mg_kg.length - 1] ?? 0;
            const cfEquilibrium = sub.equilibrium_CF_mg_kg ?? cfAtTcontact;

            legendData.push({
                name: sub.name,
                color: sub.color,
                cfTcontact: cfAtTcontact,
                cfEquilibrium: cfEquilibrium,
                sml: sub.SML_mg_kg,
            });

            // Kinetics curve
            datasets.push({
                label: sub.name,
                data: xData.map((t, j) => ({
                    x: t,
                    y: sub.CF_mg_kg[j] || 0,
                })),
                borderColor: sub.color,
                backgroundColor: sub.color + '20',
                fill: false,
                tension: 0.3,
                pointRadius: 0,
                borderWidth: 2,
            });

            // SML line for this substance (same color, dashed)
            if (chartConfig.showSML && sub.SML_mg_kg) {
                datasets.push({
                    label: `SML (${sub.name})`,
                    data: xData.map(t => ({ x: t, y: sub.SML_mg_kg })),
                    borderColor: sub.color,
                    borderDash: [5, 5],
                    pointRadius: 0,
                    fill: false,
                    borderWidth: 1.5,
                });
            }
        });
    } else {
        // Single substance (backward compatibility)
        const cfAtTcontact = results.final_CF_mg_kg ?? results.CF_mg_kg?.[results.CF_mg_kg.length - 1] ?? 0;
        const cfEquilibrium = results.equilibrium_CF_mg_kg ?? cfAtTcontact;

        legendData.push({
            name: 'Substance',
            color: 'rgb(59, 130, 246)',
            cfTcontact: cfAtTcontact,
            cfEquilibrium: cfEquilibrium,
            sml: results.SML_mg_kg,
        });

        datasets.push({
            label: 'CF (mg/kg)',
            data: xData.map((t, i) => ({
                x: t,
                y: results.CF_mg_kg?.[i] || 0,
            })),
            borderColor: 'rgb(59, 130, 246)',
            backgroundColor: 'rgba(59, 130, 246, 0.1)',
            fill: true,
            tension: 0.3,
            pointRadius: 0,
            borderWidth: 2,
        });

        if (chartConfig.showSML && results.SML_mg_kg) {
            datasets.push({
                label: 'SML',
                data: xData.map(t => ({ x: t, y: results.SML_mg_kg })),
                borderColor: 'rgb(34, 197, 94)',
                borderDash: [5, 5],
                pointRadius: 0,
                fill: false,
                borderWidth: 1.5,
            });
        }
    }

    // Configure scales with grid
    const gridConfig = chartConfig.showGrid ? {
        color: 'rgba(0, 0, 0, 0.1)',
        drawBorder: true,
    } : {
        display: false,
    };

    const yScale = {
        title: { display: true, text: 'Concentration (mg/kg food)' },
        beginAtZero: chartConfig.concScale !== 'log',
        grid: gridConfig,
    };

    if (chartConfig.concScale === 'log') {
        yScale.type = 'logarithmic';
        yScale.min = 0.001;
    }

    // Build annotations for tcontact vertical line and CF values
    const annotations = {};
    const tcontactX = chartConfig.timeScale === 'sqrt' ? Math.sqrt(tcontactDays) : tcontactDays;

    // Add step boundary lines for multi-step simulations
    const stepsSummary = results.steps_summary || [];
    if (results.multi_step && stepsSummary.length > 1) {
        let cumulativeTime = 0;
        const stepColors = ['#3B82F6', '#10B981', '#F59E0B', '#EF4444', '#8B5CF6'];

        // Convert duration to days based on unit
        const toDays = (value, unit) => {
            const conversions = {
                'minutes': 1 / 1440, 'hours': 1 / 24, 'days': 1, 'months': 30, 'weeks': 7
            };
            return value * (conversions[unit] || 1);
        };

        stepsSummary.forEach((step, i) => {
            const stepDuration = toDays(step.duration, step.duration_unit);
            const stepEndTime = cumulativeTime + stepDuration;
            const stepEndX = chartConfig.timeScale === 'sqrt' ? Math.sqrt(stepEndTime) : stepEndTime;

            // Add step boundary line (except for last step)
            if (i < stepsSummary.length - 1) {
                annotations[`stepBoundary_${i}`] = {
                    type: 'line',
                    xMin: stepEndX,
                    xMax: stepEndX,
                    borderColor: step.with_food ? 'rgba(34, 197, 94, 0.5)' : 'rgba(245, 158, 11, 0.5)',
                    borderWidth: 1.5,
                    borderDash: [3, 3],
                };
                annotations[`stepLabel_${i}`] = {
                    type: 'label',
                    xValue: stepEndX,
                    yValue: 'max',
                    yAdjust: -10 - (i % 2) * 15,  // Alternate positions to avoid overlap
                    content: `Step ${step.index} end`,
                    backgroundColor: step.with_food ? 'rgba(34, 197, 94, 0.1)' : 'rgba(245, 158, 11, 0.1)',
                    borderColor: step.with_food ? 'rgba(34, 197, 94, 0.5)' : 'rgba(245, 158, 11, 0.5)',
                    borderWidth: 1,
                    borderRadius: 3,
                    font: { size: 9 },
                    color: step.with_food ? '#166534' : '#92400e',
                    padding: 2,
                    position: 'start',
                };
            }

            cumulativeTime = stepEndTime;
        });
    }

    if (chartConfig.showTcontact && tcontactDays) {
        annotations.tcontactLine = {
            type: 'line',
            xMin: tcontactX,
            xMax: tcontactX,
            borderColor: 'rgba(100, 100, 100, 0.6)',
            borderWidth: 2,
            borderDash: [4, 4],
            label: {
                display: true,
                content: `tc = ${tcontactDays.toFixed(1)} d`,
                position: 'start',
                backgroundColor: 'rgba(100, 100, 100, 0.8)',
                font: { size: 10 },
            },
        };

        // Add point annotations for CF at tcontact for each substance
        legendData.forEach((sub, i) => {
            if (sub.cfTcontact !== undefined) {
                annotations[`cfPoint_${i}`] = {
                    type: 'point',
                    xValue: tcontactX,
                    yValue: sub.cfTcontact,
                    backgroundColor: sub.color,
                    borderColor: 'white',
                    borderWidth: 2,
                    radius: 6,
                };
                annotations[`cfLabel_${i}`] = {
                    type: 'label',
                    xValue: tcontactX,
                    yValue: sub.cfTcontact,
                    content: `${sub.cfTcontact.toFixed(1)}`,
                    position: { x: 'end', y: 'center' },
                    xAdjust: 8,
                    backgroundColor: 'rgba(255, 255, 255, 0.9)',
                    borderColor: sub.color,
                    borderWidth: 1,
                    borderRadius: 4,
                    font: { size: 10, weight: 'bold' },
                    color: sub.color,
                    padding: 3,
                };
            }
        });
    }

    state.migrationChart = new Chart(ctx, {
        type: 'line',
        data: { datasets },
        options: {
            responsive: true,
            maintainAspectRatio: false,
            parsing: false,
            interaction: {
                mode: 'nearest',
                intersect: false,
            },
            scales: {
                x: {
                    type: 'linear',
                    title: { display: true, text: xLabel },
                    grid: gridConfig,
                },
                y: yScale,
            },
            plugins: {
                legend: {
                    display: false,  // We use custom legend
                },
                annotation: {
                    annotations: annotations,
                },
                tooltip: {
                    callbacks: {
                        label: (ctx) => {
                            const label = ctx.dataset.label || '';
                            const yVal = ctx.parsed.y;
                            return `${label}: ${yVal.toFixed(2)} mg/kg`;
                        },
                    },
                },
                zoom: {
                    pan: {
                        enabled: true,
                        mode: 'xy',
                    },
                    zoom: {
                        wheel: {
                            enabled: true,
                        },
                        pinch: {
                            enabled: true,
                        },
                        mode: 'xy',
                    },
                },
            },
        },
    });

    // Update custom legend with CF values
    updateMigrationLegend(legendData, tcontactDays);
}

function updateMigrationLegend(legendData, tcontactDays) {
    const legendDiv = document.getElementById('migration-chart-legend');
    if (!legendDiv) return;

    legendDiv.innerHTML = legendData.map(sub => {
        const eqPercent = sub.cfEquilibrium > 0 ? ((sub.cfTcontact / sub.cfEquilibrium) * 100).toFixed(0) : '--';
        const smlText = sub.sml ? `SML: ${sub.sml} mg/kg` : '';
        const complianceIcon = sub.sml ? (sub.cfTcontact < sub.sml ? '✓' : '✗') : '';
        const complianceColor = sub.sml ? (sub.cfTcontact < sub.sml ? 'text-green-600' : 'text-red-600') : '';

        return `
            <div class="flex items-center gap-2 bg-gray-50 dark:bg-gray-700 px-2 py-1 rounded">
                <span class="w-3 h-3 rounded-full" style="background: ${sub.color}"></span>
                <span class="font-medium">${sub.name}</span>
                <span>CF(tc): <b>${sub.cfTcontact.toFixed(2)}</b></span>
                <span>CF(∞): ${sub.cfEquilibrium.toFixed(2)}</span>
                <span>(${eqPercent}%)</span>
                ${smlText ? `<span class="${complianceColor}">${complianceIcon} ${smlText}</span>` : ''}
            </div>
        `;
    }).join('');
}

function drawConcentrationProfileChart(results) {
    const profileContainer = document.getElementById('concentration-profile-container');
    if (!profileContainer) return;

    // Destroy existing chart
    if (state.profileChart) {
        state.profileChart.destroy();
        state.profileChart = null;
    }

    // Check for multi-substance profiles
    const substances = results.substances || [];
    const hasMultiSubstanceProfiles = substances.length > 1 && substances.some(s => s.concentration_profile?.Cx_mg_kg?.length > 0);

    // Get x data from global profiles or first substance
    const profiles = results.concentration_profiles;
    const xData = profiles?.x_um || profiles?.x_normalized;
    const xLabel = profiles?.x_um ? 'Position (µm)' : 'Position (normalized)';
    const hasRealX = !!profiles?.x_um;

    // Check if we have any profile data
    const hasGlobalProfile = profiles && xData && profiles.Cx_mg_kg?.length > 0;
    const hasAnyProfile = hasGlobalProfile || hasMultiSubstanceProfiles;

    if (!hasAnyProfile) {
        profileContainer.innerHTML = '<div class="text-gray-500 dark:text-gray-400 text-sm p-4 text-center">No concentration profile data available</div>';
        return;
    }

    // Build substance selector if multiple substances
    let substanceSelector = '';
    if (hasMultiSubstanceProfiles) {
        substanceSelector = `
            <select id="profile-substance-select" class="border dark:border-gray-600 rounded px-2 py-1 text-xs bg-white dark:bg-gray-700 mr-2">
                ${substances.map((s, i) => `<option value="${i}" style="color:${s.color}">${s.name}</option>`).join('')}
            </select>
        `;
    }

    // Build step selector for multi-step simulations
    let stepSelector = '';
    const isMultiStep = results.multi_step || (results.steps_summary && results.steps_summary.length > 1);
    const stepsSummary = results.steps_summary || [];
    if (isMultiStep && stepsSummary.length > 1) {
        stepSelector = `
            <select id="profile-step-select" class="border dark:border-gray-600 rounded px-2 py-1 text-xs bg-white dark:bg-gray-700 mr-2" title="Filter profiles by step">
                <option value="all">All Steps</option>
                ${stepsSummary.map(s => `<option value="${s.index}">Step ${s.index}: ${s.temperature_C}°C ${s.with_food ? '' : '(set-off)'}</option>`).join('')}
            </select>
        `;
    }

    profileContainer.innerHTML = `
        <div id="profile-chart-wrapper">
            <div class="chart-controls-bar flex items-center justify-between mb-2">
                <div>
                    <h4 class="font-medium text-gray-800 dark:text-gray-200 text-sm">Concentration Profiles C(x) in Packaging</h4>
                    <p class="text-xs text-gray-500 dark:text-gray-400">x = 0: food interface | x = max: back of layer</p>
                </div>
                <div class="flex items-center flex-wrap gap-1">
                    ${substanceSelector}
                    ${stepSelector}
                    <label class="flex items-center gap-1 text-xs text-gray-600 dark:text-gray-300 mr-2">
                        <input type="checkbox" id="profile-show-grid" checked class="rounded w-3 h-3">
                        Grid
                    </label>
                    <div class="chart-zoom-controls">
                        <button onclick="zoomChart('profile', 'in')" title="Zoom In">+</button>
                        <button onclick="zoomChart('profile', 'out')" title="Zoom Out">-</button>
                        <button onclick="zoomChart('profile', 'reset')" title="Reset Zoom">↺</button>
                    </div>
                    <button id="profile-fullscreen" class="ml-2 px-2 py-1 border dark:border-gray-600 rounded text-sm hover:bg-gray-100 dark:hover:bg-gray-700" title="Toggle Fullscreen">
                        ⛶
                    </button>
                </div>
            </div>
            ${isMultiStep ? `
            <div class="step-legend mb-2 flex flex-wrap gap-2 text-xs">
                ${stepsSummary.map(s => `
                    <span class="px-2 py-1 rounded ${s.with_food ? 'bg-green-100 text-green-800' : 'bg-amber-100 text-amber-800'}">
                        Step ${s.index}: ${s.temperature_C}°C, ${s.duration}${s.duration_unit?.charAt(0) || 'd'}
                        ${s.with_food ? '' : ` (${s.setoff_type || 'set-off'})`}
                    </span>
                `).join('')}
            </div>
            ` : ''}
            <div class="chart-container" style="height: 300px; position: relative;">
                <canvas id="profile-chart"></canvas>
            </div>
        </div>
    `;

    document.getElementById('profile-fullscreen').addEventListener('click', () => {
        toggleChartFullscreen('profile-chart-wrapper', 'profile');
    });

    // Store results for substance switching
    state.profileResults = results;
    state.selectedProfileStep = 'all';  // Track selected step filter

    // Add substance selector listener
    const substanceSelect = document.getElementById('profile-substance-select');
    if (substanceSelect) {
        substanceSelect.addEventListener('change', (e) => {
            renderProfileChartForSubstance(parseInt(e.target.value), state.selectedProfileStep);
        });
    }

    // Add step selector listener
    const stepSelect = document.getElementById('profile-step-select');
    if (stepSelect) {
        stepSelect.addEventListener('change', (e) => {
            state.selectedProfileStep = e.target.value;
            const substanceIdx = substanceSelect ? parseInt(substanceSelect.value) : 0;
            renderProfileChartForSubstance(substanceIdx, e.target.value);
        });
    }

    // Render initial chart (first substance or global profile)
    renderProfileChartForSubstance(0, 'all');
}

function renderProfileChartForSubstance(substanceIndex, stepFilter = 'all') {
    const results = state.profileResults;
    if (!results) return;

    const ctx = document.getElementById('profile-chart')?.getContext('2d');
    if (!ctx) return;

    // Destroy existing chart
    if (state.profileChart) {
        state.profileChart.destroy();
        state.profileChart = null;
    }

    const substances = results.substances || [];
    const globalProfiles = results.concentration_profiles;
    const xData = globalProfiles?.x_um || globalProfiles?.x_normalized || [];
    const xLabel = globalProfiles?.x_um ? 'Position (µm)' : 'Position (normalized)';
    const hasRealX = !!globalProfiles?.x_um;

    // Get profile data for selected substance or use global
    let profileData, substanceColor, substanceName;
    const substance = substances[substanceIndex];

    // Check if we need to filter by step (for multi-step simulations)
    if (stepFilter !== 'all' && substance?.step_profiles?.length > 0) {
        // Filter by specific step
        const stepIdx = parseInt(stepFilter);
        const stepProfile = substance.step_profiles.find(sp => sp.step_index === stepIdx);
        if (stepProfile) {
            profileData = {
                times_days: stepProfile.times_days,
                Cx_mg_kg: stepProfile.Cx_mg_kg,
            };
            substanceColor = substance.color;
            substanceName = `${substance.name} - Step ${stepIdx} (${stepProfile.temperature_C}°C${stepProfile.step_type?.includes('setoff') ? ', set-off' : ''})`;
        }
    } else if (substance?.concentration_profile?.Cx_mg_kg?.length > 0) {
        profileData = substance.concentration_profile;
        substanceColor = substance.color;
        substanceName = substance.name;
    } else {
        profileData = globalProfiles;
        substanceColor = substances[0]?.color || 'rgb(59, 130, 246)';
        substanceName = substances[0]?.name || 'Substance';
    }

    if (!profileData?.Cx_mg_kg?.length || !xData.length) {
        return;
    }

    // Color gradient for time progression
    const nProfiles = profileData.times_days?.length || 0;
    const baseColor = substanceColor || '#3B82F6';

    // Parse base color to RGB (handles both hex and rgb formats)
    let baseR = 59, baseG = 130, baseB = 246;  // Default blue
    if (baseColor.startsWith('#')) {
        // Hex format: #RRGGBB
        const hex = baseColor.slice(1);
        baseR = parseInt(hex.slice(0, 2), 16);
        baseG = parseInt(hex.slice(2, 4), 16);
        baseB = parseInt(hex.slice(4, 6), 16);
    } else if (baseColor.startsWith('rgb')) {
        // RGB format: rgb(r, g, b)
        const rgbMatch = baseColor.match(/\d+/g);
        if (rgbMatch) {
            baseR = parseInt(rgbMatch[0]);
            baseG = parseInt(rgbMatch[1]);
            baseB = parseInt(rgbMatch[2]);
        }
    }

    // Different color gradients for different steps when showing all
    const stepColors = ['#3B82F6', '#10B981', '#F59E0B', '#EF4444', '#8B5CF6'];

    const colorScale = (i, stepIdx = 0) => {
        const ratio = nProfiles > 1 ? i / (nProfiles - 1) : 0;
        // Lighten color for earlier times, darken for later
        const factor = 0.4 + ratio * 0.6;
        const r = Math.round(baseR * factor + (1 - factor) * 200);
        const g = Math.round(baseG * factor + (1 - factor) * 200);
        const b = Math.round(baseB * factor + (1 - factor) * 200);
        return `rgb(${r}, ${g}, ${b})`;
    };

    const datasets = (profileData.times_days || []).map((t, i) => ({
        label: `t = ${t.toFixed(2)} d`,
        data: xData.map((x, j) => ({
            x: x,
            y: profileData.Cx_mg_kg[i]?.[j] || 0,
        })),
        borderColor: colorScale(i),
        backgroundColor: 'transparent',
        pointRadius: 0,
        tension: 0.3,
        borderWidth: i === nProfiles - 1 ? 3 : 1.5,
    }));

    // Grid configuration
    const gridConfig = {
        color: 'rgba(0, 0, 0, 0.1)',
        drawBorder: true,
    };

    // X-axis configuration
    const xScaleConfig = {
        type: 'linear',
        title: { display: true, text: xLabel },
        grid: gridConfig,
    };
    if (!hasRealX) {
        xScaleConfig.min = 0;
        xScaleConfig.max = 1;
    }

    state.profileChart = new Chart(ctx, {
        type: 'line',
        data: { datasets },
        options: {
            responsive: true,
            maintainAspectRatio: false,
            parsing: false,
            interaction: {
                mode: 'nearest',
                intersect: false,
            },
            scales: {
                x: xScaleConfig,
                y: {
                    title: { display: true, text: 'Concentration (mg/kg)' },
                    beginAtZero: true,
                    grid: gridConfig,
                },
            },
            plugins: {
                legend: {
                    display: true,
                    position: 'right',
                    labels: {
                        boxWidth: 10,
                        font: { size: 10 },
                    },
                },
                tooltip: {
                    callbacks: {
                        label: (ctx) => {
                            const label = ctx.dataset.label || '';
                            const yVal = ctx.parsed.y;
                            const xVal = ctx.parsed.x;
                            return `${label}: C(${xVal.toFixed(1)}) = ${yVal.toFixed(2)} mg/kg`;
                        },
                    },
                },
                zoom: {
                    pan: {
                        enabled: true,
                        mode: 'xy',
                    },
                    zoom: {
                        wheel: {
                            enabled: true,
                        },
                        pinch: {
                            enabled: true,
                        },
                        mode: 'xy',
                    },
                },
            },
        },
    });

    // Grid toggle
    document.getElementById('profile-show-grid')?.addEventListener('change', (e) => {
        if (state.profileChart) {
            const show = e.target.checked;
            state.profileChart.options.scales.x.grid.display = show;
            state.profileChart.options.scales.y.grid.display = show;
            state.profileChart.update();
        }
    });
}

// ========== JOBS TAB ==========
async function loadJobsList() {
    const container = document.getElementById('jobs-list');
    const statusFilter = document.getElementById('jobs-filter-status').value;

    try {
        const params = new URLSearchParams();
        if (statusFilter) params.append('status', statusFilter);

        const data = await apiGet(`/jobs/list?${params.toString()}`);

        if (data.jobs.length === 0) {
            container.innerHTML = '<div class="text-gray-500 text-sm p-4 text-center">No jobs found</div>';
            return;
        }

        container.innerHTML = data.jobs.map(job => `
            <div class="border rounded-lg p-4 hover:bg-gray-50 job-item" data-job="${job.job_id}">
                <div class="flex items-center justify-between">
                    <div>
                        <div class="font-medium">${job.name}</div>
                        <div class="text-sm text-gray-500">
                            ${job.layer_count} layers, ${job.step_count} steps |
                            ${new Date(job.created_at).toLocaleString()}
                        </div>
                    </div>
                    <div class="flex items-center space-x-3">
                        <span class="px-2 py-1 rounded text-xs ${getStatusClass(job.status)}">${job.status}</span>
                        <button class="text-blue-600 hover:text-blue-800 load-job" data-job="${job.job_id}">Load</button>
                        <button class="text-red-500 hover:text-red-700 delete-job" data-job="${job.job_id}">Delete</button>
                    </div>
                </div>
            </div>
        `).join('');

        // Add listeners
        container.querySelectorAll('.load-job').forEach(btn => {
            btn.addEventListener('click', () => loadJob(btn.dataset.job));
        });

        container.querySelectorAll('.delete-job').forEach(btn => {
            btn.addEventListener('click', () => deleteJob(btn.dataset.job));
        });

        // Load stats
        loadJobStats();
    } catch (e) {
        console.error('Failed to load jobs:', e);
        container.innerHTML = '<div class="text-red-500 text-sm p-4 text-center">Failed to load jobs</div>';
    }
}

function getStatusClass(status) {
    const classes = {
        completed: 'bg-green-100 text-green-800',
        failed: 'bg-red-100 text-red-800',
        running: 'bg-blue-100 text-blue-800',
        draft: 'bg-gray-100 text-gray-800',
        pending: 'bg-yellow-100 text-yellow-800',
    };
    return classes[status] || 'bg-gray-100 text-gray-800';
}

async function loadJob(jobId) {
    try {
        const data = await apiGet(`/jobs/load/${jobId}/config`);
        if (data.config) {
            // Apply config to state
            if (data.config.layers) {
                state.layers = data.config.layers;
                setLayerCount(state.layers.length);
            }
            if (data.config.steps) state.steps = data.config.steps;
            if (data.config.food) state.food = data.config.food;
            if (data.config.geometry) state.geometry = data.config.geometry;

            document.getElementById('sim-name').value = data.name || '';
            alert(`Loaded job: ${data.name}`);

            // Switch to simulate tab
            document.querySelector('.tab-btn[data-tab="simulate"]').click();
        }
    } catch (e) {
        console.error('Failed to load job:', e);
        alert('Failed to load job');
    }
}

async function deleteJob(jobId) {
    if (!confirm('Delete this job?')) return;

    try {
        await apiDelete(`/jobs/delete/${jobId}`);
        loadJobsList();
    } catch (e) {
        console.error('Failed to delete job:', e);
    }
}

async function loadJobStats() {
    try {
        const data = await apiGet('/jobs/stats');
        const container = document.getElementById('jobs-stats');
        container.innerHTML = `
            <div class="grid grid-cols-4 gap-4 text-center">
                <div>
                    <div class="text-2xl font-bold">${data.total}</div>
                    <div class="text-xs text-gray-500">Total Jobs</div>
                </div>
                <div>
                    <div class="text-2xl font-bold text-green-600">${data.by_status?.completed || 0}</div>
                    <div class="text-xs text-gray-500">Completed</div>
                </div>
                <div>
                    <div class="text-2xl font-bold text-red-600">${data.by_status?.failed || 0}</div>
                    <div class="text-xs text-gray-500">Failed</div>
                </div>
                <div>
                    <div class="text-2xl font-bold">${data.total_size_mb}</div>
                    <div class="text-xs text-gray-500">MB Storage</div>
                </div>
            </div>
        `;
    } catch (e) {
        console.error('Failed to load stats:', e);
    }
}

function initJobsTab() {
    document.getElementById('jobs-filter-status').addEventListener('change', loadJobsList);

    document.getElementById('clear-old-jobs').addEventListener('click', async () => {
        if (!confirm('Delete jobs older than 30 days?')) return;
        try {
            const result = await apiDelete('/jobs/clear?older_than_days=30');
            alert(`Deleted ${result.deleted_count} old jobs`);
            loadJobsList();
        } catch (e) {
            console.error('Failed to clear jobs:', e);
        }
    });
}

// ========== CONFIG TAB ==========
async function loadConfig() {
    try {
        const data = await apiGet('/config/');
        const config = data.config;

        // Helper to safely set element value
        const setVal = (id, value) => {
            const el = document.getElementById(id);
            if (el) el.value = value;
        };

        // User
        setVal('config-user-name', config.user?.name || '');
        setVal('config-user-org', config.user?.organization || '');
        setVal('config-user-email', config.user?.email || '');
        setVal('config-user-role', config.user?.role || 'analyst');

        // Solver (Patankar settings)
        setVal('config-solver-nmesh', config.solver?.nmesh || 600);
        setVal('config-solver-nmeshmin', config.solver?.nmeshmin || 20);
        setVal('config-solver-ntimes', config.solver?.ntimes || 1000);
        setVal('config-solver-timescale', config.solver?.timescale || 'sqrt');
        // Advanced ODE settings
        setVal('config-solver-reltol', config.solver?.RelTol || '1e-6');
        setVal('config-solver-abstol', config.solver?.AbsTol || '1e-9');
        setVal('config-solver-maxstep', config.solver?.max_step || 'auto');

        // Store solver settings in state
        state.solverSettings = config.solver || {
            nmesh: 600, nmeshmin: 20, ntimes: 1000, timescale: 'sqrt'
        };

        // Estimators
        setVal('config-D-model', config.estimators?.D_model || 'piringer');
        setVal('config-k-model', config.estimators?.k_model || 'fhp');

        // Display
        setVal('config-unit-thickness', config.display?.units_thickness || 'um');
        setVal('config-unit-time', config.display?.units_time || 'days');
        setVal('config-decimals', config.display?.decimal_places || 3);

    } catch (e) {
        console.error('Failed to load config:', e);
    }
}

function initConfigTab() {
    // Helper to get element value safely
    const getVal = (id, defaultVal = '') => {
        const el = document.getElementById(id);
        return el ? el.value : defaultVal;
    };

    const saveBtn = document.getElementById('config-save');
    if (saveBtn) {
        saveBtn.addEventListener('click', async () => {
            const config = {
                user: {
                    name: getVal('config-user-name'),
                    organization: getVal('config-user-org'),
                    email: getVal('config-user-email'),
                    role: getVal('config-user-role', 'analyst'),
                },
                solver: {
                    // Patankar solver settings
                    nmesh: parseInt(getVal('config-solver-nmesh', '600')) || 600,
                    nmeshmin: parseInt(getVal('config-solver-nmeshmin', '20')) || 20,
                    ntimes: parseInt(getVal('config-solver-ntimes', '1000')) || 1000,
                    timescale: getVal('config-solver-timescale', 'sqrt'),
                    // Advanced ODE settings
                    RelTol: parseFloat(getVal('config-solver-reltol', '1e-6')) || 1e-6,
                    AbsTol: parseFloat(getVal('config-solver-abstol', '1e-9')) || 1e-9,
                    max_step: getVal('config-solver-maxstep', 'auto'),
                },
                estimators: {
                    D_model: getVal('config-D-model', 'piringer'),
                    k_model: getVal('config-k-model', 'fhp'),
                },
                display: {
                    units_thickness: getVal('config-unit-thickness', 'um'),
                    units_time: getVal('config-unit-time', 'days'),
                    decimal_places: parseInt(getVal('config-decimals', '3')) || 3,
                },
            };

            // Store solver settings in state for use in simulation
            state.solverSettings = config.solver;

            try {
                await apiPut('/config/', config);
                showToast('Configuration saved', 'success');
            } catch (e) {
                console.error('Failed to save config:', e);
                showToast('Failed to save configuration', 'error');
            }
        });
    }

    const resetBtn = document.getElementById('config-reset');
    if (resetBtn) {
        resetBtn.addEventListener('click', async () => {
            if (!confirm('Reset all settings to defaults?')) return;
            try {
                await apiPost('/config/reset');
                loadConfig();
                showToast('Configuration reset to defaults', 'success');
            } catch (e) {
                console.error('Failed to reset config:', e);
            }
        });
    }
}

// ========== EXPORT FUNCTIONALITY ==========
async function initExports() {
    // Check available export formats
    try {
        const data = await apiGet('/export/formats');
        const formats = data.formats || {};

        // Update button states based on availability
        const exportBtns = {
            'export-csv-btn': formats.csv,
            'export-xlsx-btn': formats.xlsx,
            'export-json-btn': formats.json,
            'export-pdf-btn': formats.pdf,
            'export-png-btn': formats.png,
            'export-svg-btn': formats.svg,
        };

        Object.entries(exportBtns).forEach(([btnId, formatInfo]) => {
            const btn = document.getElementById(btnId);
            if (btn && formatInfo && !formatInfo.available) {
                btn.disabled = true;
                btn.classList.add('opacity-50', 'cursor-not-allowed');
                btn.title = `${formatInfo.description} - Not available (missing dependency)`;
            }
        });
    } catch (e) {
        console.warn('Could not check export formats:', e);
    }

    // Export button handlers
    document.getElementById('export-csv-btn')?.addEventListener('click', () => exportResults('csv'));
    document.getElementById('export-xlsx-btn')?.addEventListener('click', () => exportResults('xlsx'));
    document.getElementById('export-json-btn')?.addEventListener('click', () => exportResults('json'));
    document.getElementById('export-pdf-btn')?.addEventListener('click', () => exportResults('pdf'));
    document.getElementById('export-png-btn')?.addEventListener('click', () => exportResults('png'));
    document.getElementById('export-svg-btn')?.addEventListener('click', () => exportResults('svg'));

    // Report modal button
    document.getElementById('show-report-btn')?.addEventListener('click', showReportModal);
}

// ========== PRINT-FRIENDLY REPORT ==========
function showReportModal() {
    if (!currentResults) {
        alert('No simulation results to display. Run a simulation first.');
        return;
    }

    const modal = document.getElementById('report-modal');
    const content = document.getElementById('report-content');
    if (!modal || !content) return;

    // Build report HTML
    const results = currentResults;
    const substances = results.substances || [];
    const tcontactDays = results.tcontact_days || 0;

    let html = `
        <div class="text-center mb-8 border-b pb-6">
            <h1 class="text-3xl font-bold text-green-600 mb-2">🍏⏩🍎 SFPPy Studio</h1>
            <p class="text-xl text-gray-600">Food Contact Migration Analysis Report</p>
            <p class="text-sm text-gray-400 mt-2">Generated: ${new Date().toLocaleString()}</p>
        </div>

        <section class="mb-8">
            <h2 class="text-xl font-bold text-gray-800 border-b pb-2 mb-4">📋 Simulation Summary</h2>
            <div class="grid grid-cols-2 gap-4 text-sm">
                <div><strong>Job ID:</strong> ${state.currentJobId || 'N/A'}</div>
                <div><strong>Contact Time:</strong> ${tcontactDays.toFixed(1)} days</div>
                <div><strong>Substances:</strong> ${substances.length}</div>
                <div><strong>Status:</strong> ${results.compliant ? '<span class="text-green-600 font-bold">✓ COMPLIANT</span>' : '<span class="text-red-600 font-bold">✗ NON-COMPLIANT</span>'}</div>
            </div>
        </section>

        <section class="mb-8">
            <h2 class="text-xl font-bold text-gray-800 border-b pb-2 mb-4">🧪 Substances & Results</h2>
            <table class="w-full text-sm border-collapse">
                <thead>
                    <tr class="bg-gray-100">
                        <th class="border px-3 py-2 text-left">Substance</th>
                        <th class="border px-3 py-2 text-right">CF @ t<sub>c</sub> (mg/kg)</th>
                        <th class="border px-3 py-2 text-right">SML (mg/kg)</th>
                        <th class="border px-3 py-2 text-right">Margin (%)</th>
                        <th class="border px-3 py-2 text-center">Status</th>
                    </tr>
                </thead>
                <tbody>
                    ${substances.map(sub => `
                        <tr>
                            <td class="border px-3 py-2">
                                <span style="display:inline-block; width:12px; height:12px; background:${sub.color}; border-radius:2px; margin-right:8px;"></span>
                                ${sub.name}
                            </td>
                            <td class="border px-3 py-2 text-right font-mono">${sub.CF_at_tcontact?.toFixed(4) || '--'}</td>
                            <td class="border px-3 py-2 text-right">${sub.SML_mg_kg?.toFixed(1) || 'N/A'}</td>
                            <td class="border px-3 py-2 text-right">${sub.margin_percent?.toFixed(1) || '--'}</td>
                            <td class="border px-3 py-2 text-center ${sub.compliant !== false ? 'text-green-600' : 'text-red-600'} font-bold">
                                ${sub.compliant !== false ? '✓' : '✗'}
                            </td>
                        </tr>
                    `).join('')}
                </tbody>
            </table>
        </section>

        <section class="mb-8">
            <h2 class="text-xl font-bold text-gray-800 border-b pb-2 mb-4">🧱 Layer Configuration</h2>
            <table class="w-full text-sm border-collapse">
                <thead>
                    <tr class="bg-gray-100">
                        <th class="border px-3 py-2 text-left">#</th>
                        <th class="border px-3 py-2 text-left">Material</th>
                        <th class="border px-3 py-2 text-right">Thickness</th>
                    </tr>
                </thead>
                <tbody>
                    ${state.layers.map((layer, i) => `
                        <tr>
                            <td class="border px-3 py-2">${i + 1}</td>
                            <td class="border px-3 py-2">${layer.material || 'Unknown'}</td>
                            <td class="border px-3 py-2 text-right">${layer.thickness || 0} ${layer.thickness_unit || 'µm'}</td>
                        </tr>
                    `).join('')}
                </tbody>
            </table>
        </section>

        <section class="mb-8">
            <h2 class="text-xl font-bold text-gray-800 border-b pb-2 mb-4">🌡️ Contact Conditions</h2>
            <table class="w-full text-sm border-collapse">
                <thead>
                    <tr class="bg-gray-100">
                        <th class="border px-3 py-2 text-left">Step</th>
                        <th class="border px-3 py-2 text-right">Temperature</th>
                        <th class="border px-3 py-2 text-right">Duration</th>
                        <th class="border px-3 py-2 text-center">Food Contact</th>
                    </tr>
                </thead>
                <tbody>
                    ${state.steps.map((step, i) => `
                        <tr>
                            <td class="border px-3 py-2">${i + 1}</td>
                            <td class="border px-3 py-2 text-right">${step.temperature || 25}°C</td>
                            <td class="border px-3 py-2 text-right">${step.duration || 10} ${step.duration_unit || 'days'}</td>
                            <td class="border px-3 py-2 text-center">${step.with_food !== false ? '✓' : '–'}</td>
                        </tr>
                    `).join('')}
                </tbody>
            </table>
        </section>

        <section class="mb-8">
            <h2 class="text-xl font-bold text-gray-800 border-b pb-2 mb-4">Food & Geometry</h2>
            <div class="grid grid-cols-2 gap-4 text-sm">
                <div><strong>Food Category:</strong> ${state.food.category || 'Not specified'}</div>
                <div><strong>Simulant:</strong> ${state.food.simulant || 'Not specified'}</div>
                <div><strong>Shape:</strong> ${state.geometry.shape || 'cylinder'}</div>
                <div><strong>Volume:</strong> ${state.geometry.volume_cm3?.toFixed(2) || '--'} cm³</div>
                <div><strong>Surface:</strong> ${state.geometry.surface_cm2?.toFixed(2) || '--'} cm²</div>
                <div><strong>V/S Ratio:</strong> ${state.geometry.vs_ratio_cm?.toFixed(3) || '--'} cm</div>
            </div>
        </section>

        <div class="text-center text-xs text-gray-400 mt-12 pt-4 border-t">
            <p>Generated by SFPPy Studio - Scientific Framework for Food Packaging Migration Prediction</p>
            <p class="mt-1">© INRAE + Generative Simulation | ${new Date().getFullYear()}</p>
        </div>
    `;

    content.innerHTML = html;
    modal.classList.remove('hidden');
    document.body.style.overflow = 'hidden';
}

function closeReportModal() {
    const modal = document.getElementById('report-modal');
    if (modal) {
        modal.classList.add('hidden');
        document.body.style.overflow = '';
    }
}

async function exportResults(format) {
    if (!state.currentJobId) {
        alert('No simulation results to export. Run a simulation first.');
        return;
    }

    const btn = document.getElementById(`export-${format}-btn`);
    const originalText = btn?.textContent;

    try {
        // Show loading state
        if (btn) {
            btn.disabled = true;
            btn.textContent = '...';
        }

        // Create download URL
        const url = `/api/export/${state.currentJobId}/${format}`;

        // Fetch and download
        const response = await fetch(url);

        if (!response.ok) {
            const errorData = await response.json().catch(() => ({}));
            throw new Error(errorData.detail || `Export failed: ${response.status}`);
        }

        // Get filename from Content-Disposition header
        const disposition = response.headers.get('Content-Disposition');
        let filename = `SFPPy_export.${format}`;
        if (disposition) {
            const match = disposition.match(/filename="?([^"]+)"?/);
            if (match) filename = match[1];
        }

        // Create blob and download
        const blob = await response.blob();
        const downloadUrl = URL.createObjectURL(blob);

        const a = document.createElement('a');
        a.href = downloadUrl;
        a.download = filename;
        document.body.appendChild(a);
        a.click();
        document.body.removeChild(a);

        URL.revokeObjectURL(downloadUrl);

        // Show success
        showExportToast(`Exported as ${format.toUpperCase()}`);

    } catch (e) {
        console.error('Export failed:', e);
        alert(`Export failed: ${e.message}`);
    } finally {
        // Restore button state
        if (btn) {
            btn.disabled = false;
            btn.textContent = originalText;
        }
    }
}

function showExportToast(message) {
    showToast(message, 'success');
}

// General toast notification function
function showToast(message, type = 'info', duration = 3000) {
    const colors = {
        success: 'bg-sfppy-green',
        error: 'bg-red-500',
        warning: 'bg-orange-500',
        info: 'bg-blue-500',
    };
    const icons = {
        success: '✓',
        error: '✗',
        warning: '⚠',
        info: 'ℹ',
    };

    const toast = document.createElement('div');
    toast.className = `fixed bottom-4 right-4 ${colors[type] || colors.info} text-white px-4 py-3 rounded-lg shadow-lg z-50 flex items-center gap-2`;
    toast.style.animation = 'slideIn 0.3s ease';
    toast.innerHTML = `<span class="text-lg">${icons[type] || icons.info}</span><span>${message}</span>`;
    document.body.appendChild(toast);

    setTimeout(() => {
        toast.style.opacity = '0';
        toast.style.transition = 'opacity 0.3s';
        setTimeout(() => toast.remove(), 300);
    }, duration);
}

// ========== FITTING STATE ==========
const fittingState = {
    dataSource: 'synthetic',  // 'synthetic', 'experimental', 'upload'
    currentJobId: null,
    data: null,  // { time_days: [], CF: [], CF_clean: [] (for synthetic) }
    trueValues: null,  // { D: ..., k: ... } for synthetic data comparison
    dataChart: null,
    resultChart: null,
    presets: [],
};

// ========== FITTING TAB ==========
async function initFittingTab() {
    // Load presets
    await loadFittingPresets();

    // Load existing fitting jobs
    await loadFittingJobs();

    // Data source selection buttons
    document.querySelectorAll('.fitting-source-btn').forEach(btn => {
        btn.addEventListener('click', () => {
            selectFittingSource(btn.dataset.source);
        });
    });

    // Generate synthetic data button
    document.getElementById('generate-synthetic-data')?.addEventListener('click', generateSyntheticData);

    // Load experimental data button
    document.getElementById('load-experimental-data')?.addEventListener('click', loadExperimentalData);

    // Upload CSV button
    document.getElementById('upload-csv-data')?.addEventListener('click', uploadCSVData);

    // Run fitting button
    document.getElementById('run-fitting')?.addEventListener('click', runFitting);
}

async function loadFittingPresets() {
    try {
        const result = await apiGet('/fitting/presets');
        if (result.success && result.presets) {
            fittingState.presets = result.presets;
            const presetsContainer = document.getElementById('fitting-presets');
            if (presetsContainer) {
                presetsContainer.innerHTML = result.presets.map(preset => `
                    <button class="preset-btn px-2 py-1 text-xs border border-gray-300 dark:border-gray-600
                                   rounded hover:bg-gray-100 dark:hover:bg-gray-700
                                   text-gray-700 dark:text-gray-300"
                            data-preset="${preset.name}"
                            onclick="applyFittingPreset('${preset.name}')">
                        ${preset.name}
                    </button>
                `).join('');
            }
        }
    } catch (e) {
        console.warn('Failed to load fitting presets:', e);
    }
}

function applyFittingPreset(presetName) {
    const preset = fittingState.presets.find(p => p.name === presetName);
    if (!preset) return;

    // Apply values to synthetic config inputs
    if (preset.layer?.thickness_um !== undefined) {
        document.getElementById('fit-thickness').value = preset.layer.thickness_um;
    }
    if (preset.layer?.C0 !== undefined) {
        document.getElementById('fit-c0').value = preset.layer.C0;
    }
    if (preset.layer?.D_true !== undefined) {
        document.getElementById('fit-d-true').value = preset.layer.D_true.toExponential(2);
    }
    if (preset.layer?.k_true !== undefined) {
        document.getElementById('fit-k-true').value = preset.layer.k_true;
    }
    if (preset.contact_time_days !== undefined) {
        document.getElementById('fit-contact-time').value = preset.contact_time_days;
    }
    if (preset.food?.volume_L !== undefined) {
        document.getElementById('fit-volume').value = preset.food.volume_L;
    }
    if (preset.food?.surface_dm2 !== undefined) {
        document.getElementById('fit-surface-area').value = preset.food.surface_dm2;
    }
    if (preset.noise_std !== undefined) {
        document.getElementById('fit-noise').value = preset.noise_std;
    }
    if (preset.n_points !== undefined) {
        document.getElementById('fit-n-points').value = preset.n_points;
    }

    // Visual feedback
    document.querySelectorAll('.preset-btn').forEach(btn => {
        if (btn.dataset.preset === presetName) {
            btn.classList.add('bg-sfppy-green', 'text-white', 'border-sfppy-green');
        } else {
            btn.classList.remove('bg-sfppy-green', 'text-white', 'border-sfppy-green');
        }
    });
}

function selectFittingSource(source) {
    fittingState.dataSource = source;

    // Update button styles
    document.querySelectorAll('.fitting-source-btn').forEach(btn => {
        if (btn.dataset.source === source) {
            btn.classList.add('border-sfppy-green', 'bg-green-50', 'dark:bg-green-900/30', 'text-sfppy-green', 'font-medium');
            btn.classList.remove('border-gray-300', 'dark:border-gray-600', 'text-gray-700', 'dark:text-gray-300');
        } else {
            btn.classList.remove('border-sfppy-green', 'bg-green-50', 'dark:bg-green-900/30', 'text-sfppy-green', 'font-medium');
            btn.classList.add('border-gray-300', 'dark:border-gray-600', 'text-gray-700', 'dark:text-gray-300');
        }
    });

    // Show/hide config sections
    document.getElementById('fitting-synthetic-config').classList.toggle('hidden', source !== 'synthetic');
    document.getElementById('fitting-experimental-config').classList.toggle('hidden', source !== 'experimental');
    document.getElementById('fitting-upload-config').classList.toggle('hidden', source !== 'upload');

    // Hide data preview and results when switching sources
    document.getElementById('fitting-data-preview').classList.add('hidden');
    document.getElementById('fitting-options').classList.add('hidden');
    document.getElementById('fitting-results').classList.add('hidden');
}

async function generateSyntheticData() {
    const btn = document.getElementById('generate-synthetic-data');
    const originalText = btn.textContent;

    try {
        btn.disabled = true;
        btn.textContent = '⏳ Generating...';

        // Collect config from inputs
        const config = {
            layer: {
                thickness_um: parseFloat(document.getElementById('fit-thickness').value) || 100,
                C0: parseFloat(document.getElementById('fit-c0').value) || 1000,
                D_true: parseFloat(document.getElementById('fit-d-true').value) || 1e-14,
                k_true: parseFloat(document.getElementById('fit-k-true').value) || 0.1,
            },
            food: {
                volume_L: parseFloat(document.getElementById('fit-volume').value) || 1,
                surface_dm2: parseFloat(document.getElementById('fit-surface-area').value) || 6,
            },
            contact_time_days: parseFloat(document.getElementById('fit-contact-time').value) || 10,
            noise_std: parseFloat(document.getElementById('fit-noise').value) || 0.01,
            n_points: parseInt(document.getElementById('fit-n-points').value) || 30,
        };

        // Add seed if provided
        const seed = document.getElementById('fit-seed').value;
        if (seed) {
            config.seed = parseInt(seed);
        }

        const result = await apiPost('/fitting/synthetic/generate', config);

        if (result.success) {
            fittingState.currentJobId = result.job_id;
            fittingState.data = {
                time_days: result.time_days,
                CF: result.CF_noisy,
                CF_clean: result.CF_clean,
            };
            fittingState.trueValues = {
                D: config.layer.D_true,
                k: config.layer.k_true,
            };

            // Show data preview
            showDataPreview(result);

            // Show fitting options
            document.getElementById('fitting-options').classList.remove('hidden');

        } else {
            throw new Error(result.detail || 'Failed to generate synthetic data');
        }

    } catch (e) {
        console.error('Error generating synthetic data:', e);
        alert(`Error: ${e.message}`);
    } finally {
        btn.disabled = false;
        btn.textContent = originalText;
    }
}

async function loadExperimentalData() {
    const btn = document.getElementById('load-experimental-data');
    const originalText = btn.textContent;

    try {
        btn.disabled = true;
        btn.textContent = '⏳ Loading...';

        // Parse textarea data
        const dataText = document.getElementById('fit-exp-data').value.trim();
        if (!dataText) {
            throw new Error('Please enter experimental data');
        }

        const thickness_um = parseFloat(document.getElementById('fit-exp-thickness').value) || 100;
        const C0 = parseFloat(document.getElementById('fit-exp-c0').value) || 1000;

        // Use the CSV endpoint which accepts plain CSV text
        const result = await apiPost('/fitting/experimental/csv', {
            data_csv: dataText,
            thickness_um: thickness_um,
            C0: C0,
        });

        if (result.success) {
            fittingState.currentJobId = result.job_id;
            fittingState.data = {
                time_days: result.time_days,
                CF: result.CF,
            };
            fittingState.trueValues = null;  // No true values for experimental data

            // Show data preview
            showDataPreview(result);

            // Show fitting options
            document.getElementById('fitting-options').classList.remove('hidden');

        } else {
            throw new Error(result.detail || 'Failed to load experimental data');
        }

    } catch (e) {
        console.error('Error loading experimental data:', e);
        alert(`Error: ${e.message}`);
    } finally {
        btn.disabled = false;
        btn.textContent = originalText;
    }
}

async function uploadCSVData() {
    const btn = document.getElementById('upload-csv-data');
    const originalText = btn.textContent;
    const fileInput = document.getElementById('fit-csv-file');

    try {
        btn.disabled = true;
        btn.textContent = '⏳ Uploading...';

        const file = fileInput.files[0];
        if (!file) {
            throw new Error('Please select a CSV file');
        }

        const thickness_um = parseFloat(document.getElementById('fit-upload-thickness').value) || 100;
        const C0 = parseFloat(document.getElementById('fit-upload-c0').value) || 1000;

        // Read file content
        const fileContent = await file.text();

        // Upload via API
        const formData = new FormData();
        formData.append('file', file);
        formData.append('thickness_um', thickness_um);
        formData.append('C0', C0);

        const response = await fetch('/api/fitting/experimental/upload', {
            method: 'POST',
            body: formData,
        });

        const result = await response.json();

        if (result.success) {
            fittingState.currentJobId = result.job_id;
            fittingState.data = {
                time_days: result.time_days,
                CF: result.CF,
            };
            fittingState.trueValues = null;

            // Show data preview
            showDataPreview(result);

            // Show fitting options
            document.getElementById('fitting-options').classList.remove('hidden');

        } else {
            throw new Error(result.detail || 'Failed to upload CSV');
        }

    } catch (e) {
        console.error('Error uploading CSV:', e);
        alert(`Error: ${e.message}`);
    } finally {
        btn.disabled = false;
        btn.textContent = originalText;
    }
}

function showDataPreview(result) {
    const previewDiv = document.getElementById('fitting-data-preview');
    const jobIdSpan = document.getElementById('fitting-job-id');
    const statsDiv = document.getElementById('fitting-data-stats');

    // If no result passed, build from fittingState
    if (!result) {
        result = {
            job_id: fittingState.currentJobId,
            time_days: fittingState.data?.time_days || [],
            CF: fittingState.data?.CF || [],
            CF_noisy: fittingState.data?.CF || [],
            CF_clean: fittingState.data?.CF_clean,
        };
    }

    // Show job ID
    jobIdSpan.textContent = `Job: ${result.job_id || fittingState.currentJobId}`;

    // Show stats
    const timeDays = result.time_days || [];
    const cfValues = result.CF_noisy || result.CF || [];
    const n = timeDays.length;

    if (n === 0 || cfValues.length === 0) {
        statsDiv.innerHTML = '<span class="text-gray-400">No data</span>';
        previewDiv.classList.remove('hidden');
        return;
    }

    const maxCF = Math.max(...cfValues);
    const minCF = Math.min(...cfValues);
    const maxTime = Math.max(...timeDays);

    let statsHtml = `<span class="mr-4">📊 ${n} points</span>`;
    statsHtml += `<span class="mr-4">⏱️ 0 - ${maxTime.toFixed(1)} days</span>`;
    statsHtml += `<span class="mr-4">📈 CF: ${minCF.toFixed(2)} - ${maxCF.toFixed(2)} mg/kg</span>`;

    if (result.CF_equilibrium !== undefined) {
        statsHtml += `<span>∞ CF_eq: ${result.CF_equilibrium.toFixed(2)} mg/kg</span>`;
    }

    statsDiv.innerHTML = statsHtml;

    // Render chart
    renderDataPreviewChart(result);

    // Show preview section
    previewDiv.classList.remove('hidden');
}

function renderDataPreviewChart(result) {
    const ctx = document.getElementById('fitting-data-chart')?.getContext('2d');
    if (!ctx) return;

    // Destroy existing chart
    if (fittingState.dataChart) {
        fittingState.dataChart.destroy();
    }

    const datasets = [];

    // Grid configuration
    const gridConfig = {
        color: 'rgba(0, 0, 0, 0.1)',
        drawBorder: true,
    };

    // Add noisy data points (if synthetic) or experimental data
    datasets.push({
        label: result.CF_noisy ? 'Data (with noise)' : 'Experimental Data',
        data: result.time_days.map((t, i) => ({
            x: t,
            y: (result.CF_noisy || result.CF)[i],
        })),
        borderColor: '#22C55E',
        backgroundColor: 'rgba(34, 197, 94, 0.5)',
        pointRadius: 5,
        pointHoverRadius: 7,
        pointStyle: 'circle',
        showLine: false,
    });

    // Add clean curve if synthetic
    if (result.CF_clean) {
        datasets.push({
            label: 'True Curve',
            data: result.time_days.map((t, i) => ({
                x: t,
                y: result.CF_clean[i],
            })),
            borderColor: '#3B82F6',
            backgroundColor: 'transparent',
            pointRadius: 0,
            showLine: true,
            borderWidth: 2,
        });
    }

    fittingState.dataChart = new Chart(ctx, {
        type: 'scatter',
        data: { datasets },
        options: {
            responsive: true,
            maintainAspectRatio: false,
            interaction: {
                mode: 'nearest',
                intersect: false,
            },
            plugins: {
                legend: {
                    display: true,
                    position: 'top',
                },
                tooltip: {
                    callbacks: {
                        label: (ctx) => {
                            const label = ctx.dataset.label || '';
                            return `${label}: t=${ctx.parsed.x.toFixed(2)}d, CF=${ctx.parsed.y.toFixed(2)} mg/kg`;
                        },
                    },
                },
                zoom: {
                    pan: {
                        enabled: true,
                        mode: 'xy',
                    },
                    zoom: {
                        wheel: {
                            enabled: true,
                        },
                        pinch: {
                            enabled: true,
                        },
                        mode: 'xy',
                    },
                },
            },
            scales: {
                x: {
                    title: {
                        display: true,
                        text: 'Time (days)',
                    },
                    type: 'linear',
                    beginAtZero: true,
                    grid: gridConfig,
                },
                y: {
                    title: {
                        display: true,
                        text: 'CF (mg/kg)',
                    },
                    beginAtZero: true,
                    grid: gridConfig,
                },
            },
        },
    });
}

async function runFitting() {
    const btn = document.getElementById('run-fitting');
    const originalText = btn.textContent;

    try {
        btn.disabled = true;
        btn.textContent = '⏳ Fitting...';

        if (!fittingState.currentJobId || !fittingState.data) {
            throw new Error('Please generate or load data first');
        }

        // Collect fitting options
        const fitD = document.getElementById('fit-D-enabled').checked;
        const fitK = document.getElementById('fit-k-enabled').checked;

        if (!fitD && !fitK) {
            throw new Error('Please enable at least one parameter to fit (D or k)');
        }

        const bounds = {
            D_min: parseFloat(document.getElementById('fit-d-min').value) || 1e-18,
            D_max: parseFloat(document.getElementById('fit-d-max').value) || 1e-10,
            k_min: parseFloat(document.getElementById('fit-k-min').value) || 0.001,
            k_max: parseFloat(document.getElementById('fit-k-max').value) || 1000,
        };

        const request = {
            job_id: fittingState.currentJobId,
            fit_D: fitD,
            fit_k: fitK,
            bounds: bounds,
        };

        // Add true values for comparison if available
        if (fittingState.trueValues) {
            request.D_true = fittingState.trueValues.D;
            request.k_true = fittingState.trueValues.k;
        }

        const result = await apiPost('/fitting/fit', request);

        if (result.success) {
            showFittingResults(result);
        } else {
            throw new Error(result.detail || 'Fitting failed');
        }

    } catch (e) {
        console.error('Error running fitting:', e);
        const errorMsg = e.message || (typeof e === 'object' ? JSON.stringify(e) : String(e));
        showToast(`Fitting error: ${errorMsg}`, 'error', 5000);
    } finally {
        btn.disabled = false;
        btn.textContent = originalText;
    }
}

function showFittingResults(result) {
    const resultsDiv = document.getElementById('fitting-results');

    // Show fitted parameters
    document.getElementById('fit-result-D').textContent = result.D_fitted?.toExponential(3) || 'N/A';
    document.getElementById('fit-result-k').textContent = result.k_fitted?.toFixed(4) || 'N/A';

    // Show fit quality
    document.getElementById('fit-result-R2').textContent = result.R2?.toFixed(4) || 'N/A';
    document.getElementById('fit-result-RMSE').textContent = result.RMSE?.toFixed(4) || 'N/A';

    // Show comparison with true values if available
    const comparisonDiv = document.getElementById('fit-comparison');
    if (result.comparison) {
        document.getElementById('fit-compare-D-true').textContent = result.comparison.D_true?.toExponential(3) || 'N/A';
        document.getElementById('fit-compare-k-true').textContent = result.comparison.k_true?.toFixed(4) || 'N/A';

        const dError = result.comparison.D_error_percent;
        const kError = result.comparison.k_error_percent;

        const dErrorSpan = document.getElementById('fit-compare-D-error');
        const kErrorSpan = document.getElementById('fit-compare-k-error');

        dErrorSpan.textContent = dError !== undefined ? `${dError > 0 ? '+' : ''}${dError.toFixed(2)}%` : 'N/A';
        kErrorSpan.textContent = kError !== undefined ? `${kError > 0 ? '+' : ''}${kError.toFixed(2)}%` : 'N/A';

        // Color code errors
        dErrorSpan.className = `font-mono ${Math.abs(dError) < 5 ? 'text-green-600' : Math.abs(dError) < 20 ? 'text-yellow-600' : 'text-red-600'}`;
        kErrorSpan.className = `font-mono ${Math.abs(kError) < 5 ? 'text-green-600' : Math.abs(kError) < 20 ? 'text-yellow-600' : 'text-red-600'}`;

        comparisonDiv.classList.remove('hidden');
    } else {
        comparisonDiv.classList.add('hidden');
    }

    // Render result chart
    renderFittingResultChart(result);

    // Show results section
    resultsDiv.classList.remove('hidden');

    // Update jobs list
    loadFittingJobs();
}

function renderFittingResultChart(result) {
    const ctx = document.getElementById('fitting-result-chart')?.getContext('2d');
    if (!ctx) return;

    // Destroy existing chart
    if (fittingState.resultChart) {
        fittingState.resultChart.destroy();
    }

    const datasets = [];

    // Grid configuration
    const gridConfig = {
        color: 'rgba(0, 0, 0, 0.1)',
        drawBorder: true,
    };

    // Experimental/noisy data
    datasets.push({
        label: 'Data',
        data: fittingState.data.time_days.map((t, i) => ({
            x: t,
            y: fittingState.data.CF[i],
        })),
        borderColor: '#22C55E',
        backgroundColor: 'rgba(34, 197, 94, 0.5)',
        pointRadius: 5,
        pointHoverRadius: 7,
        pointStyle: 'circle',
        showLine: false,
    });

    // True curve (if synthetic)
    if (fittingState.data.CF_clean) {
        datasets.push({
            label: 'True Curve',
            data: fittingState.data.time_days.map((t, i) => ({
                x: t,
                y: fittingState.data.CF_clean[i],
            })),
            borderColor: '#3B82F6',
            backgroundColor: 'transparent',
            pointRadius: 0,
            showLine: true,
            borderWidth: 2,
            borderDash: [5, 5],
        });
    }

    // Fitted curve - handle both nested (fitted_curve) and flat (CF_fitted) structures
    const fittedCF = result.fitted_curve?.CF_fitted || result.CF_fitted;
    const fittedTime = result.fitted_curve?.time_days || fittingState.data.time_days;

    if (fittedCF && fittedCF.length > 0) {
        datasets.push({
            label: 'Fitted Curve',
            data: fittedTime.map((t, i) => ({
                x: t,
                y: fittedCF[i],
            })),
            borderColor: '#F97316',
            backgroundColor: 'transparent',
            pointRadius: 0,
            showLine: true,
            borderWidth: 2.5,
        });
    }

    // Build annotations for fit quality
    const annotations = {};
    if (result.R2 !== undefined) {
        annotations.r2Label = {
            type: 'label',
            xValue: fittingState.data.time_days[Math.floor(fittingState.data.time_days.length * 0.8)],
            yValue: Math.max(...fittingState.data.CF) * 0.9,
            content: `R² = ${result.R2.toFixed(4)}`,
            backgroundColor: 'rgba(249, 115, 22, 0.9)',
            borderRadius: 4,
            font: { size: 11, weight: 'bold' },
            color: 'white',
            padding: 4,
        };
    }

    fittingState.resultChart = new Chart(ctx, {
        type: 'scatter',
        data: { datasets },
        options: {
            responsive: true,
            maintainAspectRatio: false,
            interaction: {
                mode: 'nearest',
                intersect: false,
            },
            plugins: {
                legend: {
                    display: true,
                    position: 'top',
                },
                annotation: {
                    annotations: annotations,
                },
                tooltip: {
                    callbacks: {
                        label: (ctx) => {
                            const label = ctx.dataset.label || '';
                            return `${label}: t=${ctx.parsed.x.toFixed(2)}d, CF=${ctx.parsed.y.toFixed(2)} mg/kg`;
                        },
                    },
                },
                zoom: {
                    pan: {
                        enabled: true,
                        mode: 'xy',
                    },
                    zoom: {
                        wheel: {
                            enabled: true,
                        },
                        pinch: {
                            enabled: true,
                        },
                        mode: 'xy',
                    },
                },
            },
            scales: {
                x: {
                    title: {
                        display: true,
                        text: 'Time (days)',
                    },
                    type: 'linear',
                    beginAtZero: true,
                    grid: gridConfig,
                },
                y: {
                    title: {
                        display: true,
                        text: 'CF (mg/kg)',
                    },
                    beginAtZero: true,
                    grid: gridConfig,
                },
            },
        },
    });
}

async function loadFittingJobs() {
    try {
        const result = await apiGet('/fitting/jobs');
        const listDiv = document.getElementById('fitting-jobs-list');

        if (!result.success || !result.jobs || result.jobs.length === 0) {
            listDiv.innerHTML = '<span class="text-gray-400">No fitting jobs yet.</span>';
            return;
        }

        listDiv.innerHTML = result.jobs.slice(0, 5).map(job => {
            const status = job.has_fit_result ? 'fitted' : 'pending';
            const statusColor = job.has_fit_result ? 'text-green-600' : 'text-yellow-600';
            const dataType = job.type || 'unknown';
            return `
            <div class="flex items-center justify-between py-1 border-b border-gray-100 dark:border-gray-700 last:border-0 cursor-pointer hover:bg-gray-50 dark:hover:bg-gray-700"
                 onclick="selectFittingJob('${job.job_id}')">
                <div>
                    <span class="font-mono text-xs">${job.job_id.substring(0, 8)}...</span>
                    <span class="ml-2 text-xs ${statusColor}">${status}</span>
                </div>
                <span class="text-xs text-gray-400">${dataType}</span>
            </div>`;
        }).join('');

    } catch (e) {
        console.warn('Failed to load fitting jobs:', e);
    }
}

async function selectFittingJob(jobId) {
    try {
        const result = await apiGet(`/fitting/job/${jobId}`);
        if (result.success) {
            fittingState.currentJobId = jobId;
            // API returns data at top level, not nested under 'job'
            fittingState.data = {
                time_days: result.data?.time_days || [],
                CF: result.data?.CF || result.data?.CF_noisy || result.data?.CF_exp || [],
                CF_clean: result.data?.CF_clean,
            };
            // True values from config for synthetic data
            if (result.type === 'synthetic' && result.config?.layer) {
                fittingState.trueValues = {
                    D_true: result.config.layer.D_true,
                    k_true: result.config.layer.k_true,
                };
            } else {
                fittingState.trueValues = null;
            }

            showToast(`Selected job ${jobId.substring(0, 8)}...`, 'success');

            // Show the data preview
            if (fittingState.data.time_days.length > 0) {
                showDataPreview();
                // Show fitting options so user can run a fit
                document.getElementById('fitting-options')?.classList.remove('hidden');
            }

            // Load fit results if available
            if (result.fit_result) {
                showFittingResults({
                    success: true,
                    ...result.fit_result,
                    comparison: result.type === 'synthetic' && result.config?.layer ? {
                        D_true: result.config.layer.D_true,
                        k_true: result.config.layer.k_true,
                        D_error_percent: result.fit_result.D_fitted && result.config.layer.D_true
                            ? 100 * Math.abs(result.fit_result.D_fitted - result.config.layer.D_true) / result.config.layer.D_true
                            : null,
                        k_error_percent: result.fit_result.k_fitted && result.config.layer.k_true
                            ? 100 * Math.abs(result.fit_result.k_fitted - result.config.layer.k_true) / result.config.layer.k_true
                            : null,
                    } : result.fit_result.comparison,
                });
            }
        } else {
            showToast('Failed to load job: ' + (result.error || 'Unknown error'), 'error');
        }
    } catch (e) {
        console.error('Failed to select fitting job:', e);
        showToast('Failed to load job', 'error');
    }
}

// ========== HEALTH & QUEUE MONITORING ==========
let healthCheckInterval = null;
let lastHealthStatus = null;

/**
 * Check server health and update indicator
 */
async function checkServerHealth() {
    const indicator = document.getElementById('health-indicator');
    const container = document.getElementById('server-health');
    if (!indicator || !container) return;

    try {
        const response = await fetch('/api/health', { timeout: 5000 });
        const data = await response.json();

        if (response.ok && data.status === 'healthy') {
            indicator.className = 'w-2 h-2 rounded-full bg-green-400';
            indicator.classList.remove('animate-pulse');
            container.title = `Server: Healthy | SFPPy: ${data.sfppy_available ? 'Available' : 'Mock mode'}`;
            lastHealthStatus = 'healthy';
        } else {
            indicator.className = 'w-2 h-2 rounded-full bg-yellow-400 animate-pulse';
            container.title = 'Server: Degraded';
            lastHealthStatus = 'degraded';
        }
    } catch (e) {
        indicator.className = 'w-2 h-2 rounded-full bg-red-500 animate-pulse';
        container.title = 'Server: Unreachable - Click 🔄 to retry';
        lastHealthStatus = 'offline';
    }

    // Also update queue count
    await updateQueueCount();
}

/**
 * Update pending job count in header
 */
async function updateQueueCount() {
    const countEl = document.getElementById('queue-count');
    if (!countEl) return;

    try {
        const data = await apiGet('/jobs/list');
        const jobs = data.jobs || [];
        const pendingJobs = jobs.filter(j => j.status === 'pending' || j.status === 'running');
        countEl.textContent = pendingJobs.length;
        countEl.title = `${pendingJobs.length} pending/running job(s)`;

        // Highlight if there are pending jobs
        if (pendingJobs.length > 0) {
            countEl.classList.add('text-yellow-300');
        } else {
            countEl.classList.remove('text-yellow-300');
        }
    } catch (e) {
        countEl.textContent = '?';
    }
}

/**
 * Initialize health monitoring
 */
function initHealthMonitoring() {
    // Check immediately
    checkServerHealth();

    // Check every 30 seconds
    healthCheckInterval = setInterval(checkServerHealth, 30000);

    // Restart/refresh button
    const restartBtn = document.getElementById('restart-server-btn');
    if (restartBtn) {
        restartBtn.addEventListener('click', async () => {
            restartBtn.classList.add('animate-spin');
            await checkServerHealth();
            setTimeout(() => restartBtn.classList.remove('animate-spin'), 500);

            // Show current status
            if (lastHealthStatus === 'healthy') {
                showToast('Server is healthy', 'success');
            } else if (lastHealthStatus === 'offline') {
                showToast('Server unreachable - please restart the server', 'error');
            } else {
                showToast('Server status: ' + lastHealthStatus, 'warning');
            }
        });
    }

    // Health indicator click - show detailed status
    const healthDiv = document.getElementById('server-health');
    if (healthDiv) {
        healthDiv.addEventListener('click', async () => {
            try {
                const response = await fetch('/api/info');
                const info = await response.json();
                alert(`SFPPy Studio Server Info:\n\n` +
                    `API Version: ${info.api_version || 'unknown'}\n` +
                    `SFPPy: ${info.sfppy_version || 'N/A'}\n` +
                    `Modules: ${info.modules_available?.join(', ') || 'checking...'}`);
            } catch (e) {
                alert('Could not fetch server info. Server may be offline.');
            }
        });
    }
}

// ========== INITIALIZATION ==========
document.addEventListener('DOMContentLoaded', async () => {
    initTabs();
    initHelp();
    initAssembly();
    initFoodTab();
    initSubstancesTab();
    initSimulateTab();
    initJobsTab();
    initConfigTab();
    initExports();
    initFittingTab();
    initHealthMonitoring();

    // Initialize layer visualization
    updateLayerVisualization();

    // Show Substances tab first (as it's the logical first step)
    const substancesTab = document.querySelector('.tab-btn[data-tab="substances"]');
    if (substancesTab) {
        substancesTab.click();
    }

    // Initialize backend session
    await initSession();
});
