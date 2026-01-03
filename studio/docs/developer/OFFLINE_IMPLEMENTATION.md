# SFPPy Studio - Offline Implementation Specification

**Version:** 1.0
**Date:** 2026-01-02
**Author:** Olivier Vitrac, PhD, HDR

---

## 1. Definition of "Offline"

**"Offline" means no internet access, NOT no backend.**

The FastAPI server runs on `localhost` (or LAN). The browser must function without external CDN/API calls.

```
┌─────────────────┐      ┌─────────────────┐
│     Browser     │ ←──→ │  FastAPI Server │  ← ALWAYS AVAILABLE
│  (JavaScript)   │      │   (localhost)   │
└─────────────────┘      └─────────────────┘
         ↓ BLOCKED
┌─────────────────┐
│   Internet      │  ← CDNs, PubChem live queries
│  (CDNs, APIs)   │
└─────────────────┘
```

---

## 2. Current Storage Architecture

### 2.1 Storage Locations (Current State)

| Location | Purpose | Offline Status |
|----------|---------|----------------|
| `studio/jobs/` | Simulation job results (JSON) | ✅ Local filesystem |
| `studio/examples/` | Example session files | ✅ Local filesystem |
| `studio/data/` | Static data files | ✅ Local filesystem |
| `patankar/cache.PubChem/` | Cached PubChem data | ✅ Pre-populated |
| In-memory (`SessionStore`) | Current session state | ⚠️ Lost on restart |
| `storage/` (project root) | LlamaIndex RAG index | ❌ NOT related to Studio |

### 2.2 Required Changes for Offline

**Rule 1: All Studio storage MUST be inside `studio/` directory**

```
studio/
├── jobs/              # Job results (existing)
├── examples/          # Example sessions (existing)
├── data/              # Static data (existing)
├── storage/           # NEW: Move session state here
│   ├── sessions/      # Persisted sessions
│   └── user_prefs/    # User preferences
└── static/
    └── vendor/        # NEW: Vendored CDN libraries
```

**Rule 2: Browser storage for user preferences**

Use `localStorage`/`IndexedDB` for:
- User preferences (dark mode, help bar dismissed, etc.)
- Draft session state (auto-save)
- Offline queue (pending API calls)

---

## 3. CDN Dependencies to Vendor

### 3.1 Current External Dependencies

```html
<!-- From index.html - MUST BE VENDORED -->
<script src="https://cdn.tailwindcss.com"></script>
<script src="https://cdn.jsdelivr.net/npm/chart.js"></script>
```

### 3.2 Vendor Location

```
studio/app/static/vendor/
├── tailwindcss/
│   └── tailwind.min.css    # Pre-compiled CSS (not CDN JIT)
├── chartjs/
│   └── chart.min.js        # v4.x
└── LICENSES.md             # License attribution
```

### 3.3 Version Pinning

Create `studio/app/static/vendor/manifest.json`:
```json
{
  "tailwindcss": {"version": "3.4.0", "license": "MIT"},
  "chart.js": {"version": "4.4.1", "license": "MIT"}
}
```

---

## 4. Service Worker Strategy

### 4.1 Cache Strategies by Route

| Route Pattern | Strategy | Rationale |
|---------------|----------|-----------|
| `/static/**` | **CacheFirst** | Static assets, rarely change |
| `/api/substances/**` | **CacheFirst** | PubChem data is static |
| `/api/polymers/**` | **CacheFirst** | Polymer database is static |
| `/api/simulants/**` | **CacheFirst** | Simulant data is static |
| `/api/jobs/**` | **NetworkFirst** | Jobs may update, fallback to cache |
| `/api/sessions/**` | **NetworkFirst** | Sessions are dynamic |
| `/api/simulation/run` | **NetworkOnly** | Must execute on server |
| `/api/fitting/**` | **NetworkOnly** | Must execute on server |

### 4.2 Service Worker Implementation

```javascript
// service-worker.js
const CACHE_NAME = 'sfppy-studio-v1';
const STATIC_CACHE = 'sfppy-static-v1';

const PRECACHE_URLS = [
  '/',
  '/static/js/app.js',
  '/static/vendor/tailwind.min.css',
  '/static/vendor/chart.min.js',
  '/api/substances/common',
  '/api/polymers/list',
  '/api/simulants/list'
];
```

---

## 5. Offline Behavior by Feature

### 5.1 Features That Work Offline

| Feature | Offline Behavior |
|---------|------------------|
| View saved jobs | ✅ From cache/IndexedDB |
| Load example sessions | ✅ Pre-cached |
| UI navigation | ✅ Fully cached |
| Dark mode toggle | ✅ localStorage |
| Help tips | ✅ Embedded in JS |
| View cached substance data | ✅ Pre-cached |
| Edit session (draft) | ✅ localStorage auto-save |

### 5.2 Features Requiring Network

| Feature | Offline Behavior |
|---------|------------------|
| Run simulation | 🔴 Queue request, show message |
| Search PubChem (new) | 🔴 Show "offline" message |
| Fitting computation | 🔴 Queue request, show message |
| Export PDF/CSV | 🟡 From cached results only |

### 5.3 Graceful Degradation UI

When offline:
```
┌─────────────────────────────────────────────────────────┐
│ 📴 Offline Mode                                    [x]  │
│ Some features are unavailable. Queued requests: 2       │
└─────────────────────────────────────────────────────────┘
```

---

## 6. PubChem Data Strategy

### 6.1 Pre-cached Substances

SFPPy already caches regulated substances in `patankar/cache.PubChem/`:
- All EU 10/2011 substances
- All US FDA 21 CFR substances
- All Chinese GB 9685 substances
- ~500+ substances with full properties

### 6.2 Offline Substance Search

```javascript
// Search cached substances first
const cachedResults = await searchCachedSubstances(query);
if (cachedResults.length > 0) {
  return cachedResults;
}

// If online, query PubChem
if (navigator.onLine) {
  return await queryPubChem(query);
}

// Offline and not cached
return { error: 'Substance not in offline cache. Connect to internet.' };
```

---

## 7. Session Persistence Strategy

### 7.1 Current State (In-Memory)

```python
# session.py - Current
class SessionStore:
    _sessions: Dict[str, SimulationSession] = {}  # Lost on restart
```

### 7.2 Proposed Hybrid Storage

```
┌─────────────────────────────────────────────────────────┐
│                    Session State                         │
├─────────────────────────────────────────────────────────┤
│  Browser (IndexedDB)     │  Server (studio/storage/)    │
│  ─────────────────────   │  ─────────────────────────   │
│  • Draft auto-save       │  • Completed sessions        │
│  • User preferences      │  • Job results               │
│  • Offline queue         │  • Exported files            │
└─────────────────────────────────────────────────────────┘
```

### 7.3 Auto-Save to Browser

```javascript
// Auto-save draft every 30 seconds
setInterval(() => {
  if (state.modified) {
    localStorage.setItem('sfppy_draft_session', JSON.stringify(state));
    state.modified = false;
  }
}, 30000);

// Restore on page load
const draft = localStorage.getItem('sfppy_draft_session');
if (draft) {
  const restored = JSON.parse(draft);
  if (confirm('Restore previous session?')) {
    Object.assign(state, restored);
  }
}
```

---

## 8. Implementation Checklist

### Phase 1: Vendor CDN Dependencies
- [ ] Download and vendor Tailwind CSS (pre-compiled)
- [ ] Download and vendor Chart.js
- [ ] Create `manifest.json` with versions
- [ ] Update `index.html` to use local files
- [ ] Add cache-busting query strings

### Phase 2: Service Worker
- [ ] Create `service-worker.js`
- [ ] Add SW registration in `index.html`
- [ ] Implement precache for static assets
- [ ] Implement runtime caching strategies
- [ ] Add cache versioning and cleanup

### Phase 3: Offline UI
- [ ] Add offline indicator in header
- [ ] Add request queue for offline operations
- [ ] Show appropriate messages for unavailable features
- [ ] Implement graceful degradation

### Phase 4: Session Persistence
- [ ] Move `storage/` inside `studio/` (if used)
- [ ] Implement browser auto-save (localStorage/IndexedDB)
- [ ] Add session restore on page load
- [ ] Implement offline queue with retry

### Phase 5: Testing
- [ ] Test with Chrome DevTools "Offline" mode
- [ ] Test first load online → offline reload
- [ ] Test cache invalidation on version bump
- [ ] Test offline queue retry when back online

---

## 9. Files to Create/Modify

### New Files
```
studio/app/static/vendor/tailwind.min.css
studio/app/static/vendor/chart.min.js
studio/app/static/vendor/manifest.json
studio/app/static/vendor/LICENSES.md
studio/app/static/js/service-worker.js
studio/app/static/js/offline-manager.js
```

### Modified Files
```
studio/app/templates/index.html  # Replace CDN links, add SW registration
studio/app/static/js/app.js      # Add offline detection, auto-save
studio/app/main.py               # Add caching headers middleware
```

---

## 10. Notes

### The `storage/` Directory at Project Root

The `storage/` directory at `/home/olivi/natacha/python/storage/` contains LlamaIndex RAG data (docstore, vector store, etc.). This is **NOT related to Studio** and should be:
- Either moved to its own project
- Or clearly documented as separate from Studio

**Studio should NOT use this directory.**

### Browser Storage Limits

- `localStorage`: ~5-10MB (varies by browser)
- `IndexedDB`: ~50MB+ (practically unlimited with user permission)

For large job results (some are 1.9MB), use IndexedDB, not localStorage.

---

*Document generated for SFPPy Studio offline implementation planning.*
