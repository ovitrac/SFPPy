/**
 * SFPPy Studio Service Worker
 * Enables offline operation with intelligent caching strategies.
 *
 * Cache Strategies:
 * - CacheFirst: Static assets, substance data (rarely changes)
 * - NetworkFirst: Jobs, sessions (dynamic, but fallback to cache)
 * - NetworkOnly: Simulation runs, fitting (must execute on server)
 */

const CACHE_VERSION = 'v0.3';
const STATIC_CACHE = `sfppy-static-${CACHE_VERSION}`;
const DATA_CACHE = `sfppy-data-${CACHE_VERSION}`;
const IMAGE_CACHE = `sfppy-images-${CACHE_VERSION}`;

// URLs to precache on install
const PRECACHE_URLS = [
    '/',
    '/static/js/app.js',
    '/static/css/style.css',
    '/static/vendor/chart.min.js',
    '/static/vendor/tailwind-offline.css',
    '/static/vendor/manifest.json',
    '/api/assembly/polymers',
    '/api/food/simulants',
    '/api/substances/common',
];

// URLs that should use CacheFirst strategy
const CACHE_FIRST_PATTERNS = [
    /^\/static\//,
    /^\/api\/substances\/structure\//,
    /^\/api\/substances\/common/,
    /^\/api\/polymers\//,
    /^\/api\/simulants\//,
];

// URLs that should use NetworkFirst strategy
const NETWORK_FIRST_PATTERNS = [
    /^\/api\/jobs\//,
    /^\/api\/sessions\//,
    /^\/api\/examples\//,
];

// URLs that must always go to network (no caching)
const NETWORK_ONLY_PATTERNS = [
    /^\/api\/simulation\/run/,
    /^\/api\/fitting\//,
    /^\/api\/substances\/search/,  // Live PubChem search
    /^\/api\/substances\/detail/,  // May need fresh ToxTree data
];

/**
 * Install event - precache essential assets
 */
self.addEventListener('install', (event) => {
    console.log('[SW] Installing service worker...');
    event.waitUntil(
        caches.open(STATIC_CACHE)
            .then((cache) => {
                console.log('[SW] Precaching app shell');
                return cache.addAll(PRECACHE_URLS.map(url => {
                    return new Request(url, { cache: 'reload' });
                })).catch(err => {
                    console.warn('[SW] Some precache URLs failed:', err);
                    // Continue anyway - some URLs might not exist yet
                });
            })
            .then(() => {
                console.log('[SW] Skip waiting...');
                return self.skipWaiting();
            })
    );
});

/**
 * Activate event - clean old caches
 */
self.addEventListener('activate', (event) => {
    console.log('[SW] Activating service worker...');
    event.waitUntil(
        caches.keys()
            .then((cacheNames) => {
                return Promise.all(
                    cacheNames
                        .filter((name) => {
                            // Delete old versioned caches
                            return name.startsWith('sfppy-') &&
                                   !name.endsWith(CACHE_VERSION);
                        })
                        .map((name) => {
                            console.log('[SW] Deleting old cache:', name);
                            return caches.delete(name);
                        })
                );
            })
            .then(() => {
                console.log('[SW] Claiming clients...');
                return self.clients.claim();
            })
    );
});

/**
 * Fetch event - handle requests with appropriate strategy
 */
self.addEventListener('fetch', (event) => {
    const url = new URL(event.request.url);
    const path = url.pathname;

    // Skip non-GET requests
    if (event.request.method !== 'GET') {
        return;
    }

    // Skip cross-origin requests (except PubChem images)
    if (url.origin !== self.location.origin) {
        // Handle PubChem image requests specially
        if (url.hostname === 'pubchem.ncbi.nlm.nih.gov' && path.includes('/image/')) {
            event.respondWith(handlePubChemImage(event.request, url));
            return;
        }
        return;
    }

    // Determine caching strategy based on URL pattern
    if (matchesPattern(path, NETWORK_ONLY_PATTERNS)) {
        // Network only - don't intercept
        return;
    }

    if (matchesPattern(path, CACHE_FIRST_PATTERNS)) {
        event.respondWith(cacheFirst(event.request));
        return;
    }

    if (matchesPattern(path, NETWORK_FIRST_PATTERNS)) {
        event.respondWith(networkFirst(event.request));
        return;
    }

    // Default: Network first with cache fallback
    event.respondWith(networkFirst(event.request));
});

/**
 * Check if path matches any pattern in the array
 */
function matchesPattern(path, patterns) {
    return patterns.some(pattern => pattern.test(path));
}

/**
 * Cache-first strategy: Try cache, fallback to network
 */
async function cacheFirst(request) {
    const cached = await caches.match(request);
    if (cached) {
        // Return cached response, but update cache in background
        updateCache(request);
        return cached;
    }

    try {
        const response = await fetch(request);
        if (response.ok) {
            const cache = await caches.open(STATIC_CACHE);
            cache.put(request, response.clone());
        }
        return response;
    } catch (error) {
        console.warn('[SW] Network failed for cache-first:', request.url);
        return new Response('Offline - resource not cached', {
            status: 503,
            statusText: 'Service Unavailable'
        });
    }
}

/**
 * Network-first strategy: Try network, fallback to cache
 */
async function networkFirst(request) {
    try {
        const response = await fetch(request);
        if (response.ok) {
            const cache = await caches.open(DATA_CACHE);
            cache.put(request, response.clone());
        }
        return response;
    } catch (error) {
        console.warn('[SW] Network failed, trying cache:', request.url);
        const cached = await caches.match(request);
        if (cached) {
            return cached;
        }

        // Return offline response for API calls
        if (request.url.includes('/api/')) {
            return new Response(JSON.stringify({
                success: false,
                offline: true,
                error: 'You are offline. This data is not available in cache.'
            }), {
                status: 503,
                headers: { 'Content-Type': 'application/json' }
            });
        }

        return new Response('Offline - page not cached', {
            status: 503,
            statusText: 'Service Unavailable'
        });
    }
}

/**
 * Handle PubChem image requests - PubChem primary, local SFPPy fallback
 * Strategy: Use fast PubChem CDN when online, fall back to local cache offline
 */
async function handlePubChemImage(request, url) {
    // Extract CID from PubChem image URL
    const cidMatch = url.search.match(/cid=(\d+)/);
    if (!cidMatch) {
        return fetch(request);
    }

    const cid = cidMatch[1];
    const localUrl = `/api/substances/structure/${cid}.png`;

    // Check browser cache first (fastest)
    const cached = await caches.match(request);
    if (cached) {
        console.log('[SW] Serving cached PubChem image for CID:', cid);
        return cached;
    }

    // Try remote PubChem (primary source when online)
    try {
        const response = await fetch(request);
        if (response.ok) {
            console.log('[SW] Fetched PubChem image for CID:', cid);
            // Cache for future offline use
            const cache = await caches.open(IMAGE_CACHE);
            cache.put(request, response.clone());
            return response;
        }
    } catch (error) {
        console.log('[SW] PubChem unavailable for CID:', cid, '- trying local fallback');
    }

    // Fallback to local SFPPy cache (offline mode)
    try {
        const localResponse = await fetch(localUrl);
        if (localResponse.ok) {
            console.log('[SW] Serving local SFPPy image for CID:', cid);
            return localResponse;
        }
    } catch (e) {
        // Local also not available
    }

    // Return 404 if no source available
    return new Response('', {
        status: 404,
        statusText: 'Image not available'
    });
}

/**
 * Update cache in background (stale-while-revalidate pattern)
 */
async function updateCache(request) {
    try {
        const response = await fetch(request);
        if (response.ok) {
            const cache = await caches.open(STATIC_CACHE);
            await cache.put(request, response);
        }
    } catch (error) {
        // Silent fail - we already served from cache
    }
}

/**
 * Message handler for cache management
 */
self.addEventListener('message', (event) => {
    if (event.data.type === 'SKIP_WAITING') {
        self.skipWaiting();
    }

    if (event.data.type === 'CACHE_SUBSTANCE') {
        // Pre-cache a substance's data and image
        const cid = event.data.cid;
        if (cid) {
            caches.open(IMAGE_CACHE).then(cache => {
                cache.add(`/api/substances/structure/${cid}.png`);
            });
        }
    }

    if (event.data.type === 'GET_CACHE_STATUS') {
        Promise.all([
            caches.open(STATIC_CACHE).then(c => c.keys()),
            caches.open(DATA_CACHE).then(c => c.keys()),
            caches.open(IMAGE_CACHE).then(c => c.keys()),
        ]).then(([staticKeys, dataKeys, imageKeys]) => {
            event.ports[0].postMessage({
                static: staticKeys.length,
                data: dataKeys.length,
                images: imageKeys.length,
                version: CACHE_VERSION
            });
        });
    }

    if (event.data.type === 'CLEAR_CACHE') {
        caches.keys().then(names => {
            Promise.all(names.map(name => caches.delete(name)));
        }).then(() => {
            event.ports[0].postMessage({ cleared: true });
        });
    }
});

console.log('[SW] Service worker loaded');
