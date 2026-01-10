/**
 * SFPPy Studio Offline Manager
 * Handles offline detection, request queuing, and graceful degradation.
 */

class OfflineManager {
    constructor() {
        this.isOnline = navigator.onLine;
        this.pendingRequests = [];
        this.swRegistration = null;
        this.onlineListeners = [];
        this.offlineListeners = [];

        this.init();
    }

    /**
     * Initialize offline manager
     */
    init() {
        // Listen for online/offline events
        window.addEventListener('online', () => this.handleOnline());
        window.addEventListener('offline', () => this.handleOffline());

        // Register service worker
        this.registerServiceWorker();

        // Initialize UI indicator
        this.updateUI();

        // Load any pending requests from localStorage
        this.loadPendingRequests();

        console.log('[OfflineManager] Initialized. Online:', this.isOnline);
    }

    /**
     * Register the service worker
     */
    async registerServiceWorker() {
        if ('serviceWorker' in navigator) {
            try {
                this.swRegistration = await navigator.serviceWorker.register(
                    '/service-worker.js',  // Served from root with Service-Worker-Allowed header
                    { scope: '/' }
                );

                console.log('[OfflineManager] Service worker registered:', this.swRegistration.scope);

                // Handle updates
                this.swRegistration.addEventListener('updatefound', () => {
                    const newWorker = this.swRegistration.installing;
                    newWorker.addEventListener('statechange', () => {
                        if (newWorker.state === 'installed' && navigator.serviceWorker.controller) {
                            this.showUpdateNotification();
                        }
                    });
                });
            } catch (error) {
                console.error('[OfflineManager] Service worker registration failed:', error);
            }
        } else {
            console.warn('[OfflineManager] Service workers not supported');
        }
    }

    /**
     * Handle going online
     */
    handleOnline() {
        console.log('[OfflineManager] Back online');
        this.isOnline = true;
        this.updateUI();

        // Process pending requests
        this.processPendingRequests();

        // Notify listeners
        this.onlineListeners.forEach(fn => fn());
    }

    /**
     * Handle going offline
     */
    handleOffline() {
        console.log('[OfflineManager] Gone offline');
        this.isOnline = false;
        this.updateUI();

        // Notify listeners
        this.offlineListeners.forEach(fn => fn());
    }

    /**
     * Add listener for online event
     */
    onOnline(callback) {
        this.onlineListeners.push(callback);
    }

    /**
     * Add listener for offline event
     */
    onOffline(callback) {
        this.offlineListeners.push(callback);
    }

    /**
     * Queue a request for later execution when back online
     */
    queueRequest(url, options, description) {
        const request = {
            id: Date.now() + Math.random().toString(36).substr(2, 9),
            url,
            options,
            description,
            timestamp: new Date().toISOString()
        };

        this.pendingRequests.push(request);
        this.savePendingRequests();
        this.updateUI();

        console.log('[OfflineManager] Queued request:', description);
        return request.id;
    }

    /**
     * Process pending requests when back online
     */
    async processPendingRequests() {
        if (this.pendingRequests.length === 0) return;

        console.log('[OfflineManager] Processing', this.pendingRequests.length, 'pending requests');

        const requests = [...this.pendingRequests];
        this.pendingRequests = [];
        this.savePendingRequests();

        for (const request of requests) {
            try {
                const response = await fetch(request.url, request.options);
                if (response.ok) {
                    console.log('[OfflineManager] Request succeeded:', request.description);
                    this.showToast(`✅ ${request.description}`, 'success');
                } else {
                    throw new Error(`HTTP ${response.status}`);
                }
            } catch (error) {
                console.error('[OfflineManager] Request failed:', request.description, error);
                // Re-queue failed requests
                this.pendingRequests.push(request);
                this.savePendingRequests();
            }
        }

        this.updateUI();
    }

    /**
     * Save pending requests to localStorage
     */
    savePendingRequests() {
        try {
            localStorage.setItem('sfppy_pending_requests', JSON.stringify(this.pendingRequests));
        } catch (e) {
            console.warn('[OfflineManager] Could not save pending requests:', e);
        }
    }

    /**
     * Load pending requests from localStorage
     */
    loadPendingRequests() {
        try {
            const saved = localStorage.getItem('sfppy_pending_requests');
            if (saved) {
                this.pendingRequests = JSON.parse(saved);
            }
        } catch (e) {
            console.warn('[OfflineManager] Could not load pending requests:', e);
        }
    }

    /**
     * Update UI to reflect online/offline state
     */
    updateUI() {
        // Update indicator
        const indicator = document.getElementById('offline-indicator');
        if (indicator) {
            if (this.isOnline) {
                indicator.classList.add('hidden');
            } else {
                indicator.classList.remove('hidden');
                const badge = indicator.querySelector('.offline-badge');
                const count = indicator.querySelector('.offline-count');
                if (badge) badge.textContent = '📴 Offline Mode';
                if (count && this.pendingRequests.length > 0) {
                    count.textContent = `Queued requests: ${this.pendingRequests.length}`;
                    count.classList.remove('hidden');
                } else if (count) {
                    count.classList.add('hidden');
                }
            }
        }

        // Add/remove body class for CSS targeting
        document.body.classList.toggle('is-offline', !this.isOnline);
        document.body.classList.toggle('is-online', this.isOnline);
    }

    /**
     * Show update notification
     */
    showUpdateNotification() {
        const notification = document.createElement('div');
        notification.className = 'fixed bottom-4 right-4 bg-blue-600 text-white px-4 py-3 rounded-lg shadow-lg z-50 flex items-center gap-3';
        notification.innerHTML = `
            <span>🔄 New version available</span>
            <button onclick="window.location.reload()" class="bg-white text-blue-600 px-3 py-1 rounded text-sm font-medium hover:bg-blue-50">
                Refresh
            </button>
            <button onclick="this.parentElement.remove()" class="text-blue-200 hover:text-white">✕</button>
        `;
        document.body.appendChild(notification);
    }

    /**
     * Show toast notification
     */
    showToast(message, type = 'info') {
        // Use existing toast function if available
        if (window.showToast) {
            window.showToast(message, type);
            return;
        }

        const toast = document.createElement('div');
        const bgColor = {
            success: 'bg-green-600',
            error: 'bg-red-600',
            warning: 'bg-yellow-600',
            info: 'bg-blue-600'
        }[type] || 'bg-gray-600';

        toast.className = `fixed bottom-4 left-1/2 transform -translate-x-1/2 ${bgColor} text-white px-4 py-2 rounded-lg shadow-lg z-50`;
        toast.textContent = message;
        document.body.appendChild(toast);

        setTimeout(() => toast.remove(), 3000);
    }

    /**
     * Get cache status from service worker
     */
    async getCacheStatus() {
        if (!navigator.serviceWorker.controller) {
            return null;
        }

        return new Promise((resolve) => {
            const messageChannel = new MessageChannel();
            messageChannel.port1.onmessage = (event) => {
                resolve(event.data);
            };
            navigator.serviceWorker.controller.postMessage(
                { type: 'GET_CACHE_STATUS' },
                [messageChannel.port2]
            );
        });
    }

    /**
     * Clear all caches
     */
    async clearCache() {
        if (!navigator.serviceWorker.controller) {
            return false;
        }

        return new Promise((resolve) => {
            const messageChannel = new MessageChannel();
            messageChannel.port1.onmessage = (event) => {
                resolve(event.data.cleared);
            };
            navigator.serviceWorker.controller.postMessage(
                { type: 'CLEAR_CACHE' },
                [messageChannel.port2]
            );
        });
    }

    /**
     * Pre-cache a substance
     */
    cacheSubstance(cid) {
        if (navigator.serviceWorker.controller) {
            navigator.serviceWorker.controller.postMessage({
                type: 'CACHE_SUBSTANCE',
                cid: cid
            });
        }
    }

    /**
     * Get structure image URL with fallback
     */
    getStructureImageUrl(cid) {
        if (!cid) return '';

        // Use local endpoint first - service worker handles fallback
        return `/api/substances/structure/${cid}.png`;
    }

    /**
     * Check if a feature is available offline
     */
    isFeatureAvailable(feature) {
        const offlineFeatures = [
            'viewJobs',
            'loadExamples',
            'viewCachedSubstances',
            'editSession',
            'darkMode',
            'help'
        ];

        const onlineOnlyFeatures = [
            'runSimulation',
            'searchPubChem',
            'runFitting',
            'exportPdf'
        ];

        if (this.isOnline) return true;

        return offlineFeatures.includes(feature) && !onlineOnlyFeatures.includes(feature);
    }
}

// Create global instance
window.offlineManager = new OfflineManager();

// Export for module usage
if (typeof module !== 'undefined' && module.exports) {
    module.exports = OfflineManager;
}
