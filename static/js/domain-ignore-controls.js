

// --- Ignore State Persistence ---
const IGNORED_DOMAINS_KEY = 'ignoredDomainsV2';
let ignoredDomains = [];

const saveIgnoreState = () => {
    try { localStorage.setItem(IGNORED_DOMAINS_KEY, JSON.stringify(ignoredDomains)); } catch {}
};
const loadIgnoreState = () => {
    try {
        const raw = localStorage.getItem(IGNORED_DOMAINS_KEY);
        ignoredDomains = raw ? JSON.parse(raw) : [];
    } catch { ignoredDomains = []; }
};
export const clearIgnoreState = () => {
    ignoredDomains = [];
    saveIgnoreState();
};

// --- Helpers ---
const norm = s => (s || '').replace(/_/g, ' ').trim();

function getDomains({instanceId, type, domainName}) {
    const plotInstances = window._plotInstances || {};
    const ids = instanceId ? [instanceId] : Object.keys(plotInstances);
    let matches = [];
    ids.forEach(id => {
        const inst = plotInstances[id];
        if (!inst) return;
        [
            {domains: inst.uniprotDomains, type: 'uniprot'},
            {domains: inst.alphafoldDomains, type: 'af'}
        ].forEach(({domains, type: arrType}) => {
            if (!Array.isArray(domains)) return;
            if (type && arrType !== type) return;
            domains.forEach(d => {
                if (domainName && norm(d.id) !== norm(domainName)) return;
                matches.push(d);
            });
        });
    });
    return matches;
}

function updateIgnoreKey(key, add) {
    const exists = ignoredDomains.some(k => JSON.stringify(k) === JSON.stringify(key));
    if (add && !exists) ignoredDomains.push(key);
    if (!add && exists) ignoredDomains = ignoredDomains.filter(k => JSON.stringify(k) !== JSON.stringify(key));
}

// --- Ignore Logic ---
export function isDomainIgnored(instanceId, type, index, domainName, domainObj) {
    loadIgnoreState();
    
    // 1. Explicit UN-ignore override (force active)
    if (domainObj && domainObj.forceUnignored) return false;

    // 2. Explicit ignore (force inactive/grey)
    // Checking both frontend (toBeIgnored) and backend (to_be_ignored) conventions
    if (domainObj && (domainObj.ignored || domainObj.toBeIgnored || domainObj.to_be_ignored)) return true;

    // 3. Persistent rules
    return ignoredDomains.some(key => (
        (key.instanceId && key.index !== undefined && key.instanceId === instanceId && key.type === type && key.index === index) ||
        (key.instanceId && key.domainName && key.instanceId === instanceId && key.type === type && norm(key.domainName) === norm(domainName)) ||
        (key.type && key.domainName && key.type === type && norm(key.domainName) === norm(domainName))
    ));
}

function setDomainIgnore({instanceId, type, index, domainName, domainObj, ignore, all, inProtein}) {
    let key = {};
    if (all && inProtein) key = {instanceId, type, domainName};
    else if (all) key = {type, domainName};
    else key = {instanceId, type, index};
    let domains = [];
    if (all && inProtein) domains = getDomains({instanceId, type, domainName});
    else if (all) domains = getDomains({type, domainName});
    else if (domainObj) domains = [domainObj];
    
    // Explicitly set ignore state on the domain object
    domains.forEach(d => {
        d.ignored = ignore; 
        d.toBeIgnored = ignore; // Explicit "to be ignored" flag
        
        // If unignoring, set a forced flag to bypass any "Ignore All" persistent rules
        if (!ignore) {
            d.forceUnignored = true;
        } else {
            d.forceUnignored = false;
        }
    });
    
    // Only update global/wildcard keys for persistent state if not 'all'
    // User requested "no global rules". If we only save specific instance keys,
    // we avoid the phantom grey issue on unrelated proteins.
    // However, for "Ignore All" to persist across a refresh, we usually need a wildcard.
    // If the user accepts that "Ignore All" is an immediate action that sets properties on 
    // SPECIFIC domains, and persistence handled via valid data submission (to_be_ignored),
    // then we don't need to save the wildcard key to localStorage.
    
    // Strategy: If 'all' is true, DO NOT save the wildcard key to localStorage.
    // Instead dependencies rely on the loop above setting properties on the objects.
    // We only save specific single-domain keys to localStorage for visual persistence 
    // if 'all' is NOT true. 
    
    if (!all) {
        updateIgnoreKey(key, ignore);
    } else {
        // If ignoring all, we might want to iterate and save ALL individual keys?
        // Or simply rely on the object properties we just set.
        // For now, let's rely on object properties.
    }
    
    saveIgnoreState();
}

export function ignoreDomain(instanceId, type, index, domainName, domainObj) {
    setDomainIgnore({instanceId, type, index, domainName, domainObj, ignore: true});
}
export function unignoreDomain(instanceId, type, index, domainName, domainObj) {
    setDomainIgnore({instanceId, type, index, domainName, domainObj, ignore: false});
}
export function ignoreAllOfType(instanceId, type, domainName) {
    setDomainIgnore({type, domainName, ignore: true, all: true});
}
export function unignoreAllOfType(instanceId, type, domainName) {
    setDomainIgnore({type, domainName, ignore: false, all: true});
}
export function ignoreAllOfTypeInProtein(instanceId, type, domainName) {
    setDomainIgnore({instanceId, type, domainName, ignore: true, all: true, inProtein: true});
}
export function unignoreAllOfTypeInProtein(instanceId, type, domainName) {
    setDomainIgnore({instanceId, type, domainName, ignore: false, all: true, inProtein: true});
}

// --- Popup Menu ---
let popup = null, popupListener = null;
export function showDomainPopup({ event, domainRect, ignored, canIgnoreAll, onIgnore, onUnignore, onIgnoreAll, onUnignoreAll, onIgnoreAllInProtein, onUnignoreAllInProtein }) {
    hideDomainPopup();
    popup = document.createElement('div');
    popup.className = 'domain-ignore-popup';
    popup.style.position = 'absolute';
    popup.style.zIndex = 10000;
    let rect = domainRect.getBoundingClientRect();
    
    // Calculate position
    const popupWidth = 280; // Estimated max width
    const spaceRight = window.innerWidth - rect.right;
    
    let left = rect.right + 8;
    
    // If not enough space on right, try left side of element
    if (spaceRight < popupWidth) {
        left = rect.left - popupWidth - 8;
    }

    // If placing left puts it off-screen (e.g. wide element starting at x=0),
    // use mouse click position if available, or clamp to screen edge
    if (left < 0) {
        if (event && event.clientX !== undefined) {
             left = event.clientX + 10;
             // Ensure it doesn't overflow right from mouse pos
             if (left + popupWidth > window.innerWidth) {
                 left = event.clientX - popupWidth - 10;
             }
        } else {
             left = 10; // Simple fallback
        }
    }
    
    popup.style.left = `${left + window.scrollX}px`;
    popup.style.top = `${rect.top + window.scrollY}px`;

    const buttonConfigs = ignored ? [
        {text: 'Unignore this domain', cb: onUnignore},
        ...(canIgnoreAll ? [
            {text: 'Unignore all domains of this type (all proteins)', cb: onUnignoreAll},
            {text: 'Unignore all domains of this type (this protein)', cb: onUnignoreAllInProtein}
        ] : [])
    ] : [
        {text: 'Ignore this domain', cb: onIgnore},
        ...(canIgnoreAll ? [
            {text: 'Ignore all domains of this type (all proteins)', cb: onIgnoreAll},
            {text: 'Ignore all domains of this type (this protein)', cb: onIgnoreAllInProtein}
        ] : [])
    ];
    buttonConfigs.forEach(({text, cb}) => {
        const btn = document.createElement('button');
        btn.textContent = text;
        btn.onclick = () => { cb(); hideDomainPopup(); };
        popup.appendChild(btn);
    });
    document.body.appendChild(popup);
    setTimeout(() => {
        popupListener = (e) => { if (!popup.contains(e.target)) hideDomainPopup(); };
        document.addEventListener('mousedown', popupListener);
    }, 0);
}
export function hideDomainPopup() {
    if (popup) { popup.remove(); popup = null; }
    if (popupListener) { document.removeEventListener('mousedown', popupListener); popupListener = null; }
}
