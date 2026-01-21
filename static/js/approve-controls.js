const DEFAULT_APPROVE_COLUMN_WIDTH_PX = 28;

const BULK_LABELS = {
    within: 'Approve all fragmentations within parameters',
    close: 'Approve all fragmentations close to meeting parameters',
    all: 'Approve all',
};

// module-level approval state
const approvedStates = {};

export function getApprovedState(instanceId) {
    return Boolean(approvedStates[instanceId]);
}

export function setApprovedState(instanceId, isApproved) {
    approvedStates[instanceId] = Boolean(isApproved);
}

export function createApproveCheckbox({ instanceId, checked = false, onChange }) {
    const checkbox = document.createElement('input');
    checkbox.type = 'checkbox';
    checkbox.className = 'approve-checkbox';
    checkbox.id = `${instanceId}-approve-checkbox`;
    checkbox.title = 'Approve';
    checkbox.setAttribute('aria-label', 'Approve');
    checkbox.checked = Boolean(checked);

    checkbox.addEventListener('change', () => {
        onChange?.(Boolean(checkbox.checked), checkbox);
    });

    return checkbox;
}

export function ensureBulkApproveControls({
    plotSectionsContainer,
    approveColumnWidthPx = DEFAULT_APPROVE_COLUMN_WIDTH_PX,
    getInstanceIds,
    getInstanceData,
    applyApproved,
}) {
    if (!plotSectionsContainer) return;
    if (document.getElementById('approve-controls')) return;

    const root = document.createElement('div');
    root.id = 'approve-controls';
    root.className = 'approve-controls';

    const approveCol = document.createElement('div');
    approveCol.className = 'approve-controls__col';
    approveCol.style.width = `${approveColumnWidthPx}px`;
    approveCol.style.flex = `0 0 ${approveColumnWidthPx}px`;

    const labelsCol = document.createElement('div');
    labelsCol.className = 'approve-controls__labels';

    const labelsByKey = {};

    const formatLabel = (key, counts) => {
        const base = BULK_LABELS[key];
        const n = counts?.[key];
        return Number.isFinite(n) ? `${base} (${n})` : base;
    };

    const computeCounts = () => {
        const ids = getInstanceIds?.() || [];
        let errors = 0;
        let success = 0;
        let close = 0;

        ids.forEach((instanceId) => {
            const data = getInstanceData?.(instanceId) || {};
            const status = classifyProtein(data);
            if (status === 'unknown') {
                errors += 1;
                return;
            }

            if (status === 'satisfactory' || status === 'too_short') success += 1;
            if (status === 'satisfactory' || status === 'slightly_out' || status === 'too_short') close += 1;
        });

        return {
            within: success,
            close,
            all: Math.max(0, ids.length - errors),
        };
    };

    const updateLabels = () => {
        const counts = computeCounts();
        ['within', 'close', 'all'].forEach((key) => {
            const el = labelsByKey[key];
            if (!el) return;
            el.textContent = formatLabel(key, counts);
        });
    };

    const scheduleUpdateLabels = () => {
        // Plots/instances are created right after this function returns.
        // Defer to the next paint so getInstanceIds() reflects the final list.
        requestAnimationFrame(() => updateLabels());
    };

    const addAction = (key) => {
        const cb = document.createElement('input');
        cb.type = 'checkbox';
        cb.className = 'approve-controls__action';
        cb.id = `bulk-approve-${key}`;
        cb.title = BULK_LABELS[key];
        cb.setAttribute('aria-label', BULK_LABELS[key]);

        const label = document.createElement('div');
        label.className = 'approve-controls__label';
        label.textContent = BULK_LABELS[key];
        labelsByKey[key] = label;

        cb.addEventListener('change', () => {
            if (!cb.checked) return;
            bulkApproveByRule({
                ruleKey: key,
                getInstanceIds,
                getInstanceData,
                applyApproved,
            });
            // action checkbox: reset
            cb.checked = false;
        });

        approveCol.appendChild(cb);
        labelsCol.appendChild(label);
    };

    addAction('within');
    addAction('close');
    addAction('all');

    root.appendChild(approveCol);
    root.appendChild(labelsCol);
    plotSectionsContainer.prepend(root);

    document.addEventListener('plotinstancesready', scheduleUpdateLabels);
    scheduleUpdateLabels();
}

function getFragmentParams() {
    const params = window.fragmentParams || {};
    const length = params.length || {};
    const overlap = params.overlap || {};
    const min = Number(length.min);
    const max = Number(length.max);
    const ovMin = Number(overlap.min);
    const ovMax = Number(overlap.max);

    return {
        min,
        max,
        maxSlight: Number.isFinite(max) ? max * 1.2 : NaN,
        ovMin,
        ovMax,
    };
}

function classifyProtein({ proteinLength, fragmentIndices }) {
    const { min, max, maxSlight, ovMin, ovMax } = getFragmentParams();

    const length = Number(proteinLength);
    if (!Number.isFinite(length) || !Number.isFinite(min) || !Number.isFinite(max) || !Number.isFinite(ovMin) || !Number.isFinite(ovMax)) {
        return 'unknown';
    }

    if (length <= max) return 'too_short';

    const fragments = Array.isArray(fragmentIndices) ? fragmentIndices : [];
    if (fragments.length === 0) return 'significantly_out';

    let within = true;
    let anyAboveMax = false;
    let anyAboveSlight = false;
    let anyBelowMin = false;
    let anyOverlapOutside = false;

    for (const frag of fragments) {
        if (!Array.isArray(frag) || frag.length !== 2) return 'significantly_out';
        const start = parseInt(frag[0], 10);
        const end = parseInt(frag[1], 10);
        if (!Number.isFinite(start) || !Number.isFinite(end) || end < start) return 'significantly_out';
        const fragLen = end - start + 1;

        if (fragLen < min) {
            anyBelowMin = true;
            within = false;
        }
        if (fragLen > max) {
            anyAboveMax = true;
            within = false;
        }
        if (Number.isFinite(maxSlight) && fragLen > maxSlight) {
            anyAboveSlight = true;
        }
    }

    if (fragments.length >= 2) {
        for (let i = 0; i < fragments.length - 1; i++) {
            const prevEnd = parseInt(fragments[i][1], 10);
            const nextStart = parseInt(fragments[i + 1][0], 10);
            if (!Number.isFinite(prevEnd) || !Number.isFinite(nextStart)) return 'significantly_out';
            const overlapLen = prevEnd - nextStart + 1;
            if (overlapLen < ovMin || overlapLen > ovMax) {
                anyOverlapOutside = true;
                within = false;
            }
        }
    }

    if (within) return 'satisfactory';

    // Slightly out: only issue is fragment length marginally above max (<= 1.2x max)
    if (anyAboveMax && !anyAboveSlight && !anyBelowMin && !anyOverlapOutside) {
        return 'slightly_out';
    }

    return 'significantly_out';
}

function bulkApproveByRule({ ruleKey, getInstanceIds, getInstanceData, applyApproved }) {
    const ids = getInstanceIds?.() || [];
    ids.forEach((instanceId) => {
        const data = getInstanceData?.(instanceId);
        const status = classifyProtein(data || {});

        const shouldApprove = (() => {
            if (ruleKey === 'all') return true;
            if (ruleKey === 'close') return status === 'satisfactory' || status === 'slightly_out' || status === 'too_short';
            if (ruleKey === 'within') return status === 'satisfactory' || status === 'too_short';
            return false;
        })();

        if (shouldApprove) applyApproved?.(instanceId, true);
    });
}
