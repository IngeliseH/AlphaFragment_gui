import { isDomainIgnored } from './domain-ignore-controls.js';
import { getApprovedState } from './approve-controls.js';

function parseInstanceIdFromCheckboxId(id) {
    // Expected: p<index>-approve-checkbox
    if (!id || typeof id !== 'string') return null;
    const m = id.match(/^(p\d+)-approve-checkbox$/);
    return m ? m[1] : null;
}

function parseIndexFromInstanceId(instanceId) {
    const m = String(instanceId || '').match(/^p(\d+)$/);
    if (!m) return null;
    const oneBased = parseInt(m[1], 10);
    if (!Number.isFinite(oneBased) || oneBased <= 0) return null;
    return oneBased - 1;
}

document.addEventListener('DOMContentLoaded', () => {
    const form = document.getElementById('refragment-form');
    if (!form) return;

    // No longer using separate approved_indices field.

    const finishButton = document.getElementById('finish-generate-button');

    const unapprovedButton = document.getElementById('refragment-unapproved-button');
    const allButton = document.getElementById('refragment-all-button');

    const updateRefragmentButtonWidths = () => {
        if (!unapprovedButton || !allButton) return;

        // Reset to natural size before measuring.
        unapprovedButton.style.width = '';
        allButton.style.width = '';

        // Measure after layout.
        const measure = () => {
            const w1 = unapprovedButton.getBoundingClientRect().width;
            const w2 = allButton.getBoundingClientRect().width;
            const max = Math.ceil(Math.max(w1, w2));
            if (!Number.isFinite(max) || max <= 0) return;
            unapprovedButton.style.width = `${max}px`;
            allButton.style.width = `${max}px`;
        };

        requestAnimationFrame(() => {
            requestAnimationFrame(measure);
        });
    };

    const updateFinishButtonState = () => {
        if (!finishButton) return;
        const checkboxes = document.querySelectorAll('input.approve-checkbox');
        const total = checkboxes.length;
        if (total === 0) {
            finishButton.disabled = true;
            return;
        }
        let allApproved = true;
        checkboxes.forEach((cb) => {
            if (!(cb instanceof HTMLInputElement)) return;
            if (!cb.checked) allApproved = false;
        });
        finishButton.disabled = !allApproved;
    };

    // Update as approvals change (manual checkbox clicks).
    document.addEventListener('change', (ev) => {
        const target = ev.target;
        if (target && target instanceof HTMLInputElement && target.classList.contains('approve-checkbox')) {
            updateFinishButtonState();
        }
    });

    // Update when approval state is changed programmatically (e.g. bulk-approve).
    document.addEventListener('approvalstatechange', () => {
        updateFinishButtonState();
    });


    form.addEventListener('submit', (e) => {
        const submitter = e.submitter;
        const action = submitter?.getAttribute?.('value') || '';

        // Prevent default submission to handle data processing manually
        e.preventDefault();

        // 1. Filter data and update ignore/approval states
        if (window.proteinDataMap) {
            console.log('[Refragment] Processing protein data before submission...');
            const filteredProteinData = [];
            const plotInstances = window._plotInstances || {};
            let proteinIndex = 0;

            for (const [proteinName, proteinData] of window.proteinDataMap.entries()) {
                
                let instanceId = null;
                for (const [id, inst] of Object.entries(plotInstances)) {
                    if (inst.proteinName === proteinName) {
                        instanceId = id;
                        break;
                    }
                }

                const isApproved = instanceId ? getApprovedState(instanceId) : false;

                // Prepare to process data even if approved, to capture ignore states
                if (!instanceId) {
                    console.log(`[Refragment] No instanceId found for protein: ${proteinName}`);
                    filteredProteinData.push({ ...proteinData, name: proteinName, isApproved: isApproved });
                    proteinIndex++;
                    continue;
                }

                // Function to process domains and attach flags based on current plot state
                const { alphafoldDomains = [], uniprotDomains = [] } = proteinData;
                const processDomains = (domains, type) => {
                    if (!Array.isArray(domains)) return [];
                    return domains.map((d, idx) => {
                         const ignored = isDomainIgnored(instanceId, type, idx, d.id, d);
                         
                         // Logic for APPROVED proteins:
                         // User wants: "domains in approved proteins should be able to be set 'to be ignored'"
                         // If action == refragment_unapproved:
                         //   - if 'ignored' (grey/toBeIgnored) -> keep as 'toBeIgnored' (so it stays grey)
                         //   - if 'not ignored' -> stays not ignored
                         // If action == refragment_all:
                         //   - if 'ignored' -> becomes 'currently_ignored' (hidden) via server logic
                         
                         // General logic: always capture the 'toBeIgnored' flag.
                         // The server decides what to do with it based on 'currently_ignored' transitions.
                         
                         if (ignored) {
                             return { ...d, toBeIgnored: true };
                         }
                         if (!ignored && d.currently_ignored) {
                             return { ...d, currently_ignored: false };
                         }
                         return d;
                    });
                };

                const filteredAF = processDomains(alphafoldDomains, 'af');
                const filteredUni = processDomains(uniprotDomains, 'uniprot');
                
                // If it is approved, we just push the data with updated flags.
                // The server will see it is in 'approved_set' and skip re-fragmentation,
                // BUT it will return the updated domain flags so UI reflects them.
                
                filteredProteinData.push({
                    ...proteinData,
                    name: proteinName,
                    alphafoldDomains: filteredAF,
                    uniprotDomains: filteredUni,
                    isApproved: isApproved
                });
                
                console.log(`[Refragment] Processed ${proteinName} (Approved: ${isApproved})`);
                proteinIndex++;
            }

            // Update hidden field
            const currentDataField = document.getElementById('current_protein_data');
            if (currentDataField) {
                currentDataField.value = JSON.stringify(filteredProteinData);
            }
        }

        // Ensure action is sent when submitting programmatically
        let actionInput = form.querySelector('input[name="action"][type="hidden"]');
        if (!actionInput) {
            actionInput = document.createElement('input');
            actionInput.type = 'hidden';
            actionInput.name = 'action';
            form.appendChild(actionInput);
        }
        actionInput.value = action;

        // Now submit the form programmatically
        form.submit();
    });

    // Initial state
    updateFinishButtonState();

    // Keep refragment buttons equal width (compact).
    updateRefragmentButtonWidths();
    window.addEventListener('resize', updateRefragmentButtonWidths);
});
