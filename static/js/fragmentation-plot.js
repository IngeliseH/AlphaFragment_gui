import {
    createApproveCheckbox,
    ensureBulkApproveControls,
    getApprovedState,
    setApprovedState,
} from './approve-controls.js';
import {
    isDomainIgnored,
    ignoreDomain,
    unignoreDomain,
    ignoreAllOfType,
    unignoreAllOfType,
    ignoreAllOfTypeInProtein,
    unignoreAllOfTypeInProtein,
    showDomainPopup,
    hideDomainPopup
} from './domain-ignore-controls.js';

let _plotInstances = {};
window._plotInstances = _plotInstances;

// Keep a consistent coordinate system across all proteins so plots don't
// change scale depending on page scrollbars or container timing.
const _PLOT_VIEWBOX_WIDTH = 1200;
const _PLOT_VIEWBOX_HEIGHT = 125;

// Remove grey from palette (was '#C2C1B1')
const _domainColors = ['#ED7AB0', '#8DC640', '#F68B1F', '#9F83BC', "#FFDD55", '#6DC8BF',
    '#A74399', '#A6ADD3', '#E64425', '#00A45D', '#BA836E', '#3E4291'];
const _domainColorMap = {};
let _currentColorIndex = 0;

// =============================================================================
// Public API Functions
// =============================================================================
export function initializeDomainFragmentPlots({ proteins, interactionRegions = [] }) {
    const fallbackMessage = document.getElementById('domain-fragments-fallback-message');
    if (fallbackMessage) fallbackMessage.style.display = 'none';

    // Use the proteinDataMap from the window (set up in the template)
    const proteinDataMap = window.proteinDataMap;

    const plotSectionsContainer = document.getElementById('domain-fragment-plots-container');
    ensureBulkApproveControls({
        plotSectionsContainer,
        getInstanceIds: () => Object.keys(_plotInstances),
        getInstanceData: (instanceId) => {
            const instance = _plotInstances[instanceId];
            if (!instance) return null;
            return {
                proteinLength: instance.proteinLength,
                fragmentIndices: instance.fragmentIndices,
            };
        },
        applyApproved: (instanceId, isApproved) => {
            _applyApprovedState(instanceId, isApproved);
        },
    });

    for (let i = 1; i <= proteins.length; i++) {
        const proteinName = proteins[i-1];
        const interactionRegion = (interactionRegions.length >= i) ? interactionRegions[i-1] : [];
        const proteinData = proteinDataMap ? proteinDataMap.get(proteinName) : null;
        _initializePlotInstance({
            instanceId: `p${i}`,
            proteinName,
            proteinData,
            interactionRegion,
            plotSectionsContainer
        });
    }

    document.dispatchEvent(new CustomEvent('plotinstancesready', {
        detail: { instanceIds: Object.keys(_plotInstances) },
    }));
}

// =============================================================================
// Core Logic
// =============================================================================
function _initializePlotInstance({instanceId, proteinName, proteinData, interactionRegion = [], plotSectionsContainer}) {
    const section = document.createElement('div');
    section.id = `domain-fragment-plot-${instanceId}-section`;

    // Restore approved state if passed from server
    if (proteinData && proteinData.isApproved) {
        setApprovedState(instanceId, true);
    }

    // Layout: [approve checkbox column] [content column]
    section.className = 'approve-row';

    const approveCol = document.createElement('div');
    approveCol.className = 'approve-row__approve-col';

    const contentCol = document.createElement('div');
    contentCol.className = 'approve-row__content';

    const subheading = document.createElement('h3');
    subheading.className = 'page-subtitle';
    subheading.textContent = proteinName;
    subheading.style.textAlign = 'left';
    subheading.style.margin = '0 0 8px 0';
    subheading.style.lineHeight = '1.15';

    const approveCheckbox = createApproveCheckbox({
        instanceId,
        checked: getApprovedState(instanceId),
        onChange: (checked) => {
            _applyApprovedState(instanceId, checked);
        },
    });

    const container = document.createElement('div');
    container.id = `domain-fragment-plot-container-${instanceId}`;
    container.className = 'domain-fragment-plot-container';
    container.style.minHeight = `${_PLOT_VIEWBOX_HEIGHT}px`;
    _displayInfo(container, "Loading data for domain/fragment plot...");

    approveCol.appendChild(approveCheckbox);
    contentCol.appendChild(subheading);
    contentCol.appendChild(container);
    section.appendChild(approveCol);
    section.appendChild(contentCol);
    plotSectionsContainer.appendChild(section);

    const instance = {
        proteinName,
        proteinLength: null,
        fragmentIndices: null,
        alphafoldDomains: null,
        uniprotDomains: null,
        interactionRegion: interactionRegion || [],
        isCollapsibleTableCollapsed: true,
        isApproved: getApprovedState(instanceId),
        approveCheckboxId: approveCheckbox.id,
        containerSelector: container.id,
    };
    _plotInstances[instanceId] = instance;

    if (!proteinData || proteinData.length === null || isNaN(proteinData.length)) {
        const message = proteinName ? `Length data not available or invalid for ${proteinName}.` : 'Protein not specified.';
        _displayInfo(container, message, true);
        return;
    }

    instance.proteinLength = proteinData.length;
    instance.fragmentIndices = proteinData.fragmentIndices;
    instance.alphafoldDomains = proteinData.alphafoldDomains;
    instance.uniprotDomains = proteinData.uniprotDomains;

    container.innerHTML = '';
    _renderPlot(container, instanceId);
    _renderCollapsibleTable(container, instanceId);

    // Apply initial collapsed state if previously approved.
    _setInstanceCollapsed(instanceId, instance.isApproved);
}

function _updatePlot(instanceId) {
    const instance = _plotInstances[instanceId];
    if (!instance) return;

    const container = document.getElementById(instance.containerSelector);
    if (!container) return;
    container.innerHTML = '';
    if (!instance.proteinName || !instance.proteinLength) {
        const message = instance.proteinName ? `Length data not available for ${instance.proteinName}.` : 'Protein not specified.';
        _displayInfo(container, message, true);
        return null;
    }

    _renderPlot(container, instanceId);
    _renderCollapsibleTable(container, instanceId);

    _setInstanceCollapsed(instanceId, Boolean(instance.isApproved));
}

function _setInstanceCollapsed(instanceId, shouldCollapse) {
    const instance = _plotInstances[instanceId];
    if (!instance) return;
    const container = document.getElementById(instance.containerSelector);
    if (!container) return;
    container.style.display = shouldCollapse ? 'none' : '';
}

function _applyApprovedState(instanceId, isApproved) {
    const instance = _plotInstances[instanceId];
    if (!instance) return;

    const approved = Boolean(isApproved);
    instance.isApproved = approved;
    setApprovedState(instanceId, approved);

    const checkbox = instance.approveCheckboxId ? document.getElementById(instance.approveCheckboxId) : null;
    if (checkbox) checkbox.checked = approved;

    _setInstanceCollapsed(instanceId, approved);

    document.dispatchEvent(new CustomEvent('approvalstatechange', {
        detail: { instanceId, approved },
    }));
}

function _renderPlot(container, instanceId) {
    const instance = _plotInstances[instanceId];
    const { proteinLength, fragmentIndices, alphafoldDomains, uniprotDomains } = instance;
    const margin = { top: 0, right: 60, bottom: 0, left: 60 };
    const dimensions = _calculatePlotDimensions(container, margin);

    const svg = _createSvgElement("svg", {
        "width": "100%",
        "height": dimensions.containerHeight,
        "viewBox": `0 0 ${_PLOT_VIEWBOX_WIDTH} ${dimensions.containerHeight}`
    });

    const svgGroup = _createSvgElement("g", {
        "transform": `translate(${margin.left}, ${margin.top})`
    });

    // Always render a baseline protein bar; domains/fragments layer on top.
    const rectY = dimensions.plotHeight / 2 - (dimensions.barHeight / 2);
    const rectHeight = dimensions.barHeight;
    const proteinBar = _createSvgElement("rect", {
        "x": 0,
        "y": rectY,
        "width": dimensions.plotWidth,
        "height": rectHeight,
        "fill": "#ccc"
    });
    svgGroup.appendChild(proteinBar);

    if (alphafoldDomains && alphafoldDomains.length > 0 && proteinLength) {
        // filter out currently_ignored for plot
        const visibleAF = alphafoldDomains.filter(d => !d.currently_ignored);
        if (visibleAF.length > 0) {
            _renderDomains(svgGroup, instanceId, {
                plotWidth: dimensions.plotWidth,
                yPosition: (dimensions.plotHeight / 2) - (dimensions.alphafoldDomainHeight / 2),
                type: 'af',
                domainHeight: dimensions.alphafoldDomainHeight,
                labelFontSize: 10,
                domainsOverride: visibleAF
            });
        }
    }

    if (uniprotDomains && uniprotDomains.length > 0 && proteinLength) {
        // filter out currently_ignored for plot
        const visibleUni = uniprotDomains.filter(d => !d.currently_ignored);
        if (visibleUni.length > 0) {
            _renderDomains(svgGroup, instanceId, {
                plotWidth: dimensions.plotWidth,
                yPosition: (dimensions.plotHeight / 2) - (dimensions.uniprotDomainHeight / 2),
                type: 'uniprot',
                domainHeight: dimensions.uniprotDomainHeight,
                labelFontSize: 10,
                domainsOverride: visibleUni
            });
        }
    }

    if (fragmentIndices && fragmentIndices.length > 0 && proteinLength) {
        _renderFragments(svgGroup, fragmentIndices, instanceId, {
            plotWidth: dimensions.plotWidth,
            yPosition: dimensions.plotHeight / 2,
            height: dimensions.fragmentBarHeight
        });
    }

    const startLabel = _createProteinLabel("1", -15, dimensions.plotHeight / 2, { textAnchor: "end" });
    svgGroup.appendChild(startLabel);

    const endLabel = _createProteinLabel(proteinLength, dimensions.plotWidth + 15, dimensions.plotHeight / 2, { textAnchor: "start" });
    svgGroup.appendChild(endLabel);

    svg.appendChild(svgGroup);
    container.appendChild(svg);
}

function _renderCollapsibleTable(container, instanceId) {
    const instance = _plotInstances[instanceId];
    const { alphafoldDomains, uniprotDomains } = instance;

    const domainInfoSection = document.createElement('div');
    domainInfoSection.className = 'collapsible-subsection';
    if (instance.isCollapsibleTableCollapsed) {
        domainInfoSection.classList.add('collapsed');
    }
    domainInfoSection.id = `domain-info-collapsible-section-${instanceId}`;

    const titleDiv = document.createElement('div');
    titleDiv.className = 'collapsible-subsection-title';

    const titleText = document.createElement('h4');
    titleText.textContent = 'Domain Details';
    titleDiv.appendChild(titleText);

    const icon = document.createElement('i');
    icon.className = 'fas fa-chevron-up';
    titleDiv.appendChild(icon);

    domainInfoSection.appendChild(titleDiv);

    const contentDiv = document.createElement('div');
    contentDiv.className = 'collapsible-subsection-content';

    const table = document.createElement('table');
    table.id = `domain-info-table-${instanceId}`;

    const thead = table.createTHead();
    const headerRow = thead.insertRow();

    const thName = document.createElement('th');
    thName.textContent = 'Name';
    headerRow.appendChild(thName);

    const thPosition = document.createElement('th');
    thPosition.textContent = 'Position';
    headerRow.appendChild(thPosition);

    const tbody = table.createTBody();

    const allDomainsForTable = [
        ...(alphafoldDomains || []).map((domain) => ({ ...domain, type: 'af', originalIndex: (alphafoldDomains || []).indexOf(domain), isVisible: !domain.currently_ignored })),
        ...(uniprotDomains || []).map((domain) => ({ ...domain, type: 'uniprot', originalIndex: (uniprotDomains || []).indexOf(domain), isVisible: !domain.currently_ignored })),
    ];
    
    // Filter out hidden domains from the table entirely
    const visibleDomainsForTable = allDomainsForTable.filter(d => d.isVisible);

    visibleDomainsForTable.sort((a, b) => a.start - b.start || a.end - b.end);

    visibleDomainsForTable.forEach(domainEntry => {
        const row = tbody.insertRow();
        const domainRectId = `${domainEntry.type}-domain-${instanceId}-${domainEntry.originalIndex}`;
        const startLabelId = `${domainEntry.type}-start-label-${instanceId}-${domainEntry.originalIndex}`;
        const endLabelId = `${domainEntry.type}-end-label-${instanceId}-${domainEntry.originalIndex}`;
        const baseIdLabelId = domainEntry.type === 'uniprot' ? `${domainEntry.type}-baseid-label-${instanceId}-${domainEntry.originalIndex}` : null;

        row.id = `${domainEntry.type}-row-${instanceId}-${domainEntry.originalIndex}`;
        row.dataset.domainType = domainEntry.type;
        row.dataset.domainIndex = domainEntry.originalIndex;

        const originalDomain = (domainEntry.type === 'af' ? alphafoldDomains : uniprotDomains)[domainEntry.originalIndex];
        const normalizedId = _normalizeDomainId(domainEntry.id);
        const ignored = isDomainIgnored(instanceId, domainEntry.type, domainEntry.originalIndex, normalizedId, originalDomain);


        const cellName = row.insertCell();
        cellName.textContent = domainEntry.type === 'af' ? 
            'AlphaFold' : 
            normalizedId;

        const cellPosition = row.insertCell();
        cellPosition.textContent = `${domainEntry.start}-${domainEntry.end}`;

        if (ignored) {
            row.classList.add('domain-row-ignored');
        }

        row.addEventListener("mouseover", () => {
            _handleDomainHover(
                true,
                instanceId,
                domainRectId,
                startLabelId,
                endLabelId,
                baseIdLabelId,
                row,
                domainEntry.start,
                domainEntry.end,
                ignored
            );
        });

        row.addEventListener("mouseout", () => {
            _handleDomainHover(
                false,
                instanceId,
                domainRectId,
                startLabelId,
                endLabelId,
                baseIdLabelId,
                row,
                domainEntry.start,
                domainEntry.end,
                ignored
            );
        });

        row.addEventListener('click', (evt) => {
            evt.stopPropagation();
            const canIgnoreAll = domainEntry.type !== 'af';
            let allDomains = [];
            if (window.proteinDataMap) {
                for (const proteinName of window.proteinDataMap.keys()) {
                    const pdata = window.proteinDataMap.get(proteinName);
                    if (pdata && Array.isArray(pdata.uniprotDomains)) {
                        allDomains.push(...pdata.uniprotDomains.map(d => ({...d, type: 'uniprot'})));
                    }
                    if (pdata && Array.isArray(pdata.alphafoldDomains)) {
                        allDomains.push(...pdata.alphafoldDomains.map(d => ({...d, type: 'af'})));
                    }
                }
            }
            const proteinDomains = [
                ...(uniprotDomains || []).map(d => ({...d, type: 'uniprot'})),
                ...(alphafoldDomains || []).map(d => ({...d, type: 'af'}))
            ];
            showDomainPopup({
                event: evt,
                domainRect: row,
                ignored,
                canIgnoreAll,
                onIgnore: () => {
                    ignoreDomain(instanceId, domainEntry.type, domainEntry.originalIndex, normalizedId, originalDomain);
                    _updatePlot(instanceId);
                },
                onUnignore: () => {
                    unignoreDomain(instanceId, domainEntry.type, domainEntry.originalIndex, normalizedId, originalDomain);
                    _updatePlot(instanceId);
                },
                onIgnoreAll: () => {
                    ignoreAllOfType(instanceId, domainEntry.type, normalizedId, allDomains);
                    Object.keys(_plotInstances).forEach(_updatePlot);
                },
                onUnignoreAll: () => {
                    unignoreAllOfType(instanceId, domainEntry.type, normalizedId, allDomains);
                    Object.keys(_plotInstances).forEach(_updatePlot);
                },
                onIgnoreAllInProtein: () => {
                    ignoreAllOfTypeInProtein(instanceId, domainEntry.type, normalizedId, proteinDomains);
                    _updatePlot(instanceId);
                },
                onUnignoreAllInProtein: () => {
                    unignoreAllOfTypeInProtein(instanceId, domainEntry.type, normalizedId, proteinDomains);
                    _updatePlot(instanceId);
                }
            });
        });
    });

    contentDiv.appendChild(table);
    domainInfoSection.appendChild(contentDiv);

    titleDiv.addEventListener('click', () => {
        domainInfoSection.classList.toggle('collapsed');
        instance.isCollapsibleTableCollapsed = !instance.isCollapsibleTableCollapsed;
    });

    container.appendChild(domainInfoSection);
}

// =============================================================================
// Plot Drawing
// =============================================================================
function _calculatePlotDimensions(container, margin) {
    const containerHeight = _PLOT_VIEWBOX_HEIGHT;
    const plotHeight = containerHeight - margin.top - margin.bottom;
    const plotWidth = _PLOT_VIEWBOX_WIDTH - margin.left - margin.right;

    const barHeight = 10;
    const fragmentBarHeight = barHeight * 1.5;
    const uniprotDomainHeight = fragmentBarHeight * 6;
    const alphafoldDomainHeight = fragmentBarHeight * 7;

    return {
        containerHeight,
        plotHeight,
        plotWidth,
        barHeight,
        fragmentBarHeight,
        uniprotDomainHeight,
        alphafoldDomainHeight,
    };
}

function _renderDomains(svgGroup, instanceId, config) {
    const { plotWidth, yPosition, type, domainHeight, domainsOverride } = config;
    const instance = _plotInstances[instanceId];
    const { proteinLength } = instance;
    const domains = domainsOverride || (type === 'uniprot' ? instance.uniprotDomains : instance.alphafoldDomains);

    window.domainPlot_domainBaseIdToColor = _domainColorMap;
    window.domainPlotInstancesData = _plotInstances;

    const fullList = type === 'uniprot' ? instance.uniprotDomains : instance.alphafoldDomains;

    domains.forEach((domain, visualIndex) => {
        // Find stable index from the full list
        const index = fullList.indexOf(domain);
        
        if (domain.start === undefined || domain.end === undefined || domain.start > domain.end) return;
        const start = Math.max(1, domain.start);
        const end = Math.min(proteinLength, domain.end);
        if (end < start) return;
        const normalizedId = _normalizeDomainId(domain.id);
        const denominator = proteinLength > 1 ? proteinLength - 1 : 1;
        const x1_orig = ((start - 1) / denominator) * plotWidth;
        const x2_orig = ((end - 1) / denominator) * plotWidth;
        const rect_x = x1_orig - 0.5;
        const domainWidth = Math.max(1, x2_orig - x1_orig + 1);
        // Use new attribute-based ignore logic
        // If instance is approved, we treat it as finalized. 
        // We only show 'ignored' (grey) if the domain object itself has a flag (to_be_ignored), 
        // but we ignore the global localStorage rules because this protein shouldn't change.
        // Also pass isApproved to isDomainIgnored context if needed, or simply logic here:
        
        let ignored = isDomainIgnored(instanceId, type, index, normalizedId, domain);
        if (instance.isApproved && !domain.ignored && !domain.toBeIgnored && !domain.to_be_ignored) {
            // If approved and not explicitly flagged on the object, override "global" ignores
            ignored = false;
        }

        const fillColor = ignored ? '#bdbdbd' : (type === 'uniprot' ? _assignDomainColor(normalizedId) : 'lightblue');

        const domainRect = _createSvgElement("rect", {
            "id": `${type}-domain-${instanceId}-${index}`,
            "x": rect_x,
            "y": yPosition,
            "width": domainWidth,
            "height": domainHeight,
            "fill": fillColor,
            "opacity": 0.6,
            "class": ignored ? 'domain-ignored' : ''
        });

        const labelOffsetVertical = 5;
        const labelOffsetHorizontal = 2;

        const startLabel = _createHoverLabel(domain.start, x1_orig - labelOffsetHorizontal, yPosition + domainHeight + labelOffsetVertical, {textAnchor: "end"});
        startLabel.id = `${type}-start-label-${instanceId}-${index}`;

        const endLabel = _createHoverLabel(domain.end, x2_orig + labelOffsetHorizontal, yPosition + domainHeight + labelOffsetVertical, {textAnchor: "start"});
        endLabel.id = `${type}-end-label-${instanceId}-${index}`;

        const labels = [startLabel, endLabel];
        if (type === 'uniprot') {
            const baseIdLabel = _createHoverLabel(normalizedId, x1_orig + (x2_orig - x1_orig) / 2, yPosition - labelOffsetVertical);
            baseIdLabel.id = `${type}-baseid-label-${instanceId}-${index}`;
            labels.push(baseIdLabel);
        }

        domainRect.addEventListener("mouseover", () => {
            _handleDomainHover(
                true,
                instanceId,
                `${type}-domain-${instanceId}-${index}`,
                `${type}-start-label-${instanceId}-${index}`,
                `${type}-end-label-${instanceId}-${index}`,
                type === 'uniprot' ? `${type}-baseid-label-${instanceId}-${index}` : null,
                document.getElementById(`${type}-row-${instanceId}-${index}`),
                domain.start,
                domain.end,
                ignored
            );
        });
        domainRect.addEventListener("mouseout", () => {
            _handleDomainHover(
                false,
                instanceId,
                `${type}-domain-${instanceId}-${index}`,
                `${type}-start-label-${instanceId}-${index}`,
                `${type}-end-label-${instanceId}-${index}`,
                type === 'uniprot' ? `${type}-baseid-label-${instanceId}-${index}` : null,
                document.getElementById(`${type}-row-${instanceId}-${index}`),
                domain.start,
                domain.end,
                ignored
            );
        });
        // Popup menu on click
        domainRect.addEventListener('click', (evt) => {
            evt.stopPropagation();
            const canIgnoreAll = type !== 'af';
            // Gather all domains for all proteins for global ignore
            let allDomains = [];
            if (window.proteinDataMap) {
                let i = 1;
                for (const proteinName of window.proteinDataMap.keys()) {
                    const pdata = window.proteinDataMap.get(proteinName);
                    if (pdata && Array.isArray(pdata.uniprotDomains)) {
                        allDomains.push(...pdata.uniprotDomains.map(d => ({...d, type: 'uniprot'})));
                    }
                    if (pdata && Array.isArray(pdata.alphafoldDomains)) {
                        allDomains.push(...pdata.alphafoldDomains.map(d => ({...d, type: 'af'})));
                    }
                    i++;
                }
            }
            // All domains for this protein
            const proteinDomains = [
                ...(instance.uniprotDomains || []).map(d => ({...d, type: 'uniprot'})),
                ...(instance.alphafoldDomains || []).map(d => ({...d, type: 'af'}))
            ];
            showDomainPopup({
                event: evt,
                domainRect,
                ignored,
                canIgnoreAll,
                onIgnore: () => {
                    ignoreDomain(instanceId, type, index, normalizedId, domain);
                    _updatePlot(instanceId);
                },
                onUnignore: () => {
                    unignoreDomain(instanceId, type, index, normalizedId, domain);
                    _updatePlot(instanceId);
                },
                onIgnoreAll: () => {
                    ignoreAllOfType(instanceId, type, normalizedId, allDomains);
                    // Update all plots
                    Object.keys(_plotInstances).forEach(_updatePlot);
                },
                onUnignoreAll: () => {
                    unignoreAllOfType(instanceId, type, normalizedId, allDomains);
                    Object.keys(_plotInstances).forEach(_updatePlot);
                },
                onIgnoreAllInProtein: () => {
                    ignoreAllOfTypeInProtein(instanceId, type, normalizedId, proteinDomains);
                    _updatePlot(instanceId);
                },
                onUnignoreAllInProtein: () => {
                    unignoreAllOfTypeInProtein(instanceId, type, normalizedId, proteinDomains);
                    _updatePlot(instanceId);
                }
            });
        });

        svgGroup.appendChild(domainRect);
        labels.forEach(label => svgGroup.appendChild(label));
    });
}

function _renderFragments(svgGroup, fragments, instanceId, config) {
    const { plotWidth, yPosition, height } = config;
    const instance = _plotInstances[instanceId];
    const { proteinLength } = instance;

    const lengthParams = window.fragmentParams && window.fragmentParams.length ? window.fragmentParams.length : null;
    const minLen = lengthParams ? Number(lengthParams.min) : NaN;
    const maxLen = lengthParams ? Number(lengthParams.max) : NaN;

    const getFragmentFill = (fragmentLength) => {
        // Fallback to previous color if params are missing/unparseable.
        if (!Number.isFinite(fragmentLength) || !Number.isFinite(minLen) || !Number.isFinite(maxLen)) {
            return "lightcoral";
        }
        if (fragmentLength < minLen) return "white";
        if (fragmentLength <= maxLen) return "lightgreen";
        if (fragmentLength <= maxLen * 1.2) return "lightyellow";
        return "lightcoral";
    };

    fragments.forEach((frag, i) => {
        if (!Array.isArray(frag) || frag.length !== 2) return;
        let [start, end] = frag;
        start = Math.max(1, parseInt(start));
        end = Math.min(proteinLength, parseInt(end));
        if (isNaN(start) || isNaN(end) || end < start) return;

        const fragmentLength = end - start + 1;

        const denominator = proteinLength > 1 ? proteinLength - 1 : 1;
        const x1_orig = ((start - 1) / denominator) * plotWidth;
        const x2_orig = ((end - 1) / denominator) * plotWidth;

        const rect_x = x1_orig - 0.5;
        const width = Math.max(1, x2_orig - x1_orig + 1);

        const yRect = i % 2 === 0 ? (yPosition - height) : yPosition;

        const originalFill = getFragmentFill(fragmentLength);
        const originalStroke = "dimgrey";
        const hoverStroke = "black";

        const fragmentRect = _createSvgElement("rect", {
            "x": rect_x,
            "y": yRect,
            "width": width,
            "height": height,
            "fill": originalFill,
            "opacity": "1",
            "stroke": originalStroke,
            "stroke-width": "0.5"
        });
        svgGroup.appendChild(fragmentRect);

        if (instance.interactionRegion.length > 0) {
            instance.interactionRegion.forEach(loc => {
                const overlapStart = Math.max(start, loc.start);
                const overlapEnd = Math.min(end, loc.end);
                if (overlapEnd < overlapStart) return;

                const x1_hl = ((overlapStart - 1) / denominator) * plotWidth;
                const x2_hl = ((overlapEnd - 1) / denominator) * plotWidth;
                const rect_x_hl = x1_hl - 0.5;
                const highlightWidth = Math.max(1, x2_hl - x1_hl + 1);
                const highlightRect = _createSvgElement('rect', {
                    "x": rect_x_hl,
                    "y": yRect,
                    "width": highlightWidth,
                    "height": height,
                    "fill": "#ff2a00",
                    "opacity": "1",
                    "stroke": "#b30000",
                    "stroke-width": "1.2",
                });
                svgGroup.appendChild(highlightRect);
            });
        }

        const labelYPos = yRect + height / 2;
        const labelOffset = 1;

        const fragStartLabel = _createHoverLabel(start, x1_orig - labelOffset, labelYPos, {textAnchor: "end"});
        svgGroup.appendChild(fragStartLabel);

        const fragEndLabel = _createHoverLabel(end, x2_orig + labelOffset, labelYPos, {textAnchor: "start"});
        svgGroup.appendChild(fragEndLabel);

        fragmentRect.addEventListener("mouseover", () => {
            fragmentRect.setAttribute("stroke", hoverStroke);
            fragStartLabel.setAttribute("visibility", "visible");
            fragEndLabel.setAttribute("visibility", "visible");
        });
        fragmentRect.addEventListener("mouseout", () => {
            fragmentRect.setAttribute("stroke", originalStroke);
            fragStartLabel.setAttribute("visibility", "hidden");
            fragEndLabel.setAttribute("visibility", "hidden");
        });
    });
}

// =============================================================================
// Internal Helpers
// =============================================================================
function _assignDomainColor(baseId) {
    if (!_domainColorMap[baseId]) {
        _domainColorMap[baseId] = _domainColors[_currentColorIndex % _domainColors.length];
        _currentColorIndex++;
    }
    return _domainColorMap[baseId];
}

function _normalizeDomainId(domainId) {
    let baseId = domainId;
    const underscoreIndex = baseId.lastIndexOf('_');
    if (underscoreIndex > -1 && /^\d+$/.test(baseId.substring(underscoreIndex + 1))) {
        baseId = baseId.substring(0, underscoreIndex);
    }
    return baseId.replace(/_/g, ' ');
}

function _handleDomainHover(isHovering, instanceId, domainRectId, startLabelId, endLabelId, baseIdLabelId, domainRow, start, end, isIgnored) {
    const domainRect = document.getElementById(domainRectId);
    if (domainRect) domainRect.setAttribute("opacity", isHovering ? "1.0" : "0.6");

    const labelColor = isIgnored ? "red" : "#333";

    const startLabel = document.getElementById(startLabelId);
    const endLabel = document.getElementById(endLabelId);
    if (startLabel) {
        startLabel.setAttribute("visibility", isHovering ? "visible" : "hidden");
        startLabel.setAttribute("fill", labelColor);
    }
    if (endLabel) {
        endLabel.setAttribute("visibility", isHovering ? "visible" : "hidden");
        endLabel.setAttribute("fill", labelColor);
    }

    const baseIdLabel = baseIdLabelId ? document.getElementById(baseIdLabelId) : null;
    if (baseIdLabel) {
        baseIdLabel.setAttribute("visibility", isHovering ? "visible" : "hidden");
        baseIdLabel.setAttribute("fill", labelColor);
    }

    if (domainRow) {
        if (isHovering) {
            domainRow.classList.add('domain-table-row-hover');
        } else {
            domainRow.classList.remove('domain-table-row-hover');
        }
    }
}

function _createSvgElement(tag, attributes = {}, textContent = '') {
    const svgNS = "http://www.w3.org/2000/svg";
    const element = document.createElementNS(svgNS, tag);
    for (const [key, value] of Object.entries(attributes)) {
        element.setAttribute(key, value);
    }
    if (textContent) {
        element.textContent = textContent;
    }
    return element;
}

function _displayInfo(container, message, isWarning = false) {
    const color = isWarning ? 'red' : 'grey';
    if (isWarning) console.error(message);
    if (container) {
        container.innerHTML = `<p style="color:${color}; text-align:center; padding-top: 20px;">${message}</p>`;
    }
}

function _createProteinLabel(value, x, y, { fontSize = 12, textAnchor = "middle", bold = false, angle = 0 } = {}) {
    const label = _createSvgElement("text", {
        "x": x,
        "y": y,
        "dy": "0.35em",
        "text-anchor": textAnchor,
        'dominant-baseline': 'middle',
        "font-size": `${fontSize}px`,
        "font-weight": bold ? "bold" : "normal",
    });
    label.textContent = value.toString();

    if (angle !== 0) {
        label.setAttribute('transform', `rotate(${angle} ${x} ${y})`);
    }

    return label;
}

function _createHoverLabel(value, x, y, { textAnchor = "middle", bold = false } = {}) {
    const label = _createSvgElement("text", {
        "x": x,
        "y": y,
        "dy": "0.35em",
        "text-anchor": textAnchor,
        "font-size": "10px",
        "fill": "#333",
        "visibility": "hidden",
        "font-weight": bold ? "bold" : "normal",
    });
    label.textContent = value.toString();
    return label;
}