import re
import os
import tempfile
import json

from flask import Flask, render_template, request
from werkzeug.utils import secure_filename

# Import AlphaFragment modules
from alphafragment.classes import Protein, Domain
from alphafragment.uniprot_fetch import fetch_uniprot_info
from alphafragment.domain_compilation import compile_domains
from alphafragment.fragment_protein import fragment_protein

app = Flask(__name__)
app.config['MAX_CONTENT_LENGTH'] = 5 * 1024 * 1024  # 5 MB


DEFAULT_FRAGMENT_PARAMS = {
    'length': {'min': 200, 'ideal': 384, 'max': 400},
    'overlap': {'min': 20, 'ideal': 30, 'max': 40},
}


def _compute_fragmentation_report(
    *,
    proteins_provided: int,
    protein_data_list: list[dict],
) -> dict:
    provided = int(proteins_provided or 0)
    ok = len(protein_data_list or [])
    errors = max(0, provided - ok)
    return {
        'proteins_provided': provided,
        'errors': errors,
    }


def _parse_fragmentation_params(form) -> tuple[dict, dict]:
    def _parse_int_field(field_name: str, *, default: int, min_value: int) -> int:
        raw = form.get(field_name, None)
        if raw is None or str(raw).strip() == '':
            return default
        try:
            value = int(str(raw).strip())
        except (TypeError, ValueError):
            raise ValueError(f"{field_name} must be an integer.")
        if value < min_value:
            raise ValueError(f"{field_name} must be >= {min_value}.")
        return value

    length = {
        'min': _parse_int_field('fragment_len_min', default=DEFAULT_FRAGMENT_PARAMS['length']['min'], min_value=1),
        'ideal': _parse_int_field('fragment_len_ideal', default=DEFAULT_FRAGMENT_PARAMS['length']['ideal'], min_value=1),
        'max': _parse_int_field('fragment_len_max', default=DEFAULT_FRAGMENT_PARAMS['length']['max'], min_value=1),
    }
    overlap = {
        'min': _parse_int_field('overlap_len_min', default=DEFAULT_FRAGMENT_PARAMS['overlap']['min'], min_value=0),
        'ideal': _parse_int_field('overlap_len_ideal', default=DEFAULT_FRAGMENT_PARAMS['overlap']['ideal'], min_value=0),
        'max': _parse_int_field('overlap_len_max', default=DEFAULT_FRAGMENT_PARAMS['overlap']['max'], min_value=0),
    }

    if length['min'] > length['max']:
        raise ValueError(f"Minimum fragment length ({length['min']}) must be less than maximum fragment length ({length['max']}).")
    if not (length['min'] <= length['ideal'] <= length['max']):
        raise ValueError("Ideal length must be within the min and max length bounds.")
    if overlap['min'] > overlap['max']:
        raise ValueError(f"Minimum overlap ({overlap['min']}) must be less than or equal to maximum overlap ({overlap['max']}).")
    if not (overlap['min'] <= overlap['ideal'] <= overlap['max']):
        raise ValueError("Ideal overlap must be within the min and max overlap bounds.")
    if overlap['max'] >= length['min']:
        raise ValueError(
            f"Maximum overlap ({overlap['max']}) must be less than the minimum fragment length ({length['min']}) to avoid overlap-length conflicts."
        )

    return length, overlap


def parse_uniprot_ids(raw_value: str) -> list[str]:
    """Parse a comma/whitespace separated list of UniProt IDs.

    Keeps input order and removes duplicates (case-insensitive).
    """
    raw_value = (raw_value or '').strip()
    if not raw_value:
        return []

    tokens = [t for t in re.split(r'[\s,;]+', raw_value) if t]
    ids: list[str] = []
    seen: set[str] = set()
    for token in tokens:
        cleaned = token.strip().upper()
        if not cleaned:
            continue
        # Be permissive: UniProt accessions can include isoform suffixes like P12345-2.
        if not re.fullmatch(r'[A-Z0-9][A-Z0-9-]*', cleaned):
            continue
        if cleaned in seen:
            continue
        seen.add(cleaned)
        ids.append(cleaned)
    return ids


def _is_alphafold_domain(type_value, domain_id) -> bool:
    """Return True if a domain should render as AlphaFold (light blue)."""
    if type_value is not None:
        type_str = str(type_value).strip().lower()
        # Upstream may label AlphaFold as 'alphafold' instead of the frontend's 'af'.
        return type_str in {'af', 'alphafold'}

    # Preserve previous heuristic when upstream doesn't provide a type.
    domain_id_str = str(domain_id).lower() if domain_id is not None else ''
    return 'uniprot' not in domain_id_str


# Helper to convert Protein object to dict for frontend JS
def protein_to_json(protein):
    """
    Convert Protein object to dict for frontend JS.
    Includes explicit serialization of ignore states (currently_ignored/to_be_ignored).
    """
    # Domains: split by type
    alphafold_domains = []
    uniprot_domains = []
    
    # Iterate all domains attached to the protein object
    for d in protein.domain_list:
        domain_id = getattr(d, 'id', None)
        raw_type = getattr(d, 'type', None)
        
        # Heuristics for typing if not explicit
        is_alphafold = _is_alphafold_domain(raw_type, domain_id)
        d_type = raw_type if raw_type is not None else ('af' if is_alphafold else 'uniprot')
        
        # Serialize fields, checking for our extended attributes
        currently_ignored = getattr(d, 'currently_ignored', False)
        to_be_ignored = getattr(d, 'to_be_ignored', False)
        
        d_dict = {
            'id': domain_id,
            'start': int(d.start) + 1,  # 1-based for JS
            'end': int(d.end) + 1,
            'type': d_type,
            'currently_ignored': currently_ignored,
            'to_be_ignored': to_be_ignored
        }
        
        if is_alphafold:
            alphafold_domains.append(d_dict)
        else:
            uniprot_domains.append(d_dict)
            
    # Fragments: list of [start, end], 1-based
    fragments = [[int(f[0]) + 1, int(f[1]) + 1] for f in protein.fragment_list]
    seq = getattr(protein, 'sequence', None)
    if isinstance(seq, str) and seq:
        length = len(seq)
    else:
        last_res = getattr(protein, 'last_res', None)
        length = (int(last_res) + 1) if last_res is not None else None
    
    is_approved = getattr(protein, 'is_approved', False)

    return {
        'name': protein.name,
        'accessionId': protein.accession_id,
        'length': length,
        'fragmentIndices': fragments,
        'alphafoldDomains': alphafold_domains,
        'uniprotDomains': uniprot_domains,
        'sequence': seq,
        'isApproved': is_approved
    }


class GuiDomain(Domain):
    """
    Extends the alphafragment Domain class to include GUI-specific ignore states.
    
    Attributes:
        - to_be_ignored (bool): Marked for ignoring in the next fragmentation pass (visual grey).
        - currently_ignored (bool): Already ignored in a previous pass (hidden/inactive).
    """
    def __init__(self, identifier, start, end, domain_type, to_be_ignored=False, currently_ignored=False):
        super().__init__(identifier, start, end, domain_type)
        self.to_be_ignored = to_be_ignored
        self.currently_ignored = currently_ignored
    
    def __str__(self):
        base = super().__str__()
        flags = []
        if self.to_be_ignored: flags.append("to_be_ignored")
        if self.currently_ignored: flags.append("currently_ignored")
        return f"{base} [{' '.join(flags)}]" if flags else base

    def __repr__(self):
        return f"GuiDomain(id={self.id}, start={self.start}, end={self.end}, type='{self.type}', to_be_ignored={self.to_be_ignored}, currently_ignored={self.currently_ignored})"



def _attach_domains_from_entry(protein: Protein, entry: dict) -> None:
    """Attach domains already present in current_protein_data to the Protein.

    This keeps refragment actions self-contained: no UniProt/AlphaFold refetch.
    Domain coordinates in the frontend are 1-based inclusive; Protein expects 0-based.
    
    WARNING: Modifies entry['uniprotDomains'] and entry['alphafoldDomains'] in place
    to update ignore statuses.
    """

    def _process_domains(list_value, *, fallback_type: str):
        if not isinstance(list_value, list):
            return
        
        # We iterate a copy because we might modify the dictionaries in the list
        for d in list_value:
            if not isinstance(d, dict):
                continue
            
            # --- Status transitions ---
            # 1. to_be_ignored -> currently_ignored
            if d.get('to_be_ignored') or d.get('toBeIgnored'):
                 d['currently_ignored'] = True
                 d['to_be_ignored'] = False
                 # clean up potential JS keys if we want consistency
                 if 'toBeIgnored' in d: del d['toBeIgnored']
            
            # 2. currently_ignored -> stays currently_ignored (unless unignored separately, handled by absence of flag?)
            # The JS should have set to_be_ignored=False if it was unignored, but currently_ignored might persist?
            # User requirement: "'currently ignored' is when it was 'to be ignored' and has now been refragmented"
            # If the user UN-ignores a currently_ignored domain, JS should ideally clear currently_ignored?
            # For now, let's assume 'currently_ignored' means it's OUT of the fragmentation pool.
            
            if d.get('currently_ignored'):
                # This domain is ignored. Do not yield it for the Protein object.
                continue

            raw_start = d.get('start', None)
            raw_end = d.get('end', None)
            try:
                start = int(raw_start) - 1
                end = int(raw_end) - 1
            except Exception:
                continue
            if start < 0 or end < start:
                continue
            domain_id = d.get('id', None)
            domain_type = d.get('type', None) or fallback_type
            
            # Only yield active domains
            yield GuiDomain(identifier=domain_id, start=start, end=end, domain_type=domain_type, to_be_ignored=False, currently_ignored=False)

    domains = []
    domains.extend(list(_process_domains(entry.get('alphafoldDomains'), fallback_type='af')))
    domains.extend(list(_process_domains(entry.get('uniprotDomains'), fallback_type='uniprot')))

    if not domains:
        return

    # Prefer writing directly to domain_list.
    # Some AlphaFragment versions enforce `Domain` instances in `add_domain`,
    # but the fragmentation logic only needs start/end/id/type-like attributes.
    domain_list = getattr(protein, 'domain_list', None)
    if isinstance(domain_list, list):
        domain_list.extend(domains)
        return

    # Fallback: try add_domain (may reject non-Domain instances).
    add_domain = getattr(protein, 'add_domain', None)
    if callable(add_domain):
        for d in domains:
            try:
                add_domain(d)
            except Exception:
                # If add_domain is strict, ignore and proceed without domains.
                pass
        return


@app.route('/', methods=['GET'])
def index():
    return render_template(
        'index.html',
        protein_data_list=None,
        proteins=None,
        fragment_params=DEFAULT_FRAGMENT_PARAMS,
    )


@app.route('/fragment', methods=['POST'])
def fragment():
    try:
        length_params, overlap_params = _parse_fragmentation_params(request.form)
    except Exception as e:
        return render_template(
            'index.html',
            error=str(e),
            protein_data_list=None,
            proteins=None,
            fragment_params={
                'length': {
                    'min': request.form.get('fragment_len_min', DEFAULT_FRAGMENT_PARAMS['length']['min']),
                    'ideal': request.form.get('fragment_len_ideal', DEFAULT_FRAGMENT_PARAMS['length']['ideal']),
                    'max': request.form.get('fragment_len_max', DEFAULT_FRAGMENT_PARAMS['length']['max']),
                },
                'overlap': {
                    'min': request.form.get('overlap_len_min', DEFAULT_FRAGMENT_PARAMS['overlap']['min']),
                    'ideal': request.form.get('overlap_len_ideal', DEFAULT_FRAGMENT_PARAMS['overlap']['ideal']),
                    'max': request.form.get('overlap_len_max', DEFAULT_FRAGMENT_PARAMS['overlap']['max']),
                },
            },
        )

    action = (request.form.get('action') or 'new').strip().lower()
    print(f"[DEBUG] Parsed action: {action}")


    print("[DEBUG] Received form fields:")
    for k in request.form:
        print(f"[DEBUG]   {k}: {request.form.get(k)[:200] if isinstance(request.form.get(k), str) else request.form.get(k)}")

    if action in {'refragment', 'refragment_all', 'refragment_unapproved', 'finish_generate'}:
        payload_raw = request.form.get('current_protein_data', '')
        print("[DEBUG] Received current_protein_data (raw):", payload_raw[:500])
        if not payload_raw.strip():
            print("[DEBUG] No current_protein_data received!")
            return render_template(
                'index.html',
                error="Refragment failed: no current proteins found on the page.",
                protein_data_list=None,
                proteins=None,
                fragment_params={'length': length_params, 'overlap': overlap_params},
            )

        try:
            current_list = json.loads(payload_raw)
            print(f"[DEBUG] Parsed current_protein_data: {type(current_list)} with {len(current_list) if isinstance(current_list, list) else 'N/A'} entries")
            if not isinstance(current_list, list):
                raise ValueError("current_protein_data must be a list")
        except Exception as e:
            print(f"[DEBUG] Error parsing current_protein_data: {e}")
            return render_template(
                'index.html',
                error=f"Refragment failed: could not parse current proteins ({e}).",
                protein_data_list=None,
                proteins=None,
                fragment_params={'length': length_params, 'overlap': overlap_params},
            )

        if action == 'finish_generate':
            # Placeholder until output-file generation is implemented.
            return render_template(
                'index.html',
                error="Finish and generate output files is not implemented yet.",
                protein_data_list=current_list,
                proteins=[(p.get('name') or p.get('accessionId') or 'protein') if isinstance(p, dict) else 'protein' for p in current_list],
                fragment_params={'length': length_params, 'overlap': overlap_params},
                fragmentation_report=_compute_fragmentation_report(
                    proteins_provided=len(current_list),
                    protein_data_list=[p for p in current_list if isinstance(p, dict)],
                ),
            )

        protein_data_list: list[dict] = []
        errors: list[str] = []

        # No longer using separate approved_indices list; isApproved is now a property of the protein object.

        for idx, entry in enumerate(current_list):
            try:
                if not isinstance(entry, dict):
                    raise ValueError("Invalid protein entry")

                # For 'refragment_unapproved', pass approved proteins through untouched.
                is_approved = entry.get('isApproved', False)
                if action == 'refragment_unapproved' and is_approved:
                    # Mark as approved so frontend can restore state (should already be True in entry)
                    entry['isApproved'] = True
                    protein_data_list.append(entry)
                    continue

                name = entry.get('name')
                accession = entry.get('accessionId')
                sequence = entry.get('sequence')
                if not isinstance(sequence, str) or not sequence:
                    raise ValueError("Missing sequence")

                safe_name = (name or accession or 'protein').strip() if isinstance(name or accession or 'protein', str) else 'protein'

                protein = Protein(name=safe_name, accession_id=accession or '', sequence=sequence)
                _attach_domains_from_entry(protein, entry)
                print(f"[DEBUG] {safe_name}: alphafoldDomains={len(entry.get('alphafoldDomains', []))}, uniprotDomains={len(entry.get('uniprotDomains', []))}")
                fragments = fragment_protein(protein, length=length_params.copy(), overlap=overlap_params.copy())
                fragment_indices = [[int(f[0]) + 1, int(f[1]) + 1] for f in fragments]
                print(f"[DEBUG] {safe_name}: fragmentIndices={fragment_indices}")

                updated = {
                    **entry,
                    'isApproved': False,  # Reset approval if re-fragmenting
                    'name': safe_name,
                    'length': len(sequence),
                    'fragmentIndices': fragment_indices,
                }
                protein_data_list.append(updated)
            except Exception as e:
                label = None
                if isinstance(entry, dict):
                    label = entry.get('name') or entry.get('accessionId')
                errors.append(f"{label or 'protein'}: {e}")

        fragmentation_report = _compute_fragmentation_report(
            proteins_provided=len(current_list),
            protein_data_list=protein_data_list,
        )

        print(f"[DEBUG] Rendering template with protein_data_list (len={len(protein_data_list)}): {protein_data_list}")
        return render_template(
            'index.html',
            error="; ".join(errors) if errors else None,
            protein_data_list=protein_data_list,
            proteins=[(p.get('name') or p.get('accessionId') or 'protein') for p in protein_data_list],
            fragment_params={'length': length_params, 'overlap': overlap_params},
            fragmentation_report=fragmentation_report,
        )

    uploaded_csv = request.files.get('protein_csv')
    has_csv = bool(uploaded_csv and uploaded_csv.filename and uploaded_csv.filename.strip())

    proteins_from_csv: list[Protein] = []
    proteins_provided = 0
    init_failures = 0
    if has_csv:
        filename = secure_filename(uploaded_csv.filename)
        if not filename.lower().endswith('.csv'):
            return render_template(
                'index.html',
                error="Please upload a .csv file.",
                protein_data_list=None,
                proteins=None,
                fragment_params={'length': length_params, 'overlap': overlap_params},
            )

        try:
            from alphafragment.process_proteins_csv import initialize_proteins_from_csv
        except Exception as e:
            return render_template(
                'index.html',
                error=f"CSV processing is unavailable (cannot import alphafragment.process_proteins_csv): {e}",
                protein_data_list=None,
                proteins=None,
                fragment_params={'length': length_params, 'overlap': overlap_params},
            )

        tmp_path = None
        try:
            with tempfile.NamedTemporaryFile(mode='wb', suffix='.csv', delete=False) as tmp:
                tmp_path = tmp.name
                uploaded_csv.save(tmp_path)

            proteins_from_csv, _df = initialize_proteins_from_csv(tmp_path)
            try:
                proteins_provided = int(len(_df))
            except Exception:
                proteins_provided = 0
            init_failures = max(0, proteins_provided - len(proteins_from_csv))
        except Exception as e:
            return render_template(
                'index.html',
                error=f"Error reading CSV: {e}",
                protein_data_list=None,
                proteins=None,
                fragment_params={'length': length_params, 'overlap': overlap_params},
            )
        finally:
            if tmp_path and os.path.exists(tmp_path):
                try:
                    os.remove(tmp_path)
                except OSError:
                    pass

        if not proteins_from_csv:
            return render_template(
                'index.html',
                error="No valid proteins could be initialized from the CSV.",
                protein_data_list=None,
                proteins=None,
                fragment_params={'length': length_params, 'overlap': overlap_params},
            )

    raw_uniprot_ids = request.form.get('uniprot_id', '')
    uniprot_ids = parse_uniprot_ids(raw_uniprot_ids)

    if not has_csv and not uniprot_ids:
        return render_template(
            'index.html',
            error="Please enter one or more UniProt IDs (comma-separated) or upload a CSV.",
            protein_data_list=None,
            proteins=None,
            fragment_params={'length': length_params, 'overlap': overlap_params},
        )

    protein_data_list = []
    errors: list[str] = []
    uniprot_issues = 0

    proteins_to_process: list[Protein] = []
    if has_csv:
        proteins_to_process.extend(proteins_from_csv)
    else:
        proteins_provided = len(uniprot_ids)
        for uniprot_id in uniprot_ids:
            try:
                data = fetch_uniprot_info(uniprot_id)
                sequence = data.get('sequence', '')
                if not sequence:
                    raise ValueError("No sequence found for this UniProt ID.")
                proteins_to_process.append(Protein(name=uniprot_id, accession_id=uniprot_id, sequence=sequence))
            except Exception as e:
                errors.append(f"{uniprot_id}: {str(e)}")
                uniprot_issues += 1

    for protein in proteins_to_process:
        try:
            # Only attempt domain compilation when an accession ID exists.
            if getattr(protein, 'accession_id', None):
                domains = compile_domains(protein, uniprot=True, alphafold=True, manual=False)
                for d in domains:
                    protein.add_domain(d)

            fragments = fragment_protein(protein, length=length_params.copy(), overlap=overlap_params.copy())
            for f in fragments:
                protein.add_fragment(f)

            protein_data_list.append(protein_to_json(protein))
        except Exception as e:
            label = getattr(protein, 'name', None) or getattr(protein, 'accession_id', None) or 'protein'
            errors.append(f"{label}: {str(e)}")

    if not protein_data_list:
        return render_template(
            'index.html',
            error="Error: " + ("; ".join(errors) if errors else "No valid UniProt IDs."),
            protein_data_list=None,
            proteins=None,
            fragment_params={'length': length_params, 'overlap': overlap_params},
        )

    fragmentation_report = _compute_fragmentation_report(
        proteins_provided=proteins_provided,
        protein_data_list=protein_data_list,
    )

    return render_template(
        'index.html',
        error="; ".join(errors) if errors else None,
        protein_data_list=protein_data_list,
        proteins=[p['name'] for p in protein_data_list],
        fragment_params={'length': length_params, 'overlap': overlap_params},
        fragmentation_report=fragmentation_report,
    )

if __name__ == '__main__':
    app.run(debug=True, port=5501)
