import re
import os
import tempfile
import json
import shutil
import pickle

from flask import Flask, render_template, request, send_file, redirect, url_for
from werkzeug.utils import secure_filename

# Import AlphaFragment modules
from alphafragment.classes import Protein, Domain
from alphafragment.uniprot_fetch import fetch_uniprot_info
from alphafragment.domain_compilation import compile_domains
from alphafragment.fragment_protein import fragment_protein
from alphafragment.process_proteins_csv import update_csv_with_fragments
from alphafragment.fragment_file_creation import output_fastas

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
            
            # Respect current ignore flags from the entry when attaching.
            currently_ignored_flag = bool(d.get('currently_ignored'))
            to_be_ignored_flag = bool(d.get('to_be_ignored'))

            yield GuiDomain(
                identifier=domain_id,
                start=start,
                end=end,
                domain_type=domain_type,
                to_be_ignored=to_be_ignored_flag,
                currently_ignored=currently_ignored_flag,
            )

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
    original_csv_rows = None
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
            original_csv_rows=original_csv_rows,
        )

    # If this is a refragmentation request, we may receive the original CSV
    # as JSON records in a hidden form field so it can be reused client-side.
    csv_rows_raw = request.form.get('original_csv_rows', '')
    if csv_rows_raw.strip():
        try:
            original_csv_rows = json.loads(csv_rows_raw)
        except Exception as e:
            print(f"[DEBUG] Could not parse original_csv_rows: {e}")
            original_csv_rows = None

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
                original_csv_rows=original_csv_rows,
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
                original_csv_rows=original_csv_rows,
            )

        protein_data_list: list[dict] = []
        errors: list[str] = []
        proteins_for_csv: list[Protein] = []

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
                # For fragmentation, use only domains that are not currently ignored,
                # but keep the full domain_list (including ignored ones) on the
                # Protein object so they are preserved for CSV export.
                all_domains = list(getattr(protein, 'domain_list', []) or [])
                active_domains = [
                    d for d in all_domains
                    if not getattr(d, 'currently_ignored', False)
                ]
                setattr(protein, 'domain_list', active_domains)

                fragments = fragment_protein(protein, length=length_params.copy(), overlap=overlap_params.copy())
                for f in fragments:
                    # Keep Protein.fragment_list in sync for CSV export
                    if hasattr(protein, 'add_fragment'):
                        protein.add_fragment(f)
                # Restore full domain list (including currently ignored) for CSV
                setattr(protein, 'domain_list', all_domains)
                fragment_indices = [[int(f[0]) + 1, int(f[1]) + 1] for f in fragments]
                print(f"[DEBUG] {safe_name}: fragmentIndices={fragment_indices}")

                # Collect for optional CSV export
                proteins_for_csv.append(protein)

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

        # If the user requested final output, generate a CSV file
        if action == 'finish_generate':
            import pandas as pd

            if original_csv_rows:
                df = pd.DataFrame(original_csv_rows)
            else:
                # Fallback: construct a minimal DataFrame from the current list
                rows = []
                for entry in current_list:
                    if not isinstance(entry, dict):
                        continue
                    rows.append({
                        'name': entry.get('name') or entry.get('accessionId') or 'protein',
                        'accession_id': entry.get('accessionId') or '',
                    })
                df = pd.DataFrame(rows)

            # Use a deterministic temporary file path so it can be downloaded later
            output_path = os.path.join(tempfile.gettempdir(), 'alphafragment_fragments.csv')

            # Add an extra column that contains only the domains actually used
            # in fragmentation (i.e. those that are not currently ignored).
            try:
                domains_used_dict = {}
                for protein in proteins_for_csv:
                    domain_list = getattr(protein, 'domain_list', None)
                    if domain_list:
                        domain_counts = {}
                        domain_entries = []
                        for domain in domain_list:
                            if getattr(domain, 'currently_ignored', False):
                                # Skip domains that were marked as currently ignored
                                # when computing "used for fragmentation".
                                continue
                            domain_id = domain.id
                            # Handle duplicate domain IDs by appending a counter
                            domain_counts[domain_id] = domain_counts.get(domain_id, 0) + 1
                            count = domain_counts[domain_id]
                            if count > 1:
                                domain_id = f"{domain_id}_{count}"
                            domain_entry = f"('{domain_id}', ({domain.start + 1}, {domain.end + 1}))"
                            domain_entries.append(domain_entry)
                        domains_str = '[' + ', '.join(domain_entries) + ']'
                        domains_used_dict[protein.name] = domains_str
                    else:
                        domains_used_dict[protein.name] = ''

                if 'name' in df.columns:
                    df['domains_used_for_fragmentation'] = df['name'].map(domains_used_dict)
            except Exception as e:
                print(f"[DEBUG] Failed to compute domains_used_for_fragmentation column: {e}")

            try:
                update_csv_with_fragments(df, output_path, proteins_for_csv)
            except Exception as e:
                print(f"[DEBUG] Failed to generate CSV via update_csv_with_fragments: {e}")
                # Fall back to rendering the page with an error
                fragmentation_report = _compute_fragmentation_report(
                    proteins_provided=len(current_list),
                    protein_data_list=protein_data_list,
                )
                return render_template(
                    'index.html',
                    error=f"Failed to generate CSV: {e}",
                    protein_data_list=protein_data_list,
                    proteins=[(p.get('name') or p.get('accessionId') or 'protein') for p in protein_data_list],
                    fragment_params={'length': length_params, 'overlap': overlap_params},
                    fragmentation_report=fragmentation_report,
                    original_csv_rows=original_csv_rows,
                )

            # Persist proteins so FASTA files can be generated later on demand.
            proteins_pickle_path = os.path.join(tempfile.gettempdir(), 'alphafragment_proteins_for_fastas.pkl')
            try:
                with open(proteins_pickle_path, 'wb') as f:
                    pickle.dump(proteins_for_csv, f)
            except Exception as e:
                print(f"[DEBUG] Failed to persist proteins_for_fastas: {e}")

            # CSV has been generated; redirect to a dedicated download page
            return redirect(url_for('fragments_ready'))

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
            original_csv_rows=original_csv_rows,
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
                original_csv_rows=original_csv_rows,
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
                original_csv_rows=original_csv_rows,
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

            # Convert the full CSV DataFrame to a list-of-dicts so it can be
            # serialized to JSON and stored client-side for later use.
            try:
                original_csv_rows = _df.to_dict(orient='records')
            except Exception as e:
                print(f"[DEBUG] Could not convert CSV DataFrame to dict records: {e}")
                original_csv_rows = None
        except Exception as e:
            return render_template(
                'index.html',
                error=f"Error reading CSV: {e}",
                protein_data_list=None,
                proteins=None,
                fragment_params={'length': length_params, 'overlap': overlap_params},
                original_csv_rows=original_csv_rows,
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
                original_csv_rows=original_csv_rows,
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
            original_csv_rows=original_csv_rows,
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
            original_csv_rows=original_csv_rows,
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
        original_csv_rows=original_csv_rows,
    )

@app.route('/fragments_ready', methods=['GET'])
def fragments_ready():
    """Show a page with a link to download the fragments CSV if available."""
    output_path = os.path.join(tempfile.gettempdir(), 'alphafragment_fragments.csv')
    csv_available = os.path.exists(output_path)

    fastas_zip_path = os.path.join(tempfile.gettempdir(), 'alphafragment_fastas.zip')
    fastas_available = os.path.exists(fastas_zip_path)
    return render_template(
        'fragments_ready.html',
        csv_available=csv_available,
        fastas_available=fastas_available,
    )


@app.route('/generate_fastas', methods=['POST'])
def generate_fastas():
    """Generate FASTA files and a zip archive on demand after CSV creation."""
    csv_path = os.path.join(tempfile.gettempdir(), 'alphafragment_fragments.csv')
    fastas_root = os.path.join(tempfile.gettempdir(), 'alphafragment_fastas')
    fastas_zip_path = os.path.join(tempfile.gettempdir(), 'alphafragment_fastas.zip')
    proteins_pickle_path = os.path.join(tempfile.gettempdir(), 'alphafragment_proteins_for_fastas.pkl')

    csv_available = os.path.exists(csv_path)

    if not csv_available:
        # Cannot generate FASTAs without a CSV run first
        return render_template(
            'fragments_ready.html',
            csv_available=False,
            fastas_available=False,
            error="Please generate the fragment CSV before creating FASTA files.",
        )

    # Load persisted proteins
    if not os.path.exists(proteins_pickle_path):
        return render_template(
            'fragments_ready.html',
            csv_available=csv_available,
            fastas_available=False,
            error="Protein data for FASTA generation is not available. Please rerun fragmentation.",
        )

    try:
        with open(proteins_pickle_path, 'rb') as f:
            proteins_for_fastas = pickle.load(f)
    except Exception as e:
        print(f"[DEBUG] Failed to load proteins_for_fastas: {e}")
        return render_template(
            'fragments_ready.html',
            csv_available=csv_available,
            fastas_available=False,
            error="Could not load protein data for FASTA generation.",
        )

    try:
        # Clean previous FASTA outputs
        if os.path.exists(fastas_root):
            shutil.rmtree(fastas_root)
        if os.path.exists(fastas_zip_path):
            os.remove(fastas_zip_path)

        os.makedirs(fastas_root, exist_ok=True)

        # Generate FASTA files
        output_fastas(proteins_for_fastas, save_location=fastas_root, method='all')

        # Zip the folder
        shutil.make_archive(
            base_name=os.path.splitext(fastas_zip_path)[0],
            format='zip',
            root_dir=fastas_root,
        )
    except Exception as e:
        print(f"[DEBUG] Error during FASTA generation: {e}")
        return render_template(
            'fragments_ready.html',
            csv_available=csv_available,
            fastas_available=False,
            error="Failed to generate FASTA files. Check server logs for details.",
        )

    # On success, show page with both CSV and FASTA downloads available
    return render_template(
        'fragments_ready.html',
        csv_available=True,
        fastas_available=True,
    )


@app.route('/download_fragments_csv', methods=['GET'])
def download_fragments_csv():
    """Serve the most recently generated fragments CSV, if available."""
    output_path = os.path.join(tempfile.gettempdir(), 'alphafragment_fragments.csv')
    if not os.path.exists(output_path):
        # No CSV has been generated yet; send back to the ready page with a message
        fastas_zip_path = os.path.join(tempfile.gettempdir(), 'alphafragment_fastas.zip')
        return render_template(
            'fragments_ready.html',
            csv_available=False,
            fastas_available=os.path.exists(fastas_zip_path),
            error="No fragment CSV is available to download yet.",
        )

    return send_file(
        output_path,
        mimetype='text/csv',
        as_attachment=True,
        download_name='alphafragment_fragments.csv',
    )


@app.route('/download_fastas_zip', methods=['GET'])
def download_fastas_zip():
    """Serve a zip archive containing FASTA files for fragment pairs."""
    fastas_zip_path = os.path.join(tempfile.gettempdir(), 'alphafragment_fastas.zip')
    csv_path = os.path.join(tempfile.gettempdir(), 'alphafragment_fragments.csv')

    if not os.path.exists(fastas_zip_path):
        return render_template(
            'fragments_ready.html',
            csv_available=os.path.exists(csv_path),
            fastas_available=False,
            error="No FASTA archive is available to download yet.",
        )

    return send_file(
        fastas_zip_path,
        mimetype='application/zip',
        as_attachment=True,
        download_name='alphafragment_fastas.zip',
    )

if __name__ == '__main__':
    app.run(debug=True, port=5501)
