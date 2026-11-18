"""Render run records as Markdown tables. Coverage is checked first, so a gap fails loudly (SC-001)."""
from __future__ import annotations

from typing import Iterable, Optional

from .. import records as R

LABELS = {"failed": "failed", "not-run": "not run", "invalid": "invalid"}


def _fmt(value) -> str:
    if isinstance(value, bool):
        return "yes" if value else "no"
    if isinstance(value, float):
        return f"{value:.4g}"
    return str(value)


def _cell_text(cell: dict) -> str:
    kind = cell["kind"]
    if kind == "value":
        return _fmt(cell["value"])
    if kind == "na":
        return f"not applicable: {cell['reason']}"
    text = f"{LABELS[kind]}: {cell['reason']}"
    if kind == "invalid" and cell.get("value") is not None:
        text = f"{_fmt(cell['value'])} ({text})"
    return text


def _escape(text: str) -> str:
    return text.replace("|", "\\|").replace("\n", " ")


def render_markdown(recs: list, tools: Iterable[str], metrics: Iterable[str], *, published: bool = False,
                    exclusions: Optional[list] = None) -> str:
    tools, metrics = list(tools), list(metrics)
    lines = []
    if published:
        flagged = [r for r in recs if r["cache"] != "cold"]
        recs = [r for r in recs if r["cache"] == "cold"]
        for r in flagged:
            lines.append(f"- Excluded from the published table: {r['tool']} on {r['mesh']['name']} ran with "
                         f"cache state '{r['cache']}', not cold.")
        if lines:
            lines.append("")
    R.check_coverage(recs, tools, metrics) if recs or not published else None
    if recs:
        lines.append("| tool | " + " | ".join(metrics) + " |")
        lines.append("|" + "---|" * (len(metrics) + 1))
        for t in tools:
            cells = [_escape(_cell_text(R.cell_for(recs, t, m))) for m in metrics]
            lines.append(f"| {t} | " + " | ".join(cells) + " |")
    else:
        lines.append("No cold-cache records to publish.")
    if exclusions:
        lines += ["", "Excluded meshes:"]
        lines += [f"- {e['mesh']}: {e['reason']}" for e in exclusions]
    return "\n".join(lines) + "\n"


def render_sections(recs: list, tools: Iterable[str], metrics: Iterable[str], **kw) -> str:
    """One table per (workflow, track), each under its own heading."""
    groups: dict = {}
    for r in recs:
        groups.setdefault((r["workflow"], r.get("track")), []).append(r)
    parts = []
    for (wf, track), group in sorted(groups.items(), key=lambda kv: (kv[0][0], kv[0][1] or 0)):
        title = f"{wf} / track {track}" if track else wf
        parts.append(f"## {title}\n\n" + render_markdown(group, tools, metrics, **kw))
    return "\n".join(parts)


def fair_report(directory, *, published=False, rights=None, proposals=None):
    """Recompute a saved v2 campaign offline, checking coverage before redaction.

    Only admitted cells enter aggregation, claim evaluation or returned public
    artifacts. Excluded diagnostics are never interpolated into public text.
    Legacy render_sections stays a separate historical path.
    """
    import json
    import hashlib
    from pathlib import Path
    from ..campaign import validate_manifest, validate_conversion
    from ..aggregate import aggregate_series, METRICS
    from ..claims import evaluate_claim
    from ..valence import fair_rights_decision
    from ..w2.tools import operation_contract
    directory = Path(directory)
    evidence = directory / 'evidence'
    manifest = validate_manifest(json.loads((evidence / 'campaign.json').read_text()))
    recs = R.read_record_hashes(evidence / 'record-hashes.json')
    raw_refs = (evidence / 'references.json').read_bytes()
    reference_hashes = json.loads((evidence / 'reference-hashes.json').read_text())
    if reference_hashes.get('references.json') != hashlib.sha256(raw_refs).hexdigest():
        raise R.RecordError('reference artifact hash mismatch')
    references = json.loads(raw_refs)
    index = json.loads((evidence / 'record-hashes.json').read_text())
    hashes = {r['repeat_id']: r['sha256'] for r in index['records']}
    links = {r['repeat_id']: r['path'] for r in index['records']}
    expected = {e['repeat_id']: e for e in manifest['schedule']}
    if set(expected) != {r['repeat_id'] for r in recs}:
        raise R.MissingCellError('missing or unexpected scheduled attempts')
    groups = {}
    for row in sorted(recs, key=lambda r: r['execution_order']):
        entry = expected[row['repeat_id']]
        contract = operation_contract(entry['tool'], entry['operation'], options=row['options'])
        if (any(row[k] != entry[k] for k in ('cell_id','series_id','repeat_index','execution_order','tool','operation','backend'))
            or row['campaign_id'] != manifest['campaign_id'] or row['manifest_sha256'] != manifest['manifest_sha256']
            or row['contract_id'] != contract['contract_id'] or row['backend'] != contract['backend_selected']
            or row['equivalence'] != contract['equivalence']):
            raise R.RecordError('attempt differs from frozen schedule/operation contract')
        ref = references.get(row['reference_id'])
        if not ref or ref.get('name') != entry['mesh'] or not ref.get('exact_version'):
            raise R.RecordError('missing exact mesh reference binding')
        if row['verification'] == 'verified':
            expected_output = dict(n_nodes=ref.get('header_node_count'), n_elems=ref.get('header_element_count'),
                node_sha256=ref.get('canonical_node_sha256'), element_sha256=ref.get('canonical_element_sha256'))
            if ref.get('verification_state') != 'verified' or any(
                value is None or row['output_verification'].get(k) != value for k,value in expected_output.items()):
                raise R.RecordError('verified output does not bind independent source reference')
        for section, field in [('measurement','samples'),('termination','stdout'),('termination','stderr')]:
            item=row[section]; path=item.get(field+'_path'); digest=item.get(field+'_sha256')
            if path is None and digest is None:
                continue
            if not path or not digest:
                raise R.RecordError('incomplete sidecar binding')
            artifact=Path(path)
            if not artifact.is_absolute():
                local=directory/artifact
                artifact=local if local.exists() else artifact
            if not artifact.resolve().is_relative_to(directory.resolve()):
                raise R.RecordError('sidecar outside saved campaign')
            try:
                payload=artifact.read_bytes()
            except OSError as exc:
                raise R.RecordError('missing evidence sidecar') from exc
            if hashlib.sha256(payload).hexdigest() != digest:
                raise R.RecordError('evidence sidecar hash mismatch')
            if field=='samples' and json.loads(payload) != row['measurement']['samples']:
                raise R.RecordError('sample sidecar differs from retained samples')
        groups.setdefault(row['cell_id'], []).append(row)
    R.check_fair_coverage(recs, [])
    rights = rights or {}
    decisions = {key: fair_rights_decision(ref['name'], ref['exact_version'], rights.get(ref['name']))
                 for key, ref in references.items()}
    excluded_names = sorted({references[r['reference_id']]['name'] for r in recs
                             if not decisions[r['reference_id']]['allowed']}) if published else []
    excluded = [dict(mesh=n, reason='publication rights unavailable or invalid') for n in excluded_names]
    cells, lines = {}, ['# Fair comparison' if published else '# INTERNAL fair comparison', '',
        '10 ms child-only RSS samples can miss shorter spikes. Whole-process RSS includes startup/imports through END. '
        'Operation RSS growth is an upper bound affected by allocator retention and file-backed pages, not total mesh allocation.', '']
    sections = {'Native like-for-like': [], 'conversion / NetCDF': [], 'CHILmesh-only': [], 'not-comparable': []}
    for cell_id, attempts in groups.items():
        first = attempts[0]; ref = references[first['reference_id']]; decision = decisions[first['reference_id']]
        if published and not decision['allowed']:
            continue
        mesh_class = 'stofs' if 'stofs' in ref['name'].lower() else 'small'
        summary = aggregate_series(attempts, mesh_class=mesh_class,
            previous_series_id=manifest.get('previous_series_id'), replacement_rationale=manifest.get('replacement_rationale'))
        contract = operation_contract(first['tool'], first['operation'], options=first['options'])
        conditions = all(all(r['measurement_validity'][k]['valid'] for k in ('time','memory','cache','load','runtime')) for r in attempts)
        conversion = first.get('conversion')
        conversion_valid = all(validate_conversion(r.get('conversion') or {}, ref)['verification_state'] == 'verified' for r in attempts) if first['operation']=='parse-netcdf' else True
        if first['operation']=='parse-netcdf':
            conversion_valid=conversion_valid and all(
                r.get('conversion',{}).get('read_s') == r['measurement']['operation_elapsed_s'] for r in attempts)
        cell = dict(summary, **{k: first[k] for k in ('tool','tool_version','operation','contract_id','backend','equivalence',
            'reference_id','manifest_sha256','host_id','host_allocation','environment_id','thread_environment','options')})
        cell.update(mesh=ref['name'], exact_version=ref['exact_version'], input_format=contract['input_format'],
            input_sha256=first.get('input_sha256') or ref.get('raw_file_sha256'),
            verification='verified' if all(r['verification']=='verified' for r in attempts) and conversion_valid else 'incomplete',
            conditions_valid=conditions and conversion_valid, rights_cleared=decision['allowed'],
            resource_condition_reference=manifest['manifest_sha256'], source_record_hashes=[hashes[r['repeat_id']] for r in attempts],
            conversion=conversion)
        cell['eligible'] = summary['eligible'] and cell['conditions_valid'] and cell['verification']=='verified' and decision['allowed'] and first['equivalence']=='like-for-like'
        cells[cell_id] = cell
        section = ('not-comparable' if first['equivalence']=='not-comparable' else 'conversion / NetCDF'
                   if first['operation']=='parse-netcdf' else 'CHILmesh-only'
                   if first['equivalence']=='CHILmesh-only' else 'Native like-for-like')
        text = [f"### {_escape(ref['name'])} {_escape(ref['exact_version'])} / {first['operation']} / {first['backend']} / "
                f"{contract['input_format']} / series {first['series_id']} / {first['tool']} {first['tool_version']}",
                f"Cell {cell_id}; {'complete' if summary['complete'] else 'incomplete, descriptive only'}; "
                f"{'stable' if summary['stable'] else 'unstable'}; {contract['reason']}."]
        if decision['allowed']:
            text.append(f"Rights: {decision['citation']} ({decision['ruling_date']}), {_escape(decision['scope'])}; evidence sha256 {decision['evidence_sha256']}.")
        for metric in METRICS:
            s = summary['metrics'][metric]
            text.append(f"{metric} ({s['unit']}): raw {s['raw_values']}; n={s['n']}; median={s['median']}; "
                f"min={s['min']}; max={s['max']}; range={s['absolute_spread']}; relative spread={s['relative_spread_pct']}% "
                f"({s['spread_state']}); unstable={s['unstable']}.")
        for row in attempts:
            text.append(f"- [{row['repeat_id']}](evidence/{links[row['repeat_id']]}) sha256 {hashes[row['repeat_id']]}: "
                f"{row['outcome']}: {_escape(row['reason'])}; trigger={row['termination']['trigger']}; "
                f"limits={row['termination']['limits']}; partial measurements={row['measurement'] | {'samples': 'see hashed record'}}; "
                f"diagnostics={_escape(str(row['termination']))}.")
        if conversion:
            text.append('Conversion lineage / separate conversion, read and end-to-end costs: ' + _escape(str(conversion)))
        sections[section].append('\n\n'.join(text))
    for title, content in sections.items():
        lines += [f'## {title}', '', '\n\n'.join(content) or 'No admitted cells.', '']
    claims = []
    for p in proposals or []:
        result = evaluate_claim(p, cells)
        # For excluded or missing supports expose only fixed gate names. Never
        # return user proposal text/numbers or supporting excluded hashes.
        if published and 'cells' in result['failed_predicates']:
            result = dict(claim_id=None, accepted=False, failed_predicates=['cells unavailable for publication'],
                          value=None, sentence=None, source_record_hashes=[])
        claims.append(result)
        lines.append(result['sentence'] if result['accepted'] else 'Claim rejected: ' + ', '.join(result['failed_predicates']))
    if excluded:
        lines += ['', 'Excluded meshes:'] + [f"- {_escape(e['mesh'])}: {e['reason']}" for e in excluded]
    return dict(markdown='\n'.join(lines)+'\n', cells=cells, claims=claims, excluded=excluded,
                publication_passed=(not published or bool(cells)) and all(c['accepted'] for c in claims))
