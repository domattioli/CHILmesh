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
