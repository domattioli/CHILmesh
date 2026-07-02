"""ADCIRC fort.15 (model control / run parameters) file I/O for CHILmesh.

fort.15 is version-dependent and heavily conditional (block ordering depends on
NWS / NTIP / NBFR / ... flags), so it is stored **byte-preserving**: ``read_fort15``
retains the full original text and extracts only the leading, unconditional header
fields; ``write_fort15`` emits the preserved bytes verbatim. This mirrors the
operator directive on issue #201 (fort.15 is a run-config sidecar, kept out of the
``.chil`` Identity hash, and must not be parse-and-rewritten).

Additive module on the ``fort13_io`` / ``gmsh_io`` precedent — no locked stage
module or constitution surface is touched.
"""
from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path


class Fort15ParseError(ValueError):
    """Raised when parsing a fort.15 file encounters an error."""
    pass


def _first_value(line: str) -> str:
    """Return the first whitespace-delimited token before any ``!`` comment."""
    body = line.split("!", 1)[0].strip()
    tokens = body.split()
    if not tokens:
        return ""
    return tokens[0]


@dataclass
class Fort15:
    """Container for an ADCIRC fort.15 file.

    Only the leading, unconditional header fields are extracted. The full
    original file text is retained in ``raw_text`` so ``write_fort15`` can
    reproduce the file byte-for-byte.
    """
    rundes: str          # run description (line 1)
    runid: str           # run identification (line 2)
    nfover: int
    nabout: int
    nscreen: int
    ihot: int
    ics: int
    im: int
    raw_text: str        # exact original file content (lossless write source)


def read_fort15(filename: str | Path) -> Fort15:
    """Read an ADCIRC fort.15 file, byte-preserving.

    Parses only the leading unconditional header (RUNDES, RUNID, NFOVER, NABOUT,
    NSCREEN, IHOT, ICS, IM); everything after IM is version-conditional and is
    retained verbatim in ``Fort15.raw_text`` rather than parsed.

    Parameters:
        filename: Path to the fort.15 file.

    Returns:
        Fort15 object with extracted header fields and the raw text.

    Raises:
        Fort15ParseError: If the file has fewer than 8 header lines or an
            expected integer header field cannot be parsed.
    """
    filename = Path(filename)
    with open(filename, "r", encoding="utf-8") as f:
        raw_text = f.read()

    lines = raw_text.splitlines()
    if len(lines) < 8:
        raise Fort15ParseError(
            f"fort.15 file too short: need at least 8 header lines, got {len(lines)}"
        )

    rundes = lines[0].split("!", 1)[0].rstrip()
    runid = lines[1].split("!", 1)[0].rstrip()

    int_fields = ["NFOVER", "NABOUT", "NSCREEN", "IHOT", "ICS", "IM"]
    values: list[int] = []
    for offset, name in enumerate(int_fields, start=2):
        tok = _first_value(lines[offset])
        try:
            values.append(int(tok))
        except ValueError as e:
            raise Fort15ParseError(
                f"fort.15 header field {name} (line {offset + 1}) parse error: {e}"
            )

    nfover, nabout, nscreen, ihot, ics, im = values

    return Fort15(
        rundes=rundes,
        runid=runid,
        nfover=nfover,
        nabout=nabout,
        nscreen=nscreen,
        ihot=ihot,
        ics=ics,
        im=im,
        raw_text=raw_text,
    )


def write_fort15(f15: Fort15, filename: str | Path) -> None:
    """Write a fort.15 file byte-for-byte from the preserved raw text.

    fort.15 is intentionally not regenerated from the parsed header fields
    (version-specific block ordering is not safely reproducible); the original
    bytes captured at read time are emitted verbatim.

    Parameters:
        f15: Fort15 object to write (must carry ``raw_text``).
        filename: Output path.

    Raises:
        Fort15ParseError: If ``f15.raw_text`` is empty (nothing to preserve).
    """
    if not f15.raw_text:
        raise Fort15ParseError(
            "Fort15.raw_text is empty; fort.15 is byte-preserving and cannot be "
            "regenerated from header fields alone."
        )
    filename = Path(filename)
    with open(filename, "w", encoding="utf-8", newline="") as f:
        f.write(f15.raw_text)


__all__ = ["Fort15", "read_fort15", "write_fort15", "Fort15ParseError"]
