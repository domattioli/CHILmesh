"""ADCIRC fort.13 (nodal attributes) file I/O for CHILmesh.

Supports reading and writing ADCIRC fort.13 nodal attribute files with
round-trip fidelity. Node IDs in fort.13 are 1-based; internally CHILmesh
uses 0-based indexing.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
import numpy as np


class Fort13ParseError(ValueError):
    """Raised when parsing a fort.13 file encounters an error."""
    pass


@dataclass
class NodalAttribute:
    """A single nodal attribute (e.g. manning roughness, elevation)."""
    name: str
    units: str
    values_per_node: int
    default_values: np.ndarray  # shape (values_per_node,), dtype float64
    nondefault: dict[int, np.ndarray] = field(default_factory=dict)  # 0-based node_id -> array


@dataclass
class Fort13:
    """Container for fort.13 nodal attributes."""
    grid_name: str
    num_nodes: int
    attributes: list[NodalAttribute]

    def attribute(self, name: str) -> NodalAttribute:
        """Retrieve attribute by name; raise KeyError if not found."""
        for attr in self.attributes:
            if attr.name == name:
                return attr
        raise KeyError(f"Attribute '{name}' not found in fort.13")

    def dense(self, name: str) -> np.ndarray:
        """Return dense array of attribute values, shape (num_nodes, values_per_node).

        Fills with default values, overlaid with nondefault entries.
        """
        attr = self.attribute(name)
        arr = np.tile(attr.default_values, (self.num_nodes, 1))
        for node_id, values in attr.nondefault.items():
            arr[node_id] = values
        return arr


def _parse_f13_header(lines: list) -> tuple:
    """Return ``(grid_name, num_nodes, num_attrs)`` from the first 3 lines."""
    if len(lines) < 3:
        raise Fort13ParseError("fort.13 file too short (need at least 3 lines)")
    grid_name = lines[0]
    try:
        num_nodes = int(lines[1])
        num_attrs = int(lines[2])
    except ValueError as e:
        raise Fort13ParseError(f"fort.13 header parse error: {e}")
    return grid_name, num_nodes, num_attrs


def _parse_f13_attribute_meta(lines: list, line_idx: int) -> tuple:
    """Parse one metadata entry (name, units, count, defaults).

    Returns
    -------
    tuple
        ``(NodalAttribute, next_line_index)``.
    """
    if line_idx + 3 > len(lines):
        raise Fort13ParseError("fort.13 metadata section incomplete")

    attr_name = lines[line_idx]
    units = lines[line_idx + 1]
    try:
        vpn = int(lines[line_idx + 2])
    except ValueError as e:
        raise Fort13ParseError(f"values_per_node parse error at line {line_idx + 2}: {e}")

    line_idx += 3

    if line_idx >= len(lines):
        raise Fort13ParseError("fort.13 default values line missing")

    default_tokens = lines[line_idx].split()
    if len(default_tokens) != vpn:
        raise Fort13ParseError(
            f"Attribute '{attr_name}' expects {vpn} default values, got {len(default_tokens)}"
        )

    try:
        default_values = np.array([float(tok) for tok in default_tokens], dtype=np.float64)
    except ValueError as e:
        raise Fort13ParseError(f"Default values parse error: {e}")

    attr = NodalAttribute(
        name=attr_name,
        units=units,
        values_per_node=vpn,
        default_values=default_values,
        nondefault={}
    )
    return attr, line_idx + 1


def _parse_f13_data_row(line: str, attr: NodalAttribute, num_nodes: int) -> tuple:
    """Parse one nondefault row into ``(0-based node id, values)``."""
    tokens = line.split()
    if len(tokens) != 1 + attr.values_per_node:
        raise Fort13ParseError(
            f"Data row for '{attr.name}' expects 1 + {attr.values_per_node} tokens, "
            f"got {len(tokens)}"
        )

    try:
        # Convert 1-based node id to 0-based
        node_id_1based = int(float(tokens[0]))
        node_id = node_id_1based - 1

        if not (0 <= node_id < num_nodes):
            raise Fort13ParseError(f"Node ID {node_id_1based} out of range [1, {num_nodes}]")

        values = np.array([float(tok) for tok in tokens[1:]], dtype=np.float64)
    except (ValueError, IndexError) as e:
        # Fort13ParseError is a ValueError, so the range error is re-wrapped
        # here; existing callers see the "Data row parse error" prefix.
        raise Fort13ParseError(f"Data row parse error: {e}")
    return node_id, values


def _parse_f13_data_block(
    lines: list, line_idx: int, attributes: list, num_nodes: int
) -> int:
    """Parse one attribute data block into the matching attribute.

    Returns the index of the first line after the block.
    """
    if line_idx >= len(lines):
        raise Fort13ParseError("fort.13 data section incomplete")

    attr_name = lines[line_idx]
    if attr_name not in [a.name for a in attributes]:
        raise Fort13ParseError(f"Unknown attribute '{attr_name}' in data section")

    attr_obj = next(a for a in attributes if a.name == attr_name)

    line_idx += 1
    if line_idx >= len(lines):
        raise Fort13ParseError(f"fort.13 num_nondefault line missing for '{attr_name}'")

    try:
        num_nondefault = int(lines[line_idx])
    except ValueError as e:
        raise Fort13ParseError(f"num_nondefault parse error: {e}")

    line_idx += 1

    for _ in range(num_nondefault):
        if line_idx >= len(lines):
            raise Fort13ParseError(f"fort.13 data row missing for '{attr_name}'")
        node_id, values = _parse_f13_data_row(lines[line_idx], attr_obj, num_nodes)
        attr_obj.nondefault[node_id] = values
        line_idx += 1
    return line_idx


def read_fort13(filename: str | Path) -> Fort13:
    """Read a fort.13 nodal attribute file.

    Converts 1-based node IDs to 0-based internal indexing.

    Parameters:
        filename: Path to the .13 file

    Returns:
        Fort13 object with parsed attributes

    Raises:
        Fort13ParseError: If file is malformed
    """
    filename = Path(filename)
    with open(filename, 'r', encoding='utf-8') as f:
        lines = [line.strip() for line in f]

    # Skip blank lines
    lines = [line for line in lines if line]

    grid_name, num_nodes, num_attrs = _parse_f13_header(lines)

    attributes: list[NodalAttribute] = []
    line_idx = 3
    for _ in range(num_attrs):
        attr, line_idx = _parse_f13_attribute_meta(lines, line_idx)
        attributes.append(attr)

    for _ in range(num_attrs):
        line_idx = _parse_f13_data_block(lines, line_idx, attributes, num_nodes)

    return Fort13(
        grid_name=grid_name,
        num_nodes=num_nodes,
        attributes=attributes
    )


def write_fort13(f13: Fort13, filename: str | Path) -> None:
    """Write a fort.13 nodal attribute file.

    Converts 0-based node IDs to 1-based for output (ADCIRC convention).

    Parameters:
        f13: Fort13 object to write
        filename: Output path
    """
    filename = Path(filename)
    with open(filename, 'w', encoding='utf-8') as f:
        # Write header
        f.write(f"{f13.grid_name}\n")
        f.write(f"{f13.num_nodes}\n")
        f.write(f"{len(f13.attributes)}\n")

        # Write metadata section
        for attr in f13.attributes:
            f.write(f"{attr.name}\n")
            f.write(f"{attr.units}\n")
            f.write(f"{attr.values_per_node}\n")
            # Write default values
            default_str = " ".join(repr(float(v)) for v in attr.default_values)
            f.write(f"{default_str}\n")

        # Write data section
        for attr in f13.attributes:
            f.write(f"{attr.name}\n")
            f.write(f"{len(attr.nondefault)}\n")
            # Sort by node id for consistent output
            for node_id in sorted(attr.nondefault.keys()):
                node_id_1based = node_id + 1
                values = attr.nondefault[node_id]
                values_str = " ".join(repr(float(v)) for v in values)
                f.write(f"{node_id_1based} {values_str}\n")


__all__ = ["Fort13", "NodalAttribute", "read_fort13", "write_fort13", "Fort13ParseError"]
