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

    if len(lines) < 3:
        raise Fort13ParseError("fort.13 file too short (need at least 3 lines)")

    # Parse header
    grid_name = lines[0]
    try:
        num_nodes = int(lines[1])
        num_attrs = int(lines[2])
    except ValueError as e:
        raise Fort13ParseError(f"fort.13 header parse error: {e}")

    # Parse metadata section
    attributes: list[NodalAttribute] = []
    line_idx = 3
    attr_names_in_order = []

    for _ in range(num_attrs):
        if line_idx + 3 > len(lines):
            raise Fort13ParseError("fort.13 metadata section incomplete")

        attr_name = lines[line_idx]
        units = lines[line_idx + 1]
        try:
            vpn = int(lines[line_idx + 2])
        except ValueError as e:
            raise Fort13ParseError(f"values_per_node parse error at line {line_idx + 2}: {e}")

        line_idx += 3

        # Parse default values
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

        line_idx += 1

        attributes.append(NodalAttribute(
            name=attr_name,
            units=units,
            values_per_node=vpn,
            default_values=default_values,
            nondefault={}
        ))
        attr_names_in_order.append(attr_name)

    # Parse data section
    for _ in range(num_attrs):
        if line_idx >= len(lines):
            raise Fort13ParseError("fort.13 data section incomplete")

        attr_name = lines[line_idx]
        if attr_name not in attr_names_in_order:
            raise Fort13ParseError(f"Unknown attribute '{attr_name}' in data section")

        # Find the attribute object
        attr_obj = next(a for a in attributes if a.name == attr_name)

        line_idx += 1
        if line_idx >= len(lines):
            raise Fort13ParseError(f"fort.13 num_nondefault line missing for '{attr_name}'")

        try:
            num_nondefault = int(lines[line_idx])
        except ValueError as e:
            raise Fort13ParseError(f"num_nondefault parse error: {e}")

        line_idx += 1

        # Parse nondefault rows
        for _ in range(num_nondefault):
            if line_idx >= len(lines):
                raise Fort13ParseError(f"fort.13 data row missing for '{attr_name}'")

            tokens = lines[line_idx].split()
            if len(tokens) != 1 + attr_obj.values_per_node:
                raise Fort13ParseError(
                    f"Data row for '{attr_name}' expects 1 + {attr_obj.values_per_node} tokens, "
                    f"got {len(tokens)}"
                )

            try:
                # Convert 1-based node id to 0-based
                node_id_1based = int(float(tokens[0]))
                node_id = node_id_1based - 1

                if not (0 <= node_id < num_nodes):
                    raise Fort13ParseError(f"Node ID {node_id_1based} out of range [1, {num_nodes}]")

                values = np.array([float(tok) for tok in tokens[1:]], dtype=np.float64)
                attr_obj.nondefault[node_id] = values
            except (ValueError, IndexError) as e:
                raise Fort13ParseError(f"Data row parse error: {e}")

            line_idx += 1

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
