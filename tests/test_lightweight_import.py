"""Regression: `import chilmesh` lightweight fort.14 path stays stdlib-only (#255).

`from chilmesh import read_fort14_raw` must not require numpy or matplotlib, so
downstream pure-stdlib consumers (Valence's fort.14 delegation, #214) keep
working. The stdlib check runs in a subprocess with the heavy deps blocked, so
it is independent of whatever the parent pytest process already imported.
"""
import subprocess
import sys
import textwrap


def test_fort14_raw_import_needs_no_numpy_or_matplotlib():
    code = textwrap.dedent(
        """
        import sys
        # Block the heavy stack: any import attempt now raises ModuleNotFoundError.
        for _m in ("numpy", "matplotlib", "matplotlib.pyplot", "matplotlib.cm", "scipy"):
            sys.modules[_m] = None
        from chilmesh import read_fort14_raw, write_fort14_raw, Fort14Raw
        assert callable(read_fort14_raw)
        assert callable(write_fort14_raw)
        print("LIGHTWEIGHT_OK")
        """
    )
    result = subprocess.run(
        [sys.executable, "-c", code], capture_output=True, text=True
    )
    assert result.returncode == 0, (
        "import chilmesh.read_fort14_raw pulled in a blocked heavy dep:\n"
        + result.stderr
    )
    assert "LIGHTWEIGHT_OK" in result.stdout


def test_heavy_surface_still_available():
    """The lazy names still resolve when the heavy deps ARE present."""
    import chilmesh

    assert chilmesh.Mesh is chilmesh.CHILmesh
    assert callable(chilmesh.summary)
    assert hasattr(chilmesh, "chilplotting")
