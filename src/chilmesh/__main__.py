"""Allow ``python -m chilmesh`` to invoke the CLI."""
from __future__ import annotations

from .cli import main

raise SystemExit(main())
