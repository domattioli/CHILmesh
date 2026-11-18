"""Private progress line for long CHILmesh operations.

One status line goes to stderr, redrawn in place with ``\\r`` and cleared when
the operation ends. stdout is never touched. When ``tqdm`` is importable it
draws the live bar; otherwise a small stdlib renderer is used.

Mode, highest priority first:

1. The private ``_progress`` keyword of a wired method (``True``, ``False``,
   ``None`` or one of ``"0"``, ``"1"``, ``"auto"``).
2. The ``CHILMESH_PROGRESS`` environment variable: ``0``, ``1`` or ``auto``.
3. ``auto``: draw only when stderr is a TTY and the operation has already run
   longer than ``_DELAY`` seconds.

Usage inside a wired method (tick per layer, iteration or chunk, never per
element)::

    with track("Peeling layers", "layer", past="Peeled", _progress=_progress) as t:
        for ...:
            t.tick()
"""
from __future__ import annotations

import os
import sys
import threading
import time
from typing import Any, Optional

# Seconds an operation must run before auto mode draws its first frame.
_DELAY = 1.0
# Seconds between glyph frames.
_FRAME_INTERVAL = 0.12
_GLYPHS = "✢✳✶✻✽"
_ENV = "CHILMESH_PROGRESS"


def _is_tty() -> bool:
    stream = sys.stderr
    try:
        return bool(stream is not None and stream.isatty())
    except Exception:
        return False


def _load_tqdm():
    """Return the ``tqdm`` class, or ``None`` when it is not installed."""
    try:
        from tqdm import tqdm
    except ImportError:
        return None
    return tqdm


def _resolve_mode(override: Any) -> str:
    """Return ``"0"``, ``"1"`` or ``"auto"``."""
    if override is True:
        return "1"
    if override is False:
        return "0"
    if isinstance(override, str) and override.strip().lower() in ("0", "1", "auto"):
        return override.strip().lower()
    env = os.environ.get(_ENV, "auto").strip().lower()
    return env if env in ("0", "1", "auto") else "auto"


def format_elapsed(seconds: float) -> str:
    """Format a duration as ``1.4s``, ``12s`` or ``2m 46s``."""
    if seconds < 10:
        return f"{seconds:.1f}s"
    total = int(seconds)
    if total < 60:
        return f"{total}s"
    minutes, secs = divmod(total, 60)
    if minutes < 60:
        return f"{minutes}m {secs:02d}s"
    hours, minutes = divmod(minutes, 60)
    return f"{hours}h {minutes:02d}m {secs:02d}s"


class _NullTracker:
    """No-op tracker used when progress is off."""

    __slots__ = ()

    def tick(self, n: int = 1) -> None:
        pass

    def __enter__(self) -> "_NullTracker":
        return self

    def __exit__(self, *exc: Any) -> None:
        return None


class _Tracker:
    def __init__(self, verb: str, unit: str, total: Optional[int], past: str, forced: bool):
        self.verb = verb
        self.unit = unit
        self.total = total
        self.past = past
        self.forced = forced
        self.count = 0
        self._t0 = 0.0
        self._lock = threading.Lock()
        self._stop = threading.Event()
        self._thread: Optional[threading.Thread] = None
        self._drawn = False
        self._bar = None
        self._frame = 0

    # -- public ---------------------------------------------------------
    def tick(self, n: int = 1) -> None:
        """Advance the counter by ``n``. O(1)."""
        with self._lock:
            self.count += n
            if self._bar is not None:
                self._bar.update(n)

    def __enter__(self) -> "_Tracker":
        self._t0 = time.monotonic()
        if self.forced:
            self._draw()
        self._thread = threading.Thread(target=self._run, name="chilmesh-progress", daemon=True)
        self._thread.start()
        return self

    def __exit__(self, exc_type, exc, tb) -> None:
        self._stop.set()
        if self._thread is not None:
            self._thread.join()
        with self._lock:
            self._clear()
            if exc_type is None and (self._drawn or self.forced):
                elapsed = format_elapsed(time.monotonic() - self._t0)
                noun = self.unit if self.count == 1 else self.unit + "s"
                self._write(f"✓ Chillin'. {self.past} {self.count} {noun} in {elapsed}\n")
        return None

    # -- internals ------------------------------------------------------
    def _run(self) -> None:
        while not self._stop.wait(_FRAME_INTERVAL):
            if not self._drawn and time.monotonic() - self._t0 <= _DELAY and not self.forced:
                continue
            with self._lock:
                if self._stop.is_set():
                    return
                self._draw_locked()

    def _draw(self) -> None:
        with self._lock:
            self._draw_locked()

    def _line(self) -> str:
        glyph = _GLYPHS[self._frame % len(_GLYPHS)]
        self._frame += 1
        parts = [format_elapsed(time.monotonic() - self._t0)]
        if self.total:
            pct = min(100, int(100 * self.count / self.total))
            parts.append(f"{self.unit} {self.count}/{self.total} · {pct}%")
        else:
            parts.append(f"{self.unit} {self.count}")
        return f"{glyph} {self.verb}… ({' · '.join(parts)})"

    def _draw_locked(self) -> None:
        if self._bar is None and not self._drawn:
            tqdm = _load_tqdm()
            if tqdm is not None:
                self._bar = tqdm(
                    total=self.total, initial=self.count, desc=self.verb,
                    unit=self.unit, file=sys.stderr, leave=False,
                )
                self._drawn = True
                return
        if self._bar is not None:
            return
        self._drawn = True
        self._write("\r" + self._line())

    def _clear(self) -> None:
        if self._bar is not None:
            self._bar.close()
            self._bar = None
        elif self._drawn:
            self._write("\r" + " " * 79 + "\r")

    @staticmethod
    def _write(text: str) -> None:
        stream = sys.stderr
        if stream is None:
            return
        try:
            stream.write(text)
            stream.flush()
        except Exception:
            pass


def track(verb: str, unit: str, total: Optional[int] = None, *,
          past: str = "Finished", _progress: Any = None):
    """Return a context manager that shows progress for one operation.

    Parameters
    ----------
    verb : present-tense label, for example ``"Peeling layers"``.
    unit : singular counter name, for example ``"layer"``.
    total : expected count, or ``None``. ``%`` shows only when known.
    past : past-tense verb for the done line, for example ``"Peeled"``.
    _progress : private per-call override of ``CHILMESH_PROGRESS``.
    """
    mode = _resolve_mode(_progress)
    if mode == "0" or (mode == "auto" and not _is_tty()):
        return _NullTracker()
    return _Tracker(verb, unit, total, past, forced=(mode == "1"))
