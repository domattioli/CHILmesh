"""Manim animation: quality distributions across CHILmesh's 4 built-in fixtures.

Renders a side-by-side comparison of element-quality distributions for
the annulus, donut, structured, and block_o fixtures. Each fixture
gets a horizontal lane with:

  * a histogram of element-quality values (filled bars),
  * box-plot markers (min, Q25, median, Q75, max) overlaid,
  * thumbnail of the mesh on the left (omitted for block_o due to
    polygon count).

Intent: one glance answers "how do these meshes compare in quality?"

Renders 1080p30 mp4.
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np
from manim import (
    DOWN,
    UP,
    LEFT,
    RIGHT,
    Polygon,
    Rectangle,
    Scene,
    Text,
    VGroup,
    Write,
    Create,
    FadeIn,
    FadeOut,
    Transform,
    Line,
    Dot,
    config,
)


config.frame_width = 16
config.frame_height = 9
config.pixel_width = 1920
config.pixel_height = 1080
config.frame_rate = 30
config.background_color = "#0e0e10"

DATA_PATH = Path("/tmp/quality_compare.json")

# Per-fixture color palette.
PALETTE = {
    "annulus":    "#5fb0ff",
    "donut":      "#ff9f43",
    "structured": "#2ecc71",
    "block_o":    "#e74c3c",
}
DIM     = "#666666"
BG_TILE = "#1a1a1f"
GOOD    = "#2ecc71"
TXT     = "#dddddd"
ACCENT  = "#ffd34e"


def _load():
    return json.loads(DATA_PATH.read_text())


class QualityCompare(Scene):
    def construct(self):
        data = _load()
        order = ["annulus", "donut", "structured", "block_o"]
        # Sort by median quality ascending so the eye reads worst -> best.
        order.sort(key=lambda n: data[n]["q_med"])

        title = Text("Element quality across CHILmesh fixtures",
                     font_size=42, weight="BOLD").to_edge(UP, buff=0.30)
        subtitle = Text(
            "histogram + min / Q25 / median / Q75 / max markers; one lane per fixture",
            font_size=20, color="#aaaaaa")
        subtitle.next_to(title, DOWN, buff=0.12)
        self.play(Write(title), FadeIn(subtitle), run_time=0.8)

        # Lane geometry.
        lane_y_top = 2.6
        lane_y_bot = -3.0
        n = len(order)
        lane_h = (lane_y_top - lane_y_bot) / n
        x_axis_left = -3.5
        x_axis_right = 6.0
        x_axis_w = x_axis_right - x_axis_left

        # Shared x-axis at the bottom.
        x_axis_y = lane_y_bot - 0.15
        x_axis = Line(np.array([x_axis_left, x_axis_y, 0]),
                       np.array([x_axis_right, x_axis_y, 0]),
                       color=DIM, stroke_width=1.5)
        x_ticks = VGroup()
        x_labels = VGroup()
        for v in [0.0, 0.25, 0.5, 0.75, 1.0]:
            tx = x_axis_left + v * x_axis_w
            tk = Line(np.array([tx, x_axis_y - 0.10, 0]),
                      np.array([tx, x_axis_y + 0.10, 0]),
                      color=DIM, stroke_width=1.5)
            lbl = Text(f"{v:.2f}", font_size=18, color="#bbbbbb")
            lbl.next_to(tk, DOWN, buff=0.10)
            x_ticks.add(tk); x_labels.add(lbl)
        x_title = Text("element quality (0 = degenerate, 1 = ideal)",
                        font_size=20, color=TXT)
        x_title.next_to(x_labels, DOWN, buff=0.30)
        self.play(Create(x_axis), FadeIn(x_ticks), FadeIn(x_labels),
                  FadeIn(x_title), run_time=0.7)

        def q_to_scene_x(q):
            return x_axis_left + float(q) * x_axis_w

        # Build all four lanes.
        N_BINS = 30
        bins = np.linspace(0.0, 1.0, N_BINS + 1)
        for idx, name in enumerate(order):
            color = PALETTE[name]
            d = data[name]
            qarr = np.asarray(d["qualities"])
            counts, _ = np.histogram(qarr, bins=bins)
            # Normalise per-lane so different mesh sizes share visual scale.
            max_count = max(1, counts.max())
            lane_center_y = lane_y_top - lane_h * idx - lane_h / 2
            bar_max_h = lane_h * 0.65
            baseline_y = lane_center_y - bar_max_h / 2

            # Histogram bars.
            hist_group = VGroup()
            for b in range(N_BINS):
                if counts[b] == 0:
                    continue
                w = x_axis_w / N_BINS
                h = (counts[b] / max_count) * bar_max_h
                x = q_to_scene_x((bins[b] + bins[b+1]) / 2)
                rect = Rectangle(width=w * 0.92, height=h,
                                  stroke_width=0,
                                  fill_color=color, fill_opacity=0.55)
                rect.move_to(np.array([x, baseline_y + h / 2, 0]))
                hist_group.add(rect)

            # Box-plot markers above bars.
            box_y = lane_center_y + bar_max_h / 2 + 0.18
            iqr_line = Line(
                np.array([q_to_scene_x(d["q_q25"]), box_y, 0]),
                np.array([q_to_scene_x(d["q_q75"]), box_y, 0]),
                color=color, stroke_width=8,
            )
            whisker = Line(
                np.array([q_to_scene_x(d["q_min"]), box_y, 0]),
                np.array([q_to_scene_x(d["q_max"]), box_y, 0]),
                color=color, stroke_width=2,
            )
            med_dot = Dot(np.array([q_to_scene_x(d["q_med"]), box_y, 0]),
                           color="#ffffff", radius=0.10)

            # Lane label (fixture name + counts).
            name_lbl = Text(name, font_size=26, color=color, weight="BOLD")
            stats_lbl = Text(
                f"{d['n_elems']:,} elems   med {d['q_med']:.3f}   min {d['q_min']:.3f}",
                font_size=16, color="#bbbbbb",
            )
            name_lbl.move_to(np.array([-6.6, lane_center_y + 0.20, 0]))
            stats_lbl.move_to(np.array([-6.6, lane_center_y - 0.18, 0]))
            name_lbl.align_to(np.array([-6.6, 0, 0]), LEFT)
            stats_lbl.align_to(np.array([-6.6, 0, 0]), LEFT)

            # Animate this lane in: label first, then hist, then box markers.
            self.play(FadeIn(name_lbl), FadeIn(stats_lbl), run_time=0.30)
            self.play(Create(hist_group, lag_ratio=0.02), run_time=0.85)
            self.play(Create(whisker), Create(iqr_line), FadeIn(med_dot),
                      run_time=0.45)
            self.wait(0.35)

        # Final commentary banner.
        msgs = [
            ("block_o", "highest-quality structured tri mesh"),
            ("structured", "uniform 4x4 quad grid -> tight quality cluster"),
            ("donut", "moderate quality, two boundary rings"),
            ("annulus", "ragged baseline -- long tail down to ~0"),
        ]
        # Show order matches sorted-by-median ascending.
        caption_lines = VGroup()
        for n_, msg in msgs:
            t = Text(f"{n_}: {msg}", font_size=18, color=PALETTE[n_])
            caption_lines.add(t)
        # We don't have room for all 4; instead a single summary.
        summary = Text(
            "ranking by median Q:   "
            + "  <  ".join(f"{n} ({data[n]['q_med']:.2f})" for n in order),
            font_size=20, color=TXT,
        )
        summary.to_edge(DOWN, buff=0.25)
        self.play(FadeIn(summary), run_time=0.6)
        self.wait(2.5)
