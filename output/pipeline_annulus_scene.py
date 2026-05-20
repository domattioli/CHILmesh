"""Manim animation: CHILmesh README pipeline on the annulus.

Four stages from the legacy 4-row README:
  Row 1: Raw annulus (Delaunay input)
  Row 2: ADMESH truss warm-start (spring relaxation vs SDF)
  Row 3: FEM smoother (Balendran direct method)
  Row 4: Right-iso smoother (angles -> 45/45/90)

Each stage has its own (potentially different) triangulation; Delaunay
re-triangulation runs after every smoother to repair flipped tris from
free-vertex motion. Transitions use fade-out/fade-in instead of
morphing because topology differs between rows.

Metrics: median element quality (higher better, green bar) and mean
per-triangle deviation from right-isoceles {45,45,90} angles (lower
better, orange bar).

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
    config,
)


config.frame_width = 16
config.frame_height = 9
config.pixel_width = 1920
config.pixel_height = 1080
config.frame_rate = 30
config.background_color = "#0e0e10"

DATA_PATH = Path("/tmp/readme_pipeline.json")

EDGE_COLOR = "#3a7fbf"
FILL_COLOR = "#5fb0ff"
ACCENT     = "#ff9f43"
GOOD       = "#2ecc71"
DIM        = "#666666"
BG_TILE    = "#1a1a1f"


def _load():
    return json.loads(DATA_PATH.read_text())


class ReadmePipeline(Scene):
    def construct(self):
        data = _load()
        rows = data["rows"]
        n_stages = len(rows)

        # Title
        title = Text("CHILmesh annulus pipeline", font_size=44,
                     weight="BOLD").to_edge(UP, buff=0.30)
        subtitle = Text(
            "Row 1 -> Row 2 -> Row 3 -> Row 4   (legacy README layout)",
            font_size=20, color="#aaaaaa")
        subtitle.next_to(title, DOWN, buff=0.12)
        self.play(Write(title), FadeIn(subtitle), run_time=0.8)

        # Compute common mesh scale across all stages.
        all_pts = np.vstack([np.asarray(r["points"]) for r in rows])
        cx, cy = all_pts.mean(axis=0)
        max_extent = max(np.ptp(all_pts[:, 0]), np.ptp(all_pts[:, 1]))
        mesh_scale = 5.2 / max_extent
        mesh_center = np.array([-3.5, -0.5, 0.0])

        def make_mesh(pts, elements):
            pts = np.asarray(pts)
            def to_scene(i):
                return np.array([
                    (pts[i, 0] - cx) * mesh_scale + mesh_center[0],
                    (pts[i, 1] - cy) * mesh_scale + mesh_center[1],
                    0.0,
                ])
            polys = VGroup()
            for elem in elements:
                p = Polygon(
                    *[to_scene(i) for i in elem],
                    stroke_color=EDGE_COLOR,
                    stroke_width=0.9,
                    fill_color=FILL_COLOR,
                    fill_opacity=0.22,
                )
                polys.add(p)
            return polys

        # Pipeline progress bar (top-right of mesh area).
        bar_w, bar_h = 4.6, 0.45
        bar_origin = np.array([3.6, 2.5, 0.0])
        progress_bg = Rectangle(width=bar_w, height=bar_h,
                                 stroke_color=DIM, stroke_width=1.5,
                                 fill_color=BG_TILE, fill_opacity=1.0)
        progress_bg.move_to(bar_origin)
        ticks, tick_labels = VGroup(), VGroup()
        for k in range(n_stages):
            x = bar_origin[0] - bar_w / 2 + (k + 0.5) * bar_w / n_stages
            tick = Line(np.array([x, bar_origin[1] - bar_h/2 - 0.05, 0]),
                        np.array([x, bar_origin[1] + bar_h/2 + 0.05, 0]),
                        color=DIM, stroke_width=1.5)
            ticks.add(tick)
            lbl = Text(rows[k]["stage"], font_size=16, color="#aaaaaa")
            lbl.next_to(tick, UP, buff=0.10)
            tick_labels.add(lbl)
        self.play(Create(progress_bg), FadeIn(ticks), FadeIn(tick_labels), run_time=0.6)

        stage_name_pos = np.array([3.6, 1.6, 0.0])
        algo_pos = np.array([3.6, 1.2, 0.0])
        stage_name = Text("", font_size=28, color=ACCENT, weight="BOLD").move_to(stage_name_pos)
        algo_text = Text("", font_size=20, color="#cccccc").move_to(algo_pos)
        self.add(stage_name, algo_text)

        # Metric bars.
        q_axis_origin = np.array([2.6, -0.4, 0.0])
        bar_max = 4.0
        q_label = Text("median quality (-> 1.0 better)", font_size=18, color="#cccccc")
        q_label.move_to(q_axis_origin + np.array([0, 0.45, 0])).align_to(q_axis_origin, LEFT)
        q_track = Rectangle(width=bar_max, height=0.30,
                             stroke_color=DIM, stroke_width=1.2,
                             fill_color=BG_TILE, fill_opacity=1.0)
        q_track.move_to(q_axis_origin + np.array([bar_max/2, 0, 0]))

        iso_axis_origin = np.array([2.6, -1.4, 0.0])
        iso_label = Text("right-iso angle deviation (<- lower better)",
                          font_size=18, color="#cccccc")
        iso_label.move_to(iso_axis_origin + np.array([0, 0.45, 0])).align_to(iso_axis_origin, LEFT)
        iso_track = Rectangle(width=bar_max, height=0.30,
                               stroke_color=DIM, stroke_width=1.2,
                               fill_color=BG_TILE, fill_opacity=1.0)
        iso_track.move_to(iso_axis_origin + np.array([bar_max/2, 0, 0]))

        self.play(FadeIn(q_label), FadeIn(q_track),
                  FadeIn(iso_label), FadeIn(iso_track), run_time=0.5)

        def q_bar_for(value):
            w = max(0.02, value) * bar_max
            r = Rectangle(width=w, height=0.30, stroke_width=0,
                          fill_color=GOOD, fill_opacity=0.95)
            r.move_to(q_axis_origin + np.array([w/2, 0, 0]))
            return r

        def iso_bar_for(dev):
            frac = 1.0 - min(1.0, dev / 45.0)
            w = max(0.02, frac) * bar_max
            r = Rectangle(width=w, height=0.30, stroke_width=0,
                          fill_color=ACCENT, fill_opacity=0.95)
            r.move_to(iso_axis_origin + np.array([w/2, 0, 0]))
            return r

        # Show Row 1.
        mesh = make_mesh(rows[0]["points"], rows[0]["elements"])
        self.play(Create(mesh, lag_ratio=0.003), run_time=2.0)
        new_stage = Text(rows[0]["name"], font_size=28, color=ACCENT, weight="BOLD").move_to(stage_name_pos)
        new_algo = Text(rows[0]["algo"], font_size=20, color="#cccccc").move_to(algo_pos)
        self.play(Transform(stage_name, new_stage), Transform(algo_text, new_algo), run_time=0.4)

        q_bar = q_bar_for(rows[0]["q_med"])
        iso_bar = iso_bar_for(rows[0]["iso_dev"])
        q_readout = Text(f"med Q = {rows[0]['q_med']:.3f}", font_size=22, color=GOOD)
        q_readout.next_to(q_track, RIGHT, buff=0.20)
        iso_readout = Text(f"dev = {rows[0]['iso_dev']:.1f} deg", font_size=22, color=ACCENT)
        iso_readout.next_to(iso_track, RIGHT, buff=0.20)
        self.play(FadeIn(q_bar), FadeIn(iso_bar),
                  FadeIn(q_readout), FadeIn(iso_readout), run_time=0.5)

        progress_marker = ticks[0].copy().set_color(ACCENT).set_stroke(width=3.5)
        self.add(progress_marker)
        self.wait(1.4)

        # Walk through remaining stages.
        for k in range(1, n_stages):
            new_mesh = make_mesh(rows[k]["points"], rows[k]["elements"])
            new_q_bar = q_bar_for(rows[k]["q_med"])
            new_iso_bar = iso_bar_for(rows[k]["iso_dev"])
            new_q_readout = Text(f"med Q = {rows[k]['q_med']:.3f}", font_size=22, color=GOOD)
            new_q_readout.next_to(q_track, RIGHT, buff=0.20)
            new_iso_readout = Text(f"dev = {rows[k]['iso_dev']:.1f} deg", font_size=22, color=ACCENT)
            new_iso_readout.next_to(iso_track, RIGHT, buff=0.20)
            new_stage = Text(rows[k]["name"], font_size=28, color=ACCENT, weight="BOLD").move_to(stage_name_pos)
            new_algo = Text(rows[k]["algo"], font_size=20, color="#cccccc").move_to(algo_pos)
            new_marker = ticks[k].copy().set_color(ACCENT).set_stroke(width=3.5)

            # Fade out old mesh, fade in new (topology may differ).
            self.play(FadeOut(mesh), run_time=0.35)
            self.play(
                FadeIn(new_mesh),
                Transform(q_bar, new_q_bar),
                Transform(iso_bar, new_iso_bar),
                Transform(q_readout, new_q_readout),
                Transform(iso_readout, new_iso_readout),
                Transform(stage_name, new_stage),
                Transform(algo_text, new_algo),
                Transform(progress_marker, new_marker),
                run_time=1.0,
            )
            mesh = new_mesh
            self.wait(1.5)

        delta_q = rows[2]["q_med"] - rows[0]["q_med"]
        delta_iso = rows[0]["iso_dev"] - rows[-1]["iso_dev"]
        caption = Text(
            f"truss + FEM lift Q from {rows[0]['q_med']:.2f} -> {rows[2]['q_med']:.2f}   |   "
            f"right-iso prep cuts angle dev by {delta_iso:.1f} deg",
            font_size=18, color="#cccccc"
        )
        caption.to_edge(DOWN, buff=0.3)
        self.play(FadeIn(caption), run_time=0.6)
        self.wait(2.5)
