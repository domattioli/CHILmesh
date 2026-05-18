"""Manim animation: iterative Laplacian smoothing on a perturbed annulus.

Shows mesh morphing toward equilibrium while a quality trace grows.
1080p30 mp4 for README embedding.
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
    ORIGIN,
    Polygon,
    Scene,
    Text,
    VGroup,
    Write,
    Create,
    FadeIn,
    FadeOut,
    Transform,
    Axes,
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

DATA_PATH = Path("/tmp/smooth_annulus.json")

MESH_COLOR = "#5fb0ff"
EDGE_COLOR = "#3a7fbf"
PLOT_COLOR = "#ff9f43"
PERTURB_COLOR = "#e74c3c"
SMOOTH_COLOR = "#2ecc71"


def _load():
    return json.loads(DATA_PATH.read_text())


class SmoothingAnnulus(Scene):
    def construct(self):
        data = _load()
        elements = data["elements"]
        snapshots = [np.asarray(s) for s in data["snapshots"]]
        qualities = data["qualities"]
        n_iter = len(snapshots) - 1

        # Layout: mesh on left half, quality plot on right.
        mesh_center = np.array([-3.4, -0.3, 0.0])
        all_pts = np.vstack(snapshots)
        cx, cy = all_pts.mean(axis=0)
        max_extent = max(np.ptp(all_pts[:, 0]), np.ptp(all_pts[:, 1]))
        mesh_scale = 5.4 / max_extent

        def make_mesh(pts):
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
                    stroke_width=1.0,
                    fill_color=MESH_COLOR,
                    fill_opacity=0.25,
                )
                polys.add(p)
            return polys

        # Title + labels.
        title = Text("FEM-style Laplacian smoothing", font_size=42,
                     weight="BOLD").to_edge(UP, buff=0.35)
        subtitle = Text(
            "annulus, 380 verts perturbed then iteratively averaged",
            font_size=20, color="#aaaaaa")
        subtitle.next_to(title, DOWN, buff=0.15)

        self.play(Write(title), FadeIn(subtitle), run_time=0.9)

        # Initial (perturbed) mesh.
        mesh = make_mesh(snapshots[0])
        perturb_label = Text("perturbed (random offset)", font_size=20,
                              color=PERTURB_COLOR)
        perturb_label.move_to(mesh_center + np.array([0, -3.5, 0]))
        self.play(Create(mesh, lag_ratio=0.003, run_time=2.0))
        self.play(FadeIn(perturb_label), run_time=0.3)
        self.wait(0.6)

        # Quality plot setup.
        plot_origin = np.array([4.4, -0.5, 0.0])
        axes = Axes(
            x_range=[0, n_iter, max(1, n_iter // 5)],
            y_range=[0.30, 0.75, 0.1],
            x_length=5.0,
            y_length=4.0,
            tips=False,
            axis_config={"color": "#888888", "stroke_width": 1.5},
        ).move_to(plot_origin)
        x_label = Text("iteration", font_size=20, color="#cccccc")
        x_label.next_to(axes.x_axis, DOWN, buff=0.25)
        y_label = Text("avg mesh quality", font_size=20, color="#cccccc")
        y_label.rotate(np.pi / 2).next_to(axes.y_axis, LEFT, buff=0.25)
        plot_box = VGroup(axes, x_label, y_label)
        self.play(FadeIn(plot_box), run_time=0.6)

        # Quality trace: start with single dot, extend each iteration.
        traced_points = [axes.c2p(0, qualities[0])]
        trace_dot_start = Dot(traced_points[0], color=PLOT_COLOR, radius=0.06)
        self.add(trace_dot_start)

        # Quality readout next to plot.
        q_text = Text(f"q = {qualities[0]:.3f}", font_size=24, color=PLOT_COLOR)
        q_text.next_to(axes, UP, buff=0.25)
        self.play(FadeIn(q_text), FadeOut(perturb_label), run_time=0.5)

        smooth_label = Text("smoothing...", font_size=20, color=SMOOTH_COLOR)
        smooth_label.move_to(mesh_center + np.array([0, -3.5, 0]))
        self.play(FadeIn(smooth_label), run_time=0.3)

        # Iterate: morph mesh + extend trace.
        trace_segments = VGroup()
        prev_run_time = 0.35
        for i in range(1, n_iter + 1):
            new_mesh = make_mesh(snapshots[i])
            new_point = axes.c2p(i, qualities[i])
            seg = Line(traced_points[-1], new_point,
                       stroke_color=PLOT_COLOR, stroke_width=3.0)
            dot = Dot(new_point, color=PLOT_COLOR, radius=0.05)

            # Pace: fast early (big changes), slower near end.
            if i <= 3:
                rt = 0.6
            elif i <= 10:
                rt = 0.25
            else:
                rt = 0.12

            new_q_text = Text(f"q = {qualities[i]:.3f}",
                              font_size=24, color=PLOT_COLOR)
            new_q_text.move_to(q_text.get_center())
            self.play(
                Transform(mesh, new_mesh),
                Create(seg),
                FadeIn(dot),
                Transform(q_text, new_q_text),
                run_time=rt,
            )
            traced_points.append(new_point)
            trace_segments.add(seg, dot)

        # Final state.
        done = Text("converged", font_size=22, color=SMOOTH_COLOR)
        done.move_to(smooth_label.get_center())
        delta = qualities[-1] - qualities[0]
        delta_text = Text(
            f"avg quality:  {qualities[0]:.3f}  ->  {qualities[-1]:.3f}   "
            f"(+{delta:.3f})",
            font_size=22, color="#dddddd",
        )
        delta_text.next_to(done, DOWN, buff=0.25)
        self.play(Transform(smooth_label, done), FadeIn(delta_text), run_time=0.6)
        self.wait(2.4)
