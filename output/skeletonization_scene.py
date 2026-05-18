"""Manim animation: CHILmesh skeletonization on a 7x7 toy mesh.

Visualizes the iterative layer-peeling from boundary inward:
  - faint grey: full mesh outline
  - blue: OV (outer vertices of current layer)
  - orange: OE (outer elements)
  - green: IE (inner elements)
  - title labels each layer index

Renders to /home/user/CHILmesh/output/skeletonization_toy.mp4
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np
from manim import (
    DOWN,
    UP,
    Polygon,
    Dot,
    Scene,
    Text,
    VGroup,
    Write,
    Create,
    FadeIn,
    FadeOut,
    Transform,
    config,
)


config.frame_width = 14
config.frame_height = 9
config.pixel_width = 1280
config.pixel_height = 720
config.background_color = "#0e0e10"

DATA_PATH = Path("/tmp/skeleton_data.json")

OV_COLOR = "#5fb0ff"     # outer vertices
OE_COLOR = "#ff9f43"     # outer elements
IE_COLOR = "#2ecc71"     # inner elements
DIM_COLOR = "#555555"


def _load():
    return json.loads(DATA_PATH.read_text())


def _xy_to_scene(xy, *, scale=1.0, offset=(-3.0, -3.0)):
    return np.array([xy[0] * scale + offset[0], xy[1] * scale + offset[1], 0.0])


class Skeletonization(Scene):
    def construct(self):
        data = _load()
        pts = np.asarray(data["points"])
        elements = data["elements"]
        layers = data["layers"]

        # Center the mesh in the frame.
        cx, cy = pts.mean(axis=0)
        max_extent = max(np.ptp(pts[:, 0]), np.ptp(pts[:, 1]))
        scale = 6.0 / max_extent
        offset = (-cx * scale, -cy * scale)

        def to_scene(i):
            return _xy_to_scene(pts[i], scale=scale, offset=offset)

        # Faint full-mesh outline (every element as thin polygon).
        mesh_polys = VGroup()
        for elem in elements:
            verts = [to_scene(i) for i in elem]
            poly = Polygon(
                *verts,
                stroke_color=DIM_COLOR,
                stroke_width=1.2,
                fill_opacity=0.0,
            )
            mesh_polys.add(poly)

        title = Text("CHILmesh skeletonization", font_size=32).to_edge(UP)
        subtitle = Text("7x7 structured triangulation", font_size=20, color=DIM_COLOR)
        subtitle.next_to(title, DOWN, buff=0.15)

        self.play(Write(title), FadeIn(subtitle))
        self.play(Create(mesh_polys, lag_ratio=0.02, run_time=2.0))
        self.wait(0.5)

        legend = VGroup(
            self._legend_swatch(OE_COLOR, "OE  (outer elements)"),
            self._legend_swatch(IE_COLOR, "IE  (inner elements)"),
            self._legend_swatch(OV_COLOR, "OV  (outer vertices)", dot=True),
        ).arrange(DOWN, aligned_edge=0).to_edge(DOWN).shift(UP * 0.2)
        self.play(FadeIn(legend))
        self.wait(0.4)

        layer_label = Text("", font_size=28).next_to(subtitle, DOWN, buff=0.3)
        self.add(layer_label)
        prev_highlights = VGroup()

        for k, layer in enumerate(layers):
            new_label = Text(f"Layer {k}", font_size=28, color=OV_COLOR)
            new_label.move_to(layer_label.get_center())
            self.play(Transform(layer_label, new_label), run_time=0.4)

            oe_polys = VGroup(*[
                Polygon(
                    *[to_scene(i) for i in elements[eid]],
                    stroke_color=OE_COLOR,
                    stroke_width=2.5,
                    fill_color=OE_COLOR,
                    fill_opacity=0.35,
                )
                for eid in layer["OE"]
            ])
            ie_polys = VGroup(*[
                Polygon(
                    *[to_scene(i) for i in elements[eid]],
                    stroke_color=IE_COLOR,
                    stroke_width=2.0,
                    fill_color=IE_COLOR,
                    fill_opacity=0.30,
                )
                for eid in layer["IE"]
            ])
            ov_dots = VGroup(*[
                Dot(to_scene(vid), color=OV_COLOR, radius=0.10)
                for vid in layer["OV"]
            ])

            self.play(FadeIn(oe_polys, lag_ratio=0.03), run_time=1.0)
            self.play(FadeIn(ie_polys, lag_ratio=0.03), run_time=0.9)
            self.play(Create(ov_dots, lag_ratio=0.04), run_time=0.9)
            self.wait(1.0)

            # Fade the layer down to "consumed" grey before peeling next.
            consumed = VGroup(oe_polys, ie_polys)
            self.play(
                consumed.animate.set_stroke(color=DIM_COLOR, width=1.0).set_fill(opacity=0.0),
                ov_dots.animate.set_color(DIM_COLOR).scale(0.6),
                run_time=0.6,
            )
            prev_highlights.add(consumed, ov_dots)

        done = Text(f"{len(layers)} layers peeled", font_size=26, color=OV_COLOR)
        done.move_to(layer_label.get_center())
        self.play(Transform(layer_label, done))
        self.wait(1.5)

    def _legend_swatch(self, color, label, *, dot=False):
        from manim import Square, RIGHT
        if dot:
            swatch = Dot(color=color, radius=0.10)
        else:
            swatch = Square(side_length=0.30, color=color, fill_opacity=0.4,
                            stroke_width=2.0)
        txt = Text(label, font_size=18, color="#dddddd")
        txt.next_to(swatch, RIGHT, buff=0.2)
        return VGroup(swatch, txt)
