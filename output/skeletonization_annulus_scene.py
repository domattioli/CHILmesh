"""Manim animation: CHILmesh skeletonization on the annulus fixture.

Renders 1080p mp4 suitable for embedding in README.md.

  - faint grey: full mesh
  - orange fill: OE (outer elements of current layer)
  - green fill:  IE (inner elements of current layer)
  - blue dots:   OV (outer vertices of current layer)
  - title bar:   "Layer k"
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
    Dot,
    Square,
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


config.frame_width = 16
config.frame_height = 9
config.pixel_width = 1920
config.pixel_height = 1080
config.frame_rate = 30
config.background_color = "#0e0e10"

DATA_PATH = Path("/tmp/skeleton_annulus.json")

OV_COLOR = "#5fb0ff"
OE_COLOR = "#ff9f43"
IE_COLOR = "#2ecc71"
DIM_COLOR = "#444444"
CONSUMED_COLOR = "#262629"


def _load():
    return json.loads(DATA_PATH.read_text())


class SkeletonizationAnnulus(Scene):
    def construct(self):
        data = _load()
        pts = np.asarray(data["points"])
        elements = data["elements"]
        layers = data["layers"]

        cx, cy = pts.mean(axis=0)
        max_extent = max(np.ptp(pts[:, 0]), np.ptp(pts[:, 1]))
        # Leave room for title (top) and legend (bottom).
        scale = 6.6 / max_extent

        def to_scene(i):
            return np.array([
                (pts[i, 0] - cx) * scale,
                (pts[i, 1] - cy) * scale - 0.3,
                0.0,
            ])

        # Background mesh polygons (kept on scene the whole time).
        mesh_polys = VGroup()
        for elem in elements:
            poly = Polygon(
                *[to_scene(i) for i in elem],
                stroke_color=DIM_COLOR,
                stroke_width=0.9,
                fill_opacity=0.0,
            )
            mesh_polys.add(poly)

        title = Text("CHILmesh skeletonization", font_size=44,
                     weight="BOLD").to_edge(UP, buff=0.35)
        subtitle = Text(
            f"annulus fixture  -  {len(elements)} elements  -  {len(layers)} layers",
            font_size=22, color="#aaaaaa",
        )
        subtitle.next_to(title, DOWN, buff=0.18)

        self.play(Write(title), run_time=0.8)
        self.play(FadeIn(subtitle), run_time=0.4)
        self.play(Create(mesh_polys, lag_ratio=0.003), run_time=2.5)
        self.wait(0.4)

        legend = self._legend()
        legend.to_edge(DOWN, buff=0.3)
        self.play(FadeIn(legend), run_time=0.5)

        # Layer counter (top-right of mesh area).
        layer_label = Text("", font_size=34, color=OV_COLOR)
        layer_label.to_corner(UP + RIGHT, buff=0.8).shift(DOWN * 0.4)
        self.add(layer_label)

        for k, layer in enumerate(layers):
            new_label = Text(f"Layer {k}", font_size=34, color=OV_COLOR, weight="BOLD")
            new_label.move_to(layer_label.get_center())
            self.play(Transform(layer_label, new_label), run_time=0.35)

            oe_polys = VGroup(*[
                Polygon(
                    *[to_scene(i) for i in elements[eid]],
                    stroke_color=OE_COLOR,
                    stroke_width=1.4,
                    fill_color=OE_COLOR,
                    fill_opacity=0.55,
                )
                for eid in layer["OE"]
            ])
            ie_polys = VGroup(*[
                Polygon(
                    *[to_scene(i) for i in elements[eid]],
                    stroke_color=IE_COLOR,
                    stroke_width=1.2,
                    fill_color=IE_COLOR,
                    fill_opacity=0.45,
                )
                for eid in layer["IE"]
            ])
            ov_dots = VGroup(*[
                Dot(to_scene(vid), color=OV_COLOR, radius=0.055)
                for vid in layer["OV"]
            ])

            self.play(FadeIn(oe_polys, lag_ratio=0.005), run_time=0.9)
            self.play(FadeIn(ie_polys, lag_ratio=0.005), run_time=0.8)
            self.play(Create(ov_dots, lag_ratio=0.01), run_time=0.7)
            self.wait(1.2)

            # Consume: dim everything in this layer before peeling next.
            self.play(
                oe_polys.animate.set_stroke(color=CONSUMED_COLOR, width=0.6).set_fill(
                    color=CONSUMED_COLOR, opacity=0.18),
                ie_polys.animate.set_stroke(color=CONSUMED_COLOR, width=0.6).set_fill(
                    color=CONSUMED_COLOR, opacity=0.18),
                ov_dots.animate.set_color(CONSUMED_COLOR).scale(0.55),
                run_time=0.55,
            )

        done = Text(f"{len(layers)} layers peeled",
                    font_size=34, color=OV_COLOR, weight="BOLD")
        done.move_to(layer_label.get_center())
        self.play(Transform(layer_label, done), run_time=0.4)
        self.wait(2.2)

    def _legend(self):
        items = VGroup(
            self._swatch(OE_COLOR, "OE  outer elements", fill=True),
            self._swatch(IE_COLOR, "IE  inner elements", fill=True),
            self._swatch(OV_COLOR, "OV  outer vertices", fill=False),
        ).arrange(RIGHT, buff=0.9)
        return items

    def _swatch(self, color, label, *, fill):
        if fill:
            mark = Square(side_length=0.30, color=color,
                          fill_color=color, fill_opacity=0.55, stroke_width=1.5)
        else:
            mark = Dot(color=color, radius=0.13)
        txt = Text(label, font_size=20, color="#dddddd")
        txt.next_to(mark, RIGHT, buff=0.18)
        return VGroup(mark, txt)
