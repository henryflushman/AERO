"""Create CZML and standalone CesiumJS orbit visualizations from Python data.

Cesium Cartesian positions are written in meters. Input positions to this module
use km to match the rest of the orbital mechanics helpers.
"""

from __future__ import annotations

import json
from dataclasses import dataclass, field
from datetime import datetime, timedelta, timezone
from pathlib import Path
from typing import Iterable

import numpy as np


def _utc(time: datetime) -> datetime:
    if not isinstance(time, datetime):
        raise TypeError("trajectory times must be datetime objects.")
    if time.tzinfo is None:
        return time.replace(tzinfo=timezone.utc)
    return time.astimezone(timezone.utc)


def _iso(time: datetime) -> str:
    return _utc(time).isoformat().replace("+00:00", "Z")


def _rgba(color) -> list[int]:
    values = list(color)
    if len(values) == 3:
        values.append(255)
    if len(values) != 4:
        raise ValueError("color must be RGB or RGBA.")
    values = [int(v) for v in values]
    if any(v < 0 or v > 255 for v in values):
        raise ValueError("RGBA values must lie in [0, 255].")
    return values


def trajectory_packet(
    times: Iterable[datetime],
    positions_km,
    *,
    entity_id: str,
    name: str | None = None,
    reference_frame: str = "INERTIAL",
    color=(0, 170, 255, 255),
    path_width: float = 2.0,
    point_pixel_size: float = 8.0,
) -> dict:
    """Build one time-dynamic CZML trajectory packet."""
    times = [_utc(t) for t in times]
    positions = np.asarray(positions_km, dtype=float)
    if positions.ndim != 2 or positions.shape[1] != 3:
        raise ValueError("positions_km must have shape (N, 3).")
    if len(times) != len(positions):
        raise ValueError("times and positions_km must have the same length.")
    if not times:
        raise ValueError("trajectory must contain at least one sample.")
    if any(times[i] > times[i + 1] for i in range(len(times) - 1)):
        raise ValueError("times must be sorted in ascending order.")

    frame = reference_frame.upper()
    if frame not in {"INERTIAL", "FIXED"}:
        raise ValueError("reference_frame must be 'INERTIAL' or 'FIXED'.")

    epoch = times[0]
    samples: list[float] = []
    for time, position_km in zip(times, positions):
        dt_s = (time - epoch).total_seconds()
        samples.extend([dt_s, *(1000.0 * position_km)])

    rgba = _rgba(color)
    availability = f"{_iso(times[0])}/{_iso(times[-1])}"
    return {
        "id": str(entity_id),
        "name": name or str(entity_id),
        "availability": availability,
        "position": {
            "epoch": _iso(epoch),
            "referenceFrame": frame,
            "interpolationAlgorithm": "LAGRANGE",
            "interpolationDegree": 5,
            "cartesian": samples,
        },
        "point": {
            "pixelSize": float(point_pixel_size),
            "color": {"rgba": rgba},
            "outlineColor": {"rgba": [255, 255, 255, 255]},
            "outlineWidth": 1.0,
        },
        "path": {
            "show": True,
            "width": float(path_width),
            "material": {"solidColor": {"color": {"rgba": rgba}}},
            "leadTime": 0,
            "trailTime": max(1.0, (times[-1] - times[0]).total_seconds()),
            "resolution": 30,
        },
        "label": {
            "text": name or str(entity_id),
            "font": "14px sans-serif",
            "showBackground": True,
            "pixelOffset": {"cartesian2": [12, -12]},
        },
    }


@dataclass
class CesiumScene:
    """Collect one or more orbit trajectories and write CZML/HTML."""

    name: str = "Intermediate Orbits"
    multiplier: float = 60.0
    _packets: list[dict] = field(default_factory=list, init=False, repr=False)
    _start: datetime | None = field(default=None, init=False, repr=False)
    _stop: datetime | None = field(default=None, init=False, repr=False)

    def add_trajectory(
        self,
        times: Iterable[datetime],
        positions_km,
        *,
        entity_id: str,
        name: str | None = None,
        reference_frame: str = "INERTIAL",
        color=(0, 170, 255, 255),
        path_width: float = 2.0,
        point_pixel_size: float = 8.0,
    ) -> dict:
        times = [_utc(t) for t in times]
        packet = trajectory_packet(
            times,
            positions_km,
            entity_id=entity_id,
            name=name,
            reference_frame=reference_frame,
            color=color,
            path_width=path_width,
            point_pixel_size=point_pixel_size,
        )
        self._packets.append(packet)
        self._start = times[0] if self._start is None else min(self._start, times[0])
        self._stop = times[-1] if self._stop is None else max(self._stop, times[-1])
        return packet

    @property
    def czml(self) -> list[dict]:
        if not self._packets:
            raise RuntimeError("Add at least one trajectory before requesting CZML.")
        document = {
            "id": "document",
            "name": self.name,
            "version": "1.0",
            "clock": {
                "interval": f"{_iso(self._start)}/{_iso(self._stop)}",
                "currentTime": _iso(self._start),
                "multiplier": float(self.multiplier),
                "range": "LOOP_STOP",
                "step": "SYSTEM_CLOCK_MULTIPLIER",
            },
        }
        return [document, *self._packets]

    def save_czml(self, path) -> Path:
        path = Path(path)
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(json.dumps(self.czml, indent=2), encoding="utf-8")
        return path

    def save_html(
        self,
        path,
        *,
        cesium_version: str = "latest",
        title: str | None = None,
    ) -> Path:
        """Write a standalone HTML viewer with the CZML embedded directly.

        The browser needs internet access to load CesiumJS from jsDelivr, but the
        orbit data itself is embedded in the HTML file, so no local web server is
        required.
        """
        path = Path(path)
        path.parent.mkdir(parents=True, exist_ok=True)
        title = title or self.name
        czml_json = json.dumps(self.czml).replace("</", "<\\/")
        cesium_root = f"https://cdn.jsdelivr.net/npm/cesium@{cesium_version}/Build/Cesium"

        html = f'''<!DOCTYPE html>
<html lang="en">
<head>
  <meta charset="utf-8" />
  <meta name="viewport" content="width=device-width, initial-scale=1.0" />
  <title>{title}</title>
  <script>window.CESIUM_BASE_URL = "{cesium_root}/";</script>
  <script src="{cesium_root}/Cesium.js"></script>
  <link href="{cesium_root}/Widgets/widgets.css" rel="stylesheet" />
  <style>
    html, body, #cesiumContainer {{ width: 100%; height: 100%; margin: 0; padding: 0; overflow: hidden; }}
  </style>
</head>
<body>
<div id="cesiumContainer"></div>
<script>
(async function () {{
  const viewer = new Cesium.Viewer("cesiumContainer", {{
    animation: true,
    timeline: true,
    baseLayerPicker: false,
    geocoder: false,
    navigationHelpButton: false,
    sceneModePicker: true,
    baseLayer: false,
    terrainProvider: new Cesium.EllipsoidTerrainProvider()
  }});

  try {{
    const naturalEarth = await Cesium.TileMapServiceImageryProvider.fromUrl(
      Cesium.buildModuleUrl("Assets/Textures/NaturalEarthII")
    );
    viewer.imageryLayers.addImageryProvider(naturalEarth);
  }} catch (error) {{
    console.warn("Natural Earth imagery could not be loaded:", error);
  }}

  const czml = {czml_json};
  const dataSource = await Cesium.CzmlDataSource.load(czml);
  viewer.dataSources.add(dataSource);
  await viewer.zoomTo(dataSource);
}})();
</script>
</body>
</html>
'''
        path.write_text(html, encoding="utf-8")
        return path


def times_from_seconds(epoch: datetime, times_s: Iterable[float]) -> list[datetime]:
    """Convert propagation seconds into datetimes for Cesium."""
    epoch = _utc(epoch)
    return [epoch + timedelta(seconds=float(t)) for t in times_s]


__all__ = ["trajectory_packet", "CesiumScene", "times_from_seconds"]
