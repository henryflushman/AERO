"""TLE helpers built on the ``sgp4`` Python package.

This module intentionally focuses on two-line elements and SGP4 state generation.
It does not replace PolySpace's SpaceObject / SpaceTrackJSONParser classes, which
already handle Space-Track JSON records and typed catalog objects.
"""

from __future__ import annotations

from dataclasses import dataclass
from datetime import datetime, timedelta, timezone
from typing import Iterable

import numpy as np

try:
    from sgp4.api import SGP4_ERRORS, Satrec, WGS72, jday
    from sgp4.conveniences import sat_epoch_datetime
except ImportError:
    SGP4_ERRORS = {}
    Satrec = None
    WGS72 = None
    jday = None
    sat_epoch_datetime = None


def _require_sgp4() -> None:
    if Satrec is None:
        raise ImportError(
            "TLETools requires the sgp4 package. Install it with: "
            "python -m pip install sgp4"
        )

try:
    from .OrbitalElements import ClassicalOrbitalElements, StateVector, rv_to_coe
except ImportError:
    from OrbitalElements import ClassicalOrbitalElements, StateVector, rv_to_coe


def _utc(time: datetime) -> datetime:
    if not isinstance(time, datetime):
        raise TypeError("time must be a datetime.")
    if time.tzinfo is None:
        return time.replace(tzinfo=timezone.utc)
    return time.astimezone(timezone.utc)


@dataclass(frozen=True)
class TLERecord:
    """One named or unnamed two-line element set."""

    line1: str
    line2: str
    name: str | None = None

    def __post_init__(self) -> None:
        line1 = self.line1.strip()
        line2 = self.line2.strip()
        if not line1.startswith("1 "):
            raise ValueError("line1 does not look like TLE line 1.")
        if not line2.startswith("2 "):
            raise ValueError("line2 does not look like TLE line 2.")
        object.__setattr__(self, "line1", line1)
        object.__setattr__(self, "line2", line2)
        if self.name is not None:
            object.__setattr__(self, "name", self.name.strip() or None)

    @classmethod
    def from_text(cls, text: str) -> "TLERecord":
        """Parse a 2-line or 3-line TLE text block."""
        lines = [line.strip() for line in text.splitlines() if line.strip()]
        if len(lines) == 2:
            return cls(lines[0], lines[1])
        if len(lines) == 3:
            return cls(lines[1], lines[2], name=lines[0])
        raise ValueError("TLE text must contain either 2 nonblank lines or name + 2 TLE lines.")

    @property
    def satrec(self):
        _require_sgp4()
        return Satrec.twoline2rv(self.line1, self.line2, WGS72)

    @property
    def norad_id(self) -> int:
        return int(self.satrec.satnum)

    @property
    def epoch(self) -> datetime:
        return sat_epoch_datetime(self.satrec).astimezone(timezone.utc)

    def state_at(self, time: datetime) -> StateVector:
        """Return the SGP4 TEME state at ``time`` in km and km/s."""
        time = _utc(time)
        sec = time.second + time.microsecond / 1e6
        jd, fr = jday(time.year, time.month, time.day, time.hour, time.minute, sec)
        error, r_km, v_km_s = self.satrec.sgp4(jd, fr)
        if error != 0:
            message = SGP4_ERRORS.get(error, f"SGP4 error code {error}")
            raise RuntimeError(message)
        return StateVector(np.asarray(r_km), np.asarray(v_km_s))

    def state_at_epoch(self) -> StateVector:
        return self.state_at(self.epoch)

    def coes_at(self, time: datetime) -> ClassicalOrbitalElements:
        """Convert the SGP4 TEME state at a time to osculating two-body COEs."""
        state = self.state_at(time)
        return rv_to_coe(state.r_km, state.v_km_s)

    def coes_at_epoch(self) -> ClassicalOrbitalElements:
        return self.coes_at(self.epoch)

    def propagate(
        self,
        start: datetime,
        stop: datetime,
        step_s: float,
    ) -> tuple[list[datetime], np.ndarray, np.ndarray]:
        """Sample the SGP4 TEME state on a uniform datetime grid."""
        if step_s <= 0.0:
            raise ValueError("step_s must be positive.")
        start = _utc(start)
        stop = _utc(stop)
        if stop < start:
            raise ValueError("stop must be at or after start.")

        total_s = (stop - start).total_seconds()
        count = int(np.floor(total_s / step_s)) + 1
        times = [start + timedelta(seconds=k * step_s) for k in range(count)]
        if times[-1] < stop:
            times.append(stop)

        states = [self.state_at(t) for t in times]
        r = np.vstack([state.r_km for state in states])
        v = np.vstack([state.v_km_s for state in states])
        return times, r, v


def parse_tle_file(path) -> list[TLERecord]:
    """Parse a text file containing repeated 2-line or 3-line TLE records."""
    from pathlib import Path

    lines = [line.strip() for line in Path(path).read_text(encoding="utf-8").splitlines() if line.strip()]
    records: list[TLERecord] = []
    i = 0
    while i < len(lines):
        if lines[i].startswith("1 "):
            if i + 1 >= len(lines):
                raise ValueError("TLE file ends after line 1.")
            records.append(TLERecord(lines[i], lines[i + 1]))
            i += 2
        else:
            if i + 2 >= len(lines):
                raise ValueError("Incomplete named TLE at end of file.")
            records.append(TLERecord(lines[i + 1], lines[i + 2], name=lines[i]))
            i += 3
    return records


__all__ = ["TLERecord", "parse_tle_file"]
