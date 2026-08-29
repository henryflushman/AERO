"""Small Python-facing wrapper around PolySpace's existing Space-Track pipeline.

This module deliberately does *not* reimplement Space-Track storage or query
logic. It calls the existing function-enabled PolySpace files:

- FetchGPHistory.py
- FetchSATCAT.py
- BuildRealDataset.py
- GPDataManager.py
"""

from __future__ import annotations

import sys
from pathlib import Path
from typing import Iterable


class SpaceTrackTools:
    """Convenient access to the existing PolySpace Space-Track functions."""

    def __init__(
        self,
        repo_root: str | Path | None = None,
        *,
        fetching_dir: str | Path | None = None,
        gp_history_cache: str | Path = "data/gp_history",
        satcat_cache: str | Path = "data/satcat",
    ) -> None:
        self.repo_root = Path(repo_root or Path.cwd()).resolve()
        self.fetching_dir = Path(
            fetching_dir or self.repo_root / "DevelopmentScripts" / "spacetrack"
        ).resolve()
        self.gp_history_cache = Path(gp_history_cache)
        self.satcat_cache = Path(satcat_cache)

        if str(self.fetching_dir) not in sys.path:
            sys.path.insert(0, str(self.fetching_dir))

        try:
            from FetchGPHistory import fetch_gp_history
            from FetchSATCAT import fetch_satcat
            from BuildRealDataset import build_real_dataset
            from GPDataManager import GPHistoryCache, SatCatCache, manage_caches
        except ImportError as exc:
            raise ImportError(
                "Could not import the function-enabled Space-Track files. "
                f"Expected them in: {self.fetching_dir}"
            ) from exc

        self._fetch_gp_history = fetch_gp_history
        self._fetch_satcat = fetch_satcat
        self._build_real_dataset = build_real_dataset
        self._manage_caches = manage_caches
        self.GPHistoryCache = GPHistoryCache
        self.SatCatCache = SatCatCache

    def status(self):
        """Return/print current GP_HISTORY and SATCAT cache status."""
        return self._manage_caches(
            gp_history_cache=self.gp_history_cache,
            satcat_cache=self.satcat_cache,
            status=True,
        )

    def fetch_gp_history(self, **kwargs):
        """Call PolySpace ``fetch_gp_history`` using this object's cache path."""
        kwargs.setdefault("gp_history_cache", self.gp_history_cache)
        return self._fetch_gp_history(**kwargs)

    def fetch_satcat(self, *, update_satcat: bool = False, dry_run: bool = False):
        """Create/reuse SATCAT, or force replacement when requested."""
        return self._fetch_satcat(
            satcat_cache=self.satcat_cache,
            update_satcat=update_satcat,
            dry_run=dry_run,
        )

    def build_real_dataset(
        self,
        epoch: str,
        *,
        tolerance_hours: float = 12.0,
        sample_mode: str = "nearest",
        output_dir: str | Path = "data/Real_Satellite_Data",
        output_name: str = "real.parquet",
        object_types: Iterable[str] = ("ALL",),
        **filters,
    ):
        """Build the analysis-ready real population from local caches."""
        return self._build_real_dataset(
            epoch=epoch,
            tolerance_hours=tolerance_hours,
            sample_mode=sample_mode,
            gp_history_cache=self.gp_history_cache,
            satcat_cache=self.satcat_cache,
            output_dir=output_dir,
            output_name=output_name,
            object_types=object_types,
            **filters,
        )


__all__ = ["SpaceTrackTools"]
