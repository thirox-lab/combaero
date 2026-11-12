"""Data-loading schema for the cooling validation dataset.

Each source folder carries a metadata.yaml describing every CSV. This
module parses those into typed dataclasses and exposes load_dataset().

The digitised coordinates in the CSVs are kept exactly as the digitiser
produced them. Where a figure's axis had to be rescaled to recover a
decade, the factor lives in metadata as `x_scale` rather than being
applied to the file, so the committed numbers stay a faithful record of
the measurement.
"""

from __future__ import annotations

import csv
from dataclasses import dataclass
from pathlib import Path
from typing import Literal

import yaml

DATA_ROOT = Path(__file__).parent / "data"

Kind = Literal["measured", "correlation", "frame"]
Extraction = Literal["tabulated", "figure-digitised"]
Confidence = Literal["exact", "band"]
# How much the SYMBOL-CLASS label is trusted, as distinct from the
# coordinate values. "disputed" means the evidence actively conflicts;
# such a series may be pooled but must never be used per-class.
ClassConfidence = Literal["confirmed", "provisional", "disputed"]


@dataclass(frozen=True)
class SourceMetadata:
    name: str
    citation: str
    secondary: bool


@dataclass(frozen=True)
class SeriesMetadata:
    """One digitised series: a CSV plus the provenance record for it."""

    path: Path
    source: SourceMetadata
    after: str | None  # primary paper, where the source is secondary
    page: int | None
    item: str
    series: str
    figure: str | None  # e.g. "4.46"; groups series for plotting
    panel: str | None  # "upper" / "lower" / "single"
    geometry: dict[str, float] | None
    alpha_deg: float | None
    x_axis: str
    x_scale: float
    # Like x_scale, for an ordinate printed in scaled units (Taslim & Spring
    # plot Nu x 1e-3); the CSV keeps the digitised numbers.
    y_scale: float
    y_axis: str
    kind: Kind
    extraction: Extraction
    confidence: Confidence
    uncertainty: float | None
    cross_check: str
    class_confidence: ClassConfidence
    scores: str | None  # correlation-set name, or None if not scored
    # Rib SHAPE beyond the angle: transverse, parallel, crossed, v, lambda,
    # v_broken. None means not recorded; a set that binds a non-trivial
    # shape refuses such a series rather than assuming one.
    rib_shape: str | None = None
    # The figure card: what the printed axes and equations say, read off
    # the page independently of where the digitiser put the points. See
    # validation/cooling/verify.py.
    verification: dict | None = None
    # Further sets that score this series, each reported on its own row and
    # never ranked against `scores` or each other. Declared, but not chosen:
    # the list must be EVERY eligible set, which a test derives by rule, so
    # an unflattering set cannot be left off (see conventions.eligible_sets).
    also_scores: tuple[str, ...] = ()
    # The rig as the source PRINTS it, where that is less than `geometry`
    # needs: wall count, which walls are ribbed, and null for every parameter
    # the source does not state. Kept apart from `geometry` so the rib runner,
    # which fills a missing e/D or P/e with a default, never sees it; read by
    # runners that score only what the printed rig supports
    # (rib_wall_count_runner).
    rig: dict | None = None

    @property
    def label(self) -> str:
        return f"{self.source.name}/{self.path.stem}"


@dataclass(frozen=True)
class Point:
    x: float  # already multiplied by x_scale
    y: float  # already multiplied by y_scale


def load_points(series: SeriesMetadata) -> list[Point]:
    """Read one CSV, applying the recorded x_scale."""
    points: list[Point] = []
    with open(series.path, newline="") as fh:
        for row in csv.reader(fh):
            if not row or not row[0].strip() or row[0].strip() == "x":
                continue
            points.append(
                Point(float(row[0]) * series.x_scale, float(row[1]) * series.y_scale)
            )
    if not points:
        raise ValueError(f"{series.path} contains no data rows")
    return points


def load_dataset(root: Path | None = None) -> list[SeriesMetadata]:
    """Load every series under the data root, in a stable order."""
    root = root or DATA_ROOT
    out: list[SeriesMetadata] = []
    for meta_path in sorted(root.glob("*/metadata.yaml")):
        raw = yaml.safe_load(meta_path.read_text())
        source = SourceMetadata(
            name=raw["paper"]["name"],
            citation=raw["paper"]["citation"],
            secondary=bool(raw["paper"].get("secondary", False)),
        )
        for entry in raw["files"]:
            path = meta_path.parent / entry["path"]
            if not path.exists():
                raise FileNotFoundError(f"{meta_path} names a missing file: {path}")
            out.append(
                SeriesMetadata(
                    path=path,
                    source=source,
                    after=entry.get("after"),
                    page=entry.get("page"),
                    item=entry["item"],
                    series=entry["series"],
                    figure=entry.get("figure"),
                    panel=entry.get("panel"),
                    geometry=entry.get("geometry"),
                    alpha_deg=entry.get("alpha_deg"),
                    x_axis=entry["x_axis"],
                    x_scale=float(entry.get("x_scale", 1.0)),
                    y_scale=float(entry.get("y_scale", 1.0)),
                    y_axis=entry["y_axis"],
                    kind=entry["kind"],
                    extraction=entry["extraction"],
                    confidence=entry["confidence"],
                    uncertainty=entry.get("uncertainty"),
                    cross_check=entry["cross_check"],
                    class_confidence=entry.get(
                        "class_confidence", "confirmed"
                    ),
                    scores=entry.get("scores"),
                    rib_shape=entry.get("rib_shape"),
                    verification=entry.get("verification"),
                    also_scores=tuple(entry.get("also_scores") or ()),
                    rig=entry.get("rig"),
                )
            )
    return out
