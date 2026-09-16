#!/usr/bin/env python
"""Certify centerline mileage and refuse to count the same road twice, including the copy someone digitized backwards.

The county certifies total centerline mileage to the state once a year and the
money that comes back is proportional to the number. Sixty-one segments in the
roads layer exist twice, usually digitized in opposite directions by two editors
several years apart. Nothing about the map looks wrong: the two copies sit on top
of each other, one line hidden under the other, and the total is overstated.

A spatial join on geometry equality does not find them. To ArcGIS, and to Shapely,
and to PostGIS, a LineString is an ORDERED list of positions, so the same road
recorded west to east and east to west is two different geometries. ST_Equals in
PostGIS does answer true for a reversed line, because it compares point sets, but
it costs a pairwise comparison and it is not what the tools people actually run in
an editing session compute. There is no equivalent in a Pro attribute join, and
"Find Identical" on SHAPE compares the geometry as stored.

This computes one direction-independent hash per line, so both copies land in the
same bucket in a single pass, and it measures the survivors geodesically on the
WGS84 ellipsoid. The report gives mileage before, mileage after, and the mileage
it removed, because a certification correction has to be a number somebody can
defend line by line.

    python roadmiles.py roads.geojson
    python roadmiles.py roads.geojson --by MAINT_AUTH
    python roadmiles.py roads.geojson --write-duplicates dupes.geojson --apply
    python roadmiles.py --self-test

Exit codes: 0 no duplicates, 1 duplicates were collapsed, 2 a file could not be
read or written, 64 usage error. The certified total is printed either way; exit 1
means the number changed, not that the run failed.
"""

from __future__ import print_function

import argparse
import hashlib
import json
import math
import os
import sys

# =============================================================================
# CONFIGURATION. Deliberately not flags. Change here, not at the call site.
# =============================================================================

# Decimal places coordinates are rounded to before they are hashed. Six decimal
# degrees is about 11 cm at the equator, which is finer than any road centerline
# was ever digitized and coarse enough to absorb the float noise a reprojection
# leaves behind. Raise it and two copies stop matching; lower it and two roads a
# metre apart start matching.
DEFAULT_PRECISION = 6

# WGS84, which is what GeoJSON coordinates are in by definition (RFC 7946).
WGS84_A = 6378137.0
WGS84_F = 1.0 / 298.257223563
WGS84_B = WGS84_A * (1.0 - WGS84_F)

# The international mile, exactly. The US survey mile is 1609.3472186944 m, and
# Florida's state plane feet are survey feet, so a pipeline that converts through
# them is working in the other one. The gap is 3.2 mm per mile: about a hundredth
# of a mile across a 5,000 mile county road network, which is under the third
# decimal this tool prints. It is set here so the choice is visible rather than
# buried in a literal.
METRES_PER_MILE = 1609.344

# Vincenty's inverse formula iterates. Roads converge in three or four passes;
# only a nearly antipodal pair runs out, and that is a data error, not a road.
MAX_ITERATIONS = 200
CONVERGENCE = 1e-12

# Decimal places in the printed mileage. Three is a thousandth of a mile, about
# 1.6 m, which is the resolution a certification is argued at.
REPORT_DECIMALS = 3

# =============================================================================
# End of CONFIGURATION.
# =============================================================================

# How a duplicate sits against the copy that was kept.
IDENTICAL = "identical"
REVERSED = "reversed"
ROTATED = "rotated loop"

# The bucket a feature falls in when --by names a field it does not carry.
MISSING = "(no value)"


class Part(object):
    """One LineString worth of road: where it is, how long, and what it hashes to.

    A MultiLineString feature contributes one Part per member. That is the unit
    the whole tool counts in, because a road split into two parts by an editor is
    two things on the map and a duplicate of one of them is one duplicate.
    """

    def __init__(self, index, part, coords, digest, length_mm, properties):
        self.index = index            # position of the feature in the input
        self.part = part              # position within a MultiLineString, else 0
        self.coords = coords          # the raw positions, as given
        self.digest = digest          # the direction-independent geometry hash
        self.length_mm = length_mm    # integer millimetres, see certify()
        self.properties = properties
        self.relation = None          # set on a duplicate, against the survivor

    @property
    def label(self):
        if self.part:
            return "feature %d part %d" % (self.index, self.part)
        return "feature %d" % self.index

    def __repr__(self):
        return "Part(%s, %s, %d mm)" % (self.label, self.digest[:8], self.length_mm)


class Report(object):
    """Every number the printed certification is made of."""

    def __init__(self, parts, groups, duplicates, skipped, by_field, by_totals,
                 conflicts, precision):
        self.parts = parts
        self.groups = groups              # digest -> [Part], input order
        self.duplicates = duplicates      # the Parts that were NOT counted
        self.skipped = skipped            # [(index, reason)]
        self.by_field = by_field
        self.by_totals = by_totals        # field value -> millimetres, net
        self.conflicts = conflicts        # [(digest, [differing values])]
        self.precision = precision

    @property
    def collapsed_mm(self):
        return sum(p.length_mm for p in self.duplicates)

    @property
    def net_mm(self):
        """Millimetres after collapsing: one survivor from each group."""
        return sum(members[0].length_mm for members in self.groups.values())

    @property
    def gross_mm(self):
        """Millimetres before collapsing, every part counted.

        Built by addition rather than measured separately, so that
        gross - net == collapsed holds as an identity over integers instead of
        holding to within whatever the float summation order happened to give.
        The whole point of the report is that the correction is defensible, and a
        correction that does not reconcile to the last millimetre is not.
        """
        return self.net_mm + self.collapsed_mm

    @property
    def duplicate_groups(self):
        return [m for m in self.groups.values() if len(m) > 1]


# ----------------------------------------------------------------- pure core

def geodesic_metres(lon1, lat1, lon2, lat2):
    """Distance on the WGS84 ellipsoid, by Vincenty's inverse formula.

    Planar distance over degrees is the other way this number silently comes out
    wrong. A degree of longitude is 111.3 km at the equator and 97.3 km at a central
    Florida latitude, so treating degrees as a flat grid overstates every
    east-west road by about 14 percent here, and by more the further north you
    go. There is no error message, only a bigger number.

    A crossing of the antimeridian needs no special case. The longitude
    difference arrives as -359.998 degrees rather than 0.002, and every use of it
    below goes through sin or cos, which do not care. The naive subtraction that
    would return half the planet is not in here.

    Raises ValueError rather than return a wrong answer for a nearly antipodal
    pair, where the formula does not converge. No road centerline is 20,000 km
    long, so this means a coordinate is wrong or the axes are swapped.
    """
    for name, value in (("lon1", lon1), ("lat1", lat1),
                        ("lon2", lon2), ("lat2", lat2)):
        if value is None or isinstance(value, bool) or not isinstance(
                value, (int, float)):
            raise ValueError("%s must be a number, got %r" % (name, value))
        # NaN and infinity are floats and pass every check above. Left alone, a
        # NaN runs all 200 iterations and comes back accused of being nearly
        # antipodal, and an infinity comes back as math.sin's "math domain
        # error". Both reach here from a converter that turned a null into a
        # float, and the report has to say that rather than guess.
        if not math.isfinite(value):
            raise ValueError("%s is %r, which is not a finite number. A null "
                             "coordinate became a float somewhere upstream"
                             % (name, value))
    for name, value in (("lat1", lat1), ("lat2", lat2)):
        if not -90.0 <= value <= 90.0:
            raise ValueError("%s is %r, outside -90..90. Coordinates are "
                             "longitude first in GeoJSON" % (name, value))

    lam_l = math.radians(lon2 - lon1)
    u1 = math.atan((1 - WGS84_F) * math.tan(math.radians(lat1)))
    u2 = math.atan((1 - WGS84_F) * math.tan(math.radians(lat2)))
    sin_u1, cos_u1 = math.sin(u1), math.cos(u1)
    sin_u2, cos_u2 = math.sin(u2), math.cos(u2)

    lam = lam_l
    for _ in range(MAX_ITERATIONS):
        sin_lam, cos_lam = math.sin(lam), math.cos(lam)
        sin_sigma = math.sqrt(
            (cos_u2 * sin_lam) ** 2
            + (cos_u1 * sin_u2 - sin_u1 * cos_u2 * cos_lam) ** 2)
        if sin_sigma == 0.0:
            return 0.0  # the two positions are the same point
        cos_sigma = sin_u1 * sin_u2 + cos_u1 * cos_u2 * cos_lam
        sigma = math.atan2(sin_sigma, cos_sigma)
        sin_alpha = cos_u1 * cos_u2 * sin_lam / sin_sigma
        cos_sq_alpha = 1 - sin_alpha ** 2
        if cos_sq_alpha == 0.0:
            # An exactly equatorial line has no point of closest approach to the
            # equator to measure from, and the term it feeds drops out.
            cos_2sigma_m = 0.0
        else:
            cos_2sigma_m = cos_sigma - 2 * sin_u1 * sin_u2 / cos_sq_alpha
        c = WGS84_F / 16 * cos_sq_alpha * (4 + WGS84_F * (4 - 3 * cos_sq_alpha))
        previous = lam
        lam = lam_l + (1 - c) * WGS84_F * sin_alpha * (
            sigma + c * sin_sigma * (cos_2sigma_m + c * cos_sigma
                                     * (-1 + 2 * cos_2sigma_m ** 2)))
        if abs(lam - previous) < CONVERGENCE:
            break
    else:
        raise ValueError(
            "the geodesic between (%r, %r) and (%r, %r) did not converge in %d "
            "iterations. Those two positions are nearly antipodal, which no road "
            "centerline is." % (lon1, lat1, lon2, lat2, MAX_ITERATIONS))

    u_sq = cos_sq_alpha * (WGS84_A ** 2 - WGS84_B ** 2) / WGS84_B ** 2
    big_a = 1 + u_sq / 16384 * (4096 + u_sq * (-768 + u_sq * (320 - 175 * u_sq)))
    big_b = u_sq / 1024 * (256 + u_sq * (-128 + u_sq * (74 - 47 * u_sq)))
    d_sigma = big_b * sin_sigma * (
        cos_2sigma_m + big_b / 4 * (
            cos_sigma * (-1 + 2 * cos_2sigma_m ** 2)
            - big_b / 6 * cos_2sigma_m * (-3 + 4 * sin_sigma ** 2)
            * (-3 + 4 * cos_2sigma_m ** 2)))
    return WGS84_B * big_a * (sigma - d_sigma)


def line_metres(coords):
    """Geodesic length of a list of positions, leg by leg."""
    total = 0.0
    for i in range(len(coords) - 1):
        (lon1, lat1), (lon2, lat2) = coords[i][:2], coords[i + 1][:2]
        total += geodesic_metres(lon1, lat1, lon2, lat2)
    return total


def miles(millimetres):
    """Millimetres as miles. One function, so every mileage agrees with itself."""
    return millimetres / 1000.0 / METRES_PER_MILE


def rounded_points(coords, precision=DEFAULT_PRECISION):
    """The positions at the stated precision, with repeats squeezed out.

    Two things happen here and both matter to the hash. Rounding is what lets a
    reprojected copy match the original. Dropping a vertex that is equal to the
    one before it, AFTER rounding, is what stops the same road matching or not
    matching depending on whether an editor double-clicked.

    The + 0.0 is not decoration. round(-0.0000001, 6) is -0.0, which compares
    equal to 0.0 and formats as "-0.000000", so a copy that approached the prime
    meridian from the other side hashed differently while every comparison in
    this file said the two were the same.
    """
    if precision < 0:
        raise ValueError("--precision cannot be negative")
    points = []
    for position in coords:
        point = (round(position[0], precision) + 0.0,
                 round(position[1], precision) + 0.0)
        if not points or point != points[-1]:
            points.append(point)
    return points


def canonical_points(coords, precision=DEFAULT_PRECISION):
    """One vertex list per road, whichever way round it was digitized.

    This is the whole idea and it is four lines of decision.

    An open line has exactly two readings, forwards and backwards, so the smaller
    of the two is a canonical name for the pair.

    A closed loop has two readings for every vertex it could have started at,
    because the start vertex of a ring is an arbitrary choice the editor made.
    So the candidates are every rotation of the ring in both directions, and
    again the smallest is the name. That is O(n^2) in the vertex count of a loop.
    A cul-de-sac or a roundabout has tens of vertices, so it does not matter; a
    ring with ten thousand vertices is not a road.
    """
    points = rounded_points(coords, precision)
    if len(points) >= 3 and points[0] == points[-1]:
        ring = points[:-1]
        best = min(tuple(seq[i:] + seq[:i])
                   for seq in (ring, ring[::-1])
                   for i in range(len(ring)))
        return list(best) + [best[0]]
    return list(min(tuple(points), tuple(reversed(points))))


def geometry_hash(coords, precision=DEFAULT_PRECISION):
    """A hex digest that is equal for the same road digitized either way.

    The digest is taken over the formatted decimal text, not over Python's repr
    of the floats, so the same coordinates give the same digest on every
    interpreter and the hashes in last year's report still mean something.
    """
    canon = canonical_points(coords, precision)
    blob = ";".join("%.*f,%.*f" % (precision, x, precision, y)
                    for x, y in canon)
    return hashlib.sha256(blob.encode("ascii")).hexdigest()


def duplicate_relation(kept, found, precision=DEFAULT_PRECISION):
    """How a duplicate sits against the copy that was kept.

    Worth saying out loud in the report. "Reversed" tells the editor these are
    two digitizing sessions, not one feature copied and pasted, which is a
    different conversation about how it got there.
    """
    a = rounded_points(kept, precision)
    b = rounded_points(found, precision)
    if a == b:
        return IDENTICAL
    if a == b[::-1]:
        return REVERSED
    return ROTATED


def geometry_parts(geometry):
    """Every LineString in a geometry, as a list of position lists.

    Raises ValueError with the reason when the geometry is not a centerline, so
    the caller can count it as skipped and name it instead of guessing.
    """
    if not isinstance(geometry, dict):
        raise ValueError("no geometry")
    kind = geometry.get("type")
    coords = geometry.get("coordinates")
    if kind == "LineString":
        lines = [coords]
    elif kind == "MultiLineString":
        if not isinstance(coords, list):
            raise ValueError("MultiLineString has no list of parts")
        lines = coords
    else:
        raise ValueError("geometry type %s is not a centerline" % (kind,))

    for line in lines:
        if not isinstance(line, list) or len(line) < 2:
            raise ValueError("a part has fewer than 2 positions")
        for position in line:
            if not isinstance(position, (list, tuple)) or len(position) < 2:
                raise ValueError("a position is not a coordinate pair")
            for value in position[:2]:
                if isinstance(value, bool) or not isinstance(value, (int, float)):
                    raise ValueError("a coordinate is not a number")
    return lines


def certify(features, precision=DEFAULT_PRECISION, by=None):
    """Measure the features, collapse the duplicates, and total what is left.

    Pure: dicts in, a Report out. No file, no network, no arcpy, which is why the
    self-test can certify a duplicated road network without one existing.

    Lengths are carried as integer millimetres. Floats would make the headline
    subtraction depend on the order the parts were summed in, and this number is
    handed to the state.
    """
    if precision < 0:
        raise ValueError("--precision cannot be negative")
    if not isinstance(features, list):
        raise ValueError("features must be a list, got %r" % type(features))

    parts = []
    skipped = []
    for index, feature in enumerate(features):
        if not isinstance(feature, dict):
            skipped.append((index, "not a GeoJSON object"))
            continue
        properties = feature.get("properties")
        if not isinstance(properties, dict):
            properties = {}
        try:
            lines = geometry_parts(feature.get("geometry"))
        except ValueError as exc:
            skipped.append((index, str(exc)))
            continue
        try:
            measured = [(line, line_metres(line)) for line in lines]
        except ValueError as exc:
            skipped.append((index, str(exc)))
            continue
        for number, (line, metres) in enumerate(measured):
            parts.append(Part(
                index, number if len(lines) > 1 else 0, line,
                geometry_hash(line, precision),
                int(round(metres * 1000.0)), properties))

    groups = {}
    duplicates = []
    for part in parts:
        members = groups.setdefault(part.digest, [])
        members.append(part)
        if len(members) > 1:
            part.relation = duplicate_relation(
                members[0].coords, part.coords, precision)
            duplicates.append(part)

    by_totals = {}
    conflicts = []
    if by:
        for digest in groups:
            members = groups[digest]
            survivor = members[0]
            value = survivor.properties.get(by, MISSING)
            key = MISSING if value is None else str(value)
            by_totals[key] = by_totals.get(key, 0) + survivor.length_mm
            # Two copies of one road filed under two different authorities is a
            # second defect sitting inside the first, and collapsing the geometry
            # hides it. Say so rather than silently keeping whichever came first.
            seen = []
            for member in members:
                other = member.properties.get(by, MISSING)
                other = MISSING if other is None else str(other)
                if other not in seen:
                    seen.append(other)
            if len(seen) > 1:
                conflicts.append((digest, seen))

    return Report(parts, groups, duplicates, skipped, by, by_totals, conflicts,
                  precision)


def describe(report):
    """Render a Report as the lines the CLI prints."""
    fmt = "%%.%df mi" % REPORT_DECIMALS
    measured = len(set(part.index for part in report.parts))
    out = ["%d part(s) from %d feature(s), %d skipped"
           % (len(report.parts), measured + len(report.skipped),
              len(report.skipped)),
           "rounded to %d decimal degrees, measured geodesically on WGS84"
           % report.precision,
           ""]
    for index, reason in report.skipped:
        out.append("  SKIPPED feature %d: %s" % (index, reason))
    if report.skipped:
        out.append("")

    out.append("  mileage before, every part counted   %s"
               % (fmt % miles(report.gross_mm)))
    out.append("  duplicate groups collapsed           %d group(s), %d part(s)"
               % (len(report.duplicate_groups), len(report.duplicates)))
    out.append("  mileage removed as duplicate         %s"
               % (fmt % miles(report.gross_mm - report.net_mm)))
    out.append("  " + "-" * 52)
    out.append("  CERTIFIED mileage                    %s"
               % (fmt % miles(report.net_mm)))

    if report.duplicate_groups:
        out.append("")
        out.append("duplicates, in input order:")
        for members in report.duplicate_groups:
            out.append("  %s  %d copies, %s"
                       % (members[0].digest[:12], len(members),
                          fmt % miles(members[0].length_mm)))
            out.append("      kept    %s" % members[0].label)
            for member in members[1:]:
                out.append("      removed %s, %s"
                           % (member.label, member.relation))

    if report.by_field:
        out.append("")
        out.append("certified mileage by %s:" % report.by_field)
        for key in sorted(report.by_totals):
            out.append("  %-28s %s" % (key, fmt % miles(report.by_totals[key])))
        for digest, values in report.conflicts:
            out.append("  CONFLICT %s: copies of one road disagree about %s: %s"
                       % (digest[:12], report.by_field, ", ".join(values)))

    out.append("VERDICT: %s" % ("COLLAPSED" if report.duplicates else "NO DUPLICATES"))
    return out


def duplicate_features(report):
    """The removed parts as GeoJSON features, for the editor who has to fix them.

    A duplicate part of a MultiLineString is written out as the part, not as the
    whole feature it came from. Writing the parent would hand back geometry that
    is not duplicated and is not to be deleted.

    The properties are copied, never the original dict, so the input the caller
    still holds is not quietly annotated.
    """
    out = []
    for part in report.duplicates:
        kept = report.groups[part.digest][0]
        properties = dict(part.properties)
        properties["roadmiles_group"] = part.digest[:12]
        properties["roadmiles_relation"] = part.relation
        properties["roadmiles_duplicate_of"] = kept.label
        properties["roadmiles_miles"] = round(miles(part.length_mm),
                                              REPORT_DECIMALS)
        out.append({"type": "Feature",
                    "properties": properties,
                    "geometry": {"type": "LineString",
                                 "coordinates": part.coords}})
    return out


# ------------------------------------------------------------------------ io

def read_features(path):
    """The feature list from a GeoJSON file, or from stdin when path is "-".

    Accepts a FeatureCollection, a bare Feature, a bare geometry, or a plain list
    of features, because all four turn up in an export and refusing three of them
    helps nobody.
    """
    if path == "-":
        text = sys.stdin.read()
    else:
        with open(path, "rb") as handle:
            text = handle.read().decode("utf-8-sig")
    data = json.loads(text)
    if isinstance(data, list):
        return data
    if not isinstance(data, dict):
        raise ValueError("this is not GeoJSON: the top level is a %s"
                         % type(data).__name__)
    kind = data.get("type")
    if kind == "FeatureCollection":
        features = data.get("features")
        if not isinstance(features, list):
            raise ValueError("FeatureCollection has no features array")
        return features
    if kind == "Feature":
        return [data]
    if kind in ("LineString", "MultiLineString"):
        return [{"type": "Feature", "properties": {}, "geometry": data}]
    raise ValueError("this is not a road layer: the top level type is %r"
                     % (kind,))


def write_duplicates(path, features, apply_it):
    """Write the duplicate parts as a GeoJSON FeatureCollection.

    Written as bytes with explicit newlines. Text mode turns every newline into
    CRLF on Windows, and the same command then produces two different files
    depending on which desk it was run from, which is a poor look on a file that
    exists to settle an argument.
    """
    if os.path.exists(path):
        print("error: %s already exists. Refusing to overwrite it." % path,
              file=sys.stderr)
        return 2
    if not apply_it:
        print("Would write %d duplicate part(s) to %s. Re-run with --apply."
              % (len(features), path))
        return 0
    payload = json.dumps({"type": "FeatureCollection", "features": features},
                         indent=1, sort_keys=True)
    try:
        with open(path, "wb") as handle:
            handle.write(payload.encode("utf-8"))
            handle.write(b"\n")
    except (OSError, IOError) as exc:
        print("error: could not write %s: %s" % (path, exc), file=sys.stderr)
        return 2
    print("Wrote %d duplicate part(s) to %s" % (len(features), path))
    return 0


# ------------------------------------------------------------------ self-test

def self_test():
    """Assertions over the decision core, then over the io layer.

    Nothing here reaches the network, a database or arcpy. The io half writes
    GeoJSON into a temporary directory, because a report that has never been
    produced from a real file has not been tested.
    """
    import io
    import shutil
    import tempfile

    passed = [0]
    failed = []

    def check(cond, label):
        if cond:
            passed[0] += 1
            print("PASS  %s" % label)
        else:
            failed.append(label)
            print("FAIL  %s" % label)

    def raises(fn, label):
        try:
            fn()
        except ValueError:
            check(True, label)
        except Exception as exc:
            check(False, "%s (wrong exception %r)" % (label, exc))
        else:
            check(False, "%s (no error raised)" % label)

    def feature(coords, kind="LineString", **props):
        return {"type": "Feature", "properties": dict(props),
                "geometry": {"type": kind, "coordinates": coords}}

    # A real road: a mile and a bit of SR-40 east of Riverton, west to east.
    SR40 = [[-82.10000, 29.18700], [-82.09000, 29.18710],
            [-82.08000, 29.18730], [-82.07000, 29.18700]]
    SR40_BACKWARDS = list(reversed(SR40))

    print("roadmiles self-test: no network, no database, no arcpy")
    print("-" * 68)

    # ---- the geodesic against the published WGS84 meridian arc
    # The meridian arc is the one length anyone can look up, so it is the one to
    # check the implementation against. These three are the standard figures.
    check(abs(geodesic_metres(0, 0, 0, 1) - 110574.4) < 1.0,
          "one degree of latitude at the equator is 110574.4 m, within a metre")
    check(abs(geodesic_metres(0, 44, 0, 45) - 111122.0) < 1.0,
          "the 44 to 45 degree meridian arc is 111122.0 m, within a metre")
    check(abs(geodesic_metres(0, 89, 0, 90) - 111693.9) < 1.0,
          "the 89 to 90 degree meridian arc is 111693.9 m, within a metre")
    check(abs(geodesic_metres(0, 0, 1, 0) - 111319.491) < 0.01,
          "one degree of longitude at the equator is the sphere's own 111319.491 m")
    check(geodesic_metres(0, 1, 0, 0) == geodesic_metres(0, 0, 0, 1),
          "the distance is the same measured either way round")
    check(abs(geodesic_metres(-82.0, 29.2, -81.0, 29.2) - 97250.45) < 0.01,
          "one degree of longitude at a central Florida latitude is 97250.45 m")

    # ---- planar on degrees, which is how this number goes wrong the other way
    naive = math.hypot(1.0, 0.0) * 111320.0
    real = geodesic_metres(-82.0, 29.2, -81.0, 29.2)
    check(naive - real > 14000.0,
          "planar degrees overstate that same east-west mile by over 14 km, "
          "14 percent, with no error message")
    check(abs(geodesic_metres(-82.0, 29.2, -82.0, 30.2) - 110847.4) < 1.0,
          "the same one degree north-south is 110847.4 m, so the planar error "
          "is not a constant you could divide out")

    # ---- the geodesic on the small distances a centerline is made of
    check(geodesic_metres(-82.0, 29.2, -82.0, 29.2) == 0.0,
          "a zero-length leg measures exactly zero")
    check(abs(geodesic_metres(-82.0, 29.2, -82.00001, 29.2) - 0.9725) < 0.001,
          "a one hundred-thousandth of a degree leg is 97 cm")
    check(geodesic_metres(-82.0, 29.2, -81.999999, 29.2) > 0.0,
          "a leg below the hashing precision still has a positive length")
    check(geodesic_metres(0, -33.9, 151.2, -33.9) > 0.0,
          "a long east-west geodesic is positive, the formula runs past one leg")

    # ---- THE ANTIMERIDIAN: the longitude difference arrives as -359.998
    check(abs(geodesic_metres(179.999, 0.0, -179.999, 0.0) - 222.64) < 0.01,
          "a segment crossing the antimeridian is 222 m, not half the "
          "planet  <-- pinned defect")
    check(abs(geodesic_metres(-179.999, 0.0, 179.999, 0.0) - 222.64) < 0.01,
          "and the same crossing measured the other way is the same 222 m")
    check(abs(geodesic_metres(179.99, 65.0, -179.99, 65.0) - 940.0) < 5.0,
          "an antimeridian crossing away from the equator is short too")
    check(line_metres([[179.999, 0.0], [-179.999, 0.0]])
          < line_metres([[0.0, 0.0], [0.01, 0.0]]) * 0.3,
          "so an antimeridian road is shorter than a hundredth of a degree of "
          "equator, not 180,000 times longer  <-- pinned defect")

    # ---- the geodesic refuses rather than guess
    raises(lambda: geodesic_metres(0.0, 0.0, 179.7, 0.5),
           "a nearly antipodal pair raises instead of returning a wrong number")
    raises(lambda: geodesic_metres(-82.0, 91.0, -82.0, 29.2),
           "a latitude over 90 raises, which is what swapped axes look like")
    raises(lambda: geodesic_metres(-82.0, 29.2, -82.0, -90.5),
           "a latitude under -90 raises too")
    raises(lambda: geodesic_metres("-82.0", 29.2, -82.0, 29.3),
           "a coordinate that is a string raises")
    raises(lambda: geodesic_metres(None, 29.2, -82.0, 29.3),
           "a null coordinate raises")
    raises(lambda: geodesic_metres(True, 29.2, -82.0, 29.3),
           "a boolean is not a coordinate, though python says it is a number")
    raises(lambda: geodesic_metres(float("nan"), 29.2, -82.0, 29.3),
           "a NaN coordinate raises")
    raises(lambda: geodesic_metres(float("inf"), 29.2, -82.0, 29.3),
           "an infinite coordinate raises before math.sin turns it into a "
           "domain error")
    nan_reason = ""
    try:
        geodesic_metres(float("nan"), 29.2, -82.0, 29.3)
    except ValueError as exc:
        nan_reason = str(exc)
    check("finite" in nan_reason and "antipodal" not in nan_reason,
          "and a NaN is named as not finite, not accused of being nearly "
          "antipodal  <-- pinned defect")

    # ---- adding the legs up
    check(line_metres([[-82.0, 29.2]]) == 0.0,
          "a single position has no legs and no length")
    check(abs(line_metres([[-82.0, 29.2], [-81.0, 29.2], [-80.0, 29.2]])
              - 2 * 97250.45) < 0.1,
          "two equal legs are twice one leg")
    check(line_metres(SR40) > 0.0, "a real centerline has a positive length")
    check(abs(line_metres(SR40) - line_metres(SR40_BACKWARDS)) < 1e-6,
          "and measures the same digitized backwards")
    check(line_metres([[-82.0, 29.2, 12.5], [-81.0, 29.2, 30.0]])
          == line_metres([[-82.0, 29.2], [-81.0, 29.2]]),
          "a z value is carried in the data and ignored in the length")

    # ---- rounding, which happens before anything is hashed
    check(rounded_points([[-82.123456789, 29.2]])[0][0] == -82.123457,
          "a position is rounded to the stated precision")
    check(rounded_points([[-82.1, 29.2], [-82.1, 29.2]]) == [(-82.1, 29.2)],
          "a vertex repeated exactly is squeezed out")
    check(rounded_points([[-82.1, 29.2], [-82.1000000001, 29.2]])
          == [(-82.1, 29.2)],
          "a vertex that rounds onto the one before it is squeezed out too")
    check(len(rounded_points([[-82.1, 29.2], [-82.2, 29.2], [-82.1, 29.2]])) == 3,
          "a vertex that repeats but not consecutively is kept, a loop needs it")
    check(rounded_points([[-0.0000001, 29.2]])[0][0] == 0.0,
          "a position just west of the prime meridian rounds to positive zero")
    check("%.6f" % rounded_points([[-0.0000001, 29.2]])[0][0] == "0.000000",
          "and formats without a minus sign, so it hashes with its "
          "neighbour  <-- pinned defect")
    check(geometry_hash([[-0.0000001, 29.2], [0.001, 29.2]])
          == geometry_hash([[0.0000001, 29.2], [0.001, 29.2]]),
          "so the two sides of the prime meridian give one hash  <-- pinned defect")
    check(rounded_points([[-82.1, 29.2]], precision=0) == [(-82.0, 29.0)],
          "precision 0 rounds to whole degrees")
    raises(lambda: rounded_points([[-82.1, 29.2]], precision=-1),
           "a negative precision raises")

    # ---- THE HEADLINE: a reversed line is the same road
    check(geometry_hash(SR40) == geometry_hash(SR40_BACKWARDS),
          "a coordinate list and its reverse hash EQUAL  <-- pinned defect")
    check(canonical_points(SR40) == canonical_points(SR40_BACKWARDS),
          "because both canonicalise to the same vertex list")
    check(SR40 != SR40_BACKWARDS,
          "and the two coordinate lists really are different, which is why a "
          "join on geometry equality misses them")
    check(len(geometry_hash(SR40)) == 64, "the hash is a sha256 hex digest")
    check(geometry_hash(SR40) != geometry_hash(
        [[-82.10000, 29.18700], [-82.07000, 29.18700]]),
        "a different road hashes differently")
    check(geometry_hash([[0.0, 0.0], [1.0, 1.0]])
          != geometry_hash([[0.0, 1.0], [1.0, 0.0]]),
          "two lines over the same bounding box are not the same road")
    check(geometry_hash([[0.0, 0.0], [1.0, 1.0], [2.0, 0.0]])
          != geometry_hash([[0.0, 0.0], [2.0, 0.0], [1.0, 1.0]]),
          "the same three vertices visited in another order are another road")
    check(geometry_hash(SR40) == geometry_hash(
        tuple(tuple(p) for p in SR40)),
        "the hash does not care whether positions are lists or tuples")
    check(geometry_hash([[-82.1, 29.2, 12.5], [-82.0, 29.2, 30.0]])
          == geometry_hash([[-82.1, 29.2], [-82.0, 29.2]]),
          "a z value does not change the hash, the same road surveyed twice "
          "carries two elevations")
    check(geometry_hash([[-82.1, 29.2], [-82.0, 29.2]]) == hashlib.sha256(
        "-82.100000,29.200000;-82.000000,29.200000".encode("ascii")).hexdigest(),
        "the digest is sha256 over the formatted decimal text, so a hash "
        "written into last year's report still means something")

    # ---- the rounding precision, either side of the line
    over = [[-82.1, 29.2], [-82.0, 29.2000006]]
    under = [[-82.1, 29.2], [-82.0, 29.2000004]]
    base = [[-82.1, 29.2], [-82.0, 29.2]]
    check(geometry_hash(base) != geometry_hash(over),
          "a vertex nudged PAST the rounding precision hashes differently")
    check(geometry_hash(base) == geometry_hash(under),
          "a vertex nudged UNDER the rounding precision hashes the same")
    check(geometry_hash(base, precision=4) == geometry_hash(over, precision=4),
          "a coarser --precision brings the nudged vertex back into the group")
    check(geometry_hash(base, precision=9) != geometry_hash(under, precision=9),
          "a finer --precision splits a group the default would have collapsed")
    check(geometry_hash(base, precision=0) == geometry_hash(
        [[-82.4, 29.4], [-82.0, 29.2]], precision=0),
        "precision 0 collapses two roads 40 km apart, which is why it is a "
        "knob and not a default")

    # ---- THE PINNED DEFECT: a loop is the same road from any vertex
    # A cul-de-sac digitized by two editors starts at whichever vertex each of
    # them clicked first. The rotation is the only difference and it is not one.
    LOOP = [[-82.10, 29.20], [-82.09, 29.20], [-82.09, 29.21],
            [-82.10, 29.21], [-82.10, 29.20]]
    LOOP_FROM_THIRD = [[-82.09, 29.21], [-82.10, 29.21], [-82.10, 29.20],
                       [-82.09, 29.20], [-82.09, 29.21]]
    LOOP_BACKWARDS = list(reversed(LOOP))
    LOOP_BOTH = list(reversed(LOOP_FROM_THIRD))
    check(geometry_hash(LOOP) == geometry_hash(LOOP_FROM_THIRD),
          "a closed loop starting at a different vertex is the SAME "
          "road  <-- pinned defect")
    check(geometry_hash(LOOP) == geometry_hash(LOOP_BACKWARDS),
          "a closed loop digitized the other way round is the same road")
    check(geometry_hash(LOOP) == geometry_hash(LOOP_BOTH),
          "rotated AND reversed is still the same road  <-- pinned defect")
    check(LOOP != LOOP_FROM_THIRD and LOOP != LOOP_BOTH,
          "and none of those three coordinate lists are equal to each other")
    check(canonical_points(LOOP)[0] == canonical_points(LOOP)[-1],
          "the canonical loop is still closed, so it is still a valid line")
    check(len(canonical_points(LOOP)) == 5,
          "and still has its five positions, no vertex was eaten by the rotation")
    check(canonical_points(LOOP) == [(-82.1, 29.2), (-82.1, 29.21),
                                     (-82.09, 29.21), (-82.09, 29.2),
                                     (-82.1, 29.2)],
          "a ring canonicalises to the SMALLEST rotation, so the ring digest "
          "is stable and not merely self-consistent")
    check(geometry_hash(LOOP) == hashlib.sha256(
        ("-82.100000,29.200000;-82.100000,29.210000;-82.090000,29.210000;"
         "-82.090000,29.200000;-82.100000,29.200000").encode("ascii")
    ).hexdigest(),
        "and last year's report can still be matched against this year's "
        "loops  <-- pinned defect")
    check(geometry_hash(LOOP) != geometry_hash(
        [[-82.10, 29.20], [-82.09, 29.20], [-82.09, 29.22],
         [-82.10, 29.21], [-82.10, 29.20]]),
        "a different loop over four corners is a different road")
    check(geometry_hash(LOOP) != geometry_hash(LOOP[:-1]),
          "an unclosed version of the loop is NOT the same road, one leg of it "
          "is missing")
    check(abs(line_metres(LOOP) - line_metres(LOOP_FROM_THIRD)) < 1e-6,
          "and the rotated loop measures the same length")
    # A three-position line whose ends meet is a degenerate loop. It must take
    # the ring branch without falling off the end of a one-vertex ring.
    check(geometry_hash([[0.0, 0.0], [1.0, 0.0], [0.0, 0.0]])
          == geometry_hash([[1.0, 0.0], [0.0, 0.0], [1.0, 0.0]]),
          "a two-vertex there-and-back loop rotates without falling over")

    # ---- the degenerate shapes a real layer carries
    check(len(canonical_points([[-82.1, 29.2], [-82.0, 29.2]])) == 2,
          "a two-point segment canonicalises to two positions, it is not "
          "mistaken for a ring")
    check(geometry_hash([[-82.1, 29.2], [-82.0, 29.2]])
          == geometry_hash([[-82.0, 29.2], [-82.1, 29.2]]),
          "and a two-point segment is direction-independent like any other")
    zero = [[-82.1, 29.2], [-82.1, 29.2]]
    check(line_metres(zero) == 0.0, "a zero-length segment measures zero")
    check(geometry_hash(zero) == hashlib.sha256(
        "-82.100000,29.200000".encode("ascii")).hexdigest(),
        "a zero-length segment hashes as the ONE position it squeezes down to")
    check(geometry_hash(zero) == geometry_hash([[-82.1, 29.2], [-82.1, 29.2],
                                                [-82.1, 29.2]]),
          "and a zero-length segment of three coincident vertices is the same one")
    check(geometry_hash(zero) != geometry_hash([[-82.0, 29.2], [-82.0, 29.2]]),
          "two zero-length segments in different places are not the same")

    # ---- how a duplicate sits against the copy that was kept
    check(duplicate_relation(SR40, list(SR40)) == IDENTICAL,
          "an exact copy is reported as identical")
    check(duplicate_relation(SR40, SR40_BACKWARDS) == REVERSED,
          "a backwards copy is reported as reversed, which names how it got there")
    check(duplicate_relation(LOOP, LOOP_FROM_THIRD) == ROTATED,
          "a loop from another start vertex is reported as a rotated loop")
    check(duplicate_relation(SR40, [[c[0] + 1e-9, c[1]] for c in SR40])
          == IDENTICAL,
          "a copy differing under the precision is identical, not reversed")

    # ---- reading geometry out of a feature
    check(geometry_parts({"type": "LineString", "coordinates": SR40}) == [SR40],
          "a LineString is one part")
    multi = {"type": "MultiLineString", "coordinates": [SR40, LOOP]}
    check(geometry_parts(multi) == [SR40, LOOP],
          "a MultiLineString is one part per member")
    check(len(geometry_parts({"type": "MultiLineString",
                              "coordinates": [SR40]})) == 1,
          "a MultiLineString of one member is one part")
    raises(lambda: geometry_parts(None), "a feature with no geometry raises")
    raises(lambda: geometry_parts({"type": "Polygon", "coordinates": [SR40]}),
           "a Polygon is not a centerline and raises")
    raises(lambda: geometry_parts({"type": "Point", "coordinates": [-82.0, 29.2]}),
           "a Point is not a centerline and raises")
    raises(lambda: geometry_parts({"type": "MultiLineString",
                                   "coordinates": "nope"}),
           "a MultiLineString with no list of parts raises")
    raises(lambda: geometry_parts({"type": "LineString",
                                   "coordinates": [[-82.0, 29.2]]}),
           "a one-position LineString raises, it is not a line")
    raises(lambda: geometry_parts({"type": "LineString",
                                   "coordinates": [[-82.0], [-81.0]]}),
           "a position with one number raises")
    raises(lambda: geometry_parts({"type": "LineString",
                                   "coordinates": [["a", "b"], [-81.0, 29.0]]}),
           "a position of strings raises")
    raises(lambda: geometry_parts({"type": "LineString",
                                   "coordinates": [None, [-81.0, 29.0]]}),
           "a null position raises")

    # ---- certify: the disaster, in miniature
    # SR-40 entered twice, the second time backwards, plus one unrelated road.
    OTHER = [[-82.20, 29.30], [-82.19, 29.30]]
    layer = [feature(SR40, ROADNAME="SR-40", MAINT_AUTH="STATE"),
             feature(OTHER, ROADNAME="NE 35TH ST", MAINT_AUTH="COUNTY"),
             feature(SR40_BACKWARDS, ROADNAME="SR-40", MAINT_AUTH="STATE")]
    r = certify(layer)
    check(len(r.parts) == 3, "three features give three parts")
    check(len(r.groups) == 2, "which fall into two groups of road")
    check(len(r.duplicates) == 1, "one part is a duplicate")
    check(len(r.duplicate_groups) == 1, "in one duplicate group")
    check(r.duplicates[0].index == 2,
          "the LATER copy is the one removed, the first stays")
    check(r.duplicates[0].relation == REVERSED,
          "and it is recorded as the backwards one  <-- pinned defect")
    check(r.groups[r.duplicates[0].digest][0].index == 0,
          "the survivor is the first copy in input order")
    check(r.skipped == [], "nothing was skipped")

    # ---- THE HEADLINE NUMBER: before, after, and the difference
    check(r.gross_mm == sum(p.length_mm for p in r.parts),
          "the gross is every part measured, each counted once")
    check(r.gross_mm - r.net_mm == r.collapsed_mm,
          "so the before/after difference equals the collapsed duplicate "
          "mileage EXACTLY  <-- pinned defect")
    check(miles(r.gross_mm - r.net_mm) == miles(r.collapsed_mm),
          "and it is exactly equal in miles too, not to within a rounding")
    check(r.collapsed_mm == r.parts[0].length_mm,
          "the mileage removed is exactly the length of the copy removed")
    check(r.net_mm == r.parts[0].length_mm + r.parts[1].length_mm,
          "and the certified total is the two distinct roads, each once")
    check(r.gross_mm > r.net_mm > 0, "the correction reduces the total")
    check(isinstance(r.gross_mm, int) and isinstance(r.net_mm, int),
          "every total is an integer count of millimetres, so no summation "
          "order can move it  <-- pinned defect")
    check(certify([feature([[0.0, 0.0], [0.0, 1.0]])]).parts[0].length_mm
          == 110574389,
          "a degree of latitude certifies as 110574389 mm, the meridian arc "
          "ROUNDED to the millimetre and not truncated  <-- pinned defect")
    check(certify([feature([[0.0, 44.0], [0.0, 45.0]])]).parts[0].length_mm
          == 111122008,
          "and the 44 to 45 arc as 111122008 mm, which rounds DOWN, so no "
          "millimetre is gained per part either  <-- pinned defect")
    check(abs(miles(r.net_mm) - 2.417) < 0.001,
          "the certified total for this layer is 2.417 mi")
    check(abs(miles(r.gross_mm) - 4.231) < 0.001,
          "against 4.231 mi before the duplicate was found")

    # ---- three copies collapse to one, not to two
    triple = certify([feature(SR40), feature(SR40_BACKWARDS), feature(SR40)])
    check(len(triple.groups) == 1, "three copies of one road are one group")
    check(len(triple.duplicates) == 2,
          "and TWO of them are removed, not one  <-- pinned defect")
    check(triple.net_mm == triple.parts[0].length_mm,
          "so the certified mileage is one copy of the road")
    check(triple.gross_mm == 3 * triple.parts[0].length_mm,
          "and the gross was three")
    check([p.relation for p in triple.duplicates] == [REVERSED, IDENTICAL],
          "each duplicate is described against the copy that was kept")

    # ---- two copies of one road that do not measure quite the same
    # A vertex nudged UNDER the rounding precision lands in the group and still
    # changes the length, so WHICH member the certified total is taken from is a
    # decision. Every fixture above has copies of identical length, which would
    # let "keep the last one" pass unnoticed.
    NUDGED = [[-82.10000, 29.18700], [-82.09000, 29.18710],
              [-82.08000, 29.18730], [-82.0700004, 29.18700]]
    nudge = [feature(SR40, MAINT_AUTH="STATE"),
             feature(NUDGED, MAINT_AUTH="STATE")]
    n = certify(nudge)
    check(len(n.groups) == 1,
          "a vertex nudged under the precision is still the same road")
    check(n.parts[0].length_mm != n.parts[1].length_mm,
          "and yet the two copies do not measure the same to the millimetre")
    check(n.net_mm == n.parts[0].length_mm,
          "the certified total is the FIRST copy's length, not the last "
          "one's  <-- pinned defect")
    check(n.collapsed_mm == n.parts[1].length_mm,
          "and the mileage removed is the length of the copy actually removed")
    check(certify(nudge, by="MAINT_AUTH").by_totals["STATE"]
          == n.parts[0].length_mm,
          "and --by buckets that same first copy, so the two totals agree")

    # ---- a MultiLineString
    mixed = [feature([SR40, OTHER], kind="MultiLineString"),
             feature(SR40)]
    m = certify(mixed)
    check(len(m.parts) == 3, "a two-part MultiLineString gives two parts")
    check(m.parts[0].part == 0 and m.parts[1].part == 1,
          "whose part numbers are recorded")
    check(m.parts[2].part == 0,
          "while a single-part feature is part 0 and is not labelled by part")
    check(m.parts[1].label == "feature 0 part 1",
          "so a duplicate can be pointed at by part")
    check(len(m.duplicates) == 1 and m.duplicates[0].index == 1,
          "a LineString duplicating ONE member of a MultiLineString is caught")
    check(m.duplicates[0].relation == IDENTICAL,
          "and named against the member it duplicates")
    check(m.net_mm == m.parts[0].length_mm + m.parts[1].length_mm,
          "the certified total counts each member of the MultiLineString once")
    check(certify([feature([SR40, SR40_BACKWARDS], kind="MultiLineString")]
                  ).collapsed_mm > 0,
          "a MultiLineString that duplicates itself is collapsed too")

    # ---- what gets skipped, and what does not go quietly
    messy = [feature(SR40),
             {"type": "Feature", "properties": {},
              "geometry": {"type": "Polygon", "coordinates": [SR40]}},
             {"type": "Feature", "properties": {}, "geometry": None},
             "not a feature at all",
             feature([[0.0, 0.0], [179.7, 0.5]])]
    s = certify(messy)
    check(len(s.parts) == 1, "only the one real centerline is measured")
    check(len(s.skipped) == 4, "and the other four are skipped, not counted")
    check([i for i, _ in s.skipped] == [1, 2, 3, 4],
          "each skipped feature is named by its position in the file")
    check("Polygon" in dict(s.skipped)[1], "the Polygon says what it was")
    check("no geometry" in dict(s.skipped)[2], "the null geometry says so")
    check("not a GeoJSON object" in dict(s.skipped)[3], "so does the string")
    check("antipodal" in dict(s.skipped)[4],
          "and a segment the geodesic refuses is skipped with its reason, it "
          "does not take the whole certification down  <-- pinned defect")
    check(s.net_mm == s.parts[0].length_mm,
          "a skipped feature contributes no mileage")

    # ---- NaN and infinity, which is what a null coordinate becomes in a
    # converter that writes floats. json.loads accepts both without complaint.
    broken = certify([feature([[float("nan"), 29.2], [-82.0, 29.2]]),
                      feature([[float("inf"), 29.2], [-82.0, 29.2]]),
                      feature(SR40)])
    check(len(broken.parts) == 1,
          "a NaN coordinate and an infinite one are not measured")
    check(len(broken.skipped) == 2
          and all("finite" in reason for _, reason in broken.skipped),
          "they are skipped, each naming the coordinate as not "
          "finite  <-- pinned defect")
    check(broken.net_mm == broken.parts[0].length_mm,
          "and neither contributes a millimetre to the certified total")

    # ---- properties the layer does not have
    bare = certify([{"type": "Feature",
                     "geometry": {"type": "LineString", "coordinates": SR40}}])
    check(len(bare.parts) == 1, "a feature with no properties key still counts")
    check(bare.parts[0].properties == {}, "and reads as having no properties")

    # ---- --by, the total split by maintenance authority
    b = certify(layer, by="MAINT_AUTH")
    check(sorted(b.by_totals) == ["COUNTY", "STATE"],
          "--by splits the total into one bucket per field value")
    check(sum(b.by_totals.values()) == b.net_mm,
          "and the buckets add up to the certified total EXACTLY  <-- pinned defect")
    check(b.by_totals["STATE"] == b.parts[0].length_mm,
          "the duplicated state road is in the STATE bucket once, not twice")
    check(b.conflicts == [], "with no conflict, both copies agree on the field")
    check(certify(layer).by_totals == {},
          "without --by there are no buckets at all")
    check(certify(layer, by="NO_SUCH_FIELD").by_totals.keys() == {MISSING},
          "a field the layer does not carry gives one (no value) bucket")
    nulls = certify([feature(SR40, MAINT_AUTH=None)], by="MAINT_AUTH")
    check(list(nulls.by_totals) == [MISSING],
          "a null field value is the (no value) bucket, not the string None")
    numbers = certify([feature(SR40, FUNC_CLASS=7)], by="FUNC_CLASS")
    check(list(numbers.by_totals) == ["7"],
          "a numeric field value becomes its own bucket, keyed as text")

    # ---- --by when the two copies disagree, which is the defect inside the defect
    split = [feature(SR40, MAINT_AUTH="STATE"),
             feature(SR40_BACKWARDS, MAINT_AUTH="COUNTY")]
    c = certify(split, by="MAINT_AUTH")
    check(len(c.conflicts) == 1,
          "two copies of one road filed under two authorities is reported")
    check(sorted(c.conflicts[0][1]) == ["COUNTY", "STATE"],
          "and both values are named, so somebody can decide which is right")
    check(c.by_totals == {"STATE": c.parts[0].length_mm},
          "the surviving copy's value is the one counted, in input order")
    check(sum(c.by_totals.values()) == c.net_mm,
          "and the buckets still add up to the certified total")

    # ---- certify's own guards
    check(certify([]).net_mm == 0, "an empty layer certifies zero miles")
    check(certify([]).gross_mm == 0, "with nothing collapsed")
    check(describe(certify([]))[-1] == "VERDICT: NO DUPLICATES",
          "and renders as clean rather than as an error")
    raises(lambda: certify("not a list"), "a non-list of features raises")
    raises(lambda: certify([], precision=-1), "a negative precision raises")

    # ---- the miles conversion
    check(miles(1609344) == 1.0, "1609344 mm is one international mile exactly")
    check(miles(0) == 0.0, "zero is zero")
    check(abs(miles(1000000) - 0.6213711922) < 1e-9, "a kilometre is 0.62137 mi")
    check(abs(miles(1609347) - 1.0) > 1e-6,
          "1609347 mm is measurably NOT one mile, so the international mile is "
          "in use here and not the US survey mile")
    check(abs(1609.3472186944 / METRES_PER_MILE - 1) < 3e-6,
          "which differs from the survey mile by 2 parts per million, under "
          "the third decimal this tool prints")

    # ---- rendering
    lines = describe(r)
    check(lines[-1] == "VERDICT: COLLAPSED",
          "a report with duplicates renders its verdict last")
    check(any("CERTIFIED mileage" in x for x in lines),
          "the certified total is labelled")
    check(any("mileage before" in x for x in lines),
          "the total before collapsing is printed as well")
    check("  CERTIFIED mileage                    2.417 mi" in lines,
          "the certified mileage is printed to three decimals on its own line")
    check("  mileage before, every part counted   4.231 mi" in lines,
          "and so is the gross, against its own label")
    check("  mileage removed as duplicate         1.814 mi" in lines,
          "and the labelled correction is the DIFFERENCE between them, not "
          "either total  <-- pinned defect")
    check(any(REVERSED in x for x in lines),
          "the duplicate is described as the reversed one")
    check(any("kept    feature 0" in x for x in lines),
          "the report names the copy it kept")
    check(any("removed feature 2" in x for x in lines),
          "and the copy it removed, so the correction is checkable by hand")
    check(len([x for x in describe(certify(messy)) if "SKIPPED" in x]) == 4,
          "every skipped feature is listed in the report")
    check(any("MAINT_AUTH" in x for x in describe(b)),
          "--by names the field it grouped on")
    check(any("CONFLICT" in x for x in describe(c)),
          "and a conflict between copies is printed")
    check(all("CONFLICT" not in x for x in describe(b)),
          "while a layer with no conflict prints none")
    check("Part(feature 0, " in repr(r.parts[0]),
          "a part reprs as its label, its digest and its length")

    # ---- the duplicate features it would write out
    written = duplicate_features(r)
    check(len(written) == 1, "one duplicate part is written out")
    check(written[0]["geometry"]["type"] == "LineString",
          "as a LineString, whatever it came from")
    check(written[0]["geometry"]["coordinates"] == SR40_BACKWARDS,
          "with the geometry exactly as it was digitized, not canonicalised")
    check(written[0]["properties"]["ROADNAME"] == "SR-40",
          "keeping the attributes it already had")
    check(written[0]["properties"]["roadmiles_duplicate_of"] == "feature 0",
          "and naming the feature it duplicates")
    check(written[0]["properties"]["roadmiles_relation"] == REVERSED,
          "and how")
    check(abs(written[0]["properties"]["roadmiles_miles"] - 1.814) < 0.001,
          "and what it is worth in miles")
    check("roadmiles_group" not in layer[2]["properties"],
          "the input feature is NOT annotated, only the copy  <-- pinned defect")
    check(duplicate_features(certify([feature(SR40)])) == [],
          "a layer with no duplicates writes no features")
    multi_written = duplicate_features(m)
    check(multi_written[0]["geometry"]["coordinates"] == SR40,
          "a duplicate of one MultiLineString member writes that member, not "
          "the whole parent  <-- pinned defect")
    check(json.loads(json.dumps(written)) == written,
          "everything written is JSON-serialisable")

    # ---- argument handling. ROADMILES_PRECISION supplies a default, so a
    # developer who exports it once used to see this self-test fail.
    os.environ["ROADMILES_PRECISION"] = "9"
    saved = os.environ.pop("ROADMILES_PRECISION")
    try:
        a = _parse(["roads.geojson"])
        check(a.apply is False, "--apply defaults to OFF")
        check(a.write_duplicates is None, "--write-duplicates defaults to unset")
        check(a.by is None, "--by defaults to unset")
        check(a.precision == 6,
              "--precision defaults to 6 decimal degrees")
        check(a.path == "roads.geojson", "the file is read as a positional")
        check(_parse(["--self-test"]).self_test, "--self-test parses")
        check(_parse(["r.geojson", "--precision", "4"]).precision == 4,
              "--precision is read")
        check(_parse(["r.geojson", "--by", "MAINT_AUTH"]).by == "MAINT_AUTH",
              "--by is read")
        check(_parse(["r.geojson", "--write-duplicates", "d.geojson"]
                     ).write_duplicates == "d.geojson",
              "--write-duplicates is read")
        check(_parse(["r.geojson", "--apply"]).apply is True, "--apply is read")
        check(_parse(["-"]).path == "-", "a bare - is read as stdin")
        os.environ["ROADMILES_PRECISION"] = "4"
        check(_parse(["r.geojson"]).precision == 4,
              "ROADMILES_PRECISION is read  <-- pinned defect")
        check(_parse(["r.geojson", "--precision", "7"]).precision == 7,
              "and the flag beats the environment variable")
        os.environ["ROADMILES_PRECISION"] = "six"
        quiet, sys.stderr = sys.stderr, io.StringIO()
        code = None
        try:
            try:
                _parse(["r.geojson"])
            except SystemExit as exc:
                code = exc.code
        finally:
            sys.stderr = quiet
        check(code == 2,
              "an unreadable ROADMILES_PRECISION exits 2, not a traceback")
    finally:
        os.environ["ROADMILES_PRECISION"] = saved
    check(os.environ.get("ROADMILES_PRECISION") == "9",
          "a developer's own ROADMILES_PRECISION is put back afterwards")
    os.environ.pop("ROADMILES_PRECISION", None)

    # ---- the harness itself. A check() that cannot record a failure would
    # report every defect above as a pass, which is the one failure no other
    # assertion here could see. Three deliberate failures are recorded against a
    # scratch mark and then taken back off the tally.
    quiet = sys.stdout
    sys.stdout = io.StringIO()
    mark = len(failed)
    try:
        check(False, "probe: a false condition must be recorded as a failure")
        raises(lambda: None, "probe: a call that raises nothing must fail")
        raises(lambda: [][0], "probe: the wrong exception must fail")
    finally:
        sys.stdout = quiet
    probe = failed[mark:]
    del failed[mark:]
    check(len(probe) == 3,
          "check() and raises() really do record a failure  <-- pinned defect")
    check("no error raised" in probe[1] and "wrong exception" in probe[2],
          "and say which way the call under test went wrong")

    # ---- the io layer, against real files in a temporary directory.
    #
    # Everything above this line is pure. Everything below writes GeoJSON to
    # disk and reads it back, because a report that has never been produced from
    # a file, and a duplicates layer that has never been written, have not been
    # tested.
    def capture(fn):
        """(what fn returned, everything it printed to either stream)."""
        buf = io.StringIO()
        out, err = sys.stdout, sys.stderr
        sys.stdout = sys.stderr = buf
        try:
            result = fn()
        finally:
            sys.stdout, sys.stderr = out, err
        return result, buf.getvalue()

    tmp = tempfile.mkdtemp(prefix="roadmiles-selftest-")

    def put(name, text):
        path = os.path.join(tmp, name)
        with open(path, "wb") as handle:
            handle.write(text.encode("utf-8"))
        return path

    try:
        collection = json.dumps({"type": "FeatureCollection", "features": layer})
        roads = put("roads.geojson", collection)
        check(len(read_features(roads)) == 3,
              "a FeatureCollection on disk reads back as three features")
        check(read_features(roads)[0]["geometry"]["coordinates"] == SR40,
              "with its coordinates intact")
        check(len(read_features(put("bare.geojson", json.dumps(layer)))) == 3,
              "a bare JSON array of features reads too")
        one = put("one.geojson", json.dumps(layer[0]))
        check(len(read_features(one)) == 1, "a single Feature reads as a layer of one")
        geom = put("geom.geojson", json.dumps({"type": "LineString",
                                               "coordinates": SR40}))
        check(read_features(geom)[0]["geometry"]["coordinates"] == SR40,
              "a bare geometry is wrapped into a feature")
        bom = put("bom.geojson", "\ufeff" + collection)
        check(len(read_features(bom)) == 3,
              "a file ArcGIS wrote with a UTF-8 byte order mark reads  <-- pinned defect")
        raises(lambda: read_features(put("num.geojson", "42")),
               "a JSON number is not GeoJSON and raises")
        raises(lambda: read_features(put("pt.geojson",
                                         json.dumps({"type": "Point"}))),
               "a bare Point is not a road layer and raises")
        raises(lambda: read_features(put("nofc.geojson",
                                         json.dumps({"type": "FeatureCollection"}))),
               "a FeatureCollection with no features array raises")

        # ---- main, end to end
        code, out = capture(lambda: main([roads]))
        check(code == 1, "a layer carrying a duplicate exits 1")
        check("CERTIFIED mileage" in out and "2.417 mi" in out,
              "and prints the certified total, read from a real file")
        check("4.231 mi" in out and "1.814 mi" in out,
              "with the gross and the correction beside it")
        check("VERDICT: COLLAPSED" in out, "and says a duplicate was collapsed")
        clean = put("clean.geojson",
                    json.dumps({"type": "FeatureCollection",
                                "features": [layer[0], layer[1]]}))
        code, out = capture(lambda: main([clean]))
        check(code == 0, "a layer with no duplicates exits 0")
        check("VERDICT: NO DUPLICATES" in out, "and says so")
        check("2.417 mi" in out,
              "and certifies the same total the corrected layer did  <-- pinned defect")

        code, out = capture(lambda: main([roads, "--by", "MAINT_AUTH"]))
        check(code == 1 and "certified mileage by MAINT_AUTH" in out,
              "--by prints the split end to end")
        check("STATE" in out and "COUNTY" in out, "with a line for each value")
        near = put("near.geojson", json.dumps(
            {"type": "FeatureCollection",
             "features": [feature([[-82.1, 29.2], [-82.0, 29.2]]),
                          feature([[-82.1004, 29.2], [-82.0, 29.2]])]}))
        code, out = capture(lambda: main([near]))
        check(code == 0 and "VERDICT: NO DUPLICATES" in out,
              "two roads four ten-thousandths of a degree apart are two roads "
              "at the default precision")
        code, out = capture(lambda: main([near, "--precision", "3"]))
        check(code == 1 and "VERDICT: COLLAPSED" in out,
              "and one road at --precision 3, so the flag reaches the hash and "
              "changes the answer end to end  <-- pinned defect")
        code, out = capture(lambda: main([roads, "--precision", "12"]))
        check(code == 1,
              "a precision finer than the data still finds an exactly reversed copy")

        code, out = capture(lambda: main(["no-such-file-4f2a1c.geojson"]))
        check(code == 2, "a file that does not exist exits 2")
        check("could not read" in out, "and says which file")
        code, out = capture(lambda: main([put("bad.geojson", "{not json")]))
        check(code == 2, "a file that is not JSON exits 2, not a traceback")
        code, out = capture(lambda: main([]))
        check(code == 64, "no file at all is a usage error")
        code, out = capture(lambda: main([roads, "--precision", "-1"]))
        check(code == 64, "a negative --precision is a usage error")

        # ---- reading the layer from stdin, which is how a curl gets in here
        real_stdin = sys.stdin
        try:
            sys.stdin = io.StringIO(collection)
            code, out = capture(lambda: main(["-"]))
        finally:
            sys.stdin = real_stdin
        check(code == 1 and "2.417 mi" in out,
              "the same layer piped in on stdin certifies the same total")

        # ---- --write-duplicates, which writes nothing without --apply
        target = os.path.join(tmp, "dupes.geojson")
        code, out = capture(lambda: main([roads, "--write-duplicates", target]))
        check(code == 1 and not os.path.exists(target),
              "--write-duplicates without --apply writes NOTHING  <-- pinned defect")
        check("Re-run with --apply" in out, "and says how to write it")
        code, out = capture(
            lambda: main([roads, "--write-duplicates", target, "--apply"]))
        check(code == 1 and os.path.isfile(target),
              "--write-duplicates --apply writes the file")
        with open(target, "rb") as handle:
            raw = handle.read()
        written_back = json.loads(raw.decode("utf-8"))
        check(written_back["type"] == "FeatureCollection",
              "what it wrote is a FeatureCollection")
        check(len(written_back["features"]) == 1, "holding the one duplicate part")
        check(written_back["features"][0]["properties"]["roadmiles_relation"]
              == REVERSED, "annotated with how it duplicates the survivor")
        check(written_back["features"][0]["geometry"]["coordinates"]
              == SR40_BACKWARDS,
              "and carrying the geometry as digitized, so it can be selected "
              "in Pro and deleted")
        check(b"\r\n" not in raw,
              "written with LF on every platform, so two desks produce the same "
              "file  <-- pinned defect")
        check(certify(read_features(target)).net_mm == r.collapsed_mm,
              "and certifying the duplicates file gives back exactly the mileage "
              "that was removed  <-- pinned defect")
        code, out = capture(
            lambda: main([roads, "--write-duplicates", target, "--apply"]))
        check(code == 2 and "already exists" in out,
              "a second --apply refuses to overwrite the file")
        check(json.loads(open(target, "rb").read().decode("utf-8")) == written_back,
              "and leaves what was there alone")
        nowhere = os.path.join(tmp, "no-such-dir", "d.geojson")
        code, out = capture(lambda: write_duplicates(nowhere, written, True))
        check(code == 2 and "could not write" in out,
              "a path that cannot be written exits 2, not a traceback")
        empty = os.path.join(tmp, "none.geojson")
        code, out = capture(
            lambda: main([clean, "--write-duplicates", empty, "--apply"]))
        check(code == 0 and os.path.isfile(empty),
              "a clean layer still writes an empty duplicates file, so a "
              "scheduled run has something to point at")
        check(json.loads(open(empty, "rb").read().decode("utf-8"))["features"]
              == [], "which holds no features")

        # ---- nothing secret is in the output, checked against the artefact
        check(b"password" not in raw.lower() and b"token" not in raw.lower(),
              "no credential reaches the file it writes, there is none to reach it")
    finally:
        shutil.rmtree(tmp, ignore_errors=True)
        check(not os.path.isdir(tmp),
              "the self-test leaves no temporary directory behind")

    print("-" * 68)
    total = passed[0] + len(failed)
    if failed:
        print("%d assertions, %d failed" % (total, len(failed)))
        for item in failed:
            print("  FAILED: %s" % item)
        return 1
    print("%d assertions, 0 failed" % total)
    return 0


# ----------------------------------------------------------------------- cli

def _parse(argv):
    ap = argparse.ArgumentParser(
        prog="roadmiles.py",
        description="Certify centerline mileage and refuse to count the same "
                    "road twice, including the copy someone digitized backwards.",
        epilog="Config precedence: flag > environment > the CONFIGURATION "
               "block. Nothing is written without --apply.",
    )
    ap.add_argument("path", nargs="?", metavar="FILE",
                    help="GeoJSON road layer, or - to read it from stdin")
    ap.add_argument("--precision", type=int,
                    default=os.environ.get("ROADMILES_PRECISION",
                                           DEFAULT_PRECISION),
                    help="decimal places coordinates are rounded to before "
                         "hashing (default %d). Env: ROADMILES_PRECISION"
                         % DEFAULT_PRECISION)
    ap.add_argument("--by", metavar="FIELD",
                    help="also break the certified total down by this property, "
                         "for example MAINT_AUTH")
    ap.add_argument("--write-duplicates", dest="write_duplicates",
                    metavar="PATH",
                    help="write the parts that were removed as a GeoJSON "
                         "FeatureCollection. Needs --apply to write.")
    ap.add_argument("--apply", action="store_true",
                    help="write the duplicates file. Without this nothing is "
                         "written.")
    ap.add_argument("--self-test", dest="self_test", action="store_true",
                    help="run the offline assertions and exit")
    # argparse runs type= over the default as well, so an unreadable
    # ROADMILES_PRECISION is argparse's own usage error rather than a traceback.
    return ap.parse_args(argv)


def main(argv=None):
    args = _parse(sys.argv[1:] if argv is None else argv)

    if args.self_test:
        return self_test()

    if not args.path:
        print("error: pass a GeoJSON road layer, or - for stdin. Use "
              "--self-test to verify the tool without one.", file=sys.stderr)
        return 64
    if args.precision < 0:
        print("error: --precision cannot be negative.", file=sys.stderr)
        return 64

    try:
        features = read_features(args.path)
    except (OSError, IOError, ValueError) as exc:
        print("error: could not read %s: %s" % (args.path, exc), file=sys.stderr)
        return 2

    report = certify(features, args.precision, args.by)
    print("roadmiles: %s" % args.path)
    for line in describe(report):
        print(line)

    if args.write_duplicates:
        print("")
        code = write_duplicates(args.write_duplicates,
                                duplicate_features(report), args.apply)
        if code:
            return code

    return 1 if report.duplicates else 0


if __name__ == "__main__":
    sys.exit(main())
