# roadmiles

Certify centerline mileage and refuse to count the same road twice, including the copy someone digitized backwards.

The county certifies total road centerline mileage to the state once a year, and the money that
comes back is proportional to the number. Sixty-one segments in the roads layer exist twice.
Usually two editors digitized the same road several years apart, and the second one drew it the
other way round: east to west instead of west to east.

Nothing about the map looks wrong. The two copies sit exactly on top of each other, one line
hidden under the other, and no rendering, no label and no selection colour gives it away. The
total is overstated, the certification is wrong, and the correction is not a number anybody in
the office can produce.

A spatial join on geometry equality does not find them. To ArcGIS, to Shapely and to a
geodatabase, a LineString is an **ordered** list of positions, so the same road recorded west to
east and east to west is two different geometries. Their coordinate lists are reverses of each
other and compare unequal, one field at a time, down the whole array. "Find Identical" on `SHAPE`
compares the geometry as stored and reports nothing.

This computes one direction-independent hash per line, so both copies land in the same bucket in
a single pass, and it measures the survivors geodesically on the WGS84 ellipsoid. It reports
mileage before, mileage after, and the mileage it removed, because a certification correction has
to be a number somebody can defend line by line.

```
$ python roadmiles.py --self-test
roadmiles self-test: no network, no database, no arcpy
--------------------------------------------------------------------
PASS  one degree of latitude at the equator is 110574.4 m, within a metre
PASS  the 44 to 45 degree meridian arc is 111122.0 m, within a metre
PASS  the 89 to 90 degree meridian arc is 111693.9 m, within a metre
PASS  planar degrees overstate that same east-west mile by over 14 km, 14 percent, with no error message
...
PASS  a segment crossing the antimeridian is 222 m, not half the planet  <-- pinned defect
PASS  a nearly antipodal pair raises instead of returning a wrong number
PASS  and a NaN is named as not finite, not accused of being nearly antipodal  <-- pinned defect
PASS  and formats without a minus sign, so it hashes with its neighbour  <-- pinned defect
PASS  a coordinate list and its reverse hash EQUAL  <-- pinned defect
PASS  a vertex nudged PAST the rounding precision hashes differently
PASS  a vertex nudged UNDER the rounding precision hashes the same
PASS  a closed loop starting at a different vertex is the SAME road  <-- pinned defect
PASS  rotated AND reversed is still the same road  <-- pinned defect
PASS  and last year's report can still be matched against this year's loops  <-- pinned defect
...
PASS  so the before/after difference equals the collapsed duplicate mileage EXACTLY  <-- pinned defect
PASS  every total is an integer count of millimetres, so no summation order can move it  <-- pinned defect
PASS  a degree of latitude certifies as 110574389 mm, the meridian arc ROUNDED to the millimetre and not truncated  <-- pinned defect
PASS  and TWO of them are removed, not one  <-- pinned defect
PASS  the certified total is the FIRST copy's length, not the last one's  <-- pinned defect
PASS  and a segment the geodesic refuses is skipped with its reason, it does not take the whole certification down  <-- pinned defect
PASS  and the buckets add up to the certified total EXACTLY  <-- pinned defect
PASS  and the labelled correction is the DIFFERENCE between them, not either total  <-- pinned defect
PASS  the input feature is NOT annotated, only the copy  <-- pinned defect
PASS  a duplicate of one MultiLineString member writes that member, not the whole parent  <-- pinned defect
PASS  check() and raises() really do record a failure  <-- pinned defect
...
PASS  a file ArcGIS wrote with a UTF-8 byte order mark reads  <-- pinned defect
PASS  and one road at --precision 3, so the flag reaches the hash and changes the answer end to end  <-- pinned defect
PASS  --write-duplicates without --apply writes NOTHING  <-- pinned defect
PASS  written with LF on every platform, so two desks produce the same file  <-- pinned defect
PASS  and certifying the duplicates file gives back exactly the mileage that was removed  <-- pinned defect
PASS  the self-test leaves no temporary directory behind
--------------------------------------------------------------------
244 assertions, 0 failed
```

## Requirements

Python 3.9 or newer. Nothing to install: no `arcpy`, no `shapely`, no `pyproj`, no third-party
package at all. It uses `json`, `math`, `hashlib` and `argparse`, and it runs on ArcGIS Pro's
Python and on a plain `python3` equally. The self-test needs no network, no database and no
geodatabase, and produces the same 244 assertions on Windows and on Linux.

```
git clone https://github.com/uhsear/roadmiles.git
```

## Quick start

```
python roadmiles.py --self-test
python roadmiles.py roads.geojson
```

Export the road layer from Pro as GeoJSON, or pipe one in:

```
curl -s 'https://services.arcgis.com/.../query?where=1=1&outFields=*&f=geojson' | python roadmiles.py -
```

## Usage

```
python roadmiles.py roads.geojson
python roadmiles.py roads.geojson --by MAINT_AUTH
python roadmiles.py roads.geojson --write-duplicates dupes.geojson --apply
```

| Flag | Default | What it does |
|---|---|---|
| `FILE` | none | GeoJSON road layer, or `-` to read it from stdin. Required. |
| `--precision` | `6` | Decimal places coordinates are rounded to before hashing. Env: `ROADMILES_PRECISION` |
| `--by` | none | Also break the certified total down by this property, for example `MAINT_AUTH`. |
| `--write-duplicates` | none | Write the parts that were removed as a GeoJSON FeatureCollection. |
| `--apply` | off | Write that file. Without it nothing is written. |
| `--self-test` | off | Run the offline assertions and exit. |

Exit codes: 0 no duplicates, 1 duplicates were collapsed, 2 a file could not be read or written,
64 usage error. The certified total is printed either way. **Exit 1 means the number changed, not
that the run failed**, so a scheduled check can treat it as "somebody has editing to do".

A run over a nine feature test layer, three of whose roads are entered twice:

```
$ python roadmiles.py riverton_roads.geojson --by MAINT_AUTH
roadmiles: riverton_roads.geojson
9 part(s) from 9 feature(s), 1 skipped
rounded to 6 decimal degrees, measured geodesically on WGS84

  SKIPPED feature 8: geometry type Polygon is not a centerline

  mileage before, every part counted   9.145 mi
  duplicate groups collapsed           3 group(s), 3 part(s)
  mileage removed as duplicate         3.280 mi
  ----------------------------------------------------
  CERTIFIED mileage                    5.865 mi

duplicates, in input order:
  09f5e71b0637  2 copies, 1.814 mi
      kept    feature 0
      removed feature 3, reversed
  fc98584dcbe1  2 copies, 1.208 mi
      kept    feature 1
      removed feature 6 part 1, reversed
  673a6249bb95  2 copies, 0.259 mi
      kept    feature 4
      removed feature 5, rotated loop

certified mileage by MAINT_AUTH:
  CITY                         0.605 mi
  COUNTY                       2.844 mi
  PRIVATE                      0.603 mi
  STATE                        1.814 mi
  CONFLICT fc98584dcbe1: copies of one road disagree about MAINT_AUTH: COUNTY, CITY
  CONFLICT 673a6249bb95: copies of one road disagree about MAINT_AUTH: COUNTY, CITY
VERDICT: COLLAPSED
```

Every duplicate is named by the feature it duplicates and by **how**. "Reversed" says these are
two digitizing sessions and not one feature copied and pasted, which is a different conversation
about how it got there.

## The direction-independent geometry hash

This is the reusable idea and it is small on purpose.

1. **Round** every position to `--precision` decimal places. Six is about 11 cm at the equator,
   finer than any centerline was digitized and coarse enough to absorb the float noise a
   reprojection leaves behind.
2. **Squeeze out** a vertex that is now equal to the one before it, so that a double-clicked
   vertex does not make two copies of a road different.
3. **Canonicalise the direction.** An open line has exactly two readings, forwards and backwards,
   so the smaller of the two is a canonical name for the pair. A closed loop has two readings for
   every vertex it could have started at, because the start vertex of a ring is an arbitrary
   choice the editor made, so the candidates are every rotation of the ring in both directions.
4. **Hash** the formatted decimal text of that one canonical vertex list with SHA-256.

Two copies of one road hash equal whichever way each was drawn. That is the whole trick, and it
runs in one pass over the layer rather than in a pairwise comparison.

The closed loop is the case that gets missed, and the self-test pins it. A cul-de-sac or a
roundabout digitized by two editors starts at whichever vertex each of them clicked first. Only
reversing the list, without rotating the ring, leaves every duplicated loop in the layer
uncollapsed and the mileage still overstated. Rotating without reversing misses half of them. The
tool does both, and there are assertions for rotated, for reversed and for rotated **and**
reversed, because a fixture that only ever tests one of the three leaves the other two free to
break with nothing going red.

## Why the lengths are geodesic

Planar distance over degrees is the other way this number silently comes out wrong.

A degree of longitude is 111.3 km at the equator and 97.3 km at a central Florida latitude, because
the meridians converge. Treat degrees as a flat grid, with whatever metres-per-degree constant is
lying around, and every east-west road here is overstated by about 14 percent. Every north-south
road is roughly right. There is no error message and no failed run: the certification is simply
larger than the roads are, by an amount that depends on which way each road happens to run, which
is why you cannot divide it back out afterwards.

Reprojecting to State Plane first is a correct answer to this and a different tool's job. This one
measures on the ellipsoid where the coordinates already are, with Vincenty's inverse formula
written out in about forty lines. The self-test checks it against the published WGS84 meridian
arc at three latitudes and holds it to within a metre in 110 km.

Two more cases the formula has to survive, both asserted:

- **A segment crossing the antimeridian.** The longitude difference arrives as -359.998 degrees
  rather than 0.002. Everything that uses it goes through `sin` or `cos`, which do not care, so a
  222 m road stays a 222 m road. The naive subtraction that returns half the planet is not in
  there.
- **A coordinate that is not a finite number.** `NaN` and `Infinity` are floats,
  `json.loads` accepts both without complaint, and both reach here from a converter
  that turned a null into a float. Left alone a `NaN` runs every iteration and comes
  back accused of being nearly antipodal, and an infinity comes back as `math.sin`'s
  "math domain error". Both are refused by name instead, and the feature is skipped
  with that reason.
- **A nearly antipodal pair**, where Vincenty does not converge. It raises rather than return a
  number, and the feature is reported as skipped with its reason instead of taking the whole
  certification down. No road centerline is 20,000 km long, so this means a coordinate is wrong
  or the axes are swapped.

## Why the totals are integers

Lengths are carried as integer millimetres, and miles are computed from them only for display.

The headline of the report is a subtraction: mileage before, minus mileage after, equals the
mileage removed. In floating point that identity depends on the order the parts were summed in,
so the three numbers can fail to reconcile in the last decimal for no reason anyone reading the
report could explain. Over integers it is an identity. `--by` buckets add up to the certified
total exactly, for the same reason.

The mile is the international mile, 1609.344 m exactly, set in the `CONFIGURATION` block so the
choice is visible. The US survey mile is 1609.3472186944 m, and Florida's state plane feet are
survey feet, so a pipeline that converts through them is working in the other one. The gap is
3.2 mm per mile, under the third decimal this tool prints.

## What it will not do

- **It does not decide which copy to delete.** It keeps the first in input order and writes the
  rest out with `--write-duplicates`, annotated with the feature they duplicate and how, so they
  can be selected in Pro and reviewed. Deleting is an editing session with a person in it.
- **It never edits the input.** `--apply` writes one new file, the duplicates layer, and refuses
  to overwrite a file that already exists. The features it writes are copies, so the layer the
  caller still holds is not quietly annotated.
- **It does not find near-duplicates.** Two independent digitizations of one road, a metre apart
  the whole way, are two different roads to this tool at any sensible `--precision`. That is a
  geometric similarity problem and it needs a tolerance, a Hausdorff distance and somebody's
  judgement. This finds the copy, which is the common case and the tractable one.
- **It does not know a road from a driveway.** It totals what is in the file. Filtering to the
  features that are certifiable is a `--by` report at most, and a definition query at best.
- **It does not split a road at an intersection or dissolve one across segments.** Segmentation is
  whatever the layer already has, and the total is the same either way.
- **No network and no credentials.** It reads a file or stdin. Pulling the layer out of ArcGIS
  Online is `curl` or `arcpy`, both of which already do it, and neither of which should have to be
  installed to check a number.

## Limits

- `--precision` is a real knob with two ways to be wrong. Too fine and two copies that went
  through different reprojections stop matching. Too coarse and separate roads start matching:
  at `--precision 0` every road in the county collapses into a handful of groups. Six is the
  default because it is finer than any centerline was drawn.
- An extra vertex along a straight run makes a copy hash differently, even though it is the same
  road on the ground. Only a vertex that lands on its neighbour after rounding is squeezed out.
  This is the boundary between "a copy" and "a redigitization", and this tool is on the copy side
  of it.
- A `MultiLineString` counts as one part per member, and a duplicate of one member is one
  duplicate. `--write-duplicates` writes that member out as a `LineString`, not the parent
  feature, because the parent also holds geometry that is not duplicated and is not to be deleted.
- Coordinates are read as longitude first, which is what GeoJSON means (RFC 7946). A layer written
  latitude first is caught only when a value goes outside -90 to 90, and a layer confined to
  Florida never does. The mileage comes out wrong and nothing says so.
- A z value is carried in the data, ignored in the length and ignored in the hash. Slope length is
  not what a centerline certification measures, and the same road surveyed twice carries two
  elevations.
- The whole layer is held in memory. A county road network is tens of thousands of features and
  fits easily; a statewide one might not.
- The loop canonicalisation compares every rotation of a ring in both directions, which is
  quadratic in the vertex count of that one loop. A roundabout has tens of vertices. A ring with
  ten thousand of them is not a road, and would be slow.

## Contributing

Open an issue or pull request on GitHub.

## Author

Built by [Asir Khan](https://www.linkedin.com/in/asir-khan-310317264/).

## License

MIT.
