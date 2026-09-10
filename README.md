# TriDecomp

Triangular decomposition of constraint segments inside an enclosing triangle. The **C++ engine does one thing**: `tessellate(segments, t)`. Walk, marks, I/O, and enclosing-triangle helpers are Python.

The same source still builds a **C++ CLI** that reads TIGER/ESRI shapefiles, optionally simplifies, decomposes, classifies leaves, and writes AGS transaction logs (see [C++ CLI](#c-cli-shapefile-pipeline)).

---

## Table of Contents

1. [Architecture](#architecture)
2. [Python library](#python-library)
3. [Prerequisites](#prerequisites)
4. [Building](#building)
5. [C++ CLI (shapefile pipeline)](#c-cli-shapefile-pipeline)
6. [Split Selection: Angle Method vs. Legacy Lookahead](#split-selection-angle-method-vs-legacy-lookahead)
7. [Pipeline Stages](#pipeline-stages)
8. [Output Files](#output-files)
9. [Reconstructing Transactions](#reconstructing-transactions)
10. [State FIPS Reference](#state-fips-reference)
11. [County FIPS Reference](#county-fips-reference)
12. [Examples](#examples)
13. [Troubleshooting](#troubleshooting)

---

## Architecture

| Layer | Language | Role |
|---|---|---|
| **`tessellate(segments, t, start="")`** | C++ `tessellate.cpp` (CGAL **EPEC**) | Split `t` against untyped constraint segments. Returns a list of `SplitEvent`s with tree addresses and nested line splits. |
| **`tessellate(..., engine="python")`** | Python (CGAL **EPIC**) | Same contract; A/B stand-in, not the chain-fidelity path. |
| **`edges_of`, `tridecomp`, `walk`, `tiles`, `gettri`, `read_tiger_polygon`** | Python | Polygon sugar, PerWalk marks, enclosing triangle, shapefile I/O. |

A polygon ring is just one source of segments. Delegation / refinement is the same primitive with extra segments and a non-empty `start` code.

**Invariants of a split**

- The two children partition the parent (disjoint interiors, union = parent).
- From the directed cut, **`0` = left child**, **`1` = right child**.
- A cut that **crosses** a constraint interior splits it; one stub goes to each child (a `line_splits` record).
- A cut that **contains both endpoints** of a constraint (the segment lies on the new shared side, possibly as a proper subsegment / T-junction) copies that segment into **both** children. No new vertex, no line-split record.

---

## Python library

### Install

From the repo root, after [building](#building) the extension:

```bash
make python
pip install -e . --no-build-isolation --no-deps
```

Import is then `from tridecomp import ...` (editable package under `python/tridecomp/`).

### Segments: shape `(M, 2, 2)`

Each row is one constraint edge — **two endpoints**, each `(x, y)`:

```
segments[i] = [[x0, y0],     # endpoint A
               [x1, y1]]     # endpoint B
```

That is the same four numbers, row-major, as shape `(M, 4)`: `x0, y0, x1, y1`. Tessellation treats the edge as undirected. Cycle **direction** matters later in `walk`, not here.

`edges_of(p)` builds this from a closed ring `p` of shape `(N, 2)` (CCW; a repeated closing vertex is dropped): `segments[i]` runs from `p[i]` to `p[i+1]`.

### `tessellate(segments, t, start="")`

```python
from tridecomp import tessellate, tridecomp, edges_of, gettri, walk, tiles
from tridecomp import read_tiger_polygon

p = read_tiger_polygon("tl_2025_us_county.zip", state="31", county="039")
t = gettri(p, rat=1.5)
events = tessellate(edges_of(p), t)          # or tridecomp(p, t)
# events = tessellate(edges_of(p), t, engine="python")  # EPIC stand-in
```

- **`t`**: shape `(3, 2)`. Clockwise input is reversed to CCW.
- **`start`**: tree address of `t`. Empty (`""`) for a fresh enclosing triangle. Pass a leaf code (e.g. `"01011"`) to split that tile further after adding segments inside it; child ids are `start+"0"` / `start+"1"`.
- **`engine`**: `"cpp"` (default, EPEC) or `"python"` (EPIC).

Each **SplitEvent** is a dict:

| Key | Meaning |
|---|---|
| `id` | Bitstring of the parent (`""` is the original enclosing triangle, or `start` on a refinement). |
| `child0` / `child1` | `id+"0"` (left of the cut) and `id+"1"` (right). |
| `parent` | `(3, 2)` triangle being split. |
| `cut` | `(2, 2)` split segment on the parent. |
| `child` | `(left, right)` child triangles, each `(3, 2)`. |
| `line_splits` | Nested list: constraint edges this cut **crossed** (new vertex). Each has `edge` `(2, 2)`, `point` `(2,)`, `seg` (two stubs). |

Leaves are implied: a triangle that never appears later as a `parent`. Replay with `tiles(events, t, start="")` → `{id: triangle}`.

### Refinement

```python
leaf_id = "01011"
leaf_tri = tiles(events, t)[leaf_id]
more = tessellate(extra_segments, leaf_tri, start=leaf_id)
all_events = events + more
leaves = tiles(all_events, t)          # original tree; that leaf is replaced
```

### `walk(events, t, cycle)`

Walks a directed cycle on the constraint complex (typically the same ring you tessellated).

- Tile to the **left** of each cycle edge → `INSIDE` (only if unmarked).
- Tile to the **right** → `OUTSIDE` (only if unmarked).
- Marks are **monotonic**: `UNMARKED → INSIDE/OUTSIDE` or unchanged. `INSIDE ↔ OUTSIDE` is a crossing; `simple` is then `False`.
- Remaining unmarked tiles flood from those seeds without toggling.

```python
from tridecomp import UNMARKED, INSIDE, OUTSIDE

w = walk(events, t, p)
w["simple"]       # True if the cycle did not cross itself relative to the tessellation
w["marks"]        # id -> UNMARKED / INSIDE / OUTSIDE
w["n_inside"], w["n_outside"]
```

The cycle edges must lie on the complex after line splits, so `tessellate` must have been given those segments if you want to certify the cycle.

### Helpers

| Function | Purpose |
|---|---|
| `tridecomp(p, t, start="")` | `tessellate(edges_of(p), t, start=start)` |
| `gettri(p, rat=1.5)` | Enclosing triangle `(3, 2)` for vertices `p` |
| `read_tiger_polygon(path, state=, county=)` | Exterior ring from a TIGER shapefile or zip (resolves relative names against cwd and the repo root) |

---

## Prerequisites

| Dependency | Purpose | Windows (vcpkg) | Linux (apt) |
|---|---|---|---|
| C++17 compiler | Build | Visual Studio 2022 | `g++` (13+) |
| CGAL | Exact geometry kernel | `vcpkg install cgal:x64-windows` | `libcgal-dev` (5.6) |
| shapelib | Shapefile reader (CLI) | `vcpkg install shapelib:x64-windows` | `libshp-dev` |
| GMP / MPFR | CGAL exact numbers | pulled in by vcpkg | `libgmp-dev` `libmpfr-dev` |
| Python ≥ 3.10 + numpy | Library | | conda/pip |
| CGAL Python (optional) | `engine="python"` EPIC | | conda-forge `cgal-cpp` / `import CGAL` |
| pyogrio (optional) | `read_tiger_polygon` | | conda/pip |

If you haven't set up vcpkg yet:

```powershell
git clone https://github.com/microsoft/vcpkg.git C:\vcpkg
cd C:\vcpkg
.\bootstrap-vcpkg.bat
.\vcpkg integrate install
```

Then install the dependencies:

```powershell
vcpkg install cgal:x64-windows shapelib:x64-windows
```

### Getting Shapefiles

Download TIGER/Line Shapefiles from the U.S. Census Bureau:

1. Go to https://www.census.gov/cgi-bin/geo/shapefiles/index.php
2. Select a year (e.g. 2024)
3. Select a layer type (e.g. "States", "Counties")
4. Select your state or download the national file
5. Unzip the download — you'll get `.shp`, `.dbf`, `.shx`, and other files

Keep all companion files (`.shp`, `.dbf`, `.shx`) in the same folder. The tool reads the `.shp` file and automatically finds the `.dbf` alongside it.

---

## Building

### Python extension (`make python`)

Builds `python/tridecomp/_cpp*.so` from `tessellate.cpp` (CGAL EPEC + GMP/MPFR; no shapelib). Default interpreter is `$HOME/miniconda3/envs/tridecomp-py/bin/python` (`PY=...` to override). The shapefile CLI is a separate binary (`make` / `TriDecomp.cpp`).

```bash
make python
pip install -e . --no-build-isolation --no-deps
```

### Visual Studio 2022 (CLI)

1. Double-click `TriDecomp.sln`
2. Set the configuration to **Release | x64** (toolbar dropdown)
3. Build → Build Solution (Ctrl+Shift+B)
4. The executable appears in `bin\Release\TriDecomp.exe`

### Command-Line (Developer Prompt)

```
cl /std:c++17 /EHsc /O2 /DNOMINMAX /D_CRT_SECURE_NO_WARNINGS ^
   TriDecomp.cpp /Fe:TriDecomp.exe ^
   /I"C:\vcpkg\installed\x64-windows\include" ^
   /link /LIBPATH:"C:\vcpkg\installed\x64-windows\lib" shp.lib
```

Adjust the vcpkg paths to match your installation.

> **Note:** This source compiles and links under g++ 13 with CGAL 5.6 + libshp on Linux. The intersection code uses the same `CGAL::Object` / `CGAL::object_cast` pattern throughout, so it builds consistently across both toolchains. Your authoritative Windows build remains VS2022 + vcpkg (x64).

### Linux (g++ / Make)

**Option A — system packages** (needs sudo):

```bash
sudo apt-get install -y g++ make libcgal-dev libshp-dev libgmp-dev libmpfr-dev
make
```

**Option B — no root** (extracts Ubuntu `.deb` files into `~/.local/tridecomp-deps`):

```bash
make deps
make
```

The Makefile prefers `/usr/include/CGAL` when present; otherwise it uses the local prefix and statically links shapelib, GMP, and MPFR so you do not need `LD_LIBRARY_PATH`.

Smoke-test the CLI with a synthetic 4-vertex square:

```bash
make test
```

Then run against a real shapefile the same way as on Windows, substituting `./TriDecomp` for `TriDecomp.exe`:

```bash
./TriDecomp tl_2024_us_state.shp 50 --state 53
```

---

## C++ CLI (shapefile pipeline)

**Shapefile → (optional) Polygon Reduction → CGAL Header Generation → Triangular Decomposition → PerWalk + MarkTriangles Classification → Fixed-Point Quantization → Walkable Boundary + Transaction Logging**

The CLI still runs the full AGS logging path in-process (marks, quantization, CSVs). The **library** path does not: `tessellate` returns only split events; `walk` is Python.

### Usage

```
TriDecomp.exe <shapefile.shp> [reduction_%] [options]
```

### Positional Arguments

| Argument | Description |
|---|---|
| `<shapefile.shp>` | Path to the ESRI/TIGER shapefile |
| `[reduction_%]` | *Optional.* Percentage of vertices to remove (0–100). Omit for full resolution. |

### Options

| Option | Default | Description |
|---|---|---|
| `--state <FIPS>` | *(all records)* | Filter by 2-digit state FIPS code (e.g. `53` = WA) |
| `--county <FIPS>` | *(all in state)* | Filter by 3-digit county FIPS code (e.g. `033` = King). Use with `--state`. |
| `--num-states <N>` | *(off)* | **Batch mode**: process the first `N` states found in the shapefile, each into its own `state_<FIPS>/` folder, then write a `summary.csv` and print a summary table |
| `--lookahead-range <A-B>` | *(off)* | **Depth sweep**: run each state at every look-ahead depth `A..B`, each into a `la_<D>/` folder, then print a per-state trend of how leaf count falls with depth. Combines with `--num-states`. |
| `--angle-method` | **on** | Use Dr. Ramkumar's angle method for split selection (default) |
| `--legacy-lookahead` | | Use the legacy recursive lookahead engine instead (for A/B comparison) |
| `--lookahead <N>` | `0` | Lookahead depth for **both** the angle method and the legacy engine (`0` = greedy). Higher `N` yields lower-leaf decompositions at increasing runtime cost. |
| `--no-quantize` | | Skip the fixed-point quantization / walkable-boundary publishing step |
| `--quiet` | | Progress + AGS records + CSVs only; skip per-iteration debug dumps. Logs are buffered (flushed at the end of the run). |
| `--verbose` | | Dump every cell's vertices and edges to `debug_out_lookahead.txt` |
| `--reduce <pct>` | *(positional)* | Alternative way to specify reduction % |
| `--header <path>` | `polygon_output.h` | Output CGAL header file path |
| `--points <path>` | `polygon_output.txt` | Output points file path |
| `--prefix <name>` | `simple_polygon` | C++ variable name in the header file |
| `--no-header` | | Skip CGAL header file generation |
| `--no-points` | | Skip points file generation |
| `--list-states` | | Show all state FIPS codes in the shapefile, then exit |
| `--list-counties` | | Show county FIPS codes and vertex counts (use with `--state`), then exit |
| `--help` | | Show help message |

---

## Split Selection: Angle Method vs. Legacy Lookahead

By default the pipeline uses **Dr. Ramkumar's angle method** to choose each triangle split. The legacy recursive lookahead engine is retained and selectable with `--legacy-lookahead` so you can compare the two head to head.

### How the angle method works

For each triangle, and for each of its three corners `v` (opposite side `u–w`):

1. The reference direction `v→u` defines angle `0`; the sweep extends to `v→w`.
2. Every boundary endpoint inside the triangle is converted to **one relative angle** (`|atan2(cross, dot)|`), so each polygon segment becomes an interval `[a_min, a_max]` on a 1-D angle line.
3. Candidate cuts are placed exactly at the angles of existing polygon endpoints — aiming the ray *through* an actual vertex, so the cut adds no new boundary point on that segment.
4. For each candidate angle `x`, a stab count gives `left` (segments fully ≤ x), `right` (segments fully ≥ x), and `cuts = n − left − right`.

The cost combines your four minimization heuristics:

- **Heuristic 1 & 2 (fewest crossings / vertex-through cuts):** `cuts` term — a cut through an existing vertex that crosses no other segment costs nothing.
- **Heuristic 3 (balance):** `1.5 · |left − right| / (left + right)` favors balanced left/right children.
- **Heuristic 4 (pruning):** at most 50 candidates on each side of the median angle (≤ 100 total per corner) are evaluated, discarding the rest for speed on full-state polygons.

Once the best corner + target endpoint is chosen, the actual split point is realized with **exact CGAL arithmetic** (`Line(v, endpoint) ∩ Segment(u, w)`), and the existing exact-kernel edge-partitioning routine records the line splits. So the angle method only *selects* the cut; the exact transaction logging is unchanged.

### Lookahead for the angle method (`--lookahead N`)

By default (`--lookahead 0`) the angle method is **greedy**: it realizes the single lowest-cost cut per triangle (the `O(k log k)` fast path). With `--lookahead N > 0` it instead expands its top `ANGLE_LOOKAHEAD_WIDTH` (8) scored cuts, continues each subtree with the angle method `N` levels deep, and keeps the cut that yields the **fewest total leaves** — the same minimize-leaves objective the legacy engine optimizes, now driven by angle cuts. The recursion shares the legacy engine's memoization and recursive-call cap, so it stays safe on full-state polygons. Deeper `N` produces lower-leaf decompositions at increasing runtime cost; e.g. on a reduced Arizona the interior-leaf count fell `4566 → 4437 → 4175` for `N = 0, 1, 2`. The same flag controls the legacy engine, so `--lookahead` now applies uniformly to both methods.

### Cost comparison

| | Per-split selection cost |
|---|---|
| Legacy greedy / lookahead | O(candidates × edges × 2^D) intersection tests |
| **Angle method** | **O(k log k) per triangle** (angles computed once per triangle) |

> **Can angles be computed only once?** Once **per triangle**, yes — the three corner sweeps are reused across all candidate cuts for that triangle. They cannot be reused across splits: when a triangle splits, each child is a new triangle whose corners, reference vectors, and segment subset differ (one corner is the fresh split point), so its angles must be recomputed.

---

## Pipeline Stages

The tool executes up to five stages sequentially:

### Stage 1: Read Shapefile

Reads the `.shp` file using shapelib. If `--state` is specified, reads the companion `.dbf` and skips records whose `STATEFP` field doesn't match; `--county` additionally matches `COUNTYFP`. Uses the first matching polygon record's exterior ring.

**Closed polygon handling**: If the ring's first and last points are identical (standard for shapefiles), the duplicate closing point is detected and stripped automatically. This is logged:

```
[SHAPEFILE] Closed polygon detected: first == last point,
            skipping duplicate (ring had 730 pts, using 729)
```

### Stage 2: Douglas-Peucker Reduction (optional)

If a reduction percentage is supplied, this uses the Ramer–Douglas–Peucker algorithm with binary search over the tolerance (epsilon) to hit the requested reduction percentage. Omit the percentage to keep every original vertex.

A 50% reduction means roughly half the vertices are removed while preserving the polygon's shape. The minimum output is 3 vertices (a valid polygon for CGAL).

```
[REDUCE] 729 -> 365 points (49.9% reduction)
```

### Stage 3: Generate Output Files

Two files are produced (unless disabled with `--no-header` / `--no-points`):

**CGAL Header File** (`polygon_output.h`): A complete C++ header that can be `#include`d directly in any CGAL program:

```cpp
#ifndef POLYGON_OUTPUT_H
#define POLYGON_OUTPUT_H

#include <vector>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Point_2.h>

using K     = CGAL::Exact_predicates_exact_constructions_kernel;
using Point = K::Point_2;

std::vector<Point> simple_polygon = {
    Point(-124.7305087498, 48.3954543602),
    Point(-124.6897583008, 48.3644409180),
    ...
};

#endif // POLYGON_OUTPUT_H
```

**Points File** (`polygon_output.txt`): Plain-text `x y` coordinates, one per line.

### Stage 4: Triangular Decomposition

The polygon is fed into the slopeline-based triangulation engine:

1. A convex hull is computed from the polygon points
2. A minimum enclosing triangle is found using rotating calipers
3. The enclosing triangle is recursively split using slopeline cuts
4. Each cut is chosen by the **angle method** (default) or the legacy lookahead engine
5. Every split (triangle and boundary line) is recorded as an AGS transaction

### Stage 5: Fixed-Point Quantization + Walkable Boundary

Unless `--no-quantize` is given, the pipeline publishes blockchain-ready transactions:

- **Global integer grid.** Longitude is a `uint32` with `0` at the international date line and 360° mapped across 2³²; latitude is an `int32` with `0` at the equator and ±90° across ±(2³¹−1). Resolution is ~0.93 cm (lon, equator) and ~0.47 cm (lat).
- **Quantize once per vertex.** Each distinct boundary vertex is snapped once, so split points shared by two sub-segments receive identical integers — the walk cannot develop a gap.
- **Walkability check.** The fully-subdivided boundary is rebuilt (original edges + every recorded split point, ordered along each edge) and verified to be a single **closed, connected** integer cycle.

```
[QUANTIZE] 13 boundary sub-segments, max snap 0.003 m, walkable=YES
```

Because snapping moves points by up to ~1 grid unit, the on-chain verifier should check that a split point lies **within tolerance** of a segment rather than exactly on it. `quantization_report.txt` reports the max observed deviation and a recommended tolerance.

---

### Stage 6: PerWalk + MarkTriangles (interior/exterior classification)

Every leaf triangle is classified as **interior** (inside the region) or **exterior**, using the two-part scheme from the dissertation:

- **Per-leaf interior test.** Whether a tile lies on the interior (left) side of the CCW boundary is decided by a fast double-precision point-in-polygon test on the tile centroid. This is mathematically the "left of the CCW boundary" test, but it stays correct on real TIGER polygons where the tessellation is **not** edge-conforming to the boundary — i.e. where a boundary vertex lies mid-side of a longer cut rather than at a triangle vertex. (A pure edge-adjacency walk leaks across the boundary in that case; the point-in-polygon test does not.)
- **PerWalk seeds.** Interior leaves that share an edge with an exterior leaf are exactly the tiles the boundary walk touches directly — these are the seeds.
- **MarkTriangles flood.** Deep-interior leaves that the walk never touches are reached by flooding from the seeds across interior–interior shared edges. Walls are the edges where an interior and an exterior leaf meet, so the flood cannot cross the boundary. Any interior leaf the flood cannot reach (e.g. one isolated by a non-conforming adjacency) is still classified correctly by the point-in-polygon test and reported as `unreached`.

Region codes: `1 = INSIDE`, `2 = OUTSIDE`, `0 = BOUNDARY` (reserved). They are named constants (`REGION_INSIDE` / `REGION_OUTSIDE` / `REGION_BOUNDARY`) and trivially re-mapped.

The result is verified two independent ways and logged:

```
[MARK] interior=10002 exterior=19952 (perwalk-seeds 9992 + flood 10) completeness=PASS
```

- **Area completeness** — Σ(interior leaf area) vs. polygon area, judged by a relative tolerance (0.1%). A small residual is expected on real boundaries because the tessellation is not perfectly boundary-conforming: a few leaves have the boundary passing through their interior, so their whole area counts as interior. The check still flags any gross misclassification.
- **Exact spot-check** — a capped, evenly-spaced sample of leaves is re-tested with CGAL's *exact* centroid point-in-polygon predicate; disagreements (should be 0) are logged to `triangle_decision_log.txt`.

Each interior leaf contributes one **MarkTriangles** transaction (Lᴵ total), included in the transaction counts.

**Performance.** The classification is linear-ish in the tessellation and runs in double precision, so it scales to full states: a ~10,000-vertex polygon (≈30,000 leaves) marks and writes all outputs in a few seconds. (An earlier exact-kernel formulation was orders of magnitude slower and could stall on large states before any CSV was written; that is fixed.)

---

---

## Batch Mode: Many States at Once (`--num-states N`)

Pass `--num-states N` to process the **first N states found in the shapefile**
(in record order) in a single invocation. Each state is written into its own
`state_<FIPS>/` folder containing that state's complete set of logs and CSVs
(the same files a single run produces), and a top-level `summary.csv` plus an
on-screen table compares all of them.

```
CGAL_SlopeLine_Integrated.exe tl_2024_us_state.shp --num-states 10
```

This runs the first 10 states; `state_36/`, `state_32/`, `state_02/`, … each get
their own `pipeline_log.txt`, `ags_transactions.txt`, `boundary_walk.csv`,
`triangle_marks.csv`, `triangle_splits_quantized.csv`, and so on. You can point
Part 2 (the prover–verifier) straight at one of them:

```
python merkle_prover_verifier.py --dir state_36 --lon -75.0 --lat 43.0
```

Reduction and engine options still apply to every state, e.g.
`... 50 --num-states 10 --legacy-lookahead --lookahead 2`.

### Summary table

After the batch, a table is printed and written to `summary.csv`:

```
FIPS Name                    OrigPts  UsedPts   Leaves  Interior  Walkable    Method   LA    AGS_Txns   Status
---------------------------------------------------------------------------------------------------
36   New York                  10052     2012     5635      2763       YES     angle    0        5678       OK
32   Nevada                     9538     1909     5327      2430       YES     angle    0        5354       OK
02   Alaska                      287       57      161        73       YES     angle    0         162       OK
```

Columns: state FIPS and name, original vs. used polygon-point counts, total leaf
tiles and interior tiles, perimeter-walk / walkable-chain status, split method
and lookahead depth, total AGS transaction count, and per-state status
(`OK`/`FAILED`). The `summary.csv` carries the same fields (with separate
`perimeter_walk` and `walkable_chain` columns) for spreadsheet analysis.

Single-state runs (`--state <FIPS>` or no filter) are unchanged — they write to
the current directory and create no folder.

### Look-ahead depth sweep (`--lookahead-range A-B`)

Pass `--lookahead-range A-B` to run **each state at every look-ahead depth from
A to B**, gathering results at each depth so you can see how the depth drives the
leaf count down. Each depth's outputs go into its own `la_<D>/` folder, and a
per-state trend table is printed (and folded into `summary.csv`, one row per
state × depth).

```
CGAL_SlopeLine_Integrated.exe tl_2024_us_state.shp --state 04 --lookahead-range 0-2
```

```
FIPS 04  Arizona
   LA    Leaves  Interior    AGS_Txns  Walkable    Runtime    dLeaves      %vs_LA
---------------------------------------------------------------------------------
    0      8964      4566        9023       YES       8.0s          0        0.0%
    1      8940      4437        8973       YES      17.5s        -24       -0.3%
    2      7790      4175        7811       YES     108.1s      -1174      -13.1%
```

`dLeaves` and `%vs_LA` are measured against the lowest depth in the range. It
combines with batch mode — outputs then nest as `state_<FIPS>/la_<D>/`:

```
CGAL_SlopeLine_Integrated.exe tl_2024_us_state.shp --num-states 10 --lookahead-range 0-3
```

The sweep applies to whichever engine is active (angle method by default, or the
legacy engine with `--legacy-lookahead`). Note the cost: each step up in depth is
substantially slower, and you pay it for every state, so wide ranges over many
states run long — `--lookahead-range 0-2` on a handful of states is a practical
starting point. Because look-ahead optimizes a horizon *estimate*, a deeper run
occasionally ties or ticks up by a leaf or two before improving again; the trend
table reports that faithfully.

---

## Output Files

After a run, the current directory contains these files:

### Log Files

| File | Purpose |
|---|---|
| `pipeline_log.txt` | Complete pipeline trace with all coordinates |
| `ags_transactions.txt` | Every TRIANGLE_SPLIT and LINE_SPLIT transaction (exact coordinates) |
| `ags_transactions_quantized.txt` | Quantized transactions + PerWalk listing + walk verification |
| `quantization_report.txt` | Grid parameters, max snap deviation, recommended verifier tolerance, walkability |
| `triangle_decision_log.txt` | Split candidate evaluation details |
| `debug_out_lookahead.txt` | Per-iteration debug output |

### Generated Polygon Files

| File | Purpose |
|---|---|
| `polygon_output.h` | CGAL header with polygon points |
| `polygon_output.txt` | Plain-text coordinates |

### CSV Files (for Visualization and Analysis)

| File | Columns | Use For |
|---|---|---|
| `triangles_decomposition.csv` | id, vertices, type, edge count | Visualizing final triangulation |
| `triangles_edges.csv` | triangle_id, edge endpoints | Mapping polygon edges to triangles |
| `triangle_splits.csv` | parent → split line → 2 children (exact) | Reconstructing the split tree |
| `line_splits.csv` | original edge → split point → 2 segments (exact) | Tracking boundary modifications |
| `triangle_splits_quantized.csv` | same as above, integer grid | On-chain ξT transactions |
| `line_splits_quantized.csv` | same as above, integer grid + snap deviation | On-chain ξP transactions |
| `boundary_walk.csv` | seq, exact + quantized sub-segment endpoints | Verifying / replaying the PerWalk cycle |
| `triangle_marks.csv` | leaf_id, region_code, region, 3 vertices, area | Per-leaf interior/exterior marks from PerWalk + MarkTriangles |

---

## Reconstructing Transactions

The logging is designed so you can fully reconstruct any operation after the fact.

### Reconstructing a Triangle Split

From `triangle_splits.csv`, each row gives the parent triangle, the cut line, and the two child triangles. Draw the parent, draw the split line, and verify it produces the two children.

### Reconstructing a Boundary Line Split

From `line_splits.csv`, each row gives the original edge, the split point, and the two resulting segments.

### Reconstructing the Input Polygon

From `pipeline_log.txt`, the `[POLYGON_REDUCED]` (or `[POLYGON_ORIGINAL]`) section lists every vertex; `[POLYGON_EDGES]` lists every edge; `[ENCLOSING_TRIANGLE]` gives the starting triangle; `[FINAL TRIANGLES]` lists every output triangle with its polygon-edge associations.

### Replaying the Walkable Boundary

`boundary_walk.csv` lists the published boundary in order. Each row's quantized end (`qtx,qty`) equals the next row's quantized start (`qsx,qsy`), and the final row's end equals the first row's start — a closed cycle. `ags_transactions_quantized.txt` repeats this as a human-readable PerWalk (first/last 20 sub-segments) plus the YES/NO walkability verdict.

### AGS Transaction Counting

Per AGS protocol: each triangle split = 1 transaction; each line split = 2 transactions (one for the split itself, one for the PerWalk boundary walk during merge).

---

## State FIPS Reference

Common codes for use with the `--state` option:

| FIPS | State | FIPS | State | FIPS | State |
|------|-------|------|-------|------|-------|
| 01 | Alabama | 18 | Indiana | 35 | New Mexico |
| 02 | Alaska | 19 | Iowa | 36 | New York |
| 04 | Arizona | 20 | Kansas | 37 | North Carolina |
| 05 | Arkansas | 21 | Kentucky | 38 | North Dakota |
| 06 | California | 22 | Louisiana | 39 | Ohio |
| 08 | Colorado | 23 | Maine | 40 | Oklahoma |
| 09 | Connecticut | 24 | Maryland | 41 | Oregon |
| 10 | Delaware | 25 | Massachusetts | 42 | Pennsylvania |
| 11 | District of Columbia | 26 | Michigan | 44 | Rhode Island |
| 12 | Florida | 27 | Minnesota | 45 | South Carolina |
| 13 | Georgia | 28 | Mississippi | 46 | South Dakota |
| 15 | Hawaii | 29 | Missouri | 47 | Tennessee |
| 16 | Idaho | 30 | Montana | 48 | Texas |
| 17 | Illinois | 31 | Nebraska | 49 | Utah |
|    |          | 32 | Nevada | 50 | Vermont |
|    |          | 33 | New Hampshire | 51 | Virginia |
|    |          | 34 | New Jersey | 53 | Washington |
|    |          |    |             | 54 | West Virginia |
|    |          |    |             | 55 | Wisconsin |
|    |          |    |             | 56 | Wyoming |

Full list: https://www.census.gov/library/reference/code-lists/ansi.html

Use `--list-states` to discover what's in any shapefile.

---

## County FIPS Reference

County codes are **3-digit** values within a state, used with `--county` (always together with `--state`). The same 3-digit code is reused across states (e.g. `001` exists in every state), so the state filter is required to disambiguate.

The fastest way to find a county code — and to pick a *lightweight* county for quick pipeline iteration — is `--list-counties`, which prints every county in a state ranked by polygon vertex count:

```
CGAL_SlopeLine_Integrated.exe tl_2024_us_county.shp --state 53 --list-counties
```

```
GEOID   St  County  Name                 Vertices
------  --  ------  -------------------  --------
53013   53  013     Clark                     433
53075   53  075     Columbia                  412
...

Total: 39 counties (state 53)

--- Top 10 counties with fewest polygon vertices ---
  (Fewest vertices = fastest triangulation)
  1. Columbia (GEOID 53075, state 53, county 075) - 412 vertices
  2. Clark (GEOID 53013, state 53, county 013) - 433 vertices
  3. ...
```

The full listing is sorted by vertex count (ascending), and the trailing block repeats the ten lightest counties explicitly.

Counties near the top of the "fewest vertices" list triangulate in a fraction of a second, so they are ideal for testing pipeline changes before committing to a full-state run.

A few commonly used examples:

| State | County FIPS | County |
|---|---|---|
| 53 (WA) | 033 | King |
| 53 (WA) | 061 | Snohomish |
| 06 (CA) | 037 | Los Angeles |
| 06 (CA) | 075 | San Francisco |
| 48 (TX) | 201 | Harris |
| 48 (TX) | 113 | Dallas |

Always confirm with `--list-counties` for the specific shapefile year, since codes occasionally change between vintages.

---

## Examples

### Python: county ring → tessellate → walk

```python
from tridecomp import (
    tessellate, tridecomp, edges_of, gettri, walk, tiles,
    read_tiger_polygon,
)

p = read_tiger_polygon("tl_2025_us_county.zip", state="31", county="039")
t = gettri(p, rat=1.5)
events = tridecomp(p, t)                    # tessellate(edges_of(p), t)
leaves = tiles(events, t)
w = walk(events, t, p)
print(len(events), len(leaves), w["simple"], w["n_inside"])
```

### Python: refine a leaf with extra segments

```python
code = next(iter(leaves))
more = tessellate(extra_segments, leaves[code], start=code)
combined = events + more
```

### Discover what's in a national shapefile

```
CGAL_SlopeLine_Integrated.exe tl_2024_us_county.shp --list-states
```

### List counties in a state, ranked by vertex count

```
CGAL_SlopeLine_Integrated.exe tl_2024_us_county.shp --state 53 --list-counties
```

### A single county at full resolution (angle method, quantized)

```
CGAL_SlopeLine_Integrated.exe tl_2024_us_county.shp --state 53 --county 033
```

### Washington state boundary, 50% reduction

```
CGAL_SlopeLine_Integrated.exe tl_2024_us_state.shp 50 --state 53
```

### A/B compare: legacy lookahead at depth 2 on the same county

```
CGAL_SlopeLine_Integrated.exe tl_2024_us_county.shp --state 06 --county 037 --legacy-lookahead --lookahead 2
```

### Texas with custom output names

```
CGAL_SlopeLine_Integrated.exe tl_2024_us_state.shp 40 --state 48 ^
    --header texas_boundary.h --points texas_points.txt --prefix texas_poly
```

### Skip quantization (geometry-only run)

```
CGAL_SlopeLine_Integrated.exe tl_2024_us_state.shp 50 --state 53 --no-quantize
```

### Light vs. heavy reduction

```
# faithful but slower
CGAL_SlopeLine_Integrated.exe tl_2024_us_state.shp 10 --state 53
# fast iteration
CGAL_SlopeLine_Integrated.exe tl_2024_us_state.shp 90 --state 53
```

---

## Troubleshooting

### "Polygon is not simple after reduction"

The Douglas-Peucker simplification removed too many vertices and caused a self-intersection. Try a smaller reduction percentage, or omit reduction entirely.

### "No matching polygon record found"

Either the FIPS code is wrong or the shapefile doesn't contain polygon geometry for that state/county. Run `--list-states` or `--list-counties` first to verify what's in the file.

### Walkable chain reports NO

The published integer boundary failed the closed/connected check. See `quantization_report.txt` — it reports the first break location. This usually indicates a split point recorded off its parent edge; the exact (`line_splits.csv`) and quantized (`line_splits_quantized.csv`) files let you locate the offending segment.

### "Cannot open shapefile"

Make sure all companion files (`.shp`, `.dbf`, `.shx`) are in the same folder. The `.dbf` file is required for state/county filtering and record names.

### Linker error: unresolved `SHPOpen`

Shapelib isn't linked. Make sure you ran `vcpkg install shapelib:x64-windows` and `vcpkg integrate install`.

### Linker error: unresolved CGAL symbols

CGAL isn't found. Verify `vcpkg install cgal:x64-windows` completed successfully. CGAL is mostly header-only, but some components need Boost and GMP which vcpkg handles automatically.

### Output files appear in unexpected directory

All output files are written to the **current working directory** (where you run the command from), not the shapefile's directory. Use `cd` to set your working directory before running.

### Program ends silently after the decomposition loop (no CSVs) — Windows/MSVC

If the console reaches the end of the `[progress]` lines (queue near 0) and then exits with no `[mark]` lines, no `PIPELINE COMPLETE`, and no CSVs, the process hit a **stack overflow** — historically caused by accumulating tens of thousands of exact (`Lazy_exact_nt`) triangle areas into one expression tree, whose recursive evaluation overflowed MSVC's default 1 MB thread stack (Linux/g++ defaults to 8 MB and masked it). The area sums are now done in plain `double`, which removes the deep tree; rebuild from the current source and it runs to completion (~15–20 s for a full state).


