# TriDecomp

Split an enclosing triangle against constraint segments. The C++ engine does one thing: **`tessellate(segments, t)`**. Walk, marks, enclosing-triangle helpers, and I/O are Python.

A polygon ring is just one source of segments. Refining a leaf later is the same call with extra segments and a non-empty `start` code.

## Repo layout

| Path | Role |
|---|---|
| `tessellate.cpp` | CGAL **EPEC** split engine |
| `tridecomp_api.h` | C++ result types (`TDSplitEvent`, …) |
| `python/tridecomp/_cpp.cpp` | Python/numpy wrapper around the engine |
| `python/tridecomp/` | Public package: `tessellate`, `walk`, `tiles`, `edges_of`, `gettri`, … |
| `Makefile` | Builds the extension (`.so`) |
| `pyproject.toml` | Editable pip install of the `tridecomp` package |

## Install

**Dependencies:** a C++17 compiler, CGAL headers, GMP, MPFR, Python ≥ 3.10, and numpy.

Linux (system packages):

```bash
sudo apt-get install -y g++ make libcgal-dev libgmp-dev libmpfr-dev
```

Without root, `make deps` unpacks Ubuntu `.deb` files into `~/.local/tridecomp-deps`. The Makefile uses `/usr/include/CGAL` when present; otherwise that local prefix.

From the repo root, using the same Python you will import with:

```bash
make python          # PY=... to override (default: $HOME/miniconda3/envs/tridecomp-py/bin/python)
pip install -e . --no-build-isolation --no-deps
```

`make python` compiles `tessellate.cpp` + `_cpp.cpp` into `python/tridecomp/_cpp*.so`. The pip command is **once per environment**: it points `import tridecomp` at this repo. It does not compile C++.

Optional: CGAL’s Python bindings (`engine="python"` / EPIC stand-in) and `pyogrio` (`read_tiger_polygon`).

### Jupyter

Use the **same** interpreter that built the `.so`. Register it as a kernel:

```bash
python -m ipykernel install --user --name tridecomp-py --display-name "tridecomp-py"
```

In the notebook pick that kernel (not a generic “Python 3”). After `make python`, restart the kernel so it loads the new `.so`.

## API

```python
from tridecomp import tessellate, tridecomp, edges_of, gettri, walk, tiles
from tridecomp import UNMARKED, INSIDE, OUTSIDE

p = np.array([[0.0, 0.0], [1.0, 0.0], [1.0, 1.0], [0.0, 1.0]])
t = gettri(p, rat=1.5)
events = tessellate(edges_of(p), t)   # or tridecomp(p, t)
# events = tessellate(edges_of(p), t, engine="python")  # EPIC stand-in
w = walk(events, t, p)
leaves = tiles(events, t)            # {id: triangle (3, 2)}
```

### Segments: shape `(M, 2, 2)`

Each row is one constraint edge — two endpoints `(x, y)`:

```
segments[i] = [[x0, y0],     # endpoint A
               [x1, y1]]     # endpoint B
```

That is the same four numbers, row-major, as shape `(M, 4)`: `x0, y0, x1, y1`. Tessellation treats the edge as undirected. Cycle **direction** matters later in `walk`, not here.

`edges_of(p)` builds this from a closed ring `p` of shape `(N, 2)` (CCW; a repeated closing vertex is dropped).

### `tessellate(segments, t, start="")`

- **`t`**: shape `(3, 2)`. Clockwise input is reversed to CCW.
- **`start`**: tree address of `t`. Empty (`""`) for a fresh enclosing triangle. Pass a leaf code (e.g. `"01011"`) to split that tile further; child ids are `start+"0"` / `start+"1"`.
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

Leaves are implied: a triangle that never appears later as a `parent`. Replay with `tiles(events, t, start="")`.

### Split invariants

- The two children partition the parent (disjoint interiors, union = parent).
- From the directed cut, **`0` = left child**, **`1` = right child**.
- A cut that **crosses** a constraint interior splits it; one stub goes to each child (a `line_splits` record).
- A cut that **contains both endpoints** of a constraint (the segment lies on the new shared side, possibly as a proper subsegment / T-junction) copies that segment into **both** children. No new vertex, no line-split record.

Split selection is the greedy **angle method** (lookahead 0), with a slope-line fallback if no scored cut realizes.

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
| `read_tiger_polygon(path, state=, county=)` | Exterior ring from a TIGER shapefile or zip (`pyogrio`) |

## Rebuilding

Run `make python` only after changing `tessellate.cpp`, `python/tridecomp/_cpp.cpp`, or `tridecomp_api.h`, then restart any notebook kernel. Python-only edits (`walk.py`, notebooks, …) do not need a rebuild.

```bash
make clean    # removes the .so
```
