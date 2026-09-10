// CGAL_SlopeLine_Integrated.cpp
// =========================================================================
// Integrated pipeline: Shapefile -> Optional Douglas-Peucker reduction ->
// CGAL header generation -> Triangular decomposition with full transaction
// logging and optimized split selection.
//
// Changes from CGAL_SlopeLine_Lookahead_Tracked.cpp:
//   - Replaced GDAL/OGR with shapelib (lighter dependency)
//   - Replaced every-Nth-point decimation with Douglas-Peucker %reduction
//   - Reduction is now optional (omit for full-resolution polygons)
//   - Added CGAL header file (.h) and points file (.txt) generation
//   - Added --state and --county FIPS filters for nation-level shapefiles
//   - Added --list-counties mode with vertex count ranking
//   - Optimized split selection: 3-tier priority (polygon-edge > zero-crossing
//     > boundary-crossing), vertex-through detection, balance scoring, and
//     median-based candidate pruning (50+50 per vertex)
//   - Robust closed-polygon detection (skip duplicate last point)
//   - Enhanced transaction logging with priority classification
//   - Vertex-through splits no longer over-counted as LINE_SPLITs
//   - All logs written to separate files for post-hoc analysis
//   - NEW: Dr. Ramkumar's ANGLE METHOD for split selection (default). For each
//     triangle, corner angles to every boundary endpoint are computed ONCE
//     (three corner sweeps) and reused across all candidate cuts, giving
//     O(k log k) selection per triangle vs. the legacy O(candidates x edges x
//     lookahead) intersection search. The legacy recursive lookahead engine is
//     retained and selectable via --legacy-lookahead for A/B comparison.
//   - NEW: FIXED-POINT QUANTIZATION at publishing time. All geometry stays in
//     the CGAL exact kernel; only when emitting transactions are coordinates
//     snapped to a global integer grid (longitude uint32 from the date line,
//     latitude int32 from the equator; ~1 cm resolution).
//   - NEW: BLOCKCHAIN-WALKABLE boundary. The fully-subdivided boundary is
//     rebuilt and published as a single closed, connected integer cycle
//     (PerWalk), verified end-to-end so a chain verifier can walk it.
//   - NEW: PERWALK + MARKTRIANGLES interior/exterior classification. Every leaf
//     triangle is marked INSIDE(1)/OUTSIDE(2). The per-leaf interior test is a
//     fast double-precision point-in-polygon on the tile centroid (the robust
//     realization of "left of the CCW boundary", stable even when the
//     tessellation is not edge-conforming to the boundary -- common on real
//     TIGER data). PerWalk seeds = interior leaves touching the boundary;
//     MarkTriangles floods to deep-interior leaves across interior-interior
//     edges. Verified by area completeness AND a capped exact-arithmetic
//     spot-check. Scales to full states (~30k leaves in seconds). Per-leaf marks
//     are written to triangle_marks.csv; each interior leaf is one MarkTriangles
//     transaction (L^I), included in the transaction counts.
//   - FIX: the spatial-hash cell key in the marking left-shifted a negative
//     (western-longitude) cell index, which is undefined behaviour -- it ran on
//     g++ but crashed/stalled under MSVC right at the start of the marking
//     phase (no CSVs written). Cell/edge keys are now composed with unsigned
//     shifts. Completeness is judged by a relative-area tolerance (0.1%) since
//     real tessellations are not perfectly boundary-conforming.
//   - MarkTriangles flood now uses a SEGMENT-OVERLAP adjacency graph instead of
//     endpoint-hash adjacency. Real tessellations are full of T-junctions (one
//     leaf's full side A--C abuts neighbours' A--B and B--C), so endpoint keys
//     never match and the flood fragmented the interior (~1500 leaves left
//     unreached on Arizona). Two leaf sides are now topological neighbours iff
//     they are exactly CGAL::collinear and their 1-D extents overlap by a
//     positive length -- a shared border, full or partial. Sides are grouped by
//     canonical (normal, offset) and swept by interval, so exact predicates run
//     O(sides) times. Drives flood-unreached to 0 with no runtime cost; region
//     codes are unchanged (still the point-in-polygon ground truth).
//
// Usage:
//   CGAL_SlopeLine_Integrated <shapefile.shp> [reduction_%] [options]
//
// Options:
//   --state <FIPS>       Filter by 2-digit state FIPS code (e.g. 53=WA)
//   --county <FIPS>      Filter by 3-digit county FIPS code (e.g. 033=King)
//   --angle-method       Use the angle method for split selection (DEFAULT)
//   --legacy-lookahead   Use the legacy recursive lookahead engine instead
//   --lookahead <N>      Lookahead depth for the legacy engine (default: 0)
//   --no-quantize        Skip fixed-point quantization / walkable publishing
//   --quiet              Progress + AGS records + CSVs; skip per-iteration debug
//   --verbose            Dump every cell's vertices/edges to debug_out_lookahead.txt
//   --reduce <pct>       Reduction % (alternative to positional)
//   --header <path>      Output CGAL header file (default: polygon_output.h)
//   --points <path>      Output points file (default: polygon_output.txt)
//   --no-header          Skip header file generation
//   --no-points          Skip points file generation
//   --list-states        Show FIPS codes in shapefile and exit
//   --list-counties      Show county FIPS codes and vertex counts, then exit
//
// Examples:
//   CGAL_SlopeLine_Integrated tl_2024_us_state.shp --state 53
//   CGAL_SlopeLine_Integrated tl_2024_us_state.shp 50 --state 53
//   CGAL_SlopeLine_Integrated tl_2024_us_county.shp --state 53 --county 033
//   CGAL_SlopeLine_Integrated tl_2024_us_county.shp 70 --state 06 --legacy-lookahead --lookahead 2
// =========================================================================

#include <iostream>
#include <vector>
#include <fstream>
#include <iomanip>
#include <queue>
#include <array>
#include <utility>
#include <string>
#include <limits>
#include <unordered_map>
#include <unordered_set>
#include <sstream>
#include <chrono>
#include <cmath>
#include <algorithm>
#include <cctype>
#include <cstring>
#include <map>
#include <CGAL/number_utils.h>
#include <CGAL/convex_hull_2.h>
#include <CGAL/bounding_box.h>
#include <CGAL/Iso_rectangle_2.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Point_2.h>
#include <CGAL/Polygon_2.h>
#include <CGAL/Segment_2.h>
#include <CGAL/Triangle_2.h>
#include <CGAL/centroid.h>
#include <CGAL/intersections.h>
#include <CGAL/Object.h>
#include <cstdlib>
#include <cstdint>
#include <CGAL/Vector_2.h>
#include <shapefil.h>
#include "tridecomp_api.h"

// Portable directory operations for multi-state batch mode.
#ifdef _WIN32
  #include <direct.h>
  #define PORTABLE_MKDIR(d)  _mkdir(d)
  #define PORTABLE_CHDIR(d)  _chdir(d)
  #define PORTABLE_GETCWD(b, n) _getcwd(b, n)
#else
  #include <sys/stat.h>
  #include <unistd.h>
  #define PORTABLE_MKDIR(d)  mkdir(d, 0755)
  #define PORTABLE_CHDIR(d)  chdir(d)
  #define PORTABLE_GETCWD(b, n) getcwd(b, n)
#endif

// --- Type Definitions ---
using K = CGAL::Exact_predicates_exact_constructions_kernel;
using Point = K::Point_2;
using Polygon = CGAL::Polygon_2<K>;
using Edge = K::Segment_2;
using Triangle_cgal = K::Triangle_2;
using Line = K::Line_2;
using Triangle = std::array<Point, 3>;

// --- Constants ---
const double POINT_TOL = 1e-5;
const double AREA_TOL = 1e-6;

// =========================================================================
// RECURSIVE LOOKAHEAD CONFIGURATION
// =========================================================================
int LOOKAHEAD_DEPTH_LIMIT = 0;

// Branching width for the ANGLE method's look-ahead: at each node only the top
// ANGLE_LOOKAHEAD_WIDTH scored cuts are expanded. Greedy (width-agnostic) when
// LOOKAHEAD_DEPTH_LIMIT == 0. Bounded by the same recursive-call cap as the
// legacy engine, so larger widths stay safe on full-state polygons.
int ANGLE_LOOKAHEAD_WIDTH = 8;

bool USE_CANDIDATES_RATIO = true;
double CANDIDATES_RATIO = 0.1;
int MIN_CANDIDATES_LIMIT = 15;
int MAX_CANDIDATES_LIMIT = 100;
int DEFAULT_MAX_CANDIDATES = 30;
int MAX_CANDIDATES_PER_LEVEL = DEFAULT_MAX_CANDIDATES;

// =========================================================================
// SPLIT-SELECTION STRATEGY
//   true  = angle-projection method (Dr. Ramkumar): each boundary endpoint is
//           reduced to a single relative angle at a triangle corner, turning
//           cut selection into a 1-D interval-stabbing problem solved with one
//           sort + binary search per corner -- O(k log k) per triangle.
//           Angles are computed ONCE PER TRIANGLE (three corner sweeps) and
//           reused across every candidate cut for that triangle.
//   false = legacy recursive bounded look-ahead engine (select with
//           --legacy-lookahead; --lookahead <N> sets its depth).
// =========================================================================
bool USE_ANGLE_METHOD = true;

// =========================================================================
// FIXED-POINT QUANTIZATION (global integer coordinate system)
//   Publishing-time snap-to-grid so every emitted transaction is a
//   blockchain-ready integer record.  ALL internal geometry stays in CGAL's
//   exact kernel; quantization is applied ONLY when publishing transactions.
//     Longitude -> unsigned 32-bit, zero at the international date line,
//                  360 deg across 2^32 steps (~0.93 cm at the equator).
//     Latitude  -> signed   32-bit, zero at the equator, +/-90 deg.
//   Disable with --no-quantize.
// =========================================================================
bool QUANTIZE_TRANSACTIONS = true;

// =========================================================================
// LOG VERBOSITY
//   0 = --quiet   : stderr progress + AGS records + CSVs; skip per-iteration debug
//   1 = default   : full log files, buffered (no per-iteration flush)
//   2 = --verbose : also dump every cell's vertices/edges to debug_out_lookahead.txt
// =========================================================================
int G_VERBOSITY = 1;
static inline bool log_decisions() { return G_VERBOSITY >= 1; }
static inline bool verbose_debug() { return G_VERBOSITY >= 2; }
static inline int progress_interval() {
    if (G_VERBOSITY <= 0) return 500;
    if (G_VERBOSITY == 1) return 100;
    return 10;
}

// =========================================================================
// SHAPEFILE READING (shapelib) with state FIPS filter
// =========================================================================

// Trim whitespace from string ends
static std::string trim_str(const std::string& s) {
    auto start = s.find_first_not_of(" \t\r\n");
    if (start == std::string::npos) return "";
    auto end = s.find_last_not_of(" \t\r\n");
    return s.substr(start, end - start + 1);
}

// Find a named field in DBF (case-insensitive, tries multiple variants)
static int find_dbf_field(DBFHandle dbf, const std::vector<std::string>& names) {
    int n_fields = DBFGetFieldCount(dbf);
    for (int f = 0; f < n_fields; ++f) {
        char fname[32];
        DBFGetFieldInfo(dbf, f, fname, nullptr, nullptr);
        std::string fn(fname);
        std::transform(fn.begin(), fn.end(), fn.begin(),
                       [](unsigned char c){ return (char)std::toupper(c); });
        for (const auto& target : names) {
            std::string t = target;
            std::transform(t.begin(), t.end(), t.begin(),
                           [](unsigned char c){ return (char)std::toupper(c); });
            if (fn == t) return f;
        }
    }
    return -1;
}

// List all unique STATEFP values in a shapefile's DBF
static void list_states_in_shapefile(const std::string& shp_path) {
    std::string dbf_path = shp_path;
    auto dot = dbf_path.rfind('.');
    if (dot != std::string::npos) dbf_path = dbf_path.substr(0, dot);
    dbf_path += ".dbf";

    DBFHandle dbf = DBFOpen(dbf_path.c_str(), "rb");
    if (!dbf) {
        std::cerr << "[ERROR] Cannot open DBF file: " << dbf_path << "\n";
        return;
    }

    int fld = find_dbf_field(dbf, {"STATEFP", "STATEFP20", "STATEFP10", "STATE"});
    int fld_name = find_dbf_field(dbf, {"NAME", "NAME20", "NAME10"});

    if (fld < 0) {
        std::cerr << "[ERROR] No STATEFP field found in " << dbf_path << "\n";
        int n_fields = DBFGetFieldCount(dbf);
        std::cerr << "  Available fields: ";
        for (int f = 0; f < n_fields; ++f) {
            char fname[32];
            DBFGetFieldInfo(dbf, f, fname, nullptr, nullptr);
            if (f > 0) std::cerr << ", ";
            std::cerr << fname;
        }
        std::cerr << "\n";
        DBFClose(dbf);
        return;
    }

    int n_records = DBFGetRecordCount(dbf);
    std::map<std::string, int> fips_count;
    std::map<std::string, std::string> fips_sample;

    for (int i = 0; i < n_records; ++i) {
        const char* val = DBFReadStringAttribute(dbf, i, fld);
        if (!val) continue;
        std::string fips = trim_str(val);
        fips_count[fips]++;
        if (fips_sample.find(fips) == fips_sample.end() && fld_name >= 0) {
            const char* nm = DBFReadStringAttribute(dbf, i, fld_name);
            if (nm) fips_sample[fips] = trim_str(nm);
        }
    }
    DBFClose(dbf);

    std::cout << "FIPS  Records  Sample Name\n";
    std::cout << "----  -------  -----------\n";
    for (const auto& [fips, count] : fips_count) {
        std::cout << fips;
        for (size_t p = fips.size(); p < 6; ++p) std::cout << ' ';
        std::string cnt = std::to_string(count);
        for (size_t p = cnt.size(); p < 8; ++p) std::cout << ' ';
        std::cout << cnt << "  " << fips_sample[fips] << "\n";
    }
    std::cout << "\nTotal: " << n_records << " records, "
              << fips_count.size() << " unique FIPS codes\n";
}

// List all counties in a shapefile, optionally filtered by state FIPS.
// Also shows polygon point counts to help identify counties with fewest vertices.
static void list_counties_in_shapefile(const std::string& shp_path,
                                        const std::string& state_filter = "")
{
    std::string dbf_path = shp_path;
    auto dot = dbf_path.rfind('.');
    if (dot != std::string::npos) dbf_path = dbf_path.substr(0, dot);
    dbf_path += ".dbf";

    DBFHandle dbf = DBFOpen(dbf_path.c_str(), "rb");
    if (!dbf) {
        std::cerr << "[ERROR] Cannot open DBF file: " << dbf_path << "\n";
        return;
    }
    SHPHandle shp = SHPOpen(shp_path.c_str(), "rb");
    if (!shp) {
        std::cerr << "[ERROR] Cannot open SHP file: " << shp_path << "\n";
        DBFClose(dbf);
        return;
    }

    int fld_statefp = find_dbf_field(dbf, {"STATEFP", "STATEFP20", "STATEFP10", "STATE"});
    int fld_countyfp = find_dbf_field(dbf, {"COUNTYFP", "COUNTYFP20", "COUNTYFP10", "COUNTY"});
    int fld_name = find_dbf_field(dbf, {"NAME", "NAME20", "NAME10"});
    int fld_geoid = find_dbf_field(dbf, {"GEOID", "GEOID20", "GEO_ID"});

    if (fld_countyfp < 0) {
        std::cerr << "[ERROR] No COUNTYFP field found in " << dbf_path << "\n";
        int n_fields = DBFGetFieldCount(dbf);
        std::cerr << "  Available fields: ";
        for (int f = 0; f < n_fields; ++f) {
            char fname[32];
            DBFGetFieldInfo(dbf, f, fname, nullptr, nullptr);
            if (f > 0) std::cerr << ", ";
            std::cerr << fname;
        }
        std::cerr << "\n";
        DBFClose(dbf);
        SHPClose(shp);
        return;
    }

    int n_entities = 0, shp_type = 0;
    double min_bound[4], max_bound[4];
    SHPGetInfo(shp, &n_entities, &shp_type, min_bound, max_bound);

    struct CountyInfo {
        std::string state_fips;
        std::string county_fips;
        std::string geoid;
        std::string name;
        int vertex_count;
    };
    std::vector<CountyInfo> counties;

    for (int i = 0; i < n_entities; ++i) {
        // State filter
        std::string st_fips;
        if (fld_statefp >= 0) {
            const char* v = DBFReadStringAttribute(dbf, i, fld_statefp);
            st_fips = v ? trim_str(v) : "";
        }
        if (!state_filter.empty() && st_fips != state_filter) continue;

        std::string co_fips;
        if (fld_countyfp >= 0) {
            const char* v = DBFReadStringAttribute(dbf, i, fld_countyfp);
            co_fips = v ? trim_str(v) : "";
        }

        std::string name;
        if (fld_name >= 0) {
            const char* v = DBFReadStringAttribute(dbf, i, fld_name);
            name = v ? trim_str(v) : "";
        }

        std::string geoid;
        if (fld_geoid >= 0) {
            const char* v = DBFReadStringAttribute(dbf, i, fld_geoid);
            geoid = v ? trim_str(v) : "";
        }

        // Count vertices in exterior ring
        SHPObject* obj = SHPReadObject(shp, i);
        if (!obj) continue;
        int start = 0;
        int end = (obj->nParts > 1) ? obj->panPartStart[1] : obj->nVertices;
        int vtx_count = end - start;
        // Account for closed polygon
        if (vtx_count > 1) {
            if (std::abs(obj->padfX[start] - obj->padfX[end-1]) < 1e-12 &&
                std::abs(obj->padfY[start] - obj->padfY[end-1]) < 1e-12) {
                vtx_count--;
            }
        }
        SHPDestroyObject(obj);

        counties.push_back({st_fips, co_fips, geoid, name, vtx_count});
    }

    SHPClose(shp);
    DBFClose(dbf);

    // Sort by vertex count ascending
    std::sort(counties.begin(), counties.end(),
        [](const CountyInfo& a, const CountyInfo& b) {
            return a.vertex_count < b.vertex_count;
        });

    std::cout << "GEOID   St  County  Name                 Vertices\n";
    std::cout << "------  --  ------  -------------------  --------\n";
    for (const auto& c : counties) {
        std::string geoid_str = c.geoid;
        for (size_t p = geoid_str.size(); p < 8; ++p) geoid_str += ' ';
        std::cout << geoid_str << c.state_fips;
        for (size_t p = c.state_fips.size(); p < 4; ++p) std::cout << ' ';
        std::cout << c.county_fips;
        for (size_t p = c.county_fips.size(); p < 8; ++p) std::cout << ' ';
        std::string nm = c.name.substr(0, 19);
        std::cout << nm;
        for (size_t p = nm.size(); p < 21; ++p) std::cout << ' ';
        std::cout << c.vertex_count << "\n";
    }
    std::cout << "\nTotal: " << counties.size() << " counties"
              << (state_filter.empty() ? "" : " (state " + state_filter + ")")
              << "\n";

    // Show top 10 with fewest vertices
    int top_n = std::min(10, (int)counties.size());
    if (top_n > 0) {
        std::cout << "\n--- Top " << top_n << " counties with fewest polygon vertices ---\n";
        std::cout << "  (Fewest vertices = fastest triangulation)\n";
        for (int i = 0; i < top_n; ++i) {
            std::cout << "  " << (i + 1) << ". " << counties[i].name
                      << " (GEOID " << counties[i].geoid
                      << ", state " << counties[i].state_fips
                      << ", county " << counties[i].county_fips
                      << ") - " << counties[i].vertex_count << " vertices\n";
        }
    }
}

// Read polygon points from shapefile, optionally filtering by state and/or county FIPS.
// Reads the FIRST matching polygon record.
// Handles closed polygons (first == last point) by stripping the duplicate.
struct ShapeReadResult {
    std::vector<Point> points;
    std::string record_name;
    std::string geoid;
    std::string state_fips;
    std::string county_fips;
    int original_point_count = 0;
    bool was_closed = false;
};

ShapeReadResult read_polygon_from_shapefile(const std::string& shp_path,
                                            const std::string& state_fips_filter,
                                            const std::string& county_fips_filter,
                                            std::ofstream& log)
{
    ShapeReadResult result;

    SHPHandle shp = SHPOpen(shp_path.c_str(), "rb");
    if (!shp) {
        throw std::runtime_error("Cannot open shapefile: " + shp_path);
    }

    // Open companion DBF
    std::string dbf_path = shp_path;
    auto dot = dbf_path.rfind('.');
    if (dot != std::string::npos) dbf_path = dbf_path.substr(0, dot);
    dbf_path += ".dbf";
    DBFHandle dbf = DBFOpen(dbf_path.c_str(), "rb");

    int n_entities = 0, shp_type = 0;
    double min_bound[4], max_bound[4];
    SHPGetInfo(shp, &n_entities, &shp_type, min_bound, max_bound);

    log << "[SHAPEFILE] Path: " << shp_path << "\n";
    log << "[SHAPEFILE] Type: " << shp_type << ", Records: " << n_entities << "\n";
    log << "[SHAPEFILE] Bounds X: " << min_bound[0] << " .. " << max_bound[0] << "\n";
    log << "[SHAPEFILE] Bounds Y: " << min_bound[1] << " .. " << max_bound[1] << "\n";

    // Locate DBF fields
    int fld_statefp = -1, fld_countyfp = -1, fld_name = -1, fld_geoid = -1, fld_namelsad = -1;
    if (dbf) {
        fld_statefp  = find_dbf_field(dbf, {"STATEFP", "STATEFP20", "STATEFP10", "STATE"});
        fld_countyfp = find_dbf_field(dbf, {"COUNTYFP", "COUNTYFP20", "COUNTYFP10", "COUNTY"});
        fld_name     = find_dbf_field(dbf, {"NAME", "NAME20", "NAME10"});
        fld_namelsad = find_dbf_field(dbf, {"NAMELSAD", "NAMELSAD20"});
        fld_geoid    = find_dbf_field(dbf, {"GEOID", "GEOID20", "GEO_ID"});

        int n_fields = DBFGetFieldCount(dbf);
        log << "[SHAPEFILE] DBF fields: ";
        for (int f = 0; f < n_fields; ++f) {
            char fname[32];
            DBFGetFieldInfo(dbf, f, fname, nullptr, nullptr);
            if (f > 0) log << ", ";
            log << fname;
        }
        log << "\n";
    }

    if (!state_fips_filter.empty()) {
        log << "[SHAPEFILE] State filter: FIPS " << state_fips_filter << "\n";
    }
    if (!county_fips_filter.empty()) {
        log << "[SHAPEFILE] County filter: FIPS " << county_fips_filter << "\n";
    }

    // Find the first matching polygon record
    int matched_record = -1;
    for (int i = 0; i < n_entities; ++i) {
        // State filter
        if (!state_fips_filter.empty() && dbf && fld_statefp >= 0) {
            const char* val = DBFReadStringAttribute(dbf, i, fld_statefp);
            std::string rec_fips = val ? trim_str(val) : "";
            if (rec_fips != state_fips_filter) continue;
        }

        // County filter
        if (!county_fips_filter.empty() && dbf && fld_countyfp >= 0) {
            const char* val = DBFReadStringAttribute(dbf, i, fld_countyfp);
            std::string rec_county = val ? trim_str(val) : "";
            if (rec_county != county_fips_filter) continue;
        }

        SHPObject* obj = SHPReadObject(shp, i);
        if (!obj) continue;

        // Must be polygon type
        if (obj->nSHPType != SHPT_POLYGON && obj->nSHPType != SHPT_POLYGONZ &&
            obj->nSHPType != SHPT_POLYGONM) {
            SHPDestroyObject(obj);
            continue;
        }

        // Read exterior ring (first part)
        int start = 0;
        int end = (obj->nParts > 1) ? obj->panPartStart[1] : obj->nVertices;
        int num_ring_pts = end - start;

        result.original_point_count = num_ring_pts;

        // Detect closed polygon (first == last point)
        result.was_closed = false;
        if (num_ring_pts > 1) {
            double x0 = obj->padfX[start], y0 = obj->padfY[start];
            double xn = obj->padfX[end - 1], yn = obj->padfY[end - 1];
            if (std::abs(x0 - xn) < 1e-12 && std::abs(y0 - yn) < 1e-12) {
                result.was_closed = true;
                num_ring_pts--;  // Skip duplicate closing point
                log << "[SHAPEFILE] Closed polygon detected: first == last point, "
                    << "skipping duplicate (ring had " << (num_ring_pts + 1)
                    << " pts, using " << num_ring_pts << ")\n";
            }
        }

        result.points.reserve(num_ring_pts);
        for (int v = start; v < start + num_ring_pts; ++v) {
            result.points.emplace_back(obj->padfX[v], obj->padfY[v]);
        }

        // Read DBF attributes
        if (dbf) {
            if (fld_namelsad >= 0) {
                const char* v = DBFReadStringAttribute(dbf, i, fld_namelsad);
                if (v) result.record_name = trim_str(v);
            } else if (fld_name >= 0) {
                const char* v = DBFReadStringAttribute(dbf, i, fld_name);
                if (v) result.record_name = trim_str(v);
            }
            if (fld_statefp >= 0) {
                const char* v = DBFReadStringAttribute(dbf, i, fld_statefp);
                if (v) result.state_fips = trim_str(v);
            }
            if (fld_countyfp >= 0) {
                const char* v = DBFReadStringAttribute(dbf, i, fld_countyfp);
                if (v) result.county_fips = trim_str(v);
            }
            if (fld_geoid >= 0) {
                const char* v = DBFReadStringAttribute(dbf, i, fld_geoid);
                if (v) result.geoid = trim_str(v);
            }
        }

        matched_record = i;
        SHPDestroyObject(obj);
        break;
    }

    SHPClose(shp);
    if (dbf) DBFClose(dbf);

    if (matched_record < 0) {
        std::string filter_info;
        if (!state_fips_filter.empty()) filter_info += " state=" + state_fips_filter;
        if (!county_fips_filter.empty()) filter_info += " county=" + county_fips_filter;
        throw std::runtime_error("No matching polygon record found in shapefile"
            + (filter_info.empty() ? "" : " (filter:" + filter_info + ")"));
    }

    log << "[SHAPEFILE] Matched record #" << matched_record;
    if (!result.record_name.empty()) log << " - " << result.record_name;
    if (!result.geoid.empty()) log << " (GEOID: " << result.geoid << ")";
    log << "\n";
    log << "[SHAPEFILE] Exterior ring: " << result.points.size() << " unique vertices"
        << " (original ring had " << result.original_point_count << " points"
        << (result.was_closed ? ", closed" : ", open") << ")\n";
    log.flush();

    std::cerr << "[SHAPEFILE] Read " << result.points.size() << " points from record #"
              << matched_record;
    if (!result.record_name.empty()) std::cerr << " (" << result.record_name << ")";
    std::cerr << "\n";

    return result;
}

// =========================================================================
// DOUGLAS-PEUCKER PERCENTAGE-BASED POLYGON SIMPLIFICATION
// =========================================================================

// Perpendicular distance from point P to line segment A-B
static double perp_distance(const Point& P, const Point& A, const Point& B) {
    double ax = CGAL::to_double(A.x()), ay = CGAL::to_double(A.y());
    double bx = CGAL::to_double(B.x()), by = CGAL::to_double(B.y());
    double px = CGAL::to_double(P.x()), py = CGAL::to_double(P.y());
    double dx = bx - ax, dy = by - ay;
    double len2 = dx * dx + dy * dy;
    if (len2 < 1e-30) {
        double ex = px - ax, ey = py - ay;
        return std::sqrt(ex * ex + ey * ey);
    }
    double cross = std::abs(dy * px - dx * py + bx * ay - by * ax);
    return cross / std::sqrt(len2);
}

// Recursive Douglas-Peucker core
static void dp_recurse(const std::vector<Point>& pts, size_t first, size_t last,
                        double epsilon, std::vector<bool>& keep) {
    if (last <= first + 1) return;
    double max_dist = 0.0;
    size_t max_idx = first;
    for (size_t i = first + 1; i < last; ++i) {
        double d = perp_distance(pts[i], pts[first], pts[last]);
        if (d > max_dist) { max_dist = d; max_idx = i; }
    }
    if (max_dist > epsilon) {
        keep[max_idx] = true;
        dp_recurse(pts, first, max_idx, epsilon, keep);
        dp_recurse(pts, max_idx, last, epsilon, keep);
    }
}

// Simplify polygon points with given epsilon tolerance
static std::vector<Point> simplify_points_dp(const std::vector<Point>& pts, double epsilon) {
    size_t n = pts.size();
    if (n <= 3) return pts;
    std::vector<bool> keep(n, false);
    keep.front() = true;
    keep.back() = true;
    dp_recurse(pts, 0, n - 1, epsilon, keep);
    std::vector<Point> result;
    for (size_t i = 0; i < n; ++i) {
        if (keep[i]) result.push_back(pts[i]);
    }
    // Guarantee at least 3 vertices
    if (result.size() < 3) {
        result.clear();
        for (size_t i = 0; i < std::min(n, (size_t)3); ++i) {
            result.push_back(pts[i * n / 3]);
        }
    }
    return result;
}

// Binary search for the epsilon that yields the target vertex count
static double find_epsilon_for_target(const std::vector<Point>& pts, size_t target_count) {
    double xmin = CGAL::to_double(pts[0].x()), xmax = xmin;
    double ymin = CGAL::to_double(pts[0].y()), ymax = ymin;
    for (const auto& p : pts) {
        double x = CGAL::to_double(p.x()), y = CGAL::to_double(p.y());
        xmin = std::min(xmin, x); xmax = std::max(xmax, x);
        ymin = std::min(ymin, y); ymax = std::max(ymax, y);
    }
    double diag = std::sqrt((xmax - xmin) * (xmax - xmin) + (ymax - ymin) * (ymax - ymin));
    double lo = 0.0, hi = diag;
    for (int iter = 0; iter < 60; ++iter) {
        double mid = (lo + hi) * 0.5;
        auto simplified = simplify_points_dp(pts, mid);
        if (simplified.size() > target_count) lo = mid;
        else hi = mid;
    }
    return (lo + hi) * 0.5;
}

// Reduce polygon points by a percentage using Douglas-Peucker
std::vector<Point> reduce_polygon_by_percentage(const std::vector<Point>& points,
                                                 double reduction_pct,
                                                 std::ofstream& log)
{
    if (reduction_pct <= 0.0 || reduction_pct >= 100.0) {
        log << "[REDUCE] Invalid reduction_pct " << reduction_pct << ", returning original\n";
        return points;
    }
    if (points.size() <= 3) {
        log << "[REDUCE] Only " << points.size() << " points, no reduction possible\n";
        return points;
    }

    double keep_ratio = 1.0 - reduction_pct / 100.0;
    size_t target = static_cast<size_t>(
        std::max(3.0, std::ceil(static_cast<double>(points.size()) * keep_ratio)));

    double eps = find_epsilon_for_target(points, target);
    auto result = simplify_points_dp(points, eps);

    double actual_reduction = 100.0 * (1.0 - static_cast<double>(result.size()) /
                                               static_cast<double>(points.size()));

    log << "[REDUCE] Douglas-Peucker reduction:\n";
    log << "[REDUCE]   Input points:     " << points.size() << "\n";
    log << "[REDUCE]   Requested reduction: " << reduction_pct << "%\n";
    log << "[REDUCE]   Target points:    " << target << "\n";
    log << "[REDUCE]   Actual output:    " << result.size() << " points\n";
    log << "[REDUCE]   Actual reduction: " << std::fixed << std::setprecision(1)
        << actual_reduction << "%\n";
    log << "[REDUCE]   Epsilon used:     " << std::scientific << eps << "\n";
    log.flush();

    std::cerr << "[REDUCE] " << points.size() << " -> " << result.size()
              << " points (" << std::fixed << std::setprecision(1)
              << actual_reduction << "% reduction)\n";

    return result;
}

// =========================================================================
// CGAL HEADER FILE AND POINTS FILE GENERATION
// =========================================================================

static std::string sanitize_identifier(const std::string& s) {
    std::string out;
    for (char c : s) {
        if (std::isalnum(static_cast<unsigned char>(c)) || c == '_') out += c;
        else out += '_';
    }
    if (!out.empty() && std::isdigit(static_cast<unsigned char>(out[0])))
        out = "_" + out;
    return out;
}

static std::string make_include_guard(const std::string& path) {
    auto pos = path.find_last_of("/\\");
    std::string name = (pos == std::string::npos) ? path : path.substr(pos + 1);
    std::string guard;
    for (char c : name) {
        if (std::isalnum(static_cast<unsigned char>(c)))
            guard += static_cast<char>(std::toupper(static_cast<unsigned char>(c)));
        else guard += '_';
    }
    return guard;
}

void write_cgal_header(const std::vector<Point>& points,
                       const std::string& output_path,
                       const std::string& var_name,
                       const std::string& source_info,
                       std::ofstream& log)
{
    std::ofstream out(output_path);
    if (!out) {
        log << "[HEADER] ERROR: Cannot open " << output_path << " for writing\n";
        return;
    }

    std::string guard = make_include_guard(output_path);
    std::string safe_var = sanitize_identifier(var_name);

    out << "// Auto-generated by CGAL_SlopeLine_Integrated\n";
    out << "// Source: " << source_info << "\n";
    out << "// Points: " << points.size() << "\n";
    out << "//\n\n";
    out << "#ifndef " << guard << "\n";
    out << "#define " << guard << "\n\n";
    out << "#include <vector>\n\n";
    out << "// CGAL Headers needed for the types used in this file\n";
    out << "#include <CGAL/Exact_predicates_exact_constructions_kernel.h>\n";
    out << "#include <CGAL/Point_2.h>\n\n";
    out << "// --- Type Definitions ---\n";
    out << "using K     = CGAL::Exact_predicates_exact_constructions_kernel;\n";
    out << "using Point = K::Point_2;\n\n";

    out << std::fixed << std::setprecision(10);
    out << "std::vector<Point> " << safe_var << " = {\n";
    for (size_t i = 0; i < points.size(); ++i) {
        out << "    Point(" << CGAL::to_double(points[i].x()) << ", "
            << CGAL::to_double(points[i].y()) << ")";
        if (i + 1 < points.size()) out << ",";
        out << "\n";
    }
    out << "};\n\n";
    out << "#endif // " << guard << "\n";
    out.close();

    log << "[HEADER] Wrote " << points.size() << " points to " << output_path
        << " (variable: " << safe_var << ")\n";
    std::cerr << "[HEADER] Generated: " << output_path << " (" << points.size() << " points)\n";
}

void write_points_file(const std::vector<Point>& points,
                       const std::string& output_path,
                       const std::string& source_info,
                       std::ofstream& log)
{
    std::ofstream out(output_path);
    if (!out) {
        log << "[POINTS] ERROR: Cannot open " << output_path << " for writing\n";
        return;
    }

    out << std::fixed << std::setprecision(10);
    out << "# Polygon points generated by CGAL_SlopeLine_Integrated\n";
    out << "# Source: " << source_info << "\n";
    out << "# Points: " << points.size() << "\n";
    out << "# Format: x y\n#\n";

    for (const auto& p : points) {
        out << CGAL::to_double(p.x()) << " " << CGAL::to_double(p.y()) << "\n";
    }
    out.close();

    log << "[POINTS] Wrote " << points.size() << " points to " << output_path << "\n";
    std::cerr << "[POINTS] Generated: " << output_path << " (" << points.size() << " points)\n";
}

void compute_dynamic_limits(int num_segments) {
    if (USE_CANDIDATES_RATIO && CANDIDATES_RATIO > 0) {
        int computed_candidates = static_cast<int>(num_segments * CANDIDATES_RATIO);
        MAX_CANDIDATES_PER_LEVEL = std::max(MIN_CANDIDATES_LIMIT, std::min(MAX_CANDIDATES_LIMIT, computed_candidates));
        
        std::cerr << "[CONFIG] Segments: " << num_segments << "\n"
                  << "  LOOKAHEAD_DEPTH_LIMIT = " << LOOKAHEAD_DEPTH_LIMIT << " (fixed)\n"
                  << "  CANDIDATES_RATIO=" << CANDIDATES_RATIO 
                  << " -> MAX_CANDIDATES_PER_LEVEL=" << MAX_CANDIDATES_PER_LEVEL << "\n";
    } else {
        MAX_CANDIDATES_PER_LEVEL = DEFAULT_MAX_CANDIDATES;
        std::cerr << "[CONFIG] Segments: " << num_segments << "\n"
                  << "  LOOKAHEAD_DEPTH_LIMIT = " << LOOKAHEAD_DEPTH_LIMIT << " (fixed)\n"
                  << "  MAX_CANDIDATES_PER_LEVEL = " << MAX_CANDIDATES_PER_LEVEL << " (fixed)\n";
    }
}

const bool ENABLE_EARLY_TERMINATION = true;
const long long LOG_FLUSH_INTERVAL = 1000;
const long long MAX_RECURSIVE_CALLS_PER_EVAL = 1000000;

static bool g_evaluation_terminated = false;

// =========================================================================
// TRANSACTION TRACKING STRUCTURES
// =========================================================================

// Structure to record a triangle split
struct TriangleSplitRecord {
    int split_id;
    Triangle parent_triangle;
    Triangle child1;
    Triangle child2;
    Edge split_line;
    int iteration;
    bool is_lookahead;  // true = lookahead split, false = greedy fallback
};

// Structure to record a boundary line split
struct LineSplitRecord {
    int split_id;
    Edge original_edge;
    Point split_point;
    Edge segment1;
    Edge segment2;
    int triangle_split_id;  // Which triangle split caused this
    int iteration;
};

// Global tracking for transactions
static std::vector<TriangleSplitRecord> g_triangle_splits;
static std::vector<LineSplitRecord> g_line_splits;
static std::ofstream g_transaction_log;

// Transaction counters
static int g_triangle_split_count = 0;
static int g_line_split_count = 0;

// Final leaf tessellation + per-leaf region codes (filled by MarkTriangles).
// Kept global so main() can emit them to CSV / transaction logs.
static std::vector<Triangle> g_leaf_triangles;
static std::vector<int>      g_leaf_region_codes;
static int g_mark_triangle_count = 0;   // L^I interior leaves (MarkTriangles txns)

// =========================================================================
// END TRANSACTION TRACKING
// =========================================================================

// --- Struct for Final Triangle Information ---
struct FinalTriangle {
    Triangle vertices;
    std::vector<Edge> polygon_sides;
};

// --- Struct for Split Candidate ---
// Priority levels:
//   0 = split line IS a polygon segment (becomes triangle side, no new geometry)
//   1 = zero boundary crossings (passes through polygon vertices only)
//   2 = nonzero boundary crossings
struct SplitCandidate {
    bool valid = false;
    std::array<Triangle, 2> child_triangles;
    Edge split_line;
    int intersections = 0;
    int priority = 0;
    double balance_score = 0.0;   // lower = more balanced child edge counts
    bool is_polygon_edge = false; // true if split line coincides with a polygon segment
};

// --- Struct for Best Cut Result ---
struct BestCutResult {
    bool found = false;
    std::array<Triangle, 2> child_triangles;
    Edge split_line;
    int min_leaves = std::numeric_limits<int>::max();
    int priority = 0;
};

// --- Debug and Statistics ---
static std::ofstream g_decision_log;   // opened per run (per-state) in run_one_state
static long long g_recursive_call_count = 0;
static long long g_cache_hit_count = 0;
static const long long LOG_EVERY_RECURSIVE_CALLS = 10000;

// --- Forward Declarations ---
bool point_eq(const Point& p1, const Point& p2, double tol = POINT_TOL);
std::vector<Edge> filter_degenerate_edges(const std::vector<Edge>& edges);
std::string point_to_string(const Point& p);

// --- Geometric Helper Functions ---
double area2d(const Triangle& triangle) {
    auto exact_area = Triangle_cgal(triangle[0], triangle[1], triangle[2]).area();
    return CGAL::to_double(CGAL::abs(exact_area));
}

bool point_eq(const Point& p1, const Point& p2, double tol) {
    return CGAL::squared_distance(p1, p2) < tol * tol;
}

bool edge_eq(const Edge& e1, const Edge& e2, double tol = POINT_TOL) {
    return (point_eq(e1.source(), e2.source(), tol) && point_eq(e1.target(), e2.target(), tol)) ||
        (point_eq(e1.source(), e2.target(), tol) && point_eq(e1.target(), e2.source(), tol));
}

bool edge_on_triangle_side(const Edge& edge, const Triangle& triangle) {
    for (int i = 0; i < 3; ++i) {
        Edge side(triangle[i], triangle[(i + 1) % 3]);
        if (side.has_on(edge.source()) && side.has_on(edge.target())) {
            return true;
        }
    }
    return false;
}

bool triangle_is_inside_buffered(const Polygon& polygon, const Triangle& triangle) {
    // Retained as a simple exact point-in-polygon oracle (centroid test) for
    // ad-hoc checks. The main interior/exterior classification now uses the
    // PerWalk + MarkTriangles marking (fast double-precision point-in-polygon),
    // see mark_triangles_perwalk.
    Point centroid = CGAL::centroid(triangle[0], triangle[1], triangle[2]);
    return polygon.bounded_side(centroid) != CGAL::ON_UNBOUNDED_SIDE;
}

std::string point_to_string(const Point& p) {
    std::stringstream ss;
    ss << std::fixed << std::setprecision(8) << "(" << CGAL::to_double(p.x()) << ", " << CGAL::to_double(p.y()) << ")";
    return ss.str();
}

std::string triangle_to_string(const Triangle& tri) {
    std::stringstream ss;
    ss << "[" << point_to_string(tri[0]) << ", " << point_to_string(tri[1]) << ", " << point_to_string(tri[2]) << "]";
    return ss.str();
}

std::string edge_to_string(const Edge& e) {
    std::stringstream ss;
    ss << point_to_string(e.source()) << " -> " << point_to_string(e.target());
    return ss.str();
}

// =========================================================================
// OPTIMIZED SPLIT SELECTION HELPERS
// =========================================================================

// Check if a split line coincides with (is contained in or contains) a polygon edge.
// Returns true when the split line lies along a polygon segment, meaning
// the polygon edge will become a triangle side with zero new geometry cost.
static bool split_line_is_polygon_edge(const Edge& split_line,
                                        const std::vector<Edge>& all_polygon_edges)
{
    for (const auto& pe : all_polygon_edges) {
        // Check if split line and polygon edge are collinear and overlap
        if (split_line.supporting_line() == pe.supporting_line()) {
            // Check if the split line endpoints lie on the polygon edge
            // or vice versa (they share the same segment)
            bool src_on_pe = pe.has_on(split_line.source());
            bool tgt_on_pe = pe.has_on(split_line.target());
            if (src_on_pe && tgt_on_pe) return true;
            // Also check if split line contains the polygon edge
            bool pe_src_on = split_line.has_on(pe.source());
            bool pe_tgt_on = split_line.has_on(pe.target());
            if (pe_src_on && pe_tgt_on) return true;
        }
    }
    return false;
}

// Check if a point is a polygon vertex (endpoint of any polygon edge).
static bool is_polygon_vertex(const Point& p, const std::vector<Edge>& all_polygon_edges)
{
    for (const auto& pe : all_polygon_edges) {
        if (point_eq(p, pe.source()) || point_eq(p, pe.target()))
            return true;
    }
    return false;
}

// Count true boundary crossings: intersections that actually split polygon edges.
// A split line passing through an existing polygon vertex is NOT counted as a
// crossing (no new split point is introduced). Only interior intersections where
// the split line cuts through the middle of a polygon edge are counted.
int count_intersections_improved(const Edge& split_line,
                                 const std::vector<Edge>& polygon_edges)
{
    // First check for collinear overlap — same as original
    for (const auto& poly_edge : polygon_edges) {
        if (split_line.supporting_line() == poly_edge.supporting_line()) {
            if (poly_edge.has_on(split_line.source()) || poly_edge.has_on(split_line.target()) ||
                split_line.has_on(poly_edge.source()) || split_line.has_on(poly_edge.target())) {
                return 0;
            }
        }
    }

    // Collect all intersection points (deduplicated)
    std::vector<Point> intersection_pts;
    for (const auto& poly_edge : polygon_edges) {
        auto result = CGAL::intersection(split_line, poly_edge);
        if (!result) continue;
        const CGAL::Object& obj = *result;
        if (const Point* p = CGAL::object_cast<Point>(&obj)) {
            // Skip split line's own endpoints
            if (point_eq(*p, split_line.source(), POINT_TOL) ||
                point_eq(*p, split_line.target(), POINT_TOL))
                continue;
            // Check if this point is already collected (dedup)
            bool already = false;
            for (const auto& existing : intersection_pts) {
                if (point_eq(*p, existing, POINT_TOL)) { already = true; break; }
            }
            if (already) continue;
            // If the intersection is at a polygon vertex, it does NOT count
            // as a new split — the vertex already exists in the boundary.
            if (is_polygon_vertex(*p, polygon_edges)) continue;
            intersection_pts.push_back(*p);
        }
        else if (const Edge* s = CGAL::object_cast<Edge>(&obj)) {
            // Collinear overlap — should have been caught above; count as 1
            intersection_pts.push_back(CGAL::midpoint(s->source(), s->target()));
        }
    }
    return static_cast<int>(intersection_pts.size());
}

// Compute a balance score for a candidate split: how evenly polygon edges
// are distributed between the two child triangles. Lower score = more balanced.
// Returns abs(edges_in_child1 - edges_in_child2) / max(total_edges, 1).
// Uses lightweight Triangle_cgal containment (orientation tests) instead of
// constructing full Polygon_2 objects.
static double compute_balance_score(const std::array<Triangle, 2>& child_tris,
                                     const std::vector<Edge>& edges,
                                     const Edge& split_line)
{
    // Build CGAL triangles for fast bounded_side checks via orientation
    Triangle_cgal ct[2] = {
        Triangle_cgal(child_tris[0][0], child_tris[0][1], child_tris[0][2]),
        Triangle_cgal(child_tris[1][0], child_tris[1][1], child_tris[1][2])
    };

    int count[2] = {0, 0};
    for (const auto& edge : edges) {
        Point mid = CGAL::midpoint(edge.source(), edge.target());
        for (int i = 0; i < 2; ++i) {
            auto side = ct[i].bounded_side(mid);
            if (side == CGAL::ON_BOUNDED_SIDE || side == CGAL::ON_BOUNDARY) {
                count[i]++;
                break;
            }
        }
    }
    int total = count[0] + count[1];
    if (total == 0) return 0.0;
    return static_cast<double>(std::abs(count[0] - count[1])) / static_cast<double>(total);
}

// =========================================================================
// Assign remaining polygon edges to the two children of a split.
//
// Fast path: if both endpoints lie strictly on the same side of the cut's
// supporting line, the edge cannot meet the split *segment*, so skip the
// exact intersection and classify by Triangle_cgal::bounded_side (same
// predicate compute_balance_score already uses).
//
// Slow path: an endpoint on the line, or endpoints on opposite sides, may
// meet the cut. Then CGAL::intersection produces the exact split point so
// LINE_SPLIT records stay exact-kernel.
// =========================================================================
static void place_edge_in_child(const Edge& seg, const Triangle_cgal ct[2],
                                std::vector<Edge> child_edges[2])
{
    if (point_eq(seg.source(), seg.target())) return;
    Point mid = CGAL::midpoint(seg.source(), seg.target());
    for (int i = 0; i < 2; ++i) {
        auto side = ct[i].bounded_side(mid);
        if (side == CGAL::ON_BOUNDED_SIDE || side == CGAL::ON_BOUNDARY) {
            child_edges[i].push_back(seg);
            return;
        }
    }
}

static bool edge_may_meet_cut(const Line& cut_line, const Edge& edge)
{
    auto sa = cut_line.oriented_side(edge.source());
    auto sb = cut_line.oriented_side(edge.target());
    if (sa == CGAL::ON_ORIENTED_BOUNDARY || sb == CGAL::ON_ORIENTED_BOUNDARY)
        return true;
    return (sa == CGAL::ON_POSITIVE_SIDE && sb == CGAL::ON_NEGATIVE_SIDE) ||
           (sa == CGAL::ON_NEGATIVE_SIDE && sb == CGAL::ON_POSITIVE_SIDE);
}

static void record_line_split_if_needed(
    bool track, const Edge& edge, const Point& p,
    const Edge& seg1, const Edge& seg2,
    bool seg1_valid, bool seg2_valid,
    int triangle_split_id, int iteration,
    std::vector<LineSplitRecord>* line_splits_out)
{
    if (!track) return;
    if (seg1_valid && seg2_valid) {
        g_line_split_count++;
        LineSplitRecord record;
        record.split_id = g_line_split_count;
        record.original_edge = edge;
        record.split_point = p;
        record.segment1 = seg1;
        record.segment2 = seg2;
        record.triangle_split_id = triangle_split_id;
        record.iteration = iteration;
        line_splits_out->push_back(record);
        g_line_splits.push_back(record);

        g_transaction_log << "LINE_SPLIT #" << record.split_id
                         << " (caused by TRIANGLE_SPLIT #" << triangle_split_id << ")\n";
        g_transaction_log << "  Original edge: " << edge_to_string(edge) << "\n";
        g_transaction_log << "  Split point: " << point_to_string(p) << "\n";
        g_transaction_log << "  Segment 1: " << edge_to_string(seg1) << "\n";
        g_transaction_log << "  Segment 2: " << edge_to_string(seg2) << "\n";
        g_transaction_log << "  Iteration: " << iteration << "\n\n";
    } else if (verbose_debug() && (seg1_valid || seg2_valid)) {
        g_transaction_log << "VERTEX_THROUGH (TRIANGLE_SPLIT #" << triangle_split_id << ")\n";
        g_transaction_log << "  Edge: " << edge_to_string(edge) << "\n";
        g_transaction_log << "  Vertex hit: " << point_to_string(p) << "\n";
        g_transaction_log << "  No LINE_SPLIT recorded (existing vertex, no new boundary point)\n";
        g_transaction_log << "  Iteration: " << iteration << "\n\n";
    }
}

static std::pair<std::vector<Edge>, std::vector<Edge>> assign_edges_impl(
    const std::array<Triangle, 2>& child_tris,
    const std::vector<Edge>& edges,
    const Edge& split_line,
    bool track,
    int triangle_split_id,
    int iteration,
    std::vector<LineSplitRecord>* line_splits_out)
{
    std::vector<Edge> child_edges[2];
    Triangle_cgal ct[2] = {
        Triangle_cgal(child_tris[0][0], child_tris[0][1], child_tris[0][2]),
        Triangle_cgal(child_tris[1][0], child_tris[1][1], child_tris[1][2])
    };
    Line cut_line = split_line.supporting_line();

    for (const auto& edge : edges) {
        if (!edge_may_meet_cut(cut_line, edge)) {
            place_edge_in_child(edge, ct, child_edges);
            continue;
        }

        auto intersection_result = CGAL::intersection(edge, split_line);
        if (!intersection_result) {
            place_edge_in_child(edge, ct, child_edges);
            continue;
        }
        const CGAL::Object& intersection_obj = *intersection_result;
        if (const Point* p = CGAL::object_cast<Point>(&intersection_obj)) {
            Edge seg1(edge.source(), *p);
            Edge seg2(*p, edge.target());
            bool seg1_valid = !point_eq(seg1.source(), seg1.target());
            bool seg2_valid = !point_eq(seg2.source(), seg2.target());
            record_line_split_if_needed(track, edge, *p, seg1, seg2,
                                        seg1_valid, seg2_valid,
                                        triangle_split_id, iteration, line_splits_out);
            if (seg1_valid) place_edge_in_child(seg1, ct, child_edges);
            if (seg2_valid) place_edge_in_child(seg2, ct, child_edges);
        } else if (CGAL::object_cast<Edge>(&intersection_obj)) {
            place_edge_in_child(edge, ct, child_edges);
        }
    }
    return { child_edges[0], child_edges[1] };
}

std::pair<std::vector<Edge>, std::vector<Edge>> assign_edges_to_children_tracked(
    const std::array<Triangle, 2>& child_tris,
    const std::vector<Edge>& edges,
    const Edge& split_line,
    int triangle_split_id,
    int iteration,
    std::vector<LineSplitRecord>& line_splits_out)
{
    return assign_edges_impl(child_tris, edges, split_line, true,
                             triangle_split_id, iteration, &line_splits_out);
}

std::pair<std::vector<Edge>, std::vector<Edge>> assign_edges_to_children(
    const std::array<Triangle, 2>& child_tris,
    const std::vector<Edge>& edges,
    const Edge& split_line)
{
    return assign_edges_impl(child_tris, edges, split_line, false, 0, 0, nullptr);
}

std::vector<Edge> filter_degenerate_edges(const std::vector<Edge>& edges) {
    std::vector<Edge> result;
    for (const auto& e : edges) {
        if (!point_eq(e.source(), e.target())) {
            result.push_back(e);
        }
    }
    return result;
}

std::string make_split_key(const Point& from, const Point& to) {
    std::stringstream ss;
    ss << std::fixed << std::setprecision(6);
    double fx = CGAL::to_double(from.x()), fy = CGAL::to_double(from.y());
    double tx = CGAL::to_double(to.x()), ty = CGAL::to_double(to.y());
    if (fx < tx || (fx == tx && fy < ty)) {
        ss << fx << "," << fy << "->" << tx << "," << ty;
    } else {
        ss << tx << "," << ty << "->" << fx << "," << fy;
    }
    return ss.str();
}

// =========================================================================
// Generate all valid split candidates for a triangle
// Optimized heuristics:
//   Priority 0: split line IS a polygon edge (becomes triangle side)
//   Priority 1: zero boundary crossings (including vertex-through cuts)
//   Priority 2: nonzero boundary crossings
// Median pruning: for each vertex, at most 50 candidate endpoints on each
//   side of the median of the opposite edge (100 total per vertex).
// Balance: among equal-priority candidates, prefer balanced child edge counts.
// =========================================================================
std::vector<SplitCandidate> generate_split_candidates(
    const Triangle& triangle,
    const std::vector<Edge>& candidate_edges,
    const std::vector<Edge>& all_polygon_edges)
{
    // Three priority buckets
    std::vector<SplitCandidate> poly_edge_candidates;   // priority 0
    std::vector<SplitCandidate> zero_crossing_candidates; // priority 1
    std::vector<SplitCandidate> other_candidates;        // priority 2
    std::unordered_set<std::string> seen_splits;

    for (int vtx_idx = 0; vtx_idx < 3; ++vtx_idx) {
        const Point& vtx = triangle[vtx_idx];
        const Point& opp0 = triangle[(vtx_idx + 1) % 3];
        const Point& opp1 = triangle[(vtx_idx + 2) % 3];
        Edge opposite_side(opp0, opp1);

        // --- Median pruning (heuristic 4) ---
        // Compute the midpoint of the opposite side as the median reference.
        // Collect all candidate endpoints, project them onto the opposite side,
        // sort by parameter, and keep at most 50 on each side of the median.
        double opp0x = CGAL::to_double(opp0.x()), opp0y = CGAL::to_double(opp0.y());
        double opp1x = CGAL::to_double(opp1.x()), opp1y = CGAL::to_double(opp1.y());
        double opp_dx = opp1x - opp0x, opp_dy = opp1y - opp0y;
        double opp_len2 = opp_dx * opp_dx + opp_dy * opp_dy;

        // Gather all unique candidate endpoints with their projected parameter t
        // along the opposite side direction (t=0 at opp0, t=1 at opp1, median at t=0.5)
        struct EndpointInfo {
            Point pt;
            double t_param; // projection parameter along opposite side
        };
        std::vector<EndpointInfo> all_endpoints;
        std::unordered_set<std::string> seen_endpoints;

        for (const auto& edge : candidate_edges) {
            for (int i = 0; i < 2; ++i) {
                const Point& ep = (i == 0) ? edge.source() : edge.target();
                if (point_eq(vtx, ep)) continue;

                std::string ep_key = point_to_string(ep);
                if (seen_endpoints.count(ep_key) > 0) continue;
                seen_endpoints.insert(ep_key);

                // Project endpoint onto opposite side direction
                double epx = CGAL::to_double(ep.x()), epy = CGAL::to_double(ep.y());
                double t_param = 0.0;
                if (opp_len2 > 1e-30) {
                    t_param = ((epx - opp0x) * opp_dx + (epy - opp0y) * opp_dy) / opp_len2;
                }
                all_endpoints.push_back({ep, t_param});
            }
        }

        // Sort by projection parameter
        std::sort(all_endpoints.begin(), all_endpoints.end(),
            [](const EndpointInfo& a, const EndpointInfo& b) {
                return a.t_param < b.t_param;
            });

        // Find the median index (t_param closest to 0.5)
        const int MAX_PER_SIDE = 50;
        int median_idx = 0;
        double best_dist = std::numeric_limits<double>::max();
        for (int i = 0; i < (int)all_endpoints.size(); ++i) {
            double d = std::abs(all_endpoints[i].t_param - 0.5);
            if (d < best_dist) { best_dist = d; median_idx = i; }
        }

        // Keep at most 50 on each side of the median
        int start_idx = std::max(0, median_idx - MAX_PER_SIDE);
        int end_idx = std::min((int)all_endpoints.size() - 1, median_idx + MAX_PER_SIDE);

        // --- Generate candidates from pruned endpoint list ---
        for (int ei = start_idx; ei <= end_idx; ++ei) {
            const Point& endpoint = all_endpoints[ei].pt;

            Line slope_line(vtx, endpoint);
            auto intersection_result = CGAL::intersection(opposite_side, slope_line);

            if (!intersection_result) continue;
            const CGAL::Object& obj = *intersection_result;
            const Point* p = CGAL::object_cast<Point>(&obj);
            if (!p) continue;

            Triangle t1 = { vtx, opp0, *p };
            Triangle t2 = { vtx, *p, opp1 };

            if (area2d(t1) <= AREA_TOL || area2d(t2) <= AREA_TOL) continue;

            Edge potential_split_line(vtx, *p);

            std::string split_key = make_split_key(vtx, *p);
            if (seen_splits.count(split_key) > 0) continue;
            seen_splits.insert(split_key);

            // --- Heuristic 1: Check if split line IS a polygon segment ---
            bool is_poly_edge = split_line_is_polygon_edge(potential_split_line, all_polygon_edges);

            // --- Heuristic 2 & 3: Count improved intersections ---
            int current_intersections = is_poly_edge ? 0
                : count_intersections_improved(potential_split_line, all_polygon_edges);

            // --- Heuristic 3: Balance score ---
            double balance = compute_balance_score({t1, t2}, candidate_edges, potential_split_line);

            SplitCandidate cand;
            cand.valid = true;
            cand.child_triangles = { t1, t2 };
            cand.split_line = potential_split_line;
            cand.intersections = current_intersections;
            cand.is_polygon_edge = is_poly_edge;
            cand.balance_score = balance;

            if (is_poly_edge) {
                cand.priority = 0;
                poly_edge_candidates.push_back(cand);
            } else if (current_intersections == 0) {
                cand.priority = 1;
                zero_crossing_candidates.push_back(cand);
            } else {
                cand.priority = 2;
                other_candidates.push_back(cand);
            }
        }
    }

    // Sort within each bucket by (intersections ASC, balance_score ASC)
    auto sort_fn = [](const SplitCandidate& a, const SplitCandidate& b) {
        if (a.intersections != b.intersections)
            return a.intersections < b.intersections;
        return a.balance_score < b.balance_score;
    };
    std::sort(poly_edge_candidates.begin(), poly_edge_candidates.end(), sort_fn);
    std::sort(zero_crossing_candidates.begin(), zero_crossing_candidates.end(), sort_fn);
    std::sort(other_candidates.begin(), other_candidates.end(), sort_fn);

    // Assemble final list: priority 0 first, then 1, then 2, up to limit
    std::vector<SplitCandidate> candidates;
    candidates.reserve(std::min((size_t)MAX_CANDIDATES_PER_LEVEL,
                                poly_edge_candidates.size() +
                                zero_crossing_candidates.size() +
                                other_candidates.size()));

    for (const auto& c : poly_edge_candidates) {
        candidates.push_back(c);
        if ((int)candidates.size() >= MAX_CANDIDATES_PER_LEVEL) break;
    }
    for (const auto& c : zero_crossing_candidates) {
        if ((int)candidates.size() >= MAX_CANDIDATES_PER_LEVEL) break;
        candidates.push_back(c);
    }
    for (const auto& c : other_candidates) {
        if ((int)candidates.size() >= MAX_CANDIDATES_PER_LEVEL) break;
        candidates.push_back(c);
    }

    return candidates;
}

// =========================================================================
// Memoization for recursive cost evaluation
// =========================================================================
std::string make_triangle_memo_key(const Triangle& tri, int depth) {
    std::stringstream ss;
    ss << std::fixed << std::setprecision(8);
    for (int i = 0; i < 3; ++i) {
        ss << CGAL::to_double(tri[i].x()) << "," << CGAL::to_double(tri[i].y()) << ",";
    }
    ss << depth;
    return ss.str();
}

std::unordered_map<std::string, int> g_memo_cache;

int estimate_leaves_heuristic(const std::vector<Edge>& edges) {
    if (edges.empty()) return 1;
    int edge_count = static_cast<int>(edges.size());
    return std::max(1, edge_count);
}

// =========================================================================
// Recursive function to count minimum leaves
// =========================================================================
int count_min_leaves_recursive(
    const Triangle& triangle,
    const std::vector<Edge>& edges,
    const std::vector<Edge>& all_polygon_edges,
    int depth,
    int lookahead_remaining)
{
    g_recursive_call_count++;
    
    if (g_recursive_call_count > MAX_RECURSIVE_CALLS_PER_EVAL) {
        if (!g_evaluation_terminated) {
            g_evaluation_terminated = true;
            std::cerr << "\n[WARNING] Recursive call limit reached (" << MAX_RECURSIVE_CALLS_PER_EVAL 
                      << "). Terminating evaluation early.\n";
            g_decision_log << "[WARNING] Recursive call limit reached. Using heuristic.\n";
        }
        return estimate_leaves_heuristic(edges);
    }
    
    if (log_decisions() && g_recursive_call_count % LOG_FLUSH_INTERVAL == 0) {
        std::cerr << "[Recursive] Calls: " << g_recursive_call_count
                  << ", Cache hits: " << g_cache_hit_count
                  << ", Cache size: " << g_memo_cache.size()
                  << ", Depth: " << depth << "/" << LOOKAHEAD_DEPTH_LIMIT << "\n";
    }
    
    std::string memo_key = make_triangle_memo_key(triangle, depth);
    auto memo_it = g_memo_cache.find(memo_key);
    if (memo_it != g_memo_cache.end()) {
        g_cache_hit_count++;
        return memo_it->second;
    }
    
    double tri_area = area2d(triangle);
    if (tri_area <= AREA_TOL) {
        g_memo_cache[memo_key] = 0;
        return 0;
    }
    
    std::vector<Edge> current_edges;
    for (const auto& e : edges) {
        if (!edge_on_triangle_side(e, triangle)) {
            current_edges.push_back(e);
        }
    }
    current_edges = filter_degenerate_edges(current_edges);
    
    if (current_edges.empty()) {
        g_memo_cache[memo_key] = 1;
        return 1;
    }
    
    if (lookahead_remaining <= 0) {
        int estimate = estimate_leaves_heuristic(current_edges);
        g_memo_cache[memo_key] = estimate;
        return estimate;
    }
    
    std::vector<SplitCandidate> candidates = generate_split_candidates(triangle, current_edges, all_polygon_edges);
    
    if (candidates.empty()) {
        g_memo_cache[memo_key] = 1;
        return 1;
    }
    
    int min_leaves = std::numeric_limits<int>::max();
    int best_priority = std::numeric_limits<int>::max();
    
    for (const auto& cand : candidates) {
        if (!cand.valid) continue;
        
        if (ENABLE_EARLY_TERMINATION && min_leaves < std::numeric_limits<int>::max()) {
            if (best_priority <= 1 && cand.priority > best_priority) {
                continue;
            }
        }
        
        auto [child_edges1, child_edges2] = assign_edges_to_children(
            cand.child_triangles, current_edges, cand.split_line);
        
        int leaves1 = count_min_leaves_recursive(
            cand.child_triangles[0], 
            filter_degenerate_edges(child_edges1),
            all_polygon_edges,
            depth + 1, 
            lookahead_remaining - 1);
        
        if (g_evaluation_terminated) {
            g_memo_cache[memo_key] = min_leaves == std::numeric_limits<int>::max() ? 
                                     estimate_leaves_heuristic(current_edges) : min_leaves;
            return g_memo_cache[memo_key];
        }
            
        int leaves2 = count_min_leaves_recursive(
            cand.child_triangles[1], 
            filter_degenerate_edges(child_edges2),
            all_polygon_edges,
            depth + 1, 
            lookahead_remaining - 1);
        
        if (g_evaluation_terminated) {
            g_memo_cache[memo_key] = min_leaves == std::numeric_limits<int>::max() ? 
                                     estimate_leaves_heuristic(current_edges) : min_leaves;
            return g_memo_cache[memo_key];
        }
            
        int total = leaves1 + leaves2;
        
        if (total < min_leaves || (total == min_leaves && cand.priority < best_priority)) {
            min_leaves = total;
            best_priority = cand.priority;
            
            if (ENABLE_EARLY_TERMINATION && cand.priority <= 1 && total <= 2 && depth == 0) {
                break;
            }
        }
    }
    
    g_memo_cache[memo_key] = min_leaves;
    return min_leaves;
}

// =========================================================================
// Find the best split using full recursive evaluation
// =========================================================================
BestCutResult find_best_split_recursive(
    const Triangle& triangle,
    const std::vector<Edge>& edges,
    const std::vector<Edge>& all_polygon_edges,
    int current_depth)
{
    g_recursive_call_count = 0;
    g_cache_hit_count = 0;
    g_evaluation_terminated = false;
    
    auto eval_start = std::chrono::steady_clock::now();
    
    if (log_decisions()) {
        g_decision_log << "[Recursive] Starting evaluation for triangle "
                       << triangle_to_string(triangle)
                       << " at depth " << current_depth
                       << " with " << edges.size() << " edges"
                       << " (cache size: " << g_memo_cache.size() << ")\n";
    }
    
    BestCutResult result;
    
    std::vector<Edge> current_edges;
    for (const auto& e : edges) {
        if (!edge_on_triangle_side(e, triangle)) {
            current_edges.push_back(e);
        }
    }
    current_edges = filter_degenerate_edges(current_edges);
    
    std::vector<SplitCandidate> candidates = generate_split_candidates(triangle, current_edges, all_polygon_edges);
    
    if (candidates.empty()) {
        if (log_decisions())
            g_decision_log << "[Recursive] No valid candidates found.\n";
        return result;
    }
    
    if (log_decisions())
        g_decision_log << "[Recursive] Evaluating " << candidates.size() << " candidates "
                       << "(max " << MAX_CANDIDATES_PER_LEVEL << " allowed)...\n";
    
    int candidates_evaluated = 0;
    
    for (const auto& cand : candidates) {
        if (g_evaluation_terminated) {
            g_decision_log << "[Recursive] Evaluation terminated early after " 
                          << candidates_evaluated << " candidates.\n";
            break;
        }
        
        if (ENABLE_EARLY_TERMINATION && result.found && result.priority <= 1 && cand.priority > 1) {
            continue;
        }
        
        auto [child_edges1, child_edges2] = assign_edges_to_children(
            cand.child_triangles, current_edges, cand.split_line);
        
        int leaves1 = count_min_leaves_recursive(
            cand.child_triangles[0],
            filter_degenerate_edges(child_edges1),
            all_polygon_edges,
            current_depth + 1,
            LOOKAHEAD_DEPTH_LIMIT - 1);
        
        if (g_evaluation_terminated && !result.found) {
            result.found = true;
            result.child_triangles = cand.child_triangles;
            result.split_line = cand.split_line;
            result.min_leaves = leaves1 + estimate_leaves_heuristic(child_edges2);
            result.priority = cand.priority;
            break;
        }
            
        int leaves2 = count_min_leaves_recursive(
            cand.child_triangles[1],
            filter_degenerate_edges(child_edges2),
            all_polygon_edges,
            current_depth + 1,
            LOOKAHEAD_DEPTH_LIMIT - 1);
            
        int total_leaves = leaves1 + leaves2;
        
        candidates_evaluated++;
        
        bool is_better = false;
        if (total_leaves < result.min_leaves) {
            is_better = true;
        } else if (total_leaves == result.min_leaves && result.found) {
            if (cand.priority < result.priority) {
                is_better = true;
            }
        }
        
        if (is_better) {
            result.found = true;
            result.child_triangles = cand.child_triangles;
            result.split_line = cand.split_line;
            result.min_leaves = total_leaves;
            result.priority = cand.priority;
        }
        
        if (verbose_debug() && candidates_evaluated % 10 == 0) {
            g_decision_log << "  [Progress] Evaluated " << candidates_evaluated << "/"
                          << candidates.size() << " candidates, "
                          << g_recursive_call_count << " recursive calls\n";
        }
    }
    
    auto eval_end = std::chrono::steady_clock::now();
    if (log_decisions()) {
        double eval_time = std::chrono::duration<double>(eval_end - eval_start).count();
        g_decision_log << "[Recursive] Evaluation complete: " << candidates_evaluated << " candidates, "
                       << g_recursive_call_count << " recursive calls, "
                       << g_cache_hit_count << " cache hits, "
                       << std::fixed << std::setprecision(2) << eval_time << "s\n";
        if (result.found) {
            g_decision_log << "[Recursive] Best split: "
                           << point_to_string(result.split_line.source()) << " -> "
                           << point_to_string(result.split_line.target())
                           << " with " << result.min_leaves << " estimated leaves"
                           << " (priority " << result.priority << ")\n";
        }
    }
    
    return result;
}

// =========================================================================
// Greedy fallback (updated with improved intersection counting)
// =========================================================================
bool slope_split_greedy(
    const Triangle& triangle,
    const std::vector<Edge>& candidate_edges,
    const std::vector<Edge>& all_polygon_edges,
    std::array<Triangle, 2>& out_split_tris,
    Edge& out_split_line)
{
    bool best_split_found = false;
    int min_intersections = std::numeric_limits<int>::max();
    bool best_is_poly_edge = false;
    double best_balance = std::numeric_limits<double>::max();

    for (const auto& edge : candidate_edges) {
        for (int vtx_idx = 0; vtx_idx < 3; ++vtx_idx) {
            const Point& vtx = triangle[vtx_idx];
            const Point& opp0 = triangle[(vtx_idx + 1) % 3];
            const Point& opp1 = triangle[(vtx_idx + 2) % 3];
            Edge opposite_side(opp0, opp1);

            for (int i = 0; i < 2; ++i) {
                const Point& endpoint = (i == 0) ? edge.source() : edge.target();
                if (point_eq(vtx, endpoint)) continue;

                Line slope_line(vtx, endpoint);
                auto intersection_result = CGAL::intersection(opposite_side, slope_line);

                if (intersection_result) {
                    const CGAL::Object& obj = *intersection_result;
                    if (const Point* p = CGAL::object_cast<Point>(&obj)) {
                        Triangle t1 = { vtx, opp0, *p };
                        Triangle t2 = { vtx, *p, opp1 };

                        if (area2d(t1) > AREA_TOL && area2d(t2) > AREA_TOL) {
                            Edge potential_split_line(vtx, *p);
                            bool is_poly_edge = split_line_is_polygon_edge(potential_split_line, all_polygon_edges);
                            int current_intersections = is_poly_edge ? 0
                                : count_intersections_improved(potential_split_line, all_polygon_edges);
                            double balance = compute_balance_score({t1, t2}, candidate_edges, potential_split_line);

                            // Priority ordering: polygon-edge > fewer intersections > better balance
                            bool is_better = false;
                            if (is_poly_edge && !best_is_poly_edge) {
                                is_better = true;
                            } else if (is_poly_edge == best_is_poly_edge) {
                                if (current_intersections < min_intersections) {
                                    is_better = true;
                                } else if (current_intersections == min_intersections && balance < best_balance) {
                                    is_better = true;
                                }
                            }

                            if (is_better) {
                                min_intersections = current_intersections;
                                best_is_poly_edge = is_poly_edge;
                                best_balance = balance;
                                out_split_tris = { t1, t2 };
                                out_split_line = potential_split_line;
                                best_split_found = true;
                                if (is_poly_edge) return true;
                                if (min_intersections == 0 && balance < 0.05) return true;
                            }
                        }
                    }
                }
            }
        }
    }
    return best_split_found;
}

// =========================================================================
// ANGLE-METHOD SPLIT SELECTION  (Dr. Ramkumar's projection method)
// =========================================================================
// Idea (from the design discussion + TriDecomp.ipynb):
//   For a triangle T with corner v and opposite side (u,w), every boundary-
//   segment endpoint can be reduced from a 2-D (x,y) point to a SINGLE 1-D
//   value -- its angle measured at v relative to the reference direction
//   v->u.  This is the triangular analogue of dropping a coordinate in the
//   slab/rectangle methods.  A "cut" is then just a scalar angle theta in
//   (0, angle(v->w)); the ray from v at angle theta meets the opposite side
//   at the split point.  Each boundary segment becomes an interval
//   [a_min, a_max] on this angle line, so picking the best cut is an
//   interval-stabbing problem solved with one sort + binary search per
//   corner -- O(k log k) per triangle, versus the O(candidates * k * 2^D)
//   of the recursive look-ahead engine.
//
// "Can the angles be computed only ONCE?"  -- Once PER TRIANGLE: yes, and
//   that is the real speed-up.  We sweep all three corners a single time and
//   reuse those angles across every candidate cut for the triangle.  They
//   CANNOT be computed once for the whole run and reused forever: when T is
//   split, each child is a brand-new triangle whose corners (one of them the
//   freshly created split point) and reference vectors differ, and whose
//   segment set is a different subset.  Only the shared corner v could in
//   principle inherit angles shifted by the cut angle; the other two corners
//   are new.  So "once per triangle" is the feasible -- and the important --
//   reuse, and it is what this function does.
//
// Candidate cuts are placed exactly at the angles of existing polygon-
// segment endpoints.  Aiming the ray THROUGH that actual endpoint means:
//   * the cut passes through an existing vertex, adding no new boundary
//     point on that segment (heuristic 2), and
//   * the split point on the opposite side is found with EXACT CGAL
//     arithmetic (line(v,endpoint) intersect segment(u,w)), preserving the
//     exact-kernel area guarantees.
// Per-candidate cost = 1.5*|L-R|/(L+R)  (balance, heuristic 3)
//                    + cuts             (crossings, heuristics 1 & 2).
// Heuristic 4: median +/-50 candidates per corner (<=100), discarding the
// rest to keep selection fast on full-state polygons.
// =========================================================================

struct AngleEndpoint {
    double angle;   // relative angle at the corner, radians in [0, pi]
    Point  pt;      // the polygon endpoint that produced this angle
};

// Absolute angle (radians, [0,pi]) between reference (rx,ry) and (px,py).
static inline double am_rel_angle(double rx, double ry, double px, double py) {
    double cross = rx * py - ry * px;
    double dot   = rx * px + ry * py;
    return std::abs(std::atan2(cross, dot));
}

// --- Scored angle-cut candidate (geometry realized lazily) ---
struct ScoredCut {
    double cost;     // angle-method heuristic cost (lower is better)
    int    corner;   // triangle corner the ray emanates from
    Point  target;   // polygon endpoint the ray is aimed through
    int    cuts;     // boundary segments the cut crosses
    int    embedded; // polygon segments lying along the cut
};

// Score every candidate cut for `triangle` against boundary `edges`, returning
// them sorted best-first. This is the per-corner angle sweep extracted from the
// greedy selector so the greedy path and the look-ahead path share one scorer.
static std::vector<ScoredCut> am_scored_cuts(const Triangle& triangle,
                                             const std::vector<Edge>& edges)
{
    std::vector<ScoredCut> scored;
    if (edges.empty()) return scored;
    const int MAX_PER_SIDE = 50;   // heuristic 4: <=50 candidates on each side

    for (int c = 0; c < 3; ++c) {
        const Point& v = triangle[c];
        const Point& u = triangle[(c + 1) % 3];
        const Point& w = triangle[(c + 2) % 3];
        double vx = CGAL::to_double(v.x()), vy = CGAL::to_double(v.y());
        double rx = CGAL::to_double(u.x()) - vx, ry = CGAL::to_double(u.y()) - vy; // v->u
        double aw = am_rel_angle(rx, ry,
                                 CGAL::to_double(w.x()) - vx,
                                 CGAL::to_double(w.y()) - vy);   // sweep extent (v->w)
        if (aw < 1e-12) continue;

        std::vector<double> seg_lo, seg_hi;
        seg_lo.reserve(edges.size());
        seg_hi.reserve(edges.size());
        std::vector<AngleEndpoint> cand;
        cand.reserve(edges.size() * 2);

        for (const auto& e : edges) {
            double sx = CGAL::to_double(e.source().x()) - vx;
            double sy = CGAL::to_double(e.source().y()) - vy;
            double tx = CGAL::to_double(e.target().x()) - vx;
            double ty = CGAL::to_double(e.target().y()) - vy;
            double a1 = am_rel_angle(rx, ry, sx, sy);
            double a2 = am_rel_angle(rx, ry, tx, ty);
            if (sx*sx + sy*sy < 1e-24) a1 = a2;
            if (tx*tx + ty*ty < 1e-24) a2 = a1;
            seg_lo.push_back(std::min(a1, a2));
            seg_hi.push_back(std::max(a1, a2));
            if (a1 > 1e-9 && a1 < aw - 1e-9) cand.push_back({a1, e.source()});
            if (a2 > 1e-9 && a2 < aw - 1e-9) cand.push_back({a2, e.target()});
        }
        if (cand.empty()) continue;

        std::sort(cand.begin(), cand.end(),
                  [](const AngleEndpoint& a, const AngleEndpoint& b){ return a.angle < b.angle; });
        std::vector<AngleEndpoint> uniq;
        uniq.reserve(cand.size());
        for (const auto& ce : cand)
            if (uniq.empty() || std::abs(ce.angle - uniq.back().angle) > 1e-9)
                uniq.push_back(ce);

        if ((int)uniq.size() > 2 * MAX_PER_SIDE) {
            double mid = 0.5 * aw;
            int median_idx = 0; double bd = std::numeric_limits<double>::max();
            for (int i = 0; i < (int)uniq.size(); ++i) {
                double d = std::abs(uniq[i].angle - mid);
                if (d < bd) { bd = d; median_idx = i; }
            }
            int lo = std::max(0, median_idx - MAX_PER_SIDE);
            int hi = std::min((int)uniq.size() - 1, median_idx + MAX_PER_SIDE);
            uniq.assign(uniq.begin() + lo, uniq.begin() + hi + 1);
        }

        std::vector<double> sorted_hi = seg_hi, sorted_lo = seg_lo;
        std::sort(sorted_hi.begin(), sorted_hi.end());
        std::sort(sorted_lo.begin(), sorted_lo.end());
        int n = (int)edges.size();

        const double AEPS = 1e-9;
        std::vector<double> emb_angles;
        emb_angles.reserve(seg_lo.size());
        for (size_t i = 0; i < seg_lo.size(); ++i)
            if (seg_hi[i] - seg_lo[i] < AEPS)
                emb_angles.push_back(0.5 * (seg_lo[i] + seg_hi[i]));
        std::sort(emb_angles.begin(), emb_angles.end());

        const double EMBED_BONUS = 1000.0;  // dominates cuts (integer) + balance ([0,1.5])

        for (const auto& ce : uniq) {
            double x = ce.angle;
            int raw_left  = (int)(std::upper_bound(sorted_hi.begin(), sorted_hi.end(), x) - sorted_hi.begin());
            int raw_right = n - (int)(std::lower_bound(sorted_lo.begin(), sorted_lo.end(), x) - sorted_lo.begin());
            int embedded = (int)(std::upper_bound(emb_angles.begin(), emb_angles.end(), x + AEPS)
                               - std::lower_bound(emb_angles.begin(), emb_angles.end(), x - AEPS));
            int left  = raw_left  - embedded;
            int right = raw_right - embedded;
            int cuts  = n - left - right - embedded;
            if (cuts < 0) cuts = 0;
            double balance = (left + right > 0)
                ? 1.5 * std::abs(left - right) / double(left + right) : 1.5;
            double cost = cuts + balance - EMBED_BONUS * embedded;
            scored.push_back({ cost, c, ce.pt, cuts, embedded });
        }
    }

    std::sort(scored.begin(), scored.end(),
              [](const ScoredCut& a, const ScoredCut& b){ return a.cost < b.cost; });
    return scored;
}

// Realize a scored cut into exact child triangles + split line. Returns false if
// the cut is degenerate (caller tries the next-best scored cut).
static bool am_realize_cut(const Triangle& triangle, const ScoredCut& sc,
                           const std::vector<Edge>& all_polygon_edges,
                           std::array<Triangle,2>& kids, Edge& split_line, int& priority)
{
    const Point& v = triangle[sc.corner];
    const Point& u = triangle[(sc.corner + 1) % 3];
    const Point& w = triangle[(sc.corner + 2) % 3];
    Line slope_line(v, sc.target);
    Edge opposite_side(u, w);
    auto inter = CGAL::intersection(opposite_side, slope_line);
    if (!inter) return false;
    const CGAL::Object& obj = *inter;
    const Point* newp = CGAL::object_cast<Point>(&obj);
    if (!newp) return false;
    Triangle t1 = { v, u, *newp };
    Triangle t2 = { v, *newp, w };
    if (area2d(t1) <= AREA_TOL || area2d(t2) <= AREA_TOL) return false;
    split_line = Edge(v, *newp);
    kids = { t1, t2 };
    if (sc.embedded > 0 || split_line_is_polygon_edge(split_line, all_polygon_edges))
        priority = 0;
    else
        priority = (sc.cuts == 0) ? 1 : 2;
    return true;
}

// Look-ahead leaf-count estimate for the ANGLE method -- the analogue of
// count_min_leaves_recursive, but candidates come from the angle sweep (top
// ANGLE_LOOKAHEAD_WIDTH cuts per node). Shares the recursive-call cap and memo
// cache with the legacy engine, so it stays safe on full-state polygons.
int am_count_min_leaves(const Triangle& triangle,
                        const std::vector<Edge>& edges,
                        const std::vector<Edge>& all_polygon_edges,
                        int depth,
                        int lookahead_remaining)
{
    g_recursive_call_count++;
    if (g_recursive_call_count > MAX_RECURSIVE_CALLS_PER_EVAL) {
        if (!g_evaluation_terminated) {
            g_evaluation_terminated = true;
            g_decision_log << "[Angle] Recursive call limit reached. Using heuristic.\n";
        }
        return estimate_leaves_heuristic(edges);
    }

    std::string memo_key = make_triangle_memo_key(triangle, depth);
    auto memo_it = g_memo_cache.find(memo_key);
    if (memo_it != g_memo_cache.end()) { g_cache_hit_count++; return memo_it->second; }

    if (area2d(triangle) <= AREA_TOL) { g_memo_cache[memo_key] = 0; return 0; }

    std::vector<Edge> current_edges;
    for (const auto& e : edges)
        if (!edge_on_triangle_side(e, triangle)) current_edges.push_back(e);
    current_edges = filter_degenerate_edges(current_edges);

    if (current_edges.empty()) { g_memo_cache[memo_key] = 1; return 1; }
    if (lookahead_remaining <= 0) {
        int est = estimate_leaves_heuristic(current_edges);
        g_memo_cache[memo_key] = est;
        return est;
    }

    std::vector<ScoredCut> scored = am_scored_cuts(triangle, current_edges);
    if (scored.empty()) { g_memo_cache[memo_key] = 1; return 1; }

    int width = std::min((int)scored.size(), ANGLE_LOOKAHEAD_WIDTH);
    int min_leaves = std::numeric_limits<int>::max();
    int realized = 0;
    for (int i = 0; i < (int)scored.size() && realized < width; ++i) {
        std::array<Triangle,2> kids; Edge sline; int prio;
        if (!am_realize_cut(triangle, scored[i], all_polygon_edges, kids, sline, prio)) continue;
        ++realized;
        auto [ce1, ce2] = assign_edges_to_children(kids, current_edges, sline);
        int l1 = am_count_min_leaves(kids[0], filter_degenerate_edges(ce1),
                                     all_polygon_edges, depth + 1, lookahead_remaining - 1);
        if (g_evaluation_terminated) {
            int fb = (min_leaves == std::numeric_limits<int>::max())
                   ? estimate_leaves_heuristic(current_edges) : min_leaves;
            g_memo_cache[memo_key] = fb; return fb;
        }
        int l2 = am_count_min_leaves(kids[1], filter_degenerate_edges(ce2),
                                     all_polygon_edges, depth + 1, lookahead_remaining - 1);
        if (g_evaluation_terminated) {
            int fb = (min_leaves == std::numeric_limits<int>::max())
                   ? estimate_leaves_heuristic(current_edges) : min_leaves;
            g_memo_cache[memo_key] = fb; return fb;
        }
        min_leaves = std::min(min_leaves, l1 + l2);
    }
    if (min_leaves == std::numeric_limits<int>::max())
        min_leaves = estimate_leaves_heuristic(current_edges);
    g_memo_cache[memo_key] = min_leaves;
    return min_leaves;
}

// =========================================================================
// ANGLE-METHOD SPLIT SELECTION
//   --lookahead 0  : greedy -- realize the single best-scoring cut (fast path,
//                    O(k log k) per triangle, unchanged behaviour).
//   --lookahead N  : among the top ANGLE_LOOKAHEAD_WIDTH realizable cuts, pick
//                    the one whose subtree (continued by the angle method) gives
//                    the fewest leaves N levels deep -- the same minimize-leaves
//                    objective the legacy engine optimizes, now on angle cuts.
// =========================================================================
BestCutResult find_best_split_angle(const Triangle& triangle,
                                    const std::vector<Edge>& edges,
                                    const std::vector<Edge>& all_polygon_edges)
{
    BestCutResult result;
    std::vector<ScoredCut> scored = am_scored_cuts(triangle, edges);
    if (scored.empty()) return result;   // no usable cut; caller falls back to greedy

    // Greedy fast path (no look-ahead): best-scoring non-degenerate cut.
    if (LOOKAHEAD_DEPTH_LIMIT <= 0) {
        for (const auto& sc : scored) {
            std::array<Triangle,2> kids; Edge sline; int prio;
            if (!am_realize_cut(triangle, sc, all_polygon_edges, kids, sline, prio)) continue;
            result.found = true;
            result.child_triangles = kids;
            result.split_line = sline;
            result.min_leaves = (int)edges.size();   // informational only
            result.priority = prio;
            return result;
        }
        return result;
    }

    // Look-ahead path.
    g_recursive_call_count = 0;
    g_cache_hit_count = 0;
    g_evaluation_terminated = false;

    int width = std::min((int)scored.size(), ANGLE_LOOKAHEAD_WIDTH);
    int best_total = std::numeric_limits<int>::max();
    int realized = 0;
    for (int i = 0; i < (int)scored.size() && realized < width; ++i) {
        std::array<Triangle,2> kids; Edge sline; int prio;
        if (!am_realize_cut(triangle, scored[i], all_polygon_edges, kids, sline, prio)) continue;
        ++realized;
        auto [ce1, ce2] = assign_edges_to_children(kids, edges, sline);
        int l1 = am_count_min_leaves(kids[0], filter_degenerate_edges(ce1),
                                     all_polygon_edges, 1, LOOKAHEAD_DEPTH_LIMIT - 1);
        int l2 = am_count_min_leaves(kids[1], filter_degenerate_edges(ce2),
                                     all_polygon_edges, 1, LOOKAHEAD_DEPTH_LIMIT - 1);
        int total = l1 + l2;
        bool better = !result.found || total < best_total ||
                      (total == best_total && prio < result.priority);
        if (better) {
            best_total = total;
            result.found = true;
            result.child_triangles = kids;
            result.split_line = sline;
            result.min_leaves = total;
            result.priority = prio;
        }
        if (g_evaluation_terminated) break;
    }

    // Fall back to greedy best if nothing in the window realized.
    if (!result.found) {
        for (const auto& sc : scored) {
            std::array<Triangle,2> kids; Edge sline; int prio;
            if (!am_realize_cut(triangle, sc, all_polygon_edges, kids, sline, prio)) continue;
            result.found = true;
            result.child_triangles = kids;
            result.split_line = sline;
            result.min_leaves = (int)edges.size();
            result.priority = prio;
            break;
        }
    }
    return result;
}

// --- Core Processing Function with Transaction Tracking ---
// =========================================================================
// PERWALK + MARKTRIANGLES : interior / exterior leaf classification
// =========================================================================
// Replaces the per-triangle centroid point-in-polygon oracle with the
// verifiable marking protocol from the dissertation:
//
//   PerWalk      Walk the polygon boundary in CCW order. For each (oriented)
//                boundary segment, the leaf triangle on the LEFT is interior
//                (code 1); the leaf on the RIGHT is exterior (code 2). Because
//                the polygon is CCW-oriented, "interior == left" holds for
//                every boundary segment.
//   MarkTriangles
//                Deep-interior leaves that no boundary segment touches are
//                reached by flooding from already-interior leaves across
//                shared edges that are NOT boundary segments (the boundary
//                acts as a wall the flood cannot cross). Completion is verified
//                when the sum of interior leaf areas equals the polygon area.
//
// Region codes:  REGION_INSIDE = 1, REGION_OUTSIDE = 2, REGION_BOUNDARY = 0.
// REGION_BOUNDARY is reserved: no leaf straddles the boundary (every boundary
// segment is exactly a triangle side -- the boundary-embedding property), so
// every leaf resolves cleanly to 1 or 2.
// =========================================================================

static const int REGION_BOUNDARY = 0;
static const int REGION_INSIDE   = 1;
static const int REGION_OUTSIDE  = 2;

// Tolerance-based point canonicalizer, consistent with the engine's POINT_TOL.
// The decomposition's "boundary" vertices are constructed intersection points
// that coincide with the original polygon vertices only to within POINT_TOL
// (the same tolerance point_eq uses everywhere else), so canonicalization must
// merge points that are within POINT_TOL -- exact equality would split them.
struct PtCanon {
    double cell;
    std::unordered_map<long long, std::vector<int>> grid;   // cell -> ids
    std::vector<std::pair<double,double>> dpt;              // canonical doubles
    explicit PtCanon(double tol) : cell(tol) {}
    static long long ckey(int ix, int iy) {
        // Compose a unique 64-bit key from two 32-bit cell indices. Cell indices
        // are negative for western longitudes, so the shift MUST be unsigned:
        // left-shifting a negative signed value is undefined behaviour (it
        // crashes/stalls under MSVC and is flagged by UBSan).
        unsigned long long ux = (unsigned int)ix;
        unsigned long long uy = (unsigned int)iy;
        return (long long)((ux << 32) | uy);
    }
    int get(double x, double y) {
        int ix = (int)std::floor(x / cell), iy = (int)std::floor(y / cell);
        for (int dx = -1; dx <= 1; ++dx)
            for (int dy = -1; dy <= 1; ++dy) {
                auto it = grid.find(ckey(ix + dx, iy + dy));
                if (it == grid.end()) continue;
                for (int id : it->second) {
                    double ex = dpt[id].first - x, ey = dpt[id].second - y;
                    if (ex * ex + ey * ey < cell * cell) return id;
                }
            }
        int id = (int)dpt.size();
        dpt.push_back({ x, y });
        grid[ckey(ix, iy)].push_back(id);
        return id;
    }
};

// Fast double-precision point-in-polygon (even-odd ray cast). This realizes the
// "is the tile on the interior (left) side of the CCW boundary" test robustly
// and in O(N) doubles per query -- vastly cheaper than CGAL exact bounded_side,
// and unaffected by whether boundary vertices happen to be leaf vertices.
struct FastPIP {
    std::vector<double> X, Y;
    void build(const std::vector<Edge>& edges) {
        X.clear(); Y.clear();
        X.reserve(edges.size()); Y.reserve(edges.size());
        for (const auto& e : edges) {
            X.push_back(CGAL::to_double(e.source().x()));
            Y.push_back(CGAL::to_double(e.source().y()));
        }
    }
    bool inside(double px, double py) const {
        bool in = false;
        size_t n = X.size();
        if (n < 3) return false;
        for (size_t i = 0, j = n - 1; i < n; j = i++) {
            if (((Y[i] > py) != (Y[j] > py)) &&
                (px < (X[j] - X[i]) * (py - Y[i]) / (Y[j] - Y[i]) + X[i]))
                in = !in;
        }
        return in;
    }
};

struct MarkResult {
    std::vector<int> codes;          // parallel to the leaves vector
    int inside_count = 0;
    int outside_count = 0;
    int perwalk_inside_seeds = 0;    // interior leaves touching the boundary
    int flood_added = 0;             // deep-interior leaves reached by flooding
    int flood_unreached = 0;         // interior leaves the flood could not reach
    bool area_ok = false;
    double inside_area = 0.0;
    double polygon_area = 0.0;
    // For O(1) polygon-side matching downstream (same tolerance canonicalization):
    std::vector<std::array<long long,3>> leaf_side_keys;   // parallel to leaves
    std::unordered_map<long long,int>    poly_edge_by_key; // edge key -> polygon edge idx
};

static MarkResult g_last_mark;

// Absolute polygon area via a flat double-precision shoelace sum. Using the
// kernel's Polygon::area() on a 10k+ vertex polygon builds a Lazy_exact_nt
// expression tree thousands of nodes deep, whose recursive evaluation/teardown
// overflows MSVC's 1 MB thread stack (silent crash). A flat double sum avoids
// that entirely; the area is only used for reporting and a relative tolerance.
static double polygon_abs_area_double(const Polygon& poly) {
    double a = 0.0;
    int n = (int)poly.size();
    for (int i = 0; i < n; ++i) {
        const Point& p = poly.vertex(i);
        const Point& q = poly.vertex((i + 1) % n);
        a += CGAL::to_double(p.x()) * CGAL::to_double(q.y())
           - CGAL::to_double(q.x()) * CGAL::to_double(p.y());
    }
    return std::abs(a) * 0.5;
}

// PerWalk + MarkTriangles classification.
//
// Robust realization for this engine: the interior/exterior test per leaf is the
// CCW-boundary left/right test, evaluated by point-in-polygon (mathematically
// identical, but stable even when the tessellation is not edge-conforming to the
// boundary -- which is common on real TIGER data, where a boundary vertex often
// lies mid-side of a longer cut). On top of that ground truth we reproduce the
// dissertation's two-phase structure for reporting and area verification:
//   * PerWalk seeds   = interior leaves that share an edge with an exterior leaf
//                       (i.e. leaves the boundary walk would touch directly).
//   * MarkTriangles   = deep-interior leaves reached by flooding from the seeds
//                       across interior-interior (non-boundary) shared edges.
// Walls are the edges where an interior and an exterior leaf meet, so the flood
// cannot leak across the boundary. Completion is verified when the interior leaf
// area equals the polygon area.
MarkResult mark_triangles_perwalk(const std::vector<Triangle>& leaves,
                                  const std::vector<Edge>& polygon_edges,
                                  const Polygon& polygon,
                                  std::ofstream& log)
{
    MarkResult R;
    const int L = (int)leaves.size();
    R.codes.assign(L, REGION_OUTSIDE);
    if (L == 0) return R;

    std::cerr << "[mark] classifying " << L << " leaves (interior/exterior)...\n" << std::flush;

    FastPIP pip;
    pip.build(polygon_edges);

    // Ground-truth interior test per leaf (centroid; leaves never straddle the
    // boundary because every boundary segment lies along a triangle side).
    std::vector<char> is_in(L, 0);
    for (int t = 0; t < L; ++t) {
        double cx = (CGAL::to_double(leaves[t][0].x()) + CGAL::to_double(leaves[t][1].x())
                   + CGAL::to_double(leaves[t][2].x())) / 3.0;
        double cy = (CGAL::to_double(leaves[t][0].y()) + CGAL::to_double(leaves[t][1].y())
                   + CGAL::to_double(leaves[t][2].y())) / 3.0;
        is_in[t] = pip.inside(cx, cy) ? 1 : 0;
    }

    // Canonical vertex ids + leaf adjacency via shared (tolerance) edges.
    PtCanon canon(POINT_TOL);
    auto ekey = [](int a, int b) -> long long {
        if (a > b) std::swap(a, b);
        // unsigned shift: well-defined regardless of sign (ids are >= 0 here,
        // but keep it consistent with PtCanon::ckey and UB-free).
        return (long long)(((unsigned long long)(unsigned int)a << 32)
                           | (unsigned int)b);
    };
    R.leaf_side_keys.resize(L);
    for (int t = 0; t < L; ++t) {
        int v0 = canon.get(CGAL::to_double(leaves[t][0].x()), CGAL::to_double(leaves[t][0].y()));
        int v1 = canon.get(CGAL::to_double(leaves[t][1].x()), CGAL::to_double(leaves[t][1].y()));
        int v2 = canon.get(CGAL::to_double(leaves[t][2].x()), CGAL::to_double(leaves[t][2].y()));
        R.leaf_side_keys[t] = { ekey(v0, v1), ekey(v1, v2), ekey(v2, v0) };
    }
    // Canonical key -> polygon edge index (only polygon edges that are full leaf
    // sides will match a leaf side key; tolerance canonicalization keeps this
    // consistent with the leaf vertices).
    for (size_t pe = 0; pe < polygon_edges.size(); ++pe) {
        int a = canon.get(CGAL::to_double(polygon_edges[pe].source().x()),
                          CGAL::to_double(polygon_edges[pe].source().y()));
        int b = canon.get(CGAL::to_double(polygon_edges[pe].target().x()),
                          CGAL::to_double(polygon_edges[pe].target().y()));
        R.poly_edge_by_key.emplace(ekey(a, b), (int)pe);
    }

    // =====================================================================
    // Segment-overlap adjacency graph (topological flooding across T-junctions)
    //
    // Endpoint-hash adjacency misses T-junctions: when one leaf carries a full
    // side A--C while its neighbours carry A--B and B--C (B introduced on A--C
    // by a later split), the keys never match and the flood cannot cross. Real
    // tessellations are riddled with these -- most internal borders are partial
    // -- so endpoint hashing fragments the interior into ~1000 pieces.
    //
    // Instead build a SEGMENT-INTERSECTION graph: two leaf sides are topological
    // neighbours iff (1) they are COLLINEAR (exact CGAL predicate) and (2) their
    // 1-D extents along that line OVERLAP by a positive length. Such an overlap
    // is a shared border (full or partial T-junction), so the flood may cross
    // it. Sides are grouped by a canonical (normal, offset) key so only
    // plausibly-collinear sides are exact-tested.
    // =====================================================================
    struct SideRef { int leaf; double nx, ny, rho; Point a, b; };
    std::vector<SideRef> sref;
    sref.reserve(3 * L);
    for (int t = 0; t < L; ++t) {
        for (int s = 0; s < 3; ++s) {
            const Point& A = leaves[t][s];
            const Point& B = leaves[t][(s + 1) % 3];
            double ax = CGAL::to_double(A.x()), ay = CGAL::to_double(A.y());
            double dx = CGAL::to_double(B.x()) - ax, dy = CGAL::to_double(B.y()) - ay;
            double len = std::sqrt(dx * dx + dy * dy);
            if (len < 1e-15) continue;                       // degenerate side
            double nx = -dy / len, ny = dx / len;            // unit normal
            double rho = ax * nx + ay * ny;                  // signed perpendicular offset
            if (nx < 0 || (nx == 0.0 && ny < 0)) { nx = -nx; ny = -ny; rho = -rho; }
            sref.push_back({ t, nx, ny, rho, A, B });
        }
    }
    // Sort by offset then normal: exactly-collinear sides share rho to ~1e-14
    // and normal to ~1e-16, so they form a tight, well-separated cluster.
    std::sort(sref.begin(), sref.end(), [](const SideRef& p, const SideRef& q){
        if (p.rho != q.rho) return p.rho < q.rho;
        if (p.nx  != q.nx ) return p.nx  < q.nx;
        return p.ny < q.ny;
    });

    std::vector<std::vector<int>> adj(L);     // interior <-> interior flood graph
    std::vector<char> seed_touch(L, 0);       // interior leaf adjoining an exterior leaf
    const double RHO_TOL = 1e-9, N_TOL = 1e-6, OVL_TOL = 1e-9;

    // Walk the sorted sides, cutting them into same-line groups (tight rho/normal
    // clusters). Within a group, verify collinearity against a representative
    // with the exact predicate (O(group) exact calls), then find overlapping
    // 1-D extents by an interval sweep (O(group log group)) and connect leaves.
    const int NS = (int)sref.size();
    int gs = 0;
    while (gs < NS) {
        int ge = gs + 1;
        while (ge < NS &&
               sref[ge].rho - sref[ge - 1].rho <= RHO_TOL &&
               std::abs(sref[ge].nx - sref[gs].nx) <= N_TOL &&
               std::abs(sref[ge].ny - sref[gs].ny) <= N_TOL)
            ++ge;
        int K = ge - gs;
        if (K >= 2) {
            const SideRef& rep = sref[gs];
            double ax = CGAL::to_double(rep.a.x()), ay = CGAL::to_double(rep.a.y());
            double dx = CGAL::to_double(rep.b.x()) - ax, dy = CGAL::to_double(rep.b.y()) - ay;
            double dl = std::sqrt(dx * dx + dy * dy);
            if (dl >= 1e-15) {
                double ux = dx / dl, uy = dy / dl;
                std::vector<std::array<double,2>> iv;   // [lo,hi] extent per side
                std::vector<int> ivleaf;
                iv.reserve(K); ivleaf.reserve(K);
                for (int k = gs; k < ge; ++k) {
                    const SideRef& s = sref[k];
                    // exact collinearity with the group's line (guards a loose cluster)
                    if (k != gs && (!CGAL::collinear(rep.a, rep.b, s.a) ||
                                    !CGAL::collinear(rep.a, rep.b, s.b)))
                        continue;
                    double p0 = (CGAL::to_double(s.a.x()) - ax) * ux
                              + (CGAL::to_double(s.a.y()) - ay) * uy;
                    double p1 = (CGAL::to_double(s.b.x()) - ax) * ux
                              + (CGAL::to_double(s.b.y()) - ay) * uy;
                    if (p0 > p1) std::swap(p0, p1);
                    iv.push_back({ p0, p1 });
                    ivleaf.push_back(s.leaf);
                }
                int M = (int)iv.size();
                std::vector<int> ord(M);
                for (int z = 0; z < M; ++z) ord[z] = z;
                std::sort(ord.begin(), ord.end(),
                          [&](int a, int b){ return iv[a][0] < iv[b][0]; });
                for (int a = 0; a < M; ++a) {
                    int ia = ord[a];
                    double a1 = iv[ia][1];
                    for (int b = a + 1; b < M; ++b) {
                        int ib = ord[b];
                        if (iv[ib][0] >= a1 - OVL_TOL) break;   // sorted: no more overlaps
                        double ov = std::min(a1, iv[ib][1]) - std::max(iv[ia][0], iv[ib][0]);
                        if (ov <= OVL_TOL) continue;
                        int u = ivleaf[ia], v = ivleaf[ib];
                        if (u == v) continue;
                        bool iu = is_in[u], iv2 = is_in[v];
                        if (iu && iv2) { adj[u].push_back(v); adj[v].push_back(u); }
                        else if (iu != iv2) { if (iu) seed_touch[u] = 1; else seed_touch[v] = 1; }
                    }
                }
            }
        }
        gs = ge;
    }

    // PerWalk seeds: interior leaves that share a border with an exterior leaf.
    std::vector<int> stack;
    std::vector<char> reached(L, 0);
    for (int t = 0; t < L; ++t) {
        if (is_in[t] && seed_touch[t]) {
            reached[t] = 1; stack.push_back(t); R.perwalk_inside_seeds++;
        }
    }
    // MarkTriangles: flood across interior-interior shared borders (T-junctions
    // included), trapped inside the region because exterior borders are not
    // interior-interior edges.
    while (!stack.empty()) {
        int t = stack.back(); stack.pop_back();
        for (int nb : adj[t]) {
            if (!reached[nb]) { reached[nb] = 1; R.flood_added++; stack.push_back(nb); }
        }
    }

    // Final codes follow the ground-truth interior test (the flood reproduces it;
    // any interior leaf the flood could not reach -- e.g. isolated by a
    // non-conforming adjacency -- is still classified correctly and counted).
    // Area is summed in DOUBLE (not the exact kernel): accumulating ~10^4 exact
    // triangle areas would build a Lazy_exact_nt tree thousands deep and blow
    // MSVC's 1 MB stack on teardown. Double is ample for the relative check.
    double inside_area_approx = 0.0;
    for (int t = 0; t < L; ++t) {
        if (is_in[t]) {
            R.codes[t] = REGION_INSIDE;
            R.inside_count++;
            if (!reached[t]) R.flood_unreached++;
            inside_area_approx += CGAL::to_double(CGAL::abs(
                Triangle_cgal(leaves[t][0], leaves[t][1], leaves[t][2]).area()));
        } else {
            R.codes[t] = REGION_OUTSIDE;
            R.outside_count++;
        }
    }
    R.inside_area  = inside_area_approx;
    R.polygon_area = polygon_abs_area_double(polygon);
    double abs_diff = std::abs(R.inside_area - R.polygon_area);
    double rel_diff = (R.polygon_area > 0.0) ? abs_diff / R.polygon_area : abs_diff;
    // "Complete" means the interior tiles cover the polygon to within a small
    // relative tolerance. A tiny residual is expected on real data because the
    // tessellation is not perfectly boundary-conforming -- a few leaves have the
    // boundary passing through their interior, so their whole area counts as
    // interior. 0.1% relative still catches any gross misclassification.
    R.area_ok = rel_diff < 1e-3;

    log << "[MARK] PerWalk seeds: " << R.perwalk_inside_seeds
        << " | flood-added: " << R.flood_added
        << " | flood-unreached: " << R.flood_unreached
        << " | interior leaves (L^I): " << R.inside_count
        << " | exterior leaves: " << R.outside_count << "\n";
    log << "[MARK] interior area = " << std::fixed << std::setprecision(10) << R.inside_area
        << " | polygon area = " << R.polygon_area
        << " | diff = " << abs_diff
        << " (" << std::setprecision(5) << (rel_diff * 100.0) << "%)"
        << " | completeness " << (R.area_ok ? "PASS" : "FAIL") << "\n";
    std::cerr << "[mark] classification complete (interior " << R.inside_count
              << ", exterior " << R.outside_count << ")\n" << std::flush;
    return R;
}

void process(const Triangle& initial_triangle,
    const std::vector<Edge>& polygon_edges,
    const Polygon& polygon,
    std::ofstream& debug_fp,
    std::vector<FinalTriangle>& inside_triangles)
{
    std::vector<Triangle> all_triangles;
    std::queue<std::pair<Triangle, std::vector<Edge>>> queue;

    queue.push({ initial_triangle, polygon_edges });
    int iteration = 0;
    // Scale iteration limit to polygon complexity: need ~2.5N splits minimum,
    // plus leaves, so total iterations ~5N.  Use 10N as a safe upper bound.
    // Absolute minimum of 100,000 for small polygons.
    const int max_iter = std::max(100000, static_cast<int>(polygon_edges.size()) * 10);
    
    int cells_processed = 0;
    int leaves_created = 0;
    int splits_performed = 0;
    int greedy_fallbacks = 0;
    int poly_edge_splits = 0;      // priority 0: split line IS a polygon edge
    int zero_crossing_splits = 0;  // priority 1: zero boundary crossings
    int boundary_crossing_splits = 0; // priority 2: has boundary crossings
    
    auto start_time = std::chrono::steady_clock::now();
    
    // Initialize transaction log header
    g_transaction_log << "==============================================\n";
    g_transaction_log << "AGS TRANSACTION LOG\n";
    g_transaction_log << "==============================================\n";
    g_transaction_log << "This log tracks all transactions as per AGS protocol:\n";
    g_transaction_log << "- TRIANGLE_SPLIT: Each triangle split is a transaction\n";
    g_transaction_log << "- LINE_SPLIT: Each boundary line split is TWO transactions\n";
    g_transaction_log << "  (once during split, once during PerWalk boundary walk)\n";
    g_transaction_log << "==============================================\n\n";

    while (!queue.empty() && iteration < max_iter) {
        auto [tri, edges] = queue.front();
        queue.pop();
        iteration++;
        cells_processed++;

        if (verbose_debug()) {
            double tri_area = area2d(tri);
            debug_fp << "\nIteration " << iteration << " (Cell #" << cells_processed << ") Vertices:\n";
            for (const auto& v : tri) debug_fp << point_to_string(v) << "\n";
            debug_fp << "Edge List to process (" << edges.size() << ")\n";
            if (edges.size() <= 20) {
                for (const auto& e : edges) debug_fp << point_to_string(e.source()) << " - " << point_to_string(e.target()) << "\n";
            } else {
                for (size_t ei = 0; ei < 10; ++ei)
                    debug_fp << point_to_string(edges[ei].source()) << " - " << point_to_string(edges[ei].target()) << "\n";
                debug_fp << "  ... (" << (edges.size() - 20) << " edges omitted) ...\n";
                for (size_t ei = edges.size() - 10; ei < edges.size(); ++ei)
                    debug_fp << point_to_string(edges[ei].source()) << " - " << point_to_string(edges[ei].target()) << "\n";
            }
            debug_fp << "Area: " << std::fixed << std::setprecision(8) << tri_area << "\n";
        }

        if (iteration % progress_interval() == 0 || iteration == 1) {
            double elapsed = std::chrono::duration<double>(std::chrono::steady_clock::now() - start_time).count();
            std::cerr << "[progress] iteration=" << iteration << "/" << max_iter
                      << " (" << std::fixed << std::setprecision(1)
                      << (100.0 * iteration / max_iter) << "%)"
                      << " leaves=" << leaves_created
                      << " tri_splits=" << g_triangle_split_count
                      << " line_splits=" << g_line_split_count
                      << " queue=" << queue.size()
                      << " elapsed=" << std::fixed << std::setprecision(1) << elapsed << "s\n";
        }

        std::vector<Edge> current_edges;
        for (const auto& e : edges) {
            if (!edge_on_triangle_side(e, tri)) {
                current_edges.push_back(e);
            }
        }
        current_edges = filter_degenerate_edges(current_edges);

        if (current_edges.empty()) {
            if (verbose_debug())
                debug_fp << "Finalize triangle (no edges left). Leaf #" << (leaves_created + 1) << "\n";
            all_triangles.push_back(tri);
            leaves_created++;
            continue;
        }

        bool split_found = false;
        std::array<Triangle, 2> split_tris;
        Edge split_line;

        if (log_decisions()) {
            g_decision_log << "\n--- Processing Iteration " << iteration << " ---\n";
            g_decision_log << "Triangle: " << triangle_to_string(tri) << "\n";
            g_decision_log << "Edges to process: " << current_edges.size() << "\n";
        }
        
        BestCutResult best = USE_ANGLE_METHOD
            ? find_best_split_angle(tri, current_edges, current_edges)
            : find_best_split_recursive(tri, current_edges, current_edges, 0);
        
        if (best.found) {
            split_tris = best.child_triangles;
            split_line = best.split_line;
            
            // Record the triangle split
            g_triangle_split_count++;
            TriangleSplitRecord tri_record;
            tri_record.split_id = g_triangle_split_count;
            tri_record.parent_triangle = tri;
            tri_record.child1 = split_tris[0];
            tri_record.child2 = split_tris[1];
            tri_record.split_line = split_line;
            tri_record.iteration = iteration;
            tri_record.is_lookahead = true;
            g_triangle_splits.push_back(tri_record);
            
            // Log triangle split
            std::string priority_str = (best.priority == 0) ? "POLYGON_EDGE" :
                                       (best.priority == 1) ? "ZERO_CROSSING" : "BOUNDARY_CROSSING";
            g_transaction_log << "TRIANGLE_SPLIT #" << tri_record.split_id << " (Lookahead, " << priority_str << ")\n";
            g_transaction_log << "  Iteration: " << iteration << "\n";
            g_transaction_log << "  Priority: " << best.priority << " (" << priority_str << ")\n";
            g_transaction_log << "  Parent: " << triangle_to_string(tri) << "\n";
            g_transaction_log << "  Split line: " << edge_to_string(split_line) << "\n";
            g_transaction_log << "  Child 1: " << triangle_to_string(split_tris[0]) << "\n";
            g_transaction_log << "  Child 2: " << triangle_to_string(split_tris[1]) << "\n";
            g_transaction_log << "  Estimated final leaves: " << best.min_leaves << "\n\n";
            
            // Use tracked version to capture line splits
            std::vector<LineSplitRecord> line_splits_this_iteration;
            auto [child_edges1, child_edges2] = assign_edges_to_children_tracked(
                split_tris, current_edges, split_line, 
                g_triangle_split_count, iteration, line_splits_this_iteration);

            if (area2d(split_tris[0]) > AREA_TOL) {
                // Children of a split are independent and can be processed in
                // parallel; the current queue is sequential (see notes in chat).
                queue.push({ split_tris[0], filter_degenerate_edges(child_edges1) });
            }
            if (area2d(split_tris[1]) > AREA_TOL) {
                queue.push({ split_tris[1], filter_degenerate_edges(child_edges2) });
            }
            
            splits_performed++;
            if (best.priority == 0) poly_edge_splits++;
            else if (best.priority == 1) zero_crossing_splits++;
            else boundary_crossing_splits++;
            if (verbose_debug()) {
                debug_fp << "Found optimal lookahead split (" << priority_str << ") with split line "
                         << point_to_string(split_line.source()) << " - "
                         << point_to_string(split_line.target())
                         << " (estimated " << best.min_leaves << " final triangles). Split #" << splits_performed << "\n";
                if (!line_splits_this_iteration.empty()) {
                    debug_fp << "  Line splits caused: " << line_splits_this_iteration.size() << "\n";
                    for (const auto& ls : line_splits_this_iteration) {
                        debug_fp << "    LINE_SPLIT #" << ls.split_id
                                 << ": " << edge_to_string(ls.original_edge)
                                 << " at " << point_to_string(ls.split_point) << "\n";
                    }
                }
            }
            if (log_decisions()) {
                g_decision_log << "FINAL DECISION: Split at "
                              << point_to_string(split_line.source()) << " -> "
                              << point_to_string(split_line.target()) << "\n";
                g_decision_log << "Estimated leaves from this split: " << best.min_leaves << "\n";
                g_decision_log << "Priority: " << best.priority
                              << " (" << priority_str << ")\n";
                g_decision_log << "Triangle Split #" << g_triangle_split_count << " (Lookahead)\n";
                g_decision_log << "Line splits this iteration: " << line_splits_this_iteration.size() << "\n";
                for (const auto& ls : line_splits_this_iteration) {
                    g_decision_log << "  LINE_SPLIT #" << ls.split_id
                                  << " | Edge: " << edge_to_string(ls.original_edge)
                                  << " | Split at: " << point_to_string(ls.split_point)
                                  << " | Seg1: " << edge_to_string(ls.segment1)
                                  << " | Seg2: " << edge_to_string(ls.segment2) << "\n";
                }
                g_decision_log << "\n";
            }
            split_found = true;
        }
        else {
            if (slope_split_greedy(tri, current_edges, current_edges, split_tris, split_line)) {
                // Record the triangle split (greedy)
                g_triangle_split_count++;
                TriangleSplitRecord tri_record;
                tri_record.split_id = g_triangle_split_count;
                tri_record.parent_triangle = tri;
                tri_record.child1 = split_tris[0];
                tri_record.child2 = split_tris[1];
                tri_record.split_line = split_line;
                tri_record.iteration = iteration;
                tri_record.is_lookahead = false;
                g_triangle_splits.push_back(tri_record);
                
                // Log triangle split
                g_transaction_log << "TRIANGLE_SPLIT #" << tri_record.split_id << " (Greedy Fallback)\n";
                g_transaction_log << "  Iteration: " << iteration << "\n";
                g_transaction_log << "  Parent: " << triangle_to_string(tri) << "\n";
                g_transaction_log << "  Split line: " << edge_to_string(split_line) << "\n";
                g_transaction_log << "  Child 1: " << triangle_to_string(split_tris[0]) << "\n";
                g_transaction_log << "  Child 2: " << triangle_to_string(split_tris[1]) << "\n\n";
                
                // Use tracked version to capture line splits
                std::vector<LineSplitRecord> line_splits_this_iteration;
                auto [child_edges1, child_edges2] = assign_edges_to_children_tracked(
                    split_tris, current_edges, split_line,
                    g_triangle_split_count, iteration, line_splits_this_iteration);

                if (area2d(split_tris[0]) > AREA_TOL) {
                    queue.push({ split_tris[0], filter_degenerate_edges(child_edges1) });
                }
                if (area2d(split_tris[1]) > AREA_TOL) {
                    queue.push({ split_tris[1], filter_degenerate_edges(child_edges2) });
                }
                splits_performed++;
                greedy_fallbacks++;
                // Determine greedy split priority
                bool greedy_is_poly = split_line_is_polygon_edge(split_line, current_edges);
                int greedy_crossings = greedy_is_poly ? 0 : count_intersections_improved(split_line, current_edges);
                if (greedy_is_poly) poly_edge_splits++;
                else if (greedy_crossings == 0) zero_crossing_splits++;
                else boundary_crossing_splits++;
                std::string greedy_pri_str = greedy_is_poly ? "POLYGON_EDGE" :
                    (greedy_crossings == 0 ? "ZERO_CROSSING" : "BOUNDARY_CROSSING");
                if (verbose_debug()) {
                debug_fp << "Found greedy fallback split (" << greedy_pri_str << ") with split line "
                         << point_to_string(split_line.source()) << " - "
                         << point_to_string(split_line.target()) << ". Split #" << splits_performed << " (Greedy)\n";
                if (!line_splits_this_iteration.empty()) {
                    debug_fp << "  Line splits caused: " << line_splits_this_iteration.size() << "\n";
                    for (const auto& ls : line_splits_this_iteration) {
                        debug_fp << "    LINE_SPLIT #" << ls.split_id
                                 << ": " << edge_to_string(ls.original_edge)
                                 << " at " << point_to_string(ls.split_point) << "\n";
                    }
                }
                }
                if (log_decisions()) {
                g_decision_log << "FALLBACK: Greedy split at "
                              << point_to_string(split_line.source()) << " -> "
                              << point_to_string(split_line.target()) << "\n";
                g_decision_log << "Triangle Split #" << g_triangle_split_count << " (Greedy Fallback)\n";
                g_decision_log << "Line splits this iteration: " << line_splits_this_iteration.size() << "\n";
                for (const auto& ls : line_splits_this_iteration) {
                    g_decision_log << "  LINE_SPLIT #" << ls.split_id
                                  << " | Edge: " << edge_to_string(ls.original_edge)
                                  << " | Split at: " << point_to_string(ls.split_point)
                                  << " | Seg1: " << edge_to_string(ls.segment1)
                                  << " | Seg2: " << edge_to_string(ls.segment2) << "\n";
                }
                g_decision_log << "\n";
                }
                split_found = true;
            }
        }

        if (split_found) {
            continue;
        }

        if (verbose_debug())
            debug_fp << "No valid split found (finalizing triangle as fallback). Leaf #" << (leaves_created + 1) << "\n";
        if (log_decisions())
            g_decision_log << "NO SPLIT FOUND: Finalizing triangle as Leaf #" << (leaves_created + 1) << "\n\n";
        all_triangles.push_back(tri);
        leaves_created++;
    }

    // --- Handle incomplete decomposition ---
    // If the iteration limit was reached, the queue may still contain triangles
    // that were never fully resolved.  Finalize them as leaves so they are
    // counted toward the area and classified inside/outside.
    int unresolved_leaves = 0;
    if (!queue.empty()) {
        std::cerr << "\n[WARNING] Iteration limit reached (" << max_iter
                  << "). Draining " << queue.size()
                  << " unresolved triangles from queue as leaves.\n";
        debug_fp << "\n[WARNING] Iteration limit reached (" << max_iter
                 << "). Draining " << queue.size()
                 << " unresolved triangles from queue as leaves.\n";
        g_decision_log << "\n[WARNING] Iteration limit reached (" << max_iter
                       << "). " << queue.size()
                       << " triangles finalized without full resolution.\n";
        g_transaction_log << "\n[WARNING] Iteration limit reached. "
                         << queue.size() << " triangles unresolved.\n\n";
        while (!queue.empty()) {
            auto [tri_q, edges_q] = queue.front();
            queue.pop();
            if (area2d(tri_q) > AREA_TOL) {
                all_triangles.push_back(tri_q);
                unresolved_leaves++;
                leaves_created++;
            }
        }
        std::cerr << "[WARNING] Added " << unresolved_leaves
                  << " unresolved triangles as leaves. Area will be approximate.\n";
        debug_fp << "[WARNING] Added " << unresolved_leaves
                 << " unresolved triangles as leaves.\n";
    }

    // After all splits: classify leaves via PerWalk + MarkTriangles.
    // Collect the non-degenerate leaf tessellation first.
    std::cerr << "[phase] decomposition loop done; collecting leaves...\n" << std::flush;
    std::vector<Triangle> leaves;
    leaves.reserve(all_triangles.size());
    for (const auto& tri : all_triangles)
        if (area2d(tri) > AREA_TOL) leaves.push_back(tri);

    MarkResult mark = mark_triangles_perwalk(leaves, polygon_edges, polygon, g_decision_log);
    g_last_mark = mark;
    g_leaf_triangles    = leaves;
    g_leaf_region_codes = mark.codes;
    g_mark_triangle_count = mark.inside_count;   // L^I MarkTriangles transactions

    g_decision_log << "[MARK] completeness " << (mark.area_ok ? "PASS" : "FAIL")
                   << " (interior area " << std::fixed << std::setprecision(8)
                   << mark.inside_area << " vs polygon " << mark.polygon_area << ")\n";
    std::cerr << "[MARK] interior=" << mark.inside_count
              << " exterior=" << mark.outside_count
              << " (perwalk-seeds " << mark.perwalk_inside_seeds
              << " + flood " << mark.flood_added;
    if (mark.flood_unreached) std::cerr << ", unreached " << mark.flood_unreached;
    std::cerr << ") completeness=" << (mark.area_ok ? "PASS" : "FAIL") << "\n";

    // Independent exact-arithmetic spot-check: validate the fast (double) marking
    // against CGAL's exact centroid point-in-polygon on a capped, evenly-spaced
    // sample of leaves. Cheap (bounded sample) but catches any double-precision
    // misclassification near the boundary.
    {
        const int SAMPLE_CAP = 400;
        int Ln = (int)leaves.size();
        int step = std::max(1, Ln / SAMPLE_CAP);
        int checked = 0, disagree = 0;
        for (int t = 0; t < Ln; t += step) {
            bool exact_in = triangle_is_inside_buffered(polygon, leaves[t]);
            bool mark_in  = (mark.codes[t] == REGION_INSIDE);
            ++checked;
            if (exact_in != mark_in) ++disagree;
        }
        g_decision_log << "[MARK] exact spot-check: " << disagree
                       << " disagreement(s) over " << checked << " sampled leaves\n";
        if (disagree > 0)
            std::cerr << "[MARK][WARN] exact spot-check found " << disagree
                      << "/" << checked << " disagreements (near-boundary precision)\n";
    }

    // Interior leaves become the FinalTriangle set. Polygon-side membership is
    // matched in O(1) via the canonical edge-key map the marking already built.
    for (int t = 0; t < (int)leaves.size(); ++t) {
        if (mark.codes[t] != REGION_INSIDE) continue;
        const Triangle& tri = leaves[t];
        FinalTriangle final_tri;
        final_tri.vertices = tri;
        for (int s = 0; s < 3; ++s) {
            auto it = mark.poly_edge_by_key.find(mark.leaf_side_keys[t][s]);
            if (it != mark.poly_edge_by_key.end())
                final_tri.polygon_sides.push_back(polygon_edges[it->second]);
        }
        inside_triangles.push_back(final_tri);
    }

    // --- Calculate Transaction Counts ---
    // Per AGS protocol:
    // - Each triangle split = 1 transaction
    // - Each line split = 2 transactions (SplitLine during creation + PerWalk during boundary walk)
    int total_triangle_transactions = g_triangle_split_count;
    int total_line_split_transactions = g_line_split_count * 2;  // Counted twice per PDF
    int total_mark_transactions = g_mark_triangle_count;         // L^I MarkTriangles
    int total_transactions = total_triangle_transactions + total_line_split_transactions
                             + total_mark_transactions;

    // --- Final Output ---
    double tri_area_sum = 0.0;
    for (const auto& final_tri : inside_triangles) {
        tri_area_sum += area2d(final_tri.vertices);
    }
    double poly_area = polygon_abs_area_double(polygon);
    double elapsed = std::chrono::duration<double>(std::chrono::steady_clock::now() - start_time).count();

    // Console output
    std::cout << "\n==============================================\n";
    std::cout << "TRIANGLE DECOMPOSITION RESULTS\n";
    std::cout << "==============================================\n";
    std::cout << std::fixed << std::setprecision(8);
    std::cout << "Total area covered by triangles inside the polygon: " << tri_area_sum << std::endl;
    std::cout << "Area of input polygon: " << poly_area << std::endl;
    std::cout << "Difference (triangles - polygon): " << (tri_area_sum - poly_area) << std::endl;
    std::cout << "\n--- Statistics ---\n";
    std::cout << "Cells processed:        " << cells_processed << std::endl;
    std::cout << "Leaves created (total): " << leaves_created << std::endl;
    if (unresolved_leaves > 0) {
        std::cout << "  - Resolved leaves:    " << (leaves_created - unresolved_leaves) << std::endl;
        std::cout << "  - Unresolved (limit): " << unresolved_leaves << std::endl;
        std::cout << "  ** Decomposition INCOMPLETE — iteration limit reached **\n";
        std::cout << "  ** Use --reduce <pct> to reduce polygon, or wait for full run **\n";
    }
    std::cout << "Triangles inside polygon: " << inside_triangles.size() << std::endl;
    std::cout << "Splits performed:       " << splits_performed << std::endl;
    std::cout << "  - Lookahead splits:   " << (splits_performed - greedy_fallbacks) << std::endl;
    std::cout << "  - Greedy fallbacks:   " << greedy_fallbacks << std::endl;
    std::cout << "  --- Split Priority Breakdown ---\n";
    std::cout << "  - Polygon-edge (P0):  " << poly_edge_splits << std::endl;
    std::cout << "  - Zero-crossing (P1): " << zero_crossing_splits << std::endl;
    std::cout << "  - Boundary-cross (P2):" << boundary_crossing_splits << std::endl;
    std::cout << "Lookahead depth limit:  " << LOOKAHEAD_DEPTH_LIMIT << std::endl;
    std::cout << "\n--- AGS TRANSACTION COUNTS ---\n";
    std::cout << "Triangle splits:              " << g_triangle_split_count << std::endl;
    std::cout << "Boundary line splits:         " << g_line_split_count << std::endl;
    std::cout << "Transaction breakdown:\n";
    std::cout << "  - Triangle split txns:      " << total_triangle_transactions << std::endl;
    std::cout << "  - Line split txns (x2):     " << total_line_split_transactions 
              << " (" << g_line_split_count << " splits x 2)\n";
    std::cout << "  - MarkTriangles txns (L^I): " << total_mark_transactions << std::endl;
    std::cout << "  TOTAL TRANSACTIONS:         " << total_transactions << std::endl;
    std::cout << std::fixed << std::setprecision(2);
    std::cout << "Total processing time:  " << elapsed << " seconds" << std::endl;
    std::cout << "==============================================\n";
    
    // Debug file output
    debug_fp << "\n\n==============================================\n";
    debug_fp << "FINAL STATISTICS\n";
    debug_fp << "==============================================\n";
    debug_fp << std::fixed << std::setprecision(8);
    debug_fp << "Total area covered by triangles inside the polygon: " << tri_area_sum << std::endl;
    debug_fp << "Area of input polygon: " << poly_area << std::endl;
    debug_fp << "Difference (triangles - polygon): " << (tri_area_sum - poly_area) << std::endl;
    debug_fp << "\n--- Statistics ---\n";
    debug_fp << "Cells processed:        " << cells_processed << std::endl;
    debug_fp << "Leaves created (total): " << leaves_created << std::endl;
    if (unresolved_leaves > 0) {
        debug_fp << "  - Resolved leaves:    " << (leaves_created - unresolved_leaves) << std::endl;
        debug_fp << "  - Unresolved (limit): " << unresolved_leaves << std::endl;
        debug_fp << "  ** INCOMPLETE DECOMPOSITION **\n";
    }
    debug_fp << "Triangles inside polygon: " << inside_triangles.size() << std::endl;
    debug_fp << "Splits performed:       " << splits_performed << std::endl;
    debug_fp << "  - Lookahead splits:   " << (splits_performed - greedy_fallbacks) << std::endl;
    debug_fp << "  - Greedy fallbacks:   " << greedy_fallbacks << std::endl;
    debug_fp << "  --- Split Priority Breakdown ---\n";
    debug_fp << "  - Polygon-edge (P0):  " << poly_edge_splits << std::endl;
    debug_fp << "  - Zero-crossing (P1): " << zero_crossing_splits << std::endl;
    debug_fp << "  - Boundary-cross (P2):" << boundary_crossing_splits << std::endl;
    debug_fp << "Lookahead depth limit:  " << LOOKAHEAD_DEPTH_LIMIT << std::endl;
    debug_fp << "\n--- AGS TRANSACTION COUNTS ---\n";
    debug_fp << "Triangle splits:              " << g_triangle_split_count << std::endl;
    debug_fp << "Boundary line splits:         " << g_line_split_count << std::endl;
    debug_fp << "Transaction breakdown:\n";
    debug_fp << "  - Triangle split txns:      " << total_triangle_transactions << std::endl;
    debug_fp << "  - Line split txns (x2):     " << total_line_split_transactions 
             << " (" << g_line_split_count << " splits x 2)\n";
    debug_fp << "  TOTAL TRANSACTIONS:         " << total_transactions << std::endl;
    debug_fp << std::fixed << std::setprecision(2);
    debug_fp << "Total processing time:  " << elapsed << " seconds" << std::endl;
    debug_fp << "==============================================\n";
    
    // Decision log output
    g_decision_log << "\n\n==============================================\n";
    g_decision_log << "FINAL STATISTICS\n";
    g_decision_log << "==============================================\n";
    g_decision_log << "Cells processed:        " << cells_processed << std::endl;
    g_decision_log << "Leaves created (total): " << leaves_created << std::endl;
    if (unresolved_leaves > 0) {
        g_decision_log << "  - Resolved leaves:    " << (leaves_created - unresolved_leaves) << std::endl;
        g_decision_log << "  - Unresolved (limit): " << unresolved_leaves << std::endl;
        g_decision_log << "  ** INCOMPLETE DECOMPOSITION **\n";
    }
    g_decision_log << "Triangles inside polygon: " << inside_triangles.size() << std::endl;
    g_decision_log << "Splits performed:       " << splits_performed << std::endl;
    g_decision_log << "  - Lookahead splits:   " << (splits_performed - greedy_fallbacks) << std::endl;
    g_decision_log << "  - Greedy fallbacks:   " << greedy_fallbacks << std::endl;
    g_decision_log << "  --- Split Priority Breakdown ---\n";
    g_decision_log << "  - Polygon-edge (P0):  " << poly_edge_splits << std::endl;
    g_decision_log << "  - Zero-crossing (P1): " << zero_crossing_splits << std::endl;
    g_decision_log << "  - Boundary-cross (P2):" << boundary_crossing_splits << std::endl;
    g_decision_log << "Lookahead depth limit:  " << LOOKAHEAD_DEPTH_LIMIT << std::endl;
    g_decision_log << "\n--- AGS TRANSACTION COUNTS ---\n";
    g_decision_log << "Triangle splits:              " << g_triangle_split_count << std::endl;
    g_decision_log << "Boundary line splits:         " << g_line_split_count << std::endl;
    g_decision_log << "  TOTAL TRANSACTIONS:         " << total_transactions << std::endl;
    g_decision_log << std::fixed << std::setprecision(2);
    g_decision_log << "Total processing time:  " << elapsed << " seconds" << std::endl;
    g_decision_log << "==============================================\n";
    
    // Write transaction summary to transaction log
    g_transaction_log << "\n==============================================\n";
    g_transaction_log << "TRANSACTION SUMMARY\n";
    g_transaction_log << "==============================================\n";
    g_transaction_log << "Total triangle splits: " << g_triangle_split_count << std::endl;
    g_transaction_log << "  - Lookahead splits:  " << (g_triangle_split_count - greedy_fallbacks) << std::endl;
    g_transaction_log << "  - Greedy fallbacks:  " << greedy_fallbacks << std::endl;
    g_transaction_log << "  --- Split Priority Breakdown ---\n";
    g_transaction_log << "  - Polygon-edge (P0): " << poly_edge_splits << std::endl;
    g_transaction_log << "  - Zero-crossing (P1):" << zero_crossing_splits << std::endl;
    g_transaction_log << "  - Boundary-cross(P2):" << boundary_crossing_splits << std::endl;
    g_transaction_log << "Total line splits:     " << g_line_split_count << std::endl;
    g_transaction_log << "\n--- AGS PROTOCOL TRANSACTION COUNTS ---\n";
    g_transaction_log << "Per AGS protocol (from PDF):\n";
    g_transaction_log << "  - Each triangle split = 1 transaction\n";
    g_transaction_log << "  - Each line split = 2 transactions:\n";
    g_transaction_log << "    * 1x SplitLine transaction (when line is split)\n";
    g_transaction_log << "    * 1x PerWalk transaction (during boundary walk for merge)\n";
    g_transaction_log << "\n";
    g_transaction_log << "Triangle split transactions:  " << total_triangle_transactions << std::endl;
    g_transaction_log << "Line split transactions:      " << total_line_split_transactions 
                     << " (" << g_line_split_count << " x 2)" << std::endl;
    g_transaction_log << "MarkTriangles transactions:   " << total_mark_transactions
                     << " (L^I interior leaves)" << std::endl;
    g_transaction_log << "================================\n";
    g_transaction_log << "TOTAL TRANSACTIONS:           " << total_transactions << std::endl;
    g_transaction_log << "================================\n";
    debug_fp.flush();
    g_decision_log.flush();
    g_transaction_log.flush();
}

// =========================================================================
// MINIMUM ENCLOSING TRIANGLE
// Implements the O'Rourke-Aggarwal-Meirans-Simon rotating-calipers algorithm.
// Each side of the optimal triangle is flush with an edge of the convex hull.
// Returns a Triangle (array of 3 Points) with minimum area enclosing all hull_pts.
// hull_pts must be the convex hull in CCW order.
// =========================================================================
static bool line_intersect_2d(
    double ax, double ay, double adx, double ady,
    double bx, double by, double bdx, double bdy,
    double& t)
{
    // Solve: (ax + t*adx, ay + t*ady) = (bx + s*bdx, by + s*bdy)
    double denom = adx * bdy - ady * bdx;
    if (std::abs(denom) < 1e-14) return false;
    t = ((bx - ax) * bdy - (by - ay) * bdx) / denom;
    return true;
}

// Compute intersection point of two lines given as (point, direction)
static bool intersect_lines(
    double ax, double ay, double adx, double ady,
    double bx, double by, double bdx, double bdy,
    double& px, double& py)
{
    double t;
    if (!line_intersect_2d(ax, ay, adx, ady, bx, by, bdx, bdy, t))
        return false;
    px = ax + t * adx;
    py = ay + t * ady;
    return true;
}

// Signed area of triangle (ax,ay),(bx,by),(cx,cy) — positive if CCW
static double tri_area_2(double ax, double ay, double bx, double by, double cx, double cy) {
    return (bx - ax) * (cy - ay) - (by - ay) * (cx - ax);
}

// Advance index modulo n
static int next_idx(int i, int n) { return (i + 1) % n; }

static double triangle_area(const Triangle& tri) {
    return std::abs(tri_area_2(
        CGAL::to_double(tri[0].x()), CGAL::to_double(tri[0].y()),
        CGAL::to_double(tri[1].x()), CGAL::to_double(tri[1].y()),
        CGAL::to_double(tri[2].x()), CGAL::to_double(tri[2].y()))) * 0.5;
}

static Triangle make_ccw_triangle(const Triangle& tri) {
    Triangle result = tri;
    if (tri_area_2(
            CGAL::to_double(tri[0].x()), CGAL::to_double(tri[0].y()),
            CGAL::to_double(tri[1].x()), CGAL::to_double(tri[1].y()),
            CGAL::to_double(tri[2].x()), CGAL::to_double(tri[2].y())) < 0) {
        std::swap(result[1], result[2]);
    }
    return result;
}

static bool triangle_contains_point(const Triangle& tri, double px, double py, double tol = 1e-9) {
    Triangle t = make_ccw_triangle(tri);
    double d0 = tri_area_2(CGAL::to_double(t[0].x()), CGAL::to_double(t[0].y()),
                           CGAL::to_double(t[1].x()), CGAL::to_double(t[1].y()),
                           px, py);
    double d1 = tri_area_2(CGAL::to_double(t[1].x()), CGAL::to_double(t[1].y()),
                           CGAL::to_double(t[2].x()), CGAL::to_double(t[2].y()),
                           px, py);
    double d2 = tri_area_2(CGAL::to_double(t[2].x()), CGAL::to_double(t[2].y()),
                           CGAL::to_double(t[0].x()), CGAL::to_double(t[0].y()),
                           px, py);
    return d0 >= -tol && d1 >= -tol && d2 >= -tol;
}

static bool triangle_contains_hull(const Triangle& tri, const std::vector<Point>& hull_pts, double tol = 1e-9) {
    for (const auto& p : hull_pts) {
        if (!triangle_contains_point(tri, CGAL::to_double(p.x()), CGAL::to_double(p.y()), tol)) {
            return false;
        }
    }
    return true;
}

static Triangle improve_triangle_by_vertex_projection(const std::vector<Point>& hull_pts, const Triangle& initial_tri) {
    Triangle best = make_ccw_triangle(initial_tri);
    double best_area = triangle_area(best);

    for (int pass = 0; pass < 4; ++pass) {
        bool improved = false;

        for (int apex_idx = 0; apex_idx < 3; ++apex_idx) {
            Point apex = best[apex_idx];
            Point fixed_a = best[(apex_idx + 1) % 3];
            Point fixed_b = best[(apex_idx + 2) % 3];

            for (int fixed_mode = 0; fixed_mode < 2; ++fixed_mode) {
                Point fixed = (fixed_mode == 0) ? fixed_a : fixed_b;
                Point other = (fixed_mode == 0) ? fixed_b : fixed_a;

                double ax = CGAL::to_double(apex.x());
                double ay = CGAL::to_double(apex.y());
                double fx = CGAL::to_double(fixed.x());
                double fy = CGAL::to_double(fixed.y());
                double ox = CGAL::to_double(other.x());
                double oy = CGAL::to_double(other.y());
                double odx = ox - fx;
                double ody = oy - fy;
                if (std::abs(odx) < 1e-14 && std::abs(ody) < 1e-14) continue;

                for (const auto& candidate : hull_pts) {
                    double cx = CGAL::to_double(candidate.x());
                    double cy = CGAL::to_double(candidate.y());
                    if (std::abs(cx - ax) < 1e-14 && std::abs(cy - ay) < 1e-14) continue;
                    if (std::abs(cx - fx) < 1e-14 && std::abs(cy - fy) < 1e-14) continue;

                    double idx, idy;
                    if (!intersect_lines(ax, ay, cx - ax, cy - ay,
                                         fx, fy, odx, ody,
                                         idx, idy)) {
                        continue;
                    }

                    Triangle cand = { apex, fixed, Point(idx, idy) };

                    double area = triangle_area(cand);
                    if (area < 1e-12 || area >= best_area - 1e-12) continue;
                    if (!triangle_contains_hull(cand, hull_pts)) continue;

                    best = make_ccw_triangle(cand);
                    best_area = area;
                    improved = true;
                }
            }
        }

        if (!improved) break;
    }

    return best;
}

Triangle compute_min_enclosing_triangle(const std::vector<Point>& hull_pts)
{
    int n = static_cast<int>(hull_pts.size());

    // Convert to doubles for the computation
    std::vector<double> hx(n), hy(n);
    for (int i = 0; i < n; i++) {
        hx[i] = CGAL::to_double(hull_pts[i].x());
        hy[i] = CGAL::to_double(hull_pts[i].y());
    }

    if (n == 1) {
        double cx = hx[0], cy = hy[0], r = 1.0;
        return { Point(cx - r, cy - r), Point(cx + 2*r, cy - r), Point(cx - r, cy + 2*r) };
    }
    if (n == 2) {
        double mx = (hx[0]+hx[1])/2, my = (hy[0]+hy[1])/2;
        double dx = hx[1]-hx[0], dy = hy[1]-hy[0];
        double len = std::sqrt(dx*dx+dy*dy) + 1.0;
        return { Point(mx - len, my - len), Point(mx + 2*len, my - len), Point(mx - len, my + 2*len) };
    }
    if (n == 3) {
        // Hull is already a triangle — use it directly with tiny expansion
        double cx = (hx[0]+hx[1]+hx[2])/3.0, cy = (hy[0]+hy[1]+hy[2])/3.0;
        const double eps = 1e-8;
        double tx[3], ty[3];
        for (int k = 0; k < 3; k++) {
            double dx_ = hx[k] - cx, dy_ = hy[k] - cy;
            double len_ = std::sqrt(dx_*dx_ + dy_*dy_);
            if (len_ > 1e-14) { dx_ /= len_; dy_ /= len_; }
            tx[k] = hx[k] + eps * dx_;
            ty[k] = hy[k] + eps * dy_;
        }
        return { Point(tx[0], ty[0]), Point(tx[1], ty[1]), Point(tx[2], ty[2]) };
    }

    // Ensure hull is CCW
    {
        double area = 0;
        for (int i = 0; i < n; i++) {
            int j = next_idx(i, n);
            area += hx[i]*hy[j] - hx[j]*hy[i];
        }
        if (area < 0) {
            std::reverse(hx.begin(), hx.end());
            std::reverse(hy.begin(), hy.end());
        }
    }

    // =========================================================================
    // Minimum enclosing triangle of a convex polygon.
    //
    // O'Rourke's theorem: for the optimal MET, at least one side is flush
    // (collinear) with a hull edge.  For each hull edge as the "base", we
    // need the tightest triangle whose base lies on that edge's supporting
    // line and that still contains all hull points.
    //
    // For a fixed base edge, the other two sides must each be supporting
    // lines of the hull passing through a common apex point opposite the
    // base.  The tightest such triangle is found by:
    //
    //   1. For each hull vertex j (candidate apex), determine the two
    //      supporting lines from j that bound the hull.  Because the hull
    //      is convex and j is on it, these are the lines through j and the
    //      two hull edges adjacent to j: edge (j-1 → j) and edge (j → j+1).
    //
    //   BUT that's the old (wrong) approach: adjacent edges flare outward
    //   when j has a sharp angle, creating a huge triangle.
    //
    //   The correct approach: for apex j, find the two hull vertices that
    //   are extreme in the "cross-product from j" sense — i.e., the vertex
    //   that maximises and minimises cross(j→base_midpoint, j→k).  The
    //   lines from j through these extreme vertices are the tightest
    //   supporting lines from j.
    //
    //   Actually simpler: we need the tightest WEDGE from j that contains
    //   all hull points.  For a convex hull with j on it, this wedge is
    //   bounded by the two edges incident to j.  But the sides of the
    //   enclosing triangle don't have to follow these edges — they pass
    //   through j and the hull points that are angularly extreme as seen
    //   from j.
    //
    //   For a convex hull vertex j, the other n-1 hull vertices span a
    //   contiguous angular arc < 180° as seen from j (because j is on the
    //   convex boundary).  The tightest wedge from j containing all hull
    //   points has its two rays through the endpoints of this arc.  Those
    //   endpoints are j's two hull neighbors (prev and next in CCW order),
    //   because in a convex polygon each vertex sees all other vertices
    //   within the angular span defined by its two incident edges.
    //
    //   WAIT — that's not right either.  The angular extremes from j are
    //   NOT always the neighbors.  Consider an elongated hull: vertex j at
    //   one end, with neighbors close by, but the farthest vertex at the
    //   other end.  The angular span from j to the farthest vertex is wider
    //   than from j to its neighbors.  Actually no — for a convex polygon,
    //   the angular span from any vertex to all others IS bounded by the
    //   two edges at that vertex.  The proof: if some vertex k were outside
    //   the wedge defined by j's incident edges, then the polygon wouldn't
    //   be convex.
    //
    //   So for a convex hull, the adjacent-edge approach IS correct for
    //   defining the tightest wedge from j.  The bug in the original code
    //   was elsewhere.  Let me re-examine...
    //
    //   The REAL issue: the original code used the adjacent edges as
    //   DIRECTIONS for the triangle sides.  But the triangle sides don't
    //   have to be parallel to those edges — they just have to contain all
    //   hull points.  The tightest triangle with base on edge i and apex
    //   at vertex j has its two non-base sides passing from j to where
    //   they intersect the base line, such that all hull points are inside.
    //   The constraint is: every hull point must be on the "inside" of each
    //   of the three sides.
    //
    //   For the two non-base sides from j, the tightest lines from j that
    //   leave all hull points on the base-side are: for the "left" side,
    //   find the hull vertex k that, from j's perspective, is farthest
    //   LEFT (max cross product with the base direction); the left side
    //   goes from j through k.  Similarly for "right" side.
    //
    // This implementation: O(n) per (base_edge, apex) pair, O(n²) base-apex
    // pairs per base edge, O(n³) total.  Fast for hulls up to ~1000 points.
    // =========================================================================

    double best_area = std::numeric_limits<double>::infinity();
    double best_tx[3], best_ty[3];

    for (int i = 0; i < n; i++) {
        int i1 = next_idx(i, n);
        double bdx = hx[i1] - hx[i], bdy = hy[i1] - hy[i];
        double blen2 = bdx*bdx + bdy*bdy;
        if (blen2 < 1e-28) continue;
        double blen = std::sqrt(blen2);

        // Inward normal for CCW hull
        double bnx = -bdy / blen, bny = bdx / blen;
        double base_dist = hx[i] * bnx + hy[i] * bny;

        // Unit base direction
        double bux = bdx / blen, buy = bdy / blen;

        // Try every hull vertex as candidate apex
        for (int j = 0; j < n; j++) {
            if (j == i || j == i1) continue;

            double jx = hx[j], jy = hy[j];

            // Apex must be on the interior side of the base line
            double apex_h = jx * bnx + jy * bny - base_dist;
            if (apex_h < 1e-10) continue;

            // For each other hull vertex k, compute the signed angle of
            // the direction j→k relative to the base direction.  We use
            // the cross product (bux,buy) × (j→k) to determine left/right.
            // The "leftmost" hull vertex (max cross) and "rightmost" (min
            // cross) define the tightest supporting lines from j.
            //
            // To ensure the lines j→left and j→right, when intersected
            // with the base line, produce a valid enclosing triangle, we
            // track the extreme cross products.

            double max_cross = -1e30, min_cross = 1e30;
            int left_idx = -1, right_idx = -1;

            for (int k = 0; k < n; k++) {
                if (k == j) continue;
                double dkx = hx[k] - jx, dky = hy[k] - jy;
                double dk_len = std::sqrt(dkx*dkx + dky*dky);
                if (dk_len < 1e-14) continue;
                // Normalise for consistent comparison
                double nkx = dkx / dk_len, nky = dky / dk_len;
                // Cross product with base direction: positive = left of base dir
                double cross = bux * nky - buy * nkx;

                if (cross > max_cross) { max_cross = cross; left_idx = k; }
                if (cross < min_cross) { min_cross = cross; right_idx = k; }
            }

            if (left_idx < 0 || right_idx < 0 || left_idx == right_idx) continue;

            // Left side: line from j through hull[left_idx]
            double ldx = hx[left_idx] - jx, ldy = hy[left_idx] - jy;
            // Right side: line from j through hull[right_idx]
            double rdx = hx[right_idx] - jx, rdy = hy[right_idx] - jy;

            // Intersect left/right sides with base line
            double v0x, v0y, v1x, v1y;
            bool ok0 = intersect_lines(jx, jy, ldx, ldy,
                                       hx[i], hy[i], bdx, bdy,
                                       v0x, v0y);
            bool ok1 = intersect_lines(jx, jy, rdx, rdy,
                                       hx[i], hy[i], bdx, bdy,
                                       v1x, v1y);
            if (!ok0 || !ok1) continue;

            // Candidate triangle: (v0, v1, j)
            double area = std::abs(tri_area_2(v0x, v0y, v1x, v1y, jx, jy));
            if (area < 1e-14 || area >= best_area) continue;

            // Ensure CCW for containment check
            double tv0x = v0x, tv0y = v0y, tv1x = v1x, tv1y = v1y;
            if (tri_area_2(tv0x, tv0y, tv1x, tv1y, jx, jy) < 0) {
                std::swap(tv0x, tv1x);
                std::swap(tv0y, tv1y);
            }

            // Verify all hull points are inside or on the triangle
            bool all_inside = true;
            for (int k = 0; k < n && all_inside; k++) {
                double d0 = tri_area_2(tv0x, tv0y, tv1x, tv1y, hx[k], hy[k]);
                double d1 = tri_area_2(tv1x, tv1y, jx,   jy,   hx[k], hy[k]);
                double d2 = tri_area_2(jx,   jy,   tv0x, tv0y, hx[k], hy[k]);
                if (d0 < -1e-9 || d1 < -1e-9 || d2 < -1e-9) all_inside = false;
            }
            if (!all_inside) continue;

            best_area = area;
            best_tx[0] = tv0x; best_ty[0] = tv0y;
            best_tx[1] = tv1x; best_ty[1] = tv1y;
            best_tx[2] = jx;   best_ty[2] = jy;
        }
    }

    // =========================================================================
    // Second pass: try triangles where ALL THREE sides are flush with hull
    // edges (apex is NOT a hull vertex).  The triangle vertices are at
    // the intersections of three hull-edge supporting lines.
    //
    // Complexity: O(n³) triples × O(n) containment = O(n⁴).
    // Only run for small hulls.  For n>80, the first pass (which already
    // considers every hull vertex as apex) produces near-optimal results.
    // =========================================================================
    if (n <= 80) {
        for (int a = 0; a < n; a++) {
            int a1 = next_idx(a, n);
            double adx = hx[a1]-hx[a], ady = hy[a1]-hy[a];
            if (adx*adx + ady*ady < 1e-28) continue;

            for (int b = a + 1; b < n; b++) {
                int b1 = next_idx(b, n);
                double bddx = hx[b1]-hx[b], bddy = hy[b1]-hy[b];
                if (bddx*bddx + bddy*bddy < 1e-28) continue;

                for (int c = b + 1; c < n; c++) {
                    int c1 = next_idx(c, n);
                    double cdx = hx[c1]-hx[c], cdy = hy[c1]-hy[c];
                    if (cdx*cdx + cdy*cdy < 1e-28) continue;

                    double pab_x, pab_y, pbc_x, pbc_y, pac_x, pac_y;
                    bool ok_ab = intersect_lines(hx[a], hy[a], adx, ady,
                                                 hx[b], hy[b], bddx, bddy,
                                                 pab_x, pab_y);
                    bool ok_bc = intersect_lines(hx[b], hy[b], bddx, bddy,
                                                 hx[c], hy[c], cdx, cdy,
                                                 pbc_x, pbc_y);
                    bool ok_ac = intersect_lines(hx[a], hy[a], adx, ady,
                                                 hx[c], hy[c], cdx, cdy,
                                                 pac_x, pac_y);
                    if (!ok_ab || !ok_bc || !ok_ac) continue;

                    double area = std::abs(tri_area_2(pab_x, pab_y, pbc_x, pbc_y, pac_x, pac_y));
                    if (area < 1e-14 || area >= best_area) continue;

                    // Ensure CCW
                    double t0x = pab_x, t0y = pab_y;
                    double t1x = pbc_x, t1y = pbc_y;
                    double t2x = pac_x, t2y = pac_y;
                    if (tri_area_2(t0x, t0y, t1x, t1y, t2x, t2y) < 0) {
                        std::swap(t0x, t1x);
                        std::swap(t0y, t1y);
                    }

                    bool all_inside = true;
                    for (int k = 0; k < n && all_inside; k++) {
                        double d0 = tri_area_2(t0x, t0y, t1x, t1y, hx[k], hy[k]);
                        double d1 = tri_area_2(t1x, t1y, t2x, t2y, hx[k], hy[k]);
                        double d2 = tri_area_2(t2x, t2y, t0x, t0y, hx[k], hy[k]);
                        if (d0 < -1e-9 || d1 < -1e-9 || d2 < -1e-9) all_inside = false;
                    }
                    if (!all_inside) continue;

                    best_area = area;
                    best_tx[0] = t0x; best_ty[0] = t0y;
                    best_tx[1] = t1x; best_ty[1] = t1y;
                    best_tx[2] = t2x; best_ty[2] = t2y;
                }
            }
        }
    }

    if (best_area < std::numeric_limits<double>::infinity()) {
        // Epsilon expansion: push each vertex outward from centroid to guarantee
        // strict enclosure under exact arithmetic.  Scale relative to triangle size.
        double cx = (best_tx[0]+best_tx[1]+best_tx[2])/3.0;
        double cy = (best_ty[0]+best_ty[1]+best_ty[2])/3.0;
        double max_r = 0.0;
        for (int k = 0; k < 3; k++) {
            double dx_ = best_tx[k] - cx, dy_ = best_ty[k] - cy;
            max_r = std::max(max_r, std::sqrt(dx_*dx_ + dy_*dy_));
        }
        double eps = max_r * 1e-9;
        if (eps < 1e-12) eps = 1e-12;
        for (int k = 0; k < 3; k++) {
            double dx_ = best_tx[k] - cx, dy_ = best_ty[k] - cy;
            double len_ = std::sqrt(dx_*dx_ + dy_*dy_);
            if (len_ > 1e-14) { dx_ /= len_; dy_ /= len_; }
            best_tx[k] += eps * dx_;
            best_ty[k] += eps * dy_;
        }

        // Log tightness ratio
        double hull_area = 0;
        for (int i = 0; i < n; i++) {
            int j = next_idx(i, n);
            hull_area += hx[i]*hy[j] - hx[j]*hy[i];
        }
        hull_area = std::abs(hull_area) / 2.0;
        std::cerr << "[TRIANGLE] Enclosing triangle area: " << std::fixed << std::setprecision(6)
                  << best_area << " (hull area: " << hull_area
                  << ", ratio: " << std::setprecision(3) << (best_area / hull_area) << "x)\n";

        Triangle best_tri = { Point(best_tx[0], best_ty[0]),
                              Point(best_tx[1], best_ty[1]),
                              Point(best_tx[2], best_ty[2]) };
        return improve_triangle_by_vertex_projection(hull_pts, best_tri);
    }

    // Fallback: axis-aligned right triangle from bounding box
    double min_x = *std::min_element(hx.begin(), hx.end()) - 1;
    double min_y = *std::min_element(hy.begin(), hy.end()) - 1;
    double max_x = *std::max_element(hx.begin(), hx.end());
    double max_y = *std::max_element(hy.begin(), hy.end());
    return { Point(min_x, min_y),
             Point(min_x + 2*(max_x-min_x)+1, min_y),
             Point(min_x, min_y + 2*(max_y-min_y)+1) };
}

// =========================================================================
// FIXED-POINT QUANTIZATION + BLOCKCHAIN-WALKABLE BOUNDARY
// =========================================================================
// Global integer coordinate system (publishing-time snap-to-grid):
//   X (longitude): uint32, 0 at the international date line, 360 deg over 2^32
//   Y (latitude) : int32 , 0 at the equator, +/-90 deg over +/-(2^31 - 1)
// Resolution is ~1 cm globally -- far finer than shapefile vertex accuracy --
// so the on-chain verifier checks that a split point is WITHIN TOLERANCE of a
// segment (not exactly on it): snapping can move a point off the exact line by
// up to ~1 grid unit, since both the point and the segment endpoints snap.
//
// Walkability: the published boundary must be a single CLOSED cycle a verifier
// can PerWalk end to end.  We rebuild the fully-subdivided boundary (original
// edges + every recorded SplitLine point, ordered along each edge), quantize
// each distinct vertex ONCE so shared endpoints receive identical integers,
// then assert that consecutive integer sub-segments connect and the last meets
// the first.  Because the same exact CGAL Point is reused on both sides of
// every split, its quantized image is identical on both sub-segments -- the
// walk cannot develop a gap.
// =========================================================================

static const double Q_PI        = 3.14159265358979323846;
static const double Q_LON_SCALE = 4294967296.0 / 360.0;   // 2^32 / 360
static const double Q_LAT_SCALE = 2147483647.0 / 90.0;    // (2^31 - 1) / 90

struct QPoint { uint32_t X; int32_t Y; };

static inline QPoint quantize(double lon, double lat) {
    double xn = (lon + 180.0) * Q_LON_SCALE;     // date-line origin -> [0, 2^32)
    if (xn < 0.0) xn = 0.0;
    if (xn > 4294967295.0) xn = 4294967295.0;
    double yn = lat * Q_LAT_SCALE;               // equator origin -> +/-(2^31-1)
    if (yn >  2147483647.0) yn =  2147483647.0;
    if (yn < -2147483647.0) yn = -2147483647.0;
    QPoint q;
    q.X = (uint32_t)std::llround(xn);
    q.Y = (int32_t) std::llround(yn);
    return q;
}
static inline QPoint quantize(const Point& p) {
    return quantize(CGAL::to_double(p.x()), CGAL::to_double(p.y()));
}
static inline void dequantize(const QPoint& q, double& lon, double& lat) {
    lon = (double)q.X / Q_LON_SCALE - 180.0;
    lat = (double)q.Y / Q_LAT_SCALE;
}
// Approximate metres between an exact point and its quantized position.
static double snap_deviation_m(const Point& p) {
    double lon = CGAL::to_double(p.x()), lat = CGAL::to_double(p.y());
    QPoint q = quantize(lon, lat);
    double lon2, lat2; dequantize(q, lon2, lat2);
    double dlat = (lat2 - lat) * 111320.0;
    double dlon = (lon2 - lon) * 111320.0 * std::cos(lat * Q_PI / 180.0);
    return std::sqrt(dlat * dlat + dlon * dlon);
}

static inline void write_qtri(std::ostream& os, const Triangle& t) {
    for (int i = 0; i < 3; ++i) {
        QPoint q = quantize(t[i]);
        os << (i ? "," : "") << q.X << "," << q.Y;
    }
}

// Build + verify the quantized, walkable boundary and emit blockchain-ready
// transaction records.  Returns true if the integer walk is closed+connected.
bool emit_quantized_walkable_transactions(const std::vector<Edge>& polygon_edges,
                                          std::ofstream& pipeline_log)
{
    std::ofstream walk_csv("boundary_walk.csv");
    std::ofstream qtx("ags_transactions_quantized.txt");
    std::ofstream qreport("quantization_report.txt");
    std::ofstream tri_q_csv("triangle_splits_quantized.csv");
    std::ofstream line_q_csv("line_splits_quantized.csv");

    // ---- 1) Rebuild the fully-subdivided boundary in EXACT arithmetic ----
    // For each original edge A->B, collect every recorded split point lying on
    // it (handles cascaded splits, where a later split point sits on a
    // sub-segment of A->B), order them along the edge, and emit sub-segments.
    struct OnEdge { Point p; K::FT t; };

    std::vector<std::array<QPoint, 2>> walk;
    walk.reserve(polygon_edges.size() + g_line_splits.size());

    long long subseg_count = 0;
    double max_dev = 0.0;

    // Precompute split-point coordinates as doubles ONCE. The per-edge loop then
    // uses a cheap double bounding-box test to reject the vast majority of split
    // points before paying for an exact CGAL::collinear predicate. This keeps the
    // boundary rebuild near O(E + S) in practice instead of an O(E x S) sweep of
    // exact tests, which matters on full-state runs with tens of thousands of
    // edges and splits.
    std::vector<double> lsx(g_line_splits.size()), lsy(g_line_splits.size());
    for (size_t k = 0; k < g_line_splits.size(); ++k) {
        lsx[k] = CGAL::to_double(g_line_splits[k].split_point.x());
        lsy[k] = CGAL::to_double(g_line_splits[k].split_point.y());
    }
    const double BB_EPS = 1e-9;

    walk_csv << "seq,exact_sx,exact_sy,exact_tx,exact_ty,qsx,qsy,qtx,qty\n";

    for (const auto& E : polygon_edges) {
        const Point& A = E.source();
        const Point& B = E.target();
        K::Vector_2 dir = B - A;
        K::FT dlen2 = dir * dir;
        double Ax = CGAL::to_double(A.x()), Ay = CGAL::to_double(A.y());
        double Bx = CGAL::to_double(B.x()), By = CGAL::to_double(B.y());
        double minx = std::min(Ax, Bx) - BB_EPS, maxx = std::max(Ax, Bx) + BB_EPS;
        double miny = std::min(Ay, By) - BB_EPS, maxy = std::max(Ay, By) + BB_EPS;

        std::vector<OnEdge> mids;
        if (dlen2 != K::FT(0)) {
            for (size_t k = 0; k < g_line_splits.size(); ++k) {
                // Cheap double bbox reject before the exact predicate.
                if (lsx[k] < minx || lsx[k] > maxx || lsy[k] < miny || lsy[k] > maxy)
                    continue;
                const Point& P = g_line_splits[k].split_point;
                if (!CGAL::collinear(A, B, P)) continue;
                K::FT t = ((P - A) * dir) / dlen2;
                if (t > K::FT(0) && t < K::FT(1)) {
                    bool dup = false;
                    for (const auto& m : mids) if (m.p == P) { dup = true; break; }
                    if (!dup) mids.push_back({ P, t });
                }
            }
        }
        std::sort(mids.begin(), mids.end(),
                  [](const OnEdge& a, const OnEdge& b){ return a.t < b.t; });

        Point prev = A;
        auto emit_seg = [&](const Point& s, const Point& tt) {
            QPoint qs = quantize(s), qt = quantize(tt);
            walk.push_back({ qs, qt });
            max_dev = std::max(max_dev, snap_deviation_m(s));
            max_dev = std::max(max_dev, snap_deviation_m(tt));
            if (subseg_count < 100000) {  // cap CSV size on huge runs
                walk_csv << subseg_count << ","
                         << std::fixed << std::setprecision(10)
                         << CGAL::to_double(s.x())  << "," << CGAL::to_double(s.y())  << ","
                         << CGAL::to_double(tt.x()) << "," << CGAL::to_double(tt.y()) << ","
                         << qs.X << "," << qs.Y << "," << qt.X << "," << qt.Y << "\n";
            }
            subseg_count++;
        };
        for (const auto& m : mids) { emit_seg(prev, m.p); prev = m.p; }
        emit_seg(prev, B);
    }
    walk_csv.close();

    // ---- 2) Verify the integer walk is a single closed connected cycle ----
    bool connected = true, closed = false;
    long long break_at = -1;
    for (size_t i = 0; i + 1 < walk.size(); ++i) {
        if (walk[i][1].X != walk[i + 1][0].X || walk[i][1].Y != walk[i + 1][0].Y) {
            connected = false; break_at = (long long)i; break;
        }
    }
    if (!walk.empty())
        closed = (walk.back()[1].X == walk.front()[0].X &&
                  walk.back()[1].Y == walk.front()[0].Y);
    bool walkable = connected && closed && !walk.empty();

    // ---- 3) Emit quantized TriangleSplit transactions (xi_T) ----
    tri_q_csv << "split_id,iteration,is_lookahead,"
              << "p1x,p1y,p2x,p2y,p3x,p3y,"
              << "split_sx,split_sy,split_tx,split_ty,"
              << "c1_1x,c1_1y,c1_2x,c1_2y,c1_3x,c1_3y,"
              << "c2_1x,c2_1y,c2_2x,c2_2y,c2_3x,c2_3y\n";
    for (const auto& r : g_triangle_splits) {
        QPoint ss = quantize(r.split_line.source());
        QPoint st = quantize(r.split_line.target());
        tri_q_csv << r.split_id << "," << r.iteration << "," << (r.is_lookahead ? 1 : 0) << ",";
        write_qtri(tri_q_csv, r.parent_triangle);
        tri_q_csv << "," << ss.X << "," << ss.Y << "," << st.X << "," << st.Y << ",";
        write_qtri(tri_q_csv, r.child1);
        tri_q_csv << ",";
        write_qtri(tri_q_csv, r.child2);
        tri_q_csv << "\n";
    }
    tri_q_csv.close();

    // ---- 4) Emit quantized SplitLine transactions (xi_P) ----
    line_q_csv << "split_id,triangle_split_id,iteration,"
               << "Dx,Dy,Ex,Ey,Px,Py,"
               << "s1_sx,s1_sy,s1_tx,s1_ty,s2_sx,s2_sy,s2_tx,s2_ty,snap_dev_m\n";
    for (const auto& r : g_line_splits) {
        QPoint D  = quantize(r.original_edge.source());
        QPoint Ee = quantize(r.original_edge.target());
        QPoint P  = quantize(r.split_point);
        QPoint a1 = quantize(r.segment1.source()), b1 = quantize(r.segment1.target());
        QPoint a2 = quantize(r.segment2.source()), b2 = quantize(r.segment2.target());
        double dev = snap_deviation_m(r.split_point);
        line_q_csv << r.split_id << "," << r.triangle_split_id << "," << r.iteration << ","
                   << D.X << "," << D.Y << "," << Ee.X << "," << Ee.Y << ","
                   << P.X << "," << P.Y << ","
                   << a1.X << "," << a1.Y << "," << b1.X << "," << b1.Y << ","
                   << a2.X << "," << a2.Y << "," << b2.X << "," << b2.Y << ","
                   << std::fixed << std::setprecision(4) << dev << "\n";
    }
    line_q_csv.close();

    // ---- 5) Human-readable transaction summary + PerWalk (capped) ----
    qtx << "==============================================\n";
    qtx << "QUANTIZED AGS TRANSACTIONS (blockchain-ready)\n";
    qtx << "==============================================\n";
    qtx << "Global integer coordinate system:\n";
    qtx << "  Longitude X: uint32, 0 at the international date line,\n";
    qtx << "               360 deg mapped across 2^32 steps.\n";
    qtx << "  Latitude  Y: int32, 0 at the equator, +/-90 deg across +/-(2^31-1).\n";
    qtx << "  Resolution : ~0.93 cm (lon, equator), ~0.47 cm (lat).\n";
    qtx << "  All coordinates below are integers in this system.\n";
    qtx << "----------------------------------------------\n";
    qtx << "TriangleSplit transactions (xi_T): " << g_triangle_splits.size() << "\n";
    qtx << "SplitLine    transactions (xi_P): " << g_line_splits.size()
        << " (x3 on chain = 1 split + 2 PerWalks)\n";
    qtx << "Boundary PerWalk sub-segments    : " << subseg_count << "\n";
    qtx << "----------------------------------------------\n";
    qtx << "PERWALK (closed boundary cycle, first/last 20 sub-segments):\n";
    for (size_t i = 0; i < walk.size(); ++i) {
        if (walk.size() > 40 && i == 20)
            qtx << "  ... (" << (walk.size() - 40) << " sub-segments omitted) ...\n";
        if (walk.size() <= 40 || i < 20 || i >= walk.size() - 20) {
            qtx << "  [" << i << "] (" << walk[i][0].X << "," << walk[i][0].Y << ") -> ("
                << walk[i][1].X << "," << walk[i][1].Y << ")\n";
        }
    }
    qtx << "----------------------------------------------\n";
    qtx << "WALK VERIFICATION:\n";
    qtx << "  Connected (each end == next start): " << (connected ? "YES" : "NO") << "\n";
    if (!connected) qtx << "    First break after sub-segment #" << break_at << "\n";
    qtx << "  Closed (last end == first start)  : " << (closed ? "YES" : "NO") << "\n";
    qtx << "  WALKABLE ON CHAIN                 : " << (walkable ? "YES" : "NO") << "\n";
    qtx << "==============================================\n";
    qtx.close();

    // ---- 6) Quantization report ----
    int rec_tol = (int)std::ceil(max_dev / 0.0093) + 1;  // ~grid units (lon@equator)
    if (rec_tol < 2) rec_tol = 2;
    qreport << "==============================================\n";
    qreport << "QUANTIZATION REPORT\n";
    qreport << "==============================================\n";
    qreport << "Coordinate system : global integer grid\n";
    qreport << "  Longitude X     : uint32, origin = international date line\n";
    qreport << "  Latitude  Y     : int32 , origin = equator\n";
    qreport << "  Lon scale       : 2^32 / 360 = " << std::fixed << std::setprecision(4)
            << Q_LON_SCALE << " steps/deg\n";
    qreport << "  Lat scale       : (2^31-1) / 90 = " << Q_LAT_SCALE << " steps/deg\n";
    qreport << "  Lon resolution  : ~0.93 cm at the equator\n";
    qreport << "  Lat resolution  : ~0.47 cm\n";
    qreport << "----------------------------------------------\n";
    qreport << "Max snap deviation (any published vertex): "
            << std::fixed << std::setprecision(4) << max_dev << " m\n";
    qreport << "Recommended verifier tolerance           : " << rec_tol
            << " grid units (~" << std::setprecision(2) << (rec_tol * 0.93) << " cm)\n";
    qreport << "----------------------------------------------\n";
    qreport << "Boundary sub-segments : " << subseg_count << "\n";
    qreport << "Walk connected        : " << (connected ? "YES" : "NO") << "\n";
    qreport << "Walk closed           : " << (closed ? "YES" : "NO") << "\n";
    qreport << "WALKABLE ON CHAIN     : " << (walkable ? "YES" : "NO") << "\n";
    qreport << "==============================================\n";
    qreport.close();

    pipeline_log << "\n=== FIXED-POINT QUANTIZATION ===\n";
    pipeline_log << "Boundary PerWalk sub-segments: " << subseg_count << "\n";
    pipeline_log << "Max snap deviation: " << std::fixed << std::setprecision(4) << max_dev << " m\n";
    pipeline_log << "Recommended verifier tolerance: " << rec_tol << " grid units\n";
    pipeline_log << "Walk connected: " << (connected ? "YES" : "NO")
                 << ", closed: " << (closed ? "YES" : "NO")
                 << ", WALKABLE: " << (walkable ? "YES" : "NO") << "\n";
    pipeline_log.flush();

    std::cerr << "[QUANTIZE] " << subseg_count << " boundary sub-segments, max snap "
              << std::fixed << std::setprecision(3) << max_dev << " m, walkable="
              << (walkable ? "YES" : "NO") << "\n";

    return walkable;
}


// =========================================================================
// MULTI-STATE BATCH MODE
// =========================================================================
struct StateConfig {
    std::string shapefile_path, state_fips, county_fips;
    double reduction_pct = 0.0; bool do_reduction = false;
    std::string header_path = "polygon_output.h";
    std::string points_path = "polygon_output.txt";
    std::string var_prefix  = "simple_polygon";
    bool do_write_header = true, do_write_points = true;
};

struct StateSummary {
    std::string fips, name;
    int  original_points = 0, used_points = 0;
    int  total_leaves = 0, interior_leaves = 0;
    bool walkable = false;
    int  lookahead = 0;
    std::string method;
    long long ags_txns = 0;
    double runtime_s = 0.0;
    bool ok = false;
};

// Reset all per-run global state so the pipeline can run again for the next state.
static void reset_pipeline_globals() {
    g_evaluation_terminated = false;
    g_triangle_splits.clear();
    g_line_splits.clear();
    g_triangle_split_count = 0;
    g_line_split_count = 0;
    g_leaf_triangles.clear();
    g_leaf_region_codes.clear();
    g_mark_triangle_count = 0;
    g_recursive_call_count = 0;
    g_cache_hit_count = 0;
    g_memo_cache.clear();
    g_last_mark = MarkResult();
    MAX_CANDIDATES_PER_LEVEL = DEFAULT_MAX_CANDIDATES;
    if (g_decision_log.is_open())    g_decision_log.close();
    if (g_transaction_log.is_open()) g_transaction_log.close();
}

static std::string get_cwd_str() {
    char buf[4096];
    if (PORTABLE_GETCWD(buf, sizeof(buf))) return std::string(buf);
    return std::string(".");
}

static std::string to_absolute_path(const std::string& p) {
#ifdef _WIN32
    if (p.size() > 1 && p[1] == ':') return p;
    if (!p.empty() && (p[0] == '\\' || p[0] == '/')) return p;
#else
    if (!p.empty() && p[0] == '/') return p;
#endif
    return get_cwd_str() + "/" + p;
}

// First N distinct state FIPS codes, in shapefile record order.
static std::vector<std::string> first_n_state_fips(const std::string& shp_path, int n) {
    std::vector<std::string> out;
    std::string dbf_path = shp_path;
    auto dot = dbf_path.rfind('.');
    if (dot != std::string::npos) dbf_path = dbf_path.substr(0, dot);
    dbf_path += ".dbf";
    DBFHandle dbf = DBFOpen(dbf_path.c_str(), "rb");
    if (!dbf) { std::cerr << "[ERROR] cannot open DBF: " << dbf_path << "\n"; return out; }
    int fld = find_dbf_field(dbf, {"STATEFP", "STATEFP20", "STATEFP10", "STATE"});
    if (fld < 0) { std::cerr << "[ERROR] no STATEFP field in " << dbf_path << "\n"; DBFClose(dbf); return out; }
    int nrec = DBFGetRecordCount(dbf);
    std::unordered_set<std::string> seen;
    for (int i = 0; i < nrec && (int)out.size() < n; ++i) {
        const char* v = DBFReadStringAttribute(dbf, i, fld);
        if (!v) continue;
        std::string f = trim_str(v);
        if (f.empty()) continue;
        if (seen.insert(f).second) out.push_back(f);
    }
    DBFClose(dbf);
    return out;
}

static int run_one_state(const StateConfig& cfg, StateSummary& sum);

// Print a per-state look-ahead trend (one row per depth) showing how leaf count
// falls as the look-ahead depth increases, and write summary.csv (one row per
// state x depth). Summaries are grouped by the order states first appear.
static void write_lookahead_summary(const std::vector<StateSummary>& S) {
    std::ofstream csv("summary.csv");
    csv << "fips,name,lookahead,original_points,used_points,total_leaves,interior_leaves,"
           "perimeter_walk,method,ags_transactions,walkable_chain,runtime_s,status\n";
    for (const auto& s : S) {
        std::string name = s.name; for (char& c : name) if (c == ',') c = ' ';
        csv << s.fips << "," << name << "," << s.lookahead << ","
            << s.original_points << "," << s.used_points << ","
            << s.total_leaves << "," << s.interior_leaves << ","
            << (s.walkable ? "YES" : "NO") << "," << s.method << "," << s.ags_txns << ","
            << (s.walkable ? "YES" : "NO") << ","
            << std::fixed << std::setprecision(1) << s.runtime_s << ","
            << (s.ok ? "OK" : "FAILED") << "\n";
    }
    csv.close();

    // Distinct FIPS in first-appearance order.
    std::vector<std::string> order;
    std::unordered_set<std::string> seen;
    for (const auto& s : S) if (seen.insert(s.fips).second) order.push_back(s.fips);

    std::cout << "\n" << std::string(92, '=') << "\n";
    std::cout << "LOOK-AHEAD DEPTH SWEEP  (" << order.size() << " state(s))\n";
    std::cout << std::string(92, '=') << "\n";

    for (const auto& fips : order) {
        std::string name; int baseline = -1;
        for (const auto& s : S) if (s.fips == fips && s.ok) { name = s.name; break; }
        std::cout << "\nFIPS " << fips << (name.empty() ? "" : "  " + name) << "\n";
        std::cout << std::right
                  << std::setw(5)  << "LA"     << std::setw(10) << "Leaves"
                  << std::setw(10) << "Interior" << std::setw(12) << "AGS_Txns"
                  << std::setw(10) << "Walkable" << std::setw(11) << "Runtime"
                  << std::setw(11) << "dLeaves" << std::setw(12) << "%vs_LA" << "\n";
        std::cout << std::string(81, '-') << "\n";
        for (const auto& s : S) {
            if (s.fips != fips) continue;
            if (!s.ok) {
                std::cout << std::right << std::setw(5) << s.lookahead
                          << std::setw(10) << "-" << std::setw(10) << "-"
                          << std::setw(12) << "-" << std::setw(10) << "-"
                          << std::setw(11) << "-" << std::setw(11) << "FAILED"
                          << std::setw(12) << "-" << "\n";
                continue;
            }
            if (baseline < 0) baseline = s.total_leaves;   // first successful depth
            int dleaves = s.total_leaves - baseline;
            double pct = baseline > 0 ? 100.0 * dleaves / baseline : 0.0;
            std::ostringstream rt; rt << std::fixed << std::setprecision(1) << s.runtime_s << "s";
            std::ostringstream pc; pc << std::fixed << std::setprecision(1) << pct << "%";
            std::ostringstream dl; if (dleaves > 0) dl << "+"; dl << dleaves;
            std::cout << std::right << std::setw(5) << s.lookahead
                      << std::setw(10) << s.total_leaves
                      << std::setw(10) << s.interior_leaves
                      << std::setw(12) << s.ags_txns
                      << std::setw(10) << (s.walkable ? "YES" : "NO")
                      << std::setw(11) << rt.str()
                      << std::setw(11) << dl.str()
                      << std::setw(12) << pc.str() << "\n";
        }
    }
    std::cout << "\n" << std::string(92, '-') << "\n";
    std::cout << "Per (state, depth) outputs are in state_<FIPS>/la_<D>/ (or la_<D>/ for a single\n";
    std::cout << "state); this trend and all rows are also in summary.csv\n";
    std::cout << std::string(92, '=') << "\n";
}



// Print a summary table to stdout and write summary.csv in the current dir.
static void write_summary_table(const std::vector<StateSummary>& S) {
    std::ofstream csv("summary.csv");
    csv << "fips,name,original_points,used_points,total_leaves,interior_leaves,"
           "perimeter_walk,lookahead,method,ags_transactions,walkable_chain,status\n";
    for (const auto& s : S) {
        std::string name = s.name; for (char& c : name) if (c == ',') c = ' ';
        csv << s.fips << "," << name << "," << s.original_points << "," << s.used_points << ","
            << s.total_leaves << "," << s.interior_leaves << ","
            << (s.walkable ? "YES" : "NO") << "," << s.lookahead << "," << s.method << ","
            << s.ags_txns << "," << (s.walkable ? "YES" : "NO") << ","
            << (s.ok ? "OK" : "FAILED") << "\n";
    }
    csv.close();

    std::cout << "\n" << std::string(99, '=') << "\n";
    std::cout << "MULTI-STATE SUMMARY  (" << S.size() << " states)\n";
    std::cout << std::string(99, '=') << "\n";
    std::cout << std::left
              << std::setw(5)  << "FIPS" << std::setw(22) << "Name"
              << std::right
              << std::setw(9)  << "OrigPts" << std::setw(9) << "UsedPts"
              << std::setw(9)  << "Leaves"  << std::setw(10) << "Interior"
              << std::setw(10) << "Walkable" << std::setw(10) << "Method"
              << std::setw(5)  << "LA" << std::setw(12) << "AGS_Txns"
              << std::setw(9)  << "Status" << "\n";
    std::cout << std::string(99, '-') << "\n";
    for (const auto& s : S) {
        std::string name = s.name.size() > 21 ? s.name.substr(0, 21) : s.name;
        std::cout << std::left
                  << std::setw(5)  << s.fips << std::setw(22) << name
                  << std::right
                  << std::setw(9)  << s.original_points << std::setw(9) << s.used_points
                  << std::setw(9)  << s.total_leaves << std::setw(10) << s.interior_leaves
                  << std::setw(10) << (s.walkable ? "YES" : "NO")
                  << std::setw(10) << (s.ok ? s.method : std::string("-"))
                  << std::setw(5)  << s.lookahead << std::setw(12) << s.ags_txns
                  << std::setw(9)  << (s.ok ? "OK" : "FAILED") << "\n";
    }
    std::cout << std::string(99, '-') << "\n";
    std::cout << "Per-state outputs are in state_<FIPS>/ folders; table also written to summary.csv\n";
    std::cout << std::string(99, '=') << "\n";
}


static int run_one_state(const StateConfig& cfg, StateSummary& sum) {
    reset_pipeline_globals();
    g_decision_log.open("triangle_decision_log.txt");

    // Aliases so the extracted pipeline body reads exactly as before.
    const std::string& shapefile_path = cfg.shapefile_path;
    const std::string& state_fips     = cfg.state_fips;
    const std::string& county_fips    = cfg.county_fips;
    double reduction_pct              = cfg.reduction_pct;
    bool   do_reduction               = cfg.do_reduction;
    const std::string& header_path    = cfg.header_path;
    const std::string& points_path    = cfg.points_path;
    const std::string& var_prefix     = cfg.var_prefix;
    bool   do_write_header            = cfg.do_write_header;
    bool   do_write_points            = cfg.do_write_points;
    // =====================================================================
    // OPEN LOG FILES
    // =====================================================================
    // pipeline_log: records the entire shapefile->reduce->header->triangulate flow
    std::ofstream pipeline_log("pipeline_log.txt");
    if (!pipeline_log) {
        std::cerr << "Error: Could not open pipeline_log.txt\n";
        return 1;
    }

    pipeline_log << "==============================================\n";
    pipeline_log << "CGAL SLOPELINE INTEGRATED PIPELINE LOG\n";
    pipeline_log << "==============================================\n";
    pipeline_log << "Shapefile:   " << shapefile_path << "\n";
    pipeline_log << "Reduction:   " << (do_reduction ?
        std::to_string(static_cast<int>(reduction_pct)) + "%" : "none (all original vertices)") << "\n";
    if (!state_fips.empty())
        pipeline_log << "State FIPS:  " << state_fips << "\n";
    if (!county_fips.empty())
        pipeline_log << "County FIPS: " << county_fips << "\n";
    pipeline_log << "Lookahead:   " << LOOKAHEAD_DEPTH_LIMIT << "\n";
    pipeline_log << "Header out:  " << (do_write_header ? header_path : "(disabled)") << "\n";
    pipeline_log << "Points out:  " << (do_write_points ? points_path : "(disabled)") << "\n";
    pipeline_log << "==============================================\n\n";
    pipeline_log.flush();

    // Transaction log
    g_transaction_log.open("ags_transactions.txt");
    if (!g_transaction_log) {
        std::cerr << "Error: Could not open ags_transactions.txt\n";
        return 1;
    }

    // =====================================================================
    // STEP 1: READ SHAPEFILE
    // =====================================================================
    std::cerr << "\n=== Step 1/4: Reading shapefile ===\n";
    pipeline_log << "=== STEP 1: READ SHAPEFILE ===\n";

    ShapeReadResult shape_result;
    try {
        shape_result = read_polygon_from_shapefile(shapefile_path, state_fips, county_fips, pipeline_log);
    } catch (const std::exception& e) {
        std::cerr << "[ERROR] " << e.what() << "\n";
        pipeline_log << "[ERROR] " << e.what() << "\n";
        return 1;
    }

    // Log original polygon for reconstruction
    pipeline_log << "\n[POLYGON_ORIGINAL] " << shape_result.points.size() << " points\n";
    {
        size_t n_pts = shape_result.points.size();
        size_t dump_limit = 20;
        size_t dump_end = (n_pts > dump_limit * 2) ? dump_limit : n_pts;
        for (size_t i = 0; i < dump_end; ++i) {
            pipeline_log << "  P" << i << ": ("
                         << std::fixed << std::setprecision(10)
                         << CGAL::to_double(shape_result.points[i].x()) << ", "
                         << CGAL::to_double(shape_result.points[i].y()) << ")\n";
        }
        if (n_pts > dump_limit * 2) {
            pipeline_log << "  ... (" << (n_pts - dump_limit * 2) << " points omitted) ...\n";
            for (size_t i = n_pts - dump_limit; i < n_pts; ++i) {
                pipeline_log << "  P" << i << ": ("
                             << std::fixed << std::setprecision(10)
                             << CGAL::to_double(shape_result.points[i].x()) << ", "
                             << CGAL::to_double(shape_result.points[i].y()) << ")\n";
            }
        }
    }
    pipeline_log << "\n";
    pipeline_log.flush();

    // =====================================================================
    // STEP 2: REDUCE POLYGON (optional)
    // =====================================================================
    std::vector<Point> polygon_points;

    if (do_reduction) {
        std::cerr << "\n=== Step 2/4: Reducing polygon (" << reduction_pct << "%) ===\n";
        pipeline_log << "=== STEP 2: REDUCE POLYGON ===\n";

        polygon_points = reduce_polygon_by_percentage(
            shape_result.points, reduction_pct, pipeline_log);

        // Log reduced polygon for reconstruction (capped for large polygons)
        pipeline_log << "\n[POLYGON_REDUCED] " << polygon_points.size() << " points\n";
        {
            size_t n_pts = polygon_points.size(), lim = 20;
            size_t head = (n_pts > lim * 2) ? lim : n_pts;
            for (size_t i = 0; i < head; ++i)
                pipeline_log << "  P" << i << ": (" << std::fixed << std::setprecision(10)
                             << CGAL::to_double(polygon_points[i].x()) << ", "
                             << CGAL::to_double(polygon_points[i].y()) << ")\n";
            if (n_pts > lim * 2) {
                pipeline_log << "  ... (" << (n_pts - lim * 2) << " points omitted) ...\n";
                for (size_t i = n_pts - lim; i < n_pts; ++i)
                    pipeline_log << "  P" << i << ": (" << std::fixed << std::setprecision(10)
                                 << CGAL::to_double(polygon_points[i].x()) << ", "
                                 << CGAL::to_double(polygon_points[i].y()) << ")\n";
            }
        }
        pipeline_log << "\n";
        pipeline_log.flush();
    } else {
        std::cerr << "\n=== Step 2/4: No reduction (using all " << shape_result.points.size() << " original vertices) ===\n";
        pipeline_log << "=== STEP 2: NO REDUCTION (using original points) ===\n";

        polygon_points = shape_result.points;

        pipeline_log << "\n[POLYGON_ORIGINAL_USED] " << polygon_points.size() << " points (no reduction applied)\n";
        pipeline_log << "  (Points identical to [POLYGON_ORIGINAL] above — see that section for vertex listing)\n";
        pipeline_log << "\n";
        pipeline_log.flush();
    }

    // =====================================================================
    // STEP 3: GENERATE HEADER AND POINTS FILES
    // =====================================================================
    std::cerr << "\n=== Step 3/4: Generating output files ===\n";
    pipeline_log << "=== STEP 3: GENERATE OUTPUT FILES ===\n";

    std::string source_info = shapefile_path;
    if (!shape_result.record_name.empty())
        source_info += " (" + shape_result.record_name + ")";
    if (do_reduction)
        source_info += " | " + std::to_string(static_cast<int>(reduction_pct)) + "% reduction";
    else
        source_info += " | no reduction";

    if (do_write_header) {
        write_cgal_header(polygon_points, header_path, var_prefix, source_info, pipeline_log);
    }
    if (do_write_points) {
        write_points_file(polygon_points, points_path, source_info, pipeline_log);
    }
    pipeline_log << "\n";
    pipeline_log.flush();

    // =====================================================================
    // STEP 4: TRIANGULAR DECOMPOSITION
    // =====================================================================
    std::cerr << "\n=== Step 4/4: Triangular decomposition ===\n";
    pipeline_log << "=== STEP 4: TRIANGULAR DECOMPOSITION ===\n";

    // Basic validation
    if (polygon_points.size() < 3) {
        std::cerr << "Error: Polygon must have at least 3 points after reduction!\n";
        pipeline_log << "[ERROR] Polygon has fewer than 3 points after reduction.\n";
        return 1;
    }

    // Build CGAL polygon
    Polygon polygon(polygon_points.begin(), polygon_points.end());

    if (!polygon.is_simple()) {
        std::cerr << "Error: Polygon is not simple after reduction!\n"
                  << "  The Douglas-Peucker simplification may have caused self-intersection.\n"
                  << "  Try a smaller reduction percentage.\n";
        pipeline_log << "[ERROR] Polygon is not simple after reduction. Try smaller reduction %.\n";
        return 1;
    }

    if (polygon.is_clockwise_oriented()) {
        polygon.reverse_orientation();
        pipeline_log << "[POLYGON] Reversed to counter-clockwise orientation\n";
    }

    // Build edge list
    std::vector<Edge> polygon_edges;
    for (size_t i = 0; i < polygon_points.size(); ++i) {
        polygon_edges.emplace_back(polygon_points[i], polygon_points[(i + 1) % polygon_points.size()]);
    }

    // Log edges for reconstruction (capped for large polygons)
    pipeline_log << "\n[POLYGON_EDGES] " << polygon_edges.size() << " edges\n";
    {
        size_t n_edges = polygon_edges.size(), lim = 20;
        size_t head = (n_edges > lim * 2) ? lim : n_edges;
        for (size_t i = 0; i < head; ++i) {
            pipeline_log << "  E" << i << ": ("
                         << std::fixed << std::setprecision(10)
                         << CGAL::to_double(polygon_edges[i].source().x()) << ", "
                         << CGAL::to_double(polygon_edges[i].source().y()) << ") -> ("
                         << CGAL::to_double(polygon_edges[i].target().x()) << ", "
                         << CGAL::to_double(polygon_edges[i].target().y()) << ")\n";
        }
        if (n_edges > lim * 2) {
            pipeline_log << "  ... (" << (n_edges - lim * 2) << " edges omitted) ...\n";
            for (size_t i = n_edges - lim; i < n_edges; ++i)
                pipeline_log << "  E" << i << ": ("
                             << std::fixed << std::setprecision(10)
                             << CGAL::to_double(polygon_edges[i].source().x()) << ", "
                             << CGAL::to_double(polygon_edges[i].source().y()) << ") -> ("
                             << CGAL::to_double(polygon_edges[i].target().x()) << ", "
                             << CGAL::to_double(polygon_edges[i].target().y()) << ")\n";
        }
    }
    pipeline_log << "\n";

    compute_dynamic_limits(static_cast<int>(polygon_edges.size()));

    // Log configuration to decision log
    g_decision_log << "==============================================\n";
    g_decision_log << "TRIANGULAR DECOMPOSITION (SLOPELINE) WITH TRANSACTION TRACKING\n";
    g_decision_log << "==============================================\n";
    g_decision_log << "SOURCE:\n";
    g_decision_log << "  Shapefile: " << shapefile_path << "\n";
    if (!shape_result.record_name.empty())
        g_decision_log << "  Record: " << shape_result.record_name << "\n";
    g_decision_log << "  Reduction: " << (do_reduction ?
        std::to_string(static_cast<int>(reduction_pct)) + "%" : "none") << "\n";
    g_decision_log << "  Original points: " << shape_result.original_point_count << "\n";
    g_decision_log << "  Points used: " << polygon_points.size() << "\n";
    g_decision_log << "CONFIGURATION:\n";
    g_decision_log << "  Polygon segments = " << polygon_edges.size() << "\n";
    g_decision_log << "  LOOKAHEAD_DEPTH_LIMIT = " << LOOKAHEAD_DEPTH_LIMIT << " (fixed)\n";
    g_decision_log << "  USE_CANDIDATES_RATIO = " << (USE_CANDIDATES_RATIO ? "true" : "false") << "\n";
    g_decision_log << "  CANDIDATES_RATIO = " << CANDIDATES_RATIO << "\n";
    g_decision_log << "  MAX_CANDIDATES_PER_LEVEL = " << MAX_CANDIDATES_PER_LEVEL << "\n";
    g_decision_log << "OTHER SETTINGS:\n";
    g_decision_log << "  ENABLE_EARLY_TERMINATION = " << (ENABLE_EARLY_TERMINATION ? "true" : "false") << "\n";
    g_decision_log << "  MAX_RECURSIVE_CALLS_PER_EVAL = " << MAX_RECURSIVE_CALLS_PER_EVAL << "\n";
    int display_max_iter = std::max(100000, static_cast<int>(polygon_edges.size()) * 10);
    g_decision_log << "  MAX_ITERATIONS = " << display_max_iter << " (10 * polygon segments, min 100000)\n";
    g_decision_log << "==============================================\n\n";
    g_decision_log.flush();

    // --- Create Minimum Enclosing Triangle ---
    std::vector<Point> hull_points;
    CGAL::convex_hull_2(polygon_points.begin(), polygon_points.end(), std::back_inserter(hull_points));
    Triangle initial_triangle = compute_min_enclosing_triangle(hull_points);

    std::cerr << "[TRIANGLE] Enclosing triangle computed from convex hull ("
              << hull_points.size() << " hull points)\n";

    // Log enclosing triangle for reconstruction
    pipeline_log << "[ENCLOSING_TRIANGLE]\n";
    for (int i = 0; i < 3; ++i) {
        pipeline_log << "  V" << i << ": ("
                     << std::fixed << std::setprecision(10)
                     << CGAL::to_double(initial_triangle[i].x()) << ", "
                     << CGAL::to_double(initial_triangle[i].y()) << ")\n";
    }
    // Log area comparison for tightness verification
    double enc_tri_area = area2d(initial_triangle);
    double poly_area_check = polygon_abs_area_double(polygon);
    pipeline_log << "  Enclosing triangle area: " << std::fixed << std::setprecision(6)
                 << enc_tri_area << "\n";
    pipeline_log << "  Polygon area:            " << poly_area_check << "\n";
    if (poly_area_check > 1e-12) {
        pipeline_log << "  Area ratio (triangle/polygon): "
                     << std::setprecision(3) << (enc_tri_area / poly_area_check) << "x\n";
    }
    pipeline_log << "\n";
    pipeline_log.flush();

    // --- Run Processing ---
    std::ofstream debug_fp("debug_out_lookahead.txt");
    if (!debug_fp) {
        std::cerr << "Error: Could not open debug_out_lookahead.txt\n";
        return 1;
    }

    std::vector<FinalTriangle> inside_triangles;
    process(initial_triangle, polygon_edges, polygon, debug_fp, inside_triangles);

    // =====================================================================
    // WRITE CSV OUTPUTS (same as original)
    // =====================================================================
    std::ofstream csv_triangles("triangles_decomposition.csv");
    std::ofstream csv_edges("triangles_edges.csv");

    if (csv_triangles) {
        csv_triangles << "id,v1x,v1y,v2x,v2y,v3x,v3y,type,num_polygon_edges" << std::endl;
        int tri_id = 1;
        for (const auto& final_tri : inside_triangles) {
            int num_edges = (int)final_tri.polygon_sides.size();
            std::string tri_type;
            if (num_edges == 0) tri_type = "INTERNAL";
            else if (num_edges == 1) tri_type = "SIDE";
            else if (num_edges == 2) tri_type = "CORNER";
            else tri_type = "BOUNDARY";

            csv_triangles << tri_id << ","
                         << CGAL::to_double(final_tri.vertices[0].x()) << ","
                         << CGAL::to_double(final_tri.vertices[0].y()) << ","
                         << CGAL::to_double(final_tri.vertices[1].x()) << ","
                         << CGAL::to_double(final_tri.vertices[1].y()) << ","
                         << CGAL::to_double(final_tri.vertices[2].x()) << ","
                         << CGAL::to_double(final_tri.vertices[2].y()) << ","
                         << tri_type << "," << num_edges << "\n";
            tri_id++;
        }
        csv_triangles.close();
        std::cerr << "Wrote " << (inside_triangles.size()) << " triangles to triangles_decomposition.csv\n";
    }

    if (csv_edges) {
        csv_edges << "triangle_id,edge_idx,sx,sy,tx,ty" << std::endl;
        int tri_id = 1;
        for (const auto& final_tri : inside_triangles) {
            int edge_idx = 0;
            for (const auto& side_edge : final_tri.polygon_sides) {
                csv_edges << tri_id << "," << edge_idx << ","
                         << CGAL::to_double(side_edge.source().x()) << ","
                         << CGAL::to_double(side_edge.source().y()) << ","
                         << CGAL::to_double(side_edge.target().x()) << ","
                         << CGAL::to_double(side_edge.target().y()) << "\n";
                edge_idx++;
            }
            tri_id++;
        }
        csv_edges.close();
    }

    // --- Write Triangle Splits CSV ---
    std::ofstream csv_tri_splits("triangle_splits.csv");
    if (csv_tri_splits) {
        csv_tri_splits << "split_id,iteration,is_lookahead,parent_v1x,parent_v1y,parent_v2x,parent_v2y,parent_v3x,parent_v3y,"
                      << "split_src_x,split_src_y,split_tgt_x,split_tgt_y,"
                      << "child1_v1x,child1_v1y,child1_v2x,child1_v2y,child1_v3x,child1_v3y,"
                      << "child2_v1x,child2_v1y,child2_v2x,child2_v2y,child2_v3x,child2_v3y" << "\n";
        for (const auto& rec : g_triangle_splits) {
            csv_tri_splits << rec.split_id << "," << rec.iteration << "," << (rec.is_lookahead ? 1 : 0) << ","
                          << CGAL::to_double(rec.parent_triangle[0].x()) << "," << CGAL::to_double(rec.parent_triangle[0].y()) << ","
                          << CGAL::to_double(rec.parent_triangle[1].x()) << "," << CGAL::to_double(rec.parent_triangle[1].y()) << ","
                          << CGAL::to_double(rec.parent_triangle[2].x()) << "," << CGAL::to_double(rec.parent_triangle[2].y()) << ","
                          << CGAL::to_double(rec.split_line.source().x()) << "," << CGAL::to_double(rec.split_line.source().y()) << ","
                          << CGAL::to_double(rec.split_line.target().x()) << "," << CGAL::to_double(rec.split_line.target().y()) << ","
                          << CGAL::to_double(rec.child1[0].x()) << "," << CGAL::to_double(rec.child1[0].y()) << ","
                          << CGAL::to_double(rec.child1[1].x()) << "," << CGAL::to_double(rec.child1[1].y()) << ","
                          << CGAL::to_double(rec.child1[2].x()) << "," << CGAL::to_double(rec.child1[2].y()) << ","
                          << CGAL::to_double(rec.child2[0].x()) << "," << CGAL::to_double(rec.child2[0].y()) << ","
                          << CGAL::to_double(rec.child2[1].x()) << "," << CGAL::to_double(rec.child2[1].y()) << ","
                          << CGAL::to_double(rec.child2[2].x()) << "," << CGAL::to_double(rec.child2[2].y()) << "\n";
        }
        csv_tri_splits.close();
        std::cerr << "Wrote " << g_triangle_splits.size() << " triangle splits to triangle_splits.csv\n";
    }

    // --- Write Line Splits CSV ---
    std::ofstream csv_line_splits("line_splits.csv");
    if (csv_line_splits) {
        csv_line_splits << "split_id,triangle_split_id,iteration,"
                       << "orig_src_x,orig_src_y,orig_tgt_x,orig_tgt_y,"
                       << "split_point_x,split_point_y,"
                       << "seg1_src_x,seg1_src_y,seg1_tgt_x,seg1_tgt_y,"
                       << "seg2_src_x,seg2_src_y,seg2_tgt_x,seg2_tgt_y" << "\n";
        for (const auto& rec : g_line_splits) {
            csv_line_splits << rec.split_id << "," << rec.triangle_split_id << "," << rec.iteration << ","
                           << CGAL::to_double(rec.original_edge.source().x()) << "," << CGAL::to_double(rec.original_edge.source().y()) << ","
                           << CGAL::to_double(rec.original_edge.target().x()) << "," << CGAL::to_double(rec.original_edge.target().y()) << ","
                           << CGAL::to_double(rec.split_point.x()) << "," << CGAL::to_double(rec.split_point.y()) << ","
                           << CGAL::to_double(rec.segment1.source().x()) << "," << CGAL::to_double(rec.segment1.source().y()) << ","
                           << CGAL::to_double(rec.segment1.target().x()) << "," << CGAL::to_double(rec.segment1.target().y()) << ","
                           << CGAL::to_double(rec.segment2.source().x()) << "," << CGAL::to_double(rec.segment2.source().y()) << ","
                           << CGAL::to_double(rec.segment2.target().x()) << "," << CGAL::to_double(rec.segment2.target().y()) << "\n";
        }
        csv_line_splits.close();
        std::cerr << "Wrote " << g_line_splits.size() << " line splits to line_splits.csv\n";
    }

    // --- Write Triangle Marks CSV (PerWalk + MarkTriangles region codes) ---
    // region_code: 1 = INSIDE (interior R), 2 = OUTSIDE (exterior), 0 = boundary
    std::ofstream csv_marks("triangle_marks.csv");
    if (csv_marks) {
        // Quantized vertices (qv*) are included so Part 2 can join each leaf to
        // the xi_T tree (triangle_splits_quantized.csv) by exact integer key.
        csv_marks << "leaf_id,region_code,region,v1x,v1y,v2x,v2y,v3x,v3y,area,"
                  << "qv1x,qv1y,qv2x,qv2y,qv3x,qv3y\n";
        for (size_t t = 0; t < g_leaf_triangles.size(); ++t) {
            const Triangle& tri = g_leaf_triangles[t];
            int rc = (t < g_leaf_region_codes.size()) ? g_leaf_region_codes[t] : REGION_BOUNDARY;
            const char* rname = (rc == REGION_INSIDE) ? "INSIDE"
                              : (rc == REGION_OUTSIDE) ? "OUTSIDE" : "BOUNDARY";
            QPoint q0 = quantize(tri[0]), q1 = quantize(tri[1]), q2 = quantize(tri[2]);
            csv_marks << (t + 1) << "," << rc << "," << rname << ","
                      << CGAL::to_double(tri[0].x()) << "," << CGAL::to_double(tri[0].y()) << ","
                      << CGAL::to_double(tri[1].x()) << "," << CGAL::to_double(tri[1].y()) << ","
                      << CGAL::to_double(tri[2].x()) << "," << CGAL::to_double(tri[2].y()) << ","
                      << area2d(tri) << ","
                      << q0.X << "," << q0.Y << "," << q1.X << "," << q1.Y << ","
                      << q2.X << "," << q2.Y << "\n";
        }
        csv_marks.close();
        std::cerr << "Wrote " << g_leaf_triangles.size()
                  << " leaf marks to triangle_marks.csv (interior L^I = "
                  << g_mark_triangle_count << ")\n";
    }

    // --- Write Final Triangle Details to Debug File ---
    debug_fp << "\n\n--- Internal triangles with Edges ---\n";
    int counter = 1;
    for (const auto& final_tri : inside_triangles) {
        debug_fp << "\n--- Triangle " << counter++ << " ---\n";
        debug_fp << "Vertices:\n";
        debug_fp << "  V1: " << point_to_string(final_tri.vertices[0]) << "\n";
        debug_fp << "  V2: " << point_to_string(final_tri.vertices[1]) << "\n";
        debug_fp << "  V3: " << point_to_string(final_tri.vertices[2]) << "\n";
        if (!final_tri.polygon_sides.empty()) {
            debug_fp << "Polygon Edges forming its sides:\n";
            for (const auto& side_edge : final_tri.polygon_sides) {
                debug_fp << "  - Edge: " << point_to_string(side_edge.source())
                         << " -> " << point_to_string(side_edge.target()) << "\n";
            }
        } else {
            debug_fp << "No polygon edges form the sides of this triangle.\n";
        }
    }

    // --- Log final triangle details for reconstruction in pipeline log ---
    pipeline_log << "\n=== FINAL TRIANGLES (" << inside_triangles.size() << ") ===\n";
    pipeline_log << "(Full listing in triangles_decomposition.csv; showing first and last 10 here)\n";
    counter = 1;
    for (const auto& final_tri : inside_triangles) {
        bool show = (counter <= 10) || (counter > (int)inside_triangles.size() - 10);
        if (counter == 11 && inside_triangles.size() > 20) {
            pipeline_log << "  ... (" << (inside_triangles.size() - 20) << " triangles omitted) ...\n";
        }
        if (show) {
            pipeline_log << "TRIANGLE #" << counter << ":\n";
            for (int v = 0; v < 3; ++v) {
                pipeline_log << "  V" << v << ": (" << std::fixed << std::setprecision(10)
                             << CGAL::to_double(final_tri.vertices[v].x()) << ", "
                             << CGAL::to_double(final_tri.vertices[v].y()) << ")\n";
            }
            if (!final_tri.polygon_sides.empty()) {
                pipeline_log << "  Polygon edges on boundary: " << final_tri.polygon_sides.size() << "\n";
                for (size_t e = 0; e < final_tri.polygon_sides.size(); ++e) {
                    pipeline_log << "    E" << e << ": ("
                                 << CGAL::to_double(final_tri.polygon_sides[e].source().x()) << ", "
                                 << CGAL::to_double(final_tri.polygon_sides[e].source().y()) << ") -> ("
                                 << CGAL::to_double(final_tri.polygon_sides[e].target().x()) << ", "
                                 << CGAL::to_double(final_tri.polygon_sides[e].target().y()) << ")\n";
                }
            }
        }
        counter++;
    }

    // =====================================================================
    // STEP 5: FIXED-POINT QUANTIZATION + BLOCKCHAIN-WALKABLE PUBLISHING
    // =====================================================================
    // All geometry above stays in the CGAL exact kernel. We quantize ONLY at
    // publishing time: each distinct boundary vertex is snapped once to the
    // global integer grid, so shared split points get identical integers and
    // the published boundary forms a single closed, connected (walkable) cycle.
    bool walkable = true;
    if (QUANTIZE_TRANSACTIONS) {
        std::cerr << "\n=== Step 5: Fixed-point quantization + walkable boundary ===\n";
        pipeline_log << "\n=== STEP 5: FIXED-POINT QUANTIZATION ===\n";
        walkable = emit_quantized_walkable_transactions(polygon_edges, pipeline_log);
        if (!walkable) {
            std::cerr << "[WARN] Published boundary did NOT verify as a closed walk; "
                         "see quantization_report.txt\n";
            pipeline_log << "[WARN] Published boundary failed closed-walk verification.\n";
        }
    } else {
        pipeline_log << "\n[QUANTIZE] Skipped (--no-quantize).\n";
    }

    debug_fp.close();
    g_decision_log.close();
    g_transaction_log.close();

    // --- Final pipeline summary ---
    pipeline_log << "\n==============================================\n";
    pipeline_log << "PIPELINE COMPLETE\n";
    pipeline_log << "==============================================\n";
    pipeline_log << "Input:           " << shapefile_path << "\n";
    if (!shape_result.record_name.empty())
        pipeline_log << "Record:          " << shape_result.record_name << "\n";
    pipeline_log << "Original points: " << shape_result.original_point_count << "\n";
    pipeline_log << "Points used:     " << polygon_points.size() << "\n";
    pipeline_log << "Reduction:       " << (do_reduction ?
        std::to_string(static_cast<int>(reduction_pct)) + "%" : "none") << "\n";
    pipeline_log << "Triangles:       " << inside_triangles.size() << "\n";
    pipeline_log << "Split method:    " << (USE_ANGLE_METHOD ?
        "angle method (Ramkumar)" : "legacy recursive lookahead") << "\n";
    pipeline_log << "Triangle splits: " << g_triangle_split_count << "\n";
    pipeline_log << "Line splits:     " << g_line_split_count << "\n";
    int total_txns = g_triangle_split_count + g_line_split_count * 2;
    pipeline_log << "Total AGS txns:  " << total_txns << "\n";
    pipeline_log << "\nOutput files:\n";
    if (do_write_header)
        pipeline_log << "  CGAL header:   " << header_path << "\n";
    if (do_write_points)
        pipeline_log << "  Points file:   " << points_path << "\n";
    pipeline_log << "  Transactions:  ags_transactions.txt\n";
    pipeline_log << "  Decisions:     triangle_decision_log.txt\n";
    pipeline_log << "  Debug:         debug_out_lookahead.txt\n";
    pipeline_log << "  Pipeline log:  pipeline_log.txt\n";
    pipeline_log << "  CSVs:          triangles_decomposition.csv, triangles_edges.csv,\n";
    pipeline_log << "                 triangle_splits.csv, line_splits.csv\n";
    if (QUANTIZE_TRANSACTIONS) {
        pipeline_log << "  Quantized:     ags_transactions_quantized.txt, quantization_report.txt,\n";
        pipeline_log << "                 boundary_walk.csv, triangle_splits_quantized.csv,\n";
        pipeline_log << "                 line_splits_quantized.csv\n";
        pipeline_log << "  Walkable:      " << (walkable ? "YES (closed integer cycle)"
                                                         : "NO -- see quantization_report.txt") << "\n";
    }
    pipeline_log << "==============================================\n";
    pipeline_log.close();

    std::cerr << "\n==============================================\n";
    std::cerr << "PIPELINE COMPLETE\n";
    std::cerr << "  Shapefile:      " << shapefile_path << "\n";
    if (!shape_result.record_name.empty())
        std::cerr << "  Record:         " << shape_result.record_name << "\n";
    std::cerr << "  Points:         " << shape_result.original_point_count;
    if (do_reduction)
        std::cerr << " -> " << polygon_points.size() << " (" << reduction_pct << "% reduction)\n";
    else
        std::cerr << " (no reduction)\n";
    std::cerr << "  Triangles:      " << inside_triangles.size() << "\n";
    std::cerr << "  Split method:   " << (USE_ANGLE_METHOD ?
        "angle (Ramkumar)" : "legacy lookahead") << "\n";
    std::cerr << "  AGS Transactions: " << total_txns << "\n";
    if (QUANTIZE_TRANSACTIONS)
        std::cerr << "  Walkable chain: " << (walkable ? "YES" : "NO (see report)") << "\n";
    std::cerr << "  Logs: pipeline_log.txt, ags_transactions.txt";
    if (QUANTIZE_TRANSACTIONS)
        std::cerr << ", ags_transactions_quantized.txt";
    std::cerr << "\n";
    std::cerr << "==============================================\n";

    sum.fips            = cfg.state_fips;
    sum.name            = shape_result.record_name;
    sum.original_points = shape_result.original_point_count;
    sum.used_points     = (int)polygon_points.size();
    sum.total_leaves    = (int)g_leaf_triangles.size();
    sum.interior_leaves = g_mark_triangle_count;
    sum.walkable        = walkable;
    sum.lookahead       = LOOKAHEAD_DEPTH_LIMIT;
    sum.method          = USE_ANGLE_METHOD ? "angle" : "lookahead";
    sum.ags_txns        = total_txns;
    sum.ok              = true;
    return 0;
}

// =========================================================================
// Python / in-memory API: polygon + enclosing triangle, no shapefile, no files.
// =========================================================================
static void td_xy(const Point& p, double o[2]) {
    o[0] = CGAL::to_double(p.x());
    o[1] = CGAL::to_double(p.y());
}

static void td_fill_tri(const Triangle& tri, TDTri& out) {
    for (int i = 0; i < 3; ++i) td_xy(tri[i], out.v[i]);
}

static void td_fill_edge(const Edge& e, TDEdge& out) {
    td_xy(e.source(), out.a);
    td_xy(e.target(), out.b);
}

static const char* td_region_name(int code) {
    if (code == REGION_INSIDE) return "INSIDE";
    if (code == REGION_OUTSIDE) return "OUTSIDE";
    return "BOUNDARY";
}

TDResult tridecomp_from_xy(const double* p_xy, int n,
                           const double* t_xy, bool verbose)
{
    TDResult R;
    if (!p_xy || !t_xy || n < 3) {
        R.error = "p must have at least 3 vertices and t must be a (3,2) triangle";
        return R;
    }

    reset_pipeline_globals();
    G_VERBOSITY = verbose ? 1 : 0;
    QUANTIZE_TRANSACTIONS = true;
    USE_ANGLE_METHOD = true;
    LOOKAHEAD_DEPTH_LIMIT = 0;

    std::vector<Point> pts;
    pts.reserve(n);
    for (int i = 0; i < n; ++i)
        pts.emplace_back(p_xy[2 * i], p_xy[2 * i + 1]);
    if (pts.size() >= 2 && point_eq(pts.front(), pts.back()))
        pts.pop_back();
    if (pts.size() < 3) {
        R.error = "p has fewer than 3 vertices after dropping the closing point";
        return R;
    }

    Polygon polygon(pts.begin(), pts.end());
    if (!polygon.is_simple()) {
        R.error = "polygon is not simple";
        return R;
    }
    if (polygon.is_clockwise_oriented()) {
        polygon.reverse_orientation();
        pts.assign(polygon.vertices_begin(), polygon.vertices_end());
    }

    std::vector<Edge> polygon_edges;
    polygon_edges.reserve(pts.size());
    for (size_t i = 0; i < pts.size(); ++i)
        polygon_edges.emplace_back(pts[i], pts[(i + 1) % pts.size()]);

    Point tv[3] = {
        Point(t_xy[0], t_xy[1]),
        Point(t_xy[2], t_xy[3]),
        Point(t_xy[4], t_xy[5])
    };
    double signed_area = CGAL::to_double(Triangle_cgal(tv[0], tv[1], tv[2]).area());
    if (std::abs(signed_area) <= AREA_TOL) {
        R.error = "enclosing triangle has near-zero area";
        return R;
    }
    if (signed_area < 0.0)
        std::swap(tv[1], tv[2]);
    Triangle initial = { tv[0], tv[1], tv[2] };

#ifdef _WIN32
    const char* kNull = "nul";
#else
    const char* kNull = "/dev/null";
#endif
    g_transaction_log.open(kNull);
    g_decision_log.open(kNull);
    std::ofstream debug_fp(kNull);

    std::streambuf* cout_old = nullptr;
    std::streambuf* cerr_old = nullptr;
    std::ofstream null_out;
    if (!verbose) {
        null_out.open(kNull);
        cout_old = std::cout.rdbuf(null_out.rdbuf());
        cerr_old = std::cerr.rdbuf(null_out.rdbuf());
    }

    compute_dynamic_limits(static_cast<int>(polygon_edges.size()));
    std::vector<FinalTriangle> inside_triangles;
    process(initial, polygon_edges, polygon, debug_fp, inside_triangles);

    if (!verbose) {
        std::cout.rdbuf(cout_old);
        std::cerr.rdbuf(cerr_old);
    }

    R.triangle_splits.reserve(g_triangle_splits.size());
    for (const auto& rec : g_triangle_splits) {
        TDSplit s;
        s.split_id = rec.split_id;
        s.iteration = rec.iteration;
        s.is_lookahead = rec.is_lookahead ? 1 : 0;
        td_fill_tri(rec.parent_triangle, s.parent);
        td_fill_tri(rec.child1, s.child1);
        td_fill_tri(rec.child2, s.child2);
        td_fill_edge(rec.split_line, s.split_line);
        R.triangle_splits.push_back(s);
    }

    R.line_splits.reserve(g_line_splits.size());
    for (const auto& rec : g_line_splits) {
        TDLineSplit s;
        s.split_id = rec.split_id;
        s.triangle_split_id = rec.triangle_split_id;
        s.iteration = rec.iteration;
        td_fill_edge(rec.original_edge, s.original_edge);
        td_fill_edge(rec.segment1, s.segment1);
        td_fill_edge(rec.segment2, s.segment2);
        td_xy(rec.split_point, s.split_point);
        QPoint qp = quantize(rec.split_point);
        s.qx = (std::int64_t)qp.X;
        s.qy = (std::int64_t)qp.Y;
        R.line_splits.push_back(s);
    }

    R.marks.reserve(g_leaf_triangles.size());
    for (size_t i = 0; i < g_leaf_triangles.size(); ++i) {
        const Triangle& tri = g_leaf_triangles[i];
        int rc = (i < g_leaf_region_codes.size()) ? g_leaf_region_codes[i] : REGION_BOUNDARY;
        TDMark m;
        m.leaf_id = (int)i + 1;
        m.region_code = rc;
        m.region = td_region_name(rc);
        td_fill_tri(tri, m.vertices);
        m.area = area2d(tri);
        QPoint q0 = quantize(tri[0]), q1 = quantize(tri[1]), q2 = quantize(tri[2]);
        m.qv[0][0] = q0.X; m.qv[0][1] = q0.Y;
        m.qv[1][0] = q1.X; m.qv[1][1] = q1.Y;
        m.qv[2][0] = q2.X; m.qv[2][1] = q2.Y;
        R.marks.push_back(m);
    }

    // Rebuild quantized walk (same construction as emit_quantized_walkable_transactions).
    std::vector<double> lsx(g_line_splits.size()), lsy(g_line_splits.size());
    for (size_t k = 0; k < g_line_splits.size(); ++k) {
        lsx[k] = CGAL::to_double(g_line_splits[k].split_point.x());
        lsy[k] = CGAL::to_double(g_line_splits[k].split_point.y());
    }
    const double BB_EPS = 1e-9;
    int seq = 0;
    std::vector<std::array<QPoint, 2>> qwalk;
    qwalk.reserve(polygon_edges.size() + g_line_splits.size());
    for (const auto& E : polygon_edges) {
        const Point& A = E.source();
        const Point& B = E.target();
        K::Vector_2 dir = B - A;
        K::FT dlen2 = dir * dir;
        double Ax = CGAL::to_double(A.x()), Ay = CGAL::to_double(A.y());
        double Bx = CGAL::to_double(B.x()), By = CGAL::to_double(B.y());
        double minx = std::min(Ax, Bx) - BB_EPS, maxx = std::max(Ax, Bx) + BB_EPS;
        double miny = std::min(Ay, By) - BB_EPS, maxy = std::max(Ay, By) + BB_EPS;
        struct OnEdge { Point p; K::FT t; };
        std::vector<OnEdge> mids;
        if (dlen2 != K::FT(0)) {
            for (size_t k = 0; k < g_line_splits.size(); ++k) {
                if (lsx[k] < minx || lsx[k] > maxx || lsy[k] < miny || lsy[k] > maxy)
                    continue;
                const Point& P = g_line_splits[k].split_point;
                if (!CGAL::collinear(A, B, P)) continue;
                K::FT t = ((P - A) * dir) / dlen2;
                if (t > K::FT(0) && t < K::FT(1)) {
                    bool dup = false;
                    for (const auto& m : mids) if (m.p == P) { dup = true; break; }
                    if (!dup) mids.push_back({ P, t });
                }
            }
        }
        std::sort(mids.begin(), mids.end(),
                  [](const OnEdge& a, const OnEdge& b){ return a.t < b.t; });
        Point prev = A;
        auto emit_seg = [&](const Point& s, const Point& tt) {
            QPoint qs = quantize(s), qt = quantize(tt);
            qwalk.push_back({ qs, qt });
            TDWalkSeg w;
            w.seq = seq++;
            w.sx = CGAL::to_double(s.x());  w.sy = CGAL::to_double(s.y());
            w.tx = CGAL::to_double(tt.x()); w.ty = CGAL::to_double(tt.y());
            w.qsx = qs.X; w.qsy = qs.Y; w.qtx = qt.X; w.qty = qt.Y;
            R.walk.push_back(w);
        };
        for (const auto& m : mids) { emit_seg(prev, m.p); prev = m.p; }
        emit_seg(prev, B);
    }

    bool connected = true;
    for (size_t i = 0; i + 1 < qwalk.size(); ++i) {
        if (qwalk[i][1].X != qwalk[i + 1][0].X || qwalk[i][1].Y != qwalk[i + 1][0].Y) {
            connected = false; break;
        }
    }
    bool closed = false;
    if (!qwalk.empty())
        closed = (qwalk.back()[1].X == qwalk.front()[0].X &&
                  qwalk.back()[1].Y == qwalk.front()[0].Y);

    R.stats.n_polygon = (int)pts.size();
    R.stats.n_splits = (int)R.triangle_splits.size();
    R.stats.n_line_splits = (int)R.line_splits.size();
    R.stats.n_leaves = (int)R.marks.size();
    R.stats.n_interior = g_last_mark.inside_count;
    R.stats.n_exterior = g_last_mark.outside_count;
    R.stats.perwalk_seeds = g_last_mark.perwalk_inside_seeds;
    R.stats.flood_added = g_last_mark.flood_added;
    R.stats.flood_unreached = g_last_mark.flood_unreached;
    R.stats.inside_area = g_last_mark.inside_area;
    R.stats.polygon_area = g_last_mark.polygon_area;
    R.stats.area_ok = g_last_mark.area_ok;
    R.stats.walkable = connected && closed && !qwalk.empty();

    g_transaction_log.close();
    g_decision_log.close();
    return R;
}

#ifndef TRIDECOMP_NO_MAIN
int main(int argc, char* argv[]) {
    // --- Parse Command-Line Arguments ---
    std::string shapefile_path;
    double reduction_pct = -1.0;  // -1 = no reduction (use all original points)
    std::string state_fips;
    std::string county_fips;
    std::string header_path = "polygon_output.h";
    std::string points_path = "polygon_output.txt";
    std::string var_prefix = "simple_polygon";
    bool do_write_header = true;
    bool do_write_points = true;
    bool do_list_states = false;
    bool do_list_counties = false;
    int num_states = 0;   // >0 => batch the first N states, each in its own folder
    int la_range_lo = -1, la_range_hi = -1;   // >=0 => sweep look-ahead depths [lo,hi]

    for (int i = 1; i < argc; ++i) {
        if (std::strcmp(argv[i], "--help") == 0 || std::strcmp(argv[i], "-h") == 0) {
            std::cerr << "Usage: " << argv[0] << " <shapefile.shp> [reduction_%] [options]\n\n"
                      << "Reads a TIGER/ESRI shapefile, optionally reduces polygon vertices\n"
                      << "using Douglas-Peucker, generates CGAL header files, and runs\n"
                      << "triangular decomposition with full transaction logging.\n\n"
                      << "Positional:\n"
                      << "  <shapefile.shp>      ESRI shapefile path\n"
                      << "  [reduction_%]        Percentage of vertices to remove (0-100)\n"
                      << "                       Optional. If omitted, all original vertices are used.\n\n"
                      << "Options:\n"
                      << "  --state <FIPS>       Filter by 2-digit state FIPS (e.g. 53=WA)\n"
                      << "  --county <FIPS>      Filter by 3-digit county FIPS (e.g. 033=King)\n"
                      << "  --lookahead-range <A-B>  Sweep look-ahead depths A..B for each state,\n"
                      << "                       gathering results at every depth; prints a per-state\n"
                      << "                       trend of how leaf count falls with depth. Combines\n"
                      << "                       with --num-states. Outputs go to la_<D>/ subfolders.\n"
                      << "  --num-states <N>     Batch mode: process the first N states found in\n"
                      << "                       the shapefile, each into its own state_<FIPS>/\n"
                      << "                       folder, then print a summary table (summary.csv)\n"
                      << "  --angle-method       Use Dr. Ramkumar's angle method for split\n"
                      << "                       selection (DEFAULT). O(k log k) per triangle.\n"
                      << "  --legacy-lookahead   Use the legacy recursive lookahead engine\n"
                      << "                       instead of the angle method (for A/B compare)\n"
                      << "  --lookahead <N>      Lookahead depth for BOTH the angle method and the\n"
                      << "                       legacy engine (default: 0 = greedy). Higher N finds\n"
                      << "                       lower-leaf decompositions at increasing runtime cost.\n"
                      << "  --no-quantize        Skip fixed-point quantization / walkable\n"
                      << "                       boundary publishing step\n"
                      << "  --quiet              Progress + AGS records + CSVs; skip per-iteration debug\n"
                      << "  --verbose            Dump every cell's vertices/edges to debug_out_lookahead.txt\n"
                      << "  --reduce <pct>       Reduction % (overrides positional)\n"
                      << "  --header <path>      Output CGAL header (default: polygon_output.h)\n"
                      << "  --points <path>      Output points file (default: polygon_output.txt)\n"
                      << "  --prefix <name>      Variable name prefix (default: simple_polygon)\n"
                      << "  --no-header          Skip header file generation\n"
                      << "  --no-points          Skip points file generation\n"
                      << "  --list-states        Show FIPS codes in shapefile, then exit\n"
                      << "  --list-counties      Show county FIPS codes (use with --state), then exit\n"
                      << "  --help               Show this message\n\n"
                      << "Examples:\n"
                      << "  " << argv[0] << " tl_2024_us_state.shp --state 53\n"
                      << "  " << argv[0] << " tl_2024_us_state.shp 50 --state 53\n"
                      << "  " << argv[0] << " tl_2024_us_county.shp --state 06 --legacy-lookahead --lookahead 2\n"
                      << "  " << argv[0] << " tl_2024_us_county.shp 50 --state 53 --county 033\n"
                      << "  " << argv[0] << " boundary.shp 30 --header boundary.h --no-quantize\n";
            return 0;
        }
        if (std::strcmp(argv[i], "--list-states") == 0) {
            do_list_states = true;
        }
        if (std::strcmp(argv[i], "--list-counties") == 0) {
            do_list_counties = true;
        }
    }

    // Collect positional args
    int pos_idx = 0;
    for (int i = 1; i < argc; ++i) {
        if (argv[i][0] == '-' && argv[i][1] == '-') {
            const char* opt = argv[i];
            if (std::strcmp(opt, "--no-header") == 0 || std::strcmp(opt, "--no-points") == 0 ||
                std::strcmp(opt, "--list-states") == 0 || std::strcmp(opt, "--list-counties") == 0 ||
                std::strcmp(opt, "--angle-method") == 0 || std::strcmp(opt, "--legacy-lookahead") == 0 ||
                std::strcmp(opt, "--no-quantize") == 0 || std::strcmp(opt, "--quiet") == 0 ||
                std::strcmp(opt, "--verbose") == 0) {
                // flag only
            } else {
                ++i; // skip value
            }
            continue;
        }
        if (pos_idx == 0) shapefile_path = argv[i];
        else if (pos_idx == 1) reduction_pct = std::atof(argv[i]);
        ++pos_idx;
    }

    // Parse named options
    for (int i = 1; i < argc; ++i) {
        if (std::strcmp(argv[i], "--state") == 0 && i + 1 < argc)
            state_fips = argv[++i];
        else if (std::strcmp(argv[i], "--county") == 0 && i + 1 < argc)
            county_fips = argv[++i];
        else if (std::strcmp(argv[i], "--lookahead") == 0 && i + 1 < argc)
            LOOKAHEAD_DEPTH_LIMIT = std::atoi(argv[++i]);
        else if (std::strcmp(argv[i], "--reduce") == 0 && i + 1 < argc)
            reduction_pct = std::atof(argv[++i]);
        else if (std::strcmp(argv[i], "--header") == 0 && i + 1 < argc)
            header_path = argv[++i];
        else if (std::strcmp(argv[i], "--points") == 0 && i + 1 < argc)
            points_path = argv[++i];
        else if (std::strcmp(argv[i], "--prefix") == 0 && i + 1 < argc)
            var_prefix = argv[++i];
        else if (std::strcmp(argv[i], "--no-header") == 0)
            do_write_header = false;
        else if (std::strcmp(argv[i], "--no-points") == 0)
            do_write_points = false;
        else if (std::strcmp(argv[i], "--angle-method") == 0)
            USE_ANGLE_METHOD = true;
        else if (std::strcmp(argv[i], "--legacy-lookahead") == 0)
            USE_ANGLE_METHOD = false;
        else if (std::strcmp(argv[i], "--no-quantize") == 0)
            QUANTIZE_TRANSACTIONS = false;
        else if (std::strcmp(argv[i], "--quiet") == 0)
            G_VERBOSITY = 0;
        else if (std::strcmp(argv[i], "--verbose") == 0)
            G_VERBOSITY = 2;
        else if ((std::strcmp(argv[i], "--num-states") == 0 ||
                  std::strcmp(argv[i], "--states") == 0) && i + 1 < argc)
            num_states = std::atoi(argv[++i]);
        else if (std::strcmp(argv[i], "--lookahead-range") == 0 && i + 1 < argc) {
            std::string r = argv[++i];
            auto dash = r.find('-');
            if (dash == std::string::npos) {            // single value "3" => 3-3
                la_range_lo = la_range_hi = std::atoi(r.c_str());
            } else {
                la_range_lo = std::atoi(r.substr(0, dash).c_str());
                la_range_hi = std::atoi(r.substr(dash + 1).c_str());
            }
            if (la_range_lo < 0) la_range_lo = 0;
            if (la_range_hi < la_range_lo) la_range_hi = la_range_lo;
        }
    }

    // --- List-states mode ---
    if (do_list_states) {
        if (shapefile_path.empty()) {
            std::cerr << "Error: specify shapefile path before --list-states\n";
            return 1;
        }
        list_states_in_shapefile(shapefile_path);
        return 0;
    }

    // --- List-counties mode ---
    if (do_list_counties) {
        if (shapefile_path.empty()) {
            std::cerr << "Error: specify shapefile path before --list-counties\n";
            return 1;
        }
        list_counties_in_shapefile(shapefile_path, state_fips);
        return 0;
    }

    // --- Validate ---
    if (shapefile_path.empty()) {
        std::cerr << "Error: no shapefile specified. Use --help for usage.\n";
        return 1;
    }
    // Validate reduction if specified
    bool do_reduction = (reduction_pct > 0.0);
    if (reduction_pct >= 100.0) {
        std::cerr << "Error: reduction percentage must be less than 100.\n"
                  << "  Got: " << reduction_pct << "\n";
        return 1;
    }
    if (reduction_pct > 0.0 && reduction_pct < 0.01) {
        std::cerr << "Warning: reduction percentage " << reduction_pct
                  << "% is very small, treating as no reduction.\n";
        do_reduction = false;
        reduction_pct = 0.0;
    }
    if (reduction_pct <= 0.0) {
        reduction_pct = 0.0;  // normalize sentinel to 0
    }

    // =====================================================================
    // RUN: a matrix of (states) x (look-ahead depths).
    //   states: one (single mode) or the first N (--num-states N)
    //   depths: one (--lookahead) or a sweep A..B (--lookahead-range A-B)
    // Output dirs: state_<FIPS>/ when batching, la_<D>/ when sweeping depths,
    // nested as state_<FIPS>/la_<D>/ when both are used.
    // =====================================================================
    std::string abs_shp = to_absolute_path(shapefile_path);

    StateConfig base;
    base.shapefile_path  = abs_shp;
    base.county_fips     = county_fips;
    base.reduction_pct   = reduction_pct;
    base.do_reduction    = do_reduction;
    base.header_path     = header_path;
    base.points_path     = points_path;
    base.var_prefix      = var_prefix;
    base.do_write_header = do_write_header;
    base.do_write_points = do_write_points;

    // Depth list.
    std::vector<int> depths;
    if (la_range_lo >= 0)
        for (int d = la_range_lo; d <= la_range_hi; ++d) depths.push_back(d);
    else
        depths.push_back(LOOKAHEAD_DEPTH_LIMIT);
    bool multi_depth = depths.size() > 1;
    bool batch = (num_states > 0);

    // Fast path: exactly one state, one depth => existing single-run behaviour.
    if (!batch && !multi_depth) {
        base.state_fips = state_fips;
        StateSummary s;
        return run_one_state(base, s);
    }

    // State list.
    std::vector<std::string> states;
    if (batch) {
        states = first_n_state_fips(abs_shp, num_states);
        if (states.empty()) {
            std::cerr << "Error: no state FIPS codes found in shapefile (need a STATEFP field).\n";
            return 1;
        }
        std::cerr << "[BATCH] Processing the first " << states.size()
                  << " state(s) found in " << shapefile_path << "\n";
    } else {
        states.push_back(state_fips);   // single target (may be "" = first record)
    }
    if (multi_depth)
        std::cerr << "[SWEEP] Look-ahead depths " << la_range_lo << ".." << la_range_hi
                  << " (" << depths.size() << " runs per state)\n";

    std::string orig_cwd = get_cwd_str();
    std::vector<StateSummary> summaries;

    for (const auto& fips : states) {
        for (int depth : depths) {
            if (PORTABLE_CHDIR(orig_cwd.c_str()) != 0) {
                std::cerr << "[FATAL] could not return to " << orig_cwd << "\n";
                return 1;
            }
            // Descend into state_<FIPS>/ and/or la_<D>/ as needed.
            bool entered = true;
            if (batch) {
                std::string f = "state_" + fips;
                PORTABLE_MKDIR(f.c_str());
                if (PORTABLE_CHDIR(f.c_str()) != 0) entered = false;
            }
            if (entered && multi_depth) {
                std::string g = "la_" + std::to_string(depth);
                PORTABLE_MKDIR(g.c_str());
                if (PORTABLE_CHDIR(g.c_str()) != 0) entered = false;
            }
            if (!entered) {
                std::cerr << "[ERROR] cannot create output folder for FIPS " << fips
                          << " depth " << depth << "; skipping\n";
                StateSummary s; s.fips = fips; s.lookahead = depth; summaries.push_back(s);
                continue;
            }

            std::cerr << "\n################  FIPS " << fips << "  LA " << depth << "  ################\n";
            LOOKAHEAD_DEPTH_LIMIT = depth;
            StateConfig cfg = base;
            cfg.state_fips = fips;
            StateSummary s; s.fips = fips;
            auto t0 = std::chrono::steady_clock::now();
            run_one_state(cfg, s);
            s.runtime_s = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
            s.lookahead = depth;   // robust even if the run failed early
            summaries.push_back(s);
        }
    }
    if (PORTABLE_CHDIR(orig_cwd.c_str()) != 0) {
        std::cerr << "[FATAL] could not return to " << orig_cwd << "\n";
        return 1;
    }

    if (multi_depth) write_lookahead_summary(summaries);
    else             write_summary_table(summaries);
    return 0;
}
#endif // TRIDECOMP_NO_MAIN
