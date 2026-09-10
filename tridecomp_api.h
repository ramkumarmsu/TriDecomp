#pragma once
// In-memory triangular decomposition (same contract as Python tridecomp(p, t)).
// Geometry is CGAL EPEC; coordinates in the result are doubles / integers.

#include <cstdint>
#include <string>
#include <vector>

struct TDTri {
    double v[3][2];
};

struct TDEdge {
    double a[2], b[2];
};

struct TDSplit {
    int split_id = 0;
    int iteration = 0;
    int priority = 0;       // unused in C++ records; kept for API parity
    int is_lookahead = 0;
    TDTri parent{}, child1{}, child2{};
    TDEdge split_line{};
};

struct TDLineSplit {
    int split_id = 0;
    int triangle_split_id = 0;
    int iteration = 0;
    TDEdge original_edge{}, segment1{}, segment2{};
    double split_point[2]{};
    std::int64_t qx = 0, qy = 0;
};

struct TDMark {
    int leaf_id = 0;
    int region_code = 0;
    std::string region;
    TDTri vertices{};
    double area = 0.0;
    std::int64_t qv[3][2]{};
};

struct TDWalkSeg {
    int seq = 0;
    double sx = 0, sy = 0, tx = 0, ty = 0;
    std::int64_t qsx = 0, qsy = 0, qtx = 0, qty = 0;
};

struct TDStats {
    int n_polygon = 0;
    int n_splits = 0;
    int n_line_splits = 0;
    int n_leaves = 0;
    int n_interior = 0;
    int n_exterior = 0;
    int perwalk_seeds = 0;
    int flood_added = 0;
    int flood_unreached = 0;
    double inside_area = 0.0;
    double polygon_area = 0.0;
    bool area_ok = false;
    bool walkable = false;
};

struct TDResult {
    std::vector<TDSplit> triangle_splits;
    std::vector<TDLineSplit> line_splits;
    std::vector<TDMark> marks;
    std::vector<TDWalkSeg> walk;
    TDStats stats;
    std::string error;
};

// p_xy: n points, row-major (x,y). t_xy: 3 points, row-major.
// Closing duplicate on p is dropped. Clockwise rings are reversed.
TDResult tridecomp_from_xy(const double* p_xy, int n,
                           const double* t_xy, bool verbose);
