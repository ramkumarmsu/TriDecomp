#pragma once
// In-memory tessellation: constraint segments + enclosing triangle.
// Geometry is CGAL EPEC; coordinates in the result are doubles.

#include <string>
#include <vector>

struct TDTri {
    double v[3][2];
};

struct TDEdge {
    double a[2], b[2];
};

struct TDLineSplit {
    TDEdge original_edge{};
    double split_point[2]{};
    TDEdge segment1{}, segment2{};
};

// id is the bitstring of the parent (start code of t, or that plus child bits).
// child0 = id + "0" = left of the directed cut; child1 = id + "1" = right.
struct TDSplitEvent {
    std::string id;
    std::string child0;
    std::string child1;
    TDTri parent{};
    TDEdge cut{};
    TDTri child[2]{};
    std::vector<TDLineSplit> line_splits;
};

struct TDTessellation {
    std::vector<TDSplitEvent> events;
    std::string error;
};

// seg_xy: nseg segments, each 4 doubles (x0,y0,x1,y1). t_xy: 3 points, row-major.
// start: tree address of t (empty = fresh enclosing triangle).
TDTessellation tessellate_from_xy(const double* seg_xy, int nseg,
                                  const double* t_xy, const char* start,
                                  bool verbose);
