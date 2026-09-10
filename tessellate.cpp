// tessellate.cpp — CGAL EPEC engine for tessellate(segments, t, start).
// Constraint segments + enclosing triangle. No shapefile, marks, or CLI.

#include "tridecomp_api.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <iostream>
#include <limits>
#include <queue>
#include <string>
#include <utility>
#include <vector>

#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Object.h>
#include <CGAL/centroid.h>
#include <CGAL/intersections.h>
#include <CGAL/number_utils.h>

using K = CGAL::Exact_predicates_exact_constructions_kernel;
using Point = K::Point_2;
using Edge = K::Segment_2;
using Triangle_cgal = K::Triangle_2;
using Line = K::Line_2;
using Triangle = std::array<Point, 3>;

static const double POINT_TOL = 1e-5;
static const double AREA_TOL = 1e-6;
static const int MAX_PER_SIDE = 50;
static const double EMBED_BONUS = 1000.0;
static const double AEPS = 1e-9;

struct LineSplitRecord {
    Edge original_edge;
    Point split_point;
    Edge segment1;
    Edge segment2;
};

struct SplitEvent {
    std::string id;
    Triangle parent;
    Triangle child0;  // left of directed cut
    Triangle child1;  // right
    Edge cut;
    std::vector<LineSplitRecord> line_splits;
};

struct WorkItem {
    Triangle tri;
    std::vector<Edge> edges;
    std::string addr;
};

struct BestCutResult {
    bool found = false;
    std::array<Triangle, 2> child_triangles;
    Edge split_line;
    int priority = 0;
};

struct AngleEndpoint {
    double angle;
    Point pt;
};

struct ScoredCut {
    double cost;
    int corner;
    Point target;
    int cuts;
    int embedded;
};

static double area2d(const Triangle& triangle)
{
    auto exact_area = Triangle_cgal(triangle[0], triangle[1], triangle[2]).area();
    return CGAL::to_double(CGAL::abs(exact_area));
}

static bool point_eq(const Point& p1, const Point& p2, double tol = POINT_TOL)
{
    return CGAL::squared_distance(p1, p2) < tol * tol;
}

static bool edge_on_triangle_side(const Edge& edge, const Triangle& triangle)
{
    for (int i = 0; i < 3; ++i) {
        Edge side(triangle[i], triangle[(i + 1) % 3]);
        if (side.has_on(edge.source()) && side.has_on(edge.target()))
            return true;
    }
    return false;
}

static std::vector<Edge> filter_degenerate_edges(const std::vector<Edge>& edges)
{
    std::vector<Edge> result;
    for (const auto& e : edges)
        if (!point_eq(e.source(), e.target()))
            result.push_back(e);
    return result;
}

static bool split_line_is_polygon_edge(const Edge& split_line,
                                       const std::vector<Edge>& all_polygon_edges)
{
    for (const auto& pe : all_polygon_edges) {
        if (split_line.supporting_line() != pe.supporting_line())
            continue;
        if (pe.has_on(split_line.source()) && pe.has_on(split_line.target()))
            return true;
        if (split_line.has_on(pe.source()) && split_line.has_on(pe.target()))
            return true;
    }
    return false;
}

static bool is_polygon_vertex(const Point& p, const std::vector<Edge>& all_polygon_edges)
{
    for (const auto& pe : all_polygon_edges)
        if (point_eq(p, pe.source()) || point_eq(p, pe.target()))
            return true;
    return false;
}

static int count_intersections_improved(const Edge& split_line,
                                        const std::vector<Edge>& polygon_edges)
{
    for (const auto& poly_edge : polygon_edges) {
        if (split_line.supporting_line() == poly_edge.supporting_line()) {
            if (poly_edge.has_on(split_line.source()) || poly_edge.has_on(split_line.target()) ||
                split_line.has_on(poly_edge.source()) || split_line.has_on(poly_edge.target()))
                return 0;
        }
    }

    std::vector<Point> intersection_pts;
    for (const auto& poly_edge : polygon_edges) {
        auto result = CGAL::intersection(split_line, poly_edge);
        if (!result) continue;
        const CGAL::Object& obj = *result;
        if (const Point* p = CGAL::object_cast<Point>(&obj)) {
            if (point_eq(*p, split_line.source(), POINT_TOL) ||
                point_eq(*p, split_line.target(), POINT_TOL))
                continue;
            bool already = false;
            for (const auto& existing : intersection_pts)
                if (point_eq(*p, existing, POINT_TOL)) { already = true; break; }
            if (already) continue;
            if (is_polygon_vertex(*p, polygon_edges)) continue;
            intersection_pts.push_back(*p);
        } else if (const Edge* s = CGAL::object_cast<Edge>(&obj)) {
            intersection_pts.push_back(CGAL::midpoint(s->source(), s->target()));
        }
    }
    return static_cast<int>(intersection_pts.size());
}

static double compute_balance_score(const std::array<Triangle, 2>& child_tris,
                                    const std::vector<Edge>& edges)
{
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

static bool edge_lies_on_cut(const Edge& edge, const Edge& split_line)
{
    if (point_eq(edge.source(), edge.target())) return false;
    return split_line.has_on(edge.source()) && split_line.has_on(edge.target());
}

static std::pair<std::vector<Edge>, std::vector<Edge>> assign_edges(
    const std::array<Triangle, 2>& child_tris,
    const std::vector<Edge>& edges,
    const Edge& split_line,
    std::vector<LineSplitRecord>* line_splits_out)
{
    std::vector<Edge> child_edges[2];
    Triangle_cgal ct[2] = {
        Triangle_cgal(child_tris[0][0], child_tris[0][1], child_tris[0][2]),
        Triangle_cgal(child_tris[1][0], child_tris[1][1], child_tris[1][2])
    };
    Line cut_line = split_line.supporting_line();

    for (const auto& edge : edges) {
        if (edge_lies_on_cut(edge, split_line)) {
            child_edges[0].push_back(edge);
            child_edges[1].push_back(edge);
            continue;
        }
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
            if (line_splits_out && seg1_valid && seg2_valid) {
                LineSplitRecord rec;
                rec.original_edge = edge;
                rec.split_point = *p;
                rec.segment1 = seg1;
                rec.segment2 = seg2;
                line_splits_out->push_back(rec);
            }
            if (seg1_valid) place_edge_in_child(seg1, ct, child_edges);
            if (seg2_valid) place_edge_in_child(seg2, ct, child_edges);
        } else if (CGAL::object_cast<Edge>(&intersection_obj)) {
            child_edges[0].push_back(edge);
            child_edges[1].push_back(edge);
        }
    }
    return { child_edges[0], child_edges[1] };
}

static inline double am_rel_angle(double rx, double ry, double px, double py)
{
    double cross = rx * py - ry * px;
    double dot = rx * px + ry * py;
    return std::abs(std::atan2(cross, dot));
}

static std::vector<ScoredCut> am_scored_cuts(const Triangle& triangle,
                                             const std::vector<Edge>& edges)
{
    std::vector<ScoredCut> scored;
    if (edges.empty()) return scored;

    for (int c = 0; c < 3; ++c) {
        const Point& v = triangle[c];
        const Point& u = triangle[(c + 1) % 3];
        const Point& w = triangle[(c + 2) % 3];
        double vx = CGAL::to_double(v.x()), vy = CGAL::to_double(v.y());
        double rx = CGAL::to_double(u.x()) - vx, ry = CGAL::to_double(u.y()) - vy;
        double aw = am_rel_angle(rx, ry,
                                 CGAL::to_double(w.x()) - vx,
                                 CGAL::to_double(w.y()) - vy);
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
            if (sx * sx + sy * sy < 1e-24) a1 = a2;
            if (tx * tx + ty * ty < 1e-24) a2 = a1;
            seg_lo.push_back(std::min(a1, a2));
            seg_hi.push_back(std::max(a1, a2));
            if (a1 > 1e-9 && a1 < aw - 1e-9) cand.push_back({a1, e.source()});
            if (a2 > 1e-9 && a2 < aw - 1e-9) cand.push_back({a2, e.target()});
        }
        if (cand.empty()) continue;

        std::sort(cand.begin(), cand.end(),
                  [](const AngleEndpoint& a, const AngleEndpoint& b) { return a.angle < b.angle; });
        std::vector<AngleEndpoint> uniq;
        uniq.reserve(cand.size());
        for (const auto& ce : cand)
            if (uniq.empty() || std::abs(ce.angle - uniq.back().angle) > 1e-9)
                uniq.push_back(ce);

        if ((int)uniq.size() > 2 * MAX_PER_SIDE) {
            double mid = 0.5 * aw;
            int median_idx = 0;
            double bd = std::numeric_limits<double>::max();
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

        std::vector<double> emb_angles;
        emb_angles.reserve(seg_lo.size());
        for (size_t i = 0; i < seg_lo.size(); ++i)
            if (seg_hi[i] - seg_lo[i] < AEPS)
                emb_angles.push_back(0.5 * (seg_lo[i] + seg_hi[i]));
        std::sort(emb_angles.begin(), emb_angles.end());

        for (const auto& ce : uniq) {
            double x = ce.angle;
            int raw_left = (int)(std::upper_bound(sorted_hi.begin(), sorted_hi.end(), x) - sorted_hi.begin());
            int raw_right = n - (int)(std::lower_bound(sorted_lo.begin(), sorted_lo.end(), x) - sorted_lo.begin());
            int embedded = (int)(std::upper_bound(emb_angles.begin(), emb_angles.end(), x + AEPS)
                               - std::lower_bound(emb_angles.begin(), emb_angles.end(), x - AEPS));
            int left = raw_left - embedded;
            int right = raw_right - embedded;
            int cuts = n - left - right - embedded;
            if (cuts < 0) cuts = 0;
            double balance = (left + right > 0)
                ? 1.5 * std::abs(left - right) / double(left + right) : 1.5;
            double cost = cuts + balance - EMBED_BONUS * embedded;
            scored.push_back({ cost, c, ce.pt, cuts, embedded });
        }
    }

    std::sort(scored.begin(), scored.end(),
              [](const ScoredCut& a, const ScoredCut& b) { return a.cost < b.cost; });
    return scored;
}

static bool am_realize_cut(const Triangle& triangle, const ScoredCut& sc,
                           const std::vector<Edge>& all_polygon_edges,
                           std::array<Triangle, 2>& kids, Edge& split_line, int& priority)
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

static BestCutResult find_best_split_angle(const Triangle& triangle,
                                           const std::vector<Edge>& edges,
                                           const std::vector<Edge>& all_polygon_edges)
{
    BestCutResult result;
    std::vector<ScoredCut> scored = am_scored_cuts(triangle, edges);
    for (const auto& sc : scored) {
        std::array<Triangle, 2> kids;
        Edge sline;
        int prio;
        if (!am_realize_cut(triangle, sc, all_polygon_edges, kids, sline, prio))
            continue;
        result.found = true;
        result.child_triangles = kids;
        result.split_line = sline;
        result.priority = prio;
        return result;
    }
    return result;
}

static bool slope_split_greedy(
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
                if (!intersection_result) continue;
                const CGAL::Object& obj = *intersection_result;
                const Point* p = CGAL::object_cast<Point>(&obj);
                if (!p) continue;

                Triangle t1 = { vtx, opp0, *p };
                Triangle t2 = { vtx, *p, opp1 };
                if (area2d(t1) <= AREA_TOL || area2d(t2) <= AREA_TOL) continue;

                Edge potential_split_line(vtx, *p);
                bool is_poly_edge = split_line_is_polygon_edge(potential_split_line, all_polygon_edges);
                int current_intersections = is_poly_edge ? 0
                    : count_intersections_improved(potential_split_line, all_polygon_edges);
                double balance = compute_balance_score({ t1, t2 }, candidate_edges);

                bool is_better = false;
                if (is_poly_edge && !best_is_poly_edge) {
                    is_better = true;
                } else if (is_poly_edge == best_is_poly_edge) {
                    if (current_intersections < min_intersections)
                        is_better = true;
                    else if (current_intersections == min_intersections && balance < best_balance)
                        is_better = true;
                }
                if (!is_better) continue;

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
    return best_split_found;
}

static void order_children_left_right(const Edge& cut, Triangle& left, Triangle& right)
{
    Line L = cut.supporting_line();
    Point m = CGAL::centroid(left[0], left[1], left[2]);
    if (L.oriented_side(m) == CGAL::ON_NEGATIVE_SIDE)
        std::swap(left, right);
}

static void td_xy(const Point& p, double o[2])
{
    o[0] = CGAL::to_double(p.x());
    o[1] = CGAL::to_double(p.y());
}

static void td_fill_tri(const Triangle& tri, TDTri& out)
{
    for (int i = 0; i < 3; ++i) td_xy(tri[i], out.v[i]);
}

static void td_fill_edge(const Edge& e, TDEdge& out)
{
    td_xy(e.source(), out.a);
    td_xy(e.target(), out.b);
}

static bool td_start_ok(const char* start, std::string& out, std::string& err)
{
    out = start ? start : "";
    for (char c : out) {
        if (c != '0' && c != '1') {
            err = "start must be a bitstring of 0/1 (empty for a fresh triangle)";
            return false;
        }
    }
    return true;
}

static std::vector<SplitEvent> run_tessellate(const Triangle& initial,
                                              const std::vector<Edge>& segments,
                                              const std::string& start_addr,
                                              bool verbose)
{
    std::queue<WorkItem> queue;
    queue.push({ initial, segments, start_addr });
    const int max_iter = std::max(100000, static_cast<int>(segments.size()) * 10);
    int iteration = 0;
    std::vector<SplitEvent> events;

    while (!queue.empty() && iteration < max_iter) {
        WorkItem item = queue.front();
        queue.pop();
        const Triangle tri = item.tri;
        const std::vector<Edge> edges = std::move(item.edges);
        const std::string addr = std::move(item.addr);
        iteration++;

        if (verbose && (iteration == 1 || iteration % 100 == 0))
            std::cerr << "[progress] iteration=" << iteration
                      << " splits=" << events.size()
                      << " queue=" << queue.size() << "\n";

        std::vector<Edge> current_edges;
        for (const auto& e : edges)
            if (!edge_on_triangle_side(e, tri))
                current_edges.push_back(e);
        current_edges = filter_degenerate_edges(current_edges);
        if (current_edges.empty())
            continue;

        std::array<Triangle, 2> kids;
        Edge split_line;
        BestCutResult best = find_best_split_angle(tri, current_edges, current_edges);
        if (best.found) {
            kids = best.child_triangles;
            split_line = best.split_line;
        } else if (!slope_split_greedy(tri, current_edges, current_edges, kids, split_line)) {
            continue;
        }

        order_children_left_right(split_line, kids[0], kids[1]);
        std::vector<LineSplitRecord> ls;
        auto assigned = assign_edges(kids, current_edges, split_line, &ls);

        SplitEvent ev;
        ev.id = addr;
        ev.parent = tri;
        ev.child0 = kids[0];
        ev.child1 = kids[1];
        ev.cut = split_line;
        ev.line_splits = std::move(ls);
        events.push_back(std::move(ev));

        if (area2d(kids[0]) > AREA_TOL)
            queue.push(WorkItem{ kids[0], filter_degenerate_edges(assigned.first), addr + "0" });
        if (area2d(kids[1]) > AREA_TOL)
            queue.push(WorkItem{ kids[1], filter_degenerate_edges(assigned.second), addr + "1" });
    }
    return events;
}

TDTessellation tessellate_from_xy(const double* seg_xy, int nseg,
                                  const double* t_xy, const char* start,
                                  bool verbose)
{
    TDTessellation R;
    if (!t_xy || nseg < 0) {
        R.error = "t must be a (3,2) triangle; segments must be a non-negative count";
        return R;
    }
    if (nseg > 0 && !seg_xy) {
        R.error = "segments pointer is null";
        return R;
    }
    std::string start_addr;
    if (!td_start_ok(start, start_addr, R.error))
        return R;

    std::vector<Edge> segments;
    segments.reserve((size_t)nseg);
    for (int i = 0; i < nseg; ++i) {
        Point a(seg_xy[4 * i],     seg_xy[4 * i + 1]);
        Point b(seg_xy[4 * i + 2], seg_xy[4 * i + 3]);
        if (point_eq(a, b)) continue;
        segments.emplace_back(a, b);
    }

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

    std::vector<SplitEvent> events = run_tessellate(initial, segments, start_addr, verbose);
    R.events.reserve(events.size());
    for (const auto& rec : events) {
        TDSplitEvent e;
        e.id = rec.id;
        e.child0 = rec.id + "0";
        e.child1 = rec.id + "1";
        td_fill_tri(rec.parent, e.parent);
        td_fill_edge(rec.cut, e.cut);
        td_fill_tri(rec.child0, e.child[0]);
        td_fill_tri(rec.child1, e.child[1]);
        e.line_splits.reserve(rec.line_splits.size());
        for (const auto& ls : rec.line_splits) {
            TDLineSplit s;
            td_fill_edge(ls.original_edge, s.original_edge);
            td_xy(ls.split_point, s.split_point);
            td_fill_edge(ls.segment1, s.segment1);
            td_fill_edge(ls.segment2, s.segment2);
            e.line_splits.push_back(s);
        }
        R.events.push_back(std::move(e));
    }
    return R;
}
