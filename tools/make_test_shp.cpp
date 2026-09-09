// Write a tiny ESRI shapefile (square polygon + STATEFP/COUNTYFP) for smoke tests.
#include <shapefil.h>
#include <cstdio>
#include <cstdlib>
#include <cstring>

int main(int argc, char** argv) {
    const char* path = (argc > 1) ? argv[1] : "testdata/unit_square";

    SHPHandle shp = SHPCreate(path, SHPT_POLYGON);
    if (!shp) {
        std::fprintf(stderr, "SHPCreate failed for %s\n", path);
        return 1;
    }

    // Simple CCW unit square in lon/lat (near the origin, valid WGS84-ish).
    double xs[5] = {-97.0, -96.0, -96.0, -97.0, -97.0};
    double ys[5] = { 32.0,  32.0,  33.0,  33.0,  32.0};
    int n = 5;
    SHPObject* obj = SHPCreateSimpleObject(SHPT_POLYGON, n, xs, ys, nullptr);
    if (!obj) {
        std::fprintf(stderr, "SHPCreateSimpleObject failed\n");
        SHPClose(shp);
        return 1;
    }
    SHPWriteObject(shp, -1, obj);
    SHPDestroyObject(obj);
    SHPClose(shp);

    DBFHandle dbf = DBFCreate(path);
    if (!dbf) {
        std::fprintf(stderr, "DBFCreate failed for %s\n", path);
        return 1;
    }
    int f_state  = DBFAddField(dbf, "STATEFP",  FTString, 2, 0);
    int f_county = DBFAddField(dbf, "COUNTYFP", FTString, 3, 0);
    int f_name   = DBFAddField(dbf, "NAME",     FTString, 32, 0);
    if (f_state < 0 || f_county < 0 || f_name < 0) {
        std::fprintf(stderr, "DBFAddField failed\n");
        DBFClose(dbf);
        return 1;
    }
    DBFWriteStringAttribute(dbf, 0, f_state,  "48");
    DBFWriteStringAttribute(dbf, 0, f_county, "113");
    DBFWriteStringAttribute(dbf, 0, f_name,   "TestSquare");
    DBFClose(dbf);

    std::printf("Wrote %s.shp / .shx / .dbf\n", path);
    return 0;
}
