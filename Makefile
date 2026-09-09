# Linux build for TriDecomp (CGAL + shapelib).
# Prefers system packages; falls back to a user-local prefix extracted from
# Ubuntu .deb files (see `make deps`).

CXX      ?= g++
CXXFLAGS ?= -std=c++17 -O2 -frounding-math -pthread
LDFLAGS  ?= -pthread

SYSTEM_CGAL := $(wildcard /usr/include/CGAL)
LOCAL_PREFIX ?= $(HOME)/.local/tridecomp-deps

ifeq ($(SYSTEM_CGAL),)
  INCFLAGS = -I$(LOCAL_PREFIX)/usr/include \
             -I$(LOCAL_PREFIX)/usr/include/x86_64-linux-gnu
  LIBDIR   = $(LOCAL_PREFIX)/usr/lib/x86_64-linux-gnu
  # Static GMP/MPFR/shapelib so the binary does not need LD_LIBRARY_PATH.
  LIBS     = $(LIBDIR)/libshp.a $(LIBDIR)/libmpfr.a $(LIBDIR)/libgmp.a
  SHP_LIB  = $(LIBDIR)/libshp.a
else
  INCFLAGS =
  LIBS     = -lshp -lgmp -lmpfr
  SHP_LIB  = -lshp
endif

.PHONY: all clean deps test

all: TriDecomp

TriDecomp: TriDecomp.cpp
	$(CXX) $(CXXFLAGS) $(INCFLAGS) -o $@ $< $(LDFLAGS) $(LIBS)

# Download and extract Ubuntu packages into $(LOCAL_PREFIX) (no root required).
deps:
	mkdir -p /tmp/tridecomp-debs "$(LOCAL_PREFIX)"
	cd /tmp/tridecomp-debs && apt-get download \
		libcgal-dev libgmp-dev libgmpxx4ldbl libmpfr-dev libshp-dev libshp4
	cd /tmp/tridecomp-debs && for f in *.deb; do dpkg-deb -x "$$f" "$(LOCAL_PREFIX)"; done

tools/make_test_shp: tools/make_test_shp.cpp
	mkdir -p tools
	$(CXX) -std=c++17 -O2 $(INCFLAGS) -o $@ $< $(SHP_LIB)

test: TriDecomp tools/make_test_shp
	mkdir -p testdata/run
	./tools/make_test_shp testdata/unit_square
	./TriDecomp testdata/unit_square.shp --list-states
	cd testdata/run && ../../TriDecomp ../unit_square.shp --state 48 --county 113

clean:
	rm -f TriDecomp tools/make_test_shp
	rm -rf testdata
