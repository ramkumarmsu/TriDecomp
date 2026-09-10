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

.PHONY: all clean deps test python

all: TriDecomp

TriDecomp: TriDecomp.cpp
	$(CXX) $(CXXFLAGS) $(INCFLAGS) -o $@ $< $(LDFLAGS) $(LIBS)

PY ?= $(HOME)/miniconda3/envs/tridecomp-py/bin/python
PY_INCS = $(shell $(PY) -c "import sysconfig,numpy; print('-I'+sysconfig.get_path('include')+' -I'+numpy.get_include())")
PY_EXT  = $(shell $(PY) -c "import sysconfig; print(sysconfig.get_config_var('EXT_SUFFIX'))")
PY_LIBDIR = $(shell $(PY) -c "import sysconfig; print(sysconfig.get_config_var('LIBDIR'))")
PY_LDLIB  = $(shell $(PY) -c "import sys; print('python%d.%d' % sys.version_info[:2])")

# Library engine is tessellate.cpp (CGAL + GMP/MPFR only; no shapelib).
ifeq ($(SYSTEM_CGAL),)
  PY_LIBS = $(LIBDIR)/libmpfr.a $(LIBDIR)/libgmp.a
else
  PY_LIBS = -lgmp -lmpfr
endif

python: python/tridecomp/_cpp$(PY_EXT)

python/tridecomp/_cpp$(PY_EXT): tessellate.cpp python/tridecomp/_cpp.cpp tridecomp_api.h
	$(CXX) -shared -fPIC $(CXXFLAGS) \
	    $(INCFLAGS) $(PY_INCS) -I. \
	    -o $@ python/tridecomp/_cpp.cpp tessellate.cpp \
	    $(LDFLAGS) $(PY_LIBS) -L$(PY_LIBDIR) -l$(PY_LDLIB)

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
	rm -f python/tridecomp/_cpp*.so python/tridecomp/*.pyc
	rm -rf testdata python/tridecomp/__pycache__
