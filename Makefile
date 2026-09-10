# Build the Python tessellate extension (CGAL EPEC + GMP/MPFR).
# Prefers system CGAL; falls back to a user-local prefix from Ubuntu .deb files
# (see `make deps`).

CXX      ?= g++
CXXFLAGS ?= -std=c++17 -O2 -frounding-math -pthread
LDFLAGS  ?= -pthread

SYSTEM_CGAL := $(wildcard /usr/include/CGAL)
LOCAL_PREFIX ?= $(HOME)/.local/tridecomp-deps

ifeq ($(SYSTEM_CGAL),)
  INCFLAGS = -I$(LOCAL_PREFIX)/usr/include \
             -I$(LOCAL_PREFIX)/usr/include/x86_64-linux-gnu
  LIBDIR   = $(LOCAL_PREFIX)/usr/lib/x86_64-linux-gnu
  LIBS     = $(LIBDIR)/libmpfr.a $(LIBDIR)/libgmp.a
else
  INCFLAGS =
  LIBS     = -lgmp -lmpfr
endif

PY ?= $(HOME)/miniconda3/envs/tridecomp-py/bin/python
PY_INCS = $(shell $(PY) -c "import sysconfig,numpy; print('-I'+sysconfig.get_path('include')+' -I'+numpy.get_include())")
PY_EXT  = $(shell $(PY) -c "import sysconfig; print(sysconfig.get_config_var('EXT_SUFFIX'))")
PY_LIBDIR = $(shell $(PY) -c "import sysconfig; print(sysconfig.get_config_var('LIBDIR'))")
PY_LDLIB  = $(shell $(PY) -c "import sys; print('python%d.%d' % sys.version_info[:2])")

.PHONY: all python deps clean

all: python

python: python/tridecomp/_cpp$(PY_EXT)

python/tridecomp/_cpp$(PY_EXT): tessellate.cpp python/tridecomp/_cpp.cpp tridecomp_api.h
	$(CXX) -shared -fPIC $(CXXFLAGS) \
	    $(INCFLAGS) $(PY_INCS) -I. \
	    -o $@ python/tridecomp/_cpp.cpp tessellate.cpp \
	    $(LDFLAGS) $(LIBS) -L$(PY_LIBDIR) -l$(PY_LDLIB)

# Download and extract Ubuntu packages into $(LOCAL_PREFIX) (no root required).
deps:
	mkdir -p /tmp/tridecomp-debs "$(LOCAL_PREFIX)"
	cd /tmp/tridecomp-debs && apt-get download \
		libcgal-dev libgmp-dev libgmpxx4ldbl libmpfr-dev
	cd /tmp/tridecomp-debs && for f in *.deb; do dpkg-deb -x "$$f" "$(LOCAL_PREFIX)"; done

clean:
	rm -f python/tridecomp/_cpp*.so python/tridecomp/*.pyc
	rm -rf python/tridecomp/__pycache__
