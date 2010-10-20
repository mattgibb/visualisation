AM_CPPFLAGS = $(BOOST_CPPFLAGS) -O2
AM_LDFLAGS = $(BOOST_GRAPH_LDFLAGS) $(BOOST_PROGRAM_OPTIONS_LDFLAGS)
LDADD = $(BOOST_GRAPH_LIBS) $(BOOST_PROGRAM_OPTIONS_LIBS)

#
bin_PROGRAMS = minpkjpath tetgentovtk extractpmjs genbundles
minpkjpath_SOURCES = main.cpp
tetgentovtk_SOURCES = tetgen_to_vtk.cpp
extractpmjs_SOURCES = extract_pmjs.cpp
genbundles_SOURCES = gen_bundles.cpp