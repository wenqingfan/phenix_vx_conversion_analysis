AUTOMAKE_OPTIONS = foreign

INCLUDES = -I$(includedir) -I$(OFFLINE_MAIN)/include

lib_LTLIBRARIES = \
libPhotonConversionAnalysisReco.la

noinst_HEADERS = \
ReconstructionLinkDef.h

pkginclude_HEADERS = \
Reconstruction.h

libPhotonConversionAnalysisReco_la_SOURCES = \
Reconstruction.cpp \
Reconstruction_Dict.C

libPhotonConversionAnalysisReco_la_LIBADD = \
$(ROOT_LIBS) \
-lPhotonConversionAnalysisEvent

# Rule for generating CINT dictionaries from class headers.
Reconstruction_Dict.C: \
  Reconstruction.h \
  ReconstructionLinkDef.h
	rootcint -f $@ -c $(DEFAULT_INCLUDES) $(INCLUDES) $^

clean-local:
	rm -f *Dict*

