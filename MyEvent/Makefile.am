AUTOMAKE_OPTIONS = foreign

INCLUDES = -I$(includedir) -I$(OFFLINE_MAIN)/include

lib_LTLIBRARIES = \
libPhotonConversionAnalysisEvent.la

noinst_HEADERS = \
MyEventLinkDef.h

pkginclude_HEADERS = \
MyEvent.h

libPhotonConversionAnalysisEvent_la_SOURCES = \
MyEvent.C \
MyEvent_Dict.C

libPhotonConversionAnalysisEvent_la_LIBADD = \
$(ROOT_LIBS)

# Rule for generating CINT dictionaries from class headers.
MyEvent_Dict.C: \
  MyEvent.h \
  MyEventLinkDef.h 
	rootcint -f $@ -c $(DEFAULT_INCLUDES) $(INCLUDES) $^

clean-local:
	rm -f *Dict*
