AUTOMAKE_OPTIONS = foreign

INCLUDES = -I$(includedir) 

lib_LTLIBRARIES = \
libMyPion.la

noinst_HEADERS = \
MyPionLinkDef.h

include_HEADERS = \
MyPion.h 

libMyPion_la_SOURCES = \
MyPion.C \
MyPion_Dict.C

libMyPion_la_LIBADD = \
$(ROOT_LIBS)

# Rule for generating CINT dictionaries from class headers.
MyPion_Dict.C: \
  MyPion.h \
  MyPionLinkDef.h 
	rootcint -f $@ -c $(DEFAULT_INCLUDES) $(INCLUDES) $^

clean-local:
	rm -f *Dict*
