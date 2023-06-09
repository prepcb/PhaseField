## Top-level HDF5 Makefile(.in)
##
# Copyright by The HDF Group.
# Copyright by the Board of Trustees of the University of Illinois.
# All rights reserved.
#
# This file is part of HDF5.  The full HDF5 copyright notice, including
# terms governing use, modification, and redistribution, is contained in
# the files COPYING and Copyright.html.  COPYING can be found at the root
# of the source code distribution tree; Copyright.html can be found at the
# root level of an installed copy of the electronic HDF5 document set and
# is linked from the top-level documents page.  It can also be found at
# http://hdfgroup.org/HDF5/doc/Copyright.html.  If you do not have
# access to either file, you may request a copy from help@hdfgroup.org.
## 
## This makefile mostly just reinvokes make in the various subdirectories
## but does so in the correct order.  You can alternatively invoke make from
## each subdirectory manually.
##
top_srcdir=..
top_builddir=..
srcdir=.

##------------------------------------------------------------ -*- makefile -*-
# Copyright by The HDF Group.
# Copyright by the Board of Trustees of the University of Illinois.
# All rights reserved.
#
# This file is part of HDF5.  The full HDF5 copyright notice, including
# terms governing use, modification, and redistribution, is contained in
# the files COPYING and Copyright.html.  COPYING can be found at the root
# of the source code distribution tree; Copyright.html can be found at the
# root level of an installed copy of the electronic HDF5 document set and
# is linked from the top-level documents page.  It can also be found at
# http://hdfgroup.org/HDF5/doc/Copyright.html.  If you do not have
# access to either file, you may request a copy from help@hdfgroup.org.
##
##
## The following section of this makefile comes from the
## `./config/commence' file which was generated with config.status
## from `./config/commence.in'.
##-----------------------------------------------------------------------------

## Things that Make needs
.SUFFIXES:
.SUFFIXES: .c .o .lo


## Directories to search
VPATH=$(srcdir) $(top_builddir)/src $(top_srcdir)/src

## Programs
SHELL=/bin/sh
CC=mpicc
## H5_CPPFLAGS and H5_CFLAGS are used for extra flags used to build HDF5
## (flags that don't need to be included in h5cc, for instance).
CFLAGS= -fPIC -ansi -O3 -fomit-frame-pointer -finline-functions 
CPPFLAGS=-UH5_DEBUG_API -DNDEBUG -D_LARGEFILE_SOURCE -D_LARGEFILE64_SOURCE -D_FILE_OFFSET_BITS=64 -D_POSIX_SOURCE -D_BSD_SOURCE  
H5_CFLAGS= -fPIC -ansi -O3 -fomit-frame-pointer -finline-functions
H5_CPPFLAGS=-UH5_DEBUG_API -DNDEBUG -D_LARGEFILE_SOURCE -D_LARGEFILE64_SOURCE -D_FILE_OFFSET_BITS=64 -D_POSIX_SOURCE -D_BSD_SOURCE 
LDFLAGS=
ROOT=/usr/not-backed-up/prepcb9Feb2023/FreshCampfire/hdf5-1.6.10
LIBS=-lz -lm 
AR=ar
RANLIB=ranlib
PERL=perl
RM=rm -f
CP=cp
INSTALL=/usr/bin/install -c
INSTALL_PROGRAM=${INSTALL}
INSTALL_DATA=${INSTALL} -m 644
PARALLEL=yes
RUNSERIAL=
RUNPARALLEL=/usr/lib64/openmpi/bin/mpirun -np $${NPROCS:=3}
RUNTEST=$(RUNSERIAL)
TIME=time
TRACE=:

## Installation points
## Insert DESTDIR to individual installation directory because they can be
## defined independent of prefix.
prefix=/usr/not-backed-up/prepcb9Feb2023/FreshCampfire/hdf5-1.6.10
exec_prefix=${prefix}
bindir=$(DESTDIR)${exec_prefix}/bin
libdir=$(DESTDIR)${exec_prefix}/lib
includedir=$(DESTDIR)${prefix}/include
docdir=$(DESTDIR)${prefix}/doc
PUB_LIB=$(LIB)

## Shared libraries
LT_STATIC_EXEC=
DYNAMIC_DIRS=
LT=$(top_builddir)/libtool
LT_COMPILE=$(LT) --mode=compile $(CC)
LT_LINK_LIB=$(LT) --mode=link $(CC) -rpath $(libdir) $(DYNAMIC_DIRS)
LT_LINK_EXE=$(LT) --mode=link $(CC) $(LT_STATIC_EXEC) -dlopen self -rpath $(libdir) $(DYNAMIC_DIRS)
LT_RUN=$(LT) --mode=execute
LT_INSTALL_PROG=$(LT) --mode=install $(INSTALL_PROGRAM)
LT_INSTALL_LIB=$(LT) --mode=install $(INSTALL_DATA)
LT_UNINSTALL=$(LT) --mode=uninstall $(RM)

## Optional variables. We must declare them here because Irix pmake
## complains if it sees a reference to a variable which has never been
## defined. The main makefile is free to redefine these to something else.
DOCDIR=$(docdir)
SUBDOCDIR=
EXAMPLEDIR=$(docdir)/hdf5/examples/c
LIB=
LIB_SRC=
LIB_OBJ=
PUB_HDR=
PUB_PROGS=
PUB_DOCS=
PUB_SUBDOCS=
PROGS=
TEST_PROGS=
PRIV_PROGS=
TEST_PROGS_PARA=
TEST_FLAGS=
TEST_SCRIPTS=
TEST_SCRIPTS_PARA=
AUX_LIB=
EXAMPLE_PROGS=
SUBDIRS=
LIBHDF5=

## The default is to build the library and/or programs.  We must build
## them sequentially.
all:
	$(MAKE) lib
	$(MAKE) progs
	$(MAKE) tests

## The following rules insure that the Makefile is up-to-date by rerunning
## various autoconf components (although not all versions of make assume
## that the makefile is implicitly a target).  We use time stamp files to
## keep track of the most recent update of H5config.h.in and H5config.h
## because autoheader and config.status don't update the modification time
## if the contents don't change.
## Invoke it by "gmake reconfigure".
##
## Graphically, the dependencies are:
##
##		         configure.in
##                           |  |
##         +-----------------+  +------------+
##         |                                 |
##      stamp1                           configure
##   (H5config.h.in)                         |
##       |  |                                |
##       |  +---------------+    +-----------+
##       |                  |    |
##       |               config.status
##       |                     |
##       |  +------------------+
##       |  |
##      stamp2
##    (H5config.h)                   Makefile.in et al
##         |                                 |
##         +------------------+  +-----------+
##                            |  |
##                          Makefile
##
## A side effect of updating stamp1 is to generate H5config.h.in and a
## side effect of updating stamp2 is to generate H5config.h.  When using
## a version of make that doesn't treat the makefile as the initial target
## the user may want to occassionally type `make Makefile' in any source
## directory.
## The `Makefile' target has been renamed to `reconfigure' so that the
## autoconf and make depend processes do not start up automatically.
## One must do `make reconfigure' explicitedly to start the process.
## (When srcdir is used and if more than one machines are running,
## this automatic Makefile/autoconf can get things unstable.)
##
STAMP1=$(top_builddir)/config/stamp1
STAMP2=$(top_builddir)/config/stamp2

MAKEFILE_PARTS=$(srcdir)/Makefile.in		\
	       $(top_srcdir)/config/commence.in	\
	       $(top_srcdir)/config/conclude.in	\
	       $(top_srcdir)/config/depend1.in	\
	       $(top_srcdir)/config/depend2.in	\
	       $(top_srcdir)/config/depend3.in	\
	       $(top_srcdir)/config/depend4.in	\
	       $(top_srcdir)/config/dependN.in

$(STAMP1): $(top_srcdir)/configure.in
	touch $(STAMP1)
	-cd $(top_srcdir); autoheader

$(STAMP2): $(STAMP1) $(top_builddir)/config.status
	touch $(STAMP2)
	-cd $(top_builddir); \
	   CONFIG_FILES= CONFIG_HEADERS=src/H5config.h ./config.status

$(top_srcdir)/configure: $(top_srcdir)/configure.in
	-cd $(top_srcdir); autoconf

$(top_builddir)/config.status: $(top_srcdir)/configure $(STAMP1)
	-cd $(top_builddir); ./config.status --recheck

# rerun the autoconf process if any configure components have changed. 
reconfigure: $(MAKEFILE_PARTS) $(STAMP2)
	-cd $(top_builddir); CONFIG_HEADERS= ./config.status
	@for sub in   dummy; do					\
	  if test $$sub != dummy; then					\
	    -cd $$sub; $(MAKE) reconfigure;				\
	  fi;								\
	done;

## make help target
help:
	$(top_srcdir)/bin/makehelp

##-----------------------------------------------------------------------------
## The following section of this makefile comes from the middle of
## `Makefile.in' from this directory.  It was generated by running
## `config.status'.
##-----------------------------------------------------------------------------

# Subdirectories in build-order (not including 'examples')
SUBDIRS=src test tools

##############################################################################
##		    T A R G E T S
##
## all:		Build libraries, header files, tests, and programs in the
##	 	various subdirectories but does not run tests or install the
##		library, header files, or programs.  The components can be
##		built individually with the targets lib, progs, and tests.
## check:	Test the uninstalled library to make sure it works. You may
##		also say `test' or `_test' (`test' doesn't work from the top
##		level directory for some versions of make because `test' is
##		also a directory).
## install:	Installs libraries, header files, programs, and documentation
##		in the various directories under the prefix directory (lib,
##		include, bin, man, info).  Use the `--prefix=PATH' option
##		to `configure' (or `config.status') or say `--help' for
##		other alternatives.  The default prefix is `/usr/local'.
## uninstall:	Delete all the installed files that the `install' target
##		created (but not the noninstalled files such as `make all'
##		created).
## clean:	Removes temporary files except those that record the
##		configuration and those that are part of the distribution.
## mostlyclean:	Like `clean' except it doesn't delete a few files like
##		libraries, programs, and/or generated header files because
##		regenerating them is rarely necessary and takes a lot of time.
## distclean:	Deletes all files that are created by configuring or building
##		HDF5.  If you have unpacked the source and built HDF5 without
##		creating any other files, then `make distclean' will leave
##		only the files that were in the distrubution.
## maintainer-clean:
##		Like `distclean' except it deletes more files.  It deletes
##		all generated files.  This target is not intended for normal
##		users; it deletes files that may require special tools to
##		rebuild.
## TAGS:	Updates the tags table for this program.
## dep depend:	Builds dependencies in all subdirectories.  These targets
##		might not be available on certain combinations of make
##		programs and C compilers.  At the other extreme, the GNU
##		make used in combination with gcc will maintain dependency
##		information automatically.
lib progs check check-s check-p check-clean test _test uninstall:
	@:; for d in $(SUBDIRS); do				      \
	   (cd $$d && $(MAKE) $@) || exit 1;				      \
	done

tests TAGS dep depend:
	@:; for d in $(SUBDIRS); do			              \
	   (cd $$d && $(MAKE) $@) || exit 1;				      \
	done

install-all:
	@$(MAKE) install
install:
	@:; for d in $(SUBDIRS); do				      \
	   (cd $$d && $(MAKE) $@) || exit 1;				      \
	done
	@$(MAKE) install-examples

install-examples:
	(cd examples && $(MAKE) $@) || exit 1;

uninstall-examples:
	(cd examples && $(MAKE) $@) || exit 1;

check-install:
	(cd examples && $(MAKE) check) || exit 1

.PHONY: all lib progs test _test install uninstall dep depend clean           \
	mostlyclean distclean maintainer-clean

clean mostlyclean:
	@:; for d in $(SUBDIRS) examples; do	                      \
	   (cd $$d && $(MAKE) $@);					      \
	done

distclean:
	@:; for d in $(SUBDIRS) examples; do	                      \
	   (cd $$d && $(MAKE) $@);					      \
	done
	-$(RM)  Makefile

maintainer-clean:
	@echo "This target is intended for maintainers to use;"
	@echo "it deletes files that may require special tools to rebuild."
	@:; for d in $(SUBDIRS); do		                      \
	   (cd $$d && $(MAKE) $@);					      \
	done


