#! /bin/sh
#
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


#
# Tests for the embedded library information feature.
# Part 1:
# Verify the HDF5 library does contains an exact copy of the content of the
# libhdf5.settings file.
# Part 2:
# If executable is linked with the static hdf5 library (how to determine?),
# verify an executable indeed contains an exact copy of hte content of the
# libhdf5.settings file.
#
# Programmer: Albert Cheng
#	      Sep 18, 2009

# Determine the configure options of the hdf5 library and executables.

Shared_Lib=no
Static_Lib=yes
Static_exec=no


# Print a line-line message left justified in a field of 70 characters.
#
LINEMSG() {
   SPACES="                                                               "
   echo "Check file $* $SPACES" | cut -c1-70 | tr -d '\012'
}


# Print a "SKIP" message
SKIP() {
    LINEMSG $*
    echo  " -SKIP-"
}
  
# Function definitions
CHECK_LIBINFO(){
    LINEMSG $1
    if strings $1 | grep "SUMMARY OF THE HDF5 CONFIGURATION" > /dev/null; then
	echo " PASSED"
    else
	echo " FAILED"
	nerrors=`expr $nerrors + 1`
    fi
}


# MAIN Body
nerrors=0
H5_HAVE_EMBEDDED_LIBINFO=`grep '#define H5_HAVE_EMBEDDED_LIBINFO ' ../src/H5pubconf.h`

# Skip the rest if embedded-libinfo is not enabled.
if [ -z "$H5_HAVE_EMBEDDED_LIBINFO" ]; then
    echo "embedded-libinfo is not enabled. Test skipped."
    exit 0
fi

# The location of HDF library file(s) depends on whether shared lib is
# built too.
if [ -n $Shared_Lib ]; then
   h5libdir=../src/.libs
else
   h5libdir=../src
fi 

# Different OS uses different naming for shared libs.
case `uname -s` in
    Darwin)	# MacOS
        shlibsuffix=.dylib
        break
        ;;
    *)		# default
        shlibsuffix=.so
        break
        ;;
esac

h5libsettings=../src/libhdf5.settings 

# Part 1:
# Verify the HDF5 library does contains an exact copy of the content of the
# libhdf5.settings file.
# Check dynamic library file if built.
if [ x-$Shared_Lib = x-yes ]; then
    CHECK_LIBINFO ${h5libdir}/libhdf5${shlibsuffix}
else
    SKIP ${h5libdir}/libhdf5${shlibsuffix}
fi

# Though rare, libhdf5.a may not have been built.
if [ x-$Static_Lib = x-yes ]; then
    CHECK_LIBINFO ${h5libdir}/libhdf5.a
else
    SKIP ${h5libdir}/libhdf5.a
fi

# Check if executables has the lib information only if shared lib is not
# built or static-exec is used.  (Don't care static-exec since it affects
# tools binary only.)
if [ x-$Shared_Lib != x-yes ]; then
    CHECK_LIBINFO testhdf5
else
    SKIP testhdf5
fi


if [ $nerrors -gt 0 ]; then
    echo "***$nerrors errors encountered***"
    exit 1
else
    echo "No error encountered"
    exit 0
fi
