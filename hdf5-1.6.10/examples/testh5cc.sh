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
# Tests for the h5cc compiler tool
# Created: Albert Cheng, 2007/3/13
#
# Modification:
#	Albert Cheng, 2007/4/13
#	Added -shlib tests and verbose control.
#

# Initializations
# Where the tool is installed.
prefix="${prefix:-/usr/not-backed-up/prepcb9Feb2023/FreshCampfire/hdf5-1.6.10}"
PARALLEL=yes		# Am I in parallel mode?
AR="ar"
RANLIB="ranlib"
if [ "$PARALLEL" = no ]; then
    H5TOOL="h5cc"           	# The tool name
else
    H5TOOL="h5pcc"               # The tool name
fi
H5TOOL_BIN="${prefix}/bin/${H5TOOL}"   # The path of the tool binary

CMP='cmp -s'
DIFF='diff -c'

nerrors=0
verbose=${HDF5_VERBOSE:-1}      # 0: none; 1: default; 2: chatty; 3: everything
test $verbose -gt 2 && set -x

# setup my machine information.
myos=`uname -s`
myhostnama=`uname -n`

# The build (current) directory might be different than the source directory.
if test -z "$srcdir"; then
   srcdir=.
fi

# Generate some source files and library for tests.
suffix=c		# source file suffix
hdf5main=${H5TOOL}_hdf5main.$suffix
hdf5main_o=${H5TOOL}_hdf5main.o
appmain=${H5TOOL}_appmain.$suffix
appmain_o=${H5TOOL}_appmain.o
prog1=${H5TOOL}_prog1.$suffix
prog1_o=${H5TOOL}_prog1.o
prog2=${H5TOOL}_prog2.$suffix
prog2_o=${H5TOOL}_prog2.o
applib=libapp${H5TOOL}.a

# short hands
# Caution: if some *.h5 files must be cleaned here, list them by names.
# Don't use the wildcard form of *.h5 as it will wipe out even *.h5 generated
# by otehr test programs. This will cause a racing condition error when
# parallel make (e.g., gmake -j 4) is used.
temp_SRC="$hdf5main $appmain $prog1 $prog2"
temp_OBJ=`echo $temp_SRC | sed -e "s/\.${suffix}/.o/g"`
temp_FILES="a.out $applib"

# Generate appmain:
# An application Main that calls hdf5 and application's own functions.
cat > $appmain <<EOF
#include "hdf5.h"
#define H5FILE_NAME        "tmp.h5"
int
main (void)
{
    hid_t       file;         	/* file and dataset handles */

    /*
     * Create a new file using H5F_ACC_TRUNC access,
     * default file creation properties, and default file
     * access properties.
     */
    sub1();
    sub2();
    file = H5Fcreate(H5FILE_NAME, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
    H5Fclose(file);

    printf("HDF5 C Sample program ran successfully. File %s generated.\n", H5FILE_NAME);
    remove(H5FILE_NAME);
 
    return 0;
}     
EOF

# generate prog1
cat > $prog1 <<EOF
sub1(void)
{
    printf("in sub1\n");
}
EOF

# generate prog2
cat > $prog2 <<EOF
sub2(void)
{
    printf("in sub2\n");
}
EOF

# Generate HDF5 Main Program:
# An HDF5 sample program that calls hdf5 functions.
cat > $hdf5main <<EOF
#include "hdf5.h"
#define H5FILE_NAME        "tmp.h5"
int
main (void)
{
    hid_t       file;         	/* file and dataset handles */

    /*
     * Create a new file using H5F_ACC_TRUNC access,
     * default file creation properties, and default file
     * access properties.
     */
    file = H5Fcreate(H5FILE_NAME, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
    H5Fclose(file);

    printf("HDF5 C Sample program ran successfully. File %s generated.\n", H5FILE_NAME);
    remove(H5FILE_NAME);
 
    return 0;
}     
EOF


# Parse option
#   None

# Print a line-line message left justified in a field of 70 characters
# beginning with the word "Testing".
#
TESTING() {
   SPACES="                                                               "
   echo "Testing $* $SPACES" | cut -c1-70 | tr -d '\012'
}


# Debug printing
# Change : to echo to print the debug statement
DPRINT() {
    : $*
}

# Run a test and print PASS or *FAIL*.  If a test fails then increment
# the `nerrors' global variable and (if $verbose is set) display the
# failed output.  The actual output is not removed if $HDF5_NOCLEANUP is
# defined.
#
TOOLTEST() {
    out=test_$H5TOOL_$$.out
    err=test_$H5TOOL_$$.err

    # Run test.
    TESTING $H5TOOL $@
    $H5TOOL_BIN $@ > $out 2>&1
    result=$?
    if [ $result = 0 ]; then
	echo " PASSED"
	test $verbose -gt 1 && \
	    ( echo "========== results ==========="; cat $out;
	      echo "===============================================") |sed 's/^/    /'
    else
	echo "*FAILED*"
	nerrors="`expr $nerrors + 1`"
	test $verbose -gt 0 && \
	    ( echo "========== results ==========="; cat $out;
	      echo "===============================================") |sed 's/^/    /'
    fi

    # Clean up output file
    if test -z "$HDF5_NOCLEANUP"; then
	rm -f $out
    fi
}

# Print a "SKIP" message
SKIP() {
	 TESTING $H5TOOL $@
	  echo  " -SKIP-"
}


##############################################################################
###			  T H E   T E S T S                                ###
##############################################################################
#
# Group 1: HDF5 program that calls HDF5 APIs.
echo "***"Simple Compile and Link in one step.
TOOLTEST $hdf5main
# Application program that calls HDF5 and its own functions.
TOOLTEST $appmain $prog1 $prog2
# Repeat with -shlib option
echo "***"Simple Compile and Link with -shlib in one step.
TOOLTEST -shlib $hdf5main
# Application program that calls HDF5 and its own functions.
TOOLTEST -shlib $appmain $prog1 $prog2

# Group 2: Compile, then link.
echo "***"Compile and Link in two steps.
TOOLTEST -c $hdf5main
TOOLTEST $hdf5main_o
TOOLTEST -c $appmain $prog1 $prog2
TOOLTEST $appmain_o $prog1_o $prog2_o
# Repeat with -shlib option
echo "***"Compile and Link with -shlib in two steps.
TOOLTEST -c $hdf5main
TOOLTEST -shlib $hdf5main_o
TOOLTEST -c $appmain $prog1 $prog2
TOOLTEST -shlib $appmain_o $prog1_o $prog2_o

# Group3: Build external library, then link with it.
echo "***"Build external library and link with it.
TOOLTEST -c $prog1 $prog2
$AR cru $applib $prog1_o $prog2_o
$RANLIB $applib
TOOLTEST $appmain $applib
TOOLTEST $appmain_o $applib
# Repeat with -shlib option
echo "***"Build external library and link with it using -shlib.
TOOLTEST -c $prog1 $prog2
$AR cru $applib $prog1_o $prog2_o
$RANLIB $applib
TOOLTEST -shlib $appmain $applib
TOOLTEST -shlib $appmain_o $applib

# Group 4: Just preprocess, no compile, no link.
echo "***"Just preprocess, no compile, no link.
TOOLTEST -E $hdf5main
TOOLTEST -E $appmain $prog1 $prog2

##############################################################################
# END
##############################################################################

# Clean up  file
if test -z "$HDF5_NOCLEANUP"; then
    rm -f $temp_SRC $temp_OBJ $temp_FILES
fi

if test $nerrors -eq 0 ; then
   echo "All $H5TOOL tests passed."
fi

exit $nerrors
