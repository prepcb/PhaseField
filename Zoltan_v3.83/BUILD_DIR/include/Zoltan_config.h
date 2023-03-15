/* src/include/Zoltan_config.h.  Generated from Zoltan_config.h.in by configure.  */
/* src/include/Zoltan_config.h.in.  Generated from configure.ac by autoheader.  */

/* Define to dummy `main' function (if any) required to link to the Fortran
   libraries. */
/* #undef FC_DUMMY_MAIN */

/* Define if F77 and FC dummy `main' functions are identical. */
/* #undef FC_DUMMY_MAIN_EQ_F77 */

/* Define to a macro mangling the given C identifier (in lower and upper
   case), which must not contain underscores, for linking with Fortran. */
#define FC_FUNC(name,NAME) name ## _

/* As FC_FUNC, but for C identifiers containing underscores. */
#define FC_FUNC_(name,NAME) name ## _

/* Define if want to build examples */
#define HAVE_EXAMPLES /**/

/* Define if you want to build export makefiles. */
#define HAVE_EXPORT_MAKEFILES /**/

/* Define if want to build with f90interface enabled */
#define HAVE_F90INTERFACE 1

/* Define if you are using gnumake - this will shorten your link lines. */
/* #undef HAVE_GNUMAKE */

/* Define if want to build with gzip enabled */
/* #undef ZHAVE_GZIP */

/* define if we want to use MPI */
#define HAVE_MPI /**/

/* Define if want to build with nemesis_exodus enabled */
/* #undef HAVE_NEMESIS_EXODUS */

/* Define if want to build with parmetis enabled */
/* #undef HAVE_PARMETIS */

/* Define if want to build with patoh enabled */
/* #undef HAVE_PATOH */

/* Define if want to build with scotch enabled */
/* #undef HAVE_SCOTCH */

/* Define if want to build tests */
#define HAVE_TESTS /**/

/* Define if want to build zoltan-cppdriver */
#define HAVE_ZOLTAN_CPPDRIVER /**/

/* Define if want to build zoltan-examples */
#define HAVE_ZOLTAN_EXAMPLES /**/

/* Define if want to build with octreepartitioning enabled */
/* #undef HAVE_ZOLTAN_OCT */

/* Define if want to build zoltan-tests */
#define HAVE_ZOLTAN_TESTS /**/

/* software host will be cygwin */
/* #undef HOST_CYGWIN */

/* software host will be linux */
#define HOST_LINUX 1

/* software host will be solaris */
/* #undef HOST_SOLARIS */

/* Define to 1 if your C compiler doesn't accept -c and -o together. */
/* #undef NO_MINUS_C_MINUS_O */

/* define if ZOLTAN_ID_TYPE is unsigned int */
#define UNSIGNED_INT_GLOBAL_IDS 1

/* define if ZOLTAN_ID_TYPE is unsigned long */
/* #undef UNSIGNED_LONG_GLOBAL_IDS */

/* define if ZOLTAN_ID_TYPE is unsigned long long */
/* #undef UNSIGNED_LONG_LONG_GLOBAL_IDS */
