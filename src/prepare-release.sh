#!/bin/bash
#-----------------------------------------------------------------------
#                Prepare source tree for a new release
#-----------------------------------------------------------------------
# 2015-10-14 Nongnuch Artrith (NA) and Alexander Urban (AU)
#-----------------------------------------------------------------------
usage="
 Prepare source tree for a new release.

 Usage:
    $0 NEW_VERSION

 NEW_VERSION is MAJOR.MINOR.PATCH, without a v prefix.
 Run this script from src/. Release tags use vMAJOR.MINOR.PATCH.

 The script does the following:
   1. Update the VERSION file.
   2. Update the license header in every source file.
"

if [[ $# == 1 && ( "$1" == "-h" || "$1" == "--help" ) ]]; then
    echo "${usage}"
    exit 0
fi

# Validate before writing VERSION or any source headers.
version_pattern='^(0|[1-9][0-9]*)\.(0|[1-9][0-9]*)\.(0|[1-9][0-9]*)$'
if [[ $# != 1 || ! "$1" =~ ${version_pattern} ]]; then
    echo "Error: expected MAJOR.MINOR.PATCH without a v prefix." >&2
    exit 1
fi
if [[ ! -f VERSION || ! -f license-header.txt ]]; then
    echo "Error: run prepare-release.sh from src/." >&2
    exit 1
fi
set -e

#-------------- make sure all required tools are present --------------#

for tool in sed git awk fold find
do
  if [ "$(which ${tool})" == "" ]
  then
      echo " Error: required '${tool}' command not available. Aborting."
      exit 1
  fi
done

#------------------------ collect information -------------------------#

version="$1"
current_version="$(cat VERSION)"

echo
echo " Preparing bump from version ${current_version} to version $1."
echo

#------------------------ update VERSION file -------------------------#

echo "${version}" > ./VERSION

#----------------------- update license headers -----------------------#

# Fortran files
header="$(fold -w 69 -s license-header.txt | \
  awk '{s=sprintf("!+ %s", $0); sub(/ *$/, "", s); printf("%s\\n", s);}')"
for f in *.f90 *.F90 ./tests/*.f90 ./tools/*.f90 ./ext/*.f90
do
  awk '
    BEGIN { header = 0 }
    /^!\+/{
      if (header == 0) {
        header = 1;
        printf("'"${header}"'");
      };
      next
    }
    { print }
  ' $f > $f-tmp && mv $f-tmp $f
done

# Makefiles
header="$(fold -w 69 -s license-header.txt | \
  awk '{s=sprintf("#+ %s", $0); sub(/ *$/, "", s); printf("%s\\n", s);}')"
for f in Makefile $(find . -name "Makefile.inc") $(ls makefiles/Makefile.*)
do
  awk '
    BEGIN { header = 0 }
    /^#\+/{
      if (header == 0) {
        header = 1;
        printf("'"${header}"'");
      };
      next
    }
    { print }
  ' $f > $f-tmp && mv $f-tmp $f
done

#----------------------------------------------------------------------#

echo " Tag the new release with    : git tag -a v${version}"
echo " Push the tags to origin with: git push origin v${version}"
echo

exit 0
