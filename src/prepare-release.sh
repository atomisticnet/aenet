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

 NEW_VERSION is the version number of the upcoming release.

 The script does the following:
   1. Update the VERSION file.
   2. Update the license header in every source file.
"

if [[ $# -lt 1 ]] || [ "$1" == "-h" ] || [ "$1" == "--help" ]
then
    echo "${usage}"
    exit 0
fi

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

branch="$(git rev-parse --abbrev-ref HEAD)"
version="$1"
current_version="$(git describe --abbrev=0)"

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

echo " Tag the new release with    : git tag -a ${version}"
echo " Push the tags to origin with: git push --tags origin"
echo

exit 0
