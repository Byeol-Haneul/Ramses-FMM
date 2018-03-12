#! /bin/bash

#######################################
# cr_write_patch.sh [patch dir]
# creates a .f90 file with code to
# write the patch dir content to disk
#######################################
PATCHDIRS="$1"
test -z "$2" && OUTFILE=write_patch.f90 || OUTFILE=$2

cat << EOF > $OUTFILE
subroutine output_patch(filename)
  character(LEN=80)::filename
  character(LEN=80)::fileloc
  character(LEN=30)::format
  integer::ilun

  ilun=11

  fileloc=TRIM(filename)
  format="(A)"
  open(unit=ilun,file=fileloc,form='formatted')
$(
  if test -z "$PATCHDIRS" ; then
    echo "no patches"
    echo " "
  else
    PATCHDIRS=(${PATCHDIRS//:/ })
    for PATCHDIR in "${PATCHDIRS[@]}"; do
	for filename in ${PATCHDIR}/*.f90; do
            echo "$filename" 
            cat "$filename"
	done 
    done
  fi |  sed 's/\$/ /g' |
    sed "s/\"/'/g" | cat -e | sed 's/\$/\"/' | sed 's/^/  write(ilun,format)"/'
)
  close(ilun)
end subroutine output_patch
EOF
