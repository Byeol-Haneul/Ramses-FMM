#! /bin/bash
#######################################
# cr_write_makefile.sh [makefile name]
# creates a .f90 file with code to
# write the Makefile content to disk
#######################################

test -z "$1" && exit 0
MAKEFILE=$1
test -z "$2" && OUTFILE=write_makefile.f90 || OUTFILE=$2

cat << EOF >$OUTFILE
subroutine output_makefile(filename)
  character(LEN=80)::filename
  character(LEN=80)::fileloc
  character(LEN=30)::format
  integer::ilun

  ilun=11

  fileloc=TRIM(filename)
  format="(A)"
  open(unit=ilun,file=fileloc,form='formatted')
$(sed "s/\"/\"\"/g;s/^/  write(ilun,format)\"/;s/$/\"/" ${MAKEFILE})
  close(ilun)
end subroutine output_makefile
EOF
