#!/bin/bash
GITBRANCH=$(git rev-parse --abbrev-ref HEAD)
GITHASH=$(git log --pretty=format:'%H' -n 1)
GITREPO=$(git config --get remote.origin.url)
BUILDDATE=$(date +"%D-%T")
cat <<EOF >$1
subroutine write_gitinfo
  use amr_parameters, ONLY:builddate,patchdir,gitrepo,gitbranch,githash

  builddate = "$BUILDDATE"
  patchdir  = "$PATCH"
  gitrepo   = "$GITREPO"
  gitbranch = "$GITBRANCH"
  githash   = "$GITHASH"

  write(*,*)' '
  write(*,'(" compile date = ",A)')TRIM(builddate)
  write(*,'(" patch dir    = ",A)')TRIM(patchdir)
  write(*,'(" remote repo  = ",A)')TRIM(gitrepo)
  write(*,'(" local branch = ",A)')TRIM(gitbranch)
  write(*,'(" last commit  = ",A)')TRIM(githash)
  write(*,*)' '

end subroutine write_gitinfo
EOF
