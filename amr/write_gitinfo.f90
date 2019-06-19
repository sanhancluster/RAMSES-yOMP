subroutine write_gitinfo
  use amr_commons, ONLY:builddate,patchdir,gitrepo,gitbranch,githash
  implicit none

  builddate = TRIM(BUILDDATE)
  patchdir  = TRIM(PATCH)
  gitrepo   = TRIM(GITREPO)
  gitbranch = TRIM(GITBRANCH)
  githash   = TRIM(GITHASH)

  write(*,*)' '
  write(*,'(" compile date = ",A)')builddate
  write(*,'(" patch dir    = ",A)')patchdir
  write(*,'(" remote repo  = ",A)')gitrepo
  write(*,'(" local branch = ",A)')gitbranch
  write(*,'(" last commit  = ",A)')githash
  write(*,*)' '

end subroutine write_gitinfo
