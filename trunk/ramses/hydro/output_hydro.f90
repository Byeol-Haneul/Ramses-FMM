subroutine backup_hydro(filename)
  use amr_commons
  use hydro_commons
  implicit none
  character(LEN=80)::filename

  integer::ilevel,igrid,ilun
  character(LEN=5)::nchar
  character(LEN=80)::fileloc

  if(verbose)write(*,*)'Entering backup_hydro'
  ilun=ncpu+myid+10     
  call title(myid,nchar)
  fileloc=TRIM(filename)//TRIM(nchar)
  open(unit=ilun,file=fileloc,form='unformatted')
  write(ilun)ndim
  write(ilun)levelmin
  write(ilun)nlevelmax
  write(ilun)nvar
  do ilevel=levelmin,nlevelmax
     write(ilun)noct(ilevel)
     do igrid=head(ilevel),tail(ilevel)
        write(ilun)grid(igrid)%uold
     end do
  enddo
  close(ilun)
     
end subroutine backup_hydro





