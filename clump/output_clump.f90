module output_clump_module
contains
!###################################################
!###################################################
!###################################################
!###################################################
recursive subroutine r_output_clump_field(pst,input_array,input_size,output_array,output_size)
  use mdl_module
  use amr_parameters, only: flen
  use ramses_commons, only: pst_t
  use mdl_parameters
  implicit none
  type(pst_t)::pst
  integer,VALUE::input_size
  integer::output_size
  integer,dimension(1:input_size)::input_array
  integer,dimension(1:output_size)::output_array
  
  character(LEN=flen)::filename
  integer::rID

  if(pst%nLower>0)then
     rID = mdl_send_request(pst%s%mdl,MDL_OUTPUT_CLUMP,pst%iUpper+1,input_size,output_size,input_array)
     call r_output_clump_field(pst%pLower,input_array,input_size,output_array,output_size)
     call mdl_get_reply(pst%s%mdl,rID,output_size)
  else
     filename=transfer(input_array,filename)
      call output_clump_field(s,filename)
  endif

end subroutine r_output_clump_field
!#######################################################
!#######################################################
!#######################################################
!#######################################################
subroutine output_clump_field(s,filename)
  use amr_parameters, only: ndim,twotondim,flen,dp
  use hydro_parameters, only: nvar
  use ramses_commons, only: ramses_t
  use mdl_module

  implicit none
  type(rasmes_t)::s
  character(LEN=flen)::filename

  integer::ilevel,igrid,ilun,ierr,ivar,ind
  character(LEN=5)::nchar
  character(LEN=flen)::fileloc
  logical::file_exist

  associate(g=>s%g,r=>s%r,m=>s%m,mdl=>s%mdl)
  ilun=10+mdl_core(mdl)
  call title(g%myid,nchar)
  fileloc=TRIM(filename)//TRIM(nchar)
  inquire(file=fileloc,exist=file_exist)
  if (file_exist) then
     open(unit=ilun,file=fileloc,iostat=ierr)
     close(ilun,status="delete")
  end if
  open(unit=ilun,file=fileloc,access="stream",action="write",form='unformatted')
  write(ilun)ndim
  write(ilun)nvar
  write(ilun)r%levelmin
  write(ilun)r%nlevelmax
  do ilevel=r%levelmin,r%nlevelmax
     write(ilun)m%noct(ilevel)
  enddo
  do ilevel=r%levelmin,r%nlevelmax
     do igrid=m%head(ilevel),m%tail(ilevel)
        write(ilun)m%grid(igrid)%flag1
     end do
  enddo
  close(ilun)

  end associate
     
end subroutine output_clump_field
!#######################################################
!#######################################################
!#######################################################
!#######################################################
recursive subroutine r_output_clump(pst,input_array,input_size,output_array,output_size)
  use mdl_module
  use amr_parameters, only: flen
  use ramses_commons, only: pst_t
  use mdl_parameters
  implicit none
  type(pst_t)::pst
  integer,VALUE::input_size
  integer::output_size
  integer,dimension(1:input_size)::input_array
  integer,dimension(1:output_size)::output_array

  character(LEN=flen)::filename,filename2
  integer::rID
  if(pst%nLower>0)then
    rID = mdl_send_request(pst%s%mdl,MDL_OUTPUT_CLUMP,pst%iUpper+1,input_size,output_size,input_array)
    call r_output_clump(pst%pLower,input_array,input_size,output_array,output_size)
    call mdl_get_reply(pst%s%mdl,rID,output_size)
  else
    filename=transfer(input_array,filename)
    filename2=TRIM(filename)//'clump.'
    call output_clump(s,filename2)
  endif

end subroutine r_output_clump


subroutine output_clump(s,filename)
  use amr_parameters, only: ndim,dp,i8b,flen
  use hydro_parameters, only: nvar
  use ramses_commons, only: ramses_t
  implicit none
  type(ramse_t)::s
  character(LEN=flen)::filename

  integer::i,idim,ilun,ierr
  character(LEN=flen)::fileloc
  character(LEN=5)::nchar
  real(kind=4),allocatable,dimension(:)::xsp
  integer(i8b),allocatable,dimension(:)::ii8
  integer,allocatable,dimension(:)::ll
  logical::file_exist

  associate(g=>s%g,r=>s%r,c=>s%c)

  ilun=10
  call title(g%myid,nchar)
  fileloc=TRIM(filename)//TRIM(nchar)
  inquire(file=fileloc, exist=file_exist)
  if (file_exist) then
    open(unit=ilun,file=fileloc,iostat=ierr)
    close(ilun,status="delete")
  end if
  open(unit=ilun,file=TRIM(fileloc),access="stream",action="write",form='unformatted')
  rewind(ilun)
  ! Write header
  write(ilun)ndim
  write(ilun)c%npeak


  ! Write xp,denp,levelp,idp,idc,pid,sortp,alive,peak_pos,lev_peak,ind_halo,n_cell_halo,n_cells,min_dens,max_dens,clump_vol,center_of_mass,clump_velocity,halo_mass,clump_mass,
  allocate(xsp(1:c%peak))
  ! Write peak_grid
  do i=1,c%npeak
    xsp(i)=c%peak_grid(i)
  end do
  write(ilun)xsp
  ! Write peak_cell
  do i=1,c%npeak
    xsp(i)=c%peak_cell(i)
  end do
  write(ilun)xsp
  ! Write max_dens
  do i=1,c%npeak
    xsp(i)=c%max_dens(i)
  end do
  write(ilun)xsp
  ! Write lev_peak
  allocate(ll(1:c%npeak))
  do i=1,c%npeak
    ll(i)=c%lev_peak(i)
  end do
  write(ilun)ll
  ! Write new_peak
  allocate(ii8(1:c%npeak))
  do i=1,c%npeak
    ii8(i)=c%new_peak(i)
  end do
  write(ilun)ii8
  ! Write peak_pos
  do idim=1,ndim
    do i=1,c%npeak
       xsp(i)=c%peak_pos(i,idim)
    end do
    write(ilun)xsp
  end do
  ! Write ind_halo
  do i=1,c%npeak
    ll(i)=c%ind_halo(i)
  end do
  write(ilun)ll
  ! Write n_cell_halo
  do i=1,c%npeak
    ll(i)=c%n_cell_halo(i)
  end do
  write(ilun)ll
  ! Write n_cells
  do i=1,c%npeak
    ll(i)=c%n_cells(i)
  end do
  write(ilun)ll
  ! Write min_dens
  do i=1,c%npeak
    xsp(i)=c%min_dens(i)
  end do
  write(ilun)xsp
  ! Write max_dens
  do i=1,c%npeak
    xsp(i)=c%av_dens(i)
  end do
  write(ilun)xsp
  ! Write clump_vol
  do i=1,c%npeak
    xsp(i)=c%clump_vol(i)
  end do
  write(ilun)xsp
  ! Write clump_velocity
  do idim=1,ndim
    do i=1,c%npeak
       xsp(i)=c%clump_velocity(i,idim)
    end do
    write(ilun)xsp
  end do
  ! Write center_of_mass
  do idim=1,ndim
    do i=1,c%npeak
       xsp(i)=c%center_of_mass(i,idim)
    end do
    write(ilun)xsp
  end do
  ! Write halo_mass
  do i=1,c%npeak
    xsp(i)=c%halo_mass(i)
  end do
  write(ilun)xsp
  ! Write clump_mass
  do i=1,c%npeak
    xsp(i)=c%clump_mass(i)
  end do
  write(ilun)xsp
  ! Write relevance
  do i=1,c%npeak
    xsp(i)=c%relevance(i)
  end do
  write(ilun)xsp

  deallocate(ll)
  deallocate(ii8)
  deallocate(xsp)
  close(ilun)

  end associate
end subroutine output_clump

!###################################################
!###################################################
!###################################################
!###################################################
subroutine file_descriptor_clump(r,filename,write_bkp_file)
  use amr_parameters, only: ndim,flen
  use hydro_parameters, only: nvar,nener
  use amr_commons, only: run_t
  implicit none
  type(run_t)::r
  character(LEN=flen)::filename
  logical::write_bkp_file
  
  character(LEN=flen)::fileloc
  integer::ivar,ilun

  if(r%verbose)write(*,*)'Entering file_descriptor_clump'

  ilun=11

  ! Open file
  fileloc=TRIM(filename)
  open(unit=ilun,file=fileloc,form='formatted')

  ! Write variable names in backup file
  ! Write peak_grid,peak_cell,max_dens,lev_peak,new_peak,peak_pos,ind_halo,n_cell_halo,n_cells,min_dens,av_dens,clump_vol,center_of_mass,clump_velocity,halo_mass,clump_mass,relevance
  write(ilun,'("nvar        =",I11)')nvar
  ivar=1
  write(ilun,'("variable #",I2,": peak_grid")')ivar

  ivar=2
  write(ilun,'("variable #",I2,": peak_cell")')ivar

  ivar=3
  write(ilun,'("variable #",I2,": max_dens")')ivar

  ivar=4
  write(ilun,'("variable #",I2,": lev_peak")')ivar

  ivar=5
  write(ilun,'("variable #",I2,": new_peak")')ivar

  ivar=6
  write(ilun,'("variable #",I2,": peak_pos_x")')ivar
  if(ndim>1)then
    ivar=7
    write(ilun,'("variable #",I2,": peak_pos_y")')ivar
  endif
  if(ndim>2)then
    ivar=8
    write(ilun,'("variable #",I2,": peak_pos_z")')ivar
  endif

  ivar=ndim+6
  write(ilun,'("variable #",I2,": ind_halo")')ivar

  ivar=ndim+7
  write(ilun,'("variable #",I2,": n_cell_halo")')ivar

  ivar=ndim+8
  write(ilun,'("variable #",I2,": n_cells")')ivar

  ivar=ndim+9
  write(ilun,'("variable #",I2,": min_dens")')ivar

  ivar=ndim+10
  write(ilun,'("variable #",I2,": av_dens")')ivar

  !clump_vol,center_of_mass,clump_velocity,halo_mass,clump_mass,relevance
  ivar=ndim+11
  write(ilun,'("variable #",I2,": clump_vol")')ivar

  ivar=ndim+12
  write(ilun,'("variable #",I2,": center_of_mass_x")')ivar
  if(ndim>1)then
    ivar=ndim+13
    write(ilun,'("variable #",I2,": center_of_mass_y")')ivar
  endif
  if(ndim>2)then
    ivar=ndim+14
    write(ilun,'("variable #",I2,": center_of_mass_z")')ivar
  endif

  ivar=2*ndim+12
  write(ilun,'("variable #",I2,": clump_velocity_x")')ivar
  if(ndim>1)then
    ivar=2*ndim+13
    write(ilun,'("variable #",I2,": clump_velocity_y")')ivar
  endif
  if(ndim>2)then
    ivar=2*ndim+14
    write(ilun,'("variable #",I2,": clump_velocity_z")')ivar
  endif

  ivar=3*ndim+13
  write(ilun,'("variable #",I2,": halo_mass")')ivar

  ivar=3*ndim+14
  write(ilun,'("variable #",I2,": clump_mass")')ivar

  ivar=3*ndim+15
  write(ilun,'("variable #",I2,": relevance")')ivar
  
  close(ilun)

end subroutine file_descriptor_clump

end module output_clump_module
