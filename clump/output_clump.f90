module output_clump_module
contains
!###################################################
!###################################################
!###################################################
!###################################################
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
    call output_clump_properties(pst%s,filename2)
    filename2=TRIM(filename)//'clump_field.'
    call output_clump_field(pst%s,filename2)
  endif

end subroutine r_output_clump
!###################################################
!###################################################
!###################################################
!###################################################
subroutine output_clump_properties(s,filename)
  use amr_parameters, only: ndim,dp,i8b,flen
  use hydro_parameters, only: nvar
  use ramses_commons, only: ramses_t
  use clfind_commons
  implicit none
  type(ramses_t)::s
  character(LEN=flen)::filename
  real(dp)::rel_mass,rel_mass_tot
  real(dp)::particle_mass=0
  integer::i,idim,ilun,ierr,ilun2,j,n_rel_tot,n_rel
  character(LEN=flen)::fileloc
  character(LEN=5)::nchar
  logical::file_exist

#ifndef WITHOUTMPI
  real(dp)::particle_mass_tot
#endif

  associate(r=>s%r,g=>s%g,c=>s%c,p=>s%p)

  particle_mass=MINVAL(p%mp, MASK=(p%mp > 0))
#ifndef WITHOUTMPI
  call MPI_ALLREDUCE(particle_mass,particle_mass_tot,1,MPI_DOUBLE_PRECISION,MPI_MIN,MPI_COMM_WORLD,info)
  particle_mass=particle_mass_tot
#endif
  ilun=10
  call title(g%myid,nchar)
  fileloc=TRIM(filename)//TRIM(nchar)
  inquire(file=fileloc, exist=file_exist)
  if (file_exist) then
    open(unit=ilun,file=fileloc,iostat=ierr)
    close(ilun,status="delete")
  end if

  rel_mass=0
  n_rel=0

  open(unit=ilun,file=TRIM(fileloc),access="stream",action="write",form='unformatted')
  rewind(ilun)
  ! Write header
  write(ilun)ndim
  write(ilun)c%npeak

  do j=1,c%npeak
    if (c%relevance(j) > r%relevance_threshold .and. c%halo_mass(j) > r%mass_threshold*particle_mass)then
      write(ilun,'(I8,X,I2,X,I10,X,I10,8(X,1PE18.9E2))')&
        j+c%npeak_cum(g%myid-1)&
        ,c%lev_peak(j)&
        ,c%new_peak(j)&
        ,c%n_cells(j)&
        ,c%peak_pos(j,1)&
        ,c%peak_pos(j,2)&
        ,c%peak_pos(j,3)&
        ,c%min_dens(j)&
        ,c%max_dens(j)&
        ,c%clump_mass(j)/c%clump_vol(j)&
        ,c%clump_mass(j)&
        ,c%relevance(j)
      rel_mass=rel_mass+c%clump_mass(j)
      n_rel=n_rel+1
    end if

    if(r%saddle_threshold>0)then
      if(c%ind_halo(j).EQ.j+c%npeak_cum(g%myid-1).AND.c%halo_mass(j) > r%mass_threshold*particle_mass)then
         write(ilun2,'(I10,X,I10,5(X,1PE18.9E2))')&
          j+c%npeak_cum(g%myid-1)&
          ,c%n_cells_halo(j)&
          ,c%peak_pos(j,1)&
          ,c%peak_pos(j,2)&
          ,c%peak_pos(j,3)&
          ,c%max_dens(j)&
          ,c%halo_mass(j)
      endif
   endif
  enddo
  close(ilun)
  if(r%saddle_threshold>0)then
      close(ilun2)
  endif

#ifndef WITHOUTMPI
  call MPI_ALLREDUCE(n_rel,n_rel_tot,1,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,info)
  n_rel=n_rel_tot
  call MPI_ALLREDUCE(rel_mass,rel_mass_tot,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
  rel_mass=rel_mass_tot
#else
  n_rel_tot = n_rel
  rel_mass_tot = rel_mass
#endif
  
  end associate
end subroutine output_clump_properties
!###################################################
!###################################################
!###################################################
!###################################################
subroutine output_clump_field(s,filename)
  use amr_parameters, only: ndim,twotondim,flen,dp
  use ramses_commons, only: ramses_t
  use mdl_module
  
  implicit none
  type(ramses_t)::s
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
  write(ilun)2
  write(ilun)r%levelmin
  write(ilun)r%nlevelmax
  do ilevel=r%levelmin,r%nlevelmax
     write(ilun)m%noct(ilevel)
  enddo
  do ilevel=r%levelmin,r%nlevelmax
     do igrid=m%head(ilevel),m%tail(ilevel)
        write(ilun)m%grid(igrid)%rho
        write(ilun)m%grid(igrid)%flag1
     end do
  enddo
  close(ilun)
  
end associate
end subroutine output_clump_field
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
  write(ilun,'("variable #",I2,": halo_mass")')ivar

  ivar=2*ndim+13
  write(ilun,'("variable #",I2,": clump_mass")')ivar

  ivar=2*ndim+14
  write(ilun,'("variable #",I2,": relevance")')ivar
  
  close(ilun)

end subroutine file_descriptor_clump

end module output_clump_module
