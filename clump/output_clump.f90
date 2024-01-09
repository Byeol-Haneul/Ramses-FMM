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
  
  character(LEN=flen)::filename
  integer::rID

  if(pst%nLower>0)then
    rID = mdl_send_request(pst%s%mdl,MDL_OUTPUT_CLUMP,pst%iUpper+1,input_size,output_size,input_array)
    call r_output_clump(pst%pLower,input_array,input_size,output_array,output_size)
    call mdl_get_reply(pst%s%mdl,rID,output_size)
  else
    filename=transfer(input_array,filename)
    call output_clump_properties(pst%s,filename)
    call output_clump_field(pst%s,filename)
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
  !----------------------------------------------------------------
  ! This routine output the clump properties for each processor
  !----------------------------------------------------------------
  integer::i,idim,ilun,ierr,ilun2,j
  character(LEN=flen)::fileloc
  character(LEN=5)::nchar
  logical::file_exist

  associate(r=>s%r,g=>s%g,c=>s%c,p=>s%p)

  call title(g%myid,nchar)

  !-------------------------------------
  ! Open clump file and write header
  !-------------------------------------
  ilun=10
  fileloc=TRIM(filename)//'clump.'//TRIM(nchar)
  inquire(file=fileloc, exist=file_exist)
  if (file_exist) then
     open(unit=ilun,file=fileloc,iostat=ierr)
     close(ilun,status="delete")
  end if  
  open(unit=ilun,file=TRIM(fileloc),form='formatted')
  rewind(ilun)
  write(ilun,'(144A)')'   index  halo   lev   parent      ncell    peak_x             peak_y             peak_z     '//&
          '        rho-               rho+               rho_av             mass_cl            relevance   '
  !-------------------------------------
  ! Open halo file and write header
  !-------------------------------------
  if(r%saddle_threshold>0)then
     ilun2=20
     fileloc=TRIM(filename)//'halo.'//TRIM(nchar)
     inquire(file=fileloc, exist=file_exist)
     if (file_exist) then
        open(unit=ilun2,file=fileloc,iostat=ierr)
        close(ilun2,status="delete")
     end if
     open(unit=ilun2,file=TRIM(fileloc),form='formatted')
     rewind(ilun2)
     write(ilun2,'(135A)')'     index      ncell    peak_x             peak_y             peak_z     '//&
          '        rho+               mass      '
  endif

  do j=1,c%npeak
     if (c%relevance(j) > r%relevance_threshold .and. c%halo_mass(j) > r%mass_threshold*g%mp_min)then
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
     end if
     
     if(r%saddle_threshold>0)then
        if(c%ind_halo(j).EQ.j+c%npeak_cum(g%myid-1).AND.c%halo_mass(j) > r%mass_threshold*g%mp_min)then
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
  !----------------------------------------------------------------
  ! This routine output the peak patch fields for each processor
  !----------------------------------------------------------------  
  integer::ilevel,igrid,ilun,ierr,ivar,ind
  character(LEN=5)::nchar
  character(LEN=flen)::fileloc
  logical::file_exist
  
  associate(g=>s%g,r=>s%r,m=>s%m,mdl=>s%mdl)

  ilun=10+mdl_core(mdl)
  call title(g%myid,nchar)
  fileloc=TRIM(filename)//'clump_field.'//TRIM(nchar)
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
#ifdef GRAV
  do ilevel=r%levelmin,r%nlevelmax
     do igrid=m%head(ilevel),m%tail(ilevel)
        write(ilun)m%grid(igrid)%rho
        write(ilun)m%grid(igrid)%flag1
     end do
  enddo
#endif
  close(ilun)
  
end associate

end subroutine output_clump_field
!###################################################
!###################################################
!###################################################
!###################################################
subroutine file_descriptor_clump(r,filename)
  use amr_parameters, only: ndim,flen
  use hydro_parameters, only: nvar,nener
  use amr_commons, only: run_t
  implicit none
  type(run_t)::r
  character(LEN=flen)::filename
  
  character(LEN=flen)::fileloc
  integer::ivar,ilun

  if(r%verbose)write(*,*)'Entering file_descriptor_clump'

  ilun=11
  fileloc=TRIM(filename)
  open(unit=ilun,file=fileloc,form='formatted')

  write(ilun,'("nvar        =",I11)')2

  ivar=1
  write(ilun,'("variable #",I2,": density")')ivar

  ivar=2
  write(ilun,'("variable #",I2,": clump ID")')ivar
  
  close(ilun)

end subroutine file_descriptor_clump

end module output_clump_module
