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
  
  character(LEN=flen)::filename,fileloc
  integer::rID

  if(pst%nLower>0)then
    rID = mdl_send_request(pst%s%mdl,MDL_OUTPUT_CLUMP,pst%iUpper+1,input_size,output_size,input_array)
    call r_output_clump(pst%pLower,input_array,input_size,output_array,output_size)
    call mdl_get_reply(pst%s%mdl,rID,output_size)
  else
    filename=transfer(input_array,filename)
    if(pst%s%r%output_clump)then
       call output_clump_properties(pst%s,filename)
    endif
    if(pst%s%r%output_clump_field)then
       fileloc=TRIM(filename)//'peak.'
       call output_clump_field(pst%s,fileloc)
    endif
  endif

end subroutine r_output_clump
!###################################################
!###################################################
!###################################################
!###################################################
subroutine output_clump_properties(s,filename)
  use amr_parameters, only: flen
  use ramses_commons, only: ramses_t,open_file,close_file
  use clfind_commons
  implicit none
  type(ramses_t)::s
  character(LEN=flen)::filename
  !----------------------------------------------------------------
  ! This routine output the clump properties for each processor
  !----------------------------------------------------------------
  integer::ilun,j
  character(LEN=flen)::fileloc
  integer(kind=8),dimension(s%r%levelmin:s%r%nlevelmax)::nskip

  associate(r=>s%r,g=>s%g,c=>s%c)

  ! Write clump file
  fileloc=TRIM(filename)//'clump.'

  call open_file(s,fileloc,nskip,ilun)

  do j=1,c%npeak
     if (c%relevance(j) > r%relevance_threshold.AND. &
          & c%clump_mass(j) > r%mass_threshold*g%mp_min.AND. &
          & c%halo_mass(j) > r%mass_threshold*g%mp_min)then
        write(ilun,'(I10,X,I10,1X,I2,X,I10,X,I10,8(X,1PE18.9E2))')&
             j+c%npeak_cum(g%myid-1)&
             ,c%ind_halo(j)&
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
  end do

  call close_file(s,fileloc,nskip,ilun)

  ! Write halo file
  if(r%saddle_threshold>0)then

     fileloc=TRIM(filename)//'halo.'

     call open_file(s,fileloc,nskip,ilun)

     do j=1,c%npeak
        if(c%ind_halo(j).EQ.j+c%npeak_cum(g%myid-1).AND.&
             & c%halo_mass(j) > r%mass_threshold*g%mp_min.AND. &
             & c%relevance(j) > r%relevance_threshold)then
           write(ilun,'(I10,X,I10,5(X,1PE18.9E2))')&
                j+c%npeak_cum(g%myid-1)&
                ,c%n_cells_halo(j)&
                ,c%peak_pos(j,1)&
                ,c%peak_pos(j,2)&
                ,c%peak_pos(j,3)&
                ,c%max_dens(j)&
                ,c%halo_mass(j)
        endif
     enddo

     call close_file(s,fileloc,nskip,ilun)

  end if
  
  end associate

end subroutine output_clump_properties
!###################################################
!###################################################
!###################################################
!###################################################
subroutine output_clump_field(s,filename)
  use amr_parameters, only: ndim,twotondim,flen
  use ramses_commons, only: ramses_t,open_file,close_file
  implicit none
  type(ramses_t)::s
  character(LEN=flen)::filename
  !----------------------------------------------------------------
  ! This routine output the peak patch fields for each processor
  !----------------------------------------------------------------  
  integer::ilevel,igrid,ilun
  integer(kind=8),dimension(s%r%levelmin:s%r%nlevelmax)::nskip
  real(kind=4),dimension(1:twotondim)::flg
  real(kind=4),dimension(1:twotondim)::rho

  associate(g=>s%g,r=>s%r,m=>s%m,mdl=>s%mdl)

#ifdef GRAV

  call open_file(s,filename,nskip,ilun)

  do ilevel=r%levelmin,r%nlevelmax
     write(ilun,POS=nskip(ilevel))
     do igrid=m%head(ilevel),m%tail(ilevel)
        rho=real(m%grid(igrid)%rho,kind=4)
        flg=real(m%grid(igrid)%flag1,kind=4)
        write(ilun)flg
        write(ilun)rho
     end do
  enddo

  call close_file(s,filename,nskip,ilun)

#endif

  end associate

end subroutine output_clump_field
!###################################################
!###################################################
!###################################################
!###################################################
subroutine file_descriptor_clump(r,filename)
  use amr_parameters, only: ndim,flen
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
  write(ilun,'("variable #",I2,": peak ID")')ivar

  ivar=2
  write(ilun,'("variable #",I2,": density")')ivar
  
  close(ilun)

end subroutine file_descriptor_clump

end module output_clump_module
