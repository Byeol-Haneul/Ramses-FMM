#ifdef RT
module output_rt_module
contains
!###################################################
!###################################################
!###################################################
!###################################################
recursive subroutine r_output_rt(pst,input_array,input_size,output_array,output_size)
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
     rID = mdl_send_request(pst%s%mdl,MDL_OUTPUT_RT,pst%iUpper+1,input_size,output_size,input_array)
     call r_output_rt(pst%pLower,input_array,input_size,output_array,output_size)
     call mdl_get_reply(pst%s%mdl,rID,output_size)
  else
     filename=transfer(input_array,filename)
     if(index(filename,'output')==0)then
        call backup_rt(pst%s%r,pst%s%g,pst%s%m,pst%s%mdl,filename)
     else
        call output_rt(pst%s,filename)
     endif
  endif

end subroutine r_output_rt
!###################################################
!###################################################
!###################################################
!###################################################
subroutine output_rt(s,filename)
  use amr_parameters, only: ndim,twotondim,flen,dp
  use rt_parameters, only: nrtvar, nrtgroups, rt_c
  use ramses_commons, only: ramses_t,open_file,close_file
  implicit none
  type(ramses_t)::s
  character(LEN=flen)::filename
  !-----------------------------------
  ! Output rt data in file
  !-----------------------------------
  integer::ilevel,igrid,ilun,igrp,idim,ind
  integer(kind=8),dimension(s%r%levelmin:s%r%nlevelmax)::nskip
  real(kind=4),dimension(1:twotondim,1:nrtvar)::qout
  real(dp),    dimension(1:twotondim,1:nrtvar)::qold
  real(dp),    dimension(1:twotondim,1:nrtvar)::rtuold

  associate(r=>s%r,g=>s%g,m=>s%m,mdl=>s%mdl)
  call open_file(s,filename,nskip,ilun)
  do ilevel=r%levelmin,r%nlevelmax

     write(ilun,POS=nskip(ilevel))

     do igrid=m%head(ilevel),m%tail(ilevel)

        rtuold=m%grid(igrid)%rtuold
        do ind=1,twotondim
            do igrp = 1, nrtgroups
              qold(ind,1+(igrp-1)*(ndim+1)) = rtuold(ind,1+(igrp-1)*(ndim+1)) * rt_c
              do idim = 1, ndim
                qold(ind,1+idim+(igrp-1)*(ndim+1)) = rtuold(ind,1+idim+(igrp-1)*(ndim+1))
              end do
            end do
        end do
        qout=real(qold,kind=4)
        write(ilun)qout

     end do
  end do

  call close_file(s,filename,nskip,ilun)

  end associate

end subroutine output_rt
!###################################################
!###################################################
!###################################################
!###################################################
subroutine backup_rt(r,g,m,mdl,filename)
  use amr_parameters, only: ndim,flen
  use rt_parameters, only: nrtvar
  use amr_commons, only: run_t,global_t,mesh_t
  use mdl_module
  implicit none
  type(run_t)::r
  type(global_t)::g
  type(mesh_t)::m
  type(mdl_t)::mdl
  character(LEN=flen)::filename

  integer::ilevel,igrid,ilun,ierr
  character(LEN=5)::nchar
  character(LEN=flen)::fileloc
  logical::file_exist


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
  write(ilun)nrtvar
  write(ilun)r%levelmin
  write(ilun)r%nlevelmax
  do ilevel=r%levelmin,r%nlevelmax
     write(ilun)m%noct(ilevel)
  enddo
  do ilevel=r%levelmin,r%nlevelmax
     do igrid=m%head(ilevel),m%tail(ilevel)
        write(ilun)m%grid(igrid)%uold
     end do
  enddo
  close(ilun)


end subroutine backup_rt
!###################################################
!###################################################
!###################################################
!###################################################
subroutine file_descriptor_rt(r,filename,write_bkp_file)
  use amr_parameters, only: ndim,flen
  use rt_parameters, only: nrtgroups
  use amr_commons, only: run_t
  implicit none
  type(run_t)::r
  character(LEN=flen)::filename
  logical::write_bkp_file
  character(len=1), dimension(1:3), parameter :: dim_keys = ["x", "y", "z"]
  
  character(LEN=flen)::fileloc
  integer::igrp,idim,ilun

  if(r%verbose)write(*,*)'Entering file_descriptor_rt'

  ilun=11

  ! Open file
  fileloc=TRIM(filename)
  open(unit=ilun,file=fileloc,form='formatted')

  ! Write variable names in backup file
  write(ilun,'("nvar        =",I11)')nrtgroups*(1+ndim)
  do igrp = 1, nrtgroups
     write(ilun,'("variable #",I2,": photon_density_", i0.2)')1+(igrp-1)*(ndim+1), igrp
     do idim = 1, ndim
        write(ilun,'("variable #",I2,": photon_flux_", i0.2, "_", a)')1+idim+(igrp-1)*(ndim+1), igrp, dim_keys(idim)
     end do
  end do
  close(ilun)

end subroutine file_descriptor_rt
end module output_rt_module
#endif