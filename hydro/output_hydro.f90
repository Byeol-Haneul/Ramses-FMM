module output_hydro_module
contains
!###################################################
!###################################################
!###################################################
!###################################################
recursive subroutine r_output_hydro(pst,input_array,input_size,output_array,output_size)
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
     rID = mdl_send_request(pst%s%mdl,MDL_OUTPUT_HYDRO,pst%iUpper+1,input_size,output_size,input_array)
     call r_output_hydro(pst%pLower,input_array,input_size,output_array,output_size)
     call mdl_get_reply(pst%s%mdl,rID,output_size)
  else
     filename=transfer(input_array,filename)
     if(index(filename,'output')==0)then
        call backup_hydro(pst%s%r,pst%s%g,pst%s%m,pst%s%mdl,filename)
     else
        call output_hydro(pst%s%r,pst%s%g,pst%s%m,pst%s%mdl,filename)
     endif
  endif

end subroutine r_output_hydro
!###################################################
!###################################################
!###################################################
!###################################################
subroutine output_hydro(r,g,m,mdl,filename)
  use amr_parameters, only: ndim,twotondim,flen,dp
  use hydro_parameters, only: nvar
  use amr_commons, only: run_t,global_t,mesh_t
  use mdl_module
  implicit none
  type(run_t)::r
  type(global_t)::g
  type(mesh_t)::m
  type(mdl_t)::mdl
  character(LEN=flen)::filename

  real(dp),dimension(1:twotondim,1:nvar)::uold
  real(dp)::etot,ekin,dd,pp
  real(dp),dimension(1:ndim)::vv
  real(kind=4),dimension(1:twotondim,1:nvar)::uout
  integer::ilevel,igrid,ilun,ierr,ivar,ind
  character(LEN=5)::nchar
  character(LEN=flen)::fileloc
  logical::file_exist

#ifdef HYDRO

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
        uold=m%grid(igrid)%uold
        do ind=1,twotondim
           ! Compute density
           dd=uold(ind,1)
           ! Compute velocity
           vv(1:ndim)=uold(ind,2:ndim+1)/dd
           ! Compute kinetic energy
           ekin=0.
#if NDIM>0
           ekin=ekin+0.5*dd*vv(1)**2
#endif
#if NDIM>1
           ekin=ekin+0.5*dd*vv(2)**2
#endif
#if NDIM>2
           ekin=ekin+0.5*dd*vv(3)**2
#endif
           ! Compute pressure
           etot=uold(ind,ndim+2)
           pp=(r%gamma-1)*(etot-ekin)
           ! Store as primitive variables
           uold(ind,1)=dd
           uold(ind,2:ndim+1)=vv
           uold(ind,ndim+2)=pp
#if NVAR>NDIM+2+NENER
           ! Compute passive scalars
           do ivar=ndim+3,nvar
              uold(ind,ivar)=uold(ind,ivar)/dd
           end do
#endif
        end do
        uout=real(uold,kind=4)
        write(ilun)uout
     end do
  enddo
  close(ilun)

#endif
     
end subroutine output_hydro
!###################################################
!###################################################
!###################################################
!###################################################
subroutine backup_hydro(r,g,m,mdl,filename)
  use amr_parameters, only: ndim,flen
  use hydro_parameters, only: nvar
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

#ifdef HYDRO

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
        write(ilun)m%grid(igrid)%uold
     end do
  enddo
  close(ilun)

#endif

end subroutine backup_hydro
!###################################################
!###################################################
!###################################################
!###################################################
subroutine file_descriptor_hydro(r,filename,write_bkp_file)
  use amr_parameters, only: ndim,flen
  use hydro_parameters, only: nvar,nener
  use amr_commons, only: run_t
  implicit none
  type(run_t)::r
  character(LEN=flen)::filename
  logical::write_bkp_file
  
  character(LEN=flen)::fileloc
  integer::ivar,ilun

  if(r%verbose)write(*,*)'Entering file_descriptor_hydro'

  ilun=11

  ! Open file
  fileloc=TRIM(filename)
  open(unit=ilun,file=fileloc,form='formatted')

  if(write_bkp_file)then
     ! Write variable names in backup file
     write(ilun,'("nvar        =",I11)')nvar
     ivar=1
     write(ilun,'("variable #",I2,": density")')ivar
     ivar=2
     write(ilun,'("variable #",I2,": momentum_x")')ivar
     if(ndim>1)then
        ivar=3
        write(ilun,'("variable #",I2,": momentum_y")')ivar
     endif
     if(ndim>2)then
        ivar=4
        write(ilun,'("variable #",I2,": momentum_z")')ivar
     endif
#if NENER>0
     ! Non-thermal pressures
     do ivar=ndim+2,ndim+1+nener
        write(ilun,'("variable #",I2,": non_thermal_energy_",I1)')ivar,ivar-ndim-1
     end do
#endif
     ivar=ndim+2+nener
     write(ilun,'("variable #",I2,": total_energy")')ivar
#if NVAR>NDIM+2+NENER
     ! Passive scalars
     do ivar=ndim+3+nener,nvar
        write(ilun,'("variable #",I2,": density_scalar_",I1)')ivar,ivar-ndim-2-nener
     end do
#endif
  else
     ! Write variable names in output file
     write(ilun,'("nvar        =",I11)')nvar
     ivar=1
     write(ilun,'("variable #",I2,": density")')ivar
     ivar=2
     write(ilun,'("variable #",I2,": velocity_x")')ivar
     if(ndim>1)then
        ivar=3
        write(ilun,'("variable #",I2,": velocity_y")')ivar
     endif
     if(ndim>2)then
        ivar=4
        write(ilun,'("variable #",I2,": velocity_z")')ivar
     endif
#if NENER>0
     ! Non-thermal pressures
     do ivar=ndim+2,ndim+1+nener
        write(ilun,'("variable #",I2,": non_thermal_pressure_",I1)')ivar,ivar-ndim-1
     end do
#endif
     ivar=ndim+2+nener
     write(ilun,'("variable #",I2,": thermal_pressure")')ivar
#if NVAR>NDIM+2+NENER
     ! Passive scalars
     do ivar=ndim+3+nener,nvar
        write(ilun,'("variable #",I2,": scalar_",I1)')ivar,ivar-ndim-2-nener
     end do
#endif
  endif
  close(ilun)

end subroutine file_descriptor_hydro
end module output_hydro_module
