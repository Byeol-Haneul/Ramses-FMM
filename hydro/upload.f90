!#########################################################################
!#########################################################################
!#########################################################################
!#########################################################################
subroutine m_upload_fine(r,g,m,p,mdl,ilevel)
  use amr_commons, only: run_t,global_t,mesh_t
  use pm_commons, only: part_t
  use mdl_commons, only: mdl_t
  implicit none
  type(run_t)::r
  type(global_t)::g
  type(mesh_t)::m
  type(part_t)::p
  type(mdl_t)::mdl
  integer::ilevel
  !--------------------------------------------------------------------
  ! This routine is the master procedure to upload HYDRO variables
  ! from level ilevel+1 to ilevel (averaging down or restriction).
  !--------------------------------------------------------------------
  if(ilevel==r%nlevelmax)return
  if(m%noct_tot(ilevel)==0)return
  if(m%noct_tot(ilevel+1)==0)return
  if(r%verbose)write(*,111)ilevel
111 format(' Entering upload_fine for level',i2)

  call r_upload_fine(r,g,m,p,mdl,g%ncpu,1,0,ilevel)

end subroutine m_upload_fine
!################################################################
!################################################################
!################################################################
!################################################################
recursive subroutine r_upload_fine(r,g,m,p,mdl,cpu_range,input_size,output_size,ilevel)
  use amr_commons, only: run_t,global_t,mesh_t
  use pm_commons, only: part_t
  use mdl_commons, only: mdl_t
  use mdl_parameters
  implicit none
  type(run_t)::r
  type(global_t)::g
  type(mesh_t)::m
  type(part_t)::p
  type(mdl_t)::mdl
  integer::cpu_range,input_size,output_size
  integer::ilevel

  integer::next_range,next_cpu

  next_range=cpu_range/2
  next_cpu=g%myid+next_range

  if(next_range>0)then
     call mdl_send_request(mdl,MDL_UPLOAD_FINE,next_cpu,next_range,input_size,output_size,ilevel)
     call r_upload_fine(r,g,m,p,mdl,next_range,input_size,output_size,ilevel)
  else
     call upload_fine(r,g,m,ilevel)
  endif

end subroutine r_upload_fine
!###########################################################
!########################################################### 
!###########################################################
!###########################################################
subroutine upload_fine(r,g,m,ilevel)
  use amr_parameters, only: dp,ndim,twotondim
  use amr_commons, only: run_t,global_t,mesh_t
  use hydro_parameters, only: nvar,nener
  use cache_commons
  implicit none
  type(run_t)::r
  type(global_t)::g
  type(mesh_t)::m
  integer::ilevel
  !----------------------------------------------------------------------
  ! This routine performs a restriction operation (averaging down)
  ! for the hydro variables.
  !----------------------------------------------------------------------
#if NENER>0
  integer::irad
#endif
  integer::ioct,parent_cell,get_parent_cell
  integer::ind,ivar,igrid,icell,idim
  integer(kind=8),dimension(0:ndim)::hash_key
  real(dp)::average,ekin,erad

#ifdef HYDRO

  ! Set conservative variable to zero in refined cells
  do ioct=m%head(ilevel),m%tail(ilevel)
     do ivar=1,nvar
        do ind=1,twotondim
           if(m%grid(ioct)%refined(ind))then
              m%grid(ioct)%uold(ind,ivar)=0.0
           endif
        end do
     end do
  end do

  call open_cache(r,g,m,operation_upload,domain_decompos_amr)

  ! Loop over finer level grids
  hash_key(0)=ilevel+1
  do ioct=m%head(ilevel+1),m%tail(ilevel+1)

     ! Get cell and grid index
     hash_key(1:ndim)=m%grid(ioct)%ckey(1:ndim)
     parent_cell=get_parent_cell(r,g,m,hash_key,m%grid_dict,.true.,.false.)
     igrid=(parent_cell-1)/twotondim+1
     icell=parent_cell-(igrid-1)*twotondim

     ! Average conservative variables
     do ivar=1,nvar
        average=0.0d0
        do ind=1,twotondim
           average=average+m%grid(ioct)%uold(ind,ivar)
        end do
        ! Scatter result to cell
        m%grid(igrid)%uold(icell,ivar)=average/dble(twotondim)
     end do

     ! Average internal energy instead of total energy
     if(r%interpol_var==1 .or. r%interpol_var==2)then
        average=0.0d0
        do ind=1,twotondim
           ekin=0.0d0
           do idim=1,ndim
              ekin=ekin+0.5d0*m%grid(ioct)%uold(ind,idim+1)**2/max(m%grid(ioct)%uold(ind,1),r%smallr)
           end do
           erad=0.0d0
#if NENER>0
           do irad=1,nener
              erad=erad+m%grid(ioct)%uold(ind,ndim+2+irad)
           end do
#endif
           average=average+m%grid(ioct)%uold(ind,ndim+2)-ekin-erad
        end do
        ! Scatter result to cell
        ekin=0.0d0
        do idim=1,ndim
           ekin=ekin+0.5d0*m%grid(igrid)%uold(icell,idim+1)**2/max(m%grid(igrid)%uold(icell,1),r%smallr)
        end do
        erad=0.0d0
#if NENER>0
        do irad=1,nener
           erad=erad+m%grid(igrid)%uold(icell,ndim+2+irad)
        end do
#endif
        m%grid(igrid)%uold(icell,ndim+2)=average/dble(twotondim)+ekin+erad
     endif
  end do

  call close_cache(r,g,m,m%grid_dict)

#endif

end subroutine upload_fine
!##########################################################################
!##########################################################################
!##########################################################################
!##########################################################################
