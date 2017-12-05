!################################################################
!################################################################
!################################################################
!################################################################
subroutine m_synchro_hydro_fine(r,g,m,p,mdl,ilevel,dteff)
  use amr_parameters, only: ndim,dp,twotondim
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
  real(dp)::dteff
  !--------------------------------------------------------------
  ! Add gravity source terms to uold with time step dteff.
  !--------------------------------------------------------------
  integer,dimension(1:3)::input_array
  
  if(.not. r%poisson)return
  if(m%noct_tot(ilevel)==0)return
  if(r%verbose)write(*,'("   Entering synchro_hydro_fine for level",i2," and time step dt=",1PE12.5)')ilevel,dteff

  input_array(1)=ilevel
  input_array(2:3)=transfer(dteff,input_array)
  call r_synchro_hydro_fine(r,g,m,p,mdl,mdl%ncpu,3,0,input_array)
  
end subroutine m_synchro_hydro_fine
!################################################################
!################################################################
!################################################################
!################################################################
recursive subroutine r_synchro_hydro_fine(r,g,m,p,mdl,cpu_range,input_size,output_size,input_array)
  use amr_parameters, only: dp
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
  integer,dimension(1:input_size)::input_array

  integer::next_range,next_cpu
  integer::ilevel
  real(dp)::dteff
  
  next_range=cpu_range/2
  next_cpu=g%myid+next_range

  if(next_range>0)then
     call mdl_send_request(mdl,MDL_SYNCHRO_HYDRO_FINE,next_cpu,next_range,input_size,output_size,input_array)
     call r_synchro_hydro_fine(r,g,m,p,mdl,next_range,input_size,output_size,input_array)
  else
     ilevel=input_array(1)
     dteff=transfer(input_array(2:3),dteff)
     call synchro_hydro_fine(r,m,ilevel,dteff)
  endif

end subroutine r_synchro_hydro_fine
!################################################################
!################################################################
!################################################################
!################################################################
subroutine synchro_hydro_fine(r,m,ilevel,dteff)
  use amr_parameters, only: ndim,dp,twotondim
  use amr_commons, only: run_t,mesh_t
  implicit none
  type(run_t)::r
  type(mesh_t)::m
  integer::ilevel
  real(dp)::dteff
  !--------------------------------------------------------------
  ! Add gravity source terms to uold with time step dteff.
  !--------------------------------------------------------------
  integer::igrid,ind
  integer::idim,neul=ndim+2
  real(dp)::ener

#ifdef HYDRO

  ! Loop over octs
  do igrid=m%head(ilevel),m%tail(ilevel)
     ! Loop over cells
     do ind=1,twotondim

        ! Remove kinetic energy from total energy
        ener=m%grid(igrid)%uold(ind,neul)
        do idim=1,ndim
           ener=ener-0.5d0*m%grid(igrid)%uold(ind,idim+1)**2/max(m%grid(igrid)%uold(ind,1),r%smallr)
        end do
  
        ! Update momentum
#ifdef GRAV
        do idim=1,ndim
           m%grid(igrid)%uold(ind,idim+1)=m%grid(igrid)%uold(ind,idim+1)+&
                & max(m%grid(igrid)%uold(ind,1),r%smallr)*m%grid(igrid)%f(ind,idim)*dteff
        end do
#endif
        ! Update total energy
        do idim=1,ndim
           ener=ener+0.5d0*m%grid(igrid)%uold(ind,idim+1)**2/max(m%grid(igrid)%uold(ind,1),r%smallr)
        end do
        m%grid(igrid)%uold(ind,neul)=ener

     end do
     ! End loop over cells
  end do
  ! End loop over grids

#endif

end subroutine synchro_hydro_fine
!################################################################
!################################################################
!################################################################
!################################################################
recursive subroutine r_gravity_hydro_fine(r,g,m,p,mdl,cpu_range,input_size,output_size,ilevel)
  use amr_parameters, only: dp
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
     call mdl_send_request(mdl,MDL_GRAVITY_HYDRO_FINE,next_cpu,next_range,input_size,output_size,ilevel)
     call r_gravity_hydro_fine(r,g,m,p,mdl,next_range,input_size,output_size,ilevel)
  else
     call gravity_hydro_fine(r,g,m,ilevel)
  endif

end subroutine r_gravity_hydro_fine
!################################################################
!################################################################
!################################################################
!################################################################
subroutine gravity_hydro_fine(r,g,m,ilevel)
  use amr_parameters, only: ndim,dp,twotondim
  use amr_commons, only: run_t,global_t,mesh_t
  implicit none
  type(run_t)::r
  type(global_t)::g
  type(mesh_t)::m
  integer::ilevel
  !--------------------------------------------------------------
  ! This routine adds to unew the gravity source terms to unew
  ! with only half a time step. Only the momentum and the
  ! total energy are modified in array unew.
  !--------------------------------------------------------------
  integer::igrid,ind
  real(dp)::d,u,v,w,e_kin,e_prim,d_old,fact

#ifdef HYDRO

  ! Add gravity source term at time t with half time step
  do igrid=m%head(ilevel),m%tail(ilevel)
     do ind=1,twotondim

        d=max(m%grid(igrid)%unew(ind,1),r%smallr)
        u=0.0d0; v=0.0d0; w=0.0d0
        if(ndim>0)u=m%grid(igrid)%unew(ind,2)/d
        if(ndim>1)v=m%grid(igrid)%unew(ind,3)/d
        if(ndim>2)w=m%grid(igrid)%unew(ind,4)/d
        e_kin=0.5d0*d*(u**2+v**2+w**2)
        e_prim=m%grid(igrid)%unew(ind,ndim+2)-e_kin
        d_old=max(m%grid(igrid)%uold(ind,1),r%smallr)
        fact=d_old/d*0.5d0*g%dtnew(ilevel)
#ifdef GRAV
        if(ndim>0)then
           u=u+m%grid(igrid)%f(ind,1)*fact
           m%grid(igrid)%unew(ind,2)=d*u
        endif
        if(ndim>1)then
           v=v+m%grid(igrid)%f(ind,2)*fact
           m%grid(igrid)%unew(ind,3)=d*v
        end if
        if(ndim>2)then
           w=w+m%grid(igrid)%f(ind,3)*fact
           m%grid(igrid)%unew(ind,4)=d*w
        endif
#endif
        e_kin=0.5d0*d*(u**2+v**2+w**2)
        m%grid(igrid)%unew(ind,ndim+2)=e_prim+e_kin
     end do
  end do

#endif

end subroutine gravity_hydro_fine
!###########################################################
!###########################################################
!###########################################################
!###########################################################

