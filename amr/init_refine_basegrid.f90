!#########################################################################
!#########################################################################
!#########################################################################
!#########################################################################
subroutine m_init_refine_basegrid(r,g,m,p,mdl)
  use amr_commons, only: run_t,global_t,mesh_t
  use pm_commons, only: part_t
  use mdl_commons, only: mdl_t
  implicit none
  type(run_t)::r
  type(global_t)::g
  type(mesh_t)::m
  type(part_t)::p
  type(mdl_t)::mdl
  !--------------------------------------------------------------------
  ! This routine is the master procedure to set the base grid
  ! and initialize all cell-based variables within it.
  !--------------------------------------------------------------------

  if(r%verbose)write(*,*)'Entering init_refine_basegrid'

  ! Call recursive slave routine
  call r_init_refine_basegrid(r,g,m,p,mdl,mdl%ncpu,1,0,r%levelmin)

  ! Get total, min and max grid count (only in master).
  call r_noct_tot(r,g,m,p,mdl,mdl%ncpu,1,1,r%levelmin,m%noct_tot(r%levelmin))
  call r_noct_min(r,g,m,p,mdl,mdl%ncpu,1,1,r%levelmin,m%noct_min(r%levelmin))
  call r_noct_max(r,g,m,p,mdl,mdl%ncpu,1,1,r%levelmin,m%noct_max(r%levelmin))
  call r_noct_used_max(r,g,m,p,mdl,mdl%ncpu,1,1,r%levelmin,m%noct_used_max)

  ! Initialize hydro variables on the base grid
  if(r%hydro)call m_init_flow_fine(r,g,m,p,mdl,r%levelmin)

  ! Compute total mass density from gas and particles on the base grid
  call m_rho_fine(r,g,m,p,mdl,r%levelmin)

end subroutine m_init_refine_basegrid
!###############################################
!###############################################
!###############################################
!###############################################
recursive subroutine r_collect_noct(r,g,m,p,mdl,cpu_range,input_size,output_size,ilevel,noct)
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
  integer,dimension(1:output_size)::noct

  integer::next_range,next_cpu
  integer,dimension(1:output_size)::next_noct

  next_range=cpu_range/2
  next_cpu=g%myid+next_range

  if(next_range>0)then
     call mdl_send_request(mdl,MDL_COLLECT_NOCT,next_cpu,next_range,input_size,output_size,ilevel)
     call r_collect_noct(r,g,m,p,mdl,next_range,input_size,output_size,ilevel,noct)
     call mdl_get_reply(mdl,next_cpu,output_size,next_noct)
     noct=noct+next_noct
  else
     noct=0
     noct(g%myid)=m%noct(ilevel)
  endif

end subroutine r_collect_noct
!###############################################
!###############################################
!###############################################
!###############################################
recursive subroutine r_noct_tot(r,g,m,p,mdl,cpu_range,input_size,output_size,ilevel,noct_tot)
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
  integer::ilevel,noct_tot

  integer::next_range,next_cpu
  integer::next_noct_tot

  next_range=cpu_range/2
  next_cpu=g%myid+next_range

  if(next_range>0)then
     call mdl_send_request(mdl,MDL_NOCT_TOT,next_cpu,next_range,input_size,output_size,ilevel)
     call r_noct_tot(r,g,m,p,mdl,next_range,input_size,output_size,ilevel,noct_tot)
     call mdl_get_reply(mdl,next_cpu,output_size,next_noct_tot)
     noct_tot=noct_tot+next_noct_tot
  else
     noct_tot=m%noct(ilevel)
  endif

end subroutine r_noct_tot
!###############################################
!###############################################
!###############################################
!###############################################
recursive subroutine r_noct_max(r,g,m,p,mdl,cpu_range,input_size,output_size,ilevel,noct_max)
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
  integer::ilevel,noct_max

  integer::next_range,next_cpu
  integer::next_noct_max

  next_range=cpu_range/2
  next_cpu=g%myid+next_range

  if(next_range>0)then
     call mdl_send_request(mdl,MDL_NOCT_MAX,next_cpu,next_range,input_size,output_size,ilevel)
     call r_noct_max(r,g,m,p,mdl,next_range,input_size,output_size,ilevel,noct_max)
     call mdl_get_reply(mdl,next_cpu,output_size,next_noct_max)
     noct_max=MAX(noct_max,next_noct_max)
  else
     noct_max=m%noct(ilevel)
  endif

end subroutine r_noct_max
!###############################################
!###############################################
!###############################################
!###############################################
recursive subroutine r_noct_used_max(r,g,m,p,mdl,cpu_range,input_size,output_size,ilevel,noct_used_max)
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
  integer::ilevel,noct_used_max

  integer::next_range,next_cpu
  integer::next_noct_used_max

  next_range=cpu_range/2
  next_cpu=g%myid+next_range

  if(next_range>0)then
     call mdl_send_request(mdl,MDL_NOCT_USED_MAX,next_cpu,next_range,input_size,output_size,ilevel)
     call r_noct_used_max(r,g,m,p,mdl,next_range,input_size,output_size,ilevel,noct_used_max)
     call mdl_get_reply(mdl,next_cpu,output_size,next_noct_used_max)
     noct_used_max=MAX(noct_used_max,next_noct_used_max)
  else
     noct_used_max=m%noct_used
  endif

end subroutine r_noct_used_max
!###############################################
!###############################################
!###############################################
!###############################################
recursive subroutine r_noct_min(r,g,m,p,mdl,cpu_range,input_size,output_size,ilevel,noct_min)
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
  integer::ilevel,noct_min

  integer::next_range,next_cpu
  integer::next_noct_min

  next_range=cpu_range/2
  next_cpu=g%myid+next_range

  if(next_range>0)then
     call mdl_send_request(mdl,MDL_NOCT_MIN,next_cpu,next_range,input_size,output_size,ilevel)
     call r_noct_min(r,g,m,p,mdl,next_range,input_size,output_size,ilevel,noct_min)
     call mdl_get_reply(mdl,next_cpu,output_size,next_noct_min)
     noct_min=MIN(noct_min,next_noct_min)
  else
     noct_min=m%noct(ilevel)
  endif

end subroutine r_noct_min
!###############################################
!###############################################
!###############################################
!###############################################
recursive subroutine r_gather_noct_max(r,g,m,p,mdl,cpu_range,input_size,output_size,ilevel,noct_max)
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
  integer::ilevel,noct_max

  integer::next_range,next_cpu
  integer::next_noct_max

  next_range=cpu_range/2
  next_cpu=g%myid+next_range

  if(next_range>0)then
     call mdl_send_request(mdl,MDL_GATHER_NOCT_MAX,next_cpu,next_range,input_size,output_size,ilevel)
     call r_gather_noct_max(r,g,m,p,mdl,next_range,input_size,output_size,ilevel,noct_max)
     call mdl_get_reply(mdl,next_cpu,output_size,next_noct_max)
     noct_max=MAX(noct_max,next_noct_max)
  else
     noct_max=m%noct(ilevel)
  endif

end subroutine r_gather_noct_max
!###############################################
!###############################################
!###############################################
!###############################################
recursive subroutine r_init_refine_basegrid(r,g,m,p,mdl,cpu_range,input_size,output_size,ilevel)
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
     call mdl_send_request(mdl,MDL_INIT_REFINE_BASEGRID,next_cpu,next_range,input_size,output_size,ilevel)
     call r_init_refine_basegrid(r,g,m,p,mdl,next_range,input_size,output_size,ilevel)
  else
     call init_refine_basegrid(r,g,m,p,ilevel)
  endif

end subroutine r_init_refine_basegrid
!################################################################
!################################################################
!################################################################
!################################################################
subroutine init_refine_basegrid(r,g,m,p,ilevel)
  use amr_parameters, only: nhilbert,ndim,twotondim
  use amr_commons, only: run_t,global_t,mesh_t
  use pm_commons, only: part_t
  use hilbert
  use hash, only: hash_set
  implicit none
  type(run_t)::r
  type(global_t)::g
  type(mesh_t)::m
  type(part_t)::p
  integer::ilevel
  !-------------------------------------------------------
  ! This routine builds a fully refined Cartesian grid
  ! at level ilevel. Always starts at levelmin.
  !-------------------------------------------------------
  integer::i,igrid,ioct,ilev,istart
  integer(kind=8)::ikey
  integer(kind=8),dimension(1:nhilbert)::hk=0
  integer(kind=8),dimension(1:ndim)::ix=0
  integer(kind=8),dimension(0:ndim)::hash_key
  integer(kind=8),dimension(1:nhilbert,1:r%nlevelmax)::key_ref
  integer(kind=8),dimension(1:nhilbert)::coarse_key
  integer,dimension(1:r%nlevelmax)::n_same,npatch

  ! Compute starting grid index at that level
  if(ilevel.EQ.r%levelmin)then
     istart=1
  else
     istart=m%tail(ilevel-1)+1
  endif

  ! New grid in current level
  igrid=istart-1

  ! Loop over the Cartesian grid in Hilbert order
  do ikey=m%domain(ilevel)%b(1,g%myid-1), m%domain(ilevel)%b(1,g%myid)-1
     ! Compute Cartesian index from Hilbert index
     hk(1)=ikey
     ix=hilbert_reverse(hk,ilevel-1)
     ! Insert new grid in main array
     igrid=igrid+1
     if(igrid==istart)m%head(ilevel)=istart
     m%tail(ilevel)=igrid
     m%noct(ilevel)=m%noct(ilevel)+1
     m%noct_used=m%noct_used+1
     m%grid(igrid)%lev=ilevel
     m%grid(igrid)%ckey(1:ndim)=int(ix(1:ndim),kind=4)
     m%grid(igrid)%hkey(1:nhilbert)=hk(1:nhilbert)
     m%grid(igrid)%refined(1:twotondim)=.false.
     ! Insert new grid in hash table
     hash_key(0)=ilevel
     hash_key(1:ndim)=ix(1:ndim)
     call hash_set(m%grid_dict,hash_key,igrid)
  end do

  !-----------
  ! Super-octs
  !-----------
  do i=1,ilevel
     npatch(i)=twotondim**i
  end do
  ilev=ilevel
  n_same=0
  key_ref=0
  key_ref(1,1:r%nlevelmax)=-1
  do ioct=m%head(ilev),m%tail(ilev)
     m%grid(ioct)%superoct=1
     coarse_key(1:nhilbert)=m%grid(ioct)%hkey(1:nhilbert)
     do i=1,MIN(ilev-1,r%nsuperoct)
        coarse_key(1:nhilbert)=coarsen_key(coarse_key(1:nhilbert),ilev-1) ! ilev-1 used to speed up only
        if(eq_keys(coarse_key(1:nhilbert),key_ref(1:nhilbert,i)))then
           n_same(i)=n_same(i)+1
        else
           n_same(i)=1
           key_ref(1:nhilbert,i)=coarse_key(1:nhilbert)
        endif
        if(n_same(i).EQ.npatch(i))then
           m%grid(ioct-npatch(i)+1:ioct)%superoct=npatch(i)
        endif
     end do
  end do
  
end subroutine init_refine_basegrid
!################################################################
!################################################################
!################################################################
!################################################################

