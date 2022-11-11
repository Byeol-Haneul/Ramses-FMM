module init_refine_basegrid_module
contains
!#########################################################################
!#########################################################################
!#########################################################################
!#########################################################################
subroutine m_init_refine_basegrid(pst)
  use ramses_commons, only: pst_t
#ifdef GRAV
  use rho_fine_module, only: m_rho_fine
#endif
  implicit none
  type(pst_t)::pst
  !--------------------------------------------------------------------
  ! This routine is the master procedure to set the base grid
  ! and initialize all cell-based variables within it.
  !--------------------------------------------------------------------
  integer::dummy

  associate(r=>pst%s%r,g=>pst%s%g,m=>pst%s%m,p=>pst%s%p,mdl=>pst%s%mdl)
  
  if(r%verbose)write(*,*)'Entering init_refine_basegrid'

  ! Call recursive slave routine
  call r_init_refine_basegrid(pst,r%levelmin,1)

  ! Get total, min and max grid count (only in master).
  call r_noct_tot(pst,r%levelmin,1,m%noct_tot(r%levelmin),1)
  call r_noct_min(pst,r%levelmin,1,m%noct_min(r%levelmin),1)
  call r_noct_max(pst,r%levelmin,1,m%noct_max(r%levelmin),1)
  call r_noct_used_max(pst,r%levelmin,1,m%noct_used_max,1)

  ! Initialize hydro variables on the base grid
  if(r%hydro)call m_init_flow_fine(pst,r%levelmin)

#ifdef GRAV
  ! Compute total mass density from gas and particles on the base grid
  call m_rho_fine(pst,r%levelmin)
#endif
  
  end associate

end subroutine m_init_refine_basegrid
!###############################################
!###############################################
!###############################################
!###############################################
recursive subroutine r_collect_noct(pst,ilevel,input_size,noct,output_size)
  use mdl_module
  use ramses_commons, only: pst_t
  use mdl_parameters
  implicit none
  type(pst_t)::pst
!  integer,VALUE::input_size
  integer,VALUE::input_size
  integer::output_size
  integer,dimension(1:input_size)::ilevel
  integer,dimension(1:output_size)::noct

  integer,dimension(1:output_size)::next_noct
  integer::rID

  if(pst%nLower>0)then
     rID = mdl_send_request(pst%s%mdl,MDL_COLLECT_NOCT,pst%iUpper+1,input_size,output_size,ilevel)
     call r_collect_noct(pst%pLower,ilevel,input_size,noct,output_size)
     call mdl_get_reply(pst%s%mdl,rID,output_size,next_noct)
     noct=noct+next_noct
  else
     noct=0
     noct(pst%s%g%myid)=pst%s%m%noct(ilevel(1))
  endif

end subroutine r_collect_noct
!###############################################
!###############################################
!###############################################
!###############################################
recursive subroutine r_noct_tot(pst,ilevel,input_size,noct_tot,output_size)
  use mdl_module
  use ramses_commons, only: pst_t
  use mdl_parameters
  implicit none
  type(pst_t)::pst
  integer,VALUE::input_size
  integer::output_size
  integer::ilevel
  integer::noct_tot

  integer::next_noct_tot
  integer::rID

  if(pst%nLower>0)then
     rID = mdl_send_request(pst%s%mdl,MDL_NOCT_TOT,pst%iUpper+1,input_size,output_size,ilevel)
     call r_noct_tot(pst%pLower,ilevel,input_size,noct_tot,output_size)
     call mdl_get_reply(pst%s%mdl,rID,output_size,next_noct_tot)
     noct_tot=noct_tot+next_noct_tot
  else
     noct_tot=pst%s%m%noct(ilevel)
  endif

end subroutine r_noct_tot
!###############################################
!###############################################
!###############################################
!###############################################
recursive subroutine r_noct_max(pst,ilevel,input_size,noct_max,output_size)
  use mdl_module
  use ramses_commons, only: pst_t
  use mdl_parameters
  implicit none
  type(pst_t)::pst
  integer,VALUE::input_size
  integer::output_size
  integer::ilevel
  integer::noct_max

  integer::next_noct_max
  integer::rID

  if(pst%nLower>0)then
     rID = mdl_send_request(pst%s%mdl,MDL_NOCT_MAX,pst%iUpper+1,input_size,output_size,ilevel)
     call r_noct_max(pst%pLower,ilevel,input_size,noct_max,output_size)
     call mdl_get_reply(pst%s%mdl,rID,output_size,next_noct_max)
     noct_max=MAX(noct_max,next_noct_max)
  else
     noct_max=pst%s%m%noct(ilevel)
  endif

end subroutine r_noct_max
!###############################################
!###############################################
!###############################################
!###############################################
recursive subroutine r_noct_used_max(pst,ilevel,input_size,noct_used_max,output_size)
  use mdl_module
  use ramses_commons, only: pst_t
  use mdl_parameters
  implicit none
  type(pst_t)::pst
  integer,VALUE::input_size
  integer::output_size
  integer::ilevel
  integer::noct_used_max

  integer::next_noct_used_max
  integer::rID

  if(pst%nLower>0)then
     rID = mdl_send_request(pst%s%mdl,MDL_NOCT_USED_MAX,pst%iUpper+1,input_size,output_size,ilevel)
     call r_noct_used_max(pst%pLower,ilevel,input_size,noct_used_max,output_size)
     call mdl_get_reply(pst%s%mdl,rID,output_size,next_noct_used_max)
     noct_used_max=MAX(noct_used_max,next_noct_used_max)
  else
     noct_used_max=pst%s%m%noct_used
  endif

end subroutine r_noct_used_max
!###############################################
!###############################################
!###############################################
!###############################################
recursive subroutine r_noct_min(pst,ilevel,input_size,noct_min,output_size)
  use mdl_module
  use ramses_commons, only: pst_t
  use mdl_parameters
  implicit none
  type(pst_t)::pst
  integer,VALUE::input_size
  integer::output_size
  integer::ilevel
  integer::noct_min

  integer::next_noct_min
  integer::rID

  if(pst%nLower>0)then
     rID = mdl_send_request(pst%s%mdl,MDL_NOCT_MIN,pst%iUpper+1,input_size,output_size,ilevel)
     call r_noct_min(pst%pLower,ilevel,input_size,noct_min,output_size)
     call mdl_get_reply(pst%s%mdl,rID,output_size,next_noct_min)
     noct_min=MIN(noct_min,next_noct_min)
  else
     noct_min=pst%s%m%noct(ilevel)
  endif

end subroutine r_noct_min
!###############################################
!###############################################
!###############################################
!###############################################
recursive subroutine r_gather_noct_max(pst,ilevel,input_size,noct_max,output_size)
  use mdl_module
  use ramses_commons, only: pst_t
  use mdl_parameters
  implicit none
  type(pst_t)::pst
  integer,VALUE::input_size
  integer::output_size
  integer,dimension(1:input_size)::ilevel
  integer,dimension(1:output_size)::noct_max

  integer,dimension(1:output_size)::next_noct_max
  integer::rID

  if(pst%nLower>0)then
     rID = mdl_send_request(pst%s%mdl,MDL_GATHER_NOCT_MAX,pst%iUpper+1,input_size,output_size,ilevel)
     call r_gather_noct_max(pst%pLower,ilevel,input_size,noct_max,output_size)
     call mdl_get_reply(pst%s%mdl,rID,output_size,next_noct_max)
     noct_max(1)=MAX(noct_max(1),next_noct_max(1))
  else
     noct_max(1)=pst%s%m%noct(ilevel(1))
  endif

end subroutine r_gather_noct_max
!###############################################
!###############################################
!###############################################
!###############################################
recursive subroutine r_init_refine_basegrid(pst,ilevel,input_size)
  use mdl_module
  use ramses_commons, only: pst_t
  use mdl_parameters
  implicit none
  type(pst_t)::pst
  integer,VALUE::input_size
  integer::ilevel

  integer::rID

  if(pst%nLower>0)then
     rID = mdl_send_request(pst%s%mdl,MDL_INIT_REFINE_BASEGRID,pst%iUpper+1,input_size,0,ilevel)
     call r_init_refine_basegrid(pst%pLower,ilevel,input_size)
     call mdl_get_reply(pst%s%mdl,rID,0)
  else
     call init_refine_basegrid(pst%s%r,pst%s%g,pst%s%m,ilevel)
  endif

end subroutine r_init_refine_basegrid
!################################################################
!################################################################
!################################################################
!################################################################
subroutine init_refine_basegrid(r,g,m,ilevel)
  use amr_parameters, only: nhilbert,ndim,twotondim
  use amr_commons, only: run_t,global_t,mesh_t
  use hilbert
  use hash, only: hash_setp
  implicit none
  type(run_t)::r
  type(global_t)::g
  type(mesh_t)::m
  integer::ilevel
  !-------------------------------------------------------
  ! This routine builds a fully refined Cartesian grid
  ! at level ilevel. Always starts at levelmin.
  !-------------------------------------------------------
  integer::i,igrid,ioct,ilev,istart
  integer(kind=8)::ikey
  integer(kind=8),dimension(1:nhilbert)::hk
  integer(kind=8),dimension(1:ndim)::ix
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
  hk=0
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
     call hash_setp(m%grid_dict,hash_key,m%grid(igrid))
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
end module init_refine_basegrid_module
