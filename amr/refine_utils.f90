module refine_utils
  type out_refine_fine_t
    integer::make,kill
  end type out_refine_fine_t
contains
!#########################################################################
!#########################################################################
!#########################################################################
!#########################################################################
subroutine m_refine_fine(pst,ilevel)
  use mdl_module
  use ramses_commons, only: pst_t
  use init_refine_basegrid_module, only:r_noct_max,r_noct_min,r_noct_tot,r_noct_used_max
  use load_balance_module, only: m_load_balance, r_balance_part
  implicit none
  type(pst_t)::pst
  integer::ilevel
  !--------------------------------------------------------------------
  ! This routine is the master procedure to refine the AMR grid
  ! from level ilevel to nlevelmax.
  !--------------------------------------------------------------------
  integer::ilev,dummy
  type(out_refine_fine_t)::out_refine_fine
  integer,dimension(1:2)::noct
  integer,dimension(1)::alevel

  associate(s=>pst%s)
  
  if(ilevel==s%r%nlevelmax)return
  if(s%m%noct_tot(ilevel)==0)return

  if(s%r%verbose)write(*,111)ilevel
111 format(' Entering refine_fine for level ',I2)

  ! Create new octs and destroy unecessary octs
  call r_refine_fine(pst,ilevel,1,out_refine_fine,2)

  if(s%r%verbose)write(*,112)out_refine_fine%make
112 format(' ==> Make ',i6,' sub-grids')

  if(s%r%verbose)write(*,113)out_refine_fine%kill
113 format(' ==> Kill ',i6,' sub-grids')

  ! Get total, min and max grid count (only in master)
  do ilev=ilevel+1,s%r%nlevelmax
     call r_noct_tot(pst,ilev,1,s%m%noct_tot(ilev),1)
     call r_noct_min(pst,ilev,1,s%m%noct_min(ilev),1)
     call r_noct_max(pst,ilev,1,s%m%noct_max(ilev),1)
  end do

  ! Get maximum used memory (only in master)
  call r_noct_used_max(pst,ilevel,1,s%m%noct_used_max,1)
  
  ! Load balance all levels across cpus
  call m_load_balance(pst,ilevel)

  ! Get total, min and max grid count (only in master).
  do ilev=ilevel+1,s%r%nlevelmax
     call r_noct_tot(pst,ilev,1,s%m%noct_tot(ilev),1)
     call r_noct_min(pst,ilev,1,s%m%noct_min(ilev),1)
     call r_noct_max(pst,ilev,1,s%m%noct_max(ilev),1)
  end do

  ! Get maximum used memory (only in master)
  call r_noct_used_max(pst,ilevel,1,s%m%noct_used_max,1)

  ! Balance particles across cpus
  if(mdl_threads(s%mdl)>1.AND.ilevel==s%r%levelmin)then
     if(s%r%pic.AND.mod(s%g%nstep_coarse,10)==1)then
        if(s%r%verbose)write(*,*)'Entering balance_part for level',s%r%levelmin
        call r_balance_part(pst,ilevel,1,dummy,0)
     endif
  endif

  end associate
  
end subroutine m_refine_fine
!################################################################
!################################################################
!################################################################
!################################################################
recursive subroutine r_refine_fine(pst,ilevel,input_size,output,output_size)
  use mdl_module
  use ramses_commons, only: pst_t
  use mdl_parameters
  implicit none
  type(pst_t)::pst
  integer,VALUE::input_size
  integer::output_size
  integer::ilevel
  type(out_refine_fine_t)::output

  type(out_refine_fine_t)::next_output
  integer::ncreate,nkill
  integer::rID
  
  if(pst%nLower>0)then
     rID = mdl_send_request(pst%s%mdl,MDL_REFINE_FINE,pst%iUpper+1,input_size,output_size,ilevel)
     call r_refine_fine(pst%pLower,ilevel,input_size,output,output_size)
     call mdl_get_reply(pst%s%mdl,rID,output_size,next_output)
     output%make = output%make + next_output%make
     output%kill = output%kill + next_output%kill
  else
     call refine_fine(pst%s,ilevel,output%make,output%kill)
  endif

end subroutine r_refine_fine
!###############################################################
!###############################################################
!###############################################################
!###############################################################
subroutine refine_fine(s,ilevel,ncreate,nkill)
  use amr_parameters, only: ndim,nhilbert,twotondim
  use amr_commons, only: oct
  use ramses_commons, only: ramses_t
  use marshal, only:pack_fetch_refine, unpack_fetch_refine,&
                    pack_fetch_flag, unpack_fetch_flag
  use cache_commons
  use cache, only:close_cache
  use hash
  use hilbert
  use call_back, only: cache_f
  use nbors_utils_p
  implicit none
  type(ramses_t)::s
  integer::ilevel
  integer::ncreate,nkill
  !---------------------------------------------------------
  ! This routine refines cells at level ilevel if cells
  ! are flagged for refinement and are not already refined.
  ! This routine destroys refinements for cells that are 
  ! not flagged for refinement and are already refined.
  ! For single time-stepping, numerical rules are 
  ! automatically satisfied. For adaptive time-stepping,
  ! numerical rules are checked before refining any cell.
  !---------------------------------------------------------
  integer::igrid,icell,i,j,ibit,ibucket,ilev,ind,inew,ioct
  integer::noct_zero,head_zero,indx_zero
  integer::skip_bit,ikey,true_level
  integer::ind_cell,ind_parent
  integer(kind=8),dimension(0:ndim)::hash_key
  integer(kind=8),dimension(1:nhilbert,1:s%r%nlevelmax)::key_ref
  integer(kind=8),dimension(1:nhilbert)::coarse_key
  integer,dimension(1:s%r%nlevelmax)::n_same,npatch
  integer,dimension(:),allocatable::noct_level,head_level,indx_level
  integer,dimension(:),allocatable::swap_table,swap_tmp
  integer,dimension(0:twotondim-1)::bucket_count,bucket_offset
  logical::ok
  type(oct)::oct_tmp
  type(oct),pointer::gridp
  
  associate(r=>s%r,g=>s%g,m=>s%m)

  !---------------------------------------------------
  ! Step 1: if a cell is flagged for refinement and
  ! if it is not already refined, create its son grid.
  !---------------------------------------------------        
  g%ncreate=0
  m%ifree=m%noct_used+1
  do ilev=ilevel,r%nlevelmax-1

     call open_cache(s,operation_refine,domain_decompos_amr)

     do ioct=m%head(ilev),m%tail(ilev)
        do ind=1,twotondim
           ok   = m%grid(ioct)%flag1(ind)==1 .and. &
                & .not.m%grid(ioct)%refined(ind)
           if(ok)then
              ind_parent=ioct
              ind_cell=ind
              call make_new_oct(s,m%grid(ind_parent),ind_cell,ilev+1)
              g%ncreate=g%ncreate+1
           endif
        end do
     end do

     call close_cache(s,m%grid_dict)

  end do
  ncreate=g%ncreate
  
  !----------------------------------------------------------
  ! Step 2: if the parent cell is not flagged for refinement,
  ! but it is refined, then destroy the child grid.
  !----------------------------------------------------------
  g%nkill=0
  do ilev=ilevel+1,r%nlevelmax

     call open_cache(s,operation_derefine,domain_decompos_amr)

     hash_key(0)=ilev
     do ioct=m%head(ilev),m%tail(ilev)
        hash_key(1:ndim)=m%grid(ioct)%ckey(1:ndim)
        ! Get parent cell using a read-write cache
        call get_parent_cell_p(s,hash_key,m%grid_dict,gridp,icell,flush_cache=.true.,fetch_cache=.true.)
        ok   = gridp%flag1(icell)==0 .and. &
             & gridp%refined(icell)
        if(ok)then
           ! Set grid level to zero
           m%grid(ioct)%lev=0
           ! Set parent cell to "unrefined" status
           gridp%refined(icell)=.false.
           ! Free grid from hash table
           call hash_free(m%grid_dict,hash_key)
           g%nkill=g%nkill+1
        end if
     end do

     call close_cache(s,m%grid_dict)

  end do
  nkill=g%nkill

  !-----------------------------------------------------
  ! Step 3: sort new octs and empty slots according to 
  ! their level (using counting sort algorithm).
  !-----------------------------------------------------
  allocate(noct_level(r%levelmin:r%nlevelmax))
  allocate(head_level(r%levelmin:r%nlevelmax))
  allocate(indx_level(r%levelmin:r%nlevelmax))
  ! Count number of octs per bucket
  noct_level=0
  noct_zero=0
  do ioct=m%tail(ilevel)+1,m%ifree-1
     true_level=m%grid(ioct)%lev
     if(true_level.NE.0)then
        noct_level(true_level)=noct_level(true_level)+1
     else
        noct_zero=noct_zero+1
     end if
  end do
  head_level(ilevel+1)=m%tail(ilevel)+1
  do ilev=ilevel+2,r%nlevelmax
     head_level(ilev)=head_level(ilev-1)+noct_level(ilev-1)
  end do
  head_zero=head_level(r%nlevelmax)+noct_level(r%nlevelmax)

  ! Allocate main swap table
  if(m%ifree.GT.head_level(ilevel+1))then
  allocate(swap_table(head_level(ilevel+1):m%ifree-1))

  ! Build index permutation table
  indx_level=head_level
  indx_zero=head_zero
  do ioct=m%tail(ilevel)+1,m%ifree-1
     true_level=m%grid(ioct)%lev
     if(true_level.NE.0)then
        swap_table(indx_level(true_level))=ioct
        indx_level(true_level)=indx_level(true_level)+1
     else
        swap_table(indx_zero)=ioct
        indx_zero=indx_zero+1
     end if
  end do

  !-----------------------------------------------------
  ! Step 4: sort octs level by level according to their 
  ! Hilbert key using LSD Radix Sort algorithm.
  !-----------------------------------------------------
  ! Loop over levels
  do ilev=ilevel+1,r%nlevelmax
     if(noct_level(ilev)>0)then

        ! Allocate temporary swap table just for the level
        allocate(swap_tmp(head_level(ilev):head_level(ilev)+noct_level(ilev)-1))

        ! Loop over useful bits at that level
        do ibit=ilev,1,-1

           ! Get bit and key to read from
           skip_bit=ndim*(ilev-ibit)
           ikey = skip_bit / bits_per_int(ndim) + 1
           skip_bit = mod(skip_bit, bits_per_int(ndim))

           ! Count octs in buckets
           bucket_count=0
           do inew=head_level(ilev),head_level(ilev)+noct_level(ilev)-1
              ioct=swap_table(inew)
              ibucket=int(ibits(m%grid(ioct)%hkey(ikey),skip_bit,ndim),kind=4)
              bucket_count(ibucket)=bucket_count(ibucket)+1
           end do

           ! Compute offsets
           bucket_offset(0)=head_level(ilev)
           do ibucket=1,twotondim-1
              bucket_offset(ibucket)=bucket_offset(ibucket-1)+bucket_count(ibucket-1)
           end do

           ! Sort according to Hilbert key
           do inew=head_level(ilev),head_level(ilev)+noct_level(ilev)-1
              ioct=swap_table(inew)
              ibucket=int(ibits(m%grid(ioct)%hkey(ikey),skip_bit,ndim),kind=4)
              swap_tmp(bucket_offset(ibucket))=ioct
              bucket_offset(ibucket)=bucket_offset(ibucket)+1
           end do

           ! Store permutations in swap table
           do inew=head_level(ilev),head_level(ilev)+noct_level(ilev)-1
              swap_table(inew)=swap_tmp(inew)
           end do           
        end do

        ! Deallocate tmp swap array
        deallocate(swap_tmp)
     endif
  end do

  !-----------------------------------------------------
  ! Step 5: Apply permutations directly in main memory
  ! Remember: swap_table(inew)=iold means:
  ! New data at position inew COMES FROM
  ! Old data at position iold.
  !-----------------------------------------------------
  ! Perform the swap  
  do j=head_level(ilevel+1),m%ifree-1
     if(j.NE.swap_table(j))then
        hash_key(0)=m%grid(j)%lev
        hash_key(1:ndim)=m%grid(j)%ckey(1:ndim)
        if(m%grid(j)%lev>0)call hash_free(m%grid_dict,hash_key)
        oct_tmp=m%grid(j)
        i=j
        inew=swap_table(j)
        do while(inew.NE.j)
           m%grid(i)=m%grid(inew)
           hash_key(0)=m%grid(inew)%lev
           hash_key(1:ndim)=m%grid(inew)%ckey(1:ndim)
           if(m%grid(inew)%lev>0)then
              call hash_free(m%grid_dict,hash_key)
              call hash_setp(m%grid_dict,hash_key,m%grid(i))
           endif
           swap_table(i)=i
           i=inew
           inew=swap_table(inew)
        end do
        m%grid(i)=oct_tmp
        hash_key(0)=m%grid(i)%lev
        hash_key(1:ndim)=m%grid(i)%ckey(1:ndim)
        if(m%grid(i)%lev>0)then
           call hash_setp(m%grid_dict,hash_key,m%grid(i))
        end if
        swap_table(i)=i
     endif
  end do
  deallocate(swap_table)
  endif

  !-----------------------------------------------------
  ! Step 6: Clean up final AMR structure
  !-----------------------------------------------------
  do ilev=ilevel+1,r%nlevelmax
     m%head(ilev)=head_level(ilev)
     m%tail(ilev)=head_level(ilev)+noct_level(ilev)-1
     m%noct(ilev)=noct_level(ilev)
  end do
  m%noct_used=m%tail(r%nlevelmax)
  deallocate(noct_level,head_level,indx_level)

  !-----------
  ! Super-octs
  !-----------
  do ilev=1,r%nlevelmax
     npatch(ilev)=twotondim**ilev
  end do
  do ilev=ilevel+1,r%nlevelmax
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
  end do

  end associate


end subroutine refine_fine
!###############################################################
!###############################################################
!###############################################################
!###############################################################
subroutine pack_flush_refine(grid,msg_size,msg_array)
  use amr_parameters, only: ndim,twotondim
  use hydro_parameters, only: nvar
  use amr_commons, only: oct
  use cache_commons, only: msg_large_realdp
  type(oct)::grid
  integer::msg_size
  integer,dimension(1:msg_size),optional::msg_array

  integer::ind,ivar,idim
  type(msg_large_realdp)::msg

#ifdef HYDRO
  do ivar=1,nvar
     do ind=1,twotondim
        msg%realdp_hydro(ind,ivar)=grid%uold(ind,ivar)
     end do
  end do
#endif
  
#ifdef GRAV
  do idim=1,ndim
     do ind=1,twotondim
        msg%realdp_poisson(ind,idim)=grid%f(ind,idim)
     end do
  end do
  do ind=1,twotondim
     msg%realdp_poisson(ind,ndim+1)=grid%phi(ind)
     msg%realdp_poisson(ind,ndim+2)=grid%phi_old(ind)
  end do
#endif

  msg_array=transfer(msg,msg_array)

end subroutine pack_flush_refine
!###############################################################
!###############################################################
!###############################################################
!###############################################################
subroutine unpack_flush_refine(grid,msg_size,msg_array,hash_key)
  use amr_parameters, only: ndim,twotondim
  use hydro_parameters, only: nvar
  use amr_commons, only: oct
  use cache_commons, only: msg_large_realdp
  type(oct)::grid
  integer::msg_size
  integer,dimension(1:msg_size),optional::msg_array
  integer(kind=8),dimension(0:ndim)::hash_key

  integer::ind,ivar,idim
  type(msg_large_realdp)::msg

  grid%lev=hash_key(0)
  grid%ckey(1:ndim)=hash_key(1:ndim)
  msg=transfer(msg_array,msg)

  do ind=1,twotondim
     grid%refined(ind)=.false.
  end do
  
#ifdef HYDRO
  do ind=1,twotondim
     do ivar=1,nvar
        grid%uold(ind,ivar)=msg%realdp_hydro(ind,ivar)
     end do
  end do
#endif
  
#ifdef GRAV
  do ind=1,twotondim
     do idim=1,ndim
        grid%f(ind,idim)=msg%realdp_poisson(ind,idim)
     end do
     grid%phi(ind)=msg%realdp_poisson(ind,ndim+1)
     grid%phi_old(ind)=msg%realdp_poisson(ind,ndim+2)
  end do
#endif

end subroutine unpack_flush_refine
!###############################################################
!###############################################################
!###############################################################
!###############################################################
subroutine init_flush_derefine(grid,hash_key)
  use amr_parameters, only: ndim,twotondim
  use amr_commons, only: oct
  type(oct)::grid
  integer(kind=8),dimension(0:ndim)::hash_key

  grid%lev=hash_key(0)
  grid%ckey(1:ndim)=hash_key(1:ndim)
  grid%refined(1:twotondim)=.true.
  
end subroutine init_flush_derefine
!###############################################################
!###############################################################
!###############################################################
!###############################################################
subroutine pack_flush_derefine(grid,msg_size,msg_array)
  use amr_parameters, only: twotondim
  use amr_commons, only: oct
  use cache_commons, only: msg_int4
  type(oct)::grid
  integer::msg_size
  integer,dimension(1:msg_size),optional::msg_array

  integer::ind
  type(msg_int4)::msg

  do ind=1,twotondim     
     if(grid%refined(ind))then
        msg%int4(ind)=1
     else
        msg%int4(ind)=0
     endif
  end do
  msg_array=transfer(msg,msg_array)

end subroutine pack_flush_derefine
!###############################################################
!###############################################################
!###############################################################
!###############################################################
subroutine unpack_flush_derefine(grid,msg_size,msg_array,hash_key)
  use amr_parameters, only: ndim,twotondim
  use amr_commons, only: oct
  use cache_commons, only: msg_int4
  type(oct)::grid
  integer::msg_size
  integer,dimension(1:msg_size),optional::msg_array
  integer(kind=8),dimension(0:ndim)::hash_key

  integer::ind
  type(msg_int4)::msg

  grid%lev=hash_key(0)
  grid%ckey(1:ndim)=hash_key(1:ndim)
  msg=transfer(msg_array,msg)
  do ind=1,twotondim
     if(grid%refined(ind))then
        if(msg%int4(ind)==0)then
           grid%refined(ind)=.false.
        endif
     endif
  end do

end subroutine unpack_flush_derefine
!###############################################################
!###############################################################
!###############################################################
!###############################################################
subroutine make_new_oct(s,parent,icell,ilevel)
  USE, INTRINSIC :: ISO_C_BINDING, ONLY: c_associated
  use mdl_module
  use amr_parameters, only: ndim,nhilbert,twotondim,twondim,nvector
  use amr_commons, only:nbor,oct
  use hydro_parameters, only: nvar
  use ramses_commons, only: ramses_t
  use nbors_utils_p
  use cache_commons
  use hilbert
  use hash
  implicit none
#ifndef WITHOUTMPI
  include 'mpif.h'
#endif
  type(ramses_t)::s
  integer::ilevel
  integer::icell
  type(oct)::parent
  !--------------------------------------------------------------
  ! This routine creates a children oct at level ilevel.
  ! ilevel is thus the level of the new children oct.
  ! The new oct is labelled using the new index ichild.
  ! The parent cell is labeled with the parent oct index iparent
  ! and the cell index icell (from 1 to 8).
  !--------------------------------------------------------------
  integer::idim,ivar,ind,inbor,nstride,grid_cpu
  integer(kind=8),dimension(1:nhilbert)::hk
  integer(kind=8),dimension(1:ndim)::ix
  integer(kind=8),dimension(1:ndim)::cart_key
  integer(kind=8),dimension(0:ndim)::hash_key
  integer,dimension(0:twondim)::igrid_nbor,ind_nbor
  real(dp),dimension(0:twondim,1:nvar)::u1
  real(dp),dimension(1:twotondim,1:nvar)::u2
  type(nbor),dimension(0:twondim)::grid_nbor
  type(oct),pointer::child

  associate(r=>s%r,g=>s%g,m=>s%m,mdl=>s%mdl)

#if !defined(WITHOUTMPI) && !defined(MDL2)
  ! If counter is good, check on incoming messages and perform actions
  if(mdl%mail_counter==32)then
     call check_mail(s,MPI_REQUEST_NULL,m%grid_dict)
     mdl%mail_counter=0
  endif
  mdl%mail_counter=mdl%mail_counter+1
#endif

  !=================================
  ! Create new octs into main memory
  !=================================
  ! Compute Cartesian keys of new octs
  do idim=1,ndim
     nstride=2**(idim-1)
     cart_key(idim)=2*parent%ckey(idim)+MOD((icell-1)/nstride,2)
  end do
  hash_key(0)=ilevel
  hash_key(1:ndim)=cart_key(1:ndim)

  ! Compute Hilbert keys of new octs
  ix(1:ndim)=cart_key(1:ndim)
  hk(1:nhilbert)=hilbert_key(ix,ilevel-1)

  ! Check if grid sits inside processor boundaries
  if (m%domain(ilevel)%in_rank(hk)) then

     ! Set grid index to a virtual grid in local main memory
     child => m%grid(m%ifree)

     ! Go to next main memory free line
     m%ifree=m%ifree+1
     if(m%ifree.GT.r%ngridmax)then
        write(*,*)'No more free memory'
        write(*,*)'while refining...'
        write(*,*)'Increase ngridmax'
        call mdl_abort(mdl)
     end if
     call hash_setp(m%grid_dict,hash_key,child)
  ! Otherwise, determine parent processor and use the cache
  else
     grid_cpu = m%domain(ilevel)%get_rank(hk)
     ! If next cache line is occupied, free it.
     if(m%occupied(m%free_cache))call destage(s,r%ngridmax+m%free_cache,m%grid_dict)
     ! Set grid index to a virtual grid in local cache memory
     child => m%grid(r%ngridmax+m%free_cache)
     m%occupied(m%free_cache)=.true.
     m%parent_cpu(m%free_cache)=grid_cpu
     m%dirty(m%free_cache)=.true.
     ! Go to next free cache line
     m%free_cache=m%free_cache+1
     m%ncache=m%ncache+1
     if(m%free_cache.GT.r%ncachemax)m%free_cache=1
     if(m%ncache.GT.r%ncachemax)m%ncache=r%ncachemax
     ! Insert new grid in hash table
     call hash_setp(m%grid_dict,hash_key,child)
  endif

  child%lev=ilevel
  child%ckey(1:ndim)=int(cart_key(1:ndim),kind=4)
  child%hkey(1:nhilbert)=hk(1:nhilbert)
  child%refined(1:twotondim)=.false.
  child%flag1(1:twotondim)=0
  child%flag2(1:twotondim)=0
  child%superoct=1

  ! Set status of parent cell to "refined"
  parent%refined(icell)=.true.

  !=========================================================
  ! Inject parent hydro variables into new children ones
  !=========================================================     
#ifdef HYDRO

  ! Interpolate hydro variables
  do ivar=1,nvar
     do ind=1,twotondim
        child%uold(ind,ivar)=parent%uold(icell,ivar)
     enddo
  end do
  
  ! In case one wants to interpolate using high-order schemes
  if(r%interpol_type>0)then
     
     ! Get 2ndim neighboring father cells with read-only cache
     call get_twondim_nbor_parent_cell_p(s,hash_key,m%grid_dict,grid_nbor,ind_nbor,flush_cache=.false.,fetch_cache=.true.)
     do inbor=0,twondim
        do ivar=1,nvar
           u1(inbor,ivar)=grid_nbor(inbor)%p%uold(ind_nbor(inbor),ivar)
        end do
     end do
     do inbor=1,twondim
        call unlock_cache_p(s,grid_nbor(inbor)%p)
     end do
     
     ! Interpolate
     call interpol_hydro(u1,u2,r%interpol_var,r%interpol_type,r%smallr)
     
     ! Store hydro variables
     do ivar=1,nvar
        do ind=1,twotondim
           child%uold(ind,ivar)=u2(ind,ivar)
        enddo
     end do
     
  endif
  
#endif
  
#ifdef GRAV
  
  ! Interpolate (straight injection) gravity variables
  do ind=1,twotondim
     child%f(ind,1:ndim)=parent%f(icell,1:ndim)
     child%phi(ind)=parent%phi(icell)
     child%phi_old(ind)=parent%phi_old(icell)
  enddo
  
#endif

  end associate
     
end subroutine make_new_oct
!###############################################################
!###############################################################
!###############################################################
!###############################################################
end module refine_utils
