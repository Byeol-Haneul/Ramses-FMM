!#########################################################################
!#########################################################################
!#########################################################################
!#########################################################################
subroutine m_refine_fine(pst,ilevel)
  use ramses_commons, only: pst_t
  implicit none
  type(pst_t)::pst
  integer::ilevel
  !--------------------------------------------------------------------
  ! This routine is the master procedure to refine the AMR grid
  ! from level ilevel to nlevelmax.
  !--------------------------------------------------------------------
  integer::ilev
  integer,dimension(1:2)::noct

  associate(s=>pst%s)
  
  if(ilevel==s%r%nlevelmax)return
  if(s%m%noct_tot(ilevel)==0)return

  if(s%r%verbose)write(*,111)ilevel
111 format(' Entering refine_fine for level ',I2)

  ! Create new octs and destroy unecessary octs
  call r_refine_fine(pst,1,2,ilevel,noct)

  if(s%r%verbose)write(*,112)noct(1)
112 format(' ==> Make ',i6,' sub-grids')

  if(s%r%verbose)write(*,113)noct(2)
113 format(' ==> Kill ',i6,' sub-grids')

  ! Get total, min and max grid count (only in master)
  do ilev=ilevel+1,s%r%nlevelmax
     call r_noct_tot(pst,1,1,ilev,s%m%noct_tot(ilev))
     call r_noct_min(pst,1,1,ilev,s%m%noct_min(ilev))
     call r_noct_max(pst,1,1,ilev,s%m%noct_max(ilev))
  end do

  ! Get maximum used memory (only in master)
  call r_noct_used_max(pst,1,1,ilevel,s%m%noct_used_max)

  ! Load balance all levels across cpus
  call m_load_balance(pst,ilevel)

  ! Get total, min and max grid count (only in master).
  do ilev=ilevel+1,s%r%nlevelmax
     call r_noct_tot(pst,1,1,ilev,s%m%noct_tot(ilev))
     call r_noct_min(pst,1,1,ilev,s%m%noct_min(ilev))
     call r_noct_max(pst,1,1,ilev,s%m%noct_max(ilev))
  end do

  ! Get maximum used memory (only in master)
  call r_noct_used_max(pst,1,1,ilevel,s%m%noct_used_max)

  ! Balance particles across cpus
  if(s%mdl%ncpu>1.AND.ilevel==s%r%levelmin)then
     if(s%r%pic.AND.mod(s%g%nstep_coarse,10)==1)then
        if(s%r%verbose)write(*,*)'Entering balance_part for level',s%r%levelmin
        call r_balance_part(pst,1,0,ilevel)
     endif
  endif

  end associate
  
end subroutine m_refine_fine
!################################################################
!################################################################
!################################################################
!################################################################
recursive subroutine r_refine_fine(pst,input_size,output_size,ilevel,noct)
  use ramses_commons, only: pst_t
  use mdl_parameters
  implicit none
  type(pst_t)::pst
  integer::input_size,output_size
  integer::ilevel
  integer,dimension(1:2)::noct

  integer,dimension(1:2)::next_noct
  integer::ncreate,nkill
  
  if(pst%nLower>0)then
     call mdl_send_request(pst%s%mdl,MDL_REFINE_FINE,pst%iUpper+1,input_size,output_size,ilevel)
     call r_refine_fine(pst%pLower,input_size,output_size,ilevel,noct)
     call mdl_get_reply(pst%s%mdl,pst%iUpper+1,output_size,next_noct)
     noct=noct+next_noct
  else
     call refine_fine(pst%s,ilevel,ncreate,nkill)
     noct(1)=ncreate
     noct(2)=nkill
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
  use cache_commons
  use hash
  use hilbert
  use call_back, only: cache_f
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
  integer::parent_cell
  integer,external::get_parent_cell
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
              call make_new_oct(s,ind_parent,ind_cell,ilev+1)
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
        parent_cell=get_parent_cell(s,hash_key,m%grid_dict,.true.,.true.)
        igrid=(parent_cell-1)/twotondim+1
        icell=parent_cell-(igrid-1)*twotondim
        ok   = m%grid(igrid)%flag1(icell)==0 .and. &
             & m%grid(igrid)%refined(icell)
        if(ok)then
           ! Set grid level to zero
           m%grid(ioct)%lev=0
           ! Set parent cell to "unrefined" status
           m%grid(igrid)%refined(icell)=.false.
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
              call hash_set(m%grid_dict,hash_key,i)
           endif
           swap_table(i)=i
           i=inew
           inew=swap_table(inew)
        end do
        m%grid(i)=oct_tmp
        hash_key(0)=m%grid(i)%lev
        hash_key(1:ndim)=m%grid(i)%ckey(1:ndim)
        if(m%grid(i)%lev>0)then
           call hash_set(m%grid_dict,hash_key,i)
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
subroutine pack_fetch_refine(grid,msg_size,msg_array)
  use amr_parameters, only: ndim,twotondim
  use hydro_parameters, only: nvar
  use amr_commons, only: oct
  use cache_commons, only: msg_large_realdp
  type(oct)::grid
  integer::msg_size
  integer,dimension(1:msg_size)::msg_array

  integer::ind,ivar
  type(msg_large_realdp)::msg

  do ind=1,twotondim
     if(grid%refined(ind))then
        msg%int4(ind)=1
     else
        msg%int4(ind)=0
     endif
  end do
  
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
#endif
  end do

  msg_array=transfer(msg,msg_array)

end subroutine pack_fetch_refine
!###############################################################
!###############################################################
!###############################################################
!###############################################################
subroutine unpack_fetch_refine(grid,msg_size,msg_array)
  use amr_parameters, only: ndim,twotondim
  use hydro_parameters, only: nvar
  use amr_commons, only: oct
  use cache_commons, only: msg_large_realdp
  type(oct)::grid
  integer::msg_size
  integer,dimension(1:msg_size)::msg_array

  integer::ind,ivar
  type(msg_large_realdp)::msg

  msg=transfer(msg_array,msg)

  do ind=1,twotondim
     if(msg%int4(ind)==1)then
        grid%refined(ind)=.true.
     else
        grid%refined(ind)=.false.
     endif
  end do
  
#ifdef HYDRO
  do ivar=1,nvar
     do ind=1,twotondim
        grid%uold(ind,ivar)=msg%realdp_hydro(ind,ivar)
     end do
  end do
#endif

#ifdef GRAV
  do idim=1,ndim
     do ind=1,twotondim
        grid%f(ind,idim)=msg%realdp_poisson(ind,idim)
     end do
  end do
  do ind=1,twotondim
     grid%phi(ind)=msg%realdp_poisson(ind,ndim+1)
     grid%phi_old(ind)=msg%realdp_poisson(ind,ndim+2)
  end do
#endif

end subroutine unpack_fetch_refine
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
  integer,dimension(1:msg_size)::msg_array

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
subroutine unpack_flush_refine(grid,msg_size,msg_array)
  use amr_parameters, only: ndim,twotondim
  use hydro_parameters, only: nvar
  use amr_commons, only: oct
  use cache_commons, only: msg_large_realdp
  type(oct)::grid
  integer::msg_size
  integer,dimension(1:msg_size)::msg_array

  integer::ind,ivar,idim
  type(msg_large_realdp)::msg

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
subroutine init_flush_derefine(grid,msg_size)
  use amr_parameters, only: twotondim
  use amr_commons, only: oct
  type(oct)::grid
  integer::msg_size

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
  integer,dimension(1:msg_size)::msg_array

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
subroutine unpack_flush_derefine(grid,msg_size,msg_array)
  use amr_parameters, only: twotondim
  use amr_commons, only: oct
  use cache_commons, only: msg_int4
  type(oct)::grid
  integer::msg_size
  integer,dimension(1:msg_size)::msg_array

  integer::ind
  type(msg_int4)::msg

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
subroutine make_new_oct(s,iparent,icell,ilevel)
  use amr_parameters, only: ndim,nhilbert,twotondim,twondim,nvector
  use hydro_parameters, only: nvar
  use ramses_commons, only: ramses_t
  use cache_commons
  use hilbert
  use hash
  implicit none
#ifndef WITHOUTMPI
  include 'mpif.h'
#endif
  type(ramses_t)::s
  integer::ilevel
  integer::iparent,icell
  !--------------------------------------------------------------
  ! This routine creates a children oct at level ilevel.
  ! ilevel is thus the level of the new children oct.
  ! The new oct is labelled using the new index ichild.
  ! The parent cell is labeled with the parent oct index iparent
  ! and the cell index icell (from 1 to 8).
  !--------------------------------------------------------------
  integer::idim,ivar,ichild,ind,inbor,nstride,grid_cpu
  integer(kind=8),dimension(1:nhilbert)::hk
  integer(kind=8),dimension(1:ndim)::ix
  integer(kind=8),dimension(1:ndim)::cart_key
  integer(kind=8),dimension(0:ndim)::hash_key
  integer,dimension(0:twondim)::igrid_nbor,ind_nbor
  real(dp),dimension(0:twondim,1:nvar)::u1
  real(dp),dimension(1:twotondim,1:nvar)::u2

  associate(r=>s%r,g=>s%g,m=>s%m,mdl=>s%mdl)

#ifndef WITHOUTMPI
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
     cart_key(idim)=2*m%grid(iparent)%ckey(idim)+MOD((icell-1)/nstride,2)
  end do

  ! Compute Hilbert keys of new octs
  ix(1:ndim)=cart_key(1:ndim)
  hk(1:nhilbert)=hilbert_key(ix,ilevel-1)

  ! Check if grid sits inside processor boundaries
  if (m%domain(ilevel)%in_rank(hk)) then

     ! Set grid index to a virtual grid in local main memory
     ichild=m%ifree

     ! Go to next main memory free line
     m%ifree=m%ifree+1
     if(m%ifree.GT.r%ngridmax)then
        write(*,*)'No more free memory'
        write(*,*)'while refining...'
        write(*,*)'Increase ngridmax'
        call mdl_abort
     end if

  ! Otherwise, determine parent processor and use the cache
  else
     grid_cpu = m%domain(ilevel)%get_rank(hk)
     ! If next cache line is occupied, free it.
     if(m%occupied(m%free_cache))call destage(s,r%ngridmax+m%free_cache,m%grid_dict)
     ! Set grid index to a virtual grid in local cache memory
     ichild=r%ngridmax+m%free_cache
     m%occupied(m%free_cache)=.true.
     m%parent_cpu(m%free_cache)=grid_cpu
     m%dirty(m%free_cache)=.true.
     ! Go to next free cache line
     m%free_cache=m%free_cache+1
     m%ncache=m%ncache+1
     if(m%free_cache.GT.r%ncachemax)m%free_cache=1
     if(m%ncache.GT.r%ncachemax)m%ncache=r%ncachemax
  endif

  m%grid(ichild)%lev=ilevel
  m%grid(ichild)%ckey(1:ndim)=int(cart_key(1:ndim),kind=4)
  m%grid(ichild)%hkey(1:nhilbert)=hk(1:nhilbert)
  m%grid(ichild)%refined(1:twotondim)=.false.
  m%grid(ichild)%flag1(1:twotondim)=0
  m%grid(ichild)%flag2(1:twotondim)=0
  m%grid(ichild)%superoct=1

  ! Insert new grid in hash table
  hash_key(0)=ilevel
  hash_key(1:ndim)=cart_key(1:ndim)
  call hash_set(m%grid_dict,hash_key,ichild)

  ! Set status of parent cell to "refined"
  m%grid(iparent)%refined(icell)=.true.

  !=========================================================
  ! Inject parent hydro variables into new children ones
  !=========================================================     
#ifdef HYDRO

  ! Interpolate hydro variables
  do ivar=1,nvar
     do ind=1,twotondim
        m%grid(ichild)%uold(ind,ivar)=m%grid(iparent)%uold(icell,ivar)
     enddo
  end do
  
  ! In case one wants to interpolate using high-order schemes
  if(r%interpol_type>0)then
     
     ! Get 2ndim neighboring father cells with read-only cache
     call get_twondim_nbor_parent_cell(s,hash_key,m%grid_dict,igrid_nbor,ind_nbor,.false.,.true.)
     do inbor=0,twondim
        do ivar=1,nvar
           u1(inbor,ivar)=m%grid(igrid_nbor(inbor))%uold(ind_nbor(inbor),ivar)
        end do
     end do
     do inbor=1,twondim
        call unlock_cache(s,igrid_nbor(inbor))
     end do
     
     ! Interpolate
     call interpol_hydro(u1,u2,r%interpol_var,r%interpol_type,r%smallr)
     
     ! Store hydro variables
     do ivar=1,nvar
        do ind=1,twotondim
           m%grid(ichild)%uold(ind,ivar)=u2(ind,ivar)
        enddo
     end do
     
  endif
  
#endif
  
#ifdef GRAV
  
  ! Interpolate (straight injection) gravity variables
  do ind=1,twotondim
     m%grid(ichild)%f(ind,1:ndim)=m%grid(iparent)%f(icell,1:ndim)
     m%grid(ichild)%phi(ind)=m%grid(iparent)%phi(icell)
     m%grid(ichild)%phi_old(ind)=m%grid(iparent)%phi_old(icell)
  enddo
  
#endif

  end associate
     
end subroutine make_new_oct
!###############################################################
!###############################################################
!###############################################################
!###############################################################
  
