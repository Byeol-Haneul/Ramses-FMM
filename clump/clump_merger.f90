#if NDIM==3
subroutine merge_clumps(r,g,c,p,action)
    use amr_commons
    use clfind_commons
    use multigrid_fine_coarse, only:pack_fetch_phi,unpack_fetch_phi
    use rho_fine_module, only:init_flush_rho,pack_flush_rho,unpack_flush_rho
    use godunov_fine_module, only: init_flush_godunov,pack_flush_godunov,unpack_flush_godunov
    use marshal, only: pack_fetch_refine,unpack_fetch_refine
    use boundaries, only: init_bound_refine
    use cache_commons
    use cache
    use nbors_utils
    implicit none
    character(len=9)::action
    type(run_t)::r
    type(global_t)::g
    type(peak_t)::p
    type(clump_t)::c
    type(msg_twin_realdp)::dummy_twin_realdp
    type(msg_large_realdp)::dummy_large_realdp
    type(oct),pointer::gridp,gridn,gridpm
    integer,dimension(1:ndim)::ckey,ckey_nbor
    integer(kind=8),dimension(0:ndim)::hash_cell,hash_nbor
    real(dp)::dx_loc
#ifndef WITHOUTMPI
    integer::info
#endif
    !---------------------------------------------------------------------------
    ! This routine merges the irrelevant clumps
    ! -clumps are sorted by ascending max density
    ! -irrelevent clumps are merged to most relevant neighbor
    !---------------------------------------------------------------------------
    integer::j,i,merge_to,ipart
    integer::current,nmove,ipeak,jpeak,iter
    integer::nsurvive,nzero,idepth
    integer::ilev,global_peak_id
    real(dp)::value_iij,zero=0,relevance_peak
    integer,dimension(1:npeaks_max)::alive,ind_sort
    real(dp),dimension(1:npeaks_max)::peakd
    logical::do_merge=.false.

#ifndef WITHOUTMPI
    integer::mergelevel_max_global
    integer::nmove_all,nsurvive_all,nzero_all
#endif
    if(r%verbose)then
        if(action.EQ.'relevance')then
        write(*,*)'Now merging irrelevant clumps'
        endif
        if(action.EQ.'saddleden')then
        write(*,*)'Now merging clumps into halos'
        endif
    endif
    !npeaks = npeak_cpu(g%myid)-npeak_cpu(g%myid-1)
`   ! Initialize new_peak array to global peak id
    ! All peaks are alive at the start
    npeaks = c%npart 
    do i=1,npeaks
        if(action.EQ.'relevance')then
            c%alive(i) = 1
        endif
        if(action.EQ.'saddleden')then
            if(c%relevance(i)>relevance_threshold)then
                c%alive(i) = 1
            else
                c%alive(i) = 0
            endif
        endif
    end do
    ! Sort peaks by maximum peak density in ascending order
    
    if (npeaks>0) then
        allocate(testp_sort(npeaks))
        do i=1,npeaks
            denp(i)=c%denp(i)
            testp_sort(i)=i
        end do
        call quick_sort_dp(denp,testp_sort,npeaks)
        deallocate(denp)
    endif

    do i=1,npeaks
        c%sortp(i) = testp_sort(i)
    end do
    deallocate(testp_sort)

    ! Loop over peak levels
    nzero=npeaks_tot
    idepth=0
    do while(nzero>0)

        ! Compute maximum saddle density for each clump
        do i=1,hfree-1
            call get_max(g,i,sparse_saddle_dens)
        end do
        
        ! Merge peaks
        nmove=npeak_cum(g%ncpu)
        iter=0
        do while(nmove>0)
            nmove=0
            do i=npeaks,1,-1
                ipeak=c%sortp(i)
                merge_to=c%idc(ipeak)
                if(p%alive(ipeak)>0)then
                    if(action.EQ.'relevance')then
                        if(sparse_saddle_dens%maxval(ipeak)>0)then
                            relevance_peak=c%denp(ipeak)/sparse_saddle_dens%maxval(ipeak)
                        else
                            relevance_peak=c%denp(ipeak)/density_threshold
                         end if
                        do_merge=(relevance_peak<relevance_threshold)
                    endif
                    if(action.EQ.'saddleden')then
                        do_merge=(sparse_saddle_dens%maxval(ipeak)>saddle_threshold)
                    endif
                    if(do_merge)then
                        if(sparse_saddle_dens%maxloc(ipeak)>0)then
                           call get_local_peak_id(g,sparse_saddle_dens%maxloc(ipeak),jpeak)
                           if(c%denp(jpeak)>c%denp(ipeak))then
                                merge_to=c%idc(jpeak)
                           else if(c%denp(jpeak)==c%denp(ipeak))then
                                merge_to=MIN(c%idc(ipeak),c%idc(jpeak))
                           endif
                        endif
                    endif
                endif
                if(c%idc(ipeak).NE.merge_to)then
                    nmove=nmove+1
                    c%idc(ipeak)=merge_to
                endif
            enddo
            iter=iter+1
        enddo
        ! Transfer matrix elements of merged peaks to surviving peaks
        ! Create new local duplicated peaks and update communicator
        do ipeak=1,hfree-1
            if(c%alive(ipeak)>0)then
                merge_to=c%idc(ipeak)
                if(ipeak.LE.npeaks)then
                    global_peak_id=npeak_cum(g%myid-1)+ipeak
                else
                    global_peak_id=gkey(ipeak)
                endif
                if(merge_to.NE.global_peak_id)then
                    call get_local_peak_id(merge_to,jpeak)
                    current=sparse_saddle_dens%first(ipeak) ! first element of line ipeak
                    do while(current>0) ! walk the line
                        j=sparse_saddle_dens%col(current)
                        value_iij=sparse_saddle_dens%val(current) ! value of the matrix
                        ! Copy the value of density only if larger
                        if(value_iij>get_value(jpeak,j,sparse_saddle_dens))then
                            call set_value(jpeak,j,value_iij,sparse_saddle_dens)
                            call set_value(j,jpeak,value_iij,sparse_saddle_dens)
                        end if
                        current=sparse_saddle_dens%next(current)
                    end do
                    call set_value(jpeak,jpeak,zero,sparse_saddle_dens)
                endif
            endif
        enddo
        ! Set alive to zero for newly merged peaks
        nzero=0
        nsurvive=0
        do ipeak=1,npeaks
            if(c%alive(ipeak)>0)then
                merge_to=c%idc(ipeak)
                if(merge_to.NE.(npeak_cum(g%myid-1)+ipeak))then
                  c%alive(ipeak)=0
                  c%lev_peak(ipeak)=idepth
                  nzero=nzero+1
                else
                    nsurvive=nsurvive+1
                end if
            endif
        end do

        ! Remove all matrix elements corresponding to merged peaks
        do ipeak=1,hfree-1
            current=sparse_saddle_dens%first(ipeak) ! first element of line ipeak
            do while(current>0) ! walk the line
            j=sparse_saddle_dens%col(current)
            if(c%alive(ipeak)==0 .or. c%alive(j)==0)then
                call set_value(ipeak,j,zero,sparse_saddle_dens)
            endif
            current=sparse_saddle_dens%next(current)
            end do
        end do
        if(r%verbose)write(*,*)'level=',idepth,'nmove=',nzero,'survived=',nsurvive
        idepth=idepth+1
    end do
    ! End loop over peak levels
    mergelevel_max=idepth-2 ! last level has no more clumps, also idepth=idepth+1 still happens on last level.
    
    ! Compute maximum saddle density for each surviving clump
    ! Create new local duplicated peaks and update communicator
    do i=1,hfree-1
        call get_max(g,i,sparse_saddle_dens)
    end do

    ! Count surviving peaks
    nsurvive=0
    do ipeak=1,npeaks
        if(c%alive(ipeak)>0)then
            nsurvive=nsurvive+1
        endif
    end do

    if(g%myid==1)then
        if(action.EQ.'relevance')then
           write(*,*)'Found',nsurvive,' relevant peaks'
        endif
        if(action.EQ.'saddleden')then
           write(*,*)'Found',nsurvive,' halos'
        endif
     endif
    if(action.EQ.'relevance')then
        ! Compute relevance
        do ipeak=1,npeaks
            if(c%alive(ipeak)>0)then
            if (sparse_saddle_dens%maxval(ipeak)>0)then
                relevance_peak=c%denp(ipeak)/sparse_saddle_dens%maxval(ipeak)
            else
                relevance_peak=c%denp(ipeak)/density_threshold
            end if
            c%relevance(ipeak)=relevance_peak
            else
            c%relevance(ipeak)=0
            endif
        end do
        ! Merge all peaks to deepest level
        do ilev=idepth-2,0,-1
            do ipeak=1,npeaks
                if(c%lev_peak(ipeak)==ilev)then
                    merge_to=c%idc(ipeak)
                    call get_local_peak_id(merge_to,jpeak)
                    c%idc(ipeak)=c%idc(jpeak)
                endif
            end do
        end do

        if(r%hydro)then
            call open_cache(s,table=m%grid_dict,data_size=storage_size(m%grid(1))/32,&
                    hilbert=m%domain,pack_size=storage_size(dummy_large_realdp)/32,&
                    pack=pack_fetch_refine,unpack=unpack_fetch_refine,&
                    init=init_flush_godunov, flush=pack_flush_godunov,&
                    combine=unpack_flush_godunov, bound=init_bound_refine)
        else
            call open_cache(s,table=m%grid_dict,data_size=storage_size(m%grid(1))/32,&
                    hilbert=m%domain,pack_size=storage_size(dummy_twin_realdp)/32,&
                    pack=pack_fetch_phi,unpack=unpack_fetch_phi,&
                    init=init_flush_rho, flush=pack_flush_rho, combine=unpack_flush_rho)
        endif

        ! Update flag2 field
        do ipart=1,p%naprt
            dx_loc=r%boxlen/2**p%levelp(ipart)
            ckey(1:ndim)=int(p%xp(ipart,1:ndim)/dx_loc)
            hash_cell(0)=p%levelp(ipart)+1
            hash_cell(1:ndim)=ckey(1:ndim)
            call get_parent_cell(s,hash_cell,m%grid_dict,gridp,icellp,flush_cache=.true.,fetch_cache=.true.)
            if (gridp%flag2(icellp)>0)then
                call get_local_peak_id(gridp%flag2(icellp),ipeak)
                merge_to=new_peak(ipeak)
                call get_local_peak_id(merge_to,jpeak)
                gridp%flag2(icellp)=merge_to
            end if
        end do
        call close_cache(s,m%grid_dict)
        
    endif

    if(action.EQ.'saddleden')then
        ! Compute peak index for the halo
        do ipeak=1,npeaks
            c%ind_halo(ipeak)=c%idc(ipeak)
        end do
        do ilev=idepth-2,0,-1
            do ipeak=1,npeaks
               if(c%lev_peak(ipeak)==ilev)then
                  merge_to=c%ind_halo(ipeak)
                  call get_local_peak_id(merge_to,jpeak)
                  c%ind_halo(ipeak)=c%ind_halo(jpeak)
               endif
            end do
        end do
        ! Compute halo masses
        halo_mass=0
        n_cells_halo=0
        do ipeak=1,npeaks
            merge_to=c%ind_halo(ipeak)
            call get_local_peak_id(merge_to,jpeak)
            c%halo_mass(jpeak)=c%halo_mass(jpeak)+c%clump_mass(ipeak)
            c%n_cells_halo(jpeak)=c%n_cells_halo(jpeak)+c%n_cells(ipeak)
        end do

        ! Assign back halo mass to peak
        do ipeak=1,npeaks
            merge_to=c%ind_halo(ipeak)
            call get_local_peak_id(merge_to,jpeak)
            c%halo_mass(ipeak)=c%halo_mass(jpeak)
        end do
    endif

    

end subroutine merge_clumps

!################################################################
!################################################################
!################################################################
!################################################################
subroutine get_max(g,i,mat)
    use amr_commons, only: global_t
    use sparse_matrix
    use clfind_commons,ONLY: gkey,npeak_cum
    type(sparse_mat)::mat
    type(global_t)::g
    integer::i,npeaks
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! get maximum in i-th line by walking the linked list
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    integer::current,icol

    npeaks = npeak_cpu(g%myid)-npeak_cpu(g%myid-1)
  
    mat%maxval(i)=0
    mat%maxloc(i)=0
  
    ! walk the line...
    current=mat%first(i)
    do  while( current /= 0 )
       if(mat%maxval(i)<mat%val(current))then
          mat%maxval(i)=mat%val(current)
          icol=mat%col(current)
          if(icol<=npeaks)then
             mat%maxloc(i)=npeak_cum(g%myid-1)+icol
          else
             mat%maxloc(i)=gkey(icol)
          endif
       end if
       current=mat%next(current)
    end do
  end subroutine get_max

  !################################################################
!################################################################
!################################################################
!################################################################
subroutine get_local_peak_id(g,global_peak_id,local_peak_id)
    use amr_commons
    use clfind_commons
    implicit none
    type(global_t)::g
    integer::local_peak_id
    integer::ihash,ikey,jkey
    integer::global_peak_id

    if(    global_peak_id> npeak_cum(g%myid-1) .and. &
         & global_peak_id<=npeak_cum(g%myid))then
       local_peak_id=global_peak_id-npeak_cum(g%myid-1)
    else
       ihash=MOD(global_peak_id,nhash)+1 ! compute the simple prime hash key
       if(hkey(ihash)==0)then ! hash table is empty
          hkey(ihash)=hfree
          local_peak_id=hfree
          gkey(hfree)=global_peak_id
          hfree=hfree+1
          if(hfree.eq.npeaks_max)then
             write(*,*)'Too many peaks'
             write(*,*)'Increase npeaks_max'
             stop
          endif
       else
          ikey=hkey(ihash) ! collision in the hash table
          do while(ikey>0)
             jkey=ikey
             if(gkey(ikey)==global_peak_id)exit
             ikey=nkey(ikey)
          end do
          if(ikey==0)then ! peak doesn't already exist
             nkey(jkey)=hfree
             local_peak_id=hfree
             gkey(hfree)=global_peak_id
             hfree=hfree+1
             hcollision=hcollision+1
             if(hfree.eq.npeaks_max)then
                write(*,*)'Too many peaks'
                write(*,*)'Increase npeaks_max'
                stop
             endif
          else            ! peak already exists
             local_peak_id=ikey
          end if
       end if
    end if
  
end subroutine get_local_peak_id

#endif