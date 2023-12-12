#if NDIM==3

recursive subroutine r_compute_clump_properties(s,pst,c,p,ilevel,input_size)
  use mdl_module
  use ramses_commons, only: pst_t,ramses_t
  use mdl_parameters
  use clfind_commons
  implicit none
  type(pst_t)::pst
  type(peak_t)::p
  type(clump_t)::c
  type(ramses_t)::s
  integer,VALUE::input_size
  integer::output_size
  integer::ilevel

  integer::rID

  if(pst%nLower>0)then
     rID = mdl_send_request(pst%s%mdl,MDL_RESET_RHO,pst%iUpper+1,input_size,0,ilevel)!,output_size,ilevel)
     call r_compute_clump_properties(s,pst%pLower,c,p,action,ilevel,input_size)!,ntest,output_size)
     call mdl_get_reply(pst%s%mdl,rID,0)!,output_size,next_ntest)
  else
     call compute_clump_properties(pst%s%r,pst%s%g,c,p,ilevel)!,ntest)
  endif

end subroutine r_compute_clump_properties
!################################################################
!################################################################
!################################################################
!################################################################



subroutine compute_clump_properties(r,g,p,c,ilevel)
    use amr_commons
    use clfind_commons
    !use multigrid_fine_coarse, only:pack_fetch_phi,unpack_fetch_phi
    !use rho_fine_module, only:init_flush_rho,pack_flush_rho,unpack_flush_rho
    !use godunov_fine_module, only: init_flush_godunov,pack_flush_godunov,unpack_flush_godunov
    !use marshal, only: pack_fetch_refine,unpack_fetch_refine
    !use boundaries, only: init_bound_refine
    !use cache_commons
    !use cache
    !use nbors_utils
    implicit none
#ifndef WITHOUTMPI
  integer::info
#endif
    type(run_t)::r
    type(global_t)::g
    type(peak_t)::p
    type(clump_t)::c
    !type(msg_twin_realdp)::dummy_twin_realdp
    !type(msg_large_realdp)::dummy_large_realdp
    !type(oct),pointer::gridp,gridn,gridpm
    !integer,dimension(1:ndim)::ckey,ckey_nbor
    !integer(kind=8),dimension(0:ndim)::hash_cell,hash_nbor
    real(dp)::dx_loc,vol_loc
    integer::ipart,grid,peak_nr,ilevel,global_peak_id,ipeak,plevel
    real(dp)::zero=0
    !variables needed temporarily store cell properties
    real(dp)::d=0, vol=0
    ! variables related to the size of a cell on a given level
    real(dp)::dx,dx_loc,scale,vol_loc
    real(dp),dimension(1:nlevelmax)::volume
    real(dp),dimension(1:3)::skip_loc,xcell
    real(dp),dimension(1:twotondim,1:3)::xc
    integer::nx_loc,ind,ix,iy,iz,idim
    logical,dimension(1:ndim)::period
    logical::periodic

#ifndef WITHOUTMPI
    integer::i
    real(dp)::tot_mass_tot
#endif

    period(1)=(nx==1)
    period(2)=(ny==1)
    period(3)=(nz==1)

    periodic=period(1)
    periodic=periodic.or.period(2)
    periodic=periodic.or.period(3)

    !peak-patch related arrays before sharing information with other cpus

    npeaks = c%npart
    ntest = p%npart

    min_dens=huge(zero)
    n_cells=0; n_cells_halo=0
    halo_mass=0d0; clump_mass=0d0; clump_vol=0d0
    center_of_mass=0d0; clump_velocity=0d0
    peak_pos=0d0

    if(verbose)write(*,*)'Entering compute clump properties'


    !if(r%hydro)then
    !    call open_cache(s,table=m%grid_dict,data_size=storage_size(m%grid(1))/32,&
    !            hilbert=m%domain,pack_size=storage_size(dummy_large_realdp)/32,&
    !            pack=pack_fetch_refine,unpack=unpack_fetch_refine,&
    !            init=init_flush_godunov, flush=pack_flush_godunov,&
    !            combine=unpack_flush_godunov, bound=init_bound_refine)
    !else
    !    call open_cache(s,table=m%grid_dict,data_size=storage_size(m%grid(1))/32,&
    !            hilbert=m%domain,pack_size=storage_size(dummy_twin_realdp)/32,&
    !            pack=pack_fetch_phi,unpack=unpack_fetch_phi,&
    !            init=init_flush_rho, flush=pack_flush_rho, combine=unpack_flush_rho)
    !endif

    !--------------------------------------------------------------------------
    ! loop over all cells above the threshold
    !--------------------------------------------------------------------------
    do ipart=1,ntest
        global_peak_id=p%pid(ipart)
        if (global_peak_id /=0 ) then
            call get_local_peak_id(global_peak_id,peak_nr)
                
            dx_loc=r%boxlen/2**p%levelp(ipart)
            vol_loc=dx_loc**ndim
            !ckey(1:ndim)=int(c%xp(peak_nr,1:ndim)/dx_loc)
            !hash_cell(0)=c%levelp(peak_nr)+1
            !hash_cell(1:ndim)=ckey(1:ndim)
            !call get_parent_cell(s,hash_cell,m%grid_dict,gridp,icellp,flush_cache=.true.,fetch_cache=.true.)
             
            d = p%denp(ipart)!c%denp(peak_nr)
    
            ! Cell volume
            vol=vol_loc
    
            ! Number of leaf cells per clump
            c%n_cells(peak_nr)=c%n_cells(peak_nr)+1
    
            ! Min density
            c%min_dens(peak_nr)=min(d,c%min_dens(peak_nr))
    
            ! Clump mass
            c%clump_mass(peak_nr)=c%clump_mass(peak_nr)+vol*d
    
            ! Clump volume
            c%clump_vol(peak_nr)=c%clump_vol(peak_nr)+vol
    
            ! Clump center of mass location
            c%center_of_mass(peak_nr,1:3)=c%center_of_mass(peak_nr,1:3)+vol*d*c%xp(peak_nr,1:3)
    
            ! Clump center of mass velocity
            if (r%hydro)c%clump_velocity(peak_nr,1:3)=c%clump_velocity(peak_nr,1:3)+vol*p%uold(ipart,2:4)
   
        end if
    end do
    
    !--------------------------------------------------------------------------
    ! Loop over local peaks and identify true peak positions
    !--------------------------------------------------------------------------
    do ipeak=1,npeaks
        ! Peak cell coordinates
        plevel=c%levelp(ipeak)
        dx=0.5D0**plevel
        xcell(1:ndim)=c%xp(ipeak,1:3)
        call true_max(xcell(1),xcell(2),xcell(3),plevel)
        c%peak_pos(ipeak,1:3)=xcell(1:3)
    end do

    ! Compute specific quantities
    do ipeak=1,npeaks
        if (c%relevance(ipeak)>0..and.c%n_cells(ipeak)>0)then
        c%center_of_mass(ipeak,1:3)=c%center_of_mass(ipeak,1:3)/c%clump_mass(ipeak)
        c%clump_velocity(ipeak,1:3)=c%clump_velocity(ipeak,1:3)/c%clump_mass(ipeak)
        end if
        ! Initialize halo mass to clump mass
        c%halo_mass(ipeak)=clump_mass(ipeak)
    end do

    ! Calculate total mass above threshold
    tot_mass=sum(c%clump_mass(1:npeaks))

    ! Compute further properties of the clumps
    do ipeak=1,npeaks
        if (c%relevance(ipeak)>0..and.c%n_cells(ipeak)>0)then
        av_dens(ipeak)=c%clump_mass(ipeak)/c%clump_vol(ipeak)
        end if
    end do

    ! For periodic boxes, recompute center of mass relative to peak position
    if(periodic)then
        center_of_mass=0d0;
        do ipart=1,ntest
            global_peak_id=p%pid(ipart)!flag2(icellp(ipart))
            if (global_peak_id /=0 ) then
                call get_local_peak_id(global_peak_id,peak_nr)

                dx_loc=r%boxlen/2**p%levelp(ipart)
                vol_loc=dx_loc**ndim
                xcell(1:ndim)=p%xp(ipart,1:ndim)

                do idim=1,ndim
                    if (period(idim) .and. (xcell(idim)-c%peak_pos(peak_nr,idim))>boxlen*0.5)xcell(idim)=xcell(idim)-boxlen
                    if (period(idim) .and. (xcell(idim)-c%peak_pos(peak_nr,idim))<boxlen*(-0.5))xcell(idim)=xcell(idim)+boxlen
                end do

                d = p%denp(ipart)

                ! Cell volume
                vol=vol_loc

                ! Clump center of mass location
                c%center_of_mass(peak_nr,1:3)=c%center_of_mass(peak_nr,1:3)+vol*d*xcell(1:3)

            end if
        end do

        ! Compute specific quantity
        do ipeak=1,npeaks
            if (c%relevance(ipeak)>0..and.c%n_cells(ipeak)>0)then
            c%center_of_mass(ipeak,1:3)=c%center_of_mass(ipeak,1:3)/c%clump_mass(ipeak)
            end if
        end do
    endif



end subroutine compute_clump_properties
!################################################################
!################################################################
!################################################################
!################################################################



recursive subroutine r_merge_clumps(s,pst,c,p,action,ilevel,input_size)
  use mdl_module
  use ramses_commons, only: pst_t,ramses_t
  use mdl_parameters
  use clfind_commons
  implicit none
  character(len=9)::action
  type(pst_t)::pst
  type(peak_t)::p
  type(clump_t)::c
  type(ramses_t)::s
  integer,VALUE::input_size
  integer::output_size
  integer::ilevel

  integer::rID

  if(pst%nLower>0)then
     rID = mdl_send_request(pst%s%mdl,MDL_RESET_RHO,pst%iUpper+1,input_size,0,ilevel)!,output_size,ilevel)
     call r_merge_clumps(s,pst%pLower,c,p,action,ilevel,input_size)!,ntest,output_size)
     call mdl_get_reply(pst%s%mdl,rID,0)!,output_size,next_ntest)
  else
     call merge_clumps(pst%s%r,pst%s%g,c,p,action,ilevel)!,ntest)
  endif

end subroutine r_merge_clumps


subroutine merge_clumps(r,g,c,p,action,ilevel)
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
                    call get_local_peak_id(g,merge_to,jpeak)
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
                    call get_local_peak_id(g,merge_to,jpeak)
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
                call get_local_peak_id(g,gridp%flag2(icellp),ipeak)
                merge_to=new_peak(ipeak)
                !call get_local_peak_id(g,merge_to,jpeak)
                gridp%flag2(icellp)=merge_to
                !update test particles peak
                p%pid(ipart) = merge_to
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
                  call get_local_peak_id(g,merge_to,jpeak)
                  c%ind_halo(ipeak)=c%ind_halo(jpeak)
               endif
            end do
        end do
        ! Compute halo masses
        halo_mass=0
        n_cells_halo=0
        do ipeak=1,npeaks
            merge_to=c%ind_halo(ipeak)
            call get_local_peak_id(g,merge_to,jpeak)
            c%halo_mass(jpeak)=c%halo_mass(jpeak)+c%clump_mass(ipeak)
            c%n_cells_halo(jpeak)=c%n_cells_halo(jpeak)+c%n_cells(ipeak)
        end do

        ! Assign back halo mass to peak
        do ipeak=1,npeaks
            merge_to=c%ind_halo(ipeak)
            call get_local_peak_id(g,merge_to,jpeak)
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





!################################################################
!################################################################
!################################################################
!################################################################
subroutine write_clump_properties(r,to_file)
    use amr_commons
    use pm_commons,ONLY:mp
    use hydro_commons,ONLY:mass_sph
    use clfind_commons
    use mpi_mod
    implicit none
  #ifndef WITHOUTMPI
    integer,parameter::tag=1101
    integer::dummy_io,info,info2
  #endif
    logical::to_file
    !---------------------------------------------------------------------------
    ! this routine writes the clump properties to screen and to file
    !---------------------------------------------------------------------------
  
    integer::i,j,jj,ilun,ilun2,n_rel,n_rel_tot,nx_loc
    real(dp)::rel_mass,rel_mass_tot,scale,particle_mass=0
    character(LEN=80)::fileloc,filedir
    character(LEN=5)::nchar,ncharcpu
    real(dp),dimension(1:npeaks)::peakd
    integer,dimension(1:npeaks)::ind_sort
  
  #ifndef WITHOUTMPI
    real(dp)::particle_mass_tot
  #endif
  
    if (.not. to_file)return
  
    if(r%hydro)then
        particle_mass=MINVAL(mp, MASK=(mp > 0))
    else
        particle_mass=mass_sph
    endif
  
    if(to_file)then
       ilun=20
       ilun2=22
    else
       ilun=6
       ilun2=6
    end if
  
    ! print results in descending order to screen/file
    rel_mass=0
    n_rel=0
  
    if (to_file .eqv. .true.) then
       ! first create directories
       call title(ifout,nchar)
       filedir='output_'//TRIM(nchar)
       call create_output_dirs(filedir)
       ! Wait for the token
  #ifndef WITHOUTMPI
       if(IOGROUPSIZE>0) then
          if (mod(g%myid-1,IOGROUPSIZE)/=0) then
             call MPI_RECV(dummy_io,1,MPI_INTEGER,g%myid-1-1,tag,&
                  & MPI_COMM_WORLD,MPI_STATUS_IGNORE,info2)
          end if
       endif
  #endif
  
       if(IOGROUPSIZEREP>0)then
          call title(((g%myid-1)/IOGROUPSIZEREP)+1,ncharcpu)
          fileloc=TRIM(filedir)//'/group_'//TRIM(ncharcpu)//'/clump_'//TRIM(nchar)//'.txt'
       else
          fileloc=TRIM(filedir)//'/clump_'//TRIM(nchar)//'.txt'
       endif
       call title(g%myid,nchar)
       fileloc=TRIM(fileloc)//TRIM(nchar)
       open(unit=ilun,file=fileloc,form='formatted')
  
       if(saddle_threshold>0)then
          call title(ifout,nchar)
          if(IOGROUPSIZEREP>0)then
             call title(((g%myid-1)/IOGROUPSIZEREP)+1,ncharcpu)
             fileloc=TRIM(filedir)//'/group_'//TRIM(ncharcpu)//'/halo_'//TRIM(nchar)//'.txt'
          else
             fileloc=TRIM(filedir)//'/halo_'//TRIM(nchar)//'.txt'
          endif
          call title(g%myid,nchar)
          fileloc=TRIM(fileloc)//TRIM(nchar)
          open(unit=ilun2,file=fileloc,form='formatted')
       endif
    end if
  
    if (to_file .or. myid==1)then
       write(ilun,'(144A)')'   index  halo   lev   parent      ncell    peak_x             peak_y             peak_z     '//&
            '        rho-               rho+               rho_av             mass_cl            relevance   '
       if(saddle_threshold>0)then
          write(ilun2,'(135A)')'     index      ncell    peak_x             peak_y             peak_z     '//&
               '        rho+               mass      '
       endif
    end if
  
    if (particlebased_clump_output) then ! write particle based data
      do j=npeaks,1,-1
        jj=c%sortp(j)
        if (c%relevance(jj) > relevance_threshold .and. clmp_mass_pb(jj) > mass_threshold*particle_mass)then
          write(ilun,'(I8,X,I2,X,I10,X,I10,8(X,1PE18.9E2))')&
                 jj+ipeak_start(myid)&
                 ,lev_peak(jj)&
                 ,new_peak(jj)&
                 ,n_cells(jj)&
  #ifdef UNBINDINGCOM
                 ,clmp_com_pb(jj,1)&
                 ,clmp_com_pb(jj,2)&
                 ,clmp_com_pb(jj,3)&
  #else
                 ,peak_pos(jj,1)&
                 ,peak_pos(jj,2)&
                 ,peak_pos(jj,3)&
  #endif
                 ,min_dens(jj)&
                 ,max_dens(jj)&
                 ,clmp_mass_pb(jj)/clump_vol(jj)&
                 ,clmp_mass_pb(jj)&
                 ,relevance(jj)
           rel_mass=rel_mass+clmp_mass_exclusive(jj)
           n_rel=n_rel+1
        end if
  
        if(saddle_threshold>0)then
          if(ind_halo(jj).EQ.jj+ipeak_start(myid).AND.clmp_mass_pb(jj) > mass_threshold*particle_mass)then
             write(ilun2,'(I10,X,I10,5(X,1PE18.9E2))')&
                    jj+ipeak_start(myid)&
                    ,n_cells_halo(jj)&
  #ifdef UNBINDINGCOM
                    ,clmp_com_pb(jj,1)&
                    ,clmp_com_pb(jj,2)&
                    ,clmp_com_pb(jj,3)&
  #else
                    ,peak_pos(jj,1)&
                    ,peak_pos(jj,2)&
                    ,peak_pos(jj,3)&
  #endif
                    ,max_dens(jj)&
                    ,clmp_mass_pb(jj)
          endif
        endif
      end do
  
    else ! write cell based data
  
      do j=npeaks,1,-1
         jj=ind_sort(j)
         if (relevance(jj) > relevance_threshold .and. halo_mass(jj) > mass_threshold*particle_mass)then
            write(ilun,'(I8,X,I8,1X,I2,X,I10,X,I10,8(X,1PE18.9E2))')&
                 jj+ipeak_start(myid)&
                 ,ind_halo(jj)&
                 ,lev_peak(jj)&
                 ,new_peak(jj)&
                 ,n_cells(jj)&
                 ,peak_pos(jj,1)&
                 ,peak_pos(jj,2)&
                 ,peak_pos(jj,3)&
                 ,min_dens(jj)&
                 ,max_dens(jj)&
                 ,clump_mass(jj)/clump_vol(jj)&
                 ,clump_mass(jj)&
                 ,relevance(jj)
            rel_mass=rel_mass+clump_mass(jj)
            n_rel=n_rel+1
         end if
         if(saddle_threshold>0)then
            if(ind_halo(jj).EQ.jj+ipeak_start(myid).AND.halo_mass(jj) > mass_threshold*particle_mass)then
               write(ilun2,'(I10,X,I10,5(X,1PE18.9E2))')&
                    jj+ipeak_start(myid)&
                    ,n_cells_halo(jj)&
                    ,peak_pos(jj,1)&
                    ,peak_pos(jj,2)&
                    ,peak_pos(jj,3)&
                    ,max_dens(jj)&
                    ,halo_mass(jj)
            endif
         endif
      end do
    end if
  
    if (to_file)then
       close(ilun)
       if(saddle_threshold>0)then
          close(ilun2)
       endif
    end if
  
       ! Send the token
  #ifndef WITHOUTMPI
       if(IOGROUPSIZE>0) then
          if(mod(myid,IOGROUPSIZE)/=0 .and.(myid.lt.ncpu))then
             dummy_io=1
             call MPI_SEND(dummy_io,1,MPI_INTEGER,myid-1+1,tag, &
                  & MPI_COMM_WORLD,info2)
          end if
       endif
  #endif
  
  #ifndef WITHOUTMPI
    call MPI_ALLREDUCE(n_rel,n_rel_tot,1,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,info)
    n_rel=n_rel_tot
    call MPI_ALLREDUCE(rel_mass,rel_mass_tot,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
    rel_mass=rel_mass_tot
  #else
    n_rel_tot = n_rel
    rel_mass_tot = rel_mass
  #endif
    if(myid==1)then
       if(clinfo)write(*,'(A,1PE12.5)')' Total mass [code units] above threshold =',tot_mass
       if(clinfo)write(*,'(A,I10,A,1PE12.5)')' Total mass [code units] in',n_rel_tot,' listed clumps =',rel_mass_tot
    endif
  
  end subroutine write_clump_properties

















#endif