module init_xion_module
contains
!#########################################################################
!#########################################################################
!#########################################################################
!#########################################################################
subroutine m_init_xion(pst)

! Initialize hydrogen ionization state in all cells at all levels from
! density and temperature in the cells, assuming chemical equilibrium.
!-------------------------------------------------------------------------
  use ramses_commons, only: pst_t
  use upload_module, only: m_upload_fine
  implicit none
  type(pst_t)::pst
  integer::ilevel

  if(pst%s%g%myid==1) &
            write(*,*) 'Initialising to CIE ionisation fractions'
  do ilevel=pst%s%r%nlevelmax,pst%s%r%levelmin,-1
      call r_init_xion(pst,ilevel,1)
      call m_upload_fine(pst,ilevel)
  end do

end subroutine m_init_xion

!#########################################################################
!#########################################################################
!#########################################################################
!#########################################################################
recursive subroutine r_init_xion(pst,ilevel,input_size)

! Initialize hydrogen ionization state in all cells at given level from
! density and temperature in the cells, assuming chemical equilibrium.
!-------------------------------------------------------------------------
  use mdl_module
  use ramses_commons, only: pst_t
  use mdl_parameters
  implicit none
  type(pst_t)::pst
  integer,VALUE::input_size
  integer::ilevel

  integer::rID

  if(pst%nLower>0)then
     rID = mdl_send_request( &
             pst%s%mdl,MDL_INIT_XION,pst%iUpper+1,input_size,0,ilevel)
     call r_init_xion(pst%pLower,ilevel,input_size)
     call mdl_get_reply(pst%s%mdl,rID,0)
  else
     call init_xion(pst%s%r,pst%s%g,pst%s%m,pst%s%tables,ilevel)
  endif

end subroutine r_init_xion
!#########################################################################
!#########################################################################
!#########################################################################
!#########################################################################
subroutine init_xion(r,g,m,tables,ilevel)
  use constants
  use amr_parameters, only:dp,ndim,nvector,twotondim
  use hydro_parameters, only: nener, nion
  use amr_commons, only:run_t,global_t,mesh_t
  use cooling_module, only:cooling_t
  use neq_cooling_module, only: cmp_equilibrium_abundances
  use coolrates_module, only: neq_cooling_t
  implicit none
  type(run_t)::r
  type(global_t)::g
  type(mesh_t)::m
  type(neq_cooling_t)::tables
  integer::ilevel
  integer::i,ind,igrid,idim,ngrid,nleaf
  real(kind=8)::scale_nH,scale_T2,scale_l,scale_d,scale_t,scale_v
  real(kind=8)::mu,x
  integer,dimension(1:nvector)::ind_leaf
  real(kind=8),dimension(1:nvector)::nH,T2,ekk,err,emag,Zsolar
  real(dp),dimension(nion)::phI_rates       ! Photoionization rates [s-1]
  real(dp),dimension(1:nvector, 7)::nSpec    !          Species abundances
#if NENER>0
  integer::irad
#endif

  if(r%verbose.and.g%myid==1) &
                  write(*,'("   Entering init_xion for level ",I2)')ilevel

  associate(ixHI=>r%ixHi, ixHII=>r%ixHII, ixHeII=>r%ixHeII, ixHeIII=>r%ixHeIII)

  ! Conversion factor from user units to cgs units
  call units(r,g,scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)
  pHI_rates(:)=0.0                   ! No UV background for the time being

  ! Loop over cells
  do ind=1,twotondim
     ! Loop over octs with vector sweeps
     do igrid=m%head(ilevel),m%tail(ilevel),nvector

        ! Collect a vector of leaf cells
        ngrid=MIN(nvector,m%tail(ilevel)-igrid+1)
        nleaf=0
        do i=1,ngrid
           if(.NOT. m%grid(igrid+i-1)%refined(ind))then
              nleaf=nleaf+1
              ind_leaf(nleaf)=igrid+i-1
           endif
        end do
        if(nleaf.eq.0)cycle

        ! Compute rho
        do i=1,nleaf
           nH(i)=MAX(m%grid(ind_leaf(i))%uold(ind,1),r%smallr)
        end do

        ! Compute metallicity in solar units
        if(r%metal)then
           do i=1,nleaf
              Zsolar(i)=m%grid(ind_leaf(i))%uold(ind,r%imetal)/nH(i)/0.02d0
           end do
        else
           do i=1,nleaf
              Zsolar(i)=r%z_ave
           end do
        endif

        do i=1,nleaf ! Total energy
           T2(i)=m%grid(ind_leaf(i))%uold(ind,5)
        end do

        do i=1,nleaf ! Kinetic energy
           ekk(i)=0.0d0
        end do
        do idim=1,3
           do i=1,nleaf
              ekk(i)=ekk(i)+0.5d0*m%grid(ind_leaf(i))%uold(ind,idim+1)**2/nH(i)
           end do
        end do

        do i=1,nleaf ! Non-thermal energies
           err(i)=0.0d0
        end do
#if NENER>0
        do irad=1,nener
           do i=1,nleaf
              err(i)=err(i)+m%grid(ind_leaf(i))%uold(ind,5+irad)
           end do
        end do
#endif
        do i=1,nleaf ! Magnetic energy
           emag(i)=0.0d0
        end do
#ifdef MHD
        do idim=1,3
           do i=1,nleaf
              emag(i)=emag(i)+0.125d0*(m%grid(ind_leaf(i))%bold(ind,idim) &
                   &              +m%grid(ind_leaf(i))%bold(ind,idim+3))**2
           end do
        end do
#endif
        do i=1,nleaf ! Gas thermal pressure
           T2(i)=(r%gamma-1.0d0)*(T2(i)-ekk(i)-err(i)-emag(i))
        end do

        do i=1,nleaf ! Compute T2=T/mu in Kelvin
           T2(i)=T2(i)/nH(i)*scale_T2
        end do

        do i=1,nleaf ! Compute nH in H/cc
           nH(i)=nH(i)*scale_nH
        end do

        ! Do the main computation of equilibrium abundances
        do i=1,nleaf
          call cmp_equilibrium_abundances(r, tables, &
                         T2(i),nH(i),pHI_rates(:),mu,nSpec(i,:),Zsolar(i))
        end do

        ! nspec:    1=ne, 2=nH2, 3=nHI, 4=nHII, 5=nHeI, 6=nHeII, 7=n_HEIII
        do i=1,nleaf ! Update ionization states in uold

          if(r%isH2) then
              if(r%nrestart==0 .and. r%cosmo) then      ! No primordial H2
                  x = (2.*nSpec(i,2)+nSpec(i,3)) &
                    / (2.*nSpec(i,2)+nSpec(i,3)+nspec(i,4))  ! HI fraction
              else
                  x = nSpec(i,3) &
                    /(2.*nSpec(i,2)+nSpec(i,3)+nspec(i,4))   ! HI fraction
              endif
              m%grid(ind_leaf(i))%uold(ind,r%iIons-1+ixHI) = &
                                         x*m%grid(ind_leaf(i))%uold(ind,1)
          endif

          x = nSpec(i,4) &
              /(2.*nSpec(i,2)+nSpec(i,3)+nspec(i,4))        ! HII fraction
          m%grid(ind_leaf(i))%uold(ind,r%iIons-1+ixHII) = &
                                         x*m%grid(ind_leaf(i))%uold(ind,1)

          if(r%Y_He .gt. 0d0 .and. r%isHe) then
              x = nSpec(i,6) &
                      /(nSpec(i,5)+nSpec(i,6)+nSpec(i,7)) !  HeII fraction
              m%grid(ind_leaf(i))%uold(ind,r%iIons-1+ixHeII) = &
                                         x*m%grid(ind_leaf(i))%uold(ind,1)
              x = nSpec(i,7) &
                      /(nSpec(i,5)+nSpec(i,6)+nSpec(i,7)) ! HeIII fraction
              m%grid(ind_leaf(i))%uold(ind,r%iIons-1+ixHeIII) = &
                                         x*m%grid(ind_leaf(i))%uold(ind,1)
          endif
        end do

     end do
     ! End loop over grid
  end do
  ! End loop over cells
  end associate

end subroutine init_xion
end module init_xion_module