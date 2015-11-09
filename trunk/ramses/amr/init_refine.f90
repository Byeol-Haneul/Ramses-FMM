!################################################################
!################################################################
!################################################################
!################################################################
subroutine init_refine
  use amr_commons
  use pm_commons
  implicit none
  !-------------------------------------------
  ! This routine builds the initial AMR grid
  !-------------------------------------------
  integer::ilevel

  if(myid==1)write(*,*)'Building initial AMR grid'
  init=.true.

  ! Base refinement
  do ilevel=1,levelmin
     call flag
     call refine
  end do

  ! Further refinements if necessary
  do ilevel=levelmin+1,nlevelmax
     if(initfile(levelmin).ne.' '.and.initfile(ilevel).eq.' ')exit
     if(hydro)call init_flow
     if(ivar_refine==0)call init_refmap
     call flag
     call refine
     if(nremap>0)call load_balance
     if(numbtot(1,ilevel)==0)exit
  end do 

  ! Final pass to initialize the flow
  init=.false.
  if(hydro)call init_flow

end subroutine init_refine
!################################################################
!################################################################
!################################################################
!################################################################
subroutine init_refine_2
  !--------------------------------------------------------------
  ! This routine builds additional refinements to the
  ! the initial AMR grid for filetype ne 'grafic'
  !--------------------------------------------------------------
  use amr_commons
  use hydro_commons
  use pm_commons
  use poisson_commons
  use sort, only: lsd_radix_sort_particles, apply_particle_permutation
  use hilbert, only: hilbert_for_particle
  implicit none
  integer::ilevel,i,ivar, np, offset

  if(filetype.eq.'grafic')return

  offset = part_level_offset(ilevel)
  np = npart - part_level_offset(ilevel)

  ! Compute hilbert keys (probably move outside of this routine)
  call hilbert_for_particle(offset, npart - offset, 0, ilevel)
  
  ! Compute a permutation that sorts ALL particles starting from offset
  call lsd_radix_sort_particles(offset, npart - offset, ilevel, ilevel, .true.)

  call apply_particle_permutation(offset, npart - offset, ilevel) 



  do i=levelmin,nlevelmax+1

     call refine_coarse
     do ilevel=1,nlevelmax
        call build_comm(ilevel)
        call make_virtual_fine_int(cpu_map(1),ilevel)
        call refine_fine(ilevel)
        if(hydro)call init_flow_fine(ilevel)
     end do

     if(nremap>0)call load_balance

     do ilevel=levelmin,nlevelmax
        if(poisson)call rho_fine(ilevel,2,.false.)
     end do

     do ilevel=nlevelmax,levelmin,-1
        if(hydro)then
           call upload_fine(ilevel)
           do ivar=1,nvar
              call make_virtual_fine_dp(uold(1,ivar),ilevel)
           end do
           if(simple_boundary)call make_boundary_hydro(ilevel)
        endif
     end do

     do ilevel=nlevelmax,1,-1
        call flag_fine(ilevel,2)
     end do
     call flag_coarse

  end do

end subroutine init_refine_2
!################################################################
!################################################################
!################################################################
!################################################################
