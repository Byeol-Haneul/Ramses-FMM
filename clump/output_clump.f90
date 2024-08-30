module output_clump_module
  use clump_merger_module
contains
!###################################################
!###################################################
!###################################################
!###################################################
recursive subroutine r_output_clump(pst,input_array,input_size,output_array,output_size)
  use mdl_module
  use amr_parameters, only: flen
  use ramses_commons, only: pst_t
  use mdl_parameters
  implicit none
  type(pst_t)::pst
  integer,VALUE::input_size
  integer::output_size
  integer,dimension(1:input_size)::input_array
  integer,dimension(1:output_size)::output_array
  
  character(LEN=flen)::filename,fileloc
  integer::rID,no_halo

  if(pst%nLower>0)then
     rID = mdl_send_request(pst%s%mdl,MDL_OUTPUT_CLUMP,pst%iUpper+1,input_size,output_size,input_array)
     call r_output_clump(pst%pLower,input_array,input_size,output_array,output_size)
     call mdl_get_reply(pst%s%mdl,rID,output_size)
  else
     filename=transfer(input_array,filename)
     ! Write clump and halo properties
     if(pst%s%r%output_clump)then
        call output_clump_properties(pst%s,filename)
     endif
     ! Write peak id of AMR cells
     if(pst%s%r%output_peak)then
        fileloc=TRIM(filename)//'peak.'
        call output_clump_field(pst%s,fileloc)
     endif
     ! Computing and write peak id for dark matter particles'
     if(pst%s%r%output_peak_part.and.pst%s%r%pic)then
        fileloc=TRIM(filename)//'peak_part.'
        call particle_peak_id(pst%s,pst%s%p,no_halo)
        call output_peak_part(pst%s,pst%s%p,fileloc)
     endif
     ! Computing and write peak id for star particles'
     if(pst%s%r%output_peak_star.and.pst%s%r%star)then
        fileloc=TRIM(filename)//'peak_star.'
        call particle_peak_id(pst%s,pst%s%star,no_halo)
        call output_peak_part(pst%s,pst%s%star,fileloc)
     endif
  endif

end subroutine r_output_clump
!###################################################
!###################################################
!###################################################
!###################################################
subroutine output_clump_properties(s,filename)
  use amr_parameters, only: flen, nbin, ndim, dp
  use ramses_commons, only: ramses_t,open_file,close_file
  use clfind_commons
  implicit none
  type(ramses_t)::s
  character(LEN=flen)::filename
  !----------------------------------------------------------------
  ! This routine output the clump properties for each processor
  !----------------------------------------------------------------
  integer::ilun,j
  character(LEN=flen)::fileloc
  integer(kind=8),dimension(s%r%levelmin:s%r%nlevelmax)::nskip
  real(dp)::rad,r200b,rvir,concentration
  real(dp),dimension(1:nbin)::mbin

  associate(r=>s%r,g=>s%g,c=>s%c)

  ! Write clump file
  fileloc=TRIM(filename)//'clump.'

  call open_file(s,fileloc,nskip,ilun)

  do j=1,c%npeak
     if (c%relevance(j) > c%relevance_threshold.AND. &
          & c%clump_mass(j) > c%mass_threshold.AND. &
          & c%halo_mass(j) > c%mass_threshold)then
        write(ilun,'(I10,X,I10,1X,I2,X,I10,X,I10,11(X,1PE18.9E2))')&
             j+c%npeak_cum(g%myid-1)&
             ,c%ind_halo(j)&
             ,c%lev_peak(j)&
             ,c%new_peak(j)&
             ,c%n_cells(j)&
             ,c%peak_pos(j,1),c%peak_pos(j,2),c%peak_pos(j,3)&
             ,c%peak_vel(j,1),c%peak_vel(j,2),c%peak_vel(j,3)&
             ,c%min_dens(j)&
             ,c%max_dens(j)&
             ,c%clump_mass(j)/c%clump_vol(j)&
             ,c%clump_mass(j)&
             ,c%relevance(j)
     end if
  end do

  call close_file(s,fileloc,nskip,ilun)

  ! Write halo file
  if(c%saddle_threshold>0)then

     fileloc=TRIM(filename)//'halo.'

     call open_file(s,fileloc,nskip,ilun)

     do j=1,c%npeak
        if(c%ind_halo(j).EQ.j+c%npeak_cum(g%myid-1).AND.&
             & c%halo_mass(j) > c%mass_threshold.AND. &
             & c%relevance(j) > c%relevance_threshold)then
           mbin=c%mass_bin(j,1:nbin)
           rad=2d0*(c%halo_mass(j)/4d0/3.1415926*3d0/200d0)**(1d0/3d0)
           call halo_mass_def(s,mbin,rad,r200b,rvir,concentration)
           write(ilun,'(I10,X,I10,11(X,1PE18.9E2))')&
                j+c%npeak_cum(g%myid-1)&
                ,c%n_cells_halo(j)&
                ,c%peak_pos(j,1),c%peak_pos(j,2),c%peak_pos(j,3)&
                ,c%peak_vel(j,1),c%peak_vel(j,2),c%peak_vel(j,3)&
                ,c%max_dens(j)&
                ,c%halo_mass(j)&
                ,r200b&
                ,rvir&
                ,concentration
        endif
     enddo

     call close_file(s,fileloc,nskip,ilun)

  end if

  end associate

end subroutine output_clump_properties
!###################################################
!###################################################
!###################################################
!###################################################
subroutine halo_mass_def(s,mbin,rad,r200b,rvir,c)
  use amr_parameters, only: nbin, dp
  use ramses_commons, only: ramses_t
  type(ramses_t)::s
  real(dp),dimension(1:nbin)::mbin
  real(dp)::rad,r200b,rvir,c
  real(dp)::pi,G,x,delta,deltaold,deltavir,rbin,vcirc,volbin,alpha
  real(dp)::vmax,vvir,c0,cl,cr,err,const
  integer::i
  ! Constants
  pi=ACOS(-1.0D0)
  G=1d0
  if(s%r%cosmo)G=3d0/8d0/pi*s%g%omega_m*s%g%aexp
  x=s%g%omega_m/(s%g%omega_m+s%g%omega_l*s%g%aexp**3)-1
  ! Bryan & Norman (1998)
  deltavir=(18d0*pi**2+82d0*x-39d0*x*x)/s%g%omega_m
  ! Initializations
  vmax=0
  r200b=0
  rvir=0
  delta=0
  ! Loop over radial bins
  do i=1,nbin
     if(mbin(i)==0)cycle
     rbin=dble(i)/dble(nbin)*rad
     vcirc=sqrt(G*mbin(i)/rbin)
     if(vcirc>vmax)then
        vmax=vcirc
     endif
     volbin=4d0/3d0*pi*rbin**3
     deltaold=delta
     delta=mbin(i)/volbin
     if(delta<=200.and.r200b==0)then
        if(deltaold>0)then
           alpha=log(200d0/delta)/log(deltaold/delta)
           r200b=rbin*(dble(i-1)/dble(i))**alpha
        else
           r200b=rbin
        endif
     endif
     if(delta<=deltavir.and.rvir==0)then
        if(deltaold>0)then
           alpha=log(deltavir/delta)/log(deltaold/delta)
           rvir=rbin*(dble(i-1)/dble(i))**alpha
        else
           rvir=rbin
        endif
     endif
  end do
  if(delta>200d0.and.r200b==0)then
     alpha=log(200d0/delta)/log(deltaold/delta)
     r200b=rad*(dble(nbin-1)/dble(nbin))**alpha
  endif
  if(delta>deltavir.and.rvir==0)then
     alpha=log(deltavir/delta)/log(deltaold/delta)
     rvir=rad*(dble(nbin-1)/dble(nbin))**alpha
  endif
  c0=2.1626
  vvir=sqrt(G*4d0*pi/3d0*deltavir*rvir**2)
  vmax=max(vmax,vvir)
  const=vmax**2/vvir**2*c0/(log(1d0+c0)-c0/(1d0+c0))
  cl=c0
  cr=200d0
  err=1d0
  i=0
  do while(abs(err)>1d-3.and.i<100)
     c=0.5*(cl+cr)
     err=c/(log(1d0+c)-c/(1d0+c))-const
     if(err>0)then
        cr=c
     else
        cl=c
     endif
     i=i+1
  end do
end subroutine halo_mass_def
!###################################################
!###################################################
!###################################################
!###################################################
subroutine output_clump_field(s,filename)
  use amr_parameters, only: ndim,twotondim,flen
  use ramses_commons, only: ramses_t,open_file,close_file
  implicit none
  type(ramses_t)::s
  character(LEN=flen)::filename
  !----------------------------------------------------------------
  ! This routine output the peak patch fields for each processor
  !----------------------------------------------------------------  
  integer::ilevel,igrid,ilun
  integer(kind=8),dimension(s%r%levelmin:s%r%nlevelmax)::nskip
  real(kind=4),dimension(1:twotondim)::flg
  real(kind=4),dimension(1:twotondim)::rho

  associate(g=>s%g,r=>s%r,m=>s%m,mdl=>s%mdl)

#ifdef GRAV

  call open_file(s,filename,nskip,ilun)

  do ilevel=r%levelmin,r%nlevelmax
     write(ilun,POS=nskip(ilevel))
     do igrid=m%head(ilevel),m%tail(ilevel)
        rho=real(m%grid(igrid)%rho,kind=4)
        flg=real(m%grid(igrid)%flag1,kind=4)
        write(ilun)flg
        write(ilun)rho
     end do
  enddo

  call close_file(s,filename,nskip,ilun)

#endif

  end associate

end subroutine output_clump_field
!###################################################
!###################################################
!###################################################
!###################################################
subroutine file_descriptor_clump(r,filename)
  use amr_parameters, only: ndim,flen
  use amr_commons, only: run_t
  implicit none
  type(run_t)::r
  character(LEN=flen)::filename
  
  character(LEN=flen)::fileloc
  integer::ivar,ilun

  if(r%verbose)write(*,*)'Entering file_descriptor_clump'

  ilun=11
  fileloc=TRIM(filename)
  open(unit=ilun,file=fileloc,form='formatted')

  write(ilun,'("nvar        =",I11)')2

  ivar=1
  write(ilun,'("variable #",I2,": peak ID")')ivar

  ivar=2
  write(ilun,'("variable #",I2,": density")')ivar
  
  close(ilun)

end subroutine file_descriptor_clump
!#######################################################
!#######################################################
!#######################################################
!#######################################################
subroutine output_peak_part(s,p,filename)
  use amr_parameters, only: ndim,i8b,flen
  use ramses_commons, only: ramses_t,open_part_file,close_part_file
  use pm_commons, only: part_t
  implicit none
  type(ramses_t)::s
  type(part_t)::p
  character(LEN=flen)::filename
  !--------------------------------------
  ! Output peak id of particles to file
  !--------------------------------------
  integer::i,idim,ilun,ivar
  integer(kind=8),dimension(1:p%nvaralloc+1)::nskip
  integer,allocatable,dimension(:)::ll

  associate(r=>s%r,g=>s%g)

  call open_part_file(s,p,filename,nskip,ilun)

  ! Write halo id
  allocate(ll(1:p%npart))
  do i=1,p%npart
     ll(i)=p%workp(i)
  end do
  write(ilun,POS=nskip(1))
  write(ilun)ll
  deallocate(ll)

  call close_part_file(s,p,filename,nskip,ilun)

  end associate

end subroutine output_peak_part
!#########################################################################
!#########################################################################
!#########################################################################
!#########################################################################
subroutine output_peak_header(r,g,p,filename)
  use amr_parameters, only: flen
  use amr_commons, only: run_t,global_t
  use pm_commons, only: part_t
  implicit none
  type(run_t)::r
  type(global_t)::g
  type(part_t)::p
  character(LEN=flen)::filename

  ! Local variables
  integer::ilun
  character(LEN=flen)::fileloc

  if(r%verbose)write(*,*)'Entering output_peak_header'

  ilun=10!+g%myid

  ! Open file
  fileloc=TRIM(filename)
  open(unit=ilun,file=fileloc,form='formatted')

  ! Write header information
  write(ilun,*)'Total number of particles'
  write(ilun,*)p%npart_tot

  write(ilun,*)'Total number of files'
  write(ilun,*)r%nfile

  ! Keep track of what particle fields are present
  write(ilun,*)'Particle fields'
  write(ilun,'(a)',advance='no')'peak_id'
  close(ilun)

end subroutine output_peak_header
!#######################################################
!#######################################################
!#######################################################
!#######################################################
end module output_clump_module
