!#########################################################################
!#########################################################################
!#########################################################################
!#########################################################################
subroutine dump_all
  use amr_commons
  use pm_commons
  use hydro_commons
  implicit none
#ifndef WITHOUTMPI
  include 'mpif.h'
#endif
  character(LEN=5)::nchar
  character(LEN=80)::filename,filedir,filecmd
  integer::i,itest,info,irec,ierr

  if(nstep_coarse==nstep_coarse_old.and.nstep_coarse>0)return
  if(nstep_coarse==0.and.nrestart>0)return
  if(verbose)write(*,*)'Entering dump_all'

  do i=levelmin,nlevelmax
     call write_screen(i)
  end do

  call title(ifout,nchar)
  ifout=ifout+1
  if(t>=tout(iout).or.aexp>=aout(iout))iout=iout+1
  output_done=.true.

  if(ndim>1)then
     filedir='output_'//TRIM(nchar)//'/'
     filecmd='mkdir -p '//TRIM(filedir)
#ifdef NOSYSTEM
     call PXFMKDIR(TRIM(filedir),LEN(TRIM(filedir)),O'755',info)
#else
     call system(filecmd)
#endif
#ifndef WITHOUTMPI
     call MPI_BARRIER(MPI_COMM_WORLD,info)
#endif
     ! Only master process
     if(myid==1)then
        if(pic)then
           filename=TRIM(filedir)//'part_header.txt'
           call output_header(filename)
        endif
        if(hydro)then
           filename=TRIM(filedir)//'hydro_file_descriptor.txt'
           call file_descriptor_hydro(filename)
        end if
#ifdef TOTO
        if(cooling)then
           filename=TRIM(filedir)//'cooling.out'
           call output_cool(filename)
        end if
        if(sink)then
           filename=TRIM(filedir)//'sink.info'
           call output_sink(filename)
           filename=TRIM(filedir)//'sink.csv'
           call output_sink_csv(filename)
        endif
#endif
        filename=TRIM(filedir)//'info.txt'
        call output_info(filename)
        filename=TRIM(filedir)//'makefile.txt'
        call output_makefile(filename)
        filename=TRIM(filedir)//'patches.txt'
        call output_patch(filename)
        filename=TRIM(filedir)//'namelist.txt'
        call output_namelist(filename)
        filename=TRIM(filedir)//'compilation.txt'
        call output_compil(filename)
        filename=TRIM(filedir)//'params.out'
        call output_params(filename)
     endif
     ! For each process
     filename=TRIM(filedir)//'amr.out'
     call output_amr(filename)
     if(hydro)then
        filename=TRIM(filedir)//'hydro.out'
        call output_hydro(filename)
     end if
     if(poisson)then
        filename=TRIM(filedir)//'grav.out'
        call output_poisson(filename)
     end if
#ifdef TOTO
     if(pic)then
        filename=TRIM(filedir)//'part.out'
        call output_part(filename)
     end if
     if (gadget_output) then
        filename=TRIM(filedir)//'gsnapshot_'//TRIM(nchar)
        call savegadget(filename)
     end if
#endif
  end if

end subroutine dump_all
!#########################################################################
!#########################################################################
!#########################################################################
!#########################################################################
subroutine output_namelist(filename)
  use amr_commons
  use pm_commons
  use hydro_commons
  implicit none
  character(LEN=80)::filename
  ! Copy namelist file to output directory
  character::nml_char
  integer::ierr

  open(10,file=namelist_file,access="stream",action="read")
  open(11,file=filename,access="stream",action="write")
  do 
     read(10,iostat=ierr)nml_char 
     if(ierr.NE.0)exit
     write(11)nml_char 
  end do
  close(10)
  close(11)

end subroutine output_namelist
!#########################################################################
!#########################################################################
!#########################################################################
!#########################################################################
subroutine output_compil(filename)
  use amr_commons
  use pm_commons
  use hydro_commons
  implicit none
  character(LEN=80)::filename
  ! Copy compilation details to output directory
  OPEN(UNIT=11, FILE=filename, FORM='formatted')
  write(11,'(" compile date = ",A)')TRIM(builddate)
  write(11,'(" patch dir    = ",A)')TRIM(patchdir)
  write(11,'(" remote repo  = ",A)')TRIM(gitrepo)
  write(11,'(" local branch = ",A)')TRIM(gitbranch)
  write(11,'(" last commit  = ",A)')TRIM(githash)
  CLOSE(11)
end subroutine output_compil
!#########################################################################
!#########################################################################
!#########################################################################
!#########################################################################
subroutine output_params(filename)
  use amr_commons
  use hydro_commons
  use pm_commons
  implicit none
  character(LEN=80)::filename

  integer::ilun
  integer::ilevel,ibound,istart,i,igrid,idim,ind,iskip
  character(LEN=80)::fileloc

  if(verbose)write(*,*)'Entering output_params'

  !-----------------------------------
  ! Output run parameters in file
  !-----------------------------------  
  ilun=10
  fileloc=TRIM(filename)
  open(unit=ilun,file=fileloc,access="stream"&
       & ,action="write",form='unformatted')
  ! Write grid variables
  write(ilun)ncpu
  write(ilun)ndim
  write(ilun)levelmin
  write(ilun)nlevelmax
  write(ilun)boxlen
  ! Write time variables
  write(ilun)noutput,iout,ifout
  write(ilun)tout(1:noutput)
  write(ilun)aout(1:noutput)
  write(ilun)t
  write(ilun)dtold(1:nlevelmax)
  write(ilun)dtnew(1:nlevelmax)
  write(ilun)nstep,nstep_coarse
  ! Write various constants
  write(ilun)const,mass_tot_0,rho_tot
  write(ilun)omega_m,omega_l,omega_k,omega_b,h0,aexp_ini,boxlen_ini
  write(ilun)aexp,hexp,aexp_old,epot_tot_int,epot_tot_old
  write(ilun)mass_sph
  ! Write cpu boundaries
  write(ilun)nhilbert
  do ilevel=levelmin,nlevelmax
     write(ilun)bound_key_level(1:nhilbert,0:ncpu,ilevel)
  end do
  close(ilun)

end subroutine output_params
!#########################################################################
!#########################################################################
!#########################################################################
!#########################################################################
subroutine input_params(filename,ncpu_file,levelmin_file,nlevelmax_file)
  use amr_commons
  use hydro_commons
  use pm_commons
  implicit none
  character(LEN=80)::filename
  integer::ncpu_file,levelmin_file,nlevelmax_file,i
  !-----------------------------------
  ! Read run parameters from file.
  ! Note that ncpu, levelmin and nlevelmax
  ! are allowed to vary at restart.
  !-----------------------------------  
  integer::ilun
  integer::ndim_file,noutput_file
  integer::noutput_min,nlevelmax_min
  real(dp)::mass_sph_file
  character(LEN=80)::fileloc

  if(verbose)write(*,*)'Entering input_params'

  ilun=10+myid
  fileloc=TRIM(filename)
  open(unit=ilun,file=fileloc,access="stream"&
       & ,action="read",form='unformatted')
  ! Read grid variables
  read(ilun)ncpu_file
  read(ilun)ndim_file
  read(ilun)levelmin_file
  read(ilun)nlevelmax_file
  ! Overwrite boxlen with value from file
  read(ilun)boxlen
  ! Read time variables
  read(ilun)noutput_file,iout,ifout
  noutput_min=MIN(noutput,noutput_file)
  read(ilun)tout(1:noutput_min)
  read(ilun)aout(1:noutput_min)
  read(ilun)t
  nlevelmax_min=MIN(nlevelmax,nlevelmax_file)
  read(ilun)dtold(1:nlevelmax_min)
  read(ilun)dtnew(1:nlevelmax_min)
  read(ilun)nstep,nstep_coarse
  ! Read various constants
  read(ilun)const,mass_tot_0,rho_tot
  read(ilun)omega_m,omega_l,omega_k,omega_b,h0,aexp_ini,boxlen_ini
  read(ilun)aexp,hexp,aexp_old,epot_tot_int,epot_tot_old
  read(ilun)mass_sph_file
  close(ilun)
  ! For cosmo runs only, as mass_sph is not set in the namelist
  if(cosmo)mass_sph=mass_sph_file
  ! Check dimensions
  if(ndim.NE.ndim_file)then
     if(myid==1)then
        write(*,*)'Incorrect number of space dimensions in restart file'
     endif
     call clean_stop
  endif
  ! Compute movie frame number if applicable
  if(imovout>0) then
     do i=2,imovout
        if(aendmov>0)then
           if(aexp>amovout(i-1).and.aexp<amovout(i)) then
              imov=i
           endif
        else
           if(t>tmovout(i-1).and.t<tmovout(i)) then
              imov=i
           endif
        endif
     enddo
  endif

end subroutine input_params
!#########################################################################
!#########################################################################
!#########################################################################
!#########################################################################
subroutine output_amr(filename)
  use amr_commons
  use hydro_commons
  use pm_commons
  implicit none
  character(LEN=80)::filename
  !-----------------------------------
  ! Output amr grid in file
  !-----------------------------------  
  integer::ilun,mypos
  integer::ilevel,ibound,istart,i,igrid,idim,ind,iskip
  integer,allocatable,dimension(:)::ind_grid,iig
  real(dp),allocatable,dimension(:)::xdp
  real(sp),allocatable,dimension(:)::xsp
  real(dp),dimension(1:3)::skip_loc
  character(LEN=80)::fileloc
  character(LEN=5)::nchar
  real(dp)::scale_nH,scale_T2,scale_l,scale_d,scale_t,scale_v
  real(dp)::scale
  if(verbose)write(*,*)'Entering output_amr'
  ilun=myid+10
  call title(myid,nchar)
  fileloc=TRIM(filename)//TRIM(nchar)
  open(unit=ilun,file=fileloc,access="stream"&
       & ,action="write",form='unformatted')
  write(ilun)ndim
  write(ilun)levelmin
  write(ilun)nlevelmax
  do ilevel=levelmin,nlevelmax
     write(ilun)noct(ilevel)
  end do
  do ilevel=levelmin,nlevelmax
     do igrid=head(ilevel),tail(ilevel)
        write(ilun)grid(igrid)%ckey
        write(ilun)grid(igrid)%refined
     end do
  end do
  close(ilun)  
end subroutine output_amr
!#########################################################################
!#########################################################################
!#########################################################################
!#########################################################################
subroutine output_info(filename)
  use amr_commons
  use hydro_commons
  use pm_commons
  implicit none
  character(LEN=80)::filename

  integer::ilun,icpu,idom
  real(dp)::scale_nH,scale_T2,scale_l,scale_d,scale_t,scale_v
  character(LEN=80)::fileloc
  character(LEN=5)::nchar

  if(verbose)write(*,*)'Entering output_info'

  ilun=11

  ! Conversion factor from user units to cgs units
  call units(scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)

  ! Open file
  fileloc=TRIM(filename)
  open(unit=ilun,file=fileloc,form='formatted')
  
  ! Write run parameters
  write(ilun,'("ncpu        =",I11)')ncpu
  write(ilun,'("ndim        =",I11)')ndim
  write(ilun,'("levelmin    =",I11)')levelmin
  write(ilun,'("levelmax    =",I11)')nlevelmax
  write(ilun,'("ngridmax    =",I11)')ngridmax
  write(ilun,'("nstep_coarse=",I11)')nstep_coarse
  write(ilun,*)

  ! Write physical parameters
  write(ilun,'("boxlen      =",E23.15)')boxlen
  write(ilun,'("time        =",E23.15)')t
  write(ilun,'("aexp        =",E23.15)')aexp
  write(ilun,'("H0          =",E23.15)')h0
  write(ilun,'("omega_m     =",E23.15)')omega_m
  write(ilun,'("omega_l     =",E23.15)')omega_l
  write(ilun,'("omega_k     =",E23.15)')omega_k
  write(ilun,'("omega_b     =",E23.15)')omega_b
  write(ilun,'("unit_l      =",E23.15)')scale_l
  write(ilun,'("unit_d      =",E23.15)')scale_d
  write(ilun,'("unit_t      =",E23.15)')scale_t
  write(ilun,*)
  
  close(ilun)

end subroutine output_info
!#########################################################################
!#########################################################################
!#########################################################################
!#########################################################################
subroutine output_header(filename)
  use amr_commons
  use hydro_commons
  use pm_commons
  implicit none
#ifndef WITHOUTMPI
  include 'mpif.h'
#endif
  character(LEN=80)::filename

  integer::info,ilun
  integer(i8b)::tmp_long
  character(LEN=80)::fileloc

  if(verbose)write(*,*)'Entering output_header'
  
  ilun=myid+10
  
  ! Open file
  fileloc=TRIM(filename)
  open(unit=ilun,file=fileloc,form='formatted')
  
  ! Write header information
  write(ilun,*)'Total number of particles'
  write(ilun,*)npart
  
  ! Keep track of what particle fields are present
  write(ilun,*)'Particle fields'
  write(ilun,'(a)',advance='no')'pos vel mass iord level '
#ifdef OUTPUT_PARTICLE_POTENTIAL
  write(ilun,'(a)',advance='no')'phi '
#endif
  close(ilun)

end subroutine output_header
!#########################################################################
!#########################################################################
!#########################################################################
!#########################################################################
#ifdef TOTO
subroutine savegadget(filename)
  use amr_commons
  use hydro_commons
  use pm_commons
  use gadgetreadfilemod
  implicit none
#ifndef WITHOUTMPI
  include 'mpif.h'
#endif
  character(LEN=80)::filename
  TYPE (gadgetheadertype) :: header
  real,allocatable,dimension(:,:)::pos, vel
  integer(i8b),allocatable,dimension(:)::ids
  integer::i, idim, ipart
  real:: gadgetvfact
  integer::info
  integer(i8b)::npart_tot, npart_loc
  real, parameter:: RHOcrit = 2.7755d11

#ifndef WITHOUTMPI
  npart_loc=npart
#ifndef LONGINT
  call MPI_ALLREDUCE(npart_loc,npart_tot,1,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,info)
#else
  call MPI_ALLREDUCE(npart_loc,npart_tot,1,MPI_INTEGER8,MPI_SUM,MPI_COMM_WORLD,info)
#endif
#else
  npart_tot=npart
#endif

  allocate(pos(ndim, npart), vel(ndim, npart), ids(npart))
  gadgetvfact = 100.0 * boxlen_ini / aexp / SQRT(aexp)

  header%npart = 0
  header%npart(2) = npart
  header%mass = 0
  header%mass(2) = omega_m*RHOcrit*(boxlen_ini)**3/npart_tot/1.d10
  header%time = aexp
  header%redshift = 1.d0/aexp-1.d0
  header%flag_sfr = 0
  header%nparttotal = 0
#ifndef LONGINT
  header%nparttotal(2) = npart_tot
#else
  header%nparttotal(2) = MOD(npart_tot,4294967296)
#endif
  header%flag_cooling = 0
  header%numfiles = ncpu
  header%boxsize = boxlen_ini
  header%omega0 = omega_m
  header%omegalambda = omega_l
  header%hubbleparam = h0/100.0
  header%flag_stellarage = 0
  header%flag_metals = 0
  header%totalhighword = 0
#ifndef LONGINT
  header%totalhighword(2) = 0
#else
  header%totalhighword(2) = npart_tot/4294967296
#endif
  header%flag_entropy_instead_u = 0
  header%flag_doubleprecision = 0
  header%flag_ic_info = 0
  header%lpt_scalingfactor = 0
  header%unused = ' '

  do idim=1,ndim
     ipart=0
     do i=1,npartmax
        if(levelp(i)>0)then
           ipart=ipart+1
           if (ipart .gt. npart) then
                write(*,*) myid, "Ipart=",ipart, "exceeds", npart
                call clean_stop
           endif
           pos(idim, ipart)=xp(i,idim) * boxlen_ini
           vel(idim, ipart)=vp(i,idim) * gadgetvfact
           if (idim.eq.1) ids(ipart) = idp(i)
        end if
     end do
  end do

  call gadgetwritefile(filename, myid-1, header, pos, vel, ids)
  deallocate(pos, vel, ids)

end subroutine savegadget
#endif
