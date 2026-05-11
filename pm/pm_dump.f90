! Diagnostic dump hooks for the CPU<->GPU particle-pipeline diff harness.
! Active only when built with -DPART_DUMP.
! Output is written to ${PART_DUMP_DIR} (default: ./part_dump) as
! stream-access unformatted binary files, one per dump call.

module pm_dump
  use amr_parameters, only: ndim, dp, i8b
  use amr_commons, only: mesh_t, run_t
  use pm_commons, only: part_t
  implicit none

  private
  public :: dump_part_state, dump_mesh_state, dump_set_dir

  integer, parameter :: MAGIC_PART = int(z'50415254')  ! 'PART'
  integer, parameter :: MAGIC_MESH = int(z'4D455348')  ! 'MESH'

  character(len=512), save :: dump_dir = ""
  logical, save           :: dir_resolved = .false.

contains

  subroutine resolve_dir()
    integer :: status, l
    if (dir_resolved) return
    call get_environment_variable("PART_DUMP_DIR", dump_dir, length=l, status=status)
    if (status /= 0 .or. l == 0) dump_dir = "./part_dump"
    ! Best-effort mkdir; ignore failures (dir may already exist)
    call execute_command_line("mkdir -p " // trim(dump_dir), wait=.true.)
    dir_resolved = .true.
  end subroutine resolve_dir

  subroutine dump_set_dir(d)
    character(len=*), intent(in) :: d
    dump_dir = d
    dir_resolved = .true.
    call execute_command_line("mkdir -p " // trim(dump_dir), wait=.true.)
  end subroutine dump_set_dir

  subroutine dump_part_state(p, r, tag, ilevel)
    type(part_t),     intent(in) :: p
    type(run_t),      intent(in) :: r
    character(len=*), intent(in) :: tag
    integer,          intent(in) :: ilevel

    character(len=1024) :: path
    integer :: u, npart, nlevels
    integer :: ih0, ih1, it0, it1

    call resolve_dir()
    npart   = p%npart
    nlevels = r%nlevelmax
    write(path,'(a,"/part_",a,"_lev",i2.2,".bin")') trim(dump_dir), trim(tag), ilevel
    open(newunit=u, file=trim(path), access="stream", form="unformatted", status="replace")
    write(u) MAGIC_PART
    write(u) ilevel
    write(u) ndim
    write(u) nlevels
    write(u) npart
    ! headp/tailp are allocated with lower bound levelmin or levelmin-1, not 0.
    if (allocated(p%headp)) then
      ih0 = lbound(p%headp, 1)
      ih1 = ubound(p%headp, 1)
      write(u) ih0, ih1
      write(u) p%headp(ih0:ih1)
    else
      write(u) 0, -1
    end if
    if (allocated(p%tailp)) then
      it0 = lbound(p%tailp, 1)
      it1 = ubound(p%tailp, 1)
      write(u) it0, it1
      write(u) p%tailp(it0:it1)
    else
      write(u) 0, -1
    end if
    ! Force fixed on-disk widths so the diff harness reader (int64 idp, real32 mp)
    ! is invariant to build flags (LONGINT for i8b; dp kind for real arrays).
    ! See gpu_part_diff_harness/README.md "Dump format".
    if (allocated(p%idp))    write(u) int(p%idp(1:npart), kind=8)
    if (allocated(p%levelp)) write(u) p%levelp(1:npart)
    if (allocated(p%sortp))  write(u) p%sortp(1:npart)
    if (allocated(p%xp)) write(u) p%xp(1:npart, 1:ndim)
    if (allocated(p%vp)) write(u) p%vp(1:npart, 1:ndim)
    if (allocated(p%mp)) write(u) real(p%mp(1:npart), kind=4)
    if (allocated(p%fp)) then
       write(u) 1
       write(u) p%fp(1:npart, 1:ndim)
    else
       write(u) 0
    end if
    close(u)
  end subroutine dump_part_state

  subroutine dump_mesh_state(m, tag, ilevel)
    type(mesh_t),     intent(in) :: m
    character(len=*), intent(in) :: tag
    integer,          intent(in) :: ilevel

    character(len=1024) :: path
    integer :: u, ngrid, ncell

    call resolve_dir()
    ngrid = m%ifree - 1
    ncell = 0
    write(path,'(a,"/mesh_",a,"_lev",i2.2,".bin")') trim(dump_dir), trim(tag), ilevel
    open(newunit=u, file=trim(path), access="stream", form="unformatted", status="replace")
    write(u) MAGIC_MESH
    write(u) ilevel
    write(u) ngrid
#ifdef GRAV
    if (allocated(m%rho)) then
       ncell = size(m%rho, 1)
       write(u) 1
       write(u) ncell
       write(u) m%rho(:,1:ngrid)
    else
       write(u) 0
       write(u) 0
    end if
    if (allocated(m%nref)) then
       write(u) 1
       write(u) m%nref(:,1:ngrid)
    else
       write(u) 0
    end if
#else
    write(u) 0
    write(u) 0
    write(u) 0
#endif
    close(u)
  end subroutine dump_mesh_state

end module pm_dump
