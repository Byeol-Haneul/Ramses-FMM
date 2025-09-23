module lightcone_buffer_module
  use amr_parameters, only: sp, i8b
  implicit none

  type :: lightcone_buffer
     integer(i8b), allocatable :: idp(:)
     real(sp), allocatable :: xp(:,:)
     real(sp), allocatable :: vp(:,:)
     real(sp), allocatable :: mp(:)
     integer :: nstride
     integer :: ncurrent
  end type lightcone_buffer

  type :: lightcone_buffer_grav
     real(sp), allocatable :: xc(:,:)
     real(sp), allocatable :: rho(:)
     real(sp), allocatable :: phi(:)
     real(sp), allocatable :: accel(:,:)
     real(sp), allocatable :: dphidt(:)
     integer :: nstride
     integer :: ncurrent
  end type lightcone_buffer_grav

contains

  subroutine init_lightcone_buffer(buffer, nstride)
    type(lightcone_buffer), intent(inout) :: buffer
    integer, intent(in) :: nstride

    buffer%nstride = nstride
    buffer%ncurrent = 0
    allocate(buffer%idp(nstride))
    allocate(buffer%xp(nstride, 3))
    allocate(buffer%vp(nstride, 3))
    allocate(buffer%mp(nstride))
  end subroutine init_lightcone_buffer

  subroutine init_lightcone_buffer_grav(buffer, nstride)
    type(lightcone_buffer_grav), intent(inout) :: buffer
    integer, intent(in) :: nstride

    buffer%nstride = nstride
    buffer%ncurrent = 0
    allocate(buffer%xc(nstride, 3))
    allocate(buffer%rho(nstride))
    allocate(buffer%phi(nstride))
    allocate(buffer%accel(nstride, 3))
    allocate(buffer%dphidt(nstride))
  end subroutine init_lightcone_buffer_grav

  subroutine add_to_buffer(buffer, particle_id, position, velocity, mass)
    type(lightcone_buffer), intent(inout) :: buffer
    integer(i8b), intent(in) :: particle_id
    real(sp), intent(in) :: position(3)
    real(sp), intent(in) :: velocity(3)
    real(sp), intent(in) :: mass

    if (buffer%ncurrent < buffer%nstride) then
       buffer%ncurrent = buffer%ncurrent + 1
       buffer%idp(buffer%ncurrent) = particle_id
       buffer%xp(buffer%ncurrent, :) = position(:)
       buffer%vp(buffer%ncurrent, :) = velocity(:)
       buffer%mp(buffer%ncurrent) = mass
    else
       stop 'add_to_buffer: Cannot add particle to a full buffer'
    end if
  end subroutine add_to_buffer

  subroutine add_to_buffer_grav(buffer, position, rho, phi, accel, dphidt)
    type(lightcone_buffer_grav), intent(inout) :: buffer
    real(sp), intent(in) :: position(3)
    real(sp), intent(in) :: rho, phi, dphidt
    real(sp), intent(in) :: accel(3)

    if (buffer%ncurrent < buffer%nstride) then
       buffer%ncurrent = buffer%ncurrent + 1
       buffer%xc(buffer%ncurrent, :) = position(:)
       buffer%rho(buffer%ncurrent) = rho
       buffer%phi(buffer%ncurrent) = phi
       buffer%accel(buffer%ncurrent, :) = accel(:)
       buffer%dphidt(buffer%ncurrent) = dphidt
    else
       stop 'add_to_buffer_grav: Cannot add cell to a full buffer'
    end if
  end subroutine add_to_buffer_grav

  logical function buffer_is_full(buffer)
    type(lightcone_buffer), intent(in) :: buffer

    buffer_is_full = (buffer%ncurrent >= buffer%nstride)
  end function buffer_is_full

  logical function buffer_is_empty(buffer)
    type(lightcone_buffer), intent(in) :: buffer

    buffer_is_empty = (buffer%ncurrent == 0)
  end function buffer_is_empty

  subroutine empty_buffer(buffer)
    type(lightcone_buffer), intent(inout) :: buffer

    buffer%ncurrent = 0
  end subroutine empty_buffer

  logical function buffer_grav_is_full(buffer)
    type(lightcone_buffer_grav), intent(in) :: buffer

    buffer_grav_is_full = (buffer%ncurrent >= buffer%nstride)
  end function buffer_grav_is_full

  logical function buffer_grav_is_empty(buffer)
    type(lightcone_buffer_grav), intent(in) :: buffer

    buffer_grav_is_empty = (buffer%ncurrent == 0)
  end function buffer_grav_is_empty

  subroutine empty_buffer_grav(buffer)
    type(lightcone_buffer_grav), intent(inout) :: buffer

    buffer%ncurrent = 0
  end subroutine empty_buffer_grav

end module lightcone_buffer_module
