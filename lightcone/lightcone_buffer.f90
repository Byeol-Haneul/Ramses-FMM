module lightcone_buffer_module
  use amr_parameters, only: sp
  implicit none

  type :: lightcone_buffer
     real(sp), allocatable :: xp(:,:)
     real(sp), allocatable :: vp(:,:)
     integer :: nstride
     integer :: ncurrent
  end type lightcone_buffer

contains

  subroutine init_lightcone_buffer(buffer, nstride)
    type(lightcone_buffer), intent(inout) :: buffer
    integer, intent(in) :: nstride

    buffer%nstride = nstride
    buffer%ncurrent = 0
    allocate(buffer%xp(nstride, 3))
    allocate(buffer%vp(nstride, 3))
  end subroutine init_lightcone_buffer

  subroutine add_to_buffer(buffer, position, velocity)
    type(lightcone_buffer), intent(inout) :: buffer
    real(sp), intent(in) :: position(3)
    real(sp), intent(in) :: velocity(3)

    if (buffer%ncurrent < buffer%nstride) then
       buffer%ncurrent = buffer%ncurrent + 1
       buffer%xp(buffer%ncurrent, :) = position(:)
       buffer%vp(buffer%ncurrent, :) = velocity(:)
    else
       stop 'add_to_buffer: Cannot add particle to a full buffer'
    end if
  end subroutine add_to_buffer

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

end module lightcone_buffer_module
