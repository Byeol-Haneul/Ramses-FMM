module coordinates
  
  use amr_parameters, only: ndim, nvector, dp
  implicit none

contains
  
  function grid_to_integer_nvector(xgrid, ilevel, n)
    real(dp),        intent(in), dimension(1:n, 1:ndim)  :: xgrid
    integer, intent(in) :: ilevel, n
    integer(kind=8), dimension(1:n, 1:3) :: grid_to_integer_nvector

    integer  :: idim, i
    real(dp) :: fact

    fact = 2.0_dp ** ilevel
    
    do idim = 1, ndim
       do i = 1, n
          grid_to_integer_nvector(i, idim) = floor(fact * xgrid(i, idim), kind=8)
       end do
    end do

  end function grid_to_integer_nvector

  function grid_to_integer(xgrid, ilevel)
    real(dp),        intent(in), dimension(1:ndim)  :: xgrid
    integer, intent(in) :: ilevel
    integer(kind=8), dimension(1:3) :: grid_to_integer
    
    integer  :: idim       
    do idim = 1, ndim
          grid_to_integer(idim) = floor(2.0_dp ** ilevel * xgrid(idim), kind=8)
    end do
    
  end function grid_to_integer

end module coordinates
