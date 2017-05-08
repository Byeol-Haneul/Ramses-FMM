subroutine init_hydro
  use amr_commons
  use hydro_commons
  implicit none
#ifndef WITHOUTMPI
  include 'mpif.h'
#endif

  if(verbose)write(*,*)'Entering init_hydro'

  associate(h => hydro_w)
    call h%kernel_1%init(2)
    if(nsuperoct>0) call h%kernel_2%init(4)
    if(nsuperoct>1) call h%kernel_4%init(8)
    if(nsuperoct>2) call h%kernel_8%init(16)
    if(nsuperoct>3) call h%kernel_16%init(32)
    if(nsuperoct>4) call h%kernel_32%init(64)
  end associate

end subroutine init_hydro




