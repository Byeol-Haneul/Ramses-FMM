subroutine init_hydro(r)
  use amr_commons, only: run_t
  use hydro_commons, only: hydro_w
  implicit none
#ifndef WITHOUTMPI
  include 'mpif.h'
#endif
  type(run_t)::r
  
  if(r%verbose)write(*,*)'Entering init_hydro'

  ! Initialise workspace for hydro kernels
  associate(h => hydro_w)
    call h%kernel_1%init(2)
    if(r%nsuperoct>0) call h%kernel_2%init(4)
    if(r%nsuperoct>1) call h%kernel_4%init(8)
    if(r%nsuperoct>2) call h%kernel_8%init(16)
    if(r%nsuperoct>3) call h%kernel_16%init(32)
    if(r%nsuperoct>4) call h%kernel_32%init(64)
  end associate

end subroutine init_hydro




