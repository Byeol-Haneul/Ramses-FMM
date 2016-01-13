subroutine init_hydro
  use amr_commons
  use hydro_commons
  implicit none
#ifndef WITHOUTMPI
  include 'mpif.h'
#endif
  character(LEN=80)::fileloc
  character(LEN=5)::nchar
  integer::nx

  if(verbose)write(*,*)'Entering init_hydro'

  ! Allocate work space for hydro kernel
  nx=2
  iu1=-1; iu2=nx+2
  ju1=(1-ndim/2)-1*(ndim/2); ju2=(1-ndim/2)+(nx+2)*(ndim/2)
  ku1=(1-ndim/3)-1*(ndim/3); ku2=(1-ndim/3)+(nx+2)*(ndim/3)
  if1=1; if2=nx+1
  jf1=1; jf2=(1-ndim/2)+(nx+1)*(ndim/2)
  kf1=1; kf2=(1-ndim/3)+(nx+1)*(ndim/3)

  allocate(uloc(iu1:iu2,ju1:ju2,ku1:ku2,1:nvar))
  allocate(gloc(iu1:iu2,ju1:ju2,ku1:ku2,1:ndim))
  allocate(qloc(iu1:iu2,ju1:ju2,ku1:ku2,1:nvar))
  allocate(cloc(iu1:iu2,ju1:ju2,ku1:ku2))
  allocate(flux(if1:if2,jf1:jf2,kf1:kf2,1:nvar,1:ndim))
  allocate(tmp (if1:if2,jf1:jf2,kf1:kf2,1:2   ,1:ndim))
  allocate(dq  (iu1:iu2,ju1:ju2,ku1:ku2,1:nvar,1:ndim))
  allocate(qm  (iu1:iu2,ju1:ju2,ku1:ku2,1:nvar,1:ndim))
  allocate(qp  (iu1:iu2,ju1:ju2,ku1:ku2,1:nvar,1:ndim))
  allocate(fx  (iu1:iu2,ju1:ju2,ku1:ku2,1:nvar))
  allocate(tx  (iu1:iu2,ju1:ju2,ku1:ku2,1:2   ))
  allocate(divu(if1:if2,jf1:jf2,kf1:kf2))
  allocate(okoc(iu1:iu2,ju1:ju2,ku1:ku2))

end subroutine init_hydro




