subroutine init_hydro
  use amr_commons
  use hydro_commons
  implicit none
#ifndef WITHOUTMPI
  include 'mpif.h'
#endif
  character(LEN=80)::fileloc
  character(LEN=5)::nchar

  if(verbose)write(*,*)'Entering init_hydro'

  ! Allocate work space for hydro kernel
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
  allocate(okloc(iu1:iu2,ju1:ju2,ku1:ku2))
  allocate(childloc(io1:io2,jo1:jo2,ko1:ko2))
  allocate(parentloc(io1:io2,jo1:jo2,ko1:ko2))

  if(nsuperoct>0)then
  ! Allocate work space for hydro kernel
  allocate(uloc_2(iu1_2:iu2_2,ju1_2:ju2_2,ku1_2:ku2_2,1:nvar))
  allocate(gloc_2(iu1_2:iu2_2,ju1_2:ju2_2,ku1_2:ku2_2,1:ndim))
  allocate(qloc_2(iu1_2:iu2_2,ju1_2:ju2_2,ku1_2:ku2_2,1:nvar))
  allocate(cloc_2(iu1_2:iu2_2,ju1_2:ju2_2,ku1_2:ku2_2))
  allocate(flux_2(if1_2:if2_2,jf1_2:jf2_2,kf1_2:kf2_2,1:nvar,1:ndim))
  allocate(tmp_2 (if1_2:if2_2,jf1_2:jf2_2,kf1_2:kf2_2,1:2   ,1:ndim))
  allocate(dq_2  (iu1_2:iu2_2,ju1_2:ju2_2,ku1_2:ku2_2,1:nvar,1:ndim))
  allocate(qm_2  (iu1_2:iu2_2,ju1_2:ju2_2,ku1_2:ku2_2,1:nvar,1:ndim))
  allocate(qp_2  (iu1_2:iu2_2,ju1_2:ju2_2,ku1_2:ku2_2,1:nvar,1:ndim))
  allocate(fx_2  (iu1_2:iu2_2,ju1_2:ju2_2,ku1_2:ku2_2,1:nvar))
  allocate(tx_2  (iu1_2:iu2_2,ju1_2:ju2_2,ku1_2:ku2_2,1:2   ))
  allocate(divu_2(if1_2:if2_2,jf1_2:jf2_2,kf1_2:kf2_2))
  allocate(okloc_2    (iu1_2:iu2_2,ju1_2:ju2_2,ku1_2:ku2_2))
  allocate(childloc_2 (io1_2:io2_2,jo1_2:jo2_2,ko1_2:ko2_2))
  allocate(parentloc_2(io1_2:io2_2,jo1_2:jo2_2,ko1_2:ko2_2))
  endif

  if(nsuperoct>1)then
  ! Allocate work space for hydro kernel
  allocate(uloc_4(iu1_4:iu2_4,ju1_4:ju2_4,ku1_4:ku2_4,1:nvar))
  allocate(gloc_4(iu1_4:iu2_4,ju1_4:ju2_4,ku1_4:ku2_4,1:ndim))
  allocate(qloc_4(iu1_4:iu2_4,ju1_4:ju2_4,ku1_4:ku2_4,1:nvar))
  allocate(cloc_4(iu1_4:iu2_4,ju1_4:ju2_4,ku1_4:ku2_4))
  allocate(flux_4(if1_4:if2_4,jf1_4:jf2_4,kf1_4:kf2_4,1:nvar,1:ndim))
  allocate(tmp_4 (if1_4:if2_4,jf1_4:jf2_4,kf1_4:kf2_4,1:2   ,1:ndim))
  allocate(dq_4  (iu1_4:iu2_4,ju1_4:ju2_4,ku1_4:ku2_4,1:nvar,1:ndim))
  allocate(qm_4  (iu1_4:iu2_4,ju1_4:ju2_4,ku1_4:ku2_4,1:nvar,1:ndim))
  allocate(qp_4  (iu1_4:iu2_4,ju1_4:ju2_4,ku1_4:ku2_4,1:nvar,1:ndim))
  allocate(fx_4  (iu1_4:iu2_4,ju1_4:ju2_4,ku1_4:ku2_4,1:nvar))
  allocate(tx_4  (iu1_4:iu2_4,ju1_4:ju2_4,ku1_4:ku2_4,1:2   ))
  allocate(divu_4(if1_4:if2_4,jf1_4:jf2_4,kf1_4:kf2_4))
  allocate(okloc_4    (iu1_4:iu2_4,ju1_4:ju2_4,ku1_4:ku2_4))
  allocate(childloc_4 (io1_4:io2_4,jo1_4:jo2_4,ko1_4:ko2_4))
  allocate(parentloc_4(io1_4:io2_4,jo1_4:jo2_4,ko1_4:ko2_4))
  endif

  if(nsuperoct>2)then
  ! Allocate work space for hydro kernel
  allocate(uloc_8(iu1_8:iu2_8,ju1_8:ju2_8,ku1_8:ku2_8,1:nvar))
  allocate(gloc_8(iu1_8:iu2_8,ju1_8:ju2_8,ku1_8:ku2_8,1:ndim))
  allocate(qloc_8(iu1_8:iu2_8,ju1_8:ju2_8,ku1_8:ku2_8,1:nvar))
  allocate(cloc_8(iu1_8:iu2_8,ju1_8:ju2_8,ku1_8:ku2_8))
  allocate(flux_8(if1_8:if2_8,jf1_8:jf2_8,kf1_8:kf2_8,1:nvar,1:ndim))
  allocate(tmp_8 (if1_8:if2_8,jf1_8:jf2_8,kf1_8:kf2_8,1:2   ,1:ndim))
  allocate(dq_8  (iu1_8:iu2_8,ju1_8:ju2_8,ku1_8:ku2_8,1:nvar,1:ndim))
  allocate(qm_8  (iu1_8:iu2_8,ju1_8:ju2_8,ku1_8:ku2_8,1:nvar,1:ndim))
  allocate(qp_8  (iu1_8:iu2_8,ju1_8:ju2_8,ku1_8:ku2_8,1:nvar,1:ndim))
  allocate(fx_8  (iu1_8:iu2_8,ju1_8:ju2_8,ku1_8:ku2_8,1:nvar))
  allocate(tx_8  (iu1_8:iu2_8,ju1_8:ju2_8,ku1_8:ku2_8,1:2   ))
  allocate(divu_8(if1_8:if2_8,jf1_8:jf2_8,kf1_8:kf2_8))
  allocate(okloc_8    (iu1_8:iu2_8,ju1_8:ju2_8,ku1_8:ku2_8))
  allocate(childloc_8 (io1_8:io2_8,jo1_8:jo2_8,ko1_8:ko2_8))
  allocate(parentloc_8(io1_8:io2_8,jo1_8:jo2_8,ko1_8:ko2_8))
  endif

  if(nsuperoct>3)then
  ! Allocate work space for hydro kernel
  allocate(uloc_16(iu1_16:iu2_16,ju1_16:ju2_16,ku1_16:ku2_16,1:nvar))
  allocate(gloc_16(iu1_16:iu2_16,ju1_16:ju2_16,ku1_16:ku2_16,1:ndim))
  allocate(qloc_16(iu1_16:iu2_16,ju1_16:ju2_16,ku1_16:ku2_16,1:nvar))
  allocate(cloc_16(iu1_16:iu2_16,ju1_16:ju2_16,ku1_16:ku2_16))
  allocate(flux_16(if1_16:if2_16,jf1_16:jf2_16,kf1_16:kf2_16,1:nvar,1:ndim))
  allocate(tmp_16 (if1_16:if2_16,jf1_16:jf2_16,kf1_16:kf2_16,1:2   ,1:ndim))
  allocate(dq_16  (iu1_16:iu2_16,ju1_16:ju2_16,ku1_16:ku2_16,1:nvar,1:ndim))
  allocate(qm_16  (iu1_16:iu2_16,ju1_16:ju2_16,ku1_16:ku2_16,1:nvar,1:ndim))
  allocate(qp_16  (iu1_16:iu2_16,ju1_16:ju2_16,ku1_16:ku2_16,1:nvar,1:ndim))
  allocate(fx_16  (iu1_16:iu2_16,ju1_16:ju2_16,ku1_16:ku2_16,1:nvar))
  allocate(tx_16  (iu1_16:iu2_16,ju1_16:ju2_16,ku1_16:ku2_16,1:2   ))
  allocate(divu_16(if1_16:if2_16,jf1_16:jf2_16,kf1_16:kf2_16))
  allocate(okloc_16    (iu1_16:iu2_16,ju1_16:ju2_16,ku1_16:ku2_16))
  allocate(childloc_16 (io1_16:io2_16,jo1_16:jo2_16,ko1_16:ko2_16))
  allocate(parentloc_16(io1_16:io2_16,jo1_16:jo2_16,ko1_16:ko2_16))
  endif

  if(nsuperoct>4)then
  ! Allocate work space for hydro kernel
  allocate(uloc_32(iu1_32:iu2_32,ju1_32:ju2_32,ku1_32:ku2_32,1:nvar))
  allocate(gloc_32(iu1_32:iu2_32,ju1_32:ju2_32,ku1_32:ku2_32,1:ndim))
  allocate(qloc_32(iu1_32:iu2_32,ju1_32:ju2_32,ku1_32:ku2_32,1:nvar))
  allocate(cloc_32(iu1_32:iu2_32,ju1_32:ju2_32,ku1_32:ku2_32))
  allocate(flux_32(if1_32:if2_32,jf1_32:jf2_32,kf1_32:kf2_32,1:nvar,1:ndim))
  allocate(tmp_32 (if1_32:if2_32,jf1_32:jf2_32,kf1_32:kf2_32,1:2   ,1:ndim))
  allocate(dq_32  (iu1_32:iu2_32,ju1_32:ju2_32,ku1_32:ku2_32,1:nvar,1:ndim))
  allocate(qm_32  (iu1_32:iu2_32,ju1_32:ju2_32,ku1_32:ku2_32,1:nvar,1:ndim))
  allocate(qp_32  (iu1_32:iu2_32,ju1_32:ju2_32,ku1_32:ku2_32,1:nvar,1:ndim))
  allocate(fx_32  (iu1_32:iu2_32,ju1_32:ju2_32,ku1_32:ku2_32,1:nvar))
  allocate(tx_32  (iu1_32:iu2_32,ju1_32:ju2_32,ku1_32:ku2_32,1:2   ))
  allocate(divu_32(if1_32:if2_32,jf1_32:jf2_32,kf1_32:kf2_32))
  allocate(okloc_32    (iu1_32:iu2_32,ju1_32:ju2_32,ku1_32:ku2_32))
  allocate(childloc_32 (io1_32:io2_32,jo1_32:jo2_32,ko1_32:ko2_32))
  allocate(parentloc_32(io1_32:io2_32,jo1_32:jo2_32,ko1_32:ko2_32))
  endif

end subroutine init_hydro




