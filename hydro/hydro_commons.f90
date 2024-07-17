module hydro_commons
  use amr_parameters
  use hydro_parameters
  use oct_commons

  type hydro_kernel_t
     integer::iu1,iu2,ju1,ju2,ku1,ku2
     integer::if1,if2,jf1,jf2,kf1,kf2
     integer::io1,io2,jo1,jo2,ko1,ko2
     real(dp),dimension(:,:,:,:),allocatable::uloc
     real(dp),dimension(:,:,:,:),allocatable::gloc
     real(dp),dimension(:,:,:,:),allocatable::qloc
     real(dp),dimension(:,:,:),allocatable::cloc
     real(dp),dimension(:,:,:,:,:),allocatable::flux
     real(dp),dimension(:,:,:,:,:),allocatable::tmp
     real(dp),dimension(:,:,:,:,:),allocatable::dq
     real(dp),dimension(:,:,:,:,:),allocatable::qm
     real(dp),dimension(:,:,:,:,:),allocatable::qp
     real(dp),dimension(:,:,:,:),allocatable::fx
     real(dp),dimension(:,:,:,:),allocatable::tx
     real(dp),dimension(:,:,:),allocatable::divu
     logical ,dimension(:,:,:),allocatable::okloc
     integer ,dimension(:,:,:),allocatable::cellloc
     type(nbor),dimension(:,:,:),allocatable::childloc
     type(nbor),dimension(:,:,:),allocatable::gridloc
     type(nbor),dimension(:,:,:,:),allocatable::nborloc
#ifdef MHD
     real(dp),dimension(:,:,:,:),allocatable::bloc
     real(dp),dimension(:,:,:),allocatable::emfx
     real(dp),dimension(:,:,:),allocatable::emfy
     real(dp),dimension(:,:,:),allocatable::emfz
     real(dp),dimension(:,:,:),allocatable::Ex
     real(dp),dimension(:,:,:),allocatable::Ey
     real(dp),dimension(:,:,:),allocatable::Ez
     real(dp),dimension(:,:,:,:),allocatable::bf
     real(dp),dimension(:,:,:,:,:),allocatable::dbf
     real(dp),dimension(:,:,:,:,:),allocatable::qRT
     real(dp),dimension(:,:,:,:,:),allocatable::qRB
     real(dp),dimension(:,:,:,:,:),allocatable::qLT
     real(dp),dimension(:,:,:,:,:),allocatable::qLB
     type(nbor),dimension(:,:,:,:),allocatable::nborsonloc
#endif
#ifdef RT
     real(dp),dimension(:,:,:,:),allocatable::rtuloc
     real(dp),dimension(:,:,:,:,:),allocatable::rtflux
#endif
   contains
     procedure :: init => init_hydro_kernel
     procedure :: size => size_hydro_kernel
  end type hydro_kernel_t

  type hydro_workspace_t
     type(hydro_kernel_t)::kernel_1,kernel_2,kernel_4,kernel_8,kernel_16,kernel_32
  end type hydro_workspace_t

contains
  
  subroutine init_hydro_kernel(h,nn)
    use amr_parameters, only: ndim
    use hydro_parameters, only: nvar
#ifdef RT
    use rt_parameters, only: nrtvar
#endif

    integer::nn
    class(hydro_kernel_t)::h

    h%iu1=-1; h%iu2=nn+2
    h%ju1=(1-ndim/2)-1*(ndim/2); h%ju2=(1-ndim/2)+(nn+2)*(ndim/2)
    h%ku1=(1-ndim/3)-1*(ndim/3); h%ku2=(1-ndim/3)+(nn+2)*(ndim/3)
    h%if1=1; h%if2=nn+1
    h%jf1=1; h%jf2=(1-ndim/2)+(nn+1)*(ndim/2)
    h%kf1=1; h%kf2=(1-ndim/3)+(nn+1)*(ndim/3)
    h%io1=0; h%io2=nn/2+1
    h%jo1=(1-ndim/2); h%jo2=(1-ndim/2)+(nn/2+1)*(ndim/2)
    h%ko1=(1-ndim/3); h%ko2=(1-ndim/3)+(nn/2+1)*(ndim/3)
    
    allocate(h%uloc (h%iu1:h%iu2,h%ju1:h%ju2,h%ku1:h%ku2,1:nvar))
    allocate(h%gloc (h%iu1:h%iu2,h%ju1:h%ju2,h%ku1:h%ku2,1:ndim))
    allocate(h%qloc (h%iu1:h%iu2,h%ju1:h%ju2,h%ku1:h%ku2,1:nprim))
    allocate(h%cloc (h%iu1:h%iu2,h%ju1:h%ju2,h%ku1:h%ku2))
    allocate(h%okloc(h%iu1:h%iu2,h%ju1:h%ju2,h%ku1:h%ku2))
    allocate(h%dq   (h%iu1:h%iu2,h%ju1:h%ju2,h%ku1:h%ku2,1:nprim,1:ndim))
    allocate(h%qm   (h%iu1:h%iu2,h%ju1:h%ju2,h%ku1:h%ku2,1:nprim,1:ndim))
    allocate(h%qp   (h%iu1:h%iu2,h%ju1:h%ju2,h%ku1:h%ku2,1:nprim,1:ndim))
    allocate(h%fx   (h%iu1:h%iu2,h%ju1:h%ju2,h%ku1:h%ku2,1:nprim))
    allocate(h%tx   (h%iu1:h%iu2,h%ju1:h%ju2,h%ku1:h%ku2,1:2    ))

    allocate(h%flux(h%if1:h%if2,h%jf1:h%jf2,h%kf1:h%kf2,1:nprim,1:ndim))
    allocate(h%tmp (h%if1:h%if2,h%jf1:h%jf2,h%kf1:h%kf2,1:2    ,1:ndim))
    allocate(h%divu(h%if1:h%if2,h%jf1:h%jf2,h%kf1:h%kf2))

    allocate(h%childloc(h%io1:h%io2,h%jo1:h%jo2,h%ko1:h%ko2))
    allocate(h%gridloc (h%io1:h%io2,h%jo1:h%jo2,h%ko1:h%ko2))
    allocate(h%cellloc (h%io1:h%io2,h%jo1:h%jo2,h%ko1:h%ko2))
    allocate(h%nborloc (h%io1:h%io2,h%jo1:h%jo2,h%ko1:h%ko2,1:twondim))
#ifdef MHD
    allocate(h%bloc(h%iu1:h%iu2,h%ju1:h%ju2,h%ku1:h%ku2,1:6))
    allocate(h%emfx(h%if1:h%if2,h%jf1:h%jf2,h%kf1:h%kf2))
    allocate(h%emfy(h%if1:h%if2,h%jf1:h%jf2,h%kf1:h%kf2))
    allocate(h%emfz(h%if1:h%if2,h%jf1:h%jf2,h%kf1:h%kf2))
    allocate(h%Ex(h%iu1:h%iu2,h%ju1:h%ju2,h%ku1:h%ku2))
    allocate(h%Ey(h%iu1:h%iu2,h%ju1:h%ju2,h%ku1:h%ku2))
    allocate(h%Ez(h%iu1:h%iu2,h%ju1:h%ju2,h%ku1:h%ku2))
    allocate(h%bf(h%iu1:h%iu2+1,h%ju1:h%ju2+1,h%ku1:h%ku2+1,1:3))
    allocate(h%dbf(h%iu1:h%iu2+1,h%ju1:h%ju2+1,h%ku1:h%ku2+1,1:3,1:2))
    allocate(h%qRT(h%iu1:h%iu2,h%ju1:h%ju2,h%ku1:h%ku2,1:nprim,1:3))
    allocate(h%qRB(h%iu1:h%iu2,h%ju1:h%ju2,h%ku1:h%ku2,1:nprim,1:3))
    allocate(h%qLT(h%iu1:h%iu2,h%ju1:h%ju2,h%ku1:h%ku2,1:nprim,1:3))
    allocate(h%qLB(h%iu1:h%iu2,h%ju1:h%ju2,h%ku1:h%ku2,1:nprim,1:3))
    allocate(h%nborsonloc(h%io1:h%io2,h%jo1:h%jo2,h%ko1:h%ko2,1:twondim))
#endif
#ifdef RT
    allocate(h%rtuloc (h%iu1:h%iu2,h%ju1:h%ju2,h%ku1:h%ku2,1:nrtvar))    ! Joki this needs to be eventually moved to a separate rt kernel
    allocate(h%rtflux(h%if1:h%if2,h%jf1:h%jf2,h%kf1:h%kf2,1:nrtvar,1:ndim))
#endif
  end subroutine init_hydro_kernel

  function size_hydro_kernel(h)
    use amr_parameters, only: ndim
    use hydro_parameters, only: nvar
    integer::size_hydro_kernel
    class(hydro_kernel_t)::h

    integer::nint

    nint=0
    nint=nint+size(transfer(h%uloc,(/1/)))
    nint=nint+size(transfer(h%gloc,(/1/)))
    nint=nint+size(transfer(h%qloc,(/1/)))
    nint=nint+size(transfer(h%cloc,(/1/)))
    nint=nint+size(transfer(h%okloc,(/1/)))
    nint=nint+size(transfer(h%dq   ,(/1/)))
    nint=nint+size(transfer(h%qm   ,(/1/)))
    nint=nint+size(transfer(h%qp   ,(/1/)))
    nint=nint+size(transfer(h%fx   ,(/1/)))
    nint=nint+size(transfer(h%tx   ,(/1/)))

    nint=nint+size(transfer(h%flux ,(/1/)))
    nint=nint+size(transfer(h%tmp  ,(/1/)))
    nint=nint+size(transfer(h%divu ,(/1/)))

    nint=nint+size(transfer(h%childloc,(/1/)))
    nint=nint+size(transfer(h%gridloc ,(/1/)))
    nint=nint+size(transfer(h%cellloc ,(/1/)))
    nint=nint+size(transfer(h%nborloc ,(/1/)))
#ifdef MHD
    nint=nint+size(transfer(h%bloc,(/1/)))
    nint=nint+size(transfer(h%emfx,(/1/)))
    nint=nint+size(transfer(h%emfy,(/1/)))
    nint=nint+size(transfer(h%emfz,(/1/)))
    nint=nint+size(transfer(h%Ex  ,(/1/)))
    nint=nint+size(transfer(h%Ey  ,(/1/)))
    nint=nint+size(transfer(h%Ez  ,(/1/)))
    nint=nint+size(transfer(h%bf  ,(/1/)))
    nint=nint+size(transfer(h%dbf ,(/1/)))
    nint=nint+size(transfer(h%qRT ,(/1/)))
    nint=nint+size(transfer(h%qRB ,(/1/)))
    nint=nint+size(transfer(h%qLT ,(/1/)))
    nint=nint+size(transfer(h%qLB ,(/1/)))
    nint=nint+size(transfer(h%nborsonloc,(/1/)))
#endif
#ifdef RT
    nint=nint+size(transfer(h%rtuloc,(/1/)))
    nint=nint+size(transfer(h%rtflux ,(/1/)))
#endif
    size_hydro_kernel = nint

  end function size_hydro_kernel

end module hydro_commons
