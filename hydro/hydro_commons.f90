module hydro_commons
  use amr_parameters
  use hydro_parameters

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
     integer ,dimension(:,:,:),allocatable::childloc
     integer ,dimension(:,:,:),allocatable::gridloc
     integer ,dimension(:,:,:),allocatable::cellloc
     integer ,dimension(:,:,:,:),allocatable::nborloc
   contains
     procedure :: init => init_hydro_kernel
  end type hydro_kernel_t

  type hydro_workspace_t
     type(hydro_kernel_t)::kernel_1,kernel_2,kernel_4,kernel_8,kernel_16,kernel_32
  end type hydro_workspace_t

contains
  
  subroutine init_hydro_kernel(h,nn)
    use amr_parameters, only: ndim
    use hydro_parameters, only: nvar
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
    allocate(h%qloc (h%iu1:h%iu2,h%ju1:h%ju2,h%ku1:h%ku2,1:nvar))
    allocate(h%cloc (h%iu1:h%iu2,h%ju1:h%ju2,h%ku1:h%ku2))
    allocate(h%okloc(h%iu1:h%iu2,h%ju1:h%ju2,h%ku1:h%ku2))
    allocate(h%dq   (h%iu1:h%iu2,h%ju1:h%ju2,h%ku1:h%ku2,1:nvar,1:ndim))
    allocate(h%qm   (h%iu1:h%iu2,h%ju1:h%ju2,h%ku1:h%ku2,1:nvar,1:ndim))
    allocate(h%qp   (h%iu1:h%iu2,h%ju1:h%ju2,h%ku1:h%ku2,1:nvar,1:ndim))
    allocate(h%fx   (h%iu1:h%iu2,h%ju1:h%ju2,h%ku1:h%ku2,1:nvar))
    allocate(h%tx   (h%iu1:h%iu2,h%ju1:h%ju2,h%ku1:h%ku2,1:2   ))

    allocate(h%flux(h%if1:h%if2,h%jf1:h%jf2,h%kf1:h%kf2,1:nvar,1:ndim))
    allocate(h%tmp (h%if1:h%if2,h%jf1:h%jf2,h%kf1:h%kf2,1:2   ,1:ndim))
    allocate(h%divu(h%if1:h%if2,h%jf1:h%jf2,h%kf1:h%kf2))

    allocate(h%childloc(h%io1:h%io2,h%jo1:h%jo2,h%ko1:h%ko2))
    allocate(h%gridloc (h%io1:h%io2,h%jo1:h%jo2,h%ko1:h%ko2))
    allocate(h%cellloc (h%io1:h%io2,h%jo1:h%jo2,h%ko1:h%ko2))
    allocate(h%nborloc (h%io1:h%io2,h%jo1:h%jo2,h%ko1:h%ko2,1:twondim))
    
  end subroutine init_hydro_kernel

end module hydro_commons
