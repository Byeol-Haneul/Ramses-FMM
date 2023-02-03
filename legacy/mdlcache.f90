module hash
  use amr_parameters, only: ndim
  USE, INTRINSIC :: ISO_C_BINDING, ONLY: c_int64_t, c_int32_t, c_int16_t, C_FUNPTR, C_PTR
  implicit none

  integer, dimension(0:NDIM), parameter :: constants = (/5, -1640531527, 97, 1003313/)

  type hash_table
    type(c_ptr)::mdl_cache_table
  end type hash_table

  interface
    function hash_create(nElements) BIND(C,NAME="HashCreate")
      USE, INTRINSIC :: ISO_C_BINDING, ONLY : C_INT, C_PTR
      type(c_ptr)::hash_create
      integer(c_int), VALUE :: nElements
    end function hash_create

    subroutine hash_insert(table,hash,key,data) BIND(C,NAME="HashInsert")
      USE, INTRINSIC :: ISO_C_BINDING, ONLY : C_INT64_T, C_INT32_T, C_PTR
      type(c_ptr),value                           :: table
      integer(kind=4),value                       :: hash
      integer(kind=8), dimension(0:NDIM), intent(in) :: key
      integer(c_int64_t),intent(in),value         :: data
    end subroutine hash_insert

    subroutine hash_insertp(table,hash,key,data) BIND(C,NAME="HashInsert")
      USE, INTRINSIC :: ISO_C_BINDING, ONLY : C_INT64_T, C_INT32_T, C_PTR
      type(c_ptr),value                           :: table
      integer(kind=4),value                       :: hash
      integer(kind=8), dimension(0:NDIM), intent(in) :: key
      TYPE(*),intent(in),target                   :: data
    end subroutine hash_insertp

    pure function hash_lookup(table,hash,key) BIND(C,NAME="HashLookup")
      USE, INTRINSIC :: ISO_C_BINDING, ONLY : C_INT64_T, C_INT32_T, C_PTR
      integer(c_int64_t)                          :: hash_lookup
      type(c_ptr),value                           :: table
      integer(kind=4),value                       :: hash
      integer(kind=8), dimension(0:NDIM), intent(in) :: key
    end function hash_lookup

    pure function hash_lookupp(table,hash,key) BIND(C,NAME="HashLookup")
      USE, INTRINSIC :: ISO_C_BINDING, ONLY : C_INT64_T, C_INT32_T, C_PTR
      type(c_ptr)                                 :: hash_lookupp
      type(c_ptr),value                           :: table
      integer(kind=4),value                       :: hash
      integer(kind=8), dimension(0:NDIM), intent(in) :: key
    end function hash_lookupp

    subroutine hash_remove(table,hash,key) BIND(C,NAME="HashRemove")
      USE, INTRINSIC :: ISO_C_BINDING, ONLY : C_INT64_T, C_INT32_T, C_PTR
      type(c_ptr),value                           :: table
      integer(kind=4),value                       :: hash
      integer(kind=8), dimension(0:NDIM), intent(in) :: key
    end subroutine hash_remove

    subroutine hash_clear(table) BIND(C,NAME="HashClear")
      USE, INTRINSIC :: ISO_C_BINDING, ONLY : C_INT64_T, C_INT32_T, C_PTR
      type(c_ptr),value                           :: table
    end subroutine hash_clear

    subroutine hash_statistics(table) BIND(C,NAME="HashStatistics")
      USE, INTRINSIC :: ISO_C_BINDING, ONLY : C_INT64_T, C_INT32_T, C_PTR
      type(c_ptr),value                           :: table
    end subroutine hash_statistics

     subroutine pack_function_type(grid,msg_size,msg_array) BIND(C)
       !use oct_commons, only: oct
       type(*),                     intent(in) :: grid
       integer,                       intent(in) :: msg_size
       integer,dimension(1:msg_size), intent(out):: msg_array
     end subroutine pack_function_type

     subroutine unpack_function_type(grid,msg_size,msg_array) BIND(C)
       !use oct_commons, only: oct
       type(*)                                   :: grid
       integer,                       intent(in) :: msg_size
       integer,dimension(1:msg_size), intent(in) :: msg_array
     end subroutine unpack_function_type

    subroutine ramses_cache_open(mdl,cid,hash,hash_func,nDataSize,modify,ctx,get_thread,get_tile,&
                          pack_size,pack,unpack,init,flush_size,flush,combine,create)&
                          BIND(C,NAME="ramses_cache_open")
      USE, INTRINSIC :: ISO_C_BINDING, ONLY : C_INT, C_BOOL, C_PTR, C_FUNPTR
      type(c_ptr),value::mdl
      integer(c_int), VALUE             :: cid,nDataSize,pack_size,flush_size
      logical(c_bool),value             :: modify
      type(c_ptr),value                 :: hash
      type(c_funptr), intent(in), VALUE :: hash_func
      type(*),target                    :: ctx
      type(c_funptr), intent(in), VALUE :: get_thread,get_tile,pack,unpack,init,flush,combine,create
    end subroutine ramses_cache_open

  end interface

  contains
  ! =============================================================================
  pure function hash_func(key)
    integer(kind=8), dimension(0:ndim), intent(in) :: key
    integer(kind=4)                                :: hash_func

    hash_func = dot_product(key(0:ndim), constants(0:ndim))
  end function hash_func
  ! =============================================================================

  ! =============================================================================
  subroutine init_empty_hash(htable, req_size, hash_type)
    USE, INTRINSIC :: ISO_C_BINDING, ONLY : c_associated
    implicit none
    type(hash_table), intent(inout) :: htable
    integer         , intent(in)    :: req_size
    character(6)    , intent(in)    :: hash_type

    if (hash_type .ne. 'simple') then
       print*, 'only simple hash is currently supported'
       stop
    end if

    if (c_associated(htable%mdl_cache_table)) then
       print*, 'the hash table seems to be already setup'
       stop
    end if

    htable%mdl_cache_table = hash_create(req_size)

  end subroutine init_empty_hash
  ! =============================================================================

  ! =============================================================================
  subroutine reset_entire_hash(htable, resize)
    implicit none
    type(hash_table), intent(inout) :: htable
    logical, intent(in)             :: resize
    call hash_clear(htable%mdl_cache_table)
  end subroutine reset_entire_hash
  ! =============================================================================

  ! =============================================================================
  pure function hash_get(htable, key)
    implicit none
    type(hash_table),                     intent(in) :: htable
    integer(kind=8) , dimension(0:ndim), intent(in)  :: key
    integer                                          :: hash_get

    hash_get = hash_lookup(htable%mdl_cache_table,hash_func(key),key)

  end function hash_get
  ! =============================================================================

  ! =============================================================================
  pure function hash_getp(htable, key)
    USE, INTRINSIC :: ISO_C_BINDING, ONLY : C_PTR
    implicit none
    type(hash_table),                     intent(in) :: htable
    integer(kind=8) , dimension(0:ndim), intent(in)  :: key
    type(c_ptr)                                      :: hash_getp

    hash_getp = hash_lookupp(htable%mdl_cache_table,hash_func(key),key)

  end function hash_getp
  ! =============================================================================

  ! =============================================================================
  subroutine hash_set(htable, key, val)
    implicit none
    type(hash_table),                     intent(in) :: htable
    integer(kind=8) , dimension(0:ndim), intent(in) :: key
    integer,                              intent(in)    :: val

    call hash_insert(htable%mdl_cache_table,hash_func(key),key,int(val,8))

  end subroutine hash_set
  ! =============================================================================

  ! =============================================================================
  subroutine hash_setp(htable, key, val)
    implicit none
    type(hash_table),                     intent(in) :: htable
    integer(kind=8) , dimension(0:ndim), intent(in)  :: key
    TYPE(*),intent(in),target                        :: val

    call hash_insertp(htable%mdl_cache_table,hash_func(key),key,val)

  end subroutine hash_setp
  ! =============================================================================

  ! =============================================================================
  subroutine hash_free(htable, key)
    implicit none
    type(hash_table),                     intent(inout) :: htable
    integer(kind=8) , dimension(0:ndim), intent(in)    :: key

    call hash_remove(htable%mdl_cache_table,hash_func(key),key)

  end subroutine hash_free
  ! =============================================================================

  ! =============================================================================
  subroutine hash_stats(htable)
    implicit none
    type(hash_table),                     intent(inout) :: htable

    call hash_statistics(htable%mdl_cache_table)

  end subroutine hash_stats
  ! =============================================================================

end module hash