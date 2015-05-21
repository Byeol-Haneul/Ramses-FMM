SUBROUTINE quick_sort(list, order, n)
  
  ! Quick sort routine from:
  ! Brainerd, W.S., Goldberg, C.H. & Adams, J.C. (1990) "Programmer's Guide to
  ! Fortran 90", McGraw-Hill  ISBN 0-07-000248-7, pages 149-150.
  ! Modified by Alan Miller to include an associated integer array which gives
  ! the positions of the elements in the original order.
  
  use amr_parameters, ONLY: qdp
  IMPLICIT NONE
  INTEGER :: n
  REAL(qdp), DIMENSION (1:n), INTENT(INOUT)  :: list
  INTEGER, DIMENSION (1:n), INTENT(OUT)  :: order
  
  ! Local variable
  INTEGER :: i
  
  DO i = 1, n
     order(i) = i
  END DO
  
  CALL quick_sort_1(1, n)
  
CONTAINS
  
  RECURSIVE SUBROUTINE quick_sort_1(left_end, right_end)
    
    use amr_parameters, ONLY: qdp
    INTEGER, INTENT(IN) :: left_end, right_end
    
    !     Local variables
    INTEGER             :: i, j, itemp
    REAL(qdp)              :: reference, temp
    INTEGER, PARAMETER  :: max_simple_sort_size = 6
    
    IF (right_end < left_end + max_simple_sort_size) THEN
       ! Use interchange sort for small lists
       CALL interchange_sort(left_end, right_end)
       
    ELSE
       ! Use partition ("quick") sort
       reference = list((left_end + right_end)/2)
       i = left_end - 1; j = right_end + 1
       
       DO
          ! Scan list from left end until element >= reference is found
          DO
             i = i + 1
             IF (list(i) >= reference) EXIT
          END DO
          ! Scan list from right end until element <= reference is found
          DO
             j = j - 1
             IF (list(j) <= reference) EXIT
          END DO
          
          
          IF (i < j) THEN
             ! Swap two out-of-order elements
             temp = list(i); list(i) = list(j); list(j) = temp
             itemp = order(i); order(i) = order(j); order(j) = itemp
          ELSE IF (i == j) THEN
             i = i + 1
             EXIT
          ELSE
             EXIT
          END IF
       END DO
       
       IF (left_end < j) CALL quick_sort_1(left_end, j)
       IF (i < right_end) CALL quick_sort_1(i, right_end)
    END IF
    
  END SUBROUTINE quick_sort_1
  
  
  SUBROUTINE interchange_sort(left_end, right_end)
    
    use amr_parameters, ONLY: qdp
    INTEGER, INTENT(IN) :: left_end, right_end
    
    !     Local variables
    INTEGER             :: i, j, itemp
    REAL(qdp)           :: temp
    
    DO i = left_end, right_end - 1
       DO j = i+1, right_end
          IF (list(i) > list(j)) THEN
             temp = list(i); list(i) = list(j); list(j) = temp
             itemp = order(i); order(i) = order(j); order(j) = itemp
          END IF
       END DO
    END DO
    
  END SUBROUTINE interchange_sort
  
END SUBROUTINE quick_sort
!########################################################################
!########################################################################
!########################################################################
!########################################################################
!########################################################################
SUBROUTINE quick_sort_dp(list, order, n)
  use amr_parameters, ONLY: dp
  IMPLICIT NONE
  ! Quick sort routine from:
  ! Brainerd, W.S., Goldberg, C.H. & Adams, J.C. (1990) "Programmer's Guide to
  ! Fortran 90", McGraw-Hill  ISBN 0-07-000248-7, pages 149-150.
  ! Modified by Alan Miller to include an associated integer array which gives
  ! the positions of the elements in the original order.


  INTEGER :: n
  REAL(dp), DIMENSION (1:n), INTENT(INOUT)  :: list
  INTEGER, DIMENSION (1:n), INTENT(OUT)  :: order

  ! Local variable
  INTEGER :: i

  DO i = 1, n
     order(i) = i
  END DO

  CALL quick_sort_1_dp(1, n)

CONTAINS

  RECURSIVE SUBROUTINE quick_sort_1_dp(left_end, right_end)

    INTEGER, INTENT(IN) :: left_end, right_end

    !     Local variables
    INTEGER             :: i, j, itemp
    REAL(kind=8)              :: reference, temp
    INTEGER, PARAMETER  :: max_simple_sort_size = 6

    IF (right_end < left_end + max_simple_sort_size) THEN
       ! Use interchange sort for small lists
       CALL interchange_sort_dp(left_end, right_end)

    ELSE
       ! Use partition ("quick") sort
       reference = list((left_end + right_end)/2)
       i = left_end - 1; j = right_end + 1

      DO
          ! Scan list from left end until element >= reference is found
          DO
             i = i + 1
             IF (list(i) >= reference) EXIT
          END DO
          ! Scan list from right end until element <= reference is found
          DO
             j = j - 1
             IF (list(j) <= reference) EXIT
          END DO


          IF (i < j) THEN
             ! Swap two out-of-order elements
             temp = list(i); list(i) = list(j); list(j) = temp
             itemp = order(i); order(i) = order(j); order(j) = itemp
          ELSE IF (i == j) THEN
             i = i + 1
             EXIT
          ELSE
             EXIT
          END IF
       END DO
       IF (left_end < j) CALL quick_sort_1_dp(left_end, j)
       IF (i < right_end) CALL quick_sort_1_dp(i, right_end)
    END IF

  END SUBROUTINE quick_sort_1_dp


  SUBROUTINE interchange_sort_dp(left_end, right_end)

    INTEGER, INTENT(IN) :: left_end, right_end

    !     Local variables                                                                                                                                                                           
    INTEGER             :: i, j, itemp
    REAL(kind=8)           :: temp

    DO i = left_end, right_end - 1
       DO j = i+1, right_end
          IF (list(i) > list(j)) THEN
             temp = list(i); list(i) = list(j); list(j) = temp
             itemp = order(i); order(i) = order(j); order(j) = itemp
          END IF
       END DO
    END DO

  END SUBROUTINE interchange_sort_dp

END SUBROUTINE quick_sort_dp

!########################################################################
!########################################################################
!########################################################################
!########################################################################
!########################################################################
SUBROUTINE quick_sort_keys(order, n)
 
  ! Quick sort routine from:
  ! Brainerd, W.S., Goldberg, C.H. & Adams, J.C. (1990) "Programmer's Guide to
  ! Fortran 90", McGraw-Hill  ISBN 0-07-000248-7, pages 149-150.
  ! Modified by Alan Miller to include an associated integer array which gives
  ! the positions of the elements in the original order.
  
  USE pm_parameters, only:npartmax
  USE pm_commons, only:big_hkey
  IMPLICIT NONE
  INTEGER :: n
  INTEGER, DIMENSION (1:n), INTENT(OUT)  :: order
  
  ! Local variable
  INTEGER :: i
  
  DO i = 1, n
     order(i) = i
  END DO
  
  CALL quick_sort_1_keys(1, n)
  
CONTAINS
  
  RECURSIVE SUBROUTINE quick_sort_1_keys(left_end, right_end)
    
    INTEGER, INTENT(IN) :: left_end, right_end
    !     Local variables
    INTEGER             :: i, j, itemp
    INTEGER(kind=8),DIMENSION(0:2):: reference, temp
    INTEGER, PARAMETER  :: max_simple_sort_size = 6
    
    IF (right_end < left_end + max_simple_sort_size) THEN
       ! Use interchange sort for small lists
       CALL interchange_sort_keys(left_end, right_end)       
    ELSE
       ! Use partition ("quick") sort
       reference(0:2) = big_hkey((left_end + right_end)/2,0:2)
       i = left_end - 1; j = right_end + 1
       
       DO
          ! Scan list from left end until element >= reference is found
          DO
             i = i + 1
             IF (ge_3keys(big_hkey(i,0:2),reference(0:2))) EXIT
          END DO
          ! Scan list from right end until element <= reference is found
          DO
             j = j - 1
             IF (ge_3keys(reference(0:2),big_hkey(j,0:2))) EXIT
          END DO
          
          
          IF (i < j) THEN
             ! Swap two out-of-order elements
             temp = big_hkey(i,0:2); big_hkey(i,0:2) = big_hkey(j,0:2); big_hkey(j,0:2) = temp
             itemp = order(i); order(i) = order(j); order(j) = itemp
          ELSE IF (i == j) THEN
             i = i + 1
             EXIT
          ELSE
             EXIT
          END IF
       END DO
       
       IF (left_end < j) CALL quick_sort_1_keys(left_end, j)
       IF (i < right_end) CALL quick_sort_1_keys(i, right_end)
    END IF
    
  END SUBROUTINE quick_sort_1_keys
  
  
  SUBROUTINE interchange_sort_keys(left_end, right_end)
    
    use amr_parameters, ONLY: qdp
    INTEGER, INTENT(IN) :: left_end, right_end
    !     Local variables
    INTEGER             :: i, j, itemp
    INTEGER(kind=8),DIMENSION(0:2):: temp
    
    DO i = left_end, right_end - 1
       DO j = i+1, right_end
          IF (gt_3keys(big_hkey(i,0:2),big_hkey(j,0:2))) THEN
             temp(0:2) = big_hkey(i,0:2); big_hkey(i,0:2) = big_hkey(j,0:2); big_hkey(j,0:2) = temp(0:2)
             itemp = order(i); order(i) = order(j); order(j) = itemp
          END IF
       END DO
    END DO
    
  END SUBROUTINE interchange_sort_keys

  function ge_2keys(key_a,key_b)
    integer(kind=8),dimension(0:2)::key_a,key_b
    logical::ge_2keys
    ! Function to test wether a >= b for a and b two-integer hilbert keys 
    if(key_a(1) > key_b(1))then
       ge_2keys=.true.
    elseif(key_a(1) < key_b(1))then
       ge_2keys=.false.
    elseif(key_a(0) > key_b(0))then
       ge_2keys=.true.
    elseif(key_a(0) < key_b(0))then
       ge_2keys=.false.
    else
       ge_2keys=.true.
    end if
  end function ge_2keys

  function ge_3keys(key_a,key_b)
    integer(kind=8),dimension(0:2)::key_a,key_b
    logical::ge_3keys
    ! Function to test wether a >= b for a and b two-integer hilbert keys 
    if(key_a(2) > key_b(2))then
       ge_3keys=.true.
    elseif(key_a(2) < key_b(2))then
       ge_3keys=.false.
    elseif(key_a(1) > key_b(1))then
       ge_3keys=.true.
    elseif(key_a(1) < key_b(1))then
       ge_3keys=.false.
    elseif(key_a(0) > key_b(0))then
       ge_3keys=.true.
    elseif(key_a(0) < key_b(0))then
       ge_3keys=.false.
    else
       ge_3keys=.true.
    end if
  end function ge_3keys

  function gt_3keys(key_a,key_b)
    integer(kind=8),dimension(0:2)::key_a,key_b
    logical::gt_3keys
    ! Function to test wether a >= b for a and b two-integer hilbert keys 
    if(key_a(2) > key_b(2))then
       gt_3keys=.true.
    elseif(key_a(2) < key_b(2))then
       gt_3keys=.false.
    elseif(key_a(1) > key_b(1))then
       gt_3keys=.true.
    elseif(key_a(1) < key_b(1))then
       gt_3keys=.false.
    elseif(key_a(0) > key_b(0))then
       gt_3keys=.true.
    elseif(key_a(0) < key_b(0))then
       gt_3keys=.false.
    else
       gt_3keys=.false.
    end if
  end function gt_3keys
END SUBROUTINE quick_sort_keys



