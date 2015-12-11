subroutine cic(xpart, cell_index, vol, np, cic_level, level_boundary_case)
   use amr_parameters,  only: static, dp, twotondim
   use amr_commons,     only: boxlen, nvector, ndim
   use hilbert,         only: hilbert3d
   implicit none

   integer,  intent(in)                                      :: np, cic_level, level_boundary_case
   integer(kind=4), intent(inout), dimension(1:nvector, 1:8) :: cell_index
   real(dp),        intent(inout), dimension(1:nvector, 1:8) :: vol
   real(dp), intent(in), dimension(1:nvector, 1:ndim)             :: xpart

   ! Subroutine to do the Cloud-in-Cell interpolation for nvector particle positions at level cic_level.

   ! Input variables
   ! - xpart: particle coordinates
   ! - cic_level: grid level at which interpolation takes place)
   ! - np: number of particles
   ! - level_boundary_case: 1 (default) if CIC is performed only at cic_level or 2 if
   !   interpolation should be repeated at level cic_level - 1 whenever parts of the
   !   CIC cloud fall outside the boundaries of cic_level
                            
   ! Output variables:
   ! - cell_index: indices of the cells touched by the CIC clouds
   ! - vol: Fractional volume of of a particle that falls in a given cell

   integer(kind=8), dimension(1:nvector, 0:2),    save :: cloud_hkey
   integer(kind=8), dimension(1:nvector, 1:ndim), save :: ix
   integer(kind=4), dimension(1:nvector),         save :: dummy_state
   integer(kind=4), dimension(1:nvector),         save :: cell_level
   real(dp),   dimension(1:nvector, 0:1, 1:ndim), save :: cloud_boundary
   real(dp),        dimension(1:nvector, 1:3),    save :: xpart_grid
   logical,         dimension(1:nvector),         save :: repeat_coarser
   integer,         dimension(1:ndim),            save :: ind
   real(dp),        dimension(1:ndim),            save :: delta
   integer,  save :: idim, ind_cloud, ip
   real(dp), save :: part_to_grid
   integer(kind=8), save :: grid_size

   grid_size = 2_8**cic_level
   
   if (level_boundary_case==2) repeat_coarser = .false.

   ! Convert particle coordinates (0 to boxlen)
   ! into grid-coordinates (0 to 2.**grid_level)
   part_to_grid = 2.0**cic_level / boxlen
   xpart_grid(1:np, 1:ndim) = xpart(1:np, 1:ndim) * part_to_grid

   ! Compute distances of cloud boundary from nearest "integer coordinate"
   do idim = 1, ndim       

      ! upper/right/front boundary
      do ip=1,np
         cloud_boundary(ip,1,idim) = xpart_grid(ip, idim) + 0.5D0
      end do

      ! upper/rigt/front boundary rel to nearest integer
      do ip=1,np
         cloud_boundary(ip,1,idim) = cloud_boundary(ip,1,idim) - floor(cloud_boundary(ip,1,idim), kind = 8)
      end do

      ! lower/left/back boundary rel to nearest integer
      do ip=1,np
         cloud_boundary(ip,0,idim) = 1.0D0 - cloud_boundary(ip,1,idim)
      end do
   end do

#if NDIM==1
   ! Loop cloud/cell intersections
   do ind_cloud = 0, 1
      ind(1) = ind_cloud 

      ! Compute cloud volume
      do ip=1,np
         vol(ip, ind_cloud + 1) = cloud_boundary(ip,ind(1),1) 
      end do
#endif
#if NDIM==2
   ! Loop cloud/cell intersections
   do ind_cloud = 0, 3
      ind(1) = ind_cloud/2
      ind(2) = mod(ind_cloud,2)
      
      ! Compute cloud volume
      do ip=1,np
         vol(ip, ind_cloud + 1) = cloud_boundary(ip,ind(1),1) * &
              cloud_boundary(ip,ind(2),2) 
      end do
#endif
#if NDIM==3
   ! Loop cloud/cell intersections
   do ind_cloud = 0, 7
      ind(1) = ind_cloud/4
      ind(2) = mod(ind_cloud,4)/2
      ind(3) = mod(mod(ind_cloud,4),2)
      
      ! Compute cloud volume
      do ip=1,np
         vol(ip, ind_cloud + 1) = cloud_boundary(ip,ind(1),1) * &
              cloud_boundary(ip,ind(2),2) * &
              cloud_boundary(ip,ind(3),3) 
      end do
#endif
         
      ! Compute cloud corner offset from cloud center
      delta(1:ndim) = ind(1:ndim) - 0.5D0       
      
      ! Get cell indices where the cloud corners fall into
      ! (cartesian key -> hilbert key -> cell index)
      ! TODO: WHAT IF PARTICLE SITS CLOSE TO PERIODIC BOX BOUNDARY WITH AMR  -> should be ok
      do idim = 1, ndim
         do ip = 1, np
            ix(ip,idim) = floor(xpart_grid(ip,idim) + delta(idim), kind = 8)
         end do
      end do
      do idim = 1, ndim
         do ip = 1, np
            if (ix(ip, idim) >= grid_size)then
               ix(ip, idim) = ix(ip, idim) - grid_size
            end if
            if (ix(ip, idim) < 0) then
               ix(ip, idim) = ix(ip, idim) + grid_size
            end if
         end do
      end do
      
      
      ! call hilbert3d(ix(1:np,1), ix(1:np,2), ix(1:np,3), &
      !      cloud_hkey(1:np, 2), cloud_hkey(1:np, 1), cloud_hkey(1:np, 0), &
      !      dummy_state(1:np), 0, cic_level, np)

      ! call get_cell_index_from_hilbertkey(cell_index(1:np, ind_cloud + 1), cell_level(1:np), &
      !      cloud_hkey(1:np, 2), cloud_hkey(1:np, 1), cloud_hkey(1:np, 0), np, cic_level)

      call get_cell_index_from_cartesian_hash(cell_index(1:np, ind_cloud + 1), cell_level(1:np), &
            ix(1:np, 1), ix(1:np, 2), ix(1:np, 3), cic_level, np, cic_level)  
      
      ! Exclude cloud fraction which lies in coarser level
      if (level_boundary_case == 1)then
         do ip = 1, np        
            if(cell_level(ip) < cic_level)then
               vol(ip, ind_cloud + 1) = 0.d0
            end if
         end do
      end if
      if (level_boundary_case == 2)then
         do ip = 1, np        
            if(cell_level(ip) < cic_level)then
               repeat_coarser(ip) = .true.
            end if
         end do
      end if

   end do ! end loop over cloud/cell intersections
   if (level_boundary_case == 2)then
      do ip = 1, np        
         if (repeat_coarser(ip)) then
            call cic_one(xpart(ip,1:ndim), cell_index(ip, 1:twotondim), vol(ip, 1:twotondim), cic_level - 1)
         end if
      end do
   end if
      
end subroutine cic

 subroutine cic_one(xpart, cell_index, vol, cic_level)
   use amr_parameters,  only: static, mass_cut_refine, dp, twotondim
   use amr_commons,     only: boxlen, nvector, ndim
   use hilbert,         only: hilbert3d
   implicit none
   integer,  intent(in)                           :: cic_level
   integer(kind=4), intent(inout), dimension(1:1,1:8) :: cell_index
   real(dp),        intent(inout), dimension(1:8) :: vol
   real(dp), intent(in), dimension(1:ndim)        :: xpart
   
   ! Pretty ugly copy-paste version of the cic routine, to be used
   ! for one particle only. Could be replaced with a recursive version of cic
   ! which takes a mask array as an argument. However, this might slow the
   ! original routine down. GET RID OF THIS AT SOME POINT
   
   integer(kind=8), dimension(1:1,0:2),    save :: cloud_hkey
!   integer(kind=8), dimension(1:ndim), save :: id
   integer(kind=8), dimension(1:1,1:ndim), save :: ix
   real(dp),   dimension(0:1, 1:ndim), save :: cloud_boundary
   real(dp), dimension(1:ndim), save :: xpart_grid, delta
   integer,  dimension(1:ndim), save :: ind
   integer,  save :: idim, ind_cloud
   integer, dimension(1:1), save :: cell_level
   integer,  dimension(1:1),  save ::  dummy_state
   real(dp), save :: part_to_grid
   integer(kind=8), save :: grid_size

   grid_size = 2**cic_level
   ! Convert particle coordinates (0 to boxlen)
   ! into grid-coordinates (0 to 2.**grid_level)
   part_to_grid = 2.0**cic_level / dble(boxlen)
   xpart_grid(1:ndim) = xpart(1:ndim) * part_to_grid
   
   
   ! Compute distances of cloud boundary from nearest "integer coordinate"
   do idim = 1, ndim       
      
      ! upper/right/front boundary
      cloud_boundary(1,idim) = xpart_grid(idim) + 0.5D0
      
      ! nearest integer coordinate
!      id(idim) = int(cloud_boundary(1,idim), kind=8)
      
      ! upper/rigt/front boundary rel to nearest integer
      cloud_boundary(1,idim) = cloud_boundary(1,idim) - floor(cloud_boundary(1,idim), kind=8)!id(idim)
      
      ! lower/left/back boundary rel to nearest integer
      cloud_boundary(0,idim) = 1.0D0 - cloud_boundary(1,idim)
   end do

#if NDIM==1
   ! Loop cloud/cell intersections
   do ind_cloud = 0, 1
      ind(1) = ind_cloud 

      ! Compute cloud volume
      vol(ind_cloud + 1) = cloud_boundary(ind(1),1) * &
              cloud_boundary(ind(2),2) 
#endif
#if NDIM==2
   ! Loop cloud/cell intersections
   do ind_cloud = 0, 3
      ind(1) = ind_cloud/2
      ind(2) = mod(ind_cloud,2)
      
      ! Compute cloud volume
      vol(ind_cloud + 1) = cloud_boundary(ind(1),1) * &
           cloud_boundary(ind(2),2) 
#endif
#if NDIM==3
   ! Loop cloud/cell intersections
   do ind_cloud = 0, 7
      ind(1) = ind_cloud/4
      ind(2) = mod(ind_cloud,4)/2
      ind(3) = mod(mod(ind_cloud,4),2)
      
      ! Compute cloud volume
      vol(ind_cloud + 1) = cloud_boundary(ind(1),1) * &
           cloud_boundary(ind(2),2) * &
           cloud_boundary(ind(3),3) 
#endif
         
      ! Compute cloud corner offset from cloud center
      delta(1:ndim) = ind(1:ndim) - 0.5D0       
      
      ! Get cell indices where the cloud corners fall into
      ! (cartesian key -> hilbert key -> cell index)
      do idim = 1, ndim
         ix(1,idim) = modulo(floor(xpart_grid(idim) + delta(idim), kind = 8), grid_size)
      end do

      ! call hilbert3d(ix(1,1), ix(1,2), ix(1,3), &
      !      cloud_hkey(1:1,2), cloud_hkey(1:1,1), cloud_hkey(1:1,0), &
      !      dummy_state, 0, cic_level, 1)

      ! call get_cell_index_from_hilbertkey(cell_index(1,ind_cloud + 1), cell_level(1), &
      !     cloud_hkey(1,2), cloud_hkey(1,1), cloud_hkey(1,0), 1, cic_level)

      call get_cell_index_from_cartesian_hash(cell_index(1, ind_cloud + 1), cell_level(1), &
           ix(1, 1), ix(1, 2), ix(1, 3), cic_level, 1, cic_level)  
      
   end do ! end loop over cloud/cell intersections

end subroutine cic_one
