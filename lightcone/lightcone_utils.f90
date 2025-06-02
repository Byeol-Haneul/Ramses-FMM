module lightcone_utils
    use amr_parameters, only: dp
    use amr_commons, only: run_t, global_t

    implicit none    
    contains

    function cone_to_box_coordinates(rotation_matrix, observer, x_cone) result(x_box)
        ! Converts coordinates from the cone coordinate system to the box coordinate system (code units)
        real(dp), intent(in) :: rotation_matrix(3,3), observer(3)
        real(dp), intent(in) :: x_cone(3)
        real(dp) :: x_box(3)

        x_box = matmul(rotation_matrix, x_cone)
        x_box = x_box + observer
    end function cone_to_box_coordinates

    function box_to_cone_coordinates(rotation_matrix, observer, x_box) result(x_cone)
        ! Converts coordinates from the box coordinate system to the cone coordinate system (code units)
        real(dp), intent(in) :: rotation_matrix(3,3), observer(3)
        real(dp), intent(in) :: x_box(3)
        real(dp) :: x_cone(3)

        x_cone = x_box - observer
        x_cone = matmul(rotation_matrix, x_cone)
    end function box_to_cone_coordinates

    function rotation_matrix(theta, phi) result(rot)
        ! Computes the rotation matrix for a given theta and phi (in radians)
        real(dp), intent(in) :: theta, phi
        real(dp) :: rot(3,3)
        real(dp) :: cos_theta, sin_theta, cos_phi, sin_phi

        cos_theta = cos(theta)
        sin_theta = sin(theta)
        cos_phi = cos(phi)
        sin_phi = sin(phi)

        rot(1, 1) = cos_theta * cos_phi
        rot(1, 2) = -sin_phi
        rot(1, 3) = cos_phi * sin_theta
        rot(2, 1) = sin_phi * cos_theta
        rot(2, 2) = cos_phi
        rot(2, 3) = sin_phi * sin_theta
        rot(3, 1) = -sin_theta
        rot(3, 2) = 0
        rot(3, 3) = cos_theta
    end function rotation_matrix

    function frustum_corners(x1, x2, y_angle, z_angle) result(corners)
        ! Computes the corners of a frustum delimited by the inner and outer radii x1, x2 and opening angles y_angle, z_angle (in radians)
        real(dp), intent(in) :: x1, x2, y_angle, z_angle
        real(dp) :: corners(3,8)

        ! Part of the polygon close to the observer
        corners(1, 1:4) = x1/sqrt(1.0d0 + tan(y_angle)**2 + tan(z_angle)**2)
        
        corners(2, 1) = -corners(1, 1) * tan(y_angle)
        corners(2, 2) =  corners(1, 1) * tan(y_angle) 
        corners(2, 3) = -corners(1, 1) * tan(y_angle)
        corners(2, 4) =  corners(1, 1) * tan(y_angle)
        
        corners(3, 1) = -corners(1, 1) * tan(z_angle)
        corners(3, 2) = -corners(1, 1) * tan(z_angle)
        corners(3, 3) =  corners(1, 1) * tan(z_angle)
        corners(3, 4) =  corners(1, 1) * tan(z_angle)

        ! Part of the polygon far away from the observer
        corners(1, 5:8) = x2
        
        corners(2, 5) = -x2 * tan(y_angle)
        corners(2, 6) =  x2 * tan(y_angle)
        corners(2, 7) = -x2 * tan(y_angle) 
        corners(2, 8) =  x2 * tan(y_angle)
        
        corners(3, 5) = -x2 * tan(z_angle)
        corners(3, 6) = -x2 * tan(z_angle)
        corners(3, 7) =  x2 * tan(z_angle)
        corners(3, 8) =  x2 * tan(z_angle)
    
    end function frustum_corners

    function is_in_lightcone_sector(position, r1, r2, opening_angle_y, opening_angle_z)
        ! Checks if a position is inside the lightcone sector defined by the inner and outer radii r1, r2 and angle y_angle, z_angle (in radians)
        real(dp), intent(in) :: position(3), r1, r2, opening_angle_y, opening_angle_z
        logical :: is_in_lightcone_sector
        real(dp) :: r

        r = sqrt(position(1)**2 + position(2)**2 + position(3)**2)

        is_in_lightcone_sector = (r >= r1 .and. r <= r2) .and. &
                                 (is_in_2d_sector(position(1), position(2), opening_angle_y) .and. & 
                                  is_in_2d_sector(position(1), position(3), opening_angle_z))
    end function is_in_lightcone_sector

    function is_in_2d_sector(x, y, angle)
        ! Checks if a point's angle is inside [-angle, angle]
        real(dp), intent(in) :: x, y, angle
        logical :: is_in_2d_sector
        real(dp) :: cos_xy, r

        r = sqrt(x**2 + y**2)
        if (r == 0) then
            is_in_2d_sector = .true.
            return
        end if

        cos_xy = x / r
        is_in_2d_sector = (cos_xy >= cos(angle))
    end function is_in_2d_sector

    function deg2rad(deg)
        ! Converts degrees to radians
        real(dp), intent(in) :: deg
        real(dp) :: deg2rad
        real(dp), parameter :: pi = 3.14159265358979323846_dp
        
        deg2rad = deg * pi / 180.0_dp
    end function deg2rad

    function comoving2code(g, l)
        ! Converts comoving distance (Mpc) to code units
        type(global_t), intent(in) :: g
        real(dp), intent(in) :: l
        real(dp) :: comoving2code

        comoving2code = l * g%h0 / (g%boxlen_ini * 100)
    end function comoving2code

    subroutine compute_replica_range(cone_to_box_rotation, cone_observer, angle_y, angle_z, r_inner, r_outer, &
                                     first_xreplica, last_xreplica, &
                                     first_yreplica, last_yreplica, &
                                     first_zreplica, last_zreplica)
        real(dp), intent(in) :: cone_to_box_rotation(3,3), cone_observer(3)
        real(dp), intent(in) :: angle_y, angle_z, r_inner, r_outer
        integer, intent(out) :: first_xreplica, last_xreplica, &
                                first_yreplica, last_yreplica, &
                                first_zreplica, last_zreplica
        real(dp) :: xmin, xmax, ymin, ymax, zmin, zmax
        real(dp) :: corners(3, 8)
        real(dp) :: corner(3)
        integer :: i

        ! Find the corners of the frustum in cone coordinates
        corners = frustum_corners(r_inner, r_outer, angle_y, angle_z)

        ! Transform the corners of the frustum to box coordinates
        do i = 1, 8
          corner = cone_to_box_coordinates(cone_to_box_rotation, cone_observer, corners(:,i))
          corners(:,i) = corner
        end do

        ! Find the relevant replicas range
        xmin = minval(corners(1, :))
        xmax = maxval(corners(1, :))
        ymin = minval(corners(2, :))
        ymax = maxval(corners(2, :))
        zmin = minval(corners(3, :))
        zmax = maxval(corners(3, :))

        first_xreplica = floor(xmin) ! This only works for a box of unit size
        last_xreplica = floor(xmax)
        first_yreplica = floor(ymin)
        last_yreplica = floor(ymax)
        first_zreplica = floor(zmin)
        last_zreplica = floor(zmax)

    end subroutine compute_replica_range

    !====================== The following is pasted from old ramses code - Refactor and clean up later ======================
    function comoving_distance(z, Omega0, OmegaL, OmegaR, coverH0)
        ! Computes the comoving distance at a given redshift for a given cosmology
        real(dp) :: comoving_distance, z
        real(dp) :: Omega0, OmegaL, OmegaR, coverH0
        real(dp) :: res, zz, del
        real(dp), parameter :: eps=1.0d-12
        integer :: error
        
        zz = abs(z)
        error = integrate_f(2.0d0/sqrt(1.0d0+zz), 2.0d0, eps, res, del, Omega0, OmegaL, OmegaR)
        comoving_distance = coverH0 * res
        if (z < 0) comoving_distance = -comoving_distance

    end function comoving_distance

    function f(y, Omega0, OmegaL, OmegaR)
        real(dp) :: f
        real(dp) :: y, omega0, omegaL, OmegaR

        f = 1d0/sqrt(Omega0 + OmegaR * (y/2.0d0)**2 + OmegaL * (y/2.0d0)**6)

    end function f

    function integrate_f(a, b, eps, ans, del, omega0, OmegaL, OmegaR)
        integer :: integrate_f ! This is not the value of the integral, which is in ans
        real(dp) :: a, b, eps, ans, del
        real(dp) :: Omega0, OmegaL, OmegaR  
        real(dp) :: t(0:24, 0:24)
        real(dp) :: c, d, e, s, y, x, p
        integer :: n, m, i, j, k
        
        c = 0.5d0 * (b + a)
        d = 0.5d0 * (b - a)
        n = 2
        m = 1
        e = 1.d0
        t(1, 1) = 0.d0
        t(1, 2) = 2.d0 * d * f(c, Omega0, OmegaL, OmegaR)
        t(2, 1) = 0.75d0 * t(1, 2)
        
        do while (n < 24) 
            n = n + 1
            m = m * 2
            e = e * .5d0
            s = 0.d0
            do j = 2, m, 2 
                y = dble(j-1) * e
                x = 0.5d0 * y * (3.d0 - y**2)
                s = s + (1.d0 - y**2) * (f(c-d*x, omega0, omegaL, OmegaR) &
                    &                  + f(c+d*x, omega0, omegaL, OmegaR))
            enddo
            t(n,1) = 1.5d0*s*d*e + 0.5d0*t(n-1,1)
            
            p = 1.d0
            do k = 1, n-1 
                p = p * 4.d0
                i = n + 1 - k
                t(i-1, k+1) = t(i, k) + (t(i, k) - t(i-1, k))/(p-1.d0)
            enddo
            
            ans = t(1,n)
            del = abs(t(1,n)-t(2,n-1))
            if (n >= 9) then
                if (abs(del) <= eps * abs(ans)) then
                    integrate_f = 0
                    return
                endif
            endif
        enddo
        write(*,*) 'Integration did not converge: too many iterations'
        integrate_f=1
        return

    end function integrate_f
    !=======================================================================
end module lightcone_utils

