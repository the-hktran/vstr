        SUBROUTINE calc_ch4_pot(x, v)
        implicit real*8 (a-h,o-z)

        parameter (natom = 5)
        parameter (ndim = 3*natom)
        double precision, dimension(natom,3), intent(in) :: x
        double precision, intent(out) :: v
        double precision :: bohr_to_ang
        double precision :: rad_to_deg
        double precision :: wn_to_au
        dimension rij(10)
        
        rad_to_deg = 57.2958
        bohr_to_ang = 0.529177
        wn_to_au = 4.55633D-6

        r1 = 0.
        r2 = 0.
        r3 = 0.
        r4 = 0.
        c12 = 0.
        c13 = 0.
        c14 = 0.
        c23 = 0.
        c24 = 0.
        c34 = 0.
        
        do j = 1,3
                x1 = x(1, j) - x(5, j)
                r1 = r1 + x1**2
                x2 = x(2, j) - x(5, j)
                r2 = r2 + x2**2
                x3 = x(3, j) - x(5, j)
                r3 = r3 + x3**2
                x4 = x(4, j) - x(5, j)
                r4 = r4 + x4**2
                c12 = c12 + x1 * x2
                c13 = c13 + x1 * x3
                c14 = c14 + x1 * x4
                c23 = c23 + x2 * x3
                c24 = c24 + x2 * x4
                c34 = c34 + x3 * x4
        enddo
        rij(1) = sqrt(r1) * bohr_to_ang
        rij(2) = sqrt(r2) * bohr_to_ang
        rij(3) = sqrt(r3) * bohr_to_ang
        rij(4) = sqrt(r4) * bohr_to_ang
        rij(5) = acos(c12*bohr_to_ang**2/rij(1)/rij(2)) * rad_to_deg
        rij(6) = acos(c13*bohr_to_ang**2/rij(1)/rij(3)) * rad_to_deg
        rij(7) = acos(c14*bohr_to_ang**2/rij(1)/rij(4)) * rad_to_deg
        rij(8) = acos(c23*bohr_to_ang**2/rij(2)/rij(3)) * rad_to_deg
        rij(9) = acos(c24*bohr_to_ang**2/rij(2)/rij(4)) * rad_to_deg
        rij(10) = acos(c34*bohr_to_ang**2/rij(3)/rij(4)) * rad_to_deg

        call ch4_pot(rij,v)
        v = v * wn_to_au
        end
