module abundance

    use cmbtypes
    use constants
    use precision

    implicit none

    contains

    function Yp_LnLike(CMB, Yp_mean, Yp_stdev)
        Type(CMBParams)     :: CMB
        real                :: Yp_mean, Yp_stdev
        real                :: Yp_LnLike

        Yp_LnLike = (CMB%YHe - Yp_mean)**2.00e0/(2.00e0*Yp_stdev**2.00e0)

        return
    end function Yp_LnLike

    function DoH_LnLike(CMB, DoH_mean, DoH_stdev, bbn_type)
        Type(CMBParams)     :: CMB
        real                :: DoH_mean, DoH_stdev
        real                :: DoH_LnLike
        integer             :: bbn_type

        real(dl)            :: Neff_CMB, Neff_BBN
        real(dl)            :: ombh2_CMB, ombh2_BBN
        real(dl)            :: Yp, Yp_BBN, yD
        
        Neff_CMB = CMB%nnu
        ombh2_CMB = CMB%ombh2
        Yp = dble(CMB%YHe)

        call solve_dilution(Yp, Neff_CMB, ombh2_CMB, Neff_BBN, ombh2_BBN, bbn_type)
        
        if (Neff_BBN .lt. 0.00 .and. ombh2_BBN .lt. 0.00) then
            DoH_LnLike = logZero
            return
        endif

        call deuterium_helium(ombh2_BBN, Neff_BBN, Yp_BBN, yD, 2)
        if (abs(Yp_BBN-Yp) .gt. 1e-6) then
            write(*,*) 'WARNING: Yp_BBN = ', Yp_BBN, ' Yp = ', Yp
        endif            

        DoH_LnLike = (yD - DoH_mean)**2.00e0/(2.00e0*DoH_stdev**2.00e0)
        
!        print*, 'Yp = ', Yp
!        print*, 'Neff_CMB = ', Neff_CMB
!        print*, 'ombh2_CMB = ', ombh2_CMB 
!        print*, 'Neff_BBN = ', Neff_BBN
!        print*, 'ombh2_BBN = ', ombh2_BBN
!        print*, 'log[D/H] = ', yD
!        print*, 'LogLike[D/H] = ', DoH_LnLike
        !pause

        return
    end function DoH_LnLike

    subroutine deuterium_helium(ombh2, Neff, Yp, yD, method)
        
        real(dl),   intent(in)      :: ombh2
        real(dl),   intent(in)      :: Neff
        real(dl),   intent(out)     :: Yp, yD

        integer,    intent(in)      :: method

        integer         :: m, n

        real(dl)        :: eta10, S
        real(dl)        :: anm_Yp(0:5,0:3)
        real(dl)        :: anm_yD(0:4,0:3)
        
        eta10 = 273.90d0*ombh2
        S = sqrt(1.00d0 + 7.00d0*(Neff-3.00d0)/43.00d0)

        anm_Yp(0:5,0) = (/ 0.24307d0, -14.242d0,  1418.4d0, -65863.0d0,  1.4856d6,  -1.3142d7/)
        anm_Yp(0:5,1) = (/-3.6433d-2,  14.337d0, -1375.0d0,  64741.0d0, -1.4966d6,   1.3601d7/)
        anm_Yp(0:5,2) = (/ 1.6132d-2, -4.5189d0,  444.13d0, -21353.0d0,  502610.d0, -4.6405d6/)
        anm_Yp(0:5,3) = (/-1.6279d-3, 0.43362d0, -42.850d0,  2069.40d0, -48890.d0,  452740.d0/)

        anm_yD(0:4,0) = (/ 14.892d0,  -1551.6d0,  70488.d0, -1.5390d6,   1.3164d7/)
        anm_yD(0:4,1) = (/ 6.1888d0,  -916.16d0,  56639.d0, -1.6046d6,   1.7152d7/)
        anm_yD(0:4,2) = (/-0.60319d0,  118.51d0, -8556.3d0,  267580.d0, -3.0624d6/)
        anm_yD(0:4,3) = (/ 4.5346d-2, -8.7506d0,  624.51d0, -19402.d0,  221200.d0/)

        if (method .eq. 1) then
            Yp = 0.2485d0 + 0.0016d0*(eta10-6.00d0+100.0d0*(S-1.00d0))
            yD = 2.64d0*(6.00d0/(eta10-6.00d0*(S-1.00d0)))**1.60d0
        endif

        if (method .eq. 2) then
            Yp = 0.00d0
            do n=0, 5
                do m=0, 3
                    Yp = Yp + anm_Yp(n,m) * ombh2**n * Neff**m
                enddo
            enddo

            yD = 0.00d0
            do n=0, 4
                do m=0, 3
                    yD = yD + anm_yD(n,m) * ombh2**n * Neff**m
                enddo
            enddo
        endif

        yD = dlog10(yD*1.00d-5)

        return
    end subroutine deuterium_helium


    subroutine solve_dilution(Yp, Neff_CMB, ombh2_CMB, Neff_BBN, ombh2_BBN, bbn_type)
        real(dl),   intent(in)      :: Yp
        real(dl),   intent(in)      :: Neff_CMB, ombh2_CMB
        real(dl),   intent(out)     :: Neff_BBN, ombh2_BBN
        integer,    intent(in)      :: bbn_type

        real(dl)        :: try_1, try_2, lasttry
        real(dl)        :: Yp_1, Yp_2, yD, Yp_try

        try_1 = 0.00d0
        if (bbn_type .eq. 1) then
            ombh2_BBN = (try_1/Neff_CMB)**0.75d0 * ombh2_CMB
        else if (bbn_type .eq. 2) then
            ombh2_BBN = ombh2_CMB
        else if (bbn_type .eq. 3) then
            ombh2_BBN = ombh2_CMB
            Neff_BBN = Neff_CMB
            return
        endif

        call deuterium_helium(ombh2_BBN, try_1, Yp_1, yD, 2)

        try_2 = 20.00d0
        if (bbn_type .eq. 1) then
            ombh2_BBN = (try_2/Neff_CMB)**0.75d0 * ombh2_CMB
        else if (bbn_type .eq. 2) then
            ombh2_BBN = ombh2_CMB
        endif

        call deuterium_helium(ombh2_BBN, try_2, Yp_2, yD, 2)

        if (Yp .lt. Yp_1 .or. Yp .gt. Yp_2) then
!            print, "Yp_low = ", Yp_1, "Yp_high = ", Yp_2
!            print, "Yp = ", Yp
!            print, "Yp value is not in the solvable range."
        
            ombh2_BBN = -1.0d0
            Neff_BBN = -1.0d0
            return
        endif

        lasttry = -1.00d0
        do
            Neff_BBN = (try_1 + try_2)/2.0d0
            if (bbn_type .eq. 1) then
                ombh2_BBN = (Neff_BBN/Neff_CMB)**0.75d0 * ombh2_CMB
            else if (bbn_type .eq. 2) then
                ombh2_BBN = ombh2_CMB
            endif 

            call deuterium_helium(ombh2_BBN, Neff_BBN, Yp_try, yD, 2)

            if (Yp_try .lt. Yp) then
                try_1 = (try_1 + try_2)/2.0d0
            else
                try_2 = (try_1 + try_2)/2.0d0
            endif

            if (abs(Yp_try - lasttry) .lt. 1d-8) exit
            lasttry = Yp_try 
        enddo

        if (bbn_type .eq. 1) then
            ombh2_BBN = (Neff_BBN/Neff_CMB)**0.75d0 * ombh2_CMB
        else if (bbn_type .eq. 2) then
            ombh2_BBN = ombh2_CMB
        endif

        return
    end subroutine solve_dilution


end module
