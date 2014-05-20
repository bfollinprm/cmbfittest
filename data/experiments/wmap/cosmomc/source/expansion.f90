module massive_neutrinos
    use settings
    use precision
    use constants
    
    implicit none
    private 
        
    real(dl), parameter  :: const  = 7._dl/120*pi**4 ! 5.68219698_dl
    !const = int q^3 F(q) dq = 7/120*pi^4
    real(dl), parameter  :: const2 = 5._dl/7/pi**2   !0.072372274_dl
    real(dl), parameter  :: zeta3  = 1.2020569031595942853997_dl
    real(dl), parameter  :: zeta5  = 1.0369277551433699263313_dl
    real(dl), parameter  :: zeta7  = 1.0083492773819228268397_dl

    integer, parameter  :: nrhopn=2000  
    real(dl), parameter :: am_min = 0.01_dl  !0.02_dl
    !smallest a*m_nu to integrate distribution function rather than using series
    real(dl), parameter :: am_max = 600._dl 
    !max a*m_nu to integrate
          
    real(dl),parameter  :: am_minp=am_min*1.1
    real(dl), parameter :: am_maxp=am_max*0.9
   
    real(dl) dlnam, mass_nu

    real(dl), dimension(:), allocatable ::  r1,p1,dr1,dp1,ddr1

    !Sample for massive neutrino momentum
    !These settings appear to be OK for P_k accuate at 1e-3 level
    integer, parameter :: nqmax0=80 !maximum array size of q momentum samples 
    real(dl) :: nu_q(nqmax0), nu_int_kernel(nqmax0)
 
    integer nqmax !actual number of q modes evolves
 
    public massive_nu_init, mass_nu, nu_density

    contains

    subroutine massive_nu_init(omnu, nu_massive, H0)
      
!  Initialize interpolation tables for massive neutrinos.
!  Use cubic splines interpolation of log rhonu and pnu vs. log a*m.
        
        real(dl), intent(in) :: omnu, nu_massive, H0
        integer i
        real(dl) dq,dlfdlq, q, am, rhonu,pnu
        real(dl) grhom, grhog, grhor
        real(dl) spline_data(nrhopn)
     
!  nu_masses=m_nu(i)*c**2/(k_B*T_nu0).
!  Get number density n of neutrinos from
!  rho_massless/n = int q^3/(1+e^q) / int q^2/(1+e^q)=7/180 pi^4/Zeta(3)
!  then m = Omega_nu/N_nu rho_crit /n
!  Error due to velocity < 1e-5
        
        !do i=1, CP%Nu_mass_eigenstates 
        !    nu_masses(i) = const/(1.5d0*zeta3)*grhom/grhor*CP%omegan*CP%Nu_mass_fractions(i) &
        !                   /CP%Nu_mass_degeneracies(i)
        !end do
        
        grhom = 3*H0**2/c**2*1000**2
        grhog = kappa/c**2*4*sigma_boltz/c**3*2.725_dl**4*Mpc**2
        grhor = 7.0_dl/8.0_dl*(4.0_dl/11.0_dl)**(4.0_dl/3.0_dl)*grhog
        mass_nu = const/(1.5d0*zeta3)*grhom/grhor * omnu / nu_massive

        if (allocated(r1)) return
        allocate(r1(nrhopn),p1(nrhopn),dr1(nrhopn),dp1(nrhopn),ddr1(nrhopn))

        nqmax=3
        if (AccuracyLevel >1) nqmax=4
        if (AccuracyLevel >2) nqmax=5
        if (AccuracyLevel >3) nqmax=nint(AccuracyLevel*10) 
          !note this may well be worse than the 5 optimized points

        if (nqmax > nqmax0) call MpiStop('massive_nu_init: qmax > nqmax0')

        !We evolve evolve 4F_l/dlfdlq(i), so kernel includes dlfdlnq factor
        !Integration scheme gets (Fermi-Dirac thing)*q^n exact,for n=-4, -2..2
        !see CAMB notes
        if (nqmax==3) then
          !Accurate at 2e-4 level
          nu_q(1:3) = (/0.913201, 3.37517, 7.79184/)
          nu_int_kernel(1:3) = (/0.0687359, 3.31435, 2.29911/)
          
        else if (nqmax==4) then
          !This seems to be very accurate (limited by other numerics)
           nu_q(1:4) = (/0.7, 2.62814, 5.90428, 12.0/)
           nu_int_kernel(1:4) = (/0.0200251, 1.84539, 3.52736, 0.289427/)
      
        else if (nqmax==5) then
        !exact for n=-4,-2..3 
        !This seems to be very accurate (limited by other numerics)
         nu_q(1:5) = (/0.583165, 2.0, 4.0, 7.26582, 13.0/)  
         nu_int_kernel(1:5) = (/0.0081201, 0.689407, 2.8063, 2.05156, 0.126817/) 
  
        else
         dq = (12 + nqmax/5)/real(nqmax)
         do i=1,nqmax
            q=(i-0.5d0)*dq
            nu_q(i) = q 
            dlfdlq=-q/(1._dl+exp(-q))
            nu_int_kernel(i)=dq*q**3/(exp(q)+1._dl) * (-0.25_dl*dlfdlq) !now evolve 4F_l/dlfdlq(i)
            
         end do
        end if
        nu_int_kernel=nu_int_kernel/const
        
        dlnam=-(log(am_min/am_max))/(nrhopn-1)
 

        !$OMP PARALLEL DO DEFAULT(SHARED),SCHEDULE(STATIC) &
        !$OMP & PRIVATE(am, rhonu,pnu) 
        do i=1,nrhopn
          am=am_min*exp((i-1)*dlnam)
          call nuRhoPres(am,rhonu,pnu)
          r1(i)=log(rhonu)
          p1(i)=log(pnu)
        end do
        !$OMP END PARALLEL DO


        call splini(spline_data,nrhopn)
        call splder(r1,dr1,nrhopn,spline_data)
        call splder(p1,dp1,nrhopn,spline_data)
        call splder(dr1,ddr1,nrhopn,spline_data)

    end subroutine massive_nu_init

    subroutine nuRhoPres(am,rhonu,pnu)
!  Compute the density and pressure of one eigenstate of massive neutrinos,
!  in units of the mean density of one flavor of massless neutrinos.

        real(dl),  parameter :: qmax=30._dl
        integer, parameter :: nq=100
        real(dl) dum1(nq+1),dum2(nq+1)
        real(dl), intent(in) :: am
        real(dl), intent(out) ::  rhonu,pnu
        integer i
        real(dl) q,aq,v,aqdn,adq
      

!  q is the comoving momentum in units of k_B*T_nu0/c.
!  Integrate up to qmax and then use asymptotic expansion for remainder.
        adq=qmax/nq
        dum1(1)=0._dl
        dum2(1)=0._dl
        do  i=1,nq
          q=i*adq
          aq=am/q
          v=1._dl/sqrt(1._dl+aq*aq)
          aqdn=adq*q*q*q/(exp(q)+1._dl)
          dum1(i+1)=aqdn/v
          dum2(i+1)=aqdn*v
        end do
        call splint(dum1,rhonu,nq+1)
        call splint(dum2,pnu,nq+1)
!  Apply asymptotic corrrection for q>qmax and normalize by relativistic
!  energy density.
        rhonu=(rhonu+dum1(nq+1)/adq)/const
        pnu=(pnu+dum2(nq+1)/adq)/const/3._dl   
    end subroutine nuRhoPres

    subroutine nu_density(am,rhonu)
        use precision
        !use ModelParams
        real(dl), intent(in) :: am
        real(dl), intent(out) :: rhonu

!  Compute massive neutrino density in units of the mean
!  density of one eigenstate of massless neutrinos.  Use cubic splines to
!  interpolate from a table.

        real(dl) d
        integer i
      
        if (am <= am_minp) then
          rhonu=1._dl + const2*am**2  
          return
        else if (am >= am_maxp) then
          rhonu = 3/(2*const)*(zeta3*am + (15*zeta5)/2/am)
          return
        end if
        
        d=log(am/am_min)/dlnam+1._dl
        i=int(d)
        d=d-i
       
!  Cubic spline interpolation.
        rhonu=r1(i)+d*(dr1(i)+d*(3._dl*(r1(i+1)-r1(i))-2._dl*dr1(i) &
               -dr1(i+1)+d*(dr1(i)+dr1(i+1)+2._dl*(r1(i)-r1(i+1)))))
        rhonu=exp(rhonu)
    end subroutine nu_density

end module massive_neutrinos



module expansion

    use constants
    use precision
    use settings
    use Recombination

    use cmbtypes
    use massive_neutrinos

    implicit none

    private

    real(8),    parameter   :: Tcmb = 2.725d0

    !! parameters relevant to thermo-expansion *AFTER* BBN
    type ThermExp_Params
		real(8)     :: ombh2, omb
		real(8)     :: omch2, omc
        real(8)     :: omnuh2, omnu
		real(8)     :: omk
		real(8)     :: omv, w
		real(8)     :: H0
		real(8)     :: neff
		real(8)     :: Yp
        real(8)     :: nu_massive
        real(8)     :: nu_massless
        real(8)     :: m_nu
    end	type ThermExp_Params

    type(ThermExp_Params),      private :: Params
    type(RecombinationParams),  private :: Recomb
    
    real(8)     :: w_lam = -1.00_dl
    real(8)     :: grhom_ex, grhog_ex, grhor_ex, grhoc_ex, grhob_ex, grhov_ex, grhok_ex
    real(8)     :: grhonu_massless, grhonu_massive
	real(8)     :: N0, athom
    
    public SetTEParams, z_equality, z_recomb, Theta_Sound, Theta_Damp, sum_mnu_ev

    contains

    subroutine SetTEParams(CMB, first) 

        implicit none

        type(CMBParams),    intent(in)      :: CMB
        logical,            intent(in)      :: first

        Params%ombh2  = CMB%ombh2
        Params%omch2  = CMB%omch2
        Params%omnuh2 = CMB%omnuh2
        Params%omk    = CMB%omk
        Params%w      = CMB%w
        Params%H0     = CMB%H0
        Params%neff   = CMB%nnu
        Params%Yp     = CMB%YHe
        
        if (CMB%omnu .eq. 0.0_dl) then
            Params%nu_massless = Params%neff
            Params%nu_massive  = 0.0_dl
        else
            Params%nu_massive  = Params%neff
            Params%nu_massless = 0.0_dl
        endif

        w_lam = CMB%w

        Params%omb   = Params%ombh2/(Params%H0/100.0_dl)**2.0_dl
        Params%omc   = Params%omch2/(Params%H0/100.0_dl)**2.0_dl
        Params%omnu  = Params%omnuh2/(Params%H0/100.0_dl)**2.0_dl
        Params%omv   = 1.00_dl - Params%omb - Params%omc - Params%omnu - Params%omk

        grhom_ex = 3*Params%H0**2/c**2*1000**2
        grhog_ex = kappa/c**2*4*sigma_boltz/c**3*Tcmb**4*Mpc**2
        grhor_ex = 7.0_dl/8.0_dl*(4.0_dl/11.0_dl)**(4.0_dl/3.0_dl)*grhog_ex
        grhoc_ex = grhom_ex*Params%omc
        grhob_ex = grhom_ex*Params%omb
        grhov_ex = grhom_ex*Params%omv
        grhok_ex = grhom_ex*Params%omk

        grhonu_massless = Params%nu_massless * grhor_ex
        grhonu_massive  = Params%nu_massive * grhor_ex

        mass_nu = 0.0_dl
        if (grhonu_massive .ne. 0) then
            call massive_nu_init(Params%omnu, Params%nu_massive, Params%H0)
        endif
        Params%m_nu = 1.68e-4*mass_nu !! eV mass of each species
        
        if (first) then
            call Recombination_SetDefParams(Recomb)
            N0 = Params%omb * (1.00_dl - Params%Yp) * grhom_ex*c**2/kappa/m_H/Mpc**2
            athom = sigma_thomson*N0*Mpc 
        endif

        return
    end subroutine SetTEParams


    function sum_mnu_ev(nu_massive)
        
        real(8),    intent(in)  :: nu_massive
        real(8)                 :: sum_mnu_ev

        sum_mnu_ev = nu_massive * Params%m_nu
        
        return
    end function sum_mnu_ev


    subroutine z_equality(zEQ)

		real(8),    intent(out) :: zEQ

		real(8)     :: bigH, H0p, fnu
		real(8)     :: foe
		real(8)     :: omega_radh2, ommh2, omnuh2
        real(8)     :: grhorad, grhomat
		real(8)     :: aEQ
        real(8)     :: z1, z2, avg, diff, ratio
        integer     :: i

        !bigH= 100.0D3/Mpc
        !H0p = Params%H0/100._dl*bigH
        omnuh2 = Params%omnuh2

        !fnu = Params%neff*(7.d0/8.d0)*(4.d0/11.d0)**(4.d0/3.d0)
        !omega_radh2 = (8.d0*Pi*G*a_rad*(1.d0+fnu)*Tcmb**4)/(3.d0*(c*bigH)**2)
        !ommh2 = Params%ombh2 + Params%omch2
        !aEQ = omega_radh2/ommh2
        !zEQ = 1.00d0/aEQ - 1.00d0
        !print*, "old       zeq = ", zEQ

        if (omnuh2 .eq. 0) then
            grhorad = grhog_ex + grhonu_massless
            aEQ = grhorad / (grhob_ex + grhoc_ex)
            zEQ = 1.00d0/aEQ - 1.00d0
            !print*, "massless: zeq = ", zEQ
        else
            z1 = 0.0d0
            z2 = 10000.d0

            i = 0
            diff = 10.0d0
            do while (diff .gt. 1d-8)
                i = i + 1
                if (i .eq. 200) stop 'redshift finder did not converge'

                diff = ratio_radmat_massivenu(z2) - ratio_radmat_massivenu(z1)
                avg = 0.5d0*(z2 + z1)
                ratio = ratio_radmat_massivenu(avg)
                if (ratio .gt. 1.d0) then
                    z2 = avg
                else
                    z1 = avg
                end if

                !print*, ratio
            enddo

            zEQ = avg
            !print*, "massive:  zeq = ", zEQ
        endif

        return
    end subroutine z_equality

    
    function ratio_radmat_massivenu(z)
        real(8),    intent(in)     :: z
        real(8)     :: ratio_radmat_massivenu

        real(8)     :: a, rhonu
        real(8)     :: grhomat, grhorad

        a = 1.0d0/(1.0d0+z)
        grhomat = (grhob_ex + grhoc_ex) * a
        call nu_density(a*mass_nu, rhonu)
        grhorad = rhonu*grhonu_massive + grhog_ex + grhonu_massless

        ratio_radmat_massivenu = grhorad / grhomat

        return
    end function


    subroutine z_recomb(zstar)

		real(8),    intent(out) :: zstar
		real(8)     :: try1,try2,diff,avg
        integer     :: i
        
        zstar = 0.0d0
        
        call Recombination_init(Recomb, Params%omc, Params%omb, 0.0d0, Params%omv, &
        & Params%H0, Tcmb, Params%Yp, Params%neff)
!		call Recombination_init(Recomb, Pa%OmegaC, Pa%OmegaB, 0.0d0, Pa%OmegaV, Pa%H0, &
!		& Tcmb, Pa%YHe)

        try1 = 0.d0
        try2 = 10000.d0

        i = 0
        diff = 10.d0
        do while (diff .gt. 1d-8)
            i=i+1
            if (i .eq. 100) stop 'optical depth redshift finder did not converge'

            diff = tau(try2) - tau(try1)
            avg = 0.5d0*(try2+try1)
            if (tau(avg) .gt. 1.d0) then
                try2 = avg
            else
                try1 = avg
            end if
        end do

        zstar = avg

        return
    end subroutine z_recomb


    function dtau_dz(z)

		real(dl) :: dtau_dz
		real(dl), intent(in) :: z
		real(dl) :: a

        a = 1._dl/(1._dl+z)
        !ignoring reionisation, not relevant for distance measures
        dtau_dz = Recombination_xe(a) * athom * deta_da(a)

        return
    end function dtau_dz


    function tau(z)
		real(dl) :: rombint2
		real(dl) tau
		real(dl),intent(in) :: z

        tau = rombint2(dtau_dz, 0.d0, z, 1d-6, 20, 100)

        return
    end function tau


    function deta_da(a)

		real(8),    intent(in)  :: a
		real(8)     :: deta_da
		real(8)     :: a2, grhoa2
        real(8)     :: rhonu

        a2 = a**2.d0
        grhoa2 = grhok_ex*a2 + (grhoc_ex+grhob_ex)*a + grhog_ex + grhonu_massless
        if (w_lam .eq. -1.00d0) then
            grhoa2 = grhoa2 + grhov_ex*a2**2
        else
            grhoa2 = grhoa2 + grhov_ex*a**(1.0d0-3.0d0*w_lam)
        end if

        if (grhonu_massive .ne. 0) then
            call nu_density(a*mass_nu, rhonu)
            grhoa2 = grhoa2 + rhonu*grhonu_massive
        endif

        !grhoa2 = grhoa2 + grhov*a2**2.d0
        
        deta_da = sqrt(3.d0/grhoa2)
        
        return
    end function deta_da


    function dRs_da(a)

		real(8),    intent(in)  :: a
		real(8)     :: dRs_da
		real(8)     :: R, cs

        R = 3*grhob_ex*a / (4*grhog_ex)
        cs = 1.0d0/sqrt(3.0d0*(1.0d0+R))
        dRs_da = deta_da(a) * cs
        
        return
    end function dRs_da

    
    function Theta_Sound(zstar, Rs_out, Da_out)

		real(8),    intent(in)  :: zstar
		real(8),    intent(out),optional    :: Rs_out, Da_out

		real(8)     :: astar, Theta_Sound
		real(8)     :: Rs, Da, rombint

        astar = 1.0d0/(zstar+1.0d0)
        Rs = rombint(dRs_da,1d-8,astar,1d-7)
        Da = rombint(deta_da,astar,1.0d0,1.0d-7)

        if (present(Rs_out)) Rs_out = Rs
        if (present(Da_out)) Da_out = Da

        Theta_Sound = Rs/Da

        return
    end function Theta_Sound


    function Theta_Damp(zstar, Rd_out, Da_out)

		real(8),    intent(in)  :: zstar
		real(8),    intent(out),optional    :: Rd_out, Da_out

		real(8)     :: Theta_Damp
		real(8)     :: rombint
		real(8)     :: Rd2, Rd, astar, Da

        astar = 1.0d0/(1.0d0+zstar)
        Rd2 = rombint(dRd2_da, 1d-8, astar,1d-7)
        Rd = Pi*sqrt(Rd2)
        Da = rombint(deta_da,astar,1.0d0,1.0d-7)

        if (present(Rd_out)) Rd_out = Rd
        if (present(Da_out)) Da_out = Da

        Theta_Damp = Rd/Da

        return
    end function Theta_Damp


    function dRd2_da(a)

        implicit none

		real(8)     :: a
		real(8)     :: xe
		real(8)     :: R
		real(8)     :: dRd2_da

        R = 3*grhob_ex*a / (4*grhog_ex)

        xe = Recombination_xe(a)
        dRd2_da = deta_da(a)/(xe*athom/a**2.0d0)
        dRd2_da = dRd2_da/(6.d0*(1.d0+R)**2.d0)*(R**2.d0+16.d0/15.d0*(1.d0+R))

        return    
    end function dRd2_da

end module expansion



