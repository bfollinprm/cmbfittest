!Parameterization using theta = r_s/D_a instead of H_0, and tau instead of z_re
!and log(A_s) instead of A_s
!Less general, but should give better performance
!
!The well-determined parameter A_s exp(-2tau) should be found by the covariance matrix
!parameter 3 is 100*theta, parameter 4 is tau, others same as params_H except A->log(A)
!Theta is much better constrained than H_0
!
!AL Jul 2005 - fixed bug which screwed up tau values for later importance sampling
!AL Feb 2004 - fixed compiler compatibility issue with **-number, typo in zstar
!AL Dec 2003 - renamed from params_Kowosky, changed to tau - log(A)
!AL Sept 2003 - fixed bug in declaration of dtauda routine
!AL June 2003
!Assumes prior 0.4 < h < 1

   function CMBToTheta(CMB)
     use settings
     use cmbtypes
     use ModelParams
     use CMB_Cls
     use Precision
     implicit none
     Type(CMBParams) CMB
     double precision zstar, astar, atol, rs,  DA
     double precision, external :: dsoundda, rombint
     real CMBToTheta
     integer error
  
     call InitCAMB(CMB,error,.false.)

!!From Hu & Sugiyama
       zstar =  1048*(1+0.00124*CMB%ombh2**(-0.738))*(1+ &
        (0.0783*CMB%ombh2**(-0.238)/(1+39.5*CMB%ombh2**0.763)) * &
           (CMB%omdmh2+CMB%ombh2)**(0.560/(1+21.1*CMB%ombh2**1.81)))
     
       astar = 1/(1+zstar)
       atol = 1e-6
       rs = rombint(dsoundda,1d-8,astar,atol)
       DA = AngularDiameterDistance(zstar)/astar
       CMBToTheta = rs/DA
!       print *,'z* = ',zstar, 'r_s = ',rs, 'DA = ',DA, rs/DA

  end function CMBToTheta




!Mapping between array of power spectrum parameters and CAMB
     subroutine SetCAMBInitPower(P,CMB,in)
       use camb
       use settings
       use cmbtypes
       implicit none
       type(CAMBParams)  P
       Type(CMBParams) CMB

       integer, intent(in) :: in


       if (Power_Name == 'power_tilt') then

       P%InitPower%k_0_scalar = pivot_k
       P%InitPower%k_0_tensor = pivot_k

       P%InitPower%ScalarPowerAmp(in) = cl_norm*CMB%norm(norm_As)
       P%InitPower%rat(in) = CMB%norm(norm_amp_ratio)
        
       P%InitPower%an(in) = CMB%InitPower(1)
       P%InitPower%ant(in) = CMB%InitPower(2)
       P%InitPower%n_run(in) = CMB%InitPower(3)
       if (inflation_consistency) then
         P%InitPower%ant(in) = - CMB%norm(norm_amp_ratio)/8.
          !note input n_T is ignored, so should be fixed (to anything)
       end if
       else
         stop 'params_CMB:Wrong initial power spectrum'
       end if

    end subroutine SetCAMBInitPower
 

 subroutine SetForH(Params,CMB,H0, firsttime)
     use settings
     use cmbtypes
     use CMB_Cls
     use bbn
     use abundance
     implicit none
     real Params(num_Params)
     logical, intent(in) :: firsttime
     Type(CMBParams) CMB
     real h2,H0
     real(dl) logDoH, yp

     logical    :: Yp_free = .false.

     character(1)   :: ridx
  
    CMB%H0=H0
    if (firsttime) then
        CMB%ombh2 = Params(1)    
        CMB%omdmh2 = Params(2)
        CMB%zre = Params(4) !!Not actually used.. is tau in this parameterization
        CMB%Omk = Params(5)
        CMB%nufrac = Params(6)
        CMB%w = Params(7)
     
        CMB%nnu = Params(8)
        
        if (bbn_type .eq. 0) then
            CMB%ombh2_bbn = Params(10)
            CMB%nnu_bbn = Params(11)
            
            Params(10) = CMB%ombh2_bbn
            Params(11) = CMB%nnu_bbn

            if (abs(CMB%ombh2_bbn) .lt. 1.0e-6 .and. &
            & abs(CMB%nnu_bbn) .lt. 1.0e-6) Yp_free = .true.

        else if (bbn_type .eq. 1) then
            CMB%nnu_bbn = Params(11)
            CMB%ombh2_bbn = (CMB%nnu_bbn/CMB%nnu)**0.75e0 * CMB%ombh2

            Params(10) = CMB%ombh2_bbn
            Params(11) = CMB%nnu_bbn
        else if (bbn_type .eq. 2) then
            CMB%ombh2_bbn = CMB%ombh2
            CMB%nnu_bbn = Params(11)

            Params(10) = 0.00e0
            Params(11) = CMB%nnu_bbn
        else if (bbn_type .eq. 3) then
            CMB%ombh2_bbn = CMB%ombh2
            CMB%nnu_bbn = CMB%nnu

            Params(10) = 0.00e0
            Params(11) = 0.00e0
        else
         !e.g. set from free parameter..
            write(*,*) "bbn_type set incorrectly. STOP!"
            stop
         !call MpiStop('params_CMB: YHe not free parameter in default parameterization')
        end if
!        Params(10) = CMB%ombh2_bbn
!        Params(11) = CMB%nnu_bbn
        
        if (Yp_free) then
            Yp = Params(9)
            logDoH = 0.00d0
        else
            if (CMB%nnu_bbn .lt. 10.0 .and. CMB%nnu_bbn .gt. 0.047 .and. &
            & CMB%ombh2_bbn .gt. 5.0e-3 .and. CMB%ombh2_bbn .lt. 4.0e-2) then
                Yp = yp_bbn(CMB%ombh2_bbn,CMB%nnu_bbn  - 3.046)
            else
                call deuterium_helium(dble(CMB%ombh2_bbn), dble(CMB%nnu_bbn), Yp, logDoH, bbn_fitting)
            endif
        endif

        CMB%YHe = Yp
        CMB%logDoH = logDoH
        Params(9) = CMB%YHe
        
        CMB%aeisw  = Params(12)
        CMB%alens  = Params(13)
        
        !write(*,*) Params
        !write(ridx, '(I1)') instance
        !open(unit=82+instance*10, file='Params_'//trim(chain_prefix)//'_'//ridx//'.txt', status='unknown')
        !write(82+instance*10, *) Params
        !close(82+instance*10)
        
        CMB%InitPower(1:num_initpower) = Params(index_initpower:index_initpower+num_initPower-1)
        CMB%norm(1) = exp(Params(index_norm))
        CMB%norm(2:num_norm) = Params(index_norm+1:index_norm+num_norm-1)
        CMB%nuisance(1:num_nuisance_params) = Params(index_nuisance:index_nuisance+num_nuisance_params-1)
    end if
    
    CMB%h = CMB%H0/100
    h2 = CMB%h**2
    CMB%omnuh2 = CMB%omdmh2*CMB%nufrac
    CMB%omch2 = CMB%omdmh2 - CMB%omnuh2
    CMB%omb = CMB%ombh2/h2
    CMB%omc = CMB%omch2/h2
    CMB%omnu = CMB%omnuh2/h2
    CMB%omdm = CMB%omdmh2/h2
    CMB%omv = 1- CMB%omk - CMB%omb - CMB%omdm

 end  subroutine SetForH

 subroutine ParamsToCMBParams(Params, CMB)
     use settings
     use cmbtypes
     use CMB_Cls
     use expansion

     implicit none
     real Params(num_params)
     real, save :: LastParams(num_params) = 0.
     real, save :: LastH0, Lastzre

     real(8)    :: z_eq, z_star, theta_s, theta_d

     character(1)   :: ridx

     Type(CMBParams) CMB
     real DA
     real  D_b,D_t,D_try,try_b,try_t, CMBToTheta, lasttry,tau
     external CMBToTheta

     !write(ridx, '(I1)') instance
     !open(unit=82+instance*10, file='Params_'//trim(chain_prefix)//'_'//ridx//'.txt', status='old', action='read')
     !read(82+instance*10, *) Params
     !close(82+instance*10)

     if (all(Params(1:num_hard) == Lastparams(1:num_hard))) then
       call SetForH(Params,CMB,LastH0, .true.)
       CMB%zre = Lastzre
       CMB%reserved(1) = params(4)
     else

     DA = Params(3)/100.00
     try_b = 20.0
     call SetForH(Params,CMB,try_b, .true.)
     D_b = CMBToTheta(CMB)
     try_t = 200.0
     call SetForH(Params,CMB,try_t, .false.)
     D_t = CMBToTheta(CMB)

     if (CMB%YHe .lt. 0.0 .or. CMB%YHe .gt. 1.000) then
         CMB%YHe = -1.00
         CMB%H0 = 0.00
         CMB%reserved = 0
         CMB%reserved(1) = params(4)
         return
     endif

     if (DA < D_b .or. DA > D_t) then
      cmb%H0=0 !Reject it
     else
     lasttry = -1
     do
            call SetForH(Params,CMB,(try_b+try_t)/2, .false.)
            D_try = CMBToTheta(CMB)
               if (D_try < DA) then
                  try_b = (try_b+try_t)/2
               else
                  try_t = (try_b+try_t)/2
               end if
               if (abs(D_try - lasttry)< 1e-7) exit
              lasttry = D_try
     end do

    !!call InitCAMB(CMB,error)
    tau = params(4)
    CMB%zre = GetZreFromTau(CMB, tau)
    if (CMB%zre .lt. 0) then
        cmb%H0 = 0 !! reject it
        return
    endif    
  
    LastH0 = CMB%H0
    Lastzre = CMB%zre
    LastParams = Params
    end if

    CMB%reserved = 0
    CMB%reserved(1) = params(4) !tau
  
     end if
    
    if (cmb%H0 .ge. 20.0 .and. cmb%H0 .le. 200.0) then
        call SetTEParams(CMB, .true.)
        call z_equality(z_eq)
        call z_recomb(z_star)
        theta_s = Theta_Sound(z_star)
        theta_d = Theta_Damp(z_star)
    
        CMB%sum_mnu = sum_mnu_ev(dble(CMB%nnu))
        CMB%z_star = z_star
        CMB%z_eq = z_eq
        CMB%thetaS = theta_s
        CMB%thetaD = theta_d
        CMB%DoS    = theta_d/theta_s
        !print*, '\sum m_nu = ', CMB%sum_mnu, ' eV.'
    endif
 
   end subroutine ParamsToCMBParams

   subroutine CMBParamsToParams(CMB, Params)
     use settings
     use cmbtypes
     implicit none
     real Params(num_Params)
     Type(CMBParams) CMB
     real CMBToTheta
     external CMBToTheta
 
      Params(1) = CMB%ombh2 
      Params(2) = CMB%omdmh2
 
      Params(3) = CMBToTheta(CMB)*100
      Params(4) = CMB%reserved(1)
      Params(5) = CMB%omk 
      
      Params(6) = CMB%nufrac 
      Params(7) = CMB%w
      Params(8) = CMB%nnu
      Params(9) = CMB%YHe
      Params(10)= CMB%ombh2_bbn
      Params(11)= CMB%nnu_bbn
      Params(12)= CMB%aeisw
      Params(13)= CMB%alens

      Params(index_initpower:index_initpower+num_initpower-1) =CMB%InitPower(1:num_initpower) 
      Params(index_norm) = log(CMB%norm(1))
      Params(index_norm+1:index_norm+num_norm-1) = CMB%norm(2:num_norm)
      Params(index_nuisance:index_nuisance+num_nuisance_params-1)=CMB%nuisance(1:num_nuisance_params) 

   end subroutine CMBParamsToParams

   subroutine SetParamNames(Names)
    use settings
    use ParamNames
    Type(TParamNames) :: Names
 
    if (ParamNamesFile /='') then
      call ParamNames_init(Names, ParamNamesFile)
    else
     if (generic_mcmc) then
      Names%nnames=0
      if (Feedback>0) write (*,*) 'edit SetParamNames in params_CMB.f90 if you want to use named params'
     else
       call ParamNames_init(Names, trim(LocalDir)//'params_CMB.paramnames')
    end if
    end if
   end subroutine SetParamNames

 
  function CalcDerivedParams(P, derived) result (num_derived)
     use settings
     use cmbtypes
     use ParamDef
     use Lists
     implicit none
     Type(real_pointer) :: derived
     Type(ParamSet) P
     Type(CMBParams) CMB
     real r10   
     integer num_derived 
     
     num_derived = 15
   
     allocate(Derived%P(num_derived))
   
      call ParamsToCMBParams(P%P,CMB)

      if (lmax_tensor /= 0 .and. compute_tensors) then
          r10 = P%Info%Theory%cl_tensor(10,1)/P%Info%Theory%cl(10,1)
      else
        r10 = 0
      end if

      derived%P(1) = CMB%omv
      derived%P(2) = P%Info%Theory%Age
      derived%P(3) = CMB%omdm+CMB%omb
      derived%P(4) = CMB%sum_mnu
      derived%P(5) = P%Info%Theory%Sigma_8      
      derived%P(6) = CMB%zre
      derived%P(7) = r10
      derived%P(8) = CMB%H0
      derived%P(9) = CMB%logDoH
      derived%P(10) = CMB%z_eq
      derived%P(11)= CMB%z_star
      derived%P(12)= CMB%thetaS
      derived%P(13)= CMB%thetaD
      derived%P(14)= CMB%DoS
      derived%P(15)= P%Info%Theory%cl(100,num_clsS+1)
      
!      print*, "P15 = ", derived%P(15)/((101.0**2/100.0**2)/twopi)*(2725000.0)**2

  end function CalcDerivedParams
  

  subroutine WriteParams(P, mult, like)
     use settings
     use cmbtypes
     use ParamDef
     use IO
     use Lists
     implicit none
     Type(ParamSet) P
     real, intent(in) :: mult, like
     Type(CMBParams) CMB
     real, allocatable :: output_array(:)
     Type(real_pointer) :: derived
     integer numderived 
     integer CalcDerivedParams
     external CalcDerivedParams
  
    if (outfile_handle ==0) return
  
    if (generic_mcmc) then

      call IO_OutputChainRow(outfile_handle, mult, like, P%P)
     
    else
    
      numderived = CalcDerivedParams(P, derived)

      allocate(output_array(num_real_params + numderived + nuisance_params_used ))
      output_array(1:num_real_params) =  P%P(1:num_real_params)
      output_array(num_real_params+1:num_real_params+numderived) =  derived%P
      deallocate(derived%P)

      if (nuisance_params_used>0) then
       output_array(num_real_params+numderived+1:num_real_params+numderived+nuisance_params_used) = &
        P%P(num_real_params+1:num_real_params+nuisance_params_used) 
      end if
 
      call IO_OutputChainRow(outfile_handle, mult, like, output_array)
      deallocate(output_array)           
    end if

  end  subroutine WriteParams




  subroutine WriteParamsAndDat(P, mult, like)
     use settings
     use cmbtypes
     use ParamDef
     use IO
     use Lists
     implicit none
     Type(ParamSet) P
     real, intent(in) :: mult, like
     character(LEN =30) fmt
     Type(CMBParams) CMB
     real,allocatable :: output_array(:)
     Type(real_pointer) :: derived
     integer numderived 
     integer CalcDerivedParams
     external CalcDerivedParams
         
    if (outfile_handle ==0) return

      numderived = CalcDerivedParams(P, derived)

      allocate(output_array(num_real_params + numderived + num_matter_power ))
      output_array(1:num_real_params) =  P%P(1:num_real_params)
      output_array(num_real_params+1:num_real_params+numderived) =  derived%P
      deallocate(derived%P)

      output_array(num_real_params+numderived+1:num_real_params+numderived+num_matter_power) = &
        P%Info%Theory%matter_power(:,1) 

      call IO_OutputChainRow(outfile_handle, mult, like, output_array)
      deallocate(output_array)

  end  subroutine WriteParamsAndDat


  function dsoundda(a)
          use Precision
          use ModelParams
     
          implicit none
          real(dl) dsoundda,dtauda,a,R,cs
          external dtauda

           R=3.0d4*a*CP%omegab*(CP%h0/100.0d0)**2
           cs=1.0d0/sqrt(3*(1+R))
           dsoundda=dtauda(a)*cs
        
  end function dsoundda

