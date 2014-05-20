module CalcLike
 use CMB_Cls
 use cmbtypes
 use cmbdata
 use mpk
 use Random
 use settings
 use ParamDef
 use snovae
 use WeakLen
 use Lya
 use abundance
 use SNLS  !! load SNLS module - spt.lps12

 implicit none

 logical :: Use_Age_Tophat_Prior = .true.
 logical :: Use_CMB = .true.
 logical :: Use_BBN = .false.
 logical :: Use_Clusters = .false.

 logical :: use_SZ_prior = .false.
 logical :: use_PS_prior = .false.
 logical :: use_CL_prior = .false.
 
 logical :: use_H0_prior  = .false.
 logical :: use_DoH_prior = .false.
 logical :: use_Yp_prior  = .false.
 logical :: use_DoS_prior = .false.
 logical :: use_tau_prior = .false.
 
 real   :: H0_mean, H0_stdev
 real   :: Yp_mean, Yp_stdev
 real   :: DoH_mean, DoH_stdev
 real   :: DoS_mean, DoS_stdev
 real   :: tau_mean, tau_stdev
 real   :: sz_mean, sz_stdev
 real   :: ps_mean, ps_stdev
 real   :: cl_mean, cl_stdev

 integer :: H0_min = 20, H0_max = 200
 real :: Omk_min = -0.3, Omk_max = 0.3
 real :: Use_min_zre = 0
 integer :: Age_min = 10, Age_max = 20
 real :: Temperature  = 1

 real, parameter :: sz_mean_default = 5.50,  sz_stdev_default = 3.00
 real, parameter :: ps_mean_default = 19.30, ps_stdev_default = 3.50
 real, parameter :: cl_mean_default = 5.00,  cl_stdev_default = 2.50

contains

  function GenericLikelihoodFunction(Params) 
    type(ParamSet)  Params 
    real :: GenericLikelihoodFunction
   
  
   !Used when you want to plug in your own CMB-independent likelihood function:
   !set generic_mcmc=.true. in settings.f90, then write function here returning -Ln(Likelihood)
   !Parameter array is Params%P
    !!!!
   ! GenericLikelihoodFunction = greatLike(Params%P)
   ! GenericLikelihoodFunction = LogZero 
    call MpiStop('GenericLikelihoodFunction: need to write this function!')
    GenericLikelihoodFunction=0

  end function GenericLikelihoodFunction

  
  function GetLogPrior(CMB, Info) !Get -Ln(Prior)
    real GetLogPrior
    real Age
    Type (CMBParams) CMB
    Type(ParamSetInfo) Info

    GetLogPrior = logZero
 
    if (.not. generic_mcmc) then
     if (CMB%H0 < H0_min .or. CMB%H0 > H0_max) return
     if (CMB%omk < Omk_min .or. CMB%omk > Omk_max .or. CMB%Omv < 0) return
     if (CMB%zre < Use_min_zre) return
     if (CMB%YHe < 0.00) return

     Age = GetAge(CMB, Info)
      !This sets up parameters in CAMB, so do not delete unless you are careful!
 
     if (Use_Age_Tophat_Prior .and. (Age < Age_min .or. Age > Age_max) .or. Age < 0) return
    
    end if
    GetLogPrior = 0
 
  end function GetLogPrior

  function GetLogLike(Params) !Get -Ln(Likelihood)
    type(ParamSet)  Params 
    Type (CMBParams) CMB
    real GetLogLike
    real dum(1,1)
 
    if (any(Params%P > Scales%PMax) .or. any(Params%P < Scales%PMin)) then
       GetLogLike = logZero
        return
    end if

    if (generic_mcmc) then
        GetLogLike = GenericLikelihoodFunction(Params) 
        if (GetLogLike /= LogZero) GetLogLike = GetLogLike/Temperature

    else

     call ParamsToCMBParams(Params%P,CMB)

     GetLogLike = GetLogLikePost(CMB, Params%Info,dum,.false.)
    end if 
   end function GetLogLike

    
  function GetLogLikePost(CMB, Info, inCls, HasCls) 
    use cambmain, only: initvars
    use Camb, only: CAMBParams_Set 
    type(CAMBParams)  P
    real GetLogLikePost
    Type (CMBParams) CMB
    Type(ParamSetInfo) Info
    real, intent(in):: inCls(:,:)
    logical, intent(in) :: HasCls
    real acl(lmax,num_cls_tot)
    integer error

    if (generic_mcmc) stop 'GetLogLikePost: not supported for generic'

    GetLogLikePost  = GetLogPrior(CMB, Info)
    if ( GetLogLikePost >= logZero) then
       GetLogLikePost = logZero
       
    else 
       GetLogLikePost = GetLogLikePost + sum(CMB%nuisance(1:nuisance_params_used)**2)/2
          !Unit Gaussian prior on all nuisance parameters
       if (Use_BBN) GetLogLikePost = GetLogLikePost + (CMB%ombh2 - 0.022)**2/(2*0.002**2) 
          !I'm using increased error bars here
   
       if (Use_CMB .or. Use_LSS) then
          if (HasCls) then
           acl = inCls
           error =0
          else
           call GetCls(CMB, Info, acl, error)
          end if
         if (error /= 0) then
          GetLogLikePost = logZero 
         else
          if (Use_CMB) then
            GetLogLikePost = &
            CMBLnLike(acl, CMB%norm(norm_freq_ix:norm_freq_ix+num_freq_params-1),CMB%nuisance) + GetLogLikePost
            if (use_SZ_prior) GetLogLikePost = GetLogLikePost + (((CMB%norm(norm_freq_ix)-sz_mean)/sz_stdev)**2)/2
            if (use_PS_prior) GetLogLikePost = GetLogLikePost + (((CMB%norm(norm_freq_ix+2)-ps_mean)/ps_stdev)**2)/2
            if (use_CL_prior) GetLogLikePost = GetLogLikePost + (((CMB%norm(norm_freq_ix+3)-cl_mean)/cl_stdev)**2)/2

            if (use_Yp_prior) GetLogLikePost = GetLogLikePost + Yp_LnLike(CMB, Yp_mean, Yp_stdev)
            if (use_DoH_prior) GetLogLikePost = GetLogLikePost + DoH_LnLike(CMB, DoH_mean, DoH_stdev)
            if (use_DoS_prior) GetLogLikePost = GetLogLikePost + (CMB%DoS - DoS_mean)**2.00e0/(2.00e0*DoS_stdev**2.00e0)

            if (use_H0_prior) GetLogLikePost = GetLogLikePost + (CMB%H0 - H0_mean)**2.00e0/(2.00e0*H0_stdev**2.00e0)
            if (use_tau_prior) GetLogLikePost = GetLogLikePost + (CMB%reserved(1) - tau_mean)**2.00e0/(2.00e0*tau_stdev**2.00e0)
          endif
          if (Use_mpk) GetLogLikePost = GetLogLikePost + LSSLnLike(CMB, Info%theory)
          if (Use_WeakLen) GetLogLikePost = GetLogLikePost + WeakLenLnLike(CMB, Info%theory)     
          if (Use_Lya) GetLogLikePost = GetLogLikePost +  LSS_Lyalike(CMB, Info%Theory)
          if ( GetLogLikePost >= logZero) then
            GetLogLikePost = logZero
          end if
         end if
         if (Use_SN .and. GetLogLikePost /= logZero ) then
            if (Info%Theory%SN_loglike /= 0) then
             GetLogLikePost = GetLogLikePost + Info%Theory%SN_loglike
            else
             GetLogLikePost = GetLogLikePost + SN_LnLike(CMB)
            end if
               !Assume computed only every time hard parameters change
         end if

         !! adding SNLS likelihood - spt.lps12
         IF (Use_SNLS .AND. GetLogLikePost /= LogZero) THEN
            GetLogLikePost = GetLogLikePost + snls_LnLike( CMB )
         ENDIF

         if (Use_BAO .and. GetLogLikePost /= logZero ) then
            if (Info%Theory%BAO_loglike /= 0) then
             GetLogLikePost = GetLogLikePost + Info%Theory%BAO_loglike
            else
             GetLogLikePost = GetLogLikePost + BAO_LnLike(CMB)
            end if
               !Assume computed only every time hard parameters change
         end if
         if (Use_HST .and. GetLogLikePost /= logZero) then
            if (Info%Theory%HST_loglike /= 0) then
             GetLogLikePost = GetLogLikePost + Info%Theory%HST_loglike
            else
             GetLogLikePost = GetLogLikePost + HST_LnLike(CMB)
            end if
           !!Old: GetLogLikePost = GetLogLikePost + (CMB%H0 - 72)**2/(2*8**2)  !HST 
         end if  
     
         else
           if (Use_SN) GetLogLikePost = GetLogLikePost + SN_LnLike(CMB)
           if (Use_BAO)then
                !From Jason Dosset           
                !JD had to change below so new BAO will work without CMB
                !previous way z_drag was not calculated without calling get_cls.
                if (RecomputeTransfers(CMB, Info%LastParams))  then
                    call CMBToCAMB(CMB, P)
                    call CAMBParams_Set(P)
                    call InitVars
                    Info%Theory%BAO_loglike = Bao_lnLike(CMB)
                    Info%LastParams = CMB
                end if
                GetLogLikePost = GetLogLikePost + Info%Theory%BAO_loglike
           end if
           if (Use_HST) GetLogLikePost = GetLogLikePost + HST_LnLike(CMB)
           IF (Use_SNLS) GetLogLikePost = GetLogLikePost + snls_LnLike(CMB)  !! adding SNLS likelihood - spt.lps12
         end if
         
     
      if (Use_Clusters .and. GetLogLikePost /= LogZero) then
          GetLogLikePost = GetLogLikePost + &
                 (Info%Theory%Sigma_8-0.9)**2/(2*0.05**2)
          stop 'Write your cluster prior in calclike.f90 first!'
      end if

     if (GetLogLikePost /= LogZero) GetLogLikePost = GetLogLikePost/Temperature
   

    end if

  end function GetLogLikePost

end module CalcLike
