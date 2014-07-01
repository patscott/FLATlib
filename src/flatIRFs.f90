! Fermi LAT Instrumental Response Function functional routines 
!
! These return the Fermi LAT IRFs as a function of true and 
! observed photon energies and incidence angles in the LAT. 
!
! Pat Scott, Feb 2009; pat@fysik.su.se
! Modified:	Pat Scott, Oct 2011	Added support for Pass7-type PSFs
!--------------------------------------------------------------------------------


      MODULE flatIRFs

      use flatCommon
      use flatUtils
      use flatIRFparams

      implicit none

      contains

    
      RECURSIVE FUNCTION flatIRFs_Aeff(l10E_true, IncidenceAngles_true, strip) RESULT (Aeff)
      ! Returns the effective area as a function of true photon energy and incidence angles, either
      ! in total or for the thick or thin detector strips individually.
      !
      ! Input:   l10E_true               double   log (True energy of incoming photon / GeV)
      !          IncidenceAngles_true    double   [ cos (True off-axis angle of incoming photon), 
      !                                             True azimuthal angle of incoming photon (radians) ]
      !          strip                   integer  Which type of Fermi silicon strip detector to calculate
      !                                           effective area of.  Front/thin, back/thick or both/mean.
      ! Output:                          double   Effective area (m^2) 
  
      double precision, intent(IN) :: l10E_true, IncidenceAngles_true(2)
      integer, intent(IN) :: strip
      double precision :: costheta_true, phi_true
      double precision :: Aeff, DB2VAL
      external DB2VAL


      if (strip .eq. both) then

        Aeff = flatIRFs_Aeff(l10E_true, IncidenceAngles_true, thick) + &
               flatIRFs_Aeff(l10E_true, IncidenceAngles_true, thin)
        return

      elseif (strip .ne. front .and. strip .ne. back) then

        call flatUtils_crash('Invalid strip argument given to flatIRFs_Aeff')

      endif

      costheta_true = IncidenceAngles_true(1)
      phi_true = IncidenceAngles_true(2)

      Aeff = DB2VAL(l10E_true+3.d0,costheta_true,0,0,Aeff_knots_logE(:,strip),Aeff_knots_costheta(:,strip),&
             Aeff_nEnergies,Aeff_nTheta,splineOrder_logE,splineOrder_costheta,&
             Aeff_lookup(:,:,strip),little_working)
      
      if (Aeff .le. local_prec) Aeff = 0.d0

      END FUNCTION flatIRFs_Aeff


      RECURSIVE FUNCTION flatIRFs_Aeff_mean(l10E, strip) RESULT (Aeff)
      ! Returns the effective area as a function of true photon energy only, ie averaged over incidence angles in
      ! the detector.  Can be used to calculate either the total effective area or the partial areas coming from the
      ! thick or thin detector strips.
      !
      ! Input:   l10E_true               double   log (True energy of incoming photon / GeV)
      !          strip                   integer  Which type of Fermi silicon strip detector to calculate
      !                                           effective area of.  Front/thin, back/thick or both/mean.
      ! Output:                          double   Effective area (m^2) 
  
      double precision, intent(IN) :: l10E
      integer, intent(IN) :: strip
      double precision :: Aeff, Aeff_temp(1), l10Ep3(1)
      integer :: ierr

      if (strip .eq. both) then

        Aeff = flatIRFs_Aeff_mean(l10E, thick) + &
               flatIRFs_Aeff_mean(l10E, thin)
        return

      elseif (strip .ne. front .and. strip .ne. back) then

        call flatUtils_crash('Invalid strip argument given to flatIRFs_Aeff_mean')

      endif

      l10Ep3(1) = l10E + 3.d0
      call TSVAL1(Oversampled_Aeff_nEnergies, Oversampled_Aeff_logE(:,strip), mean_Aeff_lookup(:,strip), &
           mean_Aeff_zp(:,strip), mean_Aeff_sigma(:,strip), 0, 1, l10Ep3, Aeff_temp, ierr)
      if (ierr .lt. 0) call flatUtils_crash('TSVAL1(TSPACK) returned error in mean Aeff routine.')
      Aeff = Aeff_temp(1)

      if (Aeff .le. local_prec) Aeff = 0.d0

      END FUNCTION flatIRFs_Aeff_mean


      RECURSIVE FUNCTION flatIRFs_Edisp(l10E_obs, l10E_true, IncidenceAngles_true, strip) RESULT (Edisp)
      ! Returns the LAT energy dispersion as a function of true and observed photon energies and true 
      ! incidence angles, for either the thick or thin detector strips, or the LAT as a whole.
      !
      ! Input:   l10E_obs                double   log (Reconstructed energy of incoming photon / GeV)
      !          l10E_true               double   log (True energy of incoming photon / GeV)
      !          IncidenceAngles_true    double   [ cos (True off-axis angle of incoming photon), 
      !                                             True azimuthal angle of incoming photon (radians) ]
      !          strip                   integer  Which type of Fermi silicon strip detector to calculate
      !                                           effective area of.  Front/thin, back/thick or both/mean.
      ! Output:                          double   Energy dispersion (MeV^-1) 
  
      double precision, intent(IN) :: l10E_obs, l10E_true, IncidenceAngles_true(2)
      integer, intent(IN) :: strip
      double precision :: costheta_true, pars(5), Edisp

      if (strip .eq. both) then

        Edisp = 0.5d0 * flatIRFs_Edisp(l10E_obs, l10E_true, IncidenceAngles_true, thick) + &
                0.5d0 * flatIRFs_Edisp(l10E_obs, l10E_true, IncidenceAngles_true, thin)
        return

      elseif (strip .ne. front .and. strip .ne. back) then

        call flatUtils_crash('Invalid strip argument given to flatIRFs_Edisp')

      endif

      costheta_true = IncidenceAngles_true(1)
      pars(1) = flatIRFparams_xbias(costheta_true, l10E_true, strip)
      pars(2) = flatIRFparams_LS1(costheta_true, l10E_true, strip)
      pars(3) = flatIRFparams_LS2(costheta_true, l10E_true, strip)
      pars(4) = flatIRFparams_RS1(costheta_true, l10E_true, strip)
      pars(5) = flatIRFparams_RS2(costheta_true, l10E_true, strip)


      Edisp = flatIRFs_Edisp_rawfunc(l10E_obs, l10E_true, costheta_true, pars, strip) * &
              flatIRFparams_NORM(costheta_true, l10E_true, strip)

      if (Edisp .le. local_prec) Edisp = 0.d0

      END FUNCTION flatIRFs_Edisp


      DOUBLE PRECISION FUNCTION flatIRFs_Edisp_rawfunc(l10E_obs, l10E_true, costheta_true, pars, strip)
      ! Input:    
      ! Output:   
  
      double precision, intent(IN) :: l10E_obs, l10E_true, costheta_true, pars(5)
      integer, intent(IN) :: strip
      double precision :: x, scaleFactor, l10E_true_MeV, xminusxb, S1, S2
      double precision :: xbias, LS1, LS2, RS1, RS2

      xbias = pars(1)
      LS1   = pars(2)
      LS2   = pars(3)
      RS1   = pars(4)
      RS2   = pars(5)

      l10E_true_MeV = l10E_true + 3.d0
      scaleFactor = Edisp_consts(1,strip)*l10E_true_MeV*l10E_true_MeV + Edisp_consts(2,strip)*costheta_true*costheta_true + &
                    Edisp_consts(3,strip)*l10E_true_MeV + Edisp_consts(4,strip)*costheta_true + &
                    Edisp_consts(5,strip)*costheta_true*l10E_true_MeV + Edisp_consts(6,strip)
      
      x = (10.d0**(l10E_obs - l10E_true) - 1.d0) / scaleFactor
      xminusxb = x - xbias

      if (xminusxb .le. 0.d0) then  !Should be [if (x .le. 0.d0) then] when used with original bugged IRF P6 fits, cf email with R Rando April 17 2009
        S1 = LS1
        S2 = LS2
      else
        S1 = RS1
        S2 = RS2
      endif

      if (abs(xminusxb) .le. Edisp_consts(9,strip)) then

        flatIRFs_Edisp_rawfunc = exp(-0.5d0*(abs(xminusxb/S1))**Edisp_consts(7,strip))

      else

        flatIRFs_Edisp_rawfunc = exp(0.5d0*((Edisp_consts(9,strip)/S2)**Edisp_consts(8,strip) &
                 - (Edisp_consts(9,strip)/S1)**Edisp_consts(7,strip) &
                 - (abs(xminusxb/S2))**Edisp_consts(8,strip)))

      endif

      END FUNCTION flatIRFs_Edisp_rawfunc


      RECURSIVE FUNCTION flatIRFs_Edisp_mean(l10E_obs, l10E_true, strip) RESULT (Edisp)
      ! Input:    
      ! Output:   
  
      double precision, intent(IN) :: l10E_obs, l10E_true
      integer, intent(IN) :: strip
      double precision :: Edisp

      if (strip .eq. both) then

        Edisp = 0.5d0 * flatIRFs_Edisp_mean(l10E_obs, l10E_true, thick) + &
                0.5d0 * flatIRFs_Edisp_mean(l10E_obs, l10E_true, thin)
        return

      elseif (strip .ne. front .and. strip .ne. back) then

        call flatUtils_crash('Invalid strip argument given to flatIRFs_Edisp_mean')

      endif

      Edisp = surf2(l10E_obs+3.d0,l10E_true+3.d0,Oversampled_Edisp_nEnergies,Oversampled_Edisp_nEnergies,&
                    Oversampled_Edisp_logE(:,strip),Oversampled_Edisp_logE(:,strip),mean_Edisp_lookup(:,:,strip),&
                    Oversampled_Edisp_nEnergies,mean_Edisp_zp(:,strip),Edisp_splineTension)
       
      if (Edisp .le. local_prec) Edisp = 0.d0

      END FUNCTION flatIRFs_Edisp_mean


      RECURSIVE FUNCTION flatIRFs_PSF(IncidenceAngles_obs, l10E_true, IncidenceAngles_true, strip) RESULT (PSF) 
      ! Input:    
      ! Output:   

      double precision, intent(IN) :: l10E_true, IncidenceAngles_obs(2), IncidenceAngles_true(2)
      integer, intent(IN) :: strip
      double precision :: PSF, pars(5), costheta_true

      if (strip .eq. both) then

        PSF = 0.5d0 * flatIRFs_PSF(IncidenceAngles_obs, l10E_true, IncidenceAngles_true, thick) + &
              0.5d0 * flatIRFs_PSF(IncidenceAngles_obs, l10E_true, IncidenceAngles_true, thin)
        return

      elseif (strip .ne. back .and. strip .ne. front) then

        call flatUtils_crash('Invalid strip argument given to flatIRFs_PSF')

      endif

      costheta_true = IncidenceAngles_true(1)

      select case (PSFparameterisation)
      
        case (DblKingSingleSigma)

          pars(1) = flatIRFparams_sigma(costheta_true, l10E_true, strip)
          pars(2) = flatIRFparams_gammacore(costheta_true, l10E_true, strip)
          pars(3) = flatIRFparams_gammatail(costheta_true, l10E_true, strip)

        case (DblKingDblSigma)

          pars(1) = flatIRFparams_Ntail(costheta_true, l10E_true, strip)
          pars(2) = flatIRFparams_sigma(costheta_true, l10E_true, strip)
          pars(3) = flatIRFparams_sigmatail(costheta_true, l10E_true, strip)
          pars(4) = flatIRFparams_gammacore(costheta_true, l10E_true, strip)
          pars(5) = flatIRFparams_gammatail(costheta_true, l10E_true, strip)
          
      end select

      PSF = flatIRFparams_Ncore(costheta_true, l10E_true, strip) * &
            flatIRFs_PSF_rawfunc(IncidenceAngles_obs, l10E_true, IncidenceAngles_true, pars, strip)
      
      if (PSF .le. local_prec) PSF = 0.d0

      END FUNCTION flatIRFs_PSF


      DOUBLE PRECISION FUNCTION flatIRFs_PSF_rawfunc(IncidenceAngles_obs, l10E_true, IncidenceAngles_true, pars, strip)
      ! Input:    
      ! Output:   

      double precision, intent(IN) :: l10E_true, IncidenceAngles_obs(2), IncidenceAngles_true(2), pars(5)
      integer, intent(IN) :: strip
      double precision :: costheta_obs, phi_obs
      double precision :: costheta_true, phi_true
      double precision :: x, PSF, scaleFactorSq, Ntail, sigma, sigmatail, gammacore, gammatail, theta_obs, theta_true

      scaleFactorsq = PSF_consts(1,strip)*PSF_consts(1,strip) * 10.d0**(-1.6d0 * (l10E_true + 1.d0)) &
                      + PSF_consts(2,strip)*PSF_consts(2,strip)
      costheta_obs = IncidenceAngles_obs(1)
      phi_obs = IncidenceAngles_obs(2)
      costheta_true = IncidenceAngles_true(1)
      phi_true = IncidenceAngles_true(2)
      theta_obs = acos(costheta_obs)
      theta_true = acos(costheta_true)

      x = sqrt( 2.d0 * (1.d0 - cos(theta_obs-theta_true) + sin(theta_obs)*sin(theta_true) &
          * (1.d0 - cos(phi_obs-phi_true))) / scaleFactorSq)

      select case (PSFparameterisation)
      
        case (DblKingSingleSigma)

          sigma = pars(1)
          gammacore = pars(2)
          gammatail = pars(3)
      
          PSF = flatIRFs_PSFbase(x, sigma, gammacore)
          if (abs(gammatail - 1.d0) .gt. local_prec) PSF = PSF + flatIRFs_PSFbase(x, sigma, gammatail) * &
              flatIRFs_PSFbase(xbpf*sigma, sigma, gammacore)/flatIRFs_PSFbase(xbpf*sigma, sigma, gammatail)

        case (DblKingDblSigma)

          Ntail = pars(1)
          sigma = pars(2)
          sigmatail = pars(3)
          gammacore = pars(4)
          gammatail = pars(5)
      
          PSF = flatIRFs_PSFbase(x, sigma, gammacore)
          if (abs(gammatail - 1.d0) .gt. local_prec .and. abs(Ntail) .gt. local_prec) &
            PSF = PSF + Ntail*flatIRFs_PSFbase(x, sigmatail, gammatail)

        case default

          call flatUtils_crash('Unrecognised PSF parameterisation.')

      end select

      flatIRFs_PSF_rawfunc = PSF
     
      return

      END FUNCTION flatIRFs_PSF_rawfunc


      DOUBLE PRECISION FUNCTION flatIRFs_PSFbase(x, sigma, gamm)

      double precision, intent(IN) :: x, sigma, gamm

      flatIRFs_PSFbase = (1.d0 - 1.d0/gamm) * (1.d0 + 0.5d0 / gamm * x * x / sigma / sigma) ** (-1.d0 * gamm)

      return

      END FUNCTION flatIRFs_PSFbase


      RECURSIVE FUNCTION flatIRFs_PSF_mean(Direction, l10E, strip) RESULT (PSF)
      ! Input:    
      ! Output:   

      double precision, intent(IN) :: Direction(2), l10E
      integer, intent(IN) :: strip
      double precision :: theta, unscaledx, PSF

      if (strip .eq. both) then

        PSF = 0.5d0 * flatIRFs_PSF_mean(Direction, l10E, thick) + &
              0.5d0 * flatIRFs_PSF_mean(Direction, l10E, thin)
        return

      elseif (strip .ne. back .and. strip .ne. front) then

        call flatUtils_crash('Invalid strip argument given to flatIRFs_PSF_mean')

      endif

      theta = sqrt(Direction(1)*Direction(1) + Direction(2)*Direction(2))
      unscaledx = sqrt( 2.d0 * (1.d0 - cos(theta)))
      
      PSF = surf2(l10E+3.d0,unscaledx,Oversampled_PSF_nEnergies,Oversampled_PSF_nEnergies,&
                  Oversampled_PSF_logE(:,strip),Oversampled_PSF_unscaledx(:,strip),mean_PSF_lookup(:,:,strip),&
                  Oversampled_PSF_nEnergies,mean_PSF_zp(:,strip),PSF_splineTension)

      if (PSF .le. local_prec) PSF = 0.d0

      END FUNCTION flatIRFs_PSF_mean


      END MODULE flatIRFs
