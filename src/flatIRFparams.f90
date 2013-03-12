! Fermi LAT Instrumental Response Function parameter lookup routines 
!
! These return the values of the paramters in the IRFs calculated in
! flatIRFs, using lookup tables produced in flatIRFini.
!
! In principle they could all easily be combined into a single routine,
! with the relevant arrays all given an extra dimension indexed by a
! constant integer keyword (just like back/front in the detector index).
! I have instead set them up 'the long way' as separate routines and arrays
! so that when phi-dependent IRFs become available in pass7, it will be easier
! to accomodate the situation where only some (not all) parameters are
! phi-dependent.
!
! Pat Scott, Feb 2009; pat@fysik.su.se
! Modified:	Pat Scott, Oct 2011	Added support for Pass7-type PSFs
!--------------------------------------------------------------------------------


      MODULE flatIRFparams

      use flatCommon
      use flatUtils

      implicit none

      double precision DB2VAL, surf2
      external DB2VAL, surf2

      contains


      SUBROUTINE flatIRFparams_checkBounds(value_in, value_out, array_len, array, strip, string1, string2)

      integer, intent(IN) :: strip, array_len
      double precision, intent(IN) :: value_in
      double precision, intent(IN) :: array(array_len,2)
      character (len=*), intent(IN) :: string1, string2
      double precision, intent(OUT) :: value_out

      if (value_in .lt. array(1,strip) - local_prec) then
        write(*,*) string1, ': ', value_in, '   array bottom: ', array(1,strip) - local_prec 
        write(*,*) string1, ' below permitted minimum in ', string2
        call flatUtils_crash('Exiting...')
      else if (value_in .gt. array(array_len,strip) + local_prec) then
        write(*,*) string1, ': ', value_in, '   array top: ', array(array_len,strip) + local_prec 
        write(*,*) string1,'  above permitted maximum in ', string2
        call flatUtils_crash('Exiting...')
      else
        value_out = value_in
      endif

      END SUBROUTINE

      
      DOUBLE PRECISION FUNCTION flatIRFparams_NORM(costheta_in, logE_in, strip)
      ! Input:    
      ! Output:   
   
      double precision, intent(IN) :: costheta_in, logE_in
      integer, intent(IN) :: strip
      double precision :: costheta, logE
      double precision :: NORM
      
      call flatIRFparams_checkBounds(costheta_in, costheta, Edisp_nTheta, Edisp_costheta, strip, 'Costheta', 'flatIRFparams_NORM')
      call flatIRFparams_checkBounds(logE_in+3.d0, logE, Edisp_nEnergies, Edisp_logE, strip, 'Energy', 'flatIRFparams_NORM')

      NORM = surf2(logE, costheta, Edisp_nEnergies,Edisp_nTheta,Edisp_logE(:,strip),Edisp_costheta(:,strip),&
                   NORM_lookup(:,:,strip),Edisp_nEnergies,NORM_zp(:,strip),NORM_splineTension)       
     
      if (NORM > epsilon(NORM)) then
        flatIRFparams_NORM = NORM
      else
        flatIRFparams_NORM = 0.d0
      endif

      return

      END FUNCTION flatIRFparams_NORM


      DOUBLE PRECISION FUNCTION flatIRFparams_xbias(costheta_in, logE_in, strip)
      ! Input:    
      ! Output:   

      double precision, intent(IN) :: costheta_in, logE_in
      integer, intent(IN) :: strip
      double precision :: costheta, logE
      double precision :: xbias
      

      call flatIRFparams_checkBounds(costheta_in, costheta, Edisp_nTheta, Edisp_costheta, strip, 'Costheta', 'flatIRFparams_xbias')
      call flatIRFparams_checkBounds(logE_in+3.d0, logE, Edisp_nEnergies, Edisp_logE, strip, 'Energy', 'flatIRFparams_xbias')

      xbias = DB2VAL(logE,costheta,0,0,xbias_knots_logE(:,strip),xbias_knots_costheta(:,strip),&
             Edisp_nEnergies,Edisp_nTheta,splineOrder_logE,splineOrder_costheta,&
             xbias_lookup(:,:,strip),little_working)
      
      flatIRFparams_xbias = xbias
      
      return

      END FUNCTION flatIRFparams_xbias


      DOUBLE PRECISION FUNCTION flatIRFparams_LS1(costheta_in, logE_in, strip)
      ! Input:    
      ! Output:   
   
      double precision, intent(IN) :: costheta_in, logE_in
      integer, intent(IN) :: strip
      double precision :: costheta, logE
      double precision :: LS1
      
      call flatIRFparams_checkBounds(costheta_in, costheta, Edisp_nTheta, Edisp_costheta, strip, 'Costheta', 'flatIRFparams_LS1')
      call flatIRFparams_checkBounds(logE_in+3.d0, logE, Edisp_nEnergies, Edisp_logE, strip, 'Energy', 'flatIRFparams_LS1')

      LS1 = DB2VAL(logE,costheta,0,0,LS1_knots_logE(:,strip),LS1_knots_costheta(:,strip),&
             Edisp_nEnergies,Edisp_nTheta,splineOrder_logE,splineOrder_costheta,&
             LS1_lookup(:,:,strip),little_working)
      
      if (LS1 > epsilon(LS1)) then
        flatIRFparams_LS1 = LS1
      else
        flatIRFparams_LS1 = epsilon(LS1)
      endif
      
      return

      END FUNCTION flatIRFparams_LS1


      DOUBLE PRECISION FUNCTION flatIRFparams_RS1(costheta_in, logE_in, strip)
      ! Input:    
      ! Output:   
   
      double precision, intent(IN) :: costheta_in, logE_in
      integer, intent(IN) :: strip
      double precision :: costheta, logE
      double precision :: RS1
      
      call flatIRFparams_checkBounds(costheta_in, costheta, Edisp_nTheta, Edisp_costheta, strip, 'Costheta', 'flatIRFparams_RS1')
      call flatIRFparams_checkBounds(logE_in+3.d0, logE, Edisp_nEnergies, Edisp_logE, strip, 'Energy', 'flatIRFparams_RS1')

      RS1 = DB2VAL(logE,costheta,0,0,RS1_knots_logE(:,strip),RS1_knots_costheta(:,strip),&
             Edisp_nEnergies,Edisp_nTheta,splineOrder_logE,splineOrder_costheta,&
             RS1_lookup(:,:,strip),little_working)
      
      if (RS1 > epsilon(RS1)) then
        flatIRFparams_RS1 = RS1
      else
        flatIRFparams_RS1 = epsilon(RS1)
      endif
      
      return

      END FUNCTION flatIRFparams_RS1


      DOUBLE PRECISION FUNCTION flatIRFparams_LS2(costheta_in, logE_in, strip)
      ! Input:    
      ! Output:   

      double precision, intent(IN) :: costheta_in, logE_in
      integer, intent(IN) :: strip
      double precision :: costheta, logE
      double precision :: LS2
      
      call flatIRFparams_checkBounds(costheta_in, costheta, Edisp_nTheta, Edisp_costheta, strip, 'Costheta', 'flatIRFparams_LS2')
      call flatIRFparams_checkBounds(logE_in+3.d0, logE, Edisp_nEnergies, Edisp_logE, strip, 'Energy', 'flatIRFparams_LS2')

      LS2 = DB2VAL(logE,costheta,0,0,LS2_knots_logE(:,strip),LS2_knots_costheta(:,strip),&
             Edisp_nEnergies,Edisp_nTheta,splineOrder_logE,splineOrder_costheta,&
             LS2_lookup(:,:,strip),little_working)
      
      if (LS2 > epsilon(LS2)) then
        flatIRFparams_LS2 = LS2
      else
        flatIRFparams_LS2 = epsilon(LS2)
      endif
      
      return

      END FUNCTION flatIRFparams_LS2


      DOUBLE PRECISION FUNCTION flatIRFparams_RS2(costheta_in, logE_in, strip)
      ! Input:    
      ! Output:   

      double precision, intent(IN) :: costheta_in, logE_in
      integer, intent(IN) :: strip
      double precision :: costheta, logE
      double precision :: RS2
      
      call flatIRFparams_checkBounds(costheta_in, costheta, Edisp_nTheta, Edisp_costheta, strip, 'Costheta', 'flatIRFparams_RS2')
      call flatIRFparams_checkBounds(logE_in+3.d0, logE, Edisp_nEnergies, Edisp_logE, strip, 'Energy', 'flatIRFparams_RS2')

      RS2 = DB2VAL(logE,costheta,0,0,RS2_knots_logE(:,strip),RS2_knots_costheta(:,strip),&
             Edisp_nEnergies,Edisp_nTheta,splineOrder_logE,splineOrder_costheta,&
             RS2_lookup(:,:,strip),little_working)
      
      if (RS2 > epsilon(RS2)) then
        flatIRFparams_RS2 = RS2
      else
        flatIRFparams_RS2 = epsilon(RS2)
      endif
      
      return

      END FUNCTION flatIRFparams_RS2


      DOUBLE PRECISION FUNCTION flatIRFparams_Ncore(costheta_in, logE_in, strip)
      ! Input:    
      ! Output:   
   
      double precision, intent(IN) :: costheta_in, logE_in
      integer, intent(IN) :: strip
      double precision :: costheta, logE
      double precision :: Ncore
      
      call flatIRFparams_checkBounds(costheta_in, costheta, PSF_nTheta, PSF_costheta, strip, 'Costheta', 'flatIRFparams_Ncore')
      call flatIRFparams_checkBounds(logE_in+3.d0, logE, PSF_nEnergies, PSF_logE, strip, 'Energy', 'flatIRFparams_Ncore')

      Ncore = surf2(logE, costheta, PSF_nEnergies,PSF_nTheta,PSF_logE(:,strip),PSF_costheta(:,strip),&
                    Ncore_lookup(:,:,strip),PSF_nEnergies,Ncore_zp(:,strip),Ncore_splineTension)       

      if (Ncore > epsilon(Ncore)) then
        flatIRFparams_Ncore = Ncore
      else
        flatIRFparams_Ncore = 0.d0
      endif
      
      return

      END FUNCTION flatIRFparams_Ncore


      DOUBLE PRECISION FUNCTION flatIRFparams_Ntail(costheta_in, logE_in, strip)
      ! Input:    
      ! Output:   
   
      double precision, intent(IN) :: costheta_in, logE_in
      integer, intent(IN) :: strip
      double precision :: costheta, logE
      double precision :: Ntail
      
      call flatIRFparams_checkBounds(costheta_in, costheta, PSF_nTheta, PSF_costheta, strip, 'Costheta', 'flatIRFparams_Ntail')
      call flatIRFparams_checkBounds(logE_in+3.d0, logE, PSF_nEnergies, PSF_logE, strip, 'Energy', 'flatIRFparams_Ntail')

      Ntail = surf2(logE, costheta, PSF_nEnergies,PSF_nTheta,PSF_logE(:,strip),PSF_costheta(:,strip),&
                    Ntail_lookup(:,:,strip),PSF_nEnergies,Ntail_zp(:,strip),Ntail_splineTension)       

      if (Ntail > epsilon(Ntail)) then
        flatIRFparams_Ntail = Ntail
      else
        flatIRFparams_Ntail = 0.d0
      endif
      
      return

      END FUNCTION flatIRFparams_Ntail


      DOUBLE PRECISION FUNCTION flatIRFparams_sigma(costheta_in, logE_in, strip)
      ! Input:    
      ! Output:   
   
      double precision, intent(IN) :: costheta_in, logE_in
      integer, intent(IN) :: strip
      double precision :: costheta, logE
      double precision :: sigma

      call flatIRFparams_checkBounds(costheta_in, costheta, PSF_nTheta, PSF_costheta, strip, 'Costheta', 'flatIRFparams_sigma')
      call flatIRFparams_checkBounds(logE_in+3.d0, logE, PSF_nEnergies, PSF_logE, strip, 'Energy', 'flatIRFparams_sigma')

      sigma = DB2VAL(logE,costheta,0,0,sigma_knots_logE(:,strip),sigma_knots_costheta(:,strip),&
             PSF_nEnergies,PSF_nTheta,splineOrder_logE,splineOrder_costheta,&
             sigma_lookup(:,:,strip),little_working)
      
      if (sigma > epsilon(sigma)) then
        flatIRFparams_sigma = sigma
      else
        flatIRFparams_sigma = 0.d0
      endif
      
      return

      END FUNCTION flatIRFparams_sigma


      DOUBLE PRECISION FUNCTION flatIRFparams_sigmatail(costheta_in, logE_in, strip)
      ! Input:    
      ! Output:   
   
      double precision, intent(IN) :: costheta_in, logE_in
      integer, intent(IN) :: strip
      double precision :: costheta, logE
      double precision :: sigmatail

      call flatIRFparams_checkBounds(costheta_in, costheta, PSF_nTheta, PSF_costheta, strip, 'Costheta', 'flatIRFparams_sigmatail')
      call flatIRFparams_checkBounds(logE_in+3.d0, logE, PSF_nEnergies, PSF_logE, strip, 'Energy', 'flatIRFparams_sigmatail')

      sigmatail = DB2VAL(logE,costheta,0,0,sigmatail_knots_logE(:,strip),&
             sigmatail_knots_costheta(:,strip),PSF_nEnergies,PSF_nTheta,splineOrder_logE,&
             splineOrder_costheta,sigmatail_lookup(:,:,strip),little_working)
      
      if (sigmatail > epsilon(sigmatail)) then
        flatIRFparams_sigmatail = sigmatail
      else
        flatIRFparams_sigmatail = 0.d0
      endif
      
      return

      END FUNCTION flatIRFparams_sigmatail


      DOUBLE PRECISION FUNCTION flatIRFparams_gammacore(costheta_in, logE_in, strip)
      ! Input:    
      ! Output:   
   
      double precision, intent(IN) :: costheta_in, logE_in
      integer, intent(IN) :: strip
      double precision :: costheta, logE
      double precision :: gammacore
      
      call flatIRFparams_checkBounds(costheta_in, costheta, PSF_nTheta, PSF_costheta, strip, 'Costheta', 'flatIRFparams_gammacore')
      call flatIRFparams_checkBounds(logE_in+3.d0, logE, PSF_nEnergies, PSF_logE, strip, 'Energy', 'flatIRFparams_gammacore')

      gammacore = DB2VAL(logE,costheta,0,0,gammacore_knots_logE(:,strip),gammacore_knots_costheta(:,strip),&
                  PSF_nEnergies,PSF_nTheta,splineOrder_logE,splineOrder_costheta,&
                  gammacore_lookup(:,:,strip),little_working)
      
      if (gammacore > epsilon(gammacore)) then
        flatIRFparams_gammacore = gammacore
      else
        flatIRFparams_gammacore = 0.d0
      endif
      
      return

      END FUNCTION flatIRFparams_gammacore


      DOUBLE PRECISION FUNCTION flatIRFparams_gammatail(costheta_in, logE_in, strip)
      ! Input:    
      ! Output:   
   
      double precision, intent(IN) :: costheta_in, logE_in
      integer, intent(IN) :: strip
      double precision :: costheta, logE
      double precision :: gammatail
      
      call flatIRFparams_checkBounds(costheta_in, costheta, PSF_nTheta, PSF_costheta, strip, 'Costheta', 'flatIRFparams_gammatail')
      call flatIRFparams_checkBounds(logE_in+3.d0, logE, PSF_nEnergies, PSF_logE, strip, 'Energy', 'flatIRFparams_gammatail')

      gammatail = DB2VAL(logE,costheta,0,0,gammatail_knots_logE(:,strip),gammatail_knots_costheta(:,strip),&
                  PSF_nEnergies,PSF_nTheta,splineOrder_logE,splineOrder_costheta,&
                  gammatail_lookup(:,:,strip),little_working)
      
      if (gammatail > epsilon(gammatail)) then
        flatIRFparams_gammatail = gammatail
      else
        flatIRFparams_gammatail = 0.d0
      endif
      
      return

      END FUNCTION flatIRFparams_gammatail


      END MODULE flatIRFparams
