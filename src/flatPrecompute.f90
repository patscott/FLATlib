! Fermi LAT Instrumental Response Function initialisation routines
!
! These read in the IRF parameters from FITS files and  set up
! lookup tables for the IRF parameter routines in flatIRFparams.
!
! Pat Scott, Feb 2009; pat@fysik.su.se
!--------------------------------------------------------------------------------


      MODULE flatPrecompute

      use flatCommon
      use flatUtils
      use flatIRFs
      use Precision_Model
      use CUI

      implicit none

      double precision :: unscaledx, l10E_obs, l10E_true
      integer :: strip

      contains


      SUBROUTINE flatPrecompute_PSF(IRF)

      character(len=*), intent(IN) :: IRF     
      integer :: i,j, SRgType, IFAIL, ierr, fileunit
      double precision :: SVertices1(1,2), SVertices2(2,3), SValue, SAbsErr
      double precision, allocatable :: mean_PSF_edge_derivatives(:), temp(:), mean_PSF_norm(:,:)
      character (len=strlen) :: outfilename
      
      Oversampled_PSF_nEnergies = (PSF_nEnergies-1)*PSF_OversamplingFactor + 1

      allocate(Oversampled_PSF_logE(Oversampled_PSF_nEnergies,2))
      allocate(Oversampled_PSF_unscaledx(Oversampled_PSF_nEnergies,2))
      allocate(mean_PSF_lookup(Oversampled_PSF_nEnergies,Oversampled_PSF_nEnergies,2))
      allocate(mean_PSF_zp(3*Oversampled_PSF_nEnergies*Oversampled_PSF_nEnergies,2))
      allocate(mean_PSF_edge_derivatives(Oversampled_PSF_nEnergies))
      allocate(temp(3*Oversampled_PSF_nEnergies))
      allocate(mean_PSF_norm(Oversampled_PSF_nEnergies,2))
      
      mean_PSF_edge_derivatives = 0.d0

      SVertices1(1,:) = (/acos(Fermi_costhetamax), acos(Fermi_costhetamin)/)
      SVertices2(:,1) = (/-1.d0*PSF_maxAngle, -1.d0*PSF_maxAngle/)
      SVertices2(:,2) = (/-1.d0*PSF_maxAngle, PSF_maxAngle/)
      SVertices2(:,3) = (/PSF_maxAngle, -1.d0*PSF_maxAngle/)
      SRgType = HyperQuad
      
      do strip=front,back,back-front

        Oversampled_PSF_logE(Oversampled_PSF_nEnergies,strip) = PSF_logE(PSF_nEnergies,strip)
        forall(i=1:PSF_nEnergies-1,j=1:PSF_OversamplingFactor) &
               Oversampled_PSF_logE((i-1)*PSF_OversamplingFactor+j,strip) = dble(j-1)*(PSF_logE(i+1,strip) - &
               PSF_logE(i,strip))/dble(PSF_OversamplingFactor) + PSF_logE(i,strip)
        forall(i=1:Oversampled_PSF_nEnergies) Oversampled_PSF_unscaledx(i,strip) = sqrt(2.d0)*&
              (dble(i-1)/dble(Oversampled_PSF_nEnergies-1))**2

        !Populate lookup table
        do i=1,Oversampled_PSF_nEnergies
          l10E_true = Oversampled_PSF_logE(i,strip) - 3.d0
          do j=1,Oversampled_PSF_nEnergies
            unscaledx = Oversampled_PSF_unscaledx(j,strip)
            IFAIL = 0
            call CUBATR(1,flatPrecompute_PSF_Integrand,SVertices1,SRgType,SValue,SAbsErr,IFAIL,MaxPts=50000000,&
                        EpsAbs=integratorAbsErr,EpsRel=1.d-3, Job=11)
            mean_PSF_lookup(i,j,strip) = SValue / (SVertices1(1,2) - SVertices1(1,1))
          enddo
          call CUBATR()
        enddo

        !Initialise the interpolator
        ierr = 0
        call surf1(Oversampled_PSF_nEnergies,Oversampled_PSF_nEnergies,Oversampled_PSF_logE(:,strip),&
                   Oversampled_PSF_unscaledx(:,strip),mean_PSF_lookup(:,:,strip),Oversampled_PSF_nEnergies,&
                   mean_PSF_edge_derivatives,mean_PSF_edge_derivatives,mean_PSF_edge_derivatives,&
                   mean_PSF_edge_derivatives,0.d0,0.d0,0.d0,0.d0,255,mean_PSF_zp(:,strip),temp,PSF_splineTension,ierr)
        if (ierr .ne. 0) call flatUtils_crash('surf1(FITPACK) returned error in mean PSF precomputation routine.')

        !Work out the new normalisation factors
        do i=1,Oversampled_PSF_nEnergies
          l10E_true = Oversampled_PSF_logE(i,strip) - 3.d0
          IFAIL = 0
          call CUBATR(2,flatPrecompute_PSF_normIntegrand,SVertices2,SRgType,SValue,SAbsErr,IFAIL,MaxPts=2000000000,&
                      EpsAbs=integratorAbsErr,EpsRel=1.d-4,Job=11)
          mean_PSF_norm(i,strip) = 1.d0/SValue
          call CUBATR()
        enddo

      enddo

      forall(i=1:Oversampled_PSF_nEnergies, strip=front:back:back-front) mean_PSF_lookup(i,:,strip) = &
        mean_PSF_lookup(i,:,strip) * mean_PSF_norm(i,strip)
      
      ierr=0
      call ftgiou(fileunit,ierr)
      if (ierr .ne. 0) call flatUtils_crash('No free file handle available for flatPrecompute_PSF!')

      outfilename = trim(flatlib_datadir) // trim(IRF) // '/psf_' // trim(IRF) // '_mean.dat'
      open(fileunit, STATUS='REPLACE', FILE=outfilename, IOSTAT=ierr, ACTION='WRITE')
      if (ierr .ne. 0) call flatUtils_crash('Could not open write file for flatPrecompute_PSF!')

      write(fileunit, '(I10)') Oversampled_PSF_nEnergies
      write(fileunit, '(ES13.6)') Oversampled_PSF_logE
      write(fileunit, '(ES13.6)') Oversampled_PSF_unscaledx
      write(fileunit, '(ES13.6)') mean_PSF_lookup

      close(fileunit)

      deallocate(mean_PSF_edge_derivatives, temp, mean_PSF_norm)
      deallocate(mean_PSF_lookup, Oversampled_PSF_logE, Oversampled_PSF_unscaledx, mean_PSF_zp)

      END SUBROUTINE flatPrecompute_PSF


      SUBROUTINE flatPrecompute_Edisp(IRF)

      character(len=*), intent(IN) :: IRF     
      integer :: i,j, SRgType, IFAIL, ierr, fileunit
      double precision :: SVertices1(1,2), SVertices2(1,2), SValue, SAbsErr
      double precision, allocatable :: mean_Edisp_edge_derivatives(:), temp(:), mean_Edisp_norm(:,:)
      character (len=strlen) :: outfilename
      
      Oversampled_Edisp_nEnergies = (Edisp_nEnergies-1)*Edisp_OversamplingFactor + 1

      allocate(Oversampled_Edisp_logE(Oversampled_Edisp_nEnergies,2))
      allocate(mean_Edisp_lookup(Oversampled_Edisp_nEnergies, Oversampled_Edisp_nEnergies,2))
      allocate(mean_Edisp_zp(3*Oversampled_Edisp_nEnergies*Oversampled_Edisp_nEnergies,2))
      allocate(mean_Edisp_edge_derivatives(Oversampled_Edisp_nEnergies))
      allocate(temp(3*Oversampled_Edisp_nEnergies))
      allocate(mean_Edisp_norm(Oversampled_Edisp_nEnergies,2))
     
      mean_Edisp_edge_derivatives = 0.d0

      SVertices1(1,:) = (/acos(Fermi_costhetamax), acos(Fermi_costhetamin)/)
      SVertices2(1,:) = (/log10(Fermi_Emin)+3.d0, log10(Fermi_Emax)+3.d0/)
      SRgType = HyperQuad
      
      do strip=front,back,back-front

        Oversampled_Edisp_logE(Oversampled_Edisp_nEnergies,strip) = Edisp_logE(Edisp_nEnergies,strip)
        forall(i=1:Edisp_nEnergies-1,j=1:Edisp_OversamplingFactor) &
               Oversampled_Edisp_logE((i-1)*Edisp_OversamplingFactor+j,strip) = dble(j-1)*(Edisp_logE(i+1,strip) - &
               Edisp_logE(i,strip))/dble(Edisp_OversamplingFactor) + Edisp_logE(i,strip)

        !Populate lookup table
        do i=1,Oversampled_Edisp_nEnergies
          do j=1,Oversampled_Edisp_nEnergies
            l10E_obs = Oversampled_Edisp_logE(i,strip) - 3.d0
            l10E_true = Oversampled_Edisp_logE(j,strip) - 3.d0
            IFAIL = 0
            call CUBATR(1,flatPrecompute_Edisp_Integrand,SVertices1,SRgType,SValue,SAbsErr,IFAIL,MaxPts=5000000,&
                        EpsAbs=integratorAbsErr,EpsRel=integratorRelErr,Job=11)
            mean_Edisp_lookup(i,j,strip) = SValue / (SVertices1(1,2) - SVertices1(1,1))
          enddo
          call CUBATR()
        enddo

        !Initialise the interpolator
        call surf1(Oversampled_Edisp_nEnergies,Oversampled_Edisp_nEnergies,Oversampled_Edisp_logE(:,strip),&
                   Oversampled_Edisp_logE(:,strip),mean_Edisp_lookup(:,:,strip),Oversampled_Edisp_nEnergies,&
                   mean_Edisp_edge_derivatives,mean_Edisp_edge_derivatives,mean_Edisp_edge_derivatives,&
                   mean_Edisp_edge_derivatives,0.d0,0.d0,0.d0,0.d0,255,mean_Edisp_zp(:,strip),temp,Edisp_splineTension,ierr)
        if (ierr .ne. 0) call flatUtils_crash('surf1(FITPACK) returned error in mean Edisp precomputation routine.')

        !Work out the new normalisation factors
        do i=1,Oversampled_Edisp_nEnergies
          l10E_true = Oversampled_Edisp_logE(i,strip) - 3.d0
          IFAIL = 0
          call CUBATR(1,flatPrecompute_Edisp_normIntegrand,SVertices2,SRgType,SValue,SAbsErr,IFAIL,MaxPts=500000000,&
                        EpsAbs=integratorAbsErr,EpsRel=1.d-5,Job=11)
          mean_Edisp_norm(i,strip) = 1.d0/SValue
          call CUBATR()
        enddo

      enddo

      forall(i=1:Oversampled_Edisp_nEnergies, strip=front:back:back-front) mean_Edisp_lookup(i,:,strip) = &
        mean_Edisp_lookup(i,:,strip) * mean_Edisp_norm(i,strip)

      ierr=0
      call ftgiou(fileunit,ierr)
      if (ierr .ne. 0) call flatUtils_crash('No free file handle available for flatPrecompute_Edisp!')

      outfilename = trim(flatlib_datadir) // trim(IRF) // '/edisp_' // trim(IRF) // '_mean.dat'
      open(fileunit, STATUS='REPLACE', FILE=outfilename, IOSTAT=ierr, ACTION='WRITE')
      if (ierr .ne. 0) call flatUtils_crash('Could not open write file for flatPrecompute_Edisp!')

      write(fileunit, '(I10)') Oversampled_Edisp_nEnergies
      write(fileunit, '(ES13.6)') Oversampled_Edisp_logE
      write(fileunit, '(ES13.6)') mean_Edisp_lookup

      close(fileunit)

      deallocate(mean_Edisp_edge_derivatives, temp, mean_Edisp_norm)
      deallocate(mean_Edisp_lookup, Oversampled_Edisp_logE, mean_Edisp_zp)   
      
      END SUBROUTINE flatPrecompute_Edisp


      SUBROUTINE flatPrecompute_Aeff(IRF)

      character(len=*), intent(IN) :: IRF     
      integer :: i,j, SRgType, IFAIL, ierr, fileunit
      double precision :: SVertices(1,2), SValue, SAbsErr
      character (len=strlen) :: outfilename
      
      Oversampled_Aeff_nEnergies = (Aeff_nEnergies-1)*Aeff_OversamplingFactor + 1

      allocate(Oversampled_Aeff_logE(Oversampled_Aeff_nEnergies,2))
      allocate(mean_Aeff_lookup(Oversampled_Aeff_nEnergies, 2))
      
      SVertices(1,:) = (/acos(Fermi_costhetamax), acos(Fermi_costhetamin)/)
      SRgType = HyperQuad
      
      do strip=front,back,back-front

        Oversampled_Aeff_logE(Oversampled_Aeff_nEnergies,strip) = Aeff_logE(Aeff_nEnergies,strip)
        forall(i=1:Aeff_nEnergies-1,j=1:Aeff_OversamplingFactor) &
               Oversampled_Aeff_logE((i-1)*Aeff_OversamplingFactor+j,strip) = dble(j-1)*(Aeff_logE(i+1,strip) - Aeff_logE(i,strip))/&
               Aeff_OversamplingFactor + Aeff_logE(i,strip)

        !Populate lookup table
        do i=1,Oversampled_Aeff_nEnergies
          l10E_true = Oversampled_Aeff_logE(i,strip) - 3.d0
          IFAIL = 0
          call CUBATR(1,flatPrecompute_Aeff_Integrand,SVertices,SRgType,SValue,SAbsErr,IFAIL,MaxPts=5000000,&
                      EpsAbs=integratorAbsErr,EpsRel=integratorRelErr,Job=11)
          mean_Aeff_lookup(i,strip) = SValue / radperdeg
        enddo

      enddo

      ierr=0
      call ftgiou(fileunit,ierr)
      if (ierr .ne. 0) call flatUtils_crash('No free file handle available for flatPrecompute_Aeff!')

      outfilename = trim(flatlib_datadir) // trim(IRF) // '/aeff_' // trim(IRF) // '_mean.dat'
      open(fileunit, STATUS='REPLACE', FILE=outfilename, IOSTAT=ierr, ACTION='WRITE')
      if (ierr .ne. 0) call flatUtils_crash('Could not open write file for flatPrecompute_Aeff!')

      write(fileunit, '(I10)') Oversampled_Aeff_nEnergies
      write(fileunit, '(ES13.6)') Oversampled_Aeff_logE
      write(fileunit, '(ES13.6)') mean_Aeff_lookup

      close(fileunit)    

      deallocate(mean_Aeff_lookup, Oversampled_Aeff_logE)

      END SUBROUTINE flatPrecompute_Aeff


      FUNCTION flatPrecompute_Edisp_Integrand(NumFun,X) RESULT(Value)
      ! Input:    
      ! Output:   
        
      integer, intent(IN) :: NumFun
      double precision, intent(IN) :: X(:)

      double precision :: Value(NumFun)

      Value(1) = flatIRFs_Edisp(l10E_obs, l10E_true, (/cos(X(1)), 0.d0/), strip)

      END FUNCTION flatPrecompute_Edisp_Integrand


      FUNCTION flatPrecompute_Edisp_normIntegrand(NumFun,X) RESULT(Value)
      ! Input:    
      ! Output:   
        
      integer, intent(IN) :: NumFun
      double precision, intent(IN) :: X(:)

      double precision :: Value(NumFun)
      
      Value(1) = flatIRFs_Edisp_mean(X(1)-3.d0, l10E_true, strip)*10.d0**(X(1))*log(10.d0)

      END FUNCTION flatPrecompute_Edisp_normIntegrand


      FUNCTION flatPrecompute_PSF_Integrand(NumFun,X) RESULT(Value)
      ! Input:    
      ! Output:   
        
      use flatIRFParams

      integer, intent(IN) :: NumFun
      double precision, intent(IN) :: X(:)

      double precision :: Value(NumFun), costheta_true, fullx, scaleFactorSq, sigma, gammacore, gammatail

      scaleFactorSq = PSF_consts(1,strip)*PSF_consts(1,strip) * 10.d0**(-1.6d0 * (l10E_true + 1.d0)) &
                      + PSF_consts(2,strip)*PSF_consts(2,strip)
      fullx = unscaledx / sqrt(scaleFactorSq)
      costheta_true = cos(X(1))

      sigma = flatIRFparams_sigma(costheta_true, l10E_true, strip)
      gammacore = flatIRFparams_gammacore(costheta_true, l10E_true, strip)
      gammatail = flatIRFparams_gammatail(costheta_true, l10E_true, strip)
      
      Value(1) = flatIRFs_PSFbase(fullx, sigma, gammacore)
      if (abs(gammatail - 1.d0) .gt. local_prec) Value(1) = Value(1) + flatIRFs_PSFbase(fullx, sigma, gammatail) * &
            flatIRFs_PSFbase(xbpf*sigma, sigma, gammacore)/flatIRFs_PSFbase(xbpf*sigma, sigma, gammatail) 
      Value(1) = Value(1) * flatIRFparams_Ncore(costheta_true, l10E_true, strip)

      END FUNCTION flatPrecompute_PSF_Integrand


      FUNCTION flatPrecompute_PSF_normIntegrand(NumFun,X) RESULT(Value)
      ! Input:    
      ! Output:   
        
      integer, intent(IN) :: NumFun
      double precision, intent(IN) :: X(:)

      double precision :: Value(NumFun), Direction(2)
      
      Direction = (/X(1), X(2)/)
      Value(1) = flatIRFs_PSF_mean(Direction, l10E_true, strip)

      END FUNCTION flatPrecompute_PSF_normIntegrand


      FUNCTION flatPrecompute_Aeff_Integrand(NumFun,X) RESULT(Value)
      ! Input:    
      ! Output:   
        
      integer, intent(IN) :: NumFun
      double precision, intent(IN) :: X(:)

      double precision :: Value(NumFun), IncidenceAngles_true(2)

      IncidenceAngles_true = (/cos(X(1)), 0.d0/)
      Value(1) = flatIRFs_Aeff(l10E_true, IncidenceAngles_true, strip)

      END FUNCTION flatPrecompute_Aeff_Integrand


      END MODULE flatPrecompute
