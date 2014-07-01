! Fermi LAT Instrumental Response Function initialisation routines
!
! These read in the IRF parameters from FITS files and  set up
! lookup tables for the IRF parameter routines in flatIRFparams.
!
! Pat Scott, Feb 2009; pat@fysik.su.se
! Modified:	Pat Scott, Oct 2011	Added support for Pass7-type PSFs
!--------------------------------------------------------------------------------


      MODULE flatIRFini

      use flatFFTW      
      use flatPrecompute

      implicit none

      double precision, parameter :: nulle=0.

      double precision :: Edisp_pars(5), PSF_pars(5)
      double precision :: costheta_true_local, phi_true_local
      integer :: fstatus, unit, blocksize, hdutype=0, irow, jrow
      integer :: felem, datacode, datawidth
      character (len=strlen) :: record,comment
      character (len=7) :: detnam(2)
      logical anynull

      contains


      SUBROUTINE flatIRFini_initIRFs(IRF)

      character(len=*), intent(IN) :: IRF

      call flatIRFini_Edisp(IRF)
      call flatIRFini_PSF(IRF)
      call flatIRFini_Aeff(IRF)

      END SUBROUTINE
      

      SUBROUTINE flatIRFini_initFastIRFs(IRF,balkIfMissing)

      character(len=*), intent(IN) :: IRF
      logical, intent(IN) :: balkIfMissing

      call flatIRFini_Edisp_mean(IRF,balkIfMissing)
      call flatIRFini_PSF_mean(IRF,balkIfMissing)
      call flatIRFini_Aeff_mean(IRF,balkIfMissing)

      END SUBROUTINE


      SUBROUTINE flatIRFini_cleanIRFs

      call flatIRFini_Edisp_clean
      call flatIRFini_PSF_clean
      call flatIRFini_Aeff_clean

      END SUBROUTINE


      SUBROUTINE flatIRFini_cleanFastIRFs

      call flatIRFini_Edisp_mean_clean
      call flatIRFini_PSF_mean_clean
      call flatIRFini_Aeff_mean_clean
      
      END SUBROUTINE


      SUBROUTINE flatIRFini_Edisp_clean

      deallocate(Edisp_costheta, NORM_lookup, xbias_lookup, Edisp_logE)
      deallocate(LS1_lookup, RS1_lookup, LS2_lookup, RS2_lookup)
      deallocate(NORM_zp)
      deallocate(xbias_knots_logE, xbias_knots_costheta)
      deallocate(LS1_knots_logE, LS1_knots_costheta, RS1_knots_logE)
      deallocate(RS1_knots_costheta, LS2_knots_logE, LS2_knots_costheta)
      deallocate(RS2_knots_logE, RS2_knots_costheta)

      END SUBROUTINE


      SUBROUTINE flatIRFini_PSF_clean

      deallocate(PSF_costheta, Ncore_lookup, sigma_lookup, gammacore_lookup, gammatail_lookup)
      deallocate(Ncore_zp, sigma_knots_logE, sigma_knots_costheta, gammacore_knots_logE)
      deallocate(gammacore_knots_costheta, gammatail_knots_logE, gammatail_knots_costheta, PSF_logE)
      
      END SUBROUTINE


      SUBROUTINE flatIRFini_Aeff_clean

      deallocate(Aeff_costheta, Aeff_lookup, Aeff_logE, Aeff_knots_logE, Aeff_knots_costheta)

      END SUBROUTINE


      SUBROUTINE flatIRFini_Edisp_mean_clean

        deallocate(mean_Edisp_lookup, Oversampled_Edisp_logE, mean_Edisp_zp)
      
      END SUBROUTINE
      

      SUBROUTINE flatIRFini_PSF_mean_clean

        deallocate(mean_PSF_lookup, Oversampled_PSF_logE, Oversampled_PSF_unscaledx, mean_PSF_zp)

      END SUBROUTINE


      SUBROUTINE flatIRFini_Aeff_mean_clean

        deallocate(mean_Aeff_lookup, Oversampled_Aeff_logE, mean_Aeff_zp)

      END SUBROUTINE


      SUBROUTINE flatPrecompute_IRFs(IRF)

      character(len=*), intent(IN) :: IRF     

      write(*,*) '   Initializing...'
      call flatIRFini_initIRFs(IRF)
      write(*,*) '   Starting on Effective Area...'
      call flatPrecompute_Aeff(IRF)
      write(*,*) '   Finished Effective Area.  Starting on PSF...'
      call flatPrecompute_PSF(IRF)
      write(*,*) '   Finished PSF.  Starting on Energy Dispersion...'
      call flatPrecompute_Edisp(IRF)
      write(*,*) '   Finished Energy Dispersion.  Now just cleaning up...'
      call flatIRFini_cleanIRFs

      END SUBROUTINE


      SUBROUTINE flatPrecompute_FFTs(IRF, nFFTs, RA, DEC, RApix, DECpix, energyRange_in)

      character(len=*), intent(IN) :: IRF
      integer, intent(IN) :: nFFTs, RApix(nFFTs), DECpix(nFFTs)
      double precision, intent(IN) :: RA(nFFTs,2), DEC(nFFTs,2)
      double precision, intent(IN), optional :: energyRange_in(2)
      double precision :: energyRange(2)

      if (present(energyRange_in)) then
        energyRange = energyRange_in
      else 
        energyRange = (/Fermi_Emin, Fermi_Emax/)
      endif

      call flatIRFini_PSF_mean(IRF,.false.)
      call flatFFTW_Precompile(nFFTs, RA, DEC, RApix, DECpix, energyRange)
      call flatIRFini_PSF_mean_clean

      END SUBROUTINE


      SUBROUTINE flatIRFini_Edisp(IRF)
      ! Opens the FITS file containing parameter values for the LAT
      ! energy dispersion as a function of photon energy and incidence
      ! angle with respect to the spacecraft zenith, puts them in a 
      ! lookup table and initialises the interpolator.

      character(len=*), intent(IN) :: IRF
      double precision :: SVertices(1,2), SValue, SAbsErr
      double precision, allocatable :: logE_in(:,:), costheta_in(:,:), xbias_in(:)
      double precision, allocatable :: LS1_in(:), RS1_in(:), LS2_in(:), RS2_in(:), working(:)
      double precision, allocatable :: temp(:), NORM_edgeDerivs_logE(:), NORM_edgeDerivs_costheta(:)
      double precision :: Edisp_consts_in(9,2)
      integer :: Edisp_nEnergies_raw, Edisp_nTheta_raw, iflag, i, j, SRgType, IFAIL
      character (len=strlen) :: filename(2)

      filename(front) = trim(flatlib_datadir) // trim(IRF) // '/edisp_' // trim(IRF) // '_front.fits'
      filename(back) = trim(flatlib_datadir) // trim(IRF) // '/edisp_' // trim(IRF) // '_back.fits'
      detnam(front) = 'FRONT  '
      detnam(back) = 'BACK   '

      do strip=front,back,back-front

        ! Open FITS file
        fstatus=0
        call ftgiou(unit,fstatus)
        call ftopen(unit,filename(strip),0,blocksize,fstatus)
        call ftmrhd(unit,1,hdutype,fstatus)
 
        ! Check that the header looks roughly correct
        call ftgkys(unit,'TELESCOP',record,comment,fstatus)
        if (record .ne. 'GLAST  ') call flatUtils_crash('Incorrect TELESCOP header entry in energy dispersion FITS file.')
        call ftgkys(unit,'INSTRUME',record,comment,fstatus)
        if (record .ne. 'LAT    ') call flatUtils_crash('Incorrect INSTRUME header entry in energy dispersion FITS file.')
        call ftgkys(unit,'EXTNAME',record,comment,fstatus)
        if (record .ne. 'ENERGY DISPERSION') call flatUtils_crash('Incorrect EXTNAME header entry in energy dispersion FITS file.')
        call ftgkys(unit,'DETNAM',record,comment,fstatus)
        if (record .ne. detnam(strip)) call flatUtils_crash('Incorrect DETNAM header entry in energy dispersion FITS file.')

        ! Read the data dimensions      
        call fteqty(unit,1,datacode,Edisp_nEnergies_raw,datawidth,fstatus)
        call fteqty(unit,3,datacode,Edisp_nTheta_raw,datawidth,fstatus)
        Edisp_nEnergies = Edisp_nEnergies_raw + 2
        Edisp_nTheta = Edisp_nTheta_raw + 2

        if (strip .eq. front) then
          allocate(Edisp_logE(Edisp_nEnergies,2), Edisp_costheta(Edisp_nTheta,2))
          allocate(NORM_lookup(Edisp_nEnergies,Edisp_nTheta,2), xbias_lookup(Edisp_nEnergies,Edisp_nTheta,2))
          allocate(LS1_lookup(Edisp_nEnergies,Edisp_nTheta,2), RS1_lookup(Edisp_nEnergies,Edisp_nTheta,2))
          allocate(LS2_lookup(Edisp_nEnergies,Edisp_nTheta,2), RS2_lookup(Edisp_nEnergies,Edisp_nTheta,2))
          allocate(logE_in(Edisp_nEnergies_raw,2), costheta_in(Edisp_nTheta_raw,2))
          allocate(NORM_zp(3*Edisp_nEnergies*Edisp_nTheta,2), temp(2*Edisp_nTheta+Edisp_nEnergies))
          allocate(NORM_edgeDerivs_logE(Edisp_nEnergies), NORM_edgeDerivs_costheta(Edisp_nTheta))
          allocate(xbias_in(Edisp_nTheta_raw*Edisp_nEnergies_raw))
          allocate(xbias_knots_logE(Edisp_nEnergies+splineOrder_logE,2))
          allocate(xbias_knots_costheta(Edisp_nTheta+splineOrder_costheta,2))
          allocate(LS1_in(Edisp_nTheta_raw*Edisp_nEnergies_raw))
          allocate(LS1_knots_logE(Edisp_nEnergies+splineOrder_logE,2))
          allocate(LS1_knots_costheta(Edisp_nTheta+splineOrder_costheta,2))
          allocate(RS1_in(Edisp_nTheta_raw*Edisp_nEnergies_raw))
          allocate(RS1_knots_logE(Edisp_nEnergies+splineOrder_logE,2))
          allocate(RS1_knots_costheta(Edisp_nTheta+splineOrder_costheta,2))
          allocate(LS2_in(Edisp_nTheta_raw*Edisp_nEnergies_raw))
          allocate(LS2_knots_logE(Edisp_nEnergies+splineOrder_logE,2))
          allocate(LS2_knots_costheta(Edisp_nTheta+splineOrder_costheta,2))
          allocate(RS2_in(Edisp_nTheta_raw*Edisp_nEnergies_raw))
          allocate(RS2_knots_logE(Edisp_nEnergies+splineOrder_logE,2))
          allocate(RS2_knots_costheta(Edisp_nTheta+splineOrder_costheta,2))
          allocate(working(Edisp_nEnergies*Edisp_nTheta+2*max(splineOrder_logE*(Edisp_nEnergies+1), &
                   splineOrder_costheta*(Edisp_nTheta+1))))
        endif        

        ! Read the data
        felem = 1
        call ftgcvd(unit,1,1,felem,Edisp_nEnergies_raw,nulle,logE_in(:,1),anynull,fstatus)
        call ftgcvd(unit,2,1,felem,Edisp_nEnergies_raw,nulle,logE_in(:,2),anynull,fstatus)
        call ftgcvd(unit,3,1,felem,Edisp_nTheta_raw,nulle,costheta_in(:,1),anynull,fstatus)
        call ftgcvd(unit,4,1,felem,Edisp_nTheta_raw,nulle,costheta_in(:,2),anynull,fstatus)
        !call ftgcvd(unit,5,1,felem,Edisp_nEnergies_raw*Edisp_nTheta_raw,nulle,NORM_in,anynull,fstatus)
        call ftgcvd(unit,6,1,felem,Edisp_nEnergies_raw*Edisp_nTheta_raw,nulle,LS1_in,anynull,fstatus)
        call ftgcvd(unit,7,1,felem,Edisp_nEnergies_raw*Edisp_nTheta_raw,nulle,RS1_in,anynull,fstatus)
        call ftgcvd(unit,8,1,felem,Edisp_nEnergies_raw*Edisp_nTheta_raw,nulle,xbias_in,anynull,fstatus)
        call ftgcvd(unit,9,1,felem,Edisp_nEnergies_raw*Edisp_nTheta_raw,nulle,LS2_in,anynull,fstatus)
        call ftgcvd(unit,10,1,felem,Edisp_nEnergies_raw*Edisp_nTheta_raw,nulle,RS2_in,anynull,fstatus)
        if (fstatus .gt. 0) call flatUtils_crash('Read error in energy dispersion initialisation routine.')
        
        ! Read the energy scaling parameters
        call ftmrhd(unit,1,hdutype,fstatus)
        call ftgkys(unit,'EXTNAME',record,comment,fstatus)
        if (record .ne. 'EDISP_SCALING_PARAMS') call flatUtils_crash('Bad EXTNAME in scaling parameter section of energy dispersion FITS file.')
        call ftgcvd(unit,1,1,felem,9,nulle,Edisp_consts_in(:,strip),anynull,fstatus)
        if (fstatus .gt. 0) call flatUtils_crash('Scaling parameter read error in energy dispersion initialisation routine.')

        ! Close the FITS file
        call ftclos(unit, fstatus)
        call ftfiou(unit, fstatus)
        if (fstatus .gt. 0) call flatUtils_crash('Error encountered in closing energy dispersion FITS file.')

        ! Set the scaling constants 
        Edisp_consts(:,strip) = Edisp_consts_in(:,strip)

        ! Populate the bodies of the lookup tables
        Edisp_logE(2:Edisp_nEnergies-1,strip) = 0.5*log10(logE_in(:,1)*logE_in(:,2))
        Edisp_costheta(2:Edisp_nTheta-1,strip) = 0.5*(costheta_in(:,1) + costheta_in(:,2))

        xbias_lookup(2:Edisp_nEnergies-1,2:Edisp_nTheta-1,strip) = reshape(xbias_in, [Edisp_nEnergies_raw, Edisp_nTheta_raw])
        LS1_lookup(2:Edisp_nEnergies-1,2:Edisp_nTheta-1,strip) = reshape(LS1_in, [Edisp_nEnergies_raw, Edisp_nTheta_raw])
        RS1_lookup(2:Edisp_nEnergies-1,2:Edisp_nTheta-1,strip) = reshape(RS1_in, [Edisp_nEnergies_raw, Edisp_nTheta_raw])
        LS2_lookup(2:Edisp_nEnergies-1,2:Edisp_nTheta-1,strip) = reshape(LS2_in, [Edisp_nEnergies_raw, Edisp_nTheta_raw])
        RS2_lookup(2:Edisp_nEnergies-1,2:Edisp_nTheta-1,strip) = reshape(RS2_in, [Edisp_nEnergies_raw, Edisp_nTheta_raw])

        ! Pad the edges of the lookup tables with data from the edges of the first and last bins 
        Edisp_logE(1,strip) = log10(logE_in(1,1))
        Edisp_logE(Edisp_nEnergies,strip) = log10(logE_in(Edisp_nEnergies_raw,2))
        Edisp_costheta(1,strip) = costheta_in(1,1)
        Edisp_costheta(Edisp_nTheta,strip) = costheta_in(Edisp_nTheta_raw,2)

        xbias_lookup(1,2:Edisp_nTheta-1,strip) = xbias_lookup(2,2:Edisp_nTheta-1,strip)
        xbias_lookup(Edisp_nEnergies,2:Edisp_nTheta-1,strip) = xbias_lookup(Edisp_nEnergies-1,2:Edisp_nTheta-1,strip)
        xbias_lookup(:,1,strip) = xbias_lookup(:,2,strip)
        xbias_lookup(:,Edisp_nTheta,strip) = xbias_lookup(:,Edisp_nTheta-1,strip)
        LS1_lookup(1,2:Edisp_nTheta-1,strip) = LS1_lookup(2,2:Edisp_nTheta-1,strip)
        LS1_lookup(Edisp_nEnergies,2:Edisp_nTheta-1,strip) = LS1_lookup(Edisp_nEnergies-1,2:Edisp_nTheta-1,strip)
        LS1_lookup(:,1,strip) = LS1_lookup(:,2,strip)
        LS1_lookup(:,Edisp_nTheta,strip) = LS1_lookup(:,Edisp_nTheta-1,strip)
        RS1_lookup(1,2:Edisp_nTheta-1,strip) = RS1_lookup(2,2:Edisp_nTheta-1,strip)
        RS1_lookup(Edisp_nEnergies,2:Edisp_nTheta-1,strip) = RS1_lookup(Edisp_nEnergies-1,2:Edisp_nTheta-1,strip)
        RS1_lookup(:,1,strip) = RS1_lookup(:,2,strip)
        RS1_lookup(:,Edisp_nTheta,strip) = RS1_lookup(:,Edisp_nTheta-1,strip)
        LS2_lookup(1,2:Edisp_nTheta-1,strip) = LS2_lookup(2,2:Edisp_nTheta-1,strip)
        LS2_lookup(Edisp_nEnergies,2:Edisp_nTheta-1,strip) = LS2_lookup(Edisp_nEnergies-1,2:Edisp_nTheta-1,strip)
        LS2_lookup(:,1,strip) = LS2_lookup(:,2,strip)
        LS2_lookup(:,Edisp_nTheta,strip) = LS2_lookup(:,Edisp_nTheta-1,strip)
        RS2_lookup(1,2:Edisp_nTheta-1,strip) = RS2_lookup(2,2:Edisp_nTheta-1,strip)
        RS2_lookup(Edisp_nEnergies,2:Edisp_nTheta-1,strip) = RS2_lookup(Edisp_nEnergies-1,2:Edisp_nTheta-1,strip)
        RS2_lookup(:,1,strip) = RS2_lookup(:,2,strip)
        RS2_lookup(:,Edisp_nTheta,strip) = RS2_lookup(:,Edisp_nTheta-1,strip)

        ! Work out the real normalisation factors
        SVertices(1,:) = (/log10(Fermi_Emin)+3.d0, log10(Fermi_Emax)+3.d0/)
        SRgType = HyperQuad
        do j=1,Edisp_nTheta
          costheta_true_local = Edisp_costheta(j,strip)
          do i=1,Edisp_nEnergies
            l10E_true = Edisp_logE(i,strip)
            Edisp_pars = (/xbias_lookup(i,j,strip), LS1_lookup(i,j,strip), LS2_lookup(i,j,strip), &
                    RS1_lookup(i,j,strip), RS2_lookup(i,j,strip)/)
            IFAIL = 0
            call CUBATR(1,flatIRFini_Edisp_Integrand,SVertices,SRgType,SValue,SAbsErr,IFAIL,MaxPts=50000000,&
                        EpsAbs=integratorAbsErr,EpsRel=1.d-9)
            NORM_lookup(i,j,strip) = 1.d0/SValue
          enddo
        enddo

        ! Initialise the interpolator
        iflag = 0
        NORM_edgeDerivs_logE = 0.d0
        NORM_edgeDerivs_costheta = 0.d0   

        call surf1(Edisp_nEnergies,Edisp_nTheta,Edisp_logE(:,strip),Edisp_costheta(:,strip),NORM_lookup(:,:,strip),&
                   Edisp_nEnergies,NORM_edgeDerivs_costheta,NORM_edgeDerivs_costheta,NORM_edgeDerivs_logE,&
                   NORM_edgeDerivs_logE,0.d0,0.d0,0.d0,0.d0,0,NORM_zp(:,strip),temp,NORM_splineTension,iflag)
        if (iflag .ne. 0) call flatUtils_crash('surf1(FITPACK) returned error in energy dispersion NORM initialisation routine.')
        call DB2INK(Edisp_logE(:,strip),Edisp_nEnergies,Edisp_costheta(:,strip),Edisp_nTheta,xbias_lookup(:,:,strip), &
             Edisp_nEnergies,splineOrder_logE,splineOrder_costheta,xbias_knots_logE(:,strip), &
             xbias_knots_costheta(:,strip),xbias_lookup(:,:,strip),working,iflag)
        if (iflag .ne. 1) call flatUtils_crash('DB2INK returned error in energy dispersion xbias initialisation routine.')
        iflag = 0
        call DB2INK(Edisp_logE(:,strip),Edisp_nEnergies,Edisp_costheta(:,strip),Edisp_nTheta,LS1_lookup(:,:,strip), &
             Edisp_nEnergies,splineOrder_logE,splineOrder_costheta,LS1_knots_logE(:,strip), &
             LS1_knots_costheta(:,strip),LS1_lookup(:,:,strip),working,iflag)
        if (iflag .ne. 1) call flatUtils_crash('DB2INK returned error in energy dispersion LS1 initialisation routine.')
        iflag = 0
        call DB2INK(Edisp_logE(:,strip),Edisp_nEnergies,Edisp_costheta(:,strip),Edisp_nTheta,RS1_lookup(:,:,strip), &
             Edisp_nEnergies,splineOrder_logE,splineOrder_costheta,RS1_knots_logE(:,strip), &
             RS1_knots_costheta(:,strip),RS1_lookup(:,:,strip),working,iflag)
        if (iflag .ne. 1) call flatUtils_crash('DB2INK returned error in energy dispersion RS1 initialisation routine.')
        iflag = 0
        call DB2INK(Edisp_logE(:,strip),Edisp_nEnergies,Edisp_costheta(:,strip),Edisp_nTheta,LS2_lookup(:,:,strip), &
             Edisp_nEnergies,splineOrder_logE,splineOrder_costheta,LS2_knots_logE(:,strip), &
             LS2_knots_costheta(:,strip),LS2_lookup(:,:,strip),working,iflag)
        if (iflag .ne. 1) call flatUtils_crash('DB2INK returned error in energy dispersion LS2 initialisation routine.')
        iflag = 0
        call DB2INK(Edisp_logE(:,strip),Edisp_nEnergies,Edisp_costheta(:,strip),Edisp_nTheta,RS2_lookup(:,:,strip), &
             Edisp_nEnergies,splineOrder_logE,splineOrder_costheta,RS2_knots_logE(:,strip), &
             RS2_knots_costheta(:,strip),RS2_lookup(:,:,strip),working,iflag)
        if (iflag .ne. 1) call flatUtils_crash('DB2INK returned error in energy dispersion RS2 initialisation routine.')
     
      end do     

      deallocate(logE_in, costheta_in, xbias_in, LS1_in, RS1_in, LS2_in, RS2_in, working)
      deallocate(temp, NORM_edgeDerivs_logE, NORM_edgeDerivs_costheta)

      return

      END SUBROUTINE flatIRFini_Edisp


      FUNCTION flatIRFini_Edisp_Integrand(NumFun,X) RESULT(Value)
      ! Input:    
      ! Output:   
        
      integer, intent(IN) :: NumFun
      double precision, intent(IN) :: X(:)

      double precision :: Value(NumFun)

      Value(1) = flatIRFs_Edisp_rawfunc(X(1)-3.d0, l10E_true-3.d0, costheta_true_local, Edisp_pars, strip) * 10.**(X(1)) * log(10.d0)

      END FUNCTION flatIRFini_Edisp_Integrand


      SUBROUTINE flatIRFini_PSF(IRF)
      ! Opens the FITS file containing parameter values for the LAT
      ! point spread function at various photon energies and incidence
      ! angles with respect to the spacecraft zenith, puts them in a 
      ! lookup table and initialises the interpolator.

      character(len=*), intent(IN) :: IRF
      double precision :: SVertices(2,3), SValue, SAbsErr, PSF_Phi(PSF_nPhi), normSum
      double precision, allocatable :: logE_in(:,:), costheta_in(:,:), sigma_in(:), sigmatail_in(:)
      double precision, allocatable :: Ntail_in(:), gammacore_in(:), gammatail_in(:), working(:)
      double precision, allocatable :: temp(:), Ncore_edgeDerivs_logE(:), Ncore_edgeDerivs_costheta(:)
      double precision, allocatable :: Ntail_edgeDerivs_logE(:), Ntail_edgeDerivs_costheta(:)
      double precision :: PSF_consts_in(5)
      integer :: iflag, tfields, PSF_nEnergies_raw, PSF_nTheta_raw, i, j, k, SRgType, IFAIL
      character (len=strlen) :: filename(2)

      filename(front) = trim(flatlib_datadir) // trim(IRF) // '/psf_' // trim(IRF) // '_front.fits'
      filename(back) = trim(flatlib_datadir) // trim(IRF) // '/psf_' // trim(IRF) // '_back.fits'
      detnam(front) = 'FRONT  '
      detnam(back) = 'BACK   '

      forall(i=1:PSF_nPhi) PSF_Phi(i) = dble(i-1)/dble(PSF_nPhi)*2.d0*pi_fl

      do strip=front,back,back-front

        ! Open FITS file
        fstatus=0
        call ftgiou(unit,fstatus)
        call ftopen(unit,filename(strip),0,blocksize,fstatus)
        call ftmrhd(unit,1,hdutype,fstatus)

        ! Check that the header looks roughly correct
        call ftgkys(unit,'TELESCOP',record,comment,fstatus)
        if (record .ne. 'GLAST  ') call flatUtils_crash('Incorrect TELESCOP header entry in PSF FITS file.')
        call ftgkys(unit,'INSTRUME',record,comment,fstatus)
        if (record .ne. 'LAT    ') call flatUtils_crash('Incorrect INSTRUME header entry in PSF FITS file.')
        call ftgkys(unit,'EXTNAME',record,comment,fstatus)
        if (record .ne. 'RPSF') call flatUtils_crash('Incorrect EXTNAME header entry in PSF FITS file.')
        call ftgkys(unit,'DETNAM',record,comment,fstatus)
        if (record .ne. detnam(strip)) call flatUtils_crash('Incorrect DETNAM header entry in PSF FITS file.')

        ! Work out which PSF parameterisation these IRFs use
        call ftgkyj(unit,'TFIELDS',tfields,comment,fstatus)
        select case (tfields)
          case (8)
            PSFparameterisation = DblKingSingleSigma
          case (10)
            PSFparameterisation = DblKingDblSigma
          case default
            call flatUtils_crash('Unrecognised PSF parameterisation in PSF FITS file.')
        end select

        ! Read the data dimensions      
        call fteqty(unit,1,datacode,PSF_nEnergies_raw,datawidth,fstatus)
        call fteqty(unit,3,datacode,PSF_nTheta_raw,datawidth,fstatus)
        PSF_nEnergies = PSF_nEnergies_raw + 2
        PSF_nTheta = PSF_nTheta_raw + 2

        if (strip .eq. front) then
          allocate(PSF_logE(PSF_nEnergies,2), PSF_costheta(PSF_nTheta,2))
          allocate(Ncore_lookup(PSF_nEnergies,PSF_nTheta,2), sigma_lookup(PSF_nEnergies,PSF_nTheta,2))
          allocate(gammacore_lookup(PSF_nEnergies,PSF_nTheta,2), gammatail_lookup(PSF_nEnergies,PSF_nTheta,2))
          allocate(logE_in(PSF_nEnergies_raw,2), costheta_in(PSF_nTheta_raw,2))
          allocate(Ncore_zp(3*PSF_nEnergies*PSF_nTheta,2), temp(2*PSF_nTheta+PSF_nEnergies))
          allocate(Ncore_edgeDerivs_logE(PSF_nEnergies), Ncore_edgeDerivs_costheta(PSF_nTheta))
          allocate(sigma_in(PSF_nTheta_raw*PSF_nEnergies_raw))
          allocate(sigma_knots_logE(PSF_nEnergies+splineOrder_logE,2))
          allocate(sigma_knots_costheta(PSF_nTheta+splineOrder_costheta,2))
          allocate(gammacore_in(PSF_nTheta_raw*PSF_nEnergies_raw))
          allocate(gammacore_knots_logE(PSF_nEnergies+splineOrder_logE,2))
          allocate(gammacore_knots_costheta(PSF_nTheta+splineOrder_costheta,2))
          allocate(gammatail_in(PSF_nTheta_raw*PSF_nEnergies_raw))
          allocate(gammatail_knots_logE(PSF_nEnergies+splineOrder_logE,2))
          allocate(gammatail_knots_costheta(PSF_nTheta+splineOrder_costheta,2))
          allocate(working(PSF_nEnergies*PSF_nTheta+2*max(splineOrder_logE*(PSF_nEnergies+1), &
                   splineOrder_costheta*(PSF_nTheta+1))))
          if (PSFparameterisation .eq. DblKingDblSigma) then
            allocate(Ntail_lookup(PSF_nEnergies,PSF_nTheta,2), sigmatail_lookup(PSF_nEnergies,PSF_nTheta,2))
            allocate(Ntail_zp(3*PSF_nEnergies*PSF_nTheta,2))
            allocate(Ntail_edgeDerivs_logE(PSF_nEnergies), Ntail_edgeDerivs_costheta(PSF_nTheta))
            allocate(Ntail_in(PSF_nTheta_raw*PSF_nEnergies_raw))
            allocate(sigmatail_in(PSF_nTheta_raw*PSF_nEnergies_raw))
            allocate(sigmatail_knots_logE(PSF_nEnergies+splineOrder_logE,2))
            allocate(sigmatail_knots_costheta(PSF_nTheta+splineOrder_costheta,2))
          endif 
        endif        

        ! Read the data
        felem = 1
        call ftgcvd(unit,1,1,felem,PSF_nEnergies_raw,nulle,logE_in(:,1),anynull,fstatus)
        call ftgcvd(unit,2,1,felem,PSF_nEnergies_raw,nulle,logE_in(:,2),anynull,fstatus)
        call ftgcvd(unit,3,1,felem,PSF_nTheta_raw,nulle,costheta_in(:,1),anynull,fstatus)
        call ftgcvd(unit,4,1,felem,PSF_nTheta_raw,nulle,costheta_in(:,2),anynull,fstatus)
        select case (PSFparameterisation)
          case (DblKingSingleSigma)
            call ftgcvd(unit,6,1,felem,PSF_nEnergies_raw*PSF_nTheta_raw,nulle,sigma_in,anynull,fstatus)
            call ftgcvd(unit,7,1,felem,PSF_nEnergies_raw*PSF_nTheta_raw,nulle,gammacore_in,anynull,fstatus)
            call ftgcvd(unit,8,1,felem,PSF_nEnergies_raw*PSF_nTheta_raw,nulle,gammatail_in,anynull,fstatus)
          case (DBlKingDblSigma)
            call ftgcvd(unit,6,1,felem,PSF_nEnergies_raw*PSF_nTheta_raw,nulle,Ntail_in,anynull,fstatus)
            call ftgcvd(unit,7,1,felem,PSF_nEnergies_raw*PSF_nTheta_raw,nulle,sigma_in,anynull,fstatus)
            call ftgcvd(unit,8,1,felem,PSF_nEnergies_raw*PSF_nTheta_raw,nulle,sigmatail_in,anynull,fstatus)
            call ftgcvd(unit,9,1,felem,PSF_nEnergies_raw*PSF_nTheta_raw,nulle,gammacore_in,anynull,fstatus)
            call ftgcvd(unit,10,1,felem,PSF_nEnergies_raw*PSF_nTheta_raw,nulle,gammatail_in,anynull,fstatus)
                end select 
     
        ! Read the energy scaling parameters
        call ftmrhd(unit,1,hdutype,fstatus)
        call ftgkys(unit,'EXTNAME',record,comment,fstatus)
        if (record .ne. 'PSF_SCALING_PARAMS') call flatUtils_crash('Bad EXTNAME in scaling parameter section of PSF FITS file.')
        call ftgcvd(unit,1,1,felem,5,nulle,PSF_consts_in(:),anynull,fstatus)
        if (fstatus .gt. 0) call flatUtils_crash('Scaling parameter read error in PSF initialisation routine.')
        if (strip .eq. front) then
          PSF_consts(:,front) = PSF_consts_in((/1,2,5/))
        else
          PSF_consts(:,back) = PSF_consts_in(3:5)
        endif

        ! Close the FITS file
        call ftclos(unit, fstatus)
        call ftfiou(unit, fstatus)
        if (fstatus .gt. 0) call flatUtils_crash('Error encountered in reading PSF FITS file.')

        ! Populate the bodies of the lookup tables
        PSF_logE(2:PSF_nEnergies-1,strip) = 0.5*log10(logE_in(:,1)*logE_in(:,2))
        PSF_costheta(2:PSF_nTheta-1,strip) = 0.5*(costheta_in(:,1) + costheta_in(:,2))

        sigma_lookup(2:PSF_nEnergies-1,2:PSF_nTheta-1,strip) = reshape(sigma_in, [PSF_nEnergies_raw, PSF_nTheta_raw])
        gammacore_lookup(2:PSF_nEnergies-1,2:PSF_nTheta-1,strip) = reshape(gammacore_in, [PSF_nEnergies_raw, PSF_nTheta_raw])
        gammatail_lookup(2:PSF_nEnergies-1,2:PSF_nTheta-1,strip) = reshape(gammatail_in, [PSF_nEnergies_raw, PSF_nTheta_raw])
        if (PSFparameterisation .eq. DblKingDblSigma) then
          Ntail_lookup(2:PSF_nEnergies-1,2:PSF_nTheta-1,strip) = reshape(Ntail_in, [PSF_nEnergies_raw, PSF_nTheta_raw])
          sigmatail_lookup(2:PSF_nEnergies-1,2:PSF_nTheta-1,strip) = reshape(sigmatail_in, [PSF_nEnergies_raw, PSF_nTheta_raw])
        endif
          
        ! Pad the edges of the lookup table with data from the edges of the first and last bins 
        PSF_logE(1,strip) = log10(logE_in(1,1))
        PSF_logE(PSF_nEnergies,strip) = log10(logE_in(PSF_nEnergies_raw,2))
        PSF_costheta(1,strip) = costheta_in(1,1)
        PSF_costheta(PSF_nTheta,strip) = costheta_in(PSF_nTheta_raw,2)

        sigma_lookup(1,2:PSF_nTheta-1,strip) = sigma_lookup(2,2:PSF_nTheta-1,strip)
        sigma_lookup(PSF_nEnergies,2:PSF_nTheta-1,strip) = sigma_lookup(PSF_nEnergies-1,2:PSF_nTheta-1,strip)
        sigma_lookup(:,1,strip) = sigma_lookup(:,2,strip)
        sigma_lookup(:,PSF_nTheta,strip) = sigma_lookup(:,PSF_nTheta-1,strip)
        gammacore_lookup(1,2:PSF_nTheta-1,strip) = gammacore_lookup(2,2:PSF_nTheta-1,strip)
        gammacore_lookup(PSF_nEnergies,2:PSF_nTheta-1,strip) = gammacore_lookup(PSF_nEnergies-1,2:PSF_nTheta-1,strip)
        gammacore_lookup(:,1,strip) = gammacore_lookup(:,2,strip)
        gammacore_lookup(:,PSF_nTheta,strip) = gammacore_lookup(:,PSF_nTheta-1,strip)
        gammatail_lookup(1,2:PSF_nTheta-1,strip) = gammatail_lookup(2,2:PSF_nTheta-1,strip)
        gammatail_lookup(PSF_nEnergies,2:PSF_nTheta-1,strip) = gammatail_lookup(PSF_nEnergies-1,2:PSF_nTheta-1,strip)
        gammatail_lookup(:,1,strip) = gammatail_lookup(:,2,strip)
        gammatail_lookup(:,PSF_nTheta,strip) = gammatail_lookup(:,PSF_nTheta-1,strip)
        if (PSFparameterisation .eq. DblKingDblSigma) then
          sigmatail_lookup(1,2:PSF_nTheta-1,strip) = sigmatail_lookup(2,2:PSF_nTheta-1,strip)
          sigmatail_lookup(PSF_nEnergies,2:PSF_nTheta-1,strip) = sigmatail_lookup(PSF_nEnergies-1,2:PSF_nTheta-1,strip)
          sigmatail_lookup(:,1,strip) = sigmatail_lookup(:,2,strip)
          sigmatail_lookup(:,PSF_nTheta,strip) = sigmatail_lookup(:,PSF_nTheta-1,strip)
          Ntail_lookup(1,2:PSF_nTheta-1,strip) = Ntail_lookup(2,2:PSF_nTheta-1,strip)
          Ntail_lookup(PSF_nEnergies,2:PSF_nTheta-1,strip) = Ntail_lookup(PSF_nEnergies-1,2:PSF_nTheta-1,strip)
          Ntail_lookup(:,1,strip) = Ntail_lookup(:,2,strip)
          Ntail_lookup(:,PSF_nTheta,strip) = Ntail_lookup(:,PSF_nTheta-1,strip)
        endif

        ! Work out the real normalisation factors
        SVertices(:,1) = (/acos(Fermi_costhetamax), Fermi_phimin/)
        SVertices(:,2) = (/acos(Fermi_costhetamax), Fermi_phimax/)
        SVertices(:,3) = (/acos(Fermi_costhetamin), Fermi_phimin/)
        SRgType = HyperQuad
        do j=1,PSF_nTheta
          costheta_true_local = PSF_costheta(j,strip)  
          do i=1,PSF_nEnergies
            l10E_true = PSF_logE(i,strip)
            select case (PSFparameterisation)
              case (DblKingSingleSigma)
                PSF_pars = (/sigma_lookup(i,j,strip), gammacore_lookup(i,j,strip), gammatail_lookup(i,j,strip), 0.d0, 0.d0/)
               case (DblKingDblSigma)
                PSF_pars = (/Ntail_lookup(i,j,strip), sigma_lookup(i,j,strip), &
                            sigmatail_lookup(i,j,strip), gammacore_lookup(i,j,strip), gammatail_lookup(i,j,strip)/)
            end select
            normSum = 0.d0
            do k=1,PSF_nPhi
              phi_true_local = PSF_Phi(k)
              IFAIL = 0
              call CUBATR(2,flatIRFini_PSF_Integrand,SVertices,SRgType,SValue,SAbsErr,IFAIL,MaxPts=5000000,&
                        EpsAbs=integratorAbsErr,EpsRel=1.d-5)
              normSum = normSum + SValue
            enddo
            Ncore_lookup(i,j,strip) = PSF_nPhi/normSum * radperdeg * radperdeg
          enddo
        enddo
        
        ! Initialise the interpolator
        iflag = 0
        Ncore_edgeDerivs_logE = 0.d0
        Ncore_edgeDerivs_costheta = 0.d0   
        call surf1(PSF_nEnergies,PSF_nTheta,PSF_logE(:,strip),PSF_costheta(:,strip),Ncore_lookup(:,:,strip),&
                   PSF_nEnergies,Ncore_edgeDerivs_costheta,Ncore_edgeDerivs_costheta,Ncore_edgeDerivs_logE,&
                   Ncore_edgeDerivs_logE,0.d0,0.d0,0.d0,0.d0,0,Ncore_zp(:,strip),temp,Ncore_splineTension,iflag)
        if (iflag .ne. 0) call flatUtils_crash('surf1(FITPACK) returned error in PSF Ncore initialisation routine.')
        call DB2INK(PSF_logE(:,strip),PSF_nEnergies,PSF_costheta(:,strip),PSF_nTheta,sigma_lookup(:,:,strip), &
             PSF_nEnergies,splineOrder_logE,splineOrder_costheta,sigma_knots_logE(:,strip), &
             sigma_knots_costheta(:,strip),sigma_lookup(:,:,strip),working,iflag)
        if (iflag .ne. 1) call flatUtils_crash('DB2INK returned error in PSF sigma initialisation routine.')
        iflag = 0
        call DB2INK(PSF_logE(:,strip),PSF_nEnergies,PSF_costheta(:,strip),PSF_nTheta,gammacore_lookup(:,:,strip), &
             PSF_nEnergies,splineOrder_logE,splineOrder_costheta,gammacore_knots_logE(:,strip), &
             gammacore_knots_costheta(:,strip),gammacore_lookup(:,:,strip),working,iflag)
        if (iflag .ne. 1) call flatUtils_crash('DB2INK returned error in PSF gammacore initialisation routine.')
        iflag = 0
        call DB2INK(PSF_logE(:,strip),PSF_nEnergies,PSF_costheta(:,strip),PSF_nTheta,gammatail_lookup(:,:,strip), &
             PSF_nEnergies,splineOrder_logE,splineOrder_costheta,gammatail_knots_logE(:,strip), &
             gammatail_knots_costheta(:,strip),gammatail_lookup(:,:,strip),working,iflag)
        if (iflag .ne. 1) call flatUtils_crash('DB2INK returned error in PSF gammatail initialisation routine.')

        if (PSFparameterisation .eq. DblKingDblSigma) then
          iflag = 0
          Ntail_edgeDerivs_logE = 0.d0
          Ntail_edgeDerivs_costheta = 0.d0   
          call surf1(PSF_nEnergies,PSF_nTheta,PSF_logE(:,strip),PSF_costheta(:,strip),Ntail_lookup(:,:,strip),&
                   PSF_nEnergies,Ntail_edgeDerivs_costheta,Ntail_edgeDerivs_costheta,Ntail_edgeDerivs_logE,&
                   Ntail_edgeDerivs_logE,0.d0,0.d0,0.d0,0.d0,0,Ntail_zp(:,strip),temp,Ntail_splineTension,iflag)
          if (iflag .ne. 0) call flatUtils_crash('surf1(FITPACK) returned error in PSF Ntail initialisation routine.')
          call DB2INK(PSF_logE(:,strip),PSF_nEnergies,PSF_costheta(:,strip),PSF_nTheta,sigmatail_lookup(:,:,strip), &
             PSF_nEnergies,splineOrder_logE,splineOrder_costheta,sigmatail_knots_logE(:,strip), &
             sigmatail_knots_costheta(:,strip),sigmatail_lookup(:,:,strip),working,iflag)
          if (iflag .ne. 1) call flatUtils_crash('DB2INK returned error in PSF sigmatail initialisation routine.')
        endif
     
      end do

      deallocate(logE_in, costheta_in, sigma_in, gammacore_in, gammatail_in, working)
      deallocate(temp, Ncore_edgeDerivs_logE, Ncore_edgeDerivs_costheta)
      if (PSFparameterisation .eq. DblKingDblSigma) deallocate(Ntail_in, sigmatail_in, Ntail_edgeDerivs_logE, Ntail_edgeDerivs_costheta)
     
      return

      END SUBROUTINE flatIRFini_PSF


      FUNCTION flatIRFini_PSF_Integrand(NumFun,X) RESULT(Value)
      ! Input:    
      ! Output:   

      integer, intent(IN) :: NumFun
      double precision, intent(IN) :: X(:)

      double precision :: Value(NumFun), IncidenceAngles_obs(2), IncidenceAngles_true(2)

      IncidenceAngles_obs = (/cos(X(1)), X(2)/)
      IncidenceAngles_true = (/costheta_true_local, phi_true_local/)

      Value(1) = flatIRFs_PSF_rawfunc(IncidenceAngles_obs, l10E_true - 3.d0, IncidenceAngles_true, PSF_pars, strip)

      END FUNCTION flatIRFini_PSF_Integrand


      SUBROUTINE flatIRFini_Aeff(IRF)
      ! Opens the FITS file containing the effective area of the Fermi
      ! LAT as a function of photon energy and incidence angle
      ! with respect to the spacecraft zenith, puts them in a 
      ! lookup table and initialises the interpolator.'

      character(len=*), intent(IN) :: IRF     
      double precision, allocatable :: logE_in(:,:), costheta_in(:,:), aeff_in(:), working(:)
      integer :: iflag, Aeff_nEnergies_raw, Aeff_nTheta_raw
      character (len=strlen) :: filename(2)

      filename(front) = trim(flatlib_datadir) // trim(IRF) // '/aeff_' // trim(IRF) // '_front.fits'
      filename(back) = trim(flatlib_datadir) // trim(IRF) // '/aeff_' // trim(IRF) // '_back.fits'
      detnam(front) = 'FRONT  '
      detnam(back) = 'BACK   '

      do strip=front,back,back-front

        ! Open FITS file
        fstatus=0
        call ftgiou(unit,fstatus)
        call ftopen(unit,filename(strip),0,blocksize,fstatus)
        call ftmrhd(unit,1,hdutype,fstatus)

        ! Check that the header looks roughly correct
        call ftgkys(unit,'TELESCOP',record,comment,fstatus)
        if (record .ne. 'GLAST  ') call flatUtils_crash('Incorrect TELESCOP header entry in effective area FITS file.')
        call ftgkys(unit,'INSTRUME',record,comment,fstatus)
        if (record .ne. 'LAT    ') call flatUtils_crash('Incorrect INSTRUME header entry in effective area FITS file.')
        call ftgkys(unit,'EXTNAME',record,comment,fstatus)
        if (record .ne. 'EFFECTIVE AREA') call flatUtils_crash('Incorrect EXTNAME header entry in effective area FITS file.')
        call ftgkys(unit,'DETNAM',record,comment,fstatus)
        if (record .ne. detnam(strip)) call flatUtils_crash('Incorrect DETNAM header entry in effective area FITS file.')

        ! Read the data dimensions      
        call fteqty(unit,1,datacode,Aeff_nEnergies_raw,datawidth,fstatus)
        call fteqty(unit,3,datacode,Aeff_nTheta_raw,datawidth,fstatus)
        Aeff_nEnergies = Aeff_nEnergies_raw + 2
        Aeff_nTheta = Aeff_nTheta_raw + 2

        if (strip .eq. front) then
          allocate(Aeff_logE(Aeff_nEnergies,2), Aeff_costheta(Aeff_nTheta,2), Aeff_lookup(Aeff_nEnergies,Aeff_nTheta,2))
          allocate(logE_in(Aeff_nEnergies_raw,2), costheta_in(Aeff_nTheta_raw,2))
          allocate(aeff_in(Aeff_nTheta_raw*Aeff_nEnergies_raw))
          allocate(Aeff_knots_logE(Aeff_nEnergies+splineOrder_logE,2))
          allocate(Aeff_knots_costheta(Aeff_nTheta+splineOrder_costheta,2))
          allocate(working(Aeff_nEnergies*Aeff_nTheta+2*max(splineOrder_logE*(Aeff_nEnergies+1), &
                   splineOrder_costheta*(Aeff_nTheta+1))))
        endif        

        ! Read the data
        felem = 1
        call ftgcvd(unit,1,1,felem,Aeff_nEnergies_raw,nulle,logE_in(:,1),anynull,fstatus)
        call ftgcvd(unit,2,1,felem,Aeff_nEnergies_raw,nulle,logE_in(:,2),anynull,fstatus)
        call ftgcvd(unit,3,1,felem,Aeff_nTheta_raw,nulle,costheta_in(:,1),anynull,fstatus)
        call ftgcvd(unit,4,1,felem,Aeff_nTheta_raw,nulle,costheta_in(:,2),anynull,fstatus)
        call ftgcvd(unit,5,1,felem,Aeff_nEnergies_raw*Aeff_nTheta_raw,nulle,aeff_in,anynull,fstatus)

        ! Populate the body of the lookup table
        Aeff_logE(2:Aeff_nEnergies-1,strip) = 0.5*log10(logE_in(:,1)*logE_in(:,2))
        Aeff_costheta(2:Aeff_nTheta-1,strip) = 0.5*(costheta_in(:,1) + costheta_in(:,2))
        Aeff_lookup(2:Aeff_nEnergies-1,2:Aeff_nTheta-1,strip) = reshape(aeff_in, [Aeff_nEnergies_raw, Aeff_nTheta_raw])

        ! Pad the edges of the lookup table with data from the edges of the first and last bins 
        Aeff_logE(1,strip) = log10(logE_in(1,1))
        Aeff_logE(Aeff_nEnergies,strip) = log10(logE_in(Aeff_nEnergies_raw,2))

        Aeff_costheta(1,strip) = costheta_in(1,1)
        Aeff_costheta(Aeff_nTheta,strip) = costheta_in(Aeff_nTheta_raw,2)

        Aeff_lookup(1,2:Aeff_nTheta-1,strip) = Aeff_lookup(2,2:Aeff_nTheta-1,strip)
        Aeff_lookup(Aeff_nEnergies,2:Aeff_nTheta-1,strip) = Aeff_lookup(Aeff_nEnergies-1,2:Aeff_nTheta-1,strip)
        Aeff_lookup(:,1,strip) = Aeff_lookup(:,2,strip)
        Aeff_lookup(:,Aeff_nTheta,strip) = Aeff_lookup(:,Aeff_nTheta-1,strip)

        ! Initialise the interpolator
        iflag = 0
        call DB2INK(Aeff_logE(:,strip),Aeff_nEnergies,Aeff_costheta(:,strip),Aeff_nTheta,Aeff_lookup(:,:,strip), &
             Aeff_nEnergies,splineOrder_logE,splineOrder_costheta,Aeff_knots_logE(:,strip), &
             Aeff_knots_costheta(:,strip),Aeff_lookup(:,:,strip),working,iflag)
        if (iflag .ne. 1) call flatUtils_crash('DB2INK returned error in effective area initialisation routine.')

        ! Close the FITS file
        call ftclos(unit, fstatus)
        call ftfiou(unit, fstatus)
        if (fstatus .gt. 0) call flatUtils_crash('Error encountered in reading effective area FITS file.')
     
      end do

      deallocate(logE_in, costheta_in, aeff_in, working)

      return

      END SUBROUTINE flatIRFini_Aeff


      SUBROUTINE flatIRFini_Edisp_mean(IRF,balkIfMissing)

      character(len=*), intent(IN) :: IRF
      logical, intent(IN) :: balkIfMissing
      integer :: ierr, fileunit
      double precision, allocatable :: mean_Edisp_edge_derivatives(:), temp(:)
      character (len=strlen) :: infilename   

      ierr=0
      call ftgiou(fileunit,ierr)
      if (ierr .ne. 0) call flatUtils_crash('No free file handle available for flatIRFini_Edisp_mean!')

      infilename = trim(flatlib_datadir) // trim(IRF) // '/edisp_' // trim(IRF) // '_mean.dat'
      do
        open(fileunit, STATUS='OLD', FILE=infilename, IOSTAT=ierr, ACTION='READ')
        if (ierr .ne. 0) then
          if (balkIfMissing) call flatUtils_crash('Could not open mean energy dispersion file for '// trim(IRF) //' - does it exist?')
          write(*,*) '  Damn -- mean IRFs not generated yet.  This might take a while...'
          call flatPrecompute_IRFs(IRF)
        else 
          exit
        endif
      end do

      read(fileunit, '(I10)') Oversampled_Edisp_nEnergies

      allocate(Oversampled_Edisp_logE(Oversampled_Edisp_nEnergies,2))
      allocate(mean_Edisp_lookup(Oversampled_Edisp_nEnergies, Oversampled_Edisp_nEnergies,2))
      allocate(mean_Edisp_zp(3*Oversampled_Edisp_nEnergies*Oversampled_Edisp_nEnergies,2))
      allocate(mean_Edisp_edge_derivatives(Oversampled_Edisp_nEnergies))
      allocate(temp(3*Oversampled_Edisp_nEnergies))
    
      read(fileunit, '(ES13.6)') Oversampled_Edisp_logE
      read(fileunit, '(ES13.6)') mean_Edisp_lookup

      close(fileunit)    

      mean_Edisp_edge_derivatives = 0.d0

      do strip=front,back,back-front

        !Initialise interpolator
        call surf1(Oversampled_Edisp_nEnergies,Oversampled_Edisp_nEnergies,Oversampled_Edisp_logE(:,strip),&
                   Oversampled_Edisp_logE(:,strip),mean_Edisp_lookup(:,:,strip),Oversampled_Edisp_nEnergies,&
                   mean_Edisp_edge_derivatives,mean_Edisp_edge_derivatives,mean_Edisp_edge_derivatives,&
                   mean_Edisp_edge_derivatives,0.d0,0.d0,0.d0,0.d0,255,mean_Edisp_zp(:,strip),temp,Edisp_splineTension,ierr)
        if (ierr .ne. 0) call flatUtils_crash('surf1(FITPACK) returned error in mean Edisp initialisation routine.')

      enddo

      deallocate(mean_Edisp_edge_derivatives, temp)

      END SUBROUTINE flatIRFini_Edisp_mean


      SUBROUTINE flatIRFini_PSF_mean(IRF,balkIfMissing)

      character(len=*), intent(IN) :: IRF  
      logical, intent(IN) :: balkIfMissing  
      integer :: ierr, fileunit
      double precision, allocatable :: mean_PSF_edge_derivatives(:), temp(:)
      character (len=strlen) :: infilename   

      ierr=0
      call ftgiou(fileunit,ierr)
      if (ierr .ne. 0) call flatUtils_crash('No free file handle available for flatIRFini_PSF_mean!')

      infilename = trim(flatlib_datadir) // trim(IRF) // '/psf_' // trim(IRF) // '_mean.dat'
      do
        open(fileunit, STATUS='OLD', FILE=infilename, IOSTAT=ierr, ACTION='READ')
        if (ierr .ne. 0) then
          if (balkIfMissing) call flatUtils_crash('Could not open mean PSF file for '// trim(IRF) //' - does it exist?')
          write(*,*) '  Damn -- mean IRFs not generated yet.  This might take a while...'
          call flatPrecompute_IRFs(IRF)
        else 
          exit
        endif
      end do

      read(fileunit, '(I10)') Oversampled_PSF_nEnergies

      allocate(Oversampled_PSF_logE(Oversampled_PSF_nEnergies,2))
      allocate(Oversampled_PSF_unscaledx(Oversampled_PSF_nEnergies,2))
      allocate(mean_PSF_lookup(Oversampled_PSF_nEnergies, Oversampled_PSF_nEnergies,2))
      allocate(mean_PSF_zp(3*Oversampled_PSF_nEnergies*Oversampled_PSF_nEnergies,2))
      allocate(mean_PSF_edge_derivatives(Oversampled_PSF_nEnergies))
      allocate(temp(3*Oversampled_PSF_nEnergies))
    
      read(fileunit, '(ES13.6)') Oversampled_PSF_logE
      read(fileunit, '(ES13.6)') Oversampled_PSF_unscaledx
      read(fileunit, '(ES13.6)') mean_PSF_lookup

      close(fileunit)    

      mean_PSF_edge_derivatives = 0.d0

      do strip=front,back,back-front

        !Initialise interpolator
        call surf1(Oversampled_PSF_nEnergies,Oversampled_PSF_nEnergies,Oversampled_PSF_logE(:,strip),&
                   Oversampled_PSF_unscaledx(:,strip),mean_PSF_lookup(:,:,strip),Oversampled_PSF_nEnergies,&
                   mean_PSF_edge_derivatives,mean_PSF_edge_derivatives,mean_PSF_edge_derivatives,&
                   mean_PSF_edge_derivatives,0.d0,0.d0,0.d0,0.d0,255,mean_PSF_zp(:,strip),temp,PSF_splineTension,ierr)
        if (ierr .ne. 0) call flatUtils_crash('surf1(FITPACK) returned error in mean PSF initialisation routine.')

      enddo    

      deallocate(mean_PSF_edge_derivatives, temp)

      END SUBROUTINE flatIRFini_PSF_mean


      SUBROUTINE flatIRFini_Aeff_mean(IRF,balkIfMissing)

      character(len=*), intent(IN) :: IRF
      logical, intent(IN) :: balkIfMissing
      integer :: ierr, fileunit
      double precision, allocatable :: temp(:)
      character (len=strlen) :: infilename   

      ierr=0
      call ftgiou(fileunit,ierr)
      if (ierr .ne. 0) call flatUtils_crash('No free file handle available for flatIRFini_Aeff_mean!')

      infilename = trim(flatlib_datadir) // trim(IRF) // '/aeff_' // trim(IRF) // '_mean.dat'
      do
        open(fileunit, STATUS='OLD', FILE=infilename, IOSTAT=ierr, ACTION='READ')
        if (ierr .ne. 0) then
          if (balkIfMissing)  call flatUtils_crash('Could not open mean effective area file for '// trim(IRF) //' - does it exist?')
          write(*,*) '  Damn -- mean IRFs not generated yet.  This might take a while...'
          call flatPrecompute_IRFs(IRF)
        else 
          exit
        endif
      end do

      read(fileunit, '(I10)') Oversampled_Aeff_nEnergies

      allocate(Oversampled_Aeff_logE(Oversampled_Aeff_nEnergies,2))
      allocate(mean_Aeff_lookup(Oversampled_Aeff_nEnergies,2))
      allocate(mean_Aeff_zp(Oversampled_Aeff_nEnergies,2))
      allocate(mean_Aeff_sigma(Oversampled_Aeff_nEnergies,2))
      allocate(temp(2*Oversampled_Aeff_nEnergies-2))
    
      read(fileunit, '(ES13.6)') Oversampled_Aeff_logE
      read(fileunit, '(ES13.6)') mean_Aeff_lookup

      close(fileunit)    

      do strip=front,back,back-front

        !Initialise interpolator
        call TSPSI(Oversampled_Aeff_nEnergies,Oversampled_Aeff_logE(:,strip),mean_Aeff_lookup(:,strip),2,0,.false.,.false., &
             2*Oversampled_Aeff_nEnergies-2,temp,mean_Aeff_zp(:,strip),mean_Aeff_sigma(:,strip),ierr)
        if (ierr .lt. 0) call flatUtils_crash('TSPSI(TSPACK) returned error in mean Aeff initialisation routine.')

      enddo

      deallocate(temp)

      END SUBROUTINE flatIRFini_Aeff_mean
   

      END MODULE flatIRFini
