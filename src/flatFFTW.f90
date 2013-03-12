! FFTW driver module for flatlib
!
! Pat Scott, Feb 2009; pat@fysik.su.se
!--------------------------------------------------------------------------------


      MODULE flatFFTW

      use flatCommon
      use flatIRFs
      use iso_c_binding

      implicit none

      Type(c_ptr), allocatable :: plan_source(:), plan_inverse(:)
      Type(c_ptr) :: plan_PSF
      integer, allocatable :: source_FFTsize_RA(:), source_FFTsize_DEC(:), PSF_FFTsize_RA(:), PSF_FFTsize_DEC(:)
      integer, allocatable :: PSF_pix_RA(:), PSF_pix_DEC(:), padding_RA(:), padding_DEC(:) 
      
      double precision :: TransformedPSFCube_energies(TransformedPSFCube_nEnergies)

      Type(DblPtrArray1D), allocatable :: sourceangles_RA(:), sourceangles_DEC(:)
      Type(DblPtrArray2D), allocatable :: UntransformedSource(:), TransformedInverse(:)
      Type(DblPtrArray3D), allocatable :: TransformedPSFCube(:)
      Type(DblComplexPtrArray2D), allocatable :: TransformedSource(:), UntransformedInverse(:)

      contains
      

      SUBROUTINE flatFFTW_Precompile(nFFTs,RA,DEC,source_FFTsize_RA_in,source_FFTsize_DEC_in,energyRange)
      !Sets up FFTs of PSF.
      ! Performs the FFT of the PSF and puts it into a lookup table in energy for future use in the energy integral.
      ! In principle this could all be easily done totally one-dimensionally since the PSF is circularly symmetric - but
      ! since this need only be done once, and it is easier to follow the code if the full 2D convolution is performed,
      ! the following routine treats everything as if the only symmetries to the PSF are mirror symmetry in the RA
      ! and DEC axes.  Maybe this will also be useful down the track for some instrument with an oblate PSF.

      use flatUtils

      integer, intent(IN) :: nFFTs,source_FFTsize_RA_in(nFFTs),source_FFTsize_DEC_in(nFFTs)
      double precision, intent(IN) :: RA(nFFTs,2), DEC(nFFTs,2), energyRange(2)
      double precision, allocatable :: UntransformedPSF(:,:), TransformedPSF(:,:)
      double precision :: RA_min(nFFTs), RA_max(nFFTs), DEC_min(nFFTs), DEC_max(nFFTs), Direction(2)
      integer :: i, j, k, l
      Type(DblPtrArray1D), allocatable :: PSFangles_RA(:), PSFangles_DEC(:)

      allocate(PSFangles_RA(nFFTs), PSFangles_DEC(nFFTs))
      allocate(sourceangles_RA(nFFTs), sourceangles_DEC(nFFTs))
      allocate(UntransformedSource(nFFTs), TransformedInverse(nFFTs))
      allocate(TransformedSource(nFFTs), UntransformedInverse(nFFTs))
      allocate(TransformedPSFCube(nFFTs), plan_source(nFFTs), plan_inverse(nFFTs))

      allocate(source_FFTsize_RA(nFFTs),source_FFTsize_DEC(nFFTs),PSF_FFTsize_RA(nFFTs),PSF_FFTsize_DEC(nFFTs))
      allocate(PSF_pix_RA(nFFTs),PSF_pix_DEC(nFFTs),padding_RA(nFFTs),padding_DEC(nFFTs))

      RA_min = RA(:,1)*radperdeg; RA_max = RA(:,2)*radperdeg; DEC_min = DEC(:,1)*radperdeg; DEC_max = DEC(:,2)*radperdeg
      source_FFTsize_RA = source_FFTsize_RA_in
      source_FFTsize_DEC = source_FFTsize_DEC_in
      PSF_FFTsize_RA = source_FFTsize_RA/2+1
      PSF_FFTsize_DEC = source_FFTsize_DEC/2+1
      
      do i=1,nFFTs
        allocate(sourceangles_RA(i)%ptr(source_FFTsize_RA(i)))
        allocate(sourceangles_DEC(i)%ptr(source_FFTsize_DEC(i)))
        allocate(PSFangles_RA(i)%ptr(PSF_FFTsize_RA(i)))
        allocate(PSFangles_DEC(i)%ptr(PSF_FFTsize_DEC(i)))
      enddo

      !Initialise regions over which to sample data to be FFTed
      forall (i=1:nFFTs) 
        forall (j=1:PSF_FFTsize_RA(i)) PSFangles_RA(i)%ptr(j) = &
         dble(j-1)/dble(PSF_FFTsize_RA(i)-1)*0.5d0*(RA_max(i)-RA_min(i))
        forall (j=1:PSF_FFTsize_DEC(i)) PSFangles_DEC(i)%ptr(j) = &
         dble(j-1)/dble(PSF_FFTsize_DEC(i)-1)*0.5d0*(DEC_max(i)-DEC_min(i))
        forall (j=1:source_FFTsize_RA(i)) sourceangles_RA(i)%ptr(j) = &
         RA_min(i) + (dble(j)-0.5d0)/dble(source_FFTsize_RA(i)) * (RA_max(i)-RA_min(i))
        forall (j=1:source_FFTsize_DEC(i)) sourceangles_DEC(i)%ptr(j) = &
         DEC_min(i) + (dble(j)-0.5d0)/dble(source_FFTsize_DEC(i)) * (DEC_max(i)-DEC_min(i))
      end forall
     
      !Work out the energies at which FFTs of the PSF will be tabulated
      forall (i=1:TransformedPSFCube_nEnergies) TransformedPSFCube_energies(i) = log10(energyRange(1)) + &
              dble(i-1)/dble(TransformedPSFCube_nEnergies-1) * (log10(energyRange(2)) - log10(energyRange(1)))

      !Work out how far the response needs to be tabulated, and the source padded
      forall (i=1:nFFTs)
        PSF_pix_RA(i) = ceiling(PSF_maxangle/(sourceangles_RA(i)%ptr(2) - sourceangles_RA(i)%ptr(1)))
        PSF_pix_DEC(i) = ceiling(PSF_maxangle/(sourceangles_DEC(i)%ptr(2) - sourceangles_DEC(i)%ptr(1)))
      end forall

      where (PSF_pix_RA > PSF_FFTsize_RA-1) PSF_pix_RA = PSF_FFTsize_RA-1
      where (PSF_pix_DEC > PSF_FFTsize_DEC-1) PSF_pix_DEC = PSF_FFTsize_DEC-1
      where (mod(PSF_pix_RA,2) .eq. 1) PSF_pix_RA = PSF_pix_RA+1
      where (mod(PSF_pix_DEC,2) .eq. 1) PSF_pix_DEC = PSF_pix_DEC+1

      !Add padding if the PSF is not to be periodic
      if (wrapPSFConvolution) then
        padding_RA=0; padding_DEC=0
      else
        padding_RA=PSF_pix_RA; padding_DEC=PSF_pix_DEC
      endif

      !For each ROI requested, do some more allocations and then plan the FFT for each energy
      do i = 1,nFFTs

        write(*,'(A,I2)') '      Doing FFT of PSF for ROI ',i

        allocate(UntransformedPSF(PSF_FFTsize_RA(i) + padding_RA(i)/2, PSF_FFTsize_DEC(i) + padding_DEC(i)/2))
        allocate(TransformedPSF(PSF_FFTsize_RA(i) + padding_RA(i)/2, PSF_FFTsize_DEC(i) + padding_DEC(i)/2))
        allocate(TransformedPSFCube(i)%ptr(PSF_FFTsize_RA(i) + padding_RA(i)/2, source_FFTsize_DEC(i) + padding_DEC(i), &
                                  TransformedPSFCube_nEnergies))

        allocate(UntransformedSource(i)%ptr(source_FFTsize_RA(i) + padding_RA(i),source_FFTsize_DEC(i) + padding_DEC(i)))
        allocate(TransformedSource(i)%ptr((source_FFTsize_RA(i) + padding_RA(i))/2+1,source_FFTsize_DEC(i) + padding_DEC(i)))
        allocate(TransformedInverse(i)%ptr(source_FFTsize_RA(i) + padding_RA(i),source_FFTsize_DEC(i) + padding_DEC(i)))
        allocate(UntransformedInverse(i)%ptr((source_FFTsize_RA(i) + padding_RA(i))/2+1,source_FFTsize_DEC(i) + padding_DEC(i)))

        call dfftw_plan_r2r_2d(plan_PSF, PSF_FFTsize_RA(i) + padding_RA(i)/2, &
                             PSF_FFTsize_DEC(i) + padding_DEC(i)/2, UntransformedPSF, TransformedPSF,&
                             FFTW_REDFT00, FFTW_REDFT00, FFTW_DESTROY_INPUT+FFTW_EXHAUSTIVE)

        do k=1,TransformedPSFCube_nEnergies

          ! Initialise the PSF for this ROI and energy
          do l=1,PSF_pix_RA(i)
            do j=1,PSF_pix_DEC(i)
              Direction = (/PSFangles_RA(i)%ptr(l),PSFangles_DEC(i)%ptr(j)/)
              UntransformedPSF(l,j) = flatIRFs_PSF_mean(Direction,TransformedPSFCube_energies(k), both)
            enddo
          enddo
          UntransformedPSF(PSF_pix_RA(i)+1:,:) = 0.d0
          UntransformedPSF(:,PSF_pix_DEC(i)+1:) = 0.d0

          !Normalise the PSF - requires adding up the data in the four quadrants, even though UntransformedPSF 
          !                    contains only one quadrant plus the two axes (hence avoiding overcounting 1,1 and
          !                    1,: and :,1).
          UntransformedPSF = UntransformedPSF / (sum(UntransformedPSF) + &
                             sum(UntransformedPSF(2:PSF_pix_RA(i),2:PSF_pix_DEC(i))) + &
                             sum(UntransformedPSF(2:PSF_pix_RA(i),:)) + sum(UntransformedPSF(:,2:PSF_pix_DEC(i))))

          !Do the FFT of the (quarter-plane) PSF for this ROI and energy
          call dfftw_execute_r2r(plan_PSF, UntransformedPSF, TransformedPSF)

          !Setup the mirrored FFT, adding the flipside of the DEC dimension
          TransformedPSFCube(i)%ptr(:,:PSF_FFTsize_DEC(i) + padding_DEC(i)/2,k) = TransformedPSF
          TransformedPSFCube(i)%ptr(:,PSF_FFTsize_DEC(i) + padding_DEC(i)/2+1:source_FFTsize_DEC(i)+padding_DEC(i),k) = &
                             TransformedPSF(:,PSF_FFTsize_DEC(i) + padding_DEC(i)/2-1:2:-1)

        enddo

        !Clean up after FFTs of PSF for this ROI
        deallocate(UntransformedPSF, TransformedPSF)
        call dfftw_destroy_plan(plan_PSF)

        !Plan FFTs of source model for this ROI
        call dfftw_plan_dft_r2c_2d(plan_source(i), source_FFTsize_RA(i) + padding_RA(i), source_FFTsize_DEC(i) + padding_DEC(i), &
                               UntransformedSource(i)%ptr, TransformedSource(i)%ptr, FFTW_DESTROY_INPUT+FFTW_EXHAUSTIVE)
        call dfftw_plan_dft_c2r_2d(plan_inverse(i), source_FFTsize_RA(i) + padding_RA(i), source_FFTsize_DEC(i) + padding_DEC(i), &
                               UntransformedInverse(i)%ptr, TransformedInverse(i)%ptr, FFTW_DESTROY_INPUT+FFTW_EXHAUSTIVE)

      enddo

      deallocate(PSFangles_RA, PSFangles_DEC)

      END SUBROUTINE flatFFTW_Precompile


      FUNCTION flatFFTW_TransformedPSF(ROI,E_true)

      integer,intent(IN) :: ROI
      double precision, intent(IN) :: E_true
      double precision :: flatFFTW_TransformedPSF((source_FFTsize_RA(ROI)+padding_RA(ROI))/2+1,source_FFTsize_DEC(ROI)+padding_DEC(ROI))
      double precision :: weightA, weightB, gap
      integer :: lowerindex(1)

      !Algorithm 1 - no interpolation
      !lowerindex = minloc(abs(TransformedPSFCube_energies-E_true))
      !flatFFTW_TransformedPSF = TransformedPSFCube(:,:,lowerindex(1))

      !Algorthim 2 - Parallel F90 array operations
      !lowerindex = minloc(abs(TransformedPSFCube_energies + eoshift(TransformedPSFCube_energies,1) - 2.d0*E_true))
      !gap = TransformedPSFCube_energies(lowerindex(1)+1) - TransformedPSFCube_energies(lowerindex(1)) 
      !weightA = (E_true - TransformedPSFCube_energies(lowerindex(1)))/gap
      !weightB = (TransformedPSFCube_energies(lowerindex(1)+1) - E_true)/gap
      !flatFFTW_TransformedPSF = weightB*TransformedPSFCube(:,:,lowerindex(1)) + weightA*TransformedPSFCube(:,:,lowerindex(1)+1)

      !Algorithm 3 - Old-fashioned F77 search.  Still the fastest apparently...
      if (abs(E_true - TransformedPSFCube_energies(1)) .lt. 1.d-8) then
        lowerindex(1) = 1
      else
        lowerindex(1) = 0
        do while (lowerindex(1) .ne. TransformedPSFCube_nEnergies-1)
          lowerindex(1) = lowerindex(1) + 1
          if (TransformedPSFCube_energies(lowerindex(1)) .le. E_true .and. &
                  TransformedPSFCube_energies(lowerindex(1)+1) .gt. E_true) exit
        enddo
      endif
      gap = TransformedPSFCube_energies(lowerindex(1)+1) - TransformedPSFCube_energies(lowerindex(1)) 
      weightA = (E_true - TransformedPSFCube_energies(lowerindex(1)))/gap
      weightB = (TransformedPSFCube_energies(lowerindex(1)+1) - E_true)/gap
      flatFFTW_TransformedPSF = weightB*TransformedPSFCube(ROI)%ptr(:,:,lowerindex(1)) + weightA*TransformedPSFCube(ROI)%ptr(:,:,lowerindex(1)+1)
      
      END FUNCTION flatFFTW_TransformedPSF


      SUBROUTINE flatFFTW_clean(nFFTs)

      integer, intent(IN) :: nFFTs
      integer :: i

      do i = 1,nFFTs
        call dfftw_destroy_plan(plan_source(i))
        call dfftw_destroy_plan(plan_inverse(i))
        deallocate(UntransformedSource(i)%ptr, TransformedInverse(i)%ptr, TransformedSource(i)%ptr, UntransformedInverse(i)%ptr)
        deallocate(TransformedPSFCube(i)%ptr)
      enddo
      deallocate(UntransformedSource, TransformedInverse, TransformedSource, UntransformedInverse)
      deallocate(source_FFTsize_RA, source_FFTsize_DEC, PSF_FFTsize_RA, PSF_FFTsize_DEC)
      deallocate(PSF_pix_RA, PSF_pix_DEC, padding_RA, padding_DEC, sourceangles_RA, sourceangles_DEC)
      deallocate(TransformedPSFCube, plan_source, plan_inverse)

      END SUBROUTINE flatFFTW_clean


      END MODULE flatFFTW
