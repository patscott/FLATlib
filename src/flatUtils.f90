! Utilities module for flatlib
!
! Pat Scott, Feb 2009; pat@fysik.su.se
!--------------------------------------------------------------------------------


      MODULE flatUtils

      use flatCommon

      implicit none

      contains
      

      SUBROUTINE flatUtils_crash(msg)

      character (len=*), intent(IN) :: msg

      write(*,*)
      write(*,*) msg
      write(*,*)

      stop

      END SUBROUTINE flatUtils_crash

     
      FUNCTION flatUtils_GaussianMesh(logE_obs_local, energyRange)

      double precision, intent(IN) :: logE_obs_local(:), energyRange(2)
      double precision, allocatable :: flatUtils_GaussianMesh(:,:), upperRange(:), lowerRange(:)
      double precision, allocatable :: totalUpperInt(:), totalLowerInt(:)
      double precision :: erfinv, inner
      integer, allocatable :: nUpperSamples(:), nLowerSamples(:)
      integer :: i, e, n_Energies
      external erfinv

      n_Energies = size(logE_obs_local)
      if (.not. allocated(flatUtils_GaussianMesh)) allocate(flatUtils_GaussianMesh(n_Energies, nFastSamples))
      allocate(upperRange(n_Energies), lowerRange(n_Energies))
      allocate(totalUpperInt(n_Energies), totalLowerInt(n_Energies), nUpperSamples(n_Energies), nLowerSamples(n_Energies))

      upperRange = max(log10(Fermi_Emax) - logE_obs_local,0.d0)
      lowerRange = max(logE_obs_local - log10(Fermi_Emin),0.d0)

      totalUpperInt = erf(upperRange/sigmaFastSamples/sqrt(2.d0))
      totalLowerInt = erf(lowerRange/sigmaFastSamples/sqrt(2.d0))
      nUpperSamples = floor(dble(nFastSamples) * totalUpperInt / (totalUpperInt + totalLowerInt))
      nLowerSamples = nFastSamples - nUpperSamples - 1

      do e=1,n_Energies

        flatUtils_GaussianMesh(e,1) = log10(energyRange(1))
        flatUtils_GaussianMesh(e,nLowerSamples+1) = logE_obs_local(e)
        flatUtils_GaussianMesh(e,nFastSamples) = log10(energyRange(2))

        do i=1,nUpperSamples(e)-1
          inner = dble(i)/dble(nUpperSamples(e)) * totalUpperInt(e)
          flatUtils_GaussianMesh(e,nLowerSamples(e)+1+i)= logE_obs_local(e) + erfinv(inner, 1.d0-inner) * sigmaFastSamples * sqrt(2.d0)
        enddo

        do i=1,nLowerSamples(e)-1
          inner = dble(i)/dble(nLowerSamples(e)) * totalLowerInt(e)
          flatUtils_GaussianMesh(e,nLowerSamples(e)+1-i)= logE_obs_local(e) - erfinv(inner, 1.d0-inner) * sigmaFastSamples * sqrt(2.d0)
        enddo

      enddo

      END FUNCTION flatUtils_GaussianMesh


      END MODULE flatUtils
