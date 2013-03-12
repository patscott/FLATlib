! Test program for flatlib
!
! Pat Scott, Feb 2009; pat@fysik.su.se
!--------------------------------------------------------------------------------



      MODULE testCommon

      use flatCommon
      use flatIRFini
      use flatConvolve_fast
                
      double precision, parameter :: E_min = 0.02d0, E_max = 562.d0
      double precision, parameter :: costheta_min = 0.25d0, costheta_max = 1.d0
      double precision, parameter :: phi_min = 0.d0, phi_max = 2.d0*pi_fl
      integer, parameter :: imax = 1000, jmax = 300, kmax=100
      integer :: i, j, k, counted1, counted2, countrate
      double precision :: convolved_spec(26,26,14), energies(14)

      INTERFACE
        DOUBLE PRECISION FUNCTION sourceFunction(logE_true, Direction_true)
        double precision, intent(IN) :: logE_true, Direction_true(2)
        END FUNCTION sourceFunction
      END INTERFACE
    
      END MODULE testCommon


      PROGRAM flatTest

      use testCommon

      implicit none

      integer :: numROIs = 1
      double precision :: ROIbounds_RA(2) = (/-5.,5./)
      double precision :: ROIbounds_DEC(2) = (/-5.,5./)
      integer :: bins_RA(1) = (/20/)
      integer :: bins_DEC(1) = (/20/)

      call flatConvolve_fast_init('P7CLEAN_V6',numROIs,ROIbounds_RA,ROIbounds_DEC,bins_RA,bins_DEC,.false.)

      call system_clock(counted1,countrate)

      energies = (/3.d-2, 1.d-1, 7.d-1, 12.d-1, 2.d0, 40.d0, 60.d0, 80.d0, 100.d0, 200.d0, 300.d0, 400.d0, 500.d0, 560.d0/)
      do i = 1,14
        write(*,'(A,E15.5)') 'Energy = ', energies(i)
      enddo  

      convolved_spec = flatConvolve_fast_Convolution(numROIs,energies,sourceFunction,pointingType=exposure)

      call system_clock(counted2,countrate)
      write(*,*) 'elapsed time(s): ', dble(counted2 - counted1)/dble(countrate)

      do i = 1,14
        write(*,*)
        write(*,'(A,E15.5)') 'Energy = ', energies(i)
        write(*,'(10E15.5)') convolved_spec(:,10,i)
      enddo

      call flatConvolve_fast_cleanup(1)

      write(*,*) 'FLATtest completed happily.'

      END PROGRAM flatTest

 
      DOUBLE PRECISION FUNCTION sourceFunction(logE_true, Direction_true)

      use flatCommon

      double precision, intent(IN) :: logE_true, Direction_true(2)
    
      double precision :: Direction_true_rad(2)

      Direction_true_rad = Direction_true * radperdeg

      sourceFunction = 1.d0*exp(-0.5*(logE_true-2.d0)**2/30.d0**2 - 0.5*(Direction_true_rad(1) + Direction_true_rad(2))**2/5.d0)

      END FUNCTION sourceFunction


