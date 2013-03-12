! Standalone program to precompile averaged IRFs for flatlib
!
! Pat Scott, Mar 2009; pat@fysik.su.se
!--------------------------------------------------------------------------------


      PROGRAM flatAverage

      use flatIRFini

      implicit none
      character(len=strlen) :: IRF     

      if (iargc() < 1) call flatUtils_crash('Usage: flataverage IRF_set_name')
      call getarg(1,IRF)
      if (IRF(len_trim(IRF):len_trim(IRF)) .eq. '/') IRF = IRF(1:len_trim(IRF)-1)
      write(*,*)
      call flatPrecompute_IRFs(IRF)

      END PROGRAM flatAverage
