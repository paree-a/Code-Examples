SUBROUTINE knudsen_plot

  USE PRECISIONS
  USE VARIABLES
  USE GRID
  USE CONSTANTS
  
  
  REAL(real_p) ::  vel, beta, kn, beta2
  REAL(real_p) :: tnot1, tnot2
  INTEGER(int_p) :: iband, pol
 
! OPEN(UNIT = 28 , FILE = "calculations.txt" ,STATUS = "REPLACE" )
 
!        DO i = 1,ncells  
!        write (28,*) xc(i), yc(i)
! enddo
         
         tnot1 = 200.0d0
         
         
  !Do iband = 2,2  
  Do iband = 1, nbands    !!$ loop over frequency bands
        pol = polar(iband)
        vel = gpvel(iband)

        CALL relaxtime(pol, iband, tnot1, beta)
       write (28,*) "tau vs fre", one/beta, iband
!        kn = (vel/beta)/(1.0d-06)
!        write (28,*) "kn vs  fre", kn, iband  
!         tnot2 = 100.0d0
!        CALL relaxtime(pol, iband, tnot2, beta2)
!        write (28,*) "tau vs fre", one/beta2, iband
       ! write (28,*) "vel vs fre", vel, iband    
  end do
  
    Do iband = 1, nbands    !!$ loop over frequency bands
        pol = polar(iband)
        vel = gpvel(iband)
!        tnot1 = 250.0d0
        CALL relaxtime(pol, iband, tnot1, beta)
!       write (28,*) "tau vs fre", one/beta, iband
        kn = (vel/beta)/(1.0d-06)
        write (28,*) "kn vs  fre", kn, iband  
!         tnot2 = 100.0d0
!        CALL relaxtime(pol, iband, tnot2, beta2)
!        write (28,*) "tau vs fre", one/beta2, iband
       ! write (28,*) "vel vs fre", vel, iband    
  end do
  
  
!     close (28)
  
! stop
      
 end SUBROUTINE knudsen_plot      