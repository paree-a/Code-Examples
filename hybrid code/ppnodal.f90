SUBROUTINE ppnode

  USE PRECISIONS
  USE CONSTANTS !, only : zero,two
  USE VARIABLES !, only : tnot,tnotv
  USE GRID !,only : ncells,nbcfaces,lcv,lfv,nfcell,bf_to_f,nnodes,bnode,nvcell,nvface,nodeweight,xc,xv,xf,yc,yv,yf,zc,zf,zv,temp_bc

  IMPLICIT NONE

  INTEGER(int_p) :: i,j,currv,currf
  REAL(real_p) :: disntemp
  REAL(real_p), DIMENSION (:), ALLOCATABLE :: weight

  ALLOCATE (weight(nnodes))
  
  IF(debug_level > 0) WRITE(io_debug,*) "Starting ppnodals" 
  
  weight(:) = zero

  tnotv(:) = zero
!  IF (threed) THEN          !should work for either case based on wcv and wbfv definitions
      DO i = 1,ncells
         DO j = 1,nvcell(i)
            currv = lcv(i,j)
            IF( bnode(currv) == 0 ) THEN       
               tnotv(currv) = tnotv(currv) + tnot(i)*wcv(i,j);
            END IF
         END DO
      END DO
      
    DO i = 1,nnodes
       IF(bnode(i) == 1) THEN
            tnotv(i) = 0
       END IF
    ENDDO

     
      DO i = 1,nbcfaces
         currf = bf_to_f(i)
         DO j = 1,nvface(currf)
         

            currv = lfv(currf,j)

            weight(currv) = weight(currv) + wbfv(i,j) 
            tnotv(currv) = tnotv(currv) + temp_bc(i)*wbfv(i,j)
         END DO
      END DO
     
  
  
  DO i = 1,nnodes
     IF(bnode(i) /= 1)CYCLE 
     tnotv(i) = tnotv(i)/weight(i)
  END DO

  DEALLOCATE(weight)

  IF(debug_level > 0) WRITE(io_debug,*) "Exiting ppnodals" 

END SUBROUTINE ppnode