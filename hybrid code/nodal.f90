SUBROUTINE nodintpln

  USE PRECISIONS
  USE CONSTANTS !, only : zero,two
  USE VARIABLES !, only : phiv,phi,phib
  USE GRID !,only : ncells,nbcfaces,lcv,lfv,nfcell,bf_to_f,nnodes,bnode,nvcell,nvface,nodeweight,xc,xv,xf,yc,yv,yf,zc,zf,zv !wcv,wbfv

  IMPLICIT NONE

  INTEGER(int_p) :: i,j,currv,currf,iband
  REAL(real_p) :: disntemp
  !modified the weight part
  REAL(real_p), DIMENSION (:), ALLOCATABLE :: weight
  
  ALLOCATE (weight(nnodes))
  
  IF(debug_level > 0) WRITE(io_debug,*) "Starting nodal/nodintpln" 

  phiv(:,:) = zero
  weight(:) = zero
  IF (threed) THEN 

  do iband = 1,nbands
  
        DO i = 1,ncells
           DO j = 1,nvcell(i)
              currv = lcv(i,j)
              IF( bnode(currv) == 0 ) THEN

                  phiv(currv,iband) = phiv(currv,iband) + phi(i,iband)*wcv(i,j)

              END IF
           END DO
        END DO
   End do
  ELSE 

  do iband = 1,nbands
        DO i = 1,ncells
          DO j = 1,nvcell(i)
              currv = lcv(i,j)
              IF( bnode(currv) == 0 ) THEN

                  phiv(currv,iband) = phiv(currv,iband) + phi(i,iband)*wcv(i,j)

              END IF
           END DO
        END DO
   end do
  ENDIF

  do iband = 1,nbands
    
    DO i = 1,nnodes
       IF(bnode(i) == 1) THEN
            phiv(i,iband) = 0
       END IF
     ENDDO
 
 End do


  IF (threed) THEN  

  do iband = 1,nbands
         DO i = 1,nbcfaces
            currf = bf_to_f(i)
            DO j = 1,nvface(currf)
                currv = lfv(currf,j)

                weight(currv) = weight(currv) + wcv(i,j) 
                phiv(currv,iband) = phiv(currv,iband) + phib(i,iband)*wbfv(i,j)                

            END DO
         END DO
    End do
   ELSE

  do iband = 1,nbands
         DO i = 1,nbcfaces
            currf = bf_to_f(i)
            DO j = 1,nvface(currf)
                currv = lfv(currf,j)

                 weight(currv) = weight(currv) + wcv(i,j) 
                 phiv(currv,iband) = phiv(currv,iband) + phib(i,iband)*wbfv(i,j)

            END DO
         END DO
    end DO
   ENDIF
  

  
  DEALLOCATE(weight)
  
  IF(debug_level > 0) WRITE(io_debug,*) "Exiting nodal/nodintpln" 

END SUBROUTINE nodintpln



