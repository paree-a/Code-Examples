SUBROUTINE update

  USE precisions !,ONLY:int_p
  USE grid !,ONLY: ncells,nbcfaces
  USE Variables !,ONLY:phi,phione,phitwo,gw,gwone,tnot,tnotone,phib,phibone,ndir,intensity,intensityone
  IMPLICIT NONE
  INTEGER(int_p)::i,si,ic,iband
  real(real_p):: tiny2
    tiny2 = 1.0d-09
  IF(debug_level > 0) WRITE(io_debug,*) "Starting Update" 

   
   

  

  
do ic = 1, ncells 

     do iband = 1,nbands
        gwone(ic,iband) = gw(ic,iband)         
        phitwo(ic,iband) = phione(ic,iband)
        phione(ic,iband) = phi(ic,iband)  
     enddo
     
     tnotone(ic) = tnot(ic)
       
  end do
  
  
 do iband = 1,nbands
      DO i = 1,nbcfaces    
         phibone(i,iband) = phib(i,iband)
      ENDDO
 end do
  
  !           do iband = 2,2
             ! do iband = 14,20
  do iband = 1,nbands
     do ic = 1,ncells
        do si = 1,ndir
           intensityone(ic,si,iband) = intensity(ic,si,iband)            
        enddo
     enddo
  enddo


IF(debug_level > 0) WRITE(io_debug,*) "Exiting Update" 
END SUBROUTINE update
