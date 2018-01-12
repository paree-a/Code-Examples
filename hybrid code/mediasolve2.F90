SUBROUTINE mediasolve2(pol,band,vel)

! built using the regular boundary conditions

  USE PRECISIONS
  USE VARIABLES
  USE GRID
  USE CONSTANTS

  IMPLICIT NONE

  REAL(real_p) :: de,tol_rad,alpha,tempmat(7),mat_old, local_tvel,beta1, beta2, beta3, beta, beta_old
  REAL(real_p) :: con2,gammabf,vol,farea,idtsq,netflux, guess, told, inot, inotone,iprime, phi_sum, disn1, iwallf
  INTEGER(int_p) :: i,j,currf,currcell,currbf,k,istart,ic1,ic2,  face_gnum
  INTEGER(int_p) ::tempindex(7),l,x, index_old,kspace,ierr 
  
  Integer(int_p),intent(in):: band, pol
  REAL(real_p), intent(in) :: vel
  
  real(real_p):: tiny2
  tiny2 = 1.0d-08
  
  IF(debug_level > 0) WRITE(io_debug,*) "Starting mediasolve" 
 
  idtsq = idt**2 
  CALL nodintpln     
  CALL fluxt  
 
  Tol_rad = 1.0D-8    
  
        aps(:) = zero
        anbs(:,:) = zero
        scs(:) = zero
        resid(:) = zero
        r(:) = zero
        soln(:) = zero
        
        DO i = 1,ncells
           vol = volcell(i)
   
           CALL relaxtime(pol, band, tnot(i), beta) 

           con2 = (vel)*third/beta


           DO j = 1,nfcell(i)
              currf = lcf(i,j)                          

              IF ( bface(currf) == 0 ) THEN
                 ic1 = lfc(currf,1)
                 ic2 = lfc(currf,2)

                 CALL relaxtime(pol, band, tnot(ic1), beta1) 
                 CALL relaxtime(pol, band, tnot(ic2), beta2) 
                 local_tvel = (vel*wcf(currf)/beta2 + vel*(1-wcf(currf))/beta1) !/1.0d-06
                 
                 de =  areaf(currf)*con2 /disn(currf)

                 aps(i) = aps(i) + local_tvel*de 
                 anbs(i,j) = -local_tvel*de 
                 scs(i) = scs(i) - normdir(i,j) * jt(currf,band) * de  * local_tvel
              ELSE
                 local_tvel = vel/beta 

                 de =  areaf(currf)*con2 /disn(currf) 
                 aps(i) = aps(i) + de*local_tvel
                 anbs(i,j) = -de*local_tvel
                 scs(i) = scs(i) - jt(currf,band) * de * local_tvel 

              END IF
           END DO     !ending nfcells
           
           IF (jhanda)then  
                    
               guess = tnot(i)
               told = tnotone(i)   
               CALL iwalls(pol, band, guess,inot) 
               CALL relaxtime(pol, band, tnotone(i), beta_old)    

               aps(i) = aps(i)+ (one + two*idtsq/(beta**2) + two*idt/beta  + idt*(1/beta - 1/beta_old)*idt/beta )*vol    
               CALL iwalls(pol, band, told,inotone)  
               netflux = 4*pi*inot + 4*pi*idt*(inot - inotone)/beta
            
               scs(i) = scs(i) + netflux*vol + phione(i,band)*(two*idtsq/(beta**2) +two*idt/beta)*vol + & 
                        zero*four*pi*inotone*idt*vol/beta  + &  ! made dg/dt = 0
                        idt*(1/beta - 1/beta_old)*phione(i,band)*idt*vol/beta ! extra term


           ELSE
         
             guess = tnot(i)
             told = tnotone(i)
            CALL relaxtime(pol, band, tnotone(i), beta_old)   
               CALL iwalls(pol, band, guess,inot)  
               CALL iwalls(pol, band, told,inotone)  
               CALL ipcalc(pol, band, guess,iprime) 
               netflux = 4*pi*inot + 4*pi*idt*(inot - inotone)/beta
               aps(i) = aps(i)+ (one + idtsq/(beta**2) + two*idt/beta + idt*(1/beta - 1/beta_old)*idt/beta)*vol
               scs(i) = scs(i) + netflux*vol + vol*phione(i,band)*two*(idtsq/(beta**2) + idt/beta) & 
                        - idtsq*phitwo(i,band)*vol/(beta**2) +  idt*(1/beta - 1/beta_old)*phione(i,band)*idt*vol/beta
         
           ENDIF
     
           DO j = 1,nfcell(i)
               currf = lcf(i,j)
               IF ( bface(currf) == 0 ) THEN
                   IF( lfc(currf,1) == i ) THEN
                       currcell = lfc(currf,2)
                   ELSE
                       currcell = lfc(currf,1)
                   END IF
                   resid(i) = resid(i) - anbs(i,j) * phi(currcell,band)
               ELSE
                   currbf = f_to_bf(currf)
                   resid(i) = resid(i) - anbs(i,j) * phib(currbf,band)
               END IF
           END DO
           resid(i) = resid(i) - aps(i) * phi(i,band)
        
        END DO  ! ended ncells
        
        DO i = 1,nbcfaces
            currf = bf_to_f(i)
            currcell = bf_to_c(i)
            farea = areaf(currf)

           CALL relaxtime(pol, band, temp_bc(currf), beta)
           CALL iwalls(pol, band, temp_bc(currf),iwallf) 
            gammabf = two*vel*third*farea/disn(currf)/beta !/1.0d-06
            aps(i+ncells) = (one+ idt/beta)*farea + gammabf          
            anbs(i+ncells,1) = -gammabf
 
            scs(i+ncells) = gammabf*jt(currf,band) + phibone(i,band)*idt*farea/beta  + four*pi*iwallf*farea                    
            resid(i+ncells) = resid(i+ncells) - anbs(i+ncells,1) * phi(currcell,band) - & 
                              aps(i+ncells) * phib(i,band)  

        ENDDO
        DO i = 1, (ncells+nbcfaces)

        
                      anbs(i,:) = anbs(i,:) / aps(i)
              scs(i) = scs(i)/aps(i)
              resid(i) = resid(i) / aps(i)
              aps(i) = one
           resid(i) = resid(i) + scs(i)
        ENDDO
        r(:) = resid(:)
        am(:) = zero
        ia(:) = zero
        ja(:) = zero
        k = 0
        istart = 1
        DO i = 1,ncells
           ia(istart) = k + 1  
           tempmat = zero
           tempindex = zero
           tempmat(1) = aps(i)
           tempindex(1) = i 
           l = 1
           DO j = 1,nfcell(i)
              tempmat(j+1) = anbs(i,j)
              currf = lcf(i,j)
               IF( bface(currf) == 0 ) THEN
                   IF( lfc(currf,1) == i ) THEN
                       currcell = lfc(currf,2)
                   ELSE
                       currcell = lfc(currf,1)
                   END IF
                   tempindex(j+1) = currcell
               ELSE
                   tempindex(j+1) = ncells + f_to_bf(currf)
               END IF
               l = l + 1
           END DO
           DO j = 1,l-1
              DO x = j+1,l
                 IF( tempindex(x) < tempindex(j) ) THEN
                    index_old = tempindex(x)
                    mat_old = tempmat(x)
                    tempmat(x) = tempmat(j)
                    tempindex(x) = tempindex(j)
                    tempmat(j) = mat_old
                    tempindex(j) = index_old
                 END IF
              END DO
           END DO

           am(k+1:k+l) = tempmat(:)
           ja(k+1:k+l) = tempindex(:)
           k = k + l
           istart = istart + 1
        END DO
        DO i = 1,nbcfaces
           IA(istart) = k + 1
           k = k + 1
           AM(k) = anbs(i+ncells,1)
           JA(k) = bf_to_c(i)
           k = k + 1
           AM(k) = aps(i+ncells)
           JA(k) = i + ncells
           istart = istart + 1
        END DO

        IA(istart) = k + 1
        alu(:) = zero
        jlu(:) = zero
        ju(:) = zero
        jr(:) = zero
        vv(:,:) = zero
        KSPACE = 20

        CALL ilu0(ncells+nbcfaces,AM,JA,IA,alu,jlu,ju,jr,ierr)

        CALL pgmres(ncells+nbcfaces,kspace,r,soln,vv,tol_rad,4000,0,AM,JA,IA,alu,jlu,ju,ierr)

        DO i = 1,ncells

           phi(i,band) = phi(i,band) + soln(i) 

        END DO

        DO i = 1,nbcfaces

           phib(i,band) = phib(i,band) + soln(i+ncells) ! check this term   
         
        END DO

        
  
 
 
  
  IF(debug_level > 0) WRITE(io_debug,*) "Exiting mediasolve" 
  
  END SUBROUTINE mediasolve2
    






