Subroutine tempcalc

 

  Use precisions
  use constants
  use grid
  use variables

!  use constants,only:zero,four,sigma,one,three,half,two
!  use grid,only:ncells,nbcfaces,bctype,bf_to_f,lfc,vecfx,vecfy,temp_bc,xc,xf,yc,yf,zf,zc,vecfz
!  use variables,only:gw,phi,tnot,tprime,residual,tnotone,gwone,ndir,intensity,mfp,rxi,rmu,inrxi,inrmu,phib,tbprime,ret,inret,idt,phione


  Implicit none   
  integer(int_p)::i,ifcb,currf,bcon,ic,si
  real(real_p):: sdotn,insdotn,jtotal,divb,divm,dist,mfp

  Integer(int_p)::iband,pol, factor
  !added
  Real(real_p):: guess,inot,iprime,corr,imfp,ilp,jout_totbands,diffI_temp,I_not_old
  Real(real_p):: vel,beta,pbeta,inotbeta,iprimebeta,gnbeta,gnb,gna,g_transient1,g_transient2,g_by_tau
  Real(real_p):: told,uold,unew,uprime,inotone, tempsum1, tempsum2, divergence
  REAL(real_p), ALLOCATABLE, DIMENSION(:):: jout,phib_total
  IF(debug_level > 0) WRITE(io_debug,*) "Starting tempcalc" 
  
  ALLOCATE (jout(nbands))
  ALLOCATE (phib_total(nbcfaces))
  residual = zero  
  factor = 1
  
  
  
  !************modify
  
    
  Do ic = 1,ncells
        guess = tnot(ic)
        told  = tnotone(ic)     
        gna = zero     
        gnb = zero     
        uold = zero  
        unew = zero
        uprime = zero
        tempsum1 = zero
        tempsum2 = zero  
        g_transient1 =zero
        g_transient2 = zero 
        g_by_tau =zero 
        divergence = zero

        do iband = 1,nbands
            pol = polar(iband)
            vel = gpvel(iband)   
            call relaxtime(pol, iband, guess, beta)               
            call iwalls(pol, iband, guess,inot)  
            call iwalls(pol, iband, told,inotone)  
            call ipcalc(pol, iband, guess,iprime)     
            imfp = beta/vel
            ilp  = idt/vel
            gna = gna + gw(ic,iband)*(imfp + ilp) + phi(ic,iband)*(imfp + ilp)
            
            gnb = gnb + gwone(ic,iband)*ilp + phione(ic,iband)*ilp 
            uold = uold + four*pi*ilp*inotone
            unew = unew + four*pi*(imfp+ilp)*inot
            uprime = uprime + four*pi*(imfp+ilp)*iprime   
          
            tempsum1 = tempsum1 + four*pi*inot*imfp
            g_transient1 = g_transient1 +  phi(ic,iband)*ilp + gw(ic,iband)*ilp 
            g_transient2 = g_transient2 +  phione(ic,iband)*ilp  + gwone(ic,iband)*ilp 
            g_by_tau = g_by_tau +  gw(ic,iband)*imfp +phi(ic,iband)*imfp
            
            
        enddo

         tempsum2 = g_transient1 - g_transient2 !+ g_by_tau ! - tempsum1
         corr = (uold + gna - gnb - unew)/uprime  
         divergence = divergence + gna - gnb - tempsum1
         residual  = residual + corr**2      
         tprime(ic) = corr   

  Enddo
  


 phib_total(:) = 0

  Do ifcb = 1, nbcfaces    
         bcon = bctype(ifcb)
         currf = bf_to_f(ifcb)
         ic = lfc(currf,1)
         IF (threed) THEN
            dist = sqrt((xc(ic)-xf(currf))**2 + (yc(ic)-yf(currf))**2 + (zc(ic)-zf(currf))**2)
         ELSE
            dist = sqrt((xc(ic)-xf(currf))**2 + (yc(ic)-yf(currf))**2)
         ENDIF
         Select Case (bcon)
         Case(4)
             jout(:) = zero
             jout_totbands = zero
             diffI_temp = zero
             I_not_old = zero
             jtotal = zero
 !           do iband = 2,2
             ! do iband = 14,20
             do iband = 1,nbands
                 pol = polar(iband)
                 vel = gpvel(iband)
                 IF (threed) THEN
                     Do si = 1,ndir            
                        sdotn = vecfx(currf) * rmu(si) + vecfy(currf) * rxi(si) +  vecfz(currf) * ret(si)
                        insdotn = vecfx(currf) * inrmu(si) + vecfy(currf) * inrxi(si) + vecfz(currf) * inret(si)
                        jout(iband) = jout(iband) + max(0.0d0,sdotn) *intensity(ic,si,iband)*(insdotn/sdotn)
                        enddo
                 ELSE
                     Do si = 1,ndir            
                        sdotn = vecfx(currf) * rmu(si) + vecfy(currf) * rxi(si)
                        insdotn = vecfx(currf) * inrmu(si) + vecfy(currf) * inrxi(si) 
                        jout(iband) = jout(iband) + max(0.0d0,sdotn) *intensity(ic,si,iband)*(insdotn/sdotn)
                     enddo
                 ENDIF    
                 CALL relaxtime(pol, iband, tnot(ic), beta)
!                 CALL relaxtime(pol, iband, temp_bc(ifcb), beta)
                 mfp = vel/beta !/1.0d-06
                 jout_totbands = jout_totbands + jout(iband)*exp(-dist*factor/mfp)
                 phib_total(ifcb) = phib_total(ifcb) + phib(ifcb,iband)  
  !               CALL iwalls(pol, iband, tnotone(ic), inot) 
                 CALL iwalls(pol, iband, temp_bc(ifcb), inot)  
                 I_not_old = I_not_old + inot   
                 CALL ipcalc(pol, iband, temp_bc(ifcb),iprime) 
                 diffI_temp = diffI_temp + iprime
             end do ! end bands 
             jtotal = jout_totbands + half*phib_total(ifcb)
             tbprime(ifcb) = (jtotal/pi - I_not_old )/diffI_temp
!            tbprime(ifcb) = (jtotal/pi - I_not_old )/(four*0.2826d0*temp_bc(ifcb)**3)
             !write(28,*)" tbprime(ifcb)", tbprime(ifcb)     
             residual  = residual + tbprime(ifcb)**2
        Case(3)
             tbprime(ifcb) = zero 
        end select
  enddo 

residual = sqrt(residual)
!write(28,*) "residual", residual
!write(*,*) "residual", residual

IF(debug_level > 0) WRITE(io_debug,*) "Exiting tempcalc" 

end subroutine tempcalc







        
        
        
        
