subroutine flux_calc
  
  USE precisions 
  USE variables
  use grid
  Use constants !,only:zero,pi
  Implicit none
  
  INTEGER(int_p) :: i,j,ifcb,bcon,currf,ic,si,iband,pol, top, left, right, bot, face, face_ad, face_is

  
  REAL(real_p):: ftemp, jin,jout,sdotn,insdotn,actin, vel, mfp, beta, Kn, dist, inot , value
  REAL(real_p):: jeast,jwest,emit, flux_tot, flux_tot_iso_bot, flux_tot_iso_top, flux_tot_iso_right, &
            flux_tot_iso_left, flux_net, flux_tot_isot, tiny1, farea

  character(len=6):: string
!       Allocate(flux_bot(150))
!       Allocate(flux_top(150))
!      Allocate(flux_rigt(150))
!      Allocate(flux_left(150))    

 tiny1= 1.0d-08
flux_tot =zero


  flux_tot_iso_bot = zero
    flux_tot_iso_right = zero
    flux_tot_iso_left = zero
  flux_tot_iso_top = zero
    flux_tot_isot = zero
    face_ad = 0
    face_is =0
    bot = 0
    top =0
    right = 0
    left = 0
    
    
   
    
  Do ifcb = 1, nbcfaces

      currf = bf_to_f(ifcb)
      
      farea = areaf(currf)
      ic = lfc(currf,1)    
      ftemp = temp_bc(ifcb)
      jout = zero
      jin = zero
      emit = zero
      flux_net = zero
      
      Do iband = 1, nbands
         pol = polar(iband)    
        vel = gpvel(iband)
        Call relaxtime(pol, iband, T_ref, beta)
        ! write(10,*) "tau, vel, mfp", one/beta, vel, mfp
        mfp = vel/beta   
        Kn = mfp/tau_l  
               
         if (Kn .ge. cutoff) then    ! balistic part
            
                 call iwalls(pol,iband,ftemp,actin)
                   jout = jout + pi*actin
                   value  = zero
                 Do si = 1,ndir     
                    IF (threed) THEN       
                        sdotn = vecfx(currf) * rmu(si) + vecfy(currf) * rxi(si) + vecfz(currf)*ret(si)        
                        insdotn = vecfx(currf) * inrmu(si) + vecfy(currf) * inrxi(si)+ vecfz(currf)*inret(si)
        !                dist1 = sqrt((xc(i)-xf(currf))**2 + (yc(i)-yf(currf))**2 + (zc(i)-zf(currf))**2)
                     ELSE 
                        sdotn = vecfx(currf) * rmu(si) + vecfy(currf) * rxi(si)     
                        insdotn = vecfx(currf) * inrmu(si) + vecfy(currf) * inrxi(si)
          !              dist1 = sqrt((xc(i)-xf(currf))**2 + (yc(i)-yf(currf))**2) 
                     ENDIF       
                   jin = jin + max(0.0,sdotn)*intensity(ic,si,iband)*(insdotn/sdotn)
                   value = value + max(0.0,sdotn)*intensity(ic,si,iband)*(insdotn/sdotn)
                   emit = emit + max(0.0,-sdotn)*intensity(ic,si,iband)*(insdotn/sdotn)
                ENDDO
                ! modified, sept 1
         
      if ((abs(yf(currf)-1.00E-06)<tiny1) .and. (abs(xf(currf)- 0.88500E-06)<tiny1) ) then !-(1.000000000000000E-006)).LE.(1.00E-25)) then
!        flux_tot_iso_top = flux_tot_iso_top + flux_net
!        top = top + 1
      !  WRITE (io_flux,*) "paree"
!         flux_top(top) = flux_net
!        write(28,*) "flux top" , (pi*actin-value)*farea,  iband
              
      end if
        
        else                !diffuse part
        
                
!!              jout = zero
!!              jin_b = zero
!!              jin_b1 = zero
!!              jin_d = zero
         
         
!!              flux_net = zero
                 IF (threed) THEN
                    dist = sqrt((xc(ic)-xf(currf))**2 + (yc(ic)-yf(currf))**2 + (zc(ic)-zf(currf))**2)
                 ELSE
                    dist = sqrt((xc(ic)-xf(currf))**2 + (yc(ic)-yf(currf))**2)
                 ENDIF
!!              Do iband = 1, nbands
!!                    pol = polar(iband)  
!                    vel = gpvel(iband)
!!                    CALL relaxtime(pol, iband, tnot(ic), beta)
!!                    mfp = vel/beta  
                    call iwalls(pol,iband,temp_bc(ifcb),inot)
                    jout = jout + two*pi*inot
                    jin = jin + phib(ifcb,iband)/two
!!              ENDDO
        !      flux_net = jout - (jin_d)
!!              flux_net = (jout - jin_d) !*farea
!!        !!      flux_net = (jout - jin_d)   
!!              flux_tot_isot = flux_tot_isot+ flux_net
! modified, sept 2
      if ((abs(yf(currf)-1.00E-06)<tiny1) .and. (abs(xf(currf)- 0.88500E-06)<tiny1) ) then !-(1.000000000000000E-006)).LE.(1.00E-25)) then
!        flux_tot_iso_top = flux_tot_iso_top + flux_net
!        top = top + 1
      !  WRITE (io_flux,*) "paree"
!         flux_top(top) = flux_net
!        write(28,*) "flux top" , (two*pi*inot-phib(ifcb,iband)/two)*farea,  iband
              
      end if
                
        
        endif
        
      ENDDO
       flux_net = - (jin-jout)*farea
!       write(28,*) "jout, jin", jout, jin, iband
       face = bf_to_f(ifcb)
       if (abs(xf(face)-1.00E-06)<tiny1) then !-(1.000000000000000E-006)).LE.(1.00E-25)) then
        flux_tot_iso_right = flux_tot_iso_right + flux_net
        right = right + 1
        
         flux_rigt(right) =  flux_net
        write(29,*) "flux rigt" , flux_rigt(right), yf(face)        

      !  WRITE (io_flux,*) "paree"
      end if

       if (abs(yf(face)-1.00E-06)<tiny1) then !-(1.000000000000000E-006)).LE.(1.00E-25)) then
        flux_tot_iso_top = flux_tot_iso_top + flux_net
        top = top + 1
      !  WRITE (io_flux,*) "paree"
         flux_top(top) = flux_net
        write(29,*) "flux top" , flux_top(top), xf(face)
              
      end if

      if (yf(face).EQ.(0.00E+00)) then
        flux_tot_iso_bot = flux_tot_iso_bot + flux_net
        bot = bot + 1
        flux_bot(bot) = flux_net
        write(29,*) "flux bot" , flux_bot(bot), xf(face)
!        flux_bot(bot) = flux_net
!        write(29,*) "flux bot" , flux_bot(bot), xf(face)        
 
      end if

      if (xf(face).EQ.(0.00E+00)) then

        flux_tot_iso_left = flux_tot_iso_left + flux_net
        left = left + 1
        
        flux_left(left) = flux_net
        write(29,*) "flux left" , flux_left(left), yf(face)       
                
      end if
      
      
      
!      flux = flux + (jin-jout)
!      write(io_flux,*) "bcon" ,bctype(ifcb)
!      IF (threed) THEN  
!      WRITE(io_flux,10) xf(currf),yf(currf),zf(currf),flux
!      ELSE
!      WRITE(io_flux,10) xf(currf),yf(currf),flux
!      ENDIF
!      
      select case(bctype(ifcb))
      case(ADIA)
!      flux_tot = flux_tot +flux
      end select
      
        
  ENDDO
!    write(io_flux,*) "flux_tot"
!   write(io_flux,10) flux_tot

   WRITE(io_flux,*)"flux_tot_iso_right, flux_tot_iso_left, flux_tot_iso_top, flux_tot_iso_bot, jout "
   WRITE(io_flux,10) flux_tot_iso_right/right, flux_tot_iso_left/left, flux_tot_iso_top/top, &
                        flux_tot_iso_bot/bot, jout !"flux_tot, t", flux_tot/nbcfaces, t
!   WRITE(io_flux,*) tindex
   WRITE(io_flux,*)" "
  
    WRITE(io_flux,*)"t,b,l,r, farea ", top, bot, left, right, farea
   WRITE(io_flux,*)" "
   
! WRITE(io_flux,*) " total flux is", 

10 format(4(1x,E14.6))

    close (io_output)
!    close(io_flux)
  
 ! IF(debug_level > 0) WRITE(io_debug,*) "Finishing postprocess" 
   end subroutine flux_calc

