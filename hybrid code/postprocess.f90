  
  Subroutine post_process(t)

  Use precisions !,only:int_p
  Use grid !,only:nnodes,ncells,xv,yv,zv,lcv,nvcell
  use variables !,only:tnotv,kn,dt,total_time
  USE constants
  Implicit none 

  integer(int_p),intent(in) :: t
  integer(int_p):: j,i, tindex
  
  REAL(real_p):: ftemp, jin_b,jin_b1, jin_d,jout,sdotn,insdotn,inot,&
                 flux_net, dist, vel, beta, mfp, flux_tot_adia, flux_tot_isot, flux_tot_adia_bot, flux_tot_adia_top, &
                  dr, act_hot, act_cold
    INTEGER(int_p) :: ifcb,bcon,currf,ic,si,iband,pol, face_ad, face_is, bot, top, face, right , left
    
  character(len=1024) :: filename1, filename2
  character(len=1024) :: format_string1, format_string2
  
  
  IF(debug_level > 0) WRITE(io_debug,*) "Starting post_process" 
  
  CALL ppnode


  OPEN( UNIT = io_output , FILE = "contour.dat" , STATUS = "REPLACE" )
 
  OPEN( UNIT = io_info , FILE = "info.dat" , STATUS = "REPLACE" )
 
!   from here


        if (t < 10) then
            format_string1 = "(A7,I1,A4)"
            format_string2 = "(A4,I1,A4)"
        else if (t < 100) then
            format_string1 = "(A7,I2,A4)"
            format_string2 = "(A4,I2,A4)"
        else if (t < 1000) then
            format_string1 = "(A7,I3,A4)"
            format_string2 = "(A4,I3,A4)"
        else 
            format_string1 = "(A7,I4,A4)"
            format_string2 = "(A4,I4,A4)"
        endif

        write (filename1,format_string1) "contour", t, ".dat"
        write (filename2,format_string2) "info", t, ".dat"
        print *, trim(filename1)
        print *, trim(filename2)
  
  OPEN( UNIT = io_output , FILE = filename1, STATUS = "unknown" )
  OPEN( UNIT = io_info , FILE = "info.dat" , STATUS = "unknown" )
 
 
   ! end set 4, 2015
 
 
!   
!  OPEN( UNIT = io_output , FILE = filename1, STATUS = "REPLACE" )
!  OPEN( UNIT = io_info , FILE = filename2 , STATUS = "REPLACE" )
  IF (threed) THEN
        WRITE(io_output,*) 'VARIABLES = "X","Y","Z","T"'
        WRITE(io_output,*) 'ZONE T = "BDE" , N = ',nnodes,' , E = ',ncells
        WRITE(io_output,*) 'ZONETYPE = FEBRICK, DATAPACKING = BLOCK '
        WRITE(io_output,*) 'VARLOCATION=(4=NODAL)'
        WRITE(io_output,*) (xv(j), j = 1,nnodes)
        WRITE(io_output,*) (yv(j), j = 1,nnodes)
        WRITE(io_output,*) (zv(j), j = 1,nnodes)
        WRITE(io_output,*) (tnotv(j), j = 1,nnodes)
       
        DO i = 1,ncells
            WRITE(io_output,*) ( lcv(i,j) , j = 1,nvcell(i) )
        END DO
       
  ELSE
    
        WRITE(io_output,*) 'VARIABLES = "X","Y","T"'
        WRITE(io_output,*) 'ZONE T = "BDE" , N = ',nnodes,' , E = ',ncells
        WRITE(io_output,*) 'ZONETYPE = FEQUADRILATERAL, DATAPACKING = BLOCK '
        WRITE(io_output,*) 'VARLOCATION=(3=NODAL)'
        WRITE(io_output,*) (xv(j), j = 1,nnodes)
        WRITE(io_output,*) (yv(j), j = 1,nnodes)
        WRITE(io_output,*) (tnotv(j), j = 1,nnodes)
   
        DO i = 1,ncells
            WRITE(io_output,*) ( lcv(i,j) , j = 1,nvcell(i) )
        END DO
   
   ENDIF
 
!   Write(io_info,*), "Knudsen number is = " ,vel/beta/1.0d-06
   Write(io_info,*)," time is = ",t*dt
   Write(io_info,*),"time step =",  dt
   Write(io_info,*),"time index is =",t
   Write(io_info,*),"time taken per step is =", total_time/t
   CLOSE( UNIT = io_output )
   CLOSE( UNIT = io_info )


! till here




!  OPEN( UNIT = io_output , FILE = filename1, STATUS = "REPLACE" )
!  OPEN( UNIT = io_info , FILE = filename2 , STATUS = "REPLACE" )
!!  IF (threed) THEN
!!        WRITE(io_output,*) 'VARIABLES = "X","Y","Z","T"'
!!        WRITE(io_output,*) 'ZONE T = "BDE" , N = ',nnodes,' , E = ',ncells
!!        WRITE(io_output,*) 'ZONETYPE = FEBRICK, DATAPACKING = BLOCK '
!!        WRITE(io_output,*) 'VARLOCATION=(4=NODAL)'
!!        WRITE(io_output,*) (xv(j), j = 1,nnodes)
!!        WRITE(io_output,*) (yv(j), j = 1,nnodes)
!!        WRITE(io_output,*) (zv(j), j = 1,nnodes)
!!        WRITE(io_output,*) (tnotv(j), j = 1,nnodes)
!!       
!!        DO i = 1,ncells
!!            WRITE(io_output,*) ( lcv(i,j) , j = 1,nvcell(i) )
!!        END DO
!!       
!!  ELSE
!!    
!!        WRITE(io_output,*) 'VARIABLES = "X","Y","T"'
!!        WRITE(io_output,*) 'ZONE T = "BDE" , N = ',nnodes,' , E = ',ncells
!!        WRITE(io_output,*) 'ZONETYPE = FEQUADRILATERAL, DATAPACKING = BLOCK '
!!        WRITE(io_output,*) 'VARLOCATION=(3=NODAL)'
!!        WRITE(io_output,*) (xv(j), j = 1,nnodes)
!!        WRITE(io_output,*) (yv(j), j = 1,nnodes)
!!        WRITE(io_output,*) (tnotv(j), j = 1,nnodes)
!!   
!!        DO i = 1,ncells
!!            WRITE(io_output,*) ( lcv(i,j) , j = 1,nvcell(i) )
!!        END DO
!!   
!!   ENDIF
!! 
!!!   Write(io_info,*), "Knudsen number is = " ,vel/beta/1.0d-06
!!   Write(io_info,*)," time is = ",t*dt
!!   Write(io_info,*),"time step =",  dt
!!   Write(io_info,*),"time index is =",t
!!   Write(io_info,*),"time taken per step is =", total_time/t
!!   CLOSE( UNIT = io_output )
!!   CLOSE( UNIT = io_info )


dr =zero

     Do iband = 1, nbands
            pol = polar(iband)      
            call iwalls(pol,iband,255.0d0,act_hot)
            call iwalls(pol,iband,245.0d0,act_cold)
            dr = dr + pi*(act_hot - act_cold)*1e-08
        
!            write(*,*) "dr" , dr
      enddo
          write(io_flux,*) "flux dr", dr

!end if

! OPEN( UNIT = 45 , FILE = "contour.dat" , STATUS = "REPLACE" )
!  WRITE(45,*) 'VARIABLES = "X","Y","Z","T"'
!  WRITE(45,*) 'ZONE T = "BDE" , N = ',nnodes,' , E = ',ncells
!  WRITE(45,*) 'ZONETYPE = FETETRAHEDRON, DATAPACKING = BLOCK '
!  WRITE(45,*) 'VARLOCATION=(4= NODAL)'
!  WRITE(45,*) (xv(j), j = 1,nnodes)
!  WRITE(45,*) (yv(j), j = 1,nnodes)
!  WRITE(45,*) (zv(j), j = 1,nnodes) 
!  WRITE(45,*) (tnotv(j), j = 1,nnodes) 
!  DO i = 1,ncells
!     WRITE(45,*) ( lcv(i,j) , j = 1,nvcell(i) )
!  END DO
!  CLOSE( UNIT = 45 )
 
!!! 
!!!  flux_tot_adia_bot = zero
!!!    flux_tot_iso_right = zero
!!!    flux_tot_iso_left = zero
!!!  flux_tot_adia_top = zero
!!!    flux_tot_isot = zero
!!!    face_ad = 0
!!!    face_is =0
!!!    bot = 0
!!!    top =0
!!!    right = 0
!!!    left = 0
!!!   Do ifcb = 1, nbcfaces
!!!
!!!         bcon = bctype(ifcb)
!!!
!!!        Select Case(bcon)
!!!            Case(ADIA)
!!!
!!!              face_ad = face_ad + 1
!!!
!!!              currf = bf_to_f(ifcb)
!!!              ic = lfc(currf,1)    
!!!        !     ftemp = temp_bc(ifcb)
!!!              jout = zero
!!!              jin_b = zero
!!!              jin_b1 = zero
!!!              jin_d = zero
!!!        !     emit = zero
!!!
!!!
!!!
!!!      flux_net = zero
!!!         IF (threed) THEN
!!!            dist = sqrt((xc(ic)-xf(currf))**2 + (yc(ic)-yf(currf))**2 + (zc(ic)-zf(currf))**2)
!!!         ELSE
!!!            dist = sqrt((xc(ic)-xf(currf))**2 + (yc(ic)-yf(currf))**2)
!!!         ENDIF
!!!      Do iband = 1, nbands
!!!            pol = polar(iband)  
!!!            vel = gpvel(iband)
!!!            CALL relaxtime(pol, iband, tnot(ic), beta)
!!!            mfp = vel/beta  
!!!            call iwalls(pol,iband,temp_bc(ifcb),inot)
!!!            jout = jout + pi*inot
!!!        Do si = 1,ndir     
!!!            IF (threed) THEN       
!!!                sdotn = vecfx(currf) * rmu(si) + vecfy(currf) * rxi(si) + vecfz(currf)*ret(si)        
!!!                insdotn = vecfx(currf) * inrmu(si) + vecfy(currf) * inrxi(si)+ vecfz(currf)*inret(si)
!!!!                dist1 = sqrt((xc(i)-xf(currf))**2 + (yc(i)-yf(currf))**2 + (zc(i)-zf(currf))**2)
!!!             ELSE 
!!!                sdotn = vecfx(currf) * rmu(si) + vecfy(currf) * rxi(si)     
!!!                insdotn = vecfx(currf) * inrmu(si) + vecfy(currf) * inrxi(si)
!!!  !              dist1 = sqrt((xc(i)-xf(currf))**2 + (yc(i)-yf(currf))**2) 
!!!             ENDIF       
!!!             jin_b = jin_b + max(0.0,sdotn)*intensity(ic,si,iband)*(insdotn/sdotn)*exp(-dist/mfp)
!!!             jin_b1 = jin_b1 + max(0.0,sdotn)*intensity(ic,si,iband)*(insdotn/sdotn)
!!!!           emit = emit + max(0.0,-sdotn)*intensity(ic,si,iband)*(insdotn/sdotn)
!!!        ENDDO
!!!        jin_d = jin_d + phib(ifcb,iband)/two
!!!      ENDDO
!!!      flux_net = jin_b+jin_d-jout
!!!      face = bf_to_f(ifcb)
!!!  !    if (yf(face).EQ.(1.0E-006)) then !.LE.(1.00E-25)) then
!!!        if (yf(face).GT.(0.00E+00)) then !-(1.000000000000000E-006)).LE.(1.00E-25)) then
!!!        flux_tot_adia_top = flux_tot_adia_top + flux_net
!!!        top = top + 1
!!!      !  WRITE (io_flux,*) "paree"
!!!      end if
!!!      if (yf(face).EQ.(0.00E+00)) then
!!!        flux_tot_adia_bot = flux_tot_adia_bot + flux_net
!!!        bot = bot + 1
!!!      end if
!!!      
!!!       ! WRITE(io_flux,*)"flux_tot_adia "
!!!        IF (threed) THEN  
!!!        WRITE(io_flux,*) xf(currf),yf(currf),zf(currf),flux_net,t
!!!        ELSE
!!!        WRITE(io_flux,10) xf(currf),yf(currf),flux_net,jout,jin_b,jin_d  !" xf(currf),yf(currf),flux_net,jout,jin_b,jin_d, t", &
!!!                        !  xf(currf),yf(currf),flux_net,jout,jin_b,jin_d, t
!!!    !    write (io_flux,*) "top, bot", top, bot, flux_tot_adia_top , yf(face)
!!!        ENDIF
!!!  ! ENDDO
!!!
!!!
!!!
!!!!   1 close (io_output!
!!!   case (ISOT)
!!!!  WRITE(io_flux,*)" bleh", yf(face)
!!!      face_is = face_is + 1 
!!!       currf = bf_to_f(ifcb)
!!!      ic = lfc(currf,1)    
!!!!     ftemp = temp_bc(ifcb)
!!!      jout = zero
!!!      jin_b = zero
!!!      jin_b1 = zero
!!!      jin_d = zero
!!! 
!!! 
!!!      flux_net = zero
!!!         IF (threed) THEN
!!!            dist = sqrt((xc(ic)-xf(currf))**2 + (yc(ic)-yf(currf))**2 + (zc(ic)-zf(currf))**2)
!!!         ELSE
!!!            dist = sqrt((xc(ic)-xf(currf))**2 + (yc(ic)-yf(currf))**2)
!!!         ENDIF
!!!      Do iband = 1, nbands
!!!            pol = polar(iband)  
!!!            vel = gpvel(iband)
!!!            CALL relaxtime(pol, iband, tnot(ic), beta)
!!!            mfp = vel/beta  
!!!            call iwalls(pol,iband,temp_bc(ifcb),inot)
!!!            jout = jout + pi*inot
!!!        Do si = 1,ndir     
!!!            IF (threed) THEN       
!!!                sdotn = vecfx(currf) * rmu(si) + vecfy(currf) * rxi(si) + vecfz(currf)*ret(si)        
!!!                insdotn = vecfx(currf) * inrmu(si) + vecfy(currf) * inrxi(si)+ vecfz(currf)*inret(si)
!!!!                dist1 = sqrt((xc(i)-xf(currf))**2 + (yc(i)-yf(currf))**2 + (zc(i)-zf(currf))**2)
!!!             ELSE 
!!!                sdotn = vecfx(currf) * rmu(si) + vecfy(currf) * rxi(si)     
!!!                insdotn = vecfx(currf) * inrmu(si) + vecfy(currf) * inrxi(si)
!!!  !              dist1 = sqrt((xc(i)-xf(currf))**2 + (yc(i)-yf(currf))**2) 
!!!             ENDIF       
!!!             jin_b = jin_b + max(0.0,sdotn)*intensity(ic,si,iband)*(insdotn/sdotn)*exp(-dist/mfp)
!!!            ! jin_b1 = jin_b1 + max(0.0,sdotn)*intensity(ic,si,iband)*(insdotn/sdotn)
!!!!           emit = emit + max(0.0,-sdotn)*intensity(ic,si,iband)*(insdotn/sdotn)
!!!        ENDDO
!!!        jin_d = jin_d + phib(ifcb,iband)/two
!!!      ENDDO
!!!      flux_net = jout - (jin_b+jin_d)
!!!      flux_tot_isot = flux_tot_isot+ flux_net
!!!      
!!!       face = bf_to_f(ifcb)
!!!       if (abs(xf(face)-1.00E-06)<tiny) then !-(1.000000000000000E-006)).LE.(1.00E-25)) then
!!!        flux_tot_iso_right = flux_tot_iso_right + flux_net
!!!        right = right + 1
!!!      !  WRITE (io_flux,*) "paree"
!!!      end if
!!!
!!!       if (abs(yf(face)-1.00E-06)<tiny) then !-(1.000000000000000E-006)).LE.(1.00E-25)) then
!!!        flux_tot_iso_top = flux_tot_iso_top + flux_net
!!!        top = top + 1
!!!      !  WRITE (io_flux,*) "paree"
!!!      end if
!!!
!!!      if (yf(face).EQ.(0.00E+00)) then
!!!        flux_tot_iso_bot = flux_tot_iso_bot + flux_net
!!!        bot = bot + 1
!!!      end if
!!!
!!!      if (xf(face).EQ.(0.00E+00)) then
!!!        flux_tot_iso_left = flux_tot_iso_left + flux_net
!!!        left = left + 1
!!!      end if
!!!      
!!!    !  WRITE(io_flux,*)"flux_tot_isot "
!!!!      IF (threed) THEN  
!!!!      WRITE(io_flux,*) xf(currf),yf(currf),zf(currf),flux_net,t
!!!!      ELSE
!!!!      WRITE(io_flux,10) xf(currf),yf(currf),flux_net,jout,jin_b,jin_d  !" xf(currf),yf(currf),flux_net,jout,jin_b,jin_d, t", &
!!!!                        !  xf(currf),yf(currf),flux_net,jout,jin_b,jin_d, t
!!!!      ENDIF
!!!
!!! end select
!!!
!!!   ENDDO !nbcfaces
!!!!   WRITE(io_flux,*)"flux_tot_iso_right, flux_tot_iso_left, flux_tot_iso_top, flux_tot_iso_bot, jout "
!!!   WRITE(io_flux,10) flux_tot_iso_right/right, flux_tot_iso_left/left, flux_tot_iso_top/top, &
!!!                        flux_tot_iso_bot/bot !, jout !"flux_tot, t", flux_tot/nbcfaces, t
!!!!   WRITE(io_flux,*) t
!!!!   WRITE(io_flux,*)" "
!!!  
!!!!    WRITE(io_flux,*)"t,b,l,r ", top, bot, left, right
!!!!10 format(4(1x,E14.6))
!!!
!!!!   WRITE(io_flux,*)" flux_tot_adia_top, flux_tot_adia_bot "
!!!!   WRITE(io_flux,10) flux_tot_adia_top/top, flux_tot_adia_bot/bot ! flux_tot_adia/face_ad, !"flux_tot, t", flux_tot/nbcfaces, t
!!!!   WRITE(io_flux,*) t
!!!!   WRITE(io_flux,*)" "
!!!10 format(4(1x,E14.6))
!!! 

 
 
 
 
   IF(debug_level > 0) WRITE(io_debug,*) "Exiting post_process"

 
 
end subroutine post_process
  
  
  ! ****** changing file name   *************
  


    
!        if (t < 10) then
!            format_string1 = "(A7,I1,A4)"
!            format_string2 = "(A4,I1,A4)"
!        else if (t < 100) then
!            format_string1 = "(A7,I2,A4)"
!            format_string2 = "(A4,I2,A4)"
!        else if (t < 1000) then
!            format_string1 = "(A7,I3,A4)"
!            format_string2 = "(A4,I3,A4)"
!        else 
!            format_string1 = "(A7,I4,A4)"
!            format_string2 = "(A4,I4,A4)"
!        endif
!
!        write (filename1,format_string1) "contour", t, ".dat"
!        write (filename2,format_string2) "info", t, ".dat"
!        print *, trim(filename1)
!        print *, trim(filename2)
!    
!  
!    
!  OPEN( UNIT = io_output , FILE = filename1, STATUS = "REPLACE" )
!  OPEN( UNIT = io_info , FILE = filename2 , STATUS = "REPLACE" )
!  IF (threed) THEN
!        WRITE(io_output,*) 'VARIABLES = "X","Y","Z","T"'
!        WRITE(io_output,*) 'ZONE T = "BDE" , N = ',nnodes,' , E = ',ncells
!        WRITE(io_output,*) 'ZONETYPE = FEBRICK, DATAPACKING = BLOCK ' 
!        WRITE(io_output,*) 'VARLOCATION=(4=NODAL)'
!        WRITE(io_output,*) (xv(j), j = 1,nnodes)
!        WRITE(io_output,*) (yv(j), j = 1,nnodes)
!        WRITE(io_output,*) (zv(j), j = 1,nnodes)
!        WRITE(io_output,*) (tnotv(j), j = 1,nnodes)
!        
!        DO i = 1,ncells
!            WRITE(io_output,*) ( lcv(i,j) , j = 1,nvcell(i) )
!        END DO
!        
!  ELSE
!     
!        WRITE(io_output,*) 'VARIABLES = "X","Y","T"'
!        WRITE(io_output,*) 'ZONE T = "BDE" , N = ',nnodes,' , E = ',ncells
!        WRITE(io_output,*) 'ZONETYPE = FEQUADRILATERAL, DATAPACKING = BLOCK ' 
!        WRITE(io_output,*) 'VARLOCATION=(3=NODAL)'
!        WRITE(io_output,*) (xv(j), j = 1,nnodes)
!        WRITE(io_output,*) (yv(j), j = 1,nnodes)
!        WRITE(io_output,*) (tnotv(j), j = 1,nnodes)
!    
!        DO i = 1,ncells
!            WRITE(io_output,*) ( lcv(i,j) , j = 1,nvcell(i) )
!        END DO 
!    
!   ENDIF 
!  
!  ! Write(io_info,*), "Knudsen number is = " ,Kn
!   Write(io_info,*)," time is = ",t*dt
!   Write(io_info,*),"time step =",  dt 
!   Write(io_info,*),"time index is =",t
!   Write(io_info,*),"time taken per step is =", total_time/t
!   CLOSE( UNIT = io_output )
!   CLOSE( UNIT = io_info )
!
!   IF(debug_level > 0) WRITE(io_debug,*) "Exiting post_process" 
  
!end subroutine post_process
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  