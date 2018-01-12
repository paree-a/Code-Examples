
!!$   hybrid Solver for solving the non gray boltzmann transport equation
!!$   Fluids and Thermal Analysis Lab, The Ohio State University
PROGRAM Main
  USE precisions
  USE constants !,ONLY:one,toler,zero,two,sigma
  USE grid
  USE variables
  USE ifport
  IMPLICIT NONE
  INTEGER(int_p)::i,tindex,ifcb,bcon,currf,ic,si,count,iband, pol, iter
  real(real_p)::dist,jout,sdotn,insdotn, phiguess, gval, phi_temp,mfp,vel,beta, t_norm, sum, tiny2, Kn
 Real(4) :: tstep_1,tstep_2

  OPEN(UNIT = 28 , FILE = "calculations.txt" ,STATUS = "REPLACE" )
  OPEN( UNIT = 29 , FILE = "calculations2.dat" ,STATUS = "REPLACE" )
  OPEN( UNIT = 12 , FILE = "debug_out.txt",STATUS = "REPLACE" )
  OPEN( UNIT = io_flux , FILE = "flux.dat" ,STATUS = "REPLACE" )
  OPEN( UNIT = 50 , FILE = "residual.txt" , STATUS = "REPLACE"  )
   
  CALL GRID_paree
  CALL dircalc
  Call gauss_points

  Call Initialize  
  Call avgvel 
 

  IF(debug_level > 0) WRITE(io_debug,*) "Starting initial_1" 
   
  cutoff = 0.1d0
  call knudsen_plot
  
  
          do iband = 1,nbands

            pol = polar(iband)
            vel = gpvel(iband)
            Call relaxtime(pol, iband, T_ref, beta)
            call iwalls(pol, iband, initial_temp, phiguess)  
            mfp = vel/beta   
            Kn = mfp/tau_l
            gval = four*pi*phiguess
            if (Kn .ge. cutoff) then
                gw(:,iband) = gval
                gwone(:,iband) = gval
                intensity(:,:,iband) = phiguess
                intensityone(:,:,iband) = phiguess            
            else
              phi(:,iband) = zero
              phione(:,iband) = gval
              phitwo(:,iband) = zero                 
              phibone(:,iband) = gval  
            endif
        
        enddo
  

  tiny2 = 1.0d-09
  total_time = 0.0d0
  jhanda = .true.      
     
  IF(debug_level > 0) WRITE(io_debug,*) "Starting initial_3" 
    
    time_start   = etime(ta)
  
  DO tindex = 1,tmax 

     tstep_1 = etime(ta)
     residual = one

     iter = 0
     DO WHILE (residual > toler) 
       residual = zero        
        iter = iter + 1
        do iband = 1,nbands

            pol = polar(iband)
            vel = gpvel(iband)
            Call relaxtime(pol, iband, T_ref, beta)
            mfp = vel/beta   
            Kn = mfp/tau_l

            
            if (Kn .ge. cutoff) then
               call disord(pol,iband,mfp)
            else
                call mediasolve2(pol,iband,vel)
            endif
        
        enddo
        CALL tempcalc 
        
        OPEN( UNIT = 50 , FILE = "residual.txt" , ACCESS = 'APPEND',STATUS = 'Unknown' )
        write(50,*), tindex,residual
        write(*,*), tindex,residual
        close(50)
        DO i = 1,ncells                  
           tnot(i) = tnot(i)+ tprime(i)

        ENDDO

        DO i = 1,nbcfaces                 
           temp_bc(i) = temp_bc(i)+ tbprime(i)

        ENDDO

     END DO  ! end residual check

      sum = zero

      CALL update
     
    do i= 1,ncells

      if (abs(xc(i)- 0.49500E-06)<tiny2) then
     
          WRITE(29,*)tnot(i),yc(i),tindex
          endif
          write(28,*) "tnot", tnot(i)
     enddo
     
     write(28,*)   " tindex", tindex
     write(28,*) ""
     write(28,*) ""
     
     write(29,*) " tindex", tindex
     write(29,*) ""
     write(29,*) ""
     
     tstep_2 = etime(ta)
     WRITE(*,*) "Time taken for current time step =", tstep_2-tstep_1, tindex
     WRITE(io_debug,*) "Time taken for current time step =", tstep_2-tstep_1, tindex

     CALL post_process(tindex) 
     
     call flux_calc
         
     jhanda = .false.
  ENDDO
  
  time_end = etime(ta)
  total_time =  (time_end - time_start)   
  WRITE(*,*) "Total time taken =", total_time
  WRITE(io_debug,*) "Total time taken =", total_time
  
CLOSE (io_debug)
close (28)
close (29)
close (io_flux)


END PROGRAM main
