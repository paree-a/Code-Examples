  Subroutine reflec(sin,ibface,sdotn,sout)

  Use precisions
  Use constants !,only:two
  Use grid !,only:vecfx,vecfy,vecfz
  Use variables !,only:rmu,rxi,ndir,rxi,ret, io_debug, debug_level 

  Implicit None 

  INTEGER(int_p)::i,sout
  INTEGER(int_p), INTENT(in) :: sin,ibface
  
  REAL(real_p) :: srx,sry,srz
  REAL(real_p), INTENT(in):: sdotn
  
  
 ! IF(debug_level > 0) WRITE(io_debug,*) "Starting reflection_reflec" 

  ! now the components of intensity are known the value of G( in this case theta_b) at cell centers needs to be calculated   ! calculation of G
  srx = rmu(sin) - two*sdotn*vecfx(ibface)
  sry = rxi(sin) - two*sdotn*vecfy(ibface)
  If (threed) then 
  srz = ret(sin) - two*sdotn*vecfz(ibface)
  endif
    sout = 0
  do i = 1, ndir
     if (threed) then
         if ((abs(srx - rmu(i)) .lt. 1.0d-6) .and. (abs(sry - rxi(i)).lt. 1.0d-6) .and. (abs(srz - ret(i)).lt. 1.0d-6) ) then
                sout = i
                exit       
         endif
     else
         if ((abs(srx - rmu(i)) .lt. 1.0d-6) .and. (abs(sry - rxi(i)).lt. 1.0d-6) ) then
                sout = i
                exit       
         endif
   
     endif
  end do
  if (sout==0) then
        stop
  endif
  
 ! IF(debug_level > 0) WRITE(io_debug,*) "Finishing reflection_reflec" 
 
  END subroutine reflec