!!$ Ad astra per aspera
!!$ This subroutine calculates the phonon intensity 
!!$ The inputs to this subroutine are polarization, band index and temperature

  Subroutine iwalls(pol, band, temp,actint)
  Use precisions
  Use constants !only:zero,two,one,wminla,vsla,four,cla,dirac,half,eight,pi,wminta,vsta,cta,hobol
  Use variables !,only: delta_la,delta_ta,wi,xi,nla, io_debug, debug_level
  Implicit none

  Integer(int_p):: i
  Integer(int_p), Intent(In):: pol,band
 
  Real(real_p):: fmin,fmax, a, b, c1, c2, yi, wavenum, disfun, func
  Real(real_p), Intent(In)::temp  
  Real(real_p), INTENT(out):: actint
 
  !IF(debug_level > 0) WRITE(io_debug,*) "Starting iwallcal_iwall"
   
  actint = zero
  
  Select case (pol)
    case(0) ! la phonons
        fmin = wminla + (band-1)*delta_la
        fmax = fmin + delta_la
        a = fmin
        b = fmax
        c1 = two/(b-a)
        c2 = one-c1*b
        DO i = 1,20  ! gauss-legendre quadriture
            yi = (xi(i) - c2)/c1 
            wavenum = (( -vsla + sqrt((vsla**2) + four*(yi)*cla))/(two*cla))
!            wavenum = yi/vsla
            disfun = one/(exp(hobol*yi/(temp)) - one)
            func = yi*disfun*(wavenum**2)   
            actint = actint + wi(i)*func
        ENDDO
        actint = half*actint*(b-a)
        actint = actint*dirac/(eight*(pi**3))    


!actint = 0.2826d0*(temp**4)/pi
!
        !WRITE(28,*) "c1,c2,actint la phonons", c1,c2, actint

    case(1) ! ta phonons
        fmin = wminta + (band-1-nla)*delta_ta
        fmax = fmin + delta_ta
        a = fmin
        b = fmax
        c1 = two/(b-a)
        c2 = one-c1*b     
        DO i = 1,20
            yi = (xi(i) - c2)/c1 
            wavenum = (( -vsta + sqrt((vsta**2) + four*(yi)*cta))/(two*cta))
!             wavenum = yi/vsta
            disfun = one/(exp(hobol*yi/(temp)) - one)
            func = yi*disfun*(wavenum**2)     
            actint = actint + wi(i)*two*func
        ENDDO
        actint = half*actint*(b-a)
        actint = actint*dirac/(eight*(pi**3))

        !WRITE(28,*) "c1,c2,actint ta silicon phonons", c1,c2, actint
!  
  END select
  
 ! ADD other band if neccessary--assy

 !IF(debug_level > 0) WRITE(io_debug,*) "Finishing iwallcal_iwall" 

  End subroutine iwalls










