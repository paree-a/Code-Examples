!!$ This subroutine calculates the band averaged velocity of each frequency band and allocates integer values to each of the polarizations
  Subroutine avgvel
  Use precisions,Only:real_p,int_p
  Use constants,only:wmaxla,wminta,wmaxta,wminla,vsla,vsta,cla,cta,six,four,half
  Use variables,Only:delta_la,delta_ta,nla,nta,gpvel,polar,centfreq, io_debug, debug_level
  Implicit None

  Integer(int_p)::i
  Real(real_p) :: flow,fhigh,freval
  
  IF(debug_level > 0) WRITE(io_debug,*) "Starting avgvelocity" 

!!$ the input to this subroutine is the number of LA bands NLA and number of TA bands NTA
!!$ assuming that the bins are equally spaced
  delta_la = (wmaxla - wminla)/nla
  delta_ta = (wmaxta - wminta)/nta 

!!$ Assuming quadratic dispersion relation:
!!$ LA Branch
  Do i = 1, nla 
     flow = wminla + (i-1)*delta_la    
     fhigh = flow + delta_la 
     freval = (flow + fhigh)*half
     centfreq(i) = freval
     gpvel(i) = sqrt((vsla**2.0) + four * cla * freval)    
     polar(i) = 0 !!$ Assigning integer values to various polarizations
  End Do

!!$ TA Branch
  Do i = 1, nta
     flow = wminta + (i-1)*delta_ta
     fhigh = flow + delta_ta
     freval = (flow + fhigh)*half
     centfreq(i + nla) = freval
     gpvel(nla + i) = sqrt((vsta**2.0) + four * cta * freval) 
     polar(i + nla) = 1
  End Do
!!$ Add other polarization branches if required
  
  IF(debug_level > 0) WRITE(io_debug,*) "Finished avgvelocity" 
  
  End Subroutine avgvel
