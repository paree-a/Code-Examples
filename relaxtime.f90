!!$ Ad astra per aspera
!!$ this subroutine calculates the band avearges scattering time scales for phonons 
!!$ as per the expressions given by holland
!!$ refer Mazumder and Majumdar JHT 2001, for details

  Subroutine relaxtime(pol, band, temp,beta)
  
  Use precisions !,Only:int_p,real_p 
  Use constants !,only:wminla,one,bl,wminta,two,btn,zero,wmax_half,hobol,half,btu,third
  Use variables !,only:delta_la,delta_ta,nla,xi,wi,centfreq, io_debug, debug_level
  Implicit None 

  Integer(INT_P),Intent(IN):: pol, band
  
  Real(REAL_P):: nor,umk,freval
  Real(REAL_P),Intent(IN):: temp
  Real(REAL_P),Intent(OUT):: beta 


  
    Select Case(pol) 

    Case(0) !!$ LA phonon
        freval = centfreq(band)      
        beta = bl*(T_ref**3)*(freval**2)  

    Case(1)  !!$ TA phonon     
        freval = centfreq(band) 
        nor   = btn*(T_ref**4)*freval
        umk = zero
        
        If (freval .ge. wmax_half) Then                  
            umk = btu*(freval**2)/(Sinh(hobol * freval/T_ref))             
        End If
        
        beta = nor + umk 

  End Select
  

  End Subroutine relaxtime


