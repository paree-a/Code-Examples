!**********************************************************************
subroutine dircalc



  !**********************************************************************
  ! purpose :
  !         Main Program to initiate choice of discrete ordinate 
  !         quadrature scheme
  !**********************************************************************
  use precisions
  use variables      
  use constants !,only: half,two,pi,zero
  implicit none

  real(real_p) :: theta,anphi,sum
  integer(int_p) :: i, count,j

  Allocate(Omega(Ntheta*Nphi))
  Allocate(rmu(Ntheta*Nphi),rxi(Ntheta*Nphi),ret(Ntheta*Nphi))
  Allocate(inrmu(Ntheta*Nphi),inrxi(Ntheta*Nphi),inret(Ntheta*Nphi))    
  dphi = two*pi/Nphi
  dtheta = pi/Ntheta
  ndir = int(Ntheta*Nphi)
  IF(debug_level > 0) WRITE(io_debug,*) "Starting snweights file and dircalc" 
 
  count = 0
  sum =  zero
  

  
  do i = 1, Ntheta
     theta = (i-half)*dtheta     
     do j = 1,Nphi
        anphi = (j-half)*dphi
        count = count + 1
        Omega(count) = two*sin(theta)*sin(half*dtheta)*dphi
        sum = sum + Omega(count)
        rmu(count) = sin(theta)*sin(anphi)
        rxi(count) = sin(theta)*cos(anphi)
        ret(count) = cos(theta)

        inrmu(count) = sin(anphi) * sin(dphi*half) * (dtheta - cos(two*theta)*sin(dtheta))
        inrxi(count) = cos(anphi) * sin(dphi*half) * (dtheta - cos(two*theta)*sin(dtheta))
        inret(count) = half * dphi * sin(two*theta) * sin(dtheta)

     enddo
  enddo

 IF(debug_level > 0) WRITE(io_debug,*) "Exiting snweights and dircalc" 

end subroutine dircalc
