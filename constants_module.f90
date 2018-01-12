!***********************************************************************
! Definitions of basic constants

  Module constants
  Use precisions,Only:real_p
  Implicit None
  Save

  Real(real_p), Parameter :: zero=0.0d0,one=1.0d0,two=2.0d0,three=3.0d0,four=4.0D0, & 
                             five=5.0d0, Eight = 8.0d0, nine = 9.0d0, &
                             half=0.5d0,third=0.33333333d0, fourth = 0.25D0, &
                             pi = 3.1415926535898d0, tiny =1.0D-20,tol =1.0d-8, &
                             tzero = 1.0d-15,six = 6.0d0

  Real(real_p):: dirac = 1.054571628d-34, boltz = 1.3806503d-23, &
                 hobol = 7.63822401661014d-12
                 
  ! parameters related to Silicon dispersion from page 39
!  Real(real_p),Parameter:: wminta = 0.0d0, wminla = 0.0d0, wmaxla=7.728337675901222d13,& 
!                           wmaxta = 7.728337675901222d13, cla=0.0d0,vsla=6.67d03, &
!                           cta = 0.0d0,vsta = 6.67d03,wmaxlo=9.88d13, &
!                           wmaxtos= 10.20d13,wmintos = 8.7124d13, vslos = 0.0d0, &
!                           clos = -1.6d-7, vstos = -2.57d3, ctos= 1.11d-7


  Real(real_p),Parameter:: wminta = 0.0d0, wminla = 0.0d0, wmaxla=7.728337675901222d13,& 
                           wmaxta = 2.97927417405992d13, cla=-2.0d-7,vsla=9.01d3, &
                           cta = -2.26d-7,vsta= 5.23d3,wmaxlo=9.88d13, &
                           wmaxtos= 10.20d13,wmintos = 8.7124d13, vslos = 0.0d0, &
                           clos = -1.6d-7, vstos = -2.57d3, ctos= 1.11d-7
                           
  ! relaxation time scale related parameters eqn 2.32, 2.33
  Real(real_p), PARAMETER::wmax_half = 2.417d13, bl = 2.0d-24, btu = 5.5d-18, &
                           btn = 9.3d-13, grun = 0.59d0, lattice_const = 5.43d-10, &
                           rho = 2.33d3, toler = 1.0d-04

  End Module constants
!***********************************************************************
