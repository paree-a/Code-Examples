  Module VARIABLES
  Use PRECISIONS
  Implicit None
  SAVE

  Integer(int_p), Parameter :: nphi =100, ntheta =1 ,tmax = 100 ! this is for 2d
  Integer(int_p):: ndir,nbands
  Integer(int_p):: nla, nta,nlo,nto ! No of bands in la and ta polarizations
  Integer(int_p), Allocatable ::  IA(:), JA(:)
  Integer(int_p), Allocatable, Dimension(:) :: jlu, ju, jr
  Integer(int_p), Allocatable, Dimension(:) :: polar
  REAL(real_p), Allocatable, dimension(:):: flux_bot, flux_top
  REAL(real_p), Allocatable, dimension(:):: flux_rigt, flux_left 
  
  ! File numbers used throughout the code
  INTEGER(int_p), PARAMETER :: io_ace = 11, &   ! input file "model_setup_BTE.in"
                               io_debug = 12, & ! file for debug printout
                               io_user_input = 13, & ! file with input info for simulation
                               io_residual= 14, & ! file with  count and residual values
                               io_output = 15, &     ! file for temperature contour
                               io_flux = 16, &  ! File for heat flux output
                               io_info = 17  ! File for info on time on time steps
  ! Debug level indicates granularity of printing
  ! 0 => no printing
  ! 5 => everything printed

   
  INTEGER(int_p), PARAMETER :: debug_level = 0, &  ! Debug Level
                               max_iter = 200, &   ! Maximum number of outer iterations
                               iter_gmres = 100, &    ! Maximum number of GMRES iterations
                               freq_time_output = 1 ! Output frequency
  
  REAL(real_p), PARAMETER :: relax_temp = 0.45d0, & ! Relaxation factor for Temperature
                             tol_inner = 1.0D-3, & ! GMRES tolerance
                             tol_outer = 1.0D-4  ! Outer loop tolerance
  
  Real(4) :: time_start,time_end,total_time,ta(2) 
  Real(real_p) :: residual,dt,idt, T_ref, tau_l
  Real(real_p):: dtheta,dphi,hotwall,coldwall 
  Real(real_p), Parameter :: initial_temp = 245.0d0
    Real(real_p):: dos, cutoff
  Real(real_p):: delta_la,delta_ta
  REAL(real_p), DIMENSION(:,:), POINTER :: gw,gwone
  REAL(real_p), DIMENSION(:,:), POINTER :: phib,phi,phiv,jt,tb
   REAL(real_p),DIMENSION(:,:),POINTER:: phione,phitwo,phibone
  Real(real_p), Allocatable, Dimension(:) :: alu
  Real(real_p), Allocatable, Dimension(:) :: rmu,rxi,ret,inrmu,inrxi,inret
  Real(real_p), Allocatable, Dimension(:) :: Omega
  Real(real_p), Allocatable, Dimension(:) :: gpvel
  Real(real_p), Allocatable, Dimension(:) :: centfreq
  Real(real_p), Allocatable  ::  Am(:), R(:), soln(:)
  Real(real_p), Dimension(:), Pointer :: tnot,tprime,tnotv,tbprime,tnotone ! flux
  Real(real_p), Dimension(:), Pointer :: resid,s
  Real(real_p), Dimension(:), Pointer :: scs,aps 
  Real(real_p), Dimension(20):: xi,wi   
  Real(real_p), Dimension(:,:) , Pointer :: anbs , flux
  Real(real_p), Allocatable, Dimension(:,:) :: vv
  Real(real_p), Dimension(:,:), Pointer :: gc,gcone  
  Real(real_p), Dimension(:,:,:), Pointer :: intensity,intensityone 
  LOGICAL:: jhanda
  End Module variables
