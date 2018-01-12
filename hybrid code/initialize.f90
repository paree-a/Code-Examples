Subroutine Initialize



  USE precisions 
  USE variables
  USE grid !,only:ncells,nnodes,nbcfaces,nfaces,temp_bc
  USE constants !,only:zero,sigma,pi,four,one,two
  Implicit none  
  integer(int_p) :: dim
  real(real_p)::phiguess,gval
!  nla = 1
!  nta = 0
  nla = 20
  nta = 20
!  mfp = Kn*1.0d0
!  tau = mfp/gpvel ! the relaxation time       
  dim = ncells + nbcfaces
  nbands = nla + nta
  ALLOCATE(tnot(ncells),tnotone(ncells),tprime(ncells),tnotv(nnodes))! dont change this
  ALLOCATE(gw(ncells,nbands),gwone(ncells,nbands))
  ALLOCATE(resid(dim),s(ncells) )
  AlLOCATE(r(dim),soln(dim) )
!  ALLOCATE(aps_global(dim),scs_global(dim),anbs_global(dim,6))
  ALLOCATE(am(dim*7),ja(dim*7),ia(dim+1))
!  ALLOCATE(ju_global(dim*7),vv_global(dim,21),jlu_global(dim*7),jr_global(2*dim),alu_global(dim*7))  
  ALLOCATE(intensity(ncells,ndir,nbands),intensityone(ncells,ndir,nbands))
!  ALLOCATE(jw(nbcfaces))
  Allocate(tbprime(nbcfaces),phibone(nbcfaces,nbands),phione(ncells,nbands),phitwo(ncells,nbands))
  ALLOCATE( phiv(nnodes,nbands),phib(nbcfaces,nbands),phi(ncells,nbands),jt(nfaces,nbands),tb(nbcfaces,nbands))!,jw(nbcfaces)) 
!  ALLOCATE(resid(dim))
!  AlLOCATE(r(dim),soln(dim) )
  ALLOCATE( aps(dim), scs(dim), anbs(dim,6))
!  ALLOCATE( am(dim*7),ja(dim*7),ia(dim+1))
  ALLOCATE( ju(dim*7) , vv(dim,21) , jlu(dim*7) , jr(2*dim), alu(dim*7)) 
  
  ! **********added new
  Allocate(polar(nbands))
  Allocate(gpvel(nbands))  
  Allocate(centfreq(nbands)) 
  
         Allocate(flux_bot(150))
       Allocate(flux_top(150))
      Allocate(flux_rigt(150))
      Allocate(flux_left(150))     
  
  IF(debug_level > 0) WRITE(io_debug,*) "Starting initialize" 
  



 ! phiguess = sigma*(initial_temp**4)/pi
! write(*,*), phiguess
  !gval = four*pi*phiguess
  s(:) = zero   
  phi(:,:) = zero
  phione(:,:) = zero
  phitwo(:,:) = zero
  phib(:,:) = zero
  phibone(:,:) = zero

  gw(:,:) = zero
  gwone(:,:) = zero
  intensity(:,:,:) = zero
  intensityone(:,:,:) =zero
  
  tnot(:) = initial_temp
  tnotone(:) = initial_temp   
  tprime(:) = zero
  tbprime(:) = zero
!  dt = 5.0d-12
 dt = 8.5723d-11
  idt = one/dt 
  T_ref = 250.0d0
  tau_l = 1.0d-06
   
   
   
     !****changed

  
   
   
  IF(debug_level > 0) WRITE(io_debug,*) "Exiting initialize"  
    
END subroutine Initialize
