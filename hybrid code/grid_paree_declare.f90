!***********************************************************************
! Parameters and data used ! 3d grid

  Module GRID

  Use PRECISIONS
  Implicit None
  SAVE
  ! Declare global variables
  ! USER GLOBAL VARIABLE DECLARATION BEGIN
  ! number of cells, faces etc.
  Integer(int_p) :: ncells,nfaces,nnodes,nbcfaces, ndomains, ndomains_rad
  INTEGER :: nf_max, nbc_patches       
  INTEGER(int_p), PARAMETER :: ADIA = 4, ISOT = 3

  ! arrays relating to boundary faces.
  ! bface stores bface index
  ! f_to_bf is mapping from global face index to bface
  ! bf_to_f is mapping from bface to global face index
  Integer(int_p),Dimension(:), Allocatable :: f_to_bf, bf_to_f, bf_to_c, &
                                              nfcell, bctype, bcname, nvcell
  INTEGER(int_p), DIMENSION(:), ALLOCATABLE :: bnode,bface
  INTEGER (int_p), DIMENSION(:), ALLOCATABLE :: nvface, no_f_patch 
  INTEGER (int_p), DIMENSION(:,:), ALLOCATABLE :: normdir 
  Integer(real_p), Dimension(:,:), Allocatable :: lcf, lcv, lfv, lfc, lcc 

  ! cell center coordinates, face center coordinates, vertex coordinates, face area
  REAL(real_p), Dimension(:), Allocatable :: xc, yc, zc, xf, yf, zf, xv, yv, zv, areaf, &
                                             vecfx, vecfy, vecfz, volcell,  tgtx,tgty, tgtz, nodeweight, tdotlf
  REAL(real_p), DIMENSION(:), ALLOCATABLE :: disn,phibc, wcf
  REAL(real_p), DIMENSION(:,:), ALLOCATABLE :: wcv,wbfv
  Real(real_p), Dimension(:), ALLOCATABLE :: temp_bc
  
  LOGICAL :: threed = .true.
  
  End Module GRID
    
