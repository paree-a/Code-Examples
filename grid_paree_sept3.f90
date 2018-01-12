SUBROUTINE GRID_paree
    
  USE precisions
  USE grid
  USE constants !, ONLY: zero,one
  USE variables !, ONLY: io_ace,io_debug,debug_level
    
  ! Declare global variables
  ! USER GLOBAL VARIABLE DECLARATION BEGIN

  INTEGER(int_p) :: icell, iface, inode, ibf, i, j, k, cell_1, cell_2, node_1, &
                    node_2, ifc, ic, ibc, throw_away1, throw_away2
  INTEGER(int_p) :: nv_max1,nv_max2, ierr,flag2
    
  REAL(real_p) :: xtemp1, xtemp2, disntemp
  REAL(real_p), DIMENSION(:), ALLOCATABLE :: sumweight, sumwtf
    
  LOGICAL :: exists = .true.
  LOGICAL :: flag1 = .false.
    
  IF(debug_level > 0) WRITE(io_debug,*) "Starting grid_and_boundary" 

  ncells = 0
  nfaces = 0
  nnodes = 0
  nbcfaces = 0    
   
            
!  INQUIRE(file="model_setup_BTE_100cells_120_80.in", EXIST = exists)
    INQUIRE(file="model_setup_BTE_trial_all100.in", EXIST = exists)
  IF (.not. exists) THEN
    print *, "File model_setup.in could not be found"
    STOP
  ENDIF
  
!  OPEN (UNIT=io_ace,FILE='model_setup_BTE_100x100_295k.in', status = "old")
 OPEN (UNIT=io_ace,FILE='model_setup_BTE_22k_cells_255k.in', status = "old")
!

!  INQUIRE(file="model_setup_BTE.in", EXIST = exists)
!  IF (.not. exists) THEN
!    print *, "File model_setup.in could not be found"
!    STOP
!  ENDIF
!  
!  OPEN (UNIT=io_ace,FILE='model_setup_BTE.in', status = "old")

            
  READ(io_ace,*) ncells,nfaces,nnodes,nbcfaces
             
  IF(debug_level > 0) WRITE(io_debug,*) "Finished reading global grid size data"
                             
  READ(io_ace,*) throw_away1, throw_away2
  IF(throw_away1 == 0) THEN
    threed = .false.
  ELSE
    threed = .true.
  ENDIF

  IF(throw_away2 == 0) THEN
    flag1 = .false.
  ELSE
    flag1 = .true.
  ENDIF

  IF(debug_level > 0) WRITE(io_debug,*) &
        "Finished reading other options", threed, flag1       
                                
 
  ALLOCATE(nfcell(ncells), stat=ierr)
  ALLOCATE(nvcell(ncells), stat=ierr)
  DO i=1,ncells
    READ(io_ace,*) nfcell(i),nvcell(i)
  END DO
            
  IF(debug_level > 0) WRITE(io_debug,*) &
      "Finished reading number of faces and vertices for each cell"

               
! Reading in geometric data
  IF(threed)THEN  !3D
  
    ALLOCATE (xc(ncells), stat=ierr)
    ALLOCATE (yc(ncells), stat=ierr)
    ALLOCATE (zc(ncells), stat=ierr)
    ALLOCATE (volcell(ncells), stat=ierr)
          
    DO i=1,ncells
      READ(io_ace,*) xc(i), yc(i), zc(i), volcell(i)
    END DO                 

    ALLOCATE (xf(nfaces), stat=ierr)
    ALLOCATE (yf(nfaces), stat=ierr)
    ALLOCATE (zf(nfaces), stat=ierr)
      
    ALLOCATE (vecfx(nfaces), stat=ierr)
    ALLOCATE (vecfy(nfaces), stat=ierr)
    ALLOCATE (vecfz(nfaces), stat=ierr)
    ALLOCATE (areaf(nfaces), stat=ierr)

! Allocate tangents
    ALLOCATE (tgtx(nfaces), stat=ierr)
    ALLOCATE (tgty(nfaces), stat=ierr)
    ALLOCATE (tgtz(nfaces), stat=ierr)
          
    DO i=1,nfaces
      READ(io_ace,*) xf(i),yf(i),zf(i),areaf(i),vecfx(i),vecfy(i),vecfz(i)
    END DO       

    ALLOCATE (xv(nnodes), stat=ierr)
    ALLOCATE (yv(nnodes), stat=ierr)
    ALLOCATE (zv(nnodes), stat=ierr)
    DO i=1,nnodes
      READ(io_ace,*) xv(i), yv(i), zv(i)
    END DO
          
  ELSE  ! 2D
      
  ! Insert 2D
    ALLOCATE (xc(ncells), stat=ierr)
    ALLOCATE (yc(ncells), stat=ierr)
    ALLOCATE (volcell(ncells), stat=ierr)
          
    DO i=1,ncells
      READ(io_ace,*) xc(i), yc(i), volcell(i)
    END DO                 

    ALLOCATE (xf(nfaces), stat=ierr)
    ALLOCATE (yf(nfaces), stat=ierr)
      
    ALLOCATE (vecfx(nfaces), stat=ierr)
    ALLOCATE (vecfy(nfaces), stat=ierr)
    ALLOCATE (areaf(nfaces), stat=ierr)

! Allocate tangents
    ALLOCATE (tgtx(nfaces), stat=ierr)
    ALLOCATE (tgty(nfaces), stat=ierr)
    
    ! Only for 2D case
    ALLOCATE (tdotlf(nfaces), stat=ierr)
      
    DO i=1,nfaces
      READ(io_ace,*) xf(i),yf(i),areaf(i),vecfx(i),vecfy(i)
    END DO       

    ALLOCATE (xv(nnodes), stat=ierr)
    ALLOCATE (yv(nnodes), stat=ierr)
    DO i=1,nnodes
      READ(io_ace,*) xv(i), yv(i)
    END DO
    
  ENDIF ! 2D/3D

  nf_max = -1
  DO i=1,ncells
    nf_max = MAX(nf_max,nfcell(i))
  END DO
  ALLOCATE(lcf(ncells,nf_max), stat=ierr)
  DO i=1,ncells
    READ(io_ace,*) (lcf(i,j), j=1,nfcell(i))
  END DO
  IF(debug_level > 0) WRITE(io_debug,*) &
                       "Finished reading cell to face connectivity"
  
  nv_max1 = -1
  DO i=1,ncells
    nv_max1 = MAX(nv_max1,nvcell(i))
  END DO
  ALLOCATE(lcv(ncells,nv_max1), stat=ierr)
  DO i=1,ncells
    READ(io_ace,*) (lcv(i,j), j=1,nvcell(i))
  END DO
  IF(debug_level > 0) WRITE(io_debug,*) &
                       "Finished reading cell to vertex connectivity" 
  
  ALLOCATE(lfc(nfaces,2), stat=ierr)
  DO i=1,nfaces
    READ(io_ace,*) lfc(i,1), lfc(i,2)
  END DO
  IF(debug_level > 0) WRITE(io_debug,*) &
                       "Finished reading face to cell connectivity"  
  
  ALLOCATE(nvface(nfaces), stat=ierr)
  nv_max2 = -1
  DO i=1,nfaces
    READ(io_ace,*) nvface(i)
    nv_max2 = MAX(nv_max2,nvface(i))
  END DO
  
  ALLOCATE(lfv(nfaces,nv_max2), stat=ierr)
  DO i=1,nfaces
    READ(io_ace,*) (lfv(i,j), j=1,nvface(i))
  END DO
  IF(debug_level > 0) WRITE(io_debug,*) &
                       "Finished reading face to vertex connectivity"    
 
  ALLOCATE(lcc(ncells,nf_max), stat=ierr)
  DO i=1,ncells
    READ(io_ace,*) (lcc(i,j), j=1,nfcell(i))
  END DO
  IF(debug_level > 0) WRITE(io_debug,*) &
                       "Finished reading cell type and cell to cell connectivity" 
  
! Reading in boundary data
  READ(io_ace,*) nbc_patches
  ALLOCATE(no_f_patch(nbc_patches), stat=ierr)  ! No of faces per patch
  
  ALLOCATE (phibc(nbcfaces), stat=ierr)
  ibf = 0 
  DO j = 1, nbc_patches
    READ(io_ace,*) no_f_patch(j)
    DO i = 1, no_f_patch(j)
    ! Reading values of all variables at boundaries
    ibf = ibf + 1
    READ(io_ace,*) phibc(ibf)
    ENDDO
  ENDDO
  IF(debug_level > 0) WRITE(io_debug,*) &
                        "Finished reading boundary data"

! Reading Boundary face flag and global face to boundary face conn
  ALLOCATE(bface(nfaces), stat=ierr)
  ALLOCATE(f_to_bf(nfaces), stat=ierr)
  DO i=1,nfaces
    READ(io_ace,*) bface(i),f_to_bf(i)
  END DO
  
 ! Reading Boundary type and boundary face to global face and cell 
  ALLOCATE(bf_to_f(nbcfaces), stat=ierr)
  ALLOCATE(bf_to_c(nbcfaces), stat=ierr)
  ALLOCATE(bctype(nbcfaces), stat=ierr)
  ALLOCATE(bcname(nbcfaces), stat=ierr)
  DO i=1,nbcfaces
    READ(io_ace,*) bctype(i), &
                      bf_to_f(i),bf_to_c(i)
!                      write(28,*)"bctype",bctype(i)
  END DO
  IF(debug_level > 0) WRITE(io_debug,*) &
                       "Finished reading boundary connectivity"
  
  
! Finished all Reading. Now time to process data
! Process Distance
  ALLOCATE (disn(nfaces), stat=ierr)
  IF(threed)THEN
     DO iface=1,nfaces
        cell_1 =  lfc(iface,1)
        cell_2 =  lfc(iface,2)  
        IF(cell_1 == cell_2) THEN  ! External faces
           disn(iface) = ABS((xf(iface)-xc(cell_1))*vecfx(iface) + &
                             (yf(iface)-yc(cell_1))*vecfy(iface) + &
                             (zf(iface)-zc(cell_1))*vecfz(iface))
        ELSE  ! Interior face
           disn(iface) = ABS((xc(cell_2)-xc(cell_1))*vecfx(iface) + &
                             (yc(cell_2)-yc(cell_1))*vecfy(iface) + &
                             (zc(cell_2)-zc(cell_1))*vecfz(iface))
        END IF
     END DO
  ELSE
     DO iface=1,nfaces
        cell_1 =  lfc(iface,1)
        cell_2 =  lfc(iface,2)  
        IF(cell_1 == cell_2) THEN
           disn(iface) = ABS((xf(iface)-xc(cell_1))*vecfx(iface) + &
                             (yf(iface)-yc(cell_1))*vecfy(iface))
        ELSE
           disn(iface) = ABS((xc(cell_2)-xc(cell_1))*vecfx(iface) + &
                             (yc(cell_2)-yc(cell_1))*vecfy(iface))
        END IF
     END DO
  ENDIF
     IF(debug_level > 0) WRITE(io_debug,*) &
                         "Finished calculating disn"

! ********* ASK pROF.*************
! MODIFIED OCT 3
!   Calculating the tangents used in 2D code
    IF (.not.threed) THEN
        Do i = 1,nfaces
           node_1 = lfv(i,1)
           node_2 = lfv(i,2)
           tgtx(i) = ( xv(node_2) - xv(node_1) ) / ( Sqrt( (xv(node_2) - xv(node_1))**2 + (yv(node_2) - yv(node_1))**2 ) )
           tgty(i) = ( yv(node_2) - yv(node_1) ) / ( Sqrt( (xv(node_2) - xv(node_1))**2 + (yv(node_2) - yv(node_1))**2 ) )
        End Do
        
        Do iface=1,nfaces

           cell_1 =  lfc(iface,1)
           cell_2 =  lfc(iface,2)  
           If(cell_1 .Eq. cell_2) Then
    !          disn(iface) = Abs((xf(iface)-xc(cell_1))*vecfx(iface) + (yf(iface)-yc(cell_1))*vecfy(iface))
              tdotlf(iface) = ( tgtx(iface) * ( xf(iface) - xc(cell_1) ) + tgty(iface) * ( yf(iface) - yc(cell_1) ) )          
              !disn(iface) = ABS(SQRT((xf(iface) - xc(cell_1))**2.0 + (yf(iface) - yc(cell_1))**2.0))
           Else
    !          disn(iface) = Abs((xc(cell_2)-xc(cell_1))*vecfx(iface) + (yc(cell_2)-yc(cell_1))*vecfy(iface))
              tdotlf(iface) = ( tgtx(iface) * ( xc(cell_2) - xc(cell_1) ) + tgty(iface) * ( yc(cell_2) - yc(cell_1) ) )
              !disn(iface) = ABS(SQRT((xc(cell_1) - xc(cell_2))**2.0 + (yc(cell_1) - yc(cell_2))**2.0))
           End If
        End Do
    
    ENDIF
    
! ***********  ENDS HERE ***************                   


! Finding orientation of normal at face relative to CV
  ALLOCATE (normdir(ncells,nf_max), stat=ierr)
  normdir(:,:) = 1
  DO i=1,ncells
     DO j =1,nfcell(i)
        iface = lcf(i,j)
        IF(lfc(iface,1) == i) THEN
           normdir(i,j) = 1   ! Outward already
        ELSE
           normdir(i,j) = -1  ! Flip to make outward
        END IF
     END DO
  END DO

! Compute Interpolant from cell to face
  ALLOCATE(wcf(nfaces), stat=ierr)   
  DO iface=1,nfaces
     cell_1 = lfc(iface,1)
     cell_2 = lfc(iface,2)
     xtemp1 = SQRT((xf(iface) - xc(cell_1))**2.0 + (yf(iface) - yc(cell_1))**2.0)
     IF(threed) xtemp1 = SQRT((xf(iface) - xc(cell_1))**2.0 + &
                              (yf(iface) - yc(cell_1))**2.0 + &
                              (zf(iface) - zc(cell_1))**2.0)
     xtemp2 = SQRT((xf(iface) - xc(cell_2))**2.0 + (yf(iface) - yc(cell_2))**2.0)
     IF(threed) xtemp2 = SQRT((xf(iface) - xc(cell_2))**2.0 + &
                              (yf(iface) - yc(cell_2))**2.0 + &
                              (zf(iface) - zc(cell_2))**2.0)
     wcf(iface) = xtemp2/(xtemp1 + xtemp2)
  END DO
  IF(debug_level > 0) WRITE(io_debug,*) &
                      "Finished calculating wcf"

! Compute interpolant from cell to vertexn
  ALLOCATE(wcv(ncells,nv_max1), stat=ierr) 
  ALLOCATE(sumweight(nnodes), stat=ierr)
  ALLOCATE(bnode(nnodes), stat=ierr)

  sumweight(:) = zero
  wcv(:,:) = zero
  DO i=1,ncells
     DO j=1,nvcell(i)
        inode = lcv(i,j)
        disntemp = SQRT((xv(inode) - xc(i))**2.0 + (yv(inode) - yc(i))**2.0)
        IF(threed) disntemp = SQRT((xv(inode) - xc(i))**2.0 + &
                                   (yv(inode) - yc(i))**2.0 + &
                                   (zv(inode) - zc(i))**2.0)
        sumweight(inode) = sumweight(inode) + one/disntemp
        wcv(i,j) = one/disntemp
     END DO
  END DO

    DO i=1,ncells
       DO j=1,nvcell(i)
          inode = lcv(i,j)
          wcv(i,j) = wcv(i,j)/sumweight(inode)
       END DO
    END DO
    IF(debug_level > 0) WRITE(io_debug,*) &
                        "Finished calculating wcv"
    DEALLOCATE(sumweight)

! Compute interpolant from boundary face to vertex
  ALLOCATE(wbfv(nbcfaces,nv_max2), stat=ierr) 
  ALLOCATE(sumwtf(nnodes), stat=ierr)
    
  sumwtf(:) = zero
  wbfv(:,:) = zero
  DO i=1,nbcfaces
     ifc = bf_to_f(i)
     DO j=1,nvface(ifc)
        inode = lfv(ifc,j)
        disntemp = SQRT((xv(inode) - xf(ifc))**2.0 + (yv(inode) - yf(ifc))**2.0)
        IF(threed) disntemp = SQRT((xv(inode) - xf(ifc))**2.0 + &
                                   (yv(inode) - yf(ifc))**2.0 + &
                                   (zv(inode) - zf(ifc))**2.0)
        sumwtf(inode) = sumwtf(inode) + one/disntemp
        wbfv(i,j) = one/disntemp
     END DO
  END DO

  DO i=1,nbcfaces
     ifc = bf_to_f(i)
     DO j=1,nvface(ifc)
        inode = lfv(ifc,j)
        wbfv(i,j) = wbfv(i,j)/sumwtf(inode)
     END DO
  END DO
  IF(debug_level > 0) WRITE(io_debug,*) &
                      "Finished calculating weights"
  DEALLOCATE(sumwtf)  




! Compute interpolant from cell centre to face
 
!  ALLOCATE(wcf(ncells,nf_max), stat=ierr) 
!  ALLOCATE(sumweight(nfaces), stat=ierr)
!!  ALLOCATE(bnode(nnodes), stat=ierr)
!
!  sumweight(:) = zero
!  wcf(:,:) = zero
!  DO i=1,ncells
!     DO j=1,nfcell(i)
!        iface = lcf(i,j)
!!        inode = lcv(i,j)
!        disntemp = SQRT((xf(iface) - xc(i))**2.0 + (yf(iface) - yc(i))**2.0)
!        IF(threed) disntemp = SQRT((xf(iface) - xc(i))**2.0 + &
!                                   (yf(iface) - yc(i))**2.0 + &
!                                   (zf(iface) - zc(i))**2.0)
!        sumweight(iface) = sumweight(iface) + one/disntemp
!        wcf(i,j) = one/disntemp
!     END DO
!  END DO
!
!    DO i=1,ncells
!       DO j=1,nfcell(i)
!          iface = lcf(i,j)
!          wcf(i,j) = wcf(i,j)/sumweight(iface)
!       END DO
!    END DO
!    IF(debug_level > 0) WRITE(io_debug,*) &
!                        "Finished calculating wcf"
!    DEALLOCATE(sumweight)
!!





! Dont know about this ; not there in ashraf part

!  DO i = 1,nbcfaces
!     bf_to_c(i) = lfc(bf_to_f(i),1)
!  END DO

  ! Tagging the boundary nodes
  bnode(:) = 0
  DO iface=1,nfaces
    IF (bface(iface) == 1) THEN
      DO j=1,nvface(iface)
        inode = lfv(iface,j)
        bnode(inode) = 1
      END DO
    END IF
  END DO  
                   

! Load temperature into boundary condition array
  ALLOCATE(temp_bc(nbcfaces))
  temp_bc(:) = -100.d0
  DO ibc = 1,nbcfaces
    temp_bc(ibc) = phibc(ibc)
  ENDDO
  DEALLOCATE(phibc)
  
  CLOSE(unit = io_ace)
    
  IF(debug_level > 0) WRITE(io_debug,*) "Exiting grid_and_boundary" 
    
  END SUBROUTINE GRID_paree
