SUBROUTINE fluxt

  USE PRECISIONS
  USE CONSTANTS !, only : zero,two
  USE VARIABLES !,only : jt,phiv
  USE GRID
  
!  USE GRID, only :  nfaces,lfc,lfv,xf,yf,zf,xc,yc,zc,vecfx,vecfy,vecfz &
!                        ,nvface,lfv,areaf,xv,yv,zv,


  ! Jt is calcuated such tha, tangent goes from lfv1 to lfv2  and l goes from
  ! lfc1 to lfc2 . therefore for lfc2 direction will be reversed
  IMPLICIT NONE

  !REAL(real_p),DIMENSION(:),POINTER :: jt
  INTEGER(int_p) :: i,j,v1,v2,ic1,ic2,iband
  REAL(real_p) :: phiside,sum, xv1,xv2,yv1,yv2
  REAL(real_p),DIMENSION(3) :: p1,p2,ncl,t,l
  
!  IF(debug_level > 0) WRITE(io_debug,*) "Starting fluxt" 
  IF (threed) THEN
 !             do iband = 2,2
             ! do iband = 14,20
 do iband = 1,nbands
      DO i = 1,nfaces
         ic1 = lfc(i,1)
         ic2 = lfc(i,2)
         IF( ic1 == ic2 ) THEN
            l(1) = xf(i) - xc(ic1)
            l(2) = yf(i) - yc(ic1)
            l(3) = zf(i) - zc(ic1)
         ELSE
            l(1) = xc(ic2) - xc(ic1)
            l(2) = yc(ic2) - yc(ic1)
            l(3) = zc(ic2) - zc(ic1)
         END IF
         ncl(1) = ( vecfy(i) * l(3) - vecfz(i) * l(2) )
         ncl(2) = ( vecfz(i) * l(1) - vecfx(i) * l(3) )
         ncl(3) = ( vecfx(i) * l(2) - vecfy(i) * l(1) )
         sum = zero
         DO j = 1,nvface(i)
            IF( j == nvface(i) ) THEN
               v1 = lfv(i,j)
               v2 = lfv(i,1)
            ELSE
               v1 = lfv(i,j)
               v2 = lfv(i,j+1)
            END IF
            p1(1) = xv(v1)
            p1(2) = yv(v1)
            p1(3) = zv(v1)
            p2(1) = xv(v2)
            p2(2) = yv(v2)
            p2(3) = zv(v2)
            t = p2 - p1
            phiside = ( phiv(v1,iband) + phiv(v2,iband) ) / two
    !!$        le = SQRT( ( p2(1) - p1(1) )**two + ( p2(2) - p1(2) )**two + ( p2(3) - p1(3) )**two )
            sum = sum + phiside * ( t(1) * ncl(1) + t(2) * ncl(2) + t(3) * ncl(3) )
         END DO
         jt(i,iband) = sum / areaf(i)
         WRITE(28,*) "jt(i,iband)", jt(i,iband)
      END DO
  end do

  ELSE      !2d part
  
  
  !! start
  
  
!    Do i = 1,nfaces
!       node_1 = lfv(i,1)
!       node_2 = lfv(i,2)
!       tgtx(i) = ( xv(node_2) - xv(node_1) ) / ( Sqrt( (xv(node_2) - xv(node_1))**2 + (yv(node_2) - yv(node_1))**2 ) )
!       tgty(i) = ( yv(node_2) - yv(node_1) ) / ( Sqrt( (xv(node_2) - xv(node_1))**2 + (yv(node_2) - yv(node_1))**2 ) )
!    End Do

!    Do iface=1,nfaces
!
!       cell_1 =  lfc(iface,1)
!       cell_2 =  lfc(iface,2)  
!       If(cell_1 .Eq. cell_2) Then
!          disn(iface) = Abs((xf(iface)-xc(cell_1))*vecfx(iface) + (yf(iface)-yc(cell_1))*vecfy(iface))
!          tdotlf(iface) = ( tgtx(iface) * ( xf(iface) - xc(cell_1) ) + tgty(iface) * ( yf(iface) - yc(cell_1) ) )          
!          !disn(iface) = ABS(SQRT((xf(iface) - xc(cell_1))**2.0 + (yf(iface) - yc(cell_1))**2.0))
!       Else
!          disn(iface) = Abs((xc(cell_2)-xc(cell_1))*vecfx(iface) + (yc(cell_2)-yc(cell_1))*vecfy(iface))
!          tdotlf(iface) = ( tgtx(iface) * ( xc(cell_2) - xc(cell_1) ) + tgty(iface) * ( yc(cell_2) - yc(cell_1) ) )
!          !disn(iface) = ABS(SQRT((xc(cell_1) - xc(cell_2))**2.0 + (yc(cell_1) - yc(cell_2))**2.0))
!       End If
!    End Do
 !             do iband = 2,2
             ! do iband = 14,20
 do iband = 1,nbands
       Do i = 1,nfaces
            v1 = lfv(i,1)
            v2 = lfv(i,2)
            xv1 = xv(v1)
            yv1 = yv(v1)
            xv2 = xv(v2)
            yv2 = yv(v2)
            jt(i,iband) = ( phiv(v2,iband) - phiv(v1,iband) ) * tdotlf(i) / Sqrt(( xv2 - xv1 )**2 + ( yv2 - yv1 )**2)
            !WRITE(28,*) "jt(i,iband)", jt(i,iband)
       End Do
   End do    
  ENDIF 
  
  
  !!finish
  
  
  
  
  
  
  
  
!  
!  
!  
!      DO i = 1,nfaces
!         ic1 = lfc(i,1)
!         ic2 = lfc(i,2)
!         IF( ic1 == ic2 ) THEN
!            l(1) = xf(i) - xc(ic1)
!            l(2) = yf(i) - yc(ic1)
!         ELSE
!            l(1) = xc(ic2) - xc(ic1)
!            l(2) = yc(ic2) - yc(ic1)
!         END IF
!         ncl(1) = ( vecfy(i) * l(3) - vecfz(i) * l(2) )
!         ncl(2) = ( vecfz(i) * l(1) - vecfx(i) * l(3) )
!         ncl(3) = ( vecfx(i) * l(2) - vecfy(i) * l(1) )
!         sum = zero
!         DO j = 1,nvface(i)
!            IF( j == nvface(i) ) THEN
!               v1 = lfv(i,j)
!               v2 = lfv(i,1)
!            ELSE
!               v1 = lfv(i,j)
!               v2 = lfv(i,j+1)
!            END IF
!            p1(1) = xv(v1)
!            p1(2) = yv(v1)
!            p1(3) = zv(v1)
!            p2(1) = xv(v2)
!            p2(2) = yv(v2)
!            p2(3) = zv(v2)
!            t = p2 - p1
!            phiside = ( phiv(v1) + phiv(v2) ) / two
!    !!$        le = SQRT( ( p2(1) - p1(1) )**two + ( p2(2) - p1(2) )**two + ( p2(3) - p1(3) )**two )
!            sum = sum + phiside * ( t(1) * ncl(1) + t(2) * ncl(2) + t(3) * ncl(3) )
!         END DO
!         jt(i) = sum / areaf(i)
!      END DO
!  ENDIF
  
! IF(debug_level > 0) WRITE(io_debug,*) "Exiting fluxt" 

END SUBROUTINE fluxt


