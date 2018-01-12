    Subroutine disord(pol,band,mfp)
    Use precisions
    Use constants !,Only: zero,pi,one,three,half,two,fourth,four,tiny
    Use grid !,Only:ncells,nfcell,lcf,areaf,volcell,vecfx,vecfy, &
    !vecfz,bface,f_to_bf,lfc,xc,yc,xf,yf,temp_bc,bctype, &
    !bcname,ADIA,ISOT,lcc,normdir
    Use variables
    Implicit None

    Integer(int_p) ::tempindex(7)
    Integer(int_p) :: i,j,currf,currcell,k,istart,ibface,ic
    Integer(int_p) :: x, index_old,kspace,ierr,icell
    Integer(int_p) :: si,l,bcon,sout,idir

    Real(real_p):: tempmat(7)
    Real(real_p):: sdotn,gamma,vol,insdotn,sdotn1,insdotn1
    Real(real_p):: mat_old,alpha,vomega,vel
    Real(real_p):: rmuval,retval,inrmuval,inrxival, temper, tempwall
    Real(real_p):: rxival,actin,actout,tempnum,beta,jout,inretval

   Integer(int_p),intent(in):: band, pol
    REAL(real_p), intent(in) :: mfp

    IF(debug_level > 0) WRITE(io_debug,*) "Starting bsnord_disord"

!!    Do iband = 1, nbands    !!$ loop over frequency bands
!!        pol = polar(iband)
!!        vel = gpvel(iband)
        Do si = 1, ndir     !!$ loop over the directions
            vomega = omega(si)
!            rmuval = rmu(si)
!            rxival = rxi(si)
!            retval = ret(si)
!            inrmuval = inrmu(si)
!            inrxival = inrxi(si)
!            inretval = inret(si)
            aps(:) = zero
            anbs(:,:) = zero
            scs(:) = zero
            resid(:) = zero
            r(:) = zero
            soln(:) = zero
            Do i = 1,ncells     !!$ loop over cells
                vol = volcell(i)*vomega
                temper = tnot(i)
                Call relaxtime(pol, band, temper, beta)
!!                mfp = vel/beta
                Do j = 1,nfcell(i)
                    currf = lcf(i,j)
                    gamma = mfp*areaf(currf)/vol
                    
                    ! chagned this
                    
                 IF (threed) THEN       
                     sdotn = vecfx(currf) * rmu(si) + vecfy(currf) * rxi(si) + vecfz(currf)*ret(si)        
                     insdotn = vecfx(currf) * inrmu(si) + vecfy(currf) * inrxi(si)+ vecfz(currf)*inret(si)
!                     dist1 = sqrt((xc(i)-xf(currf))**2 + (yc(i)-yf(currf))**2 + (zc(i)-zf(currf))**2)
                 ELSE 
                     sdotn = vecfx(currf) * rmu(si) + vecfy(currf) * rxi(si)     
                     insdotn = vecfx(currf) * inrmu(si) + vecfy(currf) * inrxi(si)
!                     dist1 = sqrt((xc(i)-xf(currf))**2 + (yc(i)-yf(currf))**2) 
                 ENDIF
                    
                    
                    
!                    sdotn = vecfx(currf) * rmuval + vecfy(currf) * rxival + vecfz(currf) * retval
!                    insdotn = vecfx(currf) * inrmuval + vecfy(currf) * inrxival + vecfz(currf) * inretval
                    If ( bface(currf) == 0 ) Then ! Interior face
                        alpha = normdir(i,j)
                        sdotn = alpha*sdotn
                        insdotn = alpha*insdotn
                        aps(i) = aps(i) + Max(0.0,sdotn)*gamma*((insdotn/(sdotn+tiny)))
                        anbs(i,j) = - Max(0.0,-sdotn)*gamma*((insdotn/(sdotn+tiny)))
                    Else !Boundary face
                        ibface = f_to_bf(currf)
                        bcon = bctype(ibface)
                        tempwall = temp_bc(ibface)
                        Call iwalls(pol,band,tempwall,actin)
                        aps(i) = aps(i) + Max(0.0,sdotn)*gamma*((insdotn/(sdotn+tiny)))
                        jout = zero
                        Select Case(bcon)
                        Case(ADIA)
                        
                         Do idir = 1,ndir

                                 IF (threed) THEN       
                                     sdotn1 = vecfx(currf) * rmu(idir) + vecfy(currf) * rxi(idir) + vecfz(currf)*ret(idir)        
                                     insdotn1 = vecfx(currf) * inrmu(idir) + vecfy(currf) * inrxi(idir)+ vecfz(currf)*inret(idir)
!                                     dist1 = sqrt((xc(i)-xf(currf))**2 + (yc(i)-yf(currf))**2 + (zc(i)-zf(currf))**2)
                                 ELSE 
                                     sdotn1 = vecfx(currf) * rmu(idir) + vecfy(currf) * rxi(idir)     
                                     insdotn1 = vecfx(currf) * inrmu(idir) + vecfy(currf) * inrxi(idir)
!                                     dist1 = sqrt((xc(i)-xf(currf))**2 + (yc(i)-yf(currf))**2) 
                                 ENDIF
                 

                                jout = jout + Max(0.0d0,sdotn1)*intensity(i,idir,band)*(insdotn1/ &
                                           sdotn1)
   
                            Enddo
                            Call reflec(si,ibface,sdotn,sout)
                            scs(i) = scs(i) + dos*intensity(i,sout,band) * gamma * &
                                     Max(0.0,-sdotn)*((insdotn/(sdotn+tiny))) &
                                     +((one-dos)*jout) * gamma * Max(0.0,-sdotn)* &
                                    ((insdotn/(sdotn+tiny)))/pi
                        Case(ISOT)
                            scs(i) = scs(i) + actin * gamma * Max(0.0,-sdotn)*((insdotn/(sdotn+tiny)))
                        end select

                    End If !Bounadry or Interior Face

                End Do ! Faces of cell

                Call iwalls(pol,band,temper,actout)
                aps(i) = aps(i) + idt/beta + one
                scs(i) = scs(i) + actout + idt *intensityone(i,si,band)/beta
                Do j = 1,nfcell(i)
                    currf = lcf(i,j)
                    If ( bface(currf) == 0 ) Then
                        currcell = lcc(i,j)
                        resid(i) = resid(i) - anbs(i,j) * intensity(currcell,si,band)
                    End If
                End Do
                ! print *, i, aps(i)
                resid(i) = resid(i) - aps(i) * intensity(i,si,band)
                !!$ Normalize
                anbs(i,:) = anbs(i,:) / aps(i)
                scs(i) = scs(i)/aps(i)
                resid(i) = resid(i) / aps(i)
                aps(i) = one

            End Do !End Cell Loop

            ! end of populating the matrices
            ! you might want to normalize your matrices with the diagonal here
            resid(:) = resid(:) + scs(:)
            r(:) = resid(:) ! We need r to be passed to gmres solver.
            am(:) = zero
            ia(:) = zero
            ja(:) = zero
            k = 0
            istart = 1
            Do i = 1,ncells
                ia(istart) = k + 1
                tempmat = zero
                tempindex = 0
                l = 1
                tempmat(l) = aps(i)
                tempindex(l) = i
                Do j = 1,nfcell(i)
                    currcell = lcc(i,j)
                    IF(ABS(anbs(i,j)) > tiny)THEN
                     l = l + 1
                      tempmat(l) = anbs(i,j)        !changes made on 29th may by dr mazumder
                      tempindex(l) = currcell
                    ENDIF
                ENDDO
                    
                Do j = 1,l-1
                    Do x = j+1,l
                        If( tempindex(x) < tempindex(j) ) Then
                            index_old = tempindex(x)
                            mat_old = tempmat(x)
                            tempmat(x) = tempmat(j)
                            tempindex(x) = tempindex(j)
                            tempmat(j) = mat_old
                            tempindex(j) = index_old
                        End If
                    End Do
                End Do
                am(k+1:k+l) = tempmat(1:l)
                ja(k+1:k+l) = tempindex(1:l)
                k = k + l
                istart = istart + 1
            End Do

            IA(istart) = k + 1
            alu(:) = zero
            jlu(:) = zero
            ju(:) = zero
            jr(:) = zero
            vv(:,:) = zero
            KSPACE = 20

            Call ilu0(ncells,AM,JA,IA,alu,jlu,ju,jr,ierr)
            if ( ierr .NE. 0) then
                WRITE(io_debug,*) "Error trying to do a LU decomposition. Exiting the program!"
                STOP
            ENDIF

            Call pgmres(ncells,kspace,r,soln,vv,tol_inner,iter_gmres,0,AM,JA,IA,alu,jlu,ju,ierr)

            Do icell = 1,ncells
                intensity(icell,si,band) = intensity(icell,si,band) + soln(icell)
            End Do

            ! WRITE(*,*) "Finished direction", si

        Enddo ! Directions

        !!$ calculate the integrated intensity of each spectral bin
        Do ic = 1, ncells
            tempnum = zero
            Do si = 1, ndir
                tempnum = tempnum + intensity(ic,si,band)*omega(si)
            End Do
            gw(ic,band) = tempnum
        End Do

        ! write(*,*), "band done", iband

!!    End Do ! Bands

    IF(debug_level > 0) WRITE(io_debug,*) "Finishing bsnord disord"

    End Subroutine disord
