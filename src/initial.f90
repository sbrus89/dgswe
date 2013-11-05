      SUBROUTINE initial()

      USE globals, ONLY: pres,ne,nn,ndof, &
                         ect,xy,depth, &
                         nqpta,qpta,wpta,nqpte,qpte, &
                         Hinit,Qxinit,Qyinit, &
                         phia,phil, &
                         nope,obseg,nbou,fbseg, &
                         nobed,nfbed,nobfr,nfbfr, &
                         obamp,obph,fbamp,fbph, &
                         obamp_qpt,obph_qpt,obdepth_qpt,fbamp_qpt,fbph_qpt, &
                         depth,obnds,grid_file

      IMPLICIT NONE
      INTEGER :: i,el,l,pt,dof,bfr,ed,ind,nd,k,seg
      INTEGER :: segtype
      INTEGER :: alloc_stat
      REAL(pres) :: r,s,x1,x2,x3,y1,y2,y3,x,y
      REAL(pres) :: sigma,xc,yc,h0
      REAL(pres) :: qint,f

      sigma = 10d0
      xc = 0d0 
      yc = 0d0

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! Compute initial condition
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


      DO el = 1,ne
        DO l = 1,ndof
          qint = 0d0
          DO pt = 1,nqpta
            r = qpta(pt,1)
            s = qpta(pt,2)

            x1 = xy(1,ect(1,el))
            x2 = xy(1,ect(2,el))
            x3 = xy(1,ect(3,el))

            y1 = xy(2,ect(1,el))
            y2 = xy(2,ect(2,el))
            y3 = xy(2,ect(3,el))

            x = -.5d0*((r+s)*x1 - (1d0+r)*x2 - (1d0+s)*x3)
            y = -.5d0*((r+s)*y1 - (1d0+r)*y2 - (1d0+s)*y3)

            h0 = depth(ect(1,el))*phil(1,pt) + depth(ect(2,el))*phil(2,pt) + depth(ect(3,el))*phil(3,pt)

!             f = exp(-sigma*((x-xc)**2d0+(y-yc)**2d0))+h0 
!               f = .0002*x + h0
            f = h0

!             qint = qint + wpta(pt)*f*phia(l,pt)
            qint = qint + .5d0*wpta(pt)*f*phia(l,pt)
          ENDDO
          Hinit(el,l) = qint
        ENDDO
      ENDDO

      Qxinit(:,:) = 0d0
      Qyinit(:,:) = 0d0

      ! Write initial condition
      WRITE(63,"(A)") grid_file
      WRITE(63,"(e24.17)") 0d0
      DO dof = 1,ndof
        WRITE(63,"(16000(e24.17,1x))") (Hinit(el,dof),el=1,ne)
      ENDDO

      WRITE(641,"(A)") grid_file
      WRITE(641,"(e24.17)") 0d0
      DO dof = 1,ndof
        WRITE(641,"(16000(e24.17,1x))") (Qxinit(el,dof),el=1,ne)
      ENDDO

      WRITE(642,"(A)") grid_file
      WRITE(642,"(e24.17)") 0d0
      DO dof = 1,ndof
        WRITE(642,"(16000(e24.17,1x))") (Qyinit(el,dof),el=1,ne)
      ENDDO

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! Interpolate boundary forcing data
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      ALLOCATE(obamp_qpt(nqpte*nobfr,nobed),obph_qpt(nqpte*nobfr,nobed),obdepth_qpt(nobed,nqpte),STAT=alloc_stat)
      IF(alloc_stat /= 0) THEN
        PRINT*, 'Allocation error: obamp_qpt,obph_qpt'
      ENDIF

      nd = 1
      DO pt = 1,nqpte
        DO bfr = 1,nobfr
          ed = 1
          ind = (pt-1)*nobfr + bfr
          DO seg = 1,nope
            DO k = 1,obseg(seg)-1
              obamp_qpt(ind,ed) = .5d0*( obamp(nd,bfr)*(1d0-qpte(pt)) + obamp(nd+1,bfr)*(1d0+qpte(pt)) )
              obdepth_qpt(ed,pt) =  .5d0*( depth(obnds(k,seg))*(1d0-qpte(pt)) + depth(obnds(k+1,seg))*(1d0+qpte(pt)) )
              obph_qpt(ind,ed) =  .5d0*( obph(nd,bfr)*(1d0-qpte(pt)) + obph(nd+1,bfr)*(1d0+qpte(pt)) )
              ed = ed + 1
              nd = nd + 1
            ENDDO
            nd = nd + 1
          ENDDO
          nd = 1
        ENDDO
      ENDDO

      ALLOCATE(fbamp_qpt(nqpte*nfbfr,nfbed),fbph_qpt(nqpte*nfbfr,nfbed),STAT=alloc_stat)
      IF(alloc_stat /= 0) THEN
        PRINT*, 'Allocation error: fbamp_qpt,fbph_qpt'
      ENDIF

      nd = 1
      DO pt = 1,nqpte
        DO bfr = 1,nfbfr
          ed = 1
          DO seg = 1,nbou
            segtype = fbseg(2,seg)
            IF(segtype == 2 .OR. segtype == 12 .OR. segtype == 22)THEN
              ind = (pt-1)*nfbfr + bfr 
              DO k = 1,fbseg(1,seg)-1
                fbamp_qpt(ind,ed) = .5d0*( fbamp(nd,bfr)+fbamp(nd+1,bfr) + (fbamp(nd+1,bfr)-fbamp(nd,bfr))*qpte(pt) ) 
                fbph_qpt(ind,ed) =  .5d0*( fbph(nd,bfr)+fbph(nd+1,bfr) + (fbph(nd+1,bfr)-fbph(nd,bfr))*qpte(pt) ) 
                nd = nd + 1
                ed = ed + 1
              ENDDO
              nd = nd + 1
            ENDIF
          ENDDO
          nd = 1
        ENDDO
      ENDDO

      RETURN
      END SUBROUTINE
