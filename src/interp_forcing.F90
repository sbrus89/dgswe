      SUBROUTINE interp_forcing()

      USE globals, ONLY: pres,nqpte,qpte, &
                         nope,obseg,nbou,fbseg, &
                         nobed,nfbed,nobfr,nfbfr, &
                         obamp,obph,fbamp,fbph, &
                         obamp_qpt,obph_qpt,obdepth_qpt,fbamp_qpt,fbph_qpt, &
                         depth,obnds
      USE allocation, ONLY: alloc_forcing_arrays 

      IMPLICIT NONE

      INTEGER :: pt,bfr,seg,k
      INTEGER :: nd,ed
      INTEGER :: segtype
      REAL(pres) :: xi,l1,l2
             

      ! Interpolate boundary forcing data !
      ! qpte(1:nqpte,2,1) are the 1-D edge quadrature points 
      
      CALL alloc_forcing_arrays(3)

      nd = 1
      DO pt = 1,nqpte(1)
        xi = qpte(pt,2,1)
        l1 = .5d0*(1d0-xi)
        l2 = .5d0*(1d0+xi)
        DO bfr = 1,nobfr
          ed = 1
          DO seg = 1,nope
            DO k = 1,obseg(seg)-1
              obamp_qpt(bfr,pt,ed) = obamp(nd,bfr)*l1 + obamp(nd+1,bfr)*l2
              obph_qpt(bfr,pt,ed)  =  obph(nd,bfr)*l1 +  obph(nd+1,bfr)*l2
              obdepth_qpt(ed,pt) = depth(obnds(k,seg))*l1 + depth(obnds(k+1,seg))*l2
              ed = ed + 1
              nd = nd + 1
            ENDDO
            nd = nd + 1
          ENDDO
          nd = 1
        ENDDO
      ENDDO

      nd = 1
      DO pt = 1,nqpte(1)
        xi = qpte(pt,2,1)
        l1 = .5d0*(1d0-xi)
        l2 = .5d0*(1d0+xi)
        DO bfr = 1,nfbfr
          ed = 1
          DO seg = 1,nbou
            segtype = fbseg(2,seg)
            IF(segtype == 2 .OR. segtype == 12 .OR. segtype == 22)THEN
              DO k = 1,fbseg(1,seg)-1
                fbamp_qpt(bfr,pt,ed) = fbamp(nd,bfr)*l1 + fbamp(nd+1,bfr)*l2
                fbph_qpt(bfr,pt,ed)  =  fbph(nd,bfr)*l1 + fbph(nd+1,bfr)*l2
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