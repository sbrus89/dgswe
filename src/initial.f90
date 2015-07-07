      SUBROUTINE initial()

      USE globals, ONLY: pres,ne,nn,mndof,ndof,nel_type,mnnds,nnds,np, &
                         ect,xy,depth,elxy,elhb,hbnodes, &
                         nqpta,qpta,wpta,nqpte,qpte, &
                         H,Qx,Qy, &
                         Zinit,Hinit,Qxinit,Qyinit, &
                         phia,phil, &
                         nope,obseg,nbou,fbseg, &
                         nobed,nfbed,nobfr,nfbfr, &
                         obamp,obph,fbamp,fbph, &
                         obamp_qpt,obph_qpt,obdepth_qpt,fbamp_qpt,fbph_qpt, &
                         depth,obnds, &
                         el_type,mnelnds,nelnds,psia, &
                         detJa,mmi_init
                         
      USE allocation, ONLY: alloc_sol_arrays,alloc_forcing_arrays   
      USE basis, ONLY: tri_nodes,tri_basis,quad_nodes,quad_basis
      USE read_dginp, ONLY: out_direc,grid_file

      IMPLICIT NONE
      INTEGER :: i,j,el,l,pt,dof,bfr,ed,ind,nd,k,seg,m,n,p
      INTEGER :: segtype,et
      INTEGER :: mnqpte
      INTEGER :: ipiv(mnnds,nel_type),info
      REAL(pres) :: x1,x2,x3,y1,y2,y3,x,y
      REAL(pres) :: sigma,xc,yc,h0,h00
      REAL(pres) :: qint,f,qint2
      REAL(pres) :: xn,yn
      REAL(pres),ALLOCATABLE :: rhsH(:,:),rhsH2(:,:)      
      REAL(pres) :: r(mnnds),s(mnnds),phi(mnnds*mnnds),hb(mnnds)
      REAL(pres) :: V(mnnds,mnnds,nel_type)

      sigma = 10d0
      xc = 0d0 
      yc = 0d0
      
      CALL alloc_sol_arrays()

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! Compute initial condition
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      ALLOCATE(rhsH(ne,mndof))
      ALLOCATE(rhsH2(ne,mndof))
      
      Hinit = 0d0
      
!       DO el = 1,ne
!         et = el_type(el)
!         DO l = 1,ndof(et)
!           qint = 0d0
!           qint2 = 0d0
!           DO pt = 1,nqpta(et)
!           
!             r = qpta(pt,1,et)
!             s = qpta(pt,2,et)
!             
!             x = 0d0
!             y = 0d0
!             DO nd = 1,nelnds(el)
! 
!               xn = elxy(nd,el,1)
!               yn = elxy(nd,el,2)            
!               
!               x = x + psia(nd,pt,et)*xn
!               y = y + psia(nd,pt,et)*yn
!             ENDDO
!             
! !             x1 = xy(1,ect(1,el))
! !             x2 = xy(1,ect(2,el))
! !             x3 = xy(1,ect(3,el))
! ! 
! !             y1 = xy(2,ect(1,el))
! !             y2 = xy(2,ect(2,el))
! !             y3 = xy(2,ect(3,el))
! 
! !             x = -.5d0*((r+s)*x1 - (1d0+r)*x2 - (1d0+s)*x3)
! !             y = -.5d0*((r+s)*y1 - (1d0+r)*y2 - (1d0+s)*y3)
! 
!             h0 = 0d0
!             DO nd = 1,nelnds(el)   ! This assumes there is an equal order representation between the bathymetry and the coordinate transformation
!               h0 = h0 + psia(nd,pt,et)*elhb(nd,el)
!             ENDDO
! 
! !             h00 = depth(ect(1,el))*phil(1,pt,1) + depth(ect(2,el))*phil(2,pt,1) + depth(ect(3,el))*phil(3,pt,1)
! 
! 
! !             f = exp(-sigma*((x-xc)**2d0+(y-yc)**2d0))+h0 
! !               f = .0002*x + h0
!             f = h0
!             
!             qint = qint + wpta(pt,et)*f*phia(l,pt,et)*detJa(el,pt)
!             qint2 = qint2 + 0.5d0*wpta(pt,et)*f*phia(l,pt,et)
!           ENDDO
!           rhsH(el,l) = qint
!           rhsH2(el,l) = qint2
! !           Hinit(el,1) = el ! for debugging
!         ENDDO
!         
!         m = 1
!         DO i = 1,ndof(et)        
!           DO j = 1,ndof(et)
!             Hinit(el,i) = Hinit(el,i) + mmi_init(el,m)*rhsH(el,j)
!             m = m + 1
!           ENDDO
!         ENDDO
!       ENDDO
      
      DO et = 1,nel_type
        n = nnds(et)
        p = np(et)
        IF (mod(et,2) == 1) THEN
          CALL tri_nodes(1,p,n,r,s)
          CALL tri_basis(p,n,n,r,s,phi)       
        ELSE IF (mod(et,2) == 0) THEN
          CALL quad_nodes(1,p,n,r,s)
          CALL quad_basis(p,n,n,r,s,phi)
        ENDIF
        
        DO pt = 1,n
          DO dof = 1,n
            i = (dof-1)*n + pt
            V(dof,pt,et) = phi(i)
          ENDDO
        ENDDO
        
        CALL DGETRF(n,n,V(1,1,et),mnnds,ipiv(1,et),info)        
!         DO pt = 1,n
!             PRINT("(100(e15.5))"), (V(dof,pt,et), dof = 1,n)
!         ENDDO        
!         PRINT*, " "
        
      ENDDO
            
      
      
!       DO el = 1,ne
!         et = el_type(el)
!         n = nnds(et)
!         
!         DO i = 1,n 
!           hb(i) = elhb(i,el)
!         ENDDO
!       
!         CALL DGETRS("T",n,1,V(1,1,et),mnnds,ipiv(1,et),hb,mnnds,info)
!         
!         DO i = 1,n
!           Hinit(el,i) = hb(i)
!         ENDDO
!       ENDDO
      
!       DO el = 1,ne
!         PRINT("(I5,3(e24.17),10x,3(e24.17))"), el, (rhsH2(el,l), l = 1,ndof(1)), (Hinit(el,l), l = 1,ndof(1))
!         PRINT("(I5,3(e24.17))"), el, (abs(rhsH2(el,l) - Hinit(el,l)), l = 1,ndof(1))
!       ENDDO



      Qxinit(:,:) = 0d0
      Qyinit(:,:) = 0d0
      Zinit(:,:) = 0d0
      

!       Write initial condition

      OPEN(unit=63,file=trim(out_direc) // 'solution_H.d')
      OPEN(unit=641,file=trim(out_direc) // 'solution_Qx.d')
      OPEN(unit=642,file=trim(out_direc) // 'solution_Qy.d')
      
      WRITE(63,"(A)") grid_file
      WRITE(63,"(e24.17)") 0d0
      DO dof = 1,mndof
!         WRITE(63,"(16000(e24.17,1x))") (Hinit(el,dof),el=1,ne)
        WRITE(63,"(16000(e24.17,1x))") (Zinit(el,dof),el=1,ne)
      ENDDO

      WRITE(641,"(A)") grid_file
      WRITE(641,"(e24.17)") 0d0
      DO dof = 1,mndof
        WRITE(641,"(16000(e24.17,1x))") (Qxinit(el,dof),el=1,ne)
      ENDDO

      WRITE(642,"(A)") grid_file
      WRITE(642,"(e24.17)") 0d0
      DO dof = 1,mndof
        WRITE(642,"(16000(e24.17,1x))") (Qyinit(el,dof),el=1,ne)
      ENDDO


      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! Interpolate boundary forcing data ( qpte(1:nqpte,2,1) are the 1-D edge quadrature points ) 
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      
      mnqpte = maxval(nqpte)

      CALL alloc_forcing_arrays(3)

      nd = 1
      DO pt = 1,nqpte(1)
        DO bfr = 1,nobfr
          ed = 1
          ind = (pt-1)*nobfr + bfr
          DO seg = 1,nope
            DO k = 1,obseg(seg)-1
              obamp_qpt(ind,ed) = .5d0*( obamp(nd,bfr)*(1d0-qpte(pt,2,1)) + obamp(nd+1,bfr)*(1d0+qpte(pt,2,1)) )
              obdepth_qpt(ed,pt) =  .5d0*( depth(obnds(k,seg))*(1d0-qpte(pt,2,1)) + depth(obnds(k+1,seg))*(1d0+qpte(pt,2,1)) )
              obph_qpt(ind,ed) =  .5d0*( obph(nd,bfr)*(1d0-qpte(pt,2,1)) + obph(nd+1,bfr)*(1d0+qpte(pt,2,1)) )
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
        DO bfr = 1,nfbfr
          ed = 1
          DO seg = 1,nbou
            segtype = fbseg(2,seg)
            IF(segtype == 2 .OR. segtype == 12 .OR. segtype == 22)THEN
              ind = (pt-1)*nfbfr + bfr 
              DO k = 1,fbseg(1,seg)-1
                fbamp_qpt(ind,ed) = .5d0*( fbamp(nd,bfr)+fbamp(nd+1,bfr) + (fbamp(nd+1,bfr)-fbamp(nd,bfr))*qpte(pt,2,1) ) 
                fbph_qpt(ind,ed) =  .5d0*( fbph(nd,bfr)+fbph(nd+1,bfr) + (fbph(nd+1,bfr)-fbph(nd,bfr))*qpte(pt,2,1) ) 
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
