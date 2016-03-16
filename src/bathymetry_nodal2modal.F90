      SUBROUTINE bathymetry_nodal2modal()
      
      USE globals, ONLY: rp,ne,np,norder,order,mnnds,nnds,el_type, &
                         hbm,elhb
      USE read_dginp, ONLY: out_direc
      USE vandermonde, ONLY: vandermonde_area
      USE lapack_interfaces
      
      IMPLICIT NONE
     
      INTEGER :: et,el,nd,dof
      INTEGER :: p,n,typ,eo
      INTEGER :: info
      REAL(rp) :: V(mnnds,mnnds,norder)
      INTEGER :: ipiv(mnnds,norder)     
     
     
      DO et = 1,norder
        p = np(et)
        CALL vandermonde_area(et,p,n,V(:,:,et))
        CALL DGETRF(n,n,V(1,1,et),mnnds,ipiv(1,et),info)  
      ENDDO
     
      ALLOCATE(hbm(mnnds,ne))     
      hbm = 0d0 
      DO el = 1,ne
      
        et = el_type(el)
        typ = et + 4
        eo = order(typ) 
        n = nnds(eo)    
                
        DO nd = 1,n
          hbm(nd,el) = elhb(nd,el)
        ENDDO
        
        CALL DGETRS("T",n,1,V(1,1,eo),mnnds,ipiv(1,eo),hbm(1,el),mnnds,info)                         
        
      ENDDO
      
      OPEN(unit=65,FILE=trim(out_direc) //'hb_modal.d')
      DO dof = 1,mnnds
        WRITE(65,"(16000(e24.17,1x))") (hbm(dof,el), el = 1,ne)
      ENDDO
      CLOSE(65)     
     
     
      RETURN
      END SUBROUTINE bathymetry_nodal2modal