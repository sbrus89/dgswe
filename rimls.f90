      PROGRAM rimls

      USE globals
      USE kdtree2_module
      USE evaluate
      USE basis

      IMPLICIT NONE

      INTEGER :: el,nd,pt,i,j,ed
      INTEGER :: mninds
      
      Erad = 6378206.4d0 
      
!      ! beaufort      
!       lambda0 = -76d0*deg2rad
!       phi0 = 33d0*deg2rad

!       ! shinnecock
!       lambda0 = -72.432511d0*deg2rad
!       phi0 = 40.666091d0*deg2rad
      
      ! EC2001
      lambda0 = -79d0*deg2rad
      phi0 = 35d0*deg2rad      
      
!      ! Cartesian      
!       Erad = 1d0
!       lambda0 = 0d0
!       phi0 = 0d0

      
      CALL read_input()            
      
      ndof(1) = (p+1)*(p+2)/2
      ndof(2) = (p+1)*(p+1)
      ndof(3) = ndof(1)
      ndof(4) = ndof(2)
      mndof = maxval(ndof)
      
      np(1) = 1
      np(2) = 1
      np(3) = ctp
      np(4) = ctp        
      
      nnds(1) = 3
      nnds(2) = 4
      nnds(3) = (ctp+1)*(ctp+2)/2
      nnds(4) = (ctp+1)*(ctp+1) 
      mnnds = maxval(nnds)    
      
      nverts(1) = 3
      nverts(2) = 4
      nverts(3) = 3
      nverts(4) = 4
      

      
      CALL read_grid(base)
      CALL read_grid(fine)
      
      CALL connect(base) 
      CALL connect(fine)       

      
      CALL vandermonde()  
      
      CALL transformation()      
      
      CALL normals(base)
      
      CALL coordinates(fine)
      
      ALLOCATE(xyhw(mnnds,fine%ned,3))
      
      mninds = nnds(3)-3*(np(3)-1)-3
      
!       CALL invcpp(ne,mninds,mnnds,xyhi,xyhw)
      
      OPEN(unit=9,file='interior_nodes.d')
      WRITE(9,*) fine%ne,mninds
      DO el = 1,fine%ne
        DO i = 1,mninds
!           WRITE(9,*) (xyhw(i,el,j), j = 1,3)
          WRITE(9,*) fine%xyhi(i,el,1)/(Erad*cos(phi0))+lambda0, fine%xyhi(i,el,2)/Erad, fine%xyhi(i,el,3)
        ENDDO
      ENDDO
      
      CLOSE(9)
      
!       CALL invcpp(ned,np(3)-1,mnnds,xyhe,xyhw)
      
      OPEN(unit=10,file='edge_nodes.d')
      WRITE(10,*) fine%ned,np(3)-1
      DO ed = 1,fine%ned      
        DO i = 1,np(3)-1
!           WRITE(10,*) (xyhw(i,ed,j), j = 1,3)
          WRITE(10,*) fine%xyhe(i,ed,1)/(Erad*cos(phi0))+lambda0, fine%xyhe(i,ed,2)/Erad, fine%xyhe(i,ed,3)          
        ENDDO
      ENDDO
      
      CLOSE(10)
      
      OPEN(unit=11,file='centers.d')
      WRITE(11,*) base%ne    
      DO el = 1,base%ne      
        WRITE(11,*) (base%xyhc(j,el), j = 1,3)
      ENDDO      
      
      CLOSE(11)
      
      OPEN(unit=12,file='normals.d')
      WRITE(12,*) base%ne
      DO el = 1,base%ne
        WRITE(12,*) (base%nhb(i,el),i=1,3)
      ENDDO
      
      CLOSE(12)
      
      
      OPEN(unit=13,file='boundary_nodes.d')
      WRITE(13,*) fine%ned,np(3)-1
      DO ed = 1,fine%ned      
        DO i = 1,np(3)-1
          WRITE(13,*) fine%bnd_flag(i,ed)        
        ENDDO
      ENDDO
      
      CLOSE(13)      
      
      
      ! Build kd-tree           
!       tree_xy => kdtree2_create(vxy , rearrange=.true., sort=.true.)
      tree_xy => kdtree2_create(base%xy  , rearrange=.true., sort=.true.)
      tree_c  => kdtree2_create(base%xyhc, rearrange=.true., sort=.true.)
      
      ALLOCATE(kdresults(base%ne))       
!       
      CALL grid_size(base)
      
      PRINT("(A)"), "Computing rimls surface: verticies"
      CALL rimls_surface(fine%nn,1,1,fine%xyhv)      
      PRINT("(A)"), "Computing rimls surface: edges"      
      CALL rimls_surface(fine%ned,np(3)-1,mnnds,fine%xyhe)
      PRINT("(A)"), "Computing rimls surface: interior"
      CALL rimls_surface(fine%ne,nnds(3)-3*(np(3)-1)-3,mnnds,fine%xyhi)

      
      
      OPEN(unit=9,file='rimls_interior_nodes.d')
      WRITE(9,*) fine%ne,mninds
      DO el = 1,fine%ne
        DO i = 1,mninds
          WRITE(9,*) (fine%xyhi(i,el,j), j = 1,3)
        ENDDO
      ENDDO
      
      CLOSE(9)      
      
      
      OPEN(unit=10,file='rimls_edge_nodes.d')
      WRITE(10,*) fine%ned,np(3)-1
      DO ed = 1,fine%ned      
        DO i = 1,np(3)-1
          WRITE(10,*) (fine%xyhe(i,ed,j), j = 1,3)
        ENDDO
      ENDDO
      
      CLOSE(10)
      
      OPEN(unit=11,file='rimls_vertex_nodes.d')
      WRITE(11,*) fine%nn,1
      DO nd = 1,fine%nn
          WRITE(11,*) (fine%xyhv(1,nd,j), j = 1,3)
      ENDDO
      
      CLOSE(11)      
      
     
      END PROGRAM rimls