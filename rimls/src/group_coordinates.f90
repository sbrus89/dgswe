      MODULE group_coordinates
      
      IMPLICIT NONE
      
      CONTAINS


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      SUBROUTINE separate_coordinates(mesh)
      
      USE globals, ONLY: rp,nverts,grid                                               
      
      IMPLICIT NONE
      
      TYPE(grid) :: mesh      
      
      INTEGER :: el,i,pt,led,j,ed
      INTEGER :: et,n,pn,nnd,nv,bed,edt
      
      REAL(rp) :: xpt,ypt,ytest      
           
           
      ALLOCATE(mesh%npts_vertex(mesh%nn),mesh%xyh_vertex(1,mesh%nn,3))       
      ALLOCATE(mesh%npts_edge(mesh%ned),mesh%xyh_edge(mesh%mnnds,mesh%ned,3))
      ALLOCATE(mesh%npts_interior(mesh%ne),mesh%xyh_interior(mesh%mnnds,mesh%ne,3))       
      ALLOCATE(mesh%bnd_flag(mesh%mnnds,mesh%ned))
     
      mesh%xyh_vertex = 0d0      
      mesh%xyh_edge = 0d0
      mesh%xyh_interior = 0d0
      mesh%bnd_flag = 0      
      
      PRINT "(A)", "Grouping vertex nodes..."
      
      DO i = 1,mesh%nn  
        mesh%xyh_vertex(1,i,1) = mesh%xy(1,i)
        mesh%xyh_vertex(1,i,2) = mesh%xy(2,i)
        mesh%xyh_vertex(1,i,3) = mesh%depth(i)
        mesh%npts_vertex(i) = 1
      ENDDO
      mesh%tpts_vertex = mesh%nn
      
      
      
      PRINT "(A)", "Grouping extra edge nodes..."        
      mesh%tpts_edge = 0
      DO ed = 1,mesh%ned
      
        el = mesh%ged2el(1,ed)
        led = mesh%ged2led(1,ed)
        et = mesh%el_type(el)
        nv = nverts(et)
        edt = mesh%ed_type(ed)
        
        IF (mod(et,2) == 1) THEN   
          pn = mesh%np(5)
        ELSE IF (mod(et,2) == 0) THEN
          pn = mesh%np(6)
        ENDIF    
        
        DO pt = 1,pn-1
          
          j = mod(led,nv)*pn + pt + 1          
          
          mesh%xyh_edge(pt,ed,1) = mesh%elxyh(j,el,1)
          mesh%xyh_edge(pt,ed,2) = mesh%elxyh(j,el,2)
          mesh%xyh_edge(pt,ed,3) = mesh%elhb(j,el)
          
          IF (edt /= 0) THEN
            mesh%bnd_flag(pt,ed) = 1
          ENDIF
          
        ENDDO      
        
        mesh%npts_edge(ed) = pn-1
        mesh%tpts_edge = mesh%tpts_edge + pn-1
        
!         IF (ed == 76) THEN
!           PRINT*, el
!           PRINT*, mesh%xyh_edge(1,ed,1),mesh%xyh_edge(1,ed,2)
!           STOP
!         ENDIF        
      ENDDO
      
      
      
      
      
      PRINT "(A)", "Grouping extra interior element nodes..."          
      mesh%tpts_interior = 0    
      DO el = 1,mesh%ne
      
        et = mesh%el_type(el)
        nv = nverts(et)
        
        IF (mod(et,2) == 1) THEN   
          pn = mesh%np(5)
          nnd = mesh%nnds(5)
        ELSE IF (mod(et,2) == 0) THEN
          pn = mesh%np(6)
          nnd = mesh%nnds(6)          
        ENDIF    

        pt = 0
        DO j = nv*(pn-1)+nv+1,nnd
        
          pt = pt + 1           
          
          mesh%xyh_interior(pt,el,1) = mesh%elxyh(j,el,1)
          mesh%xyh_interior(pt,el,2) = mesh%elxyh(j,el,2)
          mesh%xyh_interior(pt,el,3) = mesh%elhb(j,el)          
          
        ENDDO
        
        mesh%npts_interior(el) = pt
        mesh%tpts_interior = mesh%tpts_interior + pt        
      ENDDO
      
      RETURN
      END SUBROUTINE separate_coordinates
      
      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  

      SUBROUTINE gather_coordinates(mesh)
      
      USE globals, ONLY: grid,nverts
      
      IMPLICIT NONE
      
      TYPE(grid) :: mesh
      
      INTEGER :: el,ed,pt,i,v    
      INTEGER :: et,nv,led,n,nd,nnd,hbp

      hbp = mesh%hbp      
      
      mesh%elxyh = 0d0
      mesh%elhb = 0d0
      
      ! Verticies
      DO el = 1,mesh%ne
        et = mesh%el_type(el)
        nv = nverts(et)
        
        DO v = 1,nv
          i = (v-1)*hbp + 1
          nd = mesh%ect(v,el)
          
          mesh%elxyh(i,el,1) = mesh%xyh_vertex(1,nd,1)
          mesh%elxyh(i,el,2) = mesh%xyh_vertex(1,nd,2)
          mesh%elhb(i,el)    = mesh%xyh_vertex(1,nd,3)          
        ENDDO
      ENDDO
      
      ! Edge points
      DO ed = 1,mesh%ned
        
        el = mesh%ged2el(1,ed)
        led = mesh%ged2led(1,ed)
        et = mesh%el_type(el)
        nv = nverts(et)       
          
        DO pt = 1,hbp-1          
          i = mod(led,nv)*hbp + pt + 1

          mesh%elxyh(i,el,1) = mesh%xyh_edge(pt,ed,1)
          mesh%elxyh(i,el,2) = mesh%xyh_edge(pt,ed,2)          
          mesh%elhb(i,el)    = mesh%xyh_edge(pt,ed,3)          
        ENDDO          
          
        IF (mesh%ed_type(ed) == 0) THEN
        
          el = mesh%ged2el(2,ed)
          led = mesh%ged2led(2,ed)
          et = mesh%el_type(el)
          nv = nverts(et)                 
          
          DO pt = 1,hbp-1          
            i = mod(led,nv)*hbp + hbp - pt + 1     

            mesh%elxyh(i,el,1) = mesh%xyh_edge(pt,ed,1)
            mesh%elxyh(i,el,2) = mesh%xyh_edge(pt,ed,2)            
            mesh%elhb(i,el)    = mesh%xyh_edge(pt,ed,3)                       
          ENDDO            

        ENDIF   
        
!         PRINT*, " "
          
      ENDDO
      
      ! Interior points
      DO el = 1,mesh%ne
        et = mesh%el_type(el)
        n = mesh%nnds(et)
        nv = nverts(et)
        
        IF (mod(et,2) == 1) THEN   
          nnd = mesh%nnds(5)
        ELSE IF (mod(et,2) == 0) THEN
          nnd = mesh%nnds(6)          
        ENDIF    

        pt = 0
        DO i = nv*(hbp-1)+nv+1,nnd 
          pt = pt + 1

          mesh%elxyh(i,el,1) = mesh%xyh_interior(pt,el,1)
          mesh%elxyh(i,el,2) = mesh%xyh_interior(pt,el,2)
          mesh%elhb(i,el)    = mesh%xyh_interior(pt,el,3)           
        ENDDO
      ENDDO
      
      RETURN
      END SUBROUTINE gather_coordinates

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  

      END MODULE group_coordinates