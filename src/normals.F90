      SUBROUTINE normals()

      USE globals, ONLY: pres,ned,el_type,ged2el,ged2led,elxy, &
                         nverts,nnds,nqpta,nqpte, &
                         psia,dpsidr,dpsids, &
                         nx_pt,ny_pt,Spe,cfac

      USE transformation, ONLY: element_transformation,cpp_transformation


      IMPLICIT NONE
      
      INTEGER :: i,nd,el,pt,ed,led,edpt
      INTEGER :: et,nnd,nqa,nqe,nv    
      REAL(pres) :: xpt,ypt
      REAL(pres) :: drdx,drdy,dsdx,dsdy,detJ
      REAL(pres) :: Sp
      REAL(pres) :: nx,ny
    

      
      nx_pt = 0d0
      ny_pt = 0d0          
      
      DO ed = 1,ned
      
        el = ged2el(1,ed)
        led = ged2led(1,ed)
        
        et = el_type(el) 
        
        nnd = nnds(et)
        nqa = nqpta(et)
        nqe = nqpte(et)
        nv = nverts(et)        
        
        DO i = 1,nqe
          pt = nqa + (led-1)*nqe+i
          edpt = (led-1)*nqe+i
          

          CALL element_transformation(nnd,elxy(:,el,1),elxy(:,el,2),psia(:,pt,et),xpt,ypt, &
                                      dpsidr(:,pt,et),dpsids(:,pt,et),drdx,drdy,dsdx,dsdy,detJ)
            
          CALL cpp_transformation(ypt,Sp)
                     
          
          IF (nv == 3) THEN
              SELECT CASE(led)
                CASE(1)
                  nx = drdx + dsdx
                  ny = drdy + dsdy
                CASE(2)
                  nx = -drdx
                  ny = -drdy                
                CASE(3)
                  nx = -dsdx
                  ny = -dsdy               
              END SELECT  
              
          ELSE IF(nv == 4) THEN        
              SELECT CASE(led)
                CASE(1)
                  nx = drdx
                  ny = drdy           
                CASE(2)
                  nx = dsdx
                  ny = dsdy            
                CASE(3)
                  nx = -drdx
                  ny = -drdy             
                CASE(4)
                  nx = -dsdx
                  ny = -dsdy      
              END SELECT             
          ENDIF
          
          nx_pt(ed,i) = nx/sqrt(nx*nx+ny*ny)
          ny_pt(ed,i) = ny/sqrt(nx*nx+ny*ny)           
          
          Spe(ed,i) = Sp
          cfac(ed,i) = ny_pt(ed,i)**2+(nx_pt(ed,i)*Sp)**2
        
        ENDDO
        
      ENDDO
      
  
      RETURN
      END SUBROUTINE normals