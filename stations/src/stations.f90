      PROGRAM stations

      USE globals, ONLY:rp,ndof,mndof,np,mnnds,nnds,nverts,nel_type,ne,nn,nepn,mnepn,epn, &                                    
                        xy,ect,el_type,elxy,elhb,H,Qx,Qy,elsta,hbsta,nsta,xysta,phi_sta,t
      USE read_dginp, ONLY:read_input,p,ctp,hbp,out_direc,lines
      USE basis, ONLY: element_basis
      USE edge_connectivity_mod, ONLY: elements_per_node
      USE find_element, ONLY: find_element_init,in_element
      USE messenger2, ONLY: message_init
      USE shape_functions_mod, ONLY: shape_functions_area_eval
      USE bathymetry_interp_mod, ONLY: bathymetry_interp_eval
      

      IMPLICIT NONE

      INTEGER :: el,nd,pt,sta,i,line,dof
      INTEGER :: n1,n2,nvert,elin,leds(4),et,n,npts,nv
      REAL(rp) :: r(1),s(1),hb,rs(2)
      REAL(rp) :: Hsta,Qxsta,Qysta
      REAL(rp), ALLOCATABLE, DIMENSION(:,:) :: l,phi      


      CALL message_init()
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
                        
      CALL read_grid()                 
      CALL elements_per_node(ne,nn,nverts,el_type,ect,nepn,mnepn,epn)                   

      
      
      
      
      CALL read_stations() 

      CALL find_element_init(nel_type,nverts,np,nnds,nn,xy,nepn,epn)         

      ALLOCATE(phi_sta(mndof,nsta))      
      ALLOCATE(elsta(nsta),hbsta(nsta))
      ALLOCATE(l(mnnds,1))
      ALLOCATE(phi(mndof,1))         
      
      
      ! Find element each station is located in and determine local coordinates of station
      npts = 1
      DO sta = 1,nsta
      
        CALL in_element(xysta(:,sta),el_type,elxy,elin,leds,rs)       

        et = el_type(elin)                  
        nv = nverts(et)     
        
        r(1) = rs(1)
        s(1) = rs(2)
        
        CALL shape_functions_area_eval(nv,hbp,n,npts,r,s,l)
        CALL bathymetry_interp_eval(n,elhb(:,elin),l(:,1),hb)            
        CALL element_basis(et,p,ndof(et),npts,r,s,phi)           
        
        phi_sta(:,sta) = phi(:,1)           
        elsta(sta) = elin                                
        hbsta(sta) = hb                                         
      
      ENDDO
      
        
        
        
        

      
      ALLOCATE(H(ne,mndof),Qx(ne,mndof),Qy(ne,mndof))   

      
      OPEN(UNIT=61,FILE=trim(out_direc) //"station_H.d")
      OPEN(UNIT=621,FILE=trim(out_direc) //"station_Qx.d")     
      OPEN(UNIT=622,FILE=trim(out_direc) //"station_Qy.d")           
      
      OPEN(UNIT=63, FILE=trim(out_direc) //"solution_H.d")
      OPEN(UNIT=641,FILE=trim(out_direc) //"solution_Qx.d")     
      OPEN(UNIT=642,FILE=trim(out_direc) //"solution_Qy.d")      
      
      READ(63,*) 
      READ(641,*)
      READ(642,*)
      
      WRITE(61,*) nsta
      WRITE(621,*) nsta
      WRITE(622,*) nsta
      
      DO line = 1,lines+1
        READ(63,*) t
        READ(641,*) t
        READ(642,*) t
        
        PRINT("(A,I5)"), "Writing time snap: ", line
        
        DO dof = 1,mndof
          READ(63,*) (H(el,dof), el = 1,ne)
          READ(641,*) (Qx(el,dof), el = 1,ne)
          READ(642,*) (Qy(el,dof), el = 1,ne)
        ENDDO
        
        WRITE(61,*) t
        WRITE(621,*) t
        WRITE(622,*) t           
                
        DO sta = 1,nsta
                 
          elin = elsta(sta)                 
          et = el_type(elin)                 
                    
          Hsta = 0d0
          Qxsta = 0d0
          Qysta = 0d0
          DO dof = 1,ndof(et)
            Hsta = Hsta + H(elin,dof)*phi_sta(dof,sta)
            Qxsta = Qxsta + Qx(elin,dof)*phi_sta(dof,sta)            
            Qysta = Qysta + Qy(elin,dof)*phi_sta(dof,sta)
          ENDDO                
          
          WRITE(61,*) xysta(1,sta),xysta(2,sta),Hsta,hbsta(sta)
          WRITE(621,*) xysta(1,sta),xysta(2,sta),Qxsta
          WRITE(622,*) xysta(1,sta),xysta(2,sta),Qysta
        ENDDO
        
      ENDDO
      
      CLOSE(63)
      CLOSE(641)
      CLOSE(642)
      
      CLOSE(61)
      CLOSE(621)
      CLOSE(622)      
      END PROGRAM stations