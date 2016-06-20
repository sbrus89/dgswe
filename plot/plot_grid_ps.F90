      PROGRAM plot_grid_ps
      
      USE globals, ONLY: rp,ne,nn,nverts,ndof,mndof,np,mnp,nnds,mnnds,nel_type,ect,xy,depth, &
                         nope,neta,obseg,obnds,nvel,nbou,fbseg,fbnds,bndxy,grid_name, &
                         el_type,elxy,elhb, &
                         nepn,epn,mnepn,ned,ged2el,ged2el,ged2led,ged2nn,ed_type,recv_edge, &
                         nied,iedn,nobed,obedn,nfbed,fbedn,nnfbed,nfbedn,nfbednn, &
                         psic,psiv
      USE grid_file_mod
      USE basis, ONLY: element_nodes,element_basis
      USE read_write_output, ONLY: read_solution_full
      USE read_dginp, ONLY: read_input,out_direc,p,ctp,hbp,grid_file, &
                            bathy_file,curve_file,cb_file_exists,hb_file_exists
      USE plot_mod, ONLY: write_psheader,plot_contours,plot_mesh,close_ps
      USE triangulation, ONLY: reference_element_delaunay
      USE edge_connectivity_mod
      USE curvilinear_nodes_mod
      USE transformation
      USE shape_functions_mod
      
      IMPLICIT NONE
      
      INTEGER :: coord_sys
      REAL(rp) :: slam0,sphi0      
      REAL(rp) :: xmin,xmax
      REAL(rp) :: ymin,ymax
      REAL(rp) :: rmin,rmax
      REAL(rp) :: smin,smax
      
      INTEGER :: el,nd,pt
      INTEGER :: et,nv
      REAL(rp) :: ax,bx
      REAL(rp) :: ay,by
      
      INTEGER :: ps
      INTEGER :: pc
      INTEGER :: pplt(4)
      INTEGER :: nplt(4)
      INTEGER :: mnpp
      INTEGER :: ntri(4)
      INTEGER :: nred 
      INTEGER, ALLOCATABLE, DIMENSION(:,:,:) :: rect

      INTEGER :: i,j
      INTEGER :: nsnap
      INTEGER :: snap
      INTEGER :: space
      INTEGER :: n,ndf,npts,nnd,dof
      REAL(rp) :: xpt,ypt
      REAL(rp), ALLOCATABLE, DIMENSION(:) :: t      
      REAL(rp), ALLOCATABLE, DIMENSION(:,:,:) :: Z,Qx,Qy,hb 
      REAL(rp), ALLOCATABLE, DIMENSION(:) :: r,s
      REAL(rp), ALLOCATABLE, DIMENSION(:,:,:) :: phi   
      REAL(rp), DIMENSION(:,:,:), ALLOCATABLE :: xyplt   
      REAL(rp), DIMENSION(:,:), ALLOCATABLE :: Z_val 
      REAL(rp), DIMENSION(:,:), ALLOCATABLE :: hb_val           
      REAL(rp), DIMENSION(:,:), ALLOCATABLE :: vel_val   
      REAL(rp) :: Qx_val,Qy_val,H_val
      INTEGER :: outside
      INTEGER, DIMENSION(:), ALLOCATABLE :: el_in
      REAL(rp) :: xbox_min,xbox_max,ybox_min,ybox_max

      ps = 6
      pc = 6
      
      CALL read_input(0,"../work")
      
      ndof(1) = (p+1)*(p+2)/2
      ndof(2) = (p+1)**2
      ndof(3) = ndof(1)
      ndof(4) = ndof(2)      
      mndof = maxval(ndof)
      
      nverts(1) = 3
      nverts(2) = 4
      nverts(3) = 3
      nverts(4) = 4
      
      np(1) = 1
      np(2) = 1
      np(3) = ctp
      np(4) = ctp  
      mnp = maxval(np)+1

      nnds(1) = 3
      nnds(2) = 4
      nnds(3) = (ctp+1)*(ctp+2)/2
      nnds(4) = (ctp+1)*(ctp+1)      
      mnnds = maxval(nnds)   
      
      pplt(1) = ps
      pplt(2) = ps
      pplt(3) = pc
      pplt(4) = pc

      nplt(1) = (ps+1)*(ps+2)/2
      nplt(2) = (ps+1)*(ps+1)
      nplt(3) = (pc+1)*(pc+2)/2
      nplt(4) = (pc+1)*(pc+1)
      mnpp = maxval(nplt)

      PRINT*, grid_file
      CALL read_header(0,grid_file,grid_name,ne,nn)        
      CALL read_coords(nn,xy,depth)
      CALL read_connectivity(ne,ect,el_type) 
      CALL init_element_coordinates(ne,ctp,el_type,nverts,xy,ect,elxy)                  
      CALL read_open_boundaries(nope,neta,obseg,obnds)            
      CALL read_flow_boundaries(nbou,nvel,fbseg,fbnds)      
      CALL read_bathy_file(0,bathy_file,hbp,ne,el_type,nverts,depth,ect,elhb,hb_file_exists)                  
      CALL read_curve_file(0,curve_file,ctp,nbou,xy,bndxy,cb_file_exists)      
      CALL print_grid_info(grid_file,grid_name,ne,nn)    
      
      CALL elements_per_node(ne,nn,nverts,el_type,ect,nepn,mnepn,epn)       
      CALL find_edge_pairs(ne,nverts,el_type,ect,nepn,epn,ned,ged2el,ged2nn,ged2led)      
      CALL find_interior_edges(ned,ged2el,nied,iedn,ed_type,recv_edge)      
      CALL find_open_edges(nope,obseg,obnds,ged2nn,nobed,obedn,ed_type,recv_edge)            
      CALL find_flow_edges(nbou,fbseg,fbnds,ged2nn,nnfbed,nfbedn,nfbednn,nfbed,fbedn,recv_edge,ed_type)     
      nred = 0
      CALL print_connect_info(mnepn,ned,nied,nobed,nfbed,nnfbed,nred)

      PRINT*, "Calculating curved boundary information..."
      CALL shape_functions_linear_at_ctp(nel_type,np,psiv)                   
      CALL eval_coordinates_curved(ctp,nnds,nverts,el_type,xy,ect,fbseg,fbnds, &
                                   nnfbed,nfbedn,nfbednn,ged2el,ged2led, &
                                   psiv,bndxy,elxy)     
      
      PRINT*, "Calculating additional ploting point coordinates..."
      space = 1  
      ALLOCATE(r(mnpp),s(mnpp))      
      ALLOCATE(psic(mnnds,mnpp,nel_type))
      DO et = 1,nel_type     
        CALL element_nodes(et,space,pplt(et),npts,r,s)                  
        CALL shape_functions_area_eval(et,np(et),nnd,npts,r,s,psic(:,:,et))     
      ENDDO                                    
           
      ALLOCATE(xyplt(mnpp,ne,2))
      DO el = 1,ne      
        et = el_type(el)                          
        nnd = nnds(et)
        npts = nplt(et)
        DO pt = 1,npts              
          CALL element_transformation(nnd,elxy(:,el,1),elxy(:,el,2),psic(:,pt,et),xpt,ypt)           
          xyplt(pt,el,1) = xpt
          xyplt(pt,el,2) = ypt
        ENDDO
      ENDDO       
             
      
      PRINT*, "Evaluating reference element coordinate information..."
      ALLOCATE(phi(mndof,mnpp,nel_type))
      ALLOCATE(rect(3,3*mnpp,nel_type))
      DO et = 1,nel_type
        CALL element_nodes(et,1,pplt(et),n,r,s)
        CALL element_basis(et,p,ndf,n,r,s,phi(:,:,et))
        CALL reference_element_delaunay(n,r,s,ntri(et),rect(:,:,et))
        
!         DO i = 1,3
!           PRINT "(*(I5))", (rect(i,j,et), j = 1,ntri(et))
!         ENDDO    
!         PRINT*, ""
      ENDDO      


      
      PRINT*, "Boxing and scaling coordinates..."

      xmax = -1d10
      ymax = -1d10
      xmin = 1d10
      ymin = 1d10
      
      xbox_min = 3.2e5
      xbox_max = 3.4e5
      ybox_min = 3.24e6
      ybox_max = 3.26e6
      
!       xbox_min = -xmin
!       xbox_max = -xmax
!       ybox_min = -ymin
!       ybox_max = -ymax
      
      ALLOCATE(el_in(ne))
      el_in = 1
      
      DO el = 1,ne      
        et = el_type(el)                          
        nnd = nnds(et)
        npts = nplt(et)
        outside = 0
        DO pt = 1,npts  
        
          xpt = xyplt(pt,el,1)
          ypt = xyplt(pt,el,2)
        
          IF (xpt > xmax) THEN
            xmax = xpt
          ENDIF
          
          IF (xpt < xmin) THEN
            xmin = xpt
          ENDIF

          IF (ypt > ymax) THEN
            ymax = ypt
          ENDIF
          
          IF (ypt < ymin) THEN
            ymin = ypt
          ENDIF          
          
          IF (xpt < xbox_min .or. xpt > xbox_max) THEN
            outside = 1
            xmin = xbox_min
            xmax = xbox_max
          ENDIF
          
          IF (ypt < ybox_min .or. ypt > ybox_max) THEN
            outside = 1
            ymin = ybox_min
            ymax = ybox_max
          ENDIF

        ENDDO
        IF (outside == 1) THEN
          el_in(el) = 0
        ENDIF
      ENDDO    
      
      
      rmin = 10d0
      rmax = 602d0
      smin = 0d0
      smax = 792d0
      
      ax = (rmin/(xmin-xmax)+rmax/(xmax-xmin))
      bx = -(rmin*xmax/(xmin-xmax)+rmax*xmin/(xmax-xmin))
      
      ay = ax     ! axis equal
      by = smin-ax*ymin            
      
      DO el = 1,ne
        et = el_type(el)
        npts = nplt(et)
        nv = nverts(et)
        DO nd = 1,npts
          xyplt(nd,el,1) = ax*xyplt(nd,el,1) + bx
          xyplt(nd,el,2) = ay*xyplt(nd,el,2) + by                       
        ENDDO
      ENDDO
      
      DO nd = 1,nn
        xy(1,nd) = ax*xy(1,nd) + bx
        xy(2,nd) = ay*xy(2,nd) + by
      ENDDO
      
      
      

        
      PRINT*, "Reading solutions..."  
      CALL read_solution_full("","../work/Z.sol","N",t,Z) 
      CALL read_solution_full("","../work/Qx.sol","N",t,Qx)  
      CALL read_solution_full("","../work/Qy.sol","N",t,Qy)        
      CALL read_solution_full("","../work/hb.sol","N",t,hb)           
      nsnap = SIZE(Z,3)   
      
      PRINT*, "Evaluating solutions at additional plotting points..."
      snap = 49
      ALLOCATE(Z_val(mnpp,ne)) 
      ALLOCATE(hb_val(mnpp,ne))
      ALLOCATE(vel_val(mnpp,ne))
      DO el = 1,ne
        et = el_type(el)
        npts = nplt(et)
        ndf = ndof(et)
        DO nd = 1,npts
          Z_val(nd,el) = 0d0
          hb_val(nd,el) = 0d0
          Qx_val = 0d0
          Qy_val = 0d0
          DO dof = 1,ndf
            Z_val(nd,el) = Z_val(nd,el) + Z(dof,el,snap)*phi(dof,nd,et)
            hb_val(nd,el) = hb_val(nd,el) + hb(dof,el,1)*phi(dof,nd,et) 
            Qx_val = Qx_val + Qx(dof,el,snap)*phi(dof,nd,et)           
            Qy_val = Qy_val + Qy(dof,el,snap)*phi(dof,nd,et)            
            H_val = Z_val(nd,el) + hb_val(nd,el)
            vel_val(nd,el) = sqrt((Qx_val/H_val)**2 + (Qy_val/H_val)**2)
          ENDDO    
        ENDDO
        
      ENDDO      
      
      PRINT*, "Writing PostScript file..."
      CALL write_psheader()      
      
      CALL plot_contours(nplt,ntri,rect,ne,el_type,el_in,xyplt,vel_val)
      
!       CALL plot_mesh(ne,nverts,el_type,el_in,xy,ect)
      
      CALL close_ps()

      END PROGRAM plot_grid_ps
      
