      PROGRAM convert
      
      USE globals, ONLY: rp
      USE read_dginp
      USE grid_file_mod, ONLY:read_header,read_coords,read_connectivity
      USE edge_connectivity_mod, ONLY: elements_per_node
      USE basis, ONLY: element_nodes,element_basis
      USE version, ONLY: version_information,gitSHA
      USE write_adcirc, ONLY: write_6364_header,write_6364_snap 
      USE read_write_output, ONLY: read_solution_full,time_snaps
      
      IMPLICIT NONE
           
      INTEGER :: myrank = 0     
      CHARACTER(100) :: grid_name
      INTEGER :: ne
      INTEGER :: nn     
      INTEGER :: nsnap
      INTEGER, ALLOCATABLE, DIMENSION(:) :: el_type      
      INTEGER, ALLOCATABLE, DIMENSION(:,:) :: ect 
      REAL(rp), ALLOCATABLE, DIMENSION(:,:) :: xy
      REAL(rp), ALLOCATABLE, DIMENSION(:) :: depth  
      
      INTEGER, PARAMETER :: nel_type = 4
      INTEGER :: ndof(nel_type)
      INTEGER :: mndof
      INTEGER :: np(nel_type)
      INTEGER :: mnp
      INTEGER :: nnds(nel_type)
      INTEGER :: mnnds      
      INTEGER :: nverts(nel_type)
      
      INTEGER :: mnepn
      INTEGER, ALLOCATABLE, DIMENSION(:) :: nepn 
      INTEGER, ALLOCATABLE, DIMENSION(:,:) :: epn     

      
      REAL(rp) :: Zval,Uval,Vval
      REAL(rp), ALLOCATABLE, DIMENSION(:,:,:) :: Z,Qx,Qy,hb      
      REAL(rp), ALLOCATABLE, DIMENSION(:) :: t      
      REAL(rp), ALLOCATABLE, DIMENSION(:) :: Z_node,Qx_node,Qy_node,hb_node   
      REAL(rp) :: u_node,v_node
      REAL(rp), ALLOCATABLE, DIMENSION(:,:) :: Z_nodal,u_nodal

      REAL(rp), ALLOCATABLE, DIMENSION(:) :: r,s
      REAL(rp), ALLOCATABLE, DIMENSION(:,:) :: phi
            
      INTEGER :: i,j,dof,snap
      INTEGER :: et,el,nd,nv,ndf
      INTEGER :: vert,n,nd_check
      
      CHARACTER(32) :: rundes
      CHARACTER(24) :: runid
      CHARACTER(24) :: agrid
      CHARACTER(40) :: sha1
      INTEGER :: tskp
      INTEGER :: nout
      REAL(rp) :: time
      INTEGER :: it
      
      INTEGER :: unit63
      INTEGER :: unit64
      
      
      CALL version_information(unit=6)
      
      CALL read_input(0,".")    
      
      tf = tf*86400d0
      CALL time_snaps(sol_opt,sol_snap,tf,dt,tskp,nout)
      
      CALL sizes(ndof,mndof,nverts,np,mnp,nnds,mnnds)
      
      CALL read_header(myrank,grid_file,grid_name,ne,nn)        
      CALL read_coords(nn,xy,depth,h0)     
      CALL read_connectivity(ne,ect,el_type)     
      
      CALL read_solution_full(out_direc,"Z.sol","N",t,Z) 
      CALL read_solution_full(out_direc,"Qx.sol","N",t,Qx)  
      CALL read_solution_full(out_direc,"Qy.sol","N",t,Qy)        
      CALL read_solution_full(out_direc,"hb.sol","N",t,hb)        
      
      nsnap = SIZE(Z,3)
      
      CALL elements_per_node(ne,nn,nverts,el_type,ect,nepn,mnepn,epn)   
      
      
      ALLOCATE(r(mnnds),s(mnnds))
      ALLOCATE(phi(mndof,mnnds))
      DO et = 1,nel_type
        CALL element_nodes(et,0,np(et),n,r,s)
        CALL element_basis(et,p,ndf,n,r,s,phi)
      ENDDO
      
      ALLOCATE(Z_nodal(nn,1))
      ALLOCATE(u_nodal(nn,2))
      
      
      ALLOCATE(Z_node(mnepn))
      ALLOCATE(Qx_node(mnepn))
      ALLOCATE(Qy_node(mnepn))
      ALLOCATE(hb_node(mnepn))
      
      rundes = gitSHA
      runid = sha1("dgswe.inp","./")
      agrid = sha1(grid_file,out_direc)
      
      CALL write_6364_header("fort.63",rundes,runid,agrid,nout,nn,dt,tskp,1,unit63)
      CALL write_6364_header("fort.64",rundes,runid,agrid,nout,nn,dt,tskp,2,unit64)      
            
      
      DO snap = 1,nsnap
      
        time = t(snap)
      
        DO nd = 1,nn
          DO i = 1,nepn(nd)
          
            el = epn(i,nd)
            et = el_type(el)
            nv = nverts(et)
            ndf = ndof(et)          
            
      vnum: DO j = 1,nv
              nd_check = ect(j,el)
              IF (nd_check == nd) THEN
                vert = j
                EXIT vnum
              ENDIF
            ENDDO vnum
            
            Z_node(i) = 0d0
            Qx_node(i) = 0d0
            Qy_node(i) = 0d0
            hb_node(i) = 0d0
            DO dof = 1,ndf
              Z_node(i)  = Z_node(i)  + Z(dof,el,snap)*phi(dof,vert)
              Qx_node(i) = Qx_node(i) + Qx(dof,el,snap)*phi(dof,vert)
              Qy_node(i) = Qy_node(i) + Qy(dof,el,snap)*phi(dof,vert)            
              hb_node(i) = hb_node(i) + hb(dof,el,1)*phi(dof,vert)
            ENDDO
            
          ENDDO
          
          Zval = -999d0
          Uval = -999d0
          Vval = -999d0
          DO i = 1,nepn(nd)
          
            IF (Z_node(i) > Zval) THEN
              Zval = Z_node(i)
            ENDIF
            
            u_node = Qx_node(i)/(Z_node(i)+hb_node(i))
            IF (u_node > Uval) THEN
              Uval = u_node
            ENDIF
            
            v_node = Qy_node(i)/(Z_node(i)+hb_node(i))
            IF (v_node > Vval) THEN
              Vval = v_node
            ENDIF
                  
          ENDDO
          
  !         Zval = 0d0
  !         Uval = 0d0
  !         Vval = 0d0
  !         DO i = 1,nepn(nd)
  !         
  !           Zval = Zval + Z_node(i)          
  !           Uval = Uval + Qx_node(i)/(Z_node(i)+hb_node(i))         
  !           Vval = Vval + Qy_node(i)/(Z_node(i)+hb_node(i))
  !                 
  !         ENDDO   
  !         
  !         Zval = Zval/nepn(nd)
  !         Uval = Uval/nepn(nd)
  !         Vval = Vval/nepn(nd)
          
          Z_nodal(nd,1) = Zval
          u_nodal(nd,1) = Uval
          u_nodal(nd,2) = Vval
          
        ENDDO
        
        it = int(time/dt)
        
        CALL write_6364_snap(unit63,nn,1,time,it,Z_nodal)
        CALL write_6364_snap(unit64,nn,2,time,it,u_nodal)
      
      ENDDO 
      

      END PROGRAM convert