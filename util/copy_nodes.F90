      PROGRAM copy_nodes
      
      USE globals, ONLY: base,eval,rp
      USE kdtree2_module
      USE grid_file_mod          

      IMPLICIT NONE
      
      INTEGER :: nd
      INTEGER :: clnd
      INTEGER :: nverts(4)
      REAL(rp) :: xy(2)
      TYPE(kdtree2), POINTER :: tree_xy      
      TYPE(kdtree2_result), ALLOCATABLE, DIMENSION(:) :: closest         

      
      nverts(1) = 3
      nverts(2) = 4
      nverts(3) = 3
      nverts(4) = 4         
      
      base%grid_file = '/home/sbrus/data-drive/galveston_spline_flux_fix/grids/galveston_tri.grd'
      eval%grid_file = '/home/sbrus/data-drive/galveston_spline_flux_fix/grids/galveston_quad.grd'
      
      CALL read_grid(base)
      CALL read_grid(eval)
      
      tree_xy => kdtree2_create(base%xy(1:2,1:base%nn), rearrange=.true., sort=.true.)
      ALLOCATE(closest(base%nn))        


      
      DO nd = 1,eval%nn
      
        xy(1) = eval%xy(1,nd)
        xy(2) = eval%xy(2,nd)
        CALL kdtree2_n_nearest(tp=tree_xy,qv=xy,nn=1,results=closest)       
        clnd = closest(1)%idx 
            
        eval%xy(1,nd) = base%xy(1,clnd)
        eval%xy(2,nd) = base%xy(2,clnd)
        eval%depth(nd) = base%depth(nd)
      
      ENDDO
      
      
      
      CALL write_header('fort.14_out',eval%grid_name,eval%ne,eval%nn)      
      CALL write_coords(eval%nn,eval%xy,eval%depth)      
      CALL write_connectivity(eval%ne,eval%ect,eval%el_type,nverts)      
      CALL write_open_boundaries(eval%nope,eval%neta,eval%obseg,eval%obnds)      
      CALL write_flow_boundaries(eval%nbou,eval%nvel,eval%fbseg,eval%fbnds)
            
      
      END PROGRAM copy_nodes
      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!       
      
      SUBROUTINE read_grid(mesh)
      
      USE globals, ONLY: rp,grid
      USE grid_file_mod      
      
      IMPLICIT NONE
      
      TYPE(grid) :: mesh
      REAL(rp) :: h0 = 0d0
      
      CALL read_header(0,mesh%grid_file,mesh%grid_name,mesh%ne,mesh%nn)      
      CALL read_coords(mesh%nn,mesh%xy,mesh%depth,h0)            
      CALL read_connectivity(mesh%ne,mesh%ect,mesh%el_type)      
      CALL read_open_boundaries(mesh%nope,mesh%neta,mesh%obseg,mesh%obnds) 
      CALL read_flow_boundaries(mesh%nbou,mesh%nvel,mesh%fbseg,mesh%fbnds)       
      CALL print_grid_info(mesh%grid_file,mesh%grid_name,mesh%ne,mesh%nn)            
      
      
      END SUBROUTINE read_grid