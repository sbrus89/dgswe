      SUBROUTINE create_101btype_stations()

      USE globals, ONLY: rp,nbou,fbseg,fbnds,xy,nsta,xysta
      USE read_dginp, ONLY: out_direc

      IMPLICIT NONE
      
      INTEGER, PARAMETER :: nista = 5 ! number of stations bewteen 101 boundary nodes
      INTEGER, PARAMETER :: nacross = 40 ! number of stations across sections
      INTEGER, PARAMETER :: nsections = 3
      REAL(rp) :: r(nacross),h
      INTEGER :: sta,i,sec
      INTEGER :: bou,nd,n1,n2
      INTEGER :: nnds,btype
      LOGICAL :: file_exists
      INTEGER :: nd_across(2,nsections)
      

      !grid_file = /home/sbrus/data-drive/galveston/dgswe/quad2/galveston_quad2.grd          
      
        nd_across(1,1) = 938
        nd_across(2,1) = 857
        
        nd_across(1,2) = 1154
        nd_across(2,2) = 1081
        
        nd_across(1,3) = 233
        nd_across(2,3) = 186     
           
      
        ! Use the nodes along 101 type boundaries
        
        h = 2d0/real(nista+1,rp)
        DO i = 1,nista
          r(i) = -1d0 + h*real(i,rp)
          PRINT*, r(i)
        ENDDO        
      
        nsta = 0
      
        DO bou = 1,nbou
          nnds = fbseg(1,bou)
          btype = fbseg(2,bou)
          IF (btype == 101) THEN
            DO nd = 1,nnds
              nsta = nsta + 1            
            ENDDO
          ENDIF
        ENDDO       
      
        ALLOCATE(xysta(2,(nista+1)*nsta+nsections*nacross))
      
        nsta = 0
           
      
        DO bou = 1,nbou
          nnds = fbseg(1,bou)
          btype = fbseg(2,bou)
          PRINT*, bou,nnds,btype
          IF (btype == 101) THEN
            DO nd = 1,nnds-1
              nsta = nsta + 1            
              xysta(1,nsta) = xy(1,fbnds(nd,bou))
              xysta(2,nsta) = xy(2,fbnds(nd,bou))             
              DO i = 1,nista
                nsta = nsta + 1            
                xysta(1,nsta) = .5d0*((1d0-r(i))*xy(1,fbnds(nd,bou)) + (1d0+r(i))*xy(1,fbnds(nd+1,bou)))
                xysta(2,nsta) = .5d0*((1d0-r(i))*xy(2,fbnds(nd,bou)) + (1d0+r(i))*xy(2,fbnds(nd+1,bou)))  
              ENDDO
            ENDDO
            nsta = nsta + 1            
            xysta(1,nsta) = xy(1,fbnds(nnds,bou))
            xysta(2,nsta) = xy(2,fbnds(nnds,bou))     
            EXIT ! only do first one
          ENDIF
        ENDDO  
        
        PRINT*, "101 nsta = ",nsta
        
        ! Calculate stations along line between two points
        
        h = 2d0/real(nacross+1,rp)
        DO i = 1,nacross
          r(i) = -1d0 + h*real(i,rp)
          PRINT*, r(i)
        ENDDO             
        
    
        
        DO sec = 1,nsections
          n1 = nd_across(1,sec)
          n2 = nd_across(2,sec)
          DO i = 1,nacross
            nsta = nsta + 1            
            xysta(1,nsta) = .5d0*((1d0-r(i))*xy(1,n1) + (1d0+r(i))*xy(1,n2))
            xysta(2,nsta) = .5d0*((1d0-r(i))*xy(2,n1) + (1d0+r(i))*xy(2,n2))
          ENDDO
        ENDDO
        
      
        OPEN(UNIT=15,FILE="./stations_101.d")
        WRITE(15,*) nsta
        DO sta = 1,nsta
          WRITE(15,"(2(e24.17,2x))") xysta(1,sta), xysta(2,sta)
        ENDDO
        CLOSE(15)
      
        
        
        
      RETURN
      END SUBROUTINE create_101btype_stations