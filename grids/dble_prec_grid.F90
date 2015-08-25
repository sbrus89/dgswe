      PROGRAM dble_pres_grid

      IMPLICIT NONE
      
      INTEGER, PARAMETER :: prec = kind(1d0)
      CHARACTER(100) :: grid_file_input,grid_file_output
      CHARACTER(100) :: forcing_file_input,forcing_file_output      
      CHARACTER(100) :: grid_name
      CHARACTER(10) :: obtag,fbtag      
      INTEGER :: i,j,k,el
      INTEGER :: ne,nn,nnd
      INTEGER :: nope,neta,nbseg,obnd
      INTEGER :: nbou,nvel,btype,fbnd
      INTEGER :: nobfr    
      INTEGER :: nfbfr      
      REAL(prec) :: obfreq,obnfact,obeq
      REAL(prec) :: obamp,obph     
      REAL(prec) :: fbfreq,fbnfact,fbeq
      REAL(prec) :: fbamp,fbph       
      LOGICAL :: file_exists
      LOGICAL :: any_nfb
      
      REAL(prec) :: x,y,h
      INTEGER, DIMENSION(16) :: ect
      INTEGER, ALLOCATABLE, DIMENSION(:,:) :: fbseg
      




      grid_file_input = "inlet1.grd"
      forcing_file_input = "inlet1.bfr"
      
      grid_file_output = "inlet1_dble.grd"
      forcing_file_output = "inlet1_dble.bfr"

      
      
      
      ! open grid file
      INQUIRE(FILE=TRIM(ADJUSTL(grid_file_input)), EXIST = file_exists)
      IF(file_exists == .FALSE.) THEN
        PRINT*, "grid file does not exist"
        STOP  
      ENDIF
      
      OPEN(UNIT=14,  FILE=TRIM(ADJUSTL(grid_file_input)))     
      OPEN(UNIT=141, FILE=TRIM(ADJUSTL(grid_file_output)))   
                       
      ! name of grid
      READ(14,"(A)")   grid_name    
      WRITE(141,"(A)") grid_name                                     

      ! number of elements and number of nodes
      READ(14,*), ne, nn     
      WRITE(141,"(2(I8,1x))") ne,nn      
      
      ! node coordinates and depths
      DO i = 1,nn                                                      
        READ(14,*) j, x, y, h 
        WRITE(141,"(I8,1x,3(D24.17,1x))") j, x, y, h 
      ENDDO

      ! element connectivity
      DO i = 1,ne
        READ(14,*) el,nnd,(ect(j), j=1,nnd)       
        WRITE(141,"(5(I8,1x))") el,nnd,(ect(j), j=1,nnd) 
      ENDDO           
      
      ! open boundaries
      READ(14,*) nope             
      WRITE(141,"(I8,19x,A)") nope, "! number of open boundaries"
      READ(14,*) neta  
      WRITE(141,"(I8,19x,A)") neta, "! number of total elevation specified boundary nodes"

      DO i = 1,nope                                                     
        READ(14,*) nbseg  ! number of nodes in open boundary segment
        WRITE(141,"(I8,10x,A,1x,I8)") nbseg, "! number of nodes in open boundary", i
        DO j = 1,nbseg
          READ(14,*) obnd  ! open boundary node numbers
          WRITE(141,"(I8)") obnd
        ENDDO
      ENDDO

      ! normal flow boundaries
      READ(14,*) nbou  
      WRITE(141,"(I8,19x,A)") nbou, "! number of normal flow boundaries"
      READ(14,*) nvel  
      WRITE(141,"(I8,19x,A)") nvel, "! total number of normal flow nodes"
      
      ALLOCATE(fbseg(2,nbou))

      DO i = 1,nbou
        READ(14,*), nbseg, btype ! number of nodes in segment, boundary type
        WRITE(141,"(I8,1x,I8,10x,A,1x,I8)") nbseg, btype, "! number of nodes in normal flow boundary", i
        fbseg(1,i) = nbseg
        fbseg(2,i) = btype
        DO j = 1,nbseg
          READ(14,*), fbnd       ! normal flow boundary node numbers
           WRITE(141,"(I8)") fbnd
        ENDDO

      ENDDO

      CLOSE(14) 
      CLOSE(141)
      
      
      
      
      ! open forcing file
      INQUIRE(FILE=TRIM(ADJUSTL(forcing_file_input)), EXIST = file_exists)
      IF(file_exists == .FALSE.) THEN
        PRINT*, "forcing file does not exist"
        STOP  
      ENDIF
      
      OPEN(UNIT=15,  FILE=TRIM(ADJUSTL(forcing_file_input)))     
      OPEN(UNIT=151, FILE=TRIM(ADJUSTL(forcing_file_output)))         
      
      
      ! open boundary forcing data
      READ(15,*) nobfr
      WRITE(151,"(I8)") nobfr

      DO i = 1,nobfr
        READ(15,*) obtag
        WRITE(151,"(A)") obtag
        
        READ(15,*) obfreq,obnfact,obeq
        WRITE(151,"(3(D24.17,1x))") obfreq,obnfact,obeq
      ENDDO

      DO i = 1,nobfr
        READ(15,*) obtag
        WRITE(151,"(A)") obtag
        
        DO j = 1,neta
          READ(15,*) obamp,obph
          WRITE(151,"(2(D24.17,1x))") obamp,obph
        ENDDO
      ENDDO


      ! flow boundary forcing data
      any_nfb = .false. ! determine if there are normal flow boundaries
      DO i = 1,nbou
        btype = fbseg(2,i)
        IF(btype == 2 .OR. btype == 12 .OR. btype == 22)THEN
          any_nfb = .true.
        ENDIF        
      ENDDO      

      nfbfr = 0
      IF (any_nfb) THEN
        READ(15,*) nfbfr
        WRITE(151,"(I8)") nfbfr

        DO i = 1,nfbfr
          READ(15,*) fbtag    
          WRITE(151,"(A)") fbtag
          
          READ(15,*) fbfreq,fbnfact,fbeq
          WRITE(151,"(3(D24.17,1x))") fbfreq,fbnfact,fbeq
        ENDDO


        DO i = 1,nfbfr
          READ(15,*) fbtag 
          WRITE(151,"(A)") fbtag
          
          DO j = 1,nbou
            btype = fbseg(2,j)
            IF(btype == 2 .OR. btype == 12 .OR. btype == 22)THEN
              DO k = 1,fbseg(1,j)
                READ(15,*) fbamp,fbph
                WRITE(151,"(2(D24.17,1x))") fbamp,fbph
              ENDDO
            ENDIF
          ENDDO
        ENDDO
      ENDIF
      
      CLOSE(15)      
      CLOSE(151)

     
      END PROGRAM dble_pres_grid