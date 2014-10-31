      SUBROUTINE metis()

      USE globals, ONLY: nn,ned,ged2nn,nproc,mnepn,part

      IMPLICIT NONE
      
      INTEGER :: ed,nd,j
      INTEGER :: n1,n2,nd2
      INTEGER :: tot
      INTEGER :: tmp
      
      INTEGER, ALLOCATABLE, DIMENSION(:) :: nadjnds
      INTEGER, ALLOCATABLE, DIMENSION(:,:) :: adjnds
      
      INTEGER :: numflag,wgtflag
      INTEGER :: edgecut
      INTEGER :: nparts
      INTEGER, ALLOCATABLE, DIMENSION(:) :: xadj
      INTEGER, ALLOCATABLE, DIMENSION(:) :: adjncy
      INTEGER, ALLOCATABLE, DIMENSION(:) :: vwgt
      INTEGER, ALLOCATABLE, DIMENSION(:) :: adjwgt
      INTEGER :: options(5)
      

      
      ! Find the nodes adjacent each node
      ALLOCATE(nadjnds(nn))
      ALLOCATE(adjnds(mnepn,nn))
      
      nadjnds = 0
      DO ed = 1,ned
        n1 = ged2nn(1,ed) ! find the node numbers on each edge
        n2 = ged2nn(2,ed)
        
        nadjnds(n1) = nadjnds(n1) + 1 ! count the nodes adjacent to node n1
        nadjnds(n2) = nadjnds(n2) + 1 ! count the nodes adjacent to node n2
        
        adjnds(nadjnds(n1),n1) = n2 ! node n2 is adjacent to node n1
        adjnds(nadjnds(n2),n2) = n1 ! node n1 is adjacent to node n2       
      ENDDO
      
      ! sort the node numbers adjacent to each node 
      ! (I'm not sure why this matters, but it seems to reduce the edge cut count)
      DO nd = 1,nn
        DO ed = 1,nadjnds(nd)
          DO j = ed,nadjnds(nd)
            IF(adjnds(j,nd) < adjnds(ed,nd)) THEN
              tmp = adjnds(j,nd)
              adjnds(j,nd) = adjnds(ed,nd)
              adjnds(ed,nd) = tmp
            ENDIF
          ENDDO
        ENDDO
      ENDDO
      
!       DO nd = 1,nn
!         PRINT("(I7,5x,10(I5))"), nadjnds(nd), (adjnds(ed,nd), ed = 1,nadjnds(nd))
!       ENDDO      
      
      ! Compute the node weights
      ! (ADCPREP uses the number of adjacent nodes)
      ALLOCATE(vwgt(nn))
      
      DO nd = 1,nn
        vwgt(nd) = nadjnds(nd)
!         PRINT*, nd,vwgt(nd)
      ENDDO
      
      ! Create the adjacency arrays used by metis
      ! Also compute the edge weight 
      ! (ADCPREP uses the sum of the adjacent nodes for each node on the edge)
      ALLOCATE(xadj(nn+1))
      ALLOCATE(adjncy(2*ned),adjwgt(2*ned))
      
      xadj(1) = 1
      tot = 0
      DO nd = 1,nn
        DO ed = 1,nadjnds(nd)
          tot = tot + 1
          nd2 = adjnds(ed,nd)
          adjncy(tot) = nd2 
          adjwgt(tot) = nadjnds(nd) + nadjnds(nd2)
!           PRINT*, tot, adjncy(tot), adjwgt(tot)
        ENDDO
        xadj(nd+1) = tot+1
      ENDDO      
      
      numflag = 1 ! use fortran numbering
      wgtflag = 3 ! use edge and node weights
      
      options(1) = 1 ! use non-default options
      options(2) = 3 ! use sorted heavy edge matching type
      options(3) = 1 ! use multilevel recursive bisection algorithm during initial partitioning 
      options(4) = 3 ! minimize connectivity among the subdomains 
      options(5) = 0 ! always set to zero
      
      nparts = nproc ! number of partitions
      
      ALLOCATE(part(nn))
      
      CALL METIS_PartGraphKway(nn,xadj,adjncy,vwgt,adjwgt,wgtflag, &
                               numflag,nparts,options,edgecut,part)
      
      
      PRINT*,"edgecut = ", edgecut
      
      RETURN 
      END SUBROUTINE metis