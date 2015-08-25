      MODULE messenger2

#ifdef CMPI      
      USE mpi
#endif      

!$    USE omp_lib 

      USE globals, ONLY: pres,nrblk

      IMPLICIT NONE
      
      INTEGER :: ierr
      INTEGER :: nproc
      INTEGER :: nthreads      
      INTEGER :: myrank,myid
      INTEGER :: world_group
      INTEGER :: comp_group
      INTEGER :: comp_comm
      INTEGER :: comm_dist_graph
       
      CHARACTER(6) :: dirname
      INTEGER :: lname      
      
      INTEGER nproc_sr
      INTEGER, ALLOCATABLE, DIMENSION(:) :: proc_sr
      INTEGER, ALLOCATABLE, DIMENSION(:) :: ned_sr
      INTEGER, ALLOCATABLE, DIMENSION(:,:) :: el_sr
      INTEGER, ALLOCATABLE, DIMENSION(:,:) :: led_sr
      INTEGER, ALLOCATABLE, DIMENSION(:) :: proc_group
      
      INTEGER, ALLOCATABLE, DIMENSION(:) :: lel2gel
      
      INTEGER :: nred
      INTEGER, ALLOCATABLE, DIMENSION(:) :: redn
      
      INTEGER :: mnired
      INTEGER :: mnelred
      
      TYPE :: edge_ptr_array
        REAL(pres), POINTER :: ptr
      END TYPE edge_ptr_array
      
      TYPE(edge_ptr_array), ALLOCATABLE, DIMENSION(:,:) :: Hri !,Hre
      TYPE(edge_ptr_array), ALLOCATABLE, DIMENSION(:,:) :: Zri,Zre      
      TYPE(edge_ptr_array), ALLOCATABLE, DIMENSION(:,:) :: Qxri,Qxre
      TYPE(edge_ptr_array), ALLOCATABLE, DIMENSION(:,:) :: Qyri,Qyre   
      TYPE(edge_ptr_array), ALLOCATABLE, DIMENSION(:,:) :: xmri,ymri,xymri    
      TYPE(edge_ptr_array), ALLOCATABLE, DIMENSION(:,:) :: Zfri,Hfri,Qxfri,Qyfri
      
      REAL(pres), ALLOCATABLE, DIMENSION(:) :: xmre,ymre,xymre
      REAL(pres), ALLOCATABLE, DIMENSION(:) :: Hre
      
      REAL(pres), ALLOCATABLE, TARGET, DIMENSION(:,:) :: sol_recv
      REAL(pres), ALLOCATABLE, TARGET, DIMENSION(:,:) :: sol_send      

      TYPE(edge_ptr_array), ALLOCATABLE, DIMENSION(:,:) :: send_ptr    
      
      REAL(pres), ALLOCATABLE, DIMENSION(:,:) :: hbr      
      REAL(pres), ALLOCATABLE, DIMENSION(:,:) :: rnx,rny
      REAL(pres), ALLOCATABLE, DIMENSION(:,:) :: rcfac      
      REAL(pres), ALLOCATABLE, DIMENSION(:,:) :: detJe_recv      
      
      INTEGER, ALLOCATABLE, DIMENSION(:) :: solreq_recv
      INTEGER, ALLOCATABLE, DIMENSION(:) :: solreq_send
      INTEGER, ALLOCATABLE, DIMENSION(:) :: solreq      
      
      INTEGER, ALLOCATABLE, DIMENSION(:) :: win      
      
#ifdef ALIGN64
!DIR$ ATTRIBUTES ALIGN:64 :: Zri,Zre,Hri,Hre,Qxri,Qxre,Qyri,Qyre
!DIR$ ATTRIBUTES ALIGN:64 :: xmri,xmre,ymri,ymre,xymri,xymre
!DIR$ ATTRIBUTES ALIGN:64 :: Zfri,Hfri,Qxfri,Qyfri
!DIR$ ATTRIBUTES ALIGN:64 :: rnx,rny,len_area_recv
!DIR$ ATTRIBUTES ALIGN:64 :: sol_send,sol_recv
!DIR$ ATTRIBUTES ALIGN:64 :: send_ptr
#endif

#ifdef ALIGN32
!DIR$ ATTRIBUTES ALIGN:32 :: Zri,Zre,Hri,Hre,Qxri,Qxre,Qyri,Qyre
!DIR$ ATTRIBUTES ALIGN:32 :: xmri,xmre,ymri,ymre,xymri,xymre
!DIR$ ATTRIBUTES ALIGN:32 :: Zfri,Hfri,Qxfri,Qyfri
!DIR$ ATTRIBUTES ALIGN:32 :: rnx,rny,len_area_recv
!DIR$ ATTRIBUTES ALIGN:32 :: sol_send,sol_recv
!DIR$ ATTRIBUTES ALIGN:32 :: send_ptr
#endif 

#ifdef ALIGN16
!DIR$ ATTRIBUTES ALIGN:16 :: Zri,Zre,Hri,Hre,Qxri,Qxre,Qyri,Qyre
!DIR$ ATTRIBUTES ALIGN:16 :: xmri,xmre,ymri,ymre,xymri,xymre
!DIR$ ATTRIBUTES ALIGN:16 :: Zfri,Hfri,Qxfri,Qyfri
!DIR$ ATTRIBUTES ALIGN:16 :: rnx,rny,len_area_recv
!DIR$ ATTRIBUTES ALIGN:16 :: sol_send,sol_recv
!DIR$ ATTRIBUTES ALIGN:16 :: send_ptr
#endif

      CONTAINS
      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      SUBROUTINE directory_name()
       
      IMPLICIT NONE

#ifdef CMPI       
      dirname = "PE0000"
      lname = 6
       
      WRITE(dirname(3:lname),"(I4.4)") myrank
       
!      PRINT*, "dirname = ", dirname
#else
      dirname = "."
      lname = 1
#endif       
       
      RETURN
      END SUBROUTINE directory_name

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       
      SUBROUTINE message_init()
         
      IMPLICIT NONE
      
      CHARACTER(10) :: date,time
      
          
       
#ifdef CMPI       
      CALL MPI_INIT(ierr)
      CALL MPI_COMM_SIZE(MPI_COMM_WORLD,nproc,ierr)
      CALL MPI_COMM_RANK(MPI_COMM_WORLD,myrank,ierr)
!      CALL MPI_COMM_GROUP(MPI_COMM_WORLD,world_group,ierr)


      IF(myrank == 0) THEN       
        PRINT*, "Processor rank is ", myrank, "/", nproc
      ENDIF
       
      CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)       
      
      nrblk = 1
#else
      myrank = 0
      nrblk = 1      
#endif

#ifdef openmp      
      nthreads = omp_get_max_threads()        
      
      PRINT*, "Thread numbers :"
!$OMP  parallel private(myid)
      myid = omp_get_thread_num()
      PRINT*, myid        
!$OMP end parallel

      nrblk = nthreads  
#endif

      CALL directory_name()
      
      
      CALL DATE_AND_TIME(date,time)     
      IF (myrank == 0) THEN
        PRINT*, ' '
        PRINT*, 'Started: ',date,time
        PRINT*, ' '        
      ENDIF      
      
      CALL version()
       
      RETURN           
      END SUBROUTINE message_init
       
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!      


      SUBROUTINE read_message_files()
      
      USE globals, ONLY: pres,ne,mndof,nlines
      
      IMPLICIT NONE
      
#ifdef CMPI      
      INTEGER :: pe,ed,el,j
      LOGICAL :: file_exists
      INTEGER :: npe,tne,mne,ne2,ndof2
      REAL(pres) :: lines2
      INTEGER :: reorder,info
      
      INQUIRE(FILE=dirname(1:lname)//'/'//'fort.18', EXIST = file_exists)
      IF(file_exists == .FALSE.) THEN
        PRINT*, "fort.18 file does not exist for process", myrank
        CALL abort()
      ENDIF     
      
      OPEN(UNIT=18,FILE=dirname(1:lname)//'/'//'fort.18')
      
      READ(18,182) nproc_sr
      
      ALLOCATE(proc_sr(nproc_sr),ned_sr(nproc_sr))
      ALLOCATE(el_sr(ne,nproc_sr),led_sr(ne,nproc_sr))
      
      DO pe = 1,nproc_sr
!        j = mod(pe-1+myrank,nproc_sr)+1 ! DG ADCIRC method
!        READ(18,181) proc_sr(j),ned_sr(j) 
!        READ(18,180) (el_sr(ed,j), ed = 1,ned_sr(j))
!        READ(18,180) (led_sr(ed,j), ed = 1,ned_sr(j))
        

        READ(18,181) proc_sr(pe),ned_sr(pe) !DG ADCIRC may try to do some straggering here, not really sure about what they are doing 
        READ(18,180) (el_sr(ed,pe), ed = 1,ned_sr(pe))
        READ(18,180) (led_sr(ed,pe), ed = 1,ned_sr(pe))
      ENDDO
 
      CLOSE(18)

!      ALLOCATE(proc_group(nproc_sr+1))
!      DO pe = 1,nproc_sr
!        proc_group(pe) = proc_sr(pe)
!      ENDDO      
!      proc_group(nproc_sr+1) = myrank

!      CALL MPI_GROUP_INCL(world_group,nproc_sr+1,proc_group,comp_group,ierr)
!      CALL MPI_COMM_CREATE(MPI_COMM_WORLD,comp_group,comp_comm,ierr)

      reorder = 1
      info = MPI_INFO_NULL

      CALL MPI_Dist_graph_create_adjacent(MPI_COMM_WORLD,nproc_sr,proc_sr,ned_sr, &
                                                         nproc_sr,proc_sr,ned_sr, &
                                                         info,reorder,comm_dist_graph,ierr)      
      
      INQUIRE(FILE=dirname(1:lname)//'/'//'fort.80', EXIST = file_exists)
      IF(file_exists == .FALSE.) THEN
        PRINT*, "fort.80 file does not exist for process", myrank
        CALL abort()
      ENDIF   
      
      OPEN(UNIT=80,FILE=dirname(1:lname)//'/'//'fort.80')  
      
      READ(80,*) npe
      READ(80,*) tne
      READ(80,*) mne
      READ(80,*) ne2
      READ(80,*) ndof2
      READ(80,*) lines2
      
      ! Check the post processing information against run input 
      IF(npe /= nproc) THEN
        PRINT*, "number of processors does not agree with files"
        CALL abort()
      ENDIF
      
      IF(ne /= ne2) THEN
        PRINT*, "number of elements does not agree with files"
        CALL abort()
      ENDIF
      
      IF(ndof2 /= mndof) THEN
        PRINT*, "number of degrees of freedom does not agree with files"
        CALL abort()
      ENDIF
      
      IF(lines2 /= nlines) THEN
        PRINT*, "number of output lines does not agree with files"
        CALL abort()
      ENDIF
      
      ALLOCATE(lel2gel(ne))
      
      DO el = 1,ne
        READ(80,*) lel2gel(el)
      ENDDO
      
      CLOSE(80)
      
      
 180  FORMAT(8X,9I8)
 181  FORMAT(8X,2I8) 
 182  FORMAT(8X,I8) 
#endif   
      RETURN      
      END SUBROUTINE read_message_files
    

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


      SUBROUTINE message_setup()
      
      USE globals, ONLY: ned,mnqpte,nqpte,detJe,nx_pt,ny_pt,hbqpte, &
                         Zflux,Hflux,Qxflux,Qyflux, &
                         Zqpt,Hqpt,Qxqpt,Qyqpt, &
                         xmom,ymom,xymom, &
                         ged2el,ged2led,gel2ael, &
                         Spe,cfac
      
      IMPLICIT NONE
      
#ifdef CMPI
      INTEGER :: pe,ed,pt
      INTEGER :: el,led,gpt,ged,el_in
      INTEGER :: i,edcnt
      INTEGER :: gp_in,gp_ex
      INTEGER :: mnedsr
      INTEGER :: match
      
      mnedsr = MAXVAL(ned_sr)
      
      ALLOCATE(send_ptr(3*mnedsr*mnqpte,nproc_sr))
      ALLOCATE(sol_send(3*mnedsr*mnqpte,nproc_sr))
      
      DO pe = 1,nproc_sr
        i = 0
        DO ed = 1,ned_sr(pe)
          el = el_sr(ed,pe)
          el_in = gel2ael(el)
          led = led_sr(ed,pe)          
          DO pt = 1,nqpte(1)
            i = i+1
            
!             gpt = (led-1)*nqpte(1) + pt
            gpt = (led-1)*nqpte(1) + nqpte(1) - pt + 1            
            
!             send_ptr(i,pe)%ptr => Hqpt(el_in,gpt)
            send_ptr(i,pe)%ptr => Zqpt(el_in,gpt)
            send_ptr(ned_sr(pe)*nqpte(1) + i,pe)%ptr => Qxqpt(el_in,gpt)
            send_ptr(2*ned_sr(pe)*nqpte(1) + i,pe)%ptr => Qyqpt(el_in,gpt)
          ENDDO
        ENDDO
      ENDDO
      
      
      ALLOCATE(sol_recv(3*mnedsr*mnqpte,nproc_sr))
!       ALLOCATE(Hre(nred,mnqpte))
      ALLOCATE(Zre(nred,mnqpte),Qxre(nred,mnqpte),Qyre(nred,mnqpte))
      ALLOCATE(Zri(nred,mnqpte),Hri(nred,mnqpte),Qxri(nred,mnqpte),Qyri(nred,mnqpte))
      ALLOCATE(xmri(nred,mnqpte),ymri(nred,mnqpte),xymri(nred,mnqpte))
      ALLOCATE(xmre(nred),ymre(nred),xymre(nred))
      ALLOCATE(Hre(nred))      
      ALLOCATE(hbr(nred,mnqpte))
      ALLOCATE(Zfri(nred,mnqpte),Hfri(nred,mnqpte),Qxfri(nred,mnqpte),Qyfri(nred,mnqpte))      
      
      edcnt = 0     
      DO pe = 1,nproc_sr
        i = 0
        DO ed = 1,ned_sr(pe)
          edcnt = edcnt + 1
          el = el_sr(ed,pe)
          led = led_sr(ed,pe)
          DO pt = 1,nqpte(1)
            i = i+1
            
!             gp_ex = nqpte(1) - pt + 1
            gp_ex = pt
            
!             Hre(edcnt,gp_ex)%ptr => sol_recv(i,pe)
            Zre(edcnt,gp_ex)%ptr => sol_recv(i,pe)            
            Qxre(edcnt,gp_ex)%ptr => sol_recv(ned_sr(pe)*nqpte(1) + i,pe)
            Qyre(edcnt,gp_ex)%ptr => sol_recv(2*ned_sr(pe)*nqpte(1) + i,pe)
            
            gp_in = (led-1)*nqpte(1) + pt
            el_in = gel2ael(el)
            
            hbr(edcnt,pt) = hbqpte(el_in,gp_in) ! assumes continuous bathymetry
            
            Hri(edcnt,pt)%ptr  => Hqpt(el_in,gp_in)
            Zri(edcnt,pt)%ptr  => Zqpt(el_in,gp_in)
            Qxri(edcnt,pt)%ptr => Qxqpt(el_in,gp_in)
            Qyri(edcnt,pt)%ptr => Qyqpt(el_in,gp_in)
            
            xmri(edcnt,pt)%ptr  => xmom(el_in,gp_in)
            ymri(edcnt,pt)%ptr  => ymom(el_in,gp_in)
            xymri(edcnt,pt)%ptr => xymom(el_in,gp_in)
            
            Hfri(edcnt,pt)%ptr  => Hflux(el_in,gp_in)
            Zfri(edcnt,pt)%ptr  => Zflux(el_in,gp_in)            
            Qxfri(edcnt,pt)%ptr => Qxflux(el_in,gp_in)
            Qyfri(edcnt,pt)%ptr => Qyflux(el_in,gp_in)
          ENDDO
        ENDDO
      ENDDO
      
      ALLOCATE(rnx(nred,mnqpte),rny(nred,mnqpte),rcfac(nred,mnqpte),detJe_recv(nred,mnqpte))
      
      edcnt = 0
      DO pe = 1,nproc_sr
        DO ed = 1,ned_sr(pe)
          edcnt = edcnt+1
          
          el = el_sr(ed,pe)
          led = led_sr(ed,pe)
          
          match = 0
          DO i = 1,nred
            ged = redn(i)
            IF(el == ged2el(1,ged) .and. led == ged2led(1,ged)) THEN
            
              DO pt = 1,nqpte(1)
                rnx(edcnt,pt) = nx_pt(ged,pt)*Spe(ged,pt)
                rny(edcnt,pt) = ny_pt(ged,pt)
                
                rcfac(edcnt,pt) = cfac(ed,pt)
              
                detJe_recv(edcnt,pt) = detJe(ged,pt)
              ENDDO
              
              match = 1
              
            ENDIF
          ENDDO
          
          IF(match == 0) THEN
            PRINT*, "global edge not found for recieve edge"
            CALL abort()
          ENDIF
          
          
        ENDDO
      ENDDO
#endif       
      
      RETURN
      END SUBROUTINE message_setup
     

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


      SUBROUTINE communication_setup()
      
      USE globals, ONLY: nqpte,pres
      
      IMPLICIT NONE
#ifdef CMPI      
      
      INTEGER :: pe,tag
      
      ALLOCATE(solreq_recv(nproc_sr))
      ALLOCATE(solreq_send(nproc_sr))
      ALLOCATE(solreq(2*nproc_sr))
      
      tag = 1
      
      DO pe = 1,nproc_sr
        CALL MPI_RECV_INIT(sol_recv(1,pe),3*ned_sr(pe)*nqpte(1),MPI_DOUBLE_PRECISION,proc_sr(pe), &
                           tag,comm_dist_graph,solreq_recv(pe),ierr)

!        CALL MPI_RECV_INIT(sol_recv(1,pe),3*ned_sr(pe)*nqpte(1),MPI_DOUBLE_PRECISION,proc_sr(pe), &
!                           tag,MPI_COMM_WORLD,solreq_recv(pe),ierr)

!        CALL MPI_RECV_INIT(sol_recv(1,pe),3*ned_sr(pe)*nqpte(1),MPI_DOUBLE_PRECISION,pe, &
!                           tag,comp_comm,solreq_recv(pe),ierr)
      ENDDO
      
      DO pe = 1,nproc_sr
        CALL MPI_SEND_INIT(sol_send(1,pe),3*ned_sr(pe)*nqpte(1),MPI_DOUBLE_PRECISION,proc_sr(pe), &
                           tag,comm_dist_graph,solreq_send(pe),ierr)

!        CALL MPI_SEND_INIT(sol_send(1,pe),3*ned_sr(pe)*nqpte(1),MPI_DOUBLE_PRECISION,proc_sr(pe), &
!                           tag,MPI_COMM_WORLD,solreq_send(pe),ierr)                           

!        CALL MPI_SEND_INIT(sol_send(1,pe),3*ned_sr(pe)*nqpte(1),MPI_DOUBLE_PRECISION,pe, &
!                           tag,comp_comm,solreq_send(pe),ierr)                                 
      ENDDO

      DO pe = 1,nproc_sr
        solreq(pe) = solreq_recv(pe)
        solreq(pe+nproc_sr) = solreq_send(pe)
      ENDDO
      
!       ALLOCATE(win(nproc_sr))
!       
!       DO pe = 1,nproc_sr
!         CALL MPI_Win_create(sol_send(1,pe),3*ned_sr(pe)*nqpte(1)*pres,pres, &
!                             MPI_INFO_NULL,MPI_COMM_WORLD,win(pe))
!       ENDDO
      
#endif        
      RETURN
      END SUBROUTINE communication_setup
    

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#ifdef CMPI
      SUBROUTINE message_recieve()
      
      IMPLICIT NONE

      CALL MPI_STARTALL(nproc_sr,solreq_recv,ierr)
      
      RETURN
      END SUBROUTINE message_recieve
#endif      

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#ifdef CMPI
      SUBROUTINE message_send()
      
      USE globals, ONLY: nqpte
      
      IMPLICIT NONE
      
      INTEGER :: pe,i
      
      DO pe = 1,nproc_sr
        DO i = 1,3*ned_sr(pe)*nqpte(1)
        
          sol_send(i,pe) = send_ptr(i,pe)%ptr
          
        ENDDO
      ENDDO
      
     CALL MPI_STARTALL(nproc_sr,solreq_send,ierr)
!       CALL MPI_STARTALL(2*nproc_sr,solreq,ierr)
      
      RETURN
      END SUBROUTINE message_send
#endif      

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

       SUBROUTINE finish()

       IMPLICIT NONE
       
       
#ifdef CMPI       
       CALL MPI_FINALIZE(ierr)
       
       IF(myrank == 0) THEN
         PRINT*, "MPI terminated, status = ", ierr
       ENDIF
#endif
       STOP

       RETURN
       END SUBROUTINE finish

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

       SUBROUTINE abort()

       IMPLICIT NONE
       
       INTEGER :: errorcode
       
       errorcode = 0
        
#ifdef CMPI       
       CALL MPI_ABORT(MPI_COMM_WORLD,errorcode,ierr)
       
       PRINT*, "MPI aborted, status = ", ierr
#else      

       STOP  
       
#endif

       RETURN
       END SUBROUTINE abort

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      END MODULE messenger2
