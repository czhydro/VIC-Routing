c  SUBROUTINES FOR INITIALIZATION (roughly)
c  init_array()
c  create_vic_names()
c  search_catchment()
c

      SUBROUTINE INIT_ARRAY( A, NROW, NCOL, VALUE )

C     Initialiase float array A to VALUE

      IMPLICIT NONE

      INTEGER NCOL, NROW
      INTEGER I, J
      REAL A(NCOL,NROW)
      REAL VALUE

      DO J=1, NROW
         DO I=1, NCOL
            A(I,J)=VALUE
         END DO
      END DO

      RETURN
      END

      SUBROUTINE MATCH_GRID(INPATH,JLOC,ILOC, FluxFile)
c     calculate the shortest distance between 
c     main VIC and the grids points created 
c     by routing module and takes the VIC flux
c     with the shortest distance to read the data
c     Added by - Vikalp Mishra (NASA SERVIR)
c              - May 7 2021

      IMPLICIT NONE
         
      REAL LAT1, LON1, LAT2, LON2, rLAT2, rLON2
      REAL TMP1, TMP2, DEL_LAT, DEL_LON
      REAL R, JLOC, ILOC, MDIS, FLAT, FLON
      REAL, ALLOCATABLE::DIST(:)
      logical ex1
      INTEGER N,I,IOS,IND, stat, mloc
   
      CHARACTER*72 INPATH,LINE, cmd,FluxFile
      CHARACTER*100, ALLOCATABLE::FLIST(:)
      CHARACTER*72 S1,S2
	

      inquire(File = TRIM(ADJUSTL('list.txt')), exist = ex1)
      if (ex1) then 
         OPEN(11,FILE=TRIM(ADJUSTL('list.txt')),IOSTAT=IOS)
      else 
         WRITE(cmd,*)trim('ls '//INPATH)//trim('* > list.txt')
         STAT = SYSTEM(CMD)
         OPEN(11,FILE=TRIM(ADJUSTL('list.txt')),IOSTAT=IOS)
      endif

      !WRITE(cmd,*)trim('ls '//INPATH)//trim('* > list.txt')
      !STAT = SYSTEM(CMD)
      !OPEN(11,FILE=TRIM(ADJUSTL('list.txt')),IOSTAT=IOS)


      IF (IOS/=0) STOP 'ERROR OPENING FLUX LIST FILE'
      N = 0
      DO 
	READ(11,'(A90)', IOSTAT= IOS) LINE
	IF (IOS /= 0) EXIT
    	N = N+1
      END DO
      !WRITE(*,*)'LIST HAS TOTAL ',N, 'FLUX FILES'
      ALLOCATE(FLIST(N))
      allocate(dist(n))
      REWIND(11)
      DO I = 1,N
	READ(11,'(A90)') FLIST(I)
	IND = SCAN(FLIST(I),'_') ! SEPERATING STRIB BY UNDERSCRES
	S1 = FLIST(I)(IND+1:)
	IND = SCAN(S1,'_')
	READ(S1(1:IND-1),'(F7.4)')LAT2
	READ(S1(IND+1:),'(F7.4)')LON2
        !write(*,*) 'LAT/LON 2: ', LAT2, lon2
        R = 6371000 ! EARTH'S RADIUS IN METER
	DEL_LAT = (JLOC - LAT2)*3.14/180
	DEL_LON = (ILOC - LON2)*3.14/180
	! LAT/LON(s) IN RADIANS
     	LAT1 = JLOC/180
	LON1 = ILOC/180
	rLAT2 = LAT2/180
	rLON2 = LON2/180
      
	TMP1 = ((DEL_LAT/2)**2)+COS(LAT1)*COS(rLAT2)*(SIN(DEL_LON/2)**2)
	TMP2 = 2*ATAN2(SQRT(TMP1),SQRT(1-TMP1))
	DIST(I) = R*TMP2
      END DO
      mdis = minval(DIST)
      mloc = minval(minloc(DIST))
      !write(*,*)mdis,mloc, flist(mloc)
      write(FluxFile,*)trim(adjustl(flist(mloc)))
      close(11)

      RETURN	
      END

      SUBROUTINE CREATE_VIC_NAMES( JLOC, ILOC, EXTEN, CLEN, DPREC )


c     create string containing vic file names to be
c     appended to path given in input file

c     filenames allowed a maximum of 5 decimal places

      IMPLICIT NONE

      CHARACTER*10 JICHAR(2)
      CHARACTER*20 EXTEN
      REAL JLOC, ILOC
      INTEGER NSPACE, CLEN, CLEN_OLD, DPREC, I

      !WRITE(JICHAR(1),'(F10.5)')JLOC
      !WRITE(JICHAR(2),'(F10.5)')ILOC
      !write(*,*)JLOC
      !write(*,*)ILOC
      CLEN_OLD=1
      DO I=1,2
         NSPACE=1
 5       IF(JICHAR(I)(NSPACE:NSPACE).EQ.' ')THEN
            NSPACE=NSPACE+1
            GOTO 5
         ENDIF
         CLEN=CLEN_OLD+11-NSPACE-5+DPREC
         EXTEN(CLEN_OLD:CLEN)=JICHAR(I)(NSPACE:5+DPREC)
         IF(I.EQ.1)THEN
            EXTEN(CLEN:CLEN)='_'
         ENDIF
         CLEN_OLD=CLEN+1
      END DO

      CLEN=CLEN-1

      RETURN
      END





      SUBROUTINE SEARCH_CATCHMENT
     & (PI,PJ,DIREC,NCOL,NROW,NO_OF_BOX,CATCHIJ,PMAX,
     $  IROW,ICOL)

      IMPLICIT NONE

      INTEGER PI,PJ,I,J,NCOL,NROW,PMAX,ICOL,IROW
      INTEGER II, JJ, III, JJJ
      INTEGER DIREC(NCOL,NROW,2)
      INTEGER NO_OF_BOX
      INTEGER CATCHIJ(PMAX,2)

C****** CATCHMENTS ***************************************

      NO_OF_BOX = 0

      DO I = 1, ICOL
         DO J = 1, IROW
            !print*, pi, pj,i, j, no_of_box
            II = I
            JJ = J
 300        CONTINUE
            IF ((II .GT. ICOL) .OR. (II .LT.1) .OR.
     &          (JJ .GT. IROW) .OR. (JJ .LT.1)) THEN
               GOTO 310
            END IF
            IF ((II .EQ. PI) .AND. (JJ .EQ. PJ)) THEN
               NO_OF_BOX = NO_OF_BOX + 1
               CATCHIJ(NO_OF_BOX,1) = I
               CATCHIJ(NO_OF_BOX,2) = J
               GOTO 310
            ELSE
               IF ((DIREC(II,JJ,1).NE.0) .AND.    !check if the current
     &             (DIREC(II,JJ,2) .NE.0)) THEN   !ii,jj cell routes down
                     III = DIREC(II,JJ,1)         !to the subbasin outlet
                     JJJ = DIREC(II,JJ,2)         !point, following the
                     II  = III                    !direction of direc(,)
                     JJ  = JJJ                    !from each cell
                     !write(*,*) i,j, II, JJ
                     GOTO 300
               END IF   			  !if you get there,
            END IF                                !no_of_box increments
 310        CONTINUE                              !and you try another
         END DO                                   !cell.
      END DO

      WRITE(*,*) 'Number of grid cells upstream of present station',
     $     no_of_box

      RETURN
      END
