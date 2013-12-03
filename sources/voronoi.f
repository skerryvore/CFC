      SUBROUTINE VORONOI(NCC,LCC,N,X,Y,LIST,LPTR,LEND,NT,LCCC,ICCC,
     .                   LCT,LTRI,IER)
      INTEGER NCC, LCC(*), N, LIST(*), LPTR(*), LEND(N),
     .        NT, LCT(NCC), IER, LTRI(9,*), ICCC(6,*) 
      DOUBLE PRECISION X(*), Y(*), LCCC(14,*)
C
C     DETERMINE A VORONOI MOSAIC. TO BE USED WITH THE FUNCTIONS FROM TRIPACK.
C
C     THIS IS DONE BY DETERMINING THE CIRCUMCIRCLE CENTERS OF A 
C     DELAUNAY TRIANGULATION, USING ITS TRIANGLE NEIGBOURSHIP RELATIONS AS 
C     ADJACENCY RELATION FOR THE CIRCUMCIRCLE CENTERS.
C
C     THE TRIANGULATION PRODUCED WITH TRLIST FROM TRIPACK CAN 
C     INCLUDE EMPTY TRIANGLES WITH AREA 0 (IF COLLINEAR POINTS EXIST).
C     THESE TRIANGLES MAKE TROUBLES IF THEY OCCUR ON THE BORDER OF THE 
C     TRIANGULATION. SO WE MUST REMOVE THEM.
C     
C     A.GEBHARDT <albrecht.gebhardt@uni-klu.ac.at>
C
C     ON ENTRY:
C     ... SEE E.G. TRLIST
C
C     ON EXIT:
C     NT   - NUMBER OF TRIANGLES
C     LCCC - LIST OF CIRCUMCIRCLE CENTERS (3 by NT):
C            ROW 1,2: X, Y COORDINATES
C            ROW 3  : TRIANGLE AREA (-1 = REMOVED TRIANGLE)
C            ROW 4  : ASPECT RATIO
C            ROW 5  : CIRCUMCIRCLE RADUIS
C            ROW 6 to 14: for additive weighted voronoi diagrams
C     ICCC   LIST OF NEIGHBOUR RELATIONS (6 by NT):
C            ROW 1-3: NEIGHBOUR TRIANGLE INDICES
C            ROW 4-6: NODE INDICES OF TRIANGLE
C
C     Modules required by VORONOI:  TRLIST, CIRCUM, QSORT, RMSHNB

C     DETERMINE ALL CIRCUMCIRCLE CENTERS
      INTEGER NROW, IT, IP(3), NTRI, IT1, IT2
      DOUBLE PRECISION X1, X2, X3, Y1, Y2, Y3, XC, YC, SA, AR, 
     .     XP(3), YP(3), DEPS, CR
      LOGICAL RATIO
      
      RATIO = .TRUE.
      DEPS = 1.0D-7
      NROW = 9

      CALL TRLIST(NCC,LCC,N,LIST,LPTR,LEND,NROW,NT,LTRI,LCT,IER)
      NTRI=NT

      DO 10 IT = 1,NTRI
         X1 = X(LTRI(1,IT))
         Y1 = Y(LTRI(1,IT))
         X2 = X(LTRI(2,IT))
         Y2 = Y(LTRI(2,IT))
         X3 = X(LTRI(3,IT))
         Y3 = Y(LTRI(3,IT))
         CALL CIRCUM(X1,Y1,X2,Y2,X3,Y3,RATIO,XC,YC,CR,SA,AR)
         LCCC(1,IT) = XC
         LCCC(2,IT) = YC
         LCCC(3,IT) = SA
         LCCC(4,IT) = AR
         LCCC(5,IT) = CR
         ICCC(1,IT) = LTRI(4,IT)
         ICCC(2,IT) = LTRI(5,IT)
         ICCC(3,IT) = LTRI(6,IT)
         ICCC(4,IT) = LTRI(1,IT)
         ICCC(5,IT) = LTRI(2,IT)
         ICCC(6,IT) = LTRI(3,IT)
 10   CONTINUE

C     REMOVE TRIANGLES WITH AREA 0 FROM THE BORDER, UPDATE NEIGBOUR RELATIONS.
C     AREA=0 MEANS COLLINEAR POINTS, 3 CASES (I'M NOT SURE IF THEY ALL WILL
C     OCCUR!):
C       ALL POINTS DIFFERENT
C       2   POINTS COINCIDE
C       ALL POINTS COINCIDE
C
 11   CONTINUE
      DO 20 IT = 1,NTRI
         IF (((LCCC(3,IT).EQ.0) .OR.(LCCC(4,IT).LE.DEPS)) .AND. 
     .        (ICCC(1,IT).EQ.0 .OR. 
     .        ICCC(2,IT).EQ.0 .OR. 
     .        ICCC(3,IT).EQ.0)) THEN
C     WE HAVE A 0-TRIANGLE AT THE BORDER
            XP(1)=X(ICCC(4,IT))
            XP(2)=X(ICCC(5,IT))
            XP(3)=X(ICCC(6,IT))
            YP(1)=Y(ICCC(4,IT))
            YP(2)=Y(ICCC(5,IT))
            YP(3)=Y(ICCC(6,IT))
            CALL QSORT(3,XP,IP)
C     CHECK IF ALL X AND/OR Y COORDINATES COINCIDE
            IF (XP(IP(1)) .EQ. XP(IP(3))) THEN
               CALL QSORT(3,XP,IP)
               IF (YP(IP(1)) .EQ. YP(IP(3))) THEN
C     ALL POINTS COINCIDE, ALL NEIGHBOUR TRIANGLES HAVE AREA 0
C     AND WILL BE REMOVED LATER.
C     NOTHING TO DO (?)
                  CONTINUE
               ELSE
C     X1=X2=X3, USE Y-COORDINATES 
                  IF ((YP(IP(1)).EQ.YP(IP(2))) 
     .                 .OR. (YP(IP(2)).EQ.YP(IP(3)))) THEN
C     TWO POINTS COINCIDE
                     IF (YP(IP(1)).EQ.YP(IP(2))) THEN
                        IT1=ICCC(IP(1),IT)
                        IT2=ICCC(IP(2),IT)
                        CALL RMSHNB(IT,IT1,IT2,LCCC,ICCC)
                     ELSE
                        IT1=ICCC(IP(2),IT)
                        IT2=ICCC(IP(3),IT)
                        CALL RMSHNB(IT,IT1,IT2,LCCC,ICCC)
                     ENDIF
                  ELSE
C     ALL POINTS DIFFERENT
                     IT1=ICCC(IP(1),IT)
                     IT2=ICCC(IP(2),IT)
                     CALL RMSHNB(IT,IT1,IT2,LCCC,ICCC)
                     IT1=ICCC(IP(3),IT)
                     IT2=ICCC(IP(2),IT)
                     CALL RMSHNB(IT,IT1,IT2,LCCC,ICCC)
                  ENDIF
               ENDIF
            ELSE
C     USE X COORDINATES
               IF ((XP(IP(1)).EQ.XP(IP(2))) 
     .              .OR. (XP(IP(2)).EQ.XP(IP(3)))) THEN
C     TWO POINTS COINCIDE
                  IF (XP(IP(1)).EQ.XP(IP(2))) THEN
                     IT1=ICCC(IP(1),IT)
                     IT2=ICCC(IP(2),IT)
                     CALL RMSHNB(IT,IT1,IT2,LCCC,ICCC)
                  ELSE
                     IT1=ICCC(IP(2),IT)
                     IT2=ICCC(IP(3),IT)
                     CALL RMSHNB(IT,IT1,IT2,LCCC,ICCC)
                  ENDIF
               ELSE
C     ALL POINTS DIFFERENT
                  IT1=ICCC(IP(1),IT)
                  IT2=ICCC(IP(2),IT)
                  CALL RMSHNB(IT,IT1,IT2,LCCC,ICCC)
                  IT1=ICCC(IP(2),IT)
                  IT2=ICCC(IP(3),IT)
                  CALL RMSHNB(IT,IT1,IT2,LCCC,ICCC)
               ENDIF
            ENDIF
         ENDIF         
 20   CONTINUE
C     REPEAT THE ABOVE STEPS UNTIL NO MORE AREA 0 TRIANGLES 
C     ARE LEFT ON THE BORDER.
      DO 30 IT = 1,NTRI
         IF ((LCCC(3,IT).EQ.0) .AND. 
     .       (ICCC(1,IT).EQ.0 .OR. 
     .        ICCC(2,IT).EQ.0 .OR. 
     .        ICCC(3,IT).EQ.0)) THEN
            GO TO 11
         ENDIF
 30   CONTINUE
      RETURN

      END

      SUBROUTINE RMSHNB(I0,I1,I2,LCCC,ICCC)
      INTEGER I0,I1,I2,K,ICCC(6,*) 
      DOUBLE PRECISION LCCC(14,*)
C     REMOVE A SHARED NEIGBOUR TRIANGLE I0 FROM THE NEIGBOUR LISTS
C     OF TRIANGLE I1 AND I2 (REPLACE OCCURENCE OF I0 IN LIST OF I1 WITH I2
C     AND VICE VERSA.)
      DO 1 K=1,3
         IF(I1 .GT. 0) THEN 
            IF (ICCC(K,I1) .EQ. I0) ICCC(K,I1)=I2
         ENDIF
         IF(I2 .GT. 0) THEN 
            IF (ICCC(K,I2) .EQ. I0) ICCC(K,I2)=I1
         ENDIF
 1    CONTINUE
C     MARK TRIANGLE I0 AS REMOVED (AREA= -1)
      LCCC(3,I0)=-1
      RETURN
      END
