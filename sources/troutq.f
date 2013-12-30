      SUBROUTINE TROUTQ (NCC,LCC,N,X,Y,LIST,LPTR,LEND,LOUT,
     .                   NNABS,NPTR,NPTR1,NABOR,NBNOS,NA,NB,NT)
      INTEGER NCC, LCC(*), N, LIST(*), LPTR(*), LEND(N),LOUT,
     .        NNABS(*),NPTR(*),NPTR1(*),NABOR(*),NBNOS(*),NA,NB,NT
      DOUBLE PRECISION    X(N), Y(N)
C
C***********************************************************
C
C                                               From TRIPACK
C                                            Robert J. Renka
C                                  Dept. of Computer Science
C                                       Univ. of North Texas
C                                             (817) 565-2767
C                                                   08/22/91
C
C   modified version of TRPRNT
C
C   further modification Roger.Bivand@nhh.no 26 April 2000
C
C   Given a triangulation of a set of points in the plane,
C this subroutine returns the adjacency lists and, option-
C ally, the nodal coordinates.  The
C list of neighbors of a boundary node is followed by index
C 0.  The numbers of boundary nodes, triangles, and arcs,
C and the constraint curve starting indexes, if any, are
C also returned.
C
C
C On input:
C
C       NCC = Number of constraints.
C
C       LCC = List of constraint curve starting indexes (or
C             dummy array of length 1 if NCC = 0).
C
C       N = Number of nodes in the triangulation.
C           3 .LE. N .LE. 9999.
C
C       X,Y = Arrays of length N containing the coordinates
C             of the nodes in the triangulation -- not used
C             unless PRNTX = TRUE.
C
C       LIST,LPTR,LEND = Data structure defining the trian-
C                        gulation.  Refer to subroutine
C                        TRMESH.
C
C
C None of these parameters are altered by this routine.
C
C Output
C  
C     NNABS number of neighbours
C     NPTR start pointers
C     NPTR1 end pointers
C     NABOR list of neighbours
C     NBNOS boundary node indices
C     NB number of boundary nodes
C     NA number of arcs
C     NT number of triangles
C
C Modules required by TRPRNT:  None
C
C***********************************************************
C
      INTEGER I, INC, K, KK, LP, LPL, LUN, 
     .        ND, NL, NLMAX, NMAX, NODE, NN
      DATA  NMAX/9999/,  NLMAX/60/
C
      NN = N
      LUN = LOUT
      IF (LUN .LT. 0  .OR.  LUN .GT. 99) LUN = 6
C
C test the range of N.
C
C     IF (NN .LT. 3  .OR.  NN .GT. NMAX) THEN
C
C N is outside its valid range.
C
c        WRITE (LUN,110)
C       GO TO 5
C     ENDIF
C
C Initialize NL (the number of lines printed on the current
C   page) and NB (the number of boundary nodes encountered).
C
      NL = 6
      NB = 0

C
C Print X, Y, and LIST.
C
        KK = 0
        NPTR(1) = KK + 1
        DO 4 NODE = 1,NN
          LPL = LEND(NODE)
          LP = LPL
          K = 0
    3     K = K + 1
          KK = KK + 1
          LP = LPTR(LP)
          ND = LIST(LP)
          NABOR(KK) = ND
          IF (LP .NE. LPL) GO TO 3

          IF (ND .LE. 0) THEN
            NABOR(KK) = -ND
            K = K + 1
            NB = NB + 1
            NBNOS(NB) = NODE
          ENDIF

          INC = (K-1)/8 + 2
          NL = NL + INC
          IF (NL .GT. NLMAX) THEN
            NL = INC
          ENDIF
          NPTR1(NODE) = KK
          IF (NODE .LT. NN) NPTR(NODE+1) = NPTR1(NODE) + 1
          NNABS(NODE) = NPTR1(NODE) - NPTR(NODE) + 1
    4   CONTINUE
C
C Print NB, NA, and NT (boundary nodes, arcs, and
C   triangles).
C
      NT = 2*NN - NB - 2
      NA = NT + NN - 1
C
C Print NCC and LCC.
C
C5    RETURN
C
C Print formats:
C
C 110 FORMAT (1X,10X,'*** N IS OUTSIDE ITS VALID',
C    .        ' RANGE ***')
      END

