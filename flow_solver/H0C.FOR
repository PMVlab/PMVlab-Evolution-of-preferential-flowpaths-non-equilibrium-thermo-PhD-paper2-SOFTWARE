      program h0c
C***********************************************************************
c
c     This program solves the steady state flow in a heterogeneus aqifer
c     with a single pumping well at a particular node in the domain.
c     You can place up to 20 pumping/injection nodes in the field.
c
c     OUTPUT: head field correspondig to the given Boundary Conditons
c             and LOG-transmissivity in input.
c     The output will be a binary file of the heads at the nodes of the grid.
c
C     SUBROUTINE SOLV. SOLVES THE LINEAR SYSTEM MEMORIZED
C     AS ISOLATED ELEMENTS WITH INCOMPLETE CHOLESKY 
C     CONJUGATE GRADIENT (ICCG)
c
c ---
c     BOUNDARY CONDITIONS:
c     1) Dirichlet boundary along vertical sides (fixed head)
c     2) No-flow boundary along horizontal sides.
c ---
c
c     NN         Number of blocks with constant T along X-axis
c     MM         Number of blocks with constant T along X-axis
c     RRX        Number of segments into which a block is subdivided
c                (along X-axix)
c     RRY        Number of segments into which a block is subdivided
c                (along Y-axix) 
c     NNX        Domain sides lenght in terms of correlation length 
c                (along X-axis)
c     NNY        Domain sides lenght in terms of correlation length 
c                (along Y-axis)
c     NNI        Number of nodes along X-axis
c     NNJ        Number of nodes along Y-axis
c     NNPUNTI    Number of unknown nodes
c     NTK        Number of NON-ZERO Elements of the A matrix; A is symmetric,
c                so we store only the main diagonal and the upper part. In this case
c                NTK=1/2*(4x4+( 2*(NNI-4)+2*(NJ-2) )*6 + (NI-4)*(NJ-2)*9 + NNPUNTI)
c     TK(NTK)    Contains the symmetric part of the elements of matrix A
c     V(NNPUNTI) It is the right hand side of the system of equations
c
c ---   Grid Specifications
c
c     DX         Spacing along X-dir. (in terms of correlation scales)
c     DY         Spacing along Y-dir. (in terms of correlation scales)
c
c ---   Domain Properties
c
c     T(NN*MM)   Vector of known elements trasmissivity
c     Y(NN*MM)   Vector of known elements log-trasmissivity
c
c     As input, from the subroutine LEGGIT you read a binary file
c     which contains the logconductivity distribution. See the subroutine for
c     additional comments.
c
c --- Pumping rates and wells location
c
c     NPOZZI     it is set = 1; it is the number of wells in the domain
c     Q(20)      it the pumping rate at the specified well. +1 means
c                extraction, while -1 means injection.
c     XP(20)     X location of the pumping wells
c     YP(20)     Y location of the pumping wells
c
c     See the subroutine pump for ref.
c
c ---    HEAD SPECIFICATIONS:
c 
c     FI(NNPUNTI) Vector of unknown head values
c     FI0(NJ)     Vector of known heads along right hand side of domain
c     FI1(NJ)     Vector of known heads along right left side of domain
c
c     For the Dirichlet boundary see sub READBOUND
c     The output is given by subroutine SALVA. Ref. for furhter comments.
c
C
C  VENGONO ASSEGNATI I VALORI DEI PARAMETRI FONDAMENTALI E CALCOLATI
C  QUELLI DEI PARAMETRI DIPENDENTI; CON TALI VALORI VENGONO DIMENSIONATI
C  I VETTORI CHE SI USERANNO I SEGUITO.
C
C*************************************************************************
      PARAMETER(NN=150,MM=60,RRX=1,RRY=1,NNX=NN/5,NNY=MM/5)
      PARAMETER(NNI=NN*RRX+1,NNJ=MM*RRY+1,NNPUNTI=NNJ*(NNI-2))  
      PARAMETER(NTK=5*NNPUNTI-3*(NNI+NNJ-3)-1)
      REAL*8 TK(NTK),V(NNPUNTI),UU(NNPUNTI),VV(NNPUNTI),
     &       RR(NNPUNTI),AF(NTK),FI(NNPUNTI),ERRMAX,REALERR
      REAL*8 FI0(MM+1),FI1(MM+1)
      REAL*8 T(MM*NN),Q(20),XP(20),YP(20)
      real*8 dx,dy
c
      INTEGER RX,RY,IPUNT(NNPUNTI+1)
c
      INTEGER INDRIG(NTK),NPZ(20)
c
      CHARACTER*4 CHT
c
c --- Begin loop for well location over Grid nodes
c
      do 1111 inmc=1,100

c       iht1=0
       write(cht,'(i4.4)')inmc
       write(*,*)cht
c
c --- Grid specifications
c
      RX=RRX
      RY=RRY
      NI=NNI
      NJ=NNJ
      NPUNTI=NNPUNTI
      N=NN
      M=MM
      DX=(1.d0*NNX)/(NN*RRX)
      DY=(1.d0*NNY)/(MM*RRY)    
c      WRITE(*,*)'DX=',DX,'   DY=',DY
c
c --- reads log-transmissivity in input
c
      CALL LEGGIT(M,N,T,CHT)
c
c --- reads Dirichlet Boundary conditions along the two vertical sides 
c
      call readbound(nj,fi1,fi0,AJ,NNX)
c
c --- Reads location and strenght of pumping well
c
      call pump(inmc,jnmc,npozzi,XP,YP,Q,NPZ,DX,DY,NJ,NNX,NNY)
c
c --- Sets Error Criterion for computations
c
      ERRMAX=DX*AJ*1E-11
C        WRITE(*,*)'FI0=',FI0,'   FI1=',FI1
c      WRITE(*,*)'ERRMAX=',ERRMAX
c
c --- Builds up the matrix
c 
      CALL FAITK(NPUNTI,M,RX,RY,DX,DY,NI,NJ,FI0,FI1,T,TK,INDRIG,
     &  IPUNT,V,Q,NPZ,NPOZZI)

      MINMAT=0
      DO 5 I=1,NTK
5       MINMAT=MINMAT+ABS(TK(I))
      MINMAT=MINMAT/NTK
      REALERR=ERRMAX
c      WRITE(*,*)'REALERR=',REALERR
c
c --- Solve the system of equations
c
      CALL SOLUTION(ERR1,ERR2,TK,INDRIG,IPUNT,V,FI,NPUNTI,UU,VV,RR,
     &  AF,FI0,REALERR,ITE,NI,NJ,DX)
c
c --- Saves the head values into a binary file
c
      CALL SALVA(FI,M,N,NI,NJ,DX,DY,FI0,FI1,ERRMAX,CHT,NPOZZI,XP,YP,Q)


      WRITE(*,1)N,M,NPUNTI,ITE
C1       FORMAT(' N=',I3,'  M=',I2,'   NPUNTI=',I4,'    N. ITER.=',I4)
1       FORMAT(' N=',I3,'  M=',I3,'   NPUNTI=',I6,'    N. ITER.=',I4)
C---------------
C       OPEN(UNIT=1,FILE='FIP'//CHT//'.DAT',STATUS='UNKNOWN')
C       DO 20 I=1,NI
C       DO 20 J=1,NJ
C       X=(I-1)*DX
C       Y=(J-1)*DY
C       IF(I.EQ.1) THEN
C         WRITE(1,*) X,Y,FI0
C       ELSEIF(I.EQ.NI) THEN
C         WRITE(1,*) X,Y,FI1
C       ELSE
C         K=NJ*(I-1)+J-NJ
C         RFI=FI(K)
C         WRITE(1,*) X,Y,RFI
C       ENDIF
C20     CONTINUE
C       CLOSE(1)
C----------------

1111  continue

      STOP
      END

C*************************************************
c
      subroutine readbound(nj,fi1,fi0,AJ,NNX)
      real*8 fi1(*),fi0(*),fip1,fip0              
      
c        WRITE(*,*)'GRADIENTE MEDIO J = '
c        READ(*,*) AJ
      AJ=0.01
c     WRITE(*,*)'GRADIENTE MEDIO J = ',AJ
c
c     Left hand side boundary
c
      FIp0=1.0d2
c          1.0d-3 
c		 5.d-2
c     Right hand side boundary
c
      FIp1=0.d0
       
c        FIp1=AJ*NNX*1.d0
c        FIp0=0.d0  
c        write(*,*)fip1
c	read(*,*)
C******* APERTURA FILES FI AL CONTORNO ***********
      DO 130 KK=1,NJ
130     FI0(KK)=FIP0
C130     WRITE(*,*)'IN1=',IN1,'FI0=',FI0(IN1)
C        WRITE(*,*)'BATTI UN TASTO......'
C130     READ(*,*)
      DO 131 KK=1,NJ
131     FI1(KK)=FIP1
C        WRITE(*,*)'IM1=',IM1,'FI1=',FI1(IM1)
C        WRITE(*,*)'BATTI UN TASTO......'
C131     READ(*,*)
      return
      end
C************************************************************
      subroutine pump(inmc,jnmc,npozzi,XP,YP,Q,NPZ,DX,DY,NJ,NNX,NNY)
      real*8 XP(*),YP(*),Q(*),DX,DY
      
      integer NPZ(*)
       
      NPOZZI=1
C
C
       XP(1)=NNX/2.
       YP(1)=NNY/2.

c      XP(1)=DX*(inmc-1)
c      YP(1)=(jnmc-1)*DY
      
      write(*,*)inmc,jnmc
      write(*,*)'XP:',XP(1),'YP:',YP(1)
c      read(*,*)
      
      Q(1)=0.d0

      WRITE(*,*)'Q =',Q(1)
c      read(*,*)

      DO 10 L=1,NPOZZI
c        K=(Inmc-2)*NJ+Jnmc
        I=(XP(L)/DX+1.5)                
        J=(YP(L)/DY+1.5)                
        K=(I-2)*NJ+J
c        write(*,*)'K =',K
c        read(*,*)
        NPZ(L)=K
10    CONTINUE

      return
      end

C************************************************************
c CALL LEGGIT(M,N,T,CHT)
      SUBROUTINE LEGGIT(M,N,T,CHT1)
      REAL*8 T(*),Y(150*60)
      real*8 sumt,sumlog,PK,SK,VARK,TGEOM,P,S,YVAR
      CHARACTER*4 CHT1
c        OPEN(UNIT=11,FILE='TRS0001.txt')
        OPEN(UNIT=11,FILE='TRS'//CHT1//'.txt')
c	write(*,*)CHT1
c	write(*,*)'press a key'
c	read(*,*)
c        read(11,*)
c        read(11,*)
c        read(11,*)
c        OPEN(UNIT=11,FILE='TRS0001.BIN',STATUS='OLD',
c     &  ACCESS='DIRECT',RECL=4)
c
c      write(*,*)n,m
c
      sumlog=0.
      sumt=0.
      DO 10 IND=1,N*M
c        READ(11,REC=IND+7) YY
c        READ(11,*) it,YY    
        READ(11,*) YY
        Y(IND)=dble(yy) 
        T(IND)=DEXP(Y(IND))
c	write(*,*)ind,YY,T(ind)
c     Queste due righe servono soltanto se voglio fare un campo omogeneo in K=KG
c        T(ind)=1.d0
cc        Y(ind)=log(8.64)
ccc   &&&&&
        sumlog=sumlog+Y(IND)
        sumt=sumt+T(IND)
10    CONTINUE        
c      
      sumlog=sumlog/(N*M*1.)
      sumt=sumt/(N*M*1.)
      TGEOM=DEXP(sumlog)
      PK=0.d0
      SK=0.d0
      VARK=0.d0
      P=0.d0
      S=0.d0
      YVAR=0.d0
      do 11 i=1,N*M
        S=Y(i)-sumlog
        SK=T(i)-sumt
        P=S*S       
        PK=SK*SK
        YVAR=YVAR+P
        VARK=VARK+PK
11    continue

      YVAR=YVAR/(N*M*1.d0)
      VARK=VARK/(N*M*1.d0)
c      write(*,*)YVAR,SUMLOG
c      read(*,*)      
c      WRITE(*,*)'Y ave   NUM   = ',sumlog
c      WRITE(*,*)'K geom. NUM   = ',TGEOM
c      write(*,*)'Y VAR   NUM   = ',YVAR
c      WRITE(*,*)'K arit. NUM   = ',sumt
c      write(*,*)'K arit. TEO   = ',TGEOM*exp(YVAR/2.)
c      write(*,*)'Var. K  NUM   = ',VARK
c      write(*,*)'Var. K  TEO   = ',(DEXP(YVAR)-1.d0)*sumt*sumt
      
      CLOSE(11)

      RETURN
      END
c
C*****************************************************************
      SUBROUTINE FAITK(NPUNTI,M,RX,RY,DX,DY,NI,NJ,FI0,FI1,T,TK,
     &  INDRIG,IPUNT,V,Q,NPZ,NPOZZI)
C*****************************************************************
      REAL*8 TK(*),V(*)
      REAL*8 FI0(*),FI1(*)
      REAL*8 T(*),TE(4),TKTEMP(150,4),Q(*)
      REAL*8 DX,DY,DYX,C1,C2,C3
      INTEGER RX,RY,IPUNT(*)
      INTEGER INDRIG(*),ITEMP(150),IRIGTEMP(150,4),NPZ(20)
c
c --- These are the coefficients in the A matrix
c
      DYX=DY/DX
c
c --- This is, for example node 1 speaking with node 3 in an element
c --- - 2*C1 is the node speaking with itself
c
      C1=-(DYX+1.d0/DYX)/6.d0
c
c --- This is, for example node 1 speaking with node 4 in an element
c
      C2=(DYX-2.d0/DYX)/6.d0
c
c --- This is, for example node 1 speaking with node 2 in an element
c
      C3=(1.d0/DYX-2.d0*DYX)/6.d0
c      
c      write(*,*)'DYX:',DYX
c      write(*,*)'C1:',C1
c      write(*,*)'C2:',C2
c      write(*,*)'C3:',C3
c      read(*,*)
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&& CHECK PASSAGGIO DATI
c        IP0=0
c        IP1=0
c        DO 132 II=1,NJ
c        IP0=IP0+1
c        WRITE(*,*)'IP0=',IP0,'FI0=',FI0(IP0)
c        WRITE(*,*)'BATTI UN TASTO......'
c132     READ(*,*)
c        DO 133 KK=1,NJ
c        IP1=IP1+1
c        WRITE(*,*)'IP1=',IP1,'FI1=',FI1(IP1)
c        WRITE(*,*)'BATTI UN TASTO......'
c133     READ(*,*)
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      INDTK=0  
c      IN=0
c      IM=0
c
c --- Loop over the points of the unknowns
c      
      DO 10 K=1,NPUNTI
c
c --- Identifies column I and row J of the Kth node
c     in the original grid. For example, when K=1,
c     the node is in 2nd column, 1st row, and so on.
c     The numbering of the nodes is first by y and
c     then by x. The nodes are numbered consecutively
c --- starting from the lowest left one upward.
c
        I=(K-1)/NJ+2
        J=K-((K-1)/NJ)*NJ
        NUMNODO=(I-1)*NJ+J
c
c --- Identifies the transmissivity of the four elements
c     that share the node
c  
        CALL FAITE(I,J,M,RX,RY,TE,T)

        IF(K.EQ.1) GO TO 25

        INDR=MOD(K-1,NJ+1)+1
        ITEST=0
c
c --- Redirects the components of the T() matrix
c
        DO 20 L=1,4
C       IF(TKTEMP(INDR,L).EQ.0) GOTO 20
          IF(TKTEMP(INDR,L).LT.0.OR.TKTEMP(INDR,L).GT.0) THEN
            GOTO 15
          ELSE
            GOTO 20
          ENDIF   
15        INDTK=INDTK+1
C         INDTK=INDTK+1
        TK(INDTK)=TKTEMP(INDR,L)
        INDRIG(INDTK)=IRIGTEMP(INDR,L)
        TKTEMP(INDR,L)=0
        IF(ITEST.EQ.0) THEN
          IPUNT(K)=INDTK
          ITEST=1
         ENDIF

20      CONTINUE
        
        ITEMP(INDR)=0

25      INDTK=INDTK+1
c
c --- This is node K which speaks with itself and it is 
c     stored in the location INDTK of vector TK(INDTK)
c
        TK(INDTK)=-2*(TE(1)+TE(2)+TE(3)+TE(4))*C1
        INDRIG(INDTK)=K
        IF(ITEST.EQ.0) IPUNT(K)=INDTK
c
c --- Sets to zero the correspondent line in the right hand side vector
c
        V(K)=0.d0
c
C ****** Dirichlet Boundary *****************
c     Assignes the values of the V() vector
c --- which is the right hand side of the problem                                       
c
c --- Corner node at the bottom left
c
      IF(I.EQ.2.AND.J.EQ.1)THEN
        I1=1
        I2=2
        V(K)=-FI0(I1)*(TE(3)+TE(4))*C3
     &       -FI0(I2)*(TE(3)+TE(4))*C1      
c        WRITE(*,*)'I=',I,'I1=',I1,'I2=',I2,'V(K)=',K,V(K)
C        WRITE(*,*)'BATTI UN TASTO.....'
C        READ(*,*)
      ENDIF
c
c --- Corner node at the upper left
c
      IF(I.EQ.2.AND.J.EQ.NJ)THEN
        I1=NJ
        I2=NJ-1
        V(K)=-FI0(I1)*(TE(3)+TE(4))*C3
     &       -FI0(I2)*(TE(3)+TE(4))*C1      
c        WRITE(*,*)'I=',I,'I1=',I1,'I2=',I2,'V(K)=',K,V(K)
C        WRITE(*,*)'BATTI UN TASTO.....'
C        READ(*,*)
      ENDIF
c
c --- Corner node at the bottom right
c
      IF(I.EQ.(NI-1).and.j.eq.1)THEN
        I1=1
        I2=2
        V(K)=-FI1(I1)*(TE(1)+TE(2))*C3
     &       -FI1(I2)*(TE(1)+TE(2))*C1
c        WRITE(*,*)'I=',I,'I1=',I1,'I2=',I2,'V(K)=',K,V(K)
C        WRITE(*,*)'BATTI UN TASTO.....'
C        READ(*,*)
      ENDIF
c
c --- Corner node at the upper right
c
      IF(I.EQ.(NI-1).and.j.eq.NJ)THEN
        I1=NJ
        I2=NJ-1
        V(K)=-FI1(I1)*(TE(1)+TE(2))*C3
     &       -FI1(I2)*(TE(1)+TE(2))*C1
c        WRITE(*,*)'I=',I,'I1=',I1,'I2=',I2,'V(K)=',K,V(K)
C        WRITE(*,*)'BATTI UN TASTO.....'
C        READ(*,*)
      ENDIF    
c
c --- Nodes along the left side (Y-DIR.)
c
      if(i.eq.2.and.j.gt.1.and.j.lt.nj)then
        I1=J-1
        I2=I1+1
        I3=I2+1
        V(K)=-(C3*(TE(3)+TE(4))*FI0(I2)+
     $           C1*(TE(3)*FI0(I1)+TE(4)*FI0(I3)))
      endif
c
c --- Nodes along the right side (X-DIR.)
c
      if(i.eq.(NI-1).and.j.gt.1.and.j.lt.nj)then
        I1=J-1
        I2=I1+1
        I3=I2+1
        V(K)=-(C3*(TE(1)+TE(2))*FI1(I2)+
     $           C1*(TE(1)*FI1(I3)+TE(2)*FI1(I1)))
      endif      
c
C*****************************************************
c
c --- Take into account the presence of wells into the right hand side
c
      DO 5 L=1,NPOZZI
        IF(K.EQ.NPZ(L)) V(K)=V(K)-Q(L)
5     CONTINUE
c
c --- Compute the non-zero elements in the matrix
c     Only 4 terms are computed because of symmetry:
c     1. Node above node K
c     2. Node shifted 1 column to the right and 1 location down
c     3. Node shifted 1 column to the right at the same location as K
c     4. Node shifted 1 column to the right and 1 location up
c
c --- 1. This is the node above node K
c
        IL=K+1
        IF(IL.GT.NPUNTI) GOTO 10
        INDCL=MOD(IL-1,NJ+1)+1
        ITEMP(INDCL)=ITEMP(INDCL)+1
        TKTEMP(INDCL,ITEMP(INDCL))=(TE(1)+TE(4))*C2
        IRIGTEMP(INDCL,ITEMP(INDCL))=K
c
c --- 2. This is the node shifted 1 column to the right and 1 location down
c
        IL=K+NJ-1
        IF(IL.GT.NPUNTI) GOTO 10
        INDCL=MOD(IL-1,NJ+1)+1
        ITEMP(INDCL)=ITEMP(INDCL)+1
        TKTEMP(INDCL,ITEMP(INDCL))=TE(2)*C1
        IRIGTEMP(INDCL,ITEMP(INDCL))=K
c
c --- 3. This is the Node shifted 1 column to the right at the same location as K
c
        IL=K+NJ
        IF(IL.GT.NPUNTI) GOTO 10
        INDCL=MOD(IL-1,NJ+1)+1
        ITEMP(INDCL)=ITEMP(INDCL)+1
        TKTEMP(INDCL,ITEMP(INDCL))=(TE(1)+TE(2))*C3
        IRIGTEMP(INDCL,ITEMP(INDCL))=K
c
c --- 4. This is the node shifted 1 column to the right and 1 location up
c
        IL=K+NJ+1
        IF(IL.GT.NPUNTI) GOTO 10
        INDCL=MOD(IL-1,NJ+1)+1
        ITEMP(INDCL)=ITEMP(INDCL)+1
        TKTEMP(INDCL,ITEMP(INDCL))=TE(1)*C1
        IRIGTEMP(INDCL,ITEMP(INDCL))=K

10    CONTINUE

      IPUNT(NPUNTI+1)=INDTK+1
      
      RETURN
      END


C************************************************************** 
      SUBROUTINE FAITE(I,J,M,RX,RY,TE,T)
      REAL*8 TE(4),T(*),CI(4),CJ(4)
      INTEGER RX,RY

      DATA CI/.5,.5,-.5,-.5/
      DATA CJ/.5,-.5,-.5,.5/
c
c --- Location on the grid taking into account the
c     number of segments for each element
c
      BX=(I-1.)/RX+1. 
      BY=(J-1.)/RY+1. 
c
c --- Each node shares 4 elements.
c
      DO 10 L=1,4
        IB=BX+CI(L)/RX
        JB=BY+CJ(L)/RY
        IF(JB.GE.1.AND.JB.LE.M) THEN
          TE(L)=T((IB-1)*M+JB)
        ELSE
c
c --- nodes on the no-flux boundaries have
c --- two elements. 
c
          TE(L)=0.d0 
        ENDIF
10    CONTINUE

      RETURN
      END

C****************************************************
      SUBROUTINE SOLUTION(ERR1,ERR2,TK,INDRIG,IPUNT,V,FI,NPUNTI,
     &  UU,VV,RR,AF,FI0,ERRMAX,ITE,NI,NJ,DX)

C       REAL*8 TK(*),V(*),FI(*),ERRMAX,UU(*),VV(*),RR(*),AF(*),FIX
C     &  ,DERR2

      REAL*8 TK(*),V(*),FI(*),ERRMAX,UU(*),VV(*),RR(*),AF(*)
      REAL*8 FI0(*),DX
      INTEGER IPUNT(*)
      INTEGER INDRIG(*)

C       K=0
C       ALX=(NI-1)*DX
C
C       DO 5 I=2,NI-1
C       X=(I-1)*DX
C       FIX=FI0*(1-X/ALX)
C       DO 5 J=1,NJ
C       K=K+1
C       FI(K)=FIX
C5      CONTINUE



      CALL SOLV(TK,INDRIG,IPUNT,V,FI,UU,VV,RR,AF,NPUNTI,ERRMAX,ITE)

C       ERR1=0
C       DERR2=0
C       K=0

C       DO 10 I=2,NI-1
C       X=(I-1)*DX
C       FIX=FI0*(1-X/ALX)
C       DO 10 J=1,NJ
C       K=K+1
C       ERR1=ERR1+(FI(K)-FIX)**2
C       DERR2=DMAX1(DERR2,DABS(FI(K)-FIX))
C10     CONTINUE
C       ERR2=DERR2
C       IF(ERR1.GT.0) ERR1=SQRT(ERR1)
	RETURN
	END


C*****************************************************************
C
C     SUBROUTINE SOLV. RISOLVE IL SISTEMA LINEARE MEMORIZZATO
C     AD ELEMENTI ISOLATI CON IL METODO CHOLESKY INCOMPLETO
C     GRADIENTE CONIUGATO (ICCG)
C
	SUBROUTINE SOLV(TK,INDRIG,IPUNT,AL,ELRE,UU,VV,RR,AF,N,TESTM
     &  ,ITE)
	IMPLICIT DOUBLE PRECISION (A-H,O-Z)
	INTEGER IPUNT(*)
	INTEGER INDRIG(*)
	REAL*8  TK(*),AF(*),AL(*),ELRE(*),UU(*),VV(*),RR(*)

C        NITEM=N
	NITEM=1100
	LX=0
      
      CALL ICDCGF(TK,AF,INDRIG,IPUNT,VV,N)
      
      CALL ICDCGS(TK,AF,INDRIG,IPUNT,AL,ELRE,TESTM,NITEM,ITE,RR,UU,
     1  VV,N,LX)
      
      RETURN
      END


C********************************************************************
C
C     SUBROUTINE ICDCGF
C
	SUBROUTINE ICDCGF(AA,AF,LL,KK,VV,N)
	IMPLICIT DOUBLE PRECISION (A-H,O-Z)
	INTEGER KK(*)
	INTEGER LL(*)
	REAL*8 AA(*),AF(*),VV(*)
	K=KK(N+1)-1
	DO 1 I=1,K
	AF(I)=AA(I)
1       CONTINUE
	CALL ILDLTF(AF,LL,KK,VV,N)
	RETURN
	END
C*****************************************************************
C
C  
	SUBROUTINE ICDCGS(AA,AF,LL,KK,BB,XX,RERMAX,ITEMAX,ITE,RR,UU,
     1  VV,N,LX)
	IMPLICIT DOUBLE PRECISION (A-H,O-Z)
	INTEGER KK(*)
	INTEGER LL(*)
	REAL*8 AA(*),AF(*),BB(*),XX(*),RR(*),UU(*),VV(*)
	CALL MULMAT(AA,LL,KK,XX,RR,N)
	DO 10 I=1,N
	RR(I)=RR(I)-BB(I)
	VV(I)=RR(I)
10      CONTINUE
	CALL ILDLTS(AF,LL,KK,VV,N)
	S=0.D0
	DO 20 I=1,N
	S=S+RR(I)*VV(I)
20      CONTINUE
	RV=S
	DO 30 I=1,N
	UU(I)=VV(I)
30      CONTINUE
	ITE=0
40      IF(LX.EQ.0)GOTO 50
	R=DABS(RR(LX))
	X=DABS(XX(LX))
	GOTO 70
50      R=0.D0
	X=0.D0
	DO 60 I=1,N
	R=DMAX1(R,ABS(RR(I)))
	X=DMAX1(X,ABS(XX(I)))
60      CONTINUE
70      IF(R.LE.X*RERMAX.OR.ITE.GE.ITEMAX) RETURN
	ITE=ITE+1
C----- Check on the number of iterations
c	WRITE(*,*) 'ITE=',ITE
C-----
	CALL MULMAT(AA,LL,KK,UU,VV,N)
	S=0.D0
	DO 80 I=1,N
	S=S+UU(I)*VV(I)
80      CONTINUE
	UV=S
	RVO=RV
	SRV=-RVO/UV
	DO 90 I=1,N
	XX(I)=XX(I)+SRV*UU(I)
	RR(I)=RR(I)+SRV*VV(I)
	VV(I)=RR(I)
90      CONTINUE
	CALL ILDLTS(AF,LL,KK,VV,N)
	S=0.D0
	DO 100 I=1,N
	S=S+RR(I)*VV(I)
100     CONTINUE
	RV=S
	SRV=RV/RVO
	DO 110 I=1,N
	UU(I)=VV(I)+SRV*UU(I)
110     CONTINUE
	GO TO 40
	END

C********************************************************************
C
C     SUBROUTINE ILDLTF
C
	SUBROUTINE ILDLTF(AA,LL,KK,BB,N)
	IMPLICIT DOUBLE PRECISION (A-H,O-Z)
	INTEGER KK(*)
	REAL*8 AA(*),BB(*)
	INTEGER LL(*)
	DO 40 K=1,N
	K1=KK(K)
	K2=KK(K+1)-1
	DO 10 L=1,K
	BB(L)=0.D0
10      CONTINUE
	DO 40 L=K1,K2
	LJ=LL(L)
	KS1=KK(LJ)
	KS2=KK(LJ+1)-2
	S=AA(L)
	IF(KS1.GT.KS2)GOTO 30
	DO 20 LS=KS1,KS2
	LSRV=LL(LS)
	S=S-AA(LS)*BB(LSRV)
20      CONTINUE
30      AA(L)=S
	BB(LJ)=S/AA(KS2+1)
40      CONTINUE
	DO 50 K=1,N
	BB(K)=0.D0
50      CONTINUE
	RETURN
	END

C*************************************************************
C  
C 
	SUBROUTINE ILDLTS(AA,LL,KK,BB,N)
	IMPLICIT DOUBLE PRECISION (A-H,O-Z)
	INTEGER LL(*)
	INTEGER KK(*)
	REAL*8 AA(*),BB(*)
	DO 120 K=1,N
	K1=KK(K)
	K2=KK(K+1)-2
	S=BB(K)
	IF(K1.GT.K2)GOTO 110
	DO 100 L=K1,K2
	LJ=LL(L)
	S=S-AA(L)*BB(LJ)
100     CONTINUE
110     BB(K)=S/AA(K2+1)
120     CONTINUE
	DO 130 K=1,N
	KS=KK(K+1)-1
	BB(K)=BB(K)*AA(KS)
130     CONTINUE
	DO 150 K=N,1,-1
	K1=KK(K)
	K2=KK(K+1)-2
	BB(K)=BB(K)/AA(K2+1)
	IF(K1.GT.K2)GOTO 150
	BSRV=BB(K)
	DO 140 L=K2,K1,-1
	LJ=LL(L)
	BB(LJ)=BB(LJ)-BSRV*AA(L)
140     CONTINUE
150     CONTINUE
	RETURN
	END

C********************************************************************
C 
C     SUBROUTINE MULMAT
C
	SUBROUTINE MULMAT(AA,LL,KK,BB,CC,N)
	IMPLICIT DOUBLE PRECISION (A-H,O-Z)
	INTEGER KK(*)
C       INTEGER *2 LL(*)
	INTEGER LL(*)
	REAL*8 AA(*),BB(*),CC(*)
	DO 50 I=1,N
50      CC(I)=0.D0
	I=0
	DO 20 K=1,N
	KS=KK(K+1)-KK(K)
	J=I
	S=0.D0
	DO 10 L=1,KS
	J=J+1
	LJ=LL(J)
	S=S+AA(J)*BB(LJ)
10      CONTINUE
	CC(K)=S
	I=I+KS
20      CONTINUE
	DO 40 K=N,1,-1
	KS=KK(K+1)-KK(K)
	J=I
	I=I-KS
	IF(KS.LE.1)GOTO 40
	KS1=KS-1
	DO 30 L=KS1,1,-1
	J=J-1
	LJ=LL(J)
	CC(LJ)=CC(LJ)+AA(J)*BB(K)
30      CONTINUE
40      CONTINUE
	RETURN
	END

C***************************************************************
      SUBROUTINE SALVA(FI,M,N,NI,NJ,DX,DY,FI0,FI1,ERRMAX,CHT,
     &  NPOZZI,XX,YY,QQ)
C
C***************************************************************
      REAL*8 FI(*),ERRMAX,FITEMP
      REAL*8 FI0(*),FI1(*)
      REAL*8 XX(*),YY(*),QQ(*),DX,DY
      CHARACTER*4 CHT
      IN=0
      IM=0
	write(*,*)cht
c --- writes output to ascii file
c ---
      OPEN(UNIT=12,FILE='FIP'//CHT//'.txt',STATUS='UNKNOWN')
c ---
c	write(12,*)'Hydraulic Head field'
c      write(12,*)N,M,NI,NJ,DX,DY
c	write(12,*)'Node No.      Head'
	IND=0
	IN=0
	IM=0
      DO 10 I=1,NI
      DO 20 J=1,NJ
      IND=IND+1
      IF(I.EQ.1) THEN
        IN=IN+1
        FITEMP=FI0(IN)
c        WRITE(12,*) IND,FITEMP
        WRITE(12,*) FITEMP
        GO TO 20
      ENDIF
      IF(I.EQ.NI) THEN
        IM=IM+1
        FITEMP=FI1(IM)
c        WRITE(12,*) IND,FITEMP
        WRITE(12,*) FITEMP
      GO TO 20
      ENDIF
c      WRITE(12,*) IND-8,FI(IND-NJ-LAST)
c      WRITE(12,*) IND,FI(IND-NJ)
      WRITE(12,*) FI(IND-NJ)
20      CONTINUE
10      CONTINUE


c --- writes output to binary file
c ---

c      OPEN(UNIT=12,FILE='FIP'//CHT//'.BIN',STATUS='UNKNOWN',
c     &  ACCESS='DIRECT',RECL=8)
      
c      WRITE(12,REC=1) N
c      WRITE(12,REC=2) M
c      WRITE(12,REC=3) NI
c      WRITE(12,REC=4) NJ
c      WRITE(12,REC=5) DX
c      WRITE(12,REC=6) DY
c      RERRMAX=ERRMAX
c      WRITE(12,REC=7) RERRMAX
c      WRITE(12,REC=8) NPOZZI
c      IND=8
c      DO 1 I=1,NPOZZI
c      WRITE(12,REC=IND+1) XX(I)
c      WRITE(12,REC=IND+2) YY(I)
c      WRITE(12,REC=IND+3) QQ(I)
c1       IND=IND+3       
c      LAST=IND
c      DO 10 I=1,NI
c      DO 20 J=1,NJ
c      IND=IND+1
c      IF(I.EQ.1) THEN
c        IN=IN+1
c        FITEMP=FI0(IN)
c        WRITE(12,REC=IND) FITEMP
c        GO TO 20
c      ENDIF
c      IF(I.EQ.NI) THEN
c        IM=IM+1
c        FITEMP=FI1(IM)
c        WRITE(12,REC=IND) FITEMP
c      GO TO 20
c      ENDIF
cc       WRITE(12,REC=IND) FI(IND-NJ-7)
c      WRITE(12,REC=IND) FI(IND-NJ-LAST)
c20      CONTINUE
c10      CONTINUE
      CLOSE(12)
      
      RETURN
      END
