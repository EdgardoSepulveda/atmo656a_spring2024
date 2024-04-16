	REAL RRR(200)  
        REAL L, MU0
	common /PI/g,mu0,nnn
	g=0.85
	PI = 3.1415926
	X0 = 0.0
	Y0 = 0.0
	ZTOP = 4.0
	ZBASE = 0.0
	Z0 = ZTOP
        mm=9
c	do 1000 mm=0, 10
	w0=mm/10.
	NUMBER = 10000
	NNN = 200
c       read(*,*) number
	MU0=-0.7
	THETA0 = ACOS( MU0 )
	FAI0 = 3.0*PI / 2.0

	SECIN = 0.0
	ASEC = SECNDS(SECIN)
	II = - INT(ASEC)
	RI = RAN1( II )
c	WRITE(*,*) II, RI
        II = 0

	call scatang(NNN,RRR)
	
	NR = 0
	NT = 0
        NA=0
	DO 200 N = 1, NUMBER
	THETA = THETA0
	FAI = FAI0
	RI = RAN1( II )
c	WRITE(*,*) THETA*180/PI, FAI*180/PI
c	WRITE(*,*) N,RI

	L = - 1.0 * ALOG( 1 - RI )

	X = X0 + L * SIN( THETA ) * COS( FAI )
	Y = Y0 + L * SIN( THETA ) * SIN( FAI )
	Z = Z0 + L * COS( THETA )
c	WRITE(*,*) N, L, X, Y, Z
	kkk=0
  100   RI=RAN1(II)	
	If(w0 .GT. RI) THEN
            NA=NA+1
	    GOTO 200
        ENDIF
  
	IF( Z .GT. ZTOP ) THEN
	    NR = NR +1
	    GOTO 200
	ENDIF

	IF( Z .LT. ZBASE ) THEN
	    NT = NT +1
	    GOTO 200
	ENDIF

	RI = RAN1( II )
	ll=int(RI*200.)+1
	thetas=RRR(ll)
	RI = RAN1( II )
	FAIs = 2*PI*RI
	RI = RAN1( II )
	L = - 1.0 * ALOG( 1 - RI )

	DX = L*(SIN(THETAs)*SIN(FAIs)*COS(THETA)*COS(FAI)
     &     +    SIN(THETAs)*COS(FAIs)*SIN(FAI)
     &	   +    COS(THETAs)*SIN(THETA)*COS(FAI) )
	DY = L*(SIN(THETAs)*SIN(FAIs)*COS(THETA)*SIN(FAI)
     &     -    SIN(THETAs)*COS(FAIs)*COS(FAI)
     &	   +    COS(THETAs)*SIN(THETA)*SIN(FAI) )
	DZ = L*(COS(THETAs)*COS(THETA)
     &     -    SIN(THETAs)*SIN(FAIs)*SIN(THETA))

	X = X + DX
	Y = Y + DY
	Z = Z + DZ

c	write(*,*)  THETA*180/PI, THETAs*180/PI
 	THETA = ACOS(DZ/L)
c	write(*,*)  THETA*180/PI, DZ, L, '--------'
	IF( DX .EQ. 0.0 ) THEN
	   IF( DY .GE. 0.0 ) THEN
	      FAI = PI / 2.0
	   ELSE 
	      FAI = 3.0 * PI / 2.0
           ENDIF
       	ELSE
	    FAI = ATAN( DY / DX )
	    IF( DX .LT. 0.0 ) THEN
	       FAI = FAI + PI
	    ELSE IF( DY .LT. 0.0 ) THEN
	       FAI = FAI + 2.0*PI
   	    ENDIF
	ENDIF	
c	WRITE(*,*) THETA*180/PI,FAI*180/PI,THETAs*180/PI,FAIs*180/PI
c	WRITE(*,*) X, Y, Z
	kkk=kkk+1
	GOTO 100
 200    CONTINUE
	T = FLOAT(NT) / NUMBER
	R = FLOAT(NR) / NUMBER
	A = FLOAT(NA) / number
c	WRITE(*,*) 'NR=', NR, ' NT=', NT
	WRITE(*,*)  w0,   A,   R,   T
c 1000   continue
	STOP
	END

	FUNCTION RAN1(IDUM)
	DIMENSION R(97)
	PARAMETER  (M1=259200,IA1=7141,IC1=54773,RM1=1.0/M1)
	PARAMETER  (M2=134456,IA2=8121,IC2=28411,RM2=1.0/M2)
	PARAMETER  (M3=243000,IA3=4561,IC3=51349)
	DATA IFF/0/

	IF( IDUM.LT.0.OR.IFF.EQ.O ) THEN
	IFF = 1
	IX1 = MOD( IC1-IDUM, M1 )
	IX1 = MOD( IA1*IX1+IC1, M1 )
	IX2 = MOD( IX1, M2 )
	IX1 = MOD( IA1*IX1+IC1, M1 )
	IX3 = MOD( IX1, M3 )
	DO 10 J = 1, 97
	   IX1 = MOD( IA1*IX1+IC1, M1 )
	   IX2 = MOD( IA2*IX2+IC2, M2 )
	   R(J) = (FLOAT(IX1)+FLOAT(IX2)*RM2)*RM1
 10     CONTINUE
	  IDUM = 1
	ENDIF
C
	IX1 = MOD( IA1*IX1+IC1, M1 )
	IX2 = MOD( IA2*IX2+IC2, M2 )
	IX3 = MOD( IA3*IX3+IC3, M3 )
	J = 1 + (97*IX3)/M3
	IF( J .GT. 97 .OR. J .LT. 1 ) PAUSE
	RAN1 = R(J)
	R(J) = (FLOAT(IX1)+FLOAT(IX2)*RM2) * RM1
	RETURN
	END

	subroutine scatang(n,theta)
	REAL  theta(n), mu(201), amu(200)
	mu(1)=-1
	b=1.+g**2
	c=2.*g/(nnn*(1.-g**2))
	d=b-2.*g*mu0
	do 90 i=2, n+1
	mu(i)=(b-(1./(c+1./SQRT (d)))**2.)/(2.*g)
        amu(i-1)=(mu(i)+mu(i-1))/2.
	THETA(i-1)=acos(amu(i-1))
c 	write(*,*) i-1 , amu(i-1), theta(i-1)
 90     continue
	return
	end

