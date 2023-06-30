!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!! Author : Mohammad Razeghi          9665611105      DATE : 7/1/2018
!!!!CFD PROJECT#3...SOLVING THERMAL BOUNDRAY LAYER EQUATION
!!!!    Wall Heat Flux is Constant
!!!!	Upwind Method
!         -----------------------------
!         |							  |  H = 0.001 m
!         |							  |  fully developed
!         ----------------------------- 
!                   L = 1 m
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   VARIABLES   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!   L  = Length of the Channel
!   H  = Height of the Channel
!   W  = Width of the Channel
!   IM = Number of Grids in x_dir
!   JM = Number of Grids in y_dir
!   dx = Size of the x_dir Grid
!   dy = Size of the y_dir Grid
!   T  = Temperature
!   Um = Mean Velocity
!   Tw = Wall Temp
!   q  = Wall Flux
!   k  = Conduction Coefficient
!   ro = Density
!  alfa= Thermal Diffusion
!  RE  = REYNOLDS NUMBER
!  VI  = VISCOUSITY

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 PROGRAM UPWIND1
    IMPLICIT NONE

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  DEFINING VARIABLES  !!!!!!!!!!!!!!!!!!!!!!!!!!!!

    INTEGER,PARAMETER                  ::  IM=80,JM=80
    Real*8,PARAMETER                   ::  T_inf=300.0,RE=1000.0,PR=7.,ERROR=1E-04,Tw=450.
    Real*8,PARAMETER                   ::  alfa=1.436E-07,ro=1E+03,k=0.6,W=1E-03
    INTEGER                            ::  i,j,ITER
    REAL*8                             ::  L,H=1E-03,Um,dx,dy,err=1.,GAMA,VI
    REAL*8,Dimension(2:IM-1,2:JM-1)    ::  Fn,Fs,Fe,Fw,De,Dw,Dn,Ds,sp,su
    REAL*8,Dimension(IM,JM)            ::  T,T2,Xp,Yp,U,V,X,Y
    REAL*8,Dimension(IM)               ::  T_STAR,TM,X_STAR,HH,NU
    REAL*8,Dimension(2:IM-1,2:JM-1)    ::  an,as,ae,aw,ap,Areaw,Areas,Areae,Arean
    REAL*8,Dimension(2:JM-1)           ::  AA1,BB1,CC1,DD1
    REAL*8,Dimension(2:IM-1)           ::  AA,BB,CC,DD

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!  Finding The Average Velocity   !!!!!!!!!!!!!!!!!!!!!
	GAMA	= ALFA*RO
    VI		= PR*ALFA
	UM		= (RE*VI)/H
    !L		= H*0.05*RE*PR
    L=1
    

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   CALCULATING GRID SIZE !!!!!!!!!!!!!!!!!!!!!!!!!!!
    dx  =   (L/(IM-2))
    dy  =   (H/(JM-2))
    !!!!!!!!!!!!!!!!!!!!!!!!   DEFINING X,Y,U VECTORS FOR GRID    !!!!!!!!!!!!!!!!!!!!
    X(IM,:) = L
    Y(:,JM) = H/2.
    DO i=1,IM-1
        DO j=1,JM-1
            X(i,j)    =   (i-1) * dx
            X(i,jm)   =   (i-1) * dx
            Y(i,j)    =   -H/2. + (j-1) * dy
            Y(im,j)   =   -H/2. + (j-1) * dy
        END DO
    END DO

    Do i=2,im-1
        Do j=2,jm-1
            Yp(1,j) =(Y(1,j-1)+Y(1,j))/2.0
            Yp(im,j)=(Y(im,j-1)+Y(im,j))/2.0
            Xp(i,1) =(X(i-1,1)+X(i,1))/2.0
            Xp(i,jm)=(X(i-1,jm)+X(i,jm))/2.0
            Xp(i,j) =(X(i-1,j)+X(i,j))/2.0
            Yp(i,j) =(Y(i,j-1)+Y(i,j))/2.0
         END DO
     END DO
    Xp(1,:)=X(1,:)
    Yp(:,1)=Y(:,1)
    Xp(im,:)=X(im,:)
    Yp(:,jm)=Y(:,jm)
    DO i=1,im
        DO j=1,JM
            U(i,j)    =   3.0/2.0*um*(1-((2*Yp(i,j))/H)**2)
            V(i,j)    =   0.0
        END DO
    END DO
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!  initial Values  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    T(:,:) =   300
    T(1,2:jm-1)  =   T_inf

		T(:,1)=Tw
        T(:,jm)=Tw

    T2 = T


    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !   Defining D,F for each surface when grid is structured

    DO i=2,im-1
        DO j=2,jm-1
            Areae(i,j)= W*(Y(i,j)-Y(i,j-1))
            Areaw(i,j)= W*(Y(i-1,j)-Y(i-1,j-1))
            Arean(i,j)= W*(X(i,j)-X(i-1,j))
            Areas(i,j)= W*(X(i,j-1)-X(i-1,j-1))
         END DO
    END DO
    Do i=2,Im-1
        Do j=2,jm-1

          Fe(i,j)   =   ro*(U(i+1,j)+U(i,j))/2.0*Areae(i,j)
          Fw(i,j)   =   ro*(U(i,j)+U(i-1,j))/2.0*Areaw(i,j)
          Fn(i,j)   =   ro*(V(i,j+1)+V(i,j))/2.0*Arean(i,j)
          Fs(i,j)   =   ro*(V(i,j)+V(i,j-1))/2.0*Areas(i,j)

          De(i,j)   =   (GAMA*Areae(i,j))/(Xp(i+1,j)-Xp(i,j))
          Dw(i,j)   =   (GAMA*Areaw(i,j))/(Xp(i,j)-Xp(i-1,j))
          Dn(i,j)   =   (GAMA*Arean(i,j))/(Yp(i,j+1)-Yp(i,j))
          Ds(i,j)   =   (GAMA*Areas(i,j))/(Yp(i,j)-Yp(i,j-1))

        end do
    end do

    !!!!!!!!!!  Upwind Method   !!!!!!!

    Do i=2,im-1
        Do j=2,jm-1

            ae(i,j) =   De(i,j)     +   max(.0,-Fe(i,j))
            aw(i,j) =   Dw(i,j)     +   max(.0,Fw(i,j))
            an(i,j) =   Dn(i,j)     +   max(.0,-Fn(i,j))
            as(i,j) =   Ds(i,j)     +   max(.0,Fs(i,j))

            Su(i,j) =   0.0
            Sp(i,j) =   0.0
            
            END DO
            END DO

            
     ! Boundary Condition
     Do i=2,im-1
        Do j=2,jm-1

            !!! Left Wall
            Su(2,j) =   aw(2,j)*T_inf
            Sp(2,j) =   -aw(2,j)

            !!! Right Wall
            !Su(im-1,j)  =   0.0
            !Sp(im-1,j)  =   -ae(im-1,j)

            !!! Top Wall   Constant Wall Temperature
            Su(i,jm-1)  =   an(i,jm-1)*Tw
            Sp(i,jm-1)  =   -an(i,jm-1)

            !!! Bottom Wall  Constant Wall Temperature
            Su(i,2) =   as(i,2)*Tw
            Sp(i,2) =   -as(i,2)
		END DO
     END DO

    !!!    Boundary Condition at Corners

    !!!    Left Bottom Corner
    Sp(2,2)    =   -aw(2,2)-as(2,2)
    Su(2,2)    =   aw(2,2)*T_inf + as(2,2)*Tw

    !!!    Left Top Corner
    Sp(2,jm-1) =   -aw(2,jm-1)-an(2,jm-1)
    Su(2,jm-1) =   aw(2,jm-1)*T_inf + an(2,jm-1)*Tw

    !!!    Right Bottom Corner
    Sp(im-1,2) =   -as(im-1,2)
    Su(im-1,2) =   as(im-1,2)*Tw

    !!!    Right Top Corner
    Sp(im-1,jm-1)= -an(im-1,jm-1)
    Su(im-1,jm-1)= an(im-1,jm-1)*Tw

    ae(im-1,2:jm-1)  =   0.0
    an(2:im-1,jm-1)  =   0.0
    aw(2,2:jm-1)     =   0.0
    as(2:im-1,2)     =   0.0


    do i=2,im-1
        do j=2,jm-1
            aP(i,j)  =   as(i,j)+an(i,j)+aw(i,j)+ae(i,j)+Fe(i,j)-Fw(i,j)+Fn(i,j)-Fs(i,j)-Sp(i,j)
        end do
    end do
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   Using ADI Method  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	DO WHILE (ERR>=ERROR)
        T2(:,:)=T(:,:)
        

        ! SWEEPING WEST TO EAST
        DO I=2,IM-1
            DO J=2,JM-1
                DD1(J) = AP(I,J)
                BB1(J) = AS(I,J)
                AA1(J) = AN(I,J)
                CC1(J) = AW(I,J)*T(I-1,J)+AE(I,J)*T(I+1,J)+SU(I,J)
            END DO
            CALL TDMA(2,JM-1,BB1,DD1,AA1,CC1)
            T(I,2:JM-1)=CC1(2:JM-1)
        END DO
        
         !SWEEPING SOUTH TO NORTH
         DO J=2,JM-1
           DO I=2,IM-1
             DD(I) = AP(I,J)
             BB(I) = AW(I,J)
             AA(I) = AE(I,J)
             CC(I) = AN(I,J)*T(I,J+1)+AS(I,J)*T(I,J-1)+SU(I,J)
           END DO
           CALL TDMA(2,IM-1,BB,DD,AA,CC)
           T(2:IM-1,J) = CC(2:IM-1)
         END DO         
         ERR=SUM(ABS(T2(2:IM-1,2:JM-1)-T(2:IM-1,2:JM-1)))
         
    END DO


	DO I=1,IM
        TM(I) = (1./(UM*H))*(sum(u(i,2:jm-1)*T(i,2:jm-1)*DY))
        T_STAR(I)=(T(I,JM/2)-T(i,1))/(TM(I)-T(i,1))
    END DO


    !! NUSSELT NUMB

    DO I=1,IM
      HH(I)=(2*K*(T(I,JM)-T(I,JM-1)))/(DY*(T(i,jm)-TM(I)))
      NU(I)=(HH(I)*2*H)/K
      X_STAR(I)= (8*XP(I,JM/2))/(H*re*pr)
      
    END DO
    
  !!!!!!!!!!!!!!!!!!!!!!  Getting Outputs  !!!!!!!!!!!!!!  	

	print*,L
    print*,dx
    print*,dy
    print*,err
    
!!!!!!!!!!!!!!!!!	Opening Files    !!!!!!!!!!!!!!!!!!!!!!

    OPEN (UNIT = 11, FILE="Find T.plt", STATUS="REPLACE")
    OPEN (UNIT = 12, FILE="T in Center of Channel.plt", STATUS="REPLACE")
    OPEN (UNIT = 13, FILE="T DIMMENSSION.plt",STATUS="REPLACE")
	OPEN (UNIT = 14, FILE="NUSSELT.plt",STATUS="REPLACE")


	WRITE(11,*)'VARIABLES = X,Y,T_N'
    WRITE(11,*)'zone     i=',im-2,'         j=',jm-2,'                      F=POINT'

	Do I=2,IM-1
    	 WRITE(12,*) Xp(i,jm/2),T(i,jm/2)
         WRITE(13,*) XP(I,JM/2),T_STAR(I)
         WRITE(14,*) X_star(i),NU(i),hh(i)
    END DO
	
    DO I=2,IM-1
        DO j=2,Jm-1
            WRITE(11,*) Xp(i,j),Yp(i,j),T(i,j)
        END DO
    END DO

!!!!!!!!!!!!!!   Closing Files  !!!!!!!!!!!!!!
	CLOSE (UNIT = 11)
    CLOSE (UNIT = 12)
    CLOSE (UNIT = 13)
    CLOSE (UNIT = 14)

END PROGRAM UPWIND1

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! SUBROUTINE  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE TDMA(IL,IU,BB,DD,AA,CC)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: IL, IU
    REAL*8, DIMENSION(IL:IU), INTENT(IN) :: AA, BB
    REAL*8, DIMENSION(IL:IU), INTENT(INOUT) :: CC, DD
    !..
    !.. IL = SUBSCRIPT OF FIRST EQUATION
    !.. IU = SUBSCRIPT OF LAST EQUATION
    !.. BB = COEFFICIENT BEHIND DIAGONAL
    !.. DD = COEFFICIENT ON DIAGONAL
    !.. AA = COEFFICIENT AHEAD OF DIAGONAL
    !.. CC = ELEMENT OF CONSTANT VECTOR
    !..
    !.. ESTABLISH UPPER TRIANGULAR MATRIX
    !..
    INTEGER :: LP,i,k
    REAL*8, DIMENSION(IL:IU) :: R,CP
    R(IL)  =  AA(IL) / DD(IL)
    CP(IL) =  CC(IL) / DD(IL)

    LP = IL + 1
    DO i = LP, IU
        R(I)  =  AA(I) / (DD(I)-BB(I)*R(I-1))
        CP(I) =  (BB(I)*CP(I-1)+CC(I))/(DD(I)-BB(I)*R(I-1))
    END DO

    !..
    !.. BACK SUBSTITUTION
    !..

        CC(IU) = CP(IU)
    DO i = LP, IU
        k = IU - i + IL
        CC(k) = CC(K+1)*R(K)+CP(K)
    END DO
    !..
    !.. SOLUTION STORED IN CC
    !..

    RETURN
END SUBROUTINE TDMA
