program   A_PSOR
    implicit none
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!                                                                                         !
! AUTHOR: Mohammad Jamal Razeghi         Student Number : 9665611105
! CFD Project #2                        Date: May 2018

! Program For Solving: Lid Driven Cavity Using Vorticity-Stream function Equation
    !Using PSOR method to discretize the Poisson equation of Stream Lines
    ! Our final solution is for steady state of vorticity
    ! The error of this result is 0.001

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                            Defining Variables                                          !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!      U , V      : Velocity in X , Y Direction
!        P        : Pressure
!     Px & Py     : Pressure Derivation   in x & y direction
!     sy,sye      : Stream function in the solution
!        W        : Vorticity
!        m        : A matrix for calculation of Vorticity
!       im        : Number Of Nodes in X Direction
!       jm        : Number Of Nodes in Y Direction
!  	  dx , dy     : Step Size Of Space
!  	  a , b       : length and width
!   	dt		  : Step Size Of Time
!        Re       : Reynolds Number
!        Nu       : Kinematic Viscosity
!    err , error  : Error in the solution
!       beta      : Coefficient over Relaxation
!       res       : Residual for calculation of vorticity in steady form
!     numstep     : Number of time step for solution
!     it1         : Number of iterations for streamFunction convergence
!     it2         : Number of iterations for vorticity convergence

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                         Variables and Initial Value                                     !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    real,Parameter                   ::  a=0.0,b=1.0,vel=1.0,Re=100.0,dx=.01,dy=.01
    real,Dimension(:,:),Allocatable  ::  w,vor,m,sy,sye,u,v,p,x,y,px,py
    real                             ::  beta=1.5,error=.001,t=.0,dt,Nu,err=.0,d=.0,res=1.0
    integer                          ::  i=0,j=0,k=0,it1=0,it2=0,im,jm,numstep
    Open (UNIT = 1, FILE="ComparingPSOR&ADI.txt", STATUS="REPLACE")
    WRITE(1,*)'VARIABLES = Voriteration,VorticityRes,Streamite1,Error'
    WRITE(1,*)'zone'
    Nu=1/Re
    !dt=(0.25*(dx**2))/Nu
    dt=.001
    im=(b-a)/(dx)+1
    jm=(b-a)/(dy)+1
    numstep=5/dt+1
    write(*,*) im
    write(*,*) jm
    write(*,*) dt
    write(*,*) numstep
    Allocate (w(im,jm),vor(im,jm),m(im,jm),sy(im,jm),sye(im,jm),u(im,jm),v(im,jm))
    Allocate (p(im,jm),x(im,jm),y(im,jm),px(im,jm),py(im,jm))

!   Initial Values for Velocity Stream Function and Vorticity
!   sy is constant on the walls. We can consider that equal to zero(0)

    u(1:im,1:jm)=.0
    v(1:im,1:jm)=.0
    u(2:im-1,jm)=vel
    sy(1:im,1:jm)=.0
    w(1:im,1:jm)=.0
    m(1:im,1:jm)=.0

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                 Start the Time to Solve both Equation                           !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    main :do k=1,numstep
        t=k*dt
        err=1
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                               Solving the Stream Function Equation with PSOR                         !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        do while (err>=error)
            it1=it1+1
            sye=sy
            vor=w
            do j=2,jm-1
                do i=2,im-1
                    sy(i,j)=(1.0-beta)*sy(i,j)+(0.25*beta)*(sy(i+1,j)+sy(i-1,j)+sy(i,j+1)+sy(i,j-1)+(dx*dx)*w(i,j))
                end do
            end do
            err=0.0
            do i=1,im
                do j=1,jm
                    err=err+abs(sy(i,j)-sye(i,j))
                end do
            end do
		end do

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!							Solving Vorticity Using First Order Implicit Euler						      !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!	Vorticity in Boundary condition

        do i=2,im-1
            do j=2,jm-1
                w(i,1)=-2.0*sy(i,2)/(dx*dx)					! bottom wall
                w(i,jm)=-2.0*sy(i,jm-1)/(dx*dx)-2.0/dx		! top wall
                w(1,j)=-2.0*sy(2,j)/(dx*dx)					! right wall
                w(im,j)=-2.0*sy(im-1,j)/(dx*dx)				! left wall
            end do
		end do

!	Vorticity in Internal Nodes

        do i=2,im-1
            do j=2,jm-1
                m(i,j)=(1/(4*dx*dx))*((-sy(i,j+1)+sy(i,j-1))*(w(i+1,j)-w(i-1,j)) &
                    +(sy(i+1,j)-sy(i-1,j))*(w(i,j+1)-w(i,j-1))  &
                        +(4/Re)*(w(i+1,j)+w(i-1,j)+w(i,j+1)+w(i,j-1)-4.0*w(i,j)))
            end do
        end do
        w(2:im-1,2:jm-1)=w(2:im-1,2:jm-1)+dt*m(2:im-1,2:jm-1)
        res=0.0
        do i=1,im
            do j=1,jm
            res=res+abs(vor(i,j)-w(i,j))
            end do
        end do
        it2=it2+1
        write(1,*)    it2,res,it1,0.001

! Exit the calculation for steady state of vorticity

        if (res<=error) exit main
	end do main
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!										Finding U & V       										!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    do j=2,jm-1
        do i=2,im-1
            u(i,j)=((sy(i,j+1)-sy(i,j-1))/(2*dx))
            v(i,j)=((-sy(i+1,j)+sy(i-1,j))/(2*dx))
        end do
    end do
    do j=1,jm
        do i=1,im
            x(i,j)=dx*(i-1)
            y(i,j)=dx*(j-1)
        end do
	end do

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!									Finding Pressure												  !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    p(1:im,1:jm)=0.0
    px=p
    py=p
    do i=2,im-1
        do j=2,jm-1
            px(i,j)=(-u(i,j)*((u(i+1,j)-u(i-1,j))/(2*dx))-v(i,j)*((u(i,j+1)-u(i,j-1))/(2*dx)) &
                +(1/(Re*dx*dx))*(u(i+1,j)+u(i-1,j)+u(i,j-1)+u(i,j+1)-4*u(i,j)))
        end do
    end do
    do j=1,jm
        py(1,j)=(1/(Re*dx*dx))*(-v(4,j)+4*v(3,j)-5*v(2,j)+2*v(1,j))
    end do
    do i=1,im
        px(i,1)=(1/(Re*dx*dx))*(-u(i,4)+4*u(i,3)-5*u(i,2)+2*u(i,1))
    end do
    do i=2,im-1
        px(i,jm)=(-u(i,jm)*(u(i+1,jm)-u(i-1,jm))/(2*dx)) &
            +(1/(Re*dx*dx))*((u(i+1,jm)+u(i-1,jm)-2*u(i,jm))+2*u(i,jm)-5*u(i,jm-1)+4*u(i,jm-2)-u(i,jm-3))
    end do
    px(1,jm)=(1/(Re*dx*dx))*(2*u(1,jm)-5*u(2,jm)+4*u(3,jm)-u(4,jm))

!	Reference pressure

    p(1,1)=.0
    do j=1,jm-1
        p(1,j+1)=p(1,j)+py(1,j)*dx
    end do
    do j=1,jm
        do i=1,im-1
            p(i+1,j)=p(i,j)+px(i,j)*dx
        end do
	end do

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!									  Writing Tecplot Outputs									     !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    close (UNIT = 1)
    open (UNIT = 2, FILE="StreamLineCOUNTOUR.txt", STATUS="REPLACE")
    WRITE(2,*)'VARIABLES = X,Y,U,V,sy,w,p'
    WRITE(2,*)'ZONE     i=',im,'         j=',jm,'                      F=POINT'
    do j=1,jm
        do i=1,im
            write(2,*)    x(i,j),y(i,j),u(i,j),v(i,j),sy(i,j),w(i,j),p(i,j)
        end do
    end do
    close (UNIT = 2)
    open (UNIT = 3, FILE="U_Velocity.txt", STATUS="REPLACE")
    do j=1,jm
        write(3,*)    u(1+(im-1)/2,j),y(1+(im-1)/2,j)
    end do
    close (UNIT = 3)
    open (UNIT = 4, FILE="V_Velocity.txt", STATUS="REPLACE")
    do i=1,im
        write(4,*)    x(i,1+(jm-1)/2),v(i,1+(jm-1)/2)
    end do
    close (UNIT = 4)

END program A_PSOR
