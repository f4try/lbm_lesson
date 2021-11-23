!========================================================================
!    This code was written by Li Chen at Xi'an Jiaotong University.          
!    This code was for single-phase in a 2D channel. 
!    Pressures are known at the left inlet and right outlet,
!    and non_slip conditions at the top and bottom walls.
!    One can refer to the following papers for more details: 
!    Li Chen et al., Water Resources Research Volume: 50(12): 9343-9365, 2014
!    Li Chen: lichennht08@mail.xjtu.edu.cn. 
!========================================================================
MODULE START_L
    PARAMETER (nx=21,ny=81)   
    integer::I,J,K,LAST,ITER
    double precision,PARAMETER::XL=20.E-5,YL=80.E-5     
    double precision,dimension(nx)::X
    double precision,dimension(ny)::Y
    double precision,dimension(0:nx+1,0:ny+1)::U,V,PRE
!----------------------------------------------------
    double precision::DX,DY,DT
    double precision,PARAMETER::C=1.d0,CS2=1.d0/3.d0
    integer,dimension(0:8)::FCX=(/0,1,0,-1,0,1,-1,-1,1/)
    integer,dimension(0:8)::FCY=(/0,0,1,0,-1,1,1,-1,-1/)
    double precision,dimension(0:8)::wi
    double precision,dimension(0:2)::lambda
    !------------------------------------------------------------------ 
    double precision,dimension(0:8,0:nx+1,0:ny+1)::f1,ff1
    double precision::ftao,vmu_phy,vmu_lat,preleft,preright,feq1
!--------------------------------------------------------------------
    integer,dimension(0:nx+1,0:ny+1)::ls
    logical,dimension(0:nx+1,0:ny+1)::walls
    double precision::delta
    double precision::sumc_last,sumu
 END MODULE
!========================================================================

!========================================================================
PROGRAM MAIN
USE START_L

    CALL SOLID_STRUCTURE
    CALL INITIALIZATION

    DO iter=1,last
        CALL COLLISIONF
        CALL STREAMF
        CALL BOUNDARYF
        CALL MACROF
        if(mod(iter,1000).eq.0) CALL OUTPUT
    ENDDO  
     
END PROGRAM
!==========================================================================

!==========================================================================
SUBROUTINE INITIALIZATION
USE START_L
double precision::z1,z2

    dx=xl/float(nx-1)
    dy=dx  
    last=500000
    delta=1.d0
    lambda(0)=-5.d0/3.d0
    lambda(1)=1.0d0/3.d0
    lambda(2)=1.d0/12.d0
    wi(0)=4.d0/9.d0
    wi(1:4)=1.d0/9.d0
    wi(5:8)=1.d0/36.d0
    vmu_phy=20.e-6
    ftao=1.d0
    vmu_lat=(ftao-0.5d0)/3.d0
    scale=vmu_phy/vmu_lat
    dt=dx**2./scale
    preleft=1.0002d0
    preright=1.d0
    
    do j=1,ny
    do i=1,nx
        pre(i,j)=preleft-float(i-1)/float(nx-1)*(preleft-preright)
        u(i,j)=0.d0
        v(i,j)=0.d0
    enddo
    enddo
	 
    do j=1,ny
    do i=1,nx	
        z2=u(i,j)**2.d0+v(i,j)**2.d0
        do k=0,8 
            z1=fcx(k)*u(i,j)+fcy(k)*v(i,j)
            if(k.eq.0) then
                feq1=lambda(0)*pre(i,j)+wi(k)*(3.d0*z1+4.5d0*z1**2.-1.5d0*z2)
            elseif(k.le.4.and.k.ge.1) then
                feq1=lambda(1)*pre(i,j)+wi(k)*(3.d0*z1+4.5d0*z1**2.-1.5d0*z2)
            elseif(k.le.8.and.k.ge.5) then
                feq1=lambda(2)*pre(i,j)+wi(k)*(3.d0*z1+4.5d0*z1**2.-1.5d0*z2)
            endif				 
            f1(k,i,j)=feq1
            ff1(k,i,j)=feq1
        enddo 
    enddo
    enddo
    ! write(*,*) f1(3,2,1),ff1(3,2,1),feq1

RETURN
END SUBROUTINE
!==========================================================================

!========================================================================-=
SUBROUTINE SOLID_STRUCTURE
USE START_L
    ! ls represents the porous structure: 0 denotes nodes of void space, 1 denotes solid node.
    ls=0
    ls(:,ny:ny+1)=1
    ls(:,0:1)=1    
    walls=.false. 
    do j=0,ny+1
    do i=0,nx+1    
        if(ls(i,j).eq.1) then
            walls(i,j)=.true.
        endif
    enddo
    enddo    
RETURN
END SUBROUTINE
!==========================================================================

!==========================================================================
SUBROUTINE COLLISIONF
USE START_L
double precision::z1,z2

    do j=1,ny
    do i=1,nx
        if(.not.walls(i,j)) then
        z2=u(i,j)**2.+v(i,j)**2.
        do k=0,8
            z1=fcx(k)*u(i,j)+fcy(k)*v(i,j)
            if(k.eq.0) then
                feq1=lambda(0)*pre(i,j)+wi(k)*(3.d0*z1+4.5d0*z1**2.-1.5d0*z2)
            elseif(k.le.4.and.k.ge.1) then
                feq1=lambda(1)*pre(i,j)+wi(k)*(3.d0*z1+4.5d0*z1**2.-1.5d0*z2)
            elseif(k.le.8.and.k.ge.5) then
                feq1=lambda(2)*pre(i,j)+wi(k)*(3.d0*z1+4.5d0*z1**2.-1.5d0*z2)
            endif					  
            ff1(k,i,j)=f1(k,i,j)-1.d0/ftao*(f1(k,i,j)-feq1)
	   enddo
     endif
    enddo
    enddo 

RETURN
END SUBROUTINE
!==========================================================================

!==========================================================================
SUBROUTINE STREAMF
USE START_L

!-------------------periodic boundary along y-------------------------------------
    do j=1,ny	    
    do i=1,nx
    do k=0,8 
	     f1(k,i,j)=ff1(k,i-int(fcx(k)),j-int(fcy(k)))
    enddo
    enddo
    enddo

RETURN
END SUBROUTINE
!==========================================================================

!==========================================================================
SUBROUTINE BOUNDARYF
USE START_L

    do j=1,ny
    do i=1,nx       
	    if(walls(i,j)) then
		  ff1(1,i,j)=f1(3,i,j)
		  ff1(3,i,j)=f1(1,i,j)
		  ff1(2,i,j)=f1(4,i,j)
		  ff1(4,i,j)=f1(2,i,j)
		  ff1(5,i,j)=f1(7,i,j)
		  ff1(7,i,j)=f1(5,i,j)
		  ff1(6,i,j)=f1(8,i,j)
		  ff1(8,i,j)=f1(6,i,j)
	   endif	   
    enddo
    enddo
	
    do j=1,ny
        z1=fcx(1)*u(2,j)+fcy(1)*v(2,j)
		z2=u(2,j)**2.d0+v(2,j)**2.d0
        feq1=lambda(1)*pre(2,j)+wi(1)*(3.d0*z1+4.5d0*z1**2.-1.5d0*z2)
		f1(1,1,j)=lambda(1)*pre(1,j)+wi(1)*(3.d0*z1+4.5d0*z1**2.-1.5d0*z2) +F1(1,2,J)-FEQ1
 		z1=fcx(5)*u(2,j)+fcy(5)*v(2,j)
        feq1=lambda(2)*pre(2,j)+wi(5)*(3.d0*z1+4.5d0*z1**2.-1.5d0*z2)
		f1(5,1,j)=lambda(2)*pre(1,j)+wi(5)*(3.d0*z1+4.5d0*z1**2.-1.5d0*z2)+F1(5,2,J)-FEQ1
 		z1=fcx(8)*u(2,j)+fcy(8)*v(2,j)
        feq1=lambda(2)*pre(2,j)+wi(8)*(3.d0*z1+4.5d0*z1**2.-1.5d0*z2)
		f1(8,1,J)=lambda(2)*pre(1,j)+wi(8)*(3.d0*z1+4.5d0*z1**2.-1.5d0*z2)+F1(8,2,J)-FEQ1
  enddo
        
  do j=1,ny
		z1=fcx(3)*u(nx-1,J)+fcy(3)*v(nx-1,j)
		z2=u(nx-1,j)**2.d0+v(nx-1,j)**2.d0
        feq1=lambda(1)*pre(nx-1,j)+wi(3)*(3.d0*z1+4.5d0*z1**2.-1.5d0*z2)
		f1(3,nx,j)=lambda(1)*pre(nx,j)+wi(3)*(3.d0*z1+4.5d0*z1**2.-1.5d0*z2) +F1(3,nx-1,J)-FEQ1
		z1=fcx(6)*u(nx-1,j)+fcy(6)*v(nx-1,j)
        feq1=lambda(2)*pre(nx-1,j)+wi(6)*(3.d0*z1+4.5d0*z1**2.-1.5d0*z2)
		f1(6,nx,j)=lambda(2)*pre(nx,j)+wi(6)*(3.d0*z1+4.5d0*z1**2.-1.5d0*z2) +F1(6,nx-1,J)-FEQ1
		z1=fcx(7)*u(nx-1,j)+fcy(7)*v(nx-1,j)
        feq1=lambda(2)*pre(nx-1,j)+wi(7)*(3.d0*z1+4.5d0*z1**2.-1.5d0*z2)
		f1(7,nx,j)=lambda(2)*pre(nx,j)+wi(7)*(3.d0*z1+4.5d0*z1**2.-1.5d0*z2) +F1(7,nx-1,J)-FEQ1
   enddo
	
RETURN
END SUBROUTINE
!==========================================================================

!==========================================================================
SUBROUTINE MACROF
USE START_L

    do j=1,ny
    do i=1,nx  
        IF(.not.walls(i,j)) THEN	
            temppre=0.0d0
            tempu=0.0d0
            tempv=0.0d0
            do k=1,8
                temppre=temppre+f1(k,i,j)
                tempu=tempu+f1(k,i,j)*fcx(k)
                tempv=tempv+f1(k,i,j)*fcy(k)
            enddo
            u(i,j)=tempu
            v(i,j)=tempv
            temp1=u(i,j)**2.d0+v(i,j)**2.d0
            pre(i,j)=(temppre-2.d0/3.d0*temp1)/(-lambda(0))
            ! write(*,*) u(i,j),v(i,j),pre(i,j)
            ! write(*,*) temppre,v(i,j),f1(2,i,j)
            ! write(*,*) temppre,v(i,j),f1(2,i,j)
        elseif(walls(i,j))then
            u(i,j)=0.d0
            v(i,j)=0.d0
            pre(i,j)=0.d0
        endif             	
    enddo
    enddo

    do j=1,ny
         pre(1,j)=preleft
         pre(nx,j)=preright
    enddo

RETURN
END SUBROUTINE
!========================================================================

!========================================================================
SUBROUTINE OUTPUT
USE START_L

    sumu_last=sumu
    sumu=0.d0
    do j=1,ny
    do i=1,nx
        if(.not.walls(i,j)) sumu=sumu+u(i,j)
        ! write(*,*) u(i,j)
    enddo
    enddo             
    delta=abs(sumu_last-sumu)/abs(sumu) 
    write(*,*) iter,u(nx-10,ny/2),delta,sumu
    if(delta.le.1.e-8) then
        open(10,file="velocity_pressure.dat")
        write(10,*)'VARIABLES= X,Y,u,v,pre' 
        WRITE(10,*)'ZONE I=',nx,',J=',ny,',T=TT'
        do j=1,ny
        do i=1,nx
            write(10,*) i,j,u(i,j),v(i,j),pre(i,j)
        enddo
        enddo
        close(10)
        stop
    endif
    ! stop

RETURN
END SUBROUTINE
!========================================================================