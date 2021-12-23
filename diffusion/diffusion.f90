!========================================================================
!    This code was written by Li Chen at Xi'an Jiaotong University.          
!    This code was for pure diffusion in a 2D channel. 
!    Concentrations are known at the left inlet and right outlet,
!    and no flux boundary conditions at the top and bottom walls.
!    One can refer to the following papers for more details: 
!          Li Chen et al., International Journal of Thermal Sciences, 51:132-144,2012
!          Li Chen et al., Physical Review E, 87(4):043306,2013
!    Li Chen: lichennht08@mail.xjtu.edu.cn. 
!========================================================================
MODULE START_L
    PARAMETER (nx=61,ny=21)  
    integer*4::i,j,k,iter
    double precision,parameter::xl=60.E-6,yl=20.E-6  
!--------------------------------------------------------------------
    double precision::dx,dy,dt
    integer*4,dimension(0:4)::fcx=(/0,1,0,-1,0/)
    integer*4,dimension(0:4)::fcy=(/0,0,1,0,-1/)
!--------------------------------------------------------------------
    double precision,dimension(nx,ny)::con,u,v
    double precision,dimension(0:4,0:nx+1,0:ny+1)::g,gg
    double precision::dtao,ts,geq,diff_phy,diff_lat
    double precision::con_in,con_out
!--------------------------------------------------------------------
    logical, dimension(0:nx+1,0:ny+1)::walls
    integer*4,dimension(0:nx+1,0:ny+1)::ls
   double precision::sumc_last,sumc
    double precision::delta
!--------------------------------------------------------------
 END MODULE 
!========================================================================

!========================================================================
PROGRAM MAIN
USE START_L

    dx=xl/float(nx-1)
    dy=dx

    CALL SOLID_STRUCTURE
    CALL INITIALIZATION

    diff_phy=1.e-5 ! unit m^2/s diffusivity in physical unit
    dtao=1.d0       
    ts=0.2d0
    diff_lat=(1.d0-ts)*(dtao-0.5d0)/2.d0 !diffusivity in lattice unit
    scale=diff_phy/diff_lat
    dt=dx**2.d0/scale  ! physical time per lattice step

    last=200000
    delta=1.d0 
       	   	    
    do iter=1,last
        CALL COLLISIOND     
        CALL STREAMD               
        CALL BOUNDARYD             
        CALL MACROD
        if(mod(iter,1000).eq.0) CALL OUTPUT
    enddo	    

END PROGRAM
!========================================================================

!========================================================================
SUBROUTINE INITIALIZATION
USE START_L

    u=0.d0
    v=0.d0
    con=0.d0
    con_in=1.d0
    con_out=0.d0
    do j=1,ny
    do i=1,nx
        do k=0,4
            temp_uv=fcx(k)*u(i,j)+fcy(k)*v(i,j)
            if(k.eq.0) then
                geq=con(i,j)*(ts+0.5d0*temp_uv)
            else
                geq=con(i,j)*((1.d0-ts)/4.d0 + 0.5d0*temp_uv)         
            endif   
            gg(k,i,j)=geq
            g(k,i,j)=geq
        enddo      
    enddo
    enddo

RETURN
END SUBROUTINE
!========================================================================

!========================================================================
SUBROUTINE SOLID_STRUCTURE
USE START_L

    ! ls represents the porous structure: 0 denotes nodes of void space, 1 denotes solid node.
    ls=0   
    ls(:,0:1)=1
    ls(:,ny:ny+1)=1   
     
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
!========================================================================

!========================================================================
SUBROUTINE COLLISIOND
USE START_L

    do j=1,ny
    do i=1,nx	
        if(.not.walls(i,j)) then             
            do k=0,4
                temp_uv=fcx(k)*u(i,j)+fcy(k)*v(i,j)
                if(k.eq.0) then
                    geq=con(i,j)*(ts+0.5d0*temp_uv)
                else
                    geq=con(i,j)*((1.-ts)/4.d0 + 0.5d0*temp_uv)         
                endif   
                    gg(k,i,j)=g(k,i,j)-1.d0/dtao*(g(k,i,j)-geq)
                ! if(gg(k,i,j).gt.0.1d0) then
                !     write(*,*) i,j,k,gg(k,i,j)
                ! endif
            enddo      
        endif
    enddo
    enddo
    

RETURN
END SUBROUTINE
!========================================================================

!========================================================================
SUBROUTINE STREAMD
USE START_L
 
	do j=1,ny
	do i=1,nx
	  g(0,i,j)=gg(0,i,j)
	  g(1,i,j)=gg(1,i-1,j)
	  g(3,i,j)=gg(3,i+1,j)
	  g(2,i,j)=gg(2,i,j-1)
	  g(4,i,j)=gg(4,i,j+1)	  
	enddo
    enddo
        
RETURN 
END SUBROUTINE
!========================================================================

!========================================================================
SUBROUTINE BOUNDARYD
USE START_L

    do j=1,ny   
    do i=1,nx		 
         if(walls(i,j)) then		   
		    gg(1,i,j)=g(3,i,j)
		    gg(3,i,j)=g(1,i,j)
		    gg(2,i,j)=g(4,i,j)
		    gg(4,i,j)=g(2,i,j)
         endif		   
    enddo
    enddo

    do i=1,nx
        gg(2,i,1)=g(4,i,1)  ! bounce back for noflux boundary condition
        gg(4,i,ny)=g(2,i,ny) !bounce back for noflux boundary condition
    enddo

    do j=1,ny
        g(1,1,j)=con_in-g(0,1,j)-g(2,1,j)-g(3,1,j)-g(4,1,j) !concentration boundary condition
        g(3,nx,j)=con_out-g(0,nx,j)-g(1,nx,j)-g(2,nx,j)-g(4,nx,j) !concentration boundary condition
    enddo
                
RETURN 
END SUBROUTINE
!========================================================================

!========================================================================
SUBROUTINE MACROD
USE START_L
	
    do j=1,ny
    do i=1,nx
        if(.not.walls(i,j)) then
            con(i,j)=g(0,i,j)+g(1,i,j)+g(2,i,j)+g(3,i,j)+g(4,i,j)
        endif
    enddo
    enddo
	   
    do j=1,ny
        con(1,j)=con_in
        con(nx,j)=con_out
    enddo

RETURN
END SUBROUTINE
!========================================================================

!========================================================================
SUBROUTINE OUTPUT
USE START_L

    sumc_last=sumc
    sumc=0.d0
    do j=1,ny
    do i=1,nx
        if(.not.walls(i,j)) sumc=sumc+con(i,j)
    enddo
    enddo             
    delta=abs(sumc_last-sumc)/abs(sumc) 
    write(*,*) iter,con(nx/2,ny/2),delta,sumc
       
    open(10,file="concentration.dat")
    write(10,*)'VARIABLES= X,Y,con' 
    WRITE(10,*)'ZONE I=',nx,',J=',ny-2,',T=TT'
    do j=2,ny-1
    do i=1,nx
        write(10,*) i,j,con(i,j)
    enddo
    enddo
    close(10)

    open(10,file="concentration_x.dat")
    do i=1,nx
        write(10,*) i,con(i,ny/2)
    enddo
    close(10)

    if(delta.le.1.e-8) stop
    ! stop

RETURN
END SUBROUTINE
!========================================================================