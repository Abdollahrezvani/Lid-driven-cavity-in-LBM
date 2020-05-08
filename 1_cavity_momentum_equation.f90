


program cavity_lbm_part1
implicit none
integer, parameter :: nx=129, ny=129
real :: f(0:nx-1,0:ny-1,0:8),f_help(0:nx-1,0:ny-1,0:8) 
real :: v_x(0:8),v_y(0:8) 
real ,parameter:: delta_x=1.0, delta_t=1.0, tau=0.8
real :: visco
real :: density,wi(0:8),u_x,u_y,density_loc,p,f_eq
real :: rho, ux, uy, U_cavity, Re
integer x,y,i,t_max,time,x_e,x_w,y_n,y_s
logical Solid(0:nx-1,0:ny-1)
!-----------------
density=1.d0
visco=(tau-0.5)/3.d0
Print*, ' Input Re Number:'
read*,Re
U_cavity=visco*Re/real(nx-1)
!t_max=600000
t_max=40000
v_x(0)=0.d0; v_x(1)=1.d0 ; v_x(2)=0.d0 ; v_x(3)=-1.d0; v_x(4)=0.d0
v_y(0)=0.d0; v_y(1)=0.d0 ; v_y(2)=1.d0 ; v_y(3)=0.d0 ; v_y(4)=-1.d0
v_x(5)=1.d0; v_x(6)=-1.d0; v_x(7)=-1.d0; v_x(8)=1.d0
v_y(5)=1.d0; v_y(6)=1.d0 ; v_y(7)=-1.d0; v_y(8)=-1.d0
wi(0)=4.d0/9.d0
wi(1)=1.d0/9.d0 ;  wi(2)=1.d0/9.d0;  wi(3)=1.d0/9.d0;  wi(4)=1.d0/9.d0
wi(5)=1.d0/36.d0; wi(6)=1.d0/36.d0; wi(7)=1.d0/36.d0; wi(8)=1.d0/36.d0
!------------------ Boundary
do x=0,nx-1
	do y=0,ny-1
		solid(x,y) = .false.
	end do
end do
do x=0,nx-1
	solid(x,0) = .true.
	solid(x,ny-1) = .true.
end do
do y=0,ny-1
	solid(0,y) = .true.
	solid(nx-1,y) = .true.
end do
!---------------------------- Initial
do x=0,nx-1
	do y=0,ny-1
		do i=0,8
			f(x,y,i) = density*wi(i)
		end do
	end do
end do    
!----------------------------- Main Loop
do time = 1, t_max
     if (time .ge. 10 .and. mod(time,10*nx) .eq. 0) then      
       call check_density(nx,ny,f,time)
	   call write_results(nx,ny,solid,f,density)
	 end if
	  !------- Stream
		do x = 0, nx-1
			do y = 0, ny-1
			  y_n = mod(y+1+ny,ny)
			  x_e = mod(x+1+nx,nx)
			  y_s = mod(y-1+ny,ny)
			  x_w = mod(x-1+nx,nx) 
			  f_help(x  ,y  ,0) = f(x,y,0)
			  f_help(x_e,y  ,1) = f(x,y,1)
			  f_help(x  ,y_n,2) = f(x,y,2)
			  f_help(x_w,y  ,3) = f(x,y,3)
			  f_help(x  ,y_s,4) = f(x,y,4)
			  f_help(x_e,y_n,5) = f(x,y,5)
			  f_help(x_w,y_n,6) = f(x,y,6)
			  f_help(x_w,y_s,7) = f(x,y,7)
			  f_help(x_e,y_s,8) = f(x,y,8)
			end do
		 end do
	!---------------- Zou_He BC
	! Upper wall
	ux=U_cavity
	uy=0.d0
	y=ny-1
	do x=0,nx-1
		rho=(f_help(x,y,0)+f_help(x,y,1)+f_help(x,y,3)+2.d0*(f_help(x,y,2)+f_help(x,y,6)+f_help(x,y,5)))/(1.d0+uy)
		f_help(x,y,4)=f_help(x,y,2)-2.d0/3.d0*rho*uy
		f_help(x,y,7)=f_help(x,y,5)+0.5*(f_help(x,y,1)-f_help(x,y,3))-1.d0/6.d0*rho*uy-0.5*rho*ux
		f_help(x,y,8)=f_help(x,y,6)+0.5*(f_help(x,y,3)-f_help(x,y,1))-1.d0/6.d0*rho*uy+0.5*rho*ux
	end do
	! Lower wall
	ux=0.d0
	uy=0.d0
	y=0
	do x=0,nx-1
		rho=(f_help(x,y,0)+f_help(x,y,1)+f_help(x,y,3)+2.d0*(f_help(x,y,4)+f_help(x,y,7)+f_help(x,y,8)))/(1.d0-uy)
		f_help(x,y,2)=f_help(x,y,4)+2.d0/3.d0*rho*uy
		f_help(x,y,5)=f_help(x,y,7)-0.5*(f_help(x,y,1)-f_help(x,y,3))+1.d0/6.d0*rho*uy+0.5*rho*ux
		f_help(x,y,6)=f_help(x,y,8)-0.5*(f_help(x,y,3)-f_help(x,y,1))+1.d0/6.d0*rho*uy-0.5*rho*ux
	end do
	! Left wall
	ux=0.d0
	uy=0.d0
	x=0
	do y=0,ny-1
		rho=(f_help(x,y,0)+f_help(x,y,2)+f_help(x,y,4)+2.d0*(f_help(x,y,3)+f_help(x,y,6)+f_help(x,y,7)))/(1.d0-ux)
		f_help(x,y,1)=f_help(x,y,3)+2.d0/3.d0*rho*ux
		f_help(x,y,5)=f_help(x,y,7)-0.5*(f_help(x,y,2)-f_help(x,y,4))+1.d0/6.d0*rho*ux+0.5*rho*uy
		f_help(x,y,8)=f_help(x,y,6)+0.5*(f_help(x,y,2)-f_help(x,y,4))+1.d0/6.d0*rho*ux-0.5*rho*uy
	end do
	! Right wall
	ux=0.d0
	uy=0.d0
	x=nx-1
	do y=0,ny-1
		rho=(f_help(x,y,0)+f_help(x,y,2)+f_help(x,y,4)+2.d0*(f_help(x,y,1)+f_help(x,y,5)+f_help(x,y,8)))/(1.d0+ux)
		f_help(x,y,3)=f_help(x,y,1)-2.d0/3.d0*rho*ux
		f_help(x,y,7)=f_help(x,y,5)+0.5*(f_help(x,y,2)-f_help(x,y,4))-1.d0/6.d0*rho*ux-0.5*rho*uy
		f_help(x,y,6)=f_help(x,y,8)-0.5*(f_help(x,y,2)-f_help(x,y,4))-1.d0/6.d0*rho*ux+0.5*rho*uy
	end do
	 !--------------- Collision
		do x = 0, nx-1
			do y = 0, ny-1
				density_loc = 0.d0
				do i = 0, 8 
				  density_loc = density_loc + f_help(x,y,i)
				end do
				u_x=0.d0
				u_y=0.d0
				do i = 0, 8 
				  u_x = u_x + f_help(x,y,i) * v_x(i)
				  u_y = u_y + f_help(x,y,i) * v_y(i)
				end do			
				u_x = u_x / density_loc
				u_y = u_y / density_loc
				do i = 0, 8 
				  p = v_x(i) * u_x + v_y(i) * u_y
				  f_eq = wi(i)*density_loc*(1.d0+3*p+4.5*p*p-1.5*(u_x * u_x + u_y * u_y))
				  f(x,y,i) = f_help(x,y,i) + ( f_eq - f_help(x,y,i) ) / tau
				end do
			end do
		  end do
end do

call write_results(nx,ny,solid,f,density)
print*,'Done '
print*,'********************    end     ********************'
pause	
End program
!================== Subroutines
      subroutine check_density(nx,ny,f,time)  
      implicit none 
      integer  nx,ny,time
      real  f(0:nx-1,0:ny-1,0:8)
      integer  x,y,n
      real  n_sum
      n_sum = 0.d0
        do y = 0, ny-1
          do x = 0, nx-1
			do n = 0, 8
				n_sum = n_sum + f(x,y,n)
			end do
		  end do
		end do
      write(6,*) '*** Iteration number = ', time
      write(6,*) '*** Integral density = ', n_sum
      write(6,*) '***'
      return
      end subroutine
!==================  write
      subroutine  write_results(nx,ny,solid,f,density)
      implicit none   
      integer  nx,ny
      real ::  f(0:nx-1,0:ny-1,0:8),density
      logical  solid(0:nx-1,0:ny-1)
      integer  x,y,i,obsval
      real*8  u_x,u_y,d_loc,press,c_squ
      c_squ = 1.d0 / 3.d0
      open(11,file='Data.dat')
      write(11,*) 'VARIABLES = X, Y, Ux, Uy, PRESS, OBST' 
      write(11,*) 'ZONE I=', nx, ', J=', ny, ', F=POINT'
      do y = 0, ny-1
        do x = 0, nx-1
          if (solid(x,y)) then 
            obsval = 1
            u_x = 0.d0
            u_y = 0.d0
            press = density * c_squ
          else
            d_loc = 0.d0
            do i = 0, 8 
              d_loc = d_loc + f(x,y,i)
			end do
            u_x = (f(x,y,1) + f(x,y,5) + f(x,y,8)-(f(x,y,3) + f(x,y,6) + f(x,y,7))) / d_loc
            u_y = (f(x,y,2) + f(x,y,5) + f(x,y,6)-(f(x,y,4) + f(x,y,7) + f(x,y,8))) / d_loc
            press = d_loc * c_squ 
            obsval = 0 
          end if
          write(11,*) x, y, u_x, u_y, press, obsval
		end do
	  end do
      close(11)
      return
      pause
      end subroutine
!==================