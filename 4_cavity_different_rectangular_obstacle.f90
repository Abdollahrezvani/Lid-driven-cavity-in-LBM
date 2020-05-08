


program cavity_lbm_part4
implicit none
integer, parameter :: nx=129,ny=129
double precision f(0:nx-1,0:ny-1,0:8),f_help(0:nx-1,0:ny-1,0:8)
double precision g(0:nx-1,0:ny-1,0:8),g_help(0:nx-1,0:ny-1,0:8)
double precision v_x(0:8),v_y(0:8) 
double precision,parameter:: tau=0.8
double precision visco,Pr,th_loc,alfa,tau_g,T_initial
double precision density,wi(0:8),u_x,u_y,density_loc,t_hot,t_cold,f_eq,p,g_eq
double precision rho,ux,uy,U_cavity,Re,sum
integer x,y,i,t_max,time,x_e,x_w,y_n,y_s
logical Solid(0:nx-1,0:ny-1)
!-----------------
density=1.d0
visco=(tau-0.5)/3.d0
Print*, 'Input Re Number:'
read*,Re
Print*, 'Input Pr Number:'
read*,Pr
t_hot=1
t_cold=0
T_initial=0.0
U_cavity=visco*Re/real(nx-1)
t_max=60000
!t_max=600000
alfa=visco/Pr
tau_g=3.0d0*alfa+0.50d0
v_x(0)=0.d0; v_x(1)=1.d0 ; v_x(2)=0.d0 ; v_x(3)=-1.d0; v_x(4)=0.d0
v_y(0)=0.d0; v_y(1)=0.d0 ; v_y(2)=1.d0 ; v_y(3)=0.d0 ; v_y(4)=-1.d0
v_x(5)=1.d0; v_x(6)=-1.d0; v_x(7)=-1.d0; v_x(8)=1.d0
v_y(5)=1.d0; v_y(6)=1.d0 ; v_y(7)=-1.d0; v_y(8)=-1.d0
wi(0)=4.d0/9.d0
wi(1)=1.d0/9.d0 ;  wi(2)=1.d0/9.d0;  wi(3)=1.d0/9.d0;  wi(4)=1.d0/9.d0
wi(5)=1.d0/36.d0; wi(6)=1.d0/36.d0; wi(7)=1.d0/36.d0; wi(8)=1.d0/36.d0
!------------------ Boundary
!$$$$$$ do x=0,nx-1
!$$$$$$ do y=0,ny-1
!$$$$$$ solid(x,y)=.false.
!$$$$$$ end do
!$$$$$$ end do

do x=0,nx-1
	do y=0,ny-1
		solid(x,y) = .false.
      !solid(40:80,40:80)=.true.
      
      solid(10:30,10:30)=.true.
      solid(60:80,60:80)=.true.
      
	end do
end do

do x=0,nx-1
solid(x,0)=.true.
solid(x,ny-1)=.true.
end do
do y=0,ny-1
solid(0,y)=.true.
solid(nx-1,y)=.true.
end do
!---------------------------- Initial
do x=0,nx-1
do y=0,ny-1
do i=0,8
f(x,y,i)=density*wi(i)
g(x,y,i)=T_initial*wi(i)
end do
end do
end do    
!----------------------------- Main Loop
do time = 1, t_max
  if (time .ge. 10 .and. mod(time,10*nx) .eq. 0) then      
   call check_density(nx,ny,f,time)
   call  results(nx,ny,solid,f,density,g)
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
g_help(x ,y ,0) = g(x,y,0)
g_help(x_e,y ,1) = g(x,y,1)
g_help(x ,y_n,2) = g(x,y,2)
g_help(x_w,y ,3) = g(x,y,3)
g_help(x ,y_s,4) = g(x,y,4)
g_help(x_e,y_n,5) = g(x,y,5)
g_help(x_w,y_n,6) = g(x,y,6)
g_help(x_w,y_s,7) = g(x,y,7)
g_help(x_e,y_s,8) = g(x,y,8)
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

sum=g_help(x,y,0)+g_help(x,y,1)+g_help(x,y,2)+g_help(x,y,3)+g_help(x,y,5)+g_help(x,y,6)
g_help(x,y,4)=wi(4)/(wi(4)+wi(7)+wi(8))*(t_cold-sum)
g_help(x,y,7)=wi(7)/(wi(4)+wi(7)+wi(8))*(t_cold-sum)
g_help(x,y,8)=wi(8)/(wi(4)+wi(7)+wi(8))*(t_cold-sum)
end do

 !bounce-back
     !upper wall-barrier
     
  do x=10,30
     f_help(x,30,2)=f_help(x,30,4)
     f_help(x,30,5)=f_help(x,30,7)
     f_help(x,30,6)=f_help(x,30,8)
  end do

!upper wall barrier

do x=10,30
sum=g_help(x,30,0)+g_help(x,30,1)+g_help(x,30,4)+g_help(x,30,3)+g_help(x,30,7)+g_help(x,30,8)
g_help(x,30,2)=wi(2)/(wi(2)+wi(5)+wi(6))*(t_hot-sum)
g_help(x,30,5)=wi(5)/(wi(2)+wi(5)+wi(6))*(t_hot-sum)
g_help(x,30,6)=wi(6)/(wi(2)+wi(5)+wi(6))*(t_hot-sum)
end do

!upper wall barrier

do x=60,80
     f_help(x,80,2)=f_help(x,80,4)
     f_help(x,80,5)=f_help(x,80,7)
     f_help(x,80,6)=f_help(x,80,8)
  end do

! upper wall barrier

do x=60,80
sum=g_help(x,80,0)+g_help(x,80,1)+g_help(x,80,4)+g_help(x,80,3)+g_help(x,80,7)+g_help(x,80,8)
g_help(x,80,2)=wi(2)/(wi(2)+wi(5)+wi(6))*(t_hot-sum)
g_help(x,80,5)=wi(5)/(wi(2)+wi(5)+wi(6))*(t_hot-sum)
g_help(x,80,6)=wi(6)/(wi(2)+wi(5)+wi(6))*(t_hot-sum)
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

!$$$$$$ g_help(x,y,1)=g_help(x,y+1,1)
!$$$$$$ g_help(x,y,2)=g_help(x,y+1,2)
!$$$$$$ g_help(x,y,3)=g_help(x,y+1,3)
!$$$$$$ g_help(x,y,4)=g_help(x,y+1,4)
!$$$$$$ g_help(x,y,5)=g_help(x,y+1,5)
!$$$$$$ g_help(x,y,6)=g_help(x,y+1,6)
!$$$$$$ g_help(x,y,7)=g_help(x,y+1,7)
!$$$$$$ g_help(x,y,8)=g_help(x,y+1,8)

sum=g_help(x,y,0)+g_help(x,y,1)+g_help(x,y,3)+g_help(x,y,4)+g_help(x,y,7)+g_help(x,y,8)
g_help(x,y,2)=wi(2)/(wi(2)+wi(5)+wi(6))*(t_cold-sum)
g_help(x,y,5)=wi(5)/(wi(2)+wi(5)+wi(6))*(t_cold-sum)
g_help(x,y,6)=wi(6)/(wi(2)+wi(5)+wi(6))*(t_cold-sum)
end do

!bounce-back
     !lower wall-barrier
  do x=10,30
     f_help(x,10,4)=f_help(x,10,2)
     f_help(x,10,7)=f_help(x,10,5)
     f_help(x,10,8)=f_help(x,10,6)
  end do

!lower wall barrier

do x=10,30
sum=g_help(x,10,0)+g_help(x,10,1)+g_help(x,10,2)+g_help(x,10,3)+g_help(x,10,5)+g_help(x,10,6)
g_help(x,10,4)=wi(4)/(wi(4)+wi(7)+wi(8))*(t_hot-sum)
g_help(x,10,7)=wi(7)/(wi(4)+wi(7)+wi(8))*(t_hot-sum)
g_help(x,10,8)=wi(8)/(wi(4)+wi(7)+wi(8))*(t_hot-sum)
end do

! lower wall barrier

  do x=60,80
     f_help(x,60,4)=f_help(x,60,2)
     f_help(x,60,7)=f_help(x,60,5)
     f_help(x,60,8)=f_help(x,60,6)
  end do

! lower wall barrier

do x=60,80
sum=g_help(x,60,0)+g_help(x,60,1)+g_help(x,60,2)+g_help(x,60,3)+g_help(x,60,5)+g_help(x,60,6)
g_help(x,60,4)=wi(4)/(wi(4)+wi(7)+wi(8))*(t_hot-sum)
g_help(x,60,7)=wi(7)/(wi(4)+wi(7)+wi(8))*(t_hot-sum)
g_help(x,60,8)=wi(8)/(wi(4)+wi(7)+wi(8))*(t_hot-sum)
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

sum=g_help(x,y,0)+g_help(x,y,2)+g_help(x,y,3)+g_help(x,y,4)+g_help(x,y,6)+g_help(x,y,7)
g_help(x,y,1)=wi(1)/(wi(1)+wi(5)+wi(8))*(t_cold-sum)
g_help(x,y,5)=wi(5)/(wi(1)+wi(5)+wi(8))*(t_cold-sum)
g_help(x,y,8)=wi(8)/(wi(1)+wi(5)+wi(8))*(t_cold-sum)
end do

!bounce-back
     !left wall-barrier
  do y=10,30
     f_help(10,y,3)=f_help(10,y,1)
     f_help(10,y,6)=f_help(10,y,8)
     f_help(10,y,7)=f_help(10,y,5)
  end do

!left wall barrier

do y=10,30
sum=g_help(10,y,0)+g_help(10,y,1)+g_help(10,y,2)+g_help(10,y,4)+g_help(10,y,5)+g_help(10,y,8)
g_help(10,y,3)=wi(3)/(wi(3)+wi(6)+wi(7))*(t_hot-sum)
g_help(10,y,6)=wi(6)/(wi(3)+wi(6)+wi(7))*(t_hot-sum)
g_help(10,y,7)=wi(7)/(wi(3)+wi(6)+wi(7))*(t_hot-sum)
end do

!left wall barrier

  do y=60,80
     f_help(60,y,3)=f_help(60,y,1)
     f_help(60,y,6)=f_help(60,y,8)
     f_help(60,y,7)=f_help(60,y,5)
  end do

!left wall barrier

do y=60,80
sum=g_help(60,y,0)+g_help(60,y,1)+g_help(60,y,2)+g_help(60,y,4)+g_help(60,y,5)+g_help(60,y,8)
g_help(60,y,3)=wi(3)/(wi(3)+wi(6)+wi(7))*(t_hot-sum)
g_help(60,y,6)=wi(6)/(wi(3)+wi(6)+wi(7))*(t_hot-sum)
g_help(60,y,7)=wi(7)/(wi(3)+wi(6)+wi(7))*(t_hot-sum)
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

sum=g_help(x,y,0)+g_help(x,y,1)+g_help(x,y,2)+g_help(x,y,4)+g_help(x,y,5)+g_help(x,y,8)
g_help(x,y,3)=wi(3)/(wi(3)+wi(6)+wi(7))*(t_cold-sum)
g_help(x,y,6)=wi(6)/(wi(3)+wi(6)+wi(7))*(t_cold-sum)
g_help(x,y,7)=wi(7)/(wi(3)+wi(6)+wi(7))*(t_cold-sum)
end do

!bounce-back
     !right wall-barrier
     
  do y=10,30
     f_help(30,y,1)=f_help(30,y,3)
     f_help(30,y,5)=f_help(30,y,7)
     f_help(30,y,8)=f_help(30,y,6)
  end do

!right wall barrier

do y=10,30
sum=g_help(30,y,0)+g_help(30,y,3)+g_help(30,y,2)+g_help(30,y,4)+g_help(30,y,7)+g_help(30,y,6)
g_help(30,y,1)=wi(1)/(wi(1)+wi(5)+wi(8))*(t_hot-sum)
g_help(30,y,5)=wi(5)/(wi(1)+wi(5)+wi(8))*(t_hot-sum)
g_help(30,y,8)=wi(8)/(wi(1)+wi(5)+wi(8))*(t_hot-sum)
end do

!right wall barrier
     
  do y=60,80
     f_help(80,y,1)=f_help(80,y,3)
     f_help(80,y,5)=f_help(80,y,7)
     f_help(80,y,8)=f_help(80,y,6)
  end do
     

!right wall barrier

do y=60,80
sum=g_help(80,y,0)+g_help(80,y,3)+g_help(80,y,2)+g_help(80,y,4)+g_help(80,y,7)+g_help(80,y,6)
g_help(80,y,1)=wi(1)/(wi(1)+wi(5)+wi(8))*(t_hot-sum)
g_help(80,y,5)=wi(5)/(wi(1)+wi(5)+wi(8))*(t_hot-sum)
g_help(80,y,8)=wi(8)/(wi(1)+wi(5)+wi(8))*(t_hot-sum)
end do


     
!--------------- Collision
do x = 0, nx-1
do y = 0, ny-1
density_loc = 0.d0
th_loc=0.d0
do i = 0, 8 
density_loc = density_loc + f_help(x,y,i)
th_loc=th_loc+g_help(x,y,i)
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
g_eq = wi(i)*th_loc*(1.d0+3*p+4.50d0*p*p-1.5*(u_x * u_x + u_y * u_y))
g(x,y,i) = g_help(x,y,i) + ( g_eq - g_help(x,y,i) ) / tau_g
end do
end do
end do
end do
call  results(nx,ny,solid,f,density,g)
print*,'********************    end     ********************'
pause	
End program
!================== Subroutines
subroutine check_density(nx,ny,f,time)  
implicit none 
integer  nx,ny,time
real*8  f(0:nx-1,0:ny-1,0:8)
integer  x,y,n
real*8 n_sum
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
subroutine results(nx,ny,solid,f,density,g)
implicit none
integer nx,ny
double precision f(0:nx-1,0:ny-1,0:8),density,g(0:nx-1,0:ny-1,0:8)
logical solid(0:nx-1,0:ny-1)
integer x,y,i,obsval
real u_x,u_y,d_loc,press,c_squ,th
c_squ = 1.d0 / 3.d0
open(11,file='Data.dat')
write(11,*) 'VARIABLES = X, Y, Ux, Uy, PRESS, OBST, T'
write(11,*) 'ZONE I=', nx, ', J=', ny, ', F=POINT'
do y = 0,nx-1
do x = 0,ny-1
if (solid(x,y)) then
obsval = 1
u_x = 0.d0
u_y = 0.d0
press = density *c_squ
th=0.0d0
else
d_loc = 0.d0
th=0.0d0
do i = 0, 8
d_loc = d_loc + f(x,y,i)
th = th + g(x,y,i)
end do
u_x = (f(x,y,1) + f(x,y,5) + f(x,y,8)-(f(x,y,3) + f(x,y,6) + f(x,y,7))) / d_loc
u_y = (f(x,y,2) + f(x,y,5) + f(x,y,6)-(f(x,y,4) + f(x,y,7) + f(x,y,8))) / d_loc
press = d_loc * c_squ
obsval = 0
end if
write(11,*) x, y, u_x, u_y, press, obsval,th
end do
end do
close(11)
return
pause
end subroutine
!==================