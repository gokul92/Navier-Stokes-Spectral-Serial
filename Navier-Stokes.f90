module ns

	implicit none

	include "fftw3.f"

	contains

	subroutine fourier_derivative(n, array)						! subroutine calculates derivative
											! in fourier space. 
		implicit none

		integer :: n, i
		double precision :: n1
		double complex, dimension(0:n-1) :: array

		do i = 0, n-1	 							! calculating derivative of input				
			if (i .le. n/2 - 1) then 					! function. This is done bymu(i, j, k)ltiplying
				n1 = 1.0d0*i
				array(i) = array(i)*dcmplx(0.0d0,n1)			! i*k (i is the imaginary unit) to every
			else								! element in the array. Here, k runs from
				n1 = 1.0d0*(i-n)
				array(i) = array(i)*dcmplx(0.0d0,n1)			! -l/2 to l/2 - 1. Note that the
			end if								! array is in wrapped around format only.
		end do
	
		array(n/2) = 0.0d0							! setting oddball 	wavenumber = 0

	end subroutine


	subroutine x_derivative(input_3d, nx, ny, nz)

		implicit none

		integer :: i, j, k, nx, ny, nz		
		double complex, dimension(0:nx-1, 0:ny-1, 0:nz-1) :: input_3d
		double complex, dimension(0:nx-1) :: input_x
		integer*8 ::plan_x_fwd, plan_x_bwd

		call dfftw_plan_dft_1d(plan_x_fwd, nx, input_x, input_x, fftw_forward, fftw_measure)
		call dfftw_plan_dft_1d(plan_x_bwd, nx, input_x, input_x, fftw_backward, fftw_measure)
		

		do k = 0, nz-1
			do j = 0, ny-1										! transforming to fourier space.
				input_x = input_3d(:, j, k)							! copy to 1d array from 3d array.
				call dfftw_execute_dft(plan_x_fwd, input_x, input_x)				! call fftw 1d fft subroutine.
				call fourier_derivative(nx, input_x)						! call derivative subroutine.
				input_3d(:, j, k) = input_x							! copy back to 3d input array from
			end do											! 1d array.
		end do

		input_3d = input_3d/nx										! renormalization after x-fft.

		do k = 0, nz-1
			do j = 0, ny-1										! transforming back to physical space.
				input_x = input_3d(:, j, k)
				call dfftw_execute_dft(plan_x_bwd, input_x, input_x)
				input_3d(:, j, k) = input_x
			end do
		end do

		call dfftw_destroy_plan(plan_x_fwd)
		call dfftw_destroy_plan(plan_x_bwd)

	end subroutine


	subroutine y_derivative(input_3d, nx, ny, nz)

		implicit none

		integer :: i, j, k, nx, ny, nz		
		double complex, dimension(0:nx-1, 0:ny-1, 0:nz-1) :: input_3d
		double complex, dimension(0:ny-1) :: input_y
		integer*8 ::plan_y_fwd, plan_y_bwd

		call dfftw_plan_dft_1d(plan_y_fwd, ny, input_y, input_y, fftw_forward, fftw_measure)
		call dfftw_plan_dft_1d(plan_y_bwd, ny, input_y, input_y, fftw_backward, fftw_measure)
		

		do k = 0, nz-1
			do i = 0, nx-1
				input_y = input_3d(i, :, k)
				call dfftw_execute_dft(plan_y_fwd, input_y, input_y)
				call fourier_derivative(ny, input_y)
				input_3d(i, :, k) = input_y
			end do
		end do

		input_3d = input_3d/ny

		do k = 0, nz-1
			do i = 0, nx-1
				input_y = input_3d(i, :, k)
				call dfftw_execute_dft(plan_y_bwd, input_y, input_y)
				input_3d(i, :, k) = input_y
			end do
		end do

		call dfftw_destroy_plan(plan_y_fwd)
		call dfftw_destroy_plan(plan_y_bwd)

	end subroutine


	subroutine z_derivative(input_3d, nx, ny, nz)
		
		implicit none

		integer :: i, j, k, nx, ny, nz		
		double complex, dimension(0:nx-1, 0:ny-1, 0:nz-1) :: input_3d
		double complex, dimension(0:nz-1) :: input_z
		integer*8 ::plan_z_fwd, plan_z_bwd

		call dfftw_plan_dft_1d(plan_z_fwd, nz, input_z, input_z, fftw_forward, fftw_measure)
		call dfftw_plan_dft_1d(plan_z_bwd, nz, input_z, input_z, fftw_backward, fftw_measure)
		

		do j = 0, ny-1
			do i = 0, nx-1
				input_z = input_3d(i, j, :)
				call dfftw_execute_dft(plan_z_fwd, input_z, input_z)
				call fourier_derivative(nz, input_z)
				input_3d(i, j, :) = input_z
			end do
		end do

		input_3d = input_3d/nz

		do j = 0, ny-1
			do i = 0, nx-1
				input_z = input_3d(i, j, :)
				call dfftw_execute_dft(plan_z_bwd, input_z, input_z)
				input_3d(i, j, :) = input_z
			end do
		end do

		call dfftw_destroy_plan(plan_z_fwd)
		call dfftw_destroy_plan(plan_z_bwd)

	end subroutine
	
	subroutine navier_stokes(nx, ny, nz, gama, Re, Pr, mu, p, rho, u, v, w, et, t, f_m, f_mx, f_my, f_mz, f_e)

		implicit none

		integer, intent(in) :: nx, ny, nz
		integer :: i, j, k
		double complex, intent(in), dimension(0:nx-1, 0:ny-1, 0:nz-1) :: mu, p, rho, u, v, w, et, t
		double complex, dimension(0:nx-1, 0:ny-1, 0:nz-1) :: f_m, f_mx, f_my, f_mz, f_e, temp, temp1
		double precision :: gama, Re, Pr

!***************************************************************** X Momentum ***************************************************!

		f_mx = 0.0d0
 

		temp = u
		call x_derivative(temp, nx, ny, nz)
		temp = (mu/Re)*temp
		call x_derivative(temp, nx, ny, nz)
		f_mx = f_mx + (4.0d0/3)*temp
		
		temp = u
		call y_derivative(temp, nx, ny, nz)
		temp = (mu/Re)*temp
		call y_derivative(temp, nx, ny, nz)
		f_mx = f_mx + temp
		
		temp = u
		call z_derivative(temp, nx, ny, nz)
		temp = (mu/Re)*temp
		call z_derivative(temp, nx, ny, nz)
		f_mx = f_mx + temp
		
		temp = rho*u*u
		call x_derivative(temp, nx, ny, nz)
		f_mx = f_mx - 0.50d0*temp
		
		temp = rho*u
		call x_derivative(temp, nx, ny, nz)
		f_mx = f_mx - 0.50d0*u*temp
	
		temp = u
		call x_derivative(temp, nx, ny, nz)
		f_mx = f_mx - 0.50d0*rho*u*temp
	
		temp = rho*u*v
		call y_derivative(temp, nx, ny, nz)
		f_mx = f_mx - 0.50d0*temp
	
		temp = rho*u
		call y_derivative(temp, nx, ny, nz)
		f_mx = f_mx - 0.50d0*v*temp	

		temp = v
		call y_derivative(temp, nx, ny, nz)
		f_mx = f_mx - 0.50d0*rho*u*temp
		
		temp = rho*u*w
		call z_derivative(temp, nx, ny, nz)
		f_mx = f_mx - 0.50d0*temp
			
		temp = rho*u
		call z_derivative(temp, nx, ny, nz)
		f_mx = f_mx - 0.50d0*w*temp
	
		temp = w
		call z_derivative(temp, nx, ny, nz)
		f_mx = f_mx - 0.50d0*rho*u*temp
		
		temp = v
		call x_derivative(temp, nx, ny, nz)
		temp = (mu/Re)*temp
		call y_derivative(temp, nx, ny, nz)
		f_mx = f_mx + temp

		temp = v
		call y_derivative(temp, nx, ny, nz)
		temp = (mu/Re)*temp
		call x_derivative(temp, nx, ny, nz)
		f_mx = f_mx - (2.0d0/3.0d0)*temp
	
		temp = w
		call x_derivative(temp, nx, ny, nz)
		temp = (mu/Re)*temp
		call z_derivative(temp, nx, ny, nz)
		f_mx = f_mx + temp
	
		temp = w
		call z_derivative(temp, nx, ny, nz)
		temp = (mu/Re)*temp
		call x_derivative(temp, nx, ny, nz)
		f_mx = f_mx - (2.0d0/3.0d0)*temp
	
		temp = p
		call x_derivative(temp, nx, ny, nz)
		f_mx = f_mx - temp
		!print *, "f_mx modified due to pressure is", maxval(abs(real(f_mx)))
		

!***************************************************************** Y Momentum ***************************************************!

		f_my = 0.0
	
		temp = v
		call y_derivative(temp, nx, ny, nz)
		temp = (mu/Re)*temp
		call y_derivative(temp, nx, ny, nz)
		f_my = f_my + (4.0d0/3)*temp

		temp = w
		call y_derivative(temp, nx, ny, nz)
		temp = (mu/Re)*temp
		call z_derivative(temp, nx, ny, nz)
		f_my = f_my + temp

		temp = v
		call z_derivative(temp, nx, ny, nz)
		temp = (mu/Re)*temp
		call z_derivative(temp, nx, ny, nz)
		f_my = f_my + temp

		temp = v
		call x_derivative(temp, nx, ny, nz)
		temp = (mu/Re)*temp
		call x_derivative(temp, nx, ny, nz)
		f_my = f_my + temp

		temp = u
		call y_derivative(temp, nx, ny, nz)
		temp = (mu/Re)*temp
		call x_derivative(temp, nx, ny, nz)
		f_my = f_my + temp

		temp = u
		call x_derivative(temp, nx, ny, nz)
		temp = (mu/Re)*temp
		call y_derivative(temp, nx, ny, nz)
		f_my = f_my - (2.0d0/3.0d0)*temp

		temp = w
		call z_derivative(temp, nx, ny, nz)
		temp = (mu/Re)*temp
		call y_derivative(temp, nx, ny, nz)
		f_my = f_my - (2.0d0/3.0d0)*temp

		temp = rho*v*u
		call x_derivative(temp, nx, ny, nz)
		f_my = f_my - 0.50d0*temp
	
		temp = rho*v
		call x_derivative(temp, nx, ny, nz)
		f_my = f_my - 0.50d0*u*temp

		temp = u
		call x_derivative(temp, nx, ny, nz)
		f_my = f_my - 0.50d0*rho*v*temp
	
		temp = rho*v*v
		call y_derivative(temp, nx, ny, nz)
		f_my = f_my - 0.50d0*temp

		temp = rho*v
		call y_derivative(temp, nx, ny, nz)
		f_my = f_my - 0.50d0*v*temp

		temp = v
		call y_derivative(temp, nx, ny, nz)
		f_my = f_my - 0.50d0*rho*v*temp
	
		temp = rho*v*w
		call z_derivative(temp, nx, ny, nz)
		f_my = f_my - 0.50d0*temp

		temp = rho*v
		call z_derivative(temp, nx, ny, nz)
		f_my = f_my - 0.50d0*w*temp

		temp = w
		call z_derivative(temp, nx, ny, nz)
		f_my = f_my - 0.50d0*rho*v*temp

		temp = p
		call y_derivative(temp, nx, ny, nz)
		f_my = f_my - temp
	
		

!***************************************************************** Z Momentum ***************************************************!

		f_mz = 0.0

		temp = w	
		call z_derivative(temp, nx, ny, nz)
		temp = (mu/Re)*temp
		call z_derivative(temp, nx, ny, nz)
		f_mz = f_mz + (4.0d0/3)*temp

		temp = v
		call z_derivative(temp, nx, ny, nz)
		temp = (mu/Re)*temp
		call y_derivative(temp, nx, ny, nz)
		f_mz = f_mz + temp

		temp = w
		call y_derivative(temp, nx, ny, nz)
		temp = (mu/Re)*temp
		call y_derivative(temp, nx, ny, nz)
		f_mz = f_mz + temp

		temp = w
		call x_derivative(temp, nx, ny, nz)
		temp = (mu/Re)*temp
		call x_derivative(temp, nx, ny, nz)
		f_mz = f_mz + temp

		temp = u
		call z_derivative(temp, nx, ny, nz)
		temp = (mu/Re)*temp
		call x_derivative(temp, nx, ny, nz)
		f_mz = f_mz + temp

		temp = u
		call x_derivative(temp, nx, ny, nz)
		temp = (mu/Re)*temp
		call z_derivative(temp, nx, ny, nz)
		f_mz = f_mz - (2.0d0/3.0d0)*temp

		temp = v
		call y_derivative(temp, nx, ny, nz)
		temp = (mu/Re)*temp
		call z_derivative(temp, nx, ny, nz)
		f_mz = f_mz - (2.0d0/3.0d0)*temp

		temp = rho*w*u
		call x_derivative(temp, nx, ny, nz)
		f_mz = f_mz - 0.50d0*temp

		temp = rho*w
		call x_derivative(temp, nx, ny, nz)
		f_mz = f_mz - 0.50d0*u*temp

		temp = u
		call x_derivative(temp, nx, ny, nz)
		f_mz = f_mz - 0.50d0*rho*w*temp

		temp = rho*w*v
		call y_derivative(temp, nx, ny, nz)
		f_mz = f_mz - 0.50d0*temp

		temp = rho*w
		call y_derivative(temp, nx, ny, nz)
		f_mz = f_mz - 0.50d0*v*temp

		temp = v
		call y_derivative(temp, nx, ny, nz)
		f_mz = f_mz - 0.50d0*rho*w*temp

		temp = rho*w*w
		call z_derivative(temp, nx, ny, nz)
		f_mz = f_mz - 0.50d0*temp

		temp = rho*w
		call z_derivative(temp, nx, ny, nz)
		f_mz = f_mz - 0.50d0*w*temp

		temp = w
		call z_derivative(temp, nx, ny, nz)
		f_mz = f_mz - 0.50d0*rho*w*temp

		temp = p
		call z_derivative(temp, nx, ny, nz)
		f_mz = f_mz - temp
	
	

!***************************************************************** Mass *********************************************************!

		f_m = 0.0
	
		temp = rho*u
		call x_derivative(temp, nx, ny, nz)
		f_m = f_m - temp
	
		temp = rho*v
		call y_derivative(temp, nx, ny, nz)
		f_m = f_m - temp

		temp = rho*w
		call z_derivative(temp, nx, ny, nz)
		f_m = f_m - temp

	

!***************************************************************** Energy Equation *********************************************************!

		f_e = 0.0

		temp = rho*et*u
		call x_derivative(temp, nx, ny, nz)
		f_e = f_e - 0.50d0*temp

		temp = rho*et
		call x_derivative(temp, nx, ny, nz)
		f_e = f_e - 0.50d0*u*temp

		temp = u
		call x_derivative(temp, nx, ny, nz)
		f_e = f_e - 0.50d0*rho*et*temp
	
		temp = rho*et*v
		call y_derivative(temp, nx, ny, nz)
		f_e = f_e - 0.50d0*temp
	
		temp = rho*et
		call y_derivative(temp, nx, ny, nz)
		f_e = f_e - 0.50d0*v*temp

		temp = v
		call y_derivative(temp, nx, ny, nz)
		f_e = f_e - 0.50d0*rho*et*temp

		temp = rho*et*w
		call z_derivative(temp, nx, ny, nz)
		f_e = f_e - 0.50d0*temp

		temp = rho*et
		call z_derivative(temp, nx, ny, nz)
		f_e = f_e - 0.50d0*w*temp

		temp = w
		call z_derivative(temp, nx, ny, nz)
		f_e = f_e - 0.50d0*rho*et*temp
	
		temp = p*u
		call x_derivative(temp, nx, ny, nz)
		f_e = f_e - 0.50d0*temp	

		temp = p
		call x_derivative(temp, nx, ny, nz)
		f_e = f_e - 0.50d0*u*temp

		temp = u
		call x_derivative(temp, nx, ny, nz)
		f_e = f_e - 0.50d0*p*temp

		temp = p*v
		call y_derivative(temp, nx, ny, nz)
		f_e = f_e - 0.50d0*temp	

		temp = p
		call y_derivative(temp, nx, ny, nz)
		f_e = f_e - 0.50d0*v*temp

		temp = v
		call y_derivative(temp, nx, ny, nz)
		f_e = f_e - 0.50d0*p*temp

		temp = p*w
		call z_derivative(temp, nx, ny, nz)
		f_e = f_e - 0.50d0*temp	

		temp = p
		call z_derivative(temp, nx, ny, nz)
		f_e = f_e - 0.50d0*w*temp

		temp = w
		call z_derivative(temp, nx, ny, nz)
		f_e = f_e - 0.50d0*p*temp
 
		temp = u
		call x_derivative(temp, nx, ny, nz)
		temp = (mu/Re)*u*temp
		call x_derivative(temp, nx, ny, nz)
		f_e = f_e + (2.0d0/3.0d0)*temp

		temp = u
		call x_derivative(temp, nx, ny, nz)
		temp = (mu/Re)*temp
		call x_derivative(temp, nx, ny, nz)
		f_e = f_e + (2.0d0/3.0d0)*u*temp

		temp = u
		call x_derivative(temp, nx, ny, nz)
		f_e = f_e + (2.0d0/3.0d0)*(mu/Re)*temp*temp

		temp = v
		call y_derivative(temp, nx, ny, nz)
		temp = (mu/Re)*u*temp
		call x_derivative(temp, nx, ny, nz)
		f_e = f_e - (1.0d0/3.0d0)*temp
		
		temp = v
		call y_derivative(temp, nx, ny, nz)
		temp = (mu/Re)*temp
		call x_derivative(temp, nx, ny, nz)
		f_e = f_e - (1.0d0/3.0d0)*u*temp

		temp = u
		call x_derivative(temp, nx, ny, nz)
		temp1 = v
		call y_derivative(temp1, nx, ny, nz)
		f_e = f_e - (1.0d0/3.0d0)*(mu/Re)*temp*temp1

		temp = w
		call z_derivative(temp, nx, ny, nz)
		temp = (mu/Re)*u*temp
		call x_derivative(temp, nx, ny, nz)
		f_e = f_e - (1.0d0/3.0d0)*temp

		temp = w
		call z_derivative(temp, nx, ny, nz)
		temp = (mu/Re)*temp
		call x_derivative(temp, nx, ny, nz)
		f_e = f_e - (1.0d0/3.0d0)*u*temp

		temp = u
		call x_derivative(temp, nx, ny, nz)
		temp1 = w
		call z_derivative(temp1, nx, ny, nz)
		f_e = f_e - (1.0d0/3.0d0)*(mu/Re)*temp*temp1
	
		temp = u
		call y_derivative(temp, nx, ny, nz)
		temp = (mu/Re)*v*temp
		call x_derivative(temp, nx, ny, nz)
		f_e = f_e + 0.50d0*temp

		temp = u
		call y_derivative(temp, nx, ny, nz)
		temp = (mu/Re)*temp
		call x_derivative(temp, nx, ny, nz)
		f_e = f_e + 0.50d0*v*temp

		temp = u
		call y_derivative(temp, nx, ny, nz)
		temp1 = v
		call x_derivative(temp1, nx, ny, nz)
		f_e = f_e + 0.50d0*(mu/Re)*temp*temp1

		temp = v	
		call x_derivative(temp, nx, ny, nz)
		temp = (mu/Re)*v*temp
		call x_derivative(temp, nx, ny, nz)
		f_e = f_e + 0.50d0*temp

		temp = v
		call x_derivative(temp, nx, ny, nz)
		temp = (mu/Re)*temp
		call x_derivative(temp, nx, ny, nz)
		f_e = f_e + 0.50d0*v*temp

		temp = v
		call x_derivative(temp, nx, ny, nz)
		f_e = f_e + 0.50d0*(mu/Re)*temp*temp

		temp = u
		call z_derivative(temp, nx, ny, nz)
		temp = (mu/Re)*w*temp
		call x_derivative(temp, nx, ny, nz)
		f_e = f_e + 0.50d0*temp

		temp = u
		call z_derivative(temp, nx, ny, nz)
		temp = (mu/Re)*temp
		call x_derivative(temp, nx, ny, nz)
		f_e = f_e + 0.50d0*w*temp

		temp = u
		call z_derivative(temp, nx, ny, nz)
		temp1 = w
		call x_derivative(temp1, nx, ny, nz)
		f_e = f_e + 0.50d0*(mu/Re)*temp*temp1
		
		temp = w
		call x_derivative(temp, nx, ny, nz)
		temp = (mu/Re)*w*temp
		call x_derivative(temp, nx, ny, nz)
		f_e = f_e + 0.50d0*temp
	
		temp = w
		call x_derivative(temp, nx, ny, nz)
		temp = (mu/Re)*temp
		call x_derivative(temp, nx, ny, nz)
		f_e = f_e + 0.50d0*w*temp
	
		temp = w
		call x_derivative(temp, nx, ny, nz)
		f_e = f_e + 0.50d0*(mu/Re)*temp*temp
		
		temp = v
		call x_derivative(temp, nx, ny, nz)
		temp = (mu/Re)*u*temp
		call y_derivative(temp, nx, ny, nz)
		f_e = f_e + 0.50d0*temp
	
		temp = v
		call x_derivative(temp, nx, ny, nz)
		temp = (mu/Re)*temp
		call y_derivative(temp, nx, ny, nz)
		f_e = f_e + 0.50d0*u*temp
	
		temp = u
		call y_derivative(temp, nx, ny, nz)
		temp1 = v
		call x_derivative(temp1, nx, ny, nz)
		f_e = f_e + 0.50d0*(mu/Re)*temp*temp1
	
		temp = u
		call y_derivative(temp, nx, ny, nz)
		temp = (mu/Re)*u*temp
		call y_derivative(temp, nx, ny, nz)
		f_e = f_e + 0.50d0*temp

		temp = u
		call y_derivative(temp, nx, ny, nz)
		temp = (mu/Re)*temp
		call y_derivative(temp, nx, ny, nz)
		f_e = f_e + 0.50d0*u*temp
	
		temp = u
		call y_derivative(temp, nx, ny, nz)
		f_e = f_e + 0.50d0*(mu/Re)*temp*temp
	
		temp = v
		call y_derivative(temp, nx, ny, nz)
		temp = (mu/Re)*v*temp
		call y_derivative(temp, nx, ny, nz)
		f_e = f_e + (2.0d0/3.0d0)*temp
	
		temp = v
		call y_derivative(temp, nx, ny, nz)
		temp = (mu/Re)*temp
		call y_derivative(temp, nx, ny, nz)
		f_e = f_e + (2.0d0/3.0d0)*v*temp
	
		temp = v
		call y_derivative(temp, nx, ny, nz)
		f_e = f_e + (2.0d0/3.0d0)*(mu/Re)*temp*temp
	
		temp = u
		call x_derivative(temp, nx, ny, nz)
		temp = (mu/Re)*v*temp
		call y_derivative(temp, nx, ny, nz)
		f_e = f_e - (1.0d0/3.0d0)*temp
	
		temp = u
		call x_derivative(temp, nx, ny, nz)
		temp = (mu/Re)*temp
		call y_derivative(temp, nx, ny, nz)
		f_e = f_e - (1.0d0/3.0d0)*v*temp
	
		temp = u
		call x_derivative(temp, nx, ny, nz)
		temp1 = v
		call y_derivative(temp1, nx, ny, nz)
		f_e = f_e - (1.0d0/3.0d0)*(mu/Re)*temp*temp1
	
		temp = w
		call z_derivative(temp, nx, ny, nz)
		temp = (mu/Re)*v*temp
		call y_derivative(temp, nx, ny, nz)
		f_e = f_e - (1.0d0/3.0d0)*temp

		temp = w
		call z_derivative(temp, nx, ny, nz)
		temp = (mu/Re)*temp
		call y_derivative(temp, nx, ny, nz)
		f_e = f_e - (1.0d0/3.0d0)*v*temp
	
		temp = v
		call y_derivative(temp, nx, ny, nz)
		temp1 = w
		call z_derivative(temp1, nx, ny, nz)
		f_e = f_e - (1.0d0/3.0d0)*(mu/Re)*temp*temp1
	
		temp = v
		call z_derivative(temp, nx, ny, nz)
		temp = (mu/Re)*w*temp
		call y_derivative(temp, nx, ny, nz)
		f_e = f_e + 0.50d0*temp
	
		temp = v
		call z_derivative(temp, nx, ny, nz)
		temp = (mu/Re)*temp
		call y_derivative(temp, nx, ny, nz)
		f_e = f_e + 0.50d0*w*temp
	
		temp = v
		call z_derivative(temp, nx, ny, nz)
		temp1 = w
		call y_derivative(temp1, nx, ny, nz)
		f_e = f_e + 0.50d0*(mu/Re)*temp*temp1
		
		temp = w
		call y_derivative(temp, nx, ny, nz)
		temp = (mu/Re)*w*temp
		call y_derivative(temp, nx, ny, nz)
		f_e = f_e + 0.50d0*temp
	
		temp = w
		call y_derivative(temp, nx, ny, nz)
		temp = (mu/Re)*temp
		call y_derivative(temp, nx, ny, nz)
		f_e = f_e + 0.50d0*w*temp
	
		temp = w
		call y_derivative(temp, nx, ny, nz)
		f_e = f_e + 0.50d0*(mu/Re)*temp*temp
	
		temp = w
		call x_derivative(temp, nx, ny, nz)
		temp = (mu/Re)*u*temp
		call z_derivative(temp, nx, ny, nz)
		f_e = f_e + 0.50d0*temp

		temp = w
		call x_derivative(temp, nx, ny, nz)
		temp = (mu/Re)*temp
		call z_derivative(temp, nx, ny, nz)
		f_e = f_e + 0.50d0*u*temp

		temp = w
		call x_derivative(temp, nx, ny, nz)
		temp1 = u
		call z_derivative(temp1, nx, ny, nz)
		f_e = f_e + 0.50d0*(mu/Re)*temp*temp1

		temp = u
		call z_derivative(temp, nx, ny, nz)
		temp = (mu/Re)*u*temp
		call z_derivative(temp, nx, ny, nz)
		f_e = f_e + 0.50d0*temp

		temp = u
		call z_derivative(temp, nx, ny, nz)
		temp = (mu/Re)*temp
		call z_derivative(temp, nx, ny, nz)
		f_e = f_e + 0.50d0*u*temp

		temp = u
		call z_derivative(temp, nx, ny, nz)
		f_e = f_e + 0.50d0*(mu/Re)*temp*temp

		temp = w
		call y_derivative(temp, nx, ny, nz)
		temp = (mu/Re)*v*temp
		call z_derivative(temp, nx, ny, nz)
		f_e = f_e + 0.50d0*temp

		temp = w
		call y_derivative(temp, nx, ny, nz)
		temp = (mu/Re)*temp
		call z_derivative(temp, nx, ny, nz)
		f_e = f_e + 0.50d0*v*temp

		temp = w
		call y_derivative(temp, nx, ny, nz)
		temp1 = v
		call z_derivative(temp1, nx, ny, nz)
		f_e = f_e + 0.50d0*(mu/Re)*temp*temp1

		temp = v
		call z_derivative(temp, nx, ny, nz)
		temp = (mu/Re)*v*temp
		call z_derivative(temp, nx, ny, nz)
		f_e = f_e + 0.50d0*temp
	
		temp = v
		call z_derivative(temp, nx, ny, nz)
		temp = (mu/Re)*temp
		call z_derivative(temp, nx, ny, nz)
		f_e = f_e + 0.50d0*v*temp
	
		temp = v
		call z_derivative(temp, nx, ny, nz)
		f_e = f_e + 0.50d0*(mu/Re)*temp*temp
	
		temp = w
		call z_derivative(temp, nx, ny, nz)
		temp = (mu/Re)*w*temp
		call z_derivative(temp, nx, ny, nz)
		f_e = f_e + (2.0d0/3.0d0)*temp

		temp = w
		call z_derivative(temp, nx, ny, nz)
		temp = (mu/Re)*temp
		call z_derivative(temp, nx, ny, nz)
		f_e = f_e + (2.0d0/3.0d0)*w*temp
	
		temp = w
		call z_derivative(temp, nx, ny, nz)
		f_e = f_e + (2.0d0/3.0d0)*(mu/Re)*temp*temp
	
		temp = u
		call x_derivative(temp, nx, ny, nz)
		temp = (mu/Re)*w*temp
		call z_derivative(temp, nx, ny, nz)
		f_e = f_e - (1.0d0/3.0d0)*temp
	
		temp = u
		call x_derivative(temp, nx, ny, nz)
		temp = (mu/Re)*temp
		call z_derivative(temp, nx, ny, nz)
		f_e = f_e - (1.0d0/3.0d0)*w*temp
	
		temp = u
		call x_derivative(temp, nx, ny, nz)
		temp1 = w
		call z_derivative(temp1, nx, ny, nz)
		f_e = f_e - (1.0/3)*(mu/Re)*temp*temp1
	
		temp = v
		call y_derivative(temp, nx, ny, nz)
		temp = (mu/Re)*w*temp
		call z_derivative(temp, nx, ny, nz)
		f_e = f_e - (1.0d0/3.0d0)*temp
	
		temp = v
		call y_derivative(temp, nx, ny, nz)
		temp = (mu/Re)*temp
		call z_derivative(temp, nx, ny, nz)
		f_e = f_e - (1.0d0/3.0d0)*w*temp
	
		temp = v
		call y_derivative(temp, nx, ny, nz)
		temp1 = w
		call z_derivative(temp1, nx, ny, nz)
		f_e = f_e - (1.0d0/3.0d0)*(mu/Re)*temp*temp1
	
		temp = t
		call x_derivative(temp, nx, ny, nz)
		temp = temp*(mu/((gama-1)*Re*Pr))
		call x_derivative(temp, nx, ny, nz)
		f_e = f_e + temp

		temp = t
		call y_derivative(temp, nx, ny, nz)
		temp = temp*(mu/((gama-1)*Re*Pr))
		call y_derivative(temp, nx, ny, nz)
		f_e = f_e + temp 

		temp = t
		call z_derivative(temp, nx, ny, nz)
		temp = temp*(mu/((gama-1)*Re*Pr))
		call z_derivative(temp, nx, ny, nz)
		f_e = f_e + temp

	end subroutine

end module

program ns_main

	use ns																		! use module to take derivative
																			! and fft along x, y & z
	implicit none

	integer, parameter :: nx = 8, ny = 8, nz = 128, nt = 100											! grid points along x, y & z.
	integer :: i, j, k, time
	double complex, dimension(0:nx-1, 0:ny-1, 0:nz-1) :: mu, p, rho, u, v, w, et, t, f_m, f_mx, f_my, f_mz, f_e	
	double complex, dimension(0:nx-1, 0:ny-1, 0:nz-1) :: solution
	double complex, dimension(0:nx-1, 0:ny-1, 0:nz-1) :: m_rhs, mx_rhs, my_rhs, mz_rhs, e_rhs			
	double complex, dimension(0:nx-1, 0:ny-1, 0:nz-1) :: temp_mu, temp_p, temp_rho, temp_u, temp_v, temp_w, temp_et, temp_t
	double precision, dimension(0:nx - 1) :: x													! 1d arrays to store
	double precision, dimension(0:ny - 1) :: y													! physical space
	double precision, dimension(0:nz - 1) :: z													! points.
	double precision, parameter :: t_star = 300.0d0, gama = 1.4d0, Re = 500.0, Pr = 0.7d0	
	double precision, parameter :: rho_star = 1.2d0, rho_one = 0.3d0*rho_star, rho_two = 3.0d0*rho_star
	double precision, parameter :: R_air = 288.16, c_star = sqrt(gama*R_air*t_star), mu_star = 1.983d0*0.0001			
	double precision, parameter :: l_star = Re*mu_star/(rho_star*c_star), bval = 3.0d0/l_star
	double precision, parameter :: x_one = 2.327*l_star, x_two = 3.956*l_star
	double precision, parameter :: pi = acos(-1.0d0), dx = 2.0d0*pi/nx, dy = 2.0d0*pi/ny, dz = 2.0d0*pi/nz, dt = 0.01
	character(10) :: filename
	real :: num

	
	do i = 0, nx-1								
		x(i) = i*dx							
	end do									
										
	do j = 0, ny-1								
		y(j) = j*dy  
	end do

	do k = 0, nz-1
		z(k) = k*dz
	end do

	do k = 0, nz-1
		do j = 0, ny-1
			do i = 0, nx-1										
				rho(i, j, k) = rho_one/rho_star + 0.5d0*((rho_two-rho_one)/rho_star)*				&
						(tanh(bval*(z(k)*l_star-x_one))- tanh(bval*(z(k)*l_star-x_two)))
				u(i, j, k) = 0.0d0
				v(i, j, k) = 0.0d0 
				w(i, j, k) = 0.0d0
				t(i, j, k) = 1.0d0 
				p(i, j, k) = rho(i, j, k)*t(i, j, k)/gama
				mu(i, j, k) = t(i, j, k)**0.670d0
				!solution(i, j, k) = -(0.5d0*bval*l_star*(rho_two-rho_one)/(gama*rho_star))			&
				!			*(1.0d0/(cosh(bval*(x(i)*l_star-x_one)))**2 				&
				!			- 1.0d0/(cosh(bval*(x(i)*l_star-x_two)))**2)
			end do
		end do
	end do
	
	et = t/(gama*(gama-1.0d0)) + 0.5d0*(u**2 + v**2 + w**2)

	!open(unit = 22, file = "rho@t=0.0.dat")
	!do i = 0, nx-1
	!	write(22, *) real(rho(i, 1, 1))
	!end do
	!close(22)

	do time = 1, nt

		m_rhs = 0.0d0
		mx_rhs = 0.0d0
		my_rhs = 0.0d0
		mz_rhs = 0.0d0
		e_rhs = 0.0d0

		call navier_stokes(nx, ny, nz, gama, Re, Pr, mu, p, rho, u, v, w, et, t, f_m, f_mx, f_my, f_mz, f_e)

		m_rhs = m_rhs + (1.0d0/6.0d0)*f_m
		mx_rhs = mx_rhs + (1.0d0/6.0d0)*f_mx
		my_rhs = my_rhs + (1.0d0/6.0d0)*f_my
		mz_rhs = mz_rhs + (1.0d0/6.0d0)*f_m
		e_rhs = e_rhs + (1.0d0/6.0d0)*f_e

		temp_rho = rho + 0.5d0*dt*f_m
		temp_u = (rho*u + 0.5d0*dt*f_mx)/temp_rho
		temp_v = (rho*v + 0.5d0*dt*f_my)/temp_rho
		temp_w = (rho*w + 0.5d0*dt*f_mz)/temp_rho
		temp_et = (rho*et + 0.5d0*dt*f_e)/temp_rho
		temp_t = (temp_et - 0.5d0*(temp_u**2 + temp_v**2 + temp_w**2))*gama*(gama-1.0d0)
		temp_p = temp_rho*temp_t/gama
		temp_mu = temp_t**(0.67d0)

		call navier_stokes(nx, ny, nz, gama, Re, Pr, temp_mu, temp_p, temp_rho, temp_u, temp_v, temp_w, temp_et,		&
				 temp_t, f_m, f_mx, f_my, f_mz, f_e)

		m_rhs = m_rhs + (1.0d0/3.0d0)*f_m
		mx_rhs = mx_rhs + (1.0d0/3.0d0)*f_mx
		my_rhs = my_rhs + (1.0d0/3.0d0)*f_my
		mz_rhs = mz_rhs + (1.0d0/3.0d0)*f_mz
		e_rhs = e_rhs + (1.0d0/3.0d0)*f_e

		temp_rho = rho + 0.5d0*dt*f_m
		temp_u = (rho*u + 0.5d0*dt*f_mx)/temp_rho
		temp_v = (rho*v + 0.5d0*dt*f_my)/temp_rho
		temp_w = (rho*w + 0.5d0*dt*f_mz)/temp_rho
		temp_et = (rho*et + 0.5d0*dt*f_e)/temp_rho
		temp_t = (temp_et - 0.5d0*(temp_u**2 + temp_v**2 + temp_w**2))*gama*(gama-1.0d0)
		temp_p = temp_rho*temp_t/gama
		temp_mu = temp_t**(0.67d0)

		call navier_stokes(nx, ny, nz, gama, Re, Pr, temp_mu, temp_p, temp_rho, temp_u, temp_v, temp_w, temp_et,		&
				 temp_t, f_m, f_mx, f_my, f_mz, f_e)
	
		m_rhs = m_rhs + (1.0d0/3.0d0)*f_m
		mx_rhs = mx_rhs + (1.0d0/3.0d0)*f_mx
		my_rhs = my_rhs + (1.0d0/3.0d0)*f_my
		mz_rhs = mz_rhs + (1.0d0/3.0d0)*f_mz
		e_rhs = e_rhs + (1.0d0/3.0d0)*f_e
	
	
		temp_rho = rho + dt*f_m
		temp_u = (rho*u + dt*f_mx)/temp_rho
		temp_v = (rho*v + dt*f_my)/temp_rho
		temp_w = (rho*w + dt*f_mz)/temp_rho
		temp_et = (rho*et + dt*f_e)/temp_rho
		temp_t = (temp_et - 0.5d0*(temp_u**2 + temp_v**2 + temp_w**2))*gama*(gama-1.0d0)
		temp_p = temp_rho*temp_t/gama
		temp_mu = temp_t**(0.67d0)

		call navier_stokes(nx, ny, nz, gama, Re, Pr, temp_mu, temp_p, temp_rho, temp_u, temp_v, temp_w, temp_et,		&
				 temp_t, f_m, f_mx, f_my, f_mz, f_e)

	
		m_rhs = m_rhs + (1.0d0/6.0d0)*f_m
		mx_rhs = mx_rhs + (1.0d0/6.0d0)*f_mx
		my_rhs = my_rhs + (1.0d0/6.0d0)*f_my
		mz_rhs = mz_rhs + (1.0d0/6.0d0)*f_mz
		e_rhs = e_rhs + (1.0d0/6.0d0)*f_e
	
		temp_rho = rho		! assigning to temp_rho, the value of rho at nth iterative level.
		
		rho = temp_rho + dt*m_rhs
		u = (temp_rho*u + dt*mx_rhs)/rho
		v = (temp_rho*v + dt*my_rhs)/rho
		w = (temp_rho*w + dt*mz_rhs)/rho
		et = (temp_rho*et + dt*e_rhs)/rho
		t = (et -0.5d0*(u**2 + v**2 + w**2))*gama*(gama-1.d0)
		p = rho*t/gama
		mu = t**0.67d0

		if (mod(time,nt) .eq. 0) then	
			num = time*0.01
			write(filename,'(F10.2)') num	
			open(unit = 20, file = "rho,t,p,u@t="//filename//".dat")
			do k = 0, nz-1
				write(20, '(5E30.18)') z(k), real(rho(1, 1, k)), real(t(1, 1, k)), real(p(1, 1, k)), real(w(1, 1, k))
			end do
			close(20)		
		end if
	end do

end program





	











	















	
