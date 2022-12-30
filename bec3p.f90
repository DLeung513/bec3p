! BEC3P.f90 : Self-gravitating Bose-Einstein condensate

program main
	USE parameters3
	complex*16::psi(0:Nx, 0:Ny, 0:Nz), psi_n(0:Nx, 0:Ny, 0:Nz)
	double precision::U(0:Nx, 0:Ny, 0:Nz), density(0:Nx, 0:Ny, 0:Nz)
	double precision::phi(0:Nx, 0:Ny, 0:Nz)
	double precision::res(0:Nx, 0:Ny, 0:Nz)
	double precision::energy, init
	complex*16::dt, xi, eye
	double precision::t, mu, dx, dy, dz, x(0:Nx), y(0:Ny), z(0:Nz)
	double precision::ex, ey, ez, omega, gamma, lambda
	double precision::pi, norm, norm_n, norm1, sqnorm
	integer::i, j, k, itime
	logical bLoop
	common/timeblock/t, dt
	common/increments/dx, dy, dz
	double precision h2, GN4pi, idx2, idy2, idz2
	common/constants/eye, h2, GN4pi, idx2, idy2, idz2
	common/residuals/res

	! Fixed parameters

	pi = 4.0d0 * datan(1.0d0)
	eye = (0.0, 1.0)
	dx = (xr - xl) / Nx
	dy = (yr - yl) / Ny
	dz = (zr - zl) / Nz

	h2 = dx * dx * dy * dy * dz * dz / 2 / (dx * dx * dy * dy &
                 + dy * dy * dz * dz + dz * dz * dx * dx)
	GN4pi = 4.0d0 * pi * GN
	idx2 = 1 / dx / dx
	idy2 = 1 / dy / dy
	idz2 = 1 / dz / dz

     
	! Initial conditions

	mu = 0.5d0 * sqrt(4.0d0 * c / pi)

	psi = 0	! Zero arrays 
	psi_n = 0
	U = 0
	phi = 0
	res = 0

	dt = tau
	ex = ex0
	ey = ey0
	ez = ez0
	omega = omega0
	gamma = gamma0

	xi = (eye + gamma) / (gamma * gamma + 1.0d0)

	! Initial state

	t = 0.0
	itime = 0
	do i = 0, Nx
		do j = 0, Ny
			do k = 0, Nz
		psi(i, j, k) = init(mu, i, j, k, ex, ey, ez)
		density(i, j, k) = psi(i, j, k) * dconjg(psi(i, j, k))
			end do 
		end do 
	end do
	call movie(0, psi, phi)	! Output data for contour plots

	! Boundary conditions psi=0
	do i = 0, Nx
		do j = 0, Ny
			psi(i, j, 0) = 0
			psi(i, j, Nz) = 0
		end do 
	end do
	do j = 0, Ny
		do k = 0, Nz
			psi(0, j, k) = 0
			psi(Nx, j, k) = 0
		end do
	end do
	do k = 0, Nz
		do i = 0, Nx
			psi(i, 0, k) = 0
			psi(i, Ny, k) = 0
		end do
	end do

	! Output initial state

	call get_normsimp(psi, norm_n)

	! get initial U 

	call get_phi(phi, density)
	call get_U(U, psi, mu, ex, ey, ez, phi)

	! Time loop

	do itime = 1, time_n
		t = t + dt
		call calc_rhs_x(xi,omega, psi, psi_n, U)	! Find psi'=Rx(psi_old)
		call get_normsimp(psi, norm)
		bLoop = .true.
		do while (bLoop)	! Solve Lx(psi'')=psi' and iterate
									! due to nonlinear term in U
			call solve_x(xi, omega, psi, psi_n, U)
			call get_U(U, psi, mu, ex, ey, ez, phi)
			call get_normsimp(psi, norm1)
			if (abs(norm - norm1) < abs(1d-4 * norm)) bLoop = .false.
			norm = norm1
		end do
		call calc_rhs_y(xi, omega, psi, psi_n, U)	! Find psi'''=Ry(psi'')
		call get_normsimp(psi, norm)
		bLoop = .true.
		do while (bLoop)
			call solve_y(xi, omega, psi, psi_n, U)	! Solve Ly(psi'''')=psi'''
			call get_U(U, psi, mu, ex, ey, ez, phi)
			call get_normsimp(psi, norm1)
			if (abs(norm - norm1) < abs(1d-4 * norm)) bLoop = .false.
			norm = norm1
		end do
		call calc_rhs_z(xi, psi, psi_n, U)			! Find psi'''''=Rz(psi'''')
		call get_normsimp(psi, norm)
		bLoop = .true.
		do while (bLoop)
			call solve_z(xi, psi, psi_n, U)			! Solve Lz(psi_new)=psi'''''
			call get_U(U, psi, mu, ex, ey, ez, phi)
			call get_normsimp(psi, norm1)
			if (abs(norm - norm1) < abs(1d-4 * norm)) bLoop = .false.
			norm = norm1
		end do
		density = psi * dconjg(psi)
		call get_phi(phi, density)
		call get_U(U, psi, mu, ex, ey, ez, phi)	! Find U correspdng to psi_new
		write(6, '("N=", i6, ", t=", g8.3, ", E=", g11.4, ", P=", g11.4 &
			  )'), itime, t, energy(ex, ey, ez, psi, phi), norm !&
		flush (6)
		if (itime .gt. nstep0 .and. mod(itime, nstep2) .eq. 0) then
			call movie(itime, psi, phi)	! Output data for contour plots
		endif
	end do	! Close time loop
end    

!*********************************************************************
pure double precision function init(mu, i, j, k, ex, ey, ez)
													! Initial state
	USE parameters3
	integer, intent(in)::i, j, k
	double precision, intent(in)::ex, ey, ez, mu
	double precision::F, x, y, z, rr, rf
	double precision::dx, dy, dz
	common/increments/dx, dy, dz

	x = xl + i * dx
	y = yl + j * dy
	z = zl + k * dz
	rr = sqrt((1 + ex) * x * x + (1 + ey) * y * y + (1 + ez) * z * z)
	F = 0d0
	! The grav. potential for the n=0 Lane-Emden solution
	! rho = R^2 - rr^2.
	rf = xr
	if (rf > yr) rf = yr
	if (rf > zr) rf = zr
	rf = rf / 3d0
	if (rr * rr / R < 3d0 * rf * rf) then
		F = sqrt((rf * rf - rr * rr / R / 3d0) * N)
	endif
	init = F
	return
end

!**********************************************************************
double precision function energy(ex, ey, ez, psi, phi)	! Energy
	USE parameters3
	complex*16 psi(0:Nx, 0:Ny, 0:Nz)
	double precision, intent(in)::phi(0:Nx, 0:Ny, 0:Nz)
	double precision E, f, mdpsisq, dx, dy, dz, V, ex, ey, ez
	integer i, j, k
	common/increments/dx, dy, dz

	E = 0
	V = 0
	f = 1.0
	do k = 1, Nz - 1
		do j = 1, Ny - 1
			do i = 1, Nx - 1
				mdpsisq = psi(i, j, k) * dconjg(psi(i, j, k))
				V = V + phi(i, j, k)
				E = E + 0.25d0 * (mdpsisq - 1) * (mdpsisq - 1)
				E = E + 0.5d0 * ((dreal(psi(i + 1, j, k) - psi(i - 1, j, k)) / &
				  (2.0d0 * dx)) ** 2 + (dimag(psi(i + 1, j, k) - &
				  psi(i - 1, j, k)) / (2.0d0 * dx)) ** 2 + &
				  (dreal(psi(i, j + 1, k) - psi(i, j - 1, k)) / (2.0d0 * dy)) &
					** 2 + (dimag(psi(i, j + 1, k) - psi(i, j - 1, k)) / &
					(2.0d0 * dy)) ** 2 + (dreal(psi(i, j, k + 1) - &
					psi(i, j, k - 1)) / (2.0d0 * dz)) ** 2 + &
					(dimag(psi(i, j, k + 1) - psi(i, j, k - 1)) &
					/ (2.0d0 * dz)) ** 2)
			end do	   
		end do
	end do
	energy = E * f * dx * dy * dz + V
	return
end

!**********************************************************************
subroutine get_U(U, psi, mu, ex, ey, ez, phi)	! Find U
	USE parameters3
	complex*16 psi(0:Nx, 0:Ny, 0:Nz)
	double precision U(0:Nx, 0:Ny, 0:Nz), mu
	double precision, intent(in)::phi(0:Nx, 0:Ny, 0:Nz)
	double precision ex, ey, ez
	integer i, j, k

	U = mu - c * psi * dconjg(psi)
	U = U - phi
	return
end

!**********************************************************************
subroutine calc_rhs_x(xi, omega, psi, psi_n, U)	! Find psi_n= Rx(psi)
	USE parameters3
	complex*16 psi(0:Nx, 0:Ny, 0:Nz), psi_n(0:Nx, 0:Ny, 0:Nz)
	complex*16 eye, xi, dt
	double precision dx, dy, dz, t, y, omega, U(0:Nx, 0:Ny, 0:Nz)
	double precision idx2
	complex*16 foo1, foo2, foo3, foo4
	integer i, j, k
	common/timeblock/t, dt
	common/increments/dx, dy, dz
	common/constants/eye

	idx2 = 1.0d0 / (dx * dx)
	foo1 = eye * omega * xi * dt / 3.0d0 / (4.0d0 * dx)
	foo3 = 0.25d0 * xi * dt / 3.0d0 * idx2
	foo4 = 0.5d0 * xi * dt / 3.0d0
	do i=1,Nx-1
		do j=1,Ny-1
			y = yl + j * dy
			foo2 = foo1 * y
			do k=1,Nz-1
				psi_n(i, j, k) = (foo3 - foo2) * psi(i - 1, j, k) + &
					(1 + foo4 * (U(i, j, k) / 3.0d0 - idx2)) * &
					psi(i, j, k) + (foo3 + foo2) * psi(i + 1, j, k)
			end do
		end do
	end do
	return
end

!**********************************************************************
subroutine calc_rhs_y(xi, omega, psi, psi_n, U)	! Find psi_n=Ry(psi)
	USE parameters3
	complex*16 psi(0:Nx, 0:Ny, 0:Nz), psi_n(0:Nx, 0:Ny, 0:Nz)
	complex*16 eye, xi, dt
	double precision dx, dy, dz, t, x, omega, U(0:Nx, 0:Ny, 0:Nz)
	double precision idy2
	complex*16 foo1, foo2, foo3, foo4
	integer i, j, k
	common/timeblock/t, dt
	common/increments/dx, dy, dz
	common/constants/eye

	idy2 = 1.0d0 / (dy * dy)
	foo1 = eye * omega * xi * dt / 3.0d0 / (4.0d0 * dy)
	foo3 = 0.25d0 * xi * dt / 3.0d0 * idy2
	foo4 = 0.5d0 * xi * dt / 3.0d0
	do i = 1, Nx - 1
		x = xl + i * dx
		foo2 = foo1 * x
		do j = 1, Ny - 1
			do k = 1, Nz - 1
				psi_n(i, j, k) = (foo3 + foo2) * psi(i, j - 1, k) + &
					(1 + foo4 * (U(i, j, k) / 3.0d0 - idy2)) * psi(i, j, k) + &
					(foo3 - foo2) * psi(i, j + 1, k)
			end do
		end do
	end do
	return
end

!**********************************************************************
subroutine calc_rhs_z(xi, psi, psi_n, U)	! Find psi_n=Ry(psi)
	USE parameters3
	complex*16 psi(0:Nx, 0:Ny, 0:Nz), psi_n(0:Nx, 0:Ny, 0:Nz)
	double precision U(0:Nx, 0:Ny, 0:Nz)
	complex*16 xi, dt, aa, bb, foo4
	double precision dx, dy, dz, t
	integer i, j, k
	common/timeblock/t, dt
	common/increments/dx, dy, dz

	aa = xi * dt / 3.0d0 / (4.0d0 * dz * dz)
	bb = 1 - xi * dt / 3.0d0 / (2.0d0 * dz * dz)
	foo4 = 0.5d0 * xi * dt / 3.0d0
	do i = 1, Nx - 1
		do j = 1, Ny - 1
			do k = 1, Nz - 1
				psi_n(i, j, k) = aa * psi(i, j, k - 1) + &
					(bb - foo4 * U(i, j, k) / 3.0d0) * psi(i, j, k) + &
					aa * psi(i, j, k + 1)
			end do
		end do
	end do
	return
end

!**********************************************************************
subroutine solve_x(xi, omega, psi, psi_n, U)	! Solve Lx(psi)=psi_n for psi
	USE parameters3
	complex*16 psi(0:Nx, 0:Ny, 0:Nz), psi_n(0:Nx, 0:Ny, 0:Nz)
	double precision U(0:Nx, 0:Ny, 0:Nz), dx, dy, dz, t, omega, y
	complex*16 aa(Nx - 1), bb(Nx - 1), cc(Nx - 1), rr(Nx - 1)
	complex*16 dt, xi, eye
	complex*16 foo1, foo2, foo3, foo4
	integer i, j, k
	common/timeblock/t, dt
	common/increments/dx, dy, dz
	common/constants/eye

	foo1 = -xi * dt / 3.0d0 / (4.0d0 * dx * dx)
	foo2 = eye * omega * xi * dt / 3.0d0 / (4.0d0 * dx)
	foo3 = (1.0d0, 0.0d0) + xi * dt / 3.0d0 / (2.0d0 * dx * dx)
	foo4 = 0.5d0 * xi * dt / 3.0d0
	do k = 1, Nz - 1
		do j = 1, Ny - 1
			y = yl + j * dy
			aa(1) = foo1 + foo2 * y
			bb(1) = foo3
			cc(1) = foo1 - foo2 * y
			do i = 2, Nx - 1
				aa(i) = aa(1)
				bb(i) = bb(1) - foo4 * U(i, j, k) / 3.0d0
				cc(i) = cc(1)
				rr(i) = psi_n(i, j, k)
			end do
			bb(1) = bb(1) - foo4 * U(1, j, k) / 3.0d0
			rr(1) = psi_n(1, j, k)
			call thomas(aa, bb, cc, rr, Nx - 1)
			do i = 1, Nx - 1
				psi(i, j, k) = rr(i)
			end do
			psi(0, j, k) = 0
			psi(Nx, j, k) = 0
		end do
	end do
	return
end

!**********************************************************************
subroutine solve_y(xi, omega, psi, psi_n, U)	! Solve Ly(psi)=psi_n for psi

	USE parameters3
	complex*16 psi(0:Nx, 0:Ny, 0:Nz), psi_n(0:Nx, 0:Ny, 0:Nz)
	double precision dx, dy, dz, t, omega, x
	double precision U(0:Nx, 0:Ny, 0:Nz)
	complex*16 aa(Ny - 1), bb(Ny - 1), cc(Ny - 1), rr(Ny - 1)
	complex*16 dt, xi, eye
	complex*16 foo1, foo2, foo3, foo4
	integer i, j, k
	common/timeblock/t, dt
	common/increments/dx, dy, dz
	common/constants/eye

	foo1 = -xi * dt / 3.0d0 / (4.0d0 * dy * dy)
	foo2 = eye * omega * xi * dt / 3.0d0 / (4.0d0 * dy)
	foo3 = 1 + xi * dt / 3.0d0 / (2.0d0 * dy * dy)
	foo4 = xi * dt * 0.5d0 / 3.0d0
	do i = 1, Nx - 1
		do k = 1, Nz - 1
			x = xl + i * dx
			aa(1) = foo1 - foo2 * x
			bb(1) = foo3
			cc(1) = foo1 + foo2 * x
			do j = 2, Ny - 1
				aa(j) = aa(1)
				bb(j) = bb(1) - foo4 * U(i, j, k) / 3.0d0
				cc(j) = cc(1)
				rr(j) = psi_n(i, j, k)
			end do
			bb(1) = bb(1) - foo4 * U(i, 1, k) / 3.0d0
			rr(1) = psi_n(i, 1, k)
			call thomas(aa, bb, cc, rr, Ny - 1)
			do j = 1, Ny - 1
				psi(i, j, k) = rr(j)
			end do
			psi(i, 0, k) = 0
			psi(i, Ny, k) = 0
		end do
	end do
	return
end

!**********************************************************************
subroutine solve_z(xi, psi, psi_n, U)	! Solve Ly(psi)=psi_n for psi
	USE parameters3
	complex*16 psi(0:Nx, 0:Ny, 0:Nz), psi_n(0:Nx, 0:Ny, 0:Nz)
	double precision dx, dy, dz, t
	double precision U(0:Nx, 0:Ny, 0:Nz)
	complex*16 aa(Nz - 1), bb(Nz - 1), cc(Nz - 1), rr(Nz - 1)
	complex*16 dt, xi
	complex*16 foo1, foo2, foo4
	integer i, j, k
	common/timeblock/t, dt
	common/increments/dx, dy, dz

	foo1 = -xi * dt / 3.0d0 / (4.0d0 * dz * dz)
	foo2 = 1 + xi * dt / 3.0d0 / (2.0d0 * dz * dz)
	foo4 = xi * dt * 0.5d0 / 3.0d0
	do i = 1, Nx - 1
		do j = 1, Ny - 1
			aa(1) = foo1
			bb(1) = foo2
			cc(1) = foo1
			do k = 2, Nz - 1
				aa(k) = aa(1)
				bb(k) = bb(1) - foo4 * U(i, j, k) / 3.0d0
				cc(k) = cc(1)
				rr(k) = psi_n(i, j, k)
			end do
			bb(1) = bb(1) - foo4 * U(i, j, 1) / 3.0d0
			rr(1) = psi_n(i, j, 1)
			call thomas(aa, bb, cc, rr, Nz - 1)
			do k = 1, Nz - 1
				psi(i, j, k) = rr(k)
			end do
			psi(i, j, 0) = 0
			psi(i, j, Nz) = 0
		end do
	end do
	return
end
!**********************************************************************
subroutine thomas(l,d,u,r,m)	! Tridiagonal matrix solver with Thomas
								! algorithm. The arrays l,d,u and r contain
								! respectively the subdiagonal, the diagonal,
								! the superdiagonal and the right hand side.
								! m is the size of the system. At output the
								! solution is put in r, while the original
								! right hand side is destroyed.
	integer m, j
	complex*16 l(m), d(m), u(m), r(m), a
    
	do j = 2, m
		a = -l(j) / d(j - 1)
		d(j) = d(j) + a * u(j - 1)
		r(j) = r(j) + a * r(j - 1)
	end do
	r(m) = r(m) / d(m)
	do j = m - 1, 1, -1
		r(j) = (r(j) - u(j) * r(j + 1)) / d(j)
	end do
	return
end

!***************************************************************************
subroutine get_normsimp(psi,norm)	! Find the norm
	USE parameters3
	complex*16, intent(in)::psi(0:Nx, 0:Ny, 0:Nz)
	double precision, intent(out)::norm
	double precision::dx, dy, dz, normz(0:Nx,0:Ny), normy(0:Nx)
	integer::i, j, k
	common/increments/dx, dy, dz

	norm=0.0
	do i=0, Nx - 1
		do j=0, Ny - 1
			normz(i, j) = psi(i, j, 0) * dconjg(psi(i, j, 0)) + &
							psi(i, j, Nz - 1) * dconjg(psi(i, j, Nz - 1))
			do k = 1, Nz - 3, 2
				normz(i, j) = normz(i, j) + &
							4.0d0 * psi(i, j, k) * dconjg(psi(i, j, k)) + &
							2.0d0 * psi(i, j, k + 1) * dconjg(psi(i, j, k + 1))
			end do
		end do
	end do
	do i = 0, Nx - 1 
		normy(i) = normz(i, 0) + normz(i, Ny - 1)         
		do j = 1, Ny - 3, 2
			normy(i) = normy(i) + 4.0d0 * normz(i, j) + 2.0d0 * normz(i, j + 1)	  
		end do
	end do
	norm = normy(0) + normy(Nx - 1)
	do i = 1, Nx - 3, 2
		norm = norm + 4.0d0 * normy(i) + 2.0d0 * normy(i + 1)	
	end do
	norm = norm * dx * dy * dz / 27.0d0
	return
end subroutine get_normsimp
    
!**********************************************************************
subroutine movieZ(itime,psi, phi)	! Outputs files
	USE parameters3
	complex*16,intent(in)::psi(0:Nx, 0:Ny, 0:Nz)
	double precision,intent(in)::phi(0:Nx, 0:Ny, 0:Nz)
	double precision::x, y, z, dx, dy, dz, phpsi, depsi
	character*16 ch, cs, cu
	integer,intent(in)::itime
	integer::i, j, k
	common/increments/dx, dy, dz

	ch = 'densZ*******.dat'
	cs = 'phasZ*******.dat'
	cu = 'gravZ*******.dat'
	write(ch(6:12), "(i7.7)") itime
	write(cs(6:12), "(i7.7)") itime
	write(cu(6:12), "(i7.7)") itime
	open(unit = 8, file = ch, status = 'unknown')
	open(unit = 10, file = cs, status = 'unknown')
	open(unit = 11, file = cu, status = 'unknown')

	! Output files for contour plots

	do k = 0, Nz, 2
		do j = 0, Ny, 2
			do i=0, Nx, 2
				x = xl + i * dx
				y = yl + j * dy
				z = zl + k * dz
				depsi = zabs(psi(i, j, k)) ** 2
				phpsi = datan2(dimag(psi(i, j, k)) + 1.0d-15, &
								dreal(psi(i, j, k)))
				write(8, "(4g14.5e3)") x, y, z, depsi
				write(11, "(4g14.5e3)") x, y, z, phi(i, j, k)
				write(10, "(4g14.5e3)") x, y, z, phpsi
			end do
			write(8, "(3f14.5)") ! Leave blank line for Gnuplot before next j
			write(11, "(3f14.5)")
			write(10, "(3f14.5)") 
		end do
	end do

	close(8)
	close(10)
	close(11)
	return
end

subroutine movieX(itime,psi, phi)	! Outputs files
	USE parameters3
	complex*16,intent(in)::psi(0:Nx, 0:Ny, 0:Nz)
	double precision,intent(in)::phi(0:Nx, 0:Ny, 0:Nz)
	double precision::x, y, z, dx, dy, dz, phpsi, depsi
	character*16 ch, cs, cu
	integer,intent(in)::itime
	integer::i, j, k
	common/increments/dx, dy, dz

	ch = 'densX*******.dat'
	cs = 'phasX*******.dat'
	cu = 'gravX*******.dat'
	write(ch(6:12), "(i7.7)") itime
	write(cs(6:12), "(i7.7)") itime
	write(cu(6:12), "(i7.7)") itime
	open(unit = 8, file = ch, status = 'unknown')
	open(unit = 10, file = cs, status = 'unknown')
	open(unit = 11, file = cu, status = 'unknown')

	! Output files for contour plots

	do i=0, Nx, 2
		do j = 0, Ny, 2
			do k = 0, Nz, 2
				x = xl + i * dx
				y = yl + j * dy
				z = zl + k * dz
				depsi = zabs(psi(i, j, k)) ** 2
				phpsi = datan2(dimag(psi(i, j, k)) + 1.0d-15, &
								dreal(psi(i, j, k)))
				write(8, "(4g14.5e3)") x, y, z, depsi
				write(11, "(4g14.5e3)") x, y, z, phi(i, j, k)
				write(10, "(4g14.5e3)") x, y, z, phpsi
			end do
			write(8, "(3f14.5)") ! Leave blank line for Gnuplot before next j
			write(11, "(3f14.5)")
			write(10, "(3f14.5)") 
		end do
	end do

	close(8)
	close(10)
	close(11)
	return
end

subroutine movieY(itime,psi, phi)	! Outputs files
	USE parameters3
	complex*16,intent(in)::psi(0:Nx, 0:Ny, 0:Nz)
	double precision,intent(in)::phi(0:Nx, 0:Ny, 0:Nz)
	double precision::x, y, z, dx, dy, dz, phpsi, depsi
	character*16 ch, cs, cu
	integer,intent(in)::itime
	integer::i, j, k
	common/increments/dx, dy, dz

	ch = 'densY*******.dat'
	cs = 'phasY*******.dat'
	cu = 'gravY*******.dat'
	write(ch(6:12), "(i7.7)") itime
	write(cs(6:12), "(i7.7)") itime
	write(cu(6:12), "(i7.7)") itime
	open(unit = 8, file = ch, status = 'unknown')
	open(unit = 10, file = cs, status = 'unknown')
	open(unit = 11, file = cu, status = 'unknown')

	! Output files for contour plots

	do j = 0, Ny, 2
		do i=0, Nx, 2
			do k = 0, Nz, 2
				x = xl + i * dx
				y = yl + j * dy
				z = zl + k * dz
				depsi = zabs(psi(i, j, k)) ** 2
				phpsi = datan2(dimag(psi(i, j, k)) + 1.0d-15, &
								dreal(psi(i, j, k)))
				write(8, "(4g14.5e3)") x, y, z, depsi
				write(11, "(4g14.5e3)") x, y, z, phi(i, j, k)
				write(10, "(4g14.5e3)") x, y, z, phpsi
			end do
			write(8, "(3f14.5)") ! Leave blank line for Gnuplot before next j
			write(11, "(3f14.5)")
			write(10, "(3f14.5)") 
		end do
	end do

	close(8)
	close(10)
	close(11)
	return
end

subroutine movie(itime, psi, phi)
	USE parameters3
	complex*16,intent(in)::psi(0:Nx, 0:Ny, 0:Nz)
	double precision,intent(in)::phi(0:Nx, 0:Ny, 0:Nz)
	integer,intent(in)::itime

	call movieX(itime, psi, phi)
	call movieY(itime, psi, phi)
	call movieZ(itime, psi, phi)
	return
end

!*****************************************************************************   
subroutine get_velocity(psi, vx, vy, vz, density)	! Calculate fluid velocity
	USE parameters3    
	complex*16,intent(in)::psi(0:Nx, 0:Ny, 0:Nz)
	complex*16::eye
	double precision, intent(out)::vx(1:Nx, 1:Ny, 1:Nz), vy(1:Nx, 1:Ny, 1:Nz)
	double precision, intent(out)::vz(1:Nx, 1:Ny, 1:Nz)
	double precision, intent(in)::density(0:Nx, 0:Ny, 0:Nz)
	double precision::dx, dy, dz
	integer::i, j, k
	common/increments/dx, dy, dz
	common/constants/eye

	do k = 1, Nz
		do j = 1, Ny
			do i = 1, Nx
				if (density(i, j, k) .gt. 1e-44) then
					vx(i, j, k) = (-eye/(2.0d0 * (psi(i, j, k) * &
						dconjg(psi(i, j, k)) + 1.0d-15))) * &
						(dconjg(psi(i, j, k)) * ((psi(i+1, j, k) - &
						psi(i - 1, j, k)) / (2.0d0 * dx)) - psi(i, j, k) * &
						((dconjg(psi(i + 1, j, k)) - &
						dconjg(psi(i - 1, j, k))) / (2.0d0 * dx)))
					vy(i, j, k) = (-eye/(2.0d0 * (psi(i, j, k) * &
						dconjg(psi(i, j, k)) + 1.0d-15))) * &
						(dconjg(psi(i, j, k)) * ((psi(i, j + 1, k) - &
						psi(i, j - 1, k)) / 2.0d0 * dy) - psi(i, j, k) * &
						((dconjg(psi(i, j + 1, k)) - &
						dconjg(psi(i, j - 1, k))) / (2.0d0 * dy)))
					vz(i, j, k) = (-eye / (2.0d0 * (psi(i, j, k) * &
						dconjg(psi(i, j, k)) + 1.0d-15))) * &
						(dconjg(psi(i, j, k)) * ((psi(i, j, k + 1) - &
						psi(i, j, k - 1)) / (2.0d0 * dz)) - psi(i, j, k) * &
						((dconjg(psi(i, j, k + 1)) - &
						dconjg(psi(i, j, k - 1))) / (2.0d0 * dz)))
				endif 	   	                     
			end do        
		end do
	end do
	return
end subroutine get_velocity

!**********************************************************************
subroutine get_density(psi, density)
	USE parameters3
	complex*16,intent(in)::psi(0:Nx, 0:Ny, 0:Nz)
	double precision,intent(out)::density(0:Nx, 0:Ny, 0:Nz)
	integer::i, j, k

	density = psi * dconjg(psi)
	return                  
end subroutine get_density

!**********************************************************************
subroutine get_phi(phi, density)	! Grav. potential via Poisson's Eq.
	USE parameters3
	double precision, intent(in)::density(0:Nx, 0:nY, 0:Nz)
	double precision::phi(0:Nx, 0:Ny, 0:Nz)
	double precision::res(0:Nx, 0:Ny, 0:Nz), rtot1, rtot2
	integer::i, j, k, nn
	double precision dx, dy, dz
	common/increments/dx, dy, dz
	common/residuals/res
	complex*16::eye
	double precision h2, GN4pi, idx2, idy2, idz2
	common/constants/eye, h2, GN4pi, idx2, idy2, idz2

	rtot1 = -1
	rtot2 = 1
	nn = 0
	do while (abs(rtot2 - rtot1) > 1e-5 * abs(rtot1))
		rtot2 = rtot1
		forall (k = 1 : Nz - 1, j = 1 : Ny - 1, i = 1 : Nx - 1)
			res(i, j, k) = ((phi(i + 1, j, k) + phi(i - 1, j, k)) * idx2 + &
							(phi(i, j - 1, k) + phi(i, j + 1, k)) * idy2 + &
							(phi(i, j, k - 1) + phi(i, j, k + 1)) * idz2 - &
							GN4pi * density(i, j, k)) * h2
		end forall
		phi = res
		rtot1 = sum(phi)
		nn = nn + 1
	end do
	return
end subroutine get_phi
