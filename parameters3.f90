module parameters3  

	integer,parameter::Nx = 20, Ny = 20, Nz = 20	! Simulation grid size

	! Physical size of simulation volume
	double precision,parameter::xl = -50.0d0, yl = -50.0d0, zl = -50.0d0
	double precision,parameter::xr = 50.0d0, yr = 50.0d0, zr = 50.0d0

	! Simulation parameters
	double precision,parameter::tau= 0.01d0	! Time step
	integer,parameter::time_n=50000			! Number of iterations to run
	double precision,parameter::GN = 1d-100	! Newton's constant
	double precision,parameter::c = 0.7d-97	! BEC interaction coupling strength
	double precision,parameter::N = 5d95	! Particle number density
	double precision,parameter::R = 0.8d0	! Size of initial condensate
	double precision,parameter::ex0 = 0.03d0	! Softening parameters
	double precision,parameter::ey0 = 0.09d0
	double precision,parameter::ez0 = 4.0d0
	double precision,parameter::omega0 = 0.2d0	! Initial rotation
	double precision,parameter::gamma0 = 0.0d0	! Softening parameter

	! Output control
	integer,parameter::nstep0 = 5	! number of steps of initial transient
	integer,parameter::nstep2 = 100	! how often contour plot is output

end module parameters3
