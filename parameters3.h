#define SQ(x) ((x)*(x))

// User-configurable parameters begin here

// Grid size
#define Nx 80
#define Ny 80
#define Nz 80

// Largest grid dimension -- do not alter
#define NN ((Nx>Ny?Nx:Ny)>Nz?(Nx>Ny?Nx:Ny):Nz)

// Physical size of simulation volume
const double xl = -50.0f, yl = -50.0f, zl = -50.0f;
const double xr = 50.0f, yr = 50.0f, zr = 50.0f;

// Simulation parameters
const double tau = 0.01f;		// Time step
const int time_n = 500;		// Number of iterations to run
const double GN = 1e-100;		// Newton's constant
const double c = 0.7e-97;		// BEC interaction coupling strength
const double N = 5e95;			// Particle number density
const double R = 0.8;			// Size of initial condensate
const double ex0 = 0.03f;		// Softening parameters
const double ey0 = 0.09f;
const double ez0 = 4.0f;
const double omega0 = 0.1f;		// Initial rotation
const double gamma0 = 0.00f;		// Softening parameter

// Output control
const int nstep0 = 5;	// number of steps of initial transient without output
const int nstep2 = 100;	// every how many steps contour plot is output

