// BEC3P.cpp : Self-gravitating Bose-Einstein condensate
//

#include "stdafx.h"

#include "parameters3.h"

#pragma warning(disable:4996)

using namespace std;

const int Nn = (Nx + 1) * (Ny + 1) * (Nz + 1);

complex<double> eye;
double dt, t, dx, dy, dz, idx2, idy2, idz2;
double pi;

complex<double> psi[Nn];
complex<double> psi_n[Nn];
double dens[Nn];
double phi[Nn];
double UU[Nn];
double res[Nn];


#define ijk(i,j,k) ((((i) * (Ny + 1)) + (j)) * (Nz + 1) + (k))

#define psi(i,j,k) (psi[ijk(i,j,k)])
#define psi_n(i,j,k) (psi_n[ijk(i,j,k)])
#define density(i,j,k) (dens[ijk(i,j,k)])
#define phi(i,j,k) (phi[ijk(i,j,k)])
#define U(i,j,k) (UU[ijk(i,j,k)])

// Forward declarations
double init(double, double, int, int, int, double, double, double);
double get_normsimp();

void get_phi();
void get_U(double, double, double, double, double);
void calc_rhs_x(complex<double>, double);
void solve_x(complex<double>, double);
void calc_rhs_y(complex<double>, double);
void solve_y(complex<double>, double);
void calc_rhs_z(complex<double>);
void solve_z(complex<double>);
double energy(double, double, double);
void movie(int);
void thomas(complex<double> *, complex<double> *, complex<double> *,
			complex<double> *, int m);
void get_density();

 
int _tmain(int argc, _TCHAR* argv[])
{
	complex<double> xi; // coefficient for time derivative after softening. Here consider non-dessipative case, and xi=i
	double mu; // chemical potential
	double ex, ey, ez, omega, gamma; //Softening parameters and rotational speed
	double norm, norm_n;
	int i, j, k, itime;

	// Fixed parameters

	pi = 4 * atan((double)1);
	eye = complex<double>(0, 1);
	dx = (xr - xl) / Nx;
	dy = (yr - yl) / Ny;
	dz = (zr - zl) / Nz;

	idx2 = 1 / dx / dx;
	idy2 = 1 / dy / dy;
	idz2 = 1 / dz / dz;
     
	// Initial conditions

	mu = (double)0.5 * sqrt(4 * c / pi);

	// Zero arrays
	memset(phi, 0, sizeof(phi));
	memset(psi, 0, sizeof(psi));
	memset(psi_n, 0, sizeof(psi));
	memset(UU, 0, sizeof(UU));

	dt = tau;
	ex = ex0;
	ey = ey0;
	ez = ez0;
	omega = omega0;
	gamma = gamma0;

	xi = (eye + gamma) / (gamma * gamma + 1);

	// Initial state

	t = 0.0;
	itime = 0;
	for (k = 0; k <= Nz; k++)
		for (j = 0; j <= Ny; j++)
			for (i = 0; i <= Nx; i++)
		psi(i, j, k) = complex<double>(init(mu, c, i, j, k, ex, ey, ez), 0);
	get_density();	// We need the initial density for the grav. pot.
	movie(itime);	// Output data for contour plots

	// Boundary conditions psi=0
	for (i = 0; i <= Nx; i++)
		for (j = 0; j <= Ny; j++)
	{
		psi(i, j, 0) = 0;
		psi(i, j, Nz) = 0;
	}
	for (j = 0; j <= Ny; j++)
		for (k = 0; k <= Nz; k++)
	{
		psi(0, j, k) = 0;
		psi(Nx, j, k) = 0;
	}
	for (k = 0; k <= Nz; k++)
		for (i = 0; i <= Nx; i++)
	{
		psi(i, 0, k) = 0;
		psi(i, Ny, k) = 0;
	}

	// Output initial state

	norm_n = get_normsimp();

	// get initial U 

	get_phi();
	get_U(mu, c, ex, ey, ez);

	// Time loop

	for (itime = 1; itime <= time_n; itime++)
	{
		bool bLoop;
		int a = 0;
		t = t + dt;
		calc_rhs_x(xi, omega);	// Find psi'=Rx(psi_old)
		norm =  get_normsimp();
		
		for (bLoop = true; bLoop;)	// Solve Lx(psi'')=psi' and iterate
		{							// due to nonlinear term in U
			solve_x(xi, omega);
			get_U(mu, c, ex, ey, ez);
			double norm1 = get_normsimp();
			if (fabs(norm1 - norm) < fabs(1e-4 * norm)) bLoop = false;
			norm = norm1;
		}
		calc_rhs_y(xi, omega);		// Find psi'''=Ry(psi'')
		norm =  get_normsimp();
		for (bLoop = true; bLoop;)
		{
			solve_y(xi, omega);		// Solve Ly(psi'''')=psi'''
			get_U(mu, c, ex, ey, ez);
			double norm1 = get_normsimp();
			if (fabs(norm1 - norm) < fabs(1e-4 * norm)) bLoop = false;
			norm = norm1;
		}
		calc_rhs_z(xi);				// Find psi'''''=Rz(psi'''')
		norm =  get_normsimp();
		for (bLoop = true; bLoop;)
		{
			solve_z(xi);			// Solve Lz(psi_new)=psi'''''
			get_U(mu, c, ex, ey, ez);
			double norm1 = get_normsimp();
			if (fabs(norm1 - norm) < fabs(1e-4 * norm)) bLoop = false;
			norm = norm1;
		}
		get_density();
		get_phi();
		get_U(mu, c, ex, ey, ez);	// Find U correspdng to psi_new
		printf("N=%6d, t=%8.3lf, E=%11.4lg, P=%11.4lg"
			   "\n", itime, t, energy(ex, ey, ez), norm);
		fflush(stdout);
		if (itime > nstep0 && itime % nstep2 == 0)
			movie(itime);			// Output data for contour plots
	}	// Close time loop

	return true;
}

//*********************************************************************
double init(double mu, double c, int i, int j, int k,
			double ex, double ey, double ez)		// Initial state
{
	double F, x, y, z, r, rf;

	x = xl + i * dx;
	y = yl + j * dy;
	z = zl + k * dz;
	F = 0.0;
	r = sqrt((1 + ex) * x * x + (1 + ey) * y * y + (1 + ez) * z * z);

	// The grav. potential for the n=0 Lane-Emden solution,
	// rho = R^2 - r^2.
	rf = xr;
	if (rf > yr) rf = yr;
	if (rf > zr) rf = zr;
	rf /= 3;
	if (r * r / R < 3 * rf * rf)
		F = sqrt((rf * rf - r * r / R / 3) * N);
	return F;
}

//**********************************************************************
// Maybe the energy in lab frame, i.e. not calculating the rotational energy
//**********************************************************************
double energy(double ex, double ey, double ez)	// Energy
{
	double E, f, mdpsisq, VV;
	int i, j, k;

	E = 0;
	VV = 0;
	f = 1.0;
	for (k = 1; k < Nz; k++)
		for (j = 1; j < Ny; j++)
			for (i = 1; i < Nx; i++)
	{
				mdpsisq = SQ(real(psi(i, j, k))) + SQ(imag(psi(i, j, k)));
				VV = VV + phi(i, j, k); // shouldn't * density?
				E = E + (double)0.25 * (mdpsisq - 1) * (mdpsisq - 1); // shouldn't * coupling constant?
				E = E + (double)0.5 * (SQ(real(psi(i + 1, j, k) -
					psi(i - 1, j, k)) / (2 * dx)) + SQ(imag(psi(i + 1, j, k) -
					psi(i - 1, j, k)) / (2 * dx)) +
					SQ(real(psi(i, j + 1, k) - psi(i, j - 1, k)) / (2 * dy))
					+ SQ(imag(psi(i, j + 1, k) - psi(i, j - 1, k)) /
					(2 * dy)) + SQ(real(psi(i, j, k + 1) -
					psi(i, j, k - 1)) / (2 * dz)) +
					SQ(imag(psi(i, j, k + 1) - psi(i, j, k - 1)) / (2 * dz))); // This is dimensionless kinetic energy
	}
	return E * f * dx * dy * dz + VV;
}

//**********************************************************************
//The function get_U computes the potential term (without rotation) in 
// the Gross-Pitaevskii equation: H = kinetic - U - Omega*L
//**********************************************************************
void get_U(double mu, double c, double ex, double ey, double ez)	// Find U
{
	int i, j, k;

	for (k = 0; k <= Nz; k++)
		for (j = 0; j <= Ny; j++)
			for (i = 0; i <= Nx; i++)
		U(i, j, k) = mu - c * (SQ(real(psi(i, j, k))) + SQ(imag(psi(i, j, k))))
					- phi(i, j, k);
}

//**********************************************************************
// Calculate d_n, which is used to find psi_n+1 by Thomas algorirhm
// By ADI, we propogate along x for each y. So d/dy=0 and we fix j for each iteration
// Notedd - i*omega(xd/dy-yd/dx)/3 can be approximated by finite difference
// Define U = mu-k\\rho-phi_gravity
// b = 1-i*(dt/3)/2(V/3+1/dx^2) = 1-i*(dt/3)/2(-U/3 +1/dx^2)=1+foo4*(U/3-idx2)
// a = i*(dt/3)/(4*dx^2)-i(dt/3)*i*omega*y/3/(4dx)=foo3-foo2
// c = i*(dt/3)/(4*dx^2)-i(dt/3)*i*omega*y/3/(4dx)=foo3+foo2
//**********************************************************************
void calc_rhs_x(complex<double> xi, double omega)	// Find psi_n= Rx(psi)
{
	double y, idx2;
	complex<double> foo1, foo2, foo3, foo4;
	int i, j, k;

	idx2 = 1 / (dx * dx);
	foo1 = eye * omega * xi * dt / 3. / (4 * dx);  
	foo3 = (double)0.25 * xi * dt / 3. * idx2; 
	foo4 = (double)0.5 * xi * dt / 3.;
	for (i = 1; i < Nx; i++)
		for (j = 1; j < Ny; j++)
	{
		y = yl + j * dy; // For each y, propagate along x
		foo2 = foo1 * y; // Rotational energy (i*omega*ydx)
		for (k = 1; k < Nz; k++)
		{
			psi_n(i, j, k) = (foo3 - foo2) * psi(i - 1, j, k) +
					((double)1 + foo4 * (U(i, j, k) / 3. - idx2)) *
					psi(i, j, k) + (foo3 + foo2) * psi(i + 1, j, k);
		}
	}
}

//**********************************************************************
// Calculate d_n, which is used to find psi_n+1 by Thomas algorirhm
// By ADI, we propogate along y for each x. So d/dx=0 and we fix i for each iteration
// Notedd - i*omega(xd/dy-yd/dx)/3 can be approximated by finite difference
// Define U = mu-k\\rho-phi_gravity
// b = 1-i*(dt/3)/2(V/3+1/dy^2) = 1-i*(dt/3)/2(phi_gravity/3 +1/dy^2)=1+foo4*(U/3-idy2)
// a = i*(dt/3)/(4*dy^2)+i(dt/3)*i*omega*x/3/(4dy)=foo3+foo2
// c = i*(dt/3)/(4*dy^2)-i(dt/3)*i*omega*x/3/(4dy)=foo3-foo2
//**********************************************************************
void calc_rhs_y(complex<double> xi, double omega)		// Find psi_n=Ry(psi)
{
	double x, idy2;
	complex<double> foo1, foo2, foo3, foo4;
	int i, j, k;

	idy2 = 1 / (dy * dy);
	foo1 = eye * omega * xi * dt / 3. / (4 * dy);
	foo3 = (double)0.25 * xi * dt / 3. * idy2;
	foo4 = (double)0.5 * xi * dt / 3.;
	for (i = 1; i < Nx; i++)
	{
		x = xl + i * dx; // For each x, propagate along y
		foo2 = foo1 * x;
		for (j = 1; j < Ny; j++)
			for (k = 1; k < Nz; k++)
		{
			psi_n(i, j, k) = (foo3 + foo2) * psi(i, j - 1, k) +
				((double)1 + foo4 * (U(i, j, k) / 3. - idy2)) * psi(i, j, k)
				+ (foo3 - foo2) * psi(i, j + 1, k);
		}
	}
}

//**********************************************************************
// Calculate d_n, which is used to find psi_n+1 by Thomas algorirhm
// By ADI, we propogate along z, so d/dx = d/dy = 0
// Define U = mu-k\\rho-phi_gravity
// Question: Why U(i,j,k) this time has negative sign? I have changed it from bb - foo4 * U(i, j, k) / 3. to bb + foo4 * U(i, j, k) / 3.
//**********************************************************************
void calc_rhs_z(complex<double> xi)	// Find psi_n=Ry(psi)
{
	complex<double> aa, bb, foo4;
	int i, j, k;

	aa = xi * dt / 3. / (4 * dz * dz);
	bb = (double)1 - xi * dt / 3. / (2 * dz * dz);
	foo4 = (double)0.5 * xi * dt / 3.;
	for (i = 1; i < Nx; i++)
		for (j = 1; j < Ny; j++)
			for (k = 1; k < Nz; k++)
		psi_n(i, j, k) = aa * psi(i, j, k - 1) +
						(bb + foo4 * U(i, j, k) / 3.) * psi(i, j, k) + aa * psi(i, j, k + 1);
}

//**********************************************************************
// For each z and y, propogate along x for dt/3. Find psi_n+1 by Thomas algorirhm
//**********************************************************************
void solve_x(complex<double> xi, double omega)	// Solve Lx(psi)=psi_n for psi
{
	double y;
	complex<double> aa[Nx], bb[Nx], cc[Nx], rr[Nx];
	complex<double> foo1, foo2, foo3, foo4;
	int i, j, k;

	foo1 = -xi * dt / 3. / (4 * dx * dx);
	foo2 = eye * omega * xi * dt / 3. / (4 * dx);
	foo3 = (double)1 + xi * dt / 3. / (2 * dx * dx);
	foo4 = (double)0.5 * xi * dt / 3.;
	for (k = 1; k < Nz; k++)
		for (j = 1; j < Ny; j++)
	{
		y = yl + j * dy;
		aa[1] = foo1 + foo2 * y;
		bb[1] = foo3;
		cc[1] = foo1 - foo2 * y;
		for (i = 2; i < Nx; i++)
		{
			aa[i] = aa[1];
			bb[i] = bb[1] - foo4 * U(i, j, k) / 3.;
			cc[i] = cc[1];
			rr[i] = psi_n(i, j, k);
		}
		bb[1] = bb[1] - foo4 * U(1, j, k) / 3.;
		rr[1] = psi_n(1, j, k);
		thomas(aa, bb, cc, rr, Nx - 1);
		for (i = 1; i < Nx; i++) psi(i, j, k) = rr[i];
		psi(0, j, k) = 0;
		psi(Nx, j, k) = 0;
	}
}

//**********************************************************************
// For each z and x, propogate along y for dt/3. Find psi_n+1 by Thomas algorirhm
//**********************************************************************
void solve_y(complex<double> xi, double omega)	// Solve Ly(psi)=psi_n for psi
{
	double x;
	complex<double> aa[Ny], bb[Ny], cc[Ny], rr[Ny];
	complex<double> foo1, foo2, foo3, foo4;
	int i, j, k;

	foo1 = -xi * dt / 3. / (4 * dy * dy);
	foo2 = eye * omega * xi * dt / 3. / (4 * dy);
	foo3 = (double)1 + xi * dt / 3. / (2 * dy * dy);
	foo4 = (double)0.5 * xi * dt / 3.;
	for (i = 1; i < Nx; i++)
		for (k = 1; k < Nz; k++)
	{
		x = xl + i * dx;
		aa[1] = foo1 - foo2 * x;
		bb[1] = foo3;
		cc[1] = foo1 + foo2 * x;
		for (j = 2; j < Ny; j++)
		{
			aa[j] = aa[1];
			bb[j] = bb[1] - foo4 * U(i, j, k) / 3.;
			cc[j] = cc[1];
			rr[j] = psi_n(i, j, k);
		}
		bb[1] = bb[1] - foo4 * U(i, 1, k) / 3.;
		rr[1] = psi_n(i, 1, k);
		thomas(aa, bb, cc, rr, Ny - 1);
		for (j = 1; j < Ny; j++) psi(i, j, k) = rr[j];
		psi(i, 0, k) = 0;
		psi(i, Ny, k) = 0;
	}
}

//**********************************************************************
// For each x and y, propogate along z for dt/3. Find psi_n+1 by Thomas algorirhm
//**********************************************************************
void solve_z(complex<double> xi)	// Solve Ly(psi)=psi_n for psi
{
	complex<double> aa[Nz], bb[Nz], cc[Nz], rr[Nz];
	complex<double> foo1, foo2, foo4;
	int i, j, k;

	foo1 = -xi * dt / 3. / (4 * dz * dz);
	foo2 = (double)1 + xi * dt / 3. / (2 * dz * dz);
	foo4 = (double)0.5 * xi * dt / 3.;
	for (i = 1; i < Nx; i++)
		for (j = 1; j < Ny; j++)
	{
		aa[1] = foo1;
		bb[1] = foo2;
		cc[1] = foo1;
		for (k = 2; k < Nz; k++)
		{
			aa[k] = aa[1];
			bb[k] = bb[1] - foo4 * U(i, j, k) / 3.;
			cc[k] = cc[1];
			rr[k] = psi_n(i, j, k);
		}
		bb[1] = bb[1] - foo4 * U(i, j, 1) / 3.;
		rr[1] = psi_n(i, j, 1);
		thomas(aa, bb, cc, rr, Nz - 1);
		for (k = 1; k < Nz; k++)
		{
			psi(i, j, k) = rr[k];
		}
		psi(i, j, 0) = 0;
		psi(i, j, Nz) = 0;
	}
}

//**********************************************************************
// The function thomas solves a tridiagonal system of linear equations using the Thomas algorithm.
// The arrays l,d,u and r contain respectively the subdiagonal, the diagonal, the superdiagonal and the right hand side.
// m is the size of the system. At output the solution is put in r, while the original right hand side is destroyed.
//**********************************************************************
void thomas(complex<double> *l,	// Tridiagonal matrix solver with Thomas
			complex<double> *d,	// algorithm. The arrays l,d,u and r contain
			complex<double> *u,	// respectively the subdiagonal, the diagonal,
			complex<double> *r,	// the superdiagonal and the right hand side.
			int m)				// m is the size of the system. At output the
{								// solution is put in r, while the original
								// right hand side is destroyed.
	int j;
	complex<double> a;
    
	for (j = 2; j <= m; j++)
	{
		a = -l[j] / d[j - 1];
		d[j] = d[j] + a * u[j - 1];
		r[j] = r[j] + a * r[j - 1];
	}
	r[m] = r[m] / d[m];
	for (j = m - 1; j >= 1; j--)
		r[j] = (r[j] - u[j] * r[j + 1]) / d[j];
}

//***************************************************************************
// Find the norm using Composite Simpson's 1/3 rule
// Integrate along z and then y and then x direction
// Check https://en.wikipedia.org/wiki/Simpson%27s_rule#Simpson%27s_3/8_rule
//***************************************************************************
double get_normsimp()	
{
	static double normz[Nx + 1][Ny + 1], normy[Nx + 1], norm;
	int i, j, k;

	norm = 0.0;
	for (i = 0; i < Nx; i++)
		for (j = 0; j < Ny; j++)
	{
		normz[i][j] = SQ(real(psi(i, j, 0))) + SQ(imag(psi(i, j, 0))) +
					SQ(real(psi(i, j, Nz - 1))) + SQ(imag(psi(i, j, Nz - 1)));
		for (k = 1; k <= Nz - 3; k += 2)
			normz[i][j] += 4 * (SQ(real(psi(i, j, k))) +
								SQ(imag(psi(i, j, k)))) +
						   2 * (SQ(real(psi(i, j, k + 1))) +
								SQ(imag(psi(i, j, k + 1))));
	}
	for (i = 0; i < Nx; i++)
	{
		normy[i] = normz[i][0] + normz[i][Ny - 1];
		for (j = 1; j <= Ny - 3; j += 2)
			normy[i] += 4 * normz[i][j] + 2 * normz[i][j + 1];
	}
	norm = normy[0] + normy[Nx - 1];
	for (i = 1; i <= Nx - 3; i += 2)
		norm += 4 * normy[i] + 2 * normy[i + 1];
	return norm * dx * dy * dz / 27;
}
    
//**********************************************************************
// For each Z, print information of x, y
//**********************************************************************
void movieZ(int itime)	// Outputs files
{
	double x, y, z, phpsi, depsi;
	char ch[17], cs[17], cu[17];
	int i, j, k;
	FILE *file8, *file10, *file11;

	sprintf(ch, "densZ%07d.dat", itime);
	sprintf(cs, "phasZ%07d.dat", itime);
	sprintf(cu, "gravZ%07d.dat", itime);
	file8 = fopen(ch, "w");
	file10 = fopen(cs, "w");
	file11 = fopen(cu, "w");

	// Output files for contour plots

	for (k = 0; k <= Nz; k += 2)
		for (j = 0; j <= Ny; j += 2)
	{
		for (i = 0; i <= Nx; i += 2)
		{
			x = xl + i * dx;
			y = yl + j * dy;
			z = zl + k * dz;
			depsi = SQ(real(psi(i, j, k))) + SQ(imag(psi(i, j, k)));
			phpsi = atan2(imag(psi(i, j, k)) + 1e-15, real(psi(i, j, k)));
			fprintf(file8, "%lf %lf %lf %lg\n", x, y, z, depsi);
			fprintf(file11, "%lf %lf %lf %lg\n", x, y, z, phi(i, j, k));
			fprintf(file10, "%lf %lf %lf %lg\n", x, y, z, phpsi);
		}
		fprintf(file8, "\n");	// Leave blank line for Gnuplot before next j
		fprintf(file11, "\n");
		fprintf(file10, "\n");
	}

	fclose(file8);
	fclose(file10);
	fclose(file11);
}

//**********************************************************************
// For each X, print information of y, z
//**********************************************************************
void movieX(int itime)	// Outputs files
{
	double x, y, z, phpsi, depsi;
	char ch[17], cs[17], cu[17];
	int i, j, k;
	FILE *file8, *file10, *file11;

	sprintf(ch, "densX%07d.dat", itime);
	sprintf(cs, "phasX%07d.dat", itime);
	sprintf(cu, "gravX%07d.dat", itime);
	file8 = fopen(ch, "w");
	file10 = fopen(cs, "w");
	file11 = fopen(cu, "w");

	// Output files for contour plots

	for (i = 0; i <= Nx; i += 2)
		for (j = 0; j <= Ny; j += 2)
	{
		for (k = 0; k <= Nz; k += 2)
		{
			x = xl + i * dx;
			y = yl + j * dy;
			z = zl + k * dz;
			depsi = SQ(real(psi(i, j, k))) + SQ(imag(psi(i, j, k)));
			phpsi = atan2(imag(psi(i, j, k)) + 1e-15, real(psi(i, j, k)));
			fprintf(file8, "%lf %lf %lf %lg\n", x, y, z, depsi);
			fprintf(file11, "%lf %lf %lf %lg\n", x, y, z, phi(i, j, k));
			fprintf(file10, "%lf %lf %lf %lg\n", x, y, z, phpsi);
		}
		fprintf(file8, "\n");	// Leave blank line for Gnuplot before next j
		fprintf(file11, "\n");
		fprintf(file10, "\n");
	}

	fclose(file8);
	fclose(file10);
	fclose(file11);
}
//**********************************************************************
// For each Y, print information of x, z
//**********************************************************************
void movieY(int itime)	// Outputs files
{
	double x, y, z, phpsi, depsi;
	char ch[17], cs[17], cu[17];
	int i, j, k;
	FILE *file8, *file10, *file11;

	sprintf(ch, "densY%07d.dat", itime);
	sprintf(cs, "phasY%07d.dat", itime);
	sprintf(cu, "gravY%07d.dat", itime);
	file8 = fopen(ch, "w");
	file10 = fopen(cs, "w");
	file11 = fopen(cu, "w");

	// Output files for contour plots

	for (j = 0; j <= Ny; j += 2)
		for (i = 0; i <= Nx; i += 2)
	{
		for (k = 0; k <= Nz; k += 2)
		{
			x = xl + i * dx;
			y = yl + j * dy;
			z = zl + k * dz;
			depsi = SQ(real(psi(i, j, k))) + SQ(imag(psi(i, j, k)));
			phpsi = atan2(imag(psi(i, j, k)) + 1e-15, real(psi(i, j, k)));
			fprintf(file8, "%lf %lf %lf %lg\n", x, y, z, depsi);
			fprintf(file11, "%lf %lf %lf %lg\n", x, y, z, phi(i, j, k));
			fprintf(file10, "%lf %lf %lf %lg\n", x, y, z, phpsi);
		}
		fprintf(file8, "\n");	// Leave blank line for Gnuplot before next j
		fprintf(file11, "\n");
		fprintf(file10, "\n");
	}

	fclose(file8);
	fclose(file10);
	fclose(file11);
}

void movie(int itime)
{
	movieZ(itime);
	movieX(itime);
	movieY(itime);
}

//**********************************************************************
// The function get_phi computes the gravitational potential
// corresponding to the current values of \\psi using the relaxation
// method.
//**********************************************************************
void get_phi()	
{
	double rtot1, rtot2;
	const double h2 = dx * dx * dy * dy * dz * dz / 2 / (dx * dx * dy * dy
					+ dy * dy * dz * dz + dz * dz * dx * dx);
	int i, j, k, n;
	const double GN4pi = 4 * pi * GN;

	rtot1 = -1;
	rtot2 = 1;
	for (n = 0; fabs(rtot2 - rtot1) > 1e-5 * fabs(rtot1); n++)
	{
		rtot2 = rtot1;
		for (k = 1; k < Nz; k++)
			for (j = 1; j < Ny; j++)
				for (i = 1; i < Nx; i++)
		{
			res[ijk(i, j, k)] = ((phi(i + 1, j, k) + phi(i - 1, j, k)) * idx2 +
								 (phi(i, j - 1, k) + phi(i, j + 1, k)) * idy2 +
								 (phi(i, j, k - 1) + phi(i, j, k + 1)) * idz2 -
								 GN4pi * density(i, j, k)) * h2;
		}
		memcpy(phi, res, sizeof(phi));
		rtot1 = 0;
		for (i = 0; i <= Nx; i++)
			for (j = 0; j < Ny; j++)
				for (k = 0; k < Nz; k++)
			rtot1 += phi(i, j, k);
	}
}

void get_density()
{
	int i, j, k;

	for (i = 0; i <= Nx; i++)
		for (j = 0; j <= Ny; j++)
			for (k = 0; k <= Nz; k++)
		density(i, j, k) = SQ(real(psi(i, j, k))) + SQ(imag(psi(i, j, k)));
}
