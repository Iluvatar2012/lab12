#include <iostream>
#include <fstream>
#include <complex>
#include <sstream>
//-----------------------------------
using namespace std;
//-----------------------------------
typedef complex<double> cmplx;
//-----------------------------------
void init( cmplx* const psi0, const double eta, const double sigma, const double dx,
          const int Nx);

void writeToFile(const cmplx* const v, const string s, const double dx,
                 const int Nx, const double xmin);

void lin_step (cmplx* const psi0, const double dt, const double dx, const int Nx);
void nonlin_step (cmplx* const psi0, const double dt, const int Nx);
//-----------------------------------
int main(){

	const int Nx = 4000;
	const double L = 800;
	const double xmin = 0;
	const double Tend = 50;
	const double dx = L / (Nx - 1);
	const double dt = dx  / 10;
	const int Na = 10;
	int Nk = int(Tend / Na / dt + 0.5);

	const double eta = 0.2;

	stringstream strm;

	cmplx* psi0 = new cmplx[Nx];

	init(psi0, eta, dx, dt,Nx);

	writeToFile(psi0,"psi_0", dx,Nx,xmin);


	for (int i = 1; i <= Na; i++) {

		lin_step(psi0, dt/2., dx, Nx);

		for (int j = 1; j <= Nk-1; j++) {
			nonlin_step (psi0, dt, Nx);
			lin_step (psi0, dt, dx, Nx);
		}

		nonlin_step (psi0,  dt, Nx);
		lin_step (psi0, dt/2., dx, Nx);

		strm.str("");
		strm << "psi_" << i;
		writeToFile(psi0,strm.str(), dx,Nx,xmin);
	}

	return 0;
}
//-----------------------------------
void writeToFile(const cmplx* const v, const string s, const double dx,
                 const int Nx, const double xmin)
{
	ofstream out(s.c_str());
	for(int i=0; i<Nx; i++){
		double x = xmin + i * dx;
		out << x << "\t" << norm(v[i]) << "\t" << v[i].real() << "\t" << v[i].imag() << endl;
	}
	out.close();
}
//-----------------------------------

void init( cmplx* const psi0, const double eta,  const double dx, const double dt,
          const int Nx)
{
	const double x0 = dx*Nx * 0.5;
	const double f = sqrt(2) * eta;
	for(int i=0;i<Nx; i++){
		double x = i*dx - x0;
		psi0[i] = 2*f/cosh(eta * x);
	}
}

void lin_step (cmplx* const psi0, const double dt, const double dx, const int Nx) {

	const cmplx alpha = cmplx(0.0, -dt/(dx*dx));

	cmplx* o = new cmplx[Nx];
	cmplx* r = new cmplx[Nx];

	for (int i=0; i<Nx; i++) {
		r[i] = psi0[i];
	}

	o[0] = -alpha/(1.0+2.0*alpha);
	r[0] = r[0]/(1.0+2.0*alpha);

	for (int i=1; i<Nx; i++) {
		o[i] = -alpha/(1.0+2.0*alpha + alpha*o[i-1]);
		r[i] = (r[i] + alpha*r[i-1])/(1.0+2.0*alpha + alpha*o[i-1]);
	}

	psi0[Nx-1] = r[Nx-1];

	for (int i=2; i<=Nx; i++) {
		psi0[Nx-i] = r[Nx-i] - o[Nx-i]*psi0[Nx-i+1];
	}

	delete[] r;
	delete[] o;

}

void nonlin_step (cmplx* const psi0, const double dt, const int Nx) {

	const cmplx expt = cmplx (0.0, dt);

	for (int i=0; i<Nx; i++) {
		psi0[i] = psi0[i] * exp(expt*pow(abs(psi0[i]),2));
	}
}
