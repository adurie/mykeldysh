#include <iostream>
#include <complex>
#include <cmath>
#include <fstream>
#include <string>
#include <gsl/gsl_integration.h>
#include <eigen3/Eigen/Dense>
#include <eigen3/Eigen/Eigenvalues>
#include <vector>

using namespace Eigen;
using namespace std;
typedef complex<double> dcomp;

double fermi(double arg, double Ef){
	const double kT = 0.01;
	return 1./(1.+exp((arg-Ef)/kT));
}

vector<double> f(double x, double z, double a, double E, double Ef, int N, double theta) {
// ...NM|ins|FM(0)|NM(n)|FM(theta)...
	dcomp i;
	double F = cos(x*a)+cos(z*a);

	/* const double V = 0.2; */
	const double V = 0.0;
	const double t = 0.5;
	const double nab = -0.6795;
	Matrix2cd T, NM1, NM2, FM1, FM2, ins, I, S;
	//same hopping throughout for simplicity
	T << t,0.,0.,t;

	NM1 << -2.8, 0., 0., -2.8;
	NM2 = NM1;
	ins << 0.1, 0., 0., 0.1;
	FM1 << -1.6295 + nab, 0., 0., -1.6295 - nab;
	FM2 = FM1; 
	I << 1.,0.,0.,1.;
	S << cos(theta/2.),sin(theta/2.),-sin(theta/2.),cos(theta/2.);
	FM2 = S.inverse()*FM2*S;
	i=-1;
	i=sqrt(i);
	const dcomp im =(1e-3)*i;
//apply the bias to the LHS
	NM1 = NM1 + I*V;
	ins = ins + I*V;

	Matrix2cd OMV2=(E+im)*I-FM2-2.*T*F;

	Matrix4cd X2,O2;
	X2 << 	0,	0,	1/t,	0,
		0,	0,	0,	1/t,
		-t,	0,	OMV2(0,0)/t,OMV2(0,1)/t,
		0,	-t,	OMV2(1,0)/t,OMV2(1,1)/t;
	ComplexEigenSolver<Matrix4cd> ces2;
	ces2.compute(X2);
	O2 = ces2.eigenvectors();
	Matrix2cd b2 = O2.topRightCorner(2,2);
	Matrix2cd d2 = O2.bottomRightCorner(2,2);
	Matrix2cd GR = b2*d2.inverse();

	Matrix2cd OM = (E+im)*I - 2.*T*F;

	Matrix2cd OMV1=(E+im)*I-NM1-2.*T*F;

	double om=E-2.*t*F;

	Matrix4cd X,O,Oinv;
	X << 	0,	0,	1/t,	0,
		0,	0,	0,	1/t,
		-t,	0,	OMV1(0,0)/t,OMV1(0,1)/t,
		0,	-t,	OMV1(1,0)/t,OMV1(1,1)/t;
	ComplexEigenSolver<Matrix4cd> ces;
	ces.compute(X);
	O = ces.eigenvectors();
	Matrix2cd b = O.topRightCorner(2,2);
	Matrix2cd d = O.bottomRightCorner(2,2);
	Matrix2cd GL = b*d.inverse();

	Matrix2cd Pauli, xPau, yPau, zPau;
	xPau << 0.,1.,1.,0.;
	yPau << 0.,-i,i,0.;
	zPau << 1.,0.,0.,-1.;

	Pauli = yPau;
	double spincurrent;
	Matrix2cd A,B,TOT,GM;
//lim is thickness of layer 2
	const int lim = 8;
//build thickness of layer 2 to lim layers
	for (int it=0; it < lim; ++it){
		ins = ins - I*(V*it/(lim*1.-1));
		GL = (OM - ins -T*GL*T).inverse();
	}
//lim2 is thickness of layer 3
	const int lim2 = 16;
//build thickness of layer 3 to lim2 layers
	for (int it=0; it < lim2; ++it){

		GL = (OM - FM1 -T*GL*T).inverse();
	}
//adlayer layer 2 from layer 1 to spacer thickness, N
	vector<double> result;
	result.reserve(N);
	for (int it=0; it < N-1; ++it){
		A = (I-GR*T.adjoint()*GL*T).inverse();
		B = (I-GR.adjoint()*T.adjoint()*GL.adjoint()*T).inverse();
		TOT = (B-A)*Pauli;
		spincurrent = (1./(4.*M_PI))*real(TOT.trace())*(fermi(E,Ef+V)+fermi(E,Ef));
		result.emplace_back(spincurrent);
		GL = (OM - NM2 -T*GL*T).inverse();
	}
	return result;
}

vector<double> int_theta(double x, double z, double a, double E, double Ef, int N){
	vector<double> result;
	vector<double> integrate;
	result.reserve(N);
	integrate.reserve(N);
	for (int i = 0; i < N; i++)
		result[i] = 0.;
	double theta;

	const int n = 10;
	for (int k=0; k<n+1; k++) {
		theta = k*M_PI/n;
		integrate = f(x, z, a, E, Ef, N, theta);
		for (int i = 0; i < N; i++){
			if ((k==0)||(k==n))
				result[i] += (0.5/n)*integrate[i];
			else 
				result[i] += (1./n)*integrate[i];
		}
	}	
	return result;
}

vector<double> int_energy(double x, double z, double a, double Ef, int N){
	vector<double> result;
	vector<double> integrate;
	result.reserve(N);
	integrate.reserve(N);
	for (int i = 0; i < N; i++)
		result[i] = 0.;

	double E;
	double end = 0.5;
	double start = -2.;

	const int n = 10000;
	for (int k=0; k<n+1; k++) {
		E = start + k*(end-start)/n;
		integrate = int_theta(x, z, a, E, Ef, N);
		for (int i = 0; i < N; i++){
			if ((k==0)||(k==n))
				result[i] += (0.5/n)*integrate[i];
			else 
				result[i] += (1./n)*integrate[i];
		}
	}	
	return result;
}

vector<double> int_kpoints(double a, double Ef, int N){
	vector<double> result;
	vector<double> integrate;
	result.reserve(N);
	integrate.reserve(N);
	double x, z;
	for (int i = 0; i < N; i++)
		result[i] = 0.;

	int n = 700;
	double factor = 8./(n*n);
	for (int k = 0; k!=n+1; k++){
		if (k%2!=0){
			x = M_PI*k/n;
			for (int l = 0; l!=k+1; l++){
				if (l%2!=0){
					z = M_PI*l/n;
					/* integrate = int_theta(x, z, 1, Ef, Ef, N); */
					integrate = int_energy(x, z, 1, Ef, N);
					for (int i = 0; i < N; i++){
						if ((k==1) && (l==1))
							result[i] += factor*0.5*integrate[i];
						else if (k==l)
							result[i] += factor*0.5*integrate[i];
						else
							result[i] += factor*integrate[i];
					}
				}
			}
		}
		if (k%20 == 0)
			cout<<"     "<<k/(n*1.+1)*100.<<"% completed"<<endl;
	}
	return result;
}

int main() 
{
	//number of atomic planes
	// plot output of spincurrent against energy
	string Mydata;
	ofstream Myfile;	
	Mydata = "tmp.txt";
	/* Mydata = "eq_sc_fixed_k_no_V.txt"; */
	Myfile.open( Mydata.c_str(),ios::trunc );
	const double Ef = -0.1;

	// number of spacer layers
	int N = 30;
	vector<double> answer;
	answer.reserve(N);
	/* answer = int_theta(0, 0, 1, Ef, Ef, i); */
	/* answer = int_energy(0, 0, 1, Ef, N); */
	answer = int_kpoints(1, Ef, N);
	/* answer = f(0, 0, 1, Ef, Ef, i, 0); */
	for (int i = 0; i < N; i++){
		Myfile<<i<<" "<<answer[i]<<endl;
	}
	return 0;
}
