#include <iostream>
#include <complex>
#include <cmath>

using namespace std;

double integrand(double x, double z){
	return 1;
}

int main(){
	double x, z;
	double integrate, result;
	result = 0.;

	int n = 700;
	double factor = 8./(n*n);
	for (int k = 0; k!=n+1; k++){
		if (k%2!=0){
			x = M_PI*k/n;
			for (int l = 0; l!=k+1; l++){
				if (l%2!=0){
					z = M_PI*l/n;
					integrate = integrand(x, z);
					if ((k==1) && (l==1))
						result += factor*0.5*integrate;
					else if (k==l)
						result += factor*0.5*integrate;
					else
						result += factor*integrate;
				}
			}
		}
		if (k%20 == 0)
			cout<<"     "<<k/(n*1.+1)*100.<<"% completed"<<endl;
	}
	cout<<result<<endl;
	return 0;
}
