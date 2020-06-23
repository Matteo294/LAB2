#include <iostream>
#include <cmath>

using namespace std;

#define R 17.5e-3 / 2 

double f(double x, double d);

int main(){

    srand(time(NULL));

    const int npoints = 1e6;
    const float d = 0.100; // 10 mm
    const double mu0 = 4 * M_PI * 1e-7;
    const int N1 = 28;
    const int N2 = 28;
    const double sigma1 = M_PI * pow(R, 2);
    const double sigma2 = sigma1;
    
    double I = 0;
    double h = (double) 2*M_PI / npoints;
    double x = 0;
    for (int i = 0; i < npoints; i++){
        I += h * f(x, d);
        x += h;
    }
    cout << "Value: " << I * mu0 * N1 * N2 * pow(R, 2) / 2 << endl;
    cout << "Valore in approssimazione a dipolo: " << mu0/(4*M_PI) * 2 * N1 * N2 * sigma1 * sigma2 / pow(d, 3) << endl;

    return 0;
}

double f(double x, double d){
    return (double) cos(x) / sqrt(2 * pow(R, 2) * (1 - cos(x)) + pow(d, 2));
}