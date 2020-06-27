#include <iostream>
#include <cmath>
#include <fstream>

using namespace std;

#define R 17.5e-3 / 2 
#define DMIN 0.020 // distanza minima 10mm
#define DMAX 0.200 // distanza massima 200mm
#define NPUNTI 500 // Distanze in cui calcolare l'integrale

double f(double x, double d);

int main(){

    srand(time(NULL));

    double dstep = (DMAX - DMIN) / NPUNTI;
    double d = DMIN;
    const int npoints = 1e5; // Punti per l'integrale
    const double mu0 = 4 * M_PI * 1e-7;
    const int N1 = 30;
    const int N2 = 28;
    const double sigma1 = M_PI * pow(R, 2);
    const double sigma2 = sigma1;

    fstream myfile;
    myfile.open("induzione.csv");

    for (int j = 0; j < NPUNTI; j++){
        
        d += dstep;
        
        double I = 0;
        double h = (double) 2*M_PI / npoints;
        double x = 0;
        for (int i = 0; i < npoints; i++){
            I += h * f(x, d);
            x += h;
        }

        double val = I * mu0 * N1 * N2 * pow(R, 2) / 2;
        double approx = mu0/(4*M_PI) * 2 * N1 * N2 * sigma1 * sigma2 / pow(d, 3);
        /*cout << "Valore integrale doppio: " << val << endl;
        cout << "Valore in approssimazione a dipolo: " << approx << endl;
        cout << "\n" << endl;*/

        myfile << d << "," << val << "," << approx << endl;
    }

    myfile.close();
    return 0;
}

double f(double x, double d){
    return (double) cos(x) / sqrt(2 * pow(R, 2) * (1 - cos(x)) + pow(d, 2));
}