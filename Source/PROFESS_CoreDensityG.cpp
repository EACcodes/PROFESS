#include <cstdlib>
#include <string>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <cmath>
#include <cassert>
#include <ctime>

using namespace std;

void Simpson_Integral
(
    const int mesh,
    const double *func,
    const double *rab,
    double &asum
)
{
    /*     simpson's rule integration. On input:
    !      mesh = mhe number of grid points (should be odd)
    !      func(i)= function to be integrated
    !      rab(i) = r(i) * dr(i)/di * di
    !      For the logarithmic grid not including r=0 :
    !      r(i) = r_0*exp((i-1)*dx) ==> rab(i)=r(i)*dx
    !      For the logarithmic grid including r=0 :
    !      r(i) = a(exp((i-1)*dx)-1) ==> rab(i)=(r(i)+a)*dx
    !      Output in asum = \sum_i c_i f(i)*rab(i) = \int_0^\infty f(r) dr
    !      where c_i are alternativaly 2/3, 4/3 except c_1 = c_mesh = 1/3
    */
    //  simpson's rule integrator for function stored on the
    //  radial logarithmic mesh
    //	routine assumes that mesh is an odd number so run check
    if (mesh % 2 == 0)
    {
        cout << "\n error in subroutine simpson ";
        cout << "\n routine assumes mesh is odd but mesh = "
             << mesh << endl;
		exit(0);
    }

    asum = 0.00;
    const double r12 = 1.00 / 12.00;
    double f3 = func [0] * rab [0] * r12;
    for (int i = 1;i < mesh;i += 2)
    {
        const double f1 = f3;
        const double f2 = func [i] * rab [i] * r12;
        f3 = func [i + 1] * rab [i + 1] * r12;
        asum += 4.00 * f1 + 16.00 * f2 + 4.00 * f3;
    }
    return;
}//

int main()
{
	cout << "----------------------------------------------------------------------" << endl;
	cout << " WELCOME. " << endl;
	cout << " This code is used to transform the real space density (1D)" << endl;
	cout << " into G space density (1D), which can then be read by PROFESS" << endl;
	cout << "----------------------------------------------------------------------" << endl;

	cout << " Information for input file " << endl;
	cout << " The input file has the name: INPUT.dat " << endl;
	cout << " The first line in INPUT.dat is the number of mesh points in real space (INTEGER)" << endl;
	cout << " Begin from second line, INPUT.dat should then have two columns, the first one is r (Bohr)," << endl;
	cout << " with the second one is rho(r), be careful it's NOT rho(r)*r*r " << endl;
	cout << " The rho(G) will then be written to OUTPUT.dat" << endl;
	cout << " " << endl;

	double pi = 3.1415926535897932384626;

	// number of G points
	int nG = 30000;
	// delta G
	double deltaG = 0.002;

	cout << " (YOU CAN CHANGE THE G SPCAE PARAMETERS INSIDE THE CODE IF YOU LIKE)" << endl;
	cout << " The parameters for G space " << endl;
	cout << " The number of points are " << nG << endl;
	cout << " The interval of G points is " << deltaG << " Bohr^-1"  << endl; 
	cout << " " << endl;

	int nmesh = -1;

	// open the file
	string file = "INPUT.dat";
	ifstream ifs(file.c_str());
	if(!ifs)
	{
		cout << " Can't find the " << file << endl;
		exit(0);
	}
	else
	{
		cout << " Thank you! We find the input file " << file << endl;
	}
	ifs >> nmesh;
	cout << " nmesh=" << nmesh << endl;
	if(nmesh==0)
	{
		cout << " We need nmesh ! = 0" << endl;
		cout << " Put it in the first line of " << file << endl;
		exit(0);
	}

	double* r = new double[nmesh];
	double* rab = new double[nmesh];
	double* rho1d = new double[nmesh];
	for(int i=0; i<nmesh; ++i)
	{
	   ifs >> r[i] >> rho1d[i];

	   //************************************************************
	   // rho1d is just simply the core density, don't have r^2 term
	   //************************************************************
	   rho1d[i] = rho1d[i] * r[i] * r[i];

	   if(i>0)
	   {
			rab[i] = r[i]-r[i-1];
	   } 
	}



	cout << " " << endl;
	cout << " The final core density is printed to OUTPUT.dat" << endl;
	cout << " You can use this file in PROFESS calculations. " << endl;
	double* fun2 = new double[nmesh];
	// print out the G space form of density.	
	ofstream out("OUTPUT.dat");
	out << nG << " " << deltaG << endl;
	double g0 = 0.0;
	Simpson_Integral(nmesh,rho1d,rab,g0);
	// the first point.
	out << "0 " << g0 << endl;
	cout << " The first point rho(g)=" << g0 << endl;
	for(int ig=1; ig<nG; ++ig)
	{
//		cout << " ig=" << ig << endl;
		double g = ig * deltaG;
		for(int ir=0; ir<nmesh; ++ir)
		{
			double gr = g * r[ir];
			//-------------------------------
			// rho1d is in unit of Bohr^2
			// rho1d has r^2 already
			//-------------------------------
			fun2[ir]=rho1d[ir]*sin(gr)/gr;
		}
		double value;
		// value is density in G space
		Simpson_Integral(nmesh,fun2,rab,value);
		out << g << " " << setprecision(15) << value << endl;
	}
	return 0;
}
