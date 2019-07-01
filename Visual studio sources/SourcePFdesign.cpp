#define BZZ_COMPILER 3
#include "BzzMath.hpp"
#include "cmath"

void SysDiff(BzzVector &y, double t, BzzVector &dy);
FILE *res;

//DATA

// PHYSICAL PROPERTIES
double Tp0 = 303;										// K
double Tg0 = 403;										// K
double P = 1;											// atm
double CpL = 1;											// kcal/kg/K
double CpG = 0.25;										// kcal/kg/K
double DHev = 540;										// kcal/kg
double rhoL = 1000;										// kg/m^3
double muG = 2.3*1.0e-5;								// kg/m/s
double kG = 8.*1.0e-6;									// kcal/m/s/K
double Diff = 1.8*1.e-5;								// m^2/s

// Antoine equation
double A = 18.3036;
double B = 3816.44;
double C = -46.13;

// Process stream
double Qmilk = 1750. / 3600;							// kg/h -> kg/s
double fat = 4.76 / 100;								// kg/kg
double Qfat = Qmilk*fat;								// kg/s
double Ql = Qmilk - Qfat;								// kg/s
double Win = Ql / Qfat;									// kg/kg
double Wout = 0.005;									// kg/kg
double dp = 2 * 1.0e-4;									// m
double mp0 = rhoL*3.14 / 6 * pow(dp, 3);				// kg
double mdry = mp0 / (1 + Win);							// kg
double mout = mdry*(1 + Wout);							// kg
double nAvp = Qmilk / mp0;								// drops/s

// Gas stream
double Gdry = 20;										// kg/s
double D = 5.5;											// m


void main(void)
{	
	double t = 0.;
	double delta = 0.1;
	BzzVector y0(6, mp0, 0, Tp0, Tg0, 0, 0);
	BzzVector y, dy;
	
	BzzOdeStiff o(y0, t, SysDiff);
	
	res = fopen("res.txt", "w");
	
	

	for (int i = 0; i <= 10000; i++) //int stands for integer
	{
		t = double(i*1e-3);
		y = o(t); // Discretization of time
		
		double rhoG = P * 101325 / 8.314 / y[4] * 29 * 1e-3;
		double vg = (Gdry+y[2])*4 / (rhoG * (BZZ_PI_GRECO*pow(D, 2)));
		
		printf("\n%d\t%e\t%e", i, y[1], y[2]); //we can either use square brackets or round brackets %d stands for integer
		fprintf(res, "\n%f\t%e\t%e\t%e\t%e\t%e\t%e\t%e", i*1e-3, y[1], y[2]+Gdry,y[3],y[4],y[5],y[6],vg);
		
		if (y[1] <= mout)
		{
			break;
		}
	}
	
	fclose(res);
}

void SysDiff(BzzVector &y, double t, BzzVector &dy)
{
	// mp Gvap Tp Tg vp z

	// Pw = 29 / 18 * y[2] * P / Gdry
	// P°= exp(A-B/(y[3]+C))/760						// mmHg--->atm
	// Sp = 3.14*pow(pow((y[1]/rhoL)*6/3.14,0.3333),2)

	// Mass and heat transfer
	double rhoG = P * 101325 / 8.314 / y[4] * 29 * 1e-3;	// kg/m^3
	double vg = (Gdry+y[2])*4 / (rhoG * (BZZ_PI_GRECO*pow(D, 2)));
	double Re = rhoG*y[5]*dp / muG;							// Reynolds
	double Pr = muG*CpG / kG;								// Prandtl
	double Sc = muG / rhoG / Diff;							// Schmidt
	double Nu = 2. + 0.4*pow(Re, 0.5)*pow(Pr, 0.333);		// Nusselt
	double Sh = 2. + 0.4*pow(Re, 0.5)*pow(Sc, 0.333);		// Sherwood
	double h = Nu*kG / dp;									// kcal/m^2/s/K
	double Kc = Sh*Diff / dp;								// m/s
	double Kpw = Kc / 0.0821 / Tg0 * 18;					// kg/s/atm
	double dp = pow((y[1] / rhoL) * 6 / 3.14, 0.3333);		//m
	double Sp = 3.14*pow(dp, 2);							//m^2
	double P0 = exp(A - B / (y[3] + C)) / 760;				//atm
	double Pwa = 29 / 18 * y[2] * P / Gdry;

	dy[1] = Kpw*Sp*(Pwa-P0);
	dy[2] = -Kpw*Sp*(Pwa-P0)*nAvp;
	dy[3] = (h*Sp*(y[4]-y[3])+Kpw*Sp*(Pwa-P0)*DHev)/y[1]/CpL;
	dy[4] = (-h*Sp*(y[4] - y[3]))*nAvp / (y[2] + Gdry) / CpG;
	dy[5] = ((rhoL - rhoG)*9.81*3.14*pow(dp, 3) / 6 - 3 * muG*y[5] * 3.14*dp - y[5] * Kpw*Sp*(Pwa - P0)) / y[1];
	dy[6] = y[5] + vg;
}