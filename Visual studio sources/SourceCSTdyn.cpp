#define BZZ_COMPILER 3
#include "BzzMath.hpp"
#include "cmath"

void SysDiff(BzzVector &y, double t, BzzVector &dy);
void SysSS(BzzVector &x, BzzVector &f);
FILE *res_SS;
FILE *res;
FILE *res1;
FILE *res2;
// DATA


double Cpg = 1000;				//J/kg/°C
double Cpw = 4185;				//J/kg/K
double Cpwv = 1800;				//J/kg/°C
double Cps = 1657;				//J/kg/°C
double Dws = 6.67*1.e-10;		//m^2/s
double Dwg = 1.8*1.e-5;			//m^2/s
double eps = 0.98;				//-
double rhos = 800;				//kg/m^3
double rhowa = 1000;			//kg/m3
double R=8.314;					//J/mol/K
double P = 101325;				//Pa
double MW_air = 29*1.e-3;		//kg/mol
double MW_wa = 18*1e-3;			//kg/mol
double r = 5e-5;				//m
double a = 3 * (1 - eps) / r;	//1/m
double V = 1.26;				//m3
double Hwv = 2501*1e3;			//J/kg (DHev(0°C) considered constant with temperature change)

// Antoine
double A = 18.3036;
double B = 3816.44;
double C = -46.13;

//Flows
double Gs = 2.5e-4;				//kg/s
double Gg = 1.872e-1;			//kg/s

//Transport phenomena
double Nu = 2;
double Sh = 2;
double kcs = 3 * pow(BZZ_PI_GRECO, 2)*Dws / r;	//W/m/K
double kcg = Sh*Dwg / (2 * r);	//W/m/K

//Initial conditions
double Ts0 = 27;		//°C
double Tg0 = 160;		//°C
double Xs0 = 4;			// kg_w/kg_s
double Xg0 = 0.0196;	//kg_w/kg_g

//Other datas
double rhog = 1;//kg/m3



void main(void)
{	
	//Non linear system for steady state conditions
	BzzVector x0(4, Xs0, Xg0, Ts0, Tg0);
	BzzVector x,f;
	BzzVector Xmin(4, 0, 0, 0);
	BzzVector Xmax(4, 4, 20, 160);
	BzzNonLinearSystem nls(x0, SysSS);
	//nls.SetMinimumConstraints(Xmin);
	//nls.SetMaximumConstraints(Xmax);
	double tolA = 1.e-4, tolR = 1.e-3;
	//nls.SetTolerance(tolA, tolR);

	nls();
	nls.BzzPrint("Results");
	nls.GetSolution(&x,&f);

	res_SS=fopen("ResSS.txt","w");

	printf("Xs=%f\t%f\nXg=%f\t%f\nTs=%f\t%f\nTg=%f\t%f\n", x[1],f[1], x[2],f[2], x[3],f[3], x[4],f[4]);
	fclose(res_SS);
	system("pause");

	//ODE system for dynamic modeling
	
	BzzVector y0(4, Xs0, Xg0, Ts0, Tg0 );
	BzzVector y, dy;
	BzzVector ymin(4, 0, 0, 0, 0);
	double time = 0.;
	
	res = fopen("Res.txt", "w");
	res1 = fopen("Res1.txt", "w");
	res2 = fopen("Res2.txt", "w");

	BzzOdeStiff o(y0, 0., SysDiff);
	o.SetMinimumConstraints(&ymin);
	for (int i = 0; i <= 10000; i++)
	{
		time = double(i)*1;
		y = o(time);
		
		double aa = y[2] * MW_air / MW_wa;									//-
		double pw = aa / (1 + aa)*P;									//Pa
		double kg = 8.4044*1.e-5*(y[4] + 273.15) + 4.63*1.e-5;			//W/m/K
		double mi = 4.25*1.e-8*(y[4] + 273.15) + 5.87*1.e-6;			//kg/m/s
		double hg = Nu*kg / (2 * r);									//W/m^2/K
		double P0 = exp(A - B / (y[3] + 273.15 + C)) / 760 * 101325;	//Pa
		double j = kcg / (R*(y[4] + 273.15))*a*MW_wa*(P0 - pw)*V;		//kg/s
		double q = hg*a*(y[4] - y[3])*V;								//W

		double GMB = (Gs*(Xs0-y[1])+Gg*(Xg0-y[2]))/(Gs*Xs0+Gg*Xg0);		//global mass balance relative error on water 
		

		fprintf(res, "\n%f\t%e\t%e\t%f\t%f", time, y[1], y[2],y[3],y[4]);
		fprintf(res2, "\n%f\t%f", time, GMB);
		fprintf(res1, "\n%f\tpw=%f\tkg=%f\tmi=%f\thg=%f\tP0=%f\tq=%f\tj=%f", time, pw, kg, mi, hg, P0, q, j);		
	}
	
	fclose(res2);
	fclose(res1);
	fclose(res);
	system("pause");
}

void SysDiff(BzzVector &y, double t, BzzVector &dy)
{
	double aa= y[2]*MW_air / MW_wa;									//-
	double pw = aa / (1 + aa)*P;									//Pa
	double kg = 8.4044*1.e-5*(y[4] + 273.15) + 4.63*1.e-5;			//W/m/K
	double mi = 4.25*1.e-8*(y[4] + 273.15) + 5.87*1.e-6;			//kg/m/s
	double hg = Nu*kg / (2 * r);									//W/m^2/K
	double P0 = exp(A - B / (y[3] + 273.15 + C)) / 760 * 101325;	//Pa
	double j = kcg / (R*(y[4] + 273.15))*a*MW_wa*(P0 - pw)*V;		//kg/s
	double q = hg*a*(y[4] - y[3])*V;								//W
	
	
	dy[1] = -j / (rhos*(1 - eps)*V) - Gs*(y[1] - Xs0) / (rhos*(1 - eps)*V);
	dy[2] = j / (rhog*eps*V) - Gg*(y[2] - Xg0) / (rhog*eps*V);
	dy[3] = q / (V*rhos*(1 - eps)*(Cps + Cpw*y[1])) - j*Hwv / (V*rhos*(1 - eps)*(Cps + Cpw*y[1])) - Cpw*y[3] * dy[1] / ((Cps + Cpw*y[1])) - (Gs*((Cps + Cpw*y[1])*y[3] - (Cps + Cpw*Xs0)*Ts0)) / (rhos*(1 - eps)*V*(Cps + Cpw*y[1]));
	dy[4] = -q / (V*rhog*eps*(Cpg + Cpwv*y[2])) + j*Hwv / (eps*rhog*V*(Cpg + Cpwv*y[2])) - (Hwv + Cpwv*y[4])*dy[2] / (Cpg + Cpwv*y[2]) - (Gg*((Cpg*y[4] + (Hwv + Cpwv*y[4])*y[2]) - (Cpg*Tg0 + (Hwv + Cpwv*Tg0)*Xg0))) / (eps*rhog*V*(Cpg + Cpwv*y[2]));// -19.1*5.06*(y[4] - 25) / (eps*rhog*V*(Cpg + Cpwv*y[2]));
	
}
void SysSS(BzzVector &x, BzzVector &f)
{
	double aa = x[2] * MW_air / MW_wa;									//-
	double pw = aa / (1 + aa)*P;										//Pa						
	double kg = 8.4044*1.e-5*(x[4] + 273.15) + 4.63*1.e-5;				//W/m/K
	double mi = 4.25*1.e-8*(x[4] + 273.15) + 5.87*1.e-6;				//kg/m/s
	double hg = Nu*kg / (2 * r);										//W/m^2/K
	double P0 = exp(A - B / (x[3] + 273.15 + C)) / 760 * 101325;		//Pa
	double j = kcg / (R*(x[4] + 273.15))*a*MW_wa*(P0 - pw)*V;			//kg/s
	double q = hg*a*(x[4] - x[3])*V;									//W

	f[1] = -j - Gs*(x[1] - Xs0);
	f[2] = j  - Gg*(x[2] - Xg0);
	f[3] = q - j*Hwv - (Gs*((Cps + Cpw*x[1])*x[3] - (Cps + Cpw*Xs0)*Ts0));
	f[4] = -q + j*Hwv - (Gg*((Cpg*x[4] + (Hwv + Cpwv*x[4])*x[2]) - (Cpg*Tg0 + (Hwv + Cpwv*Tg0)*Xg0)));
}
	

