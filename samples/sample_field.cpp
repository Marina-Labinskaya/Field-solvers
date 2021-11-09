#include "elec_magn_field.h"

double f(double x, double y = 0.0, double z = 0.0) {
	estimated_area Area;
	return double(sin(2 * PI * (x - Area.a) / (Area.b - Area.a)));
}

double f_real_solution(double t, double x, double y = 0.0, double z = 0.0) {
	estimated_area Area;
	return double(sin(2 * PI * (x - Area.a - c * t) / (Area.b - Area.a)));
}


double f2( double y, double x = 0.0, double z = 0.0) {
	estimated_area Area;
	return double(sin(2 * PI * (y - Area.a) / (Area.b - Area.a)));
}

double f_real_solution2(double t,  double y, double x = 0.0, double z = 0.0) {
	estimated_area Area;
	return double(sin(2 * PI * (y - Area.a - c * t) / (Area.b - Area.a)));
}


int main()
{
	// сферическая волна

	std::ofstream fout("spherical_wave2.txt");
	const double t1 = 100.0;
	int kz = 0;

	elec_magn_field F;
	int i = 1;
	int k = int(F.T / F.steps.dt);

	F.start_FDTD();

	for (double t = F.steps.dt; t < t1; t += F.steps.dt) {
		if (i <= k)
			F.modify_Jz(i);
		else {
			F.set_zero_Jz();
		}
		F.FDTD();
		++i;
	}
	for (int i = 0; i < F.nx; ++i) {
		for (int j = 0; j < F.ny; ++j)
			fout << F.E.z[i][j][F.nz / 2] << " ";
		fout << std::endl;
	}
	
	fout.close();

	return 0;
}
