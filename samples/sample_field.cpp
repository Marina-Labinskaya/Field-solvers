#include "elec_magn_field.h"

double f(double x, double y = 0.0, double z = 0.0) {
	estimated_area Area;
	return double(sin(2 * PI * (x - Area.A) / (Area.B - Area.A)));
}

double f_real_solution(double t, double x, double y = 0.0, double z = 0.0) {
	estimated_area Area;
	return double(sin(2 * PI * (x - Area.A - c * t) / (Area.B - Area.A)));
}


double f2( double y, double x = 0.0, double z = 0.0) {
	estimated_area Area;
	return double(sin(2 * PI * (y - Area.C) / (Area.D - Area.C)));
}

double f_real_solution2(double t,  double y, double x = 0.0, double z = 0.0) {
	estimated_area Area;
	return double(sin(2 * PI * (y - Area.C - c * t) / (Area.D - Area.C)));
}


int main()
{
	/*int test_x = 4;
	int test_y = 3;
	std::ofstream fout("Ey_field.csv");
	elec_magn_field F(f);
	for (double t = F.steps.dt; t < F.T; t += F.steps.dt) {
		elec_magn_field F2 = F.FDTD();
		fout << t << ";" << F2.E.y[test_x][test_y] << ";" << f_real_solution(t, test_x * F2.steps.dx + F2.area.A) << std::endl;
	}
	fout.close();
	std::ofstream fout3("Ey_x.csv");
	elec_magn_field F1(f);
	for (double t = F1.steps.dt; t < F1.T; t += F1.steps.dt)
		F1.FDTD();
	for (int i = 0; i < F1.E.y.size(); ++i)
		fout3 << i * F1.steps.dx << ";" << F1.E.y[i][test_y] << ";" << f_real_solution(F1.steps.dt, i* F1.steps.dx + F1.area.A) << std::endl;
	fout3.close();*/

	/*std::ofstream fout2("Ex_y.csv");
	elec_magn_field F3(f2);
	for (double t = F3.steps.dt; t < F3.T; t += F3.steps.dt)
		F3.FDTD();
	std::cout << F3.E.x.size();
	for (int i = 0; i < F3.E.x[test_x].size(); ++i)
		fout2 << i * F3.steps.dy << ";" << F3.E.x[test_x][i] << ";" << f_real_solution2(F3.steps.dt, i* F3.steps.dy + F3.area.C) << std::endl;
	fout2.close();*/

	// сферическая волна

	std::ofstream fout("spherical_wave.txt");
	const double t1 = 30.0;

	field_steps st{ c, c, c };
	estimated_area ar{ (- n * c / 2), (n * c / 2 ), (- n * c / 2) , (n * c / 2) };
	elec_magn_field F(st, ar);
	int i = 0;
	int k = int(F.T / F.steps.dt);
	for (double t = 0; t < t1; t += F.steps.dt) {
		if (i <= k)
			F.modify_Jz_center(i);
		else {
			F.set_zero_Jz_center();
		}
		F.FDTD();
		++i;
	}
	for (int i = 0; i < F.nx; ++i) {
		for (int j = 0; j < F.ny; ++j)
			fout << F.E.z[i][j] << " ";
		fout << std::endl;
	}
	
	fout.close();

	std::ofstream fout2("x.txt");
	for (int i = 0; i < F.nx; ++i)
	   for (int j = 0; j < F.ny; ++j) {
			fout2 << F.area.A + i * F.steps.dx << " ";
	}

	fout2.close();

	std::ofstream fout3("y.txt");
	for (int i = 0; i < F.nx; ++i)
		for (int j = 0; j < F.ny; ++j) {
			fout3 << F.area.A + j * F.steps.dx << " ";
		}

	fout3.close();

	return 0;
}
