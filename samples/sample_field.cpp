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

	double t1;
	setlocale(LC_ALL, "Rus");
	std::cout << "Введите время t1: ";
	std::cin >> t1;
	std::cout << std::endl;
	std::ofstream fout("spherical_wave2.txt");
	

	elec_magn_field<> f;

	int i = 1;
	int k = int(f.T / f.steps.dt);

	for (double t = f.steps.dt; t <= t1; t += f.steps.dt) {
		if (i <= k) {
			f.modify_Jz(i);
		}
		else {
			f.set_zero_Jz();
		}
		f.PSTD(t);
		i++;
	}

	for (int i = 0; i < f.nx; ++i) {
		for (int j = 0; j < f.ny; ++j)
			fout << f.E.z(i, j, 0) << " ";
		fout << std::endl;
	}
	

	fout.close();
	std::cout << "Вычисления выполнены";
	return 0;
}
