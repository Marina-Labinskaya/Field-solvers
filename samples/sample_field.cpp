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
	std::ofstream fout("spherical_wave_modify_Jy.txt");
	

	elec_magn_field<> f;

	int i = 1;
	int _k = int(f.T / f.steps.dt);

	for (double t = f.steps.dt; t <= t1; t += f.steps.dt) {
		if (i <= _k) {
			f.modify_Jy(i);
		}
		else {
			f.set_zero_Jy();
		}
		f.PSTD(t);
		i++;
	}

	/*for (int j = 0; j < f.ny; ++j) {
		for (int k = 0; k < f.nz; ++k)
			fout << f.E.x(0, j, k) << " ";
		fout << std::endl;
	}*/
	
	for (int i = 0; i < f.nx; ++i) {
		for (int k = 0; k < f.nz; ++k)
			fout << f.E.y(i, 0, k) << " ";
		fout << std::endl;
	}

	fout.close();
	std::cout << "Вычисления выполнены";
	return 0;
}
