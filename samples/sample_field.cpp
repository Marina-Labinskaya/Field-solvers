#include "elec_magn_field.h"

int main()
{
	// сферическая волна

	double t1;
	setlocale(LC_ALL, "Rus");
	std::cout << "Введите время t1: ";
	std::cin >> t1;
	std::cout << std::endl;
	std::ofstream fout("spherical_wave_modify_Jy.txt");
	

	elec_magn_field<> f(0.1, c, c, c, -n * c /2, n * c /2, 0, c, - n * c / 2, n * c / 2);

	int i = 1;
	int _k = int(f.T / f.steps.dt);

	for (double t = 0; t <= t1; t += f.steps.dt) {
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
