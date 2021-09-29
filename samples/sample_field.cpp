#include "elec_magn_field.h"

double f(double x, double y = 0.0, double z = 0.0) {
	return double(sin(2 * PI * (x - A) / (B - A)));
}

double f_real_solution(double t, double x, double y = 0.0, double z = 0.0) {
	return double(sin(2 * PI * (x - A - c * t) / (B - A)));
}

int main()
{
	int test_x = 4;
	int test_y = 3;
	std::ofstream fout("Ey_field.csv");
	elec_magn_field Field(f);
	for (double t = 0; t < T; t += dt) {
		elec_magn_field F = Field.field_integrate(t);
		fout << t << ";" << F.E.Ey[test_x][test_y] << std::endl;
	}
	fout.close();
	std::ofstream fout2("Ey_real_field.csv");
	for (double t = 0; t < T; t += dt) {
		fout2 << t << ";" << f_real_solution(t, test_x * dx + A) << std::endl;
	}
	fout2.close();
	return 0;
}
