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
	for (double t = dt; t < T; t += dt) {
		elec_magn_field F = Field.field_integrate();
		fout << t << ";" << F.E.y[test_x][test_y] << ";" << f_real_solution(t, test_x * dx + A) << std::endl;
		

	}
	fout.close();
	std::ofstream fout3("Ey_x.csv");
	elec_magn_field Field2(f);
	elec_magn_field F2;
	for (double t = dt; t < T; t += dt)
		F2 = Field2.field_integrate();
	for (int i = 0; i < F2.E.y.size(); ++i)
		fout3 << i * dx << ";" << F2.E.y[i][test_y] << ";" << f_real_solution(dt, i* dx + A) << std::endl;
	fout3.close();
	/*std::ofstream fout2("Ey_real_field.csv");
	for (double t = 0; t < T; t += dt) {
		fout2 << t << ";" << f_real_solution(t, test_x * dx + A) << std::endl;
	}
	fout2.close();*/
	return 0;
}
