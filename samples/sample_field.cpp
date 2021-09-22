#include "elec_magn_field.h"

double f(double x, double y = 0.0, double z = 0.0) {
	return double(sin(2 * PI * (x - A) / (B - A)));
}

double f_real_solution(double t, double x, double y = 0.0, double z = 0.0) {
	return double(sin(2 * PI * (x - A - c * t) / (B - A)));
}

int main()
{
	double t = 0.001;
	elec_magn_field Field(f);
	elec_magn_field F = Field.field_integrate(t);
	std::cout << " real solution:" << std::endl;
	for (unsigned i = 0; i < m; ++i)
		 {
			std::cout << "i= " << i << std::endl;
			std::cout << std::endl;
			std::cout << f_real_solution(t, A + dx * i);
			std::cout << std::endl;
		}
	return 0;
}
