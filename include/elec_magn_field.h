#pragma once

#include <cmath>
#include <iostream>
#include <vector>
#include <fstream>

const double PI = 3.14;
const double c = 3e+10; // скорость света
const int n = 150;

struct field_characteristics {
	std::vector<std::vector<std::vector<double>>> x;
	std::vector<std::vector<std::vector<double>>> y;
	std::vector<std::vector<std::vector<double>>> z;
};

struct estimated_area {
	const double a = - n * c / 2; // начало сетки
	const double b = n * c / 2; // конец сетки
};

struct field_steps {
	const double dx = c;
	const double dy = c;
	const double dz = c;
	const double dt = 1;
};

class elec_magn_field {

public:
	field_characteristics E, B, J;

	estimated_area area;

	field_steps steps;

	const double T = 16;
	const double Tx = 16 * c;
	const double Ty = 16 * c;
	const double Tz = 16 * c;

	// кол-во узлов в сетке
	const int nx = int((area.b - area.a) / steps.dx);
	const int ny = int((area.b - area.a) / steps.dy);
	const int nz = int((area.b - area.a) / steps.dz);

	elec_magn_field() {

		estimated_area area{};
		field_steps steps{};

		std::vector<double> v(nz);
		std::vector<std::vector<double>> v2;
		v2.assign(ny, v);

		E.x.assign(nx, v2);
		E.y.assign(nx, v2);
		E.z.assign(nx, v2);

		B.x.assign(nx, v2);
		B.y.assign(nx, v2);
		B.z.assign(nx, v2);

		J.x.assign(nx, v2);
		J.y.assign(nx, v2);
		J.z.assign(nx, v2);

		for (int i = 0; i < nx; ++i)
			for (int j = 0; j < ny; ++j)
			    for (int k = 0; k < nz; ++k) {
				    E.x[i][j][k] = 0.0;
				    E.y[i][j][k] = 0.0;
				    E.z[i][j][k] = 0.0;

				    B.x[i][j][k] = 0.0;
				    B.y[i][j][k] = 0.0;
				    B.z[i][j][k] = 0.0;

				    J.x[i][j][k] = 0.0;
				    J.y[i][j][k] = 0.0;
				    J.z[i][j][k] = 0.0;
			}

	}


	elec_magn_field(double(*pf)(double, double, double)) {

		estimated_area area{};
		field_steps steps{};

		std::vector<double> v(nz);
		std::vector<std::vector<double>> v2;
		v2.assign(ny, v);

		E.x.assign(nx, v2);
		E.y.assign(nx, v2);
		E.z.assign(nx, v2);

		B.x.assign(nx, v2);
		B.y.assign(nx, v2);
		B.z.assign(nx, v2);

		J.x.assign(nx, v2);
		J.y.assign(nx, v2);
		J.z.assign(nx, v2);

		for (int i = 0; i < nx; ++i) {
			for (int j = 0; j < ny; ++j)
				for (int k = 0; k < nz; ++k) {
					E.x[i][j][k] = 0.0;
					E.y[i][j][k] = pf(steps.dx * i + area.a, 0.0, 0.0);
					E.z[i][j][k] = 0.0;

					B.x[i][j][k] = 0.0;
					B.y[i][j][k] = 0.0;
					B.z[i][j][k] = pf(steps.dx * i + area.a + steps.dx * 0.5, 0.0, 0.0);

					J.x[i][j][k] = 0.0;
					J.y[i][j][k] = 0.0;
					J.z[i][j][k] = 0.0;

					
				}
			}	
	}

	elec_magn_field& get_boundary_conditions_E() {

		// i == 0
		for (int j = 0; j < ny; ++j) {
			if (j != 0) {
				for (int k = 0; k < nz; ++k) {
					E.z[0][j][k] = E.z[0][j][k] - 4.0 * PI * steps.dt * J.z[0][j][k] + c * steps.dt * ((B.y[0][j][k] - B.y[nx - 1][j][k]) / (2.0 * steps.dx) - (B.x[0][j][k] - B.x[0][j - 1][k]) / (2.0 * steps.dy));
					if (k != 0) {
						E.x[0][j][k] = E.x[0][j][k] - 4.0 * PI * steps.dt * J.x[0][j][k] + c * steps.dt * ((B.z[0][j][k] - B.z[0][j - 1][k]) / (2.0 * steps.dy) - (B.y[0][j][k] - B.y[0][j][k - 1]) / (2.0 * steps.dz));
						E.y[0][j][k] = E.y[0][j][k] - 4.0 * PI * steps.dt * J.y[0][j][k] + c * steps.dt * ((B.x[0][j][k] - B.x[0][j][k - 1]) / (2.0 * steps.dz) - (B.z[0][j][k] - B.z[nx - 1][j][k]) / (2.0 * steps.dx));
					}
					if (k == 0) {
						E.x[0][j][k] = E.x[0][j][k] - 4.0 * PI * steps.dt * J.x[0][j][k] + c * steps.dt * ((B.z[0][j][k] - B.z[0][j - 1][k]) / (2.0 * steps.dy) - (B.y[0][j][k] - B.y[0][j][nz - 1]) / (2.0 * steps.dz));
						E.y[0][j][k] = E.y[0][j][k] - 4.0 * PI * steps.dt * J.y[0][j][k] + c * steps.dt * ((B.x[0][j][k] - B.x[0][j][nz - 1]) / (2.0 * steps.dz) - (B.z[0][j][k] - B.z[nx - 1][j][k]) / (2.0 * steps.dx));
					}
				}
			}
			if (j == 0) {
				for (int k = 0; k < nz; ++k) {
					E.z[0][j][k] = E.z[0][j][k] - 4.0 * PI * steps.dt * J.z[0][j][k] + c * steps.dt * ((B.y[0][j][k] - B.y[nx - 1][j][k]) / (2.0 * steps.dx) - (B.x[0][j][k] - B.x[0][ny - 1][k]) / (2.0 * steps.dy));
					if (k != 0) {
						E.x[0][j][k] = E.x[0][j][k] - 4.0 * PI * steps.dt * J.x[0][j][k] + c * steps.dt * ((B.z[0][j][k] - B.z[0][ny - 1][k]) / (2.0 * steps.dy) - (B.y[0][j][k] - B.y[0][j][k - 1]) / (2.0 * steps.dz));
						E.y[0][j][k] = E.y[0][j][k] - 4.0 * PI * steps.dt * J.y[0][j][k] + c * steps.dt * ((B.x[0][j][k] - B.x[0][j][k - 1]) / (2.0 * steps.dz) - (B.z[0][j][k] - B.z[nx - 1][j][k]) / (2.0 * steps.dx));
					}
					if (k == 0) {
						E.x[0][j][k] = E.x[0][j][k] - 4.0 * PI * steps.dt * J.x[0][j][k] + c * steps.dt * ((B.z[0][j][k] - B.z[0][ny - 1][k]) / (2.0 * steps.dy) - (B.y[0][j][k] - B.y[0][j][nz - 1]) / (2.0 * steps.dz));
						E.y[0][j][k] = E.y[0][j][k] - 4.0 * PI * steps.dt * J.y[0][j][k] + c * steps.dt * ((B.x[0][j][k] - B.x[0][j][nz - 1]) / (2.0 * steps.dz) - (B.z[0][j][k] - B.z[nx - 1][j][k]) / (2.0 * steps.dx));
					}
				}
			}
		}

		// j == 0
		for (int i = 0; i < nx; ++i) {
			if (i != 0) {
				for (int k = 0; k < nz; ++k) {
					E.z[i][0][k] = E.z[i][0][k] - 4.0 * PI * steps.dt * J.z[i][0][k] + c * steps.dt * ((B.y[i][0][k] - B.y[i - 1][0][k]) / (2.0 * steps.dx) - (B.x[i][0][k] - B.x[i][ny - 1][k]) / (2.0 * steps.dy));
					if (k != 0) {
						E.x[i][0][k] = E.x[i][0][k] - 4.0 * PI * steps.dt * J.x[i][0][k] + c * steps.dt * ((B.z[i][0][k] - B.z[i][ny - 1][k]) / (2.0 * steps.dy) - (B.y[i][0][k] - B.y[i][0][k - 1]) / (2.0 * steps.dz));
						E.y[i][0][k] = E.y[i][0][k] - 4.0 * PI * steps.dt * J.y[i][0][k] + c * steps.dt * ((B.x[i][0][k] - B.x[i][0][k - 1]) / (2.0 * steps.dz) - (B.z[i][0][k] - B.z[i - 1][0][k]) / (2.0 * steps.dx));
					}
					if (k == 0) {
						E.x[i][0][k] = E.x[i][0][k] - 4.0 * PI * steps.dt * J.x[i][0][k] + c * steps.dt * ((B.z[i][0][k] - B.z[i][ny - 1][k]) / (2.0 * steps.dy) - (B.y[i][0][k] - B.y[i][0][nz - 1]) / (2.0 * steps.dz));
						E.y[i][0][k] = E.y[i][0][k] - 4.0 * PI * steps.dt * J.y[i][0][k] + c * steps.dt * ((B.x[i][0][k] - B.x[i][0][nz - 1]) / (2.0 * steps.dz) - (B.z[i][0][k] - B.z[i - 1][0][k]) / (2.0 * steps.dx));
					}
				}
			}
			if (i == 0) {
				for (int k = 0; k < nz; ++k) {
					E.z[i][0][k] = E.z[i][0][k] - 4.0 * PI * steps.dt * J.z[i][0][k] + c * steps.dt * ((B.y[i][0][k] - B.y[nx - 1][0][k]) / (2.0 * steps.dx) - (B.x[i][0][k] - B.x[i][ny - 1][k]) / (2.0 * steps.dy));
					if (k != 0) {
						E.x[i][0][k] = E.x[i][0][k] - 4.0 * PI * steps.dt * J.x[i][0][k] + c * steps.dt * ((B.z[i][0][k] - B.z[i][ny - 1][k]) / (2.0 * steps.dy) - (B.y[i][0][k] - B.y[i][0][k - 1]) / (2.0 * steps.dz));
						E.y[i][0][k] = E.y[i][0][k] - 4.0 * PI * steps.dt * J.y[i][0][k] + c * steps.dt * ((B.x[i][0][k] - B.x[i][0][k - 1]) / (2.0 * steps.dz) - (B.z[i][0][k] - B.z[nx - 1][0][k]) / (2.0 * steps.dx));
					}
					if (k == 0) {
						E.x[i][0][k] = E.x[i][0][k] - 4.0 * PI * steps.dt * J.x[i][0][k] + c * steps.dt * ((B.z[i][0][k] - B.z[i][ny - 1][k]) / (2.0 * steps.dy) - (B.y[i][0][k] - B.y[i][0][nz - 1]) / (2.0 * steps.dz));
						E.y[i][0][k] = E.y[i][0][k] - 4.0 * PI * steps.dt * J.y[i][0][k] + c * steps.dt * ((B.x[i][0][k] - B.x[i][0][nz - 1]) / (2.0 * steps.dz) - (B.z[i][0][k] - B.z[nx - 1][0][k]) / (2.0 * steps.dx));
					}
				}
			}
		}

		// k == 0
		for (int j = 0; j < ny; ++j) {
			if (j != 0) {
				for (int i = 0; i < nx; ++i) {
					E.x[i][j][0] = E.x[i][j][0] - 4.0 * PI * steps.dt * J.x[i][j][0] + c * steps.dt * ((B.z[i][j][0] - B.z[i][j - 1][0]) / (2.0 * steps.dy) - (B.y[i][j][0] - B.y[i][j][nz - 1]) / (2.0 * steps.dz));
					if (i != 0) {
						E.y[i][j][0] = E.y[i][j][0] - 4.0 * PI * steps.dt * J.y[i][j][0] + c * steps.dt * ((B.x[i][j][0] - B.x[i][j][nz - 1]) / (2.0 * steps.dz) - (B.z[i][j][0] - B.z[i - 1][j][0]) / (2.0 * steps.dx));
						E.z[i][j][0] = E.z[i][j][0] - 4.0 * PI * steps.dt * J.z[i][j][0] + c * steps.dt * ((B.y[i][j][0] - B.y[i - 1][j][0]) / (2.0 * steps.dx) - (B.x[i][j][0] - B.x[i][j - 1][0]) / (2.0 * steps.dy));
					}
					if (i == 0) {
						E.y[i][j][0] = E.y[i][j][0] - 4.0 * PI * steps.dt * J.y[i][j][0] + c * steps.dt * ((B.x[i][j][0] - B.x[i][j][nz - 1]) / (2.0 * steps.dz) - (B.z[i][j][0] - B.z[nx - 1][j][0]) / (2.0 * steps.dx));
						E.z[i][j][0] = E.z[i][j][0] - 4.0 * PI * steps.dt * J.z[i][j][0] + c * steps.dt * ((B.y[i][j][0] - B.y[nx - 1][j][0]) / (2.0 * steps.dx) - (B.x[i][j][0] - B.x[i][j - 1][0]) / (2.0 * steps.dy));
					}
				}
			}

			if (j == 0) {
				for (int i = 0; i < nx; ++i) {
					E.x[i][j][0] = E.x[i][j][0] - 4.0 * PI * steps.dt * J.x[i][j][0] + c * steps.dt * ((B.z[i][j][0] - B.z[i][ny - 1][0]) / (2.0 * steps.dy) - (B.y[i][j][0] - B.y[i][j][nz - 1]) / (2.0 * steps.dz));
					if (i != 0) {
						E.y[i][j][0] = E.y[i][j][0] - 4.0 * PI * steps.dt * J.y[i][j][0] + c * steps.dt * ((B.x[i][j][0] - B.x[i][j][nz - 1]) / (2.0 * steps.dz) - (B.z[i][j][0] - B.z[i - 1][j][0]) / (2.0 * steps.dx));
						E.z[i][j][0] = E.z[i][j][0] - 4.0 * PI * steps.dt * J.z[i][j][0] + c * steps.dt * ((B.y[i][j][0] - B.y[i - 1][j][0]) / (2.0 * steps.dx) - (B.x[i][j][0] - B.x[i][ny - 1][0]) / (2.0 * steps.dy));
					}
					if (i == 0) {
						E.y[i][j][0] = E.y[i][j][0] - 4.0 * PI * steps.dt * J.y[i][j][0] + c * steps.dt * ((B.x[i][j][0] - B.x[i][j][nz - 1]) / (2.0 * steps.dz) - (B.z[i][j][0] - B.z[nx - 1][j][0]) / (2.0 * steps.dx));
						E.z[i][j][0] = E.z[i][j][0] - 4.0 * PI * steps.dt * J.z[i][j][0] + c * steps.dt * ((B.y[i][j][0] - B.y[nx - 1][j][0]) / (2.0 * steps.dx) - (B.x[i][j][0] - B.x[i][ny - 1][0]) / (2.0 * steps.dy));
					}
				}
			}
		}

		/*for (int i = 0; i < nx; i++)
		     { // j == 0
			E.x[i][0][k] = E.x[i][0][k] - 4.0 * PI * steps.dt * J.x[i][0][k] + c * steps.dt * (B.z[i][1][k] - B.z[i][ny - 1][k]) / (2.0 * steps.dy);
			if ((i != 0) && (i != nx - 1)) {
				E.y[i][0] = E.y[i][0] - 4.0 * PI * steps.dt * J.y[i][0] - c * steps.dt * (B.z[i + 1][0] - B.z[i - 1][0]) / (2.0 * steps.dx);
				E.z[i][0] = E.z[i][0] - 4.0 * PI * steps.dt * J.z[i][0] + c * steps.dt * ((B.y[i + 1][0] - B.y[i - 1][0]) / (2.0 * steps.dx) - \
					(B.x[i][1] - B.x[i][ny - 1]) / (2.0 * steps.dy));
			}
			if (i == 0) {
				E.y[i][0] = E.y[i][0] - 4.0 * PI * steps.dt * J.y[i][0] - c * steps.dt * (B.z[i + 1][0] - B.z[nx - 1][0]) / (2.0 * steps.dx);
				E.z[i][0] = E.z[i][0] - 4.0 * PI * steps.dt * J.z[i][0] + c * steps.dt * ((B.y[i + 1][0] - B.y[nx - 1][0]) / (2.0 * steps.dx) - \
					(B.x[i][1] - B.x[i][ny - 1]) / (2.0 * steps.dy));
			}
			if (i == nx - 1) {
				E.y[i][0] = E.y[i][0] - 4.0 * PI * steps.dt * J.y[i][0] - c * steps.dt * (B.z[0][0] - B.z[i - 1][0]) / (2.0 * steps.dx);
				E.z[i][0] = E.z[i][0] - 4.0 * PI * steps.dt * J.z[i][0] + c * steps.dt * ((B.y[0][0] - B.y[i - 1][0]) / (2.0 * steps.dx) - \
					(B.x[i][1] - B.x[i][ny - 1]) / (2.0 * steps.dy));
			}
		}

		for (int i = 0; i < nx; i++) { // j == ny-1
			E.x[i][ny - 1] = E.x[i][ny - 1] - 4.0 * PI * steps.dt * J.x[i][ny - 1] + c * steps.dt * (B.z[i][0] - B.z[i][ny - 2]) / (2.0 * steps.dy);
			if ((i != 0) && (i != nx - 1)) {
				E.y[i][ny - 1] = E.y[i][ny - 1] - 4.0 * PI * steps.dt * J.y[i][ny - 1] - c * steps.dt * (B.z[i + 1][ny - 1] - B.z[i - 1][ny - 1]) / (2.0 * steps.dx);
				E.z[i][ny - 1] = E.z[i][ny - 1] - 4.0 * PI * steps.dt * J.z[i][ny - 1] + c * steps.dt * ((B.y[i + 1][ny - 1] - B.y[i - 1][ny - 1]) / (2.0 * steps.dx) - \
					(B.x[i][0] - B.x[i][ny - 2]) / (2.0 * steps.dy));
			}
			if (i == 0) {
				E.y[i][ny - 1] = E.y[i][ny - 1] - 4.0 * PI * steps.dt * J.y[i][ny - 1] - c * steps.dt * (B.z[i + 1][ny - 1] - B.z[nx - 1][ny - 1]) / (2.0 * steps.dx);
				E.z[i][ny - 1] = E.z[i][ny - 1] - 4.0 * PI * steps.dt * J.z[i][ny - 1] + c * steps.dt * ((B.y[i + 1][ny - 1] - B.y[nx - 1][ny - 1]) / (2.0 * steps.dx) - \
					(B.x[i][0] - B.x[i][ny - 2]) / (2.0 * steps.dy));
			}
			if (i == nx - 1) {
				E.y[i][ny - 1] = E.y[i][ny - 1] - 4.0 * PI * steps.dt * J.y[i][ny - 1] - c * steps.dt * (B.z[0][ny - 1] - B.z[i - 1][ny - 1]) / (2.0 * steps.dx);
				E.z[i][ny - 1] = E.z[i][ny - 1] - 4.0 * PI * steps.dt * J.z[i][ny - 1] + c * steps.dt * ((B.y[0][ny - 1] - B.y[i - 1][ny - 1]) / (2.0 * steps.dx) - \
					(B.x[i][0] - B.x[i][ny - 2]) / (2.0 * steps.dy));
			}
		}

		for (int j = 1; j < ny - 1; j++) { // i == 0
			E.x[0][j] = E.x[0][j] - 4.0 * PI * steps.dt * J.x[0][j] + c * steps.dt * (B.z[0][j + 1] - B.z[0][j - 1]) / (2.0 * steps.dy);
			E.y[0][j] = E.y[0][j] - 4.0 * PI * steps.dt * J.y[0][j] - c * steps.dt * (B.z[1][j] - B.z[nx - 1][j]) / (2.0 * steps.dx);
			E.z[0][j] = E.z[0][j] - 4.0 * PI * steps.dt * J.z[0][j] + c * steps.dt * ((B.y[1][j] - B.y[nx - 1][j]) / (2.0 * steps.dx) - \
				(B.x[0][j + 1] - B.x[0][j - 1]) / (2.0 * steps.dy));
		}

		for (int j = 1; j < ny - 1; j++) { // i == nx-1
			E.x[nx - 1][j] = E.x[nx - 1][j] - 4.0 * PI * steps.dt * J.x[nx - 1][j] + c * steps.dt * (B.z[nx - 1][j + 1] - B.z[nx - 1][j - 1]) / (2.0 * steps.dy);
			E.y[nx - 1][j] = E.y[nx - 1][j] - 4.0 * PI * steps.dt * J.y[nx - 1][j] - c * steps.dt * (B.z[0][j] - B.z[nx - 2][j]) / (2.0 * steps.dx);
			E.z[nx - 1][j] = E.z[nx - 1][j] - 4.0 * PI * steps.dt * J.z[nx - 1][j] + c * steps.dt * ((B.y[0][j] - B.y[nx - 2][j]) / (2.0 * steps.dx) - \
				(B.x[nx - 1][j + 1] - B.x[nx - 1][j - 1]) / (2.0 * steps.dy));
		}*/

		return (*this);
	}

	elec_magn_field& get_boundary_conditions_B() {

		// i == nx - 1
		for (int j = 0; j < ny; ++j) {
			if (j != ny - 1) {
				for (int k = 0; k < nz; ++k) {
					B.z[nx - 1][j][k] = B.z[nx - 1][j][k] + c * steps.dt * ((E.x[nx - 1][j + 1][k] - E.x[nx - 1][j][k]) / (2.0 * steps.dy) - (E.y[0][j][k] - E.y[nx - 1][j][k]) / (2.0 * steps.dx));
					if (k != nz - 1) {
						B.x[nx - 1][j][k] = B.x[nx - 1][j][k] + c * steps.dt * ((E.y[nx - 1][j][k + 1] - E.y[nx - 1][j][k]) / (2.0 * steps.dz) - (E.z[nx - 1][j + 1][k] - E.z[nx - 1][j][k]) / (2.0 * steps.dy));
						B.y[nx - 1][j][k] = B.y[nx - 1][j][k] + c * steps.dt * ((E.z[0][j][k] - E.z[nx - 1][j][k]) / (2.0 * steps.dx) - (E.x[nx - 1][j][k + 1] - E.x[nx - 1][j][k]) / (2.0 * steps.dz));
					}
					if (k == nz - 1) {
						B.x[nx - 1][j][k] = B.x[nx - 1][j][k] + c * steps.dt * ((E.y[nx - 1][j][0] - E.y[nx - 1][j][k]) / (2.0 * steps.dz) - (E.z[nx - 1][j + 1][k] - E.z[nx - 1][j][k]) / (2.0 * steps.dy));
						B.y[nx - 1][j][k] = B.y[nx - 1][j][k] + c * steps.dt * ((E.z[0][j][k] - E.z[nx - 1][j][k]) / (2.0 * steps.dx) - (E.x[nx - 1][j][0] - E.x[nx - 1][j][k]) / (2.0 * steps.dz));
					}
				}
			}
			if (j == ny - 1) {
				for (int k = 0; k < nz; ++k) {
					B.z[nx - 1][j][k] = B.z[nx - 1][j][k] + c * steps.dt * ((E.x[nx - 1][0][k] - E.x[nx - 1][j][k]) / (2.0 * steps.dy) - (E.y[0][j][k] - E.y[nx - 1][j][k]) / (2.0 * steps.dx));
					if (k != nz - 1) {
						B.x[nx - 1][j][k] = B.x[nx - 1][j][k] + c * steps.dt * ((E.y[nx - 1][j][k + 1] - E.y[nx - 1][j][k]) / (2.0 * steps.dz) - (E.z[nx - 1][0][k] - E.z[nx - 1][j][k]) / (2.0 * steps.dy));
						B.y[nx - 1][j][k] = B.y[nx - 1][j][k] + c * steps.dt * ((E.z[0][j][k] - E.z[nx - 1][j][k]) / (2.0 * steps.dx) - (E.x[nx - 1][j][k + 1] - E.x[nx - 1][j][k]) / (2.0 * steps.dz));
					}
					if (k == nz - 1) {
						B.x[nx - 1][j][k] = B.x[nx - 1][j][k] + c * steps.dt * ((E.y[nx - 1][j][0] - E.y[nx - 1][j][k]) / (2.0 * steps.dz) - (E.z[nx - 1][0][k] - E.z[nx - 1][j][k]) / (2.0 * steps.dy));
						B.y[nx - 1][j][k] = B.y[nx - 1][j][k] + c * steps.dt * ((E.z[0][j][k] - E.z[nx - 1][j][k]) / (2.0 * steps.dx) - (E.x[nx - 1][j][0] - E.x[nx - 1][j][k]) / (2.0 * steps.dz));
					}
				}
			}
		}

		// j == ny - 1
		for (int i = 0; i < nx; ++i) {
			if (i != nx - 1) {
				for (int k = 0; k < nz; ++k) {
					B.z[i][ny - 1][k] = B.z[i][ny - 1][k] + c * steps.dt * ((E.x[i][0][k] - E.x[i][ny - 1][k]) / (2.0 * steps.dy) - (E.y[i + 1][ny - 1][k] - E.y[i][ny - 1][k]) / (2.0 * steps.dx));
					if (k != nz - 1) {
						B.x[i][ny - 1][k] = B.x[i][ny - 1][k] + c * steps.dt * ((E.y[i][ny - 1][k + 1] - E.y[i][ny - 1][k]) / (2.0 * steps.dz) - (E.z[i][0][k] - E.z[i][ny - 1][k]) / (2.0 * steps.dy));
						B.y[i][ny - 1][k] = B.y[i][ny - 1][k] + c * steps.dt * ((E.z[i + 1][ny - 1][k] - E.z[i][ny - 1][k]) / (2.0 * steps.dx) - (E.x[i][ny - 1][k + 1] - E.x[i][ny - 1][k]) / (2.0 * steps.dz));
					}
					if (k == nz - 1) {
						B.x[i][ny - 1][k] = B.x[i][ny - 1][k] + c * steps.dt * ((E.y[i][ny - 1][0] - E.y[i][ny - 1][k]) / (2.0 * steps.dz) - (E.z[i][0][k] - E.z[i][ny - 1][k]) / (2.0 * steps.dy));
						B.y[i][ny - 1][k] = B.y[i][ny - 1][k] + c * steps.dt * ((E.z[i + 1][ny - 1][k] - E.z[i][ny - 1][k]) / (2.0 * steps.dx) - (E.x[i][ny - 1][0] - E.x[i][ny - 1][k]) / (2.0 * steps.dz));
					}
				}
			}
			if (i == nx - 1) {
				for (int k = 0; k < nz; ++k) {
					B.z[i][ny - 1][k] = B.z[i][ny - 1][k] + c * steps.dt * ((E.x[i][0][k] - E.x[i][ny - 1][k]) / (2.0 * steps.dy) - (E.y[0][ny - 1][k] - E.y[i][ny - 1][k]) / (2.0 * steps.dx));
					if (k != nz - 1) {
						B.x[i][ny - 1][k] = B.x[i][ny - 1][k] + c * steps.dt * ((E.y[i][ny - 1][k + 1] - E.y[i][ny - 1][k]) / (2.0 * steps.dz) - (E.z[i][0][k] - E.z[i][ny - 1][k]) / (2.0 * steps.dy));
						B.y[i][ny - 1][k] = B.y[i][ny - 1][k] + c * steps.dt * ((E.z[0][ny - 1][k] - E.z[i][ny - 1][k]) / (2.0 * steps.dx) - (E.x[i][ny - 1][k + 1] - E.x[i][ny - 1][k]) / (2.0 * steps.dz));
					}
					if (k == nz - 1) {
						B.x[i][ny - 1][k] = B.x[i][ny - 1][k] + c * steps.dt * ((E.y[i][ny - 1][0] - E.y[i][ny - 1][k]) / (2.0 * steps.dz) - (E.z[i][0][k] - E.z[i][ny - 1][k]) / (2.0 * steps.dy));
						B.y[i][ny - 1][k] = B.y[i][ny - 1][k] + c * steps.dt * ((E.z[0][ny - 1][k] - E.z[i][ny - 1][k]) / (2.0 * steps.dx) - (E.x[i][ny - 1][0] - E.x[i][ny - 1][k]) / (2.0 * steps.dz));
					}
				}
			}
		}


		// k == nz - 1
		for (int i = 0; i < nx; ++i) {
			if (i != nx - 1) {
				for (int j = 0; j < ny; ++j) {
					B.y[i][j][nz - 1] = B.y[i][j][nz - 1] + c * steps.dt * ((E.z[i + 1][j][nz - 1] - E.z[i][j][nz - 1]) / (2.0 * steps.dx) - (E.x[i][j][0] - E.x[i][j][nz - 1]) / (2.0 * steps.dz));
					if (j != ny - 1) {
						B.x[i][j][nz - 1] = B.x[i][j][nz - 1] + c * steps.dt * ((E.y[i][j][0] - E.y[i][j][nz - 1]) / (2.0 * steps.dz) - (E.z[i][j + 1][nz - 1] - E.z[i][j][nz - 1]) / (2.0 * steps.dy));
						B.z[i][j][nz - 1] = B.z[i][j][nz - 1] + c * steps.dt * ((E.x[i][j + 1][nz - 1] - E.x[i][j][nz - 1]) / (2.0 * steps.dy) - (E.y[i + 1][j][nz - 1] - E.y[i][j][nz - 1]) / (2.0 * steps.dx));
					}
					if (j == ny - 1) {
						B.x[i][j][nz - 1] = B.x[i][j][nz - 1] + c * steps.dt * ((E.y[i][j][0] - E.y[i][j][nz - 1]) / (2.0 * steps.dz) - (E.z[i][0][nz - 1] - E.z[i][j][nz - 1]) / (2.0 * steps.dy));
						B.z[i][j][nz - 1] = B.z[i][j][nz - 1] + c * steps.dt * ((E.x[i][0][nz - 1] - E.x[i][j][nz - 1]) / (2.0 * steps.dy) - (E.y[i + 1][j][nz - 1] - E.y[i][j][nz - 1]) / (2.0 * steps.dx));
					}
				}
			}
			if (i == nx - 1) {
				for (int j = 0; j < ny; ++j) {
					B.y[i][j][nz - 1] = B.y[i][j][nz - 1] + c * steps.dt * ((E.z[0][j][nz - 1] - E.z[i][j][nz - 1]) / (2.0 * steps.dx) - (E.x[i][j][0] - E.x[i][j][nz - 1]) / (2.0 * steps.dz));
					if (j != ny - 1) {
						B.x[i][j][nz - 1] = B.x[i][j][nz - 1] + c * steps.dt * ((E.y[i][j][0] - E.y[i][j][nz - 1]) / (2.0 * steps.dz) - (E.z[i][j + 1][nz - 1] - E.z[i][j][nz - 1]) / (2.0 * steps.dy));
						B.z[i][j][nz - 1] = B.z[i][j][nz - 1] + c * steps.dt * ((E.x[i][j + 1][nz - 1] - E.x[i][j][nz - 1]) / (2.0 * steps.dy) - (E.y[0][j][nz - 1] - E.y[i][j][nz - 1]) / (2.0 * steps.dx));
					}
					if (j == ny - 1) {
						B.x[i][j][nz - 1] = B.x[i][j][nz - 1] + c * steps.dt * ((E.y[i][j][0] - E.y[i][j][nz - 1]) / (2.0 * steps.dz) - (E.z[i][0][nz - 1] - E.z[i][j][nz - 1]) / (2.0 * steps.dy));
						B.z[i][j][nz - 1] = B.z[i][j][nz - 1] + c * steps.dt * ((E.x[i][0][nz - 1] - E.x[i][j][nz - 1]) / (2.0 * steps.dy) - (E.y[0][j][nz - 1] - E.y[i][j][nz - 1]) / (2.0 * steps.dx));
					}
				}
			}
		}

		/*for (int i = 0; i < nx; i++) { // j == 0
			B.x[i][0] = B.x[i][0] - c * steps.dt * (E.z[i][1] - E.z[i][ny - 1]) / (2.0 * steps.dy);
			if ((i != 0) && (i != nx - 1)) {
				B.y[i][0] = B.y[i][0] + c * steps.dt * (E.z[i + 1][0] - E.z[i - 1][0]) / (2.0 * steps.dx);
				B.z[i][0] = B.z[i][0] - c * steps.dt * ((E.y[i + 1][0] - E.y[i - 1][0]) / (2.0 * steps.dx) - (E.x[i][1] - E.x[i][ny - 1]) / (2.0 * steps.dy));
			}
			if (i == 0) {
				B.y[i][0] = B.y[i][0] + c * steps.dt * (E.z[i + 1][0] - E.z[nx - 1][0]) / (2.0 * steps.dx);
				B.z[i][0] = B.z[i][0] - c * steps.dt * ((E.y[i + 1][0] - E.y[nx - 1][0]) / (2.0 * steps.dx) - (E.x[i][1] - E.x[i][ny - 1]) / (2.0 * steps.dy));
			}
			if (i == nx - 1) {
				B.y[i][0] = B.y[i][0] + c * steps.dt * (E.z[0][0] - E.z[i - 1][0]) / (2.0 * steps.dx);
				B.z[i][0] = B.z[i][0] - c * steps.dt * ((E.y[0][0] - E.y[i - 1][0]) / (2.0 * steps.dx) - (E.x[i][1] - E.x[i][ny - 1]) / (2.0 * steps.dy));
			}
		}

		for (int i = 0; i < nx; i++) { // j == ny - 1
			B.x[i][ny - 1] = B.x[i][ny - 1] - c * steps.dt * (E.z[i][0] - E.z[i][ny - 2]) / (2.0 * steps.dy);
			if ((i != 0) && (i != nx - 1)) {
				B.y[i][ny - 1] = B.y[i][ny - 1] + c * steps.dt * (E.z[i + 1][ny - 1] - E.z[i - 1][ny - 1]) / (2.0 * steps.dx);
				B.z[i][ny - 1] = B.z[i][ny - 1] - c * steps.dt * ((E.y[i + 1][ny - 1] - E.y[i - 1][ny - 1]) / (2.0 * steps.dx) - (E.x[i][0] - E.x[i][ny - 2]) / (2.0 * steps.dy));
			}
			if (i == 0) {
				B.y[i][ny - 1] = B.y[i][ny - 1] + c * steps.dt * (E.z[i + 1][ny - 1] - E.z[nx - 1][ny - 1]) / (2.0 * steps.dx);
				B.z[i][ny - 1] = B.z[i][ny - 1] - c * steps.dt * ((E.y[i + 1][ny - 1] - E.y[nx - 1][ny - 1]) / (2.0 * steps.dx) - (E.x[i][0] - E.x[i][ny - 2]) / (2.0 * steps.dy));
			}
			if (i == nx - 1) {
				B.y[i][ny - 1] = B.y[i][ny - 1] + c * steps.dt * (E.z[0][ny - 1] - E.z[i - 1][ny - 1]) / (2.0 * steps.dx);
				B.z[i][ny - 1] = B.z[i][ny - 1] - c * steps.dt * ((E.y[0][ny - 1] - E.y[i - 1][ny - 1]) / (2.0 * steps.dx) - (E.x[i][0] - E.x[i][ny - 2]) / (2.0 * steps.dy));
			}
		}

		for (int j = 1; j < ny - 1; j++) {  // i == 0
			B.x[0][j] = B.x[0][j] - c * steps.dt * (E.z[0][j + 1] - E.z[0][j - 1]) / (2.0 * steps.dy);
			B.y[0][j] = B.y[0][j] + c * steps.dt * (E.z[1][j] - E.z[nx - 1][j]) / (2.0 * steps.dx);
			B.z[0][j] = B.z[0][j] - c * steps.dt * ((E.y[1][j] - E.y[nx - 1][j]) / (2.0 * steps.dx) - (E.x[0][j + 1] - E.x[0][j - 1]) / (2.0 * steps.dy));
		}

		for (int j = 1; j < ny - 1; j++) { // i == nx-1
			B.x[nx - 1][j] = B.x[nx - 1][j] - c * steps.dt * (E.z[nx - 1][j + 1] - E.z[nx - 1][j - 1]) / (2.0 * steps.dy);
			B.y[nx - 1][j] = B.y[nx - 1][j] + c * steps.dt * (E.z[0][j] - E.z[nx - 2][j]) / (2.0 * steps.dx);
			B.z[nx - 1][j] = B.z[nx - 1][j] - c * steps.dt * ((E.y[0][j] - E.y[nx - 2][j]) / (2.0 * steps.dx) - (E.x[nx - 1][j + 1] - E.x[nx - 1][j - 1]) / (2.0 * steps.dy));
		}*/

		return (*this);
	}


   

	elec_magn_field& FDTD() {
		// обновление B

		for (int i = 0; i < nx - 1; i++)
			for (int j = 0; j < ny - 1; j++)
				for (int k = 0; k < nz - 1; k++) {
					B.x[i][j][k] = B.x[i][j][k] + c * steps.dt * ((E.y[i][j][k + 1] - E.y[i][j][k]) / (2.0 * steps.dz) - (E.z[i][j + 1][k] - E.z[i][j][k]) / (2.0 * steps.dy));
				    B.y[i][j][k] = B.y[i][j][k] + c * steps.dt * ((E.z[i + 1][j][k] - E.z[i][j][k]) / (2.0 * steps.dx) - (E.x[i][j][k + 1] - E.x[i][j][k]) / (2.0 * steps.dz));
				    B.z[i][j][k] = B.z[i][j][k] + c * steps.dt * ((E.x[i][j + 1][k] - E.x[i][j][k]) / (2.0 * steps.dy) - (E.y[i + 1][j][k] - E.y[i][j][k]) / (2.0 * steps.dx));
			}

		get_boundary_conditions_B();

		// обновление E

		for (int i = 1; i < nx; i++)
			for (int j = 1; j < ny; j++)
			    for (int k = 1; k < nz; k++) {
				    E.x[i][j][k] = E.x[i][j][k] - 4.0 * PI * steps.dt * J.x[i][j][k] + c * steps.dt * ((B.z[i][j][k] - B.z[i][j - 1][k]) / (2.0 * steps.dy) - (B.y[i][j][k] - B.y[i][j][k - 1]) / (2.0 * steps.dz));
				    E.y[i][j][k] = E.y[i][j][k] - 4.0 * PI * steps.dt * J.y[i][j][k] + c * steps.dt * ((B.x[i][j][k] - B.x[i][j][k - 1]) / (2.0 * steps.dz) - (B.z[i][j][k] - B.z[i - 1][j][k]) / (2.0 * steps.dx));
				    E.z[i][j][k] = E.z[i][j][k] - 4.0 * PI * steps.dt * J.z[i][j][k] + c * steps.dt * ((B.y[i][j][k] - B.y[i - 1][j][k]) / (2.0 * steps.dx) - (B.x[i][j][k] - B.x[i][j - 1][k]) / (2.0 * steps.dy));
			}

		get_boundary_conditions_E();

		return (*this);
	}


	elec_magn_field& modify_Jz(int m) {
		double t = steps.dt * (m - 0.5);
		double x, y, z;
		for (int i = 0; i < nx; ++i)
			for (int j = 0; j < ny; ++j)
			    for (int k = 0; k < nz; ++k) {
				    x = area.a + i * steps.dx;
				    y = area.a + j * steps.dy;
				    z = area.a + k * steps.dz;
					if ((x >= -Tx / 4.0) && (x <= Tx / 4.0) && (y >= -Ty / 4.0) && (y <= Ty / 4.0) && (z >= -Tz / 4.0) && (z <= Tz / 4.0)) {
						J.z[i][j][k] = sin(2 * PI * t / T) * pow(cos(2 * PI * x / Tx), 2) * pow(cos(2 * PI * y / Ty), 2) * pow(cos(2 * PI * z / Tz), 2);
						//std::cout << "Jz[" << i << "][" << j << "]=" << J.z[i][j][k] << std::endl;
					}
					
			}
		return (*this);
	}

	elec_magn_field& set_zero_Jz() {
		for (int i = 0; i < nx; ++i)
			for (int j = 0; j < ny; ++j)
				for (int k = 0; k < nz; ++k)
		            J.z[i][j][k] = 0.0;
		return (*this);
	}




};

