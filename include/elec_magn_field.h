#pragma once

#include <cmath>
#include <iostream>
#include <vector>
#include <fstream>

const double PI = 3.14;
const double c = 3e+10; // скорость света
const int n = 200;

struct field_characteristics {
	std::vector<std::vector<double>> x;
	std::vector<std::vector<double>> y;
	std::vector<std::vector<double>> z;
};

struct estimated_area {
	const double A, B, C, D;
	estimated_area(double a = - n * c / 2, double b = n * c / 2, double c = - n * c / 2, double d = n * c / 2) : A{ a }, B{ b }, C{ c }, D{ d } {}
};

struct field_steps {
	const double dx;
	const double dy;
	const double dz;
	const double dt;
	field_steps(double _dx = c, double _dy = c,  double _dz = c, double _dt = 1 ) : dx{ _dx }, dy{ _dy },  dz{ _dz }, dt{ _dt } {}
};

class elec_magn_field {

public:
	field_characteristics E, B, J;

	estimated_area area;

	field_steps steps;

	const double T = 16;
	const double Tx = 16 * c;
	const double Ty = 16 * c;

	// кол-во узлов в сетке
	const int nx = int((area.B - area.A) / steps.dx);
	const int ny = int((area.D - area.C) / steps.dy);
	const int nz;

	elec_magn_field() : nz{ 1 } {

		estimated_area area{};
		field_steps steps{};

		E.x.assign(nx, std::vector<double>(ny));
		E.y.assign(nx, std::vector<double>(ny));
		E.z.assign(nx, std::vector<double>(ny));

		B.x.assign(nx, std::vector<double>(ny));
		B.y.assign(nx, std::vector<double>(ny));
		B.z.assign(nx, std::vector<double>(ny));

		J.x.assign(nx, std::vector<double>(ny));
		J.y.assign(nx, std::vector<double>(ny));
		J.z.assign(nx, std::vector<double>(ny));

		for (int i = 0; i < nx; ++i)
			for (int j = 0; j < ny; ++j) {
				E.x[i][j] = 0.0;
				E.y[i][j] = 0.0;
				E.z[i][j] = 0.0;

				B.x[i][j] = 0.0;
				B.y[i][j] = 0.0;
				B.z[i][j] = 0.0;

				J.x[i][j] = 0.0;
				J.y[i][j] = 0.0;
				J.z[i][j] = 0.0;
			}

	}

	elec_magn_field(field_steps st, estimated_area& ar, int _nz = 1) : nz{ _nz} {

		field_steps steps(st.dx, st.dy, st.dz, st.dt);
		estimated_area area(ar.A, ar.B, ar.C, ar.D);

		E.x.assign(nx, std::vector<double>(ny));
		E.y.assign(nx, std::vector<double>(ny));
		E.z.assign(nx, std::vector<double>(ny));

		B.x.assign(nx, std::vector<double>(ny));
		B.y.assign(nx, std::vector<double>(ny));
		B.z.assign(nx, std::vector<double>(ny));

		J.x.assign(nx, std::vector<double>(ny));
		J.y.assign(nx, std::vector<double>(ny));
		J.z.assign(nx, std::vector<double>(ny));

		for (int i = 0; i < nx; ++i)
			for (int j = 0; j < ny; ++j) {
				E.x[i][j] = 0.0;
				E.y[i][j] = 0.0;
				E.z[i][j] = 0.0;

				B.x[i][j] = 0.0;
				B.y[i][j] = 0.0;
				B.z[i][j] = 0.0;

				J.x[i][j] = 0.0;
				J.y[i][j] = 0.0;
				J.z[i][j] = 0.0;
			}
	}

	elec_magn_field(double(*pf)(double, double, double)) : nz{ 1 } {

		estimated_area area{};
		field_steps steps{};

		E.x.assign(nx, std::vector<double>(ny));
		E.y.assign(nx, std::vector<double>(ny));
		E.z.assign(nx, std::vector<double>(ny));

		B.x.assign(nx, std::vector<double>(ny));
		B.y.assign(nx, std::vector<double>(ny));
		B.z.assign(nx, std::vector<double>(ny));

		J.x.assign(nx, std::vector<double>(ny));
		J.y.assign(nx, std::vector<double>(ny));
		J.z.assign(nx, std::vector<double>(ny));

		for (int i = 0; i < nx; ++i)
			for (int j = 0; j < ny; ++j) {

				E.x[i][j] = 0.0;
				//E.x[i][j] = pf(steps.dy * j + area.C, 0.0, 0.0);
				E.y[i][j] = pf(steps.dx * i + area.A, 0.0, 0.0);
				//E.y[i][j] = 0.0;
				E.z[i][j] = 0.0;

				B.x[i][j] = 0.0;
				B.y[i][j] = 0.0;
				B.z[i][j] = E.y[i][j];
				//B.z[i][j] = -E.x[i][j];

				J.x[i][j] = 0.0;
				J.y[i][j] = 0.0;
				J.z[i][j] = 0.0;

				//std::cout << "Ex[" << i << "][" << j << "]=" << E.x[i][j] << std::endl;
				//std::cout << "Ey[" << i << "][" << j << "]=" << E.y[i][j] << std::endl;
				//std::cout << "Bz[" << i << "][" << j << "]=" << B.z[i][j] << std::endl;
				//std::cout << std::endl;
			}
	}

	elec_magn_field& get_boundary_conditions_E() {

		for (int i = 0; i < nx; i++) { // j == 0
			E.x[i][0] = E.x[i][0] - 4.0 * PI * steps.dt * J.x[i][0] + c * steps.dt * (B.z[i][1] - B.z[i][ny - 1]) / (2.0 * steps.dy);
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
		}

		return (*this);
	}

	elec_magn_field& get_boundary_conditions_B() {

		for (int i = 0; i < nx; i++) { // j == 0
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
		}

		return (*this);
	}

    elec_magn_field& FDTD() {
		// интегрирование

		// обновление E

		   for (int i = 1; i < nx - 1; i++)
			   for (int j = 1; j < ny - 1; j++) {
				   E.x[i][j] = E.x[i][j] - 4.0 * PI * steps.dt * J.x[i][j] + c * steps.dt * (B.z[i][j + 1] - B.z[i][j - 1]) / (2.0 * steps.dy);
				   E.y[i][j] = E.y[i][j] - 4.0 * PI * steps.dt * J.y[i][j] - c * steps.dt * (B.z[i + 1][j] - B.z[i - 1][j]) / (2.0 * steps.dx);
				   E.z[i][j] = E.z[i][j] - 4.0 * PI * steps.dt * J.z[i][j] + c * steps.dt * ((B.y[i + 1][j] - B.y[i - 1][j]) / (2.0 * steps.dx) - \
					   (B.x[i][j + 1] - B.x[i][j - 1]) / (2.0 * steps.dy));
			   }

		   get_boundary_conditions_E();

		// обновление B

		   for (int i = 1; i < nx - 1; i++)
			   for (int j = 1; j < ny - 1; j++) {
				   B.x[i][j] = B.x[i][j] - c * steps.dt * (E.z[i][j + 1] - E.z[i][j - 1]) / (2.0 * steps.dy);
				   B.y[i][j] = B.y[i][j] + c * steps.dt * (E.z[i + 1][j] - E.z[i - 1][j]) / (2.0 * steps.dx);
				   B.z[i][j] = B.z[i][j] - c * steps.dt * ((E.y[i + 1][j] - E.y[i - 1][j]) / (2.0 * steps.dx) - (E.x[i][j + 1] - E.x[i][j - 1]) / (2.0 * steps.dy));
			   }

		   get_boundary_conditions_B();

		return (*this);
	}

	elec_magn_field& modify_Jz_center(int i) {
		double w1 = 4 * PI;
		double w2 = PI / 8;
		//int k = int(T / steps.dt);
		int k = int((2 * PI) / (w2 * steps.dt));
		double t = steps.dt * i;
		//J.z[nx / 2][ny / 2] = sin(2 * PI * t / T);
	    J.z[nx / 2][ny / 2] = sin(w1 * t) * pow (sin (w2 * t), 2);
		return (*this);
	}

	elec_magn_field& modify_Jz(int i) {
		double t = steps.dt * i;
		double x, y;
		for (int k = 0; k < nx; ++k)
			for (int j = 0; j < ny; ++j) {
				x = area.A + k * steps.dx;
				y = area.C + j * steps.dy;
				if ((x >= -Tx / 4.0) && (x <= Tx / 4.0) && (y >= -Ty / 4.0) && (y <= Ty / 4.0))
					J.z[k][j] = sin(2 * PI * t / T) * pow(cos(2 * PI * x / Tx), 2) * pow(cos(2 * PI * y / Ty), 2);
			}
		return (*this);
	}

	elec_magn_field& set_zero_Jz_center() {
		J.z[nx / 2][ny / 2] = 0.0;
		return (*this);
	}

	elec_magn_field& set_zero_Jz() {
		for (int i = 0; i < nx; ++i)
			for (int j = 0; j < ny; ++j)
		        J.z[i][j] = 0.0;
		return (*this);
	}




};

