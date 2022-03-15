#pragma once

#include <cmath>
#include <iostream>
#include <vector>
#include <fstream>
#include <complex.h>
#include <fftw3_mkl.h>

const double PI = 3.1415926536;
const double c = 3e+10; // скорость света
const int n = 50;


template <typename T = double>
class Array {
public:
	std::vector<T> data;
	int nx, ny, nz, size;

	Array() {
		nx = 0;
		ny = 0;
		nz = 0;
		size = 0;
	}

	void resize(int _nx, int _ny, int _nz) {
		nx = _nx;
		ny = _ny;
		nz = _nz;
		size = nx * ny * nz;
		data.assign(size, 0.0);
	}

	~Array() {
		data.clear();
	}

	double& operator() (int i, int j, int k) {
		return data[i* ny * nz + j * nz + k];
	}

	double operator() (int i, int j, int k) const {
		return data[i * ny * nz + j * nz + k];
	}

	double& Array::operator[](int i) {
		return data[i];
	}
};

template<typename T = double>
struct field_characteristics {
	Array<T> x;
	Array<T> y;
	Array<T> z;
};

struct estimated_area {
	const double a = - n * c / 2; // начало сетки
	const double b = n * c / 2; // конец сетки
};

struct point {
	double x;
	double y;
	double z;
};

struct field_steps {
	const double dx = c;
	const double dy = c;
	const double dz = c;
	const double dt = 0.1;
};

template <typename T = double>
class elec_magn_field {

public:
	field_characteristics<T> E, B, J;

	estimated_area area;

	field_steps steps;

	const double T = 16;
	const double Tx = 16 * c;
	const double Ty = 16 * c;
	const double Tz = 16 * c;

	// кол-во узлов в сетке
	const int nx = int((area.b - area.a) / steps.dx);
	const int ny = int((area.b - area.a) / steps.dy);
	//const int nz = int((area.b - area.a) / steps.dz);
	const int nz = 1;
	const int size = nx * ny * nz;

	elec_magn_field() {

		estimated_area area{};
		field_steps steps{};

		E.x.resize(nx, ny, nz);
		E.y.resize(nx, ny, nz);
		E.z.resize(nx, ny, nz);

		B.x.resize(nx, ny, nz);
		B.y.resize(nx, ny, nz);
		B.z.resize(nx, ny, nz);

		J.x.resize(nx, ny, nz);
		J.y.resize(nx, ny, nz);
		J.z.resize(nx, ny, nz);
	}

	elec_magn_field& get_boundary_conditions_E(double t) {

		// i == 0
		for (int j = 0; j < ny; ++j) {
			if (j != 0) {
				for (int k = 0; k < nz; ++k) {
					E.z[0][j][k] = E.z[0][j][k] - 4.0 * PI * t * J.z[0][j][k] + c * t * ((B.y[0][j][k] - B.y[nx - 1][j][k]) / (2.0 * steps.dx) - (B.x[0][j][k] - B.x[0][j - 1][k]) / (2.0 * steps.dy));
					if (k != 0) {
						E.x[0][j][k] = E.x[0][j][k] - 4.0 * PI * t * J.x[0][j][k] + c * t * ((B.z[0][j][k] - B.z[0][j - 1][k]) / (2.0 * steps.dy) - (B.y[0][j][k] - B.y[0][j][k - 1]) / (2.0 * steps.dz));
						E.y[0][j][k] = E.y[0][j][k] - 4.0 * PI * t * J.y[0][j][k] + c * t * ((B.x[0][j][k] - B.x[0][j][k - 1]) / (2.0 * steps.dz) - (B.z[0][j][k] - B.z[nx - 1][j][k]) / (2.0 * steps.dx));
					}
					if (k == 0) {
						E.x[0][j][k] = E.x[0][j][k] - 4.0 * PI * t * J.x[0][j][k] + c * t * ((B.z[0][j][k] - B.z[0][j - 1][k]) / (2.0 * steps.dy) - (B.y[0][j][k] - B.y[0][j][nz - 1]) / (2.0 * steps.dz));
						E.y[0][j][k] = E.y[0][j][k] - 4.0 * PI * t * J.y[0][j][k] + c * t * ((B.x[0][j][k] - B.x[0][j][nz - 1]) / (2.0 * steps.dz) - (B.z[0][j][k] - B.z[nx - 1][j][k]) / (2.0 * steps.dx));
					}
				}
			}
			if (j == 0) {
				for (int k = 0; k < nz; ++k) {
					E.z[0][j][k] = E.z[0][j][k] - 4.0 * PI * t * J.z[0][j][k] + c * t * ((B.y[0][j][k] - B.y[nx - 1][j][k]) / (2.0 * steps.dx) - (B.x[0][j][k] - B.x[0][ny - 1][k]) / (2.0 * steps.dy));
					if (k != 0) {
						E.x[0][j][k] = E.x[0][j][k] - 4.0 * PI * t * J.x[0][j][k] + c * t * ((B.z[0][j][k] - B.z[0][ny - 1][k]) / (2.0 * steps.dy) - (B.y[0][j][k] - B.y[0][j][k - 1]) / (2.0 * steps.dz));
						E.y[0][j][k] = E.y[0][j][k] - 4.0 * PI * t * J.y[0][j][k] + c * t * ((B.x[0][j][k] - B.x[0][j][k - 1]) / (2.0 * steps.dz) - (B.z[0][j][k] - B.z[nx - 1][j][k]) / (2.0 * steps.dx));
					}
					if (k == 0) {
						E.x[0][j][k] = E.x[0][j][k] - 4.0 * PI * t * J.x[0][j][k] + c * t * ((B.z[0][j][k] - B.z[0][ny - 1][k]) / (2.0 * steps.dy) - (B.y[0][j][k] - B.y[0][j][nz - 1]) / (2.0 * steps.dz));
						E.y[0][j][k] = E.y[0][j][k] - 4.0 * PI * t * J.y[0][j][k] + c * t * ((B.x[0][j][k] - B.x[0][j][nz - 1]) / (2.0 * steps.dz) - (B.z[0][j][k] - B.z[nx - 1][j][k]) / (2.0 * steps.dx));
					}
				}
			}
		}

		// j == 0
		for (int i = 0; i < nx; ++i) {
			if (i != 0) {
				for (int k = 0; k < nz; ++k) {
					E.z[i][0][k] = E.z[i][0][k] - 4.0 * PI * t * J.z[i][0][k] + c * t * ((B.y[i][0][k] - B.y[i - 1][0][k]) / (2.0 * steps.dx) - (B.x[i][0][k] - B.x[i][ny - 1][k]) / (2.0 * steps.dy));
					if (k != 0) {
						E.x[i][0][k] = E.x[i][0][k] - 4.0 * PI * t * J.x[i][0][k] + c * t * ((B.z[i][0][k] - B.z[i][ny - 1][k]) / (2.0 * steps.dy) - (B.y[i][0][k] - B.y[i][0][k - 1]) / (2.0 * steps.dz));
						E.y[i][0][k] = E.y[i][0][k] - 4.0 * PI * t * J.y[i][0][k] + c * t * ((B.x[i][0][k] - B.x[i][0][k - 1]) / (2.0 * steps.dz) - (B.z[i][0][k] - B.z[i - 1][0][k]) / (2.0 * steps.dx));
					}
					if (k == 0) {
						E.x[i][0][k] = E.x[i][0][k] - 4.0 * PI * t * J.x[i][0][k] + c * t * ((B.z[i][0][k] - B.z[i][ny - 1][k]) / (2.0 * steps.dy) - (B.y[i][0][k] - B.y[i][0][nz - 1]) / (2.0 * steps.dz));
						E.y[i][0][k] = E.y[i][0][k] - 4.0 * PI * t * J.y[i][0][k] + c * t * ((B.x[i][0][k] - B.x[i][0][nz - 1]) / (2.0 * steps.dz) - (B.z[i][0][k] - B.z[i - 1][0][k]) / (2.0 * steps.dx));
					}
				}
			}
			if (i == 0) {
				for (int k = 0; k < nz; ++k) {
					E.z[i][0][k] = E.z[i][0][k] - 4.0 * PI * t * J.z[i][0][k] + c * t * ((B.y[i][0][k] - B.y[nx - 1][0][k]) / (2.0 * steps.dx) - (B.x[i][0][k] - B.x[i][ny - 1][k]) / (2.0 * steps.dy));
					if (k != 0) {
						E.x[i][0][k] = E.x[i][0][k] - 4.0 * PI * t * J.x[i][0][k] + c * t * ((B.z[i][0][k] - B.z[i][ny - 1][k]) / (2.0 * steps.dy) - (B.y[i][0][k] - B.y[i][0][k - 1]) / (2.0 * steps.dz));
						E.y[i][0][k] = E.y[i][0][k] - 4.0 * PI * t * J.y[i][0][k] + c * t * ((B.x[i][0][k] - B.x[i][0][k - 1]) / (2.0 * steps.dz) - (B.z[i][0][k] - B.z[nx - 1][0][k]) / (2.0 * steps.dx));
					}
					if (k == 0) {
						E.x[i][0][k] = E.x[i][0][k] - 4.0 * PI * t * J.x[i][0][k] + c * t * ((B.z[i][0][k] - B.z[i][ny - 1][k]) / (2.0 * steps.dy) - (B.y[i][0][k] - B.y[i][0][nz - 1]) / (2.0 * steps.dz));
						E.y[i][0][k] = E.y[i][0][k] - 4.0 * PI * t * J.y[i][0][k] + c * t * ((B.x[i][0][k] - B.x[i][0][nz - 1]) / (2.0 * steps.dz) - (B.z[i][0][k] - B.z[nx - 1][0][k]) / (2.0 * steps.dx));
					}
				}
			}
		}

		// k == 0
		for (int j = 0; j < ny; ++j) {
			if (j != 0) {
				for (int i = 0; i < nx; ++i) {
					E.x[i][j][0] = E.x[i][j][0] - 4.0 * PI * t * J.x[i][j][0] + c * t * ((B.z[i][j][0] - B.z[i][j - 1][0]) / (2.0 * steps.dy) - (B.y[i][j][0] - B.y[i][j][nz - 1]) / (2.0 * steps.dz));
					if (i != 0) {
						E.y[i][j][0] = E.y[i][j][0] - 4.0 * PI * t * J.y[i][j][0] + c * t * ((B.x[i][j][0] - B.x[i][j][nz - 1]) / (2.0 * steps.dz) - (B.z[i][j][0] - B.z[i - 1][j][0]) / (2.0 * steps.dx));
						E.z[i][j][0] = E.z[i][j][0] - 4.0 * PI * t * J.z[i][j][0] + c * t * ((B.y[i][j][0] - B.y[i - 1][j][0]) / (2.0 * steps.dx) - (B.x[i][j][0] - B.x[i][j - 1][0]) / (2.0 * steps.dy));
					}
					if (i == 0) {
						E.y[i][j][0] = E.y[i][j][0] - 4.0 * PI * t * J.y[i][j][0] + c * t * ((B.x[i][j][0] - B.x[i][j][nz - 1]) / (2.0 * steps.dz) - (B.z[i][j][0] - B.z[nx - 1][j][0]) / (2.0 * steps.dx));
						E.z[i][j][0] = E.z[i][j][0] - 4.0 * PI * t * J.z[i][j][0] + c * t * ((B.y[i][j][0] - B.y[nx - 1][j][0]) / (2.0 * steps.dx) - (B.x[i][j][0] - B.x[i][j - 1][0]) / (2.0 * steps.dy));
					}
				}
			}

			if (j == 0) {
				for (int i = 0; i < nx; ++i) {
					E.x[i][j][0] = E.x[i][j][0] - 4.0 * PI * t * J.x[i][j][0] + c * t * ((B.z[i][j][0] - B.z[i][ny - 1][0]) / (2.0 * steps.dy) - (B.y[i][j][0] - B.y[i][j][nz - 1]) / (2.0 * steps.dz));
					if (i != 0) {
						E.y[i][j][0] = E.y[i][j][0] - 4.0 * PI * t * J.y[i][j][0] + c * t * ((B.x[i][j][0] - B.x[i][j][nz - 1]) / (2.0 * steps.dz) - (B.z[i][j][0] - B.z[i - 1][j][0]) / (2.0 * steps.dx));
						E.z[i][j][0] = E.z[i][j][0] - 4.0 * PI * t * J.z[i][j][0] + c * t * ((B.y[i][j][0] - B.y[i - 1][j][0]) / (2.0 * steps.dx) - (B.x[i][j][0] - B.x[i][ny - 1][0]) / (2.0 * steps.dy));
					}
					if (i == 0) {
						E.y[i][j][0] = E.y[i][j][0] - 4.0 * PI * t * J.y[i][j][0] + c * t * ((B.x[i][j][0] - B.x[i][j][nz - 1]) / (2.0 * steps.dz) - (B.z[i][j][0] - B.z[nx - 1][j][0]) / (2.0 * steps.dx));
						E.z[i][j][0] = E.z[i][j][0] - 4.0 * PI * t * J.z[i][j][0] + c * t * ((B.y[i][j][0] - B.y[nx - 1][j][0]) / (2.0 * steps.dx) - (B.x[i][j][0] - B.x[i][ny - 1][0]) / (2.0 * steps.dy));
					}
				}
			}
		}

		return (*this);
	}

	elec_magn_field& get_boundary_conditions_B(double t) {

		// i == nx - 1
		for (int j = 0; j < ny; ++j) {
			if (j != ny - 1) {
				for (int k = 0; k < nz; ++k) {
					B.z[nx - 1][j][k] = B.z[nx - 1][j][k] + c * t * ((E.x[nx - 1][j + 1][k] - E.x[nx - 1][j][k]) / (2.0 * steps.dy) - (E.y[0][j][k] - E.y[nx - 1][j][k]) / (2.0 * steps.dx));
					if (k != nz - 1) {
						B.x[nx - 1][j][k] = B.x[nx - 1][j][k] + c * t * ((E.y[nx - 1][j][k + 1] - E.y[nx - 1][j][k]) / (2.0 * steps.dz) - (E.z[nx - 1][j + 1][k] - E.z[nx - 1][j][k]) / (2.0 * steps.dy));
						B.y[nx - 1][j][k] = B.y[nx - 1][j][k] + c * t * ((E.z[0][j][k] - E.z[nx - 1][j][k]) / (2.0 * steps.dx) - (E.x[nx - 1][j][k + 1] - E.x[nx - 1][j][k]) / (2.0 * steps.dz));
					}
					if (k == nz - 1) {
						B.x[nx - 1][j][k] = B.x[nx - 1][j][k] + c * t * ((E.y[nx - 1][j][0] - E.y[nx - 1][j][k]) / (2.0 * steps.dz) - (E.z[nx - 1][j + 1][k] - E.z[nx - 1][j][k]) / (2.0 * steps.dy));
						B.y[nx - 1][j][k] = B.y[nx - 1][j][k] + c * t * ((E.z[0][j][k] - E.z[nx - 1][j][k]) / (2.0 * steps.dx) - (E.x[nx - 1][j][0] - E.x[nx - 1][j][k]) / (2.0 * steps.dz));
					}
				}
			}
			if (j == ny - 1) {
				for (int k = 0; k < nz; ++k) {
					B.z[nx - 1][j][k] = B.z[nx - 1][j][k] + c * t * ((E.x[nx - 1][0][k] - E.x[nx - 1][j][k]) / (2.0 * steps.dy) - (E.y[0][j][k] - E.y[nx - 1][j][k]) / (2.0 * steps.dx));
					if (k != nz - 1) {
						B.x[nx - 1][j][k] = B.x[nx - 1][j][k] + c * t * ((E.y[nx - 1][j][k + 1] - E.y[nx - 1][j][k]) / (2.0 * steps.dz) - (E.z[nx - 1][0][k] - E.z[nx - 1][j][k]) / (2.0 * steps.dy));
						B.y[nx - 1][j][k] = B.y[nx - 1][j][k] + c * t * ((E.z[0][j][k] - E.z[nx - 1][j][k]) / (2.0 * steps.dx) - (E.x[nx - 1][j][k + 1] - E.x[nx - 1][j][k]) / (2.0 * steps.dz));
					}
					if (k == nz - 1) {
						B.x[nx - 1][j][k] = B.x[nx - 1][j][k] + c * t * ((E.y[nx - 1][j][0] - E.y[nx - 1][j][k]) / (2.0 * steps.dz) - (E.z[nx - 1][0][k] - E.z[nx - 1][j][k]) / (2.0 * steps.dy));
						B.y[nx - 1][j][k] = B.y[nx - 1][j][k] + c * t * ((E.z[0][j][k] - E.z[nx - 1][j][k]) / (2.0 * steps.dx) - (E.x[nx - 1][j][0] - E.x[nx - 1][j][k]) / (2.0 * steps.dz));
					}
				}
			}
		}

		// j == ny - 1
		for (int i = 0; i < nx; ++i) {
			if (i != nx - 1) {
				for (int k = 0; k < nz; ++k) {
					B.z[i][ny - 1][k] = B.z[i][ny - 1][k] + c * t * ((E.x[i][0][k] - E.x[i][ny - 1][k]) / (2.0 * steps.dy) - (E.y[i + 1][ny - 1][k] - E.y[i][ny - 1][k]) / (2.0 * steps.dx));
					if (k != nz - 1) {
						B.x[i][ny - 1][k] = B.x[i][ny - 1][k] + c * t * ((E.y[i][ny - 1][k + 1] - E.y[i][ny - 1][k]) / (2.0 * steps.dz) - (E.z[i][0][k] - E.z[i][ny - 1][k]) / (2.0 * steps.dy));
						B.y[i][ny - 1][k] = B.y[i][ny - 1][k] + c * t * ((E.z[i + 1][ny - 1][k] - E.z[i][ny - 1][k]) / (2.0 * steps.dx) - (E.x[i][ny - 1][k + 1] - E.x[i][ny - 1][k]) / (2.0 * steps.dz));
					}
					if (k == nz - 1) {
						B.x[i][ny - 1][k] = B.x[i][ny - 1][k] + c * t * ((E.y[i][ny - 1][0] - E.y[i][ny - 1][k]) / (2.0 * steps.dz) - (E.z[i][0][k] - E.z[i][ny - 1][k]) / (2.0 * steps.dy));
						B.y[i][ny - 1][k] = B.y[i][ny - 1][k] + c * t * ((E.z[i + 1][ny - 1][k] - E.z[i][ny - 1][k]) / (2.0 * steps.dx) - (E.x[i][ny - 1][0] - E.x[i][ny - 1][k]) / (2.0 * steps.dz));
					}
				}
			}
			if (i == nx - 1) {
				for (int k = 0; k < nz; ++k) {
					B.z[i][ny - 1][k] = B.z[i][ny - 1][k] + c * t * ((E.x[i][0][k] - E.x[i][ny - 1][k]) / (2.0 * steps.dy) - (E.y[0][ny - 1][k] - E.y[i][ny - 1][k]) / (2.0 * steps.dx));
					if (k != nz - 1) {
						B.x[i][ny - 1][k] = B.x[i][ny - 1][k] + c * t * ((E.y[i][ny - 1][k + 1] - E.y[i][ny - 1][k]) / (2.0 * steps.dz) - (E.z[i][0][k] - E.z[i][ny - 1][k]) / (2.0 * steps.dy));
						B.y[i][ny - 1][k] = B.y[i][ny - 1][k] + c * t * ((E.z[0][ny - 1][k] - E.z[i][ny - 1][k]) / (2.0 * steps.dx) - (E.x[i][ny - 1][k + 1] - E.x[i][ny - 1][k]) / (2.0 * steps.dz));
					}
					if (k == nz - 1) {
						B.x[i][ny - 1][k] = B.x[i][ny - 1][k] + c * t * ((E.y[i][ny - 1][0] - E.y[i][ny - 1][k]) / (2.0 * steps.dz) - (E.z[i][0][k] - E.z[i][ny - 1][k]) / (2.0 * steps.dy));
						B.y[i][ny - 1][k] = B.y[i][ny - 1][k] + c * t * ((E.z[0][ny - 1][k] - E.z[i][ny - 1][k]) / (2.0 * steps.dx) - (E.x[i][ny - 1][0] - E.x[i][ny - 1][k]) / (2.0 * steps.dz));
					}
				}
			}
		}


		// k == nz - 1
		for (int i = 0; i < nx; ++i) {
			if (i != nx - 1) {
				for (int j = 0; j < ny; ++j) {
					B.y[i][j][nz - 1] = B.y[i][j][nz - 1] + c * t * ((E.z[i + 1][j][nz - 1] - E.z[i][j][nz - 1]) / (2.0 * steps.dx) - (E.x[i][j][0] - E.x[i][j][nz - 1]) / (2.0 * steps.dz));
					if (j != ny - 1) {
						B.x[i][j][nz - 1] = B.x[i][j][nz - 1] + c * t * ((E.y[i][j][0] - E.y[i][j][nz - 1]) / (2.0 * steps.dz) - (E.z[i][j + 1][nz - 1] - E.z[i][j][nz - 1]) / (2.0 * steps.dy));
						B.z[i][j][nz - 1] = B.z[i][j][nz - 1] + c * t * ((E.x[i][j + 1][nz - 1] - E.x[i][j][nz - 1]) / (2.0 * steps.dy) - (E.y[i + 1][j][nz - 1] - E.y[i][j][nz - 1]) / (2.0 * steps.dx));
					}
					if (j == ny - 1) {
						B.x[i][j][nz - 1] = B.x[i][j][nz - 1] + c * t * ((E.y[i][j][0] - E.y[i][j][nz - 1]) / (2.0 * steps.dz) - (E.z[i][0][nz - 1] - E.z[i][j][nz - 1]) / (2.0 * steps.dy));
						B.z[i][j][nz - 1] = B.z[i][j][nz - 1] + c * t * ((E.x[i][0][nz - 1] - E.x[i][j][nz - 1]) / (2.0 * steps.dy) - (E.y[i + 1][j][nz - 1] - E.y[i][j][nz - 1]) / (2.0 * steps.dx));
					}
				}
			}
			if (i == nx - 1) {
				for (int j = 0; j < ny; ++j) {
					B.y[i][j][nz - 1] = B.y[i][j][nz - 1] + c * t * ((E.z[0][j][nz - 1] - E.z[i][j][nz - 1]) / (2.0 * steps.dx) - (E.x[i][j][0] - E.x[i][j][nz - 1]) / (2.0 * steps.dz));
					if (j != ny - 1) {
						B.x[i][j][nz - 1] = B.x[i][j][nz - 1] + c * t * ((E.y[i][j][0] - E.y[i][j][nz - 1]) / (2.0 * steps.dz) - (E.z[i][j + 1][nz - 1] - E.z[i][j][nz - 1]) / (2.0 * steps.dy));
						B.z[i][j][nz - 1] = B.z[i][j][nz - 1] + c * t * ((E.x[i][j + 1][nz - 1] - E.x[i][j][nz - 1]) / (2.0 * steps.dy) - (E.y[0][j][nz - 1] - E.y[i][j][nz - 1]) / (2.0 * steps.dx));
					}
					if (j == ny - 1) {
						B.x[i][j][nz - 1] = B.x[i][j][nz - 1] + c * t * ((E.y[i][j][0] - E.y[i][j][nz - 1]) / (2.0 * steps.dz) - (E.z[i][0][nz - 1] - E.z[i][j][nz - 1]) / (2.0 * steps.dy));
						B.z[i][j][nz - 1] = B.z[i][j][nz - 1] + c * t * ((E.x[i][0][nz - 1] - E.x[i][j][nz - 1]) / (2.0 * steps.dy) - (E.y[0][j][nz - 1] - E.y[i][j][nz - 1]) / (2.0 * steps.dx));
					}
				}
			}
		}
		return (*this);
	}

	elec_magn_field& start_FDTD() {
		for (int i = 0; i < nx - 1; i++)
			for (int j = 0; j < ny - 1; j++)
				for (int k = 0; k < nz - 1; k++) {
					B.x[i][j][k] = B.x[i][j][k] + c * steps.dt * 0.5 * ((E.y[i][j][k + 1] - E.y[i][j][k]) / (2.0 * steps.dz) - (E.z[i][j + 1][k] - E.z[i][j][k]) / (2.0 * steps.dy));
					B.y[i][j][k] = B.y[i][j][k] + c * steps.dt * 0.5 *((E.z[i + 1][j][k] - E.z[i][j][k]) / (2.0 * steps.dx) - (E.x[i][j][k + 1] - E.x[i][j][k]) / (2.0 * steps.dz));
					B.z[i][j][k] = B.z[i][j][k] + c * steps.dt * 0.5 * ((E.x[i][j + 1][k] - E.x[i][j][k]) / (2.0 * steps.dy) - (E.y[i + 1][j][k] - E.y[i][j][k]) / (2.0 * steps.dx));
				}
		get_boundary_conditions_B(steps.dt * 0.5);

		for (int i = 1; i < nx; i++)
			for (int j = 1; j < ny; j++)
				for (int k = 1; k < nz; k++) {
					E.x[i][j][k] = E.x[i][j][k] - 4.0 * PI * steps.dt * J.x[i][j][k] + c * steps.dt * ((B.z[i][j][k] - B.z[i][j - 1][k]) / (2.0 * steps.dy) - (B.y[i][j][k] - B.y[i][j][k - 1]) / (2.0 * steps.dz));
					E.y[i][j][k] = E.y[i][j][k] - 4.0 * PI * steps.dt * J.y[i][j][k] + c * steps.dt * ((B.x[i][j][k] - B.x[i][j][k - 1]) / (2.0 * steps.dz) - (B.z[i][j][k] - B.z[i - 1][j][k]) / (2.0 * steps.dx));
					E.z[i][j][k] = E.z[i][j][k] - 4.0 * PI * steps.dt * J.z[i][j][k] + c * steps.dt * ((B.y[i][j][k] - B.y[i - 1][j][k]) / (2.0 * steps.dx) - (B.x[i][j][k] - B.x[i][j - 1][k]) / (2.0 * steps.dy));
				}

		get_boundary_conditions_E(steps.dt);

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

		get_boundary_conditions_B(steps.dt);

		// обновление E

		for (int i = 1; i < nx; i++)
			for (int j = 1; j < ny; j++)
				for (int k = 1; k < nz; k++) {
					E.x[i][j][k] = E.x[i][j][k] - 4.0 * PI * steps.dt * J.x[i][j][k] + c * steps.dt * ((B.z[i][j][k] - B.z[i][j - 1][k]) / (2.0 * steps.dy) - (B.y[i][j][k] - B.y[i][j][k - 1]) / (2.0 * steps.dz));
					E.y[i][j][k] = E.y[i][j][k] - 4.0 * PI * steps.dt * J.y[i][j][k] + c * steps.dt * ((B.x[i][j][k] - B.x[i][j][k - 1]) / (2.0 * steps.dz) - (B.z[i][j][k] - B.z[i - 1][j][k]) / (2.0 * steps.dx));
					E.z[i][j][k] = E.z[i][j][k] - 4.0 * PI * steps.dt * J.z[i][j][k] + c * steps.dt * ((B.y[i][j][k] - B.y[i - 1][j][k]) / (2.0 * steps.dx) - (B.x[i][j][k] - B.x[i][j - 1][k]) / (2.0 * steps.dy));
				}

		get_boundary_conditions_E(steps.dt);

		return (*this);
	}


	elec_magn_field& modify_Jz(int m) {
		double t = steps.dt * (m - 0.5);
		double x, y;
		int z = 0;
		for (int i = 0; i < nx; ++i)
			for (int j = 0; j < ny; ++j) {
				//for (int k = 0; k < nz; ++k) {
				x = area.a + i * steps.dx;
				y = area.a + j * steps.dy;
				// z = area.a + k * steps.dz;
				if ((x >= -Tx / 4.0) && (x <= Tx / 4.0) && (y >= -Ty / 4.0) && (y <= Ty / 4.0) && (z >= -Tz / 4.0) && (z <= Tz / 4.0)) {
					J.z(i, j, 0) = sin(2 * PI * t / T) * pow(cos(2 * PI * x / Tx), 2) * pow(cos(2 * PI * y / Ty), 2) * pow(cos(2 * PI * z / Tz), 2);
				}
			}
		return (*this);
	}

	elec_magn_field& set_zero_Jz() {
		for (int i = 0; i < nx; ++i)
			for (int j = 0; j < ny; ++j)
				for (int k = 0; k < nz; ++k)
					J.z(i, j, k) = 0.0;
		return (*this);
	}

	point w(int i, int j, int k) {
		point res;
		if (i <= nx / 2)
			res.x = (2 * PI * i / (area.b - area.a));
		else
			res.x = (2 * PI * (i - nx) / (area.b - area.a));
		if (j <= ny / 2)
			res.y = (2 * PI * j / (area.b - area.a));
		else
			res.y = (2 * PI * (j - ny) / (area.b - area.a));
		if (k <= nz / 2)
			res.z = (2 * PI * k / (area.b - area.a));
		else
			res.z = (2 * PI * (k - nz) / (area.b - area.a));

		return res;
	}

	double K(int i, int j, int k) {
		return std::sqrt(w(i, j, k).x * w(i, j, k).x + w(i, j, k).y * w(i, j, k).y + w(i, j, k).z * w(i, j, k).z);
	}

	point K_norm(int i, int j, int k) {
		point res;
		res.x = w(i, j, k).x / K(i, j, k);
		res.y = w(i, j, k).y / K(i, j, k);
		res.z = w(i, j, k).z / K(i, j, k);
		return res;
	}

	point vect_mult_K(int i, int j, int k, double x, double y, double z) {
		point res;
		res.x = K_norm(i, j, k).y * z - K_norm(i, j, k).z * y;
		res.y = K_norm(i, j, k).x * z - K_norm(i, j, k).z * x;
		res.z = K_norm(i, j, k).x * y - K_norm(i, j, k).y * x;

		return res;
	}

	double scalar_mult_K(int i, int j, int k, double x, double y, double z) {
		return (K_norm(i, j, k).x * x + K_norm(i, j, k).y * y + K_norm(i, j, k).z * z);
	}

	elec_magn_field& PSTD() {

		/*double* Ex_for_dft = (double*)malloc(sizeof(double) * nx * ny * nz);
		double* Ey_for_dft = (double*)malloc(sizeof(double) * nx * ny * nz);
		double* Ez_for_dft = (double*)malloc(sizeof(double) * nx * ny * nz);
		double* Bx_for_dft = (double*)malloc(sizeof(double) * nx * ny * nz);
		double* By_for_dft = (double*)malloc(sizeof(double) * nx * ny * nz);
		double* Bz_for_dft = (double*)malloc(sizeof(double) * nx * ny * nz);
		double* Jx_for_dft = (double*)malloc(sizeof(double) * nx * ny * nz);
		double* Jy_for_dft = (double*)malloc(sizeof(double) * nx * ny * nz);
		double* Jz_for_dft = (double*)malloc(sizeof(double) * nx * ny * nz);

		int q = 0;
		for (int i = 0; i < nx; i++) {
			for (int j = 0; j < ny; j++) {
				for (int k = 0; k < nz; k++) {
					Ex_for_dft[q] = E.x(i, j, k);
					Ez_for_dft[q] = E.y(i, j, k);
					Ez_for_dft[q] = E.z(i, j, k);
					Bx_for_dft[q] = B.x(i, j, k);
					By_for_dft[q] = B.y(i, j, k);
					Bz_for_dft[q] = B.z(i, j, k);
					Jx_for_dft[q] = J.x(i, j, k);
					Jy_for_dft[q] = J.y(i, j, k);
					Jz_for_dft[q] = J.z(i, j, k);
					q++;
				}
			}
		}*/
		//PSTD


			int q = 0;

			fftw_complex *Bx_out = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * size);
			fftw_complex *By_out = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * size);
			fftw_complex *Bz_out = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * size);
			fftw_complex *Ex_out = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * size);
			fftw_complex *Ey_out = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * size);
			fftw_complex *Ez_out = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * size);
			fftw_complex *Jx_out = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * size);
			fftw_complex *Jy_out = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * size);
			fftw_complex *Jz_out = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * size);

			fftw_plan Ex, Ey, Ez, Bx, By, Bz, Jx, Jy, Jz;
			fftw_plan Ex_, Ey_, Ez_, Bx_, By_, Bz_, Jx_, Jy_, Jz_;

			Ex = fftw_plan_dft_r2c_3d(nx, ny, nz, &E.x.data[0], Ex_out, FFTW_ESTIMATE);
			Ey = fftw_plan_dft_r2c_3d(nx, ny, nz, &E.y.data[0], Ey_out, FFTW_ESTIMATE);
			Ez = fftw_plan_dft_r2c_3d(nx, ny, nz, &E.z.data[0], Ez_out, FFTW_ESTIMATE);
			Bx = fftw_plan_dft_r2c_3d(nx, ny, nz, &B.x.data[0], Bx_out, FFTW_ESTIMATE);
			By = fftw_plan_dft_r2c_3d(nx, ny, nz, &B.y.data[0], By_out, FFTW_ESTIMATE);
			Bz = fftw_plan_dft_r2c_3d(nx, ny, nz, &B.z.data[0], Bz_out, FFTW_ESTIMATE);
			Jx = fftw_plan_dft_r2c_3d(nx, ny, nz, &J.x.data[0], Jx_out, FFTW_ESTIMATE);
			Jy = fftw_plan_dft_r2c_3d(nx, ny, nz, &J.y.data[0], Jy_out, FFTW_ESTIMATE);
			Jz = fftw_plan_dft_r2c_3d(nx, ny, nz, &J.z.data[0], Jz_out, FFTW_ESTIMATE);

			fftw_execute(Ex);
			fftw_execute(Ey);
			fftw_execute(Ez);
			fftw_execute(Bx);
			fftw_execute(By);
			fftw_execute(Bz);
			fftw_execute(Jx);
			fftw_execute(Jy);
			fftw_execute(Jz);

			/*int i = 1;
			int k = int(T / steps.dt);

			for (double t = f.steps.dt; t <= t1; t += f.steps.dt) {
				if (i <= k) {
					f.modify_Jz(i);

				}
				else {
					f.set_zero_Jz();
				}*/
				for (int _i = 0; _i < nx; _i++) {
					for (int _j = 0; _j < ny; _j++) {
						for (int _k = 0; _k < nz; _k++) {

							Bx_out[q][0] += c * 0.5 * steps.dt * vect_mult_K(_i, _j, _k, Ex_out[q][1], Ey_out[q][1], Ez_out[q][1]).x;

							Bx_out[q][1] += -c * 0.5 * steps.dt * vect_mult_K(_i, _j, _k, Ex_out[q][0], Ey_out[q][0], Ez_out[q][0]).x;

							By_out[q][0] += c * 0.5 * steps.dt * vect_mult_K(_i, _j, _k, Ex_out[q][1], Ey_out[q][1], Ez_out[q][1]).y;

							By_out[q][1] += -c * 0.5 * steps.dt * vect_mult_K(_i, _j, _k, Ex_out[q][0], Ey_out[q][0], Ez_out[q][0]).y;

							Bz_out[q][0] += c * 0.5 * steps.dt * vect_mult_K(_i, _j, _k, Ex_out[q][1], Ey_out[q][1], Ez_out[q][1]).z;

							Bz_out[q][1] += -c * 0.5 * steps.dt * vect_mult_K(_i, _j, _k, Ex_out[q][0], Ey_out[q][0], Ez_out[q][0]).z;

							Ex_out[q][0] += -c * steps.dt * vect_mult_K(_i, _j, _k, Bx_out[q][1], By_out[q][1], \
								Bz_out[q][1]).x - 4 * PI * steps.dt * Jx_out[q][0];

							Ex_out[q][1] += c * steps.dt * vect_mult_K(_i, _j, _k, Bx_out[q][0], By_out[q][0], \
								Bz_out[q][0]).x - 4 * PI * steps.dt * Jx_out[q][1];

							Ey_out[q][0] += -c * steps.dt * vect_mult_K(_i, _j, _k, Bx_out[q][1], By_out[q][1], \
								Bz_out[q][1]).y - 4 * PI * steps.dt * Jy_out[q][0];

							Ey_out[q][1] += c * steps.dt * vect_mult_K(_i, _j, _k, Bx_out[q][0], By_out[q][0], \
								Bz_out[q][0]).y - 4 * PI * steps.dt * Jy_out[q][1];

							Ez_out[q][0] += -c * steps.dt * vect_mult_K(_i, _j, _k, Bx_out[q][1], By_out[q][1], \
								Bz_out[q][1]).z - 4 * PI * steps.dt * Jz_out[q][0];

							Ez_out[q][1] += c * steps.dt * vect_mult_K(_i, _j, _k, Bx_out[q][0], By_out[q][0], \
								Bz_out[q][0]).x - 4 * PI * steps.dt * Jz_out[q][1];

							q++;

						}
					}
				}
			//}

			Ex_ = fftw_plan_dft_c2r_3d(nx, ny, nz, Ex_out, &E.x.data[0], FFTW_ESTIMATE);
			Ey_ = fftw_plan_dft_c2r_3d(nx, ny, nz, Ey_out, &E.y.data[0], FFTW_ESTIMATE);
			Ez_ = fftw_plan_dft_c2r_3d(nx, ny, nz, Ez_out, &E.z.data[0], FFTW_ESTIMATE);
			Bx_ = fftw_plan_dft_c2r_3d(nx, ny, nz, Bx_out, &B.x.data[0], FFTW_ESTIMATE);
			By_ = fftw_plan_dft_c2r_3d(nx, ny, nz, By_out, &B.y.data[0], FFTW_ESTIMATE);
			Bz_ = fftw_plan_dft_c2r_3d(nx, ny, nz, Bz_out, &B.z.data[0], FFTW_ESTIMATE);
			Jx_ = fftw_plan_dft_c2r_3d(nx, ny, nz, Jx_out, &J.x.data[0], FFTW_ESTIMATE);
			Jy_ = fftw_plan_dft_c2r_3d(nx, ny, nz, Jy_out, &J.y.data[0], FFTW_ESTIMATE);
			Jz_ = fftw_plan_dft_c2r_3d(nx, ny, nz, Jz_out, &J.z.data[0], FFTW_ESTIMATE);

			fftw_execute(Ex_);
			fftw_execute(Ey_);
			fftw_execute(Ez_);
			fftw_execute(Bx_);
			fftw_execute(By_);
			fftw_execute(Bz_);
			fftw_execute(Jx_);
			fftw_execute(Jy_);
			fftw_execute(Jz_);

			for (int _i = 0; _i < size; _i++) {
				E.x[_i] /= size;
				E.y[_i] /= size;
				E.z[_i] /= size;
				B.x[_i] /= size;
				B.y[_i] /= size;
				B.z[_i] /= size;
				J.x[_i] /= size;
				J.y[_i] /= size;
				J.z[_i] /= size;
			}

			/*q = 0;
			for (int _i = 0; _i < nx; _i++)
				for (int _j = 0; _j < ny; _j++)
					for (int _k = 0; _k < nz; _k++) {
						E.x[_i][_j][_k] = Ex_for_dft[q];
						E.y[_i][_j][_k] = Ey_for_dft[q];
						E.z[_i][_j][_k] = Ez_for_dft[q];
						B.x[_i][_j][_k] = Bx_for_dft[q];
						B.y[_i][_j][_k] = By_for_dft[q];
						B.z[_i][_j][_k] = Bz_for_dft[q];
						J.x[_i][_j][_k] = Jx_for_dft[q];
						J.y[_i][_j][_k] = Jy_for_dft[q];
						J.z[_i][_j][_k] = Jz_for_dft[q];
						q++;
					}*/

			fftw_destroy_plan(Ex);
			fftw_destroy_plan(Ey);
			fftw_destroy_plan(Ez);
			fftw_destroy_plan(Bx);
			fftw_destroy_plan(By);
			fftw_destroy_plan(Bz);
			fftw_destroy_plan(Jx);
			fftw_destroy_plan(Jy);
			fftw_destroy_plan(Jz);

			fftw_destroy_plan(Ex_);
			fftw_destroy_plan(Ey_);
			fftw_destroy_plan(Ez_);
			fftw_destroy_plan(Bx_);
			fftw_destroy_plan(By_);
			fftw_destroy_plan(Bz_);
			fftw_destroy_plan(Jx_);
			fftw_destroy_plan(Jy_);
			fftw_destroy_plan(Jz_);

			fftw_free(Ex_out);
			fftw_free(Ey_out);
			fftw_free(Ez_out);
			fftw_free(Bx_out);
			fftw_free(By_out);
			fftw_free(Bz_out);
			fftw_free(Jx_out);
			fftw_free(Jy_out);
			fftw_free(Jz_out);	

		/*free(Ex_for_dft);
		free(Ey_for_dft);
		free(Ez_for_dft);
		free(Bx_for_dft);
		free(By_for_dft);
		free(Bz_for_dft);
		free(Jx_for_dft);
		free(Jy_for_dft);
		free(Jz_for_dft);
		*/
		//for (int k = 0; k < nz; k++) {
				//	double C = std::cos(K(i, j, k) * c * steps.dt * 0.5 / 2.0);
				//	double S = std::sin(K(i, j, k) * c * steps.dt * 0.5 / 2.0);
				//	/*B_x_out[0][0] = C * B_x_out[0][0] + S * vect_mult_K(i, j, k, E_x_out[0][1], E_y_out[0][1], \
				//		E_z_out[0][1]).x - 4 * PI * (1 - C) / (K(i, j, k) * c) * vect_mult_K(i, j, k, J_x_out[0][1], J_y_out[0][1], \
				//			J_z_out[0][1]).x;
				//	B_x_out[0][1] = C * B_x_out[0][1] - S * vect_mult_K(i, j, k, E_x_out[0][0], E_y_out[0][0], \
				//		E_z_out[0][0]).x + 4 * PI * (1 - C) / (K(i, j, k) * c) * vect_mult_K(i, j, k, J_x_out[0][0], J_y_out[0][0], \
				//			J_z_out[0][0]).x;
				//	B_y_out[0][0] = C * B_y_out[0][0] + S * vect_mult_K(i, j, k, E_x_out[0][1], E_y_out[0][1], \
				//		E_z_out[0][1]).y - 4 * PI * (1 - C) / (K(i, j, k) * c) * vect_mult_K(i, j, k, J_x_out[0][1], J_y_out[0][1], \
				//			J_z_out[0][1]).y;
				//	B_y_out[0][1] = C * B_y_out[0][1] - S * vect_mult_K(i, j, k, E_x_out[0][0], E_y_out[0][0], \
				//		E_z_out[0][0]).y + 4 * PI * (1 - C) / (K(i, j, k) * c) * vect_mult_K(i, j, k, J_x_out[0][0], J_y_out[0][0], \
				//			J_z_out[0][0]).y;
				//	B_z_out[0][0] = C * B_z_out[0][0] + S * vect_mult_K(i, j, k, E_x_out[0][1], E_y_out[0][1], \
				//		E_z_out[0][1]).z - 4 * PI * (1 - C) / (K(i, j, k) * c) * vect_mult_K(i, j, k, J_x_out[0][1], J_y_out[0][1], \
				//			J_z_out[0][1]).z;
				//	B_z_out[0][1] = C * B_z_out[0][1] - S * vect_mult_K(i, j, k, E_x_out[0][0], E_y_out[0][0], \
				//		E_z_out[0][0]).z + 4 * PI * (1 - C) / (K(i, j, k) * c) * vect_mult_K(i, j, k, J_x_out[0][0], J_y_out[0][0], \
				//			J_z_out[0][0]).z;*/
				//	Bx_out[k][0] += c * 0.5 * steps.dt * vect_mult_K(i, j, k, Ex_out[k][1], Ey_out[k][1], Ez_out[k][1]).x;
				//	Bx_out[k][1] += - c * 0.5 * steps.dt * vect_mult_K(i, j, k, Ex_out[k][0], Ey_out[k][0], Ez_out[k][0]).x;
				//	By_out[k][0] += c * 0.5 * steps.dt * vect_mult_K(i, j, k, Ex_out[k][1], Ey_out[k][1], Ez_out[k][1]).y;
				//	By_out[k][1] += - c * 0.5 * steps.dt * vect_mult_K(i, j, k, Ex_out[k][0], Ey_out[k][0], Ez_out[k][0]).y;
				//	Bz_out[k][0] += c * 0.5 * steps.dt * vect_mult_K(i, j, k, Ex_out[k][1], Ey_out[k][1], Ez_out[k][1]).z;
				//	Bz_out[k][1] += - c * 0.5 * steps.dt * vect_mult_K(i, j, k, Ex_out[k][0], Ey_out[k][0], Ez_out[k][0]).z;
			//	for (int k = 0; k < nz; k++) {
			//		double C = std::cos(K(i, j, k) * c * steps.dt / 2.0);
			//		double S = std::sin(K(i, j, k) * c * steps.dt / 2.0);
			//		/*E_x_out[0][0] = C * E_x_out[0][0] - S * vect_mult_K(i, j, k, B_x_out[0][1], B_y_out[0][1], \
			//			B_z_out[0][1]).x - 4 * PI * S / (K(i, j, k) * c) * J_x_out[0][0] + (1 - C) * K_norm(i, j, k).x * scalar_mult_K(i, j, k, \
			//				E_x_out[0][0], E_y_out[0][0], E_z_out[0][0]) + 4 * PI * K_norm(i, j, k).x * scalar_mult_K(i, j, k, \
			//					J_x_out[0][0], J_y_out[0][0], J_z_out[0][0]) * (S / (K(i, j, k) * c) - steps.dt);
			//		E_x_out[0][1] = C * E_x_out[0][1] + S * vect_mult_K(i, j, k, B_x_out[0][0], B_y_out[0][0], \
			//			B_z_out[0][0]).x - 4 * PI * S / (K(i, j, k) * c) * J_x_out[0][1] + (1 - C) * K_norm(i, j, k).x * scalar_mult_K(i, j, k, \
			//				E_x_out[0][1], E_y_out[0][1], E_z_out[0][1]) + 4 * PI * K_norm(i, j, k).x * scalar_mult_K(i, j, k, \
			//					J_x_out[0][1], J_y_out[0][1], J_z_out[0][1]) * (S / (K(i, j, k) * c) - steps.dt);
			//		E_y_out[0][0] = C * E_y_out[0][0] - S * vect_mult_K(i, j, k, B_x_out[0][1], B_y_out[0][1], \
			//			B_z_out[0][1]).y - 4 * PI * S / (K(i, j, k) * c) * J_y_out[0][0] + (1 - C) * K_norm(i, j, k).y * scalar_mult_K(i, j, k, \
			//				E_x_out[0][0], E_y_out[0][0], E_z_out[0][0]) + 4 * PI * K_norm(i, j, k).y * scalar_mult_K(i, j, k, \
			//					J_x_out[0][0], J_y_out[0][0], J_z_out[0][0]) * (S / (K(i, j, k) * c) - steps.dt);
			//		E_y_out[0][1] = C * E_y_out[0][1] + S * vect_mult_K(i, j, k, B_x_out[0][0], B_y_out[0][0], \
			//			B_z_out[0][0]).y - 4 * PI * S / (K(i, j, k) * c) * J_y_out[0][1] + (1 - C) * K_norm(i, j, k).y * scalar_mult_K(i, j, k, \
			//				E_x_out[0][1], E_y_out[0][1], E_z_out[0][1]) + 4 * PI * K_norm(i, j, k).y * scalar_mult_K(i, j, k, \
			//					J_x_out[0][1], J_y_out[0][1], J_z_out[0][1]) * (S / (K(i, j, k) * c) - steps.dt);
			//		E_z_out[0][0] = C * E_z_out[0][0] - S * vect_mult_K(i, j, k, B_x_out[0][1], B_y_out[0][1], \
			//			B_z_out[0][1]).z - 4 * PI * S / (K(i, j, k) * c) * J_z_out[0][0] + (1 - C) * K_norm(i, j, k).z * scalar_mult_K(i, j, k, \
			//				E_x_out[0][0], E_y_out[0][0], E_z_out[0][0]) + 4 * PI * K_norm(i, j, k).z * scalar_mult_K(i, j, k, \
			//					J_x_out[0][0], J_y_out[0][0], J_z_out[0][0]) * (S / (K(i, j, k) * c) - steps.dt);
			//		E_z_out[0][1] = C * E_z_out[0][1] + S * vect_mult_K(i, j, k, B_x_out[0][0], B_y_out[0][0], \
			//			B_z_out[0][0]).x - 4 * PI * S / (K(i, j, k) * c) * J_z_out[0][1] + (1 - C) * K_norm(i, j, k).z * scalar_mult_K(i, j, k, \
			//				E_x_out[0][1], E_y_out[0][1], E_z_out[0][1]) + 4 * PI * K_norm(i, j, k).z * scalar_mult_K(i, j, k, \
			//					J_x_out[0][1], J_y_out[0][1], J_z_out[0][1]) * (S / (K(i, j, k) * c) - steps.dt);*/
			//		Ex_out[k][0] += - c * steps.dt * vect_mult_K(i, j, k, Bx_out[k][1], By_out[k][1], \
			//			Bz_out[k][1]).x - 4 * PI * S / (K(i, j, k) * c) * Jx_out[k][0];
			//		Ex_out[k][1] += c * steps.dt * vect_mult_K(i, j, k, Bx_out[k][0], By_out[k][0], \
			//			Bz_out[k][0]).x - 4 * PI * steps.dt * Jx_out[k][1];
			//		Ey_out[k][0] += - c * steps.dt * vect_mult_K(i, j, k, Bx_out[k][1], By_out[k][1], \
			//			Bz_out[k][1]).y - 4 * PI * steps.dt * Jy_out[k][0];
			//		Ey_out[k][1] += c * steps.dt * vect_mult_K(i, j, k, Bx_out[k][0], By_out[k][0], \
			//			Bz_out[k][0]).y - 4 * PI * steps.dt * Jy_out[k][1];
			//		Ez_out[k][0] += - c * steps.dt * vect_mult_K(i, j, k, Bx_out[k][1], By_out[k][1], \
			//			Bz_out[k][1]).z - 4 * PI * steps.dt * Jz_out[k][0];
			//		Ez_out[k][1] +=  c * steps.dt * vect_mult_K(i, j, k, Bx_out[k][0], By_out[k][0], \
			//			Bz_out[k][0]).x - 4 * PI * steps.dt * Jz_out[k][1];
		return (*this);

	}

	//elec_magn_field& PSATD() {
	//	// обновление B
	//	for (int i = 0; i < nx; i++)
	//		for (int j = 0; j < ny; j++)
	//			{
	//				fftw_complex *B_x_out;
	//				fftw_complex *B_y_out;
	//				fftw_complex *B_z_out;
	//				fftw_complex *E_x_out;
	//				fftw_complex *E_y_out;
	//				fftw_complex *E_z_out;
	//				fftw_complex *J_x_out;
	//				fftw_complex *J_y_out;
	//				fftw_complex *J_z_out;
	//				B_x_out = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * nz);
	//				B_y_out = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * nz);
	//				B_z_out = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * nz);
	//				E_x_out = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * nz);
	//				E_y_out = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * nz);
	//				E_z_out = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * nz);
	//				J_x_out = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * nz);
	//				J_y_out = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * nz);
	//				J_z_out = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * nz);
	//				fftw_plan p;
	//				p = fftw_plan_dft_r2c_1d(nz, B.x[i][j].data(), B_x_out, FFTW_ESTIMATE);
	//				fftw_execute(p);
	//				p = fftw_plan_dft_r2c_1d(nz, B.y[i][j].data(), B_y_out, FFTW_ESTIMATE);
	//				fftw_execute(p);
	//				p = fftw_plan_dft_r2c_1d(nz, B.z[i][j].data(), B_z_out, FFTW_ESTIMATE);
	//				fftw_execute(p);
	//				p = fftw_plan_dft_r2c_1d(nz, E.x[i][j].data(), E_x_out, FFTW_ESTIMATE);
	//				fftw_execute(p);
	//				p = fftw_plan_dft_r2c_1d(nz, E.y[i][j].data(), E_y_out, FFTW_ESTIMATE);
	//				fftw_execute(p);
	//				p = fftw_plan_dft_r2c_1d(nz, E.z[i][j].data(), E_z_out, FFTW_ESTIMATE);
	//				fftw_execute(p);
	//				p = fftw_plan_dft_r2c_1d(nz, J.x[i][j].data(), J_x_out, FFTW_ESTIMATE);
	//				fftw_execute(p);
	//				p = fftw_plan_dft_r2c_1d(nz, J.y[i][j].data(), J_y_out, FFTW_ESTIMATE);
	//				fftw_execute(p);
	//				p = fftw_plan_dft_r2c_1d(nz, J.z[i][j].data(), J_z_out, FFTW_ESTIMATE);
	//				fftw_execute(p);
	//				for (int k = 0; k < nz; k++) {
	//					double C = std::cos(K(i, j, k) * c * steps.dt / 2.0);
	//					double S = std::sin(K(i, j, k) * c * steps.dt / 2.0);
	//					/*B_x_out[0][0] = C * B_x_out[0][0] + S * vect_mult_K(i, j, k, E_x_out[0][1], E_y_out[0][1], \
	//						E_z_out[0][1]).x - 4 * PI * (1 - C) / (K(i, j, k) * c) * vect_mult_K(i, j, k, J_x_out[0][1], J_y_out[0][1], \
	//							J_z_out[0][1]).x;
	//					B_x_out[0][1] = C * B_x_out[0][1] - S * vect_mult_K(i, j, k, E_x_out[0][0], E_y_out[0][0], \
	//						E_z_out[0][0]).x + 4 * PI * (1 - C) / (K(i, j, k) * c) * vect_mult_K(i, j, k, J_x_out[0][0], J_y_out[0][0], \
	//							J_z_out[0][0]).x;
	//					B_y_out[0][0] = C * B_y_out[0][0] + S * vect_mult_K(i, j, k, E_x_out[0][1], E_y_out[0][1], \
	//						E_z_out[0][1]).y - 4 * PI * (1 - C) / (K(i, j, k) * c) * vect_mult_K(i, j, k, J_x_out[0][1], J_y_out[0][1], \
	//							J_z_out[0][1]).y;
	//					B_y_out[0][1] = C * B_y_out[0][1] - S * vect_mult_K(i, j, k, E_x_out[0][0], E_y_out[0][0], \
	//						E_z_out[0][0]).y + 4 * PI * (1 - C) / (K(i, j, k) * c) * vect_mult_K(i, j, k, J_x_out[0][0], J_y_out[0][0], \
	//							J_z_out[0][0]).y;
	//					B_z_out[0][0] = C * B_z_out[0][0] + S * vect_mult_K(i, j, k, E_x_out[0][1], E_y_out[0][1], \
	//						E_z_out[0][1]).z - 4 * PI * (1 - C) / (K(i, j, k) * c) * vect_mult_K(i, j, k, J_x_out[0][1], J_y_out[0][1], \
	//							J_z_out[0][1]).z;
	//					B_z_out[0][1] = C * B_z_out[0][1] - S * vect_mult_K(i, j, k, E_x_out[0][0], E_y_out[0][0], \
	//						E_z_out[0][0]).z + 4 * PI * (1 - C) / (K(i, j, k) * c) * vect_mult_K(i, j, k, J_x_out[0][0], J_y_out[0][0], \
	//							J_z_out[0][0]).z;*/
	//							// PSTD
	//					B_x_out[0][0] = B_x_out[0][0] + c * steps.dt * vect_mult_K(i, j, k, E_x_out[0][1], E_y_out[0][1], E_z_out[0][1]).x;
	//					B_x_out[0][1] = B_x_out[0][1] - c * steps.dt * vect_mult_K(i, j, k, E_x_out[0][0], E_y_out[0][0], E_z_out[0][0]).x;
	//					B_y_out[0][0] = B_y_out[0][0] + c * steps.dt * vect_mult_K(i, j, k, E_x_out[0][1], E_y_out[0][1], E_z_out[0][1]).y;
	//					B_y_out[0][1] = B_y_out[0][1] - c * steps.dt * vect_mult_K(i, j, k, E_x_out[0][0], E_y_out[0][0], E_z_out[0][0]).y;
	//					B_z_out[0][0] = B_z_out[0][0] + c * steps.dt * vect_mult_K(i, j, k, E_x_out[0][1], E_y_out[0][1], E_z_out[0][1]).z;
	//					B_z_out[0][1] = B_z_out[0][1] - c * steps.dt * vect_mult_K(i, j, k, E_x_out[0][0], E_y_out[0][0], E_z_out[0][0]).z;
	//				}
	//				p = fftw_plan_dft_c2r_1d(nz, B_x_out, B.x[i][j].data(), FFTW_ESTIMATE);
	//				fftw_execute(p);
	//				p = fftw_plan_dft_c2r_1d(nz, B_y_out, B.y[i][j].data(), FFTW_ESTIMATE);
	//				fftw_execute(p);
	//				p = fftw_plan_dft_c2r_1d(nz, B_z_out, B.z[i][j].data(), FFTW_ESTIMATE);
	//				fftw_execute(p);
	//				p = fftw_plan_dft_c2r_1d(nz, E_x_out, E.x[i][j].data(), FFTW_ESTIMATE);
	//				fftw_execute(p);
	//				p = fftw_plan_dft_c2r_1d(nz, E_y_out, E.y[i][j].data(), FFTW_ESTIMATE);
	//				fftw_execute(p);
	//				p = fftw_plan_dft_c2r_1d(nz, E_z_out, E.z[i][j].data(), FFTW_ESTIMATE);
	//				fftw_execute(p);
	//				fftw_destroy_plan(p);
	//				fftw_free(E_x_out);
	//				fftw_free(E_y_out);
	//				fftw_free(E_z_out);
	//				fftw_free(B_x_out);
	//				fftw_free(B_y_out);
	//				fftw_free(B_z_out);
	//				fftw_free(J_x_out);
	//				fftw_free(J_y_out);
	//				fftw_free(J_z_out);
	//			}
	//	// обновление E
	//	for (int i = 0; i < nx; i++)
	//		for (int j = 0; j < ny; j++)
	//			{
	//				fftw_complex *B_x_out;
	//				fftw_complex *B_y_out;
	//				fftw_complex *B_z_out;
	//				fftw_complex *E_x_out;
	//				fftw_complex *E_y_out;
	//				fftw_complex *E_z_out;
	//				fftw_complex *J_x_out;
	//				fftw_complex *J_y_out;
	//				fftw_complex *J_z_out;
	//				B_x_out = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * nz);
	//				B_y_out = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * nz);
	//				B_z_out = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * nz);
	//				E_x_out = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * nz);
	//				E_y_out = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * nz);
	//				E_z_out = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * nz);
	//				J_x_out = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * nz);
	//				J_y_out = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * nz);
	//				J_z_out = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * nz);
	//				fftw_plan p;
	//				p = fftw_plan_dft_r2c_1d(nz, B.x[i][j].data(), B_x_out, FFTW_ESTIMATE);
	//				fftw_execute(p);
	//				p = fftw_plan_dft_r2c_1d(nz, B.y[i][j].data(), B_y_out, FFTW_ESTIMATE);
	//				fftw_execute(p);
	//				p = fftw_plan_dft_r2c_1d(nz, B.z[i][j].data(), B_z_out, FFTW_ESTIMATE);
	//				fftw_execute(p);
	//				p = fftw_plan_dft_r2c_1d(nz, E.x[i][j].data(), E_x_out, FFTW_ESTIMATE);
	//				fftw_execute(p);
	//				p = fftw_plan_dft_r2c_1d(nz, E.y[i][j].data(), E_y_out, FFTW_ESTIMATE);
	//				fftw_execute(p);
	//				p = fftw_plan_dft_r2c_1d(nz, E.z[i][j].data(), E_z_out, FFTW_ESTIMATE);
	//				fftw_execute(p);
	//				p = fftw_plan_dft_r2c_1d(nz, J.x[i][j].data(), J_x_out, FFTW_ESTIMATE);
	//				fftw_execute(p);
	//				p = fftw_plan_dft_r2c_1d(nz, J.y[i][j].data(), J_y_out, FFTW_ESTIMATE);
	//				fftw_execute(p);
	//				p = fftw_plan_dft_r2c_1d(nz, J.z[i][j].data(), J_z_out, FFTW_ESTIMATE);
	//				fftw_execute(p);
	//				for (int k = 0; k < nz; k++) {
	//					double C = std::cos(K(i, j, k) * c * steps.dt / 2.0);
	//					double S = std::sin(K(i, j, k) * c * steps.dt / 2.0);
	//					/*E_x_out[0][0] = C * E_x_out[0][0] - S * vect_mult_K(i, j, k, B_x_out[0][1], B_y_out[0][1], \
	//						B_z_out[0][1]).x - 4 * PI * S / (K(i, j, k) * c) * J_x_out[0][0] + (1 - C) * K_norm(i, j, k).x * scalar_mult_K(i, j, k, \
	//							E_x_out[0][0], E_y_out[0][0], E_z_out[0][0]) + 4 * PI * K_norm(i, j, k).x * scalar_mult_K(i, j, k, \
	//								J_x_out[0][0], J_y_out[0][0], J_z_out[0][0]) * (S / (K(i, j, k) * c) - steps.dt);
	//					E_x_out[0][1] = C * E_x_out[0][1] + S * vect_mult_K(i, j, k, B_x_out[0][0], B_y_out[0][0], \
	//						B_z_out[0][0]).x - 4 * PI * S / (K(i, j, k) * c) * J_x_out[0][1] + (1 - C) * K_norm(i, j, k).x * scalar_mult_K(i, j, k, \
	//							E_x_out[0][1], E_y_out[0][1], E_z_out[0][1]) + 4 * PI * K_norm(i, j, k).x * scalar_mult_K(i, j, k, \
	//								J_x_out[0][1], J_y_out[0][1], J_z_out[0][1]) * (S / (K(i, j, k) * c) - steps.dt);
	//					E_y_out[0][0] = C * E_y_out[0][0] - S * vect_mult_K(i, j, k, B_x_out[0][1], B_y_out[0][1], \
	//						B_z_out[0][1]).y - 4 * PI * S / (K(i, j, k) * c) * J_y_out[0][0] + (1 - C) * K_norm(i, j, k).y * scalar_mult_K(i, j, k, \
	//							E_x_out[0][0], E_y_out[0][0], E_z_out[0][0]) + 4 * PI * K_norm(i, j, k).y * scalar_mult_K(i, j, k, \
	//								J_x_out[0][0], J_y_out[0][0], J_z_out[0][0]) * (S / (K(i, j, k) * c) - steps.dt);
	//					E_y_out[0][1] = C * E_y_out[0][1] + S * vect_mult_K(i, j, k, B_x_out[0][0], B_y_out[0][0], \
	//						B_z_out[0][0]).y - 4 * PI * S / (K(i, j, k) * c) * J_y_out[0][1] + (1 - C) * K_norm(i, j, k).y * scalar_mult_K(i, j, k, \
	//							E_x_out[0][1], E_y_out[0][1], E_z_out[0][1]) + 4 * PI * K_norm(i, j, k).y * scalar_mult_K(i, j, k, \
	//								J_x_out[0][1], J_y_out[0][1], J_z_out[0][1]) * (S / (K(i, j, k) * c) - steps.dt);
	//					E_z_out[0][0] = C * E_z_out[0][0] - S * vect_mult_K(i, j, k, B_x_out[0][1], B_y_out[0][1], \
	//						B_z_out[0][1]).z - 4 * PI * S / (K(i, j, k) * c) * J_z_out[0][0] + (1 - C) * K_norm(i, j, k).z * scalar_mult_K(i, j, k, \
	//							E_x_out[0][0], E_y_out[0][0], E_z_out[0][0]) + 4 * PI * K_norm(i, j, k).z * scalar_mult_K(i, j, k, \
	//								J_x_out[0][0], J_y_out[0][0], J_z_out[0][0]) * (S / (K(i, j, k) * c) - steps.dt);
	//					E_z_out[0][1] = C * E_z_out[0][1] + S * vect_mult_K(i, j, k, B_x_out[0][0], B_y_out[0][0], \
	//						B_z_out[0][0]).x - 4 * PI * S / (K(i, j, k) * c) * J_z_out[0][1] + (1 - C) * K_norm(i, j, k).z * scalar_mult_K(i, j, k, \
	//							E_x_out[0][1], E_y_out[0][1], E_z_out[0][1]) + 4 * PI * K_norm(i, j, k).z * scalar_mult_K(i, j, k, \
	//								J_x_out[0][1], J_y_out[0][1], J_z_out[0][1]) * (S / (K(i, j, k) * c) - steps.dt);
	//*/
	//// PSTD
	//					E_x_out[k][0] = E_x_out[k][0] - c * steps.dt * vect_mult_K(i, j, k, B_x_out[k][1], B_y_out[k][1], \
	//						B_z_out[k][1]).x - 4 * PI * steps.dt * J_x_out[k][0];
	//					E_x_out[k][1] = E_x_out[k][1] + c * steps.dt * vect_mult_K(i, j, k, B_x_out[k][0], B_y_out[k][0], \
	//						B_z_out[k][0]).x - 4 * PI * steps.dt * J_x_out[k][1];
	//					E_y_out[k][0] = E_y_out[k][0] - c * steps.dt * vect_mult_K(i, j, k, B_x_out[k][1], B_y_out[k][1], \
	//						B_z_out[k][1]).y - 4 * PI * steps.dt * J_y_out[k][0];
	//					E_y_out[k][1] = E_y_out[k][1] + c * steps.dt * vect_mult_K(i, j, k, B_x_out[k][0], B_y_out[k][0], \
	//						B_z_out[k][0]).y - 4 * PI * steps.dt * J_y_out[k][1];
	//					E_z_out[k][0] = E_z_out[k][0] - c * steps.dt * vect_mult_K(i, j, k, B_x_out[k][1], B_y_out[k][1], \
	//						B_z_out[k][1]).z - 4 * PI * steps.dt * J_z_out[k][0];
	//					E_z_out[k][1] = E_z_out[k][1] + c * steps.dt * vect_mult_K(i, j, k, B_x_out[k][0], B_y_out[k][0], \
	//						B_z_out[k][0]).x - 4 * PI * steps.dt * J_z_out[k][1];
	//				}
	//				p = fftw_plan_dft_c2r_1d(nz, B_x_out, B.x[i][j].data(), FFTW_ESTIMATE);
	//				fftw_execute(p);
	//				p = fftw_plan_dft_c2r_1d(nz, B_y_out, B.y[i][j].data(), FFTW_ESTIMATE);
	//				fftw_execute(p);
	//				p = fftw_plan_dft_c2r_1d(nz, B_z_out, B.z[i][j].data(), FFTW_ESTIMATE);
	//				fftw_execute(p);
	//				p = fftw_plan_dft_c2r_1d(nz, E_x_out, E.x[i][j].data(), FFTW_ESTIMATE);
	//				fftw_execute(p);
	//				p = fftw_plan_dft_c2r_1d(nz, E_y_out, E.y[i][j].data(), FFTW_ESTIMATE);
	//				fftw_execute(p);
	//				p = fftw_plan_dft_c2r_1d(nz, E_z_out, E.z[i][j].data(), FFTW_ESTIMATE);
	//				fftw_execute(p);
	//				fftw_destroy_plan(p);
	//				fftw_free(E_x_out);
	//				fftw_free(E_y_out);
	//				fftw_free(E_z_out);
	//				fftw_free(B_x_out);
	//				fftw_free(B_y_out);
	//				fftw_free(B_z_out);
	//				fftw_free(J_x_out);
	//				fftw_free(J_y_out);
	//				fftw_free(J_z_out);
	//			}
	//	return (*this);
	//}
};