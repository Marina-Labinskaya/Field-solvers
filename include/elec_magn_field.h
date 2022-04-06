#pragma once

#include <cmath>
#include <iostream>
#include <vector>
#include <fstream>
#include <complex>
#include <fftw3_mkl.h>

const double PI = 3.1415926536;
const double c = 3e+10;
const int n = 100;


template <typename T = double>
class Array {
public:
	std::vector<T> data_arr;
	int nx, ny, nz, size;

	Array() {
		nx = 0;
		ny = 0;
		nz = 0;
		size = 0;
	}

	Array(int _nx, int _ny, int _nz) {
		resize(_nx, _ny, _nz);
	}

	void resize(int _nx, int _ny, int _nz) {
		nx = _nx;
		ny = _ny;
		nz = _nz;
		size = nx * ny * nz;
		data_arr.assign(size, 0.0);
	}

	~Array() {
		data_arr.clear();
	}

	T& operator() (int i, int j, int k) {
		return data_arr[i * ny * nz + j * nz + k];
	}

	T operator() (int i, int j, int k) const {
		return data_arr[i * ny * nz + j * nz + k];
	}

	T& Array::operator[](int i) {
		return data_arr[i];
	}

	T* data() {
		return data_arr.data();
	}
};

template<typename T = double>
struct field_characteristics {
	Array<T> x;
	Array<T> y;
	Array<T> z;
};

struct estimated_area {
	const double ax;
	const double bx;
	const double ay;
	const double by;
	const double az;
	const double bz;
	estimated_area(double _ax, double _bx, double _ay, double _by, double _az, double _bz) : ax{ _ax }, bx{ _bx }, ay{ _ay }, by{ _by }, az{ _az }, bz{ _bz } {}
};

struct point {
	double x;
	double y;
	double z;
};

struct field_steps {
	const double dx;
	const double dy;
	const double dz;
	const double dt;
	field_steps(double _dx, double _dy, double _dz, double _dt) : dx{ _dx }, dy{ _dy }, dz{ _dz }, dt{ _dt } {}
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

	const int nx;
	const int ny;
	const int nz;
	const int size;

	elec_magn_field(double _dt = 0.1, double _dx = c, double _dy = c, double _dz = c, double _ax = - n * c / 2, double _bx = n * c / 2, double _ay = - n * c / 2, double _by = n * c / 2,
		double _az = -n * c / 2, double _bz = n * c / 2) : area{ _ax, _bx, _ay, _by, _az, _bz }, steps{ _dx, _dy, _dz, _dt }, nx{int((_bx - _ax)/_dx) + 1}, ny{ int((_by - _ay) / _dy) + 1 },
		nz{ int((_bz - _az) / _dz) + 1}, size{nx * ny * nz} {

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
					E.z(0, j, k) = E.z(0, j, k) - 4.0 * PI * t * J.z(0, j, k) + c * t * ((B.y(0, j, k) - B.y(nx - 1, j, k)) / (2.0 * steps.dx) - \
						(B.x(0, j, k) - B.x(0, j - 1, k)) / (2.0 * steps.dy));
					if (k != 0) {
						E.x(0, j, k) = E.x(0, j, k) - 4.0 * PI * t * J.x(0, j, k) + c * t * ((B.z(0, j, k) - B.z(0, j - 1, k)) / (2.0 * steps.dy) - \
							(B.y(0, j, k) - B.y(0, j, k - 1)) / (2.0 * steps.dz));
						E.y(0, j, k) = E.y(0, j, k) - 4.0 * PI * t * J.y(0, j, k) + c * t * ((B.x(0, j, k) - B.x(0, j, k - 1)) / (2.0 * steps.dz) - \
							(B.z(0, j, k) - B.z(nx - 1, j, k)) / (2.0 * steps.dx));
					}
					if (k == 0) {
						E.x(0, j, k) = E.x(0, j, k) - 4.0 * PI * t * J.x(0, j, k) + c * t * ((B.z(0, j, k) - B.z(0, j - 1, k)) / (2.0 * steps.dy) - \
							(B.y(0, j, k) - B.y(0, j, nz - 1)) / (2.0 * steps.dz));
						E.y(0, j, k) = E.y(0, j, k) - 4.0 * PI * t * J.y(0, j, k) + c * t * ((B.x(0, j, k) - B.x(0, j, nz - 1)) / (2.0 * steps.dz) - \
							(B.z(0, j, k) - B.z(nx - 1, j, k)) / (2.0 * steps.dx));
					}
				}
			}
			if (j == 0) {
				for (int k = 0; k < nz; ++k) {
					E.z(0, j, k) = E.z(0, j, k) - 4.0 * PI * t * J.z(0, j, k) + c * t * ((B.y(0, j, k) - B.y(nx - 1, j, k)) / (2.0 * steps.dx) - \
						(B.x(0, j, k) - B.x(0, ny - 1, k)) / (2.0 * steps.dy));
					if (k != 0) {
						E.x(0, j, k) = E.x(0, j, k) - 4.0 * PI * t * J.x(0, j, k) + c * t * ((B.z(0, j, k) - B.z(0, ny - 1, k)) / (2.0 * steps.dy) - \
							(B.y(0, j, k) - B.y(0, j, k - 1)) / (2.0 * steps.dz));
						E.y(0, j, k) = E.y(0, j, k) - 4.0 * PI * t * J.y(0, j, k) + c * t * ((B.x(0, j, k) - B.x(0, j, k - 1)) / (2.0 * steps.dz) - \
							(B.z(0, j, k) - B.z(nx - 1, j, k)) / (2.0 * steps.dx));
					}
					if (k == 0) {
						E.x(0, j, k) = E.x(0, j, k) - 4.0 * PI * t * J.x(0, j, k) + c * t * ((B.z(0, j, k) - B.z(0, ny - 1, k)) / (2.0 * steps.dy) - \
							(B.y(0, j, k) - B.y(0, j, nz - 1)) / (2.0 * steps.dz));
						E.y(0, j, k) = E.y(0, j, k) - 4.0 * PI * t * J.y(0, j, k) + c * t * ((B.x(0, j, k) - B.x(0, j, nz - 1)) / (2.0 * steps.dz) - \
							(B.z(0, j, k) - B.z(nx - 1, j, k)) / (2.0 * steps.dx));
					}
				}
			}
		}

		// j == 0
		for (int i = 0; i < nx; ++i) {
			if (i != 0) {
				for (int k = 0; k < nz; ++k) {
					E.z(i, 0, k) = E.z(i, 0, k) - 4.0 * PI * t * J.z(i, 0, k) + c * t * ((B.y(i, 0, k) - B.y(i - 1, 0, k)) / (2.0 * steps.dx) - \
						(B.x(i, 0, k) - B.x(i, ny - 1, k)) / (2.0 * steps.dy));
					if (k != 0) {
						E.x(i, 0, k) = E.x(i, 0, k) - 4.0 * PI * t * J.x(i, 0, k) + c * t * ((B.z(i, 0, k) - B.z(i, ny - 1, k)) / (2.0 * steps.dy) - \
							(B.y(i, 0, k) - B.y(i, 0, k - 1)) / (2.0 * steps.dz));
						E.y(i, 0, k) = E.y(i, 0, k) - 4.0 * PI * t * J.y(i, 0, k) + c * t * ((B.x(i, 0, k) - B.x(i, 0, k - 1)) / (2.0 * steps.dz) - \
							(B.z(i, 0, k) - B.z(i - 1, 0, k)) / (2.0 * steps.dx));
					}
					if (k == 0) {
						E.x(i, 0, k) = E.x(i, 0, k) - 4.0 * PI * t * J.x(i, 0, k) + c * t * ((B.z(i, 0, k) - B.z(i, ny - 1, k)) / (2.0 * steps.dy) - \
							(B.y(i, 0, k) - B.y(i, 0, nz - 1)) / (2.0 * steps.dz));
						E.y(i, 0, k) = E.y(i, 0, k) - 4.0 * PI * t * J.y(i, 0, k) + c * t * ((B.x(i, 0, k) - B.x(i, 0, nz - 1)) / (2.0 * steps.dz) - \
							(B.z(i, 0, k) - B.z(i - 1, 0, k)) / (2.0 * steps.dx));
					}
				}
			}
			if (i == 0) {
				for (int k = 0; k < nz; ++k) {
					E.z(i, 0, k) = E.z(i, 0, k) - 4.0 * PI * t * J.z(i, 0, k) + c * t * ((B.y(i, 0, k) - B.y(nx - 1, 0, k)) / (2.0 * steps.dx) - \
						(B.x(i, 0, k) - B.x(i, ny - 1, k)) / (2.0 * steps.dy));
					if (k != 0) {
						E.x(i, 0, k) = E.x(i, 0, k) - 4.0 * PI * t * J.x(i, 0, k) + c * t * ((B.z(i, 0, k) - B.z(i, ny - 1, k)) / (2.0 * steps.dy) - \
							(B.y(i, 0, k) - B.y(i, 0, k - 1)) / (2.0 * steps.dz));
						E.y(i, 0, k) = E.y(i, 0, k) - 4.0 * PI * t * J.y(i, 0, k) + c * t * ((B.x(i, 0, k) - B.x(i, 0, k - 1)) / (2.0 * steps.dz) - \
							(B.z(i, 0, k) - B.z(nx - 1, 0, k)) / (2.0 * steps.dx));
					}
					if (k == 0) {
						E.x(i, 0, k) = E.x(i, 0, k) - 4.0 * PI * t * J.x(i, 0, k) + c * t * ((B.z(i, 0, k) - B.z(i, ny - 1, k)) / (2.0 * steps.dy) - \
							(B.y(i, 0, k) - B.y(i, 0, nz - 1)) / (2.0 * steps.dz));
						E.y(i, 0, k) = E.y(i, 0, k) - 4.0 * PI * t * J.y(i, 0, k) + c * t * ((B.x(i, 0, k) - B.x(i, 0, nz - 1)) / (2.0 * steps.dz) - \
							(B.z(i, 0, k) - B.z(nx - 1, 0, k)) / (2.0 * steps.dx));
					}
				}
			}
		}

		// k == 0
		for (int j = 0; j < ny; ++j) {
			if (j != 0) {
				for (int i = 0; i < nx; ++i) {
					E.x(i, j, 0) = E.x(i, j, 0) - 4.0 * PI * t * J.x(i, j, 0) + c * t * ((B.z(i, j, 0) - B.z(i, j - 1, 0)) / (2.0 * steps.dy) - \
						(B.y(i, j, 0) - B.y(i, j, nz - 1)) / (2.0 * steps.dz));
					if (i != 0) {
						E.y(i, j, 0) = E.y(i, j, 0) - 4.0 * PI * t * J.y(i, j, 0) + c * t * ((B.x(i, j, 0) - B.x(i, j, nz - 1)) / (2.0 * steps.dz) - \
							(B.z(i, j, 0) - B.z(i - 1, j, 0)) / (2.0 * steps.dx));
						E.z(i, j, 0) = E.z(i, j, 0) - 4.0 * PI * t * J.z(i, j, 0) + c * t * ((B.y(i, j, 0) - B.y(i - 1, j, 0)) / (2.0 * steps.dx) - \
							(B.x(i, j, 0) - B.x(i, j - 1, 0)) / (2.0 * steps.dy));
					}
					if (i == 0) {
						E.y(i, j, 0) = E.y(i, j, 0) - 4.0 * PI * t * J.y(i, j, 0) + c * t * ((B.x(i, j, 0) - B.x(i, j, nz - 1)) / (2.0 * steps.dz) - \
							(B.z(i, j, 0) - B.z(nx - 1, j, 0)) / (2.0 * steps.dx));
						E.z(i, j, 0) = E.z(i, j, 0) - 4.0 * PI * t * J.z(i, j, 0) + c * t * ((B.y(i, j, 0) - B.y(nx - 1, j, 0)) / (2.0 * steps.dx) - \
							(B.x(i, j, 0) - B.x(i, j - 1, 0)) / (2.0 * steps.dy));
					}
				}
			}

			if (j == 0) {
				for (int i = 0; i < nx; ++i) {
					E.x(i, j, 0) = E.x(i, j, 0) - 4.0 * PI * t * J.x(i, j, 0) + c * t * ((B.z(i, j, 0) - B.z(i, ny - 1, 0)) / (2.0 * steps.dy) - \
						(B.y(i, j, 0) - B.y(i, j, nz - 1)) / (2.0 * steps.dz));
					if (i != 0) {
						E.y(i, j, 0) = E.y(i, j, 0) - 4.0 * PI * t * J.y(i, j, 0) + c * t * ((B.x(i, j, 0) - B.x(i, j, nz - 1)) / (2.0 * steps.dz) - \
							(B.z(i, j, 0) - B.z(i - 1, j, 0)) / (2.0 * steps.dx));
						E.z(i, j, 0) = E.z(i, j, 0) - 4.0 * PI * t * J.z(i, j, 0) + c * t * ((B.y(i, j, 0) - B.y(i - 1, j, 0)) / (2.0 * steps.dx) - \
							(B.x(i, j, 0) - B.x(i, ny - 1, 0)) / (2.0 * steps.dy));
					}
					if (i == 0) {
						E.y(i, j, 0) = E.y(i, j, 0) - 4.0 * PI * t * J.y(i, j, 0) + c * t * ((B.x(i, j, 0) - B.x(i, j, nz - 1)) / (2.0 * steps.dz) - \
							(B.z(i, j, 0) - B.z(nx - 1, j, 0)) / (2.0 * steps.dx));
						E.z(i, j, 0) = E.z(i, j, 0) - 4.0 * PI * t * J.z(i, j, 0) + c * t * ((B.y(i, j, 0) - B.y(nx - 1, j, 0)) / (2.0 * steps.dx) - \
							(B.x(i, j, 0) - B.x(i, ny - 1, 0)) / (2.0 * steps.dy));
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
					B.z(nx - 1, j, k) = B.z(nx - 1, j, k) + c * t * ((E.x(nx - 1, j + 1, k) - E.x(nx - 1, j, k)) / (2.0 * steps.dy) - \
						(E.y(0, j, k) - E.y(nx - 1, j, k)) / (2.0 * steps.dx));
					if (k != nz - 1) {
						B.x(nx - 1, j, k) = B.x(nx - 1, j, k) + c * t * ((E.y(nx - 1, j, k + 1) - E.y(nx - 1, j, k)) / (2.0 * steps.dz) - \
							(E.z(nx - 1, j + 1, k) - E.z(nx - 1, j, k)) / (2.0 * steps.dy));
						B.y(nx - 1, j, k) = B.y(nx - 1, j, k) + c * t * ((E.z(0, j, k) - E.z(nx - 1, j, k)) / (2.0 * steps.dx) - \
							(E.x(nx - 1, j, k + 1) - E.x(nx - 1, j, k)) / (2.0 * steps.dz));
					}
					if (k == nz - 1) {
						B.x(nx - 1, j, k) = B.x(nx - 1, j, k) + c * t * ((E.y(nx - 1, j, 0) - E.y(nx - 1, j, k)) / (2.0 * steps.dz) - \
							(E.z(nx - 1, j + 1, k) - E.z(nx - 1, j, k)) / (2.0 * steps.dy));
						B.y(nx - 1, j, k) = B.y(nx - 1, j, k) + c * t * ((E.z(0, j, k) - E.z(nx - 1, j, k)) / (2.0 * steps.dx) - \
							(E.x(nx - 1, j, 0) - E.x(nx - 1, j, k)) / (2.0 * steps.dz));
					}
				}
			}
			if (j == ny - 1) {
				for (int k = 0; k < nz; ++k) {
					B.z(nx - 1, j, k) = B.z(nx - 1, j, k) + c * t * ((E.x(nx - 1, 0, k) - E.x(nx - 1, j, k)) / (2.0 * steps.dy) - \
						(E.y(0, j, k) - E.y(nx - 1, j, k)) / (2.0 * steps.dx));
					if (k != nz - 1) {
						B.x(nx - 1, j, k) = B.x(nx - 1, j, k) + c * t * ((E.y(nx - 1, j, k + 1) - E.y(nx - 1, j, k)) / (2.0 * steps.dz) - \
							(E.z(nx - 1, 0, k) - E.z(nx - 1, j, k)) / (2.0 * steps.dy));
						B.y(nx - 1, j, k) = B.y(nx - 1, j, k) + c * t * ((E.z(0, j, k) - E.z(nx - 1, j, k)) / (2.0 * steps.dx) - \
							(E.x(nx - 1, j, k + 1) - E.x(nx - 1, j, k)) / (2.0 * steps.dz));
					}
					if (k == nz - 1) {
						B.x(nx - 1, j, k) = B.x(nx - 1, j, k) + c * t * ((E.y(nx - 1, j, 0) - E.y(nx - 1, j, k)) / (2.0 * steps.dz) - \
							(E.z(nx - 1, 0, k) - E.z(nx - 1, j, k)) / (2.0 * steps.dy));
						B.y(nx - 1, j, k) = B.y(nx - 1, j, k) + c * t * ((E.z(0, j, k) - E.z(nx - 1, j, k)) / (2.0 * steps.dx) - \
							(E.x(nx - 1, j, 0) - E.x(nx - 1, j, k)) / (2.0 * steps.dz));
					}
				}
			}
		}

		// j == ny - 1
		for (int i = 0; i < nx; ++i) {
			if (i != nx - 1) {
				for (int k = 0; k < nz; ++k) {
					B.z(i, ny - 1, k) = B.z(i, ny - 1, k) + c * t * ((E.x(i, 0, k) - E.x(i, ny - 1, k)) / (2.0 * steps.dy) - \
						(E.y(i + 1, ny - 1, k) - E.y(i, ny - 1, k)) / (2.0 * steps.dx));
					if (k != nz - 1) {
						B.x(i, ny - 1, k) = B.x(i, ny - 1, k) + c * t * ((E.y(i, ny - 1, k + 1) - E.y(i, ny - 1, k)) / (2.0 * steps.dz) - \
							(E.z(i, 0, k) - E.z(i, ny - 1, k)) / (2.0 * steps.dy));
						B.y(i, ny - 1, k) = B.y(i, ny - 1, k) + c * t * ((E.z(i + 1, ny - 1, k) - E.z(i, ny - 1, k)) / (2.0 * steps.dx) - \
							(E.x(i, ny - 1, k + 1) - E.x(i, ny - 1, k)) / (2.0 * steps.dz));
					}
					if (k == nz - 1) {
						B.x(i, ny - 1, k) = B.x(i, ny - 1, k) + c * t * ((E.y(i, ny - 1, 0) - E.y(i, ny - 1, k)) / (2.0 * steps.dz) - \
							(E.z(i, 0, k) - E.z(i, ny - 1, k)) / (2.0 * steps.dy));
						B.y(i, ny - 1, k) = B.y(i, ny - 1, k) + c * t * ((E.z(i + 1, ny - 1, k) - E.z(i, ny - 1, k)) / (2.0 * steps.dx) - \
							(E.x(i, ny - 1, 0) - E.x(i, ny - 1, k)) / (2.0 * steps.dz));
					}
				}
			}
			if (i == nx - 1) {
				for (int k = 0; k < nz; ++k) {
					B.z(i, ny - 1, k) = B.z(i, ny - 1, k) + c * t * ((E.x(i, 0, k) - E.x(i, ny - 1, k)) / (2.0 * steps.dy) - \
						(E.y(0, ny - 1, k) - E.y(i, ny - 1, k)) / (2.0 * steps.dx));
					if (k != nz - 1) {
						B.x(i, ny - 1, k) = B.x(i, ny - 1, k) + c * t * ((E.y(i, ny - 1, k + 1) - E.y(i, ny - 1, k)) / (2.0 * steps.dz) - \
							(E.z(i, 0, k) - E.z(i, ny - 1, k)) / (2.0 * steps.dy));
						B.y(i, ny - 1, k) = B.y(i, ny - 1, k) + c * t * ((E.z(0, ny - 1, k) - E.z(i, ny - 1, k)) / (2.0 * steps.dx) - \
							(E.x(i, ny - 1, k + 1) - E.x(i, ny - 1, k)) / (2.0 * steps.dz));
					}
					if (k == nz - 1) {
						B.x(i, ny - 1, k) = B.x(i, ny - 1, k) + c * t * ((E.y(i, ny - 1, 0) - E.y(i, ny - 1, k)) / (2.0 * steps.dz) - \
							(E.z(i, 0, k) - E.z(i, ny - 1, k)) / (2.0 * steps.dy));
						B.y(i, ny - 1, k) = B.y(i, ny - 1, k) + c * t * ((E.z(0, ny - 1, k) - E.z(i, ny - 1, k)) / (2.0 * steps.dx) - \
							(E.x(i, ny - 1, 0) - E.x(i, ny - 1, k)) / (2.0 * steps.dz));
					}
				}
			}
		}


		// k == nz - 1
		for (int i = 0; i < nx; ++i) {
			if (i != nx - 1) {
				for (int j = 0; j < ny; ++j) {
					B.y(i, j, nz - 1) = B.y(i, j, nz - 1) + c * t * ((E.z(i + 1, j, nz - 1) - E.z(i, j, nz - 1)) / (2.0 * steps.dx) - \
						(E.x(i, j, 0) - E.x(i, j, nz - 1)) / (2.0 * steps.dz));
					if (j != ny - 1) {
						B.x(i, j, nz - 1) = B.x(i, j, nz - 1) + c * t * ((E.y(i, j, 0) - E.y(i, j, nz - 1)) / (2.0 * steps.dz) - \
							(E.z(i, j + 1, nz - 1) - E.z(i, j, nz - 1)) / (2.0 * steps.dy));
						B.z(i, j, nz - 1) = B.z(i, j, nz - 1) + c * t * ((E.x(i, j + 1, nz - 1) - E.x(i, j, nz - 1)) / (2.0 * steps.dy) - \
							(E.y(i + 1, j, nz - 1) - E.y(i, j, nz - 1)) / (2.0 * steps.dx));
					}
					if (j == ny - 1) {
						B.x(i, j, nz - 1) = B.x(i, j, nz - 1) + c * t * ((E.y(i, j, 0) - E.y(i, j, nz - 1)) / (2.0 * steps.dz) - \
							(E.z(i, 0, nz - 1) - E.z(i, j, nz - 1)) / (2.0 * steps.dy));
						B.z(i, j, nz - 1) = B.z(i, j, nz - 1) + c * t * ((E.x(i, 0, nz - 1) - E.x(i, j, nz - 1)) / (2.0 * steps.dy) - \
							(E.y(i + 1, j, nz - 1) - E.y(i, j, nz - 1)) / (2.0 * steps.dx));
					}
				}
			}
			if (i == nx - 1) {
				for (int j = 0; j < ny; ++j) {
					B.y(i, j, nz - 1) = B.y(i, j, nz - 1) + c * t * ((E.z(0, j, nz - 1) - E.z(i, j, nz - 1)) / (2.0 * steps.dx) - \
						(E.x(i, j, 0) - E.x(i, j, nz - 1)) / (2.0 * steps.dz));
					if (j != ny - 1) {
						B.x(i, j, nz - 1) = B.x(i, j, nz - 1) + c * t * ((E.y(i, j, 0) - E.y(i, j, nz - 1)) / (2.0 * steps.dz) - \
							(E.z(i, j + 1, nz - 1) - E.z(i, j, nz - 1)) / (2.0 * steps.dy));
						B.z(i, j, nz - 1) = B.z(i, j, nz - 1) + c * t * ((E.x(i, j + 1, nz - 1) - E.x(i, j, nz - 1)) / (2.0 * steps.dy) - \
							(E.y(0, j, nz - 1) - E.y(i, j, nz - 1)) / (2.0 * steps.dx));
					}
					if (j == ny - 1) {
						B.x(i, j, nz - 1) = B.x(i, j, nz - 1) + c * t * ((E.y(i, j, 0) - E.y(i, j, nz - 1)) / (2.0 * steps.dz) - \
							(E.z(i, 0, nz - 1) - E.z(i, j, nz - 1)) / (2.0 * steps.dy));
						B.z(i, j, nz - 1) = B.z(i, j, nz - 1) + c * t * ((E.x(i, 0, nz - 1) - E.x(i, j, nz - 1)) / (2.0 * steps.dy) - \
							(E.y(0, j, nz - 1) - E.y(i, j, nz - 1)) / (2.0 * steps.dx));
					}
				}
			}
		}
		return (*this);
	}


	elec_magn_field& FDTD(double t) {
		for (int i = 0; i < nx - 1; i++)
			for (int j = 0; j < ny - 1; j++)
				for (int k = 0; k < nz - 1; k++) {
					B.x(i, j, k) = B.x(i, j, k) + c * steps.dt * ((E.y(i, j, k + 1) - E.y(i, j, k)) / (2.0 * steps.dz) - \
						(E.z(i, j + 1, k) - E.z(i, j, k)) / (2.0 * steps.dy));
					B.y(i, j, k) = B.y(i, j, k) + c * steps.dt * ((E.z(i + 1, j, k) - E.z(i, j, k)) / (2.0 * steps.dx) - \
						(E.x(i, j, k + 1) - E.x(i, j, k)) / (2.0 * steps.dz));
					B.z(i, j, k) = B.z(i, j, k) + c * steps.dt * ((E.x(i, j + 1, k) - E.x(i, j, k)) / (2.0 * steps.dy) - \
						(E.y(i + 1, j, k) - E.y(i, j, k)) / (2.0 * steps.dx));
				}
		if (t)
			get_boundary_conditions_B(steps.dt);
		else
			get_boundary_conditions_B(steps.dt / 2);

		for (int i = 1; i < nx; i++)
			for (int j = 1; j < ny; j++)
				for (int k = 1; k < nz; k++) {
					E.x(i, j, k) = E.x(i, j, k) - 4.0 * PI * steps.dt * J.x(i, j, k) + c * steps.dt * ((B.z(i, j, k) - B.z(i, j - 1, k)) / (2.0 * steps.dy) - \
						(B.y(i, j, k) - B.y(i, j, k - 1)) / (2.0 * steps.dz));
					E.y(i, j, k) = E.y(i, j, k) - 4.0 * PI * steps.dt * J.y(i, j, k) + c * steps.dt * ((B.x(i, j, k) - B.x(i, j, k - 1)) / (2.0 * steps.dz) - \
						(B.z(i, j, k) - B.z(i - 1, j, k)) / (2.0 * steps.dx));
					E.z(i, j, k) = E.z(i, j, k) - 4.0 * PI * steps.dt * J.z(i, j, k) + c * steps.dt * ((B.y(i, j, k) - B.y(i - 1, j, k)) / (2.0 * steps.dx) - \
						(B.x(i, j, k) - B.x(i, j - 1, k)) / (2.0 * steps.dy));
				}

		get_boundary_conditions_E(steps.dt);

		return (*this);
	}

	elec_magn_field& modify_Jx(int m) {
		double t = steps.dt * (m - 0.5);
		double z, y;
		int x = 0;
		for (int i = 0; i < nx; ++i) {
			for (int j = 0; j < ny; ++j) {
				for (int k = 0; k < nz; ++k) {
					z = area.az + k * steps.dz;
					y = area.ay + j * steps.dy;
					if (nx > 1)
						x = area.ax + i * steps.dx;
					if ((x >= -Tx / 4.0) && (x <= Tx / 4.0) && (y >= -Ty / 4.0) && (y <= Ty / 4.0) && (z >= -Tz / 4.0) && (z <= Tz / 4.0)) {
						J.x(0, j, k) = sin(2 * PI * t / T) * pow(cos(2 * PI * x / Tx), 2) * pow(cos(2 * PI * y / Ty), 2) * pow(cos(2 * PI * z / Tz), 2);
					}
				}
			}
		}
		return (*this);
	}

	elec_magn_field& modify_Jy(int m) {
		double t = steps.dt * (m - 0.5);
		double x, z;
		int y = 0;
		for (int i = 0; i < nx; ++i)
			for (int j = 0; j < ny; ++j) {
				for (int k = 0; k < nz; ++k) {
					x = area.ax + i * steps.dx;
					z = area.az + k * steps.dz;
					if (ny > 1)
						y = area.ay + j * steps.dy;
					if ((x >= -Tx / 4.0) && (x <= Tx / 4.0) && (y >= -Ty / 4.0) && (y <= Ty / 4.0) && (z >= -Tz / 4.0) && (z <= Tz / 4.0)) {
						J.y(i, 0, k) = sin(2 * PI * t / T) * pow(cos(2 * PI * x / Tx), 2) * pow(cos(2 * PI * y / Ty), 2) * pow(cos(2 * PI * z / Tz), 2);
					}
				}
			}
		return (*this);
	}

	elec_magn_field& modify_Jz(int m) {
		double t = steps.dt * (m - 0.5);
		double x, y;
		int z = 0;
		for (int i = 0; i < nx; ++i)
			for (int j = 0; j < ny; ++j) {
				for (int k = 0; k < nz; ++k) {
					x = area.ax + i * steps.dx;
					y = area.ay + j * steps.dy;
					if (nz > 1)
						z = area.az + k * steps.dz;
					if ((x >= -Tx / 4.0) && (x <= Tx / 4.0) && (y >= -Ty / 4.0) && (y <= Ty / 4.0) && (z >= -Tz / 4.0) && (z <= Tz / 4.0)) {
						J.z(i, j, 0) = sin(2 * PI * t / T) * pow(cos(2 * PI * x / Tx), 2) * pow(cos(2 * PI * y / Ty), 2) * pow(cos(2 * PI * z / Tz), 2);
					}
				}
			}
		return (*this);
	}

	elec_magn_field& set_zero_Jx() {
		for (int i = 0; i < nx; ++i)
			for (int j = 0; j < ny; ++j)
				for (int k = 0; k < nz; ++k)
					J.x(i, j, k) = 0.0;
		return (*this);
	}

	elec_magn_field& set_zero_Jy() {
		for (int i = 0; i < nx; ++i)
			for (int j = 0; j < ny; ++j)
				for (int k = 0; k < nz; ++k)
					J.y(i, j, k) = 0.0;
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
			res.x = (2 * PI * i / (area.bx - area.ax));
		else
			res.x = (2 * PI * (i - nx) / (area.bx - area.ax));
		if (j <= ny / 2)
			res.y = (2 * PI * j / (area.by - area.ay));
		else
			res.y = (2 * PI * (j - ny) / (area.by - area.ay));
		if (k <= nz / 2)
			res.z = (2 * PI * k / (area.bz - area.az));
		else
			res.z = (2 * PI * (k - nz) / (area.bz - area.az));

		return res;
	}

	double k_vec_len(int i, int j, int k) {
		return sqrt(w(i, j, k).x * w(i, j, k).x + w(i, j, k).y * w(i, j, k).y + w(i, j, k).z * w(i, j, k).z);
	}

	point k_norm(int i, int j, int k) {
		point res;
		if (i != 0 && j != 0 && k != 0) {
			res.x = w(i, j, k).x / k_vec_len(i, j, k);
			res.y = w(i, j, k).y / k_vec_len(i, j, k);
			res.z = w(i, j, k).z / k_vec_len(i, j, k);
		}
		else
		{
			res.x = 0.0;
			res.y = 0.0;
			res.z = 0.0;
		}

		return res;
	}


	elec_magn_field& PSTD(double t) {
		Array<std::complex<double>> Bx_out(nx, ny, nz / 2 + 1);
		Array<std::complex<double>> By_out(nx, ny, nz / 2 + 1);
		Array<std::complex<double>> Bz_out(nx, ny, nz / 2 + 1);
		Array<std::complex<double>> Ex_out(nx, ny, nz / 2 + 1);
		Array<std::complex<double>> Ey_out(nx, ny, nz / 2 + 1);
		Array<std::complex<double>> Ez_out(nx, ny, nz / 2 + 1);
		Array<std::complex<double>> Jx_out(nx, ny, nz / 2 + 1);
		Array<std::complex<double>> Jy_out(nx, ny, nz / 2 + 1);
		Array<std::complex<double>> Jz_out(nx, ny, nz / 2 + 1);

		fftw_plan Ex_plan, Ey_plan, Ez_plan, Bx_plan, By_plan, Bz_plan, Jx_plan, Jy_plan, Jz_plan;
		fftw_plan _Ex_plan, _Ey_plan, _Ez_plan, _Bx_plan, _By_plan, _Bz_plan, _Jx_plan, _Jy_plan, _Jz_plan;

		Ex_plan = fftw_plan_dft_r2c_3d(nx, ny, nz, E.x.data(), (fftw_complex *)(Ex_out.data()), FFTW_ESTIMATE);
		Ey_plan = fftw_plan_dft_r2c_3d(nx, ny, nz, E.y.data(), (fftw_complex *)(Ey_out.data()), FFTW_ESTIMATE);
		Ez_plan = fftw_plan_dft_r2c_3d(nx, ny, nz, E.z.data(), (fftw_complex *)(Ez_out.data()), FFTW_ESTIMATE);
		Bx_plan = fftw_plan_dft_r2c_3d(nx, ny, nz, B.x.data(), (fftw_complex *)(Bx_out.data()), FFTW_ESTIMATE);
		By_plan = fftw_plan_dft_r2c_3d(nx, ny, nz, B.y.data(), (fftw_complex *)(By_out.data()), FFTW_ESTIMATE);
		Bz_plan = fftw_plan_dft_r2c_3d(nx, ny, nz, B.z.data(), (fftw_complex *)(Bz_out.data()), FFTW_ESTIMATE);
		Jx_plan = fftw_plan_dft_r2c_3d(nx, ny, nz, J.x.data(), (fftw_complex *)(Jx_out.data()), FFTW_ESTIMATE);
		Jy_plan = fftw_plan_dft_r2c_3d(nx, ny, nz, J.y.data(), (fftw_complex *)(Jy_out.data()), FFTW_ESTIMATE);
		Jz_plan = fftw_plan_dft_r2c_3d(nx, ny, nz, J.z.data(), (fftw_complex *)(Jz_out.data()), FFTW_ESTIMATE);

		fftw_execute(Ex_plan);
		fftw_execute(Ey_plan);
		fftw_execute(Ez_plan);
		fftw_execute(Bx_plan);
		fftw_execute(By_plan);
		fftw_execute(Bz_plan);
		fftw_execute(Jx_plan);
		fftw_execute(Jy_plan);
		fftw_execute(Jz_plan);

		double dt_E = steps.dt;
		double dt_B = steps.dt;

		for (int _i = 0; _i < nx; _i++) {
			for (int _j = 0; _j < ny; _j++) {
				for (int _k = 0; _k < nz / 2 + 1; _k++) {

				if (t == 0)
					dt_B = steps.dt / 2;

				Bx_out(_i, _j, _k).real(Bx_out(_i, _j, _k).real() + c * dt_B * (w(_i, _j, _k).y * Ez_out(_i, _j, _k).imag() - w(_i, _j, _k).z * Ey_out(_i, _j, _k).imag()));
				Bx_out(_i, _j, _k).imag(Bx_out(_i, _j, _k).imag() - c * dt_B * (w(_i, _j, _k).y * Ez_out(_i, _j, _k).real() - w(_i, _j, _k).z * Ey_out(_i, _j, _k).real()));
				By_out(_i, _j, _k).real(By_out(_i, _j, _k).real() - c * dt_B * (w(_i, _j, _k).x * Ez_out(_i, _j, _k).imag() - w(_i, _j, _k).z * Ex_out(_i, _j, _k).imag()));
				By_out(_i, _j, _k).imag(By_out(_i, _j, _k).imag() + c * dt_B * (w(_i, _j, _k).x * Ez_out(_i, _j, _k).real() - w(_i, _j, _k).z * Ex_out(_i, _j, _k).real()));
				Bz_out(_i, _j, _k).real(Bz_out(_i, _j, _k).real() + c * dt_B * (w(_i, _j, _k).x * Ey_out(_i, _j, _k).imag() - w(_i, _j, _k).y * Ex_out(_i, _j, _k).imag()));
				Bz_out(_i, _j, _k).imag(Bz_out(_i, _j, _k).imag() - c * dt_B * (w(_i, _j, _k).x * Ey_out(_i, _j, _k).real() - w(_i, _j, _k).y * Ex_out(_i, _j, _k).real()));

				Ex_out(_i, _j, _k).real(Ex_out(_i, _j, _k).real() - c * dt_E * (w(_i, _j, _k).y * Bz_out(_i, _j, _k).imag() - w(_i, _j, _k).z * By_out(_i, _j, _k).imag()) \
					- 4 * PI * dt_E * Jx_out(_i, _j, _k).real());
				Ex_out(_i, _j, _k).imag(Ex_out(_i, _j, _k).imag() + c * dt_E * (w(_i, _j, _k).y * Bz_out(_i, _j, _k).real() - w(_i, _j, _k).z * By_out(_i, _j, _k).real()) \
					- 4 * PI * dt_E * Jx_out(_i, _j, _k).imag());
				Ey_out(_i, _j, _k).real(Ey_out(_i, _j, _k).real() + c * dt_E * (w(_i, _j, _k).x * Bz_out(_i, _j, _k).imag() - w(_i, _j, _k).z * Bx_out(_i, _j, _k).imag()) \
					- 4 * PI * dt_E * Jy_out(_i, _j, _k).real());
				Ey_out(_i, _j, _k).imag(Ey_out(_i, _j, _k).imag() - c * dt_E * (w(_i, _j, _k).x * Bz_out(_i, _j, _k).real() - w(_i, _j, _k).z * Bx_out(_i, _j, _k).real()) \
					- 4 * PI * dt_E * Jy_out(_i, _j, _k).imag());
				Ez_out(_i, _j, _k).real(Ez_out(_i, _j, _k).real() - c * dt_E * (w(_i, _j, _k).x * By_out(_i, _j, _k).imag() - w(_i, _j, _k).y * Bx_out(_i, _j, _k).imag()) \
					- 4 * PI * dt_E * Jz_out(_i, _j, _k).real());
				Ez_out(_i, _j, _k).imag(Ez_out(_i, _j, _k).imag() + c * dt_E * (w(_i, _j, _k).x * By_out(_i, _j, _k).real() - w(_i, _j, _k).y * Bx_out(_i, _j, _k).real()) \
					- 4 * PI * dt_E * Jz_out(_i, _j, _k).imag());

				dt_B = steps.dt;
				}
			}
		}

		_Ex_plan = fftw_plan_dft_c2r_3d(nx, ny, nz, (fftw_complex*)(Ex_out.data()), E.x.data(), FFTW_ESTIMATE);
		_Ey_plan = fftw_plan_dft_c2r_3d(nx, ny, nz, (fftw_complex*)(Ey_out.data()), E.y.data(), FFTW_ESTIMATE);
		_Ez_plan = fftw_plan_dft_c2r_3d(nx, ny, nz, (fftw_complex*)(Ez_out.data()), E.z.data(), FFTW_ESTIMATE);
		_Bx_plan = fftw_plan_dft_c2r_3d(nx, ny, nz, (fftw_complex*)(Bx_out.data()), B.x.data(), FFTW_ESTIMATE);
		_By_plan = fftw_plan_dft_c2r_3d(nx, ny, nz, (fftw_complex*)(By_out.data()), B.y.data(), FFTW_ESTIMATE);
		_Bz_plan = fftw_plan_dft_c2r_3d(nx, ny, nz, (fftw_complex*)(Bz_out.data()), B.z.data(), FFTW_ESTIMATE);
		_Jx_plan = fftw_plan_dft_c2r_3d(nx, ny, nz, (fftw_complex*)(Jx_out.data()), J.x.data(), FFTW_ESTIMATE);
		_Jy_plan = fftw_plan_dft_c2r_3d(nx, ny, nz, (fftw_complex*)(Jy_out.data()), J.y.data(), FFTW_ESTIMATE);
		_Jz_plan = fftw_plan_dft_c2r_3d(nx, ny, nz, (fftw_complex*)(Jz_out.data()), J.z.data(), FFTW_ESTIMATE);

		fftw_execute(_Ex_plan);
		fftw_execute(_Ey_plan);
		fftw_execute(_Ez_plan);
		fftw_execute(_Bx_plan);
		fftw_execute(_By_plan);
		fftw_execute(_Bz_plan);
		fftw_execute(_Jx_plan);
		fftw_execute(_Jy_plan);
		fftw_execute(_Jz_plan);

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

		fftw_destroy_plan(Ex_plan);
		fftw_destroy_plan(Ey_plan);
		fftw_destroy_plan(Ez_plan);
		fftw_destroy_plan(Bx_plan);
		fftw_destroy_plan(By_plan);
		fftw_destroy_plan(Bz_plan);
		fftw_destroy_plan(Jx_plan);
		fftw_destroy_plan(Jy_plan);
		fftw_destroy_plan(Jz_plan);

		fftw_destroy_plan(_Ex_plan);
		fftw_destroy_plan(_Ey_plan);
		fftw_destroy_plan(_Ez_plan);
		fftw_destroy_plan(_Bx_plan);
		fftw_destroy_plan(_By_plan);
		fftw_destroy_plan(_Bz_plan);
		fftw_destroy_plan(_Jx_plan);
		fftw_destroy_plan(_Jy_plan);
		fftw_destroy_plan(_Jz_plan);

		return (*this);
	}

	elec_magn_field& PSATD(double t) {
		Array<std::complex<double>> Bx_out(nx, ny, nz / 2 + 1);
		Array<std::complex<double>> By_out(nx, ny, nz / 2 + 1);
		Array<std::complex<double>> Bz_out(nx, ny, nz / 2 + 1);
		Array<std::complex<double>> Ex_out(nx, ny, nz / 2 + 1);
		Array<std::complex<double>> Ey_out(nx, ny, nz / 2 + 1);
		Array<std::complex<double>> Ez_out(nx, ny, nz / 2 + 1);
		Array<std::complex<double>> Jx_out(nx, ny, nz / 2 + 1);
		Array<std::complex<double>> Jy_out(nx, ny, nz / 2 + 1);
		Array<std::complex<double>> Jz_out(nx, ny, nz / 2 + 1);

		fftw_plan Ex_plan, Ey_plan, Ez_plan, Bx_plan, By_plan, Bz_plan, Jx_plan, Jy_plan, Jz_plan;
		fftw_plan _Ex_plan, _Ey_plan, _Ez_plan, _Bx_plan, _By_plan, _Bz_plan, _Jx_plan, _Jy_plan, _Jz_plan;

		Ex_plan = fftw_plan_dft_r2c_3d(nx, ny, nz, E.x.data(), (fftw_complex *)(Ex_out.data()), FFTW_ESTIMATE);
		Ey_plan = fftw_plan_dft_r2c_3d(nx, ny, nz, E.y.data(), (fftw_complex *)(Ey_out.data()), FFTW_ESTIMATE);
		Ez_plan = fftw_plan_dft_r2c_3d(nx, ny, nz, E.z.data(), (fftw_complex *)(Ez_out.data()), FFTW_ESTIMATE);
		Bx_plan = fftw_plan_dft_r2c_3d(nx, ny, nz, B.x.data(), (fftw_complex *)(Bx_out.data()), FFTW_ESTIMATE);
		By_plan = fftw_plan_dft_r2c_3d(nx, ny, nz, B.y.data(), (fftw_complex *)(By_out.data()), FFTW_ESTIMATE);
		Bz_plan = fftw_plan_dft_r2c_3d(nx, ny, nz, B.z.data(), (fftw_complex *)(Bz_out.data()), FFTW_ESTIMATE);
		Jx_plan = fftw_plan_dft_r2c_3d(nx, ny, nz, J.x.data(), (fftw_complex *)(Jx_out.data()), FFTW_ESTIMATE);
		Jy_plan = fftw_plan_dft_r2c_3d(nx, ny, nz, J.y.data(), (fftw_complex *)(Jy_out.data()), FFTW_ESTIMATE);
		Jz_plan = fftw_plan_dft_r2c_3d(nx, ny, nz, J.z.data(), (fftw_complex *)(Jz_out.data()), FFTW_ESTIMATE);

		fftw_execute(Ex_plan);
		fftw_execute(Ey_plan);
		fftw_execute(Ez_plan);
		fftw_execute(Bx_plan);
		fftw_execute(By_plan);
		fftw_execute(Bz_plan);
		fftw_execute(Jx_plan);
		fftw_execute(Jy_plan);
		fftw_execute(Jz_plan);

		double dt_E = steps.dt;
		double dt_B = steps.dt;
		std::complex<double> I(0.0, 1.0);

		for (int _i = 0; _i < nx; _i++) {
			for (int _j = 0; _j < ny; _j++) {
				for (int _k = 0; _k < nz / 2 + 1; _k++) {
					if (t == 0)
						dt_B = steps.dt / 2.0;

					double C_B = cos(k_vec_len(_i, _j, _k) * c * dt_B / 2.0);
					double S_B = sin(k_vec_len(_i, _j, _k) * c * dt_B / 2.0);
					double C_E = cos(k_vec_len(_i, _j, _k) * c * dt_E / 2.0);
					double S_E = sin(k_vec_len(_i, _j, _k) * c * dt_E / 2.0);

					if (!(_i == 0 && _j == 0 && _k == 0)) {

						/*Bx_out(_i, _j, _k).real(C_B * Bx_out(_i, _j, _k).real() + S_B * (k_norm(_i, _j, _k).y * Ez_out(_i, _j, _k).imag() - Ey_out(_i, _j, _k).imag() * k_norm(_i, _j, _k).z) \
							- 4 * PI * (1 - C_B) / (k_vec_len(_i, _j, _k) * c) * (k_norm(_i, _j, _k).y * Jz_out(_i, _j, _k).imag() - Jy_out(_i, _j, _k).imag() * k_norm(_i, _j, _k).z));

						Bx_out(_i, _j, _k).imag(C_B * Bx_out(_i, _j, _k).imag() - S_B * (k_norm(_i, _j, _k).y * Ez_out(_i, _j, _k).real() - Ey_out(_i, _j, _k).real() * k_norm(_i, _j, _k).z) \
							+ 4 * PI * (1 - C_B) / (k_vec_len(_i, _j, _k) * c) * (k_norm(_i, _j, _k).y * Jz_out(_i, _j, _k).real() - Jy_out(_i, _j, _k).real() * k_norm(_i, _j, _k).z));

						By_out(_i, _j, _k).real(C_B * By_out(_i, _j, _k).real() - S_B * (k_norm(_i, _j, _k).x * Ez_out(_i, _j, _k).imag() - Ex_out(_i, _j, _k).imag() * k_norm(_i, _j, _k).z) \
							+ 4 * PI * (1 - C_B) / (k_vec_len(_i, _j, _k) * c) * (k_norm(_i, _j, _k).x * Jz_out(_i, _j, _k).imag() - Jx_out(_i, _j, _k).imag() * k_norm(_i, _j, _k).z));

						By_out(_i, _j, _k).imag(C_B * By_out(_i, _j, _k).imag() + S_B * (k_norm(_i, _j, _k).x * Ez_out(_i, _j, _k).real() - Ex_out(_i, _j, _k).real() * k_norm(_i, _j, _k).z) \
							- 4 * PI * (1 - C_B) / (k_vec_len(_i, _j, _k) * c) * (k_norm(_i, _j, _k).x * Jz_out(_i, _j, _k).real() - Jx_out(_i, _j, _k).real() * k_norm(_i, _j, _k).z));

						Bz_out(_i, _j, _k).real(C_B * Bz_out(_i, _j, _k).real() + S_B * (k_norm(_i, _j, _k).x * Ey_out(_i, _j, _k).imag() - Ex_out(_i, _j, _k).imag() * k_norm(_i, _j, _k).y) \
							- 4 * PI * (1 - C_B) / (k_vec_len(_i, _j, _k) * c) * (k_norm(_i, _j, _k).x * Jy_out(_i, _j, _k).imag() - Jx_out(_i, _j, _k).imag() * k_norm(_i, _j, _k).y));

						Bz_out(_i, _j, _k).imag(C_B * Bz_out(_i, _j, _k).imag() - S_B * (k_norm(_i, _j, _k).x * Ey_out(_i, _j, _k).real() - Ex_out(_i, _j, _k).real() * k_norm(_i, _j, _k).y) \
							+ 4 * PI * (1 - C_B) / (k_vec_len(_i, _j, _k) * c) * (k_norm(_i, _j, _k).x * Jy_out(_i, _j, _k).real() - Jx_out(_i, _j, _k).real() * k_norm(_i, _j, _k).y));

						Ex_out(_i, _j, _k).real(C_E * Ex_out(_i, _j, _k).real() - S_E * (k_norm(_i, _j, _k).y * Bz_out(_i, _j, _k).imag() - By_out(_i, _j, _k).imag() * k_norm(_i, _j, _k).z) \
							- 4 * PI * S_E / (k_vec_len(_i, _j, _k) * c) * Jx_out(_i, _j, _k).real() + (1 - C_E) * k_norm(_i, _j, _k).x * (k_norm(_i, _j, _k).x * Ex_out(_i, _j, _k).real() \
								+ k_norm(_i, _j, _k).y * Ey_out(_i, _j, _k).real() + k_norm(_i, _j, _k).z * Ez_out(_i, _j, _k).real()) + 4 * PI * k_norm(_i, _j, _k).x * \
								(S_E / (k_vec_len(_i, _j, _k) * c) - steps.dt) * (k_norm(_i, _j, _k).x * Jx_out(_i, _j, _k).real() + k_norm(_i, _j, _k).y * Jy_out(_i, _j, _k).real() + \
									k_norm(_i, _j, _k).z * Jz_out(_i, _j, _k).real()));

						Ex_out(_i, _j, _k).imag(C_E * Ex_out(_i, _j, _k).imag() + S_E * (k_norm(_i, _j, _k).y * Bz_out(_i, _j, _k).real() - By_out(_i, _j, _k).real() * k_norm(_i, _j, _k).z) \
							- 4 * PI * S_E / (k_vec_len(_i, _j, _k) * c) * Jx_out(_i, _j, _k).imag() + (1 - C_E) * k_norm(_i, _j, _k).x * (k_norm(_i, _j, _k).x * Ex_out(_i, _j, _k).imag() \
								+ k_norm(_i, _j, _k).y * Ey_out(_i, _j, _k).imag() + k_norm(_i, _j, _k).z * Ez_out(_i, _j, _k).imag()) + 4 * PI * k_norm(_i, _j, _k).x * \
								(S_E / (k_vec_len(_i, _j, _k) * c) - steps.dt) * (k_norm(_i, _j, _k).x * Jx_out(_i, _j, _k).imag() + k_norm(_i, _j, _k).y * Jy_out(_i, _j, _k).imag() + \
									k_norm(_i, _j, _k).z * Jz_out(_i, _j, _k).imag()));

						Ey_out(_i, _j, _k).real(C_E * Ey_out(_i, _j, _k).real() + S_E * (k_norm(_i, _j, _k).x * Bz_out(_i, _j, _k).imag() - Bx_out(_i, _j, _k).imag() * k_norm(_i, _j, _k).z) \
							- 4 * PI * S_E / (k_vec_len(_i, _j, _k) * c) * Jy_out(_i, _j, _k).real() + (1 - C_E) * k_norm(_i, _j, _k).y * (k_norm(_i, _j, _k).x * Ex_out(_i, _j, _k).real() \
								+ k_norm(_i, _j, _k).y * Ey_out(_i, _j, _k).real() + k_norm(_i, _j, _k).z * Ez_out(_i, _j, _k).real()) + 4 * PI * k_norm(_i, _j, _k).y * \
								(S_E / (k_vec_len(_i, _j, _k) * c) - steps.dt) * (k_norm(_i, _j, _k).x * Jx_out(_i, _j, _k).real() + k_norm(_i, _j, _k).y * Jy_out(_i, _j, _k).real() + \
									k_norm(_i, _j, _k).z * Jz_out(_i, _j, _k).real()));

						Ey_out(_i, _j, _k).imag(C_E * Ey_out(_i, _j, _k).imag() - S_E * (k_norm(_i, _j, _k).x * Bz_out(_i, _j, _k).real() - Bx_out(_i, _j, _k).real() * k_norm(_i, _j, _k).z) \
							- 4 * PI * S_E / (k_vec_len(_i, _j, _k) * c) * Jy_out(_i, _j, _k).imag() + (1 - C_E) * k_norm(_i, _j, _k).y * (k_norm(_i, _j, _k).x * Ex_out(_i, _j, _k).imag() \
								+ k_norm(_i, _j, _k).y * Ey_out(_i, _j, _k).imag() + k_norm(_i, _j, _k).z * Ez_out(_i, _j, _k).imag()) + 4 * PI * k_norm(_i, _j, _k).y * \
								(S_E / (k_vec_len(_i, _j, _k) * c) - steps.dt) * (k_norm(_i, _j, _k).x * Jx_out(_i, _j, _k).imag() + k_norm(_i, _j, _k).y * Jy_out(_i, _j, _k).imag() + \
									k_norm(_i, _j, _k).z * Jz_out(_i, _j, _k).imag()));

						Ez_out(_i, _j, _k).real(C_E * Ez_out(_i, _j, _k).real() - S_E * (k_norm(_i, _j, _k).x * By_out(_i, _j, _k).imag() - Bx_out(_i, _j, _k).imag() * k_norm(_i, _j, _k).y) \
							- 4 * PI * S_E / (k_vec_len(_i, _j, _k) * c) * Jz_out(_i, _j, _k).real() + (1 - C_E) * k_norm(_i, _j, _k).z * (k_norm(_i, _j, _k).x * Ex_out(_i, _j, _k).real() \
								+ k_norm(_i, _j, _k).y * Ey_out(_i, _j, _k).real() + k_norm(_i, _j, _k).z * Ez_out(_i, _j, _k).real()) + 4 * PI * k_norm(_i, _j, _k).z * \
								(S_E / (k_vec_len(_i, _j, _k) * c) - steps.dt) * (k_norm(_i, _j, _k).x * Jx_out(_i, _j, _k).real() + k_norm(_i, _j, _k).y * Jy_out(_i, _j, _k).real() + \
									k_norm(_i, _j, _k).z * Jz_out(_i, _j, _k).real()));

						Ez_out(_i, _j, _k).imag(C_E * Ez_out(_i, _j, _k).imag() + S_E * (k_norm(_i, _j, _k).x * By_out(_i, _j, _k).real() - Bx_out(_i, _j, _k).real() * k_norm(_i, _j, _k).y) \
							- 4 * PI * S_E / (k_vec_len(_i, _j, _k) * c) * Jz_out(_i, _j, _k).imag() + (1 - C_E) * k_norm(_i, _j, _k).z * (k_norm(_i, _j, _k).x * Ex_out(_i, _j, _k).imag() \
								+ k_norm(_i, _j, _k).y * Ey_out(_i, _j, _k).imag() + k_norm(_i, _j, _k).z * Ez_out(_i, _j, _k).imag()) + 4 * PI * k_norm(_i, _j, _k).z * \
								(S_E / (k_vec_len(_i, _j, _k) * c) - steps.dt) * (k_norm(_i, _j, _k).x * Jx_out(_i, _j, _k).imag() + k_norm(_i, _j, _k).y * Jy_out(_i, _j, _k).imag() + \
									k_norm(_i, _j, _k).z * Jz_out(_i, _j, _k).imag()));*/
						Bx_out(_i, _j, _k) = C_B * Bx_out(_i, _j, _k) - I * S_B * (k_norm(_i, _j, _k).y * Ez_out(_i, _j, _k) - Ey_out(_i, _j, _k) * k_norm(_i, _j, _k).z) \
							+ 4 * PI * I * (1 - C_B) / (k_vec_len(_i, _j, _k) * c) * (k_norm(_i, _j, _k).y * Jz_out(_i, _j, _k) - Jy_out(_i, _j, _k) * k_norm(_i, _j, _k).z);

						By_out(_i, _j, _k) = C_B * By_out(_i, _j, _k) + I * S_B * (k_norm(_i, _j, _k).x * Ez_out(_i, _j, _k) - Ex_out(_i, _j, _k) * k_norm(_i, _j, _k).z) \
							- 4 * PI * I * (1 - C_B) / (k_vec_len(_i, _j, _k) * c) * (k_norm(_i, _j, _k).x * Jz_out(_i, _j, _k) - Jx_out(_i, _j, _k) * k_norm(_i, _j, _k).z);

						Bz_out(_i, _j, _k) = (C_B * Bz_out(_i, _j, _k) - I * S_B * (k_norm(_i, _j, _k).x * Ey_out(_i, _j, _k) - Ex_out(_i, _j, _k) * k_norm(_i, _j, _k).y) \
							+ 4 * PI * I * (1 - C_B) / (k_vec_len(_i, _j, _k) * c) * (k_norm(_i, _j, _k).x * Jy_out(_i, _j, _k) - Jx_out(_i, _j, _k) * k_norm(_i, _j, _k).y));

						Ex_out(_i, _j, _k) = C_E * Ex_out(_i, _j, _k) + I * S_E * (k_norm(_i, _j, _k).y * Bz_out(_i, _j, _k) - By_out(_i, _j, _k) * k_norm(_i, _j, _k).z) \
							- 4 * PI * S_E / (k_vec_len(_i, _j, _k) * c) * Jx_out(_i, _j, _k) + (1 - C_E) * k_norm(_i, _j, _k).x * (k_norm(_i, _j, _k).x * Ex_out(_i, _j, _k) \
								+ k_norm(_i, _j, _k).y * Ey_out(_i, _j, _k) + k_norm(_i, _j, _k).z * Ez_out(_i, _j, _k)) + 4 * PI * k_norm(_i, _j, _k).x * \
								(S_E / (k_vec_len(_i, _j, _k) * c) - steps.dt) * (k_norm(_i, _j, _k).x * Jx_out(_i, _j, _k) + k_norm(_i, _j, _k).y * Jy_out(_i, _j, _k) + \
									k_norm(_i, _j, _k).z * Jz_out(_i, _j, _k));

						Ey_out(_i, _j, _k) = C_E * Ey_out(_i, _j, _k) - I * S_E * (k_norm(_i, _j, _k).x * Bz_out(_i, _j, _k) - Bx_out(_i, _j, _k) * k_norm(_i, _j, _k).z) \
							- 4 * PI * S_E / (k_vec_len(_i, _j, _k) * c) * Jy_out(_i, _j, _k) + (1 - C_E) * k_norm(_i, _j, _k).y * (k_norm(_i, _j, _k).x * Ex_out(_i, _j, _k) \
								+ k_norm(_i, _j, _k).y * Ey_out(_i, _j, _k) + k_norm(_i, _j, _k).z * Ez_out(_i, _j, _k)) + 4 * PI * k_norm(_i, _j, _k).y * \
								(S_E / (k_vec_len(_i, _j, _k) * c) - steps.dt) * (k_norm(_i, _j, _k).x * Jx_out(_i, _j, _k) + k_norm(_i, _j, _k).y * Jy_out(_i, _j, _k) + \
									k_norm(_i, _j, _k).z * Jz_out(_i, _j, _k));

						Ez_out(_i, _j, _k) = C_E * Ez_out(_i, _j, _k) + I * S_E * (k_norm(_i, _j, _k).x * By_out(_i, _j, _k) - Bx_out(_i, _j, _k) * k_norm(_i, _j, _k).y) \
							- 4 * PI * S_E / (k_vec_len(_i, _j, _k) * c) * Jz_out(_i, _j, _k) + (1 - C_E) * k_norm(_i, _j, _k).z * (k_norm(_i, _j, _k).x * Ex_out(_i, _j, _k) \
								+ k_norm(_i, _j, _k).y * Ey_out(_i, _j, _k) + k_norm(_i, _j, _k).z * Ez_out(_i, _j, _k)) + 4 * PI * k_norm(_i, _j, _k).z * \
								(S_E / (k_vec_len(_i, _j, _k) * c) - steps.dt) * (k_norm(_i, _j, _k).x * Jx_out(_i, _j, _k) + k_norm(_i, _j, _k).y * Jy_out(_i, _j, _k) + \
									k_norm(_i, _j, _k).z * Jz_out(_i, _j, _k));

					}
				}
			}
		}

		_Ex_plan = fftw_plan_dft_c2r_3d(nx, ny, nz, (fftw_complex*)(Ex_out.data()), E.x.data(), FFTW_ESTIMATE);
		_Ey_plan = fftw_plan_dft_c2r_3d(nx, ny, nz, (fftw_complex*)(Ey_out.data()), E.y.data(), FFTW_ESTIMATE);
		_Ez_plan = fftw_plan_dft_c2r_3d(nx, ny, nz, (fftw_complex*)(Ez_out.data()), E.z.data(), FFTW_ESTIMATE);
		_Bx_plan = fftw_plan_dft_c2r_3d(nx, ny, nz, (fftw_complex*)(Bx_out.data()), B.x.data(), FFTW_ESTIMATE);
		_By_plan = fftw_plan_dft_c2r_3d(nx, ny, nz, (fftw_complex*)(By_out.data()), B.y.data(), FFTW_ESTIMATE);
		_Bz_plan = fftw_plan_dft_c2r_3d(nx, ny, nz, (fftw_complex*)(Bz_out.data()), B.z.data(), FFTW_ESTIMATE);
		_Jx_plan = fftw_plan_dft_c2r_3d(nx, ny, nz, (fftw_complex*)(Jx_out.data()), J.x.data(), FFTW_ESTIMATE);
		_Jy_plan = fftw_plan_dft_c2r_3d(nx, ny, nz, (fftw_complex*)(Jy_out.data()), J.y.data(), FFTW_ESTIMATE);
		_Jz_plan = fftw_plan_dft_c2r_3d(nx, ny, nz, (fftw_complex*)(Jz_out.data()), J.z.data(), FFTW_ESTIMATE);

		fftw_execute(_Ex_plan);
		fftw_execute(_Ey_plan);
		fftw_execute(_Ez_plan);
		fftw_execute(_Bx_plan);
		fftw_execute(_By_plan);
		fftw_execute(_Bz_plan);
		fftw_execute(_Jx_plan);
		fftw_execute(_Jy_plan);
		fftw_execute(_Jz_plan);

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

		fftw_destroy_plan(Ex_plan);
		fftw_destroy_plan(Ey_plan);
		fftw_destroy_plan(Ez_plan);
		fftw_destroy_plan(Bx_plan);
		fftw_destroy_plan(By_plan);
		fftw_destroy_plan(Bz_plan);
		fftw_destroy_plan(Jx_plan);
		fftw_destroy_plan(Jy_plan);
		fftw_destroy_plan(Jz_plan);

		fftw_destroy_plan(_Ex_plan);
		fftw_destroy_plan(_Ey_plan);
		fftw_destroy_plan(_Ez_plan);
		fftw_destroy_plan(_Bx_plan);
		fftw_destroy_plan(_By_plan);
		fftw_destroy_plan(_Bz_plan);
		fftw_destroy_plan(_Jx_plan);
		fftw_destroy_plan(_Jy_plan);
		fftw_destroy_plan(_Jz_plan);

		return (*this);
	}
};