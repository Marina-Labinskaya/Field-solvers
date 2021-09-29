#pragma once

#include <cmath>
#include <iostream>
#include <vector>
#include <fstream>

const double c = 3 * pow(10, 8); // скорость света
const double dt = 1 * pow(10, -5); // шаг по времени
const double dx = 5000.0; // шаг по координате x
const double dy = 5000.0; // шаг по координате y

const double PI = 3.14;

// границы расчетной области

// по координате x
const double A = 10000.0;
const double B = 60000.0;

// по координате y
const double C = 20000.0;
const double D = 80000.0;

// количество узлов по x и y
const int m = int((B - A) / dx + 1);
const int n = int((D - C) / dy + 1);

// по t
const double T = 1 * pow (10, -3);

struct elec_field {
	std::vector<std::vector<double>> Ex;
	std::vector<std::vector<double>> Ey;
	std::vector<std::vector<double>> Ez;
};

struct magn_field {
	std::vector<std::vector<double>> Bx;
	std::vector<std::vector<double>> By;
	std::vector<std::vector<double>> Bz;
};

struct elec_curr {
	std::vector<std::vector<double>> Jx;
	std::vector<std::vector<double>> Jy;
	std::vector<std::vector<double>> Jz;
};

class elec_magn_field {

public:
	elec_field E;
	magn_field B;
	elec_curr J;

	elec_magn_field() {

		E.Ex.assign(m, std::vector<double>(n));
		E.Ey.assign(m, std::vector<double>(n));
		E.Ez.assign(m, std::vector<double>(n));

		B.Bx.assign(m, std::vector<double>(n));
		B.By.assign(m, std::vector<double>(n));
		B.Bz.assign(m, std::vector<double>(n));

		J.Jx.assign(m, std::vector<double>(n));
		J.Jy.assign(m, std::vector<double>(n));
		J.Jz.assign(m, std::vector<double>(n));

		for (unsigned int i = 0; i < m; ++i)
			for (unsigned int j = 0; j < n; ++j) {
				E.Ex[i][j] = 0.0;
				E.Ey[i][j] = 0.0;
				E.Ez[i][j] = 0.0;

				B.Bx[i][j] = 0.0;
				B.By[i][j] = 0.0;
				B.Bz[i][j] = 0.0;

				J.Jx[i][j] = 0.0;
				J.Jy[i][j] = 0.0;
				J.Jz[i][j] = 0.0;
			}

	}

	elec_magn_field(double(*pf)(double, double, double)) {

		E.Ex.assign(m, std::vector<double>(n));
		E.Ey.assign(m, std::vector<double>(n));
		E.Ez.assign(m, std::vector<double>(n));

		B.Bx.assign(m, std::vector<double>(n));
		B.By.assign(m, std::vector<double>(n));
		B.Bz.assign(m, std::vector<double>(n));

		J.Jx.assign(m, std::vector<double>(n));
		J.Jy.assign(m, std::vector<double>(n));
		J.Jz.assign(m, std::vector<double>(n));

		for (unsigned int i = 0; i < m; ++i)
			for (unsigned int j = 0; j < n; ++j) {

				E.Ex[i][j] = 0.0;
				E.Ey[i][j] = pf(dx * i + A, 0.0, 0.0);
				E.Ez[i][j] = 0.0;

				B.Bx[i][j] = 0.0;
				B.By[i][j] = 0.0;
				B.Bz[i][j] = E.Ey[i][j];

				J.Jx[i][j] = 0.0;
				J.Jy[i][j] = 0.0;
				J.Jz[i][j] = 0.0;

				//std::cout << "Ex[" << i << "][" << j << "]=" << E.Ex[i][j] << std::endl;
				//std::cout << "Ey[" << i << "][" << j << "]=" << E.Ey[i][j] << std::endl;
				//std::cout << "Bz[" << i << "][" << j << "]=" << B.Bz[i][j] << std::endl;
				//std::cout << std::endl;
			}
	}

	   elec_magn_field& field_integrate(double t) {
		// интегрирование

		//	for (double t = dt; t < T; t += dt) {

		// обновление Ex

		for (unsigned int i = 0; i < m; ++i) {
			E.Ex[i][0] = E.Ex[i][0] - 4.0 * PI * dt * J.Jx[i][0] + c * dt * ((B.Bz[i][1] - B.Bz[i][n - 1]) / (2.0 * dy));
			E.Ex[i][n - 1] = E.Ex[i][n - 1] - 4.0 * PI * dt * J.Jx[i][n - 1] + c * dt * ((B.Bz[i][0] - B.Bz[i][n - 2]) / (2.0 * dy));
			//std::cout << "i= " << i << std::endl;
			//std::cout << std::endl;
			//std::cout << "Ex[i][0]= " << E.Ex[i][0] << std::endl;
			//std::cout << "Ex[i][n-1]= " << E.Ex[i][n - 1] << std::endl;
			for (unsigned int j = 1; j < n - 1; ++j) {
				//std::cout << "i= " << i << ", j= " << j << std::endl;
				//std::cout << std::endl;
				E.Ex[i][j] = E.Ex[i][j] - 4.0 * PI * dt * J.Jx[i][j] + c * dt * ((B.Bz[i][j + 1] - B.Bz[i][j - 1]) / (2.0 * dy));
				//std::cout << "Ex[i][j]= " << E.Ex[i][j] << std::endl;
			}
		}

		// обновление Ey

		for (unsigned int j = 0; j < n; ++j) {
			E.Ey[0][j] = E.Ey[0][j] - 4.0 * PI * dt * J.Jy[0][j] - c * dt * ((B.Bz[1][j] - B.Bz[m - 1][j]) / (2.0 * dx));
			E.Ey[m - 1][j] = E.Ey[m - 1][j] - 4.0 * PI * dt * J.Jy[m - 1][j] - c * dt * ((B.Bz[0][j] - B.Bz[m - 2][j]) / (2.0 * dx));
			//std::cout << "j= " << j << std::endl;
			//std::cout << std::endl;
			//std::cout << "Ey[0][j]= " << E.Ey[0][j] << std::endl;
			//std::cout << "Ey[m-1][j]= " << E.Ey[m - 1][j] << std::endl;
			for (unsigned int i = 1; i < m - 1; ++i) {
				//std::cout << "i= " << i << ", j= " << j << std::endl;
				//std::cout << std::endl;
				E.Ey[i][j] = E.Ey[i][j] - 4.0 * PI * dt * J.Jy[i][j] - c * dt * ((B.Bz[i + 1][j] - B.Bz[i - 1][j]) / (2.0 * dx));
				//std::cout << "Ey[i][j]= " << E.Ey[i][j] << std::endl;
			}

		}

		// обновление Ez

		E.Ez[0][0] = E.Ez[0][0] - 4.0 * PI * dt * J.Jz[0][0] + c * dt * \
			((B.By[1][0] - B.By[m - 1][0]) / (2.0 * dx) - (B.Bx[0][1] - B.Bx[0][n - 1]) / (2.0 * dy));
		E.Ez[0][n - 1] = E.Ez[0][n - 1] - 4.0 * PI * dt * J.Jz[0][n - 1] + c * dt * \
			((B.By[1][n - 1] - B.By[m - 1][n - 1]) / (2.0 * dx) - (B.Bx[0][0] - B.Bx[0][n - 2]) / (2.0 * dy));
		E.Ez[m - 1][0] = E.Ez[m - 1][0] - 4.0 * PI * dt * J.Jz[m - 1][0] + c * dt * \
			((B.By[0][0] - B.By[m - 2][0]) / (2.0 * dx) - (B.Bx[m - 1][1] - B.Bx[m - 1][n - 1]) / (2.0 * dy));
		E.Ez[m - 1][n - 1] = E.Ez[m - 1][n - 1] - 4.0 * PI * dt * J.Jz[m - 1][n - 1] + c * dt * \
			((B.By[0][n - 1] - B.By[m - 2][n - 1]) / (2.0 * dx) - (B.Bx[m - 1][0] - B.Bx[m - 1][n - 2]) / (2.0 * dy));

		for (unsigned int i = 1; i < m - 1; ++i) {
			E.Ez[i][0] = E.Ez[i][0] - 4.0 * PI * dt * J.Jz[i][0] + c * dt * \
				((B.By[i + 1][0] - B.By[i - 1][0]) / (2.0 * dx) - (B.Bx[i][1] - B.Bx[i][n - 1]) / (2.0 * dy));
			E.Ez[i][n - 1] = E.Ez[i][n - 1] - 4.0 * PI * dt * J.Jz[i][n - 1] + c * dt * \
				((B.By[i + 1][n - 1] - B.By[i - 1][n - 1]) / (2.0 * dx) - (B.Bx[i][0] - B.Bx[i][n - 2]) / (2.0 * dy));
			//std::cout << "Ez[i][0]= " << E.Ez[i][0] << std::endl;
			//std::cout << "Ez[i][n - 1]= " << E.Ez[i][n - 1] << std::endl;
			//std::cout << std::endl;
		}

		for (unsigned int j = 1; j < n - 1; ++j) {
			E.Ez[0][j] = E.Ez[0][j] - 4.0 * PI * dt * J.Jz[0][j] + c * dt * \
				((B.By[1][j] - B.By[m - 1][j]) / (2.0 * dx) - (B.Bx[0][j + 1] - B.Bx[0][j - 1]) / (2.0 * dy));
			E.Ez[m - 1][j] = E.Ez[m - 1][j] - 4.0 * PI * dt * J.Jz[m - 1][j] + c * dt * \
				((B.By[0][j] - B.By[m - 2][j]) / (2.0 * dx) - (B.Bx[m - 1][j + 1] - B.Bx[m - 1][j - 1]) / (2.0 * dy));
			//std::cout << "Ez[0][j]= " << E.Ez[0][j] << std::endl;
			//std::cout << "Ez[m-1][j]= " << E.Ez[m - 1][j] << std::endl;
			//std::cout << std::endl;

		}

		for (unsigned int i = 1; i < m - 1; ++i)
			for (unsigned int j = 1; j < n - 1; ++j) {
				E.Ez[i][j] = E.Ez[i][j] - 4.0 * PI * dt * J.Jz[i][j] + c * dt * \
					((B.By[i + 1][j] - B.By[i - 1][j]) / (2.0 * dx) - (B.Bx[i][j + 1] - B.Bx[i][j - 1]) / (2.0 * dy));
				//std::cout << "i= " << i << ", j= " << j << std::endl;
				//std::cout << std::endl;
				//std::cout << "Ez[i][j]= " << E.Ez[i][j] << std::endl;
				//std::cout << std::endl;
			}

		// обновление Bx

		for (unsigned int i = 0; i < m; ++i) {
			B.Bx[i][0] =B.Bx[i][0] - c * dt * ((E.Ez[i][1] - E.Ez[i][n - 1]) / (2.0 * dy));
			B.Bx[i][n - 1] = B.Bx[i][n - 1] - c * dt * ((E.Ez[i][0] - E.Ez[i][n - 2]) / (2.0 * dy));
			//std::cout << "i= " << i << std::endl;
			//std::cout << std::endl;
			//std::cout << "Bx[i][0]= " << B.Bx[i][0] << std::endl;
			//std::cout << "Bx[i][n-1]= " << B.Bx[i][n - 1] << std::endl;
			for (unsigned int j = 1; j < n - 1; ++j) {
				B.Bx[i][j] = B.Bx[i][j] - c * dt * ((E.Ez[i][j + 1] - E.Ez[i][j - 1]) / (2.0 * dy));
				//std::cout << "i= " << i << ", j= " << j << std::endl;
				//std::cout << std::endl;
				//std::cout << "Bx[i][j]= " << B.Bx[i][j] << std::endl;
			}
		}

		// обновление By

		for (unsigned int j = 0; j < n; ++j) {
			B.By[0][j] = B.By[0][j] + c * dt * ((E.Ez[1][j] - E.Ez[m - 1][j]) / (2.0 * dx));
			B.By[m - 1][j] = B.By[m - 1][j] + c * dt * ((E.Ez[0][j] - E.Ez[m - 2][j]) / (2.0 * dx));
			//std::cout << "j= " << j << std::endl;
			//std::cout << std::endl;
			//std::cout << "By[0][j]= " << B.By[0][j] << std::endl;
			//std::cout << "By[m-1][j]= " << B.By[m - 1][j] << std::endl;
			for (unsigned int i = 1; i < m - 1; ++i) {
				//std::cout << "i= " << i << ", j= " << j << std::endl;
				//std::cout << std::endl;
				B.By[i][j] = B.By[i][j] + c * dt * ((E.Ez[i + 1][j] - E.Ez[i - 1][j]) / (2.0 * dx));
				//std::cout << "By[i][j]= " << B.By[i][j] << std::endl;
			}
		}

		// обновление Bz

		B.Bz[0][0] = B.Bz[0][0] - c * dt * \
			((E.Ey[1][0] - E.Ey[m - 1][0]) / (2.0 * dx) - (E.Ex[0][1] - E.Ex[0][n - 1]) / (2.0 * dy));
		//std::cout << "Bz[0][0]= " << B.Bz[0][0] << std::endl;
		//std::cout << std::endl;
		B.Bz[0][n - 1] = B.Bz[0][n - 1] - c * dt * \
			((E.Ey[1][n - 1] - E.Ey[m - 1][n - 1]) / (2.0 * dx) - (E.Ex[0][0] - E.Ex[0][n - 2]) / (2.0 * dy));
		//std::cout << "Bz[0][n - 1]= " << B.Bz[0][n - 1] << std::endl;
		//std::cout << std::endl;
		B.Bz[m - 1][0] = B.Bz[m - 1][0] - c * dt * \
			((E.Ey[0][0] - E.Ey[m - 2][0]) / (2.0 * dx) - (E.Ex[m - 1][1] - E.Ex[m - 1][n - 1]) / (2.0 * dy));
		//std::cout << "Bz[m - 1][0]= " << B.Bz[m - 1][0] << std::endl;
		//std::cout << std::endl;
		B.Bz[m - 1][n - 1] = B.Bz[m - 1][n - 1] - c * dt * \
			((E.Ey[0][n - 1] - E.Ey[m - 2][n - 1]) / (2.0 * dx) - (E.Ex[m - 1][0] - E.Ex[m - 1][n - 2]) / (2.0 * dy));
		//std::cout << "Bz[m - 1][n - 1]= " << B.Bz[m - 1][n - 1] << std::endl;
		//std::cout << std::endl;

		for (unsigned int i = 1; i < m - 1; ++i) {
			B.Bz[i][0] = B.Bz[i][0] - c * dt * \
				((E.Ey[i + 1][0] - E.Ey[i - 1][0]) / (2.0 * dx) - (E.Ex[i][1] - E.Ex[i][n - 2]) / (2.0 * dy));
			B.Bz[i][n - 1] = B.Bz[i][n - 1] - c * dt * \
				((E.Ey[i + 1][n - 1] - E.Ey[i - 1][n - 1]) / (2.0 * dx) - (E.Ex[i][0] - E.Ex[i][n - 2]) / (2.0 * dy));
			//std::cout << "i= " << i << std::endl;
			//std::cout << "Bz[i][0]= " << B.Bz[i][0] << std::endl;
			//std::cout << "Bz[i][n-1]= " << B.Bz[i][n - 1] << std::endl;
			//std::cout << std::endl;
		}

		for (unsigned int j = 1; j < n - 1; ++j) {
			B.Bz[0][j] = B.Bz[0][j] - c * dt * \
				((E.Ey[1][j] - E.Ey[m - 1][j]) / (2.0 * dx) - (E.Ex[0][j + 1] - E.Ex[0][j - 1]) / (2.0 * dy));
			B.Bz[m - 1][j] = B.Bz[m - 1][j] - c * dt * \
				((E.Ey[0][j] - E.Ey[m - 2][j]) / (2.0 * dx) - (E.Ex[m - 1][j + 1] - E.Ex[m - 1][j - 1]) / (2.0 * dy));
			//std::cout << "j= " << j << std::endl;
			//std::cout << "Bz[i][j]= " << B.Bz[0][j] << std::endl;
			//std::cout << "Bz[i][j]= " << B.Bz[m - 1][j] << std::endl;
			//std::cout << std::endl;
		}

		for (unsigned int i = 1; i < m - 1; ++i)
			for (unsigned int j = 1; j < n - 1; ++j) {
				B.Bz[i][j] = B.Bz[i][j] - c * dt * \
					((E.Ey[i + 1][j] - E.Ey[i - 1][j]) / (2.0 * dx) - (E.Ex[i][j + 1] - E.Ex[i][j - 1]) / (2.0 * dy));
				//std::cout << " i= " << i << " j= " << j << std::endl;
				//std::cout << std::endl;
				//std::cout << "Bz[i][j]= " << B.Bz[i][j] << std::endl;
				//std::cout << std::endl;

			}
		//	}
		return (*this);
	}
};

