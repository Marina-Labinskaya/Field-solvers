#pragma once

#include <cmath>
#include <iostream>
#include <vector>
#include <fstream>

const double c = 3 * pow(10, 8); // скорость света
const double dx = 5 * pow(10, 8); // шаг по координате x
const double dy = 5 * pow(10, 8); // шаг по координате y
const double dt = 1 * pow(10, -3);

const double PI = 3.14;

// границы расчетной области

// по координате x
const double A = 0.0;
const double B = 10 * pow(10, 9);

// по координате y
const double C = 0.0;
const double D = 5 * pow(10, 9);

// количество узлов по x и y
const int m = int((B - A) / dx + 1);
const int n = int((D - C) / dy + 1);

// по t
const double T = 1 * pow (10, -1);

struct field_characteristics {
	std::vector<std::vector<double>> x;
	std::vector<std::vector<double>> y;
	std::vector<std::vector<double>> z;
};

class elec_magn_field {

public:
	field_characteristics E, B, J;

	elec_magn_field() {

		E.x.assign(m, std::vector<double>(n));
		E.y.assign(m, std::vector<double>(n));
		E.z.assign(m, std::vector<double>(n));

		B.x.assign(m, std::vector<double>(n));
		B.y.assign(m, std::vector<double>(n));
		B.z.assign(m, std::vector<double>(n));

		J.x.assign(m, std::vector<double>(n));
		J.y.assign(m, std::vector<double>(n));
		J.z.assign(m, std::vector<double>(n));

		for (unsigned int i = 0; i < m; ++i)
			for (unsigned int j = 0; j < n; ++j) {
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

	elec_magn_field(double(*pf)(double, double, double)) {

		E.x.assign(m, std::vector<double>(n));
		E.y.assign(m, std::vector<double>(n));
		E.z.assign(m, std::vector<double>(n));

		B.x.assign(m, std::vector<double>(n));
		B.y.assign(m, std::vector<double>(n));
		B.z.assign(m, std::vector<double>(n));

		J.x.assign(m, std::vector<double>(n));
		J.y.assign(m, std::vector<double>(n));
		J.z.assign(m, std::vector<double>(n));

		for (unsigned int i = 0; i < m; ++i)
			for (unsigned int j = 0; j < n; ++j) {

				E.x[i][j] = 0.0;
				E.y[i][j] = pf(dx * i + A, 0.0, 0.0);
				E.z[i][j] = 0.0;

				B.x[i][j] = 0.0;
				B.y[i][j] = 0.0;
				B.z[i][j] = E.y[i][j];

				J.x[i][j] = 0.0;
				J.y[i][j] = 0.0;
				J.z[i][j] = 0.0;

				//std::cout << "Ex[" << i << "][" << j << "]=" << E.x[i][j] << std::endl;
				//std::cout << "Ey[" << i << "][" << j << "]=" << E.y[i][j] << std::endl;
				//std::cout << "Bz[" << i << "][" << j << "]=" << B.z[i][j] << std::endl;
				//std::cout << std::endl;
			}
	}

	   elec_magn_field& field_integrate() {
		// интегрирование

		/*// обновление E

		   for (int i = 1; i < m - 1; i++)
			   for (int j = 1; j < n - 1; j++) {
				   E.x[i][j] = E.x[i][j] - 4.0 * PI * dt * J.x[i][j] + c * dt * (B.z[i][j + 1] - B.z[i][j - 1]) / (2.0 * dy);
				   E.y[i][j] = E.y[i][j] - 4.0 * PI * dt * J.y[i][j] - c * dt * (B.z[i + 1][j] - B.z[i - 1][j]) / (2.0 * dx);
				   E.z[i][j] = E.z[i][j] - 4.0 * PI * dt * J.z[i][j] + c * dt * ((B.y[i + 1][j] - B.y[i - 1][j]) / (2.0 * dx) - \
					   (B.x[i][j + 1] - B.x[i][j - 1]) / (2.0 * dy));
			   }

		   for (int i = 0; i < m; i++) { // j==0
			   E.x[i][0] = E.x[i][0] - 4.0 * PI * dt * J.x[i][0] + c * dt * (B.z[i][1] - B.z[i][n - 1]) / (2.0 * dy);
			   if ((i != 0) && (i != m - 1)) {
				   E.y[i][0] = E.y[i][0] - 4.0 * PI * dt * J.y[i][0] - c * dt * (B.z[i + 1][0] - B.z[i - 1][0]) / (2.0 * dx);
				   E.z[i][0] = E.z[i][0] - 4.0 * PI * dt * J.z[i][0] + c * dt * ((B.y[i + 1][0] - B.y[i - 1][0]) / (2.0 * dx) - \
					   (B.x[i][1] - B.x[i][n - 1]) / (2.0 * dy));
			   }
			   if (i == 0) {
				   E.y[i][0] = E.y[i][0] - 4.0 * PI * dt * J.y[i][0] - c * dt * (B.z[i + 1][0] - B.z[m - 1][0]) / (2.0 * dx);
				   E.z[i][0] = E.z[i][0] - 4.0 * PI * dt * J.z[i][0] + c * dt * ((B.y[i + 1][0] - B.y[m - 1][0]) / (2.0 * dx) - \
					   (B.x[i][1] - B.x[i][n - 1]) / (2.0 * dy));
			   }
			   if (i == m - 1) {
				   E.y[i][0] = E.y[i][0] - 4.0 * PI * dt * J.y[i][0] - c * dt * (B.z[0][0] - B.z[i - 1][0]) / (2.0 * dx);
				   E.z[i][0] = E.z[i][0] - 4.0 * PI * dt * J.z[i][0] + c * dt * ((B.y[0][0] - B.y[i - 1][0]) / (2.0 * dx) - \
					   (B.x[i][1] - B.x[i][n - 1]) / (2.0 * dy));
			   }
		   }

		   for (int i = 0; i < m; i++) { // j==ny-1
			   E.x[i][n - 1] = E.x[i][n - 1] - 4.0 * PI * dt * J.x[i][n - 1] + c * dt * (B.z[i][0] - B.z[i][n - 2]) / (2.0 * dy);
			   if ((i != 0) && (i != m - 1)) {
				   E.y[i][n - 1] = E.y[i][n - 1] - 4.0 * PI * dt * J.y[i][n - 1] - c * dt * (B.z[i + 1][n - 1] - B.z[i - 1][n - 1]) / (2.0 * dx);
				   E.z[i][n - 1] = E.z[i][n - 1] - 4.0 * PI * dt * J.z[i][n - 1] + c * dt * ((B.y[i + 1][n - 1] - B.y[i - 1][n - 1]) / (2.0 * dx) - \
					   (B.x[i][0] - B.x[i][n - 2]) / (2.0 * dy));
			   }
			   if (i == 0) {
				   E.y[i][n - 1] = E.y[i][n - 1] - 4.0 * PI * dt * J.y[i][n - 1] - c * dt * (B.z[i + 1][n - 1] - B.z[m - 1][n - 1]) / (2.0 * dx);
				   E.z[i][n - 1] = E.z[i][n - 1] - 4.0 * PI * dt * J.z[i][n - 1] + c * dt * ((B.y[i + 1][n - 1] - B.y[m - 1][n - 1]) / (2.0 * dx) - \
					   (B.x[i][0] - B.x[i][n - 2]) / (2.0 * dy));
			   }
			   if (i == m - 1) {
				   E.y[i][n - 1] = E.y[i][n - 1] - 4.0 * PI * dt * J.y[i][n - 1] - c * dt * (B.z[0][n - 1] - B.z[i - 1][n - 1]) / (2.0 * dx);
				   E.z[i][n - 1] = E.z[i][n - 1] - 4.0 * PI * dt * J.z[i][n - 1] + c * dt * ((B.y[0][n - 1] - B.y[i - 1][n - 1]) / (2.0 * dx) - \
					   (B.x[i][0] - B.x[i][n - 2]) / (2.0 * dy));
			   }
		   }

		   for (int j = 1; j < n - 1; j++) { // i==0
			   E.x[0][j] = E.x[0][j] - 4.0 * PI * dt * J.x[0][j] + c * dt * (B.z[0][j + 1] - B.z[0][j - 1]) / (2.0 * dy);
			   E.y[0][j] = E.y[0][j] - 4.0 * PI * dt * J.y[0][j] - c * dt * (B.z[1][j] - B.z[m - 1][j]) / (2.0 * dx);
			   E.z[0][j] = E.z[0][j] - 4.0 * PI * dt * J.z[0][j] + c * dt * ((B.y[1][j] - B.y[m - 1][j]) / (2.0 * dx) - \
				   (B.x[0][j + 1] - B.x[0][j - 1]) / (2.0 * dy));
		   }

		   for (int j = 1; j < n - 1; j++) { // i==nx-1
			   E.x[m - 1][j] = E.x[m - 1][j] - 4.0 * PI * dt * J.x[m - 1][j] + c * dt * (B.z[m - 1][j + 1] - B.z[m - 1][j - 1]) / (2.0 * dy);
			   E.y[m - 1][j] = E.y[m - 1][j] - 4.0 * PI * dt * J.y[m - 1][j] - c * dt * (B.z[0][j] - B.z[m - 2][j]) / (2.0 * dx);
			   E.z[m - 1][j] = E.z[m - 1][j] - 4.0 * PI * dt * J.z[m - 1][j] + c * dt * ((B.y[0][j] - B.y[m - 2][j]) / (2.0 * dx) - \
				   (B.x[m - 1][j + 1] - B.x[m - 1][j - 1]) / (2.0 * dy));
		   }



		   // обновление B

		   for (int i = 1; i < m - 1; i++)
			   for (int j = 1; j < n - 1; j++) {
				   B.x[i][j] = B.x[i][j] - c * dt * (E.z[i][j + 1] - E.z[i][j - 1]) / (2.0 * dy);
				   B.y[i][j] = B.y[i][j] + c * dt * (E.z[i + 1][j] - E.z[i - 1][j]) / (2.0 * dx);
				   B.z[i][j] = B.z[i][j] - c * dt * ((E.y[i + 1][j] - E.y[i - 1][j]) / (2.0 * dx) - (E.x[i][j + 1] - E.x[i][j - 1]) / (2.0 * dy));
			   }

		   for (int i = 0; i < m; i++) { // j==0
			   B.x[i][0] = B.x[i][0] - c * dt * (E.z[i][1] - E.z[i][n - 1]) / (2.0 * dy);
			   if ((i != 0) && (i != m - 1)) {
				   B.y[i][0] = B.y[i][0] + c * dt * (E.z[i + 1][0] - E.z[i - 1][0]) / (2.0 * dx);
				   B.z[i][0] = B.z[i][0] - c * dt * ((E.y[i + 1][0] - E.y[i - 1][0]) / (2.0 * dx) - (E.x[i][1] - E.x[i][n - 1]) / (2.0 * dy));
			   }
			   if (i == 0) {
				   B.y[i][0] = B.y[i][0] + c * dt * (E.z[i + 1][0] - E.z[m - 1][0]) / (2.0 * dx);
				   B.z[i][0] = B.z[i][0] - c * dt * ((E.y[i + 1][0] - E.y[m - 1][0]) / (2.0 * dx) - (E.x[i][1] - E.x[i][n - 1]) / (2.0 * dy));
			   }
			   if (i == m - 1) {
				   B.y[i][0] = B.y[i][0] + c * dt * (E.z[0][0] - E.z[i - 1][0]) / (2.0 * dx);
				   B.z[i][0] = B.z[i][0] - c * dt * ((E.y[0][0] - E.y[i - 1][0]) / (2.0 * dx) - (E.x[i][1] - E.x[i][n - 1]) / (2.0 * dy));
			   }
		   }

		   for (int i = 0; i < m; i++) { // j==ny-1
			   B.x[i][n - 1] = B.x[i][n - 1] - c * dt * (E.z[i][0] - E.z[i][n - 2]) / (2.0 * dy);
			   if ((i != 0) && (i != m - 1)) {
				   B.y[i][n - 1] = B.y[i][n - 1] + c * dt * (E.z[i + 1][n - 1] - E.z[i - 1][n - 1]) / (2.0 * dx);
				   B.z[i][n - 1] = B.z[i][n - 1] - c * dt * ((E.y[i + 1][n - 1] - E.y[i - 1][n - 1]) / (2.0 * dx) - (E.x[i][0] - E.x[i][n - 2]) / (2.0 * dy));
			   }
			   if (i == 0) {
				   B.y[i][n - 1] = B.y[i][n - 1] + c * dt * (E.z[i + 1][n - 1] - E.z[m - 1][n - 1]) / (2.0 * dx);
				   B.z[i][n - 1] = B.z[i][n - 1] - c * dt * ((E.y[i + 1][n - 1] - E.y[m - 1][n - 1]) / (2.0 * dx) - (E.x[i][0] - E.x[i][n - 2]) / (2.0 * dy));
			   }
			   if (i == m - 1) {
				   B.y[i][n - 1] = B.y[i][n - 1] + c * dt * (E.z[0][n - 1] - E.z[i - 1][n - 1]) / (2.0 * dx);
				   B.z[i][n - 1] = B.z[i][n - 1] - c * dt * ((E.y[0][n - 1] - E.y[i - 1][n - 1]) / (2.0 * dx) - (E.x[i][0] - E.x[i][n - 2]) / (2.0 * dy));
			   }
		   }

		   for (int j = 1; j < n - 1; j++) {  // i==0
			   B.x[0][j] = B.x[0][j] - c * dt * (E.z[0][j + 1] - E.z[0][j - 1]) / (2.0 * dy);
			   B.y[0][j] = B.y[0][j] + c * dt * (E.z[1][j] - E.z[m - 1][j]) / (2.0 * dx);
			   B.z[0][j] = B.z[0][j] - c * dt * ((E.y[1][j] - E.y[m - 1][j]) / (2.0 * dx) - (E.x[0][j + 1] - E.x[0][j - 1]) / (2.0 * dy));
		   }

		   for (int j = 1; j < n - 1; j++) { // i==nx-1
			   B.x[m - 1][j] = B.x[m - 1][j] - c * dt * (E.z[m - 1][j + 1] - E.z[m - 1][j - 1]) / (2.0 * dy);
			   B.y[m - 1][j] = B.y[m - 1][j] + c * dt * (E.z[0][j] - E.z[m - 2][j]) / (2.0 * dx);
			   B.z[m - 1][j] = B.z[m - 1][j] - c * dt * ((E.y[0][j] - E.y[m - 2][j]) / (2.0 * dx) - (E.x[m - 1][j + 1] - E.x[m - 1][j - 1]) / (2.0 * dy));
		   }*/

		// обновление Ex

		for (unsigned int i = 0; i < m; ++i) {
			E.x[i][0] = E.x[i][0] - 4.0 * PI * dt * J.x[i][0] + c * dt * ((B.z[i][1] - B.z[i][n - 1]) / (2.0 * dy));
			E.x[i][n - 1] = E.x[i][n - 1] - 4.0 * PI * dt * J.x[i][n - 1] + c * dt * ((B.z[i][0] - B.z[i][n - 2]) / (2.0 * dy));
			//std::cout << "i= " << i << std::endl;
			//std::cout << std::endl;
			//std::cout << "Ex[i][0]= " << E.x[i][0] << std::endl;
			//std::cout << "Ex[i][n-1]= " << E.x[i][n - 1] << std::endl;
			for (unsigned int j = 1; j < n - 1; ++j) {
				//std::cout << "i= " << i << ", j= " << j << std::endl;
				//std::cout << std::endl;
				E.x[i][j] = E.x[i][j] - 4.0 * PI * dt * J.x[i][j] + c * dt * ((B.z[i][j + 1] - B.z[i][j - 1]) / (2.0 * dy));
				//std::cout << "Ex[i][j]= " << E.Ex[i][j] << std::endl;
			}
		}

		// обновление Ey

		for (unsigned int j = 0; j < n; ++j) {
			E.y[0][j] = E.y[0][j] - 4.0 * PI * dt * J.y[0][j] - c * dt * ((B.z[1][j] - B.z[m - 1][j]) / (2.0 * dx));
			E.y[m - 1][j] = E.y[m - 1][j] - 4.0 * PI * dt * J.y[m - 1][j] - c * dt * ((B.z[0][j] - B.z[m - 2][j]) / (2.0 * dx));
			//std::cout << "j= " << j << std::endl;
			//std::cout << std::endl;
			//std::cout << "Ey[0][j]= " << E.y[0][j] << std::endl;
			//std::cout << "Ey[m-1][j]= " << E.y[m - 1][j] << std::endl;
			for (unsigned int i = 1; i < m - 1; ++i) {
				//std::cout << "i= " << i << ", j= " << j << std::endl;
				//std::cout << std::endl;
				E.y[i][j] = E.y[i][j] - 4.0 * PI * dt * J.y[i][j] - c * dt * ((B.z[i + 1][j] - B.z[i - 1][j]) / (2.0 * dx));
				//std::cout << "Ey[i][j]= " << E.y[i][j] << std::endl;
			}
		}

		// обновление Ez

		E.z[0][0] = E.z[0][0] - 4.0 * PI * dt * J.z[0][0] + c * dt * \
			((B.y[1][0] - B.y[m - 1][0]) / (2.0 * dx) - (B.x[0][1] - B.x[0][n - 1]) / (2.0 * dy));
		E.z[0][n - 1] = E.z[0][n - 1] - 4.0 * PI * dt * J.z[0][n - 1] + c * dt * \
			((B.y[1][n - 1] - B.y[m - 1][n - 1]) / (2.0 * dx) - (B.x[0][0] - B.x[0][n - 2]) / (2.0 * dy));
		E.z[m - 1][0] = E.z[m - 1][0] - 4.0 * PI * dt * J.z[m - 1][0] + c * dt * \
			((B.y[0][0] - B.y[m - 2][0]) / (2.0 * dx) - (B.x[m - 1][1] - B.x[m - 1][n - 1]) / (2.0 * dy));
		E.z[m - 1][n - 1] = E.z[m - 1][n - 1] - 4.0 * PI * dt * J.z[m - 1][n - 1] + c * dt * \
			((B.y[0][n - 1] - B.y[m - 2][n - 1]) / (2.0 * dx) - (B.x[m - 1][0] - B.x[m - 1][n - 2]) / (2.0 * dy));

		for (unsigned int i = 1; i < m - 1; ++i) {
			E.z[i][0] = E.z[i][0] - 4.0 * PI * dt * J.z[i][0] + c * dt * \
				((B.y[i + 1][0] - B.y[i - 1][0]) / (2.0 * dx) - (B.x[i][1] - B.x[i][n - 1]) / (2.0 * dy));
			E.z[i][n - 1] = E.z[i][n - 1] - 4.0 * PI * dt * J.z[i][n - 1] + c * dt * \
				((B.y[i + 1][n - 1] - B.y[i - 1][n - 1]) / (2.0 * dx) - (B.x[i][0] - B.x[i][n - 2]) / (2.0 * dy));
			//std::cout << "Ez[i][0]= " << E.Ez[i][0] << std::endl;
			//std::cout << "Ez[i][n - 1]= " << E.Ez[i][n - 1] << std::endl;
			//std::cout << std::endl;
		}

		for (unsigned int j = 1; j < n - 1; ++j) {
			E.z[0][j] = E.z[0][j] - 4.0 * PI * dt * J.z[0][j] + c * dt * \
				((B.y[1][j] - B.y[m - 1][j]) / (2.0 * dx) - (B.x[0][j + 1] - B.x[0][j - 1]) / (2.0 * dy));
			E.z[m - 1][j] = E.z[m - 1][j] - 4.0 * PI * dt * J.z[m - 1][j] + c * dt * \
				((B.y[0][j] - B.y[m - 2][j]) / (2.0 * dx) - (B.x[m - 1][j + 1] - B.x[m - 1][j - 1]) / (2.0 * dy));
			//std::cout << "Ez[0][j]= " << E.Ez[0][j] << std::endl;
			//std::cout << "Ez[m-1][j]= " << E.Ez[m - 1][j] << std::endl;
			//std::cout << std::endl;

		}

		for (unsigned int i = 1; i < m - 1; ++i)
			for (unsigned int j = 1; j < n - 1; ++j) {
				E.z[i][j] = E.z[i][j] - 4.0 * PI * dt * J.z[i][j] + c * dt * \
					((B.y[i + 1][j] - B.y[i - 1][j]) / (2.0 * dx) - (B.x[i][j + 1] - B.x[i][j - 1]) / (2.0 * dy));
				//std::cout << "i= " << i << ", j= " << j << std::endl;
				//std::cout << std::endl;
				//std::cout << "Ez[i][j]= " << E.Ez[i][j] << std::endl;
				//std::cout << std::endl;
			}

		// обновление Bx

		for (unsigned int i = 0; i < m; ++i) {
			B.x[i][0] =B.x[i][0] - c * dt * ((E.z[i][1] - E.z[i][n - 1]) / (2.0 * dy));
			B.x[i][n - 1] = B.x[i][n - 1] - c * dt * ((E.z[i][0] - E.z[i][n - 2]) / (2.0 * dy));
			//std::cout << "i= " << i << std::endl;
			//std::cout << std::endl;
			//std::cout << "Bx[i][0]= " << B.Bx[i][0] << std::endl;
			//std::cout << "Bx[i][n-1]= " << B.Bx[i][n - 1] << std::endl;
			for (unsigned int j = 1; j < n - 1; ++j) {
				B.x[i][j] = B.x[i][j] - c * dt * ((E.z[i][j + 1] - E.z[i][j - 1]) / (2.0 * dy));
				//std::cout << "i= " << i << ", j= " << j << std::endl;
				//std::cout << std::endl;
				//std::cout << "Bx[i][j]= " << B.Bx[i][j] << std::endl;
			}
		}

		// обновление By

		for (unsigned int j = 0; j < n; ++j) {
			B.y[0][j] = B.y[0][j] + c * dt * ((E.z[1][j] - E.z[m - 1][j]) / (2.0 * dx));
			B.y[m - 1][j] = B.y[m - 1][j] + c * dt * ((E.z[0][j] - E.z[m - 2][j]) / (2.0 * dx));
			//std::cout << "j= " << j << std::endl;
			//std::cout << std::endl;
			//std::cout << "By[0][j]= " << B.By[0][j] << std::endl;
			//std::cout << "By[m-1][j]= " << B.By[m - 1][j] << std::endl;
			for (unsigned int i = 1; i < m - 1; ++i) {
				//std::cout << "i= " << i << ", j= " << j << std::endl;
				//std::cout << std::endl;
				B.y[i][j] = B.y[i][j] + c * dt * ((E.z[i + 1][j] - E.z[i - 1][j]) / (2.0 * dx));
				//std::cout << "By[i][j]= " << B.By[i][j] << std::endl;
			}
		}

		// обновление Bz

		B.z[0][0] = B.z[0][0] - c * dt * \
			((E.y[1][0] - E.y[m - 1][0]) / (2.0 * dx) - (E.x[0][1] - E.x[0][n - 1]) / (2.0 * dy));
		//std::cout << "Bz[0][0]= " << B.Bz[0][0] << std::endl;
		//std::cout << std::endl;
		B.z[0][n - 1] = B.z[0][n - 1] - c * dt * \
			((E.y[1][n - 1] - E.y[m - 1][n - 1]) / (2.0 * dx) - (E.x[0][0] - E.x[0][n - 2]) / (2.0 * dy));
		//std::cout << "Bz[0][n - 1]= " << B.Bz[0][n - 1] << std::endl;
		//std::cout << std::endl;
		B.z[m - 1][0] = B.z[m - 1][0] - c * dt * \
			((E.y[0][0] - E.y[m - 2][0]) / (2.0 * dx) - (E.x[m - 1][1] - E.x[m - 1][n - 1]) / (2.0 * dy));
		//std::cout << "Bz[m - 1][0]= " << B.Bz[m - 1][0] << std::endl;
		//std::cout << std::endl;
		B.z[m - 1][n - 1] = B.z[m - 1][n - 1] - c * dt * \
			((E.y[0][n - 1] - E.y[m - 2][n - 1]) / (2.0 * dx) - (E.x[m - 1][0] - E.x[m - 1][n - 2]) / (2.0 * dy));
		//std::cout << "Bz[m - 1][n - 1]= " << B.Bz[m - 1][n - 1] << std::endl;
		//std::cout << std::endl;

		for (unsigned int i = 1; i < m - 1; ++i) {
			B.z[i][0] = B.z[i][0] - c * dt * \
				((E.y[i + 1][0] - E.y[i - 1][0]) / (2.0 * dx) - (E.x[i][1] - E.x[i][n - 2]) / (2.0 * dy));
			B.z[i][n - 1] = B.z[i][n - 1] - c * dt * \
				((E.y[i + 1][n - 1] - E.y[i - 1][n - 1]) / (2.0 * dx) - (E.x[i][0] - E.x[i][n - 2]) / (2.0 * dy));
			//std::cout << "i= " << i << std::endl;
			//std::cout << "Bz[i][0]= " << B.Bz[i][0] << std::endl;
			//std::cout << "Bz[i][n-1]= " << B.Bz[i][n - 1] << std::endl;
			//std::cout << std::endl;
		}

		for (unsigned int j = 1; j < n - 1; ++j) {
			B.z[0][j] = B.z[0][j] - c * dt * \
				((E.y[1][j] - E.y[m - 1][j]) / (2.0 * dx) - (E.x[0][j + 1] - E.x[0][j - 1]) / (2.0 * dy));
			B.z[m - 1][j] = B.z[m - 1][j] - c * dt * \
				((E.y[0][j] - E.y[m - 2][j]) / (2.0 * dx) - (E.x[m - 1][j + 1] - E.x[m - 1][j - 1]) / (2.0 * dy));
			//std::cout << "j= " << j << std::endl;
			//std::cout << "Bz[i][j]= " << B.Bz[0][j] << std::endl;
			//std::cout << "Bz[i][j]= " << B.Bz[m - 1][j] << std::endl;
			//std::cout << std::endl;
		}

		for (unsigned int i = 1; i < m - 1; ++i)
			for (unsigned int j = 1; j < n - 1; ++j) {
				B.z[i][j] = B.z[i][j] - c * dt * \
					((E.y[i + 1][j] - E.y[i - 1][j]) / (2.0 * dx) - (E.x[i][j + 1] - E.x[i][j - 1]) / (2.0 * dy));
				//std::cout << " i= " << i << " j= " << j << std::endl;
				//std::cout << std::endl;
				//std::cout << "Bz[i][j]= " << B.z[i][j] << std::endl;
				//std::cout << std::endl;

			}
		return (*this);
	}
};

