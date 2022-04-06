#include "elec_magn_field.h"

#include <gtest.h>

const double ax = 0, bx = 5E+10, ay = 0, by = 5E+10, az = 0, bz = 0;

TEST(FDTD, correct_result_Ex_for_axis_x)
{
	double real_value = 0;
	double t1 = 0.1;
	double dt = 0.001, dx = 1E+9, dy = 1E+10, dz = 1E+9;
	elec_magn_field<> f(dt, dx, dy, dz, ax, bx, ay, by, az, bz);
	for (int i = 0; i < f.nx; i++)
		for (int j = 0; j < f.ny; j++) 
			{
				f.E.y(i, j, 0) = sin(2 * PI * (f.steps.dx * i) / (f.area.bx - f.area.ax));
				f.B.z(i, j, 0) = sin(2 * PI * (f.steps.dx * i) / (f.area.bx - f.area.ax));
			}
	for (double t = 0; t < t1; t += f.steps.dt)
		f.FDTD(t);
	EXPECT_NEAR(f.E.x(3, 4, 0), real_value, 1e-3);
}

TEST(FDTD, correct_result_Ey_for_axis_x)
{
	double t1 = 0.1;
	double dt = 0.001, dx = 1E+9, dy = 1E+10, dz = 1E+9;
	elec_magn_field<> f(dt, dx, dy, dz, ax, bx, ay, by, az, bz);
	double real_value = sin(2 * PI * (3 * f.steps.dx - c * t1) / (f.area.bx - f.area.ax));
	
	for (int i = 0; i < f.nx; i++)
		for (int j = 0; j < f.ny; j++) {
			f.E.y(i, j, 0) = sin(2 * PI * (f.steps.dx * i) / (f.area.bx - f.area.ax));
			f.B.z(i, j, 0) = sin(2 * PI * (f.steps.dx * i) / (f.area.bx - f.area.ax));
		}
	for (double t = 0; t < t1; t += f.steps.dt)
		f.FDTD(t);
	EXPECT_NEAR(f.E.x(3, 4, 0), real_value, 1e-3);
}

TEST(FDTD, correct_result_Ez_for_axis_x)
{
	double real_value = 0;
	double t1 = 0.1;
	double dt = 0.001, dx = 1E+9, dy = 1E+10, dz = 1E+9;
	elec_magn_field<> f(dt, dx, dy, dz, ax, bx, ay, by, az, bz);
	for (int i = 0; i < f.nx; i++)
		for (int j = 0; j < f.ny; j++) {
			f.E.y(i, j, 0) = sin(2 * PI * (f.steps.dx * i) / (f.area.bx - f.area.ax));
			f.B.z(i, j, 0) = sin(2 * PI * (f.steps.dx * i) / (f.area.bx - f.area.ax));
		}
	for (double t = 0; t < t1; t += f.steps.dt)
		f.FDTD(t);
	EXPECT_NEAR(f.E.x(3, 4, 0), real_value, 1e-3);
}

TEST(FDTD, correct_result_Bx_for_axis_x)
{
	double dt = 0.001, t1 = 0.1, real_value, dx = 1E+9, dy = 1E+10, dz = 1E+9;
	real_value = 0;
	elec_magn_field<> f(dt, dx, dy, dz, ax, bx, ay, by, az, bz);
	for (int i = 0; i < f.nx; i++)
		for (int j = 0; j < f.ny; j++) {
			f.E.y(i, j, 0) = sin(2 * PI * (f.steps.dx * i) / (f.area.bx - f.area.ax));
			f.B.z(i, j, 0) = sin(2 * PI * (f.steps.dx * i) / (f.area.bx - f.area.ax));
		}
	for (double t = 0; t < t1; t += f.steps.dt)
		f.FDTD(t);
	EXPECT_NEAR(f.E.x(3, 4, 0), real_value, 1e-3);
}

TEST(FDTD, correct_result_By_for_axis_x)
{
	double dt = 0.001, t1 = 0.1, real_value, dx = 1E+9, dy = 1E+10, dz = 1E+9;
	real_value = 0;
	elec_magn_field<> f(dt, dx, dy, dz, ax, bx, ay, by, az, bz);
	for (int i = 0; i < f.nx; i++)
		for (int j = 0; j < f.ny; j++) {
			f.E.y(i, j, 0) = sin(2 * PI * (f.steps.dx * i) / (f.area.bx - f.area.ax));
			f.B.z(i, j, 0) = sin(2 * PI * (f.steps.dx * i) / (f.area.bx - f.area.ax));
		}
	for (double t = 0; t < t1; t += f.steps.dt)
		f.FDTD(t);
	EXPECT_NEAR(f.E.x(3, 4, 0), real_value, 1e-3);
}

TEST(FDTD, correct_result_Bz_for_axis_x)
{
	double dt = 0.001, t1 = 0.1, real_value, dx = 1E+9, dy = 1E+10, dz = 1E+9;
	elec_magn_field<> f(dt, dx, dy, dz, ax, bx, ay, by, az, bz);
	real_value = sin(2 * PI * (3 * f.steps.dx - c * t1) / (f.area.bx - f.area.ax));
	for (int i = 0; i < f.nx; i++)
		for (int j = 0; j < f.ny; j++) {
			f.E.y(i, j, 0) = sin(2 * PI * (f.steps.dx * i) / (f.area.bx - f.area.ax));
			f.B.z(i, j, 0) = sin(2 * PI * (f.steps.dx * i) / (f.area.bx - f.area.ax));
		}
	for (double t = 0; t < t1; t += f.steps.dt)
		f.FDTD(t);
	EXPECT_NEAR(f.E.x(3, 4, 0), real_value, 1e-3);
}


TEST(FDTD, correct_result_Ex_for_axis_y)
{
	double dt = 0.001, t1 = 0.1, real_value, dx = 1E+10, dy = 1E+9, dz = 1E+9;
	elec_magn_field<> f(dt, dx, dy, dz, ax, bx, ay, by, az, bz);
	real_value = sin(2.0 * PI * (4 * f.steps.dy - c * t1) / (f.area.by - f.area.ay));
	for (int i = 0; i < f.nx; i++)
		for (int j = 0; j < f.ny; j++) {
			f.E.x(i, j, 0) = sin(2 * PI * (f.steps.dy * j) / (f.area.by - f.area.ay));
			f.B.z(i, j, 0) = -sin(2 * PI * (f.steps.dy * j) / (f.area.by - f.area.ay));
		}
	for (double t = 0; t < t1; t += f.steps.dt)
		f.FDTD(t);
	EXPECT_NEAR(f.E.x(3, 4, 0), real_value, 1e-3);
}

TEST(FDTD, correct_result_Ey_for_axis_y)
{
	double dt = 0.001, t1 = 0.1, real_value, dx = 1E+10, dy = 1E+9, dz = 1E+9;
	real_value = 0;
	elec_magn_field<> f(dt, dx, dy, dz, ax, bx, ay, by, az, bz);
	for (int i = 0; i < f.nx; i++)
		for (int j = 0; j < f.ny; j++) {
			f.E.x(i, j, 0) = sin(2 * PI * (f.steps.dy * j) / (f.area.by - f.area.ay));
			f.B.z(i, j, 0) = -sin(2 * PI * (f.steps.dy * j) / (f.area.by - f.area.ay));
		}
	for (double t = 0; t < t1; t += f.steps.dt)
		f.FDTD(t);
	EXPECT_NEAR(f.E.x(3, 4, 0), real_value, 1e-3);
}

TEST(FDTD, correct_result_Ez_for_axis_y)
{
	double dt = 0.001, t1 = 0.1, real_value, dx = 1E+10, dy = 1E+9, dz = 1E+9;
	real_value = 0;
	elec_magn_field<> f(dt, dx, dy, dz, ax, bx, ay, by, az, bz);
	for (int i = 0; i < f.nx; i++)
		for (int j = 0; j < f.ny; j++) {
			f.E.x(i, j, 0) = sin(2 * PI * (f.steps.dy * j) / (f.area.by - f.area.ay));
			f.B.z(i, j, 0) = -sin(2 * PI * (f.steps.dy * j) / (f.area.by - f.area.ay));
		}
	for (double t = 0; t < t1; t += f.steps.dt)
		f.FDTD(t);
	EXPECT_NEAR(f.E.x(3, 4, 0), real_value, 1e-3);
}

TEST(FDTD, correct_result_Bx_for_axis_y)
{
	double dt = 0.001, t1 = 0.1, real_value, dx = 1E+10, dy = 1E+9, dz = 1E+9;
	real_value = 0;
	elec_magn_field<> f(dt, dx, dy, dz, ax, bx, ay, by, az, bz);
	for (int i = 0; i < f.nx; i++)
		for (int j = 0; j < f.ny; j++) {
			f.E.x(i, j, 0) = sin(2 * PI * (f.steps.dy * j) / (f.area.by - f.area.ay));
			f.B.z(i, j, 0) = -sin(2 * PI * (f.steps.dy * j) / (f.area.by - f.area.ay));
		}
	for (double t = 0; t < t1; t += f.steps.dt)
		f.FDTD(t);
	EXPECT_NEAR(f.E.x(3, 4, 0), real_value, 1e-3);
}

TEST(FDTD, correct_result_By_for_axis_y)
{
	double dt = 0.001, t1 = 0.1, real_value, dx = 1E+10, dy = 1E+9, dz = 1E+9;
	real_value = 0;
	elec_magn_field<> f(dt, dx, dy, dz, ax, bx, ay, by, az, bz);
	for (int i = 0; i < f.nx; i++)
		for (int j = 0; j < f.ny; j++) {
			f.E.x(i, j, 0) = sin(2 * PI * (f.steps.dy * j) / (f.area.by - f.area.ay));
			f.B.z(i, j, 0) = -sin(2 * PI * (f.steps.dy * j) / (f.area.by - f.area.ay));
		}
	for (double t = 0; t < t1; t += f.steps.dt)
		f.FDTD(t);
	EXPECT_NEAR(f.E.x(3, 4, 0), real_value, 1e-3);
}

TEST(FDTD, correct_result_Bz_for_axis_y)
{
	double dt = 0.001, t1 = 0.1, real_value, dx = 1E+10, dy = 1E+9, dz = 1E+9;
	elec_magn_field<> f(dt, dx, dy, dz, ax, bx, ay, by, az, bz);
	real_value = -sin(2.0 * PI * (4 * f.steps.dy - c * t1) / (f.area.by - f.area.ay));
	for (int i = 0; i < f.nx; i++)
		for (int j = 0; j < f.ny; j++) {
			f.E.x(i, j, 0) = sin(2 * PI * (f.steps.dy * j) / (f.area.by - f.area.ay));
			f.B.z(i, j, 0) = -sin(2 * PI * (f.steps.dy * j) / (f.area.by - f.area.ay));
		}
	for (double t = 0; t < t1; t += f.steps.dt)
		f.FDTD(t);
	EXPECT_NEAR(f.E.x(3, 4, 0), real_value, 1e-3);
}

//________________________

TEST(PSTD, correct_result_Ex_for_axis_x)
{
	double real_value = 0;
	double t1 = 0.1;
	double dt = 0.001, dx = 1E+9, dy = 1E+10, dz = 1E+9;
	elec_magn_field<> f(dt, dx, dy, dz, ax, bx, ay, by, az, bz);
	for (int i = 0; i < f.nx; i++)
		for (int j = 0; j < f.ny; j++)
		{
			f.E.y(i, j, 0) = sin(2 * PI * (f.steps.dx * i) / (f.area.bx - f.area.ax));
			f.B.z(i, j, 0) = sin(2 * PI * (f.steps.dx * i) / (f.area.bx - f.area.ax));
		}
	for (double t = 0; t < t1; t += f.steps.dt)
		f.PSTD(t);
	EXPECT_NEAR(f.E.x(3, 4, 0), real_value, 1e-3);
}

TEST(PSTD, correct_result_Ey_for_axis_x)
{
	double t1 = 0.1;
	double dt = 0.001, dx = 1E+9, dy = 1E+10, dz = 1E+9;
	elec_magn_field<> f(dt, dx, dy, dz, ax, bx, ay, by, az, bz);
	double real_value = sin(2 * PI * (3 * f.steps.dx - c * t1) / (f.area.bx - f.area.ax));

	for (int i = 0; i < f.nx; i++)
		for (int j = 0; j < f.ny; j++) {
			f.E.y(i, j, 0) = sin(2 * PI * (f.steps.dx * i) / (f.area.bx - f.area.ax));
			f.B.z(i, j, 0) = sin(2 * PI * (f.steps.dx * i) / (f.area.bx - f.area.ax));
		}
	for (double t = 0; t < t1; t += f.steps.dt)
		f.PSTD(t);
	EXPECT_NEAR(f.E.x(3, 4, 0), real_value, 1e-3);
}

TEST(PSTD, correct_result_Ez_for_axis_x)
{
	double real_value = 0;
	double t1 = 0.1;
	double dt = 0.001, dx = 1E+9, dy = 1E+10, dz = 1E+9;
	elec_magn_field<> f(dt, dx, dy, dz, ax, bx, ay, by, az, bz);
	for (int i = 0; i < f.nx; i++)
		for (int j = 0; j < f.ny; j++) {
			f.E.y(i, j, 0) = sin(2 * PI * (f.steps.dx * i) / (f.area.bx - f.area.ax));
			f.B.z(i, j, 0) = sin(2 * PI * (f.steps.dx * i) / (f.area.bx - f.area.ax));
		}
	for (double t = 0; t < t1; t += f.steps.dt)
		f.PSTD(t);
	EXPECT_NEAR(f.E.x(3, 4, 0), real_value, 1e-3);
}

TEST(PSTD, correct_result_Bx_for_axis_x)
{
	double dt = 0.001, t1 = 0.1, real_value, dx = 1E+9, dy = 1E+10, dz = 1E+9;
	real_value = 0;
	elec_magn_field<> f(dt, dx, dy, dz, ax, bx, ay, by, az, bz);
	for (int i = 0; i < f.nx; i++)
		for (int j = 0; j < f.ny; j++) {
			f.E.y(i, j, 0) = sin(2 * PI * (f.steps.dx * i) / (f.area.bx - f.area.ax));
			f.B.z(i, j, 0) = sin(2 * PI * (f.steps.dx * i) / (f.area.bx - f.area.ax));
		}
	for (double t = 0; t < t1; t += f.steps.dt)
		f.PSTD(t);
	EXPECT_NEAR(f.E.x(3, 4, 0), real_value, 1e-3);
}

TEST(PSTD, correct_result_By_for_axis_x)
{
	double dt = 0.001, t1 = 0.1, real_value, dx = 1E+9, dy = 1E+10, dz = 1E+9;
	real_value = 0;
	elec_magn_field<> f(dt, dx, dy, dz, ax, bx, ay, by, az, bz);
	for (int i = 0; i < f.nx; i++)
		for (int j = 0; j < f.ny; j++) {
			f.E.y(i, j, 0) = sin(2 * PI * (f.steps.dx * i) / (f.area.bx - f.area.ax));
			f.B.z(i, j, 0) = sin(2 * PI * (f.steps.dx * i) / (f.area.bx - f.area.ax));
		}
	for (double t = 0; t < t1; t += f.steps.dt)
		f.PSTD(t);
	EXPECT_NEAR(f.E.x(3, 4, 0), real_value, 1e-3);
}

TEST(PSTD, correct_result_Bz_for_axis_x)
{
	double dt = 0.001, t1 = 0.1, real_value, dx = 1E+9, dy = 1E+10, dz = 1E+9;
	elec_magn_field<> f(dt, dx, dy, dz, ax, bx, ay, by, az, bz);
	real_value = sin(2 * PI * (3 * f.steps.dx - c * t1) / (f.area.bx - f.area.ax));
	for (int i = 0; i < f.nx; i++)
		for (int j = 0; j < f.ny; j++) {
			f.E.y(i, j, 0) = sin(2 * PI * (f.steps.dx * i) / (f.area.bx - f.area.ax));
			f.B.z(i, j, 0) = sin(2 * PI * (f.steps.dx * i) / (f.area.bx - f.area.ax));
		}
	for (double t = 0; t < t1; t += f.steps.dt)
		f.PSTD(t);
	EXPECT_NEAR(f.E.x(3, 4, 0), real_value, 1e-3);
}


TEST(PSTD, correct_result_Ex_for_axis_y)
{
	double dt = 0.001, t1 = 0.1, real_value, dx = 1E+10, dy = 1E+9, dz = 1E+9;
	elec_magn_field<> f(dt, dx, dy, dz, ax, bx, ay, by, az, bz);
	real_value = sin(2.0 * PI * (4 * f.steps.dy - c * t1) / (f.area.by - f.area.ay));
	for (int i = 0; i < f.nx; i++)
		for (int j = 0; j < f.ny; j++) {
			f.E.x(i, j, 0) = sin(2 * PI * (f.steps.dy * j) / (f.area.by - f.area.ay));
			f.B.z(i, j, 0) = -sin(2 * PI * (f.steps.dy * j) / (f.area.by - f.area.ay));
		}
	for (double t = 0; t < t1; t += f.steps.dt)
		f.PSTD(t);
	EXPECT_NEAR(f.E.x(3, 4, 0), real_value, 1e-3);
}

TEST(PSTD, correct_result_Ey_for_axis_y)
{
	double dt = 0.001, t1 = 0.1, real_value, dx = 1E+10, dy = 1E+9, dz = 1E+9;
	real_value = 0;
	elec_magn_field<> f(dt, dx, dy, dz, ax, bx, ay, by, az, bz);
	for (int i = 0; i < f.nx; i++)
		for (int j = 0; j < f.ny; j++) {
			f.E.x(i, j, 0) = sin(2 * PI * (f.steps.dy * j) / (f.area.by - f.area.ay));
			f.B.z(i, j, 0) = -sin(2 * PI * (f.steps.dy * j) / (f.area.by - f.area.ay));
		}
	for (double t = 0; t < t1; t += f.steps.dt)
		f.PSTD(t);
	EXPECT_NEAR(f.E.x(3, 4, 0), real_value, 1e-3);
}

TEST(PSTD, correct_result_Ez_for_axis_y)
{
	double dt = 0.001, t1 = 0.1, real_value, dx = 1E+10, dy = 1E+9, dz = 1E+9;
	real_value = 0;
	elec_magn_field<> f(dt, dx, dy, dz, ax, bx, ay, by, az, bz);
	for (int i = 0; i < f.nx; i++)
		for (int j = 0; j < f.ny; j++) {
			f.E.x(i, j, 0) = sin(2 * PI * (f.steps.dy * j) / (f.area.by - f.area.ay));
			f.B.z(i, j, 0) = -sin(2 * PI * (f.steps.dy * j) / (f.area.by - f.area.ay));
		}
	for (double t = 0; t < t1; t += f.steps.dt)
		f.PSTD(t);
	EXPECT_NEAR(f.E.x(3, 4, 0), real_value, 1e-3);
}

TEST(PSTD, correct_result_Bx_for_axis_y)
{
	double dt = 0.001, t1 = 0.1, real_value, dx = 1E+10, dy = 1E+9, dz = 1E+9;
	real_value = 0;
	elec_magn_field<> f(dt, dx, dy, dz, ax, bx, ay, by, az, bz);
	for (int i = 0; i < f.nx; i++)
		for (int j = 0; j < f.ny; j++) {
			f.E.x(i, j, 0) = sin(2 * PI * (f.steps.dy * j) / (f.area.by - f.area.ay));
			f.B.z(i, j, 0) = -sin(2 * PI * (f.steps.dy * j) / (f.area.by - f.area.ay));
		}
	for (double t = 0; t < t1; t += f.steps.dt)
		f.PSTD(t);
	EXPECT_NEAR(f.E.x(3, 4, 0), real_value, 1e-3);
}

TEST(PSTD, correct_result_By_for_axis_y)
{
	double dt = 0.001, t1 = 0.1, real_value, dx = 1E+10, dy = 1E+9, dz = 1E+9;
	real_value = 0;
	elec_magn_field<> f(dt, dx, dy, dz, ax, bx, ay, by, az, bz);
	for (int i = 0; i < f.nx; i++)
		for (int j = 0; j < f.ny; j++) {
			f.E.x(i, j, 0) = sin(2 * PI * (f.steps.dy * j) / (f.area.by - f.area.ay));
			f.B.z(i, j, 0) = -sin(2 * PI * (f.steps.dy * j) / (f.area.by - f.area.ay));
		}
	for (double t = 0; t < t1; t += f.steps.dt)
		f.PSTD(t);
	EXPECT_NEAR(f.E.x(3, 4, 0), real_value, 1e-3);
}

TEST(PSTD, correct_result_Bz_for_axis_y)
{
	double dt = 0.001, t1 = 0.1, real_value, dx = 1E+10, dy = 1E+9, dz = 1E+9;
	elec_magn_field<> f(dt, dx, dy, dz, ax, bx, ay, by, az, bz);
	real_value = -sin(2.0 * PI * (4 * f.steps.dy - c * t1) / (f.area.by - f.area.ay));
	for (int i = 0; i < f.nx; i++)
		for (int j = 0; j < f.ny; j++) {
			f.E.x(i, j, 0) = sin(2 * PI * (f.steps.dy * j) / (f.area.by - f.area.ay));
			f.B.z(i, j, 0) = -sin(2 * PI * (f.steps.dy * j) / (f.area.by - f.area.ay));
		}
	for (double t = 0; t < t1; t += f.steps.dt)
		f.PSTD(t);
	EXPECT_NEAR(f.E.x(3, 4, 0), real_value, 1e-3);
}