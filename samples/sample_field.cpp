#include "elec_magn_field.h"

double f(double x, double y = 0.0, double z = 0.0) {
	estimated_area Area;
	return double(sin(2 * PI * (x - Area.a) / (Area.b - Area.a)));
}

double f_real_solution(double t, double x, double y = 0.0, double z = 0.0) {
	estimated_area Area;
	return double(sin(2 * PI * (x - Area.a - c * t) / (Area.b - Area.a)));
}


double f2( double y, double x = 0.0, double z = 0.0) {
	estimated_area Area;
	return double(sin(2 * PI * (y - Area.a) / (Area.b - Area.a)));
}

double f_real_solution2(double t,  double y, double x = 0.0, double z = 0.0) {
	estimated_area Area;
	return double(sin(2 * PI * (y - Area.a - c * t) / (Area.b - Area.a)));
}


int main()
{
	// сферическая волна

	double t1;
	setlocale(LC_ALL, "Rus");
	std::cout << "Введите время t1: ";
	std::cin >> t1;
	std::cout << std::endl;
	std::ofstream fout("spherical_wave2.txt");
	

	elec_magn_field<> f;

	//double* Ex_for_dft = (double*)malloc(sizeof(double) * f.nx * f.ny * f.nz);
	//double* Ey_for_dft = (double*)malloc(sizeof(double) * f.nx * f.ny * f.nz);
	//double* Ez_for_dft = (double*)malloc(sizeof(double) * f.nx * f.ny * f.nz);
	//double* Bx_for_dft = (double*)malloc(sizeof(double) * f.nx * f.ny * f.nz);
	//double* By_for_dft = (double*)malloc(sizeof(double) * f.nx * f.ny * f.nz);
	//double* Bz_for_dft = (double*)malloc(sizeof(double) * f.nx * f.ny * f.nz);
	//double* Jx_for_dft = (double*)malloc(sizeof(double) * f.nx * f.ny * f.nz);
	//double* Jy_for_dft = (double*)malloc(sizeof(double) * f.nx * f.ny * f.nz);
	//double* Jz_for_dft = (double*)malloc(sizeof(double) * f.nx * f.ny * f.nz);
	//int q = 0;
	//for (int i = 0; i < f.nx; i++) {
	//	for (int j = 0; j < f.ny; j++) {
	//		for (int k = 0; k < f.nz; k++) {
	//			Ex_for_dft[q] = 0;
	//			Ez_for_dft[q] = 0;
	//			Ez_for_dft[q] = 0;
	//			Bx_for_dft[q] = 0;
	//			By_for_dft[q] = 0;
	//			Bz_for_dft[q] = 0;
	//			Jx_for_dft[q] = 0;
	//			Jy_for_dft[q] = 0;
	//			Jz_for_dft[q] = 0;
	//			q++;
	//		}
	//	}
	//}
	//fftw_complex *Bx_out = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * f.nx * f.ny * f.nz);
	//fftw_complex *By_out = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * f.nx * f.ny * f.nz);
	//fftw_complex *Bz_out = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * f.nx * f.ny * f.nz);
	//fftw_complex *Ex_out = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * f.nx * f.ny * f.nz);
	//fftw_complex *Ey_out = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * f.nx * f.ny * f.nz);
	//fftw_complex *Ez_out = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * f.nx * f.ny * f.nz);
	//fftw_complex *Jx_out = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * f.nx * f.ny * f.nz);
	//fftw_complex *Jy_out = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * f.nx * f.ny * f.nz);
	//fftw_complex *Jz_out = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * f.nx * f.ny * f.nz);
	//fftw_plan Ex, Ey, Ez, Bx, By, Bz, Jx, Jy, Jz;
	//fftw_plan Ex_, Ey_, Ez_, Bx_, By_, Bz_, Jx_, Jy_, Jz_;
	//Ex = fftw_plan_dft_r2c_3d(f.nx, f.ny, f.nz, Ex_for_dft, Ex_out, FFTW_ESTIMATE);
	//Ey = fftw_plan_dft_r2c_3d(f.nx, f.ny, f.nz, Ey_for_dft, Ey_out, FFTW_ESTIMATE);
	//Ez = fftw_plan_dft_r2c_3d(f.nx, f.ny, f.nz, Ez_for_dft, Ez_out, FFTW_ESTIMATE);
	//Bx = fftw_plan_dft_r2c_3d(f.nx, f.ny, f.nz, Bx_for_dft, Bx_out, FFTW_ESTIMATE);
	//By = fftw_plan_dft_r2c_3d(f.nx, f.ny, f.nz, By_for_dft, By_out, FFTW_ESTIMATE);
	//Bz = fftw_plan_dft_r2c_3d(f.nx, f.ny, f.nz, Bz_for_dft, Bz_out, FFTW_ESTIMATE);
	//Jx = fftw_plan_dft_r2c_3d(f.nx, f.ny, f.nz, Jx_for_dft, Jx_out, FFTW_ESTIMATE);
	//Jy = fftw_plan_dft_r2c_3d(f.nx, f.ny, f.nz, Jy_for_dft, Jy_out, FFTW_ESTIMATE);
	//Jz = fftw_plan_dft_r2c_3d(f.nx, f.ny, f.nz, Jz_for_dft, Jz_out, FFTW_ESTIMATE);
	//fftw_execute(Ex);
	//fftw_execute(Ey);
	//fftw_execute(Ez);
	//fftw_execute(Bx);
	//fftw_execute(By);
	//fftw_execute(Bz);
	//fftw_execute(Jx);
	//fftw_execute(Jy);
	//fftw_execute(Jz);

	int i = 1;
	int k = int(f.T / f.steps.dt);

	for (double t = f.steps.dt; t <= t1; t += f.steps.dt) {
		if (i <= k) {
			f.modify_Jz(i);

		}
		else {
			f.set_zero_Jz();
		}
		f.PSTD();
		i++;
		/*q = 0;
		for (int _i = 0; _i < f.nx; _i++) {
			for (int _j = 0; _j < f.ny; _j++) {
				for (int _k = 0; _k < f.nz; _k++) {

					Bx_out[q][0] += c * 0.5 * f.steps.dt * f.vect_mult_K(_i, _j, _k, Ex_out[q][1], Ey_out[q][1], Ez_out[q][1]).x;

					Bx_out[q][1] += -c * 0.5 * f.steps.dt * f.vect_mult_K(_i, _j, _k, Ex_out[q][0], Ey_out[q][0], Ez_out[q][0]).x;

					By_out[q][0] += c * 0.5 * f.steps.dt * f.vect_mult_K(_i, _j, _k, Ex_out[q][1], Ey_out[q][1], Ez_out[q][1]).y;

					By_out[q][1] += -c * 0.5 * f.steps.dt * f.vect_mult_K(_i, _j, _k, Ex_out[q][0], Ey_out[q][0], Ez_out[q][0]).y;

					Bz_out[q][0] += c * 0.5 * f.steps.dt * f.vect_mult_K(_i, _j, _k, Ex_out[q][1], Ey_out[q][1], Ez_out[q][1]).z;

					Bz_out[q][1] += -c * 0.5 * f.steps.dt * f.vect_mult_K(_i, _j, _k, Ex_out[q][0], Ey_out[q][0], Ez_out[q][0]).z;

					Ex_out[q][0] += -c * f.steps.dt * f.vect_mult_K(_i, _j, _k, Bx_out[q][1], By_out[q][1], \
						Bz_out[q][1]).x - 4 * PI * f.steps.dt * Jx_out[q][0];

					Ex_out[q][1] += c * f.steps.dt * f.vect_mult_K(_i, _j, _k, Bx_out[q][0], By_out[q][0], \
						Bz_out[q][0]).x - 4 * PI * f.steps.dt * Jx_out[q][1];

					Ey_out[q][0] += -c * f.steps.dt * f.vect_mult_K(_i, _j, _k, Bx_out[q][1], By_out[q][1], \
						Bz_out[q][1]).y - 4 * PI * f.steps.dt * Jy_out[q][0];

					Ey_out[q][1] += c * f.steps.dt * f.vect_mult_K(_i, _j, _k, Bx_out[q][0], By_out[q][0], \
						Bz_out[q][0]).y - 4 * PI * f.steps.dt * Jy_out[q][1];

					Ez_out[q][0] += -c * f.steps.dt * f.vect_mult_K(_i, _j, _k, Bx_out[q][1], By_out[q][1], \
						Bz_out[q][1]).z - 4 * PI * f.steps.dt * Jz_out[q][0];

					Ez_out[q][1] += c * f.steps.dt * f.vect_mult_K(_i, _j, _k, Bx_out[q][0], By_out[q][0], \
						Bz_out[q][0]).x - 4 * PI * f.steps.dt * Jz_out[q][1];

					q++;

				}
			}
		}


		++i;
	}
*/
/*Ex_ = fftw_plan_dft_c2r_3d(f.nx, f.ny, f.nz, Ex_out, Ex_for_dft, FFTW_ESTIMATE);
Ey_ = fftw_plan_dft_c2r_3d(f.nx, f.ny, f.nz, Ey_out, Ey_for_dft, FFTW_ESTIMATE);
Ez_ = fftw_plan_dft_c2r_3d(f.nx, f.ny, f.nz, Ez_out, Ez_for_dft, FFTW_ESTIMATE);
Bx_ = fftw_plan_dft_c2r_3d(f.nx, f.ny, f.nz, Bx_out, Bx_for_dft, FFTW_ESTIMATE);
By_ = fftw_plan_dft_c2r_3d(f.nx, f.ny, f.nz, By_out, By_for_dft, FFTW_ESTIMATE);
Bz_ = fftw_plan_dft_c2r_3d(f.nx, f.ny, f.nz, Bz_out, Bz_for_dft, FFTW_ESTIMATE);
Jx_ = fftw_plan_dft_c2r_3d(f.nx, f.ny, f.nz, Jx_out, Jx_for_dft, FFTW_ESTIMATE);
Jy_ = fftw_plan_dft_c2r_3d(f.nx, f.ny, f.nz, Jy_out, Jy_for_dft, FFTW_ESTIMATE);
Jz_ = fftw_plan_dft_c2r_3d(f.nx, f.ny, f.nz, Jz_out, Jz_for_dft, FFTW_ESTIMATE);

fftw_execute(Ex_);
fftw_execute(Ey_);
fftw_execute(Ez_);
fftw_execute(Bx_);
fftw_execute(By_);
fftw_execute(Bz_);
fftw_execute(Jx_);
fftw_execute(Jy_);
fftw_execute(Jz_);

for (int _i = 0; _i < f.nx * f.ny * f.nz; _i++) {
	Ex_for_dft[_i] /= f.nx * f.ny * f.nz;
	Ey_for_dft[_i] /= f.nx * f.ny * f.nz;
	Ez_for_dft[_i] /= f.nx * f.ny * f.nz;
	Bx_for_dft[_i] /= f.nx * f.ny * f.nz;
	By_for_dft[_i] /= f.nx * f.ny * f.nz;
	Bz_for_dft[_i] /= f.nx * f.ny * f.nz;
	Jx_for_dft[_i] /= f.nx * f.ny * f.nz;
	Jy_for_dft[_i] /= f.nx * f.ny * f.nz;
	Jz_for_dft[_i] /= f.nx * f.ny * f.nz;
}

q = 0;
for (int _i = 0; _i < f.nx; _i++)
	for (int _j = 0; _j < f.ny; _j++)
		for (int _k = 0; _k < f.nz; _k++) {
			f.E.x[_i][_j][_k] = Ex_for_dft[q];
			f.E.y[_i][_j][_k] = Ey_for_dft[q];
			f.E.z[_i][_j][_k] = Ez_for_dft[q];
			f.B.x[_i][_j][_k] = Bx_for_dft[q];
			f.B.y[_i][_j][_k] = By_for_dft[q];
			f.B.z[_i][_j][_k] = Bz_for_dft[q];
			f.J.x[_i][_j][_k] = Jx_for_dft[q];
			f.J.y[_i][_j][_k] = Jy_for_dft[q];
			f.J.z[_i][_j][_k] = Jz_for_dft[q];
			q++;
		}
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





free(Ex_for_dft);
free(Ey_for_dft);
free(Ez_for_dft);
free(Bx_for_dft);
free(By_for_dft);
free(Bz_for_dft);
free(Jx_for_dft);
free(Jy_for_dft);
free(Jz_for_dft);*/
	}

	for (int i = 0; i < f.nx; ++i) {
		for (int j = 0; j < f.ny; ++j)
			fout << f.E.z(i, j, 0) << " ";
		fout << std::endl;
	}
	

	fout.close();
	std::cout << "Вычисления выполнены";
	return 0;
}
