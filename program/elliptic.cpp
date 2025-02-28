#include "r.h"

double T_analytic_func(double A, double B, double R1, double R2, double r, double phi)
{
	return B + B*R2/(R1*R1 + R2*R2)*(r + R1*R1/r)*cos(phi) + A*R1*R1*R1*0.5/(R1*R1*R1*R1 + R2*R2*R2*R2)*(r*r - R2*R2*R2*R2/r/r)*sin(2*phi);
}

double C_norm(vector<vector<double>> a)
{
	double norm = 0;
	for (size_t i = 0; i < a.size(); ++i)
        for (size_t j = 0; j < a[i].size(); ++j)
            norm = max(norm, abs(a[i][j]));
    return norm;
}

double C_error(vector<vector<double>> a, vector<vector<double>> b)
{
	double norm = 0;
	for (size_t i = 0; i < a.size(); ++i)
        for (size_t j = 0; j < a[i].size(); ++j)
            norm = max(norm, abs(a[i][j] - b[i][j]));
    return norm;
}

double L_norm(vector<vector<double>> a, double R1, double dr, double dphi)
{
	double norm = 0;
	for(size_t i = 0; i < a.size(); i++)
		for(size_t j = 0; j < a[i].size(); j++)
			norm += abs(a[i][j]) * (R1 + i*dr) * dr * dphi;
	return norm;
}

double L_error(vector<vector<double>> a, vector<vector<double>> b, double R1, double dr, double dphi)
{
	double norm = 0;
	for(size_t i = 0; i < a.size(); i++)
		for(size_t j = 0; j < a[i].size(); j++)
			norm += abs(a[i][j] - b[i][j]) * (R1 + i*dr) * dr * dphi;
	return norm;
}

vector<double> temperature_ring(
	int N,
	int M,
	double A,
	double B,
	double R1,
	double R2,
	double tau,
	string file_analytic,
	string file_numerical,
	string file_error
	)
{
	vector<vector<double>>
		T_prev(N+1, vector<double>(M, 0)),
		T_next(N+1, vector<double>(M, 0));

	double
		dr = (R2 - R1) / N,
		dphi = 2 * PI / M,
		t = 0,
		dt, r, phi,
		radial, polar;

	dt = 0.495 * dr * dr * R1 * R1 * dphi * dphi / (dr*dr + R1*R1*dphi*dphi);

	cout << "dt = " << dt << endl;

	vector<bool> time_flag(100, false), err_flag(10, false);

	//analytic solution

	vector<vector<double>>
	T_analytic(N+1, vector<double>(M, 0));

	for(int i = 0; i < N + 1; i++)
		for(int j = 0; j < M; j++)
		{
			r = R1 + i*dr;
			phi = j*dphi;
			T_analytic[i][j] = T_analytic_func(A, B, R1, R2, r, phi);
		}

	FILE *tmp;

	cout << "t = " << t << endl;

	cout << "C-norm of (absolute) error: " << C_error(T_analytic, T_prev) << endl;
	cout << "L-norm of (absolute) error: " << L_error(T_analytic, T_prev, R1, dr, dphi) << endl;
	cout << "C-norm of (relative) error: " << C_error(T_analytic, T_prev) / C_norm(T_analytic) << endl;
	cout << "L-norm of (relative) error: " << L_error(T_analytic, T_prev, R1, dr, dphi)
											/ L_norm(T_analytic, R1, dr, dphi) << endl;

	tmp = fopen("output_results/numerical_0.txt", "w");
	for(int i = 0; i < N + 1; i++)
	{
		for(int j = 0; j < M; j++)
		{
			fprintf(tmp, "%lf ", T_prev[i][j]);
		}
		fprintf(tmp, "\n");
	}
	fclose(tmp);

	while(t < tau)
	{
		//inner points (simple scheme)

		// #pragma omp parallel for collapse(2) private(radial, polar) shared(T_prev, T_next, dt, R1, dr, dphi, N, M)
        for(int i = 1; i < N; i++)
            for(int j = 0; j < M; j++)
            {
                radial = (R1 + (i+0.5)*dr) * (T_prev[i+1][j] - T_prev[i][j]) / dr;
                radial -= (R1 + (i-0.5)*dr) * (T_prev[i][j] - T_prev[i-1][j]) / dr;
                radial /= (R1 + i*dr) * dr;

                polar = (T_prev[i][(j+1) % M] - 2*T_prev[i][j] + T_prev[i][(j-1+M) % M])
                        / ((R1 + i*dr) * (R1 + i*dr) * dphi * dphi);

                T_next[i][j] = T_prev[i][j] + dt * (radial + polar);
            }

        //boundary conditions

        // #pragma omp parallel for shared(T_prev, T_next, M, dr, A, dphi, B)
		for(int j = 0; j < M; j++)
		{
			T_next[0][j] = T_prev[1][j] - dr*A*sin(2*j*dphi);
			T_next[N][j] = B*(1 + cos(j*dphi));
		}
		
		//vector update

		// #pragma omp parallel for collapse(2) shared(T_prev, T_next, N, M)
        for(int i = 0; i < N + 1; i++)
        {
            for(int j = 0; j < M; j++)
                T_prev[i][j] = T_next[i][j];
        }

		t += dt;

		// progress bar
		for(int l = 0; l < 100; l++)
		{
			if(t > tau * (double(l) + 1) / 100.0 && !time_flag[l])
			{
				time_flag[l] = true;
				printf("\n%d%% completed\n", l+1);
				cout << "t = " << t << endl;

				cout << "C-norm of (absolute) error: " << C_error(T_analytic, T_prev) << endl;
				cout << "L-norm of (absolute) error: " << L_error(T_analytic, T_prev, R1, dr, dphi) << endl;
				cout << "C-norm of (relative) error: " << C_error(T_analytic, T_prev) / C_norm(T_analytic) << endl;
				cout << "L-norm of (relative) error: " << L_error(T_analytic, T_prev, R1, dr, dphi)
														/ L_norm(T_analytic, R1, dr, dphi) << endl;

				tmp = fopen(("output_results/numerical_" + to_string(l+1) + ".txt").c_str(), "w");
				for(int i = 0; i < N + 1; i++)
				{
					for(int j = 0; j < M; j++)
					{
						fprintf(tmp, "%lf ", T_prev[i][j]);
					}
					fprintf(tmp, "\n");
				}
				fclose(tmp);
			}
		}
	}

	

	FILE
		*T_n = fopen(("output_results/" + file_numerical).c_str(), "w"),
		*T_a = fopen(("output_results/" + file_analytic).c_str(), "w"),
		*err = fopen(("output_results/" + file_error).c_str(), "w");
	for(int i = 0; i < N + 1; i++)
	{
		for(int j = 0; j < M; j++)
		{
			fprintf(T_n, "%lf ", T_prev[i][j]);
			fprintf(T_a, "%lf ", T_analytic[i][j]);
			fprintf(err, "%lf ", abs(T_prev[i][j] - T_analytic[i][j]));
		}
		fprintf(T_n, "\n");
		fprintf(T_a, "\n");
		fprintf(err, "\n");
	}
	fclose(T_n);
	fclose(T_a);
	fclose(err);

	vector<double> error(4, 0);
	error[0] = C_error(T_analytic, T_prev);
	error[1] = L_error(T_analytic, T_prev, R1, dr, dphi);
	error[2] = C_error(T_analytic, T_prev) / C_norm(T_analytic);
	error[3] = L_error(T_analytic, T_prev, R1, dr, dphi)
				/ L_norm(T_analytic, R1, dr, dphi);
	return error;
}