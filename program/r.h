#pragma once
#include <iostream>
#include <cmath>
#include <string>
#include <vector>
#include <ctime>
#include <omp.h>

#define PI 3.141592653589793

using namespace std;

double T_analytic_func(double A, double B, double R1, double R2, double r, double phi);
double C_norm(vector<vector<double>> a);
double C_error(vector<vector<double>> a, vector<vector<double>> b);
double L_norm(vector<vector<double>> a, double R1, double dr, double dphi);
double L_error(vector<vector<double>> a, vector<vector<double>> b, double R1, double dr, double dphi);

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
	);