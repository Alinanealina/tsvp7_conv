#include <iostream>
#include <complex>
using namespace std;
const int N = 2, M = 2 * N - 1;
int L = 0;
const float pi = 3.1415926535;
complex<float> i(0, 1);
int T = 0;
void ob(complex<float> A[N], complex<float> B[N])
{
	complex<float> C[M] = { 0 };
	int j, k, l;
	T = 0;
	cout << "\n\n Обычная свертка:";
	for (j = 0; j < M; j++)
	{
		for (k = 0; k < N; k++)
		{
			for (l = 0; l < N; l++)
			{
				if (k + l == j)
				{
					C[j] += A[k] * B[l];
					if ((A[k] * B[l]) != (float)0)
						T += 2;
				}
			}
		}
		//if (C[j].real() != 0 || C[j].imag() != 0)
			//cout << " " << C[j];
	}
	for (j = 0; j < L + N - 1; j++)
		cout << " " << C[j];
	cout << "\n Количество операций: " << T;
}

void DPF(bool f, complex<float> A[2 * N])
{
	complex<float> A1[2 * N];
	int j, k;
	float j2, k2, n = 2 * N;
	for (j = 0; j < 2 * N; j++)
	{
		A1[j] = A[j];
		A[j] = 0;
	}
	if (f)
	{
		for (j = 0; j < 2 * N; j++)
		{
			j2 = j;
			for (k = 0; k < 2 * N; k++)
			{
				k2 = k;
				A[j] += exp(-2 * pi * i * j2 * k2 / n) * A1[k];
				T += 5;
			}
			A[j] /= n;
		}
	}
	else
	{
		for (j = 0; j < 2 * N; j++)
		{
			j2 = j;
			for (k = 0; k < 2 * N; k++)
			{
				k2 = k;
				A[j] += exp(2 * pi * i * j2 * k2 / n) * A1[k];
				T += 5;
			}
		}
	}
}
void DPF_sver(complex<float> A[N], complex<float> B[N])
{
	int j;
	complex<float> A1[2 * N] = { 0 }, B1[2 * N] = { 0 }, C[2 * N] = { 0 };
	T = 0;
	for (j = 0; j < N; j++)
	{
		A1[j] = A[j];
		B1[j] = B[j];
	}
	cout << "\n ДПФ свертка:";
	DPF(true, A1);
	DPF(true, B1);
	for (j = 0; j < 2 * N; j++)
	{
		C[j] = (float)(2 * N) * A1[j] * B1[j];
		T++;
	}
	DPF(false, C);
	for (j = 0; j < L + N - 1; j++)
	{
		//if ((int)(C[j].real() * 100000) != 0 || (int)(C[j].imag() * 100000) != 0)
			printf("\n (%0.10f, %0.10f)", C[j].real(), C[j].imag());
	}
	cout << "\n\n Количество операций: " << T;
}



void PBPF(bool f, complex<float> A[2 * N])
{
	complex<float> X[2 * N], A1[N][N] = { 0 }, A2[N][N] = { 0 }, f1[N][N] = { 0 };
	int M = 0, p1 = 1, p2 = N * 2, j1, j2, k1, k2, j, k;
	float n = 2 * N, fp1, fp2, fj1, fj2, fk1, fk2;
	for (p1 = 1; p1 < p2; p1++)
	{
		if ((N * 2) % p1 == 0)
		{
			p2 = (N * 2) / p1;
			j = p1;
		}
	}
	p1 = j;
	fp1 = p1;
	fp2 = p2;
	for (j = 0; j < 2 * N; j++)
	{
		X[j] = A[j];
		A[j] = 0;
	}
	if (f)
	{
		for (k1 = 0; k1 < p1; k1++)
		{
			fk1 = k1;
			for (j2 = 0; j2 < p2; j2++)
			{
				fj2 = j2;
				for (j1 = 0; j1 < p1; j1++)
				{
					fj1 = j1;
					A1[k1][j2] += X[j2 + p2 * j1] * exp(-2 * pi * i * fj1 * fk1 / fp1);
					T++;
				}
				A1[k1][j2] /= fp1;
				//printf("\n (%0.10f, %0.10f)", A1[k1][j2].real(), A1[k1][j2].imag());
			}
		}
		for (k1 = 0; k1 < p1; k1++)
		{
			fk1 = k1;
			for (k2 = 0; k2 < p2; k2++)
			{
				fk2 = k2;
				for (j2 = 0; j2 < p2; j2++)
				{
					fj2 = j2;
					A2[k1][k2] += A1[k1][j2] * exp(-2 * pi * i * ((fk1 + fp1 * fk2) * fj2 / (fp1 * fp2)));
					T++;
				}
				A2[k1][k2] /= fp2;
			}
		}
		for (j = 0; j < p2; j++)
		{
			for (k = 0; k < p1; k++)
				A[j * p1 + k] = A2[k][j];
		}
	}
	else
	{
		for (j = 0; j < p2; j++)
		{
			for (k = 0; k < p1; k++)
				A2[k][j] = X[j * p1 + k];
		}
		for (k1 = 0; k1 < p1; k1++)
		{
			fk1 = k1;
			for (k2 = 0; k2 < p2; k2++)
			{
				fk2 = k2;
				for (j2 = 0; j2 < p2; j2++)
				{
					fj2 = j2;
					f1[k1][j2] += A2[k1][k2] * exp(2 * pi * i * ((fk1 + fp1 * fk2) * fj2 / (fp1 * fp2)));
					T++;
				}
			}
		}
		for (k1 = 0; k1 < p1; k1++)
		{
			fk1 = k1;
			for (j2 = 0; j2 < p2; j2++)
			{
				fj2 = j2;
				for (j1 = 0; j1 < p1; j1++)
				{
					fj1 = j1;
					A[k1 * p2 + j2] += f1[j1][j2] * exp(2 * pi * i * fj1 * fk1 / fp1);
					T++;
				}
			}
		}
	}
}
void PBPF_sver(complex<float> A[N], complex<float> B[N])
{
	int j;
	complex<float> A1[2 * N] = { 0 }, B1[2 * N] = { 0 }, C[2 * N] = { 0 };
	T = 0;
	for (j = 0; j < N; j++)
	{
		A1[j] = A[j];
		B1[j] = B[j];
	}
	cout << "\n ПБПФ свертка:";
	PBPF(true, A1);
	PBPF(true, B1);
	for (j = 0; j < 2 * N; j++)
	{
		C[j] = (float)(2 * N) * A1[j] * B1[j];
		T++;
	}
	PBPF(false, C);
	for (j = 0; j < L + N - 1; j++)
	{
		//if ((int)(C[j].real() * 100000) != 0 || (int)(C[j].imag() * 100000) != 0)
			printf("\n (%0.10f, %0.10f)", C[j].real(), C[j].imag());
	}
	cout << "\n\n Количество операций: " << T;
}

int main()
{
	setlocale(LC_ALL, "Russian");
	int j;
	complex<float> A[N] = { 1, 2 }, B[N] = { 3, 4 };
	for (j = 0; (A[j].real() != 0 || A[j].imag() != 0) && j < N; j++)
		L++;
	cout << " Исходные массивы:\n A:";
	for (j = 0; j < N; j++)
		cout << " " << A[j];
	cout << "\n B:";
	for (j = 0; j < N; j++)
		cout << " " << B[j];
	ob(A, B);
	cout << endl << "_______________";
	DPF_sver(A, B);
	cout << endl << "_______________";
	PBPF_sver(A, B);
	return 0;
}