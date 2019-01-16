#include <iostream>
#include <conio.h>
#include <stdio.h>
#include <fstream>

using namespace std;

/*----------------------------------------*/

double T[51][21];
double xi[51][21];
double omega[51][21];
double dt = 0.0000001;
double stepx = 0.08;
double stepy = 0.1;
double S = 0.05;

/*----------------------------------------*/

void Conditions()
{
	

	for (int j = 5; j < 15; j++)		// window
	{
		T[0][j] = -15;
	}

	for (int j = 15; j < 21; j++)		// battery
	{
		T[0][j] = 60;
	}

	for (int j = 18; j < 20; j++)		// in
	{
		T[50][j] = 20;
		omega[50][j] = 0;
		xi[50][j] = (-j / 18) * (S / 0.2);
	}

	for (int j = 0; j < 4; j++)			// out
	{
		T[50][j] = 20;
		omega[50][j] = 0;
		xi[50][j] = (50 - j)*(S / 0.2) / (50 - 3);
	}
}

/*----------------------------------------*/

void compute_T_new()
{
	double** Tnew = new double*[51];
	for (int i = 0; i < 51; i++)
		Tnew[i] = new double[21];
	for (int j = 0; j < 21; j++)
	{
		for (int i = 0; i < 51; i++)
		{
			if ((i == 0) || (i == 50) || (j == 0) || (j == 20))
			{
				Tnew[i][j] = T[i][j];
			}
			else
			{
				Tnew[i][j] = T[i][j] + dt * (
					((T[i - 1][j] - 2 * T[i][j] + T[i + 1][j]) / (stepx*stepx)) +
					((T[i][j - 1] - 2 * T[i][j] + T[i][j + 1]) / (stepy*stepy)) -
					(((xi[i][j + 1] - xi[i][j]) / stepy) * ((T[i + 1][j] - T[i][j]) / stepx)) +
					(((xi[i + 1][j] - xi[i][j]) / stepx) * ((T[i][j + 1] - T[i][j]) / stepy))
					);
			}
		}
	}
	for (int i = 0; i < 51; i++)
	{
		for (int j = 0; j < 21; j++)
		{
			T[i][j] = Tnew[i][j];
		}
	}
	for (int i = 0; i < 51; i++)
		delete[] Tnew[i];
	//delete Tnew;
}

/*----------------------------------------*/

void compute_omega_new()
{
	double** omeganew = new double*[51];
	for (int i = 0; i < 51; i++)
		omeganew[i] = new double[21];
	for (int j = 1; j < 20; j++)
	{
		for (int i = 1; i < 50; i++)
		{
			double u = (xi[i][j + 1] - xi[i][j - 1]) / (2 * stepy);
			double v = -(xi[i + 1][j] - xi[i - 1][j]) / (2 * stepx);
			double I = 10000 * (T[i][j + 1] - T[i][j - 1]) / (2 * stepx);
			double Lx = (abs(u) - u) * (omega[i + 1][j] - omega[i][j]) / (2 * stepy) -
						(abs(u) - u) * (omega[i][j] - omega[i - 1][j]) / (2 * stepx) +
						((omega[i + 1][j] - 2 * omega[i][j] + omega[i - 1][j]) / (stepx*stepx)) *
						(abs(0.71 - (abs(u) / 2)*stepx) + 0.71 - (abs(u) / 2)*stepy);
			
			double Ly = (abs(v) - v) * (omega[i][j + 1] - omega[i][j]) / (2 * stepy) -
						(abs(v) + v) * (omega[i][j] - omega[i][j - 1]) / (2 * stepy) +
						((omega[i][j+1] - 2 * omega[i][j] + omega[i][j-1]) / (stepy*stepy)) *
						(abs(0.71 - (abs(v) / 2)*stepy) + 0.71 - (abs(v) / 2)*stepy);

			omeganew[i][j] = omega[i][j] + dt*(Lx + Ly + I);
		}
	}
	for (int i = 0; i < 51; i++)
	{
		for (int j = 0; j < 21; j++)
		{
			omega[i][j] = omeganew[i][j];
		}
	}
	for (int i = 0; i < 51; i++)
		delete[] omeganew[i];
	//delete omeganew;
}

/*----------------------------------------*/

void compute_xi_new()
{
	double err;
	do
	{
		double** xinew = new double*[51];
		for (int i = 0; i < 51; i++)
			xinew[i] = new double[21];
		err = 0;
		double E = 0;
		double F = 0;
		
		for (int j = 1; j < 20; j++)
		{
			for (int i = 1; i < 50; i++)
			{
				double A = 1 / (stepx*stepx);
				double B = 2 / (stepx*stepx) + (1 / dt);
				double C = 1 / (stepx*stepx);
				double D = -(xi[i][j] / dt) - ((xi[i][j + 1] - (2 * xi[i][j]) + xi[i][j - 1]) / (stepy*stepy)) + omega[i][j];

				if (i == 1)
				{
					F = (D - A*xi[i - 1][j]) / (A * 0 - B);
					E = C / (B - A * 0);
				}
				else
				{
					F = (D - A*F) / (A*E - B);
					E = C / (B - A*E);
				}
				xinew[i][j] = E*xi[i + 1][j] + F;
				err += xi[i][j] - xinew[i][j];
				cout << '\t' << err << endl;
			}
		}

		for (int i = 0; i < 51; i++)
		{
			for (int j = 0; j < 21; j++)
			{
				xi[i][j] = xinew[i][j];
			}
		}
		//delete xinew;
		//for (int i = 0; i < 51; i++)
		//	delete[] xinew[i];
	} while (err > 0.001);
	//delete xinew;
}

/*----------------------------------------*/

void insertIntoFile(char sym)
{
	ofstream fout;
	switch (sym)
	{
	case 'T': {
		fout.open("T.txt");
		for (int j = 0; j < 21; j++)
		{
			for (int i = 0; i < 51; i++)
			{
				fout << T[i][j] << '\t';
			}
			fout << endl;
		}
		break;
	}
	case 'X': {
		fout.open("xi.txt");
		for (int j = 0; j < 21; j++)
		{
			for (int i = 0; i < 51; i++)
			{
				fout << xi[i][j] << '\t';
			}
			fout << endl;
		}
		break;
	}
	case 'W': {
		fout.open("omega.txt");
		for (int j = 0; j < 21; j++)
		{
			for (int i = 0; i < 51; i++)
			{
				fout << omega[i][j] << '\t';
			}
			fout << endl;
		}
		break;
	}
	default:
		break;
	}
}

/*----------------------------------------*/

int main()
{
	cout << "program has been started..." << endl;
	int i = 0;

	for (int i = 0; i < 51; i++)
	{
		for (int j = 0; j < 21; j++)
		{
			T[i][j] = 1;
			xi[i][j] = 0;
			omega[i][j] = 0;
		}
	}

	Conditions();

	while (i < 10000)
	{
		
		compute_T_new();
		compute_omega_new();
		compute_xi_new();
		
		i++;
		cout << i << endl;
	}

	for (int j = 0; j < 21; j++)
	{
		for (int i = 0; i < 51; i++)
		{
			cout << T[i][j] << '\t';
		}
		cout << endl;
	}

	insertIntoFile('T');
	insertIntoFile('X');
	insertIntoFile('W');
}