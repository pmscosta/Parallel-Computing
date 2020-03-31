//#include <omp.h>
#include <stdio.h>
#include <iostream>
#include <iomanip>
#include <time.h>
#include <cstdlib>
#include <fstream>
#include <papi.h>

using namespace std;

#define SYSTEMTIME clock_t

ofstream timeLogger;

void first_implementation(int m_ar, int m_br)
{

	SYSTEMTIME Time1, Time2;

	char st[100];
	double temp;
	int i, j, k;

	double *pha, *phb, *phc;

	pha = (double *)malloc((m_ar * m_ar) * sizeof(double));
	phb = (double *)malloc((m_ar * m_ar) * sizeof(double));
	phc = (double *)malloc((m_ar * m_ar) * sizeof(double));

	for (i = 0; i < m_ar; i++)
		for (j = 0; j < m_ar; j++)
			pha[i * m_ar + j] = (double)1.0;

	for (i = 0; i < m_br; i++)
		for (j = 0; j < m_br; j++)
			phb[i * m_br + j] = (double)(i + 1);

	Time1 = clock();

	for (i = 0; i < m_ar; i++)
	{
		for (j = 0; j < m_br; j++)
		{
			for (k = 0; k < m_ar; k++)
			{
				phc[i * m_ar + j] += pha[i * m_ar + k] * phb[k * m_br + j];
			}
		}
	}

	Time2 = clock();
	sprintf(st, "%4.4f", (double)(Time2 - Time1) / CLOCKS_PER_SEC);
	cout << "Time: " << st << " seconds\n";

	timeLogger << m_ar << ", " << m_br << ", " << st << ", ";

	cout << "Result matrix: " << endl;
	for (i = 0; i < 1; i++)
	{
		for (j = 0; j < min(10, m_br); j++)
			cout << phc[j] << " ";
	}
	cout << endl;

	free(pha);
	free(phb);
	free(phc);
}

void blockMult(int m_ar, int m_br, int blockSize)
{
	SYSTEMTIME Time1, Time2;

	char st[100];
	double temp;
	int i, j, k;

	double *pha, *phb, *phc;

	pha = (double *)malloc((m_ar * m_ar) * sizeof(double));
	phb = (double *)malloc((m_ar * m_ar) * sizeof(double));
	phc = (double *)malloc((m_ar * m_ar) * sizeof(double));

	for (i = 0; i < m_ar; i++)
		for (j = 0; j < m_ar; j++)
			pha[i * m_ar + j] = (double)1.0;

	for (i = 0; i < m_br; i++)
		for (j = 0; j < m_br; j++)
			phb[i * m_br + j] = (double)(i + 1);

	Time1 = clock();
	for (int i0 = 0; i0 < m_ar; i0 += blockSize)
	{
		for (int j0 = 0; j0 < m_br; j0 += blockSize)
		{
			for (int k0 = 0; k0 < m_ar; k0 += blockSize)
			{
				for (i = 0; i < blockSize; i++)
				{
					for (k = 0; k < blockSize; k++)
					{
						for (j = 0; j < blockSize; j++)
						{
							//cout << pha[(i0+i)*m_ar+k+k0] * phb[(k+k0)*m_br+j+j0] <<std::endl;
							phc[(i0 + i) * m_ar + j0 + j] += pha[(i0 + i) * m_ar + k + k0] * phb[(k + k0) * m_br + j + j0];
						}
					}
				}
			}
		}
	}

	Time2 = clock();

	sprintf(st, "%4.4f", (double)(Time2 - Time1) / CLOCKS_PER_SEC);
	cout << "Time: " << st << " seconds\n";
	timeLogger << m_ar << ",\t" << m_br << ",\t" << blockSize << ",\t" << st << ",\t";

	cout << "Result matrix: " << endl;
	for (i = 0; i < 1; i++)
	{
		for (j = 0; j < min(10, m_br); j++)
			cout << phc[j] << " ";
	}
	cout << endl;

	free(pha);
	free(phb);
	free(phc);
}

void OnMultLine(int m_ar, int m_br)
{

	SYSTEMTIME Time1, Time2;

	char st[100];
	double temp;
	int i, j, k;

	double *pha, *phb, *phc;

	pha = (double *)malloc((m_ar * m_ar) * sizeof(double));
	phb = (double *)malloc((m_ar * m_ar) * sizeof(double));
	phc = (double *)malloc((m_ar * m_ar) * sizeof(double));

	for (i = 0; i < m_ar; i++)
		for (j = 0; j < m_ar; j++)
			pha[i * m_ar + j] = (double)1.0;

	for (i = 0; i < m_br; i++)
		for (j = 0; j < m_br; j++)
			phb[i * m_br + j] = (double)(i + 1);

	Time1 = clock();

	for (i = 0; i < m_ar; i++)
	{
		for (k = 0; k < m_ar; k++)
		{
			for (j = 0; j < m_br; j++)
			{
				phc[i * m_ar + j] += pha[i * m_ar + k] * phb[k * m_br + j];
			}
		}
	}

	Time2 = clock();

	sprintf(st, "%4.4f", (double)(Time2 - Time1) / CLOCKS_PER_SEC);
	cout << "Time: " << st << " seconds\n";

	timeLogger << m_ar << ", " << m_br << ", " << st << ", ";

	cout << "Result matrix: " << endl;
	for (i = 0; i < 1; i++)
	{
		for (j = 0; j < min(10, m_br); j++)
			cout << phc[j] << " ";
	}
	cout << endl;

	free(pha);
	free(phb);
	free(phc);
}

float produtoInterno(float *v1, float *v2, int col)
{
	int i;
	float soma = 0.0;

	for (i = 0; i < col; i++)
		soma += v1[i] * v2[i];

	return (soma);
}

void handle_error(int retval)
{
	printf("PAPI error %d: %s\n", retval, PAPI_strerror(retval));
	exit(1);
}

void init_papi()
{
	int retval = PAPI_library_init(PAPI_VER_CURRENT);
	if (retval != PAPI_VER_CURRENT && retval < 0)
	{
		printf("PAPI library version mismatch!\n");
		exit(1);
	}
	if (retval < 0)
		handle_error(retval);

	std::cout << "PAPI Version Number: MAJOR: " << PAPI_VERSION_MAJOR(retval)
			  << " MINOR: " << PAPI_VERSION_MINOR(retval)
			  << " REVISION: " << PAPI_VERSION_REVISION(retval) << "\n";
}

int main(int argc, char *argv[])
{

	char c;
	int lin, col, blockSize, nt = 1;
	int op;

	int EventSet = PAPI_NULL;
	long long values[2];
	int ret;

	ret = PAPI_library_init(PAPI_VER_CURRENT);
	if (ret != PAPI_VER_CURRENT)
		std::cout << "FAIL" << endl;

	ret = PAPI_create_eventset(&EventSet);
	if (ret != PAPI_OK)
		cout << "ERRO: create eventset" << endl;

	ret = PAPI_add_event(EventSet, PAPI_L1_DCM);
	if (ret != PAPI_OK)
		cout << "ERRO: PAPI_L1_DCM" << endl;

	ret = PAPI_add_event(EventSet, PAPI_L2_DCM);
	if (ret != PAPI_OK)
		cout << "ERRO: PAPI_L2_DCM" << endl;

	op = atoi(argv[1]);
	lin = atoi(argv[2]);
	col = atoi(argv[3]);
	blockSize = atoi(argv[4]);

	timeLogger.open("timeLogger.txt", std::ios_base::app);

	// Start counting
	ret = PAPI_start(EventSet);
	if (ret != PAPI_OK)
		cout << "ERRO: Start PAPI" << endl;

	switch (op)
	{
	case 1:
		first_implementation(lin, col);
		break;
	case 2:
		OnMultLine(lin, col);
		break;
	case 3:
		cout << "Block impl with blockSize: " << blockSize << endl;
		blockMult(lin, col, blockSize);
		break;
	}

	ret = PAPI_stop(EventSet, values);
	if (ret != PAPI_OK)
		cout << "ERRO: Stop PAPI" << endl;
	printf("L1 DCM: %lld \n", values[0]);
	printf("L2 DCM: %lld \n", values[1]);

	timeLogger << values[0] << ", " << values[1] << endl;

	timeLogger.close();

	ret = PAPI_reset(EventSet);
	if (ret != PAPI_OK)
		std::cout << "FAIL reset" << endl;

	ret = PAPI_remove_event(EventSet, PAPI_L1_DCM);
	if (ret != PAPI_OK)
		std::cout << "FAIL remove event" << endl;

	ret = PAPI_remove_event(EventSet, PAPI_L2_DCM);
	if (ret != PAPI_OK)
		std::cout << "FAIL remove event" << endl;

	ret = PAPI_destroy_eventset(&EventSet);
	if (ret != PAPI_OK)
		std::cout << "FAIL destroy" << endl;
}
