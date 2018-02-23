#include<iostream>
#include<ctime>
#include<fstream>
#include<cmath>
#include<set>
#include<string>
using namespace std;
double* probability(unsigned int ***A, int x, int y, int z, double*q, double a, double b, unsigned int iteration);
double path(unsigned int ***A, int, double*, double, double, unsigned int, double*, set<string>* sequence, int* chain);
void S_count(unsigned int ***A, set<string>* sequence, int N, int*S, unsigned int iteration);
double Gyration_Radius(int* chain, int N);
int d = 3;
const int QQQ = 700;
int main()
{
	srand(time(0));
	unsigned int ***A = (unsigned int ***)malloc(QQQ*sizeof(unsigned int**));
	for (int i = 0; i < QQQ; i++)
	{
		A[i] = (unsigned int**)malloc(QQQ*sizeof(unsigned int*));
		for (int j = 0; j < QQQ; j++)
		{
			A[i][j] = (unsigned int*)malloc(QQQ*sizeof(unsigned int));
			for (int k = 0; k < QQQ; k++)
				A[i][j][k] = 0;
		}
	}
	int N, rep = 900, multiply;
	double answer = 0, ret = 0, repeats = 0, s1 = 0, sb = 0, Gradius = 0;//anwser* - R, ret - checker, repeats-rep,S1;
																		 //cin >> N;
	double*prob = (double*)malloc((2 * d + 1)*sizeof(double));
	double*ans = (double*)malloc(2 * sizeof(double));
	int*S = new int[2];
	set<string>* sequence = new set<string>;
	sequence[0].insert(to_string(QQQ / 2) + "|" + to_string(QQQ / 2) + "|" + to_string(QQQ / 2));
	ofstream f;
	f.open("out.txt");
	f << "b" << " " << "a" << " " << "R" << " " << "V" << " " << "Gradius" << " " << "s1" << " " << "sb" << " " << "s1 + sb" << " " << "N" << endl;
	bool a_c = false;
	bool b_c = true;
	for (double b = 10000; b <= 10000; b = b*((int)b_c * 5 + (int)!b_c * 2))
	{
		b_c=!b_c;
		cout << "b= " << b << endl;
		a_c = true;
		for (double a = 10000; a <= 10000;a=a*((int)a_c * 3) + a*(int)!a_c * 10/3.)
		{
			a_c = !a_c;
			cout << a << endl;
			for (N = 5000000; N <= 10000000; N += 1000000)
			{
				int* chain = new int[N];
				for (int q = 0; q < rep; q++)
				{
					multiply = q + 11000 * a;
					ret = path(A, N, prob, a, b, multiply, ans, sequence, chain);
					S_count(A, sequence, N, S, multiply);
					if (ret == -1)
					{
						cout << "Out of cell on parameters b=" << b << "  a=" << a << "  q=" << q << endl;
						q--;
						continue;
					}
					answer += ans[0] / (1.*rep);
					repeats += ans[1] / (1.*rep);
					Gradius += Gyration_Radius(chain, N) / (1.*rep);
					s1 += S[0] / (1.*rep);
					sb += S[1] / (1.*rep);
					//cout << q << endl;
					sequence[0].clear();
					sequence[0].insert(to_string(QQQ / 2) + "|" + to_string(QQQ / 2) + "|" + to_string(QQQ / 2));
				}
				f << b << " " << a << " " << answer << " " << N*1. - repeats << " " <<Gradius<<" "<< s1 << " " << sb << " " << s1 + sb << " " << N << endl;
				answer = 0;
				repeats = 0;
				Gradius = 0;
				s1 = 0;
				sb = 0;
				for (int i = 0; i < QQQ; i++)
					for (int j = 0; j < QQQ; j++)
						for (int k = 0; k < QQQ; k++)
							A[i][j][k] = 0;
				cout << N << endl;
				delete[]chain;
			}
			f << endl;
		}
	}
}
double path(unsigned int ***A, int N, double*prob, double a, double b, unsigned int iteration, double* answer, set<string>*sequence, int*chain)
{
	double sum, pb_result;
	int x, y, z, count = 0;
	x = QQQ / 2; y = QQQ / 2; z = QQQ / 2;
	A[x][y][z] = iteration;
	for (int i = 0; i < N; i++)
	{
		prob = probability(A, x, y, z, prob, a, b, iteration);
		if (prob == NULL)
			return -1;
		sum = prob[0] * 1. / prob[2 * d];
		pb_result = rand()*1. / RAND_MAX*1.;//выбросили выроятность
		for (int j = 1; j <= 2 * d; j++)// сверяем вероятности и выбираем, куда идти.
		{
			if (pb_result <= sum)
			{
				if (prob[j - 1] == a)
					count++;
				if (j - 1 == 0)
				{
					A[x - 1][y][z] = iteration;
					x -= 1;
					sequence[0].insert(to_string(x) + "|" + to_string(y) + "|" + to_string(z));
					chain[i] = 0;
				}
				else if (j - 1 == 1)
				{
					A[x + 1][y][z] = iteration;
					x += 1;
					sequence[0].insert(to_string(x) + "|" + to_string(y) + "|" + to_string(z));
					chain[i] = 1;
				}
				else if (j - 1 == 2)
				{
					A[x][y - 1][z] = iteration;
					y -= 1;
					sequence[0].insert(to_string(x) + "|" + to_string(y) + "|" + to_string(z));
					chain[i] = 2;
				}
				else if (j - 1 == 3)
				{
					A[x][y + 1][z] = iteration;
					y += 1;
					sequence[0].insert(to_string(x) + "|" + to_string(y) + "|" + to_string(z));
					chain[i] = 3;
				}
				else if (j - 1 == 4)
				{
					A[x][y][z - 1] = iteration;
					z -= 1;
					sequence[0].insert(to_string(x) + "|" + to_string(y) + "|" + to_string(z));
					chain[i] = 4;
				}
				else
				{
					A[x][y][z + 1] = iteration;
					z += 1;
					sequence[0].insert(to_string(x) + "|" + to_string(y) + "|" + to_string(z));
					chain[i] = 5;
				}
				break;
			}
			sum += (1.*prob[j] / prob[2 * d]);
		}
	}
	answer[0] = sqrt((QQQ / 2 - x)*(QQQ / 2 - x) + (QQQ / 2 - y)*(QQQ / 2 - y) + (QQQ / 2 - z)*(QQQ / 2 - z));
	answer[1] = count;
	return 1;
}
double* probability(unsigned int ***A, int x, int y, int z, double*q, double a, double b, unsigned int iteration)
{
	if ((x - 2 == 0) || (x + 2 == QQQ - 1) || (y - 2 == 0) || (y + 2 == QQQ - 1) || (z - 2 == 0) || (z + 2 == QQQ - 1))
		return NULL;
	for (int i = 0; i < 2 * d; i++)
		q[i] = 1.00001;
	if (A[x - 1][y][z] == iteration)
		q[0] = a;
	else if (A[x - 2][y][z] == iteration)
		q[0] = b;

	if (A[x + 1][y][z] == iteration)
		q[1] = a;
	else if (A[x + 2][y][z] == iteration)
		q[1] = b;

	if (A[x][y - 1][z] == iteration)
		q[2] = a;
	else if (A[x][y - 2][z] == iteration)
		q[2] = b;

	if (A[x][y + 1][z] == iteration)
		q[3] = a;
	else if (A[x][y + 2][z] == iteration)
		q[3] = b;

	if (A[x][y][z - 1] == iteration)
		q[4] = a;
	else if (A[x][y][z - 2] == iteration)
		q[4] = b;

	if (A[x][y][z + 1] == iteration)
		q[5] = a;
	else if (A[x][y][z + 2] == iteration)
		q[5] = b;

	if (A[x - 1][y - 1][z] == iteration)
	{
		if (q[0] != a)
			q[0] = b;
		if (q[2] != a)
			q[2] = b;
	}
	if (A[x - 1][y + 1][z] == iteration)
	{
		if (q[0] != a)
			q[0] = b;
		if (q[3] != a)
			q[3] = b;
	}
	if (A[x + 1][y + 1][z] == iteration)
	{
		if (q[1] != a)
			q[1] = b;
		if (q[3] != a)
			q[3] = b;
	}
	if (A[x + 1][y - 1][z] == iteration)
	{
		if (q[1] != a)
			q[1] = b;
		if (q[2] != a)
			q[2] = b;
	}


	if (A[x - 1][y][z - 1] == iteration)
	{
		if (q[0] != a)
			q[0] = b;
		if (q[4] != a)
			q[4] = b;
	}
	if (A[x - 1][y][z + 1] == iteration)
	{
		if (q[0] != a)
			q[0] = b;
		if (q[5] != a)
			q[5] = b;
	}
	if (A[x + 1][y][z + 1] == iteration)
	{
		if (q[1] != a)
			q[1] = b;
		if (q[5] != a)
			q[5] = b;
	}
	if (A[x + 1][y][z - 1] == iteration)
	{
		if (q[1] != a)
			q[1] = b;
		if (q[4] != a)
			q[4] = b;
	}


	if (A[x][y - 1][z - 1] == iteration)
	{
		if (q[2] != a)
			q[2] = b;
		if (q[4] != a)
			q[4] = b;
	}
	if (A[x][y - 1][z + 1] == iteration)
	{
		if (q[2] != a)
			q[2] = b;
		if (q[5] != a)
			q[5] = b;
	}
	if (A[x][y + 1][z + 1] == iteration)
	{
		if (q[3] != a)
			q[3] = b;
		if (q[5] != a)
			q[5] = b;
	}
	if (A[x][y + 1][z - 1] == iteration)
	{
		if (q[3] != a)
			q[3] = b;
		if (q[4] != a)
			q[4] = b;
	}
	double sum = 0;
	for (int i = 0; i < 2 * d; i++)
		sum += q[i];
	q[2 * d] = sum;
	return q;
}
void S_count(unsigned int ***A, set<string>* sequence, int N, int*S, unsigned int iteration)
{
	S[0] = S[1] = 0; //s1 & sb
	int x, y, z;
	int xi, yi, st;
	bool flag = true;
	for (auto &seq : sequence[0])
	{
		xi = yi = -1;
		st = 0;
		while (yi == -1)
		{
			if (seq[st] == '|')
				if (xi == -1)
					xi = st;
				else if (yi == -1)
					yi = st;
			st++;
		}
		x = stoi(seq.substr(0, xi));
		y = stoi(seq.substr(xi + 1, yi - xi - 1));
		z = stoi(seq.substr(yi + 1));
		if (A[x - 1][y][z] != iteration)
			if (A[x - 2][y][z] == iteration || A[x - 1][y - 1][z] == iteration || A[x - 1][y + 1][z] == iteration || A[x - 1][y][z - 1] == iteration || A[x - 1][y][z + 1] == iteration)
				S[1] += 1;
			else
				S[0] += 1;
		if (A[x + 1][y][z] != iteration)
			if (A[x + 2][y][z] == iteration || A[x + 1][y - 1][z] == iteration || A[x + 1][y + 1][z] == iteration || A[x + 1][y][z - 1] == iteration || A[x + 1][y][z + 1] == iteration)
				S[1] += 1;
			else
				S[0] += 1;
		if (A[x][y - 1][z] != iteration)
			if (A[x][y - 2][z] == iteration || A[x - 1][y - 1][z] == iteration || A[x + 1][y - 1][z] == iteration || A[x][y - 1][z - 1] == iteration || A[x][y - 1][z + 1] == iteration)
				S[1] += 1;
			else
				S[0] += 1;
		if (A[x][y + 1][z] != iteration)
			if (A[x][y + 2][z] == iteration || A[x - 1][y + 1][z] == iteration || A[x + 1][y + 1][z] == iteration || A[x][y + 1][z - 1] == iteration || A[x][y + 1][z + 1] == iteration)
				S[1] += 1;
			else
				S[0] += 1;
		if (A[x][y][z - 1] != iteration)
			if (A[x][y][z - 2] == iteration || A[x - 1][y][z - 1] == iteration || A[x + 1][y][z - 1] == iteration || A[x][y - 1][z - 1] == iteration || A[x][y + 1][z - 1] == iteration)
				S[1] += 1;
			else
				S[0] += 1;
		if (A[x][y][z + 1] != iteration)
			if (A[x][y][z + 2] == iteration || A[x - 1][y][z + 1] == iteration || A[x + 1][y][z + 1] == iteration || A[x][y - 1][z + 1] == iteration || A[x][y + 1][z + 1] == iteration)
				S[1] += 1;
			else
				S[0] += 1;
	}
}
double Gyration_Radius(int* chain, int N)
{
	int x = QQQ / 2;
	int y = QQQ / 2;
	int z = QQQ / 2;
	int vec[3];
	vec[0] = QQQ / 2;
	vec[1] = QQQ / 2;
	vec[2] = QQQ / 2;
	for (int i = 0; i < N - 1; i++)
	{
		if (chain[i] == 0)
			x--;
		else if (chain[i] == 1)
			x++;
		else if (chain[i] == 2)
			y--;
		else if (chain[i] == 3)
			y++;
		else if (chain[i] == 4)
			z--;
		else
			z++;
		vec[0] += x;
		vec[1] += y;
		vec[2] += z;
	}
	double vec0[3];
	vec0[0] = vec[0] / N*1.;
	vec0[1] = vec[1] / N*1.;
	vec0[2] = vec[2] / N*1.;

	x = QQQ / 2;
	y = QQQ / 2;
	z = QQQ / 2;
	double S2 = 0;
	S2 = (x - vec0[0])*(x - vec0[0]) + (y - vec0[1])*(y - vec0[1]) + (z - vec0[2])*(z - vec0[2]);
	for (int i = 0; i < N - 1; i++)
	{
		if (chain[i] == 0)
			x--;
		else if (chain[i] == 1)
			x++;
		else if (chain[i] == 2)
			y--;
		else if (chain[i] == 3)
			y++;
		else if (chain[i] == 4)
			z--;
		else
			z++;
		S2 += (x - vec0[0])*(x - vec0[0]) + (y - vec0[1])*(y - vec0[1]) + (z - vec0[2])*(z - vec0[2]);
	}
	return S2*1. / N;
}