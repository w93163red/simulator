// simuate_2.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"
#include <cmath>
#include <iostream>
#include <fstream>
#include <cstdlib>
#include <set>
#include <vector>
#include <algorithm>
#include <map>
#include <ctime>
#include <random>
#include <functional>
#include <assert.h>

using namespace std;
int T;
double BRT = 0.008;
int M;
double u;
double max_rate = 0;
int cache_size = 64;

class Task
{
public:
	int period;
	double c;
	set<int> UCB;
	set<int> ECB;
	Task(int p, double cc) : period(p), c(cc) {}
};

bool comp(Task a, Task b)
{
	return a.period > b.period;
}

long long gcd(long long a, long long b)
{
	return b == 0 ? a : gcd(b, a % b);
}

long long lcm(vector<int> num)
{
	long long ans = num[0] * num[1] / gcd(num[0], num[1]);
	for (int i = 0; i < num.size(); i++)
		ans = ans * num[i] / gcd(ans, num[i]);
	return ans;
}

vector<vector<double>> UUniFastDiscard(int n, double u, int nsets)
{
	vector<vector<double>> sets;
	double nextSumU;
	srand(time(NULL));
	while (sets.size() < nsets)
	{
		vector<double> utilizations;
		double sumU = u;
		bool flag = true;
		for (int i = 1; i < n; i++)
		{
			nextSumU = sumU * pow((rand() / (RAND_MAX*1.0)), (1.0 / (n - i)));
			utilizations.push_back(sumU - nextSumU);
			sumU = nextSumU;
		}
		utilizations.push_back(nextSumU);
		for (double ut : utilizations)
		{
			if (ut > 1)
			{
				flag = false;
				break;
			}
		}
		if (flag)
			sets.push_back(utilizations);
	}

	return sets;
}

vector<vector<int>> periods_generator(int n, int nsets)
{
	vector<vector<int>> ans;

	default_random_engine generator;
	uniform_int_distribution<int> distribution(1, 4);

	for (int i = 0; i < nsets; i++)
	{
		vector<int> v;
		for (int j = 0; j < n; j++)
		{
			int period = distribution(generator);
			v.push_back(2 << period);
		}
		ans.push_back(v);
	}
	return ans;
}

double cal_DBF(vector<Task> task_list, int T)
{
	double sum = 0;
	for (Task t : task_list)
	{
		if (t.period <= T)
			sum += max((double)0.0, (floor((T - t.period) / t.period) + 1) * t.c);
	}
	return sum;
}

double cal_DBF_Ju(vector<Task> task_list, int T)
{
	double sum = 0;
	for (int i = 0; i < task_list.size(); i++)
	{
		if (task_list[i].period >= T)
			continue;
		double gamma = 0;
		for (int j = i + 1; j < task_list.size(); j++)
		{
			if (task_list[i].period > task_list[j].period)
			{
				vector<int> v(128);
				vector<int>::iterator it;
				it = set_intersection(task_list[i].UCB.begin(), task_list[i].UCB.end(),
					task_list[j].ECB.begin(), task_list[j].ECB.end(), v.begin());
				v.resize(it - v.begin());
				gamma += max((double)0.0, ceil((double)(task_list[i].period - task_list[j].period) / (double)task_list[j].period)) 
					* v.size();
			}
		}
		gamma *= BRT;
		sum += max((double)0.0, floor((double)(T - task_list[i].period) / (double)task_list[i].period)) * gamma;
	}
	return sum;
}

double cal_DBF_combined(vector<Task> task_list, double &ECB, double &UCB, int T)
{
	double sum = 0;

	//ECB, UCB
	for (int i = task_list.size() - 1; i >= 0; i--)
	{
		if (task_list[i].period >= T)
			continue;
		double gamma = 0;
		long times;
		multiset<int> res;
		multiset<int> UCB_res;
		vector<int> v(128);
		vector<int>::iterator it;
		for (int k = 0; k < i; k++)
		{
			times = max(0.0, ceil((task_list[k].period - task_list[i].period) / task_list[i].period)) *
				max(0.0, 1 + floor((T - task_list[i].period) / task_list[i].period));
			set<int> med;
			for (int h = i; h < task_list.size(); h++)
			{
				v = vector<int>(128);
				it = set_union(med.begin(), med.end(), task_list[h].ECB.begin(), task_list[h].ECB.end(), v.begin());
				v.resize(it - v.begin());
				med.clear();
				med = set<int>(v.begin(), v.end());
			}
			v = vector<int>(128);
			it = set_intersection(task_list[k].UCB.begin(), task_list[k].UCB.end(), med.begin(), med.end(), v.begin());
			v.resize(it - v.begin());
			for(int n: med)
				for (int ii = 0; ii < times; ii++)
				{
					res.insert(n);
				}
			for (int n: task_list[k].UCB)
				for (int ii = 0; ii < times; ii++)
				{
					UCB_res.insert(n);
				}
		}
		int L = max(0.0, 1 + floor((T - task_list[i].period) / task_list[i].period));
		ECB += BRT * L * res.size();

		multiset<int> ECB_res;
		for (int n : task_list[i].ECB)
			for (int ii = 0; ii < L; ii++)
			{
				ECB_res.insert(n);
			}
		v = vector<int>(128);
		it = set_intersection(UCB_res.begin(), UCB_res.end(), ECB_res.begin(), ECB_res.end(), v.begin());
		v.resize(it - v.begin());
		UCB += BRT * v.size();

		sum += min(ECB, UCB);
	}
	return sum;
}

double cal_DBF_condensed(vector<Task> task_list, int T)
{
	return 0;
}

int main()
{
	ofstream myfile;
	myfile.open("result.txt");
	int n = 12;
	int nsets = 100;
	random_device rd;
	mt19937 rng(rd());
	uniform_int_distribution<int> uni(0, 10);
	vector<int> kernel = { 2 };
	//Generate Taskset	
	for (int M : kernel)
	{
		for (u = 0.5; u < M*1.0; u += 0.1*M / 2)
		{
			BRT = 0.008;
			vector<vector<double>> utils = UUniFastDiscard(n, u, nsets);
			vector<vector<int>> periods = periods_generator(n, nsets);
			vector<int> ans = { 0,0,0,0,0,0 };
			vector<Task> taskset;
			myfile << "M: " << M << " u: " << u << endl;
			cout << "M: " << M << " u: " << u << endl;
			for (int i = 0; i < periods.size(); i++)
			{
				bool DBF_pass = true, DBF_Ju_pass = true, ECB_pass = true;
				bool ucb_pass = true, DBF_combined_pass = true, DBF_condensed_pass = true;
				vector<Task> taskset;
				vector<int> task_time;
				for (int j = 0; j < periods[i].size(); j++)
				{
					Task t(periods[i][j], (double)utils[i][j] * (double)periods[i][j]);
					for (int k = j * 8; k < (j + 1) * 8; k++)
					{
						t.UCB.insert(k % cache_size);
					}
					t.ECB.insert(i * 8 % cache_size);
					t.ECB.insert((i * 8 + 1) % cache_size);
					taskset.push_back(t);
					task_time.push_back(t.period);
				}
				sort(taskset.begin(), taskset.end(), comp);
				int hp = taskset[0].period * 2;
				for (int T = hp; T >= 1; T--)
				{
					double ecb = 0, ucb = 0;
					double DBF = cal_DBF(taskset, T);
					double DBF_Ju;
					if (DBF_Ju_pass)
						DBF_Ju = cal_DBF_Ju(taskset, T);
					double DBF_combined;
					if (DBF_combined_pass)
						DBF_combined = cal_DBF_combined(taskset, ecb, ucb, T);
					double DBF_condensed;
					if (DBF_condensed_pass)
						DBF_condensed = cal_DBF_condensed(taskset, T);

					//pre-DBF test
					if (cal_DBF(taskset, T) > M*T)
					{
						DBF_pass = false;
						//cout << "DBF: " << DBF << " " << T << endl;
					}
					if (DBF_Ju_pass && DBF_Ju + DBF > M*T)
					{
						DBF_Ju_pass = false;
						//cout << "DBF_Ju: " << DBF_Ju << " " << T << endl;
					};
					if (ECB_pass > M*T)
					{
						ECB_pass = false;
						//cout << "ECB: " << ecb << " " << T << endl;
					}
					if (ucb_pass > M*T) {
						ucb_pass = false;
						//cout << "UCB: " << ucb << " " << T << endl;
					}
					if (DBF_combined_pass && DBF_combined > M*T)
					{
						DBF_combined_pass = false;
						//cout << "COMBINED: " << DBF_combined << " " << T << endl;
					}
					if (DBF_condensed_pass && DBF_condensed > M*T)
					{
						DBF_condensed_pass = false;
						//cout << "DBF_Condensed: " << DBF_condensed << " " << T << endl;
					}
				}
				if (DBF_pass) ans[0]++;
				if (DBF_Ju_pass) ans[1]++;
				if (ECB_pass) ans[2]++;
				if (ucb_pass) ans[3]++;
				if (DBF_combined_pass) ans[4]++;
				if (DBF_condensed_pass) ans[5]++;
				
			}
			for (int n : ans)
			{
				myfile << n << ",";
			}
			myfile << endl;
		}
	}
    return 0;
}

