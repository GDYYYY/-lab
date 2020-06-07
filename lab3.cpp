#include "common.h"
#include <algorithm>

using namespace std;

//You should only code here.Don't edit any other files in this
int func1(int amount, vector<int> &coins)
{
	int *dp = new int[amount + 1];
	int num = coins.size();
	// sort(coins.begin(),coins.end());
	dp[0] = 1;
	for (int i = 1; i < amount + 1; i++)
		dp[i] = 0;
	for (int i = 0; i < num; i++)
	{
		for (int j = 0; j < amount + 1; j++)
		{
			if (coins[i] <= j)
				dp[j] = dp[j] + dp[j - coins[i]];
		}
	}
	int ans = dp[amount];
	delete[] dp;
	return ans;
}

int func2(int amount, vector<vector<int>> &conquer)
{
	int **dp = new int *[amount];
	int ans = 0;
	for (int i = 0; i < amount; i++)
		dp[i] = new int[amount];
	for (int i = 0; i < amount; i++)
		for (int j = 0; j < amount; j++)
			dp[i][j] = 0;

	for (int j = 0; j < amount; j++)
	{
		dp[j][(j + 1) % amount] = conquer[j][(j + 1) % amount] ? 1 : -1;
	}
	for (int i = 2; i <= amount; i++)
	{
		for (int j = 0; j < amount; j++)
		{
			int tmp = 0;
			for (int t = 1; t < i; t++)
			{
				//两者均能与中介决斗并至少有一个能赢中介
				if ((dp[j][(j + t) % amount] == 1 || dp[(j + t) % amount][(j + i) % amount] == -1) && dp[j][(j + t) % amount] && dp[(j + t) % amount][(j + i) % amount])
				{
					tmp = conquer[j][(j + i) % amount] ? 1 : -1;
					break;
				}
			}
			dp[j][(j + i) % amount] = tmp;
		}
	}

	for (int j = 0; j < amount; j++)
	{
		if ((dp[j][j]) != 0)//可以同自己决斗（即可以赢其他所有人）
			ans++;
	}

	for (int i = 0; i < amount; i++)
		delete[] dp[i];
	delete[] dp;
	return ans;
}
//获得高斯消元后的矩阵
void precompute(double **matrix, int n, double **factor)
{
	for (int k = 0; k < n; k++)
	{
		for (int i = k + 1; i < n; i++)
		{
			factor[i][k] = matrix[i][k] / matrix[k][k]; //初等行变换的系数
			for (int j = 0; j < n; j++)
				matrix[i][j] = matrix[i][j] - factor[i][k] * matrix[k][j];
		}
	}
}
//计算求解
void compute(double **matrix, double **dp, int n, int hp, double **factor)
{
	for (int k = 0; k < n; k++)
	{
		for (int i = k + 1; i < n; i++)
		{
			matrix[i][n] = matrix[i][n] - factor[i][k] * matrix[k][n];
		}
	}
	dp[hp][n - 1] = matrix[n - 1][n] / matrix[n - 1][n - 1];
	for (int i = n - 2; i >= 0; i--)
	{
		double sum = 0;
		for (int j = i + 1; j < n; j++)
			sum += matrix[i][j] * dp[hp][j];
		dp[hp][i] = (matrix[i][n] - sum) / matrix[i][i];
	}
}

double func3(int n, int hp, vector<int> &damage, vector<int> &edges)
{
	//声明变量
	double ans = 0;
	int *numOfEdge = new int[n];
	int **graph = new int *[n];
	double **matrix = new double *[n];
	double **factor = new double *[n];
	for (int i = 0; i < n; i++)
	{
		numOfEdge[i] = 0;
		graph[i] = new int[n];
		matrix[i] = new double[n + 1];
		factor[i] = new double[n];
	}
	for (int i = 0; i < n; i++)
		for (int j = 0; j < n; j++)
			graph[i][j] = 0;
	double **dp = new double *[hp + 1];
	for (int i = 0; i <= hp; i++)
		dp[i] = new double[n];
	for (int i = 0; i <= hp; i++)
		for (int j = 0; j < n; j++)
			dp[i][j] = 0;

	int x, y, size = edges.size();
	double c;
	for (int i = 0; i < size; i += 2)
	{
		x = edges[i] - 1;
		y = edges[i + 1] - 1;
		if (x != n - 1)
			graph[x][y] = 1;
		if (y != n - 1)
			graph[y][x] = 1;
		//到达终点不会再去其他点,因此为单向边
		numOfEdge[x]++;
		numOfEdge[y]++;
	} //存图

	for (int i = hp; i >= 0; i--)
	{
		if (i == hp)
		{
			for (int j = 0; j < n; j++)
			{
				c = 0;
				for (int k = 0; k < n; k++)
				{
					matrix[j][k] = 0;
					if (graph[k][j] != 0 && j != k) //有从k通向j的路径
					{
						if (damage[j] == 0)
							matrix[j][k] = 1 / (double)numOfEdge[k];
						else if (i + damage[j] <= hp)
							//加上伤害值后血量不会超过满血（若超过说明不可能在该血量到达，即可能性应该为0）
							c += dp[i + damage[j]][k] / (double)numOfEdge[k];
					}
					else if (j == k)
						matrix[j][k] = -1;
				}
				matrix[j][n] = -c; //常数项
			}
			matrix[0][n] -= 1; //起始位置期望概率多1
			precompute(matrix, n, factor);
		}
		else
		{
			for (int j = 0; j < n; j++)
			{
				c = 0;
				for (int k = 0; k < n; k++)
				{
					if (graph[k][j] != 0 && j != k && damage[j] > 0 && i + damage[j] <= hp)
						//有从k通向j的路径
						//且加上伤害值后血量不会超过满血（若超过说明不可能在该血量到达，即可能性应该为0）
						c += dp[i + damage[j]][k] / (double)numOfEdge[k];
				}
				matrix[j][n] = -c; //常数项
			}
		}
		compute(matrix, dp, n, i, factor);
	}

	for (int i = hp; i > 0; i--)
		ans += dp[i][n - 1];

	for (int i = 0; i < hp; i++)
		delete[] dp[i];
	delete[] dp;
	for (int i = 0; i < n; i++)
	{
		delete[] graph[i];
		delete[] matrix[i];
		delete[] factor[i];
	}
	delete[] graph;
	delete[] matrix;
	delete[] factor;
	delete[] numOfEdge;
	return ans;
}
