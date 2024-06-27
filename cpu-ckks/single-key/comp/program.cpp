#include "program.h"

void upgrade_oddbaby(long n, Tree& tree)	// n should be odd
{
	long d = ceil_to_int(log(static_cast<double>(n))/log(2.0));		// required minimum depth
	long m, l;
	long min, total_min = 10000, min_m = 0, min_l = 0;
	Tree min_tree, total_min_tree;

	for(l=1; pow2(l)-1<=n; l++)
	{
		for(m=1; pow2(m-1)<n; m++)
		{
			// initialization
			vector<vector<int>> f(n+1, vector<int>(d+1,0));
			vector<vector<Tree>> G(n+1, vector<Tree>(d+1, Tree(evaltype::oddbaby)));
			f[1][1] = 0;
			for(int i=3; i<=n; i+=2) f[i][1] = 10000;

			// recursion
			for(int j=2; j<=d; j++) 
			{
				for(int i=1; i<=n; i+=2) 
				{
					if(i <= pow2(l)-1 && i <= pow2(j-1)) f[i][j] = 0;
					else
					{
						min = 10000;
						min_tree.clear();
						for(int k=1; k<=m-1 && pow2(k)<i && k<j; k++)	// g = 2^k
						{
							long g = pow2(k);
							if(f[i-g][j-1] + f[g-1][j] +1 < min) 
							{
								min = f[i-g][j-1] + f[g-1][j] +1;
								min_tree.merge(G[g-1][j],G[i-g][j-1],g);
							}
						}
						f[i][j] = min;
						G[i][j] = min_tree;
					}
				}
			}
			if(f[n][d] + pow2(l-1) + m - 2 < total_min) {
				total_min = f[n][d] + pow2(l-1) + m - 2;
				total_min_tree = G[n][d];
				min_m = m;
				min_l = l;
			}
		}
	}

	// cout << "deg " << n << ": " << total_min << endl;
	// cout << "m: " << min_m << ", l: " << min_l << endl;
	tree = total_min_tree;
	tree.m = min_m;
	tree.l = min_l;

	
}
void upgrade_baby(long n, Tree& tree)
{
	long d = ceil_to_int(log(static_cast<double>(n+1))/log(2.0));		// required minimum depth
	long m, b;
	long min, total_min = 10000, min_m = 0, min_b = 0;
	Tree min_tree, total_min_tree;
	evaltype type = evaltype::baby;

	// cout << "minimum depth: " << d << endl;

	// n==1
	if(n==1)
	{
		total_min = 0;
		total_min_tree = Tree(type);
		min_m = 1;
		min_b = 1;
	}
	
	for(b=1; b<=n; b++)
	{
		for(m=1; pow2(m-1)*b<=n; m++)
		{
		//	cout << "Stage b,m: " << b << " " << m << endl;

			// initialization
			vector<vector<int>> f(n+1, vector<int>(d+1,0));
			vector<vector<Tree>> G(n+1, vector<Tree>(d+1, Tree(type)));

			// recursion
			for(int j=1; j<=d; j++) 
			{
				for(int i=1; i<=n; i++) 
				{
					int k;
					if(i+1 > pow2(j))
					{
						f[i][j] = 10000;
						G[i][j] = Tree(type);
					}
					else if(b==1 && m>=2 && i<=2 && i<=pow2(j-1))
					{
						f[i][j] = 0;
						G[i][j] = Tree(type);
					}
					else if(i<=b && i<=pow2(j-1)) 
					{
						f[i][j] = 0;
						G[i][j] = Tree(type);
					}
					else
					{
						min = 10000;
						min_tree.clear();
						for(int k=2; k<=b; k++) 	// g = k
						{
							long g = k;
							if(g<=pow2(j-1) && 2<=g && g<i && f[i-g][j-1] + f[g-1][j] +1 < min)
							{
								min = f[i-g][j-1] + f[g-1][j] +1;
								min_tree.merge(G[g-1][j],G[i-g][j-1],g);
							} 
						}
						for(int k=0; k<=m-1; k++)	// g = 2^k b
						{
							long g = pow2(k)*b;
							if(g<=pow2(j-1) && 2<=g && g<i && f[i-g][j-1] + f[g-1][j] +1 < min)
							{
								min = f[i-g][j-1] + f[g-1][j] +1;
								min_tree.merge(G[g-1][j],G[i-g][j-1],g);
							} 
						}
						f[i][j] = min;
						G[i][j] = min_tree;
						if(min == 10000) 
						{
					//		cout << "no g found " << b << " " << m << " " << j << " " << i << endl;
					//		throw std::runtime_error("this case should not occur!");
						}
					}
				}
			}
			if(f[n][d] + m + b - 2 < total_min) {
				total_min = f[n][d] + m + b - 2;
				total_min_tree = G[n][d];
				min_m = m;
				min_b = b;
			}
		}
	}

	// cout << "deg " << n << ": " << total_min << endl;
	// cout << "m: " << min_m << ", b: " << min_b << endl;
	tree = total_min_tree;
	tree.m = min_m;
	tree.b = min_b;		

}