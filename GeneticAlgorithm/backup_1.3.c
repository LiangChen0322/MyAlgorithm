/************************************************************************
 * 
 * pow(2,10)<1024*2<pow(2,11)
 * 编码的二进制长度为22 
 ************************************************************************/
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <time.h>
#include "math.h"

#define N				200			/* 分解取值范围，精确到小数点后4位 */
#define FITNESS(k)      (-k*k)

#define SIZE			30			// population
#define MAXGEN			150			//
#define SEED			732		// initial GNG seed 
#define P_CORSS			0.75		// 
//#define P_MUTATION		0.10		//
#define P_MU_LOW		0.30
#define P_MU_HIGH		0.07

#define LEN				8			// 数据位数,也就是基因位数
#define randd()			((double)rand() /RAND_MAX)
#define randi(k)		((int)(randd() *(k) +0.5))

#define MAX(a, b)		((a)>(b)?(a):(b))
#define MIN(a, b)		((a)<(b)?(a):(b))

/******************************************************************
 * 种群个体的结构体 
 ******************************************************************/
typedef struct node{
	//char    x[LEN];				//
	unsigned int	value;			//
	long			fitness;		//
	long			fitsum;			// 
}node;

node cur[SIZE], next[SIZE], max, min;

int gen_cnt;				// 迭代次数

/****************************************************************************
 * Calculate fitness of 
 ****************************************************************************/
void cal_fitness()
{
	int i, j;
	int k;
	long d;

	for ( i=0,d=0; i<SIZE; i++ ) {
		k = (int)(cur[i].value) - 100;
		//k -= 100;
		cur[i].fitness = FITNESS(k);
		cur[i].fitsum = i>0?(cur[i].fitness+cur[i-1].fitsum):(cur[0].fitness);
	}
}

/*****************************************************************************
 * Population initial
 *****************************************************************************/
void population_init()
{
	int k;
	int i;
	int j;

	for ( i=0; i<SIZE; i++ ) {
		cur[i].value = randi(N);
		k = (int)(cur[i].value) - 100;
		cur[i].fitness = FITNESS(k);
		cur[i].fitsum = i>0?(cur[i].fitness+cur[i-1].fitsum):(cur[0].fitness);
	}
}

int sel()
{
	double p = randd();
	double sum = cur[SIZE-1].fitsum;
	int i;

	for ( i=0; i<SIZE; i++ ) {
		if ( cur[i].fitsum/sum>p ) {
			return i;
		} 
	}
}

/**********************************************************************
 * 迭代
 **********************************************************************/
void tran()
{
	int i, j;
	unsigned int pos;
	int k;
	int n_max;
	//找当前种群最优个体 
	max = cur[0];
	for ( i=1; i<SIZE-1; i++ ) {
		if ( cur[i].fitness>max.fitness ) {
			max = cur[i];
			n_max = i;
		}
	}
	for ( k=0; k<SIZE; k+=2 ) {
		int r;
		// 选择交叉个体 
		i = sel();
		j = sel();

		// 选择交叉位置 
		pos = randi(LEN-1);

		// 交叉
		if ( randd()<P_CORSS ) {
		printf("P_CORSS:%u, %u->", next[k].value, next[k+1].value);
		r = 1 <<(pos+1) -1;
		next[k].value = (cur[i].value &r)
				+ (cur[j].value >>pos) <<pos;
		next[k+1].value = cur[j].value &r
				+ (cur[i].value >>pos) <<pos;
		next[k].value = cur[k].value>N?
				cur[i].value:cur[k].value;
		next[k+1].value = cur[k+1].value>N?
				cur[j].value:cur[j].value;
		printf("%u, %u\n", next[k].value, next[k+1].value);
		} else {
			next[k].value = cur[i].value;
			next[k+1].value = cur[j].value;
		}
		// 变异
		if ( randd()<P_MU_LOW ) {
			int tmp;
			pos = randi(LEN/3+0.5-1);
			printf(">>P_MU_LOW:(%3d)%u, %u", pos, next[k].value, next[k+1].value);
			pos = 1 << pos;
			tmp = next[k].value;
			pos = next[k].value & pos;
			next[k].value ^= pos;
			printf("\np%d,t%d,v%d", pos, tmp, next[k].value);
			next[k].value = cur[k].value>N?
					tmp:next[k].value;
			printf(",%d", next[k].value);

			pos = randi(LEN/3+0.5-1);
			pos = 1 << pos;
			tmp = next[k+1].value;
			pos = next[k+1].value & pos;
			next[k+1].value ^= pos;
			next[k+1].value = cur[k+1].value>N?
					tmp:cur[k+1].value;
			printf("->%u, %u\n", next[k].value, next[k+1].value);
		}
		if ( randd()<P_MU_HIGH ) {
			pos = randi(LEN*2/3+0.5-1) + (LEN/3+0.5);
			printf(">>P_MU_HIGH:(%3d)%u, %u->", pos, next[k].value, next[k+1].value);
			r = 1 << pos;
			pos = next[k+1].value & r;
			next[k].value ^= pos;
			next[k].value = cur[k].value>N?
					pos:cur[k].value;

			pos = randi(LEN*2/3+0.5-1) + (LEN/3+0.5);
			r = 1 << pos;
			pos = next[k+1].value & r;
			next[k+1].value ^= pos;
			next[k+1].value = cur[k+1].value>N?
					pos:cur[k+1].value;
			printf("%u, %u\n", next[k].value, next[k+1].value);
		}

		r = (int)(next[k].value) - 100;
		next[k].fitness = FITNESS(r);
		r = (int)(next[k+1].value) - 100;
		next[k+1].fitness = FITNESS(r);
	}
	//找下一代的最差个体
	min = next[0];
	j = 0;
    
	for ( i=1; i<SIZE-1; i++ ) {
		if ( next[i].fitness<min.fitness ) {
			min = next[i];
			j = i;
		}
	}
	//用上一代的最优个体替换下一代的最差个体
	next[0] = max;
	next[j] = max;
	//next[n_max] = max;
	//printf("ok\n");

	memcpy(cur, next, sizeof(cur));

	cal_fitness();
}

// 打印种群
void printcur()
{
    int i;
    for ( i=0; i<SIZE; i++ ) { 
        //print(cur[i]);
        printf("NO.%2d: Value = %8u fitness = %8ld\n", 
                i, cur[i].value, cur[i].fitness);
    }
}


void GA()
{
	int     cnt=0;
	long    k;
	long    ans;
	long    val;
	int     i, j;
    
    while ( cnt++<MAXGEN ) {
        tran();
#if 1
        printf("\nmax = %ld\n",max.fitness);
        printcur();
#endif
    }
    ans = cur[0].fitness;
    cnt = 0;
#if 0
    for( i=1; i<SIZE; i++ ) {
        ans = MIN(ans, cur[i].fitness);
    }
    printf("answer = %ld\n", ans);
#else
    for( i=1; i<SIZE; i++ ) {
        if ( cur[i].fitness>ans ) {
            ans = cur[i].fitness;
            cnt = i;
        }
    }
    k = cur[cnt].value - 100;
    printf("k=%ld, anwer:%ld\n", k, ans);
#endif
}

int main()
{
	//srand((unsigned)time(NULL));
	srand(SEED);

	population_init();
	//cur[3].value = 100;
	//cur[3].fitness = 0;
	GA();

	system("pause");
	return 0;
}
