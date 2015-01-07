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
#define FITNESS(k)      (-(k)*(k)+20*k)

#define SIZE			30			// population
#define MAXGEN			150			//
#define SEED			1124		// initial GNG seed 
#define P_CORSS			0.70		// 
//#define P_MUTATION		0.10		//
#define P_MU_LOW		0.20
#define P_MU_HIGH		0.07

#define LEN				8			// 数据位数,也就是基因位数
#define randd()			((double)rand() /RAND_MAX)
#define randi(k)		((int)((double)rand() *(k)/RAND_MAX +0.5))

#define MAX(a, b)		((a)>(b)?(a):(b))
#define MIN(a, b)		((a)<(b)?(a):(b))

/******************************************************************
 * 种群个体的结构体 
 ******************************************************************/
typedef struct node{
	//char    x[LEN];				//
	unsigned int	value;			//
	double			fitness;		//
	double			fitsum;			// 
}node;

node cur[SIZE], next[SIZE], max, min;

int gen_cnt;				// 迭代次数

/****************************************************************************
 * Calculate fitness of individual
 ****************************************************************************/
void cal_fitness(struct node *pop)
{
	int i, j;
	int k;
	long d;

	k = (int)(pop[0].value) - 100;
	pop[0].fitness = FITNESS(k);
	pop[0].fitsum = pop[0].fitness;
	
	for ( i=1; i<SIZE; i++ ) {
		k = (int)(pop[i].value) - 100;
		pop[i].fitness = FITNESS(k);
		pop[i].fitsum = pop[i].fitness + pop[i-1].fitsum;
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
	}
	cal_fitness(cur);
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
	return SIZE-1;
}

/**********************************************************************
 * 迭代
 **********************************************************************/
void tran()
{
	int i, j;
	unsigned int pos;
	int k;
	int n_max = 0;

	/* 找当前种群最优个体 */
	max = cur[0];
	for ( i=1; i<SIZE-1; i++ ) {
		if ( cur[i].fitness>max.fitness ) {
			max = cur[i];
			n_max = i;
		}
	}

	for ( k=0; k<SIZE; k+=2 ) {
		int r;
		/* Sellect the individual for cross*/
		i = sel();
		j = sel();
		/* Cross */ 
		if ( randd()<P_CORSS ) {
			int tmp_0, tmp_1;

			/* 选择交叉位置 */ 
			pos = randi(LEN-1);
			//printf("P_CORSS[%2u]:(%u, %u)", pos, cur[i].value, cur[j].value);
			r = (1 <<(pos+1)) -1;
			pos = ((1<<LEN) - 1)^r;
			//printf("-%d,%d-,", r, pos);
			tmp_0 = (cur[i].value & r) | (cur[j].value & pos);
			tmp_1 = (cur[j].value & r) | (cur[i].value & pos);
			next[k].value = tmp_0>N?
					cur[i].value:tmp_0;
			next[k+1].value = tmp_1>N?
					cur[j].value:tmp_1;
			//printf("[%u, %u]->%u, %u\n", tmp_0, tmp_1, next[k].value, next[k+1].value);
		} else {
			next[k].value = cur[i].value;
			next[k+1].value = cur[j].value;
		}
		// 变异
		if ( randd()<P_MU_LOW ) {
			int tmp;
			pos = 1<<(LEN/2)-1;			
			pos = randi(pos);
			//printf(">>P_MU_LOW:(%3d)%u, %u", pos, next[k].value, next[k+1].value);
			tmp = next[k].value;
			//pos = tmp & pos;
			tmp ^= pos;
			//printf("\np%d,t%d,v%d", pos, tmp, next[k].value);
			next[k].value = tmp>N?
					next[k].value:tmp;
			//printf(",%d", next[k].value);

			pos = 1<<(LEN/2)-1;			
			pos = randi(pos);
			tmp = next[k+1].value;
			//pos = tmp & pos;
			tmp ^= pos;
			next[k+1].value = tmp>N?
					cur[k+1].value:tmp;
			//printf("->%u, %u\n", next[k].value, next[k+1].value);
		}
#if 0
		if ( randd()<P_MU_HIGH ) {
			pos = randi(LEN-LEN/2) + (LEN/2);
			printf(">>P_MU_HIGH:(%3d)%u, %u->", pos, next[k].value, next[k+1].value);
			r = 1 << pos;
			pos = next[k+1].value & r;
			next[k].value ^= pos;
			next[k].value = cur[k].value>N?
					pos:cur[k].value;

			pos = randi(LEN-LEN/2) + (LEN/2);
			r = 1 << pos;
			pos = next[k+1].value & r;
			next[k+1].value ^= pos;
			next[k+1].value = cur[k+1].value>N?
					pos:cur[k+1].value;
			printf("%u, %u\n", next[k].value, next[k+1].value);
		}
#else
		if ( randd()<P_MU_HIGH ) {
			pos = randi(LEN-1);
			//printf(">>P_MU_HIGH:(%3d)%u, %u->", pos, next[k].value, next[k+1].value);
			r = 1 << pos;
			pos = next[k+1].value & r;
			next[k].value ^= pos;
			next[k].value = cur[k].value>N?
					pos:cur[k].value;

			pos = randi(LEN-1);
			r = 1 << pos;
			pos = next[k+1].value & r;
			next[k+1].value ^= pos;
			next[k+1].value = cur[k+1].value>N?
					pos:cur[k+1].value;
			//printf("%u, %u\n", next[k].value, next[k+1].value);
		}
#endif
	}

	cal_fitness(next);		/* calculate the fitness of next generation */

	for (i=0; i<SIZE; i++) {
		/* copy the best indivadual to the next generation */
		if (cur[i].fitness == max.fitness && 
			next[i].fitness < max.fitness) {
			next[i] = max;
		}
	}
	/* search the worst individual of the new generation and replace it*/
	min = cur[0];
	j = 0;
	for ( i=1; i<SIZE; i++ ) {
		if (next[i].fitness < min.fitness) {
			min.fitness = cur[i].fitness;
			j = i;
		}
	}
	next[j] = max;
	if (next[0].fitness < max.fitness) {
		next[0] = max;
	}

	memcpy(cur, next, sizeof(cur));
}

// 打印种群
void printcur()
{
    int i;
    for ( i=0; i<SIZE; i++ ) { 
        //print(cur[i]);
        printf("NO.%2d: Value = %8u fitness = %6.3f\n", 
                i, cur[i].value, cur[i].fitness);
    }
}


void GA()
{
	int 	cnt=0;
	long	k;
	double	ans;
	long 	val;
	int 	i, j;

	while ( cnt++<MAXGEN ) {
		printf("Generation:%4d", cnt);
		tran();
#if 1
		printf("\nmax = %6.3f\n", max.fitness);
		printcur();
		printf("=======OK=======\n\n");
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
	printf("k=%ld, anwer:%6.3f\n", k, ans);
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
