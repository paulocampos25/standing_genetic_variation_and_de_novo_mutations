#include<math.h>
#include<iostream>
#include<string.h>
#include<stdio.h>
#include<stdlib.h> 
#include<iomanip> 				
#include<fstream>
#include<tuple>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

#include <unistd.h>

using namespace std;

#define PI 3.1415927

#define NIL (0)    
#define IM1 2147483563
#define IM2 2147483399
#define AM (1.0/IM1)
#define IMM1 (IM1-1)
#define IA1 40014
#define IA2 40692
#define IQ1 53668
#define IQ2 52774
#define IR1 12211
#define IR2 3791
#define NTAB 32
#define NDIV (1+IMM1/NTAB)
#define EPS 1.2e-7
#define RNMX (1.0-EPS)

long idum;
gsl_rng *gerador;


struct populacao{
  int N;
  int N0;
  double mut;
  int **step;
  int *nstep;
  int *contind;
  int ntraits;
  double sigma_r;
  double sigma_trait;
  double sigma_trait_mut;
  int n_subpopul;
  double alpha_r;
  int K;
  double migration;
  double lambda;
  int g;
  int diversity;
  double r;
};


struct otimo{
  
  double *traits;
  int tau;
  double *v;
  double *trait_shift;
  
};



struct sequencia{
  int *cont_mut;
  int **ind_mut;
  int **ind_mut_stand;
  int L;
  int ind_min;
  int *global_phase;
  int ind_max;
  int ind_max_ant;
  int ind_max_landscape;
  int dham;
  int Nmax;
  int *max;
  double roughness;
  double rough_local;
  //  int *path;
  int path_size;
  double *fitness;
  double **traits;
  double **traits_stand;
  double **traits_new;
  double W_max;
  double *trait_init;
  int *nmut;
  int *nmut_stand;
  int *nmut1;
  double *wildtype;
  int *label_newmut;
  int *label_newmut1;
};


struct mutacao{
  double **effect;
  int contlabel;
  double *position;
  double **effect_stand;
  int contlabel_stand;
  double *position_stand;
  
};

void fitnessinit(struct sequencia *sequence, struct populacao *popul, struct otimo *optimum, struct mutacao *mutation);
void fitnessinit_standing(struct sequencia *sequence, struct populacao *popul, struct otimo *optimum, struct mutacao *mutation);
void standing(struct sequencia *sequence, struct populacao *popul, struct otimo *optimum, struct mutacao *mutation);
void dynamics(struct sequencia *sequence, struct populacao *popul, struct otimo *optimum, struct mutacao *mutation);
void fitness_evaluation(struct sequencia *sequence, struct populacao *popul, struct otimo *optimum, struct mutacao *mutation);
void density(struct sequencia *sequence, struct populacao *popul, struct otimo *optimum, struct mutacao *mutation);

int main(int ac, char **av)
{
  FILE *ptt1, *ptt2, *ptt3;

  otimo optimum;

  sequencia sequence;

  populacao popul;

  mutacao mutation;

  int nwalks, max1, cont_path_total, cont_path, cont_newpath, i, j, k, cont_walk, verif, max, conf, config, t, tmax, model_dominant[2], *cont_seq, int_conf, min_step, cont_max, n_landscapes, popul_max, indice, ind_sub, ind_max, cont_extinction, cont_nonextinction, deme, sum_cont_extinction, cont_ind, cont_rescue_total, cont_rescue_partial, cont_rescue, t_eq;

  double fit_max, fit_med, soma1, cont_conf, sum_trait_2, drop_fitness, *trait, sum_t_2, distance_phen, divergence, f_max;

  long seed;  //define a seed

  char arq1[1000], arq2[1000], arq3[1000];

  /*  if (ac < 13 + atoi(av[1]) - 1) {  // Update the condition to reflect the new count
        std::cout << "Start the program like this:\n" << av[0]
                  << " <N> <U> <s> <sd> <alpha> <Nlevel> <tmax> <Tterm> <config> <Nbins> <Nmax> <semente> <K1> <K2> ... <Kn>\n"
                  << std::endl;
        return -1;
	}*/

  if (ac!=9)
    {
      cout  <<  "start the program like this:\n" << av[0]
            << " <N> <U> <s> <sd> <alpha> <Nlevel> <tmax> <Tterm> <config> <Nbins> <Nmax> <semente> \n"
            << endl;
      exit (-1);
    }

  j=0;
  
  popul.K = atoi (av[++j]);
  popul.ntraits = atoi (av[++j]);
  popul.mut = atof (av[++j]);
  n_landscapes = atoi (av[++j]);
  popul.alpha_r = atof (av[++j]);
  popul.sigma_trait_mut = atof (av[++j]);
  drop_fitness = atof(av[++j]); 
  sequence.W_max = atof (av[++j]);

    
  cout  << "#invocation: ";
  for (int i=0; i<ac; i++){
    cout << av[i] << " ";
  }
  cout << endl;

  seed = time (NULL) * getpid();  //set the seed to system time
  gerador=gsl_rng_alloc(gsl_rng_mt19937);  
  gsl_rng_set(gerador, seed); //give the seed to random generator
  
  srand(time(NULL));

  max1 = (int) pow(2.,(double)(sequence.L));
  if( sequence.L==32 )
    max1 = RAND_MAX;

  popul_max = 100000;

  sequence.traits = new double*[2*popul.K];
  for( j=0; j<(2*popul.K); j++ )
    sequence.traits[j] = new double[popul.ntraits];

  sequence.ind_mut = new int*[2*popul.K];
  for( j=0; j<(2*popul.K); j++ )
    sequence.ind_mut[j] = new int[100];

  sequence.traits_stand = new double*[2*popul.K];
  for( j=0; j<(2*popul.K); j++ )
    sequence.traits_stand[j] = new double[popul.ntraits];

  sequence.ind_mut_stand = new int*[2*popul.K];
  for( j=0; j<(2*popul.K); j++ )
    sequence.ind_mut_stand[j] = new int[100];
  
  sequence.traits_new = new double*[2*popul.K];
  for( i=0; i<(2*popul.K); i++ )
    sequence.traits_new[i] = new double[popul.ntraits];
  
  sequence.trait_init = new double[popul.ntraits];

  sequence.label_newmut = new int[2*popul.K];

  sequence.label_newmut1 = new int[2*popul.K];

  sequence.wildtype = new double[popul.ntraits];
  
  optimum.traits = new double[popul.ntraits];

  optimum.trait_shift = new double[popul.ntraits];

  sequence.fitness = new double[2*popul.K];

  sequence.nmut = new int[2*popul.K];

  sequence.nmut_stand = new int[2*popul.K];

  sequence.nmut1 = new int[2*popul.K];

  trait = new double[popul.ntraits];

  mutation.effect = new double*[10000000];
  for( i=0; i<10000000; i++ )
    mutation.effect[i] = new double[popul.ntraits];
  
  mutation.position = new double[10000000];

  mutation.effect_stand = new double*[10000000];
  for( i=0; i<10000000; i++ )
    mutation.effect_stand[i] = new double[popul.ntraits];
  
  mutation.position_stand = new double[10000000];

  config = n_landscapes;
  conf = 0;
  cont_conf = 0;

  sum_cont_extinction = 0;

  cont_extinction = 0;
  cont_nonextinction = 0;
  cont_rescue = 0;
  while( (++conf)<=config )
    {

      popul.N = popul.K;

      sum_trait_2 = -2*log((1-drop_fitness)/sequence.W_max);
      sum_t_2 = 0;
      for( k=0; k<popul.ntraits; k++ )
	{
	  trait[k] = gsl_ran_gaussian(gerador, 1.);
	  sum_t_2 += (trait[k]*trait[k]);
	}

      for( k=0; k<popul.ntraits; k++ )
	trait[k] = trait[k]*sqrt(sum_trait_2)/sqrt(sum_t_2);
      
      for( k=0; k<popul.ntraits; k++ )
	optimum.trait_shift[k] = trait[k];
            
      for( k=0; k<popul.ntraits; k++ )
	optimum.traits[k] = optimum.trait_shift[k];
	      
      fitnessinit_standing(&sequence,&popul,&optimum,&mutation);

      t = 0;
      verif = 0;
      
      while( (popul.N>0) )
	{
	  t++;

	  fitness_evaluation(&sequence,&popul,&optimum,&mutation);

	  dynamics(&sequence,&popul,&optimum,&mutation);

	  density(&sequence,&popul,&optimum,&mutation);

	  cont_ind = 0;

	  for( k=0; k<popul.N; k++ )
	    if( sequence.fitness[k]>1 )
	      cont_ind++;
	  
	  if( (cont_ind>200) )
	    {
	      cont_rescue++;
	      break;
	    }

	  if( (popul.N==0) )
	    {
	      cont_extinction++;
	      break;
	    }	  
	  
	}

    }

 
 
  sprintf(arq2,"STANDING-DeNovo-K%d-Sigma%g-ntraits%d-U%g-Wmax%g.dat",popul.K,popul.sigma_trait_mut,popul.ntraits,popul.mut,sequence.W_max);
  ptt2 = fopen(arq2,"a");
  fprintf(ptt2,"%g \t %g  \t %g  \n",drop_fitness,((double)cont_extinction/n_landscapes),((double)cont_rescue/n_landscapes));
  fclose(ptt2);

}


void standing(sequencia *sequence, populacao *popul, otimo *optimum, mutacao *mutation)
{
  int i, max, max1, ind, j, k;

  for( i=0; i<popul->N; i++ )
    {
      for( j=0; j<popul->ntraits; j++ )
	sequence->traits_stand[i][j] = sequence->traits[i][j];
      sequence->nmut_stand[i] = sequence->nmut[i];
    }
  
  mutation->contlabel_stand = mutation->contlabel;

  for( j=0; j<mutation->contlabel; j++ )
    {
      mutation->position_stand[j] = mutation->position[j];
      for( k=0; k<popul->ntraits; k++ )
	mutation->effect_stand[j][k] = mutation->effect[j][k];
    }
  
}


void fitnessinit_standing(sequencia *sequence, populacao *popul, otimo *optimum, mutacao *mutation)
{
  int i, max, max1, ind, j, k;

  double S_V_G;

  S_V_G = pow((popul->sigma_trait_mut*popul->sigma_trait_mut*popul->mut),0.25);

  for( i=0; i<popul->N; i++ )
    {
      for( j=0; j<popul->ntraits; j++ )
	sequence->traits[i][j] = gsl_ran_gaussian(gerador, S_V_G);
      sequence->nmut[i] = 0;
    }
  
}


void fitnessinit(sequencia *sequence, populacao *popul, otimo *optimum, mutacao *mutation)
{
  int i, max, max1, ind, j, k;


  for( i=0; i<popul->ntraits; i++ )
    {
      sequence->trait_init[i] = 0.;
      sequence->wildtype[i] = 0;
    }

  for( i=0; i<popul->N; i++ )
      {
	for( j=0; j<popul->ntraits; j++ )
	  sequence->traits[i][j] = sequence->trait_init[j];
	sequence->nmut[i] = 0;
      }

  mutation->contlabel = 0;;
  
}


void fitness_evaluation(sequencia *sequence, populacao *popul, otimo *optimum, mutacao *mutation)
{
  int k, indice, max1, deme;
  double delta_r[100], distance, factor;

  factor = sequence->W_max;
  for( indice=0; indice<popul->N; indice++ )
    {
      distance = 0;
      for( k=0; k<popul->ntraits; k++ )
	distance += (optimum->traits[k]-sequence->traits[indice][k])*(optimum->traits[k]-sequence->traits[indice][k]);
      
      sequence->fitness[indice] = factor*exp(-distance/(2*popul->alpha_r*popul->alpha_r));
    }
  
}


void dynamics(sequencia *sequence, populacao *popul, otimo *optimum, mutacao *mutation)
{
  int i, k, m, contbin, kinf, aleat, dig, oper, j, bit, verif, k1, cont, nlabelaux, size, indtau, *ndel, auxind, *ind1, pos, oper1, ind_pilha[2], max1, nmut, n_offspring, **pilha, seq_ind, deme, **aux_ind_mut;

  unsigned int *N_newborn;
     
  double x, soma, aux, soma1, aux1, mut, Fit, *auxF, malp, mutb, mb, r, Fmax, fitness, saux, F_tot, distance_2, med, a, factor, sum_2, y[100], norm, r_mut;

  double *prob;


  aux_ind_mut = new int*[2*popul->K];
  for( i=0; i<(2*popul->K); i++ )
    aux_ind_mut[i] = new int[100];

  cont = 0;
  for( i=0; i<popul->N; i++ )
    {
      fitness = sequence->fitness[i];
      n_offspring = gsl_ran_poisson(gerador,fitness);
      
      for( j=0; j<n_offspring; j++ )
	{
	  for( k=0; k<popul->ntraits; k++ )
	    sequence->traits_new[cont][k] = sequence->traits[i][k];
	  sequence->nmut1[cont] = sequence->nmut[i];
	      
	  r = gsl_ran_flat(gerador, 0., 1.);
	  if( (r < popul->mut) )
	    {
	      for( k=0; k<popul->ntraits; k++ )
		{
		  sequence->traits_new[cont][k] += gsl_ran_gaussian(gerador, popul->sigma_trait_mut);
		}
	      sequence->nmut1[cont]++;
	    }
	  cont++;
	}

    }

  popul->N = cont;
  for( i=0; i<popul->N; i++ )
    {
      sequence->nmut[i] = sequence->nmut1[i];
      for( j=0; j<popul->ntraits; j++ )
	sequence->traits[i][j] = sequence->traits_new[i][j];
    }
  
  
  for( i=0; i<(2*popul->K); i++ )
    delete[] aux_ind_mut[i];
  delete[] aux_ind_mut;  

}



#include <vector>
#include <algorithm>
#include <ctime>

// Assuming appropriate definitions for sequencia, populacao, and otimo

void density(sequencia *sequence, populacao *popul, otimo *optimum, mutacao *mutation)
{
    int j, deme;

    if (popul->N > popul->K)
      {
	// Create an index vector for shuffling
	std::vector<int> indices(popul->N);
	for (int i = 0; i < popul->N; ++i) indices[i] = i;
	
	// Initialize random seed
	std::srand(unsigned(std::time(0)));

	// Random shuffle the indices
	std::random_shuffle(indices.begin(), indices.end());

	// Create new fitness and traits arrays with resized population
	std::vector<double> new_fitness(popul->K);
	std::vector<int> new_nmut(popul->K);
	std::vector<std::vector<double>> new_traits(popul->K, std::vector<double>(popul->ntraits));
	
	for (int i = 0; i < popul->K; ++i)
	  {
	    new_fitness[i] = sequence->fitness[indices[i]];
	    new_nmut[i] = sequence->nmut[indices[i]];

	    for (j = 0; j < popul->ntraits; j++)
	      new_traits[i][j] = sequence->traits[indices[i]][j];
	    
	  }
	    
	    
	// Copy back the trimmed population
	for (int i = 0; i < popul->K; ++i)
	  {
	    sequence->fitness[i] = new_fitness[i];
	    for (j = 0; j < popul->ntraits; j++)
	      sequence->traits[i][j] = new_traits[i][j];
	    sequence->nmut[i] = new_nmut[i];
	  }
	    	    
	popul->N = popul->K;
	
      }
}



