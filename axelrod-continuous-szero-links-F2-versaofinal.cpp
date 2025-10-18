// compila-se da seguinte forma c++ -O3 name.cpp -o name -lm -lgsl -lgslcblas

#include<math.h>
#include<iostream>
#include<string.h>
#include<stdio.h>
#include<stdlib.h> 
#include<iomanip> 				
#include<fstream>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_sort.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_sort_double.h>
#include <unistd.h>

using namespace std;
 
#define PI 3.1415927


long idum;
gsl_rng *gerador;

struct rede{
  int **Node;
  int *C;
  int N;
  double v;
  int satisf;
  int L;
  int cont_distance;
  double sum_distance;
  double sum_nsteps;
  int nclusters;
  int largest;
  double diameter;
  double alpha;
  double d_c;
  int **pilha;
  int **links;
  int *nlinks;
  int npilha;
  int second;
};

struct individuo{
  double *posx;
  double *posy;
  int *time_refrat;
  double **state;
  double mu;
  double alpha;
  double beta;
  double R;
  double radius;
  double *heading;
  double *posx_init;
  double *posy_init;
  double *A;
  int A_min;
  int A_max;
  int *nsteps;
  int F;
  double q;
  double sigma_trait_mut;
  double d;

};



void definenet(struct rede *network);
void conf_init(struct individuo *individual, struct rede *network);
void dynamics(struct individuo *individual, struct rede *network);
void HK(struct individuo *individual, struct rede *network);
double fat(int m);

int i1=0;

int main(int ac, char **av)
{

  FILE *ptt1, *ptt2;

  rede network;

  individuo individual;

  int i, j, k, sem, conf, ind, config, verif, pilha[10000], indmaior, t, tmax, cont, cont_measures, activity;

  double x, sum_med_activity, sum_med_nclusters, sum_med_largest, sum_med_diameter, sum_activity, sum_nclusters, sum_largest, sum_diameter, sum_g, sum_g_2, var_g, sum_med_g, sum_med_g_2, med_g, sum_second;

  long seed;  //define a seed

  char arq1[200], arq2[200];

  if (ac!=5)
    {
      cout  <<  "start the program like this:\n" << av[0]
            << " <L> <U> <s> <sd> <alpha> <Nlevel> <tmax> <Tterm> <config> <Nbins> <Nmax> <semente> \n"
            << endl;
      exit (-1);
    }
  
  /**** Leitura de dados por meio de um script ****/

  j = 0;
  network.L = atoi (av[++j]);
  individual.F = atoi (av[++j]);
  individual.d = atof (av[++j]);
  config = atoi (av[++j]);
  
  cout  << "#invocation: ";
  for (int i=0; i<ac; i++){
    cout << av[i] << " ";
  }
  cout << endl;

  seed = time (NULL) * getpid();  //set the seed to system time
  gerador=gsl_rng_alloc(gsl_rng_mt19937);  
  gsl_rng_set(gerador, seed); //give the seed to random generator

 
  network.N = network.L*network.L;

  network.pilha = new int *[4*network.N];
  for( i=0; i<(4*network.N); i++ )
    network.pilha[i] = new int[2];	 
  network.C = new int[network.N];
  network.Node = new int *[network.N];
  for( i=0; i<network.N; i++ )
    network.Node[i] = new int[10];

  network.nlinks= new int[network.N];
  network.links = new int *[network.N];
  for( i=0; i<network.N; i++ )
    network.links[i] = new int[4];
  
  individual.state = new double *[network.N];
  for( i=0; i<network.N; i++ )
    individual.state[i] = new double[individual.F];

  network.sum_distance = 0;
  network.sum_nsteps = 0;
  network.cont_distance = 0;
  sum_med_g = sum_med_g_2 = 0;
  sum_med_activity = sum_med_nclusters =  sum_med_largest = sum_med_diameter = 0;

  conf = 0;
  sum_largest = 0;
  sum_nclusters = 0;
  sum_second = 0;

  definenet(&network);
  while( (++conf)<=config )
    {
      conf_init(&individual,&network);

      t = 0;
  
      while( network.npilha>0 )
	{
	  t++;

	  dynamics(&individual,&network);

	}

      HK(&individual,&network);

      sum_largest += (double)network.largest/network.N;

      sum_nclusters += (double)network.nclusters/network.N;

      sum_second += (double)network.second/network.N;
      
      sprintf(arq2,"Output-PartialOutcomes-Axelrod-L%d-F%d-sigma_zero-d%g.dat",network.L,individual.F,individual.d);
      ptt2 = fopen(arq2,"a");
      fprintf(ptt2,"%g   %g   %d  %g \n",individual.d,((double)network.largest/network.N),network.nclusters,((double)network.second/network.N));
      fclose(ptt2);
 
	  
    }

  sprintf(arq2,"Output-Optimized-Axelrod-L%d-F%d-sigma_zero.dat",network.L,individual.F);
  ptt2 = fopen(arq2,"a");
  fprintf(ptt2,"%g   %g   %g\n",individual.d,(sum_largest/config),(sum_nclusters/config));
  fclose(ptt2);

  
  for (int i = 0; i < (4 * network.N); i++)
    delete[] network.pilha[i];
  delete[] network.pilha;

  delete[] network.C;
  
  for (int i = 0; i < network.N; i++)
    delete[] network.Node[i];

  delete[] network.Node;

  delete[] network.nlinks;

  for (int i = 0; i < network.N; i++)
    delete[] network.links[i];

  delete[] network.links;

  for (int i = 0; i < network.N; i++)
    delete[] individual.state[i];
  delete[] individual.state;
  
}


void definenet( rede *network )
{
  int i, nsitios, k, cont, sitio, m, index, key1, position, position1, M, verif, j;

  double x, y;

  for( i=0; i<network->N; i++ )
    network->C[i] = 4;

  for( i=0; i<network->N; i++ )
    {
      network->Node[i][0] = ( (i<network->L) ) ? (network->N-network->L+i) : (i-network->L);
      network->Node[i][1] = ( ((i%network->L)==0 ) ) ? (i+network->L-1) : (i-1);
      network->Node[i][2] =  ( (i<(network->N-network->L)) ) ? (i+network->L) : (i%network->L);
      network->Node[i][3] =  ( ((i%network->L)==(network->L-1)) ) ? (i-network->L+1) : (i+1);
    }

}


void conf_init(individuo *individual, rede *network)
{
  int i, k, verif, check, j, *attempt, indice, ind_seq, cont, sim;

  double soma, max, x, d_x, d_y, distance, sum_cos, sum_sin, sum_cos_1, sum_sin_1, R, R_1;

  for( i=0; i<network->N; i++ )
    {
      for( j=0; j<individual->F; j++ )
	individual->state[i][j] = gsl_ran_flat(gerador, 0., 1.);
      network->nlinks[i] = 0;
    }

  network->npilha = 0;
  
  for( i=0; i<network->N; i++ )
    {
      cont = 0;
      for( j=0; j<2; j++ )
	{
	  sim = 0;
	  for( ind_seq=0; ind_seq<individual->F; ind_seq++ )
	    if(fabs(individual->state[i][ind_seq]-individual->state[network->Node[i][j]][ind_seq]) <= individual->d)
	      sim++;
	  
	  if( (sim>0) && (sim<individual->F) )
	    {
	      network->pilha[network->npilha][0] = i;
	      network->pilha[network->npilha][1] = network->Node[i][j];
	      network->links[i][network->nlinks[i]] = network->npilha;
	      network->links[network->Node[i][j]][network->nlinks[network->Node[i][j]]] = network->npilha;
	      network->nlinks[i]++;
	      network->nlinks[network->Node[i][j]]++;
	      network->npilha++;
	    }
	}

    }

}




void dynamics(individuo *individual, rede *network)
{
  int pilha[1000], i, k, m, cont, ind_neigh, cont_diff, j, verif, check, indice, sim, soma_1, soma_move, d_ham, maior_sim, ind_maior, pilha_diff[100], ind_seq, chosen, ind_feature, cont_ind, ind_pilha, pilha_aux[5], ind_pilha_aux, npilha_aux, indice_node[5], indice_link, kl, l, indice_target;
  double x, soma, r, indx, indy, p[1000], p_interact, mut_effect, distance_trait, distance;

  FILE *ptt2;

  char arq2[200];
  
  cont = 0;

  ind_pilha = (int)(network->npilha*gsl_ran_flat(gerador, 0., 1.));
  sim = 0;
  cont_diff = 0;
  for( ind_seq=0; ind_seq<individual->F; ind_seq++ )
    {
      distance = individual->d - fabs(individual->state[network->pilha[ind_pilha][0]][ind_seq]-individual->state[network->pilha[ind_pilha][1]][ind_seq]);
      if( distance>=0 )
	sim++;
      else
	{
	  pilha_diff[cont_diff] = ind_seq;
	  cont_diff++;
	}
    }

  /*** para F=2 diferem em apenas um sítio ***/

  ind_seq = pilha_diff[0];
  ind_neigh = (gsl_ran_flat(gerador, 0., 1.)<0.5) ? 0 : 1;
  mut_effect = individual->state[network->pilha[ind_pilha][ind_neigh]][ind_seq];
  individual->state[network->pilha[ind_pilha][1-ind_neigh]][ind_seq] = mut_effect;

  indice_node[0] = network->pilha[ind_pilha][1-ind_neigh];
  indice_node[1] = network->pilha[ind_pilha][ind_neigh];

  indice_target = network->pilha[ind_pilha][1-ind_neigh];

  for( k=0; k<2; k++ )
    {
      for( j=0; j<(network->nlinks[indice_node[k]]-1); j++ )
	if( network->links[indice_node[k]][j]==ind_pilha )
	  {
	    network->links[indice_node[k]][j] = network->links[indice_node[k]][network->nlinks[indice_node[k]]-1];
	    break;
	  }
      network->nlinks[indice_node[k]]--;
    }


  network->pilha[ind_pilha][0] = network->pilha[network->npilha-1][0];
  network->pilha[ind_pilha][1] = network->pilha[network->npilha-1][1];
  
  indice_node[0] = network->pilha[network->npilha-1][0];
  indice_node[1] = network->pilha[network->npilha-1][1];

  for( l=0; l<2; l++ )
    for( k=0; k<network->nlinks[indice_node[l]]; k++ )
      if( network->links[indice_node[l]][k]==(network->npilha-1) )
	{
	  network->links[indice_node[l]][k] = ind_pilha;
	  break;
	}
  network->npilha--;
  
  i = indice_target;
  
  for( j=0; j<network->C[i]; j++ )
    {
      sim = 0;
      for( ind_seq=0; ind_seq<individual->F; ind_seq++ )
	{
	  distance = individual->d - fabs(individual->state[i][ind_seq]-individual->state[network->Node[i][j]][ind_seq]);
	  if( distance>=0 )
	    sim++;
	}
 
      if( (sim!=0) && (sim!=individual->F) )
	{
	  verif = 0;
	  for( k=0; k<network->nlinks[i]; k++ )
	    if( (network->pilha[network->links[i][k]][0]==network->Node[i][j]) || (network->pilha[network->links[i][k]][1]==network->Node[i][j]) )
	      {
		verif = 1;
		break;
	      }
	  if( verif==0 )
	    {
	      network->pilha[network->npilha][0] = i;
	      network->pilha[network->npilha][1] = network->Node[i][j];
	     
	      network->links[i][network->nlinks[i]] = network->npilha;
	      network->links[network->Node[i][j]][network->nlinks[network->Node[i][j]]] = network->npilha;
	    
	      network->nlinks[i]++;
	      network->nlinks[network->Node[i][j]]++;
	      network->npilha++;
	    }	  
	}

      if(  (sim==0) || (sim==individual->F) )
	{
	  verif = 0;
	  for( k=0; k<network->nlinks[i]; k++ )
	    if( (network->pilha[network->links[i][k]][0]==network->Node[i][j]) || (network->pilha[network->links[i][k]][1]==network->Node[i][j]) )
	      {
		verif = 1;
		indice_link = network->links[i][k];
		break;
	      }
	  if( verif==1 )
	    {
	      network->pilha[indice_link][0] = network->pilha[network->npilha-1][0];
	      network->pilha[indice_link][1] = network->pilha[network->npilha-1][1];

	      indice_node[0] = network->pilha[network->npilha-1][0];
	      indice_node[1] = network->pilha[network->npilha-1][1];

	      for( l=0; l<2; l++ )
		for( k=0; k<network->nlinks[indice_node[l]]; k++ )
		  if( network->links[indice_node[l]][k]==(network->npilha-1) )
		    network->links[indice_node[l]][k] = indice_link;
	      
	      network->npilha--;

	      indice_node[0] = i;
	      indice_node[1] = network->Node[i][j];
	      
	      for( kl=0; kl<2; kl++ )
		{
		  for( l=0; l<(network->nlinks[indice_node[kl]]-1); l++ )
		    if( network->links[indice_node[kl]][l]==indice_link )
		      {
			network->links[indice_node[kl]][l] = network->links[indice_node[kl]][network->nlinks[indice_node[kl]]-1];
			break;
		      }
		  network->nlinks[indice_node[kl]]--;
		}
	      
	    }
	  
	}
      

 
      
    }
      
  
}




void HK(individuo *individual, rede *network)
{

  FILE *ptt2;
  
  int *label, *cluster, i, j, indice, *cont_clusters, cont1, labmenor, *pilha_cluster, npilha_clusters, verif, contlabel, maior_cluster, cont_state, cont_totcluster, *pilha, ind_maior, npilha, ind_seq, sim, indice_save, ref, i_ref, contador, contador_label, pilha_label[10000], target, k, second_cluster;
    
  double d_x, d_y, distance, d_max, distance2;

  char arq2[1000];
  
  
  label = new int[network->N];

  cluster = new int[network->N+1];

  cont_clusters = new int[network->N+1];

  pilha_cluster = new int[network->N+1];

  pilha = new int[network->N+1];

  
  /* for( i=0; i<network->N; i++ )
    C_ant[i] = 0;

  for( i=0; i<network->N; i++ )
    {
      for( j=0; j<network->C[i]; j++ )
	{
	  sim = 0;
	  for( ind_seq=0; ind_seq<individual->F; ind_seq++ )
	    {
	      distance = individual->d - fabs(individual->state[i][ind_seq]-individual->state[network->Node[i][j]][ind_seq]);
	      if( distance>=0 )
		sim++;
	    }
	  if( sim==individual->F )
	    {
	      Node_ant[i][C_ant[i]] = network->Node[i][j];
	      C_ant[i]++;
	    }
	}
	}*/

  for( i=0; i<network->N; i++ )
    label[i] = 0;

  contlabel = 1;

  for( i=0; i<network->N; i++ )
    {
      if( label[i]==0 )
	{
	  label[i] = contlabel;
	  contlabel++;
	}

      for( j=0; j<network->C[i]; j++ )
	{
	  if( (label[network->Node[i][j]]>0) && (label[network->Node[i][j]]!=label[i]) && (fabs(individual->state[i][0]-individual->state[network->Node[i][j]][0])<individual->d) )
	    {
	      ref = label[network->Node[i][j]];
	      target = label[i];
	      for( k=0; k<network->N; k++ )
		if( label[k]==ref )
		  label[k] = target;
	    }
	  if( (label[network->Node[i][j]]==0) && (fabs(individual->state[i][0]-individual->state[network->Node[i][j]][0])<individual->d) )
	    label[network->Node[i][j]] = label[i];
	}

    }
  


  

  for( i=0; i<contlabel; i++ )
    cont_clusters[i] = 0;
 
  npilha_clusters = 0;
  for( i=0; i<network->N; i++ )
    {
      cont_clusters[label[i]]++;
      verif = 0;
      for( j=0; j<npilha_clusters; j++ )
	if( pilha_cluster[j]==label[i] )
	  {
	    verif = 1;
	    break;
	  }
      if( verif==0 )
	{
	  pilha_cluster[npilha_clusters] = label[i];
	  npilha_clusters++;
	}
    }

  
 
  maior_cluster = 0;
  cont_totcluster = 0;
  ind_maior = 0;
  second_cluster = 0;
  for( i=0; i<npilha_clusters; i++ )
    {
      if( (cont_clusters[pilha_cluster[i]]>maior_cluster) )
	{
	  second_cluster = maior_cluster;
	  maior_cluster = cont_clusters[pilha_cluster[i]];
	  ind_maior = pilha_cluster[i];
	}
      else if( (cont_clusters[pilha_cluster[i]]>second_cluster) )
	second_cluster = cont_clusters[pilha_cluster[i]];
      
      cont_totcluster += cont_clusters[pilha_cluster[i]];
    }
  



  network->nclusters = npilha_clusters;
  network->largest =  maior_cluster;
  network->second = second_cluster;
  //  network->diameter = d_max;

  
  delete[] cluster;
  delete[] label;

  delete[] cont_clusters;
  delete[] pilha_cluster;
  delete[] pilha;

 
	  
}
	     
