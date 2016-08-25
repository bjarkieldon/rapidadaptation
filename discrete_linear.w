\pdfoutput=1
\documentclass[12pt]{cweb}%
\usepackage[english]{babel}
\usepackage[T1]{fontenc}
\usepackage{amsfonts, amsmath, amssymb}
\usepackage{tikzsymbols}
\usepackage{fullpage}
%\usepackage{bbold}%
\usepackage{bm}
\usepackage{url}
\usepackage{natbib}
\usepackage[all]{xy}%
\usepackage{color}%
\usepackage{rotating}%
\usepackage{a4wide,fullpage}%
\usepackage{setspace}%
\usepackage{enumerate}%
\usepackage{wasysym}
\usepackage{textcomp}
\usepackage{hyperref}
\setstretch{1.5}%
\newcommand{\one}[1]{\ensuremath{\mathbb{1}_{\left( #1 \right)}  } }%%
\newcommand{\g}{\,\boldsymbol{|}\,}%
\newcommand{\EE}[1]{\mathbb{E}\left[ #1 \right]}%
\newcommand{\pr}[1]{\ensuremath{\mathbb{P}\left( #1 \right) } }%%
\newcommand{\im}{\ensuremath{\imath} }%
\newcommand{\jm}{\ensuremath{\jmath} }%
\newcommand{\T}{\ensuremath{\mathcal{T}}}%
\newcommand{\be}{\begin{equation}}%
\newcommand{\ee}{\end{equation}}%
\newcommand{\tG}{\scalebox{1.2}{\tt G}}%
\newcommand{\tB}{\scalebox{1.2}{\tt B}}%
\newcommand{\tK}{\scalebox{1.2}{\tt K}}%
\newcommand{\norm}[2]{\ensuremath{\boldsymbol{|}{#1}\boldsymbol{|}_{#2} } }%
\newcommand{\In}{\ensuremath{\mathcal{I}_n} }%
\newcommand{\R}{\ensuremath{\mathbb{R}} }%
\newcommand{\bd}{\begin{displaymath}}%
\newcommand{\ed}{\end{displaymath}}%
\newcommand{\uZ}{\ensuremath{\underline{Z}}}%
\newcommand{\uR}{\ensuremath{\underline{R}}}%
\newcommand{\uX}{\ensuremath{\underline{X}}}%
\newcommand{\uY}{\ensuremath{\underline{Y}}}%
\newcommand{\cleanup}{ // clear all used memory: }%
\newcommand{\bone}[1]{\ensuremath{ \mathbb{1}\left( #1 \right) } }%
\newcommand{\EP}[2]{\ensuremath{\mathbb{E}^{#1}\left[ #2 \right] } }%
\newcommand{\hj}{\ensuremath{\hat{\jm}}}%
\title{\bf Discrete-time simulator for multi-loci selection on a trait }
\author{ CWEB technical report\\ bjarki eldon\\ Museum f\"ur Naturkunde \\ 
{Leibniz Institut f\"ur Evolutions- und
Biodiversit\"atsforschung} \\  Berlin, Germany}
\date{\today }%
\begin{document}
\maketitle

\begin{abstract}
Simulates the time (number of generations)  it takes a diploid  population to reach a predefined optimum trait
value.  Selection acts on the (multi-locus) trait,  and we compare the time to
reach optimum for  two reproduction schemes.    This CWEB \citep{knuth92literateprog,knuth1989cweb,KL2001}
technical report describes corresponding C code.  CWEB documents may be compiled with {\tt cweave} and {\tt
ctangle}.
 
\end{abstract}


\tableofcontents


@* {\bf Copyright}. 

Copyright {\copyright} 2016 Bjarki Eldon \newline

This document and any source code it contains  is distributed under the terms of the GNU General Public Licence (version $\ge 3$).  You
should have received a copy of the licence along with this file (see file COPYING).  


    The source codes  described in this document  are  free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This document and the code it contains   is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this file (see COPYING).  If not, see \url{http://www.gnu.org/licenses/}.




@* {\bf Introduction}. 

We are interested in investigating if  high fecundity and skewed distribution
of offspring (HFSDO), in addition to selection acting on a multi-loci trait,
drives rapid adaptation.  In particular,  does the combination of HFSDO and
selection acting on a multi-loci trait  enable  more  rapid  adaptation than
in the absence of  HFSDO.  Also, in the single-locus case,  in both haploid
and diploid case, we would like to know if  HFSDO enables more rapid
adaptation than a non-skewed reproduction mechanism. 



Our population size is constant at  $N$.  In each generation,  each  pair
independently  contributes offspring. We let $X_i$ denote the random number of
offspring of pair $i$, and the $X_i$ are iid, and independent between
generations.  Let $b$, $\alpha$, $\psi$ be positive constants; define
$\bone{A} :=1$ if $A$ holds, and zero otherwise; write $ [n]_0  :=  \{0,1,
\ldots, n\}$ for $n \in \mathbb{N} := \{1, 2, \ldots \}$; $[n] := \{1, 2,
\ldots, n\}$.       The law of $X_i$ is said to be  "skewed" if it is given by  
\be\label{eq:skewed}
  \pr{X_i = k} =  \frac{b^\alpha}{(b + k)^\alpha} -  \frac{b^\alpha\bone{0 \le k < \psi}}{(1+b+k)^{\alpha}}, \quad k \in [\psi]_0;
\ee 
the parameter  $b$ shifts the distribution to allow for mass at zero. The
distribution \eqref{eq:skewed} is monotonically decreasing with increasing
$k$. The law \eqref{eq:skewed} is similar in form to the model by   \cite{schweinsberg2003coalescent}.  However, the model  by  \cite{schweinsberg2003coalescent} is a limit statement, and as such can only be used  in  theoretical calculations.   




We compare the law \eqref{eq:skewed} of $X_i$ with a Poisson distribution with
parameter $\lambda > 0$. For $\alpha < 2$ the law \eqref{eq:skewed} has a much
heavier tail than the Poisson.  The population size is constant at $2N$
generations; $N$ pairs of individuals independently contribute juveniles. The
total number of juveniles generated in any given generation therefore needs to
be at least $2N$.  If the law of $X_i$ is Poisson, we take $\lambda \ge 2$;
otherwise the sampling takes a long time since the total pool of juveniles is
regenerated until the total count of juveniles is at least $2N$.


Given the allelic types  $g_{i,j}^{(1)}$ and    $g_{i,j}^{(2)}$ at locus $j
\in [L]$  of individual $i \in [2N]$,    the trait value of individual $i$ is computed as
\be\label{eq:trait} z_i = \frac{1}{L}\sum_{j \in [L]} \xi_j\left(
g_{i,j}^{(1)}  +   g_{i,j}^{(2)}  \right). \ee 
The trait value of the population is the average of the individual
trait values.   


The optimum trait value is denoted by $z_0$.  
 Given the trait value $z_i$  of individual $i$, the fitness is given by 
\be\label{eq:fitness}
     w_i  =  \exp\left( -s ( z_i - z_0 )^2  \right).   
\ee
 
Equation \eqref{eq:trait} for the trait value  allows one, for
example, to consider linear effects of allelic types.  By way of example,   set
$\xi_j = 1$ for simplicity, and let homozygotes for type 0 be the fittest
genotype, heterozygotes less fit, and homozygotes for type 1  the least fit.
Then $z_0 = 0$, and if there's one locus, the possible  fitness values are 1,
$e^{-s}$, and $e^{-2s}$. 


If, in any given generation, the number of juveniles $(\sum_i x_i)$
exceeds $2N$, we draw $\sum_i x_i$ independent exponentials $(T_i)$,
where $T_i$ has rate $w_i$.  The $2N$ juveniles with the lowest times
then form the next set of adults.  The genotypes of the parents do not
affect the number of juveniles contributed, rather indirectly the
viability of the juveniles.  The viability of a juvenile depends on
its' genotypes.   


@* {\bf Compile and run }.  

Use {\tt cweave} on the {\tt .w} file to generate {\tt .tex} file, and {\tt
ctangle} to generate a {\tt .c} file.  To compile an executable:
\begin{verbatim}
gcc <flags> discrete_linear.c -lm -lgsl -lgslcblas
\end{verbatim}




The command to run simulations (assuming {\tt a.out} is the executable) is 
\begin{verbatim}
a.out <N> <L> <a> <b> <p> <s> <z> <e> <R> <r> <f>
\end{verbatim}
where {\tt N} is number of pairs,  {\tt L} is number of loci,  {\tt a} is
parameter $\alpha$;  {\tt b} is parameter  $b$, the shift of the distribution;   {\tt p} is $\psi$,  the
offspring cutoff;   {\tt s} is the selection strength $s$;  {\tt z} the
optimum $z_0$;   {\tt e} is the allowed distance $\varepsilon$ from optimum;
{\tt R} is number of runs;  {\tt r} is a random seed,     {\tt f} is the name of the file  storing the results, the
time to reach optimum. If  $\alpha > 0$, the law of $X_i$ is assumed skewed;
otherwise set  $\alpha = 0$ to sample from the Poisson with   $b$  taken as
the parameter $(\lambda)$ for  the  Poisson. 

By way of example,   
\begin{verbatim}
a.out 50000 1 0. 2.1 10 10. 0. 0.000001 5 123 results.out
\end{verbatim}
writes the following into {\tt results.out}:
\begin{verbatim}
58
58
60
57
56
\end{verbatim}

@* {\bf Code}.

The configuration of the population is stored as a  $(2L)\times (2N)$ matrix,
with rows denoting loci, and  each column stores the types of one diploid
individual.   Rows  1 to $L$ store the types  an individual received
from one parent, and  rows $L+1$ to $2L$ store the types received from the
other parent.    Rows $j$ and $j+L$ for $j \in [L]$ hold the types at
locus $j$.   The configuration of  juveniles is stored in the same way.



The GSL library  (see \url{https://www.gnu.org/software/gsl/}) is used   for
random number generation and sorting.   If the number of  juveniles exceeds
$2N$ in any given generation,     the exponential times of the juveniles are
sorted into ascending order, and  the $2N$ juveniles with the lowest times
form the new set of adults.  We set  a maximum allowed number of juveniles for
computer memory purposes.   


@*1 {\bf Random number generator }. 

A  random number generator of choice is declaired  using the
$GSL\_RNG\_TYPE$ environment variable.  The default generator is the
`Mersenne~Twister' random number generator \cite{matsumoto1998mersenne} as
implemented in  GSL.


@<random number generator@>=@#

@t declare the random number generator $rngtype$ @>@#
gsl_rng * rngtype ; @#
@t Define the function $setup\_rng$ which initializes $rngtype$:@>@#
void setup_rng( unsigned long int seed ) @#
{@#
   @q *gsl_rng_alloc( const gsl_rng_type *T ) @>
  @q const gsl_rng_type *T ; @>
   
   @q T = gsl_rng_default ; @>
  @q rngtype = gsl_rng_alloc( T ) ; @>
  
  @t set the type as $mt19937$ @>@#
  rngtype = gsl_rng_alloc(gsl_rng_mt19937); @#
  gsl_rng_set(rngtype,  seed); 
  
      gsl_rng_env_setup();  @#
   @q gsl_rng_default_seed = rngseed ; @>
  @q  printf("seed %lu\n", gsl_rng_default_seed); @>
}


@*1 {\bf Definitions}. 

We set the maximum number of juveniles for  memory purposes.   

@<object definitions@>=@#
#define MAX_JUVENILES 10000000


@*1 {\bf Draw values for $X_i$ }. 

 Draw values for $X_i$; the diploid juveniles. We have, for $b, \alpha > 0$ some constants,  
\bd
  \pr{X_i = k} =  \frac{b^\alpha}{(b + k)^\alpha} -  \frac{b^\alpha\bone{0 \le k < \psi}}{(1+b+k)^{\alpha}}, \quad k \in [\psi]_0, 
\ed 
and we observe that $\sum_{0 \le k\le \psi } \pr{X_i = k} = 1$.   

@<initialize distribution for $X_i$@>=@#
void drawXi( int N,   int psi, double a,  double b, gsl_ran_discrete_t *
Pmass,    int * tXi,   gsl_rng * r)
{@#
          @t \newline @>
          /* $N$ is number of pairs;  $psi$ is the truncation of the skewed   offspring
          distribution \eqref{eq:skewed};  the  number of juveniles
          contributed by each pair is stored in  $tXi$;  $Pmass$ is  the distribution \eqref{eq:skewed}  */
           @t \newline @>
        int k, teljari ; @#
    tXi[0] = 0; @#
    teljari = 0 ; @#
    @t \newline @>
       /* we need the total number of juveniles to  be at least $2N$  */
         @t \newline @>
  while( (tXi[0] < 2*N)  &&  (teljari < 1000000) ){@#
  teljari = teljari + 1 ; @#
   tXi[0] = 0; @#
   for( k = 1 ; k <= N ; k++){@#
     @t \newline @>
     /* if $a > 0$ we assume the skewed law \eqref{eq:skewed}; otherwise
     Poisson with parameter $\lambda = b$ */
      @t \newline @>
     tXi[k] = (a > 0. ? (int)gsl_ran_discrete( r, Pmass) : gsl_ran_poisson(r, b) ) ; @#
     tXi[0] = tXi[0]  +  tXi[k]; }} @#


      @t \newline @>
         /* $assert$ breaks the execution if the statement does not hold */
             @t \newline @>
 assert( teljari < 1000000); @#
  
     assert( tXi[0] <= MAX_JUVENILES); @#
}




@*1 {\bf Update population}. 

Update population given values $x_i$ of juveniles generated by each
pair.  The genotypes of the parents do not affect the number of
juveniles contributed.  Genotypes are assigned to each juvenile, and
trait value computed.  If the number of juveniles exceeds the
population size $2N$, selection sorts out the $2N$ fittest (on
average) individuals.  We randomly pair the parents; since we are not
modeling fecundity selection, we can separately  draw $N$ numbers of juveniles,
and then assign to randomly formed parent pairs. 


The module $juvenile\_genotypes$ returns the trait value of the $2N$
juveniles that form the new set of  adults.  



@<assign juvenile genotypes@>=@#

double juvenile_genotypes( int N,  int L, double s,  double znull,  double epsilon, int ** Pop, int
* indindex, int * tXi,  gsl_matrix_int * tempJuve,   double *Z,
double *locuseffects,  double *etimes, size_t * aindex, int * allelefr,   gsl_rng * r )
{@#

          @t \newline @>@#
          /*  $N$ is number of pairs;  $L$ is  number of loci;  $s$ is
          selection coefficient;  $znull$ is  optimum trait value;  $epsilon$ is maximum allowed distance from optimum;  $Pop$ is the  matrix for the population configuration of types;  $indindex$ is  vector of indexes for  pairing the parents;  $tXi$ stores the $N$  numbers of juveniles;  $Z$ stores the juvenile trait values;  $tempJuve$  stores the  type configuration of the juveniles;  $locuseffects$ is the vector of the effects of loci, the      $\xi_j$ in  \eqref{eq:trait};  $etimes$ stores the exponential times; $aindex$  storest the sorted indexes of juveniles;  $allelefr$ stores  the count of an allele of some type, for consistency check;  $r$ is the random number generator    */
            @t \newline @>@#


  int i, j, k,  B, xindex; @#
  double Zbar = 0. ; @#
  allelefr[0] = 0; @#

  @t \newline @>@#
  /* $tXi[0] =  x_1 + \cdots +  x_N$ is the total number of juveniles */
    @t \newline @>@#
    xindex = 0 ; @#
  for( i = 1 ; i <= N ; i++){ @#
    @t \newline @>@#
    // check if  pair $i$ produced potential offspring, ie.\ if $X_i > 0$
      @t \newline @>@#
      
   if( tXi[i] > 0 ){ @#
    for( k = 1 ; k <= tXi[i] ; k++){ @#
        for( j = 1 ; j <= 2*L ; j++){ @#
          @t \newline @>@#
          // if $B$ takes value 0 we assign type from set $[L]$, otherwise from set $\{L+1, \ldots , 2L\}$; genotypes for range $[L]$ are assigned from  individual  $\sigma(i)$, otherwise from individual  $\sigma(N+i)$, where $\sigma$ denotes the permutation of the indexes
             @t \newline @>@#
            B =   (int)gsl_ran_bernoulli( r, .5); @#
            
            gsl_matrix_int_set( tempJuve,  j, xindex + 1, Pop[ (j > L ? (B > 0
              ? j : j -L) : (B > 0 ? j + L : j))][ (j > L ? indindex[N+i] :
              indindex[i]) ] ) ;} @#
              @t \newline @>@#
              // now have assigned genotypes for juvenile, compute trait value 
                 @t \newline @>@#
                
                Z[xindex] = 0. ; @#
              for( j = 1 ; j <= L ; j++){ @#
                
              Z[xindex]  =  Z[xindex]  +  ( locuseffects[j] * ( (double)(
                 gsl_matrix_int_get( tempJuve, j, xindex+1) +
                 gsl_matrix_int_get( tempJuve, j + L, xindex+1))) / ((double)L) ) ; } @#
                 @t \newline @>@#
                 // given trait value $z$ for juvenile, compute the fitness value $w = e^{-s(z - z_0)^2}$, and draw exponential with rate $w$
                   @t \newline @>@#
                   etimes[xindex] = gsl_ran_exponential(r,  1./gsl_sf_exp( -s*  gsl_pow_2( Z[xindex]  -  znull) ));    @#
                  
            xindex = xindex + 1; } } } @#

            assert ( xindex == tXi[0] ); @#

    @t \newline @>@#

    assert( tXi[0] >= N+N) ; @#

     if( tXi[0] > 2*N){

       @t \newline @> @# 
        /*  sort the times into ascending order, and  store the  sorted indexes in $aindex$ */
           @t \newline @> @# 
      gsl_sort_index( aindex,  etimes,  1,  tXi[0]); @#

       @t \newline @>
       /* the first $2N$ indexes in $aindex$  are the indexes of the surviving juveniles */
          @t \newline @>
           for( i = 0 ; i < N+N ; i++){ @#
          
                Zbar =  Zbar +   Z[ (int)aindex[i] ]/( (double)(N+N) ) ; @#
              
               for( j = 1 ; j <= L+L ; j++){ @#
               allelefr[0] =  allelefr[0] + (gsl_matrix_int_get( tempJuve, j, (int)aindex[i]) == 0 ? 1 : 0) ;
                 Pop[ j][ i+1] =  gsl_matrix_int_get(tempJuve, j, (int)aindex[i] ) ; } } }
                 else{ @# 
                 @t \newline @>
                   /* $tXi[0] = 2N$,  so all  juveniles survive */ 
                     @t \newline @>
                     for( i = 0 ; i < N+N ; i++){ @#
                Zbar =  Zbar +   Z[i] / ( (double)(N+N) ) ; @#
               for( j = 1 ; j <= L+L ; j++){ @#
               allelefr[0] =  allelefr[0] + (gsl_matrix_int_get( tempJuve, j, i) == 0 ? 1 : 0) ;
                 Pop[ j][ i+1]  =  gsl_matrix_int_get(tempJuve, j, i + 1 ) ; } } }

                 return( Zbar ) ; @#                 
}







@*1 {\bf Simulator}. 

A discrete-time simulator for a single run.   

@<discrete-time simulator@>=@#
int   simulator( int N,  int L, double a, double b,  int Psi, double s, double znull,  double epsilon,   
double * leffects,    gsl_ran_discrete_t * PmXi,  int ** Pop,    gsl_rng * r)
{@#

   // $P$ is the population;  $2N$ is number of diploid individuals;  $L$ is number of unlinked loci;  $a$ is  the  skewness parameter; $gamma = \boldsymbol{\gamma} = (\gamma_1, \ldots, \gamma_L)$ is the  vector of values of locus effects.  First we initialize the population $P$.
     @t \newline @>
     // $Pop$ is  the current genotypes of the $2N$ diploid individuals; $tempXi$ is the random  number of juveniles per parent pair; $pindex$ is the  index of  individuals used to form pairs; 
       @t \newline @>

    @q  gsl_matrix_int * Pop = gsl_matrix_int_calloc( L +L + 1,  N + N  + 1); @>
    

    int * tempXi = (int *)calloc( N+1, sizeof(int)); @#
    size_t * jindex = (size_t *)calloc( MAX_JUVENILES, sizeof(size_t) ) ; @#
    int * pindex = (int *)calloc( N + N, sizeof(int)); @#

    gsl_matrix_int * gjuveniles =  gsl_matrix_int_calloc( L + L + 1, MAX_JUVENILES + 1); @#
    
    double  * traits  =  (double *)calloc( MAX_JUVENILES, sizeof(double));  @#
    double  * Etimes  =  (double *)calloc( MAX_JUVENILES, sizeof(double));  @#
    
  int i, j, B,  iter ; @#
  int * allfr =  (int *)calloc(1, sizeof(int)); @#

  for(i = 1; i <= 2*N ; i++){@#
  pindex[i-1] = i ; @#
  for( j = 1 ; j <= 2*L ; j++){ @#
  @t \newline @>
    // $\pr{B = 1} = \tfrac 12$; ie.\ we assign initial alleles with equal probability 
     @t \newline @>
     B = (int)gsl_ran_bernoulli( r, .5); @#

     Pop[ j][ i] = B; }} @#
       @t \newline @>

       double zbar = 0. ; @#
       iter = 0; @#

       for( i = 1; i <= 2*N ; i++){ @#
        for( j = 1 ; j <= L ; j++){@#
          zbar = zbar  +  (leffects[j]* ( Pop[ j][i] +
         Pop[ j+L][i] ) ); }} @#
         @t \newline @>
          // print out the initial value of the mean trait 
            @t \newline @>
        zbar = zbar/( (double)(N+N)); @#
        @q  printf("initial  %g\n", zbar); @>   

      while( ((iter < 100000) && ( fabs( zbar - znull) > epsilon )) ){ @#
         @t \newline @>
  // Pair the diploid individuals, the parents: first shuffle the index values
    @t \newline @>
  gsl_ran_shuffle( r, pindex, 2*N, sizeof(int)); @#

  drawXi(N,  Psi, a, b,  PmXi,  tempXi, r); @#    
        @t \newline @>
        // now have values of $X_i$ for all $N$ pairs;  update population
          @t \newline @>
          
         @q   printf("%d  %d\n", tempXi[0], 2*N) ; @>

          zbar =   juvenile_genotypes(N, L, s, znull, epsilon, Pop, pindex, tempXi, gjuveniles,
           traits, leffects, Etimes, jindex,  allfr,   r); @#             
            iter = iter + 1 ; @#

              
            @q  printf("allfr  %d\n",  ((allfr[0]  > 0) && (allfr[0]  < 4*N)) ? 0 : 1);  @>
            
            @t \newline @>
            /* only need to check if  all alleles are 1  */
              @t \newline @>
            iter = ( (allfr[0]  >  0)  ?  iter : 1000000 ) ; @#  
           
              } @#
              printf("allfr  %d  %d\n", allfr[0], iter );  @#

            @t \newline @>
            // free used memory
             @t \newline @>
             @q gsl_matrix_int_free( Pop); @>

             free( tempXi); @#
             free(jindex); @#
             free( pindex); @#
             gsl_matrix_int_free (  gjuveniles); @#
           free(  allfr ); @#
             free( traits); @#
             free( Etimes); @#

               @t \newline @>
               /* return the count of  generations needed to reach optimum */
                  @t \newline @>

            return( iter ) ; @#
}



@*1 {\bf many runs}. 

Generate many replicates and  the time, number of generations,   to optimum.


@<many runs@>=@#
void replicates( int N, int L, double a, double b,  int psi,  double
seleccoef,  double  znull, double epsilon,     int nruns,   double *
loceffects,  gsl_rng * r ,      char skra[200] )
{ @#

   double * PXi =  (double *)calloc( 1 + psi, sizeof(double)); @#
   int ** mPop = alloc_2d_array( L, N); @#
   int k ; @#
  double mean = 0. ; @#
  for( k = 0 ; k <=  psi ; k++){ @#
    PXi[k] =  (a > 0. ? pow(b, a)*(  pow( 1./( ((double)k) + b), a) -  (k <
  psi ?  pow( 1./( ((double)(1+k))  + b), a) : 0.)) : 1. ) ; @#

    assert( PXi[k] >= 0.) ; @#

    mean = mean  +  ( ((double)k) * PXi[k] ); } @#
    
    assert( mean >= 2. ); @#
 
   gsl_ran_discrete_t *  Pmass = gsl_ran_discrete_preproc( 1 + psi, PXi); @#

   FILE * f = fopen(skra, "w");

  int rep, svar;  @#
  for( rep = 0 ; rep < nruns; rep++){ @# 
        svar = simulator( N,  L, a, b,  psi, seleccoef,  znull,  epsilon,  loceffects,  Pmass,  mPop, r); @#
        fprintf(f,  "%d\n", svar ) ; @#
        printf("%d\n", rep); @#
    }@#

     free( PXi); @#
      gsl_ran_discrete_free( Pmass); @# 
      fclose(f); @#
     
     free_2d_array( mPop,  L); @#
}


@*1 {\bf Allocate a 2d array }. 

Allocate a 2d array; use it for large population size.


@<allocate array@>=@#
int  **alloc_2d_array(  int rows, int cols)
{ @#
  
  int ** m = NULL ; @#
  m = (int **)calloc( rows + rows + 1, sizeof( int *)) ; @#
  int i ; @#
  for( i = 0 ; i <= rows + rows ; i++ ){ @#
     m[i] = (int *)calloc( cols +  cols + 1, sizeof( int ) ) ; } @#
  
  return( m ) ; @#
}


@*1 {\bf Free a 2d array }.  

Free a 2d array.

@<free 2d array@>=@#
void free_2d_array( int **a, int n )
{@#
  
  int i ; @#
  for( i = 0 ; i <= n + n; i++ ){ @#
    free( a[i] ); } @#
  free( a ) ; a = NULL ; @#
}

 
@*1 {\bf the $main$ function}. 

@C

@<Includes@>@#

@<print matrix@>@#

@<allocate array@>@#

@<free 2d array@>@#

@<random number generator@>@#

@<object definitions@>@#

@<initialize distribution for $X_i$@>@#

@<assign juvenile genotypes@>@#

@<discrete-time simulator@>@#

@<many runs@>@#

int main(int argc, char * argv[])@#
{@#

 @t initialise the random number generator @>@#
 setup_rng( (unsigned long int)atoi(argv[10]) ) ; @#

 double * leffects = (double *)calloc( atoi(argv[2]) + 1, sizeof(double)); @#
 
 int ell;@#
 for( ell = 1 ; ell <=  atoi(argv[2]); ell++){@#
   leffects[ell] = 1. ; }@#
   leffects[1] = 1. ; @#
 
 replicates( atoi(argv[1]),  atoi(argv[2]),  atof(argv[3]),  atof(argv[4]),  atoi(argv[5]), atof(argv[6]), atof(argv[7]),  atof(argv[8]),  atoi(argv[9]), leffects,    rngtype,   argv[11] );@#

  @q void simulator( int N,  int L, double a, double b,  int Psi,  double selectcoef,   double * leffects, gsl_rng * r) @>
  @q  simulator( atoi(argv[1]),  atoi(argv[2]),  atof(argv[3]),  atof(argv[4]),   atoi(argv[5]),  atof(argv[6]), leffects, rngtype); @>

@t \cleanup @>@#

gsl_rng_free( rngtype ) ; @#
free( leffects); @#

return GSL_SUCCESS;  @#

@t end of $main$ function @>@# 
}


@*1 {Print a GSL integer matrix  }. 

print matrix

@<print matrix@>=@#
void printmatrix( int rows, int cols,  gsl_matrix_int * m )
{

        int i,j; @#
        for(i  = 1 ; i <= rows; i++){ @#
          for( j = 1 ; j <= cols ; j++){ @#
            printf("%d ", gsl_matrix_int_get( m, i, j)); } @#
              printf("\n"); } @#
}


@* {\bf Includes}. 

@<Includes@>=@#
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_sf_pow_int.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_sf_elementary.h>
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_fit.h>
#include <gsl/gsl_multifit_nlin.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_sf_exp.h>
#include <gsl/gsl_sf_log.h>
#include <gsl/gsl_sf_expint.h>
#include <gsl/gsl_combination.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_combination.h>
#include <gsl/gsl_statistics_double.h>
#include <gsl/gsl_statistics_int.h>
#include <gsl/gsl_sort.h>
#include <assert.h>


@* {\bf Funding}. 

Funding provided by DFG grant STE 325/17-1 to Wolfgang Stephan through the DFG
Priority Programme SPP 1819 "Rapid Evolutionary Adaptation" \\
(\url{https://dfg-spp1819.uni-hohenheim.de/105254?L=1}).   The generous
support of Museum f\"ur Naturkunde in Berlin is  warmly acknowledged. 

@* {\bf References}. 


\bibliographystyle{plain}
\bibliography{references}


@
\end{document}
