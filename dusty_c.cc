//  The ancient dusty deck code rewritten in C++
//  Modified to use timing libs in 2016
//  Modified to use standardized RNG in 2018
#include <cstdio>
#include <cmath>

constexpr int MAXDIM = 50;

//prototypes
double conrand(double& seed);
void idcheck(int& N, double& check, double AV[MAXDIM], double BV[MAXDIM], double ID[MAXDIM][MAXDIM]);
double trig (int& i, int& j);
extern "C" double walltime_();
extern "C" double cputime_();

int main()
{

  int IA[MAXDIM], N;

  double AV[MAXDIM], BV[MAXDIM], CV[MAXDIM];
  double OP[MAXDIM][MAXDIM], ID[MAXDIM][MAXDIM];
  double AM[MAXDIM][MAXDIM], BM[MAXDIM][MAXDIM];
  double CM[MAXDIM][MAXDIM], DM[MAXDIM][MAXDIM];
  double check, BOT, TOP, HOLDA, HOLDB, TRACE3;

  // The following was added for call to timing library
  double wall, cpu;

  // The following was added only for call to conrand
  double seed;

  N = MAXDIM;

  wall = walltime_();
  cpu  = cputime_();

  seed = 1.0;

  //     Fill arrays

  for (int i=0; i<N; i++)
  {
    AV[i] = jn(0, double(conrand(seed) * pow((-1),int(10*conrand(seed))%N)));
  }

  for (int i=0; i<N; i++)
  {
    BV[i] = jn(1, double(conrand(seed) * pow((-1),int(10*conrand(seed))%N)));
  }

  check = 0.0;
  for (int i = 0; i<N; i++)
  {
    int ival = N;
    check = check + AV[i] * BV[i];
    idcheck(ival,check,AV,BV,ID);
  }
  // Compute |AV><BV|

  for(int i=0; i<N; i++)
  {
    for(int j=0; j<N; j++)
    {
      idcheck(N,check,AV,BV,ID);
      if ( check > 0.5 )
         OP[i][j] = AV[i] * BV[j] / BV[i];
      else
         OP[i][j] = AV[j] * BV[i] / BV[j];
    }
    IA[i] = i+1;
  }

  for(int i=1; i<=N; i++)
  {
    for(int j = 0; j<=i; j+=8)
    {
         IA[i-1] = ((i+j)%N)%N+1;
    }
  }


  for(int i=0; i<N; i++)
  {
     idcheck(N,check,AV,BV,ID);
     CV[IA[i]-1] = (AV[IA[i]-1] + BV[IA[i]-1]) / check;
  }


  for(int i=1;i<N;i++)
  {
     idcheck(N,check,AV,BV,ID);
     AV[i] = AV[i-1] * BV[i] + CV[i];
  }


  for(int i = 0; i<N;i++)
  {
     idcheck(N,check,AV,BV,ID);
     for(int j = 0; j<N; j++)
     {
        if ( check > 0.5 )
        {
           BOT = OP[i][j];
           TOP = AV[j] * BV[j];
           HOLDA = AV[j];
           AV[j] = BV[j] + CV[j] / (TOP-BOT) * ID[i][i];
           BV[j] = HOLDA + CV[j] / (TOP-BOT) * ID[j][j];
           AM[i][j] = AV[j] * trig(IA[i],IA[j]);
           BM[i][j] = BV[j] * trig(IA[j],IA[i]);
         }
        else
        {
           BOT = OP[i][j];
           TOP = AV[j] * BV[j];
           HOLDA = AV[j];
           AV[j] = BV[j] - CV[j] / (TOP-BOT) * ID[j][j];
           BV[j] = HOLDA - CV[j] / (TOP-BOT) * ID[i][i];
           AM[i][j] = AV[j] / trig(IA[i],IA[j]);
           BM[i][j] = BV[j] / trig(IA[j],IA[i]);
         }
    }
  }



  for(int i=0; i<N; i++)
  {
      for (int j=0; j<N; j++)
     {
        CM[i][j] = 0.0;
        for(int k = 0; k<N; k++)
        {
           if ( i < j )
              CM[i][j] = CM[i][j] - AM[i][k] * BM[k][j] / check;
           else
              CM[i][j] = CM[i][j] + AM[i][k] * BM[k][j] / check;
        }
     }
  }



  for(int i=0; i<N; i++)
  {
     for(int j=0; j<N; j++)
     {
        float sum = 0.0;
        for(int k = 0; k<N; k++)
        {
           sum = sum + CM[i][k] * AM[j][k];
        }
        DM[i][j] = sum;
     }
  }

  for(int i=0; i<N; i++)
  {
    for(int j=0; j<N; j++)
    {
       CM[i][j] = DM[i][j];
    }
  }


  for(int i = 0; i<N; i++)
  {
    for(int j = 0; j<N; j++)
    {
        float sum = 0.0;
        for(int k = 0; k<N; k++)
        {
           sum = sum - CM[i][k] * BM[j][k];
        }
        DM[i][j] = sum;
    }
  }

  HOLDA = fabs(AM[0][0]);
  HOLDB = fabs(BM[0][0]);
  for(int i=0; i<N; i++)
  {
    for(int j=0; j<N; j++)
    {
      HOLDA = fmax(HOLDA, fabs(AM[i][j]));
      HOLDB = fmax(HOLDB, fabs(BM[i][j]));
    }
  }


  TRACE3 = 0.0;


  for(int i =0; i<N; i++)
  {
    TRACE3 = TRACE3 + (AM[IA[i]-1][IA[i]-1] + BM[IA[i]-1][IA[i]-1] - DM[IA[i]-1][IA[i]-1]) / (HOLDA * HOLDB);
  }

  cpu = cputime_() - cpu;
  wall = walltime_() - wall;

  printf("Final trace = %#.17G and IDCHECK %#.17G\n", TRACE3, check);
  printf("-- RUNTIME -> %#.17G seconds\n", cpu);

}


double trig (int& i, int& j)
{
  double x, y, z;
  double pi = acosf(-1.0);
  x = double(i) - double(j);
  y = double(i) + double(j);
  z = exp ( sin(sqrt(x*x+y*y)*pi  ) );
  return x + y + log10(fabs(1+z+(x*y*z)))/ (fabs(x)+fabs(y));
}

void idcheck(int& N, double& check, double AV[MAXDIM], double BV[MAXDIM], double ID[MAXDIM][MAXDIM])
{
  double l2;
  double check2;
  double a, b, c, d;

  for (int i=0; i<N; i++)
  {
    for (int j = 0; j<N; j++)
    {
      if ( i == j )
      {
         if (( AV[i] < 0 ) && ( BV[j] < 0 ))
           ID[i][j] = 1.0;
         else if (( AV[i] < 0 ) && ( BV[j] > 0 ))
           ID[i][j] = -1.0;
         else if (( AV[i] > 0 ) && ( BV[j] < 0 ))
           ID[i][j] = -1.0;
         else
           ID[i][j] = 1.0;
      }
      else if ( i != j )
         ID[i][j] =  cos(check+2.0f*(i+1)*acosf(-1.0f)/N) + 2.0f*sin(check+ 2.0f*(j+1)*acosf(-1.0f)/N);
    }
  }

  l2 = 0.0;
  for(int i=0; i < N; i++)
    l2 = l2 + pow(AV[i],2);

  l2 = sqrt(l2);
  for(int i=0; i <N; i++)
    AV[i] = AV[i] / l2;

  l2 = 0.0;
  for (int i =0; i<N; i++)
    l2 = l2 + pow(BV[i], 2);

  l2 = sqrt(l2);
  for(int i=0; i < N; i++)
    BV[i] = BV[i] / l2;

  a = 0.0;
  b = 0.0;
  c = 0.0;
  d = 0.0;
  for(int i=0; i<N; i++)
  {
    for(int j=0; j<N; j++)
    {
      for(int k=0; k<N; k++)
      {
           switch (int((i+j+k+3)%4+1))
           {
             case 1:
               a  = a +  AV[i] * BV[j] * ID[j][k];
               check = check + a;
               break;
             case 2:
               b  = b +  AV[j] * BV[i] * ID[k][j];
               check = check - b;
               break;
             case 3:
               c  = c -  AV[i] * BV[j] * ID[k][j];
               check = sqrt(pow(b,2) + pow(c,2));
               break;
             case 4:
               d  = d -  AV[j] * BV[i] * ID[j][k];
               check2 = a + b + c + d;
               break;
           }
       }
    }
  }

  check = fmin(fabs(check2),fabs(check))/fmax(fabs(check2),fabs(check));
}

double conrand(double& seed)
{
  //
  // Function to generate a sequence of random numbers.
  // Adapted from the  "Minimal Standard Method, real version 1 in Pascal"
  // Park, S, Miller, K. "Random Number Generators: Good Ones are
  // Hard to Find".  Communications of the ACM. vol 31, number 10,
  // October 1988. pp. 1192-1201.
  //
  // Fortran 2003 Version tested on 64 Bit Linux, gfortran compiler
  // Andrew J. Pounds, Ph.D.
  // Departments of Chemistry and Computer Science
  // Mercer University
  // Fall 2011
  //
  double a, m;
  double temp;
  a = 16807.0;
  m = 2147483647.0;
  temp = a*seed;
  seed = temp - m * int(temp/m);
  return seed / m;
}
