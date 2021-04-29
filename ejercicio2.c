#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <pthread.h>

#define PI 3.14159265358979323846

#define BY_ROW(row, col, len) ((row) * (len) + (col))
#define BY_COL(row, col, len) ((row) + (col) * (len))

int N;
int BS;
int Th;

double *A;
double *B;
double *Cp;
double *Cs;
double *T;
double *R1;
double *R2;
double *M;
double *r1a;
double *r2b;

double avgR1;
double avgR2;

pthread_mutex_t mutex;
pthread_barrier_t barrier_avgr, barrier_c;

/* Time in seconds from some point in the past */
double dwalltime();

void *worker(void *ptr);

void *calculo_secuencial();

double randFP(double min, double max)
{
  double range = (max - min);
  double div = RAND_MAX / range;
  return min + (rand() / div);
}

int main(int argc, char *argv[])
{
  // Controla los argumentos al programa
  if (argc != 4 || (N = atoi(argv[1])) <= 0 || (BS = atoi(argv[2])) <= 0 || (N % BS != 0) || (Th = atoi(argv[3])) <= 0 || (Th != 1 && Th % 2 != 0) || Th > 8)
  {
    printf("\nError, modo de uso: %s N BS T (N debe ser multiplo de BS y T debe ser multiplo de 2 y menor o igual a 8)\n", argv[0]);
    printf("Parametros: Matriz de %d x %d, con bloques de %d y %d threads\n", N, N, BS, Th);

    return 0;
  }
  printf("Calculando algoritmo paralelo empleando Pthreads. Matriz de %d x %d, con bloques de %d y %d threads\n", N, N, BS, Th);

  // Reservar memoria para las matrices y auxiliares
  A = (double *)malloc(sizeof(double) * N * N);
  B = (double *)malloc(sizeof(double) * N * N);
  Cp = (double *)malloc(sizeof(double) * N * N);
  Cs = (double *)malloc(sizeof(double) * N * N);
  T = (double *)malloc(sizeof(double) * N * N);
  R1 = (double *)malloc(sizeof(double) * N * N);
  R2 = (double *)malloc(sizeof(double) * N * N);
  M = (double *)malloc(sizeof(double) * N * N);
  r1a = (double *)malloc(sizeof(double) * N * N);
  r2b = (double *)malloc(sizeof(double) * N * N);

  avgR1 = 0;
  avgR2 = 0;

  pthread_t *threads = (pthread_t *)malloc(Th * sizeof(pthread_t));

  pthread_mutex_init(&mutex, NULL);

  pthread_barrier_init(&barrier_avgr, NULL, Th);
  pthread_barrier_init(&barrier_c, NULL, Th);

  int i, j;
  int *ids = (int *)malloc(sizeof(int) * Th);

  time_t t;
  srand((unsigned)time(&t));

  for (int i = 0; i < N; i++)
  {
    for (int j = 0; j < N; j++)
    {
      // Inicializar matrices A, B, T
      A[BY_COL(i, j, N)] = randFP(0, 1);
      B[BY_COL(i, j, N)] = randFP(0, 1);
      T[BY_ROW(i, j, N)] = randFP(0, 1);

      // La matriz M tiene rango 0 a 2*PI
      M[BY_ROW(i, j, N)] = randFP(0, 2 * PI);
    }
  }

  // Llamado pthread
  double timetick = dwalltime();

  for (i = 0; i < Th; ++i)
  {
    ids[i] = i;
    pthread_create(&threads[i], NULL, worker, &ids[i]);
  }

  for (i = 0; i < Th; ++i)
  {
    pthread_join(threads[i], NULL);
  }
  double ptime = dwalltime() - timetick;
  printf("Tiempo Paralelo Pthread en segundos %f \n", ptime);

  //Llamado secuencial
  timetick = dwalltime();

  calculo_secuencial(0);

  double stime = dwalltime() - timetick;

  printf("Tiempo secuencial en segundos %f \n", stime);

  //Verifica el resultado
  int check = 1;
  for (i = 0; i < N; i++)
  {
    for (j = 0; j < N; j++)
    {
      if (fabs(Cp[BY_ROW(i, j, N)] - Cs[BY_ROW(i, j, N)]) > 0.000001)
      {
        printf("[%d,%d][%f-%f]\n", i, j, Cp[BY_ROW(i, j, N)], Cs[BY_ROW(i, j, N)]);
      }
      check = check && (fabs(Cp[BY_ROW(i, j, N)] - Cs[BY_ROW(i, j, N)]) < 0.000001);
    }
  }

  if (check)
  {
    printf("Calculo de matrices resultado correcto\n");
  }
  else
  {
    printf("Calculo de matrices resultado erroneo\n");
  }

  printf("Speedup = %f \n", stime / ptime);

  printf("Eficiencia = %f \n", (stime / ptime) / Th);

  // Limpieza
  free(A);
  free(B);
  free(Cp);
  free(Cs);
  free(T);
  free(R1);
  free(R2);
  free(M);
  free(r1a);
  free(r2b);
  free(threads);
  pthread_mutex_destroy(&mutex);
  pthread_barrier_destroy(&barrier_avgr);
  pthread_barrier_destroy(&barrier_c);
}

void *worker(void *ptr)
{
  int *p, id;
  p = (int *)ptr;
  id = *p;

  // Calcular AVGR1 y 𝑅1(𝑖,𝑗) = (1 − 𝑇(𝑖,𝑗))(1 − 𝑐𝑜𝑠𝜃(𝑖,𝑗)) + 𝑇(𝑖,𝑗) 𝑠𝑖𝑛𝜃(𝑖,𝑗)
  // y AVGR2y 𝑅2(𝑖,𝑗) = (1 − 𝑇(𝑖,𝑗))(1 − sin𝜃(𝑖,𝑗)) + 𝑇(𝑖,𝑗) cos𝜃(𝑖,𝑗)
  double lavgR1 = 0;
  double lavgR2 = 0;
  for (int i = 0; i < (N / Th); i++)
  {
    for (int j = 0; j < N; j++)
    {
      R1[BY_ROW(id + (Th * i), j, N)] = (1 - T[BY_ROW(id + (Th * i), j, N)]) * (1 - cos(M[BY_ROW(id + (Th * i), j, N)])) + T[BY_ROW(id + (Th * i), j, N)] * sin(M[BY_ROW(id + (Th * i), j, N)]);
      R2[BY_ROW(id + (Th * i), j, N)] = (1 - T[BY_ROW(id + (Th * i), j, N)]) * (1 - sin(M[BY_ROW(id + (Th * i), j, N)])) + T[BY_ROW(id + (Th * i), j, N)] * cos(M[BY_ROW(id + (Th * i), j, N)]);
      lavgR1 += R1[BY_ROW(id + (Th * i), j, N)];
      lavgR2 += R2[BY_ROW(id + (Th * i), j, N)];
    }
  }
  pthread_mutex_lock(&mutex);
  avgR1 += lavgR1;
  avgR2 += lavgR2;
  pthread_mutex_unlock(&mutex);

  pthread_barrier_wait(&barrier_avgr);

  if (id == 0)
  {
    pthread_mutex_lock(&mutex);
    avgR1 /= (N * N);
    avgR2 /= (N * N);
    pthread_mutex_unlock(&mutex);
  }

  // Multiplicación R1*A
  for (int i = 0; i < (N / Th); i += BS)
  {
    for (int j = 0; j < N; j += BS)
    {
      r1a[BY_ROW(id + (Th * i), j, N)] = 0;
      for (int k = 0; k < N; k += BS)
      {
        double *m1 = &R1[BY_ROW(id + (Th * i), k, N)];
        double *m2 = &A[BY_COL(k, j, N)];
        double *mr = &r1a[BY_ROW(id + (Th * i), k, N)];
        for (int x = 0; x < BS; x++)
        {
          for (int y = 0; y < BS; y++)
          {
            for (int z = 0; z < BS; z++)
            {
              mr[BY_ROW(x, y, N)] += m1[BY_ROW(x, z, N)] * m2[BY_COL(z, y, N)];
              printf("R1-%d,%d-%f\n", y, z, mr[BY_ROW(x, y, N)]);
            }
          }
        }
      }
    }
  }

  // Multiplicación R2*B
  for (int i = 0; i < (N / Th); i += BS)
  {
    for (int j = 0; j < N; j += BS)
    {
      r2b[BY_ROW(id + (Th * i), j, N)] = 0;
      for (int k = 0; k < N; k += BS)
      {
        double *m1 = &R2[BY_ROW(id + (Th * i), k, N)];
        double *m2 = &B[BY_COL(k, j, N)];
        double *mr = &r2b[BY_ROW(id + (Th * i), k, N)];
        for (int x = 0; x < BS; x++)
        {
          for (int y = 0; y < BS; y++)
          {
            for (int z = 0; z < BS; z++)
            {
              mr[BY_ROW(x, y, N)] += m1[BY_ROW(x, z, N)] * m2[BY_COL(z, y, N)];
              printf("R2-%d,%d-%f\n", y, z, mr[BY_ROW(x, y, N)]);
            }
          }
        }
      }
    }
  }

  pthread_barrier_wait(&barrier_c);

  // Calcular C = T + avgR1 * avgR1 * (R1 * A + R2 * B)
  for (int i = 0; i < (N / Th); i++)
  {
    for (int j = 0; j < N; j++)
    {
      Cp[BY_ROW(id + (Th * i), j, N)] = T[BY_ROW(id + (Th * i), j, N)] + avgR1 * avgR2 * (r1a[BY_ROW(id + (Th * i), j, N)] + r2b[BY_ROW(id + (Th * i), j, N)]);
    }
  }
  pthread_exit(0);
}

void *calculo_secuencial()
{
  // Calcular AVGR1 y 𝑅1(𝑖,𝑗) = (1 − 𝑇(𝑖,𝑗))(1 − 𝑐𝑜𝑠𝜃(𝑖,𝑗)) + 𝑇(𝑖,𝑗) 𝑠𝑖𝑛𝜃(𝑖,𝑗)
  // y AVGR2y 𝑅2(𝑖,𝑗) = (1 − 𝑇(𝑖,𝑗))(1 − sin𝜃(𝑖,𝑗)) + 𝑇(𝑖,𝑗) cos𝜃(𝑖,𝑗)
  double avgR1 = 0;
  double avgR2 = 0;
  for (int i = 0; i < N; i++)
  {
    for (int j = 0; j < N; j++)
    {
      R1[BY_ROW(i, j, N)] = (1 - T[BY_ROW(i, j, N)]) * (1 - cos(M[BY_ROW(i, j, N)])) + T[BY_ROW(i, j, N)] * sin(M[BY_ROW(i, j, N)]);
      R2[BY_ROW(i, j, N)] = (1 - T[BY_ROW(i, j, N)]) * (1 - sin(M[BY_ROW(i, j, N)])) + T[BY_ROW(i, j, N)] * cos(M[BY_ROW(i, j, N)]);
      avgR1 += R1[BY_ROW(i, j, N)];
      avgR2 += R2[BY_ROW(i, j, N)];
    }
  }
  avgR1 /= (N * N);
  avgR2 /= (N * N);

  // Multiplicación R1*A
  for (int i = 0; i < N; i += BS)
  {
    for (int j = 0; j < N; j += BS)
    {
      r1a[BY_ROW(i, j, N)] = 0;
      for (int k = 0; k < N; k += BS)
      {
        double *m1 = &R1[BY_ROW(i, k, N)];
        double *m2 = &A[BY_COL(k, j, N)];
        double *mr = &r1a[BY_ROW(i, k, N)];
        for (int x = 0; x < BS; x++)
        {
          for (int y = 0; y < BS; y++)
          {
            for (int z = 0; z < BS; z++)
            {
              mr[BY_ROW(x, y, N)] += m1[BY_ROW(x, z, N)] * m2[BY_COL(z, y, N)];
              printf("R1-%d,%d-%f\n", y, z, mr[BY_ROW(x, y, N)]);
            }
          }
        }
      }
    }
  }

  // Multiplicación R2*B
  for (int i = 0; i < N; i += BS)
  {
    for (int j = 0; j < N; j += BS)
    {
      r2b[BY_ROW(i, j, N)] = 0;
      for (int k = 0; k < N; k += BS)
      {
        double *m1 = &R2[BY_ROW(i, k, N)];
        double *m2 = &B[BY_COL(k, j, N)];
        double *mr = &r2b[BY_ROW(i, k, N)];
        for (int x = 0; x < BS; x++)
        {
          for (int y = 0; y < BS; y++)
          {
            for (int z = 0; z < BS; z++)
            {
              mr[BY_ROW(x, y, N)] += m1[BY_ROW(x, z, N)] * m2[BY_COL(z, y, N)];
              printf("R2-%d,%d-%f\n", y, z, mr[BY_ROW(x, y, N)]);
            }
          }
        }
      }
    }
  }

  // Calcular C = T + avgR1 * avgR1 * (R1 * A + R2 * B)
  for (int i = 0; i < N; i++)
  {
    for (int j = 0; j < N; j++)
    {
      Cs[BY_ROW(i, j, N)] = T[BY_ROW(i, j, N)] + avgR1 * avgR2 * (r1a[BY_ROW(i, j, N)] + r2b[BY_ROW(i, j, N)]);
    }
  }
}

#include <sys/time.h>

double dwalltime()
{
  double sec;
  struct timeval tv;

  gettimeofday(&tv, NULL);
  sec = tv.tv_sec + tv.tv_usec / 1000000.0;
  return sec;
}