#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <omp.h>

#define PI 3.14159265358979323846

#define BY_ROW(row, col, len) ((row) * (len) + (col))
#define BY_COL(row, col, len) ((row) + (col) * (len))

/* Time in seconds from some point in the past */
double dwalltime();

double randFP(double min, double max)
{
  double range = (max - min);
  double div = RAND_MAX / range;
  return min + (rand() / div);
}

int main(int argc, char *argv[])
{
  int N, BS, Th;

  // Controla los argumentos al programa
  if (argc != 4 || (N = atoi(argv[1])) <= 0 || (BS = atoi(argv[2])) <= 0 || (Th = atoi(argv[3])) <= 0 || (Th != 1 && Th % 2 != 0) || Th > 8 || (N % (BS * Th) != 0))
  {
    printf("\nError, modo de uso: %s N BS T (N debe ser multiplo de BS x T y T debe ser multiplo de 2 y menor o igual a 8)\n", argv[0]);
    printf("Parametros: Matriz de %d x %d, con bloques de %d y %d threads\n", N, N, BS, Th);

    return 0;
  }
  printf("Calculando algoritmo paralelo empleando OpenMP. Matriz de %d x %d, con bloques de %d y %d threads\n", N, N, BS, Th);

  omp_set_num_threads(Th);

  // Reservar memoria para las matrices y auxiliares
  double *A = (double *)malloc(sizeof(double) * N * N);
  double *B = (double *)malloc(sizeof(double) * N * N);
  double *Cp = (double *)malloc(sizeof(double) * N * N);
  double *Cs = (double *)malloc(sizeof(double) * N * N);
  double *T = (double *)malloc(sizeof(double) * N * N);
  double *R1 = (double *)malloc(sizeof(double) * N * N);
  double *R2 = (double *)malloc(sizeof(double) * N * N);
  double *M = (double *)malloc(sizeof(double) * N * N);
  double *r1ap = (double *)malloc(sizeof(double) * N * N);
  double *r2bp = (double *)malloc(sizeof(double) * N * N);
  double *r1as = (double *)malloc(sizeof(double) * N * N);
  double *r2bs = (double *)malloc(sizeof(double) * N * N);

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
  double avgR1 = 0;
  double avgR2 = 0;
  //CALCULO PARALELO
  double timetick = dwalltime();

#pragma omp parallel shared(R1, R2, A, B, T, M, r1ap, r2bp, Cp, N, BS)
  {
  // Calcular AVGR1 y ????1(????,????) = (1 ??? ????(????,????))(1 ??? ????????????????(????,????)) + ????(????,????) ????????????????(????,????)
  // y AVGR2y ????2(????,????) = (1 ??? ????(????,????))(1 ??? sin????(????,????)) + ????(????,????) cos????(????,????)

#pragma omp for reduction(+:avgR1) schedule(static) reduction(+:avgR2)
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

#pragma omp single
    {
      avgR1 /= N * N;
      avgR2 /= N * N;
    }

    // Multiplicaci??n R1*A
#pragma omp for schedule(static) nowait
    for (int i = 0; i < N; i += BS)
    {
      for (int j = 0; j < N; j += BS)
      {
        r1ap[BY_ROW(i, j, N)] = 0;
        for (int k = 0; k < N; k += BS)
        {
          double *m1 = &R1[BY_ROW(i, k, N)];
          double *m2 = &A[BY_COL(k, j, N)];
          double *mr = &r1ap[BY_ROW(i, j, N)];
          for (int x = 0; x < BS; x++)
          {
            for (int y = 0; y < BS; y++)
            {
              for (int z = 0; z < BS; z++)
              {
                mr[BY_ROW(x, y, N)] += m1[BY_ROW(x, z, N)] * m2[BY_COL(z, y, N)];
              }
            }
          }
        }
      }
    }

    // Multiplicaci??n R2*B
#pragma omp for schedule(static)
    for (int i = 0; i < N; i += BS)
    {
      for (int j = 0; j < N; j += BS)
      {
        r2bp[BY_ROW(i, j, N)] = 0;
        for (int k = 0; k < N; k += BS)
        {
          double *m1 = &R2[BY_ROW(i, k, N)];
          double *m2 = &B[BY_COL(k, j, N)];
          double *mr = &r2bp[BY_ROW(i, j, N)];
          for (int x = 0; x < BS; x++)
          {
            for (int y = 0; y < BS; y++)
            {
              for (int z = 0; z < BS; z++)
              {
                mr[BY_ROW(x, y, N)] += m1[BY_ROW(x, z, N)] * m2[BY_COL(z, y, N)];
              }
            }
          }
        }
      }
    }

    // Calcular C = T + avgR1 * avgR2 * (R1 * A + R2 * B)
#pragma omp for schedule(static)
    for (int i = 0; i < N; i++)
    {
      for (int j = 0; j < N; j++)
      {
        Cp[BY_ROW(i, j, N)] = T[BY_ROW(i, j, N)] + avgR1 * avgR2 * (r1ap[BY_ROW(i, j, N)] + r2bp[BY_ROW(i, j, N)]);
      }
    }
  }

  double ptime = dwalltime() - timetick;
  printf("Tiempo Paralelo OpenMP en segundos %f \n", ptime);

  //CALCULO SECUENCIAL
  timetick = dwalltime();

  // Calcular AVGR1 y ????1(????,????) = (1 ??? ????(????,????))(1 ??? ????????????????(????,????)) + ????(????,????) ????????????????(????,????)
  // y AVGR2y ????2(????,????) = (1 ??? ????(????,????))(1 ??? sin????(????,????)) + ????(????,????) cos????(????,????)
  avgR1 = 0;
  avgR2 = 0;
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

  avgR1 /= N * N;
  avgR2 /= N * N;

  // Multiplicaci??n R1*A
  for (int i = 0; i < N; i += BS)
  {
    for (int j = 0; j < N; j += BS)
    {
      r1as[BY_ROW(i, j, N)] = 0;
      for (int k = 0; k < N; k += BS)
      {
        double *m1 = &R1[BY_ROW(i, k, N)];
        double *m2 = &A[BY_COL(k, j, N)];
        double *mr = &r1as[BY_ROW(i, j, N)];
        for (int x = 0; x < BS; x++)
        {
          for (int y = 0; y < BS; y++)
          {
            for (int z = 0; z < BS; z++)
            {
              mr[BY_ROW(x, y, N)] += m1[BY_ROW(x, z, N)] * m2[BY_COL(z, y, N)];
            }
          }
        }
      }
    }
  }

  // Multiplicaci??n R2*B
  for (int i = 0; i < N; i += BS)
  {
    for (int j = 0; j < N; j += BS)
    {
      r2bs[BY_ROW(i, j, N)] = 0;
      for (int k = 0; k < N; k += BS)
      {
        double *m1 = &R2[BY_ROW(i, k, N)];
        double *m2 = &B[BY_COL(k, j, N)];
        double *mr = &r2bs[BY_ROW(i, j, N)];
        for (int x = 0; x < BS; x++)
        {
          for (int y = 0; y < BS; y++)
          {
            for (int z = 0; z < BS; z++)
            {
              mr[BY_ROW(x, y, N)] += m1[BY_ROW(x, z, N)] * m2[BY_COL(z, y, N)];
            }
          }
        }
      }
    }
  }

  // Calcular C = T + avgR1 * avgR2 * (R1 * A + R2 * B)
  for (int i = 0; i < N; i++)
  {
    for (int j = 0; j < N; j++)
    {
      Cs[BY_ROW(i, j, N)] = T[BY_ROW(i, j, N)] + avgR1 * avgR2 * (r1as[BY_ROW(i, j, N)] + r2bs[BY_ROW(i, j, N)]);
    }
  }

  double stime = dwalltime() - timetick;

  printf("Tiempo secuencial en segundos %f \n", stime);

  //Valido salidas
  int check = 1;
  //Verifica el resultado
  for (int i = 0; i < N; i++)
  {
    for (int j = 0; j < N; j++)
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
  free(r1ap);
  free(r2bp);
  free(r1as);
  free(r2bs);
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
