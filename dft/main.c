#include <stdlib.h>
#include <stdio.h>
#include <stdint.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <fftw3.h>

#define FFT_SIZE 44100

double coss[FFT_SIZE];

void init(void)
{
  for (size_t i = 0; i < FFT_SIZE; i++) {
    double x = (double)i / FFT_SIZE;
    coss[i] = cos(M_PI * 2 * x); // exchange cos here with e.g. triangle will get DFT by triangle instead of sin/cos :)
  }
}

void dft_sse2(double * output, double * input, double * costab, uint32_t fft_size);

void dft(double * output, double * input, double * costab, uint32_t fft_size)
{
  memset(output, 0, sizeof(double) * fft_size);
  uint32_t fft_hsize = fft_size / 2;

  for (size_t n = 0; n < fft_size; n += 4) {
    uint32_t idxs[8] = {0, n, 2 * n, 3 * n, 0, 0, 0, 0};
    idxs[4] = idxs[0] + fft_size * 1 / 4;
    idxs[5] = idxs[1] + fft_size * 1 / 4;
    idxs[6] = idxs[2] + fft_size * 1 / 4;
    idxs[7] = idxs[3] + fft_size * 1 / 4;

    double in[4] = {input[n + 0], input[n + 1], input[n + 2], input[n + 3]};

    uint32_t iidxs[4] = { idxs[0], idxs[1], idxs[2], idxs[3] };
    for (uint32_t k = 0; k < fft_hsize; k += 4) {
      if (idxs[0] >= fft_size) idxs[0] -= fft_size;
      if (idxs[1] >= fft_size) idxs[1] -= fft_size;
      if (idxs[2] >= fft_size) idxs[2] -= fft_size;
      if (idxs[3] >= fft_size) idxs[3] -= fft_size;
      if (idxs[4] >= fft_size) idxs[4] -= fft_size;
      if (idxs[5] >= fft_size) idxs[5] -= fft_size;
      if (idxs[6] >= fft_size) idxs[6] -= fft_size;
      if (idxs[7] >= fft_size) idxs[7] -= fft_size;

      if (idxs[0] >= fft_size) idxs[0] -= fft_size;
      if (idxs[1] >= fft_size) idxs[1] -= fft_size;
      if (idxs[2] >= fft_size) idxs[2] -= fft_size;
      if (idxs[3] >= fft_size) idxs[3] -= fft_size;
      if (idxs[4] >= fft_size) idxs[4] -= fft_size;
      if (idxs[5] >= fft_size) idxs[5] -= fft_size;
      if (idxs[6] >= fft_size) idxs[6] -= fft_size;
      if (idxs[7] >= fft_size) idxs[7] -= fft_size;

      if (idxs[0] >= fft_size) idxs[0] -= fft_size;
      if (idxs[1] >= fft_size) idxs[1] -= fft_size;
      if (idxs[2] >= fft_size) idxs[2] -= fft_size;
      if (idxs[3] >= fft_size) idxs[3] -= fft_size;
      if (idxs[4] >= fft_size) idxs[4] -= fft_size;
      if (idxs[5] >= fft_size) idxs[5] -= fft_size;
      if (idxs[6] >= fft_size) idxs[6] -= fft_size;
      if (idxs[7] >= fft_size) idxs[7] -= fft_size;

      if (idxs[0] >= fft_size) idxs[0] -= fft_size;
      if (idxs[1] >= fft_size) idxs[1] -= fft_size;
      if (idxs[2] >= fft_size) idxs[2] -= fft_size;
      if (idxs[3] >= fft_size) idxs[3] -= fft_size;
      if (idxs[4] >= fft_size) idxs[4] -= fft_size;
      if (idxs[5] >= fft_size) idxs[5] -= fft_size;
      if (idxs[6] >= fft_size) idxs[6] -= fft_size;
      if (idxs[7] >= fft_size) idxs[7] -= fft_size;

      double cs[8] = {
        costab[idxs[0]], costab[idxs[1]], costab[idxs[2]], costab[idxs[3]],
        costab[idxs[4]], costab[idxs[5]], costab[idxs[6]], costab[idxs[7]]};

      output[k + 0] += cs[0] * in[0];
      output[k + 1] += cs[1] * in[1];
      output[k + 2] += cs[2] * in[2];
      output[k + 3] += cs[3] * in[3];
      output[k + 0 + fft_hsize] += cs[4] * in[0];
      output[k + 1 + fft_hsize] += cs[5] * in[1];
      output[k + 2 + fft_hsize] += cs[6] * in[2];
      output[k + 3 + fft_hsize] += cs[7] * in[3];

      idxs[0] += iidxs[0];
      idxs[1] += iidxs[1];
      idxs[2] += iidxs[2];
      idxs[3] += iidxs[3];
      idxs[4] += iidxs[0];
      idxs[5] += iidxs[1];
      idxs[6] += iidxs[2];
      idxs[7] += iidxs[3];
    }
  }
}

void print_ppm(const char * output, double * input, size_t width, size_t height,
    size_t scale)
{
  size_t w = width / scale;
  double * tab = malloc(sizeof(double) * w);
  for (size_t x = 0; x < w; x++) {
    double a = 0.;
    for (size_t z = 0; z < scale; z++)
      a += input[x * scale + z];
    a /= scale;
    tab[x] = a;
  }

  double max = input[0], min = input[0];
  for (size_t i = 0; i < w; i++) {
    max = (max < tab[i]) ? tab[i] : max;
    min = (min > tab[i]) ? tab[i] : min;
  }

  for (size_t i = 0; i < w; i++) {
    tab[i] = (tab[i] - min) / (max - min);
  }

  FILE * f = fopen(output, "w");
  fprintf(f, "P1\n%zu %zu\n", w, height);

  for (size_t y = 0; y < height; y++) {
    double l = (double)(height - y) / (double)height;

    for (size_t x = 0; x < w; x++) {
      fprintf(f, "%u%c", tab[x] <= l, (x + 1 == width) ? '\n' : ' ');
    }
  }

  fclose(f);
}

uint64_t now(void)
{
  struct timespec ts;
  clock_gettime(CLOCK_MONOTONIC, &ts);
  return (uint64_t)ts.tv_sec * 1000000000 + ts.tv_nsec;
}


int main(int argc, char * argv[])
{
  init();

  double * in = fftw_malloc(sizeof(double) * FFT_SIZE);
  double * out = malloc(sizeof(double) * FFT_SIZE);
  double * out2 = malloc(sizeof(double) * FFT_SIZE);
  fftw_complex * out3 = fftw_malloc(sizeof(fftw_complex) * FFT_SIZE);
  fftw_plan p = fftw_plan_dft_r2c_1d(FFT_SIZE, in, out3, FFTW_ESTIMATE);

  for (size_t i = 0; i < FFT_SIZE; i++) {
    in[i] = sin((double)i / 500 * (M_PI * 2.0));
    in[i] += sin((double)i / 750 * (M_PI * 2.0));
  }

  // print_ppm("out1.ppm", in, FFT_SIZE, 1000, 8);
  dft(out, in, coss, FFT_SIZE);
  dft_sse2(out2, in, coss, FFT_SIZE);
  fftw_execute(p);

  uint64_t t_dft = 0, t_sse = 0, t_fftw = 0;
  for (size_t i = 0; i < 10; i++) {
    uint64_t t1 = now();
    dft(out, in, coss, FFT_SIZE);
    uint64_t t2 = now();
    dft_sse2(out2, in, coss, FFT_SIZE);
    uint64_t t3 = now();
    fftw_execute(p);
    uint64_t t4 = now();
    t_dft += t2 - t1;
    t_sse += t3 - t2;
    t_fftw += t4 - t3;
  }

  printf("%llu\n%llu\n%llu\n", t_dft, t_sse, t_fftw);

  //print_ppm("out2.ppm", out, FFT_SIZE, 1000, 8);
  for (size_t i = 0; i < FFT_SIZE; i++) {
    if (out[i] != out2[i])
      printf("%zu %f %f %e\n", i, out[i], out2[i], fabs(out[i] - out2[i]));
  }

  fftw_free(in);
  free(out);
  free(out2);
  fftw_free(out3);
  fftw_destroy_plan(p);

  return 0;
}
