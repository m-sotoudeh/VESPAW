#include "vespaw_gpu_backend.h"

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cctype>
#include <cstdint>

#include <cuda_runtime.h>
#include <cufft.h>
#include <cublas_v2.h>
#include <cusolverDn.h>
#include <cuComplex.h>

/* ========================================================================= */
/* Error-handling helpers                                                    */
/* ========================================================================= */

#define VESPAW_CUDA_CHECK(call) do {                                           \
  cudaError_t err__ = (call);                                                  \
  if (err__ != cudaSuccess) {                                                  \
    fprintf(stderr, "VESPAW CUDA error (%s:%d): %s (%d)\n",                    \
            __FILE__, __LINE__, cudaGetErrorString(err__), (int)err__);        \
    fflush(stderr);                                                            \
    exit(1);                                                                   \
  }                                                                            \
} while (0)

#define VESPAW_CUFFT_CHECK(call) do {                                          \
  cufftResult err__ = (call);                                                  \
  if (err__ != CUFFT_SUCCESS) {                                                \
    fprintf(stderr, "VESPAW cuFFT error (%s:%d): code %d\n",                   \
            __FILE__, __LINE__, (int)err__);                                   \
    fflush(stderr);                                                            \
    exit(1);                                                                   \
  }                                                                            \
} while (0)

#define VESPAW_CUBLAS_CHECK(call) do {                                         \
  cublasStatus_t err__ = (call);                                               \
  if (err__ != CUBLAS_STATUS_SUCCESS) {                                        \
    fprintf(stderr, "VESPAW cuBLAS error (%s:%d): code %d\n",                  \
            __FILE__, __LINE__, (int)err__);                                   \
    fflush(stderr);                                                            \
    exit(1);                                                                   \
  }                                                                            \
} while (0)

#define VESPAW_CUSOLVER_CHECK(call) do {                                       \
  cusolverStatus_t err__ = (call);                                             \
  if (err__ != CUSOLVER_STATUS_SUCCESS) {                                      \
    fprintf(stderr, "VESPAW cuSOLVER error (%s:%d): code %d\n",                \
            __FILE__, __LINE__, (int)err__);                                   \
    fflush(stderr);                                                            \
    exit(1);                                                                   \
  }                                                                            \
} while (0)

/* ========================================================================= */
/* Global persistent runtime state, arena, and LRU cuFFT plan cache          */
/* ========================================================================= */

namespace {

enum ArenaSlot {
  SLOT_X = 0,
  SLOT_Y,
  SLOT_Z,
  SLOT_WORK,
  SLOT_G2,
  SLOT_BUCKET,
  SLOT_POT,
  SLOT_IGTOSTRIPE,
  SLOT_ISTRIPETOYZ,
  SLOT_MINUSG,
  SLOT_OCC,
  SLOT_EI2KR,
  SLOT_GVEC,
  SLOT_VTAU,
  SLOT_RHO,
  SLOT_EIG_A,
  SLOT_EIG_B,
  SLOT_EIG_W,
  SLOT_EIG_WORK,
  SLOT_EIG_INFO,
  SLOT_PRO_TMP1,
  SLOT_PRO_TMP2,
  SLOT_INDEX_MISC,
  SLOT_CUFFT_WORK,
  ARENA_SLOT_COUNT
};

struct FftPlanEntry {
  int valid;
  int rank;
  int n0;
  int n1;
  int n2;
  int istride;
  int idist;
  int batch;
  size_t work_size;
  cufftHandle plan;
  void *bound_work;
  unsigned long long last_used;
};

static const int FFT_PLAN_CACHE_SIZE = 32;
static const int GPU_EIG_CUTOFF = 64;
static const size_t GPU_MAT_THRESHOLD = 32768ULL;

static int              g_initialized = 0;
static int              g_gpu_enabled = 0;
static int              g_dev_id = -1;
static int              g_mpirank = 0;
static cudaStream_t     g_stream = nullptr;
static cublasHandle_t   g_cublas = nullptr;
static cusolverDnHandle_t g_cusolver = nullptr;

static void  *g_arena_ptr[ARENA_SLOT_COUNT] = {nullptr};
static size_t g_arena_cap[ARENA_SLOT_COUNT] = {0};

static FftPlanEntry     g_plan_cache[FFT_PLAN_CACHE_SIZE];
static unsigned long long g_plan_clock = 0;

/* Persistent device projector cache (PROCACHEPRO) */
static void  *g_procache_dev = nullptr;
static size_t g_procache_cap = 0;
static int    g_procache_ngl = 0;
static int    g_procache_npro = 0;
static int    g_procache_nat = 0;
static char   g_procache_id[32] = {0};

static void *arena_get(int slot, size_t bytes) {
  if (bytes == 0) bytes = 16;
  if (bytes > g_arena_cap[slot]) {
    if (g_arena_ptr[slot] != nullptr) {
      VESPAW_CUDA_CHECK(cudaStreamSynchronize(g_stream));
      VESPAW_CUDA_CHECK(cudaFree(g_arena_ptr[slot]));
      g_arena_ptr[slot] = nullptr;
      g_arena_cap[slot] = 0;
    }
    VESPAW_CUDA_CHECK(cudaMalloc(&g_arena_ptr[slot], bytes));
    g_arena_cap[slot] = bytes;
  }
  return g_arena_ptr[slot];
}

static int parse_env_disabled() {
  const char *val = getenv("VESPAW_USE_GPU");
  if (!val || !val[0]) return 0;
  if (strcmp(val, "0") == 0) return 1;
  char c0 = (char)tolower((unsigned char)val[0]);
  if (c0 == 'f' || c0 == 'n') return 1;
  if (strcmp(val, "off") == 0 || strcmp(val, "OFF") == 0) return 1;
  return 0;
}

static cufftHandle get_cached_cufft_plan(int rank, int n0, int n1, int n2,
                                         int istride, int idist, int batch) {
  ++g_plan_clock;
  int free_idx = -1;
  int lru_idx = 0;
  unsigned long long min_used = ~0ULL;

  for (int i = 0; i < FFT_PLAN_CACHE_SIZE; ++i) {
    FftPlanEntry &e = g_plan_cache[i];
    if (!e.valid) {
      if (free_idx < 0) free_idx = i;
      continue;
    }
    if (e.rank == rank && e.n0 == n0 && e.n1 == n1 && e.n2 == n2 &&
        e.istride == istride && e.idist == idist && e.batch == batch) {
      e.last_used = g_plan_clock;
      void *work_ptr = arena_get(SLOT_CUFFT_WORK, e.work_size);
      if (e.bound_work != work_ptr) {
        VESPAW_CUFFT_CHECK(cufftSetWorkArea(e.plan, work_ptr));
        e.bound_work = work_ptr;
      }
      return e.plan;
    }
    if (e.last_used < min_used) {
      min_used = e.last_used;
      lru_idx = i;
    }
  }

  int slot = (free_idx >= 0) ? free_idx : lru_idx;
  if (g_plan_cache[slot].valid) {
    VESPAW_CUDA_CHECK(cudaStreamSynchronize(g_stream));
    VESPAW_CUFFT_CHECK(cufftDestroy(g_plan_cache[slot].plan));
    g_plan_cache[slot].valid = 0;
  }

  int dims[3] = {n0, n1, n2};
  cufftHandle plan;
  size_t work_size = 0;
  VESPAW_CUFFT_CHECK(cufftCreate(&plan));
  VESPAW_CUFFT_CHECK(cufftSetAutoAllocation(plan, 0));
  VESPAW_CUFFT_CHECK(cufftMakePlanMany(plan, rank, dims,
                                       nullptr, istride, idist,
                                       nullptr, istride, idist,
                                       CUFFT_Z2Z, batch, &work_size));
  VESPAW_CUFFT_CHECK(cufftSetStream(plan, g_stream));
  void *work_ptr = arena_get(SLOT_CUFFT_WORK, work_size);
  VESPAW_CUFFT_CHECK(cufftSetWorkArea(plan, work_ptr));

  FftPlanEntry &ne = g_plan_cache[slot];
  ne.valid = 1;
  ne.rank = rank;
  ne.n0 = n0;
  ne.n1 = n1;
  ne.n2 = n2;
  ne.istride = istride;
  ne.idist = idist;
  ne.batch = batch;
  ne.work_size = work_size;
  ne.plan = plan;
  ne.bound_work = work_ptr;
  ne.last_used = g_plan_clock;
  return plan;
}

static size_t block_grid_1d(size_t n, int threads) {
  return (n + (size_t)threads - 1) / (size_t)threads;
}

/* ========================================================================= */
/* CUDA Kernels                                                              */
/* ========================================================================= */

__global__ void kernel_scale_c8(cuDoubleComplex *data, size_t n, double scale) {
  size_t idx = (size_t)blockIdx.x * blockDim.x + threadIdx.x;
  if (idx < n) {
    cuDoubleComplex v = data[idx];
    data[idx] = make_cuDoubleComplex(cuCreal(v) * scale, cuCimag(v) * scale);
  }
}

__global__ void kernel_scale_r8(double *data, size_t n, double scale) {
  size_t idx = (size_t)blockIdx.x * blockDim.x + threadIdx.x;
  if (idx < n) {
    data[idx] *= scale;
  }
}

__global__ void kernel_slab_pack(int len, int nslab, int nr1l, int stride_r1,
                                 int nfft, int stride_fft,
                                 const cuDoubleComplex *src, cuDoubleComplex *dst) {
  size_t total = (size_t)len * (size_t)nslab * (size_t)nr1l * (size_t)nfft;
  size_t idx = (size_t)blockIdx.x * blockDim.x + threadIdx.x;
  if (idx >= total) return;

  size_t slab_stride = (size_t)len;
  int k = (int)(idx % slab_stride);
  size_t t1 = idx / slab_stride;
  int islab = (int)(t1 % (size_t)nslab);
  size_t t2 = t1 / (size_t)nslab;
  int ir1l = (int)(t2 % (size_t)nr1l);
  int ifft = (int)(t2 / (size_t)nr1l);

  size_t sidx = (size_t)ifft * (size_t)stride_fft +
                (size_t)ir1l * (size_t)stride_r1 +
                (size_t)islab * (size_t)len + (size_t)k;
  dst[idx] = src[sidx];
}

__global__ void kernel_slab_unpack_scaled(int len, int nslab, int nr1l, int stride_r1,
                                          int nfft, int stride_fft, double scale,
                                          const cuDoubleComplex *src, cuDoubleComplex *dst) {
  size_t total = (size_t)len * (size_t)nslab * (size_t)nr1l * (size_t)nfft;
  size_t idx = (size_t)blockIdx.x * blockDim.x + threadIdx.x;
  if (idx >= total) return;

  size_t slab_stride = (size_t)len;
  int k = (int)(idx % slab_stride);
  size_t t1 = idx / slab_stride;
  int islab = (int)(t1 % (size_t)nslab);
  size_t t2 = t1 / (size_t)nslab;
  int ir1l = (int)(t2 % (size_t)nr1l);
  int ifft = (int)(t2 / (size_t)nr1l);

  size_t sidx = (size_t)ifft * (size_t)stride_fft +
                (size_t)ir1l * (size_t)stride_r1 +
                (size_t)islab * (size_t)len + (size_t)k;
  cuDoubleComplex v = src[idx];
  dst[sidx] = make_cuDoubleComplex(cuCreal(v) * scale, cuCimag(v) * scale);
}

__global__ void kernel_scatter_stripes_to_3d(int ngl, int nr1, int nrl, int nfft,
                                             const int *igtostripe,
                                             const int *istripetoyz,
                                             const cuDoubleComplex *fofg,
                                             cuDoubleComplex *grid3d) {
  int ig = blockIdx.x * blockDim.x + threadIdx.x;
  int ifft = blockIdx.y;
  if (ig >= ngl || ifft >= nfft) return;

  int s = igtostripe[ig] - 1;
  int ir1 = s % nr1;
  int istripel = s / nr1;
  int i23 = istripetoyz[istripel] - 1;
  int idx3d = ir1 + nr1 * i23;
  grid3d[(size_t)idx3d + (size_t)ifft * (size_t)nrl] =
      fofg[(size_t)ig + (size_t)ifft * (size_t)ngl];
}

__global__ void kernel_gather_3d_to_stripes_scaled(int ngl, int nr1, int nrl, int nfft,
                                                   double scale,
                                                   const int *igtostripe,
                                                   const int *istripetoyz,
                                                   const cuDoubleComplex *grid3d,
                                                   cuDoubleComplex *fofg) {
  int ig = blockIdx.x * blockDim.x + threadIdx.x;
  int ifft = blockIdx.y;
  if (ig >= ngl || ifft >= nfft) return;

  int s = igtostripe[ig] - 1;
  int ir1 = s % nr1;
  int istripel = s / nr1;
  int i23 = istripetoyz[istripel] - 1;
  int idx3d = ir1 + nr1 * i23;
  cuDoubleComplex v = grid3d[(size_t)idx3d + (size_t)ifft * (size_t)nrl];
  fofg[(size_t)ig + (size_t)ifft * (size_t)ngl] =
      make_cuDoubleComplex(cuCreal(v) * scale, cuCimag(v) * scale);
}

__global__ void kernel_vpsi_mul_ndim1(int nrl, int nblk, const double *v,
                                      cuDoubleComplex *psiofr) {
  int ir = blockIdx.x * blockDim.x + threadIdx.x;
  int ib = blockIdx.y;
  if (ir >= nrl || ib >= nblk) return;

  double pot = v[ir];
  size_t off = (size_t)ir + (size_t)ib * (size_t)nrl;
  cuDoubleComplex val = psiofr[off];
  psiofr[off] = make_cuDoubleComplex(pot * cuCreal(val), pot * cuCimag(val));
}

__global__ void kernel_vpsi_mul_ndim2(int nrl, int nblk, const double *v,
                                      cuDoubleComplex *psiofr) {
  int ir = blockIdx.x * blockDim.x + threadIdx.x;
  int ib = blockIdx.y;
  if (ir >= nrl || ib >= nblk) return;

  double v0 = v[ir];
  double v1 = v[ir + (size_t)nrl];
  double v2 = v[ir + 2 * (size_t)nrl];
  double v3 = v[ir + 3 * (size_t)nrl];
  double vupup = v0 + v3;
  double vdndn = v0 - v3;
  cuDoubleComplex vupdn = make_cuDoubleComplex(v1, -v2);

  size_t off_up = (size_t)ir + (size_t)(2 * ib) * (size_t)nrl;
  size_t off_dn = (size_t)ir + (size_t)(2 * ib + 1) * (size_t)nrl;
  cuDoubleComplex pup = psiofr[off_up];
  cuDoubleComplex pdn = psiofr[off_dn];

  psiofr[off_up] = cuCadd(make_cuDoubleComplex(vupup * cuCreal(pup), vupup * cuCimag(pup)),
                          cuCmul(vupdn, pdn));
  psiofr[off_dn] = cuCadd(make_cuDoubleComplex(vdndn * cuCreal(pdn), vdndn * cuCimag(pdn)),
                          cuCmul(cuConj(vupdn), pup));
}

__global__ void kernel_gather_3d_vpsi_kin(int ngl, int nr1, int nrl, int nblk_fft,
                                          int tbucket, double scale,
                                          const int *igtostripe,
                                          const int *istripetoyz,
                                          const double *g2,
                                          const double *bucket,
                                          const cuDoubleComplex *grid3d,
                                          const cuDoubleComplex *psi,
                                          cuDoubleComplex *hpsi) {
  int ig = blockIdx.x * blockDim.x + threadIdx.x;
  int ifft = blockIdx.y;
  if (ig >= ngl || ifft >= nblk_fft) return;

  int s = igtostripe[ig] - 1;
  int ir1 = s % nr1;
  int istripel = s / nr1;
  int i23 = istripetoyz[istripel] - 1;
  int idx3d = ir1 + nr1 * i23;

  cuDoubleComplex vpot = grid3d[(size_t)idx3d + (size_t)ifft * (size_t)nrl];
  cuDoubleComplex p = psi[(size_t)ig + (size_t)ifft * (size_t)ngl];
  double kin = 0.5 * g2[ig] + (tbucket ? bucket[ig] : 0.0);

  hpsi[(size_t)ig + (size_t)ifft * (size_t)ngl] =
      make_cuDoubleComplex(cuCreal(vpot) * scale + kin * cuCreal(p),
                           cuCimag(vpot) * scale + kin * cuCimag(p));
}

__global__ void kernel_scatter_tsuper_mapped(int nhalf, int nglh, int nr1, int nrl,
                                             int npairs, int nb,
                                             const int *ig1_list, const int *ig2_list,
                                             const int *idx3d_map,
                                             const cuDoubleComplex *psi,
                                             cuDoubleComplex *grid3d) {
  int igh = blockIdx.x * blockDim.x + threadIdx.x;
  int ipair = blockIdx.y;
  if (igh >= nhalf || ipair >= npairs) return;

  int ib1 = 2 * ipair;
  int ib2 = ib1 + 1;
  cuDoubleComplex a = psi[(size_t)igh + (size_t)ib1 * (size_t)nglh];
  cuDoubleComplex b = (ib2 < nb) ? psi[(size_t)igh + (size_t)ib2 * (size_t)nglh]
                                 : make_cuDoubleComplex(0.0, 0.0);

  cuDoubleComplex f1 = make_cuDoubleComplex(cuCreal(a) - cuCimag(b),
                                            cuCimag(a) + cuCreal(b));
  cuDoubleComplex f2 = make_cuDoubleComplex(cuCreal(a) + cuCimag(b),
                                            -cuCimag(a) + cuCreal(b));

  int i1 = idx3d_map[ig1_list[igh]];
  int i2 = idx3d_map[ig2_list[igh]];
  size_t base = (size_t)ipair * (size_t)nrl;
  grid3d[base + (size_t)i1] = f1;
  if (i2 != i1) {
    grid3d[base + (size_t)i2] = f2;
  }
}

__global__ void kernel_accum_rho_tsuper(int nrl, int npairs, int nb,
                                        const double *occ,
                                        const cuDoubleComplex *psiofr,
                                        double *rho) {
  int ir = blockIdx.x * blockDim.x + threadIdx.x;
  if (ir >= nrl) return;
  double sum = rho[ir];
  for (int ipair = 0; ipair < npairs; ++ipair) {
    cuDoubleComplex z = psiofr[(size_t)ir + (size_t)ipair * (size_t)nrl];
    double re = cuCreal(z);
    double im = cuCimag(z);
    int ib1 = 2 * ipair;
    int ib2 = ib1 + 1;
    double f1 = occ[ib1];
    if (f1 != 0.0) sum += f1 * re * re;
    if (ib2 < nb) {
      double f2 = occ[ib2];
      if (f2 != 0.0) sum += f2 * im * im;
    }
  }
  rho[ir] = sum;
}

__global__ void kernel_accum_rho_tinv(int nrl, int nblk, int ib0, int nb,
                                      const double *occ,
                                      const cuDoubleComplex *ei2kr,
                                      const cuDoubleComplex *psiofr,
                                      double *rho) {
  int ir = blockIdx.x * blockDim.x + threadIdx.x;
  if (ir >= nrl) return;
  double sum = rho[ir];
  cuDoubleComplex phase2 = ei2kr[ir];
  for (int k = 0; k < nblk; ++k) {
    int ibh = ib0 + k;
    double f1 = occ[2 * ibh];
    double f2 = (2 * ibh + 1 < nb) ? occ[2 * ibh + 1] : 0.0;
    double s1 = 0.5 * (f1 + f2);
    double s2 = 0.5 * (f1 - f2);
    if (s1 == 0.0 && s2 == 0.0) continue;
    cuDoubleComplex p = psiofr[(size_t)ir + (size_t)k * (size_t)nrl];
    double mod2 = cuCreal(p) * cuCreal(p) + cuCimag(p) * cuCimag(p);
    double term = s1 * mod2;
    if (s2 != 0.0) {
      cuDoubleComplex p2 = cuCmul(p, p);
      cuDoubleComplex p2e = cuCmul(p2, phase2);
      term += s2 * cuCreal(p2e);
    }
    sum += term;
  }
  rho[ir] = sum;
}

__global__ void kernel_accum_rho_ndim1(int nrl, int nblk, int ib0,
                                       const double *occ,
                                       const cuDoubleComplex *psiofr,
                                       double *rho) {
  int ir = blockIdx.x * blockDim.x + threadIdx.x;
  if (ir >= nrl) return;
  double sum = rho[ir];
  for (int k = 0; k < nblk; ++k) {
    double f1 = occ[ib0 + k];
    if (f1 == 0.0) continue;
    cuDoubleComplex p = psiofr[(size_t)ir + (size_t)k * (size_t)nrl];
    sum += f1 * (cuCreal(p) * cuCreal(p) + cuCimag(p) * cuCimag(p));
  }
  rho[ir] = sum;
}

__global__ void kernel_accum_rho_ndim2(int nrl, int nblk, int ib0,
                                       const double *occ,
                                       const cuDoubleComplex *psiofr,
                                       double *rho) {
  int ir = blockIdx.x * blockDim.x + threadIdx.x;
  if (ir >= nrl) return;
  double r0 = rho[ir];
  double r1 = rho[ir + (size_t)nrl];
  double r2 = rho[ir + 2 * (size_t)nrl];
  double r3 = rho[ir + 3 * (size_t)nrl];
  for (int k = 0; k < nblk; ++k) {
    double f1 = occ[ib0 + k];
    if (f1 == 0.0) continue;
    cuDoubleComplex pu = psiofr[(size_t)ir + (size_t)(2 * k) * (size_t)nrl];
    cuDoubleComplex pd = psiofr[(size_t)ir + (size_t)(2 * k + 1) * (size_t)nrl];
    r0 += f1 * (cuCreal(pu) * cuCreal(pu) + cuCimag(pu) * cuCimag(pu));
    r3 += f1 * (cuCreal(pd) * cuCreal(pd) + cuCimag(pd) * cuCimag(pd));
    cuDoubleComplex c = cuCmul(pu, cuConj(pd));
    r1 += f1 * cuCreal(c);
    r2 += f1 * cuCimag(c);
  }
  rho[ir] = r0;
  rho[ir + (size_t)nrl] = r1;
  rho[ir + 2 * (size_t)nrl] = r2;
  rho[ir + 3 * (size_t)nrl] = r3;
}

__global__ void kernel_mul_ig_grad(int ngl, int nbh, int idir,
                                   const double *gvec,
                                   const cuDoubleComplex *psi,
                                   cuDoubleComplex *psi_d) {
  int ig = blockIdx.x * blockDim.x + threadIdx.x;
  int ib = blockIdx.y;
  if (ig >= ngl || ib >= nbh) return;
  double gd = gvec[(size_t)idir + 3 * (size_t)ig];
  cuDoubleComplex z = psi[(size_t)ig + (size_t)ib * (size_t)ngl];
  /* i * gd * (re + i*im) = (-gd*im) + i*(gd*re) */
  psi_d[(size_t)ig + (size_t)ib * (size_t)ngl] =
      make_cuDoubleComplex(-gd * cuCimag(z), gd * cuCreal(z));
}

__global__ void kernel_accum_ig_htaupsi(int ngl, int nbh, int idir,
                                        const double *gvec,
                                        const cuDoubleComplex *wofg,
                                        cuDoubleComplex *htaupsi) {
  int ig = blockIdx.x * blockDim.x + threadIdx.x;
  int ib = blockIdx.y;
  if (ig >= ngl || ib >= nbh) return;
  double gd = gvec[(size_t)idir + 3 * (size_t)ig];
  size_t off = (size_t)ig + (size_t)ib * (size_t)ngl;
  cuDoubleComplex z = wofg[off];
  cuDoubleComplex h = htaupsi[off];
  htaupsi[off] = make_cuDoubleComplex(cuCreal(h) - gd * cuCimag(z),
                                      cuCimag(h) + gd * cuCreal(z));
}

__global__ void kernel_transpose_c8(int rows, int cols,
                                    const cuDoubleComplex *in,
                                    cuDoubleComplex *out) {
  int r = blockIdx.x * blockDim.x + threadIdx.x;
  int c = blockIdx.y * blockDim.y + threadIdx.y;
  if (r < rows && c < cols) {
    out[(size_t)c + (size_t)r * (size_t)cols] =
        in[(size_t)r + (size_t)c * (size_t)rows];
  }
}

__global__ void kernel_adjust_tsuper_gweight(int n1, int n2, int ngl, int ndim,
                                             int tsuper, int ngamma, double gweight,
                                             const cuDoubleComplex *f1,
                                             const cuDoubleComplex *f2,
                                             cuDoubleComplex *mat) {
  int i1 = blockIdx.x * blockDim.x + threadIdx.x;
  int i2 = blockIdx.y * blockDim.y + threadIdx.y;
  if (i1 >= n1 || i2 >= n2) return;

  cuDoubleComplex m = mat[(size_t)i1 + (size_t)i2 * (size_t)n1];
  if (tsuper) {
    m = make_cuDoubleComplex(2.0 * cuCreal(m), 2.0 * cuCimag(m));
    if (ngamma > 0) {
      int g0 = ngamma - 1;
      for (int idim = 0; idim < ndim; ++idim) {
        cuDoubleComplex a = f1[(size_t)g0 + (size_t)ngl * ((size_t)idim + (size_t)ndim * (size_t)i1)];
        cuDoubleComplex b = f2[(size_t)g0 + (size_t)ngl * ((size_t)idim + (size_t)ndim * (size_t)i2)];
        m = cuCsub(m, cuCmul(cuConj(a), b));
      }
    }
  }
  mat[(size_t)i1 + (size_t)i2 * (size_t)n1] =
      make_cuDoubleComplex(cuCreal(m) * gweight, cuCimag(m) * gweight);
}

__global__ void kernel_dotdiag_c8(int ngl, int ndim, int n, int tsuper,
                                  int ngamma, double gweight,
                                  const cuDoubleComplex *f1,
                                  const cuDoubleComplex *f2,
                                  double *diag) {
  __shared__ double sdata[256];
  int k = blockIdx.x;
  if (k >= n) return;
  int tid = threadIdx.x;
  size_t vec_len = (size_t)ngl * (size_t)ndim;
  size_t base = (size_t)k * vec_len;

  double sum = 0.0;
  for (size_t i = (size_t)tid; i < vec_len; i += (size_t)blockDim.x) {
    cuDoubleComplex a = f1[base + i];
    cuDoubleComplex b = f2[base + i];
    sum += cuCreal(a) * cuCreal(b) + cuCimag(a) * cuCimag(b);
  }
  sdata[tid] = sum;
  __syncthreads();

  for (int stride = blockDim.x / 2; stride > 0; stride >>= 1) {
    if (tid < stride) sdata[tid] += sdata[tid + stride];
    __syncthreads();
  }

  if (tid == 0) {
    double s = sdata[0];
    if (tsuper) {
      s *= 2.0;
      if (ngamma > 0) {
        int g0 = ngamma - 1;
        for (int idim = 0; idim < ndim; ++idim) {
          cuDoubleComplex a = f1[base + (size_t)g0 + (size_t)idim * (size_t)ngl];
          cuDoubleComplex b = f2[base + (size_t)g0 + (size_t)idim * (size_t)ngl];
          s -= (cuCreal(a) * cuCreal(b) + cuCimag(a) * cuCimag(b));
        }
      }
    }
    diag[k] = s * gweight;
  }
}

__global__ void kernel_precondition_teter(int ngl, int ncols,
                                          const double *g2,
                                          const double *ekin,
                                          const cuDoubleComplex *r_in,
                                          cuDoubleComplex *d_out) {
  int ig = blockIdx.x * blockDim.x + threadIdx.x;
  int col = blockIdx.y;
  if (ig >= ngl || col >= ncols) return;

  double ek = ekin[col];
  if (ek < 0.5) ek = 0.5;
  double fact = 0.5 / ek;
  double x = g2[ig] * fact;
  double x2 = x * x;
  double x3 = x2 * x;
  double x4 = x2 * x2;
  double num = 27.0 + 18.0 * x + 12.0 * x2 + 8.0 * x3;
  double f = num / (num + 16.0 * x4);
  size_t off = (size_t)ig + (size_t)col * (size_t)ngl;
  cuDoubleComplex z = r_in[off];
  d_out[off] = make_cuDoubleComplex(f * cuCreal(z), f * cuCimag(z));
}

__global__ void kernel_expandpro_all(int ngl, int nat, int lnxx, int lmx,
                                     const int *isp_of_atom,
                                     const int *atom_ipro0,
                                     const int *sp_ibpro0,
                                     const int *lnx_of_sp,
                                     const int *lox_all,
                                     const double *r_atoms,
                                     const double *gvec,
                                     const double *barepro,
                                     const double *ylm,
                                     cuDoubleComplex *procache) {
  int ig = blockIdx.x * blockDim.x + threadIdx.x;
  int iat = blockIdx.y;
  if (ig >= ngl || iat >= nat) return;

  int isp = isp_of_atom[iat] - 1;
  int lnx = lnx_of_sp[isp];
  int ibpro0 = sp_ibpro0[isp];
  int ipro0 = atom_ipro0[iat];

  double rx = r_atoms[0 + 3 * iat];
  double ry = r_atoms[1 + 3 * iat];
  double rz = r_atoms[2 + 3 * iat];
  double gx = gvec[0 + 3 * ig];
  double gy = gvec[1 + 3 * ig];
  double gz = gvec[2 + 3 * ig];
  double phase = -(gx * rx + gy * ry + gz * rz);
  double cph = cos(phase);
  double sph = sin(phase);

  int lmn = 0;
  for (int ln = 0; ln < lnx; ++ln) {
    int l = lox_all[ln + isp * lnxx];
    int lm_base = l * l;
    double bp = barepro[(size_t)ig + (size_t)(ibpro0 + ln) * (size_t)ngl];

    /* (-i)^L */
    cuDoubleComplex cil;
    switch (l & 3) {
      case 0: cil = make_cuDoubleComplex( 1.0,  0.0); break;
      case 1: cil = make_cuDoubleComplex( 0.0, -1.0); break;
      case 2: cil = make_cuDoubleComplex(-1.0,  0.0); break;
      default: cil = make_cuDoubleComplex( 0.0,  1.0); break;
    }
    cuDoubleComplex ceigr = make_cuDoubleComplex(
        cph * cuCreal(cil) - sph * cuCimag(cil),
        cph * cuCimag(cil) + sph * cuCreal(cil));

    for (int m = 0; m < 2 * l + 1; ++m) {
      int lm = lm_base + m;
      double s = bp * ylm[(size_t)ig + (size_t)lm * (size_t)ngl];
      procache[(size_t)ig + (size_t)(ipro0 + lmn) * (size_t)ngl] =
          make_cuDoubleComplex(s * cuCreal(ceigr), s * cuCimag(ceigr));
      ++lmn;
    }
  }
}

/* Host helpers for device 3-D scatter + batched cuFFT + scale */
static void gpu_fft_3d_device_core(int dir_sign, int nr1, int nr2, int nr3,
                                   int nfft, void *d_grid) {
  cufftHandle plan = get_cached_cufft_plan(3, nr3, nr2, nr1, 1, nr1 * nr2 * nr3, nfft);
  int cufft_dir = (dir_sign < 0) ? CUFFT_FORWARD : CUFFT_INVERSE;
  VESPAW_CUFFT_CHECK(cufftExecZ2Z(plan,
                                  (cufftDoubleComplex *)d_grid,
                                  (cufftDoubleComplex *)d_grid,
                                  cufft_dir));
  if (dir_sign < 0) {
    size_t ntot = (size_t)nr1 * (size_t)nr2 * (size_t)nr3 * (size_t)nfft;
    double inv_n = 1.0 / (double)((size_t)nr1 * (size_t)nr2 * (size_t)nr3);
    int threads = 256;
    kernel_scale_c8<<<(unsigned)block_grid_1d(ntot, threads), threads, 0, g_stream>>>(
        (cuDoubleComplex *)d_grid, ntot, inv_n);
  }
}

} /* namespace */

/* ========================================================================= */
/* Public C ABI Implementation                                               */
/* ========================================================================= */

extern "C" {

void vespaw_gpu_init(int mpirank, int mpisize, int *enabled_out, int *dev_id_out) {
  (void)mpisize;
  if (g_initialized) {
    if (enabled_out) *enabled_out = g_gpu_enabled;
    if (dev_id_out) *dev_id_out = g_dev_id;
    return;
  }
  g_initialized = 1;
  g_mpirank = mpirank;
  memset(g_plan_cache, 0, sizeof(g_plan_cache));

  if (parse_env_disabled()) {
    g_gpu_enabled = 0;
    g_dev_id = -1;
    if (enabled_out) *enabled_out = 0;
    if (dev_id_out) *dev_id_out = -1;
    return;
  }

  int ndev = 0;
  cudaError_t err = cudaGetDeviceCount(&ndev);
  if (err != cudaSuccess || ndev <= 0) {
    g_gpu_enabled = 0;
    g_dev_id = -1;
    if (enabled_out) *enabled_out = 0;
    if (dev_id_out) *dev_id_out = -1;
    return;
  }

  int dev = (mpirank >= 0) ? (mpirank % ndev) : 0;
  const char *env_dev = getenv("VESPAW_GPU_ID");
  if (env_dev && env_dev[0] && mpisize <= 1) {
    int parsed = atoi(env_dev);
    if (parsed >= 0 && parsed < ndev) dev = parsed;
  }

  VESPAW_CUDA_CHECK(cudaSetDevice(dev));
  VESPAW_CUDA_CHECK(cudaStreamCreateWithFlags(&g_stream, cudaStreamNonBlocking));
  VESPAW_CUBLAS_CHECK(cublasCreate(&g_cublas));
  VESPAW_CUBLAS_CHECK(cublasSetStream(g_cublas, g_stream));
  VESPAW_CUSOLVER_CHECK(cusolverDnCreate(&g_cusolver));
  VESPAW_CUSOLVER_CHECK(cusolverDnSetStream(g_cusolver, g_stream));

  g_dev_id = dev;
  g_gpu_enabled = 1;
  if (enabled_out) *enabled_out = 1;
  if (dev_id_out) *dev_id_out = dev;
}

void vespaw_gpu_shutdown(void) {
  if (!g_initialized) return;
  if (g_gpu_enabled) {
    VESPAW_CUDA_CHECK(cudaStreamSynchronize(g_stream));
    for (int i = 0; i < FFT_PLAN_CACHE_SIZE; ++i) {
      if (g_plan_cache[i].valid) {
        VESPAW_CUFFT_CHECK(cufftDestroy(g_plan_cache[i].plan));
        g_plan_cache[i].valid = 0;
      }
    }
    for (int i = 0; i < ARENA_SLOT_COUNT; ++i) {
      if (g_arena_ptr[i]) {
        VESPAW_CUDA_CHECK(cudaFree(g_arena_ptr[i]));
        g_arena_ptr[i] = nullptr;
        g_arena_cap[i] = 0;
      }
    }
    if (g_procache_dev) {
      VESPAW_CUDA_CHECK(cudaFree(g_procache_dev));
      g_procache_dev = nullptr;
      g_procache_cap = 0;
    }
    if (g_cusolver) {
      VESPAW_CUSOLVER_CHECK(cusolverDnDestroy(g_cusolver));
      g_cusolver = nullptr;
    }
    if (g_cublas) {
      VESPAW_CUBLAS_CHECK(cublasDestroy(g_cublas));
      g_cublas = nullptr;
    }
    if (g_stream) {
      VESPAW_CUDA_CHECK(cudaStreamDestroy(g_stream));
      g_stream = nullptr;
    }
  }
  g_gpu_enabled = 0;
  g_dev_id = -1;
  g_initialized = 0;
}

int vespaw_gpu_is_enabled(void) {
  return g_gpu_enabled;
}

void vespaw_gpu_sync(void) {
  if (g_gpu_enabled && g_stream) {
    VESPAW_CUDA_CHECK(cudaStreamSynchronize(g_stream));
  }
}

void vespaw_gpu_fft_1d(int dir_sign, int len, int nfft, const void *h_in, void *h_out) {
  if (!g_gpu_enabled || len <= 0 || nfft <= 0) return;
  size_t ntot = (size_t)len * (size_t)nfft;
  size_t bytes = ntot * sizeof(cuDoubleComplex);
  void *d_x = arena_get(SLOT_X, bytes);
  void *d_y = arena_get(SLOT_Y, bytes);

  VESPAW_CUDA_CHECK(cudaMemcpyAsync(d_x, h_in, bytes, cudaMemcpyHostToDevice, g_stream));
  cufftHandle plan = get_cached_cufft_plan(1, len, 0, 0, 1, len, nfft);
  int cufft_dir = (dir_sign < 0) ? CUFFT_FORWARD : CUFFT_INVERSE;
  VESPAW_CUFFT_CHECK(cufftExecZ2Z(plan, (cufftDoubleComplex *)d_x,
                                  (cufftDoubleComplex *)d_y, cufft_dir));
  if (dir_sign < 0) {
    int threads = 256;
    kernel_scale_c8<<<(unsigned)block_grid_1d(ntot, threads), threads, 0, g_stream>>>(
        (cuDoubleComplex *)d_y, ntot, 1.0 / (double)len);
  }
  VESPAW_CUDA_CHECK(cudaMemcpyAsync(h_out, d_y, bytes, cudaMemcpyDeviceToHost, g_stream));
  VESPAW_CUDA_CHECK(cudaStreamSynchronize(g_stream));
}

void vespaw_gpu_fft_slabs(int dir_sign, int len, int nslab, int nr1l, int stride_r1,
                          int nfft, int stride_fft, int ioff, void *h_inout) {
  if (!g_gpu_enabled || len <= 0 || nslab <= 0 || nr1l <= 0 || nfft <= 0) return;

  size_t span_elems = (size_t)(nfft - 1) * (size_t)stride_fft +
                      (size_t)(nr1l - 1) * (size_t)stride_r1 +
                      (size_t)nslab * (size_t)len;
  size_t span_bytes = span_elems * sizeof(cuDoubleComplex);
  size_t pack_elems = (size_t)len * (size_t)nslab * (size_t)nr1l * (size_t)nfft;
  size_t pack_bytes = pack_elems * sizeof(cuDoubleComplex);

  cuDoubleComplex *h_base = ((cuDoubleComplex *)h_inout) + ioff;
  void *d_span = arena_get(SLOT_X, span_bytes);
  void *d_pack = arena_get(SLOT_Y, pack_bytes);

  VESPAW_CUDA_CHECK(cudaMemcpyAsync(d_span, h_base, span_bytes, cudaMemcpyHostToDevice, g_stream));

  int threads = 256;
  kernel_slab_pack<<<(unsigned)block_grid_1d(pack_elems, threads), threads, 0, g_stream>>>(
      len, nslab, nr1l, stride_r1, nfft, stride_fft,
      (const cuDoubleComplex *)d_span, (cuDoubleComplex *)d_pack);

  cufftHandle plan = get_cached_cufft_plan(1, len, 0, 0, 1, len, nslab * nr1l * nfft);
  int cufft_dir = (dir_sign < 0) ? CUFFT_FORWARD : CUFFT_INVERSE;
  VESPAW_CUFFT_CHECK(cufftExecZ2Z(plan, (cufftDoubleComplex *)d_pack,
                                  (cufftDoubleComplex *)d_pack, cufft_dir));

  double scale = (dir_sign < 0) ? (1.0 / (double)len) : 1.0;
  kernel_slab_unpack_scaled<<<(unsigned)block_grid_1d(pack_elems, threads), threads, 0, g_stream>>>(
      len, nslab, nr1l, stride_r1, nfft, stride_fft, scale,
      (const cuDoubleComplex *)d_pack, (cuDoubleComplex *)d_span);

  VESPAW_CUDA_CHECK(cudaMemcpyAsync(h_base, d_span, span_bytes, cudaMemcpyDeviceToHost, g_stream));
  VESPAW_CUDA_CHECK(cudaStreamSynchronize(g_stream));
}

void vespaw_gpu_fft_3d(int dir_sign, int n1, int n2, int n3, int nfft,
                       const void *h_in, void *h_out) {
  if (!g_gpu_enabled || n1 <= 0 || n2 <= 0 || n3 <= 0 || nfft <= 0) return;
  size_t ntot = (size_t)n1 * (size_t)n2 * (size_t)n3 * (size_t)nfft;
  size_t bytes = ntot * sizeof(cuDoubleComplex);
  void *d_grid = arena_get(SLOT_X, bytes);

  VESPAW_CUDA_CHECK(cudaMemcpyAsync(d_grid, h_in, bytes, cudaMemcpyHostToDevice, g_stream));
  gpu_fft_3d_device_core(dir_sign, n1, n2, n3, nfft, d_grid);
  VESPAW_CUDA_CHECK(cudaMemcpyAsync(h_out, d_grid, bytes, cudaMemcpyDeviceToHost, g_stream));
  VESPAW_CUDA_CHECK(cudaStreamSynchronize(g_stream));
}

void vespaw_gpu_fft_3d_mu(int dir_sign, int ngl, int nrl, int nr1, int nr2, int nr3,
                          int nfft, int nstripelx, const int *igtostripe,
                          const int *istripetoyz, void *h_fofg, void *h_fofr) {
  if (!g_gpu_enabled || ngl <= 0 || nrl <= 0 || nfft <= 0) return;
  size_t g_bytes = (size_t)ngl * (size_t)nfft * sizeof(cuDoubleComplex);
  size_t r_bytes = (size_t)nrl * (size_t)nfft * sizeof(cuDoubleComplex);
  void *d_g = arena_get(SLOT_X, g_bytes);
  void *d_r = arena_get(SLOT_Y, r_bytes);
  void *d_ig2s = arena_get(SLOT_IGTOSTRIPE, (size_t)ngl * sizeof(int));
  void *d_is2yz = arena_get(SLOT_ISTRIPETOYZ, (size_t)nstripelx * sizeof(int));

  VESPAW_CUDA_CHECK(cudaMemcpyAsync(d_ig2s, igtostripe, (size_t)ngl * sizeof(int),
                                    cudaMemcpyHostToDevice, g_stream));
  VESPAW_CUDA_CHECK(cudaMemcpyAsync(d_is2yz, istripetoyz, (size_t)nstripelx * sizeof(int),
                                    cudaMemcpyHostToDevice, g_stream));

  int threads = 256;
  dim3 grid_g((unsigned)block_grid_1d((size_t)ngl, threads), (unsigned)nfft);

  if (dir_sign > 0) {
    /* GTOR */
    VESPAW_CUDA_CHECK(cudaMemcpyAsync(d_g, h_fofg, g_bytes, cudaMemcpyHostToDevice, g_stream));
    VESPAW_CUDA_CHECK(cudaMemsetAsync(d_r, 0, r_bytes, g_stream));
    kernel_scatter_stripes_to_3d<<<grid_g, threads, 0, g_stream>>>(
        ngl, nr1, nrl, nfft, (const int *)d_ig2s, (const int *)d_is2yz,
        (const cuDoubleComplex *)d_g, (cuDoubleComplex *)d_r);
    cufftHandle plan = get_cached_cufft_plan(3, nr3, nr2, nr1, 1, nrl, nfft);
    VESPAW_CUFFT_CHECK(cufftExecZ2Z(plan, (cufftDoubleComplex *)d_r,
                                    (cufftDoubleComplex *)d_r, CUFFT_INVERSE));
    VESPAW_CUDA_CHECK(cudaMemcpyAsync(h_fofr, d_r, r_bytes, cudaMemcpyDeviceToHost, g_stream));
  } else {
    /* RTOG */
    VESPAW_CUDA_CHECK(cudaMemcpyAsync(d_r, h_fofr, r_bytes, cudaMemcpyHostToDevice, g_stream));
    cufftHandle plan = get_cached_cufft_plan(3, nr3, nr2, nr1, 1, nrl, nfft);
    VESPAW_CUFFT_CHECK(cufftExecZ2Z(plan, (cufftDoubleComplex *)d_r,
                                    (cufftDoubleComplex *)d_r, CUFFT_FORWARD));
    double scale = 1.0 / (double)nrl;
    kernel_gather_3d_to_stripes_scaled<<<grid_g, threads, 0, g_stream>>>(
        ngl, nr1, nrl, nfft, scale, (const int *)d_ig2s, (const int *)d_is2yz,
        (const cuDoubleComplex *)d_r, (cuDoubleComplex *)d_g);
    VESPAW_CUDA_CHECK(cudaMemcpyAsync(h_fofg, d_g, g_bytes, cudaMemcpyDeviceToHost, g_stream));
  }
  VESPAW_CUDA_CHECK(cudaStreamSynchronize(g_stream));
}

void vespaw_gpu_vpsi_serial(int ngl, int ndim, int nbh, int nrl, int nr1, int nr2, int nr3,
                            int nstripelx, int tsuper, int tbucket,
                            const int *igtostripe, const int *istripetoyz,
                            const int *minusg, const double *g2, const double *bucket,
                            const double *v, const void *h_psi, void *h_hpsi) {
  (void)tsuper;
  (void)minusg;
  if (!g_gpu_enabled || ngl <= 0 || nbh <= 0 || nrl <= 0) return;

  int nfft_blk_bands = (ndim == 1) ? 32 : 16;
  if (nfft_blk_bands > nbh) nfft_blk_bands = nbh;
  int max_fft = ndim * nfft_blk_bands;

  size_t psi_all_bytes = (size_t)ngl * (size_t)ndim * (size_t)nbh * sizeof(cuDoubleComplex);
  void *d_psi_all = arena_get(SLOT_X, psi_all_bytes);
  void *d_hpsi_all = arena_get(SLOT_Z, psi_all_bytes);
  void *d_grid3d = arena_get(SLOT_Y, (size_t)nrl * (size_t)max_fft * sizeof(cuDoubleComplex));
  void *d_pot = arena_get(SLOT_POT, (size_t)nrl * (size_t)(ndim * ndim) * sizeof(double));
  void *d_g2 = arena_get(SLOT_G2, (size_t)ngl * sizeof(double));
  void *d_ig2s = arena_get(SLOT_IGTOSTRIPE, (size_t)ngl * sizeof(int));
  void *d_is2yz = arena_get(SLOT_ISTRIPETOYZ, (size_t)nstripelx * sizeof(int));
  void *d_bkt = tbucket ? arena_get(SLOT_BUCKET, (size_t)ngl * sizeof(double)) : nullptr;

  VESPAW_CUDA_CHECK(cudaMemcpyAsync(d_psi_all, h_psi, psi_all_bytes, cudaMemcpyHostToDevice, g_stream));
  VESPAW_CUDA_CHECK(cudaMemcpyAsync(d_pot, v, (size_t)nrl * (size_t)(ndim * ndim) * sizeof(double),
                                    cudaMemcpyHostToDevice, g_stream));
  VESPAW_CUDA_CHECK(cudaMemcpyAsync(d_g2, g2, (size_t)ngl * sizeof(double),
                                    cudaMemcpyHostToDevice, g_stream));
  VESPAW_CUDA_CHECK(cudaMemcpyAsync(d_ig2s, igtostripe, (size_t)ngl * sizeof(int),
                                    cudaMemcpyHostToDevice, g_stream));
  VESPAW_CUDA_CHECK(cudaMemcpyAsync(d_is2yz, istripetoyz, (size_t)nstripelx * sizeof(int),
                                    cudaMemcpyHostToDevice, g_stream));
  if (tbucket) {
    VESPAW_CUDA_CHECK(cudaMemcpyAsync(d_bkt, bucket, (size_t)ngl * sizeof(double),
                                      cudaMemcpyHostToDevice, g_stream));
  }

  double scale = 1.0 / (double)nrl;
  int threads = 256;

  for (int ib0 = 0; ib0 < nbh; ib0 += nfft_blk_bands) {
    int nblk = nbh - ib0;
    if (nblk > nfft_blk_bands) nblk = nfft_blk_bands;
    int nblk_fft = ndim * nblk;

    const cuDoubleComplex *d_psi_blk =
        ((const cuDoubleComplex *)d_psi_all) + (size_t)ib0 * (size_t)ndim * (size_t)ngl;
    cuDoubleComplex *d_hpsi_blk =
        ((cuDoubleComplex *)d_hpsi_all) + (size_t)ib0 * (size_t)ndim * (size_t)ngl;

    VESPAW_CUDA_CHECK(cudaMemsetAsync(d_grid3d, 0,
                                      (size_t)nrl * (size_t)nblk_fft * sizeof(cuDoubleComplex),
                                      g_stream));
    dim3 grid_g((unsigned)block_grid_1d((size_t)ngl, threads), (unsigned)nblk_fft);
    kernel_scatter_stripes_to_3d<<<grid_g, threads, 0, g_stream>>>(
        ngl, nr1, nrl, nblk_fft, (const int *)d_ig2s, (const int *)d_is2yz,
        d_psi_blk, (cuDoubleComplex *)d_grid3d);

    cufftHandle plan = get_cached_cufft_plan(3, nr3, nr2, nr1, 1, nrl, nblk_fft);
    VESPAW_CUFFT_CHECK(cufftExecZ2Z(plan, (cufftDoubleComplex *)d_grid3d,
                                    (cufftDoubleComplex *)d_grid3d, CUFFT_INVERSE));

    dim3 grid_r((unsigned)block_grid_1d((size_t)nrl, threads), (unsigned)nblk);
    if (ndim == 1) {
      kernel_vpsi_mul_ndim1<<<grid_r, threads, 0, g_stream>>>(
          nrl, nblk, (const double *)d_pot, (cuDoubleComplex *)d_grid3d);
    } else {
      kernel_vpsi_mul_ndim2<<<grid_r, threads, 0, g_stream>>>(
          nrl, nblk, (const double *)d_pot, (cuDoubleComplex *)d_grid3d);
    }

    VESPAW_CUFFT_CHECK(cufftExecZ2Z(plan, (cufftDoubleComplex *)d_grid3d,
                                    (cufftDoubleComplex *)d_grid3d, CUFFT_FORWARD));

    kernel_gather_3d_vpsi_kin<<<grid_g, threads, 0, g_stream>>>(
        ngl, nr1, nrl, nblk_fft, tbucket, scale, (const int *)d_ig2s, (const int *)d_is2yz,
        (const double *)d_g2, (const double *)d_bkt,
        (const cuDoubleComplex *)d_grid3d, d_psi_blk, d_hpsi_blk);
  }

  VESPAW_CUDA_CHECK(cudaMemcpyAsync(h_hpsi, d_hpsi_all, psi_all_bytes,
                                    cudaMemcpyDeviceToHost, g_stream));
  VESPAW_CUDA_CHECK(cudaStreamSynchronize(g_stream));
}

void vespaw_gpu_density_serial(int ngl, int ngl_full, int ndim, int nbh, int nb, int nrl,
                               int nr1, int nr2, int nr3, int nstripelx,
                               int tsuper, int tinv, const double *occ,
                               const int *igtostripe, const int *istripetoyz,
                               const int *minusg, const void *h_ei2kr,
                               const void *h_psi, double *h_rho) {
  if (!g_gpu_enabled || ngl <= 0 || ngl_full <= 0 || nrl <= 0 || nbh <= 0) return;

  size_t rho_bytes = (size_t)nrl * (size_t)(ndim * ndim) * sizeof(double);
  void *d_rho = arena_get(SLOT_RHO, rho_bytes);
  VESPAW_CUDA_CHECK(cudaMemsetAsync(d_rho, 0, rho_bytes, g_stream));

  void *d_occ = arena_get(SLOT_OCC, (size_t)nb * sizeof(double));
  void *d_ig2s = arena_get(SLOT_IGTOSTRIPE, (size_t)ngl_full * sizeof(int));
  void *d_is2yz = arena_get(SLOT_ISTRIPETOYZ, (size_t)nstripelx * sizeof(int));
  VESPAW_CUDA_CHECK(cudaMemcpyAsync(d_occ, occ, (size_t)nb * sizeof(double),
                                    cudaMemcpyHostToDevice, g_stream));
  VESPAW_CUDA_CHECK(cudaMemcpyAsync(d_ig2s, igtostripe, (size_t)ngl_full * sizeof(int),
                                    cudaMemcpyHostToDevice, g_stream));
  VESPAW_CUDA_CHECK(cudaMemcpyAsync(d_is2yz, istripetoyz, (size_t)nstripelx * sizeof(int),
                                    cudaMemcpyHostToDevice, g_stream));

  int threads = 256;
  dim3 grid_r((unsigned)block_grid_1d((size_t)nrl, threads), 1);

  if (tsuper) {
    /* Count independent G>=0 vectors and precompute 0-based ig1, ig2, and idx3d maps */
    int nhalf = 0;
    for (int ig1 = 0; ig1 < ngl_full; ++ig1) {
      if (minusg[ig1] - 1 >= ig1) ++nhalf;
    }
    size_t idx_bytes = ((size_t)nhalf * 2 + (size_t)ngl_full) * sizeof(int);
    int *h_tmp_idx = (int *)malloc(idx_bytes);
    int *h_ig1 = h_tmp_idx;
    int *h_ig2 = h_tmp_idx + nhalf;
    int *h_idx3d = h_ig2 + nhalf;
    int igh = 0;
    for (int ig1 = 0; ig1 < ngl_full; ++ig1) {
      int ig2 = minusg[ig1] - 1;
      if (ig2 >= ig1) {
        h_ig1[igh] = ig1;
        h_ig2[igh] = ig2;
        ++igh;
      }
      int s = igtostripe[ig1] - 1;
      int ir1 = s % nr1;
      int istripel = s / nr1;
      int i23 = istripetoyz[istripel] - 1;
      h_idx3d[ig1] = ir1 + nr1 * i23;
    }

    size_t psi_bytes = (size_t)nhalf * (size_t)nb * sizeof(cuDoubleComplex);
    void *d_psi = arena_get(SLOT_X, psi_bytes);
    void *d_idx_dev = arena_get(SLOT_INDEX_MISC, idx_bytes);
    VESPAW_CUDA_CHECK(cudaMemcpyAsync(d_psi, h_psi, psi_bytes, cudaMemcpyHostToDevice, g_stream));
    VESPAW_CUDA_CHECK(cudaMemcpyAsync(d_idx_dev, h_tmp_idx, idx_bytes,
                                      cudaMemcpyHostToDevice, g_stream));

    int npairs = (nb + 1) / 2;
    int pair_blk = (npairs > 32) ? 32 : npairs;
    void *d_grid3d = arena_get(SLOT_Y, (size_t)nrl * (size_t)pair_blk * sizeof(cuDoubleComplex));

    const int *d_ig1 = (const int *)d_idx_dev;
    const int *d_ig2 = d_ig1 + nhalf;
    const int *d_idx3d = d_ig2 + nhalf;

    for (int p0 = 0; p0 < npairs; p0 += pair_blk) {
      int np = npairs - p0;
      if (np > pair_blk) np = pair_blk;
      VESPAW_CUDA_CHECK(cudaMemsetAsync(d_grid3d, 0,
                                        (size_t)nrl * (size_t)np * sizeof(cuDoubleComplex),
                                        g_stream));
      dim3 grid_half((unsigned)block_grid_1d((size_t)nhalf, threads), (unsigned)np);
      const cuDoubleComplex *d_psi_sub =
          ((const cuDoubleComplex *)d_psi) + (size_t)(2 * p0) * (size_t)nhalf;
      kernel_scatter_tsuper_mapped<<<grid_half, threads, 0, g_stream>>>(
          nhalf, nhalf, nr1, nrl, np, nb - 2 * p0,
          d_ig1, d_ig2, d_idx3d, d_psi_sub, (cuDoubleComplex *)d_grid3d);

      cufftHandle plan = get_cached_cufft_plan(3, nr3, nr2, nr1, 1, nrl, np);
      VESPAW_CUFFT_CHECK(cufftExecZ2Z(plan, (cufftDoubleComplex *)d_grid3d,
                                      (cufftDoubleComplex *)d_grid3d, CUFFT_INVERSE));
      kernel_accum_rho_tsuper<<<grid_r, threads, 0, g_stream>>>(
          nrl, np, nb - 2 * p0, ((const double *)d_occ) + 2 * p0,
          (const cuDoubleComplex *)d_grid3d, (double *)d_rho);
    }

    VESPAW_CUDA_CHECK(cudaStreamSynchronize(g_stream));
    free(h_tmp_idx);
  } else {
    size_t psi_bytes = (size_t)ngl * (size_t)ndim * (size_t)nbh * sizeof(cuDoubleComplex);
    void *d_psi = arena_get(SLOT_X, psi_bytes);
    VESPAW_CUDA_CHECK(cudaMemcpyAsync(d_psi, h_psi, psi_bytes, cudaMemcpyHostToDevice, g_stream));

    void *d_ei2kr = nullptr;
    if (tinv && h_ei2kr) {
      d_ei2kr = arena_get(SLOT_EI2KR, (size_t)nrl * sizeof(cuDoubleComplex));
      VESPAW_CUDA_CHECK(cudaMemcpyAsync(d_ei2kr, h_ei2kr,
                                        (size_t)nrl * sizeof(cuDoubleComplex),
                                        cudaMemcpyHostToDevice, g_stream));
    }

    int nfft_blk_bands = (ndim == 1) ? 32 : 16;
    if (nfft_blk_bands > nbh) nfft_blk_bands = nbh;
    int max_fft = ndim * nfft_blk_bands;
    void *d_grid3d = arena_get(SLOT_Y, (size_t)nrl * (size_t)max_fft * sizeof(cuDoubleComplex));

    for (int ib0 = 0; ib0 < nbh; ib0 += nfft_blk_bands) {
      int nblk = nbh - ib0;
      if (nblk > nfft_blk_bands) nblk = nfft_blk_bands;
      int nblk_fft = ndim * nblk;

      const cuDoubleComplex *d_psi_blk =
          ((const cuDoubleComplex *)d_psi) + (size_t)ib0 * (size_t)ndim * (size_t)ngl;

      VESPAW_CUDA_CHECK(cudaMemsetAsync(d_grid3d, 0,
                                        (size_t)nrl * (size_t)nblk_fft * sizeof(cuDoubleComplex),
                                        g_stream));
      dim3 grid_g((unsigned)block_grid_1d((size_t)ngl, threads), (unsigned)nblk_fft);
      kernel_scatter_stripes_to_3d<<<grid_g, threads, 0, g_stream>>>(
          ngl, nr1, nrl, nblk_fft, (const int *)d_ig2s, (const int *)d_is2yz,
          d_psi_blk, (cuDoubleComplex *)d_grid3d);

      cufftHandle plan = get_cached_cufft_plan(3, nr3, nr2, nr1, 1, nrl, nblk_fft);
      VESPAW_CUFFT_CHECK(cufftExecZ2Z(plan, (cufftDoubleComplex *)d_grid3d,
                                      (cufftDoubleComplex *)d_grid3d, CUFFT_INVERSE));

      if (tinv) {
        kernel_accum_rho_tinv<<<grid_r, threads, 0, g_stream>>>(
            nrl, nblk, ib0, nb, (const double *)d_occ,
            (const cuDoubleComplex *)d_ei2kr,
            (const cuDoubleComplex *)d_grid3d, (double *)d_rho);
      } else if (ndim == 1) {
        kernel_accum_rho_ndim1<<<grid_r, threads, 0, g_stream>>>(
            nrl, nblk, ib0, (const double *)d_occ,
            (const cuDoubleComplex *)d_grid3d, (double *)d_rho);
      } else {
        kernel_accum_rho_ndim2<<<grid_r, threads, 0, g_stream>>>(
            nrl, nblk, ib0, (const double *)d_occ,
            (const cuDoubleComplex *)d_grid3d, (double *)d_rho);
      }
    }
  }

  VESPAW_CUDA_CHECK(cudaMemcpyAsync(h_rho, d_rho, rho_bytes, cudaMemcpyDeviceToHost, g_stream));
  VESPAW_CUDA_CHECK(cudaStreamSynchronize(g_stream));
}

void vespaw_gpu_tau_serial(int mode, int ngl, int ndim, int nbh, int nb, int nrl,
                           int nr1, int nr2, int nr3, int nstripelx,
                           int tsuper, int tinv, const double *occ,
                           const double *gvec, const double *kvec_cart,
                           const int *igtostripe, const int *istripetoyz,
                           const int *minusg, const double *vtau,
                           const void *h_psi, double *h_tau, void *h_htaupsi) {
  (void)tsuper;
  (void)tinv;
  (void)kvec_cart;
  (void)minusg;
  if (!g_gpu_enabled || ndim != 1 || ngl <= 0 || nbh <= 0 || nrl <= 0) return;

  size_t psi_bytes = (size_t)ngl * (size_t)nbh * sizeof(cuDoubleComplex);
  void *d_psi = arena_get(SLOT_X, psi_bytes);
  void *d_psi_d = arena_get(SLOT_WORK, psi_bytes);
  void *d_grid3d = arena_get(SLOT_Y, (size_t)nrl * (size_t)nbh * sizeof(cuDoubleComplex));
  void *d_gvec = arena_get(SLOT_GVEC, 3 * (size_t)ngl * sizeof(double));
  void *d_ig2s = arena_get(SLOT_IGTOSTRIPE, (size_t)ngl * sizeof(int));
  void *d_is2yz = arena_get(SLOT_ISTRIPETOYZ, (size_t)nstripelx * sizeof(int));

  VESPAW_CUDA_CHECK(cudaMemcpyAsync(d_psi, h_psi, psi_bytes, cudaMemcpyHostToDevice, g_stream));
  VESPAW_CUDA_CHECK(cudaMemcpyAsync(d_gvec, gvec, 3 * (size_t)ngl * sizeof(double),
                                    cudaMemcpyHostToDevice, g_stream));
  VESPAW_CUDA_CHECK(cudaMemcpyAsync(d_ig2s, igtostripe, (size_t)ngl * sizeof(int),
                                    cudaMemcpyHostToDevice, g_stream));
  VESPAW_CUDA_CHECK(cudaMemcpyAsync(d_is2yz, istripetoyz, (size_t)nstripelx * sizeof(int),
                                    cudaMemcpyHostToDevice, g_stream));

  int threads = 256;
  dim3 grid_g((unsigned)block_grid_1d((size_t)ngl, threads), (unsigned)nbh);
  dim3 grid_r((unsigned)block_grid_1d((size_t)nrl, threads), (unsigned)nbh);
  cufftHandle plan = get_cached_cufft_plan(3, nr3, nr2, nr1, 1, nrl, nbh);

  if (mode == 0) {
    /* WAVES_KINETICDENSITY */
    void *d_occ = arena_get(SLOT_OCC, (size_t)nb * sizeof(double));
    void *d_tau = arena_get(SLOT_RHO, (size_t)nrl * sizeof(double));
    VESPAW_CUDA_CHECK(cudaMemcpyAsync(d_occ, occ, (size_t)nb * sizeof(double),
                                      cudaMemcpyHostToDevice, g_stream));
    VESPAW_CUDA_CHECK(cudaMemsetAsync(d_tau, 0, (size_t)nrl * sizeof(double), g_stream));

    for (int idir = 0; idir < 3; ++idir) {
      kernel_mul_ig_grad<<<grid_g, threads, 0, g_stream>>>(
          ngl, nbh, idir, (const double *)d_gvec,
          (const cuDoubleComplex *)d_psi, (cuDoubleComplex *)d_psi_d);
      VESPAW_CUDA_CHECK(cudaMemsetAsync(d_grid3d, 0,
                                        (size_t)nrl * (size_t)nbh * sizeof(cuDoubleComplex),
                                        g_stream));
      kernel_scatter_stripes_to_3d<<<grid_g, threads, 0, g_stream>>>(
          ngl, nr1, nrl, nbh, (const int *)d_ig2s, (const int *)d_is2yz,
          (const cuDoubleComplex *)d_psi_d, (cuDoubleComplex *)d_grid3d);
      VESPAW_CUFFT_CHECK(cufftExecZ2Z(plan, (cufftDoubleComplex *)d_grid3d,
                                      (cufftDoubleComplex *)d_grid3d, CUFFT_INVERSE));
      kernel_accum_rho_ndim1<<<grid_r.x, threads, 0, g_stream>>>(
          nrl, nbh, 0, (const double *)d_occ,
          (const cuDoubleComplex *)d_grid3d, (double *)d_tau);
    }

    kernel_scale_r8<<<grid_r.x, threads, 0, g_stream>>>(
        (double *)d_tau, (size_t)nrl, 0.5);
    VESPAW_CUDA_CHECK(cudaMemcpyAsync(h_tau, d_tau, (size_t)nrl * sizeof(double),
                                      cudaMemcpyDeviceToHost, g_stream));
    VESPAW_CUDA_CHECK(cudaStreamSynchronize(g_stream));
  } else {
    /* WAVES_TAUHPSI */
    void *d_vtau = arena_get(SLOT_VTAU, (size_t)nrl * sizeof(double));
    void *d_htau = arena_get(SLOT_Z, psi_bytes);
    VESPAW_CUDA_CHECK(cudaMemcpyAsync(d_vtau, vtau, (size_t)nrl * sizeof(double),
                                      cudaMemcpyHostToDevice, g_stream));
    VESPAW_CUDA_CHECK(cudaMemsetAsync(d_htau, 0, psi_bytes, g_stream));
    double scale = 1.0 / (double)nrl;

    for (int idir = 0; idir < 3; ++idir) {
      kernel_mul_ig_grad<<<grid_g, threads, 0, g_stream>>>(
          ngl, nbh, idir, (const double *)d_gvec,
          (const cuDoubleComplex *)d_psi, (cuDoubleComplex *)d_psi_d);
      VESPAW_CUDA_CHECK(cudaMemsetAsync(d_grid3d, 0,
                                        (size_t)nrl * (size_t)nbh * sizeof(cuDoubleComplex),
                                        g_stream));
      kernel_scatter_stripes_to_3d<<<grid_g, threads, 0, g_stream>>>(
          ngl, nr1, nrl, nbh, (const int *)d_ig2s, (const int *)d_is2yz,
          (const cuDoubleComplex *)d_psi_d, (cuDoubleComplex *)d_grid3d);
      VESPAW_CUFFT_CHECK(cufftExecZ2Z(plan, (cufftDoubleComplex *)d_grid3d,
                                      (cufftDoubleComplex *)d_grid3d, CUFFT_INVERSE));
      kernel_vpsi_mul_ndim1<<<grid_r, threads, 0, g_stream>>>(
          nrl, nbh, (const double *)d_vtau, (cuDoubleComplex *)d_grid3d);
      VESPAW_CUFFT_CHECK(cufftExecZ2Z(plan, (cufftDoubleComplex *)d_grid3d,
                                      (cufftDoubleComplex *)d_grid3d, CUFFT_FORWARD));
      kernel_gather_3d_to_stripes_scaled<<<grid_g, threads, 0, g_stream>>>(
          ngl, nr1, nrl, nbh, scale, (const int *)d_ig2s, (const int *)d_is2yz,
          (const cuDoubleComplex *)d_grid3d, (cuDoubleComplex *)d_psi_d);
      kernel_accum_ig_htaupsi<<<grid_g, threads, 0, g_stream>>>(
          ngl, nbh, idir, (const double *)d_gvec,
          (const cuDoubleComplex *)d_psi_d, (cuDoubleComplex *)d_htau);
    }
    size_t ntot = (size_t)ngl * (size_t)nbh;
    kernel_scale_c8<<<(unsigned)block_grid_1d(ntot, threads), threads, 0, g_stream>>>(
        (cuDoubleComplex *)d_htau, ntot, -0.5);
    VESPAW_CUDA_CHECK(cudaMemcpyAsync(h_htaupsi, d_htau, psi_bytes,
                                      cudaMemcpyDeviceToHost, g_stream));
    VESPAW_CUDA_CHECK(cudaStreamSynchronize(g_stream));
  }
}

void vespaw_gpu_procache_upload(const char *gset_id, int id_len, int nat,
                                int ngl, int npro, const double *r,
                                const void *h_procache) {
  (void)r;
  if (!g_gpu_enabled || ngl <= 0 || npro <= 0) return;
  size_t bytes = (size_t)ngl * (size_t)npro * sizeof(cuDoubleComplex);
  if (bytes > g_procache_cap) {
    if (g_procache_dev) {
      VESPAW_CUDA_CHECK(cudaStreamSynchronize(g_stream));
      VESPAW_CUDA_CHECK(cudaFree(g_procache_dev));
      g_procache_dev = nullptr;
      g_procache_cap = 0;
    }
    VESPAW_CUDA_CHECK(cudaMalloc(&g_procache_dev, bytes));
    g_procache_cap = bytes;
  }
  VESPAW_CUDA_CHECK(cudaMemcpyAsync(g_procache_dev, h_procache, bytes,
                                    cudaMemcpyHostToDevice, g_stream));
  VESPAW_CUDA_CHECK(cudaStreamSynchronize(g_stream));
  g_procache_ngl = ngl;
  g_procache_npro = npro;
  g_procache_nat = nat;
  memset(g_procache_id, 0, sizeof(g_procache_id));
  if (gset_id && id_len > 0) {
    int cpy = (id_len < 31) ? id_len : 31;
    memcpy(g_procache_id, gset_id, (size_t)cpy);
  }
}

void vespaw_gpu_projections(int ngl, int ndim, int nb, int npro,
                            int tsuper, int ngamma, double gweight1,
                            const void *h_psi, void *h_propsi) {
  if (!g_gpu_enabled || !g_procache_dev || g_procache_ngl != ngl || g_procache_npro != npro) {
    return;
  }
  if (ngl <= 0 || nb <= 0 || npro <= 0) return;

  int nb_eff = ndim * nb;
  size_t psi_bytes = (size_t)ngl * (size_t)nb_eff * sizeof(cuDoubleComplex);
  size_t out_bytes = (size_t)npro * (size_t)nb_eff * sizeof(cuDoubleComplex);
  void *d_psi = arena_get(SLOT_X, psi_bytes);
  void *d_tmp1 = arena_get(SLOT_PRO_TMP1, out_bytes);
  void *d_tmp2 = arena_get(SLOT_PRO_TMP2, out_bytes);

  VESPAW_CUDA_CHECK(cudaMemcpyAsync(d_psi, h_psi, psi_bytes, cudaMemcpyHostToDevice, g_stream));

  cuDoubleComplex one = make_cuDoubleComplex(1.0, 0.0);
  cuDoubleComplex zero = make_cuDoubleComplex(0.0, 0.0);
  VESPAW_CUBLAS_CHECK(cublasZgemm(g_cublas, CUBLAS_OP_C, CUBLAS_OP_N,
                                  npro, nb_eff, ngl,
                                  &one,
                                  (const cuDoubleComplex *)g_procache_dev, ngl,
                                  (const cuDoubleComplex *)d_psi, ngl,
                                  &zero,
                                  (cuDoubleComplex *)d_tmp1, npro));

  dim3 blk(16, 16);
  dim3 grd1((unsigned)((npro + 15) / 16), (unsigned)((nb_eff + 15) / 16));
  kernel_adjust_tsuper_gweight<<<grd1, blk, 0, g_stream>>>(
      npro, nb_eff, ngl, 1, tsuper, ngamma, gweight1,
      (const cuDoubleComplex *)g_procache_dev,
      (const cuDoubleComplex *)d_psi,
      (cuDoubleComplex *)d_tmp1);

  kernel_transpose_c8<<<grd1, blk, 0, g_stream>>>(
      npro, nb_eff, (const cuDoubleComplex *)d_tmp1, (cuDoubleComplex *)d_tmp2);

  VESPAW_CUDA_CHECK(cudaMemcpyAsync(h_propsi, d_tmp2, out_bytes,
                                    cudaMemcpyDeviceToHost, g_stream));
  VESPAW_CUDA_CHECK(cudaStreamSynchronize(g_stream));
}

void vespaw_gpu_addpro(int ngl, int ndim, int nb, int npro,
                       const void *h_propsi, void *h_psi) {
  if (!g_gpu_enabled || !g_procache_dev || g_procache_ngl != ngl || g_procache_npro != npro) {
    return;
  }
  if (ngl <= 0 || nb <= 0 || npro <= 0) return;

  int nb_eff = ndim * nb;
  size_t psi_bytes = (size_t)ngl * (size_t)nb_eff * sizeof(cuDoubleComplex);
  size_t propsi_bytes = (size_t)nb_eff * (size_t)npro * sizeof(cuDoubleComplex);
  void *d_psi = arena_get(SLOT_X, psi_bytes);
  void *d_propsi = arena_get(SLOT_PRO_TMP1, propsi_bytes);

  VESPAW_CUDA_CHECK(cudaMemcpyAsync(d_psi, h_psi, psi_bytes, cudaMemcpyHostToDevice, g_stream));
  VESPAW_CUDA_CHECK(cudaMemcpyAsync(d_propsi, h_propsi, propsi_bytes, cudaMemcpyHostToDevice, g_stream));

  cuDoubleComplex one = make_cuDoubleComplex(1.0, 0.0);
  VESPAW_CUBLAS_CHECK(cublasZgemm(g_cublas, CUBLAS_OP_N, CUBLAS_OP_T,
                                  ngl, nb_eff, npro,
                                  &one,
                                  (const cuDoubleComplex *)g_procache_dev, ngl,
                                  (const cuDoubleComplex *)d_propsi, nb_eff,
                                  &one,
                                  (cuDoubleComplex *)d_psi, ngl));

  VESPAW_CUDA_CHECK(cudaMemcpyAsync(h_psi, d_psi, psi_bytes, cudaMemcpyDeviceToHost, g_stream));
  VESPAW_CUDA_CHECK(cudaStreamSynchronize(g_stream));
}

void vespaw_gpu_scalarproduct(int is_hermitian, int ngl, int ndim,
                              int n1, const void *h_f1,
                              int n2, const void *h_f2,
                              int tsuper, int ngamma, double gweight1,
                              void *h_prod, int *done_out) {
  if (!g_gpu_enabled || ngl <= 0 || ndim <= 0 || n1 <= 0 || n2 <= 0) {
    if (done_out) *done_out = 0;
    return;
  }
  size_t work_est = (size_t)ngl * (size_t)ndim * (size_t)n1 * (size_t)n2;
  if (work_est < GPU_MAT_THRESHOLD) {
    if (done_out) *done_out = 0;
    return;
  }

  int ng = ngl * ndim;
  size_t f1_bytes = (size_t)ng * (size_t)n1 * sizeof(cuDoubleComplex);
  size_t f2_bytes = (size_t)ng * (size_t)n2 * sizeof(cuDoubleComplex);
  size_t prod_bytes = (size_t)n1 * (size_t)n2 * sizeof(cuDoubleComplex);

  void *d_f1 = arena_get(SLOT_X, f1_bytes);
  void *d_f2 = is_hermitian ? d_f1 : arena_get(SLOT_Y, f2_bytes);
  void *d_prod = arena_get(SLOT_Z, prod_bytes);

  VESPAW_CUDA_CHECK(cudaMemcpyAsync(d_f1, h_f1, f1_bytes, cudaMemcpyHostToDevice, g_stream));
  if (!is_hermitian) {
    VESPAW_CUDA_CHECK(cudaMemcpyAsync(d_f2, h_f2, f2_bytes, cudaMemcpyHostToDevice, g_stream));
  }

  cuDoubleComplex one = make_cuDoubleComplex(1.0, 0.0);
  cuDoubleComplex zero = make_cuDoubleComplex(0.0, 0.0);
  VESPAW_CUBLAS_CHECK(cublasZgemm(g_cublas, CUBLAS_OP_C, CUBLAS_OP_N,
                                  n1, n2, ng,
                                  &one,
                                  (const cuDoubleComplex *)d_f1, ng,
                                  (const cuDoubleComplex *)d_f2, ng,
                                  &zero,
                                  (cuDoubleComplex *)d_prod, n1));

  dim3 blk(16, 16);
  dim3 grd((unsigned)((n1 + 15) / 16), (unsigned)((n2 + 15) / 16));
  kernel_adjust_tsuper_gweight<<<grd, blk, 0, g_stream>>>(
      n1, n2, ngl, ndim, tsuper, ngamma, gweight1,
      (const cuDoubleComplex *)d_f1,
      (const cuDoubleComplex *)d_f2,
      (cuDoubleComplex *)d_prod);

  VESPAW_CUDA_CHECK(cudaMemcpyAsync(h_prod, d_prod, prod_bytes, cudaMemcpyDeviceToHost, g_stream));
  VESPAW_CUDA_CHECK(cudaStreamSynchronize(g_stream));
  if (done_out) *done_out = 1;
}

void vespaw_gpu_dotdiag(int ngl, int ndim, int n,
                        const void *h_f1, const void *h_f2,
                        int tsuper, int ngamma, double gweight1,
                        double *h_diag, int *done_out) {
  if (!g_gpu_enabled || ngl <= 0 || ndim <= 0 || n <= 0) {
    if (done_out) *done_out = 0;
    return;
  }
  size_t vec_tot = (size_t)ngl * (size_t)ndim * (size_t)n;
  if (vec_tot < GPU_MAT_THRESHOLD) {
    if (done_out) *done_out = 0;
    return;
  }
  size_t bytes = vec_tot * sizeof(cuDoubleComplex);
  void *d_f1 = arena_get(SLOT_X, bytes);
  void *d_f2 = (h_f1 == h_f2) ? d_f1 : arena_get(SLOT_Y, bytes);
  void *d_diag = arena_get(SLOT_RHO, (size_t)n * sizeof(double));

  VESPAW_CUDA_CHECK(cudaMemcpyAsync(d_f1, h_f1, bytes, cudaMemcpyHostToDevice, g_stream));
  if (h_f1 != h_f2) {
    VESPAW_CUDA_CHECK(cudaMemcpyAsync(d_f2, h_f2, bytes, cudaMemcpyHostToDevice, g_stream));
  }

  kernel_dotdiag_c8<<<(unsigned)n, 256, 0, g_stream>>>(
      ngl, ndim, n, tsuper, ngamma, gweight1,
      (const cuDoubleComplex *)d_f1,
      (const cuDoubleComplex *)d_f2,
      (double *)d_diag);

  VESPAW_CUDA_CHECK(cudaMemcpyAsync(h_diag, d_diag, (size_t)n * sizeof(double),
                                    cudaMemcpyDeviceToHost, g_stream));
  VESPAW_CUDA_CHECK(cudaStreamSynchronize(g_stream));
  if (done_out) *done_out = 1;
}

void vespaw_gpu_zgemm(int transa, int transb, int m, int n, int k,
                      const void *alpha, const void *h_a, int lda,
                      const void *h_b, int ldb, const void *beta,
                      void *h_c, int ldc) {
  if (!g_gpu_enabled || m <= 0 || n <= 0 || k <= 0) return;

  auto map_op = [](int t) -> cublasOperation_t {
    if (t == 1) return CUBLAS_OP_T;
    if (t == 2) return CUBLAS_OP_C;
    return CUBLAS_OP_N;
  };
  int a_cols = (transa == 0) ? k : m;
  int b_cols = (transb == 0) ? n : k;
  size_t a_bytes = (size_t)lda * (size_t)a_cols * sizeof(cuDoubleComplex);
  size_t b_bytes = (size_t)ldb * (size_t)b_cols * sizeof(cuDoubleComplex);
  size_t c_bytes = (size_t)ldc * (size_t)n * sizeof(cuDoubleComplex);

  void *d_a = arena_get(SLOT_X, a_bytes);
  void *d_b = arena_get(SLOT_Y, b_bytes);
  void *d_c = arena_get(SLOT_Z, c_bytes);

  VESPAW_CUDA_CHECK(cudaMemcpyAsync(d_a, h_a, a_bytes, cudaMemcpyHostToDevice, g_stream));
  VESPAW_CUDA_CHECK(cudaMemcpyAsync(d_b, h_b, b_bytes, cudaMemcpyHostToDevice, g_stream));
  const cuDoubleComplex *b_val = (const cuDoubleComplex *)beta;
  if (cuCreal(*b_val) != 0.0 || cuCimag(*b_val) != 0.0) {
    VESPAW_CUDA_CHECK(cudaMemcpyAsync(d_c, h_c, c_bytes, cudaMemcpyHostToDevice, g_stream));
  }

  VESPAW_CUBLAS_CHECK(cublasZgemm(g_cublas, map_op(transa), map_op(transb),
                                  m, n, k,
                                  (const cuDoubleComplex *)alpha,
                                  (const cuDoubleComplex *)d_a, lda,
                                  (const cuDoubleComplex *)d_b, ldb,
                                  b_val,
                                  (cuDoubleComplex *)d_c, ldc));

  VESPAW_CUDA_CHECK(cudaMemcpyAsync(h_c, d_c, c_bytes, cudaMemcpyDeviceToHost, g_stream));
  VESPAW_CUDA_CHECK(cudaStreamSynchronize(g_stream));
}

void vespaw_gpu_zhegvd(int n, const void *h_h, const void *h_s,
                       double *h_e, void *h_vec, int *info_out) {
  if (!g_gpu_enabled || n < GPU_EIG_CUTOFF) {
    if (info_out) *info_out = -999;
    return;
  }
  size_t mat_bytes = (size_t)n * (size_t)n * sizeof(cuDoubleComplex);
  size_t eig_bytes = (size_t)n * sizeof(double);

  void *d_a = arena_get(SLOT_EIG_A, mat_bytes);
  void *d_b = arena_get(SLOT_EIG_B, mat_bytes);
  void *d_w = arena_get(SLOT_EIG_W, eig_bytes);
  void *d_info = arena_get(SLOT_EIG_INFO, sizeof(int));

  VESPAW_CUDA_CHECK(cudaMemcpyAsync(d_a, h_h, mat_bytes, cudaMemcpyHostToDevice, g_stream));
  VESPAW_CUDA_CHECK(cudaMemcpyAsync(d_b, h_s, mat_bytes, cudaMemcpyHostToDevice, g_stream));

  int lwork = 0;
  VESPAW_CUSOLVER_CHECK(cusolverDnZhegvd_bufferSize(
      g_cusolver, CUSOLVER_EIG_TYPE_1, CUSOLVER_EIG_MODE_VECTOR, CUBLAS_FILL_MODE_UPPER,
      n, (const cuDoubleComplex *)d_a, n, (const cuDoubleComplex *)d_b, n,
      (const double *)d_w, &lwork));

  void *d_work = arena_get(SLOT_EIG_WORK, (size_t)lwork * sizeof(cuDoubleComplex));
  VESPAW_CUSOLVER_CHECK(cusolverDnZhegvd(
      g_cusolver, CUSOLVER_EIG_TYPE_1, CUSOLVER_EIG_MODE_VECTOR, CUBLAS_FILL_MODE_UPPER,
      n, (cuDoubleComplex *)d_a, n, (cuDoubleComplex *)d_b, n,
      (double *)d_w, (cuDoubleComplex *)d_work, lwork, (int *)d_info));

  int info = 0;
  VESPAW_CUDA_CHECK(cudaMemcpyAsync(&info, d_info, sizeof(int), cudaMemcpyDeviceToHost, g_stream));
  VESPAW_CUDA_CHECK(cudaMemcpyAsync(h_e, d_w, eig_bytes, cudaMemcpyDeviceToHost, g_stream));
  VESPAW_CUDA_CHECK(cudaMemcpyAsync(h_vec, d_a, mat_bytes, cudaMemcpyDeviceToHost, g_stream));
  VESPAW_CUDA_CHECK(cudaStreamSynchronize(g_stream));

  if (info_out) *info_out = info;
}

void vespaw_gpu_zheevd(int n, const void *h_h, double *h_e,
                       void *h_vec, int *info_out) {
  if (!g_gpu_enabled || n < GPU_EIG_CUTOFF) {
    if (info_out) *info_out = -999;
    return;
  }
  size_t mat_bytes = (size_t)n * (size_t)n * sizeof(cuDoubleComplex);
  size_t eig_bytes = (size_t)n * sizeof(double);

  void *d_a = arena_get(SLOT_EIG_A, mat_bytes);
  void *d_w = arena_get(SLOT_EIG_W, eig_bytes);
  void *d_info = arena_get(SLOT_EIG_INFO, sizeof(int));

  VESPAW_CUDA_CHECK(cudaMemcpyAsync(d_a, h_h, mat_bytes, cudaMemcpyHostToDevice, g_stream));

  int lwork = 0;
  VESPAW_CUSOLVER_CHECK(cusolverDnZheevd_bufferSize(
      g_cusolver, CUSOLVER_EIG_MODE_VECTOR, CUBLAS_FILL_MODE_UPPER,
      n, (const cuDoubleComplex *)d_a, n, (const double *)d_w, &lwork));

  void *d_work = arena_get(SLOT_EIG_WORK, (size_t)lwork * sizeof(cuDoubleComplex));
  VESPAW_CUSOLVER_CHECK(cusolverDnZheevd(
      g_cusolver, CUSOLVER_EIG_MODE_VECTOR, CUBLAS_FILL_MODE_UPPER,
      n, (cuDoubleComplex *)d_a, n, (double *)d_w,
      (cuDoubleComplex *)d_work, lwork, (int *)d_info));

  int info = 0;
  VESPAW_CUDA_CHECK(cudaMemcpyAsync(&info, d_info, sizeof(int), cudaMemcpyDeviceToHost, g_stream));
  VESPAW_CUDA_CHECK(cudaMemcpyAsync(h_e, d_w, eig_bytes, cudaMemcpyDeviceToHost, g_stream));
  VESPAW_CUDA_CHECK(cudaMemcpyAsync(h_vec, d_a, mat_bytes, cudaMemcpyDeviceToHost, g_stream));
  VESPAW_CUDA_CHECK(cudaStreamSynchronize(g_stream));

  if (info_out) *info_out = info;
}

void vespaw_gpu_precondition_teter(int ngl, int ndim, int nact,
                                   const double *h_g2, const double *h_ekin,
                                   const void *h_r, void *h_d, int *done_out) {
  int ncols = ndim * nact;
  if (!g_gpu_enabled || ngl <= 0 || ncols <= 0 || (size_t)ngl * (size_t)ncols < GPU_MAT_THRESHOLD) {
    if (done_out) *done_out = 0;
    return;
  }
  size_t vec_bytes = (size_t)ngl * (size_t)ncols * sizeof(cuDoubleComplex);
  void *d_g2 = arena_get(SLOT_G2, (size_t)ngl * sizeof(double));
  void *d_ek = arena_get(SLOT_OCC, (size_t)ncols * sizeof(double));
  void *d_r = arena_get(SLOT_X, vec_bytes);
  void *d_d = arena_get(SLOT_Y, vec_bytes);

  VESPAW_CUDA_CHECK(cudaMemcpyAsync(d_g2, h_g2, (size_t)ngl * sizeof(double),
                                    cudaMemcpyHostToDevice, g_stream));
  VESPAW_CUDA_CHECK(cudaMemcpyAsync(d_ek, h_ekin, (size_t)ncols * sizeof(double),
                                    cudaMemcpyHostToDevice, g_stream));
  VESPAW_CUDA_CHECK(cudaMemcpyAsync(d_r, h_r, vec_bytes, cudaMemcpyHostToDevice, g_stream));

  int threads = 256;
  dim3 grd((unsigned)block_grid_1d((size_t)ngl, threads), (unsigned)ncols);
  kernel_precondition_teter<<<grd, threads, 0, g_stream>>>(
      ngl, ncols, (const double *)d_g2, (const double *)d_ek,
      (const cuDoubleComplex *)d_r, (cuDoubleComplex *)d_d);

  VESPAW_CUDA_CHECK(cudaMemcpyAsync(h_d, d_d, vec_bytes, cudaMemcpyDeviceToHost, g_stream));
  VESPAW_CUDA_CHECK(cudaStreamSynchronize(g_stream));
  if (done_out) *done_out = 1;
}

void vespaw_gpu_expandpro_all(int ngl, int nat, int nsp, int lmnxx, int lnxx,
                              int lmx, int nbarepro, int npro,
                              const int *h_isp, const int *h_lnx,
                              const int *h_lmnx, const int *h_lox,
                              const double *h_r, const double *h_gvec,
                              const double *h_barepro, const double *h_ylm,
                              void *h_procache, int *done_out) {
  (void)lmnxx;
  if (!g_gpu_enabled || ngl <= 0 || nat <= 0 || npro <= 0) {
    if (done_out) *done_out = 0;
    return;
  }

  /* Precompute 0-based atom_ipro0(nat) and species_ibpro0(nsp) cleanly */
  size_t int_count = (size_t)nat * 2 + (size_t)nsp * 2 + (size_t)lnxx * (size_t)nsp;
  int *h_ibuf = (int *)malloc(int_count * sizeof(int));
  int *h_atom_isp = h_ibuf;
  int *h_atom_ipro = h_atom_isp + nat;
  int *h_sp_ibpro = h_atom_ipro + nat;
  int *h_sp_lnx = h_sp_ibpro + nsp;
  int *h_sp_lox = h_sp_lnx + nsp;

  int acc = 0;
  for (int isp = 0; isp < nsp; ++isp) {
    h_sp_ibpro[isp] = acc;
    h_sp_lnx[isp] = h_lnx[isp];
    acc += h_lnx[isp];
  }
  acc = 0;
  for (int iat = 0; iat < nat; ++iat) {
    int isp1 = h_isp[iat];
    h_atom_isp[iat] = isp1;
    h_atom_ipro[iat] = acc;
    acc += h_lmnx[isp1 - 1];
  }
  memcpy(h_sp_lox, h_lox, (size_t)lnxx * (size_t)nsp * sizeof(int));

  void *d_ibuf = arena_get(SLOT_INDEX_MISC, int_count * sizeof(int));
  void *d_r = arena_get(SLOT_WORK, 3 * (size_t)nat * sizeof(double));
  void *d_gvec = arena_get(SLOT_GVEC, 3 * (size_t)ngl * sizeof(double));
  void *d_bpro = arena_get(SLOT_POT, (size_t)ngl * (size_t)nbarepro * sizeof(double));
  void *d_ylm = arena_get(SLOT_G2, (size_t)ngl * (size_t)lmx * sizeof(double));

  size_t pro_bytes = (size_t)ngl * (size_t)npro * sizeof(cuDoubleComplex);
  if (pro_bytes > g_procache_cap) {
    if (g_procache_dev) {
      VESPAW_CUDA_CHECK(cudaStreamSynchronize(g_stream));
      VESPAW_CUDA_CHECK(cudaFree(g_procache_dev));
      g_procache_dev = nullptr;
      g_procache_cap = 0;
    }
    VESPAW_CUDA_CHECK(cudaMalloc(&g_procache_dev, pro_bytes));
    g_procache_cap = pro_bytes;
  }

  VESPAW_CUDA_CHECK(cudaMemcpyAsync(d_ibuf, h_ibuf, int_count * sizeof(int),
                                    cudaMemcpyHostToDevice, g_stream));
  VESPAW_CUDA_CHECK(cudaMemcpyAsync(d_r, h_r, 3 * (size_t)nat * sizeof(double),
                                    cudaMemcpyHostToDevice, g_stream));
  VESPAW_CUDA_CHECK(cudaMemcpyAsync(d_gvec, h_gvec, 3 * (size_t)ngl * sizeof(double),
                                    cudaMemcpyHostToDevice, g_stream));
  VESPAW_CUDA_CHECK(cudaMemcpyAsync(d_bpro, h_barepro,
                                    (size_t)ngl * (size_t)nbarepro * sizeof(double),
                                    cudaMemcpyHostToDevice, g_stream));
  VESPAW_CUDA_CHECK(cudaMemcpyAsync(d_ylm, h_ylm, (size_t)ngl * (size_t)lmx * sizeof(double),
                                    cudaMemcpyHostToDevice, g_stream));

  const int *d_atom_isp = (const int *)d_ibuf;
  const int *d_atom_ipro = d_atom_isp + nat;
  const int *d_sp_ibpro = d_atom_ipro + nat;
  const int *d_sp_lnx = d_sp_ibpro + nsp;
  const int *d_sp_lox = d_sp_lnx + nsp;

  int threads = 256;
  dim3 grd((unsigned)block_grid_1d((size_t)ngl, threads), (unsigned)nat);
  kernel_expandpro_all<<<grd, threads, 0, g_stream>>>(
      ngl, nat, lnxx, lmx, d_atom_isp, d_atom_ipro, d_sp_ibpro, d_sp_lnx, d_sp_lox,
      (const double *)d_r, (const double *)d_gvec,
      (const double *)d_bpro, (const double *)d_ylm,
      (cuDoubleComplex *)g_procache_dev);

  VESPAW_CUDA_CHECK(cudaMemcpyAsync(h_procache, g_procache_dev, pro_bytes,
                                    cudaMemcpyDeviceToHost, g_stream));
  VESPAW_CUDA_CHECK(cudaStreamSynchronize(g_stream));
  g_procache_ngl = ngl;
  g_procache_npro = npro;
  g_procache_nat = nat;

  free(h_ibuf);
  if (done_out) *done_out = 1;
}

void *vespaw_gpu_dev_alloc(size_t bytes) {
  if (!g_gpu_enabled || bytes == 0) return nullptr;
  void *ptr = nullptr;
  VESPAW_CUDA_CHECK(cudaMalloc(&ptr, bytes));
  return ptr;
}

void vespaw_gpu_dev_free(void *ptr) {
  if (ptr) {
    VESPAW_CUDA_CHECK(cudaFree(ptr));
  }
}

void vespaw_gpu_h2d(void *d_dst, const void *h_src, size_t bytes) {
  if (!g_gpu_enabled || !d_dst || !h_src || bytes == 0) return;
  VESPAW_CUDA_CHECK(cudaMemcpy(d_dst, h_src, bytes, cudaMemcpyHostToDevice));
}

void vespaw_gpu_d2h(void *h_dst, const void *d_src, size_t bytes) {
  if (!g_gpu_enabled || !h_dst || !d_src || bytes == 0) return;
  VESPAW_CUDA_CHECK(cudaMemcpy(h_dst, d_src, bytes, cudaMemcpyDeviceToHost));
}

} /* extern "C" */
