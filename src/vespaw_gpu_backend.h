#ifndef VESPAW_GPU_BACKEND_H
#define VESPAW_GPU_BACKEND_H

#include <stddef.h>

#ifdef __cplusplus
extern "C" {
#endif

/* Initialize / tear down persistent CUDA runtime, stream, cuBLAS, cuSOLVER, and caches. */
void vespaw_gpu_init(int mpirank, int mpisize, int *enabled_out, int *dev_id_out);
void vespaw_gpu_shutdown(void);
int  vespaw_gpu_is_enabled(void);
void vespaw_gpu_sync(void);

/* Batch 1-D, slab-1-D, and 3-D complex double cuFFT wrappers (dir_sign: +1=GTOR, -1=RTOG). */
void vespaw_gpu_fft_1d(int dir_sign, int len, int nfft, const void *h_in, void *h_out);
void vespaw_gpu_fft_slabs(int dir_sign, int len, int nslab, int nr1l, int stride_r1,
                          int nfft, int stride_fft, int ioff, void *h_inout);
void vespaw_gpu_fft_3d(int dir_sign, int n1, int n2, int n3, int nfft,
                       const void *h_in, void *h_out);

/* Serial (ntasks == 1) batched stripe <-> 3-D real-space FFT via single 3-D cuFFT plan. */
void vespaw_gpu_fft_3d_mu(int dir_sign, int ngl, int nrl, int nr1, int nr2, int nr3,
                          int nfft, int nstripelx, const int *igtostripe,
                          const int *istripetoyz, void *h_fofg, void *h_fofr);

/* Fused serial VPSI: scatter + 3-D GTOR FFT + local potential (NDIM=1,2) + RTOG + kinetic + bucket. */
void vespaw_gpu_vpsi_serial(int ngl, int ndim, int nbh, int nrl, int nr1, int nr2, int nr3,
                            int nstripelx, int tsuper, int tbucket,
                            const int *igtostripe, const int *istripetoyz,
                            const int *minusg, const double *g2, const double *bucket,
                            const double *v, const void *h_psi, void *h_hpsi);

/* GPU interstitial charge density (serial: TSUPER, TINV, complex NDIM=1, and spinor NDIM=2). */
void vespaw_gpu_density_serial(int ngl, int ngl_full, int ndim, int nbh, int nb, int nrl,
                               int nr1, int nr2, int nr3, int nstripelx,
                               int tsuper, int tinv, const double *occ,
                               const int *igtostripe, const int *istripetoyz,
                               const int *minusg, const void *h_ei2kr,
                               const void *h_psi, double *h_rho);

/* GPU meta-GGA kinetic energy density (mode=0: TAU) and H_tau|psi> (mode=1). */
void vespaw_gpu_tau_serial(int mode, int ngl, int ndim, int nbh, int nb, int nrl,
                           int nr1, int nr2, int nr3, int nstripelx,
                           int tsuper, int tinv, const double *occ,
                           const double *gvec, const double *kvec_cart,
                           const int *igtostripe, const int *istripetoyz,
                           const int *minusg, const double *vtau,
                           const void *h_psi, double *h_tau, void *h_htaupsi);

/* Persistent GPU projector cache + GPU <P|PSI> and |PSI> += |P><P|PSI>. */
void vespaw_gpu_procache_upload(const char *gset_id, int id_len, int nat,
                                int ngl, int npro, const double *r,
                                const void *h_procache);
void vespaw_gpu_projections(int ngl, int ndim, int nb, int npro,
                            int tsuper, int ngamma, double gweight1,
                            const void *h_psi, void *h_propsi);
void vespaw_gpu_addpro(int ngl, int ndim, int nb, int npro,
                       const void *h_propsi, void *h_psi);

/* GPU scalar products, batched diagonal inner products, and generic ZGEMM. */
void vespaw_gpu_scalarproduct(int is_hermitian, int ngl, int ndim,
                              int n1, const void *h_f1,
                              int n2, const void *h_f2,
                              int tsuper, int ngamma, double gweight1,
                              void *h_prod, int *done_out);
void vespaw_gpu_dotdiag(int ngl, int ndim, int n,
                        const void *h_f1, const void *h_f2,
                        int tsuper, int ngamma, double gweight1,
                        double *h_diag, int *done_out);
void vespaw_gpu_zgemm(int transa, int transb, int m, int n, int k,
                      const void *alpha, const void *h_a, int lda,
                      const void *h_b, int ldb, const void *beta,
                      void *h_c, int ldc);

/* cuSOLVER Hermitian dense eigensolvers (ZHEGVD & ZHEEVD above cutoff N >= 64). */
void vespaw_gpu_zhegvd(int n, const void *h_h, const void *h_s,
                       double *h_e, void *h_vec, int *info_out);
void vespaw_gpu_zheevd(int n, const void *h_h, double *h_e,
                       void *h_vec, int *info_out);

/* GPU batched Teter preconditioner and full-projector construction kernel. */
void vespaw_gpu_precondition_teter(int ngl, int ndim, int nact,
                                   const double *h_g2, const double *h_ekin,
                                   const void *h_r, void *h_d, int *done_out);
void vespaw_gpu_expandpro_all(int ngl, int nat, int nsp, int lmnxx, int lnxx,
                              int lmx, int nbarepro, int npro,
                              const int *h_isp, const int *h_lnx,
                              const int *h_lmnx, const int *h_lox,
                              const double *h_r, const double *h_gvec,
                              const double *h_barepro, const double *h_ylm,
                              void *h_procache, int *done_out);

/* Device memory helpers for CUDA-aware MPI or explicit staging. */
void *vespaw_gpu_dev_alloc(size_t bytes);
void  vespaw_gpu_dev_free(void *ptr);
void  vespaw_gpu_h2d(void *d_dst, const void *h_src, size_t bytes);
void  vespaw_gpu_d2h(void *h_dst, const void *d_src, size_t bytes);

#ifdef __cplusplus
}
#endif

#endif /* VESPAW_GPU_BACKEND_H */
