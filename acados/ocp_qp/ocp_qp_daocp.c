/*
 * Copyright (c) The acados authors.
 *
 * This file is part of acados.
 *
 * The 2-Clause BSD License
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *
 * 1. Redistributions of source code must retain the above copyright notice,
 * this list of conditions and the following disclaimer.
 *
 * 2. Redistributions in binary form must reproduce the above copyright notice,
 * this list of conditions and the following disclaimer in the documentation
 * and/or other materials provided with the distribution.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 * AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 * ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
 * LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
 * CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
 * SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
 * INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
 * CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
 * ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 * POSSIBILITY OF SUCH DAMAGE.;
 */


#include <assert.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

// daocp
#include "daocp/include/daocp.h"
#include "daocp/include/internal.h"

// blasfeo
#include "blasfeo/include/blasfeo.h"

// acados
#include "acados/ocp_qp/ocp_qp_common.h"
#include "acados/ocp_qp/ocp_qp_daocp.h"
#include "acados/utils/mem.h"
#include "acados/utils/print.h"
#include "acados/utils/timing.h"
#include "acados/utils/types.h"


/************************************************
 * opts
 ************************************************/

acados_size_t ocp_qp_daocp_opts_calculate_size(void *config_, void *dims_)
{
    return sizeof(ocp_qp_daocp_opts);
}


void *ocp_qp_daocp_opts_assign(void *config_, void *dims_, void *raw_memory)
{
    return raw_memory;
}


void ocp_qp_daocp_opts_initialize_default(void *config_, void *dims_, void *opts_)
{
    ocp_qp_daocp_opts* opts = opts_;
    opts->print_level = 0;
    opts->warm_start = 0;
    opts->first_run = 1;
    opts->equality_detection_tolerance = 0.0;
    daocp_args_set_default(&opts->daocp_opts);
    return;
}



void ocp_qp_daocp_opts_update(void *config_, void *dims_, void *opts_)
{
    return;
}

void ocp_qp_daocp_opts_set(void *config_, void *opts_, const char *field, void *value)
{
    ocp_qp_daocp_opts* opts = opts_;

    if (!strcmp(field, "max_iter") || !strcmp(field, "iter_max"))
    {
        int *tmp_ptr = value;
        opts->daocp_opts.max_iter = *tmp_ptr;
    }
    else if (!strcmp(field, "print_level"))
    {
        int* print_level = (int *) value;
        opts->print_level = *print_level;
    }
    else if (!strcmp(field, "warm_start"))
    {
        int* warm_start = (int *) value;
        opts->warm_start = *warm_start;
    }
    else if (!strcmp(field, "tol_ineq"))
    {
        double *tol = value;
        opts->daocp_opts.primal_tol = *tol;
    }
    else if (!strcmp(field, "selection_strategy"))
    {
        int code = *((int*) value);
        if (code == 0) opts->daocp_opts.selection = DAOCP_SELECT_GREEDY;
        else if (code == 1) opts->daocp_opts.selection = DAOCP_SELECT_MOST_VIOLATED;
        else {
            printf("\nerror: ocp_qp_daocp_opts_set: invalid selection_strategy %d\n", code);
            exit(1);
        }
    }
    else if (!strcmp(field, "equality_detection_tolerance"))
    {
        opts->equality_detection_tolerance = *((double*) value);
        if (opts->equality_detection_tolerance < 0.0)
        {
            printf("\nerror: equality_detection_tolerance must be nonnegative\n");
            exit(1);
        }
    }
    else if (!strcmp(field, "tol_stat") || !strcmp(field, "tol_eq")
             || !strcmp(field, "tol_ineq") || !strcmp(field, "tol_comp"))
    {
        /* DAOCP does not expose termination tolerances. These fields are
         * nevertheless set unconditionally by the OCP NLP option layer. */
    }
    else
    {
        printf("\nWARNING: ocp_qp_daocp_opts_set: field: %s not interfaced yet. Ignoring option and \n", field);
        exit(1);
    }

    return;
}

void ocp_qp_daocp_opts_get(void *config_, void *opts_, const char *field, void *value)
{
    ocp_qp_daocp_opts* opts = opts_;

    if (!strcmp(field, "max_iter") || !strcmp(field, "iter_max"))
        *((int*) value) = opts->daocp_opts.max_iter;
    else if (!strcmp(field, "print_level"))
        *((int*) value) = opts->print_level;
    else if (!strcmp(field, "warm_start"))
        *((int*) value) = opts->warm_start;
    else if (!strcmp(field, "selection_strategy"))
        *((int*) value) = opts->daocp_opts.selection;
    else if (!strcmp(field, "equality_detection_tolerance"))
        *((double*) value) = opts->equality_detection_tolerance;
    else
    {
        printf("\nerror: ocp_qp_daocp_opts_get: field %s not available\n", field);
        exit(1);
    }
}




/************************************************
 * memory
 ************************************************/

acados_size_t ocp_qp_daocp_memory_calculate_size(void *config_, void *dims_, void *opts_)
{
    ocp_qp_dims *dims = dims_;
    acados_size_t size = sizeof(ocp_qp_daocp_memory);
    size += sizeof(daocp_workspace);

    int N = dims->N;
    int* nx = dims->nx;
    int* nu = dims->nu;
    int* nb = dims->nb;
    int* nbu = dims->nbu;
    int* nbx = dims->nbx;
    int* ng = dims->ng;

    // daocp_qp data
    size += 4*(N+1)*sizeof(u32); // dims.ng, dims.ne, dims.nbu, dims.nbx
    size += (3*N+2)*sizeof(u32); // mem.nge, mem.nbue, mem.nbxe
    size += (3*N+2)*sizeof(u32*); // mem.idxge, mem.idxbue, mem.idxbxe
    for (u32 t=0; t<=N; ++t) {
        if (t<N) size += nbu[t]*sizeof(u32); // mem.idxbue[t]
        if (t>0) size += nbx[t]*sizeof(u32); // mem.idxbxe[t]
        size += ng[t]*sizeof(u32); // mem.idxge[t]
    }
    size += nx[0]*sizeof(f64); // x0
    u32 neq = 0;
    u32 nin = 0;
    u32 nb_tot = 0;
    for (u32 t=0; t<=N; ++t) {
        u32 ne_capacity_t = ng[t]+nbu[t]+nbx[t]*(t>0 ? 1 : 0);
        size += 2*nb[t]*sizeof(f64); // lbu/x[t], ubu/x[t]
        size += 2*nb[t]*sizeof(f64); // lbu/x_wrk[t], ubu/x_wrk[t]
        size += nb[t]*sizeof(u32); // idxbu[t]
        size += 2*ng[t]*sizeof(f64); // lg[t], ug[t]
        size += 2*ng[t]*sizeof(f64); // lg_wrk[t], ug_wrk[t]
        size += ne_capacity_t*sizeof(f64); // d[t]
        size += ng[t]*(nu[t]+nx[t])*sizeof(f64); // Cu[t], Cx[t]
        size += (t<N ? ne_capacity_t : 0)*nu[t]*sizeof(f64); // Du[t]
        size += ne_capacity_t*nx[t]*sizeof(f64); // Dx[t]

        neq += ne_capacity_t;
        nin += nb[t]+ng[t];
        nb_tot += nb[t];
    }
    size -= ng[0]*nx[0]*sizeof(f64); // Subtract off Cx[0]
    size += (2*N+1)*sizeof(f64*); // Cu[:], Cx[:]
    size += 2*(N+1)*sizeof(f64*); // Du[:], Dx[:]
    size += 4*(2*N+1)*sizeof(f64*); // lbu/x[:], ubu/x[:], lbu/x_wrk[:], ubu/x_wrk[:]
    size += 4*(N+1)*sizeof(f64*); // lg[:], ug[:], lg_wrk[:], ug_wrk[:]
    size += (N+1)*sizeof(f64*); // d[:]
    size += (2*N+1)*sizeof(u32*); // idxbu[:], idxbx[:]

    // daocp_workspace data
    int max_nx = 0; int max_nu = 0;
    for (u32 t=0; t<=N; ++t) if (nx[t] > max_nx) max_nx = nx[t];
    for (u32 t=0; t<N; ++t) if (nu[t] > max_nu) max_nu = nu[t];
    int tot_nx = nx[N];
    int tot_nu = 0;
    for (u32 t=0; t<N; ++t) {
        size += blasfeo_memsize_dmat(nx[t+1], nx[t+1]); // P
        size += blasfeo_memsize_dmat(nx[t], nu[t]); // Ku
        size += blasfeo_memsize_dmat(nx[t], nu[t]); // Ke
        size += blasfeo_memsize_dmat(nu[t], nu[t]); // Luu
        size += blasfeo_memsize_dmat(nu[t], nu[t]); // Lue
        size += blasfeo_memsize_dmat(nu[t], nu[t]); // Lee

        size += blasfeo_memsize_dvec(nu[t]+nx[t]); // ux_lqr
        size += blasfeo_memsize_dvec(nu[t]); // eta_lqr
        size += blasfeo_memsize_dvec(nu[t]); // b

        tot_nx += nx[t]; tot_nu += nu[t];
    }
    size += blasfeo_memsize_dvec(nx[N]); // ux_lqr[N]
    size += 6*N*sizeof(struct blasfeo_dmat); // P, Ku, Ke, Luu, Lue, Lee
    size += (3*N+1)*sizeof(struct blasfeo_dvec); // ux_lqr, eta_lqr, b
    size += (2*tot_nu + tot_nx - nx[0])*sizeof(f64); // x, u, eta
    size += (5*N+3)*sizeof(f64*); // x[:], u[:], eta[:], sol.lam[:], sol.eta[:]
    size += (N+1)*sizeof(daocp_constraint_type*); // contypes[:]
    size += (nin-nb_tot)*sizeof(daocp_constraint_type); // contypes
    size += 4*N*sizeof(u32); // cnu, rho, crho, m
    size += (N+1)*sizeof(u32*); // as.constraint_status[:]
    size += N*sizeof(f64*); // J
    size += 2*nin*sizeof(u32); // as.constraint_status
        
    int W_stride = nin+1;
    size += W_stride*sizeof(daocp_constraint); // as.xi2con
    size += 6*W_stride*sizeof(f64); // xi, xis, p, ps, dual_linear, dual_intermediate
    size += W_stride*sizeof(u32); // xi_sign
    size += (W_stride+1)*W_stride*sizeof(f64); // Ld
    size += 2*W_stride*tot_nu*sizeof(f64); // Mu, Me

    // We must allocate nin slots for equality duals in the solution
    // because across solves any inequality constraint may become an equality.
    size += nin*sizeof(f64); // sol.eta, sol.lam data
    // The upper bound on J data can be tighter (see daocp for remark on the topic)
    size += nin*tot_nu*sizeof(f64); // J data
    size += neq*max_nx*sizeof(f64); // H
    size += 2*neq*sizeof(f64); // h, tmp1
    size += (neq+1)*(max_nx+max_nu+1)*sizeof(f64); // GEtmp
    size += max_nx*DAOCP_MAX(max_nx, max_nu)*sizeof(f64); // ABtmp
    size += DAOCP_MAX(blasfeo_memsize_dmat(max_nx, max_nx), blasfeo_memsize_dmat(max_nu, max_nx)); // tmp2
    size += 2*blasfeo_memsize_dvec(DAOCP_MAX(max_nx, max_nu)); // costate0, costate1

    // Symmetric slack penalties for every possible inequality.
    size += 2*(3*N+2)*sizeof(f64*); // Zbu, Zbx, Zg, zbu, zbx, zg
    size += 2*nin*sizeof(f64);

    // Worst-case padding for double-alignment points and two cache-line
    // aligned BLASFEO backing-store regions.
    size += 6*8 + 2*64;

    // The enclosing xcond solver places its workspace immediately after this
    // memory block and requires the workspace structure to stay aligned.
    make_int_multiple_of(8, &size);

    return size;
}

static inline void* assign_ptr_vec(void* ptr_vec_, char* ptr, int* size_vec1, int* size_vec2, int* size_vec3, int unit_size, int n)
{
    char** ptr_vec = (char**) ptr_vec_;
    for (u32 i=0; i<n; ++i) {
        ptr_vec[i] = ptr;
        ptr += (size_vec1[i] - 
            (size_vec2 ? size_vec2[i] : 0) +
            (size_vec3 ? size_vec3[i] : 0))*unit_size;
    }
    return ptr;
}

static inline void* assign_ptr_mat(
    void* ptr_mat_, char* ptr,
    int* size_vec1, int* size_vec2, int* size_vec3, int* size_vec4, 
    int* size_vec5, int unit_size, int n)
{
    char** ptr_mat = (char**) ptr_mat_;
    for (u32 i=0; i<n; ++i) {
        ptr_mat[i] = ptr;
        ptr += (size_vec1[i] - 
            (size_vec2 ? size_vec2[i] : 0) +
            (size_vec3 ? size_vec3[i] : 0) +
            (size_vec4 ? size_vec4[i] : 0))*size_vec5[i]*unit_size;
    }
    return ptr; 
}

void *ocp_qp_daocp_memory_assign(void *config_, void *dims_, void *opts_, void *raw_memory)
{
    char* c_ptr = (char*) raw_memory;
    ocp_qp_dims* dims = dims_;
    ocp_qp_daocp_memory* mem = (ocp_qp_daocp_memory*) c_ptr;
    c_ptr += sizeof(ocp_qp_daocp_memory);
    mem->workspace = c_ptr;
    daocp_workspace* wrk = mem->workspace;

    int N = dims->N;
    int* nx = dims->nx;
    int* nu = dims->nu;
    int* nb = dims->nb;
    int* nbu = dims->nbu;
    int* nbx = dims->nbx;
    int* ng = dims->ng;

    int max_nx = 0; int max_nu = 0; int neq = 0; int nin = 0; int tot_nu = 0;
    for (u32 t=0; t<=N; ++t) if (nx[t] > max_nx) max_nx = nx[t];
    for (u32 t=0; t<N; ++t) {
        if (nu[t] > max_nu) max_nu = nu[t];
        tot_nu += nu[t];
    }
    for (u32 t=0; t<=N; ++t) neq += ng[t] + nbu[t] + nbx[t];
    neq -= nbx[0];
    for (u32 t=0; t<=N; ++t) nin += ng[t]+nb[t];

    // daocp_qp data
    mem->qp.dims.nbu = (u32*) (c_ptr += sizeof(daocp_workspace));
    mem->qp.dims.nbx = (u32*) (c_ptr += (N+1)*sizeof(u32));
    mem->qp.dims.ng = (u32*) (c_ptr += (N+1)*sizeof(u32));
    mem->qp.dims.ne = (u32*) (c_ptr += (N+1)*sizeof(u32));
    mem->nbue = (u32*) (c_ptr += (N+1)*sizeof(u32));
    mem->nbxe = (u32*) (c_ptr += N*sizeof(u32));
    mem->nge = (u32*) (c_ptr += (N+1)*sizeof(u32));
    c_ptr += (N+1)*sizeof(u32);
    align_char_to(8, &c_ptr);
    mem->idxbue = (u32**) c_ptr;
    mem->idxbxe = (u32**) (c_ptr += N*sizeof(u32*));
    mem->idxge = (u32**) (c_ptr += (N+1)*sizeof(u32*));
    c_ptr += (N+1)*sizeof(u32**);
    c_ptr = assign_ptr_vec(mem->idxbue, c_ptr, nbu, 0, 0, sizeof(u32), N);
    mem->idxbxe[0] = 0;
    c_ptr = assign_ptr_vec(mem->idxbxe+1, c_ptr, nbx+1, 0, 0, sizeof(u32), N);
    c_ptr = assign_ptr_vec(mem->idxge, c_ptr, ng, 0, 0, sizeof(u32), N+1);
    // Set number of detected equalities to 0
    memset(mem->nbue, 0, N*sizeof(u32));
    memset(mem->nbxe, 0, (N+1)*sizeof(u32));
    memset(mem->nge, 0, (N+1)*sizeof(u32));
    align_char_to(8, &c_ptr);
    mem->qp.x0 = (f64*) c_ptr;
    mem->qp.lbu = (f64**) (c_ptr += nx[0]*sizeof(f64));
    mem->qp.lbx = (f64**) (c_ptr += N*sizeof(f64*));
    wrk->lbu_wrk = (f64**) (c_ptr += (N+1)*sizeof(f64*));
    wrk->lbx_wrk = (f64**) (c_ptr += N*sizeof(f64*));
    mem->qp.ubu = (f64**) (c_ptr += (N+1)*sizeof(f64*));
    mem->qp.ubx = (f64**) (c_ptr += N*sizeof(f64*));
    wrk->ubu_wrk = (f64**) (c_ptr += (N+1)*sizeof(f64*));
    wrk->ubx_wrk = (f64**) (c_ptr += N*sizeof(f64*));
    mem->qp.cl = (f64**) (c_ptr += (N+1)*sizeof(f64*));
    mem->qp.cu = (f64**) (c_ptr += (N+1)*sizeof(f64*));
    wrk->lg_wrk = (f64**) (c_ptr += (N+1)*sizeof(f64*));
    wrk->ug_wrk = (f64**) (c_ptr += (N+1)*sizeof(f64*));
    mem->qp.d = (f64**) (c_ptr += (N+1)*sizeof(f64*));
    mem->qp.Cu = (f64**) (c_ptr += (N+1)*sizeof(f64*));
    mem->qp.Cx = (f64**) (c_ptr += N*sizeof(f64*));
    mem->qp.Du = (f64**) (c_ptr += (N+1)*sizeof(f64*));
    mem->qp.Dx = (f64**) (c_ptr += (N+1)*sizeof(f64*));
    mem->sol.lam = (f64**) (c_ptr += (N+1)*sizeof(f64*));
    mem->sol.eta = (f64**) (c_ptr += (N+1)*sizeof(f64*));
    mem->qp.idxbu = (u32**) (c_ptr += (N+1)*sizeof(f64*));
    mem->qp.idxbx = (u32**) (c_ptr += N*sizeof(u32*));
    c_ptr += (N+1)*sizeof(u32*);
    align_char_to(8, &c_ptr);
    // Only the first pointer in sol.lam is assigned, at allocation time.
    mem->sol.lam[0] = (f64*) c_ptr;
    c_ptr += nin*sizeof(f64);
    c_ptr = assign_ptr_vec(mem->qp.lbu, c_ptr, nbu, 0, 0, sizeof(f64), N);
    c_ptr = assign_ptr_vec(wrk->lbu_wrk, c_ptr, nbu, 0, 0, sizeof(f64), N);
    mem->qp.lbx[0] = 0; c_ptr = assign_ptr_vec(mem->qp.lbx+1, c_ptr, nbx+1, 0, 0, sizeof(f64), N);
    wrk->lbx_wrk[0] = 0; c_ptr = assign_ptr_vec(wrk->lbx_wrk+1, c_ptr, nbx+1, 0, 0, sizeof(f64), N);
    c_ptr = assign_ptr_vec(mem->qp.ubu, c_ptr, nbu, 0, 0, sizeof(f64), N);
    c_ptr = assign_ptr_vec(wrk->ubu_wrk, c_ptr, nbu, 0, 0, sizeof(f64), N);
    mem->qp.ubx[0] = 0; c_ptr = assign_ptr_vec(mem->qp.ubx+1, c_ptr, nbx+1, 0, 0, sizeof(f64), N);
    wrk->ubx_wrk[0] = 0; c_ptr = assign_ptr_vec(wrk->ubx_wrk+1, c_ptr, nbx+1, 0, 0, sizeof(f64), N);
    c_ptr = assign_ptr_vec(mem->qp.idxbu, c_ptr, nbu, 0, 0, sizeof(u32), N);
    mem->qp.idxbx[0] = 0; c_ptr = assign_ptr_vec(mem->qp.idxbx+1, c_ptr, nbx+1, 0, 0, sizeof(u32), N);
    align_char_to(8, &c_ptr);
    c_ptr = assign_ptr_vec(mem->qp.cl, c_ptr, ng, 0, 0, sizeof(f64), N+1);
    c_ptr = assign_ptr_vec(wrk->lg_wrk, c_ptr, ng, 0, 0, sizeof(f64), N+1);
    c_ptr = assign_ptr_vec(mem->qp.cu, c_ptr, ng, 0, 0, sizeof(f64), N+1);
    c_ptr = assign_ptr_vec(wrk->ug_wrk, c_ptr, ng, 0, 0, sizeof(f64), N+1);
    mem->qp.d[0] = (f64*) c_ptr; c_ptr += (nbu[0]+ng[0])*sizeof(f64);
    for (u32 t=1; t<=N; ++t) {
        mem->qp.d[t] = (f64*) c_ptr;
        c_ptr += (nbu[t]+nbx[t]+ng[t])*sizeof(f64);
    }
    c_ptr = assign_ptr_mat(mem->qp.Cu, c_ptr, ng, 0, 0, 0, nu, sizeof(f64), N);
    mem->qp.Cx[0]=0; c_ptr = assign_ptr_mat(mem->qp.Cx+1, c_ptr, ng+1, 0, 0, 0, nx+1, sizeof(f64), N);
    mem->qp.Du[0] = (f64*) c_ptr; c_ptr += (ng[0]+nbu[0])*nu[0]*sizeof(f64);
    for (u32 t=1; t<N; ++t) {
        mem->qp.Du[t] = (f64*) c_ptr;
        c_ptr += (ng[t]+nbu[t]+nbx[t])*nu[t]*sizeof(f64);
    }
    mem->qp.Du[N] = (f64*) c_ptr;
    mem->qp.Dx[0]= (f64*) c_ptr; c_ptr += (ng[0]+nbu[0])*nx[0]*sizeof(f64);
    for (u32 t=1; t<=N; ++t) {
        mem->qp.Dx[t] = (f64*) c_ptr;
        c_ptr += (ng[t]+nbu[t]+nbx[t])*nx[t]*sizeof(f64);
    }

    // One pair of penalties per constraint, shared by both bounds.
    mem->qp.Zbu = (f64**) c_ptr; c_ptr += N*sizeof(f64*);
    mem->qp.zbu = (f64**) c_ptr; c_ptr += N*sizeof(f64*);
    mem->qp.Zbx = (f64**) c_ptr; c_ptr += (N+1)*sizeof(f64*);
    mem->qp.zbx = (f64**) c_ptr; c_ptr += (N+1)*sizeof(f64*);
    mem->qp.Zg = (f64**) c_ptr; c_ptr += (N+1)*sizeof(f64*);
    mem->qp.zg = (f64**) c_ptr; c_ptr += (N+1)*sizeof(f64*);
    c_ptr = assign_ptr_vec(mem->qp.Zbu, c_ptr, nbu, 0, 0, sizeof(f64), N);
    c_ptr = assign_ptr_vec(mem->qp.zbu, c_ptr, nbu, 0, 0, sizeof(f64), N);
    mem->qp.Zbx[0] = mem->qp.zbx[0] = NULL;
    c_ptr = assign_ptr_vec(mem->qp.Zbx+1, c_ptr, nbx+1, 0, 0, sizeof(f64), N);
    c_ptr = assign_ptr_vec(mem->qp.zbx+1, c_ptr, nbx+1, 0, 0, sizeof(f64), N);
    c_ptr = assign_ptr_vec(mem->qp.Zg, c_ptr, ng, 0, 0, sizeof(f64), N+1);
    c_ptr = assign_ptr_vec(mem->qp.zg, c_ptr, ng, 0, 0, sizeof(f64), N+1);
    // Set inverse quadratic penalties to zero: default initialization is
    // all constraints are hard.
    for (u32 t=0; t<N; ++t) memset(mem->qp.Zbu[t], 0, nbu[t]*sizeof(f64)); 
    for (u32 t=1; t<=N; ++t) memset(mem->qp.Zbx[t], 0, nbx[t]*sizeof(f64)); 
    for (u32 t=0; t<=N; ++t) memset(mem->qp.Zg[t], 0, ng[t]*sizeof(f64)); 

    // daocp_workspace data
    wrk->P = (struct blasfeo_dmat*) c_ptr; c_ptr+=N*sizeof(struct blasfeo_dmat);
    wrk->Ku = (struct blasfeo_dmat*) c_ptr; c_ptr+=N*sizeof(struct blasfeo_dmat);
    wrk->Ke = (struct blasfeo_dmat*) c_ptr; c_ptr+=N*sizeof(struct blasfeo_dmat);
    wrk->Luu = (struct blasfeo_dmat*) c_ptr; c_ptr+=N*sizeof(struct blasfeo_dmat);
    wrk->Lue = (struct blasfeo_dmat*) c_ptr; c_ptr+=N*sizeof(struct blasfeo_dmat);
    wrk->Lee = (struct blasfeo_dmat*) c_ptr; c_ptr+=N*sizeof(struct blasfeo_dmat);
    wrk->ux_lqr = (struct blasfeo_dvec*) c_ptr; c_ptr+=(N+1)*sizeof(struct blasfeo_dvec);
    wrk->eta_lqr = (struct blasfeo_dvec*) c_ptr; c_ptr+=N*sizeof(struct blasfeo_dvec);
    wrk->b = (struct blasfeo_dvec*) c_ptr; c_ptr+=N*sizeof(struct blasfeo_dvec);

    wrk->u = (f64**) c_ptr; c_ptr += N*sizeof(f64*);
    wrk->x = (f64**) c_ptr; c_ptr += (N+1)*sizeof(f64*);
    wrk->eta = (f64**) c_ptr; c_ptr += N*sizeof(f64*);
    wrk->J = (f64**) c_ptr; c_ptr += N*sizeof(f64*);
    wrk->contypes = (daocp_constraint_type**) c_ptr; c_ptr += (N+1)*sizeof(daocp_constraint_type*);
    wrk->as.constraint_status = (u32**) c_ptr; c_ptr += (N+1)*sizeof(u32*);
    c_ptr = assign_ptr_vec(wrk->u, c_ptr, nu, 0, 0, sizeof(f64), N);
    wrk->x[0]=0; c_ptr = assign_ptr_vec(wrk->x+1, c_ptr, nx+1, 0, 0, sizeof(f64), N);
    if (N > 0) wrk->eta[0] = (f64*) c_ptr;
    c_ptr += tot_nu*sizeof(f64);
    for (u32 t=0; t<=N; ++t) {
        wrk->contypes[t] = (daocp_constraint_type*) c_ptr;
        c_ptr += ng[t]*sizeof(daocp_constraint_type);
    }
    for (u32 t=0; t<=N; ++t) {
        wrk->as.constraint_status[t] = (u32*) c_ptr;
        c_ptr += 2*(ng[t]+nb[t])*sizeof(u32);
    }
    wrk->cnu = (u32*) c_ptr; c_ptr += N*sizeof(u32);
    wrk->rho = (u32*) c_ptr; c_ptr += N*sizeof(u32);
    wrk->crho = (u32*) c_ptr; c_ptr += N*sizeof(u32);
    wrk->m = (u32*) c_ptr; c_ptr += N*sizeof(u32);

    align_char_to(64, &c_ptr);
    for (u32 t=0; t<N; ++t) {
        blasfeo_create_dmat(nx[t+1], nx[t+1], wrk->P+t, c_ptr);
        c_ptr += blasfeo_memsize_dmat(nx[t+1], nx[t+1]);
        blasfeo_create_dmat(nx[t], nu[t], wrk->Ku+t, c_ptr);
        c_ptr += blasfeo_memsize_dmat(nx[t], nu[t]);
        blasfeo_create_dmat(nx[t], nu[t], wrk->Ke+t, c_ptr);
        c_ptr += blasfeo_memsize_dmat(nx[t], nu[t]);
        blasfeo_create_dmat(nu[t], nu[t], wrk->Luu+t, c_ptr);
        c_ptr += blasfeo_memsize_dmat(nu[t], nu[t]);
        blasfeo_create_dmat(nu[t], nu[t], wrk->Lue+t, c_ptr);
        c_ptr += blasfeo_memsize_dmat(nu[t], nu[t]);
        blasfeo_create_dmat(nu[t], nu[t], wrk->Lee+t, c_ptr);
        c_ptr += blasfeo_memsize_dmat(nu[t], nu[t]);

        blasfeo_create_dvec(nu[t]+nx[t], wrk->ux_lqr+t, c_ptr);
        c_ptr += blasfeo_memsize_dvec(nu[t]+nx[t]);
        blasfeo_create_dvec(nu[t], wrk->eta_lqr+t, c_ptr);
        c_ptr += blasfeo_memsize_dvec(nu[t]);
        blasfeo_create_dvec(nu[t], wrk->b+t, c_ptr);
        c_ptr += blasfeo_memsize_dvec(nu[t]);
    }
    blasfeo_create_dvec(nx[N], wrk->ux_lqr+N, c_ptr);
    c_ptr += blasfeo_memsize_dvec(nx[N]);
    
    int W_stride = nin + 1;
    wrk->as.xi2con = (daocp_constraint*) c_ptr; c_ptr += W_stride*sizeof(daocp_constraint);
    align_char_to(8, &c_ptr);
    wrk->xi = (f64*) c_ptr; c_ptr+=W_stride*sizeof(f64);
    wrk->xis = (f64*) c_ptr; c_ptr+=W_stride*sizeof(f64);
    wrk->ps = (f64*) c_ptr; c_ptr+=W_stride*sizeof(f64);
    wrk->dual_intermediate = (f64*) c_ptr; c_ptr+=W_stride*sizeof(f64);
    wrk->p = (f64*) c_ptr; c_ptr+=W_stride*sizeof(f64);
    wrk->dual_linear = (f64*) c_ptr; c_ptr+=W_stride*sizeof(f64);
    wrk->xi_sign = (u32*) c_ptr; c_ptr+=W_stride*sizeof(u32);
    align_char_to(8, &c_ptr);
    wrk->Ld = (f64*) c_ptr; c_ptr+=(W_stride+1)*W_stride*sizeof(f64);
    wrk->Mu = (f64*) c_ptr; c_ptr+=W_stride*tot_nu*sizeof(f64);
    wrk->Me = (f64*) c_ptr; c_ptr+=W_stride*tot_nu*sizeof(f64);
    // Only the last J pointer is assigned, at allocation time.
    wrk->J[N-1] = (f64*) c_ptr; c_ptr+=tot_nu*nin*sizeof(f64);
    wrk->H = (f64*) c_ptr; c_ptr+=neq*max_nx*sizeof(f64);
    wrk->h = (f64*) c_ptr; c_ptr+=neq*sizeof(f64);
    wrk->tmp1 = (f64*) c_ptr; c_ptr+=neq*sizeof(f64);
    wrk->GEtmp = (f64*) c_ptr; c_ptr+=(neq+1)*(max_nx+max_nu+1)*sizeof(f64);
    wrk->ABtmp = (f64*) c_ptr; c_ptr+=max_nx*DAOCP_MAX(max_nx, max_nu)*sizeof(f64);
    align_char_to(64, &c_ptr);
    blasfeo_create_dmat(DAOCP_MAX(max_nx, max_nu), max_nx, &wrk->tmp2, c_ptr);
    c_ptr += blasfeo_memsize_dmat(DAOCP_MAX(max_nu, max_nx), max_nx);
    int max_nx_nu = DAOCP_MAX(max_nx, max_nu);
    blasfeo_create_dvec(max_nx_nu, &wrk->costate0, c_ptr);
    c_ptr += blasfeo_memsize_dvec(max_nx_nu);
    blasfeo_create_dvec(max_nx_nu, &wrk->costate1, c_ptr);
    c_ptr += blasfeo_memsize_dvec(max_nx_nu);

    ((ocp_qp_daocp_opts*) opts_)->first_run = 1;

    assert((char *) raw_memory + ocp_qp_daocp_memory_calculate_size(config_, dims, opts_) >= c_ptr);

    return mem;
}


void ocp_qp_daocp_memory_get(void *config_, void *mem_, const char *field, void* value)
{
    // qp_solver_config *config = config_;
    ocp_qp_daocp_memory *mem = mem_;

    if (!strcmp(field, "time_qp_solver_call"))
    {
        *((double*) value) = mem->time_qp_solver_call;
    }
    else if (!strcmp(field, "iter"))
    {
        *((int*) value) = ((daocp_workspace*) mem->workspace)->iters;
    }
    else if (!strcmp(field, "status"))
    {
        *((int*) value) = ((daocp_workspace*) mem->workspace)->status;
    }
    else
    {
        printf("\nerror: ocp_qp_daocp_memory_get: field %s not available\n", field);
        exit(1);
    }

    return;
}


void ocp_qp_daocp_memory_reset(void *config_, void *qp_in_, void *qp_out_, void *opts_, void *mem_, void *work_)
{
    ocp_qp_daocp_opts* opts = opts_;
    opts->first_run = 1;
}



/************************************************
 * workspace
 ************************************************/

acados_size_t ocp_qp_daocp_workspace_calculate_size(void *config_, void *dims_, void *opts_)
{
    return 0;
}

static u32 acados_daocp_contains_index(const int *indices, int offset, int count, int index)
{
    for (int i = 0; i < count; i++)
        if (indices[offset+i] == index)
            return 1;

    return 0;
}

/* DAOCP stores inverse quadratic penalties; zero denotes a hard constraint. */
static u32 acados_daocp_is_soft_constraint(const ocp_qp_in* qp, u32 t, u32 index)
{
    int slack = qp->idxs_rev[t][index];
    return slack >= 0 && (BLASFEO_DVECEL(qp->Z+t, slack) > 0.0
        || BLASFEO_DVECEL(qp->Z+t, qp->dim->ns[t]+slack) > 0.0);
}

static void acados_daocp_populate_slack_penalties(
    const ocp_qp_in* qp, u32 t, u32 index, f64* Z, f64* z)
{
    int slack = qp->idxs_rev[t][index];
    *Z = 0.0;
    *z = 0.0;
    if (slack >= 0 && BLASFEO_DVECEL(qp->Z+t, slack) > 0.0)
    {
        *Z = 1.0 / BLASFEO_DVECEL(qp->Z+t, slack);
        *z = BLASFEO_DVECEL(qp->rqz+t, qp->dim->nu[t]+qp->dim->nx[t]+slack);
    }
}

/* DAOCP has a single penalty per constraint, with nonnegative slacks.
 * Reject acados formulations that cannot be represented by these fields. */
static int acados_daocp_validate_slacks(const ocp_qp_in* qp)
{
    const ocp_qp_dims* dims = qp->dim;
    for (int t=0; t<=dims->N; ++t)
    {
        int ns = dims->ns[t], nb = dims->nb[t], ng = dims->ng[t];
        for (int i=0; i<nb+ng; ++i)
            if (qp->idxs_rev[t][i] >= 0 &&
                (BLASFEO_DVECEL(qp->d+t, i) > -BLASFEO_DVECEL(qp->d+t, nb+ng+i)
                 || BLASFEO_DVECEL(qp->d_mask+t, i) != 1.0
                 || BLASFEO_DVECEL(qp->d_mask+t, nb+ng+i) != 1.0))
            {
                printf("DAOCP requires ordered, unmasked bounds on soft constraints (stage %d, constraint %d).\n", t, i);
                return ACADOS_QP_FAILURE;
            }
        for (int k=0; k<ns; ++k)
        {
            double Zl = BLASFEO_DVECEL(qp->Z+t, k);
            double Zu = BLASFEO_DVECEL(qp->Z+t, ns+k);
            double zl = BLASFEO_DVECEL(qp->rqz+t, dims->nu[t]+dims->nx[t]+k);
            double zu = BLASFEO_DVECEL(qp->rqz+t, dims->nu[t]+dims->nx[t]+ns+k);
            if (!isfinite(Zl) || Zl <= 0.0 || !isfinite(1.0/Zl) || Zl != Zu
                || !isfinite(zl) || zl < 0.0 || zl != zu)
            {
                printf("DAOCP requires equal lower/upper slack penalties, Z > 0 and z >= 0 (stage %d, slack %d).\n", t, k);
                return ACADOS_QP_FAILURE;
            }
            if (BLASFEO_DVECEL(qp->d+t, 2*(nb+ng)+k) != 0.0
                || BLASFEO_DVECEL(qp->d+t, 2*(nb+ng)+ns+k) != 0.0)
            {
                printf("DAOCP requires zero lower bounds on slacks (stage %d, slack %d).\n", t, k);
                return ACADOS_QP_FAILURE;
            }
            int count = 0;
            for (int i=0; i<nb+ng; ++i)
                count += qp->idxs_rev[t][i] == k;
            if (count != 1)
            {
                printf("DAOCP requires each slack to soften exactly one constraint (stage %d, slack %d).\n", t, k);
                return ACADOS_QP_FAILURE;
            }
        }
        // x0 is a fixed parameter in DAOCP, so its fixing bounds must remain hard.
        if (t == 0)
            for (int i=0; i<dims->nbxe[0]; ++i)
                if (acados_daocp_is_soft_constraint(qp, 0, qp->idxe[0][dims->nbue[0]+i]))
                {
                    printf("DAOCP does not support softened initial-state fixing bounds.\n");
                    return ACADOS_QP_FAILURE;
                }
    }
    return ACADOS_SUCCESS;
}

static u32 acados_daocp_is_equality(const ocp_qp_daocp_memory* mem, u32 stage, u32 idx,
                                    daocp_constraint_type contype)
{
    u32* ncon;
    u32** con_idxs;
    switch (contype) {
    case DAOCP_BOUND_U:
        ncon = mem->nbue;
        con_idxs = mem->idxbue;
        break;
    case DAOCP_BOUND_X:
        ncon = mem->nbxe;
        con_idxs = mem->idxbxe;
        break;
    default:
        ncon = mem->nge;
        con_idxs = mem->idxge;
        break;
    }
    for (u32 i=0; i<ncon[stage]; ++i) if (con_idxs[stage][i] == idx) return 1;
    return 0;
}

static u32 acados_daocp_detect_equalities(
    const ocp_qp_in* qp_in, const ocp_qp_daocp_opts* opts,
    ocp_qp_daocp_memory* mem)
{
    const ocp_qp_dims* dims = qp_in->dim;
    u32 pattern_changed = 0;
    u32 old_count, detected_count;

    for (u32 t=0; t<=(u32) dims->N; ++t)
    {   
        // INPUT BOUNDS
        if (t < dims->N) {
        old_count = mem->nbue[t];
        detected_count = 0;
        for (u32 i=0; i<(u32) dims->nbu[t]; ++i)
        {
            if (acados_daocp_is_soft_constraint(qp_in, t, i)) continue;
            u32 explicit_equality = acados_daocp_contains_index(
                qp_in->idxe[t], 0, dims->nbue[t], (int) i);
            f64 lower = BLASFEO_DVECEL(qp_in->d+t, i);
            f64 upper = -BLASFEO_DVECEL(
                qp_in->d+t, dims->ng[t]+dims->nb[t]+i);
            if (explicit_equality || (DAOCP_ABS(upper-lower) <= opts->equality_detection_tolerance))
            {
                if (detected_count >= old_count
                    || mem->idxbue[t][detected_count] != i)
                    pattern_changed = 1;
                mem->idxbue[t][detected_count++] = i;
            }
        }
        if (detected_count != old_count) pattern_changed = 1;
        mem->nbue[t] = detected_count;
        }
        
        // STATE BOUNDS
        if (t > 0) {
        old_count = mem->nbxe[t];
        detected_count = 0;
        for (u32 i=0; i<(u32) dims->nbx[t]; ++i)
        {
            if (acados_daocp_is_soft_constraint(qp_in, t, dims->nbu[t]+i)) continue;
            u32 explicit_equality = acados_daocp_contains_index(
                qp_in->idxe[t], dims->nbue[t], dims->nbxe[t], dims->nbu[t]+(int)i);
            f64 lower = BLASFEO_DVECEL(qp_in->d+t, dims->nbu[t]+i);
            f64 upper = -BLASFEO_DVECEL(
                qp_in->d+t, dims->ng[t]+dims->nb[t]+dims->nbu[t]+i);
            if (explicit_equality || (DAOCP_ABS(upper-lower) <= opts->equality_detection_tolerance))
            {
                if (detected_count >= old_count
                    || mem->idxbxe[t][detected_count] != i)
                    pattern_changed = 1;
                mem->idxbxe[t][detected_count++] = i;
            }
        }
        if (detected_count != old_count) pattern_changed = 1;
        mem->nbxe[t] = detected_count;
        }

        // GENERAL CONSTRAINTS
        old_count = mem->nge[t];
        detected_count = 0;
        for (u32 i=0; i<(u32) dims->ng[t]; ++i)
        {
            if (acados_daocp_is_soft_constraint(qp_in, t, dims->nb[t]+i)) continue;
            u32 explicit_equality = acados_daocp_contains_index(
                qp_in->idxe[t], dims->nbue[t]+dims->nbxe[t], dims->nge[t], dims->nb[t]+(int)i);
            f64 lower = BLASFEO_DVECEL(qp_in->d+t, dims->nb[t]+i);
            f64 upper = -BLASFEO_DVECEL(
                qp_in->d+t, dims->ng[t]+2*dims->nb[t]+i);
            if (explicit_equality || (DAOCP_ABS(upper-lower) <= opts->equality_detection_tolerance))
            {
                if (detected_count >= old_count
                    || mem->idxge[t][detected_count] != i)
                    pattern_changed = 1;
                mem->idxge[t][detected_count++] = i;
            }
        }
        if (detected_count != old_count) pattern_changed = 1;
        mem->nge[t] = detected_count;
    }
    return pattern_changed;
}

static int acados_daocp_detect_changes_in_softening_pattern(
    const ocp_qp_daocp_memory* mem, const ocp_qp_in* qp
)
{
    const ocp_qp_dims* dims = qp->dim;
    for (int t=0; t<=dims->N; ++t)
    {
        int ns = dims->ns[t], nb = dims->nb[t], ng = dims->ng[t];
        int nbu = dims->nbu[t];
        // Get number of softened constraints for the previous problem
        int ns_prev = 0;
        for (int i=0; i<mem->qp.dims.nbu[t]; ++i) ns_prev += (mem->qp.Zbu[t][i] > 0);
        for (int i=0; i<mem->qp.dims.nbx[t]; ++i) ns_prev += (mem->qp.Zbx[t][i] > 0);
        for (int i=0; i<mem->qp.dims.ng[t]; ++i) ns_prev += (mem->qp.Zg[t][i] > 0);
        
        // If number of softened constraints is different, return
        if (ns != ns_prev) return 1;

        // Check that the softened constraints in the new problem are the
        // same as in the previous problem.
        for (int k=0; k<ns; ++k)
        {
            // Retrieve index of constraint associated to this slack
            int c_idx = -1;
            for (int i=0; (i<nb+ng) && c_idx<0; ++i) if (qp->idxs_rev[t][i] == k) c_idx = i;

            // Determine type of constraint
            daocp_constraint_type ctype;
            if (c_idx < nbu) ctype = DAOCP_BOUND_U;
            else if (c_idx < nb) ctype = DAOCP_BOUND_X;
            else ctype = DAOCP_MIXED;
            
            // Obtain constraint index in DAOCP's representation
            int c_idx_daocp;
            if (ctype == DAOCP_BOUND_U) {
                c_idx_daocp = c_idx;
                for (u32 i=0; i<mem->nbue[t] && mem->idxbue[t][i] < c_idx; ++i) c_idx_daocp -= 1; 
            } else if (ctype == DAOCP_BOUND_X) {
                c_idx_daocp = c_idx - nbu;
                for (u32 i=0; i<mem->nbxe[t] && mem->idxbxe[t][i] < c_idx-nbu; ++i) c_idx_daocp -= 1;
            } else {
                c_idx_daocp = c_idx - nb;
                for (u32 i=0; i<mem->nge[t] && mem->idxge[t][i] < c_idx-nb; ++i) c_idx_daocp -= 1;
            }

            // Check if the constraint was softened previously.
            if (ctype == DAOCP_BOUND_U) {
                if (mem->qp.Zbu[t][c_idx_daocp] == 0.0) return 1;
            } else if (ctype == DAOCP_BOUND_X) {
                if (mem->qp.Zbx[t][c_idx_daocp] == 0.0) return 1;
            } else {
                if (mem->qp.Zg[t][c_idx_daocp] == 0.0) return 1;
            } 
        }
    }
    return 0;
}

static void acados_daocp_process_constraints(
    ocp_qp_in* qp_in, ocp_qp_dims* dim, daocp_qp* qp_native,
    daocp_workspace* wrk, ocp_qp_daocp_memory* mem) {
    // Dimensions
    u32 N = qp_native->dims.N = (u32) qp_in->dim->N;
    qp_native->dims.nx = (u32*) dim->nx;
    qp_native->dims.nu = (u32*) dim->nu;
    for (u32 t=0; t<=N; ++t) {
        qp_native->dims.nbu[t] = dim->nbu[t] - mem->nbue[t];
        qp_native->dims.nbx[t] = t > 0 ? dim->nbx[t] - mem->nbxe[t] : 0;
        qp_native->dims.ng[t] = dim->ng[t] - mem->nge[t];
        qp_native->dims.ne[t] = mem->nbue[t] + mem->nbxe[t] + mem->nge[t];
    }

    // Extract x0
    for (u32 i=0; i<dim->nbxe[0]; ++i) {
        u32 bound_index = qp_in->idxe[0][dim->nbue[0]+i];
        u32 state_index = qp_in->idxb[0][bound_index] - dim->nu[0];
        qp_native->x0[state_index] = BLASFEO_DVECEL(qp_in->d, bound_index);
    }

    // u bounds
    for (u32 t=0; t<N; ++t) {
        f64* lbu = qp_native->lbu[t];
        f64* ubu = qp_native->ubu[t];
        u32* idxb = qp_native->idxbu[t];
        for (u32 i=0; i<dim->nbu[t]; ++i) {
            if (acados_daocp_is_equality(mem, t, i, DAOCP_BOUND_U)) continue;
            acados_daocp_populate_slack_penalties(qp_in, t, i,
                qp_native->Zbu[t] + (lbu-qp_native->lbu[t]),
                qp_native->zbu[t] + (lbu-qp_native->lbu[t]));
            *(lbu++) = BLASFEO_DVECEL(qp_in->d+t, i);
            *(ubu++) = -BLASFEO_DVECEL(qp_in->d+t, dim->ng[t]+dim->nb[t]+i);
            *(idxb++) = qp_in->idxb[t][i];
        }
    }
    // x bounds
    for (u32 t=1; t<=N; ++t) {
        f64* lbx = qp_native->lbx[t];
        f64* ubx = qp_native->ubx[t];
        u32* idxb = qp_native->idxbx[t];
        for (u32 i=0; i<dim->nbx[t]; ++i) {
            int bound_index = dim->nbu[t]+i;
            if (acados_daocp_is_equality(mem, t, i, DAOCP_BOUND_X)) continue;
            acados_daocp_populate_slack_penalties(qp_in, t, bound_index,
                qp_native->Zbx[t] + (lbx-qp_native->lbx[t]),
                qp_native->zbx[t] + (lbx-qp_native->lbx[t]));
            *(lbx++) = BLASFEO_DVECEL(qp_in->d+t, bound_index);
            *(ubx++) = -BLASFEO_DVECEL(qp_in->d+t, dim->ng[t]+dim->nb[t]+bound_index);
            *(idxb++) = qp_in->idxb[t][bound_index] - dim->nu[t];
        }
    }
    // general inequalities
    for (u32 t=0; t<=N; ++t) {
        f64* lb = qp_native->cl[t];
        f64* ub = qp_native->cu[t];
        f64* Cu = t < N ? qp_native->Cu[t] : 0;
        f64* Cx = t > 0 ? qp_native->Cx[t] : 0;
        daocp_constraint_type* contypes = wrk->contypes[t];
        for (u32 i=0; i<dim->ng[t]; ++i) {
            int constraint_index = dim->nb[t]+i;
            if (acados_daocp_is_equality(mem, t, i, DAOCP_MIXED))
                continue;
            acados_daocp_populate_slack_penalties(qp_in, t, constraint_index,
                qp_native->Zg[t] + (lb-qp_native->cl[t]),
                qp_native->zg[t] + (lb-qp_native->cl[t]));
            *lb = BLASFEO_DVECEL(qp_in->d+t, dim->nb[t]+i);
            *ub = -BLASFEO_DVECEL(qp_in->d+t, dim->ng[t]+2*dim->nb[t]+i);

            // Subtract Cx[0] x[0] contribution
            if (t==0) {
                for (u32 j=0; j<dim->nx[0]; ++j) {
                    *lb -= BLASFEO_DMATEL(qp_in->DCt+t, dim->nu[0]+j, i) * qp_native->x0[j];
                    *ub -= BLASFEO_DMATEL(qp_in->DCt+t, dim->nu[0]+j, i) * qp_native->x0[j];
                }
            }
            ++lb; ++ub;

            // Parse u-row and assess whether the constraint is x-only
            if (t != N) {
                for (u32 j=0; j<dim->nu[t]; ++j) Cu[j] = BLASFEO_DMATEL(qp_in->DCt+t, j, i);
                if (t == 0) *contypes = DAOCP_ONLY_U;
                else {
                    *contypes = DAOCP_MIXED;
                    u32 nonzero = 0;
                    for (u32 j=0; j<dim->nu[t]; ++j) 
                        if (DAOCP_ABS(Cu[j]) > DAOCP_ZERO_TOL) {
                            nonzero = 1;
                            break;
                        }
                    if (!nonzero) *contypes = DAOCP_ONLY_X;
                }
            }
            // Parse x-row
            if (t != 0) {
                for (u32 j=0; j<dim->nx[t]; ++j) Cx[j] = BLASFEO_DMATEL(qp_in->DCt+t, dim->nu[t]+j, i);
                if (t == N) *contypes = DAOCP_ONLY_X;
                else {
                    u32 nonzero = 0;
                    for (u32 j=0; j<dim->nx[t]; ++j) 
                        if (DAOCP_ABS(Cx[j]) > DAOCP_ZERO_TOL) {
                            nonzero = 1;
                            break;
                        }
                    if (!nonzero) *contypes = DAOCP_ONLY_U; 
                }
            }
            if (t < N) Cu += dim->nu[t];
            if (t > 0) Cx += dim->nx[t];
            ++contypes;
        }
    }
    
    // Equality constraints.
    for (u32 t=0; t<=N; ++t) {
        f64* d = qp_native->d[t];
        f64* Du = t < N ? qp_native->Du[t] : 0;
        f64* Dx = qp_native->Dx[t];

        u32 offset = 0;
        // Input bound equalities
        if (t < N) {
            for (u32 i=0; i<mem->nbue[t]; ++i, ++offset) {
                u32 bound_index = mem->idxbue[t][i];
                u32 uidx = qp_in->idxb[t][bound_index];
                d[offset] = BLASFEO_DVECEL(qp_in->d+t, bound_index);
                memset(Du + offset*dim->nu[t], 0, dim->nu[t]*sizeof(f64));
                memset(Dx + offset*dim->nx[t], 0, dim->nx[t]*sizeof(f64));
                Du[offset*dim->nu[t] + uidx] = 1.0;
            }
        }
        // State bound equalities
        if (t!=0) {
            for (u32 i=0; i<mem->nbxe[t]; ++i) {
                u32 bound_index = mem->idxbxe[t][i];
                u32 xidx = qp_in->idxb[t][dim->nbu[t]+bound_index] - dim->nu[t];
                d[offset] = BLASFEO_DVECEL(qp_in->d+t, dim->nbu[t]+bound_index);
                if (dim->nu[t] > 0)
                    memset(Du + offset*dim->nu[t], 0, dim->nu[t]*sizeof(f64));
                memset(Dx + offset*dim->nx[t], 0, dim->nx[t]*sizeof(f64));
                Dx[offset*dim->nx[t] + xidx] = 1.0;
                ++offset;
            }
        }
        // General equalities
        for (u32 i=0; i<mem->nge[t]; ++i) {
            u32 constraint_index = mem->idxge[t][i];
            d[offset] = BLASFEO_DVECEL(qp_in->d+t, dim->nb[t]+constraint_index);
            for (u32 j=0; j<dim->nu[t]; ++j)
                Du[offset*dim->nu[t]+j] = BLASFEO_DMATEL(qp_in->DCt+t, j, constraint_index);
            for (u32 j=0; j<dim->nx[t]; ++j)
                Dx[offset*dim->nx[t]+j] = BLASFEO_DMATEL(qp_in->DCt+t, dim->nu[t]+j, constraint_index);
            ++offset;
        }
        assert(offset == qp_native->dims.ne[t]);
    }
}

static void acados_daocp_init_workspace(daocp_workspace* wrk) {
    u32 acc = 0;
    for (u32 i=0; i<wrk->dims->N; ++i) {
        wrk->cnu[i] = acc;
        acc += wrk->dims->nu[i];
    }
    wrk->nu_tot = acc;
    acc = 0; for (u32 i=0; i<=wrk->dims->N; ++i) acc += wrk->dims->nbx[i]+wrk->dims->nbu[i]+wrk->dims->ng[i];
    wrk->W_stride = acc+1;
    daocp_reset_working_set(wrk);
}

static void acados_daocp_retrieve_slacks(ocp_qp_daocp_memory* mem, ocp_qp_dims* dim, ocp_qp_in* qp_in, ocp_qp_out* qp_out) {
    for (int t=0; t<=dim->N; ++t) {
        blasfeo_dvecse(2*dim->ns[t], 0.0, qp_out->ux+t, dim->nu[t]+dim->nx[t]);
        for (int i=0; i<dim->nb[t]+dim->ng[t]; ++i) {
            int slack = qp_in->idxs_rev[t][i];
            if (slack < 0) continue;
            double value = 0.0;
            if (i < dim->nb[t]) value = BLASFEO_DVECEL(qp_out->ux+t, qp_in->idxb[t][i]);
            else
                for (int j=0; j<dim->nu[t]+dim->nx[t]; ++j)
                    value += BLASFEO_DMATEL(qp_in->DCt+t, j, i-dim->nb[t]) * BLASFEO_DVECEL(qp_out->ux+t, j);
            int offset = dim->nu[t]+dim->nx[t]+slack;
            BLASFEO_DVECEL(qp_out->ux+t, offset) = fmax(0.0, BLASFEO_DVECEL(qp_in->d+t, i)-value);
            BLASFEO_DVECEL(qp_out->ux+t, offset+dim->ns[t]) = fmax(0.0, value+BLASFEO_DVECEL(qp_in->d+t, dim->nb[t]+dim->ng[t]+i));
        }
    } 
}

static void acados_daocp_retrieve_dual_sol(ocp_qp_daocp_memory* mem, ocp_qp_dims* dim, ocp_qp_in* qp_in, ocp_qp_out* qp_out, daocp_status return_status) {
    daocp_workspace* wrk = mem->workspace;

    // Dual solution is invalid if return status is INFEASIBLE. Zeroed out.
    if (return_status == DAOCP_INFEASIBLE) {
        for (int t=0; t<=dim->N; ++t) blasfeo_dvecse(2*(dim->nb[t]+dim->ng[t]+dim->ns[t]), 0.0, qp_out->lam+t, 0);
        return;
    }

    // Populate acados structs with equalities multiplers
    for (int t=0; t<=dim->N; ++t) {
        blasfeo_dvecse(2*(dim->nb[t]+dim->ng[t]+dim->ns[t]), 0.0, qp_out->lam+t, 0);
        
        for (int i=0; i<mem->nbue[t]; ++i) {
            BLASFEO_DVECEL(qp_out->lam + t, mem->idxbue[t][i]) = fmax(0.0, -mem->sol.eta[t][i]);
            BLASFEO_DVECEL(qp_out->lam + t, dim->nb[t]+dim->ng[t]+mem->idxbue[t][i]) = fmax(0.0, mem->sol.eta[t][i]);
        }
        for (int i=0; i<mem->nbxe[t]; ++i) {
            BLASFEO_DVECEL(qp_out->lam + t, dim->nbu[t] + mem->idxbxe[t][i]) = fmax(0.0, -mem->sol.eta[t][mem->nbue[t]+i]);
            BLASFEO_DVECEL(qp_out->lam + t, dim->nb[t]+dim->ng[t]+dim->nbu[t]+mem->idxbxe[t][i]) = fmax(0.0, mem->sol.eta[t][mem->nbue[t]+i]);
        }
        for (int i=0; i<mem->nge[t]; ++i) {
            BLASFEO_DVECEL(qp_out->lam + t, dim->nb[t]+mem->idxge[t][i]) = fmax(0.0, -mem->sol.eta[t][mem->nbue[t]+mem->nbxe[t]+i]);
            BLASFEO_DVECEL(qp_out->lam + t, dim->ng[t]+2*dim->nb[t]+mem->idxge[t][i]) = fmax(0.0, mem->sol.eta[t][mem->nbue[t]+mem->nbxe[t]+i]);
        }
    }

    // Populate acados structs with inequalities multipliers
    for (int i=0; i<wrk->as.n_active; ++i) {
        daocp_constraint* c = wrk->as.xi2con + i;
        u32 t = c->t;
        u32 idx = c->idx;
        u32 offset = c->is_upper ? dim->nb[t] + dim->ng[t] : 0;
        f64 scale = c->is_upper ? 1.0 : -1.0;
        if (c->type == DAOCP_BOUND_U) {
            for (u32 j=0; j<mem->nbue[t] && mem->idxbue[t][j] <= idx; ++j) ++idx;
            BLASFEO_DVECEL(qp_out->lam + t, offset+idx) = scale * wrk->xi[i];
        } else if (c->type == DAOCP_BOUND_X) {
            for (u32 j=0; j<mem->nbxe[t] && mem->idxbxe[t][j] <= idx; ++j) ++idx;
            BLASFEO_DVECEL(qp_out->lam + t, offset+dim->nbu[t]+idx) = scale * wrk->xi[i];
        } else {
            for (u32 j=0; j<mem->nge[t] && mem->idxge[t][j] <= idx; ++j) ++idx; 
            BLASFEO_DVECEL(qp_out->lam + t, offset+dim->nb[t]+idx) = scale * wrk->xi[i];
        }
    }

    // Populate slack multipliers
    for (int t=0; t<=dim->N; ++t) {
        for (int i=0; i<dim->nb[t]+dim->ng[t]; ++i) {
            int slack_idx = qp_in->idxs_rev[t][i];
            if (slack_idx < 0) continue;

            f64 sl = BLASFEO_DVECEL(qp_out->ux + t, dim->nx[t]+dim->nu[t]+slack_idx);
            f64 su = BLASFEO_DVECEL(qp_out->ux + t, dim->nx[t]+dim->nu[t]+dim->ns[t]+slack_idx);
            BLASFEO_DVECEL(
                qp_out->lam + t, 2*(dim->nb[t]+dim->ng[t]) + slack_idx
            ) = BLASFEO_DVECEL(qp_in->Z+t, slack_idx)*sl
                + BLASFEO_DVECEL(qp_in->rqz+t, dim->nx[t]+dim->nu[t]+slack_idx) 
                - BLASFEO_DVECEL(qp_out->lam+t, i);
            BLASFEO_DVECEL(
                qp_out->lam + t, 2*(dim->nb[t]+dim->ng[t]) + dim->ns[t] + slack_idx
            ) = BLASFEO_DVECEL(qp_in->Z+t, slack_idx)*su
                + BLASFEO_DVECEL(qp_in->rqz+t, dim->nx[t]+dim->nu[t]+slack_idx) 
                - BLASFEO_DVECEL(qp_out->lam+t, dim->nb[t]+dim->ng[t]+i);
        }
    }
}

static void acados_daocp_retrieve_equality_multipliers(ocp_qp_daocp_memory* mem, ocp_qp_dims* dim, ocp_qp_in* qp_in, ocp_qp_out* qp_out, daocp_status return_status) {
    // Equality multipliers must be computed only after the other dual multipliers
    // have already been translated into acados' format.
    daocp_workspace* wrk = mem->workspace;
    struct blasfeo_dvec tmp;
    int N = dim->N;

    // Initialize pi_N = P_N x_N + p_N + D_N' (\lam_{g,N}^u - \lam_{g,N}^l)
    //                   + \lam_{bx,N}^u - \lam_{bx,N}^l 
    blasfeo_dsymv_l(dim->nx[N], 1.0, qp_in->RSQrq+N, 0, 0, qp_out->ux+N, 0, 1.0, 
                    qp_in->rqz+N, 0, qp_out->pi + N-1, 0);
    
    // Don't add dual contribution if problem is infeasible.
    if (return_status != DAOCP_INFEASIBLE) {
        // General constraints contribution.
        tmp.pa = mem->sol.lam[N];
        blasfeo_daxpy(dim->ng[N], -1.0, qp_out->lam+N, dim->nb[N], qp_out->lam+N,
                    2*dim->nb[N]+dim->ng[N], &tmp, 0);
        blasfeo_dgemv_n(dim->nx[N], dim->ng[N], 1.0, qp_in->DCt+N, 0, 0, &tmp, 0, 1.0, qp_out->pi+N-1,
                    0, qp_out->pi+N-1, 0);
        // State bounds contribution
        for (int i=0; i<dim->nbx[N]; ++i) {
            int x_idx = qp_in->idxb[N][i];
            BLASFEO_DVECEL(qp_out->pi+N-1, x_idx) -= 
                BLASFEO_DVECEL(qp_out->lam+N, i)
                - BLASFEO_DVECEL(qp_out->lam+N, dim->nb[N]+dim->ng[N]+i);
        }
    }
    for (int t=N-1; t>=1; --t) {
        // pi_t = A_t' pi_t + Q_t x_t + S_t' u_t + q_t 
        //        + D_t' (\lam_{g,t}^u - \lam_{g,t}^l) + \lam_{bx,t}^u - \lam_{bx,t}^l
        blasfeo_dgemv_n(dim->nx[t], dim->nx[t+1], 1.0, qp_in->BAbt+t, dim->nu[t], 0, qp_out->pi+t,
                        0, 1.0, qp_in->rqz+t, dim->nu[t], qp_out->pi+t-1, 0);
        blasfeo_dsymv_l(dim->nx[t], 1.0, qp_in->RSQrq+t, dim->nu[t], dim->nu[t], qp_out->ux+t,
                        dim->nu[t], 1.0, qp_out->pi+t-1, 0, qp_out->pi+t-1, 0);
        blasfeo_dgemv_n(dim->nx[t], dim->nu[t], 1.0, qp_in->RSQrq+t, dim->nu[t], 0, qp_out->ux+t,
                        0, 1.0, qp_out->pi+t-1, 0, qp_out->pi+t-1, 0);
        
        // Do not add dual contributions if the problem was infeasible
        if (return_status == DAOCP_INFEASIBLE) continue;

        // Compute pi_t += D_t' (\lam_{g,t}^u - \lam_{g,t}^l).
        // Use mem->sol.lam[t] as temporary buffer to compute the difference.
        // This buffer is guaranteed to have at least ng[t] room.
        // We overwrite the dual export in daocp's format, but that's not an issue.
        tmp.pa = mem->sol.lam[t];
        blasfeo_daxpy(dim->ng[t], -1.0, qp_out->lam+t, dim->nb[t], qp_out->lam+t,
                      2*dim->nb[t]+dim->ng[t], &tmp, 0);
        blasfeo_dgemv_n(dim->nx[t], dim->ng[t], 1.0, qp_in->DCt+t, dim->nu[t], 0, &tmp, 0, 1.0, qp_out->pi+t-1,
                        0, qp_out->pi+t-1, 0);

        // State bounds contribution
        for (int i=0; i<dim->nbx[t]; ++i) {
            int x_idx = qp_in->idxb[t][dim->nbu[t]+i] - dim->nu[t];
            BLASFEO_DVECEL(qp_out->pi+t-1, x_idx) -= 
                BLASFEO_DVECEL(qp_out->lam+t, dim->nbu[t]+i)
                - BLASFEO_DVECEL(qp_out->lam+t, dim->nb[t]+dim->ng[t]+dim->nbu[t]+i);
        }
    }

    // Initial condition multiplier
    struct blasfeo_dvec* mu = &wrk->costate0;
    // -\mu = A_t' pi_1 + Q_0 x_0 + S_0' u_0 + q_0 + D_0' (\lam_{g,0}^u - \lam_{g,0}^l)
    blasfeo_dgemv_n(dim->nx[0], dim->nx[1], -1.0, qp_in->BAbt, dim->nu[0], 0, qp_out->pi,
                    0, -1.0, qp_in->rqz, dim->nu[0], mu, 0);
    blasfeo_dsymv_l(dim->nx[0], -1.0, qp_in->RSQrq, dim->nu[0], dim->nu[0], qp_out->ux,
                    dim->nu[0], 1.0, mu, 0, mu, 0);
    blasfeo_dgemv_n(dim->nx[0], dim->nu[0], -1.0, qp_in->RSQrq, dim->nu[0], 0, qp_out->ux,
                    0, 1.0, mu, 0, mu, 0);

    // Do not add dual contributions if the problem was infeasible
    if (return_status != DAOCP_INFEASIBLE) {
        tmp.pa = mem->sol.lam[0];
        blasfeo_daxpy(dim->ng[0], -1.0, qp_out->lam, dim->nb[0], qp_out->lam,
                    2*dim->nb[0]+dim->ng[0], &tmp, 0);
        blasfeo_dgemv_n(dim->nx[0], dim->ng[0], -1.0, qp_in->DCt, dim->nu[0], 0, &tmp, 0, 1.0, mu,
                            0, mu, 0);
    }
    for (int i = 0; i < dim->nbxe[0]; ++i) {
        int b = qp_in->idxe[0][dim->nbue[0] + i];
        int j = qp_in->idxb[0][b] - dim->nu[0];
        BLASFEO_DVECEL(qp_out->lam, b) = fmax(-BLASFEO_DVECEL(mu, j), 0.0);
        BLASFEO_DVECEL(qp_out->lam, dim->nb[0]+dim->ng[0]+b) = fmax(BLASFEO_DVECEL(mu, j), 0.0);
    }
}

int ocp_qp_daocp(void *config_, void *qp_in_, void *qp_out_, void *opts_, void *mem_, void *work_)
{
    ocp_qp_in *qp_in = qp_in_;
    ocp_qp_out *qp_out = qp_out_;
    ocp_qp_dims *dim = qp_in->dim;
    ocp_qp_daocp_memory* mem = mem_;
    ocp_qp_daocp_opts* opts = opts_;

    qp_info* info = (qp_info* ) qp_out->misc;
    acados_timer tot_timer, qp_timer, interface_timer, solver_call_timer;
    acados_tic(&tot_timer);
    acados_tic(&interface_timer);

    // QP validation
    if (qp_in->dim->N < 1) {
        printf("\nDAOCP requires a horizon of at least one shooting interval.\n");
        return ACADOS_QP_FAILURE;
    }
    if (qp_in->dim->nx[0] != qp_in->dim->nbxe[0]) {
        printf("\nDAOCP can only handle problems with fixed initial state.\n");
        return ACADOS_QP_FAILURE;
    }
    if (dim->nu[dim->N] != 0) {
        printf("DAOCP requires zero terminal controls.\n");
        return ACADOS_QP_FAILURE;
    }
    if (acados_daocp_validate_slacks(qp_in) != ACADOS_SUCCESS)
        return ACADOS_QP_FAILURE;
    
    // Conversion of data structures
    daocp_qp* qp_native = &mem->qp;
    // (Shadow) copy dynamics and cost
    qp_native->BAwt = qp_in->BAbt;
    qp_native->RSQrq = qp_in->RSQrq;
    // Write affine term of dynamics in qp_native->BAwt
    for (int t=0; t<dim->N; ++t)
        blasfeo_drowin(dim->nx[t+1], 1.0, qp_in->b+t, 0, qp_in->BAbt+t, dim->nu[t]+dim->nx[t], 0);
    // Write linear tems of cost in qp_native->RSQrq
    for (int t=0; t<=dim->N; ++t)
        blasfeo_drowin(dim->nu[t]+dim->nx[t], 1.0, qp_in->rqz+t, 0, qp_in->RSQrq+t, dim->nu[t]+dim->nx[t], 0);
    daocp_workspace* wrk = (daocp_workspace*) mem->workspace;
    wrk->dims = &qp_native->dims;
    // Detect changes in the pattern of softened constraints
    u32 slack_pattern_changed = 0;
    if (!opts->first_run) 
        slack_pattern_changed = acados_daocp_detect_changes_in_softening_pattern(mem, qp_in);
    // Detect equalities from problem data
    u32 equality_pattern_changed = acados_daocp_detect_equalities(qp_in, opts, mem);
    acados_daocp_process_constraints(qp_in, dim, qp_native, wrk, mem);
    /* Slack mappings and explicit equality classification can change between solves.
     * Start from an empty active set for soft QPs until these changes are tracked. */
    u32 reset_working_set = !opts->warm_start || equality_pattern_changed || slack_pattern_changed;
    if (opts->first_run || reset_working_set) acados_daocp_init_workspace(wrk);
    // Set solution pointers
    mem->sol.ux = qp_out->ux;

    // Interface work done
    info->interface_time = acados_toc(&interface_timer);

    acados_tic(&qp_timer);
    // The following calls before daocp_solve are effectively
    // part of the solver.
    acados_tic(&solver_call_timer); 
    daocp_solve_riccati(wrk, qp_native);
    daocp_solve_lqr(wrk, qp_native);
    daocp_compute_dual_linear_term(wrk, qp_native);
    // Invalidate chacked L_d y = -d solution
    wrk->as.n_valid_intermediate = 0;
    
    // We check whether we can use the active set information from 
    // the previous solve.
    if (!opts->first_run && !reset_working_set) {
        // Recompute cholesky of dual hessian, detecting singularity
        reset_working_set = daocp_compute_chol_from_scratch(wrk, qp_native);
        if (reset_working_set) daocp_reset_working_set(wrk);
    }
    daocp_solve(&opts->daocp_opts, qp_native, wrk, &mem->sol);
    opts->first_run = 0;
    mem->time_qp_solver_call = acados_toc(&solver_call_timer);

    acados_daocp_retrieve_slacks(mem, dim, qp_in, qp_out);
    acados_daocp_retrieve_dual_sol(mem, dim, qp_in, qp_out, wrk->status);
    acados_daocp_retrieve_equality_multipliers(mem, dim, qp_in, qp_out, wrk->status);
    ocp_qp_compute_t(qp_in, qp_out);

    if (opts->print_level > 0)
        printf("DA-OCP status: %d, iterations: %u, active constraints: %u\n",
               wrk->status, wrk->iters, wrk->as.n_active);

    info->solve_QP_time = acados_toc(&qp_timer);
    info->total_time = acados_toc(&tot_timer);
    info->num_iter = wrk->iters;

    // status
    int acados_status = ACADOS_QP_FAILURE; // generic QP failure
    daocp_status status = wrk->status;
    if (status==DAOCP_SOLVED)
        acados_status = ACADOS_SUCCESS;
    else if (status==DAOCP_MAX_ITER)
        acados_status = ACADOS_MAXITER;
    else if (status==DAOCP_INFEASIBLE)
        acados_status = ACADOS_INFEASIBLE;

    return acados_status;
}



void ocp_qp_daocp_eval_adj_sens(void *config_, void *param_qp_in_, void *seed, void *sens_qp_out_, void *opts_, void *mem_, void *work_)
{
    printf("\nerror: ocp_qp_daocp_eval_adj_sens: not implemented yet\n");
    exit(1);
}

void ocp_qp_daocp_eval_forw_sens(void *config_, void *param_qp_in_, void *seed, void *sens_qp_out_, void *opts_, void *mem_, void *work_)
{
    printf("\nerror: ocp_qp_daocp_eval_forw_sens: not implemented yet\n");
    exit(1);
}

void ocp_qp_daocp_solver_get(void *config_, void *qp_in_, void *qp_out_, void *opts_, void *mem_, const char *field, int stage, void* value, int size1, int size2)
{
    printf("\nerror: ocp_qp_daocp_solver_get: not implemented yet\n");
    exit(1);
}


void ocp_qp_daocp_terminate(void *config_, void *mem_, void *work_)
{
}




void ocp_qp_daocp_config_initialize_default(void *config_)
{
    qp_solver_config *config = config_;

    config->opts_calculate_size = &ocp_qp_daocp_opts_calculate_size;
    config->opts_assign = &ocp_qp_daocp_opts_assign;
    config->opts_initialize_default = &ocp_qp_daocp_opts_initialize_default;
    config->opts_update = &ocp_qp_daocp_opts_update;
    config->opts_set = &ocp_qp_daocp_opts_set;
    config->opts_get = &ocp_qp_daocp_opts_get;
    config->memory_calculate_size = &ocp_qp_daocp_memory_calculate_size;
    config->memory_assign = &ocp_qp_daocp_memory_assign;
    config->memory_get = &ocp_qp_daocp_memory_get;
    config->workspace_calculate_size = &ocp_qp_daocp_workspace_calculate_size;
    config->evaluate = &ocp_qp_daocp;
    config->terminate = &ocp_qp_daocp_terminate;
    config->eval_forw_sens = &ocp_qp_daocp_eval_forw_sens;
    config->eval_adj_sens = &ocp_qp_daocp_eval_adj_sens;
    config->memory_reset = &ocp_qp_daocp_memory_reset;
    config->solver_get = &ocp_qp_daocp_solver_get;

    return;
}
