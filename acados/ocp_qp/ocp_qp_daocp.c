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
    else
    {
        printf("\nerror: ocp_qp_daocp_opts_get: field %s not available\n", field);
        exit(1);
    }
}




/************************************************
 * memory
 ************************************************/

static inline u32 acados_daocp_num_equalities_at_stage(const ocp_qp_dims* dims, u32 stage)
{
    return dims->nge[stage] + dims->nbue[stage] + (stage > 0 ? dims->nbxe[stage] : 0);
}

static inline acados_size_t acados_daocp_num_elimination_alphas(
    acados_size_t nrows, acados_size_t rank)
{
    return rank > 0 ? rank * (2*nrows - rank - 1) / 2 : 0;
}

static void acados_daocp_equality_elimination_capacity(
    const ocp_qp_dims* dims, acados_size_t* alpha_capacity, acados_size_t* pivot_capacity)
{
    acados_size_t remaining_equalities = 0;
    *alpha_capacity = 0;
    *pivot_capacity = 0;

    for (u32 stage = 0; stage <= (u32) dims->N; stage++)
        remaining_equalities += acados_daocp_num_equalities_at_stage(dims, stage);

    for (u32 stage = 0; stage < (u32) dims->N; stage++)
    {
        acados_size_t rank = DAOCP_MIN((acados_size_t) dims->nu[stage], remaining_equalities);
        *pivot_capacity += rank;
        *alpha_capacity += acados_daocp_num_elimination_alphas(remaining_equalities, rank);
        remaining_equalities -= acados_daocp_num_equalities_at_stage(dims, stage);
    }
}

acados_size_t ocp_qp_daocp_memory_calculate_size(void *config_, void *dims_, void *opts_)
{
    ocp_qp_dims *dims = dims_;
    acados_size_t size = sizeof(ocp_qp_daocp_memory);
    size += sizeof(daocp_workspace);

    int N = dims->N;
    int* nx = dims->nx;
    int* nu = dims->nu;
    int* nb = dims->nb;
    int* ng = dims->ng; 
    int* nge = dims->nge;
    int* nbue = dims->nbue;
    int* nbxe = dims->nbxe;

    // daocp_qp data
    size += (4*N + 3)*sizeof(u32); // dims.ng, dims.ne, dims.nbu, dims.nbx
    size += nx[0]*sizeof(f64); // x0
    u32 neq = 0;
    u32 nin = 0;
    u32 nb_tot = 0;
    for (u32 t=0; t<=N; ++t) {
        u32 ne_t = nge[t]+nbue[t]+nbxe[t]*(t>0 ? 1 : 0);
        size += 2*(nb[t]-nbue[t]-nbxe[t])*sizeof(f64); // lbu/x[t], ubu/x[t]
        size += 2*(nb[t]-nbue[t]-nbxe[t])*sizeof(f64); // lbu/x_wrk[t], ubu/x_wrk[t]
        size += (nb[t]-nbue[t]-nbxe[t])*sizeof(u32); // idxbu[t]
        size += 2*(ng[t]-nge[t])*sizeof(f64); // lg[t], ug[t]
        size += 2*(ng[t]-nge[t])*sizeof(f64); // lg_wrk[t], ug_wrk[t]
        size += ne_t*sizeof(f64); // d[t]
        size += (ng[t]-nge[t])*(nu[t]+nx[t])*sizeof(f64); // Cu[t], Cx[t]
        size += (t<N ? ne_t : 0)*nu[t]*sizeof(f64); // Du[t]
        size += ne_t*nx[t]*sizeof(f64); // Dx[t]

        neq += ne_t;
        nin += nb[t]+ng[t]-(nge[t]+nbue[t]+nbxe[t]);
        nb_tot += nb[t]-nbue[t]-nbxe[t];
    }
    size -= (ng[0]-nge[0])*nx[0]*sizeof(f64); // Subtract off Cx[0]
    size += (2*N+1)*sizeof(f64*); // Cu[:], Cx[:]
    size += (2*N+1)*sizeof(f64*); // Du[:], Dx[:]
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
    size += (3*N+1)*sizeof(f64*); // x[:], u[:], eta[:]
    size += (N+1)*sizeof(daocp_constraint_type*); // contypes[:]
    size += (nin-nb_tot)*sizeof(daocp_constraint_type); // contypes
    size += 3*N*sizeof(u32); // cnu, rho, crho
    size += (N+1)*sizeof(u32*); // as.constraint_status[:]
    size += nin*sizeof(u32); // as.constraint_status
        
    u32 W_stride = DAOCP_MIN(tot_nu, nin)+1;
    size += W_stride*sizeof(daocp_constraint); // as.xi2con
    size += 3*W_stride*sizeof(f64); // xi, p, dual_linear
    size += W_stride*sizeof(u32); // xi_sign
    size += (W_stride+1)*W_stride*sizeof(f64); // Ld
    size += 2*W_stride*tot_nu*sizeof(f64); // Mu, Me

    size += neq*max_nx*sizeof(f64); // H
    size += 2*neq*sizeof(f64); // h, tmp1
    size += (neq+1)*(max_nx+max_nu+1)*sizeof(f64); // GEtmp
    size += max_nx*DAOCP_MAX(max_nx, max_nu)*sizeof(f64); // ABtmp
    size += DAOCP_MAX(blasfeo_memsize_dmat(max_nx, max_nx), blasfeo_memsize_dmat(max_nu, max_nx)); // tmp2
    size += 2*blasfeo_memsize_dvec(DAOCP_MAX(max_nx, max_nu)); // costate0, costate1

    // Row operations needed to map DA-OCP's eliminated equality multipliers
    // back to the original acados equality rows. These buffers are a tape:
    // they avoid allocations while recovering the dual solution.
    acados_size_t alpha_capacity, pivot_capacity;
    acados_daocp_equality_elimination_capacity(dims, &alpha_capacity, &pivot_capacity);
    size += alpha_capacity*sizeof(f64);
    size += pivot_capacity*sizeof(u32);

    // Worst-case padding for three double-alignment points and two cache-line
    // aligned BLASFEO backing-store regions.
    size += 3*8 + 2*64;

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
    memset(raw_memory, 0, ocp_qp_daocp_memory_calculate_size(config_, dims_, opts_));
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
    int* nge = dims->nge;
    int* nbue = dims->nbue;
    int* nbxe = dims->nbxe;

    // daocp_qp data
    mem->qp.dims.nbu = (u32*) (c_ptr += sizeof(daocp_workspace));
    mem->qp.dims.nbx = (u32*) (c_ptr += N*sizeof(u32));
    mem->qp.dims.ng = (u32*) (c_ptr += (N+1)*sizeof(u32));
    mem->qp.dims.ne = (u32*) (c_ptr += (N+1)*sizeof(u32));
    c_ptr += (N+1)*sizeof(u32);
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
    mem->qp.Dx = (f64**) (c_ptr += N*sizeof(f64*));
    mem->qp.idxbu = (u32**) (c_ptr += (N+1)*sizeof(f64*));
    mem->qp.idxbx = (u32**) (c_ptr += N*sizeof(u32*));
    c_ptr += (N+1)*sizeof(u32*);
    c_ptr = assign_ptr_vec(mem->qp.lbu, c_ptr, nbu, nbue, 0, sizeof(f64), N);
    c_ptr = assign_ptr_vec(wrk->lbu_wrk, c_ptr, nbu, nbue, 0, sizeof(f64), N);
    mem->qp.lbx[0] = 0; c_ptr = assign_ptr_vec(mem->qp.lbx+1, c_ptr, nbx+1, nbxe+1, 0, sizeof(f64), N);
    wrk->lbx_wrk[0] = 0; c_ptr = assign_ptr_vec(wrk->lbx_wrk+1, c_ptr, nbx+1, nbxe+1, 0, sizeof(f64), N);
    c_ptr = assign_ptr_vec(mem->qp.ubu, c_ptr, nbu, nbue, 0, sizeof(f64), N);
    c_ptr = assign_ptr_vec(wrk->ubu_wrk, c_ptr, nbu, nbue, 0, sizeof(f64), N);
    mem->qp.ubx[0] = 0; c_ptr = assign_ptr_vec(mem->qp.ubx+1, c_ptr, nbx+1, nbxe+1, 0, sizeof(f64), N);
    wrk->ubx_wrk[0] = 0; c_ptr = assign_ptr_vec(wrk->ubx_wrk+1, c_ptr, nbx+1, nbxe+1, 0, sizeof(f64), N);
    c_ptr = assign_ptr_vec(mem->qp.idxbu, c_ptr, nbu, nbue, 0, sizeof(u32), N);
    mem->qp.idxbx[0] = 0; c_ptr = assign_ptr_vec(mem->qp.idxbx+1, c_ptr, nbx+1, nbxe+1, 0, sizeof(u32), N);
    align_char_to(8, &c_ptr);
    c_ptr = assign_ptr_vec(mem->qp.cl, c_ptr, ng, nge, 0, sizeof(f64), N+1);
    c_ptr = assign_ptr_vec(wrk->lg_wrk, c_ptr, ng, nge, 0, sizeof(f64), N+1);
    c_ptr = assign_ptr_vec(mem->qp.cu, c_ptr, ng, nge, 0, sizeof(f64), N+1);
    c_ptr = assign_ptr_vec(wrk->ug_wrk, c_ptr, ng, nge, 0, sizeof(f64), N+1);
    mem->qp.d[0] = (f64*) c_ptr; c_ptr += (nbue[0]+nge[0])*sizeof(f64);
    for (u32 t=1; t<=N; ++t) {
        mem->qp.d[t] = (f64*) c_ptr;
        c_ptr += (nbue[t]+nbxe[t]+nge[t])*sizeof(f64);
    }
    c_ptr = assign_ptr_mat(mem->qp.Cu, c_ptr, ng, nge, 0, 0, nu, sizeof(f64), N);
    mem->qp.Cx[0]=0; c_ptr = assign_ptr_mat(mem->qp.Cx+1, c_ptr, ng+1, nge+1, 0, 0, nx+1, sizeof(f64), N);
    mem->qp.Du[0] = (f64*) c_ptr; c_ptr += (nge[0]+nbue[0])*nu[0]*sizeof(f64);
    c_ptr = assign_ptr_mat(mem->qp.Du+1, c_ptr, nge+1, 0, nbue+1, nbxe+1, nu+1, sizeof(f64), N-1);
    mem->qp.Dx[0]= (f64*) c_ptr; c_ptr += (nge[0]+nbue[0])*nx[0]*sizeof(f64);
    c_ptr = assign_ptr_mat(mem->qp.Dx+1, c_ptr, nge+1, 0, nbue+1, nbxe+1, nx+1, sizeof(f64), N);

    // daocp_workspace data
    int max_nx = 0; int max_nu = 0; int neq = 0; int nin = 0; int tot_nu = 0;
    for (u32 t=0; t<=N; ++t) if (nx[t] > max_nx) max_nx = nx[t];
    for (u32 t=0; t<N; ++t) {
        if (nu[t] > max_nu) max_nu = nu[t];
        tot_nu += nu[t];
    }
    for (u32 t=0; t<=N; ++t) neq += nge[t] + nbue[t] + nbxe[t];
    neq -= nbxe[0];
    for (u32 t=0; t<=N; ++t) nin += ng[t]+nb[t]-nge[t]-nbue[t]-nbxe[t];

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
    wrk->contypes = (daocp_constraint_type**) c_ptr; c_ptr += (N+1)*sizeof(daocp_constraint_type*);
    wrk->as.constraint_status = (u32**) c_ptr; c_ptr += (N+1)*sizeof(u32*);
    c_ptr = assign_ptr_vec(wrk->u, c_ptr, nu, 0, 0, sizeof(f64), N);
    wrk->x[0]=0; c_ptr = assign_ptr_vec(wrk->x+1, c_ptr, nx+1, 0, 0, sizeof(f64), N);
    if (N > 0) wrk->eta[0] = (f64*) c_ptr;
    c_ptr += tot_nu*sizeof(f64);
    for (u32 t=0; t<=N; ++t) {
        wrk->contypes[t] = (daocp_constraint_type*) c_ptr;
        c_ptr += (ng[t]-nge[t])*sizeof(daocp_constraint_type);
    }
    for (u32 t=0; t<=N; ++t) {
        wrk->as.constraint_status[t] = (u32*) c_ptr;
        c_ptr += (ng[t]+nb[t]-nbue[t]-nbxe[t]-nge[t])*sizeof(u32);
    }
    wrk->cnu = (u32*) c_ptr; c_ptr += N*sizeof(u32);
    wrk->rho = (u32*) c_ptr; c_ptr += N*sizeof(u32);
    wrk->crho = (u32*) c_ptr; c_ptr += N*sizeof(u32);

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
    
    u32 W_stride = DAOCP_MIN(tot_nu, nin)+1;
    wrk->as.xi2con = (daocp_constraint*) c_ptr; c_ptr += W_stride*sizeof(daocp_constraint);
    wrk->xi = (f64*) c_ptr; c_ptr+=W_stride*sizeof(f64);
    wrk->p = (f64*) c_ptr; c_ptr+=W_stride*sizeof(f64);
    wrk->dual_linear = (f64*) c_ptr; c_ptr+=W_stride*sizeof(f64);
    wrk->xi_sign = (u32*) c_ptr; c_ptr+=W_stride*sizeof(u32);
    align_char_to(8, &c_ptr);
    wrk->Ld = (f64*) c_ptr; c_ptr+=(W_stride+1)*W_stride*sizeof(f64);
    wrk->Mu = (f64*) c_ptr; c_ptr+=W_stride*tot_nu*sizeof(f64);
    wrk->Me = (f64*) c_ptr; c_ptr+=W_stride*tot_nu*sizeof(f64);
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

    acados_size_t alpha_capacity, pivot_capacity;
    acados_daocp_equality_elimination_capacity(dims, &alpha_capacity, &pivot_capacity);
    mem->equality_elimination_alphas = (f64*) c_ptr;
    c_ptr += alpha_capacity*sizeof(f64);
    mem->equality_elimination_pivots = (u32*) c_ptr;
    c_ptr += pivot_capacity*sizeof(u32);

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

static void acados_daocp_process_constraints(ocp_qp_in* qp_in, ocp_qp_dims* dim, daocp_qp* qp_native, daocp_workspace* wrk) {
    // Dimensions
    u32 N = qp_native->dims.N = (u32) qp_in->dim->N;
    qp_native->dims.nx = (u32*) dim->nx;
    qp_native->dims.nu = (u32*) dim->nu;
    for (u32 t=0; t<N; ++t)
        qp_native->dims.nbu[t] = dim->nbu[t] - dim->nbue[t];
    qp_native->dims.nbu[N] = 0;
    qp_native->dims.nbx[0] = 0;
    for (u32 t=1; t<=N; ++t)
        qp_native->dims.nbx[t] = dim->nbx[t] - qp_in->dim->nbxe[t];
    for (u32 t=0; t<=N; ++t) {
        qp_native->dims.ng[t] = dim->ng[t] - dim->nge[t];
        qp_native->dims.ne[t] = dim->nge[t] + dim->nbue[t] + dim->nbxe[t];
    }
    qp_native->dims.ne[0] -= dim->nbxe[0];

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
            if (acados_daocp_contains_index(qp_in->idxe[t], 0, dim->nbue[t], i))
                continue;
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
            if (acados_daocp_contains_index(qp_in->idxe[t], dim->nbue[t], dim->nbxe[t], bound_index))
                continue;
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
            if (acados_daocp_contains_index(qp_in->idxe[t], dim->nbue[t]+dim->nbxe[t], dim->nge[t], constraint_index))
                continue;
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

        // Input bound equalities
        for (u32 i=0; i<dim->nbue[t]; ++i) {
            int bound_index = qp_in->idxe[t][i];
            u32 uidx = qp_in->idxb[t][bound_index];
            d[i] = BLASFEO_DVECEL(qp_in->d+t, bound_index);
            if (dim->nu[t] > 0)
                memset(Du + i*dim->nu[t], 0, dim->nu[t]*sizeof(f64));
            memset(Dx + i*dim->nx[t], 0, dim->nx[t]*sizeof(f64));
            Du[i*dim->nu[t] + uidx] = 1.0;
        }
        // State bound equalities
        if (t!=0) {
            for (u32 i=0; i<dim->nbxe[t]; ++i) {
                u32 offset = dim->nbue[t] + i;
                int bound_index = qp_in->idxe[t][dim->nbue[t]+i];
                u32 xidx = qp_in->idxb[t][bound_index] - dim->nu[t];
                d[offset] = BLASFEO_DVECEL(qp_in->d+t, bound_index);
                if (dim->nu[t] > 0)
                    memset(Du + offset*dim->nu[t], 0, dim->nu[t]*sizeof(f64));
                memset(Dx + offset*dim->nx[t], 0, dim->nx[t]*sizeof(f64));
                Dx[offset*dim->nx[t] + xidx] = 1.0;
            }
        }
        // General equalities
        for (u32 i=0; i<dim->nge[t]; ++i) {
            u32 offset = dim->nbue[t] + (t > 0 ? dim->nbxe[t] : 0) + i;
            int constraint_index = qp_in->idxe[t][dim->nbue[t]+dim->nbxe[t]+i];
            int general_index = constraint_index - dim->nb[t];
            d[offset] = BLASFEO_DVECEL(qp_in->d+t, constraint_index);
            for (u32 j=0; j<dim->nu[t]; ++j)
                Du[offset*dim->nu[t]+j] = BLASFEO_DMATEL(qp_in->DCt+t, j, general_index);
            for (u32 j=0; j<dim->nx[t]; ++j)
                Dx[offset*dim->nx[t]+j] = BLASFEO_DMATEL(qp_in->DCt+t, dim->nu[t]+j, general_index);
        }
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
    wrk->W_stride = DAOCP_MIN(acc, wrk->nu_tot)+1;

    wrk->singular = 0;
    wrk->as.n_active = 0;
    wrk->as.max_t = 0;
    for (u32 t=0; t<=wrk->dims->N; ++t) {
        if (t<wrk->dims->N) for (u32 i=0; i<wrk->dims->nbu[t]; ++i) wrk->as.constraint_status[t][i] = 0;
        if (t>0) for (u32 i=0; i<wrk->dims->nbx[t]; ++i) wrk->as.constraint_status[t][wrk->dims->nbu[t] + i] = 0;
        for (u32 i=0; i<wrk->dims->ng[t]; ++i) wrk->as.constraint_status[t][wrk->dims->nbu[t]+wrk->dims->nbx[t] + i] = 0;
    }
}

static u32 acados_daocp_original_inequality_index(
    const ocp_qp_in* qp_in, u32 stage, const daocp_constraint* constraint)
{
    const ocp_qp_dims* dims = qp_in->dim;
    int first_index;
    int count;
    int equality_offset;
    int equality_count;

    switch (constraint->type)
    {
        case DAOCP_BOUND_U:
            first_index = 0;
            count = dims->nbu[stage];
            equality_offset = 0;
            equality_count = dims->nbue[stage];
            break;
        case DAOCP_BOUND_X:
            first_index = dims->nbu[stage];
            count = dims->nbx[stage];
            equality_offset = dims->nbue[stage];
            equality_count = dims->nbxe[stage];
            break;
        default:
            first_index = dims->nb[stage];
            count = dims->ng[stage];
            equality_offset = dims->nbue[stage] + dims->nbxe[stage];
            equality_count = dims->nge[stage];
            break;
    }

    u32 filtered_index = 0;
    for (int i = 0; i < count; i++)
    {
        int original_index = first_index + i;
        if (acados_daocp_contains_index(
                qp_in->idxe[stage], equality_offset, equality_count, original_index))
            continue;
        if (filtered_index == constraint->idx)
            return (u32) original_index;
        filtered_index++;
    }

    assert(0 && "DAOCP inequality index does not map to an acados constraint");
    return 0;
}

static u32 acados_daocp_original_equality_index(
    const ocp_qp_in* qp_in, u32 stage, u32 equality_row)
{
    const ocp_qp_dims* dims = qp_in->dim;

    if (equality_row < (u32) dims->nbue[stage])
        return (u32) qp_in->idxe[stage][equality_row];

    equality_row -= dims->nbue[stage];
    if (stage > 0 && equality_row < (u32) dims->nbxe[stage])
        return (u32) qp_in->idxe[stage][dims->nbue[stage] + equality_row];

    if (stage > 0)
        equality_row -= dims->nbxe[stage];
    assert(equality_row < (u32) dims->nge[stage]);
    return (u32) qp_in->idxe[stage]
        [dims->nbue[stage] + dims->nbxe[stage] + equality_row];
}

static inline void acados_daocp_store_signed_multiplier(
    struct blasfeo_dvec* lam, u32 nbg, u32 constraint_index, f64 multiplier)
{
    if (multiplier < 0.0)
        BLASFEO_DVECEL(lam, constraint_index) = -multiplier;
    else
        BLASFEO_DVECEL(lam, nbg + constraint_index) = multiplier;
}

/* Keep this threshold in sync with daocp_gaussian_elimination(). */
#define ACADOS_DAOCP_GE_ZERO_TOL 1e-7

static u32 acados_daocp_gaussian_elimination_record(
    f64* matrix, f64* row_tmp, f64* alphas, u32* pivots,
    u32 nrows, u32 ncols, u32 stride, u32 max_rank)
{
    u32 rank = 0;
    acados_size_t alpha_index = 0;

    for (u32 column = 0; column < ncols && rank < max_rank; column++)
    {
        u32 pivot_index = nrows;
        f64 pivot = ACADOS_DAOCP_GE_ZERO_TOL;
        for (u32 row = rank; row < nrows; row++)
        {
            f64 candidate = matrix[row*stride + column];
            if (DAOCP_ABS(candidate) > DAOCP_ABS(pivot))
            {
                pivot_index = row;
                pivot = candidate;
            }
        }
        if (pivot_index == nrows)
            continue;

        pivots[rank] = pivot_index;
        if (pivot_index != rank)
        {
            memcpy(row_tmp, matrix + rank*stride, stride*sizeof(f64));
            memcpy(matrix + rank*stride, matrix + pivot_index*stride, stride*sizeof(f64));
            memcpy(matrix + pivot_index*stride, row_tmp, stride*sizeof(f64));
        }

        for (u32 row = rank + 1; row < nrows; row++)
        {
            f64 alpha = matrix[row*stride + column] / pivot;
            alphas[alpha_index++] = alpha;
            for (u32 j = column; j < stride; j++)
                matrix[row*stride + j] -= alpha * matrix[rank*stride + j];
        }
        rank++;
    }

    return rank;
}

static void acados_daocp_record_equality_elimination(
    const ocp_qp_dims* dims, const daocp_qp* qp_native,
    daocp_workspace* wrk, ocp_qp_daocp_memory* mem)
{
    u32 N = (u32) dims->N;
    u32* nx = qp_native->dims.nx;
    u32* nu = qp_native->dims.nu;
    u32* ne = qp_native->dims.ne;
    f64* matrix = wrk->GEtmp;
    f64* H = wrk->H;
    f64* ABtmp = wrk->ABtmp;

    acados_size_t alpha_capacity, pivot_capacity;
    acados_daocp_equality_elimination_capacity(dims, &alpha_capacity, &pivot_capacity);
    acados_size_t alpha_base = alpha_capacity;
    acados_size_t pivot_base = pivot_capacity;
    acados_size_t capacity_rows = ne[N];

    u32 propagated_rows = ne[N];
    for (u32 row = 0; row < propagated_rows; row++)
        memcpy(H + row*nx[N], qp_native->Dx[N] + row*nx[N], nx[N]*sizeof(f64));

    for (i32 stage = (i32) N - 1; stage >= 0; stage--)
    {
        u32 nrows = ne[stage] + propagated_rows;
        u32 stride = nu[stage] + nx[stage] + 1;
        memset(matrix, 0, nrows*stride*sizeof(f64));

        for (u32 row = 0; row < ne[stage]; row++)
        {
            memcpy(matrix + row*stride,
                   qp_native->Du[stage] + row*nu[stage], nu[stage]*sizeof(f64));
            memcpy(matrix + row*stride + nu[stage],
                   qp_native->Dx[stage] + row*nx[stage], nx[stage]*sizeof(f64));
        }

        blasfeo_unpack_tran_dmat(
            nu[stage], nx[stage+1], qp_native->BAwt + stage, 0, 0,
            ABtmp, nx[stage+1]);
        daocp_fma_mm_nt(
            matrix + ne[stage]*stride, H, ABtmp,
            propagated_rows, nu[stage], nx[stage+1], stride);
        blasfeo_unpack_tran_dmat(
            nx[stage], nx[stage+1], qp_native->BAwt + stage, nu[stage], 0,
            ABtmp, nx[stage+1]);
        daocp_fma_mm_nt(
            matrix + ne[stage]*stride + nu[stage], H, ABtmp,
            propagated_rows, nx[stage], nx[stage+1], stride);

        capacity_rows += ne[stage];
        acados_size_t capacity_rank = DAOCP_MIN((acados_size_t) nu[stage], capacity_rows);
        acados_size_t stage_alpha_capacity =
            acados_daocp_num_elimination_alphas(capacity_rows, capacity_rank);
        alpha_base -= stage_alpha_capacity;
        pivot_base -= capacity_rank;

        u32 rank = acados_daocp_gaussian_elimination_record(
            matrix, matrix + nrows*stride,
            mem->equality_elimination_alphas + alpha_base,
            mem->equality_elimination_pivots + pivot_base,
            nrows, nu[stage], stride, nu[stage]);
        assert(rank == wrk->rho[stage]);

        for (u32 row = rank; row < nrows; row++)
            memcpy(H + (row-rank)*nx[stage],
                   matrix + row*stride + nu[stage], nx[stage]*sizeof(f64));
        propagated_rows = nrows - rank;
    }

    assert(alpha_base == 0);
    assert(pivot_base == 0);
}

static void acados_daocp_recover_equality_duals(
    const ocp_qp_in* qp_in, ocp_qp_out* qp_out,
    daocp_workspace* wrk, ocp_qp_daocp_memory* mem)
{
    const ocp_qp_dims* dims = qp_in->dim;
    u32 N = (u32) dims->N;
    u32* ne = mem->qp.dims.ne;
    f64* dual = wrk->tmp1;

    acados_size_t alpha_base = 0;
    acados_size_t pivot_base = 0;
    acados_size_t capacity_rows = 0;
    for (u32 stage = 0; stage <= N; stage++)
        capacity_rows += ne[stage];
    if (capacity_rows == 0)
        return;

    acados_daocp_record_equality_elimination(dims, &mem->qp, wrk, mem);

    u32 propagated_duals = wrk->nH0;
    memset(dual, 0, propagated_duals*sizeof(f64));

    for (u32 stage = 0; stage < N; stage++)
    {
        u32 rank = wrk->rho[stage];
        u32 nrows = rank + propagated_duals;
        memmove(dual + rank, dual, propagated_duals*sizeof(f64));
        for (u32 row = 0; row < rank; row++)
            dual[row] = BLASFEO_DVECEL(wrk->eta_lqr + stage, row) + wrk->eta[stage][row];

        for (i32 pivot_row = (i32) rank - 1; pivot_row >= 0; pivot_row--)
        {
            acados_size_t row_alpha_offset = alpha_base
                + (acados_size_t) pivot_row * (2*nrows - pivot_row - 1) / 2;
            f64 pivot_dual = dual[pivot_row];
            for (u32 row = (u32) pivot_row + 1; row < nrows; row++)
                pivot_dual -= mem->equality_elimination_alphas
                    [row_alpha_offset + row - pivot_row - 1] * dual[row];
            dual[pivot_row] = pivot_dual;

            u32 original_pivot = mem->equality_elimination_pivots[pivot_base + pivot_row];
            f64 tmp = dual[pivot_row];
            dual[pivot_row] = dual[original_pivot];
            dual[original_pivot] = tmp;
        }

        for (u32 row = 0; row < ne[stage]; row++)
        {
            u32 constraint_index = acados_daocp_original_equality_index(qp_in, stage, row);
            acados_daocp_store_signed_multiplier(
                qp_out->lam + stage, dims->nb[stage] + dims->ng[stage],
                constraint_index, dual[row]);
        }

        propagated_duals = nrows - ne[stage];
        memmove(dual, dual + ne[stage], propagated_duals*sizeof(f64));

        acados_size_t capacity_rank = DAOCP_MIN((acados_size_t) dims->nu[stage], capacity_rows);
        alpha_base += acados_daocp_num_elimination_alphas(capacity_rows, capacity_rank);
        pivot_base += capacity_rank;
        capacity_rows -= ne[stage];
    }

    assert(propagated_duals == ne[N]);
    for (u32 row = 0; row < ne[N]; row++)
    {
        u32 constraint_index = acados_daocp_original_equality_index(qp_in, N, row);
        acados_daocp_store_signed_multiplier(
            qp_out->lam + N, dims->nb[N] + dims->ng[N],
            constraint_index, dual[row]);
    }
}

static void acados_daocp_compute_state_cost_gradient(
    const ocp_qp_in* qp_in, const ocp_qp_out* qp_out,
    u32 stage, struct blasfeo_dvec* gradient)
{
    const ocp_qp_dims* dims = qp_in->dim;
    u32 nx = dims->nx[stage];
    u32 nu = dims->nu[stage];

    blasfeo_dveccp(nx, qp_in->rqz + stage, nu, gradient, 0);
    if (nu > 0)
        blasfeo_dgemv_n(
            nx, nu, 1.0, qp_in->RSQrq + stage, nu, 0,
            qp_out->ux + stage, 0, 1.0, gradient, 0, gradient, 0);
    blasfeo_dsymv_l(
        nx, 1.0, qp_in->RSQrq + stage, nu, nu,
        qp_out->ux + stage, nu, 1.0, gradient, 0, gradient, 0);
}

static void acados_daocp_add_state_constraint_dual_gradient(
    const ocp_qp_in* qp_in, ocp_qp_out* qp_out,
    u32 stage, struct blasfeo_dvec* gradient)
{
    const ocp_qp_dims* dims = qp_in->dim;
    u32 nx = dims->nx[stage];
    u32 nu = dims->nu[stage];
    u32 nb = dims->nb[stage];
    u32 ng = dims->ng[stage];
    u32 nbg = nb + ng;
    struct blasfeo_dvec* lam = qp_out->lam + stage;

    for (u32 i = 0; i < nb; i++)
    {
        u32 variable_index = qp_in->idxb[stage][i];
        if (variable_index >= nu)
            BLASFEO_DVECEL(gradient, variable_index - nu) +=
                BLASFEO_DVECEL(lam, nbg + i) - BLASFEO_DVECEL(lam, i);
    }

    if (ng > 0)
    {
        // t is overwritten by ocp_qp_compute_t() immediately after dual recovery.
        blasfeo_daxpy(
            ng, -1.0, lam, nb, lam, nbg + nb, qp_out->t + stage, 0);
        blasfeo_dgemv_n(
            nx, ng, 1.0, qp_in->DCt + stage, nu, 0,
            qp_out->t + stage, 0, 1.0, gradient, 0, gradient, 0);
    }
}

static void acados_daocp_compute_dual_solution(
    const ocp_qp_in* qp_in, ocp_qp_out* qp_out,
    daocp_workspace* wrk, ocp_qp_daocp_memory* mem)
{
    const ocp_qp_dims* dims = qp_in->dim;
    u32 N = (u32) dims->N;

    for (u32 stage = 0; stage <= N; stage++)
        blasfeo_dvecse(2*(dims->nb[stage] + dims->ng[stage]), 0.0, qp_out->lam + stage, 0);

    // DA-OCP stores only active inequality multipliers. Its signed convention is
    // negative for lower bounds and positive for upper bounds, matching
    // lam_upper - lam_lower in the acados stationarity equations.
    for (u32 i = 0; i < wrk->as.n_active; i++) {
        const daocp_constraint* constraint = wrk->as.xi2con + i;
        u32 stage = constraint->t;
        u32 constraint_index = 
            acados_daocp_original_inequality_index(qp_in, stage, constraint);
        u32 nbg = dims->nb[stage] + dims->ng[stage];
        if (constraint->is_upper)
            BLASFEO_DVECEL(qp_out->lam + stage, nbg + constraint_index) =
                DAOCP_MAX(wrk->xi[i], 0.0);
        else
            BLASFEO_DVECEL(qp_out->lam + stage, constraint_index) =
                DAOCP_MAX(-wrk->xi[i], 0.0);
    }

    acados_daocp_recover_equality_duals(qp_in, qp_out, wrk, mem);

    // Recover dynamics multipliers with the costate recursion. This is the
    // BLASFEO analogue of the pi recursion in example_d_riccati_recursion.c.
    if (N > 0) {
        acados_daocp_compute_state_cost_gradient(qp_in, qp_out, N, qp_out->pi + N - 1);
        acados_daocp_add_state_constraint_dual_gradient(qp_in, qp_out, N, qp_out->pi + N - 1);

        for (i32 stage = (i32) N - 1; stage > 0; stage--)
        {
            acados_daocp_compute_state_cost_gradient(
                qp_in, qp_out, (u32) stage, qp_out->pi + stage - 1);
            blasfeo_dgemv_n(
                dims->nx[stage], dims->nx[stage+1], 1.0,
                qp_in->BAbt + stage, dims->nu[stage], 0,
                qp_out->pi + stage, 0, 1.0,
                qp_out->pi + stage - 1, 0, qp_out->pi + stage - 1, 0);
            acados_daocp_add_state_constraint_dual_gradient(
                qp_in, qp_out, (u32) stage, qp_out->pi + stage - 1);
        }
    }

    // x0 is removed from DA-OCP's optimization variables. Recover the signed
    // multipliers of its fixing bounds from the remaining x0 stationarity rows.
    if (dims->nbxe[0] > 0)
    {
        struct blasfeo_dvec* residual = &wrk->costate0;
        acados_daocp_compute_state_cost_gradient(qp_in, qp_out, 0, residual);
        if (N > 0)
            blasfeo_dgemv_n(
                dims->nx[0], dims->nx[1], 1.0,
                qp_in->BAbt, dims->nu[0], 0,
                qp_out->pi, 0, 1.0, residual, 0, residual, 0);
        acados_daocp_add_state_constraint_dual_gradient(qp_in, qp_out, 0, residual);

        for (u32 i = 0; i < (u32) dims->nbxe[0]; i++)
        {
            u32 bound_index = qp_in->idxe[0][dims->nbue[0] + i];
            u32 state_index = qp_in->idxb[0][bound_index] - dims->nu[0];
            acados_daocp_store_signed_multiplier(
                qp_out->lam, dims->nb[0] + dims->ng[0],
                bound_index, -BLASFEO_DVECEL(residual, state_index));
        }
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
    for (u32 t=0; t<=qp_in->dim->N; ++t) 
        if (qp_in->dim->ns[t] != 0) {
            printf("\nDAOCP cannot support the %d slack variables at stage %u "
                   "(nx=%d, nu=%d, nb=%d, ng=%d, nbxe=%d, nbue=%d, nge=%d).\n",
                   qp_in->dim->ns[t], t, qp_in->dim->nx[t], qp_in->dim->nu[t],
                   qp_in->dim->nb[t], qp_in->dim->ng[t], qp_in->dim->nbxe[t],
                   qp_in->dim->nbue[t], qp_in->dim->nge[t]);
            return ACADOS_QP_FAILURE;
        }
    // Conversion of data structures
    daocp_qp* qp_native = &mem->qp;
    // (Shadow) copy dynamics and cost
    qp_native->BAwt = qp_in->BAbt;
    qp_native->RSQrq = qp_in->RSQrq;
    for (int t=0; t<dim->N; ++t)
        blasfeo_drowin(dim->nx[t+1], 1.0, qp_in->b+t, 0, qp_in->BAbt+t, dim->nu[t]+dim->nx[t], 0);
    for (int t=0; t<=dim->N; ++t)
        blasfeo_drowin(dim->nu[t]+dim->nx[t], 1.0, qp_in->rqz+t, 0, qp_in->RSQrq+t, dim->nu[t]+dim->nx[t], 0);
    daocp_workspace* wrk = (daocp_workspace*) mem->workspace;
    wrk->dims = &qp_native->dims;
    acados_daocp_process_constraints(qp_in, dim, qp_native, wrk);
    if (opts->first_run || !opts->warm_start) acados_daocp_init_workspace(wrk);
    // Set solution pointers
    mem->sol.ux = qp_out->ux;

    // Interface work done
    info->interface_time = acados_toc(&interface_timer);

    // TODO: Handle first_solve != 0
    acados_tic(&qp_timer);
    // The following calls before daocp_solve are effectively
    // part of the solver.
    acados_tic(&solver_call_timer); 
    daocp_solve_riccati(wrk, qp_native);
    daocp_solve_lqr(wrk, qp_native);
    // We check whether we can use the active set information from 
    // the previous solve.
    if (!opts->first_run) {
        // Recompute cholesky of dual hessian, detecting singularity
        u32 need_reset = daocp_compute_chol_from_scratch(wrk, qp_native);
        // Compute dual minimizer and check dual feasibility
        if (!need_reset) {
            daocp_solve_dual_eqcon_qp(wrk);
            need_reset = !daocp_is_dual_feasible(wrk->p, wrk->xi_sign, wrk->as.n_active);
            // If the dual minimizer is feasible, we can keep the working set
            if (!need_reset) memcpy(wrk->xi, wrk->p, wrk->as.n_active*sizeof(f64));
        }
        if (need_reset) daocp_reset_working_set(wrk);
    }
    daocp_solve(&opts->daocp_opts, qp_native, wrk, &mem->sol);
    opts->first_run = 0;
    mem->time_qp_solver_call = acados_toc(&solver_call_timer);

    /* fill qp_out */
    acados_daocp_compute_dual_solution(qp_in, qp_out, wrk, mem);
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
