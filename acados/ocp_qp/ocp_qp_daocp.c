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
    opts->first_run = 1;
    opts->daocp_opts.max_iter = 1000;
    opts->daocp_opts.selection = DAOCP_SELECT_GREEDY;
    return;
}



void ocp_qp_daocp_opts_update(void *config_, void *dims_, void *opts_)
{
    return;
}

void ocp_qp_daocp_opts_set(void *config_, void *opts_, const char *field, void *value)
{
    ocp_qp_daocp_opts* opts = opts_;

    if (!strcmp(field, "max_iter"))
    {
        int *tmp_ptr = value;
        opts->daocp_opts.max_iter = *tmp_ptr;
    }
    else if (!strcmp(field, "print_level"))
    {
        int* print_level = (int *) value;
        opts->print_level = *print_level;
    }
    else if (!strcmp(field, "selection_strategy"))
    {
        int code = *((int*) value);
        if (code == 0) opts->daocp_opts.selection = DAOCP_SELECT_GREEDY;
        if (code == 1) opts->daocp_opts.selection = DAOCP_SELECT_MOST_VIOLATED;
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
    printf("\nerror: ocp_qp_daocp_opts_get: not implemented for field %s\n", field);
    exit(1);
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
    int* nbx = dims->nbx;
    int* nbu = dims->nbu;
    int* ng = dims->ng; 
    int* ns = dims->ns;
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
    size += 2*blasfeo_memsize_dvec(max_nx); // costate0, costate1

    return size;
}

static inline void* assign_ptr_vec(char** ptr_vec, char* ptr, int* size_vec1, int* size_vec2, int* size_vec3, int unit_size, int n)
{
    for (u32 i=0; i<n; ++i) {
        ptr_vec[i] = ptr;
        ptr += (size_vec1[i] - 
            (size_vec2 ? size_vec2[i] : 0) +
            (size_vec3 ? size_vec3[i] : 0))*unit_size;
    }
    return ptr;
}

static inline void* assign_ptr_mat(
    char** ptr_mat, char* ptr, 
    int* size_vec1, int* size_vec2, int* size_vec3, int* size_vec4, 
    int* size_vec5, int unit_size, int n)
{
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
    mem->workspace = (c_ptr += sizeof(ocp_qp_daocp_memory)); 
    daocp_workspace* wrk = mem->workspace;

    int N = dims->N;
    int* nx = dims->nx;
    int* nu = dims->nu;
    int* nb = dims->nb;
    int* nbu = dims->nbu;
    int* nbx = dims->nbx;
    int* ng = dims->ng; 
    int* ns = dims->ns;
    int* nge = dims->nge;
    int* nbue = dims->nbue;
    int* nbxe = dims->nbxe;

    // daocp_qp data
    mem->qp.dims.nbu = (c_ptr += sizeof(daocp_workspace));
    mem->qp.dims.nbx = (c_ptr += N*sizeof(u32));
    mem->qp.dims.ng = (c_ptr += (N+1)*sizeof(u32));
    mem->qp.dims.ne = (c_ptr += (N+1)*sizeof(u32));
    mem->qp.x0 = (c_ptr += (N+1)*sizeof(u32));
    mem->qp.lbu = (c_ptr += nx[0]*sizeof(f64));
    mem->qp.lbx = (c_ptr += N*sizeof(f64*));
    wrk->lbu_wrk = (c_ptr += (N+1)*sizeof(f64*));
    wrk->lbx_wrk = (c_ptr += N*sizeof(f64*));
    mem->qp.ubu = (c_ptr += (N+1)*sizeof(f64*));
    mem->qp.ubx = (c_ptr += N*sizeof(f64*));
    wrk->ubu_wrk = (c_ptr += (N+1)*sizeof(f64*));
    wrk->ubx_wrk = (c_ptr += N*sizeof(f64*));
    mem->qp.cl = (c_ptr += (N+1)*sizeof(f64*));
    mem->qp.cu = (c_ptr += (N+1)*sizeof(f64*));
    wrk->lg_wrk = (c_ptr += (N+1)*sizeof(f64*));
    wrk->ug_wrk = (c_ptr += (N+1)*sizeof(f64*));
    mem->qp.d = (c_ptr += (N+1)*sizeof(f64*));
    mem->qp.Cu = (c_ptr += (N+1)*sizeof(f64*));
    mem->qp.Cx = (c_ptr += N*sizeof(f64*));
    mem->qp.Du = (c_ptr += (N+1)*sizeof(f64*));
    mem->qp.Dx = (c_ptr += N*sizeof(f64*));
    mem->qp.idxbu = (c_ptr += (N+1)*sizeof(f64*));
    mem->qp.idxbx = (c_ptr += N*sizeof(u32*));
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
    c_ptr = assign_ptr_vec(mem->qp.cl, c_ptr, ng, nge, 0, sizeof(f64), N+1);
    c_ptr = assign_ptr_vec(wrk->lg_wrk, c_ptr, ng, nge, 0, sizeof(f64), N+1);
    c_ptr = assign_ptr_vec(mem->qp.cu, c_ptr, ng, nge, 0, sizeof(f64), N+1);
    c_ptr = assign_ptr_vec(wrk->ug_wrk, c_ptr, ng, nge, 0, sizeof(f64), N+1);
    mem->qp.d[0] = c_ptr; c_ptr += (nbue[0]+nge[0])*sizeof(f64);
    for (u32 t=1; t<=N; ++t) {
        mem->qp.d[t] = c_ptr;
        c_ptr += (nbue[t]+nbxe[t]+nge[t])*sizeof(f64);
    }
    c_ptr = assign_ptr_mat(mem->qp.Cu, c_ptr, ng, nge, 0, 0, nu, sizeof(f64), N);
    mem->qp.Cx[0]=0; c_ptr = assign_ptr_mat(mem->qp.Cx+1, c_ptr, ng+1, nge+1, 0, 0, nx+1, sizeof(f64), N);
    mem->qp.Du[0] = c_ptr; c_ptr += (nge[0]+nbue[0])*nu[0]*sizeof(f64);
    c_ptr = assign_ptr_mat(mem->qp.Du+1, c_ptr, nge+1, 0, nbue+1, nbxe+1, nu+1, sizeof(f64), N-1);
    mem->qp.Dx[0]= c_ptr; c_ptr += (nge[0]+nbue[0])*nx[0]*sizeof(f64); 
    c_ptr = assign_ptr_mat(mem->qp.Dx+1, c_ptr, nge+1, 0, nbue+1, nbxe+1, nx+1, sizeof(f64), N);

    // daocp_workspace data
    int max_nx = 0; int max_nu = 0; int neq = 0; int nin = 0;
    for (u32 t=0; t<=N; ++t) if (nx[t] > max_nx) max_nx = nx[t];
    for (u32 t=0; t<N; ++t) if (nu[t] > max_nu) max_nu = nu[t];
    for (u32 t=0; t<=N; ++t) neq += nge[t] + nbue[t] + nbxe[t];
    neq -= nbxe[0];
    for (u32 t=0; t<=N; ++t) nin += ng[t]+nb[t]-nge[t]-nbue[t]-nbxe[t];

    wrk->P = c_ptr; c_ptr+=N*sizeof(struct blasfeo_dmat);
    wrk->Ku = c_ptr; c_ptr+=N*sizeof(struct blasfeo_dmat);
    wrk->Ke = c_ptr; c_ptr+=N*sizeof(struct blasfeo_dmat);
    wrk->Luu = c_ptr; c_ptr+=N*sizeof(struct blasfeo_dmat);
    wrk->Lue = c_ptr; c_ptr+=N*sizeof(struct blasfeo_dmat);
    wrk->Lee = c_ptr; c_ptr+=N*sizeof(struct blasfeo_dmat);
    wrk->ux_lqr = c_ptr; c_ptr+=(N+1)*sizeof(struct blasfeo_dvec);
    wrk->eta_lqr = c_ptr; c_ptr+=N*sizeof(struct blasfeo_dvec);
    wrk->b = c_ptr; c_ptr+=N*sizeof(struct blasfeo_dvec);

    wrk->u = c_ptr; c_ptr += N*sizeof(f64*);
    wrk->x = c_ptr; c_ptr += (N+1)*sizeof(f64*);
    wrk->eta = c_ptr; c_ptr += N*sizeof(f64*);
    wrk->contypes = c_ptr; c_ptr += (N+1)*sizeof(daocp_constraint_type*);
    wrk->as.constraint_status = c_ptr; c_ptr += (N+1)*sizeof(u32*);
    c_ptr = assign_ptr_vec(wrk->u, c_ptr, nu, 0, 0, sizeof(f64), N);
    c_ptr = assign_ptr_vec(wrk->eta, c_ptr, nu, 0, 0, sizeof(f64), N);
    wrk->x[0]=0; c_ptr = assign_ptr_vec(wrk->x+1, c_ptr, nx+1, 0, 0, sizeof(f64), N);
    for (u32 t=0; t<=N; ++t) {
        wrk->contypes[t] = c_ptr;
        c_ptr += (ng[t]-nge[t])*sizeof(daocp_constraint_type);
    }
    for (u32 t=0; t<=N; ++t) {
        wrk->as.constraint_status[t] = c_ptr;
        c_ptr += (ng[t]+nb[t]-nbue[t]-nbxe[t]-nge[t])*sizeof(u32);
    }
    wrk->cnu = c_ptr; c_ptr += N*sizeof(u32);
    wrk->rho = c_ptr; c_ptr += N*sizeof(u32);
    wrk->crho = c_ptr; c_ptr += N*sizeof(u32);
    
    int tot_nx = nx[N];
    int tot_nu = 0;
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

        tot_nx += nx[t]; tot_nu += nu[t];
    }
    blasfeo_create_dvec(nx[N], wrk->ux_lqr+N, c_ptr);
    c_ptr += blasfeo_memsize_dvec(nx[N]);
    
    u32 W_stride = DAOCP_MIN(tot_nu, nin)+1;
    wrk->as.xi2con = c_ptr; c_ptr += W_stride*sizeof(daocp_constraint);
    wrk->xi = c_ptr; c_ptr+=W_stride*sizeof(f64);
    wrk->p = c_ptr; c_ptr+=W_stride*sizeof(f64);
    wrk->dual_linear = c_ptr; c_ptr+=W_stride*sizeof(f64);
    wrk->xi_sign = c_ptr; c_ptr+=W_stride*sizeof(u32);
    wrk->Ld = c_ptr; c_ptr+=(W_stride+1)*W_stride*sizeof(f64);
    wrk->Mu = c_ptr; c_ptr+=W_stride*tot_nu*sizeof(f64);
    wrk->Me = c_ptr; c_ptr+=W_stride*tot_nu*sizeof(f64);
    wrk->H = c_ptr; c_ptr+=neq*max_nx*sizeof(f64);
    wrk->h = c_ptr; c_ptr+=neq*sizeof(f64);
    wrk->tmp1 = c_ptr; c_ptr+=neq*sizeof(f64);
    wrk->GEtmp = c_ptr; c_ptr+=(neq+1)*(max_nx+max_nu+1)*sizeof(f64);
    wrk->ABtmp = c_ptr; c_ptr+=max_nx*DAOCP_MAX(max_nx, max_nu)*sizeof(f64);
    blasfeo_create_dmat(DAOCP_MAX(max_nx, max_nu), max_nx, &wrk->tmp2, c_ptr);
    c_ptr += blasfeo_memsize_dmat(DAOCP_MAX(max_nu, max_nx), max_nx);
    blasfeo_create_dvec(max_nx, &wrk->costate0, c_ptr);
    c_ptr += blasfeo_memsize_dvec(max_nx);
    blasfeo_create_dvec(max_nx, &wrk->costate1, c_ptr);
    c_ptr += blasfeo_memsize_dvec(max_nx);

    assert((char *) raw_memory + ocp_qp_daocp_memory_calculate_size(config_, dims, opts_) == c_ptr);

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
    // ocp_qp_in *qp_in = qp_in_;
    // reset memory
    printf("acados: reset daocp_mem not implemented.\n");
    exit(1);
}



/************************************************
 * workspace
 ************************************************/

acados_size_t ocp_qp_daocp_workspace_calculate_size(void *config_, void *dims_, void *opts_)
{
    return 0;
}

static u32 acados_daocp_contains_index(const u32 *indices, u32 offset, u32 count, u32 index)
{
    for (u32 i = 0; i < count; i++)
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
    if (qp_in->dim->nx[0] != qp_in->dim->nbxe[0]) {
        printf("\nDAOCP can only handle problems with fixed initial state.\n");
        exit(1);
    }
    for (u32 t=0; t<=qp_in->dim->N; ++t) 
        if (qp_in->dim->ns[t] != 0) {
            printf("\nDAOCP cannot support slack variables yet.\n");
            exit(1);
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
    // TODO: Handle first_run != 0
    acados_daocp_process_constraints(qp_in, dim, qp_native, wrk);
    if (opts->first_run) acados_daocp_init_workspace(wrk);
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
    if (!first_run) {
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
    mem->time_qp_solver_call = acados_toc(&solver_call_timer);

    /* fill qp_out */
    ocp_qp_compute_t(qp_in, qp_out);

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

    return acados_status;
}



void ocp_qp_daocp_eval_adj_sens(void *config_, void *param_qp_in_, void *seed, void *sens_qp_out_, void *opts_, void *mem_, void *work_)
{
    printf("\nerror: ocp_qp_clarabel_eval_adj_sens: not implemented yet\n");
    exit(1);
}

void ocp_qp_daocp_eval_forw_sens(void *config_, void *param_qp_in_, void *seed, void *sens_qp_out_, void *opts_, void *mem_, void *work_)
{
    printf("\nerror: ocp_qp_clarabel_eval_forw_sens: not implemented yet\n");
    exit(1);
}

void ocp_qp_daocp_solver_get(void *config_, void *qp_in_, void *qp_out_, void *opts_, void *mem_, const char *field, int stage, void* value, int size1, int size2)
{
    printf("\nerror: ocp_qp_clarabel_solver_get: not implemented yet\n");
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
