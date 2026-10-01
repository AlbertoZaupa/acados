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


#ifndef ACADOS_OCP_QP_OCP_QP_DAOCP_H_
#define ACADOS_OCP_QP_OCP_QP_DAOCP_H_

#ifdef __cplusplus
extern "C" {
#endif

#include "daocp/include/daocp.h"

// acados
#include "acados/ocp_qp/ocp_qp_common.h"
#include "acados/utils/types.h"

typedef struct ocp_qp_daocp_opts_
{
    daocp_args daocp_opts;
    int print_level;
    int warm_start;
    int first_run;
    int detect_input_bound_equalities;
    double equality_detection_tolerance;
} ocp_qp_daocp_opts;


typedef struct ocp_qp_daocp_memory_
{
    daocp_qp qp;
    daocp_sol sol;
    void* workspace;
    // All hard input equalities, both explicit and automatically detected.
    u32** detected_input_equalities;
    u32* num_detected_input_equalities;
    double time_qp_solver_call;
} ocp_qp_daocp_memory;

acados_size_t ocp_qp_daocp_opts_calculate_size(void *config, void *dims);
//
void *ocp_qp_daocp_opts_assign(void *config, void *dims, void *raw_memory);
//
void ocp_qp_daocp_opts_initialize_default(void *config, void *dims, void *opts_);
//
void ocp_qp_daocp_opts_update(void *config, void *dims, void *opts_);
//
void ocp_qp_daocp_opts_set(void *config, void *opts_, const char *field, void *value);
//
void ocp_qp_daocp_opts_get(void *config, void *opts_, const char *field, void *value);
//
acados_size_t ocp_qp_daocp_memory_calculate_size(void *config, void *dims, void *opts_);
//
void *ocp_qp_daocp_memory_assign(void *config, void *dims, void *opts_, void *raw_memory);
//
void ocp_qp_daocp_memory_get(void *config, void *mem_, const char *field, void *value);
//
acados_size_t ocp_qp_daocp_workspace_calculate_size(void *config, void *dims, void *opts_);
//
/* Soft constraints require independent slacks with Zl == Zu > 0,
 * zl == zu >= 0, zero slack lower bounds and ordered, unmasked bounds.
 * Initial-state fixing bounds must be hard. Soft QPs reset the active set.
 * Primal slacks are returned; pi and lam are zero placeholders. */
int ocp_qp_daocp(void *config, void *qp_in, void *qp_out, void *opts_, void *mem_, void *work_);
//
void ocp_qp_daocp_memory_reset(void *config_, void *qp_in_, void *qp_out_, void *opts_, void *mem_, void *work_);
//
void ocp_qp_daocp_solver_get(void *config_, void *qp_in_, void *qp_out_, void *opts_, void *mem_, const char *field, int stage, void* value, int size1, int size2);
//
void ocp_qp_daocp_config_initialize_default(void *config);

#ifdef __cplusplus
} /* extern "C" */
#endif

#endif  // ACADOS_OCP_QP_OCP_QP_DAOCP_H_
