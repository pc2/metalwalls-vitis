const int uf = 16;
const int uf_io = 8;

const int num_max = 42508 + uf;

const int c_n = 2496;
const int SHIFT_REG_SIZE = 8 + 1;

const int max_iterations_max = 500;

void vdot(const int size, const double *a, const double *b, double *out)
{

    double shift_reg[SHIFT_REG_SIZE];
    double sum = 0.;

shift_init:
    for (int i = 0; i < SHIFT_REG_SIZE; i++)
    {
#pragma HLS unroll
        shift_reg[i] = 0;
    }

comp:
    for (int in = 0; in < size; in += uf)
    {
#pragma HLS LOOP_TRIPCOUNT min = c_n / uf max = c_n / uf

        double new_value = 0.0;
        for (int i = 0; i < uf; i++)
        {
#pragma HLS unroll
            new_value += a[in + i] * b[in + i];
        }
        shift_reg[SHIFT_REG_SIZE - 1] = shift_reg[0] + new_value;

    shift_operation:
        for (int i = 0; i < SHIFT_REG_SIZE - 1; i++)
        {
#pragma HLS unroll
            shift_reg[i] = shift_reg[i + 1];
        }
    }

final_sum:
    for (int i = 0; i < SHIFT_REG_SIZE - 1; i++)
    {
#pragma HLS unroll
        sum += shift_reg[i];
    }

writeC:
    *out = sum;
}

void daxpy(double *out, const int size, const double *y, const double a, const double *x)
{
    double local_out[num_max];
    for (int i = 0; i < size; i += uf)
    {
#pragma HLS loop_tripcount min = c_n / uf max = c_n / uf
#pragma HLS pipeline II = 1
        for (int ii = 0; ii < uf; ii++)
        {
#pragma HLS unroll
            local_out[i + ii] = y[i + ii] + a * x[i + ii];
        }
    }
    for (int i = 0; i < size; i += uf)
    {
#pragma HLS loop_tripcount min = c_n / uf max = c_n / uf
#pragma HLS pipeline II = 1
        for (int ii = 0; ii < uf; ii++)
        {
#pragma HLS unroll
            out[i + ii] = local_out[i + ii];
        }
    }
}

extern "C"
{
    void kernel_cg(const int num_atoms, const int iter, const double selfPotFactor, const double rsold,
                   const double *b_cg, const double *q_in, const double *res_in, const double *x_cg_in,
                   double *x_cg_out, double *q_out, double *res_out, double *rsnew, double *Ap)
    {
#pragma HLS INTERFACE m_axi port = b_cg offset = slave bundle = gmem1
#pragma HLS INTERFACE m_axi port = q_in offset = slave bundle = gmem2
#pragma HLS INTERFACE m_axi port = res_in offset = slave bundle = gmem3
#pragma HLS INTERFACE m_axi port = x_cg_in offset = slave bundle = gmem4
#pragma HLS INTERFACE m_axi port = x_cg_out offset = slave bundle = gmem8
#pragma HLS INTERFACE m_axi port = q_out offset = slave bundle = gmem5
#pragma HLS INTERFACE m_axi port = res_out offset = slave bundle = gmem6
#pragma HLS INTERFACE m_axi port = rsnew offset = slave bundle = gmem7
#pragma HLS INTERFACE m_axi port = Ap offset = slave bundle = gmem9

        // input buffer
        double b_cg_loc[num_max];
        double Ap_loc[num_max];
        double q_in_loc[num_max];
        double res_in_loc[num_max];
        double x_cg_in_loc[num_max];
        double rsold_loc = rsold;

        // read in
    read:
        for (int i = 0; i < num_atoms; i += uf_io)
        {
#pragma HLS loop_tripcount min = c_n / uf_io max = c_n / uf_io
#pragma HLS pipeline II = 1
            for (int ii = 0; ii < uf_io; ii++)
            {
#pragma HLS loop_tripcount min = uf_io max = uf_io
#pragma HLS unroll
                b_cg_loc[i + ii] = b_cg[i + ii];
                q_in_loc[i + ii] = q_in[i + ii];
                res_in_loc[i + ii] = res_in[i + ii];
                x_cg_in_loc[i + ii] = x_cg_in[i + ii];
                Ap_loc[i + ii] = Ap[i + ii];
            }
        }
        // output buffer
        double rsnew_loc;
        double q_out_loc[num_max];
        double res_out_loc[num_max];
        double x_cg_out_loc[num_max];

    selfPot:
        daxpy(Ap_loc, num_atoms, Ap_loc, -selfPotFactor, q_in_loc);

        if (iter == 0)
        {
            // Setup initial residual
            for (int i = 0; i < num_atoms; i += uf)
            {
#pragma HLS loop_tripcount min = c_n / uf max = c_n / uf
#pragma HLS pipeline II = 1
                for (int ii = 0; ii < uf; ii++)
                {
#pragma HLS unroll
                    res_out_loc[i + ii] = b_cg_loc[i + ii] - Ap_loc[i + ii];
                    q_out_loc[i + ii] = res_out_loc[i + ii];
                }
            }

            vdot(num_atoms, res_out_loc, res_out_loc, &rsnew_loc);
        }
        else
        {
            double pAp;
            vdot(num_atoms, q_in_loc, Ap_loc, &pAp);

            double alpha_cg = rsold_loc / pAp;

            daxpy(x_cg_out_loc, num_atoms, x_cg_in_loc, alpha_cg, q_in_loc);

            daxpy(res_out_loc, num_atoms, res_in_loc, -alpha_cg, Ap_loc);

            vdot(num_atoms, res_out_loc, res_out_loc, &rsnew_loc);

            /// Setup for next iteration
            double beta = rsnew_loc / rsold_loc;

            daxpy(q_out_loc, num_atoms, res_out_loc, beta, q_in_loc);
        }
    write:
        for (int i = 0; i < num_atoms; i += uf_io)
        {
#pragma HLS loop_tripcount min = c_n / uf_io max = c_n / uf_io
#pragma HLS pipeline II = 1
            for (int ii = 0; ii < uf_io; ii++)
            {
#pragma HLS unroll
                q_out[i + ii] = q_out_loc[i + ii];
                res_out[i + ii] = res_out_loc[i + ii];
                x_cg_out[i + ii] = x_cg_out_loc[i + ii];
                if (i == 0)
                {
                    rsnew[0] = rsnew_loc;
                }
            }
        }
    } // end of the kernel
} // end of the extern C
