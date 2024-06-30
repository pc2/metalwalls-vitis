#pragma once
#include <cmath>
#include <cstdint>

void MWPot_lrPotCPU(double param_a, double param_alpha, double param_b, double param_c, int32_t param_kxMax,
                    int32_t param_kyMax, int32_t param_kzMax, int32_t param_numWall1, int32_t param_numWall2,
                    double param_rksqmax, double *instream_cossinSkx, double *instream_cossinSky,
                    double *instream_cossinSkz, double *qWall, double *outstream_MWPotCPU);

void MWPot_srPotWallCPU(double param_a, double param_alpha, double param_b, double param_eta, int32_t param_numWall1,
                        double param_rcutsq, double *xyz, double *q, double *outstream_MWPotCPU);

void MWPot_k0PotCPU(double param_a, double param_b, double param_alpha, int32_t param_numWall1, int32_t param_numWall2,
                    double *z, double *q, double *outstream_MWPotCPU);

void init_xyz(int num, double *xyz, double *x, double *y, double *z);

void init_instream_cossinSk(int32_t numWall, int32_t kMax, int32_t ixyz, double kdim, double *xyzWall,
                            double *instream_cossinSk);