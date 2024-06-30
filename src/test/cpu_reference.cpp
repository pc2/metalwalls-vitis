#include "cpu_reference.hpp"

void MWPot_lrPotCPU(double param_a, double param_alpha, double param_b, double param_c, int32_t param_kxMax,
                    int32_t param_kyMax, int32_t param_kzMax, int32_t param_numWall1, int32_t param_numWall2,
                    double param_rksqmax, double *instream_cossinSkx, double *instream_cossinSky,
                    double *instream_cossinSkz, double *qWall, double *outstream_MWPotCPU)
{

    int32_t numWall = param_numWall1 + param_numWall2;
    double alphasq = param_alpha * param_alpha;
    double alphaconst = -1.0 / (4.0 * alphasq);
    double zweight = 2.0 * M_PI / param_c;
    double lrPotFactor = 4.0 / (param_a * param_b) * zweight; //(4.0 * M_PI) / (param_a*param_b*param_c); ->
                                                              // it's actually times 8, not 4

    /*printf("CPU %d %d %d %d  %+4.15e  %+4.15e %+4.15e %+4.15e %+4.15e
       %+4.15e\n", numWall,  param_kxMax,  param_kyMax, param_kzMax, alphaconst,
       param_rksqmax, 2 * M_PI / param_a,2 * M_PI / param_b, 2 * M_PI / param_c,
       lrPotFactor);//*/
    for (int32_t i = 0; i < numWall; ++i)
    {
        outstream_MWPotCPU[i] = 0.0;
    }
    // int count=0;
    for (int32_t kx = 0; kx <= param_kxMax; ++kx)
    {
        for (int32_t ky = (kx == 0) ? 1 : -param_kyMax; ky <= param_kyMax; ++ky)
        {
            for (int32_t kz = -param_kzMax; kz <= +param_kzMax; ++kz)
            {
                double rkx = kx * 2 * M_PI / param_a;
                double rky = ky * 2 * M_PI / param_b;
                double rkz = kz * 2 * M_PI / param_c;
                double rknormsq = rkx * rkx + rky * rky + rkz * rkz;

                if (rknormsq < param_rksqmax)
                {
                    double skCos = 0.0;
                    double skSin = 0.0;
                    for (int32_t i = 0; i < numWall; ++i)
                    {
                        double coskx = instream_cossinSkx[(i + kx * numWall) * 2 + 0];
                        double sinkx = instream_cossinSkx[(i + kx * numWall) * 2 + 1];

                        double cosky = instream_cossinSky[(i + abs(ky) * numWall) * 2 + 0];
                        double sinky = (ky >= 0) ? instream_cossinSky[(i + ky * numWall) * 2 + 1]
                                                 : -instream_cossinSky[(i - ky * numWall) * 2 + 1];

                        double coskz = instream_cossinSkz[(i + abs(kz) * numWall) * 2 + 0];
                        double sinkz = (kz >= 0) ? instream_cossinSkz[(i + kz * numWall) * 2 + 1]
                                                 : -instream_cossinSkz[(i - kz * numWall) * 2 + 1];

                        double coskxky = coskx * cosky - sinkx * sinky;
                        double sinkxky = sinkx * cosky + coskx * sinky;

                        double coskxkykz = coskxky * coskz - sinkxky * sinkz;
                        double sinkxkykz = sinkxky * coskz + coskxky * sinkz;

                        skCos += qWall[i] * coskxkykz;
                        skSin += qWall[i] * sinkxkykz;
                        /*if (kx==0 & ky==1 & kz <-param_kzMax+1 & i<4) printf("
                           lr CPU %d %d %d  %+4.15e %+4.15e
                           %+4.15e %+4.15e  \n", kx,ky,kz,coskxkykz,sinkxkykz,
                           cosky, sinky);//*/

                        // printf("%d %d %12.5E %12.5E %12.5E %12.5E %12.5E
                        // %12.5E
                        // \n",count,i,coskx,cosky,coskz,sinkx,sinky,sinkz);
                        // printf("%d %d %12.5E %12.5E  \n",count,i +
                        // abs(kz)*numWall,coskz,sinkz);
                    }

                    double alphaSk = exp(alphaconst * rknormsq) / rknormsq;

                    // printf("%d %d %d %12.5E %12.5E %12.5E  %12.5E %12.5E
                    // %12.5E %12.5E \n",kx,ky,kz,rkx,rky,rkz,rknormsq,alphaSk,
                    // skCos,skSin); if (kx==0 & ky==1 & kz <-param_kzMax+4)
                    // printf(" lr CPU %+4.15e %+4.15e
                    // %+4.15e  \n", skCos,skSin,alphaSk);

                    for (int32_t i = 0; i < numWall; ++i)
                    {
                        double coskx = instream_cossinSkx[(i + kx * numWall) * 2 + 0];
                        double sinkx = instream_cossinSkx[(i + kx * numWall) * 2 + 1];

                        double cosky = instream_cossinSky[(i + abs(ky) * numWall) * 2 + 0];
                        double sinky = (ky >= 0) ? instream_cossinSky[(i + ky * numWall) * 2 + 1]
                                                 : -instream_cossinSky[(i - ky * numWall) * 2 + 1];

                        double coskz = instream_cossinSkz[(i + abs(kz) * numWall) * 2 + 0];
                        double sinkz = (kz >= 0) ? instream_cossinSkz[(i + kz * numWall) * 2 + 1]
                                                 : -instream_cossinSkz[(i - kz * numWall) * 2 + 1];

                        double coskxky = coskx * cosky - sinkx * sinky;
                        double sinkxky = sinkx * cosky + coskx * sinky;

                        double coskxkykz = coskxky * coskz - sinkxky * sinkz;
                        double sinkxkykz = sinkxky * coskz + coskxky * sinkz;

                        outstream_MWPotCPU[i] += lrPotFactor * alphaSk * (coskxkykz * skCos + sinkxkykz * skSin);

                        // printf("%d %12.5E %12.5E %12.5E  %12.5E %12.5E %12.5E
                        // %12.5E %12.5E
                        // \n",count,coskx,cosky,coskz,sinkx,sinky,sinkz,coskxkykz,sinkxkykz
                        // );

                        // if (i==0) printf("%d %+4.15e %+4.15e  %+4.15e %+4.15e
                        // %+4.15e
                        // %+4.15e %+4.15e
                        // \n",count,alphaSk,coskxkykz,sinkxkykz,skCos ,skSin,
                        // lrPotFactor * alphaSk * (coskxkykz*skCos
                        // + sinkxkykz*skSin), outstream_MWPotCPU[i] );
                        // printf("%d %d
                        // %+4.15e %+4.15e  %+4.15e %+4.15e \n",count,i,skCos
                        // ,skSin, lrPotFactor * alphaSk * (coskxkykz*skCos +
                        // sinkxkykz*skSin), outstream_MWPotCPU[i] );

                        // printf("%d %+4.15e %+4.15e \n",count,lrPotFactor *
                        // alphaSk * (coskxkykz*skCos +
                        // sinkxkykz*skSin),outstream_MWPotCPU[i] );
                    }
                }
                // count +=1;
            }
        }
    }
    // printf(" count mode CPU %d \n",count);
}

/**
 * \brief Basic static function for the interface 'srPotWall1'.
 *
 * \param [in] param_a Interface Parameter "a".: Box length in the 1st
 * dimension \param [in] param_alpha Interface Parameter "alpha".: Ewald
 * gaussian parameter \param [in] param_b Interface Parameter "b".: Box length
 * in the 2nd dimension \param [in] param_eta Interface Parameter "eta".:
 * Electrode gaussian parameter \param [in] param_numWall1 Interface Parameter
 * "numWall1".: Number of particles in the first electrode \param [in]
 * param_rcutsq Interface Parameter "rcutsq".: Short-range cutoff squared
 * \param [in] xyzWall The stream should be of size (4 * param_numWall1 * 8)
 * bytes. \param [in] q The stream should be of size (param_numWall1 * 8)
 * bytes. \param [out] outstream_MWPotCPU The stream should be of size
 * (param_numWall1
 * * 8) bytes.
 */
void MWPot_srPotWallCPU(double param_a, double param_alpha, double param_b, double param_eta, int32_t param_numWall1,
                        double param_rcutsq, double *xyz, double *q, double *outstream_MWPotCPU)
{

    double arec = 1.0 / param_a;
    double brec = 1.0 / param_b;
    double etasqrt2 = param_eta / sqrt(2.0);

    for (int32_t i = 0; i < param_numWall1; ++i)
    {
        outstream_MWPotCPU[i] = 0.0;
        double xi = xyz[i * 4 + 0];
        double yi = xyz[i * 4 + 1];
        double zi = xyz[i * 4 + 2];

        for (int32_t j = 0; j < param_numWall1; ++j)
        {
            double xj = xyz[j * 4 + 0];
            double yj = xyz[j * 4 + 1];
            double zj = xyz[j * 4 + 2];
            double qj = q[j];

            // Compute minimum image distance
            double xij = xj - xi;
            double scalex = floor(xij * arec + 0.5);
            xij -= scalex * param_a;

            double yij = yj - yi;
            double scaley = floor(yij * brec + 0.5);
            yij -= scaley * param_b;

            double zij = zj - zi;

            double drnorm2 = xij * xij + yij * yij + zij * zij;
            if ((i != j) && (drnorm2 < param_rcutsq))
            {
                double drnorm = sqrt(drnorm2);
                double erfcalpha = erfc(param_alpha * drnorm);
                double erfceta = erfc(etasqrt2 * drnorm);

                outstream_MWPotCPU[i] += qj * (erfcalpha - erfceta) / drnorm;
            }
        }
    }
}

void MWPot_k0PotCPU(double param_a, double param_b, double param_alpha, int32_t param_numWall1, int32_t param_numWall2,
                    double *z, double *q, double *outstream_MWPotCPU)
{

    int32_t numWall = param_numWall1 + param_numWall2;
    double alphasq = param_alpha * param_alpha;
    double sqrpialpha = sqrt(M_PI) / param_alpha;
    double k0PotFactor = 2.0 / (param_a * param_b);

    for (int32_t i = 0; i < numWall; ++i)
    {
        outstream_MWPotCPU[i] = 0;
        double zi = z[i];
        for (int32_t j = 0; j < i; ++j)
        {
            double zj = z[j];
            double zij = zj - zi;
            double zijsq = zij * zij;
            double expij = exp(-zijsq * alphasq);
            double erfij = erf(zij * param_alpha);
            double k0potij = sqrpialpha * expij + M_PI * zij * erfij;
            outstream_MWPotCPU[j] -= k0PotFactor * q[i] * k0potij;
            outstream_MWPotCPU[i] -= k0PotFactor * q[j] * k0potij;
        }
    }
    for (int32_t i = 0; i < numWall; ++i)
    {
        outstream_MWPotCPU[i] -= q[i] * k0PotFactor * sqrpialpha;
    }
}

void init_xyz(int num, double *xyz, double *x, double *y, double *z)
{
    for (int32_t i = 0; i < num; ++i)
    {
        xyz[4 * i] = x[i];
        xyz[4 * i + 1] = y[i];
        xyz[4 * i + 2] = z[i];
        xyz[4 * i + 3] = 0.0;
    }
}

void init_instream_cossinSk(int32_t numWall, int32_t kMax, int32_t ixyz, double kdim, double *xyzWall,
                            double *instream_cossinSk)
{

    for (int32_t k = 0; k <= kMax; ++k)
    {
        for (int32_t i = 0; i < numWall; ++i)
        {
            instream_cossinSk[(i + k * numWall) * 2 + 0] = cos(k * xyzWall[4 * i + ixyz] * kdim);
            instream_cossinSk[(i + k * numWall) * 2 + 1] = sin(k * xyzWall[4 * i + ixyz] * kdim);
        }
    }
}