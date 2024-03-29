/* -------------------------------------------------------------------------- *
 *                             OpenMM Laboratory                              *
 *                             =================                              *
 *                                                                            *
 * A plugin for testing low-level code implementation for OpenMM.             *
 *                                                                            *
 * Copyright (c) 2023 Charlles Abreu                                          *
 * https://github.com/craabreu/openmm-lab                                     *
 * -------------------------------------------------------------------------- */

/**
 * This tests the CUDA implementation of FFT3D.
 */

#include "internal/CudaCuFFT3D.h"
#include "internal/CudaVkFFT3D.h"
#include "openmm/internal/AssertionUtilities.h"
#include "openmm/cuda/CudaArray.h"
#include "openmm/cuda/CudaContext.h"
#include "openmm/cuda/CudaSort.h"
#include "sfmt/SFMT.h"
#include "openmm/System.h"
#include <cmath>
#include <complex>
#include <set>
#ifdef _MSC_VER
  #define POCKETFFT_NO_VECTORS
#endif
#include "internal/pocketfft_hdronly.h"

using namespace OpenMMLab;
using namespace OpenMM;
using namespace std;

static CudaPlatform platform;

template <class FFT3D, typename Real, class Real2>
void testTransform(bool realToComplex, int xsize, int ysize, int zsize, int batch) {
    System system;
    system.addParticle(0.0);

    CudaPlatform::PlatformData platformData(
        NULL,
        system,
        "",
        "true",
        platform.getPropertyDefaultValue("CudaPrecision"),
        "false",
        platform.getPropertyDefaultValue(CudaPlatform::CudaTempDirectory()),
        platform.getPropertyDefaultValue(CudaPlatform::CudaDisablePmeStream()),
        "false",
        1,
        NULL
    );
    CudaContext& context = *platformData.contexts[0];
    context.initialize();
    context.setAsCurrent();
    OpenMM_SFMT::SFMT sfmt;
    init_gen_rand(0, sfmt);
    int gridSize = xsize*ysize*zsize;
    int outputZSize = (realToComplex ? zsize/2+1 : zsize);

    vector<vector<complex<double>>> reference(batch);
    for (int j = 0; j < batch; j++) {
        reference[j].resize(gridSize);
        for (int i = 0; i < gridSize; i++) {
            Real x = (float) genrand_real2(sfmt);
            Real y = realToComplex ? 0 : (float) genrand_real2(sfmt);
            reference[j][i] = complex<double>(x, y);
        }
    }

    vector<Real2> complexOriginal(gridSize*batch);
    Real* realOriginal = (Real*) &complexOriginal[0];
    for (int j = 0; j < batch; j++)
        for (int i = 0; i < gridSize; i++) {
            int offset = j*gridSize;
            if (realToComplex)
                realOriginal[offset+i] = reference[j][i].real();
            else {
                complexOriginal[offset+i].x = reference[j][i].real();
                complexOriginal[offset+i].y = reference[j][i].imag();
            }
        }

    CudaArray grid1(context, complexOriginal.size(), sizeof(Real2), "grid1");
    CudaArray grid2(context, complexOriginal.size(), sizeof(Real2), "grid2");
    grid1.upload(complexOriginal);

    CUstream stream = context.getCurrentStream();
    FFT3D fft(context, stream, xsize, ysize, zsize, batch, realToComplex, grid1, grid2);

    // Perform a forward FFT, then verify the result is correct.

    fft.execFFT(true);
    vector<Real2> result;
    grid2.download(result);

    vector<size_t> shape = {(size_t) xsize, (size_t) ysize, (size_t) zsize};
    vector<size_t> axes = {0, 1, 2};
    vector<ptrdiff_t> stride = {(ptrdiff_t) (ysize*zsize*sizeof(complex<double>)),
                                (ptrdiff_t) (zsize*sizeof(complex<double>)),
                                (ptrdiff_t) sizeof(complex<double>)};
    for (int j = 0; j < batch; j++) {
        pocketfft::c2c(shape, stride, stride, axes, true, reference[j].data(), reference[j].data(), 1.0);
        for (int x = 0; x < xsize; x++)
            for (int y = 0; y < ysize; y++)
                for (int z = 0; z < outputZSize; z++) {
                    int index1 = x*ysize*zsize + y*zsize + z;
                    int index2 = ((j*xsize + x)*ysize + y)*outputZSize + z;
                    ASSERT_EQUAL_TOL(reference[j][index1].real(), result[index2].x, 1e-3);
                    ASSERT_EQUAL_TOL(reference[j][index1].imag(), result[index2].y, 1e-3);
                }
    }

    // Perform a backward transform and see if we get the original values.

    fft.execFFT(false);
    grid1.download(result);
    double scale = 1.0/(xsize*ysize*zsize);
    Real* realResult = (Real*) &result[0];
    for (int j = 0; j < batch; j++)
        for (int i = 0; i < gridSize; i++) {
            int offset = j*gridSize;
            if (realToComplex) {
                ASSERT_EQUAL_TOL(realOriginal[offset+i], scale*realResult[offset+i], 1e-4);
            }
            else {
                ASSERT_EQUAL_TOL(complexOriginal[offset+i].x, scale*result[offset+i].x, 1e-4);
                ASSERT_EQUAL_TOL(complexOriginal[offset+i].y, scale*result[offset+i].y, 1e-4);
            }

        }
}

template <class FFT3D, typename Real, class Real2>
void executeTests(int batch) {
    testTransform<FFT3D, Real, Real2>(false, 28, 25, 30, batch);
    testTransform<FFT3D, Real, Real2>(true, 28, 25, 25, batch);
    testTransform<FFT3D, Real, Real2>(true, 25, 28, 25, batch);
    testTransform<FFT3D, Real, Real2>(true, 25, 25, 28, batch);
    testTransform<FFT3D, Real, Real2>(true, 21, 25, 27, batch);
}

int main(int argc, char* argv[]) {
    try {
        if (argc > 1)
            platform.setPropertyDefaultValue("CudaPrecision", string(argv[1]));
        if (platform.getPropertyDefaultValue("CudaPrecision") == "double") {
            executeTests<CudaCuFFT3D, double, double2>(1);
            executeTests<CudaCuFFT3D, double, double2>(2);
            executeTests<CudaCuFFT3D, double, double2>(3);
            executeTests<CudaVkFFT3D, double, double2>(1);
            executeTests<CudaVkFFT3D, double, double2>(2);
            executeTests<CudaVkFFT3D, double, double2>(3);
        }
        else {
            executeTests<CudaCuFFT3D, float, float2>(1);
            executeTests<CudaCuFFT3D, float, float2>(2);
            executeTests<CudaCuFFT3D, float, float2>(3);
            executeTests<CudaVkFFT3D, float, float2>(1);
            executeTests<CudaVkFFT3D, float, float2>(2);
            executeTests<CudaVkFFT3D, float, float2>(3);
        }
    }
    catch(const exception& e) {
        cout << "exception: " << e.what() << endl;
        return 1;
    }
    cout << "Done" << endl;
    return 0;
}