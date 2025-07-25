#include <iostream>
#include <fftw3.h>

int main() {
    // Define the function f: F_2^n -> R
    int n = 8;  // Size of the function
    

    double f[n] = {1, 0, 1, 0, 0, 1, 1, 0};
    // FFTW variables
    fftw_complex *in, *out;
    fftw_plan plan;

    // Allocate memory
    in = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * n);
    out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * n);

    // Initialize FFTW plan
    plan = fftw_plan_dft_1d(2, in, out, FFTW_FORWARD, FFTW_ESTIMATE);

    // Prepare input for FFTW (convert real numbers to complex)
    for (int i = 0; i < n; ++i) {
        in[i][0] = f[i];  // Real part
        in[i][1] = 0.0;   // Imaginary part
    }

    // Execute FFTW plan
    fftw_execute(plan);

    // Print DFT result
    std::cout << "DFT result:" << std::endl;
    for (int i = 0; i < n; ++i) {
        std::cout << "(" << out[i][0] << ", " << out[i][1] << ")" << std::endl;
    }

    // Clean up
    fftw_destroy_plan(plan);
    fftw_free(in);
    fftw_free(out);
    //delete[] f;

    return 0;
}
