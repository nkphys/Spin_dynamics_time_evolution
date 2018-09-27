#pragma once

#include <complex>
#include <vector>
#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <stdexcept>
#include "functions.h"


class Fft {


public:
	double PI_EFF ;
	
	/* 
	 * Computes the discrete Fourier transform (DFT) of the given complex vector, storing the result back into the vector.
	 * The vector can have any length. This is a wrapper function.
	 */
	void transform(std::vector<std::complex<double> > &vec);
	
	
	/* 
	 * Computes the inverse discrete Fourier transform (IDFT) of the given complex vector, storing the result back into the vector.
	 * The vector can have any length. This is a wrapper function. This transform does not perform scaling, so the inverse is not a true inverse.
	 */
	void inverseTransform(std::vector<std::complex<double> > &vec);
	
	
	/* 
	 * Computes the discrete Fourier transform (DFT) of the given complex vector, storing the result back into the vector.
	 * The vector's length must be a power of 2. Uses the Cooley-Tukey decimation-in-time radix-2 algorithm.
	 */
	void transformRadix2(std::vector<std::complex<double> > &vec);
	
	
	/* 
	 * Computes the discrete Fourier transform (DFT) of the given complex vector, storing the result back into the vector.
	 * The vector can have any length. This requires the convolution function, which in turn requires the radix-2 FFT function.
	 * Uses Bluestein's chirp z-transform algorithm.
	 */
	void transformBluestein(std::vector<std::complex<double> > &vec);
	
	
	/* 
	 * Computes the circular convolution of the given complex vectors. Each vector's length must be the same.
	 */
	void convolve(
		const std::vector<std::complex<double> > &vecx,
		const std::vector<std::complex<double> > &vecy,
		std::vector<std::complex<double> > &vecout);
	
};



void Fft::transform(vector<complex<double> > &vec) {
   int n = vec.size();
    if (n == 0)
        return;
    else if ((n & (n - 1)) == 0)  // Is power of 2
        transformRadix2(vec);
    else  // More complicated algorithm for arbitrary sizes
        transformBluestein(vec);
}


void Fft::inverseTransform(vector<complex<double> > &vec) {
    std::transform(vec.cbegin(), vec.cend(), vec.begin(),
        static_cast<complex<double> (*)(const complex<double> &)>(std::conj));
    transform(vec);
    std::transform(vec.cbegin(), vec.cend(), vec.begin(),
        static_cast<complex<double> (*)(const complex<double> &)>(std::conj));
}


void Fft::transformRadix2(vector<complex<double> > &vec) {


    // Length variables
   int n = vec.size();
    int levels = 0;  // Compute levels = floor(log2(n))
    for (int temp = n; temp > 1U; temp >>= 1)
        levels++;
    if (static_cast<int>(1U) << levels != n)
        throw std::domain_error("Length is not a power of 2");

    // Trignometric table
    vector<complex<double> > expTable(n / 2);
    for (int i = 0; i < n / 2; i++)
        expTable[i] = std::exp(complex<double>(0, -2 * PI_EFF * i / n));

    // Bit-reversed addressing permutation
    for (int i = 0; i < n; i++) {
       int j = reverseBits(i, levels);
        if (j > i)
            std::swap(vec[i], vec[j]);
    }

    // Cooley-Tukey decimation-in-time radix-2 FFT
    for (int size = 2; size <= n; size *= 2) {
       int halfsize = size / 2;
       int tablestep = n / size;
        for (int i = 0; i < n; i += size) {
            for (int j = i, k = 0; j < i + halfsize; j++, k += tablestep) {
                complex<double> temp = vec[j + halfsize] * expTable[k];
                vec[j + halfsize] = vec[j] - temp;
                vec[j] += temp;
            }
        }
        if (size == n)  // Prevent overflow in 'size *= 2'
            break;
    }
}


void Fft::transformBluestein(vector<complex<double> > &vec) {
    // Find a power-of-2 convolution length m such that m >= n * 2 + 1

    //double w_min = ;
    //double dw = ;
    //double dt = ;
   int n = vec.size();
   int m = 1;
    while (m / 2 <= n) {
        if (m > SIZE_MAX / 2)
            throw std::length_error("Vector too large");
        m *= 2;
    }

    // Trignometric table
    vector<complex<double> > expTable(n);
    for (int i = 0; i < n; i++) {
        unsigned long long temp = static_cast<unsigned long long>(i) * i;
        temp %= static_cast<unsigned long long>(n) * 2;
        double angle = PI_EFF * temp / n;
        // Less accurate alternative if long long is unavailable: double angle = PI_EFF * i * i / n;
        expTable[i] = std::exp(complex<double>(0, -angle));
    }

    // Temporary vectors and preprocessing
    vector<complex<double> > av(m);
    for (int i = 0; i < n; i++) {
        av[i] = vec[i] * expTable[i];
    }
    vector<complex<double> > bv(m);
    bv[0] = expTable[0];
    for (int i = 1; i < n; i++){
        bv[i] = bv[m - i] = std::conj(expTable[i]);
    }

    // Convolution
    vector<complex<double> > cv(m);
    convolve(av, bv, cv);

    // Postprocessing
    for (int i = 0; i < n; i++){
        vec[i] = cv[i] * expTable[i];
        }
}


void Fft::convolve(
        const vector<complex<double> > &xvec,
        const vector<complex<double> > &yvec,
        vector<complex<double> > &outvec) {

   int n = xvec.size();
    if (n != yvec.size() || n != outvec.size()){
        throw std::domain_error("Mismatched lengths");
    }
    vector<complex<double> > xv = xvec;
    vector<complex<double> > yv = yvec;
    transform(xv);
    transform(yv);
    for (int i = 0; i < n; i++){
        xv[i] *= yv[i];
    }
    inverseTransform(xv);
    for (int i = 0; i < n; i++)  // Scaling (because this FFT implementation omits it)
        outvec[i] = xv[i] / static_cast<double>(n);
}
