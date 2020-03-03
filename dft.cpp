
#include <iostream>
#include <iomanip>
#include <cmath>


template <unsigned N, class T=double> class DFT {
public:
  T sin_table[N*N];


  constexpr DFT() {
    for (unsigned i=0; i<N*N; i++) {
      sin_table[i] = sin(2*M_PI*i/N);
    }
  }


  inline void forward_one(T const *data, T *fft, unsigned l) const {
    fft[l] = 0;
    // The i-loop is sequential
    for (unsigned i=0; i<N; i++) {
      int sgn = -1;
      unsigned ind = i*(l/2);  // Defaults to sin(2*pi*i*(l/2)/N)
      // Only do integer arithmetic in the branched section
      if (l - ( (l >> 1) << 1 ) == 0) { // This is a cheap l%2 operator
        // if l is even, compute cos(2*pi*i*(l/2)/N) by shifting by pi/2
        ind = i*(l/2) + N/4;
        if (ind >= N*N) { ind -= N; }
        sgn = 1;
      }
      fft[l] += sgn*data[i]*sin_table[ind];
    }
    fft[l] /= N;
  }


  inline void inverse_one(T *data, T const *fft, unsigned l) const {
    data[l] = 0;
    // The i-loop is sequential
    for (unsigned i=0; i<N; i++) {
      T mysin = sin_table[i*l];
      unsigned ind = i*l + N/4;
      if (ind >= N*N) { ind -= N; }
      T mycos = sin_table[ind];
      // Only do integer arithmetic in the branched section
      if (i <= N/2) {
        ind = 2*i;
        mysin *= -1;
      } else {
        ind = 2*(N-i);
      }
      data[l] += fft[ind]*mycos + fft[ind+1]*mysin;
    }
  }


  inline void forward(T const *data, T *fft) const {
    for (unsigned l=0; l<N+2; l++) {
      forward_one(data, fft, l);
    }
  }


  inline void inverse(T *data, T const *fft) const {
    for (unsigned l=0; l<N; l++) {
      inverse_one(data, fft, l);
    }
  }
};






int main() {
  {
    unsigned constexpr N = 8;
    double data[N  ];
    double fft [N+2];

    DFT<N> dft;


    for (unsigned i=0; i<N; i++) {
      data[i] = pow(i+14.2,2);
    }
    for (unsigned i=0; i<N; i++) {
      std::cout << std::setprecision(15) << data[i]<< "\n";
    }
    std::cout << std::endl;
    std::cout << std::endl;


    // The l-loop is parallel
      for (unsigned l=0; l<N+2; l++) {
        dft.forward_one(data,fft,l);
      }


    for (unsigned i=0; i<N+2; i+=2) {
      std::cout << std::fixed << std::setprecision(15) << fft[i] << " + " << fft[i+1] << "i\n";
    }
    std::cout << std::endl;
    std::cout << std::endl;


    // The l-loop is parallel
      for (unsigned l=0; l<N; l++) {
        dft.inverse_one(data,fft,l);
      }


    for (unsigned i=0; i<N; i++) {
      std::cout << std::setprecision(15) << data[i]<< "\n";
    }
    std::cout << std::endl;
    std::cout << std::endl;
  }


  {
    unsigned constexpr N = 64;
    unsigned constexpr ITER = 100000;
    DFT<N> dft;
    double data[N  ];
    double ffts[N+2];

    auto t1 = std::clock();

    for (unsigned i=0; i<ITER; i++) {
      dft.forward(data,ffts);
      dft.inverse(data,ffts);
    }

    auto tm = std::clock() - t1;
    std::cout << "Cycles: " << tm << "\n";
    std::cout << data[0] << "\n";
  }
}



