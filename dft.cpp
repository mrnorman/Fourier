
#include <iostream>
#include <iomanip>
#include <cmath>

#define YAKL_INLINE inline


template <unsigned N, class T=double> class DFT {
  T cos_table[N];
  T sin_table[N];
public:
  YAKL_INLINE constexpr DFT() {
    for (unsigned i=0; i<N; i++) {
      cos_table[i] = cos(2*M_PI*i/N);
      sin_table[i] = sin(2*M_PI*i/N);
    }
  }

  static YAKL_INLINE unsigned constexpr wrap( unsigned ind ) {
    return ind - (ind/N)*N;
  }

  YAKL_INLINE void forward_one(T const *data, T *fft, unsigned l) const {
    fft[l] = 0;
    // The i-loop is sequential
    for (unsigned i=0; i<N; i++) {
      unsigned ind = wrap(i*(l/2));
      T trig;
      if (l - ( (l >> 1) << 1 ) == 0) { // This is a cheap l%2 operator
        trig = cos_table[ind];
      } else {
        trig = -sin_table[ind];
      }
      fft[l] += trig*data[i];
    }
    fft[l] /= N;
  }

  YAKL_INLINE void inverse_one(T *data, T const *fft, unsigned l) const {
    data[l] = 0;
    // The i-loop is sequential
    for (unsigned i=0; i<N; i++) {
      unsigned ind_tab = wrap(i*l);
      unsigned ind_fft = 2*i;
      int      sgn     = -1;
      if (i > N/2) {
        ind_fft = 2*(N-i);
        sgn     = 1;
      }
      data[l] +=     fft[ind_fft  ]*cos_table[ind_tab] +
                 sgn*fft[ind_fft+1]*sin_table[ind_tab];
    }
  }


  YAKL_INLINE void forward(T const *data, T *fft) const {
    for (unsigned l=0; l<N+2; l++) {
      forward_one(data, fft, l);
    }
  }


  YAKL_INLINE void inverse(T *data, T const *fft) const {
    for (unsigned l=0; l<N; l++) {
      inverse_one(data, fft, l);
    }
  }
};






int main() {
  {
    unsigned constexpr N = 12;
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
    for (unsigned i=0; i<N; i++) {
      data[i] = pow(i+14.2,2);
    }

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



