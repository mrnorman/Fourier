
#include <iostream>
#include <iomanip>
#include <cmath>


template <unsigned N, class T=double> class DFT {
public:
  T cos_table[N][N];
  T sin_table[N][N];


  constexpr DFT() {
    // InitTables<N-1,N-1,N,T>::init(cos_table,sin_table);
    for (unsigned j=0; j<N; j++) {
      for (unsigned i=0; i<N; i++) {
        cos_table[i][j] = cos(2*M_PI*i*j/N);
        sin_table[i][j] = sin(2*M_PI*i*j/N);
      }
    }
  }


  inline void forward_one(T const *data, T *fft, unsigned l) const {
    fft[l] = 0;
    // The i-loop is sequential
    for (unsigned i=0; i<N; i++) {
      if (l - ( (l >> 1) << 1 ) == 0) { // This is a cheap l%2 operator
        fft[l] +=  data[i]*cos_table[i][l/2];
      } else {
        fft[l] += -data[i]*sin_table[i][l/2];
      }
    }
    fft[l] /= N;
  }


  inline void inverse_one(T *data, T const *fft, unsigned l) const {
    data[l] = 0;
    // The i-loop is sequential
    for (unsigned i=0; i<N; i++) {
      if (i <= N/2) {
        data[l] += fft[2*i  ]*cos_table[i][l] -
                   fft[2*i+1]*sin_table[i][l];
      } else {
        data[l] += fft[2*(N-1-i+1)  ]*cos_table[i][l] +
                   fft[2*(N-1-i+1)+1]*sin_table[i][l];
      }
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





unsigned constexpr N = 8;

int main() {
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



