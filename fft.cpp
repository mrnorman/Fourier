
#include <iostream>
#include <iomanip>
#include <cmath>


template<unsigned N, typename T=double> class DanielsonLanczos {
public:
  static inline void constexpr apply(T* data) {
    DanielsonLanczos<N/2,T>::apply(data  );
    DanielsonLanczos<N/2,T>::apply(data+N);
 
    T wtemp,tempr,tempi,wr,wi,wpr,wpi;
    wtemp = sin(M_PI/N);
    wpr = -2.0*wtemp*wtemp;
    wpi = -sin(2*M_PI/N);
    wr = 1.0;
    wi = 0.0;
    for (unsigned i=0; i<N; i+=2) {
      tempr = data[i+N]*wr - data[i+N+1]*wi;
      tempi = data[i+N]*wi + data[i+N+1]*wr;
      data[i+N  ] = data[i  ]-tempr;
      data[i+N+1] = data[i+1]-tempi;
      data[i  ] += tempr;
      data[i+1] += tempi;
 
      wtemp = wr;
      wr += wr*wpr - wi   *wpi;
      wi += wi*wpr + wtemp*wpi;
    }
  }
};
template<typename T> class DanielsonLanczos<1,T> {
public:
  static inline void constexpr apply(T* data) { }
};



template <class T> void swap(T &a, T &b) {
  T tmp = a;
  a = b;
  b = tmp;
}



template <class T> void scramble(T *data , int N ) {
  int n = N<<1;
  int j=1;
  for (int i=1; i<n; i+=2) {
    if (j>i) {
      swap(data[j-1], data[i-1]);
      swap(data[j  ], data[i  ]);
    }
    int m = N;
    while (m>=2 && j>m) {
      j -= m;
      m >>= 1;
    }
    j += m;
  };
}



template <int I, int N, class T> class ProcessRealFFT {
public:
  static inline void constexpr process(T *data, T *tmp) {
    T xrp = (tmp[2*I  ] + tmp[2*(N-I)  ])*0.5;
    T xrm = (tmp[2*I  ] - tmp[2*(N-I)  ])*0.5;
    T xip = (tmp[2*I+1] + tmp[2*(N-I)+1])*0.5;
    T xim = (tmp[2*I+1] - tmp[2*(N-I)+1])*0.5;
    data[2*I  ] = ( xrp + cos(M_PI*I/N)*xip - sin(M_PI*I/N)*xrm )/(2*N);
    data[2*I+1] = ( xim - sin(M_PI*I/N)*xip - cos(M_PI*I/N)*xrm )/(2*N);
    ProcessRealFFT<I-1,N,T>::process(data,tmp);
  }
};
template <int N, class T> class ProcessRealFFT<1,N,T> {
public:
  static inline void constexpr process(T *data, T *tmp) {
    int constexpr I = 1;
    T xrp = (tmp[2*I  ] + tmp[2*(N-I)  ])*0.5;
    T xrm = (tmp[2*I  ] - tmp[2*(N-I)  ])*0.5;
    T xip = (tmp[2*I+1] + tmp[2*(N-I)+1])*0.5;
    T xim = (tmp[2*I+1] - tmp[2*(N-I)+1])*0.5;
    data[2*I  ] = ( xrp + cos(M_PI*I/N)*xip - sin(M_PI*I/N)*xrm )/(2*N);
    data[2*I+1] = ( xim - sin(M_PI*I/N)*xip - cos(M_PI*I/N)*xrm )/(2*N);
  }
};



template <int I, int N, class T> class ProcessRealInverseFFT {
public:
  static inline void constexpr process(T *data, T *tmp) {
    T xrp = (data[2*I  ] + data[2*(N-I)  ]);
    T xrm = (data[2*I  ] - data[2*(N-I)  ]);
    T xip = (data[2*I+1] + data[2*(N-I)+1]);
    T xim = (data[2*I+1] - data[2*(N-I)+1]);
    tmp[2*I  ] = xrp - cos(M_PI*I/N)*xip - sin(M_PI*I/N)*xrm;
    tmp[2*I+1] = xim - sin(M_PI*I/N)*xip + cos(M_PI*I/N)*xrm;
    ProcessRealInverseFFT<I-1,N,T>::process(data,tmp);
  }
};
template <int N, class T> class ProcessRealInverseFFT<0,N,T> {
public:
  static inline void constexpr process(T *data, T *tmp) {
    int constexpr I = 0;
    T xrp = (data[2*I  ] + data[2*(N-I)  ]);
    T xrm = (data[2*I  ] - data[2*(N-I)  ]);
    T xip = (data[2*I+1] + data[2*(N-I)+1]);
    T xim = (data[2*I+1] - data[2*(N-I)+1]);
    tmp[2*I  ] = xrp - cos(M_PI*I/N)*xip - sin(M_PI*I/N)*xrm;
    tmp[2*I+1] = xim - sin(M_PI*I/N)*xip + cos(M_PI*I/N)*xrm;
  }
};



// Calculated at compile time
constexpr unsigned nextPowerOfTwo(unsigned n) {
  unsigned count = 0;  
  // If n is zero or n is a power of 2, then return it
  if (n && !(n & (n - 1))) { return n; }
  while( n != 0) {
    n >>= 1;  
    count += 1;  
  }
  return 1 << count; 
}



template<unsigned SIZE, typename T=double> class GFFT {
  static unsigned constexpr N = nextPowerOfTwo(SIZE)/2;
public:
  void forward(T* data) {
    scramble(data,N);
    DanielsonLanczos<N,T>::apply(data);
  }
  void inverse(T* data) {
    // Multiply complex components by -1
    for (int i=0; i<2*N; i+=2) { data[i+1] = -data[i+1]; }
    forward(data);
    // Multiply complex components by -1
    for (int i=0; i<2*N; i+=2) { data[i+1] = -data[i+1]; }
  }
  void forwardReal(T *data, T *tmp) {
    // Copy to temporary buffer
    for (int i=0; i<2*N; i++) {
      tmp[i] = data[i];
    }
    // Compute FFT assuming complex #s are even,I*odd; even,I*odd
    forward(tmp);
    data[0    ] = (tmp[0] + tmp[1])/(2*N);
    data[1    ] = 0;
    data[2*N  ] = (tmp[0] - tmp[1])/(2*N);
    data[2*N+1] = 0;
    // Transform the FFT into the true FFT for the real sequence
    ProcessRealFFT<N-1,N,T>::process(data,tmp);
  }
  void inverseReal(T* data, T *tmp) {
    // Transform FFTs into something whose inverse reproduces the original real signal
    ProcessRealInverseFFT<N-1,N,T>::process(data,tmp);
    inverse(tmp);
    for (int i=0; i<2*N; i++) {
      data[i] = tmp[i];
    }
  }
};



int main() {
  int constexpr N = 8;
  GFFT<N> gfft;
  double data[N+2];
  double tmp [N];


  for (int i=0; i<N; i++) {
    data[i] = pow(i+14.2,2);
  }
  for (int i=0; i<N; i++) {
    std::cout << std::setprecision(15) << data[i]<< "\n";
  }
  std::cout << std::endl;
  std::cout << std::endl;


  // Forward FFT
  gfft.forwardReal(data,tmp);
  for (int i=0; i<N+2; i+=2) {
    std::cout << data[i] << " + " << data[i+1] << "i\n";
  }
  std::cout << std::endl;
  std::cout << std::endl;


  gfft.inverseReal(data,tmp);
  for (int i=0; i<N; i++) {
    std::cout << std::setprecision(15) << data[i] << "\n";
  }
  std::cout << std::endl;
  std::cout << std::endl;


}



