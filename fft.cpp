
#include <iostream>
#include <iomanip>
#include <cmath>


template<unsigned M, unsigned N, unsigned B, unsigned A>
struct SinCosSeries {
   static double value() {
      return 1-(A*M_PI/B)*(A*M_PI/B)/M/(M+1)
               *SinCosSeries<M+2,N,B,A>::value();
   }
};
 
template<unsigned N, unsigned B, unsigned A>
struct SinCosSeries<N,N,B,A> {
   static double value() { return 1.; }
};
 
template<unsigned B, unsigned A, typename T=double>
struct Sin;
 
template<unsigned B, unsigned A>
struct Sin<B,A,float> {
   static float value() {
      return (A*M_PI/B)*SinCosSeries<2,24,B,A>::value();
   }
};
template<unsigned B, unsigned A>
struct Sin<B,A,double> {
   static double value() {
      return (A*M_PI/B)*SinCosSeries<2,34,B,A>::value();
   }
};
 
template<unsigned B, unsigned A, typename T=double>
struct Cos;
 
template<unsigned B, unsigned A>
struct Cos<B,A,float> {
   static float value() {
      return SinCosSeries<1,23,B,A>::value();
   }
};
template<unsigned B, unsigned A>
struct Cos<B,A,double> {
   static double value() {
      return SinCosSeries<1,33,B,A>::value();
   }
};




template<unsigned N, typename T=double> class DanielsonLanczos {
public:
  static inline void constexpr apply(T* data) {
    DanielsonLanczos<N/2,T>::apply(data  );
    DanielsonLanczos<N/2,T>::apply(data+N);
 
    T wtemp,tempr,tempi,wr,wi,wpr,wpi;
    wtemp = Sin<N,1,T>::value(); // sin(M_PI/N);
    wpr = -2.0*wtemp*wtemp;
    wpi = -Sin<N,2,T>::value(); // -sin(2*M_PI/N);
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



template <class T> inline constexpr void swap(T &a, T &b) {
  T tmp = a;
  a = b;
  b = tmp;
}



template <class T> void scramble(T *data , unsigned N ) {
  unsigned n = N<<1;
  unsigned j=1;
  for (unsigned i=1; i<n; i+=2) {
    if (j>i) {
      swap(data[j-1], data[i-1]);
      swap(data[j  ], data[i  ]);
    }
    unsigned m = N;
    while (m>=2 && j>m) {
      j -= m;
      m >>= 1;
    }
    j += m;
  };
}



template <unsigned I, unsigned N, class T> class ProcessRealFFT {
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
template <unsigned N, class T> class ProcessRealFFT<1,N,T> {
public:
  static inline void constexpr process(T *data, T *tmp) {
    unsigned constexpr I = 1;
    T xrp = (tmp[2*I  ] + tmp[2*(N-I)  ])*0.5;
    T xrm = (tmp[2*I  ] - tmp[2*(N-I)  ])*0.5;
    T xip = (tmp[2*I+1] + tmp[2*(N-I)+1])*0.5;
    T xim = (tmp[2*I+1] - tmp[2*(N-I)+1])*0.5;
    data[2*I  ] = ( xrp + cos(M_PI*I/N)*xip - sin(M_PI*I/N)*xrm )/(2*N);
    data[2*I+1] = ( xim - sin(M_PI*I/N)*xip - cos(M_PI*I/N)*xrm )/(2*N);
  }
};



template <unsigned I, unsigned N, class T> class ProcessRealInverseFFT {
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
template <unsigned N, class T> class ProcessRealInverseFFT<0,N,T> {
public:
  static inline void constexpr process(T *data, T *tmp) {
    unsigned constexpr I = 0;
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
  static_assert(SIZE-N*2 == 0,"ERROR: Running GFFT with a non-power-of-two-size");
public:
  void forward(T* data) {
    scramble(data,N);
    DanielsonLanczos<N,T>::apply(data);
  }
  void inverse(T* data) {
    // Multiply complex components by -1
    for (unsigned i=0; i<2*N; i+=2) { data[i+1] = -data[i+1]; }
    forward(data);
    // Multiply complex components by -1
    for (unsigned i=0; i<2*N; i+=2) { data[i+1] = -data[i+1]; }
  }
  void forwardReal(T *data, T *tmp) {
    // Copy to temporary buffer
    for (unsigned i=0; i<2*N; i++) {
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
    for (unsigned i=0; i<2*N; i++) {
      data[i] = tmp[i];
    }
  }
};



int main() {
  unsigned constexpr N = 8;
  GFFT<N> gfft;
  double data[N+2];
  double tmp [N];


  for (unsigned i=0; i<N; i++) {
    data[i] = pow(i+14.2,2);
  }
  for (unsigned i=0; i<N; i++) {
    std::cout << std::setprecision(15) << data[i]<< "\n";
  }
  std::cout << std::endl;
  std::cout << std::endl;


  // Forward FFT
  gfft.forwardReal(data,tmp);
  for (unsigned i=0; i<N+2; i+=2) {
    std::cout << data[i] << " + " << data[i+1] << "i\n";
  }
  std::cout << std::endl;
  std::cout << std::endl;


  gfft.inverseReal(data,tmp);
  for (unsigned i=0; i<N; i++) {
    std::cout << std::setprecision(15) << data[i] << "\n";
  }
  std::cout << std::endl;
  std::cout << std::endl;
}



