
#include <iostream>
#include <iomanip>
#include <cmath>


int constexpr N = 8;

int main() {
  double data[N];
  double out [N+2];


  for (int i=0; i<N; i++) {
    data[i] = pow(i+14.2,2);
  }
  for (int i=0; i<N; i++) {
    std::cout << std::setprecision(15) << data[i]<< "\n";
  }
  std::cout << std::endl;
  std::cout << std::endl;


  // The l-loop is parallel
  for (int l=0; l<N+2; l++) {
    out[l] = 0;
    // The i-loop is sequential
    for (int i=0; i<N; i++) {
      if (l - ( (l >> 1) << 1 ) == 0) { // This is a cheap l%2 operator
        out[l] +=  data[i]*cos(2*M_PI*i*(l/2)/N);
      } else {
        out[l] += -data[i]*sin(2*M_PI*i*(l/2)/N);
      }
    }
    out[l] /= N;
  }


  for (int i=0; i<N+2; i+=2) {
    std::cout << std::fixed << std::setprecision(15) << out[i] << " + " << out[i+1] << "i\n";
  }
  std::cout << std::endl;
  std::cout << std::endl;


  // The l-loop is parallel
  for (int l=0; l<N; l++) {
    data[l] = 0;
    // The i-loop is sequential
    for (int i=0; i<N; i++) {
      if (i <= N/2) {
        data[l] += out[2*i  ]*cos(2*M_PI*i*l/N) -
                   out[2*i+1]*sin(2*M_PI*i*l/N);
      } else {
        data[l] += out[2*(N-1-i+1)  ]*cos(2*M_PI*i*l/N) +
                   out[2*(N-1-i+1)+1]*sin(2*M_PI*i*l/N);
      }
    }
  }


  for (int i=0; i<N; i++) {
    std::cout << std::setprecision(15) << data[i]<< "\n";
  }
  std::cout << std::endl;
  std::cout << std::endl;
}



