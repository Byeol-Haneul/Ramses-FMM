#include <stdint.h>

void hilbert3d_c_(uint32_t *x,uint32_t *y,uint32_t *z,double *order,uint32_t *bit_length,int32_t *npoint) {
  uint64_t s = 0;
  uint32_t m,ux,uy,uz,ut;
  int i, n = *npoint;

  for( i=0; i<n; ++i ) {
    ux = x[i];
    uy = y[i];
    uz = z[i];

    m = 1 << *bit_length;
    // e.g., m = 0x00100000;
    while (m) {
      s = s << 3;

      if (ux&m) {
	if (uy&m) {
	  if (uz&m) {
	    ut = ux;
	    ux = uy;
	    uy = ~uz;
	    uz = ~ut;
	    s |= 5;
	  }
	  else {
	    ut = uz;
	    uz = ux;
	    ux = uy;
	    uy = ut;
	    s |= 2;
	  }
	}
	else {
	  ux = ~ux;
	  uy = ~uy;
	  if (uz&m) {
	    s |= 4;
	  }
	  else {
	    s |= 3;
	  }
	}
      }
      else {
	if (uy&m) {
	  if (uz&m) {
	    ut = ux;
	    ux = uy;
	    uy = ~uz;
	    uz = ~ut;
	    s |= 6;
	  }
	  else {
	    ut = uz;
	    uz = ux;
	    ux = uy;
	    uy = ut;
	    s |= 1;
	  }
	}
	else {
	  if (uz&m) {
	    ut = uy;
	    uy = ux;
	    ux = ~uz;
	    uz = ~ut;
	    s |= 7;
	  }
	  else {
	    ut = uy;
	    uy = ux;
	    ux = uz;
	    uz = ut;
	    s |= 0;
	  }
	}
      }
      m = m >> 1;
    }
    order[i] = s;
  }
}
