#include <time.h>
#include "simpleRayTracer.h"

__host__ __device__ vector_t vectorCreate(dfloat x, dfloat y, dfloat z){
  vector_t v;
  v.x = x;
  v.y = y;
  v.z = z;
  return v;
}

/* Subtract two vectors and return the resulting vector_t */
__host__ __device__ vector_t vectorSub(const vector_t v1, const vector_t v2){
  return vectorCreate(v1.x - v2.x, v1.y - v2.y, v1.z - v2.z);
}

/* Multiply two vectors and return the resulting scalar (dot product) */
__host__ __device__ dfloat vectorDot(const vector_t v1, const vector_t v2){
  return v1.x * v2.x + v1.y * v2.y + v1.z * v2.z;
}

__host__ __device__ vector_t vectorCrossProduct(const vector_t v1, const vector_t v2){
  return vectorCreate(v1.y*v2.z-v1.z*v2.y,
		      v1.z*v2.x-v1.x*v2.z,
		      v1.x*v2.y-v1.y*v2.x);
}


/* Calculate Vector_T x Scalar and return resulting Vector*/ 
__host__ __device__ vector_t vectorScale(const dfloat c, const vector_t v){
  return vectorCreate(v.x * c, v.y * c, v.z * c);
}

/* Add two vectors and return the resulting vector_t */
__host__ __device__ vector_t vectorAdd(const vector_t v1, const vector_t v2){
  return vectorCreate(v1.x + v2.x, v1.y + v2.y, v1.z + v2.z);
}

__host__ __device__ dfloat vectorTripleProduct(const vector_t a, const vector_t b, const vector_t c){

  const vector_t aXb = vectorCrossProduct(a, b);
  
  return vectorDot(aXb, c); 
}

// assume b is unit vector
__host__ __device__ vector_t vectorOrthogonalize(const vector_t a, const vector_t b){

  dfloat adotb = vectorDot(a, b);

  return  vectorSub(a, vectorScale(adotb, b));
}

__host__ __device__ dfloat vectorNorm(const vector_t a){
  return  sqrt(vectorDot(a,a));
}

// return orthonormalized vector
__host__ __device__ vector_t vectorNormalize(const vector_t a){

  dfloat d = vectorNorm(a);
  if(d)
    return vectorScale(1./d, a);
  else
    return vectorCreate(0,0,0);
}

__host__ __device__ dfloat clamp(dfloat x, dfloat xmin, dfloat xmax){

  x = min(x, xmax);
  x = max(x, xmin);
  
  return x;
}


__host__ __device__ int iclamp(dfloat x, dfloat xmin, dfloat xmax){

  x = min(x, xmax);
  x = max(x, xmin);
  
  return floor(x);
}

// https://www.scratchapixel.com/code.php?id=10&origin=/lessons/3d-basic-rendering/minimal-ray-tracer-rendering-simple-shapes
// roots of a*t^2 + 2*b*t + c = 0
__host__ __device__ bool solveQuadratic(const dfloat a, const dfloat b, const dfloat c, dfloat *x0, dfloat *x1){
  
  dfloat discr = b * b - a * c;

  if (discr < 0) return false;
  else if (discr == 0) {
    x0[0] = x1[0] = - b / a;
  }
  else {
    dfloat sqrtdiscr = sqrt(discr);
    dfloat q = (b > 0) ?
      -(b + sqrtdiscr) :
      -(b - sqrtdiscr);
    x0[0] = q / a;
    x1[0] = c / q;
  }

  dfloat xmin = min(x0[0], x1[0]);
  dfloat xmax = max(x0[0], x1[0]);
  x0[0] = xmin;
  x1[0] = xmax;
  
  return true; 
}

