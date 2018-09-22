#include "simpleRayTracer.h"

// HW02Q2: complete this function
bool intersectRayDisk(const ray_t r, const disk_t disk, dfloat *pt_t)
{
 
  dfloat t = 1e9;
  vector_t s = r.start;
  vector_t d = r.dir;
  vector_t n = disk.normal;
  vector_t c = disk.center;

   double tol = 1e-8;
   double num = vectorDot(n,c) - vectorDot(n, s);
  double den = vectorDot(n, d);
  
  if(fabs(den) < tol)
    {
      return 0;
    }

  if(num/den < 0)
    {
      return 0;
    }

  vector_t test_pt = vectorAdd(vectorScale(1., s), vectorScale((num/den), d));

  if (vectorNorm(vectorSub(test_pt, c)) > disk.radius)
    {
      return 0;
    }

  if((num/den) < *pt_t)
	{
	  *pt_t = (num/den);
	  return true;
	}
	return false;
      

}

// HW02Q3: complete this function
bool intersectRayCylinder(const ray_t r, const cylinder_t cylinder, dfloat *pt_t){
 
  vector_t s = r.start;
  vector_t d = r.dir;

  vector_t c = cylinder.center;
  vector_t a = cylinder.axis;
  dfloat   R = cylinder.radius;
  dfloat   H = cylinder.height;

  dfloat t = 1e9;
  double tol = 1e-8;

  double ray_dist;
  vector_t h_term;
  vector_t const_term;
  vector_t projected_height;
  vector_t test_pt;

  h_term = vectorSub(vectorScale(1., d),vectorScale(vectorDot(d, a),a));
  const_term = vectorSub(vectorScale(1., vectorSub(s,c)), vectorScale(vectorDot(a,vectorSub(s,c)),a));
  
  double Aquad = vectorDot(h_term,h_term);
  double Bquad = 2*vectorDot(h_term,const_term);
  double Cquad = vectorDot(const_term,const_term) - R*R;
 
  

  if(Bquad*Bquad - 4*Aquad*Cquad < 0.)
    return false;

  if(-1*Bquad + sqrt(fabs(Bquad*Bquad - 4*Aquad*Cquad)) < tol)
    return false;

  if(-1*Bquad > sqrt(Bquad*Bquad - 4*Aquad*Cquad))
    {
    ray_dist = (-1*sqrt(Bquad*Bquad - 4*Aquad*Cquad) - Bquad)/(2*Aquad);
    test_pt = vectorAdd(vectorScale(1, s),vectorScale(ray_dist, d));
    projected_height = vectorSub(test_pt, c);

    if(vectorDot(projected_height, a) > 0 &&
       vectorDot(projected_height, a) < H)
      {
	t = ray_dist;
      }
    }

  if(t< (*pt_t)){
    *pt_t = t;
    return true;
  }

  return false;

}





  
  // B. If no intersection return false
  
  // C. If there is an intersection then check is t< (*pt_t)
  //    a. if false: return false
  //    b. if true: set *pt_t = t, and return true
  /*
  if(t< (*pt_t)){
    *pt_t = t;
    return true;
  }

  return false;
  */


// HW02Q4: complete this function
bool intersectRayCone(const ray_t r, const cone_t cone, dfloat *pt_t)
{

  vector_t s = r.start;
  vector_t d = r.dir;

  vector_t v = cone.vertex;
  vector_t a = cone.axis;
  dfloat   R = cone.radius;
  dfloat   H = cone.height;

  vector_t term = vectorSub(s, v);
  dfloat t = 1e9;
  
  

  
  double ray_dist;
  vector_t projected_height;
  vector_t test_pt;
  double tol = 1e-8;

  vector_t w = vectorSub(d, vectorScale(vectorDot(a,d), a));
  vector_t x = vectorSub(term, vectorScale(vectorDot(a, term), a));
  double y = (R/H)*(vectorDot(d,a));
  double z = (R/H)*((vectorDot(s, a) - vectorDot(v, a)));

  double Qa = vectorDot(w, w) - (y*y);
  double Qb = 2*(vectorDot(w, x) - (y*z));
  double Qc = vectorDot(x, x) - (z*z);

  if(Qb*Qb - 4*Qa*Qc < 0.)
    {
      return false;
    }

  if(-1*Qb + sqrt(fabs(Qb*Qb - 4*Qa*Qc)) < tol)
    {
      return false;
    }

  if(-1*Qb > sqrt(Qb *Qb - 4*Qa*Qc))
    {
      t = (-1*sqrt(Qb *Qb - 4*Qa*Qc) - Qb)/(2*Qa);

      test_pt = vectorAdd(s, vectorScale(t, d));

      projected_height = vectorSub(test_pt, v);

      if(vectorDot(projected_height, a) > 0 && vectorDot(projected_height, a) < H)
	{
	  if(t < (*pt_t))
	    {
	      *pt_t = t;
	      return true;
	    }
	}
      
    }
   t = (sqrt(Qb *Qb - 4*Qa*Qc) - Qb)/(2*Qa);
  test_pt = vectorAdd(s, vectorScale(t, d));
  projected_height = vectorSub(test_pt, v);
  
       if(vectorDot(projected_height, a) > 0 && vectorDot(projected_height, a) < H)
  	{
    if(t < (*pt_t))
  	    {
  	      *pt_t = t;
 	      return true;
  	    }
  	}
  
    return false;
  
    if(t <(*pt_t))
  	{
  	  *pt_t = t;
  	  return true;
  	}
      return false;
      
      
  
}
 

// Do not edit beyond here--------------------------------------------------------------------------->


/* Check if the ray and triangle intersect */
bool intersectRayTriangle(const ray_t r, const triangle_t tri, dfloat *t){

  // TW: unused fudge factor
  dfloat delta = 0; 
  
  bool retval = false;

  vector_t B1 = vectorSub(tri.vertices[2], tri.vertices[0]);
  vector_t B2 = vectorSub(tri.vertices[2], tri.vertices[1]);
  vector_t B3 = r.dir;

  vector_t R = vectorSub(tri.vertices[2], r.start);

  dfloat J = vectorTripleProduct(B2, B3, B1);
  
  dfloat L1 = vectorTripleProduct(B2, B3, R);
  if(L1<delta*J) return false;

  dfloat L2 = vectorTripleProduct(B3, B1, R);
  if(L2<delta*J || L1+L2>J*(1+delta)) return false;

  dfloat t0 = vectorTripleProduct(B1, B2, R)/J;

  /* Verify t1 larger than 0 and less than the original t */
  // TW: FUDGE FACTOR
  if((t0 > p_intersectDelta) && (t0 < *t)){
    *t = t0;
    retval = true;
  }

  return retval;
}

/* Check if the ray and triangle intersect */
bool intersectRayRectangle(const ray_t r, const rectangle_t rect, dfloat *t){

  vector_t C  = rect.center;
  vector_t A1 = rect.axis[0];
  vector_t A2 = rect.axis[1];
  dfloat   L1 = rect.length[0];
  dfloat   L2 = rect.length[1];

  // n = A1 x A2
  // (s + t*d - C).n  = 0
  // t = (C - s).n/(d.n)

  vector_t n = vectorCrossProduct(A1, A2);
  dfloat  t0 = vectorDot(vectorSub(C,r.start), n)/vectorDot(r.dir, n);

  // intersection behind start of ray
  if(t0<0 || t0>*t) return false;

  // X = s + t*d - C
  vector_t X = vectorAdd(vectorSub(r.start,C), vectorScale(t0, r.dir));
  
  dfloat h1 = vectorDot(A1, X)+0.5*L1; // shift
  if(h1<0 || h1>L1) return false;
  
  dfloat h2 = vectorDot(A2, X)+0.5*L2; // shift
  if(h2<0 || h2>L1) return false;

  // success
  *t = t0;
  
  return true;
}



/* Check if the ray and sphere intersect */
bool intersectRaySphere(const ray_t r, const sphere_t s, dfloat *t)
{
	
    bool retval = false;

  /* A = d.d, the vector_t dot product of the direction */
  dfloat A = vectorDot(r.dir, r.dir); 
	
  /* We need a vector_t representing the distance between the start of 
   * the ray and the position of the circle.
   * This is the term (p0 - c) 
   */
  vector_t dist = vectorSub(r.start, s.pos);
	
  /* 2d.(p0 - c) */  
  dfloat B = 2.f * vectorDot(r.dir, dist);
	
  /* (p0 - c).(p0 - c) - r^2 */
  dfloat C = vectorDot(dist, dist) - (s.radius * s.radius);

  /* find roots of quadratic */
  dfloat t0, t1;
  if(solveQuadratic(A,0.5*B,C,&t0,&t1)){
    if((t0 > p_intersectDelta) && (t0 < *t)){
      *t = t0;
      retval = true;
    }else
      retval = false;
  }else{
    retval = false;
  }

  return retval;
}


bool intersectPointGridCell(const grid_t   grid,
			    const vector_t p,
			    const int cellI,
			    const int cellJ,
			    const int cellK){
  
  if(p.x<=grid.xmin+(cellI  )*grid.dx) return false;
  if(p.x> grid.xmin+(cellI+1)*grid.dx) return false;

  if(p.y<=grid.ymin+(cellJ  )*grid.dy) return false;
  if(p.y> grid.ymin+(cellJ+1)*grid.dy) return false;
  
  if(p.z<=grid.zmin+(cellK  )*grid.dz) return false;
  if(p.z> grid.zmin+(cellK+1)*grid.dz) return false;  
  
  return true;
}

unsigned int intersectRayBox(ray_t *r, const bbox_t bbox){

  vector_t d = r->dir;
  vector_t s = r->start;

  dfloat mint = 20000;
  unsigned int face = 0;
  
  if(d.x>0){ // face 2
    dfloat newt = (bbox.xmax-s.x)/d.x; // d.x > 0
    if(newt>0){
      mint = min(mint, newt);
    }
  }

  if(d.x<0){ // face 4
    // s.x + newt*d.x = bbox.xmin
    dfloat newt = (bbox.xmin-s.x)/d.x;
    if(newt>0){
      mint = min(mint, newt);
    }
  }
  
  if(d.y>0){ // face 3
    dfloat newt = (bbox.ymax-s.y)/d.y;
    if(newt>0){
      mint = min(mint, newt);
    }
  }

  if(d.y<0){ // face 1
    dfloat newt = (bbox.ymin-s.y)/d.y;
    if(newt>0){
      mint = min(mint, newt);
    }
  }

  if(d.z>0){ // face 5
    dfloat newt = (bbox.zmax-s.z)/d.z;
    if(newt>0){
      mint = min(mint, newt);
    }
  }

  if(d.z<0){ // face 0
    dfloat newt = (bbox.zmin-s.z)/d.z;
    if(newt>0){
      mint = min(mint, newt);
    }
  }

  // now figure out which faces the ray passes through
  if(d.x>0){ // face 2
    dfloat newt = (bbox.xmax-s.x)/d.x;
    if(newt>0 && newt<=mint)
      face |= 4;
  }

  if(d.x<0){ // face 4
    dfloat newt = (bbox.xmin-s.x)/d.x;
    if(newt>0 && newt<=mint)
      face |= 16;
  }

  if(d.y>0){ // face 3
    dfloat newt = (bbox.ymax-s.y)/d.y;
    if(newt>0 && newt<=mint)
      face |= 8;
  }

  if(d.y<0){ // face 1
    dfloat newt = (bbox.ymin-s.y)/d.y;
    if(newt>0 && newt<=mint)
      face |= 2;
  }

  if(d.z>0){ // face 5
    dfloat newt = (bbox.zmax-s.z)/d.z;
    if(newt>0 && newt<=mint)
      face |= 32;
  }

  if(d.z<0){ // face 0
    dfloat newt = (bbox.zmin-s.z)/d.z;
    if(newt>0 && newt<=mint)
      face |= 1;
  }
  
  if(face>0){
    // update start of ray
    r->start = vectorAdd(s, vectorScale(mint, d));
  }

  return face;
}

bool intersectRayShape(const ray_t r, const shape_t s, dfloat *t){

  switch(s.type){
  case SPHERE:   return intersectRaySphere  (r, s.sphere,   t);
  case CONE:     return intersectRayCone    (r, s.cone,     t);
  case DISK:     return intersectRayDisk    (r, s.disk,     t);
  case CYLINDER: return intersectRayCylinder(r, s.cylinder, t);
  case RECTANGLE:return intersectRayRectangle(r, s.rectangle, t);
  case TRIANGLE: return intersectRayTriangle(r, s.triangle, t); 
  }

  return false;

  
}

