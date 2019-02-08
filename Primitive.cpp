#include "Primitive.hpp"

#include "polyroots.hpp"
#include "Helper.hpp"
#include <iostream>

using namespace std;

// intersect with sphere centered at 0,0,0
bool intersectSphere(const Ray &r, double &t, Vector &N, const double radius, bool isBoundingVol) {
  // NOTE: Making the assumption that c = 0
  glm::dvec3 a(r.a);
  glm::dvec3 b(r.b);
  glm::dvec3 d(r.b-r.a);
  double A = glm::dot(d, d);
  double B = 2 * glm::dot(d, a);
  double C = glm::dot(a, a) - radius * radius;

// cerr << r << endl;
//  cerr << "a " << glm::to_string(a) << endl;
//  cerr << "b " << glm::to_string(b) << endl;
  double roots[2];
  size_t numRoots = quadraticRoots(A, B, C, roots);
// cerr << numRoots << endl;
  if (numRoots == 0) return false;
  if (numRoots == 2) {
    if (roots[0] < 0 && roots[1] < 0) return false;
    else if (roots[0] < 0 || roots[1] < 0) {
      double m =  glm::max(roots[0],roots[1]);
#ifdef DRAW_BOUNDING_VOL
      //if (isBoundingVol) return false;
#endif
      t = m;
    } else {
      t = glm::min(roots[0], roots[1]); // Intersection thats closer to the origin
    }
  } else if (numRoots == 1) {
    t = roots[0]; // TODO: assuming this is right??
  }
  //cerr << roots[0] << " " << roots[1] << endl;

  Point P = r.pointAt(t);

  N = glm::normalize(P - Point(0,0,0,1)); // from center to point on the sphere
  return true;
}

// intersect with a box with corner at 0,0,0
static bool intersectBox(const Ray &r, double &t, Vector &N, const double s) {
  // Intersect with only the front face of a cube
  glm::dvec3 e(r.a);
  glm::dvec3 d(r.b - r.a);

  // x,y,z = 0,1,2
  static const int planes[6][2] = {{0,1},{0,1},{1,2},{1,2},{0,2},{0,2}};
  static const double corner[6][3] = {
    {0,0,0},{0,0,1},{1,0,0},{0,0,0},{0,0,0},{0,1,0} // TODO: save space
  };
  static const double normals[6][3] = {
    {0,0,-1},{0,0,1},{1,0,0},{-1,0,0},{0,-1,0},{0,1,0}
  };

  bool intersect = false;
  for (int i = 0; i < 6; i++) {
    glm::dvec3 n(normals[i][0], normals[i][1], normals[i][2]);
    glm::dvec3 p1(corner[i][0] * s, corner[i][1] * s, corner[i][2] * s);

    double dn = glm::dot(d, n);
  //  cerr << dn <<  endl;
  //  cerr << glm::to_string(d) << " " << glm::to_string(n) << endl;
    if (equals(dn, 0)) continue; // Lies on the plane, so just say theres no intersection

    double lt = glm::dot(p1-e, n) / dn;
    //  cerr << "pass 1" << endl;
    if (lt < 0) continue; // No intersection
   //   cerr << "pass 2" << endl;
    if (intersect && lt > t) continue; // Ray already intersected the cube at a closer point

    Point P = r.pointAt(lt);

    double pp0 = P[planes[i][0]];
    double pp1 = P[planes[i][1]];
 //   cerr << pp0 << " " << pp1 << endl;
    if (!(pp0 < 0 || pp0 > s || pp1 < 0 || pp1 > s)) {
      intersect = true;
      t = lt;
      N = Vector(n, 0);
//      cerr << "yes" << endl;
    }
  }

  return intersect;
}

Primitive::~Primitive()
{
}

bool Sphere::intersect(const Ray &r, double &t, Vector &N) { 
  return intersectSphere(r,t,N, 1);
}

Sphere::~Sphere()
{
}

Cube::~Cube()
{
}

bool Cube::intersect(const Ray &r, double &t, Vector &N) { 
  return intersectBox(r,t,N, 1);
}

NonhierSphere::~NonhierSphere()
{
}

bool NonhierSphere::intersect(const Ray &r, double &t, Vector &N) {
  return intersectSphere(r,t,N, m_radius);
}

NonhierBox::~NonhierBox()
{
}

bool NonhierBox::intersect(const Ray &r, double &t, Vector &N) { 
  return intersectBox(r,t,N, m_size);
}
