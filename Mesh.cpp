#include <iostream>
#include <fstream>

#include <glm/ext.hpp>
#include <iostream>

// #include "cs488-framework/ObjFileDecoder.hpp"
#include "Primitive.hpp"
#include "Mesh.hpp"
#include "Helper.hpp"

using namespace std;

Mesh::Mesh( const std::string& fname )
	: m_vertices()
	, m_faces(),
  bMat(glm::mat4()), invbMat(glm::mat4()), invtbMat(glm::mat4()), m_boundingSphereCenter(glm::vec3()), m_boundingRadius(0), fname(fname)
{
	std::string code;
	double vx, vy, vz;
	size_t s1, s2, s3;

	std::ifstream ifs( fname.c_str() );
	while( ifs >> code ) {
		if( code == "v" ) {
			ifs >> vx >> vy >> vz;
			m_vertices.push_back( glm::vec3( vx, vy, vz ) );
		} else if( code == "f" ) {
			ifs >> s1 >> s2 >> s3;
			m_faces.push_back( Triangle( s1 - 1, s2 - 1, s3 - 1 ) );
		}
	}

  generateBoundingVolume(); 
}

// builds sphere bv around the points of the mesh
void Mesh::generateBoundingVolume() {
  if (m_vertices.size() == 0) {
    m_boundingRadius = 0;
    return;
  }

  for (size_t vidx = 0; vidx < m_vertices.size(); vidx++) {
    m_boundingSphereCenter += m_vertices[vidx];
  }

  m_boundingSphereCenter = m_boundingSphereCenter / (double)m_vertices.size();

  for (size_t vidx = 0; vidx < m_vertices.size(); vidx++) {
    double vdistance = (double)glm::distance(m_boundingSphereCenter, m_vertices[vidx]);
    m_boundingRadius = glm::max(vdistance, m_boundingRadius);
  }
  
  cerr << glm::to_string(m_boundingSphereCenter) << endl;
  cerr << m_boundingRadius << endl;
  cerr << fname << endl;
  bMat = glm::translate(m_boundingSphereCenter);
  invbMat = glm::inverse(glm::translate(m_boundingSphereCenter));
  invtbMat = glm::dmat4(glm::inverseTranspose(glm::dmat3(bMat)));
}

// returns true if intersects w the bounding vol
bool Mesh::intersectBoundingVolume(const Ray &r, double &t, Vector &N) {
  Ray rp = r.transform(invbMat);
  if (intersectSphere(rp, t, N, m_boundingRadius, true)) {
    N = invtbMat * N;
    return true;
  }
  return false;
}

bool Mesh::intersect(const Ray &r, double &t, Vector &N) { 
#ifndef NO_BOUNDING_VOL
  if (m_vertices.size() > 0 && fname != "plane.obj") {
    double tbv = 1;
    Vector Nbv;
    bool bVolRes = intersectBoundingVolume(r,tbv,Nbv);

#ifdef DRAW_BOUNDING_VOL
    t = tbv;
    N = Nbv;
    return bVolRes;
#else
    if (!bVolRes) return false;
#endif

  }

#endif
  bool intersects = false;

  glm::vec3 ev(r.a);
  glm::vec3 dv(r.b - r.a);
  double a,b,c,d,e,f,g,h,i,j,k,l;
  for (size_t fidx = 0; fidx < m_faces.size(); ++fidx) { // taken from textbook
    glm::vec3 va(m_vertices[m_faces[fidx].v1]);
    glm::vec3 vb(m_vertices[m_faces[fidx].v2]);
    glm::vec3 vc(m_vertices[m_faces[fidx].v3]);

    a = va.x-vb.x; b = va.y-vb.y; c = va.z-vb.z;
    d = va.x-vc.x; e = va.y-vc.y; f = va.z-vc.z;
    g = dv.x; h = dv.y; i = dv.z;

    j = va.x - ev.x; k = va.y - ev.y; l = va.z - ev.z;

    double M = a*(e*i-h*f)+b*(g*f-d*i)+c*(d*h-e*g);

    double tp = - (f*(a*k-j*b)+e*(j*c-a*l)+d*(b*l-k*c)) / M;

    //cerr << tp << endl;
    if (tp < 0) continue; // Doesnt intersect this plane
    if (intersects && tp > t) continue; // We've already intersected at a nearer point

    double gam = (i*(a*k-j*b)+h*(j*c-a*l)+g*(b*l-k*c)) / M;
    if (gam < 0 || gam > 1) continue;

    double bet = (j*(e*i-h*f)+k*(g*f-d*i)+l*(d*h-e*g)) / M;

  //  cerr << fidx << " " << bet << " " << gam << endl;
    if (bet < 0 || (!equals(bet+gam, 1) && bet + gam > 1)) continue;

    t = tp;

    N = glm::vec4(glm::cross(vb - va, vc - vb), 0);
  //  cerr << glm::to_string(N) << " " << tp<< endl;
    intersects = true;
  }

  return intersects;
}

std::ostream& operator<<(std::ostream& out, const Mesh& mesh)
{
  out << "mesh {";
  /*
  
  for( size_t idx = 0; idx < mesh.m_verts.size(); ++idx ) {
  	const MeshVertex& v = mesh.m_verts[idx];
  	out << glm::to_string( v.m_position );
	if( mesh.m_have_norm ) {
  	  out << " / " << glm::to_string( v.m_normal );
	}
	if( mesh.m_have_uv ) {
  	  out << " / " << glm::to_string( v.m_uv );
	}
  }

*/
  out << "}";
  return out;
}
