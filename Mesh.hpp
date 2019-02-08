#pragma once

#include <vector>
#include <iosfwd>
#include <string>

#include <glm/glm.hpp>

#include "Primitive.hpp"

struct Triangle
{
	size_t v1;
	size_t v2;
	size_t v3;

	Triangle( size_t pv1, size_t pv2, size_t pv3 )
		: v1( pv1 )
		, v2( pv2 )
		, v3( pv3 )
	{}
};

// A polygonal mesh.
class Mesh : public Primitive {
public:
  Mesh( const std::string& fname );
  
  virtual bool intersect(const Ray &r, double &t, Vector &N);
  virtual glm::mat4 getTrans() { return glm::mat4(); }
private:
	std::vector<glm::vec3> m_vertices;
	std::vector<Triangle> m_faces;

  glm::mat4 bMat;
  glm::mat4 invbMat;
  glm::mat4 invtbMat;
  glm::vec3 m_boundingSphereCenter;
  double m_boundingRadius;

  std::string fname;
  void generateBoundingVolume();
  bool intersectBoundingVolume(const Ray &r, double &t, Vector &N);

    friend std::ostream& operator<<(std::ostream& out, const Mesh& mesh);
};
