#pragma once

#include "SceneNode.hpp"
#include "Primitive.hpp"
#include "Material.hpp"

#include "A4.hpp"

class GeometryNode : public SceneNode {
public:
	GeometryNode( const std::string & name, Primitive *prim, 
		Material *mat = nullptr );

	void setMaterial( Material *material );
  bool intersect(const Ray &r, double &t, Vector &N, Color &kd, Color &ks, double &shininess);

	Material *m_material;
	Primitive *m_primitive;
};
