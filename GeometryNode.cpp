#include "GeometryNode.hpp"

#include <iostream>

using namespace std;
//---------------------------------------------------------------------------------------
GeometryNode::GeometryNode(
	const std::string & name, Primitive *prim, Material *mat )
	: SceneNode( name )
	, m_material( mat )
	, m_primitive( prim )
{
	m_nodeType = NodeType::GeometryNode;
  //cerr << glm::to_string(trans) << endl; // TODO: this is wrong for now
  trans = prim->getTrans(); // Load in transforms for non-hier nodes
}

void GeometryNode::setMaterial( Material *mat )
{
	// Obviously, there's a potential memory leak here.  A good solution
	// would be to use some kind of reference counting, as in the 
	// C++ shared_ptr.  But I'm going to punt on that problem here.
	// Why?  Two reasons:
	// (a) In practice we expect the scene to be constructed exactly
	//     once.  There's no reason to believe that materials will be
	//     repeatedly overwritten in a GeometryNode.
	// (b) A ray tracer is a program in which you compute once, and 
	//     throw away all your data.  A memory leak won't build up and
	//     crash the program.

	m_material = mat;
}

// intersect w node, true if successful
bool GeometryNode::intersect(const Ray &r, double &t, Vector &N, Color &kd, Color &ks, double &shininess) {
  // TODO: hierarchy
 // cerr << "doing geom " << m_name << endl;
  //cerr << glm::to_string(trans) << endl;
  //cerr << glm::to_string(invttrans) << endl;
  Ray rp = r.transform(invtrans);
//  cerr << "result " << rp << endl;
  if (m_primitive->intersect(rp, t, N)) { // TODO: tangent line still counts
    N = invttrans * N;
    m_material->get(kd, ks, shininess);
    return true;
  }
  return false;
}
