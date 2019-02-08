#include "PhongMaterial.hpp"

PhongMaterial::PhongMaterial(
	const glm::vec3& kd, const glm::vec3& ks, double shininess )
	: m_kd(kd)
	, m_ks(ks)
	, m_shininess(shininess)
{}

PhongMaterial::~PhongMaterial()
{}

void PhongMaterial::get(Color &kd, Color &ks, double &shininess) {
  kd = m_kd; ks = m_ks; shininess = m_shininess;
}
