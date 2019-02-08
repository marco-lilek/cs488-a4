#pragma once

#include <glm/glm.hpp>

#include "A4.hpp"

bool intersectSphere(const Ray &r, double &t, Vector &N, const double radius, bool drawBoundingVol = false);

class Primitive {
public:
  virtual ~Primitive();
  virtual bool intersect(const Ray &r, double &t, Vector &N) = 0;
  virtual glm::mat4 getTrans() = 0;
};

class Sphere : public Primitive {
public:
  virtual ~Sphere();
  virtual bool intersect(const Ray &r, double &t, Vector &N);
  virtual glm::mat4 getTrans() { return glm::mat4(); }
};

class Cube : public Primitive {
public:
  virtual ~Cube();
  virtual bool intersect(const Ray &r, double &t, Vector &N);
  virtual glm::mat4 getTrans() { return glm::mat4(); }
};

class NonhierSphere : public Primitive {
public:
  NonhierSphere(const glm::vec3& pos, double radius)
    : m_pos(pos), m_radius(radius)
  {
  }
  virtual ~NonhierSphere();
  virtual bool intersect(const Ray &r, double &t, Vector &N);
  virtual glm::mat4 getTrans() { return glm::translate(m_pos); }

private:
  glm::vec3 m_pos;
  double m_radius;
};

class NonhierBox : public Primitive {
public:
  NonhierBox(const glm::vec3& pos, double size)
    : m_pos(pos), m_size(size)
  {
  }
  
  virtual ~NonhierBox();
  virtual bool intersect(const Ray &r, double &t, Vector &N);
  virtual glm::mat4 getTrans() { return glm::translate(m_pos); }

private:
  glm::vec3 m_pos;
  double m_size;
};
