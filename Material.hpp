#pragma once

#include "Ray.hpp"

class Material {
public:
  virtual ~Material();

  virtual void get(Color &kd, Color &ks, double &shininess) {
    kd = Color(0); 
    ks = Color(0); 
    shininess = 0;
  }

protected:
  Material();
};
