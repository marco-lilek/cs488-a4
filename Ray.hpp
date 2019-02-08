#pragma once

#include <glm/glm.hpp>
#include <glm/ext.hpp>
#include <iostream>

typedef glm::dvec4 Point;
typedef glm::dvec4 Vector;
typedef glm::dvec3 Color;

struct Ray {
  Point a; // From
  Point b; // To

  Point pointAt(double t) const {
    return a + (b - a) * t;
  }

  Ray transform(const glm::mat4 &tmat) const {
    Ray r;
    r.a = tmat * a;
    r.b = tmat * b;
    return r;
  }
};

std::ostream &operator<<(std::ostream &os, Ray const &ray);
