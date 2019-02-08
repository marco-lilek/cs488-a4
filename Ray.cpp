#include "Ray.hpp"

std::ostream &operator<<(std::ostream &os, Ray const &ray) {
  return os << "r{a=" << glm::to_string(ray.a) << ",b=" << glm::to_string(ray.b);
}

