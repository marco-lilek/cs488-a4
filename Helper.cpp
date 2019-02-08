
#include "Helper.hpp"

bool equals(double a, double b) {
  return glm::abs(a - b) < EPSILON;
}

// a >= b
bool ge(double a, double b) {
  return glm::abs(a - b) < EPSILON || a - b > 0;
}

// a <= b <= c
bool inb(double a, double b, double c) {
  return ge(b,a) && ge(c,b);
}
