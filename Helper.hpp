#pragma once

#include <glm/glm.hpp>
#include <glm/ext.hpp>

static const double EPSILON =  0.00000000001;

bool equals(double a, double b);

// a >= b
bool ge(double a, double b);

// a <= b <= c
bool inb(double a, double b, double c);
