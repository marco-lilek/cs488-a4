
#include "A4.hpp"
#include "GeometryNode.hpp"
#include "Helper.hpp"
#include <vector>

using namespace std;

// build matrix to go from xy on screen to world pos
void buildScreenToWorld(
    const double fov,
    const double nx,
    const double ny,
    const glm::vec3 &lookFrom, 
    const glm::vec3 &lookAt, 
    const glm::vec3 &up,
    glm::mat4 &stw) {
  const double d = glm::distance(lookFrom, lookAt);
  const double h = 2 * d * glm::tan(glm::radians(fov) / 2);
  const double w = ( nx /  ny) * h;

  const glm::mat4 transToView = glm::translate(glm::mat4(), glm::vec3(-(nx / 2), -(ny / 2), -d)); // -d since glm::lookAt looks down the -z axis
  const glm::mat4 aspectRatio = glm::scale(glm::mat4(), glm::vec3(w / nx, -h / ny, 1)); // LHS coord system
  const glm::mat4 invLook = glm::inverse(glm::lookAt(lookFrom, lookAt, up));

  stw = invLook * aspectRatio * transToView;
}

static const double FIRE_OFFSET = 0.0000003;

// build a ray
void buildRay(const glm::mat4 stw, const glm::vec3 &lookFrom, const glm::vec2 &sP, Ray &ray) {
  // World coord P
  ray.a = glm::vec4(lookFrom, 1);
  ray.b = stw * glm::vec4(sP,0,1);
}

// ------ RayTracing algo -----------
bool hit(const SceneNode * root, 
    const Ray &r, 
    double &t,
    Vector &N, Color &kd, Color &ks, double &shininess);

// returns true if no hits
bool castShadowRay(const SceneNode *root, const Ray &r, double &t) {
 // cerr << "shadowray" << endl;
  t = 1;
  Vector gN; Color gkd, gks; double gt, gs;
//  cerr << "shadow!" << endl;
  bool hitres = hit(root, r, gt, gN, gkd, gks, gs);
//  cerr << hitres << " " << gt << endl;
  return !hitres || ge(gt, 1);
}

// do phong liting model
Color phongLight(
    const SceneNode *root,
    const std::list<Light *> &lights,
    const Ray &r,
    const Point &P, 
    const Vector &N,
    const Color &kd,
    const Color &ks,
    const double shininess) {
  if ( !equals(glm::length(N), 1)) {
    cerr << glm::to_string(N) << endl;
    cerr << glm::length(N) << endl;
    assert(0);
  }
  Color col(0);
  for (Light * light : lights) {
    Ray sr;
    sr.a = P;
    
    sr.b = Point(light->position, 1);
    sr.a = sr.pointAt(FIRE_OFFSET * 1000);

    glm::dvec3 l = glm::dvec3(sr.b - sr.a);
    glm::dvec3 ln = glm::normalize(l);
    
    double t = 1; // NOTE: for now just assume no falloff
    if (castShadowRay(root, sr, t)) {
      double dv = glm::dot(ln, glm::dvec3(N));
      double atten = 
        (light->falloff[0] + 
         light->falloff[1] * t + 
         light->falloff[2] * t * t);
    //  cerr << glm::to_string(ln) << " " << glm::to_string(N) << endl;
   //   cerr << dv << " " << atten << endl;
      //assert(dv < 0);
      Color I = light->colour * (glm::max(0.0, dv) / atten);
      col += kd* I;
   //   cerr << glm::to_string(dcol) << endl;

      if (glm::dot(ks,ks) != 0 && !equals(dv,0)) {
        glm::dvec3 rr = glm::normalize(-ln + 2 * dv * glm::dvec3(N)); // Reflection ray
        glm::dvec3 v = glm::normalize(glm::dvec3(r.a - P));
        //cerr << "aa" << glm::to_string(rr) << endl;
        //cerr << "bb" << glm::to_string(v) << endl;
        double sfactor = glm::pow(glm::max(0.0, glm::dot(v,rr)), shininess) / dv;

        Color scol = ks * sfactor * I;
        col += scol;
    //    cerr << glm::to_string(scol) << endl;
      }

      //col = Color(1,0,1);
    }
  //  cerr << "donecol" << endl;
  }

  return col;
}

// returns true if hit, sets params
bool hit(const SceneNode * root, 
    const Ray &r, 
    double &t,
    Vector &N, Color &kd, Color &ks, double &shininess) {
  //cerr << "RR" << r << endl;
  bool intersects = false;
  // TODO assuming its flat
  t = 1;
//  cerr << "transform by " << root->m_name << endl;
  Ray cr(r.transform(root->invtrans));
//  cerr << cr << endl;
  for (SceneNode *child : root->children) {
    double ct;
    Color ckd,cks;
    Vector cN;
    double cshininess;
 //   cerr << "down " << child->m_name << endl;
    if (child->children.size() > 0 && hit(child, cr, ct, cN, ckd, cks, cshininess) && ct > 0 && (!intersects || ct < t)) {
      intersects = true;
      t = ct; N = cN; kd = ckd; ks = cks; shininess = cshininess;
    }
   // cerr << "up " << child->m_name << endl;

    //cerr << child->m_name << endl;
    if (child->m_nodeType == NodeType::GeometryNode) {
      GeometryNode *gNode = static_cast<GeometryNode *>(child);
      double gt;
      Vector gN;
      Color gkd, gks; double gshininess;
//  cerr <<cr << endl;
      if (gNode->intersect(cr, gt, gN, gkd, gks, gshininess) && gt > 0 && (!intersects || gt < t)) {
        t = gt; N = gN; kd = gkd; ks = gks; shininess = gshininess;
        intersects = true;
       // cerr << glm::to_string(N) << endl;
      }
   //   cerr << intersects << " " << gt << " " << t << endl;
      //cerr << t << " " << gt << endl;
      //cerr << intersects << endl;
    }
  }
//  cerr << "end" << endl;

  N = root->invttrans * N;
  //if (intersects)
    //cerr << glm::to_string(glm::normalize(N)) << endl;
  return intersects;
}

// set the color of the ray
bool rayColor(
    const SceneNode *root,
    const Color &ambient, 
    const std::list<Light *> &lights,
    const Ray &r, const size_t &maxHits,
    Color &resColor) {
  resColor = Color(0);
  Color kd(0),ks(0);
  double shininess(0);
  Vector N(0);
  Point P;
  double t;
//  cerr << "raycolor" << endl;
  if (hit(root, r, t, N, kd, ks, shininess)) {
    //cerr << glm::to_string(phongLight(root, lights, r.pointAt(t), N, kd, ks)) << endl;
    resColor = kd * ambient + phongLight(root, lights, r, r.pointAt(t - FIRE_OFFSET), glm::normalize(N), kd, ks, shininess);
    return true;
  }
  return false;
}

static const size_t MAXITER_MANDELBROT = 1000;
static std::vector<Color> palette;

// call before running ray tracing
void setupBG() {
  palette.resize(MAXITER_MANDELBROT);
  for (size_t i = 0; i < palette.size(); i++) {
    double v = (double) i / MAXITER_MANDELBROT;
    palette[i] = Color(0,v*0.7,0.3);
  }
}

// get background color at pos
Color getBackgroundAt(const size_t x, const size_t y, const size_t w, const size_t h) {
  double xn = (double)x / w * 0.05 - 0.5;
  double yn = (double)y / w * 0.05 + 0.58;
 // cerr << xn << " " << yn << endl;
  size_t k;
  double u(0),v(0),u2(0),v2(0);
  for (k = 1; k < MAXITER_MANDELBROT && (u2 + v2 < (1 << 16)); k++) {
    v = 2 * u * v + yn;
    u = u2 - v2 + xn;
    u2 = u * u;
    v2 = v* v;
  }

  if (k >= MAXITER_MANDELBROT) return Color(0,0,0);
  else {
    return palette[k];
  }
}

static const size_t SSREG = 4;

// does rendering
void A4_Render(
		// What to render
		SceneNode * root,

		// Image to write to, set to a given width and height
		Image & image,

		// Viewing parameters
		const glm::vec3 & eye,
		const glm::vec3 & view,
		const glm::vec3 & up,
		double fovy,

		// Lighting parameters
		const glm::vec3 & ambient,
		const std::list<Light *> & lights
) {

  // Fill in raytracing code here...

  std::cout << "Calling A4_Render(\n" <<
		  "\t" << *root <<
          "\t" << "Image(width:" << image.width() << ", height:" << image.height() << ")\n"
          "\t" << "eye:  " << glm::to_string(eye) << std::endl <<
		  "\t" << "view: " << glm::to_string(view) << std::endl <<
		  "\t" << "up:   " << glm::to_string(up) << std::endl <<
		  "\t" << "fovy: " << fovy << std::endl <<
          "\t" << "ambient: " << glm::to_string(ambient) << std::endl <<
		  "\t" << "lights{" << std::endl;

	for(const Light * light : lights) {
		std::cout << "\t\t" <<  *light << std::endl;
	}
	std::cout << "\t}" << std::endl;
	std:: cout <<")" << std::endl;

	size_t h = image.height();
	size_t w = image.width();
  glm::mat4 stw;
  buildScreenToWorld(fovy, w, h, eye, view, up, stw);
  Ray sampleR;
  buildRay(stw, eye, glm::vec2(w/2,h/2), sampleR);
  cerr << "sample " << sampleR << endl;

  root->bakeTransforms();
  setupBG();

  double progBoundary = 0.1;
  double progInc = 0.1;
  double A= w * h;

  cerr << endl << endl << "starting ray tracing!" << endl;
#ifdef SUPERSAMPLE
  cerr << "doing supersampling (" << SSREG*SSREG << " samples per pixels)" << endl;
#else
  cerr << "no supersampling" << endl;
#endif

#ifdef NO_BOUNDING_VOL
  cerr << "no bounding volumes" << endl;
#elif defined(DRAW_BOUNDING_VOL)
  cerr << "drawing bounding volumes" << endl;
#endif

  for (size_t x = 0; x < w; x++) {
    for (size_t y = 0; y < h; y++) {
 //     cerr << "------------" << x << " " << y << endl;

      size_t numsamples;
#ifdef SUPERSAMPLE
      numsamples = SSREG;
#else
      numsamples = 1;
#endif
      // Supersampling
      Color fincol(0); // avg
      bool isBg = false;
      for (size_t xr = 0; xr < numsamples; xr++) {
        for (size_t yr = 0; yr < numsamples; yr++) {
          double xray = x + (double) (2 * xr+ 1) / (2 * numsamples);
          double yray = y + (double) (2 * yr+1) / (2 * numsamples);
          Ray r;
          buildRay(stw, eye, glm::vec2(xray,yray), r);
          Color c;
          if (!rayColor(root,ambient,lights, r, 10, c)) {
            fincol = getBackgroundAt(x,y,w,h);
            isBg = true;
            break;
          } else {
            fincol += c / ((double)numsamples * numsamples);
          }
        }
        if (isBg) break;
      }

      image(x,y,0) = fincol.x;
      image(x,y,1) = fincol.y;
      image(x,y,2) = fincol.z;

      if ((double) (x*y) / A > progBoundary) {
        cerr << "progress: " << progBoundary*100 << "%" << endl;
        progBoundary += progInc;
      }

    }
  }
}
