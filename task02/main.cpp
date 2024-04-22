#include <cmath>
#include <filesystem>
#include <tuple>
#include <utility>
// #include <experimental/filesystem> // uncomment here if the <filesystem>
// cannot be included above
//
#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "Eigen/Core"
#include "stb_image_write.h"
//
#include "parse_svg.h"

/***
 * signed area of a triangle connecting points (p0, p1, p2) in counter-clockwise
 * order.
 * @param p0 1st point xy-coordinate
 * @param p1 2nd point xy-coordinate
 * @param p2 3rd point xy-coordinate
 * @return signed area (float)
 */
float area(const Eigen::Vector2f &p0, const Eigen::Vector2f &p1,
           const Eigen::Vector2f &p2) {
  const auto v01 = p1 - p0;
  const auto v02 = p2 - p0;
  // return 0.5f * (v01[0] * v02[1] - v01[1] * v02[0]); // right handed
  // coordinate
  return 0.5f * (v01[1] * v02[0] -
                 v01[0] * v02[1]); // left-handed coordinate (because pixel
                                   // y-coordinate is going down)
}

/***
 * compute number of intersection of a ray against a line segment
 * @param org ray origin
 * @param dir ray direction (unit normal)
 * @param ps one of the two end points
 * @param pe the other end point
 * @return number of intersection
 */
int number_of_intersection_ray_against_edge(const Eigen::Vector2f &org,
                                            const Eigen::Vector2f &dir,
                                            const Eigen::Vector2f &ps,
                                            const Eigen::Vector2f &pe) {
  auto a = area(org, org + dir, ps);
  auto b = area(org, pe, org + dir);
  auto c = area(org, ps, pe);
  auto d = area(dir + ps, ps, pe);
  if (a * b > 0.f && d * c < 0.f) {
    return 1;
  }
  return 0;
  // the following code was a bug
  // auto d = area(org + dir, ps, pe);
  // if (a * b > 0.f && d * c > 0.f && fabs(d) > fabs(c)) { return 1; }
}

/***
 *
 * @param org ray origin
 * @param dir ray direction (unit vector)
 * @param ps one of the two end points
 * @param pc control point
 * @param pe the other end point
 * @return the number of intersections
 */
int number_of_intersection_ray_against_quadratic_bezier(
    const Eigen::Vector2f &org, const Eigen::Vector2f &dir,
    const Eigen::Vector2f &ps, const Eigen::Vector2f &pc,
    const Eigen::Vector2f &pe) {
  // comment out below to do the assignment

  /**
 // Mathematica
 Solve[(1 - t) * (1 - t) * {x0,y0} + 2 * (1 - t) * t * {x1,y1}  + t * t *
 {x2,y2}=={xo,yo}+s*{xd,yd},{s,t}] // CForm
  */
  // use double or picture will go wrong. :(
  // main idea: solve out the s and t, check s and t are in range, a.k.a valid.

  const auto Power = [&](double x, auto p) -> double { return x * x; };
  const auto Sqrt = [&](double x) -> double { return std::sqrt(x); };
  const double x0 = ps.x();
  const double y0 = ps.y();
  const double x1 = pc.x();
  const double y1 = pc.y();
  const double x2 = pe.x();
  const double y2 = pe.y();
  const double xd = dir.x();
  const double yd = dir.y();
  const double xo = org.x();
  const double yo = org.y();
  const double s1 =
      (2 * x1 * y0 - x2 * y0 - xo * y0 - 2 * x0 * y1 + 2 * xo * y1 + x0 * y2 -
       xo * y2 -
       (2 * x1 * xd * Power(y0, 2)) /
           (xd * y0 - 2 * xd * y1 + xd * y2 - x0 * yd + 2 * x1 * yd - x2 * yd) +
       (2 * x2 * xd * Power(y0, 2)) /
           (xd * y0 - 2 * xd * y1 + xd * y2 - x0 * yd + 2 * x1 * yd - x2 * yd) +
       (2 * x0 * xd * y0 * y1) /
           (xd * y0 - 2 * xd * y1 + xd * y2 - x0 * yd + 2 * x1 * yd - x2 * yd) +
       (2 * x1 * xd * y0 * y1) /
           (xd * y0 - 2 * xd * y1 + xd * y2 - x0 * yd + 2 * x1 * yd - x2 * yd) -
       (4 * x2 * xd * y0 * y1) /
           (xd * y0 - 2 * xd * y1 + xd * y2 - x0 * yd + 2 * x1 * yd - x2 * yd) -
       (2 * x0 * xd * Power(y1, 2)) /
           (xd * y0 - 2 * xd * y1 + xd * y2 - x0 * yd + 2 * x1 * yd - x2 * yd) +
       (2 * x2 * xd * Power(y1, 2)) /
           (xd * y0 - 2 * xd * y1 + xd * y2 - x0 * yd + 2 * x1 * yd - x2 * yd) -
       (2 * x0 * xd * y0 * y2) /
           (xd * y0 - 2 * xd * y1 + xd * y2 - x0 * yd + 2 * x1 * yd - x2 * yd) +
       (2 * x1 * xd * y0 * y2) /
           (xd * y0 - 2 * xd * y1 + xd * y2 - x0 * yd + 2 * x1 * yd - x2 * yd) +
       (2 * x0 * xd * y1 * y2) /
           (xd * y0 - 2 * xd * y1 + xd * y2 - x0 * yd + 2 * x1 * yd - x2 * yd) -
       (2 * x1 * xd * y1 * y2) /
           (xd * y0 - 2 * xd * y1 + xd * y2 - x0 * yd + 2 * x1 * yd - x2 * yd) +
       (2 * x0 * x1 * y0 * yd) /
           (xd * y0 - 2 * xd * y1 + xd * y2 - x0 * yd + 2 * x1 * yd - x2 * yd) -
       (2 * Power(x1, 2) * y0 * yd) /
           (xd * y0 - 2 * xd * y1 + xd * y2 - x0 * yd + 2 * x1 * yd - x2 * yd) -
       (2 * x0 * x2 * y0 * yd) /
           (xd * y0 - 2 * xd * y1 + xd * y2 - x0 * yd + 2 * x1 * yd - x2 * yd) +
       (2 * x1 * x2 * y0 * yd) /
           (xd * y0 - 2 * xd * y1 + xd * y2 - x0 * yd + 2 * x1 * yd - x2 * yd) -
       (2 * Power(x0, 2) * y1 * yd) /
           (xd * y0 - 2 * xd * y1 + xd * y2 - x0 * yd + 2 * x1 * yd - x2 * yd) +
       (2 * x0 * x1 * y1 * yd) /
           (xd * y0 - 2 * xd * y1 + xd * y2 - x0 * yd + 2 * x1 * yd - x2 * yd) +
       (2 * x0 * x2 * y1 * yd) /
           (xd * y0 - 2 * xd * y1 + xd * y2 - x0 * yd + 2 * x1 * yd - x2 * yd) -
       (2 * x1 * x2 * y1 * yd) /
           (xd * y0 - 2 * xd * y1 + xd * y2 - x0 * yd + 2 * x1 * yd - x2 * yd) +
       (2 * Power(x0, 2) * y2 * yd) /
           (xd * y0 - 2 * xd * y1 + xd * y2 - x0 * yd + 2 * x1 * yd - x2 * yd) -
       (4 * x0 * x1 * y2 * yd) /
           (xd * y0 - 2 * xd * y1 + xd * y2 - x0 * yd + 2 * x1 * yd - x2 * yd) +
       (2 * Power(x1, 2) * y2 * yd) /
           (xd * y0 - 2 * xd * y1 + xd * y2 - x0 * yd + 2 * x1 * yd - x2 * yd) +
       x0 * yo - 2 * x1 * yo + x2 * yo +
       (x1 * y0 *
        Sqrt(Power(-2 * xd * y0 + 2 * xd * y1 + 2 * x0 * yd - 2 * x1 * yd, 2) -
             4 *
                 (xd * y0 - 2 * xd * y1 + xd * y2 - x0 * yd + 2 * x1 * yd -
                  x2 * yd) *
                 (xd * y0 - x0 * yd + xo * yd - xd * yo))) /
           (xd * y0 - 2 * xd * y1 + xd * y2 - x0 * yd + 2 * x1 * yd - x2 * yd) -
       (x2 * y0 *
        Sqrt(Power(-2 * xd * y0 + 2 * xd * y1 + 2 * x0 * yd - 2 * x1 * yd, 2) -
             4 *
                 (xd * y0 - 2 * xd * y1 + xd * y2 - x0 * yd + 2 * x1 * yd -
                  x2 * yd) *
                 (xd * y0 - x0 * yd + xo * yd - xd * yo))) /
           (xd * y0 - 2 * xd * y1 + xd * y2 - x0 * yd + 2 * x1 * yd - x2 * yd) -
       (x0 * y1 *
        Sqrt(Power(-2 * xd * y0 + 2 * xd * y1 + 2 * x0 * yd - 2 * x1 * yd, 2) -
             4 *
                 (xd * y0 - 2 * xd * y1 + xd * y2 - x0 * yd + 2 * x1 * yd -
                  x2 * yd) *
                 (xd * y0 - x0 * yd + xo * yd - xd * yo))) /
           (xd * y0 - 2 * xd * y1 + xd * y2 - x0 * yd + 2 * x1 * yd - x2 * yd) +
       (x2 * y1 *
        Sqrt(Power(-2 * xd * y0 + 2 * xd * y1 + 2 * x0 * yd - 2 * x1 * yd, 2) -
             4 *
                 (xd * y0 - 2 * xd * y1 + xd * y2 - x0 * yd + 2 * x1 * yd -
                  x2 * yd) *
                 (xd * y0 - x0 * yd + xo * yd - xd * yo))) /
           (xd * y0 - 2 * xd * y1 + xd * y2 - x0 * yd + 2 * x1 * yd - x2 * yd) +
       (x0 * y2 *
        Sqrt(Power(-2 * xd * y0 + 2 * xd * y1 + 2 * x0 * yd - 2 * x1 * yd, 2) -
             4 *
                 (xd * y0 - 2 * xd * y1 + xd * y2 - x0 * yd + 2 * x1 * yd -
                  x2 * yd) *
                 (xd * y0 - x0 * yd + xo * yd - xd * yo))) /
           (xd * y0 - 2 * xd * y1 + xd * y2 - x0 * yd + 2 * x1 * yd - x2 * yd) -
       (x1 * y2 *
        Sqrt(Power(-2 * xd * y0 + 2 * xd * y1 + 2 * x0 * yd - 2 * x1 * yd, 2) -
             4 *
                 (xd * y0 - 2 * xd * y1 + xd * y2 - x0 * yd + 2 * x1 * yd -
                  x2 * yd) *
                 (xd * y0 - x0 * yd + xo * yd - xd * yo))) /
           (xd * y0 - 2 * xd * y1 + xd * y2 - x0 * yd + 2 * x1 * yd -
            x2 * yd)) /
      (xd * y0 - 2 * xd * y1 + xd * y2 - x0 * yd + 2 * x1 * yd - x2 * yd);
  const double t1 =
      (2 * xd * y0 - 2 * xd * y1 - 2 * x0 * yd + 2 * x1 * yd -
       Sqrt(Power(-2 * xd * y0 + 2 * xd * y1 + 2 * x0 * yd - 2 * x1 * yd, 2) -
            4 *
                (xd * y0 - 2 * xd * y1 + xd * y2 - x0 * yd + 2 * x1 * yd -
                 x2 * yd) *
                (xd * y0 - x0 * yd + xo * yd - xd * yo))) /
      (2. *
       (xd * y0 - 2 * xd * y1 + xd * y2 - x0 * yd + 2 * x1 * yd - x2 * yd));

  const double s2 =
      (2 * x1 * y0 - x2 * y0 - xo * y0 - 2 * x0 * y1 + 2 * xo * y1 + x0 * y2 -
       xo * y2 -
       (2 * x1 * xd * Power(y0, 2)) /
           (xd * y0 - 2 * xd * y1 + xd * y2 - x0 * yd + 2 * x1 * yd - x2 * yd) +
       (2 * x2 * xd * Power(y0, 2)) /
           (xd * y0 - 2 * xd * y1 + xd * y2 - x0 * yd + 2 * x1 * yd - x2 * yd) +
       (2 * x0 * xd * y0 * y1) /
           (xd * y0 - 2 * xd * y1 + xd * y2 - x0 * yd + 2 * x1 * yd - x2 * yd) +
       (2 * x1 * xd * y0 * y1) /
           (xd * y0 - 2 * xd * y1 + xd * y2 - x0 * yd + 2 * x1 * yd - x2 * yd) -
       (4 * x2 * xd * y0 * y1) /
           (xd * y0 - 2 * xd * y1 + xd * y2 - x0 * yd + 2 * x1 * yd - x2 * yd) -
       (2 * x0 * xd * Power(y1, 2)) /
           (xd * y0 - 2 * xd * y1 + xd * y2 - x0 * yd + 2 * x1 * yd - x2 * yd) +
       (2 * x2 * xd * Power(y1, 2)) /
           (xd * y0 - 2 * xd * y1 + xd * y2 - x0 * yd + 2 * x1 * yd - x2 * yd) -
       (2 * x0 * xd * y0 * y2) /
           (xd * y0 - 2 * xd * y1 + xd * y2 - x0 * yd + 2 * x1 * yd - x2 * yd) +
       (2 * x1 * xd * y0 * y2) /
           (xd * y0 - 2 * xd * y1 + xd * y2 - x0 * yd + 2 * x1 * yd - x2 * yd) +
       (2 * x0 * xd * y1 * y2) /
           (xd * y0 - 2 * xd * y1 + xd * y2 - x0 * yd + 2 * x1 * yd - x2 * yd) -
       (2 * x1 * xd * y1 * y2) /
           (xd * y0 - 2 * xd * y1 + xd * y2 - x0 * yd + 2 * x1 * yd - x2 * yd) +
       (2 * x0 * x1 * y0 * yd) /
           (xd * y0 - 2 * xd * y1 + xd * y2 - x0 * yd + 2 * x1 * yd - x2 * yd) -
       (2 * Power(x1, 2) * y0 * yd) /
           (xd * y0 - 2 * xd * y1 + xd * y2 - x0 * yd + 2 * x1 * yd - x2 * yd) -
       (2 * x0 * x2 * y0 * yd) /
           (xd * y0 - 2 * xd * y1 + xd * y2 - x0 * yd + 2 * x1 * yd - x2 * yd) +
       (2 * x1 * x2 * y0 * yd) /
           (xd * y0 - 2 * xd * y1 + xd * y2 - x0 * yd + 2 * x1 * yd - x2 * yd) -
       (2 * Power(x0, 2) * y1 * yd) /
           (xd * y0 - 2 * xd * y1 + xd * y2 - x0 * yd + 2 * x1 * yd - x2 * yd) +
       (2 * x0 * x1 * y1 * yd) /
           (xd * y0 - 2 * xd * y1 + xd * y2 - x0 * yd + 2 * x1 * yd - x2 * yd) +
       (2 * x0 * x2 * y1 * yd) /
           (xd * y0 - 2 * xd * y1 + xd * y2 - x0 * yd + 2 * x1 * yd - x2 * yd) -
       (2 * x1 * x2 * y1 * yd) /
           (xd * y0 - 2 * xd * y1 + xd * y2 - x0 * yd + 2 * x1 * yd - x2 * yd) +
       (2 * Power(x0, 2) * y2 * yd) /
           (xd * y0 - 2 * xd * y1 + xd * y2 - x0 * yd + 2 * x1 * yd - x2 * yd) -
       (4 * x0 * x1 * y2 * yd) /
           (xd * y0 - 2 * xd * y1 + xd * y2 - x0 * yd + 2 * x1 * yd - x2 * yd) +
       (2 * Power(x1, 2) * y2 * yd) /
           (xd * y0 - 2 * xd * y1 + xd * y2 - x0 * yd + 2 * x1 * yd - x2 * yd) +
       x0 * yo - 2 * x1 * yo + x2 * yo -
       (x1 * y0 *
        Sqrt(Power(-2 * xd * y0 + 2 * xd * y1 + 2 * x0 * yd - 2 * x1 * yd, 2) -
             4 *
                 (xd * y0 - 2 * xd * y1 + xd * y2 - x0 * yd + 2 * x1 * yd -
                  x2 * yd) *
                 (xd * y0 - x0 * yd + xo * yd - xd * yo))) /
           (xd * y0 - 2 * xd * y1 + xd * y2 - x0 * yd + 2 * x1 * yd - x2 * yd) +
       (x2 * y0 *
        Sqrt(Power(-2 * xd * y0 + 2 * xd * y1 + 2 * x0 * yd - 2 * x1 * yd, 2) -
             4 *
                 (xd * y0 - 2 * xd * y1 + xd * y2 - x0 * yd + 2 * x1 * yd -
                  x2 * yd) *
                 (xd * y0 - x0 * yd + xo * yd - xd * yo))) /
           (xd * y0 - 2 * xd * y1 + xd * y2 - x0 * yd + 2 * x1 * yd - x2 * yd) +
       (x0 * y1 *
        Sqrt(Power(-2 * xd * y0 + 2 * xd * y1 + 2 * x0 * yd - 2 * x1 * yd, 2) -
             4 *
                 (xd * y0 - 2 * xd * y1 + xd * y2 - x0 * yd + 2 * x1 * yd -
                  x2 * yd) *
                 (xd * y0 - x0 * yd + xo * yd - xd * yo))) /
           (xd * y0 - 2 * xd * y1 + xd * y2 - x0 * yd + 2 * x1 * yd - x2 * yd) -
       (x2 * y1 *
        Sqrt(Power(-2 * xd * y0 + 2 * xd * y1 + 2 * x0 * yd - 2 * x1 * yd, 2) -
             4 *
                 (xd * y0 - 2 * xd * y1 + xd * y2 - x0 * yd + 2 * x1 * yd -
                  x2 * yd) *
                 (xd * y0 - x0 * yd + xo * yd - xd * yo))) /
           (xd * y0 - 2 * xd * y1 + xd * y2 - x0 * yd + 2 * x1 * yd - x2 * yd) -
       (x0 * y2 *
        Sqrt(Power(-2 * xd * y0 + 2 * xd * y1 + 2 * x0 * yd - 2 * x1 * yd, 2) -
             4 *
                 (xd * y0 - 2 * xd * y1 + xd * y2 - x0 * yd + 2 * x1 * yd -
                  x2 * yd) *
                 (xd * y0 - x0 * yd + xo * yd - xd * yo))) /
           (xd * y0 - 2 * xd * y1 + xd * y2 - x0 * yd + 2 * x1 * yd - x2 * yd) +
       (x1 * y2 *
        Sqrt(Power(-2 * xd * y0 + 2 * xd * y1 + 2 * x0 * yd - 2 * x1 * yd, 2) -
             4 *
                 (xd * y0 - 2 * xd * y1 + xd * y2 - x0 * yd + 2 * x1 * yd -
                  x2 * yd) *
                 (xd * y0 - x0 * yd + xo * yd - xd * yo))) /
           (xd * y0 - 2 * xd * y1 + xd * y2 - x0 * yd + 2 * x1 * yd -
            x2 * yd)) /
      (xd * y0 - 2 * xd * y1 + xd * y2 - x0 * yd + 2 * x1 * yd - x2 * yd);
  const double t2 =
      (2 * xd * y0 - 2 * xd * y1 - 2 * x0 * yd + 2 * x1 * yd +
       Sqrt(Power(-2 * xd * y0 + 2 * xd * y1 + 2 * x0 * yd - 2 * x1 * yd, 2) -
            4 *
                (xd * y0 - 2 * xd * y1 + xd * y2 - x0 * yd + 2 * x1 * yd -
                 x2 * yd) *
                (xd * y0 - x0 * yd + xo * yd - xd * yo))) /
      (2. *
       (xd * y0 - 2 * xd * y1 + xd * y2 - x0 * yd + 2 * x1 * yd - x2 * yd));

  return (s1 >= 0 && 0 <= t1 && t1 <= 1) + (s2 >= 0 && 0 <= t2 && t2 <= 1);
}

int main() {
  const auto input_file_path =
      std::filesystem::path(PROJECT_SOURCE_DIR) / ".." / "asset" / "r.svg";
  const auto [width, height, shape] =
      acg::svg_get_image_size_and_shape(input_file_path);
  if (width == 0) { // something went wrong in loading the function
    std::cout << "file open failure" << std::endl;
    abort();
  }
  const std::vector<std::string> outline_path =
      acg::svg_outline_path_from_shape(shape);
  const std::vector<std::vector<acg::Edge>> loops =
      acg::svg_loops_from_outline_path(outline_path);
  //
  std::vector<unsigned char> img_data(width * height,
                                      255); // grayscale image initialized white
  for (unsigned int ih = 0; ih < height; ++ih) {
    for (unsigned int iw = 0; iw < width; ++iw) {
      const auto org = Eigen::Vector2f(iw + 0.5, ih + 0.5); // pixel center
      const auto dir = Eigen::Vector2f(60., 20.);           // search direction
      int count_cross = 0;
      for (const auto &loop :
           loops) { // loop over loop (letter R have internal/external loops)
        for (const auto &edge : loop) { // loop over edge in the loop
          if (edge.is_bezier) { // in case the edge is a quadratic Bézier
            count_cross += number_of_intersection_ray_against_quadratic_bezier(
                org, dir, edge.ps, edge.pc, edge.pe);
          } else { // in case the edge is a line segment
            count_cross += number_of_intersection_ray_against_edge(
                org, dir, edge.ps, edge.pe);
          }
        }
      }
      if (count_cross % 2 == 1) {      // Jordan's curve theory
        img_data[ih * width + iw] = 0; // paint black if it is inside
      }
    }
  }
  const auto output_file_path =
      std::filesystem::path(PROJECT_SOURCE_DIR) / "output.png";
  stbi_write_png(output_file_path.string().c_str(), width, height, 1,
                 img_data.data(), width);
}
