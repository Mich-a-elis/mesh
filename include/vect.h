#pragma once

#include <iostream>
#include <cassert>
#include <limits>
#include <vector>
#include <complex>

namespace mesh
{
   class Vect
   {
   public:
      Vect(){};
      Vect(double x, double y, double z); // construct w/ specified coordinates

      double &operator[](const int &index);             // returns reference to specified coordinate (0-based indexing: x, y, z)
      const double &operator[](const int &index) const; // returns const reference to specified coordinate (0-based indexing: x, y, z)
      Vect operator+(const Vect &v) const;              // vector addition
      Vect operator-(const Vect &v) const;              // vector subtraction
      Vect operator-(void) const;                       // negation
      double operator*(const Vect &v) const;            // dot product
      Vect operator^(const Vect &v) const;              // cross product
      Vect operator*(const double &c) const;            // scalar product
      Vect operator/(const double &c) const;            // scalar division
      void operator+=(const Vect &v);                   // vector addition / assignment
      void operator-=(const Vect &v);                   // vector subtraction / assignment
      void operator*=(const double &c);                 // scalar product / assignment
      void operator/=(const double &c);                 // scalar division / assignment
      double norm(void) const;                          // returns Euclidean length
      Vect unit(void) const;                            // returns vector divided by norm
      void normalize(void);                             // divides by norm
      Vect proj(const Vect &unitDir) const;
      Vect rotate(const Vect &axis, double theta) const;

      double x = std::numeric_limits<double>::quiet_NaN();
      double y = std::numeric_limits<double>::quiet_NaN();
      double z = std::numeric_limits<double>::quiet_NaN(); // coordinates
   };

   Vect operator*(const double &c, const Vect &v);            // scalar product
   std::ostream &operator<<(std::ostream &os, const Vect &o); // prints coordinates

   double angle_between(const std::complex<double> &vec_a, const std::complex<double> &vec_b);
   double angle_between_vect(const Vect &vec_k, Vect vec_a, Vect vec_b);
   double angle_between_cross(const Vect &vec_k, Vect vec_a, Vect vec_b);
   std::pair<bool, Vect> in_triangle(std::vector<Vect> coords, const Vect &pt, Vect normal);
   Vect intersection_seg(const Vect &A, const Vect &B, const Vect &C, const Vect &D, const Vect &oZ, bool strict = false);

   double vector_norm(const std::vector<double> &vect);
   std::vector<double> operator*(const std::vector<std::vector<double>> &A, const std::vector<double> &B);
   std::vector<std::vector<double>> operator/(const std::vector<std::vector<double>> &A, const double &val);
   double determinant(const std::vector<std::vector<double>> &A);
   std::vector<std::vector<double>> inverse(const std::vector<std::vector<double>> &A);

   template <typename T>
   std::pair<bool, int> index_of(const std::vector<T> &vec, const T &elem)
   {
      std::pair<bool, int> result(false, -1);
      auto it = std::find(vec.begin(), vec.end(), elem);
      if (it != vec.end())
      {
         result.first = true;
         result.second = std::distance(vec.begin(), it);
      }
      return result;
   }

} // namespace mesh
