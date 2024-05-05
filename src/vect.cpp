#define _USE_MATH_DEFINES

#include <complex>
#include <iostream>
#include <math.h>
#include <vector>

#include "vect.h"

using namespace std;

namespace mesh
{
   Vect::Vect(double x0, double y0, double z0)
       : x(x0), y(y0), z(z0)
   {
   }

   double &Vect::operator[](const int &index)
   {
      if (index < 0 || index >= 3)
         throw std::out_of_range("Invalid index for Vec3");
      return (&x)[index];
   }

   const double &Vect::operator[](const int &index) const
   {
      if (index < 0 || index >= 3)
         throw std::out_of_range("Invalid index for Vec3");
      return (&x)[index];
   }

   Vect Vect::operator+(const Vect &v) const
   {
      return Vect(x + v.x, y + v.y, z + v.z);
   }

   Vect Vect::operator-(const Vect &v) const
   {
      return Vect(x - v.x, y - v.y, z - v.z);
   }

   Vect Vect::operator-(void) const
   {
      return Vect(-x, -y, -z);
   }

   double Vect::operator*(const Vect &v) const
   {
      return x * v.x + y * v.y + z * v.z;
   }

   Vect Vect::operator^(const Vect &v) const
   {
      return Vect(y * v.z - z * v.y, z * v.x - x * v.z, x * v.y - y * v.x);
   }

   Vect Vect::operator*(const double &c) const
   {
      return Vect(x * c, y * c, z * c);
   }

   Vect operator*(const double &c, const Vect &v)
   {
      return v * c;
   }

   Vect Vect::operator/(const double &c) const
   {
      return (*this) * (1. / c);
   }

   void Vect::operator+=(const Vect &v)
   {
      x += v.x;
      y += v.y;
      z += v.z;
   }

   void Vect::operator-=(const Vect &v)
   {
      x -= v.x;
      y -= v.y;
      z -= v.z;
   }

   void Vect::operator*=(const double &c)
   {
      x *= c;
      y *= c;
      z *= c;
   }

   void Vect::operator/=(const double &c)
   {
      (*this) *= (1. / c);
   }

   double Vect::norm(void) const
   {
      return sqrtf((*this) * (*this));
   }

   void Vect::normalize(void)
   {
      (*this) /= norm();
   }

   Vect Vect::unit(void) const
   {
      return (*this) / norm();
   }

   std::ostream &operator<<(std::ostream &os, const Vect &o)
   {
      os << "[" << o.x << ", " << o.y << ", " << o.z << "]";
      return os;
   }

   Vect Vect::proj(const Vect &normal) const
   {
      Vect normalU = normal.unit();
      return *this - normalU * (*this * normalU);
   }

   Vect Vect::rotate(const Vect &axis, double ang) const
   {

      double cosAngle = cos(ang);
      double sinAngle = sin(ang);
      double oneMinusCos = 1 - cosAngle;

      double axisX = axis.x;
      double axisY = axis.y;
      double axisZ = axis.z;

      double rotatedX = (cosAngle + (axisX * axisX) * oneMinusCos) * x +
                        (axisX * axisY * oneMinusCos - axisZ * sinAngle) * y +
                        (axisX * axisZ * oneMinusCos + axisY * sinAngle) * z;

      double rotatedY = (axisY * axisX * oneMinusCos + axisZ * sinAngle) * x +
                        (cosAngle + (axisY * axisY) * oneMinusCos) * y +
                        (axisY * axisZ * oneMinusCos - axisX * sinAngle) * z;

      double rotatedZ = (axisZ * axisX * oneMinusCos - axisY * sinAngle) * x +
                        (axisZ * axisY * oneMinusCos + axisX * sinAngle) * y +
                        (cosAngle + (axisZ * axisZ) * oneMinusCos) * z;

      return Vect(rotatedX, rotatedY, rotatedZ);
   }

   double angle_between(const std::complex<double> &vec_a, const std::complex<double> &vec_b)
   {
      assert(std::abs(vec_a) > 1e-14);
      assert(std::abs(vec_b) > 1e-14);
      complex<double> vec_j(-vec_a.imag(), vec_a.real());
      double c = vec_a.real() * vec_b.real() + vec_a.imag() * vec_b.imag();
      double s = vec_j.real() * vec_b.real() + vec_j.imag() * vec_b.imag();
      return atan2(s, c);
   }

   double angle_between_vect(const Vect &vec_k, Vect vec_a, Vect vec_b)
   {
      assert(vec_k.norm() > 1e-14);
      vec_a = vec_a.proj(vec_k);
      vec_b = vec_b.proj(vec_k);
      // cout << "vec_a.norm() " << vec_a.norm() << endl;
      // cout << "vec_b.norm() " << vec_b.norm() << endl;
      assert(vec_a.norm() > 1e-14);
      assert(vec_b.norm() > 1e-14);
      Vect vec_j = vec_a.rotate(vec_k, M_PI_2);
      double c = vec_a * vec_b;
      double s = vec_j * vec_b;
      return atan2(s, c);
   }

   double angle_between_cross(const Vect &vec_k, Vect vec_a, Vect vec_b)
   {
      double temp_ang = angle_between_vect(vec_k, vec_a, vec_b);
      double ang = temp_ang;
      for (int k = 1; k < 4; ++k)
      {
         double temp = atan2(sin(temp_ang + k * M_PI_2), cos(temp_ang + k * M_PI_2));
         if (std::abs(temp) < std::abs(ang))
            ang = temp;
      }
      return ang;
   }

   std::pair<bool, Vect> in_triangle(std::vector<Vect> coords, const Vect &pt, Vect normal)
   {
      if (std::isnan(normal.norm()) == true)
      {
         normal = (coords[1] - coords[0]) ^ (coords[2] - coords[0]);
         normal.normalize();
      }

      pair<bool, Vect> res = pair<bool, Vect>(false, Vect());
      double dist_pt_plan = abs(normal * (pt - coords[0]));

      // cout << "dist_pt_plan " << dist_pt_plan << endl;
      if (!(dist_pt_plan > 0 + 1e-8))
      {
         double D = normal * ((coords[1] - coords[0]) ^ (coords[2] - coords[0]));
         // cout << "D " << D << endl;
         assert(D > 0 + 1e-15);

         Vect lam;
         lam[0] = (normal * ((coords[1] - pt) ^ (coords[2] - pt))) / D;
         lam[1] = (normal * ((coords[2] - pt) ^ (coords[0] - pt))) / D;
         lam[2] = (normal * ((coords[0] - pt) ^ (coords[1] - pt))) / D;

         // cout << "lam " << lam << endl;
         vector<double> lam_vector = vector<double>({lam[0], lam[1], lam[2]});
         auto minElement = std::min_element(lam_vector.begin(), lam_vector.end());
         auto maxElement = std::max_element(lam_vector.begin(), lam_vector.end());
         // cout << "minElement " << *minElement << " maxElement " << *maxElement << endl;

         if (*minElement > 0 - 1e-8 && *maxElement < 1 + 1e-8)
         {
            res.first = true;
            res.second = lam;
         }
      }
      return res;
   }

   Vect intersection_seg(const Vect &A, const Vect &B, const Vect &C, const Vect &D, const Vect &oZ, bool strict)
   {
      Vect AB = B - A;
      Vect CD = D - C;
      assert(AB.norm() > 1e-15);
      assert(CD.norm() > 1e-15);
      // cout << "std::abs(oZ * AB) std::abs(oZ * CD) < 1e-6) " << std::abs(oZ.unit() * AB.unit()) << " " << std::abs(oZ.unit() * CD.unit()) << endl;
      assert(std::abs(oZ.unit() * AB.unit()) < 1e-6);
      assert(std::abs(oZ.unit() * CD.unit()) < 1e-6);

      Vect u = AB.rotate(oZ, M_PI_2).unit();
      Vect v = CD.rotate(oZ, M_PI_2).unit();

      Vect res = Vect();
      if (std::abs(AB.unit() * v.unit()) < 1e-12 || std::abs(CD.unit() * u.unit()) < 1e-12)
         assert(std::abs(AB.unit() * v.unit()) < 1e-12 && std::abs(CD.unit() * u.unit()) < 1e-12);
      else
      {
         Vect AC = C - A;
         double t = (AC * v) / (AB * v);
         Vect CA = A - C;
         double k = (CA * u) / (CD * u);
         double eps = 1e-10;
         if (strict == false)
         {
            if (0 - eps < t && t < 1 + eps && 0 - eps < k && k < 1 + eps)
               res = A + t * (B - A);
         }
         else
         {
            if (0 + eps < t && t < 1 - eps && 0 + eps < k && k < 1 - eps)
               res = A + t * (B - A);
         }
      }

      return res;
   }

   /////////////////////////////////////////////////////////////////////////////////////
   double vector_norm(const std::vector<double> &vect)
   {
      double sum = 0;
      for (int k = 0; k < vect.size(); ++k)
         sum += std::pow(vect[k], 2);
      return std::sqrt(sum);
   }

   std::vector<double> operator*(const std::vector<vector<double>> &A, const std::vector<double> &B)
   {
      for (int k = 0; k < A.size(); ++k)
         assert(A[k].size() == B.size());
      vector<double> res(A.size());
      for (int k = 0; k < A.size(); ++k)
      {
         res[k] = 0;
         for (int i = 0; i < B.size(); ++i)
            res[k] += A[k][i] * B[i];
      }
      return res;
   }

   std::vector<std::vector<double>> operator/(const std::vector<vector<double>> &A, const double &val)
   {
      vector<vector<double>> res(A.size());
      for (int k = 0; k < A.size(); ++k)
         res[k].resize(A[k].size());

      for (int i = 0; i < A.size(); ++i)
         for (int j = 0; j < A[i].size(); ++j)
            res[i][j] = A[i][j] / val;
      return res;
   }

   double determinant(const vector<vector<double>> &A)
   {
      for (int k = 0; k < A.size(); ++k)
         assert(A[k].size() == A.size());

      double det = 0;
      switch (A.size())
      {
      case 2:
         det = A[0][0] * A[1][1] - A[1][0] * A[0][1];
         break;
      case 3:
         det = A[0][0] * A[1][1] * A[2][2] + A[0][1] * A[1][2] * A[2][0] + A[0][2] * A[1][0] * A[2][1] -
               A[2][0] * A[1][1] * A[0][2] - A[2][1] * A[0][0] * A[1][2] - A[1][0] * A[2][2] * A[0][1];
         break;
      default:
         assert(true == false);
      }
      return det;
   }

   vector<vector<double>> inverse(const vector<vector<double>> &A)
   {
      for (int k = 0; k < A.size(); ++k)
         assert(A[k].size() == A.size());

      double det = 0;
      vector<vector<double>> inv(A.size());
      for (int k = 0; k < A.size(); ++k)
         inv[k].resize(A.size());

      switch (A.size())
      {
      case 1:
         inv[0][0] = 1;
         break;
      case 2:
         det = A[0][0] * A[1][1] - A[1][0] * A[0][1];
         inv[0][0] = A[1][1];
         inv[0][1] = -A[0][1];
         inv[1][0] = -A[1][0];
         inv[1][1] = A[0][0];
         inv = inv / det;
         break;
      case 3:
         det = A[0][0] * A[1][1] * A[2][2] + A[0][1] * A[1][2] * A[2][0] + A[0][2] * A[1][0] * A[2][1] -
               A[2][0] * A[1][1] * A[0][2] - A[2][1] * A[0][0] * A[1][2] - A[1][0] * A[2][2] * A[0][1];
         inv[0][0] =
             A[1][1] * A[2][2] - A[2][1] * A[1][2];
         inv[0][1] = A[0][2] * A[2][1] - A[2][2] * A[0][1];
         inv[0][2] =
             A[0][1] * A[1][2] - A[1][1] * A[0][2];
         inv[1][0] = A[1][2] * A[2][0] - A[2][2] * A[1][0];
         inv[1][1] =
             A[0][0] * A[2][2] - A[2][0] * A[0][2];
         inv[1][2] = A[0][2] * A[1][0] - A[1][2] * A[0][0];
         inv[2][0] =
             A[1][0] * A[2][1] - A[2][0] * A[1][1];
         inv[2][1] = A[0][1] * A[2][0] - A[2][1] * A[0][0];
         inv[2][2] =
             A[0][0] * A[1][1] - A[1][0] * A[0][1];
         inv = inv / det;
         break;
      default:
         assert(true == false);
      }

      return inv;
   }

} // namespace mesh
