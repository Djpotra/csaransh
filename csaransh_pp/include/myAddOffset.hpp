/*!
 * @file
 * class for calculating offset & closest lattice site of an atom
 * */
#ifndef MYADDOFFSET_CSARANSH_HPP
#define MYADDOFFSET_CSARANSH_HPP

#include <array>
#include <string>
#include <tuple>
#include <vector>

namespace csaransh {

using offsetCoords = std::tuple<std::array<double, 3>, double, std::array<double, 3>>;

class myAddOffset {
public:
  myAddOffset(double latConst,int ncell, std::string lattice, std::array<double, 3> origin,
  std::pair<std::vector<std::array<double, 3>>,std::vector<std::array<double, 3>>> & originAtoms);

  offsetCoords operator()(const std::array<double, 3> &coords);

private:
  bool _isUnitcell(double x, double y, double z, double l,
                   std::array<double, 3> origin);
  void _bccUnitcell();
  void _fccUnitcell();
  double _calcDist(std::array<long double, 3> a,
                         std::array<long double, 3> b);

  std::vector<std::array<long double, 3>> _sites;
  std::array<long double, 3> _origin;
  std::array<long double, 3> _base;
  std::vector<std::array<double, 3>> _originVertex;
  std::vector<std::array<double, 3>> _originCenter;
  long double _latConst;
  int _ncell;
  long double roundOffTo = 10000.0;
};
} // namespace csaransh

#endif // ADDOFFSET_CSARANSH_HPP
