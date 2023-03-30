/*!
 * @file
 * src for AddOffset class.
 * */
#include "myAddOffset.hpp"
#include <cmath>
#include <iostream>
#include <helper.hpp>

bool csaransh::myAddOffset::_isUnitcell(double x, double y, double z, double l,
                                      std::array<double, 3> origin) {
  double pos[3] = {x, y, z};
  using csaransh::invars::epsilon;
  for (int i = 0; i < 3; i++) {
    if (!(pos[i] >= (origin[i] * l - epsilon) &&
          pos[i] <= (origin[i] * l + l + epsilon)))
      return false;
  }
  return true;
}

/*
 * The bcc lattice can be seen as formed of two simple lattices interwoven.
 * The only positions we want to know are origins of these two simple lattices.
 * First one is at the origin of the bcc lattice and the other one is the body
 * centered atom of the first unitcell.
 * */
void csaransh::myAddOffset::_bccUnitcell() {
  _sites.clear();
  _sites.emplace_back(std::array<long double, 3>{{0.0, 0.0, 0.0}});
  _sites.emplace_back(std::array<long double, 3>{{0.0, 0.0, 1.0}});
  _sites.emplace_back(std::array<long double, 3>{{0.0, 1.0, 0.0}});
  _sites.emplace_back(std::array<long double, 3>{{0.0, 1.0, 1.0}});  
  _sites.emplace_back(std::array<long double, 3>{{1.0, 0.0, 0.0}});
  _sites.emplace_back(std::array<long double, 3>{{1.0, 0.0, 1.0}});
  _sites.emplace_back(std::array<long double, 3>{{1.0, 1.0, 0.0}});
  _sites.emplace_back(std::array<long double, 3>{{1.0, 1.0, 1.0}});
//   for (auto &it : _sites) {
//     for (auto &jt : it) {
//       jt *= _latConst;
//     }
//   }
}

/*
 * The fcc lattice can be seen as formed of four simple lattices interwoven.
 * The only positions we want to know are origins of these four simple lattices.
 * Since, here we are not computing these four coordinates first we can't return
 * after first four like we do in bcc. The complete unitcell is thus added.
 * TODO: Test thoroughly
 * */
void csaransh::myAddOffset::_fccUnitcell() {
  _sites.clear();
  _sites.emplace_back(std::array<long double, 3>{{0.0, 0.0, 0.0}});
  _sites.emplace_back(std::array<long double, 3>{{0.5, 0.5, 0.0}});
  _sites.emplace_back(std::array<long double, 3>{{0.5, 0.0, 0.5}});
  _sites.emplace_back(std::array<long double, 3>{{0.0, 0.5, 0.5}});
  for (auto &it : _sites) {
    for (auto &jt : it) {
      jt *= _latConst;
    }
  }
}

/**
 * calculate shorter of the distance between two points considering a periodic
 * boundary of a given size. The distance shorter of the two distances (direct
 * and mirror / through periodic boundary) is considered and information about
 * this is stored in the output parameter mirror.
 **/
double csaransh::myAddOffset::_calcDist(std::array<long double, 3> a,
                                            std::array<long double, 3> b) {
  long double res = 0;
  for (int i = 0; i < 3; i++) {
    long double dist = a[i] - b[i];
    res += dist * dist;
  }
  return std::sqrt(res);
}

csaransh::myAddOffset::  myAddOffset(double latConst,int ncell, std::string lattice, std::array<double, 3> origin,
  std::pair<std::vector<std::array<double, 3>>,std::vector<std::array<double, 3>>> & originAtoms){

  _latConst = latConst;
  _ncell = ncell;
  if (lattice[0] == 'b') {
    _bccUnitcell();
  } else if (lattice[0] == 'f') {
    _fccUnitcell();
  } else {
    // TODO log error
  }

  _originVertex = originAtoms.first;
  _originCenter = originAtoms.second;
  _origin = std::array<long double, 3>{{origin[0], origin[1], origin[2]}};
  _base = std::array<long double, 3>{{_originVertex[0][0], _originVertex[0][1], _originVertex[0][2]}};
}

template <class T> void print(std::array<T, 3> x) {
  for (auto jt : x) {
    std::cout << jt << ", ";
  }
  std::cout << '\n';
}


std::tuple<std::array<double, 3>, double, std::array<double, 3>>
csaransh::myAddOffset::operator()(const std::array<double, 3> &c) {
  std::array<long double, 3> coords{{c[0], c[1], c[2]}};
  std::array<int, 3> divCoords;
  std::array<double, 3> cellPos;
  for (int i = 0; i < 3; i++) {
    long double orig = coords[i] - (_base[i]);
    divCoords[i] = int(orig / _latConst);
  }
//   std::cout<<divCoords[0]<<" "<<divCoords[1]<<" "<<divCoords[2]<<"\n";
  auto min = -1.;
  // std::cout<<"working with coord: "<<coords[0]<<" "<<coords[1]<<" "<<coords[2]<<"\n";
  for (auto it : _sites) {
    std::array<int, 3> div{{divCoords[0]+int(it[0]),divCoords[1]+int(it[1]),divCoords[2]+int(it[2])}};
    if(div[0]<0 || div[0]>=_ncell || div[1]<0 || div[1]>=_ncell || div[2]<0 || div[2]>=_ncell) continue;
    int offset = div[0] + div[1] * _ncell + div[2] * _ncell * _ncell;
    std::array<long double, 3> vertex{{_originVertex[offset][0],_originVertex[offset][1],_originVertex[offset][2]}};
    auto dist = _calcDist(coords, vertex);
    if (dist < min || min == -1) {
      min = dist;
      for (int i = 0; i < 3; i++) {
        cellPos[i] = floorf((div[i]) * roundOffTo) / roundOffTo;
      }
      // std::cout<<"v: "<<dist<<"  "<<coords[0]<<" "<<coords[1]<<" "<<coords[2]<<" "<<vertex[0]<<" "<<vertex[1]<<" "<<vertex[2]<<"\n";
    }

    std::array<long double, 3> center{{_originCenter[offset][0],_originCenter[offset][1],_originCenter[offset][2]}};
    dist = _calcDist(coords, center);
    if(dist < min || min == -1){
        min = dist;
        for (int i = 0; i < 3; i++) {
            cellPos[i] = floorf((div[i] + 0.5) * roundOffTo) / roundOffTo;
        }
        // std::cout<<"c: "<<dist<<"  "<<coords[0]<<" "<<coords[1]<<" "<<coords[2]<<" "<<center[0]<<" "<<center[1]<<" "<<center[2]<<"\n";
    }
    // std::cout<<"-----------------compare with coord: "<<vertex[0]<<" "<<vertex[1]<<" "<<vertex[2]<<"\n";
    // std::cout<<"-----------------compare with coord: "<<center[0]<<" "<<center[1]<<" "<<center[2]<<"\n";
  }
  // std::cout<<coords[0]<<" "<<coords[1]<<" "<<coords[2]<<" "<<cellPos[0]*_latConst<<" "<<cellPos[1]*_latConst<<" "<<cellPos[2]*_latConst<<"  min:"<<min<<"\n";
  // if(min>0){
  //   for (auto it : _sites) {
  //     std::array<int, 3> div{{divCoords[0]+int(it[0]),divCoords[1]+int(it[1]),divCoords[2]+int(it[2])}};
  //     int offset = div[0] + div[1] * _ncell + div[2] * _ncell * _ncell;
  //     std::array<long double, 3> vertex{{_originVertex[offset][0],_originVertex[offset][1],_originVertex[offset][2]}};
  //     std::array<long double, 3> center{{_originCenter[offset][0],_originCenter[offset][1],_originCenter[offset][2]}};
      
  //     std::cout<<"-----------------compare with coord: "<<vertex[0]<<" "<<vertex[1]<<" "<<vertex[2]<<"\n";
  //     std::cout<<"-----------------compare with coord: "<<center[0]<<" "<<center[1]<<" "<<center[2]<<"\n";
  //   }
  // }
  return std::make_tuple(cellPos, min, c);
}