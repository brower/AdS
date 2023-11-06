// ising_rect_rad.cc

//#include <getopt.h>
//#include <cassert>
//#include <cmath>
//#include <cstdio>
#include <string>
#include <vector>
//#include <Eigen/Dense>
#include "ising-GTF.h"
//#include "statistics.h"

#include <iostream>

int main(int argc, char* argv[]) {

  std::vector<int> L = {4, 4, 4};
  std::vector<double> sc = {1., 1., 1.} ;
  // std::vector<double> sc = {0., 0., 0.} ;
  std::vector<double> fcc = {0., 0., 0., 0., 0., 0.} ;
  // std::vector<double> fcc = {1., 1., 1., 1., 1., 1.} ;
  std::vector<double> bcc = {0., 0., 0., 0.} ;
  // std::vector<double> bcc = {1., 1., 1., 1.} ;

  // initialize the lattice
  QfeLattice lattice;
  lattice.InitCubic(L, sc, fcc, bcc) ;
  // lattice.InitCubic(L, sc, fcc, bcc, true) ;

#if 0
  for (int s = 0; s < lattice.n_sites; s++) {
    int x = s % L[0] ;
    int y = (s / L[0]) % L[1] ;
    int z = (s / L[0]) / L[1] ;

    std::cout << "site " << s << " x " << x << " y " << y << " z " << z 
      << " fcc " << (x+y+z)%2 << " bcc " << ((x % 2 == y % 2) && (y % 2 == z % 2))
      << " nn " << lattice.sites[s].nn << std::endl ;
  }
#endif

#if 1
  // This will print out nice Mathematica graphics directives.
  // Will not include wrap-around links.

  std::cout << "verts = {" << std::endl ;

  for (int s = 0; s < lattice.n_sites; s++) {
    int x = s % L[0] ;
    int y = (s / L[0]) % L[1] ;
    int z = (s / L[0]) / L[1] ;

    std::string color, shape, size ;

    if ((x+y+z)%2 == 0) {
      color = "Red" ;
    } else {
      color = "Black" ;
    }

    if ((x%2 == y%2) && (y%2 == z%2)) {
      shape = "Cube" ;
      size  = "0.2" ;
    } else {
      shape = "Ball" ;
      size  = "0.1" ;
    }

    if ( s < lattice.n_sites - 1) {
      std::cout << "{" << color << ", " << shape << "[{"
        << x << "," << y << "," << z << "}, " << size << "]}," << std::endl ;
    } else {
      std::cout << "{" << color << ", " << shape << "[{"
        << x << "," << y << "," << z << "}, " << size << "]}" << std::endl << "}" << std::endl ;
    }

  } // end for (int s = 0; s < lattice.n_sites; s++)

  std::cout << "edges = {" << std::endl ;

  for (int l = 0; l < lattice.n_links; l++) {

    std::string color ;

    int ax = lattice.links[l].sites[0] % L[0] ;
    int ay = (lattice.links[l].sites[0] / L[0]) % L[1] ;
    int az = (lattice.links[l].sites[0] / L[0]) / L[1] ;
    int bx = lattice.links[l].sites[1] % L[0] ;
    int by = (lattice.links[l].sites[1] / L[0]) % L[1] ;
    int bz = (lattice.links[l].sites[1] / L[0]) / L[1] ;
    int dx = bx - ax ;
    int dy = by - ay ;
    int dz = bz - az ;

    int d2 = dx*dx + dy*dy + dz*dz ;

    switch ( d2 ) {
    case 1 :
      color = "Gray" ; break ;
    case 2 :
      color = "Magenta" ; break ;
    case 3 :
      color = "Cyan" ; break ;
    default :
      color = "Transparent" ;

    }



    if (l < lattice.n_links - 1) {
      std::cout << "{Opacity[" << lattice.links[l].wt << ", " << color << "], Thick, Line[{{"
        << ax << "," << ay << "," << az << "},{"
        << bx << "," << by << "," << bz << "}}]}," << std::endl ;
    } else {
      std::cout << "{Opacity[" << lattice.links[l].wt << ", " << color << "], Thick, Line[{{"
        << ax << "," << ay << "," << az << "},{"
        << bx << "," << by << "," << bz << "}}]}" << std::endl << "}" << std::endl ;
    }

  } // end for (int l = 0; l < lattice.n_links; l++)

#endif

  return 0;
}
