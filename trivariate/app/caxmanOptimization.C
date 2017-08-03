#include <iostream>
#include <fstream>
#include <cmath>
#include <string>
#include <vector>

#include <GoTools/utils/Point.h>
#include <GoTools/trivariate/Hexahedron.h>
#include <GoTools/trivariate/SplineVolume.h>

using namespace std;
using namespace Go;

int main (int argc, char* argv[]) {

  if (argc != 2) {
    std::cout << "Usage: ./caxmanOptimization inputparametervalues" << std::endl;
    return -1;
  }
  
  ifstream ifs(argv[1]);
  
  string time;
  int numpts, dimpts;
  ifs >> time >> numpts >> dimpts;

  std::string dummy;
  double alphadeg;                  
  double alpha;
  double delta;

  Point c(0.0,0.0,0.0);
  Point u(10.0,0.0,0.0);
  Point v(0.0,8.0/sqrt(3.0)+1,0.0);
  Point h(0.0,0.0,4.0);
  
  {
	// Compute the nominal geometry
    Point hn(0.0,0.0,3.9);
	  
	alphadeg = 60;
	alpha = alphadeg/360.0*2*M_PI;
    
	c = Point(0.05,0.05,0.0);
	u = Point(9.9,0.0,0.0);
    v = Point(0.0,8.0/sqrt(3.0)+0.9,0.0);
    h = Point(0.0,0.0,3.9);
    
	// Compute quantities that define the geometry
    double mu = 3.9/sqrt(3.0);
    Point vw = mu*v/v.length();
    Point w = c+h+vw;

    // Define points of the tetrahedron
    Point p1 = c;
    Point p2 = c+u;
    Point p3 = c+u+v;
    Point p4 = c+v;
    Point p5 = p1+h+vw;
    Point p6 = p2+h+vw;
    Point p7 = p3+h-vw;
    Point p8 = p4+h-vw;

    Hexahedron hex(p1,p2,p3,p4,p5,p6,p7,p8);
    //cout << (p1-p2).length() << "\n" << (p2-p3).length() << "\n" << (p5-p6).length() << "\n"  << (p6-p7).length() << endl;

    // Write as Hexahedron
    //hex.write(cout);
    //cout << "alpha = " << (acos((w*v)/w.length()/v.length())/(2*M_PI)*360) << endl;
    //cout << "delta = " << ((p5-p6).length()-10.0)/2.0 << endl;
    //cout << "width = " << (p5-p8).length() << endl;
    //cout << endl;
	  
	shared_ptr<SplineVolume> splvol(hex.geometryVolume());
    vector<shared_ptr<SplineSurface> > bd_sfs_spl = splvol->getBoundarySurfaces();
    for (size_t ki=0; ki<bd_sfs_spl.size(); ++ki)
      {
        bd_sfs_spl[ki]->writeStandardHeader(cout);
        bd_sfs_spl[ki]->write(cout);
		//Point normal;
		//bd_sfs_spl[ki]->normal(normal,0.5,0.5);
		//cout << normal << endl;
      }
  }

  // Check write as spline volumes
  //ofstream ofs("~/caxman_volumes.g2");
  for (int ix=0; ix!=numpts; ix++) {
    // Read parameters from file
    ifs >> dummy >> alphadeg >> delta;
    // Convert to radians
    double alpha = alphadeg/360.0*2*M_PI;
    // Compute quantities that define the geometry
    double mu = sqrt((16.0+pow(delta,2))/(1.0/pow(cos(alpha),2)-1.0));
    Point uw = delta*u/u.length();
    Point vw = mu*v/v.length();
    Point w = c+h+vw-uw;

    // Define points of the tetrahedron
    Point p1 = c;
    Point p2 = c+u;
    Point p3 = c+u+v;
    Point p4 = c+v;
    Point p5 = p1+h+vw-uw;
    Point p6 = p2+h+vw+uw;
    Point p7 = p3+h-vw+uw;
    Point p8 = p4+h-vw-uw;

    Hexahedron hex(p1,p2,p3,p4,p5,p6,p7,p8);
  
    //cout << (p1-p2).length() << "\n" << (p2-p3).length() << "\n" << (p5-p6).length() << "\n"  << (p6-p7).length() << endl;

    // Write as Hexahedron
    hex.write(cout);
    
    // Check parametrization is correct
    //for (int jx3=0;jx3!=2;++jx3) {
    //  for (int jx2=0;jx2!=2;++jx2) {
    //    for (int jx1=0;jx1!=2;++jx1) {
    //      Point pt;
    //      hex.point(pt,jx1,jx2,jx3);
    //      cout << pt << endl;
    //    }
    //  }
    //}

    // Test writing as SplineVolume
    //SplineVolume* sv = hex.geometryVolume();
    //ofs << "700 0 0 0\n";  
    //(*sv).write(ofs);
    
    // Check that the input values alpha, delta are correctly reflected in the geometry
    cout << "alpha = " << (acos((w*v)/w.length()/v.length())/(2*M_PI)*360) << endl;
    cout << "delta = " << ((p5-p6).length()-10.0)/2.0 << endl;
    cout << "width = " << (p5-p8).length() << "\n\n" << endl;

  }
  
}
