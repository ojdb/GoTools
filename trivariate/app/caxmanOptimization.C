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

// Only handles planar hexahedra at the mo
void writeToPLY(ofstream& ofs, 
		const Hexahedron& hex,
		const vector<string>& comments=vector<string>() 
		//const int resx=1, // do we need resolution for hexahedra? 
		//const int resy=1, 
		//const int resz=1
		) { 
  ofs << "ply\n"
      << "format ascii 1.0\n";
  
  for (size_t ix=0; ix!=comments.size(); ++ix) 
      ofs << "comment " << comments[ix] << "\n";
  
  ofs << "element vertex 8\n"
      // << (resx+1)*(resy+1)*(resz+1) << "\n"
      << "property float64 x\n"
      << "property float64 y\n"
      << "property float64 z\n"
      << "element face 12\n"
      // 2*(resx*resy+resx*resz+resy*resz) << "\n"
      << "property list uchar int vertex_indices\n"
      << "end_header\n";

  Point pt;
  for (int ix=0; ix!=1+1; ++ix) {
    for (int jx=0; jx!=1+1; ++jx) {
      for (int kx=0; kx!=1+1; ++kx) {
	hex.point(pt,double(ix),double(jx),double(kx));
	ofs << pt << "\n";	
      }
    }    
  } 

  ofs << "3 0 1 2\n3 1 3 2\n3 6 7 4\n3 7 5 4\n"
      << "3 4 5 0\n3 5 1 0\n3 2 3 6\n3 3 7 6\n"
      << "3 6 4 2\n3 4 0 2\n3 3 1 7\n3 1 5 7" << endl;
}

int main (int argc, char* argv[]) {

  if (argc < 2 || argc > 4) {
    std::cout << "Usage: ./caxmanOptimization input_ascii_grid [output_prefix (optional)] [output_prefix_nominal (optional)]" << std::endl;
    return -1;
  }
  
  ifstream ifs(argv[1]);
  if (!ifs.good()) { 
    cerr << "Error: Problem opening parameter file at: " << argv[1] << ", exiting...\n";
    return -1;
  }
  string prefix;
  string prefixnom;
  if (argc > 2) prefix = string(argv[2]); 
  if (argc > 3) prefixnom = string(argv[3]);

  string time;
  int numpts, dimpts;
  ifs >> time >> numpts >> dimpts;

  std::string dummy;
  double alphadeg;                  
  double alpha;
  double delta;

  Point c(0.05,0.05,0.0);
  Point u(9.9,0.0,0.0);
  Point v(0.0,8.0/sqrt(3.0)+0.9,0.0);
  Point h(0.0,0.0,3.9);
  
  // Compute the nominal geometry
  alphadeg = 60;
  alpha = alphadeg/360.0*2*M_PI;

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

  // Define the hexahedron
  Hexahedron nominal_hex(p1,p2,p3,p4,p5,p6,p7,p8);

  // Write as Hexahedron
  //hex.write(cout);
  cout << "Nominal geometry" << endl;
  cout << "alpha = " << (acos(((w-c)*v)/(w-c).length()/v.length())/(2*M_PI)*360) << endl;
  //cout << "delta = " << ((p5-p6).length()-10.0)/2.0 << endl;
  cout << "width = " << (p5-p8).length() << endl;
  cout << endl;

  // Define points of inflated geometry, on which the parametrized geometries are based 
  c = Point(0.0,0.0,0.0);
  u = Point(10,0.0,0.0);
  v = Point(0.0,8.0/sqrt(3.0)+1.0,0.0);
  h = Point(0.0,0.0,4.0);

  // Extract the spline surfaces
  shared_ptr<SplineVolume> splvol(nominal_hex.geometryVolume());
  vector<shared_ptr<SplineSurface> > bd_sfs_spl = splvol->getBoundarySurfaces();
  // Hack for incorrect orientation  
  bd_sfs_spl[0]->turnOrientation();
  bd_sfs_spl[3]->turnOrientation();
  bd_sfs_spl[4]->turnOrientation();
  ofstream ofsg2(prefixnom+string("nominal_geometry.g2"));
  for (size_t ki=0; ki<bd_sfs_spl.size(); ++ki)
    {
      bd_sfs_spl[ki]->writeStandardHeader(ofsg2);
      bd_sfs_spl[ki]->write(ofsg2);
      //Point normal;
      //bd_sfs_spl[ki]->normal(normal,0.5,0.5);
      //cout << normal << endl;
  }

  // Write tessellation to file
  ofstream ofs(prefixnom+string("nominal_geometry.ply"));
  writeToPLY(ofs,nominal_hex,vector<string>(1,string("N1")));

  // Loop through the parameter values and generate geometry for each
  for (int ix=0; ix!=numpts; ix++) {
    // Read parameters from file
    ifs >> dummy >> alphadeg >> delta;
    dummy = dummy.substr(0,dummy.size()-1);
    
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
    //hex.write(cout);
    
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
    cout << dummy << endl;
    cout << "alpha = " << (acos((w*v)/w.length()/v.length())/(2*M_PI)*360) << endl;
    cout << "delta = " << ((p5-p6).length()-10.0)/2.0 << endl;
    cout << "width = " << (p5-p8).length() << endl;

    // Write tessellation to file
    
    string outtmp = prefix+string("param_geom_")+dummy+string(".ply");
    cout << "Tessellated geometry written to: " << outtmp << "\n\n" << endl;
    ofstream ofstmp(outtmp);
    writeToPLY(ofstmp,hex,vector<string>(1,dummy));
    ofstmp.close(); 

  }
  
}
