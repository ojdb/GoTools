/*
 * Copyright (C) 1998, 2000-2007, 2010, 2011, 2012, 2013 SINTEF ICT,
 * Applied Mathematics, Norway.
 *
 * Contact information: E-mail: tor.dokken@sintef.no                      
 * SINTEF ICT, Department of Applied Mathematics,                         
 * P.O. Box 124 Blindern,                                                 
 * 0314 Oslo, Norway.                                                     
 *
 * This file is part of GoTools.
 *
 * GoTools is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Affero General Public License as
 * published by the Free Software Foundation, either version 3 of the
 * License, or (at your option) any later version. 
 *
 * GoTools is distributed in the hope that it will be useful,        
 * but WITHOUT ANY WARRANTY; without even the implied warranty of         
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the          
 * GNU Affero General Public License for more details.
 *
 * You should have received a copy of the GNU Affero General Public
 * License along with GoTools. If not, see
 * <http://www.gnu.org/licenses/>.
 *
 * In accordance with Section 7(b) of the GNU Affero General Public
 * License, a covered work must retain the producer line in every data
 * file that is created or manipulated using GoTools.
 *
 * Other Usage
 * You can be released from the requirements of the license by purchasing
 * a commercial license. Buying such a license is mandatory as soon as you
 * develop commercial activities involving the GoTools library without
 * disclosing the source code of your own applications.
 *
 * This file may be used in accordance with the terms contained in a
 * written agreement between you and SINTEF ICT. 
 */

#include "GoTools/geometry/ObjectHeader.h"
#include "GoTools/geometry/PointCloud.h"
#include "GoTools/lrsplines2D/LRSplineSurface.h"
#include "GoTools/lrsplines2D/LRSplineEvalGrid.h"
#include "GoTools/lrsplines2D/LRSplineUtils.h"
#include "GoTools/lrsplines2D/LinDepUtils.h"
#include "GoTools/lrsplines2D/LRSplinePlotUtils.h"

#include "newmat.h"
//#include "newmatio.h"

#include <iostream>
#include <fstream>
#include <vector>
#include <ctime>

using namespace Go;
using namespace std;
//using namespace NEWMAT;

typedef LRSplineSurface::ElementMap::const_iterator elem_it;
typedef std::vector<double>::const_iterator vecit;
typedef std::vector<LRBSpline2D*>::iterator lrb_it;



LRSplineSurface construct_basic_lrspline(const int deg, const int num, bool full_mult) 
{
  const int deg_u     = deg;
  const int deg_v     = deg;
  const int coefs_u   = deg+1+num;
  const int coefs_v   = deg+1+num;
  const int dimension = 1;
  const int num_knots = (deg+1)*2+num;
  double k_u[num_knots];
  double k_v[num_knots];
  double step = 1 / (double) (num+1);//1 / (double) deg;
  for (int ix=0;ix!=num_knots/2;++ix) {
    k_u[ix] = 0.0;
    k_v[ix] = 0.0;
    k_u[num_knots-1-ix] = 1.0;
    k_v[num_knots-1-ix] = 1.0;
  }
  for (int ix=1;ix!=(num+3)/2;++ix) {
    k_u[deg+ix] = 0.0+ix*step;
    k_v[deg+ix] = 0.0+ix*step;
    k_u[num_knots-deg-1-ix] = 1.0-ix*step;
    k_v[num_knots-deg-1-ix] = 1.0-ix*step;
  }

  //for (int ix=0;ix!=num_knots;++ix) {
  //  cout << k_u[ix] << " " << endl;
  //}
  const vector<double> knots_u(k_u,k_u+num_knots);
  const vector<double> knots_v(k_v,k_v+num_knots);
  LRSplineSurface lrs(deg_u, deg_v, 
		      coefs_u, coefs_v, 
		      dimension, 
		      knots_u.begin(), knots_v.begin() );
  return lrs;
  /*
  const int deg_u     = deg;
  const int deg_v     = deg;
  const int bmult_u   = full_mult ? (deg_u+1) : 0;
  const int bmult_v   = full_mult ? (deg_v+1) : 0;
  const int coefs_u   = num+bmult_u;
  const int coefs_v   = num+bmult_v;
  const int dimension = 1;
  const int num_knots_u = coefs_u + bmult_u - 1;
  const int num_knots_v = coefs_v + bmult_v - 1;
  double k_u[num_knots_u];
  double k_v[num_knots_v];
  double step = 0.1/ (double) num;
  for (int ix=0;ix!=num_knots_u/2;++ix) {
    k_u[ix] = 0.01+ix*step;
    k_u[num_knots_u-1-ix] = 0.99-ix*step;
  }
  for (int ix=0;ix!=num_knots_v/2;++ix) {
    k_v[ix] = 0.01+ix*step;
    k_v[num_knots_v-1-ix] = 0.990-ix*step;
  }
  for (int ix=0;ix!=num_knots_v;++ix) {
    cout << k_v[ix] << " " << endl;
    }

  const vector<double> knots_u(k_u,k_u+num_knots_u);
  const vector<double> knots_v(k_v,k_v+num_knots_v);
  LRSplineSurface lrs(deg_u, deg_v, 
		      coefs_u, coefs_v, 
		      dimension, 
		      knots_u.begin(), knots_v.begin());
  cout << "OK" << endl;
  return lrs;*/
}

void set_random_coefficients(LRSplineSurface& lrs) 
{
 srand(time(NULL));
 for ( LRSplineSurface::BSplineMap::const_iterator it = lrs.basisFunctionsBegin(); 
	it != lrs.basisFunctionsEnd(); ++it ) 
    {
      Point pt(1);
      pt[0] = (double) rand() / (double) (RAND_MAX);
      it->second->coefTimesGamma() = pt;
    }
}

double f(const double x, const double y) {
  return x;
  //return y;
  //return 1-x-y;
}

int main(int argc, char *argv[])
{

  // DEG1
  /*int deg = 1, num = 5; 
  const int nrefs = 4;
  double tt = 1.0/6.0;
  // Hierarchical
  LRSplineSurface::Refinement2D rfs[nrefs] = {
    {tt*2.5,tt*2,tt*4,XFIXED,1}, {tt*3.5,tt*2,tt*4,XFIXED,1},
    {tt*2.5,tt*2,tt*4,YFIXED,1}, {tt*3.5,tt*2,tt*4,YFIXED,1}
  };
  // Non-overloaded
  LRSplineSurface::Refinement2D rfs[nrefs] = {
    {tt*2.5,tt*1,tt*5,XFIXED,1}, {tt*3.5,tt*1,tt*5,XFIXED,1},
    {tt*2.5,tt*2,tt*4,YFIXED,1}, {tt*3.5,tt*2,tt*4,YFIXED,1}
    };*/
  /*
  // DEG2
  int deg = 2, num = 8; 
  const int nrefs = 6;
  double tt = 1.0/9.0;
  // Hierarchical
  LRSplineSurface::Refinement2D rfs[nrefs] = {
    {tt*3.5,tt*3,tt*6,XFIXED,1}, {tt*4.5,tt*3,tt*6,XFIXED,1}, {tt*5.5,tt*3,tt*6,XFIXED,1},
    {tt*3.5,tt*3,tt*6,YFIXED,1}, {tt*4.5,tt*3,tt*6,YFIXED,1}, {tt*5.5,tt*3,tt*6,YFIXED,1},
  };
  // Non-overloaded
  LRSplineSurface::Refinement2D rfs[nrefs] = {
    {tt*3.5,tt*1,tt*8,XFIXED,1}, {tt*4.5,tt*3,tt*6,XFIXED,1}, {tt*5.5,tt*1,tt*8,XFIXED,1},
    {tt*3.5,tt*3,tt*6,YFIXED,1}, {tt*4.5,tt*3,tt*6,YFIXED,1}, {tt*5.5,tt*3,tt*6,YFIXED,1},
    };*/

  // DEG3
  /*int deg = 3, num = 11; 
  const int nrefs = 8;
  double tt = 1.0/12.0;
  // Hierarchical
  LRSplineSurface::Refinement2D rfs[nrefs] = {
    {tt*4.5,tt*4,tt*8,XFIXED,1}, {tt*5.5,tt*4,tt*8,XFIXED,1}, {tt*6.5,tt*4,tt*8,XFIXED,1}, {tt*7.5,tt*4,tt*8,XFIXED,1},    
    {tt*4.5,tt*4,tt*8,YFIXED,1}, {tt*5.5,tt*4,tt*8,YFIXED,1}, {tt*6.5,tt*4,tt*8,YFIXED,1}, {tt*7.5,tt*4,tt*8,YFIXED,1},
  };
  // Non-overloaded
   LRSplineSurface::Refinement2D rfs[nrefs] = {
    {tt*4.5,tt*1,tt*11,XFIXED,1}, {tt*5.5,tt*3,tt*9,XFIXED,1}, {tt*6.5,tt*3,tt*9,XFIXED,1}, {tt*7.5,tt*1,tt*11,XFIXED,1},
     {tt*4.5,tt*1,tt*11,YFIXED,1}, {tt*5.5,tt*3,tt*9,YFIXED,1}, {tt*6.5,tt*3,tt*9,YFIXED,1}, {tt*7.5,tt*1,tt*11,YFIXED,1},
  };*/

  /*
  // PLAYGROUND
  int deg = 1, num = 9; 
  const int nrefs = 12;
  double tt = 1.0/(num+1.0);
  LRSplineSurface::Refinement2D rfs[nrefs] = {
    {0.25,0.2,0.4,XFIXED,1}, {0.35,0.2,0.4,XFIXED,1},
    {0.25,0.2,0.4,YFIXED,1}, {0.35,0.2,0.4,YFIXED,1},
    {0.275,0.25,0.35,XFIXED,1}, {0.325,0.25,0.35,XFIXED,1},
    {0.275,0.25,0.35,YFIXED,1}, {0.325,0.25,0.35,YFIXED,1},
    {0.2875,0.275,0.325,XFIXED,1}, {0.3125,0.275,0.325,XFIXED,1},
    {0.2875,0.275,0.325,YFIXED,1}, {0.3125,0.275,0.325,YFIXED,1},
  };*/

  /*
  // CROSS1
  int deg = 1, num = 9; 
  const int nrefs = 16;
  double tt = 1.0/(num+1.0);
  LRSplineSurface::Refinement2D rfs[nrefs] = {
    {0.15,0.0,0.4,XFIXED,1}, {0.15,0.0,0.4,YFIXED,1},
    {0.25,0.1,0.5,XFIXED,1}, {0.25,0.1,0.5,YFIXED,1},
    {0.35,0.2,0.6,XFIXED,1}, {0.35,0.2,0.6,YFIXED,1},
    {0.45,0.3,0.7,XFIXED,1}, {0.45,0.3,0.7,YFIXED,1},
    {0.55,0.4,0.8,XFIXED,1}, {0.55,0.4,0.8,YFIXED,1},
    {0.65,0.5,0.9,XFIXED,1}, {0.65,0.5,0.9,YFIXED,1},
    {0.75,0.6,1.0,XFIXED,1}, {0.75,0.6,1.0,YFIXED,1},
    {0.85,0.7,1.0,XFIXED,1}, {0.85,0.7,1.0,YFIXED,1}, 
  };
  */

  /*
  // CROSS2
  int deg = 2, num = 9; 
  const int nrefs = 16;
  double tt = 1.0/(num+1.0);
  LRSplineSurface::Refinement2D rfs[nrefs] = {
    {0.15,0.0,0.4,XFIXED,1}, {0.15,0.0,0.4,YFIXED,1},
    {0.25,0.1,0.5,XFIXED,1}, {0.25,0.1,0.5,YFIXED,1},
    {0.35,0.2,0.6,XFIXED,1}, {0.35,0.2,0.6,YFIXED,1},
    {0.45,0.3,0.7,XFIXED,1}, {0.45,0.3,0.7,YFIXED,1},
    {0.55,0.4,0.8,XFIXED,1}, {0.55,0.4,0.8,YFIXED,1},
    {0.65,0.5,0.9,XFIXED,1}, {0.65,0.5,0.9,YFIXED,1},
    {0.75,0.6,1.0,XFIXED,1}, {0.75,0.6,1.0,YFIXED,1},
    {0.85,0.7,1.0,XFIXED,1}, {0.85,0.7,1.0,YFIXED,1},
    };*/
  /*
 // CROSS3
  int deg = 3, num = 19; 
  const int nrefs = 34;
  double tt = 1.0/(num+1.0);
  LRSplineSurface::Refinement2D rfs[nrefs] = {
    {0.075,0.0,0.2,XFIXED,1}, {0.075,0.0,0.2,YFIXED,1},
    {0.125,0.05,0.25,XFIXED,1}, {0.125,0.05,0.25,YFIXED,1},
    {0.175,0.1,0.3,XFIXED,1}, {0.175,0.1,0.3,YFIXED,1}, 
    {0.225,0.15,0.35,XFIXED,1}, {0.225,0.15,0.35,YFIXED,1},
    {0.275,0.2,0.4,XFIXED,1}, {0.275,0.2,0.4,YFIXED,1}, 
    {0.325,0.25,0.45,XFIXED,1}, {0.325,0.25,0.45,YFIXED,1}, 
    {0.375,0.3,0.5,XFIXED,1}, {0.375,0.3,0.5,YFIXED,1}, 
    {0.425,0.35,0.55,XFIXED,1}, {0.425,0.35,0.55,YFIXED,1}, 
    {0.475,0.4,0.6,XFIXED,1}, {0.475,0.4,0.6,YFIXED,1}, 
    {0.525,0.45,0.65,XFIXED,1}, {0.525,0.45,0.65,YFIXED,1}, 
    {0.575,0.5,0.7,XFIXED,1}, {0.575,0.5,0.7,YFIXED,1}, 
    {0.625,0.55,0.75,XFIXED,1}, {0.625,0.55,0.75,YFIXED,1}, 
    {0.675,0.6,0.8,XFIXED,1}, {0.675,0.6,0.8,YFIXED,1}, 
    {0.725,0.65,0.85,XFIXED,1}, {0.725,0.65,0.85,YFIXED,1}, 
    {0.775,0.7,0.9,XFIXED,1}, {0.775,0.7,0.9,YFIXED,1}, 
    {0.825,0.75,0.95,XFIXED,1}, {0.825,0.75,0.95,YFIXED,1}, 
    {0.875,0.8,1.0,XFIXED,1}, {0.875,0.8,1.0,YFIXED,1}, 
  }; 
*/
  /*
  // CROSS3
  int deg = 3, num = 19; 
  const int nrefs = 34;
  double tt = 1.0/(num+1.0);
  LRSplineSurface::Refinement2D rfs[nrefs] = {
    {0.075,0.0,0.35,XFIXED,1}, {0.075,0.0,0.35,YFIXED,1},
    {0.125,0.0,0.4,XFIXED,1}, {0.125,0.0,0.4,YFIXED,1},
    {0.175,0.075,0.45,XFIXED,1}, {0.175,0.05,0.45,YFIXED,1}, 
    {0.225,0.125,0.5,XFIXED,1}, {0.225,0.125,0.5,YFIXED,1},
    {0.275,0.175,0.55,XFIXED,1}, {0.275,0.175,0.55,YFIXED,1}, 
    {0.325,0.225,0.6,XFIXED,1}, {0.325,0.225,0.6,YFIXED,1}, 
    {0.375,0.275,0.65,XFIXED,1}, {0.375,0.275,0.65,YFIXED,1}, 
    {0.425,0.325,0.7,XFIXED,1}, {0.425,0.325,0.7,YFIXED,1}, 
    {0.475,0.375,0.75,XFIXED,1}, {0.475,0.375,0.75,YFIXED,1}, 
    {0.525,0.425,0.8,XFIXED,1}, {0.525,0.425,0.8,YFIXED,1}, 
    {0.575,0.475,0.85,XFIXED,1}, {0.575,0.475,0.85,YFIXED,1}, 
    {0.625,0.525,0.9,XFIXED,1}, {0.625,0.525,0.9,YFIXED,1}, 
    {0.675,0.575,0.95,XFIXED,1}, {0.675,0.575,0.95,YFIXED,1}, 
    {0.725,0.625,1.0,XFIXED,1}, {0.725,0.625,1.0,YFIXED,1}, 
    {0.775,0.675,1.0,XFIXED,1}, {0.775,0.675,1.0,YFIXED,1}, 
    {0.825,0.725,1.0,XFIXED,1}, {0.825,0.725,1.0,YFIXED,1}, 
    {0.875,0.775,1.0,XFIXED,1}, {0.875,0.775,1.0,YFIXED,1}, 
  };
  */

  /*
  // L1
  int deg = 1, num = 9; 
  const int nrefs = 4;
  double tt = 1.0/(num+1.0);
  LRSplineSurface::Refinement2D rfs[nrefs] = {
    {0.35,0.2,0.9,XFIXED,1}, {0.35,0.2,0.9,YFIXED,1},
    {0.45,0.3,0.9,XFIXED,1}, {0.45,0.3,0.9,YFIXED,1},
    };*/
  /*
  // L2
  int deg = 2, num = 9; 
  const int nrefs = 6;
  double tt = 1.0/(num+1.0);
  LRSplineSurface::Refinement2D rfs[nrefs] = {
    {0.35,0.1,0.9,XFIXED,1}, {0.35,0.1,0.9,YFIXED,1},
    {0.45,0.3,0.9,XFIXED,1}, {0.45,0.3,0.9,YFIXED,1},
    {0.55,0.3,0.9,XFIXED,1}, {0.55,0.3,0.9,YFIXED,1},
    };*/
  /*
  // DEG3LARGE
  int deg = 3, num = 19; 
  const int nrefs = 24;
  double tt = 1.0/(num+1.0);
  LRSplineSurface::Refinement2D rfs[nrefs] = {
    {0.225,0.05,0.95,XFIXED,1}, {0.225,0.05,0.95,YFIXED,1},
    {0.275,0.15,0.85,XFIXED,1}, {0.275,0.15,0.85,YFIXED,1},
    {0.325,0.2,0.8,XFIXED,1}, {0.325,0.2,0.8,YFIXED,1},
    {0.375,0.2,0.8,XFIXED,1}, {0.375,0.2,0.8,YFIXED,1},
    {0.425,0.2,0.8,XFIXED,1}, {0.425,0.2,0.8,YFIXED,1},
    {0.475,0.2,0.8,XFIXED,1}, {0.475,0.2,0.8,YFIXED,1},
    {0.525,0.2,0.8,XFIXED,1}, {0.525,0.2,0.8,YFIXED,1},
    {0.575,0.2,0.8,XFIXED,1}, {0.575,0.2,0.8,YFIXED,1},
    {0.625,0.2,0.8,XFIXED,1}, {0.625,0.2,0.8,YFIXED,1},
    {0.675,0.2,0.8,XFIXED,1}, {0.675,0.2,0.8,YFIXED,1},
    {0.725,0.15,0.85,XFIXED,1}, {0.725,0.15,0.85,YFIXED,1},
    {0.775,0.05,0.95,XFIXED,1}, {0.775,0.05,0.95,YFIXED,1},
 };*/

  //vector<LRSplineSurface::Refinement2D> refs(rfs,rfs+nrefs);


  if (argc < 2 ) {
    cerr << "usage: ./" << argv[0] << " <input_file.txt> " << endl; 
    return -1;
  } 

  string infile = string(argv[1]);
  ifstream ifs(infile);

  int deg, num, nrefs;
  double tt;
  ifs >> deg >> num >> nrefs;

  num -= 1;
  vector<LRSplineSurface::Refinement2D> refs;

  
  for (int ix=0;ix!=nrefs;++ix){
    double fix, start, end;
    int dir_i;
    ifs >> fix >> start >> end;
    ifs >> dir_i;
    //if (dir_i==1) continue;
    LRSplineSurface::Refinement2D ref = {fix,start,end,Direction2D(dir_i),1};
    refs.push_back(ref);
  } 
  //cout << "SZ " << refs.size() << endl;
  cout << deg << " " << num << endl;
  LRSplineSurface lrs = construct_basic_lrspline(deg, num, true);

  set_random_coefficients(lrs);
  lrs.refine(refs);  

  ofstream ofs(infile.substr(0,infile.size()-4)+string(".eps"));
  
  writePostscriptMeshOverload(lrs,ofs);




  std::vector<LRBSpline2D*> bb = LinDepUtils::unpeelableBasisFunctions ( lrs );
  cout << "Is peelable: " << LinDepUtils::isPeelable(lrs) << endl;
  cout << "No. unpeelable basis functions: " << bb.size() << endl;
  return 0;
}

