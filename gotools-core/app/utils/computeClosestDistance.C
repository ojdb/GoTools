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



#include <vector>
#include <fstream>
#include "GoTools/geometry/ObjectHeader.h"
#include "GoTools/geometry/GoTools.h"
//#include "GoTools/geometry/GeomObject.h"
#include "GoTools/geometry/Factory.h"
#include "GoTools/geometry/Utils.h"
#include "GoTools/utils/RegistrationUtils.h"
#include "GoTools/utils/ClosestPointUtils.h"



using namespace Go;
using namespace std;


typedef pair<vector<vector<double> >, Point>  transformation_type;

transformation_type currentTransformation;

static bool abs_comp(float a, float b)
{
    return (std::fabs(a) < std::fabs(b));
}


int main( int argc, char* argv[] )
{
  GoTools::init();

  if (argc < 3 || argc > 4)
    {
      cout << "Usage:  " << argv[0] << " surfaceFile pointFile" << endl;
      return 1;
    }

  ifstream in_surf(argv[1]);
  ObjectHeader header;
  vector<shared_ptr<GeomObject> > surfaces;

  while (!in_surf.eof())
    {
      header.read(in_surf);
      shared_ptr<GeomObject> obj(Factory::createObject(header.classType()));
      obj->read(in_surf);
      surfaces.push_back(obj);
      Utils::eatwhite(in_surf);
    }
  in_surf.close();

  ifstream in_pts(argv[2]);
  vector<float> pts;

  while (!in_pts.eof())
    {
      for (int j = 0; j < 3; ++j)
	{
	  float f;
	  in_pts >> f;
	  pts.push_back(f);
	}
      Utils::eatwhite(in_pts);
    }
  in_pts.close();

  vector<vector<double> > startRotation;
  Point startTranslation;

  startTranslation = Point(0.0, 0.0, 0.0);
  startRotation.resize(3);
  for(int i = 0; i < 3; i++) {
	startRotation[i].resize(3);
	  for (int j = 0; j < 3; ++j)
	    startRotation[i][j] = i==j;
  }

  currentTransformation = transformation_type(startRotation, startTranslation);

  shared_ptr<boxStructuring::BoundingBoxStructure> structure = preProcessClosestVectors(surfaces, 200.0);
  
  vector<float> closest_pts = closestSignedDistances(pts, structure, currentTransformation.first, currentTransformation.second);

  auto abs_max_it = max_element(std::begin(closest_pts), std::end(closest_pts), abs_comp);
  auto abs_min_it = min_element(std::begin(closest_pts), std::end(closest_pts), abs_comp);

  std::cout << "abs min, abs max = " << *abs_min_it << ", " << *abs_max_it << std::endl;  

  auto max_it = max_element(std::begin(closest_pts), std::end(closest_pts));
  auto min_it = min_element(std::begin(closest_pts), std::end(closest_pts));

  std::cout << "min, max = " << *min_it << ", " << *max_it << std::endl;  
  
  ofstream out_str("distances.txt");
  for (int i = 0; i < closest_pts.size(); i += 1)
      out_str << closest_pts[i] << "\n";
	
  std::cout << "results have been written to distance.txt" << std::endl; 

}
