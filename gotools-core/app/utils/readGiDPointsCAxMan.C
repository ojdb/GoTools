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
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>


using namespace std;


int main( int argc, char* argv[] )
{

  if (argc < 3 || argc > 4)
    {
      cout << "Usage:  " << argv[0] << " result mesh" << endl;
      return 1;
    }

  vector<float> results;  
  vector<int> boundary_indices;

  ifstream in_res(argv[1]);

  std::string line;
  // Don't need first 3 lines
  std::getline(in_res, line);
  std::getline(in_res, line);
  std::getline(in_res, line);

  int index;
  float a,b,c;
  while (std::getline(in_res, line))
  {
    if (line == "End Values") { break; } // eof
	std::istringstream ss(line);
	ss >> index >> a >> b >> c;
	boundary_indices.push_back(index);
	results.push_back(a);
	results.push_back(b);
	results.push_back(c);
  }

  in_res.close();

  ifstream in_mesh(argv[2]);
  vector<float> pts;
  // Ignore first two lines
  std::getline(in_mesh, line);
  std::getline(in_mesh, line);
  int jndex = 0;
  while (std::getline(in_mesh, line))
  {
    if (line == "End Coordinates") { break; } // we don't need any more
    std::istringstream ss(line);
	ss >> index;
	if (index != boundary_indices[jndex]) continue;
	ss >> a >> b >> c;
	pts.push_back(1000.0*a);
	pts.push_back(1000.0*b);
	pts.push_back(1000.0*c);
    jndex++;
  }
  in_mesh.close();

  if (pts.size() == results.size()) std::cout << "files seem to correspond ok" << std::endl;
  else std::cout << "something seems to be wrong with the input files" << std::endl;

  std::vector<float> distorted(pts.size());

  for (int ix=0; ix!=pts.size(); ++ix) {
    distorted[ix] = pts[ix]+results[ix];
  }

  std::ofstream ofs("distorted_points.txt");
  for (int ix=0; ix!=pts.size()/3; ++ix) {
    ofs << distorted[3*ix] << " " << distorted[3*ix+1] << " " << distorted[3*ix+2] << "\n";
  }

}
