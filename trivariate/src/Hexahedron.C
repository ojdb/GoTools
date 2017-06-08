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

#include "GoTools/trivariate/Hexahedron.h"
#include "GoTools/trivariate/SplineVolume.h"


using std::vector;
using std::endl;


namespace Go
{


  // Constructor
  //===========================================================================
  Hexahedron::Hexahedron(Point p1, Point p2, Point p3, Point p4,
                         Point p5, Point p6, Point p7, Point p8) :
    p1_(p1), p2_(p2), p3_(p3), p4_(p4),
    p5_(p5), p6_(p6), p7_(p7), p8_(p8)
  //===========================================================================
  {
    if (p1_.dimension() != 3 ||
	p2_.dimension() != 3 ||
	p3_.dimension() != 3 ||
	p4_.dimension() != 3 ||
        p5_.dimension() != 3 ||
	p6_.dimension() != 3 ||
	p7_.dimension() != 3 ||
	p8_.dimension() != 3)
      {
	THROW("Dimension must be 3.");
	return;
      }
  }

  // Destructor
  //===========================================================================
  Hexahedron::~Hexahedron()
  //===========================================================================
  {
  }

  //===========================================================================
  void Hexahedron::read(std::istream& is)
  //===========================================================================
  {
    if (!is.good()) {
	THROW("Invalid geometry file!");
    }

    int dim;
    is >> dim;
    if (dim != 3)
	THROW("Dimension must be 3.");
    p1_.resize(dim);
    p2_.resize(dim);
    p3_.resize(dim);
    p4_.resize(dim);
    p5_.resize(dim);
    p6_.resize(dim);
    p7_.resize(dim);
    p8_.resize(dim);

    is >> p1_
       >> p2_
       >> p3_
       >> p4_
       >> p5_
       >> p6_
       >> p7_
       >> p8_;
  }

  //===========================================================================
  void Hexahedron::write(std::ostream& os) const
  //===========================================================================
  {
    os << dimension() << endl
       << p1_ << endl
       << p2_ << endl
       << p3_ << endl
       << p4_ << endl
       << p5_ << endl
       << p6_ << endl
       << p7_ << endl
       << p8_ << endl;
  }

  //===========================================================================
  int Hexahedron::dimension() const
  //===========================================================================
  {
    return p1_.dimension();
  }

  //===========================================================================
  ClassType Hexahedron::instanceType() const
  //===========================================================================
  {
    return classType();
  }


  //===========================================================================
  BoundingBox Hexahedron::boundingBox() const
  //===========================================================================
  {
    MESSAGE("boundingBox() not yet implemented");
    BoundingBox bb;
    return bb;
  }

  //===========================================================================
  DirectionCone Hexahedron::tangentCone(int pardir) const
  //===========================================================================
  {
    MESSAGE("tangentCone() not yet implemented");
    DirectionCone dc;
    return dc;
    //if (pardir == 0)
    //  return DirectionCone(dir_u_);
    //else if (pardir == 1)
    //  return DirectionCone(dir_v_);
    //else
    //  return DirectionCone(dir_w_);
  }

  //===========================================================================
  const Array<double,6> Hexahedron::parameterSpan() const
  //===========================================================================
  {
    Array<double,6> pSpan;
    pSpan[0] = pSpan[2] = pSpan[4] = 0.0;
    pSpan[1] = 1.0;
    pSpan[3] = 1.0;
    pSpan[5] = 1.0;
    return pSpan;
  }

  //===========================================================================
  void Hexahedron::point(Point& pt, double upar, double vpar, double wpar) const
  //===========================================================================
  {
    pt = p1_
      + (p2_-p1_)*upar + (p4_-p1_)*vpar + (p5_-p1_)*wpar
      + (p1_-p2_+p3_-p4_)*upar*vpar + (p1_-p2_-p5_+p6_)*upar*wpar + (p1_ - p4_ - p5_ + p8_)*vpar*wpar
      + (-p1_ + p2_ - p3_ + p4_ + p5_ - p6_ + p7_ - p8_)*upar*vpar*wpar;

      //pt =  p1_*(1-upar)*(1-vpar)*(1-wpar)
      //+ p2_*upar*(1-vpar)*(1-wpar)
      //+ p3_*upar*vpar*(1-wpar)
      //+ p4_*(1-upar)*vpar*(1-wpar)
      //+ p5_*(1-upar)*(1-vpar)*wpar
      //+ p6_*upar*(1-vpar)*wpar
      //+ p7_*upar*vpar*wpar
      //+ p8_*(1-upar)*vpar*wpar;
  }

  //===========================================================================
  void Hexahedron::point(vector<Point>& pts, 
                         double upar, double vpar, double wpar,
                         int derivs,
                         bool u_from_right,
                         bool v_from_right,
                         bool w_from_right,
                         double resolution ) const
  //===========================================================================
  {
    DEBUG_ERROR_IF(derivs < 0,
		   "Negative number of derivatives makes no sense.");
    int totpts = (derivs + 1)*(derivs + 2)*(derivs + 3)/6;
    int ptsz = (int)pts.size();
    DEBUG_ERROR_IF(ptsz< totpts,
		   "The vector of points must have sufficient size.");

    int dim = dimension();
    for (int i = 0; i < totpts; ++i) {
        if (pts[i].dimension() != dim) {
            pts[i].resize(dim);
        }
	pts[i].setValue(0.0);
    }

    point(pts[0], upar, vpar, wpar);
    if (derivs == 0)
        return;

    // Derivatives
    pts[1] = (p2_-p1_) + (p1_-p2_+p3_-p4_)*vpar + (p1_-p2_-p5_+p6_)*wpar + (-p1_+p2_-p3_+p4_+p5_-p6_+p7_-p8_)*vpar*wpar;
    pts[2] = (p4_-p1_) + (p1_-p2_+p3_-p4_)*upar + (p1_-p4_-p5_+p8_)*wpar + (-p1_+p2_-p3_+p4_+p5_-p6_+p7_-p8_)*upar*wpar;
    pts[3] = (p5_-p1_) + (p1_-p4_-p5_+p8_)*vpar + (p1_-p2_-p5_+p6_)*upar + (-p1_+p2_-p3_+p4_+p5_-p6_+p7_-p8_)*upar*vpar;
    
    // Second order and higher derivatives vanish. They are already
    // set to zero, so we return.
  }

  //===========================================================================
  double Hexahedron::nextSegmentVal(int dir, double par, bool forward, double tol) const
  //===========================================================================
  {
    MESSAGE("nextSegmentVal() not yet implemented");
    return 0.0;
  }

  //===========================================================================
  void Hexahedron::closestPoint(const Point& pt,
                                double& clo_u,
                                double& clo_v,
                                double& clo_w,
                                Point& clo_pt,
                                double& clo_dist,
                                double epsilon,
                                double *seed) const
  //===========================================================================
  {
    MESSAGE("closestPoint() not yet implemented");
  }

  //===========================================================================
  void Hexahedron::reverseParameterDirection(int pardir)
  //===========================================================================
  {
    MESSAGE("reverseParameterDirection() not yet implemented");
    //if (pardir == 0)
    //  {
    //  corner_ += length_u_ * dir_u_;
    //  dir_u_ = -dir_u_;
    //  }
    //else if (pardir == 1)
    //  {
    //  corner_ += length_v_ * dir_v_;
    //  dir_v_ = -dir_v_;
    //  }
    //else if (pardir == 2)
    //  {
    //  corner_ += length_w_ * dir_w_;
    //  dir_w_ = -dir_w_;
    //  }
  }

  //===========================================================================
  void Hexahedron::swapParameterDirection(int pardir1, int pardir2)
  //===========================================================================
  {
    MESSAGE("swapParameterDirection() not yet implemented");
    //if ((pardir1 == 0 && pardir2 == 1) || (pardir1 == 1 && pardir2 == 0))
    //  {
//	double tmp_l = length_u_;
//	length_u_ = length_v_;
//	length_v_ = tmp_l;
//	Point tmp_d = dir_u_;
//	dir_u_ = dir_v_;
//	dir_v_ = tmp_d;
  //    }
    //else if ((pardir1 == 0 && pardir2 == 2) || (pardir1 == 2 && pardir2 == 0))
    //  {
//	double tmp_l = length_u_;
//	length_u_ = length_w_;
//	length_w_ = tmp_l;
//	Point tmp_d = dir_u_;
//	dir_u_ = dir_w_;
//	dir_w_ = tmp_d;
    //  }
    //else if ((pardir1 == 1 && pardir2 == 2) || (pardir1 == 2 && pardir2 == 1))
    //  {
//	double tmp_l = length_v_;
//	length_v_ = length_w_;
//	length_w_ = tmp_l;
//	Point tmp_d = dir_v_;
//	dir_v_ = dir_w_;
//	dir_w_ = tmp_d;
    //  }
  }


  //===========================================================================
  vector<shared_ptr<ParamSurface> > Hexahedron::getAllBoundarySurfaces() const
  //===========================================================================
  {
    MESSAGE("getAllBoundarySurfaces() not implemented.");
    vector<shared_ptr<ParamSurface> > bound_surf;
    return bound_surf;
  }

  //===========================================================================
  void Hexahedron::translate(const Point& vec)
  //===========================================================================
  {
    ALWAYS_ERROR_IF(dimension() != vec.dimension(), "Volume and translation vector of different dimension");
    p1_ += vec;
    p2_ += vec;
    p3_ += vec;
    p4_ += vec;
    p5_ += vec;
    p6_ += vec;
    p7_ += vec;
    p8_ += vec;
  }

  //===========================================================================
  SplineVolume* Hexahedron::geometryVolume() const
  //===========================================================================
  {
    vector<Point> allCorners(8);
    allCorners[0] = p1_;
    allCorners[1] = p2_;
    allCorners[2] = p4_;
    allCorners[3] = p3_;
    allCorners[4] = p5_;
    allCorners[5] = p6_;
    allCorners[6] = p8_;
    allCorners[7] = p7_;

    vector<double> coefs(24);
    for (int i = 0, pos = 0; i < 8; ++i)
      for (int j = 0; j< 3; ++j, ++pos)
	coefs[pos] = allCorners[i][j];

    vector<double> knots_u(4), knots_v(4), knots_w(4);
    knots_u[0] = knots_u[1] = knots_v[0] = knots_v[1] = knots_w[0] = knots_w[1] = 0.0;
    knots_u[2] = knots_u[3] = 1.0;
    knots_v[2] = knots_v[3] = 1.0;
    knots_w[2] = knots_w[3] = 1.0;

    return new SplineVolume(2, 2, 2, 2, 2, 2, knots_u.begin(), knots_v.begin(), knots_w.begin(), coefs.begin(), 3);
  }



} // namespace Go

