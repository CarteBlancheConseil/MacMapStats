//----------------------------------------------------------------------------
// File : MMsDiscretization.h
// Project : MacMapStats
// Purpose : Header file : Discretization utils
// Author : Benoit Ogier, benoit.ogier@macmap.com
//
// Copyright (C) 1997-2015 Carte Blanche Conseil.
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.
//
// See the COPYING.lesser file for more information.
//
//----------------------------------------------------------------------------
// 
//----------------------------------------------------------------------------
// 30/03/2007 creation.
//----------------------------------------------------------------------------

#ifndef __MMsDiscretization__
#define __MMsDiscretization__

//----------------------------------------------------------------------------

#ifdef __cplusplus
extern "C" {
#endif

//----------------------------------------------------------------------------

#include <MacMapStats/MacMapStats.h>

//----------------------------------------------------------------------------

enum{
	kDiscretizationLimitColumn_		=1,
	kDiscretizationIndexColumn_		=2
};

//----------------------------------------------------------------------------

int d_estimate_class_count				(	MMsMatrix* mx);

MMsMatrix* d_equal_manning				(	MMsMatrix* mx,
											int c,
											int nbclass);
MMsMatrix* d_equal_range				(	MMsMatrix* mx,
											int c,
											int nbclass);
MMsMatrix* d_standardized				(	MMsMatrix* mx,
											int c,
											int nbclass);
MMsMatrix* d_bertin_average				(	MMsMatrix* mx,
											int c,
											int nbclass);
MMsMatrix* d_bertin_median				(	MMsMatrix* mx,
											int c,
											int nbclass);
MMsMatrix* d_arithmetic_progression		(	MMsMatrix* mx,
											int c,
											int nbclass);
MMsMatrix* d_geometric_progression		(	MMsMatrix* mx,
											int c,
											int nbclass);
MMsMatrix* d_quantile					(	MMsMatrix* mx,
											int c,
											int nbclass);
MMsMatrix* d_natural_breaks				(	MMsMatrix* mx,
											int c,
											int nbclass);

//----------------------------------------------------------------------------
								
#ifdef __cplusplus
}
#endif

//----------------------------------------------------------------------------

#endif
