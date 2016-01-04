//----------------------------------------------------------------------------
// File : MacMapStats.h
// Project : MacMapStats
// Purpose : Header file : Statitics lib
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
// 01/06/1996 creation.
// 30/03/2007 OS9 -> OSX.
//----------------------------------------------------------------------------

#ifndef __MacMapStats__
#define __MacMapStats__

//----------------------------------------------------------------------------

#ifdef __cplusplus
extern "C" {
#endif

//----------------------------------------------------------------------------

// Array kinds
enum{
	kMatrixNoKind_			= 0,
	kMatrixRectangularKind_	= 1,
	kMatrixDiagonalKind_	= 2
};

// Matrix
typedef struct MMsMatrix{
	int		flg;
	int		sign;
	int		nl;
	int		nc;
	double	f[];
}MMsMatrix;

//----------------------------------------------------------------------------
// MMsMatrix Utilities
MMsMatrix*	MMsMatrixAlloc		(	int ni,
									int nj);
MMsMatrix*	MMsDiagMatrixAlloc	(	int ni);
MMsMatrix*	MMsMatrixClone		(	MMsMatrix* mx);
MMsMatrix*	MMsUnitMatrixAlloc	(	int n);
void		MMsMatrixFree		(	MMsMatrix* mx);
MMsMatrix*	MMsMatrixAddLine	(	MMsMatrix* mx);
MMsMatrix*	MMsMatrixInsertLine	(	MMsMatrix* mx, 
									int idx);
MMsMatrix*	MMsMatrixRemoveLine	(	MMsMatrix* mx, 
									int idx);
double		MMsGetDouble		(	MMsMatrix* mx, 
									int i,
									int j);
void		MMsSetDouble		(	MMsMatrix* mx, 
									int i, 
									int j, 
									double d);
MMsMatrix*	MMsCloneColumn		(	MMsMatrix* mx,
									int nj);
void		MMsMatrixSort		(	MMsMatrix* mx);
void		MMsDumpTTxt			(	MMsMatrix* mx,
									FILE* fp);

//----------------------------------------------------------------------------
// 
double		MMsMin				(	int c, 
									MMsMatrix* mx);
double		MMsMax				(	int c, 
									MMsMatrix* mx);
double		MMsAverage			(	int c, 
									MMsMatrix* mx);
double		MMsMedian			(	int c, 
									MMsMatrix* mx);
double		MMsSum				(	int c, 
									MMsMatrix* mx);
double		MMsVariance			(	int c, 
									MMsMatrix* mx);
double		MMsStdDeviation		(	int c, 
									MMsMatrix* mx);
double		MMsSquariance		(	int c, 
									MMsMatrix* mx);
double		MMsCovariance		(	int c1,
									int c2, 
									MMsMatrix* mx);
double		MMsCorrelation		(	int c1,
									int c2, 
									MMsMatrix* mx);
void		MMsMultiply			(	MMsMatrix* mx, 
									double x);
									
MMsMatrix*	MMsMxMultiply		(	MMsMatrix* a, 
									MMsMatrix* b);
MMsMatrix*	MMsMxCorrelation	(	MMsMatrix* mx);
MMsMatrix*	MMsMxCovariance		(	MMsMatrix* mx);
MMsMatrix*	MMsTranspose		(	MMsMatrix* mx);
void		MMsEigenVector		(	MMsMatrix* A, 
									MMsMatrix* *R, 
									MMsMatrix* *V);
									
//----------------------------------------------------------------------------

void		Marginales			(	MMsMatrix* mx, 
									MMsMatrix* *pi, 
									MMsMatrix* *pj, 
									double *s);
MMsMatrix*	Centrees			(	MMsMatrix* mx);
MMsMatrix*	CentreesReduites	(	MMsMatrix* mx);
MMsMatrix*	Normees				(	MMsMatrix* mx);

MMsMatrix*	MatriceS			(	MMsMatrix* p);
MMsMatrix*	NormAFC				(	MMsMatrix* p, 
									MMsMatrix* pi, 
									MMsMatrix* pj, 
									float s);
									
MMsMatrix*	GetACPIndCoord		(	MMsMatrix* A, 
									MMsMatrix* R, 
									MMsMatrix* V, 
									char Scaled);
MMsMatrix*	GetACPVarCoord		(	MMsMatrix* R, 
									MMsMatrix* V);
MMsMatrix*	GetAFCIndCoord		(	MMsMatrix* A, 
									MMsMatrix* R);
MMsMatrix*	GetAFCIndContrib	(	MMsMatrix* A,
									MMsMatrix* V);
MMsMatrix*	GetAFCVarCoord		(	MMsMatrix* r,
									MMsMatrix* v,
									MMsMatrix* pj);
void		RmvVTrivial			(	MMsMatrix* *r);
MMsMatrix*	GetDstIndMat		(	MMsMatrix* mx);

//----------------------------------------------------------------------------

								
#ifdef __cplusplus
}
#endif

//----------------------------------------------------------------------------

#endif
