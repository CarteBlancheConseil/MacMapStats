//----------------------------------------------------------------------------
// File : MMsAnalysis.h
// Project : MacMapStats
// Purpose : Header file : Statistical analysis
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
// 29/08/2007 creation.
//----------------------------------------------------------------------------

#ifndef __MMsAnalysis__
#define __MMsAnalysis__

//----------------------------------------------------------------------------

#ifdef __cplusplus
extern "C" {
#endif

//----------------------------------------------------------------------------

#include <MacMapStats/MacMapStats.h>

//----------------------------------------------------------------------------

enum{
	_kAnalysisUndef	=0,
	_kAnalysisUNI	=1,
	_kAnalysisBIV	=2,
	_kAnalysisPC	=3,
	_kAnalysisPCn	=4,
	_kAnalysisFc	=5,
	_kAnalysisHC	=6,
	_kAnalysisIRA	=7,
	_kAnalysisMLT	=8
};

enum{
	_kMinIndex			=1,
	_kMaxIndex			=2,
	_kExtendIndex		=3,
	_kSumIndex			=4,
	_kAverageIndex		=5,
	_kMedianIndex		=6,
	_kStdDeviationIndex	=7,
	_kVarianceIndex		=8,
	_kcVariationIndex	=9,
	_kcSkewnessIndex	=10,
	_kcPearsonIndex		=11,
	_kcFisherIndex		=12,
	_kcCorrelationIndex	=13,
	_kCovarianceIndex	=14
};

typedef struct mmx_analysis{
	int			kind;
	MMsMatrix*	data;
	MMsMatrix*	indic;
}mmx_analysis;

typedef struct mmx_uni{
	int			kind;
	MMsMatrix*	data;
	MMsMatrix*	indic;
	MMsMatrix*	clss;
}mmx_uni;

typedef struct mmx_biv{
	int			kind;
	MMsMatrix*	data;
	MMsMatrix*	indic;
}mmx_biv;

typedef struct mmx_pca{
	int			kind;
	MMsMatrix*	data;
	MMsMatrix*	indic;
	MMsMatrix*	cv;		// corrélation des variables
	MMsMatrix*	chr;	// individus
	MMsMatrix*	var;	// variables
}mmx_pca;

typedef struct mmx_fa{
	int			kind;
	MMsMatrix*	data;
	MMsMatrix*	indic;
	MMsMatrix*	ctr;	// contributions
	MMsMatrix*	chr;	// individus
	MMsMatrix*	var;	// variables
}mmx_fa;

typedef struct mmx_hc{
	int			kind;
	MMsMatrix*	data;
	MMsMatrix*	indic;
	MMsMatrix*	dst;
}mmx_hc;

typedef struct mmx_ira{
	int			kind;
	MMsMatrix*	data;
	MMsMatrix*	indic;
	MMsMatrix*	center;
}mmx_ira;

typedef struct mmx_mlt{
	int			kind;
	MMsMatrix*	data;
	MMsMatrix*	indic;
}mmx_mlt;
	
//----------------------------------------------------------------------------

void uni_analysis	(	MMsMatrix *indata, // Univariate
						mmx_uni* outdata);
void biv_analysis	(	MMsMatrix *indata, // Bivariate
						mmx_biv* outdata);

void pcn_analysis	(	MMsMatrix *indata, // Principal Component Analysis (Normed)
						mmx_pca* outdata);
void pc_analysis	(	MMsMatrix *indata, // Principal Component Analysis
						mmx_pca* outdata);
void f_analysis		(	MMsMatrix *indata, // Factor Analysis
						mmx_fa* outdata);

void hc_analysis	(	MMsMatrix *indata, // Hierarchical Cluster Analysis
						mmx_hc* outdata);
void ira_analysis	(	MMsMatrix *indata, // Iterative Relocation Algorithm
						mmx_ira* outdata);

void mlt_analysis	(	MMsMatrix *indata, // Multivariate
						mmx_mlt* outdata);

void init_analysis	(	mmx_analysis* data,
						int kind);
void free_analysis	(	mmx_analysis* data);

//----------------------------------------------------------------------------

void linear_reg		(	mmx_biv* ana,
						double* a,
						double* b);
void logarithmic_reg(	mmx_biv* ana,
						double* a,
						double* b);
					
//----------------------------------------------------------------------------
								
#ifdef __cplusplus
}
#endif

//----------------------------------------------------------------------------

#endif
