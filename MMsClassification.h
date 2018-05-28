//----------------------------------------------------------------------------
// File : MMsClassification.h
// Project : MacMapStats
// Purpose : Header file : Classifications
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
// 05/05/2008 creation.
//----------------------------------------------------------------------------

#ifndef __MMsClassification__
#define __MMsClassification__

//----------------------------------------------------------------------------

#ifdef __cplusplus
extern "C" {
#endif

//----------------------------------------------------------------------------

#include <MacMapStats/MacMapStats.h>
#include <MacMapStats/MMsAnalysis.h>

//----------------------------------------------------------------------------

enum{
	kHCAMethodNone		=0,
	kHCADiameterMethod	=1,
	kHCAMinimalMethod	=2,
	kHCAAverageMethod	=3,
	kIRAStdMethod		=4
};

typedef struct hca_clss hca_clss;

struct hca_clss{
	int			cid;		// class id (in the tree)
	int			cidx;		// class index (in the tree)
	
	hca_clss*	sup;		// parent class
	hca_clss*	left;		// left branch
	hca_clss*	right;		// right branch
	
	double		height;		// class height (-1 for inital classes)
	double		midwidth;	// class midwidth
	
	int			ref;		// line index (only for initial classes, 0 for others)
	
	int			ncmp;		// nb initial classes
	hca_clss**	cmp;		// initial component
};
    
typedef struct hca_clss_io{ // Read/Write for 64 bits
    int			cid;		// class id (in the tree)
    int			cidx;		// class index (in the tree)
    
    int         sup;		// parent class
    int         left;		// left branch
    int         right;		// right branch
    
    double		height;		// class height (-1 for inital classes)
    double		midwidth;	// class midwidth
    
    int			ref;		// line index (only for initial classes, 0 for others)
    
    int			ncmp;		// nb initial classes
    int         cmp;		// initial component
}
hca_clss_io;


//----------------------------------------------------------------------------

// Hierarchical Cluster Analysis
hca_clss* hca_diameter	(	MMsMatrix* dst); 
hca_clss* hca_min		(	MMsMatrix* dst); 
hca_clss* hca_average	(	MMsMatrix* dst); 
void hca_iterate		(	hca_clss* tree, 
							void(*proc)(hca_clss*,void*),
							void* prm);
void hca_classfree		(	hca_clss* tree);

// Iterative Relocation Algorithm
hca_clss** ira_calc		(	mmx_ira* data); 
void ira_classfree		(	hca_clss* tree);

// Utils
hca_clss* hca_classnew  (   int ref,
                            int cid);
hca_clss** hca_stackinit(   int n,
                            int maj);
hca_clss** hca_stackadd	(	hca_clss** clsss,
							hca_clss* clss, 
							int n);
hca_clss** hca_stackrmv	(	hca_clss** clsss, 
							int idx, 
							int n);

//----------------------------------------------------------------------------
								
#ifdef __cplusplus
}
#endif

//----------------------------------------------------------------------------

#endif
