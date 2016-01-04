//----------------------------------------------------------------------------
// File : MMsIO.h
// Project : MacMapStats
// Purpose : Header file : Classification & analysis read/write utils
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
// 22/05/2008 creation.
//----------------------------------------------------------------------------

#ifndef __MMsIO__
#define __MMsIO__

//----------------------------------------------------------------------------

#ifdef __cplusplus
extern "C" {
#endif

//----------------------------------------------------------------------------

#include <MacMapStats/MacMapStats.h>
#include <MacMapStats/MMsAnalysis.h>
#include <MacMapStats/MMsClassification.h>

//----------------------------------------------------------------------------

int	MMsIOWriteAnalysis			(	const char* path,
									const char* name,
									mmx_analysis* ana);
int	MMsIOReadAnalysis			(	const char* path,
									const char* name,
									mmx_analysis* ana);

int	MMsIOWriteClassification	(	const char* path,
									const char* name,
									hca_clss** clsss,
									int nclss,
									int clssk);
int	MMsIOReadClassification		(	const char* path,
									const char* name,
									hca_clss*** clsss,
									int* nclss,
									int clssk);
int	MMsIOWriteParam				(	const char* path,
									const char* name,
									void* prm,
									int sz,
									const char* pname);
int	MMsIOReadParam				(	const char* path,
									const char* name,
									void** prm,
									int* sz,
									const char* pname);
									
//----------------------------------------------------------------------------
								
#ifdef __cplusplus
}
#endif

//----------------------------------------------------------------------------

#endif
