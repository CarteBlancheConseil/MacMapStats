//----------------------------------------------------------------------------
// File : MMsDiscretization.c
// Project : MacMapStats
// Purpose : C source file : Discretization utils
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

#include "MMsDiscretization.h"
#include "MMsMatrix_Utils.h"

//----------------------------------------------------------------------------
// 
//------------
int d_estimate_class_count(MMsMatrix* mx){
double	e=1.0+3.3*log10(mx->nl);
int		k=(e-0.5);
	return((k<e)?k:k-1);
}

//----------------------------------------------------------------------------
// 
//------------
MMsMatrix* d_equal_manning(MMsMatrix* mx,
							int c,
							int nbclass){
fprintf(stderr,"d_equal_manning\n");
MMsMatrix* buf=get_matrix(mx,c);
	if(buf==NULL){
		return(NULL);
	}
MMsMatrix* mxr=MMsMatrixAlloc(nbclass+1,2);
	if(mxr==NULL){
		MMsMatrixFree(buf);
		return(NULL);
	}

int		i;
int		n=buf->nl/nbclass,x=1;
	
	MMsSetDouble(mxr,1,kDiscretizationLimitColumn_,MMsGetDouble(buf,1,1));
	MMsSetDouble(mxr,1,kDiscretizationIndexColumn_,x);
	x+=n;
	for(i=2;i<=nbclass+1;i++){
		MMsSetDouble(mxr,i,kDiscretizationLimitColumn_,MMsGetDouble(buf,x,1));
		MMsSetDouble(mxr,i,kDiscretizationIndexColumn_,x);
		x+=n;
		if(x>buf->nl){
			x=buf->nl;
		}
	}

	MMsMatrixFree(buf);
	return(mxr);
}

//----------------------------------------------------------------------------
// 
//------------
MMsMatrix* d_equal_range(	MMsMatrix* mx,
							int c,
							int nbclass){
fprintf(stderr,"d_equal_range in\n");
MMsMatrix* buf=get_matrix(mx,c);
	if(buf==NULL){
		return(NULL);
	}
MMsMatrix* mxr=MMsMatrixAlloc(nbclass+1,2);
	if(mxr==NULL){
		MMsMatrixFree(buf);
		return(NULL);
	}

double	min=MMsGetDouble(buf,1,1);
double	amp=(MMsGetDouble(buf,buf->nl,1)-min)/(double)nbclass;
int		i;

	MMsSetDouble(mxr,1,kDiscretizationLimitColumn_,min);
	MMsSetDouble(mxr,nbclass+1,kDiscretizationLimitColumn_,MMsGetDouble(buf,buf->nl,1));

	for(i=2;i<=nbclass;i++){
		min+=amp;
		MMsSetDouble(mxr,i,kDiscretizationLimitColumn_,min);
	}

int		cl=2;
	min=MMsGetDouble(mxr,cl,kDiscretizationLimitColumn_);
	MMsSetDouble(mxr,1,kDiscretizationIndexColumn_,1);
	for(i=1;i<=buf->nl;i++){
		if(MMsGetDouble(buf,i,1)>=min){
			MMsSetDouble(mxr,cl,kDiscretizationIndexColumn_,i);
			while((MMsGetDouble(buf,i,1)>=min)&&(cl<mxr->nl)){
				cl++;
				min=MMsGetDouble(mxr,cl,kDiscretizationLimitColumn_);
				MMsSetDouble(mxr,cl,kDiscretizationIndexColumn_,i);
			}
		}
	}
	MMsMatrixFree(buf);
fprintf(stderr,"d_equal_range out\n");

	return(mxr);
}

//----------------------------------------------------------------------------
// 
//------------
MMsMatrix* d_standardized(	MMsMatrix* mx,
							int c,
							int nbclass){
fprintf(stderr,"d_standardized\n");
MMsMatrix* buf=get_matrix(mx,c);
	if(buf==NULL){
		return(NULL);
	}
MMsMatrix* bbuf=MMsMatrixAlloc(mx->nl,2);
	if(bbuf==NULL){
		MMsMatrixFree(buf);
		return(NULL);
	}

double	et=MMsStdDeviation(1,buf);
double	b=MMsGetDouble(buf,1,1);
double	e=MMsAverage(1,buf)-et/2.0;

	do{
		e-=et;
	}
	while(e>=b)/*was : (e>b)*/;
	if(e<b){
		e+=et;
	}

int		i=0,k=1;
	MMsSetDouble(bbuf,k,kDiscretizationLimitColumn_,b);
	MMsSetDouble(bbuf,k,kDiscretizationIndexColumn_,1);
	b=MMsGetDouble(buf,buf->nl,1);
	do{
		do{
			i++;
		}
		while((MMsGetDouble(buf,i,1)<e)&&(i<=buf->nl));
		i--;
		k++;
		MMsSetDouble(bbuf,k,kDiscretizationLimitColumn_,e);
		MMsSetDouble(bbuf,k,kDiscretizationIndexColumn_,i);		
		e+=et;
	}
	while(e<=b);
	if(e>b){
		k++;
		MMsSetDouble(bbuf,k,kDiscretizationLimitColumn_,e);
		MMsSetDouble(bbuf,k,kDiscretizationIndexColumn_,bbuf->nl);			
	}
	
MMsMatrix* mxr=MMsMatrixAlloc(k,2);
	if(mxr==NULL){
		MMsMatrixFree(bbuf);
		MMsMatrixFree(buf);
		return(NULL);
	}

	for(i=1;i<=k;i++){
		MMsSetDouble(mxr,i,kDiscretizationLimitColumn_,MMsGetDouble(bbuf,i,kDiscretizationLimitColumn_));
		MMsSetDouble(mxr,i,kDiscretizationIndexColumn_,MMsGetDouble(bbuf,i,kDiscretizationIndexColumn_));	
	}
	
	MMsMatrixFree(buf);
	MMsMatrixFree(bbuf);
	return(mxr);
}

//----------------------------------------------------------------------------
// 
//------------
MMsMatrix* d_bertin_average(	MMsMatrix* mx,
								int c,
								int nbclass){
fprintf(stderr,"d_bertin_average\n");
MMsMatrix* buf=get_matrix(mx,c);
	if(buf==NULL){
		return(NULL);
	}
MMsMatrix* mxr=MMsMatrixAlloc(nbclass+1,2);
	if(mxr==NULL){
		MMsMatrixFree(buf);
		return(NULL);
	}

double	cnt=MMsAverage(1,buf);
double	ampl=(cnt-MMsGetDouble(buf,1,1))/(nbclass/2.0);
double	amph=(MMsGetDouble(buf,buf->nl,1)-cnt)/(nbclass/2.0);	
double	e=cnt-(ampl*nbclass/2.0);
int		i;

//	e=cnt-(ampl*nbclass/2.0);
	for(i=1;i<=floor((nbclass+1)/2);i++){
		MMsSetDouble(mxr,i,kDiscretizationLimitColumn_,e);
		e+=ampl;
	}	
	e=cnt+(amph*nbclass/2.0);
	for(i=nbclass+1;i>ceil((nbclass+1)/2);i--){
		MMsSetDouble(mxr,i,kDiscretizationLimitColumn_,e);
		e-=amph;
	}
	
int		cl=2;
	e=MMsGetDouble(mxr,cl,kDiscretizationLimitColumn_);
	MMsSetDouble(mxr,1,kDiscretizationIndexColumn_,1);
	for(i=1;i<=buf->nl;i++){
		if(MMsGetDouble(buf,i,1)>=e){
			MMsSetDouble(mxr,cl,kDiscretizationIndexColumn_,i);
			while((MMsGetDouble(buf,i,1)>=e)&&(cl<mxr->nl)){
				cl++;
				e=MMsGetDouble(mxr,cl,kDiscretizationLimitColumn_);
				MMsSetDouble(mxr,cl,kDiscretizationIndexColumn_,i);
			}
		}
	}
	
	MMsMatrixFree(buf);
	return(mxr);
}

//----------------------------------------------------------------------------
// 
//------------
MMsMatrix* d_bertin_median(	MMsMatrix* mx,
							int c,
							int nbclass){
fprintf(stderr,"d_bertin_median\n");
MMsMatrix* buf=get_matrix(mx,c);
	if(buf==NULL){
		return(NULL);
	}
MMsMatrix* mxr=MMsMatrixAlloc(nbclass+1,2);
	if(mxr==NULL){
		MMsMatrixFree(buf);
		return(NULL);
	}
	
//double	cnt=(nbclass&1)?
double	cnt=(buf->nl&1)?
			MMsGetDouble(buf,(buf->nl/2)+1,1):
			(MMsGetDouble(buf,(buf->nl/2),1)-MMsGetDouble(buf,(buf->nl/2)+1,1))/2.0;
double	ampl=(cnt-MMsGetDouble(buf,1,1))/(nbclass/2.0);
double	amph=(MMsGetDouble(buf,buf->nl,1)-cnt)/(nbclass/2.0);	
double	e=cnt-(ampl*nbclass/2.0);	
int		i;

//	e=cnt-(ampl*nbclass/2.0);
	for(i=1;i<=floor((nbclass+1)/2);i++){
		MMsSetDouble(mxr,i,kDiscretizationLimitColumn_,e);
		e+=ampl;
	}	
	e=cnt+(amph*nbclass/2.0);
	for(i=nbclass+1;i>ceil((nbclass+1)/2);i--){
		MMsSetDouble(mxr,i,kDiscretizationLimitColumn_,e);
		e-=amph;
	}
	
int		cl=2;
	e=MMsGetDouble(mxr,cl,kDiscretizationLimitColumn_);
	MMsSetDouble(mxr,1,kDiscretizationIndexColumn_,1);
	for(i=1;i<=buf->nl;i++){
		if(MMsGetDouble(buf,i,1)>=e){
			MMsSetDouble(mxr,cl,kDiscretizationIndexColumn_,i);
			while((MMsGetDouble(buf,i,1)>=e)&&(cl<mxr->nl)){
				cl++;
				e=MMsGetDouble(mxr,cl,kDiscretizationLimitColumn_);
				MMsSetDouble(mxr,cl,kDiscretizationIndexColumn_,i);
			}
		}
	}
	
	MMsMatrixFree(buf);
	return(mxr);
}

//----------------------------------------------------------------------------
// 
//------------
MMsMatrix* d_arithmetic_progression(	MMsMatrix* mx,
										int c,
										int nbclass){
fprintf(stderr,"d_arithmetic_progression\n");
MMsMatrix* buf=get_matrix(mx,c);
	if(buf==NULL){
		return(NULL);
	}
MMsMatrix* mxr=MMsMatrixAlloc(nbclass+1,2);
	if(mxr==NULL){
		MMsMatrixFree(buf);
		return(NULL);
	}

double	min=MMsGetDouble(buf,1,1);
double	amp=(MMsGetDouble(buf,buf->nl,1)-min)/((double)nbclass*(double)(nbclass+1)/2.0);
int		i;

	MMsSetDouble(mxr,1,kDiscretizationLimitColumn_,min);
	MMsSetDouble(mxr,nbclass+1,kDiscretizationLimitColumn_,MMsGetDouble(buf,buf->nl,1));
	
	for(i=2;i<=nbclass;i++){
		MMsSetDouble(mxr,i,kDiscretizationLimitColumn_,min+amp*(double)(i-1));
		min=min+amp*(double)(i-1);
	}

int		cl=2;
	min=MMsGetDouble(mxr,cl,kDiscretizationLimitColumn_);
	MMsSetDouble(mxr,1,kDiscretizationIndexColumn_,1);
	for(i=1;i<=buf->nl;i++){
		if(MMsGetDouble(buf,i,1)>=min){
			MMsSetDouble(mxr,cl,kDiscretizationIndexColumn_,i);
			while((MMsGetDouble(buf,i,1)>=min)&&(cl<mxr->nl)){
				cl++;
				min=MMsGetDouble(mxr,cl,kDiscretizationLimitColumn_);
				MMsSetDouble(mxr,cl,kDiscretizationIndexColumn_,i);
			}
		}
	}

	MMsMatrixFree(buf);
	return(mxr);
}

//----------------------------------------------------------------------------
// 
//------------
MMsMatrix* d_geometric_progression(	MMsMatrix* mx,
									int c,
									int nbclass){
fprintf(stderr,"d_geometric_progression\n");
MMsMatrix* buf=get_matrix(mx,c);
	if(buf==NULL){
		return(NULL);
	}
MMsMatrix* mxr=MMsMatrixAlloc(nbclass+1,2);
	if(mxr==NULL){
		MMsMatrixFree(buf);
		return(NULL);
	}

double	min=MMsGetDouble(buf,1,1);
double	amp=pow(10,(log10(MMsGetDouble(buf,buf->nl,1))-log10(min))/(double)(nbclass+1));
int		i;

	MMsSetDouble(mxr,1,kDiscretizationLimitColumn_,min);
	MMsSetDouble(mxr,nbclass+1,kDiscretizationLimitColumn_,MMsGetDouble(buf,buf->nl,1));
	
	for(i=2;i<=nbclass;i++){
		MMsSetDouble(mxr,i,kDiscretizationLimitColumn_,min*amp);
		min=min*amp;
	}

int		cl=2;
	min=MMsGetDouble(mxr,cl,kDiscretizationLimitColumn_);
	MMsSetDouble(mxr,1,kDiscretizationIndexColumn_,1);
	for(i=1;i<=buf->nl;i++){
		if(MMsGetDouble(buf,i,1)>=min){
			MMsSetDouble(mxr,cl,kDiscretizationIndexColumn_,i);
			while((MMsGetDouble(buf,i,1)>=min)&&(cl<mxr->nl)){
				cl++;
				min=MMsGetDouble(mxr,cl,kDiscretizationLimitColumn_);
				MMsSetDouble(mxr,cl,kDiscretizationIndexColumn_,i);
			}
		}
	}

	MMsMatrixFree(buf);
	return(mxr);
}

//----------------------------------------------------------------------------
// 
//------------
MMsMatrix* d_quantile(	MMsMatrix* mx,
						int c,
						int nbclass){
fprintf(stderr,"d_quantile\n");
MMsMatrix* buf=get_matrix(mx,c);
	if(buf==NULL){
		return(NULL);
	}
MMsMatrix* mxr=MMsMatrixAlloc(nbclass+1,2);
	if(mxr==NULL){
		MMsMatrixFree(buf);
		return(NULL);
	}

double	xcil=((double)nbclass)/((double)buf->nl);
double	min/*=MMsGetDouble(buf,1,1)*/;
	MMsSetDouble(mxr,1,kDiscretizationLimitColumn_,MMsGetDouble(buf,1,1));
	MMsSetDouble(mxr,nbclass+1,kDiscretizationLimitColumn_,MMsGetDouble(buf,buf->nl,1));

int		i,k;
	for(i=1;i<=nbclass;i++){
		k=floor((double)(i-1)/xcil)+1;
		if(k<buf->nl){
			MMsSetDouble(mxr,i,kDiscretizationLimitColumn_,(MMsGetDouble(buf,k+1,1)+MMsGetDouble(buf,k,1))/2.0);
		}
		else{
			MMsSetDouble(mxr,i,kDiscretizationLimitColumn_,MMsGetDouble(buf,k,1));
		}
	}

int		cl=2;
	min=MMsGetDouble(mxr,cl,kDiscretizationLimitColumn_);
	MMsSetDouble(mxr,1,kDiscretizationIndexColumn_,1);
	for(i=1;i<=buf->nl;i++){
		if(MMsGetDouble(buf,i,1)>=min){
			MMsSetDouble(mxr,cl,kDiscretizationIndexColumn_,i);
			while((MMsGetDouble(buf,i,1)>=min)&&(cl<mxr->nl)){
				cl++;
				min=MMsGetDouble(mxr,cl,kDiscretizationLimitColumn_);
				MMsSetDouble(mxr,cl,kDiscretizationIndexColumn_,i);
			}
		}
	}

	MMsMatrixFree(buf);
	return(mxr);
}


//----------------------------------------------------------------------------
// 
//------------
typedef struct break_r{
	double	step;
	int		index;
}break_r;


//----------------------------------------------------------------------------
// 
//------------
static int step_comp(	const void* a,
						const void* b){
break_r*	oa=(break_r*)a;
break_r*	ob=(break_r*)b;
	return(ob->step-oa->step);
}

//----------------------------------------------------------------------------
// 
//------------
static int index_comp(	const void* a,
						const void* b){
break_r*	oa=(break_r*)a;
break_r*	ob=(break_r*)b;
	return(oa->index-ob->index);
}

//----------------------------------------------------------------------------
// 
//------------
MMsMatrix* d_natural_breaks(	MMsMatrix* mx,
								int c,
								int nbclass){
fprintf(stderr,"d_natural_breaks\n");
MMsMatrix* buf=get_matrix(mx,c);
	if(buf==NULL){
		return(NULL);
	}
	
MMsMatrix* mxr=MMsMatrixAlloc(nbclass+1,2);
	if(mxr==NULL){
		MMsMatrixFree(buf);
		return(NULL);
	}

break_r* rank=(break_r*)malloc((buf->nl-1)*sizeof(break_r));
	if(rank==NULL){
		MMsMatrixFree(mxr);
		MMsMatrixFree(buf);
		return(NULL);
	}
int	i;
	for(i=1;i<buf->nl;i++){
		rank[i-1].step=MMsGetDouble(buf,i+1,1)-MMsGetDouble(buf,i,1);
		rank[i-1].index=i+1;
	}
	qsort(rank,buf->nl-1,sizeof(break_r),step_comp);
	rank=(break_r*)realloc(rank,(nbclass-1)*sizeof(break_r));
	qsort(rank,(nbclass-1),sizeof(break_r),index_comp);
	
	MMsSetDouble(mxr,1,kDiscretizationLimitColumn_,MMsGetDouble(buf,1,1));
	MMsSetDouble(mxr,1,kDiscretizationIndexColumn_,1);
	for(i=2;i<=nbclass;i++){
		MMsSetDouble(mxr,i,kDiscretizationLimitColumn_,MMsGetDouble(buf,rank[i-2].index,1));
		MMsSetDouble(mxr,i,kDiscretizationIndexColumn_,rank[i-2].index);
	}
	MMsSetDouble(mxr,nbclass+1,kDiscretizationLimitColumn_,MMsGetDouble(buf,buf->nl,1));
	MMsSetDouble(mxr,nbclass+1,kDiscretizationIndexColumn_,buf->nl);

	free(rank);
	MMsMatrixFree(buf);
	return(mxr);
}
