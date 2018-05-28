//----------------------------------------------------------------------------
// File : MMsClassification.c
// Project : MacMapStats
// Purpose : C source file : Classifications
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

#include "MMsClassification.h"
#include "MMsAnalysis.h"

#pragma mark -> HCA
//----------------------------------------------------------------------------
// 
//------------
hca_clss* hca_classnew(int ref, int cid){
hca_clss*	clss=(hca_clss*)malloc(sizeof(hca_clss));
	clss->cid=cid;
	clss->cidx=cid;

	clss->sup=NULL;
	clss->left=NULL;
	clss->right=NULL;
	
	clss->height=-1;
	clss->midwidth=-1;

	clss->ref=ref;
	
	clss->ncmp=0;
	clss->cmp=NULL;
	return(clss);
}

//----------------------------------------------------------------------------
// 
//------------
hca_clss** hca_stackinit(int n, int maj){
int			i;
hca_clss**	clsss=(hca_clss**)malloc(sizeof(hca_clss*)*n);
	for(i=0;i<n;i++){
		clsss[i]=hca_classnew(i+1,i+1+maj);
	}
	return(clsss);
}

//----------------------------------------------------------------------------
// 
//------------
hca_clss** hca_stackadd(hca_clss** clsss, hca_clss* clss, int n){
hca_clss**	bclsss;
	if(clsss==NULL){
		bclsss=(hca_clss**)malloc(sizeof(hca_clss*));
	}
	else{
		bclsss=(hca_clss**)realloc(clsss,(n+1)*sizeof(hca_clss*));
	}
	if(!bclsss){
		return(NULL);
	}
	bclsss[n]=clss;
	return(bclsss);
}

//----------------------------------------------------------------------------
// 
//------------
hca_clss** hca_stackrmv(hca_clss** clsss, int idx, int n){
	if(idx<n){
		memmove((void*)&clsss[idx],(const void*)&clsss[idx+1],(size_t)(n-(idx+1))*sizeof(hca_clss*));
	}
hca_clss**	bclsss=(hca_clss**)realloc(clsss,(n-1)*sizeof(hca_clss*));
	if(!bclsss){
		return(NULL);
	}
	return(bclsss);
}

//----------------------------------------------------------------------------
// 
//------------
static void reports(hca_clss* cl, hca_clss* clc){
int	i;
	if(clc->ref>0){
		cl->cmp=hca_stackadd(cl->cmp,clc,cl->ncmp);
		cl->ncmp++;
	}
	else{
		for(i=0;i<clc->ncmp;i++){
			cl->cmp=hca_stackadd(cl->cmp,clc->cmp[i],cl->ncmp);
			cl->ncmp++;
		}
	}
}

//----------------------------------------------------------------------------
// 
//------------
static double max_dst(hca_clss* cl1, hca_clss* cl2, MMsMatrix* dst){
int			i,j;
double		d,dx=-1;
hca_clss	**cmp1,**cmp2;
int			ncmp1,ncmp2;

	if(cl1->ref==0){
		cmp1=cl1->cmp;
		ncmp1=cl1->ncmp;
	}
	else{
		cmp1=NULL;
		cmp1=hca_stackadd(cmp1,cl1,0);
		ncmp1=1;
	}
	if(cl2->ref==0){
		cmp2=cl2->cmp;
		ncmp2=cl2->ncmp;
	}
	else{
		cmp2=NULL;
		cmp2=hca_stackadd(cmp2,cl2,0);
		ncmp2=1;
	}
	
	for(i=0;i<ncmp1;i++){
		for(j=0;j<ncmp2;j++){
			d=MMsGetDouble(dst,cmp1[i]->ref,cmp2[j]->ref);
			if(d>dx){
				dx=d;
			}
		}
	}

	if(cl1->ref!=0){
		free(cmp1);
	}
	if(cl2->ref!=0){
		free(cmp2);
	}
	
	return(dx);
}

//----------------------------------------------------------------------------
// 
//------------
static double min_dst(hca_clss* cl1, hca_clss* cl2, MMsMatrix* dst){
int			i,j;
double		d,dx=LONG_MAX;
hca_clss	**cmp1,**cmp2;
int			ncmp1,ncmp2;

	if(cl1->ref==0){
		cmp1=cl1->cmp;
		ncmp1=cl1->ncmp;
	}
	else{
		cmp1=NULL;
		cmp1=hca_stackadd(cmp1,cl1,0);
		ncmp1=1;
	}
	if(cl2->ref==0){
		cmp2=cl2->cmp;
		ncmp2=cl2->ncmp;
	}
	else{
		cmp2=NULL;
		cmp2=hca_stackadd(cmp2,cl2,0);
		ncmp2=1;
	}

	for(i=0;i<ncmp1;i++){
		for(j=0;j<ncmp2;j++){
			d=MMsGetDouble(dst,cmp1[i]->ref,cmp2[j]->ref);
			if(d<dx){
				dx=d;
			}
		}
	}

	if(cl1->ref!=0){
		free(cmp1);
	}
	if(cl2->ref!=0){
		free(cmp2);
	}
	return(dx);
}

//----------------------------------------------------------------------------
// 
//------------
static double ave_dst(hca_clss* cl1, hca_clss* cl2, MMsMatrix* dst){
int			i,j;
double		d=0,n;
hca_clss	**cmp1,**cmp2;
int			ncmp1,ncmp2;

	if(cl1->ref==0){
		cmp1=cl1->cmp;
		ncmp1=cl1->ncmp;
	}
	else{
		cmp1=NULL;
		cmp1=hca_stackadd(cmp1,cl1,0);
		ncmp1=1;
	}
	if(cl2->ref==0){
		cmp2=cl2->cmp;
		ncmp2=cl2->ncmp;
	}
	else{
		cmp2=NULL;
		cmp2=hca_stackadd(cmp2,cl2,0);
		ncmp2=1;
	}

	for(i=0;i<ncmp1;i++){
		for(j=0;j<ncmp2;j++){
			d+=MMsGetDouble(dst,cmp1[i]->ref,cmp2[j]->ref);
		}
	}

	if(cl1->ref!=0){
		free(cmp1);
	}
	if(cl2->ref!=0){
		free(cmp2);
	}
	n=ncmp1*ncmp2;
	return(d/n);
}

//----------------------------------------------------------------------------
// Tree sort
//------------
static void hca_treepermute(hca_clss* tree){
hca_clss**	clsss=NULL;
hca_clss**	rev=NULL;
hca_clss*	clss;
int			n=0,nrev=0;
int			score=0;

	clsss=hca_stackadd(clsss,tree,0);
	n++;
	rev=hca_stackadd(rev,tree,0);
	nrev++;
	do{
		clss=clsss[n-1];
		clsss=hca_stackrmv(clsss,n,n);
		n--;
		if(clss->ref!=0){
			score++;
			clss->cidx=score;
			clss->midwidth=score;
		}
		if(clss->left){
			clsss=hca_stackadd(clsss,clss->left,n);
			n++;
			rev=hca_stackadd(rev,clss->left,nrev);
			nrev++;
		}
		if(clss->right){
			clsss=hca_stackadd(clsss,clss->right,n);
			n++;
			rev=hca_stackadd(rev,clss->right,nrev);
			nrev++;
		}
	}
	while(n>=1);

	for(n=nrev-1;n>=0;n--){
		clss=rev[n];
		if((clss->left)&&(clss->right)){
			clss->midwidth=(clss->left->midwidth+clss->right->midwidth)/2.0;
		}				
	}
	free(rev);

} 

//----------------------------------------------------------------------------
// 
//------------
static hca_clss* hca_calc(MMsMatrix* dst, double(*compare)(hca_clss*,hca_clss*,MMsMatrix*)){
hca_clss**	clsss=hca_stackinit(dst->nl,0);
int			itr,i,j,n=dst->nl,k1=0,k2=0,idx=0,lid=dst->nl;
hca_clss	*cl1=NULL,*cl2,*bcl1,*bcl2;
double		d,dmin=LONG_MAX;

	for(itr=1;itr<=dst->nl-1;itr++){
		dmin=LONG_MAX;
		for(i=0;i<n-1;i++){
			cl1=clsss[i];
			for(j=i+1;j<n;j++){
				cl2=clsss[j];
				d=compare(cl1,cl2,dst);
				if(d<dmin){
					dmin=d;
					bcl1=cl1;
					bcl2=cl2;
					k1=i;
					k2=j;
				}
			}
		}
		if(k1<k2){
			clsss=hca_stackrmv(clsss,k1,n);
			clsss=hca_stackrmv(clsss,k2-1,n-1);
		}
		else{
			clsss=hca_stackrmv(clsss,k2,n);
			clsss=hca_stackrmv(clsss,k1-1,n-1);		
		}
		
		n-=2;
		lid++;
		cl1=hca_classnew(0,lid);
		
		if(bcl1->height<bcl2->height){
			cl1->left=bcl1;
			cl1->right=bcl2;
		}
		else{
			cl1->left=bcl2;
			cl1->right=bcl1;
		}
		
		
		cl1->height=dmin;
		bcl1->sup=cl1;
		if(bcl1->ref!=0){
			idx++;
			bcl1->cidx=idx;
			bcl1->midwidth=idx;
		}
		bcl2->sup=cl1;
		if(bcl2->ref!=0){
			idx++;
			bcl2->cidx=idx;
			bcl2->midwidth=idx;
		}
		cl1->midwidth=(bcl1->midwidth+bcl2->midwidth)/2.0;

		reports(cl1,cl1->left);
		reports(cl1,cl1->right);
		clsss=hca_stackadd(clsss,cl1,n);
		
		n++;		
	}
	
	hca_treepermute(cl1);
	
	return(cl1);
} 

//----------------------------------------------------------------------------
// 
//------------
hca_clss* hca_diameter(MMsMatrix* dst){
	return(hca_calc(dst,max_dst));
} 

//----------------------------------------------------------------------------
// 
//------------
hca_clss* hca_min(MMsMatrix* dst){
	return(hca_calc(dst,min_dst));
} 

//----------------------------------------------------------------------------
// 
//------------
hca_clss* hca_average(MMsMatrix* dst){
	return(hca_calc(dst,ave_dst));
} 

//----------------------------------------------------------------------------
// 
//------------
void hca_iterate(hca_clss* tree, void(*proc)(hca_clss*,void*), void* prm){
hca_clss**	clsss=NULL;
hca_clss*	clss;
int			n=0;
	clsss=hca_stackadd(clsss,tree,0);
	n++;
	do{
		clss=clsss[n-1];
		clsss=hca_stackrmv(clsss,n,n);
		n--;
		proc(clss,prm);
		if(clss->left){
			clsss=hca_stackadd(clsss,clss->left,n);
			n++;
		}
		if(clss->right){
			clsss=hca_stackadd(clsss,clss->right,n);
			n++;
		}
	}
	while(n>=1);
} 

//----------------------------------------------------------------------------
// 
//------------
void hca_classfree(hca_clss* tree){
hca_clss**	clsss=NULL;
hca_clss*	clss;
int			n=0;
	clsss=hca_stackadd(clsss,tree,0);
	n++;
	do{
		clss=clsss[n-1];
		clsss=hca_stackrmv(clsss,n,n);
		n--;
		if(clss->left){
			clsss=hca_stackadd(clsss,clss->left,n);
			n++;
		}
		if(clss->right){
			clsss=hca_stackadd(clsss,clss->right,n);
			n++;
		}
		if(clss->cmp){
			free(clss->cmp);
		}
		free(clss);
	}
	while(n>=1);
} 

#pragma mark -> IRA
//----------------------------------------------------------------------------
// 
//------------
hca_clss** ira_calc(mmx_ira* data){
	if(!data->center){
		return(NULL);
	}
MMsMatrix*	dt=CentreesReduites(data->data);
MMsMatrix*	ct=CentreesReduites(data->center);
hca_clss**	clsss=hca_stackinit(data->center->nl,0);
hca_clss**	clssso=hca_stackinit(data->data->nl,data->center->nl);
hca_clss*	clss;
double		e,d,dtot,z;
int			i,j,k;

	do{
		dtot=0;
// init
		for(i=0;i<data->center->nl;i++){
			if(clsss[i]->cmp){
				free(clsss[i]->cmp);
				clsss[i]->cmp=NULL;
			}
			clsss[i]->ncmp=0;
		}
		
// put characters to centers
		for(i=1;i<=dt->nl;i++){
			e=LONG_MAX;
			clss=NULL;
			for(j=1;j<=ct->nl;j++){
				d=0.0;
				for(k=1;k<=dt->nc;k++){
					d+=sqr(MMsGetDouble(dt,i,k)-MMsGetDouble(ct,j,k));
				}
				d=sqrt(d);
				if(d<e){
					clss=clsss[j-1];
					e=d;
				}
			}
			clss->cmp=hca_stackadd(clss->cmp,clssso[i-1],clss->ncmp);
			clss->ncmp++;
		}
		
// compute new positions
		for(i=1;i<=ct->nl;i++){
			z=clsss[i-1]->ncmp;
			d=0;
			for(j=1;j<=ct->nc;j++){
				e=0;
				for(k=1;k<=clsss[i-1]->ncmp;k++){
					e+=MMsGetDouble(dt,clsss[i-1]->cmp[k-1]->ref,j);
				}
				e=(clsss[i-1]->ncmp==0)?0:(e/z);
				d+=sqr(MMsGetDouble(ct,i,j)-e);
				MMsSetDouble(ct,i,j,e);
			}
			dtot+=sqrt(d);
		}
		
	}
	while(dtot>0);

// compute centers values
	for(i=1;i<=ct->nc;i++){  
		d=MMsGetDouble(data->indic,_kAverageIndex,i);
		e=MMsGetDouble(data->indic,_kStdDeviationIndex,i);
		for(j=1;j<=ct->nl;j++){
			MMsSetDouble(data->center,j,i,MMsGetDouble(ct,j,i)*e+d);
		}
	}
	
// put group numbers
	for(i=0;i<data->center->nl;i++){
		for(j=0;j<clsss[i]->ncmp;j++){
			clsss[i]->cmp[j]->sup=clsss[i];
		}
	}

	return(clsss);
};

//----------------------------------------------------------------------------
// 
//------------
void ira_classfree(hca_clss* tree){
int	i;
	for(i=0;i<tree->ncmp;i++){
		free(tree->cmp[i]);
	}
	free(tree);
} 
