//----------------------------------------------------------------------------
// File : MacMapStats.c
// Project : MacMapStats
// Purpose : C source file : Statitics lib
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

#include "MacMapStats.h"

//----------------------------------------------------------------------------
// 
//------------
void dMMsEigenVector(	MMsMatrix* A, 
						MMsMatrix* *R, 
						MMsMatrix* *V);

//----------------------------------------------------------------------------
// 
//------------
static int MMsCompDouble(	const void* a,
							const void* b){
double*	oa=(double*)a;
double*	ob=(double*)b;
	if((*oa)>(*ob)){
		return(1);
	}
	if((*oa)<(*ob)){
		return(-1);
	}
	return(0);
}

//----------------------------------------------------------------------------
// 
//------------
static float MMsGetFloat(	MMsMatrix* mx, 
							int i, 
							int j){	
	return(MMsGetDouble(mx,i,j));
}

//----------------------------------------------------------------------------
// 
//------------
static void MMsSetFloat(	MMsMatrix* mx, 
							int i, 
							int j,
							float f){
	MMsSetDouble(mx,i,j,f);
}

//----------------------------------------------------------------------------
// 
//------------
static void MMsEigenVectorSort(	MMsMatrix* A, 
								MMsMatrix* R, 
								MMsMatrix* V){
long	i,j,k;
double	x;

	for(i=1;i<=A->nl;i++){
		for(j=i;j<=A->nl;j++){		
			if(MMsGetDouble(A,i,i)<MMsGetDouble(A,j,j)){
				x=MMsGetDouble(A,i,i);
				MMsSetDouble(A,i,i,MMsGetDouble(A,j,j)); 
				MMsSetDouble(A,j,j,x);
				for(k=1;k<=A->nl;k++){
					x=MMsGetDouble(R,i,k);
					MMsSetDouble(R,i,k,MMsGetDouble(R,j,k));
					MMsSetDouble(R,j,k,x);
				}
			}
		}
	}
	for(k=1;k<=A->nl;k++){
		MMsSetDouble(V,k,1,MMsGetDouble(A,k,k));
	}
}

//----------------------------------------------------------------------------
// 
//------------
static char MMsEigenVectorTestEnd(	int Lns, 
									int *l, 
									int *m, 
									int *Ind, 
									float *Thr, 
									float Ni){
char	end=0;

//fprintf(stderr,"MMsEigenVectorTestEnd->");
 	if((*m)!=Lns){
//fprintf(stderr,"(*m)!=Lns\n");
		(*m)++;
	}
	else{
		if((*l)!=Lns-1){
//fprintf(stderr,"((*l)!=Lns-1)\n");
			(*l)++;
			(*m)=(*l)+1;
		}
		else{
			if((*Ind)==1){
//fprintf(stderr,"((*Ind)==1)\n");
				(*Ind)=0;
				(*l)=1;
				(*m)=2;
			}
			else{
//fprintf(stderr,"Thr=*f (%f needed)\n",*Thr,(Ni*1.0E-6/(double)Lns));
				if(*Thr<=(Ni*1.0E-6/(float)Lns)){
					end=1;
				}
				else{
					(*Thr)/=(float)Lns;
				}
			}
		}
	}
 	return(end);
}

//----------------------------------------------------------------------------
// 
//------------
static void MMsEigenVectorCalcVec(	MMsMatrix* R, 
									int i, 
									int l, 
									int m, 
									float SinX, 
									float CosX){
float x;
	
	x=MMsGetFloat(R,l,i)*CosX-MMsGetFloat(R,m,i)*SinX;
	MMsSetFloat(R,m,i,(MMsGetFloat(R,l,i)*SinX+MMsGetFloat(R,m,i)*CosX));
	MMsSetFloat(R,l,i,x);
}

//----------------------------------------------------------------------------
// 
//------------
MMsMatrix* MMsMatrixAlloc(	int ni,
							int nj){
MMsMatrix*	mx;
size_t		sz=sizeof(MMsMatrix)+((ni*nj)*sizeof(double));

	mx=(MMsMatrix*)malloc(sz);
	if(mx==NULL){
		return(mx);
	}
	memset(mx,0,sz);
	mx->flg=1;	
	mx->sign=kMatrixRectangularKind_;
	mx->nl=ni;
	mx->nc=nj;
	return(mx);
}

//----------------------------------------------------------------------------
// 
//------------
MMsMatrix* MMsDiagMatrixAlloc(	int ni){
MMsMatrix*	mx;
size_t		sz=sizeof(MMsMatrix)+(((((ni*ni)-ni)/2l)+ni)*sizeof(double));
	mx=(MMsMatrix*)malloc(sz);
	if(mx==NULL){
		return(mx);
	}
	memset(mx,0,sz);
	mx->flg=1;	
	mx->sign=kMatrixDiagonalKind_;
	mx->nl=ni;
	mx->nc=ni;
	return(mx);
}

//----------------------------------------------------------------------------
// 
//------------
MMsMatrix* MMsMatrixClone(	MMsMatrix* mx){
MMsMatrix*	out=NULL;

	switch(mx->sign){
		case kMatrixRectangularKind_:
			out=MMsMatrixAlloc(mx->nl,mx->nc);
			memmove(out->f,mx->f,(mx->nl*mx->nc)*sizeof(double));
			break;
		case kMatrixDiagonalKind_:
			out=MMsDiagMatrixAlloc(mx->nl);
			memmove(out->f,mx->f,(((((mx->nl*mx->nl)-mx->nl)/2l)+mx->nl)*sizeof(double)));
			break;
	}
	return(out);
}

//----------------------------------------------------------------------------
// 
//------------
MMsMatrix* MMsUnitMatrixAlloc(int n){
long		i,j;
MMsMatrix*	out=MMsMatrixAlloc(n,n);
 
	if(out==NULL){
		return(out);
	}
	for(i=1;i<=out->nl;i++){
		for(j=1;j<=out->nc;j++){
			MMsSetDouble(out,i,j,(i!=j)?0:1);
		}
	}			
	return(out);
}

//----------------------------------------------------------------------------
// 
//------------
void MMsMatrixFree(MMsMatrix* mx){
	if(mx){
		free(mx);
	}
}

//----------------------------------------------------------------------------
// 
//------------
MMsMatrix* MMsMatrixAddLine(MMsMatrix* mx){
MMsMatrix*	res=MMsMatrixAlloc(mx->nl+1,mx->nc);
	if(!res){
		return(NULL);
	}
long	i,j;
	for(i=1;i<=mx->nl;i++){
		for(j=1;j<=mx->nc;j++){
			MMsSetDouble(res,i,j,MMsGetDouble(mx,i,j));
		}
	}
	for(j=1;j<=mx->nc;j++){
		MMsSetDouble(res,i,j,0);
	}
	MMsMatrixFree(mx);
	return(res);
}

//----------------------------------------------------------------------------
// 
//------------
MMsMatrix* MMsMatrixInsertLine(	MMsMatrix* mx, 
								int idx){
	if(idx<=0){
		return(NULL);
	}
	if(idx>mx->nl){
		return(MMsMatrixAddLine(mx));
	}
MMsMatrix*	res=MMsMatrixAlloc(mx->nl+1,mx->nc);
long		i,j;
	for(i=1;i<idx;i++){
		for(j=1;j<=mx->nc;j++){
			MMsSetDouble(res,i,j,MMsGetDouble(mx,i,j));
		}
	}
	for(j=1;j<=mx->nc;j++){
		MMsSetDouble(res,idx,j,0);
	}
	for(i=idx;i<=mx->nl;i++){
		for(j=1;j<=mx->nc;j++){
			MMsSetDouble(res,i+1,j,MMsGetDouble(mx,i,j));
		}
	}	
	MMsMatrixFree(mx);
	return(res);
}

//----------------------------------------------------------------------------
// 
//------------
MMsMatrix* MMsMatrixRemoveLine(	MMsMatrix* mx, 
								int idx){
	if(idx<=0){
		return(NULL);
	}
	if(idx>mx->nl){
		idx=mx->nl;
	}
MMsMatrix*	res=MMsMatrixAlloc(mx->nl-1,mx->nc);
long		i,j;
	for(i=1;i<idx;i++){
		for(j=1;j<=mx->nc;j++){
			MMsSetDouble(res,i,j,MMsGetDouble(mx,i,j));
		}
	}
	for(i=idx+1;i<=mx->nl;i++){
		for(j=1;j<=mx->nc;j++){
			MMsSetDouble(res,i-1,j,MMsGetDouble(mx,i,j));
		}
	}	
	MMsMatrixFree(mx);
	return(res);
}

//----------------------------------------------------------------------------
// 
//------------
double MMsGetDouble(MMsMatrix* mx, 
					int i, 
					int j){	
	switch(mx->sign){
		case kMatrixRectangularKind_:
			return(mx->f[((i-1l)*mx->nc+j)-1l]);
			break;
		case kMatrixDiagonalKind_:
			return(	(i<j)?
					mx->f[(int)((((double)mx->nc)-((double)(i/2.0)))*((double)(i-1l))+(j-1l))]:
					mx->f[(int)((((double)mx->nc)-((double)(j/2.0)))*((double)(j-1l))+(i-1l))]);
			break;
	}
	return(nan(""));
}

//----------------------------------------------------------------------------
// 
//------------
void MMsSetDouble(	MMsMatrix* mx, 
					int i, 
					int j,
					double f){
	switch(mx->sign){
		case kMatrixRectangularKind_:
			mx->f[((i-1)*mx->nc)+(j-1)]=f;
			break;
		case kMatrixDiagonalKind_:
			if(i<j){
				mx->f[(int)((((double)mx->nc)-((double)(i/2.0)))*((double)(i-1l))+(j-1l))]=f;
			}
			else{
				mx->f[(int)((((double)mx->nc)-((double)(j/2.0)))*((double)(j-1l))+(i-1l))]=f;
			}
			break;
	}
}

//----------------------------------------------------------------------------
// 
//------------
MMsMatrix* MMsCloneColumn(	MMsMatrix* mx,
							int nj){
MMsMatrix*	mxr;
size_t		sz=sizeof(MMsMatrix)+((mx->nl)*sizeof(double));
long		i;
	mxr=(MMsMatrix*)malloc(sz);
	if(mxr==NULL){
		return(mxr);
	}
	memset(mxr,0,sz);
	mxr->flg=1;
	mxr->sign=kMatrixRectangularKind_;
	mxr->nl=mx->nl;
	mxr->nc=1;
	for(i=1;i<=mx->nl;i++){
		MMsSetDouble(mxr,i,1,MMsGetDouble(mx,i,nj));
	}
	return(mxr);
}

//----------------------------------------------------------------------------
// 
//------------
void MMsMatrixSort(	MMsMatrix* mx){
	switch(mx->sign){
		case kMatrixRectangularKind_:
			qsort(mx->f,mx->nl*mx->nc,sizeof(double),MMsCompDouble);
			break;
		case kMatrixDiagonalKind_:
			qsort(mx->f,(((mx->nl*mx->nl)-mx->nl)/2)+mx->nl,sizeof(double),MMsCompDouble);
			break;
	}
}

//----------------------------------------------------------------------------
// 
//------------
void MMsDumpTTxt(	MMsMatrix* mx,
					FILE* fp){
long	i,j;
	for(i=1;i<=mx->nl;i++){
		for(j=1;j<=mx->nc;j++){
			fprintf(fp,"%f%c",MMsGetDouble(mx,i,j),(j==mx->nc)?'\n':'\t');
		}
	}
}

//----------------------------------------------------------------------------
// 
//------------
double MMsMin(	int c, 
				MMsMatrix* mx){
long	i;
double	r=MMsGetDouble(mx,1,c); 
	
	for(i=2;i<=mx->nl;i++){
		if(MMsGetDouble(mx,i,c)<r){
			r=MMsGetDouble(mx,i,c);
		}
	}
	return(r);
}

//----------------------------------------------------------------------------
// 
//------------
double MMsMax(	int c, 
				MMsMatrix* mx){
long	i;
double	r=MMsGetDouble(mx,1,c); 
	
	for(i=2;i<=mx->nl;i++){
		if(MMsGetDouble(mx,i,c)>r){
			r=MMsGetDouble(mx,i,c);
		}
	}
	return(r);
}


//----------------------------------------------------------------------------
// 
//------------
double MMsAverage(	int c, 
					MMsMatrix* mx){
long	i;
double	r=0; 

	for(i=1;i<=mx->nl;i++){
		r+=MMsGetDouble(mx,i,c);
	}
	return(r/(double)mx->nl);
}

//----------------------------------------------------------------------------
// 
//------------
double MMsMedian(	int c, 
					MMsMatrix* mx){
MMsMatrix* mxr=MMsCloneColumn(mx,c);
	if(mxr==NULL){
		return(nan(""));
	}
	MMsMatrixSort(mxr);
double	e=	(mxr->nl&1)?
			MMsGetDouble(mxr,(mxr->nl/2)+1,1):
			(MMsGetDouble(mxr,(mxr->nl/2)+1,1)+MMsGetDouble(mxr,(mxr->nl/2),1))/2.0;
	MMsMatrixFree(mxr);
	return(e);
}

//----------------------------------------------------------------------------
// 
//------------
double MMsSum(	int c, 
				MMsMatrix* mx){
long	i;
double	r=0; 

	for(i=1;i<=mx->nl;i++){
		r+=MMsGetDouble(mx,i,c);
	}
	return(r);
}

//----------------------------------------------------------------------------
// 
//------------
double MMsVariance(	int c, 
					MMsMatrix* mx){
long	i;
double	r=0,e; 

	e=MMsAverage(c,mx);
	for(i=1;i<=mx->nl;i++){
		r+=sqr(MMsGetDouble(mx,i,c)-e);
	}
	return(r/(double)mx->nl);
}

//----------------------------------------------------------------------------
// 
//------------
double MMsStdDeviation(	int c, 
						MMsMatrix* mx){
long	i;
double	r=0,e; 

	e=MMsAverage(c,mx);
	for(i=1;i<=mx->nl;i++){
		r+=sqr(MMsGetDouble(mx,i,c)-e);
	}
	return(sqrt(r/(double)mx->nl));
}

//----------------------------------------------------------------------------
// 
//------------
double MMsSquariance(	int c, 
					MMsMatrix* mx){
long	i;
double	r=0,e; 

	e=MMsAverage(c,mx);
	for(i=1;i<=mx->nl;i++){
		r+=sqr(MMsGetDouble(mx,i,c)-e);
	}
	return(r);
}

//----------------------------------------------------------------------------
// 
//------------
double MMsCovariance(	int c1, 
						int c2, 
						MMsMatrix* mx){
long	i;
double	r=0; 
double	e1=MMsAverage(c1,mx);
double	e2=MMsAverage(c2,mx);

	for(i=1;i<=mx->nl;i++){
		r+=(MMsGetDouble(mx,i,c1)-e1)*(MMsGetDouble(mx,i,c2)-e2);
	}
	return(r/(double)mx->nl);
}

//----------------------------------------------------------------------------
// 
//------------
double MMsCorrelation(	int c1, 
						int c2, 
						MMsMatrix* mx){
	return(	MMsCovariance(c1,c2,mx)/
			((sqrt(MMsSquariance(c1,mx)/(double)mx->nl)*
			sqrt(MMsSquariance(c2,mx)/(double)mx->nl)))	);
}

//----------------------------------------------------------------------------
// 
//------------
void MMsMultiply(MMsMatrix* mx, double x){
long i,j;

	for(i=1;i<=mx->nl;i++){
		for(j=1;j<=mx->nc;j++){
			MMsSetDouble(mx,i,j,(MMsGetDouble(mx,i,j)*x));
		}
	}
}

//----------------------------------------------------------------------------
// 
//------------
MMsMatrix* MMsMxMultiply(	MMsMatrix* a, 
							MMsMatrix* b){
MMsMatrix*	out=MMsMatrixAlloc(a->nl,b->nc);
long		i,j,k;
 
	if(out==NULL){
		return(NULL);
	}
	for(i=1;i<=out->nl;i++){
		for(j=1;j<=out->nc;j++){
			for(k=1;k<=a->nc;k++){
				MMsSetDouble(out,i,j,(MMsGetDouble(out,i,j)+(MMsGetDouble(a,i,k)*MMsGetDouble(b,k,j))));
			}
		}
	}
	return(out);
}

//----------------------------------------------------------------------------
// 
//------------
MMsMatrix* MMsMxCorrelation(	MMsMatrix* mx){
MMsMatrix*	out=MMsMatrixAlloc(mx->nc,mx->nc);
long		i,j;
 
	if(out==NULL){
		return(NULL);
	}
	for(i=1;i<=mx->nc;i++){
		for(j=1;j<=mx->nc;j++){
			MMsSetDouble(out,i,j,MMsCorrelation(i,j,mx));
		}
	}
	return(out);
}

//----------------------------------------------------------------------------
// 
//------------
MMsMatrix* MMsMxCovariance(	MMsMatrix* mx){
MMsMatrix*	out=MMsMatrixAlloc(mx->nc,mx->nc);
long		i,j;
 
	if(out==NULL){
		return(NULL);
	}
	for(i=1;i<=mx->nc;i++){
		for (j=1;j<=mx->nc;j++){
			MMsSetDouble(out,i,j,MMsCovariance(i,j,mx));
		}
	}
	return(out);
}

//----------------------------------------------------------------------------
// 
//------------
MMsMatrix* MMsTranspose(	MMsMatrix* mx){
MMsMatrix*	out=MMsMatrixAlloc(mx->nc,mx->nl);
long		i,j;
 
	if(out==NULL){
		return(NULL);
	}
	for(i=1;i<=mx->nl;i++){
		for(j=1;j<=mx->nc;j++){
			MMsSetDouble(out,j,i,MMsGetDouble(mx,i,j)); 
		}
	}
	return(out);
}

//----------------------------------------------------------------------------
// 
//------------
void MMsEigenVector(	MMsMatrix* A, 
						MMsMatrix* *R, 
						MMsMatrix* *V){		
long	i,j;
int		l,m,Ind,Im1,Im2,Il1,Il2;
float	Ni,Thr,SinX,SinX2,CosX,CosX2,SinOS,X,Y;
char	end=0;

	Ni=0;
	Ind=0;
	*R=MMsUnitMatrixAlloc(A->nl);
	*V=MMsMatrixAlloc(A->nl,1);
	for(i=1;i<=A->nl;i++){
		for(j=i;j<=A->nl;j++){
			if(i!=j){		
				Ni+=sqr(MMsGetFloat(A,i,j));
			}
		}
	}
	if(Ni>0){
		Ni=sqrtf((float)2.*Ni);
		Thr=Ni/(float)A->nl;
		l=1;
		m=l+1;
		do{
	 	 	if(fabs(MMsGetFloat(A,m,l))<Thr){
				end=MMsEigenVectorTestEnd(A->nl,&l,&m,&Ind,&Thr,Ni);
			}
			else{
				Ind=1;
				X=(MMsGetFloat(A,l,l)-MMsGetFloat(A,m,m))/(float)2.;
				Y=-MMsGetFloat(A,m,l)/sqrtf(sqr(MMsGetFloat(A,m,l))+sqr(X));
				if(X<0){
					Y=-Y;
				}
				SinX=Y/sqrtf((float)2.*((float)1.+sqrtf((float)1.-sqr(Y))));
				SinX2=sqr(SinX);
				CosX=sqrtf((float)1.-SinX2);
				CosX2=sqr(CosX);
				SinOS=SinX*CosX;
				for(i=1;i<=A->nl;i++){
					if(i!=l){
			 			if(i<m){
							Im1=i;
							Im2=m;
						}
						else if(i>m){
							Im1=m;
							Im2=i;
						}
						else if(i==m){
							MMsEigenVectorCalcVec(*R,i,l,m,SinX,CosX);
							continue;
						}
						if(i<l){
							Il1=i;
							Il2=l;
						}
						else{
							Il1=l;
							Il2=i;
						}
						X=MMsGetFloat(A,Il2,Il1)*CosX-MMsGetFloat(A,Im2,Im1)*SinX;
					    MMsSetFloat(A,Im2,Im1,(MMsGetFloat(A,Il2,Il1)*SinX+MMsGetFloat(A,Im2,Im1)*CosX));
						MMsSetFloat(A,Il2,Il1,X);
					    MMsEigenVectorCalcVec(*R,i,l,m,SinX,CosX);
					}
					else{
						MMsEigenVectorCalcVec(*R,i,l,m,SinX,CosX);
					}
				}
				X=(float)2.*MMsGetFloat(A,m,l)*SinOS;
				Y=MMsGetFloat(A,l,l)*CosX2+MMsGetFloat(A,m,m)*SinX2-X;
				X=MMsGetFloat(A,l,l)*SinX2+MMsGetFloat(A,m,m)*CosX2+X;
				MMsSetFloat(A,m,l,((MMsGetFloat(A,l,l)-MMsGetFloat(A,m,m))*SinOS+MMsGetFloat(A,m,l)*(CosX2-SinX2)));
				MMsSetFloat(A,l,l,Y);
				MMsSetFloat(A,m,m,X);
				end=MMsEigenVectorTestEnd(A->nl,&l,&m,&Ind,&Thr,Ni);
			}
		}
		while(!end); 
	}
	MMsEigenVectorSort(A,*R,*V);
}

//----------------------------------------------------------------------------
// 
//------------
void Marginales(MMsMatrix* mx, MMsMatrix* *pi, MMsMatrix* *pj, double *s){
long	i,j;
double	e;
  
	(*s)=0;
	(*pi)=MMsMatrixAlloc(mx->nl,1);
	(*pj)=MMsMatrixAlloc(1,mx->nc);
	for(i=1;i<=mx->nl;i++){
		for(j=1;j<=mx->nc;j++){
			e=MMsGetDouble(mx,i,j),
			(*s)+=e;
			MMsSetDouble((*pi),i,1,MMsGetDouble((*pi),i,1)+e);
			MMsSetDouble((*pj),1,j,MMsGetDouble((*pj),1,j)+e);
		}
	}
	for(i=1;i<=mx->nl;i++){
		e=MMsGetDouble((*pi),i,1);
		if(e==0){
			e=0.000000000001;
		}
		MMsSetDouble((*pi),i,1,e/(*s));
	} 
	for(j=1;j<=mx->nc;j++){ 
		e=MMsGetDouble((*pj),1,j);
		if(e==0){
			e=0.000000000001;
		}
		MMsSetDouble((*pj),1,j,e/(*s));
	}
}

//----------------------------------------------------------------------------
// 
//------------
MMsMatrix* Centrees(MMsMatrix* mx){
MMsMatrix*	out=MMsMatrixAlloc(mx->nl,mx->nc);
long		i,j;
double		a;

	if(out==NULL){
		return(NULL);
	}
	for(i=1;i<=mx->nc;i++){  
		a=MMsAverage(i,mx);
		for(j=1;j<=mx->nl;j++){
			MMsSetDouble(out,j,i,(MMsGetDouble(mx,j,i)-a));
		}
	}
	return(out);
}

//----------------------------------------------------------------------------
// 
//------------
MMsMatrix* CentreesReduites(MMsMatrix* mx){
MMsMatrix*	out=MMsMatrixAlloc(mx->nl,mx->nc);
long		i,j;
double		e,a;
 
	if(out==NULL){
		return(NULL);
	}
	for(i=1;i<=mx->nc;i++){
		a=MMsAverage(i,mx);
		e=MMsStdDeviation(i,mx);
		for(j=1;j<=mx->nl;j++){
			MMsSetDouble(out,j,i,((MMsGetDouble(mx,j,i)-a)/e));
		}
	}
	return(out);
}

//----------------------------------------------------------------------------
// 
//------------
MMsMatrix* Normees (MMsMatrix* mx){
MMsMatrix*	out=MMsMatrixAlloc(mx->nl,mx->nc);
long		i,j;
double		e,a;
 
	if(out==NULL){
		return(NULL);
	}
	for(i=1;i<=mx->nc;i++){
		a=MMsAverage(i,mx);
		e=MMsStdDeviation(i,mx);
		for(j=1;j<=mx->nl;j++){
			MMsSetDouble(out,j,i,((MMsGetDouble(mx,j,i)-a)/(sqrt(out->nl)*e)));
		}
	} 
	return(out);
}

//----------------------------------------------------------------------------
// 
//------------
MMsMatrix* MatriceS(MMsMatrix* p){
long		i,j;
MMsMatrix	*pi,*pj,*r,*tr;
double		s;
	
	Marginales(p,&pi,&pj,&s);
	r=MMsMatrixAlloc(p->nl,p->nc);
	for(i=1;i<=r->nl;i++){
		for(j=1;j<=r->nc;j++){
			MMsSetDouble(r,i,j,MMsGetDouble(p,i,j)/(sqrt(MMsGetDouble(pj,1,j)*MMsGetDouble(pi,i,1)*s)));
		}
	}
	MMsMatrixFree(pi);
	MMsMatrixFree(pj);
	tr=MMsTranspose(r);
	pi=MMsMxMultiply(tr,r);
	MMsMatrixFree(r);
	MMsMatrixFree(tr);
	for(i=1;i<=pi->nl;i++){
		for(j=1;j<=pi->nc;j++){
			MMsSetDouble(pi,i,j,MMsGetDouble(pi,i,j)/s);
		}
	}
	return(pi);
}

//----------------------------------------------------------------------------
// 
//------------
MMsMatrix* NormAFC (MMsMatrix* p, MMsMatrix* pi, MMsMatrix* pj, float s){
long		i,j;
MMsMatrix	*out,*m;

	out=MMsMatrixAlloc(p->nl,p->nc);
	if(out==NULL){
		return(NULL);
	}
	m=MMsMatrixAlloc(pj->nl,pj->nc);
	if(m==NULL){
		MMsMatrixFree(out);
		return(NULL);
	}
	for(i=1;i<=p->nl;i++){
		for(j=1;j<=p->nc;j++){
			MMsSetDouble(out,i,j,(MMsGetDouble(p,i,j)/(s*(sqrt(MMsGetDouble(pj,1,j))*MMsGetDouble(pi,i,1)))));
			MMsSetDouble(m,1,j,MMsGetDouble(m,1,j)+MMsGetDouble(out,i,j)*MMsGetDouble(pi,i,1));
		}
	}
	for(i=1;i<=p->nl;i++){
		for(j=1;j<=p->nc;j++){
			MMsSetDouble(out,i,j,(MMsGetDouble(out,i,j)-MMsGetDouble(m,1,j)));
		}
	}
	MMsMatrixFree(m);
	return(out);
}

//----------------------------------------------------------------------------
// 
//------------
MMsMatrix* GetACPIndCoord(MMsMatrix* A, MMsMatrix* R, MMsMatrix* V, char Scaled){
long		i,j,k;
double		Scale;
MMsMatrix*	C;
	
	C=MMsMatrixAlloc(A->nl,A->nc);
	for(i=1;i<=A->nl;i++){
		for(j=1;j<=A->nc;j++){
			for(k=1;k<=A->nc;k++){
				MMsSetDouble(C,i,j,(MMsGetDouble(C,i,j)+(MMsGetDouble(A,i,k)*MMsGetDouble(R,j,k))));
			}
		}
	}
	if(Scaled){
		for(i=1;i<=A->nc;i++){
			Scale=sqrt(fabs(MMsGetDouble(V,i,1))*(double)A->nc);
			for(j=1;j<=A->nl;j++){
				MMsSetDouble(C,j,i,(MMsGetDouble(C,j,i)/Scale));
			}
		}
	}
	return(C);
}

//----------------------------------------------------------------------------
// 
//------------
MMsMatrix* GetACPVarCoord(MMsMatrix* R, MMsMatrix* V){
long		i,j;
MMsMatrix*	C;
double		s=0,p=100;

	C=MMsMatrixAlloc(R->nl*4+3,R->nc);
	if(C==NULL){
		return(NULL);
	}
	for(i=1;i<=C->nc;i++){
		MMsSetDouble(C,1,i,MMsGetDouble(V,i,1));
		s+=MMsGetDouble(V,i,1);
	}
	for(i=1;i<=C->nc;i++){
		MMsSetDouble(C,2,i,MMsGetDouble(V,i,1)/s*p);
	}	
	s=0;
	for(i=1;i<=C->nc;i++){
		s+=MMsGetDouble(C,2,i);
		MMsSetDouble(C,3,i,s);
	}
	
	for(i=1;i<=C->nc;i++){
		for(j=1;j<=C->nc;j++){
			MMsSetDouble(C,i+3,j,(MMsGetDouble(R,j,i)));
			MMsSetDouble(C,i+3+R->nc,j,(MMsGetDouble(C,i+3,j)*sqrt(fabs(MMsGetDouble(V,j,1)))));
			MMsSetDouble(C,i+3+R->nc*2,j,(sqr(MMsGetDouble(C,i+3+R->nc,j))));
			MMsSetDouble(C,i+3+R->nc*3,j,(MMsGetDouble(C,i+3+R->nc*2,j)/MMsGetDouble(V,j,1)*p));
		}
	}
	return(C);
}

//----------------------------------------------------------------------------
// 
//------------
MMsMatrix* GetAFCIndCoord (MMsMatrix* A, MMsMatrix* R){
	return(MMsMxMultiply(A,R));
}

//----------------------------------------------------------------------------
// 
//------------
MMsMatrix* GetAFCIndContrib (MMsMatrix* A, MMsMatrix* V){
MMsMatrix*	out=MMsMatrixAlloc(A->nl,A->nc);
long		i,j;

	if(out==NULL){
		return(NULL);
	}
	for(i=1;i<=out->nl;i++){
		for(j=1;j<=out->nc;j++){
			MMsSetDouble(out,i,j,MMsGetDouble(A,i,j)/sqrt(fabs(MMsGetDouble(V,j+1,1))));
		}
	}
	return(out);
}

//----------------------------------------------------------------------------
// 
//------------
MMsMatrix* GetAFCVarCoord(MMsMatrix* r, MMsMatrix* v, MMsMatrix* pj){
MMsMatrix*	out=MMsMatrixAlloc(r->nl*4+3,r->nc);
long		i,j;
double		s=0,p=100;

	if(out==NULL){
		return(NULL);
	}
//==> Valeurs propres
	for(j=1;j<=out->nc;j++){
		p=MMsGetDouble(v,j+1,1);
		MMsSetDouble(out,1,j,p);
		s+=p;
	}
//==> Pourcentage d'inertie
	for(j=1;j<=out->nc;j++){
		MMsSetDouble(out,2,j,MMsGetDouble(out,1,j)/s*100.);
	}
//==> Pourcentage d'inertie cumulé
	s=0;
	for(j=1;j<=out->nc;j++){
		s+=MMsGetDouble(out,2,j);
		MMsSetDouble(out,3,j,s);
	}
	for(i=4;i<=r->nl+3;i++){
		for(j=1;j<=out->nc;j++){
//==> Vecteurs propres
			MMsSetDouble(out,i,j,MMsGetDouble(r,i-3,j));
//==> Correlations
			MMsSetDouble(out,i+r->nl*2,j,(MMsGetDouble(out,i,j)/sqrt(MMsGetDouble(pj,1,i-3))));
//==> Coordonnées
			MMsSetDouble(out,i+r->nl,j,(MMsGetDouble(out,i+r->nl*2,j)*sqrt(fabs(MMsGetDouble(v,j+1,1)))));
//==> Contributions
			MMsSetDouble(out,i+r->nl*3,j,(sqr(MMsGetDouble(out,i+r->nl,j))*MMsGetDouble(pj,1,i-3)/MMsGetDouble(v,j+1,1)));
		}
	}
	return(out);
}

//----------------------------------------------------------------------------
// 
//------------
void RmvVTrivial(MMsMatrix* *r){
MMsMatrix*	out=MMsMatrixAlloc((*r)->nl,(*r)->nc-1);
long		i,j;

	if(out==NULL){
		return;
	}
	for(i=1;i<=out->nl;i++){
		for(j=1;j<=out->nc;j++){
			MMsSetDouble(out,i,j,MMsGetDouble((*r),j+1,i));
		}
	}
	MMsMatrixFree(*r);
	(*r)=out;
}

//----------------------------------------------------------------------------
// 
//------------
MMsMatrix* GetDstIndMat(MMsMatrix* mx){
long		i,j,l;
double		d;
MMsMatrix*	out=MMsDiagMatrixAlloc(mx->nl);

	if(out==NULL){
		return(NULL);
	}
	for(i=1;i<mx->nl;i++){
		for(l=i+1;l<=mx->nl;l++){
			d=0;
			for(j=1;j<=mx->nc;j++){
				d+=sqr(MMsGetDouble(mx,i,j)-MMsGetDouble(mx,l,j));
			}
			d=sqrt(d);
			MMsSetDouble(out,i,l,d);
		}
	}
	return(out);
}
