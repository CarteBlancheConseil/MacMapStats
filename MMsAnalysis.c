//----------------------------------------------------------------------------
// File : MMsAnalysis.c
// Project : MacMapStats
// Purpose : C source file : Statistical analysis
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

#include "MMsAnalysis.h"
#include "MMsMatrix_Utils.h"

//----------------------------------------------------------------------------
// 
//------------
void uni_analysis(MMsMatrix *indata, mmx_uni* outdata){
int			i;
double		e,b;
double		min;
double		max;
double		ext;
double		sum;
double		ave;
double		med;
double		stddev;
double		var;
double		cvar;
double		cskew;
double		cpearson;
double		cfisher;

	outdata->kind=_kAnalysisUNI;
	outdata->data=get_matrix(indata,1);
	outdata->indic=MMsMatrixAlloc(_kcFisherIndex,1);
	outdata->clss=NULL;
	sum=0;
	for(i=1;i<=outdata->data->nl;i++){
		sum+=MMsGetDouble(outdata->data,i,1);
	}
	min=MMsGetDouble(outdata->data,1,1);
	max=MMsGetDouble(outdata->data,outdata->data->nl,1);
	ext=max-min;
	ave=sum/(double)outdata->data->nl;
	i=outdata->data->nl/2;
	med=(	(outdata->data->nl&1)!=0)? 
			MMsGetDouble(outdata->data,i,1): 
			(MMsGetDouble(outdata->data,i,1)+MMsGetDouble(outdata->data,i+1,1))/(double)2.0;
	var=0;
	cskew=0;
	cpearson=0;
	for(i=1;i<=outdata->data->nl;i++){
		e=MMsGetDouble(outdata->data,i,1)-ave;
		var+=(e*e);
		cskew+=(e*e*e);
		cpearson+=(e*e*e*e);
	}
//	e=var/(double)(outdata->data->nl-1);
	var/=(double)outdata->data->nl;
	stddev=sqrt(var);
	b=(double)outdata->data->nl*stddev*stddev*stddev;
	cskew/=b;
	b*=stddev;
	cpearson/=b;
	cfisher=cpearson-3.0;

	cvar=sqrt(var)/ave;
	
	MMsSetDouble(outdata->indic,_kMinIndex,1,min);
	MMsSetDouble(outdata->indic,_kMaxIndex,1,max);
	MMsSetDouble(outdata->indic,_kExtendIndex,1,ext);
	MMsSetDouble(outdata->indic,_kSumIndex,1,sum);
	MMsSetDouble(outdata->indic,_kAverageIndex,1,ave);
	MMsSetDouble(outdata->indic,_kMedianIndex,1,med);
	MMsSetDouble(outdata->indic,_kStdDeviationIndex,1,stddev);
	MMsSetDouble(outdata->indic,_kVarianceIndex,1,var);
	MMsSetDouble(outdata->indic,_kcVariationIndex,1,cvar);
	MMsSetDouble(outdata->indic,_kcSkewnessIndex,1,cskew);
	MMsSetDouble(outdata->indic,_kcPearsonIndex,1,cpearson);
	MMsSetDouble(outdata->indic,_kcFisherIndex,1,cfisher);
}

//----------------------------------------------------------------------------
// 
//------------
void biv_analysis(MMsMatrix *indata, mmx_biv* outdata){
int			i,j;
mmx_uni		uni;
MMsMatrix*	mx;

	outdata->kind=_kAnalysisBIV;
	outdata->data=MMsMatrixClone(indata);
	outdata->indic=MMsMatrixAlloc(_kCovarianceIndex,2);
	for(j=1;j<=2;j++){
		mx=MMsCloneColumn(indata,j);
		uni_analysis(mx,&uni);
		MMsMatrixFree(mx);
		MMsMatrixFree(uni.data);
		for(i=_kMinIndex;i<=_kcFisherIndex;i++){
			MMsSetDouble(outdata->indic,i,j,MMsGetDouble(uni.indic,i,1));
		}
		MMsMatrixFree(uni.indic);
	}
	MMsSetDouble(outdata->indic,_kcCorrelationIndex,1,MMsCorrelation(1,2,indata));
	MMsSetDouble(outdata->indic,_kcCorrelationIndex,2,MMsGetDouble(outdata->indic,_kcCorrelationIndex,1));
	MMsSetDouble(outdata->indic,_kCovarianceIndex,1,MMsCovariance(1,2,indata));
	MMsSetDouble(outdata->indic,_kCovarianceIndex,2,MMsGetDouble(outdata->indic,_kCovarianceIndex,1));
}

//----------------------------------------------------------------------------
// 
//------------
void pcn_analysis(	MMsMatrix *indata, 
					mmx_pca* outdata){
int			i,j;
mmx_uni		uni;
MMsMatrix	*mx,*vlp,*vcp;

	outdata->kind=_kAnalysisPCn;
	outdata->data=MMsMatrixClone(indata);
	outdata->indic=MMsMatrixAlloc(_kcFisherIndex,indata->nc);
	for(j=1;j<=indata->nc;j++){
		mx=MMsCloneColumn(indata,j);
		uni_analysis(mx,&uni);
		MMsMatrixFree(mx);
		MMsMatrixFree(uni.data);
		for(i=_kMinIndex;i<=_kcFisherIndex;i++){
			MMsSetDouble(outdata->indic,i,j,MMsGetDouble(uni.indic,i,1));
		}
		MMsMatrixFree(uni.indic);
	}
	mx=Normees(outdata->data);
	outdata->cv=MMsMxCorrelation(mx);
	MMsMatrixFree(mx);
	mx=MMsMatrixClone(outdata->cv);
	MMsEigenVector(mx,&vcp,&vlp);
	MMsMatrixFree(mx);
	mx=CentreesReduites(outdata->data);
	outdata->chr=GetACPIndCoord(mx,vcp,vlp,1);
	outdata->var=GetACPVarCoord(vcp,vlp);
	MMsMatrixFree(mx);
	MMsMatrixFree(vlp);
	MMsMatrixFree(vcp);
}

//----------------------------------------------------------------------------
// 
//------------
void pc_analysis(	MMsMatrix *indata, 
					mmx_pca* outdata){
int			i,j;
mmx_uni		uni;
MMsMatrix	*mx,*vlp,*vcp;

	outdata->kind=_kAnalysisPC;
	outdata->data=MMsMatrixClone(indata);
	outdata->indic=MMsMatrixAlloc(_kcFisherIndex,indata->nc);
	for(j=1;j<=indata->nc;j++){
		mx=MMsCloneColumn(indata,j);
		uni_analysis(mx,&uni);
		MMsMatrixFree(mx);
		MMsMatrixFree(uni.data);
		for(i=_kMinIndex;i<=_kcFisherIndex;i++){
			MMsSetDouble(outdata->indic,i,j,MMsGetDouble(uni.indic,i,1));
		}
		MMsMatrixFree(uni.indic);
	}
	outdata->cv=MMsMxCovariance(outdata->data);
	mx=MMsMatrixClone(outdata->cv);
	MMsMultiply(mx,outdata->cv->nl);
	MMsEigenVector(mx,&vcp,&vlp);
	MMsMatrixFree(mx);
	mx=Centrees(outdata->data);
	outdata->chr=GetACPIndCoord(mx,vcp,vlp,0);
	outdata->var=GetACPVarCoord(vcp,vlp);
	MMsMatrixFree(mx);
	MMsMatrixFree(vlp);
	MMsMatrixFree(vcp);
}

//----------------------------------------------------------------------------
// 
//------------
void f_analysis(MMsMatrix *indata, 
				mmx_fa* outdata){
int			i,j;
mmx_uni		uni;
MMsMatrix	*mx,*vlp,*vcp,*pi,*pj;
double		s;

	outdata->kind=_kAnalysisFc;
	outdata->data=MMsMatrixClone(indata);
	outdata->indic=MMsMatrixAlloc(_kcFisherIndex,indata->nc);
	for(j=1;j<=indata->nc;j++){
		mx=MMsCloneColumn(indata,j);
		uni_analysis(mx,&uni);
		MMsMatrixFree(mx);
		MMsMatrixFree(uni.data);
		for(i=_kMinIndex;i<=_kcFisherIndex;i++){
			MMsSetDouble(outdata->indic,i,j,MMsGetDouble(uni.indic,i,1));
		}
		MMsMatrixFree(uni.indic);
	}
	Marginales(outdata->data,&pi,&pj,&s);
	mx=MatriceS(outdata->data);
	MMsEigenVector(mx,&vcp,&vlp);
	MMsMatrixFree(mx);
	RmvVTrivial(&vcp);
	mx=NormAFC(outdata->data,pi,pj,s);
	outdata->chr=GetAFCIndCoord(mx,vcp);
	outdata->ctr=GetAFCIndContrib(outdata->chr,vlp);
	outdata->var=GetAFCVarCoord(vcp,vlp,pj);
	MMsMatrixFree(mx);
	MMsMatrixFree(vcp);
	MMsMatrixFree(vlp);
	MMsMatrixFree(pi);
	MMsMatrixFree(pj);
}

//----------------------------------------------------------------------------
// 
//------------
void hc_analysis(	MMsMatrix *indata, 
					mmx_hc* outdata){
int			i,j;
mmx_uni		uni;
MMsMatrix	*mx;

	outdata->kind=_kAnalysisHC;
	outdata->data=MMsMatrixClone(indata);
	outdata->indic=MMsMatrixAlloc(_kcFisherIndex,indata->nc);
	for(j=1;j<=indata->nc;j++){
		mx=MMsCloneColumn(indata,j);
		uni_analysis(mx,&uni);
		MMsMatrixFree(mx);
		MMsMatrixFree(uni.data);
		for(i=_kMinIndex;i<=_kcFisherIndex;i++){
			MMsSetDouble(outdata->indic,i,j,MMsGetDouble(uni.indic,i,1));
		}
		MMsMatrixFree(uni.indic);
	}
	mx=CentreesReduites(outdata->data);
	outdata->dst=GetDstIndMat(mx);
	MMsMatrixFree(mx);
}

//----------------------------------------------------------------------------
// 
//------------
void ira_analysis(	MMsMatrix *indata, 
					mmx_ira* outdata){
int			i,j;
mmx_uni		uni;
MMsMatrix	*mx;

	outdata->kind=_kAnalysisIRA;
	outdata->data=MMsMatrixClone(indata);
	outdata->indic=MMsMatrixAlloc(_kcFisherIndex,indata->nc);
	outdata->center=NULL;
	for(j=1;j<=indata->nc;j++){
		mx=MMsCloneColumn(indata,j);
		uni_analysis(mx,&uni);
		MMsMatrixFree(mx);
		MMsMatrixFree(uni.data);
		for(i=_kMinIndex;i<=_kcFisherIndex;i++){
			MMsSetDouble(outdata->indic,i,j,MMsGetDouble(uni.indic,i,1));
		}
		MMsMatrixFree(uni.indic);
	}
}

//----------------------------------------------------------------------------
// 
//------------
void mlt_analysis(MMsMatrix *indata, mmx_mlt* outdata){
int			i,j;
mmx_uni		uni;
MMsMatrix*	mx;
	
	outdata->kind=_kAnalysisMLT;
	outdata->data=MMsMatrixClone(indata);
	outdata->indic=MMsMatrixAlloc(_kcFisherIndex,outdata->data->nc);
	for(j=1;j<=outdata->data->nc;j++){
fprintf(stderr,"variable %d\n",j);
fflush(stderr);
		mx=MMsCloneColumn(indata,j);
		uni_analysis(mx,&uni);
		MMsMatrixFree(mx);
		MMsMatrixFree(uni.data);
		for(i=_kMinIndex;i<=_kcFisherIndex;i++){
fprintf(stderr,"indic %d\n",i);
fflush(stderr);
			MMsSetDouble(outdata->indic,i,j,MMsGetDouble(uni.indic,i,1));
		}
		MMsMatrixFree(uni.indic);
	}
}

//----------------------------------------------------------------------------
// 
//------------
void init_analysis(	mmx_analysis* data,
					int kind){
	data->kind=kind;
	data->data=NULL;
	data->indic=NULL;
	switch(kind){
		case _kAnalysisUndef:
			break;
		case _kAnalysisUNI:
			((mmx_uni*)data)->clss=NULL;
			break;
		case _kAnalysisBIV:
		case _kAnalysisMLT:
			break;
		case _kAnalysisPC:
		case _kAnalysisPCn:
			((mmx_pca*)data)->cv=NULL;
			((mmx_pca*)data)->chr=NULL;
			((mmx_pca*)data)->var=NULL;
			break;
		case _kAnalysisFc:
			((mmx_fa*)data)->ctr=NULL;
			((mmx_fa*)data)->chr=NULL;
			((mmx_fa*)data)->var=NULL;
			break;
		case _kAnalysisHC:
			((mmx_hc*)data)->dst=NULL;
			break;
		case _kAnalysisIRA:
			((mmx_ira*)data)->center=NULL;
			break;
		default:
			break;
	}
}

//----------------------------------------------------------------------------
// 
//------------
void free_analysis(	mmx_analysis* data){
fprintf(stderr,"Enter free_analysis\n");
fflush(stderr);
fprintf(stderr,"Free data\n");
fflush(stderr);
	MMsMatrixFree(((mmx_uni*)data)->data);
	((mmx_uni*)data)->data=NULL;
fprintf(stderr,"Free indic\n");
fflush(stderr);
	MMsMatrixFree(((mmx_uni*)data)->indic);
	((mmx_uni*)data)->indic=NULL;
	switch(data->kind){
		case _kAnalysisUndef:
			break;
		case _kAnalysisUNI:
fprintf(stderr,"_kAnalysisUNI:Free clss\n");
fflush(stderr);
			MMsMatrixFree(((mmx_uni*)data)->clss);
			((mmx_uni*)data)->clss=NULL;
			break;
		case _kAnalysisBIV:
		case _kAnalysisMLT:
			break;
		case _kAnalysisPC:
		case _kAnalysisPCn:
fprintf(stderr,"_kAnalysisPC/_kAnalysisPCn:Free cv\n");
fflush(stderr);
			MMsMatrixFree(((mmx_pca*)data)->cv);
			((mmx_pca*)data)->cv=NULL;
fprintf(stderr,"_kAnalysisPC/_kAnalysisPCn:Free chr\n");
fflush(stderr);
			MMsMatrixFree(((mmx_pca*)data)->chr);
			((mmx_pca*)data)->chr=NULL;
fprintf(stderr,"_kAnalysisPC/_kAnalysisPCn:Free var\n");
fflush(stderr);
			MMsMatrixFree(((mmx_pca*)data)->var);
			((mmx_pca*)data)->var=NULL;
			break;
		case _kAnalysisFc:
fprintf(stderr,"_kAnalysisFc:Free ctr\n");
fflush(stderr);
			MMsMatrixFree(((mmx_fa*)data)->ctr);
			((mmx_fa*)data)->ctr=NULL;
fprintf(stderr,"_kAnalysisFc:Free chr\n");
fflush(stderr);
			MMsMatrixFree(((mmx_fa*)data)->chr);
			((mmx_fa*)data)->chr=NULL;
fprintf(stderr,"_kAnalysisFc:Free var\n");
fflush(stderr);
			MMsMatrixFree(((mmx_fa*)data)->var);
			((mmx_fa*)data)->var=NULL;
			break;
		case _kAnalysisHC:
fprintf(stderr,"_kAnalysisHC:Free dst\n");
fflush(stderr);
			MMsMatrixFree(((mmx_hc*)data)->dst);
			((mmx_hc*)data)->dst=NULL;
			break;
		case _kAnalysisIRA:
fprintf(stderr,"_kAnalysisIRA:Free center\n");
fflush(stderr);
			MMsMatrixFree(((mmx_ira*)data)->center);
			((mmx_ira*)data)->center=NULL;
			break;
		default:
			break;
	}
fprintf(stderr,"Leave free_analysis\n");
fflush(stderr);
}

//----------------------------------------------------------------------------
// 
//------------
void linear_reg(mmx_biv* ana,
				double* a,
				double* b){
	*a=MMsGetDouble(ana->indic,_kCovarianceIndex,1)/MMsGetDouble(ana->indic,_kVarianceIndex,1);
	*b=(MMsGetDouble(ana->indic,_kAverageIndex,2)-(*a)*MMsGetDouble(ana->indic,_kAverageIndex,1));
}

//----------------------------------------------------------------------------
// 
//------------
void logarithmic_reg(	mmx_biv* ana,
						double* a,
						double* b){
mmx_biv		lana;
MMsMatrix*	mx=MMsMatrixClone(ana->data);
int			i;
	for(i=1;i<=mx->nl;i++){
		MMsSetDouble(mx,i,1,log(MMsGetDouble(mx,i,1)));
	}
	biv_analysis(mx,&lana);
	MMsMatrixFree(mx);
	linear_reg(&lana,a,b);
	MMsMatrixFree(lana.data);
	MMsMatrixFree(lana.indic);
}
