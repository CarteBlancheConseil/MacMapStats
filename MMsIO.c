//----------------------------------------------------------------------------
// File : MMsIO.c
// Project : MacMapStats
// Purpose : C source file : Classification & analysis read/write utils
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
// 10/10/2017 64 bits support.
//----------------------------------------------------------------------------

#include "MMsIO.h"
#include <pwd.h>
#include <unistd.h>
#include <sys/stat.h>

#define FileNotFoundError	-10000

//----------------------------------------------------------------------------
// 
//------------
typedef struct hca_parse{
	int			n;
	hca_clss**	lclsss;
}hca_parse;

//----------------------------------------------------------------------------
// endian swap proc
//------------
static void	swapword(int sz, void* word){
int				i;
unsigned char	buff;

	for(i=0;i<sz/2;i++){
		buff=((unsigned char*)word)[i];
		((unsigned char*)word)[i]=((unsigned char*)word)[sz-i-1];
		((unsigned char*)word)[sz-i-1]=buff;
	}
}

//----------------------------------------------------------------------------
// endian swap proc
//------------
static void	swap(int sz, int n, void* arr){
int	i;
    for(i=0;i<n;i++){
		swapword(sz,(void*)(((long)arr)+(i*sz)));
    }
}

//----------------------------------------------------------------------------
// 
//------------
static int cd(const char* path, char*	cwd){
	if(cwd&&!getcwd(cwd,PATH_MAX)){
		return(-1);
	}
	if(chdir(path)){
		return(-2);
	}
	return(0);
}

//----------------------------------------------------------------------------
// 
//------------
static int make_package(const char* path, const char* name, const char* ext){
mode_t	msk=S_IRWXU|S_IRWXG|S_IROTH|S_IXOTH;
char	hname[FILENAME_MAX];
char	cwd[PATH_MAX];
char	lext[64];
char*	d=strrchr(name,'.');
	
	if(!getcwd(cwd,PATH_MAX)){
		return(-1);
	}
	if(chdir(path)){
		chdir(cwd);
		return(-2);
	}
	
	if(strlen(ext)>0){
		sprintf(hname,"%s.%s",name,ext);
		sprintf(lext,"%s",ext);
	}
	else if(d!=NULL){
		sprintf(hname,"%s",name);
		d++;
		sprintf(lext,"%s",d);
	}
	else{
		sprintf(hname,"%s",name);
		sprintf(lext,"%s","???");
	}
	
	if(mkdir(hname,msk)){
	}
	
	if(chdir(hname)){
		chdir(cwd);
		return(-2);
	}
	if(mkdir("Contents",msk)){
	}
	if(chdir("Contents")){
		chdir(cwd);
		return(-2);
	}
	if(mkdir("Datas",msk)){
	}
	if(mkdir("Resources",msk)){
	}

FILE* fp=fopen("PkgInfo","w");
	fprintf(fp,"%s_MapS",lext);
	fclose(fp);
	chdir(cwd);

	return(0);
}

//----------------------------------------------------------------------------
// 
//------------
static int MMsIOKindRead(int* kind){
FILE*	fp=fopen("kind.knd","rb");
	if(!fp){
		return(FileNotFoundError);
	}
int	end;
	fread(&end,sizeof(int),1,fp);
	fread(kind,sizeof(int),1,fp);
	if(end!=1){
		swap(sizeof(int),1,kind);
	}
	fclose(fp);
	return(0);
}

//----------------------------------------------------------------------------
// 
//------------
static int MMsIOKindWrite(int kind){
FILE*	fp=fopen("kind.knd","wb");
	if(!fp){
		return(FileNotFoundError);
	}
int	end=1;
	fwrite(&end,sizeof(int),1,fp);
	fwrite(&kind,sizeof(int),1,fp);
	fclose(fp);
	return(0);
}

//----------------------------------------------------------------------------
// 
//------------
static int MMsIOMatrixRead(MMsMatrix** mx, FILE* fp){
size_t		sz=sizeof(MMsMatrix);
MMsMatrix	bf;
    *mx=NULL;
    if(fread(&bf,sz,1,fp)!=1){
fprintf(stderr,"MMsIOMatrixRead : fread (header) failed\n");
        return(-1);
    }
    fseek(fp,0,SEEK_SET);
    if(bf.flg!=1){
fprintf(stderr,"MMsIOMatrixRead : bad endian %d\n",bf.flg);
        swap(sizeof(int),3,&(bf.sign));
    }
    if(bf.sign==kMatrixRectangularKind_){
        (*mx)=MMsMatrixAlloc(bf.nl,bf.nc);
        sz+=(bf.nl*bf.nc*sizeof(double));
    }
    else if(bf.sign==kMatrixDiagonalKind_){
        (*mx)=MMsDiagMatrixAlloc(bf.nl);
        sz+=(((((bf.nl*bf.nl)-bf.nl)/2l)+bf.nl)*sizeof(double));
    }
    else{
fprintf(stderr,"MMsIOMatrixRead : bad matrix kind %d\n",bf.sign);
        return(-2);
    }
    if(fread(*mx,sz,1,fp)!=1){
fprintf(stderr,"MMsIOMatrixRead : fread (data) failed\n");
        return(-3);
    }
    if(bf.flg!=1){
        swap(sizeof(int),4,*mx);
        if(bf.sign==kMatrixRectangularKind_){
            swap(sizeof(double),(*mx)->nl*(*mx)->nc,(*mx)->f);
        }
        else if(bf.sign==kMatrixDiagonalKind_){
            swap(sizeof(double),(((((*mx)->nl*(*mx)->nl)-(*mx)->nl)/2l)+(*mx)->nl),(*mx)->f);
        }
    }
    return(0);
}

//----------------------------------------------------------------------------
// 
//------------
static int MMsIOMatrixReadFile(MMsMatrix** mx, const char* name){
int		status;
FILE*	fp=fopen(name,"rb");
	if(!fp){
		return(FileNotFoundError);
	}
	status=MMsIOMatrixRead(mx,fp);
	fclose(fp);
	return(status);
}

//----------------------------------------------------------------------------
// 
//------------
static int MMsIOMatrixWrite(MMsMatrix* mx, FILE* fp){
size_t		sz=sizeof(MMsMatrix);
	if(mx->sign==kMatrixRectangularKind_){
		sz+=(mx->nl*mx->nc*sizeof(double));
	}
	else if(mx->sign==kMatrixDiagonalKind_){
		sz+=(((((mx->nl*mx->nl)-mx->nl)/2l)+mx->nl)*sizeof(double));
	}
	else{
		return(-2);
	}
	fseek(fp,0,SEEK_SET);	
	if(fwrite(mx,sz,1,fp)!=1){
		return(-1);
	}
	return(0);
}

//----------------------------------------------------------------------------
// 
//------------
static int MMsIOMatrixWriteFile(MMsMatrix* mx, const char* name){
int		status;
FILE*	fp=fopen(name,"wb");
	if(!fp){
		return(FileNotFoundError);
	}
	status=MMsIOMatrixWrite(mx,fp);
	fclose(fp);
	return(status);
}

//----------------------------------------------------------------------------
// 
//------------
static void add_proc(hca_clss* clss, void* up){
hca_parse*	prs=(hca_parse*)up;
	prs->lclsss=hca_stackadd(prs->lclsss,clss,prs->n);
	prs->n++;
}

//----------------------------------------------------------------------------
// 
//------------
static int	compare(const void* a,const void* b){
hca_clss* ca=(*(hca_clss**)a);
hca_clss* cb=(*(hca_clss**)b);
	return(ca->cid-cb->cid);
}

//----------------------------------------------------------------------------
//
//------------
static void	to_io(hca_clss clss, hca_clss_io* ioclss){
    ioclss->cid=clss.cid;
    ioclss->cidx=clss.cidx;
    
    if(clss.left){
        ioclss->left=clss.left->cid;
    }
    else{
        ioclss->left=0;
    }
    if(clss.right){
        ioclss->right=clss.right->cid;
    }
    else{
        ioclss->right=0;
    }
    if(clss.sup){
        ioclss->sup=clss.sup->cid;
    }
    else{
        ioclss->sup=0;
    }

    ioclss->height=clss.height;
    ioclss->midwidth=clss.midwidth;
    ioclss->ref=clss.ref;
    ioclss->ncmp=clss.ncmp;
}

//----------------------------------------------------------------------------
//
//------------
static void	from_io(hca_clss_io ioclss, hca_clss* clss){
    clss->cid=ioclss.cid;
    clss->cidx=ioclss.cidx;
    clss->left=ioclss.left;
    clss->right=ioclss.right;
    clss->sup=ioclss.sup;
    clss->height=ioclss.height;
    clss->midwidth=ioclss.midwidth;
    clss->ref=ioclss.ref;
    clss->ncmp=ioclss.ncmp;
}

//----------------------------------------------------------------------------
// 
//------------
static int MMsIOClssWrite(	hca_clss** clsss,
							int nclsss,
							FILE* fp){
	
hca_parse	prs={0,NULL};
hca_clss	clss;
hca_clss_io	clssio;
int			i,j;

	if((nclsss==1)&&(clsss[0]->left!=NULL)&&(clsss[0]->right!=NULL)){
		hca_iterate(clsss[0],add_proc,&prs);
	}
	else{
		for(i=0;i<nclsss;i++){
			prs.lclsss=hca_stackadd(prs.lclsss,clsss[i],prs.n);
			prs.n++;
			for(j=0;j<clsss[i]->ncmp;j++){
				prs.lclsss=hca_stackadd(prs.lclsss,clsss[i]->cmp[j],prs.n);
				prs.n++;
			}
		}
	}	
	qsort(prs.lclsss,prs.n,sizeof(hca_clss*),compare);	
	fseek(fp,0,SEEK_SET);
	
	i=1;
	fwrite(&i,sizeof(int),1,fp);
	fwrite(&prs.n,sizeof(int),1,fp);
	for(i=0;i<prs.n;i++){
		
        clss=(*(prs.lclsss[i]));
        
        to_io(clss,&clssio);
        
		fwrite(&clssio,sizeof(hca_clss_io)-sizeof(int),1,fp);
		if(clss.cmp){
			for(j=0;j<clss.ncmp;j++){
				fwrite(&(clss.cmp[j]->cid),sizeof(int),1,fp);
			}
		}
	}
	
	if(prs.lclsss){
		free(prs.lclsss);
	}
	return(0);
}

//----------------------------------------------------------------------------
// 
//------------
static int MMsIOClssRead(	hca_clss** *clsss,
							int *nclsss,
							FILE* fp){
int			i,j,ed;
hca_clss_io	clssio;
hca_clss	clss;
hca_clss*	pclss;
hca_clss**	ppclss;

	*clsss=NULL;
	*nclsss=0;
	
	if(fread(&ed,sizeof(int),1,fp)!=1){
//fprintf(stderr,"fread error\n");
		return(-1);
	}
	if(fread(nclsss,sizeof(int),1,fp)!=1){
//fprintf(stderr,"fread error\n");
		return(-1);
	}

	if(ed!=1){
		swap(sizeof(int),1,nclsss);
	}

//fprintf(stderr,"endian=%d\n",ed);
//fprintf(stderr,"nb=%d\n",*nclsss);

	*clsss=(hca_clss**)malloc((*nclsss)*sizeof(hca_clss*));
	if((*clsss)==NULL){
		*nclsss=0;
		return(-2);
	}
	
	for(i=0;i<(*nclsss);i++){
//fprintf(stderr,"reading class %d\n",i);
		pclss=(hca_clss*)malloc(sizeof(hca_clss));
		(*clsss)[i]=pclss;
        
        if(fread(&clssio,sizeof(hca_clss_io)-sizeof(int),1,fp)!=1){
//fprintf(stderr,"fread error\n");
            (*nclsss)=i+1;
            return(-1);
        }
        
        if(ed!=1){
            swap(sizeof(int),1,&(clssio.cid));
            swap(sizeof(int),1,&(clssio.cidx));
            swap(sizeof(int),1,&(clssio.sup));
            swap(sizeof(int),1,&(clssio.left));
            swap(sizeof(int),1,&(clssio.right));
            swap(sizeof(double),1,&(clssio.height));
            swap(sizeof(double),1,&(clssio.midwidth));
            swap(sizeof(int),1,&(clssio.ref));
            swap(sizeof(int),1,&(clssio.ncmp));
        }
        
		if(clssio.ncmp>0){
			pclss->cmp=(hca_clss**)malloc(clssio.ncmp*sizeof(hca_clss*));
			if(fread(pclss->cmp,clssio.ncmp*sizeof(hca_clss*),1,fp)!=1){
				(*nclsss)=i+1;
				return(-1);
			}
			if(ed!=1){
				swap(sizeof(int),clssio.ncmp,pclss->cmp);
			}
		}
        else{
            pclss->ncmp=0;
            pclss->cmp=NULL;
        }
        from_io(clssio,pclss);
	}
	

//fprintf(stderr,"build graph\n");
	pclss=&clss;
	for(i=0;i<(*nclsss);i++){
		clss.cid=(int)(*clsss)[i]->sup;
		if(clss.cid>0){
			ppclss=bsearch(&pclss,(*clsss),(*nclsss),sizeof(hca_clss*),compare);
			(*clsss)[i]->sup=*ppclss;
			if(!(*clsss)[i]->sup){
				return(-3);
			}
		}
		else{
			(*clsss)[i]->sup=NULL;
		}
		clss.cid=(int)(*clsss)[i]->left;
		if(clss.cid>0){
			ppclss=bsearch(&pclss,(*clsss),(*nclsss),sizeof(hca_clss*),compare);
			(*clsss)[i]->left=*ppclss;
			if(!(*clsss)[i]->left){
				return(-3);
			}
		}
		else{
			(*clsss)[i]->left=NULL;
		}
		clss.cid=(int)(*clsss)[i]->right;
		if(clss.cid>0){
			ppclss=bsearch(&pclss,(*clsss),(*nclsss),sizeof(hca_clss*),compare);
			(*clsss)[i]->right=*ppclss;
			if(!(*clsss)[i]->right){
				return(-3);
			}		
		}
		else{
			(*clsss)[i]->right=NULL;
		}
		if((*clsss)[i]->ncmp>0){
			for(j=0;j<(*clsss)[i]->ncmp;j++){
				clss.cid=(int)(*clsss)[i]->cmp[j];
				ppclss=bsearch(&pclss,(*clsss),(*nclsss),sizeof(hca_clss*),compare);
				(*clsss)[i]->cmp[j]=*ppclss;
				if(!(*clsss)[i]->cmp[j]){
					return(-3);
				}
			}
		}
	}

//fprintf(stderr,"quit\n");

	return(0);
}

//----------------------------------------------------------------------------
// 
//------------
static int MMsIOClssWriteFile(	hca_clss** clsss,
								int nclsss,
								const char* name){
int		status;
FILE*	fp=fopen(name,"wb");
	if(!fp){
		return(FileNotFoundError);
	}
	status=MMsIOClssWrite(clsss,nclsss,fp);
	fclose(fp);
	return(status);
}

//----------------------------------------------------------------------------
// 
//------------
static int MMsIOClssReadFile(	hca_clss** *clsss,
								int *nclsss,
								const char* name){
int		status;
FILE*	fp=fopen(name,"rb");
	if(!fp){
		return(FileNotFoundError);
	}
	status=MMsIOClssRead(clsss,nclsss,fp);
	fclose(fp);
	return(status);
}

//----------------------------------------------------------------------------
// 
//------------
int	MMsIOWriteAnalysis	(	const char* path,
							const char* name,
							mmx_analysis* ana){
int		status=0;
char	hname[FILENAME_MAX];
char	cwd[PATH_MAX];
	
char*	d=strrchr(name,'.');
	
	switch(ana->kind){
		case _kAnalysisUndef:
			status=-1;
			break;
		case _kAnalysisUNI:
			if(d==NULL){
				status=make_package(path,name,"uni");
				sprintf(hname,"%s.uni",name);
			}
			else{
				status=make_package(path,name,"");
				sprintf(hname,"%s",name);
			}
			break;
		case _kAnalysisBIV:
			if(d==NULL){
				status=make_package(path,name,"biv");
				sprintf(hname,"%s.biv",name);
			}
			else{
				status=make_package(path,name,"");
				sprintf(hname,"%s",name);
			}
			break;
		case _kAnalysisMLT:
			if(d==NULL){
				status=make_package(path,name,"mlt");
				sprintf(hname,"%s.mlt",name);
			}
			else{
				status=make_package(path,name,"");
				sprintf(hname,"%s",name);
			}
			break;
		case _kAnalysisPC:
		case _kAnalysisPCn:
			if(d==NULL){
				status=make_package(path,name,"pca");
				sprintf(hname,"%s.pca",name);
			}
			else{
				status=make_package(path,name,"");
				sprintf(hname,"%s",name);
			}
			break;
		case _kAnalysisFc:
			if(d==NULL){
				status=make_package(path,name,"fca");
				sprintf(hname,"%s.fca",name);
			}
			else{
				status=make_package(path,name,"");
				sprintf(hname,"%s",name);
			}
			break;
		case _kAnalysisHC:
			if(d==NULL){
				status=make_package(path,name,"hca");
				sprintf(hname,"%s.hca",name);
			}
			else{
				status=make_package(path,name,"");
				sprintf(hname,"%s",name);
			}
			break;
		case _kAnalysisIRA:
			if(d==NULL){
				status=make_package(path,name,"ira");
				sprintf(hname,"%s.ira",name);
			}
			else{
				status=make_package(path,name,"");
				sprintf(hname,"%s",name);
			}
			break;
		default:
			status=-1;
			break;
	}
	if(status){
		return(status);
	}
	
	if(cd(path,cwd)){
		return(-1);
	}
	if(cd(hname,NULL)){
		cd(cwd,NULL);
		return(-1);
	}
	if(cd("Contents/Datas",NULL)){
		cd(cwd,NULL);
		return(-1);
	}

	status=MMsIOKindWrite(ana->kind);
	if(status){
		return(status);
	}
		
	if(ana->data){
		status=MMsIOMatrixWriteFile(ana->data,"data.mmx");
		if(status){
			cd(cwd,NULL);
			return(status);
		}
	}
	if(ana->indic){
		status=MMsIOMatrixWriteFile(ana->indic,"indic.mmx");
		if(status){
			cd(cwd,NULL);
			return(status);
		}
	}
	switch(ana->kind){
		case _kAnalysisUndef:
			break;
		case _kAnalysisUNI:
			if(((mmx_uni*)ana)->clss){
				status=MMsIOMatrixWriteFile(((mmx_uni*)ana)->clss,"clss.mmx");
			}
			break;
		case _kAnalysisBIV:
		case _kAnalysisMLT:
			break;
		case _kAnalysisPC:
		case _kAnalysisPCn:
			if(((mmx_pca*)ana)->cv){
				status=MMsIOMatrixWriteFile(((mmx_pca*)ana)->cv,"cv.mmx");
				if(status){
					break;
				}
			}
			if(((mmx_pca*)ana)->chr){
				status=MMsIOMatrixWriteFile(((mmx_pca*)ana)->chr,"chr.mmx");
				if(status){
					break;
				}
			}
			if(((mmx_pca*)ana)->var){
				status=MMsIOMatrixWriteFile(((mmx_pca*)ana)->var,"var.mmx");
			}
			break;
		case _kAnalysisFc:
			if(((mmx_fa*)ana)->ctr){
				status=MMsIOMatrixWriteFile(((mmx_fa*)ana)->ctr,"ctr.mmx");
				if(status){
					break;
				}
			}
			if(((mmx_fa*)ana)->chr){
				status=MMsIOMatrixWriteFile(((mmx_fa*)ana)->chr,"chr.mmx");
				if(status){
					break;
				}
			}
			if(((mmx_fa*)ana)->var){
				status=MMsIOMatrixWriteFile(((mmx_fa*)ana)->var,"var.mmx");
			}
			break;
		case _kAnalysisHC:
			if(((mmx_hc*)ana)->dst){
				status=MMsIOMatrixWriteFile(((mmx_hc*)ana)->dst,"dst.mmx");
			}
			break;
		case _kAnalysisIRA:
			if(((mmx_ira*)ana)->center){
				status=MMsIOMatrixWriteFile(((mmx_ira*)ana)->center,"center.mmx");
			}		
			break;
		default:
			break;
	}
	cd(cwd,NULL);
	return(status);
}

//----------------------------------------------------------------------------
// 
//------------
int	MMsIOReadAnalysis	(	const char* path,
							const char* name,
							mmx_analysis* ana){
int		status=0;
char	cwd[PATH_MAX];

	switch(ana->kind){
		case _kAnalysisUndef:
fprintf(stderr,"_kAnalysisUndef\n");
			return(-1);
			break;
		case _kAnalysisUNI:
		case _kAnalysisBIV:
		case _kAnalysisMLT:
		case _kAnalysisPC:
		case _kAnalysisPCn:
		case _kAnalysisFc:
		case _kAnalysisHC:
		case _kAnalysisIRA:
			break;
		default:
fprintf(stderr,"bad sign (%d)\n",ana->kind);
			return(-1);
			break;
	}

	if(cd(path,cwd)){
fprintf(stderr,"cd %s failed\n",path);
		return(-1);
	}
	if(cd(name,NULL)){
fprintf(stderr,"cd %s failed\n",name);
		cd(cwd,NULL);
		return(-1);
	}
	if(cd("Contents/Datas",NULL)){
fprintf(stderr,"cd %s failed\n","Contents/Datas");
		cd(cwd,NULL);
		return(-1);
	}

	status=MMsIOKindRead(&ana->kind);
	if(status){
fprintf(stderr,"MMsIOKindRead\n");
		return(status);
	}
fprintf(stderr,"analysis kind=%d\n",ana->kind);

fprintf(stderr,"data\n");
	status=MMsIOMatrixReadFile(&ana->data,"data.mmx");
	if(status){
fprintf(stderr,"MMsIOMatrixReadFile returns %d\n",status);
		cd(cwd,NULL);
		return(status);
	}
	
fprintf(stderr,"indic\n");
	status=MMsIOMatrixReadFile(&ana->indic,"indic.mmx");
	if(status){
		cd(cwd,NULL);
		return(status);
	}
	
    switch(ana->kind){
        case _kAnalysisUNI:
fprintf(stderr,"clss\n");
           status=MMsIOMatrixReadFile(&((mmx_uni*)ana)->clss,"clss.mmx");
            if(status==FileNotFoundError){
                status=0;
            }
            break;
        case _kAnalysisPC:
        case _kAnalysisPCn:
fprintf(stderr,"cv\n");
            status=MMsIOMatrixReadFile(&((mmx_pca*)ana)->cv,"cv.mmx");
            if(status){
                break;
            }
fprintf(stderr,"chr\n");
            status=MMsIOMatrixReadFile(&((mmx_pca*)ana)->chr,"chr.mmx");
            if(status){
                break;
            }
fprintf(stderr,"var\n");
            status=MMsIOMatrixReadFile(&((mmx_pca*)ana)->var,"var.mmx");
            break;
        case _kAnalysisFc:
fprintf(stderr,"ctr\n");
            status=MMsIOMatrixReadFile(&((mmx_fa*)ana)->ctr,"ctr.mmx");
            if(status){
                break;
            }
fprintf(stderr,"chr\n");
            status=MMsIOMatrixReadFile(&((mmx_fa*)ana)->chr,"chr.mmx");
            if(status){
                break;
            }
fprintf(stderr,"var\n");
            status=MMsIOMatrixReadFile(&((mmx_fa*)ana)->var,"var.mmx");
            break;
        case _kAnalysisHC:
fprintf(stderr,"dst\n");
            status=MMsIOMatrixReadFile(&((mmx_hc*)ana)->dst,"dst.mmx");
            break;
        case _kAnalysisIRA:
fprintf(stderr,"center\n");
            status=MMsIOMatrixReadFile(&((mmx_ira*)ana)->center,"center.mmx");
            break;
    }
    cd(cwd,NULL);
    return(status);
}

//----------------------------------------------------------------------------
// 
//------------
int	MMsIOWriteClassification	(	const char* path,
									const char* name,
									hca_clss** clsss,
									int nclss,
									int clssk){
char	cwd[PATH_MAX];
char	hname[FILENAME_MAX];

	switch(clssk){
		case kHCAMethodNone:
			return(-1);
			break;
		case kHCADiameterMethod:
			sprintf(hname,"diameter.cls");
			break;
		case kHCAMinimalMethod:
			sprintf(hname,"minimum.cls");
			break;
		case kHCAAverageMethod:
			sprintf(hname,"average.cls");
			break;
		case kIRAStdMethod:
			sprintf(hname,"relocation.cls");
			break;
		default:
			return(-1);
			break;
	}
	
	if(cd(path,cwd)){
		return(-1);
	}
	if(cd(name,NULL)){
		cd(cwd,NULL);
		return(-1);
	}
	if(cd("Contents/Datas",NULL)){
		cd(cwd,NULL);
		return(-1);
	}
int	status=MMsIOClssWriteFile(clsss,nclss,hname);
	cd(cwd,NULL);
	return(status);
}

//----------------------------------------------------------------------------
// 
//------------
int	MMsIOReadClassification	(	const char* path,
								const char* name,
								hca_clss*** clsss,
								int* nclss,
								int clssk){
char		cwd[PATH_MAX];
char		hname[FILENAME_MAX];
hca_clss**	lclsss;
int			i,lnclss;

	(*nclss)=0;
	(*clsss)=NULL;

	switch(clssk){
		case kHCAMethodNone:
			return(-1);
			break;
		case kHCADiameterMethod:
			sprintf(hname,"diameter.cls");
			break;
		case kHCAMinimalMethod:
			sprintf(hname,"minimum.cls");
			break;
		case kHCAAverageMethod:
			sprintf(hname,"average.cls");
			break;
		case kIRAStdMethod:
			sprintf(hname,"relocation.cls");
			break;
		default:
			return(-1);
			break;
	}
	
	if(cd(path,cwd)){
		return(-1);
	}
	if(cd(name,NULL)){
		cd(cwd,NULL);
		return(-1);
	}
	if(cd("Contents/Datas",NULL)){
		cd(cwd,NULL);
		return(-1);
	}
int	status=MMsIOClssReadFile(&lclsss,&lnclss,hname);
	cd(cwd,NULL);
	
	if(clssk==kIRAStdMethod){
		for(i=0;i<lnclss;i++){
			if(lclsss[i]->ncmp>0){
				(*clsss)=hca_stackadd((*clsss),lclsss[i],(*nclss));
				(*nclss)++;
			}
		}
	}
	else{
		for(i=lnclss-1;i>=0;i--){
			if(lclsss[i]->sup==NULL){
				(*clsss)=hca_stackadd((*clsss),lclsss[i],(*nclss));
				(*nclss)++;
				break;
			}
		}
	}
	return(status);
}

//----------------------------------------------------------------------------
// 
//------------
int	MMsIOWriteParam	(	const char* path,
						const char* name,
						void* prm,
						int sz,
						const char* pname){	
char	cwd[PATH_MAX];

	if(cd(path,cwd)){
		return(-1);
	}
	if(cd(name,NULL)){
		cd(cwd,NULL);
		return(-1);
	}
	if(cd("Contents/Resources",NULL)){
		cd(cwd,NULL);
		return(-1);
	}
	
FILE* fp=fopen(pname,"w");
	if(!fp){
		cd(cwd,NULL);
		return(FileNotFoundError);
	}
int status=0;
	if(fwrite(prm,sz,1,fp)!=1){
		status=-1;
	}
	fclose(fp);

	cd(cwd,NULL);
	return(status);
}

//----------------------------------------------------------------------------
// 
//------------
int	MMsIOReadParam	(	const char* path,
						const char* name,
						void** prm,
						int* sz,
						const char* pname){
char	cwd[PATH_MAX];

	if(cd(path,cwd)){
		return(-1);
	}
	if(cd(name,NULL)){
		cd(cwd,NULL);
		return(-1);
	}
	if(cd("Contents/Resources",NULL)){
		cd(cwd,NULL);
		return(-1);
	}	
	
FILE* fp=fopen(pname,"r");
	if(!fp){
		cd(cwd,NULL);
		return(FileNotFoundError);
	}
int status=0;

	fseek(fp,0,SEEK_END);
	(*sz)=ftell(fp);
	fseek(fp,0,SEEK_SET);
	(*prm)=malloc(*sz);
	if(!(*prm)){
		fclose(fp);
		cd(cwd,NULL);
		return(-2);
	}
	if(fread((*prm),(*sz),1,fp)!=1){
		status=-1;
	}
	fclose(fp);

	cd(cwd,NULL);

	return(status);
}
