#include<stdio.h>
#include<stdlib.h>
#include<time.h>
#include<math.h>

//exponential
double **e(double **ptr,int satir,int sutun,float sigma){
		int i,j,x=-(satir-1)/2,y=-(sutun-1)/2;
		double toplam=0;
		for(i=0;i<satir;i++)
			for(j=0;j<sutun;j++){
					ptr[i][j]=exp(-((x+i)*(x+i)+(y+j)*(y+j))/(2*sigma*sigma));
					toplam+=ptr[i][j];
			}
		for(i=0;i<satir;i++)
			for(j=0;j<sutun;j++)
					ptr[i][j]/=toplam;

		return ptr;
	}
void free_array(double **ptr, int satir){
     int i;
     for(i=0;i<satir;i++)
        free(ptr[i]);
     free(ptr);
	}
//memory allocation
double **yer(double **ptr,int satir,int sutun){
		int i;
		ptr=(double **)calloc(satir,sizeof(double));
		for(i=0;i<satir;i++)
			ptr[i]=(double *)calloc(sutun,sizeof(double));
		return ptr;
	}
//Sxx Syy Sxy calculation
double **S(double **ptr,double **ptr1,int satir,int sutun,int satir1,int sutun1){
			int i,j,k,l;
			double **toplam;
			toplam=yer(toplam,satir,sutun);
			for(i=0;i<satir;i++)
				for(j=0;j<sutun;j++)
					for(k=0;k<satir1;k++)
							for(l=0;l<sutun1;l++)
									toplam[i][j]+=ptr[i+k][j+l]*ptr1[k][l];
			return toplam;
	}
int main(){
		double **dizi2,**filtre,**Ix,**Iy;
		double **Ixx,**Iyy,**Ixy,**R,trash;
		double **Sxx=0,**Syy,**Sxy;
		float sigma;
		int satir1,satir2,x,sutun1,sutun2,y,i,j,k,l,satir3,sutun3;
		FILE *ptr;
		int **dizi1,maxRenk;
		char tur[3];

        if((ptr=fopen("D:\\lena.pgm","rb"))==NULL){
				printf("cannot open file");
		        return 0;
        		}

		fgets(tur,sizeof(tur),ptr);
		if((tur[0]!='P' && tur[1]!='5') || (tur[0]!='P' && tur[1]!='2')){
			printf("wrong file type");
            }
		fscanf(ptr,"%d",&satir1);
		fscanf(ptr,"%d",&sutun1);
		fscanf(ptr,"%d",&maxRenk);

		dizi1=(int **)calloc(satir1,sizeof(int));
		if(dizi1)
		for(i=0;i<satir1;i++)
			dizi1[i]=(int *)calloc(sutun1,sizeof(int));
			if(dizi1[i])
		for(i=0;i<satir1;i++)
			for(j=0;j<sutun1;j++)
				fscanf(ptr,"%d",&dizi1[i][j]);

		printf("enter row and column for exp\n");
		scanf("%d%d",&satir2,&sutun2);
		printf("enter sigma");
		scanf("%f",&sigma);
		x=(satir1-satir2)+1;
		y=(sutun1-sutun2)+1;

		dizi2=yer(dizi2,satir2,sutun2);
		if(dizi2){
			dizi2=e(dizi2,satir2,sutun2,sigma);

		}
		filtre=yer(filtre,x,y);
		if(filtre){
			for(i=0;i<x;i++)
				for(j=0;j<y;j++)
					for(k=0;k<satir2;k++)
							for(l=0;l<sutun2;l++)
									filtre[i][j]+=dizi1[i+k][j+l]*dizi2[k][l];

			}
			fclose(ptr);
//for Ix Ixx
		Ix=yer(Ix,x,y-2);
		Ixx=yer(Ixx,x,y-2);
		if(Ix){
		  for(i=0;i<x;i++)
			 for(j=0;j<y-2;j++){
				 Ix[i][j]=filtre[i][j+2]-2*filtre[i][j+1]+filtre[i][j];
				 	 if(Ixx)
						 Ixx[i][j]=pow(Ix[i][j],2);
			}
		}
//for Iy Iyy
		Iy=yer(Iy,x-2,y);
		Iyy=yer(Iyy,x-2,y);
		if(Iy){
		for(i=0;i<x-2;i++)
			for(j=0;j<y;j++){
				Iy[i][j]=filtre[i+2][j]-2*filtre[i+1][j]+filtre[i][j];
					if(Iyy)
						Iyy[i][j]=pow(Iy[i][j],2);
			}
		}
//for Ixy
		Ixy=yer(Ixy,x-2,y-2);
		if(Ixy){
			for(i=0;i<x-2;i++)
				for(j=0;j<y-2;j++)
					Ixy[i][j]=Ix[i][j]*(Iy[i][j]);
		}
//satir3 sutun3 S borders
		satir3=(x-satir2)-1;
		sutun3=(y-sutun2)-1;
//send Sxx to function
		Sxx=yer(Sxx,x-satir2+1,sutun3);
			if(Sxx)
				Sxx=S(Ixx,dizi2,x-satir2+1,sutun3,satir2,sutun2);

//send Syy to function
		Syy=yer(Syy,satir3,y-sutun2+1);
			if(Syy)
				Syy=S(Iyy,dizi2,satir3,y-sutun2+1,satir2,sutun2);

//send Sxy to function
		Sxy=yer(Sxy,satir3,sutun3);
			if(Sxy)
				Sxy=S(Ixy,dizi2,satir3,sutun3,satir2,sutun2);

		printf("enter treshold:");
		scanf("%lf",&trash);
		printf("%d %d",satir3,sutun3);
//R 
		
		R=yer(R,satir3,sutun3);
		if(R){
			for(i=0;i<satir3;i++)
				for(j=0;j<sutun3;j++){
					R[i][j]=(Sxx[i][j]*Syy[i][j]-pow(Sxy[i][j],2))-0.04*pow((Sxx[i][j]+Syy[i][j]),2);
			}
		}
		printf("%d %d",satir3,sutun3);
		for(i=0;i<satir3;i++)
			for(j=0;j<sutun3;j++)
		      if(R[i][j]>trash)
            dizi1[i+10][j+10]=255;
    printf("%d %d",satir3,sutun3);
		if((ptr=fopen("D:\\sonHalixx1.pgm","wb"))==NULL)
			printf("cannot open file.");
		else{
			fprintf(ptr,"%s ",tur);
			fprintf(ptr,"%d %d ",sutun1,satir1);
			fprintf(ptr,"%d ",maxRenk);
			for(i=0;i<satir1;i++)
				for(j=0;j<sutun1;j++)
					fprintf(ptr,"%c",(char) dizi1[i][j]);
        }
    fclose(ptr);
    
    for(i=0;i<satir1;i++)
        free(dizi1[i]);
    free(dizi1);
    free_array(dizi2,satir2);
    free_array(filtre,x);
    free_array(Ix,x);
    free_array(Ixx,x);
    free_array(Iy,x-2);
    free_array(Iyy,x-2);
    free_array(Ixy,x-2);
    free_array(Sxx,x-satir2+1);
    free_array(Syy,satir3);
    free_array(Sxy,satir3);
    free_array(R,satir3);
system("pause");
return 0;
}
