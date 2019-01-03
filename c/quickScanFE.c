#include <stdio.h>
#include <stdlib.h>
#include <string.h>


float densi[35278226];
float FE[35278226];
//float r[2];

float min(float * p, int n);

float getAve(float n , int * p);

void main (int argc, char **argv ){

       char *filename=argv[1];
       int slocal=atoi(argv[2])/2;
       int llocal=atoi(argv[3])/2;
       float thresh_hold=atof(argv[4]);
       //printf("%s %d %d\n",filename,slocal,llocal);

       //readin densi
       FILE *fp;
       fp=fopen(filename,"r");
       if(fp == NULL) perror("err when opening file");
       float buf=0;
       int cnt=0;
       while(!feof(fp)){
         //fgets(bufline,100,fp); 
         //float value = atof(bufline);  
         int ferr=fscanf(fp,"%f",&buf);
         //printf("%f\n",buf);
         //if(ferr != 1) perror("not one element\n"); ;
         densi[cnt++] = buf;  
         if(cnt % 1000000 == 0) fprintf(stderr,"#");
       }  
       fprintf(stderr,"\ntotal %d \n",cnt);
       fclose(fp);
       //printf ("%d %d\n",cnt,sizeof(densi)/sizeof(float));

       int i;
       for(i=1;i < cnt; i++){
         //printf("%f\n",densi[i]);  
         if(i % 1000000 == 0) fprintf(stderr,"#");
         //slocal window
         int slocal_left=i-slocal; 
         int slocal_right=i+slocal; 
         if(slocal_left <= 0){ slocal_left = 1; slocal_right=i+slocal*2;}
         if(slocal_right > cnt ){ slocal_left = cnt-slocal*2; slocal_right=cnt-1;}
         int j;
         float cnt_inner=0;
         float sum_inner=0;
         for(j=slocal_left; j<= slocal_right; j++){
            sum_inner += densi[j];      
            cnt_inner++;
         }
         float ave_slocal = sum_inner /  cnt_inner;         
         //printf("%f\n",ave_slocal);    

         //llocal window
         int llocal_left=i-llocal;
         int llocal_right=i+llocal;
         if(llocal_left <= 0){ llocal_left = 1; llocal_right=i+llocal*2;}
         if(llocal_right > cnt ){ llocal_left = cnt-llocal*2; llocal_right=cnt-1;}
         cnt_inner=0;
         sum_inner=0;
         for(j=llocal_left; j<= llocal_right; j++){
            sum_inner += densi[j];
            cnt_inner++;
         }
         float ave_llocal = sum_inner / cnt_inner; 
         //printf("%f\n",ave_llocal);    

         //pick the smaller one
         //float data[2]={ave_slocal,ave_llocal,};
         //float m=min(data,2);
         float m=0;
         if(ave_slocal < ave_llocal){m = ave_slocal;}else{m = ave_llocal;}
         
         FE[i]=m;
         if(m/densi[i] >= thresh_hold){printf("H1_japo\t%d\t%d\t%d\n",i,i+1,1);}else{printf("H1_japo\t%d\t%d\t%d\n",i,i+1,0);}
         //printf("H1_japo\t%d\t%d\t%f\n",i,i+1,m/densi[i]);      
       }//outter for

       //int test[3] = {10,1,2};
       //int x;
       //for(x=0;x<3;x++) printf("%d\n",test[x]);
       //int m = min(test,3);
       //printf ("min is %d\n",m);
       //float result = getAve(3,test);
       //printf("ave is %.6f\n",result);

       fprintf(stderr,"\n")     ; 

}//end of main


//###sub###

float min(float * p, int n){
   int i;
   int min= *p;// printf("within func paras is %d  %d\n",m, n);
   for(i=0; i<n; i++){
     //printf("current value is %d\n",*p);
     if(min > *p) { 
         min = *p ;
         //printf("yes:%d\n",*p); 
      }
     p++;
   }
   return ( min );
}

float getAve(float n ,int * p){
   int i;
   float ave=0;
   float sum=0;
   for(i=0;i<n;i++){
     sum+=*p;
     p++;
   }
   ave = sum/n ;
   return ( ave );
}


