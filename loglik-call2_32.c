#include <R.h>
#include <Rmath.h>
#include <Rinternals.h>
#include <math.h>
#include <stdio.h>


#define DIMS(x)							\
  INTEGER(coerceVector(getAttrib((x), R_DimSymbol), INTSXP))   \

#define e 2.718281828

/*long double factorial(long double initial,int n){   //compute initial*n!
    long double result=initial;
    for(int i=1;i<n+1;i++){
        result *= n;
            }
//    if(result <0.0){
//        printf("%s\n","overflow");
//    }

    return result;
}
*/
double logfactorial(int n){   //compute log(n!)
    double result=0.0;
    for(int i=1;i<n+1;i++){
        result += log(i);
    }
    return result;
}


void count(double *prob, int nmut,int* nunsign, int* nc){//, double *loglik,int j){
    for(int i=0;i<nmut;i++){
        if((prob[i] > 1e-6) & (prob[i]< 1-(1e-6))){
            (*nunsign)++;
        }else if(prob[i]> 1-(1e-6)){
            (*nc)++;
        }
    }
}


void causind(int *rcausal,double *prob, int nmut ,int N) /*generate causal matrix for a sample with nmut mutaitons.nunde_mut:*/
{
    int num,j;
    num=N;
    for(j=0;j<nmut;j++){
        if(prob[j]<1e-6) {rcausal[j]=0;}
        else if(prob[j]>1-(1e-6) ) {rcausal[j]=1;}
        else{
        if(num>0){
            rcausal[j]=num % 2;
            num=num-rcausal[j];
            num /= 2;
        }else{
            rcausal[j]=0;
        }
     }
   }
}


double conprob( int *rcausal, double *rprob_j, int nmut) /*caculate p(n_j|m_j),rcausal is the indicate of causal vec for a situation  （pnm is the output，rprob_j is the input probability）*/
{
    int j;
    double pnm;
    pnm=1.0;
    for(j = 0; j< nmut;j++){
      if(rcausal[j]==0){
      pnm *= (1-rprob_j[j]);
      }else{
        pnm *= rprob_j[j];
      }
    }
    return(pnm);
}

void caustable(int *rcausal, int *rY, int *newrY_j, int i, int nY0, int nY1)/*casual mut table of sample j corresponding to causal matrix(ncausal=2^n) */
/* rcausal: causal vector; rY_j: mut vector for jth sample; nY_j: length of rY_j ; newrY_j: output; 
 newY_j: first row is for damage mut and second row is for non_damage mut*/
{
  int index=0;
    for(int j=0; j < nY1; j++){
        newrY_j[j]=0;
        newrY_j[j+nY1]=0;
      for(int k=0;k<rY[i+j*nY0];k++){
          newrY_j[j] += rcausal[index];
          newrY_j[j+nY1] += 1-rcausal[index];
          index += 1;
      }  
    }
}


double sum(double *P, int *a, int *nP, int *na)
{
    int i, j;
    double tmp, sum;
    sum = 0.0;
    for (i = 0; i < na[0]; i++) {
        tmp = 1.0;
        for (j = 0; j < na[1]; j++)
            tmp *= P[a[i + j*na[0]] - 1 + j*nP[0]];
        sum += tmp;
    }
    return sum;
}

double sumnon(double *Pnon, int *anon, int *nPnon, int *nanon)
{
    int i, j;
    double tmp, sum;
    sum = 0.0;
    for (i = 0; i < nanon[0]; i++) {
        tmp = 1.0;
        for (j = 0; j < nanon[1]; j++)
        tmp *= Pnon[anon[i + j*nanon[0]] - 1 + nPnon[0]];
        sum += tmp;
    }
    return sum;
}


double partsum(double *P, int *a, int *nP, int *na, int beginrow, int endrow)  //caculate part some of multiplication of probability values correspoing to A matrix from beginrow to endrow.
{
    int i, j;
    double tmp, sum;
    sum = 0.0;
    for (i = beginrow; i < endrow; i++) {
        tmp = 1.0;
        for (j = 0; j < na[1]; j++)
            tmp *= P[a[i + j*na[0]] - 1 + j*nP[0]];
        sum += tmp;
    }
    return sum;
}



  
double logsum(double *P, int *a, int *nP, int *na)
{
  int i, j;
  double tmp, sum;
  sum = 0.0;
  for (i = 0; i < na[0]; i++) {
    tmp = 1.0;
    for (j = 0; j < na[1]; j++) 
      tmp *= P[a[i + j*na[0]] - 1 + j*nP[0]];
    sum += tmp;
  }
  return log(sum);
}



double dvecsum(double *x, int n)
{
  int i;
  double sum = 0.0;
  for (i = 0; i < n; i++)
    sum += x[i];
  return sum;    
}

int  partrowsum(int *x, int nx,int start){     //sum of a vector x from a start place , nx is the length of the part of x we want to calculate.
    int i,sum=0;
    for (i = start; i < start+nx; i++) {
        sum += x[i];
    }
    return sum;
}

int irowsum(int *x, int* nx, int row)
{
  int i, sum = 0;
  for (i = 0; i < nx[1]; i++) {
    sum += x[row + i*nx[0]];
  }
  return sum;    
}

void dsubmat(double *old, double *new, int *nold, int *nnew,
	     int* rowold, int* colold)
{
  int i, j;
  if (rowold != NULL && colold == NULL) {
    for (i = 0; i < nnew[0]; i++)
      for (j = 0; j < nnew[1]; j++)
	new[i + j*nnew[0]] = old[rowold[i] + j*nold[0]];
  } else if (rowold == NULL && colold != NULL) {
    for (i = 0; i < nnew[0]; i++)
      for (j = 0; j < nnew[1]; j++)
	new[i + j*nnew[0]] = old[i + colold[j]*nold[0]];
  } else if (rowold != NULL && colold != NULL) {
    for (i = 0; i < nnew[0]; i++)
      for (j = 0; j < nnew[1]; j++)
	new[i + j*nnew[0]] = old[rowold[i] + colold[j]*nold[0]];
  } else {
    for (i = 0; i < nnew[0]; i++)
      for (j = 0; j < nnew[1]; j++)
	new[i + j*nnew[0]] = old[i + j*nold[0]];
  }
}

void isubmat(int *old, int *new, int *nold, int *nnew,
	     int* rowold, int* colold)
{
  int i, j;
  if (rowold != NULL && colold == NULL) {
    for (i = 0; i < nnew[0]; i++)
      for (j = 0; j < nnew[1]; j++)
	new[i + j*nnew[0]] = old[rowold[i] + j*nold[0]];
  } else if (rowold == NULL && colold != NULL) {
    for (i = 0; i < nnew[0]; i++)
      for (j = 0; j < nnew[1]; j++)
	new[i + j*nnew[0]] = old[i + colold[j]*nold[0]];
  } else if (rowold != NULL && colold != NULL) {
    for (i = 0; i < nnew[0]; i++)
      for (j = 0; j < nnew[1]; j++)
	new[i + j*nnew[0]] = old[rowold[i] + colold[j]*nold[0]];
  } else {
    for (i = 0; i < nnew[0]; i++)
      for (j = 0; j < nnew[1]; j++)
	new[i + j*nnew[0]] = old[i + j*nold[0]];
  }
}


 SEXP loglik(SEXP x, SEXP AllPtrue, SEXP divsel,
	    SEXP Y, SEXP a, SEXP MAX,SEXP prob,SEXP prob_nondammut)
    
{
  int nx, ndivsel, *nAll, *nY, *nP, *na_N;
  int *rY, *rdivsel, *rMAX, *ra_N;
  int i, j, k, l, N,N1,c;

  double *rx, *rAllPtrue, *rans,*rprob_nondammut;
  double *AllP, *loglik;
  double *prob_i;
  double pweight_c;
  //  int *newY_i;
    
  
 
  SEXP ans;

  PROTECT(ans = allocVector(REALSXP, 1));
  rans = REAL(ans);
  nx	 = length(x);
  ndivsel = length(divsel);
  nAll = DIMS(AllPtrue);
  nY	 = DIMS(Y);
  nP = (int *)R_alloc(2, sizeof(int));

  rx = REAL(x);
  rAllPtrue = REAL(AllPtrue);
  rdivsel = INTEGER(divsel); 
  rY = INTEGER(Y);
  rMAX = INTEGER(MAX);
    rprob_nondammut = REAL(prob_nondammut);

   loglik = (double *)R_alloc(nY[0], sizeof(double));

  AllP = (double *)R_alloc(nAll[0]*nAll[1], sizeof(double));
  for (i = 0; i < nAll[0]; i++)
    for (j = 0; j < nAll[1]; j++)
      AllP[i + j*nAll[0]] = rAllPtrue[i + j*nAll[0]];
    
    double sumP=0.0;
  for(j=0; j<ndivsel ;j++)
     {
         for (i = 0; i < nAll[0] - 1; i++){
              AllP[i + (rdivsel[j] - 1)*nAll[0]] = exp(rx[i]);
              sumP += AllP[i + (rdivsel[j] - 1)*nAll[0]] ;
         }
     }
    
    sumP /= ndivsel;

    
    for(j=0; j<ndivsel ;j++)
    {
        for (i = 0; i < nAll[0] - 1; i++){
            AllP[i + (rdivsel[j] - 1)*nAll[0]] /= sumP+1;
            
        }
            AllP[nAll[0]-1 + (rdivsel[j]-1)*nAll[0]] = 1.0/(sumP+1);
    }

    
 //   for (i = 0; i < nAll[0] ; i++){
 //       printf("%f\t",AllP[i + (rdivsel[0] - 1)*nAll[0]]);
 //       }
  //  printf("\n");

 /*
       for(j=0; j<nAll[0] ;j++)
    {
        for (int k = 0; k < nAll[1]; k++){
           printf("%f\t",AllP[j + k*nAll[0]]);
        }
        printf("\n");
    }  */


 /*  for(j=0; j<ndivsel ;j++)
    {
      for (i = 0; i < nAll[0] - 1; i++)
	AllP[i + (rdivsel[j] - 1)*nAll[0]] = rx[i];
      AllP[i + (rdivsel[j]-1)*nAll[0]] = 1 - dvecsum(rx, nx);
    }
 */
    

  for (i = 0; i < nY[0]; i++) {
      loglik[i]=0.0;
      N = irowsum(rY, nY, i);
      prob_i=REAL(VECTOR_ELT(prob,i));
      int nc=0;
      int nunsign=0;
      count(prob_i,N,&nunsign,&nc);
      
      for(c=0;c<pow(2,nunsign);c++){
          int rcausal[N];
          causind(&rcausal[0],prob_i,N,c);
          pweight_c=conprob(&rcausal[0],prob_i,N);
          
          int newY_i[nY[1]*2];
          caustable(rcausal, rY , newY_i, i, nY[0],nY[1] );
          
//for(int d = 0;d<2*nY[1];d++){
//printf("%d\t",newY_i[d]);
//}
//printf("\n");

          N1=partrowsum(&newY_i[0],nY[1],0); //# of damage mut
          int N2=partrowsum(&newY_i[0],nY[1],nY[1]); //# of non dam mut
//printf("%d\t",N1);
//printf("%d\n",N2);          
          double logprob_nonmut_sum=0.0;
          if (N2!=0){
       /*       prob_nonmut_sum=1.0;
          }else{
              int rownon[N-N1];
              l=0;
              for (j = 0; j < nY[1]; j++)
              for (k = 0; k < newY_i[nY[1]+j]; k++) {
                  rownon[l++] = j;   //indicate for non_dammut, as rowAll for dammut
              }
              */
              for(j=0;j<nY[1];j++){
                  for(k = 0; k < newY_i[nY[1]+j]; k++){
                      if(rprob_nondammut[j]!=0.0){
                          logprob_nonmut_sum += log(rprob_nondammut[j]);
//printf("%d\t",j);
//printf("%f\n",logprob_nonmut_sum);

                      }else{

//printf("jump\t");
//printf("%d\t",k);
//printf("%f\n",rprob_nondammut[k]);
                          goto end;
                      }
                  }
              }
              logprob_nonmut_sum = logprob_nonmut_sum+logfactorial(N2);
              
              
 /*             nPnon[0] = N2;
              nPnon[1] = 1;
              double Pnon[nP[0]*nP[1]];
              dsubmat(prob_nondammut, &Pnon[0], nAll, nPnon, &rownon[0], NULL);
              prob_nonmut_sum=sumnon(P, ra_N, nP, na_N) */
          }
          
          if(N1==0){
              loglik[i] += pweight_c*pow(e,logprob_nonmut_sum);
          }else{
              nP[0] = N1;
              nP[1] = N1;
              double P[nP[0]*nP[1]];
              int rowAll[nP[0]];
              
          
              l = 0;
              for (j = 0; j < nY[1]; j++)
                  for (k = 0; k < newY_i[j]; k++) {
                      rowAll[l++] = j;
                  }
              
        
          
              dsubmat(AllP, &P[0], nAll, nP, &rowAll[0], NULL);
          
              na_N = DIMS(VECTOR_ELT(a, N1-1));
              ra_N = INTEGER(VECTOR_ELT(a, N1-1));
              
              double summation=sum(P, ra_N, nP, na_N);
              if(summation==0.0){
              goto end2;
              }


// printf("%.9f\t",summation);
// printf("%d\t",N1-rMAX[0]);
//printf("%.9f\t",logfactorial(N1-rMAX[0])+log(summation));  
//printf("%.9f\t",pweight_c);
//   printf("%f\t",logprob_nonmut_sum);
//printf("%f\t",pow(e,logprob_nonmut_sum));        
          //    double temp1=summation*pweight_c*prob_nonmut_sum;
              double logprob_c=logfactorial(N1-rMAX[0])+log(summation)+log(pweight_c)+logprob_nonmut_sum;
//printf("%s\t","temp");
//printf("%.9f\n",temp);

              double prob_c=pow(e,logprob_c);
//printf("%s\t","prob_c");
//printf("%.9f\n",prob_c);
     //         double temp = temp1*temp2*temp3;//sum(P, ra_N, nP, na_N)*pweight_c*pow(1.0/nY[1],N-N1)*factorial(N1-rMAX[0]);
             
              loglik[i] += prob_c;
  //            loglik[i] += sum(P, ra_N, nP, na_N)*pweight_c*pow(1.0/nY[1],N-N1)*factorial(N1-rMAX[0]);
               end2: ;
               }
               end: ;
//printf("jump!\n");
          }
// printf("loglik\t");
//  printf("%.10f\n",loglik[i]);

if(loglik[i]==0.0){
printf("%d\t",i);
//for(int a=0;a<N;a++){printf("%f\t",prob_i[a]);}
//printf("Y\n");
//for(int a=0;a<nY[1];a++){printf("%d\t",rY[i+a*nY[0]]);}
//printf("\n");
}else{    
      loglik[i]=log(loglik[i]);
      }
      }
  
    
  *rans = -dvecsum(loglik, nY[0]);

  UNPROTECT(1);

  return ans;
}

void chindex(int* seq, int nseq, int a, int omit)
{
  int i, j;

  j = 0;
  for (i = 0; i < nseq; i++) {
    if (i == omit)
      j++;
    seq[i] = j + a;
    j++;
  }
}


void makeimat(int *old, int *new, int *no, int *nn,
	      int radd, int cadd, int romit, int comit)
{
  int i, j, r, c;
  for (i = 0, r = 0; i < nn[0]; i++, r++) {
    if (i == romit)
      r++;
    for (j = 0, c = 0; j < nn[1]; j++, c++) {
      if (j == comit)
	c++;
      new[i + j*nn[0]] = old[r + radd + (c + cadd)*no[0]];
    }
  }
}

void rowreduce(double *x, int nx)
{
  int i;
  for (i = 0; i < nx; i++)
    x[i] -= x[nx-1];
}

/*SEXP predict(SEXP Yj, SEXP AllPesti, SEXP a,SEXP prob,SEXP MAX)  //function to caculate P(A<B),P(A>B),P(A=B)

{
    int nYj, *na_N,*nAll;
    int *rYj,*rMAX, *ra_N;
    int N,N1;
    
    double *rAllPesti, *rans, *rprob;
    long double pweight_c;
    
    SEXP ans;
    PROTECT(ans = allocVector(REALSXP, 3)); //Probability of P(A=B),P(A<B),P(A>B)
    rans = REAL(ans);
    
 
    nYj	 = length(Yj);
    rYj = INTEGER(Yj);
    N=rowsum(rYj,nYj);
    nAll = DIMS(AllPesti);
    rMAX=INTEGER(MAX);
    int nP[2];
    
    rAllPesti = REAL(AllPesti);
    rprob=REAL(prob);
    
    int nunsign=0;
    int nc=0;
    count(rprob,N,&nunsign,&nc);
    
    long double jointprob[3];  //P(Y_i,j intersection A<B)
    
    for (int i=0;i<3;i++ ) {
        jointprob[i] += 0.0;
    }

    
    for(int c=0;c<pow(2,nunsign);c++){
        int rcausal[N];
        causind(&rcausal[0],rprob,N,c);
        pweight_c=conprob(&rcausal[0],rprob,N);
        int newYj[nYj];
        caustable(rcausal, rYj , newYj, 0, 1,nYj);
        
        N1=rowsum(&newYj[0],nYj);   //N1: number of causal mut corresponding to c
        
        if (N1>0){
            nP[0] = N1;
            nP[1] = N1;
            double P[nP[0]*nP[1]];
            int rowAll[nP[0]];
            
            
            printf("rporb\t");
            for (int i=0; i<N; i++) {
                printf("%f\t",rprob[i]);
            }
            printf("\n");
            
            printf("N\t%d\n",N);
            
            
            printf("pweight\t%Lf\n",pweight_c);
            
            printf("rYj\t");
            for (int i=0; i<nYj; i++) {
                printf("%d\t",rYj[i]);
            }
            printf("\n");
            
            printf("rcausal\t");
            for (int i=0; i<N; i++) {
                printf("%d\t",rcausal[i]);
            }
            printf("\n");
            
            printf("newY\t");
            for (int i=0; i<nYj; i++) {
                printf("%d\t",newYj[i]);
            }
            printf("\n");
            
            int l = 0;
            for (int j = 0; j < nYj; j++)
                for (int k = 0; k < newYj[j]; k++) {
                    rowAll[l++] = j;
                }
            
            dsubmat(rAllPesti, P, nAll, nP, rowAll, NULL);
            
            printf("rowAll\t");
            for (int i=0; i<N1; i++) {
                printf("%d\t",rowAll[i]);
            }
            printf("\n");
            
            
            na_N = DIMS(VECTOR_ELT(a, N1-1));
            ra_N = INTEGER(VECTOR_ELT(a, N1-1));
            int arows[3];
            int partnum=na_N[0]/N1;
            arows[0]=partnum*newYj[0];
            printf("arrow\t%d\t",arows[0]);
            
            for( int i =1;i<3;i++){
                arows[i]=(na_N[0]/N1)*newYj[i]+arows[i-1];
                printf("%d\t",arows[i]);
            }
            
            
            
            long double temp[3];
            for(int i=0;i<3;i++){
                temp[i]=0.0;
            }
            
            
            
            temp[0]=partsum(P, ra_N, nP, na_N,0,arows[0])*pweight_c*pow(1.0/3,N-N1);
            
            temp[1]=partsum(P, ra_N, nP, na_N,arows[0],arows[1])*pweight_c*pow(1.0/3,N-N1);
            temp[2]=partsum(P, ra_N, nP, na_N,arows[1],arows[2])*pweight_c*pow(1.0/3,N-N1);
            printf("temp\t");
            for (int i=0; i<3; i++) {
                printf("%Lf\t",temp[i]);
            }
            printf("\n");
            
             temp[0]=factorial(temp[0],N1-rMAX[0]);
             temp[1]=factorial(temp[1],N1-rMAX[0]);
             temp[2]=factorial(temp[2],N1-rMAX[0]);
             
             printf("temp\t");
             for (int i=0; i<3; i++) {
             printf("%Lf\t",temp[i]);
             }
             printf("\n\n");
             
             for (int i=0;i<3;i++ ) {
             jointprob[i] += temp[i];
             }

        }
    }
    
    for( int i=0;i<3;i++){
        rans[i]=jointprob[i]*rAllPesti[i];
    }
    
    double sumrans=dvecsum(rans,3);
    for( int i=0;i<3;i++){
        rans[i]=rans[i]/sumrans;
    }
    
//    rans[0]=1.0;
//    rans[1]=1.0;
//    rans[2]=1.0;

    UNPROTECT(1);
    
    return ans;
    
} */

SEXP prob(SEXP Y, SEXP prob)

{
    int *nY, *nP;
    int *rY;
    int i, k, N,c;
    
    double *rans, *prob_nonmut;
    double *prob_i;
    long double pweight_c,sum_prob;
    //  int *newY_i;
    
    
    
    nY	 = DIMS(Y);
    nP = (int *)R_alloc(2, sizeof(int));
    rY = INTEGER(Y);
    
    
    SEXP ans;
    PROTECT(ans = allocVector(REALSXP,nY[1]));
    rans = REAL(ans);


    
    prob_nonmut =(double *)R_alloc(nY[1], sizeof(double));
    for (i = 0; i < nY[1]; i++){
        prob_nonmut[i]=0.0;
    }
  
    
     for (i = 0; i < nY[0]; i++) {
        N = irowsum(rY, nY, i);
        prob_i=REAL(VECTOR_ELT(prob,i));
        int nc=0;
        int nunsign=0;
        count(prob_i,N,&nunsign,&nc);
        
        for(c=0;c<pow(2,nunsign);c++){
            int rcausal[N];
            causind(&rcausal[0],prob_i,N,c);
            pweight_c=conprob(&rcausal[0],prob_i,N);
   //         printf("%Lf\t",pweight_c);

      
         int newY_i[nY[1]*2];
         caustable(rcausal, rY , newY_i, i, nY[0],nY[1] );
            
       
            for(k=0;k<nY[1];k++){
                prob_nonmut[k] += newY_i[k]*pweight_c;
            }
        }

    }
  
 
    sum_prob=0.0;
    
    for(k=0;k<nY[1];k++){
        sum_prob += prob_nonmut[k];
    }
 
 
    printf("%Lf\t",sum_prob);

    for(k=0;k<nY[1];k++){
        rans[k] = prob_nonmut[k]/sum_prob;
    }


    UNPROTECT(1);
    
    return ans;
}


/*    for(int c=0;c<pow(2,nunsign);c++){
        int rcausal[N];
        causind(&rcausal[0],rprob,N,c);
    
        
        printf("rporb\t");
        for (int i=0; i<N; i++) {
            printf("%f\t",rprob[i]);
        }
        printf("\n");
        
        printf("N\t%d\n",N);

 
        printf("pweight\t%Lf\n",pweight_c);
        
        printf("rYj\t");
        for (int i=0; i<nYj; i++) {
            printf("%d\t",rYj[i]);
        }
        printf("\n");

        printf("rcausal\t");
        for (int i=0; i<N; i++) {
            printf("%d\t",rcausal[i]);
        }
        printf("\n");

        for (int i=0; i<nYj; i++) {
            printf("%d\t",newYj[i]);
        }
        
        N1=rowsum(&newYj[0],nYj);   //N1: number of causal mut corresponding to c
        nP[0] = N1;
        nP[1] = N1;
        double P[nP[0]*nP[1]];
        int rowAll[nP[0]];
        
        int l = 0;
        for (int j = 0; j < nYj; j++)
            for (int k = 0; k < newYj[j]; k++) {
                rowAll[l++] = j;
            }
        
        dsubmat(rAllPesti, P, nAll, nP, rowAll, NULL);
        
        na_N = DIMS(VECTOR_ELT(a, N1-1));
        ra_N = INTEGER(VECTOR_ELT(a, N1-1));
        int arows[3];
        for( int i =0;i<3;i++)
            arows[i]=na_N[i]/N1*newYj[i];
        
        long double temp[3];
        
        temp[0]=partsum(P, ra_N, nP, na_N,1,arows[0])*pweight_c*pow(1.0/3,N-N1);
        temp[0]=factorial(temp[0],N1-rMAX[0]);
        temp[1]=partsum(P, ra_N, nP, na_N,arows[0],arows[1])*pweight_c*pow(1.0/3,N-N1);
        temp[1]=factorial(temp[1],N1-rMAX[0]);
        temp[2]=partsum(P, ra_N, nP, na_N,arows[1],arows[2])*pweight_c*pow(1.0/3,N-N1);
        temp[2]=factorial(temp[2],N1-rMAX[0]);
            
        double jointprob[3];  //P(Y_i,j intersection A<B)
            
        for (int i=0;i<3;i++ ) {
                jointprob[i] += temp[i];
        }
        
        
        for( int i=0;i<3;i++){
            rans[i]=rans[i]*rAllPesti[i];
        }
        for( int i=0;i<3;i++){
            rans[i]=dvecsum(rans,3);
        }

        
        
        
   }


  
    
 //   rans[0]=1.0;
 //   rans[1]=1.0;
 //   rans[2]=1.0;
    UNPROTECT(1);
    
    return ans;

    
} */




/*
SEXP gradient(SEXP x, SEXP AllPtrue, SEXP divsel,
	      SEXP Y, SEXP a, SEXP MAX, SEXP prob,SEXP aa)
{
  int nx, ndivsel, nAll[2], nY[2], nP[2], nP1[2], na_N[2], na_N1[2], ngrad,N1;
  int *rY, *rdivsel, *rMAX, *ra_N, *ra_N1;
  int i, j, k, l, cnt, N;

  double *rx, *rAllPtrue, *rans;
  double *AllP, *lik, *grad,*prob_i,pweight_c;
 
  SEXP ans;

  nx = length(x);
  ndivsel=length(divsel);
  nAll[0] = DIMS(AllPtrue)[0]; 
  nAll[1] = DIMS(AllPtrue)[1];
  nY[0] = DIMS(Y)[0]; 
  nY[1] = DIMS(Y)[1];

  rx = REAL(x);
  rAllPtrue = REAL(AllPtrue);
  rdivsel = INTEGER(divsel);
  rY = INTEGER(Y);
  rMAX = INTEGER(MAX);

  lik = (double *)R_alloc(nY[0], sizeof(double));

  AllP = (double *)R_alloc(nAll[0]*nAll[1], sizeof(double));
  for (i = 0; i < nAll[0]; i++)
    for (j = 0; j < nAll[1]; j++)
      AllP[i + j*nAll[0]] = rAllPtrue[i + j*nAll[0]];
 

  for(j=0; j<ndivsel ;j++)
    {
      for (i = 0; i < nAll[0] - 1; i++)
	AllP[i + (rdivsel[j] - 1)*nAll[0]] = rx[i];
      AllP[i + (rdivsel[j]-1)*nAll[0]] = 1 - dvecsum(rx, nx);
    }

  ngrad = nY[1];
  grad = (double *)R_alloc(ngrad, sizeof(double));
  for (i = 0; i < ngrad; i++)
    grad[i] = 0.0;

  PROTECT(ans = allocVector(REALSXP, ngrad - 1));
  rans = REAL(ans);
    


  for (i = 0; i < nY[0]; i++) {
      
      lik[i]=0.0;
      N = irowsum(rY, nY, i);
      prob_i=REAL(VECTOR_ELT(prob,i));
      int nc=0;
      int nunsign=0;
      count(prob_i,N,&nunsign,&nc);
      
      double grad_i[ngrad];
    //  grad_i=(double *)R_alloc(ngrad, sizeof(double));
      for (int m = 0; m < ngrad; m++)
          grad_i[m] = 0.0;
      
      for(int c=0;c<pow(2,nunsign);c++){
          int rcausal[N];
          causind(&rcausal[0],prob_i,N,c);
          pweight_c=conprob(&rcausal[0],prob_i,N);
          
          int newY_i[nY[1]];
          caustable(rcausal, rY , newY_i, i, nY[0],nY[1] );
          
          N1=rowsum(&newY_i[0],nY[1]);
          if(N1==0){
              lik[i] += pow(1.0/nY[1],N)*pweight_c;
          }else{
              nP[0] = N1;
              nP[1] = N1;
              double P[nP[0]*nP[1]];
              int rowAll[nP[0]];
              
              l = 0;
              for (j = 0; j < nY[1]; j++)
                  for (k = 0; k < newY_i[j]; k++) {
                      rowAll[l++] = j;
                  }
              
              dsubmat(AllP, &P[0], nAll, nP, &rowAll[0], NULL);
              
              na_N[0] = DIMS(VECTOR_ELT(a, N1-1))[0];
              na_N[1] = DIMS(VECTOR_ELT(a, N1-1))[1];
              ra_N = INTEGER(VECTOR_ELT(a, N1-1));
              lik[i] += sum(P, ra_N, nP, na_N)*pweight_c*pow(1.0/nY[1],N-N1)*factorial(N1-rMAX[0]);


    if (rdivsel[0] <= N1) {

      if (N1 > 1) {


		
	nP1[0] = nP[0] - 1;
	nP1[1] = nP[1] - 1;
	//P1 = (double *)R_alloc(nP1[0]*nP1[1], sizeof(double));
	//rowP = (int *)R_alloc(nP1[0], sizeof(int));
	//colP = (int *)R_alloc(nP1[1], sizeof(int));
          double P1[nP1[0]*nP1[1]];
          int rowP[nP1[0]];
          int colP[nP1[0]];

	cnt = 0;
	for (j = 0; j < nP[0]; j++) {
	  cnt++;

	  for(l=rdivsel[0];l<=imin2(N,rdivsel[ndivsel-1]);l++)
	    {
	      chindex(rowP, nP1[0], 0, cnt - 1);
	      chindex(colP, nP1[1], 0, l - 1);
	      dsubmat(P, P1, nP, nP1, rowP, colP);

	      if (l>rMAX[0]) {

		na_N1[0] = DIMS(VECTOR_ELT(a, N1-2))[0];
		na_N1[1] = DIMS(VECTOR_ELT(a, N1-2))[1];
		ra_N1 = INTEGER(VECTOR_ELT(a, N1-2));
		grad_i[rowAll[j]] +=
		  //exp(logsum(P1, ra_N1, nP1, na_N1) - loglik[i])/(N-rMAX[0]);
              sum(P1, ra_N1, nP1, na_N1)*pweight_c*pow(1.0/nY[1],N-N1)*factorial(N1-rMAX[0]);

	      } else
		{
		  na_N1[0] = DIMS(VECTOR_ELT(aa, N1-2))[0];
		  na_N1[1] = DIMS(VECTOR_ELT(aa, N1-2))[1];
		  ra_N1 = INTEGER(VECTOR_ELT(aa, N1-2));
		  grad_i[rowAll[j]] +=
            sum(P1, ra_N1, nP1, na_N1)*pweight_c*pow(1.0/nY[1],N-N1)*factorial(N1-rMAX[0]) ;

		}  
	    }
	}
	rowreduce(grad_i, ngrad);

      } else {
	for (j = 0; j < nY[1]; j++)
	  if (rY[i + j*nY[0]] == 1)
	    grad[j] += (1./P[0])*pweight_c*pow(1.0/nY[1],N-N1);

	rowreduce(grad_i, ngrad);
      }
    }
    
    }
    }
    for (j=0;j< ngrad -1;j++){
        grad_i[j] /=lik[i];
        grad[j] += grad_i[j];
    }
    
  }
    
    for (i = 0; i < ngrad - 1; i++)
    rans[i] = -grad[i];
    
  UNPROTECT(1);
  return ans;
}

*/

