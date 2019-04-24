#include "eMath.h"

void init_vector(eVector* v,int nelems){
  v->data = (double*) malloc(sizeof(double)*nelems);
  v->length = nelems;
}

void init_matrix(eMatrix* m,int rows,int cols){
  m->data = (eVector*) malloc(sizeof(eVector)*rows);
  int i;
  for(i=0;i<rows;i++){
    init_vector(&m->data[i],cols);
  }
  m->rows = rows;
  m->cols = cols;
}

void create_vector(eVector* v,double* data,int nelems){
  v->data = data;
  v->length = nelems;
}

void print_vector(eVector v){
  int i;
  printf("[");
  for(i=0;i<v.length;i++){
    printf("%lf ",v.data[i]);
  }
  printf("]\n");
}

void create_matrix(eMatrix* m,double** data,int rows,int cols){
  int i,j;
  eVector* v_data = (eVector*) malloc(sizeof(eVector)*rows);
  for(i=0;i<rows;i++){
    create_vector(v_data+i,data[i],cols);
  }
  m->data = v_data;
  m->rows = rows;
  m->cols = cols;
}

void print_matrix(eMatrix m){
  int i;
  printf("[\n");
  for(i=0;i<m.rows;i++){
    print_vector(m.data[i]);
  }
  printf("]\n");
}

void random_vector(eVector *v){
  srand(time(NULL));
  int i;
  for(i=0;i<v->length;i++){
    v->data[i] = rand()%1000/1000.0;
  }
}

void random_matrix(eMatrix *m){
  srand(time(NULL));
  int i,j;
  for(i=0;i<m->rows;i++){
    for(j=0;j<m->cols;j++){
      m->data[i].data[j] = rand()%1000/1000.0;
    }
  }
}

void ones_vector(eVector *v){
  int i;
  for(i=0;i<v->length;i++){
    v->data[i] = 1.0;
  }
}

void ones_matrix(eMatrix *m){
  int i,j;
  for(i=0;i<m->rows;i++){
    for(j=0;j<m->cols;j++){
      m->data[i].data[j] = 1.0;
    }
  }
}

void sum_vector(eVector* ans,eVector a,eVector b){
  int i;
  init_vector(ans,a.length);
  for(i=0;i<a.length;i++){
    ans->data[i] = a.data[i] + b.data[i];
  }
}

void sum_matrix(eMatrix* ans,eMatrix a,eMatrix b){
  int i,j;
  int rows=a.rows,cols=a.cols;
  init_matrix(ans,rows,cols);
  for(i=0;i<rows;i++){
    for(j=0;j<cols;j++){
      ans->data[i].data[j] = a.data[i].data[j] + b.data[i].data[j];
    }
  }
}

void dif_vector(eVector* ans,eVector a,eVector b){
  int i;
  init_vector(ans,a.length);
  for(i=0;i<a.length;i++){
    ans->data[i] = a.data[i] - b.data[i];
  }
}

void dif_matrix(eMatrix* ans,eMatrix a,eMatrix b){
  int i,j;
  int rows=a.rows,cols=a.cols;
  init_matrix(ans,rows,cols);
  for(i=0;i<rows;i++){
    for(j=0;j<cols;j++){
      ans->data[i].data[j] = a.data[i].data[j] - b.data[i].data[j];
    }
  }
}

void mul_num_vector(eVector *ans,double num,eVector a){
  int i;
  init_vector(ans,a.length);
  for(i=0;i<a.length;i++){
    ans->data[i] = num*a.data[i];
  }
}

void mul_num_matrix(eMatrix *ans,double num,eMatrix a){
  int i,j;
  int rows=a.rows,cols=a.cols;
  init_matrix(ans,rows,cols);
  for(i=0;i<rows;i++){
    for(j=0;j<cols;j++){
      ans->data[i].data[j] = num*a.data[i].data[j];
    }
  }
}

void mul_matrix(eMatrix* ans,eMatrix a,eMatrix b){
  int i,j;
  init_matrix(ans,a.rows,a.cols);
  for(i=0;i<a.rows;i++){
    for(j=0;j<a.cols;j++){
      ans->data[i].data[j] = a.data[i].data[j]*b.data[i].data[j];
    }
  }
}

void dot_vector(double* ans,eVector a,eVector b){
  int i;
  *ans = 0;
  for(i=0;i<a.length;i++){
    *ans = *ans + a.data[i]*b.data[i];
  }
}

void dot_matrix(eMatrix* ans,eMatrix a,eMatrix b){
  int i,j,k;
  init_matrix(ans,a.rows,b.cols);
  double s;
  for(i=0;i<a.rows;i++){
    for(j=0;j<b.cols;j++){
      s=0;
      for(k=0;k<a.cols;k++){
	s = s + a.data[i].data[k]*b.data[k].data[j];
      }
      ans->data[i].data[j] = s;
    }
  }
}

void trans_matrix(eMatrix* ans,eMatrix a){
  int i,j;
  init_matrix(ans,a.cols,a.rows);
  for(i=0;i<a.rows;i++){
    for(j=0;j<a.cols;j++){
      ans->data[j].data[i] = a.data[i].data[j];
    }
  }
}

void norm_vector(double *ans,eVector v){
  int i;
  *ans = 0;
  for(i=0;i<v.length;i++){
    *ans = *ans + v.data[i]*v.data[i];
  }
  *ans = sqrt(*ans);
}

void normalize_vector(eVector *out,eVector v){
  double n;
  norm_vector(&n,v);
  mul_num_vector(out,1/n,v);
}

void diag_matrix(eVector *out,eMatrix m){
  int i;
  init_vector(out,m.rows);
  for(i=0;i<m.rows;i++){
    out->data[i] = m.data[i].data[i];
  }
}

void LU_fact_matrix(eMatrix *L,eMatrix *U,eMatrix a){
  init_matrix(L,a.rows,a.cols);
  init_matrix(U,a.rows,a.cols);
  int i,j,k;
  double s;
  for(i=0;i<a.rows;i++){
    U->data[i].data[i] = 1;
    for(k=i;k<a.rows;k++){
      s = 0;
      for(j=0;j<i;j++){
	s = s + L->data[k].data[j]*U->data[j].data[i];
      }
      L->data[k].data[i] = a.data[k].data[i] - s;
    }
    for(k=i+1;k<a.rows;k++){
      s = 0;
      for(j=0;j<i;j++){
	s = s + L->data[i].data[j]*U->data[j].data[k];
      }
      U->data[i].data[k] = (a.data[i].data[k] - s)/L->data[i].data[i];
    }
  }
}

void QR_fact_matrix(eMatrix *Q,eMatrix *R,eMatrix a){
  init_matrix(Q,a.rows,a.cols);
  init_matrix(R,a.rows,a.cols);
  normalize_vector(&(R->data[0]), a.data[0]);
  int i,j;
  eVector aux,aux1;
  double c;
  for(i=0;i<a.rows;i++){
    init_vector(&aux,a.rows);
    init_vector(&aux1,a.rows);
    for(j=0;j<i;j++){
      dot_vector(&(Q->data[i].data[j]),a.data[i],R->data[j]);
      mul_num_vector(&aux,Q->data[i].data[j],R->data[j]);
      sum_vector(&aux1,aux,aux1);
    }
    dif_vector(&(R->data[i]),a.data[i],aux1);
    norm_vector(&(Q->data[i].data[i]),R->data[i]);
    normalize_vector(&(R->data[i]),R->data[i]);
  }
}

void det_matrix(double *ans,eMatrix m){
  eMatrix L,U;
  LU_fact_matrix(&L,&U,m);
  int i;
  *ans = 1;
  for(i=0;i<m.rows;i++){
    *ans = (*ans)*L.data[i].data[i];
  }
}

void str_2_vector(eVector *out,char* s,int length){
  init_vector(out,length);
  char *t;
  char dlims[3] = ", ";
  t = strtok(s,dlims);
  int i;
  for(i=0;i<length;i++){
    out->data[i] = atof(t);
    t = strtok(NULL,dlims);
  }
}

void str_2_matrix(eMatrix *out,char* s,int rows,int cols){
  init_matrix(out,rows,cols);
  char* t;
  char dlims[2] = ";";
  t = strsep(&s,dlims);
  int i,j;
  for(i=0;i<rows;i++){
    str_2_vector(&(out->data[i]),t,cols);
    t = strsep(&s,dlims);
  }
}

void load_vector(eVector *out,char* dir,int length){
  init_vector(out,length);
  FILE *f;
  f = fopen(dir,"r");
  char s[200];
  fgets(s,200,f);
  str_2_vector(out,s,length);
  fclose(f);
}
void load_matrix(eMatrix *out,char* dir,int rows,int cols){
  init_matrix(out,rows,cols);
  FILE *f;
  int i;
  f = fopen(dir,"r");
  char s[200];
  for(i=0;i<rows;i++){
    fgets(s,200,f);
    str_2_vector(&(out->data[i]),s,cols);
  }
  fclose(f);
}

void save_vector(char* dir,eVector v){
  int i;
  FILE* f;
  f = fopen(dir,"w");
  for(i=0;i<v.length;i++){
    fprintf(f,"%lf ",v.data[i]);
  }
  fprintf(f,"\n");
  fclose(f);
}

void save_matrix(char* dir,eMatrix m){
  int i,j;
  FILE* f;
  f = fopen(dir,"w");
  for(i=0;i<m.rows;i++){
    for(j=0;j<m.cols;j++){
      fprintf(f,"%lf ",m.data[i].data[j]);
    }
    fprintf(f,"\n");
  }
  fprintf(f,"\n");
  fclose(f);
}

void solve_linear_system(eVector* ans,eMatrix A,eVector b){
  eMatrix L,U;
  LU_fact_matrix(&L,&U,A);
  int i,j;
  double s;
  eVector z;
  init_vector(ans,b.length);
  init_vector(&z,b.length);
  for(i=0;i<b.length;i++){
    s = 0;
    for(j=0;j<i;j++){
      s = s + L.data[i].data[j]*z.data[j];
    }
    z.data[i] = (b.data[i] - s)/(L.data[i].data[i]);
  }
  for(i=b.length-1;i>=0;i--){
    s = 0;
    for(j=i+1;j<b.length;j++){
      s = s + U.data[i].data[j]*ans->data[j];
    }
    ans->data[i] = (z.data[i] - s);
  }
}

void inv_matrix(eMatrix* ans,eMatrix m){
  eMatrix L,U;
  LU_fact_matrix(&L,&U,m);
  init_matrix(ans,m.rows,m.cols);
  eVector z;
  init_vector(&z,m.rows);
  int i,j,k;
  double s;
  for(i=0;i<m.rows;i++){
   for(j=0;j<m.rows;j++){
    s = 0;
    for(k=0;k<j;k++){
      s = s + L.data[j].data[k]*z.data[k];
    }
    if(i==j){
      z.data[j] = (1 - s)/(L.data[j].data[j]);
    }
    else{
      z.data[j] = ( - s)/(L.data[j].data[j]);
    }
   }
   for(j=m.rows-1;j>=0;j--){
     s = 0;
     for(k=j+1;k<m.rows;k++){
       s = s + U.data[j].data[k]*ans->data[k].data[i];
     }
     ans->data[j].data[i] = (z.data[j] - s);
   } 
  }
}

void eigen_matrix(eVector* e_values,eMatrix m,int N){
  eMatrix Q,R;
  init_matrix(&Q,m.rows,m.cols);
  init_matrix(&R,m.rows,m.cols);
  int i;
  Q = m;
  for(i=0;i<N;i++){
    QR_fact_matrix(&Q,&R,Q);
    dot_matrix(&Q,R,Q);
  }
  diag_matrix(e_values,Q);
}

void merge_data(eVector *ans,eVector a,eVector b){
  init_vector(ans,a.length+b.length);
  int i_0=0,j_0=0,i_f=b.length - 1,j_f=a.length - 1;
  while(i_0 <= i_f){
    while(b.data[i_0] > a.data[j_0] && j_0 <= j_f){
      ans->data[i_0 + j_0] = a.data[j_0];
      j_0++;
    }
    ans->data[i_0 + j_0] = b.data[i_0];
    i_0++;
    if(i_0>i_f){
      break;
    }
    while(b.data[i_f] < a.data[j_f] && j_0 <= j_f){
      ans->data[i_f + j_f+1] = a.data[j_f];
      j_f--;
    }
    ans->data[i_f + j_f + 1] = b.data[i_f];
    i_f--;
  }
  while(j_0<=j_f){
    ans->data[i_0+j_0+1] = a.data[j_0];
    j_0++;
  }
}

void extract_data(eVector *ans,eVector a,int i,int j){
  init_vector(ans,j-i+1);
  int x;
  for(x=i;x<=j;x++){
    ans->data[x-i] = a.data[x];
  }
}
void sort_data(eVector* out,eVector v){
  int N;
  eVector aux1,aux2;
  if(v.length==1){
    *out = v;
  }
  else{
    N = v.length/2 - 1;
    extract_data(&aux1,v,0,N);
    extract_data(&aux2,v,N+1,v.length-1);
    sort_data(&aux1,aux1);
    sort_data(&aux2,aux2);
    merge_data(out,aux1,aux2);
  }
}

void average_data(double *ans,eVector v)
{
  *ans = 0;
  int i;
  for(i=0;i<v.length;i++){
    *ans = *ans + v.data[i];
  }
  *ans = (*ans)/v.length;
}

void median_data(double *ans,eVector v){
}

void mode_data(double *ans,eVector v){
}

void std_dev_data(double *ans,eVector v){
  *ans = 0;
  int i;
  double avg;
  average_data(&avg,v);
  for(i=0;i<v.length;i++){
    *ans = *ans + (v.data[i]-avg)*(v.data[i]-avg);
  }
  *ans = sqrt((*ans)/v.length); 
}

void least_square_data(eVector *ans,eVector x,eVector y){
  eMatrix A;
  eVector b;
  init_matrix(&A,2,2);
  init_vector(&b,2);
  double sx=0,sy=0;
  int i;
  for(i=0;i<x.length;i++){
    sx = sx + x.data[i];
    sy = sy + y.data[i];
  }
  dot_vector(&(A.data[0].data[0]),x,x);
  dot_vector(&(b.data[0]),x,y);
  A.data[0].data[1] = sx;
  A.data[1].data[0] = sx;
  A.data[1].data[1] = x.length;
  b.data[1] = x.length;
  solve_linear_system(ans,A,b);
}

void gaussian_dist(double *ans,double avg,double std_dev,double x){
  *ans = (1/(sqrt(2*PI)*std_dev)) * exp(-(x-avg)*(x-avg)/(2*std_dev*std_dev));
}

void diff_function(double* ans,void fun_f(double*,double),double x0,double h){  
  double x_h;
  fun_f(&x_h,x0 + h);
  fun_f(ans,x0);
  *ans = (x_h - *ans)/h;
}

void integrate_function(double *ans,void fun_f(double*,double),double x_0,double x_1,int N){
  double x,h,tmp;
  int i;
  h = (x_1 - x_0)/N;
  x = x_0;
  *ans = 0;
  for(i=0;i<N;i++){
    fun_f(&tmp,x);
    *ans = *ans + tmp*h;
    x = x + h;
  }
}

void linspace_mesh(eVector *out,double x0,double xf,int N){
  double h = (xf - x0)/(N-1);
  int i;
  init_vector(out,N);
  out->data[0] = x0;
  for(i=0;i<N;i++){
    out->data[i] = out->data[i-1] + h;
  }
}

void ode1_solver(eMatrix *out,void fun_ode(double*,double,double),double init_cond,eVector x_mesh){
  init_matrix(out,x_mesh.length,2);
  int i;
  out->data[0].data[0] = init_cond;
  fun_ode(&(out->data[0].data[1]),init_cond,x_mesh.data[0]);
  for(i=1;i<x_mesh.length;i++){
    out->data[i].data[0] = out->data[i-1].data[0] + out->data[i-1].data[1]*(x_mesh.data[i] - x_mesh.data[i-1]);
    fun_ode(&(out->data[i].data[1]),out->data[i].data[0],x_mesh.data[i]);
  }
}

