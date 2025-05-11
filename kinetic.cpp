#include<cstdio>
#include<cmath>
#include<iostream>
#include<ctime>
#include<random> 
#include<cstdlib>
#include <fstream> 
#include <omp.h>
#define OMP_NUM_THREADS 4

const int g_nNumberOfThreads = 4;
using namespace std;
/*ПАРАМЕТРЫ*/
#define a 0.
#define ab 50.
#define b 100.
#define TYPE double
#define EEE0 10.

#define po1 1. //плотность слева
#define po2 1. //плотность справа

#define B 1.
TYPE W1=1;
TYPE W2=1.5;
TYPE L=0.5;

TYPE *E_pos_e;					
TYPE *E_pos_i;
TYPE *pos_e;					
TYPE *pos_i;
TYPE *V_i;
TYPE *V_e;
TYPE *V_i_av;
TYPE *V_e_av;
TYPE *E;
TYPE *rho_i;
TYPE *phi;
TYPE *rho_e;
TYPE *num_i;
TYPE *num_e;
TYPE *net;
TYPE *net_2;
TYPE *alpha; //прогонка
TYPE *beta1;
int M;
int N;
TYPE dx;
TYPE dt;
int frame;
TYPE m_e;
TYPE m_i;
TYPE q;
TYPE T_i;
TYPE T_e;
TYPE T_max;
int v;
TYPE lambda;
string note;
TYPE rho_ei();
TYPE spl(TYPE x);
TYPE aa[]={0,0,0,0};
int flag1=1;
TYPE wt(TYPE x,TYPE x_L, TYPE x_P){
	TYPE t = (x-x_L)/(x_P-x_L);
	return (2*t*t*t-3*t*t+1);
}
TYPE *spl_x;
TYPE *spl_y;
TYPE fun(TYPE k1, TYPE k2, TYPE d, TYPE x0, TYPE x){
	TYPE Q = k1*x0+k2*(b-x0);
	TYPE c1 = k1 - b/x0;
	TYPE t = (x*b-x0+d)/(2*d);
	if(x<(x0-d)/b) return k1*x*b/Q;
	if(x>(x0+d)/b) return (c1*x+1.)*b/Q;

	return ((k1-k2)*(2*d)*((t*t*t*t)/2 -t*t*t+t)+k1*(x0-d)+2*d*k2*t)/Q;
}
TYPE distr_spl(TYPE x){
	for(int i = 0; i < N-1; ++i){
		if(x<spl_x[i+1] && x>=spl_x[i]) return (spl_x[i+1]-x)*spl_y[i]/(spl_x[i+1]-spl_x[i])+(x-spl_x[i])*spl_y[i+1]/(spl_x[i+1]-spl_x[i]);
	}
	return spl_y[N-1];
}
TYPE distr(TYPE k1, TYPE k2, TYPE d, TYPE x0, TYPE x){
	TYPE Q = k1*x0+k2*(b-x0);
	TYPE c1 = k1 - b/x0;
	if(x<k1*(x0-d)/Q) return Q*x/(k1*b);
	if(x>c1*(x0+d)/Q+b/Q) return Q*x/c1/b - 1/c1;
	if (flag1){
		spl_x = new TYPE[N];
		//TYPE dx_t = (TYPE)((x0+d)/b-(x0-d)/b)/(N-1);
		TYPE dx_t = (TYPE)((x0+d)/b-(x0-d)/b)/(N-1);
		spl_y = new TYPE[N];
		TYPE x_temp =(x0-d)/b; 
		FILE *file; 
        char filename[80];
        char path[80];
        sprintf(path, "data/v%d/kinmodver%d", v,v);
        sprintf(filename, "%sdistr.txt",path);
		file = fopen(filename, "w");
        //rho_i, E, phi
		//printf("%f\t%f\t%f",x_temp,c1*(x0+d)/Q+b/Q,dx_t);
		for (int i = 0; i < N; ++i){
			//x_temp+=dx_t;
			spl_x[i] = x_temp+fabs(x_temp-fun(k1, k2, d, x0, x_temp));
			spl_y[i] = fun(k1, k2, d, x0, x_temp)-fabs(x_temp-fun(k1, k2, d, x0, x_temp));
			//spl_x[i]=x_temp;
			//spl_y[i]=fun(k1, k2, d, x0, x_temp);
			x_temp+=dx_t;
			//x_temp=2*dx_t*i-fun(k1, k2, d, x0, dx_t*i);
			//spl_x[i] = x_temp+fabs(x_temp-fun(k1, k2, d, x0, x_temp));
			//spl_y[i] = fun(k1, k2, d, x0, x_temp)-fabs(x_temp-fun(k1, k2, d, x0, x_temp));
			//printf("%f \t%f \n", spl_x[i],spl_y[i]);
		}
		flag1=0;
		for(int i = 0; i < N; ++i){
		//x_temp=2*dx_t*i-fun(k1, k2, d, x0, dx_t*i);
		fprintf(file, "%f %f", spl_x[i], spl_y[i]); 
		fprintf(file, "\n");
	}
		fclose(file);
		FILE *file1;
		file1 = fopen("distr.gpi", "w");
		fprintf(file1, "set term \"png\"\n"); 
		char filename1[80];
		sprintf(filename1,"kinmodver%d_distr.png",v);
		fprintf(file1, "set output \"%s\"\n",filename1); 
		fprintf(file1, "plot sprintf(\"%s\") using 1:2 with points pt 1 ps 1 title \"distr\"",filename); 
		fclose(file1);
		system("gnuplot distr.gpi");
	}
	return distr_spl(x);
}
int rho_start(){
    TYPE m1 = (po1*(ab));
	TYPE m2 = (po2*(b-ab));
	TYPE m = m1+m2;
	q = m/M;
	TYPE dx_t = (TYPE)1./(M-1);

	for (int i = 0; i < M; ++i){
		pos_e[i] = b*distr(po1,po2,b/200.,b/2.,dx_t*i);
		pos_i[i] = b*distr(po1,po2,b/200.,b/2.,dx_t*i);
	}
	return 0;
}
int loadpar(){
	v=0;
	ifstream in;
	ofstream out;
	string s;
	in.open("parameters.txt");
	getline(in,s);
	in>>M;
	in>>N;
	in>>dt;
	in>>T_i;
	in>>T_e;
	in>>m_i;
	in>>m_e;
	in>>q;
	in>>frame;
	in>>T_max;
	in>>lambda;
	in.close();
	in.open("log.txt");
	while (getline(in, s)){v++;} 
	//in>>v;
	in.close();
	v-=1;
	out.open("log.txt", ios::app);
	out<<"\n"<<v<<"\t"<<M<<"\t"<<N<<"\t"<<dt<<"\t"<<T_i<<"\t"<<T_e<<"\t"<<m_i<<"\t"<<m_e<<"\t"<<q<<"\t"<<frame<<"\t"<<T_max<<"\t"<<lambda<<"\t"<<"";
	out.close();
	char s2[80];
	sprintf(s2, "data/v%d", v);
	sprintf(s2,"mkdir data\\v%d", v);
    //sprintf(temp,"mkdir data/v%d", v);
    system(s2);
	//_mkdir(s2);
	dx = (b/(N-1));
	E_pos_e = new TYPE[M];						
	E_pos_i = new TYPE[M];
	pos_e = new TYPE[M];						
	pos_i = new TYPE[M];
	V_i = new TYPE[M];
	V_e = new TYPE[M];
        
	V_i_av = new TYPE[N+1];
	V_e_av = new TYPE[N+1];
        
	E = new TYPE[N];
	phi= new TYPE[N+1];
	alpha= new TYPE[N+1];
	beta1= new TYPE[N+1];
	rho_i = new TYPE[N+1];
	rho_e = new TYPE[N+1];
        
	num_i = new TYPE[N-1];
	num_e = new TYPE[N-1];
	net =  new TYPE[N];
	net_2 =  new TYPE[N-1];
	for(int i=0;i<N;i++){net[i]=i*dx;}
	for(int i=0;i<N-1;i++){net_2[i]=(i+1)*dx-dx/2.;}
	return 0;
}
int initialization(){
	loadpar();
	rho_start();
	//normal_distribution<double> dist1(0.0, T_i);
	normal_distribution<double> dist2(0.0, T_e);
	minstd_rand0 gen(2);
	normal_distribution<double> dist1(0.0, 1./1000);
	//minstd_rand0 gen(2);
	for (int i = 0; i < M; ++i)
	{
		V_i[i]=0;//dist1(gen);
		V_e[i]=dist2(gen);
		while (fabs(V_e[i])>4*T_e) V_e[i]=dist2(gen);
	}
	double sr=0;
	for (int i = 0; i < M; ++i)
	{
		sr+= V_e[i];
	}
	sr/=M;
	for (int i = 0; i < M; ++i)
	{
		V_e[i]-=sr;
	}
	for (int i = 0; i < M; ++i)
	{if (V_e[i]*dt/dx>1) printf("%f %f\n", V_e[i]*dt,dx);
	}
	for (int i = 0; i < N; ++i)
	{
		E[i]=0;
	}
	for (int i = 0; i < N+1; ++i)
	{
        phi[i]=0;
		rho_i[i]=0;
		rho_e[i]=0;
	}
	rho_ei();
	return 0;
}
TYPE rho(TYPE *rho_m, TYPE *pos_m){
	TYPE ss;
	int j;
	#pragma omp parallel for
	for (int i = 0; i < N+1; ++i)
	{
		rho_m[i]=0;
	}
	#pragma omp parallel for
	for (int i = 0; i < M; ++i)
	{
		if (pos_m[i]<=net_2[0])
		{
			//rho_m[0]+=0.5+pos_m[i]/dx;
			//rho_m[N-2]+=0.5-pos_m[i]/dx;
                        rho_m[1]+=0.5+pos_m[i]/dx;
			rho_m[0]+=0.5-pos_m[i]/dx;
			continue;
		}

		j = (int)((pos_m[i]-dx/2.)/dx);
		ss=(pos_m[i]-dx/2.)/dx - j;

		if(pos_m[i]>=net_2[N-2])
		{
			//rho_m[N-2]+=1-ss;
			//rho_m[0]+=ss;
                        rho_m[N-1]+=1-ss;
			rho_m[N]+=ss;
			continue;
		}
		else{
			rho_m[j+1]+=1-ss;
			rho_m[j+2]+=ss;
		}
	}
	return 0;
}
TYPE V_MAX = 0;
TYPE V_av(TYPE *V_m_av, TYPE *pos_m, TYPE *V_m){
	TYPE ss;
	int j;
	TYPE v_max=0;
	#pragma omp parallel for
	for (int i = 0; i < N+1; ++i)
	{
		V_m_av[i]=0;
	}
	#pragma omp parallel for
	for (int i = 0; i < M; ++i)
	{
		if (v_max<fabs(V_m[i])) v_max=fabs(V_m[i]);
		if (pos_m[i]<net_2[0])
		{
			//V_m_av[0]+=V_m[i]*(0.5+pos_m[i]/dx);
			//V_m_av[N-2]+=V_m[i]*(0.5-pos_m[i]/dx);
            V_m_av[1]+=V_m[i]*(0.5+pos_m[i]/dx);
			V_m_av[0]+=V_m[i]*(0.5-pos_m[i]/dx);
			continue;
		}

		j = (int)((pos_m[i]-dx/2.)/dx);
		ss=(pos_m[i]-dx/2.)/dx - j;

		if(pos_m[i]>net_2[N-2])
		{
			//V_m_av[N-2]+=V_m[i]*(1-ss);
			//V_m_av[0]+=V_m[i]*ss;
                        V_m_av[N-1]+=V_m[i]*(1-ss);
			V_m_av[N]+=V_m[i]*ss;
			continue;
		}
		else{
			V_m_av[j+1]+=V_m[i]*(1-ss);
			V_m_av[j+2]+=V_m[i]*ss;
		}
	}
	//if (v_max*dt/dx>1){printf("\nKurant-error v_max = %f", v_max);}
	//else printf("\nv_max = %f\n", v_max);
	if (v_max>V_MAX)V_MAX=v_max;
	return 0;
}

TYPE rho_ei(){
	
	rho(rho_e,pos_e);
	rho(rho_i,pos_i);
	V_av(V_e_av, pos_e, V_e);
	V_av(V_i_av, pos_i, V_i);
	#pragma omp parallel for
	for (int i = 0; i < N+1; ++i)
	{
		if(rho_e[i] == 0) V_e_av[i] = 0;
		else
		V_e_av[i]/=rho_e[i];
		if(rho_i[i] == 0) V_i_av[i] = 0;
		else
		V_i_av[i]/=rho_i[i];
	}
	#pragma omp parallel for
	for (int i = 0; i < N+1; ++i)
	{
		rho_e[i]/=dx;
		rho_i[i]/=dx;
	}
	
	#pragma omp parallel for
	for (int i = 0; i < N+1; ++i)
	{
		rho_e[i]*=q;
		rho_i[i]*=q;
	}
	
	return 0;
}
TYPE spl(TYPE x){
	TYPE s;
	TYPE value=0;
	int j = (int)(x/dx);
	s=x/dx-j;
	if(j>N-1 || j<0) return 0;
	return (s)*E[j+1]+(1-s)*E[j];
	return value;
}
TYPE spl_2(TYPE x, TYPE *mas_to_spl){
	TYPE ss;
	int j;
	if (x<net_2[0])
		{
			//return mas_to_spl[0]*(0.5+x/dx)+mas_to_spl[N-2]*(0.5-x/dx);
                        return mas_to_spl[1]*(0.5+x/dx)+mas_to_spl[0]*(0.5-x/dx);
		}

	j = (int)((x-dx/2.)/dx);
	ss=(x-dx/2.)/dx - j;
	//if (i==temp) printf("\t %d %f", j, ss);
	if(x>net_2[N-2])
		{
			//return mas_to_spl[0]*ss+mas_to_spl[N-2]*(1-ss);
                        return mas_to_spl[N-1]*ss+mas_to_spl[N]*(1-ss);
		}
	else{
			//return mas_to_spl[j+1]*ss+mas_to_spl[j]*(1-ss);
                        return mas_to_spl[j+1]*ss+mas_to_spl[j]*(1-ss);
		}
	return 0;
}
TYPE CHECKINT(){
	TYPE I1=0;
	TYPE I2=0;
	for(int i = 0; i<N-1;i++){
		I1+=E[i]*E[i];
	}
	for(int i = 0; i<M;i++){
		I2+=V_e[i]*V_e[i]+(m_i/m_e)*V_i[i]*V_i[i];
	}
	return I1*dx+I2*q;
}
TYPE CHECKINT1(){
	TYPE I1=0;
	for(int i = 0; i<N-1;i++){
		I1+=E[i]*E[i];
	}
	return I1*dx;
}
TYPE CHECKINT2(){
	TYPE I2=0;
	for(int i = 0; i<M;i++){
		I2+=V_e[i]*V_e[i];
	}
	return I2*q;
}
TYPE CHECKINT3(){
	TYPE I3=0;
	for(int i = 0; i<M;i++){
		I3+=V_i[i]*V_i[i];
	}
	return I3*q*(m_i/m_e);
}
void savedraw_values(char* path){
	FILE *file; 
        char filename[80];
        char filename2[80];
        char filename3[80];
        char filename4[80];
        char filename5[80];
        char filename6[80];
        char filename7[80];
        sprintf(filename, "%srho_i.txt",path);
	file = fopen(filename, "w");
        //rho_i, E, phi
	for(int i = 0; i <= N; ++i){
		fprintf(file, "%f %f", (i-1/2.)*dx, rho_i[i]); 
		fprintf(file, "\n");
	}
	fclose(file);
        //V_e_av
	sprintf(filename2, "%sphi.txt",path);
        file = fopen(filename2, "w");
        for(int i = 0; i <= N; ++i){
		fprintf(file, "%f %f", (i-1/2.)*dx, phi[i]); 
		fprintf(file, "\n");
	}
        
        fclose(file);
        sprintf(filename3, "%sE.txt",path);
	file = fopen(filename3, "w");
	//	for(int i = 0; i <= N-1; ++i){
	//	fprintf(file, "%f %f", (i)*dx, spl_2((i)*dx, phi)); 
	//	fprintf(file, "\n");
	//}
         for(int i = 0; i < N; ++i){
		fprintf(file, "%f %f", (i)*dx, E[i]); 
		fprintf(file, "\n");
	}
        fclose(file);
        
        sprintf(filename4, "%spos_i.txt",path);
	file = fopen(filename4, "w");
    for(int i = 0; i < M; i+=1){
		fprintf(file, "%f %f", pos_i[i], V_i[i]); 
		fprintf(file, "\n");
	}
        fclose(file);
        
        sprintf(filename5, "%sV_e_av.txt",path);
        file = fopen(filename5, "w");
        //rho_i, E, phi
	for(int i = 0; i < N+1; ++i){
		fprintf(file, "%f %f", (i-1/2.)*dx, V_e_av[i]); 
		fprintf(file, "\n");
	}
	fclose(file);
        
    sprintf(filename6, "%sV_i_av.txt",path);
    file = fopen(filename6, "w");
        //rho_i, E, phi
	for(int i = 0; i < N+1; ++i){
		fprintf(file, "%f %f", (i-1/2.)*dx, V_i_av[i]); 
		fprintf(file, "\n");
	}
	fclose(file);

	sprintf(filename7, "%sEDlambdaMV_e_av.txt",path);
        file = fopen(filename7, "w");
        //rho_i, E, phi
	//for(int i = 0; i < N+1; ++i){
	//	fprintf(file, "%f %f", (i)*dx, E[i]+spl_2(i*dx,V_e_av)); 
	//	fprintf(file, "\n");
	//}
	fclose(file);
        

	FILE *file1;
	file1 = fopen("temp.gpi", "w");
	fprintf(file1, "set term \"png\"\n"); 
	char filename1[80];
	sprintf(filename1,"kinmodver%d.png",v);
	fprintf(file1, "set output \"%s\"\n",filename1); 
	fprintf(file1, "plot sprintf(\"%s\") using 1:2 with lines lw 2 title \"rho_i\", sprintf(\"%s\") using 1:2 with lines lw 2 title \"phi\", sprintf(\"%s\") using 1:2 with lines lw 2 title \"E\",sprintf(\"%s\") using 1:2 with points pt 0 title \"pos_i\",sprintf(\"%s\") using 1:2 with lines lw 2 title \"V_e_av\",sprintf(\"%s\") using 1:2 with lines lw 2 title \"V_i_av\",sprintf(\"%s\") using 1:2 with lines lw 2 title \"res\"",filename,filename2,filename3, filename4,filename5, filename6, filename7); 
	fclose(file1);
	system("gnuplot temp.gpi");

	//FILE *file1;
	file1 = fopen("pos_i.gpi", "w");
	fprintf(file1, "set term \"png\"\n"); 
	//char filename1[80];
	sprintf(filename1,"kinmodver%d_pos_i.png",v);
	fprintf(file1, "set output \"%s\"\n",filename1); 
	fprintf(file1, "plot sprintf(\"%s\") using 1:2 with points pt 0 title \"pos_i\"",filename4); 
	fclose(file1);
	system("gnuplot pos_i.gpi");

	sprintf(filename4, "%spos_e.txt",path);
	file = fopen(filename4, "w");
    for(int i = 0; i < M; i+=1){
		fprintf(file, "%f %f", pos_e[i], V_e[i]); 
		fprintf(file, "\n");
	}
    fclose(file);

    file1 = fopen("pos_e.gpi", "w");
	fprintf(file1, "set term \"png\"\n"); 
	//char filename1[80];
	sprintf(filename1,"kinmodver%d_pos_e.png",v);
	fprintf(file1, "set output \"%s\"\n",filename1); 
	fprintf(file1, "plot sprintf(\"%s\") using 1:2 with points pt 0 title \"pos_e\"",filename4); 
	fclose(file1);
	system("gnuplot pos_e.gpi");
        
}
int main(){
	FILE *file_temp; 
    char filename_temp[80];

    sprintf(filename_temp, "timesc%d.txt",v);
	file_temp = fopen(filename_temp, "w");
	int iters=0;
	initialization();
	double T=0;
	double graf_T=0;
	double graf_dt=frame*dt;
	double s=0;
    char temp[80];
    char path[80];
    printf("%d\n",v);
        
    sprintf(path, "data/v%d/kinmodver%d", v,v);
	TYPE INT1 = CHECKINT();
	printf("\n");
	TYPE temp1, temp2, kappa;
	temp2 = exp(-(lambda)*dt);
	temp1 = exp(-dt*lambda*(m_e/m_i));
	TYPE TF1, TF2;
	while(T<T_max){
		iters++;
		V_MAX = 0;
		rho_ei();
				rho_e[1]=rho_e[0]+rho_e[1];
				rho_i[1]=rho_i[0]+rho_i[1];
				rho_e[N-1]=rho_e[N]+rho_e[N-1];
				rho_i[N-1]=rho_i[N]+rho_i[N-1];
				rho_e[0]=rho_e[1];
				rho_e[N]=rho_e[N-1];
				rho_i[0]=rho_i[1];
				rho_i[N]=rho_i[N-1];

        
        //for (int i = 0; i <= N; ++i)
		//{
		//	E[i]=E[i]+(spl_2(i*dx, V_e_av)-spl_2(i*dx, V_i_av))*dt;
		//}
		E[0]=EEE0;
		for (int i = 1; i <= N; ++i)
		{
			E[i]=E[i-1]+(-rho_e[i]+rho_i[i])*dx;
		}
		//V_MAX=0;
		//for(int i = 0;i<M;i++){
		//	if (V_MAX<fabs(V_e[i])) V_MAX = fabs(V_e[i]);
		//}
		//dt = dx/V_MAX;
		
		#pragma omp parallel for
		{
		for (int i = 0; i < M; ++i)
		{
			V_i[i]=V_i[i]+(1./m_i)*spl(pos_i[i])*dt;
			V_e[i]=V_e[i]-spl(pos_e[i])*dt;

			//TF1 = spl(pos_i[i])/lambda + spl_2(pos_i[i], V_e_av);
			//TF2 = -spl(pos_e[i])/lambda + spl_2(pos_e[i], V_i_av);
			//V_i[i]=temp1*(V_i[i]-TF1)+TF1;
			//V_e[i]=temp2*(V_e[i]-TF2)+TF2;
			//V_i[i]=V_i[i]+(1./m_i)*spl(pos_i[i])*dt;
			//V_e[i]=V_e[i]-spl(pos_e[i])*dt;
			//V_e[i]=V_e[i]*temp2+(1-temp2)*(-spl(pos_e[i])/(lambda)+spl_2(pos_e[i], V_i_av));
			//V_i[i]=V_i[i]*temp1+(1-temp1)*(((1)*spl(pos_i[i]))/(lambda)+spl_2(pos_i[i], V_e_av));
		}
		}
		
		/*for (int i = 0; i < M; ++i)
		{
			//TF1 = spl(pos_i[i])/lambda + spl_2(pos_i[i], V_e_av);
			//TF2 = -spl(pos_e[i])/lambda + spl_2(pos_e[i], V_i_av);

			//pos_i[i]=pos_i[i]+dt*(TF1)+(1-temp1)*(V_i[i]-TF1)/((lambda*m_e/m_i));
			//pos_e[i]=pos_e[i]+dt*TF2+(1-temp2)*(V_e[i]-TF2)/(lambda);

                        if( pos_i[i]>b) {
                            pos_i[i]=2*(b)-pos_i[i];
                            V_i[i]=-V_i[i];
                        }
                        if( pos_i[i]<0) {
                            pos_i[i]=-pos_i[i];
                            V_i[i]=-V_i[i];
                        }
                        if( pos_e[i]>b) {
                            pos_e[i]=2*(b)-pos_e[i];
                            V_e[i]=-V_e[i];
                        }
                        if( pos_e[i]<0) {
                            pos_e[i]=-pos_e[i];
                            V_e[i]=-V_e[i];
                        }
		}*/
		if (T==0)
		{
			for (int i = 0; i < M; ++i)
		{
			pos_e[i]=pos_e[i]+dt*V_e[i]/2.;
			pos_i[i]=pos_i[i]+dt*V_i[i]/2.;

                        if( pos_i[i]>b) {
                            pos_i[i]=2*(b)-pos_i[i];
                            V_i[i]=-V_i[i];
                        }
                        if( pos_i[i]<0) {
                            pos_i[i]=-pos_i[i];
                            V_i[i]=-V_i[i];
                        }
                        if( pos_e[i]>b) {
                            pos_e[i]=2*(b)-pos_e[i];
                            V_e[i]=-V_e[i];
                        }
                        if( pos_e[i]<0) {
                            pos_e[i]=-pos_e[i];
                            V_e[i]=-V_e[i];
                        }
		}
		}
	else{
		#pragma omp parallel for
		{
		for (int i = 0; i < M; ++i)
		{
			pos_e[i]=pos_e[i]+dt*V_e[i];
			pos_i[i]=pos_i[i]+dt*V_i[i];
                        if( pos_i[i]>b) {
                            pos_i[i]=2*(b)-pos_i[i];
                            V_i[i]=-V_i[i];
                        }
                        if( pos_i[i]<0) {
                            pos_i[i]=-pos_i[i];
                            V_i[i]=-V_i[i];
                        }
                        if( pos_e[i]>b) {
                            pos_e[i]=2*(b)-pos_e[i];
                            V_e[i]=-V_e[i];
                        }
                        if( pos_e[i]<0) {
                            pos_e[i]=-pos_e[i];
                            V_e[i]=-V_e[i];
                        }
		}
		}
	}
		//printf("**");
		if (T>graf_T){
			//fprintf(file_temp, "%f %f", T, V_MAX); 
			//fprintf(file_temp, "\n");
			//v_maxx[(int)(T/dt)] = V_MAX;
			graf_T+=graf_dt;
			printf("\nend_prog V_MAX = %f", V_MAX);
			printf("\n%f %f %f %f %d\n",CHECKINT(), V_MAX,T, dt, iters);
			fprintf(file_temp,"%f %f %f %d\n",CHECKINT(), V_MAX,T, dt, iters);

		}
		T+=dt;
	}
	
	
	
        //fclose(file_temp);
	printf("\nend_prog V_MAX = %f", V_MAX);
        savedraw_values(path);
        printf("End");
        free(E_pos_e);						
	 free(E_pos_i);
	 free(pos_e);						
	 free(pos_i);
	 free(V_i);
	 free(V_e);
	 free(V_i_av);
	 free(V_e_av);
	 free(E);
	 free(phi);
	 free(alpha);
	 free(beta1);
	 free(rho_i);
	 free(rho_e);
	 free(num_i);
	 free(num_e);
	 free(net);
	 free(net_2);
	return 0;
}