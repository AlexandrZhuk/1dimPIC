#include<cmath>
#include<cstdio>
//#include<conio.h>
#include<iostream>
#include<cstdlib>  //библиотеки
#define TYPE double
#define a 0
#define ab 50
#define b 100.

#define M 200000 //количество частиц
#define N 2000 // узлы сетки
#define dx (1./(N-1)) // шаг сетки
#define dt 0.001 //шаг времени

#define po1 3 //плотность слева
#define po2 1. //плотность справа

TYPE W2 = po1;
TYPE W1 = po2;
TYPE L = 0.2;


#define eps 0.000001
#define B 0.001
 //тип данных в вычислениях

/*массивы*/
TYPE qua[M];
TYPE qua2[M]; //положения частиц
TYPE V[M]; // скорости частиц
TYPE E[N]; // напряженность
TYPE po[N+1]; // плотность
TYPE phi[N+1]; // потенциал
TYPE phi2[N+1];
TYPE q;
/*функции*/
int initialization();
TYPE nn(TYPE x); //считает плотность в точке
TYPE spl(TYPE x); // строит линейный сплайн по E[N]

void save_rho(char* filename);//сохраняет массивы в файл 
TYPE MaxPhi(){ //условие остановки в прогонке
	double max = 0;
	for (int i = 0; i < N+1; ++i)
	{
		if (max<fabs(phi[i]-phi2[i])) max = fabs(phi[i]-phi2[i]);
	}
	return max;
}
TYPE Nexp(TYPE x){
	return exp(x);
}

int main()
{
	int ms = 10;
	double T=0;
	int i=0;
	initialization();
	TYPE alpha[N]; //прогонка
	TYPE beta[N]; //прогонка
	TYPE s; //прогонка
	char filename[80]; //для вывода в файл
	sprintf(filename, "hybrid200000p00001dt1000N.txt");
	while(T<=0.10){
		/*Прогонка*/
		do
		{
			//printf("\t*");
			for (int k = N; k >= 0; k--)
			{
				phi2[k] = phi[k]; 
			}
			phi[N]=0;				//краевые условия
			alpha[0] = 0;			//краевые условия
			beta[0] = log(po1);			//краевые условия
			for (int k = 1; k < N; ++k)
			{
				s = 2*B/(dx*dx)+Nexp(phi[k])-B*alpha[k-1]/(dx*dx);
				alpha[k]=B/(dx*dx*s);
				beta[k] = ((B*beta[k-1])/(dx*dx)-Nexp(phi[k])*(1-phi[k])+po[k])/s;//po[k-1]
			}
			for (int k = N-1; k >= 0; k--)
			{
				phi[k] = alpha[k]*phi[k+1]+beta[k]; 
			}
		}
		while(MaxPhi()>eps);
	
		/*Остальные величины*/
		for (int i = 1; i <= N; ++i)
		{
			E[i-1] = (phi[i-1]-phi[i])/dx;
		}
		if(i==0){
			for (int i = 0; i < M; ++i)
		{
			//if (!(qua[i]>0 && qua[i]<1)) V[i]=0;
			V[i]=V[i]+spl(qua[i])*dt/2.;

		}}
		else{
			for (int i = 0; i < M; ++i)
		{
			//if (!(qua[i]>0 && qua[i]<1)) V[i]=0;
			V[i]=V[i]+spl(qua[i])*dt;

		}
		}
		

		for (int i = 0; i < M; ++i)
		{
			qua[i]=qua[i]+V[i]*dt;
			if( qua[i]>=b) {qua[i]-=b; V[i]=-V[i];}
			if( qua[i]<0) {qua[i]+=b;V[i]=-V[i];}
		}
		for (int i = 1; i < N; ++i)
		{
			po[i]=nn(dx*(i-1./2));
		}
		
		po[0]=po1;
		po[1]=po1;
		po[N-1]=1;
		po[N]=1;
		printf("%f \n", T);
		T+=dt;
		i++;
	}
	save_rho(filename);
	system("g.gpi");
	return 0;
}



TYPE spl(TYPE x){
	int j = floor(x/dx);
	if(j>=0 && j<N-1){
	TYPE s = x/dx- j;
	return (1-s)*E[j]+(s)*E[j+1];
	}
	return 0;
}
TYPE nn(TYPE x){
	//TYPE q = ((po1+po2)/2.)/M;

 	double t=0;
 	for (int i = 0; i < M; ++i)
 	{
 		if(fabs(x-qua[i])<=dx) t+=(1./dx)*(1-fabs(x-qua[i])/dx);
 	}
 	return q*t;
}
/*int initialization(){
	TYPE m1 = W1*(1-L);
	TYPE m2 = W2*L;
	TYPE m = m1+m2;
	q = m/M;
	TYPE n1 = m1/q;
	TYPE n2 = m2/q;
	TYPE h1 = (1-L)/n1;
	TYPE h2 = L/n2;
	qua[0] = 0;
	int part1 = (int)n1/2.;
	int i = 0;
	for (i = 1; i < part1; ++i)
	{
		qua[i] = qua[i-1]+h1;
	}
	for (i; i < part1+(int)n2; ++i)
	{
		qua[i] = qua[i-1]+h2;
	}
	for (i; i < M; ++i)
	{
		qua[i] = qua[i-1]+h1;
	}

	for (int i = 0; i < M; ++i)
	{
		V[i]=0;
	}
	for (int i = 0; i < N+1; ++i)
	{
		phi[i]=0;
	}
	for (int i = 1; i < N; ++i)
	{
		po[i]=nn(dx*(i-1./2));
	}
	po[0]=W1;
	po[1]=W1;
	po[N-1]=W1;
	po[N]=W1;
	return 0;
}*/
int initialization(){
	TYPE m1 = (po1*ab);
	TYPE m2 = (po2*(b-ab));
	TYPE m = m1+m2;
	q = m/M;
	TYPE n1 = m1/q;
	TYPE n2 = m2/q;
	TYPE h1 = ab/n1;
	TYPE h2 = (b-ab)/n2;
	qua[0] = 0;
	int i = 0;
	for (i = 1; i < n1; ++i)
	{
		qua[i] = qua[i-1]+h1;
	}
	for (i; i < M; ++i)
	{
		qua[i] = qua[i-1]+h2;
	}

	for (int i = 0; i < M; ++i)
	{
		V[i]=0;
	}
	for (int i = 0; i < N+1; ++i)
	{
		phi[i]=0;
	}
	for (int i = 1; i < N; ++i)
	{
		po[i]=nn(dx*(i-1./2));
	}
	po[0]=po1;
	po[1]=po1;
	po[N-1]=1;
	po[N]=1;
	return 0;
}
void save_rho(char* filename){
	FILE *file; 
	file = fopen(filename, "w");
	for(int i = 1; i < N; ++i){
		fprintf(file, "%f %f", (i-1/2.)*dx, po[i]); 
		fprintf(file, "\n");
	}
	fclose(file);
}
