// Neil Howarth project 4/17/14

#include <stdio.h>
#include <math.h>
#include <stdlib.h>


double mypoly (double x, double poly[])
{
	return  poly[0]*x*x*x*x + poly[1]*x*x*x +poly[2]*x*x + poly[3]*x + poly[4];
}

void check_roots(double left, double right, double poly[])
{
	double f_left, f_right;

	f_left = mypoly(left,poly);
	f_right = mypoly(right,poly);
	if (fabs(f_left) < .1e-04)
		printf("Root detected at %.3f  \n",left);
	else
		if(fabs(f_right) < .1e-04)
		;
		else
			if(f_left*f_right <0)
				printf("Root detected at %.3f  \n", (left+right)/2);
	return;
}

bool isodd(int x)
{
	if (x % 2 == 1)
		return true;
	else
		return false;
}


double TrapIntegrate(double a, double b, int N, double poly[])
{
	double dx,x,sum=0,area;
	dx=(b-a)/N;
	x=a+dx;
	for(int k=1; k<N; k++)
	{
	sum = sum + mypoly(x,poly);
	x=x+dx;
	}
	area=.5*dx*(mypoly(a,poly)+2*sum+mypoly(b,poly));
	return area;

}

double Simpsonsrule(double a, double b, int N, double poly[])
{

	double h,x,area=0,f1,f2;
	h=(b-a)/N;
	x=a+h;
	f1=mypoly(a,poly);
	f2=mypoly(b,poly);
	area = (f1+f2);

	for(int k=1;k<N;k++)
	{
		if (isodd(k))
		{
		f1=mypoly(x,poly);
		area=area +4.0*mypoly(x,poly);
		}
		else
		{
		f2=mypoly(x,poly);
		area=area+2.0*mypoly(x,poly);
		}
		x=x+h;
	}
	area=h/3.0*area;
	return area;


}

void main()
{
	//initializations
	FILE * ip;
	ip = fopen("projectdata.txt","w");

	int N, k;
	double a0, a1, a2, a3, a4, a, b, step, left, right, fx;

	
	printf("Enter coefficients a4, a3, a2, a1, a0:\n");
	scanf("%lf,%lf,%lf,%lf,%lf",&a4,&a3,&a2,&a1,&a0);
	printf("Enter interval limits a, b (a<b):  \n");
	scanf("%lf,  %lf", &a, &b);
	printf("Enter step size:  \n");
	scanf("%lf",&step);
	printf("\n\n");
	

	/*
	a4 = 0;
	a3 = 1;
	a2 = -2.125;
	a1 = -25;
	a0 = 53.125;
	a = -6;
	b = 6;
	step = .01;
	*/

	double poly[] = {a4, a3, a2, a1, a0};
	N = ceil((b-a)/step);
	
	// check subintervals for roots

	for (k = 0; k<=N-1;k++)
	{
		left = a +k*step;
		if(k==N-1)
			right = b;
		else
			right = left + step;
		check_roots(left,right,poly);
	}
	check_roots(b,b,poly);

	for(double j=a;j<=b;j=j+step)
		{
			fx = mypoly(j,poly);
			fprintf(ip,"%lf\t%lf\n",j,fx);
		}

	//matlab

	fclose(ip);
	FILE *ML;
	ML=fopen("launchML.m","w");

	fprintf(ML,"load projectdata.txt; \n");
	fprintf(ML,"x=projectdata(:,1);\n");
	fprintf(ML,"y=projectdata(:,2);\n");
	fprintf(ML,"plot(x,y)\n");
	fprintf(ML,"title('Project'); \n");
	fprintf(ML,"grid on; \n");
	fprintf(ML,"axis on; \n");

	fclose(ML);

	system("\"C:\\Program Files\\MATLAB\\R2013a\\bin\\matlab.exe\" /r launchML");

	//area
	double l,r,p;

	printf("Enter the interval for the area: ");
	scanf("%lf,  %lf", &l, &r);

	printf("Enter number of iterations for the area: ");
	scanf("%lf,  %lf", &p);

	printf("\n\nThe area between %.1lf and %.1lf are: ",l,r);

	double tarea = TrapIntegrate(l,r,p,poly);
	printf("\n\n\tArea from Trapezoidal rule = %lf\n\n",tarea);

	double sarea = Simpsonsrule(l,r,p,poly);
	printf("\tArea from Simpsons rule = %lf\n\n",sarea);

	system("pause");	
}
