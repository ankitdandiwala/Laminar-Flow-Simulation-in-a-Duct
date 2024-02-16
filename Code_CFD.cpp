#include <stdio.h>
#include<conio.h>
#include <math.h>
main()
{
      int i,j,n,m;
      float T[101][101],Tn[101][101];
      float x[101],y[101],u[101],pex[101];
      float p[101],p1[101],p2[101],p3[101],p4[101];
      float B,L,dx,dy,beta,umax,gama;
      float diff,sq,rms,sum=0;
      float H;
      int count=0;
      FILE *f1;
      printf("Enter the length of the domain\n");
      scanf("%f",&L);
      printf("Enter the height of the domain\n");
      scanf("%f",&B);
      printf("Enter the number of subdomain in i-direction\n");
      scanf("%d",&m);
      printf("Enter the number of subdomain in j-direction\n");
      scanf("%d",&n);
      dx=L/m;
      dy=B/n;
      beta=dx/dy;
      printf("Enter the maximum velocity\n");
      scanf("%f",&umax);
      printf("Enter the value of gama\n");
      scanf("%f",&gama);
      //GRID GENERATION
      x[0]=0;
      y[0]=0;
      x[1]=dx/2;
      y[1]=dy/2;
      for (i=2;i<=m;i++)
      {
          x[i]=x[i-1]+dx;
      }
      for (j=2;j<=n;j++)
      {
          y[j]=y[j-1]+dy;
      }
      x[m+1]=L;
      y[n+1]=B;
      
      for(j=0;j<=n+1;j++)
      {
           u[j]=4*umax*((y[j]/B)-((y[j]/B)*(y[j]/B)));
           pex[j]=u[j]*dx/gama;
           p[j]=2*dy-pex[j];
           p1[j]=6*dy*beta*beta+6*dy-pex[j];
           p2[j]=6*dy*beta*beta+2*dy+pex[j];
           p3[j]=4*dy*beta*beta+2*dy+pex[j];
           p4[j]=4*dy*beta*beta+6*dy-pex[j];         
      }
      
      //BOUNDARY CONDITIONS
      printf("Enter the temperature of the North & South wall\n");
      scanf("%f",&H);
      for(i=0;i<=n+1;i++)
      {
          T[i][n+1]=H;
          Tn[i][n+1]=H;
          T[i][0]=H;
          Tn[i][0]=H;
      }
      //INITIALIZATION
      for(i=0;i<=m+1;i++)
      {
          for(j=1;j<=n;j++)
          {
              T[i][j]=0;
              Tn[i][j]=0;
          }
      }
      start:
      printf ("%d\t%f\n",count,rms);
      sum=0;
      //Equations
      //Left bottom cell
      Tn[1][1]=(p[1]/p1[1])*T[2][1]-(2*dy*beta*beta/p1[1])*(T[1][2]+2*H);
      //Right bottom cell
      Tn[m][1]=(p[1]/p2[1])*T[m-1][1]-(2*dy*beta*beta/p2[1])*(T[m][2]+2*H);
      //Top left cell
      Tn[1][n]=(p[n]/p1[n])*T[2][n]-(2*dy*beta*beta/p1[n])*(T[1][n-1]+2*H);
      //Top right cell
      Tn[m][n]=(p[n]/p2[n])*T[m-1][n]-(2*dy*beta*beta/p2[n])*(T[m][n-1]+2*H);
      //The south row of cells
      for(i=2;i<=m-1;i++)
      {
           Tn[i][1]=(p[1]/(6*dy*beta*beta+4*dy))*(T[i-1][1]+T[i+1][1])-(beta*beta/(2+3*beta*beta))*(T[i][2]+2*H);
      }
      //The north row of cells
      for (i=2;i<=m-1;i++)
      {
           Tn[i][n]=(p[1]/(6*dy*beta*beta+4*dy))*(T[i-1][n]+T[i+1][n])-(beta*beta/(2+3*beta*beta))*(T[i][n-1]+2*H);
      }
      //The west row of cells
      for (j=2;j<=n-1;j++)
      {
           Tn[1][j]=(p[j]/p4[j])*T[2][j]-(2*dy*beta*beta/p4[j])*(T[1][j-1]+T[1][j+1]);
      }
      //The east row of cells
      for (j=2;j<=n-1;j++)
      {
           Tn[m][j]=(p[j]/p3[j])*T[m-1][j]-(2*dy*beta*beta/p3[j])*(T[m][j-1]+T[m][j+1]);
      }
      //Interior cells
      for (i=2;i<=m-1;i++)
      {
           for (j=2;j<=n-1;j++)
           {
                Tn[i][j]=(p[j]/(4*dy+4*dy*beta*beta))*(T[i+1][j]+T[i-1][j])-(beta*beta/(2+2*beta*beta))*(T[i][j+1]+T[i][j-1]);
           }
      }
      //Outlet
      for (j=2;j<=n-1;j++)
      {
           Tn[m+1][j]=Tn[m][j];
      }
      //Calculation of rms
      for(i=1;i<=m;i++)
      {
           for(j=1;j<=n;j++)
           {
               diff=Tn[i][j]-T[i][j];
               sq=diff*diff;
               sum=sum+sq;
           }
      }
      rms=sqrt(sum/n*m);
      if(rms<0.0001)
      {
          printf("rms=%f\n",rms);
          for(i=0;i<=m+1;i++)
          {
              for(j=0;j<=n+1;j++)
              {
                  printf("%d\t%d\t%f\n",i,j,Tn[i][j]);
              }
          }
          f1=fopen("output.dat","w");
          for(i=0;i<=m+1;i++)
          {
              for(j=0;j<=n+1;j++)
              {
                  fprintf(f1,"%f\t%f\t%f\n",x[i],y[j],Tn[i][j]);
              }
          }
          printf("Solution is converged\n");
      }
      else
      {                  
          for (i=1;i<=m;i++)
          {
              for (j=1;j<=n;j++)
              {
                  T[i][j]=Tn[i][j];
                  
              }
          }
	   count=count+1;      
       goto start;
       }
       getch();
}
