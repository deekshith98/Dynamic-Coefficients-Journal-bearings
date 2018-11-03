import java.lang.Math.*;
import java.io.*;
public class coefficientstwo
{	
  static  int size=100;//size of the grid
    static double dphi=6.6/(1.4*size);
 double  cl=0.000075,r=0.05,l=0.1,phi=0.0,dz=(l/r)/size,dphir=6.6/(1.4*size),ec=0.5,x=dphir/dz,dphiktemp,phitemp1,phitemp2,c8=15852.476,w;
//cl-clearance, r-bearing radius, l-length of the bearing, dphir-element size in phi direction in radains (3PI/(2*size), ec-eccentricity ratio, w is angular speed of rotor
 double deltaC=0.00000001,deltaPhi=0.0005,dphik,alpha=1.8,forceX,forceY,pmax=0;
 double kxx,kxy,kyx,kyy,cxx,cxy,cyx,cyy,som;
 int counter=0;
 int[] bd=new int[size];//to record boundary nodes of oil film after calculating static pressure
//to store coefficients of the equations involving pressure  
double[][] a = new double[size][size]; 
double[][] b = new double[size][size]; 
double[][] c = new double[size][size]; 
double[][] d = new double[size][size];
double[][] e = new double[size][size]; 
double[][] f = new double[size][size];  
double[][] sP = new double[size][size]; 
double[][] sQ = new double[size][size]; 
double[][] sO = new double[size][size];
//to store pressure values. In java default values are zero so no need to initialize pressure of each element to zero 
double[][] p = new double[1+size][1+size]; 
double[][] px = new double[1+size][1+size]; 
double[][] py = new double[1+size][1+size]; 
double[][] pX = new double[1+size][1+size]; 
double[][] pY = new double[1+size][1+size];
public static void main(String[] args)
{
coefficientstwo test = new coefficientstwo();
test.eqCoeffCalc(dphi);// this method is called first to calculate equations' coefficients
test.pCalc();//this method calculates static pressure
// The above two function will give attitude and angle and load capacity of bearing for the given eccentricity
System.out.println(“The load capacity for the eccentricity ”+ec+” is ”+forceY);
test.extraCoeffCalc();// calculates other coefficients which appear in f[i,j] to solve for px py px. and py.
test.coeffCalc();//to calculate dynamic coefficients 
test.print();// prints the final values 
}
public void eqCoeffCalc(double k)
{ 
 int m,n;
  for(m=1;m<size;m++)
{
 for(n=1;n<size;n++)
     { 
 a[m][n]=Math.pow((1+ec*Math.cos((m+1.0/2.0)*k)),3);
       b[m][n]=Math.pow((1+ec*Math.cos((m-1.0/2.0)*k)),3);
       c[m][n]=x*x*Math.pow((1+ec*Math.cos(m*k)),3);
       d[m][n]=x*x*Math.pow((1+ec*Math.cos(m*k)),3);
e[m][n]=a[m][n]+b[m][n]+c[m][n]+d[m][n];
f[m][n]=dphir*(ec*Math.cos((m+1.0/2.0)*k)-ec*Math.cos((m-1.0/2.0)*k));
 }}}
public void extraCoeffCalc()
{
  int m,n;
 for(n=1;n<size;n++)
  {
    for(m=1;m<bd[n];m++)
    {
sP[m][n]=Math.pow((1+ec*Math.cos(dphi*m)),2)*p[m][n-1]+Math.pow((1+ec*Math.cos(dphi*m)),2)*p[m][n+1]-2*Math.pow((1+ec*Math.cos(dphi*m)),2)*p[m][n];
sQ[m][n]=Math.pow((1+ec*Math.cos(dphi*(m+1.0/2.0))),2)*Math.cos(phi+dphi*(m+1.0/2.0))*p[m+1][n]+Math.pow((1+ec*Math.cos(dphi*(m-1.0/2.0))),2)*Math.cos(phi+dphi*(m-1.0/2.0))*p[m-1][n]-(Math.pow((1+ec*Math.cos(dphi*(m+1.0/2.0))),2)*Math.cos(phi+dphi*(m+1.0/2.0))+Math.pow((1+ec*Math.cos(dphi*(m-1.0/2.0))),2)*Math.cos(phi+dphi*(m-1.0/2.0)))*p[m][n];
sO[m][n]=Math.pow((1+ec*Math.cos(dphi*(m+1.0/2.0))),2)*Math.sin(phi+dphi*(m+1.0/2.0))*p[m+1][n]+Math.pow((1+ec*Math.cos(dphi*(m-1.0/2.0))),2)*Math.sin(phi+dphi*(m-1.0/2.0))*p[m-1][n]-(Math.pow((1+ec*Math.cos(dphi*(m+1.0/2.0))),2)*Math.sin(phi+dphi*(m+1.0/2.0))+Math.pow((1+ec*Math.cos(dphi*(m-1.0/2.0))),2)*Math.sin(phi+dphi*(m-1.0/2.0)))*p[m][n];
       }
}
}
public void pCalc()
{ 
int i,j;
double sum1,sum2; 
double[][] temp= new double[size][size];// to store pressure values of previous iteration to check for convergence
do// this loop runs until error function <deltaPhi
{
if(counter!=0)//first error function is calculated for two guess phi values .
{phitemp2=phi;
  phi=phi-dphik/((dphik-dphiktemp)/(phi-phitemp1));//calculating phi value using newton raphson formula
  phitemp1=phitemp2;
}
else
{
  do//this loop is to calculate error function for phi=0. this loop is executed only once
{
for(i=1;i<size;i++)
{
  for(j=1;j<size;j++)
  {
if(((phi+i*dphi)>1.39626&&(phi+i*dphi)<1.745329)||((phi+i*dphi)>4.537856&&(phi+i*dphi)<4.886921))// to check if the current angle lies in the groove region, 1.3962 is 80 degree in radains, 1.745329 is 100 degress in radians, 4.53 is 260 degree in radians, 4.886 is 280 degree in radians
      {
          p[i][j]=0;
 }
    else
     {temp[i][j]=p[i][j];
 p[i][j]=(1-alpha)*p[i][j]+alpha*(a[i][j]*p[i+1][j]+b[i][j]*p[i-1][j]+c[i][j]*p[i][j+1]+d[i][j]*p[i][j-1]-f[i][j])/e[i][j];
      if(p[i][j]<0)
        {
            p[i][j]=0;
        }
}
}
}
sum1=0;
sum2=0;
 for(i=1;i<size;i++)//loop to calculate the convergence value
 {
  for(j=1;j<size;j++)
  {if(((phi+i*dphi)>1.39626&&(phi+i*dphi)<1.745329)||((phi+i*dphi)>4.537856&&(phi+i*dphi)<4.886921))
      {}
      else
    {sum1=sum1+Math.abs(p[i][j]-temp[i][j]);
    sum2=sum2+Math.abs(p[i][j]);}
  }
 }   
}while((sum1/sum2)>deltaC);//condition to check for convergence
forceCalc(phi);//calls force calc method and calculates force values
phitemp1=phi;
dphiktemp=Math.atan(forceX/forceY);// error function
phi=11.0/126.0;// phi is set to 5 degrees
}
do// this loop is executed multiple times.this loop runs until roots of the equation converge for a given phi value.this loop implements over relaxation method 
{
for(i=1;i<size;i++)
{
  for(j=1;j<size;j++)
  {
if(((phi+i*dphi)>1.39626&&(phi+i*dphi)<1.745329)||((phi+i*dphi)>4.537856&&(phi+i*dphi)<4.886921))//loop to check if the current angle lies in the groove region
       {
          p[i][j]=0;
 }
    else
     {temp[i][j]=p[i][j];
p[i][j]=(1-alpha)*p[i][j]+alpha*(a[i][j]*p[i+1][j]+b[i][j]*p[i-1][j]+c[i][j]*p[i][j+1]+d[i][j]*p[i][j-1]-f[i][j])/e[i][j];
      if(p[i][j]<0)
        {
            p[i][j]=0;
        }
    }
 }
}
sum1=0;
sum2=0;
 for(i=1;i<size;i++)
 {
  for(j=1;j<size;j++)
  {
if(((phi+i*dphi)>1.39626&&(phi+i*dphi)<1.745329)||((phi+i*dphi)>4.537856&&(phi+i*dphi)<4.886921))
      {}
      else
    {sum1=sum1+Math.abs(p[i][j]-temp[i][j]);
    sum2=sum2+Math.abs(p[i][j]);}
  }
 }   
}while((sum1/sum2)>deltaC);
forceCalc(phi);
if(counter!=0)
{
    dphiktemp=dphik;
}
dphik=Math.atan(forceX/forceY);
counter++;
}while((Math.abs(dphik))>deltaPhi);//to check if the current phi value is the right attitude angle

for(j=1;j<size;j++)// loop to record oil film boundaries
{
    for(i=1;i<size;i++)
    {
if(((phi+i*dphi)>1.39626&&(phi+i*dphi)<1.745329)||((phi+i*dphi)>4.537856&&(phi+i*dphi)<4.886921))
      {}
      else
      {
        if(p[i][j]==0){
           bd[j]=i;
           break;
        }
        }
    }
}
for(i=1;i<size;i++)
{
    for(j=1;j<size;j++)
    {
    if(p[i][j]>pmax)
    {
       pmax=p[i][j];
    }
}
}
}
public void forceCalc(double k)
{
int i,j;
forceX=0;
forceY=0;
for(i=1;i<size;i++)//loop to integrate the force values
{
  for(j=1;j<size;j++)
  {
    forceX=forceX-c8*p[i][j]*dz*dphir*Math.sin(k+dphi*i);
    forceY=forceY-c8*p[i][j]*dz*dphir*Math.cos(k+dphi*i);
    }
}
}
public void coeffCalc()
{
     int i,j;
  pxCalc();//to calculate px
  pyCalc();//py
  pXCalc();//px.
  pYCalc();//py.
   for(j=1;j<size;j++)//integrating perturbed pressure to find dynamic coefficients
  {
      for(i=1;i<bd[j];i++)
     { 
         kxx=kxx-(c8/forceY)*dphir*dz*px[i][j]*Math.sin(phi+dphi*i);
         kxy=kxy-(c8/forceY)*dphir*dz*py[i][j]*Math.sin(phi+dphi*i);
         kyx=kyx-(c8/forceY)*dphir*dz*px[i][j]*Math.cos(phi+dphi*i);
         kyy=kyy-(c8/forceY)*dphir*dz*py[i][j]*Math.cos(phi+dphi*i);
         cxx=cxx-(c8/forceY)*dphir*dz*pX[i][j]*Math.sin(phi+dphi*i);
         cxy=cxy-(c8/forceY)*dphir*dz*pY[i][j]*Math.sin(phi+dphi*i);
         cyx=cyx-(c8/forceY)*dphir*dz*pX[i][j]*Math.cos(phi+dphi*i);
         cyy=cyy-(c8/forceY)*dphir*dz*pY[i][j]*Math.cos(phi+dphi*i);
}
    }
}
public void pxCalc()//px and other perturbed pressure are calculated similar to static pressure
{
  int i,j;
  double[][] temp= new double[size][size];
  double sum1,sum2;
  for(j=1;j<size;j++)// f[i,j ] is calculated for px. f is different for different perturbed pressures
  {
    for(i=1;i<bd[j];i++)
    {
      f[i][j]=dphir*dphir*Math.cos(phi+i*dphi)-3*sO[i][j]-3*x*x*Math.sin(phi+i*dphi)*sP[i][j];
    }
}
   do//loop implements over relaxation method
{
 for(j=1;j<size;j++)//calculates px only where fluid is present so i<bd[j] is used instead of i<size
  {
    for(i=1;i<bd[j];i++)
    {
if(((phi+i*dphi)>1.39626&&(phi+i*dphi)<1.745329)||((phi+i*dphi)>4.537856&&(phi+i*dphi)<4.886921))
         {}
         else{
      temp[i][j]=px[i][j];
px[i][j]=(1-alpha)*px[i][j]+alpha*(a[i][j]*px[i+1][j]+b[i][j]*px[i-1][j]+c[i][j]*px[i][j+1]+d[i][j]*px[i][j-1]-f[i][j])/e[i][j];
 }
  }
}
sum1=0;
sum2=0;
  for(j=1;j<size;j++)
  {
    for(i=1;i<bd[j];i++)
    {
if(((phi+i*dphi)>1.39626&&(phi+i*dphi)<1.745329)||((phi+i*dphi)>4.537856&&(phi+i*dphi)<4.886921)){}
        else{
    sum1=sum1+Math.abs(px[i][j]-temp[i][j]);
    sum2=sum2+Math.abs(px[i][j]);}
  }
 }   }
while((sum1/sum2)>deltaC);// checks for convergence of solution of px
}
public void pyCalc()
{
  int i,j;
 double[][] temp= new double[size][size];
  double sum1,sum2;
  for(j=1;j<size;j++)
  {
    for(i=1;i<bd[j];i++)
    {
      f[i][j]=-dphir*dphir*Math.sin(phi+i*dphi)-3*sQ[i][j]-3*x*x*Math.cos(phi+i*dphi)*sP[i][j];
    }
}
  do
{
 for(j=1;j<size;j++)
  {
    for(i=1;i<bd[j];i++)
    {
if(((phi+i*dphi)>1.39626&&(phi+i*dphi)<1.745329)||((phi+i*dphi)>4.537856&&(phi+i*dphi)<4.886921)){}
         else{
      temp[i][j]=py[i][j]; py[i][j]=(1-alpha)*py[i][j]+alpha*(a[i][j]*py[i+1][j]+b[i][j]*py[i-1][j]+c[i][j]*py[i][j+1]+d[i][j]*py[i][j-1]-f[i][j])/e[i][j];
}
  }
}
sum1=0;
sum2=0;
 for(j=1;j<size;j++)
  {
    for(i=1;i<bd[j];i++)
    {
if(((phi+i*dphi)>1.39626&&(phi+i*dphi)<1.745329)||((phi+i*dphi)>4.537856&&(phi+i*dphi)<4.886921)){}
        else{
    sum1=sum1+Math.abs(py[i][j]-temp[i][j]);
    sum2=sum2+Math.abs(py[i][j]);}
  }
 }   }
while((sum1/sum2)>deltaC);
}
public void pXCalc()
{
  int i,j;
double[][] temp= new double[size][size];
  double sum1,sum2;
    for(j=1;j<size;j++)
  {
    for(i=1;i<bd[j];i++)
    {
      f[i][j]=2*dphir*dphir*Math.sin(phi+i*dphi);
    }
 }
 do
{
for(j=1;j<size;j++)
  {
    for(i=1;i<bd[j];i++)
    {
if(((phi+i*dphi)>1.39626&&(phi+i*dphi)<1.745329)||((phi+i*dphi)>4.537856&&(phi+i*dphi)<4.886921)){}
        else{
      temp[i][j]=pX[i][j];
 pX[i][j]=(1-alpha)*pX[i][j]+alpha*(a[i][j]*pX[i+1][j]+b[i][j]*pX[i-1][j]+c[i][j]*pX[i][j+1]+d[i][j]*pX[i][j-1]-f[i][j])/e[i][j];
}
  }
}
sum1=0;
sum2=0;
   for(j=1;j<size;j++)
  {
    for(i=1;i<bd[j];i++)
    {
    sum1=sum1+Math.abs(pX[i][j]-temp[i][j]);
    sum2=sum2+Math.abs(pX[i][j]);
  }
 }   }
while((sum1/sum2)>deltaC);
}
public void pYCalc()
{
  int i,j;
 double[][] temp= new double[size][size];
  double sum1,sum2;
 for(j=1;j<size;j++)
  {
    for(i=1;i<bd[j];i++)
    {
 f[i][j]=2*dphir*dphir*Math.cos(phi+i*dphi);
    }
  }
do
{
 for(j=1;j<size;j++)
  {
    for(i=1;i<bd[j];i++)
    {
if(((phi+i*dphi)>1.39626&&(phi+i*dphi)<1.745329)||((phi+i*dphi)>4.537856&&(phi+i*dphi)<4.886921)){}
        else{
      temp[i][j]=pY[i][j];
 pY[i][j]=(1-alpha)*pY[i][j]+alpha*(a[i][j]*pY[i+1][j]+b[i][j]*pY[i-1][j]+c[i][j]*pY[i][j+1]+d[i][j]*pY[i][j-1]-f[i][j])/e[i][j];
   }
  }}
sum1=0;
sum2=0;
for(j=1;j<size;j++)
  {
    for(i=1;i<bd[j];i++)
    {
    sum1=sum1+Math.abs(pY[i][j]-temp[i][j]);
    sum2=sum2+Math.abs(pY[i][j]);
  } }   
}
while((sum1/sum2)>deltaC);
}
public void print()
{
    System.out.println(Math.toDegrees(phi));
   System.out.println("kxx="+(kxx*forceY/cl));
   System.out.println("kxy="+(kxy*forceY/cl));
   System.out.println("kyx="+(kyx*forceY/cl));
   System.out.println("kyy="+(kyy*forceY/cl));
   System.out.println("cxx="+(cxx*forceY/(cl*w)));
   System.out.println("cxy="+(cxy*forceY/(cl*w)));
   System.out.println("cyx="+(cyx*forceY/(cl*w)));
   System.out.println("cyy="+(cyy*forceY/(cl*w)));
}
}
