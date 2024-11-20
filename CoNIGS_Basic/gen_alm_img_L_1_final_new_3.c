/*  ********************************************************************************************************************************************

                                        code for REAL PART of generating map from   the  lower triangular matrix from cholesky for non zero off-diagonal terms. 

						for all values of l.  Array size is NOT a matter
					****************************************************************************


*/
#include<stdio.h>
#include<stdlib.h>
#include<math.h>

int main()
{

int lmax,n,np, chk, col, row, bcol, brow,bchk,aj,as,nb,nbp,n1,hmax,hmin,nchk,ic=12, achk;
long int sz=0,an,tp,ne;
double ob[5]={0}, lval,bval, aval,sum, ranary[100000000];
FILE *inpbip, *inpgrv, *outmat;

char fil[100], outfil[100];

//inpbip = fopen("./Cholesky_output/Planck2015TTlowP_dm_m10_0-14_lowermat_cl_L_1_lmax_1024.dat", "r");
inpbip = fopen("./Cholesky_output/Planck2015TTlowP_scale_dep_dm_m10_0-14_lowermat_cl_L_1_lmax_1024.dat", "r");

printf("enter the value of lmax");
scanf("%d", &lmax);

printf("enter the no of maps");
scanf("%d", &hmax);

printf("enter the value of hmin");
scanf("%d", &hmin);

n1= lmax*(lmax+1)+lmax-3;

for(int h=hmin; h<=hmax; h++)
{
 
 sprintf(fil,"./Random_numbers/imgarray4096_%d.dat", h);
 inpgrv= fopen(fil, "r");

 for(int mj= 1; mj<= n1; mj++)
{
	fscanf(inpgrv,"%le", &bval);
	ranary[mj]= bval;
}
 fseek(inpbip,0,SEEK_SET);

sprintf(outfil, "./Gen_alm_output/alm_Planck2015TTlowP_scale_dep_dm_m10_0-14_L_1_img_1024_%d.dat",h);

outmat=fopen(outfil,"w");
for(int v=1; v<=3;v++)
{
fprintf(outmat,"%le\n", ob[v]);
}
//for(int g=1; g<=n1; g++)
//{
//fscanf(inpgrv,"%le", &bval[g]);
//}
ic=12;
for(int i=2; i<=lmax; i++)
{
if(ic==i)
{
printf("i is %d\n", i);
ic= ic+50;
}
if((i-1)>=2)
{
 aj= i-1;
  }
   if((i-1)<2)
    {
     aj= i;
      }
/*achk=0;
fseek(inpgrv,0,SEEK_SET);
//for(int kj=aj; kj<=i; kj=kj+2)
if(i<=lmax)
{
an= aj*(aj+1)-aj -3;
while(achk<an)
{
fscanf(inpgrv,"%le", &bval);
achk= achk+1;
ranary[0]= bval;
}
ne= (6*aj)+9;
for(int kj=1; kj<ne; kj=kj+1)
{

fscanf(inpgrv,"%le", &bval);
ranary[kj]= bval;
//fscanf(inpgrv,"%le", &bval);
}
}
*/
for(int l=(-1)*i; l<=0; l++)
{
//for(int k=1;k<=n1;k++)
//{

///aval[k]=0;
//}

//if(abs(l)<=(i-2))
//{
//fseek(inpbip,sz,SEEK_SET);
//}
//printf("sz value is %ld\n", sz);
sum=0;
for(int j=aj; j<=i; j=j+1)
{

if(abs(l)<=j)
{

n= i*(i+1)+l-3;
np= j*(j+1)+l-3;
chk=0;
//printf("i is %d\n ", i);
//fseek(inpbip,0,SEEK_SET);
//fseek(inpgrv,0,SEEK_SET);
//nchk=0;
/*if(j==(i-2))
{
fseek(inpbip,sz,SEEK_SET);
//sz= ftell(inpbip);
//printf("sz value is %ld\n", sz);
}
*/
while(chk<1)
{
fscanf(inpbip,"%d %d %le", &col, &row, &lval);
/*if(j==(i-2))
{
fseek(inpbip,sz,SEEK_SET);
sz= ftell(inpbip);
//printf("sz value is %ld\n", sz);
}
*/
//fscanf(inpbip,"%d", &row);
//fscanf(inpbip,"%le", &lval);
//printf("value of col and lval %d\t %le\n ", col, lval);

if(row==n && col==np)
{
aval= lval;
//printf("value of col and lval %d\t %le\n ", col, lval);
chk=chk+1;
/*if(j==(i-2))
{
sz= ftell(inpbip);
}
*//*if(j==(i-2))
{
fseek(inpbip,sz,SEEK_SET);
sz= ftell(inpbip);
//printf("sz value is %ld\n", sz);
}
*/
}
}
//if(j==aj)
//{
//sz= ftell(inpbip);
//}
/*while(nchk<col)
{
fscanf(inpgrv,"%le", &bval);
nchk= nchk+1;
}
*/
//tp= col-an;
sum= sum+(ranary[col]*aval*0.707);


}

}
ob[4]= sum;
/*for(int s=1;s<=n1; s++)
{
ob[4]=ob[4]+ (bval[s]* aval[s]*0.707); // correction made
}
*/
fprintf(outmat,"%le\n", ob[4]);
//printf("%d\t  %le\n", n, ob[4]);
ob[4]=0;

}
}
fclose(inpgrv);

fclose(outmat);
}
printf("done");
//fclose(inpgrv);
//fclose(outmat);
fclose(inpbip);
}

