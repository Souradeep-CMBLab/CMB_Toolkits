/* *********************************************************************************************************************************************************************
                                            cholesky for A_{ll'}^{20}  for unbinned A_{ll'}^{20} for l= 512

						No much change except the file name
					  *****************************

*/
#include<stdio.h>
#include<stdlib.h>
#include<math.h>

extern "C" {
   void drc3jj_(double *fl2, double *fl3, double *fm2, double *fm3, double *flmin, double *flmax, double *fcg, int *ndh, int *fie);
     }
     
int  main ()
{
	int lmax, lowm, lowmp,n,nprime,lval=0,lvaln=0,  lvalbi=0,nd2, ie2,ul3chk,po,pn,lst,lnd,lvalbi2, kchk, shftn, shftn1, shftn2, npmax, coch, lchk, lbin=2, lbin2=2, abm, elp,y,xev,xev1;
	double cel[7000],mmh2,llm2,llh2,lmin2,ell,elm,totj, totm,fcg2,cgel[5000],zcgel[5000], dia, lterm, ltemp, tlowertrm, lowertrm, ltemp1, ltemp2, ndtrm2, lowertrm2, lowerdg,evendiag[10000]={0}, odddiag[10000]={0},  evendiagn[10000]={0}, odddiagn[10000]={0};
	
	double lfcg2,elmin, elmax,zfcg2,zlfcg2,flowertrm2,tt, clf, den, biposh1, biposh2, gar1, gar2, gar3, gar4, gar5, kv, pile, pil;

	double biposh2p[7000], biposh,mbipo[7000], bden;
	FILE *inpbipo, *inpcl, *outlow, *inpbipoll2;
	
	inpcl = fopen("bare_Planck2015TTlowP_totCls_v2.dat", "r");

	outlow= fopen("./Cholesky_output/Planck2015TTlowP_scale_dep_dm_m10_0-14_lowermat_cl_L_1_lmax_1024.dat", "w");	

	inpbipoll2=fopen("./BipoSH_spectra/Planck2015TTlowP_totCls_alpha_10_l_lp1_scale_dep_dm_0-14.dat", "r");
	
	printf("enter the value of lmax");
	scanf("%d", &lmax);
	printf("enter the value of J");
	scanf("%lf", &totj);
	printf("enter the value of M");
	scanf("%lf", &totm);
	npmax= ((lmax+1)*(lmax+1))-1;

	for(int j=2; j<=lmax; j++)
        {

	fscanf(inpcl,"%d %le", &lvaln, &cel[j]);   //change made for Planck 
               
                fscanf(inpbipoll2,"%d %le", &lvalbi2, &biposh2p[j]); //change made for saurabh
//                      fscanf(inpbipoll2,"%d %le %le %le %le", &lvalbi2, &biposh2p, &gar1, &gar2, &gar3);

              printf("value of cl, biposh and biposh2p are %le\t %le\n", cel[j], biposh2p[j]);

	}

	for(int el=2; el<=lmax; el++)
	{
		lowm= (-1)*el;
		kchk= 1;
	
//		fscanf(inpcl,"%d %le %le %le %le %le", &lval, &cel, &gar1, &gar2, &gar3, &gar4);
	//	fscanf(inpcl,"%d %le", &lvaln, &cel);   //change made for Planck 

	//		lval= lvaln-1;
          //     fscanf(inpbipo,"%d. %le", &lvalbi, &mbipo); // changes made for Planck

//			fscanf(inpbipo,"%d %le %le %le %le", &lvalbi, &mbipo, &gar1, &gar2, &gar3);
	       
	    //   if(el<=lmax-2)
	      // {
               // fscanf(inpbipoll2,"%d. %le", &lvalbi2, &biposh2p); //change made for saurabh
//			fscanf(inpbipoll2,"%d %le %le %le %le", &lvalbi2, &biposh2p, &gar1, &gar2, &gar3);
	       //}
                //printf("value if cl, biposh and biposh2p are %le %le %le", cel, mbipo, biposh2p);
		mmh2= 0;
							llm2=totj-el;
							llh2= totj+el;
    if(llm2<0)
    {
    llm2=el-totj;
    }
    if(llm2>=mmh2)
    {
    lmin2= llm2;
    }
    else
    {
    lmin2= mmh2;
    }
    nd2=(llh2-lmin2)+1;
    ell= el;
    elm=0;
    drc3jj_(&totj,&ell,&elm,&elm,&elmin,&elmax,zcgel,&nd2,&ie2);
    // printf("the values of l3min and l3max are %le \t %le \t \n",ll3min, ll3max);
    for(int l4=0;l4<nd2;l4++)
    {
    ul3chk=l4+elmin;
    if(ul3chk==el)
    {
     zfcg2= zcgel[l4];
//      printf(" 1 %le \n", zfcg2);
    }
    if(ul3chk==(el-1) && (el-1)>=2)
    {
	    zlfcg2= zcgel[l4];
//	   printf(" 2 %d \t %le \n",ul3chk, zlfcg2);
	    }
    }
		for(int m= lowm; m<= el; m++)
		{
		abm= abs(m);
		if(abm <= (el-1) && (el-1) >= 2)
		{
		elp= el-1;
	      	lowmp= (-1)*elp;
		      
		      // printf("a\n");
			       nprime= el*(el+1)+m-3;
			       n= elp*(elp+1)+m-3;
			       pn=pow(-1,m);
				 mmh2= m+totm;
							llm2=totj-el;
							llh2= totj+el;
    if(mmh2<0)
    {
    mmh2=(-1)*mmh2;
    }
    if(llm2<0)
    {
    llm2=el-totj;
    }
    if(llm2>=mmh2)
    {
    lmin2= llm2;
    }
    else
    {
    lmin2= mmh2;
    }
    nd2=(llh2-lmin2)+1;
    ell= el;
    elm=m;
    drc3jj_(&totj,&ell,&totm,&elm,&elmin,&elmax,cgel,&nd2,&ie2);
    // printf("the values of l3min and l3max are %le \t %le \t \n",ll3min, ll3max);
    for(int l4=0;l4<nd2;l4++)
    {
    ul3chk=l4+elmin;
    if(ul3chk==el)
    {
     fcg2= cgel[l4];
   //  printf(" 3 %d\t %le \n", ul3chk, fcg2);
    }
    if(ul3chk==(el-1) && (el-1)>=2)
    {
	    lfcg2= cgel[l4];
//	    printf(" 4 %d\t  %le \n",ul3chk, lfcg2);
	    
	    }
    }

		 bden= (elp*(elp+1))*0.159;
//                 biposh= biposh2p[elp]/bden;//changes made for saurabh//
                 biposh= biposh2p[elp];

		pil= sqrt((2*elp)+1);
                pile=sqrt((2*el)+1);
                lowertrm2= (biposh*pn*zlfcg2*lfcg2*1.73205*pil*pile);
		xev1= elp%2;	
		if (xev1 == 0)
                {
                 y= abs(m);
                 lowerdg=evendiag[y];
	//	printf("y value %d and lowerdg %le\n", y, lowerdg);		
                flowertrm2= lowertrm2/lowerdg;
               fprintf(outlow,"%d\t %d\t %le \n", n , nprime, flowertrm2);
                }
                else if(xev1!=0)
                {
                y= abs(m);
                 lowerdg=odddiag[y];
	//	 printf("y value %d and lowerdg %le\n", y, lowerdg);
//		printf("lowerdg %le\n", lowerdg);
	        flowertrm2= lowertrm2/lowerdg;
               fprintf(outlow,"%d\t %d\t %le \n", n , nprime, flowertrm2);
                }
               // flowertrm2= lowertrm2/lowerdg;
              // fprintf(outlow,"%d\t %d\t %le \n", n , nprime, flowertrm2);
		den= (el*(el+1))*0.159;
                clf= cel[el];
//              clf= cel/den; //change for planck
//                bden= (el*(el+1))*0.159;
//                biposh= mbipo[el]/bden; // changes made for saurabh
//              biposh= mbipo[el];

              //printf("cl and biposh %le\t %le\n",clf, biposh);
               dia= clf;
               //  printf("diag %le\n", dia);

			
		 ltemp= (flowertrm2*flowertrm2);
                 kchk= kchk+1;
		tt= dia-ltemp;
                tlowertrm= fabs(tt);
                lowertrm= sqrt(tlowertrm);
               // fseek(outlow,1, SEEK_END);
               fprintf(outlow,"%d\t %d\t %le \n", nprime ,nprime,lowertrm);
		xev= el%2;
                if (xev == 0)
                {
                 y= abs(m);
                 evendiagn[y]= lowertrm;
                }
                else if(xev!=0)
                {
                y= abs(m);
                 odddiagn[y]= lowertrm;
                }


		}
		else
		{
		
		elp=el-2;
	      	lowmp= (-1)*elp;
		      
		      // printf("a\n");
			       nprime= el*(el+1)+m-3;
			       n= elp*(elp+1)+m-3;
			       mmh2= m+totm;
							llm2=totj-el;
							llh2= totj+el;
    if(mmh2<0)
    {
    mmh2=(-1)*mmh2;
    }
    if(llm2<0)
    {
    llm2=el-totj;
    }
    if(llm2>=mmh2)
    {
    lmin2= llm2;
    }
    else
    {
    lmin2= mmh2;
    }
    nd2=(llh2-lmin2)+1;
    ell= el;
    elm=m;
    drc3jj_(&totj,&ell,&totm,&elm,&elmin,&elmax,cgel,&nd2,&ie2);
    // printf("the values of l3min and l3max are %le \t %le \t \n",ll3min, ll3max);
    for(int l4=0;l4<nd2;l4++)
    {
    ul3chk=l4+elmin;
    if(ul3chk==el)
    {
     fcg2= cgel[l4];
   //  printf(" 3 %le \n", fcg2);
    }
    if(ul3chk==(el-1) && (el-1)>=2)
    {
	    lfcg2= cgel[l4];
	   // printf(" 4  %le \n", lfcg2);
	    
	    }
    }		 pn=pow(-1,m);
		den= (lval*(lval+1))*0.159;
                clf= cel[el];
//              clf= cel/den; //change for planck
//                bden= (el*(el+1))*0.159;
//                biposh= mbipo[el]/bden; // changes made for saurabh
//              biposh= mbipo[el];

              //printf("cl and biposh %le\t %le\n",clf, biposh);
               dia= clf;
               //  printf("diag %le\n", dia);


                 ltemp= 0;
                 kchk= kchk+1;
                tt= dia;
                tlowertrm= fabs(tt);
                lowertrm= sqrt(tlowertrm);
               // fseek(outlow,1, SEEK_END);
               fprintf(outlow,"%d\t %d\t %le \n", nprime ,nprime,lowertrm);
		 xev=el%2;
                if (xev == 0)
                {
                 y= abs(m);
                 evendiagn[y]= lowertrm;
                }
                else if(xev!=0)
                {
                y= abs(m);
                 odddiagn[y]= lowertrm;
                }

		/*if(nprime==1)
		{
		lowerdg= lowertrm;	
		printf("l11 is %le\n", lowerdg);
		}


	*/	}
/*		for(int elp=el; elp<=el+2; elp=elp+2)
	     { 
	        if(elp<= lmax)
		{
	     // kchk=1;
	     	
	      	lowmp= (-1)*elp;
		      
		      // printf("a\n");
			       n= el*(el+1)+m-3;
			       nprime= elp*(elp+1)+m-3;
			       mmh2= m+totm;
							llm2=totj-el;
							llh2= totj+el;
    if(mmh2<0)
    {
    mmh2=(-1)*mmh2;
    }
    if(llm2<0)
    {
    llm2=el-totj;
    }
    if(llm2>=mmh2)
    {
    lmin2= llm2;
    }
    else
    {
    lmin2= mmh2;
    }
    nd2=(llh2-lmin2)+1;
    ell= el;
    elm=m;
    drc3jj_(&totj,&ell,&totm,&elm,&elmin,&elmax,cgel,&nd2,&ie2);
    // printf("the values of l3min and l3max are %le \t %le \t \n",ll3min, ll3max);
    for(int l4=0;l4<nd2;l4++)
    {
    ul3chk=l4+elmin;
    if(ul3chk==el)
    {
     fcg2= cgel[l4];
   //  printf(" 3 %le \n", fcg2);
    }
    if(ul3chk==(el-2) && (el-2)<=lmax)
    {
	    lfcg2= cgel[l4];
	   // printf(" 4  %le \n", lfcg2);
	    
	    }
    }
			       
			       if(n==nprime)
			       {
				       
				       pn= pow(-1,m);
				       
				       
				       
				       
				       
				       
					       
					       
					       
					       den= (lval*(lval+1))*0.159;
					       clf= cel;
//					       clf= cel/den; //change for planck
					       bden= (el*(el+1))*0.159;
					       biposh= mbipo/bden; // changes made for saurabh
//						biposh= mbipo;

					       //printf("cl and biposh %le\t %le\n",clf, biposh);
					       dia= clf+((fcg2)*zfcg2*biposh*pn*((2*el)+1)*2.2360679);
					     //  printf("diag %le\n", dia);
					       ltemp=0;
					      // kchk=1;
					      if(elp>4)
					      {
					       shftn= elp*(elp+1)-(elp+2);
					       }
					       else
					      
					       {
					       shftn= 14;
					       }
					       lchk= fabs(m);
					       
					       
					       if((n-shftn)>0 && lchk<= (elp-2))
					       {
					      //rewind(outlow);
						 fseek(outlow, 0, SEEK_SET);

					        while(kchk <= (n-shftn))
						{
						//fseek(outlow, 1, SEEK_SET);
						fscanf(outlow,"%d %d %le", &lst, &lnd, &lterm);
						//printf( " out %d, %d, %le \n", lst, lnd, lterm);
						if(lst== (((elp-2)*(elp-1)+m-3)) && lnd==nprime)
						{
						ltemp= ltemp+(lterm*lterm);
						 kchk= kchk+1;
						 // printf( " out %d, %d, %le \n", lst, lnd, lterm);
						}
						//kchk= kchk+1;
						//shftn= shftn+3*(npmax);
						
										       
						       }
						       }
						       else
						       {
							       ltemp= 0;
							       }
							      
							      							      
							       
							     tt= dia-ltemp;  
							     tlowertrm= fabs(tt);
							     lowertrm= sqrt(tlowertrm);
							    // fseek(outlow,1, SEEK_END);
							     fprintf(outlow,"%d\t %d\t %le \n", n ,nprime,lowertrm);
							     
				      }
				      if(n!= nprime)
				      {
					      ltemp1= 0;
					      ltemp2= 0;
					      ndtrm2=0;
					      pn=pow(-1,m);
					     // kchk=1;
					       //shftn1= 3*(n-1);
					     //shftn2= nprime;
					      if(elp == (el+2))
					      {
					         
						     bden= (el*(el+1))*0.159; 
						     biposh= biposh2p/bden;//changes made for saurabh//
//							biposh= biposh2p;
						    
						     
						    
					//	     printf(" bipo %le \n", biposh);
						     
						      }
						      
						      
						      else
						      {
							      biposh=0;
							      }
						/*	      if((n-1)>0)
					       {

					        while(kchk <= n-1)
						{
						//fseek(outlow, shftn1, SEEK_SET);
						fscanf(outlow,"%d %d %le", &lst, &lnd, &lterm);
						//shftn1= shftn1+3*(npmax);
						if(lst==kchk && lnd==n)
						{
						kchk= kchk+1;
						ltemp1= (lterm);
						}
						else
						{
						ltemp1=0;
						}
						//fseek(outlow, shftn2, SEEK_SET);
						fscanf(outlow,"%d %d %le", &lst, &lnd, &lterm);
						//shftn2= shftn2+3*(npmax);
						if(lst==kchk && lnd==nprime)
						{
						kchk= kchk+1;
						ltemp2= lterm;
						}
						else
						{
						ltemp2= 0;
						}
						ndtrm2= ndtrm2+(ltemp2*ltemp1);
						//kchk= kchk+1;
			       								       
						       }
						       }
						       else
						       {
							       ndtrm2=0;
							       }
							       pil= sqrt((2*elp)+1);
						               pile=sqrt((2*el)+1);
							       lowertrm2= (biposh*pn*zlfcg2*lfcg2*2.2360679*pil*pile);
							       flowertrm2= lowertrm2/lowertrm;
							       
							      
							      fprintf(outlow,"%d\t %d\t %le \n", n , nprime, flowertrm2);
					      }
					      
			
	   	}
		
		/*else
		{
		 flowertrm2=0;

		 fprintf(outlow,"%d\t %d\t %le \n", n , nprime, flowertrm2);
		 }
		       }
		*/	}
			printf("value lo l %d\n", el);
		for(int ji= 0; ji<=el; ji++)
                {
                evendiag[ji]= evendiagn[ji];
                odddiag[ji]= odddiagn[ji];
                }
		}
	//	fclose(inpbipo);
		fclose(inpcl);
		fclose(inpbipoll2);
		fclose(outlow);

		printf("done");
		return 0;
}
	
