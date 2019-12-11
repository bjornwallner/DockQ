#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include "molecule.h"

//
// Any comments and suggestion may be sent to:
// Author: Björn Wallner
// E-mail: bjorn.wallner@liu.se
//   
//



int main(int argc,char *argv[])		/* Main routine */
{
  molecule      m[1];
  molecule      model[1];
  int           i,j,k;
  double        d;
  double        mean=0;
  int           atoms=0;
  int           residues=0,residues2=0;
  //  int           contacts[MAXRES][MAXRES]={{}};
  /*  int **contacts,**contacts2;
      double **dist,**dist2;*/
  
  //  int           res_contacts[20][20]={{}};
  //  int           restype[20]={};
  //  int           tot_res_contacts=0;
  /// int           type_i,type_j;
  int           current_res_i=0;
  int           current_res_j=0;
  long           selected_resnum;
  int           first=0;
  double        cutoff=0;
  FILE          *fp;
  int error=0,error2=0;
  double temp;
  int tmp=0;
  int binary=1;
  char chain_name; 
  int nativeA[5000];
  int nativeB[5000];
  int modelA[5000];
  int modelB[5000];
  char chainA,chainB,chainA2,chainB2;
  int interface_contacts=0;
  int interface_contacts_model=0;
  int verbose=0;
  char set='a';
  char set_switch[100];
  char arg[100];
  //float        sum=0;
  /* Parse command line for PDB filename */
  
  if(argc>=4)
    {
      strcpy(model[0].filename,argv[1]);
      strcpy(m[0].filename,argv[2]);
      cutoff=strtod(argv[3],NULL); 
      cutoff=cutoff*cutoff;
      for(i=4; i<argc; i++) {
	if(strcmp(argv[i],"-cb") == 0) 
	  set='c'; //Read in CB only
	if(strcmp(argv[i],"-all") == 0) 
	  set='a'; //Read all atoms except hydrogen (this is the current default)
	if(strcmp(argv[i],"-verbose") == 0) 
	  verbose=1;
	
      }
      
      printf("%c\n",set);
  
  //if(argc>5)
      //	verbose=1;
    }
  else
    {
      printf("Usage: fnat [pdb_file] [native] cutoff [set default all-atom] [verbose (default off)]]\n");
      exit(1);
    }


  error=read_molecules(m,set);
 
  if(error==0)
    {  
      residues=m[0].residues;
      /*
      contacts = malloc(residues*sizeof(int *));
      dist = malloc(residues*sizeof(double *));
      for(i=0;i<residues;i++) {
	contacts[i] = malloc(residues*sizeof(int));
	dist[i] = malloc(residues*sizeof(double));
	}
      for(i=0;i<residues2;i++) {
	for(i=0;i<residues2;i++) {
	  contacts[i][j]=0;
	  dist[i][j]=999999;
	}
      }
      */
      current_res_i=0;
      for(i=0;i<m[0].atoms;i++)  
	{
	  if(m[0].atm[i].rescount!=current_res_i)
	    {
	      for(j=i,current_res_j=current_res_i+1;j<m[0].atoms;j++)
		{

		  if(m[0].atm[j].rescount!=current_res_j)
		    {

		      //	     {
		      //    printf("%d %d %d %d %f %f\n",m[0].atm[j].rescount,m[0].atm[i].rescount,i,j,crd(m,i,j),crd(m,j,i));
		      //printf("%d %s %s %i\n" ,m[0].atm[i].rescount,m[0].atm[i].name,m[0].atm[i].residue,m[0].atm[i].resnum);
		      //printf("%d %s %s %i\n" ,m[0].atm[j].rescount,m[0].atm[j].name,m[0].atm[j].residue,m[0].atm[j].resnum);
		      //printf("%d %d %d %d %f %f\n",m[0].atm[j].rescount,m[0].atm[i].rescount,i,j,crd(m,i,j),crd(m,j,i));
		      //if(abs(m[0].atm[i].rescount-m[0].atm[j].rescount)>5 &&
		      if(strcmp(m[0].atm[j].chain,m[0].atm[i].chain)!=0) {
			//		printf("HEJ %d %d %d %d %f %f\n",current_res_j,current_res_i,i,j,crd(m,i,j),crd(m,j,i));
			d=crd(m,i,j);
			//			printf("%d %d %d %d %s %s %f\n",m[0].atm[j].resnum,m[0].atm[i].resnum,i,j,m[0].atm[j].chain,m[0].atm[i].chain,d);
			if(verbose)
			  printf("ALLNATIVE: %d%s %d%s %f %d %d\n",m[0].atm[i].resnum,m[0].atm[i].chain,m[0].atm[j].resnum,m[0].atm[j].chain,sqrt(d),current_res_i,current_res_j); //,contacts[current_res_i][current_res_j]);
			chainA=m[0].atm[i].chain[0];
			chainB=m[0].atm[j].chain[0];
			if(/*contacts[current_res_i][current_res_j]==0 && */crd(m,i,j)<=cutoff)
			  {
			    nativeA[interface_contacts]=m[0].atm[i].resnum;
			    nativeB[interface_contacts]=m[0].atm[j].resnum;
			    interface_contacts++;

			    d=crd(m,i,j);
			    printf("NATIVE: %d%s %d%s %f\n",m[0].atm[i].resnum,m[0].atm[i].chain,m[0].atm[j].resnum,m[0].atm[j].chain,sqrt(d));

			    /*contacts[current_res_i][current_res_j]=1;
			      contacts[current_res_j][current_res_i]=1;*/
			    
			  //res_contacts[get_res(m[0].atm[i].residue)][get_res(m[0].atm[j].residue)]++;
			  //tot_res_contacts++;
			  }
			/*dist[current_res_i][current_res_j]=d;
			  dist[current_res_j][current_res_i]=d;*/
		      }
		      current_res_j++;
		    }
		}				     
	      current_res_i++;
	    }
	}
      //chain_name=m[0].atm[m[0].CA_ref[i]].chain[0];
    }


  error2=read_molecules(model,set);
  if(error2==0)
    {
 
      residues2=model[0].residues;
      /*
      contacts2 = malloc(residues2*sizeof(int *));
      for(i=0;i<residues2;i++) {
	contacts2[i] = malloc(residues2*sizeof(int));
	}
      for(i=0;i<residues2;i++) {
	for(i=0;i<residues2;i++) {
	  contacts2[i][j]=0;
	}
      }
      */  
      current_res_i=0;
      for(i=0;i<model[0].atoms;i++)   
	{
	  if(model[0].atm[i].rescount!=current_res_i)
	    {
	      for(j=i,current_res_j=current_res_i+1;j<model[0].atoms;j++)
		{
		  if(model[0].atm[j].rescount!=current_res_j)
		    {

		      //	     {
		      //    printf("%d %d %d %d %f %f\n",model[0].atm[j].rescount,model[0].atm[i].rescount,i,j,crd(m,i,j),crd(m,j,i));
		      //printf("%d %s %s %i\n" ,model[0].atm[i].rescount,model[0].atm[i].name,model[0].atm[i].residue,model[0].atm[i].resnum);
		      //printf("%d %s %s %i\n" ,model[0].atm[j].rescount,model[0].atm[j].name,model[0].atm[j].residue,model[0].atm[j].resnum);
		      //printf("%d %d %d %d %f %f\n",model[0].atm[j].rescount,model[0].atm[i].rescount,i,j,crd(m,i,j),crd(m,j,i));
		      //if(abs(model[0].atm[i].rescount-model[0].atm[j].rescount)>5 &&
		      if(strcmp(model[0].atm[j].chain,model[0].atm[i].chain)!=0) {

			//printf("HEJ %d %d %d %d %f %f\n",current_res_j,current_res_i,i,j,crd(m,i,j),crd(m,j,i));
			//d=crd(m,i,j);
			//printf("%d %d %d %d %s %s %f\n",model[0].atm[j].resnum,model[0].atm[i].resnum,i,j,model[0].atm[j].chain,model[0].atm[i].chain,d);
			//d=crd(model,i,j);
			//printf("ALLMODEL: %d%s %d%s %f %d %d\n",model[0].atm[i].resnum,model[0].atm[i].chain,model[0].atm[j].resnum,model[0].atm[j].chain,sqrt(d),current_res_i,current_res_j);
		
			chainA2=model[0].atm[i].chain[0];
			chainB2=model[0].atm[j].chain[0];
			if(/*contacts2[current_res_i][current_res_j]==0 &&*/ crd(model,i,j)<=cutoff)
			  {
			    //printf("%d %d %d %d %d %d %d\n", residues2,model[0].atm[i].resnum,model[0].atm[j].resnum,model[0].atm[i].rescount,model[0].atm[j].rescount,current_res_i,current_res_j);
			    modelA[interface_contacts_model]=model[0].atm[i].resnum;
			    modelB[interface_contacts_model]=model[0].atm[j].resnum;
			    interface_contacts_model++;
			    			    
			    if(verbose) {
			      d=crd(model,i,j);
			      printf("MODEL: %d%s %d%s %f %d %d\n",model[0].atm[i].resnum,model[0].atm[i].chain,model[0].atm[j].resnum,model[0].atm[j].chain,sqrt(d),current_res_i,current_res_j);
			    }
			    /*contacts2[current_res_i][current_res_j]=1;
			      contacts2[current_res_j][current_res_i]=1;*/
			  //res_contacts[get_res(model[0].atm[i].residue)][get_res(model[0].atm[j].residue)]++;
			  //tot_res_contacts++;
			  }
		      }
		      current_res_j++;
		    }
		}				     
	      current_res_i++;
	    }
	}
      //chain_name=model[0].atm[model[0].CA_ref[i]].chain[0];
    }

  int matches=0;
  if(chainA!= chainA2 || chainB != chainB2) {
    fprintf(stderr,"chain mismatch %c %c %c %c",chainA,chainA2,chainB,chainB2);
  } else {


    // find matches in nativeA-nativeB modelA-modelB
    for(i=0;i<interface_contacts;i++) {
      for(j=0;j<interface_contacts_model;j++) {
	if(nativeA[i] == modelA[j] &&
	   nativeB[i] == modelB[j]) {
	  if(verbose)
	    printf("OverlapMODEL-NATIVE %d%c %d%c\n",nativeA[i],chainA,nativeB[i],chainB);
	  matches++;

	}
      }

    }

  }

  if(interface_contacts!=0) { 
    printf("Fnat %d %d %f\n",matches,interface_contacts,(float)matches/(float)interface_contacts);
  } else {
    printf("Fnat %d %d %f\n",0,0,0.0);
  }
  if(interface_contacts_model!=0) { 
    printf("Fnonnat %d %d %f\n",interface_contacts_model-matches,interface_contacts_model,(float)(interface_contacts_model-matches)/(float)interface_contacts_model);
  } else {
    printf("Fnonnat %d %d %f\n",0,0,0.0);
  }
  

}

