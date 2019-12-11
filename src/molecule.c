#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include "molecule.h"
//#include "nr.h"
//#include "nrutil.h"
//#include "nets.h"


/* changed input to (molecule *m,char atomflag)
   atomflag == a -> read all atoms (except H)
   atomflag == c -> CA atoms
   atomflag == b -> backbone CA,C,N,O atoms
*/

int read_molecules(molecule *m,char atomflag)	/* Reads in molecules to be superimposed */
{
  int	i,j,k,atoms,residues;	/* Counter variables */
  char	buff[1200];	/* Input string */
  char	temp[1000];	/* Throwaway string */
  char	line_flag[11];	/* PDB file line mode */
  char  desc_flag[20];
  char	residue[4];	/* PDB atom info for output */
  char	name[4];
  int   resnum;
  char  resname[9];
  char  old_resname[9]="noname";
  char alt_loc_check[2]=" ";
  char x_temp[11];
  char y_temp[11];
  char z_temp[11];
  //char	number[8];
  int number;
  char chain[2];
  char inscode[2]=" ";
  char alt_loc[2]=" ";
  double	x,y,z;		/* Temporary coordinates values */
  FILE *fp;
  double bfactor=9.99;
  char temp_number[11];
  char temp_resnum[11];
  char temp_bfactor[11];
  int ss_flag=0;

  i=0; /* It is only one molecule to be read this was done for the moment instead of changing all "i" to 0 */
  

  //for (i=0;i<2;i++)
  //{
  //#ifdef  NOZLIB*/
  //fp=fopen(m[i].filename,"r");	/* Does file exist? */
  //#else
  //fp=gzopen(m[i].filename,"r");	/* Does file exist? */
  //#endif
  //printf("%s\n",m[0].filename);
  strcpy(m[0].method,"undef");
  m[0].rank=-1;
  m[0].score=-9999;

  fp=fopen(m[0].filename,"r");	/* Does file exist? */
  if (fp!=NULL)	/* If yes, read in coordinates */
    {
      /* Initialize things */
      //m[0].xcen=m[0].ycen=m[0].zcen=0;
      atoms=0;
      residues=0;
      //#ifdef NOZLIB
      //while(fgets(buff,255,fp)!=NULL)
      //#else
      //while(gzgets(fp,buff,255)!=NULL)
      //#endif
      while(fgets(buff,1000,fp)!=NULL)
	{
	  //sscanf(buff,"%-6s %4d  %-3s%1s%3s %1s%5s    %7.3lf %7.3lf %7.3lf",line_flag,&number,name,residue,chain,&resnum,inscode,&x,&y,&z);
	  //sscanf(buff,"%s %d %s %s %s",line_flag,&number,name,residue,resname);
	  strcpy(line_flag,"undef");
	  sscanf(buff,"%s",line_flag);
	  //printf("%s %s %s",m[0].filename,line_flag,&buff);
	  if(strcmp("REMARK",line_flag)==0)
	    {
	      sscanf(buff,"%s %s %s",line_flag,desc_flag,temp);
	      if(strcmp("SS",desc_flag)==0)
		{
		  strcpy(m[0].ss,temp);
		  ss_flag=1;
		  //printf("%s\n",m[0].ss);
		}
	      if(strcmp("METHOD",desc_flag)==0)
		{
		  strcpy(m[0].method,temp);
		  //printf("%s\n",m[0].method);
		}
	      if(strcmp("SCORE",desc_flag)==0)
		{
		  m[0].score=atof(temp);
		  //strcpy(m[0].method,temp);
		  //printf("%lf\n",m[0].score);
		}
    	    }
	  if(strcmp("MODEL",line_flag)==0)
	    {
	      sscanf(buff,"%s %s",line_flag,temp);
	      m[0].rank=atoi(temp);
	    }
	  
	  //ATOM   5307 3HD1 ILE   340     -27.009  -9.984  11.751  1.00  0.00
	  if(strcmp("ATOM",line_flag)==0 && buff[12] != 'H' && buff[13] != 'H')	/* Is it an ATOM entry? */
	    {
	      //printf("Heja: %s %s",m[0].filename,&buff);
	      strncpy_NULL(temp_number,&buff[6],5);
	      strncpy_NULL(name,&buff[13],3);
	      strncpy_NULL(residue,&buff[17],3);
	      if(atomflag == 'a' ||
		 (atomflag == 'c' && strcmp("CB ",name) == 0) ||
		 (atomflag == 'c' && strcmp("GLY",residue) == 0 && strcmp("CA ",name) == 0) || 
		 (atomflag == 'b' && (strcmp("CA ",name) == 0 || strcmp("C  ",name) == 0 || strcmp("O  ",name) == 0 || strcmp("N  ",name) ==0)))
		{
		  //printf("%s",buff);
		  // printf("%s %c %c %c %s\n",name,buff[12],buff[13],buff[24],buff);
		  strncpy_NULL(alt_loc,&buff[16],1);
		  strncpy_NULL(chain,&buff[21],1);
		  strncpy_NULL(temp_resnum,&buff[22],4);
		  strncpy_NULL(resname,&buff[22],5);
		  strncpy_NULL(x_temp,&buff[30],8);
		  strncpy_NULL(y_temp,&buff[38],8);
		  strncpy_NULL(z_temp,&buff[46],8);
		  
		  number=atoi(temp_number);
		  resnum=atoi(temp_resnum);
		  if(strlen(buff)>=66)
		    {

		      strncpy_NULL(temp_bfactor,&buff[60],6);
		      bfactor=atof(temp_bfactor);
		      //    printf("%lf\n",bfactor);
		    }

		  x=atof(x_temp);
		  y=atof(y_temp);
		  z=atof(z_temp);
		  //printf("test: %s %d %s %s %s %s %d %s %lf %lf %lf\n",line_flag,number,name,alt_loc,residue,chain,resnum,resname,x,y,z);
		  // printf("test: %s %d %s %s %s %s %d %s %lf %lf %lf\n",line_flag,number,name,alt_loc,residue,chain,resnum,resname,x,y,z);
	      //if (strcmp("N",name)==0)	   /*  Is it an N atom => new residue? */
	      //printf("%s %s\n",old_resname,resname);
		  if(strcmp(old_resname,resname)!=0)
		    {
		      m[0].sequence[residues]=aa321(residue);
		      residues++;
		      strcpy(alt_loc_check,alt_loc);
		      //printf("%s %s\n",resname,residue);
		    }
		  //sscanf(&buff[22],"%d %lf %lf %lf",&resnum,&x,&y,&z);
		  //printf("test %d %s %s %d %lf %lf %lf\n",number,name,residue,resnum,x,y,z);
		  if(strcmp(alt_loc_check,alt_loc)==0)
		    {
		      m[0].atm[atoms].x=x;
		      m[0].atm[atoms].y=y;
		      m[0].atm[atoms].z=z;
		      m[0].atm[atoms].bfactor=bfactor;
		      m[0].atm[atoms].resnum=resnum;
		      m[0].atm[atoms].number=number; //atoi(number);
		      m[0].atm[atoms].rescount=residues;
		      m[0].atm[atoms].selected=TRUE;
		      m[0].atm[atoms].deleted=FALSE;
		      strcpy(m[0].atm[atoms].name,name);
		      strcpy(m[0].atm[atoms].residue,residue);
		      strcpy(m[0].atm[atoms].resname,resname);
		      strcpy(m[0].atm[atoms].chain,chain);
		      
		      if(strcmp("CA ",name) == 0)
			{
			  //printf("%s %s %s %s %d\n",m[0].filename,resname,residue,alt_loc,residues);
			  m[0].CA_ref[residues-1]=atoms;
			}
		      //if(strcmp(old_resname,resname)!=0)
		      //	{
		      //	  m[0].sequence[residues-1]=aa321(residue);
		      //	}
		  
		      m[0].xcen+=x;
		      m[0].ycen+=y;
		      m[0].zcen+=z;
		      atoms++;
		      strcpy(old_resname,resname);
		    }
		}
	    }
	}
      m[0].sequence[residues]='\0';
      m[0].atoms=atoms;
      m[0].residues=residues;
      m[0].xcen=m[0].xcen/atoms;
      m[0].ycen=m[0].ycen/atoms;
      m[0].zcen=m[0].zcen/atoms;
      if(ss_flag == 0)
	{
	  //fprintf(stderr,"No secondary structure information in file!\n");
	}

      fclose(fp);
      //#ifdef NOZLIB
      //  fclose(fp);
      //#else
      //	  gzclose(fp);
      //#endif
	  //	  if (atoms!=m[0].atoms)		/* Are file sizes indentical? */
	  //{
	  //  printf("Inconsistent number of atoms in file %s\n",m[0].filename);
	  //  return(1);
	  //}
       
    }
  else
    {
      printf("Couldn't open file \"%s\"\n",m[0].filename);
      exit(1);
    }
  return(0);
}

int read_molecules_ca(molecule *m)	/* Reads in molecules to be superimposed */
{
  int	i,j,k,atoms,residues;	/* Counter variables */
  char	buff[1200];	/* Input string */
  char	temp[1000];	/* Throwaway string */
  char	line_flag[11];	/* PDB file line mode */
  char  desc_flag[20];
  char	residue[4];	/* PDB atom info for output */
  char	name[4];
  int   resnum;
  char  resname[9];
  char  old_resname[9]="noname";
  char alt_loc_check[2]=" ";
  char x_temp[11];
  char y_temp[11];
  char z_temp[11];
  //char	number[8];
  int number;
  char chain[2];
  char inscode[2]=" ";
  char alt_loc[2]=" ";
  double	x,y,z;		/* Temporary coordinates values */
  FILE *fp;
  char temp_number[11];
  char temp_resnum[11];
  
  i=0; /* It is only one molecule to be read this was done for the moment instead of changing all "i" to 0 */
  

  //for (i=0;i<2;i++)
  //{
  //#ifdef  NOZLIB*/
  //fp=fopen(m[i].filename,"r");	/* Does file exist? */
  //#else
  //fp=gzopen(m[i].filename,"r");	/* Does file exist? */
  //#endif
  //printf("%s\n",m[0].filename);
  strcpy(m[0].method,"undef");
  m[0].rank=-1;
  m[0].score=-9999;
  fp=fopen(m[0].filename,"r");	/* Does file exist? */
  if (fp!=NULL)	/* If yes, read in coordinates */
    {
      /* Initialize things */
      //m[0].xcen=m[0].ycen=m[0].zcen=0;
      atoms=0;
      residues=0;
      //#ifdef NOZLIB
      //while(fgets(buff,255,fp)!=NULL)
      //#else
      //while(gzgets(fp,buff,255)!=NULL)
      //#endif
      while(fgets(buff,1000,fp)!=NULL)
	{
	  //sscanf(buff,"%-6s %4d  %-3s%1s%3s %1s%5s    %7.3lf %7.3lf %7.3lf",line_flag,&number,name,residue,chain,&resnum,inscode,&x,&y,&z);
	  //sscanf(buff,"%s %d %s %s %s",line_flag,&number,name,residue,resname);
	  sscanf(buff,"%s",line_flag);
	  if(strcmp("REMARK",line_flag)==0)
	    {
	      //printf("%s",buff);
	      sscanf(buff,"%s %s %s",line_flag,desc_flag,temp);
	      if(strcmp("SS",desc_flag)==0)
		{
		  //printf("%s\n%d\n\n",buff,strlen(buff));
		  strcpy(m[0].ss,temp);
		  //printf("%s\n%d\n",m[0].ss,strlen(m[0].ss));
		}
	      if(strcmp("METHOD",desc_flag)==0)
		{
		  strcpy(m[0].method,temp);
		  //printf("%s\n",m[0].method);
		}
	      if(strcmp("SCORE",desc_flag)==0)
		{
		  m[0].score=atof(temp);
		  //strcpy(m[0].method,temp);
		  //printf("%lf\n",m[0].score);
		}
    
	    }
	  if(strcmp("MODEL",line_flag)==0)
	    {
	      sscanf(buff,"%s %s",line_flag,temp);
	      m[0].rank=atoi(temp);
	    }

	  if (strcmp("ATOM",line_flag)==0 && buff[13] != 'H')	/* Is it an ATOM entry? */
	    { 
	      strncpy_NULL(temp_number,&buff[6],5);
	      strncpy_NULL(name,&buff[13],3);
	      //printf("%s",&buff[6]);
	      if(strcmp("CA ",name) == 0)
		{
	      //printf("%s",&buff[6]);
		  strncpy_NULL(alt_loc,&buff[16],1);
		  strncpy_NULL(residue,&buff[17],3);
		  strncpy_NULL(chain,&buff[21],1);
		  strncpy_NULL(temp_resnum,&buff[22],4);
		  strncpy_NULL(resname,&buff[22],5);
		  strncpy_NULL(x_temp,&buff[30],8);
		  strncpy_NULL(y_temp,&buff[38],8);
		  strncpy_NULL(z_temp,&buff[46],8);
		  
		  number=atoi(temp_number);
		  resnum=atoi(temp_resnum);
		  x=atof(x_temp);
		  y=atof(y_temp);
		  z=atof(z_temp);
		
	      
	      //printf("test: %s %d %s %s %s %s %d %s %lf %lf %lf\n",line_flag,number,name,alt_loc,residue,chain,resnum,resname,x,y,z);
	      //if (strcmp("N",name)==0)	   /*  Is it an N atom => new residue? */
	      //printf("%s %s\n",old_resname,resname);
		  if(strcmp(old_resname,resname)!=0)
		  {
		    m[0].sequence[residues]=aa321(residue);
		    residues++;
		    strcpy(alt_loc_check,alt_loc);
		    //printf("%s %s\n",resname,residue);
		  }
	      //sscanf(&buff[22],"%d %lf %lf %lf",&resnum,&x,&y,&z);
		  //printf("test %d %s %s %d %lf %lf %lf\n",number,name,residue,resnum,x,y,z);
		  if(strcmp(alt_loc_check,alt_loc)==0)
		    {
		      m[0].atm[atoms].x=x;
		      m[0].atm[atoms].y=y;
		      m[0].atm[atoms].z=z;
		      m[0].atm[atoms].resnum=resnum;
		      m[0].atm[atoms].number=number; //atoi(number);
		      m[0].atm[atoms].rescount=residues;
		      m[0].atm[atoms].selected=TRUE;
		      strcpy(m[0].atm[atoms].name,name);
		      strcpy(m[0].atm[atoms].residue,residue);
		      //if(strcmp(old_resname,resname)!=0)
		      //	{
		      //	  m[0].sequence[residues-1]=aa321(residue);
		      //	}
		      
		      m[0].xcen+=x;
		      m[0].ycen+=y;
		      m[0].zcen+=z;
		      atoms++;
		      strcpy(old_resname,resname);
		    }
		}
	    }
	}
      m[0].sequence[residues]='\0';
      m[0].atoms=atoms;
      m[0].residues=residues;
      m[0].xcen=m[0].xcen/atoms;
      m[0].ycen=m[0].ycen/atoms;
      m[0].zcen=m[0].zcen/atoms;
      fclose(fp);
      //#ifdef NOZLIB
      //  fclose(fp);
      //#else
      //	  gzclose(fp);
      //#endif
	  //	  if (atoms!=m[0].atoms)		/* Are file sizes indentical? */
	  //{
	  //  printf("Inconsistent number of atoms in file %s\n",m[0].filename);
	  //  return(1);
	  //}
       
    }
  else
    {
      printf("Couldn't open file \"%s\"\n",m[0].filename);
      exit(1);
    }
  return(0);
}


int read_molecules_backbone(molecule *m)	/* Reads in molecules to be superimposed */
{
  int	i,j,k,atoms,residues;	/* Counter variables */
  char	buff[1200];	/* Input string */
  char	temp[1000];	/* Throwaway string */
  char	line_flag[11];	/* PDB file line mode */
  char  desc_flag[20];
  char	residue[4];	/* PDB atom info for output */
  char	name[4];
  int   resnum;
  char  resname[9];
  char  old_resname[9]="noname";
  char x_temp[11];
  char y_temp[11];
  char z_temp[11];
  //char	number[8];
  int number;
  char chain[2];
  char inscode[2]=" ";
  char alt_loc[2]=" ";
  double	x,y,z;		/* Temporary coordinates values */
  FILE *fp;
  char temp_number[11];
  char temp_resnum[11];
  i=0; /* It is only one molecule to be read this was done for the moment instead of changing all "i" to 0 */
  //printf("%s\n",m[0].filename);
  strcpy(m[0].method,"undef");
  m[0].rank=-1;
  m[0].score=-9999;

  fp=fopen(m[0].filename,"r");	/* Does file exist? */
  if (fp!=NULL)	/* If yes, read in coordinates */
    {
      /* Initialize things */
       atoms=0;
      residues=0;
       while(fgets(buff,1000,fp)!=NULL)
	{
	  //sscanf(buff,"%-6s %4d  %-3s%1s%3s %1s%5s    %7.3lf %7.3lf %7.3lf",line_flag,&number,name,residue,chain,&resnum,inscode,&x,&y,&z);
	  //sscanf(buff,"%s %d %s %s %s",line_flag,&number,name,residue,resname);
	  sscanf(buff,"%s",line_flag);
	  if(strcmp("REMARK",line_flag)==0)
	    {
	      sscanf(buff,"%s %s %s",line_flag,desc_flag,temp);
	      if(strcmp("SS",desc_flag)==0)
		{
		  //strcpy(m[0].ss,temp);
		  strncpy_NULL(m[0].ss,temp,strlen(temp));
		  //  printf("%s\n",m[0].ss);
		}
	      if(strcmp("METHOD",desc_flag)==0)
		{
		  //strcpy(m[0].method,temp);
		  strncpy_NULL(m[0].method,temp,strlen(temp));
		  //printf("%s\n",m[0].method);
		}
	      if(strcmp("SCORE",desc_flag)==0)
		{
		  m[0].score=atof(temp);
		  //strcpy(m[0].method,temp);
		  //printf("%lf\n",m[0].score);
		}
    
	    }
	  if(strcmp("MODEL",line_flag)==0)
	    {
	      sscanf(buff,"%s %s",line_flag,temp);
	      m[0].rank=atoi(temp);
	    }
	  if(strcmp("ATOM",line_flag)==0)	/* Is it an ATOM entry? */
	    {
	      //printf("%s",&buff[6]);
	      strncpy_NULL(temp_number,&buff[6],5);
	      strncpy_NULL(name,&buff[13],3);
	      //printf("%s",&buff[6]);
	      if(strcmp("CA ",name) == 0 ||
		 strcmp("C  ",name) == 0 || 
		 strcmp("O  ",name) == 0 ||
		 strcmp("N  ",name) ==0)
		{
		  //printf("%d %s",2,&buff[6]);
		  strncpy_NULL(alt_loc,&buff[16],1);
		  strncpy_NULL(residue,&buff[17],3);
		  strncpy_NULL(chain,&buff[21],1);
		  strncpy_NULL(temp_resnum,&buff[22],4);
		  strncpy_NULL(resname,&buff[22],5);
		  strncpy_NULL(x_temp,&buff[30],8);
		  strncpy_NULL(y_temp,&buff[38],8);
		  strncpy_NULL(z_temp,&buff[46],8);


		  //(strcmp("CA ",name) !=0)
		  //
		  number=atoi(temp_number);
		  resnum=atoi(temp_resnum);
		  x=atof(x_temp);
		  y=atof(y_temp);
		  z=atof(z_temp);
		  
		  //printf("test: %s %d %s %s %s %s %d %s %lf %lf %lf\n",line_flag,number,name,alt_loc,residue,chain,resnum,resname,x,y,z);
		  //if (strcmp("N",name)==0)	   /*  Is it an N atom => new residue? */
		  //printf("%s %s\n",old_resname,resname);
		  if(strcmp(old_resname,resname)!=0)
		    {
		      m[0].sequence[residues]=aa321(residue);
		      m[0].res_ref[residues]=atoms;
		      residues++;
		      //printf("%s %s\n",resname,residue);
		    }
		      //printf("%s %s\n",resname,residue);
		      //  }
		  //sscanf(&buff[22],"%d %lf %lf %lf",&resnum,&x,&y,&z);
		      //printf("test %d %s %s %d %lf %lf %lf\n",number,name,residue,resnum,x,y,z);
		  m[0].atm[atoms].x=x;
		  m[0].atm[atoms].y=y;
		  m[0].atm[atoms].z=z;
		  m[0].atm[atoms].resnum=resnum;
		  //printf("%d\n",m[0].atm[atoms].resnum);
		  m[0].atm[atoms].number=number; //atoi(number);
		  m[0].atm[atoms].rescount=residues;
		  m[0].atm[atoms].selected=TRUE;
		  strcpy(m[0].atm[atoms].name,name);
		  strcpy(m[0].atm[atoms].residue,residue);
		  if(strcmp("CA ",name) == 0)
		    {
		      //printf("%s %s %s %d\n",resname,residue,alt_loc,residues);
		      m[0].CA_ref[residues-1]=atoms;
		    }

		  //  if(strcmp(old_resname,resname)!=0)
		  // {
		  //    m[0].sequence[residues-1]=aa321(residue);
		  //  }

		  m[0].xcen+=x;
		  m[0].ycen+=y;
		  m[0].zcen+=z;
		  atoms++;
		  strcpy(old_resname,resname);
		  //}
		}
	    }
	}
       m[0].sequence[residues]='\0';
       m[0].atoms=atoms;
       m[0].residues=residues;
       m[0].xcen=m[0].xcen/atoms;
       m[0].ycen=m[0].ycen/atoms;
       m[0].zcen=m[0].zcen/atoms;


      fclose(fp);
      //#ifdef NOZLIB
      //  fclose(fp);
      //#else
      //	  gzclose(fp);
      //#endif
	  //	  if (atoms!=m[0].atoms)		/* Are file sizes indentical? */
	  //{
	  //  printf("Inconsistent number of atoms in file %s\n",m[0].filename);
	  //  return(1);
	  //}
       
    }
  else
    {
      printf("Couldn't open file \"%s\"\n",m[0].filename);
      exit(1);
    }
  return(0);
}


void strncpy_NULL(char *dest, char *src, size_t n)
{
  strncpy(dest, src, n);
  dest[n]='\0';
}






double crd(molecule *m,int atomno1, int atomno2)       /* atomnoX is the first atom of the residue */
{
  int i,j;
  double dist,lowest_dist;
  lowest_dist=999999;
  //printf("crd i: %d, j; %d\n",atomno1,atomno2);
  for(i=atomno1;m[0].atm[i].rescount == m[0].atm[atomno1].rescount;i++)  
    {
      //      if(strcmp("C  ",m[0].atm[i].name)!=0 && 
      //	 strcmp("O  ",m[0].atm[i].name)!=0 &&
      //	 strcmp("N  ",m[0].atm[i].name)!=0)
	{
	  for(j=atomno2;m[0].atm[j].rescount == m[0].atm[atomno2].rescount;j++)
	    {
	      //      if(strcmp("C  ",m[0].atm[j].name)!=0 && 
	      // strcmp("O  ",m[0].atm[j].name)!=0 &&
	      //	 strcmp("N  ",m[0].atm[j].name)!=0)
		{
		  dist=distance(m,i,j);
		  //printf("%f ",dist);
		  if(dist<lowest_dist)
		    lowest_dist=dist;
		  //printf("%s %s %f\n",m[0].atm[j].residue,m[0].atm[j].name,lowest_dist);
		}
	    }
	}
    }
  //printf("%s %s %d %d\n",m[0].atm[i].residue,m[0].atm[i].name,m[0].atm[i].rescount,m[0].atm[atomno1].rescount);
  //printf("%s %s %d %d\n",m[0].atm[j].residue,m[0].atm[j].name,m[0].atm[j].rescount,m[0].atm[atomno2].rescount);
  //  printf("lowest dist: %f\n",lowest_dist);
  return lowest_dist;
}





double distance(molecule *m,int atomno1,int atomno2)
{
  //printf("%s %s\n",m[0].atm[atomno1-1].name,m[0].atm[atomno2-1].name);

  return (m[0].atm[atomno1].x-m[0].atm[atomno2].x)*(m[0].atm[atomno1].x-m[0].atm[atomno2].x)+(m[0].atm[atomno1].y-m[0].atm[atomno2].y)*(m[0].atm[atomno1].y-m[0].atm[atomno2].y)+(m[0].atm[atomno1].z-m[0].atm[atomno2].z)*(m[0].atm[atomno1].z-m[0].atm[atomno2].z);
}


char aa321(char* res)
{
  //char check_res[4];

  //check_res[0]=res[0];
  //check_res[1]=res[1];
  //check_res[2]=res[2];
  //res[3]='\0';
  //printf("%s %s\n",res,res);
  
  if(strcmp("ALA",res)==0)
    return 'A';
  if(strcmp("ARG",res)==0)
    return 'R';
  if(strcmp("ASN",res)==0)
    return 'N';
  if(strcmp("ASP",res)==0)
    return 'D';
  if(strcmp("CYS",res)==0)
    return 'C';
  if(strcmp("GLN",res)==0)
    return 'Q';
  if(strcmp("GLU",res)==0)
    return 'E';
  if(strcmp("GLY",res)==0)
    return 'G';
  if(strcmp("HIS",res)==0)
    return 'H';
  if(strcmp("ILE",res)==0)
    return 'I';
  if(strcmp("LEU",res)==0)
    return 'L';
  if(strcmp("LYS",res)==0)
    return 'K';
  if(strcmp("MET",res)==0)
    return 'M';
  if(strcmp("PHE",res)==0)
    return 'F';
  if(strcmp("PRO",res)==0)
    return 'P';
  if(strcmp("SER",res)==0)
    return 'S';
  if(strcmp("THR",res)==0)
    return 'T';
  if(strcmp("TRP",res)==0)
    return 'W';
  if(strcmp("TYR",res)==0)
    return 'Y';
  if(strcmp("VAL",res)==0)
    return 'V';
  return 'X';
}
