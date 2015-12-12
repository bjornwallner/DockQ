#ifndef molecule_h
#define molecule_h
//#include "types/simple.h"
//#include "nets.h"
#define TRUE		1		/* Boolean definitions */
#define FALSE		0
#define	MAXATMS		25000		/* Maximum allowable atoms */
#define	MAXRES		3000	        /* Maximum allowable residues */
#define PI		3.14159265	/* Useful constant */
#undef max
#define max(a,b)    ((a) > (b) ? (a) : (b))
#undef min
#define min(a,b)    ((a) < (b) ? (a) : (b))

//#define SIZE (sizeof(a) / sizeof(a[0]))

typedef struct
 {
   double x,y,z;		/* Atomic coordinates */
   double rms;		/* RMS deviation */
   char residue[8];	/* PDB info for output */
   char name[8];
   int number;
    int resnum;
   char resname[8];
   int rescount;
   int selected;
   int deleted;
   double bfactor;
   char chain[2];
 } atm;

typedef struct {
  atm atm[MAXATMS];
  int CA_ref[MAXRES];
  int res_ref[MAXRES];
  double xcen,ycen,zcen;
  int	atoms;			/* Current # of atoms */
  int   residues;
  char  sequence[MAXRES];
  char  ss[MAXRES];
  double score;
  char method[200];
  int  rank;
  char	filename[1000];		/* filename to read molecule from */
  
} molecule;

typedef struct 
{
  double x,y,z;
}my_vector;


/* OBS!!!!!!!!!
   To relace *_ca and _backbone with one routine I
   changed input to read_molecules(molecule *m,char atomflag)
   atomflag == a -> read all atoms (except H)
   atomflag == c -> CA atoms
   atomflag == b -> backbone CA,C,N,O atoms
*/

enum aacodes
{
    ALA, ARG, ASN, ASP, CYS,
    GLN, GLU, GLY, HIS, ILE,
    LEU, LYS, MET, PHE, PRO,
    SER, THR, TRP, TYR, VAL,
    UNK
};


int    read_molecules(molecule *m,char atomflag);
int    read_molecules_ca(molecule *m);
int    read_molecules_backbone(molecule *m);

void   strncpy_NULL(char *dest, char *src, size_t n);

double distance(molecule *m,int atomno1, int atomno2);
double crd(molecule *m,int atomno1, int atomno2);   /*closest residue distance */
char   aa321(char *name);


#endif
