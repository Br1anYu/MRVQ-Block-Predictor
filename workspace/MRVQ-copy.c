/* ======================================================================== */
/*  MRVQ + Block Mean + Neighbor Block Copy (block prediction)              */
/*  Based on MRVQ-mean(p58).c                                               */
/*                                                                          */
/*  Block modes:                                                            */
/*     0 : MRVQ   (mean 8 bits + index 13 bits)                             */
/*     1 : Copy reconstructed LEFT      block                               */
/*     2 : Copy reconstructed UP        block                               */
/*     3 : Copy reconstructed UP-LEFT   block   (only when NB=4)            */
/*     4 : Copy reconstructed UP-RIGHT  block   (only when NB=4)            */
/*     5 : Block mean (all 16 pixels = mean, 8 bits)  (only when BM=1)      */
/*                                                                          */
/*  Scheme A (NB=2): MRVQ + LEFT + UP       [+ BM]   -> 3 [4] modes, 2 bits */
/*  Scheme B (NB=4): MRVQ + LEFT + UP + UL + UR [+ BM] -> 5 [6] modes, 3 b. */
/*  BM uses a flag code that is otherwise unused, so flag length is equal.  */
/*                                                                          */
/*  Mode selection: minimum  J = dist + LAMBDA * bits                       */
/*                                                                          */
/*  Usage: MRVQ-copy input_image compressed_image [NB] [LAMBDA] [BM]        */
/*         NB     : 2 (default) or 4                                        */
/*         LAMBDA : default 30                                              */
/*         BM     : 0 (default, off) or 1 (enable block-mean mode)          */
/*  (source file is ASCII-only to avoid Big5/UTF-8 trouble)                 */
/* ======================================================================== */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define imagesize  512
#define vectorsize 16
#define blkperrow  128                               /* blocks per row */
#define vectornum  (512*512/16)                      /* no. of vectors */
#define codenum    8192                              /* size of the codebook */

#define lenMEAN    8
#define lenINDEX   13

#define MODE_MRVQ  0
#define MODE_LEFT  1
#define MODE_UP    2
#define MODE_UL    3
#define MODE_UR    4
#define MODE_BM    5
#define MODENUM    6

char               MRCBfile[80]="MRCB_8192.raw";

unsigned char      image[imagesize][imagesize];      /* source image */
unsigned char      Cimage[imagesize][imagesize];     /* decoded image */
unsigned char      vectors[vectornum][vectorsize];   /* source vectors */
unsigned char      Rvectors[vectornum][vectorsize];  /* encoder-side reconstruction */
unsigned char      Cvectors[vectornum][vectorsize];  /* decoder-side reconstruction */
unsigned char      codebook[codenum][vectorsize];

int                IndexTable[vectornum];
int                blockmean[vectornum];             /* mean used by MRVQ / BM blocks */
int                EncodingType[vectornum];
int                cntMode[MODENUM];

/* ---- scheme parameters (set in main) ---- */
int                NB     = 2;        /* number of copy neighbors: 2 or 4   */
double             LAMBDA = 30.0;     /* RD multiplier                      */
int                USEBM  = 0;        /* 1: enable block-mean mode          */
int                lenFLAG;           /* mode flag bits: 2 (NB=2), 3 (NB=4) */
int                modeBits[MODENUM]; /* total bits of each mode            */

/* offset (in block index) of each copy neighbor, indexed by mode-1 */
int                nbOffset[4] = { -1, -blkperrow, -blkperrow-1, -blkperrow+1 };

unsigned long      squar[256];
double             BitRate;           /* bpp, filled by ImageEncoding       */

FILE               *imagefile,*result,*book;

int  ImageEncoding(void);
int  ImageDecoding(void);
int  Find_PSNR(void);
int  NeighborAvailable(int i, int mode);
unsigned long BlockDist(unsigned char *a, unsigned char *b);

int main(int argc, char **argv)
{
 int i,j;

 if (argc < 3 || argc > 6)
    {
      printf("\nError : The number of parameters isn't correct ! \n");
      printf("\nFormat: command input_image compressed_image [NB=2|4] [LAMBDA] [BM=0|1]\n");
      exit(1);
    }

 if (argc>=4) NB=atoi(argv[3]);
 if (NB!=2 && NB!=4)
    {
      printf("\nNB must be 2 or 4.\n");
      exit(1);
    }
 if (argc>=5) LAMBDA=atof(argv[4]);
 if (argc>=6) USEBM=atoi(argv[5]);
 if (USEBM!=0 && USEBM!=1)
    {
      printf("\nBM must be 0 or 1.\n");
      exit(1);
    }

 /* flag bits and per-mode bit cost
    (3 or 4 modes fit in 2 bits, 5 or 6 modes fit in 3 bits) */
 lenFLAG=(NB==2)?2:3;
 modeBits[MODE_MRVQ]=lenFLAG+lenMEAN+lenINDEX;
 for(i=MODE_LEFT;i<=MODE_UR;i++)
    modeBits[i]=lenFLAG;
 modeBits[MODE_BM]=lenFLAG+lenMEAN;

 if ((imagefile=fopen(argv[1],"rb"))==NULL)
    {
      printf("\nUnable to open input image %s.\n",argv[1]);
      exit(1);
    }

 if ((result=fopen(argv[2],"wb"))==NULL)
    {
      printf("\nUnable to open output compressed image %s.\n",argv[2]);
      exit(1);
    }

 if ((book=fopen(MRCBfile,"rb"))==NULL)
    {
      printf("\nUnable to open the codebook %s.\n",MRCBfile);
      exit(1);
    }

 /* get the codebook */
 for(i=0;i<codenum;i++)
    for(j=0;j<vectorsize;j++)
       {
         int c=fgetc(book);
         if (c==EOF)
            {
              printf("\nCodebook %s is too small (need %d bytes).\n",
                     MRCBfile,codenum*vectorsize);
              exit(1);
            }
         codebook[i][j]=(unsigned char)c;
       }
 fclose(book);

 /* get the input image */
 for(i=0;i<imagesize;i++)
    for(j=0;j<imagesize;j++)
       {
         int c=fgetc(imagefile);
         if (c==EOF)
            {
              printf("\nInput image %s is too small (need %d bytes).\n",
                     argv[1],imagesize*imagesize);
              exit(1);
            }
         image[i][j]=(unsigned char)c;
       }
 fclose(imagefile);

 /* image -> 4x4 vectors */
 for(i=0;i<imagesize;i++)
    for(j=0;j<imagesize;j++)
       vectors[j/4+(i/4)*blkperrow][j%4+(i%4)*4]=image[i][j];

 for(i=0;i<256;i++)
    squar[i]=(unsigned long)(i*i);

 ImageEncoding();
 ImageDecoding();
 Find_PSNR();

 return 0;
}


/* squared error between two 16-pixel blocks */
unsigned long BlockDist(unsigned char *a, unsigned char *b)
{
 int j;
 unsigned long d=0;
 for(j=0;j<vectorsize;j++)
    d+=squar[abs((int)a[j]-(int)b[j])];
 return d;
}


/* Is the neighbor of block i given by 'mode' inside the image?
   Only neighbors decoded earlier (raster order) are used. */
int NeighborAvailable(int i, int mode)
{
 int col=i%blkperrow;
 int row=i/blkperrow;

 switch(mode)
    {
      case MODE_LEFT: return (col>0);
      case MODE_UP:   return (row>0);
      case MODE_UL:   return (row>0 && col>0);
      case MODE_UR:   return (row>0 && col<blkperrow-1);
    }
 return 0;
}


int ImageEncoding(void)
{
 int           i,j,k,m;
 int           sum,mean,index,pixel,best;
 int           lastmode;
 unsigned long dist,min;
 double        cost[MODENUM],bestcost;
 unsigned char Tvector[vectorsize];   /* mean-removed vector   */
 unsigned char recMRVQ[vectorsize];   /* MRVQ reconstruction   */
 unsigned char recBM[vectorsize];     /* block-mean reconstruction */
 long          totalbits=0;

 lastmode=(NB==2)?MODE_UP:MODE_UR;    /* last copy mode in use */

 for(i=0;i<vectornum;i++)
    {
      /* ---------- block mean (used by MRVQ and BM modes) ---------- */
      sum=0;
      for(j=0;j<vectorsize;j++)
         sum+=(int)vectors[i][j];
      mean=(sum+vectorsize/2)/vectorsize;   /* rounded, not truncated */
      blockmean[i]=mean;

      /* ---------- mean-removed vector ---------- */
      for(j=0;j<vectorsize;j++)
         {
           pixel=(int)vectors[i][j]-mean+128;
           if (pixel<0) pixel=0;
           if (pixel>255) pixel=255;
           Tvector[j]=(unsigned char)pixel;
         }

      /* ---------- full search in codebook ---------- */
      min=0xFFFFFFFFUL;
      index=0;
      for(j=0;j<codenum;j++)
         {
           dist=0;
           for(k=0;k<vectorsize;k++)
              dist+=squar[abs((int)Tvector[k]-(int)codebook[j][k])];
           if (dist<min)
              {
                min=dist;
                index=j;
              }
         }
      IndexTable[i]=index;

      /* ---------- MRVQ ---------- */
      for(j=0;j<vectorsize;j++)
         {
           pixel=(int)codebook[index][j]-128+mean;
           if (pixel<0) pixel=0;
           if (pixel>255) pixel=255;
           recMRVQ[j]=(unsigned char)pixel;
         }
      cost[MODE_MRVQ]=(double)BlockDist(vectors[i],recMRVQ)
                      +LAMBDA*modeBits[MODE_MRVQ];

      /* ---------- copy modes ----------
         Use the RECONSTRUCTED blocks (Rvectors); the decoder only has these. */
      for(m=MODE_LEFT;m<=MODE_UR;m++)
         {
           cost[m]=1.0e30;
           if (m<=lastmode && NeighborAvailable(i,m))
              {
                cost[m]=(double)BlockDist(vectors[i],Rvectors[i+nbOffset[m-MODE_LEFT]])
                        +LAMBDA*modeBits[m];
              }
         }

      /* ---------- block mean mode (flat block, 8 bits) ---------- */
      for(j=0;j<vectorsize;j++)
         recBM[j]=(unsigned char)mean;
      cost[MODE_BM]=1.0e30;
      if (USEBM)
         cost[MODE_BM]=(double)BlockDist(vectors[i],recBM)
                       +LAMBDA*modeBits[MODE_BM];

      /* ---------- choose the minimum RD cost ---------- */
      best=MODE_MRVQ;
      bestcost=cost[MODE_MRVQ];
      for(m=MODE_LEFT;m<MODENUM;m++)
         if (cost[m]<bestcost)
            {
              bestcost=cost[m];
              best=m;
            }

      EncodingType[i]=best;
      cntMode[best]++;
      totalbits+=modeBits[best];

      /* store the encoder-side reconstruction for later blocks */
      for(j=0;j<vectorsize;j++)
         {
           if      (best==MODE_MRVQ) Rvectors[i][j]=recMRVQ[j];
           else if (best==MODE_BM)   Rvectors[i][j]=recBM[j];
           else                      Rvectors[i][j]=Rvectors[i+nbOffset[best-MODE_LEFT]][j];
         }
   }

 printf("\nScheme: NB=%d BM=%d (%d modes, flag=%d bits), LAMBDA=%.1f",
        NB,USEBM,((NB==2)?3:5)+USEBM,lenFLAG,LAMBDA);
 printf("\nMRVQ=%5d  LEFT=%5d  UP=%5d",
        cntMode[MODE_MRVQ],cntMode[MODE_LEFT],cntMode[MODE_UP]);
 if (NB==4)
    printf("  UP-LEFT=%5d  UP-RIGHT=%5d",cntMode[MODE_UL],cntMode[MODE_UR]);
 if (USEBM)
    printf("  BM=%5d",cntMode[MODE_BM]);
 BitRate=(double)totalbits/(imagesize*imagesize);
 printf("\nTotal bits = %ld, Bit Rate = %8.4f bpp",totalbits,BitRate);
 return 0;
}


/* The decoder uses ONLY EncodingType / blockmean / IndexTable / codebook
   and already-decoded blocks, so it verifies encoder/decoder consistency. */
int ImageDecoding(void)
{
 int i,j;
 int mean,index,pixel,mode;

 for(i=0;i<vectornum;i++)
    {
      mode=EncodingType[i];

      if (mode==MODE_MRVQ)
         {
           index=IndexTable[i];
           mean=blockmean[i];
           for(j=0;j<vectorsize;j++)
              {
                pixel=(int)codebook[index][j]-128+mean;
                if (pixel<0) pixel=0;
                if (pixel>255) pixel=255;
                Cvectors[i][j]=(unsigned char)pixel;
              }
         }
      else if (mode==MODE_BM)
         {
           for(j=0;j<vectorsize;j++)
              Cvectors[i][j]=(unsigned char)blockmean[i];
         }
      else
         {
           for(j=0;j<vectorsize;j++)
              Cvectors[i][j]=Cvectors[i+nbOffset[mode-MODE_LEFT]][j];
         }
    }

 /* vectors -> image */
 for(i=0;i<vectornum;i++)
    for(j=0;j<vectorsize;j++)
       Cimage[j/4+(i/blkperrow)*4][j%4+(i%blkperrow)*4]=Cvectors[i][j];

 fwrite(Cimage,sizeof(unsigned char),imagesize*imagesize,result);
 fclose(result);
 return 0;
}


int Find_PSNR(void)
{
 int       i,j;
 long long dist=0LL;                  /* long is only 32 bits on Windows */
 double    ave,psnr;

 for(i=0;i<imagesize;i++)
    for(j=0;j<imagesize;j++)
       dist+=squar[abs((int)image[i][j]-(int)Cimage[i][j])];

 ave=(double)dist/(imagesize*imagesize);
 psnr=(ave>0.0)?10*log10(255.0*255.0/ave):99.99;   /* guard: lossless */
 printf("\ntotal_dist=%10lld avg_dist=%8.3f  PSNR=%8.3f\n",dist,ave,psnr);

 /* one machine-readable line for batch scripts */
 printf("RESULT,%d,%d,%.2f,%d,%d,%d,%d,%d,%d,%.4f,%.3f\n",
        NB,USEBM,LAMBDA,cntMode[MODE_MRVQ],cntMode[MODE_LEFT],cntMode[MODE_UP],
        cntMode[MODE_UL],cntMode[MODE_UR],cntMode[MODE_BM],BitRate,psnr);
 return 0;
}
