#ifndef _outputOscStructure_
#define _outputOscStructure_



class outputOscStructure
{

  public :

  int         nnth;
  int         nnthb;
  int         nnm10;
  int         nnpb10;
  int         ntot;
  int         nmax;
  int         nn;
  int         nnew;
  int         nn_mctruth;
  int         nmue;
  int         nring;
  float         etim;
  float         edist;
  float         p_ring;
  float         pid;
  float         rapidity;
  float         efrac;
  float         dirn[3];
  float         ddir;
  
  
   int           ipnu;
   float         dirnu[3];
   float         pnu;
   int           mode;
   int           ip;
   float         dprob;
   float         dir[3];
   float         amom;
   float         amom_pri;
   float         amom_forcorr;
   float         path;
   float         wall;
   int           itype;
   int           muedk;
   float         flxg[3];
   float         flxgo[3];
   float         flxh[3];
   float         flxho[3];
   float         weightx;
   float         oscwgt;
   float         etotmue;
   float         fsiweight;
   float         oscweight3f;
   int           nn_selected;
   float         nn_output;

   float         ndist1;
   float         ndist2;
   float         totpid;
   float         evis;
   
   float         amom_corrected;

   
   void Zero(){

      ipnu    = 0  ;
      pnu     = 0. ;
      mode    = 0  ;
      ip      = 0  ;
      dprob   = 0. ;
      amom    = 0. ;
      amom_pri    = 0. ;
      amom_forcorr    = 0. ;
      path    = 0. ;
      wall    = 0. ;
      itype   = 0  ;
      muedk   = 0  ;
      weightx = 0. ;
      oscwgt  = 0. ;
      etotmue = 0. ;
      fsiweight = 0. ;
      oscweight3f = 0. ;
      nn_selected = 0;
      nn_output = 0.;
      amom_corrected = 0.;
 
      for( int i = 0; i < 3 ; i ++ )
      {
         dir[i]   = 0. ;
         dirnu[i] = 0. ;
         flxg[i]  = 0. ;
         flxgo[i] = 0. ;
         flxh[i]  = 0. ;
         flxho[i] = 0. ;
      }
   }


};



#endif
