#include "ESL/LikelihoodHelper.h"


LikelihoodHelper::LikelihoodHelper( DataManager * _dm , int skg )
{
   dm = _dm;
   skgen = skg;

   LoadWrappers();

   tot_energy = 0.0;
   mering = 0;
}


void LikelihoodHelper::LoadWrappers()
{
   dm->Get("nmue"    , nmue);
   dm->Get("evis"    , evis);
   dm->Get("etime"   , etime);
   dm->Get("etype"   , etype);
   dm->Get("ehit"    , ehit);
   dm->Get("egood"   , egood);

   dm->Get("evis"   , evis    );
   dm->Get("nhitac" , nhitac  );
   dm->Get("wall"   , wall    );
   dm->Get("nring"  , nring   );

   dm->Get("ip", ip);
   dm->Get("amome", amome);
   dm->Get("prmslg" , prmslg  ); 
}

bool LikelihoodHelper::FCFVCut()
{
  int nhitac_cut[] = { 10, 16, 16, 16 };

  bool kPass = false ;
  if ( evis(0) > 30.0 && wall(0) > 200.0 && nhitac(0) <  nhitac_cut[skgen] )
   kPass = true;

  return kPass;

}

bool LikelihoodHelper::SingleRingCut()
{
  return ( nring(0) == 1 ? true : false );
}


bool LikelihoodHelper::SubGeVCut()
{
  return ( evis(0) < 1330.0 ? true : false );
}

//////////
//
//  This is the Hayakawa cut
//  criteria for SK-I,-II, and -III
//  otherwise nmue is returned
////////////////
int LikelihoodHelper::muedcyBuild(  )
{

   int i , nmuemax;
   int lmuedcy = 0;

   int ehit_cut_1[4] = { 60 , 30 , 60 , 60 };
   int ehit_cut_2[4] = { 40 , 20 , 40 , 40 };



   nmuemax = ( nmue(0) > 10 ? 10 :  nmue(0) );
   if ( nmue(0) < 0) return 0;

   for( i = 0 ; i < nmuemax ; i++ )
   {
      if(skgen ==3) //for SK4 only
      {
          if(etime(i)<0.1) continue;
	  lmuedcy++;
      }
      else
      {
          if(  evis(0) <  1330.0 && etime(i) < 0.1 ) continue;
          if(  evis(0) >= 1330.0 && etime(i) < 1.2 ) continue;
          if(  etime(i) > 0.8    && etime(i) < 1.2 ) continue;

          if( etype(i) == 1 && ehit(i) >= ehit_cut_1[skgen] && egood(i) > 0.5 )
          {
             lmuedcy++;
          }
          else if( etype(i) >= 2 && etype(i) <= 4 && ehit(i) >= ehit_cut_2[skgen] )
          {
             lmuedcy++;
          }
      }
   }// end of loop on rings


   return lmuedcy;

}


int LikelihoodHelper::mipBuild( )
{

  int   i = 0 ;
  int   lmip = 0 ;
  
  tot_energy = 0.0;
  mering = 0;
  emax = 0.0;

  lmip = ip(0);

  for( i = 0 ; i < nring(0) ; i ++ )
  { 
     tot_energy += amome(i);
     if ( amome(i) > emax )  
     {  
        mering = i;
        emax = amome(i);
     }// end of mering selection
   }// end of ring loop

  // e-like
  if ( prmslg(mering,1) <=  prmslg(mering,2) ) 
     lmip=2;
  else // mu-like
     lmip=3;

  return lmip;
}
