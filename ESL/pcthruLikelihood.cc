#include "ESL/pcthruLikelihood.h"
#include <fstream>
#include "TTree.h"

using namespace std;

pcthruLikelihood::pcthruLikelihood( DataManager * _dm , CardReader * reader , int skg )
:LikelihoodReturn( _dm, reader )
{
   skgen = skg;

   LoadWrappers();
  
   std::string input;

   read->GetKey("pcthruLikelihood", input );

   TFile * file = new TFile( input.c_str()  ); 
   
      SetLikelihoodFile( file ); 
}



pcthruLikelihood::pcthruLikelihood( DataManager * _dm , CardReader * reader , int skg, const char * varname )
:LikelihoodReturn( _dm, reader )
{
   skgen = skg;

   LoadWrappers();
  
   std::string input;

   read->GetKey( varname , input );

   TFile * file = new TFile( input.c_str()  ); 
   SetLikelihoodFile( file ); 
}


// requires a DataManager to be defined
void pcthruLikelihood::LoadWrappers()
{
   dm->Get("nmue", nmue);
   dm->Get("evis", evis);
   dm->Get("ip", ip);
   dm->Get("nring", nring);
   dm->Get("amome", amome);
   dm->Get("amomm", amomm);
   dm->Get("nmue", nmue);
   dm->Get("evis", evis);
   dm->Get("etime", etime);
   dm->Get("etype", etype);
   dm->Get("ehit", ehit);
   dm->Get("egood", egood);
   dm->Get("pos", pos);
   dm->Get("dirtotmue", dirtotmue);
   dm->Get("etotmue", etotmue);

   dm->Get("msdir"  , msdir   );
   dm->Get("prmslg" , prmslg  ); 
   dm->Get("epos"   , epos    );
   dm->Get("probms" , probms  );

   dm->Get("dir"    , dir     );
   dm->Get("ipnu"   , ipnu    );
   dm->Get("mode"   , mode    ); 
   dm->Get("oscwgt" , oscwgt  );

  
}



float pcthruLikelihood::llBuild_pcthru( )
{
    float ll = 0.;


    int   mip    = mipBuild(); // must be called first!
       
    int   ndcye = muedcyBuild();  // muedcy is a class variable, and may be accessed externally   
    int   ringn  = ringnBuild(); 
    float transmom   = transmomBuild();
  
    
    int bin = GetEnergyBin( );

    // this is the actual likelihood value computation
    if( bin > 0 )
    {
//     ll += LoadLikelihood_pcthru( 340000, float(ndcye) , bin );
//     ll += LoadLikelihood_pcthru( 350000, float(ringn) , bin );
       
       //if( dpos > 1.0e-5 )
       //   ll += LoadLikelihood( 220000, sqrt(dpos)   , bin );
       if(ringn > 1)
       {
//         ll += LoadLikelihood_pcthru( 360000, transmom, bin);
       }

      /*if( kVerbose )
         std::cout << "mgmre " << ll
                   << " " <<  LoadLikelihood_numu( 210000, float(ndcye) , bin ) 
                   << " " <<  LoadLikelihood_numu( 220000, sqrt(dpos)    , bin )
                   << " " <<  LoadLikelihood_numu( 230000, ringn         , bin )
                   << " " <<  LoadLikelihood_numu( 240000, transmom      , bin )
                   << std::endl;*/

    }

   return ll;
}


void pcthruLikelihood::llBuild_pcthru_read()
{
     ofstream fout1("numu.txt", ios::app);
     ofstream fout2("numubar.txt", ios::app);
     ofstream fout3("nue.txt", ios::app);
     ofstream fout4("NC.txt", ios::app);
     
     int   mip    = mipBuild(); // must be called first!
       
     int   ndcye = muedcyBuild();  // muedcy is a class variable, and may be accessed externally
     int   ringn  = ringnBuild(); 
     float transmom   = transmomBuild();
     
     float ll_ndcye=0;
     float ll_ringn=0;
     float ll_transmom=0;
   
     
     int bin = GetEnergyBin( );
     
     
     if ( bin > 0)
     {
	
//     ll_ndcye = LoadLikelihood_pcthru( 340000, float(ndcye) , bin );
//     ll_ringn = LoadLikelihood_pcthru( 350000, float(ringn)  , bin );

       if(ringn > 1)
       {
//         ll_transmom = LoadLikelihood_pcthru( 360000, transmom      , bin );
       } 
       
     
     if(ll_ndcye != ll_ndcye || ll_ringn != ll_ringn || ll_transmom != ll_transmom)
     {
          cout << ndcye << ' ' << ll_ndcye << ' ' << ringn << ' ' << ll_ringn << ' ' << transmom << ' ' << ll_transmom << ' ' <<etotmue(0) << endl; 
     }
     
     if(ringn>1)
     {
         fout1 << ipnu(0) << ' ' << mode(0)  << ' ' << etotmue(0)  << ' ' << ll_ndcye + ll_ringn + ll_transmom << ' ' << oscwgt(0) << endl;
     }
     else
     {
         fout1 << ipnu(0) << ' ' << mode(0)  << ' ' << etotmue(0)  << ' ' << ll_ndcye + ll_ringn << ' ' << oscwgt(0) << endl;
     }
     
     /*
     if(ipnu(0)==14 && abs(mode(0))<30)
     {
          if( ringn > 1)
	  {
             fout1 << ll_ndcye  + ll_ringn + ll_transmom  << '\t' << ll_ndcye << '\t' <<  ll_ringn << '\t' << ll_transmom  <<'\t' << ndcye << '\t' <<  ringn << '\t' << transmom<< '\t' <<  etotmue(0) << '\t' << oscwgt(0) << endl;
	  }
	  else
	  {
             fout1 << ll_ndcye  + ll_ringn  << '\t' << ll_ndcye << '\t' <<  ll_ringn << '\t' << ll_transmom  <<'\t' << ndcye << '\t' <<  ringn << '\t' << transmom<< '\t' <<  etotmue(0) << '\t' << oscwgt(0) << endl;
	  
	  }
	 
     }
     else if (ipnu(0)==-14 && abs(mode(0))<30)
     {
          if(ringn > 1)
	  {
             fout2 << ll_ndcye  + ll_ringn + ll_transmom  << '\t' << ll_ndcye << '\t' <<  ll_ringn << '\t' << ll_transmom  <<'\t' << ndcye << '\t' <<  ringn << '\t' << transmom<< '\t' <<  etotmue(0) << '\t' << oscwgt(0) << endl;
	  }
	  else
	  {
             fout2 << ll_ndcye  + ll_ringn   << '\t' << ll_ndcye << '\t' <<  ll_ringn << '\t' << ll_transmom  <<'\t' << ndcye << '\t' <<  ringn << '\t' << transmom<< '\t' <<  etotmue(0) << '\t' << oscwgt(0) << endl;
	  
	  }
     }     
     else if (abs(ipnu(0)==12 && abs(mode(0))<30))
     {
          if(ringn > 1)
	  {
             fout3 << ll_ndcye  + ll_ringn + ll_transmom  << '\t' << ll_ndcye << '\t' <<  ll_ringn << '\t' << ll_transmom  <<'\t' << ndcye << '\t' <<  ringn << '\t' << transmom<< '\t' <<  etotmue(0) << '\t' << oscwgt(0) << endl;
	  }
	  else
	  {
             fout3 << ll_ndcye  + ll_ringn   << '\t' << ll_ndcye << '\t' <<  ll_ringn << '\t' << ll_transmom  <<'\t' << ndcye << '\t' <<  ringn << '\t' << transmom<< '\t' <<  etotmue(0) << '\t' << oscwgt(0) << endl;
	  
	  }
	      
     }
     else if (abs(mode(0))>30)
     {
          if(ringn > 1)
	  {
             fout4 << ll_ndcye  + ll_ringn + ll_transmom  << '\t' << ll_ndcye << '\t' <<  ll_ringn << '\t' << ll_transmom  <<'\t' << ndcye << '\t' <<  ringn << '\t' << transmom<< '\t' <<  etotmue(0) << '\t' << oscwgt(0) << endl;
	  }
	  else
	  {
             fout4 << ll_ndcye  + ll_ringn   << '\t' << ll_ndcye << '\t' <<  ll_ringn << '\t' << ll_transmom  <<'\t' << ndcye << '\t' <<  ringn << '\t' << transmom<< '\t' <<  etotmue(0) << '\t' << oscwgt(0) << endl;
	  
	  }  
	 
     }*/
     
     }
}




int pcthruLikelihood::mipBuild()
{

  int   i, lmip;
  
  emax = 0.;
  tot_energy = 0.;
  mering = 0;
  
  float p=0;

  lmip = ip(0);

  //selection taken from atmpd.kumac
  for( i = 0 ; i < nring(0) ; i ++ )
  { 
     p = ( prmslg(i,1) < prmslg(i,2) ? amome(i) : amomm(i) );
  
     tot_energy += amome(i);
     
     if ( p > emax )  
     {  
        mering = i;
        emax = p;
     }// end of mering selection
   }// end of ring loop

  // e-like
  if ( prmslg(mering,1) <  prmslg(mering,2) ) 
     lmip=2;
  else // mu-like
     lmip=3;

  return lmip;
}



//////////
//
//  This is the Hayakawa cut
//  criteria
////////////////
int pcthruLikelihood::muedcyBuild()
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
   
   muedcy=lmuedcy;

   return lmuedcy;

}



int pcthruLikelihood::ringnBuild()
{
     int lringn;
     
     lringn = nring(0);
     
     return lringn;
}


float pcthruLikelihood::transmomBuild()
{
     int i;
     
     emax=0;
     tot_energy=0;
     mering = 0;
     
     float p=0, pmax=0;
     
     for( i = 0; i < nring(0); i++)
     {
        p = ( prmslg(i,1) < prmslg(i,2) ? amome(i) : amomm(i));
	
	if (p > pmax) { pmax = p; mering = i;}
	
	tot_energy = tot_energy + p;     
     }
     
     
     float cos_angle[10], sin_angle[10];
     float ltransmom=0;
     
     for (int j=0; j< nring(0); j++)
     {
         p = ( prmslg(j,1) < prmslg(j,2) ? amome(j) : amomm(j));
         
         cos_angle[j] = dir(mering,0)*dir(j,0) + dir(mering,1)*dir(j,1) + dir(mering,2)*dir(j,2);
	 
	 if(cos_angle[j]*cos_angle[j] >=1)
	 {
	    sin_angle[j]=0;
	 }
	 else
	 {
	    sin_angle[j] = sqrt(1-cos_angle[j]*cos_angle[j]);
	 }
	 
	 ltransmom = ltransmom + p*sin_angle[j];
	 
     } 
     
     return ltransmom/tot_energy;

}


int pcthruLikelihood::GetEnergyBin( )
{

  int nbins = 6;
  float edges[] ={ 0.10, 1.330 , 2.5118, 5.0118, 10.00, 19.952, 1.0e6 };
  float E = etotmue(0)/1000.0;

  int bin = 0;

  for( int i = 1 ; i <= nbins ; i++ )
    if( E > edges[i-1] && E < edges[i] )
      bin = i;

  // adjust to match histogram numbering
  return bin;

}





