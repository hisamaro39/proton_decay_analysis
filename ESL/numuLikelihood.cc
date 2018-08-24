#include "ESL/numuLikelihood.h"
#include <fstream>
#include "TTree.h"

using namespace std;

numuLikelihood::numuLikelihood( DataManager * _dm , CardReader * reader , int skg )
:LikelihoodReturn( _dm, reader )
{
   skgen = skg;

   LoadWrappers();
  
   std::string input;

   read->GetKey("numuLikelihood", input );

   TFile * file = new TFile( input.c_str()  ); 
   
      SetLikelihoodFile( file ); 
}



numuLikelihood::numuLikelihood( DataManager * _dm , CardReader * reader , int skg, const char * varname )
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
void numuLikelihood::LoadWrappers()
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



float numuLikelihood::llBuild_numu( )
{
    float ll = 0.;


    int   mip    = mipBuild(); // must be called first!
       
    int   ndcye = muedcyBuild();  // muedcy is a class variable, and may be accessed externally
    int   ndcye_mu = ndcye_muBuild();
    int   ndcye_pi = ndcye_piBuild();
    float etime    = etimeBuild();

    float dpos   = dposBuild();
    int   ringn  = ringnBuild(); 
    float transmom   = transmomBuild();
  
    
    int bin = GetEnergyBin( );

    // this is the actual likelihood value computation
    if( bin > 0 )
    {
//     ll += LoadLikelihood_numu( 210000, float(ndcye_mu) , bin );
//     ll += LoadLikelihood_numu( 220000, float(ndcye_pi) , bin );
//     ll += LoadLikelihood_numu( 230000, ringn         , bin );
//     ll += LoadLikelihood_numu( 240000, transmom      , bin );
//     
       //if( dpos > 1.0e-5 )
       //   ll += LoadLikelihood( 220000, sqrt(dpos)   , bin );
//     if(etime > 1.0)
//     {
//         ll += LoadLikelihood_numu( 280000, etime, bin);
//     }

//    if( kVerbose )
//       std::cout << "mgmre " << ll
//                 << " " <<  LoadLikelihood_numu( 210000, float(ndcye) , bin ) 
//                 << " " <<  LoadLikelihood_numu( 220000, sqrt(dpos)    , bin )
//                 << " " <<  LoadLikelihood_numu( 230000, ringn         , bin )
//                 << " " <<  LoadLikelihood_numu( 240000, transmom      , bin )
//                 << std::endl;

    }

   return ll;
}


void numuLikelihood::llBuild_numu_read()
{
     ofstream fout1("numu.txt", ios::app);
     ofstream fout2("numubar.txt", ios::app);
     ofstream fout3("nue.txt", ios::app);
     ofstream fout4("NC.txt", ios::app);
     
     int   mip    = mipBuild(); // must be called first!
       
     int   ndcye = muedcyBuild();  // muedcy is a class variable, and may be accessed externally
     int   ndcye_mu = ndcye_muBuild();
     int   ndcye_pi = ndcye_piBuild();
     float etime    = etimeBuild();
    
     float dpos   = dposBuild();
     int   ringn  = ringnBuild(); 
     float transmom   = transmomBuild();
     
     float ll_ndcye_mu=0;
     float ll_ndcye_pi=0;
     float ll_ringn=0;
     float ll_transmom=0;
     float ll_dpos=0;
     float ll_etime=0;
     
     int bin = GetEnergyBin( );
     
     if ( bin > 0)
     {
//     ll_ndcye_mu = LoadLikelihood_numu( 210000, float(ndcye_mu) , bin );
//     ll_ndcye_pi = LoadLikelihood_numu( 220000, float(ndcye_pi) , bin );
//     ll_ringn  = LoadLikelihood_numu( 230000, ringn         , bin );
//     ll_transmom = LoadLikelihood_numu( 240000, transmom      , bin );

//     if(etime > 1.0)
//     {
//         ll_etime = LoadLikelihood_numu( 280000, etime, bin);
//     } 
     
     /*
     if(ll_ndcye_mu != ll_ndcye_mu || ll_ndcye_pi != ll_ndcye_pi || ll_ringn != ll_ringn || ll_transmom!= ll_transmom  || ll_etime != ll_etime)
     {
          cout << ndcye_mu << ' ' << ll_ndcye_mu << ' ' << etotmue(0) <<  endl;
	  cout << ndcye_pi << ' ' << ll_ndcye_pi << ' ' << endl;
	  cout << ringn <<  ' ' << ll_ringn <<  endl;
	  cout << transmom << ' ' << ll_transmom << endl;
	  cout << etime << ' ' << ll_etime <<  endl;
     }*/
     
     /*
     if(etime > 1.0)
     {
         fout1 << ipnu(0) << ' ' << mode(0) << ' ' << etotmue(0) << ' ' << ll_ndcye_mu + ll_ndcye_pi + ll_ringn + ll_transmom + ll_etime << '\t' << oscwgt(0) << endl;
     }
     else
     {
         fout1 <<  ipnu(0) << ' ' << mode(0) << ' ' << etotmue(0)<< ' ' <<  ll_ndcye_mu + ll_ndcye_pi + ll_ringn + ll_transmom << '\t' << oscwgt(0) << endl;
     }*/
     
     if(etime > 1.0 && !(ll_ndcye_mu != ll_ndcye_mu || ll_ndcye_pi != ll_ndcye_pi || ll_ringn != ll_ringn || ll_transmom != ll_transmom || ll_etime !=ll_etime))
     {
         fout1 << ll_ndcye_mu + ll_ndcye_pi + ll_ringn + ll_transmom + ll_etime << endl;
     }
     else if(etime <= 1.0 && !(ll_ndcye_mu != ll_ndcye_mu || ll_ndcye_pi != ll_ndcye_pi || ll_ringn != ll_ringn || ll_transmom != ll_transmom))
     {
         fout1 << ll_ndcye_mu + ll_ndcye_pi + ll_ringn + ll_transmom << endl;
     }
     else
     {
         fout1 << 999999 << endl;
     }
     
     /*
     if(ipnu(0)==14 && abs(mode(0))<30)
     {
          if(etime > 1.0)
	  {
             fout1 << ll_ndcye_mu + ll_ndcye_pi + ll_ringn + ll_transmom + ll_etime  << '\t' << ll_ndcye_mu << '\t' << ll_ndcye_pi << '\t' << ll_etime<<'\t' << ll_ringn << '\t' << ll_transmom  <<'\t' << ndcye_mu << '\t' << ndcye_pi << '\t' << etime << '\t' << ringn << '\t' << transmom<< '\t' <<  etotmue(0) << '\t' << oscwgt(0) << endl;
	  }
	  else
	  {
             fout1 << ll_ndcye_mu + ll_ndcye_pi + ll_ringn + ll_transmom   << '\t' << ll_ndcye_mu << '\t' << ll_ndcye_pi << '\t' << ll_etime<<'\t' << ll_ringn << '\t' << ll_transmom  <<'\t' << ndcye_mu << '\t' << ndcye_pi << '\t' << etime << '\t' << ringn << '\t' << transmom<< '\t' <<  etotmue(0) << '\t' << oscwgt(0) << endl;
	  
	  }
	 
     }
     else if (ipnu(0)==-14 && abs(mode(0))<30)
     {
          if(etime > 1.0)
	  {
            fout2 << ll_ndcye_mu + ll_ndcye_pi + ll_ringn + ll_transmom + ll_etime << '\t' << ll_ndcye_mu << '\t' << ll_ndcye_pi << '\t' << ll_etime <<'\t' << ll_ringn << '\t' << ll_transmom << '\t'<<  ndcye_mu << '\t' << ndcye_pi << '\t' << etime << '\t' << ringn << '\t' << transmom<< '\t' << etotmue(0) << '\t' << oscwgt(0) << endl;	 
	  }
	  else
	  {
             fout2 << ll_ndcye_mu + ll_ndcye_pi + ll_ringn + ll_transmom   << '\t' << ll_ndcye_mu << '\t' << ll_ndcye_pi << '\t' << ll_etime<<'\t' << ll_ringn << '\t' << ll_transmom  <<'\t' << ndcye_mu << '\t' << ndcye_pi << '\t' << etime << '\t' << ringn << '\t' << transmom<< '\t' <<  etotmue(0) << '\t' << oscwgt(0) << endl;
	  
	  }
     }     
     else if (abs(ipnu(0)==12 && abs(mode(0))<30))
     {
          if(etime > 1.0)
	  {
            fout3 << ll_ndcye_mu + ll_ndcye_pi + ll_ringn + ll_transmom + ll_etime << '\t' << ll_ndcye_mu << '\t' << ll_ndcye_pi << '\t' << ll_etime <<'\t' << ll_ringn << '\t' << ll_transmom  << '\t' << ndcye_mu << '\t' << ndcye_pi << '\t' << etime << '\t' << ringn << '\t' << transmom<< '\t' << etotmue(0) << '\t' << oscwgt(0) << endl;
	  }
	  else
	  {
             fout3 << ll_ndcye_mu + ll_ndcye_pi + ll_ringn + ll_transmom   << '\t' << ll_ndcye_mu << '\t' << ll_ndcye_pi << '\t' << ll_etime<<'\t' << ll_ringn << '\t' << ll_transmom  <<'\t' << ndcye_mu << '\t' << ndcye_pi << '\t' << etime << '\t' << ringn << '\t' << transmom<< '\t' <<  etotmue(0) << '\t' << oscwgt(0) << endl;
	  
	  }
	      
     }
     else if (abs(mode(0))>30)
     {
          if(etime > 1.0)
	  {
             fout4 << ll_ndcye_mu + ll_ndcye_pi + ll_ringn + ll_transmom + ll_etime << '\t' << ll_ndcye_mu << '\t' << ll_ndcye_pi << '\t' << ll_etime <<'\t' << ll_ringn << '\t' << ll_transmom << '\t'  << ndcye_mu << '\t' << ndcye_pi << '\t' << etime << '\t' << ringn << '\t' << transmom<< '\t' << etotmue(0) << '\t' << oscwgt(0) << endl;
	  }
	  else
	  {
             fout4 << ll_ndcye_mu + ll_ndcye_pi + ll_ringn + ll_transmom   << '\t' << ll_ndcye_mu << '\t' << ll_ndcye_pi << '\t' << ll_etime<<'\t' << ll_ringn << '\t' << ll_transmom  <<'\t' << ndcye_mu << '\t' << ndcye_pi << '\t' << etime << '\t' << ringn << '\t' << transmom<< '\t' <<  etotmue(0) << '\t' << oscwgt(0) << endl;
	  
	  }  
	 
     }*/
     
     }
}




int numuLikelihood::mipBuild()
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
int numuLikelihood::muedcyBuild()
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

float numuLikelihood::dposBuild()
{

//  return longest distance between primary vertex position
//  and decay electron position, divided by mering's energy
     
   float dpos;
   float dx;
   float x, y , z;
   int   nmuemax;

   int ehit_cut;


   if( skgen == 0 ) ehit_cut = 60;
   if( skgen == 1 ) ehit_cut = 30;
   if( skgen == 2 ) ehit_cut = 60;
   if( skgen == 3 ) ehit_cut = 60;


   nmuemax = ( nmue(0) > 10 ? 10 :  nmue(0) );
   if ( nmue(0) <= 0 ) return 0.;

   dpos=0.;
   for( int i = 0 ; i < nmuemax ; i++ )
   {
     if( skgen == 3) //for SK4 only
     {
         if( etime(i) < 0.1) continue;
	 if( etime(i) < 0.6) continue;
	 
	 if( etype(i) == 1 )
	 {
	    x = pos(0)-epos(i,0);
	    y = pos(1)-epos(i,1);
	    z = pos(2)-epos(i,2);
	    dx = sqrt(x*x + y*y + z*z);
	    
	    if( dx > dpos)
	    {
	        dpos = dx;
	    }
	 }
     }
     else
     {
   
     	 if( etime(i) <  0.1 ) continue;
     	 if( etime(i) >  0.8 && etime(i) < 1.2 ) continue;
     	 if( etype(i) == 1   && ehit(i) >= ehit_cut && egood(i) > 0.5 )
     	 {
         	x  = pos(0) - epos(i,0);
         	y  = pos(1) - epos(i,1);
         	z  = pos(2) - epos(i,2);
         	dx = sqrt( x*x + y*y + z*z );

		if( dx > dpos ) dpos = dx;
     	 }
     }
   }// end of loop on decay-e's

   return dpos / emax ;
}




int numuLikelihood::ndcye_muBuild()
{
     float rlist[36]={49.7651, 71.6595, 94.5605, 118.139, 142.227, 
                166.463, 190.205, 214.821, 238.699, 262.719,
                285.398, 310.097, 334.388, 357.353, 381.802,
                404.914, 427.625, 450.296, 474.211, 495.964,
                519.958, 542.554, 565.729, 587.981, 609.549,
                632.738, 654.881, 682.712, 694.223, 724.386,
                744.181, 769.437, 785.629, 812.260, 831.693,
                852.166}; 
		
     float pos_stop[3];
     float prange;
     float pmu, p1;
     int i1, i2;     
     
     int ehit_cut;
     if( skgen == 0 ) ehit_cut = 60;
     if( skgen == 1 ) ehit_cut = 30;
     if( skgen == 2 ) ehit_cut = 60;
     if( skgen == 3 ) ehit_cut = 60;
     
     pmu = amomm(mering);     		

     if(pmu <= 225)
     {
	i1=0;
	i2=1;
     }
     else if(pmu > 1975)
     {
	i1=34;
	i2=35;
     }
    else
    {
	i1 = int((pmu-225)/50);
	i2 = i1+1;
    }
    
    p1 =  225+(i1)*50;
    prange = rlist[i1]+(rlist[i2]-rlist[i1])/50*(pmu-p1);	
	      
    pos_stop[0]= prange*dir(mering,0)+pos(0);
    pos_stop[1]= prange*dir(mering,1)+pos(1);
    pos_stop[2]= prange*dir(mering,2)+pos(2);		
    
    float distance_tmp=0;
    int ndcye_mu=0;
    
    for(int i=0; i < muedcy; i++)
    {
    
        if(skgen == 3)
	{
           if(etime(i)<0.6) continue;
	   if(etype(i)==1)
	   {
	       distance_tmp = sqrt((epos(i,0)-pos_stop[0])*(epos(i,0)-pos_stop[0]) + (epos(i,1)-pos_stop[1])*(epos(i,1)-pos_stop[1]) +(epos(i,2)-pos_stop[2])*(epos(i,2)-pos_stop[2]));	     
	    
	       if(distance_tmp < 200)
	       {
	           ndcye_mu++;
	       }
	   }
	
        }
	else
	{
	    if( etime(i) < 0.1) continue;
	    if( etime(i) > 0.8 && etime(i) < 1.2) continue;
	    if( etype(i) == 1 && ehit(i) >= ehit_cut && egood(i) > 0.5)
	    {
	       distance_tmp = sqrt((epos(i,0)-pos_stop[0])*(epos(i,0)-pos_stop[0]) + (epos(i,1)-pos_stop[1])*(epos(i,1)-pos_stop[1]) +(epos(i,2)-pos_stop[2])*(epos(i,2)-pos_stop[2]));	     
	    
	       if(distance_tmp < 200)
	       {
	           ndcye_mu++;
	       }	        
	    } 
	}
    }
    
    return ndcye_mu;     
}


int numuLikelihood::ndcye_piBuild()
{
     float rlist[36]={49.7651, 71.6595, 94.5605, 118.139, 142.227, 
                166.463, 190.205, 214.821, 238.699, 262.719,
                285.398, 310.097, 334.388, 357.353, 381.802,
                404.914, 427.625, 450.296, 474.211, 495.964,
                519.958, 542.554, 565.729, 587.981, 609.549,
                632.738, 654.881, 682.712, 694.223, 724.386,
                744.181, 769.437, 785.629, 812.260, 831.693,
                852.166}; 
		
     float pos_stop[3];
     float prange;
     float pmu, p1;
     int i1, i2;     

     int ehit_cut;
     if( skgen == 0 ) ehit_cut = 60;
     if( skgen == 1 ) ehit_cut = 30;
     if( skgen == 2 ) ehit_cut = 60;
     if( skgen == 3 ) ehit_cut = 60;
     
     pmu = amomm(mering);     		

     if(pmu <= 225)
     {
	i1=0;
	i2=1;
     }
     else if(pmu > 1975)
     {
	i1=34;
	i2=35;
     }
    else
    {
	i1 = int((pmu-225)/50);
	i2 = i1+1;
    }
    
    p1 =  225+(i1)*50;
    prange = rlist[i1]+(rlist[i2]-rlist[i1])/50*(pmu-p1);	
	      
    pos_stop[0]= prange*dir(mering,0)+pos(0);
    pos_stop[1]= prange*dir(mering,1)+pos(1);
    pos_stop[2]= prange*dir(mering,2)+pos(2);		
    
    float distance_tmp=0;
    int ndcye_pi=0;
    
    for(int i=0; i < muedcy; i++)
    {
        if( skgen == 3)
	{
           if(etime(i)<0.6) continue;
	   if(etype(i)==1)
	   {
	       distance_tmp = sqrt((epos(i,0)-pos_stop[0])*(epos(i,0)-pos_stop[0]) +(epos(i,1)-pos_stop[1])*(epos(i,1)-pos_stop[1]) +(epos(i,2)-pos_stop[2])*(epos(i,2)-pos_stop[2]));	     
	    
	       if(distance_tmp >= 200)
	       {
	          ndcye_pi++;
	       }
	   }
	}
	else
	{
	    if( etime(i) < 0.1) continue;
	    if( etime(i) > 0.8 && etime(i) < 1.2) continue;
	    if( etype(i) == 1 && ehit(i) >= ehit_cut && egood(i) > 0.5)
	    {
	       distance_tmp = sqrt((epos(i,0)-pos_stop[0])*(epos(i,0)-pos_stop[0]) +(epos(i,1)-pos_stop[1])*(epos(i,1)-pos_stop[1]) +(epos(i,2)-pos_stop[2])*(epos(i,2)-pos_stop[2]));	     
	    
	       if(distance_tmp >= 200)
	       {
	          ndcye_pi++;
	       }	       
	    }
	}
	
	
    }
    
    return ndcye_pi;
    
}


float numuLikelihood::etimeBuild()
{
     float rlist[36]={49.7651, 71.6595, 94.5605, 118.139, 142.227, 
                166.463, 190.205, 214.821, 238.699, 262.719,
                285.398, 310.097, 334.388, 357.353, 381.802,
                404.914, 427.625, 450.296, 474.211, 495.964,
                519.958, 542.554, 565.729, 587.981, 609.549,
                632.738, 654.881, 682.712, 694.223, 724.386,
                744.181, 769.437, 785.629, 812.260, 831.693,
                852.166}; 
		
     float pos_stop[3];
     float prange;
     float pmu, p1;
     int i1, i2;     
     
     int ehit_cut;
     if( skgen == 0 ) ehit_cut = 60;
     if( skgen == 1 ) ehit_cut = 30;
     if( skgen == 2 ) ehit_cut = 60;
     if( skgen == 3 ) ehit_cut = 60;     
     
     pmu = amomm(mering);     		

     if(pmu <= 225)
     {
	i1=0;
	i2=1;
     }
     else if(pmu > 1975)
     {
	i1=34;
	i2=35;
     }
    else
    {
	i1 = int((pmu-225)/50);
	i2 = i1+1;
    }
    
    p1 =  225+(i1)*50;
    prange = rlist[i1]+(rlist[i2]-rlist[i1])/50*(pmu-p1);	
	    
    pos_stop[0]= prange*dir(mering,0)+pos(0);
    pos_stop[1]= prange*dir(mering,1)+pos(1);
    pos_stop[2]= prange*dir(mering,2)+pos(2);		
    
    float distance_tmp=0, distance_min = 100000;
    float etime_mu=0;
    
    for(int i=0; i < muedcy; i++)
    {
        if(skgen == 3)
	{
           if(etime(i)<0.6) continue;
	   if(etype(i)==1)
	   {
	       distance_tmp = sqrt((epos(i,0)-pos_stop[0])*(epos(i,0)-pos_stop[0]) +(epos(i,1)-pos_stop[1])*(epos(i,1)-pos_stop[1]) +(epos(i,2)-pos_stop[2])*(epos(i,2)-pos_stop[2]));	     
	    
	       if(distance_tmp < 200)
	       {
	           if(distance_tmp < distance_min)
		   {
		       distance_min = distance_tmp;
		       etime_mu = etime(i);		     
		   }
	       }
	   }
	}
	else
	{
	    if(etime(i) < 0.1) continue;
	    if(etime(i) > 0.8 && etime(i) < 1.2) continue;
	    if(etype(i) == 1 && ehit(i) >= ehit_cut && egood(i) > 0.5)
	    {
	       distance_tmp = sqrt((epos(i,0)-pos_stop[0])*(epos(i,0)-pos_stop[0]) +(epos(i,1)-pos_stop[1])*(epos(i,1)-pos_stop[1]) +(epos(i,2)-pos_stop[2])*(epos(i,2)-pos_stop[2]));	     
	    
	       if(distance_tmp < 200)
	       {
	           if(distance_tmp < distance_min)
		   {
		       distance_min = distance_tmp;
		       etime_mu = etime(i);		     
		   }
	       }	    
	    }
	}
	
    }
    
    return etime_mu;
    
}





int numuLikelihood::ringnBuild()
{
     int lringn;
     
     lringn = nring(0);
     
     return lringn;
}


float numuLikelihood::transmomBuild()
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


int numuLikelihood::GetEnergyBin( )
{

  int nbins = 6;
  float edges[] ={ 0.1, 1.00, 2.5118, 5.0118, 10.00, 19.952, 1.0e6 };
  float E = etotmue(0)/1000.0;

  int bin = 0;

  for( int i = 1 ; i <= nbins ; i++ )
    if( E > edges[i-1] && E < edges[i] )
      bin = i;

  // adjust to match histogram numbering
  return bin;

}





