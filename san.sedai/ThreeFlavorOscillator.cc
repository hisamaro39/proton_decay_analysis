#include "san.sedai/ThreeFlavorOscillator.h"
#include "src/OscNtupleManager.h"
#include <math.h>

#include <iomanip>

#include "skcore/Fortran++.h"

#include "san.sedai/sansedaiLibVer.h"

// deprecated
/*ThreeFlavorOscillator::ThreeFlavorOscillator(bool kMantleHeikin ):Oscillator()
{

  //bNu = new ThreeFlavorPropagator( kMantleHeikin ); 	// true means with mantle averaging!
  bNu = new BargerPropagator();
  OscFlag = 3;
  kInverted = false;

  sansedaiLibVer * sslvp = sansedaiLibVer::GetLibrary();
  LibraryVerMaster::GetMaster()->AddLibrary( sslvp ); 
}

ThreeFlavorOscillator::ThreeFlavorOscillator( CardReader * lReader )
{
  bool kGoodRead = true;
  bool klRead;
  nPaths = 20;
   
  klRead = lReader->GetKey("OscFlag", OscFlag );
  if ( !klRead ) kGoodRead = false;

  klRead = lReader->GetKey("nPathLengths", nPaths );
  klRead = lReader->GetKey("InvertedHierarchy", kInverted );
 
  if ( !kGoodRead )
  { cerr << "Failure in CardReader at ThreeFlavorOscillator::ThreeFlavorOscillator" << endl; abort(); }

  bNu = new BargerPropagator();

  sansedaiLibVer * sslvp = sansedaiLibVer::GetLibrary();
  LibraryVerMaster::GetMaster()->AddLibrary( sslvp ); 
}


double ThreeFlavorOscillator::operator()( EventParser* _E, ProfileSpace& P , int n )
{
   return -1.0;
}*/	

//double ThreeFlavorOscillator::Oscillate( int Point, EventParser* _E , int n , ProfileSpace& P )
double ThreeFlavorOscillator::Oscillate( TTree *ltree )
{

      /*switch ( OscFlag ){
        case 1: return Neighbor3D        ( Point,  _E , n, P );

        case 4: return PathAve3D         ( Point,  _E , n, P );
        case 5: return Official2D        ( Point,  _E , n, P );
        case 6: return NoAve2D           ( Point,  _E , n, P );
        case 7: return CPT2D             ( Point,  _E , n, P );
        case 8: return CPT2DAngleFix     ( Point,  _E , n, P );
        case 9: return ENWmodel          ( Point,  _E , n, P ); 
	
        default: return 0.;
      }*/
  return 0.;

}


/*double ThreeFlavorOscillator::PathAve3D( int Point, EventParser* _E , int n , ProfileSpace& P )
{
	
  SKEventParser* E = static_cast<SKEventParser*>(_E);

  // Converts PDG->internal OscType    e:12->1 , mu: 14->2, tau:16->3 //
  int    NuOscillatedTo = int( abs( E->GetPDG() )/2 - 5 );
  double Oscillated = 0;	
  double MCWeight = E->GetMCWeight();
  double FactorE = 0, FactorMu = 0;
  int    NuStart;
  int    i,j;
  double Energy;

  // for the stupid fortran
  float  XXPATH[20];
  int    IDNU2 = NuOscillatedTo;
  float  DirNeu = (float ) E->GetNuCosineZ();                  
  float  E2 = (float ) E->GetEnergy();
  double Path;
  double OldPath = 0.0;
  
  // averaging counters
  double PathLengthAverage = 0.0;
  

  // Zero out any erroneou NC Tau events that 
  // may be included
  if ( abs( E->GetMode() ) >= 30 && abs(NuOscillatedTo)==3 )
    return 0;

  // NC events only recieve their solar flux weighting
  if ( abs( E->GetMode() ) >= 30 ) 
    return MCWeight;


  if ( E->GetPDG() ==  12 ) IDNU2 = 1;
  if ( E->GetPDG() == -12 ) IDNU2 = 3;
  if ( E->GetPDG() ==  14 ) IDNU2 = 2;
  if ( E->GetPDG() == -14 ) IDNU2 = 4;
  if ( E->GetPDG() ==  16 ) IDNU2 = 2;
  if ( E->GetPDG() == -16 ) IDNU2 = 4;


  Energy   = E->GetEnergy();
  FactorE  = ( NuOscillatedTo == 1 ?  1 : E->GetHondaFluxRatio( 2 ) );  // Use Muon Ratio  (e/mu) //
  FactorMu = ( NuOscillatedTo == 2 || NuOscillatedTo == 3 ? 1 : E->GetHondaFluxRatio( 1 ) ); // Use ERatio  (mu/e)///


  float SKheight = 0.380;
  Fortran::nebaseline_h3d_(XXPATH, IDNU2, DirNeu, E2, SKheight);
 
  bNu->SetMNS( P.Get("S12") , P.Get("S13") , P.Get("S23") , 
               P.Get("M12") , ( kInverted ? -1.0 : 1.0 ) * P.Get("M23") , P.Get("CP"), Energy , true ,  E->GetPDG() );  

  bNu->DefinePath( E->GetNuCosineZ(), 15.00, true );  // test for mantle condition
  bNu->SetMatterPathLength();

  int lnPaths = nPaths; 
  double PathLength;

  double probe = 0.;
  double probm = 0.;

  for( j = 0 ; j < 20 ; j++)
     OldPath += (double) XXPATH[j] /  20.0;     

  for (i = 0; i < lnPaths ; i++ ) // energy path length aver
  {
     PathLengthAverage += 1.00;                                                                                                                       
     // Subtracts the pathlength in matter above
     bNu->SetAirPathLength( XXPATH[j]  );	 

     NuStart = 1;
     NuStart *= ( E->GetPDG() < 0 ? -1 : 1 );

     bNu->propagate( NuStart );

     Oscillated += bNu->GetProb( NuStart, NuOscillatedTo )*MCWeight*FactorE;
     probe = bNu->GetProb( NuStart, NuOscillatedTo );
 
     NuStart = 2;
     NuStart *= ( E->GetPDG() < 0 ? -1 : 1 );

     Oscillated += bNu->GetProb( NuStart, NuOscillatedTo )*MCWeight*FactorMu;
     probm = bNu->GetProb( NuStart, NuOscillatedTo );

  } // End of Pathlength averaging


  Oscillated /= ( PathLengthAverage );
  
  return Oscillated;

}



double ThreeFlavorOscillator::Official2D( int Point, EventParser* _E , int n , ProfileSpace& P )
{
        // Converts PDG->internal OscType    e:12->1 , mu: 14->2, tau:16->3 //
        SKEventParser* E = static_cast<SKEventParser*>(_E);

        int    NuOscillatedTo = int( abs( E->GetPDG() )/2 - 5 );
        double Oscillated = 0;
        double Prob = 0.;
        double MCWeight = E->GetMCWeight();
        int    NuStart;
        double Energy = E->GetEnergy();

        // for the stupid fortran
        float XXPATH[20];
        float SKHeight = 0.380;
        int   IDNU2;

        float DirNeu = (float ) E->GetNuCosineZ();
        float E2     = (float ) E->GetEnergy();

        double Path;
        double DM2  = P.Get("M23");
        double Sin2 = P.Get("S23");

        // for phase averaging
        double PhaseThreshold = 2.0;
        double PathCounter    = 0.0;
        double Phase  = ( 3.1415 * Energy)/( 1.267 *DM2 ); // oscillation phase


        //std::cout << Sin2 << ": Sin2 "  << std::endl;
        int i;

       // there should not be any NC tau events ...
       if ( abs( E->GetMode() ) >= 30 && abs(NuOscillatedTo)==3 ) return 0;


// 1-electron neutrino
// 2-muon neutrino
// 3-anti electron neutrino
// 4-anti muon neutrino

        if ( E->GetPDG() ==  12 ) IDNU2 = 1;
        if ( E->GetPDG() == -12 ) IDNU2 = 3;
        if ( E->GetPDG() ==  14 ) IDNU2 = 2;
        if ( E->GetPDG() == -14 ) IDNU2 = 4;
        if ( E->GetPDG() ==  16 ) IDNU2 = 2;
        if ( E->GetPDG() == -16 ) IDNU2 = 4;

        // NC and electron neutrinos are unoscillated
        if ( abs( E->GetMode() ) >= 30 || NuOscillatedTo == 1 ) return MCWeight;


        Fortran::nebaseline_h3d_(XXPATH,IDNU2, DirNeu ,E2,SKHeight);

        Oscillated = 0.0;
        if ( abs( E->GetMode() ) <= 30 && ( NuOscillatedTo == 2 || NuOscillatedTo == 3 ) )
        {
           
           Prob = 0.;
           for (i = 0; i < 20 ; i++ )
           {
              Path = (double) XXPATH[i];  // in km
              if ( Path > 1.0  ){
                 if ( Path /Phase  <= PhaseThreshold  )
                    Prob +=  Sin2 * sin( 1.267 * DM2*Path/Energy ) *  sin( 1.267 * DM2*Path/Energy );
                 if ( Path /Phase  > PhaseThreshold )
                    Prob += (1.0 * 0.5)*Sin2; // average out ...
                 PathCounter += 1.0;
              }
           }//end of loop over paths

           if ( PathCounter > 0.0 ) Prob /= PathCounter; 

           if ( NuOscillatedTo == 2 ) Oscillated  = ( 1. - Prob )* MCWeight; // n_mu survival prob
           if ( NuOscillatedTo == 3 ) Oscillated  =  Prob * MCWeight;        // n_tau disappearance prob
        }// end of CC and nu_mu || nu_tau

	return Oscillated;

}



double ThreeFlavorOscillator::NoAve2D( int Point, EventParser* _E , int n , ProfileSpace& P )
{
        // Converts PDG->internal OscType    e:12->1 , mu: 14->2, tau:16->3 //
        SKEventParser* E = static_cast<SKEventParser*>(_E);

        int    NuOscillatedTo = int( abs( E->GetPDG() )/2 - 5 );
        double Oscillated = 0;
        double Prob = 0.;
        double MCWeight = E->GetMCWeight();
        int    NuStart;
        double Energy = E->GetEnergy();

        // for the stupid fortran
        float XXPATH[20];
        float SKHeight = 0.380;
        int   IDNU2;

        float DirNeu = (float ) E->GetNuCosineZ();
        float E2     = (float ) E->GetEnergy();

        double Path;
        double DM2  = P.Get("M23");
        double Sin2 = P.Get("S23");

        // for phase averaging
        double PhaseThreshold = 2.0;
        double PathCounter    = 0.0;
        double Phase  = ( 3.1415 * Energy)/( 1.267 *DM2 ); // oscillation phase


        int i;

       // there should not be any NC tau events ...
       if ( abs( E->GetMode() ) >= 30 && abs(NuOscillatedTo)==3 ) return 0;

        // NC and electron neutrinos are unoscillated
        if ( abs( E->GetMode() ) >= 30 || NuOscillatedTo == 1 ) return MCWeight;

// 1-electron neutrino
// 2-muon neutrino
// 3-anti electron neutrino
// 4-anti muon neutrino

        if ( E->GetPDG() ==  12 ) IDNU2 = 1;
        if ( E->GetPDG() == -12 ) IDNU2 = 3;
        if ( E->GetPDG() ==  14 ) IDNU2 = 2;
        if ( E->GetPDG() == -14 ) IDNU2 = 4;
        if ( E->GetPDG() ==  16 ) IDNU2 = 2;
        if ( E->GetPDG() == -16 ) IDNU2 = 4;



        Fortran::nebaseline_h3d_(XXPATH,IDNU2, DirNeu ,E2,SKHeight);

        Oscillated = 0.0;
        if ( abs( E->GetMode() ) <= 30 && ( NuOscillatedTo == 2 || NuOscillatedTo == 3 ) )
        {
           
           Prob = 0.;
           for (i = 0; i < 20 ; i++ )
           {
              Path = (double) XXPATH[i];  // in km
              if ( Path > 1.0  )
              {
                 Prob +=  Sin2 * sin( 1.267 * DM2*Path/Energy ) *  sin( 1.267 * DM2*Path/Energy );
                 PathCounter += 1.0;
              }
           }//end of loop over paths

           if ( PathCounter > 0.0 ) Prob /= PathCounter; 

           if ( NuOscillatedTo == 2 ) Oscillated  = ( 1. - Prob )* MCWeight; // n_mu survival prob
           if ( NuOscillatedTo == 3 ) Oscillated  =  Prob * MCWeight;        // n_tau disappearance prob
        }// end of CC and nu_mu || nu_tau


	return Oscillated;

}



double ThreeFlavorOscillator::CPT2D( int Point, EventParser* _E , int n , ProfileSpace& P )
{
        // Converts PDG->internal OscType    e:12->1 , mu: 14->2, tau:16->3 //
        SKEventParser* E = static_cast<SKEventParser*>(_E);

        int    NuOscillatedTo = int( abs( E->GetPDG() )/2 - 5 );
        double Oscillated = 0;
        double Prob = 0.;
        double MCWeight = E->GetMCWeight();
        int    NuStart;
        double Energy = E->GetEnergy();

        // for the stupid fortran
        float XXPATH[20];
        float SKHeight = 0.380;
        int   IDNU2;

        float DirNeu = (float ) E->GetNuCosineZ();
        float E2     = (float ) E->GetEnergy();

        double Path;
        double DM2  = P.Get("M23");
        double Sin2 = P.Get("S23");

        double DM2a  = P.Get("M23a");
        double Sin2a = P.Get("S23a");

        // for phase averaging
        double PhaseThreshold = 2.0;
        double PathCounter    = 0.0;
        double Phase     = ( 3.1415 * Energy)/( 1.267 *DM2  ); // oscillation phase
        double PhaseBar  = ( 3.1415 * Energy)/( 1.267 *DM2a ); // oscillation phase


        //std::cout << Sin2 << ": Sin2 "  << std::endl;
        int i;

       // there should not be any NC tau events ...
       if ( abs( E->GetMode() ) >= 30 && abs(NuOscillatedTo)==3 ) return 0;

        // NC and electron neutrinos are unoscillated
        if ( abs( E->GetMode() ) >= 30 || NuOscillatedTo == 1 ) return MCWeight;

// 1-electron neutrino
// 2-muon neutrino
// 3-anti electron neutrino
// 4-anti muon neutrino

        if ( E->GetPDG() ==  12 ) IDNU2 = 1;
        if ( E->GetPDG() == -12 ) IDNU2 = 3;
        if ( E->GetPDG() ==  14 ) IDNU2 = 2;
        if ( E->GetPDG() == -14 ) IDNU2 = 4;
        if ( E->GetPDG() ==  16 ) IDNU2 = 2;
        if ( E->GetPDG() == -16 ) IDNU2 = 4;



        Fortran::nebaseline_h3d_(XXPATH,IDNU2, DirNeu ,E2,SKHeight);

        Oscillated = 0.0;
        if ( abs( E->GetMode() ) <= 30 && ( NuOscillatedTo == 2 || NuOscillatedTo == 3 ) )
        {
           
           Prob = 0.;
           for (i = 0; i < 20 ; i++ )
           {
              Path = (double) XXPATH[i];  // in km
              if ( Path > 1.0  )
              {
                 // neutrino mixing
                 if ( Path /Phase  <= PhaseThreshold  )
                 if( E->GetPDG()> 0 ) 
                     Prob +=  Sin2 * sin( 1.267 * DM2*Path/Energy ) *  sin( 1.267 * DM2*Path/Energy );

                 // neutrino mixing
                 if ( Path /Phase  > PhaseThreshold )
                 if( E->GetPDG() > 0 ) 
                     Prob += (1.0 * 0.5)*Sin2; // average out ...

                 // anti-neutrino mixing
                 if ( Path /PhaseBar  <= PhaseThreshold  )
                 if( E->GetPDG()< 0 ) 
                     Prob +=  Sin2a * sin( 1.267 * DM2a*Path/Energy ) *  sin( 1.267 * DM2a*Path/Energy );

                 // anti-neutrino mixing
                 if ( Path /PhaseBar  > PhaseThreshold )
                 if( E->GetPDG() < 0 ) 
                     Prob += (1.0 * 0.5)*Sin2a; // average out ...

                 PathCounter += 1.0;
              }
           }//end of loop over paths

           if ( PathCounter > 0.0 ) Prob /= PathCounter; 

           if ( NuOscillatedTo == 2 ) Oscillated  = ( 1. - Prob )* MCWeight; // n_mu survival prob
           if ( NuOscillatedTo == 3 ) Oscillated  =  Prob * MCWeight;        // n_tau disappearance prob
        }// end of CC and nu_mu || nu_tau


	return Oscillated;

}



double ThreeFlavorOscillator::CPT2DAngleFix( int Point, EventParser* _E , int n , ProfileSpace& P )
{
   // Converts PDG->internal OscType    e:12->1 , mu: 14->2, tau:16->3 //
   SKEventParser* E = static_cast<SKEventParser*>(_E);

   int    NuOscillatedTo = int( abs( E->GetPDG() )/2 - 5 );
   double Oscillated = 0;
   double Prob = 0.;
   double MCWeight = E->GetMCWeight();
   int    NuStart;
   double Energy = E->GetEnergy();

   // for the stupid fortran
   float XXPATH[20];
   float SKHeight = 0.380;
   int   IDNU2;

   float DirNeu = (float ) E->GetNuCosineZ();
   float E2     = (float ) E->GetEnergy();

   double Path;
   double DM2  = P.Get("M23");
   double Sin2 = P.Get("S23");

   double DM2a  = P.Get("M23a");
   double Sin2a = P.Get("S23a");

//////
///     This is To duplicate Wei's fixed mixing angle analyses
///
   Sin2 = Sin2a;


   // for phase averaging
   double PhaseThreshold = 2.0;
   double PathCounter    = 0.0;
   double Phase     = ( 3.1415 * Energy)/( 1.267 *DM2  ); // oscillation phase
   double PhaseBar  = ( 3.1415 * Energy)/( 1.267 *DM2a ); // oscillation phase


   //std::cout << Sin2 << ": Sin2 "  << std::endl;
   int i;

   // there should not be any NC tau events ...
   if ( abs( E->GetMode() ) >= 30 && abs(NuOscillatedTo)==3 ) return 0;


// 1-electron neutrino
// 2-muon neutrino
// 3-anti electron neutrino
// 4-anti muon neutrino

   if ( E->GetPDG() ==  12 ) IDNU2 = 1;
   if ( E->GetPDG() == -12 ) IDNU2 = 3;
   if ( E->GetPDG() ==  14 ) IDNU2 = 2;
   if ( E->GetPDG() == -14 ) IDNU2 = 4;
   if ( E->GetPDG() ==  16 ) IDNU2 = 2;
   if ( E->GetPDG() == -16 ) IDNU2 = 4;

   // NC and electron neutrinos are unoscillated
   if ( abs( E->GetMode() ) >= 30 || NuOscillatedTo == 1 ) return MCWeight;


   Fortran::nebaseline_h3d_(XXPATH,IDNU2, DirNeu ,E2,SKHeight);

   Oscillated = 0.0;
   if ( abs( E->GetMode() ) <= 30 && ( NuOscillatedTo == 2 || NuOscillatedTo == 3 ) )
   {
           
      Prob = 0.;
      for (i = 0; i < 20 ; i++ )
      {
         Path = (double) XXPATH[i];  // in km
         if ( Path > 1.0  )
         {
            // neutrino mixing
            if ( Path /Phase  <= PhaseThreshold  )
            if( E->GetPDG()> 0 ) 
               Prob +=  Sin2 * sin( 1.267 * DM2*Path/Energy ) *  sin( 1.267 * DM2*Path/Energy );

             // neutrino mixing
             if ( Path /Phase  > PhaseThreshold )
             if( E->GetPDG() > 0 ) 
               Prob += (1.0 * 0.5)*Sin2; // average out ...

              // anti-neutrino mixing
              if ( Path /PhaseBar  <= PhaseThreshold  )
              if( E->GetPDG()< 0 ) 
                 Prob +=  Sin2a * sin( 1.267 * DM2a*Path/Energy ) *  sin( 1.267 * DM2a*Path/Energy );

              // anti-neutrino mixing
              if ( Path /PhaseBar  > PhaseThreshold )
              if( E->GetPDG() < 0 ) 
                 Prob += (1.0 * 0.5)*Sin2a; // average out ...

              PathCounter += 1.0;
           }
        }//end of loop over paths

        if ( PathCounter > 0.0 ) Prob /= PathCounter; 

        if ( NuOscillatedTo == 2 ) Oscillated  = ( 1. - Prob )* MCWeight; // n_mu survival prob
        if ( NuOscillatedTo == 3 ) Oscillated  =  Prob * MCWeight;        // n_tau disappearance prob
     }// end of CC and nu_mu || nu_tau


   return Oscillated;

}



double ThreeFlavorOscillator::ENWmodel( int Point, EventParser* _E , int n , ProfileSpace& P )
{
        // Converts PDG->internal OscType    e:12->1 , mu: 14->2, tau:16->3 //
        SKEventParser* E = static_cast<SKEventParser*>(_E);
	//	LineIntegral ENWPotential;

	int Flavor3 = 0;
        int    NuOscillatedTo = int( abs(E->GetPDG())/2 - 5 );
        double Oscillated = 0;
        double Prob = 0.;
	double Probtau = 0.;
	double PathCounter = 0.;
        double MCWeight = E->GetMCWeight();
        int    NuStart;
        double Energy = E->GetEnergy();

	// ENW Variables
	double theta_34, alpha, theta_s;
	double mtwiddlesq, Mtwiddlesq;
	double Rho, V, MatterPath;
	double EarthRadius = 6371.0;
	double UnitConversion = 1.62187868;

        // for the stupid fortran
        float XXPATH[20];
        float SKHeight = 0.380;
        int   IDNU2;

        float DirNeu = (float ) E->GetNuCosineZ();
        float E2     = (float ) E->GetEnergy();

        double Path;
        double m   = P.Get("m");
        double M   = P.Get("MS");
	double gmV = P.Get("gmV");

	theta_34 = (atan(2*sqrt(M*m)/(M-m)))/2.;

	MatterPath = 0.;
	if (DirNeu < 0.)
	  {
	    MatterPath = sqrt( EarthRadius*EarthRadius - EarthRadius*EarthRadius*(1 - DirNeu * DirNeu))
	      - EarthRadius*DirNeu;
	  }

        int i;

       // there should not be any NC tau events ...
       if ( abs( E->GetMode() ) >= 30 && abs(NuOscillatedTo) == 3 ) return 0;

        // NC and electron neutrinos are unoscillated
        if ( abs( E->GetMode() ) >= 30 || abs(NuOscillatedTo) == 1 ) return MCWeight;


// 1-electron neutrino
// 2-muon neutrino
// 3-anti electron neutrino
// 4-anti muon neutrino

        if ( E->GetPDG() ==  12 ) IDNU2 = 1;
        if ( E->GetPDG() == -12 ) IDNU2 = 3;
        if ( E->GetPDG() ==  14 ) IDNU2 = 2;
        if ( E->GetPDG() == -14 ) IDNU2 = 4;
        if ( E->GetPDG() ==  16 ) IDNU2 = 2;
        if ( E->GetPDG() == -16 ) IDNU2 = 4;

        Fortran::nebaseline_h3d_(XXPATH,IDNU2, DirNeu ,E2,SKHeight);

        Oscillated = 0.0;

        if ( abs( E->GetMode() ) <= 30 && ( abs( NuOscillatedTo ) == 2 || abs( NuOscillatedTo ) == 3 ))
        {
           
           Prob = 0.;
	   Probtau = 0.;
	   PathCounter = 0.;

           for (i = 0; i < 20 ; i++ )
           {
              Path = (double) XXPATH[i];  // in km
              if ( Path > 1.0  )
              {
		
		// Calculate ENW Potential for Zenith Angle and Path Length
		// Units: [m],[M] = eV; [mtwiddlesq],[Mtwiddlesq] = eV^2
		// Units: [V] = neV
		// Units: [E2] = GeV 
		// Units: [L] = km

		if (E->GetPDG() > 0)
		  {
//    Rho = ENWPotential->SampleRhoHistogram("/net/10.10.10.218/work24/dziomba/enw_model/Osc3++/san.sedai/PotentialvsZenithAngle820.root",DirNeu,2);
		  }

		if (E->GetPDG() < 0)
		  {
//    Rho = ENWPotential->SampleRhoHistogram("/net/10.10.10.218/work24/dziomba/enw_model/Osc3++/san.sedai/PotentialvsZenithAngle820.root",DirNeu,-2);
		  } 


		V = gmV*UnitConversion*Rho*fabs(MatterPath/Path);

		// Calculates equations 15, 20, 21, 23, 26 from ENW paper

		alpha = 4.0*E2*V/(M*M-m*m);

		theta_s = (atan(sin(2*theta_34)/(cos(2*theta_34)+alpha)))/2.;

		mtwiddlesq = m*m + ((M*M - m*m)/2) * (1 + alpha - sqrt(1 + 2*alpha*cos(2*theta_34) + alpha*alpha));
		Mtwiddlesq = m*m + ((M*M - m*m)/2) * (1 + alpha + sqrt(1 + 2*alpha*cos(2*theta_34) + alpha*alpha));

		Prob += 1 - cos(theta_s)*cos(theta_s)*sin(1.27*mtwiddlesq*Path/E2)*sin(1.27*mtwiddlesq*Path/E2)
                  - sin(theta_s)*sin(theta_s)*sin(1.27*Mtwiddlesq*Path/E2)*sin(1.27*Mtwiddlesq*Path/E2)
                  - cos(theta_s)*cos(theta_s)*sin(theta_s)*sin(theta_s)*sin(1.27*(Mtwiddlesq - mtwiddlesq)*Path/E2)*sin(1.27*(Mtwiddlesq 
															     - mtwiddlesq)*Path/E2);
		Probtau += cos(theta_s)*cos(theta_s)*sin(1.27*mtwiddlesq*Path/E2)*sin(1.27*mtwiddlesq*Path/E2)
		  + sin(theta_s)*sin(theta_s)*sin(1.27*Mtwiddlesq*Path/E2)*sin(1.27*Mtwiddlesq*Path/E2) 
		  - cos(theta_s)*cos(theta_s)*sin(theta_s)*sin(theta_s)*sin(1.27*(Mtwiddlesq - mtwiddlesq)*Path/E2)*sin(1.27*(Mtwiddlesq
															      - mtwiddlesq)*Path/E2);

                 PathCounter += 1.0;
              }
           }//end of loop over paths

           if ( PathCounter > 0.0 ) 
	     {
	       Prob /= PathCounter; 
	       Probtau /= PathCounter;
	     }

           if ( NuOscillatedTo == 2 || NuOscillatedTo == -2 ) 
	   { 
	     Oscillated  =  Prob * MCWeight; // n_mu survival prob
	   }

	   if ( abs(NuOscillatedTo) == 3)
	   {
	     Oscillated =  Probtau * MCWeight; // n_tau disapperance prob 
	   }
        }// end of CC and nu_mu and nu_tau

	return Oscillated;
}


double ThreeFlavorOscillator::Neighbor3D( int Point, EventParser* _E , int n , ProfileSpace& P )
{

  SKEventParser* E = static_cast<SKEventParser*>(_E);

  // Converts PDG->internal OscType    e:12->1 , mu: 14->2, tau:16->3 //
  int    NuOscillatedTo = int( abs( E->GetPDG() )/2 - 5 );
  double Oscillated = 0;	
  double MCWeight = E->GetMCWeight();
  double FactorE = 0, FactorMu = 0;
  int    NuStart;
  int    i,j;
  double Energy;
  bool   kUseSquaredThetas = true ;
   
  // Mass Hierarchy Flag 
  // ( correction for solar term is done internally )
  double hFactor           = ( kInverted ? -1.0 : 1.0 ) ; 

  // for the fortran computation
  float  XXPATH[20];
  int    IDNU2    = NuOscillatedTo;
  float  DirNeu   = (float ) E->GetNuCosineZ();                  
  float  fortranE = (float ) E->GetEnergy();
  double Path;
  
  // averaging counters
  double PathLengthAverage = 0.0;
  double FullPathAve = 0.0;
  double rms = E->GetEAveRMS() ;


  // Zero out any erroneou NC Tau events that 
  // may be included
  if ( abs( E->GetMode() ) >= 30 && abs(NuOscillatedTo)==3 )
    return 0;

  // NC events only recieve their solar flux weighting
  if ( abs( E->GetMode() ) >= 30 ) 
  {
    return MCWeight;
  }


  if ( E->GetPDG() ==  12 ) IDNU2 = 1;
  if ( E->GetPDG() == -12 ) IDNU2 = 3;
  if ( E->GetPDG() ==  14 ) IDNU2 = 2;
  if ( E->GetPDG() == -14 ) IDNU2 = 4;
  if ( E->GetPDG() ==  16 ) IDNU2 = 2;
  if ( E->GetPDG() == -16 ) IDNU2 = 4;

  Energy   = E->GetEnergy();
  if( rms / Energy  >= 0.5 ) rms = 0.5 * Energy ;

  // Flux factor adjustment : 
  // if nue      ,  FactorE   = 1 , otherwise reweight from muonflux to electron flux
  // if numu/tau ,  FactorMu  = 1 , otherwise reweight from electron to muon flux 
  FactorE  = ( NuOscillatedTo == 1 ?  1 : E->GetHondaFluxRatio( 2 ) );  
  FactorMu = ( NuOscillatedTo == 2 || NuOscillatedTo == 3 ? 1 : E->GetHondaFluxRatio( 1 ) );


  int NEAve = E->GetNEAve();


  double EAverages [5] =  { Energy           ,
                            Energy - rms     , 
                            Energy - rms *0.5 , 
                            Energy + rms *0.5 , 
                            Energy + rms      
                          }; 


  if( NEAve == 0 )
  {
     EAverages[1] = 1.50 * Energy ;
     EAverages[2] = 1.25 * Energy ;
     EAverages[3] = 0.75 * Energy ;
     EAverages[4] = 0.50 * Energy ;
  }



  float SKheight = 0.380;
  Fortran::nebaseline_h3d_(XXPATH, IDNU2, DirNeu, fortranE, SKheight);

  for( j = 0 ; j < 20 ; j++)
     FullPathAve += (double) XXPATH[j] /  20.0;     

  double AvePath   [5] ;
  AvePath [0]  =   FullPathAve * 1.0e5;  // in [cm]
 
  // These are lifted from osc3d_map.F (arbitrary assignment)
  int PathAveIndices[4][4] = {  // Original Fortran Indexing 
                                10 , 9 , 11 , 12  , 
                                 8 , 7 , 13 , 14  ,
                                 6 , 5 , 15 , 16  ,
                                 4 , 3 , 17 , 18 
                             };

  for( i = 0 ; i < 4 ; i++ ) 
  {
    AvePath [i+1]  = 0.0 ;
    for( j = 0 ; j < 4 ; j++ ) 
    {
       // convert to c-indexing
       int index = PathAveIndices[i][j] - 1;
       AvePath [i+1]  +=   (double) XXPATH[ index ] / 4.0   ; 
    }
    AvePath [i+1] *= 1.0e5  ; // convert to cm
  }
  
  
  // Set once with arbitrary production height,
  // will be overridden with call to SetPathLength below
  // false argument stops matter profile specification 
  bNu->DefinePath( E->GetNuCosineZ(), 15.00, false ); 

  bNu->SetMNS( P.Get("S12") , P.Get("S13") , P.Get("S23") , 
               P.Get("M12") ,      hFactor * P.Get("M23") , 
               P.Get("CP")  , 
               Energy       , kUseSquaredThetas ,  E->GetPDG() );  

  int lnPaths = 5; 
  int Type = E->GetType() % global::nEventTypes;
  
  // skip energy averaging for 
  // upthrough types
  if( Type == global::UpThruNonShower_mu ||
      Type == global::UpThruShower_mu  )  
       lnPaths = 1 ;
  

  double PathLength;
  double probe = 0.;
  double probm = 0.;

  for (i = 0; i < lnPaths ; i++ ) // energy path length aver
  {
     PathLengthAverage += 1.00;                                                                                                                       
     bNu->SetPathLength( AvePath[i]    );	 
     bNu->SetEnergy    ( EAverages[i]  );

     NuStart = 1;
     NuStart *= ( E->GetPDG() < 0 ? -1 : 1 );

     // Only need to propage once 
     bNu->propagate( NuStart );

     Oscillated += bNu->GetProb( NuStart, NuOscillatedTo )*MCWeight*FactorE;
     probe       = bNu->GetProb( NuStart, NuOscillatedTo );

     NuStart = 2;
     NuStart *= ( E->GetPDG() < 0 ? -1 : 1 );

     Oscillated += bNu->GetProb( NuStart, NuOscillatedTo )*MCWeight*FactorMu;
     probm       = bNu->GetProb( NuStart, NuOscillatedTo );

  } // End of Pathlength averaging

  Oscillated /= ( PathLengthAverage );

  return Oscillated;

}*/
