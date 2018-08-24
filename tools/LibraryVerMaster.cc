#include "tools/LibraryVerMaster.h"

LibraryVerMaster * LibraryVerMaster::theMaster = 0;
std::map< std::string , LibraryVerMaster * > LibraryVerMaster::regLibs;

LibraryVerMaster::LibraryVerMaster()
{  
  Ver  = -1 ;
  Sub  = -1 ;
  Rev  = -1 ;
}

LibraryVerMaster* LibraryVerMaster::GetMaster()
{
  if ( LibraryVerMaster::theMaster == NULL ) 
  {
    LibraryVerMaster::theMaster = new LibraryVerMaster;
    return LibraryVerMaster::theMaster;
  }
  else 
       return LibraryVerMaster::theMaster;
}

void LibraryVerMaster::AddLibrary( LibraryVerMaster * ptr )
{
  std::string c = ptr->GetName();
  if ( ! regLibs.count( c ) ) 
  {
    std::cout << "LibraryVerMaster::AddLibrary " << c << " " << ptr << std::endl;
    //regLibs[ c ] = ptr;
  }
} 


void LibraryVerMaster::SaveVersionInfo( )
{
   TFile * f = new TFile("osc3pp.lib.ver.info.root", "recreate" );

   SaveVersionInfo( f );
   f->Close();

   delete f;
}

void LibraryVerMaster::SaveVersionInfo( TFile * f )
{

  if( f == 0 )
  {
      std::cout << "LibraryVerMaster::SaveVersionInfo requires an open TFile. " << std::endl; 
      std::cout << "  please pass this reoutine a non-null pointer"             << std::endl;
  }

  f->cd();
  
  TString lLabel = "unlabelled"; 
  int lVer ;
  int lSub ; 
  int lRev ;

  TTree * vertree = new TTree("osc3pp_ver_infor" , "Library Version Information for Osc3++" );

  vertree->Branch("libname"      , &lLabel       , 16000        , 0  );   
  vertree->Branch("ver"          , &lVer         , "ver/I"           );
  vertree->Branch("sub"          , &lSub         , "sub/I"           );
  vertree->Branch("rev"          , &lRev         , "rev/I"           );

  std::map< std::string , LibraryVerMaster * >::iterator _i;
  for( _i = regLibs.begin() ; _i != regLibs.end() ; _i++ )
  {	
    lLabel = _i->second->GetName(); 
    lVer   = _i->second->GetVersion();
    lSub   = _i->second->GetSubVersion();
    lRev   = _i->second->GetRevision();
    vertree->Fill();
  }

  vertree->Write();
}
