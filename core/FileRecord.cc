#include "core/FileRecord.h"
#include <iostream>

FileRecord::FileRecord( )
{
   nFiles   = 0 ;
   fMode    = 0 ;
   Instance = 0 ;
   nPoints  = 0 ;
   kUseNew=false;

}

FileRecord::FileRecord( const char * s, int i , int npoints ,int nfiles )
{

   FilePath = s; 
   Instance = i;
   nPoints  = npoints;
   nFiles   = nfiles;
   fMode = 0;
   kUseNew=false;
   FileInstance = BuildFileInstanceString( FilePath, Instance );

}



TString FileRecord::BuildFileInstanceString( const char * s , int instance )
{
   TString InstanceString = s;

   if(! InstanceString.Contains("*") )
      return InstanceString;

   std::stringstream ss; 

   char instance_char[50];
 
   if( fMode == 0 )
      sprintf( instance_char, "%06d" , instance );
   else 
      sprintf( instance_char, "%d" , instance );


   int pos;

   pos = InstanceString.Length()-7;  // take off .*.

   if( nFiles > 1 )
   {
      ss.str(""); ss << "."<< instance_char << ".root";
      InstanceString.Replace( pos, ss.str().length() , ss.str().c_str() );
   }
   else 
   {
      InstanceString.Remove( pos , 2 );
   }

   return InstanceString;

}


TString FileRecord::BuildOutputInstanceString( const char * s , int instance,  int nfiles  )
{

   TString base  = s ;
   std::stringstream ss;

   TString InstanceString;
   int pos ;

   char instance_char[50];
   if( fMode == 0 )
      sprintf( instance_char, "%06d" , instance );
   else 
      sprintf( instance_char, "%d" , instance );


   if( nfiles > 1 ) 
    if( base.EndsWith(".x") )
    {
       ss.str(""); ss << "."<< instance_char << ".x";
       InstanceString = ss.str(); 
       pos = base.Length()-2;

       base.Replace( pos, InstanceString.Length(), InstanceString );
    }
    else
    {
       ss.str(""); ss << "."<< instance_char << ".root";
       InstanceString = ss.str(); 
       pos = base.Length()-5;

       base.Replace( pos, InstanceString.Length(), InstanceString );
    }

   return base;
}


int FileRecord::GetInstanceFromPoint(  int Point , int npoints, int nfiles  )
{

   // 'l' stands for a local variable 
   // not the class one
   int lnPoints = ( npoints == -1 ? nPoints : npoints ); 
   int lnFiles  = ( nfiles  == -1 ? nFiles  : nfiles  ); 

   int lInstance; 
   int Seed; 
   int Blossom;

   // blech, brute force
   for( unsigned i = 0 ; i < lnFiles ; i++ ) 
   {
      lInstance = i;
      GetSeedAndBlossom( Seed, Blossom, lInstance , lnPoints, lnFiles );

      //      std::cout << " Seed: " << Seed << " " << lInstance << " " << Point << " " << Blossom << " " << lnPoints << " " << lnFiles << std::endl;

      if( Seed <= Point && Point < Blossom )
         return lInstance;
   } 

   return -1;
}



int FileRecord::GetPointOffsetFromInstance( int Point , int instance, int npoints, int nfiles)
{

   // 'l' stands for a local variable 
   // not the class one
   int lnPoints  = ( npoints   == -1 ? nPoints  : npoints   ); 
   int lnFiles   = ( nfiles    == -1 ? nFiles   : nfiles    ); 
   int lInstance = ( instance  == -1 ? Instance : instance  ); 

   int Seed; 
   int Blossom;
   
   GetSeedAndBlossom( Seed, Blossom, lInstance , lnPoints, lnFiles  );

   return Point - Seed;

}

void FileRecord::GetSeedAndBlossom( int & Seed , int & Blossom , int instance, int npoints, int nfiles ) const
{
  if(kUseNew)
    {
    GetSeedAndBlossomNew(Seed,Blossom,instance,npoints,nfiles);
    }
  else
    {
      GetSeedAndBlossomOld(Seed,Blossom,instance,npoints,nfiles);
    }
}

void FileRecord::GetSeedAndBlossomOld( int & Seed , int & Blossom , int instance, int npoints, int nfiles ) const
{



  int lnPoints  = ( npoints   == -1 ? nPoints  : npoints   );
  int lnFiles   = ( nfiles    == -1 ? nFiles   : nfiles    );
  int lInstance = ( instance  == -1 ? Instance : instance  );

  Seed    = lInstance*int(lnPoints/lnFiles);

  if ( Instance == lnFiles -1 )
    Blossom = lnPoints;
  else
    Blossom = ( Seed + int(lnPoints/lnFiles) >= lnPoints ?
	      lnPoints :
	      Seed+int(lnPoints/lnFiles) );

  //std::cout << Seed << " " << lnPoints << " " << lnFiles << " " << int(lnPoints/lnFiles)  << " " << lnFiles << " " << Blossom << std::\
  endl;                                                                                                                                  

}

void FileRecord::GetSeedAndBlossomNew( int & Seed , int & Blossom , int instance, int npoints, int nfiles ) const
{


  int lnPoints  = ( npoints   == -1 ? nPoints  : npoints   ); 
  int lnFiles   = ( nfiles    == -1 ? nFiles   : nfiles    ); 
  int lInstance = ( instance  == -1 ? Instance : instance  ); 

  int pointsperfile=int(lnPoints/lnFiles);
  int remainder=0;
  int npointsThisfile=pointsperfile;
  if(pointsperfile*lnFiles<lnPoints)
    {
      remainder=lnPoints-pointsperfile*lnFiles;
    }
  if(lInstance<remainder)
    {
      npointsThisfile=pointsperfile+1;
      Seed    = lInstance*npointsThisfile;
    }
  else
    {
      Seed=lInstance*npointsThisfile+remainder;
    }

  Blossom = ( Seed + npointsThisfile >= lnPoints ? 
                                lnPoints : 
                                Seed+npointsThisfile );


}


////////////
//
//  Build the input files as TChains 
//  based on the variables stored 
//  in the card reader
//
///
//TChain * FileRecord::CheckAndBuildChain( CardReader * MasterCard, std::string mode )
TChain * FileRecord::CheckAndBuildChain( CardReader * MasterCard, std::string mode, std::string instance  , std::string rep )
{
  std::cout << "CheckAndBuildChain" << std::endl;
   std::stringstream ss;  

   //// 
   //   First Build the Mother Tree
   ///

   // Get the list of files from 
   // CardReader
   std::map< std::string, std::string > file_list;
   ss.str(""); ss << "file_" << mode ;
   MasterCard->BuildListOfStrings( ss.str().c_str() , file_list );

   std::string treename;
   ss.str(""); ss << "tree_" << mode ;
   bool found = MasterCard->GetKey( ss.str().c_str() , treename );

   if( ! found )
   {
      std::cout << "FileRecord::CheckAndBuildChain " << ss.str() << " was not found, returning null. " << std::endl; 
      return 0;
   }

   std::cout << "size of file_list is " << file_list.size() << std::endl;
   std::cout << "treename is " << treename << std::endl;

   TChain * lchain = new TChain( treename.c_str() );

   // search and replace
   TString SR;

   // now add all of the files to the chain 
   int count = 0;
   std::map< std::string, std::string >::iterator item;
   for( item = file_list.begin() ; item != file_list.end() ; item ++ )
   {
      SR = item->second.c_str() ;

      if( instance !=  "-1" )
      {
          ss.str(""); ss << instance ;
          SR.ReplaceAll( rep.c_str() , ss.str().c_str() ); 
      }

      count = lchain->Add( SR.Data() );
      std::cout << "   adding: " << SR.Data()  
                << "  , " << count << " files " 
                << std::endl;
   }

   //////////////////////
   //// need to account for the possibility of 
   ///  friend trees, this is relevant for the sk1 tau processing
   ///   build them next
   TChain *friends[10],*fsifriends[10];   
   unsigned nFriends = 0 ;

   // Get the list of friends files from 
   // CardReader
   std::map< std::string, std::string > friends_list,fsi_friends_list;
   ss.str(""); ss << "friend_" << mode ;
   MasterCard->BuildListOfStrings( ss.str().c_str() , friends_list );
   ss.str(""); ss << "fsi_" << mode ;
   MasterCard->BuildListOfStrings( ss.str().c_str() , fsi_friends_list );

   // Get the name of the tree for our friends
   std::string friendTreeName;
   std::map< std::string, std::string > friendtree_list,fsi_friendtree_list;
   ss.str(""); ss << "friendtree_" << mode ;
   MasterCard->BuildListOfStrings( ss.str().c_str() , friendtree_list );
   ss.str(""); ss << "fsitree_" << mode ;
   MasterCard->BuildListOfStrings( ss.str().c_str() , fsi_friendtree_list );

   TString finder;
   std::string friendTreeKey;

   // now we loop over all of the existing friends
   std::cout << "size of friends_list is " << friends_list.size() << std::endl;
   for( item = friends_list.begin() ; item != friends_list.end() ; item ++ )
   {
      
      // get the Key
      finder = item->first.c_str();
      finder.ReplaceAll( "friend_" , "friendtree_" ); 
      friendTreeKey  = finder.Data();
      friendTreeName = friendtree_list[ friendTreeKey ];
      std::cout << "friendTreeName :" << friendTreeName << ": friendTreeKey :"  <<friendTreeKey <<":" <<  std::endl;

      friends[ nFriends ] = new TChain( friendTreeName.c_str() );

      SR = item->second.c_str() ;
      if( instance !=  "-1" )
      {
          ss.str(""); ss << instance ;
          SR.ReplaceAll( rep.c_str() , ss.str().c_str() ); 
      }
    

      // now add all of the files to the chain 
      count = friends[ nFriends ]->Add( SR.Data() );
      std::cout << "   adding: "   << SR.Data() 
                << "   as friend, " << count << " files , nFriends :" << nFriends <<std::endl; 

      lchain->AddFriend( friends[ nFriends ] );
      nFriends++;
   }

   std::cout << "size of fsi_friends_list is " << fsi_friends_list.size() << std::endl;
   // now we loop over all of the existing FSI friends
   nFriends = 0 ;
   for( item = fsi_friends_list.begin() ; item != fsi_friends_list.end() ; item ++ )
   {
      
      // get the Key
      finder = item->first.c_str();
      finder.ReplaceAll( "fsi_" , "fsitree_" ); 
      friendTreeKey  = finder.Data();
      friendTreeName = fsi_friendtree_list[ friendTreeKey ];
      std::cout << "FSI friendTreeName :" << friendTreeName << ": friendTreeKey :"  <<friendTreeKey <<":" <<  std::endl;

      fsifriends[ nFriends ] = new TChain( friendTreeName.c_str() );

      SR = item->second.c_str() ;
      if( instance !=  "-1" )
      {
          ss.str(""); ss << instance ;
          SR.ReplaceAll( rep.c_str() , ss.str().c_str() ); 
      }
    

      // now add all of the files to the chain 
      count = fsifriends[ nFriends ]->Add( SR.Data() );
      std::cout << "   adding: "   << SR.Data() 
                << "   as friend, " << count << " files , nFriends :" << nFriends <<std::endl; 

      std::cout << "   adding: "   << SR.Data() 
                << "   as friend, " << count << " files , nFriends :" << nFriends <<std::endl; 

      lchain->AddFriend( fsifriends[ nFriends ] );
      nFriends++;
   }


  return lchain;
}



