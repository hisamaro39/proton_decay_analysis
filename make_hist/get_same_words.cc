#include <vector>
#include <iostream>
using namespace std;
string get_same_words(vector<string> name){
  vector<string> list;
  vector< vector<string> > list_set;

  //name.push_back("mom_all_ring_reco_nring3_cut3_nring1_mulike1_michel0");
  //name.push_back("mom_all_ring_reco_nring3_cut3_nring1_mulike2_michel0");
  //name.push_back("mom_all_ring_reco_nring3_cut3_nring1_mulike3_michel0");

  for(unsigned int n=0;n<name.size();n++){
    int separation_point=0;
    for (unsigned int w=0;w<name[n].size();w++){
      if(name[n].at(w)=='_') {
        list.push_back(name[n].substr(separation_point,w-separation_point));
        separation_point=w+1;
      }
    }
    list_set.push_back(list);
    list.clear();
  }
  
  string final_name;
  vector<unsigned int> same_list_id;
  for(unsigned int l=0;l<list_set[0].size();l++){
    string list1 = list_set[0].at(l);
    final_name += list1;
    final_name += "_";
    for (unsigned int ls=1;ls<list_set.size();ls++){
      for(unsigned int l2=0;l2<list_set[ls].size();l2++){
        string list2 = list_set[ls].at(l2);
        if(list1==list2){
          same_list_id.push_back(l2);
          break;
        }
        else{
          bool add_name=true;
          for(unsigned int s=0;s<same_list_id.size();s++){
            if(l2==same_list_id[s]) add_name=false;
          }
          if(!add_name) continue;
          final_name += list2;
          final_name += "_";
        }
      }
    }
  }
  final_name = "output_hist/" + final_name.substr(0,final_name.size()-1) + ".pdf";
  return final_name;

}
