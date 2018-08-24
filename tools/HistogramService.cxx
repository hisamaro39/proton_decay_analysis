/**
 * @brief Histogramming class that stores the histograms for all systematics in a map.
 * @author Danilo Enoque Ferreira de Lima <dferreir@cern.ch>
 */
#include "tools/HistogramService.h"

#include <iostream>
#include "TFile.h"
#include <string>
#include <vector>
#include <map>
#include <algorithm>

#include "TH1F.h"
#include "TH2F.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TProfile.h"

HistogramService::HistogramService(const std::string &filename)
  : m_filename(filename) {

  //m_treeFile = 0;
  m_tree = 0;
  std::string fn = filename;
  size_t p = fn.find(".root");
  if (p != std::string::npos) {
    fn = fn.substr(0, p);
  }
  fn+= "_mini.root";
  //m_treeFile = new TFile(fn.c_str(), "recreate");
  //m_treeFile->cd();
  m_tree = new TTree("mini", "");

  std::cout << "m_filename=" << m_filename << std::endl;
  m_file = new TFile(m_filename.c_str(), "RECREATE");
  m_file->cd();
  m_lockTriggerAndSyst = false;
  m_dummy1D = new TH1F("dummy1D", "dummy 1d", 1, 0.0, 1.0);
  m_dummy1D->SetDirectory(0);
  m_dummy1DP = new TProfile("dummy1D Profile", "dummy 1d P", 1, 0.0, 1.0, 0.0, 1.0);
  m_dummy1DP->SetDirectory(0);  
  m_dummy2D = new TH2F("dummy2D", "dummy 2d", 1, 0.0, 1.0, 1, 0.0, 1.0);
  m_dummy2D->SetDirectory(0);
  m_dummyVS = new std::vector<std::string>;
}

std::string HistogramService::concat(const std::string &a, const std::string &s) {
  std::string out = a;
  out += s;
  return out;
}

HistogramService::~HistogramService(){};

void HistogramService::WriteOutput() {
  std::cout << "~HistogramService" << std::endl;
  delete m_dummyVS;
  delete m_dummy1D;
  delete m_dummy1DP;
  delete m_dummy2D;

  m_file->cd();
  for (std::vector<std::string>::const_iterator it = m_listOfTriggers.begin(); it != m_listOfTriggers.end(); ++it) {
    //m_file->cd(std::string(m_filename+std::string(":/")+*it).c_str());
    for (std::vector<std::string>::const_iterator jt = m_listOfSystematics.begin(); jt != m_listOfSystematics.end(); ++jt) {
      //std::cout << "Fill some histograms !!" << std::endl;
      for (std::map<std::string, TH1F *>::const_iterator kt = h_hist1D[*it][*jt].begin(); kt != h_hist1D[*it][*jt].end(); ++kt) {
        //std::cout << "Writing " << *it << " " << *jt << " " << kt->first << std::endl;
        //std::cout << "entry is " << h_hist1D[*it][*jt][kt->first]->GetEntries() << std::endl;;
        h_hist1D[*it][*jt][kt->first]->Write();
        //m_file->Add(h_hist1D[*it][*jt][kt->first]);
      }
      for (std::map<std::string, TProfile *>::const_iterator kt = h_profile1D[*it][*jt].begin(); kt != h_profile1D[*it][*jt].end(); ++kt) {
        //std::cout << "Writing " << *it << " " << *jt << " " << kt->first << std::endl;
        h_profile1D[*it][*jt][kt->first]->Write();
      }
      for (std::map<std::string, TH2F *>::const_iterator kt = h_hist2D[*it][*jt].begin(); kt != h_hist2D[*it][*jt].end(); ++kt) {
        h_hist2D[*it][*jt][kt->first]->Write();
      }
    }
  }
  //m_file->ls();
  m_file->cd("/");
  m_file->Close();
  for (std::map<std::string, std::vector<std::string> *>::const_iterator kt = m_vs.begin(); kt != m_vs.end(); ++kt) {
    gDirectory->WriteObject(m_vs[kt->first], kt->first.c_str());
  }
  /*if (m_tree && m_treeFile) {
    m_treeFile->cd();
    m_tree->Write();
    m_treeFile->Close();
    //delete m_treeFile;
  }*/
}

void HistogramService::addSystematics(const std::string &syst) {
  if (m_lockTriggerAndSyst) {
    std::cout << "HistogramService: Cannot add a new trigger and systematics after histograms are added." <<std::endl;
    std::cout << "HistogramService: This happened when adding systematics \""<<syst<<"\"." << std::endl;
    return;
  }
  if (std::find(m_listOfSystematics.begin(), m_listOfSystematics.end(), syst) == m_listOfSystematics.end())
    m_listOfSystematics.push_back(syst);
}

void HistogramService::addTrigger(const std::string &trigger) {
  if (m_lockTriggerAndSyst) {
    std::cout << "HistogramService: Cannot add a new trigger and systematics after histograms are added." <<std::endl;
    std::cout << "HistogramService: This happened when adding trigger \""<<trigger<<"\"." << std::endl;
    return;
  }
  m_file->cd();
  m_file->cd(std::string(m_filename+std::string(":/")).c_str());
  m_file->mkdir(trigger.c_str());
  m_file->cd(std::string(m_filename+std::string(":/")+trigger).c_str());
  m_listOfTriggers.push_back(trigger);
}

void HistogramService::create1D(const std::string &name, const std::string &title, int binsX, double lowX, double highX, bool allSyst) {
  m_lockTriggerAndSyst = true;
  //std::cout << "create1D!!" << std::endl;
  
  m_file->cd();
  for (std::vector<std::string>::const_iterator it = m_listOfTriggers.begin(); it != m_listOfTriggers.end(); ++it) {
    //std::cout << "HistogramService::create1D: entering directory " << std::string(m_filename+std::string(":/")+*it).c_str() << std::endl;
    m_file->cd(std::string(m_filename+std::string(":/")+*it).c_str());
    //std::cout << "allSyst=" << allSyst << std::endl;
    /*if (allSyst) {
      for (std::vector<std::string>::const_iterator jt = m_listOfSystematics.begin(); jt != m_listOfSystematics.end(); ++jt) {
        if (h_hist1D[*it][*jt].find(name) != h_hist1D[*it][*jt].end()) {
          if (name.find("Cutflow") == std::string::npos) { // don't complain about the Cutflow histograms, which will always be added more than once
            std::cout << "HistogramService: Histogram named \"" << name << "\" already exists. Giving up creating it." << std::endl;
          }
          return;
        }
        std::string hname = concat(name, *jt);
        std::cout << "Dir " << *it << ", syst " << *jt << ", name = " << name << ", hname = " << hname << std::endl;
        h_hist1D[*it][*jt][name] = new TH1F(hname.c_str(), title.c_str(), binsX, lowX, highX);
        h_hist1D[*it][*jt][name]->Sumw2();
      }
    } else { // only create nominal histogram
      if (h_hist1D[*it][""].find(name) != h_hist1D[*it][""].end()) {
        if (name.find("Cutflow") == std::string::npos) { // don't complain about the Cutflow histograms, which will always be added more than once
          std::cout << "HistogramService: Histogram named \"" << name << "\" already exists. Giving up creating it." << std::endl;
        }
        return;
      }*/
      h_hist1D[*it][""][name] = new TH1F(name.c_str(), title.c_str(), binsX, lowX, highX);
      h_hist1D[*it][""][name]->Sumw2();
      m_noSystHist1D.push_back(name);
      //std::cout << "size of m_noSystHist1D is " << m_noSystHist1D.size() << std::endl;
    //}
  }
}

void HistogramService::create1DVar(const std::string &name, const std::string &title, int binsX, double *b, bool allSyst) {
  m_lockTriggerAndSyst = true;
  
  m_file->cd();
  for (std::vector<std::string>::const_iterator it = m_listOfTriggers.begin(); it != m_listOfTriggers.end(); ++it) {
    //std::cout << "HistogramService::create1D: entering directory " << std::string(m_filename+std::string(":/")+*it).c_str() << std::endl;
    m_file->cd(std::string(m_filename+std::string(":/")+*it).c_str());
    if (allSyst) {
      for (std::vector<std::string>::const_iterator jt = m_listOfSystematics.begin(); jt != m_listOfSystematics.end(); ++jt) {
        if (h_hist1D[*it][*jt].find(name) != h_hist1D[*it][*jt].end()) {
          if (name.find("Cutflow") == std::string::npos) { // don't complain about the Cutflow histograms, which will always be added more than once
            std::cout << "HistogramService: Histogram named \"" << name << "\" already exists. Giving up creating it." << std::endl;
          }
          return;
        }
        std::string hname = concat(name, *jt);
        //std::cout << "Dir " << *it << ", syst " << *jt << ", name = " << name << ", hname = " << hname << std::endl;
        h_hist1D[*it][*jt][name] = new TH1F(hname.c_str(), title.c_str(), binsX, b);
        h_hist1D[*it][*jt][name]->Sumw2();
      }
    } else { // only create nominal histogram
      if (h_hist1D[*it][""].find(name) != h_hist1D[*it][""].end()) {
        if (name.find("Cutflow") == std::string::npos) { // don't complain about the Cutflow histograms, which will always be added more than once
          std::cout << "HistogramService: Histogram named \"" << name << "\" already exists. Giving up creating it." << std::endl;
        }
        return;
      }
      h_hist1D[*it][""][name] = new TH1F(name.c_str(), title.c_str(), binsX, b);
      h_hist1D[*it][""][name]->Sumw2();
      m_noSystHist1D.push_back(name);
    }
  }
}

void HistogramService::create1DP(const std::string &name, const std::string &title, int binsX, double lowX, double highX, double lowY, double highY, bool allSyst) {
  m_lockTriggerAndSyst = true;
  
  m_file->cd();
  for (std::vector<std::string>::const_iterator it = m_listOfTriggers.begin(); it != m_listOfTriggers.end(); ++it) {
    //std::cout << "HistogramService::create1D: entering directory " << std::string(m_filename+std::string(":/")+*it).c_str() << std::endl;
    m_file->cd(std::string(m_filename+std::string(":/")+*it).c_str());
    if (allSyst) {
      for (std::vector<std::string>::const_iterator jt = m_listOfSystematics.begin(); jt != m_listOfSystematics.end(); ++jt) {
        if (h_profile1D[*it][*jt].find(name) != h_profile1D[*it][*jt].end()) {
          if (name.find("Cutflow") == std::string::npos) { // don't complain about the Cutflow histograms, which will always be added more than once
            std::cout << "HistogramService: Histogram named \"" << name << "\" already exists. Giving up creating it." << std::endl;
          }
          return;
        }
        std::string hname = concat(name, *jt);
        //std::cout << "Dir " << *it << ", syst " << *jt << ", name = " << name << ", hname = " << hname << std::endl;
        h_profile1D[*it][*jt][name] = new TProfile(hname.c_str(), title.c_str(), binsX, lowX, highX, lowY, highY);
        h_profile1D[*it][*jt][name]->Sumw2();
      }
    } else { // only create nominal histogram
      if (h_profile1D[*it][""].find(name) != h_profile1D[*it][""].end()) {
        if (name.find("Cutflow") == std::string::npos) { // don't complain about the Cutflow histograms, which will always be added more than once
          std::cout << "HistogramService: Histogram named \"" << name << "\" already exists. Giving up creating it." << std::endl;
        }
        return;
      }
      h_profile1D[*it][""][name] = new TProfile(name.c_str(), title.c_str(), binsX, lowX, highX,lowY, highY);
      h_profile1D[*it][""][name]->Sumw2();
      m_noSystProfile1D.push_back(name);
    }
  }
}

void HistogramService::createVS(const std::string &name) {
  m_file->cd();
  if (m_vs.find(name) != m_vs.end()) {
    std::cout << "HistogramService: String vector named \"" << name << "\" already exists. Giving up creating it." << std::endl;
    return;
  }
  m_vs[name] = new std::vector<std::string>;
}

void HistogramService::create2D(const std::string &name, const std::string &title, int binsX, double lowX, double highX, int binsY, double lowY, double highY, bool allSyst) {
  m_lockTriggerAndSyst = true;
  
  m_file->cd();
  for (std::vector<std::string>::const_iterator it = m_listOfTriggers.begin(); it != m_listOfTriggers.end(); ++it) {
    m_file->cd(std::string(m_filename+std::string(":/")+*it).c_str());
    /*if (allSyst) {
      for (std::vector<std::string>::const_iterator jt = m_listOfSystematics.begin(); jt != m_listOfSystematics.end(); ++jt) {
        if (h_hist2D[*it][*jt].find(name) != h_hist2D[*it][*jt].end()) {
          std::cout << "HistogramService: 2D histogram named \"" << name << "\" already exists. Giving up creating it." << std::endl;
          return;
        }
        h_hist2D[*it][*jt][name] = new TH2F(concat(name, *jt).c_str(), title.c_str(), binsX, lowX, highX, binsY, lowY, highY);
        h_hist2D[*it][*jt][name]->Sumw2();
      }
    } else { // only for nominal
      if (h_hist2D[*it][""].find(name) != h_hist2D[*it][""].end()) {
        std::cout << "HistogramService: 2D histogram named \"" << name << "\" already exists. Giving up creating it." << std::endl;
        return;
      }*/
      h_hist2D[*it][""][name] = new TH2F(name.c_str(), title.c_str(), binsX, lowX, highX, binsY, lowY, highY);
      h_hist2D[*it][""][name]->Sumw2();
      m_noSystHist2D.push_back(name);
    //}
  }
}


void HistogramService::create2DVar(const std::string &name, const std::string &title, int binsX, double *vecX, int binsY, double *vecY, bool allSyst) {
  m_lockTriggerAndSyst = true;
  
  m_file->cd();
  for (std::vector<std::string>::const_iterator it = m_listOfTriggers.begin(); it != m_listOfTriggers.end(); ++it) {
    m_file->cd(std::string(m_filename+std::string(":/")+*it).c_str());
    if (allSyst) {
      for (std::vector<std::string>::const_iterator jt = m_listOfSystematics.begin(); jt != m_listOfSystematics.end(); ++jt) {
        if (h_hist2D[*it][*jt].find(name) != h_hist2D[*it][*jt].end()) {
          std::cout << "HistogramService: 2D histogram named \"" << name << "\" already exists. Giving up creating it." << std::endl;
          return;
        }
        h_hist2D[*it][*jt][name] = new TH2F(concat(name, *jt).c_str(), title.c_str(), binsX, vecX, binsY, vecY);
        h_hist2D[*it][*jt][name]->Sumw2();
      }
    } else { // only for nominal
      if (h_hist2D[*it][""].find(name) != h_hist2D[*it][""].end()) {
        std::cout << "HistogramService: 2D histogram named \"" << name << "\" already exists. Giving up creating it." << std::endl;
        return;
      }
      h_hist2D[*it][""][name] = new TH2F(name.c_str(), title.c_str(), binsX, vecX, binsY, vecY);
      h_hist2D[*it][""][name]->Sumw2();
      m_noSystHist2D.push_back(name);
    }
  }
}


TH1F *&HistogramService::h1D(const std::string &name, const std::string &trigger, const std::string &systematics) {
  //std::cout << "Getting histogram " << name << " " << trigger << " " << systematics<< std::endl;
  if ( (systematics != "") && (std::find(m_noSystHist1D.begin(), m_noSystHist1D.end(), name) != m_noSystHist1D.end()) )
    return m_dummy1D;
  //std::cout << "No systematics for " << name << std::endl;
  // this returns a dummy histogram if the histogram was not setup
  if (h_hist1D[trigger][systematics].find(name) == h_hist1D[trigger][systematics].end())
    return m_dummy1D;
  //std::cout << "no histogram found " << name << ", " << trigger << ", " << systematics << std::endl;

  return h_hist1D[trigger][systematics][name];
}

TProfile *&HistogramService::p1D(const std::string &name, const std::string &trigger, const std::string &systematics) {
  //std::cout << "Getting histogram " << name << " " << trigger << " " << systematics<< std::endl;
  if ( (systematics != "") && (std::find(m_noSystProfile1D.begin(), m_noSystProfile1D.end(), name) != m_noSystProfile1D.end()) )
    return m_dummy1DP;
  //std::cout << "No systematics for " << name << std::endl;
  // this returns a dummy histogram if the histogram was not setup
  if (h_profile1D[trigger][systematics].find(name) == h_profile1D[trigger][systematics].end())
    return m_dummy1DP;
  //std::cout << "no histogram found " << name << ", " << trigger << ", " << systematics << std::endl;

  return h_profile1D[trigger][systematics][name];
}

std::vector<std::string> *&HistogramService::vs(const std::string &name) {
  return m_vs[name];
}


TH2F *&HistogramService::h2D(const std::string &name, const std::string &trigger, const std::string &systematics) {
  if ( (systematics != "") && (std::find(m_noSystHist2D.begin(), m_noSystHist2D.end(), name) != m_noSystHist2D.end()) )
    return m_dummy2D;
  // this returns a dummy histogram if the histogram was not setup
  if (h_hist2D[trigger][systematics].find(name) == h_hist2D[trigger][systematics].end())
    return m_dummy2D;
  return h_hist2D[trigger][systematics][name];
}

