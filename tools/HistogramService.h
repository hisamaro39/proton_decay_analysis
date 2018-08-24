/**
 * @brief Histogramming class that stores the histograms for all systematics in a map.
 * @author Danilo Enoque Ferreira de Lima <dferreir@cern.ch>
 */
#ifndef HISTOGRAMSERVICE_H
#define HISTOGRAMSERVICE_H

#include <string>
#include "TH1D.h"
#include "TH2D.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TProfile.h"
#include <vector>
#include <map>
#include "TFile.h"
#include "TTree.h"

class HistogramService {

  public:
    HistogramService(const std::string &outputname = "out.root");
    virtual ~HistogramService();

    void addSystematics(const std::string &syst);
    void addTrigger(const std::string &trigger);
    void WriteOutput();

    void create1DVar(const std::string &name, const std::string &title, int binsX, double *bins, bool allSyst = true);
    void create1D(const std::string &name, const std::string &title, int binsX, double lowX, double highX, bool allSyst = true);
    void create2D(const std::string &name, const std::string &title, int binsX, double lowX, double highX, int binsY, double lowY, double highY, bool allSyst = true);
    void create2DVar(const std::string &name, const std::string &title, int binsX, double *vecX, int binsY, double *vecY, bool allSyst = true);

//    void create1DVarP(const std::string &name, const std::string &title, int binsX, double *bins, double lowY, double highY, bool allSyst = true);
    void create1DP(const std::string &name, const std::string &title, int binsX, double lowX, double highX, double lowY, double highY, bool allSyst = true);

    //TH1D *&h1D(const std::string &name, const std::string &trigger, const std::string &systematics);
    //TH2D *&h2D(const std::string &name, const std::string &trigger, const std::string &systematics);
    TH1F *&h1D(const std::string &name, const std::string &trigger, const std::string &systematics);
    TH2F *&h2D(const std::string &name, const std::string &trigger, const std::string &systematics);
    TProfile *&p1D(const std::string &name, const std::string &trigger, const std::string &systematics);

    std::string concat(const std::string &name, const std::string &suffix);

    void createVS(const std::string &name);
    std::vector<std::string> *&vs(const std::string &name);

    TFile *m_treeFile;
    TTree *m_tree;

  protected:

    TH1F *m_dummy1D;
    TH2F *m_dummy2D; // useless histograms which will be filled whenever histograms have no systematics to be filled
    TProfile *m_dummy1DP;
    std::vector<std::string> *m_dummyVS;

    std::string m_filename;
    TFile *m_file;


    std::vector<std::string> m_listOfSystematics;
    std::vector<std::string> m_listOfTriggers;

    bool m_lockTriggerAndSyst;
    // the first key is the systematics label
    // second key is the trigger name formed with [electron trigger]_[muon trigger]
    // the third key is the name
    //std::map<std::string, std::map<std::string, std::map<std::string, TH1D *> > > h_hist1D;
    //std::map<std::string, std::map<std::string, std::map<std::string, TH2D *> > > h_hist2D;
    std::map<std::string, std::map<std::string, std::map<std::string, TH1F *> > > h_hist1D;
    std::map<std::string, std::map<std::string, std::map<std::string, TH2F *> > > h_hist2D;
    std::map<std::string, std::map<std::string, std::map<std::string, TProfile *> > > h_profile1D;
    std::map<std::string, std::vector<std::string> *> m_vs;

    // lists of histograms with no syst. variations
    std::vector<std::string> m_noSystHist1D;
    std::vector<std::string> m_noSystHist2D;
    std::vector<std::string> m_noSystProfile1D;

  public:


    inline const std::vector<std::string> &listOfSystematics() const {return m_listOfSystematics;}
    inline const std::vector<std::string> &listOfTriggers() const {return m_listOfTriggers;}
};

#endif

