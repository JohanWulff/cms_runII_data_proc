#ifndef KINFITTER_H_
#define KINFITTER_H_

// C++
#include <iostream>
#include <string>
#include <map>
#include <vector>
#include <set>
#include <stdexcept>
#include <utility>

// ROOT
#include <Math/VectorUtil.h>
#include <Math/LorentzVector.h>
#include <Math/PtEtaPhiM4D.h>
#include <Math/PxPyPzM4D.h>
#include <TFile.h>
#include <TTree.h>
#include <TTreeReader.h>
#include <TTreeReaderValue.h>
#include <TLorentzVector.h>
#include <TMatrixD.h>

const int Z_MASS  = 91;  //GeV
const int H_MASS  = 125; //GeV

class KinFitter {
    // class to calculate KinFit mass of a particle given a mass hypothesis and its decay products
private:
    std::vector<float> kinINinfo;

    TLorentzVector tlv_l1 = TLorentzVector();
    TLorentzVector tlv_l2 = TLorentzVector();
    TLorentzVector tlv_b1 = TLorentzVector();
    TLorentzVector tlv_b2 = TLorentzVector();
    TVector2 ptmiss = TVector2();
    TMatrixD metcov = TMatrixD(2,2);

    std::pair<float,float> _fit(int mh1_hp, int mh2_hp);

public:
    KinFitter(std::vector<float> kinINinfo);
    ~KinFitter();
    std::pair<float,float> fit(std::string sgnHp);

};

#endif /* KINFITTER_H_ */