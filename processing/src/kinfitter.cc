#include "../../../HHKinFit2/HHKinFit2/interface/HHKinFitMasterHeavyHiggs.h"
#include "../../../HHKinFit2/HHKinFit2/interface/exceptions/HHInvMConstraintException.h"
#include "../../../HHKinFit2/HHKinFit2/interface/exceptions/HHEnergyRangeException.h"
#include "../../../HHKinFit2/HHKinFit2/interface/exceptions/HHEnergyConstraintException.h"
#include "../../../HHKinFit2/HHKinFit2/interface/exceptions/HHLimitSettingException.h"
#include "cms_runII_data_proc/processing/interface/kinfitter.hh"

KinFitter::KinFitter(std::vector<float> kinINinfo) {
    tlv_l1.SetPtEtaPhiM(kinINinfo[0], kinINinfo[1], kinINinfo[2],kinINinfo[3]);
    tlv_l2.SetPtEtaPhiM(kinINinfo[4], kinINinfo[5], kinINinfo[6],kinINinfo[7]);
    tlv_b1.SetPtEtaPhiM(kinINinfo[8], kinINinfo[9], kinINinfo[10],kinINinfo[11]);
    tlv_b2.SetPtEtaPhiM(kinINinfo[12], kinINinfo[13], kinINinfo[14],kinINinfo[15]);
    ptmiss.SetMagPhi(kinINinfo[16], kinINinfo[17]);    
    metcov(0,0) = static_cast<double>(kinINinfo[18]);
    metcov(1,0) = static_cast<double>(kinINinfo[19]);
    metcov(0,1) = static_cast<double>(kinINinfo[19]);
    metcov(1,1) = static_cast<double>(kinINinfo[20]);
}

KinFitter::~KinFitter() {}

std::pair<float,float> KinFitter::_fit(int mh1_hp, int mh2_hp) {
    HHKinFit2::HHKinFitMasterHeavyHiggs KinFit = HHKinFit2::HHKinFitMasterHeavyHiggs(tlv_b1, tlv_b2, tlv_l1, tlv_l2, ptmiss, metcov);
    KinFit.addHypo(mh1_hp, mh2_hp);

    float HHKmass = -999;
    float HHKChi2 = -999;
    bool DEBUG = false;

    bool wrongHHK=false;
    try { KinFit.fit(); }
    catch (HHKinFit2::HHInvMConstraintException const& e) {
        std::cout<<"THIS FIT DID NOT CONVERGE: INV MASS CONSTRAIN EXCEPTION (mh1_hp = " << mh1_hp << " , mh2_hp = " << mh2_hp << ")" << std::endl;
        if (DEBUG) {            
            std::cout<<"INVME Tau1"<<std::endl;
            std::cout<<"INVME (E,Px,Py,Pz,M) "<< tlv_l1.E() <<","<< tlv_l1.Px() <<","<< tlv_l1.Py() <<","<< tlv_l1.Pz() <<","<< tlv_l1.M() << std::endl;
            std::cout<<"INVME Tau2"<<std::endl;
            std::cout<<"INVME (E,Px,Py,Pz,M) "<< tlv_l2.E() <<","<< tlv_l2.Px() <<","<< tlv_l2.Py() <<","<< tlv_l2.Pz() <<","<< tlv_l2.M() << std::endl;
            std::cout<<"INVME B1"<<std::endl;
            std::cout<<"INVME (E,Px,Py,Pz,M) "<< tlv_b1.E() <<","<< tlv_b1.Px() <<","<< tlv_b1.Py() <<","<< tlv_b1.Pz() <<","<< tlv_b1.M() << std::endl;
            std::cout<<"INVME B2"<<std::endl;
            std::cout<<"INVME (E,Px,Py,Pz,M) "<< tlv_b2.E() <<","<< tlv_b2.Px() <<","<< tlv_b2.Py() <<","<< tlv_b2.Pz() <<","<< tlv_b2.M() << std::endl;
            std::cout<<"INVME MET"<<std::endl;
            std::cout<<"INVME (E,Px,Py,Pz,M) "<<","<< ptmiss.Px() <<","<< ptmiss.Py() << std::endl;
            std::cout<<"INVME METCOV "<<std::endl;
            std::cout<<"INVME "<< metcov (0,0) <<"  "<< metcov (0,1) << std::endl;
            std::cout<<"INVME "<< metcov (1,0) <<"  "<< metcov (1,1) << std::endl;
            std::cout<<"INVME tau1, tau2, b1, b2"<<std::endl;
            std::cout<<"INVME ";
            tlv_l1.Print();
            std::cout<<"INVME ";
            tlv_l2.Print();
            std::cout<<"INVME ";
            tlv_b1.Print();
            std::cout<<"INVME ";
            tlv_b2.Print();
        }
        wrongHHK=true;
    }
    catch (HHKinFit2::HHEnergyRangeException const& e) {
        std::cout<<"THIS FIT DID NOT CONVERGE: INV MASS CONSTRAIN EXCEPTION (mh1_hp = " << mh1_hp << " , mh2_hp = " << mh2_hp << ")" << std::endl;
        if (DEBUG) {
            std::cout<<"ERANGE Tau1"<<std::endl;
            std::cout<<"ERANGE (E,Px,Py,Pz,M) "<< tlv_l1.E() <<","<< tlv_l1.Px() <<","<< tlv_l1.Py() <<","<< tlv_l1.Pz() <<","<< tlv_l1.M() << std::endl;
            std::cout<<"ERANGE Tau2"<<std::endl;
            std::cout<<"ERANGE (E,Px,Py,Pz,M) "<< tlv_l2.E() <<","<< tlv_l2.Px() <<","<< tlv_l2.Py() <<","<< tlv_l2.Pz() <<","<< tlv_l2.M() << std::endl;
            std::cout<<"ERANGE B1"<<std::endl;
            std::cout<<"ERANGE (E,Px,Py,Pz,M) "<< tlv_b1.E() <<","<< tlv_b1.Px() <<","<< tlv_b1.Py() <<","<< tlv_b1.Pz() <<","<< tlv_b1.M() << std::endl;
            std::cout<<"ERANGE B2"<<std::endl;
            std::cout<<"ERANGE (E,Px,Py,Pz,M) "<< tlv_b2.E() <<","<< tlv_b2.Px() <<","<< tlv_b2.Py() <<","<< tlv_b2.Pz() <<","<< tlv_b2.M() << std::endl;
            std::cout<<"ERANGE MET"<<std::endl;
            std::cout<<"ERANGE (E,Px,Py,Pz,M) "<<","<< ptmiss.Px() <<","<< ptmiss.Py() <<std::endl;
            std::cout<<"ERANGE METCOV "<<std::endl;
            std::cout<<"ERANGE "<< metcov (0,0) <<"  "<< metcov (0,1) << std::endl;
            std::cout<<"ERANGE "<< metcov (1,0) <<"  "<< metcov (1,1) << std::endl;
            std::cout<<"ERANGE tau1, tau2, b1, b2"<<std::endl;
            std::cout<<"ERANGE ";
            tlv_l1.Print();
            std::cout<<"ERANGE ";
            tlv_l2.Print();
            std::cout<<"ERANGE ";
            tlv_b1.Print();
            std::cout<<"ERANGE ";
            tlv_b2.Print();
        }
        wrongHHK=true;
    }
    catch (HHKinFit2::HHEnergyConstraintException const& e) {
        std::cout<<"THIS FIT DID NOT CONVERGE: ENERGY CONSTRAIN EXCEPTION (mh1_hp = " << mh1_hp << " , mh2_hp = " << mh2_hp << ")" << std::endl;
        if (DEBUG) {
            std::cout<<"ECON Tau1"<<std::endl;
            std::cout<<"ECON (E,Px,Py,Pz,M) "<< tlv_l1.E() <<","<< tlv_l1.Px() <<","<< tlv_l1.Py() <<","<< tlv_l1.Pz() <<","<< tlv_l1.M() << std::endl;
            std::cout<<"ECON Tau2"<<std::endl;
            std::cout<<"ECON (E,Px,Py,Pz,M) "<< tlv_l2.E() <<","<< tlv_l2.Px() <<","<< tlv_l2.Py() <<","<< tlv_l2.Pz() <<","<< tlv_l2.M() << std::endl;
            std::cout<<"ECON B1"<<std::endl;
            std::cout<<"ECON (E,Px,Py,Pz,M) "<< tlv_b1.E() <<","<< tlv_b1.Px() <<","<< tlv_b1.Py() <<","<< tlv_b1.Pz() <<","<< tlv_b1.M() << std::endl;
            std::cout<<"ECON B2"<<std::endl;
            std::cout<<"ECON (E,Px,Py,Pz,M) "<< tlv_b2.E() <<","<< tlv_b2.Px() <<","<< tlv_b2.Py() <<","<< tlv_b2.Pz() <<","<< tlv_b2.M() << std::endl;
            std::cout<<"ECON MET"<<std::endl;
            std::cout<<"ECON (E,Px,Py,Pz,M) "<<","<< ptmiss.Px() <<","<< ptmiss.Py() <<std::endl;
            std::cout<<"ECON METCOV "<<std::endl;
            std::cout<<"ECON "<< metcov (0,0) <<"  "<< metcov (0,1) << std::endl;
            std::cout<<"ECON "<< metcov (1,0) <<"  "<< metcov (1,1) << std::endl;
            std::cout<<"ECON tau1, tau2, b1, b2"<<std::endl;
            std::cout<<"ECON ";
            tlv_l1.Print();
            std::cout<<"ECON ";
            tlv_l2.Print();
            std::cout<<"ECON ";
            tlv_b1.Print();
            std::cout<<"ECON ";
            tlv_b2.Print();
        }
        wrongHHK=true;
    }
    if(!wrongHHK) {
      HHKmass = KinFit.getMH(mh1_hp, mh2_hp);
      HHKChi2 = KinFit.getChi2(mh1_hp, mh2_hp);
    }
    else {
      HHKmass = -333;
      HHKChi2 = -333;
    }
    std::pair<float,float> local_result = std::pair(HHKmass, HHKChi2);

    return local_result;
}

std::pair<float,float> KinFitter::fit(std::string sgnHp) {
    int mh1_hp = H_MASS;
    int mh2_hp = H_MASS;

    if (sgnHp == "ZZ") {
        mh1_hp = Z_MASS;
        mh2_hp = Z_MASS;
    }

    else if (sgnHp == "ZH") {
        mh1_hp = Z_MASS;
        mh2_hp = H_MASS;
    }

    std::pair<float,float> result = std::pair(-999,-999);

    if (mh1_hp == mh2_hp) {
        result = _fit(mh1_hp, mh2_hp);
    }
    else {
        std::pair<float,float> right_fit = _fit(mh1_hp, mh2_hp);
        std::pair<float,float> left_fit  = _fit(mh2_hp, mh1_hp);

        if (right_fit.second != -333 && left_fit.second != -333) {
            if (right_fit.second < left_fit.second) { result = right_fit; }
            else { result = left_fit; }
        }
        else if (right_fit.second != -333 || left_fit.second != -333) {
            if (right_fit.second != -333) { result = right_fit; }
            else { result = left_fit; }
        }
        else {
            std::cout << "Neither the  mh1_hp,mh2_hp=" << mh1_hp << "," << mh2_hp << " fit nor the mh1_hp,mh2_hp=" << mh2_hp << "," << mh1_hp << " fit converged!!" << std::endl;
        }
    }

    return result;
}
