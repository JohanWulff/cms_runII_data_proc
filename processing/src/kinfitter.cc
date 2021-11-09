#include "../../../HHKinFit2/HHKinFit2/interface/HHKinFitMasterHeavyHiggs.h"
#include "cms_runII_data_proc/processing/interface/kinfitter.hh"
#include <typeinfo>

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

    /* FIXME IF WANTED BUT NOT NECESSARY
    bool wrongHHK=false;
    try { KinFit.fit() ; }
    catch (HHKinFit2::HHInvMConstraintException &e) {
        cout<<"INVME THIS EVENT WAS WRONG, INV MASS CONSTRAIN EXCEPTION"<<endl;
        cout<<"INVME masshypo1 = %i,    masshypo2 = %i" % mh1_hp % mh2_hp<<endl;
        cout<<"INVME Tau1"<<endl;
        cout<<"INVME (E,Px,Py,Pz,M) "<< tlv_l1.E() <<","<< tlv_l1.Px() <<","<< tlv_l1.Py() <<","<< tlv_l1.Pz() <<","<< tlv_l1.M() << endl;
        cout<<"INVME Tau2"<<endl;
        cout<<"INVME (E,Px,Py,Pz,M) "<< tlv_l2.E() <<","<< tlv_l2.Px() <<","<< tlv_l2.Py() <<","<< tlv_l2.Pz() <<","<< tlv_l2.M() << endl;
        cout<<"INVME B1"<<endl;
        cout<<"INVME (E,Px,Py,Pz,M) "<< tlv_b1.E() <<","<< tlv_b1.Px() <<","<< tlv_b1.Py() <<","<< tlv_b1.Pz() <<","<< tlv_b1.M() << endl;
        cout<<"INVME B2"<<endl;
        cout<<"INVME (E,Px,Py,Pz,M) "<< tlv_b2.E() <<","<< tlv_b2.Px() <<","<< tlv_b2.Py() <<","<< tlv_b2.Pz() <<","<< tlv_b2.M() << endl;
        cout<<"INVME MET"<<endl;
        cout<<"INVME (E,Px,Py,Pz,M) "<<","<< ptmiss.Px() <<","<< ptmiss.Py() <<endl;
        cout<<"INVME METCOV "<<endl;
        cout<<"INVME "<< metcov (0,0) <<"  "<< metcov (0,1) << endl;
        cout<<"INVME "<< metcov (1,0) <<"  "<< metcov (1,1) << endl;
        cout<<"INVME tau1, tau2, b1, b2"<<endl;
        cout<<"INVME ";
        tlv_l1.Print();
        cout<<"INVME ";
        tlv_l2.Print();
        cout<<"INVME ";
        tlv_b1.Print();
        cout<<"INVME ";
        tlv_b2.Print();
        wrongHHK=true;
    }
    catch (HHKinFit2::HHEnergyRangeException &e) {
        cout<<"ERANGE THIS EVENT WAS WRONG, INV MASS CONSTRAIN EXCEPTION"<<endl;
        cout<<"ERANGE masshypo1 = %i,    masshypo2 = %i" % mh1_hp % mh2_hp<<endl;
        cout<<"ERANGE Tau1"<<endl;
        cout<<"ERANGE (E,Px,Py,Pz,M) "<< tlv_l1.E() <<","<< tlv_l1.Px() <<","<< tlv_l1.Py() <<","<< tlv_l1.Pz() <<","<< tlv_l1.M() << endl;
        cout<<"ERANGE Tau2"<<endl;
        cout<<"ERANGE (E,Px,Py,Pz,M) "<< tlv_l2.E() <<","<< tlv_l2.Px() <<","<< tlv_l2.Py() <<","<< tlv_l2.Pz() <<","<< tlv_l2.M() << endl;
        cout<<"ERANGE B1"<<endl;
        cout<<"ERANGE (E,Px,Py,Pz,M) "<< tlv_b1.E() <<","<< tlv_b1.Px() <<","<< tlv_b1.Py() <<","<< tlv_b1.Pz() <<","<< tlv_b1.M() << endl;
        cout<<"ERANGE B2"<<endl;
        cout<<"ERANGE (E,Px,Py,Pz,M) "<< tlv_b2.E() <<","<< tlv_b2.Px() <<","<< tlv_b2.Py() <<","<< tlv_b2.Pz() <<","<< tlv_b2.M() << endl;
        cout<<"ERANGE MET"<<endl;
        cout<<"ERANGE (E,Px,Py,Pz,M) "<<","<< ptmiss.Px() <<","<< ptmiss.Py() <<endl;
        cout<<"ERANGE METCOV "<<endl;
        cout<<"ERANGE "<< metcov (0,0) <<"  "<< metcov (0,1) << endl;
        cout<<"ERANGE "<< metcov (1,0) <<"  "<< metcov (1,1) << endl;
        cout<<"ERANGE tau1, tau2, b1, b2"<<endl;
        cout<<"ERANGE ";
        tlv_l1.Print();
        cout<<"ERANGE ";
        tlv_l2.Print();
        cout<<"ERANGE ";
        tlv_b1.Print();
        cout<<"ERANGE ";
        tlv_b2.Print();
        wrongHHK=true;
    }
    catch (HHKinFit2::HHEnergyConstraintException &e) {
        cout<<"ECON THIS EVENT WAS WRONG, ENERGY CONSTRAIN EXCEPTION"<<endl;
        cout<<"ECON masshypo1 = %i,    masshypo2 = %i" % mh1_hp % mh2_hp<<endl;
        cout<<"ECON Tau1"<<endl;
        cout<<"ECON (E,Px,Py,Pz,M) "<< tlv_l1.E() <<","<< tlv_l1.Px() <<","<< tlv_l1.Py() <<","<< tlv_l1.Pz() <<","<< tlv_l1.M() << endl;
        cout<<"ECON Tau2"<<endl;
        cout<<"ECON (E,Px,Py,Pz,M) "<< tlv_l2.E() <<","<< tlv_l2.Px() <<","<< tlv_l2.Py() <<","<< tlv_l2.Pz() <<","<< tlv_l2.M() << endl;
        cout<<"ECON B1"<<endl;
        cout<<"ECON (E,Px,Py,Pz,M) "<< tlv_b1.E() <<","<< tlv_b1.Px() <<","<< tlv_b1.Py() <<","<< tlv_b1.Pz() <<","<< tlv_b1.M() << endl;
        cout<<"ECON B2"<<endl;
        cout<<"ECON (E,Px,Py,Pz,M) "<< tlv_b2.E() <<","<< tlv_b2.Px() <<","<< tlv_b2.Py() <<","<< tlv_b2.Pz() <<","<< tlv_b2.M() << endl;
        cout<<"ECON MET"<<endl;
        cout<<"ECON (E,Px,Py,Pz,M) "<<","<< ptmiss.Px() <<","<< ptmiss.Py() <<endl;
        cout<<"ECON METCOV "<<endl;
        cout<<"ECON "<< metcov (0,0) <<"  "<< metcov (0,1) << endl;
        cout<<"ECON "<< metcov (1,0) <<"  "<< metcov (1,1) << endl;
        cout<<"ECON tau1, tau2, b1, b2"<<endl;
        cout<<"ECON ";
        tlv_l1.Print();
        cout<<"ECON ";
        tlv_l2.Print();
        cout<<"ECON ";
        tlv_b1.Print();
        cout<<"ECON ";
        tlv_b2.Print();
        wrongHHK=true;
    }
    if(!wrongHHK) {
      HHKmass = kinFits.getMH();
      HHKChi2 = kinFits.getChi2();
    }
    else {
      HHKmass = -333;
    }
    std::pair<float,float> result = std::pair(HHKmass, HHKChi2);
    */

    KinFit.fit();
    std::pair<float,float> result = std::pair(KinFit.getMH(), KinFit.getChi2());

    return result;
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

    return _fit(mh1_hp, mh2_hp);
}
