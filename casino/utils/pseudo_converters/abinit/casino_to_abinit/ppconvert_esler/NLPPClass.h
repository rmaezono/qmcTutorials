#ifndef NLPP_CLASS_H
#define NLPP_CLASS_H

#include <map>
#include <string>
#include <vector>
//#include <Common/Splines/CubicSpline.h>
#include "CubicSpline.h"
#include "XMLWriterClass2.h"

using namespace std;


class ChannelPotentialClass {

 public:
  int l, n_principal;
  // Vl is stored in hartrees
  CubicSpline Vl;
  CubicSpline ul;
  bool HasProjector;
  double Cutoff, Occupation;
  double FindCutoff();
  void WriteChannel (XMLWriterClass &writer, bool writeVl);
  ChannelPotentialClass() {
   HasProjector= false;
   Occupation = 0.0;
  }
};


enum XCType { XC_LDA, XC_GGA, XC_HF, XC_DF, XC_NONE};


class NLPPClass {

 private:
  map<string,double> UnitToHartreeMap;
  map<string,double> UnitToBohrMap;
  map<XCType,string> XCMap;
  map<string,XCType> XCRevMap;
  map<int,string> ChannelMap;
  map<string,int> ChannelRevMap;
  map<int, string> ZToSymbolMap;
  map<string, int> SymbolToZMap;
  void SetupMaps();
  vector<ChannelPotentialClass> ChannelPotentials;
  int AtomicNumber;
  double PseudoCharge;
  string EnergyUnit, LengthUnit;
  int LocalChannel;
  // The grid is stored in bohr
  GeneralGrid PotentialGrid;
  bool Relativistic;
  XCType XC;

 public:
  int GetNumChannels();
  void ReadCASINO_PP (string fileName);
  void ReadCASINO_WF (string fileName, int l);
  void WriteXML (string fileName);
  void WriteABINIT (string fileName="");
  void WriteASCII();
  NLPPClass() : XC(XC_NONE), Relativistic(false) {
   SetupMaps();
  }
};

#endif
