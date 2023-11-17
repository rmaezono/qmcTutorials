#ifndef PARSER_CLASS_H
#define PARSER_CLASS_H

#include<vector>
#include<string>
#include<complex>
#include<cstdio>

using namespace std;


class ParserClass {
 private:
  string Buffer;
  int Pos;
 public:
  bool ReadFile (string fname);
  bool FindToken (string token);
  bool ReadInt (int &val);
  bool ReadDouble(double &val);
  bool ReadComplex(complex<double> &val);
  bool ReadWord (string &word);
  inline void Reset() { Pos = 0; }
  ParserClass() { /* do nothing for now */ }
};

#endif
