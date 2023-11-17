#ifndef PARSE_COMMAND_H
#define PARSE_COMMAND_H

#include <map>
#include <string>
#include <iostream>
#include <assert.h>
#include <list>
#include <vector>

using namespace std;


class ParamClass {
 private:
  string Arg, Name;
 public:
  bool NeedsArg, Found;
  string GetName () { return Name; }
  string GetArg  () { assert (NeedsArg); return Arg; }
  void SetArg (string arg) { assert (NeedsArg); Arg = arg; }
  ParamClass (string name, bool needsArg) {
   Name     = name;
   NeedsArg = needsArg;
   Found    = false;
  }
  ParamClass() {
   NeedsArg = false;
   Found    = false;
  }
};


class CommandLineParserClass {
 private:
  map<string, ParamClass> ArgMap;
  vector<string> Files;
 public:
  bool Parse (int argc, char **argv);
  inline bool Found (string name) { return ArgMap[name].Found; }
  inline string GetArg (string name) { return ArgMap[name].GetArg(); }
  inline int NumFiles() { return Files.size(); }
  string GetFile(int i) { return Files[i]; }
  CommandLineParserClass (list<ParamClass> &argList);
};

#endif
