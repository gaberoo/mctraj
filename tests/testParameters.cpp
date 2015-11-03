#include <iostream>
#include <fstream>
using namespace std;

#include "../Parameters.h"

int main() {
  ifstream in("testParameters.json");
  string json_input;

  in.seekg(0,ios::end);
  json_input.reserve(in.tellg());
  in.seekg(0,ios::beg);

  json_input.assign(istreambuf_iterator<char>(in), 
                    istreambuf_iterator<char>());

  rapidjson::Document jpars;
  try {
    jpars.Parse<0>(json_input.c_str());
    if (! jpars.IsObject()) throw 10;
  } catch (int e) {
    cerr << "Coudln't create JSON document:" << endl;
    cerr << json_input << endl;
  }

  Parameters p;

  cerr << "Reading parameters..." << flush;
  try {
    p.from_json(jpars);
  } catch (int e) {
    cerr << "Coudln't create JSON document:" << endl;
    cerr << json_input << endl;
  } catch (const char* str) {
    cerr << "JSON exception: " << str << endl;
  }
  cerr << "done" << endl;

  cout << p.nval("N") << " " << p.value("N") << endl;
  cout << p.nval("beta") << " " << p.value("beta") << endl;
  cout << p.nval("mu") << " " << p.value("mu") << endl;
  cout << p.nval("psi") << " " << p.value("psi") << endl;

  return 0;
}
