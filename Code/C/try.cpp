
#include <iostream>
#include <vector>

using namespace std;

int main(void){
  
  vector<double> v(10, 1.0);

  vector<double>::iterator v_iter = v.begin();

  for(int i=0; i < 20; ++i){
    *v_iter = (double)i;
    cout << *v_iter;
    ++v_iter;
  }

  cout << v.size() << "\n";

  return 0;
}
