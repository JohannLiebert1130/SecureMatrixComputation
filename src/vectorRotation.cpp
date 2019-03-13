#include <algorithm>
#include <iostream>
#include <vector>
using namespace std;

int main()
{
    vector<int> vec1{1,2,3,4,5,6,7,8,9};
    std::cout << "Old vector :"; 
    for(int i=0; i < vec1.size(); i++) 
        std::cout << " " << vec1[i]; 
    std::cout << "\n"; 

    int rotL = 3;
    rotate(vec1.begin(), vec1.begin() + rotL, vec1.end());
    std::cout << "New vector after left rotation :"; 
    for (int i=0; i < vec1.size(); i++) 
        std::cout<<" "<<vec1[i]; 
    std::cout << "\n\n"; 

    return 0;
}