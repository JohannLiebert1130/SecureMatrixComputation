#include <fstream>
#include <iostream>
using namespace std;

const int MAXN = 33000;
bool p[MAXN] = {false};     //如果ｉ为素数则p[i] 为　false 

void SavePrime(ofstream& file)
{
    if (file.is_open())
    {
        for(int i = 2; i < MAXN; i++)
        {
            if(p[i] == false)
            {
                file << i << "\n";
                for(int j = 2*i; j < MAXN; j += i)
                    p[j] = true;
            }
        }
        file.close();
    }
    else{
        cout << "Unable to open file";
    }    
}

int main()
{
    ofstream fout;
    fout.open("prime.txt");
    SavePrime(fout);
    return 0;
}