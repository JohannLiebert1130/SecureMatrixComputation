#include <NTL/lzz_pXFactoring.h>
#include<vector>
#include <fstream>
#include <sstream>
#include <sys/time.h>
#include <exception>

using namespace std;

class PTMatrix{
private:
    vector<vector<long> > _rowMatrix;
    vector<vector<long> > _diagonalMatrix;
    unsigned int _d;  // matrix dimension

public:
    PTMatrix(vector<vector<long> > rowMatrix, bool saveDiagonal=false);
    PTMatrix(unsigned int dimension, int numbersLimit, bool saveDiagonal=false);   //random matrix

    int getDimension() const;
    static vector<vector<long> > row2Diagonal(const vector<vector<long> >& rowMatrix);    // transform a matrix from row order to diagonal order
    void initDiagonal();
    vector<vector<long> > getRowMatrix() const;
    vector<vector<long> > getDiagonalMatrix() const;

    //Encrypting
    //EncryptedMatrix encrypt(const EncryptedArray& ea, const FHEPubKey& publicKey, bool saveRow=false) const;
    //EncryptedMatrix encrypt(const FHEPubKey& publicKey, bool saveRow=false) const;

    void print(string label="") const;   //prints the row-order matrix with some comment/label
    void printDiagonal(string label="") const;

    //matrices multiplication
    PTMatrix operator*(const PTMatrix& other) const;
    PTMatrix operator*=(const PTMatrix& other);
    
    //Mul by constant
    PTMatrix operator*(unsigned int num) const;
    PTMatrix operator*=(unsigned int num);
    
    //matrices addition
    PTMatrix operator+(const PTMatrix& other) const;
    PTMatrix operator+=(const PTMatrix& other);
    
    //matrices substruction
    PTMatrix operator-(const PTMatrix& other) const;
    PTMatrix operator-=(const PTMatrix& other);
    
    //Transpose
    PTMatrix transpose() const;
    
    PTMatrix operator>(const PTMatrix& other) const;
    PTMatrix operator<(const PTMatrix& other) const;
    PTMatrix operator>=(const PTMatrix& other) const;
    PTMatrix operator<=(const PTMatrix& other) const;
    
    bool operator==(const PTMatrix& other) const;
    bool operator!=(const PTMatrix& other) const;
    
    //operator%, apply %p on each element in the matrix. Use for checking if the encrypted multiplication's result is right, because it calculated under some modulo field
    PTMatrix operator%(unsigned int p) const;
    PTMatrix operator%=(unsigned int p);
    
    long& operator()(unsigned int row, unsigned int column); //return the element in [i,j] if it was regular matrix
    const long& operator()(unsigned int row, unsigned int column) const;
};

//EXCEPTIONS
class MatricesSizesNotMatch:public runtime_error{
public:
    MatricesSizesNotMatch(const int d1, const int d2) : runtime_error("Matrices sizes are not match!"), sz1(d1), sz2(d2) {}
    const char* what(){
        cnvt.str( "" );
        cnvt << runtime_error::what() << ": First matrix size: " << sz1 << "x" << sz1 << ", second matrix size: " << sz2 << "x" << sz2;
        return cnvt.str().c_str();
    }
    ~MatricesSizesNotMatch() throw() {};
private:
    int sz1, sz2;
    static ostringstream cnvt;
};