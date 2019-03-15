#include <NTL/lzz_pXFactoring.h>
#include "FHE.h"
#include "EncryptedArray.h"
#include "replicate.h"
#include <fstream>
#include <sstream>
#include <sys/time.h>
#include <exception>

using namespace std;
using namespace NTL;

class PTMatrix;
class EncryptedMatrix;


class PTMatrix{
private:
    vector<vector<long> > _rowMatrix;
    vector<vector<long> > _diagonalMatrix;
    unsigned int _d;  // matrix dimension

public:
    PTMatrix(vector<vector<long> > matrix, bool isDiagonal, bool saveRowMatrix=false);
    PTMatrix(unsigned int dimension, int numbersLimit, bool saveRowMatrix=false);   //random matrix

    int getDimension() const;
    static vector<vector<long> > Row2Diagonal(const vector<vector<long> >& rowMatrix);    // transform a matrix from row order to diagonal order
    static vector<vector<long> > Diagonal2Row(const vector<vector<long> >& diagonalMatrix);
    static vector<long> MatrixEncoding(const vector<vector<long> >& rowMatrix);


    vector<vector<long> > getRowMatrix() const;
    vector<vector<long> > getDiagonalMatrix() const;

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
    
    bool isEqual(const vector<long>& other) const;
    
    bool operator==(const PTMatrix& other) const;
    bool operator!=(const PTMatrix& other) const;
    
    //operator%, apply %p on each element in the matrix. Use for checking if the encrypted multiplication's result is right, because it calculated under some modulo field
    PTMatrix operator%(unsigned int p) const;
    PTMatrix operator%=(unsigned int p);
    
    long& operator()(unsigned int row, unsigned int column); //return the element in [i,j] if it was regular matrix
    const long& operator()(unsigned int row, unsigned int column) const;

    //Encrypting
    EncryptedMatrix encrypt(const EncryptedArray& ea, const FHEPubKey& publicKey, bool saveRow=false) const;
    EncryptedMatrix encrypt(const FHEPubKey& publicKey, bool saveRow=false) const;

    //matrices multiplication
    Ctxt operator*(const EncryptedMatrix& other) const;
    Ctxt operator*=(const EncryptedMatrix& other);

    static PTMatrix sigmaPermutation(unsigned int d);
    static PTMatrix tauPermutation(unsigned int d);
    static PTMatrix phiPermutation(unsigned int d, int k);
    static PTMatrix psiPermutation(unsigned int d, int k);

    vector<ZZX> DiagonalEncoding(const EncryptedArray& ea);

};

//This class represents an encrypted matrix.
class EncryptedMatrix{
private:
    vector<Ctxt> _diagonalMatrix;
    Ctxt _rowMatrix;
    bool _haveRowMatrix;
    int _d;
public:
    EncryptedMatrix(const vector<Ctxt>& encDiagonalMatrix,
                    const Ctxt& encRowMatrix,
                    int dimension, bool haveRowMatrix);

    int getD();
    Ctxt getRowMatrix();
    //Decrypt
    PTMatrix decrypt(const EncryptedArray& ea, const FHESecKey& secretKey) const;
    PTMatrix decrypt(const FHESecKey& secretKey) const;

    //matrix multyplication by vector
    Ctxt operator*(const Ctxt& vec) const;

    Ctxt LinTrans1(const vector<ZZX>& matrix) const;
    Ctxt LinTrans2(const vector<ZZX>& matrix) const;
    static Ctxt LinTrans3(const Ctxt& vec, const vector<ZZX>& matrix, int d, int k);
    static Ctxt LinTrans4(const Ctxt& vec, int d, int k);

    // Matrix Encoding Method,
    static Ctxt MatrixEncoding(vector<Ctxt> matrix);

    //matrices multiplication
    Ctxt operator*( EncryptedMatrix& other) ;
    Ctxt operator*=( EncryptedMatrix& other);


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