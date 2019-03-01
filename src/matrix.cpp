#include "matrix.h"

/* ----- some helping functions ------ */
long myModulu(long num, long mod)
//the regular modulo operator on minus numbers make trouble when trying to calculate enteries in the diagonal matrices.
//for example, -7%3=-1, but I want it to be 2 (-7+3=-4, -4+3=-1, -1+3= 2)
{
    //adding "mod" to "num" until it positive
    while(num<0)
        num+=mod;
    return num%mod;
}

unsigned int numDigit(long num){ return (num/10)==0 ? 1:numDigit(num/10)+1; } //How much digits in the input number. Uses for a "nice" printing

unsigned int largestNumInMatrixDigits(vector<vector<long> > matrix) //Find the longest number in a matrix and returns it length (in digits). Uses for a "nice" printing
{
    int largest = 0;
    long temp;
    for(unsigned int i=0; i < matrix.size(); i++)
        for(unsigned int j=0; j<matrix[i].size(); j++)
            if((temp = numDigit(matrix[i][j])) > largest)
                largest = temp;
    return largest;
}
void printNum(long num, int size){
    int len = numDigit(num);
    cout << num;
    for(int i=len; i < size; i++)
        cout << " ";
}

/* --------------------- PTMatrix (Plain Text Matrix) class --------------------------*/

vector<vector<long> > PTMatrix::row2Diagonal(const vector<vector<long> >& rowMatrix)
{
    vector<vector<long> > diagonalMatrix = vector<vector<long> >(rowMatrix[0].size(), vector<long>(rowMatrix.size(),0));
    for(unsigned int i=0, sz1 = diagonalMatrix.size(); i < sz1; i++)
        for(unsigned int j=0, sz2 = diagonalMatrix[i].size(); j < sz2; j++)
            diagonalMatrix[i][j] = rowMatrix[j][(i+j)%sz1];
    return diagonalMatrix;
}
PTMatrix::PTMatrix(vector<vector<long> > rowMatrix, bool saveDiagonal)
{
    _d = rowMatrix.size();
    _rowMatrix = rowMatrix;
    if(saveDiagonal)
        _diagonalMatrix = PTMatrix::row2Diagonal(_rowMatrix);
}

PTMatrix::PTMatrix(unsigned int dimension, int numbersLimit, bool saveDiagonal)
/*this constructor create a random matrix
params:
 sizes: - the size of the matrix
 numbersLimit = the values limit, that means that all the values in the matrix would be between 0 to numbersLimit (not included)
*/
{
    _d = dimension;
    _rowMatrix = vector<vector<long> >(dimension, vector<long>(dimension));
    for(unsigned int i=0; i < dimension; i++)
        for(unsigned int j=0; j< dimension; j++)
            _rowMatrix[i][j] = rand() % numbersLimit;

    if(saveDiagonal)
        _diagonalMatrix = PTMatrix::row2Diagonal(_rowMatrix);
}

void PTMatrix::initDiagonal()
{
    _diagonalMatrix = PTMatrix::row2Diagonal(_rowMatrix);
}

vector<vector<long> > PTMatrix::getRowMatrix() const
{
    return _rowMatrix;
}

vector<vector<long> > PTMatrix::getDiagonalMatrix() const
{
    return _diagonalMatrix;
}

void PTMatrix::print(string label) const{
    if(label.compare("")!=0)
        cout << label << endl;

    unsigned int cellSize = largestNumInMatrixDigits(_rowMatrix)+1; //+1 for space
    
    cout << " ";
    
    for(unsigned int i=0; i< _d*cellSize; i++)
        cout << "-";
    cout << endl;
    for(unsigned int i=0; i< _d; i++)
    {
        cout << "|";
        for(unsigned int j=0; j< _d; j++)
                printNum(_rowMatrix[i][j], cellSize);
        cout << "|" << endl;
    }
    cout << " ";
    for(unsigned int i=0; i< _d*cellSize; i++)
        cout << "-";
    cout << endl;
}

void PTMatrix::printDiagonal(string label) const
{
    if(label.compare("")!=0)
        cout << label << endl;

    unsigned int cellSize = largestNumInMatrixDigits(_diagonalMatrix)+1; //+1 for space
    
    cout << " ";
    
    for(unsigned int i=0; i< _d*cellSize; i++)
        cout << "-";
    cout << endl;
    for(unsigned int i=0; i< _d; i++)
    {
        cout << "|";
        for(unsigned int j=0; j< _d; j++)
                printNum(_diagonalMatrix[i][j], cellSize);
        cout << "|" << endl;
    }
    cout << " ";
    for(unsigned int i=0; i< _d*cellSize; i++)
        cout << "-";
    cout << endl;
}

int PTMatrix::getDimension() const
{
    return _d;
}

//operators

long& PTMatrix::operator()(unsigned int row, unsigned int column){
    if(row >= getDimension() || column >= getDimension())
    {
        cout << "Error, indices out of bound! MatSize: " << getDimension() << "*" << getDimension() <<", indices: " << row << "*" << column << endl;
        throw out_of_range("Error, indices out of bound!");
    }
    return _rowMatrix[row][column];
}

const long& PTMatrix::operator()(unsigned int row, unsigned int column) const{
    if(row >= getDimension() || column >= getDimension())
    {
        cout << "Error, indices out of bound! MatSize: " << getDimension() << "*" << getDimension() <<", indices: " << row << "*" << column << endl;
        throw out_of_range("Error, indices out of bound!");
    }
    return _rowMatrix[row][column];
}

PTMatrix PTMatrix::operator*(const PTMatrix& other) const{
    //check sizes
    if(getDimension() != other.getDimension())
        throw MatricesSizesNotMatch(getDimension(), other.getDimension());
    
    vector<vector<long> > res(getDimension(), vector<long>(getDimension(),0));
    for(unsigned int i=0; i < res.size(); i++)
        for(unsigned int j=0; j < res[i].size(); j++)
            for(unsigned int k = 0; k < other.getDimension(); k++)
                res[i][j] += (*this)(i,k)*other(k,j);
    
    return PTMatrix(res, false);
}

PTMatrix PTMatrix::operator*=(const PTMatrix& other){ return (*this) = (*this)*other; }

EncryptedMatrix PTMatrix::encrypt(const EncryptedArray& ea, const FHEPubKey& publicKey, bool saveDiagonal) const{
    vector<Ctxt> encDiagonalMatrix(_d, Ctxt(publicKey));
    vector<Ctxt> encRowMatrix(_d, Ctxt(publicKey));

    unsigned int nslots = ea.size();
    for(unsigned int i=0; i< _d; i++)
    {
        vector<long> temp = _rowMatrix[i];
        temp.resize(nslots,0);
        ea.encrypt(encRowMatrix[i], publicKey, temp);
    }
    
    if (saveDiagonal)
    {
        for(unsigned int i=0; i< _d; i++)
        {
            vector<long> temp = _diagonalMatrix[i];
            temp.resize(nslots,0);
            ea.encrypt(encDiagonalMatrix[i], publicKey, temp);
        }
    }

    return EncryptedMatrix(encRowMatrix, encDiagonalMatrix, _d, saveDiagonal);
}

EncryptedMatrix PTMatrix::encrypt(const FHEPubKey& publicKey, bool saveDiagonal) const{
    EncryptedArray ea(publicKey.getContext());
    return encrypt(ea, publicKey, saveDiagonal);
}


/* --------------------- EncryptedMatrix class -------------*/
EncryptedMatrix::EncryptedMatrix(const vector<Ctxt>& encRowMatrix,
                 const vector<Ctxt>& encDiagonalMatrix,
                 int dimension, bool haveDiagonalMatrix): 
                 _rowMatrix(encRowMatrix),
                 _diagonalMatrix(encDiagonalMatrix),
                 _d(dimension),
                 _haveDiagonalMatrix(haveDiagonalMatrix){}

PTMatrix EncryptedMatrix::decrypt(const EncryptedArray& ea, const FHESecKey& secretKey) const {
    vector<vector<long> > rowMatrix(_d);
    for(unsigned int i=0; i < _d; i++){
        ea.decrypt(_rowMatrix[i], secretKey, rowMatrix[i]);
        rowMatrix[i].resize(_d, 0);
    }
    return PTMatrix(rowMatrix, true);
}

PTMatrix EncryptedMatrix::decrypt(const FHESecKey& secretKey) const {
    EncryptedArray ea(secretKey.getContext());
    return decrypt(ea, secretKey);
}



int main()
{
    long m=0, p=257, r=3; // Native plaintext space                  // Computations will be 'modulo p'
    long L=10;          // Levels
    long c=2;           // Columns in key switching matrix
    long w=64;          // Hamming weight of secret key
    long d=0;
    long security = 128;
    ZZX G;
    m = FindM(security,L,c,p, d, 0, 0);

    FHEcontext context(m, p, r);
    // initialize context
    buildModChain(context, L, c);
    // modify the context, adding primes to the modulus chain
    FHESecKey secretKey(context);
    // construct a secret key structure
    const FHEPubKey& publicKey = secretKey;
    // an "upcast": FHESecKey is a subclass of FHEPubKey

    //if(0 == d)
    G = context.alMod.getFactorsOverZZ()[0];

    secretKey.GenSecKey(w);
    // actually generate a secret key with Hamming weight w

    addSome1DMatrices(secretKey);
    cout << "Generated key" << endl;


    EncryptedArray ea(context, G);
    cout << ea.size();
    vector<vector<long> > matrix{{1,2,3},{4,5,6},{7,8,9}};
    PTMatrix ptMatrix1(matrix,true);
    ptMatrix1.print();

    EncryptedMatrix encMatrix = ptMatrix1.encrypt(ea, publicKey, true);
    PTMatrix ptMatrix2 = encMatrix.decrypt(ea, secretKey);
    ptMatrix2.print();
    return 0;
}