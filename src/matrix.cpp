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

#define VEC_SIZE 4

// Simple class to measure time for each method
class Timer
{
public:
    void start() { m_start = my_clock(); }
    void stop() { m_stop = my_clock(); }
    double elapsed_time() const {
        return m_stop - m_start;
    }

private:
    double m_start, m_stop;
    double my_clock() const {
        struct timeval tv;
        gettimeofday(&tv, NULL);
        return tv.tv_sec + tv.tv_usec * 1e-6;
    }
};


/* --------------------- PTMatrix (Plain Text Matrix) class --------------------------*/

vector<vector<long> > PTMatrix::Row2Diagonal(const vector<vector<long> >& rowMatrix)
{
    vector<vector<long> > diagonalMatrix = vector<vector<long> >(rowMatrix[0].size(), vector<long>(rowMatrix.size(),0));
    for(unsigned int i=0, sz1 = diagonalMatrix.size(); i < sz1; i++)
        for(unsigned int j=0, sz2 = diagonalMatrix[i].size(); j < sz2; j++)
            diagonalMatrix[i][j] = rowMatrix[j][(i+j)%sz1];
    return diagonalMatrix;
}

vector<vector<long> > PTMatrix::Diagonal2Row(const vector<vector<long> >& diagonalMatrix)
{
    vector<vector<long> > rowMatrix = vector<vector<long> >(diagonalMatrix[0].size(), vector<long>(diagonalMatrix.size(),0));
    for(unsigned int i=0, sz1 = rowMatrix.size(); i < sz1; i++)
        for(unsigned int j=0, sz2 = rowMatrix[i].size(); j < sz2; j++)
            rowMatrix[i][j] = diagonalMatrix[myModulu(j-i,sz1)][i];
    return rowMatrix;
}

PTMatrix::PTMatrix(vector<vector<long> > matrix, bool isDiagonal, bool saveRowMatrix)
{   
    _d = matrix.size();
    if(isDiagonal)
    {
        _diagonalMatrix = matrix;
        if(saveRowMatrix)
            _rowMatrix = PTMatrix::Diagonal2Row(matrix);
    }
    else
    {
        _diagonalMatrix = PTMatrix::Row2Diagonal(matrix);
        if(saveRowMatrix)
            _rowMatrix = matrix;
    }
        
}

PTMatrix::PTMatrix(unsigned int dimension, int numbersLimit, bool saveRowMatrix)
/*this constructor create a random matrix
params:
 sizes: - the size of the matrix
 numbersLimit = the values limit, that means that all the values in the matrix would be between 0 to numbersLimit (not included)
*/
{
    _d = dimension;
    _diagonalMatrix = vector<vector<long> >(dimension, vector<long>(dimension));
    for(unsigned int i=0; i < dimension; i++)
        for(unsigned int j=0; j< dimension; j++)
            _diagonalMatrix[i][j] = rand() % numbersLimit;

    if(saveRowMatrix)
        _rowMatrix = PTMatrix::Diagonal2Row(_diagonalMatrix);
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
    unsigned int cellSize = largestNumInMatrixDigits(_diagonalMatrix)+1; //+1 for space
    
    cout << " ";
    
    for(unsigned int i=0; i< _d*cellSize; i++)
        cout << "-";
    cout << endl;
    for(unsigned int i=0; i< _d; i++)
    {
        cout << "|";
        for(unsigned int j=0; j< _d; j++)
                printNum((*this)(i,j), cellSize);
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
    if(row >= _d || column >= _d)
    {
        throw out_of_range("Error, indices out of bound!");
    }
    int i = row, j = column; //casting to int so the subtraction be ok
    return _diagonalMatrix[myModulu(j-i,_d)][row];
}

const long& PTMatrix::operator()(unsigned int row, unsigned int column) const{
    if(row >= _d || column >= _d)
    {
        throw out_of_range("Error, indices out of bound!");
    }
    int i = row, j = column; //casting to int so the subtraction be ok
    return _diagonalMatrix[myModulu(j-i,_d)][row];
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

vector<long> PTMatrix::MatrixEncoding(const vector<vector<long> >& rowMatrix)
{
    vector<long> v;
    for(int i=0; i < rowMatrix.size(); i++)
    {
        v.insert(v.end(), rowMatrix[i].begin(), rowMatrix[i].end());
    }
    return v;
}
EncryptedMatrix PTMatrix::encrypt(const EncryptedArray& ea, const FHEPubKey& publicKey, bool saveRow) const{
    vector<Ctxt> encDiagonalMatrix(_d, Ctxt(publicKey));
    Ctxt encRowMatrix(publicKey);

    unsigned int nslots = ea.size();
    for(unsigned int i=0; i< _d; i++)
        {
            vector<long> temp = _diagonalMatrix[i];
            temp.resize(nslots,0);
            ea.encrypt(encDiagonalMatrix[i], publicKey, temp);
        }
    
    if (saveRow)
    {
        vector<long> flatMatrix = PTMatrix::MatrixEncoding(_rowMatrix);
        flatMatrix.resize(nslots, 0);
        ea.encrypt(encRowMatrix, publicKey, flatMatrix);
    }

    return EncryptedMatrix(encDiagonalMatrix, encRowMatrix, _d, saveRow);
}

EncryptedMatrix PTMatrix::encrypt(const FHEPubKey& publicKey, bool saveRow) const{
    EncryptedArray ea(publicKey.getContext());
    return encrypt(ea, publicKey, saveRow);
}

PTMatrix PTMatrix::sigmaPermutation(unsigned int d){
    long dSquare = d * d;
    vector<vector<long> > rowMatrix(dSquare, vector<long>(dSquare, 0));
    for(int l=0; l < dSquare; l++)
        for(int i=0; i < d; i++)
            for(int j=0; j < d; j++)
            {
                if(l == d*i + myModulu(i+j,d)) rowMatrix[d*i+j][l] = 1;
            }
    return PTMatrix(rowMatrix, false, false);
}
PTMatrix PTMatrix::tauPermutation(unsigned int d){
    long dSquare = d * d;
    vector<vector<long> > rowMatrix(dSquare, vector<long>(dSquare, 0));
    for(int l=0; l < dSquare; l++)
        for(int i=0; i < d; i++)
            for(int j=0; j < d; j++)
            {
                if(l == d*myModulu(i+j,d) + j) rowMatrix[d*i+j][l] = 1;
            }
    return PTMatrix(rowMatrix, false, false);
}
PTMatrix PTMatrix::phiPermutation(unsigned int d, int k){
    long dSquare = d * d;
    vector<vector<long> > rowMatrix(dSquare, vector<long>(dSquare, 0));
    for(int l=0; l < dSquare; l++)
        for(int i=0; i < d; i++)
            for(int j=0; j < d; j++)
            {
                if(l == d*i + myModulu(k+j,d)) rowMatrix[d*i+j][l] = 1;
            }
    return PTMatrix(rowMatrix, false, false);
}
PTMatrix PTMatrix::psiPermutation(unsigned int d, int k){
    long dSquare = d * d;
    vector<vector<long> > rowMatrix(dSquare, vector<long>(dSquare, 0));
    for(int l=0; l < dSquare; l++)
        for(int i=0; i < d; i++)
            for(int j=0; j < d; j++)
            {
                if(l == d*myModulu(i+k,d) + j) rowMatrix[d*i+j][l] = 1;
            }
    return PTMatrix(rowMatrix, false, false);
}

vector<ZZX> PTMatrix::DiagonalEncoding(const EncryptedArray& ea)
{
    vector<ZZX> zzxVec(_d);
    for(int i = 0; i < _d; i++)
    {
        vector<long> temp = _diagonalMatrix[i];
        temp.resize(ea.size());
        ea.encode(zzxVec[i], temp);
    }
    return zzxVec;
}


/* --------------------- EncryptedMatrix class -------------*/
EncryptedMatrix::EncryptedMatrix(
                 const vector<Ctxt>& encDiagonalMatrix,
                 const Ctxt& encRowMatrix,
                 int dimension, bool haveRowMatrix): 
                 _rowMatrix(encRowMatrix),
                 _diagonalMatrix(encDiagonalMatrix),
                 _d(dimension),
                 _haveRowMatrix(haveRowMatrix){}

PTMatrix EncryptedMatrix::decrypt(const EncryptedArray& ea, const FHESecKey& secretKey) const {
    vector<vector<long> > diagonalMatrix(_d);
    for(unsigned int i=0; i < _d; i++){
        ea.decrypt(_diagonalMatrix[i], secretKey, diagonalMatrix[i]);
        diagonalMatrix[i].resize(_d, 0);
    }
    return PTMatrix(diagonalMatrix, true, false);
}

PTMatrix EncryptedMatrix::decrypt(const FHESecKey& secretKey) const {
    EncryptedArray ea(secretKey.getContext());
    return decrypt(ea, secretKey);
}

//matrix multyplication by vector. NOTE: this return a column vector! so don't use it to create a matrix (unless you want it to be column vectors matrix)
Ctxt EncryptedMatrix::operator*(const Ctxt& vec) const{
    EncryptedArray ea(vec.getContext());
    Ctxt result(vec.getPubKey());
    int len = _diagonalMatrix.size();
    
    //TODO: Still not perfectlly working
    
    Ctxt fixedVec = vec;
    cout << "ea.size " << ea.size() << " matrix row size " << _d << endl;
    if(ea.size() != _d) //Fix the problem that if the size of the vector is not nslots, the zero padding make the rotation push zeros to the begining of the vector
    {
        //replicate the vector to fill instead of zero padding
        for(unsigned int length = _d; length < ea.size(); length*=2){
            Ctxt copyVec = fixedVec;
            ea.shift(copyVec, length);  //shift length to right
            fixedVec+=copyVec;
        }
    }
    
    for(int i=0; i < len; i++)
    {
        Ctxt rotatedVec(fixedVec);   //copy vec
        ea.rotate(rotatedVec, -i);   //rotate it i right (-i left)
        rotatedVec *= _diagonalMatrix[i];
        result += rotatedVec;
    }
    return result;
}

Ctxt EncryptedMatrix::LinTrans1(const vector<ZZX>& matrix) const{
    EncryptedArray ea(_diagonalMatrix[0].getContext());
    Ctxt result(_diagonalMatrix[0].getPubKey());
    int len = matrix.size();
    
    //TODO: Still not perfectlly working
    
    Ctxt fixedVec = _rowMatrix;

    Timer shiftOperation;
    shiftOperation.start();
    if(ea.size() != len) //Fix the problem that if the size of the vector is not nslots, the zero padding make the rotation push zeros to the begining of the vector
    {
        //replicate the vector to fill instead of zero padding
        for(unsigned int length =len; length < ea.size(); length*=2){
            Ctxt copyVec = fixedVec;
            ea.shift(copyVec, length);  //shift length to right
            fixedVec+=copyVec;
        }
    }
    shiftOperation.stop();
    std::cout << "Time taken for the shiftOperation: " << shiftOperation.elapsed_time() << std::endl;

    Timer lintrans1Inner;
    lintrans1Inner.start();
    for(int i=-_d+1; i < _d; i++)
    {
        Ctxt rotatedVec(fixedVec);   //copy vec
        ea.rotate(rotatedVec, -i);   //rotate it i right (-i left)
        rotatedVec.multByConstant(matrix[myModulu(i, len)]);
        result += rotatedVec;
    }
    lintrans1Inner.stop();
    std::cout << "Time taken for the lintrans1Inner: " << lintrans1Inner.elapsed_time() << std::endl;

    return result;
}

Ctxt EncryptedMatrix::LinTrans2(const vector<ZZX>& matrix) const{
    EncryptedArray ea(_diagonalMatrix[0].getContext());
    Ctxt result(_diagonalMatrix[0].getPubKey());
    int len = matrix.size();
    
    //TODO: Still not perfectlly working
    
    Ctxt fixedVec = _rowMatrix;
    if(ea.size() != len) //Fix the problem that if the size of the vector is not nslots, the zero padding make the rotation push zeros to the begining of the vector
    {
        //replicate the vector to fill instead of zero padding
        for(unsigned int length =len; length < ea.size(); length*=2){
            Ctxt copyVec = fixedVec;
            ea.shift(copyVec, length);  //shift length to right
            fixedVec+=copyVec;
        }
    }
    
    for(int i=0; i < _d; i++)
    {
        Ctxt rotatedVec(fixedVec);   //copy vec
        ea.rotate(rotatedVec, -(_d*i));   //rotate it i right (-i left)
        rotatedVec.multByConstant(matrix[_d*i]);
        result += rotatedVec;
    }

    return result;
}

Ctxt EncryptedMatrix::LinTrans3(const Ctxt& vec, const vector<ZZX>& matrix, int d, int k){
    EncryptedArray ea(vec.getContext());
    Ctxt result(vec.getPubKey());
    int len = matrix.size();
    
    //TODO: Still not perfectlly working
    
    Ctxt fixedVec = vec;
    if(ea.size() != len) //Fix the problem that if the size of the vector is not nslots, the zero padding make the rotation push zeros to the begining of the vector
    {
        //replicate the vector to fill instead of zero padding
        for(unsigned int length =len; length < ea.size(); length*=2){
            Ctxt copyVec = fixedVec;
            ea.shift(copyVec, length);  //shift length to right
            fixedVec+=copyVec;
        }
    }
    

    Ctxt rotatedVec(fixedVec);   //copy vec
    ea.rotate(rotatedVec, -k);   //rotate it i right (-i left)
    rotatedVec.multByConstant(matrix[k]);
    result = rotatedVec;


    rotatedVec =fixedVec;   //copy vec
    ea.rotate(rotatedVec, d-k);   //rotate it i right (-i left)
    rotatedVec.multByConstant(matrix[myModulu(k-d, len)]);

    result += rotatedVec;
    return result;
}

Ctxt EncryptedMatrix::LinTrans4(const Ctxt& vec, int d, int k){
    EncryptedArray ea(vec.getContext());
    Ctxt result(vec.getPubKey());
    int len = d * d;
    
    //TODO: Still not perfectlly working
    
    Ctxt fixedVec = vec;
    if(ea.size() != len) //Fix the problem that if the size of the vector is not nslots, the zero padding make the rotation push zeros to the begining of the vector
    {
        //replicate the vector to fill instead of zero padding
        for(unsigned int length =len; length < ea.size(); length*=2){
            Ctxt copyVec = fixedVec;
            ea.shift(copyVec, length);  //shift length to right
            fixedVec+=copyVec;
        }
    }
    

    Ctxt rotatedVec(fixedVec);   //copy vec
    ea.rotate(rotatedVec, -(d*k));   //rotate it i right (-i left)
    result = rotatedVec;

    return result;
}

int EncryptedMatrix::getD()
{
    return _d;
}

Ctxt EncryptedMatrix::getRowMatrix()
{
    return _rowMatrix;
}

Ctxt EncryptedMatrix::operator*( EncryptedMatrix& other) 
{
    //check sizes
    if(_d != other.getD())
        throw MatricesSizesNotMatch(_d, other.getD());

    Ctxt vec = _diagonalMatrix[0]; //save it for faster vec.getPubKey() in the loop
    EncryptedArray ea(vec.getContext());

    vector<Ctxt> A(_d, Ctxt(vec.getPubKey())), B(_d, Ctxt(vec.getPubKey()));

    Timer a0Init;
	a0Init.start();
    Timer sigma;
    sigma.start();
    vector<ZZX> sigmaMatrix = PTMatrix::sigmaPermutation(_d).DiagonalEncoding(ea);
    sigma.stop();
    std::cout << "Time taken for the sigmaMatrix: " << sigma.elapsed_time() << std::endl;

    Timer lintrans1;
    lintrans1.start();
    A[0] = LinTrans1(sigmaMatrix);
    lintrans1.stop();
    std::cout << "Time taken for the LinTrans1: " << lintrans1.elapsed_time() << std::endl;

    //A[0] = PTMatrix::sigmaPermutation(_d).encrypt(vec.getPubKey()).LinTrans1(getRowMatrix(), _d);
    a0Init.stop();
	std::cout << "Time taken for the a[0]: " << a0Init.elapsed_time() << std::endl;

    Timer b0Init;
	b0Init.start();
    B[0] = other.LinTrans2(PTMatrix::tauPermutation(_d).DiagonalEncoding(ea));
    b0Init.stop();
	std::cout << "Time taken for the b[0]: " << b0Init.elapsed_time() << std::endl;

   for(int k = 1; k < _d; k++)
   {
        Timer akInit;
	    akInit.start();
        A[k] = LinTrans3(A[0], PTMatrix::phiPermutation(_d, k).DiagonalEncoding(ea), _d, k);
        akInit.stop();
	    std::cout << "Time taken for the a[" << k <<"]: " << akInit.elapsed_time() << std::endl;

         Timer bkInit;
	    bkInit.start();
        B[k] = LinTrans4(B[0], _d, k);
        bkInit.stop();
	    std::cout << "Time taken for the b[" << k <<"]: " << bkInit.elapsed_time() << std::endl;
   }

    Timer multiplySum;
	multiplySum.start();
   Ctxt result = A[0];
   result.multiplyBy(B[0]) ;
   for(int k = 1; k < _d; k++)
   {
       Ctxt temp = A[k];
       temp *= B[k];
       result += temp;
   }
   multiplySum.stop();
	std::cout << "Time taken for multiplySum " << multiplySum.elapsed_time() << std::endl;
    return result;
}


int main()
{
    long m = 0;                   // Specific modulus
	long p = 113;                 // Plaintext base [default=2], should be a prime number
	long r = 1;                   // Lifting [default=1]
	long L = 390;                  // Number of levels in the modulus chain [default=heuristic]
	long c = 3;                   // Number of columns in key-switching matrix [default=2]
	long w = 64;                  // Hamming weight of secret key
	long d = 0;                   // Degree of the field extension [default=1]
	long k = 80;                  // Security parameter [default=80] 
    long s = 0;                   // Minimum number of slots [default=0]
    
    Timer tInit;
	tInit.start();
    m = FindM(k, L, c, p, d, s, 0);           // Find a value for m given the specified values
    
    std::cout << "Initializing context... " << std::flush;
	FHEcontext context(m, p, r); 	          // Initialize context
	buildModChain(context, L, c);             // Modify the context, adding primes to the modulus chain
	std::cout << "OK!" << std::endl;

	std::cout << "Generating keys... " << std::flush;
	FHESecKey secretKey(context);                    // Construct a secret key structure
	const FHEPubKey& publicKey = secretKey;                 // An "upcast": FHESecKey is a subclass of FHEPubKey
	secretKey.GenSecKey(w);                          // Actually generate a secret key with Hamming weight
	addSome1DMatrices(secretKey);                    // Extra information for relinearization

    ZZX G;
     //G =  context.alMod.getFactorsOverZZ()[0]; 
    EncryptedArray ea(context, G);

    tInit.stop();
	std::cout << "Time taken for the initialization: " << tInit.elapsed_time() << std::endl;
    cout << "m:" << m << endl;
    cout << "nslots: " << ea.size() << endl;

    int dimension = 4;

    Timer ptMatrixInit;
	ptMatrixInit.start();
    PTMatrix ptMatrix1(dimension,3, true);
    EncryptedMatrix encMatrix1 = ptMatrix1.encrypt(publicKey, true);

    
    PTMatrix ptMatrix2(dimension, 3, true);
    EncryptedMatrix encMatrix2 = ptMatrix2.encrypt(publicKey, true);
    ptMatrixInit.stop();
    std::cout << "Time taken for the ptMatrix initialization: " << ptMatrixInit.elapsed_time() << std::endl;
    ptMatrix1.print();
    ptMatrix2.print();

    PTMatrix ptResult = ptMatrix1 * ptMatrix2;
    cout << "matrix multiplication result:" << endl;
    ptResult.print();

    Ctxt result = encMatrix1 * encMatrix2;
    vector<long> temp(ea.size(), 0);
    ea.decrypt(result, secretKey, temp);
    for(int i = 0; i < dimension * dimension; i++)
        cout << temp[i] << ' ';
    cout << endl;


    return 0;
}