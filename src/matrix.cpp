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



int main()
{
    PTMatrix ptMatrix(3,2,false);
    ptMatrix.print();
    vector<vector<long> > matrix{{1,2,3},{4,5,6},{7,8,9}};
    PTMatrix ptMatrix2(matrix,false);
    ptMatrix2.print();
    ptMatrix2.initDiagonal();
    ptMatrix2.printDiagonal();

    vector<vector<long> > matrix2 = ptMatrix2.getDiagonalMatrix();
    for(int i=0;i < 3; i++)
        for(int j = 0; j < 3; j++)
            cout << matrix2[i][j];
    cout << endl;
    cout << ptMatrix2(1,1);
    cout << endl;
    cout << ptMatrix2(2,0);
    PTMatrix newMatrix = ptMatrix * ptMatrix2;
    newMatrix.print();
    return 0;
}