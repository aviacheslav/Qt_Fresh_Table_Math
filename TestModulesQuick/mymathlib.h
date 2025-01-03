#ifndef MYMATHLIB_H
#define MYMATHLIB_H

#include<math.h>
#include<vector>
#include"mylib.h"
#include"myarraylib.h"
//#include"array2dlib.h"
//#include"serviceclasses.h"

//class MyMathLib
//{
//public:
//    MyMathLib();
//};

//const double _PI=3.1415927;

const int MatrixTypeN_P2S=121;
const int MatrixTypeN_P2L=122;
const int MatrixTypeN_P2A=123;
const int MatrixTypeN_V2S=221;
const int MatrixTypeN_V2L=222;
const int MatrixTypeN_V2A=223;
const int MatrixTypeN_P1S=111;
const int MatrixTypeN_P1L=112;
const int MatrixTypeN_P1A=113;
const int MatrixTypeN_V1S=211;
const int MatrixTypeN_V1L=212;
const int MatrixTypeN_V1A=213;
const int MatrixTypeN_Default=MatrixTypeN_V2L;


class TItersPrecision{
public:
    double precision;
    int MaxQIters;
    int Priority_1QI2Eps3And4Or;
    TItersPrecision(double precision=1E-6, int MaxQIters=10, int Priority_1QI2Eps3And4Or=4);
    bool ContinieIterations(double precision, int QIters);
};

class MathApprox{
public:
    double precision;//=1E-9;
    MathApprox(double precision=0);
    bool ET(double x1, double x2);

    bool GT(double x1, double x2);
    bool LT(double x1, double x2);
    bool NE(double x1, double x2);
    bool GE(double x1, double x2);
    bool LE(double x1, double x2);
    bool isPositive(double x);
    bool isNegative(double x);
    bool isZero(double x);
    bool isNonZero(double x);
};

bool fApprET(double x1, double x2, double precision);
bool fApprGT(double x1, double x2, double precision);
bool fApprLT(double x1, double x2, double precision);
bool fApprNE(double x1, double x2, double precision);
bool fApprGE(double x1, double x2, double precision);
bool fApprLE(double x1, double x2, double precision);
bool fApprPositive(double x, double precision);
bool fApprNegative(double x, double precision);
bool fApprZero(double x, double precision);
bool fApprNonZero(double x, double precision);

struct TPosInSucc{
     bool isLess, isWithin, isGreater;
     int EqualN, lessN;
};

TPosInSucc PosInSucc1(double*X, int Q, double x, double eps=0);
TPosInSucc PosInSucc(std::vector<double>X, double x, double eps=0);

double LinInterp(double*X, double*Y, int Q, double x, double eps=0);
double LinInterp(std::vector<double>X, std::vector<double>Y, double x, double eps=0);
double LinInterp2D(std::vector<std::vector<double>>Y, std::vector<double>Xext, std::vector<double>Xine, double xe, double xi,  double eps=0);

double Integr_Trapez_SimpleByElementsArray(std::vector<double>X, std::vector<double>Y);
double Integr_Trapez_SimpleByElementsArray(double*X, double*Y, int Q);

std::vector<double>GetSectionsBounds(double LB, double HB, int Q, double SectL=0, int SectL1Q2And3Or4=2, bool QOfSectsNotSectsBounds=true);

//void CalcDigits(int num, int digits[], /*int&order=0,*/ int SysBase=10, bool AsInNum=false);
void CalcDigits(int num, int digits[], int&order, int SysBase=10, TValsShowHide*vsh=NULL);//, bool AsInNum=false);

double IntPower(double a, int b);
int NaturalPower(int a, int b);

QString ComplexNumberToString(double ReX, double ImX);

double MyArctan(double x, double y, double eps = 0);
//
void EquationQuadratic(std::vector<double>C, std::vector<double>&ReX, std::vector<double>&ImX);
void EquationQuadratic(double a, double b, double c, std::vector<double>&ReX, std::vector<double>&ImX);
void EquationQuadratic(double C[], std::vector<double>&ReX, std::vector<double>&ImX);
void EquationQuadratic(double BToA, double CToA, std::vector<double>&ReX, std::vector<double>&ImX);

void EquationCubic(std::vector<double>C, std::vector<double>&ReX, std::vector<double>&ImX);
void EquationCubic(double C[], std::vector<double>&ReX, std::vector<double>&ImX);
void EquationCubic(double c3, double c2, double c1, double c0, std::vector<double>&ReX, std::vector<double>&ImX);
void EquationCubic(double C2ToC3, double C1ToC3, double C0ToC3, std::vector<double>&ReX, std::vector<double>&ImX);

int EquationOf4thPower(std::vector<double>C, std::vector<double>&ReX, std::vector<double>&ImX);
int EquationOf4thPower(double C[], std::vector<double>&ReX, std::vector<double>&ImX);
int EquationOf4thPower(double C3To4, double C2To4, double C1To4, double C0To4, std::vector<double>&ReX, std::vector<double>&ImX);
int EquationOf4thPower(double c4, double c3, double c2, double c1, double c0, std::vector<double>&ReX, std::vector<double>&ReY);

class PolyEq{
    //std::vector<double> C, ReX, ImX;
public:
    std::vector<double> C, ReX, ImX;
    //
    PolyEq();
    PolyEq(std::vector<double>C);
    PolyEq(const PolyEq&obj);
    ~PolyEq();
    //
    void Assign(const PolyEq&obj);
    PolyEq operator =(const PolyEq&obj);
    //
    void SetEquationAndCalcSolution(std::vector<double>C);
    void StartGenerator(std::vector<double>ReX);
    //void StartSolution()
    std::vector<double> GetReX();
    std::vector<double> GetImX();
    std::vector<double> GetCoefs();
    void SetOrderAndClear(int order);
    int GetOrder();
    QString EquationToString();
    std::string EquationToStdString();
    QString SolutionToString();
    std::string SolutionToToStdString();
};

void LinearApproximation(double*X, double*Y, int Q, double&a, double&b);
void LinearApproximation(const std::vector<double>&X, const std::vector<double>&Y, double&a, double&b);

double arctg2(double x, double y);

//------------------------------------------------------------------------------------------------------

class ComplexNum{
    double re, im;
public:
    ComplexNum(double re=0, double im=0);
    ComplexNum(const ComplexNum&obj);
    void Assign(const ComplexNum&obj);
    ~ComplexNum();
    ComplexNum& operator = (const ComplexNum&obj);
    double SetDecart(double re=0, double im=0);
    double SetPolar(double fi=0, double r=0);
    double GetRe();
    double GetIm();
    double GetFi();
    double GetR();
    ComplexNum GetConjugated();
    void GetPolar(double&r, double&fi);
    void ChangeRe(double val);
    void ChangeIm(double val);
    void ChangeFi(double val);
    void ChangeR(double val);
    ComplexNum& operator + (const ComplexNum&obj);
    ComplexNum& operator - (const ComplexNum&obj);
    ComplexNum& operator * (const ComplexNum&obj);
    ComplexNum& operator / (const ComplexNum&obj);
    ComplexNum& operator * (double val);
    bool operator ==(const ComplexNum&obj);
    bool operator !=(const ComplexNum&obj);
    bool IsConjugated(const ComplexNum&obj);
};

//====================================================================================================================
// Matrices classes: 12 types: _P2S, _P2L, _P2A; _V2S, _V2L, _V2A, _  _P1S, _P1L, _P1A; _V1S, _V1L, _V1A
// P - opinter, V - std::vector, A - array (my user's class)
// S - simple - doesn't use my array library
// L - uses my library array functions, but not library array classes
// S or L are not strictly provided, was no time, this may be corrected, especcially SetRow, AddRow, InsRow, DelRow
// 2 - 2D, 1 - 1D
//====================================================================================================================


class Matrix_Prototype{
    public:
    //
    virtual ~Matrix_Prototype()=default;
    //
    virtual void SetNull()=0;
    virtual void construct()=0;
    virtual Matrix_Prototype* CreateMatrix()=0;
    virtual Matrix_Prototype* clone()=0;
    virtual int GetTypeN() const=0 ;
    //
    void SetOne(double val=0);
    //
    virtual void SetSize(int QLins=1, int QCols=1)=0;
    //
    virtual int GetQLines() const=0;
    virtual int GetQColumns()  const=0;
    //
    virtual void Set(double**y, int QLins=1, int QCols=1);//data cor here by SetSize et SetComonent, size each sui
    virtual void Set(std::vector<std::vector<double>>y);
    virtual void Set(double*y, int QLines,  int QColumns=1);
    virtual void Set(std::vector<double>y, int QColumns=1);
    virtual void Set(double x, double y, double z);
    virtual void Set(double y, int QLines, int QColumns=1);
    //
    virtual double GetComponent(int LineN, int ColN) const=0;
    virtual void SetComponent(double val, int LineN, int ColN)=0;
    //
    bool LineNBelongsHere(int LineN) const;
    bool ColNBelongsHere(int LineN) const;
    bool NsBelongHere(int LineN, int ColN) const;
    std::vector<double> GetArrayOfLineN(int LineN) const;
    std::vector<double> GetArrayOfColN(int ColN) const;
    std::vector<double> GetContentAsArray1DByLines() const;
    std::vector<std::vector<double>> GetContentAsArray2DByLines() const;
    //
    void SetX(double val);//Set Component N
    void SetY(double val);
    void SetZ(double val);
    double GetX();//Get Component N
    double GetY();
    double GetZ();
    //
    virtual void SetLine(int N, double*x, int Q, int whatN1=1, int QDfltValsBefore=0, bool DefaulsAreZerosNotOwnVals=true)=0;
    virtual void SetColumn(int N, double*x, int Q, int whatN1=1, int QDfltValsBefore=0, bool DefaulsAreOwnValsNotZeros=true)=0;
    virtual void SetLine(int N, std::vector<double>x, int whatN1=1, int QDfltValsBefore=0, bool DefaulsAreZerosNotOwnVals=true)=0;
    virtual void SetColumn(int N, std::vector<double>x, int whatN1=1, int QDfltValsBefore=0, bool DefaulsAreOwnValsNotZeros=true)=0;
    //
    virtual void AddLine(double*x, int Q, int whatN1=1, int QDfltValsBefore=0)=0;
    virtual void AddColumn(double*x, int Q, int whatN1=1, int QDfltValsBefore=0)=0;
    virtual void AddLine(std::vector<double>x, int whatN1=1, int QDfltValsBefore=0)=0;
    virtual void AddColumn(std::vector<double>x, int whatN1=1, int QDfltValsBefore=0)=0;
    //
    virtual void InsLine(int N, double*x, int Q, int whatN1=1, int QDfltValsBefore=0)=0;
    virtual void InsColumn(int N, double*x, int Q, int whatN1=1, int QDfltValsBefore=0)=0;
    virtual void InsLine(int N, std::vector<double>x, int whatN1=1, int QDfltValsBefore=0)=0;
    virtual void InsColumn(int N, std::vector<double>x, int whatN1=1, int QDfltValsBefore=0)=0;
    //
    virtual void DelLine(int N=-1)=0;
    virtual void DelColumn(int N=-1)=0;
    //
    //double determinant() const;
    //double AlgSuppl(int LineN, int ColumnN) const;
    //
    void ConcatAdding(Matrix_Prototype*obj, int QShiftedCols=0, int FromN=1, int QAdded=0);
    void ConcatStretching(Matrix_Prototype*obj, int QShiftedLines=0, int FromN=1, int QAdded=0);
    //
    QString LineToString(int LineN, QString delim=", ", bool useBrackets=true);
    QString ColumnToString(int ColN, QString delim=", ", bool useBrackets=true);
    QString ToString(QString delimElem=", ", QString delimRow="; ", bool useBrackets=true);
    //QString GetLineAsStr(int LineN, QString delimElem=", ", bool useBrackets=true);
    //QString GetColumnAsStr(int ColN, QString delimElem=", ", bool useBrackets=true);
    //QString GetAllLinesAsStr(QString delimElem=", ", QString delimLines="; ",bool useBrackets=true);
    void StdCOut1D(QString delimElem=", ", QString delimRow="; ", bool useBrackets=true);
    void StdCOut2D(QString delimElem=", ", QString delimRow="; ", bool useBrackets=true);
};

class Matrix_P2S : public Matrix_Prototype{
    double**y;
    int QLins, QCols;
    public:
    Matrix_P2S(int QLines=1, int QColumns=1, double val=0);
    //Matrix_P2S(int QLins=1, int QCols=1);
    Matrix_P2S(double**y, int QLins=1, int QCols=1);
    Matrix_P2S(std::vector<std::vector<double>>y);
    Matrix_P2S(double*y, int QLines, int QColumns=1);
    Matrix_P2S(std::vector<double>y, int QColumns=1);
    Matrix_P2S(double x, double y, double z);
    //Matrix_P2S(int QLines, int QColumns=1, double val=0);
    //
    virtual ~Matrix_P2S();
    //
    virtual void SetNull();
    virtual void construct();
    virtual Matrix_Prototype* CreateMatrix();
    virtual Matrix_Prototype* clone();
    virtual int GetTypeN() const;
    //
    virtual void SetSize(int QLins=1, int QCols=1);
    //
    virtual int GetQLines() const;
    virtual int GetQColumns()  const;
    //
    virtual void Set(double**y, int QLins=1, int QCols=1);
    virtual void Set(std::vector<std::vector<double>>y);
    virtual void Set(double*y, int QLines,  int QColumns=1);
    virtual void Set(std::vector<double>y, int QColumns=1);
    virtual void Set(double x, double y, double z);
    virtual void Set(double y, int QLines, int QColumns=1);
    //
    //virtual void Set(double**y, int QLins=1, int QCols=1)=0;
    //virtual void Set(std::vector<std::vector<double>>y)=0;
    //virtual void Set(double*y, int QLines,  int QColumns=1)=0;
    //virtual void Set(std::vector<double>y, int QColumns=1)=0;
    //virtual void Set(double x, double y, double z)=0;
    //virtual void Set(double y, int QLines, int QColumns=1)=0;
    //
    virtual double GetComponent(int LineN, int ColN) const;
    virtual void SetComponent(double val, int LineN, int ColN);
    //
    virtual void SetLine(int N, double*x, int L=0, int FromN=1, int QDefaultBefore=0, bool DefaultValsAreZerosNotOwn=true);
    virtual void SetColumn(int N, double*x, int L=0, int FromN=1, int QDefaultBefore=0, bool DefaultValsAreOwnNotZeros=false);
    virtual void SetLine(int N, std::vector<double>x, int FromN=1, int QDefaultBefore=0, bool DefaultValsAreZerosNotOwn=true);
    virtual void SetColumn(int N, std::vector<double>x, int FromN=1, int QDefaultBefore=0, bool DefaultValsAreOwnNotZeros=false);
    //
    virtual void AddLine(double*x, int Q, int whatN1=1, int QDfltValsBefore=0);
    virtual void AddColumn(double*x, int Q, int whatN1=1, int QDfltValsBefore=0);
    virtual void AddLine(std::vector<double>x, int whatN1=1, int QDfltValsBefore=0);
    virtual void AddColumn(std::vector<double>x, int whatN1=1, int QDfltValsBefore=0);
    //
    virtual void InsLine(int N, double*x, int L=0, int FromN=1, int QDefaultBefore=0);
    virtual void InsColumn(int N, double*x, int L=0, int FromN=1, int QDefaultBefore=0);
    virtual void InsLine(int N, std::vector<double>x, int whatN1=1, int QDfltValsBefore=0);
    virtual void InsColumn(int N, std::vector<double>x, int whatN1=1, int QDfltValsBefore=0);
    //
    virtual void DelLine(int N);
    virtual void DelColumn(int N);
};


class Matrix_P2L : public Matrix_Prototype{
    double**y;
    int QLins, QCols;
    public:
    Matrix_P2L(int QLines=1, int QColumns=1, double val=0);
    //Matrix_P2L(int QLins=1, int QCols=1);
    Matrix_P2L(double**y, int QLins=1, int QCols=1);
    Matrix_P2L(std::vector<std::vector<double>>y);
    Matrix_P2L(double*y, int QLines, int QColumns=1);
    Matrix_P2L(std::vector<double>y, int QColumns=1);
    Matrix_P2L(double x, double y, double z);
    //Matrix_P2L(int QLines, int QColumns=1, double val=0);
    //
    virtual ~Matrix_P2L();
    //
    virtual void SetNull();
    virtual void construct();
    virtual Matrix_Prototype* CreateMatrix();
    virtual Matrix_Prototype* clone();
    virtual int GetTypeN() const;
    //
    virtual void SetSize(int QLins=1, int QCols=1);
    //
    virtual int GetQLines() const;
    virtual int GetQColumns()  const;
    //
    virtual void Set(double**y, int QLins=1, int QCols=1);
    virtual void Set(std::vector<std::vector<double>>y);
    virtual void Set(double*y, int QLines,  int QColumns=1);
    virtual void Set(std::vector<double>y, int QColumns=1);
    virtual void Set(double x, double y, double z);
    virtual void Set(double y, int QLines, int QColumns=1);
    //
    virtual double GetComponent(int LineN, int ColN) const;
    virtual void SetComponent(double val, int LineN, int ColN);
    //
    virtual void SetLine(int N, double*x, int L=0, int FromN=1, int QDefaultBefore=0, bool DefaultValsAreZerosNotOwn=true);
    virtual void SetColumn(int N, double*x, int L=0, int FromN=1, int QDefaultBefore=0, bool DefaultValsAreOwnNotZeros=false);
    virtual void SetLine(int N, std::vector<double>x, int FromN=1, int QDefaultBefore=0, bool DefaultValsAreZerosNotOwn=true);
    virtual void SetColumn(int N, std::vector<double>x, int FromN=1, int QDefaultBefore=0, bool DefaultValsAreOwnNotZeros=false);
    //
    virtual void AddLine(double*x, int Q, int whatN1=1, int QDfltValsBefore=0);
    virtual void AddColumn(double*x, int Q, int whatN1=1, int QDfltValsBefore=0);
    virtual void AddLine(std::vector<double>x, int whatN1=1, int QDfltValsBefore=0);
    virtual void AddColumn(std::vector<double>x, int whatN1=1, int QDfltValsBefore=0);
    //
    virtual void InsLine(int N, double*x, int L=0, int FromN=1, int QDefaultBefore=0);
    virtual void InsColumn(int N, double*x, int L=0, int FromN=1, int QDefaultBefore=0);
    virtual void InsLine(int N, std::vector<double>x, int whatN1=1, int QDfltValsBefore=0);
    virtual void InsColumn(int N, std::vector<double>x, int whatN1=1, int QDfltValsBefore=0);
    //
    virtual void DelLine(int N);
    virtual void DelColumn(int N);
};

class Matrix_V2S : public Matrix_Prototype{
    std::vector<std::vector<double>>y;
    public:
    Matrix_V2S(int QLines=1, int QColumns=1, double val=0);
    //Matrix_V2S(int QLins=1, int QCols=1);
    Matrix_V2S(double**y, int QLins=1, int QCols=1);
    Matrix_V2S(std::vector<std::vector<double>>y);
    Matrix_V2S(double*y, int QLines, int QColumns=1);
    Matrix_V2S(std::vector<double>y, int QColumns=1);
    Matrix_V2S(double x, double y, double z);
    //Matrix_V2S(int QLines, int QColumns=1, double val=0);
    //
    virtual ~Matrix_V2S();
    //
    virtual void SetNull();
    virtual void construct();
    virtual Matrix_Prototype* CreateMatrix();
    virtual Matrix_Prototype* clone();
    virtual int GetTypeN() const;
    //
    virtual void SetSize(int QLins=1, int QCols=1);
    //
    virtual int GetQLines() const;
    virtual int GetQColumns()  const;
    //
    virtual void Set(double**y, int QLins=1, int QCols=1);
    virtual void Set(std::vector<std::vector<double>>y);
    virtual void Set(double*y, int QLines,  int QColumns=1);
    virtual void Set(std::vector<double>y, int QColumns=1);
    virtual void Set(double x, double y, double z);
    virtual void Set(double y, int QLines, int QColumns=1);
    //
    virtual double GetComponent(int LineN, int ColN) const;
    virtual void SetComponent(double val, int LineN, int ColN);
    //
    virtual void SetLine(int N, double*x, int L=0, int FromN=1, int QDefaultBefore=0, bool DefaultValsAreZerosNotOwn=true);
    virtual void SetColumn(int N, double*x, int L=0, int FromN=1, int QDefaultBefore=0, bool DefaultValsAreOwnNotZeros=false);
    virtual void SetLine(int N, std::vector<double>x, int FromN=1, int QDefaultBefore=0, bool DefaultValsAreZerosNotOwn=true);
    virtual void SetColumn(int N, std::vector<double>x, int FromN=1, int QDefaultBefore=0, bool DefaultValsAreOwnNotZeros=false);
    //
    virtual void AddLine(double*x, int Q, int whatN1=1, int QDfltValsBefore=0);
    virtual void AddColumn(double*x, int Q, int whatN1=1, int QDfltValsBefore=0);
    virtual void AddLine(std::vector<double>x, int whatN1=1, int QDfltValsBefore=0);
    virtual void AddColumn(std::vector<double>x, int whatN1=1, int QDfltValsBefore=0);
    //
    virtual void InsLine(int N, double*x, int L=0, int FromN=1, int QDefaultBefore=0);
    virtual void InsColumn(int N, double*x, int L=0, int FromN=1, int QDefaultBefore=0);
    virtual void InsLine(int N, std::vector<double>x, int whatN1=1, int QDfltValsBefore=0);
    virtual void InsColumn(int N, std::vector<double>x, int whatN1=1, int QDfltValsBefore=0);
    //
    virtual void DelLine(int N);
    virtual void DelColumn(int N);
};

class Matrix_V2L : public Matrix_Prototype{
    std::vector<std::vector<double>>y;
    public:
    Matrix_V2L(int QLines=1, int QColumns=1, double val=0);
    //Matrix_V2L(int QLins=1, int QCols=1);//AMBIGOUS CO PREV - COMMENTI NG s
    Matrix_V2L(double**y, int QLins=1, int QCols=1);
    Matrix_V2L(std::vector<std::vector<double>>y);
    Matrix_V2L(double*y, int QLines, int QColumns=1);
    Matrix_V2L(std::vector<double>y, int QColumns=1);
    Matrix_V2L(double x, double y, double z);
    //Matrix_V2L(int QLines, int QColumns=1, double val=0);
    //
    virtual ~Matrix_V2L();
    //
    virtual void SetNull();
    virtual void construct();
    virtual Matrix_Prototype* CreateMatrix();
    virtual Matrix_Prototype* clone();
    virtual int GetTypeN() const;
    //
    virtual void SetSize(int QLins=1, int QCols=1);
    //
    virtual int GetQLines() const;
    virtual int GetQColumns()  const;
    //
    virtual void Set(double**y, int QLins=1, int QCols=1);
    virtual void Set(std::vector<std::vector<double>>y);
    virtual void Set(double*y, int QLines,  int QColumns=1);
    virtual void Set(std::vector<double>y, int QColumns=1);
    virtual void Set(double x, double y, double z);
    virtual void Set(double y, int QLines, int QColumns=1);
    //
    virtual double GetComponent(int LineN, int ColN) const;
    virtual void SetComponent(double val, int LineN, int ColN);
    //
    virtual void SetLine(int N, double*x, int L=0, int FromN=1, int QDefaultBefore=0, bool DefaultValsAreZerosNotOwn=true);
    virtual void SetColumn(int N, double*x, int L=0, int FromN=1, int QDefaultBefore=0, bool DefaultValsAreOwnNotZeros=false);
    virtual void SetLine(int N, std::vector<double>x, int FromN=1, int QDefaultBefore=0, bool DefaultValsAreZerosNotOwn=true);
    virtual void SetColumn(int N, std::vector<double>x, int FromN=1, int QDefaultBefore=0, bool DefaultValsAreOwnNotZeros=false);
    //
    virtual void AddLine(double*x, int Q, int whatN1=1, int QDfltValsBefore=0);
    virtual void AddColumn(double*x, int Q, int whatN1=1, int QDfltValsBefore=0);
    virtual void AddLine(std::vector<double>x, int whatN1=1, int QDfltValsBefore=0);
    virtual void AddColumn(std::vector<double>x, int whatN1=1, int QDfltValsBefore=0);
    //
    virtual void InsLine(int N, double*x, int L=0, int FromN=1, int QDefaultBefore=0);
    virtual void InsColumn(int N, double*x, int L=0, int FromN=1, int QDefaultBefore=0);
    virtual void InsLine(int N, std::vector<double>x, int whatN1=1, int QDfltValsBefore=0);
    virtual void InsColumn(int N, std::vector<double>x, int whatN1=1, int QDfltValsBefore=0);
    //
    virtual void DelLine(int N);
    virtual void DelColumn(int N);
};

class Matrix_P2A : public Matrix_Prototype{
    Array2D_PS<double> data;
    public:
    Matrix_P2A(int QLines=1, int QColumns=1, double val=0);
    //Matrix_P2A(int QLins=1, int QCols=1);
    Matrix_P2A(double**y, int QLins=1, int QCols=1);
    Matrix_P2A(std::vector<std::vector<double>>y);
    Matrix_P2A(double*y, int QLines, int QColumns=1);
    Matrix_P2A(std::vector<double>y, int QColumns=1);
    Matrix_P2A(double x, double y, double z);
    //Matrix_P2A(int QLines, int QColumns=1, double val=0);
    //
    virtual ~Matrix_P2A();
    //
    virtual void SetNull();
    virtual void construct();
    virtual Matrix_Prototype* CreateMatrix();
    virtual Matrix_Prototype* clone();
    virtual int GetTypeN() const;
    //
    virtual void SetSize(int QLins=1, int QCols=1);
    //
    virtual int GetQLines() const;
    virtual int GetQColumns()  const;
    //
    virtual void Set(double**y, int QLins=1, int QCols=1);
    virtual void Set(std::vector<std::vector<double>>y);
    virtual void Set(double*y, int QLines,  int QColumns=1);
    virtual void Set(std::vector<double>y, int QColumns=1);
    virtual void Set(double x, double y, double z);
    virtual void Set(double y, int QLines, int QColumns=1);
    //
    virtual double GetComponent(int LineN, int ColN) const;
    virtual void SetComponent(double val, int LineN, int ColN);
    //
    virtual void SetLine(int N, double*x, int L=0, int FromN=1, int QDefaultBefore=0, bool DefaultValsAreZerosNotOwn=true);
    virtual void SetColumn(int N, double*x, int L=0, int FromN=1, int QDefaultBefore=0, bool DefaultValsAreOwnNotZeros=false);
    virtual void SetLine(int N, std::vector<double>x, int FromN=1, int QDefaultBefore=0, bool DefaultValsAreZerosNotOwn=true);
    virtual void SetColumn(int N, std::vector<double>x, int FromN=1, int QDefaultBefore=0, bool DefaultValsAreOwnNotZeros=false);
    //
    virtual void AddLine(double*x, int Q, int whatN1=1, int QDfltValsBefore=0);
    virtual void AddColumn(double*x, int Q, int whatN1=1, int QDfltValsBefore=0);
    virtual void AddLine(std::vector<double>x, int whatN1=1, int QDfltValsBefore=0);
    virtual void AddColumn(std::vector<double>x, int whatN1=1, int QDfltValsBefore=0);
    //
    virtual void InsLine(int N, double*x, int L=0, int FromN=1, int QDefaultBefore=0);
    virtual void InsColumn(int N, double*x, int L=0, int FromN=1, int QDefaultBefore=0);
    virtual void InsLine(int N, std::vector<double>x, int whatN1=1, int QDfltValsBefore=0);
    virtual void InsColumn(int N, std::vector<double>x, int whatN1=1, int QDfltValsBefore=0);
    //
    virtual void DelLine(int N);
    virtual void DelColumn(int N);
};

class Matrix_V2A : public Matrix_Prototype{
    Array2D_V<double> data;
    public:
    Matrix_V2A(int QLines=1, int QColumns=1, double val=0);
    //Matrix_V2A(int QLins=1, int QCols=1);
    Matrix_V2A(double**y, int QLins=1, int QCols=1);
    Matrix_V2A(std::vector<std::vector<double>>y);
    Matrix_V2A(double*y, int QLines, int QColumns=1);
    Matrix_V2A(std::vector<double>y, int QColumns=1);
    Matrix_V2A(double x, double y, double z);
    //Matrix_V2A(int QLines, int QColumns=1, double val=0);
    virtual ~Matrix_V2A();
    //
    virtual void SetNull();
    virtual void construct();
    virtual Matrix_Prototype* CreateMatrix();
    virtual Matrix_Prototype* clone();
    virtual int GetTypeN() const;
    //
    virtual void SetSize(int QLins=1, int QCols=1);
    //
    virtual int GetQLines() const;
    virtual int GetQColumns()  const;
    //
    virtual void Set(double**y, int QLins=1, int QCols=1);
    virtual void Set(std::vector<std::vector<double>>y);
    virtual void Set(double*y, int QLines,  int QColumns=1);
    virtual void Set(std::vector<double>y, int QColumns=1);
    virtual void Set(double x, double y, double z);
    virtual void Set(double y, int QLines, int QColumns=1);
    //
    virtual double GetComponent(int LineN, int ColN) const;
    virtual void SetComponent(double val, int LineN, int ColN);
    //
    virtual void SetLine(int N, double*x, int L=0, int FromN=1, int QDefaultBefore=0, bool DefaultValsAreZerosNotOwn=true);
    virtual void SetColumn(int N, double*x, int L=0, int FromN=1, int QDefaultBefore=0, bool DefaultValsAreOwnNotZeros=false);
    virtual void SetLine(int N, std::vector<double>x, int FromN=1, int QDefaultBefore=0, bool DefaultValsAreZerosNotOwn=true);
    virtual void SetColumn(int N, std::vector<double>x, int FromN=1, int QDefaultBefore=0, bool DefaultValsAreOwnNotZeros=false);
    //
    virtual void AddLine(double*x, int Q, int whatN1=1, int QDfltValsBefore=0);
    virtual void AddColumn(double*x, int Q, int whatN1=1, int QDfltValsBefore=0);
    virtual void AddLine(std::vector<double>x, int whatN1=1, int QDfltValsBefore=0);
    virtual void AddColumn(std::vector<double>x, int whatN1=1, int QDfltValsBefore=0);
    //
    virtual void InsLine(int N, double*x, int L=0, int FromN=1, int QDefaultBefore=0);
    virtual void InsColumn(int N, double*x, int L=0, int FromN=1, int QDefaultBefore=0);
    virtual void InsLine(int N, std::vector<double>x, int whatN1=1, int QDfltValsBefore=0);
    virtual void InsColumn(int N, std::vector<double>x, int whatN1=1, int QDfltValsBefore=0);
    //
    virtual void DelLine(int N);
    virtual void DelColumn(int N);
};

class Matrix_P1S : public Matrix_Prototype{
    double*y;
    int QLins, QCols;
    public:
    Matrix_P1S(int QLines=1, int QColumns=1, double val=0);
    //Matrix_P1S(int QLins=1, int QCols=1);
    Matrix_P1S(double**y, int QLins=1, int QCols=1);
    Matrix_P1S(std::vector<std::vector<double>>y);
    Matrix_P1S(double*y, int QLines, int QColumns=1);
    Matrix_P1S(std::vector<double>y, int QColumns=1);
    Matrix_P1S(double x, double y, double z);
    //Matrix_P1S(int QLines, int QColumns=1, double val=0);
    //
    virtual ~Matrix_P1S();
    //
    virtual void SetNull();
    virtual void construct();
    virtual Matrix_Prototype* CreateMatrix();
    virtual Matrix_Prototype* clone();
    virtual int GetTypeN() const;
    //
    virtual void SetSize(int QLins=1, int QCols=1);
    //
    virtual int GetQLines() const;
    virtual int GetQColumns()  const;
    //
    virtual void Set(double**y, int QLins=1, int QCols=1);
    virtual void Set(std::vector<std::vector<double>>y);
    virtual void Set(double*y, int QLines,  int QColumns=1);
    virtual void Set(std::vector<double>y, int QColumns=1);
    virtual void Set(double x, double y, double z);
    virtual void Set(double y, int QLines, int QColumns=1);
    //
    virtual double GetComponent(int LineN, int ColN) const;
    virtual void SetComponent(double val, int LineN, int ColN);
    //
    virtual void SetLine(int N, double*x, int L=0, int FromN=1, int QDefaultBefore=0, bool DefaultValsAreZerosNotOwn=true);
    virtual void SetColumn(int N, double*x, int L=0, int FromN=1, int QDefaultBefore=0, bool DefaultValsAreOwnNotZeros=false);
    virtual void SetLine(int N, std::vector<double>x, int FromN=1, int QDefaultBefore=0, bool DefaultValsAreZerosNotOwn=true);
    virtual void SetColumn(int N, std::vector<double>x, int FromN=1, int QDefaultBefore=0, bool DefaultValsAreOwnNotZeros=false);
    //
    virtual void AddLine(double*x, int Q, int whatN1=1, int QDfltValsBefore=0);
    virtual void AddColumn(double*x, int Q, int whatN1=1, int QDfltValsBefore=0);
    virtual void AddLine(std::vector<double>x, int whatN1=1, int QDfltValsBefore=0);
    virtual void AddColumn(std::vector<double>x, int whatN1=1, int QDfltValsBefore=0);
    //
    virtual void InsLine(int N, double*x, int L=0, int FromN=1, int QDefaultBefore=0);
    virtual void InsColumn(int N, double*x, int L=0, int FromN=1, int QDefaultBefore=0);
    virtual void InsLine(int N, std::vector<double>x, int whatN1=1, int QDfltValsBefore=0);
    virtual void InsColumn(int N, std::vector<double>x, int whatN1=1, int QDfltValsBefore=0);
    //
    virtual void DelLine(int N);
    virtual void DelColumn(int N);
    //

};

class Matrix_P1L : public Matrix_Prototype{
    double*y;
    int QLins, QCols;
    public:
    Matrix_P1L(int QLines=1, int QColumns=1, double val=0);
    //Matrix_P1L(int QLins=1, int QCols=1);
    Matrix_P1L(double**y, int QLins=1, int QCols=1);
    Matrix_P1L(std::vector<std::vector<double>>y);
    Matrix_P1L(double*y, int QLines, int QColumns=1);
    Matrix_P1L(std::vector<double>y, int QColumns=1);
    Matrix_P1L(double x, double y, double z);
   // Matrix_P1L(int QLines, int QColumns=1, double val=0);
    //
    virtual ~Matrix_P1L();
    //
    virtual void SetNull();
    virtual void construct();
    virtual Matrix_Prototype* CreateMatrix();
    virtual Matrix_Prototype* clone();
    virtual int GetTypeN() const;
    //
    virtual void SetSize(int QLins=1, int QCols=1);
    //
    virtual int GetQLines() const;
    virtual int GetQColumns()  const;
    //
    virtual void Set(double**y, int QLins=1, int QCols=1);
    virtual void Set(std::vector<std::vector<double>>y);
    virtual void Set(double*y, int QLines,  int QColumns=1);
    virtual void Set(std::vector<double>y, int QColumns=1);
    virtual void Set(double x, double y, double z);
    virtual void Set(double y, int QLines, int QColumns=1);
    //
    virtual double GetComponent(int LineN, int ColN) const;
    virtual void SetComponent(double val, int LineN, int ColN);
    //
    virtual void SetLine(int N, double*x, int L=0, int FromN=1, int QDefaultBefore=0, bool DefaultValsAreZerosNotOwn=true);
    virtual void SetColumn(int N, double*x, int L=0, int FromN=1, int QDefaultBefore=0, bool DefaultValsAreOwnNotZeros=false);
    virtual void SetLine(int N, std::vector<double>x, int FromN=1, int QDefaultBefore=0, bool DefaultValsAreZerosNotOwn=true);
    virtual void SetColumn(int N, std::vector<double>x, int FromN=1, int QDefaultBefore=0, bool DefaultValsAreOwnNotZeros=false);
    //
    virtual void AddLine(double*x, int Q, int whatN1=1, int QDfltValsBefore=0);
    virtual void AddColumn(double*x, int Q, int whatN1=1, int QDfltValsBefore=0);
    virtual void AddLine(std::vector<double>x, int whatN1=1, int QDfltValsBefore=0);
    virtual void AddColumn(std::vector<double>x, int whatN1=1, int QDfltValsBefore=0);
    //
    virtual void InsLine(int N, double*x, int L=0, int FromN=1, int QDefaultBefore=0);
    virtual void InsColumn(int N, double*x, int L=0, int FromN=1, int QDefaultBefore=0);
    virtual void InsLine(int N, std::vector<double>x, int whatN1=1, int QDfltValsBefore=0);
    virtual void InsColumn(int N, std::vector<double>x, int whatN1=1, int QDfltValsBefore=0);
    //
    virtual void DelLine(int N);
    virtual void DelColumn(int N);
    //
};

class Matrix_V1S : public Matrix_Prototype{
    std::vector<double>y;
    int QCols;
    public:
    Matrix_V1S(int QLines=1, int QColumns=1, double val=0);
    //Matrix_V1S(int QLins=1, int QCols=1);
    Matrix_V1S(double**y, int QLins=1, int QCols=1);
    Matrix_V1S(std::vector<std::vector<double>>y);
    Matrix_V1S(double*y, int QLines, int QColumns=1);
    Matrix_V1S(std::vector<double>y, int QColumns=1);
    Matrix_V1S(double x, double y, double z);
    //Matrix_V1S(int QLines, int QColumns=1, double val=0);
    //
    virtual ~Matrix_V1S();
    //
    virtual void SetNull();
    virtual void construct();
    virtual Matrix_Prototype* CreateMatrix();
    virtual Matrix_Prototype* clone();
    virtual int GetTypeN() const;
    //
    virtual void SetSize(int QLins=1, int QCols=1);
    //
    virtual int GetQLines() const;
    virtual int GetQColumns()  const;
    //
    virtual void Set(double**y, int QLins=1, int QCols=1);
    virtual void Set(std::vector<std::vector<double>>y);
    virtual void Set(double*y, int QLines,  int QColumns=1);
    virtual void Set(std::vector<double>y, int QColumns=1);
    virtual void Set(double x, double y, double z);
    virtual void Set(double y, int QLines, int QColumns=1);
    //
    virtual double GetComponent(int LineN, int ColN) const;
    virtual void SetComponent(double val, int LineN, int ColN);
    //
    virtual void SetLine(int N, double*x, int L=0, int FromN=1, int QDefaultBefore=0, bool DefaultValsAreZerosNotOwn=true);
    virtual void SetColumn(int N, double*x, int L=0, int FromN=1, int QDefaultBefore=0, bool DefaultValsAreOwnNotZeros=false);
    virtual void SetLine(int N, std::vector<double>x, int FromN=1, int QDefaultBefore=0, bool DefaultValsAreZerosNotOwn=true);
    virtual void SetColumn(int N, std::vector<double>x, int FromN=1, int QDefaultBefore=0, bool DefaultValsAreOwnNotZeros=false);
    //
    virtual void AddLine(double*x, int Q, int whatN1=1, int QDfltValsBefore=0);
    virtual void AddColumn(double*x, int Q, int whatN1=1, int QDfltValsBefore=0);
    virtual void AddLine(std::vector<double>x, int whatN1=1, int QDfltValsBefore=0);
    virtual void AddColumn(std::vector<double>x, int whatN1=1, int QDfltValsBefore=0);
    //
    virtual void InsLine(int N, double*x, int L=0, int FromN=1, int QDefaultBefore=0);
    virtual void InsColumn(int N, double*x, int L=0, int FromN=1, int QDefaultBefore=0);
    virtual void InsLine(int N, std::vector<double>x, int whatN1=1, int QDfltValsBefore=0);
    virtual void InsColumn(int N, std::vector<double>x, int whatN1=1, int QDfltValsBefore=0);
    //
    virtual void DelLine(int N);
    virtual void DelColumn(int N);
};

class Matrix_V1L : public Matrix_Prototype{
    std::vector<double>y;
    int QCols;
    public:
    Matrix_V1L(int QLines=1, int QColumns=1, double val=0);
    //Matrix_V1L(int QLins=1, int QCols=1);
    Matrix_V1L(double**y, int QLins=1, int QCols=1);
    Matrix_V1L(std::vector<std::vector<double>>y);
    Matrix_V1L(double*y, int QLines, int QColumns=1);
    Matrix_V1L(std::vector<double>y, int QColumns=1);
    Matrix_V1L(double x, double y, double z);
    //Matrix_V1L(int QLines, int QColumns=1, double val=0);
    //
    virtual ~Matrix_V1L();
    //
    virtual void SetNull();
    virtual void construct();
    virtual Matrix_Prototype* CreateMatrix();
    virtual Matrix_Prototype* clone();
    virtual int GetTypeN() const;
    //
    virtual void SetSize(int QLins=1, int QCols=1);
    //
    virtual int GetQLines() const;
    virtual int GetQColumns()  const;
    //
    virtual void Set(double**y, int QLins=1, int QCols=1);
    virtual void Set(std::vector<std::vector<double>>y);
    virtual void Set(double*y, int QLines,  int QColumns=1);
    virtual void Set(std::vector<double>y, int QColumns=1);
    virtual void Set(double x, double y, double z);
    virtual void Set(double y, int QLines, int QColumns=1);
    //
    virtual double GetComponent(int LineN, int ColN) const;
    virtual void SetComponent(double val, int LineN, int ColN);
    //
    virtual void SetLine(int N, double*x, int L=0, int FromN=1, int QDefaultBefore=0, bool DefaultValsAreZerosNotOwn=true);
    virtual void SetColumn(int N, double*x, int L=0, int FromN=1, int QDefaultBefore=0, bool DefaultValsAreOwnNotZeros=false);
    virtual void SetLine(int N, std::vector<double>x, int FromN=1, int QDefaultBefore=0, bool DefaultValsAreZerosNotOwn=true);
    virtual void SetColumn(int N, std::vector<double>x, int FromN=1, int QDefaultBefore=0, bool DefaultValsAreOwnNotZeros=false);
    //
    virtual void AddLine(double*x, int Q, int whatN1=1, int QDfltValsBefore=0);
    virtual void AddColumn(double*x, int Q, int whatN1=1, int QDfltValsBefore=0);
    virtual void AddLine(std::vector<double>x, int whatN1=1, int QDfltValsBefore=0);
    virtual void AddColumn(std::vector<double>x, int whatN1=1, int QDfltValsBefore=0);
    //
    virtual void InsLine(int N, double*x, int L=0, int FromN=1, int QDefaultBefore=0);
    virtual void InsColumn(int N, double*x, int L=0, int FromN=1, int QDefaultBefore=0);
    virtual void InsLine(int N, std::vector<double>x, int whatN1=1, int QDfltValsBefore=0);
    virtual void InsColumn(int N, std::vector<double>x, int whatN1=1, int QDfltValsBefore=0);
    //
    virtual void DelLine(int N);
    virtual void DelColumn(int N);
};

class Matrix_P1A : public Matrix_Prototype{
    //arrayOf2DBy1DPointerBased;
    public:
    Matrix_P1A(int QLines=1, int QColumns=1, double val=0);
    //Matrix_P1A(int QLins=1, int QCols=1);
    Matrix_P1A(double**y, int QLins=1, int QCols=1);
    Matrix_P1A(std::vector<std::vector<double>>y);
    Matrix_P1A(double*y, int QLines, int QColumns=1);
    Matrix_P1A(std::vector<double>y, int QColumns=1);
    Matrix_P1A(double x, double y, double z);
    //Matrix_P1A(int QLines, int QColumns=1, double val=0);
    //
    virtual ~Matrix_P1A();
    //
    virtual void SetNull();
    virtual void construct();
    virtual Matrix_Prototype* CreateMatrix();
    virtual Matrix_Prototype* clone();
    virtual int GetTypeN() const;
    //
    virtual void SetSize(int QLins=1, int QCols=1);
    //
    virtual int GetQLines() const;
    virtual int GetQColumns()  const;
    //
    virtual void Set(double**y, int QLins=1, int QCols=1);
    virtual void Set(std::vector<std::vector<double>>y);
    virtual void Set(double*y, int QLines,  int QColumns=1);
    virtual void Set(std::vector<double>y, int QColumns=1);
    virtual void Set(double x, double y, double z);
    virtual void Set(double y, int QLines, int QColumns=1);
    //
    virtual double GetComponent(int LineN, int ColN) const;
    virtual void SetComponent(double val, int LineN, int ColN);
    //
    virtual void SetLine(int N, double*x, int L=0, int FromN=1, int QDefaultBefore=0, bool DefaultValsAreZerosNotOwn=true);
    virtual void SetColumn(int N, double*x, int L=0, int FromN=1, int QDefaultBefore=0, bool DefaultValsAreOwnNotZeros=false);
    virtual void SetLine(int N, std::vector<double>x, int FromN=1, int QDefaultBefore=0, bool DefaultValsAreZerosNotOwn=true);
    virtual void SetColumn(int N, std::vector<double>x, int FromN=1, int QDefaultBefore=0, bool DefaultValsAreOwnNotZeros=false);
    //
    virtual void AddLine(double*x, int Q, int whatN1=1, int QDfltValsBefore=0);
    virtual void AddColumn(double*x, int Q, int whatN1=1, int QDfltValsBefore=0);
    virtual void AddLine(std::vector<double>x, int whatN1=1, int QDfltValsBefore=0);
    virtual void AddColumn(std::vector<double>x, int whatN1=1, int QDfltValsBefore=0);
    //
    virtual void InsLine(int N, double*x, int L=0, int FromN=1, int QDefaultBefore=0);
    virtual void InsColumn(int N, double*x, int L=0, int FromN=1, int QDefaultBefore=0);
    virtual void InsLine(int N, std::vector<double>x, int whatN1=1, int QDfltValsBefore=0);
    virtual void InsColumn(int N, std::vector<double>x, int whatN1=1, int QDfltValsBefore=0);
    //
    virtual void DelLine(int N);
    virtual void DelColumn(int N);
};

class Matrix_V1A : public Matrix_Prototype{
    //arrayOf2DBy1DStdVectorBased;
    public:
    Matrix_V1A(int QLines=1, int QColumns=1, double val=0);
    //Matrix_V1A(int QLins=1, int QCols=1);
    Matrix_V1A(double**y, int QLins=1, int QCols=1);
    Matrix_V1A(std::vector<std::vector<double>>y);
    Matrix_V1A(double*y, int QLines, int QColumns=1);
    Matrix_V1A(std::vector<double>y, int QColumns=1);
    Matrix_V1A(double x, double y, double z);
    //Matrix_V1A(int QLines, int QColumns=1, double val=0);
    //
    virtual ~Matrix_V1A();
    //
    virtual void SetNull();
    virtual void construct();
    virtual Matrix_Prototype* CreateMatrix();
    virtual Matrix_Prototype* clone();
    virtual int GetTypeN() const;
    //
    virtual void SetSize(int QLins=1, int QCols=1);
    //
    virtual int GetQLines() const;
    virtual int GetQColumns()  const;
    //
    virtual void Set(double**y, int QLins=1, int QCols=1);
    virtual void Set(std::vector<std::vector<double>>y);
    virtual void Set(double*y, int QLines,  int QColumns=1);
    virtual void Set(std::vector<double>y, int QColumns=1);
    virtual void Set(double x, double y, double z);
    virtual void Set(double y, int QLines, int QColumns=1);
    //
    virtual double GetComponent(int LineN, int ColN) const;
    virtual void SetComponent(double val, int LineN, int ColN);
    //
    virtual void SetLine(int N, double*x, int L=0, int FromN=1, int QDefaultBefore=0, bool DefaultValsAreZerosNotOwn=true);
    virtual void SetColumn(int N, double*x, int L=0, int FromN=1, int QDefaultBefore=0, bool DefaultValsAreOwnNotZeros=false);
    virtual void SetLine(int N, std::vector<double>x, int FromN=1, int QDefaultBefore=0, bool DefaultValsAreZerosNotOwn=true);
    virtual void SetColumn(int N, std::vector<double>x, int FromN=1, int QDefaultBefore=0, bool DefaultValsAreOwnNotZeros=false);
    //
    virtual void AddLine(double*x, int Q, int whatN1=1, int QDfltValsBefore=0);
    virtual void AddColumn(double*x, int Q, int whatN1=1, int QDfltValsBefore=0);
    virtual void AddLine(std::vector<double>x, int whatN1=1, int QDfltValsBefore=0);
    virtual void AddColumn(std::vector<double>x, int whatN1=1, int QDfltValsBefore=0);
    //
    virtual void InsLine(int N, double*x, int L=0, int FromN=1, int QDefaultBefore=0);
    virtual void InsColumn(int N, double*x, int L=0, int FromN=1, int QDefaultBefore=0);
    virtual void InsLine(int N, std::vector<double>x, int whatN1=1, int QDfltValsBefore=0);
    virtual void InsColumn(int N, std::vector<double>x, int whatN1=1, int QDfltValsBefore=0);
    //
    virtual void DelLine(int N);
    virtual void DelColumn(int N);
};

class Matrix_P2S_v1{
    double**y;
    int QLins, QCols;
    public:
    Matrix_P2S_v1(int QLins=1, int QCols=1);
    Matrix_P2S_v1(double**y, int QLins=1, int QCols=1);
    Matrix_P2S_v1(std::vector<std::vector<double>>y);
    Matrix_P2S_v1(double*y, int Q, bool LineNotColumn=true);
    Matrix_P2S_v1(std::vector<double>y, bool LineNotColumn=true);
    Matrix_P2S_v1(double x, double y, double z);
    Matrix_P2S_v1(const Matrix_P2S_v1&obj);
    ~Matrix_P2S_v1();
    void construct();
    void SetNull();
    void SetNullIni();
    void Assign(const Matrix_P2S_v1&obj);
    Matrix_P2S_v1& operator = (const Matrix_P2S_v1&obj);
    void SetSize(int QLins=1, int QCols=1);
    void Set(double**y, int QLins=1, int QCols=1);
    void Set(std::vector<std::vector<double>>y);
    void Set(std::vector<double>y, bool ColOfVectorNotHorLine=true);
    void Set(double*y, int Q, bool ColOfVectorNotHorLine=true);
    void Set(double x, double y, double z);
    int GetQLines() const;
    int GetQColumns()  const;
    double GetComponent(int LineN, int ColN) const;
    void SetComponent(double val, int LineN, int ColN);
    //
    Matrix_P2S_v1 operator *(double k);

    QString LineToString(int LineN, QString delim=", ", bool useBrackets=true);
    void StdCOut2D(QString delimElem=", ", QString delimRow="; ", bool useBrackets=true);

};

class Matrix{
    Matrix_Prototype*data;
public:
    Matrix(int QLines=1, int QColumns=1, double val=0, int TypeN=MatrixTypeN_Default);
    Matrix(double**y, int QLines, int QColumns=1, int TypeN=MatrixTypeN_Default);
    Matrix(double*y, int QLines, int QColumns=1, int TypeN=MatrixTypeN_Default);
    Matrix(std::vector<std::vector<double>>y, int QLines, int QColumns=1, int TypeN=MatrixTypeN_Default);
    Matrix(std::vector<double>y, int QColumns=1, int TypeN=MatrixTypeN_Default);
    Matrix(const Matrix&obj);
    ~Matrix();
    //
    void construct(int TypeN=MatrixTypeN_Default);
    void SetNull();
    //void SetNullIni();
    //
    void SetTypeN(int TypeN);
    int GetTypeN() const;
    //
    void SetOne(double val=0);
    //
    void Assign(const Matrix&obj);
    Matrix& operator = (const Matrix&obj);
    //
    void Set(int QLines, int QColumns=1, double val=0, int TypeN=MatrixTypeN_Default);
    void Set(double**y, int QLines, int QColumns=1, int TypeN=MatrixTypeN_Default);
    void Set(double*y, int QLines, int QColumns=1, int TypeN=MatrixTypeN_Default);
    void Set(std::vector<std::vector<double>>y, int TypeN=MatrixTypeN_Default);
    void Set(std::vector<double>, int QColumns=1, int TypeN=MatrixTypeN_Default);

    //
    int GetQLines() const;
    int GetQColumns() const;
    //
    void SetComponent(double val, int LineN, int ColN=1);
    double GetComponent(int LineN, int ColN=1) const;
    //
    void SetSize(int QLines, int QColumns);
    //
    void SetLineN(int rowN, double*x, int QValues, int whatN1=1, int QDfltsBefore=0, bool DefaulsAreZerosNotOwnVals=true);
    void SetLineN(int rowN, std::vector<double>x, int whatN1=1, int QDfltsBefore=0, bool DefaulsAreZerosNotOwnVals=true);
    //
    void SetColumnN(int rowN, double*x, int QValues, int whatN1, int QDfltsBefore=0, bool DefaulsAreZerosNotOwnVals=true);
    void SetColumnN(int rowN, std::vector<double>x, int whatN1, int QDfltsBefore=0, bool DefaulsAreZerosNotOwnVals=true);
   //
    void AddLine(double*x, int QValues, int whatN1=1, int QDfltsBefore=0);
    void AddLine(std::vector<double>x, int whatN1=1, int QDfltsBefore=0);
    //
    void AddColumn(double*x, int QValues, int whatN1=1, int QDfltsBefore=0);
    void AddColumn(std::vector<double>x, int whatN1=1, int QDfltsBefore=0);
    //
    void InsLineN(int rowN, double*x, int QValues, int whatN1=1, int QDfltsBefore=0);
    void InsLineN(int rowN, std::vector<double>x, int whatN1=1, int QDfltsBefore=0);
    //
    void InsColumnN(int rowN, double*x, int QValues, int whatN1=1, int QDfltsBefore=0);
    void InsColumnN(int rowN, std::vector<double>x, int whatN1=1, int QDfltsBefore=0);
    //
    void DelLineN(int rowN);
    void DelColumnN(int rowN);
    //
    void SetX(double val);//Set Component N
    void SetY(double val);
    void SetZ(double val);
    double GetX();//Get Component N
    double GetY();
    double GetZ();
    //
    std::vector<double> GetLineNAsStdVector(int rowN) const;
    std::vector<double> GetColumnNAsStdVector(int rowN) const;
    std::vector<double> GetContentAsVector1DByLines() const;
    std::vector<std::vector<double>> GetContentAsVector2DByLines() const;
    Matrix Get(int Line1N=1, int Col1N=1, int Line2N=-1, int Col2N=-1) const;
    Matrix GetVectorOfColN(int ColN) const;
    Matrix GetVectorOfLineN(int LineN) const;
    Matrix GetLineN_AsHorisontalMatrix(int LineN) const;
    //
    void ConcatAdding(const Matrix&M, int QShiftedCols=0, int FromN=1, int QAdded=0);
    void ConcatStretching(const Matrix&M, int QShiftedLines=0, int FromN=1, int QAdded=0);
    //
    void Transpose();
    Matrix GetTransposed() const;
    void LeaveMinor(int LineN, int ColumnN);
    Matrix GetMinor(int LineN, int ColN) const;
    double determinant() const;
    double AlgSuppl(int LineN, int ColumnN) const;
    void Invert();
    Matrix GetInverted() const;
    Matrix GetUnionMatrix() const;
    void ToUnion();
    //
    QString LineToString(int LineN, QString delim=", ", bool useBrackets=true);
    QString ColumnToString(int ColN, QString delim=", ", bool useBrackets=true);
    QString ToString(QString delimElem=", ", QString delimRow="; ", bool useBrackets=true);
    QString GetLineAsStr(int LineN, QString delimElem=", ", bool useBrackets=true);
    QString GetColumnAsStr(int ColN, QString delimElem=", ", bool useBrackets=true);
    QString GetAllLinesAsStr(QString delimElem=", ", QString delimLines="; ",bool useBrackets=true);
    void StdCOut1D(QString delimElem=", ", QString delimRow="; ", bool useBrackets=true);
    void StdCOut2D(QString delimElem=", ", QString delimRow="; ", bool useBrackets=true);
};

Matrix operator +(const Matrix& M1, const Matrix& M2);
Matrix operator -(const Matrix& M1, const Matrix& M2);
Matrix operator *(const Matrix& M1, const Matrix& M2);
Matrix operator *(double k, const Matrix& M2);
Matrix operator *(const Matrix& M1, double k);
Matrix operator /(const Matrix& M1, const Matrix& M2);
bool operator ==(const Matrix& M1, const Matrix& M2);
bool operator !=(const Matrix& M1, const Matrix& M2);

Matrix fApproxSolveLinAlgEqsSysSeidel(const Matrix&M, const Matrix&V, Matrix X0, TItersPrecision Precision);

Matrix GivensRotationTransformMatrix(Matrix*M, int N1, int N2, Matrix*P, TValsShowHide*vsh=NULL);
void QRDecomposition(Matrix&Q, Matrix&R,  Matrix*P=NULL, TValsShowHide* vsh = NULL);

Matrix SolveLinearSys(const Matrix&A, const Matrix&F);
Matrix SolveLinearSysBySeidel(const Matrix&A, const Matrix&F, TItersPrecision Precision, Matrix*P=NULL);

Matrix CalcMatrixOfDirCossByEulersAngles(Matrix MEulerAngles);

Matrix CalcCoordsTransformFromOldCSToNew(Matrix*CoordsInOldCSParam, Matrix*EulerAnglesParam, Matrix*CSOriginOfNewCSInOldParam=NULL, TValsShowHide*vsh=NULL);
Matrix CalcCoordsTransformFromNewCSToOld(Matrix*CoordsInNewCSParam, Matrix*EulerAnglesParam, Matrix*CSOriginOfNewCSInOldParam=NULL, TValsShowHide*vsh=NULL);

Matrix ConcatAdding(const Matrix&M1, const Matrix&M2, int QShiftedCols=0, int FromN=1, int QAdded=0);
Matrix ConcatStretching(const Matrix&M1, const Matrix&M2, int QShiftedLines=0, int FromN=1, int QAdded=0);

#endif // MYMATHLIB_H
