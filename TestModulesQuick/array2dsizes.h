#ifndef ARRAY2DSIZES_H
#define ARRAY2DSIZES_H

#include"myarraylib.h"

//class Array2DSizes
//{
//public:
//    Array2DSizes();
//};

class Array2DSizePrototype
{
public:
    Array2DSizePrototype();
    //Array2DSizePrototype(const Array2DSizePrototype&obj);//no error lif't ma I nved, qly realf ce
    Array2DSizePrototype(int QExtRows, int ExtRowsLength);
    //Array2DSizePrototype(int QExtRows, int*Ls);//, int isNotVariaAndSetByMin1Max2=1);
    //Array2DSizePrototype(std::vector<int>Ls);
    //Array2DSizePrototype operator =(const Array2DSizePrototype&obj);
    virtual ~Array2DSizePrototype() = default;
    virtual Array2DSizePrototype* clone()const;
    //virtual void Assign(const Array2DSizePrototype&obj)=0;//in tal form ak I cog ce nablb realf'd
    //virtual void Assign(Array2DSizePrototype*obj)=0;
    virtual void Construct()=0;
    virtual void SetNull()=0;
    //
    virtual bool IsRectangularType() const=0;
    virtual bool IsVariableLengthType() const=0;
    virtual bool IsRectangularFact() const=0;
    virtual bool IsVariableLengthFact() const=0;
    //
    //virtual void SetVaria();
    //virtual void SetRectangular(int L_ifVaria_min_m1_max_less=-2)=0;
    //
    virtual void SetOne()=0;
    virtual bool ExtRowNBelongsTo(int extRowN)const=0;
    virtual bool PosBelongsTo(int extRowN, int ineRowN)=0;
    virtual void Set(int QExtRows, int IneRowsLength)=0;
    virtual void Set(int QExtRows,int*Ls)=0;
    virtual void Set(std::vector<int>Ls)=0;
    //virtual void Set1ExtRowVaria(int L)=0;
    //virtual void SetVaria(int QExtRows, int IneRowsLength)=0;
    //virtual void ConvertToRectIfAllLengthesAreEqual()=0;
    //virtual void Add(int L, int ifConstAndNotEqual_Ignore0AddExisting1SetAll2=2)=0;
    //virtual void Ins(int N, int L)=0;
    virtual void AddExt(int L=0)=0;
    virtual void InsExt(int N=0, int L=0)=0;
    virtual void DelExt(int N)=0;
    virtual void DelIneRow()=0;
    virtual void AddOrInsIneRow()=0;
    //virtual void AddOrInsTo(int L=0, int N_0ToAllPreserveType_GTT0AndMakeVaria=0)=0;
    //virtual void DelFrom(int L, int N_0ToAllPreserveType_GTT0AndMakeVaria=0)=0;
    virtual void AddOrInsTo(int N, int L)=0;
    virtual void DelFrom(int N, int L)=0;
    virtual int GetQExtRows() const=0;
    //virtual bool GetIfIsVariaLength() const=0;
    //virtual bool isVariaLength()=0;
    //virtual bool isRectangular() const=0;
    virtual int GetLength(int N=0) const=0;
    virtual int GetMinLength()const=0;
    virtual int GetMaxLength()const=0;
    virtual void SwapExtRows(int N1, int N2)=0;
    virtual void Reverse()=0;
    //virtual void SetQExtRows(int N)=0;
    //virtual void SetIneRowsLength(int L_Minus1MinMinus2Max=-2)=0;
    //virtual void SetLength(int L, int N_0ToAllPreserveType_GTT0AndMakeVaria_Minus1ToAllMakeNotVaria=0)=0;
    virtual void SetLength(int L, int N=0)=0;
    virtual QString GetAsQString(QString delim=", ")const=0;
    virtual std::string GetAsStdString(std::string delim=", ") const=0;
    virtual void ShowToConsole() const=0;
    //virtual bool isSame(const Array2DSizePrototype&obj)const=0;
    //virtual bool isSame(Array2DSizePrototype*obj)const=0;
    virtual int GetIneRowsLength()const=0;
    virtual void Transpose()=0;
    virtual bool IsForTranspose() const=0;
    virtual int LengthFullOfIneRowN(int N) const=0;//from start to 1st stop
    virtual std::vector<int>GetLengthes();
    virtual int GetQElements();
};

class Array2DSizeVarLenP: public Array2DSizePrototype
{
    int*Ls;
    int QExtRows;
public:
    //virtual Array2DSize();
    Array2DSizeVarLenP(const Array2DSizeVarLenP&obj);//I men no not ob in dyn polymorf S nablb realf'd
    Array2DSizeVarLenP(int QExtRows, int IneRowsLength);
    Array2DSizeVarLenP(int QExtRows, int*Ls);//, int isNotVariaAndSetByMin1Max2=1);
    Array2DSizeVarLenP(std::vector<int>Ls);
    //Array2DSizeVarLenP operator =(const Array2DSizePrototype&obj);
    virtual Array2DSizeVarLenP* clone()const ;
    virtual ~Array2DSizeVarLenP();
    //virtual void Assign(const Array2DSizePrototype&obj);//I men no not ob in dyn polymorf S nablb realf'd
    //virtual void Assign(Array2DSizeVarLenP*obj);//I hope ce abls ob Array2DSizeVarLenP s'heir of Array2DSizePrototype
    virtual void Construct();
    virtual void SetNull();
    //
    virtual bool IsRectangularType() const;
    virtual bool IsVariableLengthType() const;
    virtual bool IsRectangularFact() const;
    virtual bool IsVariableLengthFact() const;
    //
    //virtual void SetVaria();
    //virtual void SetRectangular(int L_ifVaria_min_m1_max_less=-2);
    //
    virtual void SetOne();
    virtual bool ExtRowNBelongsTo(int extRowN)const;
    virtual bool PosBelongsTo(int extRowN, int ineRowN);
    virtual void Set(int QExtRows, int IneRowsLength);
    virtual void Set(int QExtRows,int*Ls);
    virtual void Set(std::vector<int>Ls);
    //virtual void Set1ExtRowVaria(int L);
    //virtual void SetVaria(int QExtRows, int IneRowsLength);
    //virtual void ConvertToRectIfAllLengthesAreEqual();
    //virtual void Add(int L, int ifConstAndNotEqual_Ignore0AddExisting1SetAll2=2);
    //virtual void Ins(int N, int L);
    virtual void AddExt(int L=0);
    virtual void InsExt(int N=0, int L=0);
    virtual void DelExt(int N);
    virtual void DelIneRow();
    virtual  void AddOrInsIneRow();
    //virtual void AddOrInsTo(int L=0, int N_0ToAllPreserveType_GTT0AndMakeVaria=0);
    //virtual void DelFrom(int L, int N_0ToAllPreserveType_GTT0AndMakeVaria=0);
    virtual void AddOrInsTo(int N, int L);
    virtual void DelFrom(int N, int L);
    virtual int GetQExtRows() const;
    //virtual bool GetIfIsVariaLength() const;
    //virtual bool isVariaLength();
    //virtual bool isRectangular() const;
    virtual int GetLength(int N=0) const;
    virtual int GetMinLength()const;
    virtual int GetMaxLength()const;
    virtual void SwapExtRows(int N1, int N2);
    virtual void Reverse();
    //virtual void SetQExtRows(int N);
    //virtual void SetIneRowsLength(int L_Minus1MinMinus2Max=-1, bool preserveTypeIfVaria=false);
    //virtual void SetLength(int L, int N_0ToAllPreserveType_GTT0AndMakeVaria_Minus1ToAllMakeNotVaria=0);
    virtual void SetLength(int L, int N=0);
    virtual QString GetAsQString(QString delim=", ")const;
    virtual std::string GetAsStdString(std::string delim=", ") const;
    virtual void ShowToConsole() const;
    //virtual bool isSame(const Array2DSizeVarLenP&obj)const;
    //virtual bool isSame(Array2DSizeVarLenP*obj)const;
    virtual int GetIneRowsLength()const;
    virtual void Transpose();
    virtual bool IsForTranspose() const;
    virtual int LengthFullOfIneRowN(int N) const;//from start to 1st stop
    virtual std::vector<int>GetLengthes();
    virtual int GetQElements();
};

class Array2DSizeVarLenV: public Array2DSizePrototype
{
    std::vector<int>Ls;
public:
    //virtual Array2DSize();
    Array2DSizeVarLenV(const Array2DSizeVarLenV&obj);
    Array2DSizeVarLenV(int QExtRows, int IneRowsLengthe);
    Array2DSizeVarLenV(int QExtRows, int*Ls);//, int isNotVariaAndSetByMin1Max2=1);
    Array2DSizeVarLenV(std::vector<int>Ls);
    //Array2DSizeVarLenV operator =(const Array2DSizeVarLenV&obj);
    virtual ~Array2DSizeVarLenV();
    virtual Array2DSizeVarLenV* clone()const;
    //virtual void Assign(const Array2DSizeVarLenV&obj);
    //virtual void Assign(Array2DSizeVarLenV*obj);
    virtual void Construct();
    virtual void SetNull();
    //
    virtual bool IsRectangularType() const;
    virtual bool IsVariableLengthType() const;
    virtual bool IsRectangularFact() const;
    virtual bool IsVariableLengthFact() const;
    //
    //virtual void SetVaria();
    //virtual void SetRectangular(int L_ifVaria_min_m1_max_less=-2);
    //
    virtual void SetOne();
    virtual bool ExtRowNBelongsTo(int extRowN)const;
    virtual bool PosBelongsTo(int extRowN, int ineRowN);
    virtual void Set(int QExtRows, int ExtRowsLength);
    virtual void Set(int QExtRows,int*Ls);
    virtual void Set(std::vector<int>Ls);
    virtual void Set1ExtRowVaria(int L);
    virtual void SetVaria(int QExtRows, int IneRowsLength);
    virtual void ConvertToRectIfAllLengthesAreEqual();
    //virtual void Add(int L, int ifConstAndNotEqual_Ignore0AddExisting1SetAll2=2);
    //virtual void Ins(int N, int L);
    virtual void AddExt(int L=0);
    virtual void InsExt(int N=0, int L=0);
    virtual void DelExt(int N);
    virtual void DelIneRow();
    virtual  void AddOrInsIneRow();
    //virtual void AddOrInsTo(int L=0, int N_0ToAllPreserveType_GTT0AndMakeVaria=0);
    //virtual void DelFrom(int L, int N_0ToAllPreserveType_GTT0AndMakeVaria=0);
    virtual void AddOrInsTo(int N, int L);
    virtual void DelFrom(int N, int L);
    virtual int GetQExtRows() const;
    virtual bool GetIfIsVariaLength() const;
    virtual bool isVariaLength();
    virtual bool isRectangular() const;
    virtual int GetLength(int N=0) const;
    virtual int GetMinLength()const;
    virtual int GetMaxLength()const;
    virtual void SwapExtRows(int N1, int N2);
    virtual void Reverse();
    //virtual void SetQExtRows(int N);
    //virtual void SetIneRowsLength(int L_Minus1MinMinus2Max=-1, bool preserveTypeIfVaria=false);
    //virtual void SetLength(int L, int N_0ToAllPreserveType_GTT0AndMakeVaria_Minus1ToAllMakeNotVaria=0);
    virtual void SetLength(int L, int N=0);
    virtual QString GetAsQString(QString delim=", ")const;
    virtual std::string GetAsStdString(std::string delim=", ") const;
    virtual void ShowToConsole() const;
    //virtual bool isSame(const Array2DSizeVarLenV&obj)const;
    virtual bool isSame(Array2DSizeVarLenV*obj)const;
    virtual int GetIneRowsLength()const;
    virtual void Transpose();
    virtual bool IsForTranspose() const;
    virtual int LengthFullOfIneRowN(int N) const;//from start to 1st stop
    virtual std::vector<int>GetLengthes();
    virtual int GetQElements();
};


class Array2DSizeRect: public Array2DSizePrototype
{
    int QExtRows;
    int ExtRowLength;
public:
    //virtual Array2DSize();

    Array2DSizeRect(int QExtRows, int IneRowsLength);
    Array2DSizeRect(int QExtRows, int*Ls);//, int isNotVariaAndSetByMin1Max2=1);
    Array2DSizeRect(std::vector<int>Ls);
    Array2DSizeRect(const Array2DSizeRect&obj);
   // Array2DSizeRect operator =(const Array2DSizeRect&obj);
    virtual ~Array2DSizeRect();
   // virtual void Assign(const Array2DSizeRect&obj);
    virtual void Assign(Array2DSizeRect*obj);
    virtual void Construct();
    virtual void SetNull();
    //
    virtual bool IsRectangularType() const;
    virtual bool IsVariableLengthType() const;
    virtual bool IsRectangularFact() const;
    virtual bool IsVariableLengthFact() const;
    //
    //virtual void SetVaria();
    //virtual void SetRectangular(int L_ifVaria_min_m1_max_less=-2);
    //
    virtual void SetOne();
    virtual bool ExtRowNBelongsTo(int extRowN)const;
    virtual bool PosBelongsTo(int extRowN, int ineRowN);
    virtual void Set(int QExtRows, int IneRowsLength);
    virtual void Set(int QExtRows,int*Ls);
    virtual void Set(std::vector<int>Ls);
    virtual void Set1ExtRowVaria(int L);
    //virtual void SetVaria(int QExtRows, int IneRowsLength);
    //virtual void ConvertToRectIfAllLengthesAreEqual();
    //virtual void Add(int L, int ifConstAndNotEqual_Ignore0AddExisting1SetAll2=2);
    //virtual void Ins(int N, int L);
    virtual void AddExt(int L=0);
    virtual void InsExt(int N=0, int L=0);
    virtual void DelExt(int N);
    virtual void DelIneRow();
    virtual  void AddOrInsIneRow();
    //virtual void AddOrInsTo(int L=0, int N_0ToAllPreserveType_GTT0AndMakeVaria=0);
    //virtual void DelFrom(int L, int N_0ToAllPreserveType_GTT0AndMakeVaria=0);
    virtual void AddOrInsTo(int N, int L);
    virtual void DelFrom(int N, int L);
    virtual int GetQExtRows() const;
    //virtual bool GetIfIsVariaLength() const;
    //virtual bool isVariaLength();
    //virtual bool isRectangular() const;
    virtual int GetLength(int N=0) const;
    virtual int GetMinLength()const;
    virtual int GetMaxLength()const;
    virtual void SwapExtRows(int N1, int N2);
    virtual void Reverse();
    //virtual void SetQExtRows(int N);
    //virtual void SetIneRowsLength(int L_Minus1MinMinus2Max=-1, bool preserveTypeIfVaria=false);
    //virtual void SetLength(int L, int N_0ToAllPreserveType_GTT0AndMakeVaria_Minus1ToAllMakeNotVaria=0);
    virtual void SetLength(int L, int N=0);
    virtual QString GetAsQString(QString delim=", ")const;
    virtual std::string GetAsStdString(std::string delim=", ") const;
    virtual void ShowToConsole() const;
    //virtual bool isSame(const Array2DSizeRect&obj)const;
    virtual bool isSame(Array2DSizeRect*obj)const;
    virtual int GetIneRowsLength()const;
    virtual void Transpose();
    virtual bool IsForTranspose() const;
    virtual int LengthFullOfIneRowN(int N) const;//from start to 1st stop
    virtual std::vector<int>GetLengthes();
    virtual int GetQElements();
};

class Array2DSizeRectShortened: public Array2DSizePrototype
{
    int QExtRows;
    int ExtRowLength;
public:
    //virtual Array2DSize();

    Array2DSizeRectShortened(int QExtRows, int IneRowsLength);
    Array2DSizeRectShortened(int QExtRows, int*Ls);//, int isNotVariaAndSetByMin1Max2=1);
    Array2DSizeRectShortened(std::vector<int>Ls);
    Array2DSizeRectShortened(const Array2DSizeRectShortened&obj);
   // Array2DSizeRect operator =(const Array2DSizeRectShortened&obj);
    virtual ~Array2DSizeRectShortened();
   // virtual void Assign(const Array2DSizeRectShortened&obj);
    virtual void Assign(Array2DSizeRectShortened*obj);
    virtual void Construct();
    virtual void SetNull();
    //
    virtual bool IsRectangularType() const;
    virtual bool IsVariableLengthType() const;
    virtual bool IsRectangularFact() const;
    virtual bool IsVariableLengthFact() const;
    //
    //virtual void SetVaria();
    //virtual void SetRectangular(int L_ifVaria_min_m1_max_less=-2);
    //
    virtual void SetOne();
    virtual bool ExtRowNBelongsTo(int extRowN)const;
    virtual bool PosBelongsTo(int extRowN, int ineRowN);
    virtual void Set(int QExtRows, int IneRowsLength);
    virtual void Set(int QExtRows,int*Ls);
    virtual void Set(std::vector<int>Ls);
    virtual void Set1ExtRowVaria(int L);
    //virtual void SetVaria(int QExtRows, int IneRowsLength);
    //virtual void ConvertToRectIfAllLengthesAreEqual();
    //virtual void Add(int L, int ifConstAndNotEqual_Ignore0AddExisting1SetAll2=2);
    //virtual void Ins(int N, int L);
    virtual void AddExt(int L=0);
    virtual void InsExt(int N=0, int L=0);
    virtual void DelExt(int N);
    virtual void DelIneRow();
    virtual  void AddOrInsIneRow();
    //virtual void AddOrInsTo(int L=0, int N_0ToAllPreserveType_GTT0AndMakeVaria=0);
    //virtual void DelFrom(int L, int N_0ToAllPreserveType_GTT0AndMakeVaria=0);
    virtual void AddOrInsTo(int N, int L);
    virtual void DelFrom(int N, int L);
    virtual int GetQExtRows() const;
    //virtual bool GetIfIsVariaLength() const;
    //virtual bool isVariaLength();
    //virtual bool isRectangular() const;
    virtual int GetLength(int N=0) const;
    virtual int GetMinLength()const;
    virtual int GetMaxLength()const;
    virtual void SwapExtRows(int N1, int N2);
    virtual void Reverse();
    //virtual void SetQExtRows(int N);
    //virtual void SetIneRowsLength(int L_Minus1MinMinus2Max=-1, bool preserveTypeIfVaria=false);
    //virtual void SetLength(int L, int N_0ToAllPreserveType_GTT0AndMakeVaria_Minus1ToAllMakeNotVaria=0);
    virtual void SetLength(int L, int N=0);
    virtual QString GetAsQString(QString delim=", ")const;
    virtual std::string GetAsStdString(std::string delim=", ") const;
    virtual void ShowToConsole() const;
    //virtual bool isSame(const Array2DSizeRectShortened&obj)const;
    virtual bool isSame(Array2DSizeRectShortened*obj)const;
    virtual int GetIneRowsLength()const;
    virtual void Transpose();
    virtual bool IsForTranspose() const;
    virtual int LengthFullOfIneRowN(int N) const;//from start to 1st stop
    virtual std::vector<int>GetLengthes();
    virtual int GetQElements();
};

class Array2DSizeDPS //Dynamic Polymorph Shell
{
    Array2DSizePrototype*size;
public:

    Array2DSizeDPS(int QExtRows, int IneRowsLength, bool isNotVaria=true);
    Array2DSizeDPS(int QExtRows, int*Ls, bool RectIfAllSame=true);
    Array2DSizeDPS(std::vector<int>Ls, bool RectIfAllSame=true);
    Array2DSizeDPS(const Array2DSizeDPS&obj);
    ~Array2DSizeDPS();
    void Assign(const Array2DSizeDPS&obj);
    Array2DSizeDPS& operator =(const Array2DSizeDPS&obj);
    //void Construct();
    void SetNull();
    //
    bool IsRectangularType() const;
    bool IsVariableLengthType() const;
    bool IsRectangularFact() const;
    bool IsVariableLengthFact() const;
    //
    //virtual void SetVaria();
    //virtual void SetRectangular(int L_ifVaria_min_m1_max_less=-2);
    //
    void SetOne(bool isVaria=false);
    bool ExtRowNBelongsTo(int extRowN)const;
    bool PosBelongsTo(int extRowN, int ineRowN);
    void Set(int QExtRows, int IneRowsLength, bool PreserveType=true);
    void Set(int QExtRows,int*Ls, bool RectIfAllSame=true);

    //void Add(int L, int ifConstAndNotEqual_Ignore0AddExisting1SetAll2=2);
    //void Ins(int N, int L);
    void AddExt(int L=0);
    void InsExt(int N=0, int L=0);
    void DelExt(int N);
    void DelIneRow();
    void AddOrInsIneRow();
    //virtual void AddOrInsTo(int L=0, int N_0ToAllPreserveType_GTT0AndMakeVaria=0);
    //virtual void DelFrom(int L, int N_0ToAllPreserveType_GTT0AndMakeVaria=0);
    void AddOrInsTo(int N, int L);
    void DelFrom(int N, int L);
    int GetQExtRows() const;
    //bool GetIfIsVariaLength() const;
    //bool isVariaLength();
    //bool isRectangular() const;
    int GetLength(int N=0) const;
    int GetMinLength()const;
    int GetMaxLength()const;
    void SwapExtRows(int N1, int N2);
    void Reverse();
    //void SetQExtRows(int N);
    //void SetIneRowsLength(int L_MinMinus1MaxLess=-2, bool preserveTypeIfVaria=false);
    //void SetLength(int L, int N_0ToAllPreserveType_GTT0AndMakeVaria_Minus1ToAllMakeNotVaria=0);
    void SetLength(int L, int N=0);
    QString GetAsQString(QString delim=", ")const;
    std::string GetAsStdString(std::string delim=", ") const;
    void ShowToConsole() const;
    bool isSame(const Array2DSizeDPS&obj)const;
    int GetIneRowsLength(bool EvenIfVarButAllEqual=false)const;
    void Transpose();
    bool IsForTranspose() const;
    int LengthFullOfIneRowN(int N) const;//from start to 1st stop
    //
    void Set1ExtRowVaria(int L);
    void SetVectorPotentiallyVaria(int Q);
    void SetVaria(int QExtRows, int IneRowsLength);
    void ConvertToRectIfAllLengthesAreEqual();
    void ConvertToRect(int L_ifVaria_min_m1_max_less=-1);
    void ConvertToVaria();

    virtual std::vector<int>GetLengthes();
    virtual int GetQElements();
};



#endif // ARRAY2DSIZES_H
