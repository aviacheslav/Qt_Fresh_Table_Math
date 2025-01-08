#include "myarraylib.h"

//MyArrayLib::MyArrayLib()
//{
//
//}

// =======================================================================================================

//NumerSub Row - No Need, es
//template <typename T>std::vector<T> Array1DGetSubArray_byLs_PreMark_v1(T*x, int Lold, int Lreq=0, int whatFromN=1, int FirstDefaultValues=0){
//template <typename T>std::vector<T> Array1DGetSubArray_byLs_PreMark_v1(T*x, int Lold, int Lreq=0, int whatFromN=1, int FirstDefaultValues=0, T DefaultValue){
//template <typename T>std::vector<T> Array1DGetSubArray_byLs_PreMark_v1(T*x, int Lold, int Lreq=0, int whatFromN=1, int FirstDefaultValues=0, T*DefaultRow_LongEnough){
/*
std::vector<double>NumbersSubRow(std::vector<double> x, int LreqGiven, int FromN, int QDefaultBefore, double DefaultValParam){
    std::vector<double> R;
    //double*DefaultVal=new double;
    //*DefaultVal=DefaultValParam;
    //template <typename T>std::vector<T> Array1DGetSubArray_byLs_PreMark(T*x, int Lold, int Lreq=0, int whatFromN=1, int FirstDefaultValues=0, T*DfltValParam=NULL){
    //template <typename T>std::vector<T> Array1DGetSubArray_byLs_PreMark(std::vector<T>x, int Lreq=0, int whatFromN=1, int FirstDefaultValues=0, T*DfltValParam=NULL){
    //R= Array1DGetSubArray_byLs_PreMark(x, Lreq, FromN, QDefaultBefore, DefaultVal);
    //delete DefaultVal;
    //void Calc_Array1DSubArray_byLs_MarkupNs_vCycls(int&LreqCalcd, int&preN1_, int&preN2_, int&ownN1_, int&ownN2, int&ForOwnN1_, int&ForOwnN2, int&postN1, int&postN2_, int Lold, int LreqGiven, int whatFromN, int QDefaultValuesBefore){
    int LreqCalcd=0, preN1_=0, preN2_=0, ownN1_=0, ownN2=0, ForOwnN1_=0, ForOwnN2=0, postN1=0, postN2_=0, Lold=x.size();
    double CurVal;
    Calc_Array1DSubArray_byLs_MarkupNs_vFmls(LreqCalcd, preN1_, preN2_, ownN1_, ownN2, ForOwnN1_, ForOwnN2, postN1, postN2_, Lold, LreqGiven, FromN, QDefaultBefore);
    if(preN1_>0){
        for(int i=preN1_; i<=preN2_; i++){
            R.push_back(DefaultValParam);
        }
    }
    if(ownN1_>0){
        for(int i=ownN1_; i<=ownN2; i++){
            CurVal=x.at(i-1);
            R.push_back(CurVal);
        }
    }
    if(postN1>0){
        for(int i=postN1; i<=postN2_; i++){
            R.push_back(DefaultValParam);
        }
    }
    return R;
}

std::vector<double>NumbersSubRow(std::vector<double> x, int LreqGiven, int FromN, int QDefaultBefore, double*rowDflt){
    std::vector<double> R;
    int LreqCalcd=0, preN1_=0, preN2_=0, ownN1_=0, ownN2=0, ForOwnN1_=0, ForOwnN2=0, postN1=0, postN2_=0, Lold=x.size();
    double CurVal;
    Calc_Array1DSubArray_byLs_MarkupNs_vFmls(LreqCalcd, preN1_, preN2_, ownN1_, ownN2, ForOwnN1_, ForOwnN2, postN1, postN2_, Lold, LreqGiven, FromN, QDefaultBefore);
    if(preN1_>0){
        for(int i=preN1_; i<=preN2_; i++){
            CurVal=rowDflt[i-1];
            R.push_back(CurVal);
        }
    }
    if(ownN1_>0){
        for(int i=ownN1_; i<=ownN2; i++){
            CurVal=x.at(i-1);
            R.push_back(CurVal);
        }
    }
    if(postN1>0){
        for(int i=postN1; i<=postN2_; i++){
            CurVal=rowDflt[i-1];
            R.push_back(CurVal);
        }
    }
    return R;
}


std::vector<int>NumbersSubRow(std::vector<int> x, int LreqGiven, int FromN, int QDefaultBefore, int DefaultValParam){
    std::vector<int> R;
    int CurVal;
    int LreqCalcd=0, preN1_=0, preN2_=0, ownN1_=0, ownN2=0, ForOwnN1_=0, ForOwnN2=0, postN1=0, postN2_=0, Lold=x.size();
    //int, LreqGiven=0, whatFromN=0, QDefaultValuesBefore=0;
    //int*DefaultVal=new int;
    //*DefaultVal=DefaultValParam;
    //template <typename T>std::vector<T> Array1DGetSubArray_byLs_PreMark(T*x, int Lold, int Lreq=0, int whatFromN=1, int FirstDefaultValues=0, T*DfltValParam=NULL){
    //template <typename T>std::vector<T> Array1DGetSubArray_byLs_PreMark(std::vector<T>x, int Lreq=0, int whatFromN=1, int FirstDefaultValues=0, T*DfltValParam=NULL){
    //R= Array1DGetSubArray_byLs_PreMark(x,  Lreq, FromN, QDefaultBefore, DefaultVal);
    //delete DefaultVal;
    //void Calc_Array1DSubArray_byLs_MarkupNs_vCycls(int&LreqCalcd, int&preN1_, int&preN2_, int&ownN1_, int&ownN2, int&ForOwnN1_, int&ForOwnN2, int&postN1, int&postN2_, int Lold, int LreqGiven, int whatFromN, int QDefaultValuesBefore){
    Calc_Array1DSubArray_byLs_MarkupNs_vFmls(LreqCalcd, preN1_, preN2_, ownN1_, ownN2, ForOwnN1_, ForOwnN2, postN1, postN2_, Lold, LreqGiven, FromN, QDefaultBefore);
    if(preN1_>0){
        for(int i=preN1_; i<=preN2_; i++){
            R.push_back(DefaultValParam);
        }
    }
    if(ownN1_>0){
        for(int i=ownN1_; i<=ownN2; i++){
            CurVal=x.at(i-1);
            R.push_back(CurVal);
        }
    }
    if(postN1>0){
        for(int i=postN1; i<=postN2_; i++){
            R.push_back(DefaultValParam);
        }
    }
    return R;
}



std::vector<int>NumbersSubRow(std::vector<int> x, int LreqGiven, int FromN, int QDefaultBefore, int*rowDflt){
    std::vector<int> R;
    int LreqCalcd=0, preN1_=0, preN2_=0, ownN1_=0, ownN2=0, ForOwnN1_=0, ForOwnN2=0, postN1=0, postN2_=0, Lold=x.size();
    int CurVal;
    Calc_Array1DSubArray_byLs_MarkupNs_vFmls(LreqCalcd, preN1_, preN2_, ownN1_, ownN2, ForOwnN1_, ForOwnN2, postN1, postN2_, Lold, LreqGiven, FromN, QDefaultBefore);
    if(preN1_>0){
        for(int i=preN1_; i<=preN2_; i++){
            CurVal=rowDflt[i-1];
            R.push_back(CurVal);
        }
    }
    if(ownN1_>0){
        for(int i=ownN1_; i<=ownN2; i++){
            CurVal=x.at(i-1);
            R.push_back(CurVal);
        }
    }
    if(postN1>0){
        for(int i=postN1; i<=postN2_; i++){
            CurVal=rowDflt[i-1];
            R.push_back(CurVal);
        }
    }
    return R;
}
*/



Array2DSize Array2DSize_ConcatAdding(const Array2DSize& size1, const Array2DSize& size2, int&ERN1, int&IRN1, int&ERN2, int&IRN2, int shift, int Q, int FromN, int QDfltValsBefore, bool SetRect){
    int Q1=size1.GerQExtRows(), Q2=size2.GerQExtRows(), Q3, lmax1=size1.GetMaxLength(), lmax2=size2.GetMaxLength(), lmax3;
    Array2DSize size3;
    int LreqCalcd=0, preN1_=0, preN2_=0, ownN1_=0, ownN2=0, ForOwnN1_=0, ForOwnN2=0, postN1=0, postN2_=0, Lold=Q2;
    if(Q==0){
        Q=Q2-FromN+1+QDfltValsBefore;
    }
    Q3=Q1+Q;
    if(FromN<0){
        FromN+=(Q2+1);
    }
    Calc_Array1DSubArray_byLs_MarkupNs_vFmls(LreqCalcd, preN1_, preN2_, ownN1_, ownN2, ForOwnN1_, ForOwnN2, postN1, postN2_, Lold, Q, FromN, QDefaultBefore);
    if(shift>=0){
        //Q3 = Q1>shift+Q ? Q1 : shift+Q;
        //Calc_Array1DSubArray_byLs_MarkupNs_vFmls(LreqCalcd, preN1_, preN2_, ownN1_, ownN2, ForOwnN1_, ForOwnN2, postN1, postN2_, Lold, LreqGiven, FromN, QDefaultBefore);
        //Calc_Array1DSubArray_byLs_MarkupNs_vFmls(LreqCalcd, preN1_, preN2_, ownN1_, ownN2, ForOwnN1_, ForOwnN2, postN1, postN2_, Lold, Q, FromN, QDefaultBefore);
        ERN1=1;
        IRN1=1;
        if(ForOwnN2_>0){
            ERN2=Q1+ForOwnN1_;
            IRN2=shift+1;
        }else{
             ERN2=0;
             IRN2=0;
        }
        lmax3 = lmax1 >= shift + lmax2 ? lmax1 : shift + lmax2;
        size3.SetNull();
        if(SetResict){
            size3.Set(Q3, lmax3);
        }else{
            for(int i=1; i<=Q1; i++){
                size3.AddExt(size1.GetLength(i));
            }
            if(preN2_>0){
                for(int i=ownN1_; i<=ownN2; i++){
                    size3.AddExt(lmax3);
                }
            }
            if(ERN2>0){
               for(int i=ownN1_; i<=ownN2; i++){
                   size3.AddExt(shift+size2.GetLength(i));
                }
            }
            if(postN2_>0){
                for(int i=postN1; i<=postN2_; i++){
                    size3.AddExt(lmax3);
                }
            }
        }
    }else{
        //Q3 = -shift+Q1 > Q ? -shift+Q1 : Q;
        //Calc_Array1DSubArray_byLs_MarkupNs_vFmls(LreqCalcd, preN1_, preN2_, ownN1_, ownN2, ForOwnN1_, ForOwnN2, postN1, postN2_, Lold, Q, FromN, QDefaultBefore);
        ERN1=1;
        IRN1=-shift+1;
        if(ForOwnN2_>0){
            ERN2=Q1+ForOwnN1_;
            IRN2=1;
        }else{
             ERN2=0;
             IRN2=0;
        }
        lmax3 = -shift + lmax1 >= lmax2 ? -shift + lmax1 : lmax2;
        size3.SetNull();
        if(SetResict){
            size3.Set(Q3, lmax3);
        }else{
            for(int i=1; i<=Q1; i++){
                size3.AddExt(-shift+size1.GetLength(i));
            }
            if(preN2_>0){
                for(int i=ownN1_; i<=ownN2; i++){
                    size3.AddExt(lmax3);
                }
            }
            if(ERN2>0){
               for(int i=ownN1_; i<=ownN2; i++){
                   size3.AddExt(size2.GetLength(i));
                }
            }
            if(postN2_>0){
                for(int i=postN1; i<=postN2_; i++){
                    size3.AddExt(lmax3);
                }
            }
        }
    }
    return size3;
}


Array2DSize Array2DSize_ConcatStretching(const Array2DSize& size1, const Array2DSize& size2, int&ERN1, int&EIN1, int&ERN2, int&EIN2, int shift, int Q, int FromN, int QDfltValsBefore, bool Stretch1ToMax_Not_AddToLast, bool AppendToRect){
    int Q1=size1.GerQExtRows(), Q2=size2.GerQExtRows(), Q3, lmax1=size1.GetMaxLength(), lmax2=size2.GetMaxLength(), lmax3;
    Array2DSize size3;
    int LreqCalcd=0, preN1_=0, preN2_=0, ownN1_=0, ownN2=0, ForOwnN1_=0, ForOwnN2=0, postN1=0, postN2_=0, Lold=Q2;
    size3.SetNull();
    //if(Stretch1ToMax_Not_AddToLast){
    //
    //}else{
    //
    //}
    Calc_Array1DSubArray_byLs_MarkupNs_vFmls(LreqCalcd, preN1_, preN2_, ownN1_, ownN2, ForOwnN1_, ForOwnN2, postN1, postN2_, Lold, Q, FromN, QDefaultBefore);
    if(shift>=0){
        Q3 = Q1>shift+Q ? Q1 : shift+Q;
        if(Stretch1ToMax_Not_AddToLast){
            if(AppendToRect){
                lmax3 = lmax1 + lmax2;
                size3.Set(Q3, lmax3);
            }else{
                if(shift>Q1){
                    for(int i=1; i<=Q1; i++){
                        size3.AddExt(size1.SetLength(i););
                    }
                    for(int i=Q1+1; i<=shift; i++){
                        size3.AddExt(lmax1 + lmax2);
                    }
                    if(Stretch1ToMax_Not_AddToLast){
                        if(preN2_>0){
                            for(int i=preN1_; i<=preN2_; i++){
                                size3.AddExt(lmax1 + lmax2);
                            }
                        }
                        if(ownN2>0){
                            for(int i=ownN1_; i<=ownN2; i++){
                                size3.AddExt(size2.GetLength(i));
                            }
                        }
                        if(postN2_>0){
                            for(int i=postN1; i<=postN2_; i++){
                                size3.AddExt(lmax1 + lmax2);
                            }
                        }
                    }else{
                        for(int i=1; i<=shift; i++){

                        }
                    }
                }else{

                }
            }
        }else{
            if(AppendToRect){

            }else{

            }
        }
    }else{
        Q3 = -shift+Q1 > Q ? -shift+Q1 : Q;
        if(Stretch1ToMax_Not_AddToLast){
            if(AppendToRect){

            }else{

            }
        }else{
            if(AppendToRect){

            }else{

            }
        }
    }
    return size3;
}





//===================================================================================================

int Array2DSizeBasedOn1D_PosByCoords(int ExtRowN, int IneRowN, int QExtRows, int*Ls){
    int PosN=0, curL;
    if(ExtRowN>=1 && ExtRowN<=QExtRows){
        if(Ls!=NULL){
            for(int i=1; i<=ExtRowN-1; i++){
                curL=Ls[i-1];
                PosN+=curL;
            }
            PosN+=IneRowN;
        }else{
            for(int i=1; i<=ExtRowN-1; i++){
                curL=0;
                PosN+=curL;
            }
            PosN+=IneRowN;
        }
    }
    return PosN;
}
int Array2DSizeBasedOn1D_PosByCoords(int ExtRowN, int IneRowN, std::vector<int>Ls){
    return Array2DSizeBasedOn1D_PosByCoords(ExtRowN, IneRowN, Ls.size(), Ls.data());
}
void Array2DSizeBasedOn1D_CoordsByPos(int PosN, int QExtRows, int*Ls, int&ExtRowN, int&IneRowN){
    int Ncur=0, Lcur, QElements=0;
    ExtRowN=0;
    IneRowN=0;
    bool contin, reasonably;
    reasonably=(PosN>=1 && Ls!=NULL && PosN<=QElements);
    if(reasonably){
        for(int i=1; i<=QExtRows; i++){
            Lcur=Ls[i-1];
            QElements+=Lcur;
        }
    }
    contin=reasonably;
    while(contin){
        ExtRowN++;
        Lcur=Ls[ExtRowN-1];
        if(Ncur<PosN && Ncur+Lcur>=PosN){
            contin=false;
        }
    }
    if(reasonably){
        IneRowN=PosN-Ncur;
    }
}
void Array2DSizeBasedOn1D_CoordsByPos(int PosN, std::vector<int>Ls, int&ExtRowN, int&IneRowN){
    Array2DSizeBasedOn1D_CoordsByPos(PosN, Ls.size(), Ls.data(), ExtRowN, IneRowN);
}

int Array2DSizeBasedOn1D_ExtRowNByPos(int PosN, std::vector<int>Ls){
    int IneRowN=0, ExtRowN=0;
    Array2DSizeBasedOn1D_CoordsByPos(PosN, Ls, ExtRowN, IneRowN);
    return ExtRowN;
}
int Array2DSizeBasedOn1D_IneRowNByPos(int PosN, std::vector<int>Ls){
    int IneRowN=0, ExtRowN=0;
    Array2DSizeBasedOn1D_CoordsByPos(PosN, Ls, ExtRowN, IneRowN);
    return ExtRowN;
}
int Array2DSizeBasedOn1D_ExtRowNByPos(int PosN, int QExtRows, int*Ls){
    int IneRowN=0, ExtRowN=0;
    Array2DSizeBasedOn1D_CoordsByPos(PosN, QExtRows, Ls, ExtRowN, IneRowN);
    return ExtRowN;
}
int Array2DSizeBasedOn1D_IneRowNByPos(int PosN, int QExtRows, int*Ls){
    int IneRowN=0, ExtRowN=0;
    Array2DSizeBasedOn1D_CoordsByPos(PosN, QExtRows, Ls, ExtRowN, IneRowN);
    return IneRowN;
}


//===================================================================================================

Int2DVectorOfSizeOrCoords::Int2DVectorOfSizeOrCoords(int d1, int d2){
    this->d1=d1; this->d2=d2;
}

//---------------------------------------------------------------------------------------------------

Array2DSize::Array2DSize(){
    this->Construct();
}
Array2DSize::Array2DSize(const Array2DSize&obj){
    this->Construct();
    this->Assign(obj);
}
Array2DSize::Array2DSize(int QExtRows, int  IneRowsLength, bool isNotVaria){
    this->Construct();
    this->QExtRows=QExtRows;
    if(isNotVaria){
        this->IneRowsLength=IneRowsLength;
        this->Ls=NULL;
    }else{
        this->Ls=new int[this->QExtRows];
        for(int i=1; i<=this->QExtRows; i++){
             this->Ls[i-1]=IneRowsLength;
        }
        this->IneRowsLength=0;
    }
}
Array2DSize::Array2DSize(int QExtRows, int*Ls){//, int isNotVariaAndSetByMin1Max2){
    this->Construct();
    this->QExtRows=QExtRows;
    this->IneRowsLength=0;
    this->Ls=new int[this->QExtRows];
    if(Ls!=NULL){
        for(int i=1; i<=this->QExtRows; i++){
            this->Ls[i-1]=Ls[i-1];
        }
    }
}
Array2DSize Array2DSize::operator =(const Array2DSize&obj){
    this->Assign(obj);
    return *this;
}
Array2DSize::~Array2DSize(){
    this->SetNull();
}
void Array2DSize::Assign(const Array2DSize&obj){
    this->SetNull();
    this->QExtRows=obj.QExtRows;
    if(obj.Ls!=NULL){
        this->IneRowsLength=0;
        if(this->Ls!=NULL) delete [] this->Ls;
        this->Ls=new int[obj.QExtRows];
        for(int i=1; i<=obj.QExtRows; i++){
            this->Ls[i-1]=obj.Ls[i-1];
        }
    }else{
        if(this->Ls!=NULL) delete [] this->Ls;
        this->Ls=NULL;
        this->IneRowsLength=obj.IneRowsLength;
    }
    this->QExtRows=obj.QExtRows;
}
void Array2DSize::Construct(){
    this->Ls=NULL;
    this->IneRowsLength=0;
    this->QExtRows=0;
}
void Array2DSize::SetNull(){
    if(this->Ls!=NULL){
        delete[]this->Ls;
    }
    this->Ls=NULL;
    this->IneRowsLength=0;
    this->QExtRows=0;
}
void Array2DSize::SetOne(bool isVaria){
    this->SetNull();
    this->QExtRows=1;
    if(isVaria==false){
        this->IneRowsLength=1;
    }else{
        this->Ls=new int(1);
        this->Ls[1-1]=1;
    }
}
bool Array2DSize::ExtRowNBelongsTo(int extRowN)const{return (extRowN>=1 && extRowN<=this->QExtRows);}
bool Array2DSize::PosBelongsTo(int extRowN, int ineRowN){
    return (this->ExtRowNBelongsTo(extRowN) && ineRowN>=1 && ineRowN<=this->GetLength(extRowN));
}

void Array2DSize::Set(int QExtRows, int IneRowsLength, bool preserveTypeIfVaria){
    if(this->Ls!=NULL){
        if(preserveTypeIfVaria){
            for(int i=1; i<=this->QExtRows; i++){
                this->Ls[i-1]=IneRowsLength;
            }
        }else{
            this->SetNull();
            this->IneRowsLength=IneRowsLength;
        }
    }else{
         this->IneRowsLength=IneRowsLength;
    }
    this->QExtRows=QExtRows;
}
void Array2DSize::Set(int QExtRows, int*Ls){
    if(Ls!=NULL){
        this->SetNull();
        this->Ls=new int[QExtRows];
        for(int i=1; i<=QExtRows; i++){
            this->Ls[i-1]=Ls[i-1];
        }
        this->QExtRows=QExtRows;
    }
}
void Array2DSize::Set1ExtRowVaria(int L){
    this->SetNull();
    this->QExtRows=1;
    this->Ls=new int[1];
    this->Ls[1-1]=L;
    this->IneRowsLength=0;
}
void Array2DSize::SetVaria(int QExtRows, int IneRowsLength){
    this->SetNull();
    this->QExtRows=QExtRows;
    this->IneRowsLength=0;
    this->Ls=new int[this->QExtRows];
    for(int i=1; i<=this->QExtRows; i++){
        this->Ls[i-1]=IneRowsLength;
    }
}
void Array2DSize::ConvertToRectIfAllLengthesAreEqual(){
    int mn=this->GetMinLength(), mx=this->GetMaxLength();
    if(mn=mx){
        delete[]this->Ls;
        this->Ls=NULL;
        this->IneRowsLength=mn;
    }
}

/*void Array2DSize::Add(int L, int ifConstAndNotEqual_Ignore0AddExisting1SetAll2){
    if(this->Ls!=NULL){
        Array1DAddElement(this->Ls, this->QExtRows, L);
    }else{
        if(ifConstAndNotEqual_Ignore0AddExisting1SetAll2==2){
            this->IneRowsLength=L;
        }
        if(ifConstAndNotEqual_Ignore0AddExisting1SetAll2==1 || ifConstAndNotEqual_Ignore0AddExisting1SetAll2==2){
            this->QExtRows+=1;
        }
    }
}
void Array2DSize::Ins(int N, int L){

}*/
void  Array2DSize::AddExt(int L){
    int LBefore=this->GetLength(this->GetQExtRows());
    if(L==0) L=LBefore;
    if(this->Ls!=NULL){
        Array1DAddElement(this->Ls, this->QExtRows, L);
    }else{
        if(L!=this->IneRowsLength){
            this->Ls=new int[this->QExtRows+1];
            for(int i=1; i<=this->QExtRows; i++){
                this->Ls[i-1]=this->IneRowsLength;
            }
            this->Ls[this->QExtRows+1-1]=L;
            this->IneRowsLength=0;
        }
        this->QExtRows+=1;
    }
}
void  Array2DSize::InsExt(int N, int L){
    int LBefore=this->GetLength(this->QExtRows);
    if(N==0) N=this->QExtRows;
    if(L==0) L=LBefore;
    if(this->Ls!=NULL){
        Array1DInsElementToN(this->Ls, this->QExtRows, N, L);
    }else{
        if(L!=this->IneRowsLength){
            this->Ls=new int[this->QExtRows+1];
            for(int i=1; i<=N-1; i++){
                this->Ls[i-1]=this->IneRowsLength;
            }
            this->Ls[N-1]=L;
            for(int i=N+1; i<=this->QExtRows; i++){
                this->Ls[i-1]=this->IneRowsLength;
            }
            this->Ls[this->QExtRows+1-1]=L;
            this->IneRowsLength=0;
        }
        this->QExtRows+=1;
    }
}

void Array2DSize::DelExt(int N){
    if(this->ExtRowNBelongsTo(N)){
        if(this->Ls!=NULL){
            Array1DDelElementFromN(this->Ls, this->QExtRows, N);
        }
        //this->QExtRows-=1;//no need: Del mut l'QE
        if(this->QExtRows==0)this->IneRowsLength=0;//better let it be
    }
}
//void  Array2DSize::DelIne(int N){
//    int minL=this->GetMinLength();
//    if(N<1) N=N+minL;
//    if(N>=1 && N<=minL){
//        if(this->Ls!=NULL){
//            for(int i=1; i<=this->QExtRows; i++){
//                this->Ls[i-1]-=1;
//            }
//        }else{
//            this->IneRowsLength-=1;
//        }
//    }
//}
void  Array2DSize::DelIneRow(){
    int minL=this->GetMinLength();
    if(this->Ls!=NULL){
        for(int i=1; i<=this->QExtRows; i++){
            this->Ls[i-1]-=1;
        }
    }else{
        this->IneRowsLength-=1;
    }
}
void  Array2DSize::AddOrInsIneRow(){
    if(this->Ls!=NULL){
        for(int i=1; i<=this->QExtRows; i++){
            this->Ls[i-1]+=1;
        }
    }else{
        this->IneRowsLength+=1;
    }
}
//void  Array2DSize::InsIne(int N){
//    int minL=this->GetMinLength();
//    if(N<1) N=N+minL;
//    if(N>=1 && N<=minL){
//        if(this->Ls!=NULL){
//            for(int i=1; i<=this->QExtRows; i++){
//                this->Ls[i-1]+=1;
//            }
//        }else{
//            this->IneRowsLength+=1;
//        }
//    }
//}
//void  Array2DSize::AddIne(){
//    int minL=this->GetMinLength();
//    if(N<1) N=N+minL;
//    if(N>=1 && N<=minL){
//        if(this->Ls!=NULL){
//            for(int i=1; i<=this->QExtRows; i++){
//                this->Ls[i-1]+=1;
//            }
//        }else{
//            this->IneRowsLength+=1;
//        }
//    }
//}

/*void Array2DSize::AddOrInsTo(int L, int N_0ToAllPreserveType_GTT0AndMakeVaria){
    int curL;
    if(this->ExtRowNBelongsTo(N_0ToAllPreserveType_GTT0AndMakeVaria)){
        if(N_0ToAllPreserveType_GTT0AndMakeVaria>0){
            curL=this->GetLength(N_0ToAllPreserveType_GTT0AndMakeVaria);
            this->SetLength((curL+L), N_0ToAllPreserveType_GTT0AndMakeVaria);
        }else{
            if(N_0ToAllPreserveType_GTT0AndMakeVaria==0){//ToAllPreserveType
                if(this->Ls!=NULL){
                    for(int i=1; i<=this->QExtRows; i++){
                        this->Ls[i-1]+=L;
                    }
                }else{
                    this->IneRowsLength+=L;
                }
            }
        }
    }
}
void Array2DSize::DelFrom(int L, int N_0ToAllPreserveType_GTT0AndMakeVaria){
    int curLwas, curLwillbe;
    bool b;
    if(this->ExtRowNBelongsTo(N_0ToAllPreserveType_GTT0AndMakeVaria)){
        curLwas=this->GetLength(N_0ToAllPreserveType_GTT0AndMakeVaria);
        curLwillbe=curLwas-L;
        if(curLwillbe>0){
            this->SetLength(curLwillbe, N_0ToAllPreserveType_GTT0AndMakeVaria);
        }
    }else{
        if(N_0ToAllPreserveType_GTT0AndMakeVaria==0){//ToAllPreserveType
            if(this->GetMinLength()-L>=0){
                if(this->Ls!=NULL){
                    for(int i=1; i<=this->QExtRows; i++){
                        this->Ls[i-1]-=L;
                    }
                }else{
                    this->IneRowsLength-=L;
                }
            }
        }
    }
}*/
void Array2DSize::AddOrInsTo(int N, int L){
    if(this->ExtRowNBelongsTo(N) && L>=0){
        if(this->Ls==NULL && L!=0){
            this->Ls=new int[this->QExtRows];
            for(int i=1; i<=this->QExtRows; i++){
                this->Ls[i-1]=this->IneRowsLength;
            }
            this->IneRowsLength=0;
        }
        this->Ls[N-1]+=L;
    }
}
void Array2DSize::DelFrom(int N, int L){
    if(this->ExtRowNBelongsTo(N) && L>=0){
        if(this->Ls==NULL && L!=0 && (this->GetLength(N-1)-L)>=0){
            this->Ls=new int[this->QExtRows];
            for(int i=1; i<=this->QExtRows; i++) this->Ls[i-1]=this->IneRowsLength;
            this->IneRowsLength=0;
        }
        this->Ls[N-1]+=L;
    }
}

int Array2DSize::GetQExtRows()const{ return this->QExtRows; }
bool Array2DSize::GetIfIsVariaLength()const{ return (this->Ls!=NULL); }
bool Array2DSize::isVariaLength(){ return (this->Ls!=NULL);}
bool Array2DSize::isRectangular()const{ return (this->Ls==NULL); }
int Array2DSize::GetLength(int N)const{
   int R=0;
   if(this->Ls!=NULL){
       if(N==0) N=this->QExtRows;
       if( this->ExtRowNBelongsTo(N)){
           R=this->Ls[N-1];
       }
   } else{
       R=this->IneRowsLength;
   }
   return R;
}
int Array2DSize::GetMinLength()const{
   int mn, mx;
    if(this->Ls!=NULL){
        mn=this->Ls[1-1];
        mx=this->Ls[1-1];
        for(int i=1; i<=this->QExtRows; i++){
            if(i==1 || (i>1 && this->Ls[i-1]<mn)) mn=this->Ls[i-1];
            if(i==1 || (i>1 && this->Ls[i-1]>mx)) mx=this->Ls[i-1];
        }
    }else{
        mn=this->IneRowsLength;
        mx=this->IneRowsLength;
    }
    return mn;
}

int Array2DSize::GetMaxLength()const{
    int mn, mx;
    if(this->Ls!=NULL){
        mn=this->Ls[1-1];
        mx=this->Ls[1-1];
        for(int i=1; i<=this->QExtRows; i++){
            if(i==1 || (i>1 && this->Ls[i-1]<mn)) mn=this->Ls[i-1];
            if(i==1 || (i>1 && this->Ls[i-1]>mx)) mx=this->Ls[i-1];
        }
    }else{
        mn=this->IneRowsLength;
        mx=this->IneRowsLength;
    }
    return mx;
 }

void Array2DSize::SwapExtRows(int N1, int N2){
    if(this->Ls!=NULL && this->ExtRowNBelongsTo(N1) && this->ExtRowNBelongsTo(N2)){
        Array1DSwapElements(this->Ls, this->QExtRows, N1, N2);
    }
}
void Array2DSize::Reverse(){
    if(this->Ls!=NULL)  Array1DReverse(this->Ls, this->QExtRows);
}
void Array2DSize::SetQExtRows(int N){
    int*R=NULL;
    int mn;
    if(N>=1){
        if(this->Ls!=NULL){
            mn=N>=this->QExtRows?N:this->QExtRows;
            R=new int[N];
            for(int i=1; i<=mn; i++){
                R[i-1]=this->Ls[i-1];
            }
            for(int i=mn+1; i<=N; i++){
                R[i-1]=1;
            }
            delete[]this->Ls;
            this->Ls=R;
        }
        this->QExtRows=N;
    }
}
void Array2DSize::SetIneRowsLength(int L_Minus1MinMinus2Max, bool setNotVaria){
    if(L_Minus1MinMinus2Max==-1) L_Minus1MinMinus2Max=this->GetMinLength();
    else if(L_Minus1MinMinus2Max==-2 )L_Minus1MinMinus2Max=this->GetMaxLength();
    if(L_Minus1MinMinus2Max>=1){
        if(this->Ls!=NULL && setNotVaria==false){
            for(int i=1; i<=this->QExtRows; i++){
                this->Ls[i-1]=L_Minus1MinMinus2Max;
            }
        }else{
            delete[]this->Ls;
            this->IneRowsLength=L_Minus1MinMinus2Max;
        }
    }
}
/*void Array2DSize::SetLength(int L, int N_0ToAllPreserveType_GTT0AndMakeVaria_Minus1ToAllMakeNotVaria){
    if(N_0ToAllPreserveType_GTT0AndMakeVaria_Minus1ToAllMakeNotVaria==-1){
        if(this->Ls==NULL){
            delete[]this->Ls;
            this->Ls=NULL;
            this->IneRowsLength=L;
        }
    }else if(N_0ToAllPreserveType_GTT0AndMakeVaria_Minus1ToAllMakeNotVaria==0){
        if(this->Ls==NULL){
            this->IneRowsLength=L;
        }else{
            for(int i=1; i<=this->QExtRows; i++){
                this->Ls[i-1]=L;
            }
        }
    }else if(N_0ToAllPreserveType_GTT0AndMakeVaria_Minus1ToAllMakeNotVaria>1 && this->ExtRowNBelongsTo(N_0ToAllPreserveType_GTT0AndMakeVaria_Minus1ToAllMakeNotVaria)){
        if(this->Ls==NULL){
            this->Ls=new int[this->QExtRows];
            for(int i=1; i<=this->QExtRows; i++){
                this->Ls[i-1]=this->IneRowsLength;
            }
            this->IneRowsLength=0;
        }
        this->Ls[N_0ToAllPreserveType_GTT0AndMakeVaria_Minus1ToAllMakeNotVaria-1]=L;
    }
}*/
void Array2DSize::SetLength(int L, int N){
    if(N<0){
        N=this->QExtRows+N;
    }
    if(this->Ls!=NULL){
        if(N==0){
            for(int i=1; i<=this->QExtRows; i++){
                this->Ls[i-1]=L;
            }
        }else if(this->ExtRowNBelongsTo(N)){
            this->Ls[N-1]=L;
        }
    }else{
        if(this->ExtRowNBelongsTo(N) && L!=this->IneRowsLength){
            this->Ls=new int[this->QExtRows];
            for(int i=1; i<=this->QExtRows; i++){
                this->Ls[i-1]=this->IneRowsLength;
            }
            this->Ls[N-1]=L;
            this->IneRowsLength=0;
        }else{
            this->IneRowsLength=L;
        }
    }
}

bool  Array2DSize::isSame(const Array2DSize&obj)const{
    bool verdict=true;
    //int Lmy;//
    //int Ltheir;//
    if(this->QExtRows!=obj.QExtRows) verdict=false;
    //if(this->GetQExtRows()!=obj.GetQExtRows()) verdict=false;
    //else if((this->Ls==NULL  && obj.Ls!=NULL)||(this->Ls==NULL  && obj.Ls!=NULL)){
    //     verdict=false;
    //}
    else if(this->Ls==NULL  && obj.Ls==NULL){
        if(this->IneRowsLength!=obj.IneRowsLength){
            verdict=false;
        }
    }else{
       for(int i=1; i<=this->QExtRows; i++){
           //Lmy=this->GetLength(i);
           //Ltheir=obj.GetLength(i);
           //Lmy=this->Ls[i-1];//no
           //Ltheir=obj.Ls[i-1];//no
           //if(Lmy!=Ltheir) verdict=false;
           if(this->GetLength(i)!=obj.GetLength(i)) verdict=false;
       }
    }
    return verdict;
}

QString  Array2DSize::GetAsQString(QString delim)const{
    QString c, s;
    c.setNum(this->QExtRows);
    s="Rows: ";
    s=s+c;
    s=s+" Lengthes: ";
    if(this->Ls!=NULL){
        s=s+"(";
        for(int i=1; i<=this->QExtRows-1; i++){
            c.setNum(this->GetLength(i));
            s=s+c;
            s=s+delim;
        }
        c.setNum(this->GetLength(this->QExtRows));
        s=s+c;
        s=s+")";
    }else{
        //c.setNum(this->GetLength(this->QExtRows));
        c.setNum(this->IneRowsLength);
        s=s+c;
    }
    return s;
}

std::string Array2DSize::GetAsStdString(std::string delim)const{
    std::string s2;
    QString delimer=QString(delim.c_str());
    QString s1;
    s1=this->GetAsQString(delimer);
    s2=s1.toStdString();
    return s2;
}

void Array2DSize::ShowToConsole()const{
    std::cout<<this->GetAsStdString()<<std::endl;
}

int Array2DSize::GetIneRowsLength(bool EvenIfVarButAllEqual)const{
    int Lmin=this->GetMinLength(), Lmax=this->GetMaxLength(), Lrect;
    if(this->Ls!=NULL && Lmin==Lmax && EvenIfVarButAllEqual){
        Lrect=Lmin;
    }else{
        Lrect=this->IneRowsLength;
    }
    return Lrect;
}
void Array2DSize::Transpose(){
    int iL=this->GetIneRowsLength(true), Lmax=this->GetMaxLength();
    int *Ls=NULL;
    if(this->IsForTranspose()){
        if(this->IneRowsLength>0){
            this->Set(this->IneRowsLength, this->QExtRows);
        }else{
            Ls=new int[Lmax];
            for(int i=1; i<=Lmax; i++){
                iL=this->LengthFullOfIneRowN(i);
                Ls[i-1]=iL;
            }
            //if(this->Ls!=NULL){
            //    delete[]this->Ls;
            //}
            this->Set(Lmax, Ls);
            if(Ls!=NULL){
                delete[]Ls;
            }
        }
    }
}
bool Array2DSize::IsForTranspose() const{
    bool b=true;
    int L1=this->GetLength(1), L2 ;
    for(int i=2; i<=this->QExtRows; i++){
        L2=this->GetLength(i);
        if(L1<L2)b=false;
        L1=L2;
    }
    return b;
}
int Array2DSize::LengthFullOfIneRowN(int N)const{//from start to 1st stop
    int L=0;
    if(N>=1 && N<=this->GetMaxLength()){
        for(int i=1; i<=this->QExtRows; i++){
            if(this->GetLength(i)>=N){
                L++;
            }
        }
    }
    return L;
}


int Array2DSize::GetPositionByCoordsIfRealizedAs1D(int ExtRowN, int IneRowN){
    int PosN=0;
    if(ExtRowN>=1 && ExtRowN<=this->GetQExtRows() && IneRowN>=1 && IneRowN<=this->GetLength(ExtRowN)){
        for(int i=1; i<=ExtRowN-1; i++){
            PosN+=this->GetLength(i);
        }
        PosN+=IneRowN;
    }
    return PosN;
}

void Array2DSize::GetCoordsByPositionIfRealizedAs1D(int PosN, int&ExtRowN, int&IneRowN){
    int QElementsInAll=0, QExtRows=this->GetQExtRows(), curN=0, curL;
    for(int i=1; i<=QExtRows; i++){
        QElementsInAll+=this->GetLength(i);
    }
    ExtRowN=0;
    IneRowN=0;
    bool contin=(PosN>=1 && PosN<=QElementsInAll);
    while(contin){
        ExtRowN++;
        curL=this->GetLength(ExtRowN);
        if(curN+curL<PosN){
            curN+=curL;
        }else{
            IneRowN=PosN-curN;
            contin=false;
        }
    }
}

int Array2DSize::GetExtRowNByPositionIfRealizedAs1D(int PosN){
    int ExtRowN=0, IneRowN=0;
    this->GetCoordsByPositionIfRealizedAs1D(PosN, ExtRowN, IneRowN);
    return ExtRowN;
}

int Array2DSize::Array2DSize::GetIneRowNByPositionIfRealizedAs1D(int PosN){
    int ExtRowN=0, IneRowN=0;
    this->GetCoordsByPositionIfRealizedAs1D(PosN, ExtRowN, IneRowN);
    return IneRowN;
}

std::vector<int>Array2DSize::GetLengthes()const{
    std::vector<int>Ls;
    if(this->Ls!=NULL){
        for(int i=1; i<=this->QExtRows; i++){
            Ls.push_back(this->Ls[i-1]);
        }
    }else{
        for(int i=1; i<=this->QExtRows; i++){
            Ls.push_back(this->IneRowsLength);//wa this-IneRowsLength ma n'wa err!
        }
    }
    return Ls;
}
int Array2DSize::GetQElements(){
    int count=0, L;
    if(this->Ls!=NULL){
        for(int i=1; i<=this->QExtRows; i++){
            L=this->GetLength(i);
            count+=L;
        }
    }else{
        count=this->QExtRows*this->IneRowsLength;
    }
    return count;
}

//=======================================================================================================

/*

//Array2DSizePrototype*Array2DSizePrototype::clone() const{
//    return new Array2DSizePrototype(*this);
//}
//
bool Array2DSizePrototype::IsRectangularFact() const{
    return(this->GetMinLength()==this->GetMaxLength());
}
bool Array2DSizePrototype::IsVariableLengthFact() const{
    return(this->GetMinLength()!=this->GetMaxLength());
}



//class Array2DSizeVarLenP: public Array2DSizePrototype
//{
//    int*Ls;
//    int QExtRows;
//public:
    //virtual Array2DSize();
    //Array2DSizeVarLenP(const Array2DSizeVarLenP&obj);//I men no not ob in dyn polymorf S nablb realf'd
    Array2DSizeVarLenP::Array2DSizeVarLenP(int QExtRows, int ExtRowsLength): Array2DSizePrototype(QExtRows, ExtRowsLength){
        this->Construct();
        this->Set(QExtRows, ExtRowsLength);
    }
    Array2DSizeVarLenP::Array2DSizeVarLenP(int QExtRows, int*Ls) :Array2DSizePrototype(){//{//: Array2DSizePrototype(QExtRows, Ls){
        this->Construct();
        this->Set(QExtRows, Ls);
    }
    Array2DSizeVarLenP::Array2DSizeVarLenP(std::vector<int>Ls) :Array2DSizePrototype(){//{//: Array2DSizePrototype(Ls){
        this->Construct();
        this->Set(Ls);
    }
    //Array2DSizeVarLenP operator =(const Array2DSizePrototype&obj);
    //Array2DSizeVarLenP* Array2DSizeVarLenP::clone() const {
    //    return this;
    //}
    Array2DSizeVarLenP::~Array2DSizeVarLenP(){
        delete[]this->Ls;
    }
    //virtual void Assign(const Array2DSizePrototype&obj);//I men no not ob in dyn polymorf S nablb realf'd
    //void Array2DSizeVarLenP::Assign(Array2DSizeVarLenP*obj){//I hope ce abls ob Array2DSizeVarLenP s'heir of Array2DSizePrototype]
    //    this->SetNull();
    //    this->Set(obj->QExtRows, obj->Ls);
    //}
    void Array2DSizeVarLenP::Construct(){
        this->QExtRows=0;
        this->Ls=NULL;
    }

    void Array2DSizeVarLenP::SetNull(){
        if(this->Ls!=NULL){
            delete this->Ls;
        }
        this->QExtRows=0;
    }
    //
    bool Array2DSizeVarLenP::IsRectangularType() const{
        return false;
    }
    bool Array2DSizeVarLenP::IsVariableLengthType() const{
        return true;
    }
    //bool Array2DSizeVarLenP::IsRectangularFact() const{
    //    int Lmin=this->GetMinLength(), Lmax=this->GetMaxLength();
    //    return (Lmin==Lmax);
    //}
    //bool Array2DSizeVarLenP::IsVariableLengthFact() const{
    //    int Lmin=this->GetMinLength(), Lmax=this->GetMaxLength();
    //    return (Lmin!=Lmax);
    //}
    //
    //virtual void SetVaria();
    //virtual void SetRectangular(int L_ifVaria_min_m1_max_less=-2);
    //
    void Array2DSizeVarLenP::SetOne(){
        this->SetNull();
        this->Set(1, 1);
    }
    bool Array2DSizeVarLenP::ExtRowNBelongsTo(int extRowN) const{
        bool verdict;
        verdict=(extRowN>=1 && extRowN<=this->QExtRows);
        return verdict;
    }
    bool Array2DSizeVarLenP::PosBelongsTo(int extRowN, int ineRowN){
        int L;
        bool b;
        b=(this->ExtRowNBelongsTo(extRowN));
        if(b){
           L=this->GetLength(extRowN);
           b=(ineRowN>=1 && ineRowN<=L);
        }
        return b;
    }
    void Array2DSizeVarLenP::Set(int QExtRows, int IneRowsLength){
        this->SetNull();
        this->QExtRows=QExtRows;
        this->Ls=new int[this->QExtRows];
        for(int i=1; i<=this->QExtRows; i++){
            this->Ls[i-1]=IneRowsLength;
        }
    }
    void Array2DSizeVarLenP::Set(int QExtRows,int*Ls){
        this->SetNull();
        this->QExtRows=QExtRows;
        this->Ls=new int[this->QExtRows];
        if(this->Ls!=NULL){
            for(int i=1; i<this->QExtRows; i++){
                this->Ls[i-1]=Ls[i-1];
            }
        }else{
            for(int i=1; i<this->QExtRows; i++){
                this->Ls[i-1]=0;
            }
        }
    }
    void Array2DSizeVarLenP::Set(std::vector<int>Ls){
        this->SetNull();
        this->QExtRows=Ls.size();
        for(int i=1; i<this->QExtRows; i++){
            this->Ls[i-1]=Ls[i-1];
        }
    }
    //virtual void Set1ExtRowVaria(int L);
    //virtual void SetVaria(int QExtRows, int IneRowsLength);
    //virtual void ConvertToRectIfAllLengthesAreEqual();
    //virtual void Add(int L, int ifConstAndNotEqual_Ignore0AddExisting1SetAll2=2);
    //virtual void Ins(int N, int L);
    void Array2DSizeVarLenP::AddExt(int L){
        Array1DAddElement(this->Ls, this->QExtRows, L);
    }
    void  Array2DSizeVarLenP::InsExt(int N, int L){
        Array1DInsElementToN(this->Ls, this->QExtRows, N, L);
    }
    void Array2DSizeVarLenP::DelExt(int N){
        Array1DDelElementFromN(this->Ls, this->QExtRows, N);
    }
    void Array2DSizeVarLenP::DelIneRow(){
        for(int i=1; i<=this->QExtRows; i++){
            this->Ls[i-1]--;
        }
    }
    void Array2DSizeVarLenP::AddOrInsIneRow(){
        for(int i=1; i<=this->QExtRows; i++){
            this->Ls[i-1]++;
        }
    }
    //virtual void AddOrInsTo(int L=0, int N_0ToAllPreserveType_GTT0AndMakeVaria=0);
    //virtual void DelFrom(int L, int N_0ToAllPreserveType_GTT0AndMakeVaria=0);
    void Array2DSizeVarLenP::AddOrInsTo(int N, int L){
        if(N>=1 && N<=this->GetQExtRows() && L>=1){
            this->Ls+=L;
        }
    }
    void  Array2DSizeVarLenP::DelFrom(int N, int L){
        if(N>=1 && N<=this->GetQExtRows() && L>=1){
            this->Ls-=L;
        }
    }
    int Array2DSizeVarLenP::GetQExtRows() const{
        return this->QExtRows;
    }
    //virtual bool GetIfIsVariaLength() const;
    //virtual bool Array2DSizeVarLenP::isVariaLength() const{return true; }
    //virtual bool isRectangular() const;
    int Array2DSizeVarLenP::GetLength(int N) const{
        int L=0, minL, maxL, curL, QExtRows=this->GetQExtRows();
        if(N>=1 && N<=this->GetQExtRows()){
            L=this->Ls[N-1];
        }else if(N==0){
            minL=0;
            maxL=0;
            for(int i=1; i<=QExtRows; i++){
                curL=Ls[i-1];
                if(i==1 || (i>1 && curL<minL)){
                    minL=curL;
                }
                if(i==1 || (i>1 && curL>maxL)){
                    maxL=curL;
                }
            }
            if(minL==maxL){
                L=this->Ls[1-1];
            }//else remains 0
        }
        return L;
    }
    int Array2DSizeVarLenP::GetMinLength()const{
        int L=0, minL, maxL, curL, QExtRows=this->GetQExtRows();
        for(int i=1; i<=QExtRows; i++){
            curL=Ls[i-1];
            if(i==1 || (i>1 && curL<minL)){
                minL=curL;
            }
            if(i==1 || (i>1 && curL>maxL)){
                maxL=curL;
            }
        }
        return minL;
    }
    int Array2DSizeVarLenP::GetMaxLength()const{
        int L=0, minL, maxL, curL, QExtRows=this->GetQExtRows();
        for(int i=1; i<=QExtRows; i++){
            curL=Ls[i-1];
            if(i==1 || (i>1 && curL<minL)){
                minL=curL;
            }
            if(i==1 || (i>1 && curL>maxL)){
                maxL=curL;
            }
        }
        return maxL;
    }
    void Array2DSizeVarLenP::SwapExtRows(int N1, int N2){
        int buf, QExtRows=this->GetQExtRows();
        if(N1>=1 && N1<=QExtRows && N2>=1 && N2<=QExtRows){
            buf=this->Ls[N1-1];
            this->Ls[N1-1]=this->Ls[N2-1];
            this->Ls[N2-1]=buf;
        }
    }
    void Array2DSizeVarLenP::Reverse(){
        Array1DReverse(this->Ls, this->QExtRows);
    }
    //virtual void SetQExtRows(int N);
    //virtual void SetIneRowsLength(int L_Minus1MinMinus2Max=-1, bool preserveTypeIfVaria=false);
    //virtual void SetLength(int L, int N_0ToAllPreserveType_GTT0AndMakeVaria_Minus1ToAllMakeNotVaria=0);
    void Array2DSizeVarLenP::SetLength(int L, int N){
        int QExtRows=this->GetQExtRows();
        if(N<0){
            N=QExtRows+N+1;
        }
        if(N>=1 && N<=QExtRows){
            this->Ls[N-1]=L;
        }else if(N==0){
            for(int i=1; i<=QExtRows; i++){
                this->Ls[i-1]=L;
            }
        }
    }
    QString Array2DSizeVarLenP::GetAsQString(QString delim)const{
       QString s, sn;
       int QExtRows=this->GetQExtRows();
       sn.setNum(QExtRows);
       s="Size: "+sn+": [";
       for(int i=1; i<=QExtRows-1; i++){
          sn.setNum(this->Ls[i-1]);
          s+=sn;
          s+=delim;
       }
       sn.setNum(this->Ls[QExtRows-1]);
       s+=sn;
       s+="]";
       return s;
    }
    std::string Array2DSizeVarLenP::GetAsStdString(std::string delim) const{
        QString qsd=QString(delim.c_str());
        QString qs=this->GetAsQString(qsd);
        std::string ss=qs.toStdString();
        return ss;
    }
    void  Array2DSizeVarLenP::ShowToConsole() const{
        std::string ss=this->GetAsStdString();
        std::cout<<ss;
    }
    //virtual bool isSame(const Array2DSizeVarLenP&obj)const;
    //bool Array2DSizeVarLenP::isSame(Array2DSizeVarLenP*obj)const{
    //    bool b;
    //    int QExtRows=this->GetQExtRows();
    //    b=(this->GetQExtRows()==obj->GetQExtRows());
    //    if(b){
    //        for(int i=1; i<=QExtRows; i++){
    //            if(this->Ls[i-1]!=obj->Ls[i-1]){
    //                b=false;
    //            }
    //        }
    //    }
    //    return b;
    //}
    int Array2DSizeVarLenP::GetIneRowsLength()const{ return this->GetLength(0); }
    int Array2DSizeVarLenP::LengthFullOfIneRowN(int IneRowN) const{//from start to 1st stop
        int count, QExtRows=this->GetQExtRows(), Lmax=this->GetMaxLength(), Lcur, N=0;
        bool breakreached=false, contin=(N>=1 && N<=Lmax);
        while(contin){
            N++;
            Lcur=this->GetLength(N);
            if(Lcur>=IneRowN){
                if(breakreached==false){
                    count++;
                }else{
                    count=0;
                    contin=false;
                }
            }else{
                breakreached=true;
            }
        }
        return count;
    }
    bool Array2DSizeVarLenP::IsForTranspose() const{
        bool b=true;
        int count, QExtRows=this->GetQExtRows(), Lmax=this->GetMaxLength(), Lcur, Lnxt, N=0;
        for(int i=1; i<=QExtRows-1; i++){
            Lcur=this->LengthFullOfIneRowN(i);
            Lnxt=this->LengthFullOfIneRowN(i+1);
            if(Lcur>Lnxt){
                b=false;
            }
        }
        return b;
    }
    void Array2DSizeVarLenP::Transpose(){
        int*Ls=NULL, Lmax=this->GetMaxLength();
        if(this->IsForTranspose()){
            Ls=new int[Lmax];
            for(int i=1; i<=Lmax; i++){
                Ls[i-1]=this->LengthFullOfIneRowN(i);
            }
            this->Set(Lmax, Ls);
        }
    }
//};

    //class Array2DSizeVarLenV: public Array2DSizePrototype
    //{
    //    int*Ls;
    //    int QExtRows;
    //public:
        //virtual Array2DSize();
        //Array2DSizeVarLenV(const Array2DSizeVarLenV&obj);//I men no not ob in dyn polymorf S nablb realf'd
    Array2DSizeVarLenV::Array2DSizeVarLenV(int QExtRows, int ExtRowsLength): Array2DSizePrototype(){//{// : Array2DSizePrototype(QExtRows, ExtRowsLength){
            this->Construct();
            this->Set(QExtRows, ExtRowsLength);
        }
        Array2DSizeVarLenV::Array2DSizeVarLenV(int QExtRows, int*Ls){// : Array2DSizePrototype(QExtRows, Ls){
            this->Construct();
            this->Set(QExtRows, Ls);
        }
        Array2DSizeVarLenV::Array2DSizeVarLenV(std::vector<int>Ls){// : Array2DSizePrototype(Ls){
            this->Construct();
            this->Set(Ls);
        }
        //Array2DSizeVarLenV operator =(const Array2DSizePrototype&obj);
        //Array2DSizeVarLenV* Array2DSizeVarLenV::clone() const {
        //    return this;
        //}
        Array2DSizeVarLenV::~Array2DSizeVarLenV(){}
        //virtual void Assign(const Array2DSizePrototype&obj);//I men no not ob in dyn polymorf S nablb realf'd
        //void Array2DSizeVarLenV::Assign(Array2DSizeVarLenV*obj){//I hope ce abls ob Array2DSizeVarLenV s'heir of Array2DSizePrototype]
        //    this->SetNull();
        //    this->Set(obj->Ls);
        //}
        void Array2DSizeVarLenV::Construct(){
            this->SetNull();
        }

        void Array2DSizeVarLenV::SetNull(){
            this->Ls.clear();
        }
        //
        bool Array2DSizeVarLenV::IsRectangularType() const{
            return false;
        }
        bool Array2DSizeVarLenV::IsVariableLengthType() const{
            return true;
        }
        bool Array2DSizeVarLenV::IsRectangularFact() const{
            int Lmin=this->GetMinLength(), Lmax=this->GetMaxLength();
            return (Lmin==Lmax);
        }
        bool Array2DSizeVarLenV::IsVariableLengthFact() const{
            int Lmin=this->GetMinLength(), Lmax=this->GetMaxLength();
            return (Lmin!=Lmax);
        }
        //
        //virtual void SetVaria();
        //virtual void SetRectangular(int L_ifVaria_min_m1_max_less=-2);
        //
        void Array2DSizeVarLenV::SetOne(){
            this->SetNull();
            this->Set(1, 1);
        }
        bool Array2DSizeVarLenV::ExtRowNBelongsTo(int extRowN)const{
            bool verdict;
            verdict=(extRowN>=1 && extRowN<=this->GetQExtRows());
            return verdict;
        }
        bool Array2DSizeVarLenV::PosBelongsTo(int extRowN, int ineRowN){
            int L;
            bool b;
            b=(this->ExtRowNBelongsTo(extRowN));
            if(b){
               L=this->GetLength(extRowN);
               b=(ineRowN>=1 && ineRowN<=L);
            }
            return b;
        }
        void Array2DSizeVarLenV::Set(int QExtRows, int ExtRowsLength){
            this->SetNull();
            for(int i=1; i<=QExtRows; i++){
                this->Ls.push_back(ExtRowsLength);
            }
        }
        void Array2DSizeVarLenV::Set(int QExtRows, int*Ls){
            this->SetNull();
            if(Ls!=NULL){
                for(int i=1; i<QExtRows; i++){
                    this->Ls.push_back(Ls[i-1]);
                }
            }else{
                for(int i=1; i<QExtRows; i++){
                    this->Ls.push_back(0);
                }
            }
        }
        void Array2DSizeVarLenV::Set(std::vector<int>Ls){
            this->Ls=Ls;
        }
        //virtual void Set1ExtRowVaria(int L);
        //virtual void SetVaria(int QExtRows, int IneRowsLength);
        //virtual void ConvertToRectIfAllLengthesAreEqual();
        //virtual void Add(int L, int ifConstAndNotEqual_Ignore0AddExisting1SetAll2=2);
        //virtual void Ins(int N, int L);
        void Array2DSizeVarLenV::AddExt(int L){
            Array1DAddElement(this->Ls, L);
        }
        void  Array2DSizeVarLenV::InsExt(int N, int L){
            Array1DInsElementToN(this->Ls, N, L);
        }
        void Array2DSizeVarLenV::DelExt(int N){
            Array1DDelElementFromN(this->Ls, N);
        }
        void Array2DSizeVarLenV::DelIneRow(){
            int QExtRows=this->GetQExtRows();
            for(int i=1; i<=QExtRows; i++){
                this->Ls[i-1]--;
            }
        }
        void Array2DSizeVarLenV::AddOrInsIneRow(){
            int QExtRows=this->GetQExtRows();
            for(int i=1; i<=QExtRows; i++){
                this->Ls[i-1]++;
            }
        }
        //virtual void AddOrInsTo(int L=0, int N_0ToAllPreserveType_GTT0AndMakeVaria=0);
        //virtual void DelFrom(int L, int N_0ToAllPreserveType_GTT0AndMakeVaria=0);
        void Array2DSizeVarLenV::AddOrInsTo(int N, int L){
            if(N>=1 && N<=this->GetQExtRows() && L>=1){
                this->Ls[N-1]+=L;
            }
        }
        void  Array2DSizeVarLenV::DelFrom(int N, int L){
            if(N>=1 && N<=this->GetQExtRows() && L>=1){
                this->Ls[N-1]-=L;
            }
        }
        int Array2DSizeVarLenV::GetQExtRows() const{
            return this->Ls.size();
        }
        //virtual bool GetIfIsVariaLength() const;
        //virtual bool Array2DSizeVarLenV::isVariaLength() const{return true; }
        //virtual bool isRectangular() const;
        int Array2DSizeVarLenV::GetLength(int N) const{
            int L=0, minL, maxL, curL, QExtRows=this->GetQExtRows();
            if(N>=1 && N<=this->GetQExtRows()){
                L=this->Ls[N-1];
            }else if(N==0){
                minL=0;
                maxL=0;
                for(int i=1; i<=QExtRows; i++){
                    curL=Ls[i-1];
                    if(i==1 || (i>1 && curL<minL)){
                        minL=curL;
                    }
                    if(i==1 || (i>1 && curL>maxL)){
                        maxL=curL;
                    }
                }
                if(minL==maxL){
                    L=this->Ls[1-1];
                }//else remains 0
            }
            return L;
        }
        int Array2DSizeVarLenV::GetMinLength()const{
            int L=0, minL, maxL, curL, QExtRows=this->GetQExtRows();
            for(int i=1; i<=QExtRows; i++){
                curL=Ls[i-1];
                if(i==1 || (i>1 && curL<minL)){
                    minL=curL;
                }
                if(i==1 || (i>1 && curL>maxL)){
                    maxL=curL;
                }
            }
            return minL;
        }
        int Array2DSizeVarLenV::GetMaxLength()const{
            int L=0, minL, maxL, curL, QExtRows=this->GetQExtRows();
            for(int i=1; i<=QExtRows; i++){
                curL=Ls[i-1];
                if(i==1 || (i>1 && curL<minL)){
                    minL=curL;
                }
                if(i==1 || (i>1 && curL>maxL)){
                    maxL=curL;
                }
            }
            return maxL;
        }
        void Array2DSizeVarLenV::SwapExtRows(int N1, int N2){
            int buf, QExtRows=this->GetQExtRows();
            if(N1>=1 && N1<=QExtRows && N2>=1 && N2<=QExtRows){
                buf=this->Ls[N1-1];
                this->Ls[N1-1]=this->Ls[N2-1];
                this->Ls[N2-1]=buf;
            }
        }
        void Array2DSizeVarLenV::Reverse(){
            Array1DReverse(this->Ls);
        }
        //virtual void SetQExtRows(int N);
        //virtual void SetIneRowsLength(int L_Minus1MinMinus2Max=-1, bool preserveTypeIfVaria=false);
        //virtual void SetLength(int L, int N_0ToAllPreserveType_GTT0AndMakeVaria_Minus1ToAllMakeNotVaria=0);
        void Array2DSizeVarLenV::SetLength(int L, int N){
            int QExtRows=this->GetQExtRows();
            if(N<0){
                N=QExtRows+N+1;
            }
            if(N>=1 && N<=QExtRows){
                this->Ls[N-1]=L;
            }else if(N==0){
                for(int i=1; i<=QExtRows; i++){
                    this->Ls[i-1]=L;
                }
            }
        }
        QString Array2DSizeVarLenV::GetAsQString(QString delim)const{
           QString s, sn;
           int QExtRows=this->GetQExtRows();
           sn.setNum(QExtRows);
           s="Size: "+sn+": [";
           for(int i=1; i<=QExtRows-1; i++){
              sn.setNum(this->Ls[i-1]);
              s+=sn;
              s+=delim;
           }
           sn.setNum(this->Ls[QExtRows-1]);
           s+=sn;
           s+="]";
           return s;
        }
        std::string Array2DSizeVarLenV::GetAsStdString(std::string delim) const{
            QString qsd=QString(delim.c_str());
            QString qs=this->GetAsQString(qsd);
            std::string ss=qs.toStdString();
            return ss;
        }
        void  Array2DSizeVarLenV::ShowToConsole() const{
            std::string ss=this->GetAsStdString();
            std::cout<<ss;
        }
        //virtual bool isSame(const Array2DSizeVarLenV&obj)const;
        bool Array2DSizeVarLenV::isSame(Array2DSizeVarLenV*obj)const{
            bool b;
            int QExtRows=this->GetQExtRows();
            b=(this->GetQExtRows()==obj->GetQExtRows());
            if(b){
                for(int i=1; i<=QExtRows; i++){
                    if(this->Ls[i-1]!=obj->Ls[i-1]){
                        b=false;
                    }
                }
            }
            return b;
        }
        int Array2DSizeVarLenV::GetIneRowsLength()const{ return this->GetLength(0); }
        int Array2DSizeVarLenV::LengthFullOfIneRowN(int IneRowN) const{//from start to 1st stop
            int count, QExtRows=this->GetQExtRows(), Lmax=this->GetMaxLength(), Lcur, N=0;
            bool breakreached=false, contin=(N>=1 && N<=Lmax);
            while(contin){
                N++;
                Lcur=this->GetLength(N);
                if(Lcur>=IneRowN){
                    if(breakreached==false){
                        count++;
                    }else{
                        count=0;
                        contin=false;
                    }
                }else{
                    breakreached=true;
                }
            }
            return count;
        }
        bool  Array2DSizeVarLenV::IsForTranspose() const{
            bool b=true;
            int count, QExtRows=this->GetQExtRows(), Lmax=this->GetMaxLength(), Lcur, Lnxt, N=0;
            for(int i=1; i<=QExtRows-1; i++){
                Lcur=this->LengthFullOfIneRowN(i);
                Lnxt=this->LengthFullOfIneRowN(i+1);
                if(Lcur>Lnxt){
                    b=false;
                }
            }
            return b;
        }
        void Array2DSizeVarLenV::Transpose(){
            int*Ls=NULL, Lmax=this->GetMaxLength();
            if(this->IsForTranspose()){
                Ls=new int[Lmax];
                for(int i=1; i<=Lmax; i++){
                    Ls[i-1]=this->LengthFullOfIneRowN(i);
                }
                this->Set(Lmax, Ls);
            }
        }
    //};


//class Array2DSizeRect: public Array2DSizePrototype
//{
//    int QExtRows;
//    int ExtRowLength;
//public:
    //virtual Array2DSize();
        Array2DSizeRect::Array2DSizeRect(int QExtRows, int ExtRowsLength){// : Array2DSizePrototype(QExtRows, ExtRowsLength){
            this->Construct();
            this->Set(QExtRows, ExtRowsLength);
        }
        Array2DSizeRect::Array2DSizeRect(int QExtRows, int*Ls){// : Array2DSizePrototype(QExtRows, Ls){
            this->Construct();
            this->Set(QExtRows, Ls[1-1]);
        }
        Array2DSizeRect::Array2DSizeRect(std::vector<int>Ls){// : Array2DSizePrototype(Ls) {
            this->Construct();
            if(Ls.size()>0){
                this->Set(QExtRows, Ls[1-1]);
            }
        }
        //Array2DSizeRect operator =(const Array2DSizePrototype&obj);
        //Array2DSizeRect* Array2DSizeRect::clone() const {
        //    return this;
        //}
        Array2DSizeRect::~Array2DSizeRect(){}
        //virtual void Assign(const Array2DSizePrototype&obj);//I men no not ob in dyn polymorf S nablb realf'd
        void Array2DSizeRect::Assign(Array2DSizeRect*obj){//I hope ce abls ob Array2DSizeRect s'heir of Array2DSizePrototype]
            this->SetNull();
            this->Set(obj->QExtRows, obj->ExtRowLength);
        }
        void Array2DSizeRect::Construct(){
            this->SetNull();
        }

        void Array2DSizeRect::SetNull(){
            this->QExtRows=0;
            this->ExtRowLength=0;
        }
        //
        bool Array2DSizeRect::IsRectangularType() const{
            return true;
        }
        bool Array2DSizeRect::IsVariableLengthType() const{
            return false;
        }
        bool Array2DSizeRect::IsRectangularFact() const{
            int Lmin=this->GetMinLength(), Lmax=this->GetMaxLength();
            return (Lmin==Lmax);
        }
        bool Array2DSizeRect::IsVariableLengthFact() const{
            int Lmin=this->GetMinLength(), Lmax=this->GetMaxLength();
            return (Lmin!=Lmax);
        }
        //
        //virtual void SetVaria();
        //virtual void SetRectangular(int L_ifVaria_min_m1_max_less=-2);
        //
        void Array2DSizeRect::SetOne(){
            this->SetNull();
            this->Set(1, 1);
        }
        bool Array2DSizeRect::ExtRowNBelongsTo(int extRowN)const{
            bool verdict;
            verdict=(extRowN>=1 && extRowN<=this->GetQExtRows());
            return verdict;
        }
        bool Array2DSizeRect::PosBelongsTo(int extRowN, int ineRowN){
            int L;
            bool b;
            b=(this->ExtRowNBelongsTo(extRowN));
            if(b){
               L=this->ExtRowLength;
               b=(ineRowN>=1 && ineRowN<=L);
            }
            return b;
        }
        void Array2DSizeRect::Set(int QExtRows, int ExtRowLength){
            this->QExtRows=QExtRows;
            this->ExtRowLength=ExtRowLength;
        }
        void Array2DSizeRect::Set(int QExtRows, int*Ls){
            this->QExtRows=QExtRows;
            if(Ls!=NULL){
                this->ExtRowLength=Ls[1-1];
            }
        }
        void Array2DSizeRect::Set(std::vector<int>Ls){
            this->QExtRows=Ls.size();
            if(this->QExtRows>0){
                this->ExtRowLength=Ls[1-1];
            }else{
                this->ExtRowLength=0;
            }
        }
        //virtual void Set1ExtRowVaria(int L);
        //virtual void SetVaria(int QExtRows, int IneRowsLength);
        //virtual void ConvertToRectIfAllLengthesAreEqual();
        //virtual void Add(int L, int ifConstAndNotEqual_Ignore0AddExisting1SetAll2=2);
        //virtual void Ins(int N, int L);
        void Array2DSizeRect::AddExt(int L){
            this->ExtRowLength++;
        }
        void  Array2DSizeRect::InsExt(int N, int L){
            this->ExtRowLength++;
        }
        void Array2DSizeRect::DelExt(int N){
            this->ExtRowLength--;
        }
        void Array2DSizeRect::DelIneRow(){
            this->ExtRowLength--;
        }
        void Array2DSizeRect::AddOrInsIneRow(){
             this->ExtRowLength++;
        }
        //virtual void AddOrInsTo(int L=0, int N_0ToAllPreserveType_GTT0AndMakeVaria=0);
        //virtual void DelFrom(int L, int N_0ToAllPreserveType_GTT0AndMakeVaria=0);
        void Array2DSizeRect::AddOrInsTo(int N, int L){
            //NOp;
        }
        void Array2DSizeRect::DelFrom(int N, int L){
            //NOp;
        }
        int Array2DSizeRect::GetQExtRows() const{
            return this->QExtRows;
        }
        //virtual bool GetIfIsVariaLength() const;
        //virtual bool Array2DSizeRect::isVariaLength() const{return true; }
        //virtual bool isRectangular() const;
        int Array2DSizeRect::GetLength(int N) const{
            return this->ExtRowLength;
        }
        int Array2DSizeRect::GetMinLength()const{
            return this->ExtRowLength;
        }
        int Array2DSizeRect::GetMaxLength()const{
            return this->ExtRowLength;
        }
        void Array2DSizeRect::SwapExtRows(int N1, int N2){
            //NOp;
        }
        void Array2DSizeRect::Reverse(){
            //NOp;
        }
        //virtual void SetQExtRows(int N);
        //virtual void SetIneRowsLength(int L_Minus1MinMinus2Max=-1, bool preserveTypeIfVaria=false);
        //virtual void SetLength(int L, int N_0ToAllPreserveType_GTT0AndMakeVaria_Minus1ToAllMakeNotVaria=0);
        void Array2DSizeRect::SetLength(int L, int N){
            this->ExtRowLength=L;
        }
        QString Array2DSizeRect::GetAsQString(QString delim)const{
           QString s, sn;
           int QExtRows=this->GetQExtRows();
           sn.setNum(QExtRows);
           s="Size: E="+sn+" x I=";

           sn.setNum(this->ExtRowLength);
           s+=sn;
           return s;
        }
        std::string Array2DSizeRect::GetAsStdString(std::string delim) const{
            QString qsd=QString(delim.c_str());
            QString qs=this->GetAsQString(qsd);
            std::string ss=qs.toStdString();
            return ss;
        }
        void  Array2DSizeRect::ShowToConsole() const{
            std::string ss=this->GetAsStdString();
            std::cout<<ss;
        }
        //virtual bool isSame(const Array2DSizeRect&obj)const;
        bool Array2DSizeRect::isSame(Array2DSizeRect*obj)const{
            return (this->QExtRows == obj->QExtRows && this->ExtRowLength == obj->ExtRowLength);
        }
        int Array2DSizeRect::GetIneRowsLength()const{ return this->GetLength(0); }
        int Array2DSizeRect::LengthFullOfIneRowN(int IneRowN) const{//from start to 1st stop
            return this->QExtRows;
        }
        bool  Array2DSizeRect::IsForTranspose() const{
            return true;
        }
        void Array2DSizeRect::Transpose(){
            this->Set(this->ExtRowLength, this->QExtRows);
        }
    //};
//};




//class Array2DSizeDPS //Dynamic Polymorph Shell
//{
//    Array2DSizePrototype*size;
//public:
    //virtual Array2DSize();
        Array2DSizeDPS::Array2DSizeDPS(int QExtRows, int ExtRowsLength, bool isNotVaria){
            if(isNotVaria){
                this->size=new Array2DSizeRect(QExtRows, ExtRowsLength);
            }else{
                this->size=new Array2DSizeVarLenV(QExtRows, ExtRowsLength);
            }
            this->size=new Array2DSizeVarLenP(QExtRows, ExtRowsLength);
        }

        Array2DSizeDPS::Array2DSizeDPS(int QExtRows, int*Ls, bool RectIfAllSame){
            int Lmin, Lmax, Lcur;
            for(int i=1; i<=QExtRows; i++){
                if(i==1 || (i>1 && Lcur<Lmin)){
                    Lmin=Lcur;
                }
                if(i==1 || (i>1 && Lcur<Lmax)){
                    Lmax=Lcur;
                }
            }
            if(Lmin==Lmax && RectIfAllSame){
                this->size=new Array2DSizeRect(QExtRows, Lmax);
            }else{
                this->size=new Array2DSizeVarLenP(QExtRows, Ls);
            }
        }
        Array2DSizeDPS::Array2DSizeDPS(std::vector<int>Ls, bool RectIfAllSame){
            int QExtRows=Ls.size(), Lmin, Lmax, Lcur;
            for(int i=1; i<=QExtRows; i++){
                if(i==1 || (i>1 && Lcur<Lmin)){
                    Lmin=Lcur;
                }
                if(i==1 || (i>1 && Lcur<Lmax)){
                    Lmax=Lcur;
                }
            }
            if(Lmin==Lmax && RectIfAllSame){
                this->size=new Array2DSizeRect(QExtRows, Lmax);
            }else{
                this->size=new Array2DSizeVarLenV(Ls);
            }
        }
        Array2DSizeDPS::Array2DSizeDPS(const Array2DSizeDPS&obj){
            this->SetNull();
            this->Assign(obj);
        }

        Array2DSizeDPS::~Array2DSizeDPS(){
            delete this->size;
        }

        void Array2DSizeDPS::Assign(const Array2DSizeDPS&obj){
            this->size=obj.size->clone();
        }

        Array2DSizeDPS& Array2DSizeDPS::operator =(const Array2DSizeDPS&obj){
            this->Assign(obj);
            return*this;
        }

        //void Array2DSizeDPS::Construct(){}

        void Array2DSizeDPS::SetNull(){
            this->size->SetNull();
        }

        //
        bool Array2DSizeDPS::IsRectangularType()const{
            return this->size->IsRectangularType();
        }
        bool Array2DSizeDPS::IsRectangularFact()const{
             return this->size->IsRectangularFact();
        }
        bool Array2DSizeDPS::IsVariableLengthType()const{
            return this->size->IsVariableLengthType();
        }
        //bool Array2DSizeDPS::IsVariableLenghFact()const{
        //    return this->size->IsVariableLenghFact();
        //}
        //
        //virtual void SetVaria();
        //virtual void SetRectangular(int L_ifVaria_min_m1_max_less=-2);
        //
        void Array2DSizeDPS::SetOne(bool isVaria){
            this->size->SetOne();
        }
        bool Array2DSizeDPS::ExtRowNBelongsTo(int extRowN)const{
            return this->size->ExtRowNBelongsTo(extRowN);
        }
        bool Array2DSizeDPS::PosBelongsTo(int extRowN, int ineRowN){
            return this->size->PosBelongsTo(extRowN, ineRowN);
        }
        void Array2DSizeDPS::Set(int QExtRows, int ExtRowsLength, bool PreserveType){
            if(this->size->IsRectangularType()==false && PreserveType==false){
                delete this->size;
                this->size=new Array2DSizeVarLenV(QExtRows, ExtRowsLength);
            }else{
                this->size->Set(QExtRows, ExtRowsLength);
            }
        }
        void Array2DSizeDPS::Set(int QExtRows,int*Ls, bool RectIfAllSame){
            int Lmin, Lmax, Lcur;
            if(Ls!=NULL){
                for(int i=1; i<=QExtRows; i++){
                    Lcur=Ls[i-1];
                    if(i==1 || (i>1 && Lcur<Lmin)){
                        Lmin=Lcur;
                    }
                    if(i==1 || (i>1 && Lcur>Lmax)){
                        Lmax=Lcur;
                    }
                }
                if(Lmin==Lmax && this->size->IsRectangularType()==false){
                    delete this->size;
                    this->size=new Array2DSizeRect(QExtRows, Lmax);
                }else{
                    this->size->Set(QExtRows, Ls);
                }
            }
        }

        //void Add(int L, int ifConstAndNotEqual_Ignore0AddExisting1SetAll2=2);
        //void Ins(int N, int L);
        void Array2DSizeDPS::AddExt(int L){
            this->size->AddExt(L);
        }

        void Array2DSizeDPS::InsExt(int N, int L){
            this->size->InsExt(N, L);
        }
        void Array2DSizeDPS::DelExt(int N){
            this->size->DelExt(N);
        }
        void Array2DSizeDPS::DelIneRow(){ this->size->DelIneRow(); }
        void Array2DSizeDPS::AddOrInsIneRow(){ this->size->AddOrInsIneRow(); }
        //virtual void AddOrInsTo(int L=0, int N_0ToAllPreserveType_GTT0AndMakeVaria=0);
        //virtual void DelFrom(int L, int N_0ToAllPreserveType_GTT0AndMakeVaria=0);
        void Array2DSizeDPS::AddOrInsTo(int N, int L){ this->size->AddOrInsTo(N, L); }
        void Array2DSizeDPS::DelFrom(int N, int L){ this->size->DelFrom(N, L); }
        int Array2DSizeDPS::GetQExtRows() const { return this->size->GetQExtRows(); }
        //bool GetIfIsVariaLength() const;
        //bool isVariaLength();
        //bool isRectangular() const;
        int Array2DSizeDPS::GetLength(int N)const { return this->size->GetLength(N); }
        int Array2DSizeDPS::GetMinLength()const { return this->size->GetMinLength(); }
        int Array2DSizeDPS::GetMaxLength()const { return this->size->GetMaxLength(); }
        void Array2DSizeDPS::SwapExtRows(int N1, int N2) { return this->size->SwapExtRows(N1, N2); }
        void Array2DSizeDPS::Reverse() { this->size->Reverse(); }
        //void Array2DSizeDPS::SetQExtRows(int N) { this->size->SetQExtRows(N); }
        //void Array2DSizeDPS::SetIneRowsLength(int L_MinMinus1MaxLess, bool preserveTypeIfVaria){
        //    int QExtRows=this->GetQExtRows();
        //    if(L_MinMinus1MaxLess==-1){
        //        L_MinMinus1MaxLess=this->size->GetMinLength();
        //    }else if(L_MinM1MaxLess<0){
        //        L_MinMinus1MaxLess=this->size->GetMaxLength();
        //    }
        //    if(preserveTypeIfVaria=false && this->size->IsRectangularType()==false){
        //        delete this->size;
        //        this->size=new Array2DSizeRect(QExtRows, L_MinMinus1MaxLess);
        //    }else{
        //        this->size->Set(QExtRows, L_MinMinus1MaxLess);
        //    }
        //}
        //void SetLength(int L, int N_0ToAllPreserveType_GTT0AndMakeVaria_Minus1ToAllMakeNotVaria=0);
        void Array2DSizeDPS::SetLength(int L, int N){ this->size->SetLength(L, N); }

        QString Array2DSizeDPS::GetAsQString(QString delim)const{
            return this->size->GetAsQString(delim);
        }
        std::string Array2DSizeDPS::GetAsStdString(std::string delim) const{
            return this->size->GetAsStdString(delim);
        }
        void Array2DSizeDPS::ShowToConsole() const{ this->size->ShowToConsole(); }
        //bool Array2DSizeDPS::isSame(const Array2DSizeDPS&obj)const;
        int Array2DSizeDPS::GetIneRowsLength(bool EvenIfVarButAllEqual)const{
            int L;
            if(this->size->IsRectangularType()==false && (EvenIfVarButAllEqual==false || this->size->IsRectangularFact()==false)){
                L=0;
            }else{
                L=this->size->GetLength(0);
            }
            return L;
        }
        //
        int Array2DSizeDPS::LengthFullOfIneRowN(int N) const{//from start to 1st stop
            return this->size->LengthFullOfIneRowN(N);
        }
        bool Array2DSizeDPS::IsForTranspose() const { return this->size->IsForTranspose(); }
        void Array2DSizeDPS::Transpose() { this->size->Transpose(); }
        //
        void Array2DSizeDPS::Set1ExtRowVaria(int L){
            std::vector<int>Ls;
            delete this->size;
            Ls.push_back(L);
            this->size=new Array2DSizeVarLenV(Ls);
        }

        void Array2DSizeDPS::SetVectorPotentiallyVaria(int Q){
            std::vector<int>Ls;
            delete this->size;
            for(int i=1; i<=Q; i++){
                Ls.push_back(1);
            }
            this->size=new Array2DSizeVarLenV(Ls);
        }
        void Array2DSizeDPS::SetVaria(int QExtRows, int ExtRowsLength){
            std::vector<int>Ls;
            delete this->size;
            for(int i=1; i<=QExtRows; i++){
                Ls.push_back(ExtRowsLength);
            }
            this->size=new Array2DSizeVarLenV(Ls);
        }
        void Array2DSizeDPS::ConvertToRectIfAllLengthesAreEqual(){
            int Q=this->GetQExtRows(), L=this->GetMaxLength();
            if(this->size->IsRectangularType()==false && this->size->IsRectangularFact()==true){
                delete this->size;
                this->size=new Array2DSizeRect(Q, L);
            }
        }

        void Array2DSizeDPS::ConvertToRect(int L_ifVaria_min_m1_max_less){
            int Q=this->GetQExtRows();
            if(L_ifVaria_min_m1_max_less==-1){
                L_ifVaria_min_m1_max_less=this->size->GetMinLength();
            }else if(L_ifVaria_min_m1_max_less<0){
                L_ifVaria_min_m1_max_less=this->size->GetMaxLength();
            }
            delete this->size;
            this->size=new Array2DSizeRect(Q, L_ifVaria_min_m1_max_less);
        }
        void Array2DSizeDPS::ConvertToVaria(){
            std::vector<int>Ls;
            int QExtRows=this->size->GetQExtRows(), ExtRowsLength=this->size->GetLength(0);
            if(this->size->IsRectangularType()==false){
                for(int i=1; i<=QExtRows; i++){
                    Ls.push_back(ExtRowsLength);
                }
                delete this->size;
                this->size=new Array2DSizeVarLenV(Ls);
            }
        }

//};

 */


//=======================================================================================================

// ------------------------------------------------------------------------------------------------------

int CalcElementNIfNegativeOf(int N, int Q){
    int y=N;
    if(N<0){
        y=N+Q;
    }
    return y;
}


//===================================================================================================

//===================================================================================================

/*
//class Matrix_P{
//    Array2D_PS<double>data;
//public:
    Matrix_P::Matrix_P(int QLines, int QColumns, double**X){
        double*DefaultDouble=new double;
        (*(DefaultDouble))=0;
        this->data.Set(QLines, QColumns, true, X, DefaultDouble);
        delete DefaultDouble;
        //wi checked et appended by Set
    }

    Matrix_P::Matrix_P(std::vector<std::vector<double>>X){
        this->data.Set(X);//X, true, ics'dflt
    }

    Matrix_P::Matrix_P(int Q, double*X, bool VectNotLine){
        this->data.Set(Q, VectNotLine, X);//,true
    }

    Matrix_P::Matrix_P(std::vector<double>X, bool VectNotLine){
        this->data.Set(X, VectNotLine);//true
    }

    Matrix_P::Matrix_P(const Matrix_P& obj){

    }

    Matrix_P::~Matrix_P(){}

    void Matrix_P::Assign(const Matrix_P& obj){
        this->data=obj.data;
    }

    Matrix_P& Matrix_P::operator = (const Matrix_P& obj){
        this->Assign(obj);
        return *this;
    }

    void Matrix_P::SetSize(int QLines, int QColumns){
        double*DefaultDouble=new double;
        (*(DefaultDouble))=0;
        this->data.SetSize(QLines, QColumns, true, DefaultDouble);
        delete DefaultDouble;
    }

    int Matrix_P::GetQLines(){
        return this->data.GetQExtRows();
    }

    int Matrix_P::GetQColumns(){
        return this->data.GetLength(1);
    }

    void Matrix_P::Set(int QLines, int QColumns, double**X){
        double*DefaultDouble=new double;
        (*(DefaultDouble))=0;
        this->data.Set(QLines, QColumns, true, X, DefaultDouble);
        delete DefaultDouble;
    }

    void Matrix_P::Set(std::vector<std::vector<double>>X){
        this->data.Set(X);//true
        //append it with stretching!
    }

    void Matrix_P::Set(int Q, double*X, bool VectNotLine){
        double*DefaultDouble=new double;
        (*(DefaultDouble))=0;
        this->data.Set(Q, VectNotLine, X);//rect=true, defaultVal n'ved qut
        delete DefaultDouble;
    }

    void Matrix_P::Set(std::vector<double>X, bool VectNotLine){
        this->data.Set(X, VectNotLine);//,rect=true)
    }

    std::vector<double> Matrix_P::GetLineN(int LineN){

    }

    std::vector<double> Matrix_P::GetColumnN(int LineN){

    }

    double Matrix_P::GetElement(int LineN, int ColumnN){

    }

    void Matrix_P::SetElement(double x, int LineN, int ColumnN){

    }

//};

//class Matrix_V{
//    Array2D_V<double>data;
//public:
    Matrix_V::Matrix_V(int QLines, int QColumns, double**X){

    }

    Matrix_V::Matrix_V(std::vector<std::vector<double>>X){

    }

    Matrix_V::Matrix_V(int Q, double*X, bool VectNotLine){

    }

    Matrix_V::Matrix_V(std::vector<double>X, bool VectNotLine){

    }

    Matrix_V::Matrix_V(const Matrix_V& obj){

    }

    Matrix_V::~Matrix_V(){}

    void Matrix_V::Assign(const Matrix_V& obj){

    }

    Matrix_V& Matrix_V::operator = (const Matrix_V& obj){

    }

    void Matrix_V::SetSize(int QLines, int QColumns){

    }

    int Matrix_V::GetQLines(){

    }

    int Matrix_V::GetQColumns(){

    }

    void Matrix_V::Set(int QLines, int QColumns, double**X){

    }

    void Matrix_V::Set(std::vector<std::vector<double>>X){

    }

    void Matrix_V::Set(int Q, double*X, bool VectNotLine){

    }

    void Matrix_V::Set(std::vector<double>X, bool VectNotLinee){

    }

    std::vector<double> Matrix_V::GetLineN(int LineN){

    }

    std::vector<double> Matrix_V::GetColumnN(int LineN){

    }

    double Matrix_V::GetElement(int LineN, int ColumnN){

    }

    void Matrix_V::SetElement(double x, int LineN, int ColumnN){

    }
    */

//};
