#ifndef MYARRAYLIB_H
#define MYARRAYLIB_H

#include"mylib.h"

//class MyArrayLib
//{
//public:
//    MyArrayLib();
//};

//vikts! Ei fns, mostly om ine rows, ne copy'd hin. Et SetExtRow mab fa based on preMarks, mab AddExtRow ma mab no. Ins Ext Rows juq ns'fa

//=== wi 1D arrays ======================================================================================================================================

template <typename T> void Array1DShowConsole(T*x, int Q, std::string delim="; ", bool LineNotCol=true){
    if(x!=NULL){
        for(int i=1; i<=Q-1; i++){
            std::cout<<x[i-1];
            if(LineNotCol){
                std::cout<<delim;
            }else{
                std::cout<<std::endl;
            }
        }
        std::cout<<x[Q-1];
        std::cout<<std::endl;
    }
}
template <typename T> void Array1DShowConsole(std::vector<T>x, std::string delim="; ", bool LineNotCol=true){
    int Q=x.size();
    if(x.size()!=0){
       for(int i=1; i<=Q-1; i++){
           std::cout<<x[i-1];
           if(LineNotCol){
               std::cout<<delim;
           }else{
               std::cout<<std::endl;
           }
       }
       std::cout<<x[Q-1];
       std::cout<<std::endl;
   }
}

template <typename T> std::vector<T> Array1DAssign_VbyP(int Qy, T*y=NULL,  T*DfltValParam=NULL){
    std::vector<T>vec;
    T DfltVal;
    if(DfltValParam!=NULL) DfltVal=(*(DfltValParam));
    if(Qy==0){
        if(y==NULL){
            vec.clear();
        }
    }else{
        if(y==NULL){
            for(int i=1; i<=Qy; i++){
                vec.push_back(DfltVal);
            }
        }else{
            for(int i=1; i<=Qy; i++){
                vec.push_back(y[i-1]);
            }
        }
    }
    return vec;
}
template <typename T> void Array1DAssign(std::vector<T>&vec, T y[], int L=0){
    vec.clear();
    int Q=vec.size();
    if(L==0){
        vec.assign(y, y+sizeof(y)/sizeof(y[0]));
    }else{
        for(int i=1; i<=L; i++){
            vec.push_back(y[i-1]);
        }
    }
}

template <typename T> void  Array1DAssign_PbyP(T*&x, int Qx=0, T*y=NULL, int Qy=0,  T*DfltValParam=NULL){
    T DfltVal;
    if(Qy==0){
        if(y!=NULL){
            Qy=Qx;
        }
    }
    if(Qx!=Qy && Qx>0){
        if(x!=NULL){
            delete[]x;
            x=NULL;
        }
    }
    if(x==NULL && Qy>0){
        x=new T[Qy];
    }
    if(Qy>0){
        if(y!=NULL){
            for(int i=1; i<=Qy; i++){
                x[i-1]=y[i-1];
            }
        }else{
            if(DfltValParam!=NULL){
                DfltVal=(*(DfltValParam));
            }
            for(int i=1; i<=Qy; i++){
                x[i-1]=DfltVal;
            }
        }
    }
}
template <typename T> void  Array1DAssign(T*&x, int Qx, std::vector<T>y, T*DfltValParam=NULL){
    Array1DAssign(x, Qx, y.data(), y.size(), DfltValParam);
}
template <typename T> void  Array1DAssign(std::vector<T>&x, std::vector<T>y){
    x=y;
}

template <typename T> void  Array1DAssign(T*&x, T*y=NULL, int Q=0, T*DfltValParam=NULL){
    T DetfaulVal;
    if(DfltValParam!=NULL){
        DetfaulVal=(*(DfltValParam));
    }
    if(x!=NULL){
        delete[]x;
        x=NULL;
    }
    if(y!=NULL && Q>0){
        x=new T[Q];
        for(int i=1; i<=Q; i++){
            x[i-1]=y[i-1];
        }
    }else if(y==NULL && Q>0){
        for(int i=1; i<=Q; i++){
            x[i-1]=DetfaulVal;
        }
    }
}
template <typename T> void  Array1DAssign(T*&x, T DfltValParam, int Q=0, T*y=NULL){
    if(x!=NULL){
        delete[]x;
        x=NULL;
    }
    if(y!=NULL && Q>0){
        x=new T[Q];
        for(int i=1; i<=Q; i++){
            x[i-1]=y[i-1];
        }
    }else if(y==NULL && Q>0){
        for(int i=1; i<=Q; i++){
            x[i-1]=DfltValParam;
        }
    }
}
template <typename T> void  Array1DAssign(T*&x, std::vector<T>y){
    T DfltRndmVal;
    Array1DAssign(x, DfltRndmVal, y.size(), y.data());
}

template <typename T> void Array1DAssignRepeatingSingleVal(T val, T*&X, int Q=1){
    if(X!=NULL){
        delete[]X;
    }
    X=new T[Q];
    for(int i=1; i<=Q; i++){
        X[i-1]=val;
    }
}
template <typename T> void Array1DAssignRepeatingSingleVal(T val, std::vector<T>&X, int Q=1){
    X.clear();
    for(int i=1; i<=Q; i++){
        X.push_back(val);
    }
}
template <typename T> std::vector<T> Array1DAssignRepeatingSingleVal(T val, int Q=1){
    std::vector<T>X;
    for(int i=1; i<=Q; i++){
        X.push_back(val);
    }
    return X;
}

template <typename T>std::vector<T> Array1DGetSubArray(T*x, int wholeRowL, int FromN=1, int ToN=0, T*DfltValParam=NULL){
    std::vector<T>y;
    int LastN;
    T DfltVal;
    if(DfltValParam!=NULL){
        DfltVal=(*(DfltValParam));
    }
    if(ToN==0){
        ToN=wholeRowL;
    }
    if(FromN<1){
        FromN+=wholeRowL;
    }
    if(ToN<1){
        ToN+=wholeRowL;
    }
    if(FromN>0 || ToN>0){
        if(FromN<0){
            for(int i=FromN; i<1; i++){
                y.push_back(DfltVal);
            }
            if(ToN<=wholeRowL){
                for(int i=1; i<=ToN; i++){
                    y.push_back(x[i-1]);
                }
            }else{
                for(int i=1; i<=wholeRowL; i++){
                    y.push_back(x[i-1]);
                }
                for(int i=wholeRowL+1; i<=ToN; i++){
                    y.push_back(DfltVal);
                }
            }
        }else if(FromN<ToN && FromN<=wholeRowL){
            if(ToN<=wholeRowL){
                for(int i=FromN; i<=ToN; i++){
                    y.push_back(x[i-1]);
                }
            }else{
                for(int i=FromN; i<=wholeRowL; i++){
                    y.push_back(x[i-1]);
                }
                for(int i=wholeRowL+1; i<=ToN; i++){
                    y.push_back(DfltVal);
                }
            }
        }else if(ToN<FromN && ToN<=wholeRowL && FromN>=1){
            for(int i=FromN; i>=wholeRowL+1; i--){
                y.push_back(DfltVal);
            }
            if(ToN>=1){
                for(int i=wholeRowL; i>=ToN; i--){
                     y.push_back(x[i-1]);
                }
            }else{
                for(int i=wholeRowL; i>=1; i--){
                     y.push_back(x[i-1]);
                }
                for(int i=0; i>=ToN; i--){
                   y.push_back(DfltVal);
                }
            }
        }
    }
    return y;
}
template <typename T>std::vector<T> Array1DGetSubArray(std::vector<T>x,  int FromN=1, int ToN=0, T*DfltValParam=NULL){
   return  Array1DGetSubArray(x.data(), x.size(), FromN, ToN, DfltValParam);

}

template <typename T>std::vector<T> Array1DGetSubArray_byNs(T*x, int LOld, int N1=1, int N2=-1, T*DfltValParam=NULL){
    std::vector<T>y;
    T DfltVal;
    int Nmin, Nmax;
    if(DfltValParam!=NULL){
        DfltVal=(*(DfltValParam));
    }
    if(N1<0){
        N1+=(LOld+1);
    }
    if(N2<0){
        N2+=(LOld+1);
    }
    if(N1<=N2){
        if(N1<0){
            for(int i=1; i<=(-N1-1); i++){
                y.push_back(DfltVal);
            }
            N1=1;
        }
        Nmin= N2<=LOld ? N2: LOld;
        for(int i=1; i<=Nmin; i++){
            y.push_back(x[i-1]);
        }
        for(int i=LOld+1; i<=N2; i++){
            y.push_back(DfltVal);
        }
    }else{
        if(N1>LOld){
           for(int i=N1; i>=LOld+1; i--){
               y.push_back(DfltVal);
           }
           N1=LOld;
        }
        Nmax = N2>=1 ? N2: 1;
        for(int i=N1; i>=Nmax; i--){
             y.push_back(x[i-1]);
        }
        if(N2<0){
            for(int i=1; i>=(-N2); i--){
                y.push_back(DfltVal);
            }
        }
    }
    return y;
}
template <typename T>std::vector<T> Array1DGetSubArray_byNs(std::vector<T>x,  int N1=1, int N2=-1, T*DfltValParam=NULL){
    return  Array1DGetSubArray_byNs(x.data(), x.size(), N1, N2, DfltValParam);
}

template <typename T>std::vector<T> Array1DGetSubArray_byLs(T*x, int Lold, int Lnew=0, int FromN=1, int FirstDefaultValues=0, T*DfltValParam=NULL){
    std::vector<T>y;
    T DfltVal;
    bool contin;
    int whereCurN, whatCurN;
    if(DfltValParam!=NULL){
        DfltVal=(*(DfltValParam));
    }
    if(Lnew==0){
        Lnew=Lold-(FromN-1)+FirstDefaultValues;
    }
    for(int whereCurN=1; whereCurN<=FirstDefaultValues; whereCurN++){
        y.push_back(DfltVal);
    }
    whereCurN=FirstDefaultValues;
    whatCurN=FromN-1;
    contin=(Lold>0);
    //contin=(whereN>=1 && whereN<Lold && whatN<Lnew);
    if(x!=NULL){
        while(contin){
            whatCurN++;
            whereCurN++;
            y.push_back(x[whatCurN-1]);
            if(whereCurN==Lnew){
                contin=false;
            }
            if(whatCurN==Lold){
                contin=false;
            }
        }
    }else{
        while(contin){
            whatCurN++;
            whereCurN++;
            y.push_back(DfltVal);
            //contin=(whereN<Lold && whatN<Lnew);
            if(whereCurN==Lnew){
                contin=false;
            }
            if(whatCurN==Lold){
                contin=false;
            }
        }
    }
    whereCurN=Lold-(FromN-1)+FirstDefaultValues;
    if(whereCurN<Lnew){
        for(int i=whereCurN+1; i<=Lnew; i++){
            y.push_back(DfltVal);
        }
    }
    return y;
}

template <typename T>std::vector<T> Array1DGetSubArray_byLs_PreMark(T*x, int Lold, int Lreq=0, int whatFromN=1, int FirstDefaultValues=0, T*DfltValParam=NULL){
    std::vector<T>y;
    T DefaultValue, CurrentValue;
    int preN1_=0, preN2_=0, ownN1_=0, ownN2=0, ForOwnN1_=0, ForOwnN2=0, postN1=0, postN2_=0, Lnew=0;
    if(DfltValParam!=NULL)DefaultValue=*DfltValParam;
    //void Calc_Array1DSubArray_byLs_MarkupNs_vFmls(int&LreqCalcd, int&preN1_, int&preN2_, int&ownN1_, int&ownN2, int&ForOwnN1_, int&ForOwnN2, int&postN1, int&postN2_, int Lini, int LreqGiven=0, int whatFromN=1, int QDefaultValuesBefore=0);
    Calc_Array1DSubArray_byLs_MarkupNs_vFmls(Lnew, preN1_, preN2_, ownN1_, ownN2, ForOwnN1_, ForOwnN2, postN1, postN2_, Lold, Lreq, whatFromN, FirstDefaultValues);
    //wa ef corr: Calc_Array1DSubArray_byLs_MarkupNs(preN1_, preN2_, ownN1_, ownN2, ForOwnN1_, ForOwnN2, postN1, postN2_, Lold, Lnew, FromN, FirstDefaultValues);
    for(int i=1; i<=preN2_; i++){
        y.push_back(DefaultValue);
    }
    if(ForOwnN1_>0){
        for(int i=ownN1_; i<=ownN2; i++){
            CurrentValue=x[i-1];
            y.push_back(CurrentValue);
        }
    }
    if(postN1>0){
        for(int i=postN1; i<=postN2_; i++){
            //CurrentValue=x[i-1];
            y.push_back(DefaultValue);
        }
    }
    return y;
}
template <typename T>std::vector<T> Array1DGetSubArray_byLs_PreMark(std::vector<T>x, int Lreq=0, int whatFromN=1, int FirstDefaultValues=0, T*DfltValParam=NULL){
    std::vector<T>y;
    int Lold = x.size();
    T DefaultValue, CurrentValue;
    int preN1_=0, preN2_=0, ownN1_=0, ownN2=0, ForOwnN1_=0, ForOwnN2=0, postN1=0, postN2_=0, Lnew=0;
    if(DfltValParam!=NULL)DefaultValue=*DfltValParam;
    //wa ef corr: Calc_Array1DSubArray_byLs_MarkupNs(preN1_, preN2_, ownN1_, ownN2, ForOwnN1_, ForOwnN2, postN1, postN2_, Lold, Lnew, FromN, FirstDefaultValues);
    Calc_Array1DSubArray_byLs_MarkupNs_vFmls(Lnew, preN1_, preN2_, ownN1_, ownN2, ForOwnN1_, ForOwnN2, postN1, postN2_, Lold, Lreq, whatFromN, FirstDefaultValues);//taks in ef-fn, ms'arb hin
    for(int i=1; i<=preN2_; i++){
        y.push_back(DefaultValue);
    }
    if(ForOwnN1_>0){
        for(int i=ownN1_; i<=ownN2; i++){
            CurrentValue=x[i-1];
            y.push_back(CurrentValue);
        }
    }
    if(postN1>0){
        for(int i=postN1; i<=postN2_; i++){
            //CurrentValue=x[i-1];
            y.push_back(DefaultValue);
        }
    }
    return y;
}


template <typename T>std::vector<T> Array1DGetSubArray_byLs(std::vector<T>x, int Lnew=0, int FromN=1, int FirstDefaultValues=0, T*DfltValParam=NULL){
    //std::vector<T> Array1DGetSubArray_byLs_PreMark(T*x, int Lold, int Lreq=0, int FromN=1, int FirstDefaultValues=0, T*DfltValParam=NULL){
    return Array1DGetSubArray_byLs(             x.data(), x.size(),       Lnew,       FromN,       FirstDefaultValues, DfltValParam);
}




template <typename T>std::vector<T> Array1DGetSubArray_byLs_PreMark_v1(T*x, int Lold, T*DefaultRow, int Lreq=0, int whatFromN=1, int FirstDefaultValues=0, int dfltRowL_ifNegativeThanETLold=-1){
    //woi abls T*DefaultRow=NULL, ma so S'ak fn co no dfltVal et no shift
    //if DefaultRow==NULL so used random default val et dfltRowL_ifNegativeThanETLold=0, ma so - qut utf hic fn? Es fns co dflt 1 val et co no dflts
    //if DefaultRow s'satq long, dfltRowL_ifNegativeThanETLold s'am fin et by dflt <0
    //
    std::vector<T>y;
    T DefaultRandomValue, CurrentValue;
    if(DefaultRow==NULL){
        dfltRowL_ifNegativeThanETLold=0;
    }else if(dfltRowL_ifNegativeThanETLold<0){
        dfltRowL_ifNegativeThanETLold=Lold;
    }
    int preN1_=0, preN2_=0, ownN1_=0, ownN2=0, ForOwnN1_=0, ForOwnN2=0, postN1=0, postN2_=0, Lnew=0;
    //if(DfltValParam!=NULL)DefaultValue=*DfltValParam;
    //void Calc_Array1DSubArray_byLs_MarkupNs_vFmls(int&LreqCalcd, int&preN1_, int&preN2_, int&ownN1_, int&ownN2, int&ForOwnN1_, int&ForOwnN2, int&postN1, int&postN2_, int Lini, int LreqGiven=0, int whatFromN=1, int QDefaultValuesBefore=0);
    Calc_Array1DSubArray_byLs_MarkupNs_vFmls(Lnew, preN1_, preN2_, ownN1_, ownN2, ForOwnN1_, ForOwnN2, postN1, postN2_, Lold, Lreq, whatFromN, FirstDefaultValues);
    if(dfltRowL_ifNegativeThanETLold>=preN2_){
        for(int i=1; i<=preN2_; i++){
            CurrentValue=DefaultRow[i-1];
            y.push_back(CurrentValue);
        }
    }else{//if >0 && < preN2_ for cycle wi arbf S auto
        for(int i=1; i<=dfltRowL_ifNegativeThanETLold; i++){
            CurrentValue=DefaultRow[i-1];
            y.push_back(CurrentValue);
        }
        CurrentValue=DefaultRandomValue;
        for(int i=dfltRowL_ifNegativeThanETLold+1; i<=preN2_; i++){
            y.push_back(CurrentValue);
        }
    }
    if(ForOwnN1_>0){
        for(int i=ownN1_; i<=ownN2; i++){
            CurrentValue=x[i-1];
            y.push_back(CurrentValue);
        }
    }
    if(postN1>0){
        if(dfltRowL_ifNegativeThanETLold>=postN2_){
            for(int i=postN1; i<=postN2_; i++){
                CurrentValue=DefaultRow[i-1];
                y.push_back(CurrentValue);
            }
        }else if(dfltRowL_ifNegativeThanETLold>=postN1&& dfltRowL_ifNegativeThanETLold<=postN2_){
            for(int i=postN1; i<=dfltRowL_ifNegativeThanETLold; i++){
                CurrentValue=DefaultRow[i-1];
                y.push_back(CurrentValue);
            }
            CurrentValue=DefaultRandomValue;
            for(int i=dfltRowL_ifNegativeThanETLold+1; i<=postN2_; i++){
                y.push_back(CurrentValue);
            }
        }else{//dfltRowL_ifNegativeThanETLold<postN1
            CurrentValue=DefaultRandomValue;
            for(int i=postN1; i<=postN2_; i++){
                y.push_back(CurrentValue);
            }
        }
    }
    return y;
}

template <typename T>std::vector<T> Array1DGetSubArray_byLs_PreMark_v1(T*x, int Lold, T DefaultValue, int Lreq=0, int whatFromN=1, int FirstDefaultValues=0){
    std::vector<T>dfltRow;
    int LDfltMax=Lold+Lreq;
    if(whatFromN>0)LDfltMax+=whatFromN;
    if(FirstDefaultValues>0)LDfltMax+=FirstDefaultValues;
    for(int i=1; i<=LDfltMax; i++) dfltRow.push_back(DefaultValue);
    //std::vector<T> Array1DGetSubArray_byLs_PreMark_v1(T*x, int Lold, T*DefaultRow, int Lreq=0, int whatFromN=1, int FirstDefaultValues=0, int dfltRowL_ifNegativeThanETLold=-1){
    return           Array1DGetSubArray_byLs_PreMark_v1(  x,     Lold, dfltRow.data(),   Lreq,       whatFromN,       FirstDefaultValues=0, dfltRow.size());
}

template <typename T>std::vector<T> Array1DGetSubArray_byLs_PreMark_v1(T*x, int Lold, int Lreq=0, int whatFromN=1, int FirstDefaultValues=0){
    T DefaultRandomValue;
    //template <typename T>std::vector<T> Array1DGetSubArray_byLs_PreMark_v1(T*x, int Lold, T DefaultValue, int Lreq=0, int whatFromN=1, int FirstDefaultValues=0)
    return Array1DGetSubArray_byLs_PreMark_v1(                                 x,     Lold, DefaultRandomValue, Lreq,       whatFromN,       FirstDefaultValues);
    //also possible:
    //template <typename T>std::vector<T> Array1DGetSubArray_byLs_PreMark_v1(T*x, int Lold, T*DefaultRow, int Lreq=0, int whatFromN=1, int FirstDefaultValues=0, int dfltRowL_ifNegativeThanETLold=-1){
    //return                                Array1DGetSubArray_byLs_PreMark_v1(x,       Lold,    NULL,          Lreq,          whatFromN,      FirstDefaultValues,      0);
}






template <typename T>void Array1DSetSize(T*&x, int QOld, int QNew,  T*DfltValExt=NULL, int PreserveVals=true){
    T*y=NULL;
    T DefaultVal;
    if(DfltValExt!=NULL){
        DefaultVal=(*(DfltValExt));
    }
    int mn= QOld <= QNew ? QOld : QNew;
    if(x!=NULL && QOld!=QNew){
        y=new T[QNew];
        if(PreserveVals){
            for(int i=1; i<=mn; i++){
                y[i-1]=x[i-1];
            }
            for(int i=mn+1; i<=QNew; i++){
                y[i-1]=DefaultVal;
            }
        }
        delete[]x;
        x=y;
    }
}

template <typename T>void Array1DSetSize(std::vector<T>&x, int QNew,  T*DfltValExt=NULL, int PreserveVals=true){
    T*y=NULL;
    T DefaultVal;
    int QOld=x.size();
    if(DfltValExt!=NULL){
        DefaultVal=(*(DfltValExt));
    }
    int mn= QOld <= QNew ? QOld : QNew;
    if(x!=NULL && QOld!=QNew){
        y=new T[QNew];
        if(PreserveVals){
            for(int i=1; i<=mn; i++){
                y[i-1]=x[i-1];
            }
            for(int i=mn+1; i<=QNew; i++){
                y[i-1]=DefaultVal;
            }
        }
        x.clear();
        for(int i=1; i<=QNew; i++){
            x.push_back(y[i-1]);
        }
        delete[]y;
        y=NULL;
    }
}

template <typename T>void Array1DSetSize_v1(T*&x, int LOld, int LNew, T*dfltRow, int dfltRowLen_ifNegativeETLOld){
    std::vector<T>rowTmp=Array1DGetSubArray_byLs_PreMark_v1(x, LOld, dfltRow, LNew, 1, 0, dfltRowLen_ifNegativeETLOld);
    T val;
    if(x!=NULL){
        delete[]x;
        x=NULL;
    }
    if(LNew>0){
        x=new T[LNew];
        for(int i=1; i<=LNew; i++){
            val=rowTmp[i-1];
            x[i-1]=val;
            //x[i-1]=rowTmp[i-1];
        }
    }
}
template <typename T>void Array1DSetSize_v1(T*&x, int LOld, int LNew, T dfltVal){
    //std::vector<T>rowTmp=Array1DGetSubArray_byLs_PreMark_v1(x, LOld, dfltVal, LNew, 1, 0);
    std::vector<T>rowTmp;
    for(int i=1; i<=LNew; i++)rowTmp.push_back(dfltVal);
    Array1DSetSize_v1(x, LOld, LNew, rowTmp.data(), rowTmp.size());
}
template <typename T>void Array1DSetSize_v1(T*&x, int LOld, int LNew){
    T DfltRndmVal;
    Array1DSetSize_v1(x, LOld, LNew, DfltRndmVal);
}

template <typename T>void Array1DReSize(T*&x, int&Q, int dQ, T*DfltVal=NULL, int PreserveVals=true){
    int QOld=Q;
    Q+=dQ;
    Array1DSetSize(x, QOld, Q, DfltVal, PreserveVals);
}

template <typename T>void Array1DReSize(std::vector<T>&x, int dQ, T*DfltVal=NULL, int PreserveVals=true){
    int Q=x.size();
    int QOld=Q;
    Q+=dQ;
    Array1DSetSize(x, QOld, Q, DfltVal, PreserveVals);
}

template <typename T>void Array1DAddElement(T*&x, int&Q, T z){//new considers vrn x=NULL
    T*y=NULL;
    T buf;
    if(x==NULL) Q=0;
    y=new T[Q+1];
    if(x!=NULL){
        for(int i=1; i<=Q; i++){
            //y[i-1]=x[i-1];
            buf=x[i-1];
            y[i-1]=buf;
        }
    }
    y[Q+1-1]=z;
    if(x!=NULL){
        delete[]x;
    }
    x=y;
    Q++;
}
template <typename T>void Array1DAddElement(std::vector<T>&x, T z){
    x.push_back(z);
}


template <typename T>void Array1DInsElementToN(T*&x, int&Q, int N, T z, TValsShowHide*vsh=NULL/*, int PreserveVals=true*/){
    T*y;
   // int OldN, NewN;
    std::string msg="VectorInsElementToN starts working";
    writeln(vsh,msg);//"VectorInsElementToN starts working");
    if(N<1) N=Q+N;//N=CalcElementNIfNegativeOf(int N, int Q);
    if(N>=1 && N<=Q){
        msg="N is within range of arrays size";
        writeln(vsh,msg);//"N is within range of arrays size");
        y=new T[Q+1];
       // if(PreserveVals){
        for(int i=1; i<=Q+1; i++){
            writeln(vsh," i="+IntToStr(i)+" N="+IntToStr(N)+" Q="+IntToStr(Q));
            if(i<N){
                //oldN=i;
                //newN=i;
                y[i-1]=x[i-1];
                writeln(vsh,IntToStr(i)+" -> "+IntToStr(i));
            } else if(i>N){
                //oldN=i;
                //newN=i-1;
                y[i-1]=x[i-1-1];
                writeln(vsh,IntToStr(i)+" -> "+IntToStr(i-1));
            }else{
                writeln(vsh,IntToStr(i)+" =N! ");
            }
        }
        y[N-1]=z;
        delete[]x;
        x=y;
        Q++;
    }
    msg="VectorInsElementToN finishes working";
    writeln(vsh,msg);//"VectorInsElementToN finishes working");
}

template <typename T>void Array1DInsElementToN(std::vector<T>&x, int N, T z){
    x.insert(x.begin()+N-1, z);
}

template <typename T>void Array1DDelElementFromN(T*&x, int&Q, int N){
    T*y=NULL;
    if(N<1) N=Q-N;
    if(N>=1 && N<=Q){
        y=new T[Q-1];
        for(int i=1; i<=Q; i++){
            if(i<N)y[i-1]=x[i-1];
            else if(i>N) y[i-1-1]=x[i-1];
        }
        delete[]x;
        x=y;
        Q--;
    }
}

template <typename T>void Array1DDelElementFromN(std::vector<T>&x, int N){
    int Q=x.size();
    if(N>=1 && N<=Q){
        x.erase(x.begin()+N-1);
    }
}

//
template <typename T>void Array1DSwapElements(T*&x, int Q, int N1, int N2){
    T buf;
    if(N1<1) N1=Q+N1; //if(N1<1) N1=Q+N1;//N=CalcElementNIfNegativeOf(int N1, int Q);
    if(N2<1) N2=Q+N2; //if(N2<1) N2=Q+N2;//N=CalcElementNIfNegativeOf(int N2, int Q);
    if(N1>=1 && N1<=Q && N2>=1 && N2<=Q){
        buf=x[N1-1];
        x[N1-1]=x[N2-1];
        x[N2-1]=buf;
    }
}
template <typename T>void Array1DSwapElements(std::vector<T>&x, int N1, int N2){
    T buf;
    int Q=x.size();
    if(N1<1) N1=Q+N1; //if(N1<1) N1=Q+N1;//N=CalcElementNIfNegativeOf(int N1, int Q);
    if(N2<1) N2=Q+N2; //if(N2<1) N2=Q+N2;//N=CalcElementNIfNegativeOf(int N2, int Q);
    if(N1>=1 && N1<=Q && N2>=1 && N2<=Q){
        buf=x[N1-1];
        x[N1-1]=x[N2-1];
        x[N2-1]=buf;
    }
}
template <typename T>void Array1DReverse(T*&x, int Q){
    int N, N1, N2;
    if(Q%2==0){
        N=Q/2;
    }else{
        N=(Q-1)/2;
    }
    for(int i=1; i<=N; i++){
        N1=i;
        N2=Q-i+1;
        Array1DSwapElements(x, Q, N1, N2);
    }
}

template <typename T>void Array1DReverse(std::vector<T>&x){
    int N, N1, N2;
    int Q=x.size();
    if(Q%2==0){
        N=Q/2;
    }else{
        N=(Q-1)/2;
    }
    for(int i=1; i<=N; i++){
        N1=i;
        N2=Q-i+1;
        Array1DSwapElements(x, N1, N2);
    }
}

template <typename T> std::vector<int>Array1DSeek(T*data, T val, int curL, bool(*Equals)(const T& x1, const T& x2, int params)=EqualSimple, int FromN=1, int ToN=0){
    std::vector<int> Ns;
    int Q=curL;
    if(Q>0){
        if(ToN<1){
            ToN=Q;
            if(FromN>ToN){
                FromN=ToN;
            }
            if(FromN<1){
                FromN=1;
            }
            for(int i=FromN; i<=ToN; i++){
                 //if(data[i-1]==val){
                 if(Equals(data[i-1],val)){
                    Ns.push_back(i);
                }
            }
        }
    }
    return Ns;
}//fn
template <typename T> std::vector<int>Array1DSeek(const std::vector<T>&data, T val, bool(*Equals)(const T& x1, const T& x2, int params)=EqualSimple, int FromN=1, int ToN=0){
    std::vector<int> Ns;
    int curL=data.size();
    int Q=curL;
    if(Q>0){
        if(ToN<1){
            ToN=Q;
            if(FromN>ToN){
                FromN=ToN;
            }
            if(FromN<1){
                FromN=1;
            }
            for(int i=FromN; i<=ToN; i++){
                 //if(data[i-1]==val){
                 if(Equals(data[i-1],val)){
                    Ns.push_back(i);
                }
            }
        }
    }
    return Ns;
}//fn
template <typename T> std::vector<int>Array1DSeekValSimply(T*data, T val, int curL, int FromN=1, int ToN=0, int param=0){
    std::vector<int> Ns;
    int Q=curL;
    if(ToN==0 or ToN>Q) ToN=Q;
    if(FromN<1){
        //FromN= CalcElementNIfNegativeOf(FromN, curL);
        FromN=Q+FromN+1;
    }
    //if(FromN<=ToN){
    if(FromN<=ToN && FromN>=1){
        for(int i=FromN; i<=ToN; i++){
            if(data[i-1]==val){
                Ns.push_back(i);
            }
        }
    }
    //else
    //{
    //    for(int i=FromN; i>=ToN; i--){
    //        if(data[i-1]==val){
    //            Ns.push_back(i);
    //        }
    //    }
    //}
    return Ns;
}//fn
template <typename T> std::vector<int>Array1DSeekValSimply(std::vector<T>data, T val,  int FromN=1, int ToN=0, int param=0){
    std::vector<int> Ns;
    int Q=data.size();
    if(ToN==0 or ToN>Q) ToN=Q;
    if(FromN<=ToN){
        for(int i=FromN; i<=ToN; i++){
            if(data[i-1]==val){
                Ns.push_back(i);
            }
        }
    }else{
        for(int i=FromN; i>=ToN; i--){
            if(data[i-1]==val){
                Ns.push_back(i);
            }
        }
    }
    return Ns;
}//fn
template <typename T> bool Array1DSuccIsAtPos(T*data, T*succ, int L, int PosN, int succL, bool(*Equals)(const T& x1, const T& x2, int param)=EqualSimple, TValsShowHide*vsh=NULL){
    bool verdict;
    int ChkN;
    writeln(vsh,"SuccIsAtPos starts working");
    writeln(vsh,"Search: in Length: "+IntToStr(L)+" at PosN="+IntToStr(PosN)+" Length to seek: "+IntToStr(succL));
    verdict=(PosN+succL-1<=L); //check if it can contain so long part
    if(verdict) writeln(vsh,"check passed/ It can contain so long part");
    else  writeln(vsh,"check done/ It can NOT contain so long part");
    if(verdict){
        for(int i=1; i<=succL; i++){
            ChkN=PosN+i-1;
            //if(!(succ[i-1]==data[ChkN-1])){
            //if(Equals(succ[i-1], data[ChkN-1], 0)==false){
            //    verdict=false;
            //}
            if(Equals(succ[i-1], data[ChkN-1], 0)==false){
                verdict=false;
            }
        }
    }
    writeln(vsh,"SuccIsAtPos finishes working");
    return verdict;
}
template <typename T> bool Array1DSuccIsAtPos(const std::vector<T>&data, const std::vector<T>&succ, int PosN, bool(*Equals)(const T& x1, const T& x2, int param)=EqualSimple, TValsShowHide*vsh=NULL){
    int L=data.size();
    int succL=succ.size();
    bool verdict;
    int ChkN;
    writeln(vsh,"SuccIsAtPos starts working");
    writeln(vsh,"Search: in Length: "+IntToStr(L)+" at PosN="+IntToStr(PosN)+" Length to seek: "+IntToStr(succL));
    verdict=(PosN+succL-1<=L); //check if it can contain so long part
    if(verdict) writeln(vsh,"check passed/ It can contain so long part");
    else  writeln(vsh,"check done/ It can NOT contain so long part");
    if(verdict){
        for(int i=1; i<=succL; i++){
            ChkN=PosN+i-1;
            //if(!(succ[i-1]==data[ChkN-1])){
            //if(Equals(succ[i-1], data[ChkN-1], 0)==false){
            //    verdict=false;
            //}
            if(Equals(succ[i-1], data[ChkN-1], 0)==false){
                verdict=false;
            }
        }
    }
    writeln(vsh,"SuccIsAtPos finishes working");
    return verdict;
}
template <typename T> bool Array1DSuccIsAtPosSimply(T*data, T*succ, int L, int succL, int PosN, int param=0, TValsShowHide*vsh=NULL){
    bool verdict=true;
    int j;
    for(int i=1; i<=succL; i++){
        j=PosN+i-1;
        if(data[j-1]!=succ[i-1]){
            verdict=false;
        }
    }
    return verdict;
}
//new fn:
template <typename T> bool Array1DSuccIsAtPosSimply_v1(T*where, T*what, int whereL, int whatL, int posN, bool Reversed=false, int param=0, TValsShowHide*vsh=NULL){
    bool verdict=true;
    int whereCurN, whatCurN;
    if(!Reversed){
        if((posN+whatL-1)<=whereL){
            for(int i=1; i<=whatL; i++){
                whatCurN=i;
                whereCurN=posN+i-1;
                if(where[whereCurN-1]!=what[whatCurN-1]){
                    verdict=false;
                }
            }
        }
    }else{//Reversed
        if((posN-whatL+1)>=1){
            for(int i=1; i<=whatL; i++){
                whatCurN=i;
                whereCurN=posN-i+1;
                if(where[whereCurN-1]!=what[whatCurN-1]){
                    verdict=false;
                }
            }
        }
    }
    return verdict;
}
template <typename T> bool Array1DSuccIsAtPosSimply(std::vector<T>data, std::vector<T>succ, int PosN, int param=0, TValsShowHide*vsh=NULL){
    bool verdict=true;
    int j, L=data.size(), succL=succ.size();
    for(int i=1; i<=succL; i++){
        j=PosN+i-1;
        if(data[j-1]!=succ[i-1]){
            verdict=false;
        }
    }
    return verdict;
}
template <typename T> std::vector<int>Array1DSeekSubArraySimply(T*where, T*what, int whereL, int whatL, int FromN=1, int ToN=0, int param=0, TValsShowHide*vsh=NULL){
    bool found;
    std::vector<int>Ns;
    if(FromN<0){
        FromN=whereL+FromN+1;
    }
    if(ToN<0){
        ToN=whereL+ToN+1;
    }else if(ToN==0){
        ToN=whereL;
    }
    if(abs(ToN-FromN)+1>=whatL){

    }else{

    }

    if(FromN>=1 && FromN<=whereL && ToN>=1 && ToN<=whereL && where!=NULL && what!=NULL && abs(ToN-FromN)+1>=whatL){
        if(FromN<=ToN){
            //if(FromN+whatL-1<=ToN){
                //for(int N=FromN; N<=ToN-whatL+1; N++){
                for(int N=FromN; N<=ToN; N++){
                    // bool Array1DSuccIsAtPosSimply(T*data, T*succ, int L, int succL, int PosN, int param=0, TValsShowHide*vsh=NULL)
                    found = Array1DSuccIsAtPosSimply(where, what, whereL, whatL, N, param, vsh);
                    if(found){
                        Ns.push_back(N);
                    }
                }
            //}
        }else{
            //for(int N=FromN; N>=ToN+whatL-1; N--){
            for(int N=FromN; N>=ToN; N--){
                // bool Array1DSuccIsAtPosSimply(T*data, T*succ, int L, int succL, int PosN, int param=0, TValsShowHide*vsh=NULL)
                found = Array1DSuccIsAtPosSimply(where, what, whereL, whatL, N, param, vsh);
                if(found){
                    Ns.push_back(N);
                }
            }
        }
    }
    return Ns;
}

template <typename T> std::vector<int>Array1DSeek(T*data, T*succ, int succL, int curL, bool (*Equals)(const T& x1, const T& x2, int params)=EqualSimple, int FromN=1, int ToN=0, TValsShowHide*vsh=NULL){
    std::vector<int> Ns;
    int Q=succL-curL+1;
    TValsShowHide vsh1;
    if(vsh!=NULL)vsh1=*vsh;
    vsh1.EnableWriting();
    writeln(&vsh1, "Seek(1D in 1D) starts working");
    if(Q>0){
        if(ToN<1){
            ToN=Q;
        }
       if(FromN>ToN){
           FromN=ToN;
       }
       if(FromN<1){
           FromN=1;
       }
       writeln(&vsh1,"Length to seek in: "+IntToStr(Q)+" Length to seek: "+IntToStr(succL));
       for(int i=FromN; i<=ToN; i++){
           writeln(&vsh1,"PosNn: "+IntToStr(i)+"checking... ");
           if(Array1DSuccIsAtPos(data, succ, curL, i, succL, Equals, &vsh1)){
               Ns.push_back(i);
               writeln(&vsh1,"Found!");
           }else  writeln(&vsh1,"NOT found");
        }
    }
    writeln(&vsh1, "Seek(1D in 1D) finishes working");
    return Ns;
}
/*
//
template <typename T> std::vector<int>Array1DSeek(const std::vector<T>&data, const std::vector<T>&succ, bool (*Equals)(const T& x1, const T& x2, int params)=EqualSimple, int FromN=1, int ToN=0, TValsShowHide*vsh=NULL){
    std::vector<int> Ns;
    int Q=data.size();//curL;
    int succL=succ.size();
    writeln(vsh, "Seek(1D in 1D) starts working");
    if(Q>0){
        if(ToN<1){
            ToN=Q;
            if(FromN>ToN){
                FromN=ToN;
            }
            if(FromN<1){
                FromN=1;
            }
            writeln(vsh,"Length to seek in: "+IntToStr(Q)+" Length to seek: "+IntToStr(succL));
            for(int i=FromN; i<=ToN; i++){
                writeln(vsh,"PosNn: "+IntToStr(i)+"checking... ");
                if(Array1DSuccIsAtPos(data, succ, curL, i, succL, Equals, vsh)){
                    Ns.push_back(i);
                    writeln(vsh,"Found!");
                }else  writeln(vsh,"NOT found");
            }
        }
    }
    writeln(vsh, "Seek(1D in 1D) finishes working");
    return Ns;
}
template <typename T> std::vector<int>Array1DSeekSimply(T*data, T*succ, int curL, int succL, int FromN=1, int ToN=0, int param=0, TValsShowHide*vsh=NULL){
    std::vector<int> Ns;
    int Q=curL, j;
    writeln(vsh, "Seek(1D in 1D) starts working");
    if(ToN==0) ToN=Q;
    if(FromN<=ToN){
        if(ToN==Q)ToN=ToN-succL+1;
        writeln(vsh,"Length to seek in: "+IntToStr(Q)+" Length to seek: "+IntToStr(succL)+" order: direct");
        for(int i=FromN; i<=ToN; i++){
            writeln(vsh,"PosNn: "+IntToStr(i)+"checking... ");
            if(Array1DSuccIsAtPos(data, succ, curL, i, succL, Equals, vsh)){
                Ns.push_back(i);
                writeln(vsh,"Found!");
            }else  writeln(vsh,"NOT found");
        }
    }else{
        if(ToN==1) ToN=succL;
        writeln(vsh,"Length to seek in: "+IntToStr(Q)+" Length to seek: "+IntToStr(succL)+" order: reversed");
        for(int i=FromN; i>=ToN; i--){
            j=i-succL+1;
            writeln(vsh,"PosNn: "+IntToStr(j)+"("+IntToStr(i)+")"+" checking... ");
            if(Array1DSuccIsAtPos(data, succ, curL, j, succL, Equals, vsh)){
                Ns.push_back(i);//or j
                writeln(vsh,"Found!");
            }else  writeln(vsh,"NOT found");
        }
    }
    writeln(vsh, "Seek(1D in 1D) finishes working");
    return Ns;
}
template <typename T> std::vector<int>Array1DSeekSimply(const std::vector<T>&data, const std::vector<T>&succ, int FromN=1, int ToN=0, int param=0, TValsShowHide*vsh=NULL){
    std::vector<int> Ns;
    int Q=data.size(), succL=succ.size(), j;
    writeln(vsh, "Seek(1D in 1D) starts working");
    if(ToN==0) ToN=Q;
    if(FromN<=ToN){
        if(ToN==Q)ToN=ToN-succL+1;
        writeln(vsh,"Length to seek in: "+IntToStr(Q)+" Length to seek: "+IntToStr(succL)+" order: direct");
        for(int i=FromN; i<=ToN; i++){
            writeln(vsh,"PosNn: "+IntToStr(i)+"checking... ");
            if(Array1DSuccIsAtPos(data, succ, curL, i, succL, Equals, vsh)){
                Ns.push_back(i);
                writeln(vsh,"Found!");
            }else  writeln(vsh,"NOT found");
        }
    }else{
        if(ToN==1) ToN=succL;
        writeln(vsh,"Length to seek in: "+IntToStr(Q)+" Length to seek: "+IntToStr(succL)+" order: reversed");
        for(int i=FromN; i>=ToN; i--){
            j=i-succL+1;
            writeln(vsh,"PosNn: "+IntToStr(j)+"("+IntToStr(i)+")"+" checking... ");
            if(Array1DSuccIsAtPos(data, succ, curL, j, succL, Equals, vsh)){
                Ns.push_back(i);//or j
                writeln(vsh,"Found!");
            }else  writeln(vsh,"NOT found");
        }
    }
    writeln(vsh, "Seek(1D in 1D) finishes working");
    return Ns;
}
template <typename T> int Array1DSeekFirst(T*data, T val, int curL, bool (*Equals)(const T& x1, const T& x2, int params)=EqualSimple, int FromN=1, int ToN=0){
    int N=0;
    std::vector<int> Ns=Array1DSeek(data, val, curL, Equals, FromN, ToN);
    if(Ns.size()>0){
        N=Ns[1-1];
    }
    return N;
}
template <typename T> int Array1DSeekLast(T*data, T val, int curL, bool (*Equals)(const T& x1, const T& x2, int params)=EqualSimple, int FromN=1, int ToN=0){
    int N=0;
    std::vector<int> Ns=Seek(data, val, curL, Equals, FromN, ToN);
    int count=Ns.size();
    if(count>0){
        N=Ns[count-1];
    }
    return N;
}
template <typename T> int Array1DSeekFirst(T*data, T*succ, int succL, int curL, bool (*Equals)(const T& x1, const T& x2, int params)=EqualSimple, int FromN=1, int ToN=0){
    int N=0;
    std::vector<int> Ns=Seek(data, succ, succL,curL, Equals, FromN, ToN);
    if(Ns.size()>0){
        N=Ns[1-1];
    }
    return N;
}
template <typename T> int Array1DSeekLast(T*data, T*succ, int succL, int curL, bool (*Equals)(const T& x1, const T& x2, int params)=EqualSimple, int FromN=1, int ToN=0){
    int N=0;
    std::vector<int>  Ns=Seek(data, succ, succL,curL, Equals,  FromN, ToN);
    int count=Ns.size();
    if(count>0){
        N=Ns[count-1];
    }
    return N;
}

template <typename T> int Array1DSeekFirst(T*data, T val, int curL,  int FromN=1, int ToN=0, int param=0){
    int N=0;
    std::vector<int> Ns=Array1DSeek(data, val, curL, FromN, ToN, param);
    if(Ns.size()>0){
        N=Ns[1-1];
    }
    return N;
}
template <typename T> int Array1DSeekLast(T*data, T val, int curL, int FromN=1, int ToN=0, int param=0){
    int N=0;
    std::vector<int> Ns=Seek(data, val, curL, Equals, FromN, ToN, param);
    int count=Ns.size();
    if(count>0){
        N=Ns[count-1];
    }
    return N;
}
template <typename T> int Array1DSeekFirst(T*data, T*succ, int succL, int FromN=1, int ToN=0, int param=0){
    int N=0;
    std::vector<int> Ns=Seek(data, succ, succL, curL, FromN, ToN, param);
    if(Ns.size()>0){
        N=Ns[1-1];
    }
    return N;
}
template <typename T> int Array1DSeekLast(T*data, T*succ, int succL, int curL, int FromN=1, int ToN=0, int param=0){
    int N=0;
    std::vector<int>  Ns=Seek(data, succ, succL,curL, FromN, ToN, param);
    int count=Ns.size();
    if(count>0){
        N=Ns[count-1];
    }
    return N;
}

template <typename T> void Array1DShowToConsole(T*x, int L, char* delim=" \0"){
    for(int i=1; i<=L; i++){
        std::cout<<x[i-1];
        std::cout<<delim;
    }
    std::cout<<std::endl;
}
*/

//=== wa 1D arr =======================================================================================================================================

//std::vector<double>NumbersSubRow(std::vector<double> x, int LreqGiven=0, int FromN=1, int QDefaultBefore=0, double DefaultValParam=0);
//std::vector<int>NumbersSubRow(std::vector<int> x, int LreqGiven=0, int FromN=1, int QDefaultBefore=0, int DefaultValParam=0);
//std::vector<double>NumbersSubRow(std::vector<double> x, int LreqGiven, int FromN=1, int QDefaultBefore=0, double*rowDflt=NULL);
//std::vector<int>NumbersSubRow(std::vector<int> x, int LreqGiven, int FromN=1, int QDefaultBefore=0, int*rowDflt=NULL);

//-----------------------------------------------------------------------------------------------------------------------------------------------------

int Array2DSizeBasedOn1D_PosByCoords(int ExtRowN, int IneRowN, int QExtRows, int*Ls);
int Array2DSizeBasedOn1D_PosByCoords(int ExtRowN, int IneRowN, std::vector<int>Ls);
void Array2DSizeBasedOn1D_CoordsByPos(int PosN, std::vector<int>Ls, int&ExtRowN, int&IneRowN);
int Array2DSizeBasedOn1D_ExtRowNByPos(int PosN, std::vector<int>Ls);
int Array2DSizeBasedOn1D_IneRowNByPos(int PosN, std::vector<int>Ls);
int Array2DSizeBasedOn1D_ExtRowNByPos(int PosN, int QExtRows, int*Ls);
int Array2DSizeBasedOn1D_IneRowNByPos(int PosN, int QExtRows, int*Ls);

//== wa 1D arrays === wi 2D size =======================================================================================================================

class Int2DVectorOfSizeOrCoords{
public:
    int d1, d2;
    Int2DVectorOfSizeOrCoords(int d1=0, int d2=0);

};

class Array2DSize
{
    int*Ls;
    int QExtRows, IneRowsLength;
public:
    Array2DSize();
    Array2DSize(const Array2DSize&obj);
    Array2DSize(int QExtRows, int  IneRowsLength, bool isNotVaria=true);
    Array2DSize(int QExtRows, int*Ls);//, int isNotVariaAndSetByMin1Max2=1);
    Array2DSize operator =(const Array2DSize&obj);
    ~Array2DSize();
    void Assign(const Array2DSize&obj);
    void Construct();
    void SetNull();
    void SetOne(bool isVaria=false);
    bool ExtRowNBelongsTo(int extRowN)const;
    bool PosBelongsTo(int extRowN, int ineRowN);
    void Set(int QExtRows, int IneRowsLength, bool preserveType=false);
    void Set(int QExtRows,int*Ls);
    void Set1ExtRowVaria(int L);
    void SetVaria(int QExtRows, int IneRowsLength);
    void ConvertToRectIfAllLengthesAreEqual();
    //void Add(int L, int ifConstAndNotEqual_Ignore0AddExisting1SetAll2=2);
    //void Ins(int N, int L);
    void AddExt(int L=0);
    void InsExt(int N=0, int L=0);
    void DelExt(int N);
    void DelIneRow();
    void AddOrInsIneRow();
    //void AddOrInsTo(int L=0, int N_0ToAllPreserveType_GTT0AndMakeVaria=0);
    //void DelFrom(int L, int N_0ToAllPreserveType_GTT0AndMakeVaria=0);
    void AddOrInsTo(int N, int L);
    void DelFrom(int N, int L);
    int GetQExtRows() const;
    bool GetIfIsVariaLength()const;
    bool isVariaLength();
    bool isRectangular() const;
    int GetLength(int N=0) const;
    int GetMinLength()const;
    int GetMaxLength()const;
    void SwapExtRows(int N1, int N2);
    void Reverse();
    void SetQExtRows(int N);
    void SetIneRowsLength(int L_Minus1MinMinus2Max=-1, bool preserveTypeIfVaria=false);
    //void SetLength(int L, int N_0ToAllPreserveType_GTT0AndMakeVaria_Minus1ToAllMakeNotVaria=0);
    void SetLength(int L, int N=0);
    QString GetAsQString(QString delim=", ") const;
    std::string GetAsStdString(std::string delim=", ") const;
    void ShowToConsole() const;
    bool isSame(const Array2DSize&obj) const;
    int GetIneRowsLength(bool EvenIfVarButAllEqual=false) const;
    void Transpose();
    bool IsForTranspose() const;
    int LengthFullOfIneRowN(int N) const;//from start to 1st stop
    int GetPositionByCoordsIfRealizedAs1D(int ExtRowN, int IneRowN);
    void GetCoordsByPositionIfRealizedAs1D(int PosN, int&ExtRowN, int&IneRowN);
    int GetExtRowNByPositionIfRealizedAs1D(int PosN);
    int GetIneRowNByPositionIfRealizedAs1D(int PosN);
    std::vector<int>GetLengthes()const;
    int GetQElements();
};


/*

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
};
*/



//=== wa 2D size === wi some Arr Fns ========================================================================

template<typename T> Array2DSize GetSizeOf2DVector(std::vector<std::vector<T>>X, bool VarIfAllEqual=false){
    Array2DSize size;
    int Q, L, Lcur;
    int *Ls=NULL;
    bool equal=true;
    Q=X.size();
    if(Q==0){
        L=0;
        size.Set(Q, L);
    }else{
        L=X[1-1].size();
        Ls=new int[Q];
        for(int i=1; i<=Q; i++){
            Lcur=X[i-1].size();
            Ls[i-1]=Lcur;
            //Ls[i-1]=X[i-1].size();
            if(Ls[i-1]!=L) equal=false;
        }
     if(equal && !VarIfAllEqual){
            size.Set(Q, L);
        }else{
            size.Set(Q,Ls);
        }
        delete[]Ls;
    }
    return size;
}
template<typename T> int GetMinLengthOfVector2DExtRows(std::vector<std::vector<T>>X){
    int Q=X.size();
    int L, mn, mx;
    for(int i=1; i<=Q; i++){
        L=X[i-1].size();
        if(i==1 || (i>1 && L<mn)) mn=L;
        if(i==1 || (i>1 && L>mx)) mx=L;
    }
    return mn;
}
template<typename T> int GetMaxLengthOfVector2DExtRows(std::vector<std::vector<T>>X){
    int Q=X.size();
    int L, mn, mx;
    for(int i=1; i<=Q; i++){
        L=X[i-1].size();
        if(i==1 || (i>1 && L<mn)) mn=L;
        if(i==1 || (i>1 && L>mx)) mx=L;
    }
    return mx;
}
template<typename T> int GetLengthFullOfIneRowOfVector2D(std::vector<std::vector<T>>X, int N){
    int Q=X.size(), Lmax=GetMaxLengthOfVector2DExtRows(X), Lmin=GetMinLengthOfVector2DExtRows(X), Lcur, RowL=0;
    if(Q>0){
        for(int i=1; i<=Q; i++){
            Lcur=X[i-1].size();
            if(Lcur<=N){
                RowL++;
            }
        }
    }
    return RowL;
}
template<typename T> bool Vector2DIsForTranspose(std::vector<std::vector<T>>X) {
    int Q=X.size(), Lmax=GetMaxLengthOfVector2DExtRows(X), Lmin=GetMinLengthOfVector2DExtRows(X);
    int L1, L2;
    bool b=true;
    if(Q>0){
        L1=GetLengthFullOfIneRowOfVector2D(X, 1);
        for(int j=2; j<=Lmax; j++){
            L2=GetLengthFullOfIneRowOfVector2D(X, j);
            if(L2>L1){
                b=false;
            }
            L1=L2;
        }
    }
    return b;
}

//=== wa some Arr Fns === wi 2D arr ========================================================================

template<typename T> void Array2DSetSize_DfltNeeded(T**&X, const Array2DSize&size1, const Array2DSize&size2, T DfltValParam){
    int L1, L2, minOfL, maxOfL, Q1=size1.GetQExtRows(), Q2=size2.GetQExtRows(), Qmin= Q1<=Q2 ? Q1 : Q2;
    T**Y=NULL;
    if(Q2==0){
        if(X!=NULL){
            for(int i=1; i<=Q1; i++){
                delete[]X;
            }
         X=NULL;
        }
    }else{
        Y=new T*[Q2];
        for(int i=1; i<=Q2; i++){
            L2=size2.GetLength(i);
            Y[i-1]=new T[L2];
        }
        if(X!=NULL){
            for(int i=1; i<=Qmin; i++){
                L1=size1.GetLength(i);
                L2=size2.GetLength(i);
                minOfL = L1<=L2 ? L1 : L2;
                if(L1!=0 && L2!=0){
                    for(int j=1; j<=minOfL; j++){
                        Y[i-1][j-1]=X[i-1][j-1];
                    }
                    for(int j=minOfL+1; j<=L2; j++){
                        Y[i-1][j-1]=DfltValParam;
                    }
                }
                for(int i=1; i<=Q1; i++){
                    delete[]X;
                }
                X=NULL;
            }
        }else{
            for(int i=1; i<=Qmin; i++){
                L2=size2.GetLength(i);
                for(int j=1; j<=L2; j++){
                    Y[i-1][j-1]=DfltValParam;
                }
            }
        }
        for(int i=Qmin+1; i<=Q2; i++){
            L2=size2.GetLength(i);
            for(int j=1; j<=L2; j++){
                Y[i-1][j-1]=DfltValParam;
            }
        }
        //
        X=Y;
        //size1=size2;
    }
}

template<typename T> void Array2DSetSize(T**&X, Array2DSize&size1, const Array2DSize&size2, T*DfltValParam=NULL, bool PreserveVals=true){//16
    T**Y=NULL;
    //T DefaultVal=GetDefaultVal(DfltValParam);
    T DefaultVal;
    int Q1=size1.GetQExtRows(), Q2=size2.GetQExtRows(),  minL, L1, L2, minQ=Q1<=Q2?Q1:Q2;
    if(DfltValParam!=NULL){
        DefaultVal=(*(DfltValParam));
    }
    if(size1.isSame(size2)==false){
        //1 set size
        Y=new T*[Q2];
        for(int i=1; i<=Q2; i++){
            L2=size2.GetLength(i);
            Y[i-1]=new T[L2];
        }
        //2 assign, preserving vals
        if(X!=NULL && PreserveVals){
            for(int i=1; i<=minQ; i++){
                L1=size1.GetLength(i);
                L2=size2.GetLength(i);//hin wa size1 et err - acc viol
                minL=L1<=L2?L1:L2;
                for(int j=1; j<=minL; j++){
                    Y[i-1][j-1]=X[i-1][j-1];
                }
                for(int j=minL+1; j<=L2; j++){
                    Y[i-1][j-1]=DefaultVal;
                }
                //if(minL==L1){
                //    for(int j=minL+1; j<=L2; j++){
                //        Y[i-1][j-1]=X[i-1][j-1];
                //    }
                //}else{
                //    for(int j=minL+1; j<=L2; j++){
                //        Y[i-1][j-1]=DefaultVal;
                //    }
                //}
            }
            for(int i=minQ+1; i<=Q2; i++){
                 for(int j=1; j<=size2.GetLength(i); j++){
                     Y[i-1][j-1]=DefaultVal;
                }
            }
            //if(minQ==Q1){
            //    for(int i=minQ+1; i<=Q2; i++){
            //        for(int j=1; j<=size2.GetLength(i); j++){
            //             Y[i-1][j-1]=X[i-1][j-1];
            //        }
            //     }
            //}else{
            //    for(int i=minQ+1; i<=Q2; i++){
            //         for(int j=1; j<=size2.GetLength(i); j++){
            //             Y[i-1][j-1]=DefaultVal;
            //        }
            //     }
            //}
        }else if(DfltValParam!=NULL){
            for(int i=1; i<=Q2; i++){
                L2=size2.GetLength(i);
                for(int j=1; j<=L2; j++){
                    Y[i-1][j-1]=DefaultVal;
                }
            }
        }
        //3 del old
        if(X!=NULL){
            for(int i=1; i<=Q1; i++){
                delete[]X[i-1];
            }
            delete[]X;
        }
        //4 set new
        X=Y;
        size1=size2;
    }
}
template<typename T> void Array2DSetSize(T**&X, const Array2DSize&size1, int Q, int*Ls, T*DfltValParam=NULL, bool PreserveVals=true){//16
    Array2DSize size2;
    size2.Set(Q, Ls);
    Array2DSetSize(X, size1, size2, DfltValParam, PreserveVals);
}
template<typename T> void Array2DSetSize(T**&X, Array2DSize&size1, int Q, int IneRowslength, T*DfltValParam=NULL, bool PreserveVals=true){//16
    Array2DSize size2;
    size2.Set(Q, IneRowslength);
    //template<typename T> void Array2DSetSize(T**&X, Array2DSize&size1, const Array2DSize&size2, T*DfltValParam=NULL, bool PreserveVals=true){//16p
    Array2DSetSize(                               X,              size1,                   size2,   DfltValParam,           PreserveVals);
 }
template<typename T> void Array2DSetSize(std::vector<std::vector<T>>&X, const Array2DSize&size, T*DfltValParam=NULL){//16
    int Q1=X.size(), Q2=size.GetQExtRows(), L1, L2, QMin=Q1<=Q2?Q1:Q2, Lmin;
    T DefaultVal;//=GetDefaultVal(DfltValParam);
    if(DfltValParam!=NULL)DefaultVal=(*(DfltValParam));
    Array2DSize sizeOld=GetSizeOf2DVector(X);
    X.resize(Q2);
    for(int i=1; i<=Q2; i++){
        L2=size.GetLength(i);
        X[i-1].resize(L2);
    }
    if(DfltValParam!=NULL){
        for(int i=1; i<=QMin; i++){
            L1=sizeOld.GetLength(i);
            L2=size.GetLength(i);
            Lmin=L1<=L2?L1:L2;
            for(int j=Lmin+1; j<=L2; j++){
                X[i-1][j-1]=DefaultVal;
            }
        }
        for(int i=QMin+1; i<=Q2; i++){
            L2=size.GetLength(i);
            for(int j=1; j<=L2; j++){
                X[i-1][j-1]=DefaultVal;
            }
        }
    }
}
template<typename T> void Array2DSetSize(std::vector<std::vector<T>>&X, int Q, int*Ls, T*DfltValParam=NULL){//16
     Array2DSize size;
     size.Set(Q, Ls);
     Array2DSetSize(X, size, DfltValParam);
}
template<typename T> void Array2DSetSize(std::vector<std::vector<T>>&X, int Q, int IneRowslength, T*DfltValParam=NULL){//16
     Array2DSize size;
     int QOld=X.size(), LcurOld, LcurNew;
     size.Set(Q, IneRowslength);
     Array2DSetSize(X, size, DfltValParam);
}
//set ext row
/*
template<typename T> void Array2DSetExtRowSimple(T**&X, Array2DSize&size, int ExtRowN, int wholeRowL=0, T*rowParam=NULL, bool KeepIfRect=true, T*DefaultValParam=NULL){
    T**Y=NULL;
    std::vector<T>row;
    int Q=size.GetQExtRows(), LERN, Lmin=size.GetMinLength(), Lmax=size.GetMaxLength(), L1;
    int Lreq, preN1_=0, preN2_=0, ownN1_=0, ownN2=0, ForOwnN1_=0, ForOwnN2=0, postN1=0, postN2_=0, Lini=0, LreqGiven=0, whatFromN=1, QDefaultValuesBefore=0;
    if(ExtRowN<0){
        ExtRowN=Q+ExtRowN+1;
    }
    if(ExtRowN>=1 && ExtRowN<=Q){
        LERN=size.GetLength(ExtRowN);
        if(Lmin==Lmax && wholeRowL==0  && rowParam!=NULL){
            wholeRowL=Lmin;
        }
        if(Lmin==Lmax && KeepIfRect && wholeRowL>0){
            Lreq=Lmin;
            L1=Lreq;
        }else{
            Lreq=0;
            L1=wholeRowL;
        }
        //
        row=Array1DGetSubArray_byLs_PreMark(rowParam, wholeRowL, Lreq, 1, 0, DefaultValParam);
        //
        //if()
        //
        //void Calc_Array1DSubArray_byLs_MarkupNs_vFmls(int&LreqCalcd, int&preN1_, int&preN2_, int&ownN1_, int&ownN2, int&ForOwnN1_, int&ForOwnN2, int&postN1, int&postN2_, int Lini, int LreqGiven=0, int whatFromN=1, int QDefaultValuesBefore=0);

        Calc_Array1DSubArray_byLs_MarkupNs_vFmls(LreqCalcd, preN1_, preN2_, ownN1_, ownN2, ForOwnN1_, ForOwnN2, postN1, postN2_, Lini, LreqGiven, whatFromN, QDefaultValuesBefore);
        //

        //
        if(LERN!=L1){
            delete[]X[ExtRowN-1];
            X[ExtRowN]=new T[L1];
            //
            size.SetLength(L1, ExtRowN);
        }
        for(int i=1; i<=L1; i++){
             X[ExtRowN][i-1]=row[i-1];
        }
    } 
}
template<typename T> void Array2DSetExtRowSimple(T**&X, Array2DSize&size, int ExtRowN, std::vector<T>rowParam, bool KeepIfRect=true, T*DefaultValParam=NULL, bool DefaultIsSpecNotOwn=false){
    Array2DSetExtRowSimple(X, size, ExtRowN, rowParam.size(), rowParam.data(), KeepIfRect, DefaultValParam);
}
template<typename T> void Array2DSetExtRowSimple(std::vector<std::vector<T>>&X, int ExtRowN, std::vector<T>rowParam, bool KeepIfRect=true, T*DefaultValParam=NULL){
    Array2DSize size;
    size=GetSizeOf2DVector(X);
    std::vector<T>row;
    int Q=size.GetQExtRows(), Lmin=size.GetMinLength(), Lmax=size.GetMaxLength(), LERN, L1, wholeRowL=rowParam.size(), Lreq;
    if(ExtRowN<0){
        ExtRowN=Q+ExtRowN+1;
    }
    if(ExtRowN>=1 && ExtRowN<=Q){
        LERN=size.GetLength(ExtRowN);
        if(Lmin==Lmax && wholeRowL==0){
            wholeRowL=Lmin;
        }
        if(Lmin==Lmax && KeepIfRect && wholeRowL!=Lmin){
            Lreq=Lmin;
            L1=Lreq;
        }else{
            Lreq=0;
            L1=wholeRowL;
        }
        row=Array1DGetSubArray_byLs_PreMark(rowParam, Lreq, 1, 0,  DefaultValParam); //add default vals, id S'zq Lmin or cut if S'kq Lmax
        //if(LERN!=L1){
        //    delete[]X[ExtRowN-1];
        //    X[ExtRowN]=new T[L1];
        //    //
        //    size.SetLength(L1, ExtRowN);
        //}
        X[ExtRowN-1]=row;
    }
}
template<typename T> void Array2DSetExtRowSimple(std::vector<std::vector<T>>&X, int ExtRowN, int wholeRowL=0, T*rowParam=NULL, bool KeepIfRect=true, T*DefaultValParam=NULL){
    Array2DSize size;
    size=GetSizeOf2DVector(X);
    std::vector<T>row;
    int Q=size.GetQExtRows(), Lmin=size.GetMinLength(), Lmax=size.GetMaxLength(), LERN, L1, Lreq;
    if(ExtRowN<0){
        ExtRowN=Q+ExtRowN+1;
    }
    if(ExtRowN>=1 && ExtRowN<=Q){
        LERN=size.GetLength(ExtRowN);
        if(Lmin==Lmax && wholeRowL==0 && rowParam!=NULL){
            wholeRowL=Lmin;
        }
        if(Lmin==Lmax && KeepIfRect && wholeRowL>0){
            Lreq=Lmin;
            L1=Lreq;
        }else{
            Lreq=0;
            L1=wholeRowL;
        }
        row=Array1DGetSubArray_byLs_PreMark(rowParam, wholeRowL, Lreq, 1, 0, DefaultValParam); //add default vals, id S'zq Lmin or cut if S'kq Lmax
        //if(LERN!=L1){
        //    delete[]X[ExtRowN-1];
        //    X[ExtRowN]=new T[L1];
        //    //
        //    size.SetLength(L1, ExtRowN);
        //}
        X[ExtRowN-1]=row;
    }
}
*/
template<typename T> void Array2DSetRect(T**&X, Array2DSize&size, int QIneRows_m1_1st_m2_last_m3_min_Other_max=-4, int QExtRows=0, T*DefaultValParam=NULL, bool PreserveVals=true){
    int Q=size.GetQExtRows(), Lmin=size.GetMinLength(), Lmax=size.GetMaxLength(), L1=size.GetLength(1), LLast=size.GetLength(Q);
    if(QIneRows_m1_1st_m2_last_m3_min_Other_max==-1){
         QIneRows_m1_1st_m2_last_m3_min_Other_max=L1;
    }else if(QIneRows_m1_1st_m2_last_m3_min_Other_max==-2){
        QIneRows_m1_1st_m2_last_m3_min_Other_max=LLast;
   }else if(QIneRows_m1_1st_m2_last_m3_min_Other_max==-3){
        QIneRows_m1_1st_m2_last_m3_min_Other_max=Lmin;
   }else if(QIneRows_m1_1st_m2_last_m3_min_Other_max<0){
        QIneRows_m1_1st_m2_last_m3_min_Other_max=Lmax;
   }
    if(QExtRows==0){
        if(X==NULL){
            QExtRows=1;
        }else{
            QExtRows=size.GetQExtRows();
        }
        Array2DSetSize(X, size, QExtRows, QIneRows_m1_1st_m2_last_m3_min_Other_max, DefaultValParam, PreserveVals);
    }
}
template<typename T> void Array2DSetRect(std::vector<std::vector<T>>&X, int QIneRows, int QExtRows=0, T*DefaultValParam=NULL, bool PreserveVals=true){
    if(QExtRows==0){
        if(X==NULL){
            QExtRows=1;
        }else{
            QExtRows=X.size();
        }
        Array2DSetSize(X, QExtRows, QIneRows, DefaultValParam);
    }
}
template<typename T> void Array2DAssign(std::vector<std::vector<T>>&X, int L, T*data=NULL, bool QExtRowsNotIneRowLength=true, T*DefaultValParam=NULL){
    T DefaultVal, CurVal;
    std::vector<T>row;
    if(DefaultValParam!=NULL){
        DefaultVal=(*(DefaultValParam));
    }
    X.clear();
    if(L>0){
        if(QExtRowsNotIneRowLength){
            if(row!=NULL){
                for(int i=1; i<=L; i++){
                    row.clear();
                    CurVal=data[i-1];
                    row.push_back(CurVal);
                    X.push_back(row);
                }
            }else{
                for(int i=1; i<=L; i++){
                    row.clear();
                    CurVal=DefaultVal;
                    row.push_back(CurVal);
                    X.push_back(row);
                }
            }
        }else{
            if(row!=NULL){
                row.clear();
                for(int i=1; i<=L; i++){
                    CurVal=data[i-1];
                    row.push_back(CurVal);
                }
            }else{
                for(int i=1; i<=L; i++){
                    CurVal=DefaultVal;
                    row.push_back(CurVal);
                }
            }
            X.push_back(row);
        }
    }
}

//SetExtRow

template<typename T> void Array2DSetExtRowN(T**&X, Array2DSize&size, int ExtRowN, int wholeRowL=0, T*rowParam=NULL, int whatFromN=1, int QDefaultValuesBefore=0, bool KeepIfRect=true, T*DefaultValParam=NULL, bool DefaultIsSpecNotOwn=false){
    std::vector<T>rowVect;
    int preN1_, preN2_, ownN1_, ownN2, ForOwnN1_, ForOwnN2, postN1, postN2_,
        //Lold=wholeRowL,
        LreqGiven, LreqCalcd,
        //FromN,
        //DefaultValuesBefore
        //;
        minL, Q=size.GetQExtRows(), Lmin=size.GetMinLength(), Lmax=size.GetMaxLength(), LERN;
    T DefaultVal, CurVal;
    if(ExtRowN<0){
        ExtRowN=Q+ExtRowN+1;
    }
    if(ExtRowN>=1 && ExtRowN<=Q){
        LERN=size.GetLength(ExtRowN);
        if(wholeRowL==0){
            wholeRowL=LERN;
        }
        if(DefaultValParam!=NULL){
            DefaultVal=(*(DefaultValParam));
        }
        if(size.isRectangular() && KeepIfRect){
            //LreqGiven=size.GetIneRowsLength();//os irr!
            LreqGiven=Lmin;//or max, D s'egal
        }else{
            LreqGiven=0;
        }
        //void Calc_Array1DSubArray_byLs_MarkupNs_vFmls(int&LreqCalcd, int&preN1_, int&preN2_, int&ownN1_, int&ownN2, int&ForOwnN1_, int&ForOwnN2, int&postN1, int&postN2_, int Lini, int LreqGiven, int whatFromN, int QDefaultValuesBefore){
        //void Calc_Array1DSubArray_byLs_MarkupNs_vFmls(int&LreqCalcd, int&preN1_, int&preN2_, int&ownN1_, int&ownN2, int&ForOwnN1_, int&ForOwnN2, int&postN1, int&postN2_, int Lini, int LreqGiven=0, int whatFromN=1, int QDefaultValuesBefore=0);
        Calc_Array1DSubArray_byLs_MarkupNs_vFmls( LreqCalcd,     preN1_,    preN2_,     ownN1_,     ownN2,      ForOwnN1_,    ForOwnN2,      postN1,     postN2_, wholeRowL,       LreqGiven,       whatFromN,       QDefaultValuesBefore);
        //pre part
        if(preN2_>0){
            if(DefaultIsSpecNotOwn){
                minL = wholeRowL<=preN2_ ? wholeRowL : preN2_; //ab ini s'ver
                for(int i=preN1_; i<=minL; i++){
                    CurVal=X[ExtRowN-1][i-1];
                    rowVect.push_back(CurVal);
                }
                //for(int i=minL+1; i<=preN2_; i++){
                //    rowVect.push_back(DefaultVal);
                //}
            }//else{
            //    for(int i=preN1_; i<=preN2_; i++){
            //        rowVect.push_back(DefaultVal);
            //    }
            //}
        }
        //self part
        if(ownN2>0){
            for(int i=ownN1_; i<=ownN2; i++){
                CurVal=rowParam[i-1];
                rowVect.push_back(CurVal);
            }
        }
        //post part
        if(postN2_>0){
            minL = wholeRowL <= postN2_ ? wholeRowL : postN2_; // I men et hin S'ver
            if(DefaultIsSpecNotOwn){
                for(int i=postN1; i<=minL; i++){
                    CurVal=rowParam[i-1];
                    rowVect.push_back(CurVal);
                }
                //for(int i=minL+1; i<=postN2_; i++){
                //    rowVect.push_back(DefaultVal);
                //}
            }//else{//else row ha dflt vals
            //    for(int i=postN1; i<=postN2_; i++){
            //        rowVect.push_back(DefaultVal);
            //    }
            //}
            //for(int i=postN1; i<=postN2_; i++){
            //    if(DefaultIsSpecNotOwn && i<=wholeRowL){
            //        CurVal=rowParam[i-1];
            //    }else{
            //        CurVal=DefaultVal;
            //    }
            //    rowVect.push_back(CurVal);
            //}
        }
        //assign
        delete[]X[ExtRowN-1];
        X[ExtRowN-1]=new T[LreqCalcd];
        for(int i=1; i<=LreqCalcd; i++){
            X[ExtRowN-1][i-1]=rowVect[i-1];
        }
        //reflect
        size.SetLength(LreqCalcd, ExtRowN);
    }
}
template<typename T> void Array2DSetExtRowN(T**&X, Array2DSize&size, int ExtRowN, std::vector<T>rowParam, T*DefaultValParam=NULL, int FromN=1, int QDefaultValuesBefore=0, bool KeepIfRect=true, bool DefaultIsSpecNotOwn=false){
    //void Array2DSetExtRowN(T**&X, Array2DSize&size, int ExtRowN, int wholeRowL=0, T*rowParam=NULL, int whatFromN=1, int QDefaultValuesBefore=0, bool KeepIfRect=true, T*DefaultValParam=NULL, bool DefaultIsSpecNotOwn=false){
     Array2DSetExtRowN          (X,             size,     ExtRowN, rowParam.size(), rowParam.data(),         FromN,       QDefaultValuesBefore,        KeepIfRect,        DefaultValParam,           DefaultIsSpecNotOwn);
    //Array2DSetExtRowN(X, size, ExtRowN, rowParam, DefaultValParam, FromN, QDefaultValuesBefore, KeepIfRect);

}
template<typename T> void Array2DSetExtRowN(std::vector<std::vector<T>>&X, int ExtRowN, std::vector<T>rowParam, int whatFromN=1, int QDefaultValuesBefore=0, bool KeepIfRect=true, T*DefaultValParam=NULL, bool DefaultIsSpecNotOwn=false){//, int ifDiff_ByArr0ByNew1Min2Max3=0){
    std::vector<T>rowVect;
    Array2DSize size=GetSizeOf2DVector(X);
    int preN1_, preN2_, ownN1_, ownN2, ForOwnN1_, ForOwnN2, postN1, postN2_,
        //Lold=wholeRowL,
        LreqGiven, LreqCalcd,
        //FromN,
        //DefaultValuesBefore
        //;
        minL, Q=X.size(), LERN, Lmin=size.GetMinLength(), Lmax=size.GetMaxLength(), wholeRowL=rowParam.size();
    T DefaultVal, CurVal;
    if(ExtRowN<0){
        ExtRowN=Q+ExtRowN+1;
    }
    if(ExtRowN>=1 && ExtRowN<=Q){
        LERN=X[ExtRowN-1].size();
        //if(wholeRowL==0){
        //    wholeRowL=LERN;
        //}
        if(DefaultValParam!=NULL){
            DefaultVal=(*(DefaultValParam));
        }
        if(Lmin==Lmax && KeepIfRect){
            LreqGiven=Lmin;
        }else{
            LreqGiven=0;
        }
        Calc_Array1DSubArray_byLs_MarkupNs_vFmls( LreqCalcd,     preN1_,    preN2_,     ownN1_,     ownN2,      ForOwnN1_,    ForOwnN2,      postN1,     postN2_, wholeRowL,       LreqGiven,       whatFromN,       QDefaultValuesBefore);
        //pre part
        if(preN2_>0){
            if(DefaultIsSpecNotOwn){
                minL = wholeRowL<=preN2_ ? wholeRowL : preN2_;
                for(int i=preN1_; i<=minL; i++){
                    CurVal=X[ExtRowN-1][i-1];
                    rowVect.push_back(CurVal);
                }
                //for(int i=minL+1; i<=preN2_; i++){
                //    rowVect.push_back(DefaultVal);
                //}
            }//else{
            //    for(int i=preN1_; i<=preN2_; i++){
            //        rowVect.push_back(DefaultVal);
            //     }
            //}
        }
        //self part
        if(ownN2>0){
            for(int i=ownN1_; i<=ownN2; i++){
                CurVal=rowParam[i-1];
                rowVect.push_back(CurVal);
            }
        }
        //post part
        if(postN2_>0){
            minL = wholeRowL <= postN2_ ? wholeRowL : postN2_;
            if(DefaultIsSpecNotOwn){
                for(int i=postN1; i<=minL; i++){
                    CurVal=rowParam[i-1];
                    rowVect.push_back(CurVal);
                }
                //for(int i=minL+1; i<=postN2_; i++){
                //    rowVect.push_back(DefaultVal);
                //}
            }//else{
            //    for(int i=postN1; i<=postN2_; i++){
            //        rowVect.push_back(DefaultVal);
            //    }
            //}
            //for(int i=postN1; i<=postN2_; i++){
            //    if(DefaultIsSpecNotOwn && i<=wholeRowL){
            //        CurVal=rowParam[i-1];
            //    }else{
            //        CurVal=DefaultVal;
            //    }
            //    rowVect.push_back(CurVal);
            //}
        }
        //assign
        //delete[]X[ExtRowN-1];X[ExtRowN-1]=new T[LreqCalcd];
        //X.clear();
        //for(int i=1; i<=LreqCalcd; i++){
        //    CurVal=rowVect[i-1];
        //    X.push_back(CurVal);
        //}
        X[ExtRowN-1].clear();
        X[ExtRowN-1]=rowVect;
        //reflect
        //size.SetLength(LreqCalcd, ExtRowN);
    }
}
template<typename T> void Array2DSetExtRowN(std::vector<std::vector<T>>&X, int ExtRowN, int wholeRowL=0, T*rowParam=NULL, int FromN=1, int QDefaultValuesBefore=0, bool KeepIfRect=true, T*DefaultValParam=NULL, bool DefaultIsSpecNotOwn=false){//, int ifDiff_ByArr0ByNew1Min2Max3=0){
    //template <typename T> std::vector<T> Array1DAssign_VbyP(T*y=NULL, int Qy=0,  T*DfltValParam=NULL){
    std::vector<T>row=Array1DAssign_VbyP(wholeRowL, rowParam,  DefaultValParam);
    //Array2DSetExtRowN(std::vector<std::vector<T>>&X, int ExtRowN, std::vector<T>rowParam, int FromN=1, int DefaultValuesBefore=0, bool KeepIfRect=true, T*DefaultValParam=NULL, bool DefaultIsSpecNotOwn=false){
    Array2DSetExtRowN(                              X,     ExtRowN,                     row,      FromN,      QDefaultValuesBefore,           KeepIfRect,        DefaultValParam,      DefaultIsSpecNotOwn);
}
//
/*template<typename T> void Array2DSetExtRowN(T**&X, Array2DSize&size, int ExtRowN, int Lreq, int whatFromN=1, int QDefaultValuesBefore=0, T*DefaultValParam=NULL, bool KeepIfRect=true){
    std::vector<T>rowVect;
    int Q=size.GetQExtRows(), curL=size.GetLength(ExtRowN), minL, Lmin, Lmax;
    if(ExtRowN>=1 && ExtRowN<=Q){
        rowVect=Array1DGetSubArray_byLs_PreMark(X[ExtRowN-1], Lreq, whatFromN, QDefaultValuesBefore, DefaultValParam);
        //template<typename T> void Array2DSetExtRowN(T**&X, Array2DSize&size, int ExtRowN, std::vector<T>rowParam, T*DefaultValParam=NULL, int FromN=1, int QDefaultValuesBefore=0, bool KeepIfRect=true, bool DefaultIsSpecNotOwn=false){
        Array2DSetExtRowN(X, size, ExtRowN, rowVect, DefaultValParam, whatFromN, QDefaultValuesBefore, false, false);
        Lmin=size.GetMinLength(), Lmax=size.GetMaxLength();
        if(Lmin==Lmax && KeepIfRect){
            size.Set(Q, Lmin);
        }
    }
}
template<typename T> void Array2DSetExtRowN(std::vector<std::vector<T>>&X, int ExtRowN, int Lreq, int whatFromN=1, int QDefaultValuesBefore=0, T*DefaultValParam=NULL, bool KeepIfRect=true){
    std::vector<T>rowVect;
    Array2DSize size=GetSizeOf2DVector(X);
    int Q=size.GetQExtRows(), curL=size.GetLength(ExtRowN), minL, Lmin, Lmax;
    if(ExtRowN>=1 && ExtRowN<=Q){
        rowVect=Array1DGetSubArray_byLs_PreMark(X[ExtRowN-1], Lreq, whatFromN, QDefaultValuesBefore, DefaultValParam);
        //template<typename T> void Array2DSetExtRowN(std::vector<std::vector<T>>&X, int ExtRowN, std::vector<T>rowParam, int whatFromN=1, int QDefaultValuesBefore=0, bool KeepIfRect=true, T*DefaultValParam=NULL, bool DefaultIsSpecNotOwn=false){//, int ifDiff_ByArr0ByNew1Min2Max3=0){
        Array2DSetExtRowN(X, ExtRowN, rowVect, whatFromN, QDefaultValuesBefore, false, DefaultValParam, false);
    }
}*/
template<typename T> void Array2DSetExtRowN_DfltNeeded(T**&X, Array2DSize&size, int ExtRowN, T DefaultVal, int wholeRowL=0, T*rowParam=NULL, int whatFromN=1, int QDefaultValuesBefore=0, bool KeepIfRect=true, bool DefaultIsSpecNotOwn=false){
    //keep_if_rect seems exutz
    int QE=size.GetQExtRows(), Lmin=size.GetMinLength(), Lmax=size.GetMaxLength();
    std::vector<T>rowTmp;
    int Lreq=0;
    //Lreq;
    if(Lmin==Lmax && (wholeRowL==0 || wholeRowL==Lmin)){
       wholeRowL=Lmin;
    }else{
        size.SetLength(wholeRowL, ExtRowN);
    }
    if(ExtRowN>=1 && ExtRowN<=QE){
        if(wholeRowL>0){
            if(DefaultIsSpecNotOwn){
                rowTmp=Array1DGetSubArray_byLs_PreMark_v1(rowParam, wholeRowL, X[ExtRowN-1], Lreq, whatFromN, QDefaultValuesBefore, wholeRowL);
            }else{
                rowTmp=Array1DGetSubArray_byLs_PreMark_v1(rowParam, wholeRowL, DefaultVal, Lreq, whatFromN, QDefaultValuesBefore);
            }
            Array1DAssign(X[ExtRowN-1], rowTmp);
        }else{
            if(X[ExtRowN-1]!=NULL){
                delete[]X[ExtRowN-1];
                X[ExtRowN-1]=NULL;
            }
        }
    }
}
template<typename T> void Array2DSetExtRowN_v1(T**&X, Array2DSize&size, int ExtRowN, int wholeRowL=0, T*rowParam=NULL, int whatFromN=1, int QDefaultValuesBefore=0, T*DefaultValParam=NULL, bool KeepIfRect=true, bool DefaultIsSpecNotOwn=false){
    //keep_if_rect seems exutz
    int QE=size.GetQExtRows(), Lmin=size.GetMinLength(), Lmax=size.GetMaxLength();
    std::vector<T>rowTmp;
    int Lreq=0;
    T DefaultVal;
    if(*DefaultValParam!=NULL)DefaultVal=(*(DefaultValParam));
    //Lreq;
    if(Lmin==Lmax && (wholeRowL==0 || wholeRowL==Lmin)){
       wholeRowL=Lmin;
    }else{
        size.SetLength(wholeRowL, ExtRowN);
    }
    if(ExtRowN>=1 && ExtRowN<=QE){
        if(wholeRowL>0){
            if(DefaultIsSpecNotOwn){
                rowTmp=Array1DGetSubArray_byLs_PreMark_v1(rowParam, wholeRowL, X[ExtRowN-1], Lreq, whatFromN, QDefaultValuesBefore, wholeRowL);
            }else{
                rowTmp=Array1DGetSubArray_byLs_PreMark_v1(rowParam, wholeRowL, DefaultVal, Lreq, whatFromN, QDefaultValuesBefore);
            }
            Array1DAssign(X[ExtRowN-1], rowTmp);
        }else{
            if(X[ExtRowN-1]!=NULL){
                delete[]X[ExtRowN-1];
                X[ExtRowN-1]=NULL;
            }
        }
    }
}
template<typename T> void Array2DSetExtRowN_v2(T**&X, Array2DSize&size, int ExtRowN, int wholeRowL=0, T*rowParam=NULL, int whatFromN=1, int QDefaultValuesBefore=0, T*DefaultValParam=NULL, bool KeepIfRect=true, bool DefaultIsSpecNotOwn=false){
    T DefaultVal;
    if(*DefaultValParam!=NULL)DefaultVal=(*(DefaultValParam));
    Array2DSetExtRowN_DfltNeeded(X, size, ExtRowN, DefaultVal, wholeRowL, rowParam, whatFromN, QDefaultValuesBefore, KeepIfRect, DefaultIsSpecNotOwn);
}

//
template<typename T> void Array2DSwapVals(T**&X, const Array2DSize& size, int ExtRowN1, int IneRowN1, int ExtRowN2, int IneRowN2){
    int Q=size.GetQExtRows(), L1, L2, Lmin=size.GetMinLength(), Lmax=size.GetMaxLength();
    T bufVal;
    if(ExtRowN1<0){
        ExtRowN1=Q+ExtRowN1+1;
    }
    if(ExtRowN2<0){
        ExtRowN2=Q+ExtRowN2+1;
    }
    if(ExtRowN1>=1 && ExtRowN1<=Q && ExtRowN2>=1 && ExtRowN2<=Q){
        L1=size.GetLength(ExtRowN1);
        L2=size.GetLength(ExtRowN2);
        if(IneRowN1<0){
            IneRowN1=L1+IneRowN1+1;
            IneRowN1=Lmax+IneRowN1+1;
        }
        if(IneRowN2<0){
            IneRowN2=L2+IneRowN2+1;
            IneRowN2=Lmax+IneRowN2+1;
        }
        if(IneRowN1>=1 && IneRowN1<=L1 && IneRowN1>=1 && IneRowN1<=L2){
            bufVal=X[ExtRowN1-1][IneRowN1-1];
            X[ExtRowN1-1][IneRowN1-1]=X[ExtRowN2-1][IneRowN2-1];
            X[ExtRowN2-1][IneRowN2-1]=bufVal;
        }
    }
}
/*template<typename T>Array2DSwapVals(std::vector<std::vector<T>>&X,int ExtRowN1, int IneRowN1, int ExtRowN2, int IneRowN2){
    Array2DSize size=GetSizeOf2DVector(X);
    int Q=size.GetQExtRows(), L1, L2, Lmin=size.GetMinLength(), Lmax=size.GetMaxLength();
    T bufVal;
    if(ExtRowN1<0){
        ExtRowN1=Q+ExtRowN1+1;
    }
    if(ExtRowN2<0){
        ExtRowN2=Q+ExtRowN2+1;
    }
    if(ExtRowN1>=1 && ExtRowN1<=Q && ExtRowN2>=1 && ExtRowN2<=Q){
        L1=size.GetLength(ExtRowN1);
        L2=size.GetLength(ExtRowN2);
        if(IneRowN1<0){
            IneRowN1=L1+IneRowN1+1;
            IneRowN1=Lmax+IneRowN1+1;
        }
        if(IneRowN2<0){
            IneRowN2=L2+IneRowN2+1;
            IneRowN2=Lmax+IneRowN2+1;
        }
        if(IneRowN1>=1 && IneRowN1<=L1 && IneRowN1>=1 && IneRowN1<=L2){
            bufVal=X[ExtRowN1-1][IneRowN1-1];
            X[ExtRowN1-1][IneRowN1-1]=X[ExtRowN2-1][IneRowN2-1];
            X[ExtRowN2-1][IneRowN2-1]=bufVal;
        }
    }
}*/
template<typename T>Array2DSwapVals(std::vector<std::vector<T>>&X, int ExtRowN1, int IneRowN1, int ExtRowN2, int IneRowN2){
    Array2DSize size;
    size=GetSizeOf2DVector(X);
    int Q=size.GetQExtRows(), L1, L2, Lmin=size.GetMinLength(), Lmax=size.GetMaxLength();
    T bufVal;
    if(ExtRowN1<0){
        ExtRowN1=Q+ExtRowN1+1;
    }
    if(ExtRowN2<0){
        ExtRowN2=Q+ExtRowN2+1;
    }
    if(ExtRowN1>=1 && ExtRowN1<=Q && ExtRowN2>=1 && ExtRowN2<=Q){
        L1=size.GetLength(ExtRowN1);
        L2=size.GetLength(ExtRowN2);
        if(IneRowN1<0){
            IneRowN1=L1+IneRowN1+1;
            IneRowN1=Lmax+IneRowN1+1;
        }
        if(IneRowN2<0){
            IneRowN2=L2+IneRowN2+1;
            IneRowN2=Lmax+IneRowN2+1;
        }
        if(IneRowN1>=1 && IneRowN1<=L1 && IneRowN1>=1 && IneRowN1<=L2){
            bufVal=X[ExtRowN1-1][IneRowN1-1];
            X[ExtRowN1-1][IneRowN1-1]=X[ExtRowN2-1][IneRowN2-1];
            X[ExtRowN2-1][IneRowN2-1]=bufVal;
        }
    }
}
//
template<typename T>Array2DSwapExtRows(T**&X, Array2DSize& size, int ExtRowN1, int ExtRowN2){
    T*bufPtr=NULL, *ptr1=NULL, *ptr2=NULL;
    int Q=size.GetQExtRows(), L1, L2, Lmin=size.GetMinLength(), Lmax=size.GetMaxLength();
    T bufVal;
    if(ExtRowN1<0){
        ExtRowN1=Q+ExtRowN1+1;
    }
    if(ExtRowN2<0){
        ExtRowN2=Q+ExtRowN2+1;
    }
    if(ExtRowN1>=1 && ExtRowN1<=Q && ExtRowN2>=1 && ExtRowN2<=Q){
        bufPtr=X[ExtRowN1-1];
        X[ExtRowN1-1]=X[ExtRowN2-1];
        X[ExtRowN2-1]=bufPtr;
        //
        //all id below s' uz try et debug, ob hic vrn narb'te, ma nu arb
        //
        //bufPtr=X[ExtRowN1-1];
        //X[ExtRowN1-1]=NULL;
        //X[ExtRowN1-1]=X[ExtRowN2-1];
        //X[ExtRowN2-1]=NULL;
        //X[ExtRowN2-1]=bufPtr;
        //
        //ptr1=X[ExtRowN1-1];
        //ptr2=X[ExtRowN2-1];
        //bufPtr=ptr1;
        //ptr1=ptr2;
        //ptr2=bufPtr;
        //X[ExtRowN1-1]=ptr1;
        //X[ExtRowN2-1]=ptr2;
        //
        //ptr1=X[ExtRowN1-1];
        //ptr2=X[ExtRowN2-1];
        //bufPtr=ptr1;
        //ptr1=NULL;
        //ptr1=ptr2;
        //ptr2=NULL;
        //ptr2=bufPtr;
        //X[ExtRowN1-1]=NULL;
        //X[ExtRowN2-1]=NULL;
        //X[ExtRowN1-1]=ptr1;
        //X[ExtRowN2-1]=ptr2;
        //
        L1=size.GetLength(ExtRowN1);
        L2=size.GetLength(ExtRowN2);
        size.SetLength(L1, ExtRowN2);
        size.SetLength(L2, ExtRowN1);
    }
}
template<typename T>Array2DSwapExtRows(std::vector<std::vector<T>>&X, int ExtRowN1, int ExtRowN2){
    std::vector<T>bufRow;
    Array2DSize size;
    size=GetSizeOf2DVector(X);
    int Q=size.GetQExtRows(), L1, L2, Lmin=size.GetMinLength(), Lmax=size.GetMaxLength();
    T bufVal;
    if(ExtRowN1<0){
        ExtRowN1=Q+ExtRowN1+1;
    }
    if(ExtRowN2<0){
        ExtRowN2=Q+ExtRowN2+1;
    }
    if(ExtRowN1>=1 && ExtRowN1<=Q && ExtRowN2>=1 && ExtRowN2<=Q){
        bufRow=X[ExtRowN1-1];
        X[ExtRowN1-1]=X[ExtRowN2-1];
        X[ExtRowN2-1]=bufRow;
    }
}
//
template<typename T>Array2DSwapIneRows(T**&X, const Array2DSize& size, int IneRowN1, int IneRowN2){
    int Q=size.GetQExtRows(), L1, L2, Lmin, Lmax;
    if(Q>0){
        Lmin=size.GetMinLength();
        Lmax=size.GetMaxLength();
        if(IneRowN1<0){
            IneRowN1=Lmax+IneRowN1+1;
        }
        if(IneRowN2<0){
            IneRowN2=Lmax+IneRowN2+1;
        }
        if(IneRowN1>=1 && IneRowN1<=Lmin && IneRowN2>=1 && IneRowN2<=Lmin){
            for(int i=1; i<=Q; i++){
                Array2DSwapVals(X, size, i, IneRowN1, i, IneRowN2);
            }
        }
    }
}
template<typename T>Array2DSwapIneRows(std::vector<std::vector<T>>&X, int IneRowN1, int IneRowN2){
    Array2DSize size;
    size=GetSizeOf2DVector(X);
    int Q=size.GetQExtRows(), L1, L2, Lmin, Lmax;
    if(Q>0){
        Lmin=size.GetMinLength();
        Lmax=size.GetMaxLength();
        if(IneRowN1<0){
            IneRowN1=Lmax+IneRowN1+1;
        }
        if(IneRowN2<0){
            IneRowN2=Lmax+IneRowN2+1;
        }
        if(IneRowN1>=1 && IneRowN1<=Lmin && IneRowN2>=1 && IneRowN2<=Lmin){
            for(int i=1; i<=Q; i++){
                Array2DSwapVals(X, size, i, IneRowN1, i, IneRowN2);
            }
        }
    }
}
//Reverse ext rows
template<typename T>Array2DReverseExtRows(T**&X, Array2DSize& size){
    int Q=size.GetQExtRows(), N1, N2, Nmid;
    if(Q%2==0){
        Nmid=Q/2;
    }else{
        Nmid=(Q-1)/2;
    }
    for(N1=1; N1<=Nmid; N1++){
        N2=Q-N1+1;
        Array2DSwapExtRows(X, size, N1, N2);
    }
}
template<typename T>Array2DReverseExtRows(std::vector<std::vector<T>>&X){
    Array2DSize size;
    size=GetSizeOf2DVector(X);
    //int Q=size.GetQExtRows(), L1, L2, Lmin=size.GetMinLength(), Lmax=size.GetMaxLength();
    int Q=size.GetQExtRows(), N1, N2, Nmid;
    if(Q%2==0){
        Nmid=Q/2;
    }else{
        Nmid=(Q-1)/2;
    }
    for(N1=1; N1<Nmid; N1++){
        N2=Q-N1+1;
        Array2DSwapExtRows(X, N1, N2);
    }
}
//
template<typename T>Array2DReverseIneRows(T**&X, const Array2DSize& size){
    int Q=size.GetQExtRows(), N1, N2, Nmid, Lmax=size.GetMaxLength(), Lmin=size.GetMinLength();
    if(Lmin==Lmax){
        if(Lmax%2==0){
            Nmid=Q/2;
        }else{
            Nmid=(Lmax-1)/2;
        }
        for(N1=1; N1<Nmid; N1++){
            N2=Lmax-N1+1;
            for(int j=1; j<=Q; j++){
                Array2DSwapVals(X, size, j, N1, j, N2);
            }
        }
    }
}
template<typename T>Array2DReverseIneRows(std::vector<std::vector<T>>&X){
    Array2DSize size;
    size=GetSizeOf2DVector(X);
    int Q=size.GetQExtRows(), N1, N2, Nmid, Lmax=size.GetMaxLength(), Lmin=size.GetMinLength();
    if(Lmin==Lmax){
        if(Lmax%2==0){
            Nmid=Q/2;
        }else{
            Nmid=(Lmax-1)/2;
        }
        for(N1=1; N1<Nmid; N1++){
            N2=Lmax-N1+1;
            for(int j=1; j<=Q; j++){
                Array2DSwapVals(X, size, j, N1, j, N2);
            }
        }
    }
}
//
template<typename T>Array2DReverseExtRowN(T**&X, const Array2DSize& size, int ExtRowN){
    int Lcur;
    if(ExtRowN>=1 && ExtRowN<=size.GetQExtRows()){
        Lcur=size.GetLength(ExtRowN);
        Array1DReverse(X[ExtRowN-1], Lcur);
    }
}
template<typename T>Array2DReverseExtRowN(std::vector<std::vector<T>>&X, int ExtRowN){
    int Lcur;
    if(ExtRowN>=1 && ExtRowN<=X.size()){
        Lcur=X[ExtRowN-1].size();
        Array1DReverse(X[ExtRowN-1]);
    }
}
template<typename T>Array2DReverseIneRowN(T**&X, const Array2DSize& size, int IneRowN){
    int Q=size.GetQExtRows(), Lmin=size.GetMinLength(), Lmax=size.GetMaxLength(), Nmid, N1, N2;
    if(Q>0 && IneRowN>=1 && IneRowN<=Lmin){
        if(Q%2==0){
            Nmid=Q/2;
        }else{
            Nmid=(Q-1)/2;
        }
        for(N1=1; N1<=Nmid; N1++){
            N2=Q-N1+1;
            Array2DSwapVals(X, size, N1, IneRowN, N2, IneRowN);
        }
    }
}
template<typename T>Array2DReverseIneRowN(std::vector<std::vector<T>>&X, int IneRowN){
    Array2DSize size;
    size=GetSizeOf2DVector(X);
    int Q=size.GetQExtRows(), Lmin=size.GetMinLength(), Lmax=size.GetMaxLength(), Nmid, N1, N2;
    if(Q>0 && IneRowN>=1 && IneRowN<=Lmin){
        if(Q%2==0){
            Nmid=Q/2;
        }else{
            Nmid=(Q-1)/2;
        }
        for(N1=1; N1<=Nmid; N1++){
            N2=Q-N1+1;
            Array2DSwapVals(X, N1, IneRowN, N2, IneRowN);
        }
    }
}
//AddToExtRow
template<typename T>Array2DAddToExtRowN(T**&X, Array2DSize& size, int ExtRowN, T val){
    int Q=size.GetQExtRows(), Lmin=size.GetMinLength(), Lmax=size.GetMaxLength(), Lcur;
    if(ExtRowN>=1 && ExtRowN<=Q){
        Lcur=size.GetLength(ExtRowN);
        Array1DAddElement(X[ExtRowN-1], Lcur, val);
        size.SetLength(Lcur, ExtRowN);
    }
}
template<typename T>Array2DAddToExtRowN(std::vector<std::vector<T>>&X, int ExtRowN, T val){
    Array2DSize size;
    size=GetSizeOf2DVector(X);
    int Q=size.GetQExtRows(), Lmin=size.GetMinLength(), Lmax=size.GetMaxLength(), Lcur;
    if(ExtRowN>=1 && ExtRowN<=Q){
        Lcur=size.GetLength(ExtRowN);
        Array1DAddElement(X[ExtRowN-1], val);
        size.SetLength(Lcur, ExtRowN);
    }
}
//InsToExtRow
template<typename T>Array2DInsToExtRowN(T**&X, Array2DSize& size, int ExtRowN, T val, int IneRowPosN){
    int Q=size.GetQExtRows(), Lmin=size.GetMinLength(), Lmax=size.GetMaxLength(), Lcur;
    if(ExtRowN>=1 && ExtRowN<=Q){
        Lcur=size.GetLength(ExtRowN);
        Array1DInsElementToN(X[ExtRowN-1], Lcur, IneRowPosN, val);
        size.SetLength(Lcur, ExtRowN);
    }
}
template<typename T>Array2DInsToExtRowN(std::vector<std::vector<T>>&X, int ExtRowN, T val, int IneRowPosN){
    Array2DSize size;
    size=GetSizeOf2DVector(X);
    int Q=size.GetQExtRows(), Lmin=size.GetMinLength(), Lmax=size.GetMaxLength(), Lcur;
    if(ExtRowN>=1 && ExtRowN<=Q){
        Lcur=size.GetLength(ExtRowN);
        Array1DInsElementToN(X[ExtRowN-1],IneRowPosN, val);
        size.SetLength(Lcur, ExtRowN);
    }
}
//DelFromExtRow
template<typename T>Array2DDelFromExtRowN(T**&X, Array2DSize& size, int ExtRowN, int IneRowPosN){
    int Q=size.GetQExtRows(), Lmin=size.GetMinLength(), Lmax=size.GetMaxLength(), Lcur;
    if(ExtRowN>=1 && ExtRowN<=Q){
        Lcur=size.GetLength(ExtRowN);
        Array1DDelElementFromN(X[ExtRowN-1], Lcur, IneRowPosN);
        size.SetLength(Lcur, ExtRowN);
    }
}
template<typename T>Array2DDelFromExtRowN(std::vector<std::vector<T>>&X, int ExtRowN, int IneRowPosN){
    Array2DSize size;
    size=GetSizeOf2DVector(X);
    int Q=size.GetQExtRows(), Lmin=size.GetMinLength(), Lmax=size.GetMaxLength(), Lcur;
    if(ExtRowN>=1 && ExtRowN<=Q){
        Lcur=size.GetLength(ExtRowN);
        Array1DDelElementFromN(X[ExtRowN-1], IneRowPosN);
        size.SetLength(Lcur, ExtRowN);
    }
}
//
//adding ext rows
template<typename T>void Array2DAddExtRow(T**&X, Array2DSize&size, T*rowParam, int wholeRowL=0, int whatFromN=1, int QDefaultValsBefore=0, bool RectNotVar=true, T*DfltValParam=NULL, TValsShowHide*vsh=NULL){
    std::vector<T> row;
    T DefaultVal, CurVal, *TmpPtr=NULL;
    int Lreq, QTmp;
    Array2DSize NewSize=size;
    if(DfltValParam!=NULL){
        DefaultVal=(*(DfltValParam));
    }
    int vrnN=0;//0-add NULL et Set, 1-add 1-member row et Set,
               //2- add needed row et reset, 3 - 2 ac reset
               //4
    if((size.GetQExtRows()==0 || X==NULL) && wholeRowL>0){//creating new on NULL base
        QTmp=0;
        Lreq=0;
        //std::vector<T> std::vector<T> Array1DGetSubArray_byLs_PreMark(T*x, int Lold, int Lreq=0, int FromN=1, int FirstDefaultValues=0, T*DfltValParam=NULL);
        row=Array1DGetSubArray_byLs_PreMark(rowParam, wholeRowL, Lreq, whatFromN, QDefaultValsBefore, DfltValParam);
        //
        Array1DAddElement(X, QTmp, row.data());
        //
        if(RectNotVar){
            size.SetVaria(1, 1);
        }else{
            size.Set(1, 1);
        }
    }else{//creating
        if(size.isRectangular() && RectNotVar){
            Lreq=size.GetIneRowsLength();
        }else{
            // int CalcLreq_Array1DSubArray_byLs_MarkupNs_vFmls(int Lini, int whatFromN=1, int QDefaultValuesBefore=0, int LreqGiven=0);
            Lreq=CalcLreq_Array1DSubArray_byLs_MarkupNs_vFmls(wholeRowL, whatFromN, QDefaultValsBefore, 0);
        }
        //std::vector<T> Array1DGetSubArray_byLs_PreMark(T*x, int Lold, int Lreq=0, int FromN=1, int FirstDefaultValues=0, T*DfltValParam=NULL)
        row=Array1DGetSubArray_byLs_PreMark(rowParam, wholeRowL, Lreq, whatFromN, QDefaultValsBefore, DfltValParam);
        QTmp=size.GetQExtRows();
        //
        //Array1DAddElement(X, QTmp, row.data());//
        //I men S narb, et ob af fn stop arb, row, ad ic ptr zeig, deda, et ptr adbe noin, ne null, ma zeig ad anid random
        //T*TmpPtr=NULL;//narb, I men ob am mub
        //Array1DAddElement(X, QTmp, NULL);//S nint if ptr = NULL?! Ja, param tam ms'be ne unes NULL, ma ptr of T type, qof val=NULL!
        //std::cout<<"Add ext row: vrnN="<<vrnN<<std::endl;
        switch(vrnN){
            case 0:
                //LetIt remain NULL
            break;
            case 1:
                TmpPtr=new T[1];
                TmpPtr[1-1]=DefaultVal;
            break;
            case 2:
            case 3:
                TmpPtr=new T[Lreq];
                for(int i=1; i<=Lreq; i++){
                    CurVal=row[i-1];
                    TmpPtr[i-1]=CurVal;
                    //TmpPtr[i-1]=row[i-1];
                }
            break;
            case 4:
            case 5:
                //NOp, see further
                //ne tested yet
            break;
            case 6:
                NewSize.AddExt(row.size());
                Array2DSetSize(X, size, NewSize, DfltValParam, true);
                //ne tested yet
            break;
        }
        if(vrnN!=6 && vrnN!=4 && vrnN!=5){
            Array1DAddElement(X, QTmp, TmpPtr);
        }
        if( vrnN==4 || vrnN==5){
             Array1DAddElement(X, QTmp, row.data());
             //ne tested yet
        }
        //it makes QTmp++
        switch(vrnN){
            case 0:
                size.AddExt(0);
            break;
            case 1:
                size.AddExt(1);
            break;
            case 2:
            case 3:
                size.AddExt(Lreq);
            break;
        }
        //;
        //std::cout<<"Temporally arr2D's last row is:"<<std::endl;
        //void Arr1DStdCOut(T*X, int L, QString delim="; ", bool useBrackets=true){
        //Arr1DStdCOut(X[QTmp-1], Lreq, " ", false);//defined later
        //for(int i=1; i<=1; i++)std::cout<<X[QTmp-1][i-1]<<" ";
        //for(int i=1; i<=Lreq; i++)std::cout<<X[QTmp-1][i-1]<<" ";
        //std::cout<<std::endl;
        //
        //Array1DAddElement(X, QTmp,  row.data());
        //template<typename T> void Array2DSetExtRowN(T**&X, Array2DSize&size, int ExtRowN, int wholeRowL=0, T*rowParam=NULL, int whatFromN=1, int QDefaultValuesBefore=0, bool KeepIfRect=true, T*DefaultValParam=NULL, bool DefaultIsSpecNotOwn=false){
        if(vrnN!=3 &&  vrnN!=4 && vrnN!=5){
            Array2DSetExtRowN(X, size, QTmp, Lreq, row.data(), whatFromN, QDefaultValsBefore, RectNotVar, DfltValParam, true);
        }
        //
        //size.AddExt(Lreq);//zu afnot!
    }
}
template<typename T>void Array2DAddExtRow(T**&X, Array2DSize&size, std::vector<T>rowParam, int whatFromN=1, int QDefaultValsBefore=0, bool RectNotVar=true, T*DfltValParam=NULL, TValsShowHide*vsh=NULL){
    Array2DAddExtRow(X, size, rowParam.data(), rowParam.size(), whatFromN, QDefaultValsBefore, RectNotVar, DfltValParam, vsh);
}
template<typename T>void Array2DAddExtRow(std::vector<std::vector<T>>&X, std::vector<T>rowParam,  int whatFromN=1, int QDefaultValsBefore=0, bool RectNotVar=true, T*DfltValParam=NULL, TValsShowHide*vsh=NULL){
    std::vector<T> row;
    T DefaultVal, CurVal;
    int wholeRowL=rowParam.size();
    int Lreq;
    int  preN1_=0, preN2_=0, ownN1_=0, ownN2=0, ForOwnN1_=0, ForOwnN2=0, postN1=0, postN2_=0;//, Lnew=0;
    Array2DSize size;
    size=GetSizeOf2DVector(X);
    if(DfltValParam!=NULL){
        DefaultVal=(*(DfltValParam));
    }
    if(size.GetQExtRows()==0  && wholeRowL>0){
        Lreq=0;
        //template <typename T>std::vector<T> Array1DGetSubArray_byLs_PreMark(std::vector<T>x,  int Lreq=0, int FromN=1, int FirstDefaultValues=0, T*DfltValParam=NULL){
        row=Array1DGetSubArray_byLs_PreMark(rowParam,  Lreq, whatFromN, QDefaultValsBefore, DfltValParam);
        Array1DAddElement(X, row);
    }else{
        if(size.isRectangular() && RectNotVar){
            Lreq=size.GetIneRowsLength();
        }else{
            //void CalcLreq_Array1DSubArray_byLs_MarkupNs_vFmls(int Lini, int whatFromN=1, int QDefaultValuesBefore=0, int LreqGiven=0);
            Lreq=CalcLreq_Array1DSubArray_byLs_MarkupNs_vFmls(rowParam.size(), whatFromN, QDefaultValsBefore, 0);
        }
        //std::vector<T> Array1DGetSubArray_byLs_PreMark(T*x, int Lold, int Lreq=0, int FromN=1, int FirstDefaultValues=0, T*DfltValParam=NULL)
        //row=Array1DGetSubArray_byLs_PreMark(rowParam, wholeRowL, Lreq, whatFromN, QDefaultValsBefore, DfltValParam);
        row=Array1DGetSubArray_byLs_PreMark(rowParam, Lreq, whatFromN, QDefaultValsBefore, DfltValParam);
        //QTmp=size.GetQExtRows();
        Array1DAddElement(X, row);
        //
        //size.AddExt(Lreq);
    }
}
template<typename T>void Array2DAddExtRow(std::vector<std::vector<T>>&X, int wholeRowL=0, T*rowParam=NULL,  int FromN=1, int QDefaultValsBefore=0, bool RectNotVar=true, T*DfltValParam=NULL, TValsShowHide*vsh=NULL){
    std::vector<T>rowIni;
    //template <typename T> std::vector<T> Array1DAssign_VbyP(T*y=NULL, int Qy=0,  T*DfltValParam=NULL){
    rowIni=Array1DAssign_VbyP(wholeRowL, rowParam, DfltValParam);
    Array2DAddExtRow(X, rowIni, FromN, QDefaultValsBefore, RectNotVar, DfltValParam, vsh);
}
//

template<typename T>void Array2DAddExtRow_ifDfltNeeded_v1(T**&X, Array2DSize&size, T dfltVal, T*rowParam,  int wholeRowL=0, int whatFromN=1, int QDefaultValsBefore=0, bool RectIfPossible=true, TValsShowHide*vsh=NULL){
    T*y=NULL;
    std::vector<T>rowTmp;
    int Q=size.GetQExtRows();
    if(size.isRectangular() && wholeRowL==0){
        wholeRowL=size.GetLength(1);
    }
    if(wholeRowL>0){
        if(rowParam!=NULL){
            rowTmp=Array1DGetSubArray_byLs_PreMark_v1(rowParam, wholeRowL, dfltVal, wholeRowL, whatFromN, QDefaultValsBefore);
        }else{
            for(int i=1; i<=wholeRowL; i++){
                rowTmp.push_back(dfltVal);
             }
        }
        y=new T[wholeRowL];
        for(int i=1; i<=wholeRowL; i++){
            y[i-1]=rowTmp[i-1];
        }
        Array1DAddElement(X, Q, y);//Q=Q+1
        size.AddExt(wholeRowL);
        if(size.GetMinLength()==size.GetMaxLength() && RectIfPossible){
            size.Set(Q, wholeRowL);
        }
    }
}
template<typename T>void Array2DAddExtRow_ifDfltNeeded_v2(T**&X, Array2DSize&size, T DfltValParam, T*rowParam,  int wholeRowL=0, int whatFromN=1, int QDefaultValsBefore=0, bool RectIfPossible=true, TValsShowHide*vsh=NULL){
    T*y=NULL;
    y=new T;
    (*(y))=DfltValParam;
    //template<typename T>void Array2DAddExtRow(T**&X, Array2DSize&size, T*rowParam, int wholeRowL=0, int whatFromN=1, int QDefaultValsBefore=0, bool RectNotVar=true, T*DfltValParam=NULL, TValsShowHide*vsh=NULL){
    Arra2DAddExtRow(                                X,             size,   rowParam,     wholeRowL,       whatFromN,       QDefaultValsBefore,        RectIfPossible,    DfltValParam,                     vsh);
    delete y;
}
template<typename T>void Array2DAddExtRow_ifDfltNeeded_v3(T**&X, Array2DSize&size, T dfltVal, T*rowParam,  int wholeRowL=0, int whatFromN=1, int QDefaultValsBefore=0, int SizeSetMode_Varia0_Rect_ifPossible1_ByArr2ByAddedRow3ByMax4=1, TValsShowHide*vsh=NULL){
    T*y=NULL;
    int Lreq, minOfL, maxOfL, Lmax=size.GetMaxLength(), Q=size.GetQExtRows();
    std::vector<T>rowTmp;
    Array2DSize size2;
    if(size.isRectangular() && rowParam!=NULL && wholeRowL==0){
        wholeRowL=size.GetLength(1);
    }
    minOfL= wholeRowL <=Lmax ? wholeRowL : Lmax;
    maxOfL= wholeRowL >=Lmax ? wholeRowL : Lmax;
    switch(SizeSetMode_Varia0_Rect_ifPossible1_ByArr2ByAddedRow3ByMax4){
        case 0:
            Lreq=wholeRowL;
            //template<typename T>void Array2DAddExtRow_ifDfltNeeded_v1(T**&X, Array2DSize&size, T dfltVal, T*rowParam,  int wholeRowL=0, int whatFromN=1, int QDefaultValsBefore=0, bool RectIfPossible=true, TValsShowHide*vsh=NULL){
            Array2DAddExtRow_ifDfltNeeded_v1(                               X,             size,   dfltVal,   rowParam,      wholeRowL,       whatFromN,       QDefaultValsBefore,                       true,               vsh);
        break;
        case 1:
            Lreq=wholeRowL;
            //template<typename T>void Array2DAddExtRow_ifDfltNeeded_v1(T**&X, Array2DSize&size, T dfltVal, T*rowParam,  int wholeRowL=0, int whatFromN=1, int QDefaultValsBefore=0, bool RectIfPossible=true, TValsShowHide*vsh=NULL){
            Array2DAddExtRow_ifDfltNeeded_v1(                               X,             size,   dfltVal,   rowParam,      wholeRowL,       whatFromN,       QDefaultValsBefore,                       true,               vsh);
        break;
        case 2:

        break;
        case 3:

        break;
        case 4:

        break;
    }

    if(wholeRowL>0){
        if(rowParam!=NULL){
            rowTmp=Array1DGetSubArray_byLs_PreMark_v1(rowParam, wholeRowL, dfltVal, wholeRowL, whatFromN, QDefaultValsBefore);
        }else{
            for(int i=1; i<=wholeRowL; i++){
                rowTmp.push_back(dfltVal);
             }
        }
        y=new T[wholeRowL];
        for(int i=1; i<=wholeRowL; i++){
            y[i-1]=rowTmp[i-1];
        }
        Array1DAddElement(X, Q, y);//Q=Q+1
        size.AddExt(wholeRowL);
        if((size.GetMinLength()==size.GetMaxLength() && SizeSetMode_Varia0_Rect_ifPossible1_ByArr2ByAddedRow3ByMax4!=0==1) || SizeSetMode_Varia0_Rect_ifPossible1_ByArr2ByAddedRow3ByMax4>1){
            size.Set(Q, wholeRowL);
        }
    }
}

//
//Del ext row
template<typename T>void Array2DDelExtRow(T**&X, Array2DSize&size, int N){//42
    int Q=size.GetQExtRows(), L;
    if(N>=1 && N<=Q){
        Array1DDelElementFromN(X, Q, N);
        size.DelExt(N);
    }
}
template<typename T>void Array2DDelExtRow(std::vector<std::vector<T>>&X, int N){
    if(N>=1 && N<=X.size()){
        Array1DDelElementFromN(X, N);//if cany overload
        //X.erase(x.begin()+N-1);
    }
}
//
//Ins ext row
template<typename T> void Array2DInsExtRow(T**&X,  Array2DSize&size, int ExtRowN, int wholeRowL=0, T*rowParam=NULL, T*DefaultValParam=NULL, int whatFromN=1, int QDefaultValuesBefore=0, bool KeepIfRect=true, TValsShowHide*vsh=NULL){
    //writeln(&vsh,"Array2DInsExtRow(ptr, ptr) starts working");
    if(vsh!=NULL)std::cout<<" Array2DInsExtRow(ptr, ptr) starts working "<<std::endl;
    int Q=size.GetQExtRows();
    if(ExtRowN<0){
        ExtRowN=Q+ExtRowN+1;
    }
    if(ExtRowN>=1 && ExtRowN<=Q){
        //void Array2DAddExtRow(T**&X, Array2DSize&size, T*rowParam, int wholeRowL=0, int whatFromN=1, int QDefaultValsBefore=0, bool RectNotVar=true, T*DfltValParam=NULL, TValsShowHide*vsh=NULL)
        Array2DAddExtRow(           X,             size,   rowParam,     wholeRowL,       whatFromN,       QDefaultValuesBefore,      KeepIfRect,        DefaultValParam);
        //ob i add is done all complex work, et for Q=0 S narb
        Q=size.GetQExtRows();
        if(vsh!=NULL){
            std::cout<<" Seems to be added "<<std::endl;
        }
        if(vsh!=NULL){
            //void Arr2DStdCOut2D(T**X, const Array2DSize size, QString delimElem="; ", QString delimLines="; ", bool useBrackets=true){
            Arr2DStdCOut2D(X, size);
        }
        for(int i=Q-1; i>=ExtRowN; i--){
            Array2DSwapExtRows(X, size, i, i+1);
            if(vsh!=NULL){
                std::cout<<"Swapping rows "<<i<<" and "<<i+1<<std::endl;
                Arr2DStdCOut2D(X, size);
            }
        }
    }
    //writeln(&vsh,"Array2DInsExtRow(ptr, ptr) finishes working");
    if(vsh!=NULL)std::cout<<" Array2DInsExtRow(ptr, ptr) finishes working "<<std::endl;
}
template<typename T> void Array2DInsExtRow(T**&X,  Array2DSize&size, int ExtRowN, std::vector<T>rowParam, T*DefaultValParam=NULL, int whatFromN=1, int QDefaultValuesBefore=0, bool KeepIfRect=true){
    Array2DInsExtRow(X, size, ExtRowN, rowParam.size(), rowParam.data(), DefaultValParam, whatFromN, QDefaultValuesBefore, KeepIfRect);
}
template<typename T> void Array2DInsExtRow(std::vector<std::vector<T>>&X, int ExtRowN, std::vector<T>rowParam, T*DefaultValParam=NULL, int whatFromN=1, int QDefaultValuesBefore=0, bool KeepIfRect=true){
    Array2DSize size;
    size=GetSizeOf2DVector(X);
    //newSize=size;
    std::vector<T>row;
    int Q=size.GetQExtRows(), Lmin=size.GetMinLength(), Lmax=size.GetMaxLength(), Lreq=0;
    if(ExtRowN<0){
        ExtRowN=Q+ExtRowN+1;
    }
    if(ExtRowN>=1 && ExtRowN<=Q){
        if(Lmin==Lmax && KeepIfRect){
            Lreq=Lmin;//or Lmax, D s'egal
        }
        row=Array1DGetSubArray_byLs_PreMark(rowParam, Lreq, whatFromN, QDefaultValuesBefore, DefaultValParam);
        Array1DInsElementToN(X, ExtRowN, row);
        //ob co vects all s'unes
    }
}
template<typename T> void Array2DInsExtRow(std::vector<std::vector<T>>&X, int ExtRowN, int wholeRowL=0, T*rowParam=NULL, T*DefaultValParam=NULL, int whatFromN=1, int QDefaultValuesBefore=0, bool KeepIfRect=true){
    Array2DSize size;//, newSize;
    size=GetSizeOf2DVector(X);
    //newSize=size;
    std::vector<T>row;
    int Q=size.GetQExtRows(), Lmin=size.GetMinLength(), Lmax=size.GetMaxLength(), Lreq=0;
    if(ExtRowN<0){
        ExtRowN=Q+ExtRowN+1;
    }
    if(ExtRowN>=1 && ExtRowN<=Q){
        if(Lmin==Lmax && KeepIfRect){
            Lreq=Lmin;//or Lmax, D s'egal
        }
        row=Array1DGetSubArray_byLs_PreMark(rowParam, wholeRowL, Lreq, whatFromN, QDefaultValuesBefore, DefaultValParam);
        Array1DInsElementToN(X, ExtRowN, row);
        //ob co vects all s'unes
    }
}
//
//SetIneRow
template<typename T>Array2DSetIneRow(T**&X, Array2DSize& size, int IneRowN, T*rowParam=NULL, int wholeRowL=0, int whatFromN=1, int QDefaultValsBefore=0, T*DefaultValParam=NULL, bool ExistingInsteadOfDefault=true, bool ExistingOnlyForGTLminNotIgnore=false, bool LastPossIfNotRectAtPos0NotIgnore=false){
    int Q=size.GetQExtRows(), Lmin=size.GetMinLength(), Lmax=size.GetMaxLength(), Lcur, n, Lreq=Q;
    std::vector<T>row;//=Array1DGetSubArray_byLs_PreMark(rowParam, wholeRowL, Lreq, whatFromN, QDefaultValsBefore, DefaultValParam);
    T DefaultVal;
    int LreqCalcd=0, preN1_=0, preN2_=0, ownN1_=0, ownN2=0, ForOwnN1_=0, ForOwnN2=0, postN1=0, postN2_=0;
    //param Ns
    if(DefaultValParam!=NULL){
        DefaultVal=(*(DefaultValParam));
    }
    if(wholeRowL==0 && rowParam!=NULL){
        wholeRowL=Lreq;
    }
    if(IneRowN<0){
        IneRowN=Lmax+IneRowN+1;
    }
    //row
    //void Calc_Array1DSubArray_byLs_MarkupNs_vFmls(int&LreqCalcd, int&preN1_, int&preN2_, int&ownN1_, int&ownN2, int&ForOwnN1_, int&ForOwnN2, int&postN1, int&postN2_, int Lini, int LreqGiven=0, int whatFromN=1, int QDefaultValuesBefore=0);
    Calc_Array1DSubArray_byLs_MarkupNs_vFmls(LreqCalcd, preN1_, preN2_, ownN1_, ownN2, ForOwnN1_, ForOwnN2, postN1, postN2_, wholeRowL, Lreq, whatFromN, QDefaultValsBefore);
    row=Array1DGetSubArray_byLs_PreMark(rowParam, wholeRowL, Lreq, whatFromN, QDefaultValsBefore, DefaultValParam);
    if(ExistingInsteadOfDefault){
        //if(preN1_>0){//arb ac ce ob i=0 et length of 0. row =0, IneRowN>Lcur
            for(int i=preN1_; i<=preN2_; i++){
                Lcur=size.GetLength(i-1);
                if(IneRowN>=1 && IneRowN<=Lcur){
                    row[i-1]=X[i-1][IneRowN-1];
                }
            }
        //}
        //if(postN1>0){//arb ac ce ob i=0 et length of 0. row =0, IneRowN>Lcur
            for(int i=postN1; i<=postN2_; i++){
                Lcur=size.GetLength(i-1);
                if(IneRowN>=1 && IneRowN<=Lcur){
                    row[i-1]=X[i-1][IneRowN-1];
                }
            }
       //}
    }
    //set ipse
    if(IneRowN>=1 && IneRowN<=Lmax){
        if(IneRowN<=Lmin){
            for(int i=1; i<=Q; i++){
                X[i-1][IneRowN-1]=row[i-1];
            }
        }else if(ExistingOnlyForGTLminNotIgnore){
            for(int i=1; i<=Q; i++){
                Lcur=size.GetLength(i);
                if(IneRowN<=Lcur)
                X[i-1][IneRowN-1]=row[i-1];
            }
        }
    }else if(IneRowN==0 && LastPossIfNotRectAtPos0NotIgnore && Lmin!=Lmax){
        for(int i=1; i<=Q; i++){
            Lcur=size.GetLength(i);
            X[i-1][Lcur-1]=row[i-1];
        }
    }
}
template<typename T>Array2DSetIneRow(T**&X, Array2DSize& size, int IneRowN,std::vector<T>rowParam, int whatFromN=1, int QDefaultValsBefore=0, T*DefaultValParam=NULL, bool ExistingInsteadOfDefault=true, bool ExistingOnlyForGTLminNotIgnore=false, bool LastPossIfNotRectAtPos0NotIgnore=false){
    Array2DSetIneRow(X, size, IneRowN, rowParam.data(), rowParam.size(), whatFromN, QDefaultValsBefore, DefaultValParam, ExistingInsteadOfDefault, ExistingOnlyForGTLminNotIgnore, LastPossIfNotRectAtPos0NotIgnore);
}
template<typename T>Array2DSetIneRow(std::vector<std::vector<T>>&X, int IneRowN, std::vector<T>rowParam, int whatFromN=1, int QDefaultValsBefore=0, T*DefaultValParam=NULL, bool ExistingInsteadOfDefault=true, bool ExistingOnlyForGTLminNotIgnore=false, bool LastPossIfNotRectAtPos0NotIgnore=false){
    Array2DSize size;
    size=GetSizeOf2DVector(X);
    int wholeRowL=rowParam.size();
    int Q=size.GetQExtRows(), Lmin=size.GetMinLength(), Lmax=size.GetMaxLength(), Lcur, n, Lreq=Q;
    std::vector<T>row;//=Array1DGetSubArray_byLs_PreMark(rowParam, wholeRowL, Lreq, whatFromN, QDefaultValsBefore, DefaultValParam);
    T DefaultVal;
    int LreqCalcd=0, preN1_=0, preN2_=0, ownN1_=0, ownN2=0, ForOwnN1_=0, ForOwnN2=0, postN1=0, postN2_=0;
    //param Ns
    if(DefaultValParam!=NULL){
        DefaultVal=(*(DefaultValParam));
    }
    //if(wholeRowL==0 && rowParam!=NULL){
    //   wholeRowL=Lreq;
    //}
    if(IneRowN<0){
        IneRowN=Lmax+IneRowN+1;
    }
    //row
    //void Calc_Array1DSubArray_byLs_MarkupNs_vFmls(int&LreqCalcd, int&preN1_, int&preN2_, int&ownN1_, int&ownN2, int&ForOwnN1_, int&ForOwnN2, int&postN1, int&postN2_, int Lini, int LreqGiven=0, int whatFromN=1, int QDefaultValuesBefore=0);
    Calc_Array1DSubArray_byLs_MarkupNs_vFmls(LreqCalcd, preN1_, preN2_, ownN1_, ownN2, ForOwnN1_, ForOwnN2, postN1, postN2_, wholeRowL, Lreq, whatFromN, QDefaultValsBefore);
    row=Array1DGetSubArray_byLs_PreMark(rowParam, Lreq, whatFromN, QDefaultValsBefore, DefaultValParam);
    if(ExistingInsteadOfDefault){
        if(preN2_>0){
            for(int i=preN1_; i<=preN2_; i++){
                Lcur=size.GetLength(i-1);
                if(IneRowN>=1 && IneRowN<=Lcur){
                    row[i-1]=X[i-1][IneRowN-1];
                }
            }
        }
        if(postN1>0){
            for(int i=postN1; i<=postN2_; i++){
                Lcur=size.GetLength(i-1);
                if(IneRowN>=1 && IneRowN<=Lcur){
                    row[i-1]=X[i-1][IneRowN-1];
                }
            }
        }
    }
    //set ipse
    if(IneRowN>=1 && IneRowN<=Lmax){
        if(IneRowN<=Lmin){
            for(int i=1; i<=Q; i++){
                X[i-1][IneRowN-1]=row[i-1];
            }
        }else if(ExistingOnlyForGTLminNotIgnore){
            for(int i=1; i<=Q; i++){
                Lcur=size.GetLength(i);
                if(IneRowN<=Lcur)
                X[i-1][IneRowN-1]=row[i-1];
            }
        }
    }else if(IneRowN==0 && LastPossIfNotRectAtPos0NotIgnore && Lmin!=Lmax){
        for(int i=1; i<=Q; i++){
            Lcur=size.GetLength(i);
            X[i-1][Lcur-1]=row[i-1];
        }
    }
}
template<typename T>Array2DSetIneRow(std::vector<std::vector<T>>&X, int IneRowN, T*rowParam, int wholeRowL=0, int whatFromN=1, int QDefaultValsBefore=0, T*DefaultValParam=NULL, bool ExistingInsteadOfDefault=true, bool ExistingOnlyForGTLminNotIgnore=false, bool LastPossIfNotRectAtPos0NotIgnore=false){
    int Q=X.size();
    std::vector<T>row;
    //template <typename T> void Array1DAssign(std::vector<T>&vec, T y[], int L=0)
    Array1DAssign(row, rowParam, wholeRowL);
    //template<typename T>Array2DSetIneRow(std::vector<std::vector<T>>&X, int IneRowN, std::vector<T>rowParam, int whatFromN=1, int QDefaultValsBefore=0, T*DefaultValParam=NULL, bool ExistingInsteadOfDefault=true, bool ExistingOnlyForGTLminNotIgnore=false, bool LastPossIfNotRectAtPos0NotIgnore=false){
    Array2DSetIneRow(                                                  X,     IneRowN,                    row,     whatFromN,        QDefaultValsBefore,    DefaultValParam,           ExistingInsteadOfDefault,           ExistingOnlyForGTLminNotIgnore,            LastPossIfNotRectAtPos0NotIgnore);
}
//
//AddIneRow
template<typename T>Array2DAddIneRow(T**&X, Array2DSize& size, T*rowParam=NULL, int wholeRowL=0, int IfNonRectIgnore0Add1Stretch2=0, int whatFromN=1, int QDefaultValsBefore=0, T*DefaultValParam=NULL, bool RectNotVariaIfFirst=true){
    std::vector<T>row;
    T DefaultVal;
    int Q=size.GetQExtRows(), Lmin=size.GetMinLength(), Lmax=size.GetMaxLength(), Lcur, n;
    int Lreq;
    int *Ls=NULL;
    if(DefaultValParam!=NULL){
        DefaultVal=(*(DefaultValParam));
    }
    if(wholeRowL==0 && rowParam!=NULL){
        wholeRowL=Q;
    }
    if(Q==0 && wholeRowL>0){
        Lreq=0;
        row=Array1DGetSubArray_byLs_PreMark(rowParam, wholeRowL, Lreq, whatFromN, QDefaultValsBefore, DefaultValParam);
        Lcur=row.size();//in other if brunch Lcur means different thing! n not used
        //ob not from first et co Defaults before
        X=new T*[Lcur];
        for(int i=1; i<=Lcur; i++){
            X[i-1]=new T[1];
        }
        if(rowParam!=NULL){
            for(int i=1; i<=Lcur; i++){
                X[i-1][1-1]=row[i-1];
            }
        }else{
            for(int i=1; i<=Lcur; i++){
                X[i-1][1-1]=DefaultVal;
            }
        }
        //size:
        if(RectNotVariaIfFirst){
            size.Set(Lcur, 1);
        }else{
            Ls=new int[Lcur];
            for(int i=1; i<=Lcur; i++){
                Ls[i-1]=1;
            }
            size.Set(Lcur, Ls);
            delete[]Ls;
            Ls=NULL;
        }
    }else if(Q>0){
        Lmin=size.GetMinLength();
        Lmax=size.GetMaxLength();
        Lreq=Q;
        //row=Array1DGetSubArray_byLs_PreMark(rowParam, wholeRowL, Q, whatFromN, QDefaultValuesBefore, DefaultValParam);
        //ob ine row has length Q
        if(Lmin!=Lmax && IfNonRectIgnore0Add1Stretch2==2){
            for(int i=1; i<=Q; i++){
                Lcur=size.GetLength(i);
                n=Lcur;
                for(int j=Lcur+1; j<=Lmax; j++){
                    Array1DAddElement(X[i-1], n, DefaultVal);
                }
            }
            for(int i=1; i<=Q; i++){
                size.SetLength(Lmax, i);
            }
        }
        if(!((Lmin!=Lmax)&&(IfNonRectIgnore0Add1Stretch2==0))){
            if(rowParam!=NULL){
                row=Array1DGetSubArray_byLs_PreMark(rowParam, wholeRowL, Lreq, whatFromN, QDefaultValsBefore, DefaultValParam);
                for(int i=1; i<=Q; i++){
                    Lcur=size.GetLength(i);
                    n=Lcur;
                    Array1DAddElement(X[i-1], n, row[i-1]);
                }
            }else{
                for(int i=1; i<=Q; i++){
                    Lcur=size.GetLength(i);
                    n=Lcur;
                    Array1DAddElement(X[i-1], n, DefaultVal);
                }
            }
            if(Lmin!=Lmax || (Lmin==Lmax && size.GetIneRowsLength()==0)){
                for(int i=1; i<=Q; i++){
                    Lcur=size.GetLength(i);
                    size.SetLength(Lcur+1, i);
                }
            }else{
                size.SetIneRowsLength(Lmax+1);
            }
        }
    }
}//add ine row
template<typename T>Array2DAddIneRow(T**&X, Array2DSize& size, std::vector<T>rowParam, int IfNonRectIgnore0Add1Stretch2=0, int whatFromN=1, int QDefaultValsBefore=0, T*DefaultValParam=NULL){
    Array2DAddIneRow(X, size, rowParam.data(), rowParam.size(),  IfNonRectIgnore0Add1Stretch2, whatFromN, QDefaultValsBefore, DefaultValParam);
}
template<typename T>Array2DAddIneRow(std::vector<std::vector<T>>&X, std::vector<T>rowParam, int IfNonRectIgnore0Add1Stretch2=0, int whatFromN=1, int QDefaultValsBefore=0, T*DefaultValParam=NULL){
    int Q=X.size(), wholeRowL=rowParam.size();
    std::vector<T>rowM, rowCur;
    Array2DSize size;
    size=GetSizeOf2DVector(X);
    T curVal, DefaultVal;
    if(DefaultValParam!=NULL){
        DefaultVal=(*(DefaultValParam));
    }
    int Lmin=size.GetMinLength(), Lmax=size.GetMaxLength(), Lreq=0, Lcur, n;
    if(Q==0 && wholeRowL>0){
        Lreq=0;
        rowM=Array1DGetSubArray_byLs_PreMark(rowParam, Lreq, whatFromN, QDefaultValsBefore, DefaultValParam);
        Lcur=rowM.size();//in other if brunch Lcur means different thing! n not used
        //ob not from first et co Defaults before
        for(int i=1; i<=Lcur; i++){
            rowCur.clear();
            curVal=rowM[i-1];
            rowCur.push_back(curVal);
            X.push_back(rowCur);
        }
    }else{
        if(!((Lmin!=Lmax)&&(IfNonRectIgnore0Add1Stretch2==0))){
            Lreq=Q; //ob os ine row
            rowM=Array1DGetSubArray_byLs_PreMark(rowParam, Lreq, whatFromN, QDefaultValsBefore, DefaultValParam);
            if(Lmin!=Lmax && IfNonRectIgnore0Add1Stretch2==2){
                for(int i=1; i<=Q; i++){
                    Lcur=X[i-1].size();
                    //n=Lcur;
                    for(int j=Lcur+1; j<=Lmax; j++){
                        Array1DAddElement(X[i-1], DefaultVal);
                    }
                }
                for(int i=1; i<=Q; i++){
                    size.SetLength(Lmax, i);
                }
            }
            if(wholeRowL>0){
                for(int i=1; i<=Q; i++){
                    Lcur=size.GetLength(i);
                    n=Lcur;
                    Array1DAddElement(X[i-1], rowM[i-1]);
                }
            }else{
                for(int i=1; i<=Q; i++){
                    Lcur=size.GetLength(i);
                    n=Lcur;
                    Array1DAddElement(X[i-1], DefaultVal);
                }
            }
            if(Lmin==Lmax && size.GetIneRowsLength()>0){
                size.SetIneRowsLength(Lmax+1);
            }else{
                for(int i=1; i<=Q; i++){
                     Lcur=size.GetLength(i);
                     size.SetLength(Lcur+1, i);
                }
            }
        }
    }
}//add ine row
template<typename T>Array2DAddIneRow(std::vector<std::vector<T>>&X, T*rowParam=NULL, int wholeRowL=0, int IfNonRectIgnore0Add1Stretch2=0, int whatFromN=1, int QDefaultValsBefore=0, T*DefaultValParam=NULL){
    int Q=X.size(), Lreq=0;
    std::vector<T>row;
    if(wholeRowL==0 && rowParam!=NULL){
        wholeRowL=Q;
    }
    Lreq=Q;//and 0 if 0
    row=Array1DGetSubArray_byLs_PreMark(rowParam, wholeRowL, Lreq, whatFromN, QDefaultValsBefore, DefaultValParam);
    //template<typename T>Array2DAddIneRow(std::vector<std::vector<T>>&X, std::vector<T>rowParam, int IfNonRectIgnore0Add1Stretch2=0, int whatFromN=1, int QDefaultValsBefore=0, T*DefaultValParam=NULL){
    Array2DAddIneRow(X, row, IfNonRectIgnore0Add1Stretch2, whatFromN, QDefaultValsBefore, DefaultValParam);
}////add ine row
//Ins ine row
template<typename T>Array2DInsIneRow(T**&X, Array2DSize& size, int IneRowPosN, T*rowParam=NULL, int wholeRowL=0, int IfAfterLminIgnore0Stretch2=0, int whatFromN=1, int QDefaultValsBefore=0, T*DefaultValParam=NULL){
    int Q=size.GetQExtRows(), Lmin=size.GetMinLength(), Lmax=size.GetMaxLength(), Lcur, n, Lreq;
    std::vector<T>row;
    T DefaultVal;
    if(DefaultValParam!=NULL){
        DefaultVal=(*(DefaultValParam));
    }
    if(rowParam!=NULL && wholeRowL==0){
        wholeRowL=Q;
    }
    if(Q>0 && IneRowPosN>=1 && IneRowPosN<=Lmax && !(IneRowPosN>Lmin && IfAfterLminIgnore0Stretch2==0)){
        if(IneRowPosN>Lmin && IfAfterLminIgnore0Stretch2==2){
            for(int i=1; i<=Q; i++){
                Lcur=size.GetLength(i);
                n=Lcur;
                for(int j=Lcur+1; j<=Lmax; j++){
                    Array1DAddElement(X[i-1], n, DefaultVal);
                }
            }
        }
        Lreq=Q;//and 0 if 0
        row=Array1DGetSubArray_byLs_PreMark(rowParam, wholeRowL, Lreq, whatFromN, QDefaultValsBefore, DefaultValParam);
        if(rowParam!=NULL){
            for(int i=1; i<=Q; i++){
                Array2DInsToExtRowN(X, size, i, row[i-1], IneRowPosN);
            }
        }else{
            for(int i=1; i<=Q; i++){
                Array2DInsToExtRowN(X, size, i, DefaultVal, IneRowPosN);
            }
        }
    }
}
template<typename T>Array2DInsIneRow(T**&X, Array2DSize& size, int IneRowPosN, std::vector<T>rowParam, int IfAfterLminIgnore0Stretch2=0, int whatFromN=1, int QDefaultValsBefore=0, T*DefaultValParam=NULL){
    Array2DInsIneRow(X, size, IneRowPosN, rowParam.data(), rowParam.size(), IfAfterLminIgnore0Stretch2, whatFromN, QDefaultValsBefore, DefaultValParam);
}
template<typename T>Array2DInsIneRow(std::vector<std::vector<T>>&X, int IneRowPosN, std::vector<T>rowParam, int IfAfterLminIgnore0Stretch2=0, int whatFromN=1, int QDefaultValsBefore=0, T*DefaultValParam=NULL){
    std::vector<T>row;
    T DefaultVal;
    Array2DSize size;
    size=GetSizeOf2DVector(X);
    int Q=size.GetQExtRows(), Lmin=size.GetMinLength(), Lmax=size.GetMaxLength(), Lcur, n, Lreq;
    if(DefaultValParam!=NULL){
        DefaultVal=(*(DefaultValParam));
    }
    if(Q>0 && IneRowPosN>=1 && IneRowPosN<=Lmax && !(IneRowPosN>Lmin && IfAfterLminIgnore0Stretch2==0)){
        if(IneRowPosN>Lmin && IfAfterLminIgnore0Stretch2==2){
            for(int i=1; i<=Q; i++){
                Lcur=size.GetLength(i);
                //n=Lcur;
                for(int j=Lcur+1; j<=Lmax; j++){
                    Array1DAddElement(X[i-1],  DefaultVal);
                }
            }
        }
        Lreq=Q;//and 0 if 0
        row=Array1DGetSubArray_byLs_PreMark(rowParam, Lreq, whatFromN, QDefaultValsBefore, DefaultValParam);
        if(rowParam.size()!=0){
            for(int i=1; i<=Q; i++){
                Array2DInsToExtRowN(X, i, row[i-1], IneRowPosN);
            }
        }else{
            for(int i=1; i<=Q; i++){
                Array2DInsToExtRowN(X, i, DefaultVal, IneRowPosN);
            }
        }
    }
}
template<typename T>Array2DInsIneRow(std::vector<std::vector<T>>&X, int IneRowPosN, T*rowParam=NULL, int wholeRowL=0,  int IfAfterLminIgnore0Stretch2=0, int whatFromN=1, int QDefaultValsBefore=0, T*DefaultValParam=NULL){
    std::vector<T>row;
    T DefaultVal;
    Array2DSize size;
    size=GetSizeOf2DVector(X);
    int Q=size.GetQExtRows(), Lmin=size.GetMinLength(), Lmax=size.GetMaxLength(), Lcur, n, Lreq;
    //if(wholeRowL==0)wholeRowL=Lmin;
    Lreq=Q;//and 0 if 0
    //template <typename T>std::vector<T> Array1DGetSubArray_byLs_PreMark(std::vector<T>x, int Lreq=0, int whatFromN=1, int FirstDefaultValues=0, T*DfltValParam=NULL){
    //template <typename T>std::vector<T> Array1DGetSubArray_byLs_PreMark(T*x, int Lold, int Lreq=0, int whatFromN=1, int FirstDefaultValues=0, T*DfltValParam=NULL){
    row=Array1DGetSubArray_byLs_PreMark(rowParam, wholeRowL, Lreq, whatFromN, QDefaultValsBefore, DefaultValParam);
    //template<typename T>Array2DInsIneRow(std::vector<std::vector<T>>&X, int IneRowPosN, std::vector<T>rowParam, int IfAfterLminIgnore0Stretch2=0, int whatFromN=1, int QDefaultValsBefore=0, T*DefaultValParam=NULL){
    Array2DInsIneRow(X, IneRowPosN, row, IfAfterLminIgnore0Stretch2, whatFromN, QDefaultValsBefore, DefaultValParam);
}
//del ine row
template<typename T>Array2DDelIneRow(T**&X, Array2DSize& size, int IneRowPosN, bool ifPos0IgnoreNotDelLastOfVariaNotIgnore=false){
    int Q=size.GetQExtRows(), Lmin=size.GetMinLength(), Lmax=size.GetMaxLength(), Lcur, n;
    if(IneRowPosN>=1 && IneRowPosN<=Lmax){
        for(int i=1; i<=Q; i++){
            Array2DDelFromExtRowN(X, size, i, IneRowPosN);
        }
    }else if(IneRowPosN==0 && ifPos0IgnoreNotDelLastOfVariaNotIgnore==true && Lmin!=Lmax){
        for(int i=1; i<=Q; i++){
            Lcur=size.GetLength(i);
            n=Lcur;
            Array1DDelElementFromN(X[i-1], n, Lcur);
            size.SetLength(n, i);
        }
    }
}
template<typename T>Array2DDelIneRow(std::vector<std::vector<T>>&X, int IneRowPosN, bool ifPos0IgnoreNotDelLastOfVariaNotIgnore=false){
    Array2DSize size;
    size=GetSizeOf2DVector(X);
    int Q=size.GetQExtRows(), Lmin=size.GetMinLength(), Lmax=size.GetMaxLength(), Lcur;
    if(IneRowPosN>=1 && IneRowPosN<=Lmax){
        for(int i=1; i<=Q; i++){
            Array2DDelFromExtRowN(X, i, IneRowPosN);
        }
    }else if(IneRowPosN==0 && ifPos0IgnoreNotDelLastOfVariaNotIgnore==true && Lmin!=Lmax){
        for(int i=1; i<=Q; i++){
            Lcur=size.GetLength(i);
            Array1DDelElementFromN(X[i-1], Lcur);
        }
    }
}
//Assign
template<typename T>void Array2DAssign(T**&XTo, Array2DSize& sizeTo, T**XFrom, const Array2DSize& sizeFrom){
    //void Array2DSetSize(T**&X, const Array2DSize&size1, const Array2DSize&size2, T*DfltValParam=NULL, bool PreserveVals=true)
    T*DfltValParam=NULL;
    int Q=sizeFrom.GetQExtRows(), Lcur;
    Array2DSetSize(XTo, sizeTo, sizeFrom, DfltValParam, false);
    for(int i=1; i<=Q; i++){
       Lcur=sizeFrom.GetLength(i);
       for(int j=1; j<=Lcur; j++){
          XTo[i-1][j-1]=XFrom[i-1][j-1];
       }
    }
}
template<typename T>void Array2DAssign(T**&X, Array2DSize& size, std::vector<std::vector<T>>Z){
    int Q=size.GetQExtRows(), Lcur;
    //Array2DSize  size1=GetSizeOf2DVector(X);
    T val;
    if(X!=NULL){
        for(int i=1; i<=Q; i++){
            delete[]X[i-1];
        }
        delete[]X;
        size=GetSizeOf2DVector(Z);
        Q=size.GetQExtRows();
    }
    if(X.size()!=0){
        X=new T*[Q];
        for(int i=1; i<=Q; i++){
            Lcur=size.GetLength(i);
            X[i-1]=new T[Lcur];
        }
         for(int i=1; i<=Q; i++){
             Lcur=size.GetLength(i);
             for(int j=1; j<=Lcur; j++){
                 val=Z[i-1][j-1];
                 X[i-1][j-1]=val;
             }
         }
    }
}
template<typename T>void Array2DAssign(std::vector<std::vector<T>>&X, Array2DSize& size, T**Z=NULL, T*DefaultValParam=NULL){
    int Q=size.GetQExtRows(), Lcur;
    //Array2DSize  size1=GetSizeOf2DVector(X);
    T curVal, DefaultVal;
    X.clear();
    if(DefaultValParam!=NULL){
        DefaultVal=(*(DefaultValParam));
    }
    std::vector<T>row;
    if(Z!=NULL){
        for(int i=1; i<=Q; i++){
            Lcur=size.GetLength(i);
            row.clear();
            for(int j=1; j<=Lcur; j++){
                curVal=Z[i-1][j-1];
                row.push_back(curVal);
            }
            X.push_back(row);
        }
    }else{
        for(int i=1; i<=Q; i++){
            Lcur=size.GetLength(i);
            row.clear();
            for(int j=1; j<=Lcur; j++){
                curVal=DefaultVal;
                row.push_back(curVal);
            }
            X.push_back(row);
        }
    }
}
//
//get val
template<typename T> T Array2DGetVal(T**X, const Array2DSize&size, int extRowN, int ineRowN){//40
    T y;
    if(extRowN>=1 && extRowN<=size.GetQExtRows() && ineRowN>=1 && ineRowN<=size.GetLength(extRowN)){
       y=X[extRowN][ineRowN];
    }
    return y;
}
template<typename T> T* Array2DGetVal_AsPtr(T**X, const Array2DSize&size, int extRowN, int ineRowN){//40
    T y;
    if(extRowN>=1 && extRowN<=size.GetQExtRows() && ineRowN>=1 && ineRowN<=size.GetLength(extRowN)){
       y=X[extRowN][ineRowN];
    }
    return y;
}
template<typename T> T Array2DGetVal(std::vector<std::vector<T>>X, int extRowN, int ineRowN){//40
    T y;
    if(extRowN>=1 && extRowN<=X.size() && ineRowN>=1 && ineRowN<=X[extRowN-1].size){
       y=X[extRowN][ineRowN];
    }
    return y;
}
template<typename T> T* Array2DGetVal_AsPtr(std::vector<std::vector<T>>X, int extRowN, int ineRowN){//40
    T*y=NULL;
    if(extRowN>=1 && extRowN<=X.size() && ineRowN>=1 && ineRowN<=X[extRowN-1].size()){
       y=X[extRowN][ineRowN];
    }
    return y;
}
//
template<typename T> std::vector<T> Array2DGetExtRowN_AsVals(T**X, int ExtRowN){//30
    std::vector<T>y;
    int L;
    if(ExtRowN>=1 && ExtRowN<=X.size()){
        L=X[ExtRowN-1].size();
        for(int i=1; i<=L; i++){
            y.push_back(X[ExtRowN-1][i-1]);
        }
    }
    return y;
}
template<typename T> T* Array2DGetExtRowN_AsPtr(T**X, const Array2DSize&size, int ExtRowN){//30
    T*y=NULL;
    if(size.ExtRowNBelongsTo(ExtRowN)){
        y=X[ExtRowN-1];
    }
    return y;
}
template<typename T> std::vector<T> Array2DGetExtRowN(std::vector<std::vector<T>>X, int ExtRowN){//30
    std::vector<T>y;
    if(ExtRowN>=1 && ExtRowN<=X.size()){
        y=X[ExtRowN-1];
    }
    return y;
}
//transpose
template <typename T>bool SizeIsForTtanspose(std::vector<std::vector<T>>X){
    Array2DSize size=GetSizeOf2DVector(X);
    return size.IsForTranspose();
}
template<typename T>void Array2DTranspose(T**&X, Array2DSize&size){//
    T**Y=NULL;
    int Q=size.GetQExtRows(), Lmax=size.GetMaxLength(), Lmin=size.GetMinLength(), Lcur;
    //if(X!=NULL && Q>0 && L>0 && size.GetMinLength()==L){
    bool isForTranspose=size.IsForTranspose();
    if(Q>0 && Lmax>0 && Lmin>0 && isForTranspose){
    //LengthFullOfIneRowN(int N){//from start to 1st stop
        Y=new T*[Lmax];
        for(int i=1; i<=Lmax; i++){
            Lcur=size.LengthFullOfIneRowN(i);
            Y[i-1]=new T[Lcur];
        }
        //
        for(int i=1; i<=Lmax; i++){
            Lcur=size.LengthFullOfIneRowN(i);
            for(int j=1; j<=Lcur; j++){
                Y[i-1][j-1]=X[j-1][i-1];
            }
        }
        //
        for(int i=1; i<=Q; i++){
            delete[]X[i-1];
        }
        delete[]X;
        //
        X=Y;
        //
        size.Transpose();
    }
}
template<typename T>void Array2DTranspose(std::vector<std::vector<T>>&X){
    std::vector<std::vector<T>>Y;
    std::vector<T>Z;
    T val;
    Array2DSize size=GetSizeOf2DVector(X);
    int Q=size.GetQExtRows(), Lmax=size.GetMaxLength(), Lmin=size.GetMinLength(), hcur;
    bool isForTranspose=size.IsForTranspose();
    if(Q>0 && Lmax>0 && Lmin>0 && isForTranspose){
        for(int i=1; i<=Lmax; i++){
            Z.clear();
            hcur=size.LengthFullOfIneRowN(i);
            for(int j=1; j<=hcur; j++){
                val=X[j-1][i-1];
                Z.push_back(val);
            }
            Y.push_back(Z);
        }
        X=Y;
    }
    //X=Y;//may gbe this is more correct then place this operator before as is realized
}
template<typename T>void Array2DTransposeTo(std::vector<std::vector<T>>X){
    std::vector<std::vector<T>>Y;
    Y=X;
    Array2DTranspose(Y);
    return Y;
}

/*
template<typename T>std::vector<std::vector<T>> Array2DTransposeTo(T**X, const Array2DSize&size){//
    std::vector<std::vector<T>>Y;
    std::vector<T>Z;
    ///if(size.GetMaxLength()==L){
    //    for(int i=1; i<=L; i++){
    //        Z.clear();
    //        for(int j=1; j<=Q; j++){
    //            Z.push_back(X[j-1][i-1]);
    //        }
    //        Y.push_back(Z);
    //    }
    //}
    int Q=size.GetQExtRows(), Lmax=size.GetMaxLength(), Lmin=size.GetMinLength(), Lcur;
    bool isForTranspose=size.IsForTranspose();
    if(Q>0 && Lmax>0 && Lmin>0 && isForTranspose){
        for(int i=1; i<=Lmax; i++){
            Z.clear();
            Lcur=size.LengthFullOfIneRowN(i);
            for(int j=1; j<=Lcur; j++){
               Z.push_back(X[j-1][i-1]);
            }
            Y.push_back(Z);
        }
    }
    return Y;
}
template<typename T>std::vector<std::vector<T>> Array2DTransposeTo(std::vector<std::vector<T>>X){//rect only
    std::vector<std::vector<T>>Y;
    std::vector<T>Z;
    int Q=X.size(), L=X[1-1].size();
    for(int i=1; i<=L; i++){
        Z.clear();
        for(int j=1; j<=Q; j++){
           Z.push_back(X[j-1][i-1]);
       }
       Y.push_back(Z);
    }
    return Y;
}
//transose struct // not tested
template<typename T>std::vector<std::vector<T>> Array2DTransposeStructTo(T**&X, Array2DSize&size){
    std::vector<std::vector<T>>Y;
    std::vector<T> rowCur;
    int Q=size.GetQExtRows, Lmin=size.GetMinLength(), Lmax=size.GetMaxLength(), Lcur, LIneCur, L1=size.GetLength(1);
    if (size.IsForTranspose()){
        for(int i=1; i<=L1; i++){
            LIneCur=size.LengthFullOfIneRowN(i);
            rowCur.clear();
            rowCur=Array2DGetIneRowN(X, i);
            Y.push_back(rowCur);
        }
    }
    return Y;
}
template<typename T>void Array2DTransposeStruct(T**&X, Array2DSize&size){
    T**Y=NULL;
    std::vector<T> rowCur;
    int Q=size.GetQExtRows, Lmin=size.GetMinLength(), Lmax=size.GetMaxLength(), Lcur, LIneCur, L1=size.GetLength(1);
    if (size.IsForTranspose()){
        //create
        Y=new T*[L1];
        for(int i=1; i<=L1; i++){
            LIneCur=size.LengthFullOfIneRowN(i);
            Y[i-1]=new T[LIneCur];
        }
        //copy
        for(int i=1; i<=L1; i++){
            LIneCur=size.LengthFullOfIneRowN(i);
            rowCur.clear();
            rowCur=Array2DGetIneRowN(X, i);
            for(int j=1; j<=LIneCur; j++){
                Y[i-1][j-1]=X[j-1][i-1];
            }
        }
        //del
        for(int i=1; i<=Q; i++){
            delete[]X[i-1];
        }
        delete[]X;
        //assign
        X=Y;
    }
    return Y;
}
template<typename T>std::vector<std::vector<T>> Array2DTransposeStructTo(std::vector<std::vector<T>>X){
    Array2DSize size=GetSizeOf2DVector(X);
    std::vector<std::vector<T>>Y;
    std::vector<T> rowCur;
    int Q=size.GetQExtRows, Lmin=size.GetMinLength(), Lmax=size.GetMaxLength(), Lcur, LIneCur, L1=size.GetLength(1);
    if (size.IsForTranspose()){
        for(int i=1; i<=L1; i++){
            LIneCur=size.LengthFullOfIneRowN(i);
            rowCur.clear();
            rowCur=Array2DGetIneRowN(X, i);
            Y.push_back(rowCur);
        }
    }
    return Y;
}
template<typename T>void Array2DTransposeStructTo(std::vector<std::vector<T>>&X){
    Array2DSize size=GetSizeOf2DVector(X);
    std::vector<std::vector<T>>Y;
    std::vector<T> rowCur;
    int Q=size.GetQExtRows, Lmin=size.GetMinLength(), Lmax=size.GetMaxLength(), Lcur, LIneCur, L1=size.GetLength(1);
    if (size.IsForTranspose()){
        for(int i=1; i<=L1; i++){
            LIneCur=size.LengthFullOfIneRowN(i);
            rowCur.clear();
            rowCur=Array2DGetIneRowN(X, i);
            Y.push_back(rowCur);
        }
    }
    X=Y;
}
*/
//--------------------------------------------------
//Seek
template <typename T> std::vector<std::vector<int>> Array2D_SeekVal_Simply(T**X, const Array2DSize& size, const T& Val, int paramN=0){
    //std::vector<T>rowVect=Array1DAssign(row, L);
    std::vector<std::vector<int>>NsAll;
    std::vector<int>NsSnglLine;
    std::vector<int>pair;
    int Q=size.GetQExtRows(), count, N;
    int FromN=1, ToN=0;
    for(int i=1; i<=Q; i++){
        pair.clear();
        count=0;
        NsSnglLine=Array1DSeekValSimply(X[i-1], Val, size.GetLength(i), 1, 0, paramN);
        count=NsSnglLine.size();
        if(count>0){
            for(int j=1; j<=count; j++){
                N=NsSnglLine[j-1];
                pair.push_back(i);
                pair.push_back(N);
                NsAll.push_back(pair);
            }
        }
    }
    return NsAll;
}

template <typename T> std::vector<std::vector<int>> Array2D_SeekExtRow_Simply(T**X, const Array2DSize& size, T*row, int L, int paramN=0){
    std::vector<T>rowVect=Array1DAssign(row, L);
    std::vector<std::vector<int>>NsAll;
    std::vector<int>NsSnglLine;
    std::vector<int>pair;
    int Q=size.GetQExtRows(), count, N;
    int FromN=1, ToN=0;
    for(int i=1; i<=Q; i++){
        pair.clear();
        count=0;
        NsSnglLine=Array1DSeekSubArraySimply(X[i-1], row, size.GetLength(i), size.GetLength(i-1), L, FromN, ToN, paramN);
        count=NsSnglLine.size();
        if(count>0){
            for(int j=1; j<=count; j++){
                N=NsSnglLine[j-1];
                pair.push_back(i);
                pair.push_back(N);
                NsAll.push_back(pair);
            }
        }
    }
    return NsAll;
}
template <typename T> std::vector<std::vector<int>> Array2D_SeekExtRow_Simply(T**X, const Array2DSize& size, std::vector<T>row, int paramN=0){
    Array2D_SeekExtRow(row.data(), row.size(), paramN);
}


template <typename T>bool Array2D_Arrs2ndIsAtPosIn1st(T**Where, const Array2DSize& sizeWhere, T**What, const Array2DSize& sizeWhat, int ExtRowN, int IneRowN, bool IsTransposed=false, bool EachExtRowIsReversed=false){
    int whereQ=sizeWhere.GetQExtRows(), whatQ=sizeWhat.GetQExtRows(), whereL, whatL, whereH, whatH, whereLcur, whatLcur, whereExtRowN, whereIneRowN, whatExtRowN, whatIneRowN;
    bool verdict=true;
    if(IsTransposed==false){
        if(ExtRowN+whatQ>=whereQ){
            verdict=false;
        }else{
            if(EachExtRowIsReversed==false){
                // 11 12 13
                // 21 22
                // 31 32 33
                for(whatExtRowN=1; whatExtRowN<=whatQ; whatExtRowN++){
                    whereExtRowN=ExtRowN+whatExtRowN-1;
                    whereL=sizeWhere.GetLength(whereExtRowN);
                    whatL=sizeWhat.GetLength(whatExtRowN);
                    if(IneRowN-1+whatL>whereL){
                        verdict=false;
                    }
                }
                if(!verdict==false){
                    for(whatExtRowN=1; whatExtRowN<=whatQ; whatExtRowN++){
                        whereExtRowN=ExtRowN+whatExtRowN-1;
                        whereL=sizeWhere.GetLength(whereExtRowN);
                        whatL=sizeWhat.GetLength(whatExtRowN);
                        for(whatIneRowN=1; whatIneRowN<=whatL; whatIneRowN++){
                            whereIneRowN=IneRowN+whatIneRowN-1;
                            if(Where[whereExtRowN-1][whereIneRowN-1]!=What[whatExtRowN-1][whatIneRowN-1]){
                                verdict=false;
                            }
                        }
                    }
                }
            }else{//EachExtRowIsReversed==true
                // 13 12 11
                //    22 21
                // 33 32 31
                for(whatExtRowN=1; whatExtRowN<=whatQ; whatExtRowN++){
                    whereExtRowN=ExtRowN+whatExtRowN-1;
                    whereL=sizeWhere.GetLength(whereExtRowN);
                    whatL=sizeWhat.GetLength(whatExtRowN);
                    if(IneRowN-whatL+1<1){
                        verdict=false;
                    }
                }
                if(!verdict==false){
                    for(whatExtRowN=1; whatExtRowN<=whatQ; whatExtRowN++){
                        whereExtRowN=ExtRowN+whatExtRowN-1;
                        whereL=sizeWhere.GetLength(whereExtRowN);
                        whatL=sizeWhat.GetLength(whatExtRowN);
                        for(whatIneRowN=1; whatIneRowN<=whatL; whatIneRowN++){
                            whereIneRowN=IneRowN-whatIneRowN+1;
                            if(Where[whereExtRowN-1][whereIneRowN-1]!=What[whatExtRowN-1][whatIneRowN-1]){
                                verdict=false;
                            }
                        }
                    }
                }
            }
        }
    }else{//transposed //task s'harder
        if(EachExtRowIsReversed==false){
            // 11 21 31
            // 12 22 32
            // 13    33
            for(whatExtRowN=1; whatExtRowN<=whatQ; whatExtRowN++){
                whereIneRowN=IneRowN+whatExtRowN-1;
                whatL=sizeWhat.GetLength(whatExtRowN);
                for(whatIneRowN=1; whatIneRowN<=whatL; whatIneRowN++){
                    whereExtRowN=ExtRowN+whatIneRowN-1;
                    whereL=sizeWhere.GetLength(whereExtRowN);
                    if(whereExtRowN>=1 && whereExtRowN<=whereQ && whereIneRowN>=1 &&whereIneRowN<=whereL){
                        //NOp
                    }else{
                        verdict=false;
                    }
                }
            }
            if(!verdict==false){
                for(whatExtRowN=1; whatExtRowN<=whatQ; whatExtRowN++){
                    whereIneRowN=IneRowN+whatExtRowN-1;
                    whatL=sizeWhat.GetLength(whatExtRowN);
                    for(whatIneRowN=1; whatIneRowN<=whatL; whatIneRowN++){
                        whereExtRowN=ExtRowN+whatIneRowN-1;
                        whereL=sizeWhere.GetLength(whereExtRowN);
                        if(Where[whereExtRowN-1][whereIneRowN-1]!=What[whatExtRowN-1][whatIneRowN-1]){
                            verdict=false;
                        }
                    }
                }
            }
        }else{
            // 13    33
            // 12 22 32
            // 11 21 31
            for(whatExtRowN=1; whatExtRowN<=whatQ; whatExtRowN++){
                whereIneRowN=IneRowN+whatExtRowN-1;
                whatL=sizeWhat.GetLength(whatExtRowN);
                for(whatIneRowN=1; whatIneRowN<=whatL; whatIneRowN++){
                    whereExtRowN=ExtRowN-whatIneRowN+1;
                    whereL=sizeWhere.GetLength(whereExtRowN);
                    if(whereExtRowN>=1 && whereExtRowN<=whereQ && whereIneRowN>=1 &&whereIneRowN<=whereL){
                        //NOp
                    }else{
                        verdict=false;
                    }
                }
            }
            if(!verdict==false){
                for(whatExtRowN=1; whatExtRowN<=whatQ; whatExtRowN++){
                    whereIneRowN=IneRowN+whatExtRowN-1;
                    whatL=sizeWhat.GetLength(whatExtRowN);
                    for(whatIneRowN=1; whatIneRowN<=whatL; whatIneRowN++){
                        whereExtRowN=ExtRowN-whatIneRowN+1;
                        whereL=sizeWhere.GetLength(whereExtRowN);
                        if(Where[whereExtRowN-1][whereIneRowN-1]!=What[whatExtRowN-1][whatIneRowN-1]){
                            verdict=false;
                        }
                    }
                }
            }
        }
    }
    return verdict;
}
template <typename T>bool Array2D_Arrs2ndIsAtPosIn1st(std::vector<std::vector<T>>Where,  std::vector<std::vector<T>>What, int ExtRowN, int IneRowN, bool IsTransposed=false, bool EachExtRowIsReversed=false){
    Array2DSize sizeWhere=GetSizeOf2DVector(Where), sizeWhat=GetSizeOf2DVector(What);
    int whereQ=sizeWhere.GetQExtRows(), whatQ=sizeWhat.GetQExtRows(), whereL, whatL, whereH, whatH, whereLcur, whatLcur, whereExtRowN, whereIneRowN, whatExtRowN, whatIneRowN;
    bool verdict=true;
    if(IsTransposed==false){
        if(ExtRowN+whatQ>=whereQ){
            verdict=false;
        }else{
            if(EachExtRowIsReversed==false){
                // 11 12 13
                // 21 22
                // 31 32 33
                for(whatExtRowN=1; whatExtRowN<=whatQ; whatExtRowN++){
                    whereExtRowN=ExtRowN+whatExtRowN-1;
                    whereL=sizeWhere.GetLength(whereExtRowN);
                    whatL=sizeWhat.GetLength(whatExtRowN);
                    if(IneRowN-1+whatL>whereL){
                        verdict=false;
                    }
                }
                if(!verdict==false){
                    for(whatExtRowN=1; whatExtRowN<=whatQ; whatExtRowN++){
                        whereExtRowN=ExtRowN+whatExtRowN-1;
                        whereL=sizeWhere.GetLength(whereExtRowN);
                        whatL=sizeWhat.GetLength(whatExtRowN);
                        for(whatIneRowN=1; whatIneRowN<=whatL; whatIneRowN++){
                            whereIneRowN=IneRowN+whatIneRowN-1;
                            if(Where[whereExtRowN-1][whereIneRowN-1]!=What[whatExtRowN-1][whatIneRowN-1]){
                                verdict=false;
                            }
                        }
                    }
                }
            }else{//EachExtRowIsReversed==true
                // 13 12 11
                //    22 21
                // 33 32 31
                for(whatExtRowN=1; whatExtRowN<=whatQ; whatExtRowN++){
                    whereExtRowN=ExtRowN+whatExtRowN-1;
                    whereL=sizeWhere.GetLength(whereExtRowN);
                    whatL=sizeWhat.GetLength(whatExtRowN);
                    if(IneRowN-whatL+1<1){
                        verdict=false;
                    }
                }
                if(!verdict==false){
                    for(whatExtRowN=1; whatExtRowN<=whatQ; whatExtRowN++){
                        whereExtRowN=ExtRowN+whatExtRowN-1;
                        whereL=sizeWhere.GetLength(whereExtRowN);
                        whatL=sizeWhat.GetLength(whatExtRowN);
                        for(whatIneRowN=1; whatIneRowN<=whatL; whatIneRowN++){
                            whereIneRowN=IneRowN-whatIneRowN+1;
                            if(Where[whereExtRowN-1][whereIneRowN-1]!=What[whatExtRowN-1][whatIneRowN-1]){
                                verdict=false;
                            }
                        }
                    }
                }
            }
        }
    }else{//transposed //task s'harder
        if(EachExtRowIsReversed==false){
            // 11 21 31
            // 12 22 32
            // 13    33
            for(whatExtRowN=1; whatExtRowN<=whatQ; whatExtRowN++){
                whereIneRowN=IneRowN+whatExtRowN-1;
                for(whatIneRowN=1; whatIneRowN<=whatQ; whatIneRowN++){
                    whereExtRowN=ExtRowN+whatIneRowN-1;
                    whereL=sizeWhere.GetLength(whereExtRowN);
                    if(whereExtRowN>=1 && whereExtRowN<=whereQ && whereIneRowN>=1 &&whereIneRowN<=whereL){
                        verdict=false;
                    }
                }
            }
            if(!verdict==false){
                for(whatExtRowN=1; whatExtRowN<=whatQ; whatExtRowN++){
                    whereIneRowN=IneRowN+whatExtRowN-1;
                    for(whatIneRowN=1; whatIneRowN<=whatQ; whatIneRowN++){
                        whereExtRowN=ExtRowN+whatIneRowN-1;
                        whereL=sizeWhere.GetLength(whereExtRowN);
                        if(Where[whereExtRowN-1][whereIneRowN-1]!=What[whatExtRowN-1][whatIneRowN-1]){
                            verdict=false;
                        }
                    }
                }
            }
        }else{
            // 13    33
            // 12 22 32
            // 11 21 31
            for(whatExtRowN=1; whatExtRowN<=whatQ; whatExtRowN++){
                whereIneRowN=IneRowN+whatExtRowN-1;
                for(whatIneRowN=1; whatIneRowN<=whatQ; whatIneRowN++){
                    whereExtRowN=ExtRowN-whatIneRowN+1;
                    whereL=sizeWhere.GetLength(whereExtRowN);
                    if(whereExtRowN>=1 && whereExtRowN<=whereQ && whereIneRowN>=1 &&whereIneRowN<=whereL){
                        verdict=false;
                    }
                }
            }
            if(!verdict==false){
                for(whatExtRowN=1; whatExtRowN<=whatQ; whatExtRowN++){
                    whereIneRowN=IneRowN+whatExtRowN-1;
                    for(whatIneRowN=1; whatIneRowN<=whatQ; whatIneRowN++){
                        whereExtRowN=ExtRowN-whatIneRowN+1;
                        whereL=sizeWhere.GetLength(whereExtRowN);
                        if(Where[whereExtRowN-1][whereIneRowN-1]!=What[whatExtRowN-1][whatIneRowN-1]){
                            verdict=false;
                        }
                    }
                }
            }
        }
    }
    return verdict;
}

/*template <typename T> bool Array2D_SubArr2DIsAtPos_SameDataSucc(T**Where, const Array2DSize& sizeWhere, T**What, const Array2DSize& sizeWhat, int ExtRowN, int IneRowN){
    bool b=true;
    int WhereExtRowN, WhereIneRowN, WhereCurL, WhatExtRowN, WhatIneRowN, WhatCurL, WhereQ=sizeWhere.GetQExtRows(), WhatQ=sizeWhat.GetQExtRows();
    T WhereVal, WhatVal;
    for(WhatExtRowN=1; WhatExtRowN<=WhatQ; WhatExtRowN++){
        WhereExtRowN=WhatExtRowN+ExtRowN-1;
        WhatCurL=sizeWhat.GetLength(WhatExtRowN);
        WhereCurL=sizeWhere.GetLength(WhereCurL);
        for(WhatIneRowN=1; WhatIneRowN<=WhatCurL; WhatIneRowN++){
            WhereIneRowN=IneRowN+WhatIneRowN-1;

            WhereVal=Where[WhereExtRowN-1][WhereIneRowN-1];
            WhatVal=What[WhatExtRowN-1][WhatIneRowN-1];
            if(WhereVal!=WhatVal)
        }
    }
    return b;
}*/








template <typename T> std::vector<std::vector<int>> Array2D_SeekArrs2ndIn1st(T**Where, const Array2DSize& sizeWhere, T**What, const Array2DSize& sizeWhat, bool IsTransposed=false, bool EachExtRowIsReversed=false){
    std::vector<std::vector<int>>Ns;
    std::vector<int>Coords;
    bool isPresentHere;
    int whereQ=sizeWhere.GetQExtRows(), Lcur;
    for(int ExtRowN=1; ExtRowN<=whereQ; ExtRowN++){
        Lcur=sizeWhere.GetLength(ExtRowN);
        for(int IneRowN=1; IneRowN<=Lcur; IneRowN++){
            isPresentHere=Array2D_Arrs2ndIsAtPosIn1st(Where, sizeWhere, What, sizeWhat, ExtRowN, IneRowN, IsTransposed, EachExtRowIsReversed);
            if(isPresentHere){
                Coords.clear();
                Coords.push_back(ExtRowN);
                Coords.push_back(IneRowN);
                Ns.push_back(Coords);
            }
        }
    }
    return Ns;
}
template <typename T> std::vector<std::vector<int>> Array2D_SeekArrs2ndIn1st(std::vector<std::vector<T>>Where,  std::vector<std::vector<T>>What, bool IsTransposed=false, bool EachExtRowIsReversed=false){
    Array2DSize sizeWhere=GetSizeOf2DVector(Where), sizeWhat=GetSizeOf2DVector(What);
    std::vector<std::vector<int>>Ns;
    std::vector<int>Coords;
    bool isPresentHere;
    int whereQ=sizeWhere.GetQExtRows(), Lcur;
    for(int ExtRowN=1; ExtRowN<=whereQ; ExtRowN++){
        Lcur=sizeWhere.GetLength(ExtRowN);
        for(int IneRowN=1; IneRowN<=Lcur; IneRowN++){
            //isPresentHere=Array2D_Arrs2ndIsAtPosIn1st(Where, sizeWhere, What, sizeWhat, ExtRowN, IneRowN, IsTransposed, EachExtRowIsReversed);
            isPresentHere=Array2D_Arrs2ndIsAtPosIn1st(Where, What, ExtRowN, IneRowN, IsTransposed, EachExtRowIsReversed);
            if(isPresentHere){
                Coords.clear();
                Coords.push_back(ExtRowN);
                Coords.push_back(IneRowN);
                Ns.push_back(Coords);
            }
        }
    }
    return Ns;
}



template <typename T> std::vector<std::vector<int>> Array2D_SeekIneRow_Simply(T**X, const Array2DSize& size, T*row, int L, int paramN=0){
    std::vector<std::vector<T>>data;
    Array2DAssign(data, L, row, true);

}

//--------------------------------------------------
//Concat
template <typename T>void Array2DConcatAdding(T**&X1, Array2DSize& size1, T**X2, const Array2DSize& size2, int shiftVal=0, T*DefaultValParam=NULL, bool Rect=true, int QRec=0, int FromN=1){
    int L1ini=size1.GetMaxLength(), L2ini=size2.GetMaxLength(),
        //Q1ini=size1.GetQExtRows(), 
        Q1=size1.GetQExtRows(),
        Q2ini=size2.GetQExtRows(),
        //Q1fin=Q1ini,
        //Q1=Q1ini;
        Q2fin, L1fin, L2fin, L;//, Lreq;
    std::vector<T>row;
    if(QRec==0){//Q to add is Q how many is given minus FromN
        Q2fin=Q2ini-FromN+1;//corr'd
    }else if(QRec<=Q2ini){
        Q2fin=QRec-FromN+1;//corr'd
    }
    if(shiftVal>=0){//2nd rows shifted to the left
        L1fin=L1ini;
        L2fin=shiftVal+L2ini;
        L = L1fin >= L2fin ? L1fin : L2fin;
        //if(QRec==0){//Q to add is Q how many is given minus FromN
        //    Q2fin=Q2ini-FromN+1;//corr'd
        //}else if(QRec<=Q2ini){
        //    Q2fin=QRec-FromN+1;//corr'd
        //}
        if(Rect){
            for(int i=1; i<=Q1; i++){
                row=Array1DGetSubArray_byLs_PreMark(X1[i-1], size1.GetLength(i), L, 1, 0, DefaultValParam);
                //template<typename T> void Array2DSetExtRowN(T**&X, Array2DSize&size, int ExtRowN, int wholeRowL=0, T*rowParam=NULL, int whatFromN=1, int QDefaultValuesBefore=0, bool KeepIfRect=true, T*DefaultValParam=NULL, bool DefaultIsSpecNotOwn=false){
                //
                Array2DSetExtRowN(X1, size1, i, row, DefaultValParam, 1, 0, Rect);
            }
            for(int i=FromN; i<=FromN+Q2fin-1; i++){
                row=Array1DGetSubArray_byLs_PreMark(X2[i-1], size2.GetLength(i), L, 1, shiftVal, DefaultValParam);//L
                Array2DAddExtRow(X1, size1, row, 1, 0, true, DefaultValParam);
            }
        }else{
            for(int i=FromN; i<=FromN+Q2fin-1; i++){
                row=Array1DGetSubArray_byLs_PreMark(X2[i-1], size2.GetLength(i), 0, 1, shiftVal, DefaultValParam);//0
                Array2DAddExtRow(X1, size1, row, 1, 0, true, DefaultValParam);
            }
        }//gut, I men os ver
    }else{//shiftVal<0 //1st rows shifted to the left sdi
        L1fin=-shiftVal+L1ini;
        L2fin=L2ini;
        L = L1fin >= L2fin ? L1fin : L2fin;
        //
        if(Rect){
            for(int i=1; i<=Q1; i++){
                row=Array1DGetSubArray_byLs_PreMark(X1[i-1], size1.GetLength(i), L, 1, shiftVal, DefaultValParam);
                Array2DSetExtRowN(X1, size1, i, row, DefaultValParam, 1, 0, true);
            } 
            for(int i=FromN; i<= FromN+Q2fin-1; i++){
                row=Array1DGetSubArray_byLs_PreMark(X2[i-1], size2.GetLength(i), L, 1, 0, DefaultValParam);//L
                Array2DAddExtRow(X1, size1, row, 1, 0, true, DefaultValParam);
            }
        }else{
            for(int i=1; i<=Q1; i++){
                row=Array1DGetSubArray_byLs_PreMark(X1[i-1], size1.GetLength(i), 0, 1, shiftVal, DefaultValParam);
                //Array2DSetExtRowN(X1, size1, i, row, 1, 0, true, DefaultValParam, false);
                Array2DSetExtRowN(X1, size1, i, row, DefaultValParam, 1, 0, false);
            } 
            for(int i=FromN; i<= FromN+Q2fin-1; i++){
                Array2DAddExtRow(X1, size1, X2[i-1], size2.GetLength(i), 1, 0, false, DefaultValParam);
            }
        }
    }
}
template <typename T>void Array2DConcatAdding(std::vector<std::vector<T>>&X1, std::vector<std::vector<T>>X2, int shiftVal=0, T*DefaultValParam=NULL, bool Rect=true, int QRec=0, int FromN=1){
    Array2DSize size1=GetSizeOf2DVector(X1), size2=GetSizeOf2DVector(X2);
    int L1ini=size1.GetMaxLength(), L2ini=size2.GetMaxLength(),
         //Q1ini=size1.GetQExtRows(),
         Q1=size1.GetQExtRows(),
         Q2ini=size2.GetQExtRows(),
         //Q1fin=Q1ini,
         //Q1=Q1ini;
         Q2fin, L1fin, L2fin, L;//, Lreq;
    int Q;
    std::vector<T>row;
    bool contin;
    int whereCurN, whatCurN,  whereCurL, whatCurL;
    if(QRec==0){//Q to add is Q how many is given minus FromN
        Q2fin=Q2ini-FromN+1;//corr'd
    }else if(QRec<=Q2ini){
        Q2fin=QRec-FromN+1;//corr'd
    }
    Q=Q1+Q2fin;
    if(shiftVal>=0){
         L = L1ini >= shiftVal+L2ini ? L1ini : shiftVal+L2ini;
         if(Rect){
             for(int i=1; i<=Q1; i++){
                 row=Array1DGetSubArray_byLs_PreMark(X1[i-1], L, 1, 0, DefaultValParam);
                 Array2DSetExtRowN(X1, i, row);
             }
             for(int i=FromN; i<=FromN+Q2fin-1; i++){
                 row=Array1DGetSubArray_byLs_PreMark(X2[i-1], L, 1, shiftVal, DefaultValParam);
                 Array2DAddExtRow(X1, row);
             }
         }else{
             for(int i=FromN; i<=FromN+Q2fin-1; i++){
                 row=Array1DGetSubArray_byLs_PreMark(X2[i-1], 0, 1, shiftVal, DefaultValParam);
                 Array2DAddExtRow(X1, row);
             }
         }
    }else{
        L = -shiftVal+L1ini >= L2ini ? -shiftVal+L1ini : L2ini;
        if(Rect){
            for(int i=1; i<=Q1; i++){
                row=Array1DGetSubArray_byLs_PreMark(X1[i-1], L, 1, -shiftVal, DefaultValParam);
                Array2DSetExtRowN(X1, i, row);
            }
            for(int i=FromN; i<=FromN+Q2fin-1; i++){
                row=Array1DGetSubArray_byLs_PreMark(X2[i-1], L, 1, 0, DefaultValParam);
                Array2DAddExtRow(X1, row);
            }
        }else{
            for(int i=1; i<=Q1; i++){
                row=Array1DGetSubArray_byLs_PreMark(X1[i-1], 0, 1, -shiftVal, DefaultValParam);
                Array2DSetExtRowN(X1, i, row);
            }
            for(int i=FromN; i<=FromN+Q2fin-1; i++){
                row=Array1DGetSubArray_byLs_PreMark(X2[i-1], 0, 1, 0, DefaultValParam);
                Array2DAddExtRow(X1, row);
            }
        }
    }
}

template <typename T>void Array2DConcatStretching(T**&X1, Array2DSize& size1, T**X2, const Array2DSize& size2, int shiftVal=0, T*DefaultValParam=NULL, bool Rect=true, bool IfFirstNotRect_AddToLastNotSrtetchToMax=true, int whatFromN=1, int QRec=0, TValsShowHide*vsh=NULL){
    int Lmax1=size1.GetMaxLength(), Lmin1=size1.GetMinLength(), Lmax2=size2.GetMaxLength(), Lreq, Q1=size1.GetMaxLength(), Q2=size2.GetMaxLength(), Qmin, Q, whereCurN, whatCurN, whatCurL, whereCurLini, whereCurLtmp;
    std::vector<T>row;
    T*ptr=NULL;
    T CurVal, DefaultVal;
    bool contin=false;
    if(DefaultValParam!=NULL){
        DefaultVal=(*(DefaultValParam));
    }
    if(Rect){
        Lreq=Lmax1+Lmax2;
    }else{
        Lreq=0;
    }
    if(QRec==0){
        QRec=size2.GetQExtRows()-whatFromN+1;
    }
    if(!IfFirstNotRect_AddToLastNotSrtetchToMax && Lmax1!=Lmin1){
        Array2DSetRect(X1, size1, Lmax1);
        if(vsh!=NULL){
            std::cout<<"Now all rows have same length:"<<std::endl;
            Arr2DStdCOut2D(X1, size1, " ", "; ", false);
        }
    }
    if(shiftVal>=0){//=> added row shifts down
        Qmin = Q1 <= shiftVal+QRec ? Q1 : shiftVal+QRec;
        Q = Q1 >= shiftVal+QRec ? Q1 : shiftVal+QRec;
        //
        //if(Rect){
            if(shiftVal<Q1){
                //shifted rows
                if(vsh!=NULL){
                    std::cout<<"Shifted empty rows "<<std::endl;
                }
                for(int i=1; i<=shiftVal; i++){
                    row=Array1DGetSubArray_byLs_PreMark(X1[i-1], size1.GetLength(i), Lreq, 1, 0, DefaultValParam);
                    Array2DSetExtRowN(X1, size1, i, row, DefaultValParam, 1, 0, false);
                    if(vsh!=NULL){
                        std::cout<<"row "<<i<<":";
                        Array1DShowConsole(X1[i-1], size1.GetLength(i));
                        std::cout<<std::endl;
                    }
                }
                //now placing after shift, combined rows
                if(vsh!=NULL){
                    std::cout<<"Rows combined after shift"<<std::endl;
                }
                whatCurN=whatFromN;
                whereCurN=shiftVal+whatCurN-whatFromN+1;//iqv: 3+3-3+1
                whatCurN-=1;
                whereCurN-=1;
                contin=(whatCurN<whatFromN+QRec-1 && whereCurN<Q1);
                while(contin){
                    whatCurN++;
                    whereCurN++;
                    //for(whatCurN=whatFromN; whatCurN<=whatFromN+QRec-1; whatCurN++){
                    whatCurL=size2.GetLength(whatCurN);
                    whereCurN=shiftVal+whatCurN-whatFromN+1;
                    whereCurLini=size1.GetLength(whereCurN);
                    whereCurLtmp=whereCurLini;
                    for(int j=1; j<=whatCurL; j++){
                        CurVal=X2[whatCurN-1][j-1];
                        //
                        //Array1DAddElement(X1[whereCurN-1], whereCurLtmp, CurVal);//also possible ety ably
                        //Array2DAddToExtRowN(T**&X, Array2DSize& size, int ExtRowN, T val)
                        Array2DAddToExtRowN(X1, size1, whereCurN, CurVal);

                    }
                    whereCurLini=whereCurLtmp;
                    whereCurLini=size1.GetLength(whereCurN);
                    for(int j=whereCurLini+1; j<=Lreq; j++){
                        Array1DAddElement(X1[whereCurN-1], whereCurLtmp, DefaultVal);
                    }
                    if(vsh!=NULL){
                        std::cout<<"row "<<whereCurN<<":";
                        Array1DShowConsole(X1[whereCurN-1], size1.GetLength(whereCurN));
                        std::cout<<std::endl;
                    }
                    contin=(whatCurN<whatFromN+QRec-1 && whereCurN<Q1);
                }//while
                //if(Rect){
                //    size1.Set(Q1, Lreq);
                //}
                //rest
                if(whatCurN==Q2 && whereCurN<Q1){
                    for(int i=shiftVal+QRec+1; i<=Q1; i++){
                        //row=Array1DAssign()
                        row=Array1DGetSubArray_byLs_PreMark(X1[i-1], size1.GetLength(i), Lreq, 1, 0, DefaultValParam);
                        Array2DSetExtRowN(X1, size1, i, row, DefaultValParam, 1, 0, false);
                    }
                }else if(whatCurN<Q2 && whereCurN==Q1){
                    for(int i=whatCurN+1; i<=Q2; i++){
                        row=Array1DGetSubArray_byLs_PreMark(X2[i-1], size2.GetLength(i), Lreq, 1, Lmax1, DefaultValParam);
                        Array2DAddExtRow(X1, size1, row);
                    }
                }
            }else{//shiftVal>Q1
                Q = -shiftVal+Q1 >= QRec ? -shiftVal+Q1 : QRec; //n'utf'd
                for(int i=1; i<=Q1; i++){
                    row=Array1DGetSubArray_byLs_PreMark(X1[i-1], size1.GetLength(i), Lreq, 1, 0, DefaultValParam);
                    Array2DSetExtRowN(X1, size1, i, row, DefaultValParam, 1, 0, false);
                }
                for(int i=Q1+1; i<=shiftVal; i++){
                    row=Array1DGetSubArray_byLs_PreMark(ptr, 0, Lreq, 1, 0, DefaultValParam);//ety ably
                    row=Array1DGetSubArray_byLs_PreMark(ptr, 0, Lreq, 1, Lreq, DefaultValParam);//ety ably
                    Array2DAddExtRow(X1, size1, row);
                }
                for(int i=1; i<=QRec; i++){
                    row=Array1DGetSubArray_byLs_PreMark(X2[i-1], size2.GetLength(i), Lreq, 1, Lmax1, DefaultValParam);
                    Array2DAddExtRow(X1, size1, row);
                }
                if(Rect){
                    size1.Set(Q1, Lreq);
                }
            }
        //}else{//not rect. Ma exda id ov s' uz tbi Rect et ne, Lreq fybe all cass

        //}
    }else{//shiftVal<0 //=> first row shifts down, new row starts with whatFromN. row o'X2
        //ShiftVals
        if(vsh!=NULL){
            std::cout<<"Shifted rows X2"<<std::endl;
        }
        whatCurN=whatFromN-1;
        whereCurN=0;
        //contin=(whatCurN<QRec+whatFromN-1 && whatCurN<Q1);
        contin=(whatCurN<QRec+whatFromN-1 && whatCurN-whatFromN+1<-shiftVal);//QRec is jq <=Q2
        //for(whatCurN=-shiftVal; whatCurN>=1; whatCurN--){
        //for(whatCurN=whatFromN; whatCurN<=-shiftVal+whatFromN-1; whatCurN++){
        while(contin){
            whatCurN++;
            whereCurN++;
            row=Array1DGetSubArray_byLs_PreMark(X2[whatCurN-1], size2.GetLength(whatCurN), Lreq, 1, Lmax1, DefaultValParam);
            //Array2DInsExtRow(X1, size1, 1, row, DefaultValParam, 1, 0, false); //false ob L ms'be new, S n'ved l'L's val
            Array2DInsExtRow(X1, size1, whereCurN, row, DefaultValParam, 1, 0, false); //false ob L ms'be new, S n'ved l'L's val
            //Array2DInsExtRow(X1, size1, 1, size2.GetLength(whatCurN), X2[whatCurN-1], DefaultValParam, 1, Lmax1, false);//false ob L ms'be new, S n'ved S
            contin=(whatCurN<QRec+whatFromN-1 && whatCurN-whatFromN+1<-shiftVal);
            if(vsh!=NULL){
                std::cout<<"row "<<whatCurN<<": ";
                Array1DShowConsole(X1[whereCurN-1], size1.GetLength(whereCurN));
            }
        }
        Q1=size1.GetQExtRows();
        //nu row o'X1, ic wa #1., ha Nr  -shiftVal+1;
        if(-shiftVal<whatFromN+QRec-1){ //ei rows X2 extend l'rows o'X1
            contin=(whatCurN<whatFromN+QRec-1 && whereCurN<Q1);
            while(contin){
                whatCurN++;
                whereCurN++;
                contin=(whatCurN<whatFromN+QRec-1 && whereCurN<Q1);
                if(vsh!=NULL){
                    std::cout<<"Stretching X1["<<whereCurN<<"(former "<<whereCurN+shiftVal<<")] with X2["<<whatCurN<<"]"<<std::endl;
                }
                whatCurL=size2.GetLength(whatCurN);
                for(int j=1; j<=whatCurL; j++){
                    CurVal=X2[whatCurN-1][j-1];
                    Array2DAddToExtRowN(X1, size1, whereCurN, CurVal);
                }
                whereCurLini=size1.GetLength(whereCurN);
                if(Rect){
                    for(int j=whereCurLini+1; j<=Lreq; j++){
                        Array2DAddToExtRowN(X1, size1, whereCurN, DefaultVal);
                    }
                }
                whereCurLini=size1.GetLength(whereCurN);
                if(vsh!=NULL){
                    Array1DShowConsole(X1[whereCurN-1],whereCurLini);
                }
            }
            //last rows
            if(whatCurN<whatFromN+QRec-1 && whereCurN==Q1){
                if(vsh!=NULL){
                    std::cout<<"rows of X2 remaining adding:"<<std::endl;
                }
                for(int i=whatCurN+1; i<=whatFromN+QRec-1; i++){
                    whatCurL=size2.GetLength(i);
                    row=Array1DGetSubArray_byLs_PreMark(X2[i-1], size2.GetLength(i), Lreq, 1, Lmax1, DefaultValParam);
                    Array2DAddExtRow(X1, size1, row, 1, 0, false, DefaultValParam);
                    //
                    whereCurN++;
                    whereCurLini=size1.GetLength(i);
                    if(vsh!=NULL){
                        std::cout<<"row"<<whereCurN<<" (L="<<whereCurLini<<"): ";
                        Array1DShowConsole(X1[whereCurN-1],whereCurLini);
                    }
                }
            }else if(whatCurN==whatFromN+QRec-1 && whereCurN<Q1){
                Lmin1=size1.GetMinLength();
                Lmax1=size1.GetMaxLength();
                if(Rect && Lmin1!=Lmax1){
                    if(vsh!=NULL){
                        std::cout<<"rows of X1 remaining strtetchig"<<std::endl;
                    }
                    for(int i=whereCurN+1; i<=Q1; i++){
                        whereCurLtmp=size1.GetLength(i);
                       //for(int j=whereCurLtmp+1; j<=Lreq; j++){
                        for(int j=whereCurLtmp+1; j<=Lreq; j++){
                            Array2DAddToExtRowN(X1, size1, i, DefaultVal);
                        }
                        //
                        //whereCurN++;
                       // whereCurLini=size1.GetLength(i);
                        if(vsh!=NULL){
                            std::cout<<"row"<<i<<" (L="<<whereCurLini<<"): ";
                            Array1DShowConsole(X1[i-1],whereCurLini);
                        }
                    }
                }else{
                    //else NOp, tic rows stay in X1 arr co cu lenis.
                    if(vsh!=NULL){
                        std::cout<<"rows of X1 remain without changes"<<std::endl;
                    }
                }
            }//else NOp, tbi nablb zq, ob so n'fin while, so tdi s'eq, sdi tbi arrs gead fin
        }
    }//shiftvals
    if(Rect){
        size1.SetIneRowsLength(-1);
    }
    if(vsh!=NULL){
        std::cout<<"Concat stretching finishes working, now we have "<<Q1<<" rows"<<std::endl;
    }
}//fn concat stretching
template <typename T> void Array2DConcatStretching(std::vector<std::vector<T>>&X1, std::vector<std::vector<T>>X2, int shiftVal=0, T*DefaultValParam=NULL, bool Rect=true, bool IfFirstNotRect_AddToLastNotSrtetchToMax=true, int whatFromN=1, int QRec=0, TValsShowHide*vsh=NULL){
    Array2DSize size1=GetSizeOf2DVector(X1), size2=GetSizeOf2DVector(X2);
    int Q1=size1.GetQExtRows(), Q2=size2.GetQExtRows(), Lmax1=size1.GetMaxLength(), Lmax2=size2.GetMaxLength(), Lmin1=size1.GetMinLength(), Lmin2=size2.GetMinLength();
    int whatCurN, whereCurN, whatCurL, whereCurL, Lreq, minQ, maxQ, Lmin3, Lmax3;//minN;
    bool contin;
    T DefaultVal, CurVal;
    std::vector<T> row;
    if(DefaultValParam!=NULL){
        DefaultVal=(*(DefaultValParam));
    }
    if(QRec==0){
        QRec=Q2-whatFromN+1;//ZB: 5, ab 3, all 3 = 5-3+1
    }
    if(IfFirstNotRect_AddToLastNotSrtetchToMax==false && Lmin1!=Lmax1){
        for(int i=1; i<=Q1; i++){
            whereCurL=size1.GetLength(i);
            for(int j=whereCurL+1; j<=Lmax1; j++){
                Array2DAddToExtRowN(X1, i, DefaultVal);
            }
        }
    }
    if(shiftVal>=0){//-----------------------------------------------------------------------------------
        minQ = shiftVal+QRec<=Q1 ? shiftVal+QRec : Q1; //f'add
        maxQ = shiftVal+QRec>=Q1 ? shiftVal+QRec : Q1; //f'size
        if(shiftVal<=Q1){
            for(int i=1; i<=minQ; i++){
                whereCurN=i+shiftVal-1;
                whatCurN=i+whatFromN-1;
                whereCurL=size1.GetLength(whereCurN);
                whatCurL=size2.GetLength(whatCurN);
                for(int j=1; j<=whatCurL; j++){
                    CurVal=X2[whatCurN-1][j-1];
                    Array2DAddToExtRowN(X1, whereCurN, CurVal);
                }
                size1.SetLength(whereCurL+whatCurL, whereCurN);
            }
            if(minQ<QRec){
                //=> N1+QRec > Q1
                if(IfFirstNotRect_AddToLastNotSrtetchToMax==false){
                    for(int i=1; i<=QRec-minQ; i++){
                        whereCurN=Q1+i;
                        whatCurN=i+whatFromN+minQ-1;
                        whatCurL=size2.GetLength(i);
                        row=Array1DGetSubArray_byLs_PreMark_v1(X2[whatCurN-1].data(), X2[whatCurN-1].size(), DefaultVal, Lmax1+whatCurL, 1, Lmax1);
                        //ambigous fn
                        Array1DAddElement(X1, row);
                        size1.AddExt(Lmax1+whatCurL);
                    }
                }else{
                    for(int i=1; i<=QRec-minQ; i++){
                        whereCurN=Q1+i;
                        whatCurN=i+whatFromN+minQ-1;
                        whatCurL=size2.GetLength(i);
                        row=X2[whatCurN-1];
                        Array1DAddElement(X1, row);
                        size1.AddExt(whatCurL);
                    }
                }
            }else{
                maxQ=shiftVal+QRec;
                //forming  shiftVal-Q1 rows between X1 and X2
                row.clear();
                if(IfFirstNotRect_AddToLastNotSrtetchToMax==false){
                    row=Array1DAssignRepeatingSingleVal(DefaultVal, Lmax1+0);
                    for(int i=Q1+1; i<=shiftVal; i++){
                        Array1DAddElement(X1, row);
                    }
                }else{
                    for(int i=Q1+1; i<=shiftVal; i++){
                        Array1DAddElement(X1, row);//row is empty
                        size1.AddExt(row.size());
                    }
                }
                //adding QRec rows of X2
                if(IfFirstNotRect_AddToLastNotSrtetchToMax==false){
                    for(int i=1; i<=QRec; i++){
                        whatCurN=whatFromN+i-1;
                        whatCurL=size2.GetLength(whatCurN);
                        row = Array1DGetSubArray_byLs_PreMark_v1(   X2[whatCurN-1].data(), X2[whatCurN-1].size(), DefaultVal, Lmax1+whatCurL, 1, Lmax1);
                        Array1DAddElement(X1, row);
                        size1.AddExt(row.size());
                    }
                }else{
                    for(int i=1; i<=QRec; i++){
                        whatCurN=whatFromN+i-1;
                        whatCurL=size2.GetLength(whatCurN);
                        //row = Array1DGetSubArray_byLs_PreMark(   X2[whatCurN-1].data(), X2[whatCurN-1].size(), DefaultVal, Lmax1+whatCurL, 1, 0);
                        //Array1DAddElement(X1, row);
                        //size1.AddExt(row.size());
                        //also gut, or simply
                        Array1DAddElement(X1, X2[whatCurN-1]);
                        size1.AddExt(X2[whatCurN-1].size());
                    }
                }
            }
        }
    }else{//---------------------------------------------------------------------------------------------
        if(-shiftVal<=QRec){
            if(-shiftVal+QRec<=Q1){
                if(IfFirstNotRect_AddToLastNotSrtetchToMax==false){
                    for(int i=1; i<=-shiftVal; i++){
                        whatCurN=whatFromN+i-1;
                        whatCurL=size2.GetLength(whatCurN);
                        row = Array1DGetSubArray_byLs_PreMark_v1(   X2[whatCurN-1].data(), X2[whatCurN-1].size(), DefaultVal, Lmax1+whatCurL, 1, Lmax1);
                        Array1DInsElementToN(X1, i, row);
                        size1.InsExt(i, row.size());
                    }
                    for(int i=-shiftVal+1; i<=QRec; i++){

                    }
                }else{

                }
            }else{

            }
        }else{
            for(int i=1; i<=-QRec; i++){

            }
            for(int i=QRec+1; i<=-shiftVal; i++){

            }
        }
    }//--------------------------------------------------------------------------------------------------
    if(Rect){
        Lmin3=size1.GetMinLength();
        Lmax3=size1.GetMaxLength();
        if(Lmin3!=Lmax3){
            for(int i=1; i<=maxQ; i++){
                whereCurL=size1.GetLength(i);
                for(int j=whereCurL+1; j<=Lmax3; j++){
                    Array1DAddElement(X1[i-1], DefaultVal);
                }
            }
        }
        size1.Set(maxQ, Lmax3);
    }
}

/*template <typename T> void Array2DConcatStretching(std::vector<std::vector<T>>&X1, std::vector<std::vector<T>>X2, int shiftVal=0, T*DefaultValParam=NULL, bool Rect=true, bool IfFirstNotRect_AddToLastNotSrtetchToMax=true, int whatFromN=1, int QRec=0, TValsShowHide*vsh=NULL){
    Array2DSize size1=GetSizeOf2DVector(X1), size2=GetSizeOf2DVector(X2);
    int Q1=size1.GetQExtRows(), Q2=size2.GetQExtRows(), Lmax1=size1.GetMaxLength(), Lmax2=size2.GetMaxLength(), Lmin1=size1.GetMinLength(), Lmin2=size2.GetMinLength();
    int whatCurN, whereCurN, whatCurL, whereCurL, Lreq, minN;
    bool contin;
    T DefaultVal, CurVal;
    std::vector<T> row;
    if(DefaultValParam!=NULL){
        DefaultVal=(*(DefaultValParam));
    }
    if(Rect){
        Lreq=Lmax1+Lmax2;
    }else{
        Lreq=0;
    }
    if(QRec==0){
        QRec=Q2-whatFromN+1;
    }
    if(shiftVal>=0){
        //first shiftVal rows of X1
        minN= shiftVal<=Q1 ? shiftVal: Q1;
        for(int i=Q1+1; i<=shiftVal; i++){//if shiftVal>Q1, else do ce 0 times
            Array2DAddExtRow(X1, row, 1, 0, false, DefaultValParam);
        }
        //combined
        whereCurN=shiftVal+1;
        if(shiftVal<Q1){
            whatCurN=whatFromN;
            whatCurN--;
            whereCurN--;
            contin=(whatCurN<whatFromN+QRec-1 && whereCurN<Q1);
            while(contin){
                whatCurN++;
                whereCurN++;
                whatCurL=X2[whatCurN-1].size();
                for(int j=1; j<=whatCurL; j++){
                    CurVal=X2[whatCurN-1][j-1];
                    Array2DAddToExtRowN(X1, whereCurN, CurVal);
                }
            }
        }
        //after combined
        if(whatCurN<whatFromN+QRec-1 && whereCurN==Q1){
            for(int i=whatCurN+1; i<=whatFromN+QRec-1; i++){
                Array2DAddExtRow(X1, X2[i-1], 1, 0, false, DefaultValParam);
            }
        }//else NOp
        if(Rect){
            Q1=X1.size();
            for(int i=1; i<=Q1; i++){
                whereCurL=X1[i-1].size();
                if(whereCurL<Lmax1+Lmax2){
                    row=Array1DGetSubArray_byLs_PreMark(X1[i-1], Lreq, 1, 0, DefaultValParam);
                    //Array2DSetExtRowN(std::vector<std::vector<T>>&X, int ExtRowN, std::vector<T>rowParam, int whatFromN=1, int QDefaultValuesBefore=0, bool KeepIfRect=true, T*DefaultValParam=NULL, bool DefaultIsSpecNotOwn=false){//, int ifDiff_ByArr0ByNew1Min2Max3=0){
                    Array2DSetExtRowN(                             X1,           i,                    row,               1,                          0,                false,   DefaultValParam,                               false);
                }
            }
        }
    }else{//shiftVal<0
        //shifted X2
        minN = -shiftVal <= whatFromN+QRec-1  ? -shiftVal : whatFromN+QRec-1;
        for(int i=1; i<=minN; i++){
            whatCurN=whatFromN+i-1;
            row=Array1DGetSubArray_byLs_PreMark(X2[whatCurN-1], Lreq, 1, Lmax1, DefaultValParam);
            Array2DInsExtRow(X1, whatCurN, X2[whatCurN-1]);
        }
        for(int i=whatFromN+QRec-1; i<=-shiftVal; i++){//if shiftVal>whatFromN+QRec-1, else do ce 0 times
            Array2DAddExtRow(X1, row, 1, 0, false, DefaultValParam);
        }
        //combined rows
        whatCurN=-shiftVal+1;
        whereCurN=-shiftVal+1;//ob tic X2 were inserted to X1, before own rows of X1
        whatCurN--;
        whereCurN--;
        contin=(whatCurN<whatFromN+QRec-1 && whereCurN<Q1);
        while(contin){
            whatCurN++;
            whereCurN++;
            whatCurL=X2[whatCurN-1].size();
            for(int j=1; j<=whatCurL; j++){
                CurVal=X2[whatCurN-1][j-1];
                Array2DAddToExtRowN(X1, whereCurN, CurVal);
            }
            contin=(whatCurN<whatFromN+QRec-1 && whereCurN<Q1);
        }
        //rest
        if(whatCurN<whatFromN+QRec-1 && whereCurN==Q1){
             for(int i=whatCurN+1; i<=whatFromN+QRec-1; i++){
                 row=Array1DGetSubArray_byLs_PreMark(X2[i-1], Lreq, 1, Lmax1, DefaultValParam);
                 Array2DAddExtRow(X1, row);
             }
        }//else if(whatCurN==whatFromN+QRec-1 && whereCurN<Q1){ NOp
        //if Rect
        Q1=X1.size();
        if(Rect){
            for(int i=1; i<=Q1; i++){
                whereCurL=X1[i-1].size();
                if(whereCurL<Lmax1+Lmax2){
                    //Array2DSetExtRowN(std::vector<std::vector<T>>&X, int ExtRowN, std::vector<T>rowParam, int whatFromN=1, int QDefaultValuesBefore=0, bool KeepIfRect=true, T*DefaultValParam=NULL, bool DefaultIsSpecNotOwn=false){//, int ifDiff_ByArr0ByNew1Min2Max3=0){
                    Array2DSetExtRowN(                             X1,           i,  Lreq, 1, 0, DefaultValParam, false);
                }
            }
        }
    }
}*/

//-------------------------------------------------
template <typename T>void Arr1DStdCOut(T*X, int L, QString delim="; ", bool useBrackets=true){
    if(useBrackets) std::cout<<"[";
    for(int i=1; i<=L-1; i++){
        std::cout<<X[i-1];
        std::cout<<delim.toStdString().c_str();
    }
    std::cout<<X[L-1];
    if(useBrackets)std::cout<<"]";
}
template <typename T>void Arr1DStdCOut(std::vector<T>X, int L, QString delim="; ", bool useBrackets=true){
    if(useBrackets) std::cout<<"[";
    for(int i=1; i<=L-1; i++){
        std::cout<<X[i-1];
        std::cout<<delim.toStdString().c_str();
    }
    std::cout<<X[L-1];
    if(useBrackets)std::cout<<"]";
}
//------------
template <typename T>void Arr2DStdCOut2D(T**X, const Array2DSize size, QString delimElem="; ", QString delimLines="; ", bool useBrackets=true){
    int Q=size.GetQExtRows();
    size.ShowToConsole();
    if(useBrackets) std::cout<<"[";
    Arr1DStdCOut(X[1-1], size.GetLength(1), delimElem, useBrackets);
    std::cout<<delimLines.toStdString().c_str();
    std::cout<<std::endl;
    for(int i=2; i<=Q-1; i++){
        //if(useBrackets && i>1)std::cout<<" ";
        if(useBrackets)std::cout<<" ";
        Arr1DStdCOut(X[i-1], size.GetLength(i), delimElem, useBrackets);
        std::cout<<delimLines.toStdString().c_str();
        std::cout<<std::endl;
    }
    if(useBrackets)std::cout<<" ";
    Arr1DStdCOut(X[Q-1], size.GetLength(Q), delimElem, useBrackets);
    if(useBrackets)std::cout<<"]";
    std::cout<<std::endl;
}
template <typename T>void Arr2DStdCOut2D(std::vector<std::vector<T>>X, QString delimElem="; ", QString delimLines="; ", bool useBrackets=true){
    Array2DSize size=GetSizeOf2DVector(X);
    int Q=size.GetQExtRows();
    size.ShowToConsole();
    if(useBrackets) std::cout<<"[";
    Arr1DStdCOut(X[1-1], size.GetLength(1), delimElem, useBrackets);
    std::cout<<delimLines.toStdString().c_str();
    std::cout<<std::endl;
    for(int i=2; i<=Q-1; i++){
        //if(useBrackets && i>1)std::cout<<" ";
        if(useBrackets)std::cout<<" ";
        Arr1DStdCOut(X[i-1], size.GetLength(i), delimElem, useBrackets);
        std::cout<<delimLines.toStdString().c_str();
        std::cout<<std::endl;
    }
    if(useBrackets)std::cout<<" ";
    Arr1DStdCOut(X[Q-1], size.GetLength(Q), delimElem, useBrackets);
    if(useBrackets)std::cout<<"]";
    std::cout<<std::endl;
}
template <typename T>void Arr2DStdCOut1D(T**X, const Array2DSize size, QString delimElem=", ", QString delimLines="; ", bool useBrackets=true){
    int Q=size.GetQExtRows();
    if(useBrackets) std::cout<<"[";
    for(int i=1; i<=Q-1; i++){
        Arr1DStdCOut(X[i-1], size.GetLength(i), delimElem, useBrackets);
        std::cout<<delimLines.toStdString().c_str();
    }
    Arr1DStdCOut(X[Q-1], size.GetLength(Q), delimElem, useBrackets);
    if(useBrackets)std::cout<<"]";
}
template <typename T>void Arr2DStdCOut1D(std::vector<std::vector<T>>X,  QString delimElem=", ", QString delimLines="; ", bool useBrackets=true){
     Array2DSize size=GetSizeOf2DVector(X);
    int Q=size.GetQExtRows();
    if(useBrackets) std::cout<<"[";
    for(int i=1; i<=Q-1; i++){
        Arr1DStdCOut(X[i-1], size.GetLength(i), delimElem, useBrackets);
        std::cout<<delimLines.toStdString().c_str();
    }
    Arr1DStdCOut(X[Q-1], size.GetLength(Q), delimElem, useBrackets);
    if(useBrackets)std::cout<<"]";
}

//=== wa2D arr === wi 2D arrbased on 1D ========================================================================

template<typename T> void ArrayOf2DBy1DSetSize(T*X, Array2DSize&size1, const Array2DSize&size2, T*DfltValParam=NULL, bool PreserveVals=true){
    int L1, L2, QER1=size1.GetQExtRows(), QER2=size2.GetQExtRows(), Qmin=QER1<=QER2 ? QER1 : QER2, QEl1=size1.GetQElements(), PosN;
    T DefaultVal;
    if(DfltValParam!=NULL){
        DefaultVal=(*(DfltValParam));
    }
    if(QER1>0){
        if(QER2>0){
            for(int i=1; i<=Qmin; i++){
                L1=size1.GetLength(i);
                L2=size2.GetLength(i);
                if(L2>L1){
                    PosN=Array2DSizeBasedOn1D_PosByCoords(i, L1, size1.GetLengthes());//also possible
                    PosN=size1.GetPositionByCoordsIfRealizedAs1D(i, L1);
                    for(int j=1; j<=(L2-L1); j++){
                        Array1DInsElementToN(X, QEl1, PosN, DefaultVal);
                    }
                }else if(L2<L1){
                    PosN=Array2DSizeBasedOn1D_PosByCoords(i, L2, size1.GetLengthes());//also possible
                    PosN=size1.GetPositionByCoordsIfRealizedAs1D(i, L2);
                    for(int j=1; j<=(L2-L1); j++){
                        Array1DDelElementFromN(X, QEl1, PosN);
                    }
                }
            }
            if(QER2>QER1){
                for(int i=Qmin+1; i<=QER2; i++){
                    L2=size2.GetLength(i);
                    for(int j=1; j<=L2; j++){
                        Array1DAddElement(X, QEl1, DefaultVal);
                    }
                }
            }else if(QER2<QER1){
                PosN=size1.GetPositionByCoordsIfRealizedAs1D(QER1, size1.GetLength(QER1))+1;
                for(int i=Qmin+1; i<=QER2; i++){
                    L1=size1.GetLength(i);
                    for(int j=1; j<=L1; j++){
                        Array1DDelElementFromN(X, QEl1, PosN);
                    }
                }
            }
        }else{
            delete[]X;
        }
    }
    size1=size2;
}
template<typename T> void ArrayOf2DBy1DSetSize(T*&X, const Array2DSize&size1, int Q, int*Ls, T*DfltValParam=NULL, bool PreserveVals=true){//16
    Array2DSize size2;
    size2.Set(Q, Ls);
    ArrayOf2DBy1DSetSize(X, size1, size2, DfltValParam, PreserveVals);
}
template<typename T> void ArrayOf2DBy1DSetSize(T*&X, Array2DSize&size1, int Q, int ExtRowsLength, T*DfltValParam=NULL, bool PreserveVals=true){//16
    Array2DSize size2;
    size2.Set(Q, ExtRowsLength);
    //template<typename T> void ArrayOf2DBy1DSetSize(T**&X, Array2DSize&size1, const Array2DSize&size2, T*DfltValParam=NULL, bool PreserveVals=true){//16p
    ArrayOf2DBy1DSetSize(                               X,              size1,                   size2,   DfltValParam,           PreserveVals);
 }
template<typename T> void ArrayOf2DBy1DSetSize(std::vector<T>&X, Array2DSize&size1, const Array2DSize&size2, T*DfltValParam=NULL){
    int L1, L2, QER1=size1.GetQExtRows(), QER2=size2.GetQExtRows(), Qmin=QER1<=QER2 ? QER1 : QER2, QEl1=size1.GetQElements(), PosN;
    T DefaultVal;
    if(DfltValParam!=NULL){
        DefaultVal=(*(DfltValParam));
    }
    if(QER1>0){
        if(QER2>0){
            for(int i=1; i<=Qmin; i++){
                L1=size1.GetLength(i);
                L2=size2.GetLength(i);
                if(L2>L1){
                    PosN=Array2DSizeBasedOn1D_PosByCoords(i, L1, size1.GetLengthes());//also possible
                    PosN=size1.GetPositionByCoordsIfRealizedAs1D(i, L1);
                    for(int j=1; j<=(L2-L1); j++){
                        Array1DInsElementToN(X, PosN, DefaultVal);
                    }
                }else if(L2<L1){
                    PosN=Array2DSizeBasedOn1D_PosByCoords(i, L2, size1.GetLengthes());//also possible
                    PosN=size1.GetPositionByCoordsIfRealizedAs1D(i, L2);
                    for(int j=1; j<=(L2-L1); j++){
                        Array1DDelElementFromN(X, PosN);
                    }
                }
            }
            if(QER2>QER1){
                for(int i=Qmin+1; i<=QER2; i++){
                    L2=size2.GetLength(i);
                    for(int j=1; j<=L2; j++){
                        Array1DAddElement(X, DefaultVal);
                    }
                }
            }else if(QER2<QER1){
                PosN=size1.GetPositionByCoordsIfRealizedAs1D(QER1, size1.GetLength(QER1))+1;
                for(int i=Qmin+1; i<=QER2; i++){
                    L1=size1.GetLength(i);
                    for(int j=1; j<=L1; j++){
                        Array1DDelElementFromN(X, PosN);
                    }
                }
            }
        }else{
            X.clear();
        }
    }
    size1=size2;
}
template<typename T> void ArrayOf2DBy1DSetSize(std::vector<T>&X, Array2DSize&size1, int Q, int*Ls, T*DfltValParam=NULL){//16
     Array2DSize size2;
     size2.Set(Q, Ls);
     ArrayOf2DBy1DSetSize(X, size1, size2, DfltValParam);
}
template<typename T> void ArrayOf2DBy1DSetSize(std::vector<T>&X, Array2DSize&size1, int Q, int IneRowslength, T*DfltValParam=NULL){//16
     Array2DSize size2;
     int QOld=X.size(), LcurOld, LcurNew;
     size2.Set(Q, IneRowslength);
     ArrayOf2DBy1DSetSize(X, size1, size2, DfltValParam);
}

template<typename T> void ArrayOf2DBy1DSetRect(T*&X, Array2DSize&size, int QIneRows_m1_1st_m2_last_m3_min_Other_max=-4, int QExtRows=0, T*DefaultValParam=NULL, bool PreserveVals=true){
    int Q=size.GetQExtRows(), Lmin=size.GetMinLength(), Lmax=size.GetMaxLength(), L1=size.GetLength(1), LLast=size.GetLength(Q);
    if(QIneRows_m1_1st_m2_last_m3_min_Other_max==-1){
         QIneRows_m1_1st_m2_last_m3_min_Other_max=L1;
    }else if(QIneRows_m1_1st_m2_last_m3_min_Other_max==-2){
        QIneRows_m1_1st_m2_last_m3_min_Other_max=LLast;
   }else if(QIneRows_m1_1st_m2_last_m3_min_Other_max==-3){
        QIneRows_m1_1st_m2_last_m3_min_Other_max=Lmin;
   }else if(QIneRows_m1_1st_m2_last_m3_min_Other_max<0){
        QIneRows_m1_1st_m2_last_m3_min_Other_max=Lmax;
   }
    if(QExtRows==0){
        if(X!=NULL){
            QExtRows=size.GetQExtRows();
        }
        //void ArrayOf2DBy1DSetSize(T*&X, Array2DSize&size1, int Q,    int IneRowslength,                      T*DfltValParam=NULL, bool PreserveVals=true)
        ArrayOf2DBy1DSetSize          (X,              size, QExtRows, QIneRows_m1_1st_m2_last_m3_min_Other_max, DefaultValParam,        PreserveVals);
    }
}
template<typename T> void ArrayOf2DBy1DSetRect(std::vector<T>&X, Array2DSize&size, int QIneRows, int QExtRows=0, T*DefaultValParam=NULL){
    //void ArrayOf2DBy1DSetSize(std::vector<T>&X,  const Array2DSize&size, int Q,  int IneRowslength, T*DfltValParam=NULL)
    if(QExtRows==0 && QIneRows>0){
        QExtRows=size.GetQExtRows();
    }
    ArrayOf2DBy1DSetSize(                      X,                    size, QExtRows,   QIneRows,   DefaultValParam);
}
template<typename T> void ArrayOf2DBy1DAssign(std::vector<std::vector<T>>&X, int L, T*data=NULL, bool QExtRowsNotIneRowLength=true, T*DefaultValParam=NULL){
    T DefaultVal, CurVal;
    std::vector<T>row;
    if(DefaultValParam!=NULL){
        DefaultVal=(*(DefaultValParam));
    }
    X.clear();
    if(L>0){
        if(QExtRowsNotIneRowLength){
            if(row!=NULL){
                for(int i=1; i<=L; i++){
                    row.clear();
                    CurVal=data[i-1];
                    row.push_back(CurVal);
                    X.push_back(row);
                }
            }else{
                for(int i=1; i<=L; i++){
                    row.clear();
                    CurVal=DefaultVal;
                    row.push_back(CurVal);
                    X.push_back(row);
                }
            }
        }else{
            if(row!=NULL){
                row.clear();
                for(int i=1; i<=L; i++){
                    CurVal=data[i-1];
                    row.push_back(CurVal);
                }
            }else{
                for(int i=1; i<=L; i++){
                    CurVal=DefaultVal;
                    row.push_back(CurVal);
                }
            }
            X.push_back(row);
        }
    }
}

//SetExtRow

template<typename T> void ArrayOf2DBy1DSetExtRowN(T**&X, Array2DSize&size, int ExtRowN, int wholeRowL=0, T*rowParam=NULL, int whatFromN=1, int QDefaultValuesBefore=0, bool KeepIfRect=true, T*DefaultValParam=NULL, bool DefaultIsSpecNotOwn=false){
    std::vector<T>rowVect;
    int preN1_, preN2_, ownN1_, ownN2, ForOwnN1_, ForOwnN2, postN1, postN2_,
        //Lold=wholeRowL,
        LreqGiven, LreqCalcd,
        //FromN,
        //DefaultValuesBefore
        //;
        minL, Q=size.GetQExtRows(), Lmin=size.GetMinLength(), Lmax=size.GetMaxLength(), LERN;
    T DefaultVal, CurVal;
    if(ExtRowN<0){
        ExtRowN=Q+ExtRowN+1;
    }
    if(ExtRowN>=1 && ExtRowN<=Q){
        LERN=size.GetLength(ExtRowN);
        if(wholeRowL==0){
            wholeRowL=LERN;
        }
        if(DefaultValParam!=NULL){
            DefaultVal=(*(DefaultValParam));
        }
        if(size.isRectangular() && KeepIfRect){
            //LreqGiven=size.GetIneRowsLength();//os irr!
            LreqGiven=Lmin;//or max, D s'egal
        }else{
            LreqGiven=0;
        }
        //void Calc_Array1DSubArray_byLs_MarkupNs_vFmls(int&LreqCalcd, int&preN1_, int&preN2_, int&ownN1_, int&ownN2, int&ForOwnN1_, int&ForOwnN2, int&postN1, int&postN2_, int Lini, int LreqGiven, int whatFromN, int QDefaultValuesBefore){
        //void Calc_Array1DSubArray_byLs_MarkupNs_vFmls(int&LreqCalcd, int&preN1_, int&preN2_, int&ownN1_, int&ownN2, int&ForOwnN1_, int&ForOwnN2, int&postN1, int&postN2_, int Lini, int LreqGiven=0, int whatFromN=1, int QDefaultValuesBefore=0);
        Calc_Array1DSubArray_byLs_MarkupNs_vFmls( LreqCalcd,     preN1_,    preN2_,     ownN1_,     ownN2,      ForOwnN1_,    ForOwnN2,      postN1,     postN2_, wholeRowL,       LreqGiven,       whatFromN,       QDefaultValuesBefore);
        //pre part
        if(preN2_>0){
            if(DefaultIsSpecNotOwn){
                minL = wholeRowL<=preN2_ ? wholeRowL : preN2_; //ab ini s'ver
                for(int i=preN1_; i<=minL; i++){
                    CurVal=X[ExtRowN-1][i-1];
                    rowVect.push_back(CurVal);
                }
                //for(int i=minL+1; i<=preN2_; i++){
                //    rowVect.push_back(DefaultVal);
                //}
            }//else{
            //    for(int i=preN1_; i<=preN2_; i++){
            //        rowVect.push_back(DefaultVal);
            //    }
            //}
        }
        //self part
        if(ownN2>0){
            for(int i=ownN1_; i<=ownN2; i++){
                CurVal=rowParam[i-1];
                rowVect.push_back(CurVal);
            }
        }
        //post part
        if(postN2_>0){
            minL = wholeRowL <= postN2_ ? wholeRowL : postN2_; // I men et hin S'ver
            if(DefaultIsSpecNotOwn){
                for(int i=postN1; i<=minL; i++){
                    CurVal=rowParam[i-1];
                    rowVect.push_back(CurVal);
                }
                //for(int i=minL+1; i<=postN2_; i++){
                //    rowVect.push_back(DefaultVal);
                //}
            }//else{//else row ha dflt vals
            //    for(int i=postN1; i<=postN2_; i++){
            //        rowVect.push_back(DefaultVal);
            //    }
            //}
            //for(int i=postN1; i<=postN2_; i++){
            //    if(DefaultIsSpecNotOwn && i<=wholeRowL){
            //        CurVal=rowParam[i-1];
            //    }else{
            //        CurVal=DefaultVal;
            //    }
            //    rowVect.push_back(CurVal);
            //}
        }
        //assign
        delete[]X[ExtRowN-1];
        X[ExtRowN-1]=new T[LreqCalcd];
        for(int i=1; i<=LreqCalcd; i++){
            X[ExtRowN-1][i-1]=rowVect[i-1];
        }
        //reflect
        size.SetLength(LreqCalcd, ExtRowN);
    }
}
template<typename T> void ArrayOf2DBy1DSetExtRowN(T**&X, Array2DSize&size, int ExtRowN, std::vector<T>rowParam, T*DefaultValParam=NULL, int FromN=1, int QDefaultValuesBefore=0, bool KeepIfRect=true, bool DefaultIsSpecNotOwn=false){
    //void ArrayOf2DBy1DSetExtRowN(T**&X, Array2DSize&size, int ExtRowN, int wholeRowL=0, T*rowParam=NULL, int whatFromN=1, int QDefaultValuesBefore=0, bool KeepIfRect=true, T*DefaultValParam=NULL, bool DefaultIsSpecNotOwn=false){
     ArrayOf2DBy1DSetExtRowN          (X,             size,     ExtRowN, rowParam.size(), rowParam.data(),         FromN,       QDefaultValuesBefore,        KeepIfRect,        DefaultValParam,           DefaultIsSpecNotOwn);
    //ArrayOf2DBy1DSetExtRowN(X, size, ExtRowN, rowParam, DefaultValParam, FromN, QDefaultValuesBefore, KeepIfRect);

}
template<typename T> void ArrayOf2DBy1DSetExtRowN(std::vector<std::vector<T>>&X, int ExtRowN, std::vector<T>rowParam, int whatFromN=1, int QDefaultValuesBefore=0, bool KeepIfRect=true, T*DefaultValParam=NULL, bool DefaultIsSpecNotOwn=false){//, int ifDiff_ByArr0ByNew1Min2Max3=0){
    std::vector<T>rowVect;
    Array2DSize size=GetSizeOf2DVector(X);
    int preN1_, preN2_, ownN1_, ownN2, ForOwnN1_, ForOwnN2, postN1, postN2_,
        //Lold=wholeRowL,
        LreqGiven, LreqCalcd,
        //FromN,
        //DefaultValuesBefore
        //;
        minL, Q=X.size(), LERN, Lmin=size.GetMinLength(), Lmax=size.GetMaxLength(), wholeRowL=rowParam.size();
    T DefaultVal, CurVal;
    if(ExtRowN<0){
        ExtRowN=Q+ExtRowN+1;
    }
    if(ExtRowN>=1 && ExtRowN<=Q){
        LERN=X[ExtRowN-1].size();
        //if(wholeRowL==0){
        //    wholeRowL=LERN;
        //}
        if(DefaultValParam!=NULL){
            DefaultVal=(*(DefaultValParam));
        }
        if(Lmin==Lmax && KeepIfRect){
            LreqGiven=Lmin;
        }else{
            LreqGiven=0;
        }
        Calc_Array1DSubArray_byLs_MarkupNs_vFmls( LreqCalcd,     preN1_,    preN2_,     ownN1_,     ownN2,      ForOwnN1_,    ForOwnN2,      postN1,     postN2_, wholeRowL,       LreqGiven,       whatFromN,       QDefaultValuesBefore);
        //pre part
        if(preN2_>0){
            if(DefaultIsSpecNotOwn){
                minL = wholeRowL<=preN2_ ? wholeRowL : preN2_;
                for(int i=preN1_; i<=minL; i++){
                    CurVal=X[ExtRowN-1][i-1];
                    rowVect.push_back(CurVal);
                }
                //for(int i=minL+1; i<=preN2_; i++){
                //    rowVect.push_back(DefaultVal);
                //}
            }//else{
            //    for(int i=preN1_; i<=preN2_; i++){
            //        rowVect.push_back(DefaultVal);
            //     }
            //}
        }
        //self part
        if(ownN2>0){
            for(int i=ownN1_; i<=ownN2; i++){
                CurVal=rowParam[i-1];
                rowVect.push_back(CurVal);
            }
        }
        //post part
        if(postN2_>0){
            minL = wholeRowL <= postN2_ ? wholeRowL : postN2_;
            if(DefaultIsSpecNotOwn){
                for(int i=postN1; i<=minL; i++){
                    CurVal=rowParam[i-1];
                    rowVect.push_back(CurVal);
                }
                //for(int i=minL+1; i<=postN2_; i++){
                //    rowVect.push_back(DefaultVal);
                //}
            }//else{
            //    for(int i=postN1; i<=postN2_; i++){
            //        rowVect.push_back(DefaultVal);
            //    }
            //}
            //for(int i=postN1; i<=postN2_; i++){
            //    if(DefaultIsSpecNotOwn && i<=wholeRowL){
            //        CurVal=rowParam[i-1];
            //    }else{
            //        CurVal=DefaultVal;
            //    }
            //    rowVect.push_back(CurVal);
            //}
        }
        //assign
        //delete[]X[ExtRowN-1];X[ExtRowN-1]=new T[LreqCalcd];
        //X.clear();
        //for(int i=1; i<=LreqCalcd; i++){
        //    CurVal=rowVect[i-1];
        //    X.push_back(CurVal);
        //}
        X[ExtRowN-1].clear();
        X[ExtRowN-1]=rowVect;
        //reflect
        //size.SetLength(LreqCalcd, ExtRowN);
    }
}
template<typename T> void ArrayOf2DBy1DSetExtRowN(std::vector<std::vector<T>>&X, int ExtRowN, int wholeRowL=0, T*rowParam=NULL, int FromN=1, int QDefaultValuesBefore=0, bool KeepIfRect=true, T*DefaultValParam=NULL, bool DefaultIsSpecNotOwn=false){//, int ifDiff_ByArr0ByNew1Min2Max3=0){
    //template <typename T> std::vector<T> Array1DAssign_VbyP(T*y=NULL, int Qy=0,  T*DfltValParam=NULL){
    std::vector<T>row=Array1DAssign_VbyP(wholeRowL, rowParam,  DefaultValParam);
    //ArrayOf2DBy1DSetExtRowN(std::vector<std::vector<T>>&X, int ExtRowN, std::vector<T>rowParam, int FromN=1, int DefaultValuesBefore=0, bool KeepIfRect=true, T*DefaultValParam=NULL, bool DefaultIsSpecNotOwn=false){
    ArrayOf2DBy1DSetExtRowN(                              X,     ExtRowN,                     row,      FromN,      QDefaultValuesBefore,           KeepIfRect,        DefaultValParam,      DefaultIsSpecNotOwn);
}
//
template<typename T> void ArrayOf2DBy1DSetExtRowN(T**&X, Array2DSize&size, int ExtRowN, int Lreq, int whatFromN=1, int QDefaultValuesBefore=0, T*DefaultValParam=NULL, bool KeepIfRect=true){
    std::vector<T>rowVect;
    int Q=size.GetQExtRows(), curL=size.GetLength(ExtRowN), minL, Lmin, Lmax;
    if(ExtRowN>=1 && ExtRowN<=Q){
        rowVect=Array1DGetSubArray_byLs_PreMark(X[ExtRowN-1], Lreq, whatFromN, QDefaultValuesBefore, DefaultValParam);
        //template<typename T> void ArrayOf2DBy1DSetExtRowN(T**&X, Array2DSize&size, int ExtRowN, std::vector<T>rowParam, T*DefaultValParam=NULL, int FromN=1, int QDefaultValuesBefore=0, bool KeepIfRect=true, bool DefaultIsSpecNotOwn=false){
        ArrayOf2DBy1DSetExtRowN(X, size, ExtRowN, rowVect, DefaultValParam, whatFromN, QDefaultValuesBefore, false, false);
        Lmin=size.GetMinLength(), Lmax=size.GetMaxLength();
        if(Lmin==Lmax && KeepIfRect){
            size.Set(Q, Lmin);
        }
    }
}
template<typename T> void ArrayOf2DBy1DSetExtRowN(std::vector<std::vector<T>>&X, int ExtRowN, int Lreq, int whatFromN=1, int QDefaultValuesBefore=0, T*DefaultValParam=NULL, bool KeepIfRect=true){
    std::vector<T>rowVect;
    Array2DSize size=GetSizeOf2DVector(X);
    int Q=size.GetQExtRows(), curL=size.GetLength(ExtRowN), minL, Lmin, Lmax;
    if(ExtRowN>=1 && ExtRowN<=Q){
        rowVect=Array1DGetSubArray_byLs_PreMark(X[ExtRowN-1], Lreq, whatFromN, QDefaultValuesBefore, DefaultValParam);
        //template<typename T> void ArrayOf2DBy1DSetExtRowN(std::vector<std::vector<T>>&X, int ExtRowN, std::vector<T>rowParam, int whatFromN=1, int QDefaultValuesBefore=0, bool KeepIfRect=true, T*DefaultValParam=NULL, bool DefaultIsSpecNotOwn=false){//, int ifDiff_ByArr0ByNew1Min2Max3=0){
        ArrayOf2DBy1DSetExtRowN(X, ExtRowN, rowVect, whatFromN, QDefaultValuesBefore, false, DefaultValParam, false);
    }
}




//
template<typename T> void ArrayOf2DBy1DSwapVals(T**&X, const Array2DSize& size, int ExtRowN1, int IneRowN1, int ExtRowN2, int IneRowN2){
    int Q=size.GetQExtRows(), L1, L2, Lmin=size.GetMinLength(), Lmax=size.GetMaxLength();
    T bufVal;
    if(ExtRowN1<0){
        ExtRowN1=Q+ExtRowN1+1;
    }
    if(ExtRowN2<0){
        ExtRowN2=Q+ExtRowN2+1;
    }
    if(ExtRowN1>=1 && ExtRowN1<=Q && ExtRowN2>=1 && ExtRowN2<=Q){
        L1=size.GetLength(ExtRowN1);
        L2=size.GetLength(ExtRowN2);
        if(IneRowN1<0){
            IneRowN1=L1+IneRowN1+1;
            IneRowN1=Lmax+IneRowN1+1;
        }
        if(IneRowN2<0){
            IneRowN2=L2+IneRowN2+1;
            IneRowN2=Lmax+IneRowN2+1;
        }
        if(IneRowN1>=1 && IneRowN1<=L1 && IneRowN1>=1 && IneRowN1<=L2){
            bufVal=X[ExtRowN1-1][IneRowN1-1];
            X[ExtRowN1-1][IneRowN1-1]=X[ExtRowN2-1][IneRowN2-1];
            X[ExtRowN2-1][IneRowN2-1]=bufVal;
        }
    }
}
/*template<typename T>ArrayOf2DBy1DSwapVals(std::vector<std::vector<T>>&X,int ExtRowN1, int IneRowN1, int ExtRowN2, int IneRowN2){
    Array2DSize size=GetSizeOf2DVector(X);
    int Q=size.GetQExtRows(), L1, L2, Lmin=size.GetMinLength(), Lmax=size.GetMaxLength();
    T bufVal;
    if(ExtRowN1<0){
        ExtRowN1=Q+ExtRowN1+1;
    }
    if(ExtRowN2<0){
        ExtRowN2=Q+ExtRowN2+1;
    }
    if(ExtRowN1>=1 && ExtRowN1<=Q && ExtRowN2>=1 && ExtRowN2<=Q){
        L1=size.GetLength(ExtRowN1);
        L2=size.GetLength(ExtRowN2);
        if(IneRowN1<0){
            IneRowN1=L1+IneRowN1+1;
            IneRowN1=Lmax+IneRowN1+1;
        }
        if(IneRowN2<0){
            IneRowN2=L2+IneRowN2+1;
            IneRowN2=Lmax+IneRowN2+1;
        }
        if(IneRowN1>=1 && IneRowN1<=L1 && IneRowN1>=1 && IneRowN1<=L2){
            bufVal=X[ExtRowN1-1][IneRowN1-1];
            X[ExtRowN1-1][IneRowN1-1]=X[ExtRowN2-1][IneRowN2-1];
            X[ExtRowN2-1][IneRowN2-1]=bufVal;
        }
    }
}*/
template<typename T>ArrayOf2DBy1DSwapVals(std::vector<std::vector<T>>&X, int ExtRowN1, int IneRowN1, int ExtRowN2, int IneRowN2){
    Array2DSize size;
    size=GetSizeOf2DVector(X);
    int Q=size.GetQExtRows(), L1, L2, Lmin=size.GetMinLength(), Lmax=size.GetMaxLength();
    T bufVal;
    if(ExtRowN1<0){
        ExtRowN1=Q+ExtRowN1+1;
    }
    if(ExtRowN2<0){
        ExtRowN2=Q+ExtRowN2+1;
    }
    if(ExtRowN1>=1 && ExtRowN1<=Q && ExtRowN2>=1 && ExtRowN2<=Q){
        L1=size.GetLength(ExtRowN1);
        L2=size.GetLength(ExtRowN2);
        if(IneRowN1<0){
            IneRowN1=L1+IneRowN1+1;
            IneRowN1=Lmax+IneRowN1+1;
        }
        if(IneRowN2<0){
            IneRowN2=L2+IneRowN2+1;
            IneRowN2=Lmax+IneRowN2+1;
        }
        if(IneRowN1>=1 && IneRowN1<=L1 && IneRowN1>=1 && IneRowN1<=L2){
            bufVal=X[ExtRowN1-1][IneRowN1-1];
            X[ExtRowN1-1][IneRowN1-1]=X[ExtRowN2-1][IneRowN2-1];
            X[ExtRowN2-1][IneRowN2-1]=bufVal;
        }
    }
}
//
template<typename T>ArrayOf2DBy1DSwapExtRows(T**&X, Array2DSize& size, int ExtRowN1, int ExtRowN2){
    T*bufPtr=NULL, *ptr1=NULL, *ptr2=NULL;
    int Q=size.GetQExtRows(), L1, L2, Lmin=size.GetMinLength(), Lmax=size.GetMaxLength();
    T bufVal;
    if(ExtRowN1<0){
        ExtRowN1=Q+ExtRowN1+1;
    }
    if(ExtRowN2<0){
        ExtRowN2=Q+ExtRowN2+1;
    }
    if(ExtRowN1>=1 && ExtRowN1<=Q && ExtRowN2>=1 && ExtRowN2<=Q){
        bufPtr=X[ExtRowN1-1];
        X[ExtRowN1-1]=X[ExtRowN2-1];
        X[ExtRowN2-1]=bufPtr;
        //
        //all id below s' uz try et debug, ob hic vrn narb'te, ma nu arb
        //
        //bufPtr=X[ExtRowN1-1];
        //X[ExtRowN1-1]=NULL;
        //X[ExtRowN1-1]=X[ExtRowN2-1];
        //X[ExtRowN2-1]=NULL;
        //X[ExtRowN2-1]=bufPtr;
        //
        //ptr1=X[ExtRowN1-1];
        //ptr2=X[ExtRowN2-1];
        //bufPtr=ptr1;
        //ptr1=ptr2;
        //ptr2=bufPtr;
        //X[ExtRowN1-1]=ptr1;
        //X[ExtRowN2-1]=ptr2;
        //
        //ptr1=X[ExtRowN1-1];
        //ptr2=X[ExtRowN2-1];
        //bufPtr=ptr1;
        //ptr1=NULL;
        //ptr1=ptr2;
        //ptr2=NULL;
        //ptr2=bufPtr;
        //X[ExtRowN1-1]=NULL;
        //X[ExtRowN2-1]=NULL;
        //X[ExtRowN1-1]=ptr1;
        //X[ExtRowN2-1]=ptr2;
        //
        L1=size.GetLength(ExtRowN1);
        L2=size.GetLength(ExtRowN2);
        size.SetLength(L1, ExtRowN2);
        size.SetLength(L2, ExtRowN1);
    }
}
template<typename T>ArrayOf2DBy1DSwapExtRows(std::vector<std::vector<T>>&X, int ExtRowN1, int ExtRowN2){
    std::vector<T>bufRow;
    Array2DSize size;
    size=GetSizeOf2DVector(X);
    int Q=size.GetQExtRows(), L1, L2, Lmin=size.GetMinLength(), Lmax=size.GetMaxLength();
    T bufVal;
    if(ExtRowN1<0){
        ExtRowN1=Q+ExtRowN1+1;
    }
    if(ExtRowN2<0){
        ExtRowN2=Q+ExtRowN2+1;
    }
    if(ExtRowN1>=1 && ExtRowN1<=Q && ExtRowN2>=1 && ExtRowN2<=Q){
        bufRow=X[ExtRowN1-1];
        X[ExtRowN1-1]=X[ExtRowN2-1];
        X[ExtRowN2-1]=bufRow;
    }
}
//
template<typename T>ArrayOf2DBy1DSwapIneRows(T**&X, const Array2DSize& size, int IneRowN1, int IneRowN2){
    int Q=size.GetQExtRows(), L1, L2, Lmin, Lmax;
    if(Q>0){
        Lmin=size.GetMinLength();
        Lmax=size.GetMaxLength();
        if(IneRowN1<0){
            IneRowN1=Lmax+IneRowN1+1;
        }
        if(IneRowN2<0){
            IneRowN2=Lmax+IneRowN2+1;
        }
        if(IneRowN1>=1 && IneRowN1<=Lmin && IneRowN2>=1 && IneRowN2<=Lmin){
            for(int i=1; i<=Q; i++){
                ArrayOf2DBy1DSwapVals(X, size, i, IneRowN1, i, IneRowN2);
            }
        }
    }
}
template<typename T>ArrayOf2DBy1DSwapIneRows(std::vector<std::vector<T>>&X, int IneRowN1, int IneRowN2){
    Array2DSize size;
    size=GetSizeOf2DVector(X);
    int Q=size.GetQExtRows(), L1, L2, Lmin, Lmax;
    if(Q>0){
        Lmin=size.GetMinLength();
        Lmax=size.GetMaxLength();
        if(IneRowN1<0){
            IneRowN1=Lmax+IneRowN1+1;
        }
        if(IneRowN2<0){
            IneRowN2=Lmax+IneRowN2+1;
        }
        if(IneRowN1>=1 && IneRowN1<=Lmin && IneRowN2>=1 && IneRowN2<=Lmin){
            for(int i=1; i<=Q; i++){
                ArrayOf2DBy1DSwapVals(X, size, i, IneRowN1, i, IneRowN2);
            }
        }
    }
}
//Reverse ext rows
template<typename T>ArrayOf2DBy1DReverseExtRows(T**&X, Array2DSize& size){
    int Q=size.GetQExtRows(), N1, N2, Nmid;
    if(Q%2==0){
        Nmid=Q/2;
    }else{
        Nmid=(Q-1)/2;
    }
    for(N1=1; N1<=Nmid; N1++){
        N2=Q-N1+1;
        ArrayOf2DBy1DSwapExtRows(X, size, N1, N2);
    }
}
template<typename T>ArrayOf2DBy1DReverseExtRows(std::vector<std::vector<T>>&X){
    Array2DSize size;
    size=GetSizeOf2DVector(X);
    //int Q=size.GetQExtRows(), L1, L2, Lmin=size.GetMinLength(), Lmax=size.GetMaxLength();
    int Q=size.GetQExtRows(), N1, N2, Nmid;
    if(Q%2==0){
        Nmid=Q/2;
    }else{
        Nmid=(Q-1)/2;
    }
    for(N1=1; N1<Nmid; N1++){
        N2=Q-N1+1;
        ArrayOf2DBy1DSwapExtRows(X, N1, N2);
    }
}
//
template<typename T>ArrayOf2DBy1DReverseIneRows(T**&X, const Array2DSize& size){
    int Q=size.GetQExtRows(), N1, N2, Nmid, Lmax=size.GetMaxLength(), Lmin=size.GetMinLength();
    if(Lmin==Lmax){
        if(Lmax%2==0){
            Nmid=Q/2;
        }else{
            Nmid=(Lmax-1)/2;
        }
        for(N1=1; N1<Nmid; N1++){
            N2=Lmax-N1+1;
            for(int j=1; j<=Q; j++){
                ArrayOf2DBy1DSwapVals(X, size, j, N1, j, N2);
            }
        }
    }
}
template<typename T>ArrayOf2DBy1DReverseIneRows(std::vector<std::vector<T>>&X){
    Array2DSize size;
    size=GetSizeOf2DVector(X);
    int Q=size.GetQExtRows(), N1, N2, Nmid, Lmax=size.GetMaxLength(), Lmin=size.GetMinLength();
    if(Lmin==Lmax){
        if(Lmax%2==0){
            Nmid=Q/2;
        }else{
            Nmid=(Lmax-1)/2;
        }
        for(N1=1; N1<Nmid; N1++){
            N2=Lmax-N1+1;
            for(int j=1; j<=Q; j++){
                ArrayOf2DBy1DSwapVals(X, size, j, N1, j, N2);
            }
        }
    }
}
//
template<typename T>ArrayOf2DBy1DReverseExtRowN(T**&X, const Array2DSize& size, int ExtRowN){
    int Lcur;
    if(ExtRowN>=1 && ExtRowN<=size.GetQExtRows()){
        Lcur=size.GetLength(ExtRowN);
        Array1DReverse(X[ExtRowN-1], Lcur);
    }
}
template<typename T>ArrayOf2DBy1DReverseExtRowN(std::vector<std::vector<T>>&X, int ExtRowN){
    int Lcur;
    if(ExtRowN>=1 && ExtRowN<=X.size()){
        Lcur=X[ExtRowN-1].size();
        Array1DReverse(X[ExtRowN-1]);
    }
}
template<typename T>ArrayOf2DBy1DReverseIneRowN(T**&X, const Array2DSize& size, int IneRowN){
    int Q=size.GetQExtRows(), Lmin=size.GetMinLength(), Lmax=size.GetMaxLength(), Nmid, N1, N2;
    if(Q>0 && IneRowN>=1 && IneRowN<=Lmin){
        if(Q%2==0){
            Nmid=Q/2;
        }else{
            Nmid=(Q-1)/2;
        }
        for(N1=1; N1<=Nmid; N1++){
            N2=Q-N1+1;
            ArrayOf2DBy1DSwapVals(X, size, N1, IneRowN, N2, IneRowN);
        }
    }
}
template<typename T>ArrayOf2DBy1DReverseIneRowN(std::vector<std::vector<T>>&X, int IneRowN){
    Array2DSize size;
    size=GetSizeOf2DVector(X);
    int Q=size.GetQExtRows(), Lmin=size.GetMinLength(), Lmax=size.GetMaxLength(), Nmid, N1, N2;
    if(Q>0 && IneRowN>=1 && IneRowN<=Lmin){
        if(Q%2==0){
            Nmid=Q/2;
        }else{
            Nmid=(Q-1)/2;
        }
        for(N1=1; N1<=Nmid; N1++){
            N2=Q-N1+1;
            ArrayOf2DBy1DSwapVals(X, N1, IneRowN, N2, IneRowN);
        }
    }
}
//AddToExtRow
template<typename T>ArrayOf2DBy1DAddToExtRowN(T**&X, Array2DSize& size, int ExtRowN, T val){
    int Q=size.GetQExtRows(), Lmin=size.GetMinLength(), Lmax=size.GetMaxLength(), Lcur;
    if(ExtRowN>=1 && ExtRowN<=Q){
        Lcur=size.GetLength(ExtRowN);
        Array1DAddElement(X[ExtRowN-1], Lcur, val);
        size.SetLength(Lcur, ExtRowN);
    }
}
template<typename T>ArrayOf2DBy1DAddToExtRowN(std::vector<std::vector<T>>&X, int ExtRowN, T val){
    Array2DSize size;
    size=GetSizeOf2DVector(X);
    int Q=size.GetQExtRows(), Lmin=size.GetMinLength(), Lmax=size.GetMaxLength(), Lcur;
    if(ExtRowN>=1 && ExtRowN<=Q){
        Lcur=size.GetLength(ExtRowN);
        Array1DAddElement(X[ExtRowN-1], val);
        size.SetLength(Lcur, ExtRowN);
    }
}
//InsToExtRow
template<typename T>ArrayOf2DBy1DInsToExtRowN(T**&X, Array2DSize& size, int ExtRowN, T val, int IneRowPosN){
    int Q=size.GetQExtRows(), Lmin=size.GetMinLength(), Lmax=size.GetMaxLength(), Lcur;
    if(ExtRowN>=1 && ExtRowN<=Q){
        Lcur=size.GetLength(ExtRowN);
        Array1DInsElementToN(X[ExtRowN-1], Lcur, IneRowPosN, val);
        size.SetLength(Lcur, ExtRowN);
    }
}
template<typename T>ArrayOf2DBy1DInsToExtRowN(std::vector<std::vector<T>>&X, int ExtRowN, T val, int IneRowPosN){
    Array2DSize size;
    size=GetSizeOf2DVector(X);
    int Q=size.GetQExtRows(), Lmin=size.GetMinLength(), Lmax=size.GetMaxLength(), Lcur;
    if(ExtRowN>=1 && ExtRowN<=Q){
        Lcur=size.GetLength(ExtRowN);
        Array1DInsElementToN(X[ExtRowN-1],IneRowPosN, val);
        size.SetLength(Lcur, ExtRowN);
    }
}
//DelFromExtRow
template<typename T>ArrayOf2DBy1DDelFromExtRowN(T**&X, Array2DSize& size, int ExtRowN, int IneRowPosN){
    int Q=size.GetQExtRows(), Lmin=size.GetMinLength(), Lmax=size.GetMaxLength(), Lcur;
    if(ExtRowN>=1 && ExtRowN<=Q){
        Lcur=size.GetLength(ExtRowN);
        Array1DDelElementFromN(X[ExtRowN-1], Lcur, IneRowPosN);
        size.SetLength(Lcur, ExtRowN);
    }
}
template<typename T>ArrayOf2DBy1DDelFromExtRowN(std::vector<std::vector<T>>&X, int ExtRowN, int IneRowPosN){
    Array2DSize size;
    size=GetSizeOf2DVector(X);
    int Q=size.GetQExtRows(), Lmin=size.GetMinLength(), Lmax=size.GetMaxLength(), Lcur;
    if(ExtRowN>=1 && ExtRowN<=Q){
        Lcur=size.GetLength(ExtRowN);
        Array1DDelElementFromN(X[ExtRowN-1], IneRowPosN);
        size.SetLength(Lcur, ExtRowN);
    }
}
//
//adding ext rows
template<typename T>void ArrayOf2DBy1DAddExtRow(T**&X, Array2DSize&size, T*rowParam, int wholeRowL=0, int whatFromN=1, int QDefaultValsBefore=0, bool RectNotVar=true, T*DfltValParam=NULL, TValsShowHide*vsh=NULL){
    std::vector<T> row;
    T DefaultVal, CurVal, *TmpPtr=NULL;
    int Lreq, QTmp;
    Array2DSize NewSize=size;
    if(DfltValParam!=NULL){
        DefaultVal=(*(DfltValParam));
    }
    int vrnN=0;//0-add NULL et Set, 1-add 1-member row et Set,
               //2- add needed row et reset, 3 - 2 ac reset
               //4
    if((size.GetQExtRows()==0 || X==NULL) && wholeRowL>0){//creating new on NULL base
        QTmp=0;
        Lreq=0;
        //std::vector<T> std::vector<T> Array1DGetSubArray_byLs_PreMark(T*x, int Lold, int Lreq=0, int FromN=1, int FirstDefaultValues=0, T*DfltValParam=NULL);
        row=Array1DGetSubArray_byLs_PreMark(rowParam, wholeRowL, Lreq, whatFromN, QDefaultValsBefore, DfltValParam);
        //
        Array1DAddElement(X, QTmp, row.data());
        //
        if(RectNotVar){
            size.SetVaria(1, 1);
        }else{
            size.Set(1, 1);
        }
    }else{//creating
        if(size.isRectangular() && RectNotVar){
            Lreq=size.GetIneRowsLength();
        }else{
            // int CalcLreq_Array1DSubArray_byLs_MarkupNs_vFmls(int Lini, int whatFromN=1, int QDefaultValuesBefore=0, int LreqGiven=0);
            Lreq=CalcLreq_Array1DSubArray_byLs_MarkupNs_vFmls(wholeRowL, whatFromN, QDefaultValsBefore, 0);
        }
        //std::vector<T> Array1DGetSubArray_byLs_PreMark(T*x, int Lold, int Lreq=0, int FromN=1, int FirstDefaultValues=0, T*DfltValParam=NULL)
        row=Array1DGetSubArray_byLs_PreMark(rowParam, wholeRowL, Lreq, whatFromN, QDefaultValsBefore, DfltValParam);
        QTmp=size.GetQExtRows();
        //
        //Array1DAddElement(X, QTmp, row.data());//
        //I men S narb, et ob af fn stop arb, row, ad ic ptr zeig, deda, et ptr adbe noin, ne null, ma zeig ad anid random
        //T*TmpPtr=NULL;//narb, I men ob am mub
        //Array1DAddElement(X, QTmp, NULL);//S nint if ptr = NULL?! Ja, param tam ms'be ne unes NULL, ma ptr of T type, qof val=NULL!
        //std::cout<<"Add ext row: vrnN="<<vrnN<<std::endl;
        switch(vrnN){
            case 0:
                //LetIt remain NULL
            break;
            case 1:
                TmpPtr=new T[1];
                TmpPtr[1-1]=DefaultVal;
            break;
            case 2:
            case 3:
                TmpPtr=new T[Lreq];
                for(int i=1; i<=Lreq; i++){
                    CurVal=row[i-1];
                    TmpPtr[i-1]=CurVal;
                    //TmpPtr[i-1]=row[i-1];
                }
            break;
            case 4:
            case 5:
                //NOp, see further
                //ne tested yet
            break;
            case 6:
                NewSize.AddExt(row.size());
                ArrayOf2DBy1DSetSize(X, size, NewSize, DfltValParam, true);
                //ne tested yet
            break;
        }
        if(vrnN!=6 && vrnN!=4 && vrnN!=5){
            Array1DAddElement(X, QTmp, TmpPtr);
        }
        if( vrnN==4 || vrnN==5){
             Array1DAddElement(X, QTmp, row.data());
             //ne tested yet
        }
        //it makes QTmp++
        switch(vrnN){
            case 0:
                size.AddExt(0);
            break;
            case 1:
                size.AddExt(1);
            break;
            case 2:
            case 3:
                size.AddExt(Lreq);
            break;
        }
        //;
        //std::cout<<"Temporally arr2D's last row is:"<<std::endl;
        //void Arr1DStdCOut(T*X, int L, QString delim="; ", bool useBrackets=true){
        //Arr1DStdCOut(X[QTmp-1], Lreq, " ", false);//defined later
        //for(int i=1; i<=1; i++)std::cout<<X[QTmp-1][i-1]<<" ";
        //for(int i=1; i<=Lreq; i++)std::cout<<X[QTmp-1][i-1]<<" ";
        //std::cout<<std::endl;
        //
        //Array1DAddElement(X, QTmp,  row.data());
        //template<typename T> void ArrayOf2DBy1DSetExtRowN(T**&X, Array2DSize&size, int ExtRowN, int wholeRowL=0, T*rowParam=NULL, int whatFromN=1, int QDefaultValuesBefore=0, bool KeepIfRect=true, T*DefaultValParam=NULL, bool DefaultIsSpecNotOwn=false){
        if(vrnN!=3 &&  vrnN!=4 && vrnN!=5){
            ArrayOf2DBy1DSetExtRowN(X, size, QTmp, Lreq, row.data(), whatFromN, QDefaultValsBefore, RectNotVar, DfltValParam, true);
        }
        //
        //size.AddExt(Lreq);//zu afnot!
    }
}
template<typename T>void ArrayOf2DBy1DAddExtRow(T**&X, Array2DSize&size, std::vector<T>rowParam, int whatFromN=1, int QDefaultValsBefore=0, bool RectNotVar=true, T*DfltValParam=NULL, TValsShowHide*vsh=NULL){
    ArrayOf2DBy1DAddExtRow(X, size, rowParam.data(), rowParam.size(), whatFromN, QDefaultValsBefore, RectNotVar, DfltValParam, vsh);
}
template<typename T>void ArrayOf2DBy1DAddExtRow(std::vector<std::vector<T>>&X, std::vector<T>rowParam,  int whatFromN=1, int QDefaultValsBefore=0, bool RectNotVar=true, T*DfltValParam=NULL, TValsShowHide*vsh=NULL){
    std::vector<T> row;
    T DefaultVal, CurVal;
    int wholeRowL=rowParam.size();
    int Lreq;
    int  preN1_=0, preN2_=0, ownN1_=0, ownN2=0, ForOwnN1_=0, ForOwnN2=0, postN1=0, postN2_=0;//, Lnew=0;
    Array2DSize size;
    size=GetSizeOf2DVector(X);
    if(DfltValParam!=NULL){
        DefaultVal=(*(DfltValParam));
    }
    if(size.GetQExtRows()==0  && wholeRowL>0){
        Lreq=0;
        //template <typename T>std::vector<T> Array1DGetSubArray_byLs_PreMark(std::vector<T>x,  int Lreq=0, int FromN=1, int FirstDefaultValues=0, T*DfltValParam=NULL){
        row=Array1DGetSubArray_byLs_PreMark(rowParam,  Lreq, whatFromN, QDefaultValsBefore, DfltValParam);
        Array1DAddElement(X, row);
    }else{
        if(size.isRectangular() && RectNotVar){
            Lreq=size.GetIneRowsLength();
        }else{
            //void CalcLreq_Array1DSubArray_byLs_MarkupNs_vFmls(int Lini, int whatFromN=1, int QDefaultValuesBefore=0, int LreqGiven=0);
            Lreq=CalcLreq_Array1DSubArray_byLs_MarkupNs_vFmls(rowParam.size(), whatFromN, QDefaultValsBefore, 0);
        }
        //std::vector<T> Array1DGetSubArray_byLs_PreMark(T*x, int Lold, int Lreq=0, int FromN=1, int FirstDefaultValues=0, T*DfltValParam=NULL)
        //row=Array1DGetSubArray_byLs_PreMark(rowParam, wholeRowL, Lreq, whatFromN, QDefaultValsBefore, DfltValParam);
        row=Array1DGetSubArray_byLs_PreMark(rowParam, Lreq, whatFromN, QDefaultValsBefore, DfltValParam);
        //QTmp=size.GetQExtRows();
        Array1DAddElement(X, row);
        //
        //size.AddExt(Lreq);
    }
}
template<typename T>void ArrayOf2DBy1DAddExtRow(std::vector<std::vector<T>>&X, int wholeRowL=0, T*rowParam=NULL,  int FromN=1, int QDefaultValsBefore=0, bool RectNotVar=true, T*DfltValParam=NULL, TValsShowHide*vsh=NULL){
    std::vector<T>rowIni;
    //template <typename T> std::vector<T> Array1DAssign_VbyP(T*y=NULL, int Qy=0,  T*DfltValParam=NULL){
    rowIni=Array1DAssign_VbyP(wholeRowL, rowParam, DfltValParam);
    ArrayOf2DBy1DAddExtRow(X, rowIni, FromN, QDefaultValsBefore, RectNotVar, DfltValParam, vsh);
}
//
//Del ext row
template<typename T>void ArrayOf2DBy1DDelExtRow(T**&X, Array2DSize&size, int N){//42
    int Q=size.GetQExtRows(), L;
    if(N>=1 && N<=Q){
        Array1DDelElementFromN(X, Q, N);
        size.DelExt(N);
    }
}
template<typename T>void ArrayOf2DBy1DDelExtRow(std::vector<std::vector<T>>&X, int N){
    if(N>=1 && N<=X.size()){
        Array1DDelElementFromN(X, N);//if cany overload
        //X.erase(x.begin()+N-1);
    }
}
//
//Ins ext row
template<typename T> void ArrayOf2DBy1DInsExtRow(T**&X,  Array2DSize&size, int ExtRowN, int wholeRowL=0, T*rowParam=NULL, T*DefaultValParam=NULL, int whatFromN=1, int QDefaultValuesBefore=0, bool KeepIfRect=true, TValsShowHide*vsh=NULL){
    //writeln(&vsh,"ArrayOf2DBy1DInsExtRow(ptr, ptr) starts working");
    if(vsh!=NULL)std::cout<<" ArrayOf2DBy1DInsExtRow(ptr, ptr) starts working "<<std::endl;
    int Q=size.GetQExtRows();
    if(ExtRowN<0){
        ExtRowN=Q+ExtRowN+1;
    }
    if(ExtRowN>=1 && ExtRowN<=Q){
        //void ArrayOf2DBy1DAddExtRow(T**&X, Array2DSize&size, T*rowParam, int wholeRowL=0, int whatFromN=1, int QDefaultValsBefore=0, bool RectNotVar=true, T*DfltValParam=NULL, TValsShowHide*vsh=NULL)
        ArrayOf2DBy1DAddExtRow(           X,             size,   rowParam,     wholeRowL,       whatFromN,       QDefaultValuesBefore,      KeepIfRect,        DefaultValParam);
        //ob i add is done all complex work, et for Q=0 S narb
        Q=size.GetQExtRows();
        if(vsh!=NULL){
            std::cout<<" Seems to be added "<<std::endl;
        }
        if(vsh!=NULL){
            //void Arr2DStdCOut2D(T**X, const Array2DSize size, QString delimElem="; ", QString delimLines="; ", bool useBrackets=true){
            Arr2DStdCOut2D(X, size);
        }
        for(int i=Q-1; i>=ExtRowN; i--){
            ArrayOf2DBy1DSwapExtRows(X, size, i, i+1);
            if(vsh!=NULL){
                std::cout<<"Swapping rows "<<i<<" and "<<i+1<<std::endl;
                Arr2DStdCOut2D(X, size);
            }
        }
    }
    //writeln(&vsh,"ArrayOf2DBy1DInsExtRow(ptr, ptr) finishes working");
    if(vsh!=NULL)std::cout<<" ArrayOf2DBy1DInsExtRow(ptr, ptr) finishes working "<<std::endl;
}
template<typename T> void ArrayOf2DBy1DInsExtRow(T**&X,  Array2DSize&size, int ExtRowN, std::vector<T>rowParam, T*DefaultValParam=NULL, int whatFromN=1, int QDefaultValuesBefore=0, bool KeepIfRect=true){
    ArrayOf2DBy1DInsExtRow(X, size, ExtRowN, rowParam.size(), rowParam.data(), DefaultValParam, whatFromN, QDefaultValuesBefore, KeepIfRect);
}
template<typename T> void ArrayOf2DBy1DInsExtRow(std::vector<std::vector<T>>&X, int ExtRowN, std::vector<T>rowParam, T*DefaultValParam=NULL, int whatFromN=1, int QDefaultValuesBefore=0, bool KeepIfRect=true){
    Array2DSize size;
    size=GetSizeOf2DVector(X);
    //newSize=size;
    std::vector<T>row;
    int Q=size.GetQExtRows(), Lmin=size.GetMinLength(), Lmax=size.GetMaxLength(), Lreq=0;
    if(ExtRowN<0){
        ExtRowN=Q+ExtRowN+1;
    }
    if(ExtRowN>=1 && ExtRowN<=Q){
        if(Lmin==Lmax && KeepIfRect){
            Lreq=Lmin;//or Lmax, D s'egal
        }
        row=Array1DGetSubArray_byLs_PreMark(rowParam, Lreq, whatFromN, QDefaultValuesBefore, DefaultValParam);
        Array1DInsElementToN(X, ExtRowN, row);
        //ob co vects all s'unes
    }
}
template<typename T> void ArrayOf2DBy1DInsExtRow(std::vector<std::vector<T>>&X, int ExtRowN, int wholeRowL=0, T*rowParam=NULL, T*DefaultValParam=NULL, int whatFromN=1, int QDefaultValuesBefore=0, bool KeepIfRect=true){
    Array2DSize size;//, newSize;
    size=GetSizeOf2DVector(X);
    //newSize=size;
    std::vector<T>row;
    int Q=size.GetQExtRows(), Lmin=size.GetMinLength(), Lmax=size.GetMaxLength(), Lreq=0;
    if(ExtRowN<0){
        ExtRowN=Q+ExtRowN+1;
    }
    if(ExtRowN>=1 && ExtRowN<=Q){
        if(Lmin==Lmax && KeepIfRect){
            Lreq=Lmin;//or Lmax, D s'egal
        }
        row=Array1DGetSubArray_byLs_PreMark(rowParam, wholeRowL, Lreq, whatFromN, QDefaultValuesBefore, DefaultValParam);
        Array1DInsElementToN(X, ExtRowN, row);
        //ob co vects all s'unes
    }
}
//
//SetIneRow
template<typename T>ArrayOf2DBy1DSetIneRow(T**&X, Array2DSize& size, int IneRowN, T*rowParam=NULL, int wholeRowL=0, int whatFromN=1, int QDefaultValsBefore=0, T*DefaultValParam=NULL, bool ExistingInsteadOfDefault=true, bool ExistingOnlyForGTLminNotIgnore=false, bool LastPossIfNotRectAtPos0NotIgnore=false){
    int Q=size.GetQExtRows(), Lmin=size.GetMinLength(), Lmax=size.GetMaxLength(), Lcur, n, Lreq=Q;
    std::vector<T>row;//=Array1DGetSubArray_byLs_PreMark(rowParam, wholeRowL, Lreq, whatFromN, QDefaultValsBefore, DefaultValParam);
    T DefaultVal;
    int LreqCalcd=0, preN1_=0, preN2_=0, ownN1_=0, ownN2=0, ForOwnN1_=0, ForOwnN2=0, postN1=0, postN2_=0;
    //param Ns
    if(DefaultValParam!=NULL){
        DefaultVal=(*(DefaultValParam));
    }
    if(wholeRowL==0 && rowParam!=NULL){
        wholeRowL=Lreq;
    }
    if(IneRowN<0){
        IneRowN=Lmax+IneRowN+1;
    }
    //row
    //void Calc_Array1DSubArray_byLs_MarkupNs_vFmls(int&LreqCalcd, int&preN1_, int&preN2_, int&ownN1_, int&ownN2, int&ForOwnN1_, int&ForOwnN2, int&postN1, int&postN2_, int Lini, int LreqGiven=0, int whatFromN=1, int QDefaultValuesBefore=0);
    Calc_Array1DSubArray_byLs_MarkupNs_vFmls(LreqCalcd, preN1_, preN2_, ownN1_, ownN2, ForOwnN1_, ForOwnN2, postN1, postN2_, wholeRowL, Lreq, whatFromN, QDefaultValsBefore);
    row=Array1DGetSubArray_byLs_PreMark(rowParam, wholeRowL, Lreq, whatFromN, QDefaultValsBefore, DefaultValParam);
    if(ExistingInsteadOfDefault){
        //if(preN1_>0){//arb ac ce ob i=0 et length of 0. row =0, IneRowN>Lcur
            for(int i=preN1_; i<=preN2_; i++){
                Lcur=size.GetLength(i-1);
                if(IneRowN>=1 && IneRowN<=Lcur){
                    row[i-1]=X[i-1][IneRowN-1];
                }
            }
        //}
        //if(postN1>0){//arb ac ce ob i=0 et length of 0. row =0, IneRowN>Lcur
            for(int i=postN1; i<=postN2_; i++){
                Lcur=size.GetLength(i-1);
                if(IneRowN>=1 && IneRowN<=Lcur){
                    row[i-1]=X[i-1][IneRowN-1];
                }
            }
       //}
    }
    //set ipse
    if(IneRowN>=1 && IneRowN<=Lmax){
        if(IneRowN<=Lmin){
            for(int i=1; i<=Q; i++){
                X[i-1][IneRowN-1]=row[i-1];
            }
        }else if(ExistingOnlyForGTLminNotIgnore){
            for(int i=1; i<=Q; i++){
                Lcur=size.GetLength(i);
                if(IneRowN<=Lcur)
                X[i-1][IneRowN-1]=row[i-1];
            }
        }
    }else if(IneRowN==0 && LastPossIfNotRectAtPos0NotIgnore && Lmin!=Lmax){
        for(int i=1; i<=Q; i++){
            Lcur=size.GetLength(i);
            X[i-1][Lcur-1]=row[i-1];
        }
    }
}
template<typename T>ArrayOf2DBy1DSetIneRow(T**&X, Array2DSize& size, int IneRowN,std::vector<T>rowParam, int whatFromN=1, int QDefaultValsBefore=0, T*DefaultValParam=NULL, bool ExistingInsteadOfDefault=true, bool ExistingOnlyForGTLminNotIgnore=false, bool LastPossIfNotRectAtPos0NotIgnore=false){
    ArrayOf2DBy1DSetIneRow(X, size, IneRowN, rowParam.data(), rowParam.size(), whatFromN, QDefaultValsBefore, DefaultValParam, ExistingInsteadOfDefault, ExistingOnlyForGTLminNotIgnore, LastPossIfNotRectAtPos0NotIgnore);
}
template<typename T>ArrayOf2DBy1DSetIneRow(std::vector<std::vector<T>>&X, int IneRowN, std::vector<T>rowParam, int whatFromN=1, int QDefaultValsBefore=0, T*DefaultValParam=NULL, bool ExistingInsteadOfDefault=true, bool ExistingOnlyForGTLminNotIgnore=false, bool LastPossIfNotRectAtPos0NotIgnore=false){
    Array2DSize size;
    size=GetSizeOf2DVector(X);
    int wholeRowL=rowParam.size();
    int Q=size.GetQExtRows(), Lmin=size.GetMinLength(), Lmax=size.GetMaxLength(), Lcur, n, Lreq=Q;
    std::vector<T>row;//=Array1DGetSubArray_byLs_PreMark(rowParam, wholeRowL, Lreq, whatFromN, QDefaultValsBefore, DefaultValParam);
    T DefaultVal;
    int LreqCalcd=0, preN1_=0, preN2_=0, ownN1_=0, ownN2=0, ForOwnN1_=0, ForOwnN2=0, postN1=0, postN2_=0;
    //param Ns
    if(DefaultValParam!=NULL){
        DefaultVal=(*(DefaultValParam));
    }
    //if(wholeRowL==0 && rowParam!=NULL){
    //   wholeRowL=Lreq;
    //}
    if(IneRowN<0){
        IneRowN=Lmax+IneRowN+1;
    }
    //row
    //void Calc_Array1DSubArray_byLs_MarkupNs_vFmls(int&LreqCalcd, int&preN1_, int&preN2_, int&ownN1_, int&ownN2, int&ForOwnN1_, int&ForOwnN2, int&postN1, int&postN2_, int Lini, int LreqGiven=0, int whatFromN=1, int QDefaultValuesBefore=0);
    Calc_Array1DSubArray_byLs_MarkupNs_vFmls(LreqCalcd, preN1_, preN2_, ownN1_, ownN2, ForOwnN1_, ForOwnN2, postN1, postN2_, wholeRowL, Lreq, whatFromN, QDefaultValsBefore);
    row=Array1DGetSubArray_byLs_PreMark(rowParam, Lreq, whatFromN, QDefaultValsBefore, DefaultValParam);
    if(ExistingInsteadOfDefault){
        if(preN2_>0){
            for(int i=preN1_; i<=preN2_; i++){
                Lcur=size.GetLength(i-1);
                if(IneRowN>=1 && IneRowN<=Lcur){
                    row[i-1]=X[i-1][IneRowN-1];
                }
            }
        }
        if(postN1>0){
            for(int i=postN1; i<=postN2_; i++){
                Lcur=size.GetLength(i-1);
                if(IneRowN>=1 && IneRowN<=Lcur){
                    row[i-1]=X[i-1][IneRowN-1];
                }
            }
        }
    }
    //set ipse
    if(IneRowN>=1 && IneRowN<=Lmax){
        if(IneRowN<=Lmin){
            for(int i=1; i<=Q; i++){
                X[i-1][IneRowN-1]=row[i-1];
            }
        }else if(ExistingOnlyForGTLminNotIgnore){
            for(int i=1; i<=Q; i++){
                Lcur=size.GetLength(i);
                if(IneRowN<=Lcur)
                X[i-1][IneRowN-1]=row[i-1];
            }
        }
    }else if(IneRowN==0 && LastPossIfNotRectAtPos0NotIgnore && Lmin!=Lmax){
        for(int i=1; i<=Q; i++){
            Lcur=size.GetLength(i);
            X[i-1][Lcur-1]=row[i-1];
        }
    }
}
template<typename T>ArrayOf2DBy1DSetIneRow(std::vector<std::vector<T>>&X, int IneRowN, T*rowParam, int wholeRowL=0, int whatFromN=1, int QDefaultValsBefore=0, T*DefaultValParam=NULL, bool ExistingInsteadOfDefault=true, bool ExistingOnlyForGTLminNotIgnore=false, bool LastPossIfNotRectAtPos0NotIgnore=false){
    int Q=X.size();
    std::vector<T>row;
    //template <typename T> void Array1DAssign(std::vector<T>&vec, T y[], int L=0)
    Array1DAssign(row, rowParam, wholeRowL);
    //template<typename T>ArrayOf2DBy1DSetIneRow(std::vector<std::vector<T>>&X, int IneRowN, std::vector<T>rowParam, int whatFromN=1, int QDefaultValsBefore=0, T*DefaultValParam=NULL, bool ExistingInsteadOfDefault=true, bool ExistingOnlyForGTLminNotIgnore=false, bool LastPossIfNotRectAtPos0NotIgnore=false){
    ArrayOf2DBy1DSetIneRow(                                                  X,     IneRowN,                    row,     whatFromN,        QDefaultValsBefore,    DefaultValParam,           ExistingInsteadOfDefault,           ExistingOnlyForGTLminNotIgnore,            LastPossIfNotRectAtPos0NotIgnore);
}
//
//AddIneRow
template<typename T>ArrayOf2DBy1DAddIneRow(T**&X, Array2DSize& size, T*rowParam=NULL, int wholeRowL=0, int IfNonRectIgnore0Add1Stretch2=0, int whatFromN=1, int QDefaultValsBefore=0, T*DefaultValParam=NULL, bool RectNotVariaIfFirst=true){
    std::vector<T>row;
    T DefaultVal;
    int Q=size.GetQExtRows(), Lmin=size.GetMinLength(), Lmax=size.GetMaxLength(), Lcur, n;
    int Lreq;
    int *Ls=NULL;
    if(DefaultValParam!=NULL){
        DefaultVal=(*(DefaultValParam));
    }
    if(wholeRowL==0 && rowParam!=NULL){
        wholeRowL=Q;
    }
    if(Q==0 && wholeRowL>0){
        Lreq=0;
        row=Array1DGetSubArray_byLs_PreMark(rowParam, wholeRowL, Lreq, whatFromN, QDefaultValsBefore, DefaultValParam);
        Lcur=row.size();//in other if brunch Lcur means different thing! n not used
        //ob not from first et co Defaults before
        X=new T*[Lcur];
        for(int i=1; i<=Lcur; i++){
            X[i-1]=new T[1];
        }
        if(rowParam!=NULL){
            for(int i=1; i<=Lcur; i++){
                X[i-1][1-1]=row[i-1];
            }
        }else{
            for(int i=1; i<=Lcur; i++){
                X[i-1][1-1]=DefaultVal;
            }
        }
        //size:
        if(RectNotVariaIfFirst){
            size.Set(Lcur, 1);
        }else{
            Ls=new int[Lcur];
            for(int i=1; i<=Lcur; i++){
                Ls[i-1]=1;
            }
            size.Set(Lcur, Ls);
            delete[]Ls;
            Ls=NULL;
        }
    }else if(Q>0){
        Lmin=size.GetMinLength();
        Lmax=size.GetMaxLength();
        Lreq=Q;
        //row=Array1DGetSubArray_byLs_PreMark(rowParam, wholeRowL, Q, whatFromN, QDefaultValuesBefore, DefaultValParam);
        //ob ine row has length Q
        if(Lmin!=Lmax && IfNonRectIgnore0Add1Stretch2==2){
            for(int i=1; i<=Q; i++){
                Lcur=size.GetLength(i);
                n=Lcur;
                for(int j=Lcur+1; j<=Lmax; j++){
                    Array1DAddElement(X[i-1], n, DefaultVal);
                }
            }
            for(int i=1; i<=Q; i++){
                size.SetLength(Lmax, i);
            }
        }
        if(!((Lmin!=Lmax)&&(IfNonRectIgnore0Add1Stretch2==0))){
            if(rowParam!=NULL){
                row=Array1DGetSubArray_byLs_PreMark(rowParam, wholeRowL, Lreq, whatFromN, QDefaultValsBefore, DefaultValParam);
                for(int i=1; i<=Q; i++){
                    Lcur=size.GetLength(i);
                    n=Lcur;
                    Array1DAddElement(X[i-1], n, row[i-1]);
                }
            }else{
                for(int i=1; i<=Q; i++){
                    Lcur=size.GetLength(i);
                    n=Lcur;
                    Array1DAddElement(X[i-1], n, DefaultVal);
                }
            }
            if(Lmin!=Lmax || (Lmin==Lmax && size.GetIneRowsLength()==0)){
                for(int i=1; i<=Q; i++){
                    Lcur=size.GetLength(i);
                    size.SetLength(Lcur+1, i);
                }
            }else{
                size.SetIneRowsLength(Lmax+1);
            }
        }
    }
}//add ine row
template<typename T>ArrayOf2DBy1DAddIneRow(T**&X, Array2DSize& size, std::vector<T>rowParam, int IfNonRectIgnore0Add1Stretch2=0, int whatFromN=1, int QDefaultValsBefore=0, T*DefaultValParam=NULL){
    ArrayOf2DBy1DAddIneRow(X, size, rowParam.data(), rowParam.size(),  IfNonRectIgnore0Add1Stretch2, whatFromN, QDefaultValsBefore, DefaultValParam);
}
template<typename T>ArrayOf2DBy1DAddIneRow(std::vector<std::vector<T>>&X, std::vector<T>rowParam, int IfNonRectIgnore0Add1Stretch2=0, int whatFromN=1, int QDefaultValsBefore=0, T*DefaultValParam=NULL){
    int Q=X.size(), wholeRowL=rowParam.size();
    std::vector<T>rowM, rowCur;
    Array2DSize size;
    size=GetSizeOf2DVector(X);
    T curVal, DefaultVal;
    if(DefaultValParam!=NULL){
        DefaultVal=(*(DefaultValParam));
    }
    int Lmin=size.GetMinLength(), Lmax=size.GetMaxLength(), Lreq=0, Lcur, n;
    if(Q==0 && wholeRowL>0){
        Lreq=0;
        rowM=Array1DGetSubArray_byLs_PreMark(rowParam, Lreq, whatFromN, QDefaultValsBefore, DefaultValParam);
        Lcur=rowM.size();//in other if brunch Lcur means different thing! n not used
        //ob not from first et co Defaults before
        for(int i=1; i<=Lcur; i++){
            rowCur.clear();
            curVal=rowM[i-1];
            rowCur.push_back(curVal);
            X.push_back(rowCur);
        }
    }else{
        if(!((Lmin!=Lmax)&&(IfNonRectIgnore0Add1Stretch2==0))){
            Lreq=Q; //ob os ine row
            rowM=Array1DGetSubArray_byLs_PreMark(rowParam, Lreq, whatFromN, QDefaultValsBefore, DefaultValParam);
            if(Lmin!=Lmax && IfNonRectIgnore0Add1Stretch2==2){
                for(int i=1; i<=Q; i++){
                    Lcur=X[i-1].size();
                    //n=Lcur;
                    for(int j=Lcur+1; j<=Lmax; j++){
                        Array1DAddElement(X[i-1], DefaultVal);
                    }
                }
                for(int i=1; i<=Q; i++){
                    size.SetLength(Lmax, i);
                }
            }
            if(wholeRowL>0){
                for(int i=1; i<=Q; i++){
                    Lcur=size.GetLength(i);
                    n=Lcur;
                    Array1DAddElement(X[i-1], rowM[i-1]);
                }
            }else{
                for(int i=1; i<=Q; i++){
                    Lcur=size.GetLength(i);
                    n=Lcur;
                    Array1DAddElement(X[i-1], DefaultVal);
                }
            }
            if(Lmin==Lmax && size.GetIneRowsLength()>0){
                size.SetIneRowsLength(Lmax+1);
            }else{
                for(int i=1; i<=Q; i++){
                     Lcur=size.GetLength(i);
                     size.SetLength(Lcur+1, i);
                }
            }
        }
    }
}//add ine row
template<typename T>ArrayOf2DBy1DAddIneRow(std::vector<std::vector<T>>&X, T*rowParam=NULL, int wholeRowL=0, int IfNonRectIgnore0Add1Stretch2=0, int whatFromN=1, int QDefaultValsBefore=0, T*DefaultValParam=NULL){
    int Q=X.size(), Lreq=0;
    std::vector<T>row;
    if(wholeRowL==0 && rowParam!=NULL){
        wholeRowL=Q;
    }
    Lreq=Q;//and 0 if 0
    row=Array1DGetSubArray_byLs_PreMark(rowParam, wholeRowL, Lreq, whatFromN, QDefaultValsBefore, DefaultValParam);
    //template<typename T>ArrayOf2DBy1DAddIneRow(std::vector<std::vector<T>>&X, std::vector<T>rowParam, int IfNonRectIgnore0Add1Stretch2=0, int whatFromN=1, int QDefaultValsBefore=0, T*DefaultValParam=NULL){
    ArrayOf2DBy1DAddIneRow(X, row, IfNonRectIgnore0Add1Stretch2, whatFromN, QDefaultValsBefore, DefaultValParam);
}////add ine row
//Ins ine row
template<typename T>ArrayOf2DBy1DInsIneRow(T**&X, Array2DSize& size, int IneRowPosN, T*rowParam=NULL, int wholeRowL=0, int IfAfterLminIgnore0Stretch2=0, int whatFromN=1, int QDefaultValsBefore=0, T*DefaultValParam=NULL){
    int Q=size.GetQExtRows(), Lmin=size.GetMinLength(), Lmax=size.GetMaxLength(), Lcur, n, Lreq;
    std::vector<T>row;
    T DefaultVal;
    if(DefaultValParam!=NULL){
        DefaultVal=(*(DefaultValParam));
    }
    if(rowParam!=NULL && wholeRowL==0){
        wholeRowL=Q;
    }
    if(Q>0 && IneRowPosN>=1 && IneRowPosN<=Lmax && !(IneRowPosN>Lmin && IfAfterLminIgnore0Stretch2==0)){
        if(IneRowPosN>Lmin && IfAfterLminIgnore0Stretch2==2){
            for(int i=1; i<=Q; i++){
                Lcur=size.GetLength(i);
                n=Lcur;
                for(int j=Lcur+1; j<=Lmax; j++){
                    Array1DAddElement(X[i-1], n, DefaultVal);
                }
            }
        }
        Lreq=Q;//and 0 if 0
        row=Array1DGetSubArray_byLs_PreMark(rowParam, wholeRowL, Lreq, whatFromN, QDefaultValsBefore, DefaultValParam);
        if(rowParam!=NULL){
            for(int i=1; i<=Q; i++){
                ArrayOf2DBy1DInsToExtRowN(X, size, i, row[i-1], IneRowPosN);
            }
        }else{
            for(int i=1; i<=Q; i++){
                ArrayOf2DBy1DInsToExtRowN(X, size, i, DefaultVal, IneRowPosN);
            }
        }
    }
}
template<typename T>ArrayOf2DBy1DInsIneRow(T**&X, Array2DSize& size, int IneRowPosN, std::vector<T>rowParam, int IfAfterLminIgnore0Stretch2=0, int whatFromN=1, int QDefaultValsBefore=0, T*DefaultValParam=NULL){
    ArrayOf2DBy1DInsIneRow(X, size, IneRowPosN, rowParam.data(), rowParam.size(), IfAfterLminIgnore0Stretch2, whatFromN, QDefaultValsBefore, DefaultValParam);
}
template<typename T>ArrayOf2DBy1DInsIneRow(std::vector<std::vector<T>>&X, int IneRowPosN, std::vector<T>rowParam, int IfAfterLminIgnore0Stretch2=0, int whatFromN=1, int QDefaultValsBefore=0, T*DefaultValParam=NULL){
    std::vector<T>row;
    T DefaultVal;
    Array2DSize size;
    size=GetSizeOf2DVector(X);
    int Q=size.GetQExtRows(), Lmin=size.GetMinLength(), Lmax=size.GetMaxLength(), Lcur, n, Lreq;
    if(DefaultValParam!=NULL){
        DefaultVal=(*(DefaultValParam));
    }
    if(Q>0 && IneRowPosN>=1 && IneRowPosN<=Lmax && !(IneRowPosN>Lmin && IfAfterLminIgnore0Stretch2==0)){
        if(IneRowPosN>Lmin && IfAfterLminIgnore0Stretch2==2){
            for(int i=1; i<=Q; i++){
                Lcur=size.GetLength(i);
                //n=Lcur;
                for(int j=Lcur+1; j<=Lmax; j++){
                    Array1DAddElement(X[i-1],  DefaultVal);
                }
            }
        }
        Lreq=Q;//and 0 if 0
        row=Array1DGetSubArray_byLs_PreMark(rowParam, Lreq, whatFromN, QDefaultValsBefore, DefaultValParam);
        if(rowParam.size()!=0){
            for(int i=1; i<=Q; i++){
                ArrayOf2DBy1DInsToExtRowN(X, i, row[i-1], IneRowPosN);
            }
        }else{
            for(int i=1; i<=Q; i++){
                ArrayOf2DBy1DInsToExtRowN(X, i, DefaultVal, IneRowPosN);
            }
        }
    }
}
template<typename T>ArrayOf2DBy1DInsIneRow(std::vector<std::vector<T>>&X, int IneRowPosN, T*rowParam=NULL, int wholeRowL=0,  int IfAfterLminIgnore0Stretch2=0, int whatFromN=1, int QDefaultValsBefore=0, T*DefaultValParam=NULL){
    std::vector<T>row;
    T DefaultVal;
    Array2DSize size;
    size=GetSizeOf2DVector(X);
    int Q=size.GetQExtRows(), Lmin=size.GetMinLength(), Lmax=size.GetMaxLength(), Lcur, n, Lreq;
    //if(wholeRowL==0)wholeRowL=Lmin;
    Lreq=Q;//and 0 if 0
    //template <typename T>std::vector<T> Array1DGetSubArray_byLs_PreMark(std::vector<T>x, int Lreq=0, int whatFromN=1, int FirstDefaultValues=0, T*DfltValParam=NULL){
    //template <typename T>std::vector<T> Array1DGetSubArray_byLs_PreMark(T*x, int Lold, int Lreq=0, int whatFromN=1, int FirstDefaultValues=0, T*DfltValParam=NULL){
    row=Array1DGetSubArray_byLs_PreMark(rowParam, wholeRowL, Lreq, whatFromN, QDefaultValsBefore, DefaultValParam);
    //template<typename T>ArrayOf2DBy1DInsIneRow(std::vector<std::vector<T>>&X, int IneRowPosN, std::vector<T>rowParam, int IfAfterLminIgnore0Stretch2=0, int whatFromN=1, int QDefaultValsBefore=0, T*DefaultValParam=NULL){
    ArrayOf2DBy1DInsIneRow(X, IneRowPosN, row, IfAfterLminIgnore0Stretch2, whatFromN, QDefaultValsBefore, DefaultValParam);
}
//del ine row
template<typename T>ArrayOf2DBy1DDelIneRow(T**&X, Array2DSize& size, int IneRowPosN, bool ifPos0IgnoreNotDelLastOfVariaNotIgnore=false){
    int Q=size.GetQExtRows(), Lmin=size.GetMinLength(), Lmax=size.GetMaxLength(), Lcur, n;
    if(IneRowPosN>=1 && IneRowPosN<=Lmax){
        for(int i=1; i<=Q; i++){
            ArrayOf2DBy1DDelFromExtRowN(X, size, i, IneRowPosN);
        }
    }else if(IneRowPosN==0 && ifPos0IgnoreNotDelLastOfVariaNotIgnore==true && Lmin!=Lmax){
        for(int i=1; i<=Q; i++){
            Lcur=size.GetLength(i);
            n=Lcur;
            Array1DDelElementFromN(X[i-1], n, Lcur);
            size.SetLength(n, i);
        }
    }
}
template<typename T>ArrayOf2DBy1DDelIneRow(std::vector<std::vector<T>>&X, int IneRowPosN, bool ifPos0IgnoreNotDelLastOfVariaNotIgnore=false){
    Array2DSize size;
    size=GetSizeOf2DVector(X);
    int Q=size.GetQExtRows(), Lmin=size.GetMinLength(), Lmax=size.GetMaxLength(), Lcur;
    if(IneRowPosN>=1 && IneRowPosN<=Lmax){
        for(int i=1; i<=Q; i++){
            ArrayOf2DBy1DDelFromExtRowN(X, i, IneRowPosN);
        }
    }else if(IneRowPosN==0 && ifPos0IgnoreNotDelLastOfVariaNotIgnore==true && Lmin!=Lmax){
        for(int i=1; i<=Q; i++){
            Lcur=size.GetLength(i);
            Array1DDelElementFromN(X[i-1], Lcur);
        }
    }
}
//Assign
template<typename T>void ArrayOf2DBy1DAssign(T**&XTo, Array2DSize& sizeTo, T**XFrom, const Array2DSize& sizeFrom){
    //void ArrayOf2DBy1DSetSize(T**&X, const Array2DSize&size1, const Array2DSize&size2, T*DfltValParam=NULL, bool PreserveVals=true)
    T*DfltValParam=NULL;
    int Q=sizeFrom.GetQExtRows(), Lcur;
    ArrayOf2DBy1DSetSize(XTo, sizeTo, sizeFrom, DfltValParam, false);
    for(int i=1; i<=Q; i++){
       Lcur=sizeFrom.GetLength(i);
       for(int j=1; j<=Lcur; j++){
          XTo[i-1][j-1]=XFrom[i-1][j-1];
       }
    }
}
template<typename T>void ArrayOf2DBy1DAssign(T**&X, Array2DSize& size, std::vector<std::vector<T>>Z){
    int Q=size.GetQExtRows(), Lcur;
    //Array2DSize  size1=GetSizeOf2DVector(X);
    T val;
    if(X!=NULL){
        for(int i=1; i<=Q; i++){
            delete[]X[i-1];
        }
        delete[]X;
        size=GetSizeOf2DVector(Z);
        Q=size.GetQExtRows();
    }
    if(X.size()!=0){
        X=new T*[Q];
        for(int i=1; i<=Q; i++){
            Lcur=size.GetLength(i);
            X[i-1]=new T[Lcur];
        }
         for(int i=1; i<=Q; i++){
             Lcur=size.GetLength(i);
             for(int j=1; j<=Lcur; j++){
                 val=Z[i-1][j-1];
                 X[i-1][j-1]=val;
             }
         }
    }
}
template<typename T>void ArrayOf2DBy1DAssign(std::vector<std::vector<T>>&X, Array2DSize& size, T**Z=NULL, T*DefaultValParam=NULL){
    int Q=size.GetQExtRows(), Lcur;
    //Array2DSize  size1=GetSizeOf2DVector(X);
    T curVal, DefaultVal;
    X.clear();
    if(DefaultValParam!=NULL){
        DefaultVal=(*(DefaultValParam));
    }
    std::vector<T>row;
    if(Z!=NULL){
        for(int i=1; i<=Q; i++){
            Lcur=size.GetLength(i);
            row.clear();
            for(int j=1; j<=Lcur; j++){
                curVal=Z[i-1][j-1];
                row.push_back(curVal);
            }
            X.push_back(row);
        }
    }else{
        for(int i=1; i<=Q; i++){
            Lcur=size.GetLength(i);
            row.clear();
            for(int j=1; j<=Lcur; j++){
                curVal=DefaultVal;
                row.push_back(curVal);
            }
            X.push_back(row);
        }
    }
}
//
//get val
template<typename T> T ArrayOf2DBy1DGetVal(T**X, const Array2DSize&size, int extRowN, int ineRowN){//40
    T y;
    if(extRowN>=1 && extRowN<=size.GetQExtRows() && ineRowN>=1 && ineRowN<=size.GetLength(extRowN)){
       y=X[extRowN][ineRowN];
    }
    return y;
}
template<typename T> T* ArrayOf2DBy1DGetVal_AsPtr(T**X, const Array2DSize&size, int extRowN, int ineRowN){//40
    T y;
    if(extRowN>=1 && extRowN<=size.GetQExtRows() && ineRowN>=1 && ineRowN<=size.GetLength(extRowN)){
       y=X[extRowN][ineRowN];
    }
    return y;
}
template<typename T> T ArrayOf2DBy1DGetVal(std::vector<std::vector<T>>X, int extRowN, int ineRowN){//40
    T y;
    if(extRowN>=1 && extRowN<=X.size() && ineRowN>=1 && ineRowN<=X[extRowN-1].size){
       y=X[extRowN][ineRowN];
    }
    return y;
}
template<typename T> T* ArrayOf2DBy1DGetVal_AsPtr(std::vector<std::vector<T>>X, int extRowN, int ineRowN){//40
    T*y=NULL;
    if(extRowN>=1 && extRowN<=X.size() && ineRowN>=1 && ineRowN<=X[extRowN-1].size()){
       y=X[extRowN][ineRowN];
    }
    return y;
}
//
template<typename T> std::vector<T> ArrayOf2DBy1DGetExtRowN_AsVals(T**X, int ExtRowN){//30
    std::vector<T>y;
    int L;
    if(ExtRowN>=1 && ExtRowN<=X.size()){
        L=X[ExtRowN-1].size();
        for(int i=1; i<=L; i++){
            y.push_back(X[ExtRowN-1][i-1]);
        }
    }
    return y;
}
template<typename T> T* ArrayOf2DBy1DGetExtRowN_AsPtr(T**X, const Array2DSize&size, int ExtRowN){//30
    T*y=NULL;
    if(size.ExtRowNBelongsTo(ExtRowN)){
        y=X[ExtRowN-1];
    }
    return y;
}
template<typename T> std::vector<T> ArrayOf2DBy1DGetExtRowN(std::vector<std::vector<T>>X, int ExtRowN){//30
    std::vector<T>y;
    if(ExtRowN>=1 && ExtRowN<=X.size()){
        y=X[ExtRowN-1];
    }
    return y;
}
//transpose
//template <typename T>bool  ArrayOf2DBy1DSizeIsForTtanspose(std::vector<std::vector<T>>X){
//    Array2DSize size=GetSizeOf2DVector(X);
//    return size.IsForTranspose();
//}
template<typename T>void ArrayOf2DBy1DTranspose(T**&X, Array2DSize&size){//
    T**Y=NULL;
    int Q=size.GetQExtRows(), Lmax=size.GetMaxLength(), Lmin=size.GetMinLength(), Lcur;
    //if(X!=NULL && Q>0 && L>0 && size.GetMinLength()==L){
    bool isForTranspose=size.IsForTranspose();
    if(Q>0 && Lmax>0 && Lmin>0 && isForTranspose){
    //LengthFullOfIneRowN(int N){//from start to 1st stop
        Y=new T*[Lmax];
        for(int i=1; i<=Lmax; i++){
            Lcur=size.LengthFullOfIneRowN(i);
            Y[i-1]=new T[Lcur];
        }
        //
        for(int i=1; i<=Lmax; i++){
            Lcur=size.LengthFullOfIneRowN(i);
            for(int j=1; j<=Lcur; j++){
                Y[i-1][j-1]=X[j-1][i-1];
            }
        }
        //
        for(int i=1; i<=Q; i++){
            delete[]X[i-1];
        }
        delete[]X;
        //
        X=Y;
        //
        size.Transpose();
    }
}
template<typename T>void ArrayOf2DBy1DTranspose(std::vector<std::vector<T>>&X){
    std::vector<std::vector<T>>Y;
    std::vector<T>Z;
    T val;
    Array2DSize size=GetSizeOf2DVector(X);
    int Q=size.GetQExtRows(), Lmax=size.GetMaxLength(), Lmin=size.GetMinLength(), hcur;
    bool isForTranspose=size.IsForTranspose();
    if(Q>0 && Lmax>0 && Lmin>0 && isForTranspose){
        for(int i=1; i<=Lmax; i++){
            Z.clear();
            hcur=size.LengthFullOfIneRowN(i);
            for(int j=1; j<=hcur; j++){
                val=X[j-1][i-1];
                Z.push_back(val);
            }
            Y.push_back(Z);
        }
        X=Y;
    }
    //X=Y;//may gbe this is more correct then place this operator before as is realized
}
template<typename T>void ArrayOf2DBy1DTransposeTo(std::vector<std::vector<T>>X){
    std::vector<std::vector<T>>Y;
    Y=X;
    ArrayOf2DBy1DTranspose(Y);
    return Y;
}

/*
template<typename T>std::vector<std::vector<T>> ArrayOf2DBy1DTransposeTo(T**X, const Array2DSize&size){//
    std::vector<std::vector<T>>Y;
    std::vector<T>Z;
    ///if(size.GetMaxLength()==L){
    //    for(int i=1; i<=L; i++){
    //        Z.clear();
    //        for(int j=1; j<=Q; j++){
    //            Z.push_back(X[j-1][i-1]);
    //        }
    //        Y.push_back(Z);
    //    }
    //}
    int Q=size.GetQExtRows(), Lmax=size.GetMaxLength(), Lmin=size.GetMinLength(), Lcur;
    bool isForTranspose=size.IsForTranspose();
    if(Q>0 && Lmax>0 && Lmin>0 && isForTranspose){
        for(int i=1; i<=Lmax; i++){
            Z.clear();
            Lcur=size.LengthFullOfIneRowN(i);
            for(int j=1; j<=Lcur; j++){
               Z.push_back(X[j-1][i-1]);
            }
            Y.push_back(Z);
        }
    }
    return Y;
}
template<typename T>std::vector<std::vector<T>> ArrayOf2DBy1DTransposeTo(std::vector<std::vector<T>>X){//rect only
    std::vector<std::vector<T>>Y;
    std::vector<T>Z;
    int Q=X.size(), L=X[1-1].size();
    for(int i=1; i<=L; i++){
        Z.clear();
        for(int j=1; j<=Q; j++){
           Z.push_back(X[j-1][i-1]);
       }
       Y.push_back(Z);
    }
    return Y;
}
//transose struct // not tested
template<typename T>std::vector<std::vector<T>> ArrayOf2DBy1DTransposeStructTo(T**&X, Array2DSize&size){
    std::vector<std::vector<T>>Y;
    std::vector<T> rowCur;
    int Q=size.GetQExtRows, Lmin=size.GetMinLength(), Lmax=size.GetMaxLength(), Lcur, LIneCur, L1=size.GetLength(1);
    if (size.IsForTranspose()){
        for(int i=1; i<=L1; i++){
            LIneCur=size.LengthFullOfIneRowN(i);
            rowCur.clear();
            rowCur=ArrayOf2DBy1DGetIneRowN(X, i);
            Y.push_back(rowCur);
        }
    }
    return Y;
}
template<typename T>void ArrayOf2DBy1DTransposeStruct(T**&X, Array2DSize&size){
    T**Y=NULL;
    std::vector<T> rowCur;
    int Q=size.GetQExtRows, Lmin=size.GetMinLength(), Lmax=size.GetMaxLength(), Lcur, LIneCur, L1=size.GetLength(1);
    if (size.IsForTranspose()){
        //create
        Y=new T*[L1];
        for(int i=1; i<=L1; i++){
            LIneCur=size.LengthFullOfIneRowN(i);
            Y[i-1]=new T[LIneCur];
        }
        //copy
        for(int i=1; i<=L1; i++){
            LIneCur=size.LengthFullOfIneRowN(i);
            rowCur.clear();
            rowCur=ArrayOf2DBy1DGetIneRowN(X, i);
            for(int j=1; j<=LIneCur; j++){
                Y[i-1][j-1]=X[j-1][i-1];
            }
        }
        //del
        for(int i=1; i<=Q; i++){
            delete[]X[i-1];
        }
        delete[]X;
        //assign
        X=Y;
    }
    return Y;
}
template<typename T>std::vector<std::vector<T>> ArrayOf2DBy1DTransposeStructTo(std::vector<std::vector<T>>X){
    Array2DSize size=GetSizeOf2DVector(X);
    std::vector<std::vector<T>>Y;
    std::vector<T> rowCur;
    int Q=size.GetQExtRows, Lmin=size.GetMinLength(), Lmax=size.GetMaxLength(), Lcur, LIneCur, L1=size.GetLength(1);
    if (size.IsForTranspose()){
        for(int i=1; i<=L1; i++){
            LIneCur=size.LengthFullOfIneRowN(i);
            rowCur.clear();
            rowCur=ArrayOf2DBy1DGetIneRowN(X, i);
            Y.push_back(rowCur);
        }
    }
    return Y;
}
template<typename T>void ArrayOf2DBy1DTransposeStructTo(std::vector<std::vector<T>>&X){
    Array2DSize size=GetSizeOf2DVector(X);
    std::vector<std::vector<T>>Y;
    std::vector<T> rowCur;
    int Q=size.GetQExtRows, Lmin=size.GetMinLength(), Lmax=size.GetMaxLength(), Lcur, LIneCur, L1=size.GetLength(1);
    if (size.IsForTranspose()){
        for(int i=1; i<=L1; i++){
            LIneCur=size.LengthFullOfIneRowN(i);
            rowCur.clear();
            rowCur=ArrayOf2DBy1DGetIneRowN(X, i);
            Y.push_back(rowCur);
        }
    }
    X=Y;
}
*/
//--------------------------------------------------
//Seek
template <typename T> std::vector<std::vector<int>> ArrayOf2DBy1D_SeekVal_Simply(T**X, const Array2DSize& size, const T& Val, int paramN=0){
    //std::vector<T>rowVect=Array1DAssign(row, L);
    std::vector<std::vector<int>>NsAll;
    std::vector<int>NsSnglLine;
    std::vector<int>pair;
    int Q=size.GetQExtRows(), count, N;
    int FromN=1, ToN=0;
    for(int i=1; i<=Q; i++){
        pair.clear();
        count=0;
        NsSnglLine=Array1DSeekValSimply(X[i-1], Val, size.GetLength(i), 1, 0, paramN);
        count=NsSnglLine.size();
        if(count>0){
            for(int j=1; j<=count; j++){
                N=NsSnglLine[j-1];
                pair.push_back(i);
                pair.push_back(N);
                NsAll.push_back(pair);
            }
        }
    }
    return NsAll;
}

template <typename T> std::vector<std::vector<int>> ArrayOf2DBy1D_SeekExtRow_Simply(T**X, const Array2DSize& size, T*row, int L, int paramN=0){
    std::vector<T>rowVect=Array1DAssign(row, L);
    std::vector<std::vector<int>>NsAll;
    std::vector<int>NsSnglLine;
    std::vector<int>pair;
    int Q=size.GetQExtRows(), count, N;
    int FromN=1, ToN=0;
    for(int i=1; i<=Q; i++){
        pair.clear();
        count=0;
        NsSnglLine=Array1DSeekSubArraySimply(X[i-1], row, size.GetLength(i), size.GetLength(i-1), L, FromN, ToN, paramN);
        count=NsSnglLine.size();
        if(count>0){
            for(int j=1; j<=count; j++){
                N=NsSnglLine[j-1];
                pair.push_back(i);
                pair.push_back(N);
                NsAll.push_back(pair);
            }
        }
    }
    return NsAll;
}
template <typename T> std::vector<std::vector<int>> ArrayOf2DBy1D_SeekExtRow_Simply(T**X, const Array2DSize& size, std::vector<T>row, int paramN=0){
    ArrayOf2DBy1D_SeekExtRow(row.data(), row.size(), paramN);
}


template <typename T>bool ArrayOf2DBy1D_Arrs2ndIsAtPosIn1st(T**Where, const Array2DSize& sizeWhere, T**What, const Array2DSize& sizeWhat, int ExtRowN, int IneRowN, bool IsTransposed=false, bool EachExtRowIsReversed=false){
    int whereQ=sizeWhere.GetQExtRows(), whatQ=sizeWhat.GetQExtRows(), whereL, whatL, whereH, whatH, whereLcur, whatLcur, whereExtRowN, whereIneRowN, whatExtRowN, whatIneRowN;
    bool verdict=true;
    if(IsTransposed==false){
        if(ExtRowN+whatQ>=whereQ){
            verdict=false;
        }else{
            if(EachExtRowIsReversed==false){
                // 11 12 13
                // 21 22
                // 31 32 33
                for(whatExtRowN=1; whatExtRowN<=whatQ; whatExtRowN++){
                    whereExtRowN=ExtRowN+whatExtRowN-1;
                    whereL=sizeWhere.GetLength(whereExtRowN);
                    whatL=sizeWhat.GetLength(whatExtRowN);
                    if(IneRowN-1+whatL>whereL){
                        verdict=false;
                    }
                }
                if(!verdict==false){
                    for(whatExtRowN=1; whatExtRowN<=whatQ; whatExtRowN++){
                        whereExtRowN=ExtRowN+whatExtRowN-1;
                        whereL=sizeWhere.GetLength(whereExtRowN);
                        whatL=sizeWhat.GetLength(whatExtRowN);
                        for(whatIneRowN=1; whatIneRowN<=whatL; whatIneRowN++){
                            whereIneRowN=IneRowN+whatIneRowN-1;
                            if(Where[whereExtRowN-1][whereIneRowN-1]!=What[whatExtRowN-1][whatIneRowN-1]){
                                verdict=false;
                            }
                        }
                    }
                }
            }else{//EachExtRowIsReversed==true
                // 13 12 11
                //    22 21
                // 33 32 31
                for(whatExtRowN=1; whatExtRowN<=whatQ; whatExtRowN++){
                    whereExtRowN=ExtRowN+whatExtRowN-1;
                    whereL=sizeWhere.GetLength(whereExtRowN);
                    whatL=sizeWhat.GetLength(whatExtRowN);
                    if(IneRowN-whatL+1<1){
                        verdict=false;
                    }
                }
                if(!verdict==false){
                    for(whatExtRowN=1; whatExtRowN<=whatQ; whatExtRowN++){
                        whereExtRowN=ExtRowN+whatExtRowN-1;
                        whereL=sizeWhere.GetLength(whereExtRowN);
                        whatL=sizeWhat.GetLength(whatExtRowN);
                        for(whatIneRowN=1; whatIneRowN<=whatL; whatIneRowN++){
                            whereIneRowN=IneRowN-whatIneRowN+1;
                            if(Where[whereExtRowN-1][whereIneRowN-1]!=What[whatExtRowN-1][whatIneRowN-1]){
                                verdict=false;
                            }
                        }
                    }
                }
            }
        }
    }else{//transposed //task s'harder
        if(EachExtRowIsReversed==false){
            // 11 21 31
            // 12 22 32
            // 13    33
            for(whatExtRowN=1; whatExtRowN<=whatQ; whatExtRowN++){
                whereIneRowN=IneRowN+whatExtRowN-1;
                whatL=sizeWhat.GetLength(whatExtRowN);
                for(whatIneRowN=1; whatIneRowN<=whatL; whatIneRowN++){
                    whereExtRowN=ExtRowN+whatIneRowN-1;
                    whereL=sizeWhere.GetLength(whereExtRowN);
                    if(whereExtRowN>=1 && whereExtRowN<=whereQ && whereIneRowN>=1 &&whereIneRowN<=whereL){
                        //NOp
                    }else{
                        verdict=false;
                    }
                }
            }
            if(!verdict==false){
                for(whatExtRowN=1; whatExtRowN<=whatQ; whatExtRowN++){
                    whereIneRowN=IneRowN+whatExtRowN-1;
                    whatL=sizeWhat.GetLength(whatExtRowN);
                    for(whatIneRowN=1; whatIneRowN<=whatL; whatIneRowN++){
                        whereExtRowN=ExtRowN+whatIneRowN-1;
                        whereL=sizeWhere.GetLength(whereExtRowN);
                        if(Where[whereExtRowN-1][whereIneRowN-1]!=What[whatExtRowN-1][whatIneRowN-1]){
                            verdict=false;
                        }
                    }
                }
            }
        }else{
            // 13    33
            // 12 22 32
            // 11 21 31
            for(whatExtRowN=1; whatExtRowN<=whatQ; whatExtRowN++){
                whereIneRowN=IneRowN+whatExtRowN-1;
                whatL=sizeWhat.GetLength(whatExtRowN);
                for(whatIneRowN=1; whatIneRowN<=whatL; whatIneRowN++){
                    whereExtRowN=ExtRowN-whatIneRowN+1;
                    whereL=sizeWhere.GetLength(whereExtRowN);
                    if(whereExtRowN>=1 && whereExtRowN<=whereQ && whereIneRowN>=1 &&whereIneRowN<=whereL){
                        //NOp
                    }else{
                        verdict=false;
                    }
                }
            }
            if(!verdict==false){
                for(whatExtRowN=1; whatExtRowN<=whatQ; whatExtRowN++){
                    whereIneRowN=IneRowN+whatExtRowN-1;
                    whatL=sizeWhat.GetLength(whatExtRowN);
                    for(whatIneRowN=1; whatIneRowN<=whatL; whatIneRowN++){
                        whereExtRowN=ExtRowN-whatIneRowN+1;
                        whereL=sizeWhere.GetLength(whereExtRowN);
                        if(Where[whereExtRowN-1][whereIneRowN-1]!=What[whatExtRowN-1][whatIneRowN-1]){
                            verdict=false;
                        }
                    }
                }
            }
        }
    }
    return verdict;
}
template <typename T>bool ArrayOf2DBy1D_Arrs2ndIsAtPosIn1st(std::vector<std::vector<T>>Where,  std::vector<std::vector<T>>What, int ExtRowN, int IneRowN, bool IsTransposed=false, bool EachExtRowIsReversed=false){
    Array2DSize sizeWhere=GetSizeOf2DVector(Where), sizeWhat=GetSizeOf2DVector(What);
    int whereQ=sizeWhere.GetQExtRows(), whatQ=sizeWhat.GetQExtRows(), whereL, whatL, whereH, whatH, whereLcur, whatLcur, whereExtRowN, whereIneRowN, whatExtRowN, whatIneRowN;
    bool verdict=true;
    if(IsTransposed==false){
        if(ExtRowN+whatQ>=whereQ){
            verdict=false;
        }else{
            if(EachExtRowIsReversed==false){
                // 11 12 13
                // 21 22
                // 31 32 33
                for(whatExtRowN=1; whatExtRowN<=whatQ; whatExtRowN++){
                    whereExtRowN=ExtRowN+whatExtRowN-1;
                    whereL=sizeWhere.GetLength(whereExtRowN);
                    whatL=sizeWhat.GetLength(whatExtRowN);
                    if(IneRowN-1+whatL>whereL){
                        verdict=false;
                    }
                }
                if(!verdict==false){
                    for(whatExtRowN=1; whatExtRowN<=whatQ; whatExtRowN++){
                        whereExtRowN=ExtRowN+whatExtRowN-1;
                        whereL=sizeWhere.GetLength(whereExtRowN);
                        whatL=sizeWhat.GetLength(whatExtRowN);
                        for(whatIneRowN=1; whatIneRowN<=whatL; whatIneRowN++){
                            whereIneRowN=IneRowN+whatIneRowN-1;
                            if(Where[whereExtRowN-1][whereIneRowN-1]!=What[whatExtRowN-1][whatIneRowN-1]){
                                verdict=false;
                            }
                        }
                    }
                }
            }else{//EachExtRowIsReversed==true
                // 13 12 11
                //    22 21
                // 33 32 31
                for(whatExtRowN=1; whatExtRowN<=whatQ; whatExtRowN++){
                    whereExtRowN=ExtRowN+whatExtRowN-1;
                    whereL=sizeWhere.GetLength(whereExtRowN);
                    whatL=sizeWhat.GetLength(whatExtRowN);
                    if(IneRowN-whatL+1<1){
                        verdict=false;
                    }
                }
                if(!verdict==false){
                    for(whatExtRowN=1; whatExtRowN<=whatQ; whatExtRowN++){
                        whereExtRowN=ExtRowN+whatExtRowN-1;
                        whereL=sizeWhere.GetLength(whereExtRowN);
                        whatL=sizeWhat.GetLength(whatExtRowN);
                        for(whatIneRowN=1; whatIneRowN<=whatL; whatIneRowN++){
                            whereIneRowN=IneRowN-whatIneRowN+1;
                            if(Where[whereExtRowN-1][whereIneRowN-1]!=What[whatExtRowN-1][whatIneRowN-1]){
                                verdict=false;
                            }
                        }
                    }
                }
            }
        }
    }else{//transposed //task s'harder
        if(EachExtRowIsReversed==false){
            // 11 21 31
            // 12 22 32
            // 13    33
            for(whatExtRowN=1; whatExtRowN<=whatQ; whatExtRowN++){
                whereIneRowN=IneRowN+whatExtRowN-1;
                for(whatIneRowN=1; whatIneRowN<=whatQ; whatIneRowN++){
                    whereExtRowN=ExtRowN+whatIneRowN-1;
                    whereL=sizeWhere.GetLength(whereExtRowN);
                    if(whereExtRowN>=1 && whereExtRowN<=whereQ && whereIneRowN>=1 &&whereIneRowN<=whereL){
                        verdict=false;
                    }
                }
            }
            if(!verdict==false){
                for(whatExtRowN=1; whatExtRowN<=whatQ; whatExtRowN++){
                    whereIneRowN=IneRowN+whatExtRowN-1;
                    for(whatIneRowN=1; whatIneRowN<=whatQ; whatIneRowN++){
                        whereExtRowN=ExtRowN+whatIneRowN-1;
                        whereL=sizeWhere.GetLength(whereExtRowN);
                        if(Where[whereExtRowN-1][whereIneRowN-1]!=What[whatExtRowN-1][whatIneRowN-1]){
                            verdict=false;
                        }
                    }
                }
            }
        }else{
            // 13    33
            // 12 22 32
            // 11 21 31
            for(whatExtRowN=1; whatExtRowN<=whatQ; whatExtRowN++){
                whereIneRowN=IneRowN+whatExtRowN-1;
                for(whatIneRowN=1; whatIneRowN<=whatQ; whatIneRowN++){
                    whereExtRowN=ExtRowN-whatIneRowN+1;
                    whereL=sizeWhere.GetLength(whereExtRowN);
                    if(whereExtRowN>=1 && whereExtRowN<=whereQ && whereIneRowN>=1 &&whereIneRowN<=whereL){
                        verdict=false;
                    }
                }
            }
            if(!verdict==false){
                for(whatExtRowN=1; whatExtRowN<=whatQ; whatExtRowN++){
                    whereIneRowN=IneRowN+whatExtRowN-1;
                    for(whatIneRowN=1; whatIneRowN<=whatQ; whatIneRowN++){
                        whereExtRowN=ExtRowN-whatIneRowN+1;
                        whereL=sizeWhere.GetLength(whereExtRowN);
                        if(Where[whereExtRowN-1][whereIneRowN-1]!=What[whatExtRowN-1][whatIneRowN-1]){
                            verdict=false;
                        }
                    }
                }
            }
        }
    }
    return verdict;
}

/*template <typename T> bool ArrayOf2DBy1D_SubArr2DIsAtPos_SameDataSucc(T**Where, const Array2DSize& sizeWhere, T**What, const Array2DSize& sizeWhat, int ExtRowN, int IneRowN){
    bool b=true;
    int WhereExtRowN, WhereIneRowN, WhereCurL, WhatExtRowN, WhatIneRowN, WhatCurL, WhereQ=sizeWhere.GetQExtRows(), WhatQ=sizeWhat.GetQExtRows();
    T WhereVal, WhatVal;
    for(WhatExtRowN=1; WhatExtRowN<=WhatQ; WhatExtRowN++){
        WhereExtRowN=WhatExtRowN+ExtRowN-1;
        WhatCurL=sizeWhat.GetLength(WhatExtRowN);
        WhereCurL=sizeWhere.GetLength(WhereCurL);
        for(WhatIneRowN=1; WhatIneRowN<=WhatCurL; WhatIneRowN++){
            WhereIneRowN=IneRowN+WhatIneRowN-1;

            WhereVal=Where[WhereExtRowN-1][WhereIneRowN-1];
            WhatVal=What[WhatExtRowN-1][WhatIneRowN-1];
            if(WhereVal!=WhatVal)
        }
    }
    return b;
}*/








template <typename T> std::vector<std::vector<int>> ArrayOf2DBy1D_SeekArrs2ndIn1st(T**Where, const Array2DSize& sizeWhere, T**What, const Array2DSize& sizeWhat, bool IsTransposed=false, bool EachExtRowIsReversed=false){
    std::vector<std::vector<int>>Ns;
    std::vector<int>Coords;
    bool isPresentHere;
    int whereQ=sizeWhere.GetQExtRows(), Lcur;
    for(int ExtRowN=1; ExtRowN<=whereQ; ExtRowN++){
        Lcur=sizeWhere.GetLength(ExtRowN);
        for(int IneRowN=1; IneRowN<=Lcur; IneRowN++){
            isPresentHere=ArrayOf2DBy1D_Arrs2ndIsAtPosIn1st(Where, sizeWhere, What, sizeWhat, ExtRowN, IneRowN, IsTransposed, EachExtRowIsReversed);
            if(isPresentHere){
                Coords.clear();
                Coords.push_back(ExtRowN);
                Coords.push_back(IneRowN);
                Ns.push_back(Coords);
            }
        }
    }
    return Ns;
}
template <typename T> std::vector<std::vector<int>> ArrayOf2DBy1D_SeekArrs2ndIn1st(std::vector<std::vector<T>>Where,  std::vector<std::vector<T>>What, bool IsTransposed=false, bool EachExtRowIsReversed=false){
    Array2DSize sizeWhere=GetSizeOf2DVector(Where), sizeWhat=GetSizeOf2DVector(What);
    std::vector<std::vector<int>>Ns;
    std::vector<int>Coords;
    bool isPresentHere;
    int whereQ=sizeWhere.GetQExtRows(), Lcur;
    for(int ExtRowN=1; ExtRowN<=whereQ; ExtRowN++){
        Lcur=sizeWhere.GetLength(ExtRowN);
        for(int IneRowN=1; IneRowN<=Lcur; IneRowN++){
            //isPresentHere=ArrayOf2DBy1D_Arrs2ndIsAtPosIn1st(Where, sizeWhere, What, sizeWhat, ExtRowN, IneRowN, IsTransposed, EachExtRowIsReversed);
            isPresentHere=ArrayOf2DBy1D_Arrs2ndIsAtPosIn1st(Where, What, ExtRowN, IneRowN, IsTransposed, EachExtRowIsReversed);
            if(isPresentHere){
                Coords.clear();
                Coords.push_back(ExtRowN);
                Coords.push_back(IneRowN);
                Ns.push_back(Coords);
            }
        }
    }
    return Ns;
}



template <typename T> std::vector<std::vector<int>> ArrayOf2DBy1D_SeekIneRow_Simply(T**X, const Array2DSize& size, T*row, int L, int paramN=0){
    std::vector<std::vector<T>>data;
    ArrayOf2DBy1DAssign(data, L, row, true);

}

//--------------------------------------------------
//Concat
template <typename T>void ArrayOf2DBy1DConcatAdding(T**&X1, Array2DSize& size1, T**X2, const Array2DSize& size2, int shiftVal=0, T*DefaultValParam=NULL, bool Rect=true, int QRec=0, int FromN=1){
    int L1ini=size1.GetMaxLength(), L2ini=size2.GetMaxLength(),
        //Q1ini=size1.GetQExtRows(),
        Q1=size1.GetQExtRows(),
        Q2ini=size2.GetQExtRows(),
        //Q1fin=Q1ini,
        //Q1=Q1ini;
        Q2fin, L1fin, L2fin, L;//, Lreq;
    std::vector<T>row;
    if(QRec==0){//Q to add is Q how many is given minus FromN
        Q2fin=Q2ini-FromN+1;//corr'd
    }else if(QRec<=Q2ini){
        Q2fin=QRec-FromN+1;//corr'd
    }
    if(shiftVal>=0){//2nd rows shifted to the left
        L1fin=L1ini;
        L2fin=shiftVal+L2ini;
        L = L1fin >= L2fin ? L1fin : L2fin;
        //if(QRec==0){//Q to add is Q how many is given minus FromN
        //    Q2fin=Q2ini-FromN+1;//corr'd
        //}else if(QRec<=Q2ini){
        //    Q2fin=QRec-FromN+1;//corr'd
        //}
        if(Rect){
            for(int i=1; i<=Q1; i++){
                row=Array1DGetSubArray_byLs_PreMark(X1[i-1], size1.GetLength(i), L, 1, 0, DefaultValParam);
                //template<typename T> void ArrayOf2DBy1DSetExtRowN(T**&X, Array2DSize&size, int ExtRowN, int wholeRowL=0, T*rowParam=NULL, int whatFromN=1, int QDefaultValuesBefore=0, bool KeepIfRect=true, T*DefaultValParam=NULL, bool DefaultIsSpecNotOwn=false){
                //
                ArrayOf2DBy1DSetExtRowN(X1, size1, i, row, DefaultValParam, 1, 0, Rect);
            }
            for(int i=FromN; i<=FromN+Q2fin-1; i++){
                row=Array1DGetSubArray_byLs_PreMark(X2[i-1], size2.GetLength(i), L, 1, shiftVal, DefaultValParam);//L
                ArrayOf2DBy1DAddExtRow(X1, size1, row, 1, 0, true, DefaultValParam);
            }
        }else{
            for(int i=FromN; i<=FromN+Q2fin-1; i++){
                row=Array1DGetSubArray_byLs_PreMark(X2[i-1], size2.GetLength(i), 0, 1, shiftVal, DefaultValParam);//0
                ArrayOf2DBy1DAddExtRow(X1, size1, row, 1, 0, true, DefaultValParam);
            }
        }//gut, I men os ver
    }else{//shiftVal<0 //1st rows shifted to the left sdi
        L1fin=-shiftVal+L1ini;
        L2fin=L2ini;
        L = L1fin >= L2fin ? L1fin : L2fin;
        //
        if(Rect){
            for(int i=1; i<=Q1; i++){
                row=Array1DGetSubArray_byLs_PreMark(X1[i-1], size1.GetLength(i), L, 1, shiftVal, DefaultValParam);
                ArrayOf2DBy1DSetExtRowN(X1, size1, i, row, DefaultValParam, 1, 0, true);
            }
            for(int i=FromN; i<= FromN+Q2fin-1; i++){
                row=Array1DGetSubArray_byLs_PreMark(X2[i-1], size2.GetLength(i), L, 1, 0, DefaultValParam);//L
                ArrayOf2DBy1DAddExtRow(X1, size1, row, 1, 0, true, DefaultValParam);
            }
        }else{
            for(int i=1; i<=Q1; i++){
                row=Array1DGetSubArray_byLs_PreMark(X1[i-1], size1.GetLength(i), 0, 1, shiftVal, DefaultValParam);
                //ArrayOf2DBy1DSetExtRowN(X1, size1, i, row, 1, 0, true, DefaultValParam, false);
                ArrayOf2DBy1DSetExtRowN(X1, size1, i, row, DefaultValParam, 1, 0, false);
            }
            for(int i=FromN; i<= FromN+Q2fin-1; i++){
                ArrayOf2DBy1DAddExtRow(X1, size1, X2[i-1], size2.GetLength(i), 1, 0, false, DefaultValParam);
            }
        }
    }
}
template <typename T>void ArrayOf2DBy1DConcatAdding(std::vector<std::vector<T>>&X1, std::vector<std::vector<T>>X2, int shiftVal=0, T*DefaultValParam=NULL, bool Rect=true, int QRec=0, int FromN=1){
    Array2DSize size1=GetSizeOf2DVector(X1), size2=GetSizeOf2DVector(X2);
    int L1ini=size1.GetMaxLength(), L2ini=size2.GetMaxLength(),
         //Q1ini=size1.GetQExtRows(),
         Q1=size1.GetQExtRows(),
         Q2ini=size2.GetQExtRows(),
         //Q1fin=Q1ini,
         //Q1=Q1ini;
         Q2fin, L1fin, L2fin, L;//, Lreq;
    int Q;
    std::vector<T>row;
    bool contin;
    int whereCurN, whatCurN,  whereCurL, whatCurL;
    if(QRec==0){//Q to add is Q how many is given minus FromN
        Q2fin=Q2ini-FromN+1;//corr'd
    }else if(QRec<=Q2ini){
        Q2fin=QRec-FromN+1;//corr'd
    }
    Q=Q1+Q2fin;
    if(shiftVal>=0){
         L = L1ini >= shiftVal+L2ini ? L1ini : shiftVal+L2ini;
         if(Rect){
             for(int i=1; i<=Q1; i++){
                 row=Array1DGetSubArray_byLs_PreMark(X1[i-1], L, 1, 0, DefaultValParam);
                 ArrayOf2DBy1DSetExtRowN(X1, i, row);
             }
             for(int i=FromN; i<=FromN+Q2fin-1; i++){
                 row=Array1DGetSubArray_byLs_PreMark(X2[i-1], L, 1, shiftVal, DefaultValParam);
                 ArrayOf2DBy1DAddExtRow(X1, row);
             }
         }else{
             for(int i=FromN; i<=FromN+Q2fin-1; i++){
                 row=Array1DGetSubArray_byLs_PreMark(X2[i-1], 0, 1, shiftVal, DefaultValParam);
                 ArrayOf2DBy1DAddExtRow(X1, row);
             }
         }
    }else{
        L = -shiftVal+L1ini >= L2ini ? -shiftVal+L1ini : L2ini;
        if(Rect){
            for(int i=1; i<=Q1; i++){
                row=Array1DGetSubArray_byLs_PreMark(X1[i-1], L, 1, -shiftVal, DefaultValParam);
                ArrayOf2DBy1DSetExtRowN(X1, i, row);
            }
            for(int i=FromN; i<=FromN+Q2fin-1; i++){
                row=Array1DGetSubArray_byLs_PreMark(X2[i-1], L, 1, 0, DefaultValParam);
                ArrayOf2DBy1DAddExtRow(X1, row);
            }
        }else{
            for(int i=1; i<=Q1; i++){
                row=Array1DGetSubArray_byLs_PreMark(X1[i-1], 0, 1, -shiftVal, DefaultValParam);
                ArrayOf2DBy1DSetExtRowN(X1, i, row);
            }
            for(int i=FromN; i<=FromN+Q2fin-1; i++){
                row=Array1DGetSubArray_byLs_PreMark(X2[i-1], 0, 1, 0, DefaultValParam);
                ArrayOf2DBy1DAddExtRow(X1, row);
            }
        }
    }
}

template <typename T>void ArrayOf2DBy1DConcatStretching(T**&X1, Array2DSize& size1, T**X2, const Array2DSize& size2, int shiftVal=0, T*DefaultValParam=NULL, bool Rect=true, bool IfFirstNotRect_AddToLastNotSrtetchToMax=true, int whatFromN=1, int QRec=0, TValsShowHide*vsh=NULL){
    int Lmax1=size1.GetMaxLength(), Lmin1=size1.GetMinLength(), Lmax2=size2.GetMaxLength(), Lreq, Q1=size1.GetMaxLength(), Q2=size2.GetMaxLength(), Qmin, Q, whereCurN, whatCurN, whatCurL, whereCurLini, whereCurLtmp;
    std::vector<T>row;
    T*ptr=NULL;
    T CurVal, DefaultVal;
    bool contin=false;
    if(DefaultValParam!=NULL){
        DefaultVal=(*(DefaultValParam));
    }
    if(Rect){
        Lreq=Lmax1+Lmax2;
    }else{
        Lreq=0;
    }
    if(QRec==0){
        QRec=size2.GetQExtRows()-whatFromN+1;
    }
    if(!IfFirstNotRect_AddToLastNotSrtetchToMax && Lmax1!=Lmin1){
        ArrayOf2DBy1DSetRect(X1, size1, Lmax1);
        if(vsh!=NULL){
            std::cout<<"Now all rows have same length:"<<std::endl;
            Arr2DStdCOut2D(X1, size1, " ", "; ", false);
        }
    }
    if(shiftVal>=0){//=> added row shifts down
        Qmin = Q1 <= shiftVal+QRec ? Q1 : shiftVal+QRec;
        Q = Q1 >= shiftVal+QRec ? Q1 : shiftVal+QRec;
        //
        //if(Rect){
            if(shiftVal<Q1){
                //shifted rows
                if(vsh!=NULL){
                    std::cout<<"Shifted empty rows "<<std::endl;
                }
                for(int i=1; i<=shiftVal; i++){
                    row=Array1DGetSubArray_byLs_PreMark(X1[i-1], size1.GetLength(i), Lreq, 1, 0, DefaultValParam);
                    ArrayOf2DBy1DSetExtRowN(X1, size1, i, row, DefaultValParam, 1, 0, false);
                    if(vsh!=NULL){
                        std::cout<<"row "<<i<<":";
                        Array1DShowConsole(X1[i-1], size1.GetLength(i));
                        std::cout<<std::endl;
                    }
                }
                //now placing after shift, combined rows
                if(vsh!=NULL){
                    std::cout<<"Rows combined after shift"<<std::endl;
                }
                whatCurN=whatFromN;
                whereCurN=shiftVal+whatCurN-whatFromN+1;//iqv: 3+3-3+1
                whatCurN-=1;
                whereCurN-=1;
                contin=(whatCurN<whatFromN+QRec-1 && whereCurN<Q1);
                while(contin){
                    whatCurN++;
                    whereCurN++;
                    //for(whatCurN=whatFromN; whatCurN<=whatFromN+QRec-1; whatCurN++){
                    whatCurL=size2.GetLength(whatCurN);
                    whereCurN=shiftVal+whatCurN-whatFromN+1;
                    whereCurLini=size1.GetLength(whereCurN);
                    whereCurLtmp=whereCurLini;
                    for(int j=1; j<=whatCurL; j++){
                        CurVal=X2[whatCurN-1][j-1];
                        //
                        //Array1DAddElement(X1[whereCurN-1], whereCurLtmp, CurVal);//also possible ety ably
                        //ArrayOf2DBy1DAddToExtRowN(T**&X, Array2DSize& size, int ExtRowN, T val)
                        ArrayOf2DBy1DAddToExtRowN(X1, size1, whereCurN, CurVal);

                    }
                    whereCurLini=whereCurLtmp;
                    whereCurLini=size1.GetLength(whereCurN);
                    for(int j=whereCurLini+1; j<=Lreq; j++){
                        Array1DAddElement(X1[whereCurN-1], whereCurLtmp, DefaultVal);
                    }
                    if(vsh!=NULL){
                        std::cout<<"row "<<whereCurN<<":";
                        Array1DShowConsole(X1[whereCurN-1], size1.GetLength(whereCurN));
                        std::cout<<std::endl;
                    }
                    contin=(whatCurN<whatFromN+QRec-1 && whereCurN<Q1);
                }//while
                //if(Rect){
                //    size1.Set(Q1, Lreq);
                //}
                //rest
                if(whatCurN==Q2 && whereCurN<Q1){
                    for(int i=shiftVal+QRec+1; i<=Q1; i++){
                        //row=Array1DAssign()
                        row=Array1DGetSubArray_byLs_PreMark(X1[i-1], size1.GetLength(i), Lreq, 1, 0, DefaultValParam);
                        ArrayOf2DBy1DSetExtRowN(X1, size1, i, row, DefaultValParam, 1, 0, false);
                    }
                }else if(whatCurN<Q2 && whereCurN==Q1){
                    for(int i=whatCurN+1; i<=Q2; i++){
                        row=Array1DGetSubArray_byLs_PreMark(X2[i-1], size2.GetLength(i), Lreq, 1, Lmax1, DefaultValParam);
                        ArrayOf2DBy1DAddExtRow(X1, size1, row);
                    }
                }
            }else{//shiftVal>Q1
                Q = -shiftVal+Q1 >= QRec ? -shiftVal+Q1 : QRec; //n'utf'd
                for(int i=1; i<=Q1; i++){
                    row=Array1DGetSubArray_byLs_PreMark(X1[i-1], size1.GetLength(i), Lreq, 1, 0, DefaultValParam);
                    ArrayOf2DBy1DSetExtRowN(X1, size1, i, row, DefaultValParam, 1, 0, false);
                }
                for(int i=Q1+1; i<=shiftVal; i++){
                    row=Array1DGetSubArray_byLs_PreMark(ptr, 0, Lreq, 1, 0, DefaultValParam);//ety ably
                    row=Array1DGetSubArray_byLs_PreMark(ptr, 0, Lreq, 1, Lreq, DefaultValParam);//ety ably
                    ArrayOf2DBy1DAddExtRow(X1, size1, row);
                }
                for(int i=1; i<=QRec; i++){
                    row=Array1DGetSubArray_byLs_PreMark(X2[i-1], size2.GetLength(i), Lreq, 1, Lmax1, DefaultValParam);
                    ArrayOf2DBy1DAddExtRow(X1, size1, row);
                }
                if(Rect){
                    size1.Set(Q1, Lreq);
                }
            }
        //}else{//not rect. Ma exda id ov s' uz tbi Rect et ne, Lreq fybe all cass

        //}
    }else{//shiftVal<0 //=> first row shifts down, new row starts with whatFromN. row o'X2
        //ShiftVals
        if(vsh!=NULL){
            std::cout<<"Shifted rows X2"<<std::endl;
        }
        whatCurN=whatFromN-1;
        whereCurN=0;
        //contin=(whatCurN<QRec+whatFromN-1 && whatCurN<Q1);
        contin=(whatCurN<QRec+whatFromN-1 && whatCurN-whatFromN+1<-shiftVal);//QRec is jq <=Q2
        //for(whatCurN=-shiftVal; whatCurN>=1; whatCurN--){
        //for(whatCurN=whatFromN; whatCurN<=-shiftVal+whatFromN-1; whatCurN++){
        while(contin){
            whatCurN++;
            whereCurN++;
            row=Array1DGetSubArray_byLs_PreMark(X2[whatCurN-1], size2.GetLength(whatCurN), Lreq, 1, Lmax1, DefaultValParam);
            //ArrayOf2DBy1DInsExtRow(X1, size1, 1, row, DefaultValParam, 1, 0, false); //false ob L ms'be new, S n'ved l'L's val
            ArrayOf2DBy1DInsExtRow(X1, size1, whereCurN, row, DefaultValParam, 1, 0, false); //false ob L ms'be new, S n'ved l'L's val
            //ArrayOf2DBy1DInsExtRow(X1, size1, 1, size2.GetLength(whatCurN), X2[whatCurN-1], DefaultValParam, 1, Lmax1, false);//false ob L ms'be new, S n'ved S
            contin=(whatCurN<QRec+whatFromN-1 && whatCurN-whatFromN+1<-shiftVal);
            if(vsh!=NULL){
                std::cout<<"row "<<whatCurN<<": ";
                Array1DShowConsole(X1[whereCurN-1], size1.GetLength(whereCurN));
            }
        }
        Q1=size1.GetQExtRows();
        //nu row o'X1, ic wa #1., ha Nr  -shiftVal+1;
        if(-shiftVal<whatFromN+QRec-1){ //ei rows X2 extend l'rows o'X1
            contin=(whatCurN<whatFromN+QRec-1 && whereCurN<Q1);
            while(contin){
                whatCurN++;
                whereCurN++;
                contin=(whatCurN<whatFromN+QRec-1 && whereCurN<Q1);
                if(vsh!=NULL){
                    std::cout<<"Stretching X1["<<whereCurN<<"(former "<<whereCurN+shiftVal<<")] with X2["<<whatCurN<<"]"<<std::endl;
                }
                whatCurL=size2.GetLength(whatCurN);
                for(int j=1; j<=whatCurL; j++){
                    CurVal=X2[whatCurN-1][j-1];
                    ArrayOf2DBy1DAddToExtRowN(X1, size1, whereCurN, CurVal);
                }
                whereCurLini=size1.GetLength(whereCurN);
                if(Rect){
                    for(int j=whereCurLini+1; j<=Lreq; j++){
                        ArrayOf2DBy1DAddToExtRowN(X1, size1, whereCurN, DefaultVal);
                    }
                }
                whereCurLini=size1.GetLength(whereCurN);
                if(vsh!=NULL){
                    Array1DShowConsole(X1[whereCurN-1],whereCurLini);
                }
            }
            //last rows
            if(whatCurN<whatFromN+QRec-1 && whereCurN==Q1){
                if(vsh!=NULL){
                    std::cout<<"rows of X2 remaining adding:"<<std::endl;
                }
                for(int i=whatCurN+1; i<=whatFromN+QRec-1; i++){
                    whatCurL=size2.GetLength(i);
                    row=Array1DGetSubArray_byLs_PreMark(X2[i-1], size2.GetLength(i), Lreq, 1, Lmax1, DefaultValParam);
                    ArrayOf2DBy1DAddExtRow(X1, size1, row, 1, 0, false, DefaultValParam);
                    //
                    whereCurN++;
                    whereCurLini=size1.GetLength(i);
                    if(vsh!=NULL){
                        std::cout<<"row"<<whereCurN<<" (L="<<whereCurLini<<"): ";
                        Array1DShowConsole(X1[whereCurN-1],whereCurLini);
                    }
                }
            }else if(whatCurN==whatFromN+QRec-1 && whereCurN<Q1){
                Lmin1=size1.GetMinLength();
                Lmax1=size1.GetMaxLength();
                if(Rect && Lmin1!=Lmax1){
                    if(vsh!=NULL){
                        std::cout<<"rows of X1 remaining strtetchig"<<std::endl;
                    }
                    for(int i=whereCurN+1; i<=Q1; i++){
                        whereCurLtmp=size1.GetLength(i);
                       //for(int j=whereCurLtmp+1; j<=Lreq; j++){
                        for(int j=whereCurLtmp+1; j<=Lreq; j++){
                            ArrayOf2DBy1DAddToExtRowN(X1, size1, i, DefaultVal);
                        }
                        //
                        //whereCurN++;
                       // whereCurLini=size1.GetLength(i);
                        if(vsh!=NULL){
                            std::cout<<"row"<<i<<" (L="<<whereCurLini<<"): ";
                            Array1DShowConsole(X1[i-1],whereCurLini);
                        }
                    }
                }else{
                    //else NOp, tic rows stay in X1 arr co cu lenis.
                    if(vsh!=NULL){
                        std::cout<<"rows of X1 remain without changes"<<std::endl;
                    }
                }
            }//else NOp, tbi nablb zq, ob so n'fin while, so tdi s'eq, sdi tbi arrs gead fin
        }
    }//shiftvals
    if(Rect){
        size1.SetIneRowsLength(-1);
    }
    if(vsh!=NULL){
        std::cout<<"Concat stretching finishes working, now we have "<<Q1<<" rows"<<std::endl;
    }
}//fn concat stretching
template <typename T>void ArrayOf2DBy1DConcatStretching(std::vector<std::vector<T>>&X1,std::vector<std::vector<T>>X2,  int shiftVal=0, T*DefaultValParam=NULL, bool Rect=true, bool IfFirstNotRect_AddToLastNotSrtetchToMax=true, int whatFromN=1, int QRec=0, TValsShowHide*vsh=NULL){
    Array2DSize size1=GetSizeOf2DVector(X1), size2=GetSizeOf2DVector(X2);
    int Q1=size1.GetQExtRows(), Q2=size2.GetQExtRows(), Lmax1=size1.GetMaxLength(), Lmax2=size2.GetMaxLength(), Lmin1=size1.GetMinLength(), Lmin2=size2.GetMinLength();
    int whatCurN, whereCurN, whatCurL, whereCurL, Lreq, minN;
    bool contin;
    T DefaultVal, CurVal;
    std::vector<T> row;
    if(DefaultValParam!=NULL){
        DefaultVal=(*(DefaultValParam));
    }
    if(Rect){
        Lreq=Lmax1+Lmax2;
    }else{
        Lreq=0;
    }
    if(QRec==0){
        QRec=Q2-whatFromN+1;
    }
    if(shiftVal>=0){
        //first shiftVal rows of X1
        minN= shiftVal<=Q1 ? shiftVal: Q1;
        for(int i=Q1+1; i<=shiftVal; i++){//if shiftVal>Q1, else do ce 0 times
            ArrayOf2DBy1DAddExtRow(X1, row, 1, 0, false, DefaultValParam);
        }
        //combined
        whereCurN=shiftVal+1;
        if(shiftVal<Q1){
            whatCurN=whatFromN;
            whatCurN--;
            whereCurN--;
            contin=(whatCurN<whatFromN+QRec-1 && whereCurN<Q1);
            while(contin){
                whatCurN++;
                whereCurN++;
                whatCurL=X2[whatCurN-1].size();
                for(int j=1; j<=whatCurL; j++){
                    CurVal=X2[whatCurN-1][j-1];
                    ArrayOf2DBy1DAddToExtRowN(X1, whereCurN, CurVal);
                }
            }
        }
        //after combined
        if(whatCurN<whatFromN+QRec-1 && whereCurN==Q1){
            for(int i=whatCurN+1; i<=whatFromN+QRec-1; i++){
                ArrayOf2DBy1DAddExtRow(X1, X2[i-1], 1, 0, false, DefaultValParam);
            }
        }//else NOp
        if(Rect){
            Q1=X1.size();
            for(int i=1; i<=Q1; i++){
                whereCurL=X1[i-1].size();
                if(whereCurL<Lmax1+Lmax2){
                    row=Array1DGetSubArray_byLs_PreMark(X1[i-1], Lreq, 1, 0, DefaultValParam);
                    //ArrayOf2DBy1DSetExtRowN(std::vector<std::vector<T>>&X, int ExtRowN, std::vector<T>rowParam, int whatFromN=1, int QDefaultValuesBefore=0, bool KeepIfRect=true, T*DefaultValParam=NULL, bool DefaultIsSpecNotOwn=false){//, int ifDiff_ByArr0ByNew1Min2Max3=0){
                    ArrayOf2DBy1DSetExtRowN(X1, i, row, 1, 0, false, DefaultValParam, false);
                }
            }
        }
    }else{//shiftVal<0
        //shifted X2
        minN = -shiftVal <= whatFromN+QRec-1  ? -shiftVal : whatFromN+QRec-1;
        for(int i=1; i<=minN; i++){
            whatCurN=whatFromN+i-1;
            row=Array1DGetSubArray_byLs_PreMark(X2[whatCurN-1], Lreq, 1, Lmax1, DefaultValParam);
            ArrayOf2DBy1DInsExtRow(X1, whatCurN, X2[whatCurN-1]);
        }
        for(int i=whatFromN+QRec-1; i<=-shiftVal; i++){//if shiftVal>whatFromN+QRec-1, else do ce 0 times
            ArrayOf2DBy1DAddExtRow(X1, row, 1, 0, false, DefaultValParam);
        }
        //combined rows
        whatCurN=-shiftVal+1;
        whereCurN=-shiftVal+1;//ob tic X2 were inserted to X1, before own rows of X1
        whatCurN--;
        whereCurN--;
        contin=(whatCurN<whatFromN+QRec-1 && whereCurN<Q1);
        while(contin){
            whatCurN++;
            whereCurN++;
            whatCurL=X2[whatCurN-1].size();
            for(int j=1; j<=whatCurL; j++){
                CurVal=X2[whatCurN-1][j-1];
                ArrayOf2DBy1DAddToExtRowN(X1, whereCurN, CurVal);
            }
            contin=(whatCurN<whatFromN+QRec-1 && whereCurN<Q1);
        }
        //rest
        if(whatCurN<whatFromN+QRec-1 && whereCurN==Q1){
             for(int i=whatCurN+1; i<=whatFromN+QRec-1; i++){
                 row=Array1DGetSubArray_byLs_PreMark(X2[i-1], Lreq, 1, Lmax1, DefaultValParam);
                 ArrayOf2DBy1DAddExtRow(X1, row);
             }
        }//else if(whatCurN==whatFromN+QRec-1 && whereCurN<Q1){ NOp
        //if Rect
        Q1=X1.size();
        if(Rect){
            for(int i=1; i<=Q1; i++){
                whereCurL=X1[i-1].size();
                if(whereCurL<Lmax1+Lmax2){
                    ArrayOf2DBy1DSetExtRowN(X1, i, Lreq, 1, 0, DefaultValParam, false);
                }
            }
        }
    }
}

//-------------------------------------------------

template <typename T>void ArrayOf1DStdCOut(T*x, int Q,  QString delimElem="; ", QString delimLines="; ", bool useBrackets=true){
    T*ptr;
    if(useBrackets){
        std::cout<<"[";
    }
    for(int i=1; i<=Q; i++){
        ptr=(x+i-1);
        std::cout<<ptr;
        std::cout<<delimElem.toStdString().c_str();
    }
    if(useBrackets){
        std::cout<<"]";
    }
}
// ^ os ArrayOf1DStdCOut , not Array1DStdCOut - uz ArrayOf2Dby1DStdCOut...

template <typename T>void ArrayOf2Dby1DStdCOutSingleExtRow(T*X, const Array2DSize&size, int rowN, QString delimElem="; ", bool useBrackets=true){
    int PosN, L,
            Q=size.GetQExtRows();//f'ext rows no need
    T val;
    PosN=Array2DSizeBasedOn1D_PosByCoords(rowN, 1, size.GetLengthes());
    if(PosN!=0){
        L=size.GetLength(rowN);
        if(useBrackets){
            std::cout<<"[";
        }
        for(int i=1; i<=L-1; i++){
            PosN=Array2DSizeBasedOn1D_PosByCoords(rowN, i, size.GetLengthes());
            val=X[PosN-1];
            std::cout<<val;
            std::cout<<delimElem.toStdString().c_str();
        }
        PosN=Array2DSizeBasedOn1D_PosByCoords(rowN, L, size.GetLengthes());
        val=X[PosN-1];
        std::cout<<val;
        if(useBrackets){
            std::cout<<"]";
        }
    }
}
template <typename T> void ArrayOf2Dby1DStdCOutSingleExtRow1(T*X, const Array2DSize&size, int rowN, QString delimElem="; ", bool useBrackets=true){
    int PosN, L,
              Q=size.GetQExtRows();//f'ext rows no need
    T*ptr;
    PosN=Array2DSizeBasedOn1D_PosByCoords(rowN, 1, size.GetLengthes());
    if(PosN!=0){
        if(useBrackets)std::cout<<"[";
        L=size.GetLength(rowN);
        ptr=X[PosN-1];
        Array1DShowConsole(ptr, L, delimElem.toStdString());
        if(useBrackets)std::cout<<"]";
    }
}
template <typename T>void ArrayOf2Dby1DStdCOutSingleExtRow(std::vector<T>X, const Array2DSize&size, int rowN, QString delimElem="; ", bool useBrackets=true){
    ArrayOf2Dby1DStdCOutSingleExtRow(X.data(), size, rowN, delimElem, useBrackets);
}

template <typename T>void ArrayOf2Dby1DStdCOutSingleIneRow(T*X, const Array2DSize&size, int rowN, QString delimElem="; ", bool useBrackets=true){
    int PosN, L,
            Q=size.GetQExtRows(), Lmin=size.GetMinLength(), Lmax=size.GetMaxLength();
    T val;
    PosN=Array2DSizeBasedOn1D_PosByCoords(1, rowN, size.GetLengthes());
    if(PosN!=0 && rowN<=Lmin){
        L=size.GetLength(rowN);
        if(useBrackets){
            std::cout<<"[";
        }
        for(int i=1; i<=L-1; i++){
            PosN=Array2DSizeBasedOn1D_PosByCoords(i, rowN, size.GetLengthes());
            val=X[PosN-1];
            std::cout<<val;
            std::cout<<delimElem.toStdString().c_str();
        }
        PosN=Array2DSizeBasedOn1D_PosByCoords(rowN, L, size.GetLengthes());
        val=X[PosN-1];
        std::cout<<val;
        if(useBrackets){
            std::cout<<"]";
        }
    }
}
template <typename T>void ArrayOf2Dby1DStdCOutSingleIneRow(std::vector<T>X, const Array2DSize&size, int rowN, QString delimElem="; ", bool useBrackets=true){
     ArrayOf2Dby1DStdCOutSingleIneRow(X.data(),size, rowN, delimElem, useBrackets);
}

template <typename T>void ArrayOf2Dby1DStdCOut2D(T*X, const Array2DSize&size, QString delimElem="; ", QString delimLines="; ", bool useBrackets=true){
    int Q=size.GetQExtRows(), PosN, L;
    //template <typename T>void ArrayOf2Dby1DStdCOutSingleExtRow(T*X, const Array2DSize&size, int rowN, QString delimElem="; ", bool useBrackets=true){
    if(useBrackets)std::cout<<"[";
    ArrayOf2Dby1DStdCOutSingleExtRow(X, size, 1, delimElem, useBrackets);
    std::cout<<delimLines.toStdString().c_str();
    std::cout<<std::endl;
    for(int i=2; i<=Q-1; i++){
        if(useBrackets)std::cout<<" ";
        ArrayOf2Dby1DStdCOutSingleExtRow(X, size, i, delimElem, useBrackets);
        std::cout<<delimLines.toStdString().c_str();
        std::cout<<std::endl;
    }
    if(useBrackets)std::cout<<" ";
    ArrayOf2Dby1DStdCOutSingleExtRow(X, size, Q, delimElem, useBrackets);
    std::cout<<delimLines.toStdString().c_str();
    if(useBrackets)std::cout<<"]";
    std::cout<<std::endl;
}
template <typename T>void ArrayOf2Dby1DStdCOut2D(std::vector<T>X, const Array2DSize&size, QString delimElem="; ", QString delimLines="; ", bool useBrackets=true){
    ArrayOf2Dby1DStdCOut2D(X.data, size, delimElem, delimLines, useBrackets);
}

template <typename T>void ArrayOf2Dby1DStdCOut1D(T*X, const Array2DSize&size, QString delimElem="; ", QString delimLines="; ", bool useBrackets=true){
    int Q=size.GetQExtRows(), PosN, L;
    //template <typename T>void ArrayOf2Dby1DStdCOutSingleExtRow(T*X, const Array2DSize&size, int rowN, QString delimElem="; ", bool useBrackets=true){
    if(useBrackets)std::cout<<"[";
    ArrayOf2Dby1DStdCOutSingleExtRow(X, size, 1, delimElem, useBrackets);
    std::cout<<delimLines.toStdString().c_str();
    //std::cout<<std::endl;
    for(int i=2; i<=Q-1; i++){
        //if(useBrackets)std::cout<<" ";
        ArrayOf2Dby1DStdCOutSingleExtRow(X, size, i, delimElem, useBrackets);
        std::cout<<delimLines.toStdString().c_str();
        //std::cout<<std::endl;
    }
    //if(useBrackets)std::cout<<" ";
    ArrayOf2Dby1DStdCOutSingleExtRow(X, size, Q, delimElem, useBrackets);
    std::cout<<delimLines.toStdString().c_str();
    if(useBrackets)std::cout<<"]";
    //std::cout<<std::endl;
}
template <typename T>void ArrayOf2Dby1DStdCOut1D(std::vector<T>X, const Array2DSize&size, QString delimElem="; ", QString delimLines="; ", bool useBrackets=true){
    ArrayOf2Dby1DStdCOut1D(X.data(), size, delimElem, delimLines, useBrackets);
}

//end of Array2DBy1D


//-----------------------------------------------------------------------------------------------------

int CalcElementNIfNegativeOf(int N, int Q);

//-----------------------------------------------------------------------------------------------------

template<typename T> class Array2D_PS{
    Array2DSize size;
    T**data;
public:
    //Array2DPS(int QExtRows=1, int IneRowsLength=1, bool VarNotRect=true, T*DefaultValParam=NULL);

    Array2D_PS(const Array2DSize&size, T**X=NULL, T*DefaultValParam=NULL);
    Array2D_PS(std::vector<std::vector<T>>X, bool RectIfSameLengthes=true);
    Array2D_PS(int QExtRows=1, int IneRowsLength=1, bool RectNotVar=true, T**X=NULL, T*DefaultValParam=NULL);
    Array2D_PS(int Q, bool QRowsNotLengthOfOneExt=true, T*X=NULL, bool RectNotVar=true, T*DefaultValParam=NULL);
    Array2D_PS(std::vector<T>X, bool QRowsNotLengthOfOneExt=true, bool RectNotVar=true);
    Array2D_PS(const Array2D_PS&X);
    ~Array2D_PS();

    //
    void Assign(const Array2D_PS&X);
    Array2D_PS& operator = (const Array2D_PS&obj);
    void SetNull();
    //
    void SetNullIni();
    void Construct();
    //
    const std::vector<T> operator[](int index) const;
    //
    //void Set(int QExtRows=1, int IneRowsLength=1, bool VarNotRect=true, T*DefaultValParam=NULL);
    void Set(const Array2DSize&size, T**X=NULL, T*DefaultValParam=NULL);
    void Set(std::vector<std::vector<T>>X, bool RectIfSameLengthes=true);
    void Set(int QExtRows=1, int IneRowsLength=1, bool RectNotVar=true, T**X=NULL, T*DefaultValParam=NULL);
    void Set(int Q, bool QRowsNotLengthOfOneExt=true, T*X=NULL, bool RectNotVar=true, T*DefaultValParam=NULL);
    void Set(std::vector<T>X, bool QRowsNotLengthOfOneExt=true, bool RectNotVar=true);
    //
    void SetSize(const Array2DSize&newSize, T*DefaultValParam=NULL);
    void SetSize(int QExtRows=1, int IneRowsLength=1, bool VarNotRect=true, T*DefaultValParam=NULL);
    //
    std::vector<std::vector<T>> Get();
    std::vector<T> GetExtRowN_asVect(int ExtRowN);
    std::vector<T> GetIneRowN_asVect(int IneRowN);
    T*GetExtRowN_asPtr(int ExtRowN);
    //
    Array2DSize GetSize() const;
    int GetQExtRows() const;
    int GetLength(int ineRowN) const;
    int GetMinLength() const;
    int GetMaxLength() const;
    //
    //void SetVal(T val, int ExtRowN, int IneRowN=1);
    //
    T GetElement_AsVal(int ExtRowN, int IneRowN=1) const;
    T*GetElement_AsPtr(int ExtRowN, int IneRowN=1) const;
    //
    void SetElement(const T& val, int ExtRowN, int IneRowN=1);
    //
    void SwapVals(int ExtRowN1, int IneRowN1, int ExtRowN2, int IneRowN2);
    void SwapExtRows(int ExtRowN1, int ExtRowN2);
    void SwapIneRows(int IneRowN1, int IneRowN2);
    void ReverseExtRows();
    void ReverseIneRows();
    void ReverseExtRowN(int ExtRowN);
    void ReverseIneRowN(int IneRowN);
    void Transpose();
    //void TransposeStruct();
    //
    void DelExtRowN(int ExtRowN);
    void DelIneRowN(int IneRowN, bool ifPos0IgnoreNotDelLastOfVaria=true);
    //
    void SetExtRow(int ExtRowN, int wholeRowL=0, T*rowParam=NULL, int whatFromN=1, int QDefaultValuesBefore=0, bool KeepIfRect=true, T*DefaultValParam=NULL, bool DefaultIsSpecNotOwn=false);
    void SetExtRow(int ExtRowN, std::vector<T>rowParam, T*DefaultValParam=NULL, int FromN=1, int QDefaultValuesBefore=0, bool KeepIfRect=true, bool DefaultIsSpecNotOwn=false);

    void SetIneRow(int IneRowN, T*rowParam=NULL, int wholeRowL=0, int whatFromN=1, int QDefaultValsBefore=0, T*DefaultValParam=NULL, bool ExistingInsteadOfDefault=true, bool ExistingOnlyForGTLminNotIgnore=false, bool LastPossIfNotRectAtPos0NotIgnore=false);
    void SetIneRow(int IneRowN,std::vector<T>rowParam, int whatFromN=1, int QDefaultValsBefore=0, T*DefaultValParam=NULL, bool ExistingInsteadOfDefault=true, bool ExistingOnlyForGTLminNotIgnore=false, bool LastPossIfNotRectAtPos0NotIgnore=false);

    void AddExtRow(T*rowParam, int wholeRowL=0, int whatFromN=1, int QDefaultValsBefore=0, bool RectNotVar=true, T*DfltValParam=NULL, TValsShowHide*vsh=NULL);
    void AddExtRow(std::vector<T>rowParam, int whatFromN=1, int QDefaultValsBefore=0, bool RectNotVar=true, T*DfltValParam=NULL, TValsShowHide*vsh=NULL);

    void InsExtRow(int ExtRowN, int wholeRowL=0, T*rowParam=NULL, T*DefaultValParam=NULL, int whatFromN=1, int QDefaultValuesBefore=0, bool KeepIfRect=true, TValsShowHide*vsh=NULL);
    void InsExtRow(int ExtRowN, std::vector<T>rowParam, T*DefaultValParam=NULL, int whatFromN=1, int QDefaultValuesBefore=0, bool KeepIfRect=true);

    void AddIneRow(T*rowParam=NULL, int wholeRowL=0, int IfNonRectIgnore0Add1Stretch2=0, int whatFromN=1, int QDefaultValsBefore=0, T*DefaultValParam=NULL, bool RectNotVariaIfFirst=true);
    void AddIneRow(std::vector<T>rowParam, int IfNonRectIgnore0Add1Stretch2=0, int whatFromN=1, int QDefaultValsBefore=0, T*DefaultValParam=NULL);

    void InsIneRow(int IneRowPosN, T*rowParam=NULL, int wholeRowL=0, int IfAfterLminIgnore0Stretch2=0, int whatFromN=1, int QDefaultValsBefore=0, T*DefaultValParam=NULL);
    void InsIneRow(int IneRowPosN, std::vector<T>rowParam, int IfAfterLminIgnore0Stretch2=0, int whatFromN=1, int QDefaultValsBefore=0, T*DefaultValParam=NULL);
    //
    //template<typename T>Array2DAddToExtRowN(T**&X, Array2DSize& size, int ExtRowN, T val){
    void AddToExtRow(const T& val, int ExtRowN);
    void InsToExtRow(const T& val, int ExtRowN, int PosN);
    void DelFromExtRow(int ExtRowN, int PosN);
    //
    std::vector<std::vector<int>> SeekVal(const T& Val, int paramN=0);

    std::vector<std::vector<int>> SeekExtRow(T*row, int L, int paramN=0);
    std::vector<std::vector<int>> SeekExtRow(std::vector<T>row, int paramN=0);

    //std::vector<std::vector<int>> SeekIneRow(T*row, int L){
    //
    //}
    //std::vector<std::vector<int>> SeekIneRow(std::vector<T>row);
    std::vector<std::vector<int>> SeekSubArray(T**What, const Array2DSize& sizeWhat, bool IsTransposed=false, bool EachExtRowIsReversed=false);
    //std::vector<std::vector<int>> SeekSubArray(const std::vector<std::vector<T>>&data, bool Transposed, bool ExtReversed, bool IneReversed);
    std::vector<std::vector<int>>SeekSubArray(const Array2D_PS& data, bool Transposed=false, bool EachExtRowIsReversed=false);
    Array2D_PS& GetSubArray();

    //
    void ConcatAdding();

    void ConcatStretching();


    void Arr2DStdCOut2D(QString delimElem="; ", QString delimLines="; ", bool useBrackets=true);
    void Arr2DStdCOut1D(QString delimElem=", ", QString delimLines="; ", bool useBrackets=true);
};

template<typename T> Array2D_PS<T>::Array2D_PS(const Array2DSize&size, T**X, T*DefaultValParam){
    this->Construct();
    this->Set(size, X, DefaultValParam);
}
template<typename T> Array2D_PS<T>::Array2D_PS(std::vector<std::vector<T>>X, bool RectIfSameLengthes){
    this->Construct();
    this->Set(X, RectIfSameLengthes);
}
template<typename T> Array2D_PS<T>::Array2D_PS(int QExtRows, int IneRowsLength, bool RectNotVar, T**X, T*DefaultValParam){
    this->Construct();
    this->Set(QExtRows, IneRowsLength, RectNotVar, X, DefaultValParam);
}
template<typename T> Array2D_PS<T>::Array2D_PS(int Q, bool QRowsNotLengthOfOneExt, T*X, bool RectNotVar, T*DefaultValParam){
    this->Construct();
    this->Set(Q, QRowsNotLengthOfOneExt, X, RectNotVar, DefaultValParam);
}
template<typename T> Array2D_PS<T>::Array2D_PS(std::vector<T>X, bool QRowsNotLengthOfOneExt, bool RectNotVar){
    this->Construct();
    this->Set(X, QRowsNotLengthOfOneExt, RectNotVar);
}
template<typename T> Array2D_PS<T>::Array2D_PS(const Array2D_PS&obj){
    this->Construct();
    this->Assign(obj);
}

template<typename T> Array2D_PS<T>::~Array2D_PS(){
    this->SetNull();
}


template<typename T> void Array2D_PS<T>::SetNullIni(){
    this->data=NULL;
    this->size.SetNull();
}

template<typename T> void Array2D_PS<T>::Construct(){
    this->SetNullIni();
}

template<typename T> void Array2D_PS<T>::Assign(const Array2D_PS&obj){
    this->SetNull();
    this->Set(obj.size, obj.data);
}

template<typename T> Array2D_PS<T>& Array2D_PS<T>::operator = (const Array2D_PS&obj){
    this->Assign(obj);
    return *this;
}

template<typename T> void Array2D_PS<T>::SetNull(){
    int Q=this->size.GetQExtRows(), L=0;
    if(this->data!=NULL){
        for(int i=1; i<=Q; i++){
            delete[]this->data[i-1];
        }
        delete[]this->data;
        this->data=NULL;
    }
    this->size.SetNull();
}



template<typename T> const std::vector<T> Array2D_PS<T>::operator[](int index) const{
    std::vector<T>row;
    int L;
    if(this->data!=NULL && index>=0 && index<this->size.GetQExtRows() && this->data[index]!=NULL){
        L=this->size.GetLength(index+1);
        for(int i=1; i<=L; i++){
            row.push_back(this->data[index][i-1]);
        }
    }
    return row;
}

template<typename T> void Array2D_PS<T>::Set(const Array2DSize&size, T**X, T*DefaultValParam){
    int Q=size.GetQExtRows(), L;
    //
    this->SetNull();
    //
    this->data=new T*[Q];
    for(int i=1; i<=Q; i++){
        L=size.GetLength(i);
        if(L>0){
            this->data[i-1]=new T[L];
        }
    }
    if(X!=NULL){
        for(int i=1; i<=Q; i++){
            L=size.GetLength(i);
            if(L>0){
                for(int j=1; j<=L; j++){
                    this->data[i-1][j-1]=X[i-1][j-1];
                }
            }else{
                this->data[i-1]=NULL;
            }
        }
    }else{
        if(DefaultValParam!=NULL){
            for(int i=1; i<=Q; i++){
                L=size.GetLength(i);
                if(L>0){
                    for(int j=1; j<=L; j++){
                        this->data[i-1][j-1]=(*(DefaultValParam));
                    }
                }
            }
        }
    }
    this->size=size;
}
template<typename T> void Array2D_PS<T>::Set(std::vector<std::vector<T>>X, bool RectIfSameLengthes){
    int Q, L;
    this->SetNull();
    this->size=GetSizeOf2DVector(X);
    Q=this->size.GetQExtRows();
    //
    this->SetNull();
    //
    if(Q>0){
        this->data=new T*[Q];
        for(int i=1; i<=Q; i++){
            L=size.GetLength(i);
            if(L>0){
                this->data[i-1]=new T[L];
            }else{
                this->data[i-1]=NULL;
            }
        }
        //
        for(int i=1; i<=Q; i++){
            L=size.GetLength(i);
            if(L>0){
                for(int j=1; j<=L; j++){
                    this->data[i-1][j-1]=X[i-1][j-1];
                }
            }else{
                this->data[i-1]=NULL;
            }
        }
    }
}

template<typename T> void Array2D_PS<T>::Set(int QExtRows, int IneRowsLength, bool RectNotVar, T**X, T*DefaultValParam){
    this->SetNull();
    if(RectNotVar){
        this->size.Set(QExtRows, IneRowsLength);
    }else{
        this->size.SetVaria(QExtRows, IneRowsLength);
    }
    this->Set(this->size, X, DefaultValParam);
}
template<typename T> void Array2D_PS<T>::Set(int Q, bool QRowsNotLengthOfOneExt, T*X, bool RectNotVar, T*DefaultValParam){
    this->SetNull();
    if(QRowsNotLengthOfOneExt){
        Array2DAddIneRow(this->data, this->size, X, Q, 1, 1, 0, DefaultValParam, RectNotVar);
    }else{
        Array2DAddExtRow(this->data, this->size, X, Q, 1, 0, RectNotVar, DefaultValParam);
    }
}

template<typename T> void Array2D_PS<T>::Set(std::vector<T>X, bool QRowsNotLengthOfOneExt,  bool RectNotVar){
    this->SetNull();
    this->Set(X.size(), QRowsNotLengthOfOneExt, X.data(), RectNotVar);
}


template<typename T>void  Array2D_PS<T>::SetSize(const Array2DSize&newSize, T*DefaultValParam){
    Array2DSetSize(this->data, this->size, newSize, DefaultValParam);
}

template<typename T>void Array2D_PS<T>::SetSize(int QExtRows, int IneRowsLength, bool VarNotRect, T*DefaultValParam){
    Array2DSize newSize;
    if(VarNotRect){
        newSize.SetVaria(QExtRows, IneRowsLength);
    }else{
        newSize.Set(QExtRows, IneRowsLength);
    }
    Array2DSetSize(this->data, this->size, newSize, DefaultValParam);
}

template<typename T>std::vector<std::vector<T>> Array2D_PS<T>:: Get(){
    std::vector<std::vector<T>>Y;
    std::vector<T>row;
    T val;
    int Q, L;
    for(int i=1; i<=Q; i++){
        row.clear();
        L=this->size.GetLength(i);
        if(L>0){
            for(int j=1; j<=L; j++){
                val=this->data[i-1][j-1];
                row.push_back(val);
            }
        }
        Y.push_back(row);
    }
    return Y;
}

template<typename T>std::vector<T> Array2D_PS<T>::GetExtRowN_asVect(int ExtRowN){
    std::vector<T> row;
    T val;
    int Q=this->size.GetQExtRows(), L;
    if(ExtRowN>=1 && ExtRowN<=Q){
        L=this->size.GetLength(ExtRowN);
        for(int i=1; i<=L; i++){
            val=this->data[ExtRowN-1][i-1];
            row.push_back(val);
        }
    }
    return row;
}

template<typename T>std::vector<T>  Array2D_PS<T>::GetIneRowN_asVect(int IneRowN){
    std::vector<T> row;
    T val;
    int Q=this->size.GetQExtRows(), L, Lmin=size.GetMinLength(), Lmax=size.GetMaxLength();
    if(Q>0 && IneRowN>=1){
        if(IneRowN<=Lmin){
            for(int i=1; i<=Q; i++){
                val=this->data[i-1][IneRowN-1];
                row.push_back(val);
            }
        }else if(IneRowN<=Lmax){

        }
    }
    return row;
}
template<typename T>T* Array2D_PS<T>::GetExtRowN_asPtr(int ExtRowN){
    T*R=NULL;
    int Q=this->size.GetQExtRows();
    if(ExtRowN<0){
        ExtRowN=Q+ExtRowN+1;
    }
    if(ExtRowN>=1 && ExtRowN<=Q){
        R=this->data[ExtRowN-1];
    }
    return R;
}

//
template<typename T>Array2DSize  Array2D_PS<T>::GetSize() const{
    return this->size;
}

template<typename T>int  Array2D_PS<T>::GetQExtRows() const{
    return this->size.GetQExtRows();
}

template<typename T> int  Array2D_PS<T>::GetLength(int ineRowN) const{
    return this->size.GetLength(ineRowN);
}

template<typename T> int  Array2D_PS<T>::GetMinLength() const{
    return this->size.GetMinLength();
}
template<typename T> int  Array2D_PS<T>::GetMaxLength() const{
    return this->size.GetMaxLength();
}

//
//template<typename T>void  Array2D_PS<T>::SetVal(T val, int ExtRowN, int IneRowN){
//    int Q=this->size.GetQExtRows(), L=this->size.GetLength(ExtRowN);
//    if(ExtRowN>=1 && ExtRowN<=Q && IneRowN>=1 && IneRowN<=L){
//        this->data[ExtRowN-1][IneRowN-1]=val;
//    }
//}

//
template<typename T>T  Array2D_PS<T>::GetElement_AsVal(int ExtRowN, int IneRowN) const{
    T val;
    int Q=this->size.GetQExtRows(), L=this->size.GetLength(ExtRowN);
    if(ExtRowN>=1 && ExtRowN<=Q && IneRowN>=1 && IneRowN<=L){
        val=this->data[ExtRowN-1][IneRowN-1];
    }
    return val;
}
template<typename T>T* Array2D_PS<T>::GetElement_AsPtr(int ExtRowN, int IneRowN) const{
    T*ptr;
    int Q=this->size.GetQExtRows(), L=this->size.GetLength(ExtRowN);
    if(ExtRowN>=1 && ExtRowN<=Q && IneRowN>=1 && IneRowN<=L){
        ptr=this->data[ExtRowN-1][IneRowN-1];
    }
    return ptr;
}


template<typename T> void Array2D_PS<T>::SetElement(const T& val, int ExtRowN, int IneRowN){
    int Q=this->GetQExtRows(), L=this->GetLength(ExtRowN);
    if(ExtRowN<0){
        ExtRowN=ExtRowN+Q+1;
    }
    if(IneRowN<0){
        IneRowN=IneRowN+L+1;
    }
    if(ExtRowN>=1 && ExtRowN<=Q && IneRowN>=1 && IneRowN<=L){
        this->data[ExtRowN-1][IneRowN-1]=val;
    }
}

template<typename T>void Array2D_PS<T>::SwapVals(int ExtRowN1, int IneRowN1, int ExtRowN2, int IneRowN2){
    Array2DSwapVals(this->data, ExtRowN1, IneRowN1, ExtRowN2, IneRowN2);
}

template<typename T>void Array2D_PS<T>::SwapExtRows(int ExtRowN1, int ExtRowN2){
    Array2DSwapExtRows(this->data, this->size, ExtRowN1, ExtRowN2);
}

template<typename T>void Array2D_PS<T>::SwapIneRows(int IneRowN1, int IneRowN2){
    Array2DSwapIneRows(this->data, IneRowN1, IneRowN2);
}

template<typename T>void Array2D_PS<T>::ReverseExtRows(){
    Array2DReverseExtRows(this->data, this->size);
}

template<typename T>void Array2D_PS<T>::ReverseIneRows(){
    Array2DReverseIneRows(this->data, this->size);
}

template<typename T>void Array2D_PS<T>::ReverseExtRowN(int ExtRowN){
    Array2DReverseExtRowN(this->data, this->size, ExtRowN);
}

template<typename T>void Array2D_PS<T>::ReverseIneRowN(int IneRowN){
    Array2DReverseIneRowN(this->data, this->size, IneRowN);
}

template<typename T>void Array2D_PS<T>::Transpose(){
    Array2DTranspose(this->data, this->size);
}

template<typename T>void Array2D_PS<T>::DelExtRowN(int ExtRowN){
    Array2DDelExtRow(this->data, this->size, ExtRowN);
}

template<typename T>void Array2D_PS<T>::DelIneRowN(int IneRowN, bool ifPos0IgnoreNotDelLastOfVaria){
    Array2DDelIneRow(this->data, this->size, IneRowN, ifPos0IgnoreNotDelLastOfVaria);
}

template<typename T> void Array2D_PS<T>::SetExtRow(int ExtRowN, int wholeRowL, T*rowParam, int whatFromN, int QDefaultValuesBefore, bool KeepIfRect, T*DefaultValParam, bool DefaultIsSpecNotOwn){
    //template<typename T> void Array2DSetExtRowN(T**&X, Array2DSize&size, int ExtRowN, int wholeRowL=0, T*rowParam=NULL, int whatFromN=1, int QDefaultValuesBefore=0, bool KeepIfRect=true, T*DefaultValParam=NULL, bool DefaultIsSpecNotOwn=false)
    Array2DSetExtRowN(this->data, this->size, ExtRowN, wholeRowL, rowParam, whatFromN, QDefaultValuesBefore, KeepIfRect, DefaultValParam, DefaultIsSpecNotOwn);
}
template<typename T> void Array2D_PS<T>::SetExtRow(int ExtRowN, std::vector<T>rowParam, T*DefaultValParam, int FromN, int QDefaultValuesBefore, bool KeepIfRect, bool DefaultIsSpecNotOwn){
    //template<typename T> void Array2DSetExtRowN(T**&X, Array2DSize&size, int ExtRowN, std::vector<T>rowParam, T*DefaultValParam=NULL, int FromN=1, int QDefaultValuesBefore=0, bool KeepIfRect=true, bool DefaultIsSpecNotOwn=false){
    Array2DSetExtRowN(this->data, this->size, ExtRowN, rowParam, DefaultValParam, FromN, QDefaultValuesBefore, KeepIfRect, DefaultIsSpecNotOwn);
}

template<typename T> void Array2D_PS<T>::SetIneRow(int IneRowN, T*rowParam, int wholeRowL, int whatFromN, int QDefaultValsBefore, T*DefaultValParam, bool ExistingInsteadOfDefault, bool ExistingOnlyForGTLminNotIgnore, bool LastPossIfNotRectAtPos0NotIgnore){
    //template<typename T>Array2DSetIneRow(T**&X, Array2DSize& size, int IneRowN, T*rowParam=NULL, int wholeRowL=0, int whatFromN=1, int QDefaultValsBefore=0, T*DefaultValParam=NULL, bool ExistingInsteadOfDefault=true, bool ExistingOnlyForGTLminNotIgnore=false, bool LastPossIfNotRectAtPos0NotIgnore=false){
    Array2DSetIneRow(this->data, this->size, IneRowN, rowParam, wholeRowL, whatFromN, QDefaultValsBefore, DefaultValParam, ExistingInsteadOfDefault, ExistingOnlyForGTLminNotIgnore, LastPossIfNotRectAtPos0NotIgnore);
}
template<typename T> void Array2D_PS<T>::SetIneRow(int IneRowN,std::vector<T>rowParam, int whatFromN, int QDefaultValsBefore, T*DefaultValParam, bool ExistingInsteadOfDefault, bool ExistingOnlyForGTLminNotIgnore, bool LastPossIfNotRectAtPos0NotIgnore){
    //template<typename T>Array2DSetIneRow(T**&X, Array2DSize& size, int IneRowN,std::vector<T>rowParam, int whatFromN=1, int QDefaultValsBefore=0, T*DefaultValParam=NULL, bool ExistingInsteadOfDefault=true, bool ExistingOnlyForGTLminNotIgnore=false, bool LastPossIfNotRectAtPos0NotIgnore=false){
    Array2DSetIneRow(this->data, this->size, IneRowN, rowParam, whatFromN, QDefaultValsBefore, DefaultValParam, ExistingInsteadOfDefault, ExistingOnlyForGTLminNotIgnore, LastPossIfNotRectAtPos0NotIgnore);
}

template<typename T> void Array2D_PS<T>::AddExtRow(T*rowParam, int wholeRowL, int whatFromN, int QDefaultValsBefore, bool RectNotVar, T*DfltValParam, TValsShowHide*vsh){
    //template<typename T>void Array2DAddExtRow(T**&X, Array2DSize&size, T*rowParam, int wholeRowL=0, int whatFromN=1, int QDefaultValsBefore=0, bool RectNotVar=true, T*DfltValParam=NULL, TValsShowHide*vsh=NULL){
    Array2DAddExtRow(this->data, this->size, rowParam, wholeRowL, whatFromN, QDefaultValsBefore, RectNotVar, DfltValParam, vsh);
}
//
template<typename T> void Array2D_PS<T>::AddExtRow(std::vector<T>rowParam, int whatFromN, int QDefaultValsBefore, bool RectNotVar, T*DfltValParam, TValsShowHide*vsh){
    //template<typename T>void Array2DAddExtRow(T**&X, Array2DSize&size, std::vector<T>rowParam, int whatFromN=1, int QDefaultValsBefore=0, bool RectNotVar=true, T*DfltValParam=NULL, TValsShowHide*vsh=NULL){
    Array2DAddExtRow(this->data, this->size, rowParam, whatFromN, QDefaultValsBefore, RectNotVar, DfltValParam, vsh);
}

template<typename T> void Array2D_PS<T>::InsExtRow(int ExtRowN, int wholeRowL, T*rowParam, T*DefaultValParam, int whatFromN, int QDefaultValuesBefore, bool KeepIfRect, TValsShowHide*vsh){
    //template<typename T> void Array2DInsExtRow(T**&X,  Array2DSize&size, int ExtRowN, int wholeRowL=0, T*rowParam=NULL, T*DefaultValParam=NULL, int whatFromN=1, int QDefaultValuesBefore=0, bool KeepIfRect=true, TValsShowHide*vsh=NULL){
    Array2DInsExtRow(this->data, this->size, ExtRowN, wholeRowL, rowParam, DefaultValParam, whatFromN, QDefaultValuesBefore, KeepIfRect, vsh);
}
template<typename T> void Array2D_PS<T>::InsExtRow(int ExtRowN, std::vector<T>rowParam, T*DefaultValParam, int whatFromN, int QDefaultValuesBefore, bool KeepIfRect){
    //template<typename T> void Array2DInsExtRow(T**&X,  Array2DSize&size, int ExtRowN, std::vector<T>rowParam, T*DefaultValParam=NULL, int whatFromN=1, int QDefaultValuesBefore=0, bool KeepIfRect=true){
    Array2DInsExtRow(this->data, this->size, ExtRowN, rowParam, DefaultValParam, whatFromN, QDefaultValuesBefore, KeepIfRect);
}

template<typename T> void Array2D_PS<T>::AddIneRow(T*rowParam, int wholeRowL, int IfNonRectIgnore0Add1Stretch2, int whatFromN, int QDefaultValsBefore, T*DefaultValParam, bool RectNotVariaIfFirst){
    //template<typename T>Array2DAddIneRow(T**&X, Array2DSize& size, T*rowParam=NULL, int wholeRowL=0, int IfNonRectIgnore0Add1Stretch2=0, int whatFromN=1, int QDefaultValsBefore=0, T*DefaultValParam=NULL, bool RectNotVariaIfFirst=true){
    Array2DAddIneRow(this->data, this->size, rowParam, wholeRowL, IfNonRectIgnore0Add1Stretch2, whatFromN, QDefaultValsBefore, DefaultValParam, RectNotVariaIfFirst);
}
template<typename T> void Array2D_PS<T>::AddIneRow(std::vector<T>rowParam, int IfNonRectIgnore0Add1Stretch2, int whatFromN, int QDefaultValsBefore, T*DefaultValParam){
    //template<typename T>Array2DAddIneRow(T**&X, Array2DSize& size, std::vector<T>rowParam, int IfNonRectIgnore0Add1Stretch2=0, int whatFromN=1, int QDefaultValsBefore=0, T*DefaultValParam=NULL){
    Array2DAddIneRow(this->data, this->size, rowParam, IfNonRectIgnore0Add1Stretch2, whatFromN, QDefaultValsBefore, DefaultValParam);
}//hin nil om rect

template<typename T> void Array2D_PS<T>::InsIneRow(int IneRowPosN, T*rowParam, int wholeRowL, int IfAfterLminIgnore0Stretch2, int whatFromN, int QDefaultValsBefore, T*DefaultValParam){
    //template<typename T>Array2DInsIneRow(T**&X, Array2DSize& size, int IneRowPosN, T*rowParam=NULL, int wholeRowL=0, int IfAfterLminIgnore0Stretch2=0, int whatFromN=1, int QDefaultValsBefore=0, T*DefaultValParam=NULL){
    Array2DInsIneRow(this->data, this->size, IneRowPosN, rowParam, wholeRowL, IfAfterLminIgnore0Stretch2, whatFromN, QDefaultValsBefore, DefaultValParam);
}
template<typename T> void Array2D_PS<T>::InsIneRow(int IneRowPosN, std::vector<T>rowParam, int IfAfterLminIgnore0Stretch2, int whatFromN, int QDefaultValsBefore, T*DefaultValParam){
    //template<typename T>Array2DInsIneRow(T**&X, Array2DSize& size, int IneRowPosN, std::vector<T>rowParam, int IfAfterLminIgnore0Stretch2=0, int whatFromN=1, int QDefaultValsBefore=0, T*DefaultValParam=NULL){
    Array2DInsIneRow(this->data, this->size, IneRowPosN, rowParam, IfAfterLminIgnore0Stretch2, whatFromN, QDefaultValsBefore, DefaultValParam);
}
//
template<typename T> void Array2D_PS<T>::AddToExtRow(const T& val, int ExtRowN){
    //template<typename T>Array2DAddToExtRowN(T**&X, Array2DSize& size, int ExtRowN, T val){
    Array2DAddToExtRowN(this->data, this->size, ExtRowN, val);
}
template<typename T> void Array2D_PS<T>::InsToExtRow(const T& val, int ExtRowN, int PosN){
    //template<typename T>Array2DInsToExtRowN(T**&X, Array2DSize& size, int ExtRowN, T val, int IneRowPosN){
    Array2DInsToExtRowN(this->data, this->size, ExtRowN, val, PosN);
}
template<typename T> void Array2D_PS<T>::DelFromExtRow(int ExtRowN, int PosN){
    //template<typename T>Array2DDelFromExtRowN(T**&X, Array2DSize& size, int ExtRowN, int IneRowPosN){
    Array2DDelFromExtRowN(this->data, this->size, ExtRowN, PosN);
}
//
template<typename T> std::vector<std::vector<int>> Array2D_PS<T>::SeekVal(const T& Val, int paramN){
    return Array2D_SeekVal_Simply(this->data, this->size, Val, paramN);
}

template<typename T> std::vector<std::vector<int>> Array2D_PS<T>::SeekExtRow(T*row, int L, int paramN){
    return Array2D_SeekExtRow_Simply(this->data, this->size,row, L, paramN);
}
template<typename T> std::vector<std::vector<int>> Array2D_PS<T>::SeekExtRow(std::vector<T>row, int paramN){
    return this->SeekExtRow(row.data(), row.size(), paramN);
}

//std::vector<std::vector<int>> SeekIneRow(T*row, int L){
//
//}
//std::vector<std::vector<int>> SeekIneRow(std::vector<T>row);
template<typename T> std::vector<std::vector<int>> Array2D_PS<T>::SeekSubArray(T**What, const Array2DSize& sizeWhat, bool IsTransposed, bool EachExtRowIsReversed){
    //template <typename T> std::vector<std::vector<int>> Array2D_SeekArrs2ndIn1st(T**Where, const Array2DSize& sizeWhere, T**What, const Array2DSize& sizeWhat, bool IsTransposed=false, bool EachExtRowIsReversed=false){
    return Array2D_SeekArrs2ndIn1st(this->data, this->size, What, sizeWhat, IsTransposed, EachExtRowIsReversed);
}
//std::vector<std::vector<int>> SeekSubArray(const std::vector<std::vector<T>>&data, bool Transposed, bool ExtReversed, bool IneReversed);

template<typename T> std::vector<std::vector<int>> Array2D_PS<T>::SeekSubArray(const Array2D_PS& data, bool Transposed, bool EachExtRowIsReversed){
    return this->SeekSubArray(data.data, data.size, Transposed, EachExtRowIsReversed);
}



template<typename T> Array2D_PS<T>& Array2D_PS<T>::GetSubArray(){
    return*this;//do!
}

//
template<typename T> void Array2D_PS<T>::ConcatAdding(){}

template<typename T> void Array2D_PS<T>::ConcatStretching(){}


template<typename T> void Array2D_PS<T>::Arr2DStdCOut2D(QString delimElem, QString delimLines, bool useBrackets){
    //template <typename T>void Arr2DStdCOut2D(T**X, const Array2DSize size, QString delimElem="; ", QString delimLines="; ", bool useBrackets=true){
    Arr2DStdCOut2D(this->data, this->size, delimElem, delimLines, useBrackets);
}

template<typename T> void Array2D_PS<T>::Arr2DStdCOut1D(QString delimElem, QString delimLines, bool useBrackets){
    //template <typename T>void Arr2DStdCOut1D(T**X, const Array2DSize size, QString delimElem=", ", QString delimLines="; ", bool useBrackets=true){
    Arr2DStdCOut1D(this->data, this->size, delimElem, delimLines, useBrackets);
}

//template<typename T> QString ToString(QString rowDelim="; ", QString elemDelim=" ");

/*

template<typename T> class Array2DPB{
    Array2DSize*size;
    T**data;
public:
    //Array2DPB(int QExtRows=1, int IneRowsLength=1, bool VarNotRect=true, T*DefaultValParam=NULL);
    Array2DPB(int QExtRows=1, int IneRowsLength=1, bool RectNotVarNot=true, T*DefaultValParam=NULL);
    Array2DPB(Array2DSize*size, T**X=NULL, T*DefaultValparam=NULL);
    Array2DPB(int Q, bool QRowsNotLengthOfOneExt=true, T*X=NULL, T*DefaultValParam=NULL);
    Array2DPB(std::vector<std::vector<T>>X);
    Array2DPB(const Array2DPB&X);
    ~Array2DPB();
    //
    void SetSize(int QExtRows=1, int IneRowsLength=1, bool VarNotRect=true);
    void SetSize(const Array2DSize&size);
    //
    void Assign(const Array2DPB&X);
    Array2DPB& operator = (const Array2DPB&X);
    //
    T& operator[](int index);
    //
    //void Set(int QExtRows=1, int IneRowsLength=1, bool VarNotRect=true, T*DefaultValParam=NULL);
    void Set(int QExtRows=1, int IneRowsLength=1, bool VarNotRect=true, T*DefaultValParam=NULL);
    void Set(Array2DSize*size, T**X=NULL, T*DefaultValParam=NULL);
    void Set(int Q, bool QRowsNotLengthOfOneExt=true, T*X=NULL, T*DefaultValParam=NULL);
    void Set(std::vector<std::vector<T>>X);
    //
    std::vector<std::vector<T>> Get();
    std::vector<T> GetExtRowN_asVect(int ExtRowN);
    std::vector<T> GetIneRowN_asVect(int IneRowN);
    T*GetExtRowN_asPtr(int ExtRowN);
    T*GetIneRowN_asPtr(int InetRowN);
    //
    Array2DSize GetSize();
    int GetQExtRows();
    int GetLength(int ineRowN);
    //
    void SetVal(T val, int ExtRowN, int IneRowN=1);
    //
    T GetElement_AsVal(int ExtRowN, int IneRowN=1);
    T*GetElement_AsPtr(int ExtRowN, int IneRowN=1);
    //

};

*/

//------------------------------------------------------------------------------------

template<typename T> class Array2D_V{
   std::vector<std::vector<T>>data;
public:
    //Array2DPS(int QExtRows=1, int IneRowsLength=1, bool VarNotRect=true, T*DefaultValParam=NULL);
    Array2D_V(const Array2DSize&size, T**X=NULL, T*DefaultValParam=NULL);
    Array2D_V(std::vector<std::vector<T>>X);
    Array2D_V(int QExtRows=1, int IneRowsLength=1, T**X=NULL, T*DefaultValParam=NULL);
    Array2D_V(int Q, bool QRowsNotLengthOfOneExt=true, T*X=NULL, T*DefaultValParam=NULL);
    Array2D_V(std::vector<T>X, bool QRowsNotLengthOfOneExt=true);
    Array2D_V(const Array2D_V&obj);
    ~Array2D_V();
    //
    void Assign(const Array2D_V&X);
    Array2D_V& operator = (const Array2D_V&obj);
    void SetNull();
    //
    const std::vector<T> operator[](int index) const;
    //
    //void Set(int QExtRows=1, int IneRowsLength=1, bool VarNotRect=true, T*DefaultValParam=NULL);
    void Set(const Array2DSize&size, T**X=NULL, T*DefaultValParam=NULL);
    void Set(std::vector<std::vector<T>>X);
    void Set(int QExtRows=1, int IneRowsLength=1, T**X=NULL, T*DefaultValParam=NULL);
    void Set(int Q, bool QRowsNotLengthOfOneExt=true, T*X=NULL, T*DefaultValParam=NULL);
    void Set(std::vector<T>X, bool QRowsNotLengthOfOneExt=true);
    //
    void SetSize(const Array2DSize&newSize, T*DefaultValParam=NULL);
    void SetSize(int QExtRows=1, int IneRowsLength=1, T*DefaultValParam=NULL, bool RectNotVar=true);
    //
    std::vector<std::vector<T>> Get();
    std::vector<T> GetExtRowN_asVect(int ExtRowN);
    std::vector<T> GetIneRowN_asVect(int IneRowN);
    T*GetExtRowN_asPtr(int ExtRowN);
    //
    Array2DSize GetSize() const;
    int GetQExtRows() const;
    int GetLength(int ineRowN) const;
    int GetMinLength() const;
    int GetMaxLength() const;

    //
    void SetElement(const T& val, int ExtRowN, int IneRowN=1);
    //
    T GetElement_AsVal(int ExtRowN, int IneRowN=1) const;
    T*GetElement_AsPtr(int ExtRowN, int IneRowN=1);
    //
    void SwapVals(int ExtRowN1, int IneRowN1, int ExtRowN2, int IneRowN2);
    void SwapExtRows(int ExtRowN1, int ExtRowN2);
    void SwapIneRows(int IneRowN1, int IneRowN2);
    void ReverseExtRows();
    void ReverseIneRows();
    void ReverseExtRowN(int ExtRowN);
    void ReverseIneRowN(int IneRowN);
    void Transpose();
    //void TransposeStruct();
    //
    void DelExtRowN(int ExtRowN);
    void DelIneRowN(int IneRowN, bool ifPos0IgnoreNotDelLastOfVaria=true);
    //
    void SetExtRow(int ExtRowN, std::vector<T>rowParam, int whatFromN=1, int QDefaultValuesBefore=0, bool KeepIfRect=true, T*DefaultValParam=NULL, bool DefaultIsSpecNotOwn=false);
    void SetExtRow(int ExtRowN, int wholeRowL=0, T*rowParam=NULL, int FromN=1, int QDefaultValuesBefore=0, bool KeepIfRect=true, T*DefaultValParam=NULL, bool DefaultIsSpecNotOwn=false);

    void SetIneRow(int IneRowN, std::vector<T>rowParam, int whatFromN=1, int QDefaultValsBefore=0, T*DefaultValParam=NULL, bool ExistingInsteadOfDefault=true, bool ExistingOnlyForGTLminNotIgnore=false, bool LastPossIfNotRectAtPos0NotIgnore=false);
    void SetIneRow(int IneRowN, T*rowParam, int wholeRowL=0,int whatFromN=1, int QDefaultValsBefore=0, T*DefaultValParam=NULL, bool ExistingInsteadOfDefault=true, bool ExistingOnlyForGTLminNotIgnore=false, bool LastPossIfNotRectAtPos0NotIgnore=false);

    void AddExtRow(std::vector<T>rowParam,  int whatFromN=1, int QDefaultValsBefore=0, bool RectNotVar=true, T*DfltValParam=NULL, TValsShowHide*vsh=NULL);
    void AddExtRow(int wholeRowL=0, T*rowParam=NULL,  int FromN=1, int QDefaultValsBefore=0, bool RectNotVar=true, T*DfltValParam=NULL, TValsShowHide*vsh=NULL);

    void InsExtRow(int ExtRowN, std::vector<T>rowParam, T*DefaultValParam=NULL, int whatFromN=1, int QDefaultValuesBefore=0, bool KeepIfRect=true);
    void InsExtRow(int ExtRowN, int wholeRowL=0, T*rowParam=NULL, T*DefaultValParam=NULL, int whatFromN=1, int QDefaultValuesBefore=0, bool KeepIfRect=true);

    void AddIneRow(std::vector<T>rowParam, int IfNonRectIgnore0Add1Stretch2=0, int whatFromN=1, int QDefaultValsBefore=0, T*DefaultValParam=NULL);
    void AddIneRow(T*rowParam=NULL, int wholeRowL=0, int IfNonRectIgnore0Add1Stretch2=0, int whatFromN=1, int QDefaultValsBefore=0, T*DefaultValParam=NULL);

    void InsIneRow(int IneRowPosN, std::vector<T>rowParam, int IfAfterLminIgnore0Stretch2=0, int whatFromN=1, int QDefaultValsBefore=0, T*DefaultValParam=NULL);
    void InsIneRow(int IneRowPosN, T*rowParam=NULL, int wholeRowL=0,  int IfAfterLminIgnore0Stretch2=0, int whatFromN=1, int QDefaultValsBefore=0, T*DefaultValParam=NULL);
    //
    void AddToExtRow(T val, int ExtRowN);
    void InsToExtRow(T val, int ExtRowN, int IneRowPosN);
    void DelFromExtRow(int ExtRowN, int IneRowPosN);
    //
    std::vector<std::vector<int>> SeekVal(const T& Val, int paramN=0);

    //std::vector<std::vector<int>> SeekExtRow(T*row, int L, int paramN=0){
        //return Array2D_SeekExtRow_Simply(this->data, this->size,
    //}

    //std::vector<std::vector<int>> SeekExtRow(std::vector<T>row, int paramN=0){
        //return Array2D_SeekExtRow_Simply(this->data, this->size, row, paramN);
    //}

    //std::vector<std::vector<int>> SeekIneRow(T*row, int L);
    //std::vector<std::vector<int>> SeekIneRow(std::vector<T>row);
    std::vector<std::vector<int>> SeekSubArray(std::vector<std::vector<T>>What, bool IsTransposed=false, bool EachExtRowIsReversed=false);
    std::vector<std::vector<int>> SeekSubArray(const Array2D_V& data, bool IsTransposed=false, bool EachExtRowIsReversed=false);

    //std::vector<std::vector<int>> SeekSubArray(T**data, const Array2DSize& size, bool Transposed, bool ExtReversed, bool IneReversed);

    Array2D_V GetSubArray(int ExtRowN1, int IneRowN1, int ExtRowN2, int IneRowN2);
    Array2D_V GetSubArray(const Array2DSize& size);
    //
    void ConcatAdding();

    void ConcatStretching();

    void Arr2DStdCOut2D(QString delimElem="; ", QString delimLines="; ", bool useBrackets=true);
    void Arr2DStdCOut1D(QString delimElem=", ", QString delimLines="; ", bool useBrackets=true);
};

template<typename T> Array2D_V<T>::Array2D_V(const Array2DSize&size, T**X, T*DefaultValParam){
    this->Set(size, X, DefaultValParam);
}
template<typename T> Array2D_V<T>::Array2D_V(std::vector<std::vector<T>>X){
    this->Set(X);
}
template<typename T> Array2D_V<T>::Array2D_V(int QExtRows, int IneRowsLength, T**X, T*DefaultValParam){
    this->Set(QExtRows, IneRowsLength, X, DefaultValParam);
}
template<typename T> Array2D_V<T>::Array2D_V(int Q, bool QRowsNotLengthOfOneExt, T*X, T*DefaultValParamr){
    this->Set(Q, QRowsNotLengthOfOneExt, X);
}
template<typename T> Array2D_V<T>::Array2D_V(std::vector<T>X, bool QRowsNotLengthOfOneExt){
    this->Set(QRowsNotLengthOfOneExt, QRowsNotLengthOfOneExt, X);
}

template<typename T> Array2D_V<T>::Array2D_V(const Array2D_V&obj){
    this->Assign(obj);
}

template<typename T> Array2D_V<T>::~Array2D_V(){}



template<typename T> void Array2D_V<T>::Assign(const Array2D_V&obj){
    this->SetNull();
    this->Set(obj.data);
}

template<typename T> Array2D_V<T>& Array2D_V<T>::operator = (const Array2D_V&obj){
    this->Assign(obj);
    return *this;
}

template<typename T> void Array2D_V<T>::SetNull(){
    this->data.clear();
}



template<typename T> const std::vector<T> Array2D_V<T>::operator[](int index) const{
    std::vector<T>row;
    int L;
    if(this->data.size()>0 && index>=0 && index<this->data.size() && this->data[index].size()>0){
        L=this->data[index].size();
        for(int i=1; i<=L; i++){
            row.push_back(this->data[index][i-1]);
        }
    }
    return row;
}

template<typename T> void Array2D_V<T>::Set(const Array2DSize&size, T**X, T*DefaultValParam){
    int Q=size.GetQExtRows(), L;
    std::vector<T>row;
    T val;
    //
    this->SetNull();
    //
    if(X!=NULL){
        for(int i=1; i<=Q; i++){
            row.clear();
            L=size.GetLength(i);
            if(L>0){
                for(int j=1; j<=L; j++){
                    val=X[i-1][j-1];
                    row.push_back(val);
                }
            }
            this->data.push_back(row);
        }
    }else{
        if(DefaultValParam!=NULL){
            for(int i=1; i<=Q; i++){
                row.clear();
                L=size.GetLength(i);
                if(L>0){
                    for(int j=1; j<=L; j++){
                        this->data[i-1][j-1]=(*(DefaultValParam));
                    }
                }
            }
        }
    }
}
template<typename T> void Array2D_V<T>::Set(std::vector<std::vector<T>>X){
   this->data=X;
}

template<typename T> void Array2D_V<T>::Set(int QExtRows, int IneRowsLength, T**X, T*DefaultValParam){
    this->SetNull();
    Array2DSize size;
    //if(RectNotVar){
    //    size.Set(QExtRows, IneRowsLength);
    //}else{
    //    size.SetVaria(QExtRows, IneRowsLength);
    //}
    size.Set(QExtRows, IneRowsLength);
    this->Set(size, X, DefaultValParam);
}
template<typename T> void Array2D_V<T>::Set(int Q, bool QRowsNotLengthOfOneExt, T*X, T*DefaultValParam){
    this->SetNull();
    if(QRowsNotLengthOfOneExt){
        //this->Set(Q, 1, X, DefaultValParam);
        Array2DAddIneRow(this->data, X, Q, 0, 1, 0, DefaultValParam);
    }else{
        //this->Set(1, Q, X, DefaultValParam);
        Array2DAddExtRow(this->data, Q, X, 1, 0, false, DefaultValParam);
    }
}

template<typename T> void Array2D_V<T>::Set(std::vector<T>X, bool QRowsNotLengthOfOneExt){
    T*Dflt=NULL;
    this->SetNull();
    if(QRowsNotLengthOfOneExt){
        Array2DAddIneRow(this->data, X, 0, 1, 0, Dflt);
    }else{
        Array2DAddExtRow(this->data, X, 1, 0, false, Dflt);
    }
}


template<typename T>void  Array2D_V<T>::SetSize(const Array2DSize&newSize, T*DefaultValParam){
    Array2DSetSize(this->data, this->size, newSize, DefaultValParam);
}

template<typename T>void Array2D_V<T>::SetSize(int QExtRows, int IneRowsLength, T*DefaultValParam, bool RectNotVar){
    Array2DSize newSize, oldSize=GetSizeOf2DVector(this->data);
    if(RectNotVar){
        newSize.SetVaria(QExtRows, IneRowsLength);
    }else{
        newSize.SetVaria(QExtRows, IneRowsLength);
    }
    //template<typename T> void Array2DSetSize(std::vector<std::vector<T>>&X, const Array2DSize&size, T*DfltValParam=NULL){//16
    Array2DSetSize(this->data,newSize, DefaultValParam);
}

template<typename T>std::vector<std::vector<T>> Array2D_V<T>:: Get(){
    return this->data;
}

template<typename T>std::vector<T> Array2D_V<T>::GetExtRowN_asVect(int ExtRowN){
    std::vector<T> row;
    T val;
    int Q=this->size.GetQExtRows(), L;
    if(ExtRowN>=1 && ExtRowN<=Q){
        row=this->data[ExtRowN-1];
    }
    return row;
}

template<typename T>std::vector<T>  Array2D_V<T>::GetIneRowN_asVect(int IneRowN){
    std::vector<T> row;
    Array2DSize size=GetSizeOf2DVector(this->data);
    T val;
    int Q=this->size.GetQExtRows(), L, Lmin=size.GetMinLength(), Lmax=size.GetMaxLength();
    if(Q>0 && IneRowN>=1){
        if(IneRowN<=Lmin){
            for(int i=1; i<=Q; i++){
                val=this->data[i-1][IneRowN-1];
                row.push_back(val);
            }
        }else if(IneRowN<=Lmax){

        }
    }
    return row;
}
template<typename T>T* Array2D_V<T>::GetExtRowN_asPtr(int ExtRowN){
    T*R=NULL;
    int Q=this->size.GetQExtRows();
    if(ExtRowN<0){
        ExtRowN=Q+ExtRowN+1;
    }
    if(ExtRowN>=1 && ExtRowN<=Q){
        R=this->data[ExtRowN-1].data();
    }
    return R;
}
//
template<typename T>Array2DSize Array2D_V<T>::GetSize() const{
    return this->size;
}

template<typename T>int  Array2D_V<T>::GetQExtRows() const{
    Array2DSize size=GetSizeOf2DVector(this->data);
    return size.GetQExtRows();
}

template<typename T>int  Array2D_V<T>::GetLength(int ineRowN) const{
    Array2DSize size=GetSizeOf2DVector(this->data);
    return size.GetLength(ineRowN);
}
template<typename T>int Array2D_V<T>::GetMinLength() const{
    Array2DSize size=GetSizeOf2DVector(this->data);
    return size.GetMinLength();
}
template<typename T>int Array2D_V<T>::GetMaxLength() const{
    Array2DSize size=GetSizeOf2DVector(this->data);
    return size.GetMaxLength();
}
//

template<typename T>T  Array2D_V<T>::GetElement_AsVal(int ExtRowN, int IneRowN) const{
    T val;
    Array2DSize size=GetSizeOf2DVector(this->data);
    int Q=size.GetQExtRows(), L=size.GetLength(ExtRowN);
    if(ExtRowN>=1 && ExtRowN<=Q && IneRowN>=1 && IneRowN<=L){
        val=this->data[ExtRowN-1][IneRowN-1];
    }
    return val;
}
template<typename T> T* Array2D_V<T>::GetElement_AsPtr(int ExtRowN, int IneRowN){
    T*ptr;
    int Q=this->size.GetQExtRows(), L=this->size.GetLength(ExtRowN);
    if(ExtRowN>=1 && ExtRowN<=Q && IneRowN>=1 && IneRowN<=L){
        ptr=this->data[ExtRowN-1][IneRowN-1];
    }
    return ptr;
}

template<typename T> void Array2D_V<T>::SetElement(const T& val, int ExtRowN, int IneRowN){
    int Q=this->GetQExtRows(), L=this->GetLength(ExtRowN);
    if(ExtRowN<0){
        ExtRowN=ExtRowN+Q+1;
    }
    if(IneRowN<0){
        IneRowN=IneRowN+L+1;
    }
    if(ExtRowN>=1 && ExtRowN<=Q && IneRowN>=1 && IneRowN<=L){
        this->data[ExtRowN-1][IneRowN-1]=val;
    }
}

template<typename T>void Array2D_V<T>::SwapVals(int ExtRowN1, int IneRowN1, int ExtRowN2, int IneRowN2){
    Array2DSwapVals(this->data, ExtRowN1, IneRowN1, ExtRowN2, IneRowN2);
}

template<typename T>void Array2D_V<T>::SwapExtRows(int ExtRowN1, int ExtRowN2){
    Array2DSwapExtRows(this->data, ExtRowN1, ExtRowN2);
}

template<typename T>void Array2D_V<T>::SwapIneRows(int IneRowN1, int IneRowN2){
    Array2DSwapIneRows(this->data, IneRowN1, IneRowN2);
}

template<typename T>void Array2D_V<T>::ReverseExtRows(){
    Array2DReverseExtRows(this->data);
}

template<typename T>void Array2D_V<T>::ReverseIneRows(){
    Array2DReverseIneRows(this->data);
}

template<typename T>void Array2D_V<T>::ReverseExtRowN(int ExtRowN){
    Array2DReverseExtRowN(this->data, ExtRowN);
}

template<typename T>void Array2D_V<T>::ReverseIneRowN(int IneRowN){
    Array2DReverseIneRowN(this->data, IneRowN);
}

template<typename T>void Array2D_V<T>::Transpose(){
    Array2DTranspose(this->data);
}

template<typename T>void Array2D_V<T>::DelExtRowN(int ExtRowN){
    Array2DDelExtRow(this->data, ExtRowN);
}

template<typename T> void Array2D_V<T>::DelIneRowN(int IneRowN, bool ifPos0IgnoreNotDelLastOfVaria){
    Array2DDelIneRow(this->data, IneRowN, ifPos0IgnoreNotDelLastOfVaria);
}

template<typename T> void Array2D_V<T>::AddToExtRow(T val, int ExtRowN){
    //template<typename T>Array2DAddToExtRowN(std::vector<std::vector<T>>&X, int ExtRowN, T val)
    Array2DAddToExtRowN(this->data, ExtRowN, val);
}

template<typename T> void Array2D_V<T>::InsToExtRow(T val, int ExtRowN, int IneRowPosN){
    //template<typename T>Array2DInsToExtRowN(std::vector<std::vector<T>>&X, int ExtRowN, T val, int IneRowPosN){
    Array2DInsToExtRowN(this->data, ExtRowN, val, IneRowPosN);
}

template<typename T> void Array2D_V<T>::DelFromExtRow(int ExtRowN, int IneRowPosN){
    //template<typename T>Array2DDelFromExtRowN(std::vector<std::vector<T>>&X, int ExtRowN, int IneRowPosN){
    Array2DDelFromExtRowN(this->data, ExtRowN, IneRowPosN);
}


template<typename T>  void Array2D_V<T>::SetExtRow(int ExtRowN, std::vector<T>rowParam, int whatFromN, int QDefaultValuesBefore, bool KeepIfRect, T*DefaultValParam, bool DefaultIsSpecNotOwn){//, int ifDiff_ByArr0ByNew1Min2Max3=0){
    //template<typename T> void Array2DSetExtRowN(std::vector<std::vector<T>>&X, int ExtRowN, std::vector<T>rowParam, int whatFromN=1, int QDefaultValuesBefore=0, bool KeepIfRect=true, T*DefaultValParam=NULL, bool DefaultIsSpecNotOwn=false){//, int ifDiff_ByArr0ByNew1Min2Max3=0){
    Array2DSetExtRowN(this->data, ExtRowN, rowParam, whatFromN, QDefaultValuesBefore, KeepIfRect, DefaultValParam, DefaultIsSpecNotOwn);
}
template<typename T>  void Array2D_V<T>::SetExtRow(int ExtRowN, int wholeRowL, T*rowParam, int FromN, int QDefaultValuesBefore, bool KeepIfRect, T*DefaultValParam, bool DefaultIsSpecNotOwn){//, int ifDiff_ByArr0ByNew1Min2Max3=0){
    //template<typename T> void Array2DSetExtRowN(std::vector<std::vector<T>>&X, int ExtRowN, int wholeRowL=0, T*rowParam=NULL, int FromN=1, int QDefaultValuesBefore=0, bool KeepIfRect=true, T*DefaultValParam=NULL, bool DefaultIsSpecNotOwn=false){//, int ifDiff_ByArr0ByNew1Min2Max3=0){
    Array2DSetExtRowN(this->data, ExtRowN, wholeRowL, rowParam, FromN, QDefaultValuesBefore, KeepIfRect, DefaultValParam, DefaultIsSpecNotOwn);
}

template<typename T>  void Array2D_V<T>::SetIneRow(int IneRowN, std::vector<T>rowParam, int whatFromN, int QDefaultValsBefore, T*DefaultValParam, bool ExistingInsteadOfDefault, bool ExistingOnlyForGTLminNotIgnore, bool LastPossIfNotRectAtPos0NotIgnore){
    //template<typename T>Array2DSetIneRow(std::vector<std::vector<T>>&X, int IneRowN, std::vector<T>rowParam, int whatFromN=1, int QDefaultValsBefore=0, T*DefaultValParam=NULL, bool ExistingInsteadOfDefault=true, bool ExistingOnlyForGTLminNotIgnore=false, bool LastPossIfNotRectAtPos0NotIgnore=false){
    Array2DSetIneRow(this->data, IneRowN, rowParam, whatFromN, QDefaultValsBefore, DefaultValParam, ExistingInsteadOfDefault, ExistingOnlyForGTLminNotIgnore, LastPossIfNotRectAtPos0NotIgnore);
}
template<typename T>  void Array2D_V<T>::SetIneRow(int IneRowN, T*rowParam, int wholeRowL,int whatFromN, int QDefaultValsBefore, T*DefaultValParam, bool ExistingInsteadOfDefault, bool ExistingOnlyForGTLminNotIgnore, bool LastPossIfNotRectAtPos0NotIgnore){
    //template<typename T>Array2DSetIneRow(std::vector<std::vector<T>>&X, int IneRowN, T*rowParam, int wholeRowL=0,int whatFromN=1, int QDefaultValsBefore=0, T*DefaultValParam=NULL, bool ExistingInsteadOfDefault=true, bool ExistingOnlyForGTLminNotIgnore=false, bool LastPossIfNotRectAtPos0NotIgnore=false){
     Array2DSetIneRow(this->data, IneRowN, rowParam, wholeRowL,whatFromN, QDefaultValsBefore, DefaultValParam, ExistingInsteadOfDefault, ExistingOnlyForGTLminNotIgnore, LastPossIfNotRectAtPos0NotIgnore);
}

template<typename T> void Array2D_V<T>::AddExtRow(std::vector<T>rowParam,  int whatFromN, int QDefaultValsBefore, bool RectNotVar, T*DfltValParam, TValsShowHide*vsh){
    //template<typename T>void Array2DAddExtRow(std::vector<std::vector<T>>&X, std::vector<T>rowParam,  int whatFromN=1, int QDefaultValsBefore=0, bool RectNotVar=true, T*DfltValParam=NULL, TValsShowHide*vsh=NULL){
    Array2DAddExtRow(this->data, rowParam,  whatFromN, QDefaultValsBefore, RectNotVar, DfltValParam, vsh);
}
template<typename T> void Array2D_V<T>::AddExtRow(int wholeRowL, T*rowParam,  int FromN, int QDefaultValsBefore, bool RectNotVar, T*DfltValParam, TValsShowHide*vsh){
    //template<typename T>void Array2DAddExtRow(std::vector<std::vector<T>>&X, int wholeRowL=0, T*rowParam=NULL,  int FromN=1, int QDefaultValsBefore=0, bool RectNotVar=true, T*DfltValParam=NULL, TValsShowHide*vsh=NULL){
    Array2DAddExtRow(this->data, wholeRowL, rowParam,  FromN, QDefaultValsBefore, RectNotVar, DfltValParam, vsh);
}

template<typename T> void Array2D_V<T>::InsExtRow(int ExtRowN, std::vector<T>rowParam, T*DefaultValParam, int whatFromN, int QDefaultValuesBefore, bool KeepIfRect){
    //template<typename T> void Array2DInsExtRow(std::vector<std::vector<T>>&X, int ExtRowN, std::vector<T>rowParam, T*DefaultValParam=NULL, int whatFromN=1, int QDefaultValuesBefore=0, bool KeepIfRect=true){
    Array2DInsExtRow(this->data, ExtRowN, rowParam, DefaultValParam, whatFromN, QDefaultValuesBefore, KeepIfRect);
}
template<typename T>  void Array2D_V<T>::InsExtRow(int ExtRowN, int wholeRowL, T*rowParam, T*DefaultValParam, int whatFromN, int QDefaultValuesBefore, bool KeepIfRect){
    //template<typename T> void Array2DInsExtRow(std::vector<std::vector<T>>&X, int ExtRowN, int wholeRowL=0, T*rowParam=NULL, T*DefaultValParam=NULL, int whatFromN=1, int QDefaultValuesBefore=0, bool KeepIfRect=true){
    Array2DInsExtRow(this->data, ExtRowN, wholeRowL, rowParam, DefaultValParam, whatFromN, QDefaultValuesBefore, KeepIfRect);
}

template<typename T> void Array2D_V<T>::AddIneRow(std::vector<T>rowParam, int IfNonRectIgnore0Add1Stretch2, int whatFromN, int QDefaultValsBefore, T*DefaultValParam){
    //template<typename T>Array2DAddIneRow(std::vector<std::vector<T>>&X, std::vector<T>rowParam, int IfNonRectIgnore0Add1Stretch2=0, int whatFromN=1, int QDefaultValsBefore=0, T*DefaultValParam=NULL){
    Array2DAddIneRow(this->data, rowParam, IfNonRectIgnore0Add1Stretch2, whatFromN, QDefaultValsBefore, DefaultValParam);
}
template<typename T> void Array2D_V<T>::AddIneRow(T*rowParam, int wholeRowL, int IfNonRectIgnore0Add1Stretch2, int whatFromN, int QDefaultValsBefore, T*DefaultValParam){
    //template<typename T>Array2DAddIneRow(std::vector<std::vector<T>>&X, T*rowParam=NULL, int wholeRowL=0, int IfNonRectIgnore0Add1Stretch2=0, int whatFromN=1, int QDefaultValsBefore=0, T*DefaultValParam=NULL){
    Array2DAddIneRow(this->data, rowParam, wholeRowL, IfNonRectIgnore0Add1Stretch2, whatFromN, QDefaultValsBefore, DefaultValParam);
}

template<typename T> void Array2D_V<T>::InsIneRow(int IneRowPosN, std::vector<T>rowParam, int IfAfterLminIgnore0Stretch2, int whatFromN, int QDefaultValsBefore, T*DefaultValParam){
    //template<typename T>Array2DInsIneRow(std::vector<std::vector<T>>&X, int IneRowPosN, std::vector<T>rowParam, int IfAfterLminIgnore0Stretch2=0, int whatFromN=1, int QDefaultValsBefore=0, T*DefaultValParam=NULL){
    Array2DInsIneRow(this->data, IneRowPosN, rowParam, IfAfterLminIgnore0Stretch2, whatFromN, QDefaultValsBefore, DefaultValParam);
}
template<typename T> void Array2D_V<T>::InsIneRow(int IneRowPosN, T*rowParam, int wholeRowL,  int IfAfterLminIgnore0Stretch2, int whatFromN, int QDefaultValsBefore, T*DefaultValParam){
    //template<typename T>Array2DInsIneRow(std::vector<std::vector<T>>&X, int IneRowPosN, T*rowParam=NULL, int wholeRowL=0,  int IfAfterLminIgnore0Stretch2=0, int whatFromN=1, int QDefaultValsBefore=0, T*DefaultValParam=NULL){
     Array2DInsIneRow(this->data, IneRowPosN, rowParam, wholeRowL,  IfAfterLminIgnore0Stretch2, whatFromN, QDefaultValsBefore, DefaultValParam);
}

template<typename T> std::vector<std::vector<int>> Array2D_V<T>::SeekVal(const T& Val, int paramN){
    return Array2D_SeekVal_Simply(this->data, this->size, Val, paramN);
}

//template<typename T> std::vector<std::vector<int>> SeekExtRow(T*row, int L, int paramN=0){
      //return Array2D_SeekExtRow_Simply(this->data, this->size,
//}

//template<typename T> std::vector<std::vector<int>> SeekExtRow(std::vector<T>row, int paramN=0){
     //return Array2D_SeekExtRow_Simply(this->data, this->size, row, paramN);
//}

//std::vector<std::vector<int>> SeekIneRow(T*row, int L);
 //std::vector<std::vector<int>> SeekIneRow(std::vector<T>row);
template<typename T> std::vector<std::vector<int>> Array2D_V<T>::SeekSubArray(std::vector<std::vector<T>>What, bool IsTransposed, bool EachExtRowIsReversed ){
    //template <typename T> std::vector<std::vector<int>> Array2D_SeekArrs2ndIn1st(std::vector<std::vector<T>>Where,  std::vector<std::vector<T>>What, bool IsTransposed=false, bool EachExtRowIsReversed=false){
    return Array2D_SeekArrs2ndIn1st(this->data,  What, IsTransposed, EachExtRowIsReversed);
}

//template<typename T> std::vector<std::vector<int>> Array2D_V<T>::SeekSubArray(T**data, const Array2DSize& size, bool Transposed, bool ExtReversed, bool IneReversed);

template<typename T> Array2D_V<T> Array2D_V<T>::GetSubArray(int ExtRowM1, int IneRowN1, int ExtRowN2, int IneRowN2){
    Array2D_V R;
    return R;
}
template<typename T> Array2D_V<T> Array2D_V<T>::GetSubArray(const Array2DSize& size){
    Array2D_V R;
    return R;
}

    //
template<typename T> void Array2D_V<T>::ConcatAdding(){}

template<typename T> void Array2D_V<T>::ConcatStretching(){}

template<typename T> void Array2D_V<T>::Arr2DStdCOut2D(QString delimElem, QString delimLines, bool useBrackets){
    //template <typename T>void Arr2DStdCOut2D(std::vector<std::vector<T>>X, QString delimElem="; ", QString delimLines="; ", bool useBrackets=true){
    Arr2DStdCOut2D(this->data, delimElem, delimLines, useBrackets);
}

template<typename T> void Array2D_V<T>::Arr2DStdCOut1D(QString delimElem, QString delimLines, bool useBrackets){
    //template <typename T>void Arr2DStdCOut1D(std::vector<std::vector<T>>X,  QString delimElem=", ", QString delimLines="; ", bool useBrackets=true){
    Arr2DStdCOut1D(this->data, delimElem, delimLines, useBrackets);
}

//====================================================================================================================================================



//=====================================================================================================================================================

template<typename T> class Array2DRect_V1D{
   std::vector<T>data;
   int ExtRowLength;
public:
    //Array2DPS(int QExtRows=1, int IneRowsLength=1, bool VarNotRect=true, T*DefaultValParam=NULL);
    Array2DRect_V1D(const Array2DSize&size, T**X=NULL, T*DefaultValParam=NULL);
    Array2DRect_V1D(std::vector<std::vector<T>>X);
    Array2DRect_V1D(int QExtRows=1, int IneRowsLength=1, T**X=NULL, T*DefaultValParam=NULL);
    Array2DRect_V1D(int Q, bool QRowsNotLengthOfOneExt=true, T*X=NULL, T*DefaultValParam=NULL);
    Array2DRect_V1D(std::vector<T>X, bool QRowsNotLengthOfOneExt=true);
    Array2DRect_V1D(const Array2DRect_V1D&obj);
    ~Array2DRect_V1D();
    //
    void Assign(const Array2DRect_V1D&X);
    Array2DRect_V1D& operator = (const Array2DRect_V1D&obj);
    void SetNull();
    //
    const std::vector<T> operator[](int index) const;
    //
    //void Set(int QExtRows=1, int IneRowsLength=1, bool VarNotRect=true, T*DefaultValParam=NULL);
    //void Set(const Array2DSize&size, T**X=NULL, T*DefaultValParam=NULL);
    void Set(std::vector<std::vector<T>>X);
    void Set(int QExtRows=1, int IneRowsLength=1, T**X=NULL, T*DefaultValParam=NULL);
    void Set(int Q, bool QRowsNotLengthOfOneExt=true, T*X=NULL, T*DefaultValParam=NULL);
    void Set(std::vector<T>X, bool QRowsNotLengthOfOneExt=true);
    //
    //void SetSize(const Array2DSize&newSize, T*DefaultValParam=NULL);
    void SetSize(int QExtRows=1, int IneRowsLength=1, T*DefaultValParam=NULL, bool RectNotVar=true);
    //
    std::vector<std::vector<T>> Get();
    std::vector<T> GetExtRowN_asVect(int ExtRowN);
    std::vector<T> GetIneRowN_asVect(int IneRowN);
    //T*GetExtRowN_asPtr(int ExtRowN);
    //
    Array2DSize GetSize() const;
    int GetQExtRows() const;
    int GetLength(int ineRowN=1) const;
    int GetMinLength() const;
    int GetMaxLength() const;
    void GetCoordsByPosition(int PosN, int&ExtRowN, int&IneRowN);
    int GetExtRowNByPosition(int PosN);
    int GetIneRowNByPosition(int PosN);
    int GetPositionByCoords(int ExtRowN, int IneRowN);
    //
    void SetElement(const T& val, int ExtRowN, int IneRowN=1);
    //
    T GetElement_AsVal(int ExtRowN, int IneRowN=1) const;
    T*GetElement_AsPtr(int ExtRowN, int IneRowN=1);
    //
    void SwapVals(int ExtRowN1, int IneRowN1, int ExtRowN2, int IneRowN2);
    void SwapExtRows(int ExtRowN1, int ExtRowN2);
    void SwapIneRows(int IneRowN1, int IneRowN2);
    void ReverseExtRows();
    void ReverseIneRows();
    void ReverseExtRowN(int ExtRowN);
    void ReverseIneRowN(int IneRowN);
    void Transpose();
    //void TransposeStruct();
    //
    void DelExtRowN(int ExtRowN);
    void DelIneRowN(int IneRowN, bool ifPos0IgnoreNotDelLastOfVaria=true);
    //
    void SetExtRow(int ExtRowN, std::vector<T>rowParam, int whatFromN=1, int QDefaultValuesBefore=0, bool KeepIfRect=true, T*DefaultValParam=NULL, bool DefaultIsSpecNotOwn=false);
    void SetExtRow(int ExtRowN, int wholeRowL=0, T*rowParam=NULL, int FromN=1, int QDefaultValuesBefore=0, bool KeepIfRect=true, T*DefaultValParam=NULL, bool DefaultIsSpecNotOwn=false);

    void SetIneRow(int IneRowN, std::vector<T>rowParam, int whatFromN=1, int QDefaultValsBefore=0, T*DefaultValParam=NULL, bool ExistingInsteadOfDefault=true, bool ExistingOnlyForGTLminNotIgnore=false, bool LastPossIfNotRectAtPos0NotIgnore=false);
    void SetIneRow(int IneRowN, T*rowParam, int wholeRowL=0,int whatFromN=1, int QDefaultValsBefore=0, T*DefaultValParam=NULL, bool ExistingInsteadOfDefault=true, bool ExistingOnlyForGTLminNotIgnore=false, bool LastPossIfNotRectAtPos0NotIgnore=false);

    void AddExtRow(std::vector<T>rowParam,  int whatFromN=1, int QDefaultValsBefore=0, bool RectNotVar=true, T*DfltValParam=NULL, TValsShowHide*vsh=NULL);
    void AddExtRow(int wholeRowL=0, T*rowParam=NULL,  int FromN=1, int QDefaultValsBefore=0, bool RectNotVar=true, T*DfltValParam=NULL, TValsShowHide*vsh=NULL);

    void InsExtRow(int ExtRowN, std::vector<T>rowParam, T*DefaultValParam=NULL, int whatFromN=1, int QDefaultValuesBefore=0, bool KeepIfRect=true);
    void InsExtRow(int ExtRowN, int wholeRowL=0, T*rowParam=NULL, T*DefaultValParam=NULL, int whatFromN=1, int QDefaultValuesBefore=0, bool KeepIfRect=true);

    void AddIneRow(std::vector<T>rowParam, int IfNonRectIgnore0Add1Stretch2=0, int whatFromN=1, int QDefaultValsBefore=0, T*DefaultValParam=NULL);
    void AddIneRow(T*rowParam=NULL, int wholeRowL=0, int IfNonRectIgnore0Add1Stretch2=0, int whatFromN=1, int QDefaultValsBefore=0, T*DefaultValParam=NULL);

    void InsIneRow(int IneRowPosN, std::vector<T>rowParam, int IfAfterLminIgnore0Stretch2=0, int whatFromN=1, int QDefaultValsBefore=0, T*DefaultValParam=NULL);
    void InsIneRow(int IneRowPosN, T*rowParam=NULL, int wholeRowL=0,  int IfAfterLminIgnore0Stretch2=0, int whatFromN=1, int QDefaultValsBefore=0, T*DefaultValParam=NULL);
    //
    //void AddToExtRow(T val, int ExtRowN);
    //void InsToExtRow(T val, int ExtRowN, int IneRowPosN);
    //void DelFromExtRow(int ExtRowN, int IneRowPosN);
    //
    std::vector<std::vector<int>> SeekVal(const T& Val, int paramN=0);

    //std::vector<std::vector<int>> SeekExtRow(T*row, int L, int paramN=0){
        //return Array2D_SeekExtRow_Simply(this->data, this->size,
    //}

    //std::vector<std::vector<int>> SeekExtRow(std::vector<T>row, int paramN=0){
        //return Array2D_SeekExtRow_Simply(this->data, this->size, row, paramN);
    //}

    //std::vector<std::vector<int>> SeekIneRow(T*row, int L);
    //std::vector<std::vector<int>> SeekIneRow(std::vector<T>row);
    std::vector<std::vector<int>> SeekSubArray(std::vector<std::vector<T>>What, bool IsTransposed=false, bool EachExtRowIsReversed=false);
    std::vector<std::vector<int>> SeekSubArray(const Array2DRect_V1D& data, bool IsTransposed=false, bool EachExtRowIsReversed=false);

    //std::vector<std::vector<int>> SeekSubArray(T**data, const Array2DSize& size, bool Transposed, bool ExtReversed, bool IneReversed);

    Array2DRect_V1D GetSubArray(int ExtRowN1, int IneRowN1, int ExtRowN2, int IneRowN2);
    Array2DRect_V1D GetSubArray(const Array2DSize& size);
    //
    void ConcatAdding();

    void ConcatStretching();

    void Arr2DStdCOut2D(QString delimElem="; ", QString delimLines="; ", bool useBrackets=true);
    void Arr2DStdCOut1D(QString delimElem=", ", QString delimLines="; ", bool useBrackets=true);
};

template<typename T> Array2DRect_V1D<T>::Array2DRect_V1D(const Array2DSize&size, T**X, T*DefaultValParam){
    this->Set(size, X, DefaultValParam);
}
template<typename T> Array2DRect_V1D<T>::Array2DRect_V1D(std::vector<std::vector<T>>X){
    this->Set(X);
}
template<typename T> Array2DRect_V1D<T>::Array2DRect_V1D(int QExtRows, int IneRowsLength, T**X, T*DefaultValParam){
    this->Set(QExtRows, IneRowsLength, X, DefaultValParam);
}
template<typename T> Array2DRect_V1D<T>::Array2DRect_V1D(int Q, bool QRowsNotLengthOfOneExt, T*X, T*DefaultValParamr){
    this->Set(Q, QRowsNotLengthOfOneExt, X);
}
template<typename T> Array2DRect_V1D<T>::Array2DRect_V1D(std::vector<T>X, bool QRowsNotLengthOfOneExt){
    this->Set(QRowsNotLengthOfOneExt, QRowsNotLengthOfOneExt, X);
}

template<typename T> Array2DRect_V1D<T>::Array2DRect_V1D(const Array2DRect_V1D&obj){
    this->Assign(obj);
}

template<typename T> Array2DRect_V1D<T>::~Array2DRect_V1D(){}



template<typename T> void Array2DRect_V1D<T>::Assign(const Array2DRect_V1D&obj){
    this->SetNull();
    this->Set(obj.data);
}

template<typename T> Array2DRect_V1D<T>& Array2DRect_V1D<T>::operator = (const Array2DRect_V1D&obj){
    this->Assign(obj);
    return *this;
}

template<typename T> void Array2DRect_V1D<T>::SetNull(){
    this->data.clear();
}



template<typename T> const std::vector<T> Array2DRect_V1D<T>::operator[](int index) const{
    std::vector<T>row;
    int L;
    if(this->data.size()>0 && index>=0 && index<this->data.size() && this->data[index].size()>0){
        L=this->data[index].size();
        for(int i=1; i<=L; i++){
            row.push_back(this->data[index][i-1]);
        }
    }
    return row;
}

/*template<typename T> void Array2DRect_V1D<T>::Set(const Array2DSize&size, T**X, T*DefaultValParam){
    int Q=size.GetQExtRows(), L;
    std::vector<T>row;
    T val;
    //
    this->SetNull();
    //
    if(X!=NULL){
        for(int i=1; i<=Q; i++){
            row.clear();
            L=size.GetLength(i);
            if(L>0){
                for(int j=1; j<=L; j++){
                    val=X[i-1][j-1];
                    row.push_back(val);
                }
            }
            this->data.push_back(row);
        }
    }else{
        if(DefaultValParam!=NULL){
            for(int i=1; i<=Q; i++){
                row.clear();
                L=size.GetLength(i);
                if(L>0){
                    for(int j=1; j<=L; j++){
                        this->data[i-1][j-1]=(*(DefaultValParam));
                    }
                }
            }
        }
    }
}*/
template<typename T> void Array2DRect_V1D<T>::Set(std::vector<std::vector<T>>X){
    //this->data=X;
    Array2DSize size= GetSizeOf2DVector(X);
    T val;
    int minL=size.GetMinLength(), maxL=size.GetMaxLength(), L, QExtRows=X.size();
    if(minL==maxL){
        this->data.clear();
        L=minL;
        this->ExtRowLength=minL;
        for(int i=1; i<=QExtRows; i++){
            for(int j=1; j<=L; i++){
                val=X[i-1][j-1];
                this->data.push_back(val);
            }
        }
    }
}

template<typename T> void Array2DRect_V1D<T>::Set(int QExtRows, int IneRowsLength, T**X, T*DefaultValParam){
    this->SetNull();
    T val, dfltVal;
    this->data.clear();
    this->ExtRowLength=IneRowsLength;
    if(X!=NULL){
        for(int i=1; i<=QExtRows; i++){
            for(int j=1; j<=this->ExtRowLength; i++){
                val=X[i-1][j-1];
                this->data.push_back(val);
            }
        }
    }else{
        if(DefaultValParam!=NULL){
            dfltVal=(*(DefaultValParam));
        }
        for(int i=1; i<=QExtRows; i++){
            for(int j=1; j<=this->ExtRowLength; i++){
                this->data.push_back(dfltVal);
            }
        }
    }
}
template<typename T> void Array2DRect_V1D<T>::Set(int Q, bool QRowsNotLengthOfOneExt, T*X, T*DefaultValParam){
    this->SetNull();
    std::vector<std::vector<T>>data;
    if(QRowsNotLengthOfOneExt){
        //this->Set(Q, 1, X, DefaultValParam);
        Array2DAddIneRow(data, X, Q, 0, 1, 0, DefaultValParam);
    }else{
        //this->Set(1, Q, X, DefaultValParam);
        Array2DAddExtRow(data, Q, X, 1, 0, false, DefaultValParam);
    }
    this->Set(data);
}

template<typename T> void Array2DRect_V1D<T>::Set(std::vector<T>X, bool QRowsNotLengthOfOneExt){
    T*Dflt=NULL;
    std::vector<std::vector<T>>data;
    this->SetNull();
    if(QRowsNotLengthOfOneExt){
        Array2DAddIneRow(data, X, 0, 1, 0, Dflt);
    }else{
        Array2DAddExtRow(data, X, 1, 0, false, Dflt);
    }
    this->Set(data);
}


/*template<typename T>void  Array2DRect_V1D<T>::SetSize(const Array2DSize&newSize, T*DefaultValParam){
    Array2DSetSize(this->data, this->size, newSize, DefaultValParam);
}*/

template<typename T>void Array2DRect_V1D<T>::SetSize(int QExtRows, int IneRowsLength, T*DefaultValParam, bool RectNotVar){
    T curVal, dfltVal;
    int QExtRowsOld=this->GetQExtRows();
    int ExtRowsLengthOld=this->GetLength();
    int PosN;
    if(DefaultValParam!=NULL){
        dfltVal=(*(DefaultValParam));
    }
    if(IneRowsLength>this->ExtRowLength){
        for(int i=QExtRows; i>=1; i--){
            for(int j=1; j<=(IneRowsLength-this->ExtRowLength); j++){
                PosN=(i-1)*IneRowsLength+1;
                Array1DInsElementToN(this->data, PosN, dfltVal);
            }
        }
    }else if(IneRowsLength<this->ExtRowLength){
        for(int i=QExtRows; i>=1; i--){
            for(int j=1; j<=(IneRowsLength-this->ExtRowLength); j++){
                PosN=(i-1)*IneRowsLength+1;
                Array1DDelElementFromN(this->data, PosN);
            }
        }
    }
    //now all rows are formatted with same lengthes
    if(QExtRows>QExtRowsOld){
        for(int i=1; i<=QExtRows>QExtRowsOld; i++){
            for(int j=1; j<=IneRowsLength; i++){
                this->data.push_back(dfltVal);
            }
        }
    }else if(QExtRows<QExtRowsOld){
        PosN=this->data.size();
        for(int i=1; i<=QExtRows-QExtRowsOld; i++){
            for(int j=1; j<=IneRowsLength; i++){
                PosN=this->data.size();
                Array1DDelElementFromN(this->data, PosN);
            }
        }
    }
}

template<typename T>std::vector<std::vector<T>> Array2DRect_V1D<T>:: Get(){
    std::vector<std::vector<T>>data;
    std::vector<T>row;
    T val;
    int PosN;
    int QExtRows=this->GetQExtRows();
    for(int i=1; i<=QExtRows; i++){
        row.clear();
        for(int j=1; j<=this->ExtRowLength; j++){
            PosN=this->GetPositionByCoords(i, j);
            val=this->data[PosN-1];
            row.push_back(val);
        }
        data.push_back(row);
    }
    return this->data;
}

template<typename T>std::vector<T> Array2DRect_V1D<T>::GetExtRowN_asVect(int ExtRowN){
    std::vector<T> row;
    T val;
    int Q=this->size.GetQExtRows(), L, PosN;
    if(ExtRowN>=1 && ExtRowN<=Q){
        for(int j=1; j<=this->ExtRowLength; j++){
            PosN=this->GetPositionByCoords(ExtRowN, j);
            val=this->data[PosN-1];
            row.push_back(val);
        }
    }
    return row;
}

template<typename T>std::vector<T>  Array2DRect_V1D<T>::GetIneRowN_asVect(int IneRowN){
    std::vector<T> row;
    Array2DSize size=GetSizeOf2DVector(this->data);
    T val;
    int PosN;
    int Q=this->size.GetQExtRows(), L, Lmin=size.GetMinLength(), Lmax=size.GetMaxLength();
    if(Q>0 && IneRowN>=1){
        if(IneRowN<=Lmin){
            for(int i=1; i<=Q; i++){
                PosN=this->GetPositionByCoords(i, IneRowN);
                val=this->data[PosN-1];
                row.push_back(val);
            }
        }else if(IneRowN<=Lmax){

        }
    }
    return row;
}
/*template<typename T>T* Array2DRect_V1D<T>::GetExtRowN_asPtr(int ExtRowN){
    T*R=NULL;
    int Q=this->size.GetQExtRows();
    if(ExtRowN<0){
        ExtRowN=Q+ExtRowN+1;
    }
    if(ExtRowN>=1 && ExtRowN<=Q){
        R=this->data[ExtRowN-1].data();
    }
    return R;
}*/
//
template<typename T>Array2DSize Array2DRect_V1D<T>::GetSize() const{
    return this->size;
}

template<typename T>int  Array2DRect_V1D<T>::GetQExtRows() const{
    int QComponents=this->data.size();
    int QExtRows=QComponents/this->ExtRowLength;
    return QExtRows;
}

template<typename T>int  Array2DRect_V1D<T>::GetLength(int ineRowN) const{
    int QComponents=this->data.size();
    int QExtRows=QComponents/this->ExtRowLength;
    return this->ExtRowLength;
}
template<typename T>int Array2DRect_V1D<T>::GetMinLength() const{
    int QComponents=this->data.size();
    int QExtRows=QComponents/this->ExtRowLength;
    return this->ExtRowLength;
}
template<typename T>int Array2DRect_V1D<T>::GetMaxLength() const{
    int QComponents=this->data.size();
    int QExtRows=QComponents/this->ExtRowLength;
    return this->ExtRowLength;
}
//
template<typename T>void  Array2DRect_V1D<T>::GetCoordsByPosition(int PosN, int&ExtRowN, int&IneRowN){
    int QComponents=this->data.size();
    int QExtRows=QComponents/this->ExtRowLength;
    ExtRowN=PosN/this->ExtRowLength+1;
    IneRowN=PosN%this->ExtRowLength;;
    IneRowN=PosN-(ExtRowN-1)*this->ExtRowLength;
}

template<typename T>int  Array2DRect_V1D<T>::GetExtRowNByPosition(int PosN){
    int ExtRowN, IneRowN, QComponents, QExtRows;
    this->GetCoordsByPosition(PosN, ExtRowN, IneRowN);
    return ExtRowN;
}

template<typename T>int  Array2DRect_V1D<T>::GetIneRowNByPosition(int PosN){
    int ExtRowN, IneRowN, QComponents, QExtRows;
    this->GetCoordsByPosition(PosN, ExtRowN, IneRowN);
    return IneRowN;
}
template<typename T>int  Array2DRect_V1D<T>::GetPositionByCoords(int ExtRowN, int IneRowN){
    int QComponents=this->data.size();
    int QExtRows=QComponents/this->ExtRowLength;
    int PosN=(ExtRowN-1)*this->ExtRowLength+IneRowN;
    return PosN;
}

//
template<typename T>T  Array2DRect_V1D<T>::GetElement_AsVal(int ExtRowN, int IneRowN) const{
    T val;
    int PosN;
    Array2DSize size=GetSizeOf2DVector(this->data);
    int Q=size.GetQExtRows(), L=size.GetLength(ExtRowN);
    if(ExtRowN>=1 && ExtRowN<=Q && IneRowN>=1 && IneRowN<=L){
        PosN=this->GetPositionByCoords(ExtRowN, IneRowN);
        val=this->data[PosN-1];
    }
    return val;
}
template<typename T> T* Array2DRect_V1D<T>::GetElement_AsPtr(int ExtRowN, int IneRowN){
    T*ptr=NULL;
    int Q=this->size.GetQExtRows(), L=this->size.GetLength(ExtRowN), PosN;
    if(ExtRowN>=1 && ExtRowN<=Q && IneRowN>=1 && IneRowN<=L){
        PosN=this->GetPositionByCoords(ExtRowN, IneRowN);
        ptr=this->data[PosN-1];
    }
    return ptr;
}

template<typename T> void Array2DRect_V1D<T>::SetElement(const T& val, int ExtRowN, int IneRowN){
    int Q=this->GetQExtRows(), L=this->GetLength(ExtRowN), PosN;
    if(ExtRowN<0){
        ExtRowN=ExtRowN+Q+1;
    }
    if(IneRowN<0){
        IneRowN=IneRowN+L+1;
    }
    if(ExtRowN>=1 && ExtRowN<=Q && IneRowN>=1 && IneRowN<=L){
        PosN=this->GetPositionByCoords(ExtRowN, IneRowN);
        this->data[PosN-1]=val;
    }
}

template<typename T>void Array2DRect_V1D<T>::SwapVals(int ExtRowN1, int IneRowN1, int ExtRowN2, int IneRowN2){
    int Pos1N=this->GetPositionByCoords(ExtRowN1, IneRowN1);
    int Pos2N=this->GetPositionByCoords(ExtRowN2, IneRowN2);
    Array1DSwapVals(this->data, Pos1N, Pos2N);
}

template<typename T>void Array2DRect_V1D<T>::SwapExtRows(int ExtRowN1, int ExtRowN2){
    int Pos1, Pos2;
    for(int i=1; i<=this->ExtRowLength; i++){
        Pos1=this->GetPositionByCoords(ExtRowN1, i);
        Pos2=this->GetPositionByCoords(ExtRowN2, i);
        Array1DSwapVals(this->data, Pos1, Pos2);
    }
}

template<typename T>void Array2DRect_V1D<T>::SwapIneRows(int IneRowN1, int IneRowN2){
    int Pos1, Pos2, QExtRows=this->GetQExtRows();
    if(IneRowN1>=1 && IneRowN1<=this->ExtRowLength && IneRowN2>=1 && IneRowN2<=this->ExtRowLength){
        for(int i=1; i<=QExtRows; i++){
            Pos1=this->GetPositionByCoords(i, IneRowN1);
            Pos2=this->GetPositionByCoords(i, IneRowN2);
            Array1DSwapVals(this->data, Pos1, Pos2);
        }
    }
}

template<typename T>void Array2DRect_V1D<T>::ReverseExtRows(){
    //Array2DReverseExtRows(this->data);
}

template<typename T>void Array2DRect_V1D<T>::ReverseIneRows(){
    Array2DReverseIneRows(this->data);
}

template<typename T>void Array2DRect_V1D<T>::ReverseExtRowN(int ExtRowN){
    int PosErst=this->GetPositionByCoords(ExtRowN, 1),
        PosLast=this->GetPositionByCoords(ExtRowN, this->ExtRowLength),
        Q=PosLast-PosErst+1, N1, N2, Pos1, Pos2;
    if(Q%2==0){
        N1=Q/2;
    }else{
        N1=(Q-1)/2;
    }
    for(int i=1; i<=N1; i++){
        Pos1=PosErst+i-1;
        Pos2=PosLast+1-i;
        Array1DSwapElements(this->data, Pos1, Pos2);
    }

}

template<typename T>void Array2DRect_V1D<T>::ReverseIneRowN(int IneRowN){
    //Array2DReverseIneRowN(this->data, IneRowN);
}

template<typename T>void Array2DRect_V1D<T>::Transpose(){
    int Pos1, Pos2, QExtRows=this->GetQExtRows();
    for(int i=1; i<=QExtRows; i++){
        for(int j=1; j<=this->ExtRowLength; j++){
            Pos1=this->GetPositionByCoords(i, j);
            Pos2=this->GetPositionByCoords(j, i);
            Array1DSwapElements(this->data, Pos1, Pos2);
        }
    }
}

template<typename T>void Array2DRect_V1D<T>::DelExtRowN(int ExtRowN){
    int Pos1=this->GetPositionByCoords(ExtRowN, 1),
        Pos2=this->GetPositionByCoords(ExtRowN, this->ExtRowLength);
    for(int i=1; i<=Pos2-Pos1; i++){
        Array1DDelElementFromN(this->data, Pos1);
    }
}

template<typename T> void Array2DRect_V1D<T>::DelIneRowN(int IneRowN, bool ifPos0IgnoreNotDelLastOfVaria){
    int QExtRows=this->GetQExtRows(),
        PosN;
    for(int i=QExtRows; i>=1; i--){
        PosN=this->GetPositionByCoords(i, IneRowN);
        Array1DDelElementFromN(this->data, PosN);
    }
}

/*template<typename T> void Array2DRect_V1D<T>::AddToExtRow(T val, int ExtRowN){
    int QExtRows=this->GetQExtRows();
    if(ExtRowN<QExtRows){

    }
}

template<typename T> void Array2DRect_V1D<T>::InsToExtRow(T val, int ExtRowN, int IneRowPosN){
    //template<typename T>Array2DInsToExtRowN(std::vector<std::vector<T>>&X, int ExtRowN, T val, int IneRowPosN){
    Array2DInsToExtRowN(this->data, ExtRowN, val, IneRowPosN);
}

template<typename T> void Array2DRect_V1D<T>::DelFromExtRow(int ExtRowN, int IneRowPosN){
    //template<typename T>Array2DDelFromExtRowN(std::vector<std::vector<T>>&X, int ExtRowN, int IneRowPosN){
    Array2DDelFromExtRowN(this->data, ExtRowN, IneRowPosN);
}*/


template<typename T>  void Array2DRect_V1D<T>::SetExtRow(int ExtRowN, std::vector<T>rowParam, int whatFromN, int QDefaultValuesBefore, bool KeepIfRect, T*DefaultValParam, bool DefaultIsSpecNotOwn){//, int ifDiff_ByArr0ByNew1Min2Max3=0){
    //template<typename T> void Array2DSetExtRowN(std::vector<std::vector<T>>&X, int ExtRowN, std::vector<T>rowParam, int whatFromN=1, int QDefaultValuesBefore=0, bool KeepIfRect=true, T*DefaultValParam=NULL, bool DefaultIsSpecNotOwn=false){//, int ifDiff_ByArr0ByNew1Min2Max3=0){
    Array2DSetExtRowN(this->data, ExtRowN, rowParam, whatFromN, QDefaultValuesBefore, KeepIfRect, DefaultValParam, DefaultIsSpecNotOwn);
}
template<typename T>  void Array2DRect_V1D<T>::SetExtRow(int ExtRowN, int wholeRowL, T*rowParam, int FromN, int QDefaultValuesBefore, bool KeepIfRect, T*DefaultValParam, bool DefaultIsSpecNotOwn){//, int ifDiff_ByArr0ByNew1Min2Max3=0){
    //template<typename T> void Array2DSetExtRowN(std::vector<std::vector<T>>&X, int ExtRowN, int wholeRowL=0, T*rowParam=NULL, int FromN=1, int QDefaultValuesBefore=0, bool KeepIfRect=true, T*DefaultValParam=NULL, bool DefaultIsSpecNotOwn=false){//, int ifDiff_ByArr0ByNew1Min2Max3=0){
    Array2DSetExtRowN(this->data, ExtRowN, wholeRowL, rowParam, FromN, QDefaultValuesBefore, KeepIfRect, DefaultValParam, DefaultIsSpecNotOwn);
}

template<typename T>  void Array2DRect_V1D<T>::SetIneRow(int IneRowN, std::vector<T>rowParam, int whatFromN, int QDefaultValsBefore, T*DefaultValParam, bool ExistingInsteadOfDefault, bool ExistingOnlyForGTLminNotIgnore, bool LastPossIfNotRectAtPos0NotIgnore){
    //template<typename T>Array2DSetIneRow(std::vector<std::vector<T>>&X, int IneRowN, std::vector<T>rowParam, int whatFromN=1, int QDefaultValsBefore=0, T*DefaultValParam=NULL, bool ExistingInsteadOfDefault=true, bool ExistingOnlyForGTLminNotIgnore=false, bool LastPossIfNotRectAtPos0NotIgnore=false){
    Array2DSetIneRow(this->data, IneRowN, rowParam, whatFromN, QDefaultValsBefore, DefaultValParam, ExistingInsteadOfDefault, ExistingOnlyForGTLminNotIgnore, LastPossIfNotRectAtPos0NotIgnore);
}
template<typename T>  void Array2DRect_V1D<T>::SetIneRow(int IneRowN, T*rowParam, int wholeRowL,int whatFromN, int QDefaultValsBefore, T*DefaultValParam, bool ExistingInsteadOfDefault, bool ExistingOnlyForGTLminNotIgnore, bool LastPossIfNotRectAtPos0NotIgnore){
    //template<typename T>Array2DSetIneRow(std::vector<std::vector<T>>&X, int IneRowN, T*rowParam, int wholeRowL=0,int whatFromN=1, int QDefaultValsBefore=0, T*DefaultValParam=NULL, bool ExistingInsteadOfDefault=true, bool ExistingOnlyForGTLminNotIgnore=false, bool LastPossIfNotRectAtPos0NotIgnore=false){
     Array2DSetIneRow(this->data, IneRowN, rowParam, wholeRowL,whatFromN, QDefaultValsBefore, DefaultValParam, ExistingInsteadOfDefault, ExistingOnlyForGTLminNotIgnore, LastPossIfNotRectAtPos0NotIgnore);
}

template<typename T> void Array2DRect_V1D<T>::AddExtRow(std::vector<T>rowParam,  int whatFromN, int QDefaultValsBefore, bool RectNotVar, T*DfltValParam, TValsShowHide*vsh){
    //template<typename T>void Array2DAddExtRow(std::vector<std::vector<T>>&X, std::vector<T>rowParam,  int whatFromN=1, int QDefaultValsBefore=0, bool RectNotVar=true, T*DfltValParam=NULL, TValsShowHide*vsh=NULL){
    Array2DAddExtRow(this->data, rowParam,  whatFromN, QDefaultValsBefore, RectNotVar, DfltValParam, vsh);
}
template<typename T> void Array2DRect_V1D<T>::AddExtRow(int wholeRowL, T*rowParam,  int FromN, int QDefaultValsBefore, bool RectNotVar, T*DfltValParam, TValsShowHide*vsh){
    //template<typename T>void Array2DAddExtRow(std::vector<std::vector<T>>&X, int wholeRowL=0, T*rowParam=NULL,  int FromN=1, int QDefaultValsBefore=0, bool RectNotVar=true, T*DfltValParam=NULL, TValsShowHide*vsh=NULL){
    Array2DAddExtRow(this->data, wholeRowL, rowParam,  FromN, QDefaultValsBefore, RectNotVar, DfltValParam, vsh);
}

template<typename T> void Array2DRect_V1D<T>::InsExtRow(int ExtRowN, std::vector<T>rowParam, T*DefaultValParam, int whatFromN, int QDefaultValuesBefore, bool KeepIfRect){
    //template<typename T> void Array2DInsExtRow(std::vector<std::vector<T>>&X, int ExtRowN, std::vector<T>rowParam, T*DefaultValParam=NULL, int whatFromN=1, int QDefaultValuesBefore=0, bool KeepIfRect=true){
    Array2DInsExtRow(this->data, ExtRowN, rowParam, DefaultValParam, whatFromN, QDefaultValuesBefore, KeepIfRect);
}
template<typename T>  void Array2DRect_V1D<T>::InsExtRow(int ExtRowN, int wholeRowL, T*rowParam, T*DefaultValParam, int whatFromN, int QDefaultValuesBefore, bool KeepIfRect){
    //template<typename T> void Array2DInsExtRow(std::vector<std::vector<T>>&X, int ExtRowN, int wholeRowL=0, T*rowParam=NULL, T*DefaultValParam=NULL, int whatFromN=1, int QDefaultValuesBefore=0, bool KeepIfRect=true){
    Array2DInsExtRow(this->data, ExtRowN, wholeRowL, rowParam, DefaultValParam, whatFromN, QDefaultValuesBefore, KeepIfRect);
}

template<typename T> void Array2DRect_V1D<T>::AddIneRow(std::vector<T>rowParam, int IfNonRectIgnore0Add1Stretch2, int whatFromN, int QDefaultValsBefore, T*DefaultValParam){
    //template<typename T>Array2DAddIneRow(std::vector<std::vector<T>>&X, std::vector<T>rowParam, int IfNonRectIgnore0Add1Stretch2=0, int whatFromN=1, int QDefaultValsBefore=0, T*DefaultValParam=NULL){
    Array2DAddIneRow(this->data, rowParam, IfNonRectIgnore0Add1Stretch2, whatFromN, QDefaultValsBefore, DefaultValParam);
}
template<typename T> void Array2DRect_V1D<T>::AddIneRow(T*rowParam, int wholeRowL, int IfNonRectIgnore0Add1Stretch2, int whatFromN, int QDefaultValsBefore, T*DefaultValParam){
    //template<typename T>Array2DAddIneRow(std::vector<std::vector<T>>&X, T*rowParam=NULL, int wholeRowL=0, int IfNonRectIgnore0Add1Stretch2=0, int whatFromN=1, int QDefaultValsBefore=0, T*DefaultValParam=NULL){
    Array2DAddIneRow(this->data, rowParam, wholeRowL, IfNonRectIgnore0Add1Stretch2, whatFromN, QDefaultValsBefore, DefaultValParam);
}

template<typename T> void Array2DRect_V1D<T>::InsIneRow(int IneRowPosN, std::vector<T>rowParam, int IfAfterLminIgnore0Stretch2, int whatFromN, int QDefaultValsBefore, T*DefaultValParam){
    //template<typename T>Array2DInsIneRow(std::vector<std::vector<T>>&X, int IneRowPosN, std::vector<T>rowParam, int IfAfterLminIgnore0Stretch2=0, int whatFromN=1, int QDefaultValsBefore=0, T*DefaultValParam=NULL){
    Array2DInsIneRow(this->data, IneRowPosN, rowParam, IfAfterLminIgnore0Stretch2, whatFromN, QDefaultValsBefore, DefaultValParam);
}
template<typename T> void Array2DRect_V1D<T>::InsIneRow(int IneRowPosN, T*rowParam, int wholeRowL,  int IfAfterLminIgnore0Stretch2, int whatFromN, int QDefaultValsBefore, T*DefaultValParam){
    //template<typename T>Array2DInsIneRow(std::vector<std::vector<T>>&X, int IneRowPosN, T*rowParam=NULL, int wholeRowL=0,  int IfAfterLminIgnore0Stretch2=0, int whatFromN=1, int QDefaultValsBefore=0, T*DefaultValParam=NULL){
     Array2DInsIneRow(this->data, IneRowPosN, rowParam, wholeRowL,  IfAfterLminIgnore0Stretch2, whatFromN, QDefaultValsBefore, DefaultValParam);
}

template<typename T> std::vector<std::vector<int>> Array2DRect_V1D<T>::SeekVal(const T& Val, int paramN){
    return Array2D_SeekVal_Simply(this->data, this->size, Val, paramN);
}

//template<typename T> std::vector<std::vector<int>> SeekExtRow(T*row, int L, int paramN=0){
      //return Array2D_SeekExtRow_Simply(this->data, this->size,
//}

//template<typename T> std::vector<std::vector<int>> SeekExtRow(std::vector<T>row, int paramN=0){
     //return Array2D_SeekExtRow_Simply(this->data, this->size, row, paramN);
//}

//std::vector<std::vector<int>> SeekIneRow(T*row, int L);
 //std::vector<std::vector<int>> SeekIneRow(std::vector<T>row);
template<typename T> std::vector<std::vector<int>> Array2DRect_V1D<T>::SeekSubArray(std::vector<std::vector<T>>What, bool IsTransposed, bool EachExtRowIsReversed ){
    //template <typename T> std::vector<std::vector<int>> Array2D_SeekArrs2ndIn1st(std::vector<std::vector<T>>Where,  std::vector<std::vector<T>>What, bool IsTransposed=false, bool EachExtRowIsReversed=false){
    return Array2D_SeekArrs2ndIn1st(this->data,  What, IsTransposed, EachExtRowIsReversed);
}

//template<typename T> std::vector<std::vector<int>> Array2DRect_V1D<T>::SeekSubArray(T**data, const Array2DSize& size, bool Transposed, bool ExtReversed, bool IneReversed);

template<typename T> Array2DRect_V1D<T> Array2DRect_V1D<T>::GetSubArray(int ExtRowN1, int IneRowN1, int ExtRowN2, int IneRowN2){
    std::vector< std::vector<T>>data;
    data=this->Get();
    return*this;//do!
}
template<typename T> Array2DRect_V1D<T> Array2DRect_V1D<T>::GetSubArray(const Array2DSize&size){
    std::vector< std::vector<T>>data;
    data=this->Get();
    return*this;//do!
}

    //
template<typename T> void Array2DRect_V1D<T>::ConcatAdding(){}

template<typename T> void Array2DRect_V1D<T>::ConcatStretching(){}

template<typename T> void Array2DRect_V1D<T>::Arr2DStdCOut2D(QString delimElem, QString delimLines, bool useBrackets){
    //template <typename T>void Arr2DStdCOut2D(std::vector<std::vector<T>>X, QString delimElem="; ", QString delimLines="; ", bool useBrackets=true){
    std::vector< std::vector<T>>data;
    data=this->Get();
    Arr2DStdCOut2D(data, delimElem, delimLines, useBrackets);
}

template<typename T> void Array2DRect_V1D<T>::Arr2DStdCOut1D(QString delimElem, QString delimLines, bool useBrackets){
    //template <typename T>void Arr2DStdCOut1D(std::vector<std::vector<T>>X,  QString delimElem=", ", QString delimLines="; ", bool useBrackets=true){
    std::vector< std::vector<T>>data;
    data=this->Get();
    Arr2DStdCOut1D(this->data, delimElem, delimLines, useBrackets);
}


//=====================================================================================================================================================

template<typename T> class Array2D_V1D{
    std::vector<T>data;
    Array2DSize size;
    int ExtRowsLengthSame;
public:
    //Array2DPS(int QExtRows=1, int IneRowsLength=1, bool VarNotRect=true, T*DefaultValParam=NULL);
    Array2D_V1D(const Array2DSize&size, T**X=NULL, T*DefaultValParam=NULL);
    Array2D_V1D(std::vector<std::vector<T>>X);
    Array2D_V1D(int QExtRows=1, int IneRowsLength=1, T**X=NULL, T*DefaultValParam=NULL);
    Array2D_V1D(int Q, bool QRowsNotLengthOfOneExt=true, T*X=NULL, T*DefaultValParam=NULL, bool RectNotVar=true);
    Array2D_V1D(std::vector<T>X, bool QRowsNotLengthOfOneExt=true, bool RectNotVar=true);
    Array2D_V1D(const Array2D_V1D&obj);
    ~Array2D_V1D();
    //
    void Assign(const Array2D_V1D&X);
    Array2D_V1D& operator = (const Array2D_V1D&obj);
    void SetNull();
    //
    const std::vector<T> operator[](int index) const;
    //
    //void Set(int QExtRows=1, int IneRowsLength=1, bool VarNotRect=true, T*DefaultValParam=NULL);
    void Set(const Array2DSize&size, T**X=NULL, T*DefaultValParam=NULL);
    void Set(std::vector<std::vector<T>>X);
    void Set(int QExtRows=1, int IneRowsLength=1, T**X=NULL, T*DefaultValParam=NULL);
    void Set(int Q, bool QRowsNotLengthOfOneExt=true, T*X=NULL, T*DefaultValParam=NULL, bool RectNotVar=true);
    void Set(std::vector<T>X, bool QRowsNotLengthOfOneExt=true, bool RectNotVar=true);
    //
    void SetSize(const Array2DSize&newSize, T*DefaultValParam=NULL);
    void SetSize(int QExtRows=1, int IneRowsLength=1, T*DefaultValParam=NULL, bool RectNotVar=true);
    //
    std::vector<std::vector<T>> Get();
    std::vector<T> GetExtRowN_asVect(int ExtRowN);
    std::vector<T> GetIneRowN_asVect(int IneRowN);
    T*GetExtRowN_asPtr(int ExtRowN);
    //
    Array2DSize GetSize() const;
    int GetQExtRows() const;
    int GetLength(int ineRowN) const;
    int GetMinLength() const;
    int GetMaxLength() const;
    //
    int GetPositionByCoords(int ExtRowN, int IneRowN);
    void GetCoordsByPosition(int PosN, int&ExtRowN, int&IneRowN);
    int GetExtRowNByPosition(int PosN);
    int GetIneRowNByPosition(int PosN);
    bool isRectangular();
    //Array2DSize GetSize();
    //
    void SetRectangular(int L_IfNonRectMinM1MaxM2=-2);
    void SetVaria();
    //
    void SetElement(const T& val, int ExtRowN, int IneRowN=1);
    //
    T GetElement_AsVal(int ExtRowN, int IneRowN=1) const;
    T*GetElement_AsPtr(int ExtRowN, int IneRowN=1);
    //
    void SwapVals(int ExtRowN1, int IneRowN1, int ExtRowN2, int IneRowN2);
    void SwapExtRows(int ExtRowN1, int ExtRowN2);
    void SwapIneRows(int IneRowN1, int IneRowN2);
    void ReverseExtRows();
    void ReverseIneRows();
    void ReverseExtRowN(int ExtRowN);
    void ReverseIneRowN(int IneRowN);
    void Transpose();
    //void TransposeStruct();
    //
    void DelExtRowN(int ExtRowN);
    void DelIneRowN(int IneRowN, bool ifPos0IgnoreNotDelLastOfVaria=true);
    //
    void SetExtRow(int ExtRowN, std::vector<T>rowParam, int whatFromN=1, int QDefaultValuesBefore=0, int IfRectAndOtherLength_Ignore0_MakeNotRect1_AddWithArrayLength2_MakeRectWithAddedRowLength3=2, T*DefaultValParam=NULL, bool DefaultIsSpecNotOwn=false);
    void SetExtRow(int ExtRowN, int wholeRowL=0, T*rowParam=NULL, int FromN=1, int QDefaultValuesBefore=0, bool KeepIfRect=true, T*DefaultValParam=NULL, bool DefaultIsSpecNotOwn=false);

    void SetIneRow(int IneRowN, std::vector<T>rowParam, int whatFromN=1, int QDefaultValsBefore=0, T*DefaultValParam=NULL, bool ExistingInsteadOfDefault=true, bool ExistingOnlyForGTLminNotIgnore=false, bool LastPossIfNotRectAtPos0NotIgnore=false);
    void SetIneRow(int IneRowN, T*rowParam, int wholeRowL=0,int whatFromN=1, int QDefaultValsBefore=0, T*DefaultValParam=NULL, bool ExistingInsteadOfDefault=true, bool ExistingOnlyForGTLminNotIgnore=false, bool LastPossIfNotRectAtPos0NotIgnore=false);

    void AddExtRow(std::vector<T>rowParam,  int whatFromN=1, int QDefaultValsBefore=0, bool RectNotVar=true, T*DfltValParam=NULL, TValsShowHide*vsh=NULL);
    void AddExtRow(int wholeRowL=0, T*rowParam=NULL,  int FromN=1, int QDefaultValsBefore=0, bool RectNotVar=true, T*DfltValParam=NULL, TValsShowHide*vsh=NULL);

    void InsExtRow(int ExtRowN, std::vector<T>rowParam, T*DefaultValParam=NULL, int whatFromN=1, int QDefaultValuesBefore=0, bool KeepIfRect=true);
    void InsExtRow(int ExtRowN, int wholeRowL=0, T*rowParam=NULL, T*DefaultValParam=NULL, int whatFromN=1, int QDefaultValuesBefore=0, bool KeepIfRect=true);

    void AddIneRow(std::vector<T>rowParam, int IfNonRectIgnore0Add1Stretch2=0, int whatFromN=1, int QDefaultValsBefore=0, T*DefaultValParam=NULL);
    void AddIneRow(T*rowParam=NULL, int wholeRowL=0, int IfNonRectIgnore0Add1Stretch2=0, int whatFromN=1, int QDefaultValsBefore=0, T*DefaultValParam=NULL);

    void InsIneRow(int IneRowPosN, std::vector<T>rowParam, int IfAfterLminIgnore0Stretch2=0, int whatFromN=1, int QDefaultValsBefore=0, T*DefaultValParam=NULL);
    void InsIneRow(int IneRowPosN, T*rowParam=NULL, int wholeRowL=0,  int IfAfterLminIgnore0Stretch2=0, int whatFromN=1, int QDefaultValsBefore=0, T*DefaultValParam=NULL);
    //
    void AddToExtRow(T val, int ExtRowN);
    void InsToExtRow(T val, int ExtRowN, int IneRowPosN);
    void DelFromExtRow(int ExtRowN, int IneRowPosN);
    //
    std::vector<std::vector<int>> SeekVal(const T& Val, int paramN=0);

    //std::vector<std::vector<int>> SeekExtRow(T*row, int L, int paramN=0){
        //return Array2D_SeekExtRow_Simply(this->data, this->size,
    //}

    //std::vector<std::vector<int>> SeekExtRow(std::vector<T>row, int paramN=0){
        //return Array2D_SeekExtRow_Simply(this->data, this->size, row, paramN);
    //}

    //std::vector<std::vector<int>> SeekIneRow(T*row, int L);
    //std::vector<std::vector<int>> SeekIneRow(std::vector<T>row);
    std::vector<std::vector<int>> SeekSubArray(std::vector<std::vector<T>>What, bool IsTransposed=false, bool EachExtRowIsReversed=false);
    std::vector<std::vector<int>> SeekSubArray(const Array2D_V1D& data, bool IsTransposed=false, bool EachExtRowIsReversed=false);

    //std::vector<std::vector<int>> SeekSubArray(T**data, const Array2DSize& size, bool Transposed, bool ExtReversed, bool IneReversed);

    //Array2D_V1D& GetSubArray_byNs();
    //Array2D_V1D& GetSubArray_byLs();
    Array2D_V1D& GetSubArray(int ExtRowN1, int IneRowN1, int ExtRowN2, int IneRowN2);
    Array2D_V1D& GetSubArray(const Array2DSize& size);
    //
    void ConcatAdding();

    void ConcatStretching();

    void Arr2DStdCOut2D(QString delimElem="; ", QString delimLines="; ", bool useBrackets=true);
    void Arr2DStdCOut1D(QString delimElem=", ", QString delimLines="; ", bool useBrackets=true);
};

template<typename T> Array2D_V1D<T>::Array2D_V1D(const Array2DSize&size, T**X, T*DefaultValParam){
    this->Set(size, X, DefaultValParam);
}
template<typename T> Array2D_V1D<T>::Array2D_V1D(std::vector<std::vector<T>>X){
    Array2DSize size=GetSizeOf2DVector(X);
    T val;
    this->Ns.clear();
    this->data.clear();
    int minL=size.GetMinLength(), maxL=size.GetMaxLength(), QExtRows=size.GetQExtRows(), L;
    if(minL==maxL){
        this->ExtRowsLengthSame=minL;
    }else{
        this->ExtRowsLengthSame=0;
        for(int i=1; i<=QExtRows; i++){
            L=size.GetLength(i);
            this->Ns.push_back(L);
            for(int j=1; j<=L; j++){
                val=X[i-1][j-1];
                this->data.push_back(val);
            }
        }
    }
}
template<typename T> Array2D_V1D<T>::Array2D_V1D(int QExtRows, int IneRowsLength, T**X, T*DefaultValParam){
    this->Set(QExtRows, IneRowsLength, X, DefaultValParam);
}
template<typename T> Array2D_V1D<T>::Array2D_V1D(int Q, bool QRowsNotLengthOfOneExt, T*X, T*DefaultValParam, bool RectNotVar){
    this->Set(Q, QRowsNotLengthOfOneExt, X, DefaultValParam, RectNotVar);
}
template<typename T> Array2D_V1D<T>::Array2D_V1D(std::vector<T>X, bool QRowsNotLengthOfOneExt, bool RectNotVar){
    this->Set(QRowsNotLengthOfOneExt, QRowsNotLengthOfOneExt, X);
}

template<typename T> Array2D_V1D<T>::Array2D_V1D(const Array2D_V1D&obj){
    this->Assign(obj);
}

template<typename T> Array2D_V1D<T>::~Array2D_V1D(){}

template<typename T> void Array2D_V1D<T>::Assign(const Array2D_V1D&obj){
    this->SetNull();
    this->Set(obj.data);
}

template<typename T> Array2D_V1D<T>& Array2D_V1D<T>::operator = (const Array2D_V1D&obj){
    this->Assign(obj);
    return *this;
}

template<typename T> void Array2D_V1D<T>::SetNull(){
    this->data.clear();
}

template<typename T> const std::vector<T> Array2D_V1D<T>::operator[](int index) const{
    std::vector<T>row;
    T val;
    int L;
    if(this->data.size()>0 && index>=0 && index<this->data.size() && this->data[index].size()>0){
        L=this->data[index].size();
        for(int i=1; i<=L; i++){
            val=this->GetElement_AsVal(index+1, i);
            //val=this->data[index][i-1];
            row.push_back(val);
        }
    }
    return row;
}

template<typename T> void Array2D_V1D<T>::Set(const Array2DSize&size, T**X, T*DefaultValParam){
    int Q=size.GetQExtRows(), L;
    std::vector<T>row;
    T val;
    //
    this->SetNull();
    //
    if(X!=NULL){
        for(int i=1; i<=Q; i++){
            row.clear();
            L=size.GetLength(i);
            if(L>0){
                for(int j=1; j<=L; j++){
                    val=X[i-1][j-1];
                    row.push_back(val);
                }
            }
            this->data.push_back(row);
        }
    }else{
        if(DefaultValParam!=NULL){
            for(int i=1; i<=Q; i++){
                row.clear();
                L=size.GetLength(i);
                if(L>0){
                    for(int j=1; j<=L; j++){
                        this->data[i-1][j-1]=(*(DefaultValParam));
                    }
                }
            }
        }
        //this->SetSize(size, DefaultValParam)//no, no need to preserrve vals
    }
}
template<typename T> void Array2D_V1D<T>::Set(std::vector<std::vector<T>>X){
    std::vector<T>row;
    T val;
    int Q=X.size(), L;
    this->SetNull();
    if(Q>0){
        for(int i=1; i<=Q; i++){
            row.clear();
            L=X[i-1].size();
            if(L>0){
                for(int j=1; j<=L; j++){
                    val=X[i-1][j-1];
                    row.push_back(val);
                }
            }
            this->data.push_back(row);
        }
    }
}

template<typename T> void Array2D_V1D<T>::Set(int QExtRows, int IneRowsLength, T**X, T*DefaultValParam){
    this->SetNull();
    Array2DSize size;
    size.Set(QExtRows, IneRowsLength);
    this->Set(size, X, DefaultValParam);
}
template<typename T> void Array2D_V1D<T>::Set(int Q, bool QRowsNotLengthOfOneExt, T*X, T*DefaultValParam, bool RectNotVar){
    this->SetNull();
    T CurVal, DfltVal;
    if(DefaultValParam!=NULL){
        DfltVal=(*(DefaultValParam));
    }
    if(X!=NULL){
        for(int i=1; i<=Q; i++){
            CurVal=X[i-1];
            this->data.push_back(CurVal);
        }
    }else{
        for(int i=1; i<=Q; i++){
            this->data.push_back(CurVal);
        }
    }
    this->Ls.clear();
    if(QRowsNotLengthOfOneExt){
        if(RectNotVar){
            this->ExtRowsLengthSame=1;
        }else{
            this->ExtRowsLengthSame=0;
            for(int i=1; i<=QRowsNotLengthOfOneExt; i++){
                this->Ls.push_back(1);
            }
        }
    }else{
        if(RectNotVar){
            this->ExtRowsLengthSame=QRowsNotLengthOfOneExt;
        }else{
            this->ExtRowsLengthSame=0;
            for(int i=1; i<=1; i++){
                this->Ls.push_back(QRowsNotLengthOfOneExt);
            }
        }
    }
}

template<typename T> void Array2D_V1D<T>::Set(std::vector<T>X, bool QRowsNotLengthOfOneExt, bool RectNotVar){
    T*Dflt=NULL;
    this->Set(X.size(), QRowsNotLengthOfOneExt, X.data(), Dflt, RectNotVar);
}


template<typename T>void  Array2D_V1D<T>::SetSize(const Array2DSize&newSize, T*DefaultValParam){
    T curVal, dfltVal;
    if(DefaultValParam!=NULL){
        dfltVal=(*(DefaultValParam));
    }
    int QExtRowsOld = this->GetQExtRows(), QExtRowsNew = newSize.GetQExtRows(),
        QExtRowsMin = QExtRowsOld<=QExtRowsNew ? QExtRowsOld : QExtRowsNew,
        L, QAll, Lold, Lnew, PosN, minL;
    if(QExtRowsNew>QExtRowsOld){
        for(int i=QExtRowsOld+1; i<=QExtRowsNew; i++){
           L=newSize.GetLength(i);
           for(int j=1; j<=L; j++){
               this->data.push_back(dfltVal);
           }
           this->Ls.push_back(L);
        }
    }else if(QExtRowsNew<QExtRowsOld){
        for(int i=QExtRowsNew+1; i<=QExtRowsOld; i++){
           L=this->Ls[i-1];
           for(int j=1; j<=L; j++){
               QAll=this->data.size();
               Array1DDelElementFromN(this->data, QAll);
           }
           this->Ls.push_back(L);
        }
    }
    //now ha same ext rows quantity, correcting ext rows lengthes
    for(int i=1; i<=QExtRowsMin; i++){
        Lold=this->Ls[i-1];
        Lnew=newSize.GetLength(i);
        minL = Lnew<=Lold ? Lnew : Lold;
        PosN=this->GetPositionByCoords(i, minL+1);
        if(Lnew>Lold){
            //for(int i=1; i<=QExtRowsNew; i++){
                        for(int j=1; j<=Lnew-Lold; j++){
                Array1DInsElementToN(this->data, PosN, dfltVal);
            }
            this->Ls[i-1]+=(Lnew-Lold);
            //}
        }else if(Lnew<Lold){
            //for(int i=1; i<=QExtRowsNew; i++){
            for(int j=1; j<=Lold-Lnew; j++){
                Array1DDelElementFromN(this->data, PosN);
            }
            //}
            this->Ls[i-1]-=(Lold-Lnew);
        }
    }
}

template<typename T>void Array2D_V1D<T>::SetSize(int QExtRows, int IneRowsLength, T*DefaultValParam, bool RectNotVar){
    Array2DSize newSize, oldSize=GetSizeOf2DVector(this->data);
    if(RectNotVar){
        newSize.SetVaria(QExtRows, IneRowsLength);
    }else{
        newSize.SetVaria(QExtRows, IneRowsLength);
    }
    //template<typename T> void Array2DSetSize(std::vector<std::vector<T>>&X, const Array2DSize&size, T*DfltValParam=NULL){//16
    Array2DSetSize(this->data,newSize, DefaultValParam);
}

template<typename T>std::vector<std::vector<T>> Array2D_V1D<T>:: Get(){//Get_AsVect2D(){
    std::vector<T> row;
    std::vector<std::vector<T>>data;
    T val;
    int Q=this->GetQExtRows(), L;
    for(int i=1; i<=Q; i++){
        row.clear();
        L=this->GetLength(i);
        for(int j=1; j<=L; j++){
            val=this->GetElement_AsVal(i, j);
            row.push_back(val);
        }
        data.push_back(row);
    }
    return data;
}

template<typename T>std::vector<T> Array2D_V1D<T>::GetExtRowN_asVect(int ExtRowN){
    std::vector<T> row;
    T val;
    int Q=this->size.GetQExtRows(), L, PosN;
    if(ExtRowN>=1 && ExtRowN<=Q){
        L=this->GetLength(ExtRowN);
        for(int i=1; i<=L; i++){
            PosN=this->GetPositionByCoords(ExtRowN, i);
            val=this->data[PosN-1];
            row.push_back(val);
        }
    }
    return row;
}

template<typename T>std::vector<T>  Array2D_V1D<T>::GetIneRowN_asVect(int IneRowN){
    std::vector<T> row;
    Array2DSize size=GetSizeOf2DVector(this->data);
    T val;
    int Q=this->size.GetQExtRows(), L, Lmin=size.GetMinLength(), Lmax=size.GetMaxLength();
    if(Q>0 && IneRowN>=1){
        if(IneRowN<=Lmin){
            for(int i=1; i<=Q; i++){
                val=this->data[i-1][IneRowN-1];
                row.push_back(val);
            }
        }else if(IneRowN<=Lmax){

        }
    }
    return row;
}
template<typename T> T* Array2D_V1D<T>::GetExtRowN_asPtr(int ExtRowN){
    T*R=NULL;
    int PosN;
    int Q=this->size.GetQExtRows();
    if(ExtRowN<0){
        ExtRowN=Q+ExtRowN+1;
    }
    if(ExtRowN>=1 && ExtRowN<=Q){
        PosN=this->GetPositionByCoords(ExtRowN, 1);
        R=this->data[PosN-1];
    }
    return R;
}
//
template<typename T> Array2DSize Array2D_V1D<T>::GetSize() const{
    Array2DSize size;
    if(this->isRectangular()){
        size.Set(this->GetQExtRows(), this->ExtRowsLengthSame);
    }else{
        size.Set(this->GetQExtRows(), this->Ls.data());
    }
    return size;
}

/*template<typename T> Array2DSize Array2D_V1D<T>::GetSize(){
    Array2DSize size;
    int Q=this->GetQExtRows(), L;
    if(isRectangular()){
        size.Set(Q, this->ExtRowsLengthSame);
    }else{
        size.SetNull();
        for(int i=1; i<=Q; i++){
            L=this->GetLength(i);
            size.AddExt(L);
        }
    }
    return size;
}*/


template<typename T>int  Array2D_V1D<T>::GetQExtRows() const{
    int QExtRows;
    if(this->Ls.size()==0){
        if(this->ExtRowsLengthSame>0){
            QExtRows=this->data.size()/this->ExtRowsLengthSame;
        }
    }else{
        QExtRows=this->Ls.size();
    }
    return QExtRows;
}

template<typename T>int  Array2D_V1D<T>::GetLength(int RowN) const{
    Array2DSize size=this->GetSize();
    return size.GetLength(RowN);
}
template<typename T>int Array2D_V1D<T>::GetMinLength() const{
    Array2DSize size=this->GetSize();
    return size.GetMinLength();
}
template<typename T>int Array2D_V1D<T>::GetMaxLength() const{
    Array2DSize size=this->GetSize();
    return size.GetMaxLength();
}
//
template<typename T> int Array2D_V1D<T>::GetPositionByCoords(int ExtRowN, int IneRowN){
    int Q=this->GetQExtRows(), L, L1, PosN=0;
    if(ExtRowN>=1 && ExtRowN<=Q){
        L1=this->GetLength(ExtRowN);
        if(IneRowN>=1 && IneRowN<=L1){
            for(int i=1; i<=ExtRowN-1; i++){
                L=this->GetLength(i);
                PosN+=L;
            }
            PosN+=IneRowN;
        }
    }
    return PosN;
}

template<typename T> void Array2D_V1D<T>::GetCoordsByPosition(int PosN, int&ExtRowN, int&IneRowN){
    int QER=this->GetQExtRows(), L, N=0;
    ExtRowN=0;
    bool contin=(QER>0 && PosN<=this->data.size() && N<PosN);
    while(contin){
        ExtRowN++;
        L=this->GetLength(ExtRowN);
        if(N<PosN && N+L>=PosN){
            contin=false;
            IneRowN=PosN-N;
        }else{
            N+=L;
        }
    }
}

template<typename T> int Array2D_V1D<T>::GetExtRowNByPosition(int PosN){
    int ExtRowN=0, IneRowN=0;
    this->GetCoordsByPosition(PosN, ExtRowN, IneRowN);
    return ExtRowN;
}

template<typename T> int Array2D_V1D<T>::GetIneRowNByPosition(int PosN){
        int ExtRowN=0, IneRowN=0;
        this->GetCoordsByPosition(PosN, ExtRowN, IneRowN);
        return IneRowN;
}

template<typename T> bool Array2D_V1D<T>::isRectangular(){
    bool  b=(this->Ls.size()==0);
    return b;
}


//
template<typename T> void  Array2D_V1D<T>::SetRectangular(int L_IfNonRectMinM1MaxM2){
    int Lmax=this->GetMaxLength(), Lmin=this->GetMinLength();
    if(Lmin==Lmax){
        this->ExtRowsLengthSame=Lmin;
        this->Ls.clear();
    }else{
        if(L_IfNonRectMinM1MaxM2==-1){
            L_IfNonRectMinM1MaxM2=Lmin;
        }else if(L_IfNonRectMinM1MaxM2<0){//==-1{
            L_IfNonRectMinM1MaxM2=Lmax;
        }//else L_IfNonRectMinM1MaxM2 keeps its val
        this->SetSize(this->GetQExtRows(), L_IfNonRectMinM1MaxM2 );
    }
}

template<typename T> void Array2D_V1D<T>::SetVaria(){
    int Q=this->GetQExtRows();
    this->Ls.clear();
    for(int i=1; i<=Q; i++){
        this->Ls.push_back(this->ExtRowsLengthSame);
    }
    this->ExtRowsLengthSame=0;
}

//
template<typename T> T Array2D_V1D<T>::GetElement_AsVal(int ExtRowN, int IneRowN) const{
    T val;
    int PosN;
    Array2DSize size=GetSizeOf2DVector(this->data);
    int Q=size.GetQExtRows(), L=size.GetLength(ExtRowN);
    if(ExtRowN>=1 && ExtRowN<=Q && IneRowN>=1 && IneRowN<=L){
        PosN=this->GetPositionByCoords(ExtRowN, IneRowN);
        val=this->data[PosN-1];
    }
    return val;
}
template<typename T> T* Array2D_V1D<T>::GetElement_AsPtr(int ExtRowN, int IneRowN){
    T*ptr=NULL;
    int PosN;
    Array2DSize size=GetSizeOf2DVector(this->data);
    int Q=size.GetQExtRows(), L=size.GetLength(ExtRowN);
    if(ExtRowN>=1 && ExtRowN<=Q && IneRowN>=1 && IneRowN<=L){
        PosN=this->GetPositionByCoords(ExtRowN, IneRowN);
        ptr=this->data[PosN-1];
    }
    return ptr;
}

template<typename T> void Array2D_V1D<T>::SetElement(const T& val, int ExtRowN, int IneRowN){
    int Q=this->GetQExtRows(), L=this->GetLength(ExtRowN), PosN;
    if(ExtRowN<0){
        ExtRowN=ExtRowN+Q+1;
    }
    if(IneRowN<0){
        IneRowN=IneRowN+L+1;
    }
    if(ExtRowN>=1 && ExtRowN<=Q && IneRowN>=1 && IneRowN<=L){
        PosN=this->GetPositionByCoords(ExtRowN, IneRowN);
        this->data[PosN-1]=val;
    }
}

template<typename T>void Array2D_V1D<T>::SwapVals(int ExtRowN1, int IneRowN1, int ExtRowN2, int IneRowN2){
    int Pos1, Pos2, Pos1Ini, Pos1Fin, Pos2Ini, Pos2Fin;
    Pos1=this->GetPositionByCoords(ExtRowN1, IneRowN1);
    Pos2=this->GetPositionByCoords(ExtRowN2, IneRowN2);
    Array1DSwapElements(this->data, Pos1, Pos2);
}

template<typename T>void Array2D_V1D<T>::SwapExtRows(int ExtRowN1, int ExtRowN2){
    int Pos1, Pos2, Pos1Ini, Pos1Fin, Pos2Ini, Pos2Fin, L1, L2, Lmin;
    std::vector<T>rest;
    T val;
    Pos1Ini=this->GetPositionByCoords(ExtRowN1, 1);
    Pos1Fin=this->GetPositionByCoords(ExtRowN1, this->GetLength(ExtRowN1));
    Pos2Ini=this->GetPositionByCoords(ExtRowN2, 1);
    Pos2Fin=this->GetPositionByCoords(ExtRowN2, this->GetLength(ExtRowN2));
    L1=this->GetLength(ExtRowN1);
    L2=this->GetLength(ExtRowN2);
    Lmin= L1<=L2 ? L1 : L2;
    //swapping part of min length
    for(int i=1; i<=Lmin; i++){
        Pos1=Pos1Ini+i-1;
        Pos2=Pos2Ini+i-1;
        Array1DSwapElements(this->data, Pos1, Pos2);
    }
    //moving rest
    if(L1>L2){
        for(int i=L1-L2; i>=1; i--){
            Pos1=Pos1Fin+i-1;
            Pos2=Pos2Fin+1;
            val=this->data[Pos1-1];
            Array1DInsElementToN(this->data, Pos2, val);
        }
        for(int i=1; i<=L1-L2; i++){
            Array1DDelElementFromN(this->data, Pos1Fin+1);
        }
    }else if(L2>L1){
        for(int i=L2-L1; i>=1; i--){
            Pos1=Pos1Fin+1;
            Pos2=Pos2Fin+i-1;
            val=this->data[Pos2-1];
            Array1DInsElementToN(this->data, Pos1, val);
        }
        for(int i=1; i<=L2-L1; i++){
            Array1DDelElementFromN(this->data, Pos2Fin+1);
        }
    }
    //wrting into lengthes array
    Array1DSwapElements(this->Ls, Pos1, Pos2);
}

template<typename T>void Array2D_V1D<T>::SwapIneRows(int IneRowN1, int IneRowN2){
    int minL=this->GetMinLength(), QExtRows=this->GetQExtRows(), Pos1, Pos2;
    if(IneRowN1>=1 && IneRowN1<=minL && IneRowN2>=1 && IneRowN2<=minL){
        for(int i=1; i<=QExtRows; i++){
            Pos1=this->GetPositionByCoords(i, IneRowN1);
            Pos2=this->GetPositionByCoords(i, IneRowN2);
            Array1DSwapElements(this->data, Pos1, Pos2);
        }
    }
}

template<typename T>void Array2D_V1D<T>::ReverseExtRows(){
    int PosIni, PosFin, Pos1, Pos2, PosMed, QExtRows=this->GetQExtRows(), L;
    for(int i=1; i<=QExtRows; i++){
        L=this->GetLength(i);
        PosIni= this->GetPositionByCoords(i, 1);
        PosFin= this->GetPositionByCoords(i, L);
        if(L%2==0){
            PosMed=this->GetPositionByCoords(i, L/2);
        }else{
            PosMed=this->GetPositionByCoords(i, (L-1)/2);
        }
        for(int j=PosIni; j<=PosMed; j++){
            Pos1=j;
            Pos2=PosFin-j+1;
            Array1DSwapElements(this->data, Pos1, Pos2);
        }
    }
}

template<typename T>void Array2D_V1D<T>::ReverseIneRows(){
    Array2DSize size=this->GetSize();
    int QExtRows=this->GetQExtRows(), Lmin=this->GetMinLength(), Lmax=this->GetMaxLength(), PosIni, PosFin, Pos1, Pos2, PosMed, Hmed;
    Hmed = (QExtRows%2==0) ? QExtRows/2 : (QExtRows-1)/2;
    if(Lmin==Lmax){
        for(int j=1; j<=Lmin; j++){
            PosIni=this->GetPositionByCoords(1, j);
            PosFin=this->GetPositionByCoords(QExtRows, j);
            if(QExtRows%2==0){
                PosMed=this->GetPositionByCoords(QExtRows/2, j);
            }else{
                PosMed=this->GetPositionByCoords((QExtRows-1)/2, j);
            }
            for(int i=PosIni; i<=PosMed; i++){
                Pos1=i;
                Pos2=PosFin-i+1;
                Array1DSwapElements(this->data, Pos1, Pos2);
            }
        }
    }
}

template<typename T>void Array2D_V1D<T>::ReverseExtRowN(int ExtRowN){
    int PosIni, PosFin, Pos1, Pos2, PosMed, QExtRows=this->GetQExtRows(), L;
    if(ExtRowN>=1 && ExtRowN<=QExtRows){
    //for(int i=1; i<=QExtRows; i++){
        L=this->GetLength(ExtRowN);
        PosIni= this->GetPositionByCoords(ExtRowN, 1);
        PosFin= this->GetPositionByCoords(ExtRowN, L);
        if(L%2==0){
            PosMed=this->GetPositionByCoords(ExtRowN, L/2);
        }else{
            PosMed=this->GetPositionByCoords(ExtRowN, (L-1)/2);
        }
        for(int i=PosIni; i<=PosMed; i++){
            Pos1=i;
            Pos2=PosFin-i+1;
            Array1DSwapElements(this->data, Pos1, Pos2);
        }
    }
}

template<typename T>void Array2D_V1D<T>::ReverseIneRowN(int IneRowN){
    Array2DSize size=this->GetSize();
    int QExtRows=this->GetQExtRows(), Lmin=this->GetMinLength(), Lmax=this->GetMaxLength(), PosIni, PosFin, Pos1, Pos2, PosMed, Hmed;
    Hmed = (QExtRows%2==0) ? QExtRows/2 : (QExtRows-1)/2;
    if(Lmin==Lmax){
        if(IneRowN>=1 && IneRowN<=Lmin){
        //for(int j=1; j<=Lmin; j++){
            PosIni=this->GetPositionByCoords(1, IneRowN);
            PosFin=this->GetPositionByCoords(QExtRows, IneRowN);
            if(QExtRows%2==0){
                PosMed=this->GetPositionByCoords(QExtRows/2, IneRowN);
            }else{
                PosMed=this->GetPositionByCoords((QExtRows-1)/2, IneRowN);
            }
            for(int i=PosIni; i<=PosMed; i++){
                Pos1=i;
                Pos2=PosFin-i+1;
                Array1DSwapElements(this->data, Pos1, Pos2);
            }
        }
    }
}
template<typename T>void Array2D_V1D<T>::Transpose(){
    Array2DSize size=this->GetSize();
    std::vector<T>data;
    T val;
    int PosN, L, QExtRows=this->GetQExtRows(), Lmin=this->GetMinLength(), Lmax=this->GetMaxLength(), H;
    if(size.IsForTranspose()){
        for(int i=1; i<=Lmax; i++){
            H=size.LengthFullOfIneRowN(i);
            for(int j=1; j<=H; j++){
                val=this->GetElement_AsVal(j, i);
                data.push_back(val);
            }
        }
        //assigning data
        this->data=data;
        //now writing size
        if(this->isRectangular()==false){
            size.Transpose();
            this->Ls.clear();
            QExtRows=size.GetQExtRows();
            for(int i=1; i<=QExtRows; i++){
                L=size.GetLength(i);
                this->Ls.push_back(L);
            }
        }else{
            this->ExtRowsLengthSame=QExtRows;
        }
    }
}

template<typename T>void Array2D_V1D<T>::DelExtRowN(int ExtRowN){
    int L, N1;
    if(ExtRowN>=1 && ExtRowN<=this->GetQExtRows()){
        L=this->GetLength(ExtRowN);
        N1=this->GetPositionByCoords(ExtRowN, 1);
        for(int i=1; i<=L; i++){
            Array1DDelElementFromN(this->data, N1);
        }
        if(this->isRectangular()==false){
            Array1DDelElementFromN(this->Ls, ExtRowN);
        }
    }
}

template<typename T> void Array2D_V1D<T>::DelIneRowN(int IneRowN, bool ifPos0IgnoreNotDelLastOfVaria){
    Array2DSize size = this->GetSize();
    int H = size.LengthFullOfIneRowN(IneRowN), Lmin = size.GetMinLength(), PosN;
    int QExtRows=this->GetQExtRows(), L;
    if(IneRowN>1 && IneRowN<=Lmin){
        for(int i=1; i<=H; i++){
            PosN=this->GetPositionByCoords(i, IneRowN);
            Array1DDelElementFromN(this->data, PosN);
        }
        if(this->isRectangular()==false){
            for(int i=1; i<=H; i++){
                this->Ls[i-1]--;
            }
        }else{
            this->ExtRowsLengthSame--;
        }
    }else if(IneRowN==0 && ifPos0IgnoreNotDelLastOfVaria==false){
        for(int i=1; i<=QExtRows; i++){
            L=this->GetLength(i);
            if(L>0){
                PosN=this->GetPositionByCoords(i, L);
                Array1DDelElementFromN(this->data, PosN);
            }
        }
        if(this->isRectangular()==false){
            for(int i=1; i<=H; i++){
                this->Ls[i-1]--;
            }
        }else{
            this->ExtRowsLengthSame--;
        }
    }
}

template<typename T> void Array2D_V1D<T>::AddToExtRow(T val, int ExtRowN){
    int QExtRows=this->GetQExtRows(), L, PosN;
    if(ExtRowN>=1 && ExtRowN<=QExtRows && this->isRectangular()==false){
        L=this->GetLength(ExtRowN);
        PosN=this->GetPositionByCoords(ExtRowN, L);
        Array1DInsElementToN(this->data, PosN, val);
        this->Ls[ExtRowN-1]++;
    }
}

template<typename T> void Array2D_V1D<T>::InsToExtRow(T val, int ExtRowN, int IneRowPosN){
    int QExtRows=this->GetQExtRows(), L, PosN;
    if(ExtRowN>=1 && ExtRowN<=QExtRows && this->isRectangular()==false){
        L=this->GetLength(ExtRowN);
        if(IneRowPosN>=1 && IneRowPosN<=L){
            PosN=this->GetPositionByCoords(ExtRowN, IneRowPosN);
            Array1DInsElementToN(this->data, PosN, val);
            this->Ls[ExtRowN-1]++;
        }
    }
}

template<typename T> void Array2D_V1D<T>::DelFromExtRow(int ExtRowN, int IneRowPosN){
    int QExtRows=this->GetQExtRows(), L, PosN;
    if(ExtRowN>=1 && ExtRowN<=QExtRows && this->isRectangular()==false){
        L=this->GetLength(ExtRowN);
        if(IneRowPosN>=1 && IneRowPosN<=L){
            PosN=this->GetPositionByCoords(ExtRowN, IneRowPosN);
            Array1DDelElementFromN(this->data, PosN);
            this->Ls[ExtRowN-1]--;
        }
    }
}

template<typename T>  void Array2D_V1D<T>::SetExtRow(int ExtRowN, std::vector<T>rowParam, int whatFromN, int QDefaultValuesBefore, int IfRectAndOtherLength_Ignore0_MakeNotRect1_AddWithArrayLength2_MakeRectWithAddedRowLength3, T*DefaultValParam, bool DefaultIsSpecNotOwn){//, int ifDiff_ByArr0ByNew1Min2Max3=0){
    T val;
    std::vector<T>row;
    int PosN, QExtRows=this->GetQExtRows(), Lrow, Lext=rowParam.size(), minL, Lreq;
    if(ExtRowN>=1 && ExtRowN<=this->GetQExtRows()){
        PosN=this->GetPositionByCoords(ExtRowN, 1);
        Lrow=this->GetLength(ExtRowN);
        minL = Lrow <= Lext ? Lrow : Lext;
        if(this->isRectangular()){
            if(rowParam.size()!=this->ExtRowsLengthSame){
                switch(IfRectAndOtherLength_Ignore0_MakeNotRect1_AddWithArrayLength2_MakeRectWithAddedRowLength3){
                    case 0:
                        //NOp;
                    break;
                    case 1:
                         Lreq=0;
                         row=Array1DGetSubArray_byLs_PreMark(rowParam, Lreq, whatFromN, QDefaultValuesBefore, DefaultValParam);
                         for(int i=1; i<=minL; i++){
                             val=row[i-1];
                             this->data[PosN+i-1]=val;
                         }
                    break;
                    case2:
                        Lreq=Lrow;
                        row=Array1DGetSubArray_byLs_PreMark(rowParam, Lreq, whatFromN, QDefaultValuesBefore, DefaultValParam);
                    break;
                    case 3:
                        Lreq=0;
                        row=Array1DGetSubArray_byLs_PreMark(rowParam, Lreq, whatFromN, QDefaultValuesBefore, DefaultValParam);
                    break;

                    default:
                        //NOp;
                    break;
                }//sw

            }
        }else{

        }
    }
    //row=Array1DGetSubArray_byLs_PreMark(rowParam,
}
template<typename T>  void Array2D_V1D<T>::SetExtRow(int ExtRowN, int wholeRowL, T*rowParam, int FromN, int QDefaultValuesBefore, bool KeepIfRect, T*DefaultValParam, bool DefaultIsSpecNotOwn){//, int ifDiff_ByArr0ByNew1Min2Max3=0){
    //template<typename T> void Array2DSetExtRowN(std::vector<std::vector<T>>&X, int ExtRowN, int wholeRowL=0, T*rowParam=NULL, int FromN=1, int QDefaultValuesBefore=0, bool KeepIfRect=true, T*DefaultValParam=NULL, bool DefaultIsSpecNotOwn=false){//, int ifDiff_ByArr0ByNew1Min2Max3=0){
    Array2DSetExtRowN(this->data, ExtRowN, wholeRowL, rowParam, FromN, QDefaultValuesBefore, KeepIfRect, DefaultValParam, DefaultIsSpecNotOwn);
}

template<typename T>  void Array2D_V1D<T>::SetIneRow(int IneRowN, std::vector<T>rowParam, int whatFromN, int QDefaultValsBefore, T*DefaultValParam, bool ExistingInsteadOfDefault, bool ExistingOnlyForGTLminNotIgnore, bool LastPossIfNotRectAtPos0NotIgnore){
    //template<typename T>Array2DSetIneRow(std::vector<std::vector<T>>&X, int IneRowN, std::vector<T>rowParam, int whatFromN=1, int QDefaultValsBefore=0, T*DefaultValParam=NULL, bool ExistingInsteadOfDefault=true, bool ExistingOnlyForGTLminNotIgnore=false, bool LastPossIfNotRectAtPos0NotIgnore=false){
    Array2DSetIneRow(this->data, IneRowN, rowParam, whatFromN, QDefaultValsBefore, DefaultValParam, ExistingInsteadOfDefault, ExistingOnlyForGTLminNotIgnore, LastPossIfNotRectAtPos0NotIgnore);
}
template<typename T>  void Array2D_V1D<T>::SetIneRow(int IneRowN, T*rowParam, int wholeRowL,int whatFromN, int QDefaultValsBefore, T*DefaultValParam, bool ExistingInsteadOfDefault, bool ExistingOnlyForGTLminNotIgnore, bool LastPossIfNotRectAtPos0NotIgnore){
    //template<typename T>Array2DSetIneRow(std::vector<std::vector<T>>&X, int IneRowN, T*rowParam, int wholeRowL=0,int whatFromN=1, int QDefaultValsBefore=0, T*DefaultValParam=NULL, bool ExistingInsteadOfDefault=true, bool ExistingOnlyForGTLminNotIgnore=false, bool LastPossIfNotRectAtPos0NotIgnore=false){
     Array2DSetIneRow(this->data, IneRowN, rowParam, wholeRowL,whatFromN, QDefaultValsBefore, DefaultValParam, ExistingInsteadOfDefault, ExistingOnlyForGTLminNotIgnore, LastPossIfNotRectAtPos0NotIgnore);
}

template<typename T> void Array2D_V1D<T>::AddExtRow(std::vector<T>rowParam,  int whatFromN, int QDefaultValsBefore, bool RectNotVar, T*DfltValParam, TValsShowHide*vsh){
    //template<typename T>void Array2DAddExtRow(std::vector<std::vector<T>>&X, std::vector<T>rowParam,  int whatFromN=1, int QDefaultValsBefore=0, bool RectNotVar=true, T*DfltValParam=NULL, TValsShowHide*vsh=NULL){
    Array2DAddExtRow(this->data, rowParam,  whatFromN, QDefaultValsBefore, RectNotVar, DfltValParam, vsh);
}
template<typename T> void Array2D_V1D<T>::AddExtRow(int wholeRowL, T*rowParam,  int FromN, int QDefaultValsBefore, bool RectNotVar, T*DfltValParam, TValsShowHide*vsh){
    //template<typename T>void Array2DAddExtRow(std::vector<std::vector<T>>&X, int wholeRowL=0, T*rowParam=NULL,  int FromN=1, int QDefaultValsBefore=0, bool RectNotVar=true, T*DfltValParam=NULL, TValsShowHide*vsh=NULL){
    Array2DAddExtRow(this->data, wholeRowL, rowParam,  FromN, QDefaultValsBefore, RectNotVar, DfltValParam, vsh);
}

template<typename T> void Array2D_V1D<T>::InsExtRow(int ExtRowN, std::vector<T>rowParam, T*DefaultValParam, int whatFromN, int QDefaultValuesBefore, bool KeepIfRect){
    //template<typename T> void Array2DInsExtRow(std::vector<std::vector<T>>&X, int ExtRowN, std::vector<T>rowParam, T*DefaultValParam=NULL, int whatFromN=1, int QDefaultValuesBefore=0, bool KeepIfRect=true){
    Array2DInsExtRow(this->data, ExtRowN, rowParam, DefaultValParam, whatFromN, QDefaultValuesBefore, KeepIfRect);
}
template<typename T>  void Array2D_V1D<T>::InsExtRow(int ExtRowN, int wholeRowL, T*rowParam, T*DefaultValParam, int whatFromN, int QDefaultValuesBefore, bool KeepIfRect){
    //template<typename T> void Array2DInsExtRow(std::vector<std::vector<T>>&X, int ExtRowN, int wholeRowL=0, T*rowParam=NULL, T*DefaultValParam=NULL, int whatFromN=1, int QDefaultValuesBefore=0, bool KeepIfRect=true){
    Array2DInsExtRow(this->data, ExtRowN, wholeRowL, rowParam, DefaultValParam, whatFromN, QDefaultValuesBefore, KeepIfRect);
}

template<typename T> void Array2D_V1D<T>::AddIneRow(std::vector<T>rowParam, int IfNonRectIgnore0Add1Stretch2, int whatFromN, int QDefaultValsBefore, T*DefaultValParam){
    //template<typename T>Array2DAddIneRow(std::vector<std::vector<T>>&X, std::vector<T>rowParam, int IfNonRectIgnore0Add1Stretch2=0, int whatFromN=1, int QDefaultValsBefore=0, T*DefaultValParam=NULL){
    Array2DAddIneRow(this->data, rowParam, IfNonRectIgnore0Add1Stretch2, whatFromN, QDefaultValsBefore, DefaultValParam);
}
template<typename T> void Array2D_V1D<T>::AddIneRow(T*rowParam, int wholeRowL, int IfNonRectIgnore0Add1Stretch2, int whatFromN, int QDefaultValsBefore, T*DefaultValParam){
    //template<typename T>Array2DAddIneRow(std::vector<std::vector<T>>&X, T*rowParam=NULL, int wholeRowL=0, int IfNonRectIgnore0Add1Stretch2=0, int whatFromN=1, int QDefaultValsBefore=0, T*DefaultValParam=NULL){
    Array2DAddIneRow(this->data, rowParam, wholeRowL, IfNonRectIgnore0Add1Stretch2, whatFromN, QDefaultValsBefore, DefaultValParam);
}

template<typename T> void Array2D_V1D<T>::InsIneRow(int IneRowPosN, std::vector<T>rowParam, int IfAfterLminIgnore0Stretch2, int whatFromN, int QDefaultValsBefore, T*DefaultValParam){
    //template<typename T>Array2DInsIneRow(std::vector<std::vector<T>>&X, int IneRowPosN, std::vector<T>rowParam, int IfAfterLminIgnore0Stretch2=0, int whatFromN=1, int QDefaultValsBefore=0, T*DefaultValParam=NULL){
    Array2DInsIneRow(this->data, IneRowPosN, rowParam, IfAfterLminIgnore0Stretch2, whatFromN, QDefaultValsBefore, DefaultValParam);
}
template<typename T> void Array2D_V1D<T>::InsIneRow(int IneRowPosN, T*rowParam, int wholeRowL,  int IfAfterLminIgnore0Stretch2, int whatFromN, int QDefaultValsBefore, T*DefaultValParam){
    //template<typename T>Array2DInsIneRow(std::vector<std::vector<T>>&X, int IneRowPosN, T*rowParam=NULL, int wholeRowL=0,  int IfAfterLminIgnore0Stretch2=0, int whatFromN=1, int QDefaultValsBefore=0, T*DefaultValParam=NULL){
     Array2DInsIneRow(this->data, IneRowPosN, rowParam, wholeRowL,  IfAfterLminIgnore0Stretch2, whatFromN, QDefaultValsBefore, DefaultValParam);
}

template<typename T> std::vector<std::vector<int>> Array2D_V1D<T>::SeekVal(const T& Val, int paramN){
    return Array2D_SeekVal_Simply(this->data, this->size, Val, paramN);
}

//template<typename T> std::vector<std::vector<int>> SeekExtRow(T*row, int L, int paramN=0){
      //return Array2D_SeekExtRow_Simply(this->data, this->size,
//}

//template<typename T> std::vector<std::vector<int>> SeekExtRow(std::vector<T>row, int paramN=0){
     //return Array2D_SeekExtRow_Simply(this->data, this->size, row, paramN);
//}

//std::vector<std::vector<int>> SeekIneRow(T*row, int L);
 //std::vector<std::vector<int>> SeekIneRow(std::vector<T>row);
template<typename T> std::vector<std::vector<int>> Array2D_V1D<T>::SeekSubArray(std::vector<std::vector<T>>What, bool IsTransposed, bool EachExtRowIsReversed ){
    //template <typename T> std::vector<std::vector<int>> Array2D_SeekArrs2ndIn1st(std::vector<std::vector<T>>Where,  std::vector<std::vector<T>>What, bool IsTransposed=false, bool EachExtRowIsReversed=false){
    return Array2D_SeekArrs2ndIn1st(this->data,  What, IsTransposed, EachExtRowIsReversed);
}

//template<typename T> std::vector<std::vector<int>> Array2D_V1D<T>::SeekSubArray(T**data, const Array2DSize& size, bool Transposed, bool ExtReversed, bool IneReversed);

template<typename T> Array2D_V1D<T>& Array2D_V1D<T>::GetSubArray(int ExtRowN1, int IneRowN1, int ExtRowN2, int IneRowN2){
    Array2D_V1D R;
    return R;
}

template<typename T> Array2D_V1D<T>& Array2D_V1D<T>::GetSubArray(const Array2DSize&size){
    Array2D_V1D R;
    return R;
}

    //
template<typename T> void Array2D_V1D<T>::ConcatAdding(){}

template<typename T> void Array2D_V1D<T>::ConcatStretching(){}

template<typename T> void Array2D_V1D<T>::Arr2DStdCOut2D(QString delimElem, QString delimLines, bool useBrackets){
    //template <typename T>void Arr2DStdCOut2D(std::vector<std::vector<T>>X, QString delimElem="; ", QString delimLines="; ", bool useBrackets=true){
    std::vector< std::vector<T>>data;
    data=this->Get();
    Arr2DStdCOut2D(data, delimElem, delimLines, useBrackets);
}

template<typename T> void Array2D_V1D<T>::Arr2DStdCOut1D(QString delimElem, QString delimLines, bool useBrackets){
    //template <typename T>void Arr2DStdCOut1D(std::vector<std::vector<T>>X,  QString delimElem=", ", QString delimLines="; ", bool useBrackets=true){
    std::vector< std::vector<T>>data;
    data=this->Get();
    Arr2DStdCOut2D(data, delimElem, delimLines, useBrackets);
}

//=================================================================================================

const int Array2D_TypeN_P2S=121;
const int Array2D_TypeN_P2L=122;
const int Array2D_TypeN_V2S=221;
const int Array2D_TypeN_V2L=222;
const int Array2D_TypeN_P1S=111;
const int Array2D_TypeN_P1L=112;
const int Array2D_TypeN_V1S=211;
const int Array2D_TypeN_V1L=212;

const int Array2D_TypeN_P2S_Rect=1211;
const int Array2D_TypeN_P2L_Rect=1221;
const int Array2D_TypeN_V2S_Rect=2211;
const int Array2D_TypeN_V2L_Rect=2221;
const int Array2D_TypeN_P1S_Rect=1111;
const int Array2D_TypeN_P1L_Rect=1121;
const int Array2D_TypeN_V1S_Rect=2111;
const int Array2D_TypeN_V1L_Rect=2121;

//const int Array2D_TypeN_P2S_VLO=1212;
//const int Array2D_TypeN_P2L_VLO=1222;
//const int Array2D_TypeN_V2S_VLO=2212;
//const int Array2D_TypeN_V2L_VLO=2222;
//const int Array2D_TypeN_P1S_VLO=1112;
//const int Array2D_TypeN_P1L_VLO=1122;
//const int Array2D_TypeN_V1S_VLO=2112;
//const int Array2D_TypeN_V1L_VLO=2122;

template <class T>class Array2D_Prototype{
    //
    virtual ~Array2D_Prototype()=default;
    //
    virtual void SetNull()=0;
    virtual void construct()=0;
    virtual Array2D_Prototype* CreateArray()=0;
    virtual Array2D_Prototype* clone()=0;
    virtual int GetTypeN() const=0 ;
    //
    void SetOne(double val=0);
    //
    virtual void SetSize(const Array2DSize&size)=0;
    virtual void SetSize(int QExtRows=1, int ExtRowsLength=1)=0;
    virtual void SetSize(int QExtRows=1, int*ExtRowsLengthes=NULL)=0;//if Null => L=0
    virtual void SetSize(std::vector<int>Lengthes)=0;
    //
    virtual int GetQExtRows() const=0;
    virtual int GetLength(int ExtRowN=1) const=0;
    virtual int GetMinLength() const=0;
    virtual int GetMaxLength() const=0;
    //
    bool LineNBelongsHere(int LineN);
    bool ColNBelongsHere(int ColN);
    bool NsBelongHere(int LineN, int ColN);
    //
    std::vector<T> GetArrayOfExtRowN(int LineN);
    std::vector<T> GetArrayOfIneRowN(int ColN);
    std::vector<T> GetContentAsArray1D();
    std::vector<std::vector<T>> GetContentAsArray2D();
    //
    //
    virtual T GetComponent_AsVal(int LineN, int ColN) const=0;
    virtual T* GetComponent_AsPtr(int LineN, int ColN) const=0;
    virtual void SetComponent(const T& val, int LineN, int ColN)=0;
    //
    virtual void AddToExtRowN(int RowN, const T&val)=0;
    virtual void InsToExtRowN(int RowN, int PosN, const T&val)=0;
    virtual void DelFromExtRowN(int RowN, int PosN)=0;
    //
    virtual void SetExtRowN(int N, T*rowParam, int Q, int whatN1=1, int QDfltValsBefore=0, bool DefaulsAreZerosNotOwnVals=true)=0;
    virtual void SetExtRowN(int N, std::vector<T>rowParam, int whatN1=1, int QDfltValsBefore=0, bool DefaulsAreZerosNotOwnVals=true)=0;
    virtual void SetExtRowN(int N, T*rowParam, int Q, const T& defaultVal, int whatN1=1, int QDfltValsBefore=0, bool DefaulsAreZerosNotOwnVals=true)=0;
    //
    virtual void SetIneRowN(int N, T*rowParam, int Q, int whatN1=1, int QDfltValsBefore=0, bool DefaulsAreOwnValsNotZeros=true)=0;
    virtual void SetIneRowN(int N, std::vector<T>rowParam, int whatN1=1, int QDfltValsBefore=0, bool DefaulsAreOwnValsNotZeros=true)=0;
    //
    virtual void AddExtRowN(T*rowParam, int Q, int whatN1=1, int QDfltValsBefore=0)=0;
    virtual void AddExtRowN(std::vector<T>rowParam, int whatN1=1, int QDfltValsBefore=0)=0;
    //
    virtual void AddIneRowN(T*rowParam, int Q, int whatN1=1, int QDfltValsBefore=0)=0;
    virtual void AddIneRowN(std::vector<T>rowParam, int whatN1=1, int QDfltValsBefore=0)=0;
    //
    virtual void InsExtRowN(int N, T*rowParam, int Q, int whatN1=1, int QDfltValsBefore=0)=0;
    virtual void InsExtRowN(int N, std::vector<T>rowParam, int whatN1=1, int QDfltValsBefore=0)=0;
    //
    virtual void InsIneRowN(int N, T*rowParam, int Q, int whatN1=1, int QDfltValsBefore=0)=0;
    virtual void InsIneRowN(int N, std::vector<T>rowParam, int whatN1=1, int QDfltValsBefore=0)=0;
    //
    virtual void DelExtRowN(int N=-1)=0;
    virtual void DelIneRowN(int N=-1)=0;
    //
    virtual void Set(T**y, const Array2DSize& size, bool WriteSize=true)=0;
    virtual void Set(T**y, std::vector<int>ExtRowsLengthes, bool RectIfPossible=true, bool WriteSize=true)=0;
    virtual void Set(T**y, int QExtRows, int*ExtRowsLengthes, bool RectIfPossible=true, bool WriteSize=true)=0;
    virtual void Set(T**y, int QExtRows, int ExtRowsLength=1, bool WriteSize=true)=0;
    virtual void Set(T*y, const Array2DSize& size, bool WriteSize=true)=0;
    virtual void Set(T*y, std::vector<int>ExtRowsLengthes, bool WriteSize=true)=0;
    virtual void Set(T*y, int QExtRows, int*ExtRowsLengthe, bool WriteSize=true)=0;
    virtual void Set(T*y, int QExtRows, int ExtRowsLength=1, bool sizeVaria=false, bool WriteSize=true)=0;
    virtual void Set(std::vector<std::vector<T>>y, bool RectIfPossible=true, bool WriteSize=true)=0;
    virtual void Set(std::vector<T>y, const Array2DSize& size, bool WriteSize=true)=0;
    virtual void Set(std::vector<T>y, std::vector<int>ExtRowsLengthes, bool RectIfPossible=true, bool WriteSize=true)=0;
    virtual void Set(std::vector<T>y, int QExtRows, int*ExtRowsLengthes, bool RectIfPossible=true, bool WriteSize=true)=0;
    virtual void Set(std::vector<T>y, int QExtRows, int ExtRowsLength=1, bool sizeVaria=false, bool WriteSize=true)=0;
    virtual void Set(const T&val, const Array2DSize& size, bool WriteSize=true)=0;
    virtual void Set(const T&val, std::vector<int>ExtRowsLengthes, bool RectIfPossible=true, bool WriteSize=true)=0;
    virtual void Set(const T&val, int QExtRows, int*ExtRowsLengthes, bool RectIfPossible=true, bool WriteSize=true)=0;
    virtual void Set(const T&val, int QExtRows, int ExtRowsLength=1, bool WriteSize=true)=0;
    //
    std::vector<T> GetExtRowN_AsArray(int rowN) const;
    T*GetExtRowN_AsPtr(int ExtRowN) const;
    std::vector<T> GetIneRowN_AsArray(int rowN) const;
    std::vector<T> GetContent_AsArray1D() const;
    std::vector<std::vector<T>> GetContent_AsArray2D() const;
    //
    void SwapVals(int ExtRowN1, int IneRowN1, int ExtRowN2, int IneRowN2);
    //
    virtual void SwapExtRows(int N1, int N2)=0;
    virtual void SwapIneRows(int N1, int N2)=0;
    //
    virtual void ReverseExtRowN(int rowN)=0;
    virtual void ReverseIneRowN(int rowN)=0;
    virtual void ReverseExtRowsSuccession()=0;
    virtual void ReverseIneRowsSuccession()=0;
    //
    bool ValIsAtPos(const T&val, int ExtRowN, int IneRowN) const;
    std::vector<std::vector<int>> SeekVal(const T&val) const;
    bool ExtRowIsAtPos(T*row, int L, int ExtRowN, int IneRowN) const;
    std::vector<std::vector<int>> SeekExtRow(T*row, int L) const;
    bool IneRowIsAtPos(T*row, int L, int ExtRowN, int IneRowN) const;
    std::vector<std::vector<int>> SeekIneRow(T*row, int L) const;
    bool Arr2DIsAtPos(std::vector<std::vector<T>>arr, int ExtRowN, int IneRowN, bool transposed=false, bool ExtRowsReversed=false, bool IneRowsReversed=false) const;
    std::vector<std::vector<int>> SeekArr2D(std::vector<std::vector<T>>arr, bool transposed=false, bool ExtRowsReversed=false, bool IneRowsReversed=false) const;
    //
    void ConcatAdding(Array2D_Prototype*obj, int QShiftedCols=0, int FromN=1, int QAdded=0);
    void ConcatStretching(Array2D_Prototype*obj, int QShiftedExtRowNs=0, int FromN=1, int QAdded=0);
};

//template <class T>class Array2D_P2S
//template <class T>class Array2D_P2L
//template <class T>class Array2D_V2S
//template <class T>class Array2D_V2L
//template <class T>class Array2D_P1S
//template <class T>class Array2D_P1L
//template <class T>class Array2D_V1S
//template <class T>class Array2D_V1L

#endif // MYARRAYLIB_H
