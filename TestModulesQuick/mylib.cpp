#include "mylib.h"

//MyLib::MyLib()
//{
//
//}

void writeln(TValsShowHide*w, /*const*/ char*pch){
    //QString *s=new QString(pch);
    QString s=pch;
    QString qs=QString(pch);
    if(w!=NULL){
        if(w->Show1Hide0!=0){
            if(w->ConsoleInterface)printf("%s\n",pch);
            if(w->f!=NULL)fprintf(w->f,"%s\n",pch);
            //if(w->TxtFld!=NULL)w->TxtFld->append(s);
            //if(w->TxtFld!=NULL)w->TxtFld->append(qs);
        }
    }
    //delete s;
}
void writeln(TValsShowHide*w, std::string s){
    const char*pch=s.c_str();
    QString qs=QString(pch);
    if(w!=NULL){
        if(w->Show1Hide0!=0){
            if(w->ConsoleInterface)printf("%s\n",pch);
            if(w->f!=NULL)fprintf(w->f,"%s\n",pch);
            //if(w->TxtFld!=NULL)w->TxtFld->append(s);
            //if(w->TxtFld!=NULL)w->TxtFld->append(qs);
        }
    }
}
void writeln(TValsShowHide*w, QString s){
   // const char*pch=(s.toStdString()).c_str();
    //std::string ss=s.toStdString();
    //char charr1[400], charr2[400];
    //ss=s.toStdString();
    ////charr1=s.toStdString().c_str();//incomp types
    ////charr2=ss.c_str();//incomp types
    //const char*pch1=ss.c_str();
    std::string ss=s.toStdString();
    const char*pch=ss.c_str();
    if(w!=NULL){
        if(w->Show1Hide0!=0){
            if(w->ConsoleInterface)printf("%s\n",pch);
            if(w->f!=NULL)fprintf(w->f,"%s\n",pch);
            //if(w->TxtFld!=NULL)w->TxtFld->append(s);
        }
    }
}

int BoolToInt(bool x){
    int y;
    bool b=BoolValByDefault;
    if(b){
        if(!x) y=0;
        else y=1;
    }else{
        if(x) y=1;
        else y=0;
    }
    return y;
}

int IntOfBool(bool x){return BoolToInt(x);}

bool IntToBool(int x){
    int y;
    bool b=BoolValByDefault;
    if(b){
        y=true;
        if(x==0) y=false;
    }else{
        y=false;
        if(x==1) y=true;
    }
    return y;
}

bool BoolOfInt(int x){return IntToBool(x);}

QString FloatToStr(float x){
    QString str;
    str.setNum(x);
    return str;
}
QString FloatToStr(double x){
    QString str;
    str.setNum(x);
    return str;
}
QString IntToStr(int x){
    QString str;
    str.setNum(x);
    return str;
}

QString UniqueStrValGenerator(std::vector<QString>vals, QString proposedVal, int N, QString bef, QString aft){
    QString ValFin, ValTmp, ValCur, sn;
    sn.setNum(N);
    ValTmp=bef+proposedVal+aft;//+sn;
    bool contin=true, isUnique=true;
    int Q=vals.size();
    while(contin){
        for(int i=1; i<=Q; i++){
            ValCur=vals[i-1];
            if(ValCur==ValTmp){
                isUnique=false;
            }
        }
        if(!isUnique){
           ValTmp=ValTmp+"-N"+sn;
        }
        contin=!isUnique;
    }
    ValFin=ValTmp;
    return ValFin;
}

void Calc_Array1DSubArray_byLs_MarkupNs_vFmls(int&LreqCalcd, int&preN1_, int&preN2_, int&ownN1_, int&ownN2, int&ForOwnN1_, int&ForOwnN2, int&postN1, int&postN2_, int Lini, int LreqGiven, int whatFromN, int QDefaultValuesBefore){
    int LownMax;//, LreqGiven;
    if(whatFromN<0){
        whatFromN=Lini+whatFromN+1;
    }
    //DefaultValuesBefore
    if(QDefaultValuesBefore>0){
        preN1_=1;
        if(LreqGiven <= 1 || LreqGiven >= QDefaultValuesBefore){
            preN2_=QDefaultValuesBefore;
        }else{
            preN2_=LreqGiven;
        }
    }else{
        preN1_=0;
        preN2_=0;
    }
    //Own values
    if(whatFromN <= Lini && (LreqGiven <= 1 || (LreqGiven >= 1 && LreqGiven >= QDefaultValuesBefore))) {
        ownN1_=whatFromN;
        //ForOwnN1_ = preN2_+ownN1_;
        ForOwnN1_=preN2_+1;
        //ForOwnN2 wi def'd af
    }else{
        ownN1_=0;
        ownN2=0;
        ForOwnN1_=0;
        ForOwnN2=0;//tgo S ablb def'd nu
    }
    //Lown(Max), et, if D ns'def'd, ownN2 et ForOwnN2
    if(ownN1_>=1){
        LownMax=Lini-whatFromN+1;
        if(LreqGiven <= 1 || LreqGiven >= QDefaultValuesBefore + LownMax){//unlimited
            ownN2=Lini;
            ForOwnN2 = ForOwnN1_ + LownMax - 1;//corr'd
        }else{//limited
            ForOwnN2=LreqGiven;
            ownN2=ownN1_+(ForOwnN2-ForOwnN1_+1)-1;//vikts! fml qvi'dd!
        }
    }else{
        LownMax=0;
    }
    //LreqCalcd
    //if(LreqGiven<=0){
    //    LreqCalcd=QDefaultValuesBefore+LownMax;//LownMax==0 if whatFromN>Lini
    //}else if(LreqGiven>=)
    //post
    if(ForOwnN2>0 && ForOwnN2<LreqGiven){
        postN1=ForOwnN2+1;
        postN2_=LreqGiven;
    }else{
        postN1=0;
        postN2_=0;
    }
    //LreqCalcd
    if(LreqGiven>0){
        LreqCalcd=LreqGiven;
    }else{
         LreqCalcd=QDefaultValuesBefore+LownMax;//LownMax==0 if whatFromN>Lini
    }
}
int CalcLreq_Array1DSubArray_byLs_MarkupNs_vFmls(int Lini, int whatFromN, int QDefaultValuesBefore, int LreqGiven){
    int LreqCalcd;
    int preN1_, preN2_, ownN1_, ownN2, ForOwnN1_, ForOwnN2, postN1, postN2_;
    //void Calc_Array1DSubArray_byLs_MarkupNs_vFmls(int&LreqCalcd, int&preN1_, int&preN2_, int&ownN1_, int&ownN2, int&ForOwnN1_, int&ForOwnN2, int&postN1, int&postN2_, int Lini, int LreqGiven, int whatFromN, int QDefaultValuesBefore){
    Calc_Array1DSubArray_byLs_MarkupNs_vFmls(LreqCalcd, preN1_, preN2_, ownN1_, ownN2, ForOwnN1_, ForOwnN2, postN1, postN2_, Lini, LreqGiven, whatFromN, QDefaultValuesBefore);
    return LreqCalcd;
}

void Calc_Array1DSubArray_byLs_MarkupNs_vCycls(int&LreqCalcd, int&preN1_, int&preN2_, int&ownN1_, int&ownN2, int&ForOwnN1_, int&ForOwnN2, int&postN1, int&postN2_, int Lold, int LreqGiven, int whatFromN, int QDefaultValuesBefore){
    bool contin;
    int whereCurN, whatCurN, whereQ, whatQ;
    int Lmin, LownFact;
    int LownMax;
    int Lreq=LreqGiven;
    if(Lreq<1){
        //pre
        if(QDefaultValuesBefore>0){
            preN1_=1;
            preN2_=QDefaultValuesBefore;//ob Ns s'nat'l, so Nlast=Q/L
        }else{
            preN1_=0;
            preN2_=0;
        }
        whereCurN=preN2_;
        //own
        if(whatFromN<=Lold){
            ownN1_=whatFromN;
            LownMax=Lold-whatFromN+1;
            ownN2=LownMax;//ob Ns s'nat'l, so Nlast=Q/L
            ForOwnN1_=whereCurN+1;
            ForOwnN2=Lold;
            LownFact=ownN2-ownN1_+1;
        }else{
            ownN1_=0;
            ownN2=0;
            ForOwnN1_==0;
            ForOwnN2=0;
            LownFact=0;
        }
        //post
        postN1=0;
        postN2_=0;
        //Lreq
        LreqCalcd=QDefaultValuesBefore+LownFact;
    }else{//Lreq fixed, known prerviously = Lreq s'ved ef
        //pre
        whereCurN=0;
        whatCurN=0;
        whereQ=Lreq;
        whatQ=Lold;
        contin=(whereCurN<whereQ && whatCurN<whatQ);
        //pre
        whereQ=Lreq;
        whatQ=QDefaultValuesBefore;
        contin=(whereCurN<whereQ && whatCurN<whatQ);
        if(contin){
            preN1_=whereCurN+1;
        }else{
           preN1_=0;
        }
        preN2_=preN1_;
        while(contin){
            whereCurN++;
            whatCurN++;
            contin=(whereCurN<whereQ && whatCurN<whatQ);
        }
        preN2_=whereCurN;
        //own
        whatCurN=0;
        whatQ=Lold;//whereQ=Lreq;remains
        contin=(whereCurN<whereQ && whatCurN<whatQ);
        if(contin){
           ownN1_=1;
        }else{
           ownN1_=0;
           ownN2=0;
        }
        ownN2=ownN1_;
        while(contin){
            whereCurN++;
            whatCurN++;
            contin=(whereCurN<whereQ && whatCurN<whatQ);
        }
        ownN2=whatCurN;
        if(ownN1_==0){
            LownFact=0;
        }else{
            LownFact=ownN2-ownN1_+1;
        }
        //own global
        if(LownFact==0){
            ForOwnN1_=0;
            ForOwnN2=0;
        }else{
            ForOwnN1_=preN2_+1;
            ForOwnN2=ForOwnN2+LownFact-1;
        }
        //post
        if(whereCurN<Lreq){
            postN1=whereCurN+1;
            postN2_=Lreq;
        }else{
            postN1=0;
            postN2_=0;
        }
        //Lreq ns'mut'd
        LreqCalcd=Lreq;
    }
}

/*int CalcLreq_Array1DSubArray_byLs_MarkupNs_vCycls(int Lini, int whatFromN, int QDefaultValuesBefore, int LreqGiven){

    int LreqCalcd;//, LreqGiven;
    //LreqGiven=Lnew
    ;
    return LreqCalcd;
}*/


template<typename T> std::vector<double> GetNumericSubarray_byLs(std::vector<double> x, int FromN, int QDefaultValuesBefore, int LreqGiven){
    std::vector<double>subArr;
    int preN1_=0, preN2_=0, ownN1_=0, ownN2=0, ForOwnN1_=0, ForOwnN2=0, postN1=0, postN2_=0,
        LreqCalcd=0;
    int Lini=x.size();
    Calc_Array1DSubArray_byLs_MarkupNs_vCycls(LreqCalcd, preN1_,preN2_, ownN1_, ownN2, ForOwnN1_, ForOwnN2, postN1, postN2_, Lini, LreqGiven, FromN, QDefaultValuesBefore);
    if(preN1_!=0 && preN2_!=0){
        for(int i=preN1_; i<=preN2_; i++){
            subArr.push_back(0);
        }
    }
    if(ForOwnN1_!=0 && ForOwnN2!=0){
        for(int i=ownN1_; i<=ownN2; i++){
            subArr.push_back(x[i-1]);
        }
    }
    if(postN1!=0 && postN2_!=0){
        for(int i=postN1; i<=postN2_; i++){
            subArr.push_back(0);
        }
    }
    return  subArr;
}

std::vector<double> GetNumericSubarray_byLs(std::vector<double> x, int FromN, int QDefaultValuesBefore, int LreqGiven){
    std::vector<double>subArr;
    int preN1_=0, preN2_=0, ownN1_=0, ownN2=0, ForOwnN1_=0, ForOwnN2=0, postN1=0, postN2_=0,
        LreqCalcd=0;
    int Lini=x.size();
    Calc_Array1DSubArray_byLs_MarkupNs_vCycls(LreqCalcd, preN1_,preN2_, ownN1_, ownN2, ForOwnN1_, ForOwnN2, postN1, postN2_, Lini, LreqGiven, FromN, QDefaultValuesBefore);
    if(preN1_!=0 && preN2_!=0){
        for(int i=preN1_; i<=preN2_; i++){
            subArr.push_back(0);
        }
    }
    if(ForOwnN1_!=0 && ForOwnN2!=0){
        for(int i=ownN1_; i<=ownN2; i++){
            subArr.push_back(x[i-1]);
        }
    }
    if(postN1!=0 && postN2_!=0){
        for(int i=postN1; i<=postN2_; i++){
            subArr.push_back(0);
        }
    }
    return  subArr;
}
std::vector<double> GetNumericSubarray_byLs(std::vector<double> x, std::vector<double> ArrOfDfltVals, int FromN, int QDefaultValuesBefore, int LreqGiven){
    std::vector<double>subArr;
    int preN1_=0, preN2_=0, ownN1_=0, ownN2=0, ForOwnN1_=0, ForOwnN2=0, postN1=0, postN2_=0,
        LreqCalcd=0;
    int Lini=x.size();
    Calc_Array1DSubArray_byLs_MarkupNs_vCycls(LreqCalcd, preN1_,preN2_, ownN1_, ownN2, ForOwnN1_, ForOwnN2, postN1, postN2_, Lini, LreqGiven, FromN, QDefaultValuesBefore);
    if(preN1_!=0 && preN2_!=0){
        for(int i=preN1_; i<=preN2_; i++){
            if(ArrOfDfltVals.size()<=i){
                subArr.push_back(ArrOfDfltVals[i-1]);
            }else{
                subArr.push_back(0);
            }
        }
    }
    if(ForOwnN1_!=0 && ForOwnN2!=0){
        for(int i=ownN1_; i<=ownN2; i++){
            subArr.push_back(x[i-1]);
        }
    }
    if(postN1!=0 && postN2_!=0){
        for(int i=postN1; i<=postN2_; i++){
            if(ArrOfDfltVals.size()<=i){
                subArr.push_back(ArrOfDfltVals[i-1]);
            }else{
                subArr.push_back(0);
            }
        }
    }
    return  subArr;
}


QString MySubString_byLs_Formatting_PreMark(QString Sini, int LreqGiven, int FromN, int QDefaultBefore, QString Default){
    QString R="", Scur;
    int Lold=Sini.length();
    int preN1_=0, preN2_=0, ownN1_=0, ownN2=0, ForOwnN1_=0, ForOwnN2=0, postN1=0, postN2_=0, LreqCalcd=0;
    //void Calc_Array1DSubArray_byLs_MarkupNs_vCycls(int&LreqCalcd, int&preN1_, int&preN2_, int&ownN1_, int&ownN2, int&ForOwnN1_, int&ForOwnN2, int&postN1, int&postN2_, int Lold, int LreqGiven, int whatFromN, int QDefaultValuesBefore){
    Calc_Array1DSubArray_byLs_MarkupNs_vCycls(LreqCalcd, preN1_, preN2_, ownN1_, ownN2, ForOwnN1_, ForOwnN2, postN1, postN2_, Lold, LreqGiven, FromN, QDefaultBefore);
    if(preN2_>0){
        for(int i=preN1_; i<=preN2_; i++){
            R=R+Default;
       }
    }
    //if(ownN1_>0){
    //   for(int i=ownN1_; i<=ownN2; i++){
    //       Scur=MySubString_L1(Sini, i);
    //       R=R+Scur;
    //  }
    //}//worked gut ma so zu irrational
    Scur=Sini.mid(ownN1_-1, ownN2-ownN1_+1);
    R=R+Scur;
    if(postN2_>0){
        for(int i=postN1; i<=postN2_; i++){
             R=R+Default;
        }
    }
    return R;
    //std::string ss, ss1;
    //ss1=ss.substr(ss.begin()+1,)
}

//ne tested
void Calc_Array1DSubArray_byLs_MarkupNs_vFml_FixedDfltsAfter(int&LreqCalcd, int&preN1, int&preN2, int&ownN1, int&ownN2, int&ForOwnN1, int&ForOwnN2, int&postN1, int&postN2, int Lini, int LreqGiven, int whatFromN, int QDefaultValuesAfter){
    int LownMax=Lini-whatFromN+1, QDefaultsBefore, LownFin;
    if(LreqGiven<=0){//=>unlimited
        if(ForOwnN1<=Lini){
            QDefaultsBefore=0;
            LreqCalcd=QDefaultValuesAfter+LownMax;
            //
            preN1=0;
            preN2=0;
            ownN1=whatFromN;
            ownN2=Lini;
            ForOwnN1=1;
            ForOwnN2=preN2+LownMax;
        }else{//ma nout
            ownN1=0;
            ownN2=0;
            ForOwnN1=0;
            ForOwnN2=0;
            LreqCalcd=QDefaultValuesAfter;
        }
        postN1=ForOwnN2+1;
        postN2=postN1+QDefaultValuesAfter-1;
    }else{//LreqGiven>0
        LreqCalcd=LreqGiven;
        if(LownMax+QDefaultValuesAfter<=LreqGiven){
            LownFin=LownMax;
            QDefaultsBefore=LreqGiven-LownMax-QDefaultValuesAfter;
        }else{
            if(QDefaultValuesAfter<LreqGiven){
                LownFin=LreqGiven-QDefaultValuesAfter;
            }else{
                LownFin=0;
            }
            QDefaultsBefore=0;
        }
        if(QDefaultsBefore>0){
            preN1=1;
            preN2=QDefaultsBefore;
        }else{
            preN1=0;
            preN2=0;
        }
        if(LownFin>0){
            ownN1=whatFromN;
            ownN2=ownN1+LownFin-1;
            ForOwnN1=preN2+1;
            ForOwnN2=ForOwnN1+LownFin-1;
        }else{
            ownN1=0;
            ownN2=0;
            ForOwnN1=0;
            ForOwnN2=0;
        }
        if(QDefaultsBefore>0){
            postN1=ForOwnN2+1;
            if(QDefaultsBefore>LreqGiven){
                postN2=postN1+LreqGiven-1;
            }else{
                postN2=postN1+QDefaultsBefore-1;
            }
        }else{
           postN1=0;
           postN2=0;
        }
    }
}


std::vector<double>NumbersSubRowSafe(std::vector<double> x,int LreqGiven=0, int FromN=1, int QDefaultBefore=0, double DefaultVal=0){
    std::vector<double> R;
    double CurVal;
    int Lold=x.size();
    int preN1_=0, preN2_=0, ownN1_=0, ownN2=0, ForOwnN1_=0, ForOwnN2=0, postN1=0, postN2_=0;
    int LreqCalcd=0;
    //void Calc_Array1DSubArray_byLs_MarkupNs_vFmls(int&LreqCalcd, int&preN1_, int&preN2_, int&ownN1_, int&ownN2, int&ForOwnN1_, int&ForOwnN2, int&postN1, int&postN2_, int Lini, int LreqGiven=0, int whatFromN=1, int QDefaultValuesBefore=0);
    Calc_Array1DSubArray_byLs_MarkupNs_vFmls(LreqCalcd, preN1_, preN2_, ownN1_, ownN2, ForOwnN1_, ForOwnN2, postN1, postN2_, Lold, LreqGiven, FromN, QDefaultBefore);
    if(preN2_>0){
        for(int i=preN1_; i<=preN2_; i++){
            R.push_back(DefaultVal);
       }
    }
    if(ownN1_>0){
        for(int i=ownN1_; i<=ownN2; i++){
            CurVal=x[i-1];
            R.push_back(CurVal);
       }
    }
    if(postN2_>0){
        for(int i=postN1; i<=postN2_; i++){
              R.push_back(DefaultVal);
        }
    }
    return R;
}

std::vector<int>NumbersSubRowSafe(std::vector<int> x,int LreqGiven=0, int FromN=1, int QDefaultBefore=0, int DefaultVal=0, std::vector<double>*ArrOfDfltVals=NULL){
    std::vector<int> R;//, dfltArr;
    double val;
    int  CurVal;
    int Lold=x.size();
    int preN1_=0, preN2_=0, ownN1_=0, ownN2=0, ForOwnN1_=0, ForOwnN2=0, postN1=0, postN2_=0;
    int LreqCalcd=0;
    //void Calc_Array1DSubArray_byLs_MarkupNs_vFmls(int&LreqCalcd, int&preN1_, int&preN2_, int&ownN1_, int&ownN2, int&ForOwnN1_, int&ForOwnN2, int&postN1, int&postN2_, int Lini, int LreqGiven=0, int whatFromN=1, int QDefaultValuesBefore=0);
    Calc_Array1DSubArray_byLs_MarkupNs_vFmls(LreqCalcd, preN1_, preN2_, ownN1_, ownN2, ForOwnN1_, ForOwnN2, postN1, postN2_, Lold, LreqGiven, FromN, QDefaultBefore);
    if(preN2_>0){
        if(ArrOfDfltVals!=NULL && ArrOfDfltVals->size()>=preN2_){
            for(int i=preN1_; i<=preN2_; i++){
                val=ArrOfDfltVals->at(i-1);
                R.push_back(val);
            }
        }else{
            for(int i=preN1_; i<=preN2_; i++){
                R.push_back(DefaultVal);
            }
       }
    }
    if(ownN1_>0){
        if(ArrOfDfltVals!=NULL && ArrOfDfltVals->size()>=ownN2){
            for(int i=ownN1_; i<=ownN2; i++){
                val=ArrOfDfltVals->at(i-1);
                R.push_back(val);
            }
        }else{
            for(int i=ownN1_; i<=ownN2; i++){
                CurVal=x[i-1];
                R.push_back(CurVal);
            }
       }
    }
    if(postN2_>0){
        for(int i=postN1; i<=postN2_; i++){
              R.push_back(DefaultVal);
        }
    }
    return R;
}


QString MySubString_byNs(QString where, int N1, int N2){
    int L=where.length();
    int ownN1, ownN2, ForOwnN1, ForOwnN2, PreN1, PreN2, PostN1, PostN2;
    QString sr="", s1;
    std::vector<QString>iniS, newS;
    for(int i=1; i<=L; i++){
        s1=MySubString_L1(where, N1);
        iniS.push_back(s1);
    }
    if(N1<0)N1=L+N1+1;
    if(N2<0)N2=L+N2+1;
    if(N1>=1 && N1<=L && N2>=1 && N2<=L){

    }
    return sr;
}

//QString MySubString_byLs(QString where, int N1, int L);
//QString MySubString_byLs_Formatting(QString where, int N1, int L, QString defaultVal){}
QString MySubString_L1(QString where, int N1){ return where.mid(N1-1, 1);}

QString MySubstring_Formatting_PreMark(QString Sini, int LreqGiven, int FromN, int QDefaultBefore, QString Default){
    //int preN1_=0, preN2_=0, ownN1_=0, ownN2=0, ForOwnN1_=0, ForOwnN2=0, postN1=0, postN2_=0,
    //LreqCalcd=0;
    //int Lini=Sini.length();
    //Calc_Array1DSubArray_byLs_MarkupNs_vCycls(LreqCalcd, preN1_,preN2_, ownN1_, ownN2, ForOwnN1_, ForOwnN2, postN1, postN2_, Lini, LreqGiven, whatFromN, QDefaultValuesBefore);
    //
    QString SCur, SR="";
    int Lold=Sini.length();
    int preN1_=0, preN2_=0, ownN1_=0, ownN2=0, ForOwnN1_=0, ForOwnN2=0, postN1=0, postN2_=0;
    int LreqCalcd=0;
    Calc_Array1DSubArray_byLs_MarkupNs_vFmls(LreqCalcd, preN1_, preN2_, ownN1_, ownN2, ForOwnN1_, ForOwnN2, postN1, postN2_, Lold, LreqGiven, FromN, QDefaultBefore);
    if(preN1_!=0 && preN2_!=0){
        for(int i=preN1_; i<=preN2_; i++){
            SCur=Default;
            SR+=SCur;
        }
    }
    if(ownN1_!=0 && ownN2!=0){
        for(int i=ownN1_; i<=ownN2; i++){
            SCur=MySubString_L1(Sini, i);
            SR+=SCur;
        }
    }
    if(postN1!=0 && postN2_!=0){
        for(int i=postN1; i<=postN2_; i++){
            SCur=Default;
            SR+=SCur;
        }
    }
    return SR;
}

void CalcNsOfSubArray(int whereL, int whatL, int&FromN1, int&ToN1, int FromN0, int ToN0){
    if(whereL>0 && whatL>0 && whereL>=whatL){
        if(FromN0<0){
            FromN1=FromN0+whereL+1;
        }
        if(ToN0<0){
            ToN1=ToN0+whereL+1;
        }else if(ToN0==0){
            ToN1=whereL;
        }
        if(FromN1>=1 && FromN1<=whereL && ToN1>=1 && ToN1<=whereL){
            if(FromN1<=ToN1){
                if(ToN1>=whereL-whatL+1){
                    ToN1=whereL-whatL+1;
                }
            }else{
                if(ToN1<=whatL){
                    ToN1=whatL;
                }
            }
        }else{
            FromN1=0;
            ToN1=0;
        }
    }else{
        FromN1=0;
        ToN1=0;
    }
}

QString ToLowerCase(QString str){ return str; }
QString ToUpperCase(QString str){ return str; }
bool SubstringIsAtPos(QString whereParam, QString whatParam, int N, bool MatchCase, bool direct){
    QString whereCur, whatCur, where, what;
    int whereL=where.length(), whatL=what.length(), whereN, whatN;
    bool b=true;
    if(N<0){
        N=whereL+N+1;
    }
    if(MatchCase){
        where=whereParam;
        what=whatParam;
    }else{
        where=ToLowerCase(whereParam);
        what=ToLowerCase(whatParam);
    }
    whereL=where.length();
    whatL=what.length();
    if(whereL>0 && whatL>0 && whereL>=whatL && N>=1 && N<=whereL){
        if(direct && N<=whereL-whatL+1){
            for(int i=1; i<=whatL; i++){
                whereCur=MySubString_L1(where, N+i-1);
                whatCur=MySubString_L1(what, i);
                if(whereCur!=whatCur){
                    b=false;
                }
            }
        }else if(direct==false && N>=whatL){
            for(int i=1; i<=whatL; i++){
                whereCur=MySubString_L1(where, N-whatL+i);
                whatCur=MySubString_L1(what, i);
                if(whereCur!=whatCur){
                    b=false;
                }
            }
        }else b=false;
    }else b=false;
    return b;
}

std::vector<int>SeekSubstring(QString where, QString what, int FromN, int ToN, bool MatchCase){
    std::vector<int> Ns;
    int whereL=where.length(), whatL=what.length();
    bool direct;
    if(FromN<0){
        FromN=FromN+whereL+1;
    }
    if(ToN<0){
        ToN=ToN+whereL+1;
    }else{
        ToN=whereL;
    }
    direct=(FromN<=ToN);
    if(FromN>=1 && FromN<=whereL && ToN>=1 && ToN<=whereL && whereL>0 && whatL>0 && whereL>=whatL){
        for(int N=FromN; N<=ToN; N++){
            if(SubstringIsAtPos(where, what, N, MatchCase, direct)==true){
                Ns.push_back(N);
            }
        }
    }
    return Ns;
}


bool IsTrueWord(QString word){
    bool b=false;
    int L=word.length();
    if(L==1){
         if(word.mid(1-1, 1)=="1"){
             b=true;
         }
         else if(word.mid(1-1, 1)=="+"){
             b=true;
         }
    }
    if(L==3){
         if(
                 (word.mid(1-1, 1)=="y" || word.mid(1-1, 1)=="Y")
                 &&
                 (word.mid(2-1, 1)=="e" || word.mid(2-1, 1)=="E")
                 &&
                 (word.mid(3-1, 1)=="s" || word.mid(3-1, 1)=="S")
         ){
             b=true;
         }
    }
    if(L==4){
         if(
                 (word.mid(1-1, 1)=="t" || word.mid(1-1, 1)=="T")
                 &&
                 (word.mid(2-1, 1)=="r" || word.mid(2-1, 1)=="R")
                 &&
                 (word.mid(3-1, 1)=="u" || word.mid(3-1, 1)=="U")
                 &&
                 (word.mid(4-1, 1)=="e" || word.mid(4-1, 1)=="E")
         ){
             b=true;
         }
    }
    return b;
}
bool IsTrueWord(char word[]){
    bool b=false;
    int L=strlen(word);
    if(L==1){
         if(word[1-1]=='1'){
             b=true;
         }
         else if(word[1-1]=='+'){
             b=true;
         }
    }
    if(L==3){
         if(
                 (word[1-1]=='y' || word[1-1]=='Y')
                 &&
                 (word[2-1]=='e' || word[2-1]=='E')
                 &&
                 (word[3-1]=='s' || word[3-1]=='S')
         ){
             b=true;
         }
    }
    if(L==4){
         if(
                 (word[1-1]=='t' || word[1-1]=='T')
                 &&
                 (word[2-1]=='r' || word[2-1]=='"')
                 &&
                 (word[3-1]=='u' || word[3-1]=='U')
                 &&
                 (word[4-1]=='e' || word[4-1]=='E')
         ){
             b=true;
         }
    }
    return b;
}
bool IsFalseWord(QString word){
    bool b=false;
    int L=word.length();
    if(L==1){
         if(word.mid(1-1, 1)=="0"){
             b=true;
         }
         else if(word.mid(1-1, 1)=="-"){
             b=true;
         }
    }
    if(L==2){
         if(
                 (word.mid(1-1, 1)=="n" || word.mid(1-1, 1)=="N")
                 &&
                 (word.mid(1-1, 1)=="o" || word.mid(1-1, 1)=="O")
         ){
             b=true;
         }
    }
    if(L==5){
         if(
                 (word.mid(1-1, 1)=="f" || word.mid(1-1, 1)=="F")
                 &&
                 (word.mid(1-1, 1)=="a" || word.mid(1-1, 1)=="A")
                 &&
                 (word.mid(1-1, 1)=="l" || word.mid(1-1, 1)=="L")
                 &&
                 (word.mid(1-1, 1)=="s" || word.mid(1-1, 1)=="S")
                 &&
                 (word.mid(1-1, 1)=="e" || word.mid(1-1, 1)=="E")
         ){
             b=true;
         }
    }
    return b;
}

bool IsFalseWord(char word[]){
    bool b=false;
    int L=strlen(word);
    if(L==1){
         if(word[1-1]=='0'){
             b=true;
         }
         else if(word[1-1]=='-'){
             b=true;
         }
    }
    if(L==2){
         if(
                 (word[1-1]=='n' || word[1-1]=='N')
                 &&
                 (word[2-1]=='o' || word[2-1]=='O')
         ){
             b=true;
         }
    }
    if(L==5){
         if(
                 (word[1-1]=='f' || word[1-1]=='F')
                 &&
                 (word[2-1]=='a' || word[2-1]=='A')
                 &&
                 (word[3-1]=='l' || word[3-1]=='L')
                 &&
                 (word[4-1]=='s' || word[4-1]=='S')
                 &&
                 (word[5-1]=='e' || word[5-1]=='E')
         ){
             b=true;
         }
    }
    return b;
}

std::vector<int>GetListOfNsByRange(int N1, int N2, int L){
    std::vector<int>Ns;
    int N2min;
    if(N1==0)N1=L;
    if(N2==0) N2=N1;
    if(N1<0)N1=L+N1+1;
    if(N2<0) N2=L+N2+1;
    if(L==0) L=N2;
    if(N2<=L)N2min=N2;
    else N2min=L;
    //
    if(N2>0){
       if(N1<0){
           for(int i=1; i<=-N1; i++){
               Ns.push_back(0);
           }
       }
       if(N2>=N1){
           for(int i=N1; i<=N2min; i++){
               Ns.push_back(i);
           }
           for(int i=N2min+1; i<=L; i++){
                Ns.push_back(0);
           }
       }else{
           for(int i=N1; i>=N2; i--){
               Ns.push_back(i);
           }
       }
    }
    return Ns;
}

TValsShowHide::TValsShowHide(){
    DisableWriting();
    this->ConsoleInterface=true;
    this->f=NULL;
}
TValsShowHide::~TValsShowHide(){
    //if(this->f!=NULL)delete this->f;//no need: V ne fy new oper
}
void TValsShowHide::Assign(const TValsShowHide&obj){
    this->ConsoleInterface=obj.ConsoleInterface;
    this->Show1Hide0=obj.Show1Hide0;
    this->f=obj.f;
    //
}
//void TValsShowHide::Assign(const TValsShowHide&obj){
void TValsShowHide::SetShow1Hide0(int Show1Hide0){
    this->Show1Hide0=Show1Hide0;
}
void TValsShowHide::EnableWriting(){
    this->Show1Hide0=1;
}
void TValsShowHide::DisableWriting(){
    this->Show1Hide0=0;
}




/*
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
*/
