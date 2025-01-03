#include "tests.h"
#include"mylib.h"

//Tests::Tests()
//{
//
//}
void TestQStringFormatting(){
    QString sarr[]={"Chapter One", "Chapter Two", "Chapter Three", "Chapter Four", "Chapter Five", "Chapter Six", "Chapter Seven", "Chapter Eight", "Chapter Nine"};
    QString qs, dfltStr, hilfsStr, qscur="";
    std::string ss;
    int QDefaultBefore, FromN, existingL, requiredL, LreqGiven, QItems=9;
    //
    //QItems=sarr.size();
    //QItems=sarr.length();
    //vrn1
    QDefaultBefore=0;
    LreqGiven=15;
    FromN=1;
    dfltStr=".";
    //int zero=0;
    //qs=MySubSttring_byLs_Formatting_PreMark(sarr[1-1], LreqGiven, FromN, QDefaultBefore, dfltStr);//hin same arr as in cycle
    for(int i=1; i<=QItems; i++){
        //QString MySubSttring_byLs_Formatting_PreMark(QString Sini, int Lreq=0, int FromN=1, int QDefaultBefore=0, QString Default=" ");
        qscur=sarr[i-1];
        //qs=MySubSttring_byLs_Formatting_PreMark(sarr[i-1], LreqGiven, FromN, QDefaultBefore, dfltStr);
        qs=MySubString_byLs_Formatting_PreMark(qscur, LreqGiven, FromN, QDefaultBefore, dfltStr);
        //hilfsStr.setNum(i);
        hilfsStr=IntToStr(i);
        //qs=qs+hilfsStr;
        //ss=qs.toStdString();
        std::cout<<qs.toStdString().c_str()<<hilfsStr.toStdString().c_str()<<std::endl;
    }
    qs="";
    QString stringOfNum="12.34";
    int LsToShow[]={2, 3, 4, 5, 6, 7, 8};
    int NsToShow[]={1, 2, 3, 4, 5, 6, 7};
    int countZerosBefore, LToReq, LNumOwnTmp, LownIni=stringOfNum.length();
    for(int i=1; i<=7; i++){
        for(int j=1; j<=7; j++){
            LNumOwnTmp=LownIni-NsToShow[j-1]+1;
            LToReq=LsToShow[i-1];
            if(LNumOwnTmp<0)LNumOwnTmp=0;
            if(LNumOwnTmp<=LsToShow[i-1]){
                countZerosBefore=LsToShow[i-1]-LNumOwnTmp;
                //LToReq=LsToShow[i-1];
            }else{
                countZerosBefore=0;
                //LToReq=0;
            }
            std::cout<<"12.34, L="<<LsToShow[i-1]<<" N1="<<NsToShow[j-1]<<" LfactToShow="<<LNumOwnTmp<<" Zeros before: "<<countZerosBefore<<" : ";
            qs=MySubString_byLs_Formatting_PreMark(stringOfNum, LsToShow[i-1], NsToShow[j-1], countZerosBefore, "0");
            std::cout<<qs.toStdString().c_str()<<std::endl;
        }
    }
}

void Test1DLib(){
    TValsShowHide vsh;
    vsh.EnableWriting();
    vsh.ConsoleInterface=true;
    //vsh.
    int line1[]={11, 10, 20, 30, 15, 16, 17, 18}, line2[]={10, 20, 30};
    int Q=8;
    int* X=line1;
    int val1=55;
    std::cout<<"Val to add: "<<val1<<std::endl;
    std::cout<<"Array bef add: "<<std::endl;
    for(int i=1; i<=Q; i++){ std::cout<<X[i-1]<<"; "; }
    std::cout<<std::endl;
    Array1DAddElement(X, Q, val1);
    std::cout<<"Array af add: "<<std::endl;
    for(int i=1; i<=Q; i++){
        std::cout<<X[i-1]<<"; ";
        //std::cout<<"X["<<i<<"]="<<X[i-1]<<"; ";
    }
    std::cout<<std::endl;
    std::cout<<"Array aft del last (N="<<Q<<"):"<<std::endl;
    Array1DDelElementFromN(X, Q, Q);
    for(int i=1; i<=Q; i++){
        std::cout<<X[i-1]<<"; ";
        //std::cout<<"X["<<i<<"]="<<X[i-1]<<"; ";
    }
    std::cout<<std::endl;
    std::cout<<"Array aft ins to first (N="<<1<<"):"<<std::endl;
    Array1DInsElementToN(X, Q, 1, val1);
    for(int i=1; i<=Q; i++){
        std::cout<<X[i-1]<<"; ";
        //std::cout<<"X["<<i<<"]="<<X[i-1]<<"; ";
    }
    std::cout<<std::endl;
    std::cout<<"Array aft del first (N="<<1<<"):"<<std::endl;
    Array1DDelElementFromN(X, Q, 1);
    for(int i=1; i<=Q; i++){
        std::cout<<X[i-1]<<"; ";
        //std::cout<<"X["<<i<<"]="<<X[i-1]<<"; ";
    }
    std::cout<<std::endl;
    std::cout<<"Array aft ins to  N="<<3<<":"<<std::endl;
    Array1DInsElementToN(X, Q, 3, val1);
    for(int i=1; i<=Q; i++){
        std::cout<<X[i-1]<<"; ";
        //std::cout<<"X["<<i<<"]="<<X[i-1]<<"; ";
    }
    std::cout<<std::endl;
    std::cout<<"Array aft del from: N="<<3<<":"<<std::endl;
    Array1DDelElementFromN(X, Q, 3);
    for(int i=1; i<=Q; i++){
        std::cout<<X[i-1]<<"; ";
        //std::cout<<"X["<<i<<"]="<<X[i-1]<<"; ";
    }
    std::cout<<std::endl;
    std::cout<<"Array aft swap N1 "<<3<<" and "<<5<<":"<<std::endl;
    Array1DSwapElements(X, Q, 3, 5);
    for(int i=1; i<=Q; i++){
        std::cout<<X[i-1]<<"; ";
        //std::cout<<"X["<<i<<"]="<<X[i-1]<<"; ";
    }
    std::cout<<std::endl;
    std::cout<<"Array aft swap N2 "<<3<<" and "<<5<<":"<<std::endl;
    Array1DSwapElements(X, Q, 3, 5);
    for(int i=1; i<=Q; i++){
        std::cout<<X[i-1]<<"; ";
        //std::cout<<"X["<<i<<"]="<<X[i-1]<<"; ";
    }
    std::cout<<std::endl;
    std::cout<<"Array aft Reverse N1:"<<std::endl;
    Array1DReverse(X, Q);
    for(int i=1; i<=Q; i++){
        std::cout<<X[i-1]<<"; ";
        //std::cout<<"X["<<i<<"]="<<X[i-1]<<"; ";
    }
    std::cout<<std::endl;
    std::cout<<"Array aft Reverse N2:"<<std::endl;
    Array1DReverse(X, Q);
    for(int i=1; i<=Q; i++){
        std::cout<<X[i-1]<<"; ";
        //std::cout<<"X["<<i<<"]="<<X[i-1]<<"; ";
    }
    std::cout<<std::endl;
    std::cout<<"Array aft add :"<<std::endl;
    Array1DAddElement(X, Q, 10);
    Array1DAddElement(X, Q, 20);
    Array1DAddElement(X, Q, 30);
    Array1DAddElement(X, Q, 21);
    Array1DAddElement(X, Q, 10);
    Array1DAddElement(X, Q, 20);
    Array1DAddElement(X, Q, 22);
    for(int i=1; i<=Q; i++){
        std::cout<<X[i-1]<<"; ";
        //std::cout<<"X["<<i<<"]="<<X[i-1]<<"; ";
    }
    std::cout<<"Now Q="<<Q<<std::endl;
    std::vector<int>Ns;
    val1=20;
    std::cout<<"Seek val in whole arr, val="<<val1<<std::endl;
    //template <typename T> std::vector<int>Array1DSeekValSimply(T*data, T val, int curL, int FromN=1, int ToN=0, int param=0){
    Ns=Array1DSeekValSimply(X, val1, Q);
    for(int i=1; i<=Ns.size(); i++){
        std::cout<<Ns[i-1]<<": ";
        for(int j=Ns[i-1]; j<=Ns[i-1]+1-1; j++){
            std::cout<<X[j-1]<<" "<<std::endl;
        }
        //std::cout<<"N["<<i<<"]="<<Ns[i-1]<<"; ";
    }
    std::cout<<std::endl;
    int sQ=3,  FromN=1, ToN=0;
    std::cout<<"Seek sub-arr in whole arr"<<std::endl;
    //template <typename T> std::vector<int>Array1DSeekValSimply(T*data, T val, int curL, int FromN=1, int ToN=0, int param=0){
    Ns=Array1DSeekSubArraySimply(X, line2, Q, 3);
    for(int i=1; i<=Ns.size(); i++){
        std::cout<<Ns[i-1]<<": ";
        for(int j=Ns[i-1]; j<=Ns[i-1]+sQ-1; j++){
            std::cout<<X[j-1]<<" ";
        }
       std::cout<<std::endl;
    }
    std::cout<<std::endl;
    std::cout<<"Seek sub-arr in whole arr (1...0)"<<std::endl;
    //template <typename T> std::vector<int>Array1DSeekValSimply(T*data, T val, int curL, int FromN=1, int ToN=0, int param=0){
    Ns=Array1DSeekSubArraySimply(X, line2, Q, 3, FromN, ToN);
    for(int i=1; i<=Ns.size(); i++){
        std::cout<<Ns[i-1]<<": ";
        for(int j=Ns[i-1]; j<=Ns[i-1]+sQ-1; j++){
            std::cout<<X[j-1]<<" ";
        }
       std::cout<<std::endl;
    }
    std::cout<<std::endl;
    std::cout<<"Seek sub-arr in whole arr (1...0, visa versa order)"<<std::endl;
    FromN=14; ToN=2;
    //template <typename T> std::vector<int>Array1DSeekValSimply(T*data, T val, int curL, int FromN=1, int ToN=0, int param=0){
    Ns=Array1DSeekSubArraySimply(X, line2, Q, 3, FromN, ToN);
    for(int i=1; i<=Ns.size(); i++){
        std::cout<<Ns[i-1]<<": ";
        for(int j=Ns[i-1]; j<=Ns[i-1]+sQ-1; j++){
            std::cout<<X[j-1]<<" ";
        }
       std::cout<<std::endl;
    }
    std::cout<<std::endl;
    FromN=3; ToN=0;
    std::cout<<"Seek sub-arr in whole arr (3...0)-("<<FromN<<"..."<<ToN<<")"<<std::endl;
    //template <typename T> std::vector<int>Array1DSeekValSimply(T*data, T val, int curL, int FromN=1, int ToN=0, int param=0){
    Ns=Array1DSeekSubArraySimply(X, line2, Q, sQ, FromN, ToN);
    for(int i=1; i<=Ns.size(); i++){
        std::cout<<Ns[i-1]<<": ";
        for(int j=Ns[i-1]; j<=Ns[i-1]+sQ-1; j++){
            std::cout<<X[j-1]<<" ";
        }
       std::cout<<std::endl;
    }
    std::cout<<std::endl;
    FromN=2; ToN=7;
    std::cout<<"Seek sub-arr in whole arr (2...7)-("<<FromN<<"..."<<ToN<<")"<<std::endl;
    //template <typename T> std::vector<int>Array1DSeekValSimply(T*data, T val, int curL, int FromN=1, int ToN=0, int param=0){
    Ns=Array1DSeekSubArraySimply(X, line2, Q, sQ, FromN, ToN);
    for(int i=1; i<=Ns.size(); i++){
        std::cout<<Ns[i-1]<<": ";
        for(int j=Ns[i-1]; j<=Ns[i-1]+sQ-1; j++){
            std::cout<<X[j-1]<<" ";
        }
       std::cout<<std::endl;
    }
    std::cout<<std::endl;
    //
    int LToReq, countZerosBefore;
    int *DfltIntVal=new int;
    *DfltIntVal=0;
    //
    FromN=1; ToN=0; LToReq=0; countZerosBefore=0;
    std::cout<<"Getting Sub-Array"<<std::endl;
    //<T>std::vector<T> Array1DGetSubArray_byLs_PreMark(T*x, int Lold, int Lreq=0, int whatFromN=1, int FirstDefaultValues=0, T*DfltValParam=NULL){
    if(LToReq==0)std::cout<<"Final length unlimited"<<std::endl;
    else std::cout<<"Final length="<<LToReq<<std::endl;
    std::cout<<"FromN="<<FromN<<" default vals before: "<<countZerosBefore<<" iniL="<<Q<<std::endl;
    Ns=Array1DGetSubArray_byLs_PreMark(X, Q, LToReq, FromN, countZerosBefore, DfltIntVal);
    for(int i=1; i<=Ns.size(); i++){
        std::cout<<Ns[i-1]<<" ";
    }
    std::cout<<std::endl;
    FromN=3; ToN=0; LToReq=0; countZerosBefore=2;
    std::cout<<"Getting Sub-Array"<<std::endl;
    //<T>std::vector<T> Array1DGetSubArray_byLs_PreMark(T*x, int Lold, int Lreq=0, int whatFromN=1, int FirstDefaultValues=0, T*DfltValParam=NULL){
    if(LToReq==0)std::cout<<"Final length unlimited"<<std::endl;
    else std::cout<<"Final length="<<LToReq<<std::endl;
    std::cout<<"FromN="<<FromN<<" default vals before: "<<countZerosBefore<<" iniL="<<Q<<std::endl;
    Ns=Array1DGetSubArray_byLs_PreMark(X, Q, LToReq, FromN, countZerosBefore, DfltIntVal);
    for(int i=1; i<=Ns.size(); i++){
        std::cout<<Ns[i-1]<<" ";
    }
    std::cout<<std::endl;
    FromN=3; ToN=0; LToReq=8; countZerosBefore=2;
    std::cout<<"Getting Sub-Array"<<std::endl;
    //<T>std::vector<T> Array1DGetSubArray_byLs_PreMark(T*x, int Lold, int Lreq=0, int whatFromN=1, int FirstDefaultValues=0, T*DfltValParam=NULL){
    if(LToReq==0)std::cout<<"Final length unlimited"<<std::endl;
    else std::cout<<"Final length="<<LToReq<<std::endl;
    std::cout<<"FromN="<<FromN<<" default vals before: "<<countZerosBefore<<" iniL="<<Q<<std::endl;
    Ns=Array1DGetSubArray_byLs_PreMark(X, Q, LToReq, FromN, countZerosBefore, DfltIntVal);
    for(int i=1; i<=Ns.size(); i++){
        std::cout<<Ns[i-1]<<" ";
    }
    std::cout<<std::endl;
    FromN=3; ToN=0; LToReq=19; countZerosBefore=2;
    std::cout<<"Getting Sub-Array"<<std::endl;
    //<T>std::vector<T> Array1DGetSubArray_byLs_PreMark(T*x, int Lold, int Lreq=0, int whatFromN=1, int FirstDefaultValues=0, T*DfltValParam=NULL){
    if(LToReq==0)std::cout<<"Final length unlimited"<<std::endl;
    else std::cout<<"Final length="<<LToReq<<std::endl;
    std::cout<<"FromN="<<FromN<<" default vals before: "<<countZerosBefore<<" iniL="<<Q<<std::endl;
    Ns=Array1DGetSubArray_byLs_PreMark(X, Q, LToReq, FromN, countZerosBefore, DfltIntVal);
    for(int i=1; i<=Ns.size(); i++){
        std::cout<<Ns[i-1]<<" ";
    }
     std::cout<<std::endl;
    //
    delete DfltIntVal;
    std::cout<<std::endl;
}//Test1DLib()

void TestArr2DLib(){
    //
    TValsShowHide vsh;
    vsh.EnableWriting();
    //vsh.ConsoleInterface();
    int*DfltInt=new int;
    *DfltInt=0;
    int ExtRowN, IneRowN;
    int *Ls, *L1s;//, *DfltInt;
    int Q=9, LiniExN1=9, LcurRow;// Lini1=9, Lini2=7, Lini3=8, Lini4=8, Lini5=7, Lini6=9, Lini7=6, Lini8=5, Lini9=8;
    Array2DSize size1, size2;
    //
    //----------------------------------------------------------------------------------------------------------------------------------------------
    //
    /*
    int arr1DLine1[]={11, 12, 13, 14, 15, 16, 17, 18, 19};
    int arr1DLine2[]={21, 22, 23, 24, 25, 26, 27, 28, 29};
    int arr1DLine3[]={31, 32, 33, 34, 35, 36, 37, 38, 39};
    int arr1DLine4[]={41, 42, 43, 44, 45, 46, 47, 48, 49};
    int arr1DLine5[]={51, 52, 53, 54, 55, 56, 57, 58, 59};
    int arr1DLine6[]={61, 62, 63, 64, 65, 66, 67, 68, 69};
    int arr1DLine7[]={71, 72, 73, 74, 75, 76, 77, 78, 79};
    int arr1DLine8[]={81, 82, 83, 84, 85, 86, 87, 88, 89};
    int arr1DLine9[]={91, 92, 93, 94, 95, 96, 97, 98, 99};
    */
    //
    /*
    int arr1DLine1[]={11, 12, 13, 14, 15, 16, 17, 18, 19};
    int arr1DLine2[]={21, 100,200,300,25, 26, 27, 28, 29};
    int arr1DLine3[]={31, 400,33, 34, 35, 36, 37, 38, 39};
    int arr1DLine4[]={41, 500,600,44, 45, 46, 47, 48, 49};
    int arr1DLine5[]={51, 52, 53, 54, 55, 56, 57, 58, 59};
    int arr1DLine6[]={61, 62, 63, 64, 65, 66, 67, 68, 69};
    int arr1DLine7[]={71, 72, 73, 74, 75, 76, 77, 78, 79};
    int arr1DLine8[]={81, 82, 83, 84, 85, 86, 87, 88, 89};
    int arr1DLine9[]={91, 92, 93, 94, 95, 96, 97, 98, 99};
    */
    //
    /*
    int arr1DLine1[]={11, 12, 13, 14, 15, 16, 17, 18, 19};
    int arr1DLine2[]={300,200,100,24, 25, 26, 27, 28, 29};
    int arr1DLine3[]={31, 32, 400,34, 35, 36, 37, 38, 39};
    int arr1DLine4[]={41, 600,500,44, 45, 46, 47, 48, 49};
    int arr1DLine5[]={51, 52, 53, 54, 55, 56, 57, 58, 59};
    int arr1DLine6[]={61, 62, 63, 64, 65, 66, 67, 68, 69};
    int arr1DLine7[]={71, 72, 73, 74, 75, 76, 77, 78, 79};
    int arr1DLine8[]={81, 82, 83, 84, 85, 86, 87, 88, 89};
    int arr1DLine9[]={91, 92, 93, 94, 95, 96, 97, 98, 99};
    */
    //
    /*
    int arr1DLine1[]={11, 12, 13, 14, 15, 16, 17, 18, 19};
    int arr1DLine2[]={21, 500,400,100,25, 26, 27, 28, 29};
    int arr1DLine3[]={31, 600,33, 200,35, 36, 37, 38, 39};
    int arr1DLine4[]={41, 42, 43, 300,45, 46, 47, 48, 49};
    int arr1DLine5[]={51, 52, 53, 54, 55, 56, 57, 58, 59};
    int arr1DLine6[]={61, 62, 63, 64, 65, 66, 67, 68, 69};
    int arr1DLine7[]={71, 72, 73, 74, 75, 76, 77, 78, 79};
    int arr1DLine8[]={81, 82, 83, 84, 85, 86, 87, 88, 89};
    int arr1DLine9[]={91, 92, 93, 94, 95, 96, 97, 98, 99};
    */
    //
    /*
    int arr1DLine1[]={11, 12, 13, 14, 15, 16, 17, 18, 19};
    int arr1DLine2[]={21, 22, 23, 100,400,500,27, 28, 29};
    int arr1DLine3[]={31, 32, 33, 200,35, 600,37, 38, 39};
    int arr1DLine4[]={41, 42, 43, 300,45, 46, 47, 48, 49};
    int arr1DLine5[]={51, 52, 53, 54, 55, 56, 57, 58, 59};
    int arr1DLine6[]={61, 62, 63, 64, 65, 66, 67, 68, 69};
    int arr1DLine7[]={71, 72, 73, 74, 75, 76, 77, 78, 79};
    int arr1DLine8[]={81, 82, 83, 84, 85, 86, 87, 88, 89};
    int arr1DLine9[]={91, 92, 93, 94, 95, 96, 97, 98, 99};
    */
    //
    /*
    int arr1DLine1[]={11, 12, 13, 100,15, 16, 17, 18, 19};
    int arr1DLine2[]={21, 22, 23, 200,25, 600,27, 28, 29};
    int arr1DLine3[]={31, 32, 33, 300,400,500,37, 38, 39};
    int arr1DLine4[]={41, 42, 43, 44, 45, 46, 47, 48, 49};
    int arr1DLine5[]={51, 52, 53, 54, 55, 56, 57, 58, 59};
    int arr1DLine6[]={61, 62, 63, 64, 65, 66, 67, 68, 69};
    int arr1DLine7[]={71, 72, 73, 74, 75, 76, 77, 78, 79};
    int arr1DLine8[]={81, 82, 83, 84, 85, 86, 87, 88, 89};
    int arr1DLine9[]={91, 92, 93, 94, 95, 96, 97, 98, 99};
    */
    //
    /*
    int arr1DLine1[]={11, 12, 13, 100,15, 16, 17, 18, 19};
    int arr1DLine2[]={21, 600,23, 200,25, 26, 27, 28, 29};
    int arr1DLine3[]={31, 500,400,300,35, 26, 37, 38, 39};
    int arr1DLine4[]={41, 42, 43, 44, 45, 46, 47, 48, 49};
    int arr1DLine5[]={51, 52, 53, 54, 55, 56, 57, 58, 59};
    int arr1DLine6[]={61, 62, 63, 64, 65, 66, 67, 68, 69};
    int arr1DLine7[]={71, 72, 73, 74, 75, 76, 77, 78, 79};
    int arr1DLine8[]={81, 82, 83, 84, 85, 86, 87, 88, 89};
    int arr1DLine9[]={91, 92, 93, 94, 95, 96, 97, 98, 99};
    */
    //
    /*
    int arr1DLine1[]={11, 12, 13, 14, 15, 16, 17, 18, 19};
    int arr1DLine2[]={21, 22, 23, 24, 25, 26, 27, 28, 29};
    int arr1DLine3[]={31, 32, 33, 34, 35, 36, 37, 38, 39};
    int arr1DLine4[]={41, 42, 43, 44, 45, 46, 47, 48, 49};
    int arr1DLine5[]={51, 52, 53, 54, 55, 56, 57, 58, 59};
    int arr1DLine6[]={61, 62, 63, 64, 65, 66, 67, 68, 69};
    int arr1DLine7[]={71, 72, 73, 74, 75, 76, 77, 78, 79};
    int arr1DLine8[]={81, 82, 83, 84, 85, 86, 87, 88, 89};
    int arr1DLine9[]={91, 92, 93, 94, 95, 96, 97, 98, 99};
    */
    //
    //----
    //
    /*
    int arr1DLine1[]={11, 12, 13, 14, 15, 16, 17, 18, 19};
    int arr1DLine2[]={21, 22, 23, 24, 25, 26, 27};
    int arr1DLine3[]={31, 32, 33, 34, 35, 36, 37, 38};
    int arr1DLine4[]={41, 42, 43, 44, 45, 46, 47, 48};
    int arr1DLine5[]={51, 52, 53, 54, 55, 56, 57};
    int arr1DLine6[]={61, 62, 63, 64, 65, 66, 67, 68, 69};
    int arr1DLine7[]={71, 72, 73, 74, 75, 76};
    int arr1DLine8[]={81, 82, 83, 84, 85};
    int arr1DLine9[]={91, 92, 93, 94, 95, 96, 97, 98};
    */
    //
    //
    //
    //int arr1DLine1[]={11, 12, 13, 14, 15, 16, 17, 18, 19};
    //int arr1DLine2[]={21, 100,200,300,25, 26, 27};
    //int arr1DLine3[]={31, 400,33, 34, 35, 36, 37, 38};
    //int arr1DLine4[]={41, 500,600,44, 45, 46, 47, 48};
    //int arr1DLine5[]={51, 52, 53, 54, 55, 56, 57};
    //int arr1DLine6[]={61, 62, 63, 64, 65, 66, 67, 68, 69};
    //int arr1DLine7[]={71, 72, 73, 74, 75, 76};
    //int arr1DLine8[]={81, 82, 83, 84, 85};
    //int arr1DLine9[]={91, 92, 93, 94, 95, 96, 97, 98};
    //
    Ls=new int[Q];
    Ls[1-1]=9;//Lini1;
    Ls[2-1]=7;//Lini2;
    Ls[3-1]=8;//Lini3;
    Ls[4-1]=8;//Lini4;
    Ls[5-1]=7;//Lini5;
    Ls[6-1]=9;//Lini6;
    Ls[7-1]=8;//=Lini7;
    Ls[8-1]=6;//Lini8;
    Ls[9-1]=8;//Lini9;
    size1.Set(Q, Ls);
    //
    int arr1DLine1[]={11, 12,  13,  14,  15, 16, 17, 18, 19};
    int arr1DLine2[]={21, 100, 200, 300, 25, 26, 27};
    int arr1DLine3[]={31, 400, 33,  300, 200, 100, 37, 38};
    int arr1DLine4[]={41, 500, 600, 44,  45,  400, 47, 48};
    int arr1DLine5[]={51, 300, 53,  54,  600, 500, 57};
    int arr1DLine6[]={61, 200, 63,  600, 100, 400, 500, 68, 69};
    int arr1DLine7[]={71, 100, 400, 500, 200, 76,  600, 78};
    int arr1DLine8[]={81, 82,  83,  84,  300, 86};
    int arr1DLine9[]={91, 92,  93,  94,  95, 96, 97, 98};
    //
    int arr1DN2Line1[]={100, 200, 300};
    int arr1DN2Line2[]={400};
    int arr1DN2Line3[]={500, 600};
    //
    /*
    int arr1DLine1[]={11, 12, 13, 14, 15, 16, 17, 18, 19};
    int arr1DLine2[]={300,200,100,24, 25, 26, 27};
    int arr1DLine3[]={31, 32, 400,34, 35, 36, 37, 38};
    int arr1DLine4[]={41, 600,500,44, 45, 46, 47, 48};
    int arr1DLine5[]={51, 52, 53, 54, 55, 56, 57};
    int arr1DLine6[]={61, 62, 63, 64, 65, 66, 67, 68, 69};
    int arr1DLine7[]={71, 72, 73, 74, 75, 76};
    int arr1DLine8[]={81, 82, 83, 84, 85};
    int arr1DLine9[]={91, 92, 93, 94, 95, 96, 97, 98};
    */
    //
    /*
    int arr1DLine1[]={11, 12, 13, 14, 15, 16, 17, 18, 19};
    int arr1DLine2[]={21, 500,400,100,25, 26, 27};
    int arr1DLine3[]={31, 600,33, 200,35, 36, 37, 38};
    int arr1DLine4[]={41, 42, 43, 300,45, 46, 47, 48};
    int arr1DLine5[]={51, 52, 53, 54, 55, 56, 57};
    int arr1DLine6[]={61, 62, 63, 64, 65, 66, 67, 68, 69};
    int arr1DLine7[]={71, 72, 73, 74, 75, 76};
    int arr1DLine8[]={81, 82, 83, 84, 85};
    int arr1DLine9[]={91, 92, 93, 94, 95, 96, 97, 98};
    */
    //
    /*
    int arr1DLine1[]={11, 12, 13, 14, 15, 16, 17, 18, 19};
    int arr1DLine2[]={21, 22, 23, 100,400,500,27};
    int arr1DLine3[]={31, 32, 33, 200,35, 600,37, 38};
    int arr1DLine4[]={41, 42, 43, 300,45, 46, 47, 48};
    int arr1DLine5[]={51, 52, 53, 54, 55, 56, 57};
    int arr1DLine6[]={61, 62, 63, 64, 65, 66, 67, 68, 69};
    int arr1DLine7[]={71, 72, 73, 74, 75, 76};
    int arr1DLine8[]={81, 82, 83, 84, 85};
    int arr1DLine9[]={91, 92, 93, 94, 95, 96, 97, 98};
    */
    //
    /*
    int arr1DLine1[]={11, 12, 13, 100,15, 16, 17, 18, 19};
    int arr1DLine2[]={21, 22, 23, 200,25, 600,27};
    int arr1DLine3[]={31, 32, 33, 300,400,500,37, 38};
    int arr1DLine4[]={41, 42, 43, 44, 45, 46, 47, 48};
    int arr1DLine5[]={51, 52, 53, 54, 55, 56, 57};
    int arr1DLine6[]={61, 62, 63, 64, 65, 66, 67, 68, 69};
    int arr1DLine7[]={71, 72, 73, 74, 75, 76};
    int arr1DLine8[]={81, 82, 83, 84, 85};
    int arr1DLine9[]={91, 92, 93, 94, 95, 96, 97, 98};
    */
    //
    /*
    int arr1DLine1[]={11, 12, 13, 100,15, 16, 17, 18, 19};
    int arr1DLine2[]={21, 600,23, 200,25, 26, 27};
    int arr1DLine3[]={31, 500,400,300,35, 26, 37, 38};
    int arr1DLine4[]={41, 42, 43, 44, 45, 46, 47, 48};
    int arr1DLine5[]={51, 52, 53, 54, 55, 56, 57};
    int arr1DLine6[]={61, 62, 63, 64, 65, 66, 67, 68, 69};
    int arr1DLine7[]={71, 72, 73, 74, 75, 76};
    int arr1DLine8[]={81, 82, 83, 84, 85};
    int arr1DLine9[]={91, 92, 93, 94, 95, 96, 97, 98};
    */
    //
    /*
    int arr1DLine1[]={11, 12, 13, 14, 15, 16, 17, 18, 19};
    int arr1DLine2[]={21, 22, 23, 24, 25, 26, 27};
    int arr1DLine3[]={31, 32, 33, 34, 35, 36, 37, 38};
    int arr1DLine4[]={41, 42, 43, 44, 45, 46, 47, 48};
    int arr1DLine5[]={51, 52, 53, 54, 55, 56, 57};
    int arr1DLine6[]={61, 62, 63, 64, 65, 66, 67, 68, 69};
    int arr1DLine7[]={71, 72, 73, 74, 75, 76};
    int arr1DLine8[]={81, 82, 83, 84, 85};
    int arr1DLine9[]={91, 92, 93, 94, 95, 96, 97, 98};
    */
    int arr1DN1[]={101, 102, 103, 104, 105, 106, 107, 108, 109};
    int arr1DN2[]={201, 202, 203, 204, 205, 206, 207, 208, 209};
    //
    int *arr1DLine1Ptr=arr1DLine1,
        *arr1DLine2Ptr=arr1DLine2,
        *arr1DLine3Ptr=arr1DLine3,
        *arr1DLine4Ptr=arr1DLine4,
        *arr1DLine5Ptr=arr1DLine5,
        *arr1DLine6Ptr=arr1DLine6,
        *arr1DLine7Ptr=arr1DLine7,
        *arr1DLine8Ptr=arr1DLine8,
        *arr1DLine9Ptr=arr1DLine9,
        //
        *arr1DN2Line1Ptr=arr1DN2Line1,
        *arr1DN2Line2Ptr=arr1DN2Line2,
        *arr1DN2Line3Ptr=arr1DN2Line3,
        //
        *arr1DN1Ptr=arr1DN1,
        *arr1DN2Ptr=arr1DN2;
    int **arr2DN1Ptr, **arr2DN2Ptr;


    //
    arr2DN1Ptr=new int*[Q];
    arr2DN1Ptr[1-1]=arr1DLine1Ptr;
    arr2DN1Ptr[2-1]=arr1DLine2Ptr;
    arr2DN1Ptr[3-1]=arr1DLine3Ptr;
    arr2DN1Ptr[4-1]=arr1DLine4Ptr;
    arr2DN1Ptr[5-1]=arr1DLine5Ptr;
    arr2DN1Ptr[6-1]=arr1DLine6Ptr;
    arr2DN1Ptr[7-1]=arr1DLine7Ptr;
    arr2DN1Ptr[8-1]=arr1DLine8Ptr;
    arr2DN1Ptr[9-1]=arr1DLine9Ptr;

    //
    L1s=new int[3];
    L1s[1-1]=3;
    L1s[2-1]=1;
    L1s[3-1]=2;
    arr2DN2Ptr=new int*[3];
    arr2DN2Ptr[1-1]=arr1DN2Line1Ptr;
    arr2DN2Ptr[2-1]=arr1DN2Line2Ptr;
    arr2DN2Ptr[3-1]=arr1DN2Line3Ptr;
    //
    size2.Set(3, L1s);
    //
    std::cout<<"Ini 2D array:"<<std::endl;
    Arr2DStdCOut1D(arr2DN1Ptr, size1);
    std::cout<<std::endl;
    std::cout<<std::endl;
    Arr2DStdCOut1D(arr2DN1Ptr, size1, " ", "; ", true);
    std::cout<<std::endl;
    std::cout<<std::endl;
    Arr2DStdCOut1D(arr2DN1Ptr, size1, " ", "", true);
    std::cout<<std::endl;
    std::cout<<std::endl;
    Arr2DStdCOut1D(arr2DN1Ptr, size1, " ", "; ", false);
    std::cout<<std::endl;
    std::cout<<std::endl;
    Arr2DStdCOut1D(arr2DN1Ptr, size1, ", ", "; ", false);
    std::cout<<std::endl;
    std::cout<<std::endl;
    Arr2DStdCOut2D(arr2DN1Ptr, size1);
    std::cout<<std::endl;
    std::cout<<std::endl;
    Arr2DStdCOut2D(arr2DN1Ptr, size1, " ", "; ", true);
    std::cout<<std::endl;
    Arr2DStdCOut2D(arr2DN1Ptr, size1, " ", "", true);
    std::cout<<std::endl;
    Arr2DStdCOut2D(arr2DN1Ptr, size1, " ", "; ", false);
    std::cout<<std::endl;
    //
    bool  IsTransposed=false, EachExtRowIsReversed=false;
    std::vector<std::vector<int>> Ns;
    //
    std::cout<<" Testing search ";
    //
    IsTransposed=false;
    EachExtRowIsReversed=false;
    if(IsTransposed)std::cout<<" Row to seek IS transposed ";
    else std::cout<<" Row to seek is NOT transposed ";
    std::cout<<";";
    if(EachExtRowIsReversed)std::cout<<" Each Ext Row IS Reversed ";
    else std::cout<<" Each Ext Row is NOT Reversed ";
    std::cout<<std::endl;
    Ns= Array2D_SeekArrs2ndIn1st(arr2DN1Ptr, size1, arr2DN2Ptr, size2, IsTransposed, EachExtRowIsReversed);
    std::cout<<"Found: "<<Ns.size()<<std::endl;
    for(int i=1; i<=Ns.size(); i++){
        std::cout<<"ExtRowN: "<<Ns[i-1][1-1]<<" IneRowN: "<<Ns[i-1][2-1]<<std::endl;
    }
    std::cout<<"Must be correct answer: Found 1 ExtRowN: 2 IneRowN: 2"<<std::endl;
    std::cout<<std::endl;
    //
    IsTransposed=false;
    EachExtRowIsReversed=true;
    if(IsTransposed)std::cout<<" Row to seek IS transposed ";
    else std::cout<<" Row to seek is NOT transposed ";
    std::cout<<";";
    if(EachExtRowIsReversed)std::cout<<" Each Ext Row IS Reversed ";
    else std::cout<<" Each Ext Row is NOT Reversed ";
    std::cout<<std::endl;
    Ns= Array2D_SeekArrs2ndIn1st(arr2DN1Ptr, size1, arr2DN2Ptr, size2, IsTransposed, EachExtRowIsReversed);
    std::cout<<"Found: "<<Ns.size()<<std::endl;
    for(int i=1; i<=Ns.size(); i++){
        std::cout<<"ExtRowN: "<<Ns[i-1][1-1]<<" IneRowN: "<<Ns[i-1][2-1]<<std::endl;
    }
    std::cout<<"Must be correct answer: Found 1 ExtRowN: 3 IneRowN: 6"<<std::endl;
    std::cout<<std::endl;
    //
    IsTransposed=true;
    EachExtRowIsReversed=false;
    if(IsTransposed)std::cout<<" Row to seek IS transposed ";
    else std::cout<<" Row to seek is NOT transposed ";
    std::cout<<";";
    if(EachExtRowIsReversed)std::cout<<" Each Ext Row IS Reversed ";
    else std::cout<<" Each Ext Row is NOT Reversed ";
    std::cout<<std::endl;
    Ns= Array2D_SeekArrs2ndIn1st(arr2DN1Ptr, size1, arr2DN2Ptr, size2, IsTransposed, EachExtRowIsReversed);
    std::cout<<"Found: "<<Ns.size()<<std::endl;
    for(int i=1; i<=Ns.size(); i++){
        std::cout<<"ExtRowN: "<<Ns[i-1][1-1]<<" IneRowN: "<<Ns[i-1][2-1]<<std::endl;
    }
    std::cout<<"Must be correct answer: Found 6 ExtRowN: 6 IneRowN: 5"<<std::endl;
    std::cout<<std::endl;
    //
    IsTransposed=true;
    EachExtRowIsReversed=true;
    if(IsTransposed)std::cout<<" Row to seek IS transposed ";
    else std::cout<<" Row to seek is NOT transposed ";
    std::cout<<";";
    if(EachExtRowIsReversed)std::cout<<" Each Ext Row IS Reversed ";
    else std::cout<<" Each Ext Row is NOT Reversed ";
    std::cout<<std::endl;
    Ns= Array2D_SeekArrs2ndIn1st(arr2DN1Ptr, size1, arr2DN2Ptr, size2, IsTransposed, EachExtRowIsReversed);
    std::cout<<"Found: "<<Ns.size()<<std::endl;
    for(int i=1; i<=Ns.size(); i++){
        std::cout<<"ExtRowN: "<<Ns[i-1][1-1]<<" IneRowN: "<<Ns[i-1][2-1]<<std::endl;
    }
    std::cout<<"Must be correct answer: Found 1 ExtRowN: 7 IneRowN: 2"<<std::endl;
    std::cout<<std::endl;
    //
    std::cout<<"Press any key"<<std::endl;
    getch();
    //
    std::cout<<"Transposing:"<<std::endl;
    //
    Array2DTranspose(arr2DN1Ptr, size1);
    Arr2DStdCOut2D(arr2DN1Ptr, size1, " ", "; ", false);
    std::cout<<std::endl;
    std::cout<<"Size now: "<<size1.GetAsStdString().c_str()<<std::endl;
    Q=size1.GetQExtRows();
    std::cout<<"Transposing again:"<<std::endl;
    Array2DTranspose(arr2DN1Ptr, size1);
    Arr2DStdCOut2D(arr2DN1Ptr, size1, " ", "; ", false);
    std::cout<<std::endl;
    std::cout<<"Size now: "<<size1.GetAsStdString().c_str()<<std::endl;
    Q=size1.GetQExtRows();
    //
    //int Lcur;
    std::cout<<"Deleting from ext rows:"<<std::endl;
    for(int i=2; i<=Q; i++){
        for(int j=Q; j>=Q+1-i+1; j--){
            //Lcur=size
            Array2DDelFromExtRowN(arr2DN1Ptr, size1, i, j);//size1.GetMaxLength()-j+1);
        }
    }
    Arr2DStdCOut2D(arr2DN1Ptr, size1, " ", "; ", false);
    std::cout<<std::endl;
    std::cout<<"Size now: "<<size1.GetAsStdString().c_str()<<std::endl;
    Q=size1.GetQExtRows();
    //
    //
    std::cout<<"Adding line:" <<std::endl;
    //template<typename T>void Array2DAddExtRow(T**&X, Array2DSize&size, T*rowParam, int wholeRowL=0, int whatFromN=1, int QDefaultValsBefore=0, bool RectNotVar=true, T*DfltValParam=NULL, TValsShowHide*vsh=NULL){
    Array2DAddExtRow(arr2DN1Ptr, size1, arr1DN1Ptr, LiniExN1);
    Arr2DStdCOut2D(arr2DN1Ptr, size1, " ", "; ", false);
    std::cout<<std::endl;
    Q=size1.GetQExtRows();
    //LcurRow=9;
    //Arr1DStdCOut(arr2DN1Ptr[Q-1], LcurRow);
    std::cout<<std::endl;
    //
    std::cout<<"Deleting line from "<<Q<<" position :"<<std::endl;
    Array2DDelExtRow(arr2DN1Ptr, size1, Q);
    Arr2DStdCOut2D(arr2DN1Ptr, size1, " ", "; ", false);
    std::cout<<std::endl;
    std::cout<<"Size now: "<<size1.GetAsStdString().c_str()<<std::endl;
    Q=size1.GetQExtRows();
    //
    ExtRowN=3;
    std::cout<<"Inserting line to "<<ExtRowN<<" position :"<<std::endl;
    //Array2DInsExtRow(arr2DN1Ptr, size1, ExtRowN, 9, arr1DN1Ptr, DfltInt, 1, 0, true, &vsh);
    Array2DInsExtRow(arr2DN1Ptr, size1, ExtRowN, 9, arr1DN1Ptr, DfltInt, 1, 0, true);
    Arr2DStdCOut2D(arr2DN1Ptr, size1, " ", "; ", false);
    std::cout<<std::endl;
    Q=size1.GetQExtRows();
    //
    std::cout<<"Deleting line from "<<ExtRowN<<" position :"<<std::endl;
    Array2DDelExtRow(arr2DN1Ptr, size1, ExtRowN);
    Arr2DStdCOut2D(arr2DN1Ptr, size1, " ", "; ", false);
    std::cout<<std::endl;
    std::cout<<"Size now: "<<size1.GetAsStdString().c_str()<<std::endl;
    Q=size1.GetQExtRows();
    //
    ExtRowN=1;
    std::cout<<"Inserting line to "<<ExtRowN<<" position :"<<std::endl;
    //Array2DInsExtRow(arr2DN1Ptr, size1, ExtRowN, 9, arr1DN1Ptr, DfltInt, 1, 0, true, &vsh);
    Array2DInsExtRow(arr2DN1Ptr, size1, ExtRowN, 9, arr1DN1Ptr, DfltInt, 1, 0, true);
    Arr2DStdCOut2D(arr2DN1Ptr, size1, " ", "; ", false);
    std::cout<<std::endl;
    Q=size1.GetQExtRows();
    //
    std::cout<<"Setting line at "<<ExtRowN<<" position :"<<std::endl;
    //void Array2DSetExtRowN(T**&X, Array2DSize&size, int ExtRowN, int wholeRowL=0, T*rowParam=NULL, int whatFromN=1, int QDefaultValuesBefore=0, bool KeepIfRect=true, T*DefaultValParam=NULL, bool DefaultIsSpecNotOwn=false){
    Array2DSetExtRowN(arr2DN1Ptr, size1, ExtRowN, 9, arr1DN2Ptr);//, int whatFromN=1, int QDefaultValuesBefore=0, bool KeepIfRect=true, T*DefaultValParam=NULL, bool DefaultIsSpecNotOwn=false){
    Arr2DStdCOut2D(arr2DN1Ptr, size1, " ", "; ", false);
    std::cout<<std::endl;
    Q=size1.GetQExtRows();
    std::cout<<"Press any key:"<<std::endl;
    getch();
    //
    std::cout<<"Deleting line from "<<ExtRowN<<" position :"<<std::endl;
    Array2DDelExtRow(arr2DN1Ptr, size1, ExtRowN);
    Arr2DStdCOut2D(arr2DN1Ptr, size1, " ", "; ", false);
    std::cout<<std::endl;
    std::cout<<"Size now: "<<size1.GetAsStdString().c_str()<<std::endl;
    Q=size1.GetQExtRows();
    //
    ExtRowN=Q;
    std::cout<<"Inserting line to (last) position :"<<std::endl;
    //Array2DInsExtRow(arr2DN1Ptr, size1, ExtRowN, 9, arr1DN1Ptr, DfltInt, 1, 0, true, &vsh);
    Array2DInsExtRow(arr2DN1Ptr, size1, ExtRowN, 9, arr1DN1Ptr, DfltInt, 1, 0, true);
    Arr2DStdCOut2D(arr2DN1Ptr, size1, " ", "; ", false);
    std::cout<<std::endl;
    Q=size1.GetQExtRows();
    //
    std::cout<<"Deleting line from "<<ExtRowN<<" position :"<<std::endl;
    Array2DDelExtRow(arr2DN1Ptr, size1, ExtRowN);
    Arr2DStdCOut2D(arr2DN1Ptr, size1, " ", "; ", false);
    std::cout<<std::endl;
    std::cout<<"Size now: "<<size1.GetAsStdString().c_str()<<std::endl;
    Q=size1.GetQExtRows();
    //
    // Add, Ins, Del Swap ext rows arb, et se arb
    //
    std::cout<<"Adding column (defaults):"<<std::endl;
    //template<typename T>void Array2DAddExtRow(T**&X, Array2DSize&size, T*rowParam, int wholeRowL=0, int whatFromN=1, int QDefaultValsBefore=0, bool RectNotVar=true, T*DfltValParam=NULL, TValsShowHide*vsh=NULL){
    Array2DAddIneRow(arr2DN1Ptr, size1, arr1DN1Ptr);
    Arr2DStdCOut2D(arr2DN1Ptr, size1, " ", "; ", false);
    std::cout<<std::endl;
    Q=size1.GetQExtRows();
    //LcurRow=9;
    //Arr1DStdCOut(arr2DN1Ptr[Q-1], LcurRow);
    std::cout<<"defaults mean do nil if varia L"<<std::endl;
    std::cout<<std::endl;
    //
    std::cout<<"Adding column (do nil if varia L):"<<std::endl;
    //template<typename T>void Array2DAddExtRow(T**&X, Array2DSize&size, T*rowParam, int wholeRowL=0, int whatFromN=1, int QDefaultValsBefore=0, bool RectNotVar=true, T*DfltValParam=NULL, TValsShowHide*vsh=NULL){
    Array2DAddIneRow(arr2DN1Ptr, size1, arr1DN1Ptr, 9, 0);
    Arr2DStdCOut2D(arr2DN1Ptr, size1, " ", "; ", false);
    std::cout<<std::endl;
    Q=size1.GetQExtRows();
    //LcurRow=9;
    //Arr1DStdCOut(arr2DN1Ptr[Q-1], LcurRow);
    std::cout<<std::endl;
    //
    std::cout<<"Adding column (add to last pos of varia L lines):"<<std::endl;
    //template<typename T>void Array2DAddExtRow(T**&X, Array2DSize&size, T*rowParam, int wholeRowL=0, int whatFromN=1, int QDefaultValsBefore=0, bool RectNotVar=true, T*DfltValParam=NULL, TValsShowHide*vsh=NULL){
    Array2DAddIneRow(arr2DN1Ptr, size1, arr1DN1Ptr, 9, 1);
    Arr2DStdCOut2D(arr2DN1Ptr, size1, " ", "; ", false);
    std::cout<<std::endl;
    Q=size1.GetQExtRows();
    //LcurRow=9;
    //Arr1DStdCOut(arr2DN1Ptr[Q-1], LcurRow);
    std::cout<<"Size now: "<<size1.GetAsStdString().c_str()<<std::endl;
    std::cout<<std::endl;
    //
    std::cout<<"Deleting line from last position of each var L line):"<<std::endl;
    Array2DDelIneRow(arr2DN1Ptr, size1, 0, true);
    Arr2DStdCOut2D(arr2DN1Ptr, size1, " ", "; ", false);
    std::cout<<std::endl;
    std::cout<<"Size now: "<<size1.GetAsStdString().c_str()<<std::endl;
    Q=size1.GetQExtRows();
    //
    std::cout<<"Adding column (stretch to Lmax if varia Ls):"<<std::endl;
    //template<typename T>void Array2DAddExtRow(T**&X, Array2DSize&size, T*rowParam, int wholeRowL=0, int whatFromN=1, int QDefaultValsBefore=0, bool RectNotVar=true, T*DfltValParam=NULL, TValsShowHide*vsh=NULL){
    //Array2DAddIneRow(arr2DN1Ptr, size1, arr1DN1Ptr, 9, 2);
    Array2DAddIneRow(arr2DN1Ptr, size1, arr1DN1Ptr, 9, 2, 1, 0, DfltInt);
    Arr2DStdCOut2D(arr2DN1Ptr, size1, " ", "; ", false);
    std::cout<<std::endl;
    Q=size1.GetQExtRows();
    //LcurRow=9;
    //Arr1DStdCOut(arr2DN1Ptr[Q-1], LcurRow);
    std::cout<<"Size now: "<<size1.GetAsStdString().c_str()<<std::endl;
    std::cout<<std::endl;
    //
    std::cout<<"Setting column (N Lmax):"<<std::endl;
    //template<typename T>Array2DSetIneRow(T**&X, Array2DSize& size, int IneRowN, T*rowParam=NULL, int wholeRowL=0, int whatFromN=1, int QDefaultValsBefore=0, T*DefaultValParam=NULL, bool ExistingInsteadOfDefault=true, bool ExistingOnlyForGTLminNotIgnore=false, bool LastPossIfNotRectAtPos0NotIgnore=false){
    Array2DSetIneRow(arr2DN1Ptr, size1, size1.GetMaxLength(), arr1DN2Ptr);
    Arr2DStdCOut2D(arr2DN1Ptr, size1, " ", "; ", false);
    std::cout<<std::endl;
    Q=size1.GetQExtRows();
    std::cout<<"Size now: "<<size1.GetAsStdString().c_str()<<std::endl;
    std::cout<<std::endl;
    //
    std::cout<<"Deleting line from "<<size1.GetMaxLength()<<" position :"<<std::endl;
    Array2DDelIneRow(arr2DN1Ptr, size1, size1.GetMaxLength());
    Arr2DStdCOut2D(arr2DN1Ptr, size1, " ", "; ", false);
    std::cout<<std::endl;
    std::cout<<"Size now: "<<size1.GetAsStdString().c_str()<<std::endl;
    Q=size1.GetQExtRows();
    ////tried ady del lmax-1 - arb gut. Nu testing ins
    //Array2DInsIneRow(T**&X, Array2DSize& size, int IneRowPosN, T*rowParam=NULL, int wholeRowL=0, int IfAfterLminIgnore0Stretch2=0, int whatFromN=1, int QDefaultValsBefore=0, T*DefaultValParam=NULL)
    //
    IneRowN=1;
    std::cout<<"Inserting col to "<<IneRowN<<" position :"<<std::endl;
    Array2DInsIneRow(arr2DN1Ptr, size1, IneRowN, arr1DN1Ptr, 9, 2);//, int whatFromN=1, int QDefaultValsBefore=0, T*DefaultValParam=NULL)
    Arr2DStdCOut2D(arr2DN1Ptr, size1, " ", "; ", false);
    std::cout<<std::endl;
    std::cout<<"Size now: "<<size1.GetAsStdString().c_str()<<std::endl;
    //arr1DN1Ptr
    //
    std::cout<<"Deleting line from "<<IneRowN<<" position :"<<std::endl;
    Array2DDelIneRow(arr2DN1Ptr, size1, IneRowN);
    Arr2DStdCOut2D(arr2DN1Ptr, size1, " ", "; ", false);
    std::cout<<std::endl;
    std::cout<<"Size now: "<<size1.GetAsStdString().c_str()<<std::endl;
    Q=size1.GetQExtRows();
    ////tried ady del lmax-1 - arb gut. Nu testing ins
    //Array2DInsIneRow(T**&X, Array2DSize& size, int IneRowPosN, T*rowParam=NULL, int wholeRowL=0, int IfAfterLminIgnore0Stretch2=0, int whatFromN=1, int QDefaultValsBefore=0, T*DefaultValParam=NULL)
    //
    IneRowN=2;
    std::cout<<"Inserting col to "<<IneRowN<<" position :"<<std::endl;
    Array2DInsIneRow(arr2DN1Ptr, size1, IneRowN, arr1DN1Ptr, 9, 2);//, int whatFromN=1, int QDefaultValsBefore=0, T*DefaultValParam=NULL)
    Arr2DStdCOut2D(arr2DN1Ptr, size1, " ", "; ", false);
    std::cout<<std::endl;
    std::cout<<"Size now: "<<size1.GetAsStdString().c_str()<<std::endl;
    //arr1DN1Ptr
    //
    std::cout<<"Deleting line from "<<IneRowN<<" position :"<<std::endl;
    Array2DDelIneRow(arr2DN1Ptr, size1, IneRowN);
    Arr2DStdCOut2D(arr2DN1Ptr, size1, " ", "; ", false);
    std::cout<<std::endl;
    std::cout<<"Size now: "<<size1.GetAsStdString().c_str()<<std::endl;
    Q=size1.GetQExtRows();
    ////tried ady del lmax-1 - arb gut. Nu testing ins
    //Array2DInsIneRow(T**&X, Array2DSize& size, int IneRowPosN, T*rowParam=NULL, int wholeRowL=0, int IfAfterLminIgnore0Stretch2=0, int whatFromN=1, int QDefaultValsBefore=0, T*DefaultValParam=NULL)
    //
    IneRowN=size1.GetMaxLength();
    std::cout<<"Inserting col to "<<IneRowN<<" position :"<<std::endl;
    Array2DInsIneRow(arr2DN1Ptr, size1, IneRowN, arr1DN1Ptr, 9, 2);//, int whatFromN=1, int QDefaultValsBefore=0, T*DefaultValParam=NULL)
    Arr2DStdCOut2D(arr2DN1Ptr, size1, " ", "; ", false);
    std::cout<<std::endl;
    std::cout<<"Size now: "<<size1.GetAsStdString().c_str()<<std::endl;
    //arr1DN1Ptr
    //
    std::cout<<"Deleting line from "<<IneRowN<<" position :"<<std::endl;
    Array2DDelIneRow(arr2DN1Ptr, size1, IneRowN);
    Arr2DStdCOut2D(arr2DN1Ptr, size1, " ", "; ", false);
    std::cout<<std::endl;
    std::cout<<"Size now: "<<size1.GetAsStdString().c_str()<<std::endl;
    Q=size1.GetQExtRows();
    ////tried ady del lmax-1 - arb gut. Nu testing ins
    //Array2DInsIneRow(T**&X, Array2DSize& size, int IneRowPosN, T*rowParam=NULL, int wholeRowL=0, int IfAfterLminIgnore0Stretch2=0, int whatFromN=1, int QDefaultValsBefore=0, T*DefaultValParam=NULL)
    //
    // add, del, ins col arb gut, et set arb
    //
    //
    //Search, Concat
    // std::vector<std::vector<int>> Array2D_SeekArrs2ndIn1st(T**Where, const Array2DSize& sizeWhere, T**What, const Array2DSize& sizeWhat, bool IsTransposed=false, bool EachExtRowIsReversed=false){
    //
    std::cout<<" Testing search ";
    //
    IsTransposed=false;
    EachExtRowIsReversed=false;
    if(IsTransposed)std::cout<<" Row to seek IS transposed ";
    else std::cout<<" Row to seek is NOT transposed ";
    std::cout<<";";
    if(EachExtRowIsReversed)std::cout<<" Each Ext Row IS Reversed ";
    else std::cout<<" Each Ext Row is NOT Reversed ";
    std::cout<<std::endl;
    Ns= Array2D_SeekArrs2ndIn1st(arr2DN1Ptr, size1, arr2DN2Ptr, size2, IsTransposed, EachExtRowIsReversed);
    std::cout<<"Found: "<<Ns.size()<<std::endl;
    for(int i=1; i<=Ns.size(); i++){
        std::cout<<"ExtRowN: "<<Ns[i-1][1-1]<<" IneRowN: "<<Ns[i-1][2-1]<<std::endl;
    }
    std::cout<<"Must be correct answer: Found 1 ExtRowN: 2 IneRowN: 2"<<std::endl;
    std::cout<<std::endl;
    //
    IsTransposed=false;
    EachExtRowIsReversed=true;
    if(IsTransposed)std::cout<<" Row to seek IS transposed ";
    else std::cout<<" Row to seek is NOT transposed ";
    std::cout<<";";
    if(EachExtRowIsReversed)std::cout<<" Each Ext Row IS Reversed ";
    else std::cout<<" Each Ext Row is NOT Reversed ";
    std::cout<<std::endl;
    Ns= Array2D_SeekArrs2ndIn1st(arr2DN1Ptr, size1, arr2DN2Ptr, size2, IsTransposed, EachExtRowIsReversed);
    std::cout<<"Found: "<<Ns.size()<<std::endl;
    for(int i=1; i<=Ns.size(); i++){
        std::cout<<"ExtRowN: "<<Ns[i-1][1-1]<<" IneRowN: "<<Ns[i-1][2-1]<<std::endl;
    }
    std::cout<<"Must be correct answer: Found 1 ExtRowN: 3 IneRowN: 6"<<std::endl;
    std::cout<<std::endl;
    //
    IsTransposed=true;
    EachExtRowIsReversed=false;
    if(IsTransposed)std::cout<<" Row to seek IS transposed ";
    else std::cout<<" Row to seek is NOT transposed ";
    std::cout<<";";
    if(EachExtRowIsReversed)std::cout<<" Each Ext Row IS Reversed ";
    else std::cout<<" Each Ext Row is NOT Reversed ";
    std::cout<<std::endl;
    Ns= Array2D_SeekArrs2ndIn1st(arr2DN1Ptr, size1, arr2DN2Ptr, size2, IsTransposed, EachExtRowIsReversed);
    std::cout<<"Found: "<<Ns.size()<<std::endl;
    for(int i=1; i<=Ns.size(); i++){
        std::cout<<"ExtRowN: "<<Ns[i-1][1-1]<<" IneRowN: "<<Ns[i-1][2-1]<<std::endl;
    }
    std::cout<<"Must be correct answer: Found 6 ExtRowN: 6 IneRowN: 5"<<std::endl;
    std::cout<<std::endl;
    //
    IsTransposed=true;
    EachExtRowIsReversed=true;
    if(IsTransposed)std::cout<<" Row to seek IS transposed ";
    else std::cout<<" Row to seek is NOT transposed ";
    std::cout<<";";
    if(EachExtRowIsReversed)std::cout<<" Each Ext Row IS Reversed ";
    else std::cout<<" Each Ext Row is NOT Reversed ";
    std::cout<<std::endl;
    Ns = Array2D_SeekArrs2ndIn1st(arr2DN1Ptr, size1, arr2DN2Ptr, size2, IsTransposed, EachExtRowIsReversed);
    std::cout<<"Found: "<<Ns.size()<<std::endl;
    for(int i=1; i<=Ns.size(); i++){
        std::cout<<"ExtRowN: "<<Ns[i-1][1-1]<<" IneRowN: "<<Ns[i-1][2-1]<<std::endl;
    }
    std::cout<<"Must be correct answer: Found 1 ExtRowN: 7 IneRowN: 2"<<std::endl;
    std::cout<<std::endl;
    //
    std::cout<<"Press any key"<<std::endl;
    getch();
    //
    std::cout<<std::endl;
    //
    std::cout<<"Transposing:"<<std::endl;
    //
    Array2DTranspose(arr2DN1Ptr, size1);
    Arr2DStdCOut2D(arr2DN1Ptr, size1, " ", "; ", false);
    std::cout<<std::endl;
    std::cout<<"Size now: "<<size1.GetAsStdString().c_str()<<std::endl;
    Q=size1.GetQExtRows();
    std::cout<<"Transposing again:"<<std::endl;
    Array2DTranspose(arr2DN1Ptr, size1);
    Arr2DStdCOut2D(arr2DN1Ptr, size1, " ", "; ", false);
    std::cout<<std::endl;
    std::cout<<"Size now: "<<size1.GetAsStdString().c_str()<<std::endl;
    Q=size1.GetQExtRows();
    //
    //
    std::cout<<"Press any key"<<std::endl;
    getch();
    //
    //
    std::cout<<"Testing assign and concat"<<std::endl;
    //
    Array2DAssign(arr2DN2Ptr, size2, arr2DN1Ptr, size1);
    //
    Q=size1.GetQExtRows();
    std::cout<<"Copied row:"<<std::endl;
    Arr2DStdCOut2D(arr2DN2Ptr, size2, " ", "; ", false);
    std::cout<<std::endl;
    //
    int shiftQ=3;//3;//0;//works well by 0 and 3 and rect=true, ma exda
    bool rect=true;
    int QRec=0; //Q records to add ab added array
    int FromN=3; //Erst added N ab added array
    std::cout<<"Adding 2D arr:"<<std::endl;
    //
    std::cout<<"Now let us try concat adding..."<<std::endl;
    Array2DConcatAdding(arr2DN1Ptr, size1, arr2DN2Ptr, size2, shiftQ, DfltInt, rect, QRec, FromN);
    std::cout<<"Concatenated 2D arr (adding):"<<std::endl;
    //Arr2DStdCOut2D(arr2DN2Ptr, size2, " ", "; ", false);
    Arr2DStdCOut2D(arr2DN1Ptr, size1, " ", "; ", false);
    std::cout<<std::endl;
    //works well  - tested for rect arrays,
    //rect for em n'vikts, shiftQ = 0 et 3, FromN = 1 et 3, QRec=0 et 2 //wa errs, fa corrs
    //
    bool Stretching_ifFirstNotRect_AddToCurNotStretchToMax=false;
    shiftQ=-3;//-3;//3
    FromN=1;
    //t T void Array2DConcatStretching(T**&X1, Array2DSize& size1, T**X2, const Array2DSize& size2, int shiftVal=0, T*DefaultValParam=NULL, bool Rect=true, bool IfFirstNotRect_AddToLastNotSrtetchToMax=true, int whatFromN=1, int QRec=0){
    Array2DConcatStretching(arr2DN1Ptr, size1, arr2DN2Ptr, size2, shiftQ, DfltInt, rect, Stretching_ifFirstNotRect_AddToCurNotStretchToMax, FromN, QRec, &vsh);
    std::cout<<"Concatenated 2D arr (stretching):"<<std::endl;
    //Arr2DStdCOut2D(arr2DN2Ptr, size2, " ", "; ", false);
    Arr2DStdCOut2D(arr2DN1Ptr, size1, " ", "; ", false);
    std::cout<<std::endl;
    std::cout<<"trying to repeat..."<<std::endl;
    for(int i=1; i<=10; i++){
        Arr2DStdCOut2D(arr2DN1Ptr, size1, " ", "; ", false);
        std::cout<<std::endl;
    }
    std::cout<<"let us go further..."<<std::endl;
    //works well for rect arrays FromN=3 and 1, ShiftQ=3, 0, -3, QRecs=0 and 7  (which is same) => last row whole of X2
    //and works well if nothing else
    //

    std::cout<<"Now let us try concat adding..."<<std::endl;
    Array2DConcatAdding(arr2DN1Ptr, size1, arr2DN2Ptr, size2, shiftQ, DfltInt, rect, QRec, FromN);
    std::cout<<"Concatenated 2D arr (adding):"<<std::endl;
    //Arr2DStdCOut2D(arr2DN2Ptr, size2, " ", "; ", false);
    Arr2DStdCOut2D(arr2DN1Ptr, size1, " ", "; ", false);
    std::cout<<std::endl;
    //works well  - tested for rect arrays,
    //rect for em n'vikts, shiftQ = 0 et 3, FromN = 1 et 3, QRec=0 et 2 //wa errs, fa corrs
    //

    /*
    //
    std::cout<<"Cutting to Rect:"<<std::endl;
    Array2DSetRect(arr2DN1Ptr, size1, 4);
    Arr2DStdCOut2D(arr2DN1Ptr, size1, " ", "; ", false);
    std::cout<<std::endl;
    std::cout<<"Stretching to Rect:"<<std::endl;
    //template<typename T> void Array2DSetRect(T**&X, Array2DSize&size, int QIneRows, int QExtRows=0, T*DefaultValParam=NULL, bool PreserveVals=true){
    //Array2DSetRect(arr2DN1Ptr, size1, 8);
    Array2DSetRect(arr2DN1Ptr, size1, 8, 0, DfltInt, true);
    Arr2DStdCOut2D(arr2DN1Ptr, size1, " ", "; ", false);
    std::cout<<std::endl;
    //works well  - tested for rect arrays
    */
    //
    //-----------------------
    //Fin. Collecting Garbage
    std::cout<<"Fin. Collecting Garbage"<<std::endl;
    //
    delete[]Ls;
    delete[]L1s;
    //
    if(arr2DN1Ptr!=NULL){
        Q=size1.GetQExtRows();
        for(int i=1; i<=Q; i++){
            if(arr2DN1Ptr[i-1]!=NULL){
                delete[]arr2DN1Ptr[i-1];
            }
        }
    }
    delete[]arr2DN1Ptr;
    //
    for(int i=1; i<=3; i++){
        if(arr2DN2Ptr[i-1]!=NULL){
            delete[]arr2DN2Ptr[i-1];
        }
    }
    delete[]arr2DN2Ptr;
    //
    delete DfltInt;
    //
}//TestArr2DLib
//-------------------------------------------------------------------------------------------------------------
void PtrTest(float*X, int Q){
    for(int i=1; i<=Q; i++){
        std::cout<<X[i-1]<<std::endl;
    }
}
void TestPtr1(){
    float X1[]={1.1, 1.3, 1.3};
    std::vector<float>X2;
    X2.push_back(2.1);
    X2.push_back(2.2);
    X2.push_back(2.3);
    PtrTest(X1, 3);
    PtrTest(X2.data(), X2.size());
}
//-----------------------------------------------------------------------------------------------------------
void TestArr2DLib_VB(){
    //
    TValsShowHide vsh;
    vsh.EnableWriting();
    //vsh.ConsoleInterface();
    int*DfltInt=new int;
    *DfltInt=0;
    int ExtRowN, IneRowN;
    int *Ls, *L1s;//, *DfltInt;
    int Q=9, LiniExN1=9, LcurRow;// Lini1=9, Lini2=7, Lini3=8, Lini4=8, Lini5=7, Lini6=9, Lini7=6, Lini8=5, Lini9=8;
    Array2DSize size1, size2;
    //
    //----------------------------------------------------------------------------------------------------------------------------------------------
    //
    /*
    int arr1DLine1[]={11, 12, 13, 14, 15, 16, 17, 18, 19};
    int arr1DLine2[]={21, 22, 23, 24, 25, 26, 27, 28, 29};
    int arr1DLine3[]={31, 32, 33, 34, 35, 36, 37, 38, 39};
    int arr1DLine4[]={41, 42, 43, 44, 45, 46, 47, 48, 49};
    int arr1DLine5[]={51, 52, 53, 54, 55, 56, 57, 58, 59};
    int arr1DLine6[]={61, 62, 63, 64, 65, 66, 67, 68, 69};
    int arr1DLine7[]={71, 72, 73, 74, 75, 76, 77, 78, 79};
    int arr1DLine8[]={81, 82, 83, 84, 85, 86, 87, 88, 89};
    int arr1DLine9[]={91, 92, 93, 94, 95, 96, 97, 98, 99};
    */
    //
    /*
    int arr1DLine1[]={11, 12, 13, 14, 15, 16, 17, 18, 19};
    int arr1DLine2[]={21, 100,200,300,25, 26, 27, 28, 29};
    int arr1DLine3[]={31, 400,33, 34, 35, 36, 37, 38, 39};
    int arr1DLine4[]={41, 500,600,44, 45, 46, 47, 48, 49};
    int arr1DLine5[]={51, 52, 53, 54, 55, 56, 57, 58, 59};
    int arr1DLine6[]={61, 62, 63, 64, 65, 66, 67, 68, 69};
    int arr1DLine7[]={71, 72, 73, 74, 75, 76, 77, 78, 79};
    int arr1DLine8[]={81, 82, 83, 84, 85, 86, 87, 88, 89};
    int arr1DLine9[]={91, 92, 93, 94, 95, 96, 97, 98, 99};
    */
    //
    /*
    int arr1DLine1[]={11, 12, 13, 14, 15, 16, 17, 18, 19};
    int arr1DLine2[]={300,200,100,24, 25, 26, 27, 28, 29};
    int arr1DLine3[]={31, 32, 400,34, 35, 36, 37, 38, 39};
    int arr1DLine4[]={41, 600,500,44, 45, 46, 47, 48, 49};
    int arr1DLine5[]={51, 52, 53, 54, 55, 56, 57, 58, 59};
    int arr1DLine6[]={61, 62, 63, 64, 65, 66, 67, 68, 69};
    int arr1DLine7[]={71, 72, 73, 74, 75, 76, 77, 78, 79};
    int arr1DLine8[]={81, 82, 83, 84, 85, 86, 87, 88, 89};
    int arr1DLine9[]={91, 92, 93, 94, 95, 96, 97, 98, 99};
    */
    //
    /*
    int arr1DLine1[]={11, 12, 13, 14, 15, 16, 17, 18, 19};
    int arr1DLine2[]={21, 500,400,100,25, 26, 27, 28, 29};
    int arr1DLine3[]={31, 600,33, 200,35, 36, 37, 38, 39};
    int arr1DLine4[]={41, 42, 43, 300,45, 46, 47, 48, 49};
    int arr1DLine5[]={51, 52, 53, 54, 55, 56, 57, 58, 59};
    int arr1DLine6[]={61, 62, 63, 64, 65, 66, 67, 68, 69};
    int arr1DLine7[]={71, 72, 73, 74, 75, 76, 77, 78, 79};
    int arr1DLine8[]={81, 82, 83, 84, 85, 86, 87, 88, 89};
    int arr1DLine9[]={91, 92, 93, 94, 95, 96, 97, 98, 99};
    */
    //
    /*
    int arr1DLine1[]={11, 12, 13, 14, 15, 16, 17, 18, 19};
    int arr1DLine2[]={21, 22, 23, 100,400,500,27, 28, 29};
    int arr1DLine3[]={31, 32, 33, 200,35, 600,37, 38, 39};
    int arr1DLine4[]={41, 42, 43, 300,45, 46, 47, 48, 49};
    int arr1DLine5[]={51, 52, 53, 54, 55, 56, 57, 58, 59};
    int arr1DLine6[]={61, 62, 63, 64, 65, 66, 67, 68, 69};
    int arr1DLine7[]={71, 72, 73, 74, 75, 76, 77, 78, 79};
    int arr1DLine8[]={81, 82, 83, 84, 85, 86, 87, 88, 89};
    int arr1DLine9[]={91, 92, 93, 94, 95, 96, 97, 98, 99};
    */
    //
    /*
    int arr1DLine1[]={11, 12, 13, 100,15, 16, 17, 18, 19};
    int arr1DLine2[]={21, 22, 23, 200,25, 600,27, 28, 29};
    int arr1DLine3[]={31, 32, 33, 300,400,500,37, 38, 39};
    int arr1DLine4[]={41, 42, 43, 44, 45, 46, 47, 48, 49};
    int arr1DLine5[]={51, 52, 53, 54, 55, 56, 57, 58, 59};
    int arr1DLine6[]={61, 62, 63, 64, 65, 66, 67, 68, 69};
    int arr1DLine7[]={71, 72, 73, 74, 75, 76, 77, 78, 79};
    int arr1DLine8[]={81, 82, 83, 84, 85, 86, 87, 88, 89};
    int arr1DLine9[]={91, 92, 93, 94, 95, 96, 97, 98, 99};
    */
    //
    /*
    int arr1DLine1[]={11, 12, 13, 100,15, 16, 17, 18, 19};
    int arr1DLine2[]={21, 600,23, 200,25, 26, 27, 28, 29};
    int arr1DLine3[]={31, 500,400,300,35, 26, 37, 38, 39};
    int arr1DLine4[]={41, 42, 43, 44, 45, 46, 47, 48, 49};
    int arr1DLine5[]={51, 52, 53, 54, 55, 56, 57, 58, 59};
    int arr1DLine6[]={61, 62, 63, 64, 65, 66, 67, 68, 69};
    int arr1DLine7[]={71, 72, 73, 74, 75, 76, 77, 78, 79};
    int arr1DLine8[]={81, 82, 83, 84, 85, 86, 87, 88, 89};
    int arr1DLine9[]={91, 92, 93, 94, 95, 96, 97, 98, 99};
    */
    //
    /*
    int arr1DLine1[]={11, 12, 13, 14, 15, 16, 17, 18, 19};
    int arr1DLine2[]={21, 22, 23, 24, 25, 26, 27, 28, 29};
    int arr1DLine3[]={31, 32, 33, 34, 35, 36, 37, 38, 39};
    int arr1DLine4[]={41, 42, 43, 44, 45, 46, 47, 48, 49};
    int arr1DLine5[]={51, 52, 53, 54, 55, 56, 57, 58, 59};
    int arr1DLine6[]={61, 62, 63, 64, 65, 66, 67, 68, 69};
    int arr1DLine7[]={71, 72, 73, 74, 75, 76, 77, 78, 79};
    int arr1DLine8[]={81, 82, 83, 84, 85, 86, 87, 88, 89};
    int arr1DLine9[]={91, 92, 93, 94, 95, 96, 97, 98, 99};
    */
    //
    //----
    //
    /*
    int arr1DLine1[]={11, 12, 13, 14, 15, 16, 17, 18, 19};
    int arr1DLine2[]={21, 22, 23, 24, 25, 26, 27};
    int arr1DLine3[]={31, 32, 33, 34, 35, 36, 37, 38};
    int arr1DLine4[]={41, 42, 43, 44, 45, 46, 47, 48};
    int arr1DLine5[]={51, 52, 53, 54, 55, 56, 57};
    int arr1DLine6[]={61, 62, 63, 64, 65, 66, 67, 68, 69};
    int arr1DLine7[]={71, 72, 73, 74, 75, 76};
    int arr1DLine8[]={81, 82, 83, 84, 85};
    int arr1DLine9[]={91, 92, 93, 94, 95, 96, 97, 98};
    */
    //
    //
    //
    //int arr1DLine1[]={11, 12, 13, 14, 15, 16, 17, 18, 19};
    //int arr1DLine2[]={21, 100,200,300,25, 26, 27};
    //int arr1DLine3[]={31, 400,33, 34, 35, 36, 37, 38};
    //int arr1DLine4[]={41, 500,600,44, 45, 46, 47, 48};
    //int arr1DLine5[]={51, 52, 53, 54, 55, 56, 57};
    //int arr1DLine6[]={61, 62, 63, 64, 65, 66, 67, 68, 69};
    //int arr1DLine7[]={71, 72, 73, 74, 75, 76};
    //int arr1DLine8[]={81, 82, 83, 84, 85};
    //int arr1DLine9[]={91, 92, 93, 94, 95, 96, 97, 98};
    //
    Ls=new int[Q];
    Ls[1-1]=9;//Lini1;
    Ls[2-1]=7;//Lini2;
    Ls[3-1]=8;//Lini3;
    Ls[4-1]=8;//Lini4;
    Ls[5-1]=7;//Lini5;
    Ls[6-1]=9;//Lini6;
    Ls[7-1]=8;//=Lini7;
    Ls[8-1]=6;//Lini8;
    Ls[9-1]=8;//Lini9;
    size1.Set(Q, Ls);
    //
    int arr1DLine1[]={11, 12,  13,  14,  15, 16, 17, 18, 19};
    int arr1DLine2[]={21, 100, 200, 300, 25, 26, 27};
    int arr1DLine3[]={31, 400, 33,  300, 200, 100, 37, 38};
    int arr1DLine4[]={41, 500, 600, 44,  45,  400, 47, 48};
    int arr1DLine5[]={51, 300, 53,  54,  600, 500, 57};
    int arr1DLine6[]={61, 200, 63,  600, 100, 400, 500, 68, 69};
    int arr1DLine7[]={71, 100, 400, 500, 200, 76,  600, 78};
    int arr1DLine8[]={81, 82,  83,  84,  300, 86};
    int arr1DLine9[]={91, 92,  93,  94,  95, 96, 97, 98};
    //
    int arr1DN2Line1[]={100, 200, 300};
    int arr1DN2Line2[]={400};
    int arr1DN2Line3[]={500, 600};
    //
    /*
    int arr1DLine1[]={11, 12, 13, 14, 15, 16, 17, 18, 19};
    int arr1DLine2[]={300,200,100,24, 25, 26, 27};
    int arr1DLine3[]={31, 32, 400,34, 35, 36, 37, 38};
    int arr1DLine4[]={41, 600,500,44, 45, 46, 47, 48};
    int arr1DLine5[]={51, 52, 53, 54, 55, 56, 57};
    int arr1DLine6[]={61, 62, 63, 64, 65, 66, 67, 68, 69};
    int arr1DLine7[]={71, 72, 73, 74, 75, 76};
    int arr1DLine8[]={81, 82, 83, 84, 85};
    int arr1DLine9[]={91, 92, 93, 94, 95, 96, 97, 98};
    */
    //
    /*
    int arr1DLine1[]={11, 12, 13, 14, 15, 16, 17, 18, 19};
    int arr1DLine2[]={21, 500,400,100,25, 26, 27};
    int arr1DLine3[]={31, 600,33, 200,35, 36, 37, 38};
    int arr1DLine4[]={41, 42, 43, 300,45, 46, 47, 48};
    int arr1DLine5[]={51, 52, 53, 54, 55, 56, 57};
    int arr1DLine6[]={61, 62, 63, 64, 65, 66, 67, 68, 69};
    int arr1DLine7[]={71, 72, 73, 74, 75, 76};
    int arr1DLine8[]={81, 82, 83, 84, 85};
    int arr1DLine9[]={91, 92, 93, 94, 95, 96, 97, 98};
    */
    //
    /*
    int arr1DLine1[]={11, 12, 13, 14, 15, 16, 17, 18, 19};
    int arr1DLine2[]={21, 22, 23, 100,400,500,27};
    int arr1DLine3[]={31, 32, 33, 200,35, 600,37, 38};
    int arr1DLine4[]={41, 42, 43, 300,45, 46, 47, 48};
    int arr1DLine5[]={51, 52, 53, 54, 55, 56, 57};
    int arr1DLine6[]={61, 62, 63, 64, 65, 66, 67, 68, 69};
    int arr1DLine7[]={71, 72, 73, 74, 75, 76};
    int arr1DLine8[]={81, 82, 83, 84, 85};
    int arr1DLine9[]={91, 92, 93, 94, 95, 96, 97, 98};
    */
    //
    /*
    int arr1DLine1[]={11, 12, 13, 100,15, 16, 17, 18, 19};
    int arr1DLine2[]={21, 22, 23, 200,25, 600,27};
    int arr1DLine3[]={31, 32, 33, 300,400,500,37, 38};
    int arr1DLine4[]={41, 42, 43, 44, 45, 46, 47, 48};
    int arr1DLine5[]={51, 52, 53, 54, 55, 56, 57};
    int arr1DLine6[]={61, 62, 63, 64, 65, 66, 67, 68, 69};
    int arr1DLine7[]={71, 72, 73, 74, 75, 76};
    int arr1DLine8[]={81, 82, 83, 84, 85};
    int arr1DLine9[]={91, 92, 93, 94, 95, 96, 97, 98};
    */
    //
    /*
    int arr1DLine1[]={11, 12, 13, 100,15, 16, 17, 18, 19};
    int arr1DLine2[]={21, 600,23, 200,25, 26, 27};
    int arr1DLine3[]={31, 500,400,300,35, 26, 37, 38};
    int arr1DLine4[]={41, 42, 43, 44, 45, 46, 47, 48};
    int arr1DLine5[]={51, 52, 53, 54, 55, 56, 57};
    int arr1DLine6[]={61, 62, 63, 64, 65, 66, 67, 68, 69};
    int arr1DLine7[]={71, 72, 73, 74, 75, 76};
    int arr1DLine8[]={81, 82, 83, 84, 85};
    int arr1DLine9[]={91, 92, 93, 94, 95, 96, 97, 98};
    */
    //
    /*
    int arr1DLine1[]={11, 12, 13, 14, 15, 16, 17, 18, 19};
    int arr1DLine2[]={21, 22, 23, 24, 25, 26, 27};
    int arr1DLine3[]={31, 32, 33, 34, 35, 36, 37, 38};
    int arr1DLine4[]={41, 42, 43, 44, 45, 46, 47, 48};
    int arr1DLine5[]={51, 52, 53, 54, 55, 56, 57};
    int arr1DLine6[]={61, 62, 63, 64, 65, 66, 67, 68, 69};
    int arr1DLine7[]={71, 72, 73, 74, 75, 76};
    int arr1DLine8[]={81, 82, 83, 84, 85};
    int arr1DLine9[]={91, 92, 93, 94, 95, 96, 97, 98};
    */
    int arr1DN1[]={101, 102, 103, 104, 105, 106, 107, 108, 109};
    int arr1DN2[]={201, 202, 203, 204, 205, 206, 207, 208, 209};
    //
    int *arr1DLine1Ptr=arr1DLine1,
        *arr1DLine2Ptr=arr1DLine2,
        *arr1DLine3Ptr=arr1DLine3,
        *arr1DLine4Ptr=arr1DLine4,
        *arr1DLine5Ptr=arr1DLine5,
        *arr1DLine6Ptr=arr1DLine6,
        *arr1DLine7Ptr=arr1DLine7,
        *arr1DLine8Ptr=arr1DLine8,
        *arr1DLine9Ptr=arr1DLine9,
        //
        *arr1DN2Line1Ptr=arr1DN2Line1,
        *arr1DN2Line2Ptr=arr1DN2Line2,
        *arr1DN2Line3Ptr=arr1DN2Line3,
        //
        *arr1DN1Ptr=arr1DN1,
        *arr1DN2Ptr=arr1DN2;
    int **arr2DN1Ptr, **arr2DN2Ptr;


    //
    arr2DN1Ptr=new int*[Q];
    arr2DN1Ptr[1-1]=arr1DLine1Ptr;
    arr2DN1Ptr[2-1]=arr1DLine2Ptr;
    arr2DN1Ptr[3-1]=arr1DLine3Ptr;
    arr2DN1Ptr[4-1]=arr1DLine4Ptr;
    arr2DN1Ptr[5-1]=arr1DLine5Ptr;
    arr2DN1Ptr[6-1]=arr1DLine6Ptr;
    arr2DN1Ptr[7-1]=arr1DLine7Ptr;
    arr2DN1Ptr[8-1]=arr1DLine8Ptr;
    arr2DN1Ptr[9-1]=arr1DLine9Ptr;
    //
    L1s=new int[3];
    L1s[1-1]=3;
    L1s[2-1]=1;
    L1s[3-1]=2;
    arr2DN2Ptr=new int*[3];
    arr2DN2Ptr[1-1]=arr1DN2Line1Ptr;
    arr2DN2Ptr[2-1]=arr1DN2Line2Ptr;
    arr2DN2Ptr[3-1]=arr1DN2Line3Ptr;
    //
    size2.Set(3, L1s);
    //
    std::cout<<"Ini 2D array:"<<std::endl;
    Arr2DStdCOut1D(arr2DN1Ptr, size1);
    std::cout<<std::endl;
    std::cout<<std::endl;
    Arr2DStdCOut1D(arr2DN1Ptr, size1, " ", "; ", true);
    std::cout<<std::endl;
    std::cout<<std::endl;
    Arr2DStdCOut1D(arr2DN1Ptr, size1, " ", "", true);
    std::cout<<std::endl;
    std::cout<<std::endl;
    Arr2DStdCOut1D(arr2DN1Ptr, size1, " ", "; ", false);
    std::cout<<std::endl;
    std::cout<<std::endl;
    Arr2DStdCOut1D(arr2DN1Ptr, size1, ", ", "; ", false);
    std::cout<<std::endl;
    std::cout<<std::endl;
    Arr2DStdCOut2D(arr2DN1Ptr, size1);
    std::cout<<std::endl;
    std::cout<<std::endl;
    Arr2DStdCOut2D(arr2DN1Ptr, size1, " ", "; ", true);
    std::cout<<std::endl;
    Arr2DStdCOut2D(arr2DN1Ptr, size1, " ", "", true);
    std::cout<<std::endl;
    Arr2DStdCOut2D(arr2DN1Ptr, size1, " ", "; ", false);
    std::cout<<std::endl;
    //
    bool  IsTransposed=false, EachExtRowIsReversed=false;
    std::vector<std::vector<int>> Ns;
    //
    std::vector<int> arr1DVBN1, arr1DVBN2, curRow;
    std::vector<std::vector<int>> arr2DVBN1, arr2DVBN2, arr2DVBN3;
    //
    Array2DAssign(arr2DVBN1, size1, arr2DN1Ptr, DfltInt);
    Array2DAssign(arr2DVBN2, size2, arr2DN2Ptr, DfltInt);
    Array1DAssign(arr1DVBN1, arr1DN1, 9);
    Array1DAssign(arr1DVBN2, arr1DN2, 9);
    //--------------------------------------------------------------
    //Collecting Garbage initial
    std::cout<<"Fin. Collecting Garbage"<<std::endl;
    //
    delete[]Ls;
    delete[]L1s;
    //
    if(arr2DN1Ptr!=NULL){
        Q=size1.GetQExtRows();
        for(int i=1; i<=Q; i++){
            if(arr2DN1Ptr[i-1]!=NULL){
                delete[]arr2DN1Ptr[i-1];
            }
        }
    }
    delete[]arr2DN1Ptr;
    //
    for(int i=1; i<=3; i++){
        if(arr2DN2Ptr[i-1]!=NULL){
            delete[]arr2DN2Ptr[i-1];
        }
    }
    delete[]arr2DN2Ptr;
    //
   // delete DfltInt;//NO at the end
    //
    //--------------------------------------------------------------
    //
    std::cout<<" Testing search ";
    //
    IsTransposed=false;
    EachExtRowIsReversed=false;
    if(IsTransposed)std::cout<<" Row to seek IS transposed ";
    else std::cout<<" Row to seek is NOT transposed ";
    std::cout<<";";
    if(EachExtRowIsReversed)std::cout<<" Each Ext Row IS Reversed ";
    else std::cout<<" Each Ext Row is NOT Reversed ";
    std::cout<<std::endl;
    Ns= Array2D_SeekArrs2ndIn1st(arr2DVBN1, arr2DVBN2, IsTransposed, EachExtRowIsReversed);
    
    for(int i=1; i<=Ns.size(); i++){
        std::cout<<"ExtRowN: "<<Ns[i-1][1-1]<<" IneRowN: "<<Ns[i-1][2-1]<<std::endl;
    }
    std::cout<<"Must be correct answer: Found 1 ExtRowN: 2 IneRowN: 2"<<std::endl;
    std::cout<<std::endl;
    //
    IsTransposed=false;
    EachExtRowIsReversed=true;
    if(IsTransposed)std::cout<<" Row to seek IS transposed ";
    else std::cout<<" Row to seek is NOT transposed ";
    std::cout<<";";
    if(EachExtRowIsReversed)std::cout<<" Each Ext Row IS Reversed ";
    else std::cout<<" Each Ext Row is NOT Reversed ";
    std::cout<<std::endl;
    Ns= Array2D_SeekArrs2ndIn1st(arr2DVBN1, arr2DVBN2, IsTransposed, EachExtRowIsReversed);
    std::cout<<"Found: "<<Ns.size()<<std::endl;
    for(int i=1; i<=Ns.size(); i++){
        std::cout<<"ExtRowN: "<<Ns[i-1][1-1]<<" IneRowN: "<<Ns[i-1][2-1]<<std::endl;
    }
    std::cout<<"Must be correct answer: Found 1 ExtRowN: 3 IneRowN: 6"<<std::endl;
    std::cout<<std::endl;
    //
    IsTransposed=true;
    EachExtRowIsReversed=false;
    if(IsTransposed)std::cout<<" Row to seek IS transposed ";
    else std::cout<<" Row to seek is NOT transposed ";
    std::cout<<";";
    if(EachExtRowIsReversed)std::cout<<" Each Ext Row IS Reversed ";
    else std::cout<<" Each Ext Row is NOT Reversed ";
    std::cout<<std::endl;
    Ns= Array2D_SeekArrs2ndIn1st(arr2DVBN1, arr2DVBN2, IsTransposed, EachExtRowIsReversed);
    std::cout<<"Found: "<<Ns.size()<<std::endl;
    for(int i=1; i<=Ns.size(); i++){
        std::cout<<"ExtRowN: "<<Ns[i-1][1-1]<<" IneRowN: "<<Ns[i-1][2-1]<<std::endl;
    }
    std::cout<<"Must be correct answer: Found 6 ExtRowN: 6 IneRowN: 5"<<std::endl;
    std::cout<<std::endl;
    //
    IsTransposed=true;
    EachExtRowIsReversed=true;
    if(IsTransposed)std::cout<<" Row to seek IS transposed ";
    else std::cout<<" Row to seek is NOT transposed ";
    std::cout<<";";
    if(EachExtRowIsReversed)std::cout<<" Each Ext Row IS Reversed ";
    else std::cout<<" Each Ext Row is NOT Reversed ";
    std::cout<<std::endl;
    Ns= Array2D_SeekArrs2ndIn1st(arr2DVBN1, arr2DVBN2, IsTransposed, EachExtRowIsReversed);
    std::cout<<"Found: "<<Ns.size()<<std::endl;
    for(int i=1; i<=Ns.size(); i++){
        std::cout<<"ExtRowN: "<<Ns[i-1][1-1]<<" IneRowN: "<<Ns[i-1][2-1]<<std::endl;
    }
    std::cout<<"Must be correct answer: Found 1 ExtRowN: 7 IneRowN: 2"<<std::endl;
    std::cout<<std::endl;
    //
    std::cout<<"Press any key"<<std::endl;
    getch();
    //
    std::cout<<"Transposing:"<<std::endl;
    //
    Array2DTranspose(arr2DVBN1);
    Arr2DStdCOut2D(arr2DVBN1, " ", "; ", false);
    std::cout<<std::endl;
    size1=GetSizeOf2DVector(arr2DVBN1);
    std::cout<<"Size now: "<<size1.GetAsStdString().c_str()<<std::endl;
    Q=size1.GetQExtRows();
    std::cout<<"Transposing again:"<<std::endl;
    Array2DTranspose(arr2DVBN1);
    Arr2DStdCOut2D(arr2DVBN1, " ", "; ", false);
    std::cout<<std::endl;
    size1=GetSizeOf2DVector(arr2DVBN1);
    std::cout<<"Size now: "<<size1.GetAsStdString().c_str()<<std::endl;
    Q=size1.GetQExtRows();
    //
    std::cout<<"Press any key"<<std::endl;
    getch();
    //
    std::cout<<"Deleting from ext rows:"<<std::endl;
    for(int i=2; i<=Q; i++){
        for(int j=Q; j>=Q+1-i+1; j--){
            //Lcur=size
            Array2DDelFromExtRowN(arr2DVBN1, i, j);//size1.GetMaxLength()-j+1);
        }
    }
    Arr2DStdCOut2D(arr2DVBN1, " ", "; ", false);
    std::cout<<std::endl;
    size1=GetSizeOf2DVector(arr2DVBN1);
    std::cout<<"Size now: "<<size1.GetAsStdString().c_str()<<std::endl;
    Q=size1.GetQExtRows();
    //
    //
    std::cout<<"Adding line:" <<std::endl;
    //template<typename T>void Array2DAddExtRow(T**&X, Array2DSize&size, T*rowParam, int wholeRowL=0, int whatFromN=1, int QDefaultValsBefore=0, bool RectNotVar=true, T*DfltValParam=NULL, TValsShowHide*vsh=NULL){
    Array2DAddExtRow(arr2DVBN1, arr1DVBN1, LiniExN1);
    Arr2DStdCOut2D(arr2DVBN1, " ", "; ", false);
    std::cout<<std::endl;
    size1=GetSizeOf2DVector(arr2DVBN1);
    Q=size1.GetQExtRows();
    //LcurRow=9;
    //Arr1DStdCOut(arr2DN1Ptr[Q-1], LcurRow);
    std::cout<<std::endl;
    //
    std::cout<<"Deleting line from "<<Q<<" position :"<<std::endl;
    Array2DDelExtRow(arr2DVBN1, Q);
    Arr2DStdCOut2D(arr2DVBN1, " ", "; ", false);
    std::cout<<std::endl;
    size1=GetSizeOf2DVector(arr2DVBN1);
    std::cout<<"Size now: "<<size1.GetAsStdString().c_str()<<std::endl;
    Q=size1.GetQExtRows();
    //
    ExtRowN=3;
    std::cout<<"Inserting line to "<<ExtRowN<<" position :"<<std::endl;
    //Array2DInsExtRow(arr2DN1Ptr, size1, ExtRowN, 9, arr1DN1Ptr, DfltInt, 1, 0, true, &vsh);
    //template<typename T> void Array2DInsExtRow(std::vector<std::vector<T>>&X, int ExtRowN, std::vector<T>rowParam, T*DefaultValParam=NULL, int whatFromN=1, int QDefaultValuesBefore=0, bool KeepIfRect=true){
    Array2DInsExtRow(arr2DVBN1, ExtRowN, arr1DVBN1, DfltInt, 1, 0, true);
    Arr2DStdCOut2D(arr2DVBN1, " ", "; ", false);
    std::cout<<std::endl;
    Q=size1.GetQExtRows();
    //
    std::cout<<"Deleting line from "<<ExtRowN<<" position :"<<std::endl;
    Array2DDelExtRow(arr2DVBN1, ExtRowN);
    Arr2DStdCOut2D(arr2DVBN1, " ", "; ", false);
    std::cout<<std::endl;
    size1=GetSizeOf2DVector(arr2DVBN1);
    std::cout<<"Size now: "<<size1.GetAsStdString().c_str()<<std::endl;
    Q=size1.GetQExtRows();
    //
    ExtRowN=1;
    std::cout<<"Inserting line to "<<ExtRowN<<" position :"<<std::endl;
    //Array2DInsExtRow(arr2DN1Ptr, size1, ExtRowN, 9, arr1DN1Ptr, DfltInt, 1, 0, true, &vsh);
    //template<typename T> void Array2DInsExtRow(std::vector<std::vector<T>>&X, int ExtRowN, std::vector<T>rowParam, T*DefaultValParam=NULL, int whatFromN=1, int QDefaultValuesBefore=0, bool KeepIfRect=true){
    Array2DInsExtRow(arr2DVBN1,  ExtRowN,  arr1DVBN1, DfltInt, 1, 0, true);
    Arr2DStdCOut2D(arr2DVBN1, " ", "; ", false);
    std::cout<<std::endl;
    Q=size1.GetQExtRows();
    //
    std::cout<<"Setting line at "<<ExtRowN<<" position :"<<std::endl;
    //void Array2DSetExtRowN(T**&X, Array2DSize&size, int ExtRowN, int wholeRowL=0, T*rowParam=NULL, int whatFromN=1, int QDefaultValuesBefore=0, bool KeepIfRect=true, T*DefaultValParam=NULL, bool DefaultIsSpecNotOwn=false){
    Array2DSetExtRowN(arr2DVBN1, ExtRowN, arr1DVBN2);//, int whatFromN=1, int QDefaultValuesBefore=0, bool KeepIfRect=true, T*DefaultValParam=NULL, bool DefaultIsSpecNotOwn=false){
    Arr2DStdCOut2D(arr2DVBN1, " ", "; ", false);
    std::cout<<std::endl;
    size1=GetSizeOf2DVector(arr2DVBN1);
    Q=size1.GetQExtRows();
    std::cout<<"Press any key:"<<std::endl;
    getch();
    //
    std::cout<<"Deleting line from "<<ExtRowN<<" position :"<<std::endl;
    Array2DDelExtRow(arr2DVBN1, ExtRowN);
    Arr2DStdCOut2D(arr2DVBN1, " ", "; ", false);
    std::cout<<std::endl;
    size1=GetSizeOf2DVector(arr2DVBN1);
    std::cout<<"Size now: "<<size1.GetAsStdString().c_str()<<std::endl;
    Q=size1.GetQExtRows();
    //
    ExtRowN=Q;
    std::cout<<"Inserting line to (last) position :"<<std::endl;
    //Array2DInsExtRow(arr2DN1Ptr, size1, ExtRowN, 9, arr1DN1Ptr, DfltInt, 1, 0, true, &vsh);
    Array2DInsExtRow(arr2DVBN1, ExtRowN, arr1DVBN1, DfltInt, 1, 0, true);
    Arr2DStdCOut2D(arr2DVBN1, " ", "; ", false);
    std::cout<<std::endl;
    Q=size1.GetQExtRows();
    //
    std::cout<<"Deleting line from "<<ExtRowN<<" position :"<<std::endl;
    Array2DDelExtRow(arr2DVBN1, ExtRowN);
    Arr2DStdCOut2D(arr2DVBN1, " ", "; ", false);
    std::cout<<std::endl;
    size1=GetSizeOf2DVector(arr2DVBN1);
    std::cout<<"Size now: "<<size1.GetAsStdString().c_str()<<std::endl;
    Q=size1.GetQExtRows();
    //
    // Add, Ins, Del Swap ext rows arb, et se arb
    //
    std::cout<<"Adding column (defaults):"<<std::endl;
    //template<typename T>void Array2DAddExtRow(T**&X, Array2DSize&size, T*rowParam, int wholeRowL=0, int whatFromN=1, int QDefaultValsBefore=0, bool RectNotVar=true, T*DfltValParam=NULL, TValsShowHide*vsh=NULL){
    Array2DAddIneRow(arr2DVBN1, arr1DVBN1);
    Arr2DStdCOut2D(arr2DVBN1, " ", "; ", false);
    std::cout<<std::endl;
    size1=GetSizeOf2DVector(arr2DVBN1);
    Q=size1.GetQExtRows();
    //LcurRow=9;
    //Arr1DStdCOut(arr2DN1Ptr[Q-1], LcurRow);
    std::cout<<"defaults mean do nil if varia L"<<std::endl;
    std::cout<<std::endl;
    //
    std::cout<<"Adding column (do nil if varia L):"<<std::endl;
    //template<typename T>void Array2DAddExtRow(T**&X, Array2DSize&size, T*rowParam, int wholeRowL=0, int whatFromN=1, int QDefaultValsBefore=0, bool RectNotVar=true, T*DfltValParam=NULL, TValsShowHide*vsh=NULL){
    //template<typename T>Array2DAddIneRow(std::vector<std::vector<T>>&X, std::vector<T>rowParam, int IfNonRectIgnore0Add1Stretch2=0, int whatFromN=1, int QDefaultValsBefore=0, T*DefaultValParam=NULL){
    Array2DAddIneRow(arr2DVBN1, arr1DVBN1, 0, 1, 0, DfltInt);
    Arr2DStdCOut2D(arr2DVBN1, " ", "; ", false);
    std::cout<<std::endl;
    size1=GetSizeOf2DVector(arr2DVBN1);
    Q=size1.GetQExtRows();
    //LcurRow=9;
    //Arr1DStdCOut(arr2DN1Ptr[Q-1], LcurRow);
    std::cout<<std::endl;
    //
    std::cout<<"Adding column (add to last pos of varia L lines):"<<std::endl;
    //template<typename T>void Array2DAddExtRow(T**&X, Array2DSize&size, T*rowParam, int wholeRowL=0, int whatFromN=1, int QDefaultValsBefore=0, bool RectNotVar=true, T*DfltValParam=NULL, TValsShowHide*vsh=NULL){
    Array2DAddIneRow(arr2DVBN1, arr1DVBN1, 1);
    Arr2DStdCOut2D(arr2DVBN1, " ", "; ", false);
    std::cout<<std::endl;
    size1=GetSizeOf2DVector(arr2DVBN1);
    Q=size1.GetQExtRows();
    //LcurRow=9;
    //Arr1DStdCOut(arr2DN1Ptr[Q-1], LcurRow);
    std::cout<<"Size now: "<<size1.GetAsStdString().c_str()<<std::endl;
    std::cout<<std::endl;
    //
    std::cout<<"Deleting line from last position of each var L line):"<<std::endl;
    Array2DDelIneRow(arr2DVBN1, 0, true);
    Arr2DStdCOut2D(arr2DVBN1, " ", "; ", false);
    std::cout<<std::endl;
    size1=GetSizeOf2DVector(arr2DVBN1);
    std::cout<<"Size now: "<<size1.GetAsStdString().c_str()<<std::endl;
    Q=size1.GetQExtRows();
    //
    std::cout<<"Adding column (stretch to Lmax if varia Ls):"<<std::endl;
    //template<typename T>void Array2DAddExtRow(T**&X, Array2DSize&size, T*rowParam, int wholeRowL=0, int whatFromN=1, int QDefaultValsBefore=0, bool RectNotVar=true, T*DfltValParam=NULL, TValsShowHide*vsh=NULL){
    //Array2DAddIneRow(arr2DN1Ptr, size1, arr1DN1Ptr, 9, 2);
    Array2DAddIneRow(arr2DVBN1, arr1DVBN1, 2, 1, 0, DfltInt);
    Arr2DStdCOut2D(arr2DVBN1, " ", "; ", false);
    std::cout<<std::endl;
    size1=GetSizeOf2DVector(arr2DVBN1);
    Q=size1.GetQExtRows();
    //LcurRow=9;
    //Arr1DStdCOt(arr2DN1Ptr[Q-1], LcurRow);
    std::cout<<"Size now: "<<size1.GetAsStdString().c_str()<<std::endl;
    std::cout<<std::endl;
    //
    std::cout<<"Setting column (N Lmax):"<<std::endl;
    //template<typename T>Array2DSetIneRow(T**&X, Array2DSize& size, int IneRowN, T*rowParam=NULL, int wholeRowL=0, int whatFromN=1, int QDefaultValsBefore=0, T*DefaultValParam=NULL, bool ExistingInsteadOfDefault=true, bool ExistingOnlyForGTLminNotIgnore=false, bool LastPossIfNotRectAtPos0NotIgnore=false){
    size1=GetSizeOf2DVector(arr2DVBN1);
    Array2DSetIneRow(arr2DVBN1, size1.GetMaxLength(), arr1DVBN2);
    Arr2DStdCOut2D(arr2DVBN1, " ", "; ", false);
    std::cout<<std::endl;
    size1=GetSizeOf2DVector(arr2DVBN1);
    Q=size1.GetQExtRows();
    std::cout<<"Size now: "<<size1.GetAsStdString().c_str()<<std::endl;
    std::cout<<std::endl;
    //
    std::cout<<"Deleting line from "<<size1.GetMaxLength()<<" position :"<<std::endl;
    size1=GetSizeOf2DVector(arr2DVBN1);
    Array2DDelIneRow(arr2DVBN1, size1.GetMaxLength());
    size1=GetSizeOf2DVector(arr2DVBN1);
    Arr2DStdCOut2D(arr2DVBN1, " ", "; ", false);
    std::cout<<std::endl;
    size1=GetSizeOf2DVector(arr2DVBN1);
    std::cout<<"Size now: "<<size1.GetAsStdString().c_str()<<std::endl;
    Q=size1.GetQExtRows();
    ////tried ady del lmax-1 - arb gut. Nu testing ins
    //Array2DInsIneRow(T**&X, Array2DSize& size, int IneRowPosN, T*rowParam=NULL, int wholeRowL=0, int IfAfterLminIgnore0Stretch2=0, int whatFromN=1, int QDefaultValsBefore=0, T*DefaultValParam=NULL)
    //
    IneRowN=1;
    std::cout<<"Inserting col to "<<IneRowN<<" position :"<<std::endl;
    Array2DInsIneRow(arr2DVBN1, IneRowN, arr1DVBN1, 2);//, int whatFromN=1, int QDefaultValsBefore=0, T*DefaultValParam=NULL)
    Arr2DStdCOut2D(arr2DVBN1, " ", "; ", false);
    std::cout<<std::endl;
    size1=GetSizeOf2DVector(arr2DVBN1);
    std::cout<<"Size now: "<<size1.GetAsStdString().c_str()<<std::endl;
    //arr1DN1Ptr
    //
    std::cout<<"Deleting line from "<<IneRowN<<" position :"<<std::endl;
    Array2DDelIneRow(arr2DVBN1, IneRowN);
    Arr2DStdCOut2D(arr2DVBN1, " ", "; ", false);
    std::cout<<std::endl;
    size1=GetSizeOf2DVector(arr2DVBN1);
    std::cout<<"Size now: "<<size1.GetAsStdString().c_str()<<std::endl;
    Q=size1.GetQExtRows();
    ////tried ady del lmax-1 - arb gut. Nu testing ins
    //Array2DInsIneRow(T**&X, Array2DSize& size, int IneRowPosN, T*rowParam=NULL, int wholeRowL=0, int IfAfterLminIgnore0Stretch2=0, int whatFromN=1, int QDefaultValsBefore=0, T*DefaultValParam=NULL)
    //
    IneRowN=2;
    std::cout<<"Inserting col to "<<IneRowN<<" position :"<<std::endl;
    //Array2DInsIneRow(std::vector<std::vector<T>>&X, int IneRowPosN, std::vector<T>rowParam, int IfAfterLminIgnore0Stretch2=0, int whatFromN=1, int QDefaultValsBefore=0, T*DefaultValParam=NULL){
    Array2DInsIneRow(arr2DVBN1, IneRowN, arr1DVBN1, 2);//, int whatFromN=1, int QDefaultValsBefore=0, T*DefaultValParam=NULL)
    Arr2DStdCOut2D(arr2DVBN1, " ", "; ", false);
    std::cout<<std::endl;
    size1=GetSizeOf2DVector(arr2DVBN1);
    std::cout<<"Size now: "<<size1.GetAsStdString().c_str()<<std::endl;
    //arr1DN1Ptr
    //
    std::cout<<"Deleting line from "<<IneRowN<<" position :"<<std::endl;
    Array2DDelIneRow(arr2DVBN1, IneRowN);
    size1=GetSizeOf2DVector(arr2DVBN1);
    Arr2DStdCOut2D(arr2DVBN1, " ", "; ", false);
    std::cout<<std::endl;
    size1=GetSizeOf2DVector(arr2DVBN1);
    std::cout<<"Size now: "<<size1.GetAsStdString().c_str()<<std::endl;
    Q=size1.GetQExtRows();
    ////tried ady del lmax-1 - arb gut. Nu testing ins
    //Array2DInsIneRow(T**&X, Array2DSize& size, int IneRowPosN, T*rowParam=NULL, int wholeRowL=0, int IfAfterLminIgnore0Stretch2=0, int whatFromN=1, int QDefaultValsBefore=0, T*DefaultValParam=NULL)
    //
    IneRowN=size1.GetMaxLength();
    std::cout<<"Inserting col to "<<IneRowN<<" position :"<<std::endl;
    Array2DInsIneRow(arr2DVBN1, IneRowN, arr1DN1Ptr, 2);//, int whatFromN=1, int QDefaultValsBefore=0, T*DefaultValParam=NULL)
    Arr2DStdCOut2D(arr2DVBN1, " ", "; ", false);
    std::cout<<std::endl;
    size1=GetSizeOf2DVector(arr2DVBN1);
    std::cout<<"Size now: "<<size1.GetAsStdString().c_str()<<std::endl;
    //arr1DN1Ptr
    //
    std::cout<<"Deleting line from "<<IneRowN<<" position :"<<std::endl;
    Array2DDelIneRow(arr2DVBN1, IneRowN);
    Arr2DStdCOut2D(arr2DVBN1, " ", "; ", false);
    std::cout<<std::endl;
    size1=GetSizeOf2DVector(arr2DVBN1);
    std::cout<<"Size now: "<<size1.GetAsStdString().c_str()<<std::endl;
    Q=size1.GetQExtRows();
    ////tried ady del lmax-1 - arb gut. Nu testing ins
    //Array2DInsIneRow(T**&X, Array2DSize& size, int IneRowPosN, T*rowParam=NULL, int wholeRowL=0, int IfAfterLminIgnore0Stretch2=0, int whatFromN=1, int QDefaultValsBefore=0, T*DefaultValParam=NULL)
    //
    // add, del, ins col arb gut, et set arb
    //
    //
    //Search, Concat
    // std::vector<std::vector<int>> Array2D_SeekArrs2ndIn1st(T**Where, const Array2DSize& sizeWhere, T**What, const Array2DSize& sizeWhat, bool IsTransposed=false, bool EachExtRowIsReversed=false){
    //
    std::cout<<" Testing search ";
    //
    IsTransposed=false;
    EachExtRowIsReversed=false;
    if(IsTransposed)std::cout<<" Row to seek IS transposed ";
    else std::cout<<" Row to seek is NOT transposed ";
    std::cout<<";";
    if(EachExtRowIsReversed)std::cout<<" Each Ext Row IS Reversed ";
    else std::cout<<" Each Ext Row is NOT Reversed ";
    std::cout<<std::endl;
    Ns = Array2D_SeekArrs2ndIn1st(arr2DVBN1, arr2DVBN2, IsTransposed, EachExtRowIsReversed);
    std::cout<<"Found: "<<Ns.size()<<std::endl;
    for(int i=1; i<=Ns.size(); i++){
        std::cout<<"ExtRowN: "<<Ns[i-1][1-1]<<" IneRowN: "<<Ns[i-1][2-1]<<std::endl;
    }
    std::cout<<"Must be correct answer: Found 1 ExtRowN: 2 IneRowN: 2"<<std::endl;
    std::cout<<std::endl;
    //
    IsTransposed=false;
    EachExtRowIsReversed=true;
    if(IsTransposed)std::cout<<" Row to seek IS transposed ";
    else std::cout<<" Row to seek is NOT transposed ";
    std::cout<<";";
    if(EachExtRowIsReversed)std::cout<<" Each Ext Row IS Reversed ";
    else std::cout<<" Each Ext Row is NOT Reversed ";
    std::cout<<std::endl;
    Ns = Array2D_SeekArrs2ndIn1st(arr2DVBN1, arr2DVBN2, IsTransposed, EachExtRowIsReversed);
    std::cout<<"Found: "<<Ns.size()<<std::endl;
    for(int i=1; i<=Ns.size(); i++){
        std::cout<<"ExtRowN: "<<Ns[i-1][1-1]<<" IneRowN: "<<Ns[i-1][2-1]<<std::endl;
    }
    std::cout<<"Must be correct answer: Found 1 ExtRowN: 3 IneRowN: 6"<<std::endl;
    std::cout<<std::endl;
    //
    IsTransposed=true;
    EachExtRowIsReversed=false;
    if(IsTransposed)std::cout<<" Row to seek IS transposed ";
    else std::cout<<" Row to seek is NOT transposed ";
    std::cout<<";";
    if(EachExtRowIsReversed)std::cout<<" Each Ext Row IS Reversed ";
    else std::cout<<" Each Ext Row is NOT Reversed ";
    std::cout<<std::endl;
    Ns = Array2D_SeekArrs2ndIn1st(arr2DVBN1, arr2DVBN2, IsTransposed, EachExtRowIsReversed);
    std::cout<<"Found: "<<Ns.size()<<std::endl;
    for(int i=1; i<=Ns.size(); i++){
        std::cout<<"ExtRowN: "<<Ns[i-1][1-1]<<" IneRowN: "<<Ns[i-1][2-1]<<std::endl;
    }
    std::cout<<"Must be correct answer: Found 6 ExtRowN: 6 IneRowN: 5"<<std::endl;
    std::cout<<std::endl;
    //
    IsTransposed=true;
    EachExtRowIsReversed=true;
    if(IsTransposed)std::cout<<" Row to seek IS transposed ";
    else std::cout<<" Row to seek is NOT transposed ";
    std::cout<<";";
    if(EachExtRowIsReversed)std::cout<<" Each Ext Row IS Reversed ";
    else std::cout<<" Each Ext Row is NOT Reversed ";
    std::cout<<std::endl;
    Ns = Array2D_SeekArrs2ndIn1st(arr2DVBN1, arr2DVBN2, IsTransposed, EachExtRowIsReversed);
    std::cout<<"Found: "<<Ns.size()<<std::endl;
    for(int i=1; i<=Ns.size(); i++){
        std::cout<<"ExtRowN: "<<Ns[i-1][1-1]<<" IneRowN: "<<Ns[i-1][2-1]<<std::endl;
    }
    std::cout<<"Must be correct answer: Found 1 ExtRowN: 7 IneRowN: 2"<<std::endl;
    std::cout<<std::endl;
    //
    std::cout<<"Press any key"<<std::endl;
    getch();
    //
    std::cout<<std::endl;
    //
    std::cout<<"Transposing:"<<std::endl;
    //
    Array2DTranspose(arr2DVBN1);
    Arr2DStdCOut2D(arr2DVBN1, " ", "; ", false);
    std::cout<<std::endl;
    size1=GetSizeOf2DVector(arr2DVBN1);
    std::cout<<"Size now: "<<size1.GetAsStdString().c_str()<<std::endl;
    Q=size1.GetQExtRows();
    std::cout<<"Transposing again:"<<std::endl;
    Array2DTranspose(arr2DVBN1);
    Arr2DStdCOut2D(arr2DVBN1, " ", "; ", false);
    std::cout<<std::endl;
    std::cout<<"Size now: "<<size1.GetAsStdString().c_str()<<std::endl;
    Q=size1.GetQExtRows();
    //
    //
    std::cout<<"Press any key"<<std::endl;
    getch();
    //
    //
    std::cout<<"Testing assign and concat"<<std::endl;
    //
    //Array2DAssign(arr2DVBN2, arr2DVBN1);
    arr2DVBN2= arr2DVBN1;
    //
    size1=GetSizeOf2DVector(arr2DVBN1);
    Q=size1.GetQExtRows();
    std::cout<<"Copied row:"<<std::endl;
    Arr2DStdCOut2D(arr2DVBN2, " ", "; ", false);
    std::cout<<std::endl;
    //
    int shiftQ=3;//3;//0;//works well by 0 and 3 and rect=true, ma exda
    bool rect=true;
    int QRec=0; //Q records to add ab added array
    int FromN=3; //Erst added N ab added array
    std::cout<<"Adding 2D arr:"<<std::endl;
    //
    std::cout<<"Now let us try concat adding..."<<std::endl;
    Array2DConcatAdding(arr2DVBN1, arr2DVBN2, shiftQ, DfltInt, rect, QRec, FromN);
    std::cout<<"Concatenated 2D arr (adding):"<<std::endl;
    //Arr2DStdCOut2D(arr2DN2Ptr, size2, " ", "; ", false);
    Arr2DStdCOut2D(arr2DVBN1, " ", "; ", false);
    std::cout<<std::endl;
    //works well  - tested for rect arrays,
    //rect for em n'vikts, shiftQ = 0 et 3, FromN = 1 et 3, QRec=0 et 2 //wa errs, fa corrs
    //
    bool Stretching_ifFirstNotRect_AddToCurNotStretchToMax=false;
    shiftQ=-3;//-3;//3
    FromN=1;
    //t T void Array2DConcatStretching(T**&X1, Array2DSize& size1, T**X2, const Array2DSize& size2, int shiftVal=0, T*DefaultValParam=NULL, bool Rect=true, bool IfFirstNotRect_AddToLastNotSrtetchToMax=true, int whatFromN=1, int QRec=0){
    Array2DConcatStretching(arr2DVBN1, arr2DVBN2, shiftQ, DfltInt, rect, Stretching_ifFirstNotRect_AddToCurNotStretchToMax, FromN, QRec, &vsh);
    std::cout<<"Concatenated 2D arr (stretching):"<<std::endl;
    //Arr2DStdCOut2D(arr2DN2Ptr, size2, " ", "; ", false);
    Arr2DStdCOut2D(arr2DVBN1, " ", "; ", false);
    std::cout<<std::endl;
    std::cout<<"trying to repeat..."<<std::endl;
    for(int i=1; i<=10; i++){
        Arr2DStdCOut2D(arr2DVBN1, " ", "; ", false);
        std::cout<<std::endl;
    }
    std::cout<<"let us go further..."<<std::endl;
    //works well for rect arrays FromN=3 and 1, ShiftQ=3, 0, -3, QRecs=0 and 7  (which is same) => last row whole of X2
    //and works well if nothing else
    //

    std::cout<<"Now let us try concat adding..."<<std::endl;
    Array2DConcatAdding(arr2DVBN1, arr2DVBN2, shiftQ, DfltInt, rect, QRec, FromN);
    std::cout<<"Concatenated 2D arr (adding):"<<std::endl;
    //Arr2DStdCOut2D(arr2DN2Ptr, size2, " ", "; ", false);
    Arr2DStdCOut2D(arr2DVBN1, " ", "; ", false);
    std::cout<<std::endl;
    //works well  - tested for rect arrays,
    //rect for em n'vikts, shiftQ = 0 et 3, FromN = 1 et 3, QRec=0 et 2 //wa errs, fa corrs
    //

    /*
    //
    std::cout<<"Cutting to Rect:"<<std::endl;
    Array2DSetRect(arr2DN1Ptr, size1, 4);
    Arr2DStdCOut2D(arr2DN1Ptr, size1, " ", "; ", false);
    std::cout<<std::endl;
    std::cout<<"Stretching to Rect:"<<std::endl;
    //template<typename T> void Array2DSetRect(T**&X, Array2DSize&size, int QIneRows, int QExtRows=0, T*DefaultValParam=NULL, bool PreserveVals=true){
    //Array2DSetRect(arr2DN1Ptr, size1, 8);
    Array2DSetRect(arr2DN1Ptr, size1, 8, 0, DfltInt, true);
    Arr2DStdCOut2D(arr2DN1Ptr, size1, " ", "; ", false);
    std::cout<<std::endl;
    //works well  - tested for rect arrays
    */
    //Collecting  garbage final
     delete DfltInt;

}//TestArr2DLib


void TestKynemParams(){
    QString str;
    int DoFs[]={0,1,1,0,0,0};
    //int DoFs[]={1,1,0,0,0,0};
    double val;
    std::cout<<"TestKynemParams starts working"<<std::endl;
    //KynematicParams_vb kyn(NULL);//gut
    //KynematicParams_vb kyn();//so S ne int ze kyn s'member o'KynematicParams_vb
    KynematicParams_vb kyn;//gut
    //kyn.Set(1,1,0,0,0,0, 2,1);
    val=111;
    std::cout<<"x[1]="<<val<<std::endl;
    kyn.set_x(val);
    val=kyn.get_x(1);
    std::cout<<"x[1]="<<val<<std::endl;
    val=121;
    std::cout<<"vx[1]="<<val<<std::endl;
    kyn.set_vx(val);
    val=kyn.get_vx();
    std::cout<<"vx[1]="<<val<<std::endl;
    val=131;
    std::cout<<"ax[1]="<<val<<std::endl;
    kyn.set_ax(val);
    val=kyn.get_ax();
    std::cout<<"ax[1]="<<val<<std::endl;
    //
    std::cout<<"Kynematic params: after Vals of all derivs of X DoF of 1st point is given (ma sha nur 0th deriv)"<<std::endl;
    kyn.ShowConsole();
    //
    std::cout<<"Kynematic params: DoFs are given: x - no (but was), y - yes (and was), z - yes (but was not)"<<std::endl;
    std::cout<<"Kynematic params: "<<std::endl;
    kyn.SetParams(DoFs, 2);
    std::cout<<"Kynematic params: "<<std::endl;
    kyn.ShowConsole();
    //
    val=121;
    std::cout<<"y[1]="<<val<<std::endl;
    kyn.set_y(val);
    val=kyn.get_y(1);
    std::cout<<"y[1]="<<val<<std::endl;
    val=151;
    std::cout<<"vy[1]="<<val<<std::endl;
    kyn.set_vy(val);
    val=kyn.get_vy();
    std::cout<<"vy[1]="<<val<<std::endl;
    val=181;
    std::cout<<"ay[1]="<<val<<std::endl;
    kyn.set_ay(val);
    val=kyn.get_ay();
    std::cout<<"ay[1]="<<val<<std::endl;
    std::cout<<"Kynematic params: "<<std::endl;
    kyn.ShowConsole();
    //
    val=131;
    std::cout<<"z[1]="<<val<<std::endl;
    kyn.set_z(val);
    val=kyn.get_z(1);
    std::cout<<"z[1]="<<val<<std::endl;
    val=161;
    std::cout<<"vz[1]="<<val<<std::endl;
    kyn.set_vz(val);
    val=kyn.get_vz();
    std::cout<<"vz[1]="<<val<<std::endl;
    val=191;
    std::cout<<"az[1]="<<val<<std::endl;
    kyn.set_az(val);
    val=kyn.get_az();
    std::cout<<"az[1]="<<val<<std::endl;
    std::cout<<"Kynematic params: "<<std::endl;
    kyn.ShowConsole();
    //
    std::vector<double>row;
    row=kyn.GetPointN(1);
    std::cout<<"Point 1:"<<std::endl;
    Array1DShowConsole(row);
    row=kyn.Get_all_ax();
    std::cout<<"all ax:"<<std::endl;
    Array1DShowConsole(row);
    std::cout<<"or all ax:"<<std::endl;
    str=kyn.get_all_x_as_String();
    std::cout<<"or all ax:"<<std::endl;
    std::cout<<str.toStdString()<<std::endl;
    //
    std::cout<<"Point 1:"<<std::endl;
    Array1DShowConsole(row);
    row=kyn.Get_all_ay();
    std::cout<<"all ay:"<<std::endl;
    Array1DShowConsole(row);
    std::cout<<"or all ay:"<<std::endl;
    str=kyn.get_all_ay_as_String();
    std::cout<<str.toStdString()<<std::endl;
    //
    std::cout<<"Kynematic params: (excluded X, added y and z, no more else):"<<std::endl;
    kyn.ShowConsole();
    std::cout<<"Adding DoF FiX"<<std::endl;
    kyn.AddDoFFiX();
    kyn.set_fix(122);
    kyn.set_wx(152);;
    kyn.set_ex(162);
    std::cout<<"Kynematic params: "<<std::endl;
    kyn.ShowConsole();
    std::cout<<"Excluding DoF Y"<<std::endl;
    kyn.ExcludeDoFY();
    std::cout<<"Kynematic params: "<<std::endl;
    kyn.ShowConsole();
    //
    KynematicCharacteristics1 kynChars1;
    kynChars1.Set(1,0,1,1,0,0,2,1);
    KynematicParams_vb kyn1;
    kyn1=kyn.GetSubStr(1,&kynChars1);
    std::cout<<"Kynematic params of deriv -  with X, and No Y:"<<std::endl;
    kyn1.ShowConsole();
    std::cout<<"TestKynemParams finishes working"<<std::endl;
    double point1[]={131, 231, 331, 431, 531, 631, 731, 831, 931},
           point2[]={132, 232, 332, 432, 532, 632, 732, 832, 932},
           point3[]={133, 233, 333, 433, 533, 633, 733, 833, 933},
           point4[]={134, 234, 334, 434, 534, 634, 734, 834, 934};
    std::cout<<"Trying to add point:"<<std::endl;
    kyn1.AddPoint(point1);
    std::cout<<"Kynematic params af adding point:"<<std::endl;
    kyn1.ShowConsole();
    std::cout<<"Trying to ins point 2:"<<std::endl;
    kyn1.InsPoint(2, point2);
    std::cout<<"Kynematic params af ins point to N2:"<<std::endl;
    kyn1.ShowConsole();
    std::cout<<"Trying to add point:"<<std::endl;
    kyn1.AddPoint(point3);
    std::cout<<"Kynematic params af adding point:"<<std::endl;
    kyn1.ShowConsole();
    std::cout<<"Trying to ins point to 1:"<<std::endl;
    kyn1.InsPoint(1, point4);
    std::cout<<"Kynematic params af ins point to 1:"<<std::endl;
    kyn1.ShowConsole();
    std::cout<<"Trying to del point 4 (ne last):"<<std::endl;
    kyn1.DelPoint(4);
    std::cout<<"Kynematic params af del  point 4 (ne last):"<<std::endl;
    kyn1.ShowConsole();
    std::cout<<"Trying to del point 4 (last):"<<std::endl;
    kyn1.DelPoint(4);
    std::cout<<"Kynematic params af del  point 4  last):"<<std::endl;
    kyn1.ShowConsole();
    std::cout<<"Trying to del point 1:"<<std::endl;
    kyn1.DelPoint(1);
    std::cout<<"Kynematic params af del  point 1:"<<std::endl;
    kyn1.ShowConsole();
}

void TestMath(){
    PolyEq eq;
    double roots4[]={1, 2, 3, 4};
    std::string equation, solution;
    std::vector<double>vroots;
    for(int i=1; i<=4; i++){
        vroots.push_back(roots4[i-1]);
    }
    eq.StartGenerator(vroots);
    equation=eq.EquationToStdString();
    solution=eq.SolutionToToStdString();
    printf("Equation: %s\n",equation.c_str());
    printf("Solution: %s\n",solution.c_str());
    std::cout<<"Equation: "<<equation<<std::endl;
    std::cout<<"Solution: "<<solution<<std::endl;

}
void TestMiscel(){
    QString s0="12345";
    QString ss=MySubString_L1(s0, 3);
    std::cout<<"was: "<<s0.toStdString()<<" 3th char is: "<<ss.toStdString()<<std::endl;
    QString qsarr[]={"Chapter One", "Chapter Two", "Chapter Three", "Chapter Four", "Chapter Five", "Chapter Six", "Chapter Seven", "Chapter Eight", "Chapter Nine"};

}
void TestMatrices(){
    double Coords1[]={50, 30, 70},
           EulerAnglesDeg[]={10, 30, 70},
           OriginCoords[]={0, 0, 100};
    TValsShowHide vsh;
    vsh.EnableWriting();
    vsh.ConsoleInterface=true;
    //
    Matrix MCoords1(Coords1, 3, 1), MEulerAnglesDegS(EulerAnglesDeg, 3, 1), MOriginCoords(OriginCoords, 3, 1), MEulerAnglesRad;
    Matrix MCoords2, MCoords3;
    Matrix MEulerAnglesDeg;
    //
    std::cout<<"Matrices task"<<std::endl;
    //
    std::cout<<"Initial Data:"<<std::endl;
    //
    std::cout<<"Initial coords:"<<std::endl;
    std::cout<<"Matrix:"<<std::endl;
    MCoords1.StdCOut2D();
    //
    std::cout<<"Euler angles (deg):"<<std::endl;
    //
    std::cout<<"Matrix:"<<std::endl;
    MEulerAnglesDeg.StdCOut2D();
    //
    std::cout<<"Origin coords):"<<std::endl;
    std::cout<<"Matrix:"<<std::endl;
    MOriginCoords.StdCOut2D();
    //
    //
    std::cout<<std::endl;
    //
    //
    std::cout<<"Calculation: multiply matrix to number"<<std::endl;
    MEulerAnglesRad = MEulerAnglesDeg *(M_PI/180);
    MEulerAnglesRad.StdCOut2D();
    //
    //
    std::cout<<std::endl;
    //
    //
    std::cout<<"Calculation: multiply matrix to number"<<std::endl;
    std::cout<<"Matrices type: Matrix_PS1"<<std::endl;
    MEulerAnglesRad = MEulerAnglesDeg *(M_PI/180);
    std::cout<<"Euler angles in rad:"<<std::endl;
    MEulerAnglesRad.StdCOut2D();
    //
    //
    std::cout<<std::endl;
    //
    //
    std::cout<<"Calculation:Coord transforms: add, subtract, multiply and divide matrices"<<std::endl;
    MCoords2 = CalcCoordsTransformFromOldCSToNew(&MCoords1, &MEulerAnglesRad, &MOriginCoords, &vsh);
    std::cout<<"Answer: Coords in new Sys:"<<std::endl;
    MCoords2.StdCOut2D();
    std::cout<<std::endl;
    //

    std::cout<<"finish with matrices"<<std::endl;
}


/*
void TestMatrixSimplePS(){
    double Coords1[]={50, 30, 70},
           EulerAnglesDeg[]={10, 30, 70},
           OriginCoords[]={0, 0, 100};
    Matrix_P2S MCoords1_PS(Coords1, 3, 1), MEulerAnglesDeg_PS(EulerAnglesDeg, 3, 1), MOriginCoords_PS(OriginCoords, 3, 1), MEulerAnglesRad_PS;
    Matrix_P2S MCoords2_PS, MCoords3_PS;

    //
    std::cout<<"Matrices task"<<std::endl;
    //
    std::cout<<"Initial Data:"<<std::endl;
    //
    std::cout<<"Initial coords:"<<std::endl;
    std::cout<<"PS type:"<<std::endl;
    MCoords1_PS.StdCOut2D();

    //
    std::cout<<"Euler angles (deg):"<<std::endl;
    std::cout<<"PS type:"<<std::endl;
    MEulerAnglesDeg_PS.StdCOut2D();

    //
    std::cout<<"Origin coords):"<<std::endl;
    std::cout<<"PS type:"<<std::endl;
    MOriginCoords_PS.StdCOut2D();

    //
    //
    std::cout<<"Calculation: multiply matrix to number"<<std::endl;
    //
    //

    std::cout<<"Matrices type: Matrix_PS"<<std::endl;
    MEulerAnglesRad_PS = MEulerAnglesDeg_PS *(M_PI/180);
    std::cout<<"Euler angles in rad:"<<std::endl;
    MEulerAnglesRad_PS.StdCOut2D();

    //
    //
    std::cout<<std::endl;
    std::cout<<"Calculation:Coord transforms: add, subtract, multiply and divide matrices"<<std::endl;
    std::cout<<"Matrices type: Matrix_PS"<<std::endl;
    MCoords2_PS = MCoords2_PS.CalcCoordsTransformFromOldCSToNew(&MCoords1_PS, &MEulerAnglesRad_PS, &MOriginCoords_PS);
    std::cout<<"Answer: Coords in new Sys:"<<std::endl;
    MCoords2_PS.StdCOut2D();
    std::cout<<std::endl;
    //

    //


    std::cout<<"finish with matrices"<<std::endl;
}
*/

/*

void TestMatrixSimplePS1(){
    double Coords1[]={50, 30, 70},
           EulerAnglesDeg[]={10, 30, 70},
           OriginCoords[]={0, 0, 100};
    Matrix_P2S_v1 MCoords1_PS1(Coords1, 3, 1), MEulerAnglesDeg_PS1(EulerAnglesDeg, 3, 1), MOriginCoords_PS1(OriginCoords, 3, 1), MEulerAnglesRad_PS1;
    Matrix_P2S_v1 MCoords2_PS1, MCoords3_PS1;

    //
    std::cout<<"Matrices task"<<std::endl;
    //
    std::cout<<"Initial Data:"<<std::endl;
    //
    std::cout<<"Initial coords:"<<std::endl;
    std::cout<<"PS1 type:"<<std::endl;
    MCoords1_PS1.StdCOut2D();

    //
    std::cout<<"Euler angles (deg):"<<std::endl;
    std::cout<<"PS1 type:"<<std::endl;
    MEulerAnglesDeg_PS1.StdCOut2D();

    //
    std::cout<<"Origin coords):"<<std::endl;
    std::cout<<"P1S type:"<<std::endl;
    MOriginCoords_PS1.StdCOut2D();

    //
    //
    std::cout<<"Calculation: multiply matrix to number"<<std::endl;
    //
    //

    std::cout<<"Matrices type: Matrix_PS1"<<std::endl;
    MEulerAnglesRad_PS1 = MEulerAnglesDeg_PS1 *(M_PI/180);
    std::cout<<"Euler angles in rad:"<<std::endl;
    MEulerAnglesRad_PS1.StdCOut2D();

    //
    //
    std::cout<<std::endl;
    //
    //

    //


    std::cout<<"finish with matrices"<<std::endl;
}
*/


/*void TestMatrixSimpleA2P(){
    double Coords1[]={50, 30, 70},
           EulerAnglesDeg[]={10, 30, 70},
           OriginCoords[]={0, 0, 100};
    Matrix_P2A MCoords1_A2P(Coords1, 3, 1), MEulerAnglesDeg_A2P(EulerAnglesDeg, 3, 1), MOriginCoords_A2P(OriginCoords, 3, 1), MEulerAnglesRad_A2P;
    Matrix_P2A MCoords2_A2P, MCoords3_A2P;

    //
    std::cout<<"Matrices task"<<std::endl;
    //
    std::cout<<"Initial Data:"<<std::endl;
    //
    std::cout<<"Initial coords:"<<std::endl;
    std::cout<<"A2P type:"<<std::endl;
    MCoords1_A2P.StdCOut2D();

    //
    std::cout<<"Euler angles (deg):"<<std::endl;
    std::cout<<"A2P type:"<<std::endl;
    MEulerAnglesDeg_A2P.StdCOut2D();

    //
    std::cout<<"Origin coords):"<<std::endl;
    std::cout<<"A2P type:"<<std::endl;
    MOriginCoords_A2P.StdCOut2D();

    //
    //
    std::cout<<"Calculation: multiply matrix to number"<<std::endl;
    //
    //

    std::cout<<"Matrices type: Matrix_A2P"<<std::endl;
    MEulerAnglesRad_A2P = MEulerAnglesDeg_A2P *(M_PI/180);
    std::cout<<"Euler angles in rad:"<<std::endl;
    MEulerAnglesRad_A2P.StdCOut2D();

    //
    //
    std::cout<<std::endl;
    std::cout<<"Calculation:Coord transforms: add, subtract, multiply and divide matrices"<<std::endl;
    std::cout<<"Matrices type: Matrix_A2P"<<std::endl;
    MCoords2_A2P = MCoords2_A2P.CalcCoordsTransformFromOldCSToNew(&MCoords1_A2P, &MEulerAnglesRad_A2P, &MOriginCoords_A2P);
    std::cout<<"Answer: Coords in new Sys:"<<std::endl;
    MCoords2_A2P.StdCOut2D();
    std::cout<<std::endl;
    //
    //

    //


    std::cout<<"finish with matrices"<<std::endl;
}
*/


void TestLinearApprox(){
    TValsShowHide vsh;
    vsh.ConsoleInterface=true;
    vsh.Show1Hide0=1;
    //            1      2       3   4   5  6        8     9   10     11   12  13  14  15  16
    double X[]={-15,    -10,   -5, 0,   5, 10,   15, 20,   25, 30,   35, 40,   45, 50, 55, 60};
    double Y[]={-22.5,  -15, -7.5, 0, 7.5, 15, 22.5, 30, 37.5, 45, 52.5, 60, 67.5, 70, 65, 60};
    double a, b, ya;
    int Q=16;
    LinearApproximation(X, Y, Q, a, b);
    writeln(&vsh, "a="+FloatToStr(a)+" b="+FloatToStr(b));
    for(int i=1; i<=Q; i++){
        ya=a*X[i-1]+b;
        writeln(&vsh,"X="+FloatToStr(X[i-1])+" Yp="+FloatToStr(Y[i-1])+" Ya="+FloatToStr(ya));
    }
}

void TestPtrNotFirst(){
    int arr[10];int shiftN=2;
    int Q=3;
    int*ptr1=(arr+shiftN-1), *ptr2;
    for(int i=1; i<=10; i++){
        arr[i-1]=i;
    }
    for(int i=1; i<=Q; i++){
        ptr2=(ptr1+i-1);
        printf("%d\n",ptr2);
    }
}
void TestSeidel(){
    TItersPrecision Precision;
    Matrix M(3, 3), V1(3, 1), V2(3, 1), X0(3,1), XT(3, 1), XS, MI, MT;
    M.SetComponent(11, 1, 1);
    M.SetComponent(12, 1, 2);
    M.SetComponent(13, 1, 3);
    M.SetComponent(21, 2, 1);
    M.SetComponent(22, 2, 2);
    M.SetComponent(23, 2, 3);
    M.SetComponent(31, 3, 1);
    M.SetComponent(32, 3, 2);
    M.SetComponent(33, 3, 3);
    V1.SetComponent(74, 1, 1);
    V1.SetComponent(134, 1, 2);
    V1.SetComponent(194, 1, 3);
    V2.SetComponent(74, 1, 1);
    V2.SetComponent(134, 1, 2);
    V2.SetComponent(190, 1, 3);
    for(int i=1; i<=3; i++)XT.SetComponent(i, i, 1);
    Precision.MaxQIters=10;
    Precision.precision=1E-6;
    Precision.Priority_1QI2Eps3And4Or=4;
    std::cout<<"Haben matrices:"<<std::endl;
    std::cout<<"M :"<<std::endl;
    M.StdCOut2D();
    std::cout<<"Testing multiplication:"<<std::endl;
    XS=M*XT;
    std::cout<<"Result:"<<std::endl;
    XS.StdCOut2D();
    std::cout<<"M again:"<<std::endl;
    M.StdCOut2D();
    M.SetComponent(51, 1, 1);
    M.SetComponent(52, 2, 2);
    M.SetComponent(53, 3, 3);
    std::cout<<"M with new vals:"<<std::endl;
    M.StdCOut2D();
    std::cout<<"M inverted:"<<std::endl;
    MI=M.GetInverted();
    MI.StdCOut2D();
    std::cout<<"M transposed:"<<std::endl;
    MT=M.GetTransposed();
    MT.StdCOut2D();
    std::cout<<"Matrix precise method solution:"<<std::endl;
    //XS=V1/M;
    MI=M.GetInverted();
    XS=MI*V1;
    std::cout<<"Solution:"<<std::endl;
    XS.StdCOut2D();
    std::cout<<"Seidel approximate method solution:"<<std::endl;
    //Matrix_VL fApproxSolveLinAlgEqsSysSeidel(const Matrix_VL&M, const Matrix_VL&V, Matrix_VL X0, TItersPrecision Precision){
    XS=fApproxSolveLinAlgEqsSysSeidel(M, V2, X0, Precision);
    XS.StdCOut2D();
}

void TestMatrix(){
    //
    int MatrixTypeN=MatrixTypeN_V2A;
    //MatrixTypeN=MatrixTypeN_V2L;
    TItersPrecision Precision;
    Matrix M(3, 3),
           V1(3, 1),
           V2(3, 1),
            X0(3,1),
            XT(3, 1),
            XS,
            MI,
            MT;
        M.SetComponent(11, 1, 1);
        M.SetComponent(12, 1, 2);
        M.SetComponent(13, 1, 3);
        M.SetComponent(21, 2, 1);
        M.SetComponent(22, 2, 2);
        M.SetComponent(23, 2, 3);
        M.SetComponent(31, 3, 1);
        M.SetComponent(32, 3, 2);
        M.SetComponent(33, 3, 3);
        V1.SetComponent(74, 1, 1);
        V1.SetComponent(134, 1, 2);
        V1.SetComponent(194, 1, 3);
        V2.SetComponent(74, 1, 1);
        V2.SetComponent(134, 1, 2);
        V2.SetComponent(190, 1, 3);
        for(int i=1; i<=3; i++)XT.SetComponent(i, i, 1);
        Precision.MaxQIters=10;
        Precision.precision=1E-6;
        Precision.Priority_1QI2Eps3And4Or=4;
        std::cout<<"Haben matrices:"<<std::endl;
        std::cout<<"M :"<<std::endl;
        M.StdCOut2D();
        std::cout<<"Testing multiplication:"<<std::endl;
        XS=M*XT;
        std::cout<<"Result:"<<std::endl;
        XS.StdCOut2D();
        std::cout<<"M again:"<<std::endl;
        M.StdCOut2D();
        M.SetComponent(51, 1, 1);
        M.SetComponent(52, 2, 2);
        M.SetComponent(53, 3, 3);
        std::cout<<"M with new vals:"<<std::endl;
        M.StdCOut2D();
        std::cout<<"M inverted:"<<std::endl;
        MI=M.GetInverted();
        MI.StdCOut2D();
        std::cout<<"M transposed:"<<std::endl;
        MT=M.GetTransposed();
        MT.StdCOut2D();
        std::cout<<"Matrix precise method solution:"<<std::endl;
        //XS=V1/M;
        MI=M.GetInverted();
        XS=MI*V1;
        std::cout<<"Solution:"<<std::endl;
        XS.StdCOut2D();
        std::cout<<"Seidel approximate method solution:"<<std::endl;
        //Matrix_VL fApproxSolveLinAlgEqsSysSeidel(const Matrix_VL&M, const Matrix_VL&V, Matrix_VL X0, TItersPrecision Precision){
        XS=fApproxSolveLinAlgEqsSysSeidel(M, V2, X0, Precision);
        XS.StdCOut2D();


    double Coords1[]={50, 30, 70},
           EulerAnglesDeg[]={10, 30, 70},
           OriginCoords[]={0, 0, 100};
    Matrix MCoords1(Coords1, 3, 1), MEulerAnglesDeg(EulerAnglesDeg, 3, 1), MOriginCoords(OriginCoords, 3, 1), MEulerAnglesRad;
    Matrix MCoords2, MCoords3;
    //
    std::cout<<"Matrices task"<<std::endl;
    //
    std::cout<<"Initial Data:"<<std::endl;
    //
    std::cout<<"Initial coords:"<<std::endl;
    std::cout<<"A2P type:"<<std::endl;
    MCoords1.StdCOut2D();
    //
    std::cout<<"Euler angles (deg):"<<std::endl;
    std::cout<<"A2P type:"<<std::endl;
    MEulerAnglesDeg.StdCOut2D();
    //
    std::cout<<"Origin coords):"<<std::endl;
    std::cout<<"A2P type:"<<std::endl;
    MOriginCoords.StdCOut2D();
    //
    //
    std::cout<<"Calculation: multiply matrix to number"<<std::endl;
    //
    //
    std::cout<<"Matrices type: Matrix_A2P"<<std::endl;
    MEulerAnglesRad = MEulerAnglesDeg *(M_PI/180);
    std::cout<<"Euler angles in rad:"<<std::endl;
    MEulerAnglesRad.StdCOut2D();
    //
    //
    std::cout<<std::endl;
    std::cout<<"Calculation:Coord transforms: add, subtract, multiply and divide matrices"<<std::endl;
    std::cout<<"Matrices type: Matrix_A2P"<<std::endl;
    //Matrix CalcCoordsTransformFromOldCSToNew(Matrix*CoordsInOldCSParam, Matrix*EulerAnglesParam, Matrix*CSOriginOfNewCSInOldParam=NULL, TValsShowHide*vsh=NULL);
    //Matrix CalcCoordsTransformFromNewCSToOld(Matrix*CoordsInNewCSParam, Matrix*EulerAnglesParam, Matrix*CSOriginOfNewCSInOldParam=NULL, TValsShowHide*vsh=NULL);
    MCoords2 = CalcCoordsTransformFromOldCSToNew(&MCoords1, &MEulerAnglesRad, &MOriginCoords);
    std::cout<<"Answer: Coords in new Sys:"<<std::endl;
    MCoords2.StdCOut2D();
    std::cout<<std::endl;
    //
    //

    //


    std::cout<<"finish with matrices"<<std::endl;
}
