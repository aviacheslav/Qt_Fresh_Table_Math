#include "testmatrix.h"


// class Matrix_PS1{
//     double**y;
//     int QLins, QCols;
//     public:
Matrix_PS1::Matrix_PS1(int QLins, int QCols){
 std::cout<<"Matrix_PS constr (QL, QC) works"<<std::endl;
 this->Construct();
 this->SetSize(QLins, QCols);
}
Matrix_PS1::Matrix_PS1(double**y, int QLins, int QCols){
 std::cout<<"Matrix_PS constr (**dbl, QL, QC) works"<<std::endl;
 this->Construct();
 this->Set(y, QLins, QCols);
}
Matrix_PS1::Matrix_PS1(std::vector<std::vector<double>>y){
 std::cout<<"Matrix_PS constr (vector<<dbl>>, QL, QC) works"<<std::endl;
 this->Construct();
 this->Set(y);
}

Matrix_PS1::Matrix_PS1(double*y, int Q, bool ColOfVectorNotHorLine){
 std::cout<<"Matrix_PS constr (*dbl, Q, bool) works"<<std::endl;
 this->Construct();
 this->Set(y, Q, ColOfVectorNotHorLine);
}
Matrix_PS1::Matrix_PS1(std::vector<double>y, bool ColOfVectorNotHorLine){
 std::cout<<"Matrix_PS constr (vector<dbl>,  bool) works"<<std::endl;
 this->Construct();
 this->Set(y, ColOfVectorNotHorLine);
}
Matrix_PS1::Matrix_PS1(double x, double y, double z){
 std::cout<<"Matrix_PS constr (x, y, z) works"<<std::endl;
 this->Construct();
 this->Set(x, y, z);
}
Matrix_PS1::Matrix_PS1(const Matrix_PS1&obj){
 std::cout<<"Matrix_PS COPY constr works"<<std::endl;
 this->Construct();
 this->Assign(obj);
}
void Matrix_PS1::Construct(){ this->SetNullIni(); }
void Matrix_PS1::SetNull(){
 if(this->y!=NULL){
     for(int i=1; i<=this->QLins; i++){
         delete[] this->y[i-1];
     }
     delete[] this->y;
 }
 this->QLins=0;
 this->QCols=0;
}
void Matrix_PS1::SetNullIni(){ this->y=NULL; this->QLins=0; this->QCols=0; }
Matrix_PS1::~Matrix_PS1(){
 std::cout<<"Matrix_PS Destructor works"<<std::endl;
 this->SetNull();
}
void Matrix_PS1::Assign(const Matrix_PS1&obj){
 this->SetNull();
 this->Set(obj.y, obj.QLins, obj.QCols);
}
Matrix_PS1& Matrix_PS1::operator = (const Matrix_PS1&obj){
 this->Assign(obj);
 return*this;
}
void Matrix_PS1::SetSize(int QLins, int QCols){
 double**z=NULL;
 int QLinsMin, QColsMin;
 if(this->QLins!=QLins || this->QCols!=QCols){
     z=new double*[QLins];
     for(int i=1; i<=QLins; i++){
         z[i-1]=new double[QCols];
     }
     //
     if(this->y!=NULL){
         QLinsMin=QLins<=this->QLins ? QLins : this->QLins;
         QColsMin=QLins<=this->QCols ? QCols : this->QCols;
         for(int i=1; i<=QLinsMin; i++){
             for(int j=1; j<=QColsMin; j++){
                 z[i-1][j-1]=y[i-1][j-1];
             }
             for(int j=QColsMin+1; j<=QCols; j++){
                 z[i-1][j-1]=0;
             }
         }
         for(int i=QLinsMin+1; i<=QLins; i++){
             for(int j=1; j<=QCols; j++){
                 z[i-1][j-1]=0;
             }
         }
     }
     //
     if(this->y!=NULL){
         for(int i=1; i<=this->QLins; i++){
             delete[] this->y[i-1];
         }
         delete[] this->y;
     }
     //
     this->y=z;
     //
     this->QLins=QLins;
     this->QCols=QCols;
 }
}
void Matrix_PS1::Set(double**y, int QLins, int QCols){
 //this->SetNull();
 this->SetSize(QLins, QCols);
 if(y!=NULL){
     for(int i=1; i<=QLins; i++){
         for(int j=1; j<=QCols; j++){
             this->y[i-1][j-1]=y[i-1][j-1];
         }
     }
 }else{
     for(int i=1; i<=QLins; i++){
         for(int j=1; j<=QCols; j++){
             this->y[i-1][j-1]=0;
         }
     }
 }
}
void Matrix_PS1::Set(std::vector<std::vector<double>>y){
 //this->SetNull();
 int QLins=y.size(), QCols=y[1-1].size();
 this->SetSize(QLins, QCols);
 for(int i=1; i<=QLins; i++){
     for(int j=1; j<=QCols; j++){
         this->y[i-1][j-1]=y[i-1][j-1];
     }
 }
}
void Matrix_PS1::Set(double*y, int Q, bool ColOfVectorNotHorLine){
 if(ColOfVectorNotHorLine){
     this->SetSize(Q, 1);
     for(int i=1; i<=Q; i++){
         for(int j=1; j<=1; j++){
             this->y[i-1][j-1]=y[i-1];
         }
     }
     this->QLins=Q;
     this->QCols=1;
 }else{
     this->SetSize(1, Q);
     for(int i=1; i<=1; i++){
         for(int j=1; j<=Q; j++){
             this->y[i-1][j-1]=y[j-1];
         }
     }
     this->QLins=Q;
     this->QCols=1;
 }
}
void  Matrix_PS1::Set(std::vector<double>y, bool ColOfVectorNotHorLine){
 this->Set(y.data(), y.size(), ColOfVectorNotHorLine);
}
void Matrix_PS1::Set(double x, double y, double z){
  this->SetNull();
  this->SetSize(3, 1);
   this->y[1-1][1-1]=x;
   this->y[2-1][1-1]=y;
   this->y[3-1][1-1]=z;
}
int Matrix_PS1::GetQLines() const { return this->QLins; }
int Matrix_PS1::GetQColumns() const { return this->QCols; }
double Matrix_PS1::GetComponent(int LineN, int ColN) const{
 double y=0;
 if(LineN>=1 && LineN<=this->QLins && ColN>=1 && ColN<=this->QCols){
     y=this->y[LineN-1][ColN-1];
 }
 return y;
}
void Matrix_PS1::SetComponent(int LineN, int ColN, double val){
 if(LineN>=1 && LineN<=this->QLins && ColN>=1 && ColN<=this->QCols){
     this->y[LineN-1][ColN-1]=val;
 }
}

Matrix_PS1 Matrix_PS1::operator *(double x){
    double cur;
    for(int i=1; i<=this->QLins; i++){
        for(int j=1; j<=this->QCols; j++){
           cur=this->y[i-1][j-1]*x;
           this->y[i-1][j-1]=cur;
        }
    }
    return *this;
}


QString Matrix_PS1::LineToString(int LineN, QString delim, bool useBrackets){
    int QC=this->GetQColumns();
    QString s="", sn;
    if(LineN>=1 && LineN<=this->GetQLines()){
        if(useBrackets){
            s=s+"[";
        }
        for(int i=1; i<=QC-1; i++){
            sn.setNum(this->GetComponent(LineN, i));
            s=s+sn;
            s=s+delim;
        }
        sn.setNum(this->GetComponent(LineN, QC));
        s=s+sn;
        if(useBrackets){
            s=s+"]";
        }
    }
    return s;
}





void Matrix_PS1::StdCOut2D(QString delimElem, QString delimRow, bool useBrackets){
    int QL=this->GetQLines();
    std::string s, c;
    QString qs;
    s="";
    if(useBrackets){
        s="[";
    }
    qs=this->LineToString(1, delimElem, useBrackets);
    c=qs.toStdString();
    s=s+c;
    if(QL>1){
        s=s+delimRow.toStdString();
    }
    std::cout<<s<<std::endl;
    for(int i=2; i<=QL-1; i++){
        if(useBrackets){
            s=" ";
        }
        qs=this->LineToString(i, delimElem, useBrackets);
        c=qs.toStdString();
        s=s+c;
        s=s+delimRow.toStdString();
        std::cout<<s<<std::endl;
    }
    if(useBrackets){
        s=" ";
    }
    qs=LineToString(QL, delimElem, useBrackets);
    c=qs.toStdString();
    s=s+c;
    if(useBrackets){
        s+="]";
    }
    std::cout<<s<<std::endl;
}//fn
