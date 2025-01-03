#include "mymathlib.h"

//MyMathLib::MyMathLib()


TPosInSucc PosInSucc1(double*X, int Q, double x, double eps){
    TPosInSucc pos;
    MathApprox appr(eps);
    pos.isLess=false;
    pos.isGreater=false;
    pos.isWithin-false;
    pos.EqualN=0;
    pos.lessN=0;
    //
    if(x<X[1-1]){
        pos.isLess=true;
    }else if(x>X[Q-1]){
       pos.isGreater=true;
    }else{
        pos.isWithin=true;
        for(int i=1; i<=Q; i++){
            if(appr.ET(X[i-1], x)){
                pos.EqualN=i;
                break;
            }
        }
        if(pos.EqualN==0){
            for(int i=1; i<=Q-1; i++){
                if(appr.LT(X[i-1], x) && appr.GT(X[i+1-1], x)){
                    pos.lessN=i;
                    break;
                }
            }
        }
    }
    return pos;
}

TPosInSucc PosInSucc(std::vector<double>X, double x, double eps){
    TPosInSucc pos;
    MathApprox appr(eps);
    int Q=X.size();
    pos.isLess=false;
    pos.isGreater=false;
    pos.isWithin-false;
    pos.EqualN=0;
    pos.lessN=0;
    //
    if(x<X[1-1]){
        pos.isLess=true;
    }else if(x>X[Q-1]){
       pos.isGreater=true;
    }else{
        pos.isWithin=true;
        for(int i=1; i<=Q; i++){
            if(appr.ET(X[i-1], x)){
                pos.EqualN=i;
                break;
            }
        }
        if(pos.EqualN==0){
            for(int i=1; i<=Q-1; i++){
                if(appr.LT(X[i-1], x) && appr.GT(X[i+1-1], x)){
                    pos.lessN=i;
                    break;
                }
            }
        }
    }
    return pos;
}

double LinInterp(double*X, double*Y, int Q, double x, double eps){
    TPosInSucc pos=PosInSucc1(X, Q, x, eps);
    double y=0, y1, y2, x1, x2, k;
    if(pos.EqualN>0){
        y=Y[pos.EqualN-1];
    }else{
        if(pos.isLess){
            x1=X[1-1];
            x2=X[2-1];
            y1=Y[1-1];
            y2=Y[2-1];
        }else if(pos.isGreater){
            x1=X[Q-1-1];
            x2=X[Q-1];
            y1=Y[Q-1-1];
            y2=Y[Q-1];
        }else if(pos.lessN>0){
            x1=X[pos.lessN-1];
            x2=X[pos.lessN+1-1];
            y1=X[pos.lessN-1];
            y2=X[pos.lessN+1-1];
        }
        k=(y2-y1)/(x2-x1);
        y=k*(x-x1)+y1;
    }
    return y;
}

double LinInterp(std::vector<double>X, std::vector<double>Y, double x, double eps){
    int Q=X.size();
    TPosInSucc pos=PosInSucc(X, x, eps);
    double y=0, y1, y2, x1, x2, k;
    if(pos.EqualN>0){
        y=Y[pos.EqualN-1];
    }else{
        if(pos.isLess){
            x1=X[1-1];
            x2=X[2-1];
            y1=Y[1-1];
            y2=Y[2-1];
        }else if(pos.isGreater){
            x1=X[Q-1-1];
            x2=X[Q-1];
            y1=Y[Q-1-1];
            y2=Y[Q-1];
        }else if(pos.lessN>0){
            x1=X[pos.lessN-1];
            x2=X[pos.lessN+1-1];
            y1=X[pos.lessN-1];
            y2=X[pos.lessN+1-1];
        }
        k=(y2-y1)/(x2-x1);
        y=k*(x-x1)+y1;
    }
    return y;
}

double LinInterp2D(std::vector<std::vector<double>>Arr2D, std::vector<double>Xext, std::vector<double>Xine, double xe, double xi, double eps){
    TPosInSucc posExt=PosInSucc(Xext, xe, eps), posIne=PosInSucc(Xine, xi, eps);
    double y0, y;
    std::vector<double>X0, X, Y0, Y;
    int Q=Arr2D.size(), L=Arr2D[1-1].size();
    X=Xine;
    if(posExt.EqualN>0){
        Y=Arr2D[posExt.EqualN-1];
    }else{
        if(posExt.isLess){
            X0.clear();
            X0.push_back(Xext[1-1]);
            X0.push_back(Xext[2-1]);
            for(int i=1; i<=L; i++){
                Y0.clear();
                Y0.push_back(Arr2D[1-1][i-1]);
                Y0.push_back(Arr2D[2-1][i-1]);
                //
                y0=LinInterp(X0, Y0, xe, eps);
                Y.push_back(y0);
            }
        }else if(posExt.isGreater){
            X0.clear();
            X0.push_back(Xext[Q-1-1]);
            X0.push_back(Xext[Q-1]);
            for(int i=1; i<=L; i++){
                Y0.clear();
                Y0.push_back(Arr2D[Q-1-1][i-1]);
                Y0.push_back(Arr2D[Q-1][i-1]);
                //
                y0=LinInterp(X0, Y0, xe, eps);
                Y.push_back(y0);
            }
        }else if(posExt.lessN>0){
            X0.clear();
            X0.push_back(Xext[posExt.lessN-1]);
            X0.push_back(Xext[posExt.lessN+1-1]);
            for(int i=1; i<=L; i++){
                Y0.clear();
                Y0.push_back(Arr2D[posExt.lessN-1][i-1]);
                Y0.push_back(Arr2D[posExt.lessN+1-1][i-1]);
                //
                y0=LinInterp(X0, Y0, xe, eps);
                Y.push_back(y0);
            }
        }
        y=LinInterp(X, Y, xi, eps);
    }
    return y;
}

std::vector<double>GetSectionsBounds(double LB, double HB, int Q, bool QOfSectsNotSectsBounds=true){
    std::vector<double>Y;
    int QSects, QSectsBounds;
    double SectL, y;
    QSects = QOfSectsNotSectsBounds==true ? Q : Q-1;
    QSectsBounds=QSects+1;
    SectL=(HB-LB)/QSects;
    Y.push_back(LB);
    for(int SectBoundN=1; SectBoundN<=QSectsBounds; SectBoundN++){
        y=LB+SectL*(SectBoundN-1);
        if(SectBoundN>1 && SectBoundN<QSectsBounds){
            Y.push_back(y);
        }
    }
    Y.push_back(HB);
    return Y;
}


double Integr_Trapez_SimpleByElementsArray(std::vector<double>X, std::vector<double>Y){
    double ys=0, c1, c2, dx;
    int Q=Y.size();
    for(int i=1; i<=Q-1; i++){
        dx=X[i+1-1]-X[i-1];
        c1=(Y[i+1-1]+Y[i-1])/2;
        c2=c1*dx;
        ys+=c2;
    }
    return ys;
}

double Integr_Trapez_SimpleByElementsArray(double*X, double*Y, int Q){
    double ys=0, c1, c2, dx;
    //int Q=Y.size();
    for(int i=1; i<=Q-1; i++){
        dx=X[i+1-1]-X[i-1];
        c1=(Y[i+1-1]+Y[i-1])/2;
        c2=c1*dx;
        ys+=c2;
    }
    return ys;
}

//void CalcDigits(int num, int digits[], /*int&order,*/ int SysBase, bool AsInNum){
void CalcDigits(int num, int digits[], int&order, int SysBase, TValsShowHide*vsh){//, bool AsInNum){
    int digit, z;
    z=num>=0?num:-num;
    //if(z>SysBase)
    order=-1;
    writeln(vsh,"CalcDigits starts working. Try number "+IntToStr(z));
    while(z>=SysBase){
        order++;
        digit=z%SysBase;
        digits[order]=digit;
        writeln(vsh,"z="+IntToStr(z)+" > Sys base="+IntToStr(SysBase)+" => digit="+IntToStr(digit));
        z=(z-digit)/SysBase;
    }
    order++;
    digit=z;
    digits[order]=digit;
    writeln(vsh,"z="+IntToStr(z)+" > Sys base="+IntToStr(SysBase)+" => digit="+IntToStr(digit)+" order="+IntToStr(order));
}

double IntPower(double a, int b){
    double y=1;
    if(b>=0){
        for(int i=1; i<=b; i++){
            y*=a;
        }
    }else{
        for(int i=1; i<=b; i++){
            y/=a;
        }
    }
    return y;
}

int NaturalPower(int a, int b){
   int y=1;
   for(int i=1; i<=b; i++){
       y*=a;
   }
   return y;
}

QString ComplexNumberToString(double ReX, double ImX){
    QString s, sRe, sIm;
    sRe.setNum(ReX);
    sIm.setNum(ImX);
    if(ImX!=0){
        s+="(";
        s+=sRe;
        if(ImX>=0){
            s+=" + ";
        }
        s+=sIm;
        s+="*i)";
    }else{
        s=sRe;
    }
    return s;
}

//class TItersPrecision{
//public:
//    double precision;
//    int MaxQIters;
//    int Priority_1QI2Eps3And4Or;
TItersPrecision::TItersPrecision(double precision, int MaxQIters, int Priority_1QI2Eps3And4Or){
    this->precision=precision;
    this->MaxQIters=MaxQIters;
    this->Priority_1QI2Eps3And4Or=Priority_1QI2Eps3And4Or;
}
bool TItersPrecision::ContinieIterations(double precision, int QIters){
    bool contin=true;
    switch(this->Priority_1QI2Eps3And4Or){
        case 1:
            contin=(QIters>=this->MaxQIters);
        break;
        case 2:
            contin=(precision<=this->precision);
        break;
        case 3:
            contin=(QIters>=this->MaxQIters && QIters>=this->MaxQIters);
        break;
        case 4:
            contin=(QIters>=this->MaxQIters || QIters>=this->MaxQIters);
        break;
    }
    return contin;
}
//};

//class MathApprox{
//public:
//    double precision=1E-9;
    MathApprox::MathApprox(double precision){
        this->precision=precision;
    }
    bool MathApprox::ET(double x1, double x2){
        bool verdict=false;
        if(fabs(x1-x2)<this->precision){
            verdict=true;
        }
        return verdict;
    }
    bool MathApprox::GT(double x1, double x2){
        bool verdict=false;
        if(x1>x2 && fabs(x1-x2)>this->precision){
            verdict=true;
        }
        return verdict;
    }
    bool MathApprox::LT(double x1, double x2){
        bool verdict=false;
        if(x1<x2 && fabs(x2-x1)>this->precision){
            verdict=true;
        }
        return verdict;
    }
    bool MathApprox::NE(double x1, double x2){ return(!(this->ET(x1, x2)));}
    bool MathApprox::GE(double x1, double x2){ return(!(this->LT(x1, x2)));}
    bool MathApprox::LE(double x1, double x2){ return(!(this->GT(x1, x2)));}

    bool MathApprox::isPositive(double x){ return this->GT(x, 0); }
    bool MathApprox::isNegative(double x){ return this->LT(x, 0); }
    bool MathApprox::isZero(double x){ return this->ET(x, 0); }
    bool MathApprox::isNonZero(double x){ return this->NE(x, 0); }

    //

    bool fApprET(double x1, double x2, double precision){
        bool verdict=false;
        if(fabs(x1-x2)<precision){
            verdict=true;
        }
        return verdict;
    }
    bool fApprGT(double x1, double x2, double precision){
        bool verdict=false;
        if(x1>x2 && fabs(x1-x2)>precision){
            verdict=true;
        }
        return verdict;
    }
    bool fApprLT(double x1, double x2, double precision){
        bool verdict=false;
        if(x1<x2 && fabs(x2-x1)>precision){
            verdict=true;
        }
        return verdict;
    }
    bool fApprNE(double x1, double x2, double precision){ return(!(fApprET(x1, x2, precision))); }
    bool fApprGE(double x1, double x2, double precision){ return(!(fApprLT(x1, x2, precision))); }
    bool fApprLE(double x1, double x2, double precision){ return(!(fApprGT(x1, x2, precision))); }
    //
    bool fApprPositive(double x, double precision){ return(!(fApprGT(x, 0, precision))); }
    bool fApprNegative(double x, double precision){ return(!(fApprLT(x, 0, precision))); }
    bool fApprZero(double x, double precision){ return(!(fApprET(x, 0, precision))); }
    bool fApprNonZero(double x, double precision){ return(!(fApprNE(x, 0, precision))); }


    void EquationQuadratic(std::vector<double>C, std::vector<double>&ReX, std::vector<double>&ImX){
        double a, b, c, d, re_x, im_x;
        ReX.clear();
        ImX.clear();
        a=C[2];
        b=C[1];
        c=C[0];
        d=b*b-4*a*c;
        if(a!=0){
            if(d<0){
                re_x=-b/2/a;
                ReX.push_back(re_x);
                ReX.push_back(re_x);
                im_x=sqrt(-d)/2/a;
                ImX.push_back(-im_x);
                ImX.push_back(im_x);
            }
        }
    }
    void EquationQuadratic(double C[], std::vector<double>&ReX, std::vector<double>&ImX){
        double a, b, c;//, d, re_x, im_x;
        ReX.clear();
        ImX.clear();
        a=C[2];
        b=C[1];
        c=C[0];
        EquationQuadratic(a, b, c, ReX, ImX);
    }
    void EquationQuadratic(double BToA, double CToA, std::vector<double>&ReX, std::vector<double>&ImX){
        double a, b, c;//, d, re_x, im_x;
        ReX.clear();
        ImX.clear();
        a=1;
        b=BToA;
        c=CToA;
        EquationQuadratic(a, b, c, ReX, ImX);
    }
    void EquationQuadratic(double a, double b, double c, std::vector<double>&ReX, std::vector<double>&ImX){
        //double a, b, c;
        double d, re_x, im_x;
        ReX.clear();
        ImX.clear();
        //a=C[2];
        //b=C[1];
        //c=C[0];
        d=b*b-4*a*c;
        if(a!=0){
            if(d<0){
                re_x=-b/2/a;
                ReX.push_back(re_x);
                ReX.push_back(re_x);
                im_x=sqrt(-d)/2/a;
                ImX.push_back(-im_x);
                ImX.push_back(im_x);
            }else{
                re_x=(-b-sqrt(d))/2/a;
                ReX.push_back(re_x);
                re_x=(-b+sqrt(d))/2/a;
                ReX.push_back(re_x);
                ImX.push_back(0);
                ImX.push_back(0);
            }
        }
    }
    void EquationCubic(std::vector<double>C, std::vector<double>&ReX, std::vector<double>&ImX){
        double c3, c2, c1, c0;
        c0=C[0];
        c1=C[1];
        c2=C[2];
        c3=C[3];
        EquationCubic(c3, c2, c1, c0, ReX, ImX);
    }
    void EquationCubic(double C[], std::vector<double>&ReX, std::vector<double>&ImX){
        double c3, c2, c1, c0;
        c0=C[0];
        c1=C[1];
        c2=C[2];
        c3=C[3];
        EquationCubic(c3, c2, c1, c0, ReX, ImX);
    }
void EquationCubic(double c3, double c2, double c1, double c0, std::vector<double>&ReX, std::vector<double>&ImX){
    double d, re_x, im_x, Re1, Re2, Im2, Re3, Im3;
    //s=new string[6];
    //s[1-1]="c3="+c3.ToString()+" c2="+c2.ToString()+" c1="+c1.ToString()+" c0="+c0.ToString();
    double a, b, c, A, B, p, q, Q, alfa, CosAlfa, ReY1, ReY2, ReY3;
    a=c2/c3; b=c1/c3; c=c0/c3;
    p=-a*a/3+b;
    q=2*a*a*a/27-a*b/3+c;
    Q=p*p*p/27+q*q/4;
    B = 0;
    CosAlfa = 0;
    alfa = 0;
    A = 0;
    //s[2-1]="p="+p.ToString()+" q="+q.ToString()+" Q="+Q.ToString();
    if(Q<0){
       // s[3-1]="Equation has 3 real not equal roots";
        CosAlfa=-q/(2*sqrt(-(p/3)*(p/3)*(p/3)));
        alfa=acos(CosAlfa);
        ReY1 = 2.0 * sqrt(-(p / 3.0)) * cos(alfa / 3.0);
        ReY2 = 2 * sqrt(-(p / 3)) * cos(alfa / 3 + 2 *M_PI / 3);
        ReY3 = 2 * sqrt(-(p / 3)) * cos(alfa / 3 - 2 *M_PI / 3);
        Im2=0;
        Im3=0;
    }else{
        //if(Q>0)s[3-1]="Equation has 1 real and 2 complex roots";
        //else  s[3-1]="Equation has 3 real roots, min 2 of them are equal";
        if ((-q / 2 + sqrt(Q)) >= 0){
            A = pow((-q / 2 + sqrt(Q)), 1 / 3.0);
        }else{
            A = -pow(-(-q / 2 + sqrt(Q)), 1 / 3.0);
        }
        if ((-q / 2 - sqrt(Q)) >= 0){
            B = pow((-q / 2 - sqrt(Q)), 1 / 3.0);
        }else{
            B = -pow(-(-q / 2 - sqrt(Q)), 1 / 3.0);
        }
        ReY1=A+B;
        //Im1=0;
        ReY2=-(A+B)/2;
        Im2 = (A - B) / 2.0 * sqrt(3);
        ReY3=ReY2;
        Im3=-Im2;
    }
    //s[4-1]="-(p/3)*(p/3)*(p/3)="+(-(p/3)*(p/3)*(p/3)).ToString()+" ";
    //if(!(Q<0))s[4-1]+=" A="+A.ToString()+" B="+B.ToString();
    //else s[4-1]+=" can't calculate A & B ";
    //if (!(Q > 0)) s[4 - 1] += " Cos(alfa)= " + (CosAlfa).ToString() + " alfa=" + (alfa).ToString();
    //else s[4-1]+=" can't calculate cos & alfa ";
    //s[5 - 1] += " ReY1=" + (ReY1).ToString() + " ReY2=" + (ReY2).ToString() + " ReY3=" + (ReY3).ToString();
    Re1=ReY1-a/3;
    Re2=ReY2-a/3;
    Re3=ReY3-a/3;
    //s[6 - 1] += " Re1=" + (Re1).ToString() + " Re2=" + (Re2).ToString() + " Re3=" + (Re3).ToString();
    ReX.push_back(Re1);
    ImX.push_back(0);
    ReX.push_back(Re2);
    ImX.push_back(Im2);
    ReX.push_back(Re3);
    ImX.push_back(Im3);
}
void EquationCubic(double C2ToC3, double C1ToC3, double C0ToC3, std::vector<double>&ReX, std::vector<double>&ImX){

}

int EquationOf4thPower(std::vector<double>C, std::vector<double>&ReX, std::vector<double>&ImX){
    double c4, c3, c2, c1, c0;
    c0=C[0];
    c1=C[1];
    c2=C[2];
    c3=C[3];
    c4=C[4];
    int ErrN;
    ErrN= EquationOf4thPower(c4, c3, c2, c1, c0, ReX, ImX);
    return ErrN;
}
int EquationOf4thPower(double C[], std::vector<double>&ReX, std::vector<double>&ImX){
    double c4, c3, c2, c1, c0;
    c0=C[0];
    c1=C[1];
    c2=C[2];
    c3=C[3];
    c4=C[4];
    int ErrN;
    ErrN= EquationOf4thPower(c4, c3, c2, c1, c0, ReX, ImX);
    return ErrN;
}
int EquationOf4thPower(double C3To4, double C2To4, double C1To4, double C0To4, std::vector<double>&ReX, std::vector<double>&ImX){
    double c4, c3, c2, c1, c0;
    c0=C0To4;
    c1=C1To4;
    c2=C2To4;
    c3=C3To4;
    c4=1;
    int ErrN;
    ErrN= EquationOf4thPower(c4, c3, c2, c1, c0, ReX, ImX);
    return ErrN;
}
int EquationOf4thPower(double c4, double c3, double c2, double c1, double c0, std::vector<double>&ReX, std::vector<double>&ImX){
            double e2pD;
            bool AllRight=true;
            int ErrN=0;
            std::vector<double>ResolvRe, ResolvIm;
            std::vector<double> EqPow2N1Re, EqPow2N1Im, EqPow2N2Re,EqPow2N2Im;
            //string[] si=new string[27];
            //string []ResolveStr=null;
            double Re1, Im1, Re2, Im2, Re3, Im3, Re4, Im4;
            double ResolvC3, ResolvC2, ResolvC1, ResolvC0;//, CubicReRoot1;
            double Eq1C2,  Eq1C1, Eq1C0, Eq2C2,  Eq2C1, Eq2C0;
            double Root1Eq1Re1, Root1Eq1Im1, Root2Eq1Re2, Root2Eq1Im2,
               Root3Eq2Re1, Root3Eq2Im1, Root4Eq2Re2, Root4Eq2Im2,
               ResolvRe1, ResolvRe2, ResolvIm2, ResolvRe3, ResolvIm3;
            double p,q,r;
            double a=c3/c4, b=c2/c4, c=c1/c4, d=c0/c4;
            //
            e2pD = 0; ResolvRe1 = 0; ResolvRe2 = 0; ResolvIm2 = 0; ResolvRe3 = 0;
            ResolvIm3 = 0; Root3Eq2Re1 = 0; Root3Eq2Im1 = 0; Root4Eq2Re2 = 0; Root4Eq2Im2 = 0;
            Root1Eq1Re1 = 0; Root1Eq1Im1 = 0; Root2Eq1Re2 = 0; Root2Eq1Im2 = 0;
            //Root1Eq2Re1 = 0; Root1Eq2Im1 = 0; Root2Eq2Re2 = 0; Root2Eq2Im2 = 0;
            //
            p=b-3*a*a/8;
            q=a*a*a/8-a*b/2+c;
            r=-3*a*a*a*a/256+a*a*b/16-c*a/4+d;
            ResolvC3=2;
            ResolvC2=-p;
            ResolvC1=-2*r;
            ResolvC0=r*p-q*q/4;
            //si[1-1]="Supplementary Coefficients p="+p.ToString()+" q="+q.ToString()+" r="+r.ToString();
           // si[2-1]="Cubic resolvent: "+ResolvC3.ToString()+"*y^3+"+ResolvC2.ToString()+"*y^2+"+ResolvC1.ToString()+"*y+"+ResolvC0.ToString()+"=0";
            //CubicEquation(ResolvC3, ResolvC2, ResolvC1, ResolvC0, ref ReolvRe1, ref ResolvRe2, ref ResolvIm2, ref ResolvRe3, ref ResolvIm3, ref ResolveStr);
            EquationCubic(ResolvC3, ResolvC2, ResolvC1, ResolvC0, ResolvRe, ResolvIm);
            ResolvRe1=ResolvRe[1-1];
            ResolvRe2=ResolvRe[2-1];
            ResolvIm2=ResolvIm[2-1];
            ResolvRe3=ResolvRe[3-1];
            ResolvIm3=ResolvIm[3-1];
            //si[3-1]="1st Solution of Cubic Resolvent: "+ResolvRe1.ToString();
            //si[4-1]="Checking: "+(ResolvC3*ResolvRe1*ResolvRe1*ResolvRe1+ResolvC2*ResolvRe1*ResolvRe1+ResolvC1*ResolvRe1+ResolvC0).ToString();
            //si[5-1]=ResolveStr[1-1];
            //si[6-1]=ResolveStr[2-1];
            //    si[7-1]=ResolveStr[3-1];
            //si[8-1]=ResolveStr[4-1];
            //si[9-1]=ResolveStr[5-1];
            //    si[10-1]=ResolveStr[6-1];
            //si[11-1]="It was info from cubic resolvent. Now again equation of 4th power";
            if(ResolvRe1<=p/2){
                //MessageBox.Show("Equation can't be replaced by equivalent quadratic equations with real coefs!");
                AllRight=false;
                Root1Eq1Re1=666; Root1Eq1Im1=666; Root2Eq1Re2=666; Root2Eq1Im2=666;
                Root3Eq2Re1=666; Root3Eq2Im1=666; Root4Eq2Re2=666; Root4Eq2Im2=666;
                //si[12-1]="can't calculate the real coefs of equivalent cuadratic equations";
                //si[13-1]="2y1-p=2*"+ResolvRe1.ToString()+"-"+p.ToString()+"="+(2*ResolvRe1-p).ToString();
                //for(int i=14; i<=27; i++){
                //    si[14-1]=" ";
                //}
            }else{
                //si[12-1]="No problem with solution";
                Eq1C2=1;
                Eq1C1=-sqrt(2*ResolvRe1-p);
                Eq1C0 = q / (2 * sqrt(2 * ResolvRe1 - p)) + ResolvRe1;
                //si[13-1]="First quadratic equation: "+Eq1C2.ToString()+"*x^2+"+Eq1C1.ToString()+"*x+"+Eq1C0.ToString()+"=0";
                Eq2C2=Eq1C2;
                Eq2C1=-Eq1C1;
                Eq2C0=-q/(2*sqrt(2*ResolvRe1-p))+ResolvRe1;
                //si[14-1]="Second quadratic equation: "+Eq2C2.ToString()+"*x^2+"+Eq2C1.ToString()+"*x+"+Eq2C0.ToString()+"=0";
                //QuadraticEquation(Eq1C2, Eq1C1, Eq1C0, ref e2pD,  ref Root1Eq1Re1, ref Root1Eq1Im1, ref Root2Eq1Re2, ref Root2Eq1Im2);
                //QuadraticEquation(Eq2C2, Eq2C1, Eq2C0, ref e2pD,  ref Root3Eq2Re1, ref Root3Eq2Im1, ref Root4Eq2Re2, ref Root4Eq2Im2);
                EquationQuadratic(Eq1C2, Eq1C1, Eq1C0, EqPow2N1Re, EqPow2N1Im);
                EquationQuadratic(Eq2C2, Eq2C1, Eq2C0, EqPow2N2Re, EqPow2N2Im);
                Root1Eq1Re1=EqPow2N1Re[1-1];
                Root1Eq1Im1=EqPow2N1Im[1-1];
                Root2Eq1Re2=EqPow2N1Re[2-1];
                Root2Eq1Im2=EqPow2N1Im[2-1];
                //
                Root3Eq2Re1=EqPow2N2Re[1-1];
                Root3Eq2Im1=EqPow2N2Im[1-1];
                Root4Eq2Re2=EqPow2N2Re[2-1];
                Root4Eq2Im2=EqPow2N2Im[2-1];
                //si[15-1]="ANSWER:";
                //si[16-1]="1st Y: ";
                //if(Root1Eq1Im1==0){
                //    si[16-1]+=Root1Eq1Re1.ToString();
                //}else{
                //    si[16-1]+=Root1Eq1Re1.ToString()+"+i*"+Root1Eq1Im1.ToString();
                //}
                //si[17-1]="2nd Y: ";
                //if(Root2Eq1Im2==0){
                //    si[17-1]+=Root2Eq1Re2.ToString();
                //}else{
                //    si[17-1]+=Root2Eq1Re2.ToString()+"+i*"+Root2Eq1Im2.ToString();
                //}
                //si[18-1]="3rd Y: ";
                //if (Root3Eq2Im1 == 0)
                //{
                //    si[18 - 1] += Root3Eq2Re1.ToString();
                //}
                //else
                //{
                //    si[18 - 1] += Root3Eq2Re1.ToString() + "+i*" + Root3Eq2Im1.ToString();
                //}
                //si[19-1]="4th Y: ";
                //if(Root4Eq2Im2==0){
                //    si[19-1]+=Root4Eq2Re2.ToString();
                //}else{
                //    si[19-1]+=Root4Eq2Re2.ToString()+"+i*"+Root4Eq2Im2.ToString();
                //}
                //
                Root1Eq1Re1-=a/4;
                Root2Eq1Re2-=a/4;
                Root3Eq2Re1-=a/4;
                Root4Eq2Re2-=a/4;
                //si[20-1]="1st Root: ";
                //if(Root1Eq1Im1==0){
                //    si[20-1]+=Root1Eq1Re1.ToString();
                //}else{
                //    si[20-1]+=Root1Eq1Re1.ToString()+"+i*"+Root1Eq1Im1.ToString();
                //}
                //si[21-1]="2nd Root: ";
                //if(Root2Eq1Im2==0){
                //    si[21-1]+=Root2Eq1Re2.ToString();
                //}else{
                //    si[21 - 1] += Root2Eq1Re2.ToString() + "+i*" + Root2Eq1Im2.ToString();
                //}
                //si[22-1]="3rd Root: ";
                //if(Root3Eq2Im1==0){
                //    si[22 - 1] += Root3Eq2Re1.ToString();
                //}else{
                //    si[22-1]+=Root3Eq2Re1.ToString()+"+i*"+Root3Eq2Im1.ToString();
                //}
                //si[23-1]="4th Root: ";
                //if(Root4Eq2Im2==0){
                //    si[23-1]+=Root4Eq2Re2.ToString();
                //}else{
                //    si[23 - 1] += Root4Eq2Re2.ToString() + "+i*" + Root4Eq2Im2.ToString();
                //}
                //
                //if(Root1Eq1Im1==0){
                //    si[24-1]="Checking Root 1: "+(Root1Eq1Re1*Root1Eq1Re1*Root1Eq1Re1*Root1Eq1Re1+a*Root1Eq1Re1*Root1Eq1Re1*Root1Eq1Re1+b*Root1Eq1Re1*Root1Eq1Re1+c*Root1Eq1Re1+d).ToString();
                //}else{
                //    si[24-1]="  ";
                //}
                //if(Root2Eq1Im2==0){
                //    si[25 - 1] = "Checking Root 2: " + (Root2Eq1Re2 * Root2Eq1Re2 * Root2Eq1Re2 * Root2Eq1Re2 + a * Root2Eq1Re2 * Root2Eq1Re2 * Root2Eq1Re2 + b * Root2Eq1Re2 * Root2Eq1Re2 + c * Root2Eq1Re2 + d).ToString();
                //}else{
                //    si[25-1]="  ";
                //}
                //if(Root3Eq2Im1==0){
                //    si[26 - 1] = "Checking Root 3: " + (Root3Eq2Re1 * Root3Eq2Re1 * Root3Eq2Re1 * Root3Eq2Re1 + a * Root3Eq2Re1 * Root3Eq2Re1 * Root3Eq2Re1 + b * Root3Eq2Re1 * Root3Eq2Re1 + c * Root3Eq2Re1 + d).ToString();
                //}else{
                //    si[26-1]="  ";
                //}
                //if(Root4Eq2Im2==0){
                //    si[27 - 1] = "Checking Root 4: " + (Root4Eq2Re2 * Root4Eq2Re2 * Root4Eq2Re2 * Root4Eq2Re2 + a * Root4Eq2Re2 * Root4Eq2Re2 * Root4Eq2Re2 + b * Root4Eq2Re2 * Root4Eq2Re2 + c * Root4Eq2Re2 + d).ToString();
                //}else{
                //    si[27-1]="  ";
                //}
            }
            Re1=Root1Eq1Re1; Im1=Root1Eq1Im1; Re2=Root2Eq1Re2; Im2=Root2Eq1Im2;
            Re3=Root3Eq2Re1; Im3=Root3Eq2Im1; Re4=Root4Eq2Re2; Im4=Root4Eq2Im2;
            //s=si;
            ReX.push_back(Re1);
            ReX.push_back(Re2);
            ReX.push_back(Re3);
            ReX.push_back(Re4);
            ImX.push_back(Im1);
            ImX.push_back(Im2);
            ImX.push_back(Im3);
            ImX.push_back(Im4);
            if(AllRight)ErrN=0;
            return ErrN;
        }

//class PolyEq{
//       std::vector<double>C, ReX, ImX;
  //
PolyEq::PolyEq(){}
PolyEq::PolyEq(std::vector<double>C){
    this->SetEquationAndCalcSolution(C);
}
PolyEq::PolyEq(const PolyEq&obj){
    this->C=obj.C;
}
PolyEq::~PolyEq(){}
void PolyEq::Assign(const PolyEq&obj){
    this->C=obj.C;
    this->ReX=obj.ReX;
    this->ImX=obj.ImX;
}
PolyEq PolyEq::operator =(const PolyEq&obj){
    this->Assign(obj);
    return *this;
}
void PolyEq::SetEquationAndCalcSolution(std::vector<double>C){
    int Q=C.size(), order=Q-1;
    this->ReX.clear();
    this->ImX.clear();
    double val;
    int ErrN;
    switch(order){
        case 1:
            if(this->C[1]!=0){
                val=-this->C[1]/this->C[0];
                ReX.push_back(val);
                ImX.push_back(0);
            }
        break;
        case 2:
            EquationQuadratic(this->C, this->ReX, this->ImX);
        break;
        case 3:
            EquationCubic(this->C, this->ReX, this->ImX);
        break;
        case 4:
            ErrN=EquationOf4thPower(this->C, this->ReX, this->ImX);
        break;
    }
}
void PolyEq::StartGenerator(std::vector<double>ReX){
    int Q=ReX.size(), order=Q-0;
    double K=1,
           c11, c10,
           c22, c21, c20,
           c33, c32, c31, c30,
           c44, c43, c42, c41, c40,
           c55, c54, c53, c52, c51, c50;
    this->ReX=ReX;
    this->C.clear();
    this->ImX.clear();
    for(int i=1; i<=order; i++){
       this->ImX.push_back(0);
    }
    c11=K;
    c10=-ReX[1-1];
    //
    c20=ReX[1-1]*ReX[2-1];
    c21=-(ReX[1-1]+ReX[2-1]);
    c22=K;
    //
    c33 = K * 1;
    c32 = -K * (ReX[1 - 1] + ReX[2 - 1] + ReX[3 - 1]);
    c31 = K * (ReX[1 - 1] + ReX[2 - 1]) * +ReX[3 - 1] + ReX[1 - 1] * ReX[2 - 1];
    c30 = -K * (ReX[1 - 1] * ReX[2 - 1] * ReX[3 - 1]);
    //
    c44 = c33;
    c43 = c32 - ReX[4 - 1] * c33;
    c43 = -K * (ReX[1 - 1] + ReX[2 - 1] + ReX[3 - 1] + ReX[4 - 1]);
    c42 = c31 - ReX[4 - 1] * c32;
    c42 = K * ((ReX[1 - 1] + ReX[2 - 1]) * ReX[3 - 1] + ReX[1 - 1] * ReX[2 - 1] + ReX[4 - 1] * (ReX[1 - 1] + ReX[2 - 1] + ReX[3 - 1]));
    c41 = c30 - ReX[4 - 1] * c31;
    c40 = -ReX[4 - 1] * c30;
    c40 = K * ReX[1 - 1] * ReX[2 - 1] * ReX[3 - 1] * ReX[4 - 1];
    //try by analogy
    c55=c44;
    c54=c43-ReX[5-1]*c44;
    c53=c42-ReX[5-1]*c43;
    c52=c41-ReX[5-1]*c42;
    c51=c40-ReX[5-1]*c41;
    c50= K * ReX[1 - 1] * ReX[2 - 1] * ReX[3 - 1] * ReX[4 - 1]* ReX[5 - 1];
    //
    switch(order){
        case 1:
            if(c11!=0){
                this->ReX.push_back(-c10/c11);
                //this->ReX.push_back(c11);
                this->C.push_back(c10);
                this->C.push_back(c11);
            }
        break;
        case 2:
            if(c22>0){
                this->C.push_back(c20);
                this->C.push_back(c21);
                this->C.push_back(c22);
            }
        break;
        case 3:
            if(c33>0){
                this->C.push_back(c30);
                this->C.push_back(c31);
                this->C.push_back(c32);
                this->C.push_back(c33);
            }
        break;
        case 4:
            if(c44>0){
                this->C.push_back(c40);
                this->C.push_back(c41);
                this->C.push_back(c42);
                this->C.push_back(c43);
                this->C.push_back(c44);
            }
        break;
        case 5:

        break;
    }
    this->SetEquationAndCalcSolution(this->C);
}
std::vector<double> PolyEq::GetReX(){return this->ReX;}
std::vector<double> PolyEq::GetImX(){return this->ImX;}
std::vector<double> PolyEq::GetCoefs(){return this->C;}
void PolyEq::SetOrderAndClear(int order){
    this->ReX.clear();
    this->ImX.clear();
    this->C.clear();
    for(int i=1; i<order; i++){
        this->C.push_back(0);
        this->ReX.push_back(0);
        this->ImX.push_back(0);
    }
    this->C.push_back(0);
}
int PolyEq::GetOrder(){
    int order=0;
    order=this->ReX.size();
    return order;
}
QString PolyEq::EquationToString(){
    int powN, order=this->ReX.size(), Q=order+1;
    QString s="", sC, sPow, snum;
    double val, absVal;
    for(int i=0; i<=order; i++){
        powN=order-i;
        val=this->C[powN];
        absVal=fabs(val);
        //snum.setNum(absVal);
        snum.setNum(val);
        sPow.setNum(powN);
        if(val!=0){
            if(powN>1){//if(powN>0){
                if(powN<order){
                    if(val>=0){
                        s+="+";
                    }
                    s+=snum;
                    s+="*X^";
                    s+=sPow;
                }else{
                    s+=snum;
                    s+="*X^";
                    s+=sPow;
                }
        }else if(powN==1){
            if(powN<order){
                if(val>=0){
                    s+="+";
                }
                s+=snum;
                s+="*X";
                //s+="^";
                 //s+=sPow;
            }else{
                s+=snum;
                s+="*X";
                //s+="^";
                //s+=sPow;
            }
        }else{//nur free member
            if(val>=0){
                s+="+";
            }
            s+=snum;
            //s+="*X^";
            //s+=sPow;
            }
        }
    }
    s+="=0";
    return s;
}
std::string PolyEq::EquationToStdString(){
    std::string ss;
    QString qs;
    qs=this->EquationToString();
    ss=qs.toStdString();
    return ss;
}
QString PolyEq::SolutionToString(){
    QString ss, scur, sn;
    int order=this->ReX.size();
    for(int i=1; i<=order-1; i++){
        scur="";
        sn.setNum(i);
        scur="X["+sn+"]=";
        scur+=ComplexNumberToString(this->ReX[i-1], this->ImX[i-1]);
        ss+=scur;
        ss+="; ";
    }
    sn.setNum(order);
    scur="X["+sn+"]=";
    scur+=ComplexNumberToString(this->ReX[order-1], this->ImX[order-1]);
    ss+=scur;
    return ss;
}
std::string PolyEq::SolutionToToStdString(){
    std::string ss;
    QString qs;
    qs=this->SolutionToString();
    ss=qs.toStdString();
    return ss;
}


void LinearApproximation(double*X, double*Y, int Qini, double&a, double&b){
    double meanY=0, sumSquaredDiff=0, stdDevY;
    int Q;
    std::vector<double> newX, newY;
    double sumX=0;
    double sumY=0;
    double sumXY = 0.0;
    double sumX2 = 0.0;
    for(int i=1; i<=Qini; i++){
        meanY+=Y[i-1];
    }
    meanY/=Q;
    for(int i=1; i<=Qini; i++){
        sumSquaredDiff+=(Y[i-1]-meanY)*(Y[i-1]-meanY);
    }
    stdDevY = sqrt(sumSquaredDiff / Qini);
    const double threshold = 2.0 * stdDevY;
    for(int i=1; i<=Qini; i++){
        if(abs(Y[i] - meanY) <= threshold){
            newX.push_back(X[i-1]);
            newY.push_back(Y[i-1]);
        }
    }
    Q=newX.size();
    for(int i=1; i<=Q; i++){
        sumX+=newX[i-1];
        sumY+=newY[i-1];
        sumXY+=newX[i-1]*newY[i-1];
        sumX2+=newX[i-1]*newX[i-1];
    }
    a = (Q * sumXY - sumX * sumY) / (Q * sumX2 - sumX * sumX);
    b = (sumY - a * sumX) / Q;

}

void LinearApproximation(const std::vector<double>&X, const std::vector<double>&Y, double&a, double&b){
    double meanY, sumSquaredDiff;
    std::vector<double>X1=X, Y1=Y;
    int Q=X.size();
    if(X.size()!=Y.size() || X.empty()){
        //NOp
        a=0;
        b=0;
    }else{
        //void LinearApproximation(double*&X, double*Y, int Qini, double&a, double&b)
        LinearApproximation(X1.data(), Y1.data(), X.size(), a, b);
    }
}


double arctg2(double x, double y){
    double R=0, absval;
    int QuadrantN=0;
    if(x=0 && y!=0){

    }else if(y==0 &&x!=0){
        R=0;
    }else{
        if(x>0 && y>0){
            QuadrantN=1;
        }else if(x<0 && y>0){
            QuadrantN=2;
        }else if(x<0 && y<0){
            QuadrantN=3;
        }else if(x>0 && y<0){
            QuadrantN=4;
        }
        absval=fabs(y/x);
        switch(QuadrantN){
            case 1:
                R=absval;
            break;
            case 2:
                R=M_PI-absval;
            break;
            case 3:
                R=M_PI+absval;
            break;
            case 4:
                R=2*M_PI-absval;
            break;
        }
    }
    return R;
}

//-------------------------------------------------------------------------------------------------------------


//class ComplexNum{
//    double re, im;
//public:
    ComplexNum::ComplexNum(double re, double im){
        this->SetDecart(re, im);
    }

    ComplexNum::ComplexNum(const ComplexNum&obj){
        this->re=re;
        this->im=im;
    }

    void ComplexNum::Assign(const ComplexNum&obj){
        this->SetDecart(re, im);
    }

    ComplexNum::~ComplexNum(){}
    ComplexNum& ComplexNum::operator = (const ComplexNum&obj){
        this->Assign(obj);
        return(*(this));
    }

    double ComplexNum::SetDecart(double re, double im){
        this->re=re;
        this->im=im;
    }

    double ComplexNum::SetPolar(double fi, double r){
        this->re=r*cos(fi);
        this->im=r*sin(fi);
    }

    double ComplexNum::GetRe(){
        return this->re;
    }

    double ComplexNum::GetIm(){
        return this->im;
    }
    void ComplexNum::GetPolar(double&r, double&fi){
        r=sqrt(this->re*this->re + this->im*this->im);
        fi=arctg2(this->re, this->im);
    }

    double ComplexNum::GetFi(){
        double r=0, fi=0;
        this->GetPolar(r, fi);
        return fi;
    }

    double ComplexNum::GetR(){
        double r=0, fi=0;
        this->GetPolar(r, fi);
        return r;
    }
    void ComplexNum::ChangeRe(double val){
        this->re=val;
    }

    void ComplexNum::ChangeIm(double val){
        this->im=val;
    }
    void ComplexNum::ChangeFi(double val){
        double r=GetR();
        this->SetPolar(val, r);
    }
    ComplexNum ComplexNum::GetConjugated(){
        ComplexNum y(this->re, -this->im);
        return y;
    }

    void ComplexNum::ChangeR(double val){
        double fi=GetFi();
        this->SetPolar(fi, val);
    }
    ComplexNum& ComplexNum::operator + (const ComplexNum&obj){
        this->re+=obj.re;
        this->im+=obj.im;
        return*this;
    }

    ComplexNum& ComplexNum::operator - (const ComplexNum&obj){
        this->re-=obj.re;
        this->im-=obj.im;
        return*this;
    }
    ComplexNum& ComplexNum::operator * (const ComplexNum&obj){
        this->re = this->re * obj.re + this->im * obj.im;
        this->im = this->re * obj.im - this->im * obj.re;
        return*this;
    }
    ComplexNum& ComplexNum::operator / (const ComplexNum&obj){
        double denominator=this->re*this->re+this->im*this->im;
        this->re = this->re * obj.re - this->im * obj.im;
        this->im = this->re * obj.im + this->im * obj.re;
        this->re/=denominator;
        this->im/=denominator;
        return*this;
    }
    ComplexNum& ComplexNum::operator * (double val){
        this->re*=val;
        this->im*=val;
        return (*(this));
    }

    bool ComplexNum::operator == (const ComplexNum&obj){
        return (this->re == obj.re && this->im == obj.im);
    }

    bool ComplexNum::operator != (const ComplexNum&obj){
        return !(this->re == obj.re && this->im == obj.im);
    }

    bool ComplexNum::IsConjugated(const ComplexNum&obj){
        return(this->re==obj.re && this->im==-obj.im);
    }

//};



//============================================================================================================================

    void Matrix_Prototype::SetOne(double val){
        this->SetNull();
        this->SetSize(1, 1);
        this->SetComponent(val, 1, 1);
    }
    //
    bool Matrix_Prototype::LineNBelongsHere(int LineN) const{
        int QLines=this->GetQLines(), QColumns=this->GetQColumns();
        bool b=(LineN >=1 && LineN <=QLines);
        return b;
    }
    bool Matrix_Prototype::ColNBelongsHere(int ColN) const{
        int QLines=this->GetQLines(), QColumns=this->GetQColumns();
        bool b=(ColN>=1 && ColN<=QColumns);
        return b;
    }
    bool Matrix_Prototype::NsBelongHere(int LineN, int ColN) const{
        int QLines=this->GetQLines(), QColumns=this->GetQColumns();
        bool b=(LineN >=1 && LineN <=QLines && ColN>=1 && ColN<=QColumns);
        return b;
    }



    std::vector<double> Matrix_Prototype::GetArrayOfLineN(int LineN) const{
        std::vector<double> row;
        double val;
        int QLines=this->GetQLines(), QColumns=this->GetQColumns();
        if(this->LineNBelongsHere(LineN)){
            for(int i=1; i<=QColumns; i++){
                val=this->GetComponent(LineN, i);
                row.push_back(val);
            }
        }
        return row;
    }

    std::vector<double> Matrix_Prototype::GetArrayOfColN(int ColN) const{
        std::vector<double> row;
        double val;
        int QLines=this->GetQLines(), QColumns=this->GetQColumns();
        if(this->LineNBelongsHere(ColN)){
            for(int i=1; i<=QLines; i++){
                val=this->GetComponent(i, ColN);
                row.push_back(val);
            }
        }
        return row;
    }

    std::vector<double> Matrix_Prototype::GetContentAsArray1DByLines() const{
        double val;
        int QLins=GetQLines(), QCols=GetQColumns();
        std::vector<double>data;
        for(int i=1; i<=QLins; i++){
            for(int j=1; j<=QCols; j++){
                val=this->GetComponent(i, j);
                data.push_back(val);
            }
        }
        return data;
    }
    std::vector<std::vector<double>> Matrix_Prototype::GetContentAsArray2DByLines() const{
        int QLins=GetQLines(), QCols=GetQColumns();
        std::vector<double>row;
        std::vector<std::vector<double>>R;
        double val;
        if(QCols>-0 && QLins>0){
            for(int i=1; i<=QLins; i++){
                row.clear();
                for(int j=1; j<=QLins; j++){
                    val=this->GetComponent(i, j);
                    row.push_back(val);
                }
                R.push_back(row);
            }
        }
        return R;
    }


    void Matrix_Prototype::Set(double**y, int QLins, int QCols){
        //this->SetNull();
        this->SetSize(QLins, QCols);
        if(y!=NULL){
            for(int i=1; i<=QLins; i++){
                for(int j=1; j<=QCols; j++){
                    this->SetComponent(y[i-1][j-1], i, j);
                }
            }
        }else{
            for(int i=1; i<=QLins; i++){
                for(int j=1; j<=QCols; j++){
                    this->SetComponent(0, i, j);
                }
            }
        }
        //size
    }
    void Matrix_Prototype::Set(std::vector<std::vector<double>>y){
        //this->SetNull();
        int QLins=y.size(), QCols=y[1-1].size();
        this->SetSize(QLins, QCols);
        for(int i=1; i<=QLins; i++){
            for(int j=1; j<=QCols; j++){
                //this->y[i-1][j-1]=y[i-1][j-1];
                this->SetComponent(y[i-1][j-1], i, j);
            }
        }
        //size
    }
    void Matrix_Prototype::Set(double*y, int QLins, int QCols){
        int N=0;
        //if(ColOfVectorNotHorLine){
            this->SetSize(QLins, QCols);
            for(int i=1; i<=QLins; i++){
                for(int j=1; j<=QCols; j++){
                    N++;
                    this->SetComponent(y[N-1], i, j);
                }
            }
            //this->QLins=Qlins;
            //this->QCols=QCols;
            //size - this is for concrete type!
        //}
    }
    void  Matrix_Prototype::Set(std::vector<double>y, int QCols){
        this->Set(y.data(), y.size()/QCols, QCols);
        //size
    }
    void Matrix_Prototype::Set(double x, double y, double z){
         this->SetNull();
         this->SetSize(3, 1);
          this->SetComponent(x, 1, 1);
          this->SetComponent(y, 1, 2);
          this->SetComponent(z, 1, 3);
         //size
    }
    void Matrix_Prototype::Set(double y, int QLins, int QCols){
        std::vector<double>z;
        this->SetSize(QLins, QCols);
        for(int i=1; i<=QLins; i++){
            for(int j=1; j<=QCols; j++){
                this->SetComponent(y, i, j);
            }
        }
    }

    void Matrix_Prototype::SetX(double val){ this->SetComponent(val, 1, 1); }
    void Matrix_Prototype::SetY(double val){ this->SetComponent(val, 2, 1); }
    void Matrix_Prototype::SetZ(double val){ this->SetComponent(val, 3, 1); }
    double Matrix_Prototype::GetX(){ return this->GetComponent(1, 1); }
    double Matrix_Prototype::GetY(){ return this->GetComponent(2, 1); }
    double Matrix_Prototype::GetZ(){ return this->GetComponent(3, 1); }
    //
    void Matrix_Prototype::SetLine(int N, double*x, int Q, int whatN1, int QDfltValsBefore, bool DefaulsAreZerosNotOwnVals){
        int LreqCalcd=0, preN1_=0, preN2_=0, ownN1_=0, ownN2=0, ForOwnN1_=0, ForOwnN2=0, postN1=0, postN2_=0, Lold=Q;//x.size();
        int QLines=this->GetQLines(), QColumns=this->GetQColumns();
        int LreqGiven=QColumns;
        double val;
        std::vector<double>R;
        if(N>=1 && N<=QLines){
            Calc_Array1DSubArray_byLs_MarkupNs_vFmls(LreqCalcd, preN1_, preN2_, ownN1_, ownN2, ForOwnN1_, ForOwnN2, postN1, postN2_, Lold, LreqGiven, whatN1, QDfltValsBefore);
            if(preN1_>0){
                if(DefaulsAreZerosNotOwnVals){
                    for(int i=preN1_; i<=preN2_; i++){
                        R.push_back(0);
                    }
                }else{
                    for(int i=preN1_; i<=preN2_; i++){
                        val=this->GetComponent(N, i);
                        R.push_back(val);
                    }
                }
            }
            if(ownN1_>0){
                for(int i=ownN1_; i<=ownN2; i++){
                    val=x[i-1];
                    R.push_back(val);
                }
            }
            if(postN1>0){
                if(DefaulsAreZerosNotOwnVals){
                    for(int i=postN1; i<=postN2_; i++){
                        R.push_back(0);
                    }
                }else{
                    for(int i=postN1; i<=postN2_; i++){
                        val=this->GetComponent(N, i);
                        R.push_back(val);
                    }
                }
            }
            for(int i=1; i<=QColumns; i++){
                val=R[i-1];
                this->SetComponent(val, N, i);
            }
        }
    }
    void Matrix_Prototype::SetColumn(int N, double*x, int Q, int whatN1, int QDfltValsBefore, bool DefaulsAreZerosNotOwnVals){
        int LreqCalcd=0, preN1_=0, preN2_=0, ownN1_=0, ownN2=0, ForOwnN1_=0, ForOwnN2=0, postN1=0, postN2_=0, Lold=Q;//=x.size();
        int QLines=this->GetQLines(), QColumns=this->GetQColumns();
        int LreqGiven=QLines;
        double val;
        std::vector<double>R;
        if(N>=1 && N<=QColumns){
            Calc_Array1DSubArray_byLs_MarkupNs_vFmls(LreqCalcd, preN1_, preN2_, ownN1_, ownN2, ForOwnN1_, ForOwnN2, postN1, postN2_, Lold, LreqGiven, whatN1, QDfltValsBefore);
            if(preN1_>0){
                if(DefaulsAreZerosNotOwnVals){
                    for(int i=preN1_; i<=preN2_; i++){
                        R.push_back(0);
                    }
                }else{
                    for(int i=preN1_; i<=preN2_; i++){
                        val=this->GetComponent(N, i);
                        R.push_back(val);
                    }
                }
            }
            if(ownN1_>0){
                for(int i=ownN1_; i<=ownN2; i++){
                    val=x[i-1];
                    R.push_back(val);
                }
            }
            if(postN1>0){
                if(DefaulsAreZerosNotOwnVals){
                    for(int i=postN1; i<=postN2_; i++){
                        R.push_back(0);
                    }
                }else{
                    for(int i=postN1; i<=postN2_; i++){
                        val=this->GetComponent(N, i);
                        R.push_back(val);
                    }
                }
            }
            for(int i=1; i<=QLines; i++){
                val=R[i-1];
                this->SetComponent(val, i, N);
            }
        }
    }
    void  Matrix_Prototype::SetLine(int N, std::vector<double>x, int whatN1, int QDfltValsBefore, bool DefaulsAreZerosNotOwnVals){
        this->SetLine(N, x.data(), x.size(), whatN1, QDfltValsBefore, DefaulsAreZerosNotOwnVals);
    }
    void  Matrix_Prototype::SetColumn(int N, std::vector<double>x, int whatN1, int QDfltValsBefore, bool DefaulsAreZerosNotOwnVals){
        this->SetColumn(N, x.data(), x.size(), whatN1, QDfltValsBefore, DefaulsAreZerosNotOwnVals);
    }
    //


    //
    void Matrix_Prototype::ConcatAdding(Matrix_Prototype*obj, int QShiftedCols, int FromN, int QAdded){
        std::vector<double>row;
        double val;
        int QLinesThisWas=this->GetQLines(),  QLinesObj=obj->GetQLines(),  QLinesRslt,
            QColsThisWas=this->GetQColumns(), QColsObj=obj->GetQColumns(), QColsRslt;
        if(QAdded==0){
            QAdded=QLinesObj-FromN+1;
        }
        if(QShiftedCols>0){
            QColsRslt = QColsThisWas >= QShiftedCols+QColsObj ? QColsThisWas : QShiftedCols+QColsObj;
        }else{
            QColsRslt = -QShiftedCols+QColsObj >= QColsThisWas ? -QShiftedCols+QColsObj : QColsThisWas;
        }
        QLinesRslt=QLinesThisWas+QAdded;
        if(QShiftedCols<0){
            for(int i=1; i<=-QShiftedCols; i++){
                this->InsColumn(1, NULL, this->GetQLines());//(1, 0);
            }
        }
        if(QShiftedCols>=0){
            for(int i=FromN; i<=FromN-1+QAdded; i++){
                for(int j=1; j<=QColsObj; j++){
                    val=obj->GetComponent(i, j);
                    this->SetComponent(val, QLinesThisWas+i, QShiftedCols+j);
                }
            }
        }else{
            for(int i=FromN; i<=FromN-1+QAdded; i++){
                for(int j=1; j<=QColsObj; j++){
                    val=obj->GetComponent(i, j);
                    this->SetComponent(val, QLinesThisWas+i, j);
                }
            }
        }
    }

    void Matrix_Prototype::ConcatStretching(Matrix_Prototype*obj, int QShiftedLines, int FromN, int QAdded){
        std::vector<double>row;
        double val;
        int QLinesThisWas=this->GetQLines(),  QLinesObj=obj->GetQLines(),  QLinesRslt,
            QColsThisWas=this->GetQColumns(), QColsObj=obj->GetQColumns(), QColsRslt;
        if(QAdded==0){
            QAdded=QColsObj-FromN+1;
        }
        if(QShiftedLines>0){
            QLinesRslt = QLinesThisWas >= QShiftedLines+QLinesObj ? QLinesThisWas : QShiftedLines+QLinesObj;
        }else{
            QLinesRslt = -QShiftedLines+QLinesThisWas >= QLinesObj ? -QShiftedLines+QLinesThisWas : QLinesObj;
        }
        QColsRslt=QColsThisWas+QAdded;
        if(QShiftedLines<0){
            for(int i=1; i<=-QShiftedLines; i++){
                this->InsLine(1, NULL, this->GetQColumns());//(1, 0);
            }
        }
        if(QShiftedLines>=0){
            for(int i=1; i<=QLinesObj; i++){
                for(int j=FromN; j<=FromN-1+QAdded; j++){
                    val=obj->GetComponent(i, j);
                    this->SetComponent(val, QShiftedLines+i, QColsThisWas+j);
                }
            }
        }else{
            for(int i=1; i<=QLinesObj; i++){
                for(int j=FromN; j<=FromN-1+QAdded; j++){
                    val=obj->GetComponent(i, j);
                    this->SetComponent(val, i, QColsThisWas+j);                }
            }
        }
    }




    QString Matrix_Prototype::LineToString(int LineN, QString delimElem, bool useBrackets){
        double val;
        QString s="", sc;
        int QLines=this->GetQLines(), QColumns=this->GetQColumns();
        if(useBrackets)s+="[";
        if(LineN>=1 && LineN<=QLines && QColumns>0){
            for(int i=1; i<=QColumns-1; i++){
                val=this->GetComponent(LineN, i);
                sc.setNum(val);
                s+=sc;
                s+=delimElem;
            }
            val=this->GetComponent(LineN, QColumns);
            sc.setNum(val);
            s+=sc;
            if(useBrackets)s+="]";
        }
        return s;
    }
    QString Matrix_Prototype::ColumnToString(int ColN, QString delimElem, bool useBrackets){
        double val;
        QString s="", sc;
        int QLines=this->GetQLines(), QColumns=this->GetQColumns();
        if(useBrackets)s+="[";
        if(ColN>=1 && ColN<=QColumns && QLines>0){
            for(int i=1; i<=QLines-1; i++){
                val=this->GetComponent(i, ColN);
                sc.setNum(val);
                s+=sc;
                s+=delimElem;
            }
            val=this->GetComponent(QLines, ColN);
            sc.setNum(val);
            s+=sc;
            if(useBrackets)s+="]";
        }
        return s;
    }
    QString Matrix_Prototype::ToString(QString delimElem, QString delimLines, bool useBrackets){
        int QLines=this->GetQLines(), QColumns=this->GetQColumns();
        QString s="", sc;
        if(useBrackets)s+="[";
        if(QLines>0 && QColumns>0){
            for(int i=1; i<=QColumns-1; i++){
                //sc=this->GetLineAsStr(i, delimElem, useBrackets);
                sc=this->LineToString(i, delimElem, useBrackets);
                s+=sc;
                s+=delimLines;
            }
            //sc=this->GetLineAsStr(QLines, delimElem, useBrackets);
            sc=this->LineToString(QLines, delimElem, useBrackets);
            s+=sc;
            if(useBrackets)s+="]";
        }
        return s;
    }
    void Matrix_Prototype::StdCOut1D(QString delimElem, QString delimRow, bool useBrackets){
        //QString qs=this->GetAllLinesAsStr(delimElem, delimRow, useBrackets);
        QString qs=this->ToString(delimElem, delimRow, useBrackets);
        std::string ss=qs.toStdString();
        std::cout<<ss.c_str();
    }

    void Matrix_Prototype::StdCOut2D(QString delimElem, QString delimRow, bool useBrackets){
        int QLines=this->GetQLines(), QColumns=this->GetQColumns();
        QString qs;
        std::string ss;
        if(useBrackets)std::cout<<"[";
        if(this->GetQLines()>0 && this->GetQColumns()>0){
            //qs=this->GetLineAsStr(1, delimElem, useBrackets);
            qs=this->LineToString(1, delimElem, useBrackets);
            ss=qs.toStdString();
            std::cout<<ss.c_str()<<std::endl;
            for(int i=2; i<=QLines-1; i++){
                if(useBrackets)std::cout<<" ";
                //qs=this->GetLineAsStr(i, delimElem, useBrackets);
                qs=this->LineToString(i, delimElem, useBrackets);
                ss=qs.toStdString();
                std::cout<<ss.c_str()<<std::endl;
            }
            if(useBrackets)std::cout<<" ";
            //s=this->GetLineAsStr(QLines, delimElem, useBrackets);
            qs=this->LineToString(QLines, delimElem, useBrackets);
            ss=qs.toStdString();
            std::cout<<ss.c_str()<<std::endl;
        }
        if(useBrackets)std::cout<<"]";
    }//eofn
    
//};

   // class Matrix_P2S{
   //     double**y;
   //     int QLins, QCols;
   //     public:
Matrix_P2S::Matrix_P2S(int QLines, int QColumns, double val){
    std::cout<<"Matrix_P2S constr (dbl, QL, QC) works"<<std::endl;
    this->construct();
    //this->Set(x, QLines, QColumns);
    this->SetSize(QLines, QColumns);
    for(int i=1; i<=QLines; i++){
        for(int j=1; j<=QColumns; j++){
            this->SetComponent(val, i, j);
        }
    }
}
//Matrix_P2S::Matrix_P2S(int QLins, int QCols){
//    std::cout<<"Matrix_P2S constr (QL, QC) works"<<std::endl;
//    this->construct();
//    this->SetSize(QLins, QCols);
//}
Matrix_P2S::Matrix_P2S(double**y, int QLins, int QCols){
    std::cout<<"Matrix_P2S constr (**dbl, QL, QC) works"<<std::endl;
    this->construct();
    this->Set(y, QLins, QCols);
}
Matrix_P2S::Matrix_P2S(std::vector<std::vector<double>>y){
    std::cout<<"Matrix_P2S constr (vector<<dbl>>, QL, QC) works"<<std::endl;
    this->construct();
    this->Set(y);
}
Matrix_P2S::Matrix_P2S(double*y, int QLines, int QColumns){
    std::cout<<"Matrix_P2S constr (*dbl, QL, QC) works"<<std::endl;
    this->construct();
    this->Set(y, QLines, QColumns);
}
Matrix_P2S::Matrix_P2S(std::vector<double>y, int QColumns){
    std::cout<<"Matrix_P2S constr (vector<dbl>,  bool) works"<<std::endl;
    this->construct();
    this->Set(y, QColumns);
}
Matrix_P2S::Matrix_P2S(double x, double y, double z){
    std::cout<<"Matrix_P2S constr (x, y, z) works"<<std::endl;
    this->construct();
    this->Set(x, y, z);
}
//Matrix_P2S::Matrix_P2S(double x, int QLines, int QColumns){
//    std::cout<<"Matrix_P2S constr (dbl, QL, QC) works"<<std::endl;
//    this->construct();
//    this->Set(x, QLines, QColumns);
//}
//
Matrix_P2S::~Matrix_P2S(){
    std::cout<<"Matrix_P2S Destructor works"<<std::endl;
    this->SetNull();
}
//
void Matrix_P2S::SetNull(){
    if(this->y!=NULL){
        for(int i=1; i<=this->QLins; i++){
            delete[] this->y[i-1];
        }
        delete[] this->y;
    }
    this->y=NULL;
    this->QLins=0;
    this->QCols=0;
}
void Matrix_P2S::construct(){
    this->y=NULL;
    this->QLins=0;
    this->QCols=0;
    //
    SetOne(0);
}
Matrix_Prototype* Matrix_P2S::CreateMatrix(){
    return new Matrix_P2S();
}
Matrix_Prototype* Matrix_P2S::clone(){
    return this;
}
int Matrix_P2S::GetTypeN() const{
    return 121;//P/V, 1/2, S/L/A
}

//
void Matrix_P2S::SetSize(int QLins, int QCols){
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
//
int Matrix_P2S::GetQLines() const { return this->QLins; }
int Matrix_P2S::GetQColumns() const { return this->QCols; }


void Matrix_P2S::Set(double**y, int QLins, int QCols){
    Matrix_Prototype::Set(y, QLins, QCols);
    this->QLins=QLins;
    this->QCols=QCols;
}
void Matrix_P2S::Set(std::vector<std::vector<double>>y){
    this->QLins=y.size(), this->QCols=y[1-1].size();
    Matrix_Prototype::Set(y);
}
void Matrix_P2S::Set(double*y, int QLines,  int QColumns){
    Matrix_Prototype::Set(y, QLins, QCols);
    this->QLins=QLins;
    this->QCols=QCols;
}
void Matrix_P2S::Set(std::vector<double>y, int QColumns){
    Matrix_Prototype::Set(y, QColumns);
    this->QCols=QColumns;
    this->QLins=y.size()/QColumns;
}
void Matrix_P2S::Set(double x, double y, double z){
    Matrix_Prototype::Set(x, y, z);
    this->QLins=3;
    this->QCols=1;
}
void Matrix_P2S::Set(double y, int QLines, int QColumns){
    Matrix_Prototype::Set(y, QLines, QColumns);
    this->QCols=QColumns;
    this->QLins=QLines;
}


double Matrix_P2S::GetComponent(int LineN, int ColN) const{
    double y=0;
    if(LineN>=1 && LineN<=this->QLins && ColN>=1 && ColN<=this->QCols){
        y=this->y[LineN-1][ColN-1];
    }
    return y;
}
void Matrix_P2S::SetComponent(double val, int LineN, int ColN){
    if(LineN>=1 && LineN<=this->QLins && ColN>=1 && ColN<=this->QCols){
        this->y[LineN-1][ColN-1]=val;
    }
}
//
void Matrix_P2S::SetLine(int N, double*x, int L, int FromN, int QDefaultBefore, bool DefaultValsAreZerosNotOwn){
    //std::vector<double>row;
    int QL=this->GetQLines(), QC=this->GetQColumns();
    Array2DSize size(QL, QC);
    //int preN1_, preN2_, ownN1_, ownN2, ForOwnN1_, ForOwnN2, postN1, postN2_;
    //int LreqCalcd;
    //void Calc_Array1DSubArray_byLs_MarkupNs_vFmls(int&LreqCalcd, int&preN1_, int&preN2_, int&ownN1_, int&ownN2, int&ForOwnN1_, int&ForOwnN2, int&postN1, int&postN2_, int Lini, int LreqGiven, int whatFromN, int QDefaultValuesBefore){
    //Calc_Array1DSubArray_byLs_MarkupNs_vFmls(LreqCalcd, preN1_, preN2_, ownN1_, ownN2, ForOwnN1_, ForOwnN2, postN1, postN2_, QC, 0, whatFromN, QDefaultValuesBefore);
    double *dflt=new double;
    (*(dflt))=0;
    //template<typename T> void Array2DSetExtRowN(T**&X, Array2DSize&size, int ExtRowN, int wholeRowL=0, T*rowParam=NULL, int whatFromN=1, int QDefaultValuesBefore=0, bool KeepIfRect=true, T*DefaultValParam=NULL, bool DefaultIsSpecNotOwn=false){
    Array2DSetExtRowN(                          this->y,             size,           N,               L,               x,           FromN,             QDefaultBefore,                false,                   dflt,       DefaultValsAreZerosNotOwn);
    delete dflt;
}
void  Matrix_P2S::SetLine(int N, std::vector<double>x, int FromN, int QDfltValsBefore, bool DefaultValsAreOwnNotZeros){
    //template<typename T> void Array2DSetExtRowN(std::vector<std::vector<T>>X, int ExtRowN, std::vector<T>rowParam, T*DefaultValParam=NULL, int FromN=1, int DefaultValuesBefore=0, bool KeepIfRect=false)
    //template<typename T> void Array2DSetExtRowN(T**&X, const Array2DSize&size, int ExtRowN, std::vector<T>rowParam, T*DefaultValParam=NULL, int FromN=1, int DefaultValuesBefore=0, bool KeepIfRect=false)
    //template<typename T> void Array2DSetExtRowN(T**&X, const Array2DSize&size, int ExtRowN, int wholeRowL=0, T*rowParam=NULL, T*DefaultValParam=NULL, int FromN=1, int DefaultValuesBefore=0, bool KeepIfRect=false)
    Array2DSize size(this->GetQLines(), this->GetQColumns());
    double*defaultVal=new double;
    *defaultVal=0;
    //Array2DSetExtRowN(T**&X,  const Array2DSize&size, int ExtRowN, std::vector<T>rowParam, T*DefaultValParam=NULL, int FromN=1, int DefaultValuesBefore=0, bool KeepIfRect=false)
    Array2DSetExtRowN(this->y,                    size,           N,                      x,             defaultVal,      FromN,           QDfltValsBefore, true);
    delete defaultVal;
}
void  Matrix_P2S::SetColumn(int N, double*x, int L, int FromN, int QDefaultBefore, bool DefaultValsAreOwnNotZeros){
    int QL=this->GetQLines(), QC=this->GetQColumns();
    Array2DSize size(QL,  QC);
    double*dflt=new double;
    (*(dflt))=0;
    //template<typename T>Array2DSetIneRow(T**&X, Array2DSize& size, int IneRowN, T*rowParam=NULL, int wholeRowL=0, int whatFromN=1, int QDefaultValsBefore=0, T*DefaultValParam=NULL, bool ExistingInsteadOfDefault=true, bool ExistingOnlyForGTLminNotIgnore=false, bool LastPossIfNotRectAtPos0NotIgnore=false){
    Array2DSetIneRow(                    this->y,              size,           N,               x,               L,           FromN,       QDefaultBefore,             dflt,                    DefaultValsAreOwnNotZeros);//            ExistingOnlyForGTLminNotIgnore,      LastPossIfNotRectAtPos0NotIgnore);
    delete dflt;
}
void  Matrix_P2S::SetColumn(int N, std::vector<double>x, int whatN1, int QDfltValsBefore, bool DefaultValsAreOwnNotZeros){
    //template<typename T> void Array2DSetIneRowN(T**&X, Array2DSize&size, int IneRowN, int wholeRowL=0, T*rowParam=NULL, T*DefaultValParam=NULL, int FromN=1, int DefaultValuesBefore=0, int IfNGTLmin_Ignore0DoNil1Add2Srtetch3=1)
    //template<typename T> void Array2DSetIneRowN(T**&X, Array2DSize&size, int IneRowN, std::vector<T>rowParam, T*DefaultValExt=NULL, int FromN=1, int DefaultValuesBefore=0, int IfNGTLmin_Ignore0DoNil1Add2Srtetch3=1)
    //template<typename T> void Array2DSetIneRowN(std::vector<std::vector<T>>&X, std::vector<T>rowParam, int IneRowN, T*DefaultValExt=NULL, int FromN=1, int DefaultValuesBefore=0, int IfNGTLmin_Ignore0DoNil1Add2Srtetch3=1)
    Array2DSize size(this->GetQLines(), this->GetQColumns());
    double*defaultVal=new double;
    *defaultVal=0;
    //template<typename T>Array2DSetIneRow(T**&X, Array2DSize& size, int IneRowN, std::vector<T>rowParam, int whatFromN=1, int QDefaultValsBefore=0, T*DefaultValParam=NULL, bool ExistingInsteadOfDefault=true, bool ExistingOnlyForGTLminNotIgnore=false, bool LastPossIfNotRectAtPos0NotIgnore=false){
    Array2DSetIneRow(                    this->y,              size,           N,                      x, whatN1,                   QDfltValsBefore, defaultVal);
    delete defaultVal;
}
//
void Matrix_P2S::AddLine(double*x, int L, int whatN1, int QDefaultValsBefore){
    //template<typename T>void Array2DAddExtRow(T**&X, Array2DSize&size, T*rowParam, int wholeRowL=0, int whereFromN=1, int whatFromN=1, bool RectNotVar=false, T*DfltValParam=NULL, TValsShowHide*vsh=NULL)
    int QLines=this->GetQLines(), QColumns=this->GetQColumns();
    Array2DSize size(QLines, QColumns);
    double*defaultVal=new double;
    *defaultVal=0;
    if((QLines==1 && QColumns==1 && L>1) || ((QLines==0 || QColumns==0)&&L>0)){
        this->SetNull();
    }
    Array2DAddExtRow(this->y, size, x, L, 1, whatN1, true, defaultVal);
    delete defaultVal;
}
void Matrix_P2S::AddLine(std::vector<double>x, int whatN1, int QDefaultValsBefore){
   int QLines=this->GetQLines(), QColumns=this->GetQColumns(), L=x.size();
    Array2DSize size(QLines, QColumns);
    //template<typename T>void Array2DAddExtRow(T**&X, Array2DSize&size, T*rowParam, int wholeRowL=0, int whereFromN=1, int whatFromN=1, bool RectNotVar=false, T*DfltValParam=NULL, TValsShowHide*vsh=NULL){//44 //better but not checked
    double*defaultVal=new double;
    *defaultVal=0;
    if((QLines==1 && QColumns==1 && L>1) || ((QLines==0 || QColumns==0)&&L>0)){
        this->SetNull();
    }
    Array2DAddExtRow(this->y, size, x.data(), x.size(), 1, whatN1, true, defaultVal);
    delete defaultVal;
}
void Matrix_P2S::AddColumn(double*x, int L, int whatN1, int QDefaultValsBefore){
    //template<typename T>void Array2DAddIneRow(T**&X, Array2DSize&size, T*rowParam, int wholeRowL=0, T*DfltValParam=NULL, bool AddToEachNotStretchToMax=true, int whatFromN=1, int whereFromN=1
    int QLines=this->GetQLines(), QColumns=this->GetQColumns();
    Array2DSize size(QLines, QColumns);
    double*defaultVal=new double;
    *defaultVal=0;
    if((QLines==1 && QColumns==1 && L>1) || ((QLines==0 || QColumns==0)&&L>0)){
        this->SetNull();
    }
    Array2DAddExtRow(this->y, size, x, L, 1, whatN1, true, defaultVal);
    delete defaultVal;
}
void Matrix_P2S::AddColumn(std::vector<double>x, int FromN, int QDefaultValsBefore){
    //template<typename T>void Array2DAddIneRow(T**&X, Array2DSize&size, T*rowParam, int wholeRowL=0, T*DfltValParam=NULL, bool AddToEachNotStretchToMax=true, int whatFromN=1, int whereFromN=1){
    int QLines=this->GetQLines(), QColumns=this->GetQColumns(), L=x.size();
    Array2DSize size(QLines, QColumns);
    double*defaultVal=new double;
    *defaultVal=0;
    if( (QLines==1 && QColumns==1 && L>1) || ((QLines==0||QColumns==0)&&L>0) ){
        this->SetNull();
    }
    //template<typename T>Array2DAddIneRow(T**&X, Array2DSize& size, std::vector<T>rowParam, int IfNonRectIgnore0Add1Stretch2=0, int whatFromN=1, int QDefaultValsBefore=0, T*DefaultValParam=NULL){
    Array2DAddIneRow(                    this->y,              size,                      x,                                  0,           FromN,                        0,  defaultVal);
    delete defaultVal;
}

void Matrix_P2S::InsLine(int N, double*x, int L, int FromN, int QDefaultBefore){
    Array2DSize size(this->GetQLines(), this->GetQColumns());
    double*defaultVal=new double;
    *defaultVal=0;
    //template<typename T> void Array2DInsExtRow(T**&X,  Array2DSize&size, int ExtRowN, int wholeRowL=0, T*rowParam=NULL, T*DefaultValParam=NULL, int whatFromN=1, int QDefaultValuesBefore=0, bool KeepIfRect=true, TValsShowHide*vsh=NULL){
    Array2DInsExtRow(                          this->y,              size,           N,             L,                 x,             defaultVal,           FromN,             QDefaultBefore,                 true);
    delete defaultVal;
}
void Matrix_P2S::InsLine(int N, std::vector<double>x, int whatFromN, int QDefaultValsBefore){
     Array2DSize size(this->GetQLines(), this->GetQColumns());
     double*defaultVal=new double;
     *defaultVal=0;
     //template<typename T> void Array2DInsExtRow(T**&X,  Array2DSize&size, int ExtRowN, std::vector<T>rowParam, T*DefaultValParam=NULL, int whatFromN=1, int QDefaultValuesBefore=0, bool KeepIfRect=true){
     Array2DInsExtRow(                          this->y,              size,           N,                      x,             defaultVal,       whatFromN,        QDefaultValsBefore,                  true);
     delete defaultVal;
}
void Matrix_P2S::InsColumn(int N, double*x, int L, int FromN, int QDefaultBefore){
    Array2DSize size(this->GetQLines(), this->GetQColumns());
    double*defaultVal=new double;
    *defaultVal=0;
    //template<typename T>Array2DInsIneRow(T**&X, Array2DSize& size, int IneRowPosN, T*rowParam=NULL, int wholeRowL=0, int IfAfterLminIgnore0Stretch2=0, int whatFromN=1, int QDefaultValsBefore=0, T*DefaultValParam=NULL){
    Array2DInsIneRow(                    this->y,              size,              N,               x,               L,                                0,           FromN,          QDefaultBefore,     defaultVal);
    delete defaultVal;
}
void Matrix_P2S::InsColumn(int N, std::vector<double>x, int whatFromN, int QDefaultValsBefore){
    Array2DSize size(this->GetQLines(), this->GetQColumns());
    double*defaultVal=new double;
    *defaultVal=0;
    //template<typename T>Array2DInsIneRow(T**&X, Array2DSize& size, int IneRowPosN, std::vector<T>rowParam, int IfAfterLminIgnore0Stretch2=0, int whatFromN=1, int QDefaultValsBefore=0, T*DefaultValParam=NULL){
    Array2DInsIneRow(                     this->y,              size,              N,                      x,                                0,       whatFromN,       QDefaultValsBefore,             defaultVal);
}
void Matrix_P2S::DelLine(int N){
    Array2DSize size(this->GetQLines(), this->GetQColumns());
    //Array1DDelElementFromN(this->y, QLines, N);
    Array2DDelExtRow(this->y, size, N);
}
void Matrix_P2S::DelColumn(int N){
    Array2DSize size(this->GetQLines(), this->GetQColumns());
    //template<typename T>void Array2DDelIneRow(T**&X, Array2DSize&size, int N)
    Array2DDelIneRow(this->y, size, N);
}






// };


//-------------------------------------------------------------------------------------------------------------------------

   // class Matrix_P2L{
   //     double**y;
   //     int QLins, QCols;
   //     public:
Matrix_P2L::Matrix_P2L(int QLines, int QColumns, double val){
    std::cout<<"Matrix_P2L constr (dbl, QL, QC) works"<<std::endl;
    this->construct();
    this->SetSize(QLines, QColumns);
    for(int i=1; i<=QLines; i++){
        for(int j=1; j<=QColumns; j++){
            this->SetComponent(val, i, j);
        }
    }
}
Matrix_P2L::Matrix_P2L(double**y, int QLins, int QCols){
    std::cout<<"Matrix_P2L constr (**dbl, QL, QC) works"<<std::endl;
    this->construct();
    this->Set(y, QLins, QCols);
}
Matrix_P2L::Matrix_P2L(std::vector<std::vector<double>>y){
    std::cout<<"Matrix_P2L constr (vector<<dbl>>, QL, QC) works"<<std::endl;
    this->construct();
    this->Set(y);
}
Matrix_P2L::Matrix_P2L(double*y, int QLines, int QColumns){
    std::cout<<"Matrix_P2L constr (*dbl, QL, QC) works"<<std::endl;
    this->construct();
    this->Set(y, QLines, QColumns);
}
Matrix_P2L::Matrix_P2L(std::vector<double>y, int QColumns){
    std::cout<<"Matrix_P2L constr (vector<dbl>,  bool) works"<<std::endl;
    this->construct();
    this->Set(y, QColumns);
}
Matrix_P2L::Matrix_P2L(double x, double y, double z){
    std::cout<<"Matrix_P2L constr (x, y, z) works"<<std::endl;
    this->construct();
    this->Set(x, y, z);
}
//Matrix_P2L::Matrix_P2L(double x, int QLines, int QColumns){
//    std::cout<<"Matrix_P2L constr (dbl, QL, QC) works"<<std::endl;
//    this->construct();
//    this->Set(x, QLines, QColumns);
//
//}
//
Matrix_P2L::~Matrix_P2L(){
    std::cout<<"Matrix_P2L Destructor works"<<std::endl;
    this->SetNull();
}
//
void Matrix_P2L::SetNull(){
    if(this->y!=NULL){
        for(int i=1; i<=this->QLins; i++){
            delete[] this->y[i-1];
        }
        delete[] this->y;
    }
    this->y=NULL;
    this->QLins=0;
    this->QCols=0;
}
void Matrix_P2L::construct(){
    this->y=NULL;
    this->QLins=0;
    this->QCols=0;
    //
    this->SetOne(0);
}
Matrix_Prototype* Matrix_P2L::CreateMatrix(){
    return new Matrix_P2L();
}
Matrix_Prototype* Matrix_P2L::clone(){
    return this;
}
int Matrix_P2L::GetTypeN() const{
    return 121;//P/V, 1/2, S/L/A
}
//
void Matrix_P2L::SetSize(int QLins, int QCols){
    Array2DSize size1(this->GetQLines(), this->GetQColumns());
    Array2DSize size2(QLins, QCols);
    Array2DSetSize(this->y, size1, size2);//*TDflt, preserve=true
    this->QLins=QLins;
    this->QCols=QCols;
}
//
int Matrix_P2L::GetQLines() const { return this->QLins; }
int Matrix_P2L::GetQColumns() const { return this->QCols; }


void Matrix_P2L::Set(double**y, int QLins, int QCols){
    Matrix_Prototype::Set(y, QLins, QCols);
    this->QLins=QLins;
    this->QCols=QCols;
}
void Matrix_P2L::Set(std::vector<std::vector<double>>y){
    this->QLins=y.size(), this->QCols=y[1-1].size();
    Matrix_Prototype::Set(y);
}
void Matrix_P2L::Set(double*y, int QLines,  int QColumns){
    Matrix_Prototype::Set(y, QLins, QCols);
    this->QLins=QLins;
    this->QCols=QCols;
}
void Matrix_P2L::Set(std::vector<double>y, int QColumns){
    Matrix_Prototype::Set(y, QColumns);
    this->QCols=QColumns;
    this->QLins=y.size()/QColumns;
}
void Matrix_P2L::Set(double x, double y, double z){
    Matrix_Prototype::Set(x, y, z);
    this->QLins=3;
    this->QCols=1;
}
void Matrix_P2L::Set(double y, int QLines, int QColumns){
    Matrix_Prototype::Set(y, QLines, QColumns);
    this->QCols=QColumns;
    this->QLins=QLines;
}

double Matrix_P2L::GetComponent(int LineN, int ColN) const{
    double y=0;
    if(LineN>=1 && LineN<=this->QLins && ColN>=1 && ColN<=this->QCols){
        y=this->y[LineN-1][ColN-1];
    }
    return y;
}
void Matrix_P2L::SetComponent(double val, int LineN, int ColN){
    if(LineN>=1 && LineN<=this->QLins && ColN>=1 && ColN<=this->QCols){
        this->y[LineN-1][ColN-1]=val;
    }
}
//
void Matrix_P2L::SetLine(int N, double*x, int L, int FromN, int QDefaultBefore, bool DefaultValsAreZerosNotOwn){
    int QLins=this->GetQLines(), QCols=this->GetQColumns();
    std::vector<double>row;
    if(N>=1 && N<=QLins){
        if(DefaultValsAreZerosNotOwn){
            //
            row=Array1DGetSubArray_byLs_PreMark_v1(x, L, 0.0, QCols, FromN, QDefaultBefore);
        }else{
            row=Array1DGetSubArray_byLs_PreMark_v1(x, L, this->y[N-1], QCols, FromN, QDefaultBefore, QCols);
        }
        for(int i=1; i<=QCols; i++){
            this->y[N-1][i-1]=row[i-1];
        }
    }
}
void  Matrix_P2L::SetLine(int N, std::vector<double>x, int FromN, int QDfltValsBefore, bool DefaultValsAreOwnNotZeros){
    ////template<typename T> void Array2DSetExtRowN(std::vector<std::vector<T>>X, int ExtRowN, std::vector<T>rowParam, T*DefaultValParam=NULL, int FromN=1, int DefaultValuesBefore=0, bool KeepIfRect=false)
    ////template<typename T> void Array2DSetExtRowN(T**&X, const Array2DSize&size, int ExtRowN, std::vector<T>rowParam, T*DefaultValParam=NULL, int FromN=1, int DefaultValuesBefore=0, bool KeepIfRect=false)
    ////template<typename T> void Array2DSetExtRowN(T**&X, const Array2DSize&size, int ExtRowN, int wholeRowL=0, T*rowParam=NULL, T*DefaultValParam=NULL, int FromN=1, int DefaultValuesBefore=0, bool KeepIfRect=false)
    //Array2DSize size(this->GetQLines(), this->GetQColumns());
    //double*defaultVal=new double;
    //*defaultVal=0;
    ////Array2DSetExtRowN(T**&X,  const Array2DSize&size, int ExtRowN, std::vector<T>rowParam, T*DefaultValParam=NULL, int FromN=1, int DefaultValuesBefore=0, bool KeepIfRect=false)
    //Array2DSetExtRowN(this->y,                    size,           N,                      x,             defaultVal,      FromN,           QDfltValsBefore, true);
    //delete defaultVal;
    SetLine(N, x.data(), x.size(), FromN, QDfltValsBefore, DefaultValsAreOwnNotZeros);
}
void  Matrix_P2L::SetColumn(int N, double*x, int L, int FromN, int QDefaultBefore, bool DefaultValsAreOwnNotZeros){
    int QLins=this->GetQLines(), QCols=this->GetQColumns();
    std::vector<double>row;
    if(N>=1 && N<=QCols){
        if(DefaultValsAreOwnNotZeros=false){
            row=Array1DGetSubArray_byLs_PreMark_v1(x, L, 0.0, QCols, FromN, QDefaultBefore);//vikts: if 0 alt 0.0, S int S ak numoz et di error no fn co tal params oq
        }else{
            row=Array1DGetSubArray_byLs_PreMark_v1(x, L, this->y[N-1], QLins, FromN, QDefaultBefore, QLins);
        }
        for(int i=1; i<=QCols; i++){
            this->y[i-1][N-1]=row[i-1];
        }
    }
}
void  Matrix_P2L::SetColumn(int N, std::vector<double>x, int FromN, int QDefaultBefore, bool DefaultValsAreOwnNotZeros){
    ////template<typename T> void Array2DSetIneRowN(T**&X, Array2DSize&size, int IneRowN, int wholeRowL=0, T*rowParam=NULL, T*DefaultValParam=NULL, int FromN=1, int DefaultValuesBefore=0, int IfNGTLmin_Ignore0DoNil1Add2Srtetch3=1)
    ////template<typename T> void Array2DSetIneRowN(T**&X, Array2DSize&size, int IneRowN, std::vector<T>rowParam, T*DefaultValExt=NULL, int FromN=1, int DefaultValuesBefore=0, int IfNGTLmin_Ignore0DoNil1Add2Srtetch3=1)
    ////template<typename T> void Array2DSetIneRowN(std::vector<std::vector<T>>&X, std::vector<T>rowParam, int IneRowN, T*DefaultValExt=NULL, int FromN=1, int DefaultValuesBefore=0, int IfNGTLmin_Ignore0DoNil1Add2Srtetch3=1)
    //Array2DSize size(this->GetQLines(), this->GetQColumns());
    //double*defaultVal=new double;
    //*defaultVal=0;
    ////template<typename T>Array2DSetIneRow(T**&X, Array2DSize& size, int IneRowN, std::vector<T>rowParam, int whatFromN=1, int QDefaultValsBefore=0, T*DefaultValParam=NULL, bool ExistingInsteadOfDefault=true, bool ExistingOnlyForGTLminNotIgnore=false, bool LastPossIfNotRectAtPos0NotIgnore=false){
    //Array2DSetIneRow(                    this->y,              size,           N,                      x, whatN1,                   QDfltValsBefore, defaultVal);
    //delete defaultVal;
    SetColumn(N, x.data(), x.size(), FromN, QDefaultBefore, DefaultValsAreOwnNotZeros);
}
//
void Matrix_P2L::AddLine(double*x, int L, int whatN1, int QDefaultValsBefore){
    int QLins=this->GetQLines(), QCols=this->GetQColumns();
    Array2DSize size(QLins, QCols);
    std::vector<double>row;
    //double *line=NULL, *ptr=NULL, Qv=0;
    bool isRect=true;
    row=Array1DGetSubArray_byLs_PreMark_v1(x, L, 0.0, QCols, whatN1, QDefaultValsBefore);
    Array2DAddExtRow(this->y, size, row.data(), row.size(), whatN1, QDefaultValsBefore, isRect);
    this->QLins=size.GetQExtRows();
    //this->QCols=remains
}
void Matrix_P2L::AddLine(std::vector<double>x, int whatN1, int QDefaultValsBefore){
   int QLines=this->GetQLines(), QColumns=this->GetQColumns(), L=x.size();
    Array2DSize size(QLines, QColumns);
    //template<typename T>void Array2DAddExtRow(T**&X, Array2DSize&size, T*rowParam, int wholeRowL=0, int whereFromN=1, int whatFromN=1, bool RectNotVar=false, T*DfltValParam=NULL, TValsShowHide*vsh=NULL){//44 //better but not checked
    double*defaultVal=new double;
    *defaultVal=0;
    if((QLines==1 && QColumns==1 && L>1) || ((QLines==0 || QColumns==0)&&L>0)){
        this->SetNull();
    }
    Array2DAddExtRow(this->y, size, x.data(), x.size(), 1, whatN1, true, defaultVal);
    delete defaultVal;
}
void Matrix_P2L::AddColumn(double*x, int L, int whatN1, int QDefaultValsBefore){
    //template<typename T>void Array2DAddIneRow(T**&X, Array2DSize&size, T*rowParam, int wholeRowL=0, T*DfltValParam=NULL, bool AddToEachNotStretchToMax=true, int whatFromN=1, int whereFromN=1
    int QLines=this->GetQLines(), QColumns=this->GetQColumns();
    Array2DSize size(QLines, QColumns);
    double*defaultVal=new double;
    *defaultVal=0;
    if((QLines==1 && QColumns==1 && L>1) || ((QLines==0 || QColumns==0)&&L>0)){
        this->SetNull();
    }
    Array2DAddExtRow(this->y, size, x, L, 1, whatN1, true, defaultVal);
    delete defaultVal;
}
void Matrix_P2L::AddColumn(std::vector<double>x, int FromN, int QDefaultValsBefore){
    //template<typename T>void Array2DAddIneRow(T**&X, Array2DSize&size, T*rowParam, int wholeRowL=0, T*DfltValParam=NULL, bool AddToEachNotStretchToMax=true, int whatFromN=1, int whereFromN=1){
    int QLines=this->GetQLines(), QColumns=this->GetQColumns(), L=x.size();
    Array2DSize size(QLines, QColumns);
    double*defaultVal=new double;
    *defaultVal=0;
    if( (QLines==1 && QColumns==1 && L>1) || ((QLines==0||QColumns==0)&&L>0) ){
        this->SetNull();
    }
    //template<typename T>Array2DAddIneRow(T**&X, Array2DSize& size, std::vector<T>rowParam, int IfNonRectIgnore0Add1Stretch2=0, int whatFromN=1, int QDefaultValsBefore=0, T*DefaultValParam=NULL){
    Array2DAddIneRow(                    this->y,              size,                      x,                                  0,           FromN,                        0,  defaultVal);
    delete defaultVal;
}

void Matrix_P2L::InsLine(int N, double*x, int L, int FromN, int QDefaultBefore){
    Array2DSize size(this->GetQLines(), this->GetQColumns());
    double*defaultVal=new double;
    *defaultVal=0;
    //template<typename T> void Array2DInsExtRow(T**&X,  Array2DSize&size, int ExtRowN, int wholeRowL=0, T*rowParam=NULL, T*DefaultValParam=NULL, int whatFromN=1, int QDefaultValuesBefore=0, bool KeepIfRect=true, TValsShowHide*vsh=NULL){
    Array2DInsExtRow(                          this->y,              size,           N,             L,                 x,             defaultVal,           FromN,             QDefaultBefore,                 true);
    delete defaultVal;
}
void Matrix_P2L::InsLine(int N, std::vector<double>x, int whatFromN, int QDefaultValsBefore){
     Array2DSize size(this->GetQLines(), this->GetQColumns());
     double*defaultVal=new double;
     *defaultVal=0;
     //template<typename T> void Array2DInsExtRow(T**&X,  Array2DSize&size, int ExtRowN, std::vector<T>rowParam, T*DefaultValParam=NULL, int whatFromN=1, int QDefaultValuesBefore=0, bool KeepIfRect=true){
     Array2DInsExtRow(                          this->y,              size,           N,                      x,             defaultVal,       whatFromN,        QDefaultValsBefore,                  true);
     delete defaultVal;
}
void Matrix_P2L::InsColumn(int N, double*x, int L, int FromN, int QDefaultBefore){
    Array2DSize size(this->GetQLines(), this->GetQColumns());
    double*defaultVal=new double;
    *defaultVal=0;
    //template<typename T>Array2DInsIneRow(T**&X, Array2DSize& size, int IneRowPosN, T*rowParam=NULL, int wholeRowL=0, int IfAfterLminIgnore0Stretch2=0, int whatFromN=1, int QDefaultValsBefore=0, T*DefaultValParam=NULL){
    Array2DInsIneRow(                    this->y,              size,              N,               x,               L,                                0,           FromN,          QDefaultBefore,     defaultVal);
    delete defaultVal;
}
void Matrix_P2L::InsColumn(int N, std::vector<double>x, int whatFromN, int QDefaultValsBefore){
    Array2DSize size(this->GetQLines(), this->GetQColumns());
    double*defaultVal=new double;
    *defaultVal=0;
    //template<typename T>Array2DInsIneRow(T**&X, Array2DSize& size, int IneRowPosN, std::vector<T>rowParam, int IfAfterLminIgnore0Stretch2=0, int whatFromN=1, int QDefaultValsBefore=0, T*DefaultValParam=NULL){
    Array2DInsIneRow(                     this->y,              size,              N,                      x,                                0,       whatFromN,       QDefaultValsBefore,             defaultVal);
}
void Matrix_P2L::DelLine(int N){
    Array2DSize size(this->GetQLines(), this->GetQColumns());
    //Array1DDelElementFromN(this->y, QLines, N);
    Array2DDelExtRow(this->y, size, N);
}
void Matrix_P2L::DelColumn(int N){
    Array2DSize size(this->GetQLines(), this->GetQColumns());
    //template<typename T>void Array2DDelIneRow(T**&X, Array2DSize&size, int N)
    Array2DDelIneRow(this->y, size, N);
}



// }


//----------------------------------------------------------------------------
//class Matrix_V2S{
//    std::vector<std::vector<double>>y;
//    public:
Matrix_V2S::Matrix_V2S(int QLins, int QCols, double val){
    std::cout<<"Matrix_V2S constr (QL, QC) works"<<std::endl;
    this->construct();
    //this->SetSize(QLins, QCols);
    this->Set(val, QLins, QCols);
}
Matrix_V2S::Matrix_V2S(double**y, int QLins, int QCols){
    std::cout<<"Matrix_V2S constr (**dbl, QL, QC) works"<<std::endl;
    this->construct();
    this->Set(y, QLins, QCols);
}
Matrix_V2S::Matrix_V2S(std::vector<std::vector<double>>y){
    std::cout<<"Matrix_V2S constr (vector<<dbl>>, QL, QC) works"<<std::endl;
    this->construct();
    this->Set(y);
}
Matrix_V2S::Matrix_V2S(double*y, int QLines, int QColumns){
    std::cout<<"Matrix_V2S constr (*dbl, QL, QC) works"<<std::endl;
    this->construct();
    this->Set(y, QLines, QColumns);
}
Matrix_V2S::Matrix_V2S(std::vector<double>y, int QColumns){
    std::cout<<"Matrix_V2S constr (vector<dbl>,  bool) works"<<std::endl;
    this->construct();
    this->Set(y, QColumns);
}
Matrix_V2S::Matrix_V2S(double x, double y, double z){
    std::cout<<"Matrix_V2S constr (x, y, z) works"<<std::endl;
    this->construct();
    this->Set(x, y, z);
}
//
Matrix_V2S::~Matrix_V2S(){
    //this->SetNull();//possible but no need
}
//
void Matrix_V2S::SetNull(){
    this->y.clear();
}
void ::Matrix_V2S::construct(){
    this->y.clear();
    //
    this->SetOne(0);
}
Matrix_Prototype* Matrix_V2S::CreateMatrix(){
    return new Matrix_V2S();
}

Matrix_Prototype* Matrix_V2S::clone(){
    return this;
}
int Matrix_V2S::GetTypeN() const{
    return 221;//P/V, 1/2, S/L/A
}


void Matrix_V2S::SetSize(int QLins, int QCols){
    //Array2DSize oldSize, newSize;
    //newSize.Set(QLins, QCols);
    Array2DSetSize(this->y, QLins, QCols);
}

int Matrix_V2S::GetQLines() const { return this->y.size(); }
int Matrix_V2S::GetQColumns() const { return this->y[1-1].size(); }
//
void Matrix_V2S::Set(double**y, int QLins, int QCols){
    Matrix_Prototype::Set(y, QLins, QCols);
    //this->QLins=QLins;
    //this->QCols=QCols;
}
void Matrix_V2S::Set(std::vector<std::vector<double>>y){
    //this->QLins=y.size(), this->QCols=y[1-1].size();
    Matrix_Prototype::Set(y);
}
void Matrix_V2S::Set(double*y, int QLines,  int QColumns){
    Matrix_Prototype::Set(y, QLines, QColumns);
    //this->QLins=QLins;
    //this->QCols=QCols;
}
void Matrix_V2S::Set(std::vector<double>y, int QColumns){
    Matrix_Prototype::Set(y, QColumns);
    //this->QCols=QColumns;
    //this->QLins=y.size()/QColumns;
}
void Matrix_V2S::Set(double x, double y, double z){
    Matrix_Prototype::Set(x, y, z);
    //this->QLins=3;
    //this->QCols=1;
}
void Matrix_V2S::Set(double y, int QLines, int QColumns){
    Matrix_Prototype::Set(y, QLines, QColumns);
    //this->QCols=QColumns;
    //this->QLins=QLines;
}

//
double Matrix_V2S::GetComponent(int LineN, int ColN) const{
    double y=0;
    int QLins=this->GetQLines(), QCols=this->GetQColumns();
    if(LineN>=1 && LineN<=QLins && ColN>=1 && ColN<=QCols){
        y=this->y[LineN-1][ColN-1];
    }
    return y;
}
void Matrix_V2S::SetComponent(double val, int LineN, int ColN){
    int QLins=this->GetQLines(), QCols=this->GetQColumns();
    if(LineN>=1 && LineN<=QLins && ColN>=1 && ColN<=QCols){
        this->y[LineN-1][ColN-1]=val;
    }
}
//
void  Matrix_V2S::SetLine(int N, double*x, int Q, int whatN1, int QDfltValsBefore, bool DefaulsAreZerosNotOwnVals){
    double*defaultVal=NULL;  //defaultVal=new double; (*(defaultVal))=0; // no need let em be random, af D wi be 0
    int QLines=this->GetQLines(), QColumns=this->GetQColumns();//, Ladded;
    int LreqCalcd=0, preN1_=0, preN2_=0, ownN1_=0, ownN2=0, ForOwnN1_=0, ForOwnN2=0, postN1=0, postN2_=0, QDefaultValuesBefore=0;
    std::vector<double> row1;
    Array1DAssign(row1, x, Q);
    //std::vector<double> GetNumericSubarray_byLs(std::vector<double> x, int FromN=1, int QDefaultValuesBefore=0, int LreqGiven=0);
    std::vector<double> row2;
    std::vector<double>* VP=NULL;
    if(N>=1 && N<=QLines){
        if(DefaulsAreZerosNotOwnVals){
            row2= GetNumericSubarray_byLs(row1, whatN1, QDfltValsBefore, QColumns);
        }else{
            row2 = GetNumericSubarray_byLs(row1, this->y[N-1], whatN1, QDfltValsBefore, QColumns);
        }
        Array2DSetExtRowN(this->y, N, row2);
    }
}
void  Matrix_V2S::SetLine(int N, std::vector<double>x, int whatN1, int QDfltValsBefore, bool DefaulsAreZerosNotOwnVals){
    this->SetLine(N, x.data(), x.size(), whatN1, QDfltValsBefore, DefaulsAreZerosNotOwnVals);
}
void  Matrix_V2S::SetColumn(int N, double*x, int Q, int whatN1, int QDfltValsBefore, bool DefaulsAreOwnValsNotZeros){
    std::vector<double> row1;
    Array1DAssign(row1, x, Q);
    std::vector<double> row2;// = DefaulsAreZerosNotOwnVals==true? GetNumericSubarray_byLs(row1, whatN1, QDfltValsBefore, this->GetQLines()); : GetNumericSubarray_byLs(row1, this->y[N-1], whatN1, QDfltValsBefore, this->GetQLines());
    std::vector<double> ColumnN=this->GetArrayOfColN(N);
    int QColumns=this->GetQColumns();
    if(N>=1 && N<=QColumns){
        if(DefaulsAreOwnValsNotZeros){
            row2 = GetNumericSubarray_byLs(row1, whatN1, QDfltValsBefore, this->GetQLines());
        }else{
            row2 = GetNumericSubarray_byLs(row1, this->y[N-1], whatN1, QDfltValsBefore, this->GetQLines());
        }
        //template<typename T>Array2DSetIneRow(std::vector<std::vector<T>>&X, int IneRowN, std::vector<T>rowParam, int whatFromN=1, int QDefaultValsBefore=0, T*DefaultValParam=NULL, bool ExistingInsteadOfDefault=true, bool ExistingOnlyForGTLminNotIgnore=false, bool LastPossIfNotRectAtPos0NotIgnore=false){
        Array2DSetIneRow(this->y, N, row2);
    }
}
void  Matrix_V2S::SetColumn(int N, std::vector<double>x, int whatN1, int QDfltValsBefore, bool DefaulsAreOwnValsNotZeros){
    this->SetColumn(N, x.data(), x.size(),  whatN1, QDfltValsBefore, DefaulsAreOwnValsNotZeros);
}
//
void Matrix_V2S::AddLine(double*x, int L, int whatN1, int QDefaultValsBefore){
    int QLines=this->GetQLines(), QColumns=this->GetQColumns();
    if((QLines==1 && QColumns==1 && L>1) || ((QLines==0 || QColumns==0)&&L>0)){
        this->SetNull();
    }
    std::vector<double>row0, row;
    Array1DAssign(row0, x, L);
    row=GetNumericSubarray_byLs(row0, whatN1, QDefaultValsBefore, this->GetQColumns());
    Array2DAddExtRow(this->y, row);
}
void Matrix_V2S::AddLine(std::vector<double>x, int whatN1, int QDefaultValsBefore){
    int QLines=this->GetQLines(), QColumns=this->GetQColumns(), L=x.size();
    if((QLines==1 && QColumns==1 && L>1) || ((QLines==0 || QColumns==0)&&L>0)){
        this->SetNull();
    }
    std::vector<double> row;
    row=GetNumericSubarray_byLs(x, whatN1, QDefaultValsBefore, this->GetQColumns());
    Array2DAddExtRow(this->y, row);
}
void Matrix_V2S::AddColumn(double*x, int L, int whatN1, int QDefaultValsBefore){
    int QLines=this->GetQLines(), QColumns=this->GetQColumns();
    if((QLines==1 && QColumns==1 && L>1) || ((QLines==0 || QColumns==0)&&L>0)){
        this->SetNull();
    }
    std::vector<double>row0, row;
    Array1DAssign(row0, x, L);
    row=GetNumericSubarray_byLs(row0, whatN1, QDefaultValsBefore, this->GetQLines());
    Array2DAddIneRow(this->y, row);
}
void Matrix_V2S::AddColumn(std::vector<double>x, int whatN1, int QDefaultValsBefore){
    int QLines=this->GetQLines(), QColumns=this->GetQColumns(), L=x.size();
    if((QLines==1 && QColumns==1 && L>1) || ((QLines==0 || QColumns==0)&&L>0)){
        this->SetNull();
    }
    std::vector<double>row;
    row=GetNumericSubarray_byLs(x, whatN1, QDefaultValsBefore, this->GetQLines());
    Array2DAddIneRow(this->y, row);
}

void Matrix_V2S::InsLine(int N, double*x, int Q, int whatN1, int QDefaultValsBefore){
    int QLines=this->GetQLines(), QColumns=this->GetQColumns();
    std::vector<double> row, row0;
    Array1DAssign(row0, x, Q);
    row=GetNumericSubarray_byLs(row0, whatN1, QDefaultValsBefore, this->GetQColumns());
    Array2DInsExtRow(this->y, N, row);
}
void Matrix_V2S::InsLine(int N, std::vector<double>x, int whatN1, int QDefaultValsBefore){
    int QLines=this->GetQLines(), QColumns=this->GetQColumns(), L=x.size();
    std::vector<double> row;//, row0;
    //Array1DAssign(row0, x, Q);
    row=GetNumericSubarray_byLs(x, whatN1, QDefaultValsBefore, this->GetQColumns());
    Array2DInsExtRow(this->y, N, row);
}
void Matrix_V2S::InsColumn(int N, double*x, int Q, int whatN1, int QDefaultValsBefore){
    int QLines=this->GetQLines(), QColumns=this->GetQColumns();
    std::vector<double> row, row0;
    Array1DAssign(row0, x, Q);
    row=GetNumericSubarray_byLs(row0, whatN1, QDefaultValsBefore, this->GetQLines());
    Array2DInsIneRow(this->y, N, row);
}
void Matrix_V2S::InsColumn(int N, std::vector<double>x, int whatN1, int QDefaultValsBefore){
    int QLines=this->GetQLines(), QColumns=this->GetQColumns(), L=x.size();
    std::vector<double> row;
    row=GetNumericSubarray_byLs(x, whatN1, QDefaultValsBefore, this->GetQLines());
    Array2DInsIneRow(this->y, N, row);
}
void Matrix_V2S::DelLine(int N){
   Array2DDelExtRow(this->y, N);
}
void Matrix_V2S::DelColumn(int N){
    Array2DDelIneRow(this->y, N);
}

//
/*
void Matrix_V2S::Transpose(){
   Array2DTranspose(this->y);
}
Matrix_V2S Matrix_V2S::TransposeTo(){
    Matrix_V2S M=*this;
    M.Transpose();
    return M;
}
//
void Matrix_V2S::SetX(double val){
    if(this->GetQLines()>=1){
        this->SetComponent(val, 1, 1);
    }
}

void Matrix_V2S::SetY(double val){
    if(this->GetQLines()>=2){
        this->SetComponent(val, 2, 1);
    }
}
void Matrix_V2S::SetZ(double val){
    if(this->GetQLines()>=3){
        this->SetComponent(val, 3, 1);
    }
}
double Matrix_V2S::GetX(){
    return this->y[1-1][1-1];
}

double Matrix_V2S::GetY(){
   double y;
   if(this->GetQLines()>=2){
       y=this->y[2-1][1-1];
   }
   return y;
}

double Matrix_V2S::GetZ(){
    double y;
    if(this->GetQLines()>=3){
        y=this->y[3-1][1-1];
    }
    return y;
}


Matrix_V2S Matrix_V2S::SubMatrix(int Line1N, int Col1N, int Line2N, int Col2N){
    Matrix_V2S M;
    int QLins=this->GetQLines(), QCols=this->GetQColumns(), QLinsNew=0, QColsNew=0, LineN, ColN;
    double val;
    if(Line1N>=1 && Line1N<=QLins && Line2N>=1 && Line2N<=QLins && Col1N>=1 && Col1N<=QCols && Col2N>=1 && Col2N<=QCols){
        if(Line1N<=Line2N  && Col1N<=Col2N){
            QLinsNew=Line2N-Line1N+1;
            QColsNew=Col2N-Col1N+1;
            M.SetSize(QLinsNew, QColsNew);
            //for(int i=Line1N; i<=Line2N; i++){
            for(int i=1; i<=QLinsNew; i++){
                //LineN=i-Line1N+1;
                LineN=i+Line1N-1;
                for(int j=Col1N; j<=Col2N; j++){
                    ColN=Col1N+j-1;
                    //ColN=j-Col1N+1;
                    val=this->GetComponent(LineN, ColN);
                    //val=this->GetElement(i, j);
                    M.SetComponent(val, i, j);
                    //M.SetElement(val, LineN, ColN);
                }
            }
        }else if(Line1N<=Line2N  && Col1N>Col2N){
            QLinsNew=Line2N-Line1N+1;
            QColsNew=Col1N-Col2N+1;
            M.SetSize(QLinsNew, QColsNew);
            for(int i=1; i<=QLinsNew; i++){
                LineN=i+Line1N-1;
                for(int j=1; j>=QColsNew; j++){
                    ColN=Col2N+j-1;//(changes to Col1N)
                    val=this->GetComponent(LineN, ColN);
                    M.SetComponent(val, i, j);
                }
            }
        }else if(Line1N>=Line2N  && Col1N<=Col2N){
            QLinsNew=Line1N-Line2N+1;
            QColsNew=Col2N-Col1N+1;
            M.SetSize(QLinsNew, QColsNew);
            for(int i=1; i<=QLinsNew; i++){
                LineN=Line2N+i-1;//(changes to Line1N)
                for(int j=1; j>=QColsNew; j++){
                    ColN=Col1N+j-1;
                    val=this->GetComponent(LineN, ColN);
                    M.SetComponent(val, i, j);
                }
            }
        }else{//both 1 > 2
            QLinsNew=Line1N-Line2N+1;
            QColsNew=Col1N-Col2N+1;
            M.SetSize(QLinsNew, QColsNew);
            for(int i=1; i<=QLinsNew; i++){
                LineN=Line2N+i-1;//(changes to Line1N)
                for(int j=1; j>=QColsNew; j++){
                    ColN=Col2N+j-1;//(changes to Col1N)
                    val=this->GetComponent(LineN, ColN);
                    M.SetComponent(val, i, j);
                }
            }
        }
    }
    return M;
}


Matrix_V2S Matrix_V2S::MinorTo(int LineN, int ColN) const
       {
           Matrix_V2S M;
           M=*this;
           // M = this.DelLineTo(LineN);
           // M = M.DelColTo(ColN);
           M.DelLine(LineN);
           M.DelColumn(ColN);
           return M;
       }

       void Matrix_V2S::Minor(int LineN, int ColumnN)
       {
           this->DelColumn(ColumnN);
           this->DelLine(LineN);
       }
       double Matrix_V2S::AlgSuppl(int LineN, int ColumnN) const
       {
           double d, r;
           int QLines=this->GetQLines(), QColumns=this->GetQColumns();
           Matrix_V2S M;// in |C# wa M = new Matrix_VB();
           if ((LineN < 1) || (LineN > QLines) || (ColumnN < 1) || (ColumnN > QColumns))
           {
               r = 666;
           }
           else
           {
               M = *this;
               M = MinorTo(LineN, ColumnN);
               d = M.Determinant();
               r = d;
               if ((LineN + ColumnN)%2 != 0)
               {
                   r *= (-1);
               }
           }
           return r;
       }
        double Matrix_V2S::Determinant() const
       {
           double r, d, AlgSpl;
           int RowN;
           int QColumns=this->GetQColumns(), QLines=this->GetQLines();
           if (QColumns != QLines) return 666;
           else
           {
               if (QLines == 1)
               {
                   r = this->y[1 - 1][1 - 1];
               }
               else if (QLines == 2)
               {
                   r = this->y[1 - 1][1 - 1] * this->y[2-1][2-1] - this->y[1 - 1][2 - 1] * this->y[2 - 1][1 - 1];
               }
               else if (QLines == 3)
               {
                   r = this->y[1-1][1-1] * this->y[2-1][2-1] * this->y[3-1][3-1] + this->y[1-1][3-1] * this->y[2-1][1-1] * this->y[3-1][2-1] + this->y[3-1][1-1] * this->y[1-1][2-1] * this->y[2-1][3-1] -
                       this->y[1-1][3-1] * this->y[2-1][2-1] * this->y[3-1][1-1] - this->y[1-1][1-1] * this->y[2-1][3-1] * this->y[3-1][2-1] - this->y[3-1][3-1] * this->y[1-1][2-1] * this->y[2-1][1-1];
               }
               else
               {
                   RowN = 1;
                   r = 0;
                   for (int j = 1; j <= QColumns; j++)
                   {
                       //r += GetComponent(RowN, j) * AlgSuppl(RowN, j);
                       r += GetComponent(RowN, j) * AlgSuppl(RowN, j);
                   }
               }
           }
           return r;
       }

       Matrix_V2S Matrix_V2S::UnionMatrix() const
       {
           int QLines = this->GetQLines(), QColumns = this->GetQColumns();
           int TQLines = QColumns, TQColumns = QLines;
           //Matrix_VB MT = new Matrix_VB();
           //MT = this.TransposeTo();
           Matrix_V2S MUR;//in C# = new Matrix_VB();
           //double[][] c;
           double CurAlgSuppl;//, det = this.Determinant();
          // c = new double[TQLines][];
           //for (int i = 1; i <= TQLines; i++)
           //{
           //    c[i - 1] = new double[TQColumns];
           //}
           MUR.SetSize(TQLines, TQColumns);
           //
           for (int i = 1; i <= QLines; i++)
           {
               for (int j = 1; j <= QColumns; j++)
               {
                   CurAlgSuppl = this->AlgSuppl(i, j);
                   //c[j - 1][i - 1] = CurAlgSuppl;
                   MUR.SetComponent(CurAlgSuppl, j, i);
               }
           }
           //
           //MUR = new Matrix_VB(c, TQLines, TQColumns);
           return MUR;
       }
       Matrix_V2S Matrix_V2S::InverseMatrix() const
       {
           Matrix_V2S MIR;
           double det = this->Determinant();
           MIR = this->UnionMatrix();
           if (det != 0) MIR = MIR * (1 / det);
           return MIR;
       }
*/


/*QString Matrix_V::ToString(){
    QString str="", cur;
    str="[[";
    int QLins=this->GetQLines(), QCols=this->GetQColumns();
    for(int i=1; i<=QLins-1; i++){
        str+="[";
        for(int j=1; j<=QCols-1; j++){
            cur.setNum(this->y[i-1][j-1]);
            str+=cur;
            str+=", ";
            cur="";
        }
        cur.setNum(this->y[i-1][QCols-1]);
        str+=cur;
        str+="];";
    }
    str+="[";
    //LastLine
    for(int j=1; j<=QCols-1; j++){
        cur.setNum(this->y[QLins-1][j-1]);
        str+=cur;
        str+=", ";
        cur="";
    }
    //LastElt
    cur.setNum(this->y[QLins-1-1][QCols-1]);
    str+=cur;
    str+="]]";
    return str;
}
std::string Matrix_V::ToStdString(){
    std::string ss;
    QString qs;
    qs=this->ToString();
    ss=qs.toStdString();
    return ss;
}*/

/*
Matrix_V2S Matrix_V2S::GivensRotationTransformMatrix(int N1, int N2, Matrix_V2S*M, TValsShowHide*vsh){
    int QL;
    Matrix_V2S T(QL,QL);
    if(M==NULL)M=this;
    QL=M->GetQLines();
    double a11=M->GetComponent(N1, N1), a21=M->GetComponent(N2, N1), c12, s12, denom;
    for(int i=1; i<=QL; i++){
        for(int j=1; j<=QL; j++){
            T.SetComponent(0, i, j);
        }
    }
    for(int i=1; i<=QL; i++){
        T.SetComponent(1, i, i);
    }
    if(a21==0){
        c12=1;
        s12=0;
    }else{
        denom=sqrt(a11*a11+a21*a21);
        c12=a11/denom;
        s12=a21/denom;
    }
    writeln(vsh, " T=" +M->ToString());
    writeln(vsh, " a11=" +FloatToStr(a11) + " a21=" +FloatToStr(a21)+ " c12=" + FloatToStr(c12) + " s12=" + FloatToStr(s12));
    T.SetComponent(c12, N1, N1);
    T.SetComponent(c12, N2, N2);
    T.SetComponent(s12, N1, N2);
    T.SetComponent(-s12, N2, N1);
    return T;
}
//void QRDecomposition(Matrix_VB&Q, Matrix_VB&R,  Matrix_VB*M=NULL, TValsShowHide* vsh = NULL);
void Matrix_V2S::QRDecomposition(Matrix_V2S&Q, Matrix_V2S&R,  Matrix_V2S*M, TValsShowHide* vsh)
{
    writeln(vsh, "QRDecomposition starts working");
    if (M == NULL) M = this;
    int QL = M->GetQLines();
    std::vector<int>ns;// = new int[QL-1];
    //int countTsInLine;
    int jn, ii, jj;
    Q.SetSize(QL, QL);
    for (int i = 1; i <= QL; i++)
    {
        Q.SetComponent(1, i, i);
    }
    std::vector<std::vector<Matrix_V2S>>Ts;
    std::vector<Matrix_V2S>T1;
    Matrix_V2S M1;
    for (int i = 1; i <= QL - 1; i++)
    {
        ns.push_back(QL - 1 - i + 1);
        T1.clear();
        for(int j=1; j<=ns[i - 1]; j++){
            T1.push_back(M1);
        }
        Ts.push_back(T1);
    }
    std::vector<Matrix_V2S>TL;
    Matrix_V2S T;//, R=new Matrix(QL, QL);
    R.Assign(*M);
    writeln(vsh, "M="+M->ToString());
    writeln(vsh, "ini Q=" + Q.ToString());
    writeln(vsh, "ini R=" + R.ToString());
    for (int i = 1; i <= QL - 1; i++)
    {
        //countTsInLine = 0;
        for (int j = i + 1; j <= QL; j++)
        {
            //countTsInLine++; //in AddToVector
            writeln(vsh, "i=" + IntToStr(i)+" j="+IntToStr(j));
            T = R.GivensRotationTransformMatrix(i, j, NULL, vsh);
            writeln(vsh, "T=" + T.ToString());
            R = T * R;
            writeln(vsh, "cur R=" + R.ToString());
            jn = j - i;
            Ts[i - 1][jn - 1] = T;
        }
    }
    //R is completed, beginning wit Q
    for(int i=1; i<=QL-1; i++){
        ii = QL - 1 - i + 1;
        for (int j = 1; j <= ns[ii - 1]; j++)
        {
            jj = ns[ii - 1] - j + 1;
            T = Ts[ii - 1][jj - 1];
            Q = Q * T;
        }
    }
    Q.Transpose();
    //
    writeln(vsh, "Answer:");
    writeln(vsh,"Q="+Q.ToString()+" R="+R.ToString());
    writeln(vsh, "QRDecomposition finishes working");

}//fn QR decomp


Matrix_V2S Matrix_V2S::CalcMatrixOfDirCossByEulersAngles (Matrix_V2S MEulerAngles){
    Matrix_V2S MR(3, 3);
    double gamma = MEulerAngles.GetX() * M_PI / 180,
                       psi = MEulerAngles.GetY() * M_PI / 180,
                       theta = MEulerAngles.GetZ() * M_PI / 180;
                MR.SetComponent(cos(psi) * cos(theta), 1, 1);
                MR.SetComponent(sin(psi) * sin(gamma) - cos(psi) * sin(theta) * cos(gamma), 1, 2);
                MR.SetComponent(sin(psi) * cos(gamma) + cos(psi) * sin(theta) * sin(gamma), 1, 3);
                MR.SetComponent(sin(theta), 2, 1);
                MR.SetComponent(cos(theta) * cos(gamma), 2, 2);
                MR.SetComponent(-cos(theta) * sin(gamma), 2, 3);
                MR.SetComponent(-sin(psi) * cos(theta), 3, 1);
                MR.SetComponent(cos(psi) * sin(gamma) + sin(psi) * sin(theta) * cos(gamma), 3, 2);
                MR.SetComponent(cos(psi) * cos(gamma) - sin(psi) * sin(theta) * sin(gamma), 3, 3);
                return MR;
}

Matrix_V2S Matrix_V2S::CalcCoordsTransformFromOldCSToNew(Matrix_V2S*CoordsInOldCSParam, Matrix_V2S*EulerAnglesParam, Matrix_V2S*CSOriginOfNewCSInOldParam, TValsShowHide*vsh){
    Matrix_V2S MCoordsInOldCS(3, 1), MEulerAngles (3, 1),MCSOriginOfNewCSInOld(3, 1);
    Matrix_V2S MDirCoss(3, 3), MCoordsInNewCS(3, 1), M1(3, 3), M2(3, 3), M3(3, 3);//, M4(3, 3);
    writeln(vsh, "CalcCoordsTransformFromOldCSToNew starts working");
    if(CoordsInOldCSParam!=NULL){
        MCoordsInOldCS=(*(CoordsInOldCSParam));
    }
    writeln(vsh, "MCoordsInOldCS");
    MCoordsInOldCS.StdCOut2D();
    if(EulerAnglesParam!=NULL){
        MEulerAngles=(*(EulerAnglesParam));
    }
    writeln(vsh, "MEulerAngles");
    MCoordsInOldCS.StdCOut2D();
    if(CSOriginOfNewCSInOldParam!=NULL){
        MCSOriginOfNewCSInOld=(*(CSOriginOfNewCSInOldParam));
    }
    writeln(vsh, "MCSOriginOfNewCSInOld");
    MCSOriginOfNewCSInOld.StdCOut2D();
    MDirCoss=this->CalcMatrixOfDirCossByEulersAngles(MEulerAngles);
    writeln(vsh, "MDirCoss");
    MDirCoss.StdCOut2D();
    M1=MDirCoss.InverseMatrix();
    writeln(vsh, "MDirCoss Inverted");
    M1.StdCOut2D();
    M2=MCoordsInOldCS-MCSOriginOfNewCSInOld;
    writeln(vsh, "M2=MCoordsInOldCS-MCSOriginOfNewCSInOld");
    M2.StdCOut2D();
    writeln(vsh, "M3=M1*M2");
    M3=M1*M2;
    M3.StdCOut2D();
    MCoordsInNewCS=M3;
    writeln(vsh, "CalcCoordsTransformFromOldCSToNew finishes working");
    return MCoordsInNewCS;
    //
    //    Matrix M6, M7 = null;
    //    M1.SetSize(3, 1, 1);
    //    M3.SetSize(3, 1, 1);
    //    M2.SetSize(3, 3, 1);
    //    //M7 = M2.TransposeTo();
    //    M7 = M2.InverseMatrix();
    //    M6 = M1 - M3;
    //    M4 = M7 * M6;
    //
}

Matrix_V2S Matrix_V2S::CalcCoordsTransformFromNewCSToOld(Matrix_V2S*CoordsInNewCSParam, Matrix_V2S*EulerAnglesParam, Matrix_V2S*CSOriginOfNewCSInOldParam, TValsShowHide*vsh){
    Matrix_V2S MCoordsInNewCS(3, 1), MEulerAngles (3, 1), MCSOriginOfNewCSInOld(3, 1);
    Matrix_V2S MDirCoss(3, 3), MCoordsInOldCS(3, 1), M1(3, 3), M2(3, 3), M3(3, 3);//, M4(3, 3);
    if(CoordsInNewCSParam!=NULL){
        MCoordsInNewCS=(*(CoordsInNewCSParam));
    }
    if(EulerAnglesParam!=NULL){
        MEulerAngles=(*(EulerAnglesParam));
    }
    if(CSOriginOfNewCSInOldParam!=NULL){
        MCSOriginOfNewCSInOld=(*(CSOriginOfNewCSInOldParam));
    }
    // Matrix      CalcMatrixOfDirCossByEulersAngles(Matrix MEulerAngles);
    MDirCoss=this->CalcMatrixOfDirCossByEulersAngles(       MEulerAngles);
    M1=MDirCoss;
    M2=M1*MCoordsInNewCS;
    MCoordsInOldCS=M2+MCSOriginOfNewCSInOld;
    return MCoordsInOldCS;
    //
    //    Matrix M6 = null;
    //    M1.SetSize(3, 1, 1);
    //    M3.SetSize(3, 1, 1);
    //    M2.SetSize(3, 3, 1);
    //    M6 = M2 * M1;
    //    M4 = M6 + M3;
    //
}

*/

   // };
    //};
//};




//----------------------------------------------------------------------------
//class Matrix_VL{
//    std::vector<std::vector<double>>y;
//    public:
Matrix_V2L::Matrix_V2L(int QLins, int QCols, double val){
    //this-construct();
    //this->SetSize(QLins, QCols);
    this->Set(val, QLins, QCols);
}
Matrix_V2L::Matrix_V2L(double**y, int QLins, int QCols){
    //this->construct();
    this->Set(y, QLins, QCols);
}
Matrix_V2L::Matrix_V2L(std::vector<std::vector<double>>y){
    //this->construct();
    this->Set(y);
}
Matrix_V2L::Matrix_V2L(double*y, int QLines, int QColumns){
    //this->construct();
    this->Set(y, QLines, QColumns);
}
Matrix_V2L::Matrix_V2L(std::vector<double>y, int QColumns){
    //this->construct();
    this->Set(y, QColumns);
}
Matrix_V2L::Matrix_V2L(double x, double y, double z){
    //this->construct();
    this->Set(x, y, z);
}
//
//Matrix_V2L::Matrix_V2L(const Matrix_V2L&obj){
//    //this->construct();
//    this->Assign(obj);
//}
//
Matrix_V2L::~Matrix_V2L(){
    //this->SetNull();//possible but no need
}

void Matrix_V2L::SetNull(){
    int QLines=this->GetQLines();
    //for(int i=1; i<=QLines; i++){//possible but no need
    //    this->y[i-1].clear();
   // }
    this->y.clear();
}
void Matrix_V2L::construct(){ this->SetOne(); }//NOp;
//void Matrix_VB::SetNullIni(){ this->y=NULL; this->QLins=0; this->QCols=0; }

Matrix_Prototype* Matrix_V2L::CreateMatrix(){
     return new Matrix_V2L();
}
 Matrix_Prototype* Matrix_V2L::clone(){
    return this;
}
int Matrix_V2L::GetTypeN() const{ return MatrixTypeN_V2L; }

/*
//mab no not ob S nablb virt'l
void Matrix_V2L::Assign(const Matrix_V2L&obj){
    this->SetNull();
    std::vector<double>curRow;
    double val;
    int QLins=obj.GetQLines(), QCols=obj.GetQColumns();
    for(int i=1; i<=QLins; i++){
        curRow.clear();
        for(int j=1; j<=QCols; j++){
            val=obj.GetComponent(i, j);
            curRow.push_back(val);
        }
        this->y.push_back(curRow);
    }
}
Matrix_V2L Matrix_V2L::operator = (const Matrix_V2L&obj){
    this->Assign(obj);
    return*this;
}
*/
void Matrix_V2L::SetSize(int QLins, int QCols){
    Array2DSetSize(this->y, QLins, QCols);
}


int Matrix_V2L::GetQLines() const { return this->y.size(); }
int Matrix_V2L::GetQColumns() const { return this->y[1-1].size(); }
//
void Matrix_V2L::Set(double**y, int QLins, int QCols){
    Matrix_Prototype::Set(y, QLins, QCols);
    //this->QLins=QLins;
    //this->QCols=QCols;
}
void Matrix_V2L::Set(std::vector<std::vector<double>>y){
    //this->QLins=y.size(), this->QCols=y[1-1].size();
    Matrix_Prototype::Set(y);
}
void Matrix_V2L::Set(double*y, int QLines,  int QColumns){
    Matrix_Prototype::Set(y, QLines, QColumns);
    //this->QLins=QLins;
    //this->QCols=QCols;
}
void Matrix_V2L::Set(std::vector<double>y, int QColumns){
    Matrix_Prototype::Set(y, QColumns);
    //this->QCols=QColumns;
    //this->QLins=y.size()/QColumns;
}
void Matrix_V2L::Set(double x, double y, double z){
    Matrix_Prototype::Set(x, y, z);
    //this->QLins=3;
    //this->QCols=1;
}
void Matrix_V2L::Set(double y, int QLines, int QColumns){
    Matrix_Prototype::Set(y, QLines, QColumns);
    //this->QCols=QColumns;
    //this->QLins=QLines;
}
//
double Matrix_V2L::GetComponent(int LineN, int ColN) const{
    double y=0;
    int QLins=this->GetQLines(), QCols=this->GetQColumns();
    if(LineN>=1 && LineN<=QLins && ColN>=1 && ColN<=QCols){
        y=this->y[LineN-1][ColN-1];
    }
    return y;
}
void Matrix_V2L::SetComponent(double val, int LineN, int ColN){
    int QLins=this->GetQLines(), QCols=this->GetQColumns();
    if(LineN>=1 && LineN<=QLins && ColN>=1 && ColN<=QCols){
        this->y[LineN-1][ColN-1]=val;
    }
}


/*
Matrix_V2L Matrix_V2L::operator +(const Matrix_V2L&obj){
    double x;
    int QLins=this->GetQLines(), QCols=this->GetQColumns();
    if(QLins==obj.GetQLines() && QCols==obj.GetQColumns()){
        for(int i=1; i<=QLins; i++){
            for(int j=1; j<=QCols; j++){
                x=this->y[i-1][j-1]+obj.y[i-1][j-1];
                this->y[i-1][j-1]=x;
            }
        }
    }
    return *this;
}
Matrix_V2L Matrix_V2L::operator -(const Matrix_V2L&obj){
    double x;
    int QLins=this->GetQLines(), QCols=this->GetQColumns();
    if(QLins==obj.GetQLines() && QCols==obj.GetQColumns()){
        for(int i=1; i<=QLins; i++){
            for(int j=1; j<=QCols; j++){
                x=this->y[i-1][j-1]-obj.y[i-1][j-1];
                this->y[i-1][j-1]=x;
            }
        }
    }
    return *this;
}
Matrix_V2L Matrix_V2L::operator *(const Matrix_V2L&obj){
    //double x;
    //int QM;
    //int QLins=this->GetQLines(), QCols=this->GetQColumns();
    //if(QCols==obj.GetQLines()){
    //    QM=QCols;
    //    for(int i=1; i<=QLins; i++){
   //         for(int j=1; j<=QCols; j++){
   //             for(int k=1; k<=QM; k++){
    //                x=this->y[i-1][k-1]-obj.y[k-1][j-1];
    //                this->y[i-1][j-1]=x;
    //            }
    //        }
    //    }
   // }
    //return *this;
    double x;
    Matrix_V2L MT=obj;
    Matrix_V2L M1=(*(this));
    //Matrix *M2=&obj;//invalid conversion from const Matrix_PS* to Matrix_PS* [-fpermissive]
    //Matrix *M2=obj;//cannot convert const Matrix_PS to Matrix_PS* in initialization
    Matrix_V2L *M2;
    //M2=&obj;
    //M2=obj;//cannot convert const Matrix to Matrix* in assignment
    //M2=MT;//cannot convert const Matrix to Matrix* in assignment
    M2=&MT;
    //Matrix M3;
    int QL1, QL2, QL3, QC1, QC2, QC3, QM;
    QL1=this->GetQLines(), QL2=obj.GetQLines(), QC1=this->GetQColumns(), QC2=obj.GetQColumns();
    double x1, x2, x3, mp;
    if(QC1==QL2){
        QM=QC1;
        QM=QL2;
        QL3=QL1;
        QC3=QC2;
        this->SetSize(QL3, QC3);
        //M3.SetSize(QL3, QC3);
        for(int i=1; i<=QL3; i++){
            for(int j=1; j<=QC3; j++){
                x3=0;
                for(int k=1; k<=QM; k++){
                    x1=M1.GetComponent(i, k);
                    x2=M2->GetComponent(k, j);
                    mp=x1*x2;
                    x3+=mp;
                }
                //M3.SetComponent(x3, i, j);
                this->SetComponent(x3, i, j);
            }
        }
    }
    //(*(this))=M3;
    return *this;
}
Matrix_V2L Matrix_V2L::operator /(const Matrix_V2L&obj){
    Matrix_V2L MR, MI;
    MI = obj.InverseMatrix();
    MR = MI * (*(this));
    return MR;
}

Matrix_V2L Matrix_V2L::operator *(double x){
    double cur;
    int QLins=this->GetQLines(), QCols=this->GetQColumns();
    for(int i=1; i<=QLins; i++){
        for(int j=1; j<=QCols; j++){
           cur=this->y[i-1][j-1]*x;
           this->y[i-1][j-1]=cur;
        }
    }
    return *this;
}
*/

void Matrix_V2L::SetLine(int N, double*x, int Q, int whatN1, int QDfltValsBefore, bool DefaulsAreZerosNotOwnVals){
    double*defaultVal=NULL;  //defaultVal=new double; (*(defaultVal))=0; // no need let em be random, af D wi be 0
    int QLines=this->GetQLines(), QColumns=this->GetQColumns();//, Ladded;
    int LreqCalcd=0, preN1_=0, preN2_=0, ownN1_=0, ownN2=0, ForOwnN1_=0, ForOwnN2=0, postN1=0, postN2_=0, QDefaultValuesBefore=0;
    //std::vector<double> row1;
    //Array1DAssign(row1, x, Q);
    //std::vector<double> GetNumericSubarray_byLs(std::vector<double> x, int FromN=1, int QDefaultValuesBefore=0, int LreqGiven=0);
    std::vector<double> row2;
    //std::vector<double>* VP=NULL;
    if(N>=1 && N<=QLines){
        if(DefaulsAreZerosNotOwnVals){
            //template <typename T>std::vector<T> Array1DGetSubArray_byLs_PreMark_v1(T*x, int Lold, T DefaultValue, int Lreq=0, int whatFromN=1, int FirstDefaultValues=0){
            row2= Array1DGetSubArray_byLs_PreMark_v1(x, Q, 0.0, QColumns, whatN1, QDfltValsBefore);
        }else{
            //template <typename T>std::vector<T> Array1DGetSubArray_byLs_PreMark_v1(T*x, int Lold,        T*DefaultRow, int Lreq=0, int whatFromN=1, int FirstDefaultValues=0, int dfltRowL_ifNegativeThanETLold=-1){
            row2 = Array1DGetSubArray_byLs_PreMark_v1(                                 x,        Q, this->y[N-1].data(),   QColumns,        whatN1,          QDfltValsBefore);
        }
        Array2DSetExtRowN(this->y, N, row2);
    }
}
void Matrix_V2L::SetLine(int N, std::vector<double>x, int whatN1, int QDfltValsBefore, bool DefaulsAreZerosNotOwnVals){
    this->SetLine(N, x.data(), x.size(), whatN1, QDfltValsBefore, DefaulsAreZerosNotOwnVals);
}
void Matrix_V2L::SetColumn(int N, double*x, int Q, int whatN1, int QDfltValsBefore, bool DefaulsAreOwnValsNotZeros){
    //std::vector<double> row1;
    //Array1DAssign(row1, x, Q);
    std::vector<double> row2;// = DefaulsAreZerosNotOwnVals==true? GetNumericSubarray_byLs(row1, whatN1, QDfltValsBefore, this->GetQLines()); : GetNumericSubarray_byLs(row1, this->y[N-1], whatN1, QDfltValsBefore, this->GetQLines());
    std::vector<double> ColumnN=this->GetArrayOfColN(N);
    int QColumns=this->GetQColumns();
    if(N>=1 && N<=QColumns){
        if(DefaulsAreOwnValsNotZeros){
            //template <typename T>std::vector<T> Array1DGetSubArray_byLs_PreMark_v1(T*x, int Lold,   T*DefaultRow, int Lreq=0, int whatFromN=1, int FirstDefaultValues=0, int dfltRowL_ifNegativeThanETLold=-1){
            row2 = Array1DGetSubArray_byLs_PreMark_v1(                                 x,        Q, ColumnN.data(),   QColumns,          whatN1,        QDfltValsBefore);
        }else{
            //template <typename T>std::vector<T> Array1DGetSubArray_byLs_PreMark_v1(T*x, int Lold, T DefaultValue, int Lreq=0, int whatFromN=1, int FirstDefaultValues=0){
            row2= Array1DGetSubArray_byLs_PreMark_v1(                                  x,        Q,            0.0,   QColumns,          whatN1,          QDfltValsBefore);
        }
        //template<typename T>Array2DSetIneRow(std::vector<std::vector<T>>&X, int IneRowN, std::vector<T>rowParam, int whatFromN=1, int QDefaultValsBefore=0, T*DefaultValParam=NULL, bool ExistingInsteadOfDefault=true, bool ExistingOnlyForGTLminNotIgnore=false, bool LastPossIfNotRectAtPos0NotIgnore=false){
        Array2DSetIneRow(this->y, N, row2);
    }
}
void Matrix_V2L::SetColumn(int N, std::vector<double>x, int whatN1, int QDfltValsBefore, bool DefaulsAreOwnValsNotZeros){
    this->SetColumn(N, x.data(), x.size(),  whatN1, QDfltValsBefore, DefaulsAreOwnValsNotZeros);
}
//
void Matrix_V2L::AddLine(double*x, int L, int whatN1, int QDefaultValsBefore){
    int QLines=this->GetQLines(), QColumns=this->GetQColumns();
    if((QLines==1 && QColumns==1 && L>1) || ((QLines==0 || QColumns==0)&&L>0)){
        this->SetNull();
    }
    std::vector<double>row0, row;
    Array1DAssign(row0, x, L);
    //row=GetNumericSubarray_byLs(row0, whatN1, QDefaultValsBefore, this->GetQColumns());
    //template <typename T>std::vector<T> Array1DGetSubArray_byLs_PreMark_v1(T*x, int Lold, T DefaultValue, int Lreq=0, int whatFromN=1, int FirstDefaultValues=0){
    row= Array1DGetSubArray_byLs_PreMark_v1(                                   x,        L,             0.0,   QColumns,          whatN1,       QDefaultValsBefore);
    Array2DAddExtRow(this->y, row);
}
void Matrix_V2L::AddLine(std::vector<double>x, int whatN1, int QDefaultValsBefore){
    this->AddLine(x.data(), x.size(), whatN1, QDefaultValsBefore);
}
void Matrix_V2L::AddColumn(double*x, int L, int whatN1, int QDefaultValsBefore){
    int QLines=this->GetQLines(), QColumns=this->GetQColumns();
    if((QLines==1 && QColumns==1 && L>1) || ((QLines==0 || QColumns==0)&&L>0)){
        this->SetNull();
    }
    std::vector<double>row0, row;
    Array1DAssign(row0, x, L);
    row=GetNumericSubarray_byLs(row0, whatN1, QDefaultValsBefore, this->GetQLines());
    Array2DAddIneRow(this->y, row);
}
void Matrix_V2L::AddColumn(std::vector<double>x, int whatN1, int QDefaultValsBefore){
    this->AddColumn(x.data(), x.size(), whatN1, QDefaultValsBefore);
}

void Matrix_V2L::InsLine(int N, double*x, int Q, int whatN1, int QDefaultValsBefore){
    int QLines=this->GetQLines(), QColumns=this->GetQColumns();
    std::vector<double> row;//, row0;
    //Array1DAssign(row0, x, Q);
    //row=GetNumericSubarray_byLs(row0, whatN1, QDefaultValsBefore, this->GetQColumns());
    //template <typename T>std::vector<T> Array1DGetSubArray_byLs_PreMark_v1(T*x, int Lold, T DefaultValue, int Lreq=0, int whatFromN=1, int FirstDefaultValues=0){
    row= Array1DGetSubArray_byLs_PreMark_v1(                                   x,        Q,            0.0,   QColumns,          whatN1,       QDefaultValsBefore);
    Array2DInsExtRow(this->y, N, row);
}
void Matrix_V2L::InsLine(int N, std::vector<double>x, int whatN1, int QDefaultValsBefore){
    this->InsLine(N, x.data(), x.size(), whatN1, QDefaultValsBefore);
}
void Matrix_V2L::InsColumn(int N, double*x, int Q, int whatN1, int QDefaultValsBefore){
    int QLines=this->GetQLines(), QColumns=this->GetQColumns();
    std::vector<double> row;//, row0;
    //Array1DAssign(row0, x, Q);
    //row=GetNumericSubarray_byLs(row0, whatN1, QDefaultValsBefore, this->GetQLines());
    //template <typename T>std::vector<T> Array1DGetSubArray_byLs_PreMark_v1(T*x, int Lold, T DefaultValue, int Lreq=0, int whatFromN=1, int FirstDefaultValues=0){
    row= Array1DGetSubArray_byLs_PreMark_v1(                                   x,        Q,            0.0,     QLines,          whatN1,       QDefaultValsBefore);
    Array2DInsIneRow(this->y, N, row);
}
void Matrix_V2L::InsColumn(int N, std::vector<double>x, int whatN1, int QDefaultValsBefore){
    this->InsColumn(N, x.data(), x.size(), whatN1,  QDefaultValsBefore);
}
void Matrix_V2L::DelLine(int N){
   Array2DDelExtRow(this->y, N);
}
void Matrix_V2L::DelColumn(int N){
    Array2DDelIneRow(this->y, N);
}
/*
void Matrix_V2L::Transpose(){
   Array2DTranspose(this->y);
}
Matrix_V2L Matrix_V2L::TransposeTo(){
    Matrix_V2L M=*this;
    M.Transpose();
    return M;
}
//
void Matrix_V2L::SetX(double val){
    if(this->GetQLines()>=1){
        this->SetComponent(val, 1, 1);
    }
}

void Matrix_V2L::SetY(double val){
    if(this->GetQLines()>=2){
        this->SetComponent(val, 2, 1);
    }
}
void Matrix_V2L::SetZ(double val){
    if(this->GetQLines()>=3){
        this->SetComponent(val, 3, 1);
    }
}
double Matrix_V2L::GetX(){
    return this->y[1-1][1-1];
}

double Matrix_V2L::GetY(){
   double y;
   if(this->GetQLines()>=2){
       y=this->y[2-1][1-1];
   }
   return y;
}

double Matrix_V2L::GetZ(){
    double y;
    if(this->GetQLines()>=3){
        y=this->y[3-1][1-1];
    }
    return y;
 }
//
Matrix_V2L Matrix_V2L::SubMatrix(int Line1N, int Col1N, int Line2N, int Col2N){
    Matrix_V2L M;
    int QLins=this->GetQLines(), QCols=this->GetQColumns(), QLinsNew=0, QColsNew=0, LineN, ColN;
    double val;
    if(Line1N>=1 && Line1N<=QLins && Line2N>=1 && Line2N<=QLins && Col1N>=1 && Col1N<=QCols && Col2N>=1 && Col2N<=QCols){
        if(Line1N<=Line2N  && Col1N<=Col2N){
            QLinsNew=Line2N-Line1N+1;
            QColsNew=Col2N-Col1N+1;
            M.SetSize(QLinsNew, QColsNew);
            //for(int i=Line1N; i<=Line2N; i++){
            for(int i=1; i<=QLinsNew; i++){
                //LineN=i-Line1N+1;
                LineN=i+Line1N-1;
                for(int j=Col1N; j<=Col2N; j++){
                    ColN=Col1N+j-1;
                    //ColN=j-Col1N+1;
                    val=this->GetComponent(LineN, ColN);
                    //val=this->GetElement(i, j);
                    M.SetComponent(val, i, j);
                    //M.SetElement(val, LineN, ColN);
                }
            }
        }else if(Line1N<=Line2N  && Col1N>Col2N){
            QLinsNew=Line2N-Line1N+1;
            QColsNew=Col1N-Col2N+1;
            M.SetSize(QLinsNew, QColsNew);
            for(int i=1; i<=QLinsNew; i++){
                LineN=i+Line1N-1;
                for(int j=1; j>=QColsNew; j++){
                    ColN=Col2N+j-1;//(changes to Col1N)
                    val=this->GetComponent(LineN, ColN);
                    M.SetComponent(val, i, j);
                }
            }
        }else if(Line1N>=Line2N  && Col1N<=Col2N){
            QLinsNew=Line1N-Line2N+1;
            QColsNew=Col2N-Col1N+1;
            M.SetSize(QLinsNew, QColsNew);
            for(int i=1; i<=QLinsNew; i++){
                LineN=Line2N+i-1;//(changes to Line1N)
                for(int j=1; j>=QColsNew; j++){
                    ColN=Col1N+j-1;
                    val=this->GetComponent(LineN, ColN);
                    M.SetComponent(val, i, j);
                }
            }
        }else{//both 1 > 2
            QLinsNew=Line1N-Line2N+1;
            QColsNew=Col1N-Col2N+1;
            M.SetSize(QLinsNew, QColsNew);
            for(int i=1; i<=QLinsNew; i++){
                LineN=Line2N+i-1;//(changes to Line1N)
                for(int j=1; j>=QColsNew; j++){
                    ColN=Col2N+j-1;//(changes to Col1N)
                    val=this->GetComponent(LineN, ColN);
                    M.SetComponent(val, i, j);
                }
            }
        }
    }
    return M;
}


Matrix_V2L Matrix_V2L::MinorTo(int LineN, int ColN) const
       {
           Matrix_V2L M;
           M=*this;
           // M = this.DelLineTo(LineN);
           // M = M.DelColTo(ColN);
           M.DelLine(LineN);
           M.DelColumn(ColN);
           return M;
       }

       void Matrix_V2L::Minor(int LineN, int ColumnN)
       {
           this->DelColumn(ColumnN);
           this->DelLine(LineN);
       }
       double Matrix_V2L::AlgSuppl(int LineN, int ColumnN) const
       {
           double d, r;
           int QLines=this->GetQLines(), QColumns=this->GetQColumns();
           Matrix_V2L M;// in |C# wa M = new Matrix_VB();
           if ((LineN < 1) || (LineN > QLines) || (ColumnN < 1) || (ColumnN > QColumns))
           {
               r = 666;
           }
           else
           {
               M = *this;
               M = MinorTo(LineN, ColumnN);
               d = M.Determinant();
               r = d;
               if ((LineN + ColumnN)%2 != 0)
               {
                   r *= (-1);
               }
           }
           return r;
       }
        double Matrix_V2L::Determinant() const
       {
           double r, d, AlgSpl;
           int RowN;
           int QColumns=this->GetQColumns(), QLines=this->GetQLines();
           if (QColumns != QLines) return 666;
           else
           {
               if (QLines == 1)
               {
                   r = this->y[1 - 1][1 - 1];
               }
               else if (QLines == 2)
               {
                   r = this->y[1 - 1][1 - 1] * this->y[2-1][2-1] - this->y[1 - 1][2 - 1] * this->y[2 - 1][1 - 1];
               }
               else if (QLines == 3)
               {
                   r = this->y[1-1][1-1] * this->y[2-1][2-1] * this->y[3-1][3-1] + this->y[1-1][3-1] * this->y[2-1][1-1] * this->y[3-1][2-1] + this->y[3-1][1-1] * this->y[1-1][2-1] * this->y[2-1][3-1] -
                       this->y[1-1][3-1] * this->y[2-1][2-1] * this->y[3-1][1-1] - this->y[1-1][1-1] * this->y[2-1][3-1] * this->y[3-1][2-1] - this->y[3-1][3-1] * this->y[1-1][2-1] * this->y[2-1][1-1];
               }
               else
               {
                   RowN = 1;
                   r = 0;
                   for (int j = 1; j <= QColumns; j++)
                   {
                       //r += GetComponent(RowN, j) * AlgSuppl(RowN, j);
                       r += GetComponent(RowN, j) * AlgSuppl(RowN, j);
                   }
               }
           }
           return r;
       }

       Matrix_V2L Matrix_V2L::UnionMatrix() const
       {
           int QLines = this->GetQLines(), QColumns = this->GetQColumns();
           int TQLines = QColumns, TQColumns = QLines;
           //Matrix_VB MT = new Matrix_VB();
           //MT = this.TransposeTo();
           Matrix_V2L MUR;//in C# = new Matrix_VB();
           //double[][] c;
           double CurAlgSuppl;//, det = this.Determinant();
          // c = new double[TQLines][];
           //for (int i = 1; i <= TQLines; i++)
           //{
           //    c[i - 1] = new double[TQColumns];
           //}
           MUR.SetSize(TQLines, TQColumns);
           //
           for (int i = 1; i <= QLines; i++)
           {
               for (int j = 1; j <= QColumns; j++)
               {
                   CurAlgSuppl = this->AlgSuppl(i, j);
                   //c[j - 1][i - 1] = CurAlgSuppl;
                   MUR.SetComponent(CurAlgSuppl, j, i);
               }
           }
           //
           //MUR = new Matrix_VB(c, TQLines, TQColumns);
           return MUR;
       }
       Matrix_V2L Matrix_V2L::InverseMatrix() const
       {
           Matrix_V2L MIR;
           double det = this->Determinant();
           MIR = this->UnionMatrix();
           if (det != 0) MIR = MIR * (1 / det);
           return MIR;
       }
*/


/*QString Matrix_VL::ToString(){
    QString str="", cur;
    str="[[";
    int QLins=this->GetQLines(), QCols=this->GetQColumns();
    for(int i=1; i<=QLins-1; i++){
        str+="[";
        for(int j=1; j<=QCols-1; j++){
            cur.setNum(this->y[i-1][j-1]);
            str+=cur;
            str+=", ";
            cur="";
        }
        cur.setNum(this->y[i-1][QCols-1]);
        str+=cur;
        str+="];";
    }
    str+="[";
    //LastLine
    for(int j=1; j<=QCols-1; j++){
        cur.setNum(this->y[QLins-1][j-1]);
        str+=cur;
        str+=", ";
        cur="";
    }
    //LastElt
    cur.setNum(this->y[QLins-1-1][QCols-1]);
    str+=cur;
    str+="]]";
    return str;
}
std::string Matrix_VL::ToStdString(){
    std::string ss;
    QString qs;
    qs=this->ToString();
    ss=qs.toStdString();
    return ss;
}
*/
//Matrix_VB GivensRotationTransformMatrix(int N1, int N2, Matrix_VB*M=NULL, TValsShowHide*vsh=NULL);
/*
Matrix_V2L Matrix_V2L::GivensRotationTransformMatrix(int N1, int N2, Matrix_V2L*M, TValsShowHide*vsh){
    int QL;
    Matrix_V2L T(QL,QL);
    if(M==NULL)M=this;
    QL=M->GetQLines();
    double a11=M->GetComponent(N1, N1), a21=M->GetComponent(N2, N1), c12, s12, denom;
    for(int i=1; i<=QL; i++){
        for(int j=1; j<=QL; j++){
            T.SetComponent(0, i, j);
        }
    }
    for(int i=1; i<=QL; i++){
        T.SetComponent(1, i, i);
    }
    if(a21==0){
        c12=1;
        s12=0;
    }else{
        denom=sqrt(a11*a11+a21*a21);
        c12=a11/denom;
        s12=a21/denom;
    }
    writeln(vsh, " T=" +M->ToString());
    writeln(vsh, " a11=" +FloatToStr(a11) + " a21=" +FloatToStr(a21)+ " c12=" + FloatToStr(c12) + " s12=" + FloatToStr(s12));
    T.SetComponent(c12, N1, N1);
    T.SetComponent(c12, N2, N2);
    T.SetComponent(s12, N1, N2);
    T.SetComponent(-s12, N2, N1);
    return T;
}
//void QRDecomposition(Matrix_VB&Q, Matrix_VB&R,  Matrix_VB*M=NULL, TValsShowHide* vsh = NULL);
void Matrix_V2L::QRDecomposition(Matrix_V2L&Q, Matrix_V2L&R,  Matrix_V2L*M, TValsShowHide* vsh)
{
    writeln(vsh, "QRDecomposition starts working");
    if (M == NULL) M = this;
    int QL = M->GetQLines();
    std::vector<int>ns;// = new int[QL-1];
    //int countTsInLine;
    int jn, ii, jj;
    Q.SetSize(QL, QL);
    for (int i = 1; i <= QL; i++)
    {
        Q.SetComponent(1, i, i);
    }
    std::vector<std::vector<Matrix_V2L>>Ts;
    std::vector<Matrix_V2L>T1;
    Matrix_V2L M1;
    for (int i = 1; i <= QL - 1; i++)
    {
        ns.push_back(QL - 1 - i + 1);
        T1.clear();
        for(int j=1; j<=ns[i - 1]; j++){
            T1.push_back(M1);
        }
        Ts.push_back(T1);
    }
    std::vector<Matrix_V2L>TL;
    Matrix_V2L T;//, R=new Matrix(QL, QL);
    R.Assign(*M);
    writeln(vsh, "M="+M->ToString());
    writeln(vsh, "ini Q=" + Q.ToString());
    writeln(vsh, "ini R=" + R.ToString());
    for (int i = 1; i <= QL - 1; i++)
    {
        //countTsInLine = 0;
        for (int j = i + 1; j <= QL; j++)
        {
            //countTsInLine++; //in AddToVector
            writeln(vsh, "i=" + IntToStr(i)+" j="+IntToStr(j));
            T = R.GivensRotationTransformMatrix(i, j, NULL, vsh);
            writeln(vsh, "T=" + T.ToString());
            R = T * R;
            writeln(vsh, "cur R=" + R.ToString());
            jn = j - i;
            Ts[i - 1][jn - 1] = T;
        }
    }
    //R is completed, beginning wit Q
    for(int i=1; i<=QL-1; i++){
        ii = QL - 1 - i + 1;
        for (int j = 1; j <= ns[ii - 1]; j++)
        {
            jj = ns[ii - 1] - j + 1;
            T = Ts[ii - 1][jj - 1];
            Q = Q * T;
        }
    }
    Q.Transpose();
    //
    writeln(vsh, "Answer:");
    writeln(vsh,"Q="+Q.ToString()+" R="+R.ToString());
    writeln(vsh, "QRDecomposition finishes working");

}//fn QR decomp

Matrix_V2L Matrix_V2L::CalcMatrixOfDirCossByEulersAngles (Matrix_V2L MEulerAngles){
    Matrix_V2L MR(3, 3);
    double gamma = MEulerAngles.GetX() * M_PI / 180,
                       psi = MEulerAngles.GetY() * M_PI / 180,
                       theta = MEulerAngles.GetZ() * M_PI / 180;
                MR.SetComponent(cos(psi) * cos(theta), 1, 1);
                MR.SetComponent(sin(psi) * sin(gamma) - cos(psi) * sin(theta) * cos(gamma), 1, 2);
                MR.SetComponent(sin(psi) * cos(gamma) + cos(psi) * sin(theta) * sin(gamma), 1, 3);
                MR.SetComponent(sin(theta), 2, 1);
                MR.SetComponent(cos(theta) * cos(gamma), 2, 2);
                MR.SetComponent(-cos(theta) * sin(gamma), 2, 3);
                MR.SetComponent(-sin(psi) * cos(theta), 3, 1);
                MR.SetComponent(cos(psi) * sin(gamma) + sin(psi) * sin(theta) * cos(gamma), 3, 2);
                MR.SetComponent(cos(psi) * cos(gamma) - sin(psi) * sin(theta) * sin(gamma), 3, 3);
                return MR;
}

Matrix_V2L Matrix_V2L::CalcCoordsTransformFromOldCSToNew(Matrix_V2L*CoordsInOldCSParam, Matrix_V2L*EulerAnglesParam, Matrix_V2L*CSOriginOfNewCSInOldParam, TValsShowHide*vsh){
    Matrix_V2L MCoordsInOldCS(3, 1), MEulerAngles (3, 1),MCSOriginOfNewCSInOld(3, 1);
    Matrix_V2L MDirCoss(3, 3), MCoordsInNewCS(3, 1), M1(3, 3), M2(3, 3), M3(3, 3);//, M4(3, 3);
    writeln(vsh, "CalcCoordsTransformFromOldCSToNew starts working");
    if(CoordsInOldCSParam!=NULL){
        MCoordsInOldCS=(*(CoordsInOldCSParam));
    }
    writeln(vsh, "MCoordsInOldCS");
    MCoordsInOldCS.StdCOut2D();
    if(EulerAnglesParam!=NULL){
        MEulerAngles=(*(EulerAnglesParam));
    }
    writeln(vsh, "MEulerAngles");
    MCoordsInOldCS.StdCOut2D();
    if(CSOriginOfNewCSInOldParam!=NULL){
        MCSOriginOfNewCSInOld=(*(CSOriginOfNewCSInOldParam));
    }
    writeln(vsh, "MCSOriginOfNewCSInOld");
    MCSOriginOfNewCSInOld.StdCOut2D();
    MDirCoss=this->CalcMatrixOfDirCossByEulersAngles(MEulerAngles);
    writeln(vsh, "MDirCoss");
    MDirCoss.StdCOut2D();
    M1=MDirCoss.InverseMatrix();
    writeln(vsh, "MDirCoss Inverted");
    M1.StdCOut2D();
    M2=MCoordsInOldCS-MCSOriginOfNewCSInOld;
    writeln(vsh, "M2=MCoordsInOldCS-MCSOriginOfNewCSInOld");
    M2.StdCOut2D();
    writeln(vsh, "M3=M1*M2");
    M3=M1*M2;
    M3.StdCOut2D();
    MCoordsInNewCS=M3;
    writeln(vsh, "CalcCoordsTransformFromOldCSToNew finishes working");
    return MCoordsInNewCS;
    //
    //    Matrix M6, M7 = null;
    //    M1.SetSize(3, 1, 1);
    //    M3.SetSize(3, 1, 1);
    //    M2.SetSize(3, 3, 1);
    //    //M7 = M2.TransposeTo();
    //    M7 = M2.InverseMatrix();
    //    M6 = M1 - M3;
    //    M4 = M7 * M6;
    //
}

Matrix_V2L Matrix_V2L::CalcCoordsTransformFromNewCSToOld(Matrix_V2L*CoordsInNewCSParam, Matrix_V2L*EulerAnglesParam, Matrix_V2L*CSOriginOfNewCSInOldParam, TValsShowHide*vsh){
    Matrix_V2L MCoordsInNewCS(3, 1), MEulerAngles (3, 1), MCSOriginOfNewCSInOld(3, 1);
    Matrix_V2L MDirCoss(3, 3), MCoordsInOldCS(3, 1), M1(3, 3), M2(3, 3), M3(3, 3);//, M4(3, 3);
    if(CoordsInNewCSParam!=NULL){
        MCoordsInNewCS=(*(CoordsInNewCSParam));
    }
    if(EulerAnglesParam!=NULL){
        MEulerAngles=(*(EulerAnglesParam));
    }
    if(CSOriginOfNewCSInOldParam!=NULL){
        MCSOriginOfNewCSInOld=(*(CSOriginOfNewCSInOldParam));
    }
    // Matrix      CalcMatrixOfDirCossByEulersAngles(Matrix MEulerAngles);
    MDirCoss=this->CalcMatrixOfDirCossByEulersAngles(       MEulerAngles);
    M1=MDirCoss;
    M2=M1*MCoordsInNewCS;
    MCoordsInOldCS=M2+MCSOriginOfNewCSInOld;
    return MCoordsInOldCS;
    //
   //    Matrix M6 = null;
   //     M1.SetSize(3, 1, 1);
   //     M3.SetSize(3, 1, 1);
   //     M2.SetSize(3, 3, 1);
   //     M6 = M2 * M1;
  //      M4 = M6 + M3;
   //
}
*/


   // };
    //};
//};

// class Matrix_A2P{
//     double**y;
//     int QLins, QCols;
//     public:
Matrix_P2A::Matrix_P2A(int QLins, int QCols, double val){
    this->construct();
    this->SetSize(QLins, QCols);
    for(int i=1; i<=QLins; i++){
        for(int j=1; j<=QCols; j++){
            this->SetComponent(val, i, j);
        }
    }
}
Matrix_P2A::Matrix_P2A(double**y, int QLins, int QCols){
    this->construct();
    this->Set(y, QLins, QCols);
}
Matrix_P2A::Matrix_P2A(std::vector<std::vector<double>>y){
    this->construct();
    this->Set(y);
}

Matrix_P2A::Matrix_P2A(double*y, int QLins, int QCols){
    this->construct();
    this->Set(y, QLins, QCols);
}
Matrix_P2A::Matrix_P2A(std::vector<double>y, int QCols){
    this->construct();
    this->Set(y, QCols);
}
Matrix_P2A::Matrix_P2A(double x, double y, double z){
    this->construct();
    this->Set(x, y, z);
}
//Matrix_P2A::Matrix_P2A(const Matrix_P2A&obj){
//    this->construct();
//    this->Assign(obj);
//}
void Matrix_P2A::construct(){
    this->data.SetNullIni();
    this->SetOne(0);
}
void Matrix_P2A::SetNull(){
    this->data.SetNull();
}
//void Matrix_P2A::SetNullIni(){
//    this->data.SetNull();
//}
Matrix_P2A::~Matrix_P2A(){
    this->SetNull();
}
//void Matrix_P2A::Assign(const Matrix_P2A&obj){
//    this->data=obj.data;
//}
//Matrix_P2A& Matrix_P2A::operator = (const Matrix_P2A&obj){
//    this->Assign(obj);
//    return*this;
//}
void Matrix_P2A::SetSize(int QLins, int QCols){
    double*dflt=new double;
    (*(dflt))=0;
    this->data.SetSize(QLins, QCols, true, dflt);
    delete dflt;
}


int Matrix_P2A::GetQLines() const { return this->data.GetQExtRows(); }
int Matrix_P2A::GetQColumns() const { return this->data.GetMaxLength(); }
//
void Matrix_P2A::Set(double**y, int QLins, int QCols){
    Matrix_Prototype::Set(y, QLins, QCols);
}

void Matrix_P2A::Set(std::vector<std::vector<double>>y){
    Matrix_Prototype::Set(y);
}
void Matrix_P2A::Set(double*y, int QLines,  int QColumns){
    Matrix_Prototype::Set(y, QLines, QColumns);
}
void Matrix_P2A::Set(std::vector<double>y, int QColumns){
    Matrix_Prototype::Set(y, QColumns);
}
void Matrix_P2A::Set(double x, double y, double z){
    Matrix_Prototype::Set(x, y, z);
}
void Matrix_P2A::Set(double y, int QLines, int QColumns){
    Matrix_Prototype::Set(y, QLines, QColumns);
}
//
double Matrix_P2A::GetComponent(int LineN, int ColN) const{
    double y=0;
    int QLins=this->GetQLines(), QCols=this->GetQColumns();
    if(LineN>=1 && LineN<=QLins && ColN>=1 && ColN<=QCols){
        y=this->data.GetElement_AsVal(LineN, ColN);
    }
    return y;
}
void Matrix_P2A::SetComponent(double val, int LineN, int ColN){
    int QLins=this->GetQLines(), QCols=this->GetQColumns();
    if(LineN>=1 && LineN<=QLins && ColN>=1 && ColN<=QCols){
        this->data.SetElement(val, LineN, ColN);
    }
}

Matrix_Prototype* Matrix_P2A::CreateMatrix(){
    return new Matrix_P2A();
}

Matrix_Prototype* Matrix_P2A::clone(){
    return this;
}
int Matrix_P2A::GetTypeN() const{
    return MatrixTypeN_P2A;//P/V, 1/2, S/L/A
}


/*
Matrix_P2A Matrix_P2A::VectorOfLineN(int LineN){
 Matrix_P2A M;
 bool ColOfVectorNotHorLine=true;
 std::vector<double>row;
 int QLins=this->GetQLines(), QCols=this->GetQColumns();
 if(LineN>=1 && LineN<=QLins){
     for(int i=1; i<=QCols; i++){
         row.push_back(this->GetComponent(LineN, i));
     }
 }
 M.Set(row, ColOfVectorNotHorLine);
 return M;
}
Matrix_P2A Matrix_P2A::MatrixLineN(int LineN){
 Matrix_P2A M;
 bool ColOfVectorNotHorLine=false;
 std::vector<double>row;
 int QLins=GetQLines(), QCols=GetQColumns();
 if(LineN>=1 && LineN<=QLins){
     for(int i=1; i<=QCols; i++){
         row.push_back(this->GetComponent(LineN, i));
     }
 }
 M.Set(row, ColOfVectorNotHorLine);
 return M;
}
Matrix_P2A Matrix_P2A::VectorOfColN(int ColN){
 Matrix_P2A M;
 bool ColOfVectorNotHorLine=true;
 std::vector<double>row;
 int QLins=GetQLines(), QCols=GetQColumns();
 if(ColN>=1 && ColN<=QCols){
     for(int i=1; i<=QLins; i++){
         row.push_back(this->GetComponent(i, ColN));
     }
 }
 M.Set(row, ColOfVectorNotHorLine);
 return M;
}
std::vector<double> Matrix_P2A::GetLineN(int LineN){
 std::vector<double>R;
 int QLins=GetQLines(), QCols=GetQColumns();
 if(LineN>=1 && LineN<=QLins){
     for(int i=1; i<=QCols; i++){
         R.push_back(this->GetComponent(LineN, i));
     }
 }
 return R;
}
std::vector<double>Matrix_P2A::GetColN(int ColN){
    std::vector<double>R;
    int QLins=GetQLines(), QCols=GetQColumns();
    if(ColN>=1 && ColN<=QCols){
    for(int i=1; i<=QLins; i++){
        R.push_back(this->GetComponent(i, ColN));
    }
 }
 return R;
}
Matrix_P2A Matrix_P2A::operator +(const Matrix_P2A&obj){
    double x;
    int QLins=GetQLines(), QCols=GetQColumns();
    if(this->GetQLines()==obj.GetQLines() && this->GetQColumns()==obj.GetQColumns()){
        for(int i=1; i<=QLins; i++){
            for(int j=1; j<=QCols; j++){
                x=this->GetComponent(i, j)+obj.GetComponent(i, j);
                 this->SetComponent(x, i, j);
            }
        }
    }
    return *this;
}
Matrix_P2A Matrix_P2A::operator -(const Matrix_P2A&obj){
    double x;
    int QLins=GetQLines(), QCols=GetQColumns();
    if(this->GetQLines()==obj.GetQLines() && this->GetQColumns()==obj.GetQColumns()){
        for(int i=1; i<=QLins; i++){
            for(int j=1; j<=QCols; j++){
                x=this->GetComponent(i, j)-obj.GetComponent(i, j);
                this->SetComponent(x, i, j);
            }
        }
    }
    return *this;
}
Matrix_P2A Matrix_P2A::operator *(const Matrix_P2A&obj){
    //double x;
    // int QM, QLins=GetQLines(), QCols=GetQColumns();
    // if(this->GetQColumns()==obj.GetQLines()){
    //    QM=this->GetQColumns();
    //    for(int i=1; i<=QLins; i++){
    //        for(int j=1; j<=QCols; j++){
    //            for(int k=1; k<=QM; k++){
    //                x=this->GetComponent(i, k)-obj.GetComponent(k, j);
    //                this->SetComponent(x, i, j);
    //            }
    //        }
    //    }
    //}
    //return *this;
    double x;
    Matrix_P2A MT=obj;
    Matrix_P2A M1=(*(this));
    //Matrix *M2=&obj;//invalid conversion from const Matrix_PS* to Matrix_PS* [-fpermissive]
    //Matrix *M2=obj;//cannot convert const Matrix_PS to Matrix_PS* in initialization
    Matrix_P2A *M2;
    //M2=&obj;
    //M2=obj;//cannot convert const Matrix to Matrix* in assignment
    //M2=MT;//cannot convert const Matrix to Matrix* in assignment
    M2=&MT;
    //Matrix M3;
    int QL1, QL2, QL3, QC1, QC2, QC3, QM;
    QL1=this->GetQLines(), QL2=obj.GetQLines(), QC1=this->GetQColumns(), QC2=obj.GetQColumns();
    double x1, x2, x3, mp;
    if(QC1==QL2){
        QM=QC1;
        QM=QL2;
        QL3=QL1;
        QC3=QC2;
        this->SetSize(QL3, QC3);
        //M3.SetSize(QL3, QC3);
        for(int i=1; i<=QL3; i++){
            for(int j=1; j<=QC3; j++){
                x3=0;
                for(int k=1; k<=QM; k++){
                    x1=M1.GetComponent(i, k);
                    x2=M2->GetComponent(k, j);
                    mp=x1*x2;
                    x3+=mp;
                }
                //M3.SetComponent(x3, i, j);
                this->SetComponent(x3, i, j);
            }
        }
    }
    //(*(this))=M3;
    return *this;
}
Matrix_P2A Matrix_P2A::operator /(const Matrix_P2A&obj){
    Matrix_P2A MR, MI;
    MI = obj.InverseMatrix();
    MR = MI * (*(this));
    return MR;
}

Matrix_P2A Matrix_P2A::operator *(double x){
    double cur;
    int QLins=GetQLines(), QCols=GetQColumns();
    for(int i=1; i<=QLins; i++){
        for(int j=1; j<=QCols; j++){
            cur=this->GetComponent(i, j)*x;
            this->SetComponent(cur, i, j);
        }
    }
    return *this;
}


void Matrix_P2A::Transpose(){
    this->data.Transpose();
}
Matrix_P2A Matrix_P2A::TransposeTo(){
    Matrix_P2A M=*this;
    M.Transpose();
    return M;
}

void Matrix_P2A::SetX(double val){
    if(this->GetQLines()>=1){
        this->SetComponent(val, 1, 1);
    }
}

void Matrix_P2A::SetY(double val){
    if(this->GetQLines()>=2){
        this->SetComponent(val, 2, 1);
    }
}
void Matrix_P2A::SetZ(double val){
    if(this->GetQLines()>=3){
        this->SetComponent(val, 3, 1);
    }
}
double Matrix_P2A::GetX(){
    return this->GetComponent(1, 1);
}

double Matrix_P2A::GetY(){
   double y;
   if(this->GetQLines()>=2){
       y=this->GetComponent(2, 1);
   }
   return y;
}

double Matrix_P2A::GetZ(){
    double y;
    if(this->GetQLines()>=3){
        y=this->GetComponent(3, 1);
    }
    return y;
}
*/


/*void Matrix::AddLine(std::vector<double>x, int N1){
 int L, Lmin, whereN, whatN;
 this->SetSize(this->QLins+1, this->QCols);
 L=x.size();
 //Lmin = L <= this->QCols ? L : this->QCols;//uz N1=1
 Lmin = L <= this->QCols-N1+1 ? L : this->QCols-N1+1;
 for(int i=1; i<=N1-1; i++){
     this->y[this->QLins-1][i-1]=0;
 }
 for(int i=1; i<=Lmin; i++){
     this->y[this->QLins-1][N1-1+i-1]=x[i-1];
 }
 for(int i=Lmin+1; i<=this->QCols; i++){
     this->y[this->QLins-1][i-1]=0;
 }
 //return *this;
}
void Matrix::AddColumn(std::vector<double>x, int N1){
 int L, Lmin, whereN, whatN;
 this->SetSize(this->QLins, this->QCols+1);
 L=x.size();
 //Lmin = L <= this->QCols ? L : this->QCols;//uz N1=1
 Lmin = L <= this->QLins-N1+1 ? L : this->QLins-N1+1;
 for(int i=1; i<=N1-1; i++){
     this->y[i-1][this->QCols-1]=0;
 }
 for(int i=1; i<=Lmin; i++){
     this->y[N1-1+i-1][this->QCols-1]=x[i-1];
 }
 for(int i=Lmin+1; i<=this->QLins; i++){
     this->y[i-1][this->QCols-1]=0;
 }
 //return *this;
}
void Matrix::InsLine(int N, std::vector<double>x, int N1){
 int L, Lmin;//, whereN, whatN;
 if(N>=1 && N<=this->QLins){
     this->SetSize(this->QLins+1, this->QCols);
     for(int i=QLins-1; i>=N; i--){
         for(int j=1; j<=this->QCols; j++){
             y[i+1-1][j-1]= y[i-1][j-1];
         }
     }
     L=x.size();
     //Lmin = L <= this->QCols ? L : this->QCols;//uz N1=1
     Lmin = L <= this->QCols-N1+1 ? L : this->QCols-N1+1;
     for(int i=1; i<=N1-1; i++){
         this->y[N-1][i-1]=0;
     }
     for(int i=1; i<=Lmin; i++){
         this->y[N-1][N1-1+i-1]=x[i-1];
     }
     for(int i=Lmin+1; i<=this->QCols; i++){
         this->y[N-1][i-1]=0;
     }
 }
 //return *this;
}
void Matrix::InsColumn(int N, std::vector<double>x, int N1){
 int L, Lmin;//, whereN, whatN;
 if(N>=1 && N<=this->QCols){
     this->SetSize(this->QLins, this->QCols+1);
     for(int i=QCols-1; i>=N; i--){
         for(int j=1; j<=this->QLins; j++){
             y[j-1][i+1-1]= y[j-1][i-1];
         }
     }
     L=x.size();
     //Lmin = L <= this->QCols ? L : this->QCols;//uz N1=1
     Lmin = L <= this->QCols-N1+1 ? L : this->QCols-N1+1;
     for(int i=1; i<=N1-1; i++){
         this->y[N-1][i-1]=0;
     }
     for(int i=1; i<=Lmin; i++){
         this->y[N-1][N1-1+i-1]=x[i-1];
     }
     for(int i=Lmin+1; i<=this->QCols; i++){
         this->y[N-1][i-1]=0;
     }
 }
 //return *this;
}*/
//
void Matrix_P2A::SetLine(int N, double*x, int L, int FromN, int QDefaultBefore, bool DefaultValsAreZerosNotOwn){
    //std::vector<double>row;
    int QL=this->GetQLines(), QC=this->GetQColumns();
    Array2DSize size(QL, QC);
    //int preN1_, preN2_, ownN1_, ownN2, ForOwnN1_, ForOwnN2, postN1, postN2_;
    //int LreqCalcd;
    //void Calc_Array1DSubArray_byLs_MarkupNs_vFmls(int&LreqCalcd, int&preN1_, int&preN2_, int&ownN1_, int&ownN2, int&ForOwnN1_, int&ForOwnN2, int&postN1, int&postN2_, int Lini, int LreqGiven, int whatFromN, int QDefaultValuesBefore){
    //Calc_Array1DSubArray_byLs_MarkupNs_vFmls(LreqCalcd, preN1_, preN2_, ownN1_, ownN2, ForOwnN1_, ForOwnN2, postN1, postN2_, QC, 0, whatFromN, QDefaultValuesBefore);
    double *dflt=new double;
    (*(dflt))=0;
    //template<typename T> void Array2DSetExtRowN(T**&X, Array2DSize&size, int ExtRowN, int wholeRowL=0, T*rowParam=NULL, int whatFromN=1, int QDefaultValuesBefore=0, bool KeepIfRect=true, T*DefaultValParam=NULL, bool DefaultIsSpecNotOwn=false){
    this->data.SetExtRow(                          N,               L,               x,           FromN,             QDefaultBefore,                false,                   dflt,       DefaultValsAreZerosNotOwn);
    delete dflt;
}
void  Matrix_P2A::SetLine(int N, std::vector<double>x, int FromN, int QDfltValsBefore, bool DefaultValsAreOwnNotZeros){
 //template<typename T> void Array2DSetExtRowN(std::vector<std::vector<T>>X, int ExtRowN, std::vector<T>rowParam, T*DefaultValParam=NULL, int FromN=1, int DefaultValuesBefore=0, bool KeepIfRect=false)
 //template<typename T> void Array2DSetExtRowN(T**&X, const Array2DSize&size, int ExtRowN, std::vector<T>rowParam, T*DefaultValParam=NULL, int FromN=1, int DefaultValuesBefore=0, bool KeepIfRect=false)
 //template<typename T> void Array2DSetExtRowN(T**&X, const Array2DSize&size, int ExtRowN, int wholeRowL=0, T*rowParam=NULL, T*DefaultValParam=NULL, int FromN=1, int DefaultValuesBefore=0, bool KeepIfRect=false)
 Array2DSize size(this->GetQLines(), this->GetQColumns());
 double*defaultVal=new double;
 *defaultVal=0;
 //Array2DSetExtRowN(T**&X,  const Array2DSize&size, int ExtRowN, std::vector<T>rowParam, T*DefaultValParam=NULL, int FromN=1, int DefaultValuesBefore=0, bool KeepIfRect=false)
 this->data.SetExtRow(N,                      x,             defaultVal,      FromN,           QDfltValsBefore, true);
 delete defaultVal;
}
void  Matrix_P2A::SetColumn(int N, double*x, int L, int FromN, int QDefaultBefore, bool DefaultValsAreOwnNotZeros){
 int QL=this->GetQLines(), QC=this->GetQColumns();
 Array2DSize size(QL,  QC);
 double*dflt=new double;
 (*(dflt))=0;
 //template<typename T>Array2DSetIneRow(T**&X, Array2DSize& size, int IneRowN, T*rowParam=NULL, int wholeRowL=0, int whatFromN=1, int QDefaultValsBefore=0, T*DefaultValParam=NULL, bool ExistingInsteadOfDefault=true, bool ExistingOnlyForGTLminNotIgnore=false, bool LastPossIfNotRectAtPos0NotIgnore=false){
 this->data.SetIneRow(           N,               x,               L,           FromN,       QDefaultBefore,             dflt,                    DefaultValsAreOwnNotZeros);//            ExistingOnlyForGTLminNotIgnore,      LastPossIfNotRectAtPos0NotIgnore);
 delete dflt;
}
void  Matrix_P2A::SetColumn(int N, std::vector<double>x, int whatN1, int QDfltValsBefore, bool DefaultValsAreOwnNotZeros){
 //template<typename T> void Array2DSetIneRowN(T**&X, Array2DSize&size, int IneRowN, int wholeRowL=0, T*rowParam=NULL, T*DefaultValParam=NULL, int FromN=1, int DefaultValuesBefore=0, int IfNGTLmin_Ignore0DoNil1Add2Srtetch3=1)
 //template<typename T> void Array2DSetIneRowN(T**&X, Array2DSize&size, int IneRowN, std::vector<T>rowParam, T*DefaultValExt=NULL, int FromN=1, int DefaultValuesBefore=0, int IfNGTLmin_Ignore0DoNil1Add2Srtetch3=1)
 //template<typename T> void Array2DSetIneRowN(std::vector<std::vector<T>>&X, std::vector<T>rowParam, int IneRowN, T*DefaultValExt=NULL, int FromN=1, int DefaultValuesBefore=0, int IfNGTLmin_Ignore0DoNil1Add2Srtetch3=1)
 Array2DSize size(this->GetQLines(), this->GetQColumns());
 double*defaultVal=new double;
 *defaultVal=0;
 //template<typename T>Array2DSetIneRow(T**&X, Array2DSize& size, int IneRowN, std::vector<T>rowParam, int whatFromN=1, int QDefaultValsBefore=0, T*DefaultValParam=NULL, bool ExistingInsteadOfDefault=true, bool ExistingOnlyForGTLminNotIgnore=false, bool LastPossIfNotRectAtPos0NotIgnore=false){
 this->data.SetIneRow(           N,                      x, whatN1,                   QDfltValsBefore, defaultVal);
 delete defaultVal;
}
//
void Matrix_P2A::AddLine(double*x, int L, int whatN1, int QDefaultValsBefore){
 //template<typename T>void Array2DAddExtRow(T**&X, Array2DSize&size, T*rowParam, int wholeRowL=0, int whereFromN=1, int whatFromN=1, bool RectNotVar=false, T*DfltValParam=NULL, TValsShowHide*vsh=NULL)
 int QLines=this->GetQLines(), QColumns=this->GetQColumns();
 Array2DSize size(QLines, QColumns);
 double*defaultVal=new double;
 *defaultVal=0;
 if((QLines==1 && QColumns==1 && L>1) || ((QLines==0 || QColumns==0)&&L>0)){
     this->SetNull();
 }
 this->data.AddExtRow(x, L, 1, whatN1, true, defaultVal);
 delete defaultVal;
}
void Matrix_P2A::AddLine(std::vector<double>x, int whatN1, int QDefaultValsBefore){
int QLines=this->GetQLines(), QColumns=this->GetQColumns(), L=x.size();
 Array2DSize size(QLines, QColumns);
 //template<typename T>void Array2DAddExtRow(T**&X, Array2DSize&size, T*rowParam, int wholeRowL=0, int whereFromN=1, int whatFromN=1, bool RectNotVar=false, T*DfltValParam=NULL, TValsShowHide*vsh=NULL){//44 //better but not checked
 double*defaultVal=new double;
 *defaultVal=0;
 if((QLines==1 && QColumns==1 && L>1) || ((QLines==0 || QColumns==0)&&L>0)){
     this->SetNull();
 }
 this->data.AddExtRow(x.data(), x.size(), 1, whatN1, true, defaultVal);
 delete defaultVal;
}
void Matrix_P2A::AddColumn(double*x, int L, int whatN1, int QDefaultValsBefore){
 //template<typename T>void Array2DAddIneRow(T**&X, Array2DSize&size, T*rowParam, int wholeRowL=0, T*DfltValParam=NULL, bool AddToEachNotStretchToMax=true, int whatFromN=1, int whereFromN=1
 int QLines=this->GetQLines(), QColumns=this->GetQColumns();
 Array2DSize size(QLines, QColumns);
 double*defaultVal=new double;
 *defaultVal=0;
 if((QLines==1 && QColumns==1 && L>1) || ((QLines==0 || QColumns==0)&&L>0)){
     this->SetNull();
 }
 this->data.AddIneRow(x, L, 1, whatN1, true, defaultVal);
 delete defaultVal;
}
void Matrix_P2A::AddColumn(std::vector<double>x, int FromN, int QDefaultValsBefore){
 //template<typename T>void Array2DAddIneRow(T**&X, Array2DSize&size, T*rowParam, int wholeRowL=0, T*DfltValParam=NULL, bool AddToEachNotStretchToMax=true, int whatFromN=1, int whereFromN=1){
 int QLines=this->GetQLines(), QColumns=this->GetQColumns(), L=x.size();
 Array2DSize size(QLines, QColumns);
 double*defaultVal=new double;
 *defaultVal=0;
 if( (QLines==1 && QColumns==1 && L>1) || ((QLines==0||QColumns==0)&&L>0) ){
     this->SetNull();
 }
 //template<typename T>Array2DAddIneRow(T**&X, Array2DSize& size, std::vector<T>rowParam, int IfNonRectIgnore0Add1Stretch2=0, int whatFromN=1, int QDefaultValsBefore=0, T*DefaultValParam=NULL){
 this->data.AddIneRow(      x,                                  0,           FromN,                        0,  defaultVal);
 delete defaultVal;
}
void Matrix_P2A::InsLine(int N, double*x, int L, int FromN, int QDefaultBefore){
 Array2DSize size(this->GetQLines(), this->GetQColumns());
 double*defaultVal=new double;
 *defaultVal=0;
 //template<typename T> void Array2DInsExtRow(T**&X,  Array2DSize&size, int ExtRowN, int wholeRowL=0, T*rowParam=NULL, T*DefaultValParam=NULL, int whatFromN=1, int QDefaultValuesBefore=0, bool KeepIfRect=true, TValsShowHide*vsh=NULL){
 this->data.InsExtRow(         N,             L,                 x,             defaultVal,           FromN,             QDefaultBefore,                 true);
 delete defaultVal;
}
void Matrix_P2A::InsLine(int N, std::vector<double>x, int whatFromN, int QDefaultValsBefore){
  Array2DSize size(this->GetQLines(), this->GetQColumns());
  double*defaultVal=new double;
  *defaultVal=0;
  //template<typename T> void Array2DInsExtRow(T**&X,  Array2DSize&size, int ExtRowN, std::vector<T>rowParam, T*DefaultValParam=NULL, int whatFromN=1, int QDefaultValuesBefore=0, bool KeepIfRect=true){
  this->data.InsExtRow(      N,                      x,             defaultVal,       whatFromN,        QDefaultValsBefore,                  true);
  delete defaultVal;
}
void Matrix_P2A::InsColumn(int N, double*x, int L, int FromN, int QDefaultBefore){
 Array2DSize size(this->GetQLines(), this->GetQColumns());
 double*defaultVal=new double;
 *defaultVal=0;
 //template<typename T>Array2DInsIneRow(T**&X, Array2DSize& size, int IneRowPosN, T*rowParam=NULL, int wholeRowL=0, int IfAfterLminIgnore0Stretch2=0, int whatFromN=1, int QDefaultValsBefore=0, T*DefaultValParam=NULL){
 this->data.InsIneRow(       N,               x,               L,                                0,           FromN,          QDefaultBefore,     defaultVal);
 delete defaultVal;
}
void Matrix_P2A::InsColumn(int N, std::vector<double>x, int whatFromN, int QDefaultValsBefore){
    Array2DSize size(this->GetQLines(), this->GetQColumns());
    double*defaultVal=new double;
    *defaultVal=0;
    //void InsIneRow(int IneRowPosN, std::vector<T>rowParam, int IfAfterLminIgnore0Stretch2=0, int whatFromN=1, int QDefaultValsBefore=0, T*DefaultValParam=NULL);
    this->data.InsIneRow(         N,                       x,                                0,      whatFromN,        QDefaultValsBefore,   defaultVal);
    delete defaultVal;

}

void  Matrix_P2A::DelLine(int N){
    this->data.DelExtRowN(N);
}
void Matrix_P2A::DelColumn(int N){
    this->data.DelIneRowN(N);
}



/*
Matrix_P2A Matrix_P2A::SubMatrix(int Line1N, int Col1N, int Line2N, int Col2N){
    Matrix_P2A M;*this;
    return M;
}

Matrix_P2A Matrix_P2A::MinorTo(int LineN, int ColN) const
    {
        Matrix_P2A M;
        M=*this;
        // M = this.DelLineTo(LineN);
        // M = M.DelColTo(ColN);
        M.DelLine(LineN);
        M.DelColumn(ColN);
        return M;
    }

    void Matrix_P2A::Minor(int LineN, int ColumnN)
    {
        this->DelColumn(ColumnN);
        this->DelLine(LineN);
    }
    double Matrix_P2A::AlgSuppl(int LineN, int ColumnN)  const
    {
        double d, r;
        int QLines=this->GetQLines(), QColumns=this->GetQColumns();
        Matrix_P2A M;// in |C# wa M = new Matrix();
        if ((LineN < 1) || (LineN > QLines) || (ColumnN < 1) || (ColumnN > QColumns))
        {
            r = 666;
        }
        else
        {
            M = *this;
            M = MinorTo(LineN, ColumnN);
            d = M.Determinant();
            r = d;
            if ((LineN + ColumnN)%2 != 0)
            {
                r *= (-1);
            }
        }
        return r;
    }
     double Matrix_P2A::Determinant() const
    {
        double r, d, AlgSpl;
        int RowN;
        int QColumns=this->GetQColumns(), QLines=this->GetQLines();
        if (QColumns != QLines) return 666;
        else
        {
            if (QLines == 1)
            {
                r = this->GetComponent(1, 1);
            }
            else if (QLines == 2)
            {
                r = this->GetComponent(1, 1) * this->GetComponent(2, 2) - this->GetComponent(1, 2) * this->GetComponent(2, 1);
            }
            else if (QLines == 3)
            {
                //r = this->y[1-1][1-1]       *      this->y[2-1][2-1]  *       this->y[3-1][3-1]   +   this->y[1-1][3-1]     *      this->y[2-1][1-1]   *        this->y[3-1][2-1] +  this->y[3-1][1-1]       *  this->y[1-1][2-1]       *  this->y[2-1][3-1] -
                //    this->y[1-1][3-1]       *      this->y[2-1][2-1]  *       this->y[3-1][1-1]   -   this->y[1-1][1-1]     *        this->y[2-1][3-1] *        this->y[3-1][2-1] - this->y[3-1][3-1]        * this->y[1-1][2-1]        * this->y[2-1][1-1];
                r = this->GetComponent(1, 1) * this->GetComponent(2, 2) * this->GetComponent(3, 3) + this->GetComponent(1, 3) * this->GetComponent(2, 1) * this->GetComponent(3, 2) + this->GetComponent(3, 1) * this->GetComponent(1, 2) * this->GetComponent(2, 3) -
                    this->GetComponent(1, 3) * this->GetComponent(2, 2) * this->GetComponent(3, 1) - this->GetComponent(1, 1) * this->GetComponent(2, 3) * this->GetComponent(3, 2) - this->GetComponent(3, 3) * this->GetComponent(1, 2) * this->GetComponent(2, 1);
            }
            else
            {
                RowN = 1;
                r = 0;
                for (int j = 1; j <= QColumns; j++)
                {
                    //r += GetComponent(RowN, j) * AlgSuppl(RowN, j);
                    r += GetComponent(RowN, j) * AlgSuppl(RowN, j);
                }
            }
        }
        return r;
    }

    Matrix_P2A Matrix_P2A::UnionMatrix() const
    {
        int QLines = this->GetQLines(), QColumns = this->GetQColumns();
        int TQLines = QColumns, TQColumns = QLines;
        //Matrix MT = new Matrix();
        //MT = this.TransposeTo();
        Matrix_P2A MUR;//in C# = new Matrix();
        //double[][] c;
        double CurAlgSuppl;//, det = this.Determinant();
       // c = new double[TQLines][];
        //for (int i = 1; i <= TQLines; i++)
        //{
        //    c[i - 1] = new double[TQColumns];
        //}
        MUR.SetSize(TQLines, TQColumns);
        //
        for (int i = 1; i <= QLines; i++)
        {
            for (int j = 1; j <= QColumns; j++)
            {
                CurAlgSuppl = this->AlgSuppl(i, j);
                //c[j - 1][i - 1] = CurAlgSuppl;
                MUR.SetComponent(CurAlgSuppl, j, i);
            }
        }
        //
        //MUR = new Matrix(c, TQLines, TQColumns);
        return MUR;
    }
    Matrix_P2A Matrix_P2A::InverseMatrix() const
    {
        Matrix_P2A MIR;
        double det = this->Determinant();
        MIR = this->UnionMatrix();
        if (det != 0) MIR = MIR * (1 / det);
        return MIR;
    }
*/

/*QString Matrix_A2P::ToString(){
 QString str="", cur;
 str="[[";
 int QLins=this->GetQLines(), QCols=this->GetQColumns();
 for(int i=1; i<=QLins-1; i++){
     str+="[";
     for(int j=1; j<=QCols-1; j++){
         cur.setNum(this->GetComponent(i, j));
         str+=cur;
         str+=", ";
         cur="";
     }
     cur.setNum(this->GetComponent(i, QCols));
     str+=cur;
     str+="];";
 }
 str+="[";
 //LastLine
 for(int j=1; j<=QCols-1; j++){
     cur.setNum(this->GetComponent(QLins, j));
     str+=cur;
     str+=", ";
     cur="";
 }
 //LastElt
 cur.setNum(this->GetComponent(QLins-1, QCols));
 str+=cur;
 str+="]]";
 return str;
}
std::string Matrix_A2P::ToStdString(){
 std::string ss;
 QString qs;
 qs=this->ToString();
 ss=qs.toStdString();
 return ss;
}*/


/*
Matrix_P2A Matrix_P2A::GivensRotationTransformMatrix(int N1, int N2, Matrix_P2A*M, TValsShowHide*vsh){
 int QL;
 Matrix_P2A T(QL,QL);
 if(M==NULL)M=this;
 QL=M->GetQLines();
 double a11=M->GetComponent(N1, N1), a21=M->GetComponent(N2, N1), c12, s12, denom;
 for(int i=1; i<=QL; i++){
     for(int j=1; j<=QL; j++){
         T.SetComponent(0, i, j);
     }
 }
 for(int i=1; i<=QL; i++){
     T.SetComponent(1, i, i);
 }
 if(a21==0){
     c12=1;
     s12=0;
 }else{
     denom=sqrt(a11*a11+a21*a21);
     c12=a11/denom;
     s12=a21/denom;
 }
 writeln(vsh, " T=" +M->ToString());
 writeln(vsh, " a11=" +FloatToStr(a11) + " a21=" +FloatToStr(a21)+ " c12=" +FloatToStr(c12) + " s12=" + FloatToStr(s12));
 T.SetComponent(c12, N1, N1);
 T.SetComponent(c12, N2, N2);
 T.SetComponent(s12, N1, N2);
 T.SetComponent(-s12, N2, N1);
 return T;
}
void Matrix_P2A::QRDecomposition(Matrix_P2A&Q, Matrix_P2A&R,  Matrix_P2A*M, TValsShowHide* vsh)
{
 writeln(vsh, "QRDecomposition starts working");
 if (M == NULL) M = this;
 int QL = M->GetQLines();
 int*ns=new int[QL-1];
 //int countTsInLine;
 int jn, ii, jj;
 Q.SetSize(QL, QL);
 for (int i = 1; i <= QL; i++){
     Q.SetComponent(1, i, i);
 }
 Matrix_P2A**Ts=new Matrix_P2A*[QL-1];
 for (int i = 1; i <= QL - 1; i++)
 {
     ns[i-1] = QL - 1 - i + 1;
     Ts[i - 1] = new Matrix_P2A[ns[i - 1]];
 }
 Matrix_P2A*TL=NULL;
 Matrix_P2A T;//, R=new Matrix(QL, QL);
 R.Assign(*M);
 writeln(vsh, "M="+M->ToString());
 writeln(vsh, "ini Q=" + Q.ToString());
 writeln(vsh, "ini R=" + R.ToString());
 for (int i = 1; i <= QL - 1; i++)
 {
     //countTsInLine = 0;
     for (int j = i + 1; j <= QL; j++)
     {
         //countTsInLine++; //in AddToVector
         writeln(vsh, "i=" + IntToStr(i)+" j="+IntToStr(j));
         //Matrix_VB GivensRotationTransformMatrix(int N1, int N2, Matrix_VB*M=NULL, TValsShowHide*vsh=NULL);
         T = R.GivensRotationTransformMatrix(i, j, NULL, vsh);
         writeln(vsh, "T=" + T.ToString());
         R = T * R;
         writeln(vsh, "cur R=" + R.ToString());
         jn = j - i;
         Ts[i - 1][jn - 1] = T;
     }
 }
 //R is completed, beginning wit Q
 for(int i=1; i<=QL-1; i++){
     ii = QL - 1 - i + 1;
     for (int j = 1; j <= ns[ii - 1]; j++){
         jj = ns[ii - 1] - j + 1;
         T = Ts[ii - 1][jj - 1];
         Q = Q * T;
     }
 }
 Q.Transpose();
 //
 delete[]ns;
 if(TL!=NULL)delete[]TL;
 for(int i=1; i<=QL-1; i++){
     delete[]Ts[i-1];
 }
 delete[]Ts;
 //
 writeln(vsh, "Answer:");
 writeln(vsh,"Q="+Q.ToString()+" R="+R.ToString());
 writeln(vsh, "QRDecomposition finishes working");
 //
}//fn QR decomp

Matrix_P2A Matrix_P2A::CalcMatrixOfDirCossByEulersAngles (Matrix_P2A MEulerAngles){
    Matrix_P2A MR(3, 3);
    double gamma = MEulerAngles.GetX(), //* _PI / 180,
                       psi = MEulerAngles.GetY(),// * _PI / 180,
                       theta = MEulerAngles.GetZ();// * _PI / 180;
                MR.SetComponent(cos(psi) * cos(theta), 1, 1);
                MR.SetComponent(sin(psi) * sin(gamma) - cos(psi) * sin(theta) * cos(gamma), 1, 2);
                MR.SetComponent(sin(psi) * cos(gamma) + cos(psi) * sin(theta) * sin(gamma), 1, 3);
                MR.SetComponent(sin(theta), 2, 1);
                MR.SetComponent(cos(theta) * cos(gamma), 2, 2);
                MR.SetComponent(-cos(theta) * sin(gamma), 2, 3);
                MR.SetComponent(-sin(psi) * cos(theta), 3, 1);
                MR.SetComponent(cos(psi) * sin(gamma) + sin(psi) * sin(theta) * cos(gamma), 3, 2);
                MR.SetComponent(cos(psi) * cos(gamma) - sin(psi) * sin(theta) * sin(gamma), 3, 3);
                return MR;
}

Matrix_P2A Matrix_P2A::CalcCoordsTransformFromOldCSToNew(Matrix_P2A*CoordsInOldCSParam, Matrix_P2A*EulerAnglesParam, Matrix_P2A*CSOriginOfNewCSInOldParam, TValsShowHide*vsh){
    Matrix_P2A MCoordsInOldCS(3, 1), MEulerAngles (3, 1),MCSOriginOfNewCSInOld(3, 1);
    Matrix_P2A MDirCoss(3, 3), MCoordsInNewCS(3, 1), M1(3, 3), M2(3, 3), M3(3, 3);//, M4(3, 3);
    writeln(vsh, "CalcCoordsTransformFromOldCSToNew starts working");
    if(CoordsInOldCSParam!=NULL){
        MCoordsInOldCS=(*(CoordsInOldCSParam));
    }
    writeln(vsh, "MCoordsInOldCS");
    MCoordsInOldCS.StdCOut2D();
    if(EulerAnglesParam!=NULL){
        MEulerAngles=(*(EulerAnglesParam));
    }
    writeln(vsh, "MEulerAngles");
    MCoordsInOldCS.StdCOut2D();
    if(CSOriginOfNewCSInOldParam!=NULL){
        MCSOriginOfNewCSInOld=(*(CSOriginOfNewCSInOldParam));
    }
    writeln(vsh, "MCSOriginOfNewCSInOld");
    MCSOriginOfNewCSInOld.StdCOut2D();
    MDirCoss=this->CalcMatrixOfDirCossByEulersAngles(MEulerAngles);
    writeln(vsh, "MDirCoss");
    MDirCoss.StdCOut2D();
    M1=MDirCoss.InverseMatrix();
    writeln(vsh, "MDirCoss Inverted");
    M1.StdCOut2D();
    M2=MCoordsInOldCS-MCSOriginOfNewCSInOld;
    writeln(vsh, "M2=MCoordsInOldCS-MCSOriginOfNewCSInOld");
    M2.StdCOut2D();
    writeln(vsh, "M3=M1*M2");
    M3=M1*M2;
    M3.StdCOut2D();
    MCoordsInNewCS=M3;
    writeln(vsh, "CalcCoordsTransformFromOldCSToNew finishes working");
    return MCoordsInNewCS;
    //
    //    Matrix M6, M7 = null;
    //    M1.SetSize(3, 1, 1);
    //    M3.SetSize(3, 1, 1);
    //    M2.SetSize(3, 3, 1);
    //    //M7 = M2.TransposeTo();
    //    M7 = M2.InverseMatrix();
    //    M6 = M1 - M3;
    //    M4 = M7 * M6;
    //
}

Matrix_P2A Matrix_P2A::CalcCoordsTransformFromNewCSToOld(Matrix_P2A*CoordsInNewCSParam, Matrix_P2A*EulerAnglesParam, Matrix_P2A*CSOriginOfNewCSInOldParam, TValsShowHide*vsh){
    Matrix_P2A MCoordsInNewCS(3, 1), MEulerAngles (3, 1), MCSOriginOfNewCSInOld(3, 1);
    Matrix_P2A MDirCoss(3, 3), MCoordsInOldCS(3, 1), M1(3, 3), M2(3, 3), M3(3, 3);//, M4(3, 3);
    if(CoordsInNewCSParam!=NULL){
        MCoordsInNewCS=(*(CoordsInNewCSParam));
    }
    if(EulerAnglesParam!=NULL){
        MEulerAngles=(*(EulerAnglesParam));
    }
    if(CSOriginOfNewCSInOldParam!=NULL){
        MCSOriginOfNewCSInOld=(*(CSOriginOfNewCSInOldParam));
    }
    // Matrix      CalcMatrixOfDirCossByEulersAngles(Matrix MEulerAngles);
    MDirCoss=this->CalcMatrixOfDirCossByEulersAngles(       MEulerAngles);
    M1=MDirCoss;
    M2=M1*MCoordsInNewCS;
    MCoordsInOldCS=M2+MCSOriginOfNewCSInOld;
    return MCoordsInOldCS;
    //
    //    Matrix M6 = null;
    //    M1.SetSize(3, 1, 1);
    //    M3.SetSize(3, 1, 1);
    //    M2.SetSize(3, 3, 1);
    //    M6 = M2 * M1;
    //    M4 = M6 + M3;
    //
}
*/




// }




// class Matrix_A2V{
//     double**y;
//     int QLins, QCols;
//     public:
Matrix_V2A::Matrix_V2A(int QLins, int QCols, double val){
    //this->construct();
    this->SetSize(QLins, QCols);
    this->Set(val,QLins, QCols);
}
Matrix_V2A::Matrix_V2A(double**y, int QLins, int QCols){
     //this->construct();
     this->Set(y, QLins, QCols);
}
Matrix_V2A::Matrix_V2A(std::vector<std::vector<double>>y){
     //this->construct();
     this->Set(y);
}

Matrix_V2A::Matrix_V2A(double*y, int QLines, int QColumns){
    //this->construct();
    this->Set(y, QLines, QColumns);
}
Matrix_V2A::Matrix_V2A(std::vector<double>y, int QColumns){
    //this->construct();
    this->Set(y,  QColumns);
}
Matrix_V2A::Matrix_V2A(double x, double y, double z){
    //this->construct();
    this->Set(x, y, z);
}
/*Matrix_V2A::Matrix_V2A(const Matrix_V2A&obj){
 //this->construct();
 this->Assign(obj);
}*/
//void Matrix_A2V::construct(){ this->SetNullIni(); }
void Matrix_V2A::SetNull(){
    this->data.SetNull();
}
//void Matrix_A2V::SetNullIni(){ this->data.SetNull(); }
Matrix_V2A::~Matrix_V2A(){
    this->SetNull();
}
/*void Matrix_V2A::Assign(const Matrix_V2A&obj){
    this->data=obj.data;
}
Matrix_V2A& Matrix_V2A::operator = (const Matrix_V2A&obj){
    this->Assign(obj);
    return*this;
}*/
void Matrix_V2A::SetSize(int QLins, int QCols){
    int QLinsOld=this->GetQLines(), QColsOld=this->GetQColumns();
    //void SetSize(const Array2DSize&newSize, T*DefaultValParam=NULL);
    //void SetSize(int QExtRows=1, int IneRowsLength=1, T*DefaultValParam=NULL, bool RectNotVar=true);
    this->data.SetSize(QLins, QCols, NULL, true);
    //
    if(QLins>QLinsOld || QCols>QColsOld){
        for(int i=1; i<=QLins; i++){
            for(int j=QColsOld+1; j<=QCols; j++){
                this->SetComponent(0, i, j);
            }
        }
        for(int i=QLinsOld+1; i<=QLins; i++){
            for(int j=1; j<=QCols; j++){
                this->SetComponent(0, i, j);
            }
        }
    }
}

//int GetQExtRows() const;
//template<typename T>int  Array2D_PS<T>::GetQExtRows() const
//int GetQLines() const;
int Matrix_V2A::GetQLines() const { return this->data.GetQExtRows(); }
int Matrix_V2A::GetQColumns() const { return this->data.GetMaxLength(); }

void Matrix_V2A::Set(double**y, int QLines,  int QColumns){
    Matrix_Prototype::Set(y, QLines,  QColumns);
}
void Matrix_V2A::Set(std::vector<std::vector<double>>y){
    Matrix_Prototype::Set(y);
}
void Matrix_V2A::Set(double*y, int QLines,  int QColumns){
    Matrix_Prototype::Set(y, QLines, QColumns);
}
void Matrix_V2A::Set(std::vector<double>y, int QColumns){
    Matrix_Prototype::Set(y, QColumns);
}
void Matrix_V2A::Set(double x, double y, double z){
    Matrix_Prototype::Set(x, y, z);
}
void Matrix_V2A::Set(double y, int QLines, int QColumns){
    Matrix_Prototype::Set(y, QLines, QColumns);
}

double Matrix_V2A::GetComponent(int LineN, int ColN) const{
    double y=0;
    int QLins=this->GetQLines(), QCols=this->GetQColumns();
    if(LineN>=1 && LineN<=QLins && ColN>=1 && ColN<=QCols){
        y=this->data.GetElement_AsVal(LineN, ColN);
    }
    return y;
}
void Matrix_V2A::SetComponent(double val, int LineN, int ColN){
    int QLins=this->GetQLines(), QCols=this->GetQColumns();
    if(LineN>=1 && LineN<=QLins && ColN>=1 && ColN<=QCols){
        this->data.SetElement(val, LineN, ColN);
    }
}
void Matrix_V2A::construct(){
    //NOp;
}
Matrix_Prototype* Matrix_V2A::CreateMatrix(){
    return new Matrix_V2A();
}

Matrix_Prototype* Matrix_V2A::clone(){
    return this;
}
int Matrix_V2A::GetTypeN() const{
    return MatrixTypeN_V2A;
}

/*
Matrix_V2A Matrix_V2A::VectorOfLineN(int LineN){
 Matrix_V2A M;
 bool ColOfVectorNotHorLine=true;
 std::vector<double>row;
 int QLins=this->GetQLines(), QCols=this->GetQColumns();
 if(LineN>=1 && LineN<=QLins){
     for(int i=1; i<=QCols; i++){
         row.push_back(this->GetComponent(LineN, i));
     }
 }
 M.Set(row, ColOfVectorNotHorLine);
 return M;
}
Matrix_V2A Matrix_V2A::MatrixLineN(int LineN){
 Matrix_V2A M;
 bool ColOfVectorNotHorLine=false;
 std::vector<double>row;
 int QLins=GetQLines(), QCols=GetQColumns();
 if(LineN>=1 && LineN<=QLins){
     for(int i=1; i<=QCols; i++){
         row.push_back(this->GetComponent(LineN, i));
     }
 }
 M.Set(row, ColOfVectorNotHorLine);
 return M;
}
Matrix_V2A Matrix_V2A::VectorOfColN(int ColN){
 Matrix_V2A M;
 bool ColOfVectorNotHorLine=true;
 std::vector<double>row;
 int QLins=GetQLines(), QCols=GetQColumns();
 if(ColN>=1 && ColN<=QCols){
     for(int i=1; i<=QLins; i++){
         row.push_back(this->GetComponent(i, ColN));
     }
 }
 M.Set(row, ColOfVectorNotHorLine);
 return M;
}
std::vector<double> Matrix_V2A::GetLineN(int LineN){
 std::vector<double>R;
 int QLins=GetQLines(), QCols=GetQColumns();
 if(LineN>=1 && LineN<=QLins){
     for(int i=1; i<=QCols; i++){
         R.push_back(this->GetComponent(LineN, i));
     }
 }
 return R;
}
std::vector<double>Matrix_V2A::GetColN(int ColN){
    std::vector<double>R;
    int QLins=GetQLines(), QCols=GetQColumns();
    if(ColN>=1 && ColN<=QCols){
    for(int i=1; i<=QLins; i++){
        R.push_back(this->GetComponent(i, ColN));
    }
 }
 return R;
}
Matrix_V2A Matrix_V2A::operator +(const Matrix_V2A&obj){
    double x;
    int QLins=GetQLines(), QCols=GetQColumns();
    if(this->GetQLines()==obj.GetQLines() && this->GetQColumns()==obj.GetQColumns()){
        for(int i=1; i<=QLins; i++){
            for(int j=1; j<=QCols; j++){
                x=this->GetComponent(i, j)+obj.GetComponent(i, j);
                 this->SetComponent(x, i, j);
            }
        }
    }
    return *this;
}
Matrix_V2A Matrix_V2A::operator -(const Matrix_V2A&obj){
    double x;
    int QLins=GetQLines(), QCols=GetQColumns();
    if(this->GetQLines()==obj.GetQLines() && this->GetQColumns()==obj.GetQColumns()){
        for(int i=1; i<=QLins; i++){
            for(int j=1; j<=QCols; j++){
                x=this->GetComponent(i, j)-obj.GetComponent(i, j);
                this->SetComponent(x, i, j);
            }
        }
    }
    return *this;
}
Matrix_V2A Matrix_V2A::operator *(const Matrix_V2A&obj){
    //double x;
    //int QM, QLins=GetQLines(), QCols=GetQColumns();
    //if(this->GetQColumns()==obj.GetQLines()){
    //    QM=this->GetQColumns();
    //    for(int i=1; i<=QLins; i++){
    //        for(int j=1; j<=QCols; j++){
    //            for(int k=1; k<=QM; k++){
    //                x=this->GetComponent(i, k)-obj.GetComponent(k, j);
    //                this->SetComponent(x, i, j);
    //            }
    //        }
    //    }
    //}
    //return *this;
    double x;
    Matrix_V2A MT=obj;
    Matrix_V2A M1=(*(this));
    //Matrix *M2=&obj;//invalid conversion from const Matrix_PS* to Matrix_PS* [-fpermissive]
    //Matrix *M2=obj;//cannot convert const Matrix_PS to Matrix_PS* in initialization
    Matrix_V2A *M2;
    //M2=&obj;
    //M2=obj;//cannot convert const Matrix to Matrix* in assignment
    //M2=MT;//cannot convert const Matrix to Matrix* in assignment
    M2=&MT;
    //Matrix M3;
    int QL1, QL2, QL3, QC1, QC2, QC3, QM;
    QL1=this->GetQLines(), QL2=obj.GetQLines(), QC1=this->GetQColumns(), QC2=obj.GetQColumns();
    double x1, x2, x3, mp;
    if(QC1==QL2){
        QM=QC1;
        QM=QL2;
        QL3=QL1;
        QC3=QC2;
        this->SetSize(QL3, QC3);
        //M3.SetSize(QL3, QC3);
        for(int i=1; i<=QL3; i++){
            for(int j=1; j<=QC3; j++){
                x3=0;
                for(int k=1; k<=QM; k++){
                    x1=M1.GetComponent(i, k);
                    x2=M2->GetComponent(k, j);
                    mp=x1*x2;
                    x3+=mp;
                }
                //M3.SetComponent(x3, i, j);
                this->SetComponent(x3, i, j);
            }
        }
    }
    //(*(this))=M3;
    return *this;

}
Matrix_V2A Matrix_V2A::operator /(const Matrix_V2A&obj){
    Matrix_V2A MR, MI;
    MI = obj.InverseMatrix();
    MR = MI * (*(this));
    return MR;
}

Matrix_V2A Matrix_V2A::operator *(double x){
    double cur;
    int QLins=GetQLines(), QCols=GetQColumns();
    for(int i=1; i<=QLins; i++){
        for(int j=1; j<=QCols; j++){
            cur=this->GetComponent(i, j)*x;
            this->SetComponent(cur, i, j);
        }
    }
    return *this;
}

void Matrix_V2A::Transpose(){
    this->data.Transpose();
}
Matrix_V2A Matrix_V2A::TransposeTo(){
    Matrix_V2A M=*this;
    M.Transpose();
    return M;
}

void Matrix_V2A::SetX(double val){
    if(this->GetQLines()>=1){
        this->SetComponent(val, 1, 1);
    }
}

void Matrix_V2A::SetY(double val){
    if(this->GetQLines()>=2){
        this->SetComponent(val, 2, 1);
    }
}
void Matrix_V2A::SetZ(double val){
    if(this->GetQLines()>=3){
        this->SetComponent(val, 3, 1);
    }
}
double Matrix_V2A::GetX(){
    return this->GetComponent(1, 1);
}

double Matrix_V2A::GetY(){
   double y;
   if(this->GetQLines()>=2){
       y=this->GetComponent(2, 1);
   }
   return y;
}

double Matrix_V2A::GetZ(){
    double y;
    if(this->GetQLines()>=3){
        y=this->GetComponent(3, 1);
    }
    return y;
}
*/

/*void Matrix::AddLine(std::vector<double>x, int N1){
 int L, Lmin, whereN, whatN;
 this->SetSize(this->QLins+1, this->QCols);
 L=x.size();
 //Lmin = L <= this->QCols ? L : this->QCols;//uz N1=1
 Lmin = L <= this->QCols-N1+1 ? L : this->QCols-N1+1;
 for(int i=1; i<=N1-1; i++){
     this->y[this->QLins-1][i-1]=0;
 }
 for(int i=1; i<=Lmin; i++){
     this->y[this->QLins-1][N1-1+i-1]=x[i-1];
 }
 for(int i=Lmin+1; i<=this->QCols; i++){
     this->y[this->QLins-1][i-1]=0;
 }
 //return *this;
}
void Matrix::AddColumn(std::vector<double>x, int N1){
 int L, Lmin, whereN, whatN;
 this->SetSize(this->QLins, this->QCols+1);
 L=x.size();
 //Lmin = L <= this->QCols ? L : this->QCols;//uz N1=1
 Lmin = L <= this->QLins-N1+1 ? L : this->QLins-N1+1;
 for(int i=1; i<=N1-1; i++){
     this->y[i-1][this->QCols-1]=0;
 }
 for(int i=1; i<=Lmin; i++){
     this->y[N1-1+i-1][this->QCols-1]=x[i-1];
 }
 for(int i=Lmin+1; i<=this->QLins; i++){
     this->y[i-1][this->QCols-1]=0;
 }
 //return *this;
}
void Matrix::InsLine(int N, std::vector<double>x, int N1){
 int L, Lmin;//, whereN, whatN;
 if(N>=1 && N<=this->QLins){
     this->SetSize(this->QLins+1, this->QCols);
     for(int i=QLins-1; i>=N; i--){
         for(int j=1; j<=this->QCols; j++){
             y[i+1-1][j-1]= y[i-1][j-1];
         }
     }
     L=x.size();
     //Lmin = L <= this->QCols ? L : this->QCols;//uz N1=1
     Lmin = L <= this->QCols-N1+1 ? L : this->QCols-N1+1;
     for(int i=1; i<=N1-1; i++){
         this->y[N-1][i-1]=0;
     }
     for(int i=1; i<=Lmin; i++){
         this->y[N-1][N1-1+i-1]=x[i-1];
     }
     for(int i=Lmin+1; i<=this->QCols; i++){
         this->y[N-1][i-1]=0;
     }
 }
 //return *this;
}
void Matrix::InsColumn(int N, std::vector<double>x, int N1){
 int L, Lmin;//, whereN, whatN;
 if(N>=1 && N<=this->QCols){
     this->SetSize(this->QLins, this->QCols+1);
     for(int i=QCols-1; i>=N; i--){
         for(int j=1; j<=this->QLins; j++){
             y[j-1][i+1-1]= y[j-1][i-1];
         }
     }
     L=x.size();
     //Lmin = L <= this->QCols ? L : this->QCols;//uz N1=1
     Lmin = L <= this->QCols-N1+1 ? L : this->QCols-N1+1;
     for(int i=1; i<=N1-1; i++){
         this->y[N-1][i-1]=0;
     }
     for(int i=1; i<=Lmin; i++){
         this->y[N-1][N1-1+i-1]=x[i-1];
     }
     for(int i=Lmin+1; i<=this->QCols; i++){
         this->y[N-1][i-1]=0;
     }
 }
 //return *this;
}*/
//
void Matrix_V2A::SetLine(int N, double*x, int L, int FromN, int QDefaultBefore, bool DefaultValsAreZerosNotOwn){
 //std::vector<double>row;
 int QL=this->GetQLines(), QC=this->GetQColumns();
 Array2DSize size(QL, QC);
 //int preN1_, preN2_, ownN1_, ownN2, ForOwnN1_, ForOwnN2, postN1, postN2_;
 //int LreqCalcd;
 //void Calc_Array1DSubArray_byLs_MarkupNs_vFmls(int&LreqCalcd, int&preN1_, int&preN2_, int&ownN1_, int&ownN2, int&ForOwnN1_, int&ForOwnN2, int&postN1, int&postN2_, int Lini, int LreqGiven, int whatFromN, int QDefaultValuesBefore){
 //Calc_Array1DSubArray_byLs_MarkupNs_vFmls(LreqCalcd, preN1_, preN2_, ownN1_, ownN2, ForOwnN1_, ForOwnN2, postN1, postN2_, QC, 0, whatFromN, QDefaultValuesBefore);
 double *dflt=new double;
 (*(dflt))=0;
 //template<typename T> void Array2DSetExtRowN(T**&X, Array2DSize&size, int ExtRowN, int wholeRowL=0, T*rowParam=NULL, int whatFromN=1, int QDefaultValuesBefore=0, bool KeepIfRect=true, T*DefaultValParam=NULL, bool DefaultIsSpecNotOwn=false){
 this->data.SetExtRow(                          N,               L,               x,           FromN,             QDefaultBefore,                false,                   dflt,       DefaultValsAreZerosNotOwn);
 delete dflt;
}
void  Matrix_V2A::SetLine(int N, std::vector<double>x, int FromN, int QDfltValsBefore, bool DefaultValsAreOwnNotZeros){
 //template<typename T> void Array2DSetExtRowN(std::vector<std::vector<T>>X, int ExtRowN, std::vector<T>rowParam, T*DefaultValParam=NULL, int FromN=1, int DefaultValuesBefore=0, bool KeepIfRect=false)
 //template<typename T> void Array2DSetExtRowN(T**&X, const Array2DSize&size, int ExtRowN, std::vector<T>rowParam, T*DefaultValParam=NULL, int FromN=1, int DefaultValuesBefore=0, bool KeepIfRect=false)
 //template<typename T> void Array2DSetExtRowN(T**&X, const Array2DSize&size, int ExtRowN, int wholeRowL=0, T*rowParam=NULL, T*DefaultValParam=NULL, int FromN=1, int DefaultValuesBefore=0, bool KeepIfRect=false)
 Array2DSize size(this->GetQLines(), this->GetQColumns());
 //std::vector<double> GetNumericSubarray_byLs(std::vector<double> x, int FromN=1, int QDefaultValuesBefore=0, int LreqGiven=0);
 //std::vector<double> GetNumericSubarray_byLs(std::vector<double> x, std::vector<double> ArrOfDfltVals, int FromN=1, int QDefaultValuesBefore=0, int LreqGiven=0);
 std::vector<double>row=GetNumericSubarray_byLs(x, FromN, QDfltValsBefore, this->GetQColumns());
 //double*defaultVal=new double;
 //*defaultVal=0;
 //Array2DSetExtRowN(T**&X,  const Array2DSize&size, int ExtRowN, std::vector<T>rowParam, T*DefaultValParam=NULL, int FromN=1, int DefaultValuesBefore=0, bool KeepIfRect=false)
 //void SetExtRow(int ExtRowN, std::vector<T>rowParam, int whatFromN=1, int QDefaultValuesBefore=0, bool KeepIfRect=true, T*DefaultValParam=NULL, bool DefaultIsSpecNotOwn=false);
 this->data.SetExtRow(      N,                    row,         FromN,       QDfltValsBefore,                        true);
 //if(x.siz)
 //delete defaultVal;
}
void  Matrix_V2A::SetColumn(int N, double*x, int L, int FromN, int QDefaultBefore, bool DefaultValsAreOwnNotZeros){
 int QL=this->GetQLines(), QC=this->GetQColumns();
 Array2DSize size(QL,  QC);
 double*dflt=new double;
 (*(dflt))=0;
 //template<typename T>Array2DSetIneRow(T**&X, Array2DSize& size, int IneRowN, T*rowParam=NULL, int wholeRowL=0, int whatFromN=1, int QDefaultValsBefore=0, T*DefaultValParam=NULL, bool ExistingInsteadOfDefault=true, bool ExistingOnlyForGTLminNotIgnore=false, bool LastPossIfNotRectAtPos0NotIgnore=false){
 this->data.SetIneRow(           N,               x,               L,           FromN,       QDefaultBefore,             dflt,                    DefaultValsAreOwnNotZeros);//            ExistingOnlyForGTLminNotIgnore,      LastPossIfNotRectAtPos0NotIgnore);
 delete dflt;
}
void  Matrix_V2A::SetColumn(int N, std::vector<double>x, int whatN1, int QDfltValsBefore, bool DefaultValsAreOwnNotZeros){
 //template<typename T> void Array2DSetIneRowN(T**&X, Array2DSize&size, int IneRowN, int wholeRowL=0, T*rowParam=NULL, T*DefaultValParam=NULL, int FromN=1, int DefaultValuesBefore=0, int IfNGTLmin_Ignore0DoNil1Add2Srtetch3=1)
 //template<typename T> void Array2DSetIneRowN(T**&X, Array2DSize&size, int IneRowN, std::vector<T>rowParam, T*DefaultValExt=NULL, int FromN=1, int DefaultValuesBefore=0, int IfNGTLmin_Ignore0DoNil1Add2Srtetch3=1)
 //template<typename T> void Array2DSetIneRowN(std::vector<std::vector<T>>&X, std::vector<T>rowParam, int IneRowN, T*DefaultValExt=NULL, int FromN=1, int DefaultValuesBefore=0, int IfNGTLmin_Ignore0DoNil1Add2Srtetch3=1)
 Array2DSize size(this->GetQLines(), this->GetQColumns());
 double*defaultVal=new double;
 *defaultVal=0;
 //template<typename T>Array2DSetIneRow(T**&X, Array2DSize& size, int IneRowN, std::vector<T>rowParam, int whatFromN=1, int QDefaultValsBefore=0, T*DefaultValParam=NULL, bool ExistingInsteadOfDefault=true, bool ExistingOnlyForGTLminNotIgnore=false, bool LastPossIfNotRectAtPos0NotIgnore=false){
 this->data.SetIneRow(           N,                      x, whatN1,                   QDfltValsBefore, defaultVal);
 delete defaultVal;
}
//
void Matrix_V2A::AddLine(double*x, int L, int whatN1, int QDefaultValsBefore){
 int QLines=this->GetQLines(), QColumns=this->GetQColumns();
 //std::vector<double> GetNumericSubarray_byLs(std::vector<double> x, int FromN=1, int QDefaultValuesBefore=0, int LreqGiven=0);
 //std::vector<double> GetNumericSubarray_byLs(std::vector<double> x, std::vector<double> ArrOfDfltVals, int FromN=1, int QDefaultValuesBefore=0, int LreqGiven=0);
 std::vector<double>row0, row;
 Array1DAssign(row0, x, L);
 row=GetNumericSubarray_byLs(row0, whatN1, QDefaultValsBefore, QColumns);
 if((QLines==1 && QColumns==1 && L>1) || ((QLines==0 || QColumns==0)&&L>0)){
     this->SetNull();
 }
 this->data.AddExtRow(row);
}
void Matrix_V2A::AddLine(std::vector<double>x, int whatN1, int QDefaultValsBefore){
int QLines=this->GetQLines(), QColumns=this->GetQColumns(), L=x.size();
 Array2DSize size(QLines, QColumns);
 //template<typename T>void Array2DAddExtRow(T**&X, Array2DSize&size, T*rowParam, int wholeRowL=0, int whereFromN=1, int whatFromN=1, bool RectNotVar=false, T*DfltValParam=NULL, TValsShowHide*vsh=NULL){//44 //better but not checked
 double*defaultVal=new double;
 *defaultVal=0;
 if((QLines==1 && QColumns==1 && L>1) || ((QLines==0 || QColumns==0)&&L>0)){
     this->SetNull();
 }
 //template<typename T> void Array2D_PS<T>::AddExtRow(std::vector<T>rowParam, int whatFromN, int QDefaultValsBefore, bool RectNotVar, T*DfltValParam, TValsShowHide*vsh)
 this->data.AddExtRow(x, 1, whatN1, true, defaultVal);
 delete defaultVal;
}
void Matrix_V2A::AddColumn(double*x, int L, int whatN1, int QDefaultValsBefore){
 //template<typename T>void Array2DAddIneRow(T**&X, Array2DSize&size, T*rowParam, int wholeRowL=0, T*DfltValParam=NULL, bool AddToEachNotStretchToMax=true, int whatFromN=1, int whereFromN=1
 int QLines=this->GetQLines(), QColumns=this->GetQColumns();
 Array2DSize size(QLines, QColumns);
 double*defaultVal=new double;
 *defaultVal=0;
 if((QLines==1 && QColumns==1 && L>1) || ((QLines==0 || QColumns==0)&&L>0)){
     this->SetNull();
 }
 this->data.AddIneRow(x, L, 1, whatN1, true, defaultVal);
 delete defaultVal;
}
void Matrix_V2A::AddColumn(std::vector<double>x, int FromN, int QDefaultValsBefore){
 //template<typename T>void Array2DAddIneRow(T**&X, Array2DSize&size, T*rowParam, int wholeRowL=0, T*DfltValParam=NULL, bool AddToEachNotStretchToMax=true, int whatFromN=1, int whereFromN=1){
 int QLines=this->GetQLines(), QColumns=this->GetQColumns(), L=x.size();
 Array2DSize size(QLines, QColumns);
 double*defaultVal=new double;
 *defaultVal=0;
 if( (QLines==1 && QColumns==1 && L>1) || ((QLines==0||QColumns==0)&&L>0) ){
     this->SetNull();
 }
 //template<typename T>Array2DAddIneRow(T**&X, Array2DSize& size, std::vector<T>rowParam, int IfNonRectIgnore0Add1Stretch2=0, int whatFromN=1, int QDefaultValsBefore=0, T*DefaultValParam=NULL){
 this->data.AddIneRow(      x,                                  0,           FromN,                        0,  defaultVal);
 delete defaultVal;
}

void Matrix_V2A::InsLine(int N, double*x, int L, int FromN, int QDefaultBefore){
 Array2DSize size(this->GetQLines(), this->GetQColumns());
 double*defaultVal=new double;
 *defaultVal=0;
 //template<typename T> void Array2DInsExtRow(T**&X,  Array2DSize&size, int ExtRowN, int wholeRowL=0, T*rowParam=NULL, T*DefaultValParam=NULL, int whatFromN=1, int QDefaultValuesBefore=0, bool KeepIfRect=true, TValsShowHide*vsh=NULL){
 this->data.InsExtRow(         N,             L,                 x,             defaultVal,           FromN,             QDefaultBefore,                 true);
 delete defaultVal;
}
void Matrix_V2A::InsLine(int N, std::vector<double>x, int whatFromN, int QDefaultValsBefore){
  Array2DSize size(this->GetQLines(), this->GetQColumns());
  double*defaultVal=new double;
  *defaultVal=0;
  //template<typename T> void Array2DInsExtRow(T**&X,  Array2DSize&size, int ExtRowN, std::vector<T>rowParam, T*DefaultValParam=NULL, int whatFromN=1, int QDefaultValuesBefore=0, bool KeepIfRect=true){
  this->data.InsExtRow(      N,                      x,             defaultVal,       whatFromN,        QDefaultValsBefore,                  true);
  delete defaultVal;
}
void Matrix_V2A::InsColumn(int N, double*x, int L, int FromN, int QDefaultBefore){
 Array2DSize size(this->GetQLines(), this->GetQColumns());
 double*defaultVal=new double;
 *defaultVal=0;
 //template<typename T>Array2DInsIneRow(T**&X, Array2DSize& size, int IneRowPosN, T*rowParam=NULL, int wholeRowL=0, int IfAfterLminIgnore0Stretch2=0, int whatFromN=1, int QDefaultValsBefore=0, T*DefaultValParam=NULL){
 this->data.InsIneRow(       N,               x,               L,                                0,           FromN,          QDefaultBefore,     defaultVal);
 delete defaultVal;
}
void Matrix_V2A::InsColumn(int N, std::vector<double>x, int whatFromN, int QDefaultValsBefore){
 Array2DSize size(this->GetQLines(), this->GetQColumns());
 double*defaultVal=new double;
 *defaultVal=0;
 //template<typename T>Array2DInsIneRow(T**&X, Array2DSize& size, int IneRowPosN, std::vector<T>rowParam, int IfAfterLminIgnore0Stretch2=0, int whatFromN=1, int QDefaultValsBefore=0, T*DefaultValParam=NULL){
 //void InsIneRow(int IneRowPosN, T*rowParam=NULL, int wholeRowL=0, int IfAfterLminIgnore0Stretch2=0, int whatFromN=1, int QDefaultValsBefore=0, T*DefaultValParam=NULL);
 //void InsIneRow(T**&X, Array2DSize& size, int IneRowPosN, std::vector<T>rowParam, int IfAfterLminIgnore0Stretch2=0, int whatFromN=1, int QDefaultValsBefore=0, T*DefaultValParam=NULL);
 this->data.InsIneRow(      N,                      x,                                0,       whatFromN,       QDefaultValsBefore,             defaultVal);
}
void Matrix_V2A::DelLine(int N){
 Array2DSize size(this->GetQLines(), this->GetQColumns());
 //Array1DDelElementFromN(this->y, QLines, N);
 //
 this->data.DelExtRowN(N);
}
void Matrix_V2A::DelColumn(int N){
 Array2DSize size(this->GetQLines(), this->GetQColumns());
 //template<typename T>void Array2DDelIneRow(T**&X, Array2DSize&size, int N)
 //
 this->data.DelIneRowN(N);
}

/*
Matrix_V2A Matrix_V2A::SubMatrix(int Line1N, int Col1N, int Line2N, int Col2N){
 Matrix_V2A M;
 int QLins=this->GetQLines(), QCols=this->GetQColumns(), QLinsNew=0, QColsNew=0, LineN, ColN;
 double val;
 if(Line1N>=1 && Line1N<=QLins && Line2N>=1 && Line2N<=QLins && Col1N>=1 && Col1N<=QCols && Col2N>=1 && Col2N<=QCols){
     if(Line1N<=Line2N  && Col1N<=Col2N){
         QLinsNew=Line2N-Line1N+1;
         QColsNew=Col2N-Col1N+1;
         M.SetSize(QLinsNew, QColsNew);
         //for(int i=Line1N; i<=Line2N; i++){
         for(int i=1; i<=QLinsNew; i++){
             //LineN=i-Line1N+1;
             LineN=i+Line1N-1;
             for(int j=Col1N; j<=Col2N; j++){
                 ColN=Col1N+j-1;
                 //ColN=j-Col1N+1;
                 val=this->GetComponent(LineN, ColN);
                 //val=this->GetElement(i, j);
                 M.SetComponent(val, i, j);
                 //M.SetElement(val, LineN, ColN);
             }
         }
     }else if(Line1N<=Line2N  && Col1N>Col2N){
         QLinsNew=Line2N-Line1N+1;
         QColsNew=Col1N-Col2N+1;
         M.SetSize(QLinsNew, QColsNew);
         for(int i=1; i<=QLinsNew; i++){
             LineN=i+Line1N-1;
             for(int j=1; j>=QColsNew; j++){
                 ColN=Col2N+j-1;//(changes to Col1N)
                 val=this->GetComponent(LineN, ColN);
                 M.SetComponent(val, i, j);
             }
         }
     }else if(Line1N>=Line2N  && Col1N<=Col2N){
         QLinsNew=Line1N-Line2N+1;
         QColsNew=Col2N-Col1N+1;
         M.SetSize(QLinsNew, QColsNew);
         for(int i=1; i<=QLinsNew; i++){
             LineN=Line2N+i-1;//(changes to Line1N)
             for(int j=1; j>=QColsNew; j++){
                 ColN=Col1N+j-1;
                 val=this->GetComponent(LineN, ColN);
                 M.SetComponent(val, i, j);
             }
         }
     }else{//both 1 > 2
         QLinsNew=Line1N-Line2N+1;
         QColsNew=Col1N-Col2N+1;
         M.SetSize(QLinsNew, QColsNew);
         for(int i=1; i<=QLinsNew; i++){
             LineN=Line2N+i-1;//(changes to Line1N)
             for(int j=1; j>=QColsNew; j++){
                 ColN=Col2N+j-1;//(changes to Col1N)
                 val=this->GetComponent(LineN, ColN);
                 M.SetComponent(val, i, j);
             }
         }
     }
 }
 return M;
}



Matrix_V2A Matrix_V2A::MinorTo(int LineN, int ColN) const
    {
        Matrix_V2A M;
        M=*this;
        // M = this.DelLineTo(LineN);
        // M = M.DelColTo(ColN);
        M.DelLine(LineN);
        M.DelColumn(ColN);
        return M;
    }

    void Matrix_V2A::Minor(int LineN, int ColumnN)
    {
        this->DelColumn(ColumnN);
        this->DelLine(LineN);
    }
    double Matrix_V2A::AlgSuppl(int LineN, int ColumnN)  const
    {
        double d, r;
        int QLines=this->GetQLines(), QColumns=this->GetQColumns();
        Matrix_V2A M;// in |C# wa M = new Matrix();
        if ((LineN < 1) || (LineN > QLines) || (ColumnN < 1) || (ColumnN > QColumns))
        {
            r = 666;
        }
        else
        {
            M = *this;
            M = MinorTo(LineN, ColumnN);
            d = M.Determinant();
            r = d;
            if ((LineN + ColumnN)%2 != 0)
            {
                r *= (-1);
            }
        }
        return r;
    }
     double Matrix_V2A::Determinant() const
    {
        double r, d, AlgSpl;
        int RowN;
        int QColumns=this->GetQColumns(), QLines=this->GetQLines();
        if (QColumns != QLines) return 666;
        else
        {
            if (QLines == 1)
            {
                r = this->GetComponent(1, 1);
            }
            else if (QLines == 2)
            {
                r = this->GetComponent(1, 1) * this->GetComponent(2, 2) - this->GetComponent(1, 2) * this->GetComponent(2, 1);
            }
            else if (QLines == 3)
            {
                //r = this->y[1-1][1-1]       *      this->y[2-1][2-1]  *       this->y[3-1][3-1]   +   this->y[1-1][3-1]     *      this->y[2-1][1-1]   *        this->y[3-1][2-1] +  this->y[3-1][1-1]       *  this->y[1-1][2-1]       *  this->y[2-1][3-1] -
                //    this->y[1-1][3-1]       *      this->y[2-1][2-1]  *       this->y[3-1][1-1]   -   this->y[1-1][1-1]     *        this->y[2-1][3-1] *        this->y[3-1][2-1] - this->y[3-1][3-1]        * this->y[1-1][2-1]        * this->y[2-1][1-1];
                r = this->GetComponent(1, 1) * this->GetComponent(2, 2) * this->GetComponent(3, 3) + this->GetComponent(1, 3) * this->GetComponent(2, 1) * this->GetComponent(3, 2) + this->GetComponent(3, 1) * this->GetComponent(1, 2) * this->GetComponent(2, 3) -
                    this->GetComponent(1, 3) * this->GetComponent(2, 2) * this->GetComponent(3, 1) - this->GetComponent(1, 1) * this->GetComponent(2, 3) * this->GetComponent(3, 2) - this->GetComponent(3, 3) * this->GetComponent(1, 2) * this->GetComponent(2, 1);
            }
            else
            {
                RowN = 1;
                r = 0;
                for (int j = 1; j <= QColumns; j++)
                {
                    //r += GetComponent(RowN, j) * AlgSuppl(RowN, j);
                    r += GetComponent(RowN, j) * AlgSuppl(RowN, j);
                }
            }
        }
        return r;
    }

    Matrix_V2A Matrix_V2A::UnionMatrix() const
    {
        int QLines = this->GetQLines(), QColumns = this->GetQColumns();
        int TQLines = QColumns, TQColumns = QLines;
        //Matrix MT = new Matrix();
        //MT = this.TransposeTo();
        Matrix_V2A MUR;//in C# = new Matrix();
        //double[][] c;
        double CurAlgSuppl;//, det = this.Determinant();
       // c = new double[TQLines][];
        //for (int i = 1; i <= TQLines; i++)
        //{
        //    c[i - 1] = new double[TQColumns];
        //}
        MUR.SetSize(TQLines, TQColumns);
        //
        for (int i = 1; i <= QLines; i++)
        {
            for (int j = 1; j <= QColumns; j++)
            {
                CurAlgSuppl = this->AlgSuppl(i, j);
                //c[j - 1][i - 1] = CurAlgSuppl;
                MUR.SetComponent(CurAlgSuppl, j, i);
            }
        }
        //
        //MUR = new Matrix(c, TQLines, TQColumns);
        return MUR;
    }
    Matrix_V2A Matrix_V2A::InverseMatrix() const
    {
        Matrix_V2A MIR;
        double det = this->Determinant();
        MIR = this->UnionMatrix();
        if (det != 0) MIR = MIR * (1 / det);
        return MIR;
    }

*/

/*QString Matrix_A2V::ToString(){
 QString str="", cur;
 str="[[";
 int QLins=this->GetQLines(), QCols=this->GetQColumns();
 for(int i=1; i<=QLins-1; i++){
     str+="[";
     for(int j=1; j<=QCols-1; j++){
         cur.setNum(this->GetComponent(i, j));
         str+=cur;
         str+=", ";
         cur="";
     }
     cur.setNum(this->GetComponent(i, QCols));
     str+=cur;
     str+="];";
 }
 str+="[";
 //LastLine
 for(int j=1; j<=QCols-1; j++){
     cur.setNum(this->GetComponent(QLins, j));
     str+=cur;
     str+=", ";
     cur="";
 }
 //LastElt
 cur.setNum(this->GetComponent(QLins-1, QCols));
 str+=cur;
 str+="]]";
 return str;
}
std::string Matrix_A2V::ToStdString(){
 std::string ss;
 QString qs;
 qs=this->ToString();
 ss=qs.toStdString();
 return ss;
}
*/

/*
Matrix_V2A Matrix_V2A::GivensRotationTransformMatrix(int N1, int N2, Matrix_V2A*M, TValsShowHide*vsh){
 int QL;
 Matrix_V2A T(QL,QL);
 if(M==NULL)M=this;
 QL=M->GetQLines();
 double a11=M->GetComponent(N1, N1), a21=M->GetComponent(N2, N1), c12, s12, denom;
 for(int i=1; i<=QL; i++){
     for(int j=1; j<=QL; j++){
         T.SetComponent(0, i, j);
     }
 }
 for(int i=1; i<=QL; i++){
     T.SetComponent(1, i, i);
 }
 if(a21==0){
     c12=1;
     s12=0;
 }else{
     denom=sqrt(a11*a11+a21*a21);
     c12=a11/denom;
     s12=a21/denom;
 }
 writeln(vsh, " T=" +M->ToString());
 writeln(vsh, " a11=" +FloatToStr(a11) + " a21=" +FloatToStr(a21)+ " c12=" +FloatToStr(c12) + " s12=" + FloatToStr(s12));
 T.SetComponent(c12, N1, N1);
 T.SetComponent(c12, N2, N2);
 T.SetComponent(s12, N1, N2);
 T.SetComponent(-s12, N2, N1);
 return T;
}
void Matrix_V2A::QRDecomposition(Matrix_V2A&Q, Matrix_V2A&R,  Matrix_V2A*M, TValsShowHide* vsh)
{
 writeln(vsh, "QRDecomposition starts working");
 if (M == NULL) M = this;
 int QL = M->GetQLines();
 int*ns=new int[QL-1];
 //int countTsInLine;
 int jn, ii, jj;
 Q.SetSize(QL, QL);
 for (int i = 1; i <= QL; i++){
     Q.SetComponent(1, i, i);
 }
 Matrix_V2A**Ts=new Matrix_V2A*[QL-1];
 for (int i = 1; i <= QL - 1; i++)
 {
     ns[i-1] = QL - 1 - i + 1;
     Ts[i - 1] = new Matrix_V2A[ns[i - 1]];
 }
 Matrix_V2A*TL=NULL;
 Matrix_V2A T;//, R=new Matrix(QL, QL);
 R.Assign(*M);
 writeln(vsh, "M="+M->ToString());
 writeln(vsh, "ini Q=" + Q.ToString());
 writeln(vsh, "ini R=" + R.ToString());
 for (int i = 1; i <= QL - 1; i++)
 {
     //countTsInLine = 0;
     for (int j = i + 1; j <= QL; j++)
     {
         //countTsInLine++; //in AddToVector
         writeln(vsh, "i=" + IntToStr(i)+" j="+IntToStr(j));
         //Matrix_VB GivensRotationTransformMatrix(int N1, int N2, Matrix_VB*M=NULL, TValsShowHide*vsh=NULL);
         T = R.GivensRotationTransformMatrix(i, j, NULL, vsh);
         writeln(vsh, "T=" + T.ToString());
         R = T * R;
         writeln(vsh, "cur R=" + R.ToString());
         jn = j - i;
         Ts[i - 1][jn - 1] = T;
     }
 }
 //R is completed, beginning wit Q
 for(int i=1; i<=QL-1; i++){
     ii = QL - 1 - i + 1;
     for (int j = 1; j <= ns[ii - 1]; j++){
         jj = ns[ii - 1] - j + 1;
         T = Ts[ii - 1][jj - 1];
         Q = Q * T;
     }
 }
 Q.Transpose();
 //
 delete[]ns;
 if(TL!=NULL)delete[]TL;
 for(int i=1; i<=QL-1; i++){
     delete[]Ts[i-1];
 }
 delete[]Ts;
 //
 writeln(vsh, "Answer:");
 writeln(vsh,"Q="+Q.ToString()+" R="+R.ToString());
 writeln(vsh, "QRDecomposition finishes working");
 //
}//fn QR decomp


Matrix_V2A Matrix_V2A::CalcMatrixOfDirCossByEulersAngles (Matrix_V2A MEulerAngles){
    Matrix_V2A MR(3, 3);
    double gamma = MEulerAngles.GetX() ,//* _PI / 180,
                       psi = MEulerAngles.GetY(),// * _PI / 180,
                       theta = MEulerAngles.GetZ();// * _PI / 180;
                MR.SetComponent(cos(psi) * cos(theta), 1, 1);
                MR.SetComponent(sin(psi) * sin(gamma) - cos(psi) * sin(theta) * cos(gamma), 1, 2);
                MR.SetComponent(sin(psi) * cos(gamma) + cos(psi) * sin(theta) * sin(gamma), 1, 3);
                MR.SetComponent(sin(theta), 2, 1);
                MR.SetComponent(cos(theta) * cos(gamma), 2, 2);
                MR.SetComponent(-cos(theta) * sin(gamma), 2, 3);
                MR.SetComponent(-sin(psi) * cos(theta), 3, 1);
                MR.SetComponent(cos(psi) * sin(gamma) + sin(psi) * sin(theta) * cos(gamma), 3, 2);
                MR.SetComponent(cos(psi) * cos(gamma) - sin(psi) * sin(theta) * sin(gamma), 3, 3);
                return MR;
}

Matrix_V2A Matrix_V2A::CalcCoordsTransformFromOldCSToNew(Matrix_V2A*CoordsInOldCSParam, Matrix_V2A*EulerAnglesParam, Matrix_V2A*CSOriginOfNewCSInOldParam, TValsShowHide*vsh){
    Matrix_V2A MCoordsInOldCS(3, 1), MEulerAngles (3, 1),MCSOriginOfNewCSInOld(3, 1);
    Matrix_V2A MDirCoss(3, 3), MCoordsInNewCS(3, 1), M1(3, 3), M2(3, 3), M3(3, 3);//, M4(3, 3);
    writeln(vsh, "CalcCoordsTransformFromOldCSToNew starts working");
    if(CoordsInOldCSParam!=NULL){
        MCoordsInOldCS=(*(CoordsInOldCSParam));
    }
    writeln(vsh, "MCoordsInOldCS");
    MCoordsInOldCS.StdCOut2D();
    if(EulerAnglesParam!=NULL){
        MEulerAngles=(*(EulerAnglesParam));
    }
    writeln(vsh, "MEulerAngles");
    MCoordsInOldCS.StdCOut2D();
    if(CSOriginOfNewCSInOldParam!=NULL){
        MCSOriginOfNewCSInOld=(*(CSOriginOfNewCSInOldParam));
    }
    writeln(vsh, "MCSOriginOfNewCSInOld");
    MCSOriginOfNewCSInOld.StdCOut2D();
    MDirCoss=this->CalcMatrixOfDirCossByEulersAngles(MEulerAngles);
    writeln(vsh, "MDirCoss");
    MDirCoss.StdCOut2D();
    M1=MDirCoss.InverseMatrix();
    writeln(vsh, "MDirCoss Inverted");
    M1.StdCOut2D();
    M2=MCoordsInOldCS-MCSOriginOfNewCSInOld;
    writeln(vsh, "M2=MCoordsInOldCS-MCSOriginOfNewCSInOld");
    M2.StdCOut2D();
    writeln(vsh, "M3=M1*M2");
    M3=M1*M2;
    M3.StdCOut2D();
    MCoordsInNewCS=M3;
    writeln(vsh, "CalcCoordsTransformFromOldCSToNew finishes working");
    return MCoordsInNewCS;
    //
    //    Matrix M6, M7 = null;
    //    M1.SetSize(3, 1, 1);
    //    M3.SetSize(3, 1, 1);
    //    M2.SetSize(3, 3, 1);
    //    //M7 = M2.TransposeTo();
    //    M7 = M2.InverseMatrix();
    //    M6 = M1 - M3;
    //    M4 = M7 * M6;
    //
}

Matrix_V2A Matrix_V2A::CalcCoordsTransformFromNewCSToOld(Matrix_V2A*CoordsInNewCSParam, Matrix_V2A*EulerAnglesParam, Matrix_V2A*CSOriginOfNewCSInOldParam, TValsShowHide*vsh){
    Matrix_V2A MCoordsInNewCS(3, 1), MEulerAngles (3, 1), MCSOriginOfNewCSInOld(3, 1);
    Matrix_V2A MDirCoss(3, 3), MCoordsInOldCS(3, 1), M1(3, 3), M2(3, 3), M3(3, 3);//, M4(3, 3);
    if(CoordsInNewCSParam!=NULL){
        MCoordsInNewCS=(*(CoordsInNewCSParam));
    }
    if(EulerAnglesParam!=NULL){
        MEulerAngles=(*(EulerAnglesParam));
    }
    if(CSOriginOfNewCSInOldParam!=NULL){
        MCSOriginOfNewCSInOld=(*(CSOriginOfNewCSInOldParam));
    }
    // Matrix      CalcMatrixOfDirCossByEulersAngles(Matrix MEulerAngles);
    MDirCoss=this->CalcMatrixOfDirCossByEulersAngles(       MEulerAngles);
    M1=MDirCoss;
    M2=M1*MCoordsInNewCS;
    MCoordsInOldCS=M2+MCSOriginOfNewCSInOld;
    return MCoordsInOldCS;
    //
    //    Matrix M6 = null;
    //    M1.SetSize(3, 1, 1);
    //    M3.SetSize(3, 1, 1);
    //    M2.SetSize(3, 3, 1);
    //    M6 = M2 * M1;
    //    M4 = M6 + M3;
    //
}
*/



// }


// class Matrix_PS1{
//     double**y;
//     int QLins, QCols;
//     public:
/*

Matrix_P2S_v1::Matrix_P2S_v1(int QLins, int QCol){
 std::cout<<"Matrix_PS1 constr (QL, QC) works"<<std::endl;
 this->construct();
 this->SetSize(QLins, QCols);

}
Matrix_P2S_v1::Matrix_P2S_v1(double**y, int QLins, int QCols){
 std::cout<<"Matrix_PS1 constr (**dbl, QL, QC) works"<<std::endl;
 this->construct();
 this->Set(y, QLins, QCols);
}
Matrix_P2S_v1::Matrix_P2S_v1(std::vector<std::vector<double>>y){
 std::cout<<"Matrix_PS1 constr (vector<<dbl>>, QL, QC) works"<<std::endl;
 this->construct();
 this->Set(y);
}

Matrix_P2S_v1::Matrix_P2S_v1(double*y, int Q, bool ColOfVectorNotHorLine){
 std::cout<<"Matrix_PS1 constr (*dbl, Q, bool) works"<<std::endl;
 this->construct();
 this->Set(y, Q, ColOfVectorNotHorLine);
}
Matrix_P2S_v1::Matrix_P2S_v1(std::vector<double>y, bool ColOfVectorNotHorLine){
 std::cout<<"Matrix_PS1 constr (vector<dbl>,  bool) works"<<std::endl;
 this->construct();
 this->Set(y, ColOfVectorNotHorLine);
}
Matrix_P2S_v1::Matrix_P2S_v1(double x, double y, double z){
 std::cout<<"Matrix_PS1 constr (x, y, z) works"<<std::endl;
 this->construct();
 this->Set(x, y, z);
}
Matrix_P2S_v1::Matrix_P2S_v1(const Matrix_P2S_v1&obj){
 std::cout<<"Matrix_PS1 COPY constr works"<<std::endl;
 this->construct();
 this->Assign(obj);
}
void Matrix_P2S_v1::construct(){ this->SetNullIni(); }
void Matrix_P2S_v1::SetNull(){
 if(this->y!=NULL){
     for(int i=1; i<=this->QLins; i++){
         delete[] this->y[i-1];
     }
     delete[] this->y;
 }
 this->QLins=0;
 this->QCols=0;
}
void Matrix_P2S_v1::SetNullIni(){ this->y=NULL; this->QLins=0; this->QCols=0; }
Matrix_P2S_v1::~Matrix_P2S_v1(){
 std::cout<<"Matrix_PS1 Destructor works"<<std::endl;
 this->SetNull();
}
void Matrix_P2S_v1::Assign(const Matrix_P2S_v1&obj){
 this->SetNull();
 this->Set(obj.y, obj.QLins, obj.QCols);
}
Matrix_P2S_v1& Matrix_P2S_v1::operator = (const Matrix_P2S_v1&obj){
 this->Assign(obj);
 return*this;
}
void Matrix_P2S_v1::SetSize(int QLins, int QCols){
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
     //if(this->y!=NULL){
     //    for(int i=1; i<=this->QLins; i++){
     //        delete[] this->y[i-1];
     //    }
     //    delete[] this->y;
     //}
     this->SetNull();
     //
     this->y=z;
     //
     this->QLins=QLins;
     this->QCols=QCols;
 }
}

int Matrix_P2S_v1::GetQLines() const { return this->QLins; }
int Matrix_P2S_v1::GetQColumns() const { return this->QCols; }
double Matrix_P2S_v1::GetComponent(int LineN, int ColN) const{
 double y=0;
 if(LineN>=1 && LineN<=this->QLins && ColN>=1 && ColN<=this->QCols){
     y=this->y[LineN-1][ColN-1];
 }
 return y;
}
void Matrix_P2S_v1::SetComponent(double val, int LineN, int ColN){
 if(LineN>=1 && LineN<=this->QLins && ColN>=1 && ColN<=this->QCols){
     this->y[LineN-1][ColN-1]=val;
 }
}

Matrix_P2S_v1 Matrix_P2S_v1::operator *(double x){
    double cur;
    for(int i=1; i<=this->QLins; i++){
        for(int j=1; j<=this->QCols; j++){
           cur=this->y[i-1][j-1]*x;
           this->y[i-1][j-1]=cur;
        }
    }
    return *this;
}


QString Matrix_P2S_v1::LineToString(int LineN, QString delim, bool useBrackets){
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



void Matrix_P2S_v1::Set(std::vector<std::vector<double>>y){
    Matrix_Prototype::Set(y);
    thios-
}
void Matrix_P2S_v1::Set(double*y, int QLines,  int QColumns){
    Matrix_P2S_v1::Set(y, QLins, QCols);
}
void Matrix_P2S_v1::Set(std::vector<double>y, int QColumns){
    Matrix_P2S_v1::Set(y, QColumns);
}
void Matrix_P2S_v1::Set(double x, double y, double z){
    Matrix_P2S_v1::Set(x, y, z);
}
void Matrix_P2S_v1::Set(double y, int QLines, int QColumns){
    Matrix_P2S_v1::Set(y, QLines, QColumns);
}




void Matrix_P2S_v1::StdCOut2D(QString delimElem, QString delimRow, bool useBrackets){
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

*/




//class Matrix{
//    Matrix_Prototype*data;
//public:
    Matrix::Matrix(int QLines, int QColumns, double val, int TypeN){
        switch(TypeN){
            case MatrixTypeN_P2S:

            break;

        }
    }
    Matrix::Matrix(double**y, int QLines, int QColumns, int TypeN){
        this->construct(TypeN);
        this->data->Set(y, QLines, QColumns);
    }
    Matrix::Matrix(double*y, int QLines, int QColumns, int TypeN){
        this->construct(TypeN);
        this->data->Set(y, QLines, QColumns);
    }
    Matrix::Matrix(std::vector<std::vector<double>>y, int QLines, int QColumns, int TypeN){
        this->construct(TypeN);
        this->data->Set(y);
    }
    Matrix::Matrix(std::vector<double>y, int QColumns, int TypeN){
        this->construct(TypeN);
        this->data->Set(y, QColumns);
    }
    //
    Matrix::Matrix(const Matrix&obj){
        this->construct(obj.GetTypeN());
        this->Assign(obj);
    }
    //
    Matrix::~Matrix(){
        delete this->data;
    }

    //
    void Matrix::construct(int TypeN){
        switch(TypeN){
            case MatrixTypeN_P2S:
                this->data=new Matrix_P2S();
            break;
            case MatrixTypeN_P2L:
                this->data=new Matrix_P2L();
            break;
            case MatrixTypeN_P2A:
                this->data=new Matrix_P2A();
            break;
            case MatrixTypeN_V2S:
                this->data=new Matrix_V2S();
            break;
            case MatrixTypeN_V2L:
                this->data=new Matrix_V2L();
            break;
            case MatrixTypeN_V2A:
                this->data=new Matrix_V2A();
            break;
            //rest not written yet
            //case MatrixTypeN_P1S:
            //    this->data=new Matrix_P1S();
            //break;
            //case MatrixTypeN_P1L:
            //    this->data=new Matrix_P1L();
            //break;
            //case MatrixTypeN_P1A:
            //    this->data=new Matrix_P1A();
            //break;
            //case MatrixTypeN_V1S:
            //    this->data=new Matrix_V1S();
            //break;
            //case MatrixTypeN_V1L:
            //    this->data=new Matrix_V1L();
            //break;
            //case MatrixTypeN_V1A:
            //    this->data=new Matrix_V1A();
            //break;
        }
    }//construct

    void Matrix::SetNull(){
        this->data->SetNull();
        delete this->data;
    }

    //void Matrix::SetNullIni();//no need
    void  Matrix::SetOne(double val){
        this->data->SetOne(val);
    }

    //
    void  Matrix::Assign(const Matrix&obj){
        std::vector<std::vector<double>>content2D=obj.GetContentAsVector2DByLines();
        std::vector<double>content1D=obj.GetContentAsVector1DByLines();
        int QLines1, QLines2, QColumns1, QColumns2, TypeN1, TypeN2;
        TypeN1=this->GetTypeN();
        TypeN2=obj.GetTypeN();
        QLines1=this->GetQLines();
        QLines2=obj.GetQLines();
        QColumns1=this->GetQColumns();
        QColumns2=obj.GetQColumns();
        if(QLines1==QLines2 && QColumns1==QColumns2 && TypeN1==TypeN2){
            //NOp;
        }else{
            delete this->data;
        }
        this->Set(content1D, QColumns2);
    }

    Matrix& Matrix::operator = (const Matrix&obj){
        this->Assign(obj);
        return *this;
    }

    //


    void Matrix::Set(int QLines, int QColumns, double val, int TypeN){
        this->data->Set(val, QLines, QColumns);
        if(TypeN!=this->GetTypeN()){
            this->SetTypeN(TypeN);
        }
    }

    void Matrix::Set(double**y, int QLines, int QColumns, int TypeN){
        this->data->Set(y, QLines, QColumns);
        if(TypeN!=this->GetTypeN()){
            this->SetTypeN(TypeN);
        }
    }

    void Matrix::Set(double*y, int QLines, int QColumns, int TypeN){
        this->data->Set(y, QLines, QColumns);
        if(TypeN!=this->GetTypeN()){
            this->SetTypeN(TypeN);
        }
    }

    void Matrix::Set(std::vector<std::vector<double>>y, int TypeN){
        this->data->Set(y);
        if(TypeN!=this->GetTypeN()){
            this->SetTypeN(TypeN);
        }
    }

    void Matrix::Set(std::vector<double>y, int QColumns, int TypeN){
        this->data->Set(y, QColumns);
        if(TypeN!=this->GetTypeN()){
            this->SetTypeN(TypeN);
        }
    }

    //
    int Matrix::GetQLines() const { return this->data->GetQLines(); }
    int Matrix::GetQColumns() const { return this->data->GetQColumns(); }
    //
    void Matrix::SetComponent(double val, int LineN, int ColN){
        this->data->SetComponent(val, LineN, ColN);
    }
    double Matrix::GetComponent(int LineN, int ColN) const {
        return this->data->GetComponent(LineN, ColN);
    }
    //
    void  Matrix::SetSize(int QLines, int QColumns){
        this->data->SetSize(QLines, QColumns);
    }

    void Matrix::SetLineN(int rowN, double*x, int QValues, int whatN1, int QDfltsBefore, bool DefaulsAreZerosNotOwnVals){
        this->data->SetLine(rowN, x, QValues, whatN1, QDfltsBefore, DefaulsAreZerosNotOwnVals);
    }

    void Matrix::SetLineN(int rowN, std::vector<double>x, int whatN1, int QDfltsBefore, bool DefaulsAreZerosNotOwnVals){
        this->data->SetLine(rowN, x, whatN1, QDfltsBefore, DefaulsAreZerosNotOwnVals);
    }
    void Matrix::SetColumnN(int rowN,double*x, int QValues, int whatN1, int QDfltsBefore, bool DefaulsAreZerosNotOwnVals){
        this->data->SetColumn(rowN, x, QValues, whatN1, QDfltsBefore, DefaulsAreZerosNotOwnVals);
    }

    void Matrix::SetColumnN(int rowN,std::vector<double>x, int whatN1, int QDfltsBefore, bool DefaulsAreZerosNotOwnVals){
        this->data->SetColumn(rowN, x, whatN1, QDfltsBefore, DefaulsAreZerosNotOwnVals);
    }

    void Matrix::AddLine(double*x, int QValues, int whatN1, int QDfltsBefore){
        this->data->AddLine(x, QValues, whatN1, QDfltsBefore);
    }

    void Matrix::AddLine(std::vector<double>x, int whatN1, int QDfltsBefore){
        this->data->AddLine(x, whatN1, QDfltsBefore);
    }

    void Matrix::AddColumn(double*x, int QValues, int whatN1, int QDfltsBefore){
        //virtual void AddColumn(double*x, int Q, int whatN1=1, int QDfltValsBefore=0)=0;
        this->data->AddColumn(x, QValues, whatN1, QDfltsBefore);
    }

    void Matrix::AddColumn(std::vector<double>x, int whatN1, int QDfltsBefore){
        this->data->AddColumn(x, whatN1, QDfltsBefore);
    }
    void Matrix::InsLineN(int rowN, double*x, int QValues, int whatN1, int QDfltsBefore){
        this->data->InsLine(rowN, x, QValues, whatN1, QDfltsBefore);
    }

    void Matrix::InsLineN(int rowN, std::vector<double>x, int whatN1, int QDfltsBefore){
        this->data->InsLine(rowN, x, whatN1, QDfltsBefore);
    }

    void Matrix::InsColumnN(int rowN, double*x, int QValues, int whatN1, int QDfltsBefore){
        this->data->InsColumn(rowN, x, QValues, whatN1, QDfltsBefore);
    }

    void Matrix::InsColumnN(int rowN, std::vector<double>x, int whatN1, int QDfltsBefore){
        this->data->InsColumn(rowN, x, whatN1, QDfltsBefore);
    }

    void Matrix::DelLineN(int rowN){
        this->data->DelLine(rowN);
    }

    void Matrix::DelColumnN(int rowN){
        this->data->DelColumn(rowN);
    }
    //
    void Matrix::SetX(double val){
        this->data->SetX(val);
    }
    void Matrix::SetY(double val){
        this->data->SetY(val);
    }
    void Matrix::SetZ(double val){
        this->data->SetZ(val);
    }
    double Matrix::GetX(){
        return this->data->GetX();
    }
    double Matrix::GetY(){
       return this->data->GetY();
    }
    double Matrix::GetZ(){
       return this->data->GetZ();
    }
    //
    void Matrix::SetTypeN(int TypeN){
            Matrix_Prototype*M;
            int QLines=this->GetQLines(), QColumns=this->GetQColumns();
            switch(TypeN){
                case MatrixTypeN_P2S:
                    M=new Matrix_P2S();
                    //this->data=new Matrix_P2S();
                break;
                case MatrixTypeN_P2L:
                    M=new Matrix_P2L();
                break;
                case MatrixTypeN_P2A:
                    M=new Matrix_P2A();
                break;
                case MatrixTypeN_V2S:
                    M=new Matrix_V2S();
                break;
                case MatrixTypeN_V2L:
                    M=new Matrix_V2L();
                break;
                case MatrixTypeN_V2A:
                    M=new Matrix_V2A();
                break;
                //rest not written yet
                //case MatrixTypeN_P1S:
                //    M=new Matrix_P1S();
                //break;
                //case MatrixTypeN_P1L:
                //    M=new Matrix_P1L();
                //break;
                //case MatrixTypeN_P1A:
                //    M=new Matrix_P1A();
                //break;
                //case MatrixTypeN_V1S:
                //    M=new Matrix_V1S();
                //break;
                //case MatrixTypeN_V1L:
                //    M=new Matrix_V1L();
                //break;
                //case MatrixTypeN_V1A:
                //    M=new Matrix_V1A();
                //break;
            }
            //M->SetSize(this->GetQLines(), this->GetQColumns());
            M->Set(this->GetContentAsVector1DByLines(), this->GetQColumns());
            delete this->data;
            this->data=M;
        }

        int Matrix::GetTypeN() const{
            int TypeN=0;
            if(this->data!=NULL)TypeN=this->data->GetTypeN();
            return TypeN;
        }
        //
        std::vector<double> Matrix::GetLineNAsStdVector(int rowN) const{
            return this->data->GetArrayOfLineN(rowN);
        }
        std::vector<double> Matrix::GetColumnNAsStdVector(int rowN) const{
            return this->data->GetArrayOfColN(rowN);
        }
        std::vector<double> Matrix::GetContentAsVector1DByLines() const{
            return this->data->GetContentAsArray1DByLines();
        }

        std::vector<std::vector<double>> Matrix::GetContentAsVector2DByLines() const{
            return this->data->GetContentAsArray2DByLines();
        }
        Matrix Matrix::Get(int Line1N, int Col1N, int Line2N, int Col2N) const{//ne finished!
            Matrix R;
            std::vector<double>row;
            int QLinesIni=this->GetQLines(), QColumnsIni=this->GetQColumns(),
                QLinesFin, QColumnsFin;
            if(Line1N<0) Line1N+=QLinesIni+1;
            if(Line2N<0) Line2N+=QLinesIni+1;
            if(Col1N<0) Col1N+=QColumnsIni+1;
            if(Col2N<0) Col2N+=QColumnsIni+1;
            if(Line2N>=Line1N){
                QLinesFin=Line2N-Line1N;
                if(Col2N>=Col1N){
                    QColumnsFin=Col2N-Col1N;
                    R.SetSize(QLinesFin, QColumnsFin);
                    for(int i=1; i<=QColumnsFin; i++){
                        row.push_back(0);
                    }
                    if(Line1N<0){
                        for(int i=Line1N; i<=0; i++){//works nur if Line1n<0
                            R.AddLine(row);
                        }
                        for(int i=1; i<=Line2N; i++){

                        }
                    }
                }else{

                }
            }else{
                if(Col2N>=Col1N){

                }else{

                }
            }
        }
        Matrix Matrix::GetVectorOfColN(int ColN) const{
            std::vector<double>ColVct=this->GetColumnNAsStdVector(ColN);
            Matrix M(ColVct, 1, this->GetTypeN());
            return M;
        }
        Matrix Matrix::GetVectorOfLineN(int LineN) const{
            std::vector<double>LinVct=this->GetLineNAsStdVector(LineN);
            Matrix M(LinVct, 1, this->GetTypeN());
            return M;
        }

        Matrix Matrix::GetLineN_AsHorisontalMatrix(int LineN) const{
            std::vector<double>LinVct=this->GetLineNAsStdVector(LineN);
            Matrix M(LinVct, this->GetQColumns(), this->GetTypeN());
            return M;
        }
        //
        //void ConcatAdding(const Matrix&M, int QShiftedCols=0, int FromN=1, int QAdded=0);
        //void ConcatStretching(const Matrix&M, int QShiftedLines=0, int FromN=1, int QAdded=0);
        //
        void Matrix::ConcatAdding(const Matrix&M, int QShiftedCols, int FromN, int QAdded){
            //this->data->ConcatAdding(&M.data, QShiftedCols, FromN, QAdded); //so error: no match - qob?
            Matrix_Prototype*MP;
            MP=M.data->clone();
            this->data->ConcatAdding(MP, QShiftedCols, FromN, QAdded);
        }

        void Matrix::ConcatStretching(const Matrix&M, int QShiftedLines, int FromN, int QAdded){
            //this->data->ConcatStretching(&(M.data), QShiftedLines, FromN, QAdded); //so error: no match - qob?
            Matrix_Prototype*MP;
            MP=M.data;
            this->data->ConcatStretching(MP, QShiftedLines, FromN, QAdded);
        }

        //
        void Matrix::Transpose(){
            int QLines=this->GetQLines(), QColumns=this->GetQColumns();
            //Matrix_Prototype*M=new Matrix(QColumns, QLines, 0, this->GetTypeN());
            Matrix_Prototype*M1=new Matrix_V2L(QColumns, QLines);//vikts: 0.0, if 0 ecri no match!
            Matrix M2(QColumns, QLines, 0, this->GetTypeN());
            double val;
            for(int i=1; i<=QLines; i++){
                for(int j=1; j<=QColumns; j++){
                    val=this->GetComponent(i, j);
                    //M->SetComponent(val, j, i);
                    M1->SetComponent(val, j, i);
                    M2.SetComponent(val, j, i);
                }
            }
            delete this->data;
            this->data=M1;
        }

        Matrix Matrix::GetTransposed() const{
            int QLines=this->GetQLines(), QColumns=this->GetQColumns();
            Matrix M(QColumns, QLines, 0, this->GetTypeN());
            double val;
            for(int i=1; i<=QLines; i++){
                for(int j=1; j<=QColumns; j++){
                    val=this->GetComponent(i, j);
                    M.SetComponent(val, j, i);
                }
            }
            return M;
        }


        Matrix Matrix::GetMinor(int LineN, int ColN) const
        {
            Matrix M;
            M=*this;
            // M = this.DelLineTo(LineN);
            // M = M.DelColTo(ColN);
            M.DelLineN(LineN);
            M.DelColumnN(ColN);
            return M;
        }

        void Matrix::LeaveMinor(int LineN, int ColumnN)
        {
             this->DelColumnN(ColumnN);
             this->DelLineN(LineN);
        }

        double Matrix::determinant() const{
            double r=0, d, AlgSpl;
            int RowN;
            int QColumns=this->GetQColumns(), QLines=this->GetQLines();
            if(QColumns != QLines) return 666;
            else{
                if (QLines == 1){
                    r = this->GetComponent(1, 1);
                }else if (QLines == 2){
                    r = this->GetComponent(1, 1) * this->GetComponent(2, 2) - this->GetComponent(1, 2) * this->GetComponent(2, 1);
                }else if (QLines == 3){
                    r = this->GetComponent(1, 1) * this->GetComponent(2, 2) * this->GetComponent(3, 3) + this->GetComponent(1, 3) * this->GetComponent(2, 1) * this->GetComponent(3, 2) + this->GetComponent(3, 1) * this->GetComponent(1, 2) * this->GetComponent(2, 3) -
                        this->GetComponent(1, 3) * this->GetComponent(2, 2) * this->GetComponent(3, 1) - this->GetComponent(1, 1) * this->GetComponent(2, 3) * this->GetComponent(3, 2) - this->GetComponent(3, 3) * this->GetComponent(1, 2) * this->GetComponent(2, 1);
                }else{
                    RowN = 1;
                    r = 0;
                    for(int j = 1; j <= QColumns; j++){
                        r += GetComponent(RowN, j) * AlgSuppl(RowN, j);
                    }
                }
            }
            return r;
        }
        double Matrix::AlgSuppl(int LineN, int ColumnN) const
        {
             double d, r;
             int QLines=this->GetQLines(), QColumns=this->GetQColumns();
             Matrix M;// in |C# wa M = new Matrix();
             if ((LineN < 1) || (LineN > QLines) || (ColumnN < 1) || (ColumnN > QColumns))
             {
                 r = 666;
             }
             else
             {
                 M = *this;
                 M = GetMinor(LineN, ColumnN);
                 d = M.determinant();
                 r = d;
                 if ((LineN + ColumnN)%2 != 0)
                 {
                     r *= (-1);
                 }
             }
             return r;
        }

        Matrix Matrix::GetInverted() const
        {
            Matrix MIR;
            double det = this->determinant();
            MIR = this->GetUnionMatrix();
            if (det != 0) MIR = MIR * (1 / det);
            return MIR;
        }
        void Matrix::Invert()
        {
            Matrix MIR;
            double det = this->determinant();
            MIR = this->GetUnionMatrix();
            if (det != 0) MIR = MIR * (1 / det);
            this->Assign(MIR);
        }

        Matrix Matrix::GetUnionMatrix() const
        {
             int QLines = this->GetQLines(), QColumns = this->GetQColumns();
             int TQLines = QColumns, TQColumns = QLines;
             //Matrix MT = new Matrix();
             //MT = this.TransposeTo();
             Matrix MUR;//in C# = new Matrix();
             //double[][] c;
             double CurAlgSuppl;//, det = this.Determinant();
             // c = new double[TQLines][];
             //for (int i = 1; i <= TQLines; i++)
             //{
             //    c[i - 1] = new double[TQColumns];
             //}
             MUR.SetSize(TQLines, TQColumns);
             //
             for (int i = 1; i <= QLines; i++)
             {
                  for (int j = 1; j <= QColumns; j++)
                  {
                       CurAlgSuppl = this->AlgSuppl(i, j);
                       //c[j - 1][i - 1] = CurAlgSuppl;
                       MUR.SetComponent(CurAlgSuppl, j, i);
                  }
             }
             //
             //MUR = new Matrix(c, TQLines, TQColumns);
             return MUR;
        }

        void Matrix::ToUnion(){
            int QLines = this->GetQLines(), QColumns = this->GetQColumns();
            int TQLines = QColumns, TQColumns = QLines;
            //Matrix MT = new Matrix();
            //MT = this.TransposeTo();
            Matrix MUR;//in C# = new Matrix();
            //double[][] c;
            double CurAlgSuppl;//, det = this.Determinant();
            // c = new double[TQLines][];
            //for (int i = 1; i <= TQLines; i++)
            //{
            //    c[i - 1] = new double[TQColumns];
            //}
            MUR.SetSize(TQLines, TQColumns);
            //
            for (int i = 1; i <= QLines; i++)
            {
                 for (int j = 1; j <= QColumns; j++)
                 {
                      CurAlgSuppl = this->AlgSuppl(i, j);
                      //c[j - 1][i - 1] = CurAlgSuppl;
                      MUR.SetComponent(CurAlgSuppl, j, i);
                 }
            }
            this->Assign(MUR);
        }

        //
        QString Matrix::LineToString(int LineN, QString delim, bool useBrackets){
            QString R;
            R=this->data->LineToString(LineN, delim, useBrackets);
            return R;
        }

        QString Matrix::ColumnToString(int ColN, QString delim, bool useBrackets){
            QString R;
            R=this->data->ColumnToString(ColN, delim, useBrackets);
            return R;
        }

        QString Matrix::ToString(QString delimElem, QString delimRow, bool useBrackets){
            QString R;
            R=this->data->ToString(delimElem, delimRow, useBrackets);
            return R;
        }

        QString Matrix::GetLineAsStr(int LineN, QString delimElem, bool useBrackets){
            QString R;
            //R=this->data->GetLineAsStr(LineN, delimElem, useBrackets);
            R=this->data->LineToString(LineN, delimElem, useBrackets);
            return R;
        }
        QString Matrix::GetColumnAsStr(int ColN, QString delimElem, bool useBrackets){
            QString R;
            //=this->data->GetColumnAsStr(ColN, delimElem, useBrackets);
            R=this->data->ColumnToString(ColN, delimElem, useBrackets);
            return R;
        }

        QString Matrix::GetAllLinesAsStr(QString delimElem, QString delimLines, bool useBrackets){
            QString R;
            //R=this->data->GetAllLinesAsStr(delimElem, delimLines, useBrackets);
            R=this->data->ToString(delimElem, delimLines, useBrackets);
            return R;
        }

        void Matrix::StdCOut1D(QString delimElem, QString delimRow, bool useBrackets){
            this->data->StdCOut1D(delimElem, delimRow, useBrackets);
        }

        void Matrix::StdCOut2D(QString delimElem, QString delimRow, bool useBrackets){
            this->data->StdCOut2D(delimElem, delimRow, useBrackets);
        }

    //};


        Matrix operator +(const Matrix& M1, const Matrix& M2){
            double x3, x1, x2;
            Matrix M3;
            int QL1=M1.GetQLines(), QL2=M2.GetQLines(), QC1=M1.GetQColumns(), QC2=M2.GetQColumns(), QL3, QC3, QM;
            if(QL1==0 || QC1==0 || QL1!=QL2 || QC1!=QC2){
                //NOP;
            }else{
                QL3=QL1;
                QC3=QC2;
                M3.SetSize(QL3, QC3);
                for(int i=1; i<=QL3; i++){
                    for(int j=1; j<=QC3; j++){
                        x1=M1.GetComponent(i, j);
                        x2=M2.GetComponent(i, j);
                        x3=x1+x2;
                        M3.SetComponent(x3, i, j);
                    }
                }
            }
            return M3;
        }


        Matrix operator -(const Matrix& M1, const Matrix& M2){
            double x3, x1, x2;
            Matrix M3;
            int QL1=M1.GetQLines(), QL2=M2.GetQLines(), QC1=M1.GetQColumns(), QC2=M2.GetQColumns(), QL3, QC3, QM;
            if(QL1==0 || QC1==0 || QL1!=QL2 || QC1!=QC2){
                //NOP;
            }else{
                QL3=QL1;
                QC3=QC2;
                M3.SetSize(QL3, QC3);
                for(int i=1; i<=QL3; i++){
                    for(int j=1; j<=QC3; j++){
                        x1=M1.GetComponent(i, j);
                        x2=M2.GetComponent(i, j);
                        x3=x1-x2;
                        M3.SetComponent(x3, i, j);
                    }
                }
            }
            return M3;
        }


    Matrix operator *(const Matrix& M1, const Matrix& M2){
        double x;
        Matrix M3;
        int QL1=M1.GetQLines(), QL2=M2.GetQLines(), QC1=M1.GetQColumns(), QC2=M2.GetQColumns(), QL3, QC3, QM;
        double x1, x2, x3, mp;
        if(QC1==QL2){
            QM=QC1;
            QM=QL2;
            QL3=QL1;
            QC3=QC2;
            M3.SetSize(QL3, QC3);
            for(int i=1; i<=QL3; i++){
                for(int j=1; j<=QC3; j++){
                    x3=0;
                    for(int k=1; k<=QM; k++){
                        x1=M1.GetComponent(i, k);
                        x2=M2.GetComponent(k, j);
                        mp=x1*x2;
                        x3+=mp;
                    }
                    M3.SetComponent(x3, i, j);
                }
            }
        }
        return M3;
    }


    Matrix operator *(double k, const Matrix& M){
        Matrix R;
        double x;
        int QL=M.GetQLines(), QC=M.GetQColumns();
        R.SetSize(QL, QC);
        for(int i=1; i<=QL; i++){
            for(int j=1; j<=QC; j++){
                x=M.GetComponent(i, j);
                x*=k;
                R.SetComponent(x, i, j);
            }
        }
        return R;
    }
    Matrix operator *(const Matrix& M, double k){
        Matrix R;
        double x;
        int QL=M.GetQLines(), QC=M.GetQColumns();
        R.SetSize(QL, QC);
        for(int i=1; i<=QL; i++){
            for(int j=1; j<=QC; j++){
                x=M.GetComponent(i, j);
                x*=k;
                R.SetComponent(x, i, j);
            }
        }
        return R;
    }

    Matrix operator /(const Matrix& M1, const Matrix& M2){
        double x;
        Matrix M3, MR;
        int QL1=M1.GetQLines(), QL2=M2.GetQLines(), QC1=M1.GetQColumns(), QC2=M2.GetQColumns(), QL3, QC3, QM;
        if(QL1==QC1 && QL1==QL2){
            M3=M1.GetUnionMatrix();
            MR=M1*M2;
        }
        return MR;
    }


    Matrix fApproxSolveLinAlgEqsSysSeidel(const Matrix&M, const Matrix&V, Matrix X0, TItersPrecision Precision){
        int QLines=M.GetQLines(), QColumns=M.GetQColumns(), countIters=0;
        double delta,
               det=M.determinant(),
               val, valBef, valAft, val1, CoefErst, sum_bef_N, sum_aft_N, deltaSumQuadr;
        bool contin=true, cond_QIters, cond_delta;
        Matrix X(QLines, 1);
        if(X0.GetQColumns()!=1 || X0.GetQLines()!=M.GetQLines()){
            X0.SetSize(M.GetQLines(), 1);
        }
        //if(fApprNonZero(det, Precision.precision)){

        //}else{
            while(contin){
                countIters++;
                for(int i=1; i<=QLines; i++){
                    CoefErst=1/M.GetComponent(i, i);// CoefErst was not declared in this scope - ob wa det=M.Determinant() - big li alt mic
                    sum_bef_N=0; // sum_bef_N was not declared in this scope - big li alt mic
                    sum_aft_N=0; // sum_aft_N was not declared in this scope - big li alt mic
                    if(i==1){
                        for(int j=i+1; j<=i-1; j++){
                            valAft=M.GetComponent(i, j)*X0.GetComponent(i, 1); // valAft was not declared in this scope - big li alt mic
                            sum_aft_N+=valAft;
                        }
                    }else if(i==QLines){
                        for(int j=1; j<=i-1; j++){
                            valBef=M.GetComponent(i, j)*X0.GetComponent(i, 1); // valBef was not declared in this scope - big li alt mic
                            sum_bef_N+=valBef;
                        }
                    }else{
                        for(int j=i+1; j<=i-1; j++){
                            valAft=M.GetComponent(i, j)*X0.GetComponent(i, 1); // valAft was not declared in this scope - big li alt mic
                            sum_aft_N+=valAft; //no error shown
                        }
                        for(int j=1; j<=i-1; j++){
                            valBef=M.GetComponent(i, j)*X0.GetComponent(i, 1); // valBef was not declared in this scope - big li alt mic
                            sum_bef_N+=valBef; //no error shown
                        }
                    }
                    val1=V.GetComponent(i, 1); // val1 was not declared in this scope - big li alt mic
                    val=CoefErst*(val1-(sum_bef_N+sum_aft_N)); // val was not declared in this scope - big li alt mic
                    X.SetComponent(val, i, 1);
                    deltaSumQuadr=sqrt(val*val-val1*val1); //deltaSumQuadr was not declared in this scope - big li alt mic
                    X0.SetComponent(val, i, 1);//ors X0[i]=X[i]
                }//for
                //
                delta=sqrt(val*val-val1*val1); //delta, val, val1 was not declared in this scope - big li alt mic
                //cond_QIters=(countIters<Precision.MaxQIters);
                //cond_delta=(delta<=Precision.precision);
                contin=Precision.ContinieIterations(delta, countIters);
            }//while
        //}

        return X;
    }



    /*
    Matrix Matrix::SubMatrix(int Line1N, int Col1N, int Line2N, int Col2N){
        Matrix M;
        int QLins=this->GetQLines(), QCols=this->GetQColumns(), QLinsNew=0, QColsNew=0, LineN, ColN;
        double val;
        if(Line1N>=1 && Line1N<=QLins && Line2N>=1 && Line2N<=QLins && Col1N>=1 && Col1N<=QCols && Col2N>=1 && Col2N<=QCols){
            if(Line1N<=Line2N  && Col1N<=Col2N){
                QLinsNew=Line2N-Line1N+1;
                QColsNew=Col2N-Col1N+1;
                M.SetSize(QLinsNew, QColsNew);
                //for(int i=Line1N; i<=Line2N; i++){
                for(int i=1; i<=QLinsNew; i++){
                    //LineN=i-Line1N+1;
                    LineN=i+Line1N-1;
                    for(int j=Col1N; j<=Col2N; j++){
                        ColN=Col1N+j-1;
                        //ColN=j-Col1N+1;
                        val=this->GetComponent(LineN, ColN);
                        //val=this->GetElement(i, j);
                        M.SetComponent(val, i, j);
                        //M.SetElement(val, LineN, ColN);
                    }
                }
            }else if(Line1N<=Line2N  && Col1N>Col2N){
                QLinsNew=Line2N-Line1N+1;
                QColsNew=Col1N-Col2N+1;
                M.SetSize(QLinsNew, QColsNew);
                for(int i=1; i<=QLinsNew; i++){
                    LineN=i+Line1N-1;
                    for(int j=1; j>=QColsNew; j++){
                        ColN=Col2N+j-1;//(changes to Col1N)
                        val=this->GetComponent(LineN, ColN);
                        M.SetComponent(val, i, j);
                    }
                }
            }else if(Line1N>=Line2N  && Col1N<=Col2N){
                QLinsNew=Line1N-Line2N+1;
                QColsNew=Col2N-Col1N+1;
                M.SetSize(QLinsNew, QColsNew);
                for(int i=1; i<=QLinsNew; i++){
                    LineN=Line2N+i-1;//(changes to Line1N)
                    for(int j=1; j>=QColsNew; j++){
                        ColN=Col1N+j-1;
                        val=this->GetComponent(LineN, ColN);
                        M.SetComponent(val, i, j);
                    }
                }
            }else{//both 1 > 2
                QLinsNew=Line1N-Line2N+1;
                QColsNew=Col1N-Col2N+1;
                M.SetSize(QLinsNew, QColsNew);
                for(int i=1; i<=QLinsNew; i++){
                    LineN=Line2N+i-1;//(changes to Line1N)
                    for(int j=1; j>=QColsNew; j++){
                        ColN=Col2N+j-1;//(changes to Col1N)
                        val=this->GetComponent(LineN, ColN);
                        M.SetComponent(val, i, j);
                    }
                }
            }
        }
        return M;
    }
    */






    Matrix GivensRotationTransformMatrix(Matrix*M, int N1, int N2, Matrix*P, TValsShowHide*vsh){
        int QL;
        Matrix T(QL,QL);
        if(P==NULL)P=M;
        QL=P->GetQLines();
        double a11=P->GetComponent(N1, N1), a21=P->GetComponent(N2, N1), c12, s12, denom;
        for(int i=1; i<=QL; i++){
            for(int j=1; j<=QL; j++){
                T.SetComponent(0, i, j);
            }
        }
        for(int i=1; i<=QL; i++){
            T.SetComponent(1, i, i);
        }
        if(a21==0){
            c12=1;
            s12=0;
        }else{
            denom=sqrt(a11*a11+a21*a21);
            c12=a11/denom;
            s12=a21/denom;
        }
        writeln(vsh, " T=" +P->ToString());
        writeln(vsh, " a11=" +FloatToStr(a11) + " a21=" +FloatToStr(a21)+ " c12=" +FloatToStr(c12) + " s12=" + FloatToStr(s12));
        T.SetComponent(c12, N1, N1);
        T.SetComponent(c12, N2, N2);
        T.SetComponent(s12, N1, N2);
        T.SetComponent(-s12, N2, N1);
        return T;
    }
    /*
    void QRDecomposition(Matrix*M, Matrix&Q, Matrix&R,  Matrix*P, TValsShowHide* vsh)
    {
        writeln(vsh, "QRDecomposition starts working");
        if (M == NULL) M = this;
        int QL = M->GetQLines();
        int*ns=new int[QL-1];
        //int countTsInLine;
        int jn, ii, jj;
        Q.SetSize(QL, QL);
        for (int i = 1; i <= QL; i++){
            Q.SetComponent(1, i, i);
        }
        Matrix**Ts=new Matrix*[QL-1];
        for (int i = 1; i <= QL - 1; i++)
        {
            ns[i-1] = QL - 1 - i + 1;
            Ts[i - 1] = new Matrix[ns[i - 1]];
        }
        Matrix*TL=NULL;
        Matrix T;//, R=new Matrix(QL, QL);
        R.Assign(*M);
        writeln(vsh, "M="+M->ToString());
        writeln(vsh, "ini Q=" + Q.ToString());
        writeln(vsh, "ini R=" + R.ToString());
        for (int i = 1; i <= QL - 1; i++)
        {
            //countTsInLine = 0;
            for (int j = i + 1; j <= QL; j++)
            {
                //countTsInLine++; //in AddToVector
                writeln(vsh, "i=" + IntToStr(i)+" j="+IntToStr(j));
                //Matrix_VB GivensRotationTransformMatrix(int N1, int N2, Matrix_VB*M=NULL, TValsShowHide*vsh=NULL);
                T = R.GivensRotationTransformMatrix(i, j, NULL, vsh);
                writeln(vsh, "T=" + T.ToString());
                R = T * R;
                writeln(vsh, "cur R=" + R.ToString());
                jn = j - i;
                Ts[i - 1][jn - 1] = T;
            }
        }
        //R is completed, beginning wit Q
        for(int i=1; i<=QL-1; i++){
            ii = QL - 1 - i + 1;
            for (int j = 1; j <= ns[ii - 1]; j++){
                jj = ns[ii - 1] - j + 1;
                T = Ts[ii - 1][jj - 1];
                Q = Q * T;
            }
        }
        Q.Transpose();
        //
        delete[]ns;
        if(TL!=NULL)delete[]TL;
        for(int i=1; i<=QL-1; i++){
            delete[]Ts[i-1];
        }
        delete[]Ts;
        //
        writeln(vsh, "Answer:");
        writeln(vsh,"Q="+Q.ToString()+" R="+R.ToString());
        writeln(vsh, "QRDecomposition finishes working");
        //
    }//fn QR decomp
    */

    Matrix CalcMatrixOfDirCossByEulersAngles (Matrix MEulerAngles){
        Matrix MR(3, 3);
        double gamma = MEulerAngles.GetX() * M_PI / 180,
                           psi = MEulerAngles.GetY() * M_PI / 180,
                           theta = MEulerAngles.GetZ() * M_PI / 180;
                    MR.SetComponent(cos(psi) * cos(theta), 1, 1);
                    MR.SetComponent(sin(psi) * sin(gamma) - cos(psi) * sin(theta) * cos(gamma), 1, 2);
                    MR.SetComponent(sin(psi) * cos(gamma) + cos(psi) * sin(theta) * sin(gamma), 1, 3);
                    MR.SetComponent(sin(theta), 2, 1);
                    MR.SetComponent(cos(theta) * cos(gamma), 2, 2);
                    MR.SetComponent(-cos(theta) * sin(gamma), 2, 3);
                    MR.SetComponent(-sin(psi) * cos(theta), 3, 1);
                    MR.SetComponent(cos(psi) * sin(gamma) + sin(psi) * sin(theta) * cos(gamma), 3, 2);
                    MR.SetComponent(cos(psi) * cos(gamma) - sin(psi) * sin(theta) * sin(gamma), 3, 3);
                    return MR;
    }

    Matrix CalcCoordsTransformFromOldCSToNew(Matrix*CoordsInOldCSParam, Matrix*EulerAnglesParam, Matrix*CSOriginOfNewCSInOldParam, TValsShowHide*vsh){
        Matrix MCoordsInOldCS(3, 1), MEulerAngles (3, 1),MCSOriginOfNewCSInOld(3, 1);
        Matrix MDirCoss(3, 3), MCoordsInNewCS(3, 1), M1(3, 3), M2(3, 3), M3(3, 3);//, M4(3, 3);
        writeln(vsh, "CalcCoordsTransformFromOldCSToNew starts working");
        if(CoordsInOldCSParam!=NULL){
            MCoordsInOldCS=(*(CoordsInOldCSParam));
        }
        writeln(vsh, "MCoordsInOldCS");
        MCoordsInOldCS.StdCOut2D();
        if(EulerAnglesParam!=NULL){
            MEulerAngles=(*(EulerAnglesParam));
        }
        writeln(vsh, "MEulerAngles");
        MCoordsInOldCS.StdCOut2D();
        if(CSOriginOfNewCSInOldParam!=NULL){
            MCSOriginOfNewCSInOld=(*(CSOriginOfNewCSInOldParam));
        }
        writeln(vsh, "MCSOriginOfNewCSInOld");
        MCSOriginOfNewCSInOld.StdCOut2D();
        MDirCoss=CalcMatrixOfDirCossByEulersAngles(MEulerAngles);
        writeln(vsh, "MDirCoss");
        MDirCoss.StdCOut2D();
        M1=MDirCoss.GetInverted();
        writeln(vsh, "MDirCoss Inverted");
        M1.StdCOut2D();
        M2=MCoordsInOldCS-MCSOriginOfNewCSInOld;
        writeln(vsh, "M2=MCoordsInOldCS-MCSOriginOfNewCSInOld");
        M2.StdCOut2D();
        writeln(vsh, "M3=M1*M2");
        M3=M1*M2;
        M3.StdCOut2D();
        MCoordsInNewCS=M3;
        writeln(vsh, "CalcCoordsTransformFromOldCSToNew finishes working");
        return MCoordsInNewCS;
        /*
            Matrix M6, M7 = null;
            M1.SetSize(3, 1, 1);
            M3.SetSize(3, 1, 1);
            M2.SetSize(3, 3, 1);
            //M7 = M2.TransposeTo();
            M7 = M2.InverseMatrix();
            M6 = M1 - M3;
            M4 = M7 * M6;
        */
    }

    Matrix CalcCoordsTransformFromNewCSToOld(Matrix*CoordsInNewCSParam, Matrix*EulerAnglesParam, Matrix*CSOriginOfNewCSInOldParam, TValsShowHide*vsh){
        Matrix MCoordsInNewCS(3, 1), MEulerAngles (3, 1), MCSOriginOfNewCSInOld(3, 1);
        Matrix MDirCoss(3, 3), MCoordsInOldCS(3, 1), M1(3, 3), M2(3, 3), M3(3, 3);//, M4(3, 3);
        if(CoordsInNewCSParam!=NULL){
            MCoordsInNewCS=(*(CoordsInNewCSParam));
        }
        if(EulerAnglesParam!=NULL){
            MEulerAngles=(*(EulerAnglesParam));
        }
        if(CSOriginOfNewCSInOldParam!=NULL){
            MCSOriginOfNewCSInOld=(*(CSOriginOfNewCSInOldParam));
        }
        // Matrix      CalcMatrixOfDirCossByEulersAngles(Matrix MEulerAngles);
        MDirCoss=CalcMatrixOfDirCossByEulersAngles(       MEulerAngles);//qs'ver?
        M1=MDirCoss;
        M2=M1*MCoordsInNewCS;
        MCoordsInOldCS=M2+MCSOriginOfNewCSInOld;
        return MCoordsInOldCS;
        //
        //    Matrix M6 = null;
        //    M1.SetSize(3, 1, 1);
        //    M3.SetSize(3, 1, 1);
        //    M2.SetSize(3, 3, 1);
        //    M6 = M2 * M1;
        //    M4 = M6 + M3;
        //
    }



    Matrix ConcatAdding(const Matrix&M1, const Matrix&M2, int QShiftedCols, int FromN, int QAdded){
        Matrix MR;
        std::vector<double>row;
        double val;
        int QLines1=M1.GetQLines(),  QLines2=M2.GetQLines(),  QLinesRslt,
            QCols1=M1.GetQColumns(), QCols2=M2.GetQColumns(), QColsRslt;
        if(QAdded==0){
            QAdded=QLines2-FromN+1;
        }
        if(QShiftedCols>0){
            QColsRslt = QCols1 >= QShiftedCols+QCols2 ? QCols1 : QShiftedCols+QCols2;
        }else{
            QColsRslt = -QShiftedCols+QCols2 >= QCols1 ? -QShiftedCols+QCols2 :  QCols1;
        }
        QLinesRslt=QLines1+QAdded;
        if(FromN>=1 && FromN<=QLines2 && QAdded>0){
            MR.SetSize(QLinesRslt, QColsRslt);
            if(QShiftedCols>0){
                for(int i=1; i<=QLines1; i++){
                    for(int j=1; j<=QCols1; j++){
                        val=M1.GetComponent(i, j);
                        MR.SetComponent(i, j);
                    }
                }
                for(int i=FromN; i<=FromN-1+QAdded; i++){
                    for(int j=1; j<=QCols2; j++){
                        val=M2.GetComponent(i, j);
                        MR.SetComponent(val, QLines1+i-FromN+1, QShiftedCols+j);
                    }
                }
            }else{
                for(int i=1; i<=QLines1; i++){
                    for(int j=1; j<=QCols1; j++){
                        val=M1.GetComponent(i, j);
                        MR.SetComponent(i, j+QShiftedCols);
                    }
                }
                for(int i=FromN; i<=FromN-1+QAdded; i++){
                    for(int j=1; j<=QCols2; j++){
                        val=M2.GetComponent(i, j);
                        MR.SetComponent(val, QLines1+i-FromN+1, j);
                    }
                }
            }
        }
        return MR;
    }
    //
    Matrix ConcatStretching(const Matrix&M1, const Matrix&M2, int QShiftedLines, int FromN, int QAdded){
        Matrix MR;
        std::vector<double>row;
        double val;
        int QLines1=M1.GetQLines(),  QLines2=M2.GetQLines(),  QLinesRslt,
            QCols1=M1.GetQColumns(), QCols2=M2.GetQColumns(), QColsRslt;
        if(QAdded==0){
            QAdded=QCols2-FromN+1;
        }
        if(QShiftedLines>0){
            QLinesRslt = QLines1 >= QShiftedLines+QLines2 ? QLines1 : QShiftedLines+QLines2;
        }else{
            QLinesRslt = -QShiftedLines+QLines1 >= QLines2 ? -QShiftedLines+QLines1 : QLines2;
        }
        QColsRslt=QCols1+QAdded;
        if(FromN>=1 && FromN<=QCols2 && QAdded>0){
            MR.SetSize(QLinesRslt, QColsRslt);
            if(QShiftedLines>0){
                for(int i=1; i<=QLines1; i++){
                    for(int j=1; j<=QCols1; j++){
                        val=M1.GetComponent(i, j);
                        MR.SetComponent(i, j);
                    }
                }
                for(int i=1; i<=QLines2; i++){
                    for(int j=FromN; j<=FromN-1+QAdded; j++){
                        val=M2.GetComponent(i, j);
                        MR.SetComponent(val, QShiftedLines+i, QCols1+j);
                    }
                }
            }else{
                for(int i=1; i<=QLines1; i++){
                    for(int j=1; j<=QCols1; j++){
                        val=M1.GetComponent(i, j);
                        MR.SetComponent(-QShiftedLines+i, j);
                    }
                }
                for(int i=1; i<=QLines2; i++){
                    for(int j=FromN; j<=FromN-1+QAdded; j++){
                        val=M2.GetComponent(i, j);
                        MR.SetComponent(val, i, QCols1+j);
                    }
                }
            }
        }
        return MR;
    }






    std::vector<double> CharEqCoefsCoefsByKrylov(const Matrix& M){
        std::vector<double> C;
        int QRows=M.GetQLines(), QColumns=M.GetQColumns();
        Matrix vect(QRows, 1);
        vect.SetX(1);//rest wizeros by default
        for(int i=1; i<=QRows-1; i++){
            vect=M*vect;
        }
        return C;
    }

    bool RowsAreLinearlyDependent(std::vector<double>row1, std::vector<double>row2, double precision){
        bool b=true, contin;
        double k1, kcur;
        int Nc, L1=row1.size(), L2=row2.size();
        if(L1>0 || L1!=L2){
            for(int i=1; i<=L1; i++){
                //bool fApprZero(double x, double precision);
                if(fApprZero(row1[i-1], precision) || fApprZero(row2[i-1], precision)){
                    b=false;
                }
            }
            if(b){
                k1=row1[1-1]/row2[1-1];
                for(int i=2; i<=L1; i++){
                    kcur=row1[i-1]/row2[i-1];
                    //bool fApprNE(double x1, double x2, double precision);
                    if(fApprNE(kcur, k1, precision)){
                        b=false;
                    }
                }
            }
        }
        return b;
    }
