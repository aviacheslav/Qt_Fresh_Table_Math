#include "kynemparam.h"



//class DoFVals{
////public:
    //double displacement, velocity, acceleration;
    DoFVals::DoFVals(double displacement, double velocity, double acceleration){
        this->displacement=displacement;
        this->velocity=velocity;
        this->acceleration=acceleration;
    }
    void DoFVals::set(std::vector<double>v, bool ifEmptyLeaveNotZeros){
        int Q=v.size();
        if(Q>=0){
            this->displacement=v[0];
            if(Q>=1){
                this->velocity=v[1];
                if(Q>=2){
                    this->acceleration=v[2];
                }
            }
        }else{
            if(ifEmptyLeaveNotZeros==false){
                this->displacement=0;
                this->velocity=0;
                this->acceleration=0;
            }
        }
    }
    std::vector<double>DoFVals::get(){
        std::vector<double>R;
        R.push_back(this->displacement);
        R.push_back(this->velocity);
        R.push_back(this->acceleration);
        return R;
    }
//};




    //class KynematicCharacteristics1{
    //public:
       // int DoFsCharNum, MaxDerivN;
             //DoFAndDerivsCharNum,
//KynematicCharacteristics1::KynematicCharacteristics1(){
//    this->construct();
//}
        KynematicCharacteristics1::KynematicCharacteristics1(int X, int Y, int Z, int FiX, int FiY, int FiZ, int MaxDerivN, int QPoints){
            this->Set(X, Y, Z, FiX, FiY, FiZ, MaxDerivN, QPoints);
            //std::cout<<" DoFsCharNum="<<DoFsCharNum<<" DoFsAndDerivsCharNum="<<this->DoFsAndDerivsCharNum<<" QPoints="<<QPoints<<std::endl;
            std::cout<<"Conswtr fin arb"<<std::endl;
        }
        KynematicCharacteristics1::KynematicCharacteristics1(int*DoFsPE, int MaxDerivN, int QPoints){
            this->Set(DoFsPE,  MaxDerivN,  QPoints);
            std::cout<<"Conswtr fin arb"<<std::endl;
        }
        //KynematicCharacteristics1::KynematicCharacteristics1(int DoFsCharNums=0, int MaxDerivN=0, int QPoints=0){//do U see default params? Ce arb, ninty qy, mab ob no dflt vals at
        //    this->DoFsAndDerivsCharNum=DoFsCharNums+64*MaxDerivN;
        //    this->QPoints=QPoints;
        //}//ce work, ma let be ne so

        KynematicCharacteristics1::KynematicCharacteristics1(int DoFsCharNums, int MaxDerivN, int QPoints){//do U see default params? Ce arb, ninty qy, mab ob no dflt vals at
            this->DoFsAndDerivsCharNum=DoFsCharNums+64*MaxDerivN;
            this->QPoints=QPoints;
        }

        //void  KynematicCharacteristics1::construct(){
        //    int X=1, Y=0, Z=0, FiX=0, FiY=0, FiZ=0;
        //    int DoFs[]={0, 0, 0, 0, 0, 0};
        //    DoFs[1-1]=X;
        //    DoFs[2-1]=Y;
        //    DoFs[3-1]=Z;
        //    DoFs[4-1]=FiX;
        //    DoFs[5-1]=FiY;
        //    DoFs[6-1]=FiZ;
        //    int DoFsCharN=0;
        //    for(int i=1; i<=6; i++){
        //        DoFsCharN+=DoFs[i-1];
        //  }
        //    int MaxDerivN=2;
        //    this->DoFsAndDerivsCharNum=DoFsCharN+64*(MaxDerivN-1);
        //    //this->DoFsAndDerivsCharNum=DoFsCharN+NaturalPower(64, (MaxDerivN-1));
        //    this->QPoints=1;
        //}

        int KynematicCharacteristics1::GetQPoints()const{ return this->QPoints; }
        int KynematicCharacteristics1::GetQDoFs()const{
            int DoFs[]={0, 0, 0, 0, 0, 0, 0};
            int order=0, ArrayN=0, DoFInArrayN, MaxDerivN=0, QPoints=0, QDoFs=0, DoFsCharNum=0;
            this->Get(DoFs, MaxDerivN, QPoints, QDoFs, DoFsCharNum);
            return QDoFs;
        }
        int KynematicCharacteristics1::ArrayNByDofNDerivN(int DoFN, int DerivN){
            int DoFs[]={0, 0, 0, 0, 0, 0, 0};
            int order=0, ArrayN=0, DoFInArrayN, MaxDerivN=0, QPoints=0, QDoFs=0, DoFsCharNum=0;
            this->Get(DoFs, MaxDerivN, QPoints, QDoFs, DoFsCharNum);
            if(DoFs[DoFN-1]==1){
                DoFInArrayN=DoFN;
                for(int i=DoFN-1; i>=1; i--){
                    if(DoFs[i-1]==0){
                        DoFInArrayN-=1;
                    }
                }
                ArrayN=QDoFs*DerivN+DoFInArrayN;
            }
            return ArrayN;
        }
        void KynematicCharacteristics1::DoFNAndDerivNByArrayN(int ArrayN, int&DoFN, int&DerivN){
            int DoFs[]={0, 0, 0, 0, 0, 0, 0};
            int order=0, DoFInArrayN, MaxDerivN=0, QPoints=0, QDoFs=0, DoFsCharNum=0;
            int count=0, count1=0;
            bool contin;
            this->Get(DoFs, MaxDerivN, QPoints, QDoFs, DoFsCharNum);
            ArrayN=0;
            //
            DerivN=ArrayN/QDoFs;
            //
            if(DerivN==0){
                DoFInArrayN=ArrayN;
            }else{
                DoFInArrayN = ArrayN % (QDoFs*DerivN);
                if(DoFInArrayN==0){
                    DoFInArrayN=6;
                }
            }
            //
            contin=true;
            while(contin){
                count++;
                if(DoFs[count-1]==1){
                    count1++;
                }
                if(count1==DoFInArrayN){
                    DoFN=count;
                    contin=false;
                }
                if(count>=6){
                    contin=false;
                }
            }
        }
        void KynematicCharacteristics1::SetMaxDerivN(int val){
            int MaxDerivN=this->DoFsAndDerivsCharNum/64;
            int DoFsAndDerivs=this->DoFsAndDerivsCharNum%64;
            this->DoFsAndDerivsCharNum=MaxDerivN+64*val;
        }
        int KynematicCharacteristics1::GetMaxDerivN()const{
            return this->DoFsAndDerivsCharNum/64;//int div
        }
        int KynematicCharacteristics1::GetArrayLengthOfPointKynematicParams()const{
            int QDoFs=this->GetQDoFs(),
                MaxDerovN=this->GetMaxDerivN();
            return ( (this->GetQDoFs()) * (this->GetMaxDerivN()+1) );
        }
        int KynematicCharacteristics1::GetDoFsCharNum(){
            return this->DoFsAndDerivsCharNum%64;
        }

        void KynematicCharacteristics1::Set(int X, int Y, int Z, int FiX, int FiY, int FiZ, int MaxDerivN, int QPoints){
            int DoFs[]={0, 0, 0, 0, 0, 0, 0};
            int DoFsCharNum=0;
            DoFs[1-1]=X;
            DoFs[2-1]=Y;
            DoFs[3-1]=Z;
            DoFs[4-1]=FiX;
            DoFs[5-1]=FiY;
            DoFs[6-1]=FiZ;
            //Set(DoFs, MaxDerivN);
            for(int i=1; i<=6; i++){
                DoFsCharNum+=DoFs[i-1]*NaturalPower(2, i-1);
            }
            this->DoFsAndDerivsCharNum=DoFsCharNum+64*MaxDerivN;
            this->QPoints=QPoints;
        }
        void KynematicCharacteristics1::Set(int*DoFs, int MaxDerivN, int QPoints){
            int DoFsCharNum=0;
            int DoFs1[]={0, 0, 0, 0, 0, 0, 0};
            if(DoFs!=NULL){
                for(int i=1; i<=6; i++){
                   DoFs1[i-1]=DoFs[i-1];
                }
            }
            for(int i=1; i<=6; i++){
                DoFsCharNum+=DoFs1[i-1]*NaturalPower(2, i-1);
            }
            this->DoFsAndDerivsCharNum=DoFsCharNum+64*MaxDerivN;
            this->QPoints=QPoints;
        }

        void KynematicCharacteristics1::Set(int DoFsCharNum, int MaxDerivN, int QPoints){
            this->DoFsAndDerivsCharNum=DoFsCharNum+64*MaxDerivN;
            this->QPoints=QPoints;

        }
        void KynematicCharacteristics1::Set(int DoFsAndDerivsCharNum, int QPoints){
            this->DoFsAndDerivsCharNum=DoFsAndDerivsCharNum;
            this->QPoints=QPoints;
        }

        void KynematicCharacteristics1::Set(std::vector<int>DoFs, int MaxDerivN, int QPoints){
            int DoFsCharNum=0;
            int DoFs1[]={0, 0, 0, 0, 0, 0, 0};
            if(DoFs.size()>0){
                for(int i=1; i<=6; i++){
                   DoFs1[i-1]=DoFs[i-1];
                }
            }
            for(int i=1; i<=6; i++){
                DoFsCharNum+=DoFs1[i-1]*NaturalPower(2, i-1);
            }
            this->DoFsAndDerivsCharNum=DoFsCharNum+64*MaxDerivN;
            this->QPoints=QPoints;
        }


        void KynematicCharacteristics1::Get(int DoFs[], int&MaxDerivN, int&QPoints, int&QDoFs, int &DoFsCharNum) const {
            int order=0;
            //QDoFs=this->DoFsAndDerivsCharNum/64;
            QDoFs=0;
            MaxDerivN=this->DoFsAndDerivsCharNum/64;
            DoFsCharNum=this->DoFsAndDerivsCharNum-64*MaxDerivN;//QDoFs;
            CalcDigits(DoFsCharNum, DoFs, order, 2, NULL);
            for(int i=1; i<=6; i++){
                QDoFs+=DoFs[i-1];
            }
            QPoints=this->QPoints;
            //DoFsAndDerivsCharNum=this->DoFsAndDerivsCharNum;
        }

        void KynematicCharacteristics1::get(int&X, int&Y, int&Z, int&FiX, int&FiY, int&FiZ, int&MaxDerivN, int&QPoints){
            int DoFs[]={0, 0, 0, 0, 0, 0, 0};
            int DoFsCharNum=0, QDoFs=0;
            MaxDerivN=0; QPoints=0;
            this->Get(DoFs, MaxDerivN, QPoints, QDoFs, DoFsCharNum);
            X=DoFs[1-1];
            Y=DoFs[2-1];
            Z=DoFs[3-1];
            FiX=DoFs[4-1];
            FiY=DoFs[5-1];
            FiZ=DoFs[6-1];
        }
        void KynematicCharacteristics1::SetDoFs(int X, int Y, int Z, int FiX, int FiY, int FiZ){
            int DoFs[]={0, 0, 0, 0, 0, 0, 0};
            int DoFsCharNum=0;
            DoFs[1-1]=X;
            DoFs[2-1]=Y;
            DoFs[3-1]=Z;
            DoFs[4-1]=FiX;
            DoFs[5-1]=FiY;
            DoFs[6-1]=FiZ;
            for(int i=1; i<=6; i++){
                DoFsCharNum+=DoFs[i-1]*NaturalPower(2, i);
            }
            int MaxDerivN=this->DoFsAndDerivsCharNum/64;
            this->DoFsAndDerivsCharNum=DoFsCharNum+64*MaxDerivN;
        }
        void KynematicCharacteristics1::SetDoFs(int*DoFs){
            int DoFsArr[]={0, 0, 0, 0, 0, 0, 0};
            int DoFsCharNum=0;
            if(DoFs!=NULL){
                for(int i=1; i<6; i++){
                    DoFsArr[i-1]=DoFs[i-1];
                }
            }
            for(int i=1; i<=6; i++){
                DoFsCharNum+=DoFs[i-1]*NaturalPower(2, i);
            }
            int MaxDerivN=this->DoFsAndDerivsCharNum/64;
            this->DoFsAndDerivsCharNum=DoFsCharNum+64*MaxDerivN;
        }
        void KynematicCharacteristics1::SetDoFs(int DoFsCharNum){
            int maxDerivN=this->DoFsAndDerivsCharNum%64;
            this->DoFsAndDerivsCharNum=DoFsCharNum+64*maxDerivN;
        }
        void KynematicCharacteristics1::SetIfDoFNExists(bool val, int N){
            int DoFs[]={0, 0, 0, 0, 0, 0, 0};
            int MaxDerivN=0, QPoints=0,  QDoFs=0, DoFsCharNum=0;
            this->Get(DoFs, MaxDerivN, QPoints, QDoFs, DoFsCharNum);
            if(N>=1 && N<=6){
                DoFs[N-1]=BoolToInt(val);
           }
            this->SetDoFs(DoFs);
        }
        void KynematicCharacteristics1::SetIfDoFXExists(bool val){
           this->SetIfDoFNExists(val, 1);
        }
        void KynematicCharacteristics1::SetIfDoFYExists(bool val){
            this->SetIfDoFNExists(val, 2);
        }
        void KynematicCharacteristics1::SetIfDoFZExists(bool val){
            this->SetIfDoFNExists(val, 3);
        }
        void KynematicCharacteristics1::SetIfDoFFiXExists(bool val){
            this->SetIfDoFNExists(val, 4);
        }
        void KynematicCharacteristics1::SetIfDoFFiYExists(bool val){
            this->SetIfDoFNExists(val, 5);
        }
        void KynematicCharacteristics1::SetIfDoFFiZExists(bool val){
            this->SetIfDoFNExists(val, 6);
        }
        bool KynematicCharacteristics1::GetIfDoFNExists(int DoFN){
            int DoFs[]={0, 0, 0, 0, 0, 0, 0};
            int order=0, QPoints=0, QDoFs=0, MaxDerivN=0;
            int X=0, Y=0, Z=0, FiX=0, FiY=0, FiZ=0;
            this->get(X, Y, Z, FiX, FiY, FiZ, QPoints, MaxDerivN);
            DoFs[1-1]=X;
            DoFs[2-1]=Y;
            DoFs[3-1]=Z;
            DoFs[4-1]=FiX;
            DoFs[5-1]=FiY;
            DoFs[6-1]=FiZ;
            return (DoFs[DoFN-1]==1);
        }
        bool KynematicCharacteristics1:: GetIfDoFXExists(){
            return this->GetIfDoFNExists(1);
        }
        bool KynematicCharacteristics1:: GetIfDoFYExists(){
            return this->GetIfDoFNExists(2);
        }
        bool KynematicCharacteristics1::GetIfDoFZExists(){
            return this->GetIfDoFNExists(3);
        }
        bool KynematicCharacteristics1::GetIfDoFFiXExists(){
            return this->GetIfDoFNExists(4);
        }
        bool KynematicCharacteristics1::GetIfDoFFiYExists(){
            return this->GetIfDoFNExists(5);
        }
        bool KynematicCharacteristics1::GetIfDoFFiZExists(){
            return this->GetIfDoFNExists(6);
        }

        void KynematicCharacteristics1::AddOrInsPoint(){
            this->QPoints++;
        }
        void KynematicCharacteristics1::DelPoint(){
            this->QPoints--;
        }
        void KynematicCharacteristics1::SetQPoints(int QPoints){
            this->QPoints=QPoints;
        }

        //
        QString KynematicCharacteristics1::ToString(){
            QString s="";
            int DoFs[]={0, 0, 0, 0, 0, 0};
            int QPoints, MaxDerivN, DoFsCharNum, QDoFs=0;
            this->Get(DoFs, MaxDerivN, QPoints, QDoFs, DoFsCharNum);
            s=s+" DoFs: ";
            s=s+" Quantity: ";
            s=s+IntToStr(QDoFs);
            s=s+" Characteristics: ";
            s=s+IntToStr(DoFsCharNum);
            s=s+" X: ";
            s=s+IntToStr(DoFs[1-1]);
            s=s+" Y: ";
            s=s+IntToStr(DoFs[2-1]);
            s=s+" Z: ";
            s=s+IntToStr(DoFs[3-1]);
            s=s+" FiX: ";
            s=s+IntToStr(DoFs[4-1]);
            s=s+" FiY: ";
            s=s+IntToStr(DoFs[5-1]);
            s=s+" FiZ: ";
            s=s+IntToStr(DoFs[6-1]);
            s=s+" MaxDerivN=";
            s=s+IntToStr(MaxDerivN);
            s=s+" QPoints=";
            s=s+IntToStr(QPoints);
            return s;
        }

    //};




    //*****************************************************************************



        //class KynematicParams2
        //{
        //public:
        //    //std::vector<PointKynematicParams>data;
        ///    std::vector<double*>data;
        //    KynematicCharacteristics kynChars;
        //    //

        KynematicParams_vb::KynematicParams_vb(int  *DoFsParam, int MaxDerivN, int QPoints){
            int DoFs[]={0, 0, 0, 0, 0, 0}, arrayL;
            KynematicCharacteristics1 kynChars;
            if(DoFsParam!=NULL){
                for(int i=1; i<=6; i++) DoFs[i-1]=DoFsParam[i-1];
            }
            double*dfltVal=new double;
            kynChars.Set(DoFs, MaxDerivN, QPoints);
            arrayL=kynChars.GetArrayLengthOfPointKynematicParams();
            std::vector<std::vector<double>>data;
            Array2DSetSize(data, QPoints, arrayL, dfltVal);
            this->Set(data, &kynChars);
            delete dfltVal;
        }
        KynematicParams_vb::KynematicParams_vb(KynematicCharacteristics1 *kynCharsParam){
                //int DoFsCharNum=1, QPoints=1, MaxDerivN=2, arrayL;
                int QPoints=1, MaxDerivN=2, arrayL;
                //this->data=NULL;
                if(kynCharsParam!=NULL){
                    this->kynChars=(*(kynCharsParam));
                }//else{
                 //   this->kynChars.DoFsAndDerivsCharNum=DoFsCharNum+64*MaxDerivN;
                //}
                arrayL=this->kynChars.GetArrayLengthOfPointKynematicParams();
                std::vector<std::vector<double>>data;
                double*dfltVal=new double;
                dfltVal=0;
                Array2DSetSize(data, QPoints, arrayL, dfltVal);
                this->Set(data, &kynChars);
                delete dfltVal;
            }

            KynematicParams_vb::KynematicParams_vb(std::vector<std::vector<double>>data,  KynematicCharacteristics1 *kynChars){
                int DoFsCharNum=0, QPoints=1, MaxDerivN=2, arrayL;
                //this->data=NULL;
                if(kynChars!=NULL){
                    this->kynChars=(*(kynChars));
                }else{
                    this->kynChars.DoFsAndDerivsCharNum=DoFsCharNum+64*MaxDerivN;
                }
                arrayL=this->kynChars.GetArrayLengthOfPointKynematicParams();
                this->Set(data, kynChars);
            }
            KynematicParams_vb::KynematicParams_vb(const KynematicParams_vb&obj){
                this->Assign(obj);
            }

            KynematicParams_vb::~KynematicParams_vb(){
                ///int Q=this->kynChars.GetQPoints();
                //if(this->data!=NULL){
                //    for(int i=1; i<=Q; i++){
                //        delete [] this->data[i-1];
                //    }
                //    delete [] this->data;
               // }
            }

            void KynematicParams_vb::Assign(const KynematicParams_vb&obj){
                //in class:
                //std::vector<std::vector<double>>data;
                //KynematicCharacteristics1 kynChars;
                // void Set(std::vector<std::vector<double>>data, KynematicCharacteristics1 *kynChars);
                KynematicCharacteristics1 kynCharsNew=obj.kynChars;
                //this->Set(obj.data, &(obj.kynChars));
                this->Set(obj.data, &kynCharsNew);
            }

           KynematicParams_vb KynematicParams_vb:: operator = (const KynematicParams_vb&obj){
                this->Assign(obj);
                return *this;
           }

            void KynematicParams_vb::DelAllPoints(){
                this->data.clear();
                this->kynChars.SetQPoints(0);
            }

            KynematicCharacteristics1 KynematicParams_vb::GetKynematicCharacteristics(KynematicCharacteristics1 *kynChars){
                KynematicCharacteristics1 kynCharResult;
                if(kynChars!=NULL){
                    kynCharResult=(*(kynChars));
                }else{
                    kynCharResult=this->kynChars;
                }
                 return kynCharResult;
            }





            //
            void KynematicParams_vb::Set(std::vector<std::vector<double>>data, KynematicCharacteristics1 *kynChars){
                double*dfltVal=new double;
                dfltVal=0;
                int arrayL=this->kynChars.GetArrayLengthOfPointKynematicParams();
                int QPointsOld=this->kynChars.GetQPoints(), QPointsNew;
                KynematicCharacteristics1 kynCharsNew=this->GetKynematicCharacteristics(kynChars);
                QPointsNew=kynCharsNew.GetQPoints();
                arrayL=this->kynChars.GetArrayLengthOfPointKynematicParams();
                //if(this->data!=NULL){
                if(this->data.size()>0){
                    //for(int i=1; i<=QPointsOld; i++){
                    //    delete[]this->data[i-1];
                    //}
                    //delete[]this->data;
                    //this->data=NULL;
                    this->data.clear();
                }
                if(QPointsNew>0){
                    if(data.size()==0){
                        Array2DSetSize(this->data, QPointsNew, arrayL, dfltVal);
                    }else{
                        this->data=data;
                    }
                }
                this->kynChars=kynCharsNew;
                //
                delete dfltVal;
            }
            void KynematicParams_vb::SetMaxDerivN(int val){
                //KynematicCharacteristics kynCharActual;
                int QDoFs=this->GetQDoFs(), MaxDerivN=this->GetMaxDerivN(), QPoints=this->GetQPoints(), L=this->kynChars.GetArrayLengthOfPointKynematicParams(), n, DoFsCharNum=this->kynChars.GetDoFsCharNum();
                if(val>=0 && val<=2){
                    if(val>this->kynChars.GetMaxDerivN()){
                        n=QDoFs*(MaxDerivN+1);
                        for(int PointN=1; PointN<=QPoints; PointN++){
                            for(int i=1; i<=val-MaxDerivN; i++){
                                for(int j=1; j<=QDoFs; j++){
                                    Array1DAddElement(this->data[PointN-1],  0.0);
                                }
                            }
                        }
                        this->kynChars.Set(DoFsCharNum, val, QPoints);
                    }else if(val<this->kynChars.GetMaxDerivN()){
                        n=QDoFs*(val-MaxDerivN);
                        for(int PointN=1; PointN<=QPoints; PointN++){
                            for(int i=1; i<=n; i++){
                                Array1DDelElementFromN(this->data[PointN-1], L);
                            }
                        }
                        this->kynChars.Set(DoFsCharNum, val, QPoints);
                    }//else NOp;
                }//else NOp;
            }
            int KynematicParams_vb::GetMaxDerivN(){ return this->kynChars.GetMaxDerivN(); }
            int KynematicParams_vb::GetQDoFs(){ return this->kynChars.GetQDoFs(); }

            void KynematicParams_vb::GetParams(int&X, int&Y, int&Z, int&FiX, int&FiY, int&FiZ, int&MaxDerivN, int&QPoints){
                QPoints=0;
                int DoFsCharNum=0;
                MaxDerivN=0;
                this->kynChars.get(X, Y, Z, FiX, FiY, FiZ, MaxDerivN, QPoints);
            }

            void KynematicParams_vb::GetParams(int DoFs[], int&MaxDerivN, int&QPoints){
                int DoFsCharNum=0, QDoFs;
                MaxDerivN=0;
                QPoints=0;
                this->kynChars.Get(DoFs, MaxDerivN, QPoints, QDoFs, DoFsCharNum);
            }

            void KynematicParams_vb::GetDoFsAndDerivs(int&X, int&Y, int&Z, int&FiX, int&FiY, int&FiZ, int&MaxDerivN){
                int QPoints=0, DoFsCharNum=0;
                MaxDerivN=0;
                this->kynChars.get(X, Y, Z, FiX, FiY, FiZ, MaxDerivN, QPoints);
            }

            void KynematicParams_vb::GetDoFsAndDerivs(int DoFs[], int&MaxDerivN){
                int QPoints=0, DoFsCharNum=0, QDoFs=0;
                MaxDerivN=0;
                this->kynChars.Get(DoFs, MaxDerivN, QPoints, QDoFs, DoFsCharNum);
            }

            void KynematicParams_vb::SetParams(KynematicCharacteristics1*kynCharsNew){
                int DoFsOld[]={0, 0, 0, 0, 0, 0}, DoFsNew[]={0, 0, 0, 0, 0, 0};
                int MaxDerivNOld=0, QPointsOld=0, QDoFsOld=0, DoFsCharNumOld=0,
                    MaxDerivNNew=0, QPointsNew=0, QDoFsNew=0, DoFsCharNumNew=0,
                    LOld, LNew;
                double*dfltVal=new double;
                dfltVal=0;
                this->kynChars.Get(DoFsOld, MaxDerivNOld, QPointsOld, QDoFsOld, DoFsCharNumOld);
                kynCharsNew->Get(DoFsNew, MaxDerivNNew, QPointsNew, QDoFsNew, DoFsCharNumNew);
                this->SetMaxDerivN(MaxDerivNNew);
                if(QPointsOld>0){
                    for(int i=1; i<=6; i++){
                        if(DoFsNew[i-1]==1 && DoFsOld[i-1]==0){
                            this->AddDoFN(i);
                        }else if(DoFsNew[i-1]==0 && DoFsOld[i-1]==1){
                            this->ExcludeDoFN(i);
                        }
                    }
                }else if(QPointsNew>0){ // => this->data = NULL or point to nothing
                    LNew=kynCharsNew->GetArrayLengthOfPointKynematicParams();
                    //this->data=new double*[QPointsNew];
                    //for(int i=1; i<=QPointsNew; i++){
                    //    this->data[i-1]=new double[LNew];
                    //}
                    //for(int i=1; i<=QPointsNew; i++){
                    //    for(int j=1; j<=LNew; j++){
                    //        this->data[i-1][j-1]=0;
                    //    }
                    //}
                    Array2DSetSize(this->data, QPointsNew, LNew, dfltVal);
                }
                //this->kynChars.Set(DoFsNew, MaxDerivN);
                this->kynChars=(*(kynCharsNew));
                delete dfltVal;
            }
            void KynematicParams_vb::SetParams(int X, int Y, int Z, int FiX, int FiY, int FiZ, int MaxDerivN, int QPoints){
                int Xcheck, Ycheck, Zcheck, FiXcheck, FiYcheck, FiZcheck, MaxDerivN_check, QPoints_check;
                KynematicCharacteristics1 kynCharsNew;
                if(QPoints<0)QPoints=this->kynChars.GetQPoints();
                kynCharsNew.Set(X,  Y,  Z,  FiX,  FiY,  FiZ,  MaxDerivN, QPoints);
                kynCharsNew.get(Xcheck, Ycheck, Zcheck, FiXcheck, FiYcheck, FiZcheck, MaxDerivN_check, QPoints_check);
                this->SetParams(&kynCharsNew);
            }
            void KynematicParams_vb::SetParams(int *DoFs, int MaxDerivN, int QPoints){
                KynematicCharacteristics1 kynCharsNew;
                if(QPoints<0){
                    QPoints=this->GetQPoints();
                }
                kynCharsNew.Set(DoFs, MaxDerivN, QPoints);
                this->SetParams(&kynCharsNew);
            }


            void KynematicParams_vb::SetDoFs(int*DoFs){
                int MaxDerivN=this->GetMaxDerivN(), QPoints=this->GetQPoints();
                KynematicCharacteristics1 kynCharsNew;
                kynCharsNew.Set(DoFs, MaxDerivN, QPoints);
                this->SetParams(&kynCharsNew);
            }

            void KynematicParams_vb::SetDoFs(int X, int Y, int Z, int FiX, int FiY, int FiZ){
                int MaxDerivN=this->GetMaxDerivN(), QPoints=this->GetQPoints();
                //KynematicCharacteristics1 kynCharsNew;
                //kynCharsNew.Set(X, Y, Z, FiX, FiY, FiZ, MaxDerivN, QPoints);
                this->SetParams(X, Y, Z, FiX, FiY, FiZ, MaxDerivN, QPoints);
            }

            void KynematicParams_vb::SetDoFs(int DoFsCharNums ){
                int QPoints=this->GetQPoints(), MaxDerivN=this->GetMaxDerivN();
                KynematicCharacteristics1 kynCharsNew;
                kynCharsNew.Set(DoFsCharNums, MaxDerivN, QPoints);
                this->SetParams(&kynCharsNew);
            }
            //void KynematicParams_vb::SetDoFsAndDerivs(int DoFsAndDerivsCharNums){
            //    int QPoints=this->GetQPoints();
            //    KynematicCharacteristics1 kynCharsNew;
            //    kynCharsNew.Set()
            //}

            void KynematicParams_vb::SetDoFsAndDerivs(int*DoFs, int MaxDerivN){
                int QPoints=this->GetQPoints();
                this->SetParams(DoFs, MaxDerivN, QPoints);
            }

            void KynematicParams_vb::SetDoFsAndDerivs(int DoFsCharNums, int MaxDerivN){
                int QPoints=this->GetQPoints();
                KynematicCharacteristics1 kynCharsNew;
                kynCharsNew.Set(DoFsCharNums, MaxDerivN, QPoints);
                this->SetParams(&kynCharsNew);
            }

            void KynematicParams_vb::set_x(double val, int PointN){
                int ArrayN, DoFN=1, DerivN=0, MaxDerivN=this->GetMaxDerivN(), allL=this->GetArrayLengthOfPointKynematicParams();
                if(PointN<0)PointN+=(this->GetQPoints()+1);
                if(PointN>=1 && PointN<=this->GetQPoints()){
                    ArrayN=this->kynChars.ArrayNByDofNDerivN(DoFN, DerivN);
                    this->data[PointN-1][ArrayN-1]=val;
                }
            }
            void KynematicParams_vb::set_y(double val, int PointN){
                int ArrayN, DoFN=2, DerivN=0;
                if(PointN<0)PointN+=(this->GetQPoints()+1);
                if(PointN>=1 && PointN<=this->GetQPoints()){
                    ArrayN=this->kynChars.ArrayNByDofNDerivN(DoFN, DerivN);
                    this->data[PointN-1][ArrayN-1]=val;
                }
            }
            void KynematicParams_vb::set_z(double val, int PointN){
                int ArrayN, DoFN=3, DerivN=0;
                if(PointN<0)PointN+=(this->GetQPoints()+1);
                if(PointN>=1 && PointN<=this->GetQPoints()){
                    ArrayN=this->kynChars.ArrayNByDofNDerivN(DoFN, DerivN);
                    this->data[PointN-1][ArrayN-1]=val;
                }
            }
            void KynematicParams_vb::set_vx(double val, int PointN){
                int ArrayN, DoFN=1, DerivN=1;
                if(PointN<0)PointN+=(this->GetQPoints()+1);
                if(PointN>=1 && PointN<=this->GetQPoints() && DerivN<=this->GetMaxDerivN()){
                    ArrayN=this->kynChars.ArrayNByDofNDerivN(DoFN, DerivN);
                    this->data[PointN-1][ArrayN-1]=val;
                }
            }
            void KynematicParams_vb::set_vy(double val, int PointN){
                int ArrayN, DoFN=2, DerivN=1;
                if(PointN<0)PointN+=(this->GetQPoints()+1);
                if(PointN>=1 && PointN<=this->GetQPoints() && DerivN<=this->GetMaxDerivN()){
                    ArrayN=this->kynChars.ArrayNByDofNDerivN(DoFN, DerivN);
                    this->data[PointN-1][ArrayN-1]=val;
                }
            }
            void KynematicParams_vb::set_vz(double val, int PointN){
                int ArrayN, DoFN=3, DerivN=1;
                if(PointN<0)PointN+=(this->GetQPoints()+1);
                if(PointN>=1 && PointN<=this->GetQPoints() && DerivN<=this->GetMaxDerivN()){
                    ArrayN=this->kynChars.ArrayNByDofNDerivN(DoFN, DerivN);
                    this->data[PointN-1][ArrayN-1]=val;
                }
            }
            void KynematicParams_vb::set_ax(double val, int PointN){
                int ArrayN, DoFN=1, DerivN=2;
                if(PointN<0)PointN+=(this->GetQPoints()+1);
                if(PointN>=1 && PointN<=this->GetQPoints() && DerivN<=this->GetMaxDerivN()){
                    ArrayN=this->kynChars.ArrayNByDofNDerivN(DoFN, DerivN);
                    this->data[PointN-1][ArrayN-1]=val;
                }
            }
            void KynematicParams_vb::set_ay(double val, int PointN){
                int ArrayN, DoFN=2, DerivN=2,  allL=this->GetArrayLengthOfPointKynematicParams();
                if(PointN<0)PointN+=(this->GetQPoints()+1);
                if(PointN>=1 && PointN<=this->GetQPoints() && DerivN<=this->GetMaxDerivN()){
                    ArrayN=this->kynChars.ArrayNByDofNDerivN(DoFN, DerivN);
                    this->data[PointN-1][ArrayN-1]=val;
                }
            }
            void KynematicParams_vb::set_az(double val, int PointN){
                int ArrayN, DoFN=3, DerivN=2,  allL=this->GetArrayLengthOfPointKynematicParams();
                if(PointN<0)PointN+=(this->GetQPoints()+1);
                if(PointN>=1 && PointN<=this->GetQPoints()){
                    ArrayN=this->kynChars.ArrayNByDofNDerivN(DoFN, DerivN);
                    this->data[PointN-1][ArrayN-1]=val;
                }
            }
            void KynematicParams_vb::set_fix(double val, int PointN){
                int ArrayN, DoFN=4, DerivN=0;
                if(PointN<0)PointN+=(this->GetQPoints()+1);
                if(PointN>=1 && PointN<=this->GetQPoints()){
                    ArrayN=this->kynChars.ArrayNByDofNDerivN(DoFN, DerivN);
                    this->data[PointN-1][ArrayN-1]=val;
                }
            }
            void KynematicParams_vb::set_fiy(double val, int PointN){
                int ArrayN, DoFN=5, DerivN=0;
                if(PointN<0)PointN+=(this->GetQPoints()+1);
                if(PointN>=1 && PointN<=this->GetQPoints()){
                    ArrayN=this->kynChars.ArrayNByDofNDerivN(DoFN, DerivN);
                    this->data[PointN-1][ArrayN-1]=val;
                }
            }
            void KynematicParams_vb::set_fiz(double val, int PointN){
                int ArrayN, DoFN=6, DerivN=0;
                if(PointN<0)PointN+=(this->GetQPoints()+1);
                if(PointN>=1 && PointN<=this->GetQPoints()){
                    ArrayN=this->kynChars.ArrayNByDofNDerivN(DoFN, DerivN);
                    this->data[PointN-1][ArrayN-1]=val;
                }
            }
            void KynematicParams_vb::set_wx(double val, int PointN){
                int ArrayN, DoFN=4, DerivN=1;
                if(PointN<0)PointN+=(this->GetQPoints()+1);
                if(PointN>=1 && PointN<=this->GetQPoints() && DerivN<=this->GetMaxDerivN()){
                    ArrayN=this->kynChars.ArrayNByDofNDerivN(DoFN, DerivN);
                    this->data[PointN-1][ArrayN-1]=val;
                }
            }
            void KynematicParams_vb::set_wy(double val, int PointN){
                int ArrayN, DoFN=5, DerivN=1;
                if(PointN<0)PointN+=(this->GetQPoints()+1);
                if(PointN>=1 && PointN<=this->GetQPoints() && DerivN<=this->GetMaxDerivN()){
                    ArrayN=this->kynChars.ArrayNByDofNDerivN(DoFN, DerivN);
                    this->data[PointN-1][ArrayN-1]=val;
                }
            }
            void KynematicParams_vb::set_wz(double val, int PointN){
                int ArrayN, DoFN=6, DerivN=1;
                if(PointN<0)PointN+=(this->GetQPoints()+1);
                if(PointN>=1 && PointN<=this->GetQPoints() && DerivN<=this->GetMaxDerivN()){
                    ArrayN=this->kynChars.ArrayNByDofNDerivN(DoFN, DerivN);
                    this->data[PointN-1][ArrayN-1]=val;
                }
            }
            void KynematicParams_vb::set_ex(double val, int PointN){
                int ArrayN, DoFN=4, DerivN=2;
                if(PointN<0)PointN+=(this->GetQPoints()+1);
                if(PointN>=1 && PointN<=this->GetQPoints() && DerivN<=this->GetMaxDerivN()){
                    ArrayN=this->kynChars.ArrayNByDofNDerivN(DoFN, DerivN);
                    this->data[PointN-1][ArrayN-1]=val;
                }
            }
            void KynematicParams_vb::set_ey(double val, int PointN){
                int ArrayN, DoFN=5, DerivN=2;
                if(PointN<0)PointN+=(this->GetQPoints()+1);
                if(PointN>=1 && PointN<=this->GetQPoints() && DerivN<=this->GetMaxDerivN()){
                    ArrayN=this->kynChars.ArrayNByDofNDerivN(DoFN, DerivN);
                    this->data[PointN-1][ArrayN-1]=val;
                }
            }
            void KynematicParams_vb::set_ez(double val, int PointN){
                int ArrayN, DoFN=6, DerivN=2;
                if(PointN<0)PointN+=(this->GetQPoints()+1);
                if(PointN>=1 && PointN<=this->GetQPoints() && DerivN<=this->GetMaxDerivN()){
                    ArrayN=this->kynChars.ArrayNByDofNDerivN(DoFN, DerivN);
                    this->data[PointN-1][ArrayN-1]=val;
                }
            }
            void KynematicParams_vb::SetIfDoFXExists(int DoFN, bool val){
                if(this->kynChars.GetIfDoFNExists(DoFN)==false && val==true){
                    this->AddDoFN(DoFN);
                }else if(this->kynChars.GetIfDoFNExists(DoFN)==true && val==false){
                    this->ExcludeDoFN(DoFN);
                }
            }
            void KynematicParams_vb::SetIfDoFXExists(bool val){
                int DoFN=1;
                this->SetIfDoFXExists(DoFN, val);
            }
            void KynematicParams_vb::SetIfDoFYExists(bool val){
                int DoFN=2;
                this->SetIfDoFXExists(DoFN, val);
            }
            void KynematicParams_vb::SetIfDoFZExists(bool val){
                int DoFN=3;
                this->SetIfDoFXExists(DoFN, val);
            }
            void KynematicParams_vb::SetIfDoFFiXExists(bool val){
                int DoFN=4;
                this->SetIfDoFXExists(DoFN, val);
            }
            void KynematicParams_vb::SetIfDoFFiYExists(bool val){
                int DoFN=5;
                this->SetIfDoFXExists(DoFN, val);
            }
            void KynematicParams_vb::SetIfDoFFiZExists(bool val){
                int DoFN=6;
                this->SetIfDoFXExists(DoFN, val);
            }
            bool KynematicParams_vb::GetIfDoFNExists(int DoFN){
                return this->kynChars.GetIfDoFNExists(DoFN);
            }
            bool KynematicParams_vb:: GetIfDoFXExists(){
                return this->kynChars.GetIfDoFXExists();
            }
            bool KynematicParams_vb::GetIfDoFYExists(){
                return this->kynChars.GetIfDoFYExists();
            }
            bool KynematicParams_vb::GetIfDoFZExists(){
                return this->kynChars.GetIfDoFZExists();
            }
            bool KynematicParams_vb::GetIfDoFFiXExists(){
                return this->kynChars.GetIfDoFFiXExists();
            }
            bool KynematicParams_vb::GetIfDoFFiYExists(){
                return this->kynChars.GetIfDoFFiYExists();
            }
            bool KynematicParams_vb::GetIfDoFFiZExists(){
                return this->kynChars.GetIfDoFFiZExists();
            }
            void KynematicParams_vb::ExcludeDoFN(int ExcludedDoFN){
                int DoFs[]={0, 0, 0, 0, 0, 0};
                //std::vector<int>Ns;
                int MaxDerivN=this->GetMaxDerivN(), N, QPoints=this->GetQPoints(),
                    L=this->kynChars.GetArrayLengthOfPointKynematicParams();
                this->GetParams(DoFs, MaxDerivN, QPoints);
                if(DoFs[ExcludedDoFN-1]==1){
                    for(int PointN=1; PointN<=QPoints; PointN++){
                        for(int DerivN=MaxDerivN; DerivN>=0; DerivN--){
                            N=this->kynChars.ArrayNByDofNDerivN(ExcludedDoFN, DerivN);
                            Array1DDelElementFromN(this->data[PointN-1], N);
                        }
                    }
                }
                DoFs[ExcludedDoFN-1]=0;
                this->kynChars.Set(DoFs, MaxDerivN, QPoints);
               // this->kynChars.SetMaxDerivN();
            }
            void KynematicParams_vb::ExcludeDoFX(){
                this->ExcludeDoFN(1);
            }
            void KynematicParams_vb::ExcludeDoFY(){
                this->ExcludeDoFN(2);
            }
            void KynematicParams_vb::ExcludeDoFZ(){
                this->ExcludeDoFN(3);
            }
            void KynematicParams_vb::ExcludeDoFFiX(){
                this->ExcludeDoFN(4);
            }
            void KynematicParams_vb::ExcludeDoFFiY(){
                this->ExcludeDoFN(5);
            }
            void KynematicParams_vb::ExcludeDoFFiZ(){
                this->ExcludeDoFN(6);
            }
            void KynematicParams_vb::AddDoFN(int AddedDoFN){
                KynematicCharacteristics1 kynCharsNew;
                int order=0, MaxDerivN=this->kynChars.GetMaxDerivN(), NToInsTo, QPoints=this->kynChars.GetQPoints(), countAdded=0, LOld=this->kynChars.GetArrayLengthOfPointKynematicParams(), L=LOld;
                int DoFsOld[]={0, 0, 0, 0, 0, 0}, DoFsNew[]={0, 0, 0, 0, 0, 0};
                if(AddedDoFN>=1 && AddedDoFN<=6){
                    this->GetParams(DoFsOld, MaxDerivN, QPoints);
                    if(DoFsOld[AddedDoFN-1]==0){
                        for(int i=1; i<=6; i++){
                            DoFsNew[i-1]=DoFsOld[i-1];
                        }
                        DoFsNew[AddedDoFN-1]=1;
                        //void Set(int*DoFs, int MaxDerivN, int QPoints)
                        kynCharsNew.Set(DoFsNew, MaxDerivN, QPoints);
                        if(DoFsOld[AddedDoFN-1]==0){
                            kynCharsNew.Set(DoFsNew, MaxDerivN, QPoints);
                        }
                        //for(int DerivN=MaxDerivN; DerivN>=0; DerivN--){
                        for(int DerivN=0; DerivN<=MaxDerivN; DerivN++){
                            NToInsTo=kynCharsNew.ArrayNByDofNDerivN(AddedDoFN, DerivN);
                            //Ns.push_back(NToInsTo);
                            if(L<NToInsTo){//ce can be nur at last DerivN, et ce Add next
                                if(NToInsTo-L>1){
                                    std::cout<<"error adding DoF "<<AddedDoFN<<" DerivN="<<DerivN<<"N="<<NToInsTo<<" L="<<L<<std::endl;

                                }
                                for(int PointN=1; PointN<=QPoints; PointN++){
                                    Array1DAddElement(this->data[PointN-1],  0.0);
                                }
                                L++;
                            }else{
                                for(int PointN=1; PointN<=QPoints; PointN++){
                                    Array1DInsElementToN(this->data[PointN-1],  NToInsTo, 0.0);
                                }
                                L++;
                            }
                            //countAdded++;
                        }
                    }
                }this->kynChars=kynCharsNew;
            }
            void KynematicParams_vb::AddDoFX(){
                this->AddDoFN(1);
            }
            void KynematicParams_vb::AddDoFY(){
                 this->AddDoFN(2);
            }
            void KynematicParams_vb::AddDoFZ(){
                 this->AddDoFN(3);
            }
            void KynematicParams_vb::AddDoFFiX(){
                 this->AddDoFN(4);
            }
            void KynematicParams_vb::AddDoFFiY(){
                 this->AddDoFN(5);
            }
            void KynematicParams_vb::AddDoFFiZ(){
                 this->AddDoFN(6);
            }

            double KynematicParams_vb::get_x(int PointN){
                double val=0;
                int ArrayN=0;
                int DoFN=1, DerivN=0;
                if(PointN<0)PointN+=(this->GetQPoints()+1);
                //if(this->kynChars.GetIfDoFNExists(DoFN) && PointN>=1 && PointN<=this->data.size() && DerivN<this->kynChars.GetMaxDerivN()){
                bool DoFNExists=this->kynChars.GetIfDoFNExists(DoFN);
                if(this->kynChars.GetIfDoFNExists(DoFN) && PointN>=1 && PointN<=this->kynChars.GetQPoints() && DerivN<=this->GetMaxDerivN()){
                    ArrayN=this->kynChars.ArrayNByDofNDerivN(DoFN, DerivN);
                    val=this->data[PointN-1][ArrayN-1];
                }
                return val;
            }
            double KynematicParams_vb::get_y(int PointN){
                double val=0;
                int ArrayN=0;
                int DoFN=2, DerivN=0;
                if(PointN<0)PointN+=(this->GetQPoints()+1);
                //if(this->kynChars.GetIfDoFNExists(DoFN) && PointN>=1 && PointN<=this->data.size() && DerivN<this->kynChars.GetMaxDerivN()){
                if(this->kynChars.GetIfDoFNExists(DoFN) && PointN>=1 && PointN<=this->kynChars.GetQPoints() && DerivN<=this->kynChars.GetMaxDerivN()){
                    ArrayN=this->kynChars.ArrayNByDofNDerivN(DoFN, DerivN);
                    val=this->data[PointN-1][ArrayN-1];
                }
                return val;
            }
            double KynematicParams_vb::get_z(int PointN){
                double val=0;
                int ArrayN=0;
                int DoFN=3, DerivN=0;
                if(PointN<0)PointN+=(this->GetQPoints()+1);
                //if(this->kynChars.GetIfDoFNExists(DoFN) && PointN>=1 && PointN<=this->data.size() && DerivN<this->kynChars.GetMaxDerivN()){
                if(this->kynChars.GetIfDoFNExists(DoFN) && PointN>=1 && PointN<=this->kynChars.GetQPoints() && DerivN<=this->kynChars.GetMaxDerivN()){
                    ArrayN=this->kynChars.ArrayNByDofNDerivN(DoFN, DerivN);
                    val=this->data[PointN-1][ArrayN-1];
                }
                return val;
            }
            double KynematicParams_vb::get_vx(int PointN){
                double val=0;
                int ArrayN=0;
                int DoFN=1, DerivN=1;
                if(PointN<0)PointN+=(this->GetQPoints()+1);
                //if(this->kynChars.GetIfDoFNExists(DoFN) && PointN>=1 && PointN<=this->data.size() && DerivN<this->kynChars.GetMaxDerivN()){
                if(this->kynChars.GetIfDoFNExists(DoFN) && PointN>=1 && PointN<=this->kynChars.GetQPoints() && DerivN<=this->kynChars.GetMaxDerivN()){
                    ArrayN=this->kynChars.ArrayNByDofNDerivN(DoFN, DerivN);
                    val=this->data[PointN-1][ArrayN-1];
                }
                return val;
            }
            double KynematicParams_vb::get_vy(int PointN){
                double val=0;
                int ArrayN=0;
                int DoFN=2, DerivN=1;
                if(PointN<0)PointN+=(this->GetQPoints()+1);
                //if(this->kynChars.GetIfDoFNExists(DoFN) && PointN>=1 && PointN<=this->data.size() && DerivN<this->kynChars.GetMaxDerivN()){
                if(this->kynChars.GetIfDoFNExists(DoFN) && PointN>=1 && PointN<=this->kynChars.GetQPoints() && DerivN<=this->kynChars.GetMaxDerivN()){
                    ArrayN=this->kynChars.ArrayNByDofNDerivN(DoFN, DerivN);
                    val=this->data[PointN-1][ArrayN-1];
                }
                return val;
            }
            double KynematicParams_vb::get_vz(int PointN){
                double val=0;
                int ArrayN=0;
                int DoFN=3, DerivN=1;
                if(PointN<0)PointN+=(this->GetQPoints()+1);
                //if(this->kynChars.GetIfDoFNExists(DoFN) && PointN>=1 && PointN<=this->data.size() && DerivN<this->kynChars.GetMaxDerivN()){
                if(this->kynChars.GetIfDoFNExists(DoFN) && PointN>=1 && PointN<=this->kynChars.GetQPoints() && DerivN<=this->kynChars.GetMaxDerivN()){
                    ArrayN=this->kynChars.ArrayNByDofNDerivN(DoFN, DerivN);
                    val=this->data[PointN-1][ArrayN-1];
                }
                return val;
            }
            double KynematicParams_vb::get_ax(int PointN){
                double val=0;
                int ArrayN=0;
                int DoFN=1, DerivN=2;
                if(PointN<0)PointN+=(this->GetQPoints()+1);
                //if(this->kynChars.GetIfDoFNExists(DoFN) && PointN>=1 && PointN<=this->data.size() && DerivN<this->kynChars.GetMaxDerivN()){
                if(this->kynChars.GetIfDoFNExists(DoFN) && PointN>=1 && PointN<=this->kynChars.GetQPoints() && DerivN<=this->kynChars.GetMaxDerivN()){
                    ArrayN=this->kynChars.ArrayNByDofNDerivN(DoFN, DerivN);
                    val=this->data[PointN-1][ArrayN-1];
                }
                return val;
            }
            double KynematicParams_vb::get_ay(int PointN){
                double val=0;
                int ArrayN=0;
                int DoFN=2, DerivN=2,  allL=this->GetArrayLengthOfPointKynematicParams(), maxDerivN=this->GetMaxDerivN();
                if(PointN<0)PointN+=(this->GetQPoints()+1);
                //if(this->kynChars.GetIfDoFNExists(DoFN) && PointN>=1 && PointN<=this->data.size() && DerivN<this->kynChars.GetMaxDerivN()){
                if(this->kynChars.GetIfDoFNExists(DoFN) && PointN>=1 && PointN<=this->kynChars.GetQPoints() && DerivN<=this->kynChars.GetMaxDerivN()){
                    ArrayN=this->kynChars.ArrayNByDofNDerivN(DoFN, DerivN);
                    val=this->data[PointN-1][ArrayN-1];
                }
                return val;
            }
            double KynematicParams_vb::get_az(int PointN){
                double val=0;
                int ArrayN=0, DoFN=3, DerivN=2,  allL=this->GetArrayLengthOfPointKynematicParams(), maxDerivN=this->GetMaxDerivN();
                if(PointN<0)PointN+=(this->GetQPoints()+1);
                //if(this->kynChars.GetIfDoFNExists(DoFN) && PointN>=1 && PointN<=this->data.size() && DerivN<this->kynChars.GetMaxDerivN()){
                if(this->kynChars.GetIfDoFNExists(DoFN) && PointN>=1 && PointN<=this->kynChars.GetQPoints() && DerivN<=this->kynChars.GetMaxDerivN()){
                    ArrayN=this->kynChars.ArrayNByDofNDerivN(DoFN, DerivN);
                    val=this->data[PointN-1][ArrayN-1];
                }
                return val;
            }
            double KynematicParams_vb::get_fiX(int PointN){
                double val=0;
                int ArrayN=0;
                int DoFN=4, DerivN=0;
                if(PointN<0)PointN+=(this->GetQPoints()+1);
                //if(this->kynChars.GetIfDoFNExists(DoFN) && PointN>=1 && PointN<=this->data.size() && DerivN<this->kynChars.GetMaxDerivN()){
                if(this->kynChars.GetIfDoFNExists(DoFN) && PointN>=1 && PointN<=this->kynChars.GetQPoints() && DerivN<=this->kynChars.GetMaxDerivN()){
                    ArrayN=this->kynChars.ArrayNByDofNDerivN(DoFN, DerivN);
                    val=this->data[PointN-1][ArrayN-1];
                }
                return val;
            }
            double KynematicParams_vb::get_fiY(int PointN){
                double val=0;
                int ArrayN=0;
                int DoFN=5, DerivN=0;
                if(PointN<0)PointN+=(this->GetQPoints()+1);
                //if(this->kynChars.GetIfDoFNExists(DoFN) && PointN>=1 && PointN<=this->data.size() && DerivN<this->kynChars.GetMaxDerivN()){
                if(this->kynChars.GetIfDoFNExists(DoFN) && PointN>=1 && PointN<=this->kynChars.GetQPoints() && DerivN<=this->kynChars.GetMaxDerivN()){
                    ArrayN=this->kynChars.ArrayNByDofNDerivN(DoFN, DerivN);
                    val=this->data[PointN-1][ArrayN-1];
                }
                return val;
            }
            double KynematicParams_vb::get_fiZ(int PointN){
                double val=0;
                int ArrayN=0;
                int DoFN=6, DerivN=0;
                if(PointN<0)PointN+=(this->GetQPoints()+1);
                //if(this->kynChars.GetIfDoFNExists(DoFN) && PointN>=1 && PointN<=this->data.size() && DerivN<this->kynChars.GetMaxDerivN()){
                if(this->kynChars.GetIfDoFNExists(DoFN) && PointN>=1 && PointN<=this->kynChars.GetQPoints() && DerivN<=this->kynChars.GetMaxDerivN()){
                    ArrayN=this->kynChars.ArrayNByDofNDerivN(DoFN, DerivN);
                    val=this->data[PointN-1][ArrayN-1];
                }
                return val;
            }
            double KynematicParams_vb::get_wx(int PointN){
                double val=0;
                int ArrayN=0;
                int DoFN=4, DerivN=1;
                if(PointN<0)PointN+=(this->GetQPoints()+1);
                //if(this->kynChars.GetIfDoFNExists(DoFN) && PointN>=1 && PointN<=this->data.size() && DerivN<this->kynChars.GetMaxDerivN()){
                if(this->kynChars.GetIfDoFNExists(DoFN) && PointN>=1 && PointN<=this->kynChars.GetQPoints() && DerivN<=this->kynChars.GetMaxDerivN()){
                    ArrayN=this->kynChars.ArrayNByDofNDerivN(DoFN, DerivN);
                    val=this->data[PointN-1][ArrayN-1];
                }
                return val;
            }
            double KynematicParams_vb::get_wy(int PointN){
                double val=0;
                int ArrayN=0;
                int DoFN=5, DerivN=1;
                if(PointN<0)PointN+=(this->GetQPoints()+1);
                //if(this->kynChars.GetIfDoFNExists(DoFN) && PointN>=1 && PointN<=this->data.size() && DerivN<this->kynChars.GetMaxDerivN()){
                if(this->kynChars.GetIfDoFNExists(DoFN) && PointN>=1 && PointN<=this->kynChars.GetQPoints() && DerivN<=this->kynChars.GetMaxDerivN()){
                    ArrayN=this->kynChars.ArrayNByDofNDerivN(DoFN, DerivN);
                    val=this->data[PointN-1][ArrayN-1];
                }
                return val;
            }
            double KynematicParams_vb::get_wz(int PointN){
                double val=0;
                int ArrayN=0;
                int DoFN=6, DerivN=1;
                if(PointN<0)PointN+=(this->GetQPoints()+1);
                //if(this->kynChars.GetIfDoFNExists(DoFN) && PointN>=1 && PointN<=this->data.size() && DerivN<this->kynChars.GetMaxDerivN()){
                if(this->kynChars.GetIfDoFNExists(DoFN) && PointN>=1 && PointN<=this->kynChars.GetQPoints() && DerivN<=this->kynChars.GetMaxDerivN()){
                    ArrayN=this->kynChars.ArrayNByDofNDerivN(DoFN, DerivN);
                    val=this->data[PointN-1][ArrayN-1];
                }
                return val;
            }
            double KynematicParams_vb::get_ex(int PointN){
                double val=0;
                int ArrayN=0;
                int DoFN=4, DerivN=2;
                if(PointN<0)PointN+=(this->GetQPoints()+1);
                //if(this->kynChars.GetIfDoFNExists(DoFN) && PointN>=1 && PointN<=this->data.size() && DerivN<this->kynChars.GetMaxDerivN()){
                if(this->kynChars.GetIfDoFNExists(DoFN) && PointN>=1 && PointN<=this->kynChars.GetQPoints() && DerivN<=this->kynChars.GetMaxDerivN()){
                    ArrayN=this->kynChars.ArrayNByDofNDerivN(DoFN, DerivN);
                    val=this->data[PointN-1][ArrayN-1];
                }
                return val;
            }
            double KynematicParams_vb::get_ey(int PointN){
                double val=0;
                int ArrayN=0;
                int DoFN=5, DerivN=2;
                if(PointN<0)PointN+=(this->GetQPoints()+1);
                //if(this->kynChars.GetIfDoFNExists(DoFN) && PointN>=1 && PointN<=this->data.size() && DerivN<this->kynChars.GetMaxDerivN()){
                if(this->kynChars.GetIfDoFNExists(DoFN) && PointN>=1 && PointN<=this->kynChars.GetQPoints() && DerivN<=this->kynChars.GetMaxDerivN()){
                    ArrayN=this->kynChars.ArrayNByDofNDerivN(DoFN, DerivN);
                    val=this->data[PointN-1][ArrayN-1];
                }
                return val;
            }
            double KynematicParams_vb::get_ez(int PointN){
                double val=0;
                int ArrayN=0;
                int DoFN=6, DerivN=2;
                if(PointN<0)PointN+=(this->GetQPoints()+1);
                //if(this->kynChars.GetIfDoFNExists(DoFN) && PointN>=1 && PointN<=this->data.size() && DerivN<this->kynChars.GetMaxDerivN()){
                if(this->kynChars.GetIfDoFNExists(DoFN) && PointN>=1 && PointN<=this->kynChars.GetQPoints() && DerivN<=this->kynChars.GetMaxDerivN()){
                    ArrayN=this->kynChars.ArrayNByDofNDerivN(DoFN, DerivN);
                    val=this->data[PointN-1][ArrayN-1];
                }
                return val;
            }
            //
            std::vector<double>KynematicParams_vb::Get_all_X(){
                std::vector<double>data;
                double val;
                int QDoFs=this->kynChars.GetQDoFs(),
                    MaxDerivN=this->kynChars.GetMaxDerivN(),
                    DoFN=1, DerivN=0,
                    //QPoints=this->data.size(),
                    QPoints=this->kynChars.GetQPoints(),
                    ArrayN =this->kynChars.ArrayNByDofNDerivN(DoFN, DerivN);
                if(this->kynChars.GetIfDoFNExists(DoFN) && DerivN<=MaxDerivN){
                    for(int PointN=1; PointN<=QPoints; PointN++){
                        val=this->data[PointN-1][ArrayN-1];
                        data.push_back(val);
                    }
                }
                return data;
            }
            std::vector<double>KynematicParams_vb::Get_all_Y(){
                std::vector<double>data;
                double val;
                int QDoFs=this->kynChars.GetQDoFs(),
                    MaxDerivN=this->kynChars.GetMaxDerivN(),
                    DoFN=2, DerivN=0,
                        //QPoints=this->data.size(),
                        QPoints=this->kynChars.GetQPoints(),
                    ArrayN =this->kynChars.ArrayNByDofNDerivN(DoFN, DerivN);
                if(this->kynChars.GetIfDoFNExists(DoFN) && DerivN<=MaxDerivN){
                    for(int PointN=1; PointN<=QPoints; PointN++){
                        val=this->data[PointN-1][ArrayN-1];
                        data.push_back(val);
                    }
                }
                return data;
            }
            std::vector<double>KynematicParams_vb::Get_all_Z(){
                std::vector<double>data;
                double val;
                int QDoFs=this->kynChars.GetQDoFs(),
                    MaxDerivN=this->kynChars.GetMaxDerivN(),
                    DoFN=3, DerivN=0,
                    //QPoints=this->data.size(),
                    QPoints=this->kynChars.GetQPoints(),
                    ArrayN =this->kynChars.ArrayNByDofNDerivN(DoFN, DerivN);
                if(this->kynChars.GetIfDoFNExists(DoFN) && DerivN<=MaxDerivN){
                    for(int PointN=1; PointN<=QPoints; PointN++){
                        val=this->data[PointN-1][ArrayN-1];
                        data.push_back(val);
                    }
                }
                return data;
            }
            std::vector<double>KynematicParams_vb::Get_all_Vx(){
                std::vector<double>data;
                double val;
                int QDoFs=this->kynChars.GetQDoFs(),
                    MaxDerivN=this->kynChars.GetMaxDerivN(),
                    DoFN=1, DerivN=1,
                    //QPoints=this->data.size(),
                    QPoints=this->kynChars.GetQPoints(),
                    ArrayN =this->kynChars.ArrayNByDofNDerivN(DoFN, DerivN);
                if(this->kynChars.GetIfDoFNExists(DoFN) && DerivN<=MaxDerivN){
                    for(int PointN=1; PointN<=QPoints; PointN++){
                        val=this->data[PointN-1][ArrayN-1];
                        data.push_back(val);
                    }
                }
                return data;
            }
            std::vector<double>KynematicParams_vb::Get_all_Vy(){
                std::vector<double>data;
                double val;
                int QDoFs=this->kynChars.GetQDoFs(),
                    MaxDerivN=this->kynChars.GetMaxDerivN(),
                    DoFN=2, DerivN=1,
                    //QPoints=this->data.size(),
                    QPoints=this->kynChars.GetQPoints(),
                    ArrayN =this->kynChars.ArrayNByDofNDerivN(DoFN, DerivN);
                if(this->kynChars.GetIfDoFNExists(DoFN) && DerivN<=MaxDerivN){
                    for(int PointN=1; PointN<=QPoints; PointN++){
                        val=this->data[PointN-1][ArrayN-1];
                        data.push_back(val);
                    }
                }
                return data;
            }
            std::vector<double>KynematicParams_vb::Get_all_Vz(){
                std::vector<double>data;
                double val;
                int QDoFs=this->kynChars.GetQDoFs(),
                    MaxDerivN=this->kynChars.GetMaxDerivN(),
                    DoFN=3, DerivN=1,
                    //QPoints=this->data.size(),
                    QPoints=this->kynChars.GetQPoints(),
                    ArrayN =this->kynChars.ArrayNByDofNDerivN(DoFN, DerivN);
                if(this->kynChars.GetIfDoFNExists(DoFN) && DerivN<=MaxDerivN){
                    for(int PointN=1; PointN<=QPoints; PointN++){
                        val=this->data[PointN-1][ArrayN-1];
                        data.push_back(val);
                    }
                }
                return data;
            }
            std::vector<double>KynematicParams_vb::Get_all_ax(){
                std::vector<double>data;
                double val;
                int QDoFs=this->kynChars.GetQDoFs(),
                    MaxDerivN=this->kynChars.GetMaxDerivN(),
                    DoFN=1, DerivN=2,
                    //QPoints=this->data.size(),
                    QPoints=this->kynChars.GetQPoints(),
                    ArrayN =this->kynChars.ArrayNByDofNDerivN(DoFN, DerivN);
                if(this->kynChars.GetIfDoFNExists(DoFN) && DerivN<=MaxDerivN){
                    for(int PointN=1; PointN<=QPoints; PointN++){
                        val=this->data[PointN-1][ArrayN-1];
                        data.push_back(val);
                    }
                }
                return data;
            }
            std::vector<double>KynematicParams_vb::Get_all_ay(){
                std::vector<double>data;
                double val;
                int QDoFs=this->kynChars.GetQDoFs(),
                    MaxDerivN=this->kynChars.GetMaxDerivN(),
                    DoFN=2, DerivN=2,
                    //QPoints=this->data.size(),
                    QPoints=this->kynChars.GetQPoints(),
                    ArrayN =this->kynChars.ArrayNByDofNDerivN(DoFN, DerivN);
                if(this->kynChars.GetIfDoFNExists(DoFN) && DerivN<=MaxDerivN){
                    for(int PointN=1; PointN<=QPoints; PointN++){
                        val=this->data[PointN-1][ArrayN-1];
                        data.push_back(val);
                    }
                }
                return data;
            }
            std::vector<double>KynematicParams_vb::Get_all_az(){
                std::vector<double>data;
                double val;
                int QDoFs=this->kynChars.GetQDoFs(),
                    MaxDerivN=this->kynChars.GetMaxDerivN(),
                    DoFN=3, DerivN=2,
                    //QPoints=this->data.size(),
                    QPoints=this->kynChars.GetQPoints(),
                    ArrayN =this->kynChars.ArrayNByDofNDerivN(DoFN, DerivN);
                if(this->kynChars.GetIfDoFNExists(DoFN) && DerivN<=MaxDerivN){
                    for(int PointN=1; PointN<=QPoints; PointN++){
                        val=this->data[PointN-1][ArrayN-1];
                        data.push_back(val);
                    }
                }
                return data;
            }
            std::vector<double>KynematicParams_vb::Get_all_FiX(){
                std::vector<double>data;
                double val;
                int QDoFs=this->kynChars.GetQDoFs(),
                    MaxDerivN=this->kynChars.GetMaxDerivN(),
                    DoFN=4, DerivN=0,
                    //QPoints=this->data.size(),
                    QPoints=this->kynChars.GetQPoints(),
                    ArrayN =this->kynChars.ArrayNByDofNDerivN(DoFN, DerivN);
                if(this->kynChars.GetIfDoFNExists(DoFN) && DerivN<=MaxDerivN){
                    for(int PointN=1; PointN<=QPoints; PointN++){
                        val=this->data[PointN-1][ArrayN-1];
                        data.push_back(val);
                    }
                }
                return data;
            }
            std::vector<double>KynematicParams_vb::Get_all_FiY(){
                std::vector<double>data;
                double val;
                int QDoFs=this->kynChars.GetQDoFs(),
                    MaxDerivN=this->kynChars.GetMaxDerivN(),
                    DoFN=5, DerivN=0,
                    //QPoints=this->data.size(),
                    QPoints=this->kynChars.GetQPoints(),
                    ArrayN =this->kynChars.ArrayNByDofNDerivN(DoFN, DerivN);
                if(this->kynChars.GetIfDoFNExists(DoFN) && DerivN<=MaxDerivN){
                    for(int PointN=1; PointN<=QPoints; PointN++){
                        val=this->data[PointN-1][ArrayN-1];
                        data.push_back(val);
                    }
                }
                return data;
            }
            std::vector<double>KynematicParams_vb::Get_all_FiZ(){
                std::vector<double>data;
                double val;
                int QDoFs=this->kynChars.GetQDoFs(),
                    MaxDerivN=this->kynChars.GetMaxDerivN(),
                    DoFN=6, DerivN=0,
                    //QPoints=this->data.size(),
                    QPoints=this->kynChars.GetQPoints(),
                    ArrayN =this->kynChars.ArrayNByDofNDerivN(DoFN, DerivN);
                if(this->kynChars.GetIfDoFNExists(DoFN) && DerivN<=MaxDerivN){
                    for(int PointN=1; PointN<=QPoints; PointN++){
                        val=this->data[PointN-1][ArrayN-1];
                        data.push_back(val);
                    }
                }
                return data;
            }
            std::vector<double>KynematicParams_vb::Get_all_wx(){
                std::vector<double>data;
                double val;
                int QDoFs=this->kynChars.GetQDoFs(),
                    MaxDerivN=this->kynChars.GetMaxDerivN(),
                    DoFN=4, DerivN=1,
                    //QPoints=this->data.size(),
                    QPoints=this->kynChars.GetQPoints(),
                    ArrayN =this->kynChars.ArrayNByDofNDerivN(DoFN, DerivN);
                if(this->kynChars.GetIfDoFNExists(DoFN) && DerivN<=MaxDerivN){
                    for(int PointN=1; PointN<=QPoints; PointN++){
                        val=this->data[PointN-1][ArrayN-1];
                        data.push_back(val);
                    }
                }
                return data;
            }
            std::vector<double>KynematicParams_vb::Get_all_wy(){
                std::vector<double>data;
                double val;
                int QDoFs=this->kynChars.GetQDoFs(),
                    MaxDerivN=this->kynChars.GetMaxDerivN(),
                    DoFN=5, DerivN=1,
                    //QPoints=this->data.size(),
                    QPoints=this->kynChars.GetQPoints(),
                    ArrayN =this->kynChars.ArrayNByDofNDerivN(DoFN, DerivN);
                if(this->kynChars.GetIfDoFNExists(DoFN) && DerivN<=MaxDerivN){
                    for(int PointN=1; PointN<=QPoints; PointN++){
                        val=this->data[PointN-1][ArrayN-1];
                        data.push_back(val);
                    }
                }
                return data;
            }
            std::vector<double>KynematicParams_vb::Get_all_wz(){
                std::vector<double>data;
                double val;
                int QDoFs=this->kynChars.GetQDoFs(),
                    MaxDerivN=this->kynChars.GetMaxDerivN(),
                    DoFN=6, DerivN=1,
                    //QPoints=this->data.size(),
                    QPoints=this->kynChars.GetQPoints(),
                    ArrayN =this->kynChars.ArrayNByDofNDerivN(DoFN, DerivN);
                if(this->kynChars.GetIfDoFNExists(DoFN) && DerivN<=MaxDerivN){
                    for(int PointN=1; PointN<=QPoints; PointN++){
                        val=this->data[PointN-1][ArrayN-1];
                        data.push_back(val);
                    }
                }
                return data;
            }
            std::vector<double>KynematicParams_vb::Get_all_ex(){
                std::vector<double>data;
                double val;
                int QDoFs=this->kynChars.GetQDoFs(),
                    MaxDerivN=this->kynChars.GetMaxDerivN(),
                    DoFN=4, DerivN=2,
                    //QPoints=this->data.size(),
                    QPoints=this->kynChars.GetQPoints(),
                    ArrayN =this->kynChars.ArrayNByDofNDerivN(DoFN, DerivN);
                if(this->kynChars.GetIfDoFNExists(DoFN) && DerivN<=MaxDerivN){
                    for(int PointN=1; PointN<=QPoints; PointN++){
                        val=this->data[PointN-1][ArrayN-1];
                        data.push_back(val);
                    }
                }
                return data;
            }
            std::vector<double>KynematicParams_vb::Get_all_ey(){
                std::vector<double>data;
                double val;
                int QDoFs=this->kynChars.GetQDoFs(),
                    MaxDerivN=this->kynChars.GetMaxDerivN(),
                    DoFN=5, DerivN=2,
                    //QPoints=this->data.size(),
                    QPoints=this->kynChars.GetQPoints(),
                    ArrayN =this->kynChars.ArrayNByDofNDerivN(DoFN, DerivN);
                if(this->kynChars.GetIfDoFNExists(DoFN) && DerivN<=MaxDerivN){
                    for(int PointN=1; PointN<=QPoints; PointN++){
                        val=this->data[PointN-1][ArrayN-1];
                        data.push_back(val);
                    }
                }
                return data;
            }
            std::vector<double>KynematicParams_vb::Get_all_ez(){
                std::vector<double>data;
                double val;
                int QDoFs=this->kynChars.GetQDoFs(),
                    MaxDerivN=this->kynChars.GetMaxDerivN(),
                    DoFN=6, DerivN=2,
                    //QPoints=this->data.size(),
                    QPoints=this->kynChars.GetQPoints(),
                    ArrayN =this->kynChars.ArrayNByDofNDerivN(DoFN, DerivN);
                if(this->kynChars.GetIfDoFNExists(DoFN) && DerivN<=MaxDerivN){
                    for(int PointN=1; PointN<=QPoints; PointN++){
                        val=this->data[PointN-1][ArrayN-1];
                        data.push_back(val);
                    }
                }
                return data;
            }
            //
            void KynematicParams_vb::Set_all_X( std::vector<double>data ){
                int //QPoints=this->data.size(),
                    QPoints=this->kynChars.GetQPoints(),
                    QPointsAdded=data.size(),
                    Qmin=QPoints<=QPointsAdded?QPoints:QPointsAdded,
                    DoFN=1, DerivN=0,
                    MaxDerivN=this->kynChars.GetMaxDerivN(),
                    ArrayN=this->kynChars.ArrayNByDofNDerivN(DoFN, DerivN);
                //if(this->kynChars.GetIfDoFXExists() && this->data.size()>0 ){
                if(this->kynChars.GetIfDoFXExists() && this->kynChars.GetQPoints()>0 ){
                    for(int PointN=1; PointN<=Qmin; PointN++){
                        this->data[PointN-1][ArrayN-1]=data[PointN-1];
                    }
                    for(int PointN=Qmin+1; PointN<=QPoints; PointN++){
                        this->data[PointN-1][ArrayN-1]=data[PointN-1];
                    }
                }
            }
            void KynematicParams_vb::Set_all_Y( std::vector<double>data ){
                int //QPoints=this->data.size(),
                    QPoints=this->kynChars.GetQPoints(),
                    QPointsAdded=data.size(),
                    Qmin=QPoints<=QPointsAdded?QPoints:QPointsAdded,
                    DoFN=2, DerivN=0,
                    MaxDerivN=this->kynChars.GetMaxDerivN(),
                    ArrayN=this->kynChars.ArrayNByDofNDerivN(DoFN, DerivN);
                //if(this->kynChars.GetIfDoFXExists() && this->data.size()>0 ){
                if(this->kynChars.GetIfDoFXExists() && this->kynChars.GetQPoints()>0 ){
                    for(int PointN=1; PointN<=Qmin; PointN++){
                        this->data[PointN-1][ArrayN-1]=data[PointN-1];
                    }
                    for(int PointN=Qmin+1; PointN<=QPoints; PointN++){
                        this->data[PointN-1][ArrayN-1]=data[PointN-1];
                    }
                }
            }
            void KynematicParams_vb::Set_all_Z( std::vector<double>data ){
                int //QPoints=this->data.size(),
                    QPoints=this->kynChars.GetQPoints(),
                    QPointsAdded=data.size(),
                    Qmin=QPoints<=QPointsAdded?QPoints:QPointsAdded,
                    DoFN=3, DerivN=0,
                    MaxDerivN=this->kynChars.GetMaxDerivN(),
                    ArrayN=this->kynChars.ArrayNByDofNDerivN(DoFN, DerivN);
                //if(this->kynChars.GetIfDoFXExists() && this->data.size()>0 ){
                if(this->kynChars.GetIfDoFXExists() && this->kynChars.GetQPoints()>0 ){
                    for(int PointN=1; PointN<=Qmin; PointN++){
                        this->data[PointN-1][ArrayN-1]=data[PointN-1];
                    }
                    for(int PointN=Qmin+1; PointN<=QPoints; PointN++){
                        this->data[PointN-1][ArrayN-1]=data[PointN-1];
                    }
                }
            }
            void KynematicParams_vb::Set_all_Vx( std::vector<double>data ){
                int //QPoints=this->data.size(),
                    QPoints=this->kynChars.GetQPoints(),
                    QPointsAdded=data.size(),
                    Qmin=QPoints<=QPointsAdded?QPoints:QPointsAdded,
                    DoFN=1, DerivN=1,
                    MaxDerivN=this->kynChars.GetMaxDerivN(),
                    ArrayN=this->kynChars.ArrayNByDofNDerivN(DoFN, DerivN);
                //if(this->kynChars.GetIfDoFXExists() && this->data.size()>0 ){
                if(this->kynChars.GetIfDoFXExists() && this->kynChars.GetQPoints()>0 ){
                    for(int PointN=1; PointN<=Qmin; PointN++){
                        this->data[PointN-1][ArrayN-1]=data[PointN-1];
                    }
                    for(int PointN=Qmin+1; PointN<=QPoints; PointN++){
                        this->data[PointN-1][ArrayN-1]=data[PointN-1];
                    }
                }
            }
            void KynematicParams_vb::Set_all_Vy( std::vector<double>data ){
                int //QPoints=this->data.size(),
                    QPoints=this->kynChars.GetQPoints(),
                    QPointsAdded=data.size(),
                    Qmin=QPoints<=QPointsAdded?QPoints:QPointsAdded,
                    DoFN=2, DerivN=1,
                    MaxDerivN=this->kynChars.GetMaxDerivN(),
                    ArrayN=this->kynChars.ArrayNByDofNDerivN(DoFN, DerivN);
                //if(this->kynChars.GetIfDoFXExists() && this->data.size()>0 ){
                if(this->kynChars.GetIfDoFXExists() && this->kynChars.GetQPoints()>0 ){
                    for(int PointN=1; PointN<=Qmin; PointN++){
                        this->data[PointN-1][ArrayN-1]=data[PointN-1];
                    }
                    for(int PointN=Qmin+1; PointN<=QPoints; PointN++){
                        this->data[PointN-1][ArrayN-1]=data[PointN-1];
                    }
                }
            }
            void KynematicParams_vb::Set_all_Vz( std::vector<double>data ){
                int //QPoints=this->data.size(),
                    QPoints=this->kynChars.GetQPoints(),
                    QPointsAdded=data.size(),
                    Qmin=QPoints<=QPointsAdded?QPoints:QPointsAdded,
                    DoFN=3, DerivN=1,
                    MaxDerivN=this->kynChars.GetMaxDerivN(),
                    ArrayN=this->kynChars.ArrayNByDofNDerivN(DoFN, DerivN);
                //if(this->kynChars.GetIfDoFXExists() && this->data.size()>0 ){
                if(this->kynChars.GetIfDoFXExists() && this->kynChars.GetQPoints()>0 ){
                    for(int PointN=1; PointN<=Qmin; PointN++){
                        this->data[PointN-1][ArrayN-1]=data[PointN-1];
                    }
                    for(int PointN=Qmin+1; PointN<=QPoints; PointN++){
                        this->data[PointN-1][ArrayN-1]=data[PointN-1];
                    }
                }
            }
            void KynematicParams_vb::Set_all_ax( std::vector<double>data ){
                int //QPoints=this->data.size(),
                    QPoints=this->kynChars.GetQPoints(),
                    QPointsAdded=data.size(),
                    Qmin=QPoints<=QPointsAdded?QPoints:QPointsAdded,
                    DoFN=1, DerivN=2,
                    MaxDerivN=this->kynChars.GetMaxDerivN(),
                    ArrayN=this->kynChars.ArrayNByDofNDerivN(DoFN, DerivN);
                //if(this->kynChars.GetIfDoFXExists() && this->data.size()>0 ){
                if(this->kynChars.GetIfDoFXExists() && this->kynChars.GetQPoints()>0 ){
                    for(int PointN=1; PointN<=Qmin; PointN++){
                        this->data[PointN-1][ArrayN-1]=data[PointN-1];
                    }
                    for(int PointN=Qmin+1; PointN<=QPoints; PointN++){
                        this->data[PointN-1][ArrayN-1]=data[PointN-1];
                    }
                }
            }
            void KynematicParams_vb::Set_all_ay( std::vector<double>data ){
                int //QPoints=this->data.size(),
                    QPoints=this->kynChars.GetQPoints(),
                    QPointsAdded=data.size(),
                    Qmin=QPoints<=QPointsAdded?QPoints:QPointsAdded,
                    DoFN=2, DerivN=2,
                    MaxDerivN=this->kynChars.GetMaxDerivN(),
                    ArrayN=this->kynChars.ArrayNByDofNDerivN(DoFN, DerivN);
                //if(this->kynChars.GetIfDoFXExists() && this->data.size()>0 ){
                if(this->kynChars.GetIfDoFXExists() && this->kynChars.GetQPoints()>0 ){
                    for(int PointN=1; PointN<=Qmin; PointN++){
                        this->data[PointN-1][ArrayN-1]=data[PointN-1];
                    }
                    for(int PointN=Qmin+1; PointN<=QPoints; PointN++){
                        this->data[PointN-1][ArrayN-1]=data[PointN-1];
                    }
                }
            }
            void KynematicParams_vb::Set_all_az( std::vector<double>data ){
                int //QPoints=this->data.size(),
                    QPoints=this->kynChars.GetQPoints(),
                    QPointsAdded=data.size(),
                    Qmin=QPoints<=QPointsAdded?QPoints:QPointsAdded,
                    DoFN=3, DerivN=2,
                    MaxDerivN=this->kynChars.GetMaxDerivN(),
                    ArrayN=this->kynChars.ArrayNByDofNDerivN(DoFN, DerivN);
                //if(this->kynChars.GetIfDoFXExists() && this->data.size()>0 ){
                if(this->kynChars.GetIfDoFXExists() && this->kynChars.GetQPoints()>0 ){
                    for(int PointN=1; PointN<=Qmin; PointN++){
                        this->data[PointN-1][ArrayN-1]=data[PointN-1];
                    }
                    for(int PointN=Qmin+1; PointN<=QPoints; PointN++){
                        this->data[PointN-1][ArrayN-1]=data[PointN-1];
                    }
                }
            }
            void KynematicParams_vb::Set_all_FiX( std::vector<double>data ){
                int //QPoints=this->data.size(),
                    QPoints=this->kynChars.GetQPoints(),
                    QPointsAdded=data.size(),
                    Qmin=QPoints<=QPointsAdded?QPoints:QPointsAdded,
                    DoFN=4, DerivN=0,
                    MaxDerivN=this->kynChars.GetMaxDerivN(),
                    ArrayN=this->kynChars.ArrayNByDofNDerivN(DoFN, DerivN);
                //if(this->kynChars.GetIfDoFXExists() && this->data.size()>0 ){
                if(this->kynChars.GetIfDoFXExists() && this->kynChars.GetQPoints()>0 ){
                    for(int PointN=1; PointN<=Qmin; PointN++){
                        this->data[PointN-1][ArrayN-1]=data[PointN-1];
                    }
                    for(int PointN=Qmin+1; PointN<=QPoints; PointN++){
                        this->data[PointN-1][ArrayN-1]=data[PointN-1];
                    }
                }
            }
            void KynematicParams_vb::Set_all_FiY( std::vector<double>data ){
                int //QPoints=this->data.size(),
                    QPoints=this->kynChars.GetQPoints(),
                    QPointsAdded=data.size(),
                    Qmin=QPoints<=QPointsAdded?QPoints:QPointsAdded,
                    DoFN=5, DerivN=0,
                    MaxDerivN=this->kynChars.GetMaxDerivN(),
                    ArrayN=this->kynChars.ArrayNByDofNDerivN(DoFN, DerivN);
                //if(this->kynChars.GetIfDoFXExists() && this->data.size()>0 ){
                if(this->kynChars.GetIfDoFXExists() && this->kynChars.GetQPoints()>0 ){
                    for(int PointN=1; PointN<=Qmin; PointN++){
                        this->data[PointN-1][ArrayN-1]=data[PointN-1];
                    }
                    for(int PointN=Qmin+1; PointN<=QPoints; PointN++){
                        this->data[PointN-1][ArrayN-1]=data[PointN-1];
                    }
                }
            }
            void KynematicParams_vb::Set_all_FiZ( std::vector<double>data ){
                int //QPoints=this->data.size(),
                    QPoints=this->kynChars.GetQPoints(),
                    QPointsAdded=data.size(),
                    Qmin=QPoints<=QPointsAdded?QPoints:QPointsAdded,
                    DoFN=6, DerivN=0,
                    MaxDerivN=this->kynChars.GetMaxDerivN(),
                    ArrayN=this->kynChars.ArrayNByDofNDerivN(DoFN, DerivN);
                //if(this->kynChars.GetIfDoFXExists() && this->data.size()>0 ){
                if(this->kynChars.GetIfDoFXExists() && this->kynChars.GetQPoints()>0 ){
                    for(int PointN=1; PointN<=Qmin; PointN++){
                        this->data[PointN-1][ArrayN-1]=data[PointN-1];
                    }
                    for(int PointN=Qmin+1; PointN<=QPoints; PointN++){
                        this->data[PointN-1][ArrayN-1]=data[PointN-1];
                    }
                }
            }
            void KynematicParams_vb::Set_all_wx( std::vector<double>data ){
                int //QPoints=this->data.size(),
                    QPoints=this->kynChars.GetQPoints(),
                    QPointsAdded=data.size(),
                    Qmin=QPoints<=QPointsAdded?QPoints:QPointsAdded,
                    DoFN=4, DerivN=1,
                    MaxDerivN=this->kynChars.GetMaxDerivN(),
                    ArrayN=this->kynChars.ArrayNByDofNDerivN(DoFN, DerivN);
                //if(this->kynChars.GetIfDoFXExists() && this->data.size()>0 ){
                if(this->kynChars.GetIfDoFXExists() && this->kynChars.GetQPoints()>0 ){
                    for(int PointN=1; PointN<=Qmin; PointN++){
                        this->data[PointN-1][ArrayN-1]=data[PointN-1];
                    }
                    for(int PointN=Qmin+1; PointN<=QPoints; PointN++){
                        this->data[PointN-1][ArrayN-1]=data[PointN-1];
                    }
                }
            }
            void KynematicParams_vb::Set_all_wy( std::vector<double>data ){
                int //QPoints=this->data.size(),
                    QPoints=this->kynChars.GetQPoints(),
                    QPointsAdded=data.size(),
                    Qmin=QPoints<=QPointsAdded?QPoints:QPointsAdded,
                    DoFN=5, DerivN=1,
                    MaxDerivN=this->kynChars.GetMaxDerivN(),
                    ArrayN=this->kynChars.ArrayNByDofNDerivN(DoFN, DerivN);
                //if(this->kynChars.GetIfDoFXExists() && this->data.size()>0 ){
                if(this->kynChars.GetIfDoFXExists() && this->kynChars.GetQPoints()>0 ){
                    for(int PointN=1; PointN<=Qmin; PointN++){
                        this->data[PointN-1][ArrayN-1]=data[PointN-1];
                    }
                    for(int PointN=Qmin+1; PointN<=QPoints; PointN++){
                        this->data[PointN-1][ArrayN-1]=data[PointN-1];
                    }
                }
            }
            void KynematicParams_vb::Set_all_wz( std::vector<double>data ){
                int //QPoints=this->data.size(),
                    QPoints=this->kynChars.GetQPoints(),
                    QPointsAdded=data.size(),
                    Qmin=QPoints<=QPointsAdded?QPoints:QPointsAdded,
                    DoFN=6, DerivN=1,
                    MaxDerivN=this->kynChars.GetMaxDerivN(),
                    ArrayN=this->kynChars.ArrayNByDofNDerivN(DoFN, DerivN);
                //if(this->kynChars.GetIfDoFXExists() && this->data.size()>0 ){
                if(this->kynChars.GetIfDoFXExists() && this->kynChars.GetQPoints()>0 ){
                    for(int PointN=1; PointN<=Qmin; PointN++){
                        this->data[PointN-1][ArrayN-1]=data[PointN-1];
                    }
                    for(int PointN=Qmin+1; PointN<=QPoints; PointN++){
                        this->data[PointN-1][ArrayN-1]=data[PointN-1];
                    }
                }
            }
            void KynematicParams_vb::Set_all_ex( std::vector<double>data ){
                int //QPoints=this->data.size(),
                    QPoints=this->kynChars.GetQPoints(),
                    QPointsAdded=data.size(),
                    Qmin=QPoints<=QPointsAdded?QPoints:QPointsAdded,
                    DoFN=4, DerivN=2,
                    MaxDerivN=this->kynChars.GetMaxDerivN(),
                    ArrayN=this->kynChars.ArrayNByDofNDerivN(DoFN, DerivN);
                //if(this->kynChars.GetIfDoFXExists() && this->data.size()>0 ){
                if(this->kynChars.GetIfDoFXExists() && this->kynChars.GetQPoints()>0 ){
                    for(int PointN=1; PointN<=Qmin; PointN++){
                        this->data[PointN-1][ArrayN-1]=data[PointN-1];
                    }
                    for(int PointN=Qmin+1; PointN<=QPoints; PointN++){
                        this->data[PointN-1][ArrayN-1]=data[PointN-1];
                    }
                }
            }
            void KynematicParams_vb::Set_all_ey( std::vector<double>data ){
                int //QPoints=this->data.size(),
                    QPoints=this->kynChars.GetQPoints(),
                    QPointsAdded=data.size(),
                    Qmin=QPoints<=QPointsAdded?QPoints:QPointsAdded,
                    DoFN=5, DerivN=2,
                    MaxDerivN=this->kynChars.GetMaxDerivN(),
                    ArrayN=this->kynChars.ArrayNByDofNDerivN(DoFN, DerivN);
                //if(this->kynChars.GetIfDoFXExists() && this->data.size()>0 ){
                if(this->kynChars.GetIfDoFXExists() && this->kynChars.GetQPoints()>0 ){
                    for(int PointN=1; PointN<=Qmin; PointN++){
                        this->data[PointN-1][ArrayN-1]=data[PointN-1];
                    }
                    for(int PointN=Qmin+1; PointN<=QPoints; PointN++){
                        this->data[PointN-1][ArrayN-1]=data[PointN-1];
                    }
                }
            }
            void KynematicParams_vb::Set_all_ez( std::vector<double>data ){
                int //QPoints=this->data.size(),
                    QPoints=this->kynChars.GetQPoints(),
                    QPointsAdded=data.size(),
                    Qmin=QPoints<=QPointsAdded?QPoints:QPointsAdded,
                    DoFN=6, DerivN=2,
                    MaxDerivN=this->kynChars.GetMaxDerivN(),
                    ArrayN=this->kynChars.ArrayNByDofNDerivN(DoFN, DerivN);
                //if(this->kynChars.GetIfDoFXExists() && this->data.size()>0 ){
                if(this->kynChars.GetIfDoFXExists() && this->kynChars.GetQPoints()>0 ){
                    for(int PointN=1; PointN<=Qmin; PointN++){
                        this->data[PointN-1][ArrayN-1]=data[PointN-1];
                    }
                    for(int PointN=Qmin+1; PointN<=QPoints; PointN++){
                        this->data[PointN-1][ArrayN-1]=data[PointN-1];
                    }
                }
            }
            //
            void KynematicParams_vb::SetValsByDoFNDerivNPointN(double val, int DoFN, int DerivN, int PointN){
                int N;
                if(PointN<0){
                    PointN+=(this->GetQPoints()+1);
                }
                if(PointN>=1 && PointN<=this->GetQPoints() && DerivN>=0 && DerivN<=this->GetMaxDerivN() && this->kynChars.GetIfDoFNExists(DoFN)){
                    N=this->kynChars.ArrayNByDofNDerivN(DoFN, DerivN);
                    this->data[PointN-1][N-1]=val;
                }
            }

            double KynematicParams_vb::GetValsByDoFNDerivNPointN(int DoFN, int DerivN, int PointN){
                double val=0;
                int N;
                if(PointN<0){
                    PointN+=(this->GetQPoints()+1);
                }
                if(PointN>=1 && PointN<=this->GetQPoints() && DerivN>=0 && DerivN<=this->GetMaxDerivN() && this->kynChars.GetIfDoFNExists(DoFN)){
                    N=this->kynChars.ArrayNByDofNDerivN(DoFN, DerivN);
                    val=this->data[PointN-1][N-1];
                }
                return val;
            }

            //
            //void KynematicParams_vb::SetPointN(int N, double*data){}
            void KynematicParams_vb::SetPointN(int N, std::vector<double>data){
                if(N>=1 && N<=this->GetQPoints()){
                    this->data[N-1]=data;
                }
            }
            std::vector<double>KynematicParams_vb::GetPointN(int N)const{
                std::vector<double>point;
                int L=this->GetArrayLengthOfPointKynematicParams();
                for(int i=1; i<=L; i++){
                    point.push_back(0.0);
                }
                if(N>=1 && N<=this->GetQPoints()){
                    for(int i=1; i<=L; i++){
                        point[i-1]=this->data[N-1][i-1];
                    }
                }
                return point;
            }
            //
            int KynematicParams_vb::GetQPoints()const{
                //return this->data.size();
                return this->kynChars.GetQPoints();
            }


            int KynematicParams_vb::GetArrayLengthOfPointKynematicParams()const{
                return this->kynChars.GetArrayLengthOfPointKynematicParams();
            }

            //void KynematicParams2::AddPoint(std::vector<double>data){
                //this->data.push_back(data);
                //int Q=this->GetQPoints(), L=this->kynChars.GetArrayLengthOfPointKynematicParams();
                //std::vector<double>point;//=new double[L];
                ////double**y;
                //double val;
                //y=new double*[Q+1];
                //for(int i=1; i<=Q+1; i++){
                //    y[i-1]=new double[L];
                //}
                //
                //if(data!=NULL){
                //if(data.size()!=0){
                //    for(int i=1; i<=L; i++){
                //        val=data[i-1];
                //        //point[i-1]=val;
                //        point.push_back(val);
                //        //point[i-1]=data[i-1];
                //    }
                //}else{
                //    for(int i=1; i<=L; i++){
                //        point[i-1]=0;
                //    }
                //}
                //
                //if(this->data!=NULL){
                //if(this->data.size()!=0){
                //    for(int i=1; i<=Q; i++){
                //        for(int j=1; j<=L; j++){
                //            val=this->data[i-1][j-1];
                //            y[i-1][j-1]=val;
                //            //y[i-1][j-1]=this->data[i-1][j-1];
                //        }
                //    }
                //}
                //
                //for(int i=1; i<=L; i++){
                //    y[Q+1-1][i-1]=point[i-1];
                //}
                //
                //if(this->data!=NULL){
               //     for(int i=1; i<=Q+1; i++){
               //         delete[]this->data[i-1];
               //     }
                //    delete[]this->data;
                //}
                //this->data=y;
                ////
                //delete[]point;
                //point=NULL;
                //
                //this->kynChars.QPoints++;
                //
                //int LastN=this->GetQPoints();
                //LastN=this->data.size();
                //LastN=this->kynChars.GetQPoints();
                //QString sn;
                //std::cout<<"Last val in last point = "<<this->data[LastN-1][Q-1]<<std::endl;
                //sn=FloatToStr(this->data[LastN-1][Q-1]);
                //std::cout<<"Last val in last point = "<<sn.toStdString()<<std::endl;
               // std::cout<<"Last val in last point = "<<sn.toStdString().c_str()<<std::endl;
                //for(int i=1; i<=Q; i++){
                //    std::cout<<this->data[LastN-1][i-1]<<"; ";
                //}
                //std::cout<<std::endl;
                //std::cout<<"This poinbt was added. Add finishes working"<<std::endl;
            //}
            //
            std::vector<double>KynematicParams_vb::GetRowVectorOfPointDataByDoFsAsPtrs(double*X, double*Y, double*Z, double*FiX, double*FiY, double*FiZ){
                std::vector<double> row;
                int MaxDerivN=this->GetMaxDerivN();
                //0
                if(this->GetIfDoFXExists()){
                    if(X!=NULL){
                        row.push_back(X[0]);
                    }else{
                        row.push_back(0);
                    }
                }
                if(this->GetIfDoFYExists()){
                    if(Y!=NULL){
                        row.push_back(Y[0]);
                    }else{
                        row.push_back(0);
                    }
                }
                if(this->GetIfDoFZExists()){
                    if(Z!=NULL){
                        row.push_back(Z[0]);
                    }else{
                        row.push_back(0);
                    }
                }
                if(this->GetIfDoFFiXExists()){
                    if(FiX!=NULL){
                        row.push_back(FiX[0]);
                    }else{
                        row.push_back(0);
                    }
                }
                if(this->GetIfDoFFiYExists()){
                    if(FiY!=NULL){
                        row.push_back(FiY[0]);
                    }else{
                        row.push_back(0);
                    }
                }
                if(this->GetIfDoFFiZExists()){
                    if(FiZ!=NULL){
                        row.push_back(FiZ[0]);
                    }else{
                        row.push_back(0);
                    }
                }
                if(MaxDerivN>=1){
                    //1
                    if(this->GetIfDoFXExists()){
                        if(X!=NULL){
                            row.push_back(X[1]);
                        }else{
                            row.push_back(0);
                        }
                    }
                    if(this->GetIfDoFYExists()){
                        if(Y!=NULL){
                            row.push_back(Y[1]);
                        }else{
                            row.push_back(0);
                        }
                    }
                    if(this->GetIfDoFZExists()){
                        if(Z!=NULL){
                            row.push_back(Z[1]);
                        }else{
                            row.push_back(0);
                        }
                    }
                    if(this->GetIfDoFFiXExists()){
                        if(FiX!=NULL){
                            row.push_back(FiX[1]);
                        }else{
                            row.push_back(0);
                        }
                    }
                    if(this->GetIfDoFFiYExists()){
                        if(FiY!=NULL){
                            row.push_back(FiY[1]);
                        }else{
                            row.push_back(0);
                        }
                    }
                    if(this->GetIfDoFFiZExists()){
                        if(FiZ!=NULL){
                            row.push_back(FiZ[1]);
                        }else{
                            row.push_back(0);
                        }
                    }
                    if(MaxDerivN>=2){
                        //2
                        if(this->GetIfDoFXExists()){
                            if(X!=NULL){
                                row.push_back(X[2]);
                            }else{
                                row.push_back(0);
                            }
                        }
                        if(this->GetIfDoFYExists()){
                            if(Y!=NULL){
                                row.push_back(Y[2]);
                            }else{
                                row.push_back(0);
                            }
                        }
                        if(this->GetIfDoFZExists()){
                            if(Z!=NULL){
                                row.push_back(Z[2]);
                            }else{
                                row.push_back(0);
                            }
                        }
                        if(this->GetIfDoFFiXExists()){
                            if(FiX!=NULL){
                                row.push_back(FiX[2]);
                            }else{
                                row.push_back(0);
                            }
                        }
                        if(this->GetIfDoFFiYExists()){
                            if(FiY!=NULL){
                                row.push_back(FiY[2]);
                            }else{
                                row.push_back(0);
                            }
                        }
                        if(this->GetIfDoFFiZExists()){
                            if(FiZ!=NULL){
                                row.push_back(FiZ[2]);
                            }else{
                                row.push_back(0);
                            }
                        }
                    }
                }
                return row;
            }
            std::vector<double>KynematicParams_vb::GetRowVectorOfPointDataByDoFsAsVecs(std::vector<double>X, std::vector<double>Y, std::vector<double>Z, std::vector<double>FiX, std::vector<double>FiY, std::vector<double>FiZ){
                std::vector<double> row;
                int MaxDerivN=this->GetMaxDerivN();
                //0
                if(this->GetIfDoFXExists()){
                    if(X.size()!=0){
                        row.push_back(X[0]);
                    }else{
                        row.push_back(0);
                    }
                }
                if(this->GetIfDoFYExists()){
                    if(Y.size()!=0){
                        row.push_back(Y[0]);
                    }else{
                        row.push_back(0);
                    }
                }
                if(this->GetIfDoFZExists()){
                    if(Z.size()!=0){
                        row.push_back(Z[0]);
                    }else{
                        row.push_back(0);
                    }
                }
                if(this->GetIfDoFFiXExists()){
                    if(FiX.size()!=0){
                        row.push_back(FiX[0]);
                    }else{
                        row.push_back(0);
                    }
                }
                if(this->GetIfDoFFiYExists()){
                    if(FiY.size()!=0){
                        row.push_back(FiY[0]);
                    }else{
                        row.push_back(0);
                    }
                }
                if(this->GetIfDoFFiZExists()){
                    if(FiZ.size()!=0){
                        row.push_back(FiZ[0]);
                    }else{
                        row.push_back(0);
                    }
                }
                if(MaxDerivN>=1){
                    //1
                    if(this->GetIfDoFXExists()){
                        if(X.size()>=2){
                            row.push_back(X[1]);
                        }else{
                            row.push_back(0);
                        }
                    }
                    if(this->GetIfDoFYExists()){
                        if(Y.size()>=2){
                            row.push_back(Y[1]);
                        }else{
                            row.push_back(0);
                        }
                    }
                    if(this->GetIfDoFZExists()){
                        if(Z.size()>=2){
                            row.push_back(Z[1]);
                        }else{
                            row.push_back(0);
                        }
                    }
                    if(this->GetIfDoFFiXExists()){
                        if(FiX.size()>=2){
                            row.push_back(FiX[1]);
                        }else{
                            row.push_back(0);
                        }
                    }
                    if(this->GetIfDoFFiYExists()){
                        if(FiY.size()>=1){
                            row.push_back(FiY[1]);
                        }else{
                            row.push_back(0);
                        }
                    }
                    if(this->GetIfDoFFiZExists()){
                        if(FiZ.size()>=2){
                            row.push_back(FiZ[1]);
                        }else{
                            row.push_back(0);
                        }
                    }
                    if(MaxDerivN>=2){
                        //2
                        if(this->GetIfDoFXExists()){
                            if(X.size()>=3){
                                row.push_back(X[2]);
                            }else{
                                row.push_back(0);
                            }
                        }
                        if(this->GetIfDoFYExists()){
                            if(Y.size()>=3){
                                row.push_back(Y[2]);
                            }else{
                                row.push_back(0);
                            }
                        }
                        if(this->GetIfDoFZExists()){
                            if(Z.size()>=3){
                                row.push_back(Z[2]);
                            }else{
                                row.push_back(0);
                            }
                        }
                        if(this->GetIfDoFFiXExists()){
                            if(FiX.size()>=3){
                                row.push_back(FiX[2]);
                            }else{
                                row.push_back(0);
                            }
                        }
                        if(this->GetIfDoFFiYExists()){
                            if(FiY.size()>=3){
                                row.push_back(FiY[2]);
                            }else{
                                row.push_back(0);
                            }
                        }
                        if(this->GetIfDoFFiZExists()){
                            if(FiZ.size()>=3){
                                row.push_back(FiZ[2]);
                            }else{
                                row.push_back(0);
                            }
                        }
                    }
                }
                return row;
            }
            //
            void KynematicParams_vb::AddPoint(double*data, int LGiven){
                double *DfltVal=new double, valCur;
                int L=this->GetArrayLengthOfPointKynematicParams(), QPoints=this->GetQPoints();
                std::vector<double>PointToAddOrIns;
                DfltVal=0;
                //KynematicCharacteristics1 kynAdded
                if(LGiven==0 && data!=NULL)LGiven=L;
                //template <typename T>std::vector<T> Array1DGetSubArray(T*x, int Lold, int Lnew=0, int FromN=1, int FirstDefaultValues=0, T*DfltValParam=NULL)
                PointToAddOrIns=Array1DGetSubArray_byLs(data, LGiven, this->GetArrayLengthOfPointKynematicParams(), 1, 0, DfltVal);
                Array1DAddElement(this->data, PointToAddOrIns);
                this->kynChars.AddOrInsPoint();
                delete DfltVal;
            }
            void KynematicParams_vb::AddPoint(std::vector<double>data){
                double*DfltVal=new double;
                DfltVal=0;
                int L=this->GetArrayLengthOfPointKynematicParams();
                std::vector<double>PointToAddOrIns;
                PointToAddOrIns=Array1DGetSubArray_byLs(data, L, 1, 0, DfltVal);
                this->data.push_back(PointToAddOrIns);
                this->kynChars.AddOrInsPoint();
                delete DfltVal;
            }
            void  KynematicParams_vb::AddPoint(double*data, const KynematicCharacteristics1 kynCharsParam){
                int DoFsOld[]={0, 0, 0, 0, 0, 0}, DoFsNew[]={0, 0, 0, 0, 0, 0}, DoFsCharNumOld, DoFsCharNumNew, QDoFsOld, QDoFsNew, MaxDerivNOld, MaxDerivNNew, QPointsOld, QPointsNew;
                std::vector<double>PointToAddOrIns;
                int L=this->GetArrayLengthOfPointKynematicParams(), Lold=kynCharsParam.GetArrayLengthOfPointKynematicParams();
                //void Get(int DoFs[], int&MaxDerivN, int&QPoints, int&QDoFs, int &DoFsCharNum) const ;
                kynCharsParam.Get(DoFsOld, MaxDerivNOld, QPointsOld, QDoFsOld, DoFsCharNumOld);
                this->kynChars.Get(DoFsNew, MaxDerivNNew, QPointsNew, QDoFsNew, DoFsCharNumNew);
                if(DoFsCharNumOld==DoFsCharNumNew && MaxDerivNOld==MaxDerivNNew){
                    this->AddPoint(data);
                }else{
                    //Array1DGetSubArray_byLs
                    PointToAddOrIns=Array1DGetSubArray_byLs(data, Lold, L);
                    this->AddPoint(PointToAddOrIns);
                }
            }

            void KynematicParams_vb::AddPoint(double*X, double*Y, double*Z, double*FiX, double*FiY, double*FiZ){
                std::vector<double>PointToAddOrIns=GetRowVectorOfPointDataByDoFsAsPtrs(X, Y, Z, FiX, FiY, FiZ);
                this->data.push_back(PointToAddOrIns);
                this->kynChars.AddOrInsPoint();
            }
            void KynematicParams_vb::AddPoint(std::vector<double>X, std::vector<double>Y, std::vector<double>Z, std::vector<double>FiX, std::vector<double>FiY, std::vector<double>FiZ){
                std::vector<double>PointToAddOrIns=GetRowVectorOfPointDataByDoFsAsVecs(X, Y, Z, FiX, FiY, FiZ);
                this->data.push_back(PointToAddOrIns);
                this->kynChars.AddOrInsPoint();
            }
            //
            void KynematicParams_vb::InsPoint(int N, double*given, int LGiven){
                double *DfltVal, valCur;
                int L=this->GetArrayLengthOfPointKynematicParams(), QPoints=this->GetQPoints();
                std::vector<double>PointToAddOrIns;
                //KynematicCharacteristics1 kynAdded
                QPoints=this->GetQPoints();
                if(N<0)N+=QPoints;
                if(N>=1 && N<=QPoints){
                    DfltVal=new double;
                    DfltVal=0;
                    if(LGiven==0 && given!=NULL)LGiven=L;
                    PointToAddOrIns=Array1DGetSubArray_byLs(given, LGiven, L, 1, 0, DfltVal);
                    //
                    Array1DInsElementToN(this->data, N,PointToAddOrIns);
                    //
                    this->kynChars.AddOrInsPoint();
                    //
                    delete DfltVal;
                }
            }

            void KynematicParams_vb::InsPoint(int N, std::vector<double>data){
                int QPoints=this->GetQPoints();
                double *DfltVal;
                if(N<0)N+=QPoints;
                if(N>=1 && N<=QPoints){
                    DfltVal=new double;
                    DfltVal=0;
                    std::vector<double>PointToAddOrIns=Array1DGetSubArray_byLs(data, this->GetArrayLengthOfPointKynematicParams(), 1, 0, DfltVal);
                    Array1DInsElementToN(this->data, N, PointToAddOrIns);
                    this->kynChars.AddOrInsPoint();
                    delete DfltVal;
                }
            }
            //





            void KynematicParams_vb::InsPoint(int N, double*X, double*Y, double*Z, double*FiX, double*FiY, double*FiZ){
                std::vector<double>PointToAddOrIns=GetRowVectorOfPointDataByDoFsAsPtrs(X, Y, Z, FiX, FiY, FiZ);
                if(N>=1 && N<=this->GetQPoints()){
                    this->data.push_back(PointToAddOrIns);
                    this->kynChars.AddOrInsPoint();
                }
            }
            void KynematicParams_vb::InsPoint(int N, std::vector<double>X, std::vector<double>Y, std::vector<double>Z, std::vector<double>FiX, std::vector<double>FiY, std::vector<double>FiZ){
                std::vector<double>PointToAddOrIns=GetRowVectorOfPointDataByDoFsAsVecs(X, Y, Z, FiX, FiY, FiZ);
                if(N>=1 && N<=this->GetQPoints()){
                    Array1DInsElementToN(this->data, N, PointToAddOrIns);
                    this->kynChars.AddOrInsPoint();
                }
            }
            //
            void KynematicParams_vb::DelPoint(int N){
                int QPoints=this->GetQPoints();
                if(N<0)N+=QPoints;
                Array1DDelElementFromN(this->data, N);
                this->kynChars.DelPoint();
            }

            //

            void KynematicParams_vb::SetPoint(int N, std::vector<double>data){
                if(N>=1 && N<=this->data.size()){
                    this->data[N-1]=data;
                }
            }
            //

            QString KynematicParams_vb::PointNToString(int N, bool showAll){
                QString s="", sn;
                double val;
                int MaxDerivN=this->kynChars.GetMaxDerivN();
                if(N>=1 && N<=this->GetQPoints()){
                   //0
                   if(this->kynChars.GetIfDoFXExists() || showAll){
                       s=s+" x=";
                       val=this->get_x(N);
                       sn=FloatToStr(val);
                       //s=s+FloatToStr(this->get_x(N));
                       s=s+sn;
                   }
                   if(this->kynChars.GetIfDoFYExists() || showAll){
                       s=s+" y=";
                       val=this->get_y(N);
                       sn=FloatToStr(val);
                       //s=s+FloatToStr(this->get_y(N));
                       s=s+sn;
                   }
                   if(this->kynChars.GetIfDoFZExists() || showAll){
                       s=s+" z=";
                       val=this->get_z(N);
                       sn=FloatToStr(val);
                       //s=s+FloatToStr(this->get_z(N));
                       s=s+sn;
                   }
                   if(this->kynChars.GetIfDoFFiXExists() || showAll){
                       s=s+" fiX=";
                       val=this->get_fiX(N);
                       sn=FloatToStr(val);
                       //s=s+FloatToStr(this->get_fiX(N));
                       s=s+sn;
                   }
                   if(this->kynChars.GetIfDoFFiYExists() || showAll){
                       s=s+" fiY=";
                       val=this->get_fiY(N);
                       sn=FloatToStr(val);
                       //s=s+FloatToStr(this->get_fiY(N));
                       s=s+sn;
                   }
                   if(this->kynChars.GetIfDoFFiZExists() || showAll){
                       s=s+" fiZ=";
                       val=this->get_fiZ(N);
                       sn=FloatToStr(val);
                       //s=s+FloatToStr(this->get_fiZ(N));
                       s=s+sn;
                   }
                   //1
                   if((this->kynChars.GetIfDoFXExists() && MaxDerivN>=1) || showAll){
                       s=s+" vx=";
                       val=this->get_vx(N);
                       sn=FloatToStr(val);
                       //s=s+FloatToStr(this->get_vx(N));
                       s=s+sn;
                   }
                   if((this->kynChars.GetIfDoFYExists() && MaxDerivN>=1) || showAll){
                       s=s+" vy=";
                       val=this->get_vy(N);
                       sn=FloatToStr(val);
                       //s=s+FloatToStr(this->get_vy(N));
                       s=s+sn;
                   }
                   if((this->kynChars.GetIfDoFZExists() && MaxDerivN>=1) || showAll){
                       s=s+" vz=";
                       val=this->get_vz(N);
                       sn=FloatToStr(val);
                       //s=s+FloatToStr(this->get_vz(N));
                       s=s+sn;
                   }
                   if((this->kynChars.GetIfDoFFiXExists() && MaxDerivN>=1) || showAll){
                       s=s+" wx=";
                       val=this->get_wx(N);
                       sn=FloatToStr(val);
                       //s=s+FloatToStr(this->get_wx(N));
                       s=s+sn;
                   }
                   if((this->kynChars.GetIfDoFFiYExists() && MaxDerivN>=1) || showAll){
                       s=s+" wy=";
                       s=s+FloatToStr(this->get_wy(N));
                   }
                   if((this->kynChars.GetIfDoFFiZExists() && MaxDerivN>1) || showAll){
                       s=s+" wz=";
                       s=s+FloatToStr(this->get_wz(N));
                   }
                   //2
                   if((this->kynChars.GetIfDoFXExists() && MaxDerivN>=2) || showAll){
                       s=s+" ax=";
                       s=s+FloatToStr(this->get_ax(N));
                   }
                   if((this->kynChars.GetIfDoFYExists() && MaxDerivN>=2) || showAll){
                       s=s+" ay=";
                       s=s+FloatToStr(this->get_ay(N));
                   }
                   if((this->kynChars.GetIfDoFZExists() && MaxDerivN>=2) || showAll){
                       s=s+" az=";
                       s=s+FloatToStr(this->get_az(N));
                   }
                   if((this->kynChars.GetIfDoFFiXExists() && MaxDerivN>=2) || showAll){
                       s=s+" ex=";
                       s=s+FloatToStr(this->get_ex(N));
                   }
                   if((this->kynChars.GetIfDoFFiYExists() && MaxDerivN>=2) || showAll){
                       s=s+" ey=";
                       s=s+FloatToStr(this->get_ey(N));
                   }
                   if((this->kynChars.GetIfDoFFiZExists() && MaxDerivN>2) || showAll){
                       s=s+" ez=";
                       s=s+FloatToStr(this->get_ez(N));
                   }
                }
                return s;
            }

            QString KynematicParams_vb::DoFNDerivNToString(int DoFN, int DerivN, QString delim) {
                QString s="", sn, sname;
                double val;
                int MaxDerivN=this->kynChars.GetMaxDerivN(), QPoints=this->GetQPoints();
                int N=this->kynChars.ArrayNByDofNDerivN(DoFN, DerivN);
                if(this->kynChars.GetIfDoFNExists(DoFN) && this->kynChars.GetMaxDerivN()>=DerivN){
                    switch(DoFN){
                        case 1:
                            switch(DerivN){
                                case 0:
                                    sname="x";
                                break;
                                case 1:
                                    sname="vx";
                                break;
                                case 2:
                                    sname="ax";
                                break;
                            }
                        break;
                        case 2:
                        switch(DerivN){
                            case 0:
                                sname="y";
                            break;
                            case 1:
                                sname="vy";
                            break;
                            case 2:
                                sname="ay";
                            break;
                        }
                        break;
                        case 3:
                        switch(DerivN){
                            case 0:
                                sname="z";
                            break;
                            case 1:
                                sname="vz";
                            break;
                            case 2:
                                sname="az";
                            break;
                        }
                        break;
                        case 4:
                        switch(DerivN){
                            case 0:
                                sname="fix";
                            break;
                            case 1:
                                sname="wx";
                            break;
                            case 2:
                                sname="ex";
                            break;
                        }
                        break;
                        case 5:
                        switch(DerivN){
                            case 0:
                                sname="fiy";
                            break;
                            case 1:
                                sname="wy";
                            break;
                            case 2:
                                sname="ey";
                            break;
                        }
                        break;
                        case 6:
                        switch(DerivN){
                            case 0:
                                sname="fiz";
                            break;
                            case 1:
                                sname="wz";
                            break;
                            case 2:
                                sname="ez";
                            break;
                        }//switch ine
                        break;
                   }//switch ext
                   for(int i=1; i<=QPoints-1; i++){
                       s+=sname;
                       s+=IntToStr(i);
                       s+="=";
                       val=this->GetValsByDoFNDerivNPointN(DoFN, DerivN, i);
                       sn=FloatToStr(val);
                       s+=sn;
                       s+=delim;
                   }
                   s+=sname;
                   s+=IntToStr(QPoints);
                   s+="=";
                   val=this->GetValsByDoFNDerivNPointN(DoFN, DerivN, QPoints);
                   sn=FloatToStr(val);
                   s+=sn;
                }//switch ext
                return s;
            }

            QString KynematicParams_vb::get_all_x_as_String(QString delim){ return DoFNDerivNToString(1, 0, delim); }
            QString KynematicParams_vb::get_all_vx_as_String(QString delim){ return DoFNDerivNToString(1, 1, delim); }
            QString KynematicParams_vb::get_all_ax_as_String(QString delim){ return DoFNDerivNToString(1, 2, delim); }
            QString KynematicParams_vb::get_all_fix_as_String(QString delim){ return DoFNDerivNToString(4, 0, delim); }
            QString KynematicParams_vb::get_all_wx_as_String(QString delim){ return DoFNDerivNToString(4, 1, delim); }
            QString KynematicParams_vb::get_all_ex_as_String(QString delim){ return DoFNDerivNToString(4, 2, delim); }

            QString KynematicParams_vb::get_all_y_as_String(QString delim){ return DoFNDerivNToString(2, 0, delim); }
            QString KynematicParams_vb::get_all_vy_as_String(QString delim){ return DoFNDerivNToString(2, 1, delim); }
            QString KynematicParams_vb::get_all_ay_as_String(QString delim){ return DoFNDerivNToString(2, 2, delim); }
            QString KynematicParams_vb::get_all_fiy_as_String(QString delim){ return DoFNDerivNToString(5, 0, delim); }
            QString KynematicParams_vb::get_all_wy_as_String(QString delim){ return DoFNDerivNToString(5, 1, delim); }
            QString KynematicParams_vb::get_all_ey_as_String(QString delim){ return DoFNDerivNToString(5, 2, delim); }

            QString KynematicParams_vb::get_all_z_as_String(QString delim){ return DoFNDerivNToString(3, 0, delim); }
            QString KynematicParams_vb::get_all_vz_as_String(QString delim){ return DoFNDerivNToString(3, 1, delim); }
            QString KynematicParams_vb::get_all_az_as_String(QString delim){ return DoFNDerivNToString(3, 2, delim); }
            QString KynematicParams_vb::get_all_fiz_as_String(QString delim){ return DoFNDerivNToString(6, 0, delim); }
            QString KynematicParams_vb::get_all_wz_as_String(QString delim){ return DoFNDerivNToString(6, 1, delim); }
            QString KynematicParams_vb::get_all_ez_as_String(QString delim){ return DoFNDerivNToString(6, 2, delim); }

            QString KynematicParams_vb::ParamsToString(){ return this->kynChars.ToString(); }


                        void KynematicParams_vb::ShowConsole(){
                            int QPoints=this->GetQPoints();
                            QString qs;
                            qs=this->ParamsToString();
                            std::cout<<"Parameters: "<<qs.toStdString()<<std::endl;
                            qs="";
                            for(int i=1; i<=QPoints; i++){
                                qs=this->PointNToString(i);
                                std::cout<<"Point "<<i<<": "<<qs.toStdString()<<std::endl;
                            }
                        }

                    //};
                        
    std::vector<double>KynematicParams_vb::GetRowForNewChars(std::vector<double>given, int DoFsCharNumOld, int MaxDerivNOld, int DoFsCharNumNew, int MaxDerivNNew){
        std::vector<double>R;
        int DoFsOld[]={0, 0, 0, 0, 0, 0},  DoFsNew[]={0, 0, 0, 0, 0, 0}, order=0, QPointsOld=1, QPointsNew=1, arrNOld, DoFN, DerivN, QDoFsOld=0, QDoFsNew=0;
        double val;
        KynematicCharacteristics1 kynCharsOld, kynCharsNew;
        kynCharsOld.Set(DoFsCharNumOld, MaxDerivNOld, QPointsOld);
        kynCharsNew.Set(DoFsCharNumNew, MaxDerivNNew, QPointsNew);
        kynCharsOld.Get(DoFsOld, MaxDerivNOld, QPointsOld, QDoFsOld, DoFsCharNumOld);
        kynCharsNew.Get(DoFsNew, MaxDerivNNew, QPointsNew, QDoFsNew, DoFsCharNumNew);
        //
        for(DerivN=0; DerivN<=MaxDerivNNew; DerivN++){
            for(DoFN=1; DoFN<=6; DoFN++){
                if(DoFsNew[DoFN-1]==1){
                    arrNOld=kynCharsOld.ArrayNByDofNDerivN(DoFN, DerivN);
                    if(arrNOld>0){
                        val=given[arrNOld-1];
                        R.push_back(val);
                    }else{
                        R.push_back(0);
                    }
                }
            }
        }
        return R;
    }

    KynematicParams_vb KynematicParams_vb::AddSmartTo(const KynematicParams_vb& obj, KynematicCharacteristics1*kynCharsNewParam, int FromN){
        KynematicParams_vb R, x1, x2;
        KynematicCharacteristics1 kynCharsNew, kynCharsNew1, kynCharsNew2;
        int QPoints1=this->GetQPoints(), QPoints2=obj.GetQPoints(), QPointsNew;
        std::vector<double>point;
        if(kynCharsNewParam!=NULL){
           kynCharsNew=(*(kynCharsNewParam));
        }else{
            kynCharsNew=this->kynChars;
        }
        QPointsNew=kynCharsNew.GetQPoints();
        kynCharsNew1=kynCharsNew;
        kynCharsNew1.SetQPoints(QPoints1);
        kynCharsNew2=kynCharsNew;
        kynCharsNew2.SetQPoints(QPointsNew);//if so
        x1=this->GetSubStr(1, &kynCharsNew1);
        x2=this->GetSubStr(FromN, &kynCharsNew2);
        R.DelAllPoints();
        for(int i=1; i<=QPoints1; i++){
            point=x1.GetPointN(i);
            R.AddPoint(point);
        }
        for(int i=1; i<=QPointsNew; i++){
            point=x2.GetPointN(i);
            R.AddPoint(point);
        }
        R.kynChars=kynCharsNew;
        return R;
    }


    KynematicParams_vb KynematicParams_vb::AddSmartTo(const KynematicParams_vb& obj, int IfNotSameStr_by1st1_by2nd2_DoFsAndDerivsAnd3_DoFsAndDerivsOr4){
        int  MaxDerivN1, QPoints1, QDoFs1, DoFsCharNum1,  MaxDerivN2, QPoints2, QDoFs2, DoFsCharNum2, MaxDerivN3, DoFsCharNum3;
        int DoFs1[]={0, 0, 0, 0, 0, 0},  DoFs2[]={0, 0, 0, 0, 0, 0}, DoFs3[]={0, 0, 0, 0, 0, 0};
        this->kynChars.Get(DoFs1, MaxDerivN1, QPoints1, QDoFs1, DoFsCharNum1);
        obj.kynChars.Get(DoFs2, MaxDerivN2, QPoints2, QDoFs2, DoFsCharNum2);
        KynematicCharacteristics1 kynChars3;
        KynematicParams_vb R;
        std::vector<double>pointToAdd, pointWas;
        if(DoFsCharNum1==DoFsCharNum2 && MaxDerivN1==MaxDerivN2){
            for(int i=1; i<=QPoints2; i++){
                pointWas=obj.GetPointN(i);
                pointToAdd=pointWas;
                this->AddPoint(pointToAdd);
            }
        }else{
            switch(IfNotSameStr_by1st1_by2nd2_DoFsAndDerivsAnd3_DoFsAndDerivsOr4){
                case 1:
                    for(int i=1; i<=6; i++){
                        DoFs3[i-1]=DoFs1[i-1];
                    }
                    MaxDerivN3=MaxDerivN1;
                break;
                case 2:
                    for(int i=1; i<=6; i++){
                        DoFs3[i-1]=DoFs2[i-1];
                    }
                    MaxDerivN3=MaxDerivN2;
                break;
                case 3:
                    for(int i=1; i<=6; i++){
                        DoFs3[i-1]=DoFs1[i-1]*DoFs2[i-1];
                    }
                    MaxDerivN3 = MaxDerivN1<=MaxDerivN2 ? MaxDerivN1 : MaxDerivN2;
                break;
                case 4:
                    for(int i=1; i<=6; i++){
                        DoFs3[i-1]= DoFs1[i-1]+DoFs2[i-1] >0 ? 1 : 0;
                    }
                    MaxDerivN3 = MaxDerivN1>=MaxDerivN2 ? MaxDerivN1 : MaxDerivN2;
                break;
            }
            kynChars3.Set(DoFs3, MaxDerivN3, QPoints1+QPoints2);
            DoFsCharNum3=kynChars3.GetDoFsCharNum();
            R.DelAllPoints();
            for(int i=1; i<=QPoints1; i++){
                pointWas=this->GetPointN(i);
                pointToAdd=this->GetRowForNewChars(pointWas, DoFsCharNum1,  MaxDerivN1, DoFsCharNum3,  MaxDerivN3);
                R.AddPoint(pointToAdd);
            }
            for(int i=1; i<=QPoints2; i++){
                pointWas=obj.GetPointN(i);
                pointToAdd=this->GetRowForNewChars(pointWas, DoFsCharNum2,  MaxDerivN2, DoFsCharNum3,  MaxDerivN3);
                R.AddPoint(pointToAdd);
            }
        }
        return R;
    }

    KynematicParams_vb KynematicParams_vb::GetSubStr(int FromN, KynematicCharacteristics1*kynCharsParam){
        int  MaxDerivNOld, QPointsOld, QDoFsOld, DoFsCharNumOld,  MaxDerivNNew, DoFsCharNumNew,QDoFsNew, QPointsNew, QPointsMin, QPointsReMax,  LastNRe, LNew;
        int DoFsOld[]={0, 0, 0, 0, 0, 0},   DoFsNew[]={0, 0, 0, 0, 0, 0};
        this->kynChars.Get(DoFsOld, MaxDerivNOld, QPointsOld, QDoFsOld, DoFsCharNumOld);
        //obj.kynChars.Get(DoFs2, MaxDerivN2, QPoints2, QDoFs2, DoFsCharNum2);
        KynematicCharacteristics1 kynCharsNew, kynCharsTmp;
        KynematicParams_vb R;
        std::vector<double>pointToAdd, pointWas, pointEmpty;
        if(kynCharsParam!=NULL){
            kynCharsNew=(*(kynCharsParam));
        }else{
           kynCharsNew=this->kynChars;
        }
        this->kynChars.Get(DoFsOld, MaxDerivNOld, QPointsOld, QDoFsOld, DoFsCharNumOld);
        kynCharsNew.Get(DoFsNew, MaxDerivNNew, QPointsNew, QDoFsNew, DoFsCharNumNew);
        QPointsReMax=QPointsOld-FromN+1;
        if(QPointsNew<0)QPointsNew=QPointsReMax;
        QPointsMin=QPointsNew <= QPointsReMax ? QPointsNew : QPointsReMax;
        LastNRe=FromN-1+QPointsMin;
        if(DoFsCharNumOld==DoFsCharNumNew && MaxDerivNOld==MaxDerivNNew){
            for(int ptN=FromN; ptN<=QPointsOld; ptN++){
                pointWas=this->GetPointN(ptN);
                pointToAdd=pointWas;
                //this->AddPoint(pointToAdd);
                R.AddPoint(pointToAdd);
            }
        }else{
            R.DelAllPoints();
            kynCharsTmp=kynCharsNew;
            kynCharsTmp.SetQPoints(0);
            //R.SetParams(&kynCharsNew);
            R.SetParams(&kynCharsTmp);
            LNew=kynCharsNew.GetArrayLengthOfPointKynematicParams();
            for(int i=1; i<=LNew; i++)pointEmpty.push_back(0);
            for(int ptN=FromN; ptN<=QPointsOld; ptN++){
                pointWas=this->GetPointN(ptN);
                pointToAdd=this->GetRowForNewChars(pointWas, DoFsCharNumOld,  MaxDerivNOld, DoFsCharNumNew,  MaxDerivNNew);
                R.AddPoint(pointToAdd);
            }
            for(int i=1; i<=QPointsNew-QPointsMin; i++){
                R.AddPoint(pointEmpty);
            }
        }
        //kynCharsNew.SetQPoints(QPointsNew);
        R.kynChars=kynCharsNew;
        return R;
    }
