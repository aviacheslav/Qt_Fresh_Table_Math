#ifndef KYNEMPARAM_H
#define KYNEMPARAM_H

#include"mymathlib.h"
//#include"array2dlib.h"
#include"myarraylib.h"


class DoFVals{
public:
    double displacement, velocity, acceleration;
    DoFVals(double displacement=0, double velocity=0, double acceleration=0);
    void set(std::vector<double>v, bool ifEmptyNotZeros=true);
    std::vector<double>get();
};

class KynematicCharacteristics1{
public:
    int DoFsAndDerivsCharNum, QPoints;
    //KynematicCharacteristics1();
    //KynematicCharacteristics1(int X, int Y, int Z, int FiX, int FiY, int FiZ, int MaxDerivN, int QPoints);
    KynematicCharacteristics1(int X=1, int Y=0, int Z=0, int FiX=0, int FiY=0, int FiZ=0, int MaxDerivN=0, int QPoints=1);
    KynematicCharacteristics1(int*DoFsPE, int MaxDerivN, int QPoints);
    KynematicCharacteristics1(int DoFsCharNums, int MaxDerivN, int QPoints);
    //void construct();//not used
    void SetMaxDerivN(int val);
    int GetQDoFs()const;
    int GetDoFsCharNum();
    int GetQPoints()const;
    int ArrayNByDofNDerivN(int DoFN, int DerivN);
    void DoFNAndDerivNByArrayN(int ArrayN, int&DoFN, int&DerivN);
    int GetMaxDerivN()const;
    int GetArrayLengthOfPointKynematicParams()const;
    void Set(int X, int Y, int Z, int FiX, int FiY, int FiZ, int MaxDerivN, int QPoints);
    void Set(int*DoFs, int MaxDerivN, int QPoints);
    void Set(std::vector<int>DoFs, int MaxDerivN, int QPoints);
    void Set(int DoFsCharNums, int MaxDerivN, int QPoints);
    void Set(int DoFsAndDerivsCharNum, int QPoints);
    void Get(int DoFs[], int&MaxDerivN, int&QPoints, int&QDoFs, int &DoFsCharNum) const ;
    void get(int&X, int&Y, int&Z, int&FiX, int&FiY, int&FiZ, int&MaxDerivN, int&QPoints);
    void SetDoFs(int X, int Y, int Z, int FiX, int FiY, int FiZ);
    void SetDoFs(int*DoFs);
    void SetDoFs(int DoFsCharNums);
    void SetIfDoFNExists(bool val, int N);
    void SetIfDoFXExists(bool val);
    void SetIfDoFYExists(bool val);
    void SetIfDoFZExists(bool val);
    void SetIfDoFFiXExists(bool val);
    void SetIfDoFFiYExists(bool val);
    void SetIfDoFFiZExists(bool val);
    bool GetIfDoFNExists(int DoFN);
    bool GetIfDoFXExists();
    bool GetIfDoFYExists();
    bool GetIfDoFZExists();
    bool GetIfDoFFiXExists();
    bool GetIfDoFFiYExists();
    bool GetIfDoFFiZExists();
    //void AddPoint(std::vector<double>point);
    void AddOrInsPoint();
    void DelPoint();
    void SetQPoints(int QPoints);
    QString ToString();

};

//********************************************************************************************************************************


class KynematicParams_vb
{
public:
    //std::vector<PointKynematicParams>data;
    std::vector<std::vector<double>>data;
    KynematicCharacteristics1 kynChars;
    //
    KynematicParams_vb(int  *DoFsParam, int MaxDerivN, int QPoints=1);
    KynematicParams_vb(KynematicCharacteristics1 *kynCharsParam=NULL);
    KynematicParams_vb(std::vector<std::vector<double>>data,  KynematicCharacteristics1 *kynChars=NULL);
    KynematicParams_vb(const KynematicParams_vb&obj);
    //
    ~KynematicParams_vb();
    //
    KynematicParams_vb operator = (const KynematicParams_vb&obj);
    void Assign(const KynematicParams_vb&obj);
    //
    KynematicCharacteristics1 GetKynematicCharacteristics( KynematicCharacteristics1 *kynChars);
    //
    void DelAllPoints();//ne tested
    //
    void SetMaxDerivN(int val);
    int GetMaxDerivN();
    int GetQDoFs();
    void GetParams(int&X, int&Y, int&Z, int&FiX, int&FiY, int&FiZ, int&MaxDerivN, int&QPoints);
    void GetParams(int DoFs[], int&MaxDerivN, int&QPoints);
    void GetDoFsAndDerivs(int&X, int&Y, int&Z, int&FiX, int&FiY, int&FiZ, int&MaxDerivN);
    void GetDoFsAndDerivs(int DoFs[], int&MaxDerivN);
    //void SetParams(const KynematicCharacteristics&kynChars);
    void Set(std::vector<std::vector<double>>data, KynematicCharacteristics1 *kynChars);
    void SetParams(KynematicCharacteristics1*kynChars);
    void SetParams(int X, int Y, int Z, int FiX, int FiY, int FiZ, int MaxDerivN, int QPoints=-1);
    void SetParams(int *DoFs, int MaxDerivN, int QPoints=-1);
    void SetDoFs(int X, int Y, int Z, int FiX, int FiY, int FiZ);
    void SetDoFs(int*DoFs);
    void SetDoFs(int DoFsCharNum);
    void SetDoFsAndDerivs(int*DoFs, int MaxDerivN);
    void SetDoFsAndDerivs(int DoFsCharNums, int MaxDerivN);
    //void SetDoFsAndDerivs(int DoFsAndDerivsCharNums);
    void set_x(double val, int PointN=-1);
    void set_y(double val, int PointN=-1);
    void set_z(double val, int PointN=-1);
    void set_vx(double val, int PointN=-1);
    void set_vy(double val, int PointN=-1);
    void set_vz(double val, int PointN=-1);
    void set_ax(double val, int PointN=-1);
    void set_ay(double val, int PointN=-1);
    void set_az(double val, int PointN=-1);
    void set_fix(double val, int PointN=-1);
    void set_fiy(double val, int PointN=-1);
    void set_fiz(double val, int PointN=-1);
    void set_wx(double val, int PointN=-1);
    void set_wy(double val, int PointN=-1);
    void set_wz(double val, int PointN=-1);
    void set_ex(double val, int PointN=-1);
    void set_ey(double val, int PointN=-1);
    void set_ez(double val, int PointN=-1);
    void SetValsByDoFNDerivNPointN(double val, int DoFN, int DerivN, int PointN=-1);
    double GetValsByDoFNDerivNPointN(int DoFN, int DerivN, int PointN=-1);
    void SetIfDoFXExists(int DoFN, bool val);
    void SetIfDoFXExists(bool val);
    void SetIfDoFYExists(bool val);
    void SetIfDoFZExists(bool val);
    void SetIfDoFFiXExists(bool val);
    void SetIfDoFFiYExists(bool val);
    void SetIfDoFFiZExists(bool val);
    bool GetIfDoFNExists(int DoFN);
    bool GetIfDoFXExists();
    bool GetIfDoFYExists();
    bool GetIfDoFZExists();
    bool GetIfDoFFiXExists();
    bool GetIfDoFFiYExists();
    bool GetIfDoFFiZExists();
    void ExcludeDoFN(int ExcludedDoFN);
    void ExcludeDoFX();
    void ExcludeDoFY();
    void ExcludeDoFZ();
    void ExcludeDoFFiX();
    void ExcludeDoFFiY();
    void ExcludeDoFFiZ();
    void AddDoFN(int AddedDoFN);
    void AddDoFX();
    void AddDoFY();
    void AddDoFZ();
    void AddDoFFiX();
    void AddDoFFiY();
    void AddDoFFiZ();
    double get_x(int PointN=-1);
    double get_y(int PointN=-1);
    double get_z(int PointN=-1);
    double get_vx(int PointN=-1);
    double get_vy(int PointN=-1);
    double get_vz(int PointN=-1);
    double get_ax(int PointN=-1);
    double get_ay(int PointN=-1);
    double get_az(int PointN=-1);
    double get_fiX(int PointN=-1);
    double get_fiY(int PointN=-1);
    double get_fiZ(int PointN=-1);
    double get_wx(int PointN=-1);
    double get_wy(int PointN=-1);
    double get_wz(int PointN=-1);
    double get_ex(int PointN=-1);
    double get_ey(int PointN=-1);
    double get_ez(int PointN=-1);
    //
    std::vector<double>Get_all_X();
    std::vector<double>Get_all_Y();
    std::vector<double>Get_all_Z();
    std::vector<double>Get_all_Vx();
    std::vector<double>Get_all_Vy();
    std::vector<double>Get_all_Vz();
    std::vector<double>Get_all_ax();
    std::vector<double>Get_all_ay();
    std::vector<double>Get_all_az();
    std::vector<double>Get_all_FiX();
    std::vector<double>Get_all_FiY();
    std::vector<double>Get_all_FiZ();
    std::vector<double>Get_all_wx();
    std::vector<double>Get_all_wy();
    std::vector<double>Get_all_wz();
    std::vector<double>Get_all_ex();
    std::vector<double>Get_all_ey();
    std::vector<double>Get_all_ez();
    //
    void Set_all_X( std::vector<double>data );
    void Set_all_Y( std::vector<double>data );
    void Set_all_Z( std::vector<double>data );
    void Set_all_Vx( std::vector<double>data );
    void Set_all_Vy( std::vector<double>data );
    void Set_all_Vz( std::vector<double>data );
    void Set_all_ax( std::vector<double>data );
    void Set_all_ay( std::vector<double>data );
    void Set_all_az( std::vector<double>data );
    void Set_all_FiX( std::vector<double>data );
    void Set_all_FiY( std::vector<double>data );
    void Set_all_FiZ( std::vector<double>data );
    void Set_all_wx( std::vector<double>data );
    void Set_all_wy( std::vector<double>data );
    void Set_all_wz( std::vector<double>data );
    void Set_all_ex( std::vector<double>data );
    void Set_all_ey( std::vector<double>data );
    void Set_all_ez( std::vector<double>data );
    //
    void SetPointN(int N, double*data);
    void SetPointN(int N, std::vector<double>data);
    std::vector<double>GetPointN(int N) const;
    //
    int GetQPoints()const;
    int GetArrayLengthOfPointKynematicParams()const;
    //
    //void AddPoint(double*data=NULL);
    //void AddPoint(double data[]);

    void SetPoint(int N, double*data=NULL);
    //void SetPoint(int N, double data[]);
    void SetPoint(int N, std::vector<double>data);
    void SetPoint(int N, double*XDoF,double*YDoF,double*ZDoF,double*FiXDoF,double*FiYDoF,double*FiZDoF);
    //
    std::vector<double>GetRowVectorOfPointDataByDoFs(double*X, double*Y, double*Z, double*FiX, double*FiY, double*FiZ);
    std::vector<double>GetRowVectorOfPointDataByDoFs( std::vector<double>X, std::vector<double>Y, std::vector<double>Z, std::vector<double>FiX, std::vector<double>FiY, std::vector<double>FiZ);
    //
    std::vector<double>GetRowVectorOfPointDataByDoFsAsPtrs(double*X, double*Y, double*Z, double*FiX, double*FiY, double*FiZ);
    std::vector<double>GetRowVectorOfPointDataByDoFsAsVecs(std::vector<double>X, std::vector<double>Y, std::vector<double>Z, std::vector<double>FiX, std::vector<double>FiY, std::vector<double>FiZ);
    std::vector<double>GetRowVectorOfPointDataByDoFsAsPtrs(DoFVals*X=NULL, DoFVals*Y=NULL, DoFVals*Z=NULL, DoFVals*FiX=NULL, DoFVals*FiY=NULL, DoFVals*FiZ=NULL);
    //
    void AddPoint(double*data, int LGiven=0);
    void AddPoint(double*data, const KynematicCharacteristics1 kynCharsParam);
    void AddPoint(std::vector<double>data);//ne tested
    void AddPoint(double*X, double*Y, double*Z, double*FiX, double*FiY, double*FiZ);//ne tested
    void AddPoint(DoFVals*X=NULL, DoFVals*Y=NULL, DoFVals*Z=NULL, DoFVals*FiX=NULL, DoFVals*FiY=NULL, DoFVals*FiZ=NULL);
    void AddPoint(std::vector<double>X, std::vector<double>Y, std::vector<double>Z, std::vector<double>FiX, std::vector<double>FiY, std::vector<double>FiZ);//ne tested
    //
    void InsPoint(int N, std::vector<double>data);//ne tested
    void InsPoint(int N, double*given, int LGiven=0);
    void InsPoint(int N, double*X, double*Y, double*Z, double*FiX, double*FiY, double*FiZ);//ne tested
    void InsPoint(int N, DoFVals*X=NULL, DoFVals*Y=NULL, DoFVals*Z=NULL, DoFVals*FiX=NULL, DoFVals*FiY=NULL, DoFVals*FiZ=NULL);
    void InsPoint(int N, std::vector<double>X, std::vector<double>Y, std::vector<double>Z, std::vector<double>FiX, std::vector<double>FiY, std::vector<double>FiZ);//ne tested
    //
    void DelPoint(int N=-1);
    //
    QString PointNToString(int N, bool showAll=false);
    QString DoFNDerivNToString(int DoFN, int DerivN, QString delim="; ");

    QString get_all_x_as_String(QString delim="; ");
    QString get_all_vx_as_String(QString delim="; ");
    QString get_all_ax_as_String(QString delim="; ");
    QString get_all_fix_as_String(QString delim="; ");
    QString get_all_wx_as_String(QString delim="; ");
    QString get_all_ex_as_String(QString delim="; ");

    QString get_all_y_as_String(QString delim="; ");
    QString get_all_vy_as_String(QString delim="; ");
    QString get_all_ay_as_String(QString delim="; ");
    QString get_all_fiy_as_String(QString delim="; ");
    QString get_all_wy_as_String(QString delim="; ");
    QString get_all_ey_as_String(QString delim="; ");

    QString get_all_z_as_String(QString delim="; ");
    QString get_all_vz_as_String(QString delim="; ");
    QString get_all_az_as_String(QString delim="; ");
    QString get_all_fiz_as_String(QString delim="; ");
    QString get_all_wz_as_String(QString delim="; ");
    QString get_all_ez_as_String(QString delim="; ");

    QString ParamsToString();

    void ShowConsole();

    std::vector<double>GetRowForNewChars(std::vector<double>given, int DoFsCharNumOld, int MaxDerivNOld, int DoFsCharNumNew, int MaxDerivNNew);
    
    //void AddSmart(const KynematicParams_vb& obj, int IfNotSameStr_by1st1_by2nd2_DoFsAndDerivsAnd3_DoFsAndDerivsOr4=1);
    KynematicParams_vb AddSmartTo(const KynematicParams_vb& obj, KynematicCharacteristics1*kynCharsNewParam=NULL, int FromN=1);//ne tested
    KynematicParams_vb AddSmartTo(const KynematicParams_vb& obj, int IfNotSameStr_by1st1_by2nd2_DoFsAndDerivsAnd3_DoFsAndDerivsOr4=1);//ne tested
    KynematicParams_vb GetSubStr(int FromN=1, KynematicCharacteristics1*kynChars=NULL);//ne tested
};

#endif // KYNEMPARAM_H
