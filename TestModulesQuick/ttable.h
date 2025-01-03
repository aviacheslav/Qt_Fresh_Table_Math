#ifndef TTABLE_H
#define TTABLE_H

#include"mylib.h"

const int DataCellTypeN_Double=1;
const int DataCellTypeN_Float=2;
const int DataCellTypeN_Int=3;
const int DataCellTypeN_Int_UniqueValKeeper=10;
const int DataCellTypeN_Int_HeaderItemN=11;
const int DataCellTypeN_Bool=4;
const int DataCellTypeN_String=5;
const int DataCellTypeN_StdString=6;
const int DataCellTypeN_ComboboxOrMemo=7;
const int DataCellTypeN_DBFieldData=8;
const int DataCellTypeN_DBTableData=9;
//
class Table2DIntVector{
public:
    int d1, d2;
    Table2DIntVector(int d1=0, int d2=0);
};

class TDBTableData{
public:
    TDBTableData();
    TDBTableData(const TDBTableData&obj);
    ~TDBTableData();//no need
    //
    void Assign(const TDBTableData&obj);
    TDBTableData operator = (const TDBTableData&obj);
    //
    void SetDBTableName(QString name);
    void SetDBFileFullName(QString FullName);
    void SetDBTypeName(QString name);
    void SetDBTypeN(int TypeN);
    void SetDBName(QString name);
    //
    QString GetDBTableName(QString name);
    QString GetDBFileFullName(QString FullName);
    QString GetDBTypeName(QString name);
    int GetDBTypeN();
    QString GetDBName(QString name);
    //
    QString ToStringFull();
    std::vector<QString>GetAsStringArray();
    void SetFromStringArray(QString*arr, int Q, int N1);
};

class TDBTableHeader{
public:
    QString DBTableName;
    TDBTableData*TBTableData;
    //
    TDBTableHeader();
    TDBTableHeader(const TDBTableHeader&obj);
    ~TDBTableHeader();//necessary
    //
    void Assign(const TDBTableHeader&obj);
    TDBTableHeader operator = (const TDBTableHeader&obj);
    //
    QString ToStringFull();
    std::vector<QString>GetAsStringArray();
    void SetFromStringArray(QString*arr, int Q, int N1);
    //
    void SetTBTablName(QString name);
    QString GetTBTableName();
    void SetTBTableData(TDBTableData*TBTableData);
    TDBTableData* GetTBTableData();
    //
    void SetDBTableName(QString name);
    void SetDBFileFullName(QString FullName);
    void SetDBTypeName(QString name);
    void SetDBTypeN(int TypeN);
    void SetDBName(QString name);
    //
    QString GetDBTableName(QString name);
    QString GetDBFileFullName(QString FullName);
    QString GetDBTypeName(QString name);
    int GetDBTypeN();
    QString GetDBName(QString name);
};

class TConnectedTableData{
public:
    QString ConnectedTableName;
    QString ConnectedTableContentFieldName0;
    QString ConnectedTableKeyFieldName;
    //
    TConnectedTableData();
    TConnectedTableData(const TConnectedTableData&obj);
    ~TConnectedTableData();
    //
    void Assign(const TConnectedTableData&obj);
    TConnectedTableData operator = (const TConnectedTableData&obj);
    //
    QString ToStringFull();
    std::vector<QString>GetAsStringArray();
    void SetFromStringArray(QString*arr, int Q, int N1);
};

class TDBFieldData{
public:
    //QString ColHeaderName()=0;
    //QString DBFieldHeaderName()=0;
    QString DBFieldTypeName;
    int DBFieldTypeN;
    int DBFieldWidth;
    bool IfAutoIncement;
    bool IfKeyField;
    bool IfFieldValIsObligatory;
    int DBFieldCharacteristics;
    //
    TDBFieldData();
    TDBFieldData(const TDBFieldData&obj);
    ~TDBFieldData();
    //
    void Assign(const TDBFieldData&obj);
    TDBFieldData operator = (const TDBFieldData&obj);
    //
    QString ToStringFull();
    std::vector<QString>GetAsStringArray();
    void SetFromStringArray(QString*arr, int Q, int N1);
    //
    void SetColHeaderName(QString name);
    void SetDBFieldHeaderName(QString name);
    void SetDBFieldTypeName(QString name);
    void SetDBFieldTypeN(int TypeN);
    void SetDBFieldWidth(int width);
    void SetIfAutoIncement(bool IsAutoIncement);
    void SetIfKeyField(bool IsKeyField);
    void SetIfFieldValIsObligatory(bool ValIsObligatory);
    void SetDBFieldCharacteristics(int DBFieldCharacteristics);
    //
    QString GetColHeaderName();
    QString GetDBFieldHeaderName();
    QString GetDBFieldTypeName();
    int GetDBFieldTypeN();
    int GetDBFieldWidth();
    bool GetIfAutoIncement();
    bool GetIfKeyField();
    bool GetIfFieldValIsObligatory();
    int GetDBFieldCharacteristics();
};
class TDBFieldHeader{

};

class TDataCell{
public:
    virtual void SetByValDouble(double val, int whereToSetN=0)=0;
    virtual void SetByValFloat(float val, int whereToSetN=0)=0;
    virtual void SetByValInt(int val, int whereToSetN=0)=0;
    virtual void SetByValBool(bool val, int whereToSetN=0)=0;
    virtual void SetByValString(QString val, int whereToSetN=0)=0;
    virtual void SetByValStdString(std::string val, int whereToSetN=0)=0;
    //
    virtual double GetDoubleVal(int N=0)=0;
    virtual float GetFloatVal(int N=0)=0;
    virtual int GetIntVal(int N=0)=0;
    virtual bool GetBoolVal(int N=0)=0;
    virtual QString GetStringVal(int N=0);
    virtual std::string GetStdStringVal(int N=0);
    //
    //virtual void SetActiveN(int NewActiveN=0)=0;
    virtual void SetActiveItemN(int NewActiveItemN)=0;
    virtual void SetName1(QString Name1)=0;
    virtual void SetName2(QString Name2)=0;
    virtual void SetItem(int itemN, QString item)=0;
    virtual void AddItem(QString item)=0;
    virtual void InsItemN(int itemN, QString item)=0;
    virtual void DelItemN(int itemN)=0;
    virtual std::vector<QString>GetItems()=0;
    //
    virtual int GetTypeN()=0;
    virtual int GetActiveItemN()=0;
    //virtual int GetActiveN()=0;
    //
    virtual int GetItemsListLength()=0;
    virtual QString GetName1()=0;
    virtual QString GetName2()=0;
    virtual QString GetActiveItem()=0;
    virtual QString GetItemN(int itemN)=0;
    //-------------------------------------------------------------
    virtual QString ToStringFull()=0;
    //
    virtual std::vector<QString>GetAsStringArray()=0;
    //
    virtual void SetFromStringArray(QString*arr, int Q, int N1)=0;
    //-------------------------------------------------------------
    virtual void SetColHeaderName(QString name)=0;
    virtual void SetDBFieldHeaderName(QString name)=0;
    virtual void SetDBFieldTypeName(QString name)=0;
    virtual void SetDBFieldTypeN(int TypeN)=0;
    virtual void SetDBFieldWidth(int width)=0;
    virtual void SetIfAutoIncement(bool IsAutoIncement)=0;
    virtual void SetIfKeyField(bool IsKeyField)=0;
    virtual void SetIfFieldValIsObligatory(bool ValIsObligatory)=0;
    virtual void SetDBFieldCharacteristics(int DBFieldCharacteristics)=0;
    //
    virtual QString GetColHeaderName()=0;
    virtual QString GetDBFieldHeaderName()=0;
    virtual QString GetDBFieldTypeName()=0;
    virtual int GetDBFieldTypeN()=0;
    virtual int GetDBFieldWidth()=0;
    virtual bool GetIfAutoIncement()=0;
    virtual bool GetIfKeyField()=0;
    virtual bool GetIfFieldValIsObligatory()=0;
    virtual int GetDBFieldCharacteristics()=0;
    //-----------------------------------------------------------------
    virtual void SetDBTableName(QString name)=0;
    virtual void SetDBFileFullName(QString FullName)=0;
    virtual void SetDBTypeName(QString name)=0;
    virtual void SetDBTypeN(int TypeN)=0;
    virtual void SetDBName(QString name)=0;
    //
    virtual QString GetDBTableName(QString name)=0;
    virtual QString GetDBFileFullName(QString FullName)=0;
    virtual QString GetDBTypeName(QString name)=0;
    virtual int GetDBTypeN()=0;
    virtual QString GetDBName(QString name)=0;
    //----------------------------------------------------------------
    virtual void SetConnectedTableName(QString name)=0;
    virtual void SetConnectedTableContentFieldName(QString name)=0;
    virtual void SetConnectedTableKeyFieldName(QString name)=0;
    //
    virtual QString GetConnectedTableName(QString name)=0;
    virtual QString GetConnectedTableContentFieldName(QString name)=0;
    virtual QString GetConnectedTableKeyFieldName(QString name)=0;
};//TDataCellAbstr

class TDEataCell_Double{
double val;
public:
    virtual void SetByValDouble(double val, int whereToSetN=0);
    virtual void SetByValFloat(float val, int whereToSetN=0);
    virtual void SetByValInt(int val, int whereToSetN=0);
    virtual void SetByValBool(bool val, int whereToSetN=0);
    virtual void SetByValString(QString val, int whereToSetN=0);
    virtual void SetByValStdString(std::string val, int whereToSetN=0);
    //
    virtual double GetDoubleVal(int whereToGetN=0);
    virtual float GetFloatVal(int whereToGetN=0);
    virtual int GetIntVal(int whereToGetN=0);
    virtual bool GetBoolVal(int whereToGetN=0);
    virtual QString GetStringVal(int whereToGetN=0);
    virtual std::string GetStdStringVal(int whereToGetN=0);
    //
    //virtual void SetActiveN(int NewActiveN);
    virtual void SetActiveItemN(int NewActiveItemN);
    virtual void SetName1(QString Name1);
    virtual void SetName2(QString Name2);
    virtual void SetItem(int itemN, QString item);
    virtual void AddItem(QString item);
    virtual void InsItemN(int itemN, QString item);
    virtual void DelItemN(int itemN);
    virtual std::vector<QString>GetItems();
    //
    virtual int GetTypeN();
    virtual int GetActiveItemN();
    //virtual int GetActiveN();
    //
    virtual int GetItemsListLength();
    virtual QString GetName1();
    virtual QString GetName2();
    virtual QString GetActiveItem();
    virtual QString GetItemN(int itemN);
    //-------------------------------------------------------------
    virtual QString ToStringFull();
    //
    virtual std::vector<QString>GetAsStringArray();
    //
    virtual void SetFromStringArray(QString*arr, int Q, int N1);
    //-------------------------------------------------------------
    virtual void SetColHeaderName(QString name);
    virtual void SetDBFieldHeaderName(QString name);
    virtual void SetDBFieldTypeName(QString name);
    virtual void SetDBFieldTypeN(int TypeN);
    virtual void SetDBFieldWidth(int width);
    virtual void SetIfAutoIncement(bool IsAutoIncement);
    virtual void SetIfKeyField(bool IsKeyField);
    virtual void SetIfFieldValIsObligatory(bool ValIsObligatory);
    virtual void SetDBFieldCharacteristics(int DBFieldCharacteristics);
    //
    virtual QString GetColHeaderName();
    virtual QString GetDBFieldHeaderName();
    virtual QString GetDBFieldTypeName();
    virtual int GetDBFieldTypeN();
    virtual int GetDBFieldWidth();
    virtual bool GetIfAutoIncement();
    virtual bool GetIfKeyField();
    virtual bool GetIfFieldValIsObligatory();
    virtual int GetDBFieldCharacteristics();
    //-----------------------------------------------------------------
    virtual void SetDBTableName(QString name);
    virtual void SetDBFileFullName(QString FullName);
    virtual void SetDBTypeName(QString name);
    virtual void SetDBTypeN(int TypeN);
    virtual void SetDBName(QString name);
    //
    virtual QString GetDBTableName(QString name);
    virtual QString GetDBFileFullName(QString FullName);
    virtual QString GetDBTypeName(QString name);
    virtual int GetDBTypeN();
    virtual QString GetDBName(QString name);
    //----------------------------------------------------------------
    virtual void SetConnectedTableName(QString name);
    virtual void SetConnectedTableContentFieldName(QString name);
    virtual void SetConnectedTableKeyFieldName(QString name);
    //
    virtual QString GetConnectedTableName(QString name);
    virtual QString GetConnectedTableContentFieldName(QString name);
    virtual QString GetConnectedTableKeyFieldName(QString name);
};//TDEataCell_Double

class TDEataCell_String{
QString val;
public:
    virtual void SetByValDouble(double val, int whereToSetN=0);
    virtual void SetByValFloat(float val, int whereToSetN=0);
    virtual void SetByValInt(int val, int whereToSetN=0);
    virtual void SetByValBool(bool val, int whereToSetN=0);
    virtual void SetByValString(QString val, int whereToSetN=0);
    virtual void SetByValStdString(std::string val, int whereToSetN=0);
    //
    virtual double GetDoubleVal(int whereToGetN=0);
    virtual float GetFloatVal(int whereToGetN=0);
    virtual int GetIntVal(int whereToGetN=0);
    virtual bool GetBoolVal(int whereToGetN=0);
    virtual QString GetStringVal(int whereToGetN=0);
    virtual std::string GetStdStringVal(int whereToGetN=0);
    //
    //virtual void SetActiveN(int NewActiveN);
    virtual void SetActiveItemN(int NewActiveItemN);
    virtual void SetName1(QString Name1);
    virtual void SetName2(QString Name2);
    virtual void SetItem(int itemN, QString item);
    virtual void AddItem(QString item);
    virtual void InsItemN(int itemN, QString item);
    virtual void DelItemN(int itemN);
    virtual std::vector<QString>GetItems();
    //
    virtual int GetTypeN();
    virtual int GetActiveItemN();
    //virtual int GetActiveN();
    //
    virtual int GetItemsListLength();
    virtual QString GetName1();
    virtual QString GetName2();
    virtual QString GetActiveItem();
    virtual QString GetItemN(int itemN);
    //-------------------------------------------------------------
    virtual QString ToStringFull();
    //
    virtual std::vector<QString>GetAsStringArray();
    //
    virtual void SetFromStringArray(QString*arr, int Q, int N1);
    //-------------------------------------------------------------
    virtual void SetColHeaderName(QString name);
    virtual void SetDBFieldHeaderName(QString name);
    virtual void SetDBFieldTypeName(QString name);
    virtual void SetDBFieldTypeN(int TypeN);
    virtual void SetDBFieldWidth(int width);
    virtual void SetIfAutoIncement(bool IsAutoIncement);
    virtual void SetIfKeyField(bool IsKeyField);
    virtual void SetIfFieldValIsObligatory(bool ValIsObligatory);
    virtual void SetDBFieldCharacteristics(int DBFieldCharacteristics);
    //
    virtual QString GetColHeaderName();
    virtual QString GetDBFieldHeaderName();
    virtual QString GetDBFieldTypeName();
    virtual int GetDBFieldTypeN();
    virtual int GetDBFieldWidth();
    virtual bool GetIfAutoIncement();
    virtual bool GetIfKeyField();
    virtual bool GetIfFieldValIsObligatory();
    virtual int GetDBFieldCharacteristics();
    //-----------------------------------------------------------------
    virtual void SetDBTableName(QString name);
    virtual void SetDBFileFullName(QString FullName);
    virtual void SetDBTypeName(QString name);
    virtual void SetDBTypeN(int TypeN);
    virtual void SetDBName(QString name);
    //
    virtual QString GetDBTableName(QString name);
    virtual QString GetDBFileFullName(QString FullName);
    virtual QString GetDBTypeName(QString name);
    virtual int GetDBTypeN();
    virtual QString GetDBName(QString name);
    //----------------------------------------------------------------
    virtual void SetConnectedTableName(QString name);
    virtual void SetConnectedTableContentFieldName(QString name);
    virtual void SetConnectedTableKeyFieldName(QString name);
    //
    virtual QString GetConnectedTableName(QString name);
    virtual QString GetConnectedTableContentFieldName(QString name);
    virtual QString GetConnectedTableKeyFieldName(QString name);
};//TDEataCell_String

class TDEataCell_ComboboxOrMemo{
std::vector<QString>items;
int ActiveN;
public:
    virtual void SetByValDouble(double val, int whereToSetN=0);
    virtual void SetByValFloat(float val, int whereToSetN=0);
    virtual void SetByValInt(int val, int whereToSetN=0);
    virtual void SetByValBool(bool val, int whereToSetN=0);
    virtual void SetByValString(QString val, int whereToSetN=0);
    virtual void SetByValStdString(std::string val, int whereToSetN=0);
    //
    virtual double GetDoubleVal(int whereToGetN=0);
    virtual float GetFloatVal(int whereToGetN=0);
    virtual int GetIntVal(int whereToGetN=0);
    virtual bool GetBoolVal(int whereToGetN=0);
    virtual QString GetStringVal(int whereToGetN=0);
    virtual std::string GetStdStringVal(int whereToGetN=0);
    //
    //virtual void SetActiveN(int NewActiveN);
    virtual void SetActiveItemN(int NewActiveItemN);
    virtual void SetName1(QString Name1);
    virtual void SetName2(QString Name2);
    virtual void SetItem(int itemN, QString item);
    virtual void AddItem(QString item);
    virtual void InsItemN(int itemN, QString item);
    virtual void DelItemN(int itemN);
    virtual std::vector<QString>GetItems();
    //
    virtual int GetTypeN();
    virtual int GetActiveItemN();
    //virtual int GetActiveN();
    //
    virtual int GetItemsListLength();
    virtual QString GetName1();
    virtual QString GetName2();
    virtual QString GetActiveItem();
    virtual QString GetItemN(int itemN);
    //-------------------------------------------------------------
    virtual QString ToStringFull();
    //
    virtual std::vector<QString>GetAsStringArray();
    //
    virtual void SetFromStringArray(QString*arr, int Q, int N1);
    //-------------------------------------------------------------
    virtual void SetColHeaderName(QString name);
    virtual void SetDBFieldHeaderName(QString name);
    virtual void SetDBFieldTypeName(QString name);
    virtual void SetDBFieldTypeN(int TypeN);
    virtual void SetDBFieldWidth(int width);
    virtual void SetIfAutoIncement(bool IsAutoIncement);
    virtual void SetIfKeyField(bool IsKeyField);
    virtual void SetIfFieldValIsObligatory(bool ValIsObligatory);
    virtual void SetDBFieldCharacteristics(int DBFieldCharacteristics);
    //
    virtual QString GetColHeaderName();
    virtual QString GetDBFieldHeaderName();
    virtual QString GetDBFieldTypeName();
    virtual int GetDBFieldTypeN();
    virtual int GetDBFieldWidth();
    virtual bool GetIfAutoIncement();
    virtual bool GetIfKeyField();
    virtual bool GetIfFieldValIsObligatory();
    virtual int GetDBFieldCharacteristics();
    //-----------------------------------------------------------------
    virtual void SetDBTableName(QString name);
    virtual void SetDBFileFullName(QString FullName);
    virtual void SetDBTypeName(QString name);
    virtual void SetDBTypeN(int TypeN);
    virtual void SetDBName(QString name);
    //
    virtual QString GetDBTableName(QString name);
    virtual QString GetDBFileFullName(QString FullName);
    virtual QString GetDBTypeName(QString name);
    virtual int GetDBTypeN();
    virtual QString GetDBName(QString name);
    //----------------------------------------------------------------
    virtual void SetConnectedTableName(QString name);
    virtual void SetConnectedTableContentFieldName(QString name);
    virtual void SetConnectedTableKeyFieldName(QString name);
    //
    virtual QString GetConnectedTableName(QString name);
    virtual QString GetConnectedTableContentFieldName(QString name);
    virtual QString GetConnectedTableKeyFieldName(QString name);
};//TDEataCell_ComboboxOrMemo


/*
    virtual void SetValAndTypeDouble(double val)=0;
    virtual void SetValAndTypeFloat(float val)=0;
    virtual void SetValAndTypeInt(int val)=0;
    virtual void SetValAndTypeBool(bool val)=0;
    virtual void SetValAndTypeStr(QString val)=0;
    virtual void SetComboboxOrMemo(std::vector<QString>items, int ActiveItemN=1)=0;
*/

class DataCell{
TDataCell*cell;
public:
    DataCell();
    DataCell(double val);
    DataCell(float val);
    DataCell(int val);
    DataCell(bool val);
    DataCell(QString val);
    DataCell(std::string val);
    DataCell(std::vector<QString>items, int ActiveN=1);
    DataCell(QString ColNameToDisplay, std::vector<QString>items, TDBFieldData*TDBFieldData=NULL, TConnectedTableData*ConnectedTableData=NULL);
    DataCell(QString TableNameToDisplay, TDBFieldHeader*DBFieldHeader=NULL);
    DataCell(const DataCell&obj);
    ~DataCell();
    //
    void Assign(const DataCell&obj);
    DataCell operator = (const DataCell&obj);
    //
    DataCell operator = (double val);
    DataCell operator = (float val);
    DataCell operator = (int val);
    DataCell operator = (bool val);
    DataCell operator = (QString val);
    DataCell operator = (std::string val);
    DataCell operator = (std::vector<QString>items);
    //
    void SetValAndTypeDouble(double val);
    void SetValAndTypeFloat(float val);
    void SetValAndTypeInt(int val);
    void SetValAndTypeBool(bool val);
    void SetValAndTypeString(QString val);
    void SetValAndTypeStdString(std::string val);
    void SetValAndTypeComboboxOrMemo(std::vector<QString>items, int ActiveN=1);
    void SetValAndTypeDBFieldData(QString ColNameToDisplay, std::vector<QString>items, TDBFieldData*TDBFieldData=NULL, TConnectedTableData*ConnectedTableData=NULL);
    void SetValAndTypeDBTableData(QString TableNameToDisplay, TDBFieldHeader*DBFieldHeader=NULL);
    //
    void SetByValDouble(double val, int whereToSetN=0);
    void SetByValFloat(float val, int whereToSetN=0);
    void SetByValInt(int val, int whereToSetN=0);
    void SetByValBool(bool val, int whereToSetN=0);
    void SetByValString(QString val, int whereToSetN=0);
    void SetByValStdString(std::string val, int whereToSetN=0);
    //
    double GetDoubleVal(int whereToGetN=0);
    float GetFloatVal(int whereToGetN=0);
    int GetIntVal(int whereToGetN=0);
    bool GetBoolVal(int whereToGetN=0);
    QString GetStringVal(int whereToGetN=0);
    std::string GetStdStringVal(int whereToGetN=0);
    //
    //void SetActiveN(int NewActiveN);
    void SetActiveItemN(int NewActiveItemN);
    void SetName1(QString Name1);
    void SetName2(QString Name2);
    void SetItem(int itemN, QString item);
    void AddItem(QString item);
    void InsItemN(int itemN, QString item);
    void DelItemN(int itemN);
    std::vector<QString>GetItems();
    //
    int GetTypeN();
    int GetActiveItemN();
    //int GetActiveN();
    //
    int GetItemsListLength();
    QString GetName1();
    QString GetName2();
    QString GetActiveItem();
    QString GetItemN(int itemN);
    //-------------------------------------------------------------
    QString ToStringFull();
    //
    std::vector<QString>GetAsStringArray();
    //
    void SetFromStringArray(QString*arr, int Q, int N1);
    //-------------------------------------------------------------
    void SetColHeaderName(QString name);
    void SetDBFieldHeaderName(QString name);
    void SetDBFieldTypeName(QString name);
    void SetDBFieldTypeN(int TypeN);
    void SetDBFieldWidth(int width);
    void SetIfAutoIncement(bool IsAutoIncement);
    void SetIfKeyField(bool IsKeyField);
    void SetIfFieldValIsObligatory(bool ValIsObligatory);
    void SetDBFieldCharacteristics(int DBFieldCharacteristics);
    //
    QString GetColHeaderName();
    QString GetDBFieldHeaderName();
    QString GetDBFieldTypeName();
    int GetDBFieldTypeN();
    int GetDBFieldWidth();
    bool GetIfAutoIncement();
    bool GetIfKeyField();
    bool GetIfFieldValIsObligatory();
    int GetDBFieldCharacteristics();
    //-----------------------------------------------------------------
    void SetDBTableName(QString name);
    void SetDBFileFullName(QString FullName);
    void SetDBTypeName(QString name);
    void SetDBTypeN(int TypeN);
    void SetDBName(QString name);
    //
    QString GetDBTableName(QString name);
    QString GetDBFileFullName(QString FullName);
    QString GetDBTypeName(QString name);
    int GetDBTypeN();
    QString GetDBName(QString name);
    //----------------------------------------------------------------
    void SetConnectedTableName(QString name);
    void SetConnectedTableContentFieldName(QString name);
    void SetConnectedTableKeyFieldName(QString name);
    //
    QString GetConnectedTableName(QString name);
    QString GetConnectedTableContentFieldName(QString name);
    QString GetConnectedTableKeyFieldName(QString name);
};//TDEataCell_ComboboxOrMemo

class TTable
{
public:
    TTable();
};

#endif // TTABLE_H
