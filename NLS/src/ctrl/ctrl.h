// -------------------------------------------------------------------------
//
// Ctrl.h - This header file contains the public data structure
//          definitions for the control class.
//
// -------------------------------------------------------------------------

#ifndef _CTRL_H_
#define _CTRL_H_

class cModel; class cLinSys;

typedef enum _updatetype
{
 STANDARD,
 MODIFIED
} eUpdate;

typedef enum _ctrltype
{
 CONSTANT,
 VARIABLE
} eCtrlType;

typedef struct _control
{
 int       CtrlEq;
 double    CtrlFactor;
 double    CtrlIniFactor;
 eCtrlType CtrlType;
 int       NumMaxIte;
 int       NumMaxStep;
 double    Tol;
 eUpdate   UpdateType;
} sControl;

// -------------------------------------------------------------------------
// Control class:

class cControl
{
 protected:
  sControl  _sControl;
  int       _iConvFlag;
  int       _iCurrIte;
  int       _iCurrStep;
  cModel   *_pcModel;
  int       _iNumEq;
  cLinSys  *_pcLinSys;
  double   *_pdReference;
  double    _dTotFactor;

 public:
           cControl    ( cModel *, sControl *, cLinSys * );
  virtual ~cControl    ( void );
  virtual  void Solver ( void ) = 0;
};

#endif
