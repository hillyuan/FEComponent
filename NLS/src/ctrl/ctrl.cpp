// --------------------------------------------------------------------------
//
// Ctrl.cpp - Control class.
//
// --------------------------------------------------------------------------

#include "ctrl.h"
#include "mempack/mempack.h"
#include "model/model.h"
#include "linsys/linearsystem.h"

// --------------------------------------------------------------------------
// Public functions:

// ============================= cControl ===================================

cControl :: cControl ( cModel *pcModel, sControl *psControl, cLinSys *pcLinSys )
{
 _pcModel  = pcModel;
 _sControl = *psControl;
 _pcLinSys  = pcLinSys;

 _iNumEq      = pcModel->NumEq( );
 _pdReference = (double *) MemAlloc( _iNumEq, sizeof(double) );
 pcModel->Reference( _pdReference );
 _pcLinSys->Init(pcModel);
}

// ============================= ~cControl ==================================

cControl :: ~cControl ( void )
{
 MemFree( _pdReference );
}

