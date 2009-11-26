#ifndef ALITRDDIGITSPARAM_H
#define ALITRDDIGITSPARAM_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// Class containing parameters for digits                                    //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

/* $Id: AliTRDdigitsParam.h 34070 2009-08-04 15:34:53Z cblume $ */

#include "TObject.h"

class AliTRDdigitsParam : public TObject
{

  public:
    
    AliTRDdigitsParam();
    AliTRDdigitsParam(const AliTRDdigitsParam &p);   
    AliTRDdigitsParam &operator=(const AliTRDdigitsParam &p); 
    virtual        ~AliTRDdigitsParam();

    virtual void    Copy(TObject &p) const;

            void    SetNTimeBins(Int_t ntb)                    { fNTimeBins            = ntb;      }     
            void    SetPretiggerPhase(Int_t det, UInt_t phase) { fPretriggerPhase[det] = phase;    }
            void    SetADCbaseline(Int_t baseline)             { fADCbaseline          = baseline; }

            Int_t   GetNTimeBins() const                       { return fNTimeBins;                }
            UInt_t  GetPretriggerPhase(Int_t det) const        { return fPretriggerPhase[det];     }
            Int_t   GetADCbaseline() const                     { return fADCbaseline;              }

  protected:

            Int_t   fNTimeBins;            //  Number of timebins
            UInt_t  fPretriggerPhase[540]; //  Pretrigger phase for each detector
            Int_t   fADCbaseline;          //  ADC baseline, given in ADC channels

    ClassDef(AliTRDdigitsParam,4)          //  The parameters for digits

};
#endif
