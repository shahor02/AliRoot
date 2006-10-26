#ifndef ALIMUONTRACKPARAM_H
#define ALIMUONTRACKPARAM_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/*$Id$*/
// Revision of includes 07/05/2004

/// \ingroup rec
/// \class AliMUONTrackParam
/// \brief Track parameters in ALICE dimuon spectrometer
///
////////////////////////////////////////////////////
/// Track parameters in ALICE dimuon spectrometer
////////////////////////////////////////////////////

#include <TObject.h>
#include "AliMUONHitForRec.h"

class AliESDMuonTrack;
class AliMagF;

class AliMUONTrackParam : public TObject 
{
 public:
  AliMUONTrackParam(); // Constructor
  virtual ~AliMUONTrackParam(); // Destructor
  
  AliMUONTrackParam(const AliMUONTrackParam& theMUONTrackParam);
  AliMUONTrackParam& operator=(const  AliMUONTrackParam& theMUONTrackParam);

  void GetParamFrom(const AliESDMuonTrack& esdMuonTrack);
  void SetParamFor(AliESDMuonTrack& esdMuonTrack);


  // Get and Set methods for data
	/// return inverse bending momentum (GeV/c ** -1) times the charge (assumed forward motion)
  Double_t GetInverseBendingMomentum(void) const {return fInverseBendingMomentum;}
	/// set inverse bending momentum (GeV/c ** -1) times the charge (assumed forward motion)
  void     SetInverseBendingMomentum(Double_t InverseBendingMomentum) {fInverseBendingMomentum = InverseBendingMomentum;}
	/// return bending slope (cm ** -1)
  Double_t GetBendingSlope(void) const {return fBendingSlope;}
	/// set bending slope (cm ** -1)
  void     SetBendingSlope(Double_t BendingSlope) {fBendingSlope = BendingSlope;}
	/// return non bending slope (cm ** -1)
  Double_t GetNonBendingSlope(void) const {return fNonBendingSlope;}
	/// set non bending slope (cm ** -1)
  void     SetNonBendingSlope(Double_t NonBendingSlope) {fNonBendingSlope = NonBendingSlope;}
	/// return Z coordinate (cm)
  Double_t GetZ(void) const {return fZ;}
	/// set Z coordinate (cm)
  void     SetZ(Double_t Z) {fZ = Z;}
	/// return bending coordinate (cm)
  Double_t GetBendingCoor(void) const {return fBendingCoor;}
	/// set bending coordinate (cm)
  void     SetBendingCoor(Double_t BendingCoor) {fBendingCoor = BendingCoor;}
	/// return non bending coordinate (cm)
  Double_t GetNonBendingCoor(void) const {return fNonBendingCoor;}
	/// set non bending coordinate (cm)
  void     SetNonBendingCoor(Double_t NonBendingCoor) {fNonBendingCoor = NonBendingCoor;}
  void              SetTrackParam(AliMUONTrackParam& TrackParam);
  AliMUONHitForRec* GetHitForRecPtr(void) const;
	/// set pointeur to associated HitForRec if any
  void              SetHitForRecPtr(AliMUONHitForRec* HitForRec) {fHitForRecPtr = HitForRec;}
  
  Double_t Px() const;  // return px
  Double_t Py() const;  // return py
  Double_t Pz() const;  // return pz
  Double_t P()  const;  // return total momentum

	/// necessary for sorting TClonesArray of TrackHit's
  Bool_t IsSortable () const {return kTRUE;}
	/// "Compare" function for sorting
  Int_t Compare(const TObject* TrackParam) const;

  void ExtrapToZ(Double_t Z);
  void ExtrapToStation(Int_t Station, AliMUONTrackParam *TrackParam);
  void ExtrapToVertex(Double_t xVtx, Double_t yVtx, Double_t zVtx);  // extrapolation to vertex through the absorber (with true vertex) 
   void BransonCorrection(Double_t xVtx, Double_t yVtx, Double_t zVtx); // makes Branson correction with true vertex  
   // returns total momentum after energy loss correction in the absorber
  Double_t TotalMomentumEnergyLoss(Double_t thetaLimit, Double_t pTotal, Double_t theta);
  void FieldCorrection(Double_t Z); // makes simple magnetic field correction through the absorber 

  void ExtrapOneStepHelix(Double_t charge, Double_t step, 
			  Double_t *vect, Double_t *vout) const;
  void ExtrapOneStepHelix3(Double_t field, Double_t step, 
			   Double_t *vect, Double_t *vout) const;

  void ExtrapOneStepRungekutta(Double_t charge, Double_t step, 
			       Double_t* vect, Double_t* vout) const;
  
  virtual void Print(Option_t* opt="") const;
 
	/// set field map
  void SetField(const AliMagF* magField) {fkField = magField;}


 protected:
 private:
  Double_t fInverseBendingMomentum; ///< Inverse bending momentum (GeV/c ** -1) times the charge (assumed forward motion)
  Double_t fBendingSlope; ///< Bending slope (cm ** -1)
  Double_t fNonBendingSlope; ///< Non bending slope (cm ** -1)
  Double_t fZ; ///< Z coordinate (cm)
  Double_t fBendingCoor; ///< bending coordinate (cm)
  Double_t fNonBendingCoor; ///< non bending coordinate (cm)

  const AliMagF* fkField;     //!< field map

  void SetGeant3Parameters(Double_t *VGeant3, Double_t ForwardBackward);
  void GetFromGeant3Parameters(Double_t *VGeant3, Double_t Charge);

  void GetField(Double_t *Position, Double_t *Field) const;
  
  AliMUONHitForRec *fHitForRecPtr; //!< Pointer to associated HitForRec if any
  
  ClassDef(AliMUONTrackParam, 2) // Track parameters in ALICE dimuon spectrometer
    };
	
#endif
