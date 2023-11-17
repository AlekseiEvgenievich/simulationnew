// **************************************************************
//
//      ------------ CrystalEyeAnticoincidenceHit  ------
//           by Libo Wu (20 Sep 2021)
//           Email: libo.wu@gssi.it
//
// **************************************************************

#ifndef CrystalEyeAnticoincidenceHit_h
#define CrystalEyeAnticoincidenceHit_h 1

#include "G4VHit.hh"
#include "G4THitsCollection.hh"
#include "G4Allocator.hh"
#include "G4ThreeVector.hh"

class CrystalEyeAnticoincidenceHit : public G4VHit
{
		public:
				CrystalEyeAnticoincidenceHit();
				virtual ~CrystalEyeAnticoincidenceHit();
				CrystalEyeAnticoincidenceHit(const CrystalEyeAnticoincidenceHit&);

				const CrystalEyeAnticoincidenceHit& operator=(const CrystalEyeAnticoincidenceHit&);
				G4bool operator==(const CrystalEyeAnticoincidenceHit&) const;

				//  Memory allocation and de-allocation
				inline void* operator new(size_t);
				inline void  operator delete(void*);

				//methods from base class
				virtual void Draw();
				virtual void Print();

 private:
				G4double fEdep;
				G4double fEdepDelay;
				G4int    fAcdID;
				G4int    fOrder;
				G4ThreeVector pos;

 public:

				inline virtual void AddEnergy(G4double de) { fEdep += de;};
				inline virtual void AddEnergyDelay(G4double de) { fEdepDelay += de;};
				inline virtual void SetACDID(G4int i)      { fAcdID = i;};
				inline virtual void SetACDOrder(G4int i)   { fOrder = i;};
				inline virtual void SetPos(G4ThreeVector xyz) { pos=xyz; };

				inline G4double GetEdep()     { return fEdep;};
				inline G4double GetEdepDelay(){ return fEdepDelay;};
				inline G4int    GetACDID()    { return fAcdID;};
				inline G4int    GetACDOrder() { return fOrder;};
				inline G4ThreeVector GetPos() { return pos;};

};
//---------------------------------------------

//using CrystalEyeAnticoincidenceHitsCollection = G4THitsCollection<CrystalEyeAnticoincidenceHit>;
typedef G4THitsCollection<CrystalEyeAnticoincidenceHit> CrystalEyeAnticoincidenceHitsCollection;
extern G4ThreadLocal G4Allocator<CrystalEyeAnticoincidenceHit>* CrystalEyeAnticoincidenceHitAllocator;

//---------------------------------------------
inline void* CrystalEyeAnticoincidenceHit::operator new(size_t)
{
		if (!CrystalEyeAnticoincidenceHitAllocator) {
				CrystalEyeAnticoincidenceHitAllocator = new G4Allocator<CrystalEyeAnticoincidenceHit>;
		}
		void *hit;
		hit = (void *) CrystalEyeAnticoincidenceHitAllocator->MallocSingle();
		return hit;
}

inline void CrystalEyeAnticoincidenceHit::operator delete(void* hit)
{
		if (!CrystalEyeAnticoincidenceHitAllocator) {
				CrystalEyeAnticoincidenceHitAllocator = new G4Allocator<CrystalEyeAnticoincidenceHit>;
		}
		CrystalEyeAnticoincidenceHitAllocator->FreeSingle((CrystalEyeAnticoincidenceHit*) hit);
}

//---------------------------------------------

#endif
