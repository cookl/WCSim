#ifndef WCTEHardScatterProcess_h
#define WCTEHardScatterProcess_h 1

#include "G4VPhysicsConstructor.hh"
#include "globals.hh"
#include "G4VDiscreteProcess.hh"
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
#include "G4Electron.hh"
#include "G4MuonMinus.hh"

class WCTEHardScatterProcess : public G4VDiscreteProcess
{

public:

    explicit WCTEHardScatterProcess(const G4String& processName= "WCTEHardScatterProcess", G4ProcessType type = fUserDefined);
    virtual ~WCTEHardScatterProcess();


public:

    // virtual G4bool IsApplicable(const G4ParticleDefinition& aParticleType) override;

    // virtual void BuildPhysicsTable(const G4ParticleDefinition& aParticleType) override;

    // virtual G4double GetMeanFreePath(const G4Track& aTrack,
    //                                  G4double, 
    //                                  G4ForceCondition*) override;
    
    // virtual G4VParticleChange* PostStepDoIt(const G4Track& aTrack,
    //                                         const G4Step&  aStep) override;


    virtual G4bool IsApplicable(const G4ParticleDefinition& aParticleType);

    virtual void BuildPhysicsTable(const G4ParticleDefinition& aParticleType);

    virtual G4double GetMeanFreePath(const G4Track& aTrack,
                                     G4double, 
                                     G4ForceCondition*);
    
    virtual G4VParticleChange* PostStepDoIt(const G4Track& aTrack,
                                            const G4Step&  aStep);



    virtual G4PhysicsTable* GetPhysicsTable() const;

    virtual void DumpPhysicsTable() const;

protected:
    G4PhysicsTable* thePhysicsTable;

private:

    // WCTEHardScatterProcess(const WCTEHardScatterProcess &right) = delete;
    // WCTEHardScatterProcess& operator=(const WCTEHardScatterProcess &right) = delete;

    WCTEHardScatterProcess(const WCTEHardScatterProcess &right);
    WCTEHardScatterProcess& operator=(const WCTEHardScatterProcess &right);

    G4PhysicsOrderedFreeVector*
    CalculateMeanFreePaths( const G4Material* material) const;

private:
  G4double rhov;

};



inline
G4bool WCTEHardScatterProcess::IsApplicable(const G4ParticleDefinition& aParticleType)
{
  // return (&aParticleType == G4Electron::Electron());
  return (&aParticleType == G4MuonMinus::MuonMinus());
}

inline
void WCTEHardScatterProcess::DumpPhysicsTable() const
{
  G4int PhysicsTableSize = thePhysicsTable->entries();
  G4PhysicsOrderedFreeVector *v;

  for (G4int i = 0; i < PhysicsTableSize; ++i)
  {
    v = (G4PhysicsOrderedFreeVector*)(*thePhysicsTable)[i];
    v->DumpValues();
  }
}

inline G4PhysicsTable* WCTEHardScatterProcess::GetPhysicsTable() const
{
  return thePhysicsTable;
}


#endif /* G4OpRayleigh_h */
