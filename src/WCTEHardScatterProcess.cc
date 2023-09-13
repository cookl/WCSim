// Author : Miao Yu
// Mail : miaoyu@whu.edu.cn
// Description : Rayleigh scattering of anisotropic molecules
// Date : March, 2022


#include "WCTEHardScatterProcess.hh"

#include "G4ios.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "G4OpProcessSubType.hh"
#include "Randomize.hh"
#include "G4String.hh"
#include "G4ProcessType.hh"
#include "G4Event.hh"
#include "G4RunManager.hh"
// #include "MatrixCalc.hh"



// #include "WCTEHardScatterProcess.hh"

// class WCTEHardScatterProcess;
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

WCTEHardScatterProcess::WCTEHardScatterProcess(const G4String& processName, G4ProcessType type)
   : G4VDiscreteProcess(processName, type)
{
  SetProcessSubType(fUserDefined);

  thePhysicsTable = NULL;

  if (verboseLevel > 0) {
    G4cout << GetProcessName() << " is created " << G4endl;
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

WCTEHardScatterProcess::~WCTEHardScatterProcess()
{
  if (thePhysicsTable) {
    thePhysicsTable->clearAndDestroy();
    delete thePhysicsTable;
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


G4VParticleChange*
WCTEHardScatterProcess::PostStepDoIt(const G4Track& aTrack, const G4Step& aStep)
{
    const G4Event* g4event = G4RunManager::GetRunManager()->GetCurrentEvent();
    int eventID = g4event->GetEventID(); 
    std::cout << "***************Hard scatter has occured " << eventID << std::endl;

    const G4DynamicParticle* aParticle = aTrack.GetDynamicParticle();
    double inEnergy = aParticle->GetKineticEnergy();
    // std::cout << "Hard Scatter " << inEnergy/CLHEP::GeV << std::endl;

    aParticleChange.Initialize(aTrack);

    // if ( verboseLevel > 0) {
        // G4cout << "WCTE Hard Scatter" << G4endl;
        // G4cout << "Old Momentum Direction: "
        //        << aParticle->GetMomentumDirection() << G4endl;
        // G4cout << "Old Polarization: "
        //        << aParticle->GetPolarization() << G4endl;
    // int ijk =0;
    // std::cin >> ijk;
    // }

   //  const G4Material* aMaterial = aTrack.GetMaterial();
   //  G4MaterialPropertiesTable* aMaterialPropertiesTable = 
   //      aMaterial->GetMaterialPropertiesTable();
   //  if (aMaterialPropertiesTable->ConstPropertyExists("RHOV"))
   //      rhov = aMaterialPropertiesTable->GetConstProperty("RHOV");
   //  else
   //      rhov = 0.0;

    G4double CosThetar, SinThetar, Phir;
    G4double unit_x, unit_y, unit_z;
    G4ThreeVector NewMomentumDirection, OldMomentumDirection, MidMomentumDirection;
   //  G4ThreeVector NewPolarizationDirection, OldPolarizationDirection, MidPolarizationDirection;
   //  G4double Pr = 1;
   //  G4double Prob = 0;
    
   //  OldMomentumDirection = aParticle->GetMomentumDirection();
   //  OldPolarizationDirection = aParticle->GetPolarization();
   //  OldMomentumDirection = OldMomentumDirection.unit();
   //  OldPolarizationDirection = OldPolarizationDirection.unit();
   //  ///G4cout << "OldPolarizationDirection : " << OldPolarizationDirection.x() << " " << OldPolarizationDirection.y() << " " << OldPolarizationDirection.z() << G4endl;
   //  ///G4cout << "OldMomentumDirection : " << OldMomentumDirection.x() << " " << OldMomentumDirection.y() << " " << OldMomentumDirection.z() << G4endl;

   //  do {
   //      // STEP 1. sample out momentum *********************8
     CosThetar = G4UniformRand() * 2 - 1;
     SinThetar = std::sqrt(1 - CosThetar*CosThetar);

     Phir = twopi * G4UniformRand();
     G4double SinPhir = std::sin(Phir);
     G4double CosPhir = std::cos(Phir);

     unit_x = SinThetar * CosPhir;
     unit_y = SinThetar * SinPhir;
     unit_z = CosThetar;
     NewMomentumDirection.set (unit_x,unit_y,unit_z);

     //uniformly distibuted between 5 MeV and inEnergy
     G4double energy = 0.005*CLHEP::GeV + (inEnergy-0.005*CLHEP::GeV)*(G4UniformRand());

   //      // STEP 2. in polarization vector rotated by the polarization tensor *********************8
   //      G4double oldPolVec[3] = {OldPolarizationDirection.x(), OldPolarizationDirection.y(), OldPolarizationDirection.z()};
   //      G4double* PolTensor = MatrixCalc::rotatePolTensor(rhov);
   //      G4double* midPolVec = MatrixCalc::rotatePolVector(PolTensor, oldPolVec);
   //      G4double* normMidPolVec = MatrixCalc::norm(midPolVec);
   //      MidPolarizationDirection = G4ThreeVector(normMidPolVec[0], normMidPolVec[1], normMidPolVec[2]);

        
   //      // STEP 3.decide the out polarization vector by out momentum and mid polarization vector
   //      G4double cosAng = NewMomentumDirection.x() * MidPolarizationDirection.x() 
   //                      + NewMomentumDirection.y() * MidPolarizationDirection.y()
   //                      + NewMomentumDirection.z() * MidPolarizationDirection.z();
        
   //      if (std::abs(cosAng -1 ) < 1e-5)
   //      {
   //          //G4ThreeVector l1 = perpendicular_vector(MidPolarizationDirection);
            
   //          G4double* tmp_l1 = MatrixCalc::perpendicular_vector(normMidPolVec);
   //          G4double* norm_l1 = MatrixCalc::norm(tmp_l1);
   //          G4ThreeVector l1 = G4ThreeVector(norm_l1[0], norm_l1[1], norm_l1[2]);

   //          G4ThreeVector l2 = MidPolarizationDirection.cross(l1);
   //          G4double beta = G4UniformRand() * twopi;
   //          G4ThreeVector tmp = std::cos(beta) * l1 + std::sin(beta) * l2;
   //          NewPolarizationDirection = G4ThreeVector(tmp.x(), tmp.y(), tmp.z());
   //          delete[] tmp_l1;
   //          delete[] norm_l1;
   //      } else {
   //          G4ThreeVector tmp = MidPolarizationDirection - MidPolarizationDirection.mag() * cosAng * NewMomentumDirection;
   //          tmp = tmp.unit();
   //          NewPolarizationDirection = G4ThreeVector(tmp.x(), tmp.y(), tmp.z());
   //      }

   //      // STEP 4. calculate probability to decide if accept this sampling
   //      Prob = NewPolarizationDirection.x() * MidPolarizationDirection.x()
   //           + NewPolarizationDirection.y() * MidPolarizationDirection.y()
   //           + NewPolarizationDirection.z() * MidPolarizationDirection.z();
        
   //      Pr = G4UniformRand();
        
   //      delete[] PolTensor;
   //      delete[] midPolVec;
   //      delete[] normMidPolVec;
    
   //  } while(Pr > Prob);

   // aParticleChange.ProposePolarization(NewPolarizationDirection);
  // G4ThreeVector* oldMomDirection = aParticleChange.GetMomentumDirection();
   // std::cout << "Before change GetMomentumDirection " << aParticleChange.GetMomentumDirection()->x() << " " << aParticleChange.GetMomentumDirection()->y()  << " " << aParticleChange.GetMomentumDirection()->z()<< std::endl;
   aParticleChange.ProposeMomentumDirection(NewMomentumDirection);
  // G4ThreeVector* newMomDirection = aParticleChange.GetMomentumDirection();
   // std::cout << "After change GetMomentumDirection " << aParticleChange.GetMomentumDirection()->x() << " " << aParticleChange.GetMomentumDirection()->y()  << " " << aParticleChange.GetMomentumDirection()->z()<< std::endl;


   // std::cout << "Before change GetEnergy () " << aParticleChange.GetEnergy() << std::endl;
   
   //NB this is the kinetic energy
   aParticleChange.ProposeEnergy(energy);
   

   // std::cout << "After change GetEnergy " << aParticleChange.GetEnergy() << std::endl;
   // if (verboseLevel > 0) {
   //   G4cout << "New Polarization: "
   //        << NewPolarizationDirection << G4endl;
   //   G4cout << "Polarization Change: "
   //        << *(aParticleChange.GetPolarization()) << G4endl;
   //   G4cout << "New Momentum Direction: "
   //        << NewMomentumDirection << G4endl;
   //   G4cout << "Momentum Change: "
   //        << *(aParticleChange.GetMomentumDirection()) << G4endl;
   // }

    // int ijk =0;
    // std::cin >> ijk;

   // std::cout << "After Hard scatter energy " << inEnergy <<" -> " <<  energy << std::endl;


   return G4VDiscreteProcess::PostStepDoIt(aTrack, aStep);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


void WCTEHardScatterProcess::BuildPhysicsTable(const G4ParticleDefinition&)
{ 

  if (thePhysicsTable) {
    thePhysicsTable->clearAndDestroy();
    delete thePhysicsTable;
    thePhysicsTable = NULL;
  }

  const G4MaterialTable* theMaterialTable = G4Material::GetMaterialTable();
  const G4int numOfMaterials = G4Material::GetNumberOfMaterials();

  thePhysicsTable = new G4PhysicsTable(numOfMaterials);

  // for (G4int iMaterial = 0; iMaterial < numOfMaterials; ++iMaterial)
  // {
  //   G4Material* material = (*theMaterialTable)[iMaterial];
  //   G4MaterialPropertiesTable* materialProperties =
  //                                    material->GetMaterialPropertiesTable();
  //   G4PhysicsOrderedFreeVector* rayleigh = NULL;
  //   if (materialProperties) {
  //     rayleigh = materialProperties->GetProperty(kRAYLEIGH);
  //     if (rayleigh == NULL) rayleigh = CalculateMeanFreePaths(material);
  //   }
  //   thePhysicsTable->insertAt(iMaterial, rayleigh);
  // }


  //laurence simple version 
  for (G4int iMaterial = 0; iMaterial < numOfMaterials; ++iMaterial)
  {
    G4Material* material = (*theMaterialTable)[iMaterial];
    // G4MaterialPropertiesTable* materialProperties =
                                     // material->GetMaterialPropertiesTable();
    // G4PhysicsOrderedFreeVector* rayleigh = NULL;
    // if (materialProperties) {
    //   rayleigh = materialProperties->GetProperty(kRAYLEIGH);
    //   if (rayleigh == NULL) rayleigh = CalculateMeanFreePaths(material);
    // }
    G4PhysicsOrderedFreeVector* rayleigh = CalculateMeanFreePaths(material);                                
    thePhysicsTable->insertAt(iMaterial, rayleigh);
  }



}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double WCTEHardScatterProcess::GetMeanFreePath(const G4Track& aTrack,
                                       G4double ,
                                       G4ForceCondition*)
{
  const G4DynamicParticle* particle = aTrack.GetDynamicParticle();
  const G4double particleMomentum = particle->GetTotalMomentum();
  const G4Material* material = aTrack.GetMaterial();

  G4PhysicsOrderedFreeVector* meanFreePathVector =
                              static_cast<G4PhysicsOrderedFreeVector*>
                              ((*thePhysicsTable)(material->GetIndex()));


  G4double meanFreePath = DBL_MAX;
  if (meanFreePathVector) meanFreePath = meanFreePathVector->Value(particleMomentum);
  
  // std::cout << "Laurence Double Check " << material->GetName() << " has index " << material->GetIndex() << " mfp " << meanFreePath<< std::endl;                            

  return meanFreePath;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
G4PhysicsOrderedFreeVector*
WCTEHardScatterProcess::CalculateMeanFreePaths(const G4Material* material) const
{
  // G4MaterialPropertiesTable* materialProperties =
  //                                      material->GetMaterialPropertiesTable();

  // Retrieve the beta_T or isothermal compressibility value. For backwards
  // compatibility use a constant if the material is "Water". If the material
  // doesn't have an ISOTHERMAL_COMPRESSIBILITY constant then return
  // G4double betat;
  // if (material->GetName() == "Water") {
  //   betat = 7.658e-23*m3/MeV;
  // }
  // else if (materialProperties->ConstPropertyExists("ISOTHERMAL_COMPRESSIBILITY")) {
  //   betat = materialProperties->GetConstProperty(kISOTHERMAL_COMPRESSIBILITY);
  // }
  // else {
  //   return NULL;
  // }

  // // If the material doesn't have a RINDEX property vector then return
  // G4MaterialPropertyVector* rIndex = materialProperties->GetProperty(kRINDEX);
  // if (rIndex == NULL) return NULL;

  // // Retrieve the optional scale factor, (this just scales the scattering length
  // G4double scaleFactor = 1.0;
  // if (materialProperties->ConstPropertyExists("RS_SCALE_FACTOR")) {
  //   scaleFactor = materialProperties->GetConstProperty(kRS_SCALE_FACTOR);
  // }

  // // Retrieve the material temperature. For backwards compatibility use a
  // // constant if the material is "Water"
  // G4double temperature;
  // if (material->GetName() == "Water") {
  //   temperature = 283.15*kelvin; // Temperature of water is 10 degrees celsius
  // }
  // else {
  //   temperature = material->GetTemperature();
  // }

  G4PhysicsOrderedFreeVector* meanFreePaths =
                                             new G4PhysicsOrderedFreeVector();
  // This calculates the meanFreePath via the Einstein-Smoluchowski formula
  // const G4double c1 = scaleFactor * betat * temperature * k_Boltzmann /
  //                     ( 6.0 * pi );
  
  // double energyThreshold = 450.0*CLHEP::MeV;
  // double mfp_water = 40.0; 
  
  double energyThreshold = 0.0*CLHEP::MeV;
  double mfp_water = 400.0; 
  int nEntries =2000;                                           
  for (int index = 0; index < nEntries; ++index)
  {
    const G4double energy = index + 0.1;
  //   const G4double rIndexSquared = (*rIndex)[uRIndex] * (*rIndex)[uRIndex];
  //   const G4double xlambda = h_Planck * c_light / energy;
  //   const G4double c2 = std::pow(twopi/xlambda,4);
  //   const G4double c3 =
  //                  std::pow(((rIndexSquared-1.0)*(rIndexSquared+2.0 )/3.0),2);

  //   const G4double meanFreePath = 1.0 / ( c1 * c2 * c3 );

    // if( verboseLevel > 0) {
    G4double meanFreePath;
    if (material->GetName() == "Water" || material->GetName() == "G4_WATER"){
      if(energy>energyThreshold){
        meanFreePath= mfp_water*CLHEP::cm;
      }else{
        meanFreePath = DBL_MAX;
      }
    }else{
      meanFreePath = DBL_MAX;
    }

    if( true && index%500==0 ) {
      G4cout << material->GetName() << " " <<energy << "MeV\t" << meanFreePath << "mm" << G4endl;
    }

    meanFreePaths->InsertValues(energy, meanFreePath);
  }

  return meanFreePaths;
}
