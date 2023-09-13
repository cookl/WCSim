#include "WCSimPhysicsListFactoryMessenger.hh"
#include "WCSimPhysicsListFactory.hh"

#include "G4UIdirectory.hh"
#include "G4ios.hh"
#include "globals.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithABool.hh"

WCSimPhysicsListFactoryMessenger::WCSimPhysicsListFactoryMessenger(WCSimPhysicsListFactory* WCSimPhysFactory, G4String inValidListsString)
  :thisWCSimPhysicsListFactory(WCSimPhysFactory), ValidListsString(inValidListsString)
{
 
  G4String defaultList="FTFP_BERT";
  
  physListCmd = new G4UIcmdWithAString("/WCSim/physics/list",this);
  G4String cmd_hint = "Available options: " + ValidListsString;
  physListCmd->SetGuidance(cmd_hint);
  physListCmd->SetGuidance("See http://geant4.cern.ch/support/proc_mod_catalog/physics_lists/useCases.shtml");
  physListCmd->SetGuidance("    http://geant4.cern.ch/support/proc_mod_catalog/physics_lists/referencePL.shtml");
  physListCmd->SetGuidance("Note: Physics list is locked-in after initialization");
  
  physListCmd->SetDefaultValue(defaultList);
  physListCmd->SetCandidates(ValidListsString);  // ToDo get list of physics lists from G4PhysicsListFactory

  SetNewValue(physListCmd, defaultList);

  G4String captureModelsString = "Default HP Rad GLG4Sim";
  G4String captureModelGuidance = "Available options: " + captureModelsString;
  nCaptureModelCmd = new G4UIcmdWithAString("/WCSim/physics/nCapture",this);
  nCaptureModelCmd->SetGuidance(captureModelGuidance);
  nCaptureModelCmd->SetDefaultValue("Default");
  nCaptureModelCmd->SetCandidates(captureModelsString);

  //Laurence add messenger for new hard scatter process
  custNuclearModelCmd = new G4UIcmdWithABool("/WCSim/physics/customNuclearScatter",this);
  custNuclearModelCmd->SetGuidance("Choose whether to include the custom nuclear scatter");
  custNuclearModelCmd->SetDefaultValue(false);

}

WCSimPhysicsListFactoryMessenger::~WCSimPhysicsListFactoryMessenger()
{
  delete physListCmd;
  delete nCaptureModelCmd;
  delete custNuclearModelCmd;
  //delete WCSimDir;
}

void WCSimPhysicsListFactoryMessenger::SetNewValue(G4UIcommand* command, G4String newValue)
{
  if (command == physListCmd)
    thisWCSimPhysicsListFactory->SetList(newValue);
  else if (command == nCaptureModelCmd){
    thisWCSimPhysicsListFactory->SetnCaptModel(newValue);
  }else if(command == custNuclearModelCmd){
    bool value = command->ConvertToBool(newValue);
    thisWCSimPhysicsListFactory->SetCustNuclScatter(value);
  }
}
