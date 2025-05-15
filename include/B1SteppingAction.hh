//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
// $Id: B1SteppingAction.hh 74483 2013-10-09 13:37:06Z gcosmo $
//
/// \file B1SteppingAction.hh
/// \brief Definition of the B1SteppingAction class

#ifndef B1SteppingAction_h
#define B1SteppingAction_h 1

#include "G4UserSteppingAction.hh"
#include "globals.hh"
#include "G4SystemOfUnits.hh"
#include <fstream>

class B1EventAction;

class G4LogicalVolume;

/// Stepping action class
/// 

class B1SteppingAction : public G4UserSteppingAction
{
  public:
    B1SteppingAction(B1EventAction* eventAction);
    virtual ~B1SteppingAction();

    // method from the base class
    virtual void UserSteppingAction(const G4Step*);

  private:
    B1EventAction*  fEventAction;
    G4LogicalVolume* fScoringVolume1;
    G4LogicalVolume* fScoringVolume2;
    G4LogicalVolume* fScoringVolume3;
    G4LogicalVolume* fScoringVolume4;
    G4LogicalVolume* fScoringVolume5;
const G4double sipmEfficiency[77] = {
  0.03106050,
  0.03701624,
  0.04395153,
  0.05127865,
  0.05872332,
  0.06640308,
  0.07482731,
  0.08336909,
  0.09213421,
  0.10140086,
  0.11098838,
  0.12042264,
  0.12976113,
  0.13940491,
  0.14856381,
  0.15783047,
  0.16742038,
  0.17674090,
  0.18634621,
  0.19534604,
  0.20533104,
  0.21466953,
  0.22404393,
  0.23311663,
  0.24225398,
  0.25159247,
  0.26061489,
  0.26975224,
  0.28057052,
  0.28940616,
  0.29975752,
  0.31092778,
  0.32048177,
  0.33075411,
  0.34095461,
  0.35035057,
  0.35963159,
  0.36952423,
  0.37900638,
  0.38703319,
  0.39318484,
  0.39835692,
  0.39553579,
  0.38651206,
  0.37771336,
  0.36782072,
  0.35897071,
  0.34983336,
  0.33994893,
  0.32972149,
  0.31698710,
  0.30458240,
  0.29281750,
  0.28011078,
  0.26768340,
  0.25583071,
  0.24397801,
  0.23190981,
  0.22059588,
  0.21044566,
  0.19829126,
  0.18536105,
  0.17544788,
  0.16503188,
  0.15282001,
  0.14381196,
  0.13381260,
  0.12381323,
  0.11284410,
  0.10303869,
  0.09305369,
  0.08321236,
  0.07362245,
  0.06414030,
  0.05088683,
  0.03434693,
  0.02090105
};
const G4double photonEnergy[77] = {
  1.39130000,
  1.41509700,
  1.43972200,
  1.46522000,
  1.49163700,
  1.51902400,
  1.54743600,
  1.57693000,
  1.60615200,
  1.63500800,
  1.66339800,
  1.69121800,
  1.71998500,
  1.74806700,
  1.77534700,
  1.80349300,
  1.83254600,
  1.86064500,
  1.88766100,
  1.91345800,
  1.93790600,
  1.95876000,
  1.97577000,
  1.99526300,
  2.01961600,
  2.04687100,
  2.07250800,
  2.09637800,
  2.12080500,
  2.14580800,
  2.17400000,
  2.20561300,
  2.23541000,
  2.26602400,
  2.29748700,
  2.32685800,
  2.36004600,
  2.40050900,
  2.44566600,
  2.50628200,
  2.58457800,
  2.66792300,
  2.75682200,
  2.84738900,
  2.92518300,
  2.98760100,
  3.04253600,
  3.09425900,
  3.13692100,
  3.17522800,
  3.21286000,
  3.24890800,
  3.28408000,
  3.32002100,
  3.35181300,
  3.38800000,
  3.42497600,
  3.45882000,
  3.49200000,
  3.51492700,
  3.53953400,
  3.57428800,
  3.58838100,
  3.59310400,
  3.61690400,
  3.62410600,
  3.63859500,
  3.65320100,
  3.67533200,
  3.69023500,
  3.70525900,
  3.72040600,
  3.73567800,
  3.75107500,
  3.77441100,
  3.80598100,
  3.85200700
};
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
