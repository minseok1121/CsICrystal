//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under the terms  and *
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
// $Id: B1SteppingAction.cc 74483 2013-10-09 13:37:06Z gcosmo $
//
/// \file B1SteppingAction.cc
/// \brief Implementation of the B1SteppingAction class

#include "B1SteppingAction.hh"
#include "B1EventAction.hh"
#include "B1DetectorConstruction.hh"
#include "B1Analysis.hh"

#include "G4Step.hh"
#include "G4Event.hh"
#include "G4RunManager.hh"
#include "G4LogicalVolume.hh"
#include "G4SystemOfUnits.hh"
#include "G4INCLGlobals.hh"
#include "G4UImanager.hh"
#include <fstream>
#include <iomanip>
#include <cstdlib> // getenv 사용을 위한 헤더 추가
#include "G4Threading.hh"  // 스레드 ID 가져오기 위한 헤더 추가
#include "G4EventManager.hh"
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B1SteppingAction::B1SteppingAction(B1EventAction* eventAction)
: G4UserSteppingAction(),
  fEventAction(eventAction),
  fScoringVolume1(0),
  fScoringVolume2(0),
  fScoringVolume3(0),
  fScoringVolume4(0),
  fScoringVolume5(0)
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B1SteppingAction::~B1SteppingAction()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void B1SteppingAction::UserSteppingAction(const G4Step* step)
{
  auto& outFile = fEventAction->GetOutputFile();
  auto analysisManager = G4RootAnalysisManager::Instance();
  
  // Get the preStepPoint
  G4StepPoint* preStepPoint = step->GetPreStepPoint();
  G4Track* tr = step->GetTrack();
  //G4cout << preStepPoint->GetGlobalTime() / ns << G4endl;
  // Get volume of the current step
  if (!fScoringVolume1) { 
    const B1DetectorConstruction* detectorConstruction
      = static_cast<const B1DetectorConstruction*>
        (G4RunManager::GetRunManager()->GetUserDetectorConstruction());
    fScoringVolume1 = detectorConstruction->GetScoringVolume1();   
    fScoringVolume2 = detectorConstruction->GetScoringVolume2();   
    fScoringVolume3 = detectorConstruction->GetScoringVolume3();   
    fScoringVolume4 = detectorConstruction->GetScoringVolume4();   
    fScoringVolume5 = detectorConstruction->GetScoringVolume5();   
  }
  
  G4LogicalVolume* volume 
    = preStepPoint->GetTouchableHandle()->GetVolume()->GetLogicalVolume();
      
  // Check if we are in scoring volume
  G4String particleName = tr->GetDefinition()->GetParticleName();
  G4UImanager* UI = G4UImanager::GetUIpointer();

  G4LogicalVolume* preVolume = 
      preStepPoint->GetTouchableHandle()->GetVolume()->GetLogicalVolume();
  G4LogicalVolume* postVolume = nullptr;
  if(step->GetPostStepPoint()->GetTouchableHandle() && 
     step->GetPostStepPoint()->GetTouchableHandle()->GetVolume()) {
      postVolume = step->GetPostStepPoint()->GetTouchableHandle()->GetVolume()->GetLogicalVolume();
  }

  // Check scoring volume 1
  if (volume == fScoringVolume1) {
    auto x = preStepPoint->GetPosition()[0] / cm;
    auto y = preStepPoint->GetPosition()[1] / cm;
    auto z = preStepPoint->GetPosition()[2] / cm;
    analysisManager->FillH2(0, z, y);
    analysisManager->FillH2(3, z, x);
    if (tr->GetParticleDefinition()->GetPDGEncoding() == 11) analysisManager->FillH2(1, z, y);
    if (tr->GetParticleDefinition()->GetPDGEncoding() == 22) analysisManager->FillH2(2, z, y);
    if (tr->GetParticleDefinition()->GetPDGEncoding() == 2112) analysisManager->FillH2(4, z, y);
    auto ed = step->GetTotalEnergyDeposit() / MeV;
    analysisManager->FillH1(20, 0, ed);
  }

  // Check scoring volume 2
  if (volume == fScoringVolume2) {
    B1EventAction* eventAction = (B1EventAction*)G4EventManager::GetEventManager()->GetUserEventAction();
    if (tr->GetParticleDefinition()->GetParticleName() == "opticalphoton" &&
    step->GetPostStepPoint()->GetTouchableHandle()->GetVolume()->GetName() == "SiPM") {
analysisManager->FillH1(30,(CLHEP::h_Planck * CLHEP::c_light / step->GetPreStepPoint()->GetKineticEnergy()) / CLHEP::nanometer);
    G4double time = step->GetPostStepPoint()->GetGlobalTime() / ns;
    G4double energy = step->GetPreStepPoint()->GetKineticEnergy() / eV;

if (energy < 4.6 && energy > 3.9) {
  if (G4UniformRand() < 0.35) {
    analysisManager->FillH1(33, time);
    //eventAction->CountDetectedPhoton();
    eventAction->CountDetectedPhoton2();
    //G4cout << "finally?" << G4endl;
    analysisManager->FillH1(32,(CLHEP::h_Planck * CLHEP::c_light / step->GetPreStepPoint()->GetKineticEnergy()) / CLHEP::nanometer);
}
}

if (energy < photonEnergy[0] || energy > photonEnergy[76]) {
  step->GetTrack()->SetTrackStatus(fStopAndKill);
  return;
}
    auto it = std::min_element(std::begin(photonEnergy), std::end(photonEnergy),
    [energy](G4double a, G4double b) {
        return std::abs(a - energy) < std::abs(b - energy);
    });

int idx = std::distance(std::begin(photonEnergy), it);
G4double efficiency = sipmEfficiency[idx];

    // 해당 인덱스의 efficiency
    double eff = sipmEfficiency[idx];

    // 난수 발생 후 efficiency보다 작으면 기록
    if (G4UniformRand() < eff) {
        analysisManager->FillH1(21, time);
        eventAction->CountDetectedPhoton();
        //G4cout << "finally?" << G4endl;
        analysisManager->FillH1(29,(CLHEP::h_Planck * CLHEP::c_light / step->GetPreStepPoint()->GetKineticEnergy()) / CLHEP::nanometer);
    }

    // 감지 여부와 관계없이 Photon 종료
    //G4cout << "killed" << G4endl;
    tr->SetTrackStatus(fStopAndKill);
    return;
}
  G4Track* track = step->GetTrack();

    if (tr->GetParticleDefinition()->GetParticleName() == "opticalphoton" && track->GetTrackStatus() == fStopAndKill) analysisManager->FillH1(35,(CLHEP::h_Planck * CLHEP::c_light / step->GetPreStepPoint()->GetKineticEnergy()) / CLHEP::nanometer);
    // Optical photon 생성 카운트
    if (tr->GetParticleDefinition()->GetParticleName() == "opticalphoton" &&
        tr->GetCurrentStepNumber() == 1) {
      G4double time = step->GetPostStepPoint()->GetGlobalTime() / ns;
      G4double energy = step->GetPostStepPoint()->GetKineticEnergy() / MeV;
      analysisManager->FillH1(22, time);
        eventAction->CountCreatedPhoton();
        analysisManager->FillH1(28,(CLHEP::h_Planck * CLHEP::c_light / step->GetPreStepPoint()->GetKineticEnergy()) / CLHEP::nanometer);
    }
    auto x = preStepPoint->GetPosition()[0] / cm;
    auto y = preStepPoint->GetPosition()[1] / cm;
    auto z = preStepPoint->GetPosition()[2] / cm;
/*
    if (tr->GetParticleDefinition()->GetParticleName() != "opticalphoton"){
  const char* jobId = std::getenv("B1_JOB_ID");
  if (jobId) {
    std::string fileName = "Little_" + std::string(jobId) + "_" +std::to_string(tr->GetParticleDefinition()->GetPDGEncoding())+".txt";
    std::ofstream outFile(fileName, std::ios::app);
    //G4cout << "hmm.." << G4endl;
    outFile << tr->GetTrackID() << "     " << x << "     " << y << "     " << z << std::endl;
  }
}
  */
    auto r = sqrt(x * x + y * y);
    auto theta = preStepPoint->GetMomentumDirection().theta();
    auto e = preStepPoint->GetKineticEnergy() / MeV;
    //G4cout << tr->GetPosition()[2] / cm << "=====" << z <<"     "<< tr->GetKineticEnergy() / MeV << "     " << e << G4endl;
    auto ed = step->GetTotalEnergyDeposit() / MeV;
    if(ed > 0) eventAction->TotalED(ed);
    analysisManager->FillH1(20, 1, ed);
    //if(ed > 1) G4cout << tr->GetParticleDefinition()->GetPDGEncoding() <<"   "<< z <<"   "<< y <<"   "<< ed << G4endl;
if (tr->GetParticleDefinition()->GetPDGEncoding() == 11) {
      auto secondaries = step->GetSecondaryInCurrentStep();
    if (!secondaries || secondaries->empty()) return;
    auto process = step->GetPostStepPoint()->GetProcessDefinedStep();
    if (!process) return;
    G4String procName = process ? process->GetProcessName() : "unknown_process";

    G4double totalOpticalE = 0.0;
    int opticalCount = 0;

    std::vector<std::pair<G4String, G4double>> nonOpticalSecondaries;

    for (auto* sec : *secondaries) {
        G4String secName = sec->GetDefinition()->GetParticleName();
        G4double secE = sec->GetKineticEnergy();

        if (secName == "opticalphoton") {
            totalOpticalE += secE;
            opticalCount++;
        } else {
            nonOpticalSecondaries.emplace_back(secName, secE);
        }
    }

    // 총 opticalphoton 정보 출력
    if (opticalCount > 0) {
      if(step->GetTrack()->GetParentID()==1) outFile << "[HIGH-e-from-gamma] ";
        outFile << procName << " [e- secondaries] TrackLength = " << tr->GetTrackLength() / mm
               << " mm, Total Optical Energy = " << totalOpticalE / eV
               << " eV, Optical Count = " << opticalCount << std::endl;
    }

    // 나머지 secondaries 출력
    for (const auto& [name, energy] : nonOpticalSecondaries) {
      if(step->GetTrack()->GetParentID()==1) outFile << "[HIGH-e-from-gamma] ";
        outFile << procName << " [e- secondaries] -> secondary = " << name
               << ", E = " << energy / MeV << " MeV" << std::endl;
    }
      analysisManager->FillH2(5, z+1, abs(6+y), ed);
      analysisManager->FillH2(6, abs(x+1), abs(6+y), ed);
      analysisManager->FillH2(7, e, theta);
      analysisManager->FillH2(8, e, r);
    }
    if (tr->GetParticleDefinition()->GetPDGEncoding() == 22) {
      auto process = step->GetPostStepPoint()->GetProcessDefinedStep();
    if (!process) return;
    G4String procName = process ? process->GetProcessName() : "unknown_process";
        G4double trackLength = tr->GetTrackLength(); // 지금까지의 감마선 트랙 길이
        G4StepPoint* postPoint = step->GetPostStepPoint();

        const auto& secondaries = *step->GetSecondaryInCurrentStep();
        for (auto* sec : secondaries) {
            if (sec->GetDefinition()->GetParticleName() == "e-") {
                G4double e_energy = sec->GetKineticEnergy(); // 생성된 전자의 에너지
                // 에너지가 0.1 MeV 이상이면 특별 표시
          //if (e_energy > 0.1 * MeV) {
            if (tr->GetTrackID()==1) {
              outFile << "[HIGH-e-from-gamma] ";
          }
                // 로그 출력 (또는 파일 기록 가능)
                outFile << procName << "  TrackLength = " << trackLength / mm
                       << " mm, e- energy = " << e_energy / MeV << " MeV" << std::endl;
          }
        }
      analysisManager->FillH2(9, z+1, abs(6+y), ed);
      analysisManager->FillH2(10, abs(x+1), abs(6+y), ed);
      analysisManager->FillH2(11, e, theta);
      analysisManager->FillH2(12, e, r);
      analysisManager->FillH1(31, z+0.15, e);
    }
    if (tr->GetParticleDefinition()->GetPDGEncoding() == 2212) {
      analysisManager->FillH2(13, z, abs(y));
      analysisManager->FillH2(14, abs(x), abs(y));
      analysisManager->FillH2(15, e, theta);
      analysisManager->FillH2(16, e, r);
    }
    if (tr->GetParticleDefinition()->GetPDGEncoding() == 2112) {
      analysisManager->FillH2(17, z, abs(y));
      analysisManager->FillH2(18, abs(x), abs(y));
      analysisManager->FillH2(19, e, theta);
      analysisManager->FillH2(20, e, r);
    }
    if(postVolume != fScoringVolume2){
      if (tr->GetParticleDefinition()->GetPDGEncoding() == 11) analysisManager->FillH2(63, z+1, abs(6+y), e);
      if (tr->GetParticleDefinition()->GetPDGEncoding() == -11) analysisManager->FillH2(64, e, theta);
      if (tr->GetParticleDefinition()->GetPDGEncoding() == 22) analysisManager->FillH2(65, z+1, abs(6+y), e);
      if (tr->GetParticleDefinition()->GetPDGEncoding() == 2112) analysisManager->FillH2(66, e, theta);
      if (tr->GetParticleDefinition()->GetPDGEncoding() == 2212) analysisManager->FillH2(67, e, theta);
    }
  }

  if(volume != fScoringVolume2 && postVolume == fScoringVolume2){
G4int threadID = G4Threading::G4GetThreadId();
G4int eventID = G4RunManager::GetRunManager()->GetCurrentEvent()->GetEventID();
G4int trackID = tr->GetTrackID();
auto x = step->GetPostStepPoint()->GetPosition()[0] / cm;
auto y = step->GetPostStepPoint()->GetPosition()[1] / cm;
auto z = step->GetPostStepPoint()->GetPosition()[2] / cm;
auto r = sqrt(x * x + y * y);
auto e = step->GetPostStepPoint()->GetKineticEnergy() / MeV;
auto t = step->GetPostStepPoint()->GetGlobalTime() / ns; 
        if (tr->GetParticleDefinition()->GetPDGEncoding() == 11) analysisManager->FillH2(68, e, t);
        if (tr->GetParticleDefinition()->GetPDGEncoding() == -11) analysisManager->FillH2(69, e, t);
        if (tr->GetParticleDefinition()->GetPDGEncoding() == 22) analysisManager->FillH2(70, e, t);
        if (tr->GetParticleDefinition()->GetPDGEncoding() == 2112) analysisManager->FillH2(71, e, t);
    if (tr->GetParticleDefinition()->GetPDGEncoding() == 2212) analysisManager->FillH2(72, e, t);
        if (tr->GetParticleDefinition()->GetPDGEncoding() == 11) analysisManager->FillH2(78, z, y);
        if (tr->GetParticleDefinition()->GetPDGEncoding() == -11) analysisManager->FillH2(79, z, y);
        if (tr->GetParticleDefinition()->GetPDGEncoding() == 22) analysisManager->FillH2(80, z, y);
        if (tr->GetParticleDefinition()->GetPDGEncoding() == 2112) analysisManager->FillH2(81, z, y);
    if (tr->GetParticleDefinition()->GetPDGEncoding() == 2212) analysisManager->FillH2(82, z, y);
        if (tr->GetParticleDefinition()->GetPDGEncoding() == 11) analysisManager->FillH1(23, t, e);
        if (tr->GetParticleDefinition()->GetPDGEncoding() == -11) analysisManager->FillH1(24, t, e);
        if (tr->GetParticleDefinition()->GetPDGEncoding() == 22) analysisManager->FillH1(25, t, e);
        if (tr->GetParticleDefinition()->GetPDGEncoding() == 2112) analysisManager->FillH1(26, t, e);
    if (tr->GetParticleDefinition()->GetPDGEncoding() == 2212) analysisManager->FillH1(27, t, e);
}
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

