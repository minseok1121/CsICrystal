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
// $Id: B1PrimaryGeneratorAction.cc 69565 2013-05-08 12:35:31Z gcosmo $
//
/// \file B1PrimaryGeneratorAction.cc
/// \brief Implementation of the B1PrimaryGeneratorAction class

#include "B1PrimaryGeneratorAction.hh"
#include "B1EventAction.hh"

#include "G4LogicalVolumeStore.hh"
#include "G4LogicalVolume.hh"
#include "G4Box.hh"
#include "G4RunManager.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "G4SystemOfUnits.hh"
#include "Randomize.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B1PrimaryGeneratorAction::B1PrimaryGeneratorAction()
: G4VUserPrimaryGeneratorAction(),
  fParticleGun(0), 
  fEnvelopeBox(0)
{
  G4int n_particle = 1;
  fParticleGun  = new G4ParticleGun(n_particle);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B1PrimaryGeneratorAction::~B1PrimaryGeneratorAction()
{
  delete fParticleGun;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void B1PrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{
    // 입자 생성기 초기화
    G4ParticleDefinition* particle = G4ParticleTable::GetParticleTable()->FindParticle("gamma");
    fParticleGun->SetParticleDefinition(particle);
    G4double theta_min = 0.0;
    G4double theta_max = 45.0 * deg;
    
    // 1. Uniform in cos(theta)
    // 세타, 파이 계산
    G4double costheta = std::cos(0.0) - G4UniformRand() * (std::cos(theta_min) - std::cos(theta_max));
    G4double sintheta = std::sqrt(1.0 - costheta * costheta);
    G4double theta = std::acos(costheta);
    G4double phi = 2.0 * CLHEP::pi * G4UniformRand();
    
    // 3. Spherical to Cartesian
    /*
    G4double dirX = sintheta * std::cos(phi);
    G4double dirZ = sintheta * std::sin(phi);
    G4double dirY = costheta;
    */
    G4double dirX = sintheta * std::cos(phi);
    G4double dirY = sintheta * std::sin(phi);
    G4double dirZ = costheta;
    // ---------- 위치 설정 (Z = -1cm, 반지름 1.2cm 원 안) ----------
    //G4double radius = 1.2 * cm;
    G4double radius = 0.1 * cm;
    G4double y_range = 2.95 * cm;
    G4double x = (2.0 * G4UniformRand() - 1.0) * radius;  // x in [-0.1 cm, 0.1 cm]
    G4double y = (2.0 * G4UniformRand() - 1.0) * y_range; // y in [-2.95 cm, 2.95 cm]
    G4ThreeVector position(x, y, -0.25 * cm);
    //G4ThreeVector position(0, -3.1*cm, 0);
    fParticleGun->SetParticlePosition(position);
    /*
    G4double radius = 0.1 * cm;
    double x, z;
    do {
        x = (2.0 * G4UniformRand() - 1.0) * radius;
        z = (2.0 * G4UniformRand() - 1.0) * radius;
    } while (x*x + z*z > radius*radius);
    //G4ThreeVector position(x, -6.1*cm, z);
    G4ThreeVector position(x, -3.1*cm, z);
    //G4ThreeVector position(0, 0, -0.5 * cm);
    fParticleGun->SetParticlePosition(position);
    */
/*
    // ---------- 방향 설정 (Z > 0 방향으로만, 랜덤한 단위 벡터) ----------
    G4double dx, dy, dz;
    do {
        dx = 2.0 * G4UniformRand() - 1.0;
        dy = 2.0 * G4UniformRand() - 1.0;
        dz = G4UniformRand(); // Z > 0만 허용
    } while (dx*dx + dy*dy + dz*dz > 1.0);  // 단위벡터 구체 내에 있는지 확인

    G4ThreeVector direction(dx, dy, dz);
    direction = direction.unit();
    fParticleGun->SetParticleMomentumDirection(direction);
    */

// 이벤트 액션 가져와서 theta/phi 전달
auto eventAction = const_cast<B1EventAction*>(
    static_cast<const B1EventAction*>(
        G4RunManager::GetRunManager()->GetUserEventAction()
    )
);
eventAction->SetThetaPhi(theta, phi, x, y);
// 4. Set momentum direction
fParticleGun->SetParticleMomentumDirection(G4ThreeVector(dirX, dirY, dirZ));

    //fParticleGun->SetParticleMomentumDirection(G4ThreeVector(0, 1, 0));

    // ---------- 에너지 설정 (Co-60: 대표적 감마선 1.17 MeV & 1.33 MeV 중 택일) ----------
    G4double energies[2] = {1.17 * MeV, 1.33 * MeV};
    G4int index = G4RandFlat::shootInt(2);
    //fParticleGun->SetParticleEnergy(energies[index]);
    fParticleGun->SetParticleEnergy(1.33*MeV);

    // 이벤트에 입자 쏘기
    fParticleGun->GeneratePrimaryVertex(anEvent);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

