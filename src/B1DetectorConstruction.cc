#include "B1DetectorConstruction.hh"
#include "B1ScintSD.hh"

#include "G4RunManager.hh"
#include "G4NistManager.hh"
#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4Cons.hh"
#include "G4Orb.hh"
#include "G4Sphere.hh"
#include "G4Trd.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4SystemOfUnits.hh"
#include "G4SDManager.hh"

#include "OrganicScintillatorFactory.hh"
#include "InorganicScintillatorFactory.hh"
#include "G4Material.hh"
#include "G4OpticalSurface.hh"
#include "G4LogicalBorderSurface.hh"
#include "G4OpticalParameters.hh"
#include "G4LogicalSkinSurface.hh"

B1DetectorConstruction::B1DetectorConstruction()
: G4VUserDetectorConstruction(),
  fScoringVolume1(0),
  fScoringVolume2(0),
  fScoringVolume3(0),
  fScoringVolume4(0),
  fScoringVolume5(0)

{ }

B1DetectorConstruction::~B1DetectorConstruction()
{ }

G4VPhysicalVolume* B1DetectorConstruction::Construct()
{  
  G4String name;
  G4double density;
  G4int nel, natoms;

  G4NistManager* nist = G4NistManager::Instance();
  G4Material* det_mat = nist->FindOrBuildMaterial("G4_CESIUM_IODIDE");
//G4OpticalParameters::Instance()->SetScintByParticleType(true);
// EJ301 (det_mat) 정의
G4double photonEnergy[48] = {
  2.06855400 * eV,
  2.10825500 * eV,
  2.14950900 * eV,
  2.19241000 * eV,
  2.23705900 * eV,
  2.28356300 * eV,
  2.33204300 * eV,
  2.38262500 * eV,
  2.43545000 * eV,
  2.49067100 * eV,
  2.54845400 * eV,
  2.60898200 * eV,
  2.67245500 * eV,
  2.73909300 * eV,
  2.80914000 * eV,
  2.88286300 * eV,
  2.96056000 * eV,
  3.04256200 * eV,
  3.12923500 * eV,
  3.22099100 * eV,
  3.31829000 * eV,
  3.42165200 * eV,
  3.51283600 * eV,
  3.57961100 * eV,
  3.63890100 * eV,
  3.67953200 * eV,
  3.71409000 * eV,
  3.74930300 * eV,
  3.77795800 * eV,
  3.81439900 * eV,
  3.84032900 * eV,
  3.86661300 * eV,
  3.92027700 * eV,
  4.02809300 * eV,
  4.14200800 * eV,
  4.20364900 * eV,
  4.23516200 * eV,
  4.25796300 * eV,
  4.29029900 * eV,
  4.32313100 * eV,
  4.37091400 * eV,
  4.41976500 * eV,
  4.47984800 * eV,
  4.55728800 * eV,
  4.62119500 * eV,
  4.70364300 * eV,
  4.85971100 * eV,
  5.08465800 * eV
};

G4double rIndex[48] = {
  1.7371, 1.7382, 1.7392, 1.7397, 1.7412, 1.7434, 1.7457, 1.7479, 1.7502, 1.7526,
  1.7554, 1.7592, 1.7642, 1.7689, 1.7735, 1.7790, 1.7846, 1.7903, 1.7963, 1.8030,
  1.8111, 1.8201, 1.8301, 1.8377, 1.8442, 1.8483, 1.8518, 1.8552, 1.8580, 1.8615,
  1.8639, 1.8663, 1.8735, 1.8884, 1.9038, 1.9154, 1.9212, 1.9254, 1.9312, 1.9370,
  1.9453, 1.9536, 1.9693, 1.9899, 2.0065, 2.0272, 2.0878, 2.1917
};
G4double absLength[48] = {
  1e6 * m, 1e6 * m, 1e6 * m, 1e6 * m, 1e6 * m, 1e6 * m, 1e6 * m, 1e6 * m,
  1e6 * m, 1e6 * m, 1e6 * m, 1e6 * m, 1e6 * m, 1e6 * m, 1e6 * m, 1e6 * m,
  1e6 * m, 1e6 * m, 1e6 * m, 1e6 * m, 1e6 * m, 1e6 * m, 1e6 * m, 1e6 * m,
  1e6 * m, 1e6 * m, 1e6 * m, 1e6 * m, 1e6 * m, 1e6 * m, 1e6 * m, 1e6 * m,
  1e6 * m, 1e6 * m, 1e6 * m, 1e6 * m, 1e6 * m, 1e6 * m, 1e6 * m, 1e6 * m,
  1e6 * m, 1e6 * m, 1e6 * m, 1e6 * m, 1e6 * m, 1e6 * m, 1e6 * m, 1e6 * m
};


G4double scintComp[48] = {
  0.00607500,
  0.01955900,
  0.03124500,
  0.04203300,
  0.05057300,
  0.06046100,
  0.06720300,
  0.07349600,
  0.07619300,
  0.07709200,
  0.07529400,
  0.07349600,
  0.07124900,
  0.06855200,
  0.06540600,
  0.06225900,
  0.06091100,
  0.06225900,
  0.07079900,
  0.08518200,
  0.10630800,
  0.15125500,
  0.22069900,
  0.29913200,
  0.38048700,
  0.46408900,
  0.54184800,
  0.60657200,
  0.67669000,
  0.74231300,
  0.80928500,
  0.87625700,
  0.95446500,
  0.98997400,
  0.95266700,
  0.86816600,
  0.78501400,
  0.69197200,
  0.61106700,
  0.53555500,
  0.46453800,
  0.38453200,
  0.30452600,
  0.21732800,
  0.14855800,
  0.08248600,
  0.03079600,
  0.00517600
};
// Material Properties Table 생성
G4MaterialPropertiesTable* mpt = new G4MaterialPropertiesTable();
mpt->AddProperty("RINDEX", photonEnergy, rIndex, 48);
mpt->AddProperty("ABSLENGTH", photonEnergy, absLength, 48);
mpt->AddProperty("SCINTILLATIONCOMPONENT1", photonEnergy, scintComp, 48);

// 상수 특성 추가
mpt->AddConstProperty("SCINTILLATIONYIELD", 2000.0 / MeV);           // 5~6% of NaI(Tl)
mpt->AddConstProperty("SCINTILLATIONTIMECONSTANT1", 10.0 * ns);      // 0.01 µs
mpt->AddConstProperty("RESOLUTIONSCALE", 1.0);
mpt->AddConstProperty("SCINTILLATIONYIELD1", 1.0);                   // 단일 구성 요소 사용 시

// Material에 할당
det_mat->SetMaterialPropertiesTable(mpt);

  G4Material* world_mat = nist->FindOrBuildMaterial("G4_AIR");

  G4bool checkOverlaps = true;

  G4RotationMatrix* Rotation = new G4RotationMatrix();
  Rotation->rotateX(90*deg);
  Rotation->rotateY(0*deg);
  Rotation->rotateZ(0*deg);

  G4Box* solidWorld =    
    new G4Box("World",                       //its name
       3*m, 3*m, 3*m);     //its size
      
  G4LogicalVolume* logicWorld =                         
    new G4LogicalVolume(solidWorld,          //its solid
                        world_mat,           //its material
                        "World");            //its name
                                   
  G4VPhysicalVolume* physWorld = 
    new G4PVPlacement(0,                     //no rotation
                      G4ThreeVector(),       //at (0,0,0)
                      logicWorld,            //its logical volume
                      "World",               //its name
                      0,                     //its mother  volume
                      false,                 //no boolean operation
                      0,                     //copy number
                      checkOverlaps);        //overlaps checking
  /*
  G4double tarCenter =   12.0*0.5*cm;
  G4Box* solidTarget =    
    new G4Box("Target",                    //its name
        0.5*cm, 0.5*cm, tarCenter); //its size
        */
    G4double tarCenter =   6.0*0.5*cm;
    G4double dx = 0.15*cm;
    G4double dy = tarCenter;  // height → z방향 회전돼서 이제 y축이 됨
    G4double dz = 0.15*cm;

        G4Box* solidTarget =    
          new G4Box("Target",                    //its name
              dx, dz, tarCenter); //its size
  G4LogicalVolume* logicTarget =                         
    new G4LogicalVolume(solidTarget,            //its solid
                        det_mat,             //its material
                        "Target");         //its name
  G4VPhysicalVolume* physTarget =
  new G4PVPlacement(Rotation,                       //no rotation
                    G4ThreeVector(0,0,0),         //at (0,0,0)
                    logicTarget,                //its logical volume
                    "Target",              //its name
                    logicWorld,              //its mother  volume
                    false,                   //no boolean operation
                    0,                       //copy number
                    checkOverlaps);          //overlaps checking
    // ---- Teflon 덮개 Optical Surface 설정 ----
    G4OpticalSurface* teflonSurface = new G4OpticalSurface("TeflonSurface");
    teflonSurface->SetType(dielectric_dielectric);
    teflonSurface->SetFinish(groundfrontpainted);
    teflonSurface->SetModel(unified);
  
    G4MaterialPropertiesTable* teflonMPT = new G4MaterialPropertiesTable();
    const G4int num = 2;
    G4double ephoton[num] = {1.0 * eV, 6.0 * eV};
    G4double reflectivity[num] = {0.93, 0.93};
    //G4double reflectivity[num] = {1, 1};
    teflonMPT->AddProperty("REFLECTIVITY", ephoton, reflectivity, num);
    teflonSurface->SetMaterialPropertiesTable(teflonMPT);

    G4Material* teflon = nist->FindOrBuildMaterial("G4_TEFLON");
    G4double teflonThickness = 0.1;
  
    G4Box* teflonLeft   = new G4Box("TeflonLeft",    teflonThickness/2, dy, dz); // -x
    G4Box* teflonRight  = new G4Box("TeflonRight",   teflonThickness/2, dy, dz); // +x
    G4Box* teflonBottom = new G4Box("TeflonBottom",  dx, teflonThickness/2, dz); // -y
    G4Box* teflonBack   = new G4Box("TeflonBack",    dx, dy, teflonThickness/2); // -z
    G4Box* teflonFront  = new G4Box("TeflonFront",   dx, dy, teflonThickness/2); // +z
    G4LogicalVolume* logicLeft   = new G4LogicalVolume(teflonLeft,   teflon, "TeflonLeft");
    G4LogicalVolume* logicRight  = new G4LogicalVolume(teflonRight,  teflon, "TeflonRight");
    G4LogicalVolume* logicBottom = new G4LogicalVolume(teflonBottom, teflon, "TeflonBottom");
    G4LogicalVolume* logicBack   = new G4LogicalVolume(teflonBack,   teflon, "TeflonBack");
    G4LogicalVolume* logicFront  = new G4LogicalVolume(teflonFront,  teflon, "TeflonFront");
    G4ThreeVector posLeft   = G4ThreeVector(-dx - teflonThickness/2, 0, 0);
    G4ThreeVector posRight  = G4ThreeVector( dx + teflonThickness/2, 0, 0);
    G4ThreeVector posBottom = G4ThreeVector(0, -dy - teflonThickness/2, 0);
    G4ThreeVector posBack   = G4ThreeVector(0, 0, -dz - teflonThickness/2);
    G4ThreeVector posFront  = G4ThreeVector(0, 0,  dz + teflonThickness/2);
    G4VPhysicalVolume* physTefLeft = new G4PVPlacement(0, posLeft,   logicLeft,   "TeflonLeft",   logicWorld, false, 0, checkOverlaps);
    G4VPhysicalVolume* physTefRight =new G4PVPlacement(0, posRight,  logicRight,  "TeflonRight",  logicWorld, false, 0, checkOverlaps);
    G4VPhysicalVolume* physTefBottom = new G4PVPlacement(0, posBottom, logicBottom, "TeflonBottom", logicWorld, false, 0, checkOverlaps);
    G4VPhysicalVolume* physTefBack = new G4PVPlacement(0, posBack,   logicBack,   "TeflonBack",   logicWorld, false, 0, checkOverlaps);
    G4VPhysicalVolume* physTefFront = new G4PVPlacement(0, posFront,  logicFront,  "TeflonFront",  logicWorld, false, 0, checkOverlaps);
    new G4LogicalBorderSurface("SurfaceLeft",   physTarget, physTefLeft,   teflonSurface);
    new G4LogicalBorderSurface("SurfaceRight",  physTarget, physTefRight,  teflonSurface);
    new G4LogicalBorderSurface("SurfaceBottom", physTarget, physTefBottom,  teflonSurface);
    new G4LogicalBorderSurface("SurfaceBack",   physTarget, physTefBack,  teflonSurface);
    new G4LogicalBorderSurface("SurfaceFront",  physTarget, physTefFront,  teflonSurface);

    G4Material* aluminum = nist->FindOrBuildMaterial("G4_Al");
    G4double alThickness = 50 * micrometer;
    G4Box* alLeft   = new G4Box("AlLeft",   alThickness/2, dy, dz);
    G4Box* alRight  = new G4Box("AlRight",  alThickness/2, dy, dz);
    G4Box* alBottom = new G4Box("AlBottom", dx, alThickness/2, dz);
    G4Box* alBack   = new G4Box("AlBack",   dx, dy, alThickness/2);
    G4Box* alFront  = new G4Box("AlFront",  dx, dy, alThickness/2);
    G4LogicalVolume* logicAlLeft   = new G4LogicalVolume(alLeft,   aluminum, "AlLeft");
    G4LogicalVolume* logicAlRight  = new G4LogicalVolume(alRight,  aluminum, "AlRight");
    G4LogicalVolume* logicAlBottom = new G4LogicalVolume(alBottom, aluminum, "AlBottom");
    G4LogicalVolume* logicAlBack   = new G4LogicalVolume(alBack,   aluminum, "AlBack");
    G4LogicalVolume* logicAlFront  = new G4LogicalVolume(alFront,  aluminum, "AlFront");
    G4ThreeVector posAlLeft   = G4ThreeVector(-dx - teflonThickness - alThickness/2, 0, 0);
    G4ThreeVector posAlRight  = G4ThreeVector( dx + teflonThickness + alThickness/2, 0, 0);
    G4ThreeVector posAlBottom = G4ThreeVector(0, -dy - teflonThickness - alThickness/2, 0);
    G4ThreeVector posAlBack   = G4ThreeVector(0, 0, -dz - teflonThickness - alThickness/2);
    G4ThreeVector posAlFront  = G4ThreeVector(0, 0,  dz + teflonThickness + alThickness/2);
    new G4PVPlacement(0, posAlLeft,   logicAlLeft,   "AlLeft",   logicWorld, false, 0, checkOverlaps);
    new G4PVPlacement(0, posAlRight,  logicAlRight,  "AlRight",  logicWorld, false, 0, checkOverlaps);
    new G4PVPlacement(0, posAlBottom, logicAlBottom, "AlBottom", logicWorld, false, 0, checkOverlaps);
    new G4PVPlacement(0, posAlBack,   logicAlBack,   "AlBack",   logicWorld, false, 0, checkOverlaps);
    new G4PVPlacement(0, posAlFront,  logicAlFront,  "AlFront",  logicWorld, false, 0, checkOverlaps);                    

      // ---- SiPM 크기 설정 ----
  G4double sipmSizeXY = 0.6 * cm;
  G4double sipmThickness = 0.1 * cm;
  
  // ---- SiPM Material ----
  G4Material* silicon = nist->FindOrBuildMaterial("G4_Si");

  // ---- SiPM Solid / Logic ----
  G4Box* solidSiPM = new G4Box("SiPM", sipmSizeXY/2, sipmThickness/2, sipmSizeXY/2);
  G4LogicalVolume* logicSiPM = new G4LogicalVolume(solidSiPM, silicon, "SiPM");

// ---- PDE 적용을 위한 MaterialPropertiesTable ----
//G4MaterialPropertiesTable* sipmMPT = new G4MaterialPropertiesTable();
//sipmMPT->AddProperty("EFFICIENCY", photonEnergy, sipmEfficiency, 38);

  // ---- SiPM 위치 (윗면 중앙) ----
  // 원래 타겟은 z 방향으로 길고, x축 기준으로 90도 회전되었으므로
  // y+ 방향 (회전 전 기준 z+) 으로 sipmThickness/2만큼 올려서 붙임
  G4ThreeVector posSiPM = G4ThreeVector(0, dy + sipmThickness/2, 0);
  G4VPhysicalVolume* physSiPM = new G4PVPlacement(0, posSiPM, logicSiPM, "SiPM", logicWorld, false, 0, checkOverlaps);

  fScoringVolume1 = logicWorld;
  fScoringVolume2 = logicTarget;
  fScoringVolume3 = logicSiPM;

  return physWorld;
}

