#include "DetectorConstruction.hh"
#include "DetectorMessenger.hh"

#include "G4Material.hh"

#include "G4Tubs.hh"
#include "G4Box.hh"
#include "G4Sphere.hh"
#include "G4Trd.hh"
#include "G4SubtractionSolid.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"

#include "G4GeometryManager.hh"
#include "G4PhysicalVolumeStore.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4SolidStore.hh"
#include "G4AssemblyVolume.hh"

#include "G4VisAttributes.hh"
#include "G4Colour.hh"

#include "DetectionSystemSceptar2.hh"

#include "G4SystemOfUnits.hh" // new version geant4.10 requires units

// Class initialisation
DetectionSystemSceptar2::DetectionSystemSceptar2() {
    errors = 0;
    rotateNull = new G4RotationMatrix;      // Zero rotation matrix
    moveNull = G4ThreeVector(0.0,0.0,0.0);  // Zero translation vector

    // parameters:
    scintWidth = 10*mm;
    scintThick = 1.6*mm;
    strucFrontGap = 0.5*mm;
    mountThick = 12.0*mm;
    paintThick = 100.0*um;
    // paintThick = 0.5*mm;
    
    scintDist = (paintThick*2+scintWidth+strucFrontGap)/tan(2*M_PI/16);

    materialStructureDefined = false;
    buildStrucMount = false;

    SetStructureThickness(scintThick);
    SetStructureMaterial("Mylar");
    DefineCoords();
}

// Class destruction
DetectionSystemSceptar2::~DetectionSystemSceptar2() {
    delete fExpHallLog;
    delete fSipmShroudLog;
    delete fSipmActiveLog;
    delete fSipmConductLog;
    delete fScintillatorUpStrLog;
    delete fScintillatorDnStrLog;
    delete fStructureLog;
    delete fStructureMountLog;
    delete fUpAssembly;
    delete fDnAssembly;
    delete fStructureAssembly;
    delete fStructureMountAssembly;
}

void DetectionSystemSceptar2::Debrief() {
    G4cout << "SCEPTAR II built and placed with: " << G4endl;
    G4cout << " >> buildScint = " << buildScint << G4endl;
    G4cout << " >> buildSipm = " << buildSipm << G4endl;
    G4cout << " >> buildStruc = " << buildStruc << G4endl;
    G4cout << " >> strucThick = " << strucThick << " mm" << G4endl;
    G4cout << " >> fMaterialStructure = " << fMaterialStructure->GetName() << G4endl;
    G4cout << " >> errors = " << errors << G4endl;
}

void DetectionSystemSceptar2::SetConfig(G4int config, G4bool struc, G4bool sipm) {

    /*  config: 0 
                1 downstream + zero degree
                2 upstream
                3 upstream + downstream + zero degree
    */

    //                      0  1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 uS dS
    int allBuildPos[20]  = {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1}; // all positions
    int downBuildPos[20] = {0, 0, 0, 0, 3, 3, 3, 3, 3, 3, 3, 3, 1, 1, 1, 1, 0, 1}; // downstream
    int upBuildPos[20]   = {1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 2, 2, 0, 0, 0, 0, 0, 0}; // upstream
    int testPos[20]      = {1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}; // zero-degree downstream only

    // if (abs(config) == 1) { for (int i = 0; i < 20; i++) fBuildPos[i] = 0; } 
    if (abs(config) == 1) { for (int i = 0; i < 20; i++) fBuildPos[i] = upBuildPos[i]; } 
    if (abs(config) == 2) { for (int i = 0; i < 20; i++) fBuildPos[i] = downBuildPos[i]; } 
    if (abs(config) == 3) { for (int i = 0; i < 20; i++) fBuildPos[i] = allBuildPos[i]; } 
    if (abs(config) == 4) { for (int i = 0; i < 20; i++) fBuildPos[i] = testPos[i]; }
    if (config ==0) { buildScint = false; buildSipm = false; buildStruc = false; }
    if (config > 0) { buildScint = true; buildSipm = sipm; buildStruc = struc; }
    if (config < 0) { buildScint = false; buildSipm = false; buildStruc = true; }

    buildStrucMount = buildStruc;
}

void DetectionSystemSceptar2::BuildDetector() {

    fUpAssembly = new G4AssemblyVolume();
    fDnAssembly = new G4AssemblyVolume();
    fStructureAssembly = new G4AssemblyVolume();
    fHalfStructureAssembly = new G4AssemblyVolume();
    fStructureMountAssembly = new G4AssemblyVolume();

    if (buildScint)         ConstructScintillator();
    if (buildStruc)         ConstructStructure();
    if (buildStrucMount)    ConstructStructureMount();
    if (buildSipm)          ConstructSipm();
}

void DetectionSystemSceptar2::PlaceDetector(G4LogicalVolume* expHallLog) {

    G4int positionNumber = 0;                // Griffin position number
    G4double x, y, z;                        // cartesian coordinates
    G4double theta, phi, alpha, beta, gamma; // coordinates for rotation matrix

    // GRIFFIN POSITIONS: 

    for ( int i = 0; i <= 17; i++ ) {
        positionNumber = i;

        if (fBuildPos[positionNumber] == 0) continue;

        theta  = fCoords[positionNumber][0]*deg;
        phi    = fCoords[positionNumber][1]*deg;
        alpha  = fCoords[positionNumber][2]*deg; // yaw
        beta   = fCoords[positionNumber][3]*deg; // pitch
        gamma  = fCoords[positionNumber][4]*deg; // roll

        rotate = new G4RotationMatrix;
        rotate->rotateX(M_PI/2.0);
        rotate->rotateX(alpha);
        rotate->rotateY(beta);
        rotate->rotateZ(gamma);

        x = -scintWidth/2 - paintThick; y = scintWidth/2 + paintThick; z = scintThick/2 + paintThick + scintDist;
        moveB = G4ThreeVector(TransX(x,y,z,theta,phi), TransY(x,y,z,theta,phi), TransZ(x,z,theta));
        x = scintWidth/2 + paintThick; y = scintWidth/2 + paintThick; z = scintThick/2 + paintThick + scintDist;
        moveG = G4ThreeVector(TransX(x,y,z,theta,phi), TransY(x,y,z,theta,phi), TransZ(x,z,theta));
        x = scintWidth/2 + paintThick; y = -scintWidth/2 - paintThick; z = scintThick/2 + paintThick + scintDist;
        moveR = G4ThreeVector(TransX(x,y,z,theta,phi), TransY(x,y,z,theta,phi), TransZ(x,z,theta));
        x = -scintWidth/2 - paintThick; y = -scintWidth/2 - paintThick; z = scintThick/2 + paintThick + scintDist;
        moveW = G4ThreeVector(TransX(x,y,z,theta,phi), TransY(x,y,z,theta,phi), TransZ(x,z,theta));

        if ( buildScint || buildSipm ) {

            // fDet = {0,1,2,3}. Pure upstream.
            if ( positionNumber < 4 && fBuildPos[positionNumber] == 1 ) {
                if (fBuildPos[positionNumber] == 1) {
                    fUpAssembly->MakeImprint(expHallLog, moveB, rotate, 0);
                    fUpAssembly->MakeImprint(expHallLog, moveG, rotate, 0);
                    fUpAssembly->MakeImprint(expHallLog, moveR, rotate, 0);
                    fUpAssembly->MakeImprint(expHallLog, moveW, rotate, 0);
                }
            }

            // fDet = {4,5,6,7,8,9,10,11}. Mixed upstream/downstream.
            if ( positionNumber > 3 && positionNumber < 12 ) {
                if (fBuildPos[positionNumber] != 2) {
                    fDnAssembly->MakeImprint(expHallLog, moveG, rotate, 0);
                    fDnAssembly->MakeImprint(expHallLog, moveR, rotate, 0);
                }
                if (fBuildPos[positionNumber] != 3) {
                    fUpAssembly->MakeImprint(expHallLog, moveB, rotate, 0);
                    fUpAssembly->MakeImprint(expHallLog, moveW, rotate, 0);
                }
            }

            // fDet = {12,13,14,15}. Pure downstream.
            if ( positionNumber > 11 && positionNumber < 16 ) {
                if (fBuildPos[positionNumber] == 1) {
                    fDnAssembly->MakeImprint(expHallLog, moveB, rotate, 0);
                    fDnAssembly->MakeImprint(expHallLog, moveG, rotate, 0);
                    fDnAssembly->MakeImprint(expHallLog, moveR, rotate, 0);
                    fDnAssembly->MakeImprint(expHallLog, moveW, rotate, 0);
                }
            }

            // fDet = {17}. Zero-degree downstream.
            if ( positionNumber == 17 ) {
                if (fBuildPos[positionNumber] == 1) {
                    fDnAssembly->MakeImprint(expHallLog, moveB, rotate, 0);
                    fDnAssembly->MakeImprint(expHallLog, moveG, rotate, 0);
                    fDnAssembly->MakeImprint(expHallLog, moveR, rotate, 0);
                    fDnAssembly->MakeImprint(expHallLog, moveW, rotate, 0);
                }
            }

        }

        if (buildStruc) {
            if (fBuildPos[positionNumber] == 1) {
                x = 0.0; y = 0.0; z = strucThick/2 + scintDist;
                move = G4ThreeVector(TransX(x,y,z,theta,phi), TransY(x,y,z,theta,phi), TransZ(x,z,theta));
                fStructureAssembly->MakeImprint(expHallLog, move, rotate, 0);
                if (positionNumber == 17) {
                    x = 0.0; y = 0.0; z = strucThick + scintDist + mountThick/2;
                    move = G4ThreeVector(TransX(x,y,z,theta,phi), TransY(x,y,z,theta,phi), TransZ(x,z,theta));
                    fStructureMountAssembly->MakeImprint(expHallLog, move, rotate, 0);
                }
                if (positionNumber == 16) {
                    x = 0.0; y = 0.0; z = strucThick + scintDist + mountThick/2;
                    move = G4ThreeVector(TransX(x,y,z,theta,phi), TransY(x,y,z,theta,phi), TransZ(x,z,theta));
                    fStructureMountAssembly->MakeImprint(expHallLog, move, rotate, 0);
                }
            }
            if (fBuildPos[positionNumber] == 2) {
                x = 0.0; y = 0.0; z = strucThick/2 + scintDist;
                move = G4ThreeVector(TransX(x,y,z,theta,phi), TransY(x,y,z,theta,phi), TransZ(x,z,theta));
                fHalfStructureAssembly->MakeImprint(expHallLog, move, rotate, 0);
            }
            if (fBuildPos[positionNumber] == 3) {
                rotate = new G4RotationMatrix;
                rotate->rotateX(-M_PI/2.0);
                rotate->rotateX(alpha);
                rotate->rotateY(beta);
                rotate->rotateZ(gamma);
                x = 0.0; y = 0.0; z = strucThick/2 + scintDist;
                move = G4ThreeVector(TransX(x,y,z,theta,phi), TransY(x,y,z,theta,phi), TransZ(x,z,theta));
                fHalfStructureAssembly->MakeImprint(expHallLog, move, rotate, 0);
            }
        }
    }
}

G4double DetectionSystemSceptar2::TransX(G4double x, G4double y, G4double z, G4double theta, G4double phi) {
    return (x*cos(theta)+z*sin(theta))*cos(phi)-y*sin(phi);
}

G4double DetectionSystemSceptar2::TransY(G4double x, G4double y, G4double z, G4double theta, G4double phi) {
    return (x*cos(theta)+z*sin(theta))*sin(phi)+y*cos(phi);
}

G4double DetectionSystemSceptar2::TransZ(G4double x, G4double z, G4double theta) {
    return -x*sin(theta)+z*cos(theta);
}

void DetectionSystemSceptar2::ConstructScintillator() {

    G4VisAttributes* scintillatorVisAtt = new G4VisAttributes(G4Colour(0.6,0.6,0.6));
    scintillatorVisAtt->SetVisibility(true);

    fMaterialScint = G4Material::GetMaterial("BC404");
    CheckMaterial(fMaterialScint,"BC404");
    fMaterialScintPaint = G4Material::GetMaterial("Titanium");
    CheckMaterial(fMaterialScintPaint,"Titanium");

    G4Box* scintillator = SquareScintillator();
    fScintillatorUpStrLog = new G4LogicalVolume(scintillator, fMaterialScint, "sceptar2ScintillatorUpStrLog", 0, 0, 0); // upstream
    fScintillatorDnStrLog = new G4LogicalVolume(scintillator, fMaterialScint, "sceptar2ScintillatorDnStrLog", 0, 0, 0); // downstream
    fScintillatorUpStrLog->SetVisAttributes(scintillatorVisAtt);
    fScintillatorDnStrLog->SetVisAttributes(scintillatorVisAtt);
    fUpAssembly->AddPlacedVolume(fScintillatorUpStrLog, moveNull, rotateNull);
    fDnAssembly->AddPlacedVolume(fScintillatorDnStrLog, moveNull, rotateNull);

    // Paint:
    G4VisAttributes* scintPaintVisAtt = new G4VisAttributes(G4Colour(0.8,0.8,0.8));
    scintPaintVisAtt->SetVisibility(true);
    G4VSolid* scintPaint = ScintPaint();
    fScintPaintUpStrLog = new G4LogicalVolume(scintPaint, fMaterialScintPaint, "sceptar2ScintPaintUpStrLog", 0, 0, 0); // upstream
    fScintPaintDnStrLog = new G4LogicalVolume(scintPaint, fMaterialScintPaint, "sceptar2ScintPaintDnStrLog", 0, 0, 0); // downstream
    fScintPaintUpStrLog->SetVisAttributes(scintPaintVisAtt);
    fScintPaintDnStrLog->SetVisAttributes(scintPaintVisAtt);
    fUpAssembly->AddPlacedVolume(fScintPaintUpStrLog, moveNull, rotateNull);
    fDnAssembly->AddPlacedVolume(fScintPaintDnStrLog, moveNull, rotateNull);
}

void DetectionSystemSceptar2::ConstructSipm() {

    G4VisAttributes* sipmShroudVisAtt = new G4VisAttributes(G4Colour(0.6,0.2,0.2));
    sipmShroudVisAtt->SetVisibility(true);
    G4VisAttributes* sipmActiveVisAtt = new G4VisAttributes(G4Colour(0.2,0.2,0.6));
    sipmActiveVisAtt->SetVisibility(true);
    G4VisAttributes* sipmConductVisAtt = new G4VisAttributes(G4Colour(0.6,0.4,0.0));
    sipmConductVisAtt->SetVisibility(true);

    // G4double sipmPlasticThick = 0.21*mm;
    G4double sipmConductThick = 0.01*mm;
    G4double sipmActiveThick  = 0.20*mm;
    G4double sipmActiveDim = 1.00*mm;
    G4double sipmDim1 = 1.5*mm;
    G4double sipmDim2 = 1.8*mm;
    G4double sipmThick = 0.65*mm;

    G4VSolid* sipmPlastic = new G4Box("sipmPlastic",sipmThick/2, sipmDim1/2, sipmDim2/2);
    G4VSolid* sipmActive = new G4Box("sipmActive", sipmActiveThick/2, sipmActiveDim/2, sipmActiveDim/2);
    G4VSolid* sipmConduct = new G4Box("sipmConduct", sipmConductThick/2, sipmDim1/2, sipmDim2/2);
    G4VSolid* sipmShroud = new G4SubtractionSolid("sipmShroud", sipmPlastic, sipmActive, 0, moveNull);

    G4Material* fMaterialSipmShroud = G4Material::GetMaterial("Mylar");
    CheckMaterial(fMaterialSipmShroud,"Mylar");

    G4Material* fMaterialSipmActive = G4Material::GetMaterial("Silicon");
    CheckMaterial(fMaterialSipmActive,"Silicon");

    G4Material* fMaterialSipmConduct = G4Material::GetMaterial("Copper");
    CheckMaterial(fMaterialSipmConduct,"Copper");

    move = G4ThreeVector(scintThick/2+sipmThick/2,0.0,0.0);

    fSipmShroudLog = new G4LogicalVolume(sipmShroud, fMaterialSipmShroud, "sceptar2SipmShroudLog", 0, 0, 0);
    fSipmShroudLog->SetVisAttributes(sipmShroudVisAtt);
    fUpAssembly->AddPlacedVolume(fSipmShroudLog, move, rotateNull);
    fDnAssembly->AddPlacedVolume(fSipmShroudLog, move, rotateNull);

    fSipmActiveLog = new G4LogicalVolume(sipmActive, fMaterialSipmActive, "sceptar2SipmActiveLog", 0, 0, 0);
    fSipmActiveLog->SetVisAttributes(sipmActiveVisAtt);
    fUpAssembly->AddPlacedVolume(fSipmActiveLog, move, rotateNull);
    fDnAssembly->AddPlacedVolume(fSipmActiveLog, move, rotateNull);
   
    move = G4ThreeVector(scintThick/2+sipmThick+sipmConductThick,0.0,0.0);

    fSipmConductLog = new G4LogicalVolume(sipmConduct, fMaterialSipmConduct, "sceptar2SipmConductLog", 0, 0, 0);
    fSipmConductLog->SetVisAttributes(sipmConductVisAtt);
    fUpAssembly->AddPlacedVolume(fSipmConductLog, move, rotateNull);
    fDnAssembly->AddPlacedVolume(fSipmConductLog, move, rotateNull);
}

void DetectionSystemSceptar2::ConstructStructure() {

    if (!materialStructureDefined) {
        G4cout << " SCEPTAR 2 Structure Material Undefined." << G4endl;
        errors++;
        return;
    }

    G4VisAttributes* structureVisAtt = new G4VisAttributes(G4Colour(0.6,0.2,0.2));
    structureVisAtt->SetVisibility(true);

    G4VSolid* squareStruct = SquareStructure();
    G4VSolid* halfStruct = HalfStructure();

    // Structure Unit
    fStructureLog = new G4LogicalVolume(squareStruct, fMaterialStructure, "sceptar2StructureLog", 0, 0, 0);
    fStructureLog->SetVisAttributes(structureVisAtt);
    fStructureAssembly->AddPlacedVolume(fStructureLog, moveNull, rotateNull);

    // Structure Half Unit
    fHalfStructureLog = new G4LogicalVolume(halfStruct, fMaterialStructure, "sceptar2halfStructureLog", 0, 0, 0);
    fHalfStructureLog->SetVisAttributes(structureVisAtt);
    fHalfStructureAssembly->AddPlacedVolume(fHalfStructureLog, moveNull, rotateNull);
}

void DetectionSystemSceptar2::ConstructStructureMount() {

    if (!materialStructureDefined) {
        G4cout << "@@@ SCEPTAR 2 Structure Material Undefined." << G4endl;
        errors++;
        return;
    }

    G4VisAttributes* structureVisAtt = new G4VisAttributes(G4Colour(0.6,0.2,0.2));
    structureVisAtt->SetVisibility(true);

    // G4double mountBoxWidth = scintWidth + strucGap + strucThick*tan(2*M_PI/16);
    G4double mountBoxWidth = scintWidth + strucGap/2 + strucThick*tan(2*M_PI/16);
    G4VSolid* mountBox = new G4Box("mountBox", mountThick/2, mountBoxWidth, mountBoxWidth);
    G4VSolid* mountCut = new G4Box("mountCut", mountThick, mountBoxWidth*1.2, scintWidth);
    
    G4VSolid* mountStruc = new G4SubtractionSolid("mountStruc", mountBox, mountCut, rotateNull, moveNull);
    rotate = new G4RotationMatrix;
    rotate->rotateX(90.*deg);
    G4VSolid* mountStruc2 = new G4SubtractionSolid("mountStruc", mountStruc, mountCut, rotate, moveNull);

    fStructureMountLog = new G4LogicalVolume(mountStruc2, fMaterialStructure, "sceptar2StructureMountLog", 0, 0, 0);
    fStructureMountLog->SetVisAttributes(structureVisAtt);
    fStructureMountAssembly->AddPlacedVolume(fStructureMountLog, moveNull, rotateNull);
}

G4Box* DetectionSystemSceptar2::SquareScintillator() {
    G4double boxZ = scintWidth;
    G4double boxY = scintWidth;
    G4double boxX = scintThick;
    G4Box* squareScintillator = new G4Box("squareScintillator", boxX/2, boxY/2, boxZ/2);
    return squareScintillator;
}

G4VSolid* DetectionSystemSceptar2::ScintPaint() {
    G4Box* scintillator = SquareScintillator();
    G4Box* paintBox = new G4Box("paintBox",scintThick/2+paintThick,scintWidth/2+paintThick,scintWidth/2+paintThick);
    G4VSolid* paintShroud = new G4SubtractionSolid("paintShroud",paintBox,scintillator,0,moveNull);
    return paintShroud;
}

G4VSolid* DetectionSystemSceptar2::SquareStructure() {
    G4VSolid* strucPlane = new G4Box("strucPlane",strucThick/2, scintWidthPlusPaint*4/3, scintWidthPlusPaint*4/3);
    G4VSolid* scintHole = new G4Box("scintHole",strucThick, scintWidthPlusPaint, scintWidthPlusPaint);
    G4VSolid* struc1 = new G4SubtractionSolid("struc1", strucPlane, scintHole, 0, moveNull);
    G4VSolid* strucCut = new G4Box("strucCut", cutWidth, cutWidth, scintWidthPlusPaint*2);

    rotate = new G4RotationMatrix;
    rotate->rotateX(90.*deg);
    rotate->rotateZ(2*M_PI/16);
    move = G4ThreeVector(0.0,0.0,cutDist);
    G4VSolid* struc2 = new G4SubtractionSolid("struc2", struc1, strucCut, rotate, move);

    rotate = new G4RotationMatrix;
    rotate->rotateX(90.*deg);
    rotate->rotateZ(-2*M_PI/16);
    move = G4ThreeVector(0.0,0.0,-cutDist);
    G4VSolid* struc3 = new G4SubtractionSolid("struc3", struc2, strucCut, rotate, move);

    rotate = new G4RotationMatrix;
    rotate->rotateZ(-2*M_PI/16);
    move = G4ThreeVector(0.0,cutDist,0.0);
    G4VSolid* struc4 = new G4SubtractionSolid("struc4", struc3, strucCut, rotate, move);

    rotate = new G4RotationMatrix;
    rotate->rotateZ(2*M_PI/16);
    move = G4ThreeVector(0.0,-cutDist,0.0);
    G4VSolid* struc5 = new G4SubtractionSolid("struc5", struc4, strucCut, rotate, move);

    return struc5;
}

G4VSolid* DetectionSystemSceptar2::HalfStructure() {
    G4VSolid* scintPlane = new G4Box("scintPlane",strucThick/2, scintWidthPlusPaint*4/3, scintWidthPlusPaint*4/3);
    G4VSolid* scintHole = new G4Box("scintHole",strucThick, scintWidthPlusPaint, scintWidthPlusPaint);
    G4VSolid* struc = new G4SubtractionSolid("struc", scintPlane, scintHole, 0, moveNull);
    G4VSolid* strucCut = new G4Box("scintCut", cutWidth, cutWidth, scintWidthPlusPaint*2);

    rotate = new G4RotationMatrix;
    rotate->rotateX(90.*deg);
    rotate->rotateZ(2*M_PI/16);
    move = G4ThreeVector(0.0,0.0,cutDist);
    G4VSolid* struc2 = new G4SubtractionSolid("struc2", struc, strucCut, rotate, move);

    rotate = new G4RotationMatrix;
    rotate->rotateX(90.*deg);
    rotate->rotateZ(-2*M_PI/16);
    move = G4ThreeVector(0.0,0.0,-cutDist);
    G4VSolid* struc3 = new G4SubtractionSolid("struc3", struc2, strucCut, rotate, move);

    rotate = new G4RotationMatrix;
    rotate->rotateZ(-2*M_PI/16);
    move = G4ThreeVector(0.0,cutDist,0.0);
    G4VSolid* struc4 = new G4SubtractionSolid("struc4", struc3, strucCut, rotate, move);

    rotate = new G4RotationMatrix;
    rotate->rotateZ(2*M_PI/16);
    move = G4ThreeVector(0.0,-cutDist,0.0);
    G4VSolid* struc5 = new G4SubtractionSolid("struc5", struc4, strucCut, rotate, move);

    G4VSolid* halfCut = new G4Box("halfCut", strucThick, scintWidthPlusPaint, scintWidthPlusPaint*4/3);
    move = G4ThreeVector(0.0,-scintWidthPlusPaint,0.0);
    rotate = new G4RotationMatrix;
    G4VSolid* struc6 = new G4SubtractionSolid("struc6", struc5, halfCut, rotate, move);

    return struc6;
}

void DetectionSystemSceptar2::SetStructureMaterial(G4String material) {

    fMaterialStructure = G4Material::GetMaterial(material);
    CheckMaterial(fMaterialStructure,material);
    materialStructureDefined = true;
}

void DetectionSystemSceptar2::SetStructureThickness(G4double thickness) {
    strucThick = thickness;

    if ( strucThick < 1.1*mm ) {
        strucThick = 1.1*mm;
        G4cout << "@@@ SCEPTAR2 Structure Thickness out of range. Defaulted to minimum 1.1mm." << G4endl;
        errors++;
    }
    if ( strucThick > 6.8*mm ) {
        strucThick = 6.8*mm;
        G4cout << "@@@ SCEPTAR2 Structure Thickness out of range. Defaulted to minimum 6.8mm." << G4endl;
        errors++;
    }

    cutWidth = 1.2*strucThick;
    strucGap = strucFrontGap + tan(2*M_PI/16)*strucThick/2;
    scintWidthPlusPaint = scintWidth + 2*paintThick; 
    cutDist = scintWidthPlusPaint + strucGap + cutWidth/cos(2*M_PI/16);
}

void DetectionSystemSceptar2::CheckMaterial(G4Material* material, G4String name) {
    if( !material ) {
        G4cout << " ----> Material " << name << " not found, cannot set detector structure material! " << G4endl;
        errors++;
        return;
    }
}

void DetectionSystemSceptar2::DefineCoords() {
    // GRIFFIN Coordinates, from DetectionSystemGriffinSuppressed.cc

    // theta
    fCoords[0][0] 	= 45.0;
    fCoords[1][0] 	= 45.0;
    fCoords[2][0] 	= 45.0;
    fCoords[3][0] 	= 45.0;
    fCoords[4][0] 	= 90.0;
    fCoords[5][0] 	= 90.0;
    fCoords[6][0] 	= 90.0;
    fCoords[7][0] 	= 90.0;
    fCoords[8][0] 	= 90.0;
    fCoords[9][0] 	= 90.0;
    fCoords[10][0] 	= 90.0;
    fCoords[11][0] 	= 90.0;
    fCoords[12][0] 	= 135.0;
    fCoords[13][0] 	= 135.0;
    fCoords[14][0] 	= 135.0;
    fCoords[15][0] 	= 135.0;
    // phi
    fCoords[0][1] 	= 67.5;
    fCoords[1][1] 	= 157.5;
    fCoords[2][1] 	= 247.5;
    fCoords[3][1] 	= 337.5;
    fCoords[4][1] 	= 22.5;
    fCoords[5][1] 	= 67.5;
    fCoords[6][1] 	= 112.5;
    fCoords[7][1] 	= 157.5;
    fCoords[8][1] 	= 202.5;
    fCoords[9][1] 	= 247.5;
    fCoords[10][1] 	= 292.5;
    fCoords[11][1] 	= 337.5;
    fCoords[12][1] 	= 67.5;
    fCoords[13][1] 	= 157.5;
    fCoords[14][1] 	= 247.5;
    fCoords[15][1] 	= 337.5;
    // yaw (alpha)
    fCoords[0][2] 	= 0.0;
    fCoords[1][2] 	= 0.0;
    fCoords[2][2] 	= 0.0;
    fCoords[3][2] 	= 0.0;
    fCoords[4][2] 	= 0.0;
    fCoords[5][2] 	= 0.0;
    fCoords[6][2] 	= 0.0;
    fCoords[7][2] 	= 0.0;
    fCoords[8][2] 	= 0.0;
    fCoords[9][2] 	= 0.0;
    fCoords[10][2] 	= 0.0;
    fCoords[11][2] 	= 0.0;
    fCoords[12][2] 	= 0.0;
    fCoords[13][2] 	= 0.0;
    fCoords[14][2] 	= 0.0;
    fCoords[15][2] 	= 0.0;
    // pitch (beta)
    fCoords[0][3] 	= -45.0;
    fCoords[1][3] 	= -45.0;
    fCoords[2][3] 	= -45.0;
    fCoords[3][3] 	= -45.0;
    fCoords[4][3] 	= 0.0;
    fCoords[5][3] 	= 0.0;
    fCoords[6][3] 	= 0.0;
    fCoords[7][3] 	= 0.0;
    fCoords[8][3] 	= 0.0;
    fCoords[9][3] 	= 0.0;
    fCoords[10][3] 	= 0.0;
    fCoords[11][3] 	= 0.0;
    fCoords[12][3] 	= 45.0;
    fCoords[13][3] 	= 45.0;
    fCoords[14][3] 	= 45.0;
    fCoords[15][3] 	= 45.0;
    // roll (gamma)
    fCoords[0][4] 	= 67.5;
    fCoords[1][4] 	= 157.5;
    fCoords[2][4] 	= 247.5;
    fCoords[3][4] 	= 337.5;
    fCoords[4][4] 	= 22.5;
    fCoords[5][4] 	= 67.5;
    fCoords[6][4] 	= 112.5;
    fCoords[7][4] 	= 157.5;
    fCoords[8][4] 	= 202.5;
    fCoords[9][4] 	= 247.5;
    fCoords[10][4] 	= 292.5;
    fCoords[11][4] 	= 337.5;
    fCoords[12][4] 	= 67.5;
    fCoords[13][4] 	= 157.5;
    fCoords[14][4] 	= 247.5;
    fCoords[15][4] 	= 337.5;
    // downstream zero degree position:
    fCoords[16][0] = 0.0;
    fCoords[16][1] = 67.5;
    fCoords[16][2] = 0.0;
    fCoords[16][3] = -90.0;
    fCoords[16][4] = 67.5;
    // upstream zero degree position:
    fCoords[17][0] = 180.0;
    fCoords[17][1] = 67.5;
    fCoords[17][2] = 0.0;
    fCoords[17][3] = 90.0;
    fCoords[17][4] = 67.5;
}