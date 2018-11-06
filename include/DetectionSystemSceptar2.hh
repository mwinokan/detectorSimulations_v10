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
//
// $Id: DetectorConstruction.hh,v 1.1 2010-10-18 15:56:17 maire Exp $
// GEANT4 tag $Name: geant4-09-04-patch-02 $
//
// 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#ifndef DETECTIONSYSTEMSCEPTAR2_HH
#define DETECTIONSYSTEMSCEPTAR2_HH

#include "G4VUserDetectorConstruction.hh"
#include "globals.hh"

class DetectionSystemSceptar2   
{
public:
    DetectionSystemSceptar2();                                                 // Initialisation
    ~DetectionSystemSceptar2();                                                // Destruction
    void SetConfig(G4int configuration, G4bool buildStruc, G4bool buildSipm);  // Variable passing
    void SetStructureMaterial(G4String material);                              // Variable passing
    void SetStructureThickness(G4double thickness);                            // Variable passing
    void BuildDetector();                                                      // Build the detector
    void PlaceDetector(G4LogicalVolume* expHallLog);                           // Place the detector in logical volume
    void Debrief();
private:
    G4LogicalVolume* fExpHallLog;               // World
    G4LogicalVolume* fSipmShroudLog;            // SiPM shroud
    G4LogicalVolume* fSipmActiveLog;            // SiPM active area
    G4LogicalVolume* fSipmConductLog;           // SiPM conductor
    G4LogicalVolume* fScintillatorUpStrLog;     // Upstream Scintillator
    G4LogicalVolume* fScintillatorDnStrLog;     // Downstream Scintillator
    G4LogicalVolume* fStructureLog;             // Holding structure unit
    G4LogicalVolume* fHalfStructureLog;         // Holding structure half unit
    G4LogicalVolume* fStructureMountLog;        // Holding structure mount

    G4AssemblyVolume* fUpAssembly;              // Upstream Scintillators + SiPMs
    G4AssemblyVolume* fDnAssembly;              // Downstream Scintillators + SiPMs
    G4AssemblyVolume* fStructureAssembly;       // Holding structure Unit
    G4AssemblyVolume* fHalfStructureAssembly;   // Holding structure Half Unit
    G4AssemblyVolume* fStructureMountAssembly;  // Holding structure mount

    void DefineCoords();            // Define GRIFFIN coordinates
    void ConstructScintillator();   // Construct scintillator and add to assembly
    void ConstructSipm();           // Construct SiPM and add to assembly
    void ConstructStructure();      // Construct structure and add to assembly
    void ConstructStructureMount(); // Construct structure mount and add to assembly
    
    void CheckMaterial(G4Material* material, G4String name);

    G4Box* SquareScintillator();    // Square scintillator
    G4VSolid* SquareStructure();    // Structure unit
    G4VSolid* HalfStructure();      // Half Structure unit

    G4Material* fMaterialScint;
    G4Material* fMaterialStructure;

    // translation functions for GRIFFIN coordinates
    G4double TransX(G4double x, G4double y, G4double z, G4double theta, G4double phi);
    G4double TransY(G4double x, G4double y, G4double z, G4double theta, G4double phi);
    G4double TransZ(G4double x, G4double z, G4double theta);

    // GRIFFIN coordinate array
    double fCoords[20][5];
    int fBuildPos[20];
    
    G4int errors;

    G4double scintThick;
    G4double scintWidth;
    G4double scintDist;
    G4double strucGap;
    G4double strucFrontGap;
    G4double cutWidth;
    G4double strucThick;
    G4double cutDist;
    G4double mountThick;
    
    bool buildScint;
    bool buildScintTriangles;
    bool buildSipm;
    bool buildStruc;
    bool buildStrucMount;
    bool materialStructureDefined;

    // rotation matrices
    G4RotationMatrix* rotateNull;
    G4RotationMatrix* rotate;

    // translation vectors
    G4ThreeVector moveNull;
    G4ThreeVector move;
    G4ThreeVector moveB;
    G4ThreeVector moveG;
    G4ThreeVector moveR;
    G4ThreeVector moveW;
};

#endif