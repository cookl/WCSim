#ifndef INTERACTION_H
#define INTERACTION_H
//This defines a class object to contain interaction details to be associated with WCSimTrajectory
//this is a non-root version of WCSimRootInteraction defined to be used before the writing-out portion 

#include<string>
#include<vector>
class WCSimInteraction{

public:
	WCSimInteraction() {}

	WCSimInteraction(int inPiD, 
      float inDir[3], 
      float inMom[3], 
      int inTrackID, 
      float intVertex[3], 
      int intCode, 
      std::string intName)
{

  fInPid = inPiD;      
  fInTrackID = inTrackID;
  fIntCode = intCode;
  fIntName = intName;
  fNDaughters =0;
  // fNDaughters = nDaughters;

  for(int i=0; i<3; i++){
    fInDir[i] = inDir[i];
    fInMom[i] = inMom[i];
    fIntVertex[i] = intVertex[i];

  }
  return;
}

void AddDaughter(int pID, 
        int trackID, 
        float mom[3]){
  
  fNDaughters++;
  fDaughterPID.push_back(pID);
  fDaughterTrackID.push_back(trackID);
  fDaughterMomX.push_back(mom[0]);
  fDaughterMomY.push_back(mom[1]);
  fDaughterMomZ.push_back(mom[2]);
  return;

}

	int   getInPid(){return fInPid;} //input particle ID       
	float getInDir(int i){ return fInDir[i];} //input particle direction
	float getInMom(int i){ return fInMom[i];} //input particle Momentum
	float getInTrackID(){return fInTrackID;} //input particle track ID

	//interaction info
	float getIntVertex(int i){return fIntVertex[i];} //interaction vertex
	int   getIntCode(){return fIntCode;} //interaction code
	std::string getIntName(){return fIntName;} //interaction name
	int getNDaughters(){return fNDaughters;}
	
	int  GetDaughterPID(int i=0) const { return (i<fNDaughters) ? fDaughterPID[i] : 0;} 
  int  GetDaughterTrackID(int i=0) const { return (i<fNDaughters) ? fDaughterTrackID[i] : 0;} 
  float  GetDaughterMomX(int i=0) const { return (i<fNDaughters) ? fDaughterMomX[i] : 0;} 
  float  GetDaughterMomY(int i=0) const { return (i<fNDaughters) ? fDaughterMomY[i] : 0;} 
  float  GetDaughterMomZ(int i=0) const { return (i<fNDaughters) ? fDaughterMomZ[i] : 0;} 


private:
  //input particle info
  int   fInPid; //input particle ID       
  float fInDir[3]; //input particle direction
  float fInMom[3]; //input particle Momentum
  int fInTrackID; //input particle track ID

  //interaction info
  float fIntVertex[3]; //interaction vertex
  int fIntCode; //interaction code
  std::string fIntName; //interaction name

  //output daughters info
  int fNDaughters;
  std::vector<int> fDaughterPID;
  std::vector<int> fDaughterTrackID;
  std::vector<float> fDaughterMomX;
  std::vector<float> fDaughterMomY;
  std::vector<float> fDaughterMomZ;


};

#endif