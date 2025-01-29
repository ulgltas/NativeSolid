#include "../include/config.h"
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <algorithm>

using namespace std;

Config::Config(string filename):ConfigFileName(filename)
{
    // Setting defaults here
    NUMBER_DESIGN_VARIABLES = 1;
}

Config::~Config()
{

}

Config* Config::GetAddress()
{
    return this;
}

void Config::ReadConfig()
{
    string delimiter = "=";
    size_t pos;
    string text_line, option;
    ifstream InputFile;
    InputFile.open(ConfigFileName.c_str(), ios::in);
    if(InputFile.fail()){
      cerr << "Invalid configuration file name : " << ConfigFileName << endl;
      throw(-1);
    }
    char caract;
    while (getline(InputFile,text_line))
    {
        caract = text_line[0];
        if (caract == '%')
	{
	}
	else
	{
	pos = text_line.find(delimiter);
	option = text_line.substr(0,pos);
        text_line.erase(0,pos+delimiter.length());
        option.erase(remove(option.begin(), option.end(), ' '), option.end());
	text_line.erase(remove(text_line.begin(), text_line.end(),' '),text_line.end());
	if (option == "CSD_SOLVER") CSD_SOLVER = text_line;
    else if (option == "KIND_PROBLEM") KIND_PROBLEM = text_line;
    else if (option == "OBJ_FUNCTION") OBJ_FUNCTION = text_line;
	else if (option == "MESH_FILE") MESH_FILE = text_line;
	else if (option == "UNSTEADY_SIMULATION") UNSTEADY_SIMULATION = text_line;
	else if (option == "STRUCT_TYPE") STRUCT_TYPE = text_line;
        else if (option == "LINEARIZE") LINEARIZE = text_line;
        else if (option == "INTEGRATION_ALGO") INTEGRATION_ALGO = text_line;
        else if (option == "RESTART_SOL") RESTART_SOL = text_line;
    else if (option == "FIXED_DOF") FIXED_DOF = text_line;
	else if (option == "RESTART_FILE") RESTART_FILE = text_line;
	else if (option == "MOVING_MARKER") MOVING_MARKER = text_line;
	else if (option == "SPRING_MASS") SPRING_MASS = atof(text_line.c_str());
	else if (option == "INERTIA_CG") INERTIA_CG = atof(text_line.c_str());
	else if (option == "INERTIA_FLEXURAL") INERTIA_FLEXURAL = atof(text_line.c_str());
	else if (option == "SPRING_STIFFNESS") SPRING_STIFFNESS = atof(text_line.c_str());
	else if (option == "SPRING_DAMPING") SPRING_DAMPING = atof(text_line.c_str());
	else if (option == "TORSIONAL_STIFFNESS") TORSIONAL_STIFFNESS = atof(text_line.c_str());
	else if (option == "TORSIONAL_DAMPING") TORSIONAL_DAMPING = atof(text_line.c_str());
	else if (option == "CORD") CORD = atof(text_line.c_str());
	else if (option == "FLEXURAL_AXIS") FLEXURAL_AXIS = atof(text_line.c_str());
	else if (option == "GRAVITY_CENTER") GRAVITY_CENTER = atof(text_line.c_str());
        else if (option == "INITIAL_DISP") INITIAL_DISP = atof(text_line.c_str());
        else if (option == "INITIAL_ANGLE") INITIAL_ANGLE = atof(text_line.c_str());
	else if (option == "START_TIME") START_TIME = atof(text_line.c_str());
	else if (option == "DELTA_T") DELTA_T = atof(text_line.c_str());
	else if (option == "DELTA_ITER_WRITE") DELTAITERWRITE = atol(text_line.c_str());
	else if (option == "STOP_TIME") STOP_TIME = atof(text_line.c_str());
    else if (option == "OMEGA") OMEGA = atof(text_line.c_str());
    else if (option == "NUMBER_HARMONICS") NUMBER_HARMONICS = atol(text_line.c_str());
	else if (option == "RHO") RHO = atof(text_line.c_str());
    else if (option == "DESIGN_VARIABLE_KIND") DESIGN_VARIABLE_KIND = text_line;
    else if (option == "NUMBER_DESIGN_VARIABLES") NUMBER_DESIGN_VARIABLES = atol(text_line.c_str());
        else cout << "The option " + option + " is not recognized !" << endl;
	}
    }
    InputFile.close();

    if (CSD_SOLVER == "NATIVE") cout << "The Native solver has been chosen" << endl;
    else cout << "Cannot run the solver with other value than NATIVE for CSD_SOLVER option !" << endl;

    if (UNSTEADY_SIMULATION == "YES") cout << "Dynamic structure computation" << endl;
    else if (UNSTEADY_SIMULATION == "HARMONIC") cout << "Harmonic balance structure computation" << endl;
    else cout << "Static structure computation" << endl;

    if (STRUCT_TYPE == "SPRING_HOR" || STRUCT_TYPE == "SPRING_VER") cout << "Structural model is a plunging spring" << endl;
    else if (STRUCT_TYPE == "AIRFOIL") cout << "Structural model is a pitching plunging airfoil" << endl;
    else cout << "The specified structural model is not recognized or implemented yet !" << endl;
}

std::string Config::GetMeshFile()
{
    return MESH_FILE;
}

std::string Config::GetUnsteady()
{
    return UNSTEADY_SIMULATION;
}

std::string Config::GetCSDSolver()
{
    return CSD_SOLVER;
}

std::string Config::GetKindProblem()
{
    return KIND_PROBLEM;
}

std::string Config::GetFixedDof()
{
    return FIXED_DOF;
}

std::string Config::GetObjFunction()
{
    return OBJ_FUNCTION;
}

std::string Config::GetStructType()
{
    return STRUCT_TYPE;
}

std::string Config::GetLinearize(){
    return LINEARIZE;
}

std::string Config::GetIntegrationAlgo()
{
    return INTEGRATION_ALGO;
}

std::string Config::GetRestartSol()
{
    return RESTART_SOL;
}

std::string Config::GetRestartFile()
{
    return RESTART_FILE;
}

std::string Config::GetMovingMarker()
{
    return MOVING_MARKER;
}

std::string Config::GetDesignVariableKind()
{
    return DESIGN_VARIABLE_KIND;
}

double Config::GetStartTime()
{
    return START_TIME;
}

double Config::GetDeltaT(){
    return DELTA_T;
}

unsigned long Config::GetDeltaIterWrite(){
    return DELTAITERWRITE;
}

double Config::GetStopTime(){
    return STOP_TIME;
}

double Config::GetSpringStiffness(){
    return SPRING_STIFFNESS;
}

double Config::GetSpringMass(){
    return SPRING_MASS;
}

double Config::GetInertiaCG(){
    return INERTIA_CG;
}

double Config::GetInertiaFlexural(){
    return INERTIA_FLEXURAL;
}

double Config::GetSpringDamping(){
    return SPRING_DAMPING;
}

double Config::GetTorsionalStiffness(){
  return TORSIONAL_STIFFNESS;
}

double Config::GetTorsionalDamping(){
  return TORSIONAL_DAMPING;
}

double Config::GetCord(){
  return CORD;
}

double Config::GetFlexuralAxis(){
  return FLEXURAL_AXIS;
}

double Config::GetGravityCenter(){
    return GRAVITY_CENTER;
}

double Config::GetInitialDisp(){
  return INITIAL_DISP;
}

double Config::GetInitialAngle(){
  return INITIAL_ANGLE;
}

double Config::GetOmega(){
    return OMEGA;
}

unsigned long Config::GetNumberHarmonics(){
    return NUMBER_HARMONICS;
}

double Config::GetRho(){
    return RHO;
}

unsigned long Config::GetNumberDesignVariables()
{
    return NUMBER_DESIGN_VARIABLES;
}
