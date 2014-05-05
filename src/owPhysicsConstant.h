#ifndef OW_PHYSICS_CONSTANT_H
#define OW_PHYSICS_CONSTANT_H

#include "owOpenCLConstant.h"
#include <math.h>

#ifndef M_PI
#define M_PI 3.1415927f
#endif
//Sizes of the box containing simulated 'world'
//Sizes choice is realized this way because it should be proportional to smoothing radius h
const float h = 3.34f;
const float r0 = 0.5f * h; // distance between two boundary particle
#define XMIN 0
#define XMAX 116.9//51.0042063076//10.f * r0//33.4//24.64278211//24.64278211//33.4000081813//88.844//120.24f //11.69f
#define YMIN 0
#define YMAX 751.5//515.573339198//120.f * r0//320.198396832//368.398493//67.798493//33.4000081813//88.844//80.16f //11.69f
#define ZMIN 0
#define ZMAX 116.9//51.0042063076///10.f * r0//33.4//24.64278211//24.6427721146//66.7999918187//88.844//182.03f //11.69f

const float rho0 = 1000.0f;
const float stiffness = 0.75f;
const float hashGridCellSize = 2.0f * h;
const float hashGridCellSizeInv = 1.0f / hashGridCellSize;
const float mass = 3.25e-14f;//3.25e-14f;//0.0003f;//3.25e-14f;
const float simulationScale = 0.004f*pow(mass,1.f/3.f)/pow(0.00025f,1.f/3.f);
const float h___ = simulationScale * YMAX;
const float simulationScaleInv = 1.0f / simulationScale;
const float mu = 0.89e-4f;//7.0f;//;//4.0f;//why this value? Dynamic viscosity of water at 25 C = 0.89e-3 Pa*s
const float timeStep = 5.0e-06f;//5.0e-06f;//0.001f;//0.0005f;//0.0042f;// ATTENTION you should remember about time step: if it is larger 0.001 it can lead to 'explosion' of elastic matter objects
const float CFLLimit = 100.0f;
const float damping = 0.75f;

const float beta = timeStep*timeStep*mass*mass*2/(rho0*rho0);// B. Solenthaler's dissertation, formula 3.6 (end of page 30)
const float betaInv = 1.f/beta;
//const float Wpoly6Coefficient = 315.0f / ( 64.0f * M_PI * pow( h * simulationScale, 9.0f ) );
const double Wpoly6Coefficient = 315.0 / ( 64.0 * M_PI * pow( (double)(h*simulationScale), 9.0 ) );
//const float gradWspikyCoefficient= -45.0f / ( M_PI * pow( h * simulationScale, 6.0f ) );
const double gradWspikyCoefficient= -45.0 / ( M_PI * pow( (double)(h*simulationScale), 6.0 ) );
const double del2WviscosityCoefficient = -gradWspikyCoefficient;
const float gravity_x = 0.0f;
const float gravity_y = -9.8f;
const float gravity_z = 0.0f;
extern const float delta;
const int maxIteration = 3;
const int ELASTIC_CONNECTIONS_COUNT = 0;
const float surfTensCoeff = 0.0f;//-4.5e-10f;//-4.5e-10f;//-1.5e-15f;//-4.5e-10f;//-0.0013f;//
const float elasticityCoeff = 1.95e-05f;//10000.0f * mass

#endif // #ifndef OW_PHYSICS_CONSTANT_H
