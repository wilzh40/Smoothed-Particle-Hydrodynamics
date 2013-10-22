#ifndef OW_PHYSICS_SIMULATOR_H
#define OW_PHYSICS_SIMULATOR_H

#include "owPhysicsConstant.h"
#include "owHelper.h"
#include "owOpenCLSolver.h"

class owPhysicsFluidSimulator
{
public:
	owPhysicsFluidSimulator(void);
	owPhysicsFluidSimulator(owHelper * helper);
	~owPhysicsFluidSimulator(void);
	float * getPositionBuffer() {  return positionBuffer; };
	float * getVelocityBuffer() { ocl_solver->read_velocity_b( velocityBuffer ) ; return velocityBuffer; };
	float * getDensityBuffer() { ocl_solver->read_density_b( densityBuffer ); return densityBuffer; };
	void putVelocityBuffer(float * v_b);
	unsigned int * getParticleIndexBuffer() { ocl_solver->read_particleIndex_b( particleIndexBuffer ); return particleIndexBuffer; };
	//TODO helper functions delete after fix!!
	float * getElasticConnections() { return elasticConnections; };
	double simulationStep();
private:
	owOpenCLSolver * ocl_solver;
	float * positionBuffer;
	float * velocityBuffer;
	float * elasticConnections;
	//Helper buffers
	float * densityBuffer;
	unsigned int * particleIndexBuffer;
	float * accelerationBuffer;//TODO REMOVE after fixing
	owHelper * helper;
};

#endif //OW_PHYSICS_SIMULATOR_H