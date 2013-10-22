#include "owPhysicTests.h"
#include "owPhysicsFluidSimulator.h"
#include <iostream>

extern const int tested_id = 100;
extern int PARTICLE_COUNT;
int zero_point = 0;
void get_buffs(owPhysicsFluidSimulator * fluid_simulation, float * init_p, float * init_v, int);
float get_dist(float * from, float * to);
void zero_vel_buff( owPhysicsFluidSimulator * fluid_simulation);

void physic_friction_test(){
	int iterationNum = 300;
	int startIteration = iterationNum - 100;
	owPhysicsFluidSimulator * fluid_simulation;
	owHelper * helper;
	int i = 0;
	float * initial_position = new float[4];
	float * initial_velocity = new float[4];
	float * end_position = new float[4];
	float * end_velocity = new float[4];
	float * zero_position = new float[4];
	float * zero_velocity = new float[4];
	helper = new owHelper();
	fluid_simulation = new owPhysicsFluidSimulator(helper);
	while( i <= iterationNum ){
		get_buffs(fluid_simulation, zero_position, zero_velocity, zero_point);
		if(zero_velocity[1] > 0.f){
			std::cout <<"velocity:" << zero_velocity[0] << "\t" << zero_velocity[1] << "\t" << zero_velocity[2] << "\t" << std::endl;
			zero_vel_buff(fluid_simulation);
			
			i = 0;
		}
		if(i == 200){
			get_buffs(fluid_simulation, initial_position, initial_velocity, tested_id);
		}
		fluid_simulation->simulationStep();
		helper->refreshTime();
		i++;
	}
	get_buffs(fluid_simulation, end_position, end_velocity, tested_id);
	float result = get_dist(initial_position,end_position);
	std::cout << "Distance :" << result * simulationScale << " m" << std::endl;
	std::cout << "Time :" << (float)(100) * timeStep << " s" << std::endl;
	delete fluid_simulation;
	delete helper;
}
void zero_vel_buff( owPhysicsFluidSimulator * fluid_simulation){
	float * v_b = fluid_simulation->getVelocityBuffer();
	for(int i = 0; i < PARTICLE_COUNT; i++){
		if(int(v_b[4 * i + 3]) != BOUNDARY_PARTICLE){
			v_b[4 * i + 0] = 0.f;
			v_b[4 * i + 1] = 0.f;
			v_b[4 * i + 2] = 0.f;
			v_b[4 * i + 3] = v_b[4 * i + 3];
		}
	}
	fluid_simulation->putVelocityBuffer(v_b);
	//delete v_b;
}
float get_dist(float * from, float * to){
	return sqrt ( pow(to[0] - from[0], 2.f) + pow(to[1] - from[1], 2.f) + pow(to[2] - from[2], 2.f) + pow(to[3] - from[3], 2.f) );
}

void get_buffs(owPhysicsFluidSimulator * fluid_simulation, float * init_p, float * init_v, int id){
	float * p_b;
	float * v_b;
	p_b = fluid_simulation->getPositionBuffer();
	v_b = fluid_simulation->getVelocityBuffer();
	init_p[0] = p_b[ id * 4 + 0 ];
	init_p[1] = p_b[ id * 4 + 1 ];
	init_p[2] = p_b[ id * 4 + 2 ];
	init_p[3] = p_b[ id * 4 + 3 ];
	init_v[0] = v_b[ id * 4 + 0 ];
	init_v[1] = v_b[ id * 4 + 1 ];
	init_v[2] = v_b[ id * 4 + 2 ];
	init_v[3] = v_b[ id * 4 + 3 ];
	//delete v_b;
	//delete p_b;
}