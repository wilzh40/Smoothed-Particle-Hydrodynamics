#include "owPhysicTests.h"
#include "owPhysicsFluidSimulator.h"
#include <iostream>
#include <fstream>
#include <string>

extern const int tested_id = 100;
extern int PARTICLE_COUNT;
int zero_point = 0;
void get_buffs(owPhysicsFluidSimulator * fluid_simulation, float * init_p, float * init_v, int);
float get_dist(float * from, float * to);
void zero_vel_buff( owPhysicsFluidSimulator * fluid_simulation);
void calc_mass_center(owPhysicsFluidSimulator * fluid_simulator, float *);
void calc_mass_center_v(owPhysicsFluidSimulator * fluid_simulator, float * v);
float get_dist_3(float * from, float * to);
float mass_ = 0.f;
float true_time = 0.f;
float gravity = 9.81;
float eps = 0.0001f;
void physic_friction_test(){
	int start_v_iteration = 550;
	int end_iteration = start_v_iteration + 200;
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
	bool first = true;
	while(i <= end_iteration){
		get_buffs(fluid_simulation, zero_position, zero_velocity, zero_point);
		if(zero_velocity[1] > 0.f && i < start_v_iteration){
			std::cout <<"velocity:" << zero_velocity[0] << "\t" << zero_velocity[1] << "\t" << zero_velocity[2] << "\t" << std::endl;
			zero_vel_buff(fluid_simulation);
		}else{
			if(i > start_v_iteration && first){
				get_buffs(fluid_simulation, initial_position, initial_velocity, tested_id);
				first = false;
			}
		}
		fluid_simulation->simulationStep();
		helper->refreshTime();
		i++;
	}
	get_buffs(fluid_simulation, end_position, end_velocity, tested_id);
	float result = get_dist(initial_position,end_position);
	std::cout << "Distance :" << result * simulationScale << " m" << std::endl;
	std::cout << "Time :" << (float)(end_iteration - start_v_iteration) * timeStep << " s" << std::endl;
	delete fluid_simulation;
	delete helper;
}
void physic_suvat_test(){
	int iteration_start = 0;
	int iteration_stop = iteration_start + 100;
	owPhysicsFluidSimulator * fluid_simulation;
	owHelper * helper;
	helper = new owHelper();
	fluid_simulation = new owPhysicsFluidSimulator(helper);
	float * initial_position = new float[3];
	float * end_position = new float[3];
	float * end_velocity = new float[3];
	int i = iteration_start;
	calc_mass_center(fluid_simulation, initial_position);
	while(i < iteration_stop ){
		fluid_simulation->simulationStep();
		helper->refreshTime();
		i++;
	}
	calc_mass_center(fluid_simulation, end_position);
	calc_mass_center_v(fluid_simulation, end_velocity);
	float result = get_dist_3(initial_position,end_position);
	calc_mass_center_v(fluid_simulation, end_velocity);
	std::cout << "Distance :" << result * simulationScale << " m" << std::endl;
	std::cout << "Time :" << (float)(iteration_stop - iteration_start) * timeStep << " s" << std::endl;
	std::cout << "Mass :" << mass_ << " kg" << std::endl;
	std::cout << "Velocity :" << sqrt ( pow( end_velocity [0], 2.f ) + pow ( end_velocity [1], 2.f ) + pow ( end_velocity[2], 2.f ) ) << " m/s" << std::endl;
	delete fluid_simulation;
	delete helper;
	return;
}

void gravity_test_1(){
	int iteration_start = 0;
	int iteration_stop = iteration_start + 100;
	owPhysicsFluidSimulator * fluid_simulation;
	owHelper * helper;
	helper = new owHelper();
	fluid_simulation = new owPhysicsFluidSimulator(helper);
	float * initial_position = new float[3];
	float * end_position = new float[3];
	float * end_velocity = new float[3];
	float * prev_velocity = new float[3];
	prev_velocity[0] = 0.f;
	prev_velocity[1] = -9.81f;
	prev_velocity[2] = 0.f;
	int i = iteration_start;
	calc_mass_center(fluid_simulation, initial_position);
	float s = 0.f;
	float s_b = 0.f;

	std::string file_name = "./test_graphs/gravity_test_single_particle.txt";
	std::ofstream out_f (file_name.c_str(), std::ofstream::out);
	if( out_f.is_open() )
	{
		out_f << "Distance" << "\t" << "Exact value of Time" << "\t"<<"Simulation Time\n" ;
		float scal_velocity_val = 0.f;
		while(s * simulationScale <= 0.15){
			s_b = s;
			fluid_simulation->simulationStep();
			helper->refreshTime();
			calc_mass_center(fluid_simulation, end_position);
			s = get_dist_3(initial_position,end_position);
			i++;
			true_time = sqrt( ( 2 * s * simulationScale ) / gravity );
			out_f << s * simulationScale << "\t" << true_time << "\t"<< (float)i * timeStep << "\n" ;
			calc_mass_center_v(fluid_simulation, end_velocity);
			float scal_new = end_velocity[0] * prev_velocity[0] + end_velocity[1] * prev_velocity[1] + end_velocity[2] * prev_velocity[2];
			if( scal_velocity_val * scal_new < 0){
				std::cout << "Velocity :" << scal_new << " m/s" << std::endl;
				break;
			}
			scal_velocity_val = scal_new;
			prev_velocity[0] = end_velocity[0];
			prev_velocity[1] = end_velocity[1];
			prev_velocity[2] = end_velocity[2];
		}

	}
	out_f.close();

	//std::cout << "Distance :" << s_b * simulationScale << " m" << std::endl;
	std::cout << "Distance :" << s * simulationScale << " m" << std::endl;
	std::cout << "Time :" << (float)i * timeStep << " s" << std::endl;
	std::cout << "Mass :" << mass_ << " kg" << std::endl;
	std::cout << "Velocity :" << sqrt ( pow( end_velocity [0], 2.f ) + pow ( end_velocity [1], 2.f ) + pow ( end_velocity[2], 2.f ) ) << " m/s" << std::endl;
	true_time = sqrt( ( 2 * s * simulationScale ) / gravity );
	if( abs(true_time - (float)i * timeStep) <= eps ){
		std::cerr << "TEST GRAVITY 1 PASS" << std::endl;
	}else{
		std::cerr << "TEST GRAVITY 1 FAIL " << "difference btween real t and calculated is = " << abs(true_time - (float)i *timeStep) << std::endl;
	}
	delete fluid_simulation;
	delete helper;
	return;
}
void calc_mass_center(owPhysicsFluidSimulator * fluid_simulator, float * mass_center){
	float * p_b;
	p_b = fluid_simulator->getPositionBuffer();
	mass_center[0] = 0.f;
	mass_center[1] = 0.f;
	mass_center[2] = 0.f;
	float p_count = 0.f;
	for(int i = 0; i < PARTICLE_COUNT; i++){
		if(int(p_b[ 4 * i + 3 ]) != BOUNDARY_PARTICLE){
			mass_center[0] += p_b[4 * i + 0];
			mass_center[1] += p_b[4 * i + 1];
			mass_center[2] += p_b[4 * i + 2];
			p_count++;
		}
	}
	mass_center[0] /= p_count;
	mass_center[1] /= p_count;
	mass_center[2] /= p_count;
	mass_ = mass * p_count;
}
void calc_mass_center_v(owPhysicsFluidSimulator * fluid_simulator, float * v){
	float * v_b;
	v_b = fluid_simulator->getVelocityBuffer();
	v[0] = 0.f;
	v[1] = 0.f;
	v[2] = 0.f;
	float p_count = 0.f;
	for(int i = 0; i < PARTICLE_COUNT; i++){
		if(int(v_b[ 4 * i + 3 ]) != BOUNDARY_PARTICLE){
			v[0] += v_b[4 * i + 0];
			v[1] += v_b[4 * i + 1];
			v[2] += v_b[4 * i + 2];
			p_count++;
		}
	}
	v[0] /= p_count;
	v[1] /= p_count;
	v[2] /= p_count;
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
float get_dist_3(float * from, float * to){
	return sqrt ( pow(to[0] - from[0], 2.f) + pow(to[1] - from[1], 2.f) + pow(to[2] - from[2], 2.f) );
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
