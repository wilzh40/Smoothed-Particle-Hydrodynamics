#include "owWorldSimulation.h"
#include "owPhysicTests.h"
#include <stdio.h>

int main(int argc, char **argv)
{
	if(argc == 1)
		run( argc, argv);
	else{
		for(int i = 1; i<argc; i++){
			if(strncmp(argv[i], "-no_g", 5) == 0)//run without graphics
				run( argc, argv, false);
			if(strncmp(argv[i], "-test", 5) == 0)
				//gravity_test_1();
				physic_friction_test();
		}
	}
	return 0;
}
