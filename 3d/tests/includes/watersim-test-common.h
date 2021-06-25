/*
 * A header file to contain functions common to all WaterSim tests.
 * Put any utility functions you may need to share between tests here.
 */
#ifndef WATERSIM_WATERSIM_TEST_COMMON_H
#define WATERSIM_WATERSIM_TEST_COMMON_H

#include <cassert>
#include <string>
#include <limits>
const double interpolation_tolerance = std::numeric_limits<double>::epsilon();

std::string validation_data_dir = "/validation_data/";
#ifdef CMAKE_SOURCE_DIR
	// Macros to convert preprocessor variable to string, see
	//  https://stackoverflow.com/questions/240353/convert-a-preprocessor-token-to-a-string
	#define STRINGIFY(x) #x
	#define TOSTRING(x) STRINGIFY(x)

	std::string validation_data_path = std::string(TOSTRING(CMAKE_SOURCE_DIR)).append(validation_data_dir);
#else
	std::string validation_data_path = validation_data_dir;
#endif

std::string validation_data_ref = validation_data_path + "ref.nc";
std::string validation_data_cfg = validation_data_path + "validation-config.json";

#endif //WATERSIM_WATERSIM_TEST_COMMON_H
