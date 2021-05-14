#ifndef NCWRITER_H
#define NCWRITER_H

#include <iostream>
#include <netcdf>
#include <string>
#include <vector>
#include "Particles.h"
#include "Mac3d.h"

using namespace netCDF;
using namespace netCDF::exceptions;


class NcWriter {
	
	public:
	
	/** Constructor
	 * Params:
	 * - filePath is a string containing the path to the file to write.
	 * - break_points is the total number of break points where the data needs 
	 *   to be logged
	 * - num_particles is the total number of particles
	 * - n is the x-dimension of the MAC grid
	 * - m is the y-dimension of the MAC grid
	 * - l is the z-dimension of the MAC grid
	 */
	NcWriter( const std::string& filePath, 
			  const unsigned break_points, 
			  const unsigned num_particles, 
			  const unsigned n, 
			  const unsigned m, 
			  const unsigned l );
	
	
	/** Destructor
	 */	
	~NcWriter();


	/** Add a variable to the reference file
	 * Params:
	 * - varName is a string containing the name of the variable to add
	 * - dimName is a string containing the name of the dimension of the 
	 *   variable to add
	 * - ncType is the NcType of the specified variable
	 */
	void addVar( const std::string& varName, const std::string& dimName, NcType ncType );


	/** Write the value of a scalar variable to the reference file
	 * Params:
	 * - varName is a string containing the name of the desired variable
	 * - value is the scalar value of the specified variable
	 * - ncType is the NcType of the specified variable
	 */
	template<typename T>
	void writeScalar( const std::string& varName, const T* const value, NcType ncType ){

		try{
			
			// Open the file to write
			NcFile dataFile(_filePath, NcFile::write);
			
			// Add the new variable to the file
			NcVar var = dataFile.addVar(varName, ncType);

			// Write the value to the file
			var.putVar(value);
			
			// Define the attributes of the new variable
			// var.putAtt("units", "m");
			
		} catch(NcException& e){
			std::cout << "--- WRITE SCALAR FAILURE: " + varName + " ---" << std::endl;
			std::cout << e.what() << std::endl;
		}
	}


	/** Create a new variable named varName and write the corresponding data 
	 *  into it.
	 * Params:
	 * - data is a pointer to the first element of an array containing the data 
	 *   of the new variable.
	 */
	template<typename T>
	void write( const unsigned breakPt, const std::string& varName, const T* const data ){
		
		try{
			
			// Open the file to write
			NcFile dataFile(_filePath, NcFile::write);
			
			// Get the desired variable in the file
			NcVar var = dataFile.getVar(varName);

			std::size_t n = var.getDim(1).getSize();
			
			std::vector<std::size_t> start;
			std::vector<std::size_t> count;
			start.push_back(breakPt);
			start.push_back(0);
			count.push_back(1);
			count.push_back(n);

			// Store the data in the new variable
			var.putVar(start, count, data);
			
		} catch(NcException& e){
			std::cout << "--- WRITE FAILURE: " + varName + " ---" << std::endl;
			std::cout << e.what() << std::endl;
		}
	}
	

	/** Write all the relevant arrays of the FLIP::step_FLIP() function to the 
	 *  reference file for a specific break point
	 * - breakPt is the number of the desired break point (see FLIP::step_FLIP 
	 *   method)
	 * - The remaining arguments contain the linear arrays containing the data 
	 *   to be written
	 */
	void writeAll( const unsigned breakPt, 
				   const Particles& particles, 
				   const Mac3d* const MACGrid );
	

	private:
	
	// String containing the path to the file to write
	const std::string _filePath;
};

#endif
