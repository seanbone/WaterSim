#ifndef NCREADER_H
#define NCREADER_H

#include <iostream>
#include <netcdf>
#include <string>
#include <vector>
#include <algorithm>
#include <cmath>
#include "Particles.h"
#include "Mac3d.h"
#include "SimConfig.h"

using namespace netCDF;
using namespace netCDF::exceptions;


class NcReader {
	
	public:
	
	/** Constructor
	 * Params:
	 * - filePath is a string containing the path to the file to read
	 * - cfgPath is a string containing the path to the config.json file to be 
	 *   used
	 */
	NcReader( const std::string& filePath, const std::string& cfgPath );
	
	
	/** Destructor
	 */	
	~NcReader();
	

	/** Read the specified scalar variable from _filePath 
	 * Params:
     * - varName is a string containing the name of the variable to be read 
	 */
	template<typename T>
	T readScalar( const std::string& varName ){

		try{

			// Open the file to read
			NcFile dataFile(_filePath, NcFile::read);
			
			// Read the variable
			NcVar var = dataFile.getVar(varName);
			
			// Check if the variable was read correctly (isNull() returns 
			// true if the NcVar object is empty)
			if (var.isNull()){
				throw ("NcVar object " + varName + " is empty");
			}
			
			// Store data into destination
			T destination;
			var.getVar(&destination);

			return destination;

		} catch(NcException& e){
			std::cout << "--- READ SCALAR FAILURE: " + varName + " ---" << std::endl;
			std::cout << e.what() << std::endl;

			return 0;
		}
	}
	

	/** Read the specified variable from _filePath 
	 * Params:
	 * - varName is a string containing the name of the variable to be read 
	 * - destination is a pointer to an array of double where the data will be 
	 *   stored
	 */
	template<typename T>
	void read( const unsigned breakPt, const std::string& varName, T* destination ){
		
		try{
			
			// Open the file to read
			NcFile dataFile(_filePath, NcFile::read);
			
			// Read the variable
			NcVar var = dataFile.getVar(varName);

			std::size_t n = var.getDim(1).getSize();
			
			// Check if the variable was read correctly (isNull() returns 
			// true if the NcVar object is empty)
			if (var.isNull()){
				throw ("NcVar object " + varName + " is empty");
			}
			
			// Fill the array with read data
			std::vector<std::size_t> start, count;
			start.push_back(breakPt);
			start.push_back(0);
			count.push_back(1);
			count.push_back(n);
			
			// Store data into destination
			var.getVar(start, count, destination);

		} catch(NcException& e){
			std::cout << "--- READ FAILURE: " + varName + " ---" << std::endl;
			std::cout << e.what() << std::endl;
		}
	}


	/** Convert the linear arrays read from _filePath to the data structures of 
	 *  the FLIP algorithm implementation (Particles and Mac3d)
	 */
	void toFlipStructures();


	/** Read all the arrays present in the reference file for a specific break 
	 *  point and store them in the relative members
	 * Params:
	 * - breakPt is the number of the desired break point (see FLIP::step_FLIP 
	 *   method)
	 */
	void readAll( unsigned breakPt );


	/** Check if the model state validates
	 */
	void validate();


	/** Get the size of a dimension
	 * Params:
	 * - dimName is a string containing the name of the desired dimension
	 */
	unsigned getDim( const std::string& dimName );


	// Dimensions
	unsigned num_particles;
	unsigned timestep;

	// Reference linear arrays
	double* uMAC;
	double* vMAC;
	double* wMAC;
	double* uStar;
	double* vStar;
	double* wStar;
	double* pMAC;
	bool* fluid_cells;
	bool* solid_cells;

	Particles* particles;
	Particles* referenceParticles;

	// Pointers to old data structures used for the "test simulation"
	Mac3d* MACGrid;

	// Initial configuration object
	SimConfig cfg;
	
	
	private:

	/** Compute the relative error between a reference and a given value
	 * Params:
	 * - ref is the reference value
	 * - value is the value resulting from the "test simulation"
	 */
	double rErr( const double ref, const double value );


	/** Output an error message specifying which variable does not validate and 
	 *  the relative error
	 * Params:
	 * - flag determines if the output message should be outputted
	 * - varName is a string containing the name of the non-validating variable
	 * - maxErr is the maximum relative error of the specified variable
	 * - tol is the tolerance to achieve validation against the reference values
	 */
	void outputMessage( const bool flag, const std::string& varName, const double maxErr, const double tol );


    // String containing the path to the file to read
    const std::string _filePath;
	
	// MAC grid dimensions
	int _n;
	int _m;
	int _l;

	// MAC grid dimensions (in meters)
	double _dx;
	double _dy;
	double _dz;
};

#endif