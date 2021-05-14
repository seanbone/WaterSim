#include "NcWriter.h"

NcWriter::NcWriter( const std::string& filePath, 
					const unsigned break_points, 
					const unsigned num_particles, 
					const unsigned n, 
					const unsigned m, 
					const unsigned l ):
					_filePath(filePath)
{

	try{
		
		NcFile dataFile(_filePath, NcFile::replace);

		dataFile.addDim("break_points", break_points);
		dataFile.addDim("num_particles", num_particles);
		dataFile.addDim("x_mac_faces", (n+1) *  m    *  l   );
		dataFile.addDim("y_mac_faces",  n    * (m+1) *  l   );
		dataFile.addDim("z_mac_faces",  n    *  m    * (l+1));
		dataFile.addDim("mac_centers",  n * m * l);

		dataFile.close();

	} catch(NcException& e){
		std::cout << "--- NCWRITER CONSTRUCTOR FAILURE ---" << std::endl;
		std::cout << e.what() << std::endl;
	}

}


NcWriter::~NcWriter(){
	
}


void NcWriter::addVar( const std::string& varName, const std::string& dimName, NcType ncType ){
	
	try{
		
		// Open the file to write
		NcFile dataFile(_filePath, NcFile::write);
		
		// Get the dimensions
		NcDim bpDim = dataFile.getDim("break_points");
		NcDim nDim  = dataFile.getDim(dimName);
		
		// Store the dimensions into a vector in the correct order
		std::vector<NcDim> dims;
		dims.push_back(bpDim);
		dims.push_back(nDim);
		
		// Add the new variable to the file
		NcVar var = dataFile.addVar(varName, ncType, dims);
		
		// Define the attributes of the new variable
		// var.putAtt("units", "m");
		
	} catch(NcException& e){
		std::cout << "--- ADD VAR FAILURE: " + varName + " ---" << std::endl;
		std::cout << e.what() << std::endl;
	}
}


void NcWriter::writeAll( const unsigned breakPt, 
						 const Particles& particles, 
						 const Mac3d* const MACGrid )
{

	write(breakPt, "x", particles.x);
	write(breakPt, "y", particles.y);
	write(breakPt, "z", particles.z);
	write(breakPt, "u", particles.u);
	write(breakPt, "v", particles.v);
	write(breakPt, "w", particles.w);
	write(breakPt, "uMAC", MACGrid->pu_);
	write(breakPt, "vMAC", MACGrid->pv_);
	write(breakPt, "wMAC", MACGrid->pw_);
	write(breakPt, "uStar", MACGrid->pu_star_);
	write(breakPt, "vStar", MACGrid->pv_star_);
	write(breakPt, "wStar", MACGrid->pw_star_);
	write(breakPt, "pMAC", MACGrid->ppressure_);
	write(breakPt, "fluid_cells", MACGrid->pfluid_);
	write(breakPt, "solid_cells", MACGrid->psolid_);

}
