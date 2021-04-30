#include <iostream>
#include <netcdf>

// We are writing 2D data, a 6 x 12 grid
constexpr int nx = 6;
constexpr int ny = 12;

// Return this in event of a problem
constexpr int nc_err = 2;

int main() {
	// The default behavior of the C++ API is to throw an exception if
	// an error occurs
	try {
		// This is the data array we will write. It will just be filled
		// with a progression of numbers for this example.
		int dataOut[nx][ny];

		// Create some pretend data. If this wasn't an example program, we
		// would have some real data to write, for example, model output.
		for (int i = 0; i < nx; i++) {
			for (int j = 0; j < ny; j++) {
				dataOut[i][j] = i * ny + j;
			}
		}

		// Create the file. The Replace parameter tells netCDF to overwrite
		// this file, if it already exists.
		netCDF::NcFile dataFile("simple_xy.nc", netCDF::NcFile::replace);

		// Create netCDF dimensions
		auto xDim = dataFile.addDim("x", nx);
		auto yDim = dataFile.addDim("y", ny);

		// Define the variable. The type of the variable in this case is
		// ncInt (32-bit integer)
		auto data = dataFile.addVar("data", netCDF::ncInt, {xDim, yDim});

		// Write the data to the file. Although netCDF supports reading
		// and writing subsets of data, in this case we write all the data
		// in one operation.
		data.putVar(dataOut);

		// The file will be automatically close when the NcFile object goes
		// out of scope. This frees up any internal netCDF resources
		// associated with the file, and flushes any buffers.
	} catch (netCDF::exceptions::NcException &e) {
		std::cout << e.what() << std::endl;
		return nc_err;
	}

	// Now read the data back in
	try {
		// This is the array we will read into
		int dataIn[nx][ny];

		// Open the file for read access
		netCDF::NcFile dataFile("simple_xy.nc", netCDF::NcFile::read);

		// Retrieve the variable named "data"
		auto data = dataFile.getVar("data");
		if (data.isNull())
			return nc_err;
		data.getVar(dataIn);

		// Check the values.
		for (int i = 0; i < nx; i++) {
			for (int j = 0; j < ny; j++) {
				if (dataIn[i][j] != i * ny + j) {
					return nc_err;
				}
			}
		}
	} catch (netCDF::exceptions::NcException &e) {
		std::cout << e.what() << std::endl;
		return nc_err;
	}
}
