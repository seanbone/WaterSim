#include <iostream>

#include "SimConfig.h"
#include "WaterSim.h"


class BasicSimRunner {
	SimConfig cfg_;
	WaterSim waterSim_;

public:
	explicit BasicSimRunner(std::string cfgFile) : waterSim_(cfg_) {
		if (!cfg_.readFile(cfgFile)) {
			throw std::runtime_error("Could not find configuration file. Aborting.");
		}

		waterSim_.updateParams(cfg_);
		waterSim_.resetMembers();
		std::cout << cfg_.toString() << std::endl;
	}

	void run() {
		std::cout << "\n\n****** Starting simulation..." << std::endl;
		std::cout << waterSim_.getNumParticles() << " particles." << std::endl;

		int nx, ny, nz;
		cfg_.getGridResolution(nx, ny, nz);
		std::printf("%d cells (%d x %d x %d)\n", nx*ny*nz, nx, ny, nz);

		int max = cfg_.getMaxSteps();
		for (int step = 0; max < 0 || step < max; step++) {
			waterSim_.advance();
		}
	}
};

int main(int argc, char* argv[]) {
	bool autostart = false;

	std::string defaultConfigFile =  "config.json";
	std::string configFile = defaultConfigFile;
	if (argc > 1) {
		std::string tmp = argv[1];
		if (tmp == "-h") {
			std::cout << "Usage: ./watersim-cli [-h] [-y] [<config>]" << std::endl;
			std::cout << "\n-h\t\t Print this help message and exit.\n";
			std::cout << "-y\t\t Autostart simulation.\n";
			std::cout << "<config>\t Location of JSON configuration file. Default = '"
					  << defaultConfigFile << "'." << std::endl;
			return 0;
		} else if (tmp == "-y")
			autostart = true;
		else
			configFile = argv[1];
		if (argc > 2)
			configFile = argv[2];
	}

	auto simRunner = BasicSimRunner(configFile);

	if (!autostart) {
		std::cout << "\n\nStart simulation with above configuration? [y/n] ";
		char ans;
		std::cin >> ans;
		if (ans != 'y') {
			std::cout << "Aborting." << std::endl;
			return 0;
		}
	}

	simRunner.run();

	return 0;
}
