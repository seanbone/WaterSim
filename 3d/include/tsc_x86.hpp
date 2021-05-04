/* ==================== Assume GCC under Linux ===================== */
#include <cassert>
#include <map>
#include <string>
#include <iostream>
#include <fstream>
#include <nlohmann/json.hpp>
namespace tsc {
	typedef unsigned long long myInt64;
	typedef unsigned int INT32;

	/* This is the RDTSC timer.
	 * RDTSC is an instruction on several Intel and compatible CPUs that Reads the
	 * Time Stamp Counter. The Intel manuals contain more information.
	 */

#define COUNTER_LO(a) ((a).int32.lo)
#define COUNTER_HI(a) ((a).int32.hi)
#define COUNTER_VAL(a) ((a).int64)

#define COUNTER(a) \
		((unsigned long long)COUNTER_VAL(a))

#define COUNTER_DIFF(a,b) \
		(COUNTER(a)-COUNTER(b))

		typedef union
		{       myInt64 int64;
				struct {INT32 lo, hi;} int32;
		} tsc_counter;

#define RDTSC(cpu_c) \
	  __asm__ __volatile__ ("rdtsc" : "=a" ((cpu_c).int32.lo), "=d"((cpu_c).int32.hi))
#define CPUID() \
	  __asm__ __volatile__ ("cpuid" : : "a" (0) : "bx", "cx", "dx" )


	inline myInt64 start(void) {
		tsc_counter start;
		CPUID();
		RDTSC(start);
		return COUNTER_VAL(start);
	}

	inline myInt64 stop(myInt64 start) {
		tsc_counter end;
		RDTSC(end);
		CPUID();
		return COUNTER_VAL(end) - start;
	}
	inline std::string thousandsep(std::string str, char character = '\'') {
		for (auto it = str.rbegin() + 3; it < str.rend(); it+=3) {
			str.insert(it.base(), character);
		}
		return str;
	}

	using json = nlohmann::json;
	class TSCTimer {
		/* TSC Timer class
		 * This class provides a high-resolution timer based on the x86 Time Step Counter.
		 * There is no public constructor for TSCTimer objects; instead users are expected
		 * to obtain an instance reference by calling TSCTimer::get_timer(filename).
		 * This allows flexible use of a single timer within a project, without passing around
		 * TSCTimer objects.
		 * 
		 * Timings are performed in the following fashion:
		 *  TSCTimer::init_timer("timings.json", 1400);  // optional, set initial timer step
		 * 	 for (...) {
		 *   	TSCTimer& timer = TSCTimer::get_timer("timings.json");
		 *   	timer.start_timing("lu-solver");
		 *   	x = lu_solver.solve (A, b);
		 *   	timer.stop_timing("lu-solver", true, "some-custom-tag");
		 *   	// more stuff
		 *
		 *   	timer.step();
		 *   }
		 * The timer may be deactivated to reduce overhead, and re-activated later:
		 *   timer.deactivate();
		 *   timer.activate();
		 *
		 */
	private:
		const std::string filename;
		using time_map_t = std::map<std::string, myInt64>;
	    time_map_t start_times;
		int time_step;
		// Sorry for Singleton
		TSCTimer (const std::string& filename, int initial_time_step) : filename(filename), time_step{initial_time_step} { }

		using timer_map_t = std::map<std::string, TSCTimer>;
		inline static timer_map_t timers;
		bool active = true;

	public:
		TSCTimer(const TSCTimer& t) = delete;
		TSCTimer(TSCTimer&& t) = default;
		static TSCTimer& init_timer(const std::string& filename, int initial_time_step=0) {
			// Ensure timer has not been initialized
			// improvement: index by file, not by file path
			std::cout << "Init timer " << filename << std::endl;
			assert (timers.find(filename) == timers.end());
			auto timer_it_pair = timers.emplace(filename, TSCTimer(filename, initial_time_step));
			return timer_it_pair.first->second;
		}

		static TSCTimer& get_timer(const std::string& filename) {
			// improvement: index by file, not by file path
			if (timers.find(filename) == timers.end()) {
				return init_timer(filename);
			}
			return timers.find(filename)->second;
		}

		void start_timing (const std::string& section_name) {
			if (!active) return;
			// start a timing
			// ensure that timer isn't still running
			assert (start_times.find(section_name) == start_times.end());
			myInt64& start_time = start_times[section_name];
			start_time = start();
		}

		// stop the timing of a particular section and return the duration
		myInt64 stop_timing (const std::string& section_name, bool print_info, const std::string& tag) {
			if (!active) return 0;
			myInt64 stop_time = start();
			assert (start_times.find(section_name) != start_times.end());

			myInt64 start_time = start_times[section_name];
			myInt64 duration = stop_time - start_time;

			// write to json file
			json summary;
			summary["step"] = time_step;
			summary["name"] = section_name;
			summary["start"] = start_time;
			summary["stop"] = stop_time;
			summary["duration"] = duration;
			summary["tag"] = tag;
			std::ofstream(filename, std::ios::app) << summary << std::endl;
			// step_no, section_name, start_time, stop_time, duration, tag

			start_times.erase(section_name);
			if (print_info) {
				std::cout << "Cycles spent in " << section_name << ": " << thousandsep(std::to_string(duration)) << std::endl;
			}
			return duration;
		}

		// increase the step counter by 1
		int step() {
			return ++time_step;
		}
		void activate() {
			active = true;
		}
		void deactivate() {
			active = false;
		}
	};
}
